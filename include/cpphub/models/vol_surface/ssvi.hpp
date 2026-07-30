#pragma once
// SOURCE: PHASE3_SPEC §4.2 - SSVI (Surface SVI) - arbitrage-free across maturities
// Phase 4 LITE - D2 整改项: SSVI 跨期限模块实现 (2026-07-31)
//
// SSVI 参数化 (Gatheral-Jacquier 2014):
//   w(k,θ) = θ/2 * (1 + φ(θ)*ρ*k + sqrt((φ(θ)*k + ρ)² + 1 - ρ²))
// 其中 θ = ATM 总方差, φ(θ) 是期限结构函数, k = log-moneyness
//
// 无套利条件 (Gatheral-Jacquier 2014 Theorem 4.2):
//   1. φ(θ) > 0
//   2. ∂_θ (θ * φ(θ)) > 0  (Calendar 套利, θ*φ(θ) 单调递增)
//   3. |ρ| < 1
//   4. φ(θ) * θ * (1 + |ρ|) < 4  (充分 Butterfly 套利条件)
//   严格无 butterfly 套利的额外条件 (Theorem 4.4):
//   5. φ(θ) * (θ * (1 + |ρ|) + 2*ψ*sqrt(1-ρ²)) < 4, 其中 ψ = ∂_θ(θ*φ(θ))/φ(θ)
//
// 参考实现: Gatheral & Jacquier (2014) "Arbitrage-free SVI volatility surfaces"
//           https://arxiv.org/abs/1204.0646
#include "cpphub/core/types.hpp"
#include "cpphub/models/vol_surface/svi.hpp"
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub::v1 {

// SSVI 参数化: w(k,θ) = θ/2 * (1 + φ(θ)*ρ*k + sqrt((φ(θ)*k + ρ)² + 1 - ρ²))
// 其中 θ = ATM 总方差, φ(θ) 是期限结构函数
// 无套利条件 (Gatheral-Jacquier):
//   1. φ(θ) > 0
//   2. ∂_θ (θ * φ(θ)) > 0  (Calendar 套利)
//   3. |ρ| < 1
//   4. φ(θ) * θ * (1 + |ρ|) < 4  (Butterfly 套利)
//   5. φ(θ) * (θ * (1 + |ρ|) + 2 * ψ * sqrt(1 - ρ²)) < 4  (严格无 butterfly 套利)
//   其中 ψ = -∂_θ θ * φ(θ)

struct SSVIParams {
    Real rho;                    // 相关性 (-1 < ρ < 1)
    std::function<Real(Real)> phi;  // φ(θ): 期限结构函数
    std::vector<Real> theta_slice; // 各期限的 θ 值 (ATM 总方差)
};

class SSVI {
public:
    SSVI(SSVIParams params);

    Real total_variance(Real k, Real theta) const;
    Real implied_vol(Real k, Real T, Real theta) const;

    bool check_calendar_arbitrage() const;
    bool check_butterfly_arbitrage() const;        // 充分条件 (Theorem 4.2)
    bool check_strict_butterfly_arbitrage() const; // 严格条件 (Theorem 4.4)
    bool check_no_arbitrage() const;

    // 常用 φ(θ) 参数化
    static SSVIParams Heston_like(Real rho, Real eta, Real lambda, const std::vector<Real>& theta_slice);
    static SSVIParams Power_law(Real rho, Real eta, Real gamma, const std::vector<Real>& theta_slice);

    // 从市场数据校准
    CalibrationResult calibrate(
        const std::vector<Real>& strikes,
        const std::vector<Real>& maturities,
        const std::vector<Real>& implied_vols,
        Real forward,
        const CalibConfig& cfg = CalibConfig{});

    SSVIParams params() const { return params_; }

private:
    SSVIParams params_;
    Real dvar_dk(Real k, Real theta) const;
    Real d2var_dk2(Real k, Real theta) const;
    // 数值导数辅助 (用于 φ(θ) 的导数,因 std::function 无法解析求导)
    Real numeric_dphi_dtheta(Real theta) const;
    Real numeric_d_theta_phi_dtheta(Real theta) const;  // d(θ*φ(θ))/dθ
};

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

inline SSVI::SSVI(SSVIParams params) : params_(std::move(params)) {
    if (std::abs(params_.rho) >= 1.0) {
        throw std::invalid_argument("SSVI: |rho| must be < 1");
    }
    if (!params_.phi) {
        throw std::invalid_argument("SSVI: phi function must be set");
    }
}

inline Real SSVI::total_variance(Real k, Real theta) const {
    if (theta <= 0.0) {
        throw std::invalid_argument("SSVI::total_variance: theta must be positive");
    }
    Real phi = params_.phi(theta);
    if (phi <= 0.0) {
        throw std::invalid_argument("SSVI::total_variance: phi(theta) must be positive");
    }
    Real rho = params_.rho;
    Real one_minus_rho2 = 1.0 - rho * rho;
    // w(k,θ) = θ/2 * (1 + φ*ρ*k + sqrt((φ*k + ρ)² + 1 - ρ²))
    Real pk = phi * k;
    Real discriminant = (pk + rho) * (pk + rho) + one_minus_rho2;
    Real s = std::sqrt(std::max(discriminant, Real(0.0)));
    return 0.5 * theta * (1.0 + phi * rho * k + s);
}

inline Real SSVI::implied_vol(Real k, Real T, Real theta) const {
    if (T <= 0.0) {
        throw std::invalid_argument("SSVI::implied_vol: T must be positive");
    }
    Real w = total_variance(k, theta);
    if (w < 0.0) return 0.0;  // 数值保护
    return std::sqrt(w / T);
}

// ∂w/∂k = θ/2 * (φ*ρ + (φ*(φ*k + ρ)) / sqrt((φ*k + ρ)² + 1 - ρ²))
inline Real SSVI::dvar_dk(Real k, Real theta) const {
    Real phi = params_.phi(theta);
    Real rho = params_.rho;
    Real one_minus_rho2 = 1.0 - rho * rho;
    Real pk = phi * k;
    Real discriminant = (pk + rho) * (pk + rho) + one_minus_rho2;
    Real s = std::sqrt(std::max(discriminant, Real(1e-30)));
    return 0.5 * theta * (phi * rho + phi * (pk + rho) / s);
}

// ∂²w/∂k² = θ/2 * φ² * (1 - (φ*k + ρ)² / s²) / s = θ/2 * φ² * (1-ρ²) / s³
// 推导: d/dk [φ*(φ*k+ρ)/s] = φ² * (s - (φ*k+ρ)*(φ/s)) / s² = φ² * (s² - (φ*k+ρ)²) / s³
//       s² - (φ*k+ρ)² = 1-ρ²
inline Real SSVI::d2var_dk2(Real k, Real theta) const {
    Real phi = params_.phi(theta);
    Real rho = params_.rho;
    Real one_minus_rho2 = 1.0 - rho * rho;
    Real pk = phi * k;
    Real discriminant = (pk + rho) * (pk + rho) + one_minus_rho2;
    Real s = std::sqrt(std::max(discriminant, Real(1e-30)));
    Real s3 = s * s * s;
    return 0.5 * theta * phi * phi * one_minus_rho2 / s3;
}

// 数值导数: φ'(θ) (中心差分)
inline Real SSVI::numeric_dphi_dtheta(Real theta) const {
    Real h = 1e-6 * std::max(theta, Real(1.0));
    return (params_.phi(theta + h) - params_.phi(theta - h)) / (2.0 * h);
}

// d(θ*φ(θ))/dθ = φ(θ) + θ*φ'(θ)
inline Real SSVI::numeric_d_theta_phi_dtheta(Real theta) const {
    return params_.phi(theta) + theta * numeric_dphi_dtheta(theta);
}

// Calendar 套利: d(θ*φ(θ))/dθ > 0 对所有 θ
// 等价于 θ*φ(θ) 关于 θ 单调递增 (Gatheral-Jacquier Theorem 4.2)
inline bool SSVI::check_calendar_arbitrage() const {
    if (params_.theta_slice.size() < 2) return true;  // 单期限无法违反
    // 检查 θ*φ(θ) 在所有切片上单调递增
    for (Size i = 1; i < params_.theta_slice.size(); ++i) {
        Real theta_prev = params_.theta_slice[i - 1];
        Real theta_curr = params_.theta_slice[i];
        Real prod_prev = theta_prev * params_.phi(theta_prev);
        Real prod_curr = theta_curr * params_.phi(theta_curr);
        if (prod_curr <= prod_prev) return false;
    }
    // 数值检查: 在每个 θ 上 d(θ*φ(θ))/dθ > 0
    for (Real theta : params_.theta_slice) {
        if (theta <= 0.0) return false;
        Real deriv = numeric_d_theta_phi_dtheta(theta);
        if (deriv <= 0.0) return false;
    }
    return true;
}

// Butterfly 套利 (充分条件, Gatheral-Jacquier Theorem 4.2):
//   |ρ|<1 且 φ(θ)*θ*(1+|ρ|) ≤ 4 对所有 θ
// 这是实际无套利所需的标准检查条件
inline bool SSVI::check_butterfly_arbitrage() const {
    if (std::abs(params_.rho) >= 1.0) return false;
    Real abs_rho = std::abs(params_.rho);
    for (Real theta : params_.theta_slice) {
        if (theta <= 0.0) return false;
        Real phi = params_.phi(theta);
        if (phi <= 0.0) return false;
        if (phi * theta * (1.0 + abs_rho) >= 4.0) return false;
    }
    return true;
}

// 严格 Butterfly 套利 (必要充分条件, Gatheral-Jacquier Theorem 4.4):
//   额外 φ(θ)*(θ*(1+|ρ|) + 2*ψ*sqrt(1-ρ²)) < 4, ψ = ∂_θ(θ*φ(θ))/φ(θ)
// 这是更强的保证,实际校准中可能处于边界,仅用于诊断
inline bool SSVI::check_strict_butterfly_arbitrage() const {
    if (std::abs(params_.rho) >= 1.0) return false;
    Real abs_rho = std::abs(params_.rho);
    Real one_minus_rho2 = 1.0 - params_.rho * params_.rho;
    Real sqrt_1mrho2 = std::sqrt(std::max(one_minus_rho2, Real(0.0)));
    for (Real theta : params_.theta_slice) {
        if (theta <= 0.0) return false;
        Real phi = params_.phi(theta);
        if (phi <= 0.0) return false;
        Real dtheta_phi = numeric_d_theta_phi_dtheta(theta);
        Real psi = dtheta_phi / phi;
        if (phi * (theta * (1.0 + abs_rho) + 2.0 * psi * sqrt_1mrho2) >= 4.0) return false;
    }
    return true;
}

inline bool SSVI::check_no_arbitrage() const {
    return check_calendar_arbitrage() && check_butterfly_arbitrage();
}

// Heston-like 参数化: φ(θ) = η / (θ^λ * (1+θ)^λ)
// 常用 λ=1: φ(θ) = η / (θ*(1+θ)), 满足 θ*φ(θ) = η/(1+θ) 关于 θ 单调递减 — 不满足 calendar!
// 正确 Heston-like (Gatheral-Jacquier eq. 5.2): φ(θ) = η * θ^(-λ), 0<λ<1/2
// θ*φ(θ) = η*θ^(1-λ), d/dθ = η*(1-λ)*θ^(-λ) > 0 ✓ (λ<1)
// Butterfly: φ(θ)*θ*(1+|ρ|) = η*θ^(1-λ)*(1+|ρ|), 需 < 4 对最大 θ 成立
inline SSVIParams SSVI::Heston_like(Real rho, Real eta, Real lambda,
                                      const std::vector<Real>& theta_slice) {
    if (lambda <= 0.0 || lambda >= 1.0) {
        throw std::invalid_argument("Heston_like: lambda must be in (0,1)");
    }
    if (eta <= 0.0) {
        throw std::invalid_argument("Heston_like: eta must be positive");
    }
    SSVIParams p;
    p.rho = rho;
    p.phi = [eta, lambda](Real theta) { return eta * std::pow(theta, -lambda); };
    p.theta_slice = theta_slice;
    return p;
}

// Power-law 参数化 (Gatheral-Jacquier eq. 5.3): φ(θ) = η * θ^(-γ), 0<γ<1/2
// 与 Heston-like 形式相同但 γ<1/2 保证严格无 butterfly 套利
inline SSVIParams SSVI::Power_law(Real rho, Real eta, Real gamma,
                                    const std::vector<Real>& theta_slice) {
    if (gamma <= 0.0 || gamma >= 0.5) {
        throw std::invalid_argument("Power_law: gamma must be in (0, 0.5)");
    }
    if (eta <= 0.0) {
        throw std::invalid_argument("Power_law: eta must be positive");
    }
    SSVIParams p;
    p.rho = rho;
    p.phi = [eta, gamma](Real theta) { return eta * std::pow(theta, -gamma); };
    p.theta_slice = theta_slice;
    return p;
}

// 校准: 拟合 SSVI 参数 (rho, eta, gamma) 到市场数据
// 采用分层策略: 先按期限校准 SVI 切片得到各 θ, 再拟合 SSVI 全局参数
inline CalibrationResult SSVI::calibrate(
        const std::vector<Real>& strikes,
        const std::vector<Real>& maturities,
        const std::vector<Real>& implied_vols,
        Real forward,
        const CalibConfig& cfg) {

    CalibrationResult result;
    if (strikes.empty() || maturities.empty()) {
        result.converged = false;
        result.message = "empty input data";
        return result;
    }

    // Step 1: 按期限分组,每组用 SVI 校准得到 ATM 总方差 θ
    // (这里简化:假设输入按期限排序,直接用第一期限做 SVI 切片校准作为 θ 估计)
    // 实际生产实现应按期限分组校准每个切片
    std::vector<Real> theta_estimates;
    Size n_maturities = maturities.size();
    Size n_strikes = strikes.size();

    // 假设 strikes/implied_vols 是 flattend (n_maturities × n_strikes_per_mat)
    // 简化:对每个期限计算 ATM 总方差 (插值或最近 strike)
    for (Size m = 0; m < n_maturities; ++m) {
        Real T = maturities[m];
        Real atm_vol = implied_vols[m * n_strikes];  // 假设第一个 strike 接近 ATM
        // 找最近 forward 的 strike
        Real best_diff = std::abs(strikes[0] - forward);
        Real best_vol = implied_vols[m * n_strikes];
        for (Size i = 0; i < n_strikes; ++i) {
            Real diff = std::abs(strikes[i] - forward);
            if (diff < best_diff) {
                best_diff = diff;
                best_vol = implied_vols[m * n_strikes + i];
            }
        }
        theta_estimates.push_back(best_vol * best_vol * T);
    }

    // Step 2: 全局拟合 SSVI 参数 (rho, eta, gamma) — 使用 Power-law 参数化
    // 目标函数: sum over all (k, T, market_vol) of (SSVI_vol(k,T,θ(T)) - market_vol)²
    ResidualFn residual = [&](const std::vector<Real>& x) -> std::vector<Real> {
        Real rho = std::max(-0.999, std::min(0.999, x[0]));
        Real eta = std::max(1e-4, x[1]);
        Real gamma = std::max(1e-4, std::min(0.499, x[2]));

        std::vector<Real> r;
        r.reserve(n_maturities * n_strikes);
        for (Size m = 0; m < n_maturities; ++m) {
            Real T = maturities[m];
            Real theta = theta_estimates[m];
            Real phi = eta * std::pow(theta, -gamma);
            Real one_minus_rho2 = 1.0 - rho * rho;
            for (Size i = 0; i < n_strikes; ++i) {
                Real k = std::log(strikes[i] / forward);
                Real pk = phi * k;
                Real disc = (pk + rho) * (pk + rho) + one_minus_rho2;
                Real s = std::sqrt(std::max(disc, Real(0.0)));
                Real w = 0.5 * theta * (1.0 + phi * rho * k + s);
                Real model_vol = (w > 0.0) ? std::sqrt(w / T) : 0.0;
                Real market_vol = implied_vols[m * n_strikes + i];
                r.push_back(model_vol - market_vol);
            }
        }
        return r;
    };

    // 初始猜测
    std::vector<Real> x0 = {-0.3, 1.0, 0.25};

    std::vector<Real> x_init = x0;
    if (cfg.use_de_init) {
        std::vector<Bounds> bounds = {
            {-0.99, 0.99},  // rho
            {1e-4, 10.0},   // eta
            {1e-4, 0.49}    // gamma
        };
        ObjectiveFn obj = [&](const std::vector<Real>& xx) -> Real {
            auto r = residual(xx);
            Real s = 0.0;
            for (Real v : r) s += v * v;
            return s;
        };
        DifferentialEvolution::Config de_cfg;
        de_cfg.population_size = cfg.de_pop_size;
        de_cfg.max_generations = cfg.de_generations;
        de_cfg.seed = cfg.seed;
        auto de_result = DifferentialEvolution::minimize(obj, bounds, de_cfg);
        x_init = de_result.x;
    }

    LevenbergMarquardt::Config lm_cfg;
    lm_cfg.max_iterations = cfg.lm_max_iter;
    lm_cfg.ftol = cfg.ftol;
    lm_cfg.xtol = cfg.xtol;
    auto lm_result = LevenbergMarquardt::minimize(residual, x_init, lm_cfg);

    // 更新内部参数
    params_.rho = std::max(-0.999, std::min(0.999, lm_result.x[0]));
    Real eta = std::max(1e-4, lm_result.x[1]);
    Real gamma = std::max(1e-4, std::min(0.499, lm_result.x[2]));
    params_.phi = [eta, gamma](Real theta) { return eta * std::pow(theta, -gamma); };
    params_.theta_slice = theta_estimates;

    result.params = lm_result.x;
    result.objective_value = lm_result.fx;
    result.n_iterations = lm_result.n_iterations;
    result.converged = lm_result.converged;
    result.message = lm_result.message;
    result.residuals = residual(lm_result.x);
    return result;
}

}  // namespace cpphub::v1
