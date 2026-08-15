// =============================================================================
// egarch_model.hpp - EGARCH(1,1) QMLE 估计 (spec §2.0.3)
//
// Phase 7B v1.6 M1 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 教材锚点: Nelson 1991 / Tsay 3ed Ch 6 / arch arch/univariate/volatility.py
// 优化器: SLSQP (ADR-018, calibration/optimizer.hpp)
//
// 幻觉点防护 (spec §6.1 + 实测):
//   G5:  E|z| = √(2/π) ≈ 0.7978845608 (非 2/π)
//   G6:  非对称项用标准化残差 zₜ₋₁ = εₜ₋₁/√hₜ₋₁ (非未标准化 εₜ₋₁)
//   G23: 参数顺序与 arch 相反: spec.alpha = 非对称 z 项 = arch gamma[o];
//        spec.gamma = 对称 (|z|-E|z|) 项 = arch alpha[p] (Nelson 1991 约定)
//   G-q (实测, 源码核对 2026-08-15): arch EGARCH.backcast() 重写为
//        log(super().backcast(resids)) — backcast 返回 **log 尺度**;
//        egarch_recursion_python t=0: lnsigma2[0] = ω + β·backcast(log 尺度)
//        即 ln h₁ = ω + β·log(EWMA 方差 backcast)。verify_egarch.py 基准
//        backcast=-0.4659 (负值) 印证 log 尺度。
//
// 无内部尺度变换: arch 的 t=0 约定 (β·bc 非线性进入 log 递归) 破坏 MLE 尺度
// 等变性, 缩放后估计的是不同模型 (t=1 项偏移 β·(s²-1)·bc - 2·ln s 量级)。
// 故 EGARCH 在原始去均值尺度直接估计 (与 arch 逐位一致)。数据尺度需使
// ln(无条件方差) ∈ (-10, 10) (spec 边界), 日频/分钟频均满足。
// =============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <algorithm>

#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/calibration/optimizer.hpp"
#include "cpphub/timeseries/garch/garch_model.hpp"      // GarchConfig/detail
#include "cpphub/timeseries/garch/garch_distribution.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace garch {

/// E|z| 常数 (G5, 公开以便测试)
inline constexpr Real EGARCH_E_ABS_Z = 0.7978845608028654;  // √(2/π)

// EGARCH 参数 (Nelson 1991, G23 参数映射)
struct EGarchParams {
    Real omega;   ///< 常数项 ω
    Real alpha;   ///< 非对称系数 (θ in Nelson; 乘 zₜ₋₁; 对应 arch gamma[o])
    Real beta;    ///< 持续性 β (|β| < 1 平稳)
    Real gamma;   ///< 对称 size 系数 (γ in Nelson; 乘 (|zₜ₋₁|-E|z|); arch alpha[p])
    Real nu;      ///< t/GED 自由度 (Normal 时忽略)
};

// EGARCH 估计结果 (结构同 GarchResult, 参数替换为 EGarchParams)
struct EGarchResult {
    EGarchParams params;
    std::vector<Real> conditional_variances;
    std::vector<Real> residuals;
    std::vector<Real> std_residuals;
    Real log_likelihood;
    Real aic;
    Real bic;
    std::vector<std::vector<Real>> vcov;
    std::vector<Real> std_errors;
    bool converged;
    Size n_iterations;
    std::string message;
};

/// @brief EGARCH(1,1) 方差递归 (公开以便测试, 与 arch egarch_recursion 一致)
///
/// ln h₁ = ω + β·bc_log  (G-q: bc_log 为 **log 尺度** backcast, 即
///        log(EWMA 方差), arch EGARCH.backcast 重写返回 log, 可为负)
/// ln hₜ = ω + β·ln hₜ₋₁ + α·zₜ₋₁ + γ·(|zₜ₋₁| - E|z|)  (t ≥ 2, G6: z 标准化)
inline std::vector<Real> filter_egarch(const EGarchParams& params,
                                       const std::vector<Real>& residuals,
                                       Real bc_log) {
    const Size T = residuals.size();
    if (!std::isfinite(bc_log)) {
        throw std::invalid_argument("filter_egarch: backcast must be finite");
    }
    std::vector<Real> h(T);
    Real ln_h_prev = params.omega + params.beta * bc_log;  // G-q
    for (Size t = 0; t < T; ++t) {
        if (t > 0) {
            const Real z_prev = residuals[t - 1] / std::sqrt(h[t - 1]);
            ln_h_prev = params.omega + params.beta * ln_h_prev
                        + params.alpha * z_prev
                        + params.gamma * (std::abs(z_prev) - EGARCH_E_ABS_Z);
        }
        h[t] = std::exp(ln_h_prev);
        if (!std::isfinite(h[t]) || !(h[t] > 0.0)) {
            throw std::runtime_error("filter_egarch: invalid variance");
        }
    }
    return h;
}

/// @brief EGARCH(1,1) QMLE 估计 (原尺度, spec §2.0.3)
///
/// @throws std::invalid_argument 若 T < 10, 含 NaN, 或零方差
inline EGarchResult estimate_egarch(const std::vector<Real>& data,
                                    const GarchConfig& config = GarchConfig{}) {
    const Size T = data.size();
    if (T < 10) {
        throw std::invalid_argument("estimate_egarch: need at least 10 observations");
    }
    Real mean = 0.0;
    for (Real v : data) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("estimate_egarch: data contains NaN/Inf");
        }
        mean += v;
    }
    mean /= static_cast<Real>(T);

    // G2: 残差 = 数据 - 样本均值
    std::vector<Real> eps(T);
    for (Size t = 0; t < T; ++t) eps[t] = data[t] - mean;

    Real var = 0.0;
    for (Real e : eps) var += e * e;
    var /= static_cast<Real>(T);
    Real scale = 0.0;
    for (Real v : data) scale = std::max(scale, std::abs(v));
    const Real rel_sd = (scale > 0.0) ? std::sqrt(var) / scale : 0.0;
    if (!(var > 0.0) || rel_sd < 1e-13) {
        throw std::invalid_argument("estimate_egarch: zero variance data");
    }

    const bool has_nu = config.dist != GarchDist::Normal;
    const Size k = has_nu ? 5u : 4u;

    // G1: backcast (log 尺度, G-q: arch EGARCH.backcast = log(EWMA 方差))
    const Real bc_log = std::log(backcast_variance(eps, config.backcast_lambda,
                                                    config.backcast_window));

    // 逐观测似然 (原尺度递归)
    detail::PerObsLL pll = [&](const std::vector<Real>& theta) {
        const Real nu = has_nu ? theta[4] : 0.0;
        const EGarchParams p{theta[0], theta[1], theta[2], theta[3], nu};
        const std::vector<Real> h = filter_egarch(p, eps, bc_log);
        std::vector<Real> ll(T);
        for (Size t = 0; t < T; ++t) {
            ll[t] = log_likelihood_term(eps[t], h[t], config.dist, nu);
        }
        return ll;
    };

    // 约束 + 边界 (spec §2.0.3 Step 3): β ∈ [0,1), 其余 [-10,10]/[-1,1]
    std::vector<Bounds> bnds = {{-10.0, 10.0}, {-1.0, 1.0},
                                {0.0, 1.0 - 1e-6}, {-1.0, 1.0}};
    if (has_nu) bnds.push_back({2.05, 500.0});
    const ConstraintFn ineq = [](const std::vector<Real>& x) {
        return std::vector<Real>{1.0 - x[2] - 1e-8};  // |β| < 1 (β≥0 由 bounds)
    };

    // G16 多起始 (原尺度): HR/ML/用户/随机扰动
    const Real log_var = std::log(var);
    const Real nu0 = 8.0;
    std::vector<std::vector<Real>> starts;
    auto push_start = [&](Real om, Real al, Real be, Real ga) {
        std::vector<Real> st = {om, al, be, ga};
        if (has_nu) st.push_back(nu0);
        starts.push_back(std::move(st));
    };
    push_start(0.10 * log_var, -0.05, 0.90, 0.15);   // HR (含杠杆方向)
    push_start(0.05 * log_var, 0.00, 0.95, 0.05);    // ML
    if (config.initial_params.size() >= 4) {
        std::vector<Real> st(config.initial_params.begin(),
                             config.initial_params.begin() + 4);
        if (has_nu) {
            st.push_back(config.initial_params.size() >= 5
                             ? config.initial_params[4] : nu0);
        }
        starts.push_back(std::move(st));
    }
    {  // 随机扰动 (HR 基础 + 10% 高斯, Philox 确定性)
        Philox4x64 rng(43, 11);
        std::vector<Real> st = {0.10 * log_var, -0.05, 0.90, 0.15};
        if (has_nu) st.push_back(nu0);
        for (Size j = 0; j < st.size(); ++j) {
            const uint64_t r1 = rng();
            const uint64_t r2 = rng();
            const Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            const Real u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            const auto [z1, z2] = box_muller(u1, u2);
            (void)z2;
            st[j] += 0.10 * std::abs(st[j]) * z1;
        }
        starts.push_back(std::move(st));
    }
    if (!config.use_multistart) starts.resize(1);

    const detail::MultistartResult fit =
        detail::run_multistart(pll, starts, bnds, ineq, config.optimizer_config);

    // 最优参数处重组 (原尺度, 无映射)
    const EGarchParams phat{fit.x[0], fit.x[1], fit.x[2], fit.x[3],
                            has_nu ? fit.x[4] : 0.0};
    const std::vector<Real> h = filter_egarch(phat, eps, bc_log);
    Real ll = 0.0;
    for (Size t = 0; t < T; ++t) {
        ll += log_likelihood_term(eps[t], h[t], config.dist, phat.nu);
    }

    EGarchResult result{};
    result.params = phat;
    result.conditional_variances = h;
    result.residuals = eps;
    result.std_residuals.resize(T);
    for (Size t = 0; t < T; ++t) {
        result.std_residuals[t] = eps[t] / std::sqrt(h[t]);
    }
    result.log_likelihood = ll;
    result.aic = -2.0 * ll + 2.0 * static_cast<Real>(k);
    result.bic = -2.0 * ll + static_cast<Real>(k) * std::log(static_cast<Real>(T));

    // G9: sandwich V = H⁻¹·S·H⁻¹ (arch 约定: 中心化 OPG × n/(n-1))
    if (config.compute_sandwich) {
        try {
            const auto H = detail::numerical_hessian(pll, fit.x);
            const Size npar = fit.x.size();
            const auto G = detail::opg_scores(pll, fit.x);
            const Size Tn = G.empty() ? 0 : G[0].size();
            std::vector<Real> gbar(npar, 0.0);
            for (Size a = 0; a < npar; ++a) {
                Real acc = 0.0;
                for (Size t = 0; t < Tn; ++t) acc += G[a][t];
                gbar[a] = acc / static_cast<Real>(Tn);
            }
            std::vector<std::vector<Real>> S(npar, std::vector<Real>(npar, 0.0));
            const Real scl = static_cast<Real>(Tn) / static_cast<Real>(Tn - 1);
            for (Size t = 0; t < Tn; ++t)
                for (Size a = 0; a < npar; ++a)
                    for (Size b = 0; b < npar; ++b)
                        S[a][b] += (G[a][t] - gbar[a]) * (G[b][t] - gbar[b]);
            for (Size a = 0; a < npar; ++a)
                for (Size b = 0; b < npar; ++b) S[a][b] *= scl;
            std::vector<std::vector<Real>> Hinv;
            if (!detail::invert_matrix(H, Hinv)) {
                result.message = "sandwich: Hessian inversion failed";
            } else {
                std::vector<std::vector<Real>> HS(npar, std::vector<Real>(npar, 0.0));
                for (Size a = 0; a < npar; ++a)
                    for (Size b = 0; b < npar; ++b) {
                        Real acc = 0.0;
                        for (Size m = 0; m < npar; ++m) acc += Hinv[a][m] * S[m][b];
                        HS[a][b] = acc;
                    }
                result.vcov.assign(npar, std::vector<Real>(npar, 0.0));
                for (Size a = 0; a < npar; ++a)
                    for (Size b = 0; b < npar; ++b) {
                        Real acc = 0.0;
                        for (Size m = 0; m < npar; ++m) acc += HS[a][m] * Hinv[m][b];
                        result.vcov[a][b] = acc;
                    }
                result.std_errors.resize(npar);
                for (Size a = 0; a < npar; ++a) {
                    result.std_errors[a] =
                        (result.vcov[a][a] > 0.0) ? std::sqrt(result.vcov[a][a]) : 0.0;
                }
            }
        } catch (const std::exception& e) {
            result.message = std::string("sandwich failed: ") + e.what();
        }
    }

    result.converged = fit.converged;
    result.n_iterations = fit.n_iterations;
    if (result.message.empty()) result.message = fit.message;
    return result;
}

}  // namespace garch
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
