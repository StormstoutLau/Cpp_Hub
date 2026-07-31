#pragma once
// SOURCE: El Euch, Rosenbaum (2018) "The characteristic function of rough
//         Heston model" — Mathematical Finance, 28(1), 131-158.
// SOURCE: Smith (2017) "Option pricing under the rough Heston model" —
//         fractional CIR Euler discretisation.
//
// 模块: Rough Heston 过程层 — 分数阶 CIR (fractional CIR) + GBM 价格过程
//
// ==================== Rough Heston 模型数学 ====================
//
// 方差过程服从分数阶 CIR (Volterra 型积分形式):
//   v(t) = v(0) + (1/Γ(α)) ∫₀^t (t-s)^{α-1} [κ(θ - v(s)) ds + σ√v(s) dW(s)]
//   α = H + 1/2 ∈ (0.5, 1),  H ∈ (0, 0.5)
//
// 价格过程 (风险中性, 杠杆相关):
//   dS/S = (r - q) dt + √v(t) dZ(t),  dZ = ρ dW + √(1-ρ²) dW⊥
//
// ==================== 分数阶 Euler 离散化 (Smith 2017) ====================
//
// 均匀网格 t_j = j·Δt, j = 0..N, Δt = T/N.
// 确定性 Volterra 核:
//   K_{j,i} = ∫_{t_i}^{t_{i+1}} (t_{j+1} - s)^{α-1} ds
//           = [(t_{j+1} - t_i)^α - (t_{j+1} - t_{i+1})^α] / α
//
// 方差更新 (Full Truncation, 用 v_i⁺ = max(v_i, 0)):
//   v_{j+1} = v(0) + (κθ/Γ(α))·t_{j+1}^α
//             - (κ/Γ(α))·Σᵢ₌₀ʲ K_{j,i}·vᵢ⁺
//             + (σ/Γ(α))·Σᵢ₌₀ʲ √(K_{j,i}·vᵢ⁺)·Zᵢ
//   (随机积分用 Itô 等距近似 √(K_{j,i}·vᵢ⁺)·Zᵢ, 并截断 v_{j+1} ≥ 0)
//
// 价格更新 (log-Euler, 用 v_j⁺):
//   S_{j+1} = S_j·exp((r-q-0.5·v_j⁺)·Δt + √(v_j⁺·Δt)·(ρ·Z_j + √(1-ρ²)·Z⊥_j))
//   其中 Z_j 为方差过程噪声 (杠杆相关), Z⊥_j 为独立正交噪声.
//
// ==================== H = 0.5 (α = 1) 退化 ====================
// α=1 时 K_{j,i} = Δt (常数), 方差更新退化为标准 CIR Euler (Full Truncation),
// 与 models/diffusion/heston.hpp 的 Euler scheme 逐路径一致 (给定相同噪声序列).

#include <cmath>
#include <vector>
#include <span>
#include <stdexcept>
#include <algorithm>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/models/diffusion/process.hpp"

namespace cpphub {
inline namespace v1 {

struct RoughHestonParams {
    Real H;       // Hurst 指数, (0, 0.5]; α = H + 0.5. H=0.5 退化为标准 Heston
    Real kappa;   // 均值回归速度 > 0
    Real theta;   // 长期方差 > 0
    Real sigma;   // vol of vol > 0
    Real rho;     // 价格-方差相关性 ∈ (-1, 1)
    Real v0;      // 初始方差 > 0
    Real S0;      // 标的现价 > 0
    Real r;       // 利率
    Real q;       // 股息率
};

inline void validate_rough_heston_params(const RoughHestonParams& p) {
    if (p.H <= 0.0 || p.H > 0.5)
        throw std::invalid_argument("RoughHeston: H must be in (0, 0.5]");
    if (p.kappa <= 0.0)
        throw std::invalid_argument("RoughHeston: kappa must be positive");
    if (p.theta <= 0.0)
        throw std::invalid_argument("RoughHeston: theta must be positive");
    if (p.sigma <= 0.0)
        throw std::invalid_argument("RoughHeston: sigma must be positive");
    if (p.rho <= -1.0 || p.rho >= 1.0)
        throw std::invalid_argument("RoughHeston: rho must be in (-1, 1)");
    if (p.v0 <= 0.0)
        throw std::invalid_argument("RoughHeston: v0 must be positive");
    if (p.S0 <= 0.0)
        throw std::invalid_argument("RoughHeston: S0 must be positive");
}

// 预计算 Volterra 核系数 K_{j,i} (j ≥ i 下三角)
// K[j][i] = [(t_{j+1} - t_i)^α - (t_{j+1} - t_{i+1})^α] / α
// 返回 N×N 矩阵 (j=0..N-1, i=0..j; i>j 为 0)
std::vector<std::vector<Real>> rough_heston_kernel(Real T, Size n_steps, Real alpha) {
    if (T <= 0.0) throw std::invalid_argument("rough_heston_kernel: T must be positive");
    if (n_steps == 0) throw std::invalid_argument("rough_heston_kernel: n_steps must be positive");
    if (alpha <= 0.0) throw std::invalid_argument("rough_heston_kernel: alpha must be positive");

    const Real dt = T / static_cast<Real>(n_steps);
    std::vector<std::vector<Real>> K(n_steps, std::vector<Real>(n_steps, 0.0));
    for (Size j = 0; j < n_steps; ++j) {
        const Real t_jp1 = static_cast<Real>(j + 1) * dt;
        for (Size i = 0; i <= j; ++i) {
            const Real a = std::pow(t_jp1 - static_cast<Real>(i) * dt, alpha);
            const Real b = std::pow(t_jp1 - static_cast<Real>(i + 1) * dt, alpha);
            K[j][i] = (a - b) / alpha;
        }
    }
    return K;
}

// Rough Heston 过程 (fractional CIR + GBM 价格)
class RoughHestonProcess : public StochasticProcess {
public:
    RoughHestonProcess(const RoughHestonParams& params, Real T, Size n_steps)
        : params_(params), T_(T), n_steps_(n_steps) {
        validate_rough_heston_params(params_);
        if (T_ <= 0.0) throw std::invalid_argument("RoughHestonProcess: T must be positive");
        if (n_steps_ == 0) throw std::invalid_argument("RoughHestonProcess: n_steps must be positive");

        const Real dt = T_ / static_cast<Real>(n_steps_);
        t_grid_.resize(n_steps_ + 1);
        for (Size j = 0; j <= n_steps_; ++j) t_grid_[j] = static_cast<Real>(j) * dt;
        K_ = rough_heston_kernel(T_, n_steps_, alpha());
    }

    Size dimension() const override { return 2; }  // [S, v]
    Real spot() const override { return params_.S0; }

    // 生成一条路径: path[0..n_steps] = S 路径,
    // 若 buffer 足够大, path[n_steps+1..2*(n_steps+1)-1] = v 路径.
    void generate_path(Real T, Size n_steps, std::span<Real> path,
                       Philox4x64& rng) const override {
        if (n_steps == 0) {
            path[0] = params_.S0;
            return;
        }
        if (path.size() < n_steps + 1)
            throw std::invalid_argument("RoughHestonProcess: path buffer too small");

        std::vector<std::vector<Real>> local;
        const std::vector<std::vector<Real>>* K = &K_;
        if (!(n_steps == n_steps_ && std::abs(T - T_) < 1e-12)) {
            local = rough_heston_kernel(T, n_steps, alpha());
            K = &local;
        }
        simulate(T, n_steps, *K, path, rng);
    }

    // 便捷接口: 返回 [S_path, v_path], 各为 n_steps+1 长度
    std::vector<std::vector<Real>> generate_path(Philox4x64& rng,
                                                 Real dt_scale = 1.0) const {
        const Real T_eff = T_ * dt_scale;
        std::vector<std::vector<Real>> local;
        const std::vector<std::vector<Real>>* K = &K_;
        if (std::abs(dt_scale - 1.0) >= 1e-12) {
            local = rough_heston_kernel(T_eff, n_steps_, alpha());
            K = &local;
        }

        std::vector<Real> flat(2 * (n_steps_ + 1));
        simulate(T_eff, n_steps_, *K, flat, rng);

        std::vector<std::vector<Real>> out(2, std::vector<Real>(n_steps_ + 1, 0.0));
        for (Size i = 0; i <= n_steps_; ++i) {
            out[0][i] = flat[i];
            out[1][i] = flat[n_steps_ + 1 + i];
        }
        return out;
    }

    const RoughHestonParams& params() const { return params_; }
    Real T() const { return T_; }
    Size n_steps() const { return n_steps_; }
    Real alpha() const { return params_.H + 0.5; }

private:
    RoughHestonParams params_;
    Real T_;
    Size n_steps_;
    std::vector<Real> t_grid_;
    std::vector<std::vector<Real>> K_;  // 预计算核 (下三角, j ≥ i)

    void simulate(Real T, Size n_steps, const std::vector<std::vector<Real>>& K,
                  std::span<Real> path, Philox4x64& rng) const {
        const Real alpha = params_.H + 0.5;
        const Real dt = T / static_cast<Real>(n_steps);
        const Real sqrt_dt = std::sqrt(dt);
        const Real rho = params_.rho;
        const Real sqrt_1m_rho2 = std::sqrt(1.0 - rho * rho);
        const Real inv_gamma = 1.0 / std::tgamma(alpha);
        const Real kt_gamma = params_.kappa * params_.theta * inv_gamma;
        const Real k_gamma = params_.kappa * inv_gamma;
        const Real s_gamma = params_.sigma * inv_gamma;

        // 噪声: Z[i] 为方差过程噪声, Z_perp[i] 为正交噪声 (与 Heston Euler 噪声约定一致)
        std::vector<Real> Z(n_steps);
        std::vector<Real> Z_perp(n_steps);
        for (Size i = 0; i < n_steps; ++i) {
            uint64_t r1 = rng();
            uint64_t r2 = rng();
            double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            auto [w1, w2] = box_muller(u1, u2);
            Z[i] = rho * w1 + sqrt_1m_rho2 * w2;
            Z_perp[i] = sqrt_1m_rho2 * w1 - rho * w2;
        }

        std::vector<Real> v(n_steps + 1, 0.0);
        v[0] = params_.v0;
        path[0] = params_.S0;
        Real S = params_.S0;

        for (Size j = 0; j < n_steps; ++j) {
            const Real t_jp1 = static_cast<Real>(j + 1) * dt;
            const auto& Krow = K[j];
            Real drift_sum = 0.0;
            Real diff_sum = 0.0;
            for (Size i = 0; i <= j; ++i) {
                const Real vpos = (v[i] > 0.0) ? v[i] : 0.0;
                const Real Kji = Krow[i];
                drift_sum += Kji * vpos;
                diff_sum += std::sqrt(Kji * vpos) * Z[i];
            }

            // 分数阶 Euler + Full Truncation
            Real v_next = params_.v0 + kt_gamma * std::pow(t_jp1, alpha)
                          - k_gamma * drift_sum + s_gamma * diff_sum;
            if (v_next < 0.0) v_next = 0.0;
            v[j + 1] = v_next;

            // 价格更新 (log-Euler, 用 v_j⁺, 与方差噪声 Z[j] 相关)
            const Real v_j = (v[j] > 0.0) ? v[j] : 0.0;
            const Real dZ = rho * Z[j] + sqrt_1m_rho2 * Z_perp[j];
            S *= std::exp((params_.r - params_.q - 0.5 * v_j) * dt
                          + std::sqrt(v_j) * dZ * sqrt_dt);
            path[j + 1] = S;
        }

        if (path.size() >= 2 * (n_steps + 1)) {
            for (Size i = 0; i <= n_steps; ++i) path[n_steps + 1 + i] = v[i];
        }
    }
};

// ============ MC 定价 ============
struct RoughHestonMCResult {
    Real price;
    Real std_error;
    Real ci_lower;
    Real ci_upper;
    Size n_paths;
};

// 欧式 call/put MC 定价
RoughHestonMCResult rough_heston_price_european(
    const RoughHestonParams& params,
    Real T, Real K, bool is_call,
    Size n_steps, Size n_paths,
    uint64_t seed) {
    validate_rough_heston_params(params);
    if (T <= 0.0) throw std::invalid_argument("rough_heston MC: T must be positive");
    if (K <= 0.0) throw std::invalid_argument("rough_heston MC: K must be positive");
    if (n_steps == 0) throw std::invalid_argument("rough_heston MC: n_steps must be positive");
    if (n_paths == 0) throw std::invalid_argument("rough_heston MC: n_paths must be positive");

    RoughHestonProcess process(params, T, n_steps);

    Real sum = 0.0;
    Real sum2 = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(seed + j);
        auto paths = process.generate_path(rng, 1.0);
        const Real ST = paths[0].back();
        const Real payoff = is_call ? std::max(ST - K, 0.0) : std::max(K - ST, 0.0);
        sum += payoff;
        sum2 += payoff * payoff;
    }

    const Real mean = sum / static_cast<Real>(n_paths);
    Real var = sum2 / static_cast<Real>(n_paths) - mean * mean;
    if (var < 0.0) var = 0.0;
    const Real se = std::sqrt(var / static_cast<Real>(n_paths));
    const Real df = std::exp(-params.r * T);

    RoughHestonMCResult res;
    res.price = df * mean;
    res.std_error = df * se;
    res.ci_lower = res.price - 1.96 * res.std_error;
    res.ci_upper = res.price + 1.96 * res.std_error;
    res.n_paths = n_paths;
    return res;
}

}  // namespace v1
}  // namespace cpphub
