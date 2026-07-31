#pragma once
// SOURCE: Hagan, Kumar, Lesniewski, Woodward (2002) "Managing Smile Risk"
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" Ch.3
//
// 模块: SABR 随机过程 (Stochastic Alpha Beta Rho)
//
// ==================== SABR 过程数学 ====================
//
// SABR 模型 SDE:
//   dF(t) = σ(t) · F(t)^β · dW_1(t)
//   dσ(t) = ν · σ(t) · dW_2(t)
//   dW_1 dW_2 = ρ dt
//   F(0) = F_0 > 0,  σ(0) = α > 0
//
// 离散化方案 (Euler + log-Euler 混合):
//   - σ(t) 严格为正 (CEV 模型要求), 用 log-Euler 保证:
//       σ_{t+dt} = σ_t · exp(-0.5·ν²·dt + ν·√dt·Z_2)
//   - F(t) 用 Euler:
//       F_{t+dt} = F_t + σ_t · F_t^β · √dt · Z_1
//     β ∈ (0, 1) 时 F 可能为负, 加吸收壁 F < 0 → F = 0
//   - β = 1 (对数正态极限): 用 log-Euler 更稳定
//       F_{t+dt} = F_t · exp(-0.5·σ_t²·dt + σ_t·√dt·Z_1)
//
// 相关性:
//   Z_1 = W_1
//   Z_2 = ρ·W_1 + √(1-ρ²)·W_2   (W_1, W_2 独立标准正态)
//
// 维度: 2 (F, σ)
// 注: SABR 无简单闭式特征函数, characteristic_function() 返回默认值

#include <cmath>
#include <span>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/models/diffusion/process.hpp"
#include "cpphub/pricing/analytic/sabr_hagan.hpp"

namespace cpphub {
inline namespace v1 {

// ============ SABR 过程模拟方案 ============
enum class SABRScheme {
    EulerAbsorbing,   // F 用 Euler + 吸收壁, σ 用 log-Euler (默认, 通用)
    LogEuler          // F 和 σ 都用 log-Euler (仅 β=1 严格正确, β<1 为近似)
};

// ============ SABR 过程类 ============
class SABRProcess : public StochasticProcess {
public:
    explicit SABRProcess(SABRParams p, Real F0, Real r = 0.0, Real q = 0.0,
                          SABRScheme scheme = SABRScheme::EulerAbsorbing)
        : params_(p), F0_(F0), r_(r), q_(q), scheme_(scheme) {
        validate_sabr_params(p);
        if (F0 <= 0.0) {
            throw std::invalid_argument("SABRProcess: F0 must be positive");
        }
    }

    Size dimension() const override { return 2; }
    Real spot() const override { return F0_; }

    // SABR 无简单闭式 CF, 返回默认 (子类可重写)
    Complex characteristic_function(Complex u, Real tau) const override {
        (void)u; (void)tau;
        return Complex{0.0, 0.0};
    }

    void generate_path(Real T, Size n_steps,
                       std::span<Real> path,
                       Philox4x64& rng) const override;

    // 一步演化 (供 MC 引擎/路径相关期权使用)
    // 输入: 当前 F, 当前 sigma, 时间步长 dt, 两个相关正态 Z1, Z2
    // 输出: 下一步的 F (sigma 通过引用更新)
    Real evolve(Real F, Real& sigma, Real dt, Real Z1, Real Z2) const;

    const SABRParams& params() const { return params_; }
    Real forward0() const { return F0_; }
    Real r() const { return r_; }
    Real q() const { return q_; }

private:
    SABRParams params_;
    Real F0_;
    Real r_;
    Real q_;
    SABRScheme scheme_;

    // 生成相关正态对 (Z1, Z2)
    void generate_correlated_normals(Philox4x64& rng, Real& Z1, Real& Z2) const;
};

// ============ 相关正态生成 ============
inline void SABRProcess::generate_correlated_normals(Philox4x64& rng,
                                                       Real& Z1, Real& Z2) const {
    uint64_t r1 = rng();
    uint64_t r2 = rng();
    double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
    double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
    auto [w1, w2] = box_muller(u1, u2);
    Real rho = params_.rho;
    Z1 = w1;
    Z2 = rho * w1 + std::sqrt(1.0 - rho * rho) * w2;
}

// ============ 一步演化 ============
inline Real SABRProcess::evolve(Real F, Real& sigma, Real dt,
                                  Real Z1, Real Z2) const {
    const Real sqrt_dt = std::sqrt(dt);
    const Real nu = params_.nu;
    const Real beta = params_.beta;

    // σ 用 log-Euler: σ_{t+dt} = σ_t · exp(-0.5·ν²·dt + ν·√dt·Z_2)
    sigma = sigma * std::exp(-0.5 * nu * nu * dt + nu * sqrt_dt * Z2);

    // F 用 Euler (吸收壁) 或 log-Euler
    if (scheme_ == SABRScheme::LogEuler) {
        // log-Euler: ln F_{t+dt} = ln F_t - 0.5·σ²·F^(2β-2)·dt + σ·F^(β-1)·√dt·Z_1
        // 仅 β=1 严格, β<1 为近似 (忽略 β 对 F 的影响)
        if (F <= 0.0) return 0.0;
        Real vol_F = sigma * std::pow(F, beta - 1.0);
        Real drift = -0.5 * vol_F * vol_F * dt;
        F = F * std::exp(drift + vol_F * sqrt_dt * Z1);
    } else {
        // Euler: F_{t+dt} = F_t + σ_t·F_t^β·√dt·Z_1
        if (F <= 0.0) return 0.0;
        Real dF = sigma * std::pow(F, beta) * sqrt_dt * Z1;
        F = F + dF;
        if (F < 0.0) F = 0.0;  // 吸收壁
    }
    return F;
}

// ============ 完整路径生成 ============
inline void SABRProcess::generate_path(Real T, Size n_steps,
                                        std::span<Real> path,
                                        Philox4x64& rng) const {
    if (n_steps == 0) {
        throw std::invalid_argument("SABRProcess: n_steps must be positive");
    }
    if (path.size() < n_steps + 1) {
        throw std::invalid_argument("SABRProcess: path buffer too small");
    }

    Real dt = T / static_cast<Real>(n_steps);
    path[0] = F0_;

    Real F = F0_;
    Real sigma = params_.alpha;

    for (Size i = 0; i < n_steps; ++i) {
        Real Z1, Z2;
        generate_correlated_normals(rng, Z1, Z2);
        F = evolve(F, sigma, dt, Z1, Z2);
        path[i + 1] = F;
    }
}

}  // namespace v1
}  // namespace cpphub
