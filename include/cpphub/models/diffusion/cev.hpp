#pragma once
// SOURCE: Cox (1975) "Notes on Option Pricing I: Constant Elasticity of Variance Diffusion"
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" §3.2-3.4
// SOURCE: Schroder (1989) "Computing the CEV Option Pricing Formula"
//
// 模块: CEV (Constant Elasticity of Variance) 随机过程
//
// ==================== CEV 过程数学 ====================
//
// CEV SDE (risk-neutral):
//   dS(t) = (r - q) S(t) dt + σ S(t)^β dW(t)
//   S(0) = S_0 > 0
//
// 参数:
//   σ (sigma) — 波动率尺度, > 0
//   β (beta)  — 弹性系数, β < 1: 吸收壁在 0; β = 1: GBM; β > 1: 可能爆炸 (本版未支持)
//
// 离散化方案:
//   (A) Euler + 吸收壁 (默认, 通用):
//       S_{t+dt} = S_t + (r-q)·S_t·dt + σ·S_t^β·√dt·Z
//       若 S_{t+dt} < 0 则 S_{t+dt} = 0  (β < 1 的吸收壁)
//
//   (B) Log-Euler (仅 β = 1 严格正确):
//       S_{t+dt} = S_t · exp( (r-q-0.5·σ²)·dt + σ·√dt·Z )
//       对 β != 1 为近似 (忽略 β 对局部波动率的影响)
//
// 吸收壁处理 (β < 1):
//   - 当 S → 0 时, dS → σ·0^β·dW = 0 (因 β > 0), S 在 0 处被吸收
//   - Euler 离散可能让 S 跨过 0, 需 floor 至 0
//   - β ≤ 0 时 0^β 发散, 不支持 (本实现要求 β ∈ [0, 1])
//
// 维度: 1
// 注: CEV 无简单闭式特征函数 (Schroder 解是经过变量替换到非中心卡方, 非标准 CF)
//
// ==================== 数值稳定性要点 ====================
//
// 1. β 接近 1 时, Euler 离散与 GBM 误差小; Log-Euler 在 β=1 严格
// 2. S^β 当 S → 0 且 β < 1 时趋于 0 (良性); β = 0 时为常数 σ (Bachelier-like)
// 3. drift 项 (r-q)·S·dt 在 S=0 处为 0, 与吸收壁一致
// 4. 大步长 dt 下 Euler 可能数值不稳定 (S 跨越 0 后吸收), 建议步数 ≥ 50/年

#include <cmath>
#include <span>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/models/diffusion/process.hpp"
#include "cpphub/pricing/analytic/cev_analytic.hpp"

namespace cpphub {
inline namespace v1 {

// ============ CEV 过程模拟方案 ============
enum class CEVScheme {
    EulerAbsorbing,  // Euler + 吸收壁 (默认, 通用 β ∈ [0, 1])
    LogEuler         // log-Euler (仅 β=1 严格, β<1 为近似)
};

// ============ CEV 过程类 ============
class CEVProcess : public StochasticProcess {
public:
    explicit CEVProcess(CEVParams p, Real S0, Real r = 0.0, Real q = 0.0,
                         CEVScheme scheme = CEVScheme::EulerAbsorbing)
        : params_(p), S0_(S0), r_(r), q_(q), scheme_(scheme) {
        validate_cev_params(p);
        if (S0 <= 0.0) {
            throw std::invalid_argument("CEVProcess: S0 must be positive");
        }
        if (p.beta < 0.0) {
            throw std::invalid_argument("CEVProcess: beta < 0 not supported (singular at S=0)");
        }
    }

    Size dimension() const override { return 1; }
    Real spot() const override { return S0_; }

    // CEV 无标准闭式 CF (Schroder 解经非中心卡方变换, 非标准 CF 形式)
    // 子类或外部函数 cev_call_price 提供解析定价
    Complex characteristic_function(Complex u, Real tau) const override {
        (void)u; (void)tau;
        return Complex{0.0, 0.0};
    }

    void generate_path(Real T, Size n_steps,
                       std::span<Real> path,
                       Philox4x64& rng) const override;

    // 一步演化 (供 MC 引擎/路径相关期权使用)
    // 输入: 当前 S, 时间步长 dt, 标准正态 Z
    // 输出: 下一步的 S
    Real evolve(Real S, Real dt, Real Z) const noexcept;

    const CEVParams& params() const { return params_; }
    Real S0() const { return S0_; }
    Real r() const { return r_; }
    Real q() const { return q_; }
    CEVScheme scheme() const { return scheme_; }

private:
    CEVParams params_;
    Real S0_;
    Real r_;
    Real q_;
    CEVScheme scheme_;
};

// ============ 一步演化 ============
inline Real CEVProcess::evolve(Real S, Real dt, Real Z) const noexcept {
    const Real sqrt_dt = std::sqrt(dt);
    const Real sigma = params_.sigma;
    const Real beta = params_.beta;
    const Real drift_rate = r_ - q_;

    if (scheme_ == CEVScheme::LogEuler) {
        // log-Euler: 仅 β=1 严格; β<1 时近似 (用瞬时波动率 σ·S^(β-1))
        if (S <= 0.0) return 0.0;  // 已被吸收
        Real vol_inst = sigma * std::pow(S, beta - 1.0);
        Real drift = (drift_rate - 0.5 * vol_inst * vol_inst) * dt;
        return S * std::exp(drift + vol_inst * sqrt_dt * Z);
    }

    // Euler + 吸收壁
    if (S <= 0.0) return 0.0;
    Real diffusion = sigma * std::pow(S, beta) * sqrt_dt * Z;
    Real drift = drift_rate * S * dt;
    Real S_next = S + drift + diffusion;
    if (S_next < 0.0) S_next = 0.0;  // 吸收壁
    return S_next;
}

// ============ 完整路径生成 ============
inline void CEVProcess::generate_path(Real T, Size n_steps,
                                       std::span<Real> path,
                                       Philox4x64& rng) const {
    if (n_steps == 0) {
        throw std::invalid_argument("CEVProcess: n_steps must be positive");
    }
    if (path.size() < n_steps + 1) {
        throw std::invalid_argument("CEVProcess: path buffer too small");
    }

    Real dt = T / static_cast<Real>(n_steps);
    path[0] = S0_;

    Real S = S0_;
    for (Size i = 0; i < n_steps; ++i) {
        // 生成标准正态 (Box-Muller)
        uint64_t r1 = rng();
        uint64_t r2 = rng();
        double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
        double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
        auto [z1, z2] = box_muller(u1, u2);
        Real Z = (i < n_steps - 1) ? z1 : z1;  // 用 z1 (z2 可留作 antithetic)
        S = evolve(S, dt, Z);
        path[i + 1] = S;
    }
}

}  // namespace v1
}  // namespace cpphub
