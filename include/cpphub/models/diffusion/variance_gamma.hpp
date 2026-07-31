#pragma once
// SOURCE: Madan, Carr & Chang (1998) "The Variance Gamma Process and Option Pricing"
// SOURCE: Cont & Tankov (2004) "Financial Modelling with Jump Processes" §4.3
//
// 模块: Variance Gamma (VG) 随机过程模拟
//
// ==================== 过程模拟 ====================
//
// VG 是纯跳跃 Levy 过程, 有独立平稳增量. 模拟方式:
//
// 方法: 按时间步 Δt 采样独立增量
//   ΔG ~ Gamma(Δt/ν, ν)  (shape=Δt/ν, scale=ν, 均值=Δt, 方差=Δt·ν)
//   ΔW ~ Normal(0, ΔG)    (Brownian 在随机时间 ΔG 上)
//   ΔX = θ·ΔG + σ·ΔW      (VG 增量)
//
// 价格更新:
//   S_{t+Δt} = S_t · exp((r - q + ω)Δt + ΔX)
//   ω = (1/ν) ln(1 - θν - σ²ν/2)  (鞅修正)
//
// 注: Gamma 采样使用 std::gamma_distribution (Marsaglia-Tsang 算法),
//     Philox4x64 满足 UniformRandomBitGenerator 概念可直接传入.

#include <cmath>
#include <complex>
#include <random>
#include <span>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/models/diffusion/process.hpp"
#include "cpphub/pricing/analytic/vg_analytic.hpp"

namespace cpphub {
inline namespace v1 {

// ============ VG 过程类 ============
class VarianceGammaProcess : public StochasticProcess {
public:
    explicit VarianceGammaProcess(VGParams p, Real S0, Real r = 0.0, Real q = 0.0)
        : params_(p), S0_(S0), r_(r), q_(q),
          omega_(vg_omega(p.sigma, p.nu, p.theta)) {
        validate_vg_params(p);
        if (S0 <= 0.0) throw std::invalid_argument("VGProcess: S0 must be positive");
    }

    Size dimension() const override { return 1; }
    Real spot() const override { return S0_; }

    Complex characteristic_function(Complex u, Real tau) const override;
    void generate_path(Real T, Size n_steps,
                       std::span<Real> path, Philox4x64& rng) const override;

    // 单步增量 (供测试验证)
    // 输入: dt, 已采样的 Gamma 增量 dG, 正态 Z
    // 输出: VG 增量 ΔX = θ·ΔG + σ·sqrt(ΔG)·Z
    Real vg_increment(Real dt, Real dG, Real Z) const noexcept;

    // 采样 Gamma 增量 ΔG ~ Gamma(dt/ν, ν)
    Real sample_gamma_increment(Real dt, Philox4x64& rng) const;

    const VGParams& params() const { return params_; }
    Real S0() const { return S0_; }
    Real r() const { return r_; }
    Real q() const { return q_; }
    Real omega() const { return omega_; }

private:
    VGParams params_;
    Real S0_;
    Real r_;
    Real q_;
    Real omega_;  // 鞅修正项
};

inline Real VarianceGammaProcess::sample_gamma_increment(Real dt, Philox4x64& rng) const {
    // Gamma(shape=dt/ν, scale=ν) 增量
    Real shape = dt / params_.nu;
    Real scale = params_.nu;
    std::gamma_distribution<Real> dist(shape, scale);
    return dist(rng);
}

inline Real VarianceGammaProcess::vg_increment(Real dt, Real dG, Real Z) const noexcept {
    (void)dt;
    // ΔX = θ·ΔG + σ·sqrt(ΔG)·Z
    return params_.theta * dG + params_.sigma * std::sqrt(dG) * Z;
}

inline void VarianceGammaProcess::generate_path(Real T, Size n_steps,
                                                  std::span<Real> path,
                                                  Philox4x64& rng) const {
    Real dt = T / static_cast<Real>(n_steps);
    path[0] = S0_;

    Real S = S0_;
    Real drift_rate = r_ - q_ + omega_;

    for (Size i = 0; i < n_steps; ++i) {
        // 采样 Gamma 增量
        Real dG = sample_gamma_increment(dt, rng);

        // 采样标准正态 (用于 Brownian 部分)
        uint64_t r1 = rng();
        uint64_t r2 = rng();
        double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
        double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
        auto [z1, z2] = box_muller(u1, u2);
        (void)z2;

        // VG 增量
        Real dX = vg_increment(dt, dG, z1);

        // 价格更新
        Real drift = drift_rate * dt;
        S = S * std::exp(drift + dX);
        path[i + 1] = S;
    }
}

inline Complex VarianceGammaProcess::characteristic_function(Complex u, Real tau) const {
    return vg_characteristic_function(u, tau, S0_, r_, q_, params_);
}

}  // namespace v1
}  // namespace cpphub
