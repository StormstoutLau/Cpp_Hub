#pragma once
// SOURCE: Bates (1996) "Jumps and Stochastic Volatility: Exchange Rate Processes
//         Implicit in Deutsche Mark Options"
// SOURCE: Merton (1976) "Option Pricing When Underlying Stock Returns Are Discontinuous"
// SOURCE: Bakshi, Cao, Chen (1997) "Empirical Performance of Alternative Option Pricing Models"
//
// 模块: Bates 模型特征函数 (Heston 随机波动率 + Merton 对数正态跳跃)
//
// ==================== Bates 模型数学 ====================
//
// 风险中性 SDE:
//   dS/S = (r - q - λm) dt + sqrt(v) dW^S + dJ
//   dv   = κ(θ - v) dt + σ sqrt(v) dW^v
//   dW^S dW^v = ρ dt
//   J = compound Poisson(λ), 跳跃幅度 J_k ~ LogNormal(μ_J, σ_J²)
//   m = E[J_k - 1] = exp(μ_J + σ_J²/2) - 1  (跳跃补偿, 保证风险中性)
//
// 特征函数 (Bates 1996):
//   φ_Bates(u, τ) = φ_Heston(u, τ; r̃) × φ_Jump(u, τ)
//
// 其中:
//   r̃ = r - λm  (调整 drift 以补偿跳跃均值)
//   φ_Heston(u, τ; r̃) 为 Heston CF (用 r̃ 替代 r)
//   φ_Jump(u, τ) = exp(λτ(exp(iμ_J u - σ_J²u²/2) - 1))  (Merton 跳跃 CF)
//
// 边界:
//   λ = 0: 退化为 Heston (φ_Jump = 1)
//   σ_J = 0: 跳跃幅度恒定 exp(μ_J), 仍为 Merton 形式
//
// 注: Heston CF 的 Little Trap 修正 (Albrecher 2007) 通过 heston_cf.hpp 继承

#include <complex>
#include <cmath>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include "cpphub/pricing/fourier/characteristic_functions.hpp"

namespace cpphub {
inline namespace v1 {

// ============ Bates CF 参数结构 ============
struct BatesCFParams {
    // Heston 部分
    Real v0;       // 初始方差 > 0
    Real kappa;    // 均值回归速度 > 0
    Real theta;    // 长期方差 > 0
    Real sigma;    // 波动率波动率 > 0
    Real rho;      // 相关性 [-1, 1]
    // Merton 跳跃部分
    Real lambda;   // 跳跃强度 >= 0
    Real mu_J;     // log(J) 的均值
    Real sigma_J;  // log(J) 的标准差 > 0
    // 市场参数
    Real r;        // 无风险利率
    Real q;        // 股息率
};

// ============ 参数验证 ============
inline void validate_bates_cf_params(const BatesCFParams& p) {
    if (p.v0 <= 0.0) throw std::invalid_argument("Bates CF: v0 must be positive");
    if (p.kappa <= 0.0) throw std::invalid_argument("Bates CF: kappa must be positive");
    if (p.theta <= 0.0) throw std::invalid_argument("Bates CF: theta must be positive");
    if (p.sigma <= 0.0) throw std::invalid_argument("Bates CF: sigma must be positive");
    if (p.rho < -1.0 || p.rho > 1.0) throw std::invalid_argument("Bates CF: rho must be in [-1, 1]");
    if (p.lambda < 0.0) throw std::invalid_argument("Bates CF: lambda must be non-negative");
    if (p.sigma_J <= 0.0) throw std::invalid_argument("Bates CF: sigma_J must be positive");
}

// ============ 跳跃补偿 m = E[J-1] = exp(μ_J + σ_J²/2) - 1 ============
inline Real bates_jump_compensation(Real mu_J, Real sigma_J) {
    return std::exp(mu_J + 0.5 * sigma_J * sigma_J) - 1.0;
}

// ============ Merton 跳跃特征函数 ============
// φ_J(u, τ) = exp(λτ(exp(iμ_J u - σ_J²u²/2) - 1))
inline Complex merton_jump_cf(Complex u, Real tau,
                                Real lambda, Real mu_J, Real sigma_J) {
    if (lambda == 0.0) return Complex(1.0, 0.0);
    const Complex i(0.0, 1.0);
    Complex iu = i * u;
    // exponent = iμ_J u - σ_J²u²/2  (注意 u² 为复数平方)
    Complex exponent = iu * mu_J - 0.5 * sigma_J * sigma_J * (u * u);
    return std::exp(lambda * tau * (std::exp(exponent) - 1.0));
}

// ============ Bates 特征函数 ============
inline Complex bates_characteristic_function(Complex u, Real tau,
                                               Real S0, const BatesCFParams& p) {
    validate_bates_cf_params(p);
    if (S0 <= 0.0) throw std::invalid_argument("Bates CF: S0 must be positive");
    if (tau < 0.0) throw std::invalid_argument("Bates CF: tau must be non-negative");
    if (tau == 0.0) return Complex(1.0, 0.0);

    if (std::abs(u) < Real(1e-15)) return Complex(1.0, 0.0);

    // 跳跃补偿: 调整 drift
    Real m = bates_jump_compensation(p.mu_J, p.sigma_J);
    Real adjusted_r = p.r - p.lambda * m;

    // Heston CF (用调整后的 r)
    HestonCFParams hp{p.v0, p.kappa, p.theta, p.sigma, p.rho, adjusted_r, p.q};
    Complex heston_cf = heston_characteristic_function(u, tau, S0, hp);

    // λ=0 时跳跃 CF = 1, 直接返回 Heston
    if (p.lambda == 0.0) return heston_cf;

    // Merton 跳跃 CF
    Complex jump_cf = merton_jump_cf(u, tau, p.lambda, p.mu_J, p.sigma_J);

    return heston_cf * jump_cf;
}

// ============ Bates CF 工厂 (CharFn 闭包, 用于 COS/FFT 引擎) ============
inline CharFn make_bates_cf(Real S0, Real r, Real q,
                              const BatesCFParams& p, Real T) {
    if (S0 <= 0.0) throw std::invalid_argument("bates_cf: S0 must be positive");
    if (T <= 0.0) throw std::invalid_argument("bates_cf: T must be positive");
    BatesCFParams params = p;
    params.r = r;
    params.q = q;
    return [S0, T, params](Complex u) -> Complex {
        return bates_characteristic_function(u, T, S0, params);
    };
}

}  // namespace v1
}  // namespace cpphub
