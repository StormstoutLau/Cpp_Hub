#pragma once
// SOURCE: Madan, Carr & Chang (1998) "The Variance Gamma Process and Option Pricing"
// SOURCE: Madan & Seneta (1990) "The Variance Gamma (V.G.) Model for Share Market Returns"
// SOURCE: Cont & Tankov (2004) "Financial Modelling with Jump Processes" §4.3
//
// 模块: Variance Gamma (VG) 模型解析层 (特征函数 + 累积量 + omega)
//
// ==================== VG 模型数学 ====================
//
// VG 过程 X(t) 是一个纯跳跃 Levy 过程, 由 Gamma 时变的 Brownian 运动构造:
//   X(t) = θ·G(t) + σ·W(G(t))
//   G(t) ~ Gamma(t/ν, 1/ν)  (shape=t/ν, rate=1/ν, 即 scale=ν)
//   W 为标准 Brownian 运动, 与 G 独立
//
// 风险中性价格过程:
//   S_T = S_0 · exp((r - q + ω)T + X_T)
//   ω = (1/ν) · ln(1 - θν - σ²ν/2)  (鞅修正, 保证 E[S_T] = S_0 e^{(r-q)T})
//
// 特征函数 (ln S_T):
//   φ(u, τ) = exp(iu(ln S_0 + (r-q+ω)τ)) · (1 - iuθν + σ²νu²/2)^{-τ/ν}
//
// 累积量 (对 X_T):
//   E[X_T] = θT
//   Var[X_T] = (σ² + θ²ν)T
//   偏度 = (2θ³ν² + 3σ²θν)T / Var^{3/2}
//   峰度 = 3 + (3σ⁴ν + 12θ²σ²ν² + 2θ⁴ν³)T / Var²
//
// 边界:
//   ν → 0: VG 退化为 Black-Scholes (Gamma 时变 → 确定性时变)
//   θ = 0: 对称 VG (无偏度, 仅厚尾)
//   σ = 0: 纯跳跃 VG (X_T = θ·G_T)

#include <complex>
#include <cmath>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/pricing/fourier/characteristic_functions.hpp"

namespace cpphub {
inline namespace v1 {

// ============ VG 参数结构 ============
struct VGParams {
    Real sigma;  // Brownian 波动率 > 0
    Real nu;     // Gamma 过程方差率 > 0 (G(t) ~ Gamma(t/ν, ν))
    Real theta;  // drift 参数 (VG 偏度来源)
};

// ============ 参数验证 ============
inline void validate_vg_params(const VGParams& p) {
    if (p.sigma <= 0.0) throw std::invalid_argument("VG: sigma must be positive");
    if (p.nu <= 0.0) throw std::invalid_argument("VG: nu must be positive");
    // Feller 条件: 1 - θν - σ²ν/2 > 0 (保证 ω 有定义)
    Real condition = 1.0 - p.theta * p.nu - 0.5 * p.sigma * p.sigma * p.nu;
    if (condition <= 0.0) {
        throw std::invalid_argument("VG: Feller condition violated (1 - theta*nu - sigma^2*nu/2 > 0)");
    }
}

// ============ 鞅修正 ω = (1/ν) ln(1 - θν - σ²ν/2) ============
inline Real vg_omega(Real sigma, Real nu, Real theta) {
    Real condition = 1.0 - theta * nu - 0.5 * sigma * sigma * nu;
    if (condition <= 0.0) {
        throw std::invalid_argument("VG omega: Feller condition violated");
    }
    return std::log(condition) / nu;
}

// ============ VG 特征函数 ============
// φ(u, τ) = exp(iu(ln S_0 + (r-q+ω)τ)) · (1 - iuθν + σ²νu²/2)^{-τ/ν}
inline Complex vg_characteristic_function(Complex u, Real tau,
                                            Real S0, Real r, Real q,
                                            const VGParams& p) {
    validate_vg_params(p);
    if (S0 <= 0.0) throw std::invalid_argument("VG: S0 must be positive");
    if (tau < 0.0) throw std::invalid_argument("VG: tau must be non-negative");
    if (tau == 0.0) return Complex(1.0, 0.0);

    if (std::abs(u) < Real(1e-15)) return Complex(1.0, 0.0);

    Real omega = vg_omega(p.sigma, p.nu, p.theta);
    Real ln_S0 = std::log(S0);
    Real drift = (r - q + omega) * tau;

    const Complex i(0.0, 1.0);
    Complex iu = i * u;

    // 相位项: exp(iu(ln S_0 + drift))
    Complex phase = std::exp(iu * (ln_S0 + drift));

    // 基础项: (1 - iuθν + σ²νu²/2)^{-τ/ν}
    // 注意: u² 为复数平方
    Complex u_sq = u * u;
    Complex base = Complex(1.0, 0.0) - iu * (p.theta * p.nu)
                 + Complex(0.5 * p.sigma * p.sigma * p.nu * u_sq.real(),
                           0.5 * p.sigma * p.sigma * p.nu * u_sq.imag());
    Real exponent = -tau / p.nu;
    return phase * std::pow(base, exponent);
}

// ============ VG CF 工厂 (CharFn 闭包, 用于 COSEngine) ============
// 注意: 与 characteristic_functions.hpp 中的 make_vg_cf 功能一致,
//       此处提供独立接口便于 VGProcess 类内部使用
inline CharFn make_vg_cf_direct(Real S0, Real r, Real q,
                                  const VGParams& p, Real T) {
    if (S0 <= 0.0) throw std::invalid_argument("vg_cf: S0 must be positive");
    if (T <= 0.0) throw std::invalid_argument("vg_cf: T must be positive");
    VGParams params = p;
    return [S0, T, r, q, params](Complex u) -> Complex {
        return vg_characteristic_function(u, T, S0, r, q, params);
    };
}

// ============ VG 累积量 (用于 COS 截断和测试验证) ============
inline Real vg_cumulant_mean(Real tau, const VGParams& p) {
    // E[X_T] = θT
    return p.theta * tau;
}

inline Real vg_cumulant_variance(Real tau, const VGParams& p) {
    // Var[X_T] = (σ² + θ²ν)T
    return (p.sigma * p.sigma + p.theta * p.theta * p.nu) * tau;
}

inline Real vg_cumulant_skewness(Real tau, const VGParams& p) {
    // 三阶累积量 κ_3 = (2θ³ν² + 3σ²θν)T
    Real kappa3 = (2.0 * p.theta * p.theta * p.theta * p.nu * p.nu
                   + 3.0 * p.sigma * p.sigma * p.theta * p.nu) * tau;
    Real var = vg_cumulant_variance(tau, p);
    return kappa3 / std::pow(var, 1.5);
}

inline Real vg_cumulant_kurtosis_excess(Real tau, const VGParams& p) {
    // 四阶累积量 κ_4 = (3σ⁴ν + 12θ²σ²ν² + 2θ⁴ν³)T
    Real kappa4 = (3.0 * p.sigma * p.sigma * p.sigma * p.sigma * p.nu
                   + 12.0 * p.theta * p.theta * p.sigma * p.sigma * p.nu * p.nu
                   + 2.0 * p.theta * p.theta * p.theta * p.theta * p.nu * p.nu * p.nu) * tau;
    Real var = vg_cumulant_variance(tau, p);
    return kappa4 / (var * var);  // 超额峰度
}

}  // namespace v1
}  // namespace cpphub
