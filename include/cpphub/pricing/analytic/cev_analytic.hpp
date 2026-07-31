#pragma once
// SOURCE: Cox (1975) "Notes on Option Pricing I: Constant Elasticity of Variance Diffusion"
// SOURCE: Schroder (1989) "Computing the Constant Elasticity of Variance Option Pricing Formula"
// SOURCE: Davydov & Linetsky (2001) "Pricing and Hedging Path-Dependent Options under the CEV Process"
// SOURCE: Jeanblanc, Yor, Chesney (2009) "Mathematical Methods for Financial Markets" §9.3.4
//
// 模块: CEV (Constant Elasticity of Variance) 模型解析定价
//
// ==================== CEV 模型数学 ====================
//
// CEV SDE (risk-neutral):
//   dS = (r - q) S dt + σ S^β dW
//
// 参数:
//   σ (sigma) — 波动率尺度
//   β (beta)  — 弹性系数, β < 1: 吸收壁在 0; β = 1: GBM; β > 1: 可能爆炸
//
// 闭式解 (Schroder 1989, β < 1):
//   Call(S, K, T) = S e^{-qT} P_1 - K e^{-rT} P_2
//
//   其中:
//     P_1 = 1 - χ'²(c; 2/(1-β), d)         (生存概率, 即 S_T > K 的风险中性概率)
//     P_2 = χ'²(d; 2/(1-β) - 2, c)         (Schroder 对偶: 交换参数与非中心参数)
//
//   c (K 相关) 为 chi2 求值点, d (S 相关) 为非中心参数
//
//   对于 r = q (零 drift, μ→0 极限):
//     c = K^{2(1-β)} / (σ² (1-β)² T)
//     d = S^{2(1-β)} / (σ² (1-β)² T)
//
//   对于 r ≠ q (CIR 变换 Y = S^{2(1-β)}/[(1-β)²σ²], 漂移参数 b = -2(1-β)μ):
//     c = 2(r-q) K^{2(1-β)} / (σ²(1-β)(e^{2(1-β)(r-q)T} - 1))
//     d = 2(r-q) S^{2(1-β)} e^{2(1-β)(r-q)T} / (σ²(1-β)(e^{2(1-β)(r-q)T} - 1))
//
//   χ'²(z; k, λ) = 非中心卡方 CDF, 自由度 k, 非中心参数 λ, 在 z 处求值
//
// 极限:
//   β = 1: 退化为 GBM, 用 Black-Scholes 公式 (公式在 β=1 处奇异)
//   β = 0: dS = (r-q)S dt + σ dW (r=q 时为 Bachelier)
//   β → -∞: 趋向于对数正态的反向
//
// 注: β > 1 (爆炸情形) 此版本未实现, 需用对偶非中心卡方 (Schroder 1989)

#include <cmath>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/math.hpp"

namespace cpphub {
inline namespace v1 {

// ============ CEV 参数结构 ============
struct CEVParams {
    Real sigma;  // 波动率尺度, > 0
    Real beta;   // 弹性系数, < 1 (本实现), β=1 退化为 GBM
};

// ============ 参数验证 ============
inline void validate_cev_params(const CEVParams& p) {
    if (p.sigma <= 0.0) {
        throw std::invalid_argument("CEV: sigma must be positive");
    }
    if (p.beta > 1.0) {
        throw std::invalid_argument("CEV: beta > 1 not supported in this implementation");
    }
}

// ============ CEV 欧式 Call 定价 ============
inline Real cev_call_price(Real S, Real K, Real T, Real r, Real q, const CEVParams& p) {
    validate_cev_params(p);
    if (S <= 0.0) throw std::invalid_argument("CEV: spot S must be positive");
    if (K <= 0.0) throw std::invalid_argument("CEV: strike K must be positive");
    if (T < 0.0) throw std::invalid_argument("CEV: T must be non-negative");

    if (T == 0.0) {
        return std::max(S - K, 0.0);
    }

    const Real beta = p.beta;
    const Real sigma = p.sigma;

    // β = 1: GBM 极限, 用 Black-Scholes
    if (std::abs(beta - 1.0) < 1e-10) {
        // d1, d2
        Real sqrt_T = std::sqrt(T);
        Real vol = sigma;  // β=1 时 σ 即 BS 波动率
        Real d1 = (std::log(S / K) + (r - q + 0.5 * vol * vol) * T) / (vol * sqrt_T);
        Real d2 = d1 - vol * sqrt_T;
        return S * std::exp(-q * T) * normal_cdf(d1) - K * std::exp(-r * T) * normal_cdf(d2);
    }

    // β < 1: 非中心卡方公式
    const Real one_minus_beta = 1.0 - beta;
    const Real two_omb = 2.0 * one_minus_beta;

    Real c, d;  // 非中心卡方参数 (c 关联 K, d 关联 S)
    Real drift = r - q;

    if (std::abs(drift) < 1e-10) {
        // r = q (μ→0 极限): c, d 退化为 K^{2(1-β)}, S^{2(1-β)} 归一化
        // 由一般分支取 μ→0: 2μ/(σ²(1-β)(e^{2(1-β)μT}-1)) → 1/(σ²(1-β)²T)
        Real denom = sigma * sigma * one_minus_beta * one_minus_beta * T;
        c = std::pow(K, two_omb) / denom;
        d = std::pow(S, two_omb) / denom;
    } else {
        // r ≠ q: 含 drift 调整
        // CIR 变换 Y = S^{2(1-β)}/[(1-β)²σ²] 的漂移参数 b = -2(1-β)μ
        // 故指数为 2(1-β)μT (非 2μT)
        Real exp_arg = 2.0 * one_minus_beta * drift * T;
        Real exp_2dt = std::exp(exp_arg);
        Real denom = sigma * sigma * one_minus_beta * std::expm1(exp_arg);
        // denom = σ²(1-β)(e^{2(1-β)μT} - 1)
        c = 2.0 * drift * std::pow(K, two_omb) / denom;
        d = 2.0 * drift * std::pow(S, two_omb) * exp_2dt / denom;
    }

    // 自由度 (Schroder 1989): P_1 与 P_2 的 df 相差 2
    // ν_1 = 2/(1-β), ν_2 = 2/(1-β) - 2 = 2β/(1-β)
    Real nu_1 = 2.0 / one_minus_beta;         // P_1 用
    Real nu_2 = 2.0 / one_minus_beta - 2.0;   // P_2 用

    // P_1 = 1 - χ'²(c; nu_1, d)  (生存概率, S_T > K)
    //   c (K 相关) 为 chi2 参数, d (S 相关) 为非中心参数
    // P_2 = χ'²(d; nu_2, c)       (Schroder 对偶: 交换参数与非中心参数)
    //   d (S 相关) 为 chi2 参数, c (K 相关) 为非中心参数
    Real chi2_1 = noncentral_chi2_cdf(c, nu_1, d);
    Real chi2_2 = noncentral_chi2_cdf(d, nu_2, c);

    Real P1 = 1.0 - chi2_1;
    Real P2 = chi2_2;

    // Call = S e^{-qT} P_1 - K e^{-rT} P_2
    return S * std::exp(-q * T) * P1 - K * std::exp(-r * T) * P2;
}

// ============ CEV 欧式 Put 定价 (Call-Put Parity) ============
// CEV Call-Put Parity: C - P = S e^{-qT} - K e^{-rT}
// (仅在 β < 1 且吸收壁在 0 时严格成立; β > 1 时 parity 可能失效)
inline Real cev_put_price(Real S, Real K, Real T, Real r, Real q, const CEVParams& p) {
    Real C = cev_call_price(S, K, T, r, q, p);
    return C - S * std::exp(-q * T) + K * std::exp(-r * T);
}

}  // namespace v1
}  // namespace cpphub
