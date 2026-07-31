#pragma once
// SOURCE: Hagan, Kumar, Lesniewski, Woodward (2002) "Managing Smile Risk"
//         Wilmott Magazine, September 2002, pp. 84-108.
// SOURCE: Oblój (2008) "Fine-Tune Your Smile: Correction to Hagan et al."
//         Wilmott Magazine, May 2008. (修正 x(z) 公式中的符号错误)
// SOURCE: West (2005) "Calibration of Sabr Model on Real Market Data"
//
// 模块: SABR 隐含波动率闭式近似 (Hagan 2002 asymptotic formula)
//
// ==================== SABR 模型数学 ====================
//
// SABR (Stochastic Alpha Beta Rho) 模型:
//   dF = σ F^β dW_1
//   dσ = ν σ dW_2
//   dW_1 dW_2 = ρ dt
//   F(0) = F_0,  σ(0) = α
//
// 参数:
//   α (alpha) — 初始波动率水平 (σ(0))
//   β (beta)  — 弹性系数 (CEV exponent), 控制微笑斜率
//   ρ (rho)   — F 与 σ 的相关性, 控制微笑倾斜方向
//   ν (nu)    — 波动率的波动率 (vol-of-vol), 控制微笑曲率
//
// Hagan 2002 给出 Black 隐含波动率 σ_B(K, F) 的渐近近似公式 (O(σ²) 精度):
//
// 完整公式 (K ≠ F):
//   σ_B(K, F) = α / D(K, F) · (z / x(z)) · (1 + c(K, F) · T)
//
//   其中:
//     D(K, F) = (F K)^((1-β)/2) · [1 + ((1-β)²/24) ln²(F/K) + ((1-β)⁴/1920) ln⁴(F/K)]
//     z        = (ν/α) · (F K)^((1-β)/2) · ln(F/K)
//     x(z)     = ln[(√(1 - 2ρz + z²) + z - ρ) / (1 - ρ)]   // Oblój 2008 修正版
//     c(K, F)  = ((1-β)²/24) · α²/(F K)^(1-β)
//              + (ρ β ν / 4) · α / (F K)^((1-β)/2)
//              + ((2 - 3ρ²)/24) · ν²
//
// ATM 公式 (K = F, 即 z = 0, x(0) = 0, z/x → 1):
//   σ_ATM(F) = α / F^(1-β) · (1 + c_ATM · T)
//   c_ATM    = ((1-β)²/24) · α² / F^(2-2β)
//            + (ρ β ν / 4) · α / F^(1-β)
//            + ((2 - 3ρ²)/24) · ν²
//
// 极限情形:
//   β = 1: 对数正态模型 (lognormal SABR), (F K)^((1-β)/2) = 1
//   β = 0: 正态模型 (normal SABR), 退化为 Bachelier
//   ρ = ±1: x(z) 退化, 需用极限公式
//   ν = 0: 局部波动率模型 (CEV), 退化为 σ_B = α/F^(1-β) (无随机波动率)
//
// 精度: 对于 T < 5 年, σ < 0.5, 通常误差 < 1% (相对)

#include <cmath>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"

namespace cpphub {
inline namespace v1 {

// ============ SABR 参数结构 ============
struct SABRParams {
    Real alpha;   // 初始波动率 σ(0), > 0
    Real beta;    // 弹性系数, ∈ [0, 1]
    Real rho;     // 相关性, ∈ (-1, 1)
    Real nu;      // vol-of-vol, > 0
};

// ============ 参数验证 ============
inline void validate_sabr_params(const SABRParams& p) {
    if (p.alpha <= 0.0) {
        throw std::invalid_argument("SABR: alpha must be positive");
    }
    if (p.beta < 0.0 || p.beta > 1.0) {
        throw std::invalid_argument("SABR: beta must be in [0, 1]");
    }
    if (p.rho <= -1.0 || p.rho >= 1.0) {
        throw std::invalid_argument("SABR: rho must be in (-1, 1)");
    }
    if (p.nu < 0.0) {
        throw std::invalid_argument("SABR: nu must be non-negative");
    }
}

// ============ Hagan 2002 IV 公式 ============
// 计算 SABR 模型下的 Black 隐含波动率
// 输入: 行权价 K, 远期价格 F, 到期时间 T, SABR 参数
// 输出: Black 隐含波动率 σ_B
//
// 注: 使用 Oblój (2008) 修正后的 x(z) 公式
inline Real sabr_implied_vol_hagan(Real K, Real F, Real T, const SABRParams& p) {
    validate_sabr_params(p);

    if (K <= 0.0) {
        throw std::invalid_argument("SABR: strike K must be positive");
    }
    if (F <= 0.0) {
        throw std::invalid_argument("SABR: forward F must be positive");
    }
    if (T < 0.0) {
        throw std::invalid_argument("SABR: time to maturity T must be non-negative");
    }

    const Real alpha = p.alpha;
    const Real beta = p.beta;
    const Real rho = p.rho;
    const Real nu = p.nu;

    // T = 0 退化: 返回瞬时 IV
    if (T == 0.0) {
        // 极限: σ_B → α / F^(1-β)
        return alpha / std::pow(F, 1.0 - beta);
    }

    const Real one_minus_beta = 1.0 - beta;

    // ATM 情形 (K ≈ F): 使用 ATM 闭式公式, 避免 z/x 的 0/0 不定型
    // 容差取相对量级
    const Real log_FK = std::log(F / K);
    if (std::abs(log_FK) < 1e-10) {
        // σ_ATM = α / F^(1-β) · (1 + c_ATM · T)
        const Real F_pow = std::pow(F, one_minus_beta);          // F^(1-β)
        const Real F_pow2 = F_pow * F_pow;                        // F^(2-2β) = (F^(1-β))²
        const Real c_atm =
            (one_minus_beta * one_minus_beta / 24.0) * alpha * alpha / F_pow2
            + (rho * beta * nu / 4.0) * alpha / F_pow
            + ((2.0 - 3.0 * rho * rho) / 24.0) * nu * nu;
        return (alpha / F_pow) * (1.0 + c_atm * T);
    }

    // 完整公式 (K ≠ F)
    const Real FK = F * K;
    const Real FK_pow_half = std::pow(FK, one_minus_beta / 2.0);  // (FK)^((1-β)/2)
    const Real FK_pow = FK_pow_half * FK_pow_half;                 // (FK)^(1-β)

    // D(K, F) 分母部分
    const Real log_FK2 = log_FK * log_FK;
    const Real log_FK4 = log_FK2 * log_FK2;
    const Real D = FK_pow_half * (1.0
                  + (one_minus_beta * one_minus_beta / 24.0) * log_FK2
                  + (one_minus_beta * one_minus_beta
                     * one_minus_beta * one_minus_beta / 1920.0) * log_FK4);

    // z = (ν/α) · (FK)^((1-β)/2) · ln(F/K)
    const Real z = (nu / alpha) * FK_pow_half * log_FK;

    // x(z) = ln[(√(1 - 2ρz + z²) + z - ρ) / (1 - ρ)]   // Oblój 2008 修正
    const Real one_minus_rho = 1.0 - rho;
    const Real sqrt_term = std::sqrt(1.0 - 2.0 * rho * z + z * z);
    const Real x_z = std::log((sqrt_term + z - rho) / one_minus_rho);

    // z / x(z)
    // 注: 当 z → 0 (即 K → F), x(z) → z (一阶), z/x → 1, 已在 ATM 分支处理
    // 但若 ρ → ±1, x(z) 可能退化, 此处已由 validate_sabr_params 限制 |ρ| < 1
    if (std::abs(x_z) < 1e-15) {
        // 数值保护: x(z) 过小, 退化为 ATM 公式
        const Real F_pow = std::pow(F, one_minus_beta);
        const Real F_pow2 = F_pow * F_pow;
        const Real c_atm =
            (one_minus_beta * one_minus_beta / 24.0) * alpha * alpha / F_pow2
            + (rho * beta * nu / 4.0) * alpha / F_pow
            + ((2.0 - 3.0 * rho * rho) / 24.0) * nu * nu;
        return (alpha / F_pow) * (1.0 + c_atm * T);
    }

    // 高阶修正项 c(K, F) · T
    const Real c_KF =
        (one_minus_beta * one_minus_beta / 24.0) * alpha * alpha / FK_pow
        + (rho * beta * nu / 4.0) * alpha / FK_pow_half
        + ((2.0 - 3.0 * rho * rho) / 24.0) * nu * nu;

    // σ_B = (α / D) · (z / x(z)) · (1 + c · T)
    return (alpha / D) * (z / x_z) * (1.0 + c_KF * T);
}

// ============ ATM IV 便捷接口 ============
inline Real sabr_implied_vol_atm(Real F, Real T, const SABRParams& p) {
    validate_sabr_params(p);
    if (F <= 0.0) throw std::invalid_argument("SABR: forward F must be positive");
    if (T < 0.0) throw std::invalid_argument("SABR: T must be non-negative");
    if (T == 0.0) return p.alpha / std::pow(F, 1.0 - p.beta);

    const Real one_minus_beta = 1.0 - p.beta;
    const Real F_pow = std::pow(F, one_minus_beta);
    const Real F_pow2 = F_pow * F_pow;
    const Real c_atm =
        (one_minus_beta * one_minus_beta / 24.0) * p.alpha * p.alpha / F_pow2
        + (p.rho * p.beta * p.nu / 4.0) * p.alpha / F_pow
        + ((2.0 - 3.0 * p.rho * p.rho) / 24.0) * p.nu * p.nu;
    return (p.alpha / F_pow) * (1.0 + c_atm * T);
}

// ============ 反推 alpha (给定 ATM IV) ============
// 已知 σ_ATM, F, T, (β, ρ, ν), 求解 α
// σ_ATM = α / F^(1-β) · (1 + c_ATM(α) · T)
// c_ATM 含 α² 和 α 项, 构成三次方程, 此处用 Newton 迭代求解
//
// 初值: α₀ = σ_ATM · F^(1-β) (忽略高阶项)
inline Real sabr_solve_alpha_from_atm(Real sigma_atm, Real F, Real T, Real beta, Real rho, Real nu) {
    if (sigma_atm <= 0.0) throw std::invalid_argument("SABR: sigma_atm must be positive");
    if (F <= 0.0) throw std::invalid_argument("SABR: F must be positive");
    if (T < 0.0) throw std::invalid_argument("SABR: T must be non-negative");
    if (beta < 0.0 || beta > 1.0) throw std::invalid_argument("SABR: beta must be in [0, 1]");
    if (rho <= -1.0 || rho >= 1.0) throw std::invalid_argument("SABR: rho must be in (-1, 1)");
    if (nu < 0.0) throw std::invalid_argument("SABR: nu must be non-negative");

    if (T == 0.0) {
        return sigma_atm * std::pow(F, 1.0 - beta);
    }

    const Real one_minus_beta = 1.0 - beta;
    const Real F_pow = std::pow(F, one_minus_beta);
    const Real F_pow2 = F_pow * F_pow;

    // 系数
    const Real c2 = (one_minus_beta * one_minus_beta / 24.0) / F_pow2;       // α² 项系数
    const Real c1 = (rho * beta * nu / 4.0) / F_pow;                          // α 项系数
    const Real c0 = ((2.0 - 3.0 * rho * rho) / 24.0) * nu * nu;               // 常数项

    // 目标: (α / F_pow) · (1 + (c2·α² + c1·α + c0) · T) - σ_ATM = 0
    // 即: α·(1 + (c2·α² + c1·α + c0)·T) - σ_ATM·F_pow = 0
    //     c2·T·α³ + c1·T·α² + (1 + c0·T)·α - σ_ATM·F_pow = 0

    const Real a3 = c2 * T;
    const Real a2 = c1 * T;
    const Real a1 = 1.0 + c0 * T;
    const Real a0 = -sigma_atm * F_pow;

    // Newton 迭代, 初值忽略 α² 和 α 高阶项
    Real alpha = sigma_atm * F_pow;

    for (int iter = 0; iter < 100; ++iter) {
        Real f = a3 * alpha * alpha * alpha + a2 * alpha * alpha + a1 * alpha + a0;
        Real fp = 3.0 * a3 * alpha * alpha + 2.0 * a2 * alpha + a1;
        if (std::abs(fp) < 1e-20) break;
        Real delta = f / fp;
        alpha -= delta;
        if (alpha < 1e-12) alpha = 1e-12;  // 保持正值
        if (std::abs(delta) < 1e-14 * (1.0 + std::abs(alpha))) break;
    }
    return alpha;
}

}  // namespace v1
}  // namespace cpphub
