// =============================================================================
// garch_distribution.hpp - GARCH 分布族 (Normal / Student-t / GED) 与似然函数
//
// Phase 7B v1.6 M1 (PHASE7B_FINANCIAL_TS_SPEC.md §2.0.1)
//
// 教材锚点: Bollerslev 1986 (Normal) / Bollerslev 1987 RES (t) / Nelson 1991 (GED)
// 对照库: Python arch `arch/univariate/distribution.py`
//
// 幻觉点防护 (spec §6.1):
//   G3: Normal 似然必须含 -0.5·log(2π) 常数项
//   G14: t 分布自由度 ν 与 (ω,α,β) 联合 QMLE (由 estimate_* 承担, 此处提供 ℓₜ(ν))
//   G-GED1: GED 指数为 ν (非 2ν); 缩放常数 c 含 2^(-1/ν) 因子;
//           常数项为 log(ν) - log(c) - log(Γ(1/ν)) - (1+1/ν)·log(2)
//
// 仅依赖 core/ (ADR-017 依赖图: garch_distribution.hpp 独立)
// =============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>

#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace garch {

// GARCH 分布族 (G-ADR5: t 分布与 (ω,α,β) 联合 QMLE)
enum class GarchDist {
    Normal,   ///< 标准正态: ℓ = -0.5·[log(2π) + log(hₜ) + ε²ₜ/hₜ]
    StudentT, ///< Student-t (ν 自由度): 含 Γ((ν+1)/2)/Γ(ν/2) 项
    GED       ///< 广义误差分布 (Nelson 1991): 含 Γ(1/ν)/Γ(3/ν) 项
};

// 分布参数个数 (AIC/BIC 的 k 计数用)
inline Size num_dist_params(GarchDist dist) {
    return dist == GarchDist::Normal ? 0u : 1u;
}

/// @brief 单观测对数似然
///
/// 幻觉点 G3: Normal 必须含 -0.5·log(2π) 常数项。
/// StudentT (Bollerslev 1987, 标准化到单位方差):
///   ℓₜ = log(Γ((ν+1)/2)) - log(Γ(ν/2)) - 0.5·log(ν-2) - 0.5·log(π)
///        - 0.5·log(hₜ) - ((ν+1)/2)·log(1 + ε²ₜ/((ν-2)·hₜ))
/// GED (Nelson 1991, arch distribution.py 实测):
///   c = 2^(-1/ν) · sqrt(Γ(1/ν)/Γ(3/ν))
///   ℓₜ = log(ν) - log(c) - log(Γ(1/ν)) - (1+1/ν)·log(2)
///        - 0.5·log(hₜ) - 0.5·|εₜ/(√hₜ·c)|^ν
///
/// @throws std::invalid_argument 若 variance <= 0, 或 (t/GED 时) nu <= 2.05
inline Real log_likelihood_term(Real residual, Real variance,
                                GarchDist dist, Real nu = 0.0) {
    if (!(variance > 0.0)) {
        throw std::invalid_argument(
            "log_likelihood_term: variance must be positive");
    }

    switch (dist) {
        case GarchDist::Normal: {
            // G3: 完整常数项 (arch: lls = -0.5*(log(2*pi)+log(sigma2)+resids**2/sigma2))
            const Real e2h = residual * residual / variance;
            return -0.5 * (std::log(2.0 * PI) + std::log(variance) + e2h);
        }
        case GarchDist::StudentT: {
            // arch 下界 2.05 (StudentsT.bounds), ν>2 保证方差存在
            if (!(nu > 2.05)) {
                throw std::invalid_argument(
                    "log_likelihood_term: Student-t requires nu > 2.05");
            }
            const Real half_nu_p1 = 0.5 * (nu + 1.0);
            const Real half_nu = 0.5 * nu;
            const Real u = residual * residual / ((nu - 2.0) * variance);
            return std::lgamma(half_nu_p1) - std::lgamma(half_nu)
                   - 0.5 * std::log(nu - 2.0) - 0.5 * std::log(PI)
                   - 0.5 * std::log(variance)
                   - half_nu_p1 * std::log(1.0 + u);
        }
        case GarchDist::GED: {
            // GED 形状参数数学上 ν>0 均有效 (ν=1 Laplace, ν=2 Normal);
            // 下界 1.0 排除数值退化区 (排幻觉点 G-nu-bound: 2.05 是 StudentsT
            // 的估计下界, 非 GED 似然下界 — spec 测试矩阵自身要求 GED ν=2 可用)
            if (!(nu > 1.0)) {
                throw std::invalid_argument(
                    "log_likelihood_term: GED requires nu > 1.0");
            }
            // G-GED1 (arch arch/univariate/distribution.py 实测):
            //   c = 2^(-1/ν) · sqrt(Γ(1/ν)/Γ(3/ν)); 指数为 ν (非 2ν)
            const Real c = std::pow(2.0, -1.0 / nu)
                           * std::sqrt(std::tgamma(1.0 / nu) / std::tgamma(3.0 / nu));
            const Real abs_std = std::abs(residual) / (std::sqrt(variance) * c);
            return std::log(nu) - std::log(c) - std::lgamma(1.0 / nu)
                   - (1.0 + 1.0 / nu) * std::log(2.0)
                   - 0.5 * std::log(variance)
                   - 0.5 * std::pow(abs_std, nu);
        }
    }
    throw std::invalid_argument("log_likelihood_term: unknown distribution");
}

/// @brief 批量对数似然 (逐项求和)
/// @throws std::invalid_argument 若两序列长度不一致或任一 variance <= 0
inline Real log_likelihood(const std::vector<Real>& residuals,
                           const std::vector<Real>& variances,
                           GarchDist dist, Real nu = 0.0) {
    if (residuals.size() != variances.size()) {
        throw std::invalid_argument(
            "log_likelihood: residuals and variances size mismatch");
    }
    Real total = 0.0;
    for (Size i = 0; i < residuals.size(); ++i) {
        total += log_likelihood_term(residuals[i], variances[i], dist, nu);
    }
    return total;
}

namespace detail {

// digamma ψ(x): 递推 + 渐近级数 (x >= 6 时 ψ(x) ≈ ln x - 1/(2x) - 1/(12x²)
//   + 1/(120x⁴) - 1/(252x⁶)); 小参数用 ψ(x) = ψ(x+1) - 1/x 递推到 x >= 6
inline Real digamma(Real x) {
    if (!(x > 0.0)) {
        throw std::invalid_argument("digamma: requires x > 0");
    }
    Real result = 0.0;
    while (x < 6.0) {
        result -= 1.0 / x;
        x += 1.0;
    }
    const Real inv = 1.0 / x;
    const Real inv2 = inv * inv;
    result += std::log(x) - 0.5 * inv
              - inv2 * (1.0 / 12.0 - inv2 * (1.0 / 120.0 - inv2 / 252.0));
    return result;
}

}  // namespace detail

/// @brief 似然梯度 (关于分布参数 ν 的解析得分; 其余分量为 0)
///
/// 说明: (residuals, variances) 给定时, ℓ 仅通过 ν 依赖参数 — (ω,α,β) 的
/// 得分需要方差递归路径 ∂hₜ/∂θⱼ, 无法从静态序列解析求出, 由估计器内部
/// 用中心差分处理 (spec §2.0.1: "未提供时 SLSQP 用中心差分")。
///
/// @return {dℓ/dω, dℓ/dα, dℓ/dβ, dℓ/dν} — 前三项恒为 0, Normal 时第 4 项亦为 0
inline std::vector<Real> log_likelihood_gradient(
    const std::vector<Real>& residuals,
    const std::vector<Real>& variances,
    GarchDist dist, Real nu = 0.0) {
    if (residuals.size() != variances.size()) {
        throw std::invalid_argument(
            "log_likelihood_gradient: residuals and variances size mismatch");
    }

    std::vector<Real> grad(4, 0.0);
    if (dist == GarchDist::Normal) {
        return grad;  // G3: Normal 无分布参数
    }

    if (dist == GarchDist::StudentT) {
        if (!(nu > 2.05)) {
            throw std::invalid_argument(
                "log_likelihood_gradient: Student-t requires nu > 2.05");
        }
        // 逐观测得分 (排幻觉: 常数项 0.5ψ((ν+1)/2)-0.5ψ(ν/2)-1/(2(ν-2))
        // 每个 ℓₜ 都含一份, 必须逐观测累加, 非总和后加一次):
        //   dℓₜ/dν = 0.5ψ((ν+1)/2) - 0.5ψ(ν/2) - 1/(2(ν-2))
        //             - 0.5·log(1+uₜ) + ((ν+1)/2)·uₜ/((ν-2)(1+uₜ))
        // 其中 uₜ = ε²ₜ/((ν-2)hₜ) (对 ν 求导时为常数)
        const Real const_part = 0.5 * detail::digamma(0.5 * (nu + 1.0))
                                - 0.5 * detail::digamma(0.5 * nu)
                                - 0.5 / (nu - 2.0);
        Real d = 0.0;
        for (Size i = 0; i < residuals.size(); ++i) {
            if (!(variances[i] > 0.0)) {
                throw std::invalid_argument(
                    "log_likelihood_gradient: variance must be positive");
            }
            const Real u = residuals[i] * residuals[i] / ((nu - 2.0) * variances[i]);
            const Real s = 1.0 + u;
            d += const_part - 0.5 * std::log(s)
                 + 0.5 * (nu + 1.0) * u / ((nu - 2.0) * s);
        }
        grad[3] = d;
        return grad;
    }

    // GED: ν 进入 c, Γ(1/ν), Γ(3/ν) 与指数 — 解析式冗长易错,
    // 用中心差分 (排幻觉: 不硬造解析公式)
    if (!(nu > 1.0)) {
        throw std::invalid_argument(
            "log_likelihood_gradient: GED requires nu > 1.0");
    }
    const Real delta = 1e-6 * std::max<Real>(1.0, std::abs(nu));
    const Real nup = nu + delta;
    const Real num = nu - delta;
    grad[3] = (log_likelihood(residuals, variances, dist, nup)
               - log_likelihood(residuals, variances, dist, num))
              / (2.0 * delta);
    return grad;
}

}  // namespace garch
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
