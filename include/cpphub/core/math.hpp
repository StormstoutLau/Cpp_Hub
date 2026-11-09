#pragma once
#include <cmath>
#include <limits>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"

namespace cpphub {
inline namespace v1 {

inline Real erf(Real x) noexcept {
    return std::erf(x);
}

inline Real erfc(Real x) noexcept {
    return std::erfc(x);
}

inline Real normal_pdf(Real x) noexcept {
    return INV_SQRT_2PI * std::exp(-x * x / Real(2));
}

inline Real normal_cdf(Real x) noexcept {
    return Real(0.5) * std::erfc(-x / SQRT_2);
}

inline Real bessel_i0(Real x) noexcept {
    Real ax = std::abs(x);
    Real x2 = ax * 0.5;
    Real s = x2 * x2;
    Real sum = 1.0;
    Real term = 1.0;
    int kmax = (ax < 3.75) ? 15 : 50;
    for (int k = 1; k <= kmax; ++k) {
        term *= s / (Real(k) * Real(k));
        Real prev = sum;
        sum += term;
        if (sum == prev) break;
    }
    return sum;
}

inline Real bessel_i1(Real x) noexcept {
    Real ax = std::abs(x);
    Real x2 = ax * 0.5;
    Real s = x2 * x2;
    Real sum = 1.0;
    Real term = 1.0;
    int kmax = (ax < 3.75) ? 15 : 50;
    for (int k = 1; k <= kmax; ++k) {
        term *= s / (Real(k) * Real(k + 1));
        Real prev = sum;
        sum += term;
        if (sum == prev) break;
    }
    Real y = x2 * sum;
    return (x >= 0) ? y : -y;
}

inline Real inv_normal_cdf(Real p) {
    if (p <= 0) return -std::numeric_limits<Real>::infinity();
    if (p >= 1) return std::numeric_limits<Real>::infinity();

    bool negate = (p < 0.5);
    Real q = negate ? p : (1 - p);
    Real t = std::sqrt(-2 * std::log(q));

    static const Real c0 = 2.515517;
    static const Real c1 = 0.802853;
    static const Real c2 = 0.010328;
    static const Real d0 = 1.432788;
    static const Real d1 = 0.189269;
    static const Real d2 = 0.001308;

    Real x = t - (c0 + c1 * t + c2 * t * t) /
                 (1.0 + d0 * t + d1 * t * t + d2 * t * t * t);
    if (negate) x = -x;

    for (int i = 0; i < 3; ++i) {
        Real f = normal_cdf(x) - p;
        Real fp = normal_pdf(x);
        Real denom = 2 * fp + f * x;
        if (std::abs(denom) > std::numeric_limits<Real>::min()) {
            x -= 2 * f / denom;
        }
    }

    return x;
}

// ==================== 不完全 Gamma 函数 ====================
// SOURCE: Numerical Recipes (3rd ed.) §6.2
// 正则化下不完全 Gamma 函数 P(a, x) = γ(a, x) / Γ(a)
// 满足 P(a, 0) = 0, P(a, ∞) = 1
//
// RISK-036: P 与 Q 严格互补 — P 走级数 (x < a+1), Q 走连分式 (x >= a+1),
// 各自在自己支路内相对精确; 任何一侧用减法从另一侧推导都会在尾部灾难性失精度。
//
// 实现:
//   x < a+1: 级数展开 γ(a,x) = x^a e^{-x} Σ_{n=0}^∞ x^n / (a(a+1)...(a+n))
//   x >= a+1: 连分式展开 Q(a,x) = 1 - P(a,x) = e^{-x} x^a / Γ(a) · CF(x, a)
//
// 精度: 接近机器精度 (~1e-14)
namespace detail {
// Lentz 连分式直接计算 Q(a, x) = Γ(a,x)/Γ(a) (仅用于 x >= a+1 支路;
// 单一来源供 regularized_lower/upper_gamma 共用, 保证 P/Q 逐位互补)
inline Real gamma_q_continued_fraction(Real a, Real x) {
    Real b = x + 1.0 - a;
    Real c = 1.0 / std::numeric_limits<Real>::min();
    Real d = 1.0 / b;
    Real h = d;
    for (int i = 1; i <= 1000; ++i) {
        Real an = -static_cast<Real>(i) * (static_cast<Real>(i) - a);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < std::numeric_limits<Real>::min()) d = std::numeric_limits<Real>::min();
        c = b + an / c;
        if (std::abs(c) < std::numeric_limits<Real>::min()) c = std::numeric_limits<Real>::min();
        d = 1.0 / d;
        Real del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < std::numeric_limits<Real>::epsilon()) break;
    }
    return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
}
}  // namespace detail

inline Real regularized_lower_gamma(Real a, Real x) {
    if (x < 0.0 || a <= 0.0) {
        throw std::invalid_argument("regularized_lower_gamma: a > 0 and x >= 0 required");
    }
    if (x == 0.0) return 0.0;

    if (x < a + 1.0) {
        // 级数展开
        Real ap = a;
        Real sum = 1.0 / a;
        Real del = sum;
        for (int n = 0; n < 1000; ++n) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::abs(del) < std::abs(sum) * std::numeric_limits<Real>::epsilon()) break;
        }
        return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
    }
    // 连分式: Q 在此支路直接算出, P = 1 - Q (x >= a+1 区间 P 接近 1, 减法无损)
    return 1.0 - detail::gamma_q_continued_fraction(a, x);
}

// 正则化上不完全 Gamma 函数 Q(a, x) = Γ(a,x) / Γ(a)
//
// RISK-036: 原实现 1 - P 在 Q < ~1.1e-16 (k=1/LR≳70 的 χ² 尾部) 时舍入为 0
// (Goldberg 1991 §2.4 cancellation; NR §6.2 gammq 在连分式支路直接返回 Q)。
// 重写后: x >= a+1 直接返回连分式 Q (全尾精度到 underflow); x < a+1 区间
// Q >= ~0.083, 1 - P 相对误差 <= 3e-15, 安全。
inline Real regularized_upper_gamma(Real a, Real x) {
    if (x < 0.0 || a <= 0.0) {
        throw std::invalid_argument("regularized_upper_gamma: a > 0 and x >= 0 required");
    }
    if (x == 0.0) return 1.0;
    if (x < a + 1.0) {
        return 1.0 - regularized_lower_gamma(a, x);
    }
    return detail::gamma_q_continued_fraction(a, x);
}

// χ²(k) 上尾生存函数 sf(x) = P(χ²_k > x) = Q(k/2, x/2)
//
// RISK-036: backtesting.hpp (Kupiec/Christoffersen) 的 p 值入口; 旧
// 1 - chi2_cdf 路径在 p < ~1e-16 时归零。x <= 0 → 1 (无证据 → 不拒绝,
// 保持"负 LR → p=1"语义)。
inline Real chi2_sf(Real x, Real k) {
    if (!(k > 0.0)) {
        throw std::invalid_argument("chi2_sf: k must be positive");
    }
    if (x <= 0.0) return 1.0;
    return regularized_upper_gamma(k / Real(2), x / Real(2));
}

// ==================== 非中心卡方分布 CDF / SF ====================
// SOURCE: Noncentral chi-squared distribution
//   X ~ χ'²(k, λ), 自由度 k > 0, 非中心参数 λ >= 0
//   CDF: P(X ≤ x) = e^{-λ/2} Σ_{j=0}^∞ (λ/2)^j / j! · P(k/2 + j, x/2)
//   SF : P(X > x) = e^{-λ/2} Σ_{j=0}^∞ (λ/2)^j / j! · Q(k/2 + j, x/2)
//   其中 P/Q 是正则化下/上不完全 Gamma 函数
//
// RISK-039 (CEV 专项批):
//   旧实现大 λ 分支固定 5σ Poisson 窗口, CDF 绝对误差 ~1e-7; 且 CEV 的
//   P1 = 1 - CDF 在深 OTM 暴露该误差 (P1 地板 3.17e-7, 价格地板 3.17e-5).
//   重写: ① 新增 noncentral_chi2_sf 上尾直算; ② 大 λ 窗口加宽到 ±10σ,
//   并扩展覆盖 Gamma 因子有效区 (a = k/2+j 直到 x/2 + 10√(x/2)),
//   因为 Q(a, x/2) 在 a ≈ x/2 附近达到峰值, 可能远离 Poisson 众数.
//
// 精度: 小 λ 接近机器精度; 大 λ ~1e-13 (Poisson 10σ 尾部 + Gamma 相对项)
namespace detail {
// 非中心卡方 Poisson 混合, upper=true 返回上尾 Q 的混合, false 返回下尾 P
inline Real ncchi2_mixture(Real x, Real k, Real lambda, bool upper) {
    const Real lambda_half = lambda / 2.0;
    const Real x_half = x / 2.0;

    auto gamma = [upper](Real a, Real z) -> Real {
        return upper ? regularized_upper_gamma(a, z)
                     : regularized_lower_gamma(a, z);
    };

    // ---- 小 λ (≤ 50): 直接级数求和 ----
    if (lambda <= 50.0) {
        Real poisson_weight = std::exp(-lambda_half);
        Real sum = poisson_weight * gamma(k / 2.0, x_half);
        for (int j = 1; j < 400; ++j) {
            poisson_weight *= lambda_half / static_cast<Real>(j);
            Real term = poisson_weight * gamma(k / 2.0 + j, x_half);
            sum += term;
            if (poisson_weight < 1e-17 && std::abs(term) < 1e-17 * sum) break;
            if (poisson_weight < 1e-300) break;
        }
        return sum;
    }

    // ---- 大 λ: 对数空间 Poisson 权重 + 自适应窗口 ----
    const Real log_lambda_half = std::log(lambda_half);
    const Real sqrt_lambda_half = std::sqrt(lambda_half);

    // 窗口下界: Poisson 主体 -10σ (更低 j 的权重 < ~1e-23)
    long long j_lo = static_cast<long long>(std::floor(
        std::max(Real(0.0), lambda_half - 10.0 * sqrt_lambda_half)));
    long long j_hi_poisson = static_cast<long long>(std::floor(
        lambda_half + 10.0 * sqrt_lambda_half + 10.0));
    // Gamma 因子有效区: a = k/2+j 到 x/2 + 10√(x/2) (Q 的峰值在 a≈x/2)
    long long j_hi_gamma = static_cast<long long>(std::floor(
        x_half + 10.0 * std::sqrt(x_half + 1.0) + 10.0));
    long long j_hi = std::max(j_hi_poisson, j_hi_gamma);
    if (j_hi > 2000000) j_hi = 2000000;  // 病态参数保护

    Real log_w = -lambda_half + static_cast<Real>(j_lo) * log_lambda_half
                 - std::lgamma(static_cast<Real>(j_lo) + 1.0);

    Real sum = 0.0;
    for (long long j = j_lo; j <= j_hi; ++j) {
        if (j > j_lo) {
            log_w += log_lambda_half - std::log(static_cast<Real>(j));
        }
        if (log_w > -745.0) {
            Real w = std::exp(log_w);
            sum += w * gamma(k / 2.0 + static_cast<Real>(j), x_half);
        }
    }
    return sum;
}
}  // namespace detail

inline Real noncentral_chi2_cdf(Real x, Real k, Real lambda) {
    if (k <= 0.0) throw std::invalid_argument("noncentral_chi2_cdf: k must be positive");
    if (lambda < 0.0) throw std::invalid_argument("noncentral_chi2_cdf: lambda must be non-negative");
    if (x <= 0.0) return 0.0;
    if (lambda == 0.0) {
        return regularized_lower_gamma(k / 2.0, x / 2.0);
    }
    Real sum = detail::ncchi2_mixture(x, k, lambda, false);
    if (sum < 0.0) return 0.0;
    if (sum > 1.0) return 1.0;
    return sum;
}

// 非中心卡方上尾 SF(x; k, λ) = P(X > x) — 深尾直算, 避免 1−CDF 抵消
inline Real noncentral_chi2_sf(Real x, Real k, Real lambda) {
    if (k <= 0.0) throw std::invalid_argument("noncentral_chi2_sf: k must be positive");
    if (lambda < 0.0) throw std::invalid_argument("noncentral_chi2_sf: lambda must be non-negative");
    if (x <= 0.0) return 1.0;
    if (lambda == 0.0) {
        return regularized_upper_gamma(k / 2.0, x / 2.0);
    }
    Real sum = detail::ncchi2_mixture(x, k, lambda, true);
    if (sum < 0.0) return 0.0;
    if (sum > 1.0) return 1.0;
    return sum;
}

}  // namespace v1
}  // namespace cpphub
