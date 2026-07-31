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
// 实现:
//   x < a+1: 级数展开 γ(a,x) = x^a e^{-x} Σ_{n=0}^∞ x^n / (a(a+1)...(a+n))
//   x >= a+1: 连分式展开 Q(a,x) = 1 - P(a,x) = e^{-x} x^a / Γ(a) · CF(x, a)
//
// 精度: 接近机器精度 (~1e-14)
inline Real regularized_lower_gamma(Real a, Real x) {
    if (x < 0.0 || a <= 0.0) {
        throw std::invalid_argument("regularized_lower_gamma: a > 0 and x >= 0 required");
    }
    if (x == 0.0) return 0.0;

    const Real LOG_EPS = std::log(std::numeric_limits<Real>::epsilon());

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
    } else {
        // 连分式 (Lentz 算法), 计算 Q(a,x) = 1 - P(a,x)
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
        Real Q = std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
        return 1.0 - Q;
    }
}

// 正则化上不完全 Gamma 函数 Q(a, x) = 1 - P(a, x) = Γ(a,x) / Γ(a)
inline Real regularized_upper_gamma(Real a, Real x) {
    return 1.0 - regularized_lower_gamma(a, x);
}

// ==================== 非中心卡方分布 CDF ====================
// SOURCE: Noncentral chi-squared distribution
//   X ~ χ'²(k, λ), 自由度 k > 0, 非中心参数 λ >= 0
//   CDF: P(X ≤ x) = e^{-λ/2} Σ_{j=0}^∞ (λ/2)^j / j! · P(k/2 + j, x/2)
//   其中 P(a, x) 是正则化下不完全 Gamma 函数
//
// 实现:
//   - 小 λ (≤ 50): 直接级数求和, exp(-λ/2) 不下溢
//   - 大 λ (> 50): 对数空间递推 Poisson 权重, 避免 exp(-λ/2) 下溢
//     Poisson 众数 j* ≈ λ/2, 在 [j* - 5√(λ/2), j* + 5√(λ/2)] 范围求和
//     覆盖 > 99.9999% 概率质量
//
// 精度: 小 λ 接近机器精度 (~1e-12); 大 λ ~1e-6 (Poisson 尾部截断)
inline Real noncentral_chi2_cdf(Real x, Real k, Real lambda) {
    if (k <= 0.0) throw std::invalid_argument("noncentral_chi2_cdf: k must be positive");
    if (lambda < 0.0) throw std::invalid_argument("noncentral_chi2_cdf: lambda must be non-negative");
    if (x <= 0.0) return 0.0;

    // 中心卡方 (λ=0) 直接用 Gamma CDF
    if (lambda == 0.0) {
        return regularized_lower_gamma(k / 2.0, x / 2.0);
    }

    const Real lambda_half = lambda / 2.0;

    // ---- 小 λ (≤ 50): 直接级数求和 ----
    if (lambda <= 50.0) {
        Real poisson_weight = std::exp(-lambda_half);  // j=0: e^{-λ/2}
        Real sum = poisson_weight * regularized_lower_gamma(k / 2.0, x / 2.0);
        for (int j = 1; j < 200; ++j) {
            poisson_weight *= lambda_half / static_cast<Real>(j);
            Real term = poisson_weight * regularized_lower_gamma(k / 2.0 + j, x / 2.0);
            sum += term;
            if (poisson_weight < 1e-16 && term < 1e-16 * sum) break;
        }
        if (sum < 0.0) return 0.0;
        if (sum > 1.0) return 1.0;
        return sum;
    }

    // ---- 大 λ (> 50): 对数空间递推 ----
    // log_w_j = -λ/2 + j·log(λ/2) - lgamma(j+1)
    // 递推: log_w_{j+1} = log_w_j + log(λ/2) - log(j+1)
    const Real log_lambda_half = std::log(lambda_half);
    const Real sqrt_lambda_half = std::sqrt(lambda_half);

    // 求和范围: Poisson 众数 ± 5σ (覆盖 > 99.9999%)
    int j_start = static_cast<int>(std::max(Real(0.0), lambda_half - 5.0 * sqrt_lambda_half));
    int j_end = static_cast<int>(lambda_half + 5.0 * sqrt_lambda_half + 10.0);

    // 初始 log_w_{j_start}
    Real log_w = -lambda_half
                 + static_cast<Real>(j_start) * log_lambda_half
                 - std::lgamma(static_cast<Real>(j_start) + 1.0);

    Real sum = 0.0;
    for (int j = j_start; j <= j_end; ++j) {
        if (j > j_start) {
            log_w += log_lambda_half - std::log(static_cast<Real>(j));
        }
        // w_j = exp(log_w); 对极小值用 0 (下溢保护)
        if (log_w > -700.0) {
            Real w = std::exp(log_w);
            Real term = w * regularized_lower_gamma(k / 2.0 + static_cast<Real>(j), x / 2.0);
            sum += term;
        }
    }
    if (sum < 0.0) return 0.0;
    if (sum > 1.0) return 1.0;
    return sum;
}

}  // namespace v1
}  // namespace cpphub
