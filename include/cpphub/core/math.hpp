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

Real inv_normal_cdf(Real p) {
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

}  // namespace v1
}  // namespace cpphub
