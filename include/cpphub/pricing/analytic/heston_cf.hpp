#pragma once
#include <complex>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"

namespace cpphub {
inline namespace v1 {

struct HestonCFParams {
    Real v0;
    Real kappa;
    Real theta;
    Real sigma;
    Real rho;
    Real r;
    Real q;
};

Complex heston_characteristic_function(Complex u, Real tau,
                                        Real S0, const HestonCFParams& p);

Complex heston_variance_cf(Complex u, Real tau, const HestonCFParams& p);

inline Complex heston_characteristic_function(Complex u, Real tau,
                                               Real S0, const HestonCFParams& p) {
    if (std::abs(u) < Real(1e-15)) {
        return Complex(1, 0);
    }

    const Complex i(0, 1);

    Complex x = p.kappa - p.rho * p.sigma * u * i;

    Complex d = std::sqrt(x * x + p.sigma * p.sigma * (u * u + i * u));

    Complex g = (x - d) / (x + d);

    Complex e_dt = std::exp(-d * tau);

    Complex D = (x - d) / (p.sigma * p.sigma) * (Real(1) - e_dt) / (Real(1) - g * e_dt);

    // Direct log form: log(1 - g*e^{-dτ}) - log(1 - g).
    // When Feller condition (2κθ > σ²) holds, |g| < 1 for real u, so both
    // 1 - g*e^{-dτ} and 1 - g lie in the right half-plane (Re > 0), and the
    // principal log is continuous — no branch-cut discontinuity.
    // The previous "Little Trap" rewrite log(1/g - e^{-dτ}) - log(1/g - 1)
    // introduced a spurious branch jump when Im(log(g)) crossed zero.
    Complex log_term = std::log(Real(1) - g * e_dt) - std::log(Real(1) - g);

    Complex C = p.kappa * p.theta / (p.sigma * p.sigma) * ((x - d) * tau - Real(2) * log_term);
    C += i * u * (std::log(S0) + (p.r - p.q) * tau);

    return std::exp(C + D * p.v0);
}

inline Complex heston_variance_cf(Complex u, Real tau, const HestonCFParams& p) {
    if (std::abs(u) < Real(1e-15)) {
        return Complex(1, 0);
    }

    const Complex i(0, 1);

    Complex x = p.kappa - p.rho * p.sigma * u * i;

    Complex d = std::sqrt(x * x + p.sigma * p.sigma * (u * u + i * u));

    Complex g = (x - d) / (x + d);

    Complex e_dt = std::exp(-d * tau);

    Complex D = (x - d) / (p.sigma * p.sigma) * (Real(1) - e_dt) / (Real(1) - g * e_dt);

    // Direct log form (see heston_characteristic_function for rationale)
    Complex log_term = std::log(Real(1) - g * e_dt) - std::log(Real(1) - g);

    Complex C = p.kappa * p.theta / (p.sigma * p.sigma) * ((x - d) * tau - Real(2) * log_term);

    return std::exp(C + D * p.v0);
}

}  // namespace v1
}  // namespace cpphub
