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

    // Log-of-ratio form: log((1 - g*e^{-dτ}) / (1 - g)).
    // This is mathematically equivalent to log(1 - g*e^{-dτ}) - log(1 - g) but
    // avoids branch-cut discontinuities that arise when the two individual logs
    // land on different branches (which happens for complex u in the P1 Carr-Madan
    // integral, where u_shifted = u_real - i).
    // For real u with Feller satisfied, |g| < 1, so the ratio is real and positive,
    // and the principal log is continuous — matching the RISK-015 direct form.
    Complex log_term = std::log((Real(1) - g * e_dt) / (Real(1) - g));

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

    // Log-of-ratio form (see heston_characteristic_function for rationale)
    Complex log_term = std::log((Real(1) - g * e_dt) / (Real(1) - g));

    Complex C = p.kappa * p.theta / (p.sigma * p.sigma) * ((x - d) * tau - Real(2) * log_term);

    return std::exp(C + D * p.v0);
}

}  // namespace v1
}  // namespace cpphub
