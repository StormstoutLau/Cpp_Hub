#pragma once
#include <cmath>
#include <complex>
#include <span>
#include "cpphub/models/diffusion/process.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include "cpphub/core/linalg.hpp"
#include "cpphub/core/constants.hpp"

namespace cpphub {
inline namespace v1 {

struct HestonParams {
    Real S0;
    Real v0;
    Real kappa;
    Real theta;
    Real sigma;
    Real rho;
    Real r;
    Real q;
};

enum class HestonScheme {
    Euler,
    FullTruncation,
    QE_M,
    Exact
};

class Heston : public StochasticProcess {
public:
    explicit Heston(HestonParams p, HestonScheme scheme = HestonScheme::FullTruncation)
        : params_(p), scheme_(scheme) {}

    Size dimension() const override { return 2; }
    Real spot() const override { return params_.S0; }

    Complex characteristic_function(Complex u, Real tau) const override;

    void generate_path(Real T, Size n_steps,
                       std::span<Real> path,
                       Philox4x64& rng) const override;

    const HestonParams& params() const { return params_; }

private:
    HestonParams params_;
    HestonScheme scheme_;

    void euler_step(Real dt, Real& S, Real& v,
                    Real z1, Real z2) const;
    void full_truncation_step(Real dt, Real& S, Real& v,
                              Real z1, Real z2) const;
    void generate_correlated_normals(Philox4x64& rng,
                                      Real& z1, Real& z2) const;
};

inline void Heston::generate_correlated_normals(Philox4x64& rng,
                                                  Real& z1, Real& z2) const {
    uint64_t r1 = rng();
    uint64_t r2 = rng();
    double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
    double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
    auto [w1, w2] = box_muller(u1, u2);
    Real rho = params_.rho;
    z1 = w1;
    z2 = rho * w1 + std::sqrt(1.0 - rho * rho) * w2;
}

inline void Heston::euler_step(Real dt, Real& S, Real& v,
                                Real z1, Real z2) const {
    Real v_old = v;
    Real v_pos = (v_old > 0.0) ? v_old : 0.0;
    Real sqrt_v = std::sqrt(v_pos);

    Real kappa = params_.kappa;
    Real theta = params_.theta;
    Real sigma = params_.sigma;
    Real r = params_.r;
    Real q = params_.q;

    v = v_old + kappa * (theta - v_old) * dt + sigma * sqrt_v * z2 * std::sqrt(dt);
    if (v < 0.0) v = 0.0;

    S = S * std::exp((r - q - 0.5 * v_pos) * dt + sqrt_v * z1 * std::sqrt(dt));
}

inline void Heston::full_truncation_step(Real dt, Real& S, Real& v,
                                          Real z1, Real z2) const {
    Real v_old = v;
    Real v_pos = (v_old > 0.0) ? v_old : 0.0;
    Real sqrt_v = std::sqrt(v_pos);

    Real kappa = params_.kappa;
    Real theta = params_.theta;
    Real sigma = params_.sigma;
    Real r = params_.r;
    Real q = params_.q;

    v = v_old + kappa * (theta - v_old) * dt + sigma * sqrt_v * z2 * std::sqrt(dt);
    if (v < 0.0) v = 0.0;

    S = S * std::exp((r - q - 0.5 * v_pos) * dt + sqrt_v * z1 * std::sqrt(dt));
}

inline void Heston::generate_path(Real T, Size n_steps,
                                   std::span<Real> path,
                                   Philox4x64& rng) const {
    Real dt = T / static_cast<Real>(n_steps);
    path[0] = params_.S0;

    if (scheme_ == HestonScheme::QE_M) {
        return;
    }

    Real S = params_.S0;
    Real v = params_.v0;

    for (Size i = 0; i < n_steps; ++i) {
        Real z1, z2;
        generate_correlated_normals(rng, z1, z2);

        switch (scheme_) {
            case HestonScheme::Euler:
                euler_step(dt, S, v, z1, z2);
                break;
            case HestonScheme::FullTruncation:
                full_truncation_step(dt, S, v, z1, z2);
                break;
            default:
                break;
        }
        path[i + 1] = S;
    }
}

inline Complex Heston::characteristic_function(Complex u, Real tau) const {
    // Delegate to Little Trap (Albrecher 2007) implementation to avoid the
    // branch-cut discontinuity present in the original Heston (1993) form.
    // The original form computes log((1 - g*e^{-dτ})/(1 - g)) directly, which
    // jumps to the wrong Riemann sheet when the argument crosses the negative
    // real axis — producing incorrect CF values for certain (u, τ) pairs.
    // See: Albrecher et al. (2007) "The Little Trap of Heston's Stochastic
    // Volatility Model", and Schoutens' corrected reference table.
    HestonCFParams p{params_.v0, params_.kappa, params_.theta,
                     params_.sigma, params_.rho, params_.r, params_.q};
    return heston_characteristic_function(u, tau, params_.S0, p);
}

}  // namespace v1
}  // namespace cpphub
