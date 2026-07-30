#pragma once
#include <cmath>
#include "cpphub/models/diffusion/heston.hpp"
#include "cpphub/core/math.hpp"

namespace cpphub {
inline namespace v1 {

class HestonQE : public Heston {
public:
    explicit HestonQE(HestonParams p) : Heston(p, HestonScheme::QE_M) {}

    void generate_path(Real T, Size n_steps,
                       std::span<Real> path,
                       Philox4x64& rng) const override;

private:
    void qe_step(Real dt, Real& v, Real z2) const;
    void apply_martingale_correction(Real& S, Real v_old, Real v_new,
                                      Real dt, Real z1, Real z2) const;
    void compute_moments(Real v, Real dt,
                          Real& m, Real& s2) const;

    static constexpr Real PSI_C = 1.5;
};

inline void HestonQE::compute_moments(Real v, Real dt,
                                       Real& m, Real& s2) const {
    Real kappa = params().kappa;
    Real theta = params().theta;
    Real sigma = params().sigma;

    Real exp_kdt = std::exp(-kappa * dt);

    m = theta + (v - theta) * exp_kdt;

    Real term1 = v * sigma * sigma * exp_kdt / kappa * (1.0 - exp_kdt);
    Real term2 = theta * sigma * sigma / (2.0 * kappa) * (1.0 - exp_kdt) * (1.0 - exp_kdt);
    s2 = term1 + term2;
    if (s2 < 0.0) s2 = 0.0;
}

inline void HestonQE::qe_step(Real dt, Real& v, Real z2) const {
    Real m, s2;
    compute_moments(v, dt, m, s2);

    Real psi = s2 / (m * m + 1e-30);

    if (psi <= PSI_C) {
        Real b_sq = 2.0 / psi - 1.0 + std::sqrt(2.0 / psi * (2.0 / psi - 1.0));
        Real b = std::sqrt(b_sq);
        Real a = m / (1.0 + b_sq);
        Real x = b + z2;
        v = a * x * x;
    } else {
        Real p = (psi - 1.0) / (psi + 1.0);
        Real beta = 2.0 / (m * psi + p);
        if (m * psi + p <= 0.0) {
            beta = 1.0 / m;
        }

        Real u = normal_cdf(z2);
        u = std::min(std::max(u, 1e-15), 1.0 - 1e-15);

        if (u <= p) {
            v = 0.0;
        } else {
            v = std::log((1.0 - p) / (1.0 - u)) / beta;
        }
    }
}

inline void HestonQE::apply_martingale_correction(Real& S, Real v_old, Real v_new,
                                                    Real dt, Real z1, Real z2) const {
    Real r = params().r;
    Real q = params().q;
    Real rho = params().rho;
    Real rho2 = rho * rho;
    Real sqrt_1_minus_rho2 = std::sqrt(1.0 - rho2);

    Real v = v_new;
    if (v < 0.0) v = 0.0;

    Real w = (z1 - rho * z2) / sqrt_1_minus_rho2;

    Real veff = v * (1.0 - rho2);
    S = S * std::exp((r - q - 0.5 * veff) * dt + std::sqrt(veff * dt) * w);
}

inline void HestonQE::generate_path(Real T, Size n_steps,
                                     std::span<Real> path,
                                     Philox4x64& rng) const {
    Real dt = T / static_cast<Real>(n_steps);
    path[0] = params().S0;

    Real S = params().S0;
    Real v = params().v0;
    Real rho = params().rho;
    Real sqrt_1_minus_rho2 = std::sqrt(1.0 - rho * rho);

    for (Size i = 0; i < n_steps; ++i) {
        uint64_t r1 = rng();
        uint64_t r2 = rng();
        Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
        Real u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
        auto [w1, w2] = box_muller(u1, u2);

        Real z1 = w1;
        Real z2 = rho * w1 + sqrt_1_minus_rho2 * w2;

        Real v_old = v;
        qe_step(dt, v, z2);

        apply_martingale_correction(S, v_old, v, dt, z1, z2);

        path[i + 1] = S;
    }
}

}  // namespace v1
}  // namespace cpphub
