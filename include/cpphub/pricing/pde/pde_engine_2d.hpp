// 2D PDE engine for Heston model using ADI schemes (Craig-Sneyd / Hundsdorfer-Verwer)
//
// Solves the Heston PDE in log-spot / variance space:
//   dV/dtau = 0.5 v d2V/dx2 + (r-q-0.5v) dV/dx
//           + 0.5 xi^2 v d2V/dv2 + kappa(theta-v) dV/dv
//           + rho xi v d2V/(dx dv) - r V
// where x = ln(S/K), v = variance, tau = T - t (time to maturity).
//
// Boundary conditions:
//   x_min (S -> 0):   V = 0 (call), V = K*exp(-r*tau) (put)
//   x_max (S -> inf): V = K*exp(x - q*tau) - K*exp(-r*tau) (call), V = 0 (put)
//   v_min (v -> 0):   V = degenerate PDE solution (BSM-like with v=0)
//   v_max (v -> inf): linear extrapolation (Neumann-like)
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/pricing/pde/fdm_2d.hpp"

namespace cpphub {
inline namespace v1 {

enum class ADISchemeType {
    CraigSneyd,
    HundsdorferVerwer
};

struct PDEEngine2DConfig {
    Size n_x = 100;         // log-S grid points
    Size n_v = 80;          // variance grid points
    Size n_time = 200;      // time steps
    Real x_range = 4.0;     // absolute half-width in log-S/K space (x_min = x0 - x_range, x_max = x0 + x_range)
    Real v_range = 5.0;     // v_max = v0 + v_range * sqrt(v0) (heuristic)
    Real v_min_val = 0.0;   // v_min (typically 0)
    ADISchemeType scheme = ADISchemeType::CraigSneyd;
    Real theta = 0.5;       // ADI theta (0.5 = second-order)
    bool is_call = true;
};

class PDEEngine2D {
public:
    explicit PDEEngine2D(PDEEngine2DConfig config = PDEEngine2DConfig{})
        : config_(config) {}

    // Price European option under Heston model via ADI
    Real price(const HestonPDEParams& hp) const {
        // Build grid
        Real x0 = std::log(hp.S0 / hp.K);
        Real sigma0 = std::sqrt(hp.v0);
        Real x_min = x0 - config_.x_range;
        Real x_max = x0 + config_.x_range;
        Real v_max = hp.v0 + config_.v_range * sigma0;  // heuristic upper bound
        if (v_max <= hp.v0) v_max = hp.v0 * 5.0 + 0.01;

        FDMGrid2D grid(config_.n_x, config_.n_v, x_min, x_max,
                       config_.v_min_val, v_max);
        auto coeffs = compute_heston_coeffs(grid, hp);

        // Initial condition: payoff at maturity (tau = 0)
        // V(x, v, 0) = max(K*exp(x) - K, 0) for call (x = ln(S/K) => S = K*exp(x))
        //            = max(K - K*exp(x), 0) for put
        Size total = grid.size();
        std::vector<Real> V(total), V_new(total);
        for (Size j = 0; j < grid.n_v(); ++j) {
            for (Size i = 0; i < grid.n_x(); ++i) {
                Real S = hp.K * std::exp(grid.x(i));
                V[grid.idx(i, j)] = config_.is_call
                    ? std::max(S - hp.K, 0.0)
                    : std::max(hp.K - S, 0.0);
            }
        }

        // Time integration (tau from 0 to T)
        Real dt = hp.T / static_cast<Real>(config_.n_time);
        auto scheme = create_scheme();

        for (Size step = 0; step < config_.n_time; ++step) {
            // Apply boundary conditions before step
            apply_boundary_conditions(V, grid, hp,
                                       static_cast<Real>(step) * dt);
            scheme->step(V, V_new, dt, grid, coeffs);
            // Apply boundary conditions after step (at tau_{n+1})
            apply_boundary_conditions(V_new, grid, hp,
                                       static_cast<Real>(step + 1) * dt);
            V.swap(V_new);
        }

        // Interpolate at (x0, v0)
        return interpolate_at(V, grid, x0, hp.v0);
    }

    const PDEEngine2DConfig& config() const { return config_; }

private:
    PDEEngine2DConfig config_;

    std::unique_ptr<ADISchemeBase> create_scheme() const {
        switch (config_.scheme) {
            case ADISchemeType::HundsdorferVerwer:
                return std::make_unique<HundsdorferVerwerScheme>(config_.theta);
            case ADISchemeType::CraigSneyd:
            default:
                return std::make_unique<CraigSneydScheme>(config_.theta);
        }
    }

    // Boundary conditions
    void apply_boundary_conditions(std::vector<Real>& V, const FDMGrid2D& grid,
                                    const HestonPDEParams& hp, Real tau) const {
        Size nx = grid.n_x(), nv = grid.n_v();
        Real r = hp.r, q = hp.q, K = hp.K;
        Real disc_r = std::exp(-r * tau);
        Real disc_q = std::exp(-q * tau);

        // x boundaries (x_min: S->0, x_max: S->inf)
        for (Size j = 0; j < nv; ++j) {
            if (config_.is_call) {
                V[grid.idx(0, j)] = 0.0;  // S -> 0 => call = 0
                Real S_max = K * std::exp(grid.x_max());
                V[grid.idx(nx - 1, j)] = S_max * disc_q - K * disc_r;
            } else {
                V[grid.idx(0, j)] = K * disc_r;  // S -> 0 => put = K*exp(-r*tau)
                V[grid.idx(nx - 1, j)] = 0.0;    // S -> inf => put = 0
            }
        }

        // v boundaries
        // v_min (=0): degenerate PDE. For Heston, dV/dv term vanishes.
        // Simple approximation: V(x, 0, tau) = BSM with sigma=0 (deterministic drift)
        //   S(tau) = S0 * exp((r-q)*tau) => call payoff at tau:
        //   V(x, 0, tau) = max(K*exp(x + (r-q)*tau) - K, 0) * exp(-r*tau) for call
        //                = max(K - K*exp(x + (r-q)*tau), 0) * exp(-r*tau) for put
        if (grid.v_min() == 0.0) {
            for (Size i = 0; i < nx; ++i) {
                Real x = grid.x(i);
                Real S_tau = K * std::exp(x + (r - q) * tau);
                if (config_.is_call) {
                    V[grid.idx(i, 0)] = std::max(S_tau - K, 0.0) * disc_r;
                } else {
                    V[grid.idx(i, 0)] = std::max(K - S_tau, 0.0) * disc_r;
                }
            }
        }

        // v_max: linear extrapolation (Neumann-like)
        if (nv >= 3) {
            for (Size i = 0; i < nx; ++i) {
                Real v0 = V[grid.idx(i, nv - 3)];
                Real v1 = V[grid.idx(i, nv - 2)];
                Real v2 = V[grid.idx(i, nv - 1)];
                // Linear extrapolation: V[nv-1] = 2*V[nv-2] - V[nv-3]
                // But this overrides the corner. Keep it simple.
                V[grid.idx(i, nv - 1)] = 2.0 * v1 - v0;
                (void)v2;
            }
        }
    }

    Real interpolate_at(const std::vector<Real>& V, const FDMGrid2D& grid,
                         Real x0, Real v0) const {
        // Bilinear interpolation at (x0, v0)
        Size nx = grid.n_x(), nv = grid.n_v();
        // Find i such that grid.x(i) <= x0 < grid.x(i+1)
        Size i = 0;
        for (Size k = 0; k < nx - 1; ++k) {
            if (grid.x(k) <= x0 && x0 <= grid.x(k + 1)) { i = k; break; }
            if (k == nx - 2) i = k;
        }
        Size j = 0;
        for (Size k = 0; k < nv - 1; ++k) {
            if (grid.v(k) <= v0 && v0 <= grid.v(k + 1)) { j = k; break; }
            if (k == nv - 2) j = k;
        }

        Real x_l = grid.x(i), x_r = grid.x(i + 1);
        Real v_b = grid.v(j), v_t = grid.v(j + 1);
        Real wx = (x0 - x_l) / (x_r - x_l);
        Real wv = (v0 - v_b) / (v_t - v_b);

        Real V_bl = V[grid.idx(i, j)];
        Real V_br = V[grid.idx(i + 1, j)];
        Real V_tl = V[grid.idx(i, j + 1)];
        Real V_tr = V[grid.idx(i + 1, j + 1)];

        Real V_b = V_bl * (1.0 - wx) + V_br * wx;
        Real V_t = V_tl * (1.0 - wx) + V_tr * wx;
        return V_b * (1.0 - wv) + V_t * wv;
    }
};

}  // namespace v1
}  // namespace cpphub
