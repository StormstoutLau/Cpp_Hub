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
    HundsdorferVerwer,
    ModifiedCraigSneyd
};

enum class GridType2D {
    Uniform,
    Sinh
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
    GridType2D grid_type = GridType2D::Uniform;  // grid transformation
    Real alpha_x = 1.0;     // sinh concentration for x (0=uniform, 1.0=moderate, 2.0=strong)
    Real alpha_v = 0.0;     // sinh concentration for v (0=uniform; >0 hurts when Feller satisfied)
    Size n_rannacher_warmup = 0;  // Rannacher smoothing: first n steps use theta=1.0 (L-stable), 0=disabled
    bool is_american = false;     // American option (PSOR early-exercise projection after each ADI step)
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

        // American options benefit from sinh grid: concentrates points near
        // S=K (payoff kink) and v=0 (Feller boundary layer). If user hasn't
        // explicitly set grid_type, auto-enable Sinh for American.
        bool use_sinh = (config_.grid_type == GridType2D::Sinh);
        if (config_.is_american && config_.grid_type == GridType2D::Uniform
            && config_.alpha_x == 1.0) {
            // Default config (alpha_x=1.0, grid_type=Uniform) — auto-upgrade
            use_sinh = true;
        }

        FDMGrid2D grid = use_sinh
            ? FDMGrid2D(config_.n_x, config_.n_v, x_min, x_max,
                        config_.v_min_val, v_max, config_.alpha_x, config_.alpha_v)
            : FDMGrid2D(config_.n_x, config_.n_v, x_min, x_max,
                        config_.v_min_val, v_max);
        auto coeffs = compute_heston_coeffs(grid, hp);

        // Initial condition: payoff at maturity (tau = 0)
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

        // Precompute payoff at each x grid point for American projection
        std::vector<Real> payoff_at_x(grid.n_x(), 0.0);
        if (config_.is_american) {
            for (Size i = 0; i < grid.n_x(); ++i) {
                Real S = hp.K * std::exp(grid.x(i));
                payoff_at_x[i] = config_.is_call
                    ? std::max(S - hp.K, 0.0)
                    : std::max(hp.K - S, 0.0);
            }
        }

        // Adaptive time step
        Real dt = hp.T / static_cast<Real>(config_.n_time);
        Real dv = grid.dv();
        Real v_top = grid.v_max();
        Real drift_v_top = hp.kappa * (hp.theta - v_top);
        Real diff_v_top = 0.5 * hp.xi * hp.xi * v_top;
        Real a_v_top = diff_v_top / (dv * dv) - drift_v_top / (2.0 * dv);
        Real b_v_top = -hp.xi * hp.xi * v_top / (dv * dv);
        Real c_v_top = diff_v_top / (dv * dv) + drift_v_top / (2.0 * dv);
        Real spec_v = std::abs(b_v_top)
                      + 2.0 * std::sqrt(std::abs(a_v_top * c_v_top));
        const Real C_MAX = 4.0;
        Size n_time_actual = config_.n_time;
        if (dt * spec_v > C_MAX) {
            n_time_actual = static_cast<Size>(std::ceil(
                static_cast<Real>(config_.n_time) * dt * spec_v / C_MAX));
            dt = hp.T / static_cast<Real>(n_time_actual);
        }

        auto scheme = create_scheme();
        scheme->reset();

        for (Size step = 0; step < n_time_actual; ++step) {
            apply_boundary_conditions(V, grid, hp,
                                       static_cast<Real>(step) * dt);
            if (config_.is_american) {
                // Ikonen-Toivanen (2004): project after each implicit sub-step
                // (Y1, Y2, ytilde1, ytilde2, V_new) — not just at the end.
                // This reduces splitting error from O(dt) large-constant to
                // O(dt) small-constant, significantly improving accuracy.
                scheme->step_american(V, V_new, dt, grid, coeffs, payoff_at_x);
            } else {
                scheme->step(V, V_new, dt, grid, coeffs);
            }
            apply_boundary_conditions(V_new, grid, hp,
                                       static_cast<Real>(step + 1) * dt);
            // Final safety projection for American (redundant with step_american
            // but guards against boundary-condition interference)
            if (config_.is_american) {
                apply_early_exercise(V_new, grid, hp);
            }
            V.swap(V_new);
        }

        last_n_time_used_ = n_time_actual;
        return interpolate_at(V, grid, x0, hp.v0);
    }

    const PDEEngine2DConfig& config() const { return config_; }

    // Returns the actual number of time steps used in the last price() call.
    // May exceed config_.n_time if adaptive time stepping was triggered.
    Size last_n_time_used() const { return last_n_time_used_; }

private:
    PDEEngine2DConfig config_;
    mutable Size last_n_time_used_ = 0;

    std::unique_ptr<ADISchemeBase> create_scheme() const {
        auto make_inner = [](ADISchemeType type, Real theta)
                -> std::unique_ptr<ADISchemeBase> {
            switch (type) {
                case ADISchemeType::HundsdorferVerwer:
                    return std::make_unique<HundsdorferVerwerScheme>(theta);
                case ADISchemeType::ModifiedCraigSneyd:
                    return std::make_unique<ModifiedCraigSneydScheme>(theta);
                case ADISchemeType::CraigSneyd:
                default:
                    return std::make_unique<CraigSneydScheme>(theta);
            }
        };

        if (config_.n_rannacher_warmup > 0) {
            // Warmup: theta=1.0 (L-stable, strong damping of high-freq oscillations)
            // Main: configured theta (typically 0.5 for second-order accuracy)
            auto warmup = make_inner(config_.scheme, 1.0);
            auto main = make_inner(config_.scheme, config_.theta);
            return std::make_unique<RannacherSmoothing2D>(
                std::move(warmup), std::move(main), config_.n_rannacher_warmup);
        }
        return make_inner(config_.scheme, config_.theta);
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

        // v_max: zero-flux Neumann boundary (dV/dv = 0)
        // Use constant extrapolation V[nv-1] = V[nv-2] instead of linear
        // extrapolation. Linear extrapolation (2*V[nv-2]-V[nv-3]) amplifies
        // boundary errors on fine grids where dv is small, causing numerical
        // blow-up. Constant extrapolation is more diffusive but stable.
        // Reference: in 't Hout & Foulon (2010) §3.3 recommend Neumann at v_max.
        if (nv >= 2) {
            for (Size i = 0; i < nx; ++i) {
                V[grid.idx(i, nv - 1)] = V[grid.idx(i, nv - 2)];
            }
        }
    }

    // American option early-exercise constraint (L1 projection / PSOR projection).
    // Enforces V(x, v) >= payoff(S(x)) at all interior grid points.
    // The payoff depends only on S = K*exp(x), not on v, so the projection
    // is applied uniformly across all v-rows for each x-column.
    //
    // This is the simplest operator-splitting scheme for American options
    // (Ikonen & Toivanen 2004). After each ADI time step, project the solution
    // onto the constraint set {V >= payoff}. The scheme is O(dt) accurate
    // due to the splitting, but converges to the true American price as dt->0.
    //
    // Note: Boundary points are also projected to ensure consistency
    // (e.g., deep ITM put at S->0 has payoff = K, V should be >= K*exp(-r*tau)).
    void apply_early_exercise(std::vector<Real>& V, const FDMGrid2D& grid,
                               const HestonPDEParams& hp) const {
        Size nx = grid.n_x(), nv = grid.n_v();
        Real K = hp.K;
        for (Size i = 0; i < nx; ++i) {
            Real S = K * std::exp(grid.x(i));
            Real payoff = config_.is_call
                ? std::max(S - K, 0.0)
                : std::max(K - S, 0.0);
            for (Size j = 0; j < nv; ++j) {
                Size k = grid.idx(i, j);
                if (V[k] < payoff) V[k] = payoff;
            }
        }
    }

    Real interpolate_at(const std::vector<Real>& V, const FDMGrid2D& grid,
                         Real x0, Real v0) const {
        // Hybrid interpolation: linear in x, quadratic (3-point Lagrange) in v.
        // Uses general Lagrange formula (works for both uniform and non-uniform grids).
        // Falls back to bilinear (linear in v) when v0 is at the v_min boundary
        // (j=0) where 3-point stencil is unavailable.
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
        Real wx = (x0 - x_l) / (x_r - x_l);

        // Helper: linear interpolation in x at a given v-row j_row
        auto interp_x = [&](Size j_row) -> Real {
            return V[grid.idx(i, j_row)] * (1.0 - wx)
                 + V[grid.idx(i + 1, j_row)] * wx;
        };

        // Quadratic Lagrange in v when 3-point stencil is available
        // (j >= 1 ensures we can use points j-1, j, j+1)
        if (j >= 1) {
            Real v_jm1 = grid.v(j - 1);
            Real v_j   = grid.v(j);
            Real v_jp1 = grid.v(j + 1);
            // General 3-point Lagrange (works for uniform and non-uniform)
            Real L0 = (v0 - v_j) * (v0 - v_jp1) / ((v_jm1 - v_j) * (v_jm1 - v_jp1));
            Real L1 = (v0 - v_jm1) * (v0 - v_jp1) / ((v_j - v_jm1) * (v_j - v_jp1));
            Real L2 = (v0 - v_jm1) * (v0 - v_j) / ((v_jp1 - v_jm1) * (v_jp1 - v_j));
            Real V_jm1 = interp_x(j - 1);
            Real V_j   = interp_x(j);
            Real V_jp1 = interp_x(j + 1);
            return V_jm1 * L0 + V_j * L1 + V_jp1 * L2;
        }

        // Fallback: bilinear (linear in v) at v_min boundary
        Real V_b = interp_x(j);
        Real V_t = interp_x(j + 1);
        Real v_b = grid.v(j), v_t = grid.v(j + 1);
        Real wv = (v0 - v_b) / (v_t - v_b);
        return V_b * (1.0 - wv) + V_t * wv;
    }
};

}  // namespace v1
}  // namespace cpphub
