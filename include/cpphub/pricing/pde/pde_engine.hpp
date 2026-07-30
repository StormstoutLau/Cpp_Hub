#pragma once
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <string>
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/pricing/pde/fdm_grid.hpp"
#include "cpphub/pricing/pde/fdm_scheme.hpp"
#include "cpphub/pricing/pde/thomas_solver.hpp"
#include "cpphub/instruments/payoff/payoff.hpp"

namespace cpphub {
inline namespace v1 {

struct PDEEngineConfig {
    Size n_spatial = 400;
    Size n_time = 1000;
    Real alpha = 0.2;
    FDMSchemeType scheme = FDMSchemeType::CrankNicolson;
    Real s_multiplier = 5.0;
};

class PDEEngine {
public:
    explicit PDEEngine(PDEEngineConfig config = PDEEngineConfig{})
        : config_(config) {}

    Real price_european(const PayOff& payoff, Real S0, Real K, Real T,
                         Real r, Real q, Real sigma) const {
        PDEParams params{r, q, sigma, T, K, S0};
        FDMGrid grid(config_.n_spatial, S0, K, sigma, T, config_.alpha);
        TimeGrid time(config_.n_time, T);
        Size n = grid.size();
        Size n_steps = config_.n_time;
        Real dt = time.dt();

        std::vector<Real> V(n), V_new(n);
        for (Size i = 0; i < n; ++i) {
            V[i] = payoff(grid.s(i));
        }

        auto scheme = create_scheme();

        for (Size step = 0; step < n_steps; ++step) {
            Real tau = static_cast<Real>(step + 1) * dt;
            V_new = V;
            set_boundary(V_new, tau, grid, payoff, r, q, K, sigma, false);
            scheme->step(V, V_new, dt, grid, params);
            V = V_new;
        }

        return interpolate_at_S0(V, grid, S0);
    }

    Real price_american(const PayOff& payoff, Real S0, Real K, Real T,
                         Real r, Real q, Real sigma) const {
        PDEParams params{r, q, sigma, T, K, S0};
        FDMGrid grid(config_.n_spatial, S0, K, sigma, T, config_.alpha);
        TimeGrid time(config_.n_time, T);
        Size n = grid.size();
        Size n_steps = config_.n_time;
        Real dt = time.dt();

        std::vector<Real> V(n), V_new(n);
        for (Size i = 0; i < n; ++i) {
            V[i] = payoff(grid.s(i));
        }

        Real omega = 1.5;
        Real tol = 1e-8;
        Size max_iter = 2000;

        for (Size step = 0; step < n_steps; ++step) {
            Real tau = static_cast<Real>(step + 1) * dt;

            std::vector<Real> payoff_vals(n);
            for (Size i = 0; i < n; ++i) {
                payoff_vals[i] = payoff(grid.s(i));
            }

            std::vector<Real> a, b, c, d;
            build_matrix(grid, params, dt, a, b, c);

            d.assign(n, 0.0);
            if (config_.scheme == FDMSchemeType::ImplicitEuler) {
                for (Size i = 1; i < n - 1; ++i) d[i] = V[i];
            } else {
                std::vector<Real> a_op, b_op, c_op;
                compute_operator_coeffs(grid, params, a_op, b_op, c_op);
                Real dt_half = dt * 0.5;
                for (Size i = 1; i < n - 1; ++i) {
                    d[i] = V[i] + dt_half * (a_op[i] * V[i - 1] +
                                              b_op[i] * V[i] +
                                              c_op[i] * V[i + 1]);
                }
            }

            V_new = V;
            set_boundary(V_new, tau, grid, payoff, r, q, K, sigma, true);
            psor_solve(V_new, payoff_vals, a, b, c, d, omega, tol, max_iter);
            V = V_new;
        }

        return interpolate_at_S0(V, grid, S0);
    }

    struct Greeks {
        Real delta;
        Real gamma;
        Real theta;
    };

    Greeks greeks(const PayOff& payoff, Real S0, Real K, Real T,
                   Real r, Real q, Real sigma, bool american = false) const {
        Real h = 0.01;
        Real V0 = price_either(payoff, S0, K, T, r, q, sigma, american);
        Real Vp = price_either(payoff, S0 + h, K, T, r, q, sigma, american);
        Real Vm = price_either(payoff, S0 - h, K, T, r, q, sigma, american);

        Real delta = (Vp - Vm) / (2.0 * h);
        Real gamma = (Vp - 2.0 * V0 + Vm) / (h * h);

        Real dt_small = 1.0 / 365.0;
        Real V_forward = price_either(payoff, S0, K, T + dt_small, r, q, sigma, american);
        Real theta = -(V_forward - V0) / dt_small;

        return {delta, gamma, theta};
    }

    std::string name() const { return "PDEEngine"; }
    const PDEEngineConfig& config() const { return config_; }

private:
    PDEEngineConfig config_;

    Real price_either(const PayOff& payoff, Real S0, Real K, Real T,
                       Real r, Real q, Real sigma, bool american) const {
        return american ? price_american(payoff, S0, K, T, r, q, sigma)
                        : price_european(payoff, S0, K, T, r, q, sigma);
    }

    static Real bsm_call(Real S, Real K, Real tau, Real r, Real q, Real sigma) {
        if (tau <= 0.0) return std::max(S - K, 0.0);
        Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * tau) / (sigma * std::sqrt(tau));
        Real d2 = d1 - sigma * std::sqrt(tau);
        return S * std::exp(-q * tau) * normal_cdf(d1) - K * std::exp(-r * tau) * normal_cdf(d2);
    }

    static Real bsm_put(Real S, Real K, Real tau, Real r, Real q, Real sigma) {
        if (tau <= 0.0) return std::max(K - S, 0.0);
        Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * tau) / (sigma * std::sqrt(tau));
        Real d2 = d1 - sigma * std::sqrt(tau);
        return K * std::exp(-r * tau) * normal_cdf(-d2) - S * std::exp(-q * tau) * normal_cdf(-d1);
    }

    void set_boundary(std::vector<Real>& V, Real tau, const FDMGrid& grid,
                       const PayOff& payoff, Real r, Real q, Real K, Real sigma,
                       bool) const {
        Size n = V.size();
        Real S_min = grid.s_min();
        Real S_max = grid.s_max();
        if (payoff.name() == "Call") {
            V[0] = bsm_call(S_min, K, tau, r, q, sigma);
            V[n - 1] = bsm_call(S_max, K, tau, r, q, sigma);
        } else if (payoff.name() == "Put") {
            V[0] = bsm_put(S_min, K, tau, r, q, sigma);
            V[n - 1] = bsm_put(S_max, K, tau, r, q, sigma);
        } else {
            Real df = std::exp(-r * tau);
            V[0] = payoff(S_min) * df;
            V[n - 1] = payoff(S_max) * df;
        }
    }

    void psor_solve(std::vector<Real>& V, const std::vector<Real>& payoff,
                     const std::vector<Real>& a, const std::vector<Real>& b,
                     const std::vector<Real>& c, const std::vector<Real>& d,
                     Real omega, Real tol, Size max_iter) const {
        Size n = V.size();
        for (Size iter = 0; iter < max_iter; ++iter) {
            Real max_diff = 0.0;
            for (Size i = 1; i < n - 1; ++i) {
                Real old_val = V[i];
                Real y = (d[i] - a[i] * V[i - 1] - c[i] * V[i + 1]) / b[i];
                V[i] = std::max(payoff[i], old_val + omega * (y - old_val));
                max_diff = std::max(max_diff, std::abs(V[i] - old_val));
            }
            if (max_diff < tol) break;
        }
    }

    void build_matrix(const FDMGrid& grid, const PDEParams& params, Real dt,
                       std::vector<Real>& a, std::vector<Real>& b,
                       std::vector<Real>& c) const {
        Size n = grid.size();
        std::vector<Real> a_op, b_op, c_op;
        compute_operator_coeffs(grid, params, a_op, b_op, c_op);

        a.assign(n, 0.0); b.assign(n, 0.0); c.assign(n, 0.0);

        Real theta = (config_.scheme == FDMSchemeType::ImplicitEuler) ? 1.0 : 0.5;
        Real dt_theta = dt * theta;

        for (Size i = 1; i < n - 1; ++i) {
            a[i] = -dt_theta * a_op[i];
            b[i] = 1.0 - dt_theta * b_op[i];
            c[i] = -dt_theta * c_op[i];
        }

        b[0] = 1.0; c[0] = 0.0;
        a[n - 1] = 0.0; b[n - 1] = 1.0;
    }

    std::unique_ptr<FDMScheme> create_scheme() const {
        switch (config_.scheme) {
            case FDMSchemeType::ExplicitEuler:
                return std::make_unique<ExplicitEuler>();
            case FDMSchemeType::ImplicitEuler:
                return std::make_unique<ImplicitEuler>();
            case FDMSchemeType::CrankNicolson:
            default:
                return std::make_unique<CrankNicolson>(0.5);
        }
    }

    Real interpolate_at_S0(const std::vector<Real>& V, const FDMGrid& grid, Real S0) const {
        Size n = grid.size();
        for (Size i = 0; i < n - 1; ++i) {
            if (grid.s(i) <= S0 && S0 <= grid.s(i + 1)) {
                Real w = (S0 - grid.s(i)) / (grid.s(i + 1) - grid.s(i));
                return (1.0 - w) * V[i] + w * V[i + 1];
            }
        }
        return V[n - 1];
    }
};

}  // namespace v1
}  // namespace cpphub
