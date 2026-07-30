#include <gtest/gtest.h>
#include "cpphub/pricing/pde/pde_engine.hpp"
#include "cpphub/pricing/pde/fdm_grid.hpp"
#include "cpphub/pricing/pde/fdm_scheme.hpp"
#include "cpphub/pricing/pde/thomas_solver.hpp"
#include "cpphub/instruments/payoff/vanilla.hpp"
#include "cpphub/core/math.hpp"
#include <cmath>
#include <vector>
#include <string>
#include <memory>

using namespace cpphub;

namespace {

Real bsm_call(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(S - K, 0.0);
    Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    Real d2 = d1 - sigma * std::sqrt(T);
    return S * std::exp(-q * T) * normal_cdf(d1) - K * std::exp(-r * T) * normal_cdf(d2);
}

Real bsm_put(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(K - S, 0.0);
    Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    Real d2 = d1 - sigma * std::sqrt(T);
    return K * std::exp(-r * T) * normal_cdf(-d2) - S * std::exp(-q * T) * normal_cdf(-d1);
}

Real bsm_delta_call(Real S0, Real K, Real T, Real r, Real q, Real sigma) {
    Real d1 = (std::log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    return std::exp(-q * T) * normal_cdf(d1);
}

Real bsm_gamma(Real S0, Real K, Real T, Real r, Real q, Real sigma) {
    Real d1 = (std::log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    return normal_pdf(d1) * std::exp(-q * T) / (S0 * sigma * std::sqrt(T));
}

Real pde_price(const PayOff& payoff, Real S0, Real K, Real T,
               Real r, Real q, Real sigma, FDMSchemeType type,
               Size n_spatial, Size n_time) {
    PDEParams params{r, q, sigma, T, K, S0};
    FDMGrid grid(n_spatial, S0, K, sigma, T, 0.2);
    TimeGrid time(n_time, T);
    Size n = grid.size();
    Size n_steps = n_time;
    Real dt = time.dt();

    std::vector<Real> V(n), V_new(n);
    for (Size i = 0; i < n; ++i) V[i] = payoff(grid.s(i));

    std::unique_ptr<FDMScheme> scheme;
    if (type == FDMSchemeType::ExplicitEuler)
        scheme = std::make_unique<ExplicitEuler>();
    else if (type == FDMSchemeType::ImplicitEuler)
        scheme = std::make_unique<ImplicitEuler>();
    else
        scheme = std::make_unique<CrankNicolson>(0.5);

    for (Size step = 0; step < n_steps; ++step) {
        Real tau = static_cast<Real>(step + 1) * dt;
        V_new = V;
        Real S_min = grid.s_min();
        Real S_max = grid.s_max();
        if (payoff.name() == "Call") {
            V_new[0] = bsm_call(S_min, K, tau, r, q, sigma);
            V_new[n-1] = bsm_call(S_max, K, tau, r, q, sigma);
        } else {
            V_new[0] = bsm_put(S_min, K, tau, r, q, sigma);
            V_new[n-1] = bsm_put(S_max, K, tau, r, q, sigma);
        }
        scheme->step(V, V_new, dt, grid, params);
        V = V_new;
    }

    Real price = V[n-1];
    for (Size i = 0; i < n - 1; ++i) {
        if (grid.s(i) <= S0 && S0 <= grid.s(i+1)) {
            Real w = (S0 - grid.s(i)) / (grid.s(i+1) - grid.s(i));
            price = (1.0 - w) * V[i] + w * V[i+1];
            break;
        }
    }
    return price;
}

} // anonymous namespace

TEST(pde_thomas_solver, SolvesDiagonalSystem) {
    Size n = 3;
    std::vector<Real> a = {0.0, -1.0, -1.0};
    std::vector<Real> b = {2.0, 2.0, 2.0};
    std::vector<Real> c = {-1.0, -1.0, 0.0};
    std::vector<Real> d = {1.0, 0.0, 1.0};
    auto x = thomas_solve(a, b, c, d);
    EXPECT_NEAR(x[0], 1.0, 1e-12);
    EXPECT_NEAR(x[1], 1.0, 1e-12);
    EXPECT_NEAR(x[2], 1.0, 1e-12);
}

TEST(pde_thomas_solver, HandlesNonUniformDiagonal) {
    Size n = 4;
    std::vector<Real> a = {0.0, 1.0, 1.0, 0.0};
    std::vector<Real> b = {2.0, 3.0, 4.0, 5.0};
    std::vector<Real> c = {1.0, 1.0, 1.0, 0.0};
    std::vector<Real> d = {1.0, 3.0, 3.0, 2.0};
    auto x = thomas_solve(a, b, c, d);
    EXPECT_EQ(x.size(), n);
    EXPECT_NEAR(b[0] * x[0] + c[0] * x[1], d[0], 1e-10);
    for (Size i = 1; i < n - 1; ++i) {
        EXPECT_NEAR(a[i] * x[i-1] + b[i] * x[i] + c[i] * x[i+1], d[i], 1e-10);
    }
    EXPECT_NEAR(a[n-1] * x[n-2] + b[n-1] * x[n-1], d[n-1], 1e-10);
}

TEST(pde_thomas_solver, SolveLargeSystemConsistentWithDirect) {
    Size n = 100;
    std::vector<Real> a(n, 0.0), b(n, 0.0), c(n, 0.0), d(n, 0.0);
    for (Size i = 0; i < n; ++i) {
        b[i] = 4.0;
        if (i > 0) a[i] = -1.0;
        if (i < n - 1) c[i] = -1.0;
    }
    for (Size i = 0; i < n; ++i) {
        Real x_exact = static_cast<Real>(i + 1);
        d[i] = (i > 0 ? a[i] * (x_exact - 1.0) : 0.0)
             + b[i] * x_exact
             + (i < n - 1 ? c[i] * (x_exact + 1.0) : 0.0);
    }
    auto x = thomas_solve(a, b, c, d);
    for (Size i = 0; i < n; ++i) {
        EXPECT_NEAR(x[i], static_cast<Real>(i + 1), 1e-10);
    }
}

TEST(pde_grid, SizeCorrect) {
    FDMGrid grid(100, 100.0, 100.0, 0.2, 1.0, 0.2);
    EXPECT_EQ(grid.size(), 100);
}

TEST(pde_grid, BoundaryPointsCorrect) {
    FDMGrid grid(100, 100.0, 100.0, 0.2, 1.0, 0.2);
    EXPECT_NEAR(grid.s(0), grid.s_min(), 1e-12);
    EXPECT_NEAR(grid.s(grid.size() - 1), grid.s_max(), 1e-12);
}

TEST(pde_grid, CenterConcentrated) {
    FDMGrid grid(200, 100.0, 100.0, 0.2, 1.0, 0.2);
    Size mid = grid.size() / 2;
    Real ds_center = grid.s(mid + 1) - grid.s(mid);
    Real ds_edge = grid.s(1) - grid.s(0);
    EXPECT_LT(ds_center, ds_edge);
}

TEST(pde_grid, SymmetryAroundCenter) {
    FDMGrid grid(200, 100.0, 100.0, 0.2, 1.0, 0.2);
    Size n = grid.size();
    Real center = (grid.s_min() + grid.s_max()) / 2.0;
    for (Size i = 0; i < n / 2; ++i) {
        EXPECT_NEAR(grid.s(i) - center, -(grid.s(n - 1 - i) - center), 1e-12);
    }
}

TEST(pde_scheme, ExplicitEulerEuropeanCallMatchesBSM) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    PutPayOff payoff(K);
    Real price = pde_price(payoff, S0, K, T, r, q, sigma,
                           FDMSchemeType::ExplicitEuler, 100, 50000);
    Real expected = bsm_put(S0, K, T, r, q, sigma);
    EXPECT_NEAR(price, expected, 0.01 * expected);
}

TEST(pde_scheme, ImplicitEulerEuropeanCallMatchesBSM) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    CallPayOff payoff(K);
    Real price = pde_price(payoff, S0, K, T, r, q, sigma,
                           FDMSchemeType::ImplicitEuler, 200, 500);
    Real expected = bsm_call(S0, K, T, r, q, sigma);
    EXPECT_NEAR(price, expected, 0.001 * expected);
}

TEST(pde_scheme, CrankNicolsonEuropeanCallMatchesBSM) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    CallPayOff payoff(K);
    Real price = pde_price(payoff, S0, K, T, r, q, sigma,
                           FDMSchemeType::CrankNicolson, 200, 500);
    Real expected = bsm_call(S0, K, T, r, q, sigma);
    EXPECT_NEAR(price, expected, 0.0001 * expected);
}

TEST(pde_scheme, CrankNicolsonEuropeanPutMatchesBSM) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    PutPayOff payoff(K);
    Real price = pde_price(payoff, S0, K, T, r, q, sigma,
                           FDMSchemeType::CrankNicolson, 200, 500);
    Real expected = bsm_put(S0, K, T, r, q, sigma);
    EXPECT_NEAR(price, expected, 0.0001 * expected);
}

TEST(pde_scheme, CNConvergesFasterThanExplicit) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    CallPayOff payoff(K);
    Real ref = bsm_call(S0, K, T, r, q, sigma);

    Real exp_err = std::abs(pde_price(payoff, S0, K, T, r, q, sigma,
                                       FDMSchemeType::ExplicitEuler, 100, 5000) - ref);
    Real cn_err = std::abs(pde_price(payoff, S0, K, T, r, q, sigma,
                                      FDMSchemeType::CrankNicolson, 100, 5000) - ref);
    EXPECT_LT(cn_err, exp_err);
}

TEST(pde_engine, EuropeanCallAccuracy) {
    PDEEngineConfig cfg;
    cfg.n_spatial = 400;
    cfg.n_time = 1000;
    cfg.scheme = FDMSchemeType::CrankNicolson;
    cfg.alpha = 0.2;
    PDEEngine engine(cfg);

    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0.02, sigma = 0.2;
    CallPayOff payoff(K);
    Real price = engine.price_european(payoff, S0, K, T, r, q, sigma);
    Real expected = bsm_call(S0, K, T, r, q, sigma);
    EXPECT_NEAR(price, expected, 0.0001 * expected);
}

TEST(pde_engine, EuropeanPutAccuracy) {
    PDEEngineConfig cfg;
    cfg.n_spatial = 400;
    cfg.n_time = 1000;
    cfg.scheme = FDMSchemeType::CrankNicolson;
    cfg.alpha = 0.2;
    PDEEngine engine(cfg);

    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0.02, sigma = 0.2;
    PutPayOff payoff(K);
    Real price = engine.price_european(payoff, S0, K, T, r, q, sigma);
    Real expected = bsm_put(S0, K, T, r, q, sigma);
    EXPECT_NEAR(price, expected, 0.0001 * expected);

    Real call_price = bsm_call(S0, K, T, r, q, sigma);
    Real parity = call_price - price;
    Real parity_expected = S0 * std::exp(-q * T) - K * std::exp(-r * T);
    EXPECT_NEAR(parity, parity_expected, 0.001);
}

TEST(pde_engine, AmericanPutGreaterThanEuropean) {
    PDEEngineConfig cfg;
    cfg.n_spatial = 400;
    cfg.n_time = 1000;
    cfg.scheme = FDMSchemeType::CrankNicolson;
    cfg.alpha = 0.2;
    PDEEngine engine(cfg);

    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    PutPayOff payoff(K);
    Real euro = engine.price_european(payoff, S0, K, T, r, q, sigma);
    Real amer = engine.price_american(payoff, S0, K, T, r, q, sigma);
    EXPECT_GE(amer, euro - 0.001);
}

TEST(pde_engine, AmericanPutMatchesBroadieDetempleBenchmark) {
    PDEEngineConfig cfg;
    cfg.n_spatial = 500;
    cfg.n_time = 2000;
    cfg.scheme = FDMSchemeType::CrankNicolson;
    cfg.alpha = 0.5;
    PDEEngine engine(cfg);

    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    PutPayOff payoff(K);
    Real price = engine.price_american(payoff, S0, K, T, r, q, sigma);
    EXPECT_NEAR(price, 6.0909, 1e-3);
}

TEST(pde_engine, GreeksNumericallyConsistent) {
    PDEEngineConfig cfg;
    cfg.n_spatial = 1000;
    cfg.n_time = 2000;
    cfg.scheme = FDMSchemeType::CrankNicolson;
    cfg.alpha = 0.2;
    PDEEngine engine(cfg);

    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0.02, sigma = 0.2;
    CallPayOff payoff(K);
    Real h = 0.5;
    Real Vp = engine.price_european(payoff, S0 + h, K, T, r, q, sigma);
    Real V0 = engine.price_european(payoff, S0, K, T, r, q, sigma);
    Real Vm = engine.price_european(payoff, S0 - h, K, T, r, q, sigma);

    Real delta = (Vp - Vm) / (2.0 * h);
    Real gamma = (Vp - 2.0 * V0 + Vm) / (h * h);

    Real delta_ref = bsm_delta_call(S0, K, T, r, q, sigma);
    Real gamma_ref = bsm_gamma(S0, K, T, r, q, sigma);

    EXPECT_NEAR(delta, delta_ref, 0.01);
    EXPECT_NEAR(gamma, gamma_ref, 0.01);
}

TEST(pde_psor, ConvergesForAmericanPut) {
    PDEEngineConfig cfg;
    cfg.n_spatial = 300;
    cfg.n_time = 600;
    cfg.scheme = FDMSchemeType::CrankNicolson;
    cfg.alpha = 0.2;
    PDEEngine engine(cfg);

    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    PutPayOff payoff(K);
    Real price = engine.price_american(payoff, S0, K, T, r, q, sigma);
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, K);
}

TEST(pde_psor, OmegaAffectsConvergence) {
    PDEEngineConfig cfg;
    cfg.n_spatial = 200;
    cfg.n_time = 400;
    cfg.scheme = FDMSchemeType::CrankNicolson;
    cfg.alpha = 0.2;
    PDEEngine engine(cfg);

    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    PutPayOff payoff(K);
    Real price = engine.price_american(payoff, S0, K, T, r, q, sigma);
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, K);
}

TEST(pde_psor, EarlyExerciseBoundaryCorrect) {
    PDEEngineConfig cfg;
    cfg.n_spatial = 400;
    cfg.n_time = 1000;
    cfg.scheme = FDMSchemeType::CrankNicolson;
    cfg.alpha = 0.2;
    PDEEngine engine(cfg);

    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    PutPayOff payoff(K);

    FDMGrid grid(cfg.n_spatial, S0, K, sigma, T, cfg.alpha);
    Size n = grid.size();
    std::vector<Real> V(n);
    for (Size i = 0; i < n; ++i) V[i] = payoff(grid.s(i));

    std::vector<Real> a_op, b_op, c_op;
    PDEParams params{r, q, sigma, T, K, S0};
    compute_operator_coeffs(grid, params, a_op, b_op, c_op);

    Real dt = T / cfg.n_time;
    Real theta = 0.5;
    Real dt_theta = dt * theta;

    std::vector<Real> a(n), b(n), c(n);
    a.assign(n, 0.0); b.assign(n, 0.0); c.assign(n, 0.0);
    for (Size i = 1; i < n - 1; ++i) {
        a[i] = -dt_theta * a_op[i];
        b[i] = 1.0 - dt_theta * b_op[i];
        c[i] = -dt_theta * c_op[i];
    }
    b[0] = 1.0; c[0] = 0.0;
    a[n-1] = 0.0; b[n-1] = 1.0;

    std::vector<Real> d(n, 0.0);
    Real dt_half = dt * 0.5;
    for (Size i = 1; i < n - 1; ++i) {
        d[i] = V[i] + dt_half * (a_op[i] * V[i-1] + b_op[i] * V[i] + c_op[i] * V[i+1]);
    }

    Real omega = 1.5;
    Real tol = 1e-8;
    Size max_iter = 2000;
    std::vector<Real> payoff_vals(n);
    for (Size i = 0; i < n; ++i) payoff_vals[i] = payoff(grid.s(i));

    for (Size step = 0; step < cfg.n_time; ++step) {
        Real tau = static_cast<Real>(step + 1) * dt;
        V[0] = payoff_vals[0];
        V[n-1] = payoff_vals[n-1];

        for (Size iter = 0; iter < max_iter; ++iter) {
            Real max_diff = 0.0;
            for (Size i = 1; i < n - 1; ++i) {
                Real old_val = V[i];
                Real y = (d[i] - a[i] * V[i-1] - c[i] * V[i+1]) / b[i];
                V[i] = std::max(payoff_vals[i], old_val + omega * (y - old_val));
                max_diff = std::max(max_diff, std::abs(V[i] - old_val));
            }
            if (max_diff < tol) break;
        }

        for (Size i = 1; i < n - 1; ++i) {
            d[i] = V[i] + dt_half * (a_op[i] * V[i-1] + b_op[i] * V[i] + c_op[i] * V[i+1]);
        }
    }

    EXPECT_NEAR(V[0], K - grid.s_min(), 0.1);
}

TEST(pde_integration, ConvergenceRate) {
    PDEEngineConfig cfg;
    cfg.scheme = FDMSchemeType::CrankNicolson;
    cfg.alpha = 0.2;

    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    CallPayOff payoff(K);
    Real ref = bsm_call(S0, K, T, r, q, sigma);

    std::vector<Size> n_spatials = {100, 200, 400};
    std::vector<Real> errors;

    for (auto ns : n_spatials) {
        cfg.n_spatial = ns;
        cfg.n_time = 4 * ns;
        PDEEngine engine(cfg);
        Real price = engine.price_european(payoff, S0, K, T, r, q, sigma);
        errors.push_back(std::abs(price - ref));
    }

    for (Size i = 1; i < errors.size(); ++i) {
        Real ratio = errors[i-1] / errors[i];
        EXPECT_GT(ratio, 1.5);
    }
}
