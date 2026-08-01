// Chebyshev collocation spectral method tests for 1D BSM PDE
//
// The Chebyshev collocation method discretizes the spatial domain using
// Chebyshev-Gauss-Lobatto points and approximates derivatives via the
// Chebyshev differentiation matrix. This yields spectral (exponential)
// convergence for smooth solutions, far superior to FD methods.
//
// For BSM PDE in log-price x = ln(S/K):
//   dV/dtau = 0.5*sigma^2 * d2V/dx2 + (r-q-0.5*sigma^2) * dV/dx - r*V
// with V(x, 0) = payoff(K*exp(x)) and far-field boundary conditions.
//
// IMPORTANT: The payoff kink at S=K limits spectral convergence to
// algebraic (not exponential) for vanilla options. The domain half-width
// x_range must be chosen carefully: too large => poor resolution at S=K
// (Chebyshev points are sparsest at the domain center); too small =>
// boundary condition inaccuracy. A good rule of thumb is
// x_range ≈ 5-10 * sigma * sqrt(T).
//
// Reference: Boyd (2001) "Chebyshev and Fourier Spectral Methods",
//            Trefethen (2000) "Spectral Methods in MATLAB".
#include <gtest/gtest.h>
#include "cpphub/pricing/pde/chebyshev_pde.hpp"
#include "cpphub/core/math.hpp"
#include <cmath>

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

}  // anonymous namespace

// T2.2: Smoke test — Chebyshev PDE produces positive finite price
TEST(ChebyshevPDE, ProducesPositivePrice) {
    ChebyshevPDEConfig cfg;
    cfg.n_points = 64;
    cfg.n_time = 200;
    cfg.x_range = 5.0;
    cfg.is_call = true;
    ChebyshevPDEEngine engine(cfg);

    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Real price = engine.price(S0, K, T, r, q, sigma);
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, S0);
    EXPECT_TRUE(std::isfinite(price));
}

// T2.2: European call matches BSM analytic within 2%
// N=80, x_range=2.0 (10*sigma*sqrt(T)) gives ~0.9% relative error.
// The payoff kink at S=K limits convergence; 2% is a realistic bar.
TEST(ChebyshevPDE, EuropeanCallMatchesBSM) {
    ChebyshevPDEConfig cfg;
    cfg.n_points = 80;
    cfg.n_time = 300;
    cfg.x_range = 2.0;  // 10 * sigma * sqrt(T) for sigma=0.2, T=1
    cfg.is_call = true;
    ChebyshevPDEEngine engine(cfg);

    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Real price = engine.price(S0, K, T, r, q, sigma);
    Real analytic = bsm_call(S0, K, T, r, q, sigma);
    Real rel_err = std::abs(price - analytic) / analytic;
    EXPECT_LT(rel_err, 0.02) << "price=" << price << " BSM=" << analytic;
}

// T2.2: European put matches BSM analytic within 2%
// Put has slightly larger error than call due to the S->0 boundary
// (Dirichlet V=K*exp(-r*tau) is approximate at finite S_min).
TEST(ChebyshevPDE, EuropeanPutMatchesBSM) {
    ChebyshevPDEConfig cfg;
    cfg.n_points = 80;
    cfg.n_time = 300;
    cfg.x_range = 2.0;
    cfg.is_call = false;
    ChebyshevPDEEngine engine(cfg);

    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Real price = engine.price(S0, K, T, r, q, sigma);
    Real analytic = bsm_put(S0, K, T, r, q, sigma);
    Real rel_err = std::abs(price - analytic) / analytic;
    EXPECT_LT(rel_err, 0.02) << "price=" << price << " BSM=" << analytic;
}

// T2.2: Multiple strikes match BSM within 2%
TEST(ChebyshevPDE, MultipleStrikesMatchBSM) {
    ChebyshevPDEConfig cfg;
    cfg.n_points = 80;
    cfg.n_time = 300;
    cfg.x_range = 2.0;
    cfg.is_call = true;
    ChebyshevPDEEngine engine(cfg);

    Real S0 = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    for (Real K : {80.0, 90.0, 100.0, 110.0, 120.0}) {
        Real price = engine.price(S0, K, T, r, q, sigma);
        Real analytic = bsm_call(S0, K, T, r, q, sigma);
        Real rel_err = std::abs(price - analytic) / analytic;
        EXPECT_LT(rel_err, 0.02) << "K=" << K << " price=" << price << " BSM=" << analytic;
    }
}

// T2.2: Spectral convergence — error decreases as N increases.
// With x_range=2.0 and fine time grid, the spatial error dominates and
// decreases algebraically (not exponentially, due to payoff kink).
// N=64 should be significantly better than N=24.
TEST(ChebyshevPDE, SpectralConvergenceWithN) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Real analytic = bsm_call(S0, K, T, r, q, sigma);

    std::vector<Size> Ns = {24, 40, 64, 96};
    std::vector<Real> errors;
    for (Size N : Ns) {
        ChebyshevPDEConfig cfg;
        cfg.n_points = N;
        cfg.n_time = 400;  // fine time grid to isolate spatial error
        cfg.x_range = 2.0;
        cfg.is_call = true;
        ChebyshevPDEEngine engine(cfg);
        Real price = engine.price(S0, K, T, r, q, sigma);
        Real err = std::abs(price - analytic);
        errors.push_back(err);
    }

    // All errors should be reasonable (< 1.5 for N>=24 with x_range=2)
    for (Size i = 0; i < errors.size(); ++i) {
        EXPECT_LT(errors[i], 1.50) << "N=" << Ns[i] << " err=" << errors[i];
    }
    // Monotonic decrease: spectral convergence property
    for (Size i = 1; i < errors.size(); ++i) {
        EXPECT_LT(errors[i], errors[i - 1])
            << "Spectral convergence violated: N=" << Ns[i]
            << " err=" << errors[i] << " >= N=" << Ns[i - 1]
            << " err=" << errors[i - 1];
    }
    // N=96 should achieve < 0.10 absolute error
    EXPECT_LT(errors[3], 0.10)
        << "N=96 err=" << errors[3] << " expected < 0.10";
}

// T2.2: Put-call parity holds
// call - put = S0*exp(-q*T) - K*exp(-r*T)
TEST(ChebyshevPDE, PutCallParityHolds) {
    ChebyshevPDEConfig cfg;
    cfg.n_points = 80;
    cfg.n_time = 300;
    cfg.x_range = 2.0;
    cfg.is_call = true;
    ChebyshevPDEEngine engine_call(cfg);
    cfg.is_call = false;
    ChebyshevPDEEngine engine_put(cfg);

    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.02, sigma = 0.2;
    Real call = engine_call.price(S0, K, T, r, q, sigma);
    Real put = engine_put.price(S0, K, T, r, q, sigma);

    Real parity = call - put;
    Real parity_expected = S0 * std::exp(-q * T) - K * std::exp(-r * T);
    EXPECT_NEAR(parity, parity_expected, 0.05);
}

// T2.2: Works with different maturities
// x_range is scaled with sqrt(T) to maintain consistent resolution
// relative to the diffusion scale sigma*sqrt(T).
TEST(ChebyshevPDE, DifferentMaturitiesMatchBSM) {
    Real S0 = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.2;
    for (Real T : {0.25, 0.5, 1.0, 2.0, 5.0}) {
        ChebyshevPDEConfig cfg;
        cfg.n_points = 80;
        cfg.n_time = static_cast<Size>(300 * T);
        cfg.x_range = 2.0 * std::sqrt(T);  // scale with diffusion length
        cfg.is_call = true;
        ChebyshevPDEEngine eng(cfg);
        Real price = eng.price(S0, K, T, r, q, sigma);
        Real analytic = bsm_call(S0, K, T, r, q, sigma);
        Real rel_err = std::abs(price - analytic) / analytic;
        EXPECT_LT(rel_err, 0.02) << "T=" << T << " price=" << price << " BSM=" << analytic;
    }
}
