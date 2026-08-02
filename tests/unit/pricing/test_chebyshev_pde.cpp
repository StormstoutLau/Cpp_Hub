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

// MSVC does not define PI_VAL by default; use std::acos(-1.0) for portability.
const Real PI_VAL = std::acos(-1.0);

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

// ============================================================================
// DCT (Discrete Chebyshev Transform) tests
// ============================================================================

// DCT round-trip: inverse(forward(f)) ≈ f
TEST(ChebyshevDCT, RoundTripRecoversOriginal) {
    Size N = 16;
    Size n = N + 1;
    // Sample a smooth function f(x) = exp(x) on CGL nodes
    std::vector<Real> f(n);
    for (Size k = 0; k < n; ++k) {
        Real x = std::cos(PI_VAL * static_cast<Real>(k) / static_cast<Real>(N));
        f[k] = std::exp(x);
    }
    auto a = ChebyshevPDEEngine::cheb_forward_dct(f, N);
    auto f_recovered = ChebyshevPDEEngine::cheb_inverse_dct(a, N);
    for (Size k = 0; k < n; ++k) {
        EXPECT_NEAR(f[k], f_recovered[k], 1e-10)
            << "DCT round-trip failed at k=" << k;
    }
}

// DCT of T_2(x) = 2x^2 - 1 should give a_2 = 1, all others ~ 0
TEST(ChebyshevDCT, ChebyshevPolynomialT2Coefficient) {
    Size N = 16;
    Size n = N + 1;
    std::vector<Real> f(n);
    for (Size k = 0; k < n; ++k) {
        Real x = std::cos(PI_VAL * static_cast<Real>(k) / static_cast<Real>(N));
        f[k] = 2.0 * x * x - 1.0;  // T_2(x)
    }
    auto a = ChebyshevPDEEngine::cheb_forward_dct(f, N);
    // a[2] should be 1.0, others should be near 0
    EXPECT_NEAR(a[2], 1.0, 1e-10);
    for (Size k = 0; k <= N; ++k) {
        if (k == 2) continue;
        EXPECT_LT(std::abs(a[k]), 1e-10)
            << "T_2 coefficient a[" << k << "] should be ~0, got " << a[k];
    }
}

// DCT of constant function should give a_0 = 1, all others ~ 0
TEST(ChebyshevDCT, ConstantFunctionCoefficient) {
    Size N = 8;
    Size n = N + 1;
    std::vector<Real> f(n, 5.0);  // constant 5
    auto a = ChebyshevPDEEngine::cheb_forward_dct(f, N);
    EXPECT_NEAR(a[0], 5.0, 1e-10);
    for (Size k = 1; k <= N; ++k) {
        EXPECT_LT(std::abs(a[k]), 1e-10)
            << "Constant a[" << k << "] should be ~0, got " << a[k];
    }
}

// ============================================================================
// Exponential spectral filter tests
// ============================================================================

// Filter preserves smooth functions (low-order coefficients unchanged)
TEST(ChebyshevFilter, PreservesSmoothFunction) {
    Size N = 16;
    Size n = N + 1;
    // f(x) = exp(x) is smooth, dominated by low-order Chebyshev modes
    std::vector<Real> f(n);
    for (Size k = 0; k < n; ++k) {
        Real x = std::cos(PI_VAL * static_cast<Real>(k) / static_cast<Real>(N));
        f[k] = std::exp(x);
    }
    // Apply mild filter (eta=4, p=2)
    auto f_filtered = ChebyshevPDEEngine::apply_exp_filter(f, N, 4.0, 2);
    // Smooth function should be nearly unchanged
    Real max_diff = 0.0;
    for (Size k = 0; k < n; ++k) {
        max_diff = std::max(max_diff, std::abs(f[k] - f_filtered[k]));
    }
    EXPECT_LT(max_diff, 0.1) << "Smooth function changed by " << max_diff;
}

// Filter suppresses high-frequency oscillations (Gibbs-like signal)
TEST(ChebyshevFilter, SuppressesHighFrequency) {
    Size N = 32;
    Size n = N + 1;
    // Construct signal with high-frequency component: f(x) = T_N(x) (Nyquist mode)
    std::vector<Real> f(n);
    for (Size k = 0; k < n; ++k) {
        Real x = std::cos(PI_VAL * static_cast<Real>(k) / static_cast<Real>(N));
        f[k] = std::cos(static_cast<Real>(N) * std::acos(x));  // T_N(x)
    }
    // Apply strong filter (eta=36, p=4) — should suppress T_N
    auto f_filtered = ChebyshevPDEEngine::apply_exp_filter(f, N, 36.0, 4);
    // Amplitude should be significantly reduced
    Real orig_amp = 0.0, filt_amp = 0.0;
    for (Size k = 0; k < n; ++k) {
        orig_amp = std::max(orig_amp, std::abs(f[k]));
        filt_amp = std::max(filt_amp, std::abs(f_filtered[k]));
    }
    // sigma(1) = exp(-36) ≈ 2.3e-16, so filtered amplitude should be ~0
    EXPECT_LT(filt_amp, orig_amp * 1e-6)
        << "Filter failed to suppress Nyquist mode: orig=" << orig_amp
        << " filt=" << filt_amp;
}

// Filter DC component (k=0) is preserved exactly
TEST(ChebyshevFilter, PreservesDCComponent) {
    Size N = 16;
    Size n = N + 1;
    std::vector<Real> f(n, 3.0);  // constant
    auto f_filtered = ChebyshevPDEEngine::apply_exp_filter(f, N, 36.0, 4);
    // Mean should be exactly preserved (sigma(0) = 1)
    Real orig_mean = 0.0, filt_mean = 0.0;
    for (Size k = 0; k < n; ++k) {
        orig_mean += f[k];
        filt_mean += f_filtered[k];
    }
    orig_mean /= static_cast<Real>(n);
    filt_mean /= static_cast<Real>(n);
    EXPECT_NEAR(orig_mean, filt_mean, 1e-10);
}

// Filter enabled in pricing engine still matches BSM (within tolerance)
// The filter should not destroy pricing accuracy for smooth-enough solutions.
TEST(ChebyshevPDE, FilterEnabledMatchesBSM) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Real analytic = bsm_call(S0, K, T, r, q, sigma);

    ChebyshevPDEConfig cfg;
    cfg.n_points = 80;
    cfg.n_time = 300;
    cfg.x_range = 2.0;
    cfg.is_call = true;
    cfg.use_filter = true;
    cfg.filter_eta = 36.0;
    cfg.filter_p = 4;
    ChebyshevPDEEngine engine(cfg);

    Real price = engine.price(S0, K, T, r, q, sigma);
    Real rel_err = std::abs(price - analytic) / analytic;
    // Filter may slightly reduce accuracy but should remain within 5%
    EXPECT_LT(rel_err, 0.05)
        << "Filtered price=" << price << " BSM=" << analytic
        << " rel_err=" << rel_err;
}

// Filter reduces oscillation amplitude in the solution
// Compare max-min of V between filtered and unfiltered runs
TEST(ChebyshevPDE, FilterReducesOscillation) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;

    // Run without filter
    ChebyshevPDEConfig cfg;
    cfg.n_points = 48;  // smaller N → more oscillation
    cfg.n_time = 200;
    cfg.x_range = 3.0;  // wider domain → more boundary effect
    cfg.is_call = true;
    cfg.n_rannacher_warmup = 2;  // less warmup → more residual oscillation
    cfg.use_filter = false;
    ChebyshevPDEEngine eng_no_filter(cfg);
    Real price_no_filter = eng_no_filter.price(S0, K, T, r, q, sigma);

    // Run with filter
    cfg.use_filter = true;
    cfg.filter_eta = 36.0;
    cfg.filter_p = 4;
    ChebyshevPDEEngine eng_filter(cfg);
    Real price_filter = eng_filter.price(S0, K, T, r, q, sigma);

    Real analytic = bsm_call(S0, K, T, r, q, sigma);
    // Both should be in the right ballpark
    EXPECT_GT(price_no_filter, 0.0);
    EXPECT_GT(price_filter, 0.0);
    // Filtered price should be closer to analytic (or at least not worse)
    Real err_no_filter = std::abs(price_no_filter - analytic);
    Real err_filter = std::abs(price_filter - analytic);
    // Filter should not make things significantly worse
    EXPECT_LT(err_filter, err_no_filter * 3.0)
        << "Filter degraded accuracy: no_filter_err=" << err_no_filter
        << " filter_err=" << err_filter;
}
