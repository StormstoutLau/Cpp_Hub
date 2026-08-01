// RBF-FD (Radial Basis Function - Finite Difference) tests for multi-asset BSM PDE
//
// Test strategy (TDD, ascending complexity):
//   1. Smoke test — engine returns finite positive price
//   2. 1D BSM limit — d=1, European call/put vs analytic BSM (rel err < 5%)
//   3. 2D Margrabe exchange option — max(S1 - S2, 0), analytic formula
//   4. 2D Basket option — equal-weight basket call vs Monte Carlo
//
// References:
//   Margrabe (1978) "The Value of an Option to Exchange One Asset for Another"
//   Wright & Fornberg (2016) on RBF-FD stability and convergence
#include <gtest/gtest.h>
#include "cpphub/pricing/pde/rbf_fd.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_payoffs.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/rng.hpp"
#include <cmath>
#include <vector>
#include <numeric>

using namespace cpphub;

namespace {

// ============ BSM analytic (1D) ============
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

// ============ Margrabe (1978) exchange option analytic ============
// Payoff: max(S1(T) - S2(T), 0)
// V = S1*exp(-q1*T)*N(d1) - S2*exp(-q2*T)*N(d2)
// sigma_eff^2 = sigma1^2 + sigma2^2 - 2*rho*sigma1*sigma2
// d1 = (ln(S1/S2) + (q2 - q1 + 0.5*sigma_eff^2)*T) / (sigma_eff*sqrt(T))
// d2 = d1 - sigma_eff*sqrt(T)
Real margrabe_exchange(Real S1, Real S2,
                        Real q1, Real q2,
                        Real sigma1, Real sigma2,
                        Real rho, Real T) {
    Real sigma_eff = std::sqrt(sigma1 * sigma1 + sigma2 * sigma2
                                - 2.0 * rho * sigma1 * sigma2);
    if (sigma_eff < 1e-12) {
        // Degenerate: deterministic spread, exercise if S1 > S2 at T (in expectation)
        // Just return discounted intrinsic at-expectation
        Real fwd1 = S1 * std::exp((0.0 - q1) * T);
        Real fwd2 = S2 * std::exp((0.0 - q2) * T);
        return std::max(fwd1 - fwd2, 0.0) * std::exp(-0.0 * T);
    }
    Real d1 = (std::log(S1 / S2) + (q2 - q1 + 0.5 * sigma_eff * sigma_eff) * T)
              / (sigma_eff * std::sqrt(T));
    Real d2 = d1 - sigma_eff * std::sqrt(T);
    Real r_eff = 0.0;  // Margrabe uses r=0 in d1/d2 (the discount is on q1, q2)
    // Note: Standard Margrabe (no drift in d1 beyond q2-q1) is consistent with
    // the change-of-numeraire derivation; here we use r=0 implicitly since the
    // exchange option is invariant to interest rate.
    return S1 * std::exp(-q1 * T) * normal_cdf(d1)
         - S2 * std::exp(-q2 * T) * normal_cdf(d2);
}

// ============ Monte Carlo basket pricer (2 assets) ============
// Returns discounted average payoff: exp(-r*T) * E[max(w1*S1 + w2*S2 - K, 0)]
Real mc_basket_call(Real S1, Real S2,
                     Real sigma1, Real sigma2, Real rho,
                     Real r, Real q1, Real q2, Real T, Real K,
                     Real w1, Real w2,
                     Size n_paths, unsigned long seed) {
    MultiAssetGBMConfig cfg;
    cfg.S0 = {S1, S2};
    cfg.sigma = {sigma1, sigma2};
    cfg.q = {q1, q2};
    cfg.r = r;
    cfg.T = T;
    cfg.n_steps = 1;  // terminal payoff only
    cfg.correlation = {{1.0, rho}, {rho, 1.0}};
    MultiAssetGBMPathGenerator gen(cfg);

    Philox4x64 rng(seed, 0);
    Real sum = 0.0;
    Real sum2 = 0.0;
    for (Size p = 0; p < n_paths; ++p) {
        auto paths = gen.generate_path(rng);
        Real basket = w1 * paths[0].back() + w2 * paths[1].back();
        Real payoff = std::max(basket - K, 0.0);
        sum += payoff;
        sum2 += payoff * payoff;
    }
    Real mean = sum / static_cast<Real>(n_paths);
    Real variance = (sum2 - sum * sum / static_cast<Real>(n_paths))
                    / static_cast<Real>(n_paths - 1);
    Real std_err = std::sqrt(variance / static_cast<Real>(n_paths));
    (void)std_err;  // not asserted tightly; MC just a cross-check
    return std::exp(-r * T) * mean;
}

}  // anonymous namespace

// ============================================================
// 1. Smoke test: engine returns finite positive price
// ============================================================
TEST(RBFFDEngine, SmokeTestProducesFinitePositivePrice) {
    RBFFDConfig cfg;
    cfg.n_per_dim = 15;
    cfg.stencil_size = 12;
    cfg.n_time = 50;
    cfg.theta = 0.5;
    cfg.domain_width = 5.0;
    RBFFDEngine engine(cfg);

    MultiAssetBSMParams params;
    params.S0 = {100.0};
    params.sigma = {0.2};
    params.q = {0.0};
    params.corr = {{1.0}};
    params.r = 0.05;
    params.T = 1.0;
    params.K = {100.0};

    auto payoff = [](const std::vector<Real>& S) -> Real {
        return std::max(S[0] - 100.0, 0.0);
    };

    Real price = engine.price_european(payoff, params);
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, 200.0);
    EXPECT_TRUE(std::isfinite(price));
}

// ============================================================
// 2. 1D BSM limit: European call matches BSM analytic
// ============================================================
TEST(RBFFDEngine, OneDimEuropeanCallMatchesBSM) {
    RBFFDConfig cfg;
    cfg.n_per_dim = 31;
    cfg.stencil_size = 15;
    cfg.n_time = 100;
    cfg.theta = 0.5;
    cfg.domain_width = 5.0;
    cfg.rbf = RBFType::Polyharmonic5;
    cfg.poly_degree = 2;
    RBFFDEngine engine(cfg);

    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    MultiAssetBSMParams params;
    params.S0 = {S0};
    params.sigma = {sigma};
    params.q = {q};
    params.corr = {{1.0}};
    params.r = r;
    params.T = T;
    params.K = {K};

    auto payoff = [K](const std::vector<Real>& S) -> Real {
        return std::max(S[0] - K, 0.0);
    };

    Real price = engine.price_european(payoff, params);
    Real analytic = bsm_call(S0, K, T, r, q, sigma);
    Real rel_err = std::abs(price - analytic) / analytic;
    // RBF-FD with PHS r^5 + poly_degree=2 should achieve < 5% on 1D BSM
    EXPECT_LT(rel_err, 0.05) << "RBF-FD=" << price << " BSM=" << analytic
                              << " rel_err=" << rel_err;
}

// ============================================================
// 3. 1D BSM put matches BSM analytic
// ============================================================
TEST(RBFFDEngine, OneDimEuropeanPutMatchesBSM) {
    RBFFDConfig cfg;
    cfg.n_per_dim = 31;
    cfg.stencil_size = 15;
    cfg.n_time = 100;
    cfg.theta = 0.5;
    cfg.domain_width = 5.0;
    RBFFDEngine engine(cfg);

    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    MultiAssetBSMParams params;
    params.S0 = {S0};
    params.sigma = {sigma};
    params.q = {q};
    params.corr = {{1.0}};
    params.r = r;
    params.T = T;
    params.K = {K};

    auto payoff = [K](const std::vector<Real>& S) -> Real {
        return std::max(K - S[0], 0.0);
    };

    Real price = engine.price_european(payoff, params);
    Real analytic = bsm_put(S0, K, T, r, q, sigma);
    Real rel_err = std::abs(price - analytic) / analytic;
    EXPECT_LT(rel_err, 0.05) << "RBF-FD=" << price << " BSM=" << analytic
                              << " rel_err=" << rel_err;
}

// ============================================================
// 4. 1D BSM with dividend yield
// ============================================================
TEST(RBFFDEngine, OneDimEuropeanCallWithDividend) {
    RBFFDConfig cfg;
    cfg.n_per_dim = 31;
    cfg.stencil_size = 15;
    cfg.n_time = 100;
    cfg.theta = 0.5;
    cfg.domain_width = 5.0;
    RBFFDEngine engine(cfg);

    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.03, sigma = 0.25;
    MultiAssetBSMParams params;
    params.S0 = {S0};
    params.sigma = {sigma};
    params.q = {q};
    params.corr = {{1.0}};
    params.r = r;
    params.T = T;
    params.K = {K};

    auto payoff = [K](const std::vector<Real>& S) -> Real {
        return std::max(S[0] - K, 0.0);
    };

    Real price = engine.price_european(payoff, params);
    Real analytic = bsm_call(S0, K, T, r, q, sigma);
    Real rel_err = std::abs(price - analytic) / analytic;
    EXPECT_LT(rel_err, 0.05) << "RBF-FD=" << price << " BSM=" << analytic
                              << " rel_err=" << rel_err;
}

// ============================================================
// 5. 2D Margrabe exchange option (analytic cross-check)
// ============================================================
TEST(RBFFDEngine, TwoDimMargrabeExchangeOption) {
    RBFFDConfig cfg;
    cfg.n_per_dim = 21;
    cfg.stencil_size = 18;
    cfg.n_time = 80;
    cfg.theta = 0.5;
    cfg.domain_width = 4.0;
    RBFFDEngine engine(cfg);

    Real S1 = 100.0, S2 = 100.0;
    Real sigma1 = 0.2, sigma2 = 0.25;
    Real rho = 0.5;
    Real r = 0.05, q1 = 0.0, q2 = 0.0, T = 1.0;

    MultiAssetBSMParams params;
    params.S0 = {S1, S2};
    params.sigma = {sigma1, sigma2};
    params.q = {q1, q2};
    params.corr = {{1.0, rho}, {rho, 1.0}};
    params.r = r;
    params.T = T;
    params.K = {S1, S2};  // log-S centered at S0 (ATM)

    auto payoff = [](const std::vector<Real>& S) -> Real {
        return std::max(S[0] - S[1], 0.0);
    };

    Real price = engine.price_european(payoff, params);
    Real analytic = margrabe_exchange(S1, S2, q1, q2, sigma1, sigma2, rho, T);
    Real rel_err = (analytic > 1e-6)
                    ? std::abs(price - analytic) / analytic
                    : std::abs(price - analytic);
    // 2D RBF-FD with PHS+poly should achieve < 8% on Margrabe (sparse grid)
    EXPECT_LT(rel_err, 0.08) << "RBF-FD=" << price << " Margrabe=" << analytic
                              << " rel_err=" << rel_err;
}

// ============================================================
// 6. 2D Margrabe with rho = 0 (uncorrelated assets)
// ============================================================
TEST(RBFFDEngine, TwoDimMargrabeUncorrelated) {
    RBFFDConfig cfg;
    cfg.n_per_dim = 21;
    cfg.stencil_size = 18;
    cfg.n_time = 80;
    cfg.theta = 0.5;
    cfg.domain_width = 4.0;
    RBFFDEngine engine(cfg);

    Real S1 = 100.0, S2 = 110.0;
    Real sigma1 = 0.2, sigma2 = 0.3;
    Real rho = 0.0;
    Real r = 0.05, q1 = 0.0, q2 = 0.0, T = 1.0;

    MultiAssetBSMParams params;
    params.S0 = {S1, S2};
    params.sigma = {sigma1, sigma2};
    params.q = {q1, q2};
    params.corr = {{1.0, rho}, {rho, 1.0}};
    params.r = r;
    params.T = T;
    params.K = {S1, S2};

    auto payoff = [](const std::vector<Real>& S) -> Real {
        return std::max(S[0] - S[1], 0.0);
    };

    Real price = engine.price_european(payoff, params);
    Real analytic = margrabe_exchange(S1, S2, q1, q2, sigma1, sigma2, rho, T);
    Real rel_err = (analytic > 1e-6)
                    ? std::abs(price - analytic) / analytic
                    : std::abs(price - analytic);
    EXPECT_LT(rel_err, 0.08) << "RBF-FD=" << price << " Margrabe=" << analytic
                              << " rel_err=" << rel_err;
}

// ============================================================
// 7. 2D Basket call vs Monte Carlo
// ============================================================
TEST(RBFFDEngine, TwoDimBasketCallVsMonteCarlo) {
    RBFFDConfig cfg;
    cfg.n_per_dim = 21;
    cfg.stencil_size = 18;
    cfg.n_time = 80;
    cfg.theta = 0.5;
    cfg.domain_width = 4.0;
    RBFFDEngine engine(cfg);

    Real S1 = 100.0, S2 = 100.0;
    Real sigma1 = 0.2, sigma2 = 0.25;
    Real rho = 0.3;
    Real r = 0.05, q1 = 0.0, q2 = 0.0, T = 1.0;
    Real K = 100.0;
    Real w1 = 0.5, w2 = 0.5;

    MultiAssetBSMParams params;
    params.S0 = {S1, S2};
    params.sigma = {sigma1, sigma2};
    params.q = {q1, q2};
    params.corr = {{1.0, rho}, {rho, 1.0}};
    params.r = r;
    params.T = T;
    params.K = {S1, S2};

    auto payoff = [K, w1, w2](const std::vector<Real>& S) -> Real {
        Real basket = w1 * S[0] + w2 * S[1];
        return std::max(basket - K, 0.0);
    };

    Real pde_price = engine.price_european(payoff, params);
    Real mc_price = mc_basket_call(S1, S2, sigma1, sigma2, rho,
                                     r, q1, q2, T, K, w1, w2,
                                     200000, 42);
    // Allow 5% relative tolerance (MC std err ~0.5%, PDE error ~3%)
    Real rel_diff = std::abs(pde_price - mc_price) / mc_price;
    EXPECT_LT(rel_diff, 0.05) << "PDE=" << pde_price << " MC=" << mc_price
                               << " rel_diff=" << rel_diff;
}

// ============================================================
// 8. Convergence: increasing n_per_dim improves accuracy
// ============================================================
TEST(RBFFDEngine, ConvergenceWithGridRefinement) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Real analytic = bsm_call(S0, K, T, r, q, sigma);

    MultiAssetBSMParams params;
    params.S0 = {S0};
    params.sigma = {sigma};
    params.q = {q};
    params.corr = {{1.0}};
    params.r = r;
    params.T = T;
    params.K = {K};

    auto payoff = [K](const std::vector<Real>& S) -> Real {
        return std::max(S[0] - K, 0.0);
    };

    // Coarse grid
    RBFFDConfig cfg_coarse;
    cfg_coarse.n_per_dim = 15;
    cfg_coarse.stencil_size = 10;
    cfg_coarse.n_time = 50;
    cfg_coarse.theta = 0.5;
    cfg_coarse.domain_width = 5.0;
    RBFFDEngine engine_coarse(cfg_coarse);
    Real err_coarse = std::abs(engine_coarse.price_european(payoff, params) - analytic);

    // Fine grid
    RBFFDConfig cfg_fine;
    cfg_fine.n_per_dim = 41;
    cfg_fine.stencil_size = 20;
    cfg_fine.n_time = 200;
    cfg_fine.theta = 0.5;
    cfg_fine.domain_width = 5.0;
    RBFFDEngine engine_fine(cfg_fine);
    Real err_fine = std::abs(engine_fine.price_european(payoff, params) - analytic);

    EXPECT_LT(err_fine, err_coarse) << "Fine grid err=" << err_fine
                                     << " >= Coarse grid err=" << err_coarse;
}

// ============================================================
// 9. Validation: invalid params throw
// ============================================================
TEST(RBFFDEngine, InvalidParamsThrow) {
    RBFFDConfig cfg;
    cfg.n_per_dim = 3;  // too small
    EXPECT_THROW(RBFFDEngine engine(cfg), std::invalid_argument);

    cfg.n_per_dim = 25;
    cfg.theta = 1.5;  // out of range
    EXPECT_THROW(RBFFDEngine engine(cfg), std::invalid_argument);

    cfg.theta = 0.5;
    cfg.poly_degree = 3;  // too high
    EXPECT_THROW(RBFFDEngine engine(cfg), std::invalid_argument);
}

// ============================================================
// 10. Validation: dimension mismatch throws
// ============================================================
TEST(RBFFDEngine, DimensionMismatchThrows) {
    RBFFDConfig cfg;
    cfg.n_per_dim = 15;
    cfg.stencil_size = 12;
    cfg.n_time = 50;
    RBFFDEngine engine(cfg);

    // sigma size mismatch
    MultiAssetBSMParams params;
    params.S0 = {100.0, 100.0};
    params.sigma = {0.2};  // wrong size
    params.q = {0.0, 0.0};
    params.corr = {{1.0, 0.0}, {0.0, 1.0}};
    params.r = 0.05;
    params.T = 1.0;

    auto payoff = [](const std::vector<Real>& S) -> Real {
        return std::max(S[0] - 100.0, 0.0);
    };
    EXPECT_THROW(engine.price_european(payoff, params), std::invalid_argument);
}
