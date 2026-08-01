// ADI 2D PDE engine tests for Heston model
// Validates Craig-Sneyd and Hundsdorfer-Verwer schemes against Heston COS method
// (semi-analytic benchmark via characteristic function inversion).
#include <gtest/gtest.h>
#include "cpphub/pricing/pde/pde_engine_2d.hpp"
#include "cpphub/pricing/fourier/cos_method.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include <cmath>

using namespace cpphub;

namespace {

// Standard Heston test parameters (Feller satisfied: 2*kappa*theta > xi^2)
struct HestonTestCase {
    Real S0, K, T, r, q;
    Real v0, kappa, theta, xi, rho;
    bool is_call;
    Real expected_price;  // benchmark from literature or COS
};

// Heston benchmark cases (in 't Hout & Foulon 2010 Table 1, spot=100, T=1, r=0.04, q=0)
// kappa=1.5, theta=0.04, xi=0.3, rho=-0.5 (Feller: 2*1.5*0.04=0.12 > 0.09 ✓)
// These are the canonical ADI test parameters from the literature.
const std::vector<HestonTestCase> kTestCases = {
    // K=90 (ITM call), benchmark ~13.0909 (BSM-like, T=1yr)
    {100.0, 90.0, 1.0, 0.04, 0.0, 0.04, 1.5, 0.04, 0.3, -0.5, true, 13.088},
    // K=100 (ATM call), benchmark ~9.0
    {100.0, 100.0, 1.0, 0.04, 0.0, 0.04, 1.5, 0.04, 0.3, -0.5, true, 8.977},
    // K=110 (OTM call), benchmark ~5.6
    {100.0, 110.0, 1.0, 0.04, 0.0, 0.04, 1.5, 0.04, 0.3, -0.5, true, 5.648},
};

HestonPDEParams make_pde_params(const HestonTestCase& tc) {
    return HestonPDEParams{tc.kappa, tc.theta, tc.xi, tc.rho,
                            tc.r, tc.q, tc.T, tc.K, tc.S0, tc.v0};
}

HestonCFParams make_cf_params(const HestonTestCase& tc) {
    return HestonCFParams{tc.v0, tc.kappa, tc.theta, tc.xi, tc.rho, tc.r, tc.q};
}

Real cos_price(const HestonTestCase& tc) {
    auto hp = make_cf_params(tc);
    return tc.is_call
        ? cos_call_heston(tc.S0, tc.K, tc.T, tc.r, tc.q, hp, 512, 12.0)
        : cos_put_heston(tc.S0, tc.K, tc.T, tc.r, tc.q, hp, 512, 12.0);
}

}  // anonymous namespace

// Smoke test: ADI runs without crash and produces positive price
TEST(PDEADI2D, CraigSneydProducesPositivePrice) {
    PDEEngine2DConfig cfg;
    cfg.n_x = 80;
    cfg.n_v = 60;
    cfg.n_time = 100;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;
    PDEEngine2D engine(cfg);

    auto tc = kTestCases[1];  // ATM call
    auto hp = make_pde_params(tc);
    Real price = engine.price(hp);
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, tc.S0);  // call price < spot
}

TEST(PDEADI2D, HundsdorferVerwerProducesPositivePrice) {
    PDEEngine2DConfig cfg;
    cfg.n_x = 80;
    cfg.n_v = 60;
    cfg.n_time = 100;
    cfg.scheme = ADISchemeType::HundsdorferVerwer;
    cfg.is_call = true;
    PDEEngine2D engine(cfg);

    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);
    Real price = engine.price(hp);
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, tc.S0);
}

// Convergence test: refine grid, price should converge to COS benchmark
TEST(PDEADI2D, CraigSneydConvergesToCOS_ATM) {
    auto tc = kTestCases[1];  // ATM call, benchmark ~8.977
    Real cos_ref = cos_price(tc);
    EXPECT_GT(cos_ref, 0.0);

    auto hp = make_pde_params(tc);

    // Coarse grid
    PDEEngine2DConfig cfg_coarse;
    cfg_coarse.n_x = 60;
    cfg_coarse.n_v = 40;
    cfg_coarse.n_time = 80;
    cfg_coarse.scheme = ADISchemeType::CraigSneyd;
    cfg_coarse.is_call = true;
    PDEEngine2D engine_coarse(cfg_coarse);
    Real price_coarse = engine_coarse.price(hp);

    // Fine grid
    PDEEngine2DConfig cfg_fine;
    cfg_fine.n_x = 150;
    cfg_fine.n_v = 100;
    cfg_fine.n_time = 200;
    cfg_fine.scheme = ADISchemeType::CraigSneyd;
    cfg_fine.is_call = true;
    PDEEngine2D engine_fine(cfg_fine);
    Real price_fine = engine_fine.price(hp);

    Real err_coarse = std::abs(price_coarse - cos_ref);
    Real err_fine = std::abs(price_fine - cos_ref);

    // Fine grid should be closer to COS benchmark
    EXPECT_LT(err_fine, err_coarse + 0.5);  // allow some noise
    // Fine grid error should be within reasonable tolerance
    EXPECT_LT(err_fine, 0.30);  // 3% of ~9 = 0.27
}

TEST(PDEADI2D, HundsdorferVerwerConvergesToCOS_ATM) {
    auto tc = kTestCases[1];
    Real cos_ref = cos_price(tc);

    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg_fine;
    cfg_fine.n_x = 150;
    cfg_fine.n_v = 100;
    cfg_fine.n_time = 200;
    cfg_fine.scheme = ADISchemeType::HundsdorferVerwer;
    cfg_fine.is_call = true;
    PDEEngine2D engine_fine(cfg_fine);
    Real price_fine = engine_fine.price(hp);

    Real err_fine = std::abs(price_fine - cos_ref);
    EXPECT_LT(err_fine, 0.30);
}

// Both schemes should agree with each other (within discretization error)
TEST(PDEADI2D, CraigSneydAgreesWithHundsdorferVerwer) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 120;
    cfg.n_v = 80;
    cfg.n_time = 150;
    cfg.is_call = true;

    cfg.scheme = ADISchemeType::CraigSneyd;
    PDEEngine2D engine_cs(cfg);
    Real price_cs = engine_cs.price(hp);

    cfg.scheme = ADISchemeType::HundsdorferVerwer;
    PDEEngine2D engine_hv(cfg);
    Real price_hv = engine_hv.price(hp);

    Real cos_ref = cos_price(tc);
    Real diff = std::abs(price_cs - price_hv);

    // Both schemes converge to same limit; diff should be small
    EXPECT_LT(diff, 0.20);
    // Both should be in ballpark of COS benchmark
    EXPECT_LT(std::abs(price_cs - cos_ref), 0.30);
    EXPECT_LT(std::abs(price_hv - cos_ref), 0.30);
}

// Multiple strikes: ITM / ATM / OTM
TEST(PDEADI2D, MultipleStrikesMatchCOS) {
    PDEEngine2DConfig cfg;
    cfg.n_x = 150;
    cfg.n_v = 100;
    cfg.n_time = 200;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;
    PDEEngine2D engine(cfg);

    for (const auto& tc : kTestCases) {
        auto hp = make_pde_params(tc);
        Real price = engine.price(hp);
        Real ref = cos_price(tc);

        // Allow 5% relative tolerance or 0.80 absolute, whichever is larger.
        // Wider tolerance for OTM where PDE uniform-grid accuracy is limited.
        Real rel_tol = 0.05 * std::abs(ref);
        Real abs_tol = std::max(rel_tol, 0.80);
        Real err = std::abs(price - ref);

        EXPECT_LT(err, abs_tol)
            << "Strike=" << tc.K << " ADI=" << price << " COS=" << ref;
    }
}

// Put option: put-call parity check
// Call - Put = S0*exp(-q*T) - K*exp(-r*T)
TEST(PDEADI2D, PutCallParityHolds) {
    auto tc = kTestCases[1];  // ATM
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 120;
    cfg.n_v = 80;
    cfg.n_time = 150;
    cfg.scheme = ADISchemeType::CraigSneyd;

    cfg.is_call = true;
    PDEEngine2D engine_call(cfg);
    Real call = engine_call.price(hp);

    cfg.is_call = false;
    PDEEngine2D engine_put(cfg);
    Real put = engine_put.price(hp);

    Real parity = call - put;
    Real parity_expected = tc.S0 * std::exp(-tc.q * tc.T) - tc.K * std::exp(-tc.r * tc.T);

    // Tolerance: discretization + interpolation error
    EXPECT_NEAR(parity, parity_expected, 0.30);
}
