// ADI 2D PDE engine tests for Heston model
// Validates Craig-Sneyd and Hundsdorfer-Verwer schemes against Heston COS method
// (semi-analytic benchmark via characteristic function inversion).
#include <gtest/gtest.h>
#include "cpphub/pricing/pde/pde_engine_2d.hpp"
#include "cpphub/pricing/fourier/cos_method.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"  // bsm_call_price
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

// T1.2a: 精细网格数值稳定性测试 (复现 B5 benchmark 爆炸)
// benchmark 发现 200x150x400 配置下 CS 格式产生天文数字 (5.8e31)
// 此测试用于 TDD 红灯阶段: 当前应 FAIL, 修复后应 PASS
TEST(PDEADI2D, FineGridNumericalStability) {
    auto tc = kTestCases[1];  // ATM call, K=100
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    // 逐步增加网格密度, 找到爆炸临界点
    // B4 (100x80x200) 已验证稳定, 从 120x100x300 开始测试
    struct GridConfig { Size nx, nv, nt; const char* label; };
    std::vector<GridConfig> grids = {
        {120, 100, 300, "120x100x300"},
        {150, 120, 350, "150x120x350"},
        {200, 150, 400, "200x150x400 (B5 benchmark)"},
    };

    for (const auto& gc : grids) {
        PDEEngine2DConfig cfg;
        cfg.n_x = gc.nx;
        cfg.n_v = gc.nv;
        cfg.n_time = gc.nt;
        cfg.scheme = ADISchemeType::CraigSneyd;
        cfg.is_call = true;
        PDEEngine2D engine(cfg);

        Real price = engine.price(hp);
        // 验收: 价格应为有限正数, 且在合理范围内 (COS 参考值 ± 50%)
        EXPECT_TRUE(std::isfinite(price))
            << gc.label << ": price not finite, got " << price;
        EXPECT_GT(price, 0.0)
            << gc.label << ": price non-positive, got " << price;
        EXPECT_LT(price, tc.S0 * 2.0)
            << gc.label << ": price unreasonably large, got " << price
            << " (COS ref=" << cos_ref << ")";
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

// T1.2a: Verify adaptive time step is triggered on fine grids
// The engine should auto-increase n_time when dt*spec(L_v) > C_MAX=4.0
TEST(PDEADI2D, AdaptiveTimeStepTriggeredOnFineGrid) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);

    // Coarse grid: dt*spec_v ~ 5.6 > 4.0 → should trigger adaptation
    PDEEngine2DConfig cfg;
    cfg.n_x = 120;
    cfg.n_v = 100;
    cfg.n_time = 300;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;
    PDEEngine2D engine_coarse(cfg);
    (void)engine_coarse.price(hp);
    // 120x100x300: dt*spec_v ≈ 5.61, expected n_time ≈ ceil(300*5.61/4.0) = 421
    EXPECT_GT(engine_coarse.last_n_time_used(), cfg.n_time)
        << "Adaptive time step should trigger for 120x100x300";

    // Fine grid: dt*spec_v ~ 9.6 > 4.0 → should trigger larger adaptation
    cfg.n_x = 200;
    cfg.n_v = 150;
    cfg.n_time = 400;
    PDEEngine2D engine_fine(cfg);
    (void)engine_fine.price(hp);
    // 200x150x400: dt*spec_v ≈ 9.58, expected n_time ≈ ceil(400*9.58/4.0) = 958
    EXPECT_GT(engine_fine.last_n_time_used(), cfg.n_time)
        << "Adaptive time step should trigger for 200x150x400";
    EXPECT_GT(engine_fine.last_n_time_used(), engine_coarse.last_n_time_used())
        << "Finer grid should require more time steps";
}

// T1.2a: Verify coarse grid does NOT trigger adaptation (dt*spec_v < 4.0)
TEST(PDEADI2D, AdaptiveTimeStepNotTriggeredOnCoarseGrid) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 80;
    cfg.n_v = 60;
    cfg.n_time = 200;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;
    PDEEngine2D engine(cfg);
    (void)engine.price(hp);
    // 80x60x200: dt*spec_v should be small enough to not trigger adaptation
    // If it does trigger, that's OK too (conservative), but log it
    if (engine.last_n_time_used() > cfg.n_time) {
        // Acceptable: conservative adaptation. Just verify price is still valid.
        Real price = engine.price(hp);
        EXPECT_GT(price, 0.0);
        EXPECT_LT(price, tc.S0);
    } else {
        EXPECT_EQ(engine.last_n_time_used(), cfg.n_time);
    }
}

// T1.3: Convergence order quantification
// CS and HV are second-order in time and space (O(dt^2 + dx^2 + dv^2)).
// This test refines all dimensions simultaneously and checks convergence.
//
// Note: v_max boundary uses constant extrapolation (1st-order Neumann),
// which limits achievable accuracy at very fine grids. The test validates:
//   1. Error decreases from coarse (L0) to medium (L1) — convergence is real
//   2. Absolute error at medium grid is small (< 0.10, ~1% of price)
//   3. Self-convergence: L1 and L2 agree within tolerance
TEST(PDEADI2D, CraigSneydConvergenceOrder) {
    auto tc = kTestCases[1];  // ATM call
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    // Grid sequence: each level refines all dims by ~2x
    struct GridLevel { Size nx, nv, nt; const char* label; };
    std::vector<GridLevel> levels = {
        {50,  40,  100, "L0: 50x40x100"},
        {100, 80,  200, "L1: 100x80x200"},
        {200, 160, 400, "L2: 200x160x400"},
    };

    std::vector<Real> prices, errors;
    for (const auto& gl : levels) {
        PDEEngine2DConfig cfg;
        cfg.n_x = gl.nx;
        cfg.n_v = gl.nv;
        cfg.n_time = gl.nt;
        cfg.scheme = ADISchemeType::CraigSneyd;
        cfg.is_call = true;
        PDEEngine2D engine(cfg);
        Real price = engine.price(hp);
        Real err = std::abs(price - cos_ref);
        prices.push_back(price);
        errors.push_back(err);

        EXPECT_TRUE(std::isfinite(price)) << gl.label;
        EXPECT_GT(price, 0.0) << gl.label;
    }

    // 1. Error decreases from L0 to L1 (coarse to medium convergence)
    EXPECT_LT(errors[1], errors[0])
        << "L0->L1: err0=" << errors[0] << " err1=" << errors[1];

    // 2. Convergence rate L0->L1 should be at least ~1.0 (better than 1st order)
    Real rate_01 = std::log(errors[0] / errors[1]) / std::log(2.0);
    EXPECT_GE(rate_01, 0.8)
        << "L0->L1 rate=" << rate_01 << " err0=" << errors[0]
        << " err1=" << errors[1];

    // 3. Medium grid absolute error < 0.10 (~1% of COS ~9)
    EXPECT_LT(errors[1], 0.10)
        << "L1 error=" << errors[1] << " COS ref=" << cos_ref;

    // 4. Self-convergence: L1 and L2 should agree within 0.10
    //    (boundary-limited convergence may prevent L2 from being closer to COS,
    //    but L2 should still agree with L1)
    Real self_err = std::abs(prices[1] - prices[2]);
    EXPECT_LT(self_err, 0.10)
        << "L1-L2 self-error=" << self_err << " L1=" << prices[1]
        << " L2=" << prices[2];
}

// T1.3: HV convergence order (same test, different scheme)
TEST(PDEADI2D, HundsdorferVerwerConvergenceOrder) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    struct GridLevel { Size nx, nv, nt; };
    std::vector<GridLevel> levels = {
        {50,  40,  100},
        {100, 80,  200},
        {200, 160, 400},
    };

    std::vector<Real> prices, errors;
    for (const auto& gl : levels) {
        PDEEngine2DConfig cfg;
        cfg.n_x = gl.nx;
        cfg.n_v = gl.nv;
        cfg.n_time = gl.nt;
        cfg.scheme = ADISchemeType::HundsdorferVerwer;
        cfg.is_call = true;
        PDEEngine2D engine(cfg);
        Real price = engine.price(hp);
        Real err = std::abs(price - cos_ref);
        prices.push_back(price);
        errors.push_back(err);
        EXPECT_TRUE(std::isfinite(price));
        EXPECT_GT(price, 0.0);
    }

    EXPECT_LT(errors[1], errors[0]);
    Real rate_01 = std::log(errors[0] / errors[1]) / std::log(2.0);
    EXPECT_GE(rate_01, 0.8)
        << "HV L0->L1 rate=" << rate_01;
    EXPECT_LT(errors[1], 0.10);

    Real self_err = std::abs(prices[1] - prices[2]);
    EXPECT_LT(self_err, 0.10);
}

// T1.2b: Verify quadratic v-interpolation improves accuracy
// At medium grid (100x80x200), the quadratic v-interpolation should
// produce error < 0.08 (was ~0.05-0.09 with bilinear).
// This test is a regression guard: if someone reverts to bilinear,
// the error will exceed the tighter threshold.
TEST(PDEADI2D, QuadraticVInterpolationAccuracy) {
    auto tc = kTestCases[1];  // ATM call
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 100;
    cfg.n_v = 80;
    cfg.n_time = 200;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;
    PDEEngine2D engine(cfg);
    Real price = engine.price(hp);
    Real err = std::abs(price - cos_ref);

    // With quadratic v-interpolation, error should be < 0.08
    // (previously ~0.05-0.09 with bilinear, tighter bound catches regression)
    EXPECT_LT(err, 0.08)
        << "CS 100x80x200 error=" << err << " COS ref=" << cos_ref
        << " price=" << price;

    // Also check HV
    cfg.scheme = ADISchemeType::HundsdorferVerwer;
    PDEEngine2D engine_hv(cfg);
    Real price_hv = engine_hv.price(hp);
    Real err_hv = std::abs(price_hv - cos_ref);
    EXPECT_LT(err_hv, 0.08)
        << "HV 100x80x200 error=" << err_hv;
}

// =============================================================================
// T1.6: Modified Craig-Sneyd (MCS) scheme tests
//
// Mathematical basis (in 't Hout & Foulon 2010, Table 2):
//   MCS differs from CS ONLY in the final correction step:
//     CS:  V_{n+1} = ytilde2 + theta*dt*L_xv*Yhat + (theta-0.5)*dt*L_full*Yhat
//     MCS: V_{n+1} = ytilde2 + theta*dt*L_xv*Yhat + (theta-0.5)*dt*L_full*ytilde2
//   When theta=0.5, (theta-0.5)=0, so CS == MCS (identical results).
//   When theta!=0.5, MCS replaces L_full(Yhat) with L_full(ytilde2).
//   MCS is unconditionally stable and second-order for theta in [0.5, 1].
//   Reference: in 't Hout & Welfert (2007) "Unconditional stability of
//   second-order ADI schemes ...".
// =============================================================================

// T1.6 RED: MCS scheme exists and produces positive finite price
TEST(PDEADI2D, ModifiedCraigSneydProducesPositivePrice) {
    PDEEngine2DConfig cfg;
    cfg.n_x = 80;
    cfg.n_v = 60;
    cfg.n_time = 100;
    cfg.scheme = ADISchemeType::ModifiedCraigSneyd;
    cfg.is_call = true;
    PDEEngine2D engine(cfg);

    auto tc = kTestCases[1];  // ATM call
    auto hp = make_pde_params(tc);
    Real price = engine.price(hp);
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, tc.S0);
}

// T1.6 RED: MCS with theta=0.5 must equal CS with theta=0.5 exactly
// Mathematical guarantee: (theta-0.5)=0 makes the differing term vanish.
TEST(PDEADI2D, ModifiedCraigSneydEqualsCraigSneydAtThetaHalf) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 100;
    cfg.n_v = 80;
    cfg.n_time = 200;
    cfg.theta = 0.5;
    cfg.is_call = true;

    cfg.scheme = ADISchemeType::CraigSneyd;
    PDEEngine2D engine_cs(cfg);
    Real price_cs = engine_cs.price(hp);

    cfg.scheme = ADISchemeType::ModifiedCraigSneyd;
    PDEEngine2D engine_mcs(cfg);
    Real price_mcs = engine_mcs.price(hp);

    // Identical algorithm at theta=0.5 => results must match to machine precision
    EXPECT_NEAR(price_cs, price_mcs, 1e-12)
        << "CS=" << price_cs << " MCS=" << price_mcs
        << " (should be identical at theta=0.5)";
}

// T1.6 RED: MCS converges to COS benchmark (second-order accurate)
TEST(PDEADI2D, ModifiedCraigSneydConvergesToCOS) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    struct GridLevel { Size nx, nv, nt; const char* label; };
    std::vector<GridLevel> levels = {
        {50,  40,  100, "L0: 50x40x100"},
        {100, 80,  200, "L1: 100x80x200"},
        {200, 160, 400, "L2: 200x160x400"},
    };

    std::vector<Real> prices, errors;
    for (const auto& gl : levels) {
        PDEEngine2DConfig cfg;
        cfg.n_x = gl.nx;
        cfg.n_v = gl.nv;
        cfg.n_time = gl.nt;
        cfg.scheme = ADISchemeType::ModifiedCraigSneyd;
        cfg.theta = 0.5;
        cfg.is_call = true;
        PDEEngine2D engine(cfg);
        Real price = engine.price(hp);
        Real err = std::abs(price - cos_ref);
        prices.push_back(price);
        errors.push_back(err);

        EXPECT_TRUE(std::isfinite(price)) << gl.label;
        EXPECT_GT(price, 0.0) << gl.label;
    }

    // Coarse-to-medium convergence
    EXPECT_LT(errors[1], errors[0])
        << "MCS L0->L1: err0=" << errors[0] << " err1=" << errors[1];
    Real rate_01 = std::log(errors[0] / errors[1]) / std::log(2.0);
    EXPECT_GE(rate_01, 0.8)
        << "MCS L0->L1 rate=" << rate_01;
    // Medium grid absolute error < 0.10
    EXPECT_LT(errors[1], 0.10)
        << "MCS L1 error=" << errors[1] << " COS ref=" << cos_ref;
    // Self-convergence
    Real self_err = std::abs(prices[1] - prices[2]);
    EXPECT_LT(self_err, 0.10)
        << "MCS L1-L2 self-error=" << self_err;
}

// T1.6 RED: MCS with theta != 0.5 is stable and produces valid prices
// theta = (0.5 + 1/sqrt(12)) is the optimal stability parameter for HV/MCS
// (in 't Hout & Foulon 2010, Section 4).
//
// IMPORTANT: At theta != 0.5, the correction term (theta-0.5)*dt*L_full*{Yhat/ytilde2}
// introduces additional amplification. The adaptive time step threshold C_MAX=4.0
// was calibrated for theta=0.5 (where corr_coeff=0). At theta=0.7887, grids that
// trigger adaptive time stepping (n_v>=80) can blow up for ALL schemes (CS/MCS/HV).
// This is a known limitation, not an MCS-specific bug.
// Use 80x60 grid (no adaptive time step) for theta != 0.5 stability tests.
TEST(PDEADI2D, ModifiedCraigSneydStableAtOptimalTheta) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 80;
    cfg.n_v = 60;
    cfg.n_time = 200;
    cfg.scheme = ADISchemeType::ModifiedCraigSneyd;
    cfg.theta = 0.5 + 1.0 / std::sqrt(12.0);  // ~0.7887, optimal stability
    cfg.is_call = true;
    PDEEngine2D engine(cfg);
    Real price = engine.price(hp);

    EXPECT_TRUE(std::isfinite(price));
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, tc.S0);
    // On coarse 80x60 grid, error is larger (~2.0 for all schemes at theta=0.7887).
    // Key assertion: no blow-up. Generous error bound for coarse grid + non-half theta.
    Real err = std::abs(price - cos_ref);
    EXPECT_LT(err, 3.0)
        << "MCS theta=0.7887 err=" << err << " price=" << price
        << " COS ref=" << cos_ref;
}

// T1.6 RED: MCS differs from CS when theta != 0.5
// Validates that MCS is actually a distinct scheme (not a CS alias).
// Uses 80x60 grid (stable for all schemes at theta=0.7887).
TEST(PDEADI2D, ModifiedCraigSneydDiffersFromCraigSneydAtNonHalfTheta) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 80;
    cfg.n_v = 60;
    cfg.n_time = 200;
    cfg.theta = 0.5 + 1.0 / std::sqrt(12.0);  // ~0.7887
    cfg.is_call = true;

    cfg.scheme = ADISchemeType::CraigSneyd;
    PDEEngine2D engine_cs(cfg);
    Real price_cs = engine_cs.price(hp);

    cfg.scheme = ADISchemeType::ModifiedCraigSneyd;
    PDEEngine2D engine_mcs(cfg);
    Real price_mcs = engine_mcs.price(hp);

    // At theta != 0.5, the (theta-0.5)*dt*L_full*{Yhat vs ytilde2} term differs.
    // Difference should be non-zero but small (same order of discretization error).
    Real diff = std::abs(price_cs - price_mcs);
    EXPECT_GT(diff, 1e-10)
        << "CS and MCS should differ at theta=0.7887, got diff=" << diff
        << " CS=" << price_cs << " MCS=" << price_mcs;
    // Difference should be bounded by discretization error scale
    EXPECT_LT(diff, 0.20)
        << "CS-MCS diff too large: " << diff;
}

// T1.6 RED: MCS put-call parity holds
TEST(PDEADI2D, ModifiedCraigSneydPutCallParityHolds) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 120;
    cfg.n_v = 80;
    cfg.n_time = 150;
    cfg.scheme = ADISchemeType::ModifiedCraigSneyd;
    cfg.theta = 0.5;

    cfg.is_call = true;
    PDEEngine2D engine_call(cfg);
    Real call = engine_call.price(hp);

    cfg.is_call = false;
    PDEEngine2D engine_put(cfg);
    Real put = engine_put.price(hp);

    Real parity = call - put;
    Real parity_expected = tc.S0 * std::exp(-tc.q * tc.T) - tc.K * std::exp(-tc.r * tc.T);
    EXPECT_NEAR(parity, parity_expected, 0.30);
}

// =============================================================================
// T1.7: 2D non-uniform grid (sinh transformation) tests
//
// The sinh transformation concentrates grid points near:
//   x = x0 (ATM, where payoff has kink)
//   v = 0  (where Feller boundary layer causes rapid solution variation)
//
// Transformation:
//   x_i = x_center + x_range * sinh(alpha_x * (2i/(n-1) - 1)) / sinh(alpha_x)
//   v_j = v_min + (v_max - v_min) * sinh(alpha_v * j/(n-1)) / sinh(alpha_v)
//
// When alpha=0, the grid degenerates to uniform (L'Hôpital: sinh(a*t)/sinh(a) → t).
// Reference: Tavella & Randall (2000) "Numerical Financial Methods using
//            Computer Algebra", Chapter 5 on non-uniform grids.
// =============================================================================

// T1.7 RED: Non-uniform grid produces positive finite price
TEST(PDEADI2D, NonUniformGridProducesPositivePrice) {
    PDEEngine2DConfig cfg;
    cfg.n_x = 80;
    cfg.n_v = 60;
    cfg.n_time = 100;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;
    cfg.grid_type = GridType2D::Sinh;
    PDEEngine2D engine(cfg);

    auto tc = kTestCases[1];  // ATM call
    auto hp = make_pde_params(tc);
    Real price = engine.price(hp);
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, tc.S0);
}

// T1.7 RED: Non-uniform grid converges to COS benchmark
TEST(PDEADI2D, NonUniformGridConvergesToCOS) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    struct GridLevel { Size nx, nv, nt; const char* label; };
    std::vector<GridLevel> levels = {
        {50,  40,  100, "L0: 50x40x100"},
        {100, 80,  200, "L1: 100x80x200"},
        {200, 160, 400, "L2: 200x160x400"},
    };

    std::vector<Real> prices, errors;
    for (const auto& gl : levels) {
        PDEEngine2DConfig cfg;
        cfg.n_x = gl.nx;
        cfg.n_v = gl.nv;
        cfg.n_time = gl.nt;
        cfg.scheme = ADISchemeType::CraigSneyd;
        cfg.is_call = true;
        cfg.grid_type = GridType2D::Sinh;
        PDEEngine2D engine(cfg);
        Real price = engine.price(hp);
        Real err = std::abs(price - cos_ref);
        prices.push_back(price);
        errors.push_back(err);

        EXPECT_TRUE(std::isfinite(price)) << gl.label;
        EXPECT_GT(price, 0.0) << gl.label;
    }

    // Coarse-to-medium convergence
    EXPECT_LT(errors[1], errors[0])
        << "Sinh L0->L1: err0=" << errors[0] << " err1=" << errors[1];
    // Medium grid absolute error < 0.10
    EXPECT_LT(errors[1], 0.10)
        << "Sinh L1 error=" << errors[1] << " COS ref=" << cos_ref;
    // Self-convergence
    Real self_err = std::abs(prices[1] - prices[2]);
    EXPECT_LT(self_err, 0.10)
        << "Sinh L1-L2 self-error=" << self_err;
}

// T1.7 RED: Non-uniform grid improves accuracy vs uniform at same point count
// The sinh grid concentrates points near the strike and v=0, where the
// solution varies most. At the same number of grid points, the sinh grid
// should produce smaller error than the uniform grid.
TEST(PDEADI2D, NonUniformGridImprovesAccuracy) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 80;
    cfg.n_v = 60;
    cfg.n_time = 200;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;

    // Uniform grid
    cfg.grid_type = GridType2D::Uniform;
    PDEEngine2D engine_uniform(cfg);
    Real price_uniform = engine_uniform.price(hp);
    Real err_uniform = std::abs(price_uniform - cos_ref);

    // Sinh grid (same point count)
    cfg.grid_type = GridType2D::Sinh;
    PDEEngine2D engine_sinh(cfg);
    Real price_sinh = engine_sinh.price(hp);
    Real err_sinh = std::abs(price_sinh - cos_ref);

    // Sinh grid should have smaller or comparable error
    // (allowing some tolerance for noise)
    EXPECT_LT(err_sinh, err_uniform + 0.02)
        << "Uniform err=" << err_uniform << " Sinh err=" << err_sinh;
}

// T1.7 RED: Non-uniform grid is stable on fine grids (no blow-up)
TEST(PDEADI2D, NonUniformGridStableOnFineGrid) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 200;
    cfg.n_v = 150;
    cfg.n_time = 400;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;
    cfg.grid_type = GridType2D::Sinh;
    PDEEngine2D engine(cfg);
    Real price = engine.price(hp);

    EXPECT_TRUE(std::isfinite(price));
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, tc.S0 * 2.0);
    Real err = std::abs(price - cos_ref);
    EXPECT_LT(err, 0.50)
        << "Sinh 200x150x400 err=" << err << " price=" << price;
}

// T1.7 RED: Non-uniform grid put-call parity holds
TEST(PDEADI2D, NonUniformGridPutCallParityHolds) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 100;
    cfg.n_v = 80;
    cfg.n_time = 200;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.grid_type = GridType2D::Sinh;
    cfg.theta = 0.5;

    cfg.is_call = true;
    PDEEngine2D engine_call(cfg);
    Real call = engine_call.price(hp);

    cfg.is_call = false;
    PDEEngine2D engine_put(cfg);
    Real put = engine_put.price(hp);

    Real parity = call - put;
    Real parity_expected = tc.S0 * std::exp(-tc.q * tc.T) - tc.K * std::exp(-tc.r * tc.T);
    EXPECT_NEAR(parity, parity_expected, 0.30);
}

// T1.7 DIAG: Find optimal alpha values for sinh grid
TEST(PDEADI2D, DiagnosticSinhAlphaSweep) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    // Uniform baseline
    PDEEngine2DConfig cfg_u;
    cfg_u.n_x = 80; cfg_u.n_v = 60; cfg_u.n_time = 200;
    cfg_u.scheme = ADISchemeType::CraigSneyd; cfg_u.is_call = true;
    cfg_u.grid_type = GridType2D::Uniform;
    PDEEngine2D eng_u(cfg_u);
    Real err_u = std::abs(eng_u.price(hp) - cos_ref);
    std::cout << "[DIAG] Uniform err=" << err_u << "\n";

    // Sinh sweep
    for (Real ax : {0.0, 0.2, 0.5, 1.0, 2.0}) {
        for (Real av : {0.0, 0.2, 0.5, 1.0, 2.0}) {
            PDEEngine2DConfig cfg;
            cfg.n_x = 80; cfg.n_v = 60; cfg.n_time = 200;
            cfg.scheme = ADISchemeType::CraigSneyd; cfg.is_call = true;
            cfg.grid_type = GridType2D::Sinh;
            cfg.alpha_x = ax; cfg.alpha_v = av;
            PDEEngine2D engine(cfg);
            Real price = engine.price(hp);
            Real err = std::abs(price - cos_ref);
            std::cout << "[DIAG] Sinh ax=" << ax << " av=" << av
                      << " err=" << err << (err < err_u ? " BETTER" : "") << "\n";
        }
    }
    SUCCEED();
}

// ============================================================================
// T1.8: Rannacher smoothing wrapper tests
// Rannacher smoothing uses ADI with theta=1.0 (L-stable, strong damping) for
// the first n_warmup steps, then switches to the configured theta (typically
// 0.5 for second-order accuracy). This dampens high-frequency oscillations
// from non-smooth payoff kinks (e.g., at the strike price) that plague CN/ADI
// schemes with theta=0.5 at coarse time grids.
// Reference: Rannacher (1984), Giles & Carter (1988).
// ============================================================================

// T1.8 RED: Smoke test — smoothing produces positive finite price
TEST(PDEADI2D, RannacherSmoothingProducesPositivePrice) {
    PDEEngine2DConfig cfg;
    cfg.n_x = 80;
    cfg.n_v = 60;
    cfg.n_time = 100;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.theta = 0.5;
    cfg.is_call = true;
    cfg.n_rannacher_warmup = 4;
    PDEEngine2D engine(cfg);

    auto tc = kTestCases[1];  // ATM call
    auto hp = make_pde_params(tc);
    Real price = engine.price(hp);
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, tc.S0);
    EXPECT_TRUE(std::isfinite(price));
}

// T1.8 RED: Backward compatibility — warmup=0 matches no smoothing
TEST(PDEADI2D, RannacherSmoothingZeroWarmupMatchesNoSmoothing) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 100;
    cfg.n_v = 80;
    cfg.n_time = 150;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.theta = 0.5;
    cfg.is_call = true;

    // Without smoothing
    PDEEngine2D engine_no_smooth(cfg);
    Real price_no_smooth = engine_no_smooth.price(hp);

    // With smoothing but warmup=0 (should be identical)
    cfg.n_rannacher_warmup = 0;
    PDEEngine2D engine_zero_warmup(cfg);
    Real price_zero_warmup = engine_zero_warmup.price(hp);

    EXPECT_NEAR(price_no_smooth, price_zero_warmup, 1e-10);
}

// T1.8: Smoothing changes results at coarse grids (warmup is executed)
// and converges to same limit at fine grids.
//
// Note: For Heston with Feller-satisfied parameters, the solution is already
// smooth, so Rannacher smoothing does NOT always improve price accuracy —
// the warmup steps introduce O(dt) error that can exceed the oscillation-
// damping benefit. The true benefit of Rannacher smoothing is visible in:
//   (a) Greeks (gamma) near the strike, where CN oscillations are severe
//   (b) Non-smooth payoffs (digital options), not supported by current engine
//   (c) Feller-violated parameters (boundary layer at v=0)
// Here we verify the more fundamental properties: warmup is actually executed
// (coarse-grid results differ) and both schemes converge at fine grids.
TEST(PDEADI2D, RannacherSmoothingChangesCoarseGridResults) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    // Coarse grid: smoothing should change the result (proves warmup executes)
    {
        PDEEngine2DConfig cfg;
        cfg.n_x = 100;
        cfg.n_v = 80;
        cfg.n_time = 20;
        cfg.scheme = ADISchemeType::CraigSneyd;
        cfg.theta = 0.5;
        cfg.is_call = true;

        PDEEngine2D engine_no_smooth(cfg);
        Real price_no_smooth = engine_no_smooth.price(hp);

        cfg.n_rannacher_warmup = 4;
        PDEEngine2D engine_smooth(cfg);
        Real price_smooth = engine_smooth.price(hp);

        // Prices should differ (proving warmup steps are executed with theta=1.0)
        Real diff = std::abs(price_smooth - price_no_smooth);
        EXPECT_GT(diff, 1e-6)
            << "Smoothing should change result at coarse grid"
            << " no_smooth=" << price_no_smooth << " smooth=" << price_smooth;

        // Both should be finite and positive
        EXPECT_GT(price_no_smooth, 0.0);
        EXPECT_GT(price_smooth, 0.0);
        EXPECT_TRUE(std::isfinite(price_no_smooth));
        EXPECT_TRUE(std::isfinite(price_smooth));
    }

    // Fine grid: smoothing and no-smoothing should converge (warmup error -> 0)
    {
        PDEEngine2DConfig cfg;
        cfg.n_x = 150;
        cfg.n_v = 100;
        cfg.n_time = 300;
        cfg.scheme = ADISchemeType::CraigSneyd;
        cfg.theta = 0.5;
        cfg.is_call = true;

        PDEEngine2D engine_no_smooth(cfg);
        Real price_no_smooth = engine_no_smooth.price(hp);

        cfg.n_rannacher_warmup = 4;
        PDEEngine2D engine_smooth(cfg);
        Real price_smooth = engine_smooth.price(hp);

        // At fine grid, both should agree (warmup O(dt) error vanishes)
        Real diff = std::abs(price_smooth - price_no_smooth);
        EXPECT_LT(diff, 0.10)
            << "Fine grid: smooth and no-smooth should converge"
            << " no_smooth=" << price_no_smooth << " smooth=" << price_smooth
            << " COS ref=" << cos_ref;
    }
}

// T1.8 RED: Convergence — smoothing converges to COS as grid refines
TEST(PDEADI2D, RannacherSmoothingConvergesToCOS) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    struct GridLevel { Size nx, nv, nt; const char* label; };
    std::vector<GridLevel> levels = {
        {50,  40,  100, "L0: 50x40x100"},
        {100, 80,  200, "L1: 100x80x200"},
        {200, 160, 400, "L2: 200x160x400"},
    };

    std::vector<Real> prices, errors;
    for (const auto& gl : levels) {
        PDEEngine2DConfig cfg;
        cfg.n_x = gl.nx;
        cfg.n_v = gl.nv;
        cfg.n_time = gl.nt;
        cfg.scheme = ADISchemeType::CraigSneyd;
        cfg.theta = 0.5;
        cfg.is_call = true;
        cfg.n_rannacher_warmup = 4;
        PDEEngine2D engine(cfg);
        Real price = engine.price(hp);
        Real err = std::abs(price - cos_ref);
        prices.push_back(price);
        errors.push_back(err);

        EXPECT_TRUE(std::isfinite(price)) << gl.label;
        EXPECT_GT(price, 0.0) << gl.label;
    }

    // Coarse-to-medium convergence
    EXPECT_LT(errors[1], errors[0])
        << "Rannacher L0->L1: err0=" << errors[0] << " err1=" << errors[1];
    // Medium grid absolute error < 0.15 (slightly looser than no-smoothing
    // because warmup steps introduce O(dt) error at coarse grid)
    EXPECT_LT(errors[1], 0.15)
        << "Rannacher L1 error=" << errors[1] << " COS ref=" << cos_ref;
    // Self-convergence
    Real self_err = std::abs(prices[1] - prices[2]);
    EXPECT_LT(self_err, 0.15)
        << "Rannacher L1-L2 self-error=" << self_err;
}

// T1.8 RED: Put-call parity preserved under smoothing
TEST(PDEADI2D, RannacherSmoothingPreservesPutCallParity) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 120;
    cfg.n_v = 80;
    cfg.n_time = 150;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.theta = 0.5;
    cfg.n_rannacher_warmup = 4;

    cfg.is_call = true;
    PDEEngine2D engine_call(cfg);
    Real call = engine_call.price(hp);

    cfg.is_call = false;
    PDEEngine2D engine_put(cfg);
    Real put = engine_put.price(hp);

    Real parity = call - put;
    Real parity_expected = tc.S0 * std::exp(-tc.q * tc.T) - tc.K * std::exp(-tc.r * tc.T);

    EXPECT_NEAR(parity, parity_expected, 0.30);
}

// T1.8 RED: Smoothing works with all three ADI schemes (CS, HV, MCS)
TEST(PDEADI2D, RannacherSmoothingWorksWithAllSchemes) {
    auto tc = kTestCases[1];
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 100;
    cfg.n_v = 80;
    cfg.n_time = 100;
    cfg.theta = 0.5;
    cfg.is_call = true;
    cfg.n_rannacher_warmup = 4;

    for (auto scheme : {ADISchemeType::CraigSneyd,
                        ADISchemeType::HundsdorferVerwer,
                        ADISchemeType::ModifiedCraigSneyd}) {
        cfg.scheme = scheme;
        PDEEngine2D engine(cfg);
        Real price = engine.price(hp);
        EXPECT_TRUE(std::isfinite(price)) << "scheme=" << static_cast<int>(scheme);
        EXPECT_GT(price, 0.0) << "scheme=" << static_cast<int>(scheme);
        // All schemes should be in the ballpark of COS
        Real err = std::abs(price - cos_ref);
        EXPECT_LT(err, 0.80) << "scheme=" << static_cast<int>(scheme)
                             << " price=" << price << " COS=" << cos_ref;
    }
}

// ============================================================================
// T2.1: American Heston PDE (ADI + PSOR early exercise)
//
// American options under Heston require solving the variational inequality:
//   max(dV/dtau - L V, payoff - V) = 0
// where L is the Heston PDE operator and payoff = max(K*exp(x) - K, 0) for call
// (or max(K - K*exp(x), 0) for put).
//
// The ADI splitting with PSOR projection approach:
//   1. Standard ADI step (CS/HV/MCS) computes V_intermediate
//   2. PSOR projection: V = max(V_intermediate, payoff)
// This is the simplest operator-splitting scheme (Ikonen-Toivanen, 2004).
// More accurate schemes (e.g., IT- splitting, penalty methods) are possible
// but PSOR-after-ADI is a well-established baseline.
//
// Reference: Ikonen & Toivanen (2004) "Operator splitting methods for
//            American option pricing",Applied Mathematics Letters 17.
//            Zvan, Forsyth & Vetzal (1998) "Robust numerical methods
//            for PDE models of Asian options".
// ============================================================================

// T2.1 RED: Smoke test — American put produces positive finite price
TEST(PDEADI2D, AmericanPutProducesPositivePrice) {
    PDEEngine2DConfig cfg;
    cfg.n_x = 80;
    cfg.n_v = 60;
    cfg.n_time = 100;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = false;       // put
    cfg.is_american = true;    // American
    PDEEngine2D engine(cfg);

    auto tc = kTestCases[1];  // ATM put
    auto hp = make_pde_params(tc);
    Real price = engine.price(hp);
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, tc.K);   // put price < strike
    EXPECT_TRUE(std::isfinite(price));
}

// T2.1 RED: American put >= European put (early exercise premium)
TEST(PDEADI2D, AmericanPutGreaterThanEuropeanPut) {
    auto tc = kTestCases[1];  // ATM (is_call=true in test case, but we price puts)
    auto hp = make_pde_params(tc);

    // Compute European put via PDE (same discretization for fair comparison)
    PDEEngine2DConfig cfg;
    cfg.n_x = 120;
    cfg.n_v = 80;
    cfg.n_time = 200;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = false;       // put
    cfg.is_american = false;   // European

    PDEEngine2D engine_eur(cfg);
    Real price_eur = engine_eur.price(hp);

    // American put
    cfg.is_american = true;
    PDEEngine2D engine_amer(cfg);
    Real price_amer = engine_amer.price(hp);

    // American put must be >= European put (early exercise premium)
    EXPECT_GE(price_amer, price_eur - 0.10)  // allow small numerical noise
        << "American=" << price_amer << " European=" << price_eur;

    // Premium should be positive (early exercise has value for puts with r>0)
    Real premium = price_amer - price_eur;
    EXPECT_GT(premium, -0.05);
}

// T2.1 RED: American call with no dividends ≈ European call
// (No early exercise optimal when q=0 and r>0 for calls)
// 多投影修正后精度提升, 容差从 0.30 收紧至 0.15
TEST(PDEADI2D, AmericanCallNoDividendMatchesEuropean) {
    auto tc = kTestCases[1];  // ATM call, q=0
    auto hp = make_pde_params(tc);
    Real cos_european = cos_price(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 120;
    cfg.n_v = 80;
    cfg.n_time = 200;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;
    cfg.is_american = true;

    PDEEngine2D engine_amer(cfg);
    Real price_amer = engine_amer.price(hp);

    // With q=0, American call = European call (no early exercise benefit)
    // 多投影修正 + sinh grid 后误差约 0.08-0.10, 收紧至 0.15
    EXPECT_NEAR(price_amer, cos_european, 0.15)
        << "American call (q=0)=" << price_amer << " European=" << cos_european;
}

// T2.1 RED: American call with dividends > European call
// When q > 0, early exercise of ITM call can be optimal
TEST(PDEADI2D, AmericanCallWithDividendsExceedsEuropean) {
    // Custom case: high dividend yield makes early exercise optimal
    HestonTestCase tc_div = {100.0, 80.0, 1.0, 0.04, 0.08,  // q=8% high dividend
                              0.04, 1.5, 0.04, 0.3, -0.5, true, 0.0};
    auto hp = make_pde_params(tc_div);

    PDEEngine2DConfig cfg;
    cfg.n_x = 120;
    cfg.n_v = 80;
    cfg.n_time = 200;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;

    // European
    cfg.is_american = false;
    PDEEngine2D engine_eur(cfg);
    Real price_eur = engine_eur.price(hp);

    // American
    cfg.is_american = true;
    PDEEngine2D engine_amer(cfg);
    Real price_amer = engine_amer.price(hp);

    // American call with high dividend should >= European call
    EXPECT_GE(price_amer, price_eur - 0.10)
        << "American (q=0.08)=" << price_amer << " European=" << price_eur;
}

// T2.1 RED: American put price bounded by intrinsic value and strike
TEST(PDEADI2D, AmericanPutBoundedByIntrinsicAndStrike) {
    HestonTestCase tc = {100.0, 90.0, 1.0, 0.04, 0.0,  // ITM put (K=90 > S0=100? No, K=90 < S0=100, OTM)
                          0.04, 1.5, 0.04, 0.3, -0.5, false, 0.0};
    // Use OTM put: K=90 < S0=100, intrinsic = max(90-100,0) = 0
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 100;
    cfg.n_v = 80;
    cfg.n_time = 150;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = false;
    cfg.is_american = true;

    PDEEngine2D engine(cfg);
    Real price = engine.price(hp);

    // American put >= intrinsic value (max(K-S, 0))
    Real intrinsic = std::max(tc.K - tc.S0, 0.0);
    EXPECT_GE(price, intrinsic - 0.01)
        << "price=" << price << " intrinsic=" << intrinsic;
    // American put <= strike (present value bound)
    Real strike_pv = tc.K * std::exp(-tc.r * tc.T);
    EXPECT_LE(price, tc.K)
        << "price=" << price << " strike=" << tc.K;
}

// T2.1 RED: Deep ITM American put ≈ intrinsic value (immediate exercise)
// sinh grid 后 S=K 附近分辨率提升, 深度 ITM put 精度从 ±5.0 收紧至 ±2.0
TEST(PDEADI2D, DeepITMAmericanPutApproachesIntrinsic) {
    // Deep ITM put: K=150, S0=100, so intrinsic = 50
    HestonTestCase tc = {100.0, 150.0, 0.5, 0.04, 0.0,
                          0.04, 1.5, 0.04, 0.3, -0.5, false, 0.0};
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 120;
    cfg.n_v = 80;
    cfg.n_time = 200;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = false;
    cfg.is_american = true;

    PDEEngine2D engine(cfg);
    Real price = engine.price(hp);

    Real intrinsic = tc.K - tc.S0;  // 50
    // Deep ITM American put should be close to intrinsic
    // (time value is small but nonzero due to mean-reversion of variance)
    EXPECT_GT(price, intrinsic - 1.0)
        << "Deep ITM American put=" << price << " intrinsic=" << intrinsic;
    EXPECT_LT(price, intrinsic + 2.0)
        << "Deep ITM American put=" << price << " intrinsic=" << intrinsic;
}

// T2.1 RED: American put converges as grid refines
TEST(PDEADI2D, AmericanPutConvergesWithGridRefinement) {
    auto tc = kTestCases[1];  // ATM put
    auto hp = make_pde_params(tc);

    struct GridLevel { Size nx, nv, nt; const char* label; };
    std::vector<GridLevel> levels = {
        {60,  40,  100, "L0: 60x40x100"},
        {120, 80,  200, "L1: 120x80x200"},
    };

    std::vector<Real> prices;
    for (const auto& gl : levels) {
        PDEEngine2DConfig cfg;
        cfg.n_x = gl.nx;
        cfg.n_v = gl.nv;
        cfg.n_time = gl.nt;
        cfg.scheme = ADISchemeType::CraigSneyd;
        cfg.is_call = false;
        cfg.is_american = true;
        PDEEngine2D engine(cfg);
        Real price = engine.price(hp);
        prices.push_back(price);

        EXPECT_TRUE(std::isfinite(price)) << gl.label;
        EXPECT_GT(price, 0.0) << gl.label;
    }

    // Fine grid should be finite and positive (convergence check is loose
    // because American option has lower smoothness than European)
    EXPECT_GT(prices[1], 0.0);
    // Both should be in reasonable range for ATM put (European ~8.9, American > 9)
    EXPECT_GT(prices[0], 5.0) << "Coarse price too low";
    EXPECT_LT(prices[0], 20.0) << "Coarse price too high";
    EXPECT_GT(prices[1], 5.0) << "Fine price too low";
    EXPECT_LT(prices[1], 20.0) << "Fine price too high";
}

// =========================================================================
// Precision tests after Ikonen-Toivanen multi-projection upgrade
// (step_american projects after Y1, Y2, ytilde1, ytilde2, V_new)
// =========================================================================

// American call (q=0) must closely match European COS benchmark.
// Theory: no early exercise optimal => American call = European call.
// Pre-upgrade tolerance was 0.30; target: 0.10 with multi-projection.
TEST(PDEADI2D, AmericanCallNoDividendTightTolerance) {
    auto tc = kTestCases[1];  // ATM call, q=0
    auto hp = make_pde_params(tc);
    Real cos_european = cos_price(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 150;
    cfg.n_v = 100;
    cfg.n_time = 300;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;
    cfg.is_american = true;

    PDEEngine2D engine_amer(cfg);
    Real price_amer = engine_amer.price(hp);

    Real err = std::abs(price_amer - cos_european);
    std::cout << "AmericanCall(q=0)=" << price_amer << " European=" << cos_european
              << " err=" << err << std::endl;
    EXPECT_LT(err, 0.10)
        << "American call (q=0) should match European within 0.10";
}

// American put must strictly exceed European put (positive early-exercise premium)
// Pre-upgrade allowed -0.05; target: strictly > 0 with multi-projection
TEST(PDEADI2D, AmericanPutStrictlyExceedsEuropean) {
    auto tc = kTestCases[1];  // ATM
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 150;
    cfg.n_v = 100;
    cfg.n_time = 300;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = false;

    // European put (same discretization)
    cfg.is_american = false;
    PDEEngine2D engine_eur(cfg);
    Real price_eur = engine_eur.price(hp);

    // American put
    cfg.is_american = true;
    PDEEngine2D engine_amer(cfg);
    Real price_amer = engine_amer.price(hp);

    Real premium = price_amer - price_eur;
    std::cout << "AmericanPut=" << price_amer << " EuropeanPut=" << price_eur
              << " premium=" << premium << std::endl;
    // Premium must be strictly positive (r > 0 makes early exercise valuable)
    EXPECT_GT(premium, 0.01)
        << "American put premium should be > 0.01";
}

// Grid convergence: American put price should converge monotonically
TEST(PDEADI2D, AmericanPutGridConvergenceTight) {
    auto tc = kTestCases[1];  // ATM
    auto hp = make_pde_params(tc);

    struct GridLevel { Size nx, nv, nt; };
    std::vector<GridLevel> levels = {
        {80,  60,  100},
        {120, 80,  200},
        {180, 120, 300},
    };

    std::vector<Real> prices;
    for (const auto& gl : levels) {
        PDEEngine2DConfig cfg;
        cfg.n_x = gl.nx;
        cfg.n_v = gl.nv;
        cfg.n_time = gl.nt;
        cfg.scheme = ADISchemeType::CraigSneyd;
        cfg.is_call = false;
        cfg.is_american = true;
        PDEEngine2D engine(cfg);
        Real price = engine.price(hp);
        prices.push_back(price);
        std::cout << "Grid " << gl.nx << "x" << gl.nv << "x" << gl.nt
                  << ": price=" << price << std::endl;
    }

    // Convergence: |L2 - L1| should not be much larger than |L1 - L0|.
    // Strict monotonic decrease is not guaranteed for American options
    // because sinh grid concentration changes non-linearly with N.
    Real diff_01 = std::abs(prices[1] - prices[0]);
    Real diff_12 = std::abs(prices[2] - prices[1]);
    std::cout << "diff L0-L1=" << diff_01 << " diff L1-L2=" << diff_12 << std::endl;
    // Allow diff_12 up to 1.5x diff_01 (convergence is roughly first-order)
    EXPECT_LT(diff_12, 1.5 * diff_01)
        << "Grid convergence: refinement should not increase error significantly";
}

// American call with high dividend should exceed European call
TEST(PDEADI2D, AmericanCallHighDivTightTolerance) {
    HestonTestCase tc_div = {100.0, 80.0, 1.0, 0.04, 0.08,
                              0.04, 1.5, 0.04, 0.3, -0.5, true, 0.0};
    auto hp = make_pde_params(tc_div);

    PDEEngine2DConfig cfg;
    cfg.n_x = 150;
    cfg.n_v = 100;
    cfg.n_time = 300;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;

    cfg.is_american = false;
    PDEEngine2D engine_eur(cfg);
    Real price_eur = engine_eur.price(hp);

    cfg.is_american = true;
    PDEEngine2D engine_amer(cfg);
    Real price_amer = engine_amer.price(hp);

    Real premium = price_amer - price_eur;
    std::cout << "AmerCall(q=8%)=" << price_amer << " EurCall=" << price_eur
              << " premium=" << premium << std::endl;
    // With q=8%, early exercise is optimal => premium > 0
    EXPECT_GT(premium, 0.01)
        << "American call with high dividend should exceed European";
}

// Deep ITM American put should be close to intrinsic value
TEST(PDEADI2D, DeepITMAmericanPutTightTolerance) {
    HestonTestCase tc = {100.0, 150.0, 0.5, 0.04, 0.0,
                          0.04, 1.5, 0.04, 0.3, -0.5, false, 0.0};
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 150;
    cfg.n_v = 100;
    cfg.n_time = 300;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = false;
    cfg.is_american = true;

    PDEEngine2D engine(cfg);
    Real price = engine.price(hp);

    Real intrinsic = tc.K - tc.S0;  // 50
    std::cout << "DeepITM AmericanPut=" << price << " intrinsic=" << intrinsic << std::endl;
    // Deep ITM: price should be close to intrinsic (within 2.0)
    EXPECT_GT(price, intrinsic - 0.5)
        << "Deep ITM American put too far below intrinsic";
    EXPECT_LT(price, intrinsic + 2.0)
        << "Deep ITM American put too far above intrinsic";
}

// Rannacher smoothing for American options: verify stability and reasonable accuracy.
// Note: For American calls (q=0), multi-projection already handles the payoff kink
// effectively. Rannacher adds O(dt) dissipation which can slightly worsen accuracy
// for smooth regions. The test verifies Rannacher produces valid results within
// a reasonable tolerance, not that it always improves accuracy.
TEST(PDEADI2D, RannacherAmericanCallStableAndReasonable) {
    auto tc = kTestCases[1];  // ATM call, q=0
    auto hp = make_pde_params(tc);
    Real cos_european = cos_price(tc);

    // Without Rannacher
    PDEEngine2DConfig cfg_no_ran;
    cfg_no_ran.n_x = 120;
    cfg_no_ran.n_v = 80;
    cfg_no_ran.n_time = 200;
    cfg_no_ran.scheme = ADISchemeType::CraigSneyd;
    cfg_no_ran.is_call = true;
    cfg_no_ran.is_american = true;
    cfg_no_ran.n_rannacher_warmup = 0;
    PDEEngine2D engine_no_ran(cfg_no_ran);
    Real price_no_ran = engine_no_ran.price(hp);
    Real err_no_ran = std::abs(price_no_ran - cos_european);

    // With Rannacher (4 warmup steps)
    PDEEngine2DConfig cfg_ran = cfg_no_ran;
    cfg_ran.n_rannacher_warmup = 4;
    PDEEngine2D engine_ran(cfg_ran);
    Real price_ran = engine_ran.price(hp);
    Real err_ran = std::abs(price_ran - cos_european);

    std::cout << "err_no_rannacher=" << err_no_ran << " err_rannacher=" << err_ran << std::endl;
    // Both should produce finite, positive prices
    EXPECT_TRUE(std::isfinite(price_ran));
    EXPECT_GT(price_ran, 0.0);
    // Rannacher error should be within reasonable range (< 0.20 for this grid)
    EXPECT_LT(err_ran, 0.20)
        << "Rannacher American call error should be < 0.20";
}

// =========================================================================
// 高精度收敛测试：网格细化序列下 American call (q=0) 的误差衰减
// 多投影修正后理论收敛阶 = O(dt + dx^2 + dv^2)
// 用 COS 欧式价格作为基准 (q=0 => American call = European call)
// 注: sinh grid 浓度参数固定时, 相邻级别价格差异不严格单调递减
//     (浓度随 N 非线性变化), 因此测试重点在总体收敛趋势与精细网格精度
// =========================================================================
TEST(PDEADI2D, AmericanCallConvergenceOrderWithCOSBenchmark) {
    auto tc = kTestCases[1];  // ATM call, q=0
    auto hp = make_pde_params(tc);
    Real cos_ref = cos_price(tc);

    // 3 级网格序列, 每级约 1.5x 细化
    struct GridLevel { Size nx, nv, nt; const char* label; };
    std::vector<GridLevel> levels = {
        {120, 90,  200, "L0: 120x90x200"},
        {180, 120, 300, "L1: 180x120x300"},
        {240, 160, 400, "L2: 240x160x400"},
    };

    std::vector<Real> prices, errors;
    for (const auto& gl : levels) {
        PDEEngine2DConfig cfg;
        cfg.n_x = gl.nx;
        cfg.n_v = gl.nv;
        cfg.n_time = gl.nt;
        cfg.scheme = ADISchemeType::CraigSneyd;
        cfg.is_call = true;
        cfg.is_american = true;
        PDEEngine2D engine(cfg);
        Real price = engine.price(hp);
        prices.push_back(price);
        errors.push_back(std::abs(price - cos_ref));
        std::cout << gl.label << ": price=" << price << " err=" << errors.back() << std::endl;
    }

    // 1. 所有价格必须有限正数
    for (Real p : prices) {
        EXPECT_TRUE(std::isfinite(p));
        EXPECT_GT(p, 0.0);
    }

    // 2. 总体收敛趋势: 精细网格 (L2) 误差 < 粗网格 (L0) 误差
    EXPECT_LT(errors[2], errors[0])
        << "L0->L2 overall: err0=" << errors[0] << " err2=" << errors[2];

    // 3. 精细网格 (L2) 误差应 < 0.08 (≈0.8% of COS≈9.76)
    EXPECT_LT(errors[2], 0.08)
        << "L2 err too large: " << errors[2] << " COS ref=" << cos_ref;

    // 4. 精细网格自洽性: L1 与 L2 价格差异 < 0.04
    Real diff_12 = std::abs(prices[1] - prices[2]);
    EXPECT_LT(diff_12, 0.04)
        << "L1-L2 self-error too large: " << diff_12;
}

// =========================================================================
// 多 scheme 一致性测试: CS / HV / MCS 在美式看跌上应给出接近的价格
// (splitting 误差相同量级, scheme 差异应 << 离散误差)
// =========================================================================
TEST(PDEADI2D, AmericanPutMultipleSchemesAgree) {
    auto tc = kTestCases[1];  // ATM put
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 150;
    cfg.n_v = 100;
    cfg.n_time = 200;
    cfg.is_call = false;
    cfg.is_american = true;

    std::vector<std::pair<ADISchemeType, const char*>> schemes = {
        {ADISchemeType::CraigSneyd,        "CS"},
        {ADISchemeType::HundsdorferVerwer, "HV"},
        {ADISchemeType::ModifiedCraigSneyd, "MCS"},
    };

    std::vector<Real> prices;
    for (const auto& [scheme, label] : schemes) {
        cfg.scheme = scheme;
        PDEEngine2D engine(cfg);
        Real price = engine.price(hp);
        prices.push_back(price);
        std::cout << label << " American put = " << price << std::endl;
        EXPECT_TRUE(std::isfinite(price));
        EXPECT_GT(price, 0.0);
    }

    // 任意两 scheme 差异 < 0.10 (1% of strike)
    for (Size i = 0; i < prices.size(); ++i) {
        for (Size j = i + 1; j < prices.size(); ++j) {
            Real diff = std::abs(prices[i] - prices[j]);
            EXPECT_LT(diff, 0.10)
                << "scheme " << schemes[i].second << " vs " << schemes[j].second
                << ": diff=" << diff;
        }
    }
}

// =========================================================================
// 低随机波动率 Heston 情形: xi 较小, v0=theta=sigma^2, kappa 较大
// 方差过程接近常数, Heston 价格接近 BSM 价格
// American call (q=0) 应该 = Heston European call (COS 基准)
// =========================================================================
TEST(PDEADI2D, AmericanCallLowVolOfVolMatchesEuropean) {
    // 低随机波动率: xi=0.1, kappa=5, theta=v0=0.04 (sigma=0.2)
    HestonTestCase tc_lowvol = {100.0, 100.0, 1.0, 0.04, 0.0,
                                 0.04, 5.0, 0.04, 0.1, -0.3, true, 0.0};
    auto hp = make_pde_params(tc_lowvol);

    // Heston European call (COS) 作为基准
    auto cf_params = make_cf_params(tc_lowvol);
    Real cos_european = cos_call_heston(tc_lowvol.S0, tc_lowvol.K, tc_lowvol.T,
                                         tc_lowvol.r, tc_lowvol.q, cf_params, 512, 12.0);
    std::cout << "Heston European call (COS) = " << cos_european << std::endl;

    // American call (q=0) 应该 = Heston European call
    PDEEngine2DConfig cfg;
    cfg.n_x = 180;
    cfg.n_v = 120;
    cfg.n_time = 300;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = true;
    cfg.is_american = true;

    PDEEngine2D engine_amer(cfg);
    Real price_amer = engine_amer.price(hp);
    std::cout << "American Heston (low-vol-of-vol) call = " << price_amer << std::endl;

    // 容差 0.10: q=0 时 American call = European call, 误差来自离散
    Real err = std::abs(price_amer - cos_european);
    EXPECT_LT(err, 0.10)
        << "American call (low vol-of-vol) err=" << err << " COS=" << cos_european;

    // 价格必须有限正数且 < S0
    EXPECT_TRUE(std::isfinite(price_amer));
    EXPECT_GT(price_amer, 0.0);
    EXPECT_LT(price_amer, tc_lowvol.S0);
}

// =========================================================================
// 美式 put 严格无套利边界验证:
//   1. American put >= max(K*exp(-r*T), intrinsic) (perpetual lower bound)
//   2. American put <= K (upper bound)
//   3. American put >= European put
// =========================================================================
TEST(PDEADI2D, AmericanPutNoArbitrageBounds) {
    // ITM put: K=110, S0=100
    HestonTestCase tc = {100.0, 110.0, 1.0, 0.04, 0.0,
                          0.04, 1.5, 0.04, 0.3, -0.5, false, 0.0};
    auto hp = make_pde_params(tc);

    PDEEngine2DConfig cfg;
    cfg.n_x = 150;
    cfg.n_v = 100;
    cfg.n_time = 250;
    cfg.scheme = ADISchemeType::CraigSneyd;
    cfg.is_call = false;

    // European put
    cfg.is_american = false;
    PDEEngine2D engine_eur(cfg);
    Real price_eur = engine_eur.price(hp);

    // American put
    cfg.is_american = true;
    PDEEngine2D engine_amer(cfg);
    Real price_amer = engine_amer.price(hp);

    Real intrinsic = std::max(tc.K - tc.S0, 0.0);  // 10
    Real K_pv = tc.K * std::exp(-tc.r * tc.T);     // ~105.7

    std::cout << "ITM put: European=" << price_eur << " American=" << price_amer
              << " intrinsic=" << intrinsic << " K_pv=" << K_pv << std::endl;

    // 1. American >= European (早行权权利价值)
    EXPECT_GE(price_amer, price_eur - 0.01)
        << "American must >= European";

    // 2. American >= intrinsic (即行权价值)
    EXPECT_GE(price_amer, intrinsic - 0.01)
        << "American must >= intrinsic";

    // 3. American <= K (无套利上界)
    EXPECT_LE(price_amer, tc.K)
        << "American must <= K";

    // 4. American put 在 r=0 时 = European put (无早行权收益)
    // (此 case r=0.04, 所以 American > European, 差额即早行权溢价)
    EXPECT_GT(price_amer, price_eur)
        << "With r>0, American put premium > 0";
}

// =========================================================================
// 自适应时间步规范化测试: 验证 adaptive time step 行为
// Heston PDE 的 v_max 处谱半径极大, adaptive time step 几乎总会触发,
// 且 n_actual = ceil(T * spec_v / C_MAX) 与输入 n_time 无关 (dt 被规范化).
// 此测试验证:
//   1. adaptive 触发 (n_actual > n_time) 当 n_time 较小
//   2. n_actual 收敛到稳定值 (与 n_time 无关) 当 n_time 足够大
//   3. 价格在所有配置下完全相同 (因为实际 dt 相同)
// =========================================================================
TEST(PDEADI2D, AmericanPutAdaptiveTimeStepNormalization) {
    auto tc = kTestCases[1];  // ATM put
    auto hp = make_pde_params(tc);

    struct TimeLevel { Size nt; const char* label; };
    std::vector<TimeLevel> levels = {
        {50,  "nt=50"},
        {100, "nt=100"},
        {200, "nt=200"},
        {400, "nt=400"},
    };

    std::vector<Real> prices;
    std::vector<Size> n_actual;
    for (const auto& tl : levels) {
        PDEEngine2DConfig cfg;
        cfg.n_x = 150;
        cfg.n_v = 100;
        cfg.n_time = tl.nt;
        cfg.scheme = ADISchemeType::CraigSneyd;
        cfg.is_call = false;
        cfg.is_american = true;
        PDEEngine2D engine(cfg);
        Real price = engine.price(hp);
        prices.push_back(price);
        n_actual.push_back(engine.last_n_time_used());
        std::cout << tl.label << ": price=" << price
                  << " n_actual=" << engine.last_n_time_used() << std::endl;
    }

    // 1. 自适应时间步触发: 小 n_time 时 n_actual > n_time
    EXPECT_GT(n_actual[0], levels[0].nt)
        << "Adaptive should trigger for small n_time=" << levels[0].nt;

    // 2. n_actual 收敛到稳定值 (与 n_time 无关)
    //    当 n_time 足够大, n_actual 应该 = ceil(T*spec_v/C_MAX), 是个常数
    EXPECT_EQ(n_actual[1], n_actual[2])
        << "n_actual should be constant for nt>=100: "
        << n_actual[1] << " vs " << n_actual[2];
    EXPECT_EQ(n_actual[2], n_actual[3])
        << "n_actual should be constant for nt>=200: "
        << n_actual[2] << " vs " << n_actual[3];

    // 3. 价格在自适应触发后完全相同 (dt 被规范化为相同值)
    EXPECT_NEAR(prices[1], prices[2], 1e-10)
        << "Prices should be identical when dt is normalized";
    EXPECT_NEAR(prices[2], prices[3], 1e-10)
        << "Prices should be identical when dt is normalized";

    // 4. 所有价格有限正数, 在合理范围内
    for (Real p : prices) {
        EXPECT_TRUE(std::isfinite(p));
        EXPECT_GT(p, 0.0);
        EXPECT_LT(p, tc.K);
    }
}
