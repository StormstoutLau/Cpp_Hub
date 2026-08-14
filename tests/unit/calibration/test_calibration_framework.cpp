#include <gtest/gtest.h>
#include "cpphub/calibration/optimizer.hpp"
#include "cpphub/calibration/objective.hpp"
#include "cpphub/calibration/calibrator.hpp"
#include "cpphub/models/vol_surface/svi.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"  // bsm_call_price
#include <cmath>
#include <vector>

using namespace cpphub::v1;

TEST(M3CompileCheck, LMQuadraticResidual) {
    // Quadratic residuals: r = [x0-2, x1-3, x0+x1-5], min at (2,3)
    ResidualFn r = [](const std::vector<Real>& x) -> std::vector<Real> {
        return {x[0] - 2.0, x[1] - 3.0, x[0] + x[1] - 5.0};
    };
    LevenbergMarquardt::Config cfg;
    cfg.max_iterations = 200;
    auto result = LevenbergMarquardt::minimize(r, {0.0, 0.0}, cfg);
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(result.x[0], 2.0, 1e-6);
    EXPECT_NEAR(result.x[1], 3.0, 1e-6);
}

TEST(M3CompileCheck, NelderMeadQuadratic) {
    ObjectiveFn f = [](const std::vector<Real>& x) -> Real {
        return (x[0] - 3.0)*(x[0] - 3.0) + (x[1] + 2.0)*(x[1] + 2.0);
    };
    NelderMead::Config cfg;
    cfg.ftol = 1e-12;
    cfg.xtol = 1e-12;
    auto result = NelderMead::minimize(f, {0.0, 0.0}, cfg);
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(result.x[0], 3.0, 1e-6);
    EXPECT_NEAR(result.x[1], -2.0, 1e-6);
}

TEST(M3CompileCheck, DEGlobalSearch) {
    ObjectiveFn f = [](const std::vector<Real>& x) -> Real {
        return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    };
    std::vector<Bounds> bounds = {{-5, 5}, {-5, 5}, {-5, 5}};
    DifferentialEvolution::Config cfg;
    cfg.max_generations = 100;
    cfg.tol = 1e-10;
    auto result = DifferentialEvolution::minimize(f, bounds, cfg);
    EXPECT_LT(std::abs(result.x[0]), 1e-3);
    EXPECT_LT(std::abs(result.x[1]), 1e-3);
    EXPECT_LT(std::abs(result.x[2]), 1e-3);
}

TEST(M3CompileCheck, SVITotalVariance) {
    SVIParams p{0.04, 0.4, 0.0, 0.1, 0.0};
    SVI s(p);
    EXPECT_NEAR(s.total_variance(0.0), 0.08, 1e-12);
    EXPECT_NEAR(s.implied_vol(0.0, 1.0), std::sqrt(0.08), 1e-12);
}

TEST(M3CompileCheck, SVINoArbitrage) {
    SVIParams p{0.04, 0.4, 0.0, 0.1, 0.0};
    SVI s(p);
    EXPECT_TRUE(s.check_butterfly_arbitrage());
    EXPECT_TRUE(s.check_butterfly_at(0.0));
}

TEST(M3CompileCheck, SVICalibration) {
    SVIParams true_p{0.04, 0.4, -0.3, 0.1, 0.0};
    SVI true_svi(true_p);
    Real T = 1.0;
    Real F = 100.0;
    std::vector<Real> strikes, maturities, ivs;
    for (int i = 0; i < 11; ++i) {
        Real k = -0.4 + 0.08 * i;
        Real K = F * std::exp(k);
        Real w = true_svi.total_variance(k);
        Real iv = std::sqrt(w / T);
        strikes.push_back(K);
        maturities.push_back(T);
        ivs.push_back(iv);
    }
    SVI calib_svi(SVIParams{0.0, 0.1, 0.0, 0.1, 0.0});
    CalibConfig cfg;
    // RISK-013 修复: DE 100 代 × 50 种群对 5 维 SVI 不足,需增加到 500 代 × 100 种群
    // 才能可靠找到全局最优 (DE 收敛后 LM 1 次迭代即满足 gtol)
    cfg.de_pop_size = 100;
    cfg.de_generations = 500;
    cfg.lm_max_iter = 2000;
    cfg.ftol = 1e-14;
    cfg.xtol = 1e-14;
    auto result = calib_svi.calibrate(strikes, maturities, ivs, F, cfg);
    EXPECT_TRUE(result.converged);
    // Check fitted total variance matches at sample points (functional fit)
    SVI fitted(SVIParams{result.params[0], result.params[1], result.params[2],
                         result.params[3], result.params[4]});
    for (int i = 0; i < 11; ++i) {
        Real k = -0.4 + 0.08 * i;
        Real w_fit = fitted.total_variance(k);
        Real w_true = true_svi.total_variance(k);
        EXPECT_NEAR(w_fit, w_true, 1e-4) << "k=" << k;
    }
}

// ===========================================================================
// ObjectiveFunction tests
// ===========================================================================

TEST(ObjectiveFunction, PriceWeightedResiduals) {
    // Two quotes, equal model=market → zero residuals
    std::vector<MarketQuote> quotes = {
        {100.0, 1.0, 10.0, 0.2, 50.0},
        {110.0, 1.0, 5.0, 0.18, 30.0}
    };
    PerQuoteModelFn model = [](const std::vector<Real>&, const MarketQuote& q) -> Real {
        return q.market_price;  // identical → zero residual
    };
    ObjectiveFunction obj(model, quotes, WeightingScheme::PriceWeighted);
    auto r = obj.residuals({});
    EXPECT_EQ(r.size(), 2u);
    EXPECT_NEAR(r[0], 0.0, 1e-12);
    EXPECT_NEAR(r[1], 0.0, 1e-12);
    EXPECT_NEAR(obj.evaluate({}), 0.0, 1e-12);
}

TEST(ObjectiveFunction, RelativeErrorWeighting) {
    // Model = market + 1.0, weights = 1/market, residual_i = w_i * diff_i = 1.0 / market_i
    std::vector<MarketQuote> quotes = {
        {100.0, 1.0, 10.0, 0.2, 0.0},
        {110.0, 1.0, 5.0, 0.18, 0.0}
    };
    PerQuoteModelFn model = [](const std::vector<Real>&, const MarketQuote& q) -> Real {
        return q.market_price + 1.0;
    };
    ObjectiveFunction obj(model, quotes, WeightingScheme::RelativeError);
    auto r = obj.residuals({});
    EXPECT_NEAR(r[0], 1.0 / 10.0, 1e-12);
    EXPECT_NEAR(r[1], 1.0 / 5.0, 1e-12);
    // evaluate = 0.5 * sum( w_i * diff_i^2 ) = 0.5 * (0.1*1 + 0.2*1) = 0.15
    EXPECT_NEAR(obj.evaluate({}), 0.5 * (0.1 * 1.0 + 0.2 * 1.0), 1e-12);
}

TEST(ObjectiveFunction, IVObjectiveComparesIV) {
    // make_iv_objective must compare model output against implied_vol, not market_price
    std::vector<MarketQuote> quotes = {
        {100.0, 1.0, 10.0, 0.20, 0.0},
        {110.0, 1.0, 5.0, 0.18, 0.0}
    };
    // Model returns exactly the market IV → zero residual
    auto iv_fn = [](const std::vector<Real>&, Real K, Real T) -> Real {
        (void)K; (void)T;
        return 0.19;  // midpoint
    };
    auto obj = ObjectiveFunction::make_iv_objective(iv_fn, quotes,
                WeightingScheme::RelativeError);
    auto r = obj.residuals({});
    EXPECT_NEAR(r[0], (0.19 - 0.20) / 0.20, 1e-12);
    EXPECT_NEAR(r[1], (0.19 - 0.18) / 0.18, 1e-12);
}

TEST(ObjectiveFunction, ToObjectiveFnMinimization) {
    // 1D: model = x, market = 5, PriceWeighted: weight = 1/sqrt(|market|)
    std::vector<MarketQuote> quotes = {{100.0, 1.0, 5.0, 0.2, 0.0}};
    PerQuoteModelFn model = [](const std::vector<Real>& x, const MarketQuote&) -> Real {
        return x[0];
    };
    ObjectiveFunction obj(model, quotes, WeightingScheme::PriceWeighted);
    // objective(x) = 0.5 * (1/sqrt(5)) * (x-5)^2
    auto f = obj.to_objective_fn();
    Real w = 1.0 / std::sqrt(5.0);
    Real expected = 0.5 * w * (3.0 - 5.0) * (3.0 - 5.0);
    EXPECT_NEAR(f({3.0}), expected, 1e-12);
    // Minimizer: x = 5
    Real fmin = f({5.0});
    EXPECT_NEAR(fmin, 0.0, 1e-12);
}

TEST(ObjectiveFunction, EmptyQuotesThrows) {
    EXPECT_THROW(
        ObjectiveFunction obj(
            [](const std::vector<Real>&, const MarketQuote&) -> Real { return 0.0; },
            {}, WeightingScheme::PriceWeighted),
        std::invalid_argument);
}

// ===========================================================================
// bsm_implied_vol tests
// ===========================================================================

TEST(BsmImpliedVol, RoundTrip) {
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real C = bsm_call_price(S, K, T, r, q, sigma);
    Real iv = bsm_implied_vol(C, S, K, T, r, q, true);
    EXPECT_NEAR(iv, sigma, 1e-8);

    // OTM call
    K = 110.0;
    C = bsm_call_price(S, K, T, r, q, sigma);
    iv = bsm_implied_vol(C, S, K, T, r, q, true);
    EXPECT_NEAR(iv, sigma, 1e-8);

    // ITM call (higher price, but Newton-Raphson should still converge)
    K = 90.0;
    C = bsm_call_price(S, K, T, r, q, sigma);
    iv = bsm_implied_vol(C, S, K, T, r, q, true);
    EXPECT_NEAR(iv, sigma, 1e-8);
}

TEST(BsmImpliedVol, PutRoundTrip) {
    Real S = 100.0, K = 110.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.25;
    Real P = bsm_put_price(S, K, T, r, q, sigma);
    Real iv = bsm_implied_vol(P, S, K, T, r, q, false);
    EXPECT_NEAR(iv, sigma, 1e-8);
}

// ===========================================================================
// HestonCalibrator tests
// ===========================================================================

TEST(HestonCalibrator, FellerConditionCheck) {
    HestonParams ok{0.04, 2.0, 0.04, 0.2, -0.5};
    EXPECT_TRUE(HestonCalibrator::check_feller(ok));   // 2*2*0.04=0.16 > 0.04
    HestonParams bad{0.04, 0.5, 0.04, 0.3, -0.5};
    EXPECT_FALSE(HestonCalibrator::check_feller(bad)); // 2*0.5*0.04=0.04 < 0.09
}

TEST(HestonCalibrator, DefaultBounds) {
    auto b = HestonCalibrator::default_bounds();
    ASSERT_EQ(b.size(), 5u);
    EXPECT_GT(b[0].upper, b[0].lower);  // v0
    EXPECT_GT(b[1].upper, b[1].lower);  // kappa
    EXPECT_GT(b[2].upper, b[2].lower);  // theta
    EXPECT_GT(b[3].upper, b[3].lower);  // sigma_v
    EXPECT_GT(b[4].upper, b[4].lower);  // rho
    EXPECT_LE(b[4].lower, -0.9);
    EXPECT_GE(b[4].upper, 0.9);
}

TEST(HestonCalibrator, EmptyQuotesFails) {
    HestonCalibrator cal;
    cal.set_market(100.0, 0.05, 0.0);
    CalibConfig cfg;
    cfg.use_de_init = false;
    auto result = cal.calibrate({}, cfg);
    EXPECT_FALSE(result.converged);
}

TEST(HestonCalibrator, RoundTripOnSyntheticQuotes) {
    // Generate call prices from a known Heston model, then invert to IVs and
    // calibrate. We check the recovered parameters are in the right ballpark.
    Real S = 100.0, r = 0.05, q = 0.0;
    HestonParams true_p{0.04, 2.0, 0.04, 0.3, -0.5};
    // Sanity: Feller satisfied (2*2*0.04=0.16 > 0.09)
    ASSERT_TRUE(HestonCalibrator::check_feller(true_p));

    std::vector<MarketQuote> quotes;
    std::vector<Real> strikes = {80.0, 90.0, 100.0, 110.0, 120.0};
    Real T = 1.0;
    for (Real K : strikes) {
        Real price = detail::heston_call_price_cf(S, K, T, r, q, true_p);
        Real iv = bsm_implied_vol(price, S, K, T, r, q, true);
        quotes.push_back({K, T, price, iv, 0.0});
    }

    HestonCalibrator cal;
    cal.set_market(S, r, q);
    CalibConfig cfg;
    cfg.use_de_init = false;  // skip DE for speed in unit test
    cfg.lm_max_iter = 200;
    auto result = cal.calibrate(quotes, cfg);
    // Convergence is not guaranteed from a single LM start; the test passes if
    // either we converged, or at least the residuals are small (functional fit).
    if (result.converged) {
        HestonParams recovered = cal.extract_params(result.params);
        // Check ballpark: v0 should be near 0.04, theta near 0.04
        EXPECT_GT(recovered.v0, 0.0);
        EXPECT_GT(recovered.theta, 0.0);
        EXPECT_GT(recovered.kappa, 0.0);
        EXPECT_GT(recovered.sigma_v, 0.0);
        EXPECT_LT(std::abs(recovered.rho), 1.0);
    }
    // Regardless of convergence, objective should be small (functional match)
    Real max_resid = 0.0;
    for (Real r_i : result.residuals) {
        max_resid = std::max(max_resid, std::abs(r_i));
    }
    // IV residuals use RelativeError weighting (residual = (iv_model-iv_mkt)/iv_mkt)
    EXPECT_LT(max_resid, 0.05) << "max |residual|=" << max_resid;
}

// ===========================================================================
// SABRCalibrator tests
// ===========================================================================

TEST(SABRCalibrator, DefaultBounds) {
    auto b = SABRCalibrator::default_bounds();
    ASSERT_EQ(b.size(), 4u);
    EXPECT_GT(b[0].upper, b[0].lower);  // alpha
    EXPECT_GE(b[1].lower, 0.0);          // beta >= 0
    EXPECT_LE(b[1].upper, 1.0);          // beta <= 1
    EXPECT_GT(b[2].upper, b[2].lower);  // nu
    EXPECT_LE(b[3].lower, -0.9);         // rho
    EXPECT_GE(b[3].upper, 0.9);
}

TEST(SABRCalibrator, EmptyQuotesFails) {
    SABRCalibrator cal;
    cal.set_market(100.0, 0.05, 0.0);
    CalibConfig cfg;
    cfg.use_de_init = false;
    auto result = cal.calibrate({}, cfg);
    EXPECT_FALSE(result.converged);
}

TEST(SABRCalibrator, HaganFormulaATM) {
    // ATM Hagan formula has a closed form; check it returns a sensible value
    Real F = 100.0, K = 100.0, T = 1.0;
    SABRParams sp{0.2, 0.5, 0.3, -0.2};
    Real iv = detail::sabr_implied_vol_hagan(F, K, T, sp);
    EXPECT_GT(iv, 0.0);
    EXPECT_LT(iv, 1.0);
    // ATM vol should be roughly alpha / F^(1-beta) for short maturity
    Real iv_atm_approx = sp.alpha / std::pow(F, 1.0 - sp.beta);
    EXPECT_NEAR(iv, iv_atm_approx, 0.05);  // + higher-order T correction
}

TEST(SABRCalibrator, HaganFormulaSmileShape) {
    // IV should be decreasing in K for negative rho (skew)
    Real F = 100.0, T = 1.0;
    SABRParams sp{0.2, 0.5, 0.3, -0.3};
    Real iv_OTM = detail::sabr_implied_vol_hagan(F, 110.0, T, sp);
    Real iv_ATM = detail::sabr_implied_vol_hagan(F, 100.0, T, sp);
    Real iv_ITM = detail::sabr_implied_vol_hagan(F, 90.0, T, sp);
    // For negative rho: OTM call (K>F) has lower IV than ITM call (K<F)
    EXPECT_LT(iv_OTM, iv_ATM);
    EXPECT_GT(iv_ITM, iv_ATM);
}

TEST(SABRCalibrator, RoundTripOnSyntheticQuotes) {
    // Generate IVs from known SABR params, calibrate, check functional match.
    Real F = 100.0, r = 0.0, q = 0.0;
    SABRParams true_p{0.25, 0.5, 0.4, -0.3};
    std::vector<Real> strikes = {80.0, 90.0, 100.0, 110.0, 120.0};
    Real T = 1.0;
    std::vector<MarketQuote> quotes;
    for (Real K : strikes) {
        Real iv = detail::sabr_implied_vol_hagan(F * std::exp(r * T), K, T, true_p);
        // market_price field is irrelevant for IV calibration, fill with 0
        quotes.push_back({K, T, 0.0, iv, 0.0});
    }

    SABRCalibrator cal;
    cal.set_market(F, r, q);
    CalibConfig cfg;
    cfg.use_de_init = false;
    cfg.lm_max_iter = 500;
    auto result = cal.calibrate(quotes, cfg);
    // Check functional fit: residuals (relative error) should be small
    Real max_resid = 0.0;
    for (Real r_i : result.residuals) {
        max_resid = std::max(max_resid, std::abs(r_i));
    }
    EXPECT_LT(max_resid, 0.05) << "max |residual|=" << max_resid;
    // Recovered params should be in valid region
    if (result.params.size() == 4) {
        SABRParams rec = cal.extract_params(result.params);
        EXPECT_GT(rec.alpha, 0.0);
        EXPECT_GE(rec.beta, 0.0);
        EXPECT_LE(rec.beta, 1.0);
        EXPECT_GT(rec.nu, 0.0);
        EXPECT_LT(std::abs(rec.rho), 1.0);
    }
}

// ===========================================================================
// SABRCalibrator fixed-beta mode tests
// ===========================================================================

TEST(SABRCalibrator, FixedBetaState) {
    SABRCalibrator cal;
    EXPECT_FALSE(cal.has_fixed_beta());
    cal.set_fixed_beta(0.5);
    EXPECT_TRUE(cal.has_fixed_beta());
    EXPECT_NEAR(cal.fixed_beta(), 0.5, 1e-12);
    cal.clear_fixed_beta();
    EXPECT_FALSE(cal.has_fixed_beta());
}

TEST(SABRCalibrator, FixedBetaInvalidRejects) {
    SABRCalibrator cal;
    EXPECT_THROW(cal.set_fixed_beta(-0.1), std::invalid_argument);
    EXPECT_THROW(cal.set_fixed_beta(1.5), std::invalid_argument);
}

TEST(SABRCalibrator, FixedBetaBoundsReduced) {
    auto b_full = SABRCalibrator::default_bounds();
    auto b_fix  = SABRCalibrator::default_bounds_fixed_beta();
    ASSERT_EQ(b_full.size(), 4u);
    ASSERT_EQ(b_fix.size(), 3u);
    // Reduced bounds should drop the beta slot (index 1)
    EXPECT_NEAR(b_fix[0].lower, b_full[0].lower, 1e-12);  // alpha
    EXPECT_NEAR(b_fix[0].upper, b_full[0].upper, 1e-12);
    EXPECT_NEAR(b_fix[1].lower, b_full[2].lower, 1e-12);  // nu (skip beta)
    EXPECT_NEAR(b_fix[1].upper, b_full[2].upper, 1e-12);
    EXPECT_NEAR(b_fix[2].lower, b_full[3].lower, 1e-12);  // rho
    EXPECT_NEAR(b_fix[2].upper, b_full[3].upper, 1e-12);
}

TEST(SABRCalibrator, FixedBetaRoundTrip) {
    // Generate IVs from known SABR params with beta=0.5, fix beta=0.5, calibrate (alpha,nu,rho).
    Real F = 100.0, r = 0.0, q = 0.0;
    SABRParams true_p{0.25, 0.5, 0.4, -0.3};
    std::vector<Real> strikes = {80.0, 90.0, 100.0, 110.0, 120.0};
    Real T = 1.0;
    std::vector<MarketQuote> quotes;
    for (Real K : strikes) {
        Real iv = detail::sabr_implied_vol_hagan(F * std::exp(r * T), K, T, true_p);
        quotes.push_back({K, T, 0.0, iv, 0.0});
    }

    SABRCalibrator cal;
    cal.set_market(F, r, q);
    cal.set_fixed_beta(0.5);  // fix beta to true value
    CalibConfig cfg;
    cfg.use_de_init = false;
    cfg.lm_max_iter = 500;
    auto result = cal.calibrate(quotes, cfg);

    EXPECT_TRUE(result.converged) << result.message;
    ASSERT_EQ(result.params.size(), 3u);  // 3-param form
    SABRParams rec = cal.extract_params(result.params);
    // With beta fixed to the true value, alpha/nu/rho should recover tightly
    EXPECT_NEAR(rec.alpha, true_p.alpha, 0.05);
    EXPECT_NEAR(rec.beta,   0.5,        1e-12);  // exactly fixed
    EXPECT_NEAR(rec.nu,    true_p.nu,   0.05);
    EXPECT_NEAR(rec.rho,   true_p.rho,  0.10);
    // Functional fit should be tight
    Real max_resid = 0.0;
    for (Real r_i : result.residuals) max_resid = std::max(max_resid, std::abs(r_i));
    EXPECT_LT(max_resid, 0.01) << "max |residual|=" << max_resid;
}

TEST(SABRCalibrator, FixedBetaEquityConvention) {
    // Equity convention: beta=0.5. Synthetic smile should calibrate cleanly.
    Real F = 100.0, r = 0.0, q = 0.0;
    SABRParams true_p{0.20, 0.5, 0.3, -0.2};
    std::vector<Real> strikes = {85.0, 95.0, 100.0, 105.0, 115.0};
    Real T = 0.5;
    std::vector<MarketQuote> quotes;
    for (Real K : strikes) {
        Real iv = detail::sabr_implied_vol_hagan(F * std::exp(r * T), K, T, true_p);
        quotes.push_back({K, T, 0.0, iv, 0.0});
    }

    SABRCalibrator cal;
    cal.set_market(F, r, q);
    cal.set_fixed_beta(0.5);  // equity convention
    CalibConfig cfg;
    cfg.use_de_init = false;
    cfg.lm_max_iter = 500;
    auto result = cal.calibrate(quotes, cfg);

    EXPECT_TRUE(result.converged);
    SABRParams rec = cal.extract_params(result.params);
    EXPECT_NEAR(rec.beta, 0.5, 1e-12);
    EXPECT_GT(rec.alpha, 0.0);
    EXPECT_GT(rec.nu, 0.0);
    EXPECT_LT(std::abs(rec.rho), 1.0);
}

TEST(SABRCalibrator, FixedBetaFXNormalVol) {
    // FX normal-vol convention: beta=0.0. Generate smile and fix beta=0.0.
    Real F = 1.10, r = 0.0, q = 0.0;  // EUR/USD-like forward
    SABRParams true_p{0.0080, 0.0, 0.20, -0.1};
    std::vector<Real> strikes = {1.05, 1.08, 1.10, 1.12, 1.15};
    Real T = 1.0;
    std::vector<MarketQuote> quotes;
    for (Real K : strikes) {
        Real iv = detail::sabr_implied_vol_hagan(F * std::exp(r * T), K, T, true_p);
        quotes.push_back({K, T, 0.0, iv, 0.0});
    }

    SABRCalibrator cal;
    cal.set_market(F, r, q);
    cal.set_fixed_beta(0.0);  // FX normal vol
    CalibConfig cfg;
    cfg.use_de_init = false;
    cfg.lm_max_iter = 500;
    auto result = cal.calibrate(quotes, cfg);

    EXPECT_TRUE(result.converged);
    SABRParams rec = cal.extract_params(result.params);
    EXPECT_NEAR(rec.beta, 0.0, 1e-12);
    EXPECT_GT(rec.alpha, 0.0);
}

// ===========================================================================
// v1.1 Task 5: 正则化与早停测试
// ===========================================================================

// --- LM 正则化基础: r=[x-2], 真值 x=2, prior=0, lambda=10 → 结果向 0 偏移 ---
TEST(LMRegularization, PullsTowardPrior) {
    // 残差: r(x) = x - 2, 真值 x*=2
    // 无正则化: min 0.5*(x-2)^2 → x=2
    // 有正则化 (lambda=10, prior=0): min 0.5*(x-2)^2 + 0.5*10*x^2
    //   d/dx = (x-2) + 10*x = 0 → 11x = 2 → x = 2/11 ≈ 0.1818
    ResidualFn r = [](const std::vector<Real>& x) -> std::vector<Real> {
        return {x[0] - 2.0};
    };
    LevenbergMarquardt::Config cfg;
    cfg.max_iterations = 200;
    cfg.lambda_reg = 10.0;
    cfg.params_prior = {0.0};
    auto result = LevenbergMarquardt::minimize(r, {1.0}, cfg);
    EXPECT_TRUE(result.converged);
    // 解析解: x = 2/(1+10) = 2/11
    EXPECT_NEAR(result.x[0], 2.0 / 11.0, 1e-6);
}

// --- LM 正则化 lambda=0 等价于无正则化 ---
TEST(LMRegularization, LambdaZeroEquivalentToNoReg) {
    ResidualFn r = [](const std::vector<Real>& x) -> std::vector<Real> {
        return {x[0] - 2.0, x[1] - 3.0};
    };
    LevenbergMarquardt::Config cfg_no_reg;
    cfg_no_reg.max_iterations = 200;
    auto result_no_reg = LevenbergMarquardt::minimize(r, {0.0, 0.0}, cfg_no_reg);

    LevenbergMarquardt::Config cfg_reg;
    cfg_reg.max_iterations = 200;
    cfg_reg.lambda_reg = 0.0;  // 禁用正则化
    cfg_reg.params_prior = {5.0, 5.0};  // 即使设了 prior, lambda=0 也不影响
    auto result_reg = LevenbergMarquardt::minimize(r, {0.0, 0.0}, cfg_reg);

    EXPECT_NEAR(result_reg.x[0], result_no_reg.x[0], 1e-10);
    EXPECT_NEAR(result_reg.x[1], result_no_reg.x[1], 1e-10);
}

// --- LM 正则化: prior.size() != x.size() 时自动禁用 ---
TEST(LMRegularization, SizeMismatchDisablesReg) {
    ResidualFn r = [](const std::vector<Real>& x) -> std::vector<Real> {
        return {x[0] - 2.0};
    };
    LevenbergMarquardt::Config cfg;
    cfg.max_iterations = 200;
    cfg.lambda_reg = 10.0;
    cfg.params_prior = {0.0, 0.0};  // size=2 != x.size()=1, 应禁用
    auto result = LevenbergMarquardt::minimize(r, {1.0}, cfg);
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(result.x[0], 2.0, 1e-6);  // 无正则化, 回到真值
}

// --- LM 早停: 设置较大 RMSE 阈值, 验证提前停止 ---
TEST(LMEarlyStop, StopsWhenRMSEBelowThreshold) {
    // r(x) = x - 2, 真值 x=2, RMSE = |x-2|
    // 设置 early_stop_rmse = 0.5, 当 |x-2| < 0.5 时停止
    ResidualFn r = [](const std::vector<Real>& x) -> std::vector<Real> {
        return {x[0] - 2.0};
    };

    // 启用早停: RMSE < 0.5 时停止 (宽松阈值, 应在第 1 次迭代就停止)
    LevenbergMarquardt::Config cfg_es;
    cfg_es.max_iterations = 200;
    cfg_es.early_stop_rmse = 0.5;
    auto result_es = LevenbergMarquardt::minimize(r, {0.0}, cfg_es);

    EXPECT_TRUE(result_es.converged);
    // 早停应触发, message 包含 "early_stop"
    EXPECT_NE(result_es.message.find("early_stop"), std::string::npos)
        << "expected early_stop in message, got: " << result_es.message;
    EXPECT_LT(result_es.n_iterations, cfg_es.max_iterations);
}

// --- DE 正则化基础: f(x) = x^2, prior=3, lambda=10 → 结果向 3 偏移 ---
TEST(DERegularization, PullsTowardPrior) {
    // f(x) = x[0]^2 + x[1]^2, 真值 (0,0)
    // 有正则化 (lambda=10, prior=(3,3)):
    //   min x^2+y^2 + 0.5*10*((x-3)^2+(y-3)^2)
    //   d/dx = 2x + 10(x-3) = 0 → 12x = 30 → x = 2.5
    ObjectiveFn f = [](const std::vector<Real>& x) -> Real {
        return x[0] * x[0] + x[1] * x[1];
    };
    std::vector<Bounds> bounds = {{-5, 5}, {-5, 5}};
    DifferentialEvolution::Config cfg;
    cfg.max_generations = 100;
    cfg.tol = 1e-10;
    cfg.lambda_reg = 10.0;
    cfg.params_prior = {3.0, 3.0};
    auto result = DifferentialEvolution::minimize(f, bounds, cfg);
    // 解析解: x = 30/12 = 2.5
    EXPECT_NEAR(result.x[0], 2.5, 0.05);
    EXPECT_NEAR(result.x[1], 2.5, 0.05);
    // 返回的 fx 应为原始目标 (不含正则化项)
    EXPECT_NEAR(result.fx, 2.5 * 2.5 + 2.5 * 2.5, 0.5);
}

// --- DE 正则化 lambda=0 等价于无正则化 ---
TEST(DERegularization, LambdaZeroEquivalentToNoReg) {
    ObjectiveFn f = [](const std::vector<Real>& x) -> Real {
        return (x[0] - 3.0) * (x[0] - 3.0) + (x[1] + 2.0) * (x[1] + 2.0);
    };
    std::vector<Bounds> bounds = {{-5, 5}, {-5, 5}};

    DifferentialEvolution::Config cfg_no_reg;
    cfg_no_reg.max_generations = 100;
    cfg_no_reg.tol = 1e-10;
    auto result_no_reg = DifferentialEvolution::minimize(f, bounds, cfg_no_reg);

    DifferentialEvolution::Config cfg_reg;
    cfg_reg.max_generations = 100;
    cfg_reg.tol = 1e-10;
    cfg_reg.lambda_reg = 0.0;
    cfg_reg.params_prior = {0.0, 0.0};
    auto result_reg = DifferentialEvolution::minimize(f, bounds, cfg_reg);

    EXPECT_NEAR(result_reg.x[0], result_no_reg.x[0], 1e-3);
    EXPECT_NEAR(result_reg.x[1], result_no_reg.x[1], 1e-3);
}

// --- DE 早停: 当原始目标 < threshold^2 时停止 ---
TEST(DEEarlyStop, StopsWhenObjectiveBelowThreshold) {
    // f(x) = x^2, 真值 0, early_stop_rmse=0.1 → f < 0.01 时停止
    ObjectiveFn f = [](const std::vector<Real>& x) -> Real {
        return x[0] * x[0] + x[1] * x[1];
    };
    std::vector<Bounds> bounds = {{-5, 5}, {-5, 5}};

    DifferentialEvolution::Config cfg;
    cfg.max_generations = 200;
    cfg.tol = 1e-12;  // 严格收敛容差, 确保早停先于 stagnation
    cfg.early_stop_rmse = 0.1;  // f < 0.01 时停止
    auto result = DifferentialEvolution::minimize(f, bounds, cfg);
    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.message.find("early_stop"), 0u)
        << "expected early_stop message, got: " << result.message;
}

// --- SABR 校准器端到端正则化测试 ---
TEST(CalibRegularizationE2E, SABRRegularizationPullsTowardPrior) {
    // 生成 SABR 合成数据 (beta 固定模式), 然后用正则化校准
    // 验证: 当数据噪声大时, 正则化使参数向先验偏移
    Real F = 100.0, r = 0.0, q = 0.0;
    SABRParams true_p{0.20, 0.5, 0.3, -0.2};
    std::vector<Real> strikes = {85.0, 95.0, 100.0, 105.0, 115.0};
    Real T = 0.5;
    std::vector<MarketQuote> quotes;
    for (Real K : strikes) {
        Real iv = detail::sabr_implied_vol_hagan(F * std::exp(r * T), K, T, true_p);
        quotes.push_back({K, T, 0.0, iv, 0.0});
    }

    // 无正则化校准
    SABRCalibrator cal_no_reg;
    cal_no_reg.set_market(F, r, q);
    cal_no_reg.set_fixed_beta(0.5);
    CalibConfig cfg_no_reg;
    cfg_no_reg.use_de_init = false;
    cfg_no_reg.lm_max_iter = 200;
    auto result_no_reg = cal_no_reg.calibrate(quotes, cfg_no_reg);

    // 有正则化校准 (prior 远离真值, lambda 较大)
    SABRCalibrator cal_reg;
    cal_reg.set_market(F, r, q);
    cal_reg.set_fixed_beta(0.5);
    CalibConfig cfg_reg;
    cfg_reg.use_de_init = false;
    cfg_reg.lm_max_iter = 200;
    cfg_reg.lambda_reg = 100.0;
    // SABR 固定 beta 模式: params = [alpha, nu, rho]
    // 先验设为远离真值的值
    cfg_reg.params_prior = {0.5, 1.0, 0.5};  // 远离 true (0.2, 0.3, -0.2)
    auto result_reg = cal_reg.calibrate(quotes, cfg_reg);

    // 正则化结果应向先验偏移: alpha_reg 比 alpha_no_reg 更接近先验 0.5
    Real alpha_no_reg = result_no_reg.params[0];
    Real alpha_reg = result_reg.params[0];
    EXPECT_LT(std::abs(alpha_reg - 0.5), std::abs(alpha_no_reg - 0.5))
        << "regularized alpha should be closer to prior 0.5";
}

// --- CalibConfig 新字段默认值验证 ---
TEST(CalibConfigDefaults, NewFieldsHaveCorrectDefaults) {
    CalibConfig cfg;
    EXPECT_DOUBLE_EQ(cfg.lambda_reg, 0.0);
    EXPECT_TRUE(cfg.params_prior.empty());
    EXPECT_DOUBLE_EQ(cfg.early_stop_rmse, 0.0);
}

// ===========================================================================
// SLSQP tests (ADR-018 implementation boundary)
// 12 test cases covering P0 basic functionality, GARCH application, and
// numerical stability. References:
//   - scipy.optimize.minimize(method='SLSQP') for numerical baselines
//   - arch package GARCH constraint conventions
// ===========================================================================

// P0: Unconstrained quadratic minimization
// min (x-1)^2 + (y-2)^2,  solution x=1, y=2
TEST(M3CompileCheck, SLSQPUnconstrainedQuadratic) {
    ObjectiveFn f = [](const std::vector<Real>& x) {
        Real d1 = x[0] - 1.0;
        Real d2 = x[1] - 2.0;
        return d1 * d1 + d2 * d2;
    };
    std::vector<Bounds> bounds = {{-10.0, 10.0}, {-10.0, 10.0}};
    SLSQP::Config cfg;
    cfg.max_iterations = 100;
    auto r = SLSQP::minimize(f, {0.0, 0.0}, bounds, {}, {}, cfg);
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.x[0], 1.0, 1e-4);
    EXPECT_NEAR(r.x[1], 2.0, 1e-4);
    EXPECT_NEAR(r.fx, 0.0, 1e-8);
}

// P0: Bound-constrained quadratic
// min (x-5)^2 s.t. x in [0, 2],  solution x=2 (boundary)
TEST(M3CompileCheck, SLSQPBoundConstraint) {
    ObjectiveFn f = [](const std::vector<Real>& x) {
        Real d = x[0] - 5.0;
        return d * d;
    };
    std::vector<Bounds> bounds = {{0.0, 2.0}};
    auto r = SLSQP::minimize(f, {1.0}, bounds);
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.x[0], 2.0, 1e-4);
    EXPECT_NEAR(r.fx, 9.0, 1e-4);
}

// P0: Inequality constraint c(x) >= 0
// min (x-3)^2 s.t. x >= 1  (constraint: x - 1 >= 0), solution x=3 (interior)
TEST(M3CompileCheck, SLSQPInequalityConstraint) {
    ObjectiveFn f = [](const std::vector<Real>& x) {
        Real d = x[0] - 3.0;
        return d * d;
    };
    std::vector<Bounds> bounds = {{-10.0, 10.0}};
    // Inequality: x - 1 >= 0
    ConstraintFn ineq = [](const std::vector<Real>& x) {
        return std::vector<Real>{x[0] - 1.0};
    };
    auto r = SLSQP::minimize(f, {0.5}, bounds, {ineq}, {});
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.x[0], 3.0, 1e-4);
    // Verify constraint satisfied
    EXPECT_GE(r.x[0] - 1.0, -1e-6);
}

// P0: Equality constraint c(x) = 0
// min x^2 + y^2 s.t. x + y = 2,  solution x=1, y=1
TEST(M3CompileCheck, SLSQPEqualityConstraint) {
    ObjectiveFn f = [](const std::vector<Real>& x) {
        return x[0] * x[0] + x[1] * x[1];
    };
    std::vector<Bounds> bounds = {{-10.0, 10.0}, {-10.0, 10.0}};
    // Equality: x + y - 2 = 0
    ConstraintFn eq = [](const std::vector<Real>& x) {
        return std::vector<Real>{x[0] + x[1] - 2.0};
    };
    SLSQP::Config cfg;
    cfg.max_iterations = 200;
    auto r = SLSQP::minimize(f, {0.5, 0.5}, bounds, {}, {eq}, cfg);
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.x[0], 1.0, 1e-3);
    EXPECT_NEAR(r.x[1], 1.0, 1e-3);
    EXPECT_NEAR(r.x[0] + r.x[1], 2.0, 1e-4);
}

// P0: Mixed constraints (equality + inequality + bounds)
// min (x-2)^2 + (y-1)^2 s.t. x + y <= 2 (ineq: 2 - x - y >= 0), x - y = 1 (eq)
// Solution: x=1.5, y=0.5 (active inequality + equality)
TEST(M3CompileCheck, SLSQPMixedConstraints) {
    ObjectiveFn f = [](const std::vector<Real>& x) {
        Real d1 = x[0] - 2.0;
        Real d2 = x[1] - 1.0;
        return d1 * d1 + d2 * d2;
    };
    std::vector<Bounds> bounds = {{-10.0, 10.0}, {-10.0, 10.0}};
    ConstraintFn ineq = [](const std::vector<Real>& x) {
        return std::vector<Real>{2.0 - x[0] - x[1]};
    };
    ConstraintFn eq = [](const std::vector<Real>& x) {
        return std::vector<Real>{x[0] - x[1] - 1.0};
    };
    SLSQP::Config cfg;
    cfg.max_iterations = 200;
    auto r = SLSQP::minimize(f, {0.0, 0.0}, bounds, {ineq}, {eq}, cfg);
    EXPECT_TRUE(r.converged);
    // x + y = 2 (active inequality), x - y = 1 => x=1.5, y=0.5
    EXPECT_NEAR(r.x[0], 1.5, 1e-3);
    EXPECT_NEAR(r.x[1], 0.5, 1e-3);
    EXPECT_GE(2.0 - r.x[0] - r.x[1], -1e-6);
    EXPECT_NEAR(r.x[0] - r.x[1], 1.0, 1e-4);
}

// P0: GARCH parameter constraints (non-negativity + stationarity)
// Constraints: omega > 0, alpha >= 0, beta >= 0, alpha + beta < 1
// We test constraint satisfaction, not actual GARCH calibration.
TEST(M3CompileCheck, SLSQPGARCHConstraints) {
    // Dummy objective: minimize distance to (0.1, 0.05, 0.9)
    // But alpha + beta = 0.95 < 1, so unconstrained optimum is feasible.
    ObjectiveFn f = [](const std::vector<Real>& x) {
        Real d1 = x[0] - 0.1;
        Real d2 = x[1] - 0.05;
        Real d3 = x[2] - 0.9;
        return d1 * d1 + d2 * d2 + d3 * d3;
    };
    std::vector<Bounds> bounds = {
        {1e-8, 10.0},   // omega > 0
        {0.0, 10.0},    // alpha >= 0
        {0.0, 10.0}     // beta >= 0
    };
    // Stationarity: 1 - alpha - beta >= 0  => alpha + beta <= 1
    ConstraintFn stationarity = [](const std::vector<Real>& x) {
        return std::vector<Real>{1.0 - x[1] - x[2]};
    };
    auto r = SLSQP::minimize(f, {0.05, 0.02, 0.5}, bounds, {stationarity}, {});
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.x[0], 0.1, 1e-3);
    EXPECT_NEAR(r.x[1], 0.05, 1e-3);
    EXPECT_NEAR(r.x[2], 0.9, 1e-3);
    // Constraints satisfied
    EXPECT_GT(r.x[0], 0.0);
    EXPECT_GE(r.x[1], 0.0);
    EXPECT_GE(r.x[2], 0.0);
    EXPECT_GE(1.0 - r.x[1] - r.x[2], -1e-6);
}

// P0: GARCH-style calibration with infeasible start
// Start at (0.3, 0.7, 0.5) where alpha+beta=1.2 > 1 (infeasible)
// Optimizer must recover to feasible region
TEST(M3CompileCheck, SLSQPGARCHCalibration) {
    // Objective: mimic GARCH negative log-likelihood shape (quadratic approx)
    // Target: omega=0.1, alpha=0.1, beta=0.85 (alpha+beta=0.95 < 1)
    ObjectiveFn f = [](const std::vector<Real>& x) {
        Real d1 = x[0] - 0.1;
        Real d2 = x[1] - 0.1;
        Real d3 = x[2] - 0.85;
        return d1 * d1 + d2 * d2 + d3 * d3;
    };
    std::vector<Bounds> bounds = {
        {1e-8, 10.0},
        {0.0, 10.0},
        {0.0, 10.0}
    };
    ConstraintFn stationarity = [](const std::vector<Real>& x) {
        return std::vector<Real>{1.0 - x[1] - x[2]};
    };
    SLSQP::Config cfg;
    cfg.max_iterations = 200;
    // Infeasible start: alpha + beta = 0.7 + 0.5 = 1.2 > 1
    auto r = SLSQP::minimize(f, {0.3, 0.7, 0.5}, bounds, {stationarity}, {}, cfg);
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.x[0], 0.1, 1e-3);
    EXPECT_NEAR(r.x[1], 0.1, 1e-3);
    EXPECT_NEAR(r.x[2], 0.85, 1e-3);
    EXPECT_GE(1.0 - r.x[1] - r.x[2], -1e-6);
}

// P1: Infeasible start recovery
// min x^2 s.t. x >= 5, start at x=0 (infeasible)
TEST(M3CompileCheck, SLSQPInfeasibleStart) {
    ObjectiveFn f = [](const std::vector<Real>& x) {
        return x[0] * x[0];
    };
    std::vector<Bounds> bounds = {{-100.0, 100.0}};
    ConstraintFn ineq = [](const std::vector<Real>& x) {
        return std::vector<Real>{x[0] - 5.0};  // x >= 5
    };
    SLSQP::Config cfg;
    cfg.max_iterations = 200;
    auto r = SLSQP::minimize(f, {0.0}, bounds, {ineq}, {}, cfg);
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.x[0], 5.0, 1e-3);
    EXPECT_GE(r.x[0] - 5.0, -1e-6);
}

// P1: Rosenbrock function (classic unconstrained test)
// min (1-x)^2 + 100*(y - x^2)^2,  solution x=1, y=1
TEST(M3CompileCheck, SLSQPRosenbrock) {
    ObjectiveFn f = [](const std::vector<Real>& x) {
        Real d1 = 1.0 - x[0];
        Real d2 = x[1] - x[0] * x[0];
        return d1 * d1 + 100.0 * d2 * d2;
    };
    std::vector<Bounds> bounds = {{-5.0, 5.0}, {-5.0, 5.0}};
    SLSQP::Config cfg;
    cfg.max_iterations = 500;
    cfg.ftol = 1e-10;
    cfg.gtol = 1e-8;
    auto r = SLSQP::minimize(f, {-1.2, 1.0}, bounds, {}, {}, cfg);
    // Rosenbrock is a hard test; relax tolerance
    EXPECT_NEAR(r.x[0], 1.0, 1e-2);
    EXPECT_NEAR(r.x[1], 1.0, 1e-2);
}

// P1: Max iterations reached (non-converging case)
TEST(M3CompileCheck, SLSQPMaxIterations) {
    // Highly oscillatory objective that won't converge
    Size eval_count = 0;
    ObjectiveFn f = [&eval_count](const std::vector<Real>& x) {
        ++eval_count;
        return std::sin(x[0] * 100.0) * std::cos(x[0] * 100.0) + x[0] * x[0];
    };
    std::vector<Bounds> bounds = {{-2.0, 2.0}};
    SLSQP::Config cfg;
    cfg.max_iterations = 3;  // very low to force max-iter
    auto r = SLSQP::minimize(f, {0.5}, bounds, {}, {}, cfg);
    // Should not crash; either converged or max iterations
    EXPECT_LE(r.n_iterations, 3);
}

// P1: Numerical baseline vs scipy SLSQP (basic quadratic)
// scipy: minimize(lambda x: (x[0]-1.5)**2 + (x[1]-2.5)**2, [0,0],
//                 method='SLSQP', bounds=[(-5,5),(-5,5)])
// Expected: x=[1.5, 2.5], fun=0
TEST(M3CompileCheck, SLSQPVsScipyBasic) {
    ObjectiveFn f = [](const std::vector<Real>& x) {
        Real d1 = x[0] - 1.5;
        Real d2 = x[1] - 2.5;
        return d1 * d1 + d2 * d2;
    };
    std::vector<Bounds> bounds = {{-5.0, 5.0}, {-5.0, 5.0}};
    auto r = SLSQP::minimize(f, {0.0, 0.0}, bounds);
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.x[0], 1.5, 1e-4);
    EXPECT_NEAR(r.x[1], 2.5, 1e-4);
    EXPECT_NEAR(r.fx, 0.0, 1e-8);
}

// P1: Numerical baseline vs scipy SLSQP (constrained GARCH-like)
// scipy: minimize(lambda x: (x[0]-0.1)**2 + (x[1]-0.1)**2 + (x[2]-0.85)**2,
//                 [0.3, 0.7, 0.5], method='SLSQP',
//                 bounds=[(1e-8,10),(0,10),(0,10)],
//                 constraints={'type':'ineq','fun': lambda x: 1-x[1]-x[2]})
// Expected: x=[0.1, 0.1, 0.85], fun=0
TEST(M3CompileCheck, SLSQPVsScipyGARCH) {
    ObjectiveFn f = [](const std::vector<Real>& x) {
        Real d1 = x[0] - 0.1;
        Real d2 = x[1] - 0.1;
        Real d3 = x[2] - 0.85;
        return d1 * d1 + d2 * d2 + d3 * d3;
    };
    std::vector<Bounds> bounds = {
        {1e-8, 10.0},
        {0.0, 10.0},
        {0.0, 10.0}
    };
    ConstraintFn stationarity = [](const std::vector<Real>& x) {
        return std::vector<Real>{1.0 - x[1] - x[2]};
    };
    SLSQP::Config cfg;
    cfg.max_iterations = 200;
    auto r = SLSQP::minimize(f, {0.3, 0.7, 0.5}, bounds, {stationarity}, {}, cfg);
    EXPECT_TRUE(r.converged);
    // GARCH parameter-level tolerance: 1e-4 (vs scipy baseline)
    EXPECT_NEAR(r.x[0], 0.1, 1e-4);
    EXPECT_NEAR(r.x[1], 0.1, 1e-4);
    EXPECT_NEAR(r.x[2], 0.85, 1e-4);
    EXPECT_GE(1.0 - r.x[1] - r.x[2], -1e-6);
}
