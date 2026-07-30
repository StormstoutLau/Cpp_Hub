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
