// v1.2 Batch 8: 篮子 CDS (Basket CDS / Nth-to-Default) 单元测试
// 覆盖: BasketCDSConfig / BasketCDSResult
//       + BasketCDSMCPricer (Copula MC 模拟)
//       + NthToDefaultHomogeneousPricer (单因子 Gaussian Copula 半解析)
//       + make_homogeneous_credit_curves / make_heterogeneous_credit_curves
// 测试维度: 配置验证/概率单调性/相关性效应/期限效应/MC vs 半解析一致性/
//          par spread/Greeks 边界
#include <gtest/gtest.h>
#include "cpphub/instruments/credit/basket_cds.hpp"
#include "cpphub/instruments/credit/copula.hpp"
#include "cpphub/instruments/credit/credit_curve.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

using namespace cpphub;

namespace {
// 平坦零息曲线
ZeroCurve flat_zero_curve(Real rate, Real max_T = 30.0) {
    std::vector<Real> Ts, rs;
    for (Size i = 1; i <= 30; ++i) {
        Real T = static_cast<Real>(i);
        if (T > max_T) break;
        Ts.push_back(T);
        rs.push_back(rate);
    }
    return ZeroCurve(Ts, rs, ZeroCurve::InterpType::LinearZero);
}
}  // namespace

// ============================================================
// 1. BasketCDSConfig 基本验证
// ============================================================
TEST(BasketCDSConfigTest, DefaultValues) {
    BasketCDSConfig cfg;
    EXPECT_NEAR(cfg.notional, 1.0, 1e-15);
    EXPECT_NEAR(cfg.spread, 0.01, 1e-15);
    EXPECT_NEAR(cfg.maturity, 5.0, 1e-15);
    EXPECT_EQ(cfg.n_premiums, 20u);
    EXPECT_EQ(cfg.nth_default, 1u);
    EXPECT_TRUE(cfg.is_buyer);
}

TEST(BasketCDSConfigTest, GenerateSchedule) {
    BasketCDSConfig cfg;
    cfg.maturity = 2.0;
    cfg.n_premiums = 8;
    cfg.generate_schedule();
    ASSERT_EQ(cfg.payment_times.size(), 8u);
    ASSERT_EQ(cfg.year_fractions.size(), 8u);
    EXPECT_NEAR(cfg.payment_times[0], 0.25, 1e-12);
    EXPECT_NEAR(cfg.payment_times[7], 2.0, 1e-12);
    for (Real tau : cfg.year_fractions) {
        EXPECT_NEAR(tau, 0.25, 1e-12);
    }
}

TEST(BasketCDSConfigTest, ValidationThrowsOnInvalid) {
    BasketCDSConfig cfg;
    cfg.notional = -1.0;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.notional = 1.0;
    cfg.spread = -0.01;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.spread = 0.01;
    cfg.maturity = 0.0;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.maturity = 5.0;
    cfg.n_premiums = 0;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.n_premiums = 20;
    cfg.nth_default = 0;  // 必须为正
    EXPECT_THROW(cfg.validate(), std::invalid_argument);
}

// ============================================================
// 2. NthToDefaultHomogeneousPricer 基本行为
// ============================================================
TEST(NthToDefaultPricerTest, ConstructorValidates) {
    auto discount = flat_zero_curve(0.03);
    EXPECT_THROW(NthToDefaultHomogeneousPricer(-0.1, 0.05, 0.6, 5, discount, 5.0),
                  std::invalid_argument);
    EXPECT_THROW(NthToDefaultHomogeneousPricer(1.0, 0.05, 0.6, 5, discount, 5.0),
                  std::invalid_argument);
    EXPECT_THROW(NthToDefaultHomogeneousPricer(0.3, -0.1, 0.6, 5, discount, 5.0),
                  std::invalid_argument);
    EXPECT_THROW(NthToDefaultHomogeneousPricer(0.3, 1.5, 0.6, 5, discount, 5.0),
                  std::invalid_argument);
    EXPECT_THROW(NthToDefaultHomogeneousPricer(0.3, 0.05, -0.1, 5, discount, 5.0),
                  std::invalid_argument);
    EXPECT_THROW(NthToDefaultHomogeneousPricer(0.3, 0.05, 1.5, 5, discount, 5.0),
                  std::invalid_argument);
    EXPECT_THROW(NthToDefaultHomogeneousPricer(0.3, 0.05, 0.6, 0, discount, 5.0),
                  std::invalid_argument);
}

TEST(NthToDefaultPricerTest, Accessors) {
    auto discount = flat_zero_curve(0.03);
    NthToDefaultHomogeneousPricer p(0.3, 0.10, 0.6, 10, discount, 5.0);
    EXPECT_NEAR(p.rho(), 0.3, 1e-15);
    EXPECT_EQ(p.n_names(), 10u);
    EXPECT_NEAR(p.pd_at_maturity(), 0.10, 1e-15);
    EXPECT_NEAR(p.lgd(), 0.6, 1e-15);
}

TEST(NthToDefaultPricerTest, ZeroPDGivesZeroProb) {
    auto discount = flat_zero_curve(0.03);
    NthToDefaultHomogeneousPricer p(0.3, 0.0, 0.6, 5, discount, 5.0);
    EXPECT_NEAR(p.prob_nth_default_by(5.0, 1), 0.0, 1e-15);
    EXPECT_NEAR(p.prob_nth_default_by(5.0, 3), 0.0, 1e-15);
    EXPECT_NEAR(p.survival_prob(5.0, 1), 1.0, 1e-15);
}

TEST(NthToDefaultPricerTest, ProbNthDefaultMonotonicInTime) {
    // P(τ_(N) ≤ t) 随 t 单调递增
    auto discount = flat_zero_curve(0.03);
    NthToDefaultHomogeneousPricer p(0.3, 0.10, 0.6, 10, discount, 5.0);
    Real prev = -1.0;
    for (Size k = 0; k <= 10; ++k) {
        Real t = 0.5 * k;  // 0 → 5
        Real prob = p.prob_nth_default_by(t, 2);
        EXPECT_GE(prob, prev);
        prev = prob;
    }
}

TEST(NthToDefaultPricerTest, ProbNthDefaultDecreasesWithN) {
    // 给定 t, P(τ_(N) ≤ t) 随 N 单调递减 (第 N 个违约比第 1 个违约更难)
    auto discount = flat_zero_curve(0.03);
    NthToDefaultHomogeneousPricer p(0.3, 0.10, 0.6, 5, discount, 5.0);
    Real prev = 2.0;
    for (Size N = 1; N <= 5; ++N) {
        Real prob = p.prob_nth_default_by(5.0, N);
        EXPECT_LE(prob, prev);
        prev = prob;
    }
}

TEST(NthToDefaultPricerTest, ProbNthDefaultIncreasesWithPD) {
    // PD 越高, P(τ_(N) ≤ t) 越大
    auto discount = flat_zero_curve(0.03);
    Real prev = -1.0;
    for (Size k = 0; k <= 5; ++k) {
        Real pd = 0.02 + 0.02 * k;
        NthToDefaultHomogeneousPricer p(0.3, pd, 0.6, 5, discount, 5.0);
        Real prob = p.prob_nth_default_by(5.0, 2);
        EXPECT_GE(prob, prev);
        prev = prob;
    }
}

TEST(NthToDefaultPricerTest, ProbNthDefaultBounds) {
    // 概率应在 [0, 1] 之间
    auto discount = flat_zero_curve(0.03);
    NthToDefaultHomogeneousPricer p(0.3, 0.30, 0.6, 5, discount, 5.0);
    for (Size N = 1; N <= 5; ++N) {
        Real prob = p.prob_nth_default_by(5.0, N);
        EXPECT_GE(prob, 0.0);
        EXPECT_LE(prob, 1.0);
    }
    // N > n_names → 概率为 0
    EXPECT_NEAR(p.prob_nth_default_by(5.0, 6), 0.0, 1e-15);
}

TEST(NthToDefaultPricerTest, HigherCorrelationReducesFTDSpread) {
    // FTD: 高相关性 → 多名字同时违约 → 第 1 个违约概率反而降低
    // (相关性强 → 一旦触发, 多个一起违约, 但触发概率降低)
    auto discount = flat_zero_curve(0.03);
    NthToDefaultHomogeneousPricer p_low(0.1, 0.10, 0.6, 10, discount, 5.0);
    NthToDefaultHomogeneousPricer p_high(0.6, 0.10, 0.6, 10, discount, 5.0);

    BasketCDSConfig cfg;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 1;

    auto res_low = p_low.price(cfg);
    auto res_high = p_high.price(cfg);

    EXPECT_LT(res_high.par_spread, res_low.par_spread);
}

TEST(NthToDefaultPricerTest, ParSpreadZeroesOutPV) {
    auto discount = flat_zero_curve(0.03);
    NthToDefaultHomogeneousPricer p(0.3, 0.10, 0.6, 5, discount, 5.0);

    BasketCDSConfig cfg;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 1;
    auto res = p.price(cfg);
    Real ps = res.par_spread;

    BasketCDSConfig cfg_par = cfg;
    cfg_par.spread = ps;
    auto res_par = p.price(cfg_par);
    EXPECT_NEAR(res_par.pv, 0.0, 1e-8);
}

TEST(NthToDefaultPricerTest, PriceReturnsValidResult) {
    auto discount = flat_zero_curve(0.03);
    NthToDefaultHomogeneousPricer p(0.3, 0.10, 0.6, 5, discount, 5.0);

    BasketCDSConfig cfg;
    cfg.notional = 10.0;
    cfg.spread = 0.02;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 2;
    auto res = p.price(cfg);

    EXPECT_EQ(res.survival_prob_at_t.size(), 20u);
    EXPECT_EQ(res.default_prob_at_t.size(), 20u);
    EXPECT_GE(res.prob_nth_default_by_maturity, 0.0);
    EXPECT_LE(res.prob_nth_default_by_maturity, 1.0);
    EXPECT_GT(res.risky_pv01, 0.0);
    EXPECT_GE(res.par_spread, 0.0);
}

TEST(NthToDefaultPricerTest, BuyerSellerPVFlipped) {
    auto discount = flat_zero_curve(0.03);
    NthToDefaultHomogeneousPricer p(0.3, 0.10, 0.6, 5, discount, 5.0);

    BasketCDSConfig cfg;
    cfg.spread = 0.02;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 1;
    cfg.is_buyer = true;
    auto res_buyer = p.price(cfg);

    cfg.is_buyer = false;
    auto res_seller = p.price(cfg);

    EXPECT_NEAR(res_buyer.pv, -res_seller.pv, 1e-12);
}

// ============================================================
// 3. BasketCDSMCPricer 基本行为
// ============================================================
TEST(BasketCDSMCPricerTest, ConstructorValidates) {
    auto discount = flat_zero_curve(0.03);
    std::vector<CreditCurve> empty_curves;
    EXPECT_THROW(BasketCDSMCPricer(empty_curves, discount),
                  std::invalid_argument);
}

TEST(BasketCDSMCPricerTest, Accessors) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.05, 0.4);
    BasketCDSMCPricer pricer({cc, cc, cc}, discount);
    EXPECT_EQ(pricer.n_names(), 3u);
}

TEST(BasketCDSMCPricerTest, PriceWithGaussianCopula) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.05, 0.4);
    std::vector<CreditCurve> curves(5, cc);
    BasketCDSMCPricer pricer(curves, discount);

    GaussianCopula copula(make_equicorrelation(5, 0.3));
    BasketCDSConfig cfg;
    cfg.notional = 10.0;
    cfg.spread = 0.02;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 1;

    auto result = pricer.price(copula, cfg, 10000, 42);
    EXPECT_EQ(result.survival_prob_at_t.size(), 20u);
    EXPECT_EQ(result.default_prob_at_t.size(), 20u);
    EXPECT_GE(result.prob_nth_default_by_maturity, 0.0);
    EXPECT_LE(result.prob_nth_default_by_maturity, 1.0);
    EXPECT_GT(result.risky_pv01, 0.0);
    EXPECT_GE(result.par_spread, 0.0);
}

TEST(BasketCDSMCPricerTest, CopulaNameSizeMismatchThrows) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.05, 0.4);
    std::vector<CreditCurve> curves(5, cc);
    BasketCDSMCPricer pricer(curves, discount);

    GaussianCopula copula(make_equicorrelation(3, 0.3));  // 3 vs 5
    BasketCDSConfig cfg;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 1;
    EXPECT_THROW(pricer.price(copula, cfg, 100, 42), std::invalid_argument);
}

TEST(BasketCDSMCPricerTest, NthDefaultExceedsNNamesThrows) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.05, 0.4);
    std::vector<CreditCurve> curves(5, cc);
    BasketCDSMCPricer pricer(curves, discount);

    GaussianCopula copula(make_equicorrelation(5, 0.3));
    BasketCDSConfig cfg;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 6;  // > 5
    EXPECT_THROW(pricer.price(copula, cfg, 100, 42), std::invalid_argument);
}

TEST(BasketCDSMCPricerTest, ZeroPDGivesZeroProtection) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.0, 0.4);  // zero hazard
    std::vector<CreditCurve> curves(5, cc);
    BasketCDSMCPricer pricer(curves, discount);

    GaussianCopula copula(make_equicorrelation(5, 0.3));
    BasketCDSConfig cfg;
    cfg.spread = 0.02;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 1;

    auto result = pricer.price(copula, cfg, 5000, 42);
    EXPECT_NEAR(result.pv_protection_leg, 0.0, 1e-10);
    EXPECT_NEAR(result.prob_nth_default_by_maturity, 0.0, 1e-10);
    EXPECT_NEAR(result.par_spread, 0.0, 1e-10);
    EXPECT_GT(result.pv_premium_leg, 0.0);
    EXPECT_GT(result.risky_pv01, 0.0);
}

TEST(BasketCDSMCPricerTest, ParSpreadZeroesOutPV) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.08, 0.4);
    std::vector<CreditCurve> curves(5, cc);
    BasketCDSMCPricer pricer(curves, discount);

    GaussianCopula copula(make_equicorrelation(5, 0.3));
    BasketCDSConfig cfg;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 1;

    auto result = pricer.price(copula, cfg, 20000, 42);
    Real ps = result.par_spread;

    BasketCDSConfig cfg_par = cfg;
    cfg_par.spread = ps;
    auto result_par = pricer.price(copula, cfg_par, 20000, 42);
    // MC 误差允许
    Real tol = std::max(1e-4, 0.02 * std::abs(result_par.pv_premium_leg));
    EXPECT_NEAR(result_par.pv, 0.0, tol);
}

TEST(BasketCDSMCPricerTest, FTDHigherSpreadThanSecondToDefault) {
    // 第 1 个违约比第 2 个违约更易触发 → FTD par spread 应高于 second
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.05, 0.4);
    std::vector<CreditCurve> curves(5, cc);
    BasketCDSMCPricer pricer(curves, discount);

    GaussianCopula copula(make_equicorrelation(5, 0.3));
    BasketCDSConfig cfg;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 1;
    auto res_ftd = pricer.price(copula, cfg, 20000, 42);

    cfg.nth_default = 2;
    auto res_std = pricer.price(copula, cfg, 20000, 42);

    EXPECT_GT(res_ftd.par_spread, res_std.par_spread);
}

TEST(BasketCDSMCPricerTest, HigherPDGivesHigherSpread) {
    // PD 越高, par spread 越高
    auto discount = flat_zero_curve(0.03);
    GaussianCopula copula(make_equicorrelation(5, 0.3));
    BasketCDSConfig cfg;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 1;

    auto cc_low = CreditCurve::flat(0.03, 0.4);
    std::vector<CreditCurve> curves_low(5, cc_low);
    BasketCDSMCPricer pricer_low(curves_low, discount);
    auto res_low = pricer_low.price(copula, cfg, 20000, 42);

    auto cc_high = CreditCurve::flat(0.10, 0.4);
    std::vector<CreditCurve> curves_high(5, cc_high);
    BasketCDSMCPricer pricer_high(curves_high, discount);
    auto res_high = pricer_high.price(copula, cfg, 20000, 42);

    EXPECT_GT(res_high.par_spread, res_low.par_spread);
}

TEST(BasketCDSMCPricerTest, BuyerSellerPVFlipped) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.05, 0.4);
    std::vector<CreditCurve> curves(5, cc);
    BasketCDSMCPricer pricer(curves, discount);

    GaussianCopula copula(make_equicorrelation(5, 0.3));
    BasketCDSConfig cfg;
    cfg.spread = 0.02;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 1;
    cfg.is_buyer = true;
    auto res_buyer = pricer.price(copula, cfg, 5000, 42);

    cfg.is_buyer = false;
    auto res_seller = pricer.price(copula, cfg, 5000, 42);

    EXPECT_NEAR(res_buyer.pv, -res_seller.pv, 1e-10);
}

// ============================================================
// 4. MC vs 半解析一致性 (同质组合)
// ============================================================
TEST(BasketCDSMCvsSemiTest, MCMatchesSemiAnalyticForHomogeneous) {
    // 同质组合 + 单因子 Gaussian Copula
    // MC 模拟结果应接近半解析结果
    auto discount = flat_zero_curve(0.03);
    Real h = 0.05;
    Real R = 0.4;
    Real pd_5y = 1.0 - std::exp(-h * 5.0);
    Real lgd = 1.0 - R;
    Real rho = 0.3;
    Size n_names = 5;
    Real maturity = 5.0;

    NthToDefaultHomogeneousPricer semi(rho, pd_5y, lgd, n_names, discount, maturity);
    BasketCDSConfig cfg;
    cfg.maturity = maturity;
    cfg.n_premiums = 20;
    cfg.nth_default = 1;
    auto res_semi = semi.price(cfg);

    auto cc = CreditCurve::flat(h, R);
    std::vector<CreditCurve> curves(n_names, cc);
    BasketCDSMCPricer mc(curves, discount);
    OneFactorGaussianCopula copula(rho, n_names);
    auto res_mc = mc.price(copula, cfg, 50000, 42);

    // MC 应接近半解析 (10% 相对容差 + MC 误差)
    Real rel_tol = 0.10;
    EXPECT_NEAR(res_mc.prob_nth_default_by_maturity,
                 res_semi.prob_nth_default_by_maturity,
                 rel_tol * res_semi.prob_nth_default_by_maturity + 1e-3);
    EXPECT_NEAR(res_mc.par_spread, res_semi.par_spread,
                 rel_tol * res_semi.par_spread + 1e-4);
}

TEST(BasketCDSMCvsSemiTest, MCMatchesForSecondToDefault) {
    // 第 2 个违约情形
    auto discount = flat_zero_curve(0.03);
    Real h = 0.05;
    Real R = 0.4;
    Real pd_5y = 1.0 - std::exp(-h * 5.0);
    Real lgd = 1.0 - R;
    Real rho = 0.3;
    Size n_names = 5;
    Real maturity = 5.0;

    NthToDefaultHomogeneousPricer semi(rho, pd_5y, lgd, n_names, discount, maturity);
    BasketCDSConfig cfg;
    cfg.maturity = maturity;
    cfg.n_premiums = 20;
    cfg.nth_default = 2;
    auto res_semi = semi.price(cfg);

    auto cc = CreditCurve::flat(h, R);
    std::vector<CreditCurve> curves(n_names, cc);
    BasketCDSMCPricer mc(curves, discount);
    OneFactorGaussianCopula copula(rho, n_names);
    auto res_mc = mc.price(copula, cfg, 50000, 42);

    Real rel_tol = 0.15;
    EXPECT_NEAR(res_mc.prob_nth_default_by_maturity,
                 res_semi.prob_nth_default_by_maturity,
                 rel_tol * res_semi.prob_nth_default_by_maturity + 1e-3);
}

// ============================================================
// 5. BasketCDSMCPricer 与 t-Copula
// ============================================================
TEST(BasketCDSMCPricerTest, WorksWithTCopula) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.05, 0.4);
    std::vector<CreditCurve> curves(5, cc);
    BasketCDSMCPricer pricer(curves, discount);

    TCopula copula(make_equicorrelation(5, 0.3), 5.0);
    BasketCDSConfig cfg;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.nth_default = 1;

    EXPECT_NO_THROW({
        auto result = pricer.price(copula, cfg, 5000, 42);
        EXPECT_GE(result.par_spread, 0.0);
    });
}

// ============================================================
// 6. 便捷工厂函数
// ============================================================
TEST(BasketCDSFactoryTest, MakeHomogeneousCreditCurves) {
    auto curves = make_homogeneous_credit_curves(5, 0.05, 0.4);
    ASSERT_EQ(curves.size(), 5u);
    for (const auto& cc : curves) {
        EXPECT_NEAR(cc.hazard_rate(1.0), 0.05, 1e-12);
        EXPECT_NEAR(cc.recovery_rate(), 0.4, 1e-15);
    }
}

TEST(BasketCDSFactoryTest, MakeHeterogeneousCreditCurves) {
    std::vector<Real> hrs = {0.02, 0.04, 0.06, 0.08, 0.10};
    auto curves = make_heterogeneous_credit_curves(hrs, 0.4);
    ASSERT_EQ(curves.size(), 5u);
    for (Size i = 0; i < 5; ++i) {
        EXPECT_NEAR(curves[i].hazard_rate(1.0), hrs[i], 1e-12);
        EXPECT_NEAR(curves[i].recovery_rate(), 0.4, 1e-15);
    }
}
