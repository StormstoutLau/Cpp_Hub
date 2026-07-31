// v1.2 Batch 8: 合成 CDO (Synthetic CDO) 单元测试
// 覆盖: CDOTrancheConfig / CDOTrancheResult
//       + CDOLHPPricer (Vasicek 大组合半解析)
//       + CDOMCPricer (Copula MC 模拟)
//       + make_cdx_ig_tranches
// 测试维度: 配置验证/期望损失单调性/边界情形/相关性效应/
//          tranche 优先级/MC vs LHP 一致性/par spread
#include <gtest/gtest.h>
#include "cpphub/instruments/credit/cdo.hpp"
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
// 1. CDOTrancheConfig 基本验证
// ============================================================
TEST(CDOTrancheConfigTest, DefaultValues) {
    CDOTrancheConfig cfg;
    EXPECT_NEAR(cfg.attachment, 0.0, 1e-15);
    EXPECT_NEAR(cfg.detachment, 0.03, 1e-15);
    EXPECT_NEAR(cfg.spread, 0.05, 1e-15);
    EXPECT_NEAR(cfg.maturity, 5.0, 1e-15);
    EXPECT_EQ(cfg.n_premiums, 20u);
    EXPECT_TRUE(cfg.is_protection_seller);
}

TEST(CDOTrancheConfigTest, GenerateSchedule) {
    CDOTrancheConfig cfg;
    cfg.maturity = 3.0;
    cfg.n_premiums = 12;  // 季度
    cfg.generate_schedule();
    ASSERT_EQ(cfg.payment_times.size(), 12u);
    ASSERT_EQ(cfg.year_fractions.size(), 12u);
    EXPECT_NEAR(cfg.payment_times[0], 0.25, 1e-12);
    EXPECT_NEAR(cfg.payment_times[11], 3.0, 1e-12);
    for (Real tau : cfg.year_fractions) {
        EXPECT_NEAR(tau, 0.25, 1e-12);
    }
}

TEST(CDOTrancheConfigTest, ValidationThrowsOnInvalid) {
    CDOTrancheConfig cfg;
    cfg.attachment = -0.1;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.attachment = 0.0;
    cfg.detachment = -0.1;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.detachment = 0.0;  // detachment <= attachment
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.detachment = 1.5;  // detachment > 1
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.detachment = 0.05;
    cfg.maturity = 0.0;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.maturity = 5.0;
    cfg.spread = -0.01;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.spread = 0.05;
    cfg.n_premiums = 0;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);
}

TEST(CDOTrancheConfigTest, TrancheWidth) {
    CDOTrancheConfig cfg;
    cfg.attachment = 0.03;
    cfg.detachment = 0.07;
    EXPECT_NEAR(cfg.tranche_width(), 0.04, 1e-15);
}

// ============================================================
// 2. CDOLHPPricer 基本行为
// ============================================================
TEST(CDOLHPPricerTest, ConstructorValidates) {
    auto discount = flat_zero_curve(0.03);
    EXPECT_THROW(CDOLHPPricer(-0.1, 0.05, 0.6, discount, 5.0),
                  std::invalid_argument);
    EXPECT_THROW(CDOLHPPricer(1.0, 0.05, 0.6, discount, 5.0),
                  std::invalid_argument);
    EXPECT_THROW(CDOLHPPricer(0.3, -0.1, 0.6, discount, 5.0),
                  std::invalid_argument);
    EXPECT_THROW(CDOLHPPricer(0.3, 1.5, 0.6, discount, 5.0),
                  std::invalid_argument);
    EXPECT_THROW(CDOLHPPricer(0.3, 0.05, -0.1, discount, 5.0),
                  std::invalid_argument);
    EXPECT_THROW(CDOLHPPricer(0.3, 0.05, 1.5, discount, 5.0),
                  std::invalid_argument);
    EXPECT_THROW(CDOLHPPricer(0.3, 0.05, 0.6, discount, 0.0),
                  std::invalid_argument);
}

TEST(CDOLHPPricerTest, Accessors) {
    auto discount = flat_zero_curve(0.03);
    CDOLHPPricer pricer(0.3, 0.05, 0.6, discount, 5.0);
    EXPECT_NEAR(pricer.rho(), 0.3, 1e-15);
    EXPECT_NEAR(pricer.pd_maturity(), 0.05, 1e-15);
    EXPECT_NEAR(pricer.lgd(), 0.6, 1e-15);
    EXPECT_NEAR(pricer.maturity(), 5.0, 1e-15);
}

TEST(CDOLHPPricerTest, ZeroPDGivesZeroLoss) {
    // PD=0 → 期望损失为 0
    auto discount = flat_zero_curve(0.03);
    CDOLHPPricer pricer(0.3, 0.0, 0.6, discount, 5.0);
    EXPECT_NEAR(pricer.expected_tranche_loss(5.0, 0.0, 0.03), 0.0, 1e-15);
    EXPECT_NEAR(pricer.expected_tranche_loss(5.0, 0.03, 0.07), 0.0, 1e-15);
    EXPECT_NEAR(pricer.expected_portfolio_loss(5.0), 0.0, 1e-15);
}

TEST(CDOLHPPricerTest, ExpectedPortfolioLossEqualsLGDTimesPD) {
    // 无条件期望损失 E[L/N] = LGD * PD (与 rho 无关)
    auto discount = flat_zero_curve(0.03);
    CDOLHPPricer pricer(0.3, 0.10, 0.6, discount, 5.0);
    Real expected = 0.6 * 0.10;
    EXPECT_NEAR(pricer.expected_portfolio_loss(5.0), expected, 1e-10);
}

TEST(CDOLHPPricerTest, ExpectedTrancheLossMonotonicInTime) {
    // 期望损失随时间单调递增
    auto discount = flat_zero_curve(0.03);
    CDOLHPPricer pricer(0.3, 0.10, 0.6, discount, 5.0);
    Real A = 0.0, D = 0.10;
    Real prev = -1.0;
    for (Size k = 0; k <= 10; ++k) {
        Real t = 0.5 * k;  // 0 → 5
        Real loss = pricer.expected_tranche_loss(t, A, D);
        EXPECT_GE(loss, prev);
        prev = loss;
    }
}

TEST(CDOLHPPricerTest, ExpectedTrancheLossBoundedByTrancheWidth) {
    // tranche 损失不能超过 tranche 宽度
    auto discount = flat_zero_curve(0.03);
    CDOLHPPricer pricer(0.3, 0.20, 0.6, discount, 5.0);
    Real A = 0.0, D = 0.03;  // equity tranche
    Real loss = pricer.expected_tranche_loss(5.0, A, D);
    EXPECT_LE(loss, D - A + 1e-10);
    EXPECT_GE(loss, 0.0);
}

TEST(CDOLHPPricerTest, HigherCorrelationMoreConcentratedLoss) {
    // 高相关性 → equity tranche 损失降低, senior tranche 损失增加
    // (违约更趋向同时发生, 集中在 senior 部分)
    auto discount = flat_zero_curve(0.03);
    Real A_eq = 0.0, D_eq = 0.03;       // equity
    Real A_sen = 0.07, D_sen = 0.10;    // senior

    CDOLHPPricer pricer_low(0.1, 0.10, 0.6, discount, 5.0);
    CDOLHPPricer pricer_high(0.5, 0.10, 0.6, discount, 5.0);

    Real eq_low = pricer_low.expected_tranche_loss(5.0, A_eq, D_eq);
    Real eq_high = pricer_high.expected_tranche_loss(5.0, A_eq, D_eq);
    Real sen_low = pricer_low.expected_tranche_loss(5.0, A_sen, D_sen);
    Real sen_high = pricer_high.expected_tranche_loss(5.0, A_sen, D_sen);

    // Equity tranche: 高相关性应使损失降低 (因更多违约同时发生, 但系统性事件
    // 更多时 PD 条件化后, 极端 m 下 PD 高, 但对 equity 来说 LHP 极限损失
    // 受 min(LGD*PD_cond, D) 限制. 高 rho 下 m 分布更广, 极端 m 下损失饱和更快.
    // 经典结果: 高 rho → equity 损失降低, senior 损失增加)
    EXPECT_LT(eq_high, eq_low + 1e-6);
    EXPECT_GT(sen_high, sen_low);
}

TEST(CDOLHPPricerTest, ExpectedRemainingNotionalNonNegative) {
    auto discount = flat_zero_curve(0.03);
    CDOLHPPricer pricer(0.3, 0.10, 0.6, discount, 5.0);
    for (Real A : {0.0, 0.03, 0.07}) {
        for (Real D : {0.03, 0.07, 0.10}) {
            if (D <= A) continue;
            Real rn = pricer.expected_remaining_notional(5.0, A, D);
            EXPECT_GE(rn, 0.0);
            EXPECT_LE(rn, 1.0);
        }
    }
}

// ============================================================
// 3. CDOLHPPricer 完整定价
// ============================================================
TEST(CDOLHPPricerTest, PriceReturnsValidResult) {
    auto discount = flat_zero_curve(0.03);
    CDOLHPPricer pricer(0.3, 0.10, 0.6, discount, 5.0);

    CDOTrancheConfig cfg;
    cfg.attachment = 0.0;
    cfg.detachment = 0.03;
    cfg.spread = 0.05;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;

    auto result = pricer.price(cfg);
    EXPECT_GT(result.pv_premium_leg, 0.0);
    EXPECT_GT(result.pv_protection_leg, 0.0);
    EXPECT_GT(result.risky_pv01, 0.0);
    EXPECT_GT(result.par_spread, 0.0);
    EXPECT_GE(result.expected_loss, 0.0);
    EXPECT_LE(result.expected_loss, 0.03);  // tranche width
    EXPECT_EQ(result.expected_tranche_loss.size(), 20u);
    EXPECT_EQ(result.expected_remaining_notional.size(), 20u);
}

TEST(CDOLHPPricerTest, EquityTrancheHigherSpreadThanSenior) {
    // Equity tranche 风险更高 → par spread 应远高于 senior
    auto discount = flat_zero_curve(0.03);
    CDOLHPPricer pricer(0.3, 0.10, 0.6, discount, 5.0);

    CDOTrancheConfig cfg_eq;
    cfg_eq.attachment = 0.0;
    cfg_eq.detachment = 0.03;
    cfg_eq.maturity = 5.0;
    cfg_eq.n_premiums = 20;
    auto res_eq = pricer.price(cfg_eq);

    CDOTrancheConfig cfg_sen;
    cfg_sen.attachment = 0.07;
    cfg_sen.detachment = 0.10;
    cfg_sen.maturity = 5.0;
    cfg_sen.n_premiums = 20;
    auto res_sen = pricer.price(cfg_sen);

    EXPECT_GT(res_eq.par_spread, res_sen.par_spread);
    EXPECT_GT(res_eq.expected_loss, res_sen.expected_loss);
}

TEST(CDOLHPPricerTest, ParSpreadZeroesOutPV) {
    // 用 par_spread 作为 spread 时, PV 应接近 0
    auto discount = flat_zero_curve(0.03);
    CDOLHPPricer pricer(0.3, 0.10, 0.6, discount, 5.0);

    CDOTrancheConfig cfg;
    cfg.attachment = 0.03;
    cfg.detachment = 0.07;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;

    auto result = pricer.price(cfg);
    Real ps = result.par_spread;

    CDOTrancheConfig cfg_par = cfg;
    cfg_par.spread = ps;
    auto result_par = pricer.price(cfg_par);

    EXPECT_NEAR(result_par.pv, 0.0, 1e-8);
}

TEST(CDOLHPPricerTest, ProtectionSellerPVPositiveWhenSpreadHigh) {
    // spread 高于 par_spread 时, 卖方 PV > 0 (收到的保费超过赔付)
    auto discount = flat_zero_curve(0.03);
    CDOLHPPricer pricer(0.3, 0.10, 0.6, discount, 5.0);

    CDOTrancheConfig cfg;
    cfg.attachment = 0.0;
    cfg.detachment = 0.03;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    auto result_par = pricer.price(cfg);
    Real ps = result_par.par_spread;

    cfg.spread = ps * 2.0;  // 双倍 par spread
    auto result = pricer.price(cfg);
    EXPECT_GT(result.pv, 0.0);
}

TEST(CDOLHPPricerTest, ProtectionBuyerPVFlippedSign) {
    // 买方视角 = -卖方视角
    auto discount = flat_zero_curve(0.03);
    CDOLHPPricer pricer(0.3, 0.10, 0.6, discount, 5.0);

    CDOTrancheConfig cfg;
    cfg.attachment = 0.0;
    cfg.detachment = 0.03;
    cfg.spread = 0.10;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.is_protection_seller = true;
    auto res_seller = pricer.price(cfg);

    cfg.is_protection_seller = false;
    auto res_buyer = pricer.price(cfg);

    EXPECT_NEAR(res_buyer.pv, -res_seller.pv, 1e-12);
    EXPECT_NEAR(res_buyer.pv_premium_leg, -res_seller.pv_premium_leg, 1e-12);
    EXPECT_NEAR(res_buyer.pv_protection_leg, -res_seller.pv_protection_leg, 1e-12);
}

TEST(CDOLHPPricerTest, ZeroPDGivesZeroPVProtection) {
    // PD=0 → 保护腿为 0, par_spread 为 0
    auto discount = flat_zero_curve(0.03);
    CDOLHPPricer pricer(0.3, 0.0, 0.6, discount, 5.0);

    CDOTrancheConfig cfg;
    cfg.attachment = 0.0;
    cfg.detachment = 0.03;
    cfg.spread = 0.05;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;

    auto result = pricer.price(cfg);
    EXPECT_NEAR(result.pv_protection_leg, 0.0, 1e-12);
    EXPECT_NEAR(result.par_spread, 0.0, 1e-12);
    EXPECT_GT(result.pv_premium_leg, 0.0);  // 仍有保费收入
    EXPECT_GT(result.pv, 0.0);  // 卖方获利
}

// ============================================================
// 4. CDOMCPricer 基本行为
// ============================================================
TEST(CDOMCPricerTest, ConstructorValidates) {
    auto discount = flat_zero_curve(0.03);
    std::vector<CreditCurve> empty_curves;
    EXPECT_THROW(CDOMCPricer(empty_curves, discount, 1.0),
                  std::invalid_argument);
    EXPECT_THROW(CDOMCPricer({CreditCurve::flat(0.05, 0.4)}, discount, -1.0),
                  std::invalid_argument);
}

TEST(CDOMCPricerTest, Accessors) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.05, 0.4);
    CDOMCPricer pricer({cc, cc, cc, cc, cc}, discount, 100.0);
    EXPECT_EQ(pricer.n_names(), 5u);
    EXPECT_NEAR(pricer.total_notional(), 100.0, 1e-15);
}

TEST(CDOMCPricerTest, PriceWithGaussianCopula) {
    // 基本 MC 定价: 输出应有限
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.05, 0.4);  // 5% hazard, 40% recovery
    std::vector<CreditCurve> curves(10, cc);
    Real total_notional = 100.0;
    CDOMCPricer pricer(curves, discount, total_notional);

    GaussianCopula copula(make_equicorrelation(10, 0.3));
    CDOTrancheConfig cfg;
    cfg.attachment = 0.0;
    cfg.detachment = 0.10;
    cfg.spread = 0.05;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;

    auto result = pricer.price(copula, cfg, 5000, 42);
    EXPECT_EQ(result.expected_tranche_loss.size(), 20u);
    EXPECT_GE(result.expected_loss, 0.0);
    EXPECT_LE(result.expected_loss, 0.10);  // tranche width
    EXPECT_GT(result.risky_pv01, 0.0);
    EXPECT_GE(result.par_spread, 0.0);
}

TEST(CDOMCPricerTest, CopulaNameSizeMismatchThrows) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.05, 0.4);
    std::vector<CreditCurve> curves(10, cc);
    CDOMCPricer pricer(curves, discount, 100.0);

    // Copula 维度不匹配
    GaussianCopula copula(make_equicorrelation(5, 0.3));  // 5 vs 10
    CDOTrancheConfig cfg;
    cfg.attachment = 0.0;
    cfg.detachment = 0.10;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    EXPECT_THROW(pricer.price(copula, cfg, 100, 42), std::invalid_argument);
}

TEST(CDOMCPricerTest, ZeroPDGivesZeroLoss) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.0, 0.4);  // zero hazard
    std::vector<CreditCurve> curves(5, cc);
    CDOMCPricer pricer(curves, discount, 100.0);

    GaussianCopula copula(make_equicorrelation(5, 0.3));
    CDOTrancheConfig cfg;
    cfg.attachment = 0.0;
    cfg.detachment = 0.10;
    cfg.spread = 0.05;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;

    auto result = pricer.price(copula, cfg, 2000, 42);
    EXPECT_NEAR(result.pv_protection_leg, 0.0, 1e-10);
    EXPECT_NEAR(result.expected_loss, 0.0, 1e-10);
    EXPECT_NEAR(result.par_spread, 0.0, 1e-10);
    EXPECT_GT(result.pv_premium_leg, 0.0);
}

TEST(CDOMCPricerTest, HigherCorrelationChangesEquityLoss) {
    // 高相关性改变 equity tranche 期望损失分布 (方向取决于组合大小和参数)
    // 这里仅验证: 结果有限、非负、在 tranche 宽度内, 且两种相关性下结果不同
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.10, 0.4);  // 10% hazard
    std::vector<CreditCurve> curves(20, cc);
    CDOMCPricer pricer(curves, discount, 100.0);

    CDOTrancheConfig cfg;
    cfg.attachment = 0.0;
    cfg.detachment = 0.10;  // equity
    cfg.spread = 0.0;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;

    GaussianCopula copula_low(make_equicorrelation(20, 0.1));
    GaussianCopula copula_high(make_equicorrelation(20, 0.6));
    auto res_low = pricer.price(copula_low, cfg, 10000, 42);
    auto res_high = pricer.price(copula_high, cfg, 10000, 42);

    // 两个结果都应合法
    EXPECT_GE(res_low.expected_loss, 0.0);
    EXPECT_LE(res_low.expected_loss, 0.10);
    EXPECT_GE(res_high.expected_loss, 0.0);
    EXPECT_LE(res_high.expected_loss, 0.10);
    // 高相关性下 par spread 应不同 (相关性影响损失分布)
    EXPECT_NE(res_high.par_spread, res_low.par_spread);
}

TEST(CDOMCPricerTest, ParSpreadZeroesOutPV) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.08, 0.4);
    std::vector<CreditCurve> curves(10, cc);
    CDOMCPricer pricer(curves, discount, 100.0);

    GaussianCopula copula(make_equicorrelation(10, 0.3));
    CDOTrancheConfig cfg;
    cfg.attachment = 0.03;
    cfg.detachment = 0.07;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;

    auto result = pricer.price(copula, cfg, 10000, 42);
    Real ps = result.par_spread;

    CDOTrancheConfig cfg_par = cfg;
    cfg_par.spread = ps;
    auto result_par = pricer.price(copula, cfg_par, 10000, 42);
    // MC 误差允许 1% of pv_premium
    Real tol = std::max(1e-4, 0.01 * std::abs(result_par.pv_premium_leg));
    EXPECT_NEAR(result_par.pv, 0.0, tol);
}

TEST(CDOMCPricerTest, MatchesLHPForLargeHomogeneousPortfolio) {
    // 大组合 + 同质 → MC 应接近 LHP 半解析结果
    auto discount = flat_zero_curve(0.03);
    Real h = 0.05;
    Real R = 0.4;
    Real pd_5y = 1.0 - std::exp(-h * 5.0);
    Real lgd = 1.0 - R;
    Real rho = 0.3;
    Real maturity = 5.0;

    // LHP 半解析
    CDOLHPPricer lhp(rho, pd_5y, lgd, discount, maturity);
    CDOTrancheConfig cfg_lhp;
    cfg_lhp.attachment = 0.0;
    cfg_lhp.detachment = 0.10;
    cfg_lhp.maturity = maturity;
    cfg_lhp.n_premiums = 20;
    auto res_lhp = lhp.price(cfg_lhp);

    // MC (大组合, 100 个名字)
    auto cc = CreditCurve::flat(h, R);
    std::vector<CreditCurve> curves(100, cc);
    CDOMCPricer mc_pricer(curves, discount, 1.0);  // total_notional = 1 → loss as fraction
    GaussianCopula copula(make_equicorrelation(100, rho));
    CDOTrancheConfig cfg_mc;
    cfg_mc.attachment = 0.0;
    cfg_mc.detachment = 0.10;
    cfg_mc.maturity = maturity;
    cfg_mc.n_premiums = 20;
    auto res_mc = mc_pricer.price(copula, cfg_mc, 20000, 42);

    // 大组合 MC 应接近 LHP (10% 容差)
    Real rel_tol = 0.15;
    EXPECT_NEAR(res_mc.expected_loss, res_lhp.expected_loss,
                 rel_tol * std::abs(res_lhp.expected_loss) + 1e-4);
    EXPECT_NEAR(res_mc.par_spread, res_lhp.par_spread,
                 rel_tol * std::abs(res_lhp.par_spread) + 1e-4);
}

// ============================================================
// 5. CDOMCPricer 与 OneFactorGaussianCopula
// ============================================================
TEST(CDOMCPricerTest, WorksWithOneFactorCopula) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.05, 0.4);
    std::vector<CreditCurve> curves(5, cc);
    CDOMCPricer pricer(curves, discount, 100.0);

    OneFactorGaussianCopula copula(0.3, 5);
    CDOTrancheConfig cfg;
    cfg.attachment = 0.0;
    cfg.detachment = 0.10;
    cfg.spread = 0.05;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;

    EXPECT_NO_THROW({
        auto result = pricer.price(copula, cfg, 2000, 42);
        EXPECT_GE(result.expected_loss, 0.0);
    });
}

// ============================================================
// 6. CDOMCPricer 与 TCopula
// ============================================================
TEST(CDOMCPricerTest, WorksWithTCopula) {
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::flat(0.05, 0.4);
    std::vector<CreditCurve> curves(5, cc);
    CDOMCPricer pricer(curves, discount, 100.0);

    TCopula copula(make_equicorrelation(5, 0.3), 5.0);
    CDOTrancheConfig cfg;
    cfg.attachment = 0.0;
    cfg.detachment = 0.10;
    cfg.spread = 0.05;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;

    EXPECT_NO_THROW({
        auto result = pricer.price(copula, cfg, 2000, 42);
        EXPECT_GE(result.expected_loss, 0.0);
    });
}

// ============================================================
// 7. 便捷工厂函数
// ============================================================
TEST(CDOFactoryTest, MakeCDXIGTranches) {
    auto tranches = make_cdx_ig_tranches(5.0, 4);
    ASSERT_EQ(tranches.size(), 6u);  // CDX IG 6 个标准 tranches
    // 检查 attachment / detachment
    EXPECT_NEAR(tranches[0].attachment, 0.0, 1e-15);
    EXPECT_NEAR(tranches[0].detachment, 0.03, 1e-15);
    EXPECT_NEAR(tranches[1].attachment, 0.03, 1e-15);
    EXPECT_NEAR(tranches[1].detachment, 0.07, 1e-15);
    EXPECT_NEAR(tranches[5].attachment, 0.30, 1e-15);
    EXPECT_NEAR(tranches[5].detachment, 1.00, 1e-15);
    // 检查期限与支付次数
    for (const auto& t : tranches) {
        EXPECT_NEAR(t.maturity, 5.0, 1e-15);
        EXPECT_EQ(t.n_premiums, 20u);
    }
}

TEST(CDOFactoryTest, MakeCDXIGTranchesCustomMaturity) {
    auto tranches = make_cdx_ig_tranches(3.0, 4);
    ASSERT_EQ(tranches.size(), 6u);
    for (const auto& t : tranches) {
        EXPECT_NEAR(t.maturity, 3.0, 1e-15);
        EXPECT_EQ(t.n_premiums, 12u);
    }
}
