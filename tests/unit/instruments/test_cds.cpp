// v1.2 Batch 4: 信用违约互换 (CDS) 单元测试
// 覆盖: CreditCurve / CDSConfig / CreditDefaultSwap / CDSCurveBootstrapper
#include <gtest/gtest.h>
#include "cpphub/instruments/credit/credit_curve.hpp"
#include "cpphub/instruments/credit/cds.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

using namespace cpphub;

// ============================================================
// 辅助: 构造平坦零息曲线 (单一零利率)
// ============================================================
namespace {
ZeroCurve flat_zero_curve(Real rate, Real max_T = 30.0) {
    std::vector<Real> Ts;
    std::vector<Real> rs;
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
// CreditCurve 测试
// ============================================================
TEST(CreditCurveTest, FromFlatCDSSpread) {
    Real s = 0.02, R = 0.40;
    auto cc = CreditCurve::from_flat_cds(s, R);
    // h = s / (1-R) = 0.02/0.6 = 0.0333...
    Real h_expected = s / (1.0 - R);
    EXPECT_NEAR(cc.hazard_rate(5.0), h_expected, 1e-12);
    EXPECT_NEAR(cc.survival_prob(5.0), std::exp(-h_expected * 5.0), 1e-12);
    EXPECT_NEAR(cc.recovery_rate(), R, 0.0);
    EXPECT_NEAR(cc.lgd(), 1.0 - R, 1e-15);
}

TEST(CreditCurveTest, ExpectedLoss) {
    // h = 0.02 flat, R = 0.4 → LGD = 0.6, PD(0,5) = 1 - exp(-0.1)
    auto cc = CreditCurve::flat(0.02, 0.40);
    Real pd5 = 1.0 - std::exp(-0.02 * 5.0);
    EXPECT_NEAR(cc.expected_loss(5.0), 0.6 * pd5, 1e-12);
    // 区间预期损失
    Real pd_1_5 = std::exp(-0.02 * 1.0) - std::exp(-0.02 * 5.0);
    EXPECT_NEAR(cc.expected_loss(1.0, 5.0), 0.6 * pd_1_5, 1e-12);
}

TEST(CreditCurveTest, InvalidRecoveryRate) {
    auto pd = PDCurve::flat(0.02);
    EXPECT_THROW(CreditCurve(pd, -0.1), std::invalid_argument);
    EXPECT_THROW(CreditCurve(pd, 1.0), std::invalid_argument);
    EXPECT_THROW(CreditCurve(pd, 1.5), std::invalid_argument);
}

TEST(CreditCurveTest, ZeroHazardRateNoDefault) {
    auto cc = CreditCurve::flat(0.0, 0.40);
    EXPECT_NEAR(cc.survival_prob(100.0), 1.0, 1e-15);
    EXPECT_NEAR(cc.default_prob(100.0), 0.0, 1e-15);
    EXPECT_NEAR(cc.expected_loss(100.0), 0.0, 1e-15);
}

// ============================================================
// CDSConfig / Schedule 生成测试
// ============================================================
TEST(CDSConfigTest, GenerateScheduleQuarterly) {
    CDSConfig cfg;
    cfg.notional = 1.0;
    cfg.spread = 0.01;
    cfg.start_time = 0.0;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.generate_schedule();
    EXPECT_EQ(cfg.payment_times.size(), 20u);
    EXPECT_EQ(cfg.year_fractions.size(), 20u);
    EXPECT_NEAR(cfg.payment_times[0], 0.25, 1e-12);
    EXPECT_NEAR(cfg.payment_times[19], 5.0, 1e-12);
    EXPECT_NEAR(cfg.year_fractions[0], 0.25, 1e-12);
    // 严格递增
    for (Size i = 1; i < cfg.payment_times.size(); ++i) {
        EXPECT_GT(cfg.payment_times[i], cfg.payment_times[i - 1]);
    }
}

TEST(CDSConfigTest, GenerateScheduleCustom) {
    // 外部指定 payment_times (不应被覆盖)
    CDSConfig cfg;
    cfg.notional = 1.0;
    cfg.spread = 0.01;
    cfg.start_time = 0.0;
    cfg.maturity = 3.0;
    cfg.n_premiums = 6;
    cfg.payment_times = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
    cfg.year_fractions = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    cfg.generate_schedule();
    EXPECT_EQ(cfg.payment_times.size(), 6u);
    EXPECT_NEAR(cfg.payment_times[0], 0.5, 1e-15);
}

TEST(CDSConfigTest, ZeroPremiumsThrows) {
    CDSConfig cfg;
    cfg.n_premiums = 0;
    EXPECT_THROW(cfg.generate_schedule(), std::invalid_argument);
}

TEST(CreditDefaultSwapTest, InvalidParameters) {
    CDSConfig cfg;
    cfg.notional = 0.0;  // 非法
    cfg.spread = 0.01;
    cfg.maturity = 5.0;
    cfg.n_premiums = 20;
    cfg.payment_times.clear();
    cfg.year_fractions.clear();
    EXPECT_THROW(CreditDefaultSwap cds(cfg), std::invalid_argument);

    cfg.notional = 1.0;
    cfg.spread = -0.01;  // 非法
    EXPECT_THROW(CreditDefaultSwap cds(cfg), std::invalid_argument);

    cfg.spread = 0.01;
    cfg.maturity = 0.0;  // maturity <= start_time
    EXPECT_THROW(CreditDefaultSwap cds(cfg), std::invalid_argument);
}

// ============================================================
// CDS 定价 - 解析对比 (r=0, h flat)
// ============================================================
TEST(CDSPricingTest, PremiumLegAnalyticRZeroHFlat) {
    // r=0 (P=1), h=0.02 flat, R=0.4, N=1, s=0.01, 5Y quarterly
    // PV_premium = s * tau * Σ (Q_i + 0.5 * (Q_{i-1} - Q_i)) = s * tau * Σ 0.5*(Q_{i-1}+Q_i)
    auto discount = flat_zero_curve(0.0);
    auto cc = CreditCurve::flat(0.02, 0.40);
    auto cfg = make_cds(1.0, 0.01, 5.0, 4, true);
    CreditDefaultSwap cds(cfg);

    Real s = 0.01, tau = 0.25;
    Real expected = 0.0;
    for (Size i = 1; i <= 20; ++i) {
        Real t_prev = static_cast<Real>(i - 1) * tau;
        Real t_i = static_cast<Real>(i) * tau;
        Real Q_prev = std::exp(-0.02 * t_prev);
        Real Q_i = std::exp(-0.02 * t_i);
        expected += s * tau * 0.5 * (Q_prev + Q_i);
    }
    EXPECT_NEAR(cds.premium_leg_pv(cc, discount), expected, 1e-12);
}

TEST(CDSPricingTest, ProtectionLegAnalyticRZeroHFlat) {
    // r=0, h=0.02 flat, R=0.4 → LGD=0.6
    // PV_protection = LGD * Σ PD(t_{i-1}, t_i) * P(t_mid) = 0.6 * Σ (Q_{i-1} - Q_i) * 1
    auto discount = flat_zero_curve(0.0);
    auto cc = CreditCurve::flat(0.02, 0.40);
    auto cfg = make_cds(1.0, 0.01, 5.0, 4, true);
    CreditDefaultSwap cds(cfg);

    Real tau = 0.25;
    Real expected = 0.0;
    for (Size i = 1; i <= 20; ++i) {
        Real t_prev = static_cast<Real>(i - 1) * tau;
        Real t_i = static_cast<Real>(i) * tau;
        Real Q_prev = std::exp(-0.02 * t_prev);
        Real Q_i = std::exp(-0.02 * t_i);
        expected += 0.6 * (Q_prev - Q_i);
    }
    EXPECT_NEAR(cds.protection_leg_pv(cc, discount), expected, 1e-12);
}

TEST(CDSPricingTest, RiskyPV01Analytic) {
    // Risky_pv01 = Σ tau * (Q_i + 0.5*(Q_{i-1}-Q_i)) = tau * Σ 0.5*(Q_{i-1}+Q_i)
    auto discount = flat_zero_curve(0.0);
    auto cc = CreditCurve::flat(0.02, 0.40);
    auto cfg = make_cds(1.0, 0.01, 5.0, 4, true);
    CreditDefaultSwap cds(cfg);

    Real tau = 0.25;
    Real expected = 0.0;
    for (Size i = 1; i <= 20; ++i) {
        Real t_prev = static_cast<Real>(i - 1) * tau;
        Real t_i = static_cast<Real>(i) * tau;
        Real Q_prev = std::exp(-0.02 * t_prev);
        Real Q_i = std::exp(-0.02 * t_i);
        expected += tau * 0.5 * (Q_prev + Q_i);
    }
    EXPECT_NEAR(cds.risky_pv01(cc, discount), expected, 1e-12);
}

TEST(CDSPricingTest, PVIdentity) {
    // PV (buyer) = PV_protection - PV_premium
    auto discount = flat_zero_curve(0.05);
    auto cc = CreditCurve::flat(0.02, 0.40);
    auto cfg = make_cds(1.0, 0.01, 5.0, 4, true);
    CreditDefaultSwap cds(cfg);
    Real prot = cds.protection_leg_pv(cc, discount);
    Real prem = cds.premium_leg_pv(cc, discount);
    EXPECT_NEAR(cds.pv(cc, discount), prot - prem, 1e-12);
}

TEST(CDSPricingTest, BuyerSellerSignFlip) {
    auto discount = flat_zero_curve(0.05);
    auto cc = CreditCurve::flat(0.02, 0.40);
    auto cfg_buyer = make_cds(1.0, 0.01, 5.0, 4, true);
    auto cfg_seller = make_cds(1.0, 0.01, 5.0, 4, false);
    CreditDefaultSwap cds_b(cfg_buyer);
    CreditDefaultSwap cds_s(cfg_seller);
    Real pv_b = cds_b.pv(cc, discount);
    Real pv_s = cds_s.pv(cc, discount);
    EXPECT_NEAR(pv_b, -pv_s, 1e-12);
}

TEST(CDSPricingTest, ParSpreadRelation) {
    // par_spread = PV_protection / Risky_pv01
    // 用此 spread 重构 CDS, PV 应 ≈ 0
    auto discount = flat_zero_curve(0.05);
    auto cc = CreditCurve::flat(0.02, 0.40);
    auto cfg = make_cds(1.0, 0.01, 5.0, 4, true);
    CreditDefaultSwap cds(cfg);
    Real s_par = cds.par_spread(cc, discount);
    Real prot = cds.protection_leg_pv(cc, discount);
    Real rpv01 = cds.risky_pv01(cc, discount);
    EXPECT_NEAR(s_par, prot / rpv01, 1e-12);

    // 用 s_par 作为 spread, PV 应 ≈ 0
    auto cfg_par = make_cds(1.0, s_par, 5.0, 4, true);
    CreditDefaultSwap cds_par(cfg_par);
    EXPECT_NEAR(cds_par.pv(cc, discount), 0.0, 1e-10);
}

TEST(CDSPricingTest, ParSpreadApproxFromFlatCDS) {
    // from_flat_cds: h = s/(1-R) 单点近似 (忽略折现/应计)
    // 真实 par_spread 应接近 s 但不完全相等
    // 当 h 小 (s 小), par_spread → s (线性近似成立)
    Real s = 0.001;  // 10bp, 小 spread
    Real R = 0.40;
    auto discount = flat_zero_curve(0.03);
    auto cc = CreditCurve::from_flat_cds(s, R);
    auto cfg = make_cds(1.0, s, 5.0, 4, true);
    CreditDefaultSwap cds(cfg);
    Real s_par = cds.par_spread(cc, discount);
    // 小 spread 近似: |s_par - s| / s < 5%
    Real rel_err = std::abs(s_par - s) / s;
    EXPECT_LT(rel_err, 0.05);
}

TEST(CDSPricingTest, ZeroHazardRateNoDefault) {
    // h=0 → 无违约, protection leg = 0, PV (buyer) = -premium leg ≤ 0
    auto discount = flat_zero_curve(0.05);
    auto cc = CreditCurve::flat(0.0, 0.40);
    auto cfg = make_cds(1.0, 0.01, 5.0, 4, true);
    CreditDefaultSwap cds(cfg);
    EXPECT_NEAR(cds.protection_leg_pv(cc, discount), 0.0, 1e-15);
    EXPECT_NEAR(cds.par_spread(cc, discount), 0.0, 1e-15);
    // buyer PV = -premium (因为付保费但无赔付)
    Real prem = cds.premium_leg_pv(cc, discount);
    EXPECT_NEAR(cds.pv(cc, discount), -prem, 1e-12);
    EXPECT_LT(cds.pv(cc, discount), 0.0);
}

TEST(CDSPricingTest, HigherSpreadReducesBuyerPV) {
    // spread 越高 (付越多), buyer PV 越低
    auto discount = flat_zero_curve(0.05);
    auto cc = CreditCurve::flat(0.02, 0.40);
    auto cfg_low = make_cds(1.0, 0.005, 5.0, 4, true);
    auto cfg_high = make_cds(1.0, 0.05, 5.0, 4, true);
    CreditDefaultSwap cds_low(cfg_low);
    CreditDefaultSwap cds_high(cfg_high);
    EXPECT_GT(cds_low.pv(cc, discount), cds_high.pv(cc, discount));
}

TEST(CDSPricingTest, NotionalScaling) {
    // PV 与 notional 线性相关
    auto discount = flat_zero_curve(0.05);
    auto cc = CreditCurve::flat(0.02, 0.40);
    auto cfg1 = make_cds(1.0, 0.01, 5.0, 4, true);
    auto cfg10 = make_cds(10.0, 0.01, 5.0, 4, true);
    CreditDefaultSwap cds1(cfg1);
    CreditDefaultSwap cds10(cfg10);
    EXPECT_NEAR(cds10.pv(cc, discount), 10.0 * cds1.pv(cc, discount), 1e-10);
}

// ============================================================
// CDS Curve Bootstrap 测试
// ============================================================
TEST(CDSBootstrapTest, SingleQuoteSelfConsistency) {
    // 1 个 5Y quote, bootstrap 后用该曲线对该 CDS 定价, PV ≈ 0
    Real s = 0.02, R = 0.40;
    auto discount = flat_zero_curve(0.03);
    std::vector<CDSCurveBootstrapper::CDSQuote> quotes = {
        {5.0, s}
    };
    CDSCurveBootstrapper boot(quotes, R, discount, 4);
    auto cc = boot.bootstrap();

    // 用该曲线对 5Y CDS 定价, PV 应 ≈ 0
    auto cfg = make_cds(1.0, s, 5.0, 4, true);
    CreditDefaultSwap cds(cfg);
    Real pv = cds.pv(cc, discount);
    EXPECT_NEAR(pv, 0.0, 1e-8);

    // par_spread 应 ≈ s
    Real s_par = cds.par_spread(cc, discount);
    EXPECT_NEAR(s_par, s, 1e-8);
}

TEST(CDSBootstrapTest, MultiQuoteSelfConsistency) {
    // 多 quote bootstrap, 每个 quote 期限的 CDS PV ≈ 0
    Real R = 0.40;
    auto discount = flat_zero_curve(0.03);
    std::vector<CDSCurveBootstrapper::CDSQuote> quotes = {
        {1.0, 0.01},
        {3.0, 0.015},
        {5.0, 0.02},
        {7.0, 0.025},
        {10.0, 0.03}
    };
    CDSCurveBootstrapper boot(quotes, R, discount, 4);
    auto cc = boot.bootstrap();

    // 对每个 quote 期限的 CDS, PV ≈ 0, par_spread ≈ quote spread
    for (const auto& q : quotes) {
        auto cfg = make_cds(1.0, q.spread, q.maturity, 4, true);
        CreditDefaultSwap cds(cfg);
        Real pv = cds.pv(cc, discount);
        EXPECT_NEAR(pv, 0.0, 1e-8)
            << "Maturity=" << q.maturity << " spread=" << q.spread;
        Real s_par = cds.par_spread(cc, discount);
        EXPECT_NEAR(s_par, q.spread, 1e-8)
            << "Maturity=" << q.maturity << " spread=" << q.spread;
    }
}

TEST(CDSBootstrapTest, FlatCurveProducesFlatHazard) {
    // 所有 quote spread 相同, bootstrap 后各段 hazard rate 应接近
    Real R = 0.40;
    Real s_flat = 0.02;
    auto discount = flat_zero_curve(0.03);
    std::vector<CDSCurveBootstrapper::CDSQuote> quotes = {
        {1.0, s_flat},
        {3.0, s_flat},
        {5.0, s_flat},
        {7.0, s_flat},
        {10.0, s_flat}
    };
    CDSCurveBootstrapper boot(quotes, R, discount, 4);
    auto cc = boot.bootstrap();

    // 各段 hazard rate 应接近 s_flat / (1-R) = 0.0333...
    Real h_approx = s_flat / (1.0 - R);
    const auto& hrs = cc.pd_curve().hazard_rates();
    for (Size i = 0; i < hrs.size(); ++i) {
        Real rel_err = std::abs(hrs[i] - h_approx) / h_approx;
        EXPECT_LT(rel_err, 0.05)  // 容差 5% (折现与应计引入小偏差
            << "Segment " << i << " h=" << hrs[i];
    }
}

TEST(CDSBootstrapTest, IncreasingSpreadProducesIncreasingHazard) {
    // 期限结构向上倾斜 (长端 spread 大), hazard rate 应单调 (大致) 递增
    Real R = 0.40;
    auto discount = flat_zero_curve(0.03);
    std::vector<CDSCurveBootstrapper::CDSQuote> quotes = {
        {1.0, 0.005},
        {3.0, 0.01},
        {5.0, 0.02},
        {7.0, 0.03},
        {10.0, 0.04}
    };
    CDSCurveBootstrapper boot(quotes, R, discount, 4);
    auto cc = boot.bootstrap();

    // 平均 hazard rate 应随期限递增 (检查累积 PD 的凸性)
    Real pd_1 = cc.default_prob(1.0);
    Real pd_5 = cc.default_prob(5.0);
    Real pd_10 = cc.default_prob(10.0);
    EXPECT_GT(pd_5, pd_1);
    EXPECT_GT(pd_10, pd_5);

    // 长端瞬时 hazard 应大于短端
    Real h_short = cc.hazard_rate(0.5);
    Real h_long = cc.hazard_rate(9.5);
    EXPECT_GT(h_long, h_short);
}

TEST(CDSBootstrapTest, InvalidParameters) {
    auto discount = flat_zero_curve(0.03);
    // 空 quotes
    EXPECT_THROW(CDSCurveBootstrapper({}, 0.4, discount, 4), std::invalid_argument);
    // 非正 maturity
    std::vector<CDSCurveBootstrapper::CDSQuote> bad1 = {{0.0, 0.01}};
    EXPECT_THROW(CDSCurveBootstrapper(bad1, 0.4, discount, 4), std::invalid_argument);
    // 负 spread
    std::vector<CDSCurveBootstrapper::CDSQuote> bad2 = {{5.0, -0.01}};
    EXPECT_THROW(CDSCurveBootstrapper(bad2, 0.4, discount, 4), std::invalid_argument);
    // 非递增 maturity
    std::vector<CDSCurveBootstrapper::CDSQuote> bad3 = {{5.0, 0.01}, {3.0, 0.02}};
    EXPECT_THROW(CDSCurveBootstrapper(bad3, 0.4, discount, 4), std::invalid_argument);
    // recovery 越界
    std::vector<CDSCurveBootstrapper::CDSQuote> ok_q = {{5.0, 0.01}};
    EXPECT_THROW(CDSCurveBootstrapper(ok_q, -0.1, discount, 4), std::invalid_argument);
    EXPECT_THROW(CDSCurveBootstrapper(ok_q, 1.0, discount, 4), std::invalid_argument);
    // freq = 0
    EXPECT_THROW(CDSCurveBootstrapper(ok_q, 0.4, discount, 0), std::invalid_argument);
}

// ============================================================
// CDS Par Spread 与 from_flat_cds 近似关系
// ============================================================
TEST(CDSPricingTest, ParSpreadFromFlatCDSApproximation) {
    // 经典近似: s ≈ h * (1-R), 当 h 小且利率低时近似好
    // 检查 |par_spread - s_input| / s_input < 1% (小 spread)
    Real R = 0.40;
    auto discount = flat_zero_curve(0.01);  // 低利率减小折现效应
    Real s = 0.001;  // 10bp
    auto cc = CreditCurve::from_flat_cds(s, R);
    auto cfg = make_cds(1.0, s, 5.0, 4, true);
    CreditDefaultSwap cds(cfg);
    Real s_par = cds.par_spread(cc, discount);
    Real rel_err = std::abs(s_par - s) / s;
    EXPECT_LT(rel_err, 0.01) << "s_input=" << s << " s_par=" << s_par;
}

// ============================================================
// 端到端: bootstrap + CDS 定价 + par spread 验证
// ============================================================
TEST(CDSEndToEndTest, BootstrappedCurvePricesInputQuotes) {
    // 真实场景: 从市场 CDS quotes bootstrap 曲线, 再用该曲线定价
    Real R = 0.40;
    auto discount = flat_zero_curve(0.025);  // 2.5% OIS
    std::vector<CDSCurveBootstrapper::CDSQuote> market_quotes = {
        {1.0, 0.0080},   // 80bp
        {2.0, 0.0120},   // 120bp
        {3.0, 0.0160},   // 160bp
        {5.0, 0.0220},   // 220bp
        {7.0, 0.0270},   // 270bp
        {10.0, 0.0320},  // 320bp
    };
    CDSCurveBootstrapper boot(market_quotes, R, discount, 4);
    auto cc = boot.bootstrap();

    // 验证 1: 每个 market quote 的 CDS PV ≈ 0
    for (const auto& q : market_quotes) {
        auto cfg = make_cds(10.0, q.spread, q.maturity, 4, true);
        CreditDefaultSwap cds(cfg);
        Real pv = cds.pv(cc, discount);
        EXPECT_NEAR(pv, 0.0, 1e-7)
            << "Maturity=" << q.maturity << " spread=" << q.spread;
    }

    // 验证 2: 用该曲线定价一个非 quote 期限的 CDS (如 4Y)
    // 应得到介于 3Y (160bp) 和 5Y (220bp) 之间的 par spread
    auto cfg_4y = make_cds(1.0, 0.018, 4.0, 4, true);  // spread 任意, 用 par_spread
    CreditDefaultSwap cds_4y(cfg_4y);
    Real s_par_4y = cds_4y.par_spread(cc, discount);
    // 4Y par spread 应介于 3Y 和 5Y quote 之间
    EXPECT_GT(s_par_4y, 0.0160);
    EXPECT_LT(s_par_4y, 0.0220);
}

TEST(CDSEndToEndTest, CDSOnBootstrappedCurveMatchesInputSpread) {
    // 单点验证: bootstrap 单 quote 后, par_spread 精确恢复输入
    Real R = 0.40;
    Real s_input = 0.0150;  // 150bp
    auto discount = flat_zero_curve(0.02);
    std::vector<CDSCurveBootstrapper::CDSQuote> quotes = {{5.0, s_input}};
    CDSCurveBootstrapper boot(quotes, R, discount, 4);
    auto cc = boot.bootstrap();

    auto cfg = make_cds(1.0, s_input, 5.0, 4, true);
    CreditDefaultSwap cds(cfg);
    Real s_par = cds.par_spread(cc, discount);
    EXPECT_NEAR(s_par, s_input, 1e-9);

    // buyer PV (spread = par) ≈ 0
    EXPECT_NEAR(cds.pv(cc, discount), 0.0, 1e-10);

    // 若 spread > par, buyer PV < 0 (付太多)
    auto cfg_high = make_cds(1.0, s_input + 0.005, 5.0, 4, true);
    CreditDefaultSwap cds_high(cfg_high);
    EXPECT_LT(cds_high.pv(cc, discount), 0.0);

    // 若 spread < par, buyer PV > 0 (保护便宜)
    auto cfg_low = make_cds(1.0, s_input - 0.005, 5.0, 4, true);
    CreditDefaultSwap cds_low(cfg_low);
    EXPECT_GT(cds_low.pv(cc, discount), 0.0);
}
