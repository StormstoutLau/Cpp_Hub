// v1.2 Batch 6: Total Return Swap (TRS) 单元测试
// 覆盖: TRSConfig/TRSResult + price_trs (确定性远期定价) + trs_par_spread
// 测试维度: 基本定价/远期价格/信用损失/par spread/零违约/对称性/边界
#include <gtest/gtest.h>
#include "cpphub/instruments/credit/trs.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include "cpphub/instruments/credit/credit_curve.hpp"
#include <cmath>
#include <vector>

using namespace cpphub;

namespace {
// 平坦零息曲线辅助函数
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
// 1. TRSConfig 基本验证
// ============================================================
TEST(TRSConfigTest, DefaultValues) {
    TRSConfig cfg;
    EXPECT_EQ(cfg.notional, 1.0);
    EXPECT_EQ(cfg.S0, 100.0);
    EXPECT_EQ(cfg.asset_yield, 0.0);
    EXPECT_EQ(cfg.funding_spread, 0.005);
    EXPECT_EQ(cfg.maturity, 5.0);
    EXPECT_EQ(cfg.n_payments, 20u);
    EXPECT_TRUE(cfg.is_receiver);
}

TEST(TRSConfigTest, GenerateSchedule) {
    TRSConfig cfg;
    cfg.maturity = 2.0;
    cfg.n_payments = 8;  // 季度
    cfg.generate_schedule();
    ASSERT_EQ(cfg.payment_times.size(), 8u);
    EXPECT_NEAR(cfg.payment_times[0], 0.25, 1e-10);
    EXPECT_NEAR(cfg.payment_times[7], 2.0, 1e-10);
}

TEST(TRSConfigTest, ValidationThrowsOnInvalid) {
    TRSConfig cfg;
    cfg.notional = -1.0;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.notional = 1.0;
    cfg.S0 = -100.0;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.S0 = 100.0;
    cfg.maturity = 0.0;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);
}

// ============================================================
// 2. 零违约 TRS (无信用风险) — 资产远期回报 ≈ 融资成本
// ============================================================
TEST(TRSPriceTest, ZeroDefaultRiskGivesForwardReturnMinusFunding) {
    // 无信用风险: hazard = 0, recovery = 0
    // 远期回报 ≈ r*tau, 融资 ≈ (r + spread)*tau
    // PV ≈ -spread * risky_annuity (receiver 支付 spread)
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.0, 30.0, 100);  // zero hazard

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.S0 = 100.0;
    cfg.asset_yield = 0.0;
    cfg.funding_spread = 0.01;  // 100bp
    cfg.maturity = 5.0;
    cfg.n_payments = 20;  // 季度

    auto result = price_trs(cfg, discount, credit);

    // 无信用风险: PV_asset_return ≈ r * annuity, PV_funding ≈ (r + spread) * annuity
    // PV ≈ -spread * annuity
    // 粗略检查: PV 应为负 (receiver 支付 spread), 且量级 ~ spread * N * T
    EXPECT_LT(result.pv, 0.0);
    // spread=1%, N=100, T=5: PV ≈ -0.01 * 100 * 5 = -5 (粗略)
    EXPECT_GT(result.pv, -10.0);
    EXPECT_LT(result.pv, -1.0);
    // 信用损失应为 0
    EXPECT_NEAR(result.pv_credit_loss, 0.0, 1e-10);
}

TEST(TRSPriceTest, ZeroSpreadGivesZeroPVWithZeroDefault) {
    // spread=0 且无违约: PV ≈ 0 (资产回报 = 融资成本)
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.0, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.S0 = 100.0;
    cfg.asset_yield = 0.0;
    cfg.funding_spread = 0.0;  // 无利差
    cfg.maturity = 5.0;
    cfg.n_payments = 20;

    auto result = price_trs(cfg, discount, credit);
    // 远期回报 ≈ r*tau, 融资 ≈ r*tau, PV ≈ 0
    EXPECT_NEAR(result.pv, 0.0, 0.5);  // 离散化误差
}

// ============================================================
// 3. 信用风险对 PV 的影响
// ============================================================
TEST(TRSPriceTest, CreditRiskReducesReceiverPV) {
    // 增加信用风险 → receiver PV 降低 (信用损失增加)
    auto discount = flat_zero_curve(0.05);

    CreditCurve credit_safe = CreditCurve::flat(0.0, 0.4, 30.0, 100);     // 无违约
    CreditCurve credit_risky = CreditCurve::flat(0.02, 0.4, 30.0, 100);   // 2% hazard

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.S0 = 100.0;
    cfg.asset_yield = 0.0;
    cfg.funding_spread = 0.01;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;

    auto result_safe = price_trs(cfg, discount, credit_safe);
    auto result_risky = price_trs(cfg, discount, credit_risky);

    // 有信用风险: receiver 额外承担信用损失, PV 应更低
    EXPECT_LT(result_risky.pv, result_safe.pv);
    // 信用损失项应为负 (receiver 支出)
    EXPECT_LT(result_risky.pv_credit_loss, 0.0);
    EXPECT_NEAR(result_safe.pv_credit_loss, 0.0, 1e-10);
}

TEST(TRSPriceTest, HigherHazardGivesMoreCreditLoss) {
    auto discount = flat_zero_curve(0.05);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;
    cfg.funding_spread = 0.01;

    Real prev_loss = 0.0;
    for (Real h : {0.01, 0.02, 0.05, 0.10}) {
        CreditCurve credit = CreditCurve::flat(h, 0.4, 30.0, 100);
        auto result = price_trs(cfg, discount, credit);
        EXPECT_LT(result.pv_credit_loss, prev_loss);  // 更负
        prev_loss = result.pv_credit_loss;
    }
}

TEST(TRSPriceTest, HigherRecoveryReducesCreditLoss) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit_low_R = CreditCurve::flat(0.03, 0.20, 30.0, 100);
    CreditCurve credit_high_R = CreditCurve::flat(0.03, 0.80, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;

    auto result_low_R = price_trs(cfg, discount, credit_low_R);
    auto result_high_R = price_trs(cfg, discount, credit_high_R);

    // 高回收率 → 信用损失更小 (绝对值)
    EXPECT_LT(result_low_R.pv_credit_loss, result_high_R.pv_credit_loss);
}

// ============================================================
// 4. 资产 yield 对 PV 的影响
// ============================================================
TEST(TRSPriceTest, AssetYieldApproximatelyCancelsInPV) {
    // 远期价格 F(t) = S0*exp((r-q)*t) 已含 yield 调整:
    //   total_return = exp((r-q)*tau) - 1 + q*tau ≈ r*tau (一阶与 q 无关)
    // 二阶效应: d(total_return)/dq = tau*(1 - exp((r-q)*tau)) < 0 (当 r>q)
    // 因此 asset_yield 上升 → PV 略降 (而非上升)
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;
    cfg.funding_spread = 0.01;

    cfg.asset_yield = 0.0;
    auto result_no_yield = price_trs(cfg, discount, credit);

    cfg.asset_yield = 0.03;  // 3% yield
    auto result_with_yield = price_trs(cfg, discount, credit);

    // 一阶效应相消, 二阶效应使 PV 略降 (r > q 时)
    EXPECT_LT(result_with_yield.pv, result_no_yield.pv);
    // 变化量应较小 (远小于 yield * N * T = 0.03 * 100 * 5 = 15)
    EXPECT_GT(std::abs(result_with_yield.pv - result_no_yield.pv), 0.0);
    EXPECT_LT(std::abs(result_with_yield.pv - result_no_yield.pv), 1.0);
}

// ============================================================
// 5. Receiver vs Payer 对称性
// ============================================================
TEST(TRSPriceTest, PayerIsNegativeOfReceiver) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;
    cfg.funding_spread = 0.01;

    cfg.is_receiver = true;
    auto result_receiver = price_trs(cfg, discount, credit);

    cfg.is_receiver = false;
    auto result_payer = price_trs(cfg, discount, credit);

    EXPECT_NEAR(result_payer.pv, -result_receiver.pv, 1e-10);
    EXPECT_NEAR(result_payer.pv_asset_leg, -result_receiver.pv_asset_leg, 1e-10);
    EXPECT_NEAR(result_payer.pv_funding_leg, -result_receiver.pv_funding_leg, 1e-10);
}

// ============================================================
// 6. Notional 线性性
// ============================================================
TEST(TRSPriceTest, NotionalLinearity) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    TRSConfig cfg;
    cfg.S0 = 100.0;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;
    cfg.funding_spread = 0.01;
    cfg.notional = 100.0;

    auto result_100 = price_trs(cfg, discount, credit);

    cfg.notional = 200.0;
    auto result_200 = price_trs(cfg, discount, credit);

    EXPECT_NEAR(result_200.pv, 2.0 * result_100.pv, 1e-8);
    EXPECT_NEAR(result_200.pv_asset_leg, 2.0 * result_100.pv_asset_leg, 1e-8);
    EXPECT_NEAR(result_200.pv_funding_leg, 2.0 * result_100.pv_funding_leg, 1e-8);
}

// ============================================================
// 7. Par Spread 测试
// ============================================================
TEST(TRSParSpreadTest, ZeroDefaultGivesZeroParSpread) {
    // 无信用风险 + asset_yield=0 → par_spread ≈ 0
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.0, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.asset_yield = 0.0;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;

    Real ps = trs_par_spread(cfg, discount, credit);
    // 无违约, 无 yield: 资产远期回报 = 融资成本, par_spread ≈ 0
    EXPECT_NEAR(ps, 0.0, 1e-6);
}

TEST(TRSParSpreadTest, CreditRiskGivesNegativeParSpread) {
    // 约定: receiver 承担参考资产信用损失 (见 CreditRiskReducesReceiverPV)
    // 有信用风险 → pv_asset 含信用损失项 (负) → par_spread = (pv_asset - pv_funding_base)/annuity < 0
    // par_spread < 0 表示 receiver 应获得信用补偿 (即支付负 spread)
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.asset_yield = 0.0;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;

    Real ps = trs_par_spread(cfg, discount, credit);
    EXPECT_LT(ps, 0.0);
    // 粗略: hazard=2%, LGD=0.6 → 年化预期损失 ≈ 1.2%, par_spread ≈ -1.2%
    EXPECT_LT(ps, -0.005);
    EXPECT_GT(ps, -0.025);
}

TEST(TRSParSpreadTest, ParSpreadZeroesOutPV) {
    // 用 par_spread 作为 funding_spread, PV 应 ≈ 0
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.asset_yield = 0.0;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;

    Real ps = trs_par_spread(cfg, discount, credit);
    cfg.funding_spread = ps;
    auto result = price_trs(cfg, discount, credit);
    EXPECT_NEAR(result.pv, 0.0, 1e-8);
}

TEST(TRSParSpreadTest, HigherHazardGivesMoreNegativeParSpread) {
    // receiver 承担信用损失: hazard 越高 → 信用损失越大 → par_spread 越负
    auto discount = flat_zero_curve(0.05);
    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;

    Real prev_ps = 0.0;  // h=0 时 par_spread ≈ 0
    for (Real h : {0.01, 0.02, 0.05, 0.10}) {
        CreditCurve credit = CreditCurve::flat(h, 0.4, 30.0, 100);
        Real ps = trs_par_spread(cfg, discount, credit);
        EXPECT_LT(ps, prev_ps);  // 更负
        prev_ps = ps;
    }
}

TEST(TRSParSpreadTest, AssetYieldReducesParSpread) {
    // 资产 yield 增加 → receiver 收入更多 → 需要的 par_spread 降低
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;

    cfg.asset_yield = 0.0;
    Real ps_no_yield = trs_par_spread(cfg, discount, credit);

    cfg.asset_yield = 0.03;
    Real ps_with_yield = trs_par_spread(cfg, discount, credit);

    EXPECT_LT(ps_with_yield, ps_no_yield);
}

// ============================================================
// 8. 远期价格逻辑验证
// ============================================================
TEST(TRSPriceTest, ForwardPriceGrowthMatchesRate) {
    // 无 yield, 无违约: asset leg PV ≈ (forward return) * annuity
    // forward return = exp(r*tau) - 1 ≈ r*tau
    // asset leg ≈ r * Σ tau * N * P(0,t) (无违约)
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.0, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.asset_yield = 0.0;
    cfg.funding_spread = 0.0;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;

    auto result = price_trs(cfg, discount, credit);
    // asset leg PV ≈ funding leg PV (无利差, 无违约)
    EXPECT_NEAR(result.pv_asset_leg, result.pv_funding_leg, 0.5);
}

// ============================================================
// 9. 现金流明细验证
// ============================================================
TEST(TRSPriceTest, CashflowDetailsSize) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.maturity = 3.0;
    cfg.n_payments = 12;

    auto result = price_trs(cfg, discount, credit);
    EXPECT_EQ(result.asset_leg_cashflows.size(), 12u);
    EXPECT_EQ(result.funding_leg_cashflows.size(), 12u);
    EXPECT_EQ(result.payment_times.size(), 12u);
}

TEST(TRSPriceTest, FundingCashflowsPositive) {
    // receiver 视角: funding leg 现金流应为正 (支出)
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.maturity = 3.0;
    cfg.n_payments = 12;
    cfg.funding_spread = 0.01;

    auto result = price_trs(cfg, discount, credit);
    for (Real cf : result.funding_leg_cashflows) {
        EXPECT_GT(cf, 0.0);  // receiver 支出 = 正
    }
}

// ============================================================
// 10. 自定义支付时间表
// ============================================================
TEST(TRSPriceTest, CustomPaymentSchedule) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.maturity = 2.0;
    cfg.n_payments = 4;
    cfg.payment_times = {0.5, 1.0, 1.5, 2.0};

    auto result = price_trs(cfg, discount, credit);
    EXPECT_EQ(result.payment_times.size(), 4u);
    EXPECT_NEAR(result.payment_times[0], 0.5, 1e-10);
    EXPECT_NEAR(result.payment_times[3], 2.0, 1e-10);
}

// ============================================================
// 11. 长期 TRS 信用损失递增
// ============================================================
TEST(TRSPriceTest, LongerMaturityGivesMoreCreditLoss) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.n_payments = 20;
    cfg.funding_spread = 0.01;

    Real prev_loss = 0.0;
    for (Real T : {1.0, 2.0, 5.0, 10.0}) {
        cfg.maturity = T;
        cfg.n_payments = static_cast<Size>(T * 4);  // 季度
        cfg.payment_times.clear();
        auto result = price_trs(cfg, discount, credit);
        // 长期 → 更多违约概率 → 信用损失更负
        EXPECT_LT(result.pv_credit_loss, prev_loss - 1e-10);
        prev_loss = result.pv_credit_loss;
    }
}

// ============================================================
// 12. 边界情况
// ============================================================
TEST(TRSPriceTest, SinglePayment) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.maturity = 1.0;
    cfg.n_payments = 1;

    auto result = price_trs(cfg, discount, credit);
    EXPECT_EQ(result.payment_times.size(), 1u);
    // 单期 TRS 仍应有合理 PV
    EXPECT_TRUE(std::isfinite(result.pv));
}

TEST(TRSPriceTest, ZeroHazardZeroRecoveryNoCreditLoss) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.0, 30.0, 100);

    TRSConfig cfg;
    cfg.notional = 100.0;
    cfg.maturity = 5.0;
    cfg.n_payments = 20;

    auto result = price_trs(cfg, discount, credit);
    EXPECT_NEAR(result.pv_credit_loss, 0.0, 1e-10);
}
