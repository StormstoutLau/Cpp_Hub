// v1.2 Batch 7: 信用利差期权 (Credit Spread Option) 单元测试
// 覆盖: CSOConfig/CSOResult + price_credit_spread_option (Black 1976)
//       + credit_spread_implied_vol + credit_spread_parity_check
// 测试维度: 基本定价/平价关系/敲出效应/Greeks/隐含波动率/边界/退化情形
#include <gtest/gtest.h>
#include "cpphub/instruments/credit/credit_spread_option.hpp"
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
// 1. CSOConfig 基本验证
// ============================================================
TEST(CSOConfigTest, DefaultValues) {
    CSOConfig cfg;
    EXPECT_EQ(cfg.notional, 1.0);
    EXPECT_EQ(cfg.strike_spread, 0.01);
    EXPECT_EQ(cfg.option_maturity, 1.0);
    EXPECT_EQ(cfg.cds_tenor, 5.0);
    EXPECT_EQ(cfg.forward_spread, 0.01);
    EXPECT_EQ(cfg.spread_volatility, 0.3);
    EXPECT_TRUE(cfg.is_payer);
    EXPECT_TRUE(cfg.is_knock_out);
}

TEST(CSOConfigTest, ValidationThrowsOnInvalid) {
    CSOConfig cfg;
    cfg.notional = -1.0;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.notional = 1.0;
    cfg.strike_spread = -0.01;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.strike_spread = 0.01;
    cfg.option_maturity = -1.0;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.option_maturity = 1.0;
    cfg.forward_spread = -0.01;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    cfg.forward_spread = 0.01;
    cfg.spread_volatility = -0.5;
    EXPECT_THROW(cfg.validate(), std::invalid_argument);
}

// ============================================================
// 2. 基本定价 (ATM, ITM, OTM)
// ============================================================
TEST(CSOPriceTest, ATMCallHasPositiveValue) {
    // ATM: F = K
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.4, 30.0, 100);  // 无违约

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.01;       // 100bp
    cfg.forward_spread = 0.01;      // ATM
    cfg.option_maturity = 1.0;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.5;
    cfg.is_payer = true;
    cfg.is_knock_out = false;       // 先测试非敲出

    auto result = price_credit_spread_option(cfg, discount, credit);
    EXPECT_GT(result.pv, 0.0);

    // ATM call 解析: PV = df * F * [2N(0.5σ√T) - 1] * N * τ
    // d1 = 0.5*σ*√T = 0.25, N(0.25) 由 normal_cdf 计算
    Real d1 = 0.25;
    Real N_d1 = normal_cdf(d1);
    Real expected = std::exp(-0.05) * 0.01 * (2.0 * N_d1 - 1.0) * 100 * 5;
    EXPECT_NEAR(result.pv, expected, 1e-10);
}

TEST(CSOPriceTest, ITMCallGreaterThanATM) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.4, 30.0, 100);

    CSOConfig cfg_atm;
    cfg_atm.notional = 100.0;
    cfg_atm.strike_spread = 0.01;
    cfg_atm.forward_spread = 0.01;
    cfg_atm.option_maturity = 1.0;
    cfg_atm.cds_tenor = 5.0;
    cfg_atm.spread_volatility = 0.3;
    cfg_atm.is_knock_out = false;
    auto result_atm = price_credit_spread_option(cfg_atm, discount, credit);

    CSOConfig cfg_itm = cfg_atm;
    cfg_itm.forward_spread = 0.015;  // ITM (F > K)
    auto result_itm = price_credit_spread_option(cfg_itm, discount, credit);

    EXPECT_GT(result_itm.pv, result_atm.pv);
    EXPECT_GT(result_itm.intrinsic_value, 0.0);
    EXPECT_NEAR(result_atm.intrinsic_value, 0.0, 1e-12);  // ATM 内在价值 = 0
}

TEST(CSOPriceTest, OTMCallLessThanATM) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.4, 30.0, 100);

    CSOConfig cfg_atm;
    cfg_atm.notional = 100.0;
    cfg_atm.strike_spread = 0.01;
    cfg_atm.forward_spread = 0.01;
    cfg_atm.option_maturity = 1.0;
    cfg_atm.cds_tenor = 5.0;
    cfg_atm.spread_volatility = 0.3;
    cfg_atm.is_knock_out = false;
    auto result_atm = price_credit_spread_option(cfg_atm, discount, credit);

    CSOConfig cfg_otm = cfg_atm;
    cfg_otm.forward_spread = 0.008;  // OTM (F < K)
    auto result_otm = price_credit_spread_option(cfg_otm, discount, credit);

    EXPECT_LT(result_otm.pv, result_atm.pv);
    EXPECT_NEAR(result_otm.intrinsic_value, 0.0, 1e-12);  // OTM 内在价值 = 0
}

TEST(CSOPriceTest, PutCallMonotonicityInForward) {
    // Call PV 随 F 单调递增; Put PV 随 F 单调递减
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.01;
    cfg.option_maturity = 1.0;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.3;
    cfg.is_knock_out = false;

    Real prev_call_pv = -1.0;
    Real prev_put_pv = 1e10;
    for (Real F : {0.005, 0.008, 0.010, 0.012, 0.015, 0.020}) {
        cfg.forward_spread = F;
        cfg.is_payer = true;
        auto call = price_credit_spread_option(cfg, discount, credit);
        EXPECT_GT(call.pv, prev_call_pv);
        prev_call_pv = call.pv;

        cfg.is_payer = false;
        auto put = price_credit_spread_option(cfg, discount, credit);
        EXPECT_LT(put.pv, prev_put_pv);
        prev_put_pv = put.pv;
    }
}

// ============================================================
// 3. Put-Call Parity 验证
// ============================================================
TEST(CSOParityTest, PutCallParityHolds) {
    // PV_call - PV_put = Q * df * (F - K) * N * τ
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.012;
    cfg.forward_spread = 0.010;  // OTM call / ITM put
    cfg.option_maturity = 2.0;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.4;

    // 敲出
    cfg.is_knock_out = true;
    Real diff_ko = credit_spread_parity_check(cfg, discount, credit);
    EXPECT_NEAR(diff_ko, 0.0, 1e-10);

    // 非敲出
    cfg.is_knock_out = false;
    Real diff_nko = credit_spread_parity_check(cfg, discount, credit);
    EXPECT_NEAR(diff_nko, 0.0, 1e-10);
}

TEST(CSOParityTest, ParityHoldsAcrossStrikes) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.forward_spread = 0.01;
    cfg.option_maturity = 1.0;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.3;
    cfg.is_knock_out = true;

    for (Real K : {0.005, 0.008, 0.010, 0.012, 0.015, 0.020}) {
        cfg.strike_spread = K;
        Real diff = credit_spread_parity_check(cfg, discount, credit);
        EXPECT_NEAR(diff, 0.0, 1e-10);
    }
}

// ============================================================
// 4. 敲出效应 (Knock-out)
// ============================================================
TEST(CSOKnockOutTest, KnockOutReducesValue) {
    // 有信用风险时, 敲出期权价值 < 非敲出期权价值
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.05, 0.4, 30.0, 100);  // 5% hazard

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.01;
    cfg.forward_spread = 0.01;
    cfg.option_maturity = 3.0;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.4;

    cfg.is_knock_out = true;
    auto result_ko = price_credit_spread_option(cfg, discount, credit);

    cfg.is_knock_out = false;
    auto result_nko = price_credit_spread_option(cfg, discount, credit);

    EXPECT_LT(result_ko.pv, result_nko.pv);
    EXPECT_LT(result_ko.survival_prob, 1.0);
    EXPECT_NEAR(result_nko.survival_prob, 1.0, 1e-12);
}

TEST(CSOKnockOutTest, NoDefaultRiskGivesEqualKOandNKO) {
    // 无信用风险时, 敲出与非敲出价值相同
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.01;
    cfg.forward_spread = 0.01;
    cfg.option_maturity = 2.0;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.3;

    cfg.is_knock_out = true;
    auto result_ko = price_credit_spread_option(cfg, discount, credit);

    cfg.is_knock_out = false;
    auto result_nko = price_credit_spread_option(cfg, discount, credit);

    EXPECT_NEAR(result_ko.pv, result_nko.pv, 1e-12);
}

TEST(CSOKnockOutTest, HigherHazardReducesKOValue) {
    // 敲出期权: hazard 越高 → Q 越低 → PV 越低
    auto discount = flat_zero_curve(0.05);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.01;
    cfg.forward_spread = 0.01;
    cfg.option_maturity = 3.0;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.3;
    cfg.is_knock_out = true;

    Real prev_pv = 1e10;
    for (Real h : {0.0, 0.01, 0.02, 0.05, 0.10}) {
        CreditCurve credit = CreditCurve::flat(h, 0.4, 30.0, 100);
        auto result = price_credit_spread_option(cfg, discount, credit);
        EXPECT_LT(result.pv, prev_pv + 1e-10);  // 单调递减 (允许相等)
        EXPECT_NEAR(result.survival_prob, std::exp(-h * 3.0), 1e-6);
        prev_pv = result.pv;
    }
}

// ============================================================
// 5. Greeks 验证
// ============================================================
TEST(CSOGreeksTest, DeltaMatchesFiniteDifference) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.012;
    cfg.forward_spread = 0.010;
    cfg.option_maturity = 1.5;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.4;
    cfg.is_knock_out = true;

    // Call delta
    cfg.is_payer = true;
    auto result = price_credit_spread_option(cfg, discount, credit);

    Real eps = 1e-6;
    CSOConfig cfg_up = cfg;
    cfg_up.forward_spread += eps;
    CSOConfig cfg_dn = cfg;
    cfg_dn.forward_spread -= eps;
    Real pv_up = price_credit_spread_option(cfg_up, discount, credit).pv;
    Real pv_dn = price_credit_spread_option(cfg_dn, discount, credit).pv;
    Real fd_delta = (pv_up - pv_dn) / (2.0 * eps);

    EXPECT_NEAR(result.delta, fd_delta, 1e-4);

    // Put delta
    cfg.is_payer = false;
    result = price_credit_spread_option(cfg, discount, credit);
    cfg_up.is_payer = false;
    cfg_dn.is_payer = false;
    pv_up = price_credit_spread_option(cfg_up, discount, credit).pv;
    pv_dn = price_credit_spread_option(cfg_dn, discount, credit).pv;
    fd_delta = (pv_up - pv_dn) / (2.0 * eps);
    EXPECT_NEAR(result.delta, fd_delta, 1e-4);
}

TEST(CSOGreeksTest, GammaMatchesFiniteDifference) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.012;
    cfg.forward_spread = 0.010;
    cfg.option_maturity = 1.5;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.4;
    cfg.is_knock_out = true;
    cfg.is_payer = true;

    auto result = price_credit_spread_option(cfg, discount, credit);

    Real eps = 1e-4;
    CSOConfig cfg_up = cfg;
    cfg_up.forward_spread += eps;
    CSOConfig cfg_dn = cfg;
    cfg_dn.forward_spread -= eps;
    Real pv_up = price_credit_spread_option(cfg_up, discount, credit).pv;
    Real pv_dn = price_credit_spread_option(cfg_dn, discount, credit).pv;
    Real pv_0 = result.pv;
    Real fd_gamma = (pv_up - 2.0 * pv_0 + pv_dn) / (eps * eps);

    // gamma 量级 ~36000 (因 F=0.01 较小), 使用相对容差
    Real rel_err = std::abs(result.gamma - fd_gamma) / std::abs(result.gamma);
    EXPECT_LT(rel_err, 1e-3);
}

TEST(CSOGreeksTest, VegaMatchesFiniteDifference) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.012;
    cfg.forward_spread = 0.010;
    cfg.option_maturity = 1.5;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.4;
    cfg.is_knock_out = true;
    cfg.is_payer = true;

    auto result = price_credit_spread_option(cfg, discount, credit);

    Real eps = 1e-6;
    CSOConfig cfg_up = cfg;
    cfg_up.spread_volatility += eps;
    CSOConfig cfg_dn = cfg;
    cfg_dn.spread_volatility -= eps;
    Real pv_up = price_credit_spread_option(cfg_up, discount, credit).pv;
    Real pv_dn = price_credit_spread_option(cfg_dn, discount, credit).pv;
    Real fd_vega = (pv_up - pv_dn) / (2.0 * eps);

    EXPECT_NEAR(result.vega, fd_vega, 1e-4);
}

TEST(CSOGreeksTest, GammaAndVegaEqualForCallAndPut) {
    // 相同参数下, call 和 put 的 gamma/vega 相等
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.012;
    cfg.forward_spread = 0.010;
    cfg.option_maturity = 1.5;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.4;
    cfg.is_knock_out = true;

    cfg.is_payer = true;
    auto result_call = price_credit_spread_option(cfg, discount, credit);

    cfg.is_payer = false;
    auto result_put = price_credit_spread_option(cfg, discount, credit);

    EXPECT_NEAR(result_call.gamma, result_put.gamma, 1e-10);
    EXPECT_NEAR(result_call.vega, result_put.vega, 1e-10);
}

// ============================================================
// 6. 隐含波动率
// ============================================================
TEST(CSOImpliedVolTest, RecoversOriginalVolatility) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.012;
    cfg.forward_spread = 0.010;
    cfg.option_maturity = 1.5;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.45;
    cfg.is_knock_out = true;

    cfg.is_payer = true;
    auto result = price_credit_spread_option(cfg, discount, credit);
    Real iv = credit_spread_implied_vol(cfg, discount, credit, result.pv);
    EXPECT_NEAR(iv, 0.45, 1e-6);

    cfg.is_payer = false;
    result = price_credit_spread_option(cfg, discount, credit);
    iv = credit_spread_implied_vol(cfg, discount, credit, result.pv);
    EXPECT_NEAR(iv, 0.45, 1e-6);
}

TEST(CSOImpliedVolTest, ThrowsOnNoArbitrageViolation) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.012;
    cfg.forward_spread = 0.010;
    cfg.option_maturity = 1.5;
    cfg.cds_tenor = 5.0;
    cfg.is_knock_out = true;

    // 负价格
    EXPECT_THROW(credit_spread_implied_vol(cfg, discount, credit, -1.0), std::invalid_argument);

    // 过高价格 (超过上界)
    Real Q = credit.survival_prob(cfg.option_maturity);
    Real df = discount.discount_factor(cfg.option_maturity);
    Real upper = Q * df * cfg.forward_spread * cfg.notional * cfg.cds_tenor;
    EXPECT_THROW(credit_spread_implied_vol(cfg, discount, credit, upper * 2.0), std::runtime_error);
}

// ============================================================
// 7. 边界与退化情形
// ============================================================
TEST(CSOEdgeCaseTest, ZeroVolGivesIntrinsicValue) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.008;
    cfg.forward_spread = 0.012;  // ITM call
    cfg.option_maturity = 1.0;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.0;  // 无波动
    cfg.is_knock_out = false;

    auto result = price_credit_spread_option(cfg, discount, credit);
    Real df = std::exp(-0.05 * 1.0);
    Real expected = df * (0.012 - 0.008) * 100 * 5;
    EXPECT_NEAR(result.pv, expected, 1e-10);
}

TEST(CSOEdgeCaseTest, ZeroVolOTMCallGivesZero) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.012;
    cfg.forward_spread = 0.008;  // OTM call
    cfg.option_maturity = 1.0;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.0;
    cfg.is_knock_out = false;

    auto result = price_credit_spread_option(cfg, discount, credit);
    EXPECT_NEAR(result.pv, 0.0, 1e-12);
}

TEST(CSOEdgeCaseTest, NotionalLinearity) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.012;
    cfg.forward_spread = 0.010;
    cfg.option_maturity = 1.5;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.4;
    cfg.is_knock_out = true;

    auto result_100 = price_credit_spread_option(cfg, discount, credit);

    cfg.notional = 200.0;
    auto result_200 = price_credit_spread_option(cfg, discount, credit);

    EXPECT_NEAR(result_200.pv, 2.0 * result_100.pv, 1e-8);
    EXPECT_NEAR(result_200.delta, 2.0 * result_100.delta, 1e-8);
    EXPECT_NEAR(result_200.vega, 2.0 * result_100.vega, 1e-8);
}

TEST(CSOEdgeCaseTest, TenorLinearity) {
    // PV 对 cds_tenor 线性 (因 payoff = spread * N * τ)
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.012;
    cfg.forward_spread = 0.010;
    cfg.option_maturity = 1.5;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.4;
    cfg.is_knock_out = true;

    auto result_5y = price_credit_spread_option(cfg, discount, credit);

    cfg.cds_tenor = 10.0;
    auto result_10y = price_credit_spread_option(cfg, discount, credit);

    EXPECT_NEAR(result_10y.pv, 2.0 * result_5y.pv, 1e-8);
}

TEST(CSOEdgeCaseTest, HigherVolGivesHigherValue) {
    // 对 call 和 put, 波动率越高 → 期权价值越高
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.012;
    cfg.forward_spread = 0.010;
    cfg.option_maturity = 1.5;
    cfg.cds_tenor = 5.0;
    cfg.is_knock_out = true;

    // Call
    cfg.is_payer = true;
    Real prev_call_pv = -1.0;
    for (Real sigma : {0.1, 0.2, 0.3, 0.5, 0.8}) {
        cfg.spread_volatility = sigma;
        auto result = price_credit_spread_option(cfg, discount, credit);
        EXPECT_GT(result.pv, prev_call_pv);
        prev_call_pv = result.pv;
    }

    // Put
    cfg.is_payer = false;
    Real prev_put_pv = -1.0;
    for (Real sigma : {0.1, 0.2, 0.3, 0.5, 0.8}) {
        cfg.spread_volatility = sigma;
        auto result = price_credit_spread_option(cfg, discount, credit);
        EXPECT_GT(result.pv, prev_put_pv);
        prev_put_pv = result.pv;
    }
}

TEST(CSOEdgeCaseTest, LongerMaturityGivesHigherValue) {
    // 对 ATM 期权, 期限越长 → 价值越高 (theta 为正? 实际 ATM 期限效应通常为正)
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.4, 30.0, 100);  // 无违约, 简化

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.010;
    cfg.forward_spread = 0.010;  // ATM
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.3;
    cfg.is_knock_out = false;  // 无违约, 敲出无影响
    cfg.is_payer = true;

    Real prev_pv = -1.0;
    for (Real T : {0.25, 0.5, 1.0, 2.0, 5.0}) {
        cfg.option_maturity = T;
        auto result = price_credit_spread_option(cfg, discount, credit);
        // ATM call 随 T 单调递增 (因 PV = df * F * [2N(0.5σ√T) - 1] * N * τ)
        EXPECT_GT(result.pv, prev_pv);
        prev_pv = result.pv;
    }
}

// ============================================================
// 8. 对称性: Payer vs Receiver
// ============================================================
TEST(CSOSymmetryTest, PayerReceiverSymmetryAtParity) {
    // 当 F = K (ATM) 时, call PV = put PV
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.02, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.010;
    cfg.forward_spread = 0.010;  // ATM
    cfg.option_maturity = 1.5;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.4;
    cfg.is_knock_out = true;

    cfg.is_payer = true;
    auto result_call = price_credit_spread_option(cfg, discount, credit);

    cfg.is_payer = false;
    auto result_put = price_credit_spread_option(cfg, discount, credit);

    EXPECT_NEAR(result_call.pv, result_put.pv, 1e-10);
}

// ============================================================
// 9. 数值精度验证
// ============================================================
TEST(CSONumericsTest, HandComputedATMCallPrice) {
    // 手算验证: ATM call, 无违约, 平坦利率
    // F = K = 0.01, σ = 0.3, T = 1, N = 100, τ = 5, r = 0.05
    // d1 = 0.5 * σ * √T = 0.15
    // d2 = -0.15
    // N(d1) = N(0.15) ≈ 0.55962
    // N(d2) = N(-0.15) ≈ 0.44038
    // PV = exp(-0.05) * 0.01 * (0.55962 - 0.44038) * 100 * 5
    //    = 0.95123 * 0.01 * 0.11924 * 500
    //    = 0.56715
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.0, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.01;
    cfg.forward_spread = 0.01;
    cfg.option_maturity = 1.0;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.3;
    cfg.is_knock_out = false;
    cfg.is_payer = true;

    auto result = price_credit_spread_option(cfg, discount, credit);

    Real d1_expected = 0.5 * 0.3 * 1.0;  // 0.15
    Real d2_expected = -0.15;
    EXPECT_NEAR(result.d1, d1_expected, 1e-10);
    EXPECT_NEAR(result.d2, d2_expected, 1e-10);

    Real N_d1 = 0.5596177;  // N(0.15)
    Real N_d2 = 0.4403823;  // N(-0.15)
    Real df = std::exp(-0.05);
    Real expected = df * 0.01 * (N_d1 - N_d2) * 100 * 5;
    EXPECT_NEAR(result.pv, expected, 1e-6);
    EXPECT_GT(result.pv, 0.5);
    EXPECT_LT(result.pv, 0.6);
}

TEST(CSONumericsTest, DiscountFactorAndSurvivalStored) {
    auto discount = flat_zero_curve(0.05);
    CreditCurve credit = CreditCurve::flat(0.03, 0.4, 30.0, 100);

    CSOConfig cfg;
    cfg.notional = 100.0;
    cfg.strike_spread = 0.01;
    cfg.forward_spread = 0.01;
    cfg.option_maturity = 2.0;
    cfg.cds_tenor = 5.0;
    cfg.spread_volatility = 0.3;
    cfg.is_knock_out = true;

    auto result = price_credit_spread_option(cfg, discount, credit);
    EXPECT_NEAR(result.discount_factor, std::exp(-0.05 * 2.0), 1e-10);
    EXPECT_NEAR(result.survival_prob, std::exp(-0.03 * 2.0), 1e-6);
}
