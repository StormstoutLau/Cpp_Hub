// =============================================================================
// test_liquidity_measures.cpp
// Phase 5 v1.4.3 - HFE Liquidity Measures (23 种度量 + getTradeDirection)
//
// 对标 R highfrequency 1.0.3:
//   getLiquidityMeasures(tqData, win=300)  — liquidityMeasures.R L231
//   getTradeDirection(tqData)              — liquidityMeasures.R L346
//
// 排幻觉点 (R 源码实测 2026-08-03):
//   D1: getTradeDirection 是 tick rule + midpoint 覆盖混合, 非纯 Lee-Ready
//       - 首元素 = TRUE → 1 (buy), 不是 0
//       - rets==0 → NA + LOCF (Last Observation Carried Forward)
//       - midpoint 覆盖在 LOCF 之后: price<mid→-1, price>mid→1, price==mid→保留
//   D2: realizedSpread 用 lead shift mid[t+win], 越界为 NaN
//   D3: depthImbalanceRatio = (direction*OFRSIZ/BIDSIZ)^direction
//       D=-1 时结果为**负倒数** (-BIDSIZ/OFRSIZ), 非 +BIDSIZ/OFRSIZ
//
// 容差: 1e-10 (R 对标)
//
// SOURCE:
//   [Hasbrouck 2009] Hasbrouck, "Trading Costs and Returns for U.S. Equities"
//   [Lee-Ready 1991] Lee, Ready, "Inferring trade direction from intraday data"
// R 源码: tests/fixtures/hfe/hf_src/highfrequency/R/liquidityMeasures.R
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/liquidity/liquidity_measures.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace cpphub::v1::hfecon;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
// 容差
constexpr Real TOL = 1e-10;

// 辅助: 检查 NaN
bool is_nan(Real x) { return std::isnan(x); }

// 辅助: 比较 two vectors with NaN awareness
void expect_near_or_nan(const std::vector<Real>& actual,
                         const std::vector<Real>& expected,
                         Real tol, const std::string& name) {
    ASSERT_EQ(actual.size(), expected.size())
        << "size mismatch in " << name;
    for (Size i = 0; i < actual.size(); ++i) {
        if (is_nan(expected[i])) {
            EXPECT_TRUE(is_nan(actual[i]))
                << name << "[" << i << "] expected NaN";
        } else {
            EXPECT_NEAR(actual[i], expected[i], tol)
                << name << "[" << i << "] expected " << expected[i];
        }
    }
}
} // namespace

// =============================================================================
// TEST 1: getTradeDirection — tick rule + LOCF (price == mid, 不触发覆盖)
// R: c(TRUE, fifelse(rets>0,TRUE,fifelse(rets<0,FALSE,NA)))*2-1, then nafill locf
//
// 数据: price == mid 全部相等, midpoint 覆盖不触发 (price==mid → 保留 tick rule)
//   price = [100, 101, 100, 100, 102]
//   bid   = [99,  100, 99,  99,  101]
//   ofr   = [101, 102, 101, 101, 103]
//   mid   = [100, 101, 100, 100, 102]
//   rets  = [1, -1, 0, 2]
//   tick rule: [TRUE, TRUE, FALSE, NA, TRUE] → [1, 1, -1, NA, 1]
//   LOCF:      [1, 1, -1, -1, 1]  (NA → 前值 -1)
//   midpoint:  all equal → 保留 tick = [1, 1, -1, -1, 1]
// =============================================================================
TEST(LiquidityMeasuresTest, GetTradeDirectionTickRule) {
    std::vector<Real> price = {100, 101, 100, 100, 102};
    std::vector<Real> bid   = {99,  100, 99,  99,  101};
    std::vector<Real> ofr   = {101, 102, 101, 101, 103};

    auto dir = get_trade_direction(price, bid, ofr);
    ASSERT_EQ(dir.size(), 5u);
    EXPECT_EQ(dir[0], 1);   // 首元素 = buy (TRUE → 1)
    EXPECT_EQ(dir[1], 1);   // uptick → buy
    EXPECT_EQ(dir[2], -1);  // downtick → sell
    EXPECT_EQ(dir[3], -1);  // zero-tick → LOCF (前值 -1)
    EXPECT_EQ(dir[4], 1);   // uptick → buy
}

// =============================================================================
// TEST 2: getTradeDirection — midpoint 覆盖 (price != mid, 覆盖触发)
// R: fifelse(price<midpoints, -1, fifelse(price>midpoints, 1, buys))
//
// 数据: 构造 tick rule 与 midpoint 冲突, 验证 midpoint 优先
//   price = [100, 101, 99]
//   bid   = [100, 100, 99]
//   ofr   = [101, 101, 100]
//   mid   = [100.5, 100.5, 99.5]
//   rets  = [1, -2]
//   tick rule: [TRUE, TRUE, FALSE] → [1, 1, -1]
//   midpoint 覆盖:
//     i=0: price=100 < mid=100.5 → -1 (覆盖 tick 的 1!)
//     i=1: price=101 > mid=100.5 → 1 (与 tick 一致)
//     i=2: price=99 < mid=99.5 → -1 (与 tick 一致)
//   direction = [-1, 1, -1]
// =============================================================================
TEST(LiquidityMeasuresTest, GetTradeDirectionMidpointOverride) {
    std::vector<Real> price = {100, 101, 99};
    std::vector<Real> bid   = {100, 100, 99};
    std::vector<Real> ofr   = {101, 101, 100};

    auto dir = get_trade_direction(price, bid, ofr);
    ASSERT_EQ(dir.size(), 3u);
    EXPECT_EQ(dir[0], -1);  // price < mid → 覆盖为 sell
    EXPECT_EQ(dir[1], 1);   // price > mid → buy
    EXPECT_EQ(dir[2], -1);  // price < mid → sell
}

// =============================================================================
// TEST 3: getLiquidityMeasures — 23 种度量手算验证 (n=3, win=1)
// R 源码: liquidityMeasures.R L268-312
//
// 数据:
//   price=[100,101,99], bid=[99,100,98], ofr=[101,102,100]
//   size=[100,200,150], ofrsiz=[500,400,600], bidsiz=[300,500,400]
//   mid=[100,101,99]
//   rets=[1,-2], tick=[1,1,-1], midpoint all equal → dir=[1,1,-1]
// =============================================================================
TEST(LiquidityMeasuresTest, All23Measures) {
    std::vector<Real> price   = {100, 101, 99};
    std::vector<Real> bid     = {99,  100, 98};
    std::vector<Real> ofr     = {101, 102, 100};
    std::vector<Real> size    = {100, 200, 150};
    std::vector<Real> ofrsiz  = {500, 400, 600};
    std::vector<Real> bidsiz  = {300, 500, 400};
    int win = 1;

    auto m = get_liquidity_measures(price, bid, ofr, size, ofrsiz, bidsiz,
                                     std::nullopt, win);

    // --- 1. effectiveSpread = 2*D*(P-mid) ---
    // mid=[100,101,99], D=[1,1,-1], P-mid=[0,0,0] → [0,0,0]
    expect_near_or_nan(m.effectiveSpread, {0.0, 0.0, 0.0}, TOL, "effectiveSpread");

    // --- 2. realizedSpread = 2*D*(P - mid[t+win]), win=1 ---
    // t=0: 2*1*(100-101)=-2, t=1: 2*1*(101-99)=4, t=2: 越界 NaN
    expect_near_or_nan(m.realizedSpread, {-2.0, 4.0, std::numeric_limits<Real>::quiet_NaN()},
                        TOL, "realizedSpread");

    // --- 3. valueTrade = SIZE*P ---
    expect_near_or_nan(m.valueTrade, {10000.0, 20200.0, 14850.0}, TOL, "valueTrade");

    // --- 4. signedValueTrade = D*valueTrade ---
    expect_near_or_nan(m.signedValueTrade, {10000.0, 20200.0, -14850.0}, TOL,
                        "signedValueTrade");

    // --- 5. depthImbalanceDifference = D*(OFRSIZ-BIDSIZ)/(OFRSIZ+BIDSIZ) ---
    // t=0: 1*(500-300)/(500+300)=0.25
    // t=1: 1*(400-500)/(400+500)=-1/9
    // t=2: -1*(600-400)/(600+400)=-0.2
    expect_near_or_nan(m.depthImbalanceDifference,
                        {0.25, -1.0/9.0, -0.2}, TOL, "depthImbalanceDifference");

    // --- 6. depthImbalanceRatio = (D*OFRSIZ/BIDSIZ)^D ---
    // t=0: (1*500/300)^1 = 5/3
    // t=1: (1*400/500)^1 = 0.8
    // t=2: (-1*600/400)^(-1) = (-1.5)^(-1) = -2/3  (D3: 负倒数!)
    expect_near_or_nan(m.depthImbalanceRatio,
                        {5.0/3.0, 0.8, -2.0/3.0}, TOL, "depthImbalanceRatio");

    // --- 7. proportionalEffectiveSpread = effectiveSpread/mid ---
    expect_near_or_nan(m.proportionalEffectiveSpread, {0.0, 0.0, 0.0}, TOL,
                        "proportionalEffectiveSpread");

    // --- 8. proportionalRealizedSpread = realizedSpread/mid ---
    // t=0: -2/100=-0.02, t=1: 4/101, t=2: NaN
    expect_near_or_nan(m.proportionalRealizedSpread,
                        {-0.02, 4.0/101.0, std::numeric_limits<Real>::quiet_NaN()},
                        TOL, "proportionalRealizedSpread");

    // --- 9. priceImpact = (effectiveSpread - realizedSpread)/2 ---
    // t=0: (0-(-2))/2=1, t=1: (0-4)/2=-2, t=2: NaN
    expect_near_or_nan(m.priceImpact,
                        {1.0, -2.0, std::numeric_limits<Real>::quiet_NaN()},
                        TOL, "priceImpact");

    // --- 10. proportionalPriceImpact = priceImpact/mid ---
    // t=0: 1/100=0.01, t=1: -2/101, t=2: NaN
    expect_near_or_nan(m.proportionalPriceImpact,
                        {0.01, -2.0/101.0, std::numeric_limits<Real>::quiet_NaN()},
                        TOL, "proportionalPriceImpact");

    // --- 11. halfTradedSpread = D*(P-mid) ---
    expect_near_or_nan(m.halfTradedSpread, {0.0, 0.0, 0.0}, TOL, "halfTradedSpread");

    // --- 12. proportionalHalfTradedSpread = halfTradedSpread/mid ---
    expect_near_or_nan(m.proportionalHalfTradedSpread, {0.0, 0.0, 0.0}, TOL,
                        "proportionalHalfTradedSpread");

    // --- 13. squaredLogReturn = (log(P)-log(P[t-1]))^2 ---
    // t=0: NaN, t=1: (log(101)-log(100))^2, t=2: (log(99)-log(101))^2
    Real lr1 = std::log(101.0) - std::log(100.0);
    Real lr2 = std::log(99.0) - std::log(101.0);
    expect_near_or_nan(m.squaredLogReturn,
                        {std::numeric_limits<Real>::quiet_NaN(), lr1*lr1, lr2*lr2},
                        TOL, "squaredLogReturn");

    // --- 14. absLogReturn = |log(P)-log(P[t-1])| ---
    expect_near_or_nan(m.absLogReturn,
                        {std::numeric_limits<Real>::quiet_NaN(), std::fabs(lr1), std::fabs(lr2)},
                        TOL, "absLogReturn");

    // --- 15. quotedSpread = OFR-BID ---
    expect_near_or_nan(m.quotedSpread, {2.0, 2.0, 2.0}, TOL, "quotedSpread");

    // --- 16. proportionalQuotedSpread = quotedSpread/mid ---
    expect_near_or_nan(m.proportionalQuotedSpread,
                        {2.0/100.0, 2.0/101.0, 2.0/99.0}, TOL,
                        "proportionalQuotedSpread");

    // --- 17. logQuotedSpread = log(OFR/BID) ---
    expect_near_or_nan(m.logQuotedSpread,
                        {std::log(101.0/99.0), std::log(102.0/100.0), std::log(100.0/98.0)},
                        TOL, "logQuotedSpread");

    // --- 18. logQuotedSize = log(OFRSIZ)+log(BIDSIZ) ---
    expect_near_or_nan(m.logQuotedSize,
                        {std::log(500.0)+std::log(300.0),
                         std::log(400.0)+std::log(500.0),
                         std::log(600.0)+std::log(400.0)}, TOL, "logQuotedSize");

    // --- 19. quotedSlope = quotedSpread/logQuotedSize ---
    Real qs0 = 2.0 / (std::log(500.0)+std::log(300.0));
    Real qs1 = 2.0 / (std::log(400.0)+std::log(500.0));
    Real qs2 = 2.0 / (std::log(600.0)+std::log(400.0));
    expect_near_or_nan(m.quotedSlope, {qs0, qs1, qs2}, TOL, "quotedSlope");

    // --- 20. logQSlope = logQuotedSpread/logQuotedSize ---
    Real lqs0 = std::log(101.0/99.0) / (std::log(500.0)+std::log(300.0));
    Real lqs1 = std::log(102.0/100.0) / (std::log(400.0)+std::log(500.0));
    Real lqs2 = std::log(100.0/98.0) / (std::log(600.0)+std::log(400.0));
    expect_near_or_nan(m.logQSlope, {lqs0, lqs1, lqs2}, TOL, "logQSlope");

    // --- 21. midQuoteSquaredReturn = (log(mid)-log(mid[t-1]))^2 ---
    // mid=[100,101,99], 同 squaredLogReturn
    expect_near_or_nan(m.midQuoteSquaredReturn,
                        {std::numeric_limits<Real>::quiet_NaN(), lr1*lr1, lr2*lr2},
                        TOL, "midQuoteSquaredReturn");

    // --- 22. midQuoteAbsReturn = |log(mid)-log(mid[t-1])| ---
    expect_near_or_nan(m.midQuoteAbsReturn,
                        {std::numeric_limits<Real>::quiet_NaN(), std::fabs(lr1), std::fabs(lr2)},
                        TOL, "midQuoteAbsReturn");

    // --- 23. signedTradeSize = D*SIZE ---
    expect_near_or_nan(m.signedTradeSize, {100.0, 200.0, -150.0}, TOL, "signedTradeSize");
}

// =============================================================================
// TEST 4: 用户 DIRECTION 输入 — 跳过 getTradeDirection, 直接使用
// R: if('DIRECTION' %in% colnames(tqData)) direction := DIRECTION
// =============================================================================
TEST(LiquidityMeasuresTest, UserDirectionInput) {
    std::vector<Real> price  = {100, 101, 102};
    std::vector<Real> bid    = {99,  100, 101};
    std::vector<Real> ofr    = {101, 102, 103};
    std::vector<Real> size   = {100, 200, 150};
    std::vector<Real> ofrsiz = {500, 400, 600};
    std::vector<Real> bidsiz = {300, 500, 400};
    std::vector<int> userDir = {-1, 1, -1};  // 用户指定方向

    auto m = get_liquidity_measures(price, bid, ofr, size, ofrsiz, bidsiz,
                                     userDir, 1);

    // effectiveSpread = 2*D*(P-mid), mid=[100,101,102]
    // t=0: 2*(-1)*(100-100)=0, t=1: 2*1*(101-101)=0, t=2: 2*(-1)*(102-102)=0
    expect_near_or_nan(m.effectiveSpread, {0.0, 0.0, 0.0}, TOL, "effectiveSpread");

    // signedTradeSize = D*SIZE = [-100, 200, -150]
    expect_near_or_nan(m.signedTradeSize, {-100.0, 200.0, -150.0}, TOL,
                        "signedTradeSize");

    // signedValueTrade = D*SIZE*P = [-1*100*100, 1*200*101, -1*150*102]
    expect_near_or_nan(m.signedValueTrade, {-10000.0, 20200.0, -15300.0}, TOL,
                        "signedValueTrade");
}

// =============================================================================
// TEST 5: realizedSpread 越界 NaN — win >= n 时全部越界
// R: shift(midpoints, win, type="lead") 越界 = NA → NaN
// =============================================================================
TEST(LiquidityMeasuresTest, RealizedSpreadOutOfBoundsNaN) {
    std::vector<Real> price  = {100, 101, 102};
    std::vector<Real> bid    = {99,  100, 101};
    std::vector<Real> ofr    = {101, 102, 103};
    std::vector<Real> size   = {100, 200, 150};
    std::vector<Real> ofrsiz = {500, 400, 600};
    std::vector<Real> bidsiz = {300, 500, 400};

    // win=3 (n=3), 全部越界 → realizedSpread 全 NaN
    auto m = get_liquidity_measures(price, bid, ofr, size, ofrsiz, bidsiz,
                                     std::nullopt, 3);
    for (Size i = 0; i < 3; ++i) {
        EXPECT_TRUE(is_nan(m.realizedSpread[i]))
            << "realizedSpread[" << i << "] should be NaN (win=3 >= n=3)";
    }
}

// =============================================================================
// TEST 6: depthImbalanceRatio D=-1 负倒数 (D3 排幻觉)
// R: (direction * OFRSIZ / BIDSIZ) ^ direction
// D=-1: (-OFRSIZ/BIDSIZ)^(-1) = -BIDSIZ/OFRSIZ (负倒数, 非 +BIDSIZ/OFRSIZ)
// =============================================================================
TEST(LiquidityMeasuresTest, DepthImbalanceRatioNegativeReciprocal) {
    // 构造 D=-1 的场景: price < mid → 强制 sell
    std::vector<Real> price  = {98};       // price < mid → D=-1
    std::vector<Real> bid    = {100};
    std::vector<Real> ofr    = {102};      // mid=101
    std::vector<Real> size   = {100};
    std::vector<Real> ofrsiz = {200};
    std::vector<Real> bidsiz = {100};
    // D=-1: depthImbalanceRatio = (-1*200/100)^(-1) = (-2)^(-1) = -0.5

    auto m = get_liquidity_measures(price, bid, ofr, size, ofrsiz, bidsiz,
                                     std::nullopt, 1);
    ASSERT_EQ(m.depthImbalanceRatio.size(), 1u);
    // 负倒数 -0.5, 非 +0.5
    EXPECT_NEAR(m.depthImbalanceRatio[0], -0.5, TOL);
}

// =============================================================================
// TEST 7: 异常处理 — 长度不一致抛 invalid_argument
// =============================================================================
TEST(LiquidityMeasuresTest, LengthMismatchThrows) {
    std::vector<Real> price  = {100, 101};
    std::vector<Real> bid    = {99};  // 长度不一致
    std::vector<Real> ofr    = {101};
    std::vector<Real> size   = {100};
    std::vector<Real> ofrsiz = {500};
    std::vector<Real> bidsiz = {300};

    EXPECT_THROW(
        get_liquidity_measures(price, bid, ofr, size, ofrsiz, bidsiz,
                                std::nullopt, 1),
        std::invalid_argument);
}

// =============================================================================
// TEST 8: 用户 DIRECTION 含非法值 (非 ±1) 抛 invalid_argument
// R: if(any(!(tqData$DIRECTION %in% c(-1,1)))) stop(...)
// =============================================================================
TEST(LiquidityMeasuresTest, InvalidUserDirectionThrows) {
    std::vector<Real> price  = {100, 101};
    std::vector<Real> bid    = {99,  100};
    std::vector<Real> ofr    = {101, 102};
    std::vector<Real> size   = {100, 200};
    std::vector<Real> ofrsiz = {500, 400};
    std::vector<Real> bidsiz = {300, 500};
    std::vector<int> badDir  = {1, 0};  // 0 非法

    EXPECT_THROW(
        get_liquidity_measures(price, bid, ofr, size, ofrsiz, bidsiz,
                                badDir, 1),
        std::invalid_argument);
}
