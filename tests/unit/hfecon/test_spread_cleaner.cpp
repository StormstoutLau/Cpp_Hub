// =============================================================================
// test_spread_cleaner.cpp
// Phase 5 v1.4.3 - HFE Spread Cleaner (rmLargeSpread / rmNegativeSpread / spreadPrices)
//
// 对标 R highfrequency 1.0.3:
//   rmLargeSpread(qData, maxi = 50)  — dataHandling.R L1617
//   rmNegativeSpread(qData)          — dataHandling.R L1670
//   spreadPrices(data)               — dataHandling.R L3019
//
// 容差: 1e-12 (无浮点累积, 纯筛选/转换)
//
// SOURCE:
//   [BKS 2022] Boudt, Kleen, Sjoerup, JSS 104(8), 1-36, doi:10.18637/jss.v104.i08
// R 源码实测: 2026-08-03
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/liquidity/spread_cleaner.hpp"
#include <vector>
#include <cmath>

using namespace cpphub::v1::hfecon;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
// 构造报价序列辅助函数
std::vector<QuoteObservation> make_quotes(
    const std::vector<Real>& bids,
    const std::vector<Real>& asks,
    const std::vector<int64_t>& ts) {
    std::vector<QuoteObservation> q;
    for (Size i = 0; i < bids.size(); ++i) {
        q.push_back({ts[i], bids[i], asks[i]});
    }
    return q;
}
} // namespace

// =============================================================================
// TEST 1: rm_negative_spread — 保留 OFR > BID (严格大于)
// R: qData[OFR > BID]
// =============================================================================
TEST(SpreadCleanerTest, RmNegativeSpreadBasic) {
    auto quotes = make_quotes(
        {100.0, 100.0, 100.0, 100.0},  // BID
        {100.1, 99.9,  100.0, 100.2},  // OFR (第2条 OFR<BID, 第3条 OFR==BID)
        {1, 2, 3, 4});
    auto result = rm_negative_spread(quotes);
    // 期望: 保留第1条 (100.1>100) 和第4条 (100.2>100), 删除第2/3条
    ASSERT_EQ(result.size(), 2u);
    EXPECT_DOUBLE_EQ(result[0].bid, 100.0);
    EXPECT_DOUBLE_EQ(result[0].ask, 100.1);
    EXPECT_DOUBLE_EQ(result[1].bid, 100.0);
    EXPECT_DOUBLE_EQ(result[1].ask, 100.2);
}

// =============================================================================
// TEST 2: rm_large_spread — 每日 SPREAD 中位数筛选, 保留 SPREAD < SPREAD_MEDIAN * maxi
// R: qData[, SPREAD := OFR - BID][, SPREAD_MEDIAN := median(SPREAD), by = "DATE"][SPREAD < (SPREAD_MEDIAN * maxi)]
// =============================================================================
TEST(SpreadCleanerTest, RmLargeSpreadBasic) {
    // 构造单日数据: 大部分 spread=0.1, 一条异常 spread=10.0
    // median(0.1, 0.1, 0.1, 0.1, 10.0) = 0.1
    // maxi=50 -> 阈值 = 0.1 * 50 = 5.0
    // 保留 spread < 5.0 的 4 条, 删除 spread=10.0 的 1 条
    std::vector<QuoteObservation> quotes;
    for (int i = 0; i < 4; ++i) {
        quotes.push_back({static_cast<int64_t>(i), 100.0, 100.1});  // spread=0.1
    }
    quotes.push_back({4, 100.0, 110.0});  // spread=10.0 (异常)

    auto result = rm_large_spread(quotes, 50.0);
    ASSERT_EQ(result.size(), 4u);
    for (const auto& q : result) {
        EXPECT_NEAR(q.ask - q.bid, 0.1, 1e-12);
    }
}

// =============================================================================
// TEST 3: rm_large_spread — 多日分组, 每日独立计算中位数
// R: by = "DATE" 分组
// =============================================================================
TEST(SpreadCleanerTest, RmLargeSpreadMultiDay) {
    // Day 1 (ts 0-3): spreads = {0.1, 0.1, 10.0, 0.1}, median=0.1, 阈值=5.0, 保留 3 条
    // Day 2 (ts 86400*1e9 + 0-3): spreads = {1.0, 1.0, 1.0}, median=1.0, 阈值=50.0, 保留 3 条
    const int64_t one_day_ns = 86400LL * 1000000000LL;
    std::vector<QuoteObservation> quotes = {
        {0,             100.0, 100.1},   // day1 spread=0.1
        {1,             100.0, 100.1},   // day1 spread=0.1
        {2,             100.0, 110.0},   // day1 spread=10.0 (异常)
        {3,             100.0, 100.1},   // day1 spread=0.1
        {one_day_ns + 0, 200.0, 201.0},  // day2 spread=1.0
        {one_day_ns + 1, 200.0, 201.0},  // day2 spread=1.0
        {one_day_ns + 2, 200.0, 201.0},  // day2 spread=1.0
    };
    auto result = rm_large_spread(quotes, 50.0);
    // Day1: 保留 3 条 (删除 spread=10.0), Day2: 保留 3 条
    ASSERT_EQ(result.size(), 6u);
    // 验证 day1 保留的 spread 都是 0.1
    for (Size i = 0; i < 3; ++i) {
        EXPECT_NEAR(result[i].ask - result[i].bid, 0.1, 1e-12);
    }
    // 验证 day2 保留的 spread 都是 1.0
    for (Size i = 3; i < 6; ++i) {
        EXPECT_NEAR(result[i].ask - result[i].bid, 1.0, 1e-12);
    }
}

// =============================================================================
// TEST 4: spread_prices — 长格式 (DT, SYMBOL, PRICE) → 宽格式
// R: split by SYMBOL -> outer join on DT
// =============================================================================
TEST(SpreadCleanerTest, SpreadPricesLongToWide) {
    // 构造长格式: 2 个 SYMBOL, 部分时间点不重叠
    std::vector<int64_t> dt = {
        0, 1, 2, 3,        // ETF: 全部 4 个时间点
        1, 2, 3            // AAA: 时间点 1,2,3 (时间点 0 缺失)
    };
    std::vector<std::string> symbols = {
        "ETF", "ETF", "ETF", "ETF",
        "AAA", "AAA", "AAA"
    };
    std::vector<Real> prices = {
        100.0, 101.0, 102.0, 103.0,  // ETF
        50.0, 51.0, 52.0             // AAA
    };

    auto result = spread_prices(dt, symbols, prices);
    // 期望: 宽格式行数 = 4 (dt=0,1,2,3), 列数 = 3 (DT + ETF + AAA)
    // dt=0: ETF=100.0, AAA=NaN
    // dt=1: ETF=101.0, AAA=50.0
    // dt=2: ETF=102.0, AAA=51.0
    // dt=3: ETF=103.0, AAA=52.0
    ASSERT_EQ(result.size(), 4u);
    EXPECT_EQ(result[0].dt, 0);
    EXPECT_NEAR(result[0].prices[0], 100.0, 1e-12);
    EXPECT_TRUE(std::isnan(result[0].prices[1]));

    EXPECT_EQ(result[1].dt, 1);
    EXPECT_NEAR(result[1].prices[0], 101.0, 1e-12);
    EXPECT_NEAR(result[1].prices[1], 50.0, 1e-12);

    EXPECT_EQ(result[3].dt, 3);
    EXPECT_NEAR(result[3].prices[0], 103.0, 1e-12);
    EXPECT_NEAR(result[3].prices[1], 52.0, 1e-12);
}
