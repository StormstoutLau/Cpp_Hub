// =============================================================================
// test_amihud.cpp
// Phase 5 v1.4.3 - HFE Amihud Illiquidity (Amihud 2002, JFM 52(5))
//
// 公式: ILLIQ = (1/N) * sum(|r_t| / DVOL_t)
//   r_t     = 日对数收益率
//   DVOL_t  = 日美元成交额 (dollar volume, 单位: 货币)
//   N       = 观测数
//
// R 对照: highfrequency 无直接实现; 公式来自 Amihud (2002)
// 容差: 1e-12 (纯算术, 无浮点累积)
//
// SOURCE:
//   [Amihud 2002] Yakov Amihud, J. Financial Markets 52(5), 2002
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/liquidity/amihud.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace cpphub::v1::hfecon;
using cpphub::v1::Real;
using cpphub::v1::Size;

// =============================================================================
// TEST 1: 基本计算 — 手算验证
// 输入:
//   dailyReturns      = [0.01, -0.02, 0.03]  (3 天对数收益率)
//   dailyDollarVolume = [1e7,  2e7,  3e7]    (3 天美元成交额)
// 手算:
//   |r_1|/DVOL_1 = 0.01 / 1e7 = 1e-9
//   |r_2|/DVOL_2 = 0.02 / 2e7 = 1e-9
//   |r_3|/DVOL_3 = 0.03 / 3e7 = 1e-9
//   sum          = 3e-9
//   ILLIQ        = (1/3) * 3e-9 = 1e-9
// =============================================================================
TEST(AmihudTest, BasicCalculation) {
    std::vector<Real> dailyReturns = {0.01, -0.02, 0.03};
    std::vector<Real> dailyDollarVolume = {1e7, 2e7, 3e7};

    Real illiq = amihud_illiquidity(dailyReturns, dailyDollarVolume);
    EXPECT_NEAR(illiq, 1e-9, 1e-12);
}

// =============================================================================
// TEST 2: 零成交额异常 — 任一 DVOL_t = 0 时应抛 invalid_argument
// Amihud 公式含除法, 零成交额会导致除零. R 实现中通常用 na.rm 但 C++ 显式抛异常更安全.
// =============================================================================
TEST(AmihudTest, ZeroDollarVolumeThrows) {
    std::vector<Real> dailyReturns = {0.01, -0.02, 0.03};
    std::vector<Real> dailyDollarVolume = {1e7, 0.0, 3e7};  // 第 2 天成交额为 0

    EXPECT_THROW(
        amihud_illiquidity(dailyReturns, dailyDollarVolume),
        std::invalid_argument);
}

// =============================================================================
// TEST 3: 长度不一致异常 — returns 和 dollarVolume 长度不同时抛异常
// =============================================================================
TEST(AmihudTest, LengthMismatchThrows) {
    std::vector<Real> dailyReturns = {0.01, -0.02, 0.03};
    std::vector<Real> dailyDollarVolume = {1e7, 2e7};  // 长度不一致

    EXPECT_THROW(
        amihud_illiquidity(dailyReturns, dailyDollarVolume),
        std::invalid_argument);
}

// =============================================================================
// TEST 4: 空输入异常 — 空序列时抛 invalid_argument
// =============================================================================
TEST(AmihudTest, EmptyInputThrows) {
    std::vector<Real> dailyReturns;
    std::vector<Real> dailyDollarVolume;

    EXPECT_THROW(
        amihud_illiquidity(dailyReturns, dailyDollarVolume),
        std::invalid_argument);
}

// =============================================================================
// TEST 5: 单日输入 — N=1, 验证无累积错误
// =============================================================================
TEST(AmihudTest, SingleDayInput) {
    std::vector<Real> dailyReturns = {0.05};
    std::vector<Real> dailyDollarVolume = {5e6};

    Real illiq = amihud_illiquidity(dailyReturns, dailyDollarVolume);
    // (1/1) * |0.05| / 5e6 = 1e-8
    EXPECT_NEAR(illiq, 1e-8, 1e-12);
}
