// =============================================================================
// test_aj_jump_test.cpp
// Phase 5 v1.4.3 - HFE AJ Jump Test (Aït-Sahalia & Jacod 2009)
//
// 对标 R highfrequency 1.0.3:
//   AJjumpTest(pData, p=4, k=2, alignBy, alignPeriod, alphaMultiplier=4, alpha=0.975)
//   — jumpTests.R L106 + internalJumpTests.R
//
// 排幻觉点 (R 源码实测 2026-08-03):
//   D4: alpha 动态 = alphaMultiplier * sqrt(RV), 覆盖入参 alpha (R bug: critical value 用动态 alpha)
//   D5: seq(1,N,h) 整数步长抽样, h = alignPeriod * scale(alignBy)
//   D6: rse = abs(makeReturns(pData[selection])), 过滤价格后重算收益 (非 r[selection])
//   D7: Ap = (1/N)^(1-p/2)/mup * sum(rse^p)
//   D8: calculateNpk 含 fmupk(p,k) 查表项
//   D9: fmupk 硬编码表 (p=2,3,4 × k=2,3,4), 其他 (p,k) 走 MC
//
// 容差: 1e-10 (R 对标)
//
// SOURCE:
//   [AS-J 2009] Aït-Sahalia & Jacod, Annals of Statistics 37(1), 184-222
// R 源码: tests/fixtures/hfe/hf_src/highfrequency/R/jumpTests.R L106
//         tests/fixtures/hfe/hf_src/highfrequency/R/internalJumpTests.R
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/tests/aj_jump_test.hpp"
#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>

using namespace cpphub::v1::hfecon;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
constexpr Real TOL = 1e-10;
} // namespace

// =============================================================================
// TEST 1: fmupk 查表 — 9 个硬编码值 (D9)
// R: internalJumpTests.R L17-48
// =============================================================================
TEST(AJJumpTest, FmupkLookupTable) {
    // p=2
    EXPECT_NEAR(detail::fmupk(2, 2), 4.00, TOL);
    EXPECT_NEAR(detail::fmupk(2, 3), 5.00, TOL);
    EXPECT_NEAR(detail::fmupk(2, 4), 6.00, TOL);
    // p=3
    EXPECT_NEAR(detail::fmupk(3, 2), 24.07, TOL);
    EXPECT_NEAR(detail::fmupk(3, 3), 33.63, TOL);
    EXPECT_NEAR(detail::fmupk(3, 4), 43.74, TOL);
    // p=4
    EXPECT_NEAR(detail::fmupk(4, 2), 204.04, TOL);
    EXPECT_NEAR(detail::fmupk(4, 3), 320.26, TOL);
    EXPECT_NEAR(detail::fmupk(4, 4), 455.67, TOL);
}

// =============================================================================
// TEST 2: mup 解析公式 — 标准正态 p 阶绝对矩
// R: mup = 2^(p/2) * gamma((p+1)/2) / gamma(0.5)
// 手算:
//   mup(2) = 2^1 * gamma(1.5)/gamma(0.5) = 2 * 0.5*sqrt(pi)/sqrt(pi) = 1.0
//   mup(4) = 2^2 * gamma(2.5)/gamma(0.5) = 4 * 0.75*sqrt(pi)/sqrt(pi) = 3.0
// =============================================================================
TEST(AJJumpTest, MupAnalyticFormula) {
    EXPECT_NEAR(detail::mup(2), 1.0, TOL);
    EXPECT_NEAR(detail::mup(4), 3.0, TOL);
    // mup(1) = 2^0.5 * gamma(1) / gamma(0.5) = sqrt(2) * 1 / sqrt(pi) = sqrt(2/pi)
    EXPECT_NEAR(detail::mup(1), std::sqrt(2.0 / std::acos(-1.0)), TOL);
}

// =============================================================================
// TEST 3: calculateNpk(4, 2) — 手算验证 (D8)
// R: internalJumpTests.R L52-55
//
// 手算:
//   mup(4) = 3.0
//   mu2p(4) = 2^4 * gamma(4.5)/gamma(0.5) = 16 * 6.5625 = 105.0
//   fmupk(4,2) = 204.04
//   npk = (1/mup^2) * (k^(p-2)*(1+k)*mu2p + k^(p-2)*(k-1)*mup^2 - 2*k^(p/2-1)*fmupk)
//       = (1/9) * (4*3*105 + 4*1*9 - 2*2*204.04)
//       = (1/9) * (1260 + 36 - 816.16)
//       = (1/9) * 479.84
//       = 53.315555...
// =============================================================================
TEST(AJJumpTest, CalculateNpkP4K2) {
    Real npk = detail::calculate_npk(4, 2);
    Real expected = (1.0 / 9.0) * (4.0 * 3.0 * 105.0 + 4.0 * 1.0 * 9.0 - 2.0 * 2.0 * 204.04);
    EXPECT_NEAR(npk, expected, TOL);
    EXPECT_NEAR(npk, 53.31555555555556, TOL);
}

// =============================================================================
// TEST 4: calculateV — 已知 rse 序列数值验证 (D7)
// R: internalJumpTests.R L60-65
//
// 构造 rse = [0.01, 0.02, 0.01] (3 个小收益), p=4, k=2, N=100
// 手算:
//   mup(4) = 3.0, mu2p(4) = 105.0
//   sum(rse^4) = 0.01^4 + 0.02^4 + 0.01^4 = 1e-8 + 16e-8 + 1e-8 = 18e-8
//   sum(rse^8) = 0.01^8 + 0.02^8 + 0.01^8 = 1e-16 + 256e-16 + 1e-16 = 258e-16
//   Ap = (1/100)^(1-2) / 3.0 * 18e-8 = 100 / 3.0 * 18e-8 = 600 * 1e-8 = 6e-6
//   A2p = (1/100)^(1-4) / 105.0 * 258e-16 = 1e6 / 105.0 * 258e-16 = 9523.81 * 258e-16
//       = 2.45714e-12
//   npk = calculateNpk(4,2) = 53.315555...
//   V = npk * A2p / (N * Ap^2) = 53.3156 * 2.45714e-12 / (100 * (6e-6)^2)
//     = 53.3156 * 2.45714e-12 / (100 * 3.6e-11)
//     = 1.3100e-10 / 3.6e-9
//     = 0.03639...
// =============================================================================
TEST(AJJumpTest, CalculateVKnownRse) {
    std::vector<Real> rse = {0.01, 0.02, 0.01};
    Real p = 4.0, k = 2.0;
    Size N = 100;

    Real V = detail::calculate_v(rse, p, k, N);

    // 手算中间量
    Real mup = 3.0;
    Real mu2p = 105.0;
    Real sum_rse_p = std::pow(0.01, 4) + std::pow(0.02, 4) + std::pow(0.01, 4);
    Real sum_rse_2p = std::pow(0.01, 8) + std::pow(0.02, 8) + std::pow(0.01, 8);
    Real Ap = std::pow(1.0 / N, 1 - p / 2) / mup * sum_rse_p;
    Real A2p = std::pow(1.0 / N, 1 - p) / mu2p * sum_rse_2p;
    Real npk = detail::calculate_npk(4, 2);
    Real V_expected = npk * A2p / (N * Ap * Ap);

    EXPECT_NEAR(V, V_expected, TOL);
    // V 应为正数
    EXPECT_GT(V, 0.0);
}

// =============================================================================
// TEST 5: aj_jump_test 端到端 — p=4, k=2 默认参数 (D4-D6 完整流程)
// 构造无跳跃的 GBM 价格序列, 验证 ztest 有限且量级合理 (|ztest| < 10)
// =============================================================================
TEST(AJJumpTest, EndToEndP4K2NoJumps) {
    // 构造 GBM 价格序列 (无跳跃, seed=42, n=200, sigma=0.001)
    std::mt19937_64 gen(42);
    std::normal_distribution<Real> dist(0.0, 0.001);
    std::vector<Real> prices;
    prices.reserve(200);
    Real log_p = std::log(100.0);
    prices.push_back(100.0);
    for (int i = 1; i < 200; ++i) {
        log_p += dist(gen);
        prices.push_back(std::exp(log_p));
    }

    auto result = aj_jump_test(prices, 4, 2, "seconds", 1, 4.0, 0.975);

    // ztest 应有限
    EXPECT_TRUE(std::isfinite(result.ztest));
    // 无跳跃时 |ztest| 应在合理范围 (< 10)
    EXPECT_LT(std::fabs(result.ztest), 10.0);
    // pvalue 应在 [0, 1]
    EXPECT_GE(result.pvalue, 0.0);
    EXPECT_LE(result.pvalue, 1.0);
}

// =============================================================================
// TEST 6: aj_jump_test — alignBy="minutes", alignPeriod=5 (D5 整数抽样)
// 验证不同 alignBy/alignPeriod 下 ztest 有限
// =============================================================================
TEST(AJJumpTest, AlignByMinutesAlignPeriod5) {
    // 构造 1000 个价格点 (模拟 1 秒间隔)
    std::mt19937_64 gen(123);
    std::normal_distribution<Real> dist(0.0, 0.0005);
    std::vector<Real> prices;
    prices.reserve(1000);
    Real log_p = std::log(50.0);
    prices.push_back(50.0);
    for (int i = 1; i < 1000; ++i) {
        log_p += dist(gen);
        prices.push_back(std::exp(log_p));
    }

    // alignBy="minutes", alignPeriod=5 → h = 5*60 = 300
    // 注意: h=300, N=999 → 子采样仅 4 个价格 (3 收益) / 2 个价格 (1 收益)
    // 稀疏子采样使 |z| 极大, pvalue 可能下溢为 0 (浮点正常行为)
    auto result = aj_jump_test(prices, 4, 2, "minutes", 5, 4.0, 0.975);
    EXPECT_TRUE(std::isfinite(result.ztest));
    EXPECT_GE(result.pvalue, 0.0);
    EXPECT_LE(result.pvalue, 1.0);
}

// =============================================================================
// TEST 7: 异常处理 — 空输入或 n < 2 抛 invalid_argument
// =============================================================================
TEST(AJJumpTest, EmptyInputThrows) {
    std::vector<Real> empty_prices;
    EXPECT_THROW(
        aj_jump_test(empty_prices, 4, 2, "seconds", 1, 4.0, 0.975),
        std::invalid_argument);
}

// =============================================================================
// TEST 8: 异常处理 — 未知 alignBy 抛 invalid_argument
// =============================================================================
TEST(AJJumpTest, UnknownAlignByThrows) {
    std::vector<Real> prices = {100, 101, 102, 103};
    EXPECT_THROW(
        aj_jump_test(prices, 4, 2, "ticks", 1, 4.0, 0.975),
        std::invalid_argument);
}
