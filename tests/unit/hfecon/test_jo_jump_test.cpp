// =============================================================================
// test_jo_jump_test.cpp
// Phase 5 v1.4.3 - HFE JO Jump Test (Jiang & Oomen 2008)
//
// 对标 R highfrequency 1.0.3:
//   JOjumpTest(pData, power=4, alignBy, alignPeriod, alpha=0.975)
//   — jumpTests.R L446 + internalJumpTests.R + internals.cpp L207
//
// 排幻觉点 (R 源码实测 2026-08-03):
//   D10: R = simre(pData) 简单收益率, r = makeReturns(pData) 对数收益率
//        SwV = 2*sum(R - r)  (论文用同一收益率, R 源码用两种)
//   D11: SwV = 2*sum(R-r), 不是 sum(2*(R-r)) — 数值相同但语义明确
//   D12: mu1 = 2^(6/2)*gamma(3.5)/gamma(0.5) = μ₆ = 15.0 (6 阶矩, 非 power 阶)
//   D13: rollApplyProdWrapper: m = m-1 后窗口 m 个元素, 输出长度 n-m+1
//        internals.cpp L207-216 实测: out[i] = prod(x[i:i+m]), 共 m+1 个元素
//
// 容差: 1e-10 (R 对标)
//
// SOURCE:
//   [JO 2008] Jiang & Oomen, Mathematical Finance 18(3), doi:10.1111/j.1467-9965.2008.00343.x
// R 源码: tests/fixtures/hfe/hf_src/highfrequency/R/jumpTests.R L446
//         tests/fixtures/hfe/hf_src/highfrequency/R/internalJumpTests.R (simre, RBPVar)
//         tests/fixtures/hfe/hf_src/highfrequency/src/internals.cpp L207 (rollApplyProdWrapper)
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/tests/jo_jump_test.hpp"
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
// TEST 1: simre 简单收益率计算 (D10)
// R: internalJumpTests.R L79-101
//   R[0] = 0
//   R[i] = P[i]/P[i-1] - 1
//
// 手算: P = [100, 110, 105]
//   R[0] = 0
//   R[1] = 110/100 - 1 = 0.1
//   R[2] = 105/110 - 1 = -5/110 = -0.04545454545454545...
// =============================================================================
TEST(JOJumpTest, SimreSimpleReturns) {
    std::vector<Real> prices = {100.0, 110.0, 105.0};
    auto R = detail::simre(prices);

    ASSERT_EQ(R.size(), 3u);
    EXPECT_NEAR(R[0], 0.0, TOL);
    EXPECT_NEAR(R[1], 0.1, TOL);
    EXPECT_NEAR(R[2], -5.0 / 110.0, TOL);
}

// =============================================================================
// TEST 2: make_returns 对数收益率 (D10)
// R: makeReturns(pData)
//   r[0] = 0
//   r[i] = log(P[i]) - log(P[i-1])
//
// 手算: P = [100, 110, 105]
//   r[0] = 0
//   r[1] = log(110) - log(100) = log(1.1) = 0.09531017980432493
//   r[2] = log(105) - log(110) = log(105/110) = -0.04652001563082960
// =============================================================================
TEST(JOJumpTest, MakeReturnsLogReturns) {
    std::vector<Real> prices = {100.0, 110.0, 105.0};
    auto r = detail::make_returns(prices);

    ASSERT_EQ(r.size(), 3u);
    EXPECT_NEAR(r[0], 0.0, TOL);
    EXPECT_NEAR(r[1], std::log(1.1), TOL);
    EXPECT_NEAR(r[2], std::log(105.0 / 110.0), TOL);
}

// =============================================================================
// TEST 3: roll_apply_prod_wrapper 窗口验证 (D13)
// R: internals.cpp L207-216
//   m = m - 1 (关键偏移!)
//   out 长度 = n - m (即 n - (orig_m - 1) = n - orig_m + 1)
//   out[i] = prod(x[i : i+m]), 共 m+1 = orig_m 个元素
//
// 手算: r = [1, 2, 3, 4, 5], orig_m = 4
//   m = 4 - 1 = 3
//   out 长度 = 5 - 3 = 2
//   out[0] = prod(r[0:3]) = 1*2*3*4 = 24
//   out[1] = prod(r[1:4]) = 2*3*4*5 = 120
// =============================================================================
TEST(JOJumpTest, RollApplyProdWrapperWindow) {
    std::vector<Real> r = {1.0, 2.0, 3.0, 4.0, 5.0};
    auto out = detail::roll_apply_prod_wrapper(r, 4);

    ASSERT_EQ(out.size(), 2u);
    EXPECT_NEAR(out[0], 24.0, TOL);  // 1*2*3*4
    EXPECT_NEAR(out[1], 120.0, TOL); // 2*3*4*5
}

// =============================================================================
// TEST 4: rbp_var 双幂变差 (RBPVar)
// R: internalRealizedMeasures.R L256-262
//   bv = (pi/2) * sum(|r[0:n-2]| * |r[1:n-1]|)
//   注意: r[0]=0, 所以第一项 |r[0]|*|r[1]| = 0
//
// 手算: r = [0, 0.01, 0.02, -0.01, 0.03] (n=5)
//   sum = |0|*|0.01| + |0.01|*|0.02| + |0.02|*|-0.01| + |-0.01|*|0.03|
//       = 0 + 0.0002 + 0.0002 + 0.0003 = 0.0007
//   bv = (pi/2) * 0.0007
// =============================================================================
TEST(JOJumpTest, RbpVarBipowerVariance) {
    std::vector<Real> r = {0.0, 0.01, 0.02, -0.01, 0.03};
    Real bv = detail::rbp_var(r);
    Real pi = std::acos(-1.0);
    Real expected = (pi / 2.0) * 0.0007;
    EXPECT_NEAR(bv, expected, TOL);
}

// =============================================================================
// TEST 5: mu1 6 阶矩解析公式 (D12)
// R: jumpTests.R L495
//   mu1 = 2^(6/2) * gamma(1/2*(6+1)) / gamma(1/2)
//       = 2^3 * gamma(3.5) / gamma(0.5)
//       = 8 * (15/8 * sqrt(pi)) / sqrt(pi)
//       = 15.0
//
// 手算: gamma(3.5) = 2.5*1.5*0.5*sqrt(pi) = 1.875*sqrt(pi)
//       gamma(0.5) = sqrt(pi)
//       mu1 = 8 * 1.875 = 15.0
// =============================================================================
TEST(JOJumpTest, Mu1SixthMoment) {
    Real mu1 = detail::mu1_jo();
    EXPECT_NEAR(mu1, 15.0, TOL);
}

// =============================================================================
// TEST 6: jo_jump_test 端到端 power=4 默认参数 (D10-D13 完整流程)
// 构造无跳跃的 GBM 价格序列, 验证 ztest 有限且量级合理
// =============================================================================
TEST(JOJumpTest, EndToEndPower4NoJumps) {
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

    auto result = jo_jump_test(prices, 4, "seconds", 1, 0.975);

    EXPECT_TRUE(std::isfinite(result.ztest));
    EXPECT_LT(std::fabs(result.ztest), 20.0);
    EXPECT_GE(result.pvalue, 0.0);
    EXPECT_LE(result.pvalue, 1.0);
    EXPECT_LT(result.criticalLower, 0.0);
    EXPECT_GT(result.criticalUpper, 0.0);
}

// =============================================================================
// TEST 7: jo_jump_test power=6 分支 (D10-D13 完整流程, power=6)
// 构造无跳跃的 GBM 价格序列, 验证 ztest 有限且量级合理
// =============================================================================
TEST(JOJumpTest, EndToEndPower6NoJumps) {
    std::mt19937_64 gen(123);
    std::normal_distribution<Real> dist(0.0, 0.001);
    std::vector<Real> prices;
    prices.reserve(300);
    Real log_p = std::log(100.0);
    prices.push_back(100.0);
    for (int i = 1; i < 300; ++i) {
        log_p += dist(gen);
        prices.push_back(std::exp(log_p));
    }

    auto result = jo_jump_test(prices, 6, "seconds", 1, 0.975);

    EXPECT_TRUE(std::isfinite(result.ztest));
    EXPECT_LT(std::fabs(result.ztest), 20.0);
    EXPECT_GE(result.pvalue, 0.0);
    EXPECT_LE(result.pvalue, 1.0);
}

// =============================================================================
// TEST 8: 异常处理 — 空输入抛 invalid_argument
// =============================================================================
TEST(JOJumpTest, EmptyInputThrows) {
    std::vector<Real> empty_prices;
    EXPECT_THROW(
        jo_jump_test(empty_prices, 4, "seconds", 1, 0.975),
        std::invalid_argument);
}

// =============================================================================
// TEST 9: 异常处理 — N 不足 (N <= power+1 导致分母 N-power-1 <= 0)
// power=4 → 需要 N > 5, 即至少 7 个价格点
// =============================================================================
TEST(JOJumpTest, InsufficientNThrows) {
    std::vector<Real> prices = {100.0, 101.0, 102.0, 103.0, 104.0, 105.0};
    EXPECT_THROW(
        jo_jump_test(prices, 4, "seconds", 1, 0.975),
        std::invalid_argument);
}

// =============================================================================
// TEST 10: 异常处理 — 未知 power 抛 invalid_argument
// R 源码只支持 power=4 和 power=6
// =============================================================================
TEST(JOJumpTest, UnknownPowerThrows) {
    std::vector<Real> prices = {100.0, 101.0, 102.0, 103.0, 104.0,
                                105.0, 106.0, 107.0, 108.0, 109.0};
    EXPECT_THROW(
        jo_jump_test(prices, 5, "seconds", 1, 0.975),
        std::invalid_argument);
}
