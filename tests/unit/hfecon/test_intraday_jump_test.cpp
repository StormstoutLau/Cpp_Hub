// =============================================================================
// test_intraday_jump_test.cpp
// Phase 5 v1.4.3 - HFE Intraday Jump Test (Lee & Mykland 2008)
//
// 对标 R highfrequency 1.0.3:
//   intradayJumpTest(pData, volEstimator="RM", driftEstimator="none", alpha=0.95, ...)
//   — jumpTests.R L583 (RM 模式简化版, PARM 推迟 v1.4.4)
//
// 排幻觉点 (R 源码实测 2026-08-03):
//   D14: vol 调整 = sqrt(vol^2 / (lookBackPeriod-2))  (Lee-Mykland 原文无此调整)
//   D15: Cn 无 sqrt(2/pi) 常数  (R 去除常数使 L~N(0,1), Lee-Mykland Eq.12 有常数)
//   D16: n = NROW(pData) 原始观测数  (非对齐后观测数)
//
// 简化假设 (v1.4.3):
//   - 输入为单日等间隔价格 vector (无需 aggregatePrice/按日分组)
//   - 仅支持 rBPCov RM 估计器 (rMinRVar/rMedRVar 推迟 v1.4.4)
//   - drift = 0 (driftEstimator="none")
//
// 容差: 1e-8 (R 对标, 滚动窗口浮点累积)
//
// SOURCE:
//   [LM 2008] Lee & Mykland, JFE 6(5)
//   [COP 2014] Christensen, Oomen, Podolskij, JFE 144, 576-599
// R 源码: tests/fixtures/hfe/hf_src/highfrequency/R/jumpTests.R L583
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/tests/intraday_jump_test.hpp"
#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>

using namespace cpphub::v1::hfecon;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
constexpr Real TOL = 1e-8;
} // namespace

// =============================================================================
// TEST 1: rolling_rbp_var 滚动窗口 rBPCov (D14 前置)
// R: spotVol RM 模式 = 滚动窗口 rBPCov
//   rbp_var(r) = (pi/2) * sum(|r[0:n-2]| * |r[1:n-1]|)
//   滚动: vol[j] = rbp_var(r[j-lookBack+1 : j])
//
// 手算: r = [0.01, 0.02, -0.01, 0.03, 0.02], lookBack = 3
//   窗口 0 [0.01, 0.02, -0.01]:
//     sum = |0.01|*|0.02| + |0.02|*|-0.01| = 0.0002 + 0.0002 = 0.0004
//     vol[0] = (pi/2) * 0.0004
//   窗口 1 [0.02, -0.01, 0.03]:
//     sum = |0.02|*|-0.01| + |-0.01|*|0.03| = 0.0002 + 0.0003 = 0.0005
//     vol[1] = (pi/2) * 0.0005
//   窗口 2 [-0.01, 0.03, 0.02]:
//     sum = |-0.01|*|0.03| + |0.03|*|0.02| = 0.0003 + 0.0006 = 0.0009
//     vol[2] = (pi/2) * 0.0009
//   输出长度 = 5 - 3 + 1 = 3
// =============================================================================
TEST(IntradayJumpTest, RollingRbpVar) {
    std::vector<Real> r = {0.01, 0.02, -0.01, 0.03, 0.02};
    auto vol = detail::rolling_rbp_var(r, 3);

    ASSERT_EQ(vol.size(), 3u);
    Real pi = std::acos(-1.0);
    EXPECT_NEAR(vol[0], (pi / 2.0) * 0.0004, TOL);
    EXPECT_NEAR(vol[1], (pi / 2.0) * 0.0005, TOL);
    EXPECT_NEAR(vol[2], (pi / 2.0) * 0.0009, TOL);
}

// =============================================================================
// TEST 2: lee_mykland_critical_value 临界值 (D15, D16)
// R: jumpTests.R L707-722
//   n = NROW(pData)  (D16: 原始观测数)
//   Cn = sqrt(2*log(n)) - (log(pi)+log(log(n)))/(2*sqrt(2*log(n)))  (D15: 无 sqrt(2/pi))
//   Sn = 1/sqrt(2*log(n))
//   betastar = -log(-log(1-alpha))
//   criticalValue = Cn + Sn*betastar
//
// Python 精确计算: n=100, alpha=0.95
//   log(100)         = 4.6051701859880918
//   sqrt(2*log(100)) = 3.0348542587702929
//   log(pi)          = 1.1447298858494002
//   log(log(100))    = 1.5271796258079011
//   Cn               = 2.5946503339960034
//   Sn               = 0.3295051144911304
//   betastar         = -1.0971887003649483
//   criticalValue    = 2.2331210456638764
//
// 排幻觉 (2026-08-03 review):
//   原测试期望值 2.2331421269504335 基于手算 (log(log(100)) 精度不足 1.527179736
//   vs 精确值 1.5271796258...), 导致 Cn 偏差 2.6e-5. Python 验证确认正确值为
//   2.2331210456638764, C++ 实现与 Python 一致.
// =============================================================================
TEST(IntradayJumpTest, LeeMyklandCriticalValue) {
    Real cv = detail::lee_mykland_critical_value(100, 0.95);
    EXPECT_NEAR(cv, 2.2331210456638764, TOL);
}

// =============================================================================
// TEST 3: lee_mykland_critical_value — 不同 n 和 alpha
// R: betastar = -log(-log(1-alpha))
//   alpha=0.95 → 1-alpha=0.05 → betastar=-1.097
//   alpha=0.99 → 1-alpha=0.01 → betastar=-1.527 (更负)
//   cv = Cn + Sn*betastar, Sn > 0, 所以 alpha↑ → betastar↓ → cv↓
//
// 排幻觉 (2026-08-03 review):
//   原测试断言 cv99 > cv95 错误. 实际 alpha 越大 → betastar 越负 → cv 越小.
//   R 文档称 alpha 为 "significance level" 但公式用 1-alpha, 语义为 confidence
//   level: alpha↑ = 更确信无跳跃 = 临界值↓ = 更容易检测跳跃.
// =============================================================================
TEST(IntradayJumpTest, LeeMyklandCriticalValueN1000A99) {
    Real cv95 = detail::lee_mykland_critical_value(1000, 0.95);
    Real cv99 = detail::lee_mykland_critical_value(1000, 0.99);
    // alpha=0.99 的临界值应小于 alpha=0.95 (betastar 更负)
    EXPECT_LT(cv99, cv95);
    // 两者都应为正数
    EXPECT_GT(cv95, 0.0);
    EXPECT_GT(cv99, 0.0);
}

// =============================================================================
// TEST 4: intraday_jump_test 端到端 — 无跳跃 GBM 序列 (D14-D16 完整流程)
// 构造无跳跃的 GBM 价格序列, 验证 ztest 有限且量级合理
// =============================================================================
TEST(IntradayJumpTest, EndToEndNoJumps) {
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

    auto result = intraday_jump_test(prices, 10, 0.95, "rBPCov");

    // ztest 数量 = (n-1) - lookBackPeriod + 1 = n - lookBackPeriod
    // (returns 长度 = n-1, 去除 r[0]=0; vol 长度 = returns - lookBack + 1)
    EXPECT_EQ(result.ztest.size(), 200u - 10u);
    EXPECT_EQ(result.spotVol.size(), 200u - 10u);

    // 所有 ztest 应有限
    for (Real z : result.ztest) {
        EXPECT_TRUE(std::isfinite(z));
    }

    // criticalValue 应为正数
    EXPECT_GT(result.criticalValue, 0.0);

    // n 应等于原始观测数 (D16)
    EXPECT_EQ(result.n, 200);
}

// =============================================================================
// TEST 5: intraday_jump_test — vol 调整验证 (D14)
// R: vol$spot = sqrt(spot^2 / (lookBackPeriod-2)) where spot = sqrt(RBPVar)
//    → vol$spot = sqrt(RBPVar / (lookBackPeriod-2))
// 构造已知收益率序列, 验证 spotVol[0] = sqrt(rbp_var / (lookBack-2))
//
// 排幻觉 (2026-08-03 review):
//   原测试错误写成 sqrt(rbp_var^2 / (lookBack-2)), 混淆了 vol(=sqrt(RBPVar))
//   与 rbp_var(=RBPVar). R 源码 vol$spot^2 = RBPVar (非 RBPVar^2).
//   另: 原测试用 r[0]=0 构造 r 向量, 但实现用 returns (去除 r[0]=0),
//       预期值应基于 returns[0:lookBack-1] 而非 r[0:lookBack-1].
// =============================================================================
TEST(IntradayJumpTest, VolAdjustmentD14) {
    // 构造已知价格序列, 使收益率已知
    // P = [100, 101, 102, 101, 103, 102, 104, 105, 103, 106, 107, 108]
    std::vector<Real> prices = {100, 101, 102, 101, 103, 102, 104, 105, 103, 106, 107, 108};
    int lookBack = 5;

    auto result = intraday_jump_test(prices, lookBack, 0.95, "rBPCov");

    // 手动计算 returns (与实现一致: returns[i] = log(P[i+1]/P[i]), 无 r[0]=0)
    std::vector<Real> returns(prices.size() - 1);
    for (Size i = 0; i + 1 < prices.size(); ++i) {
        returns[i] = std::log(prices[i + 1] / prices[i]);
    }

    // 第一个窗口 returns[0:4] (lookBack=5 个元素)
    Real sum_rbp = 0.0;
    for (Size i = 0; i + 1 < static_cast<Size>(lookBack); ++i) {
        sum_rbp += std::fabs(returns[i]) * std::fabs(returns[i + 1]);
    }
    Real pi = std::acos(-1.0);
    Real rbp_var_first = (pi / 2.0) * sum_rbp;

    // D14 (修正): vol = sqrt(RBPVar / (lookBack-2))
    // R: vol$spot = sqrt((sqrt(RBPVar))^2 / (K-2)) = sqrt(RBPVar / (K-2))
    Real expected_vol = std::sqrt(rbp_var_first / (lookBack - 2));

    ASSERT_FALSE(result.spotVol.empty());
    EXPECT_NEAR(result.spotVol[0], expected_vol, TOL);
}

// =============================================================================
// TEST 6: 异常处理 — 空输入抛 invalid_argument
// =============================================================================
TEST(IntradayJumpTest, EmptyInputThrows) {
    std::vector<Real> empty_prices;
    EXPECT_THROW(
        intraday_jump_test(empty_prices, 10, 0.95, "rBPCov"),
        std::invalid_argument);
}

// =============================================================================
// TEST 7: 异常处理 — lookBackPeriod 不足 (n <= lookBackPeriod)
// =============================================================================
TEST(IntradayJumpTest, InsufficientLookBackThrows) {
    std::vector<Real> prices = {100.0, 101.0, 102.0, 103.0, 104.0};
    EXPECT_THROW(
        intraday_jump_test(prices, 10, 0.95, "rBPCov"),
        std::invalid_argument);
}

// =============================================================================
// TEST 8: 异常处理 — lookBackPeriod <= 2 (分母 lookBack-2 = 0)
// =============================================================================
TEST(IntradayJumpTest, LookBackLe2Throws) {
    std::vector<Real> prices = {100.0, 101.0, 102.0, 103.0, 104.0};
    EXPECT_THROW(
        intraday_jump_test(prices, 2, 0.95, "rBPCov"),
        std::invalid_argument);
}
