// =============================================================================
// test_jump_test_diagnostics.cpp - Phase 7A Wave 3 跳跃检验多重修正测试
//
// 10 用例:
//   Bonferroni 修正 (3) + Benjamini-Hochberg 修正 (3) + 跳跃检验联合 (4)
//
// 排幻觉点覆盖:
//   H23 (BH 修正控制 FDR, 非 FWER; Bonferroni 控制 FWER)
//
// 教材锚点: Bonferroni 1936, Benjamini-Hochberg 1995
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/hfecon/tests/jump_test_diagnostics.hpp"

using cpphub::v1::hfecon::multiple_test_correction;
using cpphub::v1::hfecon::jump_test_diagnostics;
using cpphub::v1::hfecon::MultipleTestCorrectionResult;
using cpphub::v1::hfecon::JumpTestDiagnosticsResult;
using cpphub::v1::Real;
using cpphub::v1::Size;

// =============================================================================
// Bonferroni 修正 (3 用例)
// =============================================================================

// --- Bonferroni 1: 基本修正公式 adjusted_p = min(m * p, 1.0) ---
TEST(BonferroniCorrection, BasicFormula) {
    std::vector<Real> p = {0.01, 0.04, 0.03};
    auto res = multiple_test_correction(
        p, MultipleTestCorrectionResult::Method::Bonferroni, 0.05);

    ASSERT_EQ(res.adjusted_p_values.size(), 3u);
    // m=3: adjusted = min(3*p, 1)
    EXPECT_NEAR(res.adjusted_p_values[0], 0.03, 1e-10);  // 3*0.01
    EXPECT_NEAR(res.adjusted_p_values[1], 0.12, 1e-10);  // 3*0.04
    EXPECT_NEAR(res.adjusted_p_values[2], 0.09, 1e-10);  // 3*0.03

    // reject if adjusted < alpha=0.05
    EXPECT_TRUE(res.reject_null[0]);   // 0.03 < 0.05
    EXPECT_FALSE(res.reject_null[1]);  // 0.12 > 0.05
    EXPECT_FALSE(res.reject_null[2]);  // 0.09 > 0.05
    EXPECT_EQ(res.n_rejections, 1u);
}

// --- Bonferroni 2: 截断到 1.0 ---
TEST(BonferroniCorrection, TruncationToOne) {
    std::vector<Real> p = {0.5, 0.8};
    auto res = multiple_test_correction(
        p, MultipleTestCorrectionResult::Method::Bonferroni, 0.05);

    // m=2: adjusted = min(2*p, 1)
    EXPECT_NEAR(res.adjusted_p_values[0], 1.0, 1e-10);  // min(1.0, 1.0)
    EXPECT_NEAR(res.adjusted_p_values[1], 1.0, 1e-10);  // min(1.6, 1.0)
    EXPECT_EQ(res.n_rejections, 0u);
}

// --- Bonferroni 3: 全拒绝 ---
TEST(BonferroniCorrection, AllRejected) {
    std::vector<Real> p = {0.001, 0.002, 0.003, 0.004};
    auto res = multiple_test_correction(
        p, MultipleTestCorrectionResult::Method::Bonferroni, 0.05);

    // m=4: adjusted = 4*p = {0.004, 0.008, 0.012, 0.016}, 全 < 0.05
    for (Size i = 0; i < 4; ++i) {
        EXPECT_TRUE(res.reject_null[i]) << "i=" << i;
    }
    EXPECT_EQ(res.n_rejections, 4u);
}

// =============================================================================
// Benjamini-Hochberg 修正 (3 用例)
// 排幻觉点 H23: BH 控制 FDR, 非 FWER
// =============================================================================

// --- BH 1: 基本BH修正, 验证排序和 adjusted p ---
TEST(BHCorrection, BasicFormula) {
    // p 值故意乱序: {0.04, 0.01, 0.03}
    // 排序后: p_(1)=0.01, p_(2)=0.03, p_(3)=0.04
    // BH adjusted (从大到小累积最小):
    //   k=3: m*p_(3)/3 = 3*0.04/3 = 0.04, cum_min=0.04
    //   k=2: m*p_(2)/2 = 3*0.03/2 = 0.045, cum_min=min(0.04, 0.045)=0.04
    //   k=1: m*p_(1)/1 = 3*0.01/1 = 0.03, cum_min=min(0.04, 0.03)=0.03
    // adjusted = {0.04 (for p=0.04), 0.03 (for p=0.01), 0.04 (for p=0.03)}
    std::vector<Real> p = {0.04, 0.01, 0.03};
    auto res = multiple_test_correction(
        p, MultipleTestCorrectionResult::Method::BenjaminiHochberg, 0.05);

    ASSERT_EQ(res.adjusted_p_values.size(), 3u);
    // p=0.04 (rank 3): adjusted = 0.04
    EXPECT_NEAR(res.adjusted_p_values[0], 0.04, 1e-10);
    // p=0.01 (rank 1): adjusted = 0.03
    EXPECT_NEAR(res.adjusted_p_values[1], 0.03, 1e-10);
    // p=0.03 (rank 2): adjusted = 0.04 (累积最小值)
    EXPECT_NEAR(res.adjusted_p_values[2], 0.04, 1e-10);
}

// --- BH 2: BH 比 Bonferroni 更强力 (更多拒绝) ---
// 排幻觉点 H23: BH 控制 FDR, 允许更多拒绝
TEST(BHCorrection, MorePowerfulThanBonferroni) {
    std::vector<Real> p = {0.01, 0.02, 0.03, 0.04};
    auto res_bh = multiple_test_correction(
        p, MultipleTestCorrectionResult::Method::BenjaminiHochberg, 0.05);
    auto res_bonf = multiple_test_correction(
        p, MultipleTestCorrectionResult::Method::Bonferroni, 0.05);

    // BH 拒绝数 >= Bonferroni 拒绝数
    EXPECT_GE(res_bh.n_rejections, res_bonf.n_rejections);

    // 具体验证:
    // Bonferroni: adjusted = 4*p = {0.04, 0.08, 0.12, 0.16}, 拒绝 {0.04<0.05}
    // BH: 排序 p_(1)=0.01, p_(2)=0.02, p_(3)=0.03, p_(4)=0.04
    //   k=4: 4*0.04/4=0.04, cum_min=0.04
    //   k=3: 4*0.03/3=0.04, cum_min=0.04
    //   k=2: 4*0.02/2=0.04, cum_min=0.04
    //   k=1: 4*0.01/1=0.04, cum_min=0.04
    //   adjusted all = 0.04 < 0.05 → 全拒绝
    EXPECT_EQ(res_bonf.n_rejections, 1u);  // 仅 p=0.01
    EXPECT_EQ(res_bh.n_rejections, 4u);    // 全拒绝
}

// --- BH 3: 边界情况 - 单个 p 值 ---
TEST(BHCorrection, SinglePValue) {
    std::vector<Real> p = {0.03};
    auto res = multiple_test_correction(
        p, MultipleTestCorrectionResult::Method::BenjaminiHochberg, 0.05);

    // m=1: adjusted = 1*0.03/1 = 0.03
    EXPECT_NEAR(res.adjusted_p_values[0], 0.03, 1e-10);
    EXPECT_TRUE(res.reject_null[0]);
    EXPECT_EQ(res.n_rejections, 1u);
}

// =============================================================================
// 跳跃检验联合诊断 (4 用例)
// =============================================================================

// --- 联合 1: 无跳跃 (所有统计量小) → consensus_jumps=0 ---
TEST(JumpTestJoint, NoJumpsDetected) {
    const Size n_days = 5;
    std::vector<Real> bns(n_days, 0.5), aj(n_days, 0.3),
        jo(n_days, 0.4), rank(n_days, 0.2);

    auto res = jump_test_diagnostics(bns, aj, jo, rank, 0.05);

    // 所有 |Z| < 1.96 → p > 0.05, 不拒绝
    EXPECT_EQ(res.consensus_jumps, 0u);
    // 全不拒绝 → consistent=true
    EXPECT_TRUE(res.consistent);
}

// --- 联合 2: 全部跳跃 (所有统计量大) → consensus_jumps=n_days ---
TEST(JumpTestJoint, AllJumpsDetected) {
    const Size n_days = 3;
    std::vector<Real> bns(n_days, 5.0), aj(n_days, 4.5),
        jo(n_days, 5.5), rank(n_days, 4.8);

    auto res = jump_test_diagnostics(bns, aj, jo, rank, 0.05);

    // 所有 |Z| > 3.5 → p < 0.001, Bonferroni 修正后仍拒绝
    EXPECT_EQ(res.consensus_jumps, n_days);
    // 全拒绝 → consistent=true
    EXPECT_TRUE(res.consistent);
}

// --- 联合 3: 不一致 (部分检验拒绝) → consistent=false ---
TEST(JumpTestJoint, InconsistentTests) {
    const Size n_days = 2;
    // Day 1: BNS/AJ 大 (拒绝), JO/Rank 小 (不拒绝) → 不一致
    // Day 2: 全小 (不拒绝) → 一致
    std::vector<Real> bns = {5.0, 0.5}, aj = {4.5, 0.3},
        jo = {0.4, 0.4}, rank = {0.2, 0.2};

    auto res = jump_test_diagnostics(bns, aj, jo, rank, 0.05);

    // Day 1: 2/4 拒绝 → 非 0 或 4 → 不一致
    // Day 2: 0/4 拒绝 → 一致
    // 总体: 不一致
    EXPECT_FALSE(res.consistent);
}

// --- 联合 4: 异常情况 - 空输入抛异常 ---
TEST(JumpTestJoint, EmptyInputThrows) {
    std::vector<Real> empty;
    EXPECT_THROW(jump_test_diagnostics(empty, empty, empty, empty),
                 std::invalid_argument);
}

// =============================================================================
// 异常情况测试
// =============================================================================

// --- 异常 1: 空 p 值向量 ---
TEST(MultipleTestCorrectionExceptions, EmptyPValuesThrows) {
    std::vector<Real> empty;
    EXPECT_THROW(
        multiple_test_correction(empty,
            MultipleTestCorrectionResult::Method::Bonferroni),
        std::invalid_argument);
}

// --- 异常 2: alpha 越界 ---
TEST(MultipleTestCorrectionExceptions, InvalidAlphaThrows) {
    std::vector<Real> p = {0.01, 0.05};
    EXPECT_THROW(
        multiple_test_correction(p,
            MultipleTestCorrectionResult::Method::Bonferroni, 0.0),
        std::invalid_argument);
    EXPECT_THROW(
        multiple_test_correction(p,
            MultipleTestCorrectionResult::Method::Bonferroni, 1.0),
        std::invalid_argument);
}

// --- 异常 3: 跳跃检验大小不匹配 ---
TEST(JumpTestJointExceptions, SizeMismatchThrows) {
    std::vector<Real> bns(3, 1.0), aj(2, 1.0), jo(3, 1.0), rank(3, 1.0);
    EXPECT_THROW(jump_test_diagnostics(bns, aj, jo, rank),
                 std::invalid_argument);
}
