// =============================================================================
// test_adf_test.cpp - ADF 单位根检验测试 (20 用例, spec §3.1 测试矩阵)
//
// 基准: tests/unit/timeseries/unit_root_baseline.inc (arch 8.0.0 自动生成,
//       脚本 tests/fixtures/timeseries/verify_unit_root.py, 硬编码策略)
//
// 容差:
//   - statistic/nobs/n_lags: 1e-10 (同算法对照, 实际逐位一致)
//   - 临界值: 1e-12 (MacKinnon response surface, 纯多项式)
//   - p 值: 1e-12 (erfc vs scipy ndtr, 中心区 ~1e-16)
// =============================================================================
#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <cmath>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/adf_test.hpp"
#include "unit_root_baseline.inc"

namespace ur = cpphub::v1::timeseries::unit_root;
using cpphub::Real;
using cpphub::Size;

static std::vector<Real> y_rw() {
    return {ur::baseline::Y_RW, ur::baseline::Y_RW + ur::baseline::T};
}
static std::vector<Real> y_ar() {
    return {ur::baseline::Y_AR, ur::baseline::Y_AR + ur::baseline::T};
}
static std::vector<Real> y_tr() {
    return {ur::baseline::Y_TR, ur::baseline::Y_TR + ur::baseline::T};
}

// 固定 lag=16 案例统一断言: stat/p/cv/nobs/n_lags
static void check_adf_case(const ur::ADFResult& r,
                           const ur::baseline::AdfCase& c) {
    EXPECT_NEAR(r.statistic, c.stat, 1e-10);
    EXPECT_NEAR(r.p_value, c.p, 1e-12);
    EXPECT_NEAR(r.critical_value_1pct, c.cv1, 1e-12);
    EXPECT_NEAR(r.critical_value_5pct, c.cv5, 1e-12);
    EXPECT_NEAR(r.critical_value_10pct, c.cv10, 1e-12);
    EXPECT_EQ(r.n_obs, c.nobs);
    EXPECT_EQ(r.n_lags, c.n_lags);
}

// ---------------------------------------------------------------------------
// 1-6: ADF_FIXED — 250 样本, 固定 lag=16 (Schwert(250)), 1e-10 对照
// ---------------------------------------------------------------------------
TEST(AdfTest, FixedLagRwN) {
    check_adf_case(ur::adf_test(y_rw(), "n", 16), ur::baseline::ADF_FIXED[0]);
}

TEST(AdfTest, FixedLagRwC) {
    check_adf_case(ur::adf_test(y_rw(), "c", 16), ur::baseline::ADF_FIXED[1]);
}

TEST(AdfTest, FixedLagRwCt) {
    check_adf_case(ur::adf_test(y_rw(), "ct", 16), ur::baseline::ADF_FIXED[2]);
}

TEST(AdfTest, FixedLagArN) {
    check_adf_case(ur::adf_test(y_ar(), "n", 16), ur::baseline::ADF_FIXED[3]);
}

TEST(AdfTest, FixedLagArC) {
    check_adf_case(ur::adf_test(y_ar(), "c", 16), ur::baseline::ADF_FIXED[4]);
}

TEST(AdfTest, FixedLagArCt) {
    check_adf_case(ur::adf_test(y_ar(), "ct", 16), ur::baseline::ADF_FIXED[5]);
}

// ---------------------------------------------------------------------------
// 7-10: ADF_AIC — arch lags=None, method='aic' (固定截断评估样本)
// ---------------------------------------------------------------------------
TEST(AdfTest, AicRwN) {
    const auto r = ur::adf_test(y_rw(), "n", 0, true);
    check_adf_case(r, ur::baseline::ADF_AIC[0]);
}

TEST(AdfTest, AicRwC) {
    check_adf_case(ur::adf_test(y_rw(), "c", 0, true),
                   ur::baseline::ADF_AIC[1]);
}

TEST(AdfTest, AicArN) {
    check_adf_case(ur::adf_test(y_ar(), "n", 0, true),
                   ur::baseline::ADF_AIC[2]);
}

TEST(AdfTest, AicArC) {
    check_adf_case(ur::adf_test(y_ar(), "c", 0, true),
                   ur::baseline::ADF_AIC[3]);
}

// ---------------------------------------------------------------------------
// 11: 默认 lag=0 => Schwert 自动 + 上限保护 (U1)
// ---------------------------------------------------------------------------
TEST(AdfTest, DefaultLagIsSchwert) {
    const auto r = ur::adf_test(y_rw(), "n");
    EXPECT_EQ(r.n_lags, ur::baseline::SCHWERT_250);
    EXPECT_NEAR(r.statistic, ur::baseline::ADF_FIXED[0].stat, 1e-10);
}

// ---------------------------------------------------------------------------
// 12: 拒绝方向 — H0: 单位根, 左尾拒绝 (U4)
// ---------------------------------------------------------------------------
TEST(AdfTest, RejectNullDirection) {
    // rw/n: stat=+0.58 > cv5 → 不拒绝
    const auto rw = ur::adf_test(y_rw(), "n", 16);
    EXPECT_FALSE(rw.reject_null);
    EXPECT_GT(rw.statistic, rw.critical_value_5pct);
    // ar/n: stat=-2.78 < cv5=-1.94 → 拒绝
    const auto ar = ur::adf_test(y_ar(), "n", 16);
    EXPECT_TRUE(ar.reject_null);
    EXPECT_LT(ar.statistic, ar.critical_value_5pct);
}

// ---------------------------------------------------------------------------
// 13: trend auto — 确定性趋势序列应选 "ct" (U2, 方向性验证)
// ---------------------------------------------------------------------------
TEST(AdfTest, TrendAutoSelectsCtForTrendSeries) {
    const auto r = ur::adf_test(y_tr(), "auto", 4);
    EXPECT_EQ(r.trend_spec, "ct");
    // auto 结果与显式 ct 一致 (数值路径相同)
    const auto explicit_ct = ur::adf_test(y_tr(), "ct", 4);
    EXPECT_NEAR(r.statistic, explicit_ct.statistic, 1e-12);
}

// ---------------------------------------------------------------------------
// 14-16: 异常输入防护
// ---------------------------------------------------------------------------
TEST(AdfTest, InvalidTrendThrows) {
    EXPECT_THROW(ur::adf_test(y_rw(), "xyz", 4), std::invalid_argument);
}

TEST(AdfTest, SampleTooSmallThrows) {
    EXPECT_THROW(ur::adf_test({1.0, 2.0, 3.0}, "c", 0),
                 std::invalid_argument);
}

TEST(AdfTest, LagTooLargeThrows) {
    EXPECT_THROW(ur::adf_test(y_rw(), "c", 300), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 17: 系数布局 — [trend..., γ, δ₁..δ_p], γ 位置 = trend_cols
// ---------------------------------------------------------------------------
TEST(AdfTest, CoefficientsLayout) {
    const Size p = 2;
    const auto r = ur::adf_test(y_ar(), "c", p);
    EXPECT_EQ(r.coefficients.size(), 1 + 1 + p);  // [const, y_{t-1}, dy1, dy2]
    EXPECT_EQ(r.std_errors.size(), r.coefficients.size());
    EXPECT_EQ(r.n_obs, ur::baseline::T - 1 - p);
    // 平稳 AR: γ < 0 (U4 方向)
    EXPECT_LT(r.coefficients[1], 0.0);
    // 自身一致性: statistic = γ / SE(γ)
    EXPECT_NEAR(r.statistic, r.coefficients[1] / r.std_errors[1], 1e-10);
}

// ---------------------------------------------------------------------------
// 18: AIC 评估样本固定截断语义 — 最终回归 nobs = (T-1) - p, 非 max_lag
// ---------------------------------------------------------------------------
TEST(AdfTest, AicSampleTruncationSemantics) {
    // rw/n: AIC 选 p=5 → nobs=244 (非评估样本 233)
    EXPECT_EQ(ur::baseline::ADF_AIC[0].nobs, ur::baseline::T - 1 - 5);
    EXPECT_EQ(ur::baseline::ADF_AIC[0].n_lags, 5);
    // ar/n: AIC 选 p=0 → nobs=249
    EXPECT_EQ(ur::baseline::ADF_AIC[2].nobs, ur::baseline::T - 1);
}

// ---------------------------------------------------------------------------
// 19: "n" 与 "nc" 别名等价
// ---------------------------------------------------------------------------
TEST(AdfTest, NcAliasMatchesN) {
    const auto a = ur::adf_test(y_rw(), "n", 16);
    const auto b = ur::adf_test(y_rw(), "nc", 16);
    EXPECT_NEAR(a.statistic, b.statistic, 0.0);
    EXPECT_NEAR(a.p_value, b.p_value, 0.0);
    EXPECT_EQ(b.trend_spec, "nc");
}

// ---------------------------------------------------------------------------
// 20: summary 明确 H0/H1 方向
// ---------------------------------------------------------------------------
TEST(AdfTest, SummaryMentionsNullHypothesis) {
    const auto r = ur::adf_test(y_rw(), "c", 4);
    EXPECT_NE(r.summary.find("unit root"), std::string::npos);
    EXPECT_NE(r.summary.find("stationary"), std::string::npos);
}
