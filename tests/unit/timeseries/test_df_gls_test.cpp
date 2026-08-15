// =============================================================================
// test_df_gls_test.cpp - DF-GLS 检验测试 (18 用例, spec §3.4 测试矩阵)
//
// 基准: tests/unit/timeseries/unit_root_baseline.inc (arch 8.0.0 自动生成,
//       脚本 tests/fixtures/timeseries/verify_unit_root.py, 硬编码策略)
//
// 容差:
//   - statistic/nobs/lags: 1e-10; 临界值 1e-12; p 值 1e-12
// =============================================================================
#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <cmath>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/df_gls_test.hpp"
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

static void check_dfgls_case(const ur::DFGlsResult& r,
                             const ur::baseline::DfglsCase& c) {
    EXPECT_NEAR(r.statistic, c.stat, 1e-10);
    EXPECT_NEAR(r.p_value, c.p, 1e-12);
    EXPECT_NEAR(r.critical_value_1pct, c.cv1, 1e-12);
    EXPECT_NEAR(r.critical_value_5pct, c.cv5, 1e-12);
    EXPECT_NEAR(r.critical_value_10pct, c.cv10, 1e-12);
    EXPECT_EQ(r.n_obs, c.nobs);
    EXPECT_EQ(r.n_lags, c.lags);
}

// ---------------------------------------------------------------------------
// 1-4: arch 默认 AIC lag (rw→5, ar→0), 250 样本 c/ct (1e-10)
// ---------------------------------------------------------------------------
TEST(DfGlsTest, RwC) {
    check_dfgls_case(ur::df_gls_test(y_rw(), "c"), ur::baseline::DFGLS_CASES[0]);
}

TEST(DfGlsTest, RwCt) {
    check_dfgls_case(ur::df_gls_test(y_rw(), "ct"),
                     ur::baseline::DFGLS_CASES[1]);
}

TEST(DfGlsTest, ArC) {
    check_dfgls_case(ur::df_gls_test(y_ar(), "c"), ur::baseline::DFGLS_CASES[2]);
}

TEST(DfGlsTest, ArCt) {
    check_dfgls_case(ur::df_gls_test(y_ar(), "ct"),
                     ur::baseline::DFGLS_CASES[3]);
}

// ---------------------------------------------------------------------------
// 5-6: c̄ 与 ρ̄ 精确值 (U8/U9-修正)
// ---------------------------------------------------------------------------
TEST(DfGlsTest, CBarValues) {
    const auto c = ur::df_gls_test(y_rw(), "c");
    EXPECT_DOUBLE_EQ(c.c_bar, -7.0);
    const auto ct = ur::df_gls_test(y_rw(), "ct");
    EXPECT_DOUBLE_EQ(ct.c_bar, -13.5);
}

TEST(DfGlsTest, RhoBarValues) {
    // ρ̄ = 1 + c̄/T (非 c̄ 本身, U9-修正)
    const auto c = ur::df_gls_test(y_rw(), "c");
    EXPECT_NEAR(c.rho_bar, 1.0 - 7.0 / 250.0, 1e-15);   // 0.972
    EXPECT_GT(c.rho_bar, 0.9);
    EXPECT_LT(c.rho_bar, 1.0);
    const auto ct = ur::df_gls_test(y_rw(), "ct");
    EXPECT_NEAR(ct.rho_bar, 1.0 - 13.5 / 250.0, 1e-15);  // 0.946
}

// ---------------------------------------------------------------------------
// 7: lag 选择 = OLS-detrended + trend "n" + AIC (U20, 非 MAIC/Schwert)
// ---------------------------------------------------------------------------
TEST(DfGlsTest, LagSelectionIsAicOnOlsDetrended) {
    // rw → 5, ar → 0 (arch 默认 method="aic")
    EXPECT_EQ(ur::df_gls_test(y_rw(), "c").n_lags, 5);
    EXPECT_EQ(ur::df_gls_test(y_ar(), "c").n_lags, 0);
    EXPECT_EQ(ur::df_gls_test(y_rw(), "ct").n_lags, 5);
    EXPECT_EQ(ur::df_gls_test(y_ar(), "ct").n_lags, 0);
}

// ---------------------------------------------------------------------------
// 8: 显式 max_lag 覆盖自动选择
// ---------------------------------------------------------------------------
TEST(DfGlsTest, ExplicitMaxLag) {
    const auto r = ur::df_gls_test(y_rw(), "c", 4);
    EXPECT_EQ(r.n_lags, 4);
    EXPECT_EQ(r.n_obs, ur::baseline::T - 1 - 4);
}

// ---------------------------------------------------------------------------
// 9: 自动选择与显式最优 lag 数值一致 (同回归路径)
// ---------------------------------------------------------------------------
TEST(DfGlsTest, AutoMatchesExplicitOptimalLag) {
    const auto auto_r = ur::df_gls_test(y_rw(), "c");    // 自动选 5
    const auto manual_r = ur::df_gls_test(y_rw(), "c", 5);
    EXPECT_DOUBLE_EQ(auto_r.statistic, manual_r.statistic);
}

// ---------------------------------------------------------------------------
// 10: 回归样本语义 — n_obs = (T-1) - n_lags
// ---------------------------------------------------------------------------
TEST(DfGlsTest, RegressionNobsSemantics) {
    const auto r = ur::df_gls_test(y_rw(), "c");
    EXPECT_EQ(r.n_obs, ur::baseline::T - 1 - r.n_lags);
    // CV 用回归 nobs (244), 非原始 T (250) — 与基准 DFGLS_CASES[0].nobs 一致
    EXPECT_EQ(ur::baseline::DFGLS_CASES[0].nobs, 244);
}

// ---------------------------------------------------------------------------
// 11-13: 异常输入防护
// ---------------------------------------------------------------------------
TEST(DfGlsTest, InvalidTrendThrows) {
    EXPECT_THROW(ur::df_gls_test(y_rw(), "n"), std::invalid_argument);
    EXPECT_THROW(ur::df_gls_test(y_rw(), "xyz"), std::invalid_argument);
}

TEST(DfGlsTest, SampleTooSmallThrows) {
    EXPECT_THROW(ur::df_gls_test({1.0, 2.0, 3.0, 4.0}, "c"),
                 std::invalid_argument);
}

TEST(DfGlsTest, LagTooLargeThrows) {
    EXPECT_THROW(ur::df_gls_test(y_rw(), "c", 300), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 14: 拒绝方向 — H0: 单位根, 左尾拒绝
// ---------------------------------------------------------------------------
TEST(DfGlsTest, RejectNullDirection) {
    // rw/c: stat=0.51 > cv5 → 不拒绝单位根
    EXPECT_FALSE(ur::df_gls_test(y_rw(), "c").reject_null);
    // ar/c: stat=-8.20 < cv5 → 拒绝单位根 (平稳)
    const auto ar = ur::df_gls_test(y_ar(), "c");
    EXPECT_TRUE(ar.reject_null);
    EXPECT_LT(ar.statistic, ar.critical_value_5pct);
}

// ---------------------------------------------------------------------------
// 15: GLS detrending ≠ OLS detrending — DF-GLS 统计量与 ADF 不同 (U9)
// ---------------------------------------------------------------------------
TEST(DfGlsTest, GlsDiffersFromOlsAdf) {
    const auto dfgls = ur::df_gls_test(y_rw(), "c");
    const auto adf = ur::adf_test(y_rw(), "c", dfgls.n_lags);
    // DF-GLS 与基准一致; 同 lag 的 OLS-ADF 统计量应明显不同
    EXPECT_NEAR(dfgls.statistic, ur::baseline::DFGLS_CASES[0].stat, 1e-10);
    EXPECT_GT(std::fabs(dfgls.statistic - adf.statistic), 0.05);
}

// ---------------------------------------------------------------------------
// 16: 趋势序列 ct → 强拒绝单位根 (围绕趋势平稳)
// ---------------------------------------------------------------------------
TEST(DfGlsTest, TrendSeriesRejectsWithCt) {
    const auto r = ur::df_gls_test(y_tr(), "ct");
    EXPECT_TRUE(r.reject_null);
    EXPECT_LT(r.statistic, r.critical_value_1pct);
}

// ---------------------------------------------------------------------------
// 17: 默认参数 = trend "ct" + AIC 自动
// ---------------------------------------------------------------------------
TEST(DfGlsTest, DefaultParameters) {
    const auto r = ur::df_gls_test(y_rw());
    EXPECT_EQ(r.trend_spec, "ct");
    check_dfgls_case(r, ur::baseline::DFGLS_CASES[1]);
}

// ---------------------------------------------------------------------------
// 18: summary 方向
// ---------------------------------------------------------------------------
TEST(DfGlsTest, SummaryMentionsNullHypothesis) {
    const auto r = ur::df_gls_test(y_ar(), "c");
    EXPECT_NE(r.summary.find("unit root"), std::string::npos);
    EXPECT_NE(r.summary.find("stationary"), std::string::npos);
}
