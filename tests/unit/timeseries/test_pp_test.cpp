// =============================================================================
// test_pp_test.cpp - Phillips-Perron 单位根检验测试 (15 用例, spec §3.2 测试矩阵)
//
// 基准: tests/unit/timeseries/unit_root_baseline.inc (arch 8.0.0 自动生成,
//       脚本 tests/fixtures/timeseries/verify_unit_root.py, 硬编码策略)
//
// 容差:
//   - statistic/nobs/bandwidth: 1e-10 (同算法对照, 实际逐位一致)
//   - 临界值: 1e-12; p 值: 1e-12
//   - tiny 手算分量: 1e-15 (整数数据, 精确算术)
// =============================================================================
#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <cmath>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/pp_test.hpp"
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
static std::vector<Real> tiny_y() {
    return {ur::baseline::TINY_Y, ur::baseline::TINY_Y + ur::baseline::TINY_N};
}

static void check_pp_case(const ur::PPResult& r,
                          const ur::baseline::PpCase& c) {
    EXPECT_NEAR(r.statistic, c.stat, 1e-10);
    EXPECT_NEAR(r.p_value, c.p, 1e-12);
    EXPECT_NEAR(r.critical_value_1pct, c.cv1, 1e-12);
    EXPECT_NEAR(r.critical_value_5pct, c.cv5, 1e-12);
    EXPECT_NEAR(r.critical_value_10pct, c.cv10, 1e-12);
    EXPECT_EQ(r.n_obs, c.nobs);
    EXPECT_EQ(r.bandwidth, c.lags);
}

// ---------------------------------------------------------------------------
// 1-4: 默认带宽 = Schwert(250)=16, 250 样本 rw/ar × c/ct (1e-10)
// ---------------------------------------------------------------------------
TEST(PpTest, RwC) {
    check_pp_case(ur::pp_test(y_rw(), "c"), ur::baseline::PP_CASES[0]);
}

TEST(PpTest, RwCt) {
    check_pp_case(ur::pp_test(y_rw(), "ct"), ur::baseline::PP_CASES[1]);
}

TEST(PpTest, ArC) {
    check_pp_case(ur::pp_test(y_ar(), "c"), ur::baseline::PP_CASES[2]);
}

TEST(PpTest, ArCt) {
    check_pp_case(ur::pp_test(y_ar(), "ct"), ur::baseline::PP_CASES[3]);
}

// ---------------------------------------------------------------------------
// 5: 显式带宽 lags=8 (U-ADR8)
// ---------------------------------------------------------------------------
TEST(PpTest, ExplicitBandwidth8) {
    check_pp_case(ur::pp_test(y_rw(), "c", 8), ur::baseline::PP_CASES[4]);
}

// ---------------------------------------------------------------------------
// 6: tiny 手算 — y=[1..10], trend n, lags=2, 全链路统计量 (1e-12)
// ---------------------------------------------------------------------------
TEST(PpTest, TinyHandComputedStatistic) {
    const auto r = ur::pp_test(tiny_y(), "n", 2);
    EXPECT_NEAR(r.statistic, ur::baseline::PP_TINY_STAT, 1e-12);
    EXPECT_EQ(r.n_obs, ur::baseline::TINY_N - 1);
    EXPECT_EQ(r.bandwidth, 2);
}

// ---------------------------------------------------------------------------
// 7: tiny 分量白盒 — rho/bse/gamma0/lam2/s2 与 arch 逐位一致 (U5-sigma2:
//    第一项 gamma0=SSR/n 无 df 修正, 第二项 s2=SSR/(n-k) df 修正)
// ---------------------------------------------------------------------------
TEST(PpTest, TinyComponentsMatchArch) {
    const auto y = tiny_y();
    const Size n = ur::baseline::TINY_N - 1;
    std::vector<Real> lhs(n);
    std::vector<std::vector<Real>> X(n, std::vector<Real>(1));
    for (Size i = 0; i < n; ++i) {
        lhs[i] = y[i + 1] - y[i];
        X[i][0] = y[i];
    }
    const auto ols = ur::detail::ols_fit(lhs, X);
    // rho = 1 + γ̂
    EXPECT_NEAR(1.0 + ols.beta[0], ur::baseline::PP_TINY_RHO, 1e-15);
    EXPECT_NEAR(ols.bse[0], ur::baseline::PP_TINY_BSE, 1e-15);
    // gamma0 = SSR/n (无 df 修正)
    EXPECT_NEAR(ols.ssr / static_cast<Real>(n), ur::baseline::PP_TINY_GAMMA0,
                1e-15);
    // s2 = SSR/(n-k) (df 修正)
    EXPECT_NEAR(ols.ssr / static_cast<Real>(n - ols.n_params),
                ur::baseline::PP_TINY_S2, 1e-15);
    // lam2 = NW Bartlett 长期方差, bandwidth=2
    EXPECT_NEAR(ur::long_run_variance(ols.resid, 2), ur::baseline::PP_TINY_LAM2,
                1e-15);
    // 两套方差定义必须不同 (排 U5-sigma2: 混用单一定义则失败)
    EXPECT_NE(ur::baseline::PP_TINY_GAMMA0, ur::baseline::PP_TINY_S2);
}

// ---------------------------------------------------------------------------
// 8: 默认带宽 = Schwert 规则 (U5, 与 ADF 同公式)
// ---------------------------------------------------------------------------
TEST(PpTest, DefaultBandwidthIsSchwert) {
    const auto r = ur::pp_test(y_rw(), "c");
    EXPECT_EQ(r.bandwidth, ur::baseline::SCHWERT_250);
}

// ---------------------------------------------------------------------------
// 9: PP 无 ADF 式上限保护 — T=10: Schwert=7 直接生效
//    (ADF 同样本会被 cap 到 (T-1)/2-1-tc=3)
//    注: trend 用 "n" — tiny_y=[1..10] 为完美线性, "c" 下 SSR=0 触发退化防护
// ---------------------------------------------------------------------------
TEST(PpTest, NoUpperCapUnlikeAdf) {
    const auto r = ur::pp_test(tiny_y(), "n");
    EXPECT_EQ(r.bandwidth, ur::schwert_lag(ur::baseline::TINY_N));  // 7
    EXPECT_EQ(ur::schwert_lag(ur::baseline::TINY_N), 7);
}

// ---------------------------------------------------------------------------
// 10: 拒绝方向 — H0: 单位根, 左尾拒绝
// ---------------------------------------------------------------------------
TEST(PpTest, RejectNullDirection) {
    EXPECT_FALSE(ur::pp_test(y_rw(), "c").reject_null);
    const auto ar = ur::pp_test(y_ar(), "c");
    EXPECT_TRUE(ar.reject_null);
    EXPECT_LT(ar.statistic, ar.critical_value_5pct);
}

// ---------------------------------------------------------------------------
// 11-13: 异常输入防护
// ---------------------------------------------------------------------------
TEST(PpTest, InvalidTrendThrows) {
    EXPECT_THROW(ur::pp_test(y_rw(), "xyz"), std::invalid_argument);
}

TEST(PpTest, SampleTooSmallThrows) {
    EXPECT_THROW(ur::pp_test({1.0, 2.0, 3.0}, "c"), std::invalid_argument);
}

TEST(PpTest, BandwidthTooLargeThrows) {
    EXPECT_THROW(ur::pp_test(y_rw(), "c", 300), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 14: 常数序列退化防护 — 残差全零 → runtime_error
//     (显式 bandwidth=2 避开默认 Schwert(5)=6 >= nobs 的前置 invalid_argument)
// ---------------------------------------------------------------------------
TEST(PpTest, DegenerateConstantSeriesThrows) {
    EXPECT_THROW(ur::pp_test({5.0, 5.0, 5.0, 5.0, 5.0}, "n", 2),
                 std::runtime_error);
}

// ---------------------------------------------------------------------------
// 15: "n" 与 "nc" 别名等价
// ---------------------------------------------------------------------------
TEST(PpTest, NcAliasMatchesN) {
    const auto a = ur::pp_test(y_rw(), "n", 8);
    const auto b = ur::pp_test(y_rw(), "nc", 8);
    EXPECT_DOUBLE_EQ(a.statistic, b.statistic);
    EXPECT_EQ(b.trend_spec, "nc");
}
