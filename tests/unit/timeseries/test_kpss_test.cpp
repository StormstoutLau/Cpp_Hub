// =============================================================================
// test_kpss_test.cpp - KPSS 平稳性检验测试 (15 用例, spec §3.3 测试矩阵)
//
// 基准: tests/unit/timeseries/unit_root_baseline.inc (arch 8.0.0 自动生成,
//       脚本 tests/fixtures/timeseries/verify_unit_root.py, 硬编码策略)
//
// 容差:
//   - statistic/lags: 1e-10 (同算法对照, 实际逐位一致)
//   - 临界值: 1e-12; p 值: 1e-12 (插值方向性由 test_mackinnon_cv 覆盖)
//   - tiny 手算: 1e-15 / Hobijn 手算: 精确
// =============================================================================
#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <cmath>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/kpss_test.hpp"
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

static void check_kpss_case(const ur::KPSSResult& r,
                            const ur::baseline::KpssCase& c) {
    EXPECT_NEAR(r.statistic, c.stat, 1e-10);
    EXPECT_NEAR(r.p_value, c.p, 1e-12);
    EXPECT_NEAR(r.critical_value_1pct, c.cv1, 1e-12);
    EXPECT_NEAR(r.critical_value_5pct, c.cv5, 1e-12);
    EXPECT_NEAR(r.critical_value_10pct, c.cv10, 1e-12);
    EXPECT_EQ(r.bandwidth, c.lags);
}

// ---------------------------------------------------------------------------
// 1-4: Hobijn 自动带宽 (rw→10, ar→6), 250 样本 (1e-10)
// ---------------------------------------------------------------------------
TEST(KpssTest, RwC) {
    check_kpss_case(ur::kpss_test(y_rw(), "c"), ur::baseline::KPSS_CASES[0]);
}

TEST(KpssTest, RwCt) {
    check_kpss_case(ur::kpss_test(y_rw(), "ct"), ur::baseline::KPSS_CASES[1]);
}

TEST(KpssTest, ArC) {
    check_kpss_case(ur::kpss_test(y_ar(), "c"), ur::baseline::KPSS_CASES[2]);
}

TEST(KpssTest, ArCt) {
    check_kpss_case(ur::kpss_test(y_ar(), "ct"), ur::baseline::KPSS_CASES[3]);
}

// ---------------------------------------------------------------------------
// 5: legacy 模式 = 显式 Schwert 带宽 16 (U11)
// ---------------------------------------------------------------------------
TEST(KpssTest, LegacySchwertBandwidth16) {
    check_kpss_case(ur::kpss_test(y_rw(), "c", ur::baseline::SCHWERT_250),
                    ur::baseline::KPSS_CASES[4]);
}

// ---------------------------------------------------------------------------
// 6: tiny 手算 — y=[1..10], trend c, lags=2 (1e-12)
// ---------------------------------------------------------------------------
TEST(KpssTest, TinyHandComputed) {
    const auto r = ur::kpss_test(tiny_y(), "c", 2);
    EXPECT_NEAR(r.statistic, ur::baseline::KPSS_TINY_STAT, 1e-12);
    EXPECT_EQ(r.bandwidth, 2);
    EXPECT_EQ(r.n_obs, ur::baseline::TINY_N);
}

// ---------------------------------------------------------------------------
// 7: tiny LAM 分量白盒 — 退化均值后 Bartlett LRV 逐位一致 (1e-15)
// ---------------------------------------------------------------------------
TEST(KpssTest, TinyLamComponent) {
    const auto y = tiny_y();
    const Size n = ur::baseline::TINY_N;
    Real mean = 0.0;
    for (Real v : y) mean += v;
    mean /= static_cast<Real>(n);
    std::vector<Real> u(n);
    for (Size i = 0; i < n; ++i) u[i] = y[i] - mean;
    // 18.22 量级下 1 ULP ≈ 3.55e-15, 1e-13 ≈ 28 ULP 仍远严于 spec 1e-10
    EXPECT_NEAR(ur::long_run_variance(u, 2), ur::baseline::KPSS_TINY_LAM,
                1e-13);
}

// ---------------------------------------------------------------------------
// 8: Hobijn 手算小样本 — u=[1,1,1,1]: covlags=1, s0=2.5, s1=1.5,
//    γ̂=1.1447·0.6^(2/3), lags=int(γ̂·4^(1/3))=1 (逐位可手算)
// ---------------------------------------------------------------------------
TEST(KpssTest, HobijnHandComputedSmallCase) {
    EXPECT_EQ(ur::hobijn_bandwidth({1.0, 1.0, 1.0, 1.0}), 1);
    // 交替序列 [1,-1,...]: s0 = 1 + (-7/4) < 0 → 退化抛异常
    EXPECT_THROW(ur::hobijn_bandwidth({1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0,
                                       -1.0}),
                 std::runtime_error);
    EXPECT_THROW(ur::hobijn_bandwidth({}), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 9: Hobijn 与 legacy 带宽不同 → 统计量不同 (U11 方向性)
// ---------------------------------------------------------------------------
TEST(KpssTest, HobijnDiffersFromLegacy) {
    const auto hobijn = ur::kpss_test(y_rw(), "c");
    const auto legacy = ur::kpss_test(y_rw(), "c", ur::baseline::SCHWERT_250);
    EXPECT_NE(hobijn.bandwidth, legacy.bandwidth);
    EXPECT_NE(hobijn.statistic, legacy.statistic);
}

// ---------------------------------------------------------------------------
// 10: 拒绝方向 — H0: 平稳性, 右尾拒绝 (U12, 与 ADF 相反)
// ---------------------------------------------------------------------------
TEST(KpssTest, RejectNullDirection) {
    // rw: LM=1.31 > cv5=0.4614 → 拒绝平稳性 (非平稳)
    const auto rw = ur::kpss_test(y_rw(), "c");
    EXPECT_TRUE(rw.reject_null);
    EXPECT_GT(rw.statistic, rw.critical_value_5pct);
    // ar: LM=0.114 < cv5 → 不拒绝平稳性
    const auto ar = ur::kpss_test(y_ar(), "c");
    EXPECT_FALSE(ar.reject_null);
    EXPECT_LT(ar.statistic, ar.critical_value_5pct);
}

// ---------------------------------------------------------------------------
// 11-13: 异常输入防护
// ---------------------------------------------------------------------------
TEST(KpssTest, InvalidTrendThrows) {
    // KPSS 不支持 "n"/"nc"/"ctt" — 只有 c/ct
    EXPECT_THROW(ur::kpss_test(y_rw(), "n"), std::invalid_argument);
    EXPECT_THROW(ur::kpss_test(y_rw(), "xyz"), std::invalid_argument);
}

TEST(KpssTest, SampleTooSmallThrows) {
    EXPECT_THROW(ur::kpss_test({1.0, 2.0}, "c"), std::invalid_argument);
}

TEST(KpssTest, BandwidthTooLargeThrows) {
    EXPECT_THROW(ur::kpss_test(y_rw(), "c", 300), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 14: 常数序列退化防护 — 残差全零 → Hobijn s0≤0 → runtime_error
// ---------------------------------------------------------------------------
TEST(KpssTest, DegenerateConstantSeriesThrows) {
    EXPECT_THROW(ur::kpss_test({5.0, 5.0, 5.0, 5.0, 5.0}, "c"),
                 std::runtime_error);
}

// ---------------------------------------------------------------------------
// 15: summary 明确 H0/H1 方向 (U12)
// ---------------------------------------------------------------------------
TEST(KpssTest, SummaryMentionsDirection) {
    const auto r = ur::kpss_test(y_ar(), "ct");
    EXPECT_NE(r.summary.find("stationary"), std::string::npos);
    EXPECT_NE(r.summary.find("unit root"), std::string::npos);
    // 默认 trend = "c"
    EXPECT_EQ(ur::kpss_test(y_ar()).trend_spec, "c");
}
