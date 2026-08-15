// =============================================================================
// test_unit_root_common.cpp - 单位根共享工具测试 (12 用例, spec §3.0.1)
//
// Phase 7B v1.6 M2 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 基准: tests/unit/timeseries/unit_root_baseline.inc (arch 8.0.0 自动生成,
//       脚本 tests/fixtures/timeseries/verify_unit_root.py)
//
// 幻觉点覆盖:
//   U1: Schwert 规则 ceil(12·(T/100)^0.25) — T=250 → 16 (基准 SCHWERT_250)
//   U5/U11: cov_nw 复刻 γ_j 分母恒为 n — L0/L2 手算基准
//   AIC lag 选择 vs arch ADF(lags=None, method='aic') 4 组精确对照
//   U2: select_trend_spec 方向性验证 (线性 → ct, RW → 非 ct)
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/timeseries/unit_root/unit_root_common.hpp"
#include "unit_root_baseline.inc"

namespace ur = cpphub::v1::timeseries::unit_root;
using cpphub::Real;
using cpphub::Size;

// 基准序列 (namespace ...::unit_root::baseline)
static std::vector<Real> make_rw() {
    return {ur::baseline::Y_RW, ur::baseline::Y_RW + ur::baseline::T};
}
static std::vector<Real> make_ar() {
    return {ur::baseline::Y_AR, ur::baseline::Y_AR + ur::baseline::T};
}

// ---------------------------------------------------------------------------
// 1. Schwert lag: T=250 → 16 (基准), 手算补充 T=100/50 (U1 向上取整)
// ---------------------------------------------------------------------------
TEST(UnitRootCommon, SchwertLagMatchesArch) {
    EXPECT_EQ(ur::schwert_lag(250), 16u);              // 基准 SCHWERT_250
    EXPECT_EQ(ur::schwert_lag(ur::baseline::T), 16u);
    // 手算: 12·(100/100)^0.25 = 12 → 12
    EXPECT_EQ(ur::schwert_lag(100), 12u);
    // 手算: 12·(50/100)^0.25 = 12·0.8409 = 10.09 → ceil = 11 (非 floor=10)
    EXPECT_EQ(ur::schwert_lag(50), 11u);
}

// ---------------------------------------------------------------------------
// 2. Schwert 边界: T=0 抛异常
// ---------------------------------------------------------------------------
TEST(UnitRootCommon, SchwertLagThrowsOnEmpty) {
    EXPECT_THROW(ur::schwert_lag(0), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 3. LRV L=0: u=[1..5] → γ₀ = (1+4+9+16+25)/5 = 11.0 (COVNW_U5_L0)
// ---------------------------------------------------------------------------
TEST(UnitRootCommon, LongRunVarianceL0HandComputed) {
    const std::vector<Real> u = {1, 2, 3, 4, 5};
    const Real lrv = ur::long_run_variance(u, 0);
    EXPECT_NEAR(lrv, ur::baseline::COVNW_U5_L0, 1e-12);
    EXPECT_DOUBLE_EQ(lrv, 11.0);
}

// ---------------------------------------------------------------------------
// 4. LRV L=2 Bartlett: 手算 25.133333333333336 (COVNW_U5_L2)
//    γ₁ = (2+6+12+20)/5 = 8, w₁ = 2/3; γ₂ = (3+8+15)/5 = 5.2, w₂ = 1/3
//    λ = 11 + 2·(2/3·8 + 1/3·5.2) = 11 + 2·7.0667 = 25.1333
// ---------------------------------------------------------------------------
TEST(UnitRootCommon, LongRunVarianceL2HandComputed) {
    const std::vector<Real> u = {1, 2, 3, 4, 5};
    const Real lrv = ur::long_run_variance(u, 2);
    EXPECT_NEAR(lrv, ur::baseline::COVNW_U5_L2, 1e-12);
}

// ---------------------------------------------------------------------------
// 5. LRV 非 Bartlett 核抛异常 (U6/U11-kernel: v1.6 仅 Bartlett)
// ---------------------------------------------------------------------------
TEST(UnitRootCommon, LongRunVarianceThrowsOnNonBartlett) {
    const std::vector<Real> u = {1, 2, 3, 4, 5};
    EXPECT_THROW(ur::long_run_variance(u, 2, "QS"),
                 std::invalid_argument);
    EXPECT_THROW(ur::long_run_variance(u, 2, "Parzen"),
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 6. LRV 空序列抛异常
// ---------------------------------------------------------------------------
TEST(UnitRootCommon, LongRunVarianceThrowsOnEmpty) {
    const std::vector<Real> u;
    EXPECT_THROW(ur::long_run_variance(u, 0), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 7. AIC lag 选择 vs arch (4 组精确对照, ADF_AIC 基准)
//    rw/n→5, rw/c→5, ar/n→0, ar/c→0 (Python 复刻预验证一致)
// ---------------------------------------------------------------------------
TEST(UnitRootCommon, SelectLagByIcAicMatchesArch) {
    const auto rw = make_rw();
    const auto ar = make_ar();
    EXPECT_EQ(ur::select_lag_by_ic(rw, "n", 0, "aic"), 5u);
    EXPECT_EQ(ur::select_lag_by_ic(rw, "c", 0, "aic"), 5u);
    EXPECT_EQ(ur::select_lag_by_ic(ar, "n", 0, "aic"), 0u);
    EXPECT_EQ(ur::select_lag_by_ic(ar, "c", 0, "aic"), 0u);
}

// ---------------------------------------------------------------------------
// 8. BIC lag 选择: 惩罚 log(n)·p > 2·p → 选择的 lag ≤ AIC (方向性)
// ---------------------------------------------------------------------------
TEST(UnitRootCommon, SelectLagByIcBicNotLargerThanAic) {
    const auto rw = make_rw();
    const auto ar = make_ar();
    EXPECT_LE(ur::select_lag_by_ic(rw, "c", 0, "bic"),
              ur::select_lag_by_ic(rw, "c", 0, "aic"));
    EXPECT_LE(ur::select_lag_by_ic(ar, "c", 0, "bic"),
              ur::select_lag_by_ic(ar, "c", 0, "aic"));
}

// ---------------------------------------------------------------------------
// 9. 非法 criterion 抛异常
// ---------------------------------------------------------------------------
TEST(UnitRootCommon, SelectLagByIcThrowsOnBadCriterion) {
    const auto rw = make_rw();
    EXPECT_THROW(ur::select_lag_by_ic(rw, "c", 0, "hqic"),
                 std::invalid_argument);
    EXPECT_THROW(ur::select_lag_by_ic(rw, "bad"), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 10. trend spec: 完美线性趋势 → "ct" (U2 方向性)
// ---------------------------------------------------------------------------
TEST(UnitRootCommon, SelectTrendSpecDetectsLinearTrend) {
    const std::vector<Real> lin = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    EXPECT_EQ(ur::select_trend_spec(lin), "ct");
}

// ---------------------------------------------------------------------------
// 11. trend spec: 零均值白噪声 → "nc"; 非零均值 → "c" (U2 方向性)
// ---------------------------------------------------------------------------
TEST(UnitRootCommon, SelectTrendSpecOnNoiseAndConstant) {
    // 均值 0 白噪声 (固定种子手工序列, 常数项不显著)
    const std::vector<Real> noise = {0.31, -0.21, 0.04, -0.38, 0.19,
                                     -0.11, 0.27, -0.30, 0.09, -0.16,
                                     0.22, -0.05, 0.13, -0.27, 0.01};
    EXPECT_EQ(ur::select_trend_spec(noise), "nc");
    // 非零均值 (常数项显著) → "c"
    const std::vector<Real> shifted = {10.31, 9.79, 10.04, 9.62, 10.19,
                                       9.89, 10.27, 9.70, 10.09, 9.84,
                                       10.22, 9.95, 10.13, 9.73, 10.01};
    EXPECT_EQ(ur::select_trend_spec(shifted), "c");
}

// ---------------------------------------------------------------------------
// 12. trend spec: 小样本抛异常
// ---------------------------------------------------------------------------
TEST(UnitRootCommon, SelectTrendSpecThrowsOnSmallSample) {
    const std::vector<Real> small = {1, 2, 3, 4};
    EXPECT_THROW(ur::select_trend_spec(small), std::invalid_argument);
}
