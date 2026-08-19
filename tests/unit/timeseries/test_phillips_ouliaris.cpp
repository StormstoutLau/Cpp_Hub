// =============================================================================
// test_phillips_ouliaris.cpp - Phillips-Ouliaris 协整检验 (Phase 7C M3, 10 用例)
//
// 基准: coint_baseline.inc (urca 1.3-4 ca.po, 1e-8)
// 幻觉点: CI12 (Pu 残差基方向依赖 / Pz 方向无关, Pz 优先)
// =============================================================================

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/timeseries/cointegration/phillips_ouliaris.hpp"
#include "coint_baseline.inc"

namespace cb = cpphub::v1::timeseries::coint_baseline;
namespace co = cpphub::v1::timeseries::cointegration;
using cpphub::v1::Real;

namespace {

std::vector<Real> y1() {
    return {std::begin(cb::Y1), std::end(cb::Y1)};
}
std::vector<Real> y2() {
    return {std::begin(cb::Y2), std::end(cb::Y2)};
}
std::vector<Real> y3() {
    return {std::begin(cb::Y3), std::end(cb::Y3)};
}

}  // namespace

// 1. Pu vs urca (协整对, demean × short 带宽), 1e-8
TEST(PhillipsOuliaris, PuVsUrca) {
    const auto r0 = co::phillips_ouliaris(y1(), y2(), "Pu", "none", 0);
    EXPECT_NEAR(r0.statistic, cb::PO_Y1Y2_PU_NONE_SHORT[0], 1e-8);
    const auto r1 = co::phillips_ouliaris(y1(), y2(), "Pu", "constant", 0);
    EXPECT_NEAR(r1.statistic, cb::PO_Y1Y2_PU_CONSTANT_SHORT[0], 1e-8);
    const auto r2 = co::phillips_ouliaris(y1(), y2(), "Pu", "trend", 0);
    EXPECT_NEAR(r2.statistic, cb::PO_Y1Y2_PU_TREND_SHORT[0], 1e-8);
    EXPECT_EQ(r0.type, "Pu");
    EXPECT_EQ(r0.n_obs, 249u);
}

// 2. Pz vs urca (协整对, demean × short), 1e-8
TEST(PhillipsOuliaris, PzVsUrca) {
    const auto r0 = co::phillips_ouliaris(y1(), y2(), "Pz", "none", 0);
    EXPECT_NEAR(r0.statistic, cb::PO_Y1Y2_PZ_NONE_SHORT[0], 1e-8);
    const auto r1 = co::phillips_ouliaris(y1(), y2(), "Pz", "constant", 0);
    EXPECT_NEAR(r1.statistic, cb::PO_Y1Y2_PZ_CONSTANT_SHORT[0], 1e-8);
    const auto r2 = co::phillips_ouliaris(y1(), y2(), "Pz", "trend", 0);
    EXPECT_NEAR(r2.statistic, cb::PO_Y1Y2_PZ_TREND_SHORT[0], 1e-8);
    EXPECT_EQ(r2.type, "Pz");
}

// 3. long 带宽 (数值 lag=15 ↔ urca "long"), Pu/Pz, 1e-8
TEST(PhillipsOuliaris, LongBandwidth) {
    const auto pu = co::phillips_ouliaris(y1(), y2(), "Pu", "none", 15);
    EXPECT_NEAR(pu.statistic, cb::PO_Y1Y2_PU_NONE_LONG[0], 1e-8);
    EXPECT_EQ(pu.lag, 15u);
    const auto pz = co::phillips_ouliaris(y1(), y2(), "Pz", "constant", 15);
    EXPECT_NEAR(pz.statistic, cb::PO_Y1Y2_PZ_CONSTANT_LONG[0], 1e-8);
}

// 4. CI12: Pz 方向无关 (swap 两参统计量不变; 双向均对 urca 基准)
//    容差 1e-9: 数学上精确不变, 数值上受求和顺序影响 — urca 自身两个方向
//    基准相差 2.45e-12 (227.29907868302834 vs 227.29907868303079), 1e-12 不可达
TEST(PhillipsOuliaris, PzDirectionInvariance) {
    const auto fwd = co::phillips_ouliaris(y1(), y2(), "Pz", "none", 0);
    const auto rev = co::phillips_ouliaris(y2(), y1(), "Pz", "none", 0);
    EXPECT_NEAR(fwd.statistic, rev.statistic, 1e-9);  // CI12 核心断言
    EXPECT_NEAR(rev.statistic, cb::PO_Y2Y1_PZ_NONE_SHORT[0], 1e-8);
}

// 5. CI12: Pu 方向依赖 (对调列统计量改变; 双向均对 urca 基准)
TEST(PhillipsOuliaris, PuDirectionDependence) {
    const auto fwd = co::phillips_ouliaris(y1(), y2(), "Pu", "none", 0);
    const auto rev = co::phillips_ouliaris(y2(), y1(), "Pu", "none", 0);
    EXPECT_NEAR(fwd.statistic, cb::PO_Y1Y2_PU_NONE_SHORT[0], 1e-8);
    EXPECT_NEAR(rev.statistic, cb::PO_Y2Y1_PU_NONE_SHORT[0], 1e-8);
    EXPECT_NE(fwd.statistic, rev.statistic);  // CI12 核心断言
}

// 6. short 带宽公式: lmax = trunc(4·(nobs/100)^0.25) = 5 (nobs=249)
TEST(PhillipsOuliaris, ShortBandwidthFormula) {
    const auto r = co::phillips_ouliaris(y1(), y2(), "Pz", "none", 0);
    EXPECT_EQ(r.lag, 5u);
    EXPECT_EQ(r.lag,
              static_cast<cpphub::v1::Size>(
                  std::trunc(4.0 * std::pow(249.0 / 100.0, 0.25))));
}

// 7. CV 表锚 vs urca cval (10/5/1%; 全 demean × 双类型)
TEST(PhillipsOuliaris, CriticalValueAnchors) {
    struct Row { const char* tag; const double* ref; };
    const Row rows[] = {
        {"PO_Y1Y2_PU_NONE_SHORT", cb::PO_Y1Y2_PU_NONE_SHORT_CVAL},
        {"PO_Y1Y2_PU_CONSTANT_SHORT", cb::PO_Y1Y2_PU_CONSTANT_SHORT_CVAL},
        {"PO_Y1Y2_PU_TREND_SHORT", cb::PO_Y1Y2_PU_TREND_SHORT_CVAL},
        {"PO_Y1Y2_PZ_NONE_SHORT", cb::PO_Y1Y2_PZ_NONE_SHORT_CVAL},
        {"PO_Y1Y2_PZ_CONSTANT_SHORT", cb::PO_Y1Y2_PZ_CONSTANT_SHORT_CVAL},
        {"PO_Y1Y2_PZ_TREND_SHORT", cb::PO_Y1Y2_PZ_TREND_SHORT_CVAL},
    };
    const char* dets[] = {"none", "constant", "trend"};
    for (int i = 0; i < 6; ++i) {
        const auto r = co::phillips_ouliaris(
            y1(), y2(), i < 3 ? "Pu" : "Pz", dets[i % 3], 0);
        EXPECT_DOUBLE_EQ(r.cv_5pct, rows[i].ref[1]) << rows[i].tag;
    }
    // 表锚 (转录保真, phillips_ouliaris.hpp static_assert 同源)
    EXPECT_DOUBLE_EQ(cb::PO_Y1Y2_PU_NONE_SHORT_CVAL[0], 20.3933);
    EXPECT_DOUBLE_EQ(cb::PO_Y1Y2_PU_NONE_SHORT_CVAL[2], 38.3413);
    EXPECT_DOUBLE_EQ(cb::PO_Y1Y2_PZ_NONE_SHORT_CVAL[1], 40.8217);
}

// 8. 拒绝语义: 协整对拒绝 (stat > cv5), 非协整对不拒绝
TEST(PhillipsOuliaris, RejectSemantics) {
    const auto coint = co::phillips_ouliaris(y1(), y2(), "Pz", "none", 0);
    EXPECT_TRUE(coint.reject_null);
    EXPECT_EQ(coint.reject_null, coint.statistic > coint.cv_5pct);
    const auto nc = co::phillips_ouliaris(y1(), y3(), "Pz", "none", 0);
    EXPECT_NEAR(nc.statistic, cb::PO_Y1Y3_PZ_NONE_SHORT[0], 1e-8);
    EXPECT_EQ(nc.reject_null, nc.statistic > nc.cv_5pct);
}

// 9. 非协整对 vs urca (Pu/Pz), 1e-8
TEST(PhillipsOuliaris, NoCointegrationPair) {
    const auto pu = co::phillips_ouliaris(y1(), y3(), "Pu", "constant", 0);
    EXPECT_NEAR(pu.statistic, cb::PO_Y1Y3_PU_CONSTANT_SHORT[0], 1e-8);
    EXPECT_FALSE(pu.reject_null);
    const auto pz = co::phillips_ouliaris(y1(), y3(), "Pz", "trend", 0);
    EXPECT_NEAR(pz.statistic, cb::PO_Y1Y3_PZ_TREND_SHORT[0], 1e-8);
}

// 10. 输入校验
TEST(PhillipsOuliaris, InputValidation) {
    EXPECT_THROW(co::phillips_ouliaris(y1(), {1.0, 2.0}),
                 std::invalid_argument);  // 长度不齐
    auto nan_y = y2();
    nan_y[3] = std::nan("");
    EXPECT_THROW(co::phillips_ouliaris(y1(), nan_y),
                 std::invalid_argument);  // NaN
    EXPECT_THROW(co::phillips_ouliaris(y1(), y2(), "zz"),
                 std::invalid_argument);  // type 非法
    EXPECT_THROW(co::phillips_ouliaris(y1(), y2(), "Pz", "zz"),
                 std::invalid_argument);  // demean 非法
    EXPECT_THROW(co::phillips_ouliaris({1.0, 2.0}, {1.0, 2.0}),
                 std::invalid_argument);  // T < 5
}
