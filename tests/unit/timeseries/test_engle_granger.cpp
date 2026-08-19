// =============================================================================
// test_engle_granger.cpp - Engle-Granger 两步法 (Phase 7C M3, 14 用例)
//
// 基准: coint_baseline.inc (statsmodels 0.14.4, 1e-10)
// 幻觉点: CI1 (N=2 协整响应面) / CI2 (p 与 cv 不同源, 分列断言) /
//         CI3 (方向依赖)
// =============================================================================

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/timeseries/cointegration/engle_granger.hpp"
#include "cpphub/timeseries/cointegration/mackinnon_coint_cv.hpp"
#include "coint_baseline.inc"

namespace cb = cpphub::v1::timeseries::coint_baseline;
namespace co = cpphub::v1::timeseries::cointegration;
using cpphub::v1::Real;
using cpphub::v1::Size;

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

const double* eg_arr(const char* pair, const char* trend, const char* what) {
    // 命名拼接 EG_<PAIR>_<TREND>_<WHAT> (全大写) — 仅用于已知存在的基准
    static std::string buf;
    buf = std::string("EG_") + pair + "_" + trend + "_" + what;
    for (auto& ch : buf) ch = static_cast<char>(::toupper(ch));
    if (buf == "EG_Y1_Y2_C_T") return cb::EG_Y1_Y2_C_T;
    if (buf == "EG_Y1_Y2_C_P") return cb::EG_Y1_Y2_C_P;
    if (buf == "EG_Y1_Y2_C_CV") return cb::EG_Y1_Y2_C_CV;
    if (buf == "EG_Y1_Y2_CT_T") return cb::EG_Y1_Y2_CT_T;
    if (buf == "EG_Y1_Y2_CT_P") return cb::EG_Y1_Y2_CT_P;
    if (buf == "EG_Y1_Y2_CT_CV") return cb::EG_Y1_Y2_CT_CV;
    if (buf == "EG_Y1_Y2_CTT_T") return cb::EG_Y1_Y2_CTT_T;
    if (buf == "EG_Y1_Y2_CTT_P") return cb::EG_Y1_Y2_CTT_P;
    if (buf == "EG_Y1_Y2_CTT_CV") return cb::EG_Y1_Y2_CTT_CV;
    if (buf == "EG_Y1_Y2_N_T") return cb::EG_Y1_Y2_N_T;
    if (buf == "EG_Y1_Y2_N_P") return cb::EG_Y1_Y2_N_P;
    if (buf == "EG_Y2_Y1_C_T") return cb::EG_Y2_Y1_C_T;
    if (buf == "EG_Y2_Y1_C_P") return cb::EG_Y2_Y1_C_P;
    if (buf == "EG_Y2_Y1_C_CV") return cb::EG_Y2_Y1_C_CV;
    if (buf == "EG_Y1_Y3_C_T") return cb::EG_Y1_Y3_C_T;
    if (buf == "EG_Y1_Y3_C_P") return cb::EG_Y1_Y3_C_P;
    if (buf == "EG_Y1_Y3_C_CV") return cb::EG_Y1_Y3_C_CV;
    ADD_FAILURE() << "unknown baseline: " << buf;
    return nullptr;
}

}  // namespace

// 1. 统计量 vs statsmodels — 协整对全趋势 (CI1: 协整临界面, 非 ADF 表)
TEST(EngleGranger, StatisticAllTrends) {
    const char* trends[] = {"n", "c", "ct", "ctt"};
    for (const char* tr : trends) {
        const auto r = co::engle_granger(y1(), y2(), tr);
        const double ref = eg_arr("Y1_Y2", tr, "T")[0];
        EXPECT_NEAR(r.statistic, ref, 1e-10) << "trend=" << tr;
        EXPECT_EQ(r.trend, tr);
        EXPECT_EQ(r.n_obs, 250u);
    }
}

// 2. p 值 vs statsmodels (CI2: 1994 渐近, 与 cv 不同源)
TEST(EngleGranger, PValueVsStatsmodels) {
    const char* trends[] = {"n", "c", "ct", "ctt"};
    for (const char* tr : trends) {
        const auto r = co::engle_granger(y1(), y2(), tr);
        const double ref = eg_arr("Y1_Y2", tr, "P")[0];
        EXPECT_NEAR(r.p_value, ref, 1e-10) << "trend=" << tr;
    }
}

// 3. 临界值 vs statsmodels (2010 响应面, nobs−1 修正; c/ct/ctt)
TEST(EngleGranger, CriticalValuesVsStatsmodels) {
    const char* trends[] = {"c", "ct", "ctt"};
    for (const char* tr : trends) {
        const auto r = co::engle_granger(y1(), y2(), tr);
        const double* ref = eg_arr("Y1_Y2", tr, "CV");
        EXPECT_NEAR(r.cv_1pct, ref[0], 1e-10) << "trend=" << tr;
        EXPECT_NEAR(r.cv_5pct, ref[1], 1e-10) << "trend=" << tr;
        EXPECT_NEAR(r.cv_10pct, ref[2], 1e-10) << "trend=" << tr;
    }
}

// 4. trend="n": cv = NaN×3 (statsmodels 无 2010 n 表), p 值仍有效 (CI2)
TEST(EngleGranger, TrendNNaNcriticalValues) {
    const auto r = co::engle_granger(y1(), y2(), "n");
    EXPECT_TRUE(std::isnan(r.cv_1pct));
    EXPECT_TRUE(std::isnan(r.cv_5pct));
    EXPECT_TRUE(std::isnan(r.cv_10pct));
    EXPECT_FALSE(r.reject_null);  // NaN cv → 不拒绝
    EXPECT_NEAR(r.p_value, eg_arr("Y1_Y2", "N", "P")[0], 1e-10);
    EXPECT_NEAR(r.statistic, eg_arr("Y1_Y2", "N", "T")[0], 1e-10);
}

// 5. CI3 方向依赖: LHS 选择改变统计量; 双方向各对基准
TEST(EngleGranger, DirectionDependence) {
    const auto fwd = co::engle_granger(y1(), y2(), "c");
    const auto rev = co::engle_granger(y2(), y1(), "c");
    EXPECT_NEAR(fwd.statistic, eg_arr("Y1_Y2", "C", "T")[0], 1e-10);
    EXPECT_NEAR(rev.statistic, eg_arr("Y2_Y1", "C", "T")[0], 1e-10);
    EXPECT_NE(fwd.statistic, rev.statistic);  // CI3 核心断言
    // cv 与方向无关 (同一 N=2 表)
    EXPECT_NEAR(fwd.cv_5pct, rev.cv_5pct, 0.0);
}

// 6. 非协整对: 统计量远不及临界值, 不拒绝
TEST(EngleGranger, NoCointegrationPair) {
    const auto r = co::engle_granger(y1(), y3(), "c");
    EXPECT_NEAR(r.statistic, eg_arr("Y1_Y3", "C", "T")[0], 1e-10);
    EXPECT_NEAR(r.p_value, eg_arr("Y1_Y3", "C", "P")[0], 1e-10);
    EXPECT_FALSE(r.reject_null);
    EXPECT_GT(r.p_value, 0.30);
}

// 7. 协整对拒绝语义: reject_null == (statistic < cv_5pct)
TEST(EngleGranger, RejectSemantics) {
    for (const char* tr : {"c", "ct", "ctt"}) {
        const auto r = co::engle_granger(y1(), y2(), tr);
        EXPECT_TRUE(r.reject_null);
        EXPECT_EQ(r.reject_null, r.statistic < r.cv_5pct);
    }
}

// 8. 输入校验 (§1.4-5)
TEST(EngleGranger, InputValidation) {
    EXPECT_THROW(co::engle_granger(y1(), {1.0, 2.0}, "c"),
                 std::invalid_argument);                       // 长度不齐
    EXPECT_THROW(co::engle_granger({1.0, 2.0, 3.0}, {1.0, 2.0, 3.0}, "c"),
                 std::invalid_argument);                       // T < 6
    auto nan_y = y2();
    nan_y[5] = std::nan("");
    EXPECT_THROW(co::engle_granger(y1(), nan_y, "c"),
                 std::invalid_argument);                       // NaN
    EXPECT_THROW(co::engle_granger(y1(), y2(), "zz"),
                 std::invalid_argument);                       // trend 非法
}

// 9. MacKinnon 1994 p 值边界 (CI2 机制复刻)
TEST(EngleGranger, MackinnonPValueBoundaries) {
    // N=2, c: tau_max=0.92 → stat 1.0 ⇒ p=1; tau_min=−18.86 → stat −20 ⇒ p=0
    EXPECT_EQ(co::mackinnon_coint_p_value(1.0, "c", 2), 1.0);
    EXPECT_EQ(co::mackinnon_coint_p_value(-20.0, "c", 2), 0.0);
    // N=1, n: tau_max=inf (1e300 编码) — 任意有限 stat 不触发右尾截断
    const Real p = co::mackinnon_coint_p_value(5.0, "n", 1);
    EXPECT_GT(p, 0.0);
    EXPECT_LE(p, 1.0);
    EXPECT_THROW(co::mackinnon_coint_p_value(0.0, "c", 7),
                 std::invalid_argument);                       // N>6 表外
    EXPECT_THROW(co::mackinnon_coint_p_value(0.0, "x", 2),
                 std::invalid_argument);
}

// 10. 2010 响应面: 渐近值 = β∞ (表锚)
TEST(EngleGranger, MackinnonCVAsymptoticAnchors) {
    const auto c2 = co::mackinnon_coint_critical_values_asymptotic(2, "c");
    EXPECT_DOUBLE_EQ(c2[0], -3.89644);
    EXPECT_DOUBLE_EQ(c2[1], -3.33613);
    EXPECT_DOUBLE_EQ(c2[2], -3.04445);
    const auto ct2 = co::mackinnon_coint_critical_values_asymptotic(2, "ct");
    EXPECT_DOUBLE_EQ(ct2[0], -4.32762);
    const auto nans = co::mackinnon_coint_critical_values_asymptotic(2, "n");
    EXPECT_TRUE(std::isnan(nans[0]));
}

// 11. 2010 响应面公式: CV(T) = β∞ + β1/T + β2/T² + β3/T³ 手算复现 (c, N=2)
TEST(EngleGranger, MackinnonCVFiniteSampleFormula) {
    const Size T = 249;
    const auto cv = co::mackinnon_coint_critical_values(2, "c", T);
    const double inv = 1.0 / 249.0;
    // tau_c_2010 N=2: 1% [−3.89644, −10.9519, −33.527, 0] 等 (表锚)
    const double refs[3][4] = {{-3.89644, -10.9519, -33.527, 0.0},
                               {-3.33613, -6.1101, -6.823, 0.0},
                               {-3.04445, -4.2412, -2.720, 0.0}};
    for (int lvl = 0; lvl < 3; ++lvl) {
        const double want = refs[lvl][0] + refs[lvl][1] * inv +
                            refs[lvl][2] * inv * inv +
                            refs[lvl][3] * inv * inv * inv;
        EXPECT_NEAR(cv[lvl], want, 1e-12) << "level " << lvl;
    }
    // 与 statsmodels coint 基准一致 (同一 T=249)
    const auto r = co::engle_granger(y1(), y2(), "c");
    EXPECT_NEAR(r.cv_5pct, eg_arr("Y1_Y2", "C", "CV")[1], 1e-10);
}

// 12. 响应面参数校验
TEST(EngleGranger, MackinnonCVValidation) {
    EXPECT_THROW(co::mackinnon_coint_critical_values(13, "c", 100),
                 std::invalid_argument);
    EXPECT_THROW(co::mackinnon_coint_critical_values(2, "c", 0),
                 std::invalid_argument);
    EXPECT_THROW(co::mackinnon_coint_critical_values(2, "z", 100),
                 std::invalid_argument);
}

// 13. 共线性保护 (statsmodels: R² ≥ 1−100√ε → t = −∞, p = 0)
TEST(EngleGranger, CollinearityGuard) {
    std::vector<Real> doubled;
    for (Real v : y1()) doubled.push_back(2.0 * v);
    const auto r = co::engle_granger(doubled, y1(), "c");  // R² = 1 精确共线
    EXPECT_EQ(r.statistic, -std::numeric_limits<Real>::infinity());
    EXPECT_EQ(r.p_value, 0.0);
    EXPECT_TRUE(r.reject_null);
}

// 14. summary 回显
TEST(EngleGranger, SummaryFields) {
    const auto r = co::engle_granger(y1(), y2(), "ct");
    EXPECT_NE(r.summary.find("Engle-Granger"), std::string::npos);
    EXPECT_NE(r.summary.find("ct"), std::string::npos);
}
