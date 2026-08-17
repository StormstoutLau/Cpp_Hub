// =============================================================================
// test_zivot_andrews.cpp - Zivot-Andrews 检验测试 (15 用例, spec §1.2/§9.5)
//
// 基准: tests/unit/timeseries/za_baseline.inc (verify_za.py 自动生成:
//       statsmodels 0.14.6 Baum 模式 + verify_za.R urca 1.3-4 固定 lag 转录)
//
// 容差 (spec §1.3):
//   - vs statsmodels (baum 模式): stat 1e-10 / p 1e-12 / lag,bp 精确
//   - vs urca (固定 lag=1, trim 放开): stat 1e-8 (实测逐位 ~1e-12)
//   - 论文临界值: EXPECT_DOUBLE_EQ (精确)
//
// 幻觉点覆盖 (spec §9.5):
//   ZA1 双模式 lag (用例 15) / ZA2 trim 参数化 (9-11) / ZA3 DU/DT 构造经
//   双库基准隐含 + 断点定位 (7) / ZA4 双临界值表 (13-14) / ZA5 min 取向 (8)
// =============================================================================
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"
#include "cpphub/timeseries/unit_root/zivot_andrews_test.hpp"
#include "za_baseline.inc"

namespace ur = cpphub::v1::timeseries::unit_root;
using ur::ZAModel;
using ur::zivot_andrews_test;
using cpphub::Real;
using cpphub::Size;

static std::vector<Real> y_za() {
    return {ur::za_baseline::Y_ZA, ur::za_baseline::Y_ZA + ur::za_baseline::T};
}

// ---------------------------------------------------------------------------
// 1-3: statsmodels 基准 (Baum 模式: baum_preselect=true, trim=0.15) — 1e-10/1e-12
// ---------------------------------------------------------------------------
TEST(ZivotAndrewsTest, SmBaselineModelA) {
    const auto r = zivot_andrews_test(y_za(), ZAModel::A, 0, true, 0.15);
    EXPECT_NEAR(r.statistic, ur::za_baseline::SM_A.stat, 1e-10);
    EXPECT_NEAR(r.p_value_baum, ur::za_baseline::SM_A.p, 1e-12);
    EXPECT_EQ(r.break_index, ur::za_baseline::SM_A.bpidx);
    EXPECT_EQ(r.n_lags, ur::za_baseline::SM_A.lag);
}

TEST(ZivotAndrewsTest, SmBaselineModelB) {
    const auto r = zivot_andrews_test(y_za(), ZAModel::B, 0, true, 0.15);
    EXPECT_NEAR(r.statistic, ur::za_baseline::SM_B.stat, 1e-10);
    EXPECT_NEAR(r.p_value_baum, ur::za_baseline::SM_B.p, 1e-12);
    EXPECT_EQ(r.n_lags, ur::za_baseline::SM_B.lag);
    // ZA1 已知对照边界: statsmodels "t" 模型 DT 用 cutoff−1 边界 ⇒ bpidx 差 1;
    // C++ 取 urca 形 (主对照), break_index = SM_B.bpidx − 1 (stat/p 已 1e-14 对齐)
    EXPECT_EQ(r.break_index, ur::za_baseline::SM_B.bpidx - 1);
}

TEST(ZivotAndrewsTest, SmBaselineModelC) {
    const auto r = zivot_andrews_test(y_za(), ZAModel::C, 0, true, 0.15);
    EXPECT_NEAR(r.statistic, ur::za_baseline::SM_C.stat, 1e-10);
    EXPECT_NEAR(r.p_value_baum, ur::za_baseline::SM_C.p, 1e-12);
    EXPECT_EQ(r.break_index, ur::za_baseline::SM_C.bpidx);
    EXPECT_EQ(r.n_lags, ur::za_baseline::SM_C.lag);
}

// ---------------------------------------------------------------------------
// 4-6: urca 基准 (主模式: fixed_lag=1, trim=0.001 放开网格) — 1e-8
// ---------------------------------------------------------------------------
TEST(ZivotAndrewsTest, UrcaFixedLagIntercept) {
    const auto r = zivot_andrews_test(y_za(), ZAModel::A, 1, false, 0.001);
    EXPECT_NEAR(r.statistic, ur::za_baseline::URCA_A_STAT, 1e-8);
    EXPECT_EQ(r.break_index + 1, ur::za_baseline::URCA_A_BPOINT);  // 1-based bpoint
    EXPECT_EQ(r.n_lags, 1);
}

TEST(ZivotAndrewsTest, UrcaFixedLagTrend) {
    const auto r = zivot_andrews_test(y_za(), ZAModel::B, 1, false, 0.001);
    EXPECT_NEAR(r.statistic, ur::za_baseline::URCA_B_STAT, 1e-8);
    EXPECT_EQ(r.break_index + 1, ur::za_baseline::URCA_B_BPOINT);
}

TEST(ZivotAndrewsTest, UrcaFixedLagBoth) {
    const auto r = zivot_andrews_test(y_za(), ZAModel::C, 1, false, 0.001);
    EXPECT_NEAR(r.statistic, ur::za_baseline::URCA_C_STAT, 1e-8);
    EXPECT_EQ(r.break_index + 1, ur::za_baseline::URCA_C_BPOINT);
}

// ---------------------------------------------------------------------------
// 7: 断点搜索定位 — 崩溃均值构造断点 0-based 70, Model A/C 应在 ±5 内 (ZA3 方向)
// ---------------------------------------------------------------------------
TEST(ZivotAndrewsTest, BreakSearchNearTrueBreak) {
    const auto ra = zivot_andrews_test(y_za(), ZAModel::A, 0, true, 0.15);
    const auto rc = zivot_andrews_test(y_za(), ZAModel::C, 0, true, 0.15);
    // Size 无符号, 用 double 绝对差防下溢
    EXPECT_NEAR(static_cast<double>(ra.break_index),
                static_cast<double>(ur::za_baseline::TRUE_BREAK), 5.0);
    EXPECT_NEAR(static_cast<double>(rc.break_index),
                static_cast<double>(ur::za_baseline::TRUE_BREAK), 5.0);
}

// ---------------------------------------------------------------------------
// 8: ZA5 — 统计量 = min_{Tb} t(α̂) (最负), path 全有限
// ---------------------------------------------------------------------------
TEST(ZivotAndrewsTest, StatisticIsMinOfPath) {
    const auto r = zivot_andrews_test(y_za(), ZAModel::C, 0, true, 0.15);
    ASSERT_FALSE(r.t_stats_path.empty());
    for (Real v : r.t_stats_path) {
        EXPECT_TRUE(std::isfinite(v));
    }
    EXPECT_NEAR(r.statistic, *std::min_element(r.t_stats_path.begin(),
                                               r.t_stats_path.end()),
                1e-15);
}

// ---------------------------------------------------------------------------
// 9: ZA2 — trim 网格参数化 (0.15 → 84 候选; 0.25 → 60)
// ---------------------------------------------------------------------------
TEST(ZivotAndrewsTest, TrimGridSize) {
    const auto r15 = zivot_andrews_test(y_za(), ZAModel::C, 0, true, 0.15);
    EXPECT_EQ(r15.t_stats_path.size(), ur::za_baseline::GRID_015);
    EXPECT_DOUBLE_EQ(r15.trim, 0.15);
    const auto r25 = zivot_andrews_test(y_za(), ZAModel::C, 0, true, 0.25);
    // trimcnt=30 ⇒ b in [31, 90] ⇒ 60 候选
    EXPECT_EQ(r25.t_stats_path.size(), 60);
    EXPECT_DOUBLE_EQ(r25.trim, 0.25);
}

// ---------------------------------------------------------------------------
// 10: ZA2 — trim 越界 (≤0 或 >1/3) 抛 invalid_argument
// ---------------------------------------------------------------------------
TEST(ZivotAndrewsTest, TrimBoundsThrow) {
    EXPECT_THROW(zivot_andrews_test(y_za(), ZAModel::C, 0, true, 0.0),
                 std::invalid_argument);
    EXPECT_THROW(zivot_andrews_test(y_za(), ZAModel::C, 0, true, -0.1),
                 std::invalid_argument);
    EXPECT_THROW(zivot_andrews_test(y_za(), ZAModel::C, 0, true, 0.34),
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 11: 网格空 (大 lag 内收 b_lo 越过 b_hi) 抛 invalid_argument
// ---------------------------------------------------------------------------
TEST(ZivotAndrewsTest, EmptyGridThrows) {
    std::vector<Real> y30(30);
    for (Size i = 0; i < 30; ++i) {
        y30[i] = static_cast<Real>(i);  // 非常数
    }
    // T=30, k=20, trim=0.3: b_lo = max(10, 22) = 22 > b_hi = min(21, 29) = 21
    EXPECT_THROW(zivot_andrews_test(y30, ZAModel::C, 20, false, 0.3),
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 12: 异常输入 (T<10 / NaN / 零方差)
// ---------------------------------------------------------------------------
TEST(ZivotAndrewsTest, InvalidInputsThrow) {
    EXPECT_THROW(zivot_andrews_test({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0},
                                    ZAModel::C),
                 std::invalid_argument);
    auto y_nan = y_za();
    y_nan[5] = std::nan("");
    EXPECT_THROW(zivot_andrews_test(y_nan, ZAModel::C), std::invalid_argument);
    EXPECT_THROW(zivot_andrews_test(std::vector<Real>(50, 3.14), ZAModel::C),
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 13: ZA4 — ZA1992 论文表临界值 (主) 精确相等 + reject 语义
// ---------------------------------------------------------------------------
TEST(ZivotAndrewsTest, PaperCriticalValues) {
    const auto ra = zivot_andrews_test(y_za(), ZAModel::A, 0, true, 0.15);
    EXPECT_DOUBLE_EQ(ra.critical_1pct, -5.34);
    EXPECT_DOUBLE_EQ(ra.critical_5pct, -4.80);
    EXPECT_DOUBLE_EQ(ra.critical_10pct, -4.58);
    const auto rb = zivot_andrews_test(y_za(), ZAModel::B, 0, true, 0.15);
    EXPECT_DOUBLE_EQ(rb.critical_1pct, -4.93);
    EXPECT_DOUBLE_EQ(rb.critical_5pct, -4.42);
    EXPECT_DOUBLE_EQ(rb.critical_10pct, -4.11);
    const auto rc = zivot_andrews_test(y_za(), ZAModel::C, 0, true, 0.15);
    EXPECT_DOUBLE_EQ(rc.critical_1pct, -5.57);
    EXPECT_DOUBLE_EQ(rc.critical_5pct, -5.08);
    EXPECT_DOUBLE_EQ(rc.critical_10pct, -4.82);
    // H0: 带单断点单位根; stat > cv5 ⇒ 不拒绝 (左尾)
    EXPECT_FALSE(rc.reject_null);
    EXPECT_EQ(rc.reject_null, rc.statistic < rc.critical_5pct);
}

// ---------------------------------------------------------------------------
// 14: ZA4 — MC 表 p 值插值 (statsmodels L2529 语义: 节点命中/中点/端点 clamp)
// ---------------------------------------------------------------------------
TEST(ZivotAndrewsTest, McPvalueInterpolation) {
    // 表节点精确命中 (Model A: 1%/5%/10%)
    EXPECT_NEAR(ur::detail::za_mc_pvalue(-5.27644, ZAModel::A), 0.01, 1e-15);
    EXPECT_NEAR(ur::detail::za_mc_pvalue(-4.81067, ZAModel::A), 0.05, 1e-15);
    EXPECT_NEAR(ur::detail::za_mc_pvalue(-4.56618, ZAModel::A), 0.10, 1e-15);
    // 线性插值中点 (0.9% 与 1.0% 节点间: -5.30294/-5.27644)
    EXPECT_NEAR(ur::detail::za_mc_pvalue((-5.30294 + -5.27644) / 2.0, ZAModel::A),
                0.0095, 1e-12);
    // 端点 clamp
    EXPECT_DOUBLE_EQ(ur::detail::za_mc_pvalue(-99.0, ZAModel::A), 0.00001);
    EXPECT_DOUBLE_EQ(ur::detail::za_mc_pvalue(1.0, ZAModel::A), 0.999);
    // B/C 1% 锚 + ZA4 trap: MC c 1% = −5.27644 (非论文 −5.34, 非 0.1% 分位 −5.83192)
    EXPECT_NEAR(ur::detail::za_mc_pvalue(-5.03421, ZAModel::B), 0.01, 1e-15);
    EXPECT_NEAR(ur::detail::za_mc_pvalue(-5.57556, ZAModel::C), 0.01, 1e-15);
    EXPECT_DOUBLE_EQ(ur::detail::ZA_MC_C[10].cv, -5.27644);
}

// ---------------------------------------------------------------------------
// 15: ZA1 — 双模式 lag 语义 (baum 预选 vs Schwert 自动 vs 用户固定) + 默认参数
// ---------------------------------------------------------------------------
TEST(ZivotAndrewsTest, DualModeLagSemantics) {
    // Baum 对照模式: AIC 一次性预选 (该数据选 0, 与 statsmodels 一致)
    const auto rb = zivot_andrews_test(y_za(), ZAModel::C, 0, true, 0.15);
    EXPECT_EQ(rb.n_lags, 0);
    // 主模式默认 (fixed_lag=0, baum=false): Schwert 自动 k=13 ≠ Baum 的 0
    const auto rd = zivot_andrews_test(y_za(), ZAModel::C);  // 全默认
    EXPECT_EQ(rd.n_lags, ur::schwert_lag(ur::za_baseline::T));
    EXPECT_EQ(rd.n_lags, 13);
    // 两模式统计量显著不同 (lag 不同 ⇒ 不同回归)
    EXPECT_GT(std::fabs(rd.statistic - rb.statistic), 1e-3);
    // 用户固定 lag 回显
    const auto rf = zivot_andrews_test(y_za(), ZAModel::C, 3, false, 0.15);
    EXPECT_EQ(rf.n_lags, 3);
    // 默认参数回显 (spec §2.2: model=C, trim=0.15)
    EXPECT_EQ(rd.model, ZAModel::C);
    EXPECT_DOUBLE_EQ(rd.trim, 0.15);
    EXPECT_NE(rd.summary.find("Zivot-Andrews"), std::string::npos);
}
