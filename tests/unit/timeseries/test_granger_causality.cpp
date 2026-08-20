// =============================================================================
// test_granger_causality.cpp - 线性 Granger 因果检验 (Phase 7C M1, 16 用例)
//
// 基准: granger_baseline.inc (statsmodels 0.14.4 grangercausalitytests +
//   OLS f_test 增广 Wald + 显式 NW 三明治, gen_phase7c_granger_baseline.py)
//   GR_*   10 值: [f, f_p, params_f, params_f_p, chi2, chi2_p, lr, lr_p,
//                  df_denom, df_num]                 容差 1e-10
//   GRTY_*  6 值: [wald, chi2_p, sm_f, sm_f_p, k, df_denom]  容差 1e-8/1e-10
//   GRHAC_* 6 值: [wald, chi2_p, sm_f, sm_f_p, p, df_denom]  容差 1e-8/1e-10
//
// 数据: X→Y 真因果 (DGP 一阶) / Z 独立 / (Y1,Y2) I(1) 协整对水平值
//
// 幻觉点 (spec §9.1):
//   GR1 (df 公式)      — df2 = nobs−(2p+1), nobs = T−p; 四统计量内部一致性
//   GR2 (TY df=k)      — 错误自由度 k+d_max 无法复现基准 p 值
//   GR3 (d_max 外部)   — d_max 不干涉标准路径; NaN 政策
//   GR4 (增广阶不进约束) — p2d2 基准: 估计 p+d 阶, 约束仅前 k 阶
//   GR5 (非稳健默认)   — HAC 独立数值路径, 与标准 χ² 版差异断言
//   GR6 (方向)         — 显式 (cause,effect) 形参; fwd/rev 各对基准
//   GR7 (I(1) 失效)    — I(1) 水平值统计量入基准 (集成场景 2 做对比)
//
// HAC 带宽语义: 公开 API hac_bandwidth=0 → NW 规则默认
//   floor(4·(249/100)^(2/9)) = 4 (GRHAC_FWD_P1_L4 锁定);
//   White 退化 (L=0) 仅 detail::hac_wald_statistic 可达 (GRHAC_FWD_P1_L0)
// =============================================================================

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cpphub/econometrics/core/special_functions.hpp"
#include "cpphub/timeseries/arima/granger_test.hpp"
#include "granger_baseline.inc"

namespace gb = cpphub::v1::timeseries::granger_baseline;
namespace gr = cpphub::v1::timeseries::granger;
namespace ed = cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

std::vector<Real> vx() {
    return {std::begin(gb::X), std::end(gb::X)};
}
std::vector<Real> vy() {
    return {std::begin(gb::Y), std::end(gb::Y)};
}
std::vector<Real> vz() {
    return {std::begin(gb::Z), std::end(gb::Z)};
}
std::vector<Real> vy1() {
    return {std::begin(gb::Y1), std::end(gb::Y1)};
}
std::vector<Real> vy2() {
    return {std::begin(gb::Y2), std::end(gb::Y2)};
}

/// 标准四统计量 + df 全字段断言 (GR_* 10 值布局, GR1/GR6)
void expect_standard(const gr::GrangerResult& r, const double* ref) {
    EXPECT_NEAR(r.f_stat, ref[0], 1e-10);
    EXPECT_NEAR(r.f_p, ref[1], 1e-10);
    EXPECT_NEAR(r.params_f_stat, ref[2], 1e-10);
    EXPECT_NEAR(r.params_f_p, ref[3], 1e-10);
    EXPECT_NEAR(r.chi2_stat, ref[4], 1e-10);
    EXPECT_NEAR(r.chi2_p, ref[5], 1e-10);
    EXPECT_NEAR(r.lr_stat, ref[6], 1e-10);
    EXPECT_NEAR(r.lr_p, ref[7], 1e-10);
    EXPECT_EQ(r.df2, static_cast<Size>(ref[8]));
    EXPECT_EQ(r.df1, static_cast<Size>(ref[9]));
}

/// TY/HAC Wald + p 断言 (6 值布局: [wald, chi2_p, ...])
void expect_wald(Real stat, Real p, const double* ref) {
    EXPECT_NEAR(stat, ref[0], 1e-8);
    EXPECT_NEAR(p, ref[1], 1e-10);
}

}  // namespace

// 1. 标准四统计量: x→y 真因果 (p=1,2,4) — statsmodels 1e-10 (GR1 df 公式)
TEST(GrangerCausality, StandardStatsTrueCausality) {
    const auto r1 = gr::granger_test(vx(), vy(), 1, 0, false);
    expect_standard(r1, gb::GR_FWD_P1);
    const auto r2 = gr::granger_test(vx(), vy(), 2, 0, false);
    expect_standard(r2, gb::GR_FWD_P2);
    const auto r4 = gr::granger_test(vx(), vy(), 4, 0, false);
    expect_standard(r4, gb::GR_FWD_P4);
    // DGP: y_t = 0.3 + 0.5 y_{t−1} + 0.4 x_{t−1} + e ⇒ 一阶真因果显著
    EXPECT_LT(r1.f_p, 0.01);
}

// 2. 标准四统计量: 反向 y→x (GR6 方向; 反向无因果)
TEST(GrangerCausality, StandardStatsReverseDirection) {
    const auto r1 = gr::granger_test(vy(), vx(), 1, 0, false);
    expect_standard(r1, gb::GR_REV_P1);
    const auto r2 = gr::granger_test(vy(), vx(), 2, 0, false);
    expect_standard(r2, gb::GR_REV_P2);
    EXPECT_GT(r1.f_p, 0.05);  // 反向不显著
}

// 3. 标准四统计量: 独立 z→y 零假设
TEST(GrangerCausality, StandardStatsIndependentNull) {
    const auto r1 = gr::granger_test(vz(), vy(), 1, 0, false);
    expect_standard(r1, gb::GR_NULL_P1);
    const auto r2 = gr::granger_test(vz(), vy(), 2, 0, false);
    expect_standard(r2, gb::GR_NULL_P2);
    EXPECT_GT(r1.f_p, 0.05);
    EXPECT_GT(r2.f_p, 0.05);
}

// 4. 标准四统计量: I(1) 协整对水平值 (GR7 场景数据; 对比留给集成场景 2)
TEST(GrangerCausality, StandardStatsI1Levels) {
    const auto r1 = gr::granger_test(vy2(), vy1(), 1, 0, false);
    expect_standard(r1, gb::GR_I1_P1);
    const auto r2 = gr::granger_test(vy2(), vy1(), 2, 0, false);
    expect_standard(r2, gb::GR_I1_P2);
}

// 5. GR6 方向: 显式 (cause,effect) 形参; fwd ≠ rev 各对自身基准;
//    statsmodels 输入 [effect,cause] 的方向陷阱由形参顺序消灭
TEST(GrangerCausality, DirectionExplicitCauseEffect) {
    const auto fwd = gr::granger_test(vx(), vy(), 1, 0, false);  // x→y
    const auto rev = gr::granger_test(vy(), vx(), 1, 0, false);  // y→x
    EXPECT_NEAR(fwd.f_stat, gb::GR_FWD_P1[0], 1e-10);
    EXPECT_NEAR(rev.f_stat, gb::GR_REV_P1[0], 1e-10);
    EXPECT_NE(fwd.f_stat, rev.f_stat);  // 方向核心断言: 交换 ≠ 对称
    // params_ftest 与 ssr_ftest 数学等价、数值独立路径 (两列分录)
    EXPECT_NEAR(fwd.params_f_stat, fwd.f_stat, 1e-9);
    EXPECT_NEAR(rev.params_f_stat, rev.f_stat, 1e-9);
}

// 6. GR1 df 公式 + 四统计量内部一致性:
//    χ² = F·p·nobs/df2 (同源 SSR 代数恒等式); LR ≤ χ² (对数不等式)
TEST(GrangerCausality, DfFormulaInternalConsistency) {
    struct Row {
        const double* ref;
        Size p;
    };
    const Row rows[] = {{gb::GR_FWD_P1, 1}, {gb::GR_FWD_P2, 2},
                        {gb::GR_FWD_P4, 4}};
    for (const auto& row : rows) {
        const auto r = gr::granger_test(vx(), vy(), row.p, 0, false);
        EXPECT_EQ(r.df1, row.p);
        EXPECT_EQ(r.df2, static_cast<Size>(row.ref[8]));
        const Size nobs = 250 - row.p;  // T − p
        // χ² = nobs·dSSR/SSR_u = F·p·nobs/df2
        EXPECT_NEAR(r.chi2_stat,
                    r.f_stat * static_cast<Real>(row.p) *
                        static_cast<Real>(nobs) / static_cast<Real>(r.df2),
                    1e-8)
            << "p=" << row.p;
        // ln(1+x) ≤ x ⇒ LR = nobs·ln(SSR_r/SSR_u) ≤ χ²
        EXPECT_LE(r.lr_stat, r.chi2_stat) << "p=" << row.p;
    }
}

// 7. TY 增广 Wald: 真因果 (p1d1/p2d1/p2d2);
//    GR4: p2d2 估计 4 阶但约束仅前 k=2 阶 cause 滞后 (基准锁定)
TEST(GrangerCausality, TyWaldTrueCausality) {
    const auto r1 = gr::granger_test(vx(), vy(), 1, 1, false);
    EXPECT_TRUE(r1.has_ty);
    expect_wald(r1.ty_wald_stat, r1.ty_wald_p, gb::GRTY_FWD_P1_D1);
    const auto r2 = gr::granger_test(vx(), vy(), 2, 1, false);
    expect_wald(r2.ty_wald_stat, r2.ty_wald_p, gb::GRTY_FWD_P2_D1);
    const auto r3 = gr::granger_test(vx(), vy(), 2, 2, false);
    expect_wald(r3.ty_wald_stat, r3.ty_wald_p, gb::GRTY_FWD_P2_D2);
    // 增广阶改变估计样本 ⇒ wald 随 d_max 变化 (非平凡差异)
    EXPECT_NE(r2.ty_wald_stat, r3.ty_wald_stat);
    EXPECT_LT(r1.ty_wald_p, 0.01);  // 真因果 TY 亦显著 (水平值合法推断)
}

// 8. GR2: TY Wald 自由度 = k (非 k+d_max):
//    基准 p 值 = χ²(k) 上尾; 错误自由度 χ²(k+d) 差 ~1e-6/1e-7 ≫ 容差
TEST(GrangerCausality, TyWaldDfEqualsK) {
    const auto r1 = gr::granger_test(vx(), vy(), 1, 1, false);
    EXPECT_NEAR(r1.ty_wald_p, gb::GRTY_FWD_P1_D1[1], 1e-10);
    const Real wrong1 = ed::detail::chi2_sf(2.0, r1.ty_wald_stat);  // k+d=2
    EXPECT_GT(std::abs(wrong1 - gb::GRTY_FWD_P1_D1[1]), 1e-9);

    const auto r2 = gr::granger_test(vx(), vy(), 2, 1, false);
    EXPECT_NEAR(r2.ty_wald_p, gb::GRTY_FWD_P2_D1[1], 1e-10);
    const Real wrong2 = ed::detail::chi2_sf(3.0, r2.ty_wald_stat);  // k+d=3
    EXPECT_GT(std::abs(wrong2 - gb::GRTY_FWD_P2_D1[1]), 1e-9);
}

// 9. TY: 反向 + I(1) 对 (I(1) 水平值 TY 为 GR7 合法推断)
TEST(GrangerCausality, TyWaldReverseAndI1) {
    const auto rev = gr::granger_test(vy(), vx(), 1, 1, false);
    expect_wald(rev.ty_wald_stat, rev.ty_wald_p, gb::GRTY_REV_P1_D1);
    const auto i1a = gr::granger_test(vy2(), vy1(), 1, 1, false);
    expect_wald(i1a.ty_wald_stat, i1a.ty_wald_p, gb::GRTY_I1_P1_D1);
    const auto i1b = gr::granger_test(vy2(), vy1(), 2, 1, false);
    expect_wald(i1b.ty_wald_stat, i1b.ty_wald_p, gb::GRTY_I1_P2_D1);
}

// 10. NaN 政策 (§1.4-5): d_max=0 ⇒ TY NaN; with_hac=false ⇒ HAC NaN;
//     标准字段恒有效
TEST(GrangerCausality, NanPolicy) {
    const auto r = gr::granger_test(vx(), vy(), 1, 0, false);
    EXPECT_FALSE(r.has_ty);
    EXPECT_TRUE(std::isnan(r.ty_wald_stat));
    EXPECT_TRUE(std::isnan(r.ty_wald_p));
    EXPECT_FALSE(r.has_hac);
    EXPECT_TRUE(std::isnan(r.hac_wald_stat));
    EXPECT_TRUE(std::isnan(r.hac_wald_p));
    EXPECT_NEAR(r.f_stat, gb::GR_FWD_P1[0], 1e-10);  // 标准字段不受影响
    EXPECT_FALSE(r.summary.empty());
}

// 11. HAC-Wald 显式带宽 L=3: fwd p1/p2 + rev + null + i1 五场景
TEST(GrangerCausality, HacWaldExplicitBandwidth) {
    const auto r1 = gr::granger_test(vx(), vy(), 1, 0, true, 3);
    EXPECT_TRUE(r1.has_hac);
    expect_wald(r1.hac_wald_stat, r1.hac_wald_p, gb::GRHAC_FWD_P1_L3);
    const auto r2 = gr::granger_test(vx(), vy(), 2, 0, true, 3);
    expect_wald(r2.hac_wald_stat, r2.hac_wald_p, gb::GRHAC_FWD_P2_L3);
    const auto r3 = gr::granger_test(vy(), vx(), 1, 0, true, 3);
    expect_wald(r3.hac_wald_stat, r3.hac_wald_p, gb::GRHAC_REV_P1_L3);
    const auto r4 = gr::granger_test(vz(), vy(), 1, 0, true, 3);
    expect_wald(r4.hac_wald_stat, r4.hac_wald_p, gb::GRHAC_NULL_P1_L3);
    const auto r5 = gr::granger_test(vy2(), vy1(), 1, 0, true, 3);
    expect_wald(r5.hac_wald_stat, r5.hac_wald_p, gb::GRHAC_I1_P1_L3);
}

// 12. HAC 默认带宽 = NW 经验法则 floor(4·(T/100)^(2/9)):
//     nobs=249/248 → 4; GRHAC_FWD_P1_L4 锁定 hac_bandwidth=0 的默认行为
TEST(GrangerCausality, HacWaldDefaultBandwidthRule) {
    EXPECT_EQ(gr::detail::nw_default_bandwidth(249), 4u);
    EXPECT_EQ(gr::detail::nw_default_bandwidth(248), 4u);
    // 默认 (hac_bandwidth=0) 解析到 L=4 → 与显式 4 同一基准
    const auto rdef = gr::granger_test(vx(), vy(), 1, 0, true, 0);
    expect_wald(rdef.hac_wald_stat, rdef.hac_wald_p, gb::GRHAC_FWD_P1_L4);
    const auto rexpl = gr::granger_test(vx(), vy(), 1, 0, true, 4);
    EXPECT_NEAR(rexpl.hac_wald_stat, gb::GRHAC_FWD_P1_L4[0], 1e-8);
    EXPECT_NEAR(rexpl.hac_wald_stat, rdef.hac_wald_stat, 1e-12);
    // 默认 ≠ L=3 (证明默认确实解析为 4, 而非其他带宽)
    EXPECT_GT(std::abs(rdef.hac_wald_stat - gb::GRHAC_FWD_P1_L3[0]), 1e-3);
}

// 13. HAC White 退化 (L=0: Ω = Σ xu xu', 无自协方差项):
//     公开 API 的 0 语义是 NW 默认, White 仅 detail::hac_wald_statistic 可达
TEST(GrangerCausality, HacWaldWhiteDegenerate) {
    const Real w0 = gr::detail::hac_wald_statistic(vx(), vy(), 1, 0);
    EXPECT_NEAR(w0, gb::GRHAC_FWD_P1_L0[0], 1e-8);
    EXPECT_NEAR(ed::detail::chi2_sf(1.0, w0), gb::GRHAC_FWD_P1_L0[1],
                1e-10);
    // 带宽语义真实改变 Ω: White ≠ NW(L=3) ≠ NW(L=4)
    EXPECT_GT(std::abs(w0 - gb::GRHAC_FWD_P1_L3[0]), 1e-3);
    EXPECT_GT(std::abs(w0 - gb::GRHAC_FWD_P1_L4[0]), 1e-3);
    EXPECT_THROW(gr::detail::hac_wald_statistic(vx(), vy(), 1, 300),
                 std::invalid_argument);
}

// 14. GR5: 标准四统计量均非异方差稳健; HAC 为独立稳健数值路径:
//     (a) HAC-Wald 与标准 χ² 版差异显著 (b) HAC 开关不污染标准统计量
TEST(GrangerCausality, HacVsStandardDifference) {
    const auto std_only = gr::granger_test(vx(), vy(), 1, 0, false);
    const auto with_h = gr::granger_test(vx(), vy(), 1, 0, true, 3);
    EXPECT_NEAR(with_h.hac_wald_stat, gb::GRHAC_FWD_P1_L3[0], 1e-8);
    // 稳健修正产生实质性差异 (35.28 vs 38.45)
    EXPECT_GT(std::abs(with_h.hac_wald_stat - std_only.chi2_stat), 1.0);
    // 标准统计量与 HAC 开关完全无关 (同一路径)
    EXPECT_NEAR(std_only.f_stat, with_h.f_stat, 0.0);
    EXPECT_NEAR(std_only.chi2_stat, with_h.chi2_stat, 0.0);
    EXPECT_NEAR(std_only.lr_stat, with_h.lr_stat, 0.0);
    EXPECT_NEAR(std_only.params_f_stat, with_h.params_f_stat, 0.0);
}

// 15. GR3 API 契约: d_max 纯外部给定, 不做模型内单位根自适应 —
//     标准 F/χ²/LR 与 d_max 无关 (逐位一致); 仅 TY 分支随 d_max 变化
TEST(GrangerCausality, ApiContractDmaxExternal) {
    const auto r0 = gr::granger_test(vx(), vy(), 2, 0, false);
    const auto r1 = gr::granger_test(vx(), vy(), 2, 1, false);
    const auto r2 = gr::granger_test(vx(), vy(), 2, 2, false);
    EXPECT_NEAR(r0.f_stat, r1.f_stat, 0.0);
    EXPECT_NEAR(r0.chi2_stat, r2.chi2_stat, 0.0);
    EXPECT_NEAR(r0.lr_stat, r1.lr_stat, 0.0);
    EXPECT_EQ(r0.df1, 2u);
    EXPECT_EQ(r1.df1, 2u);  // lag 不因 d_max 改变 (无自适应)
    EXPECT_NEAR(r1.ty_wald_stat, gb::GRTY_FWD_P2_D1[0], 1e-8);
    EXPECT_NEAR(r2.ty_wald_stat, gb::GRTY_FWD_P2_D2[0], 1e-8);
}

// 16. 输入校验: 长度不齐/空序列/lag=0/NaN/样本不足 (标准+TY)/带宽越界
TEST(GrangerCausality, InputValidation) {
    const std::vector<Real> a(50, 1.0);
    const std::vector<Real> b(49, 1.0);
    EXPECT_THROW(gr::granger_test(a, b, 1), std::invalid_argument);
    EXPECT_THROW(gr::granger_test({}, {}), std::invalid_argument);
    EXPECT_THROW(gr::granger_test(a, a, 0), std::invalid_argument);
    std::vector<Real> bad(50, 1.0);
    bad[7] = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(gr::granger_test(bad, a, 1), std::invalid_argument);
    EXPECT_THROW(gr::granger_test(a, bad, 1), std::invalid_argument);
    // 标准样本不足: T=7, p=2 ⇒ nobs=5 ≤ k_u=5
    std::vector<Real> tiny(7);
    for (Size i = 0; i < tiny.size(); ++i) {
        tiny[i] = std::sin(static_cast<Real>(i));
    }
    EXPECT_THROW(gr::granger_test(tiny, tiny, 2), std::invalid_argument);
    // TY 样本不足: T=10, p=2, d=1 ⇒ n_ty=7 ≤ k_ty=7 (标准路径先通过)
    std::vector<Real> e10(10), c10(10);
    for (Size i = 0; i < 10; ++i) {
        e10[i] = std::cos(0.7 * static_cast<Real>(i));
        c10[i] = std::sin(0.9 * static_cast<Real>(i));
    }
    EXPECT_THROW(gr::granger_test(c10, e10, 2, 1, false),
                 std::invalid_argument);
    // 带宽越界: nobs=249, L=300
    EXPECT_THROW(gr::granger_test(vx(), vy(), 1, 0, true, 300),
                 std::invalid_argument);
}
