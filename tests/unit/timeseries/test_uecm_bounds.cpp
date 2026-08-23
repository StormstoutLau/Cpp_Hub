// =============================================================================
// test_uecm_bounds.cpp - UECM 重参数化 + PSS 边限检验 (Phase 8 M1 Step 4,
// 18 用例)
//
// 基准: ardl_baseline.inc UECM_*/BOUNDS_*/DANISH_* (statsmodels 0.14.6
// UECM.from_ardl().fit() + bounds_test 主基准 1e-10; R ARDL 0.2.5
// bounds_f/t_test 交叉 1e-8 — verify_ardl.py [S4] 三方自证:
// B 变换 ≤4.5e-15, F ≤8.5e-14, t ≤7.1e-15) +
// pss_t_bounds.inc (PSS 2001 Table CII, R ARDL 转录与原文逐位一致)
//
// 幻觉点: A4 (t 仅 Case I/III/V, statsmodels 无 t) / PSS6 (k=force 数) /
//   PSS7 (F=Wald/q) / PSS8 (UECM 样本=ARDL 同一样本) / PSS9 (F 表档
//   {10,5,1,0.1}% vs t 表档 {10,5,2.5,1}% 对偶) / PSS10 (statsmodels
//   bounds_test 查表 k 偏移 +1 — danish 回归断言) / case-trend 相容
//   (R ARDL parse_case 镜像) / ctt 无 case 映射
// =============================================================================

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/ardl/ardl_model.hpp"
#include "cpphub/timeseries/ardl/uecm_bounds.hpp"
#include "ardl_baseline.inc"

namespace ardl = cpphub::v1::timeseries::ardl;
namespace ab = cpphub::v1::ardl_baseline;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

template <Size N>
std::vector<Real> vec(const double (&a)[N]) {
    return std::vector<Real>(std::begin(a), std::end(a));
}

std::vector<std::vector<Real>> exog_k1() {
    return {vec(ab::DATA_K1_X1)};
}

std::vector<std::vector<Real>> exog_k2() {
    return {vec(ab::DATA_K2_X1), vec(ab::DATA_K2_X2)};
}

/// k1 夹具 → ARDL → UECM (trend 由 case 对应, 对拍前提)
ardl::UecmResult uecm_k1(Size p, const std::vector<Size>& q,
                         const std::string& trend) {
    ardl::ArdlSpec s;
    s.p = p;
    s.q = q;
    s.trend = trend;
    return ardl::uecm_from_ardl(
        ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s));
}

ardl::UecmResult uecm_k2(Size p, const std::vector<Size>& q,
                         const std::string& trend) {
    ardl::ArdlSpec s;
    s.p = p;
    s.q = q;
    s.trend = trend;
    return ardl::uecm_from_ardl(
        ardl::ardl_fit(vec(ab::DATA_K2_Y), exog_k2(), s));
}

/// UECM params/bse 对照 (vs statsmodels UECM.from_ardl().fit() 1e-10)
template <Size NP, Size NB>
void check_uecm(const ardl::UecmResult& u, const double (&params)[NP],
                const double (&bse)[NB], Size nobs) {
    ASSERT_EQ(u.params.size(), NP);
    ASSERT_EQ(u.vcov.size(), NP);
    for (Size i = 0; i < NP; ++i) {
        EXPECT_NEAR(u.params[i], params[i], 1e-10) << "param " << i;
        EXPECT_NEAR(std::sqrt(u.vcov[i][i]), bse[i], 1e-10) << "bse " << i;
    }
    EXPECT_EQ(u.n_obs_used, nobs);
}

/// 手工构造最小 UECM (结论分支测试; vcov = 单位阵)
ardl::UecmResult manual_uecm(Real pi_y, Real pi_1) {
    ardl::UecmResult u;
    u.source_spec.p = 2;
    u.source_spec.q = {1};
    u.source_spec.trend = "c";
    // 列序 [const, π_y, π_1, ψ_1, γ_{1,0}] (case III: kt=1)
    u.params = {0.0, pi_y, pi_1, 0.0, 0.0};
    const Size n = 5;
    u.vcov.assign(n, std::vector<Real>(n, 0.0));
    for (Size i = 0; i < n; ++i) u.vcov[i][i] = 1.0;
    u.pi = {pi_y, pi_1};
    u.pi_vcov = {{1.0, 0.0}, {0.0, 1.0}};
    u.n_obs_used = 100;
    return u;
}

}  // namespace

// 1. B 变换 k1 五 case: params/bse vs statsmodels UECM 1e-10
TEST(UecmBounds, Reparam_K1_FiveCases) {
    check_uecm(uecm_k1(2, {1}, "n"), ab::UECM_C1_K1_PARAMS,
               ab::UECM_C1_K1_BSE, 248);
    check_uecm(uecm_k1(2, {1}, "c"), ab::UECM_C2_K1_PARAMS,
               ab::UECM_C2_K1_BSE, 248);
    check_uecm(uecm_k1(2, {1}, "c"), ab::UECM_C3_K1_PARAMS,
               ab::UECM_C3_K1_BSE, 248);
    check_uecm(uecm_k1(2, {1}, "ct"), ab::UECM_C4_K1_PARAMS,
               ab::UECM_C4_K1_BSE, 248);
    check_uecm(uecm_k1(2, {1}, "ct"), ab::UECM_C5_K1_PARAMS,
               ab::UECM_C5_K1_BSE, 248);
}

// 2. B 变换 k2 (逐变量阶 q=[1,2]): params/bse vs statsmodels 1e-10
TEST(UecmBounds, Reparam_K2_TwoCases) {
    check_uecm(uecm_k2(2, {1, 2}, "c"), ab::UECM_C3_K2_PARAMS,
               ab::UECM_C3_K2_BSE, 248);
    check_uecm(uecm_k2(2, {1, 2}, "ct"), ab::UECM_C5_K2_PARAMS,
               ab::UECM_C5_K2_BSE, 248);
}

// 3. π 块提取: pi = params[kt..kt+k], pi_vcov 子块一致
TEST(UecmBounds, PiBlockExtraction) {
    const auto u = uecm_k1(2, {1}, "c");  // case III 布局
    ASSERT_EQ(u.pi.size(), 2u);           // π_y + π_1
    EXPECT_EQ(u.pi[0], u.params[1]);      // kt=1 → π_y = params[1]
    EXPECT_EQ(u.pi[1], u.params[2]);
    EXPECT_EQ(u.pi_vcov[0][0], u.vcov[1][1]);
    EXPECT_EQ(u.pi_vcov[0][1], u.vcov[1][2]);
    EXPECT_EQ(u.pi_vcov[1][1], u.vcov[2][2]);
    // σ² 与样本透传 (PSS8: UECM 样本 = ARDL 同一样本)
    const auto u2 = uecm_k2(2, {1, 2}, "c");
    EXPECT_EQ(u2.n_obs_used, 248u);  // T=250 − hb=max(2,2)=2 (PSS8)
}

// 4. F 五 case k1: stat vs statsmodels 1e-10 + R 1e-8; p 双界 1e-10
TEST(UecmBounds, F_K1_FiveCases_VsStatsmodelsAndR) {
    {
        auto r = ardl::bounds_f_test(uecm_k1(2, {1}, "n"), ardl::PssCase::I);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C1_K1_STAT[0], 1e-10);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C1_K1_STAT[1], 1e-8);
        EXPECT_NEAR(r.p_values[0], ab::BOUNDS_C1_K1_STAT[2], 1e-10);
        EXPECT_NEAR(r.p_values[1], ab::BOUNDS_C1_K1_STAT[3], 1e-10);
        EXPECT_FALSE(r.is_t_test);
        EXPECT_EQ(r.k, 1u);
    }
    {
        auto r = ardl::bounds_f_test(uecm_k1(2, {1}, "c"), ardl::PssCase::II);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C2_K1_STAT[0], 1e-10);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C2_K1_STAT[1], 1e-8);
        EXPECT_NEAR(r.p_values[0], ab::BOUNDS_C2_K1_STAT[2], 1e-10);
        EXPECT_NEAR(r.p_values[1], ab::BOUNDS_C2_K1_STAT[3], 1e-10);
    }
    {
        auto r = ardl::bounds_f_test(uecm_k1(2, {1}, "c"), ardl::PssCase::III);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C3_K1_STAT[0], 1e-10);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C3_K1_STAT[1], 1e-8);
        EXPECT_NEAR(r.p_values[0], ab::BOUNDS_C3_K1_STAT[2], 1e-10);
        EXPECT_NEAR(r.p_values[1], ab::BOUNDS_C3_K1_STAT[3], 1e-10);
    }
    {
        auto r = ardl::bounds_f_test(uecm_k1(2, {1}, "ct"), ardl::PssCase::IV);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C4_K1_STAT[0], 1e-10);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C4_K1_STAT[1], 1e-8);
    }
    {
        auto r = ardl::bounds_f_test(uecm_k1(2, {1}, "ct"), ardl::PssCase::V);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C5_K1_STAT[0], 1e-10);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C5_K1_STAT[1], 1e-8);
        EXPECT_NEAR(r.p_values[0], ab::BOUNDS_C5_K1_STAT[2], 1e-10);
        EXPECT_NEAR(r.p_values[1], ab::BOUNDS_C5_K1_STAT[3], 1e-10);
    }
}

// 5. F k2 (k_PSS=2): stat + p
TEST(UecmBounds, F_K2_TwoCases) {
    {
        auto r = ardl::bounds_f_test(uecm_k2(2, {1, 2}, "c"),
                                     ardl::PssCase::III);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C3_K2_STAT[0], 1e-10);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C3_K2_STAT[1], 1e-8);
        EXPECT_NEAR(r.p_values[0], ab::BOUNDS_C3_K2_STAT[2], 1e-10);
        EXPECT_NEAR(r.p_values[1], ab::BOUNDS_C3_K2_STAT[3], 1e-10);
        EXPECT_EQ(r.k, 2u);  // PSS6: k = force 变量数
    }
    {
        auto r = ardl::bounds_f_test(uecm_k2(2, {1, 2}, "ct"),
                                     ardl::PssCase::V);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C5_K2_STAT[0], 1e-10);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C5_K2_STAT[1], 1e-8);
    }
}

// 6. F 约束数五 case: I/III/V = k+1, II/IV = k+2 (PSS7 配套)
TEST(UecmBounds, F_RestrictionCounts) {
    EXPECT_EQ(ardl::bounds_f_test(uecm_k1(2, {1}, "n"), ardl::PssCase::I)
                  .n_restrictions, 2u);   // k+1 = 2
    EXPECT_EQ(ardl::bounds_f_test(uecm_k1(2, {1}, "c"), ardl::PssCase::II)
                  .n_restrictions, 3u);   // k+2 = 3
    EXPECT_EQ(ardl::bounds_f_test(uecm_k1(2, {1}, "c"), ardl::PssCase::III)
                  .n_restrictions, 2u);
    EXPECT_EQ(ardl::bounds_f_test(uecm_k1(2, {1}, "ct"), ardl::PssCase::IV)
                  .n_restrictions, 3u);
    EXPECT_EQ(ardl::bounds_f_test(uecm_k1(2, {1}, "ct"), ardl::PssCase::V)
                  .n_restrictions, 2u);
    EXPECT_EQ(ardl::bounds_f_test(uecm_k2(2, {1, 2}, "c"), ardl::PssCase::III)
                  .n_restrictions, 3u);   // k=2 → k+1 = 3
}

// 7. F CV: 本库键 (k_PSS, case) 直查 vs 基准 (statsmodels 表, 键直查非
//    bounds_test 返回值 — PSS10)
TEST(UecmBounds, F_CriticalValues_KeySemantics) {
    const auto check = [](const ardl::BoundsTestResult& r,
                          const double cv[]) {
        for (Size i = 0; i < 4; ++i) {
            EXPECT_NEAR(r.crit_vals[0][i], cv[i], 1e-12) << "I0 " << i;
            EXPECT_NEAR(r.crit_vals[1][i], cv[4 + i], 1e-12) << "I1 " << i;
        }
    };
    check(ardl::bounds_f_test(uecm_k1(2, {1}, "n"), ardl::PssCase::I),
          ab::BOUNDS_C1_K1_CV);
    check(ardl::bounds_f_test(uecm_k1(2, {1}, "c"), ardl::PssCase::III),
          ab::BOUNDS_C3_K1_CV);
    check(ardl::bounds_f_test(uecm_k1(2, {1}, "ct"), ardl::PssCase::V),
          ab::BOUNDS_C5_K1_CV);
    check(ardl::bounds_f_test(uecm_k2(2, {1, 2}, "c"), ardl::PssCase::III),
          ab::BOUNDS_C3_K2_CV);  // k=2 键
}

// 8. t 三 case × k1/k2: vs R ARDL 1e-8 + 构造式 (π_y/bse) 1e-10
TEST(UecmBounds, T_Case135_VsR) {
    {
        auto r = ardl::bounds_t_test(uecm_k1(2, {1}, "n"), ardl::PssCase::I);
        EXPECT_TRUE(r.is_t_test);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C1_K1_T[0], 1e-8);  // R
        EXPECT_NEAR(r.stat, ab::BOUNDS_C1_K1_T[1], 1e-10); // statsmodels 构造式
    }
    {
        auto r = ardl::bounds_t_test(uecm_k1(2, {1}, "c"), ardl::PssCase::III);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C3_K1_T[0], 1e-8);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C3_K1_T[1], 1e-10);
    }
    {
        auto r = ardl::bounds_t_test(uecm_k1(2, {1}, "ct"), ardl::PssCase::V);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C5_K1_T[0], 1e-8);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C5_K1_T[1], 1e-10);
    }
    {
        auto r = ardl::bounds_t_test(uecm_k2(2, {1, 2}, "c"), ardl::PssCase::III);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C3_K2_T[0], 1e-8);
    }
    {
        auto r = ardl::bounds_t_test(uecm_k2(2, {1, 2}, "ct"), ardl::PssCase::V);
        EXPECT_NEAR(r.stat, ab::BOUNDS_C5_K2_T[0], 1e-8);
    }
}

// 9. t Case II/IV → throw (PSS Table CII 边界, A4)
TEST(UecmBounds, T_Case24_Throw) {
    const auto u_c = uecm_k1(2, {1}, "c");
    const auto u_ct = uecm_k1(2, {1}, "ct");
    EXPECT_THROW(ardl::bounds_t_test(u_c, ardl::PssCase::II),
                 std::invalid_argument);
    EXPECT_THROW(ardl::bounds_t_test(u_ct, ardl::PssCase::IV),
                 std::invalid_argument);
}

// 10. t CV = PSS 2001 Table CII (R ARDL 转录, 与原文逐位一致);
//     ⚠️ 档位 {10,5,2.5,1}% ≠ F 表 {10,5,1,0.1}% (PSS9 对偶)
TEST(UecmBounds, T_CriticalValues_PSS2001TableCII) {
    // case I k=1 全四档 (原文 Table CII(i))
    auto r1 = ardl::bounds_t_test(uecm_k1(2, {1}, "n"), ardl::PssCase::I);
    const double c1_i0[] = {-1.62, -1.95, -2.24, -2.58};
    const double c1_i1[] = {-2.28, -2.60, -2.90, -3.22};
    for (Size i = 0; i < 4; ++i) {
        EXPECT_NEAR(r1.crit_vals[0][i], c1_i0[i], 1e-12);
        EXPECT_NEAR(r1.crit_vals[1][i], c1_i1[i], 1e-12);
    }
    // case III k=1 (原文锚 −2.86/−3.22)
    auto r3 = ardl::bounds_t_test(uecm_k1(2, {1}, "c"), ardl::PssCase::III);
    EXPECT_NEAR(r3.crit_vals[0][1], -2.86, 1e-12);
    EXPECT_NEAR(r3.crit_vals[1][1], -3.22, 1e-12);
    // case V k=2 5% (原文 Table CII(v))
    auto r5 = ardl::bounds_t_test(uecm_k2(2, {1, 2}, "ct"), ardl::PssCase::V);
    EXPECT_NEAR(r5.crit_vals[0][1], -3.41, 1e-12);
    EXPECT_NEAR(r5.crit_vals[1][1], -3.95, 1e-12);
    // case III k=2 5%
    EXPECT_NEAR(ardl::bounds_t_test(uecm_k2(2, {1, 2}, "c"),
                                    ardl::PssCase::III).crit_vals[1][1],
                -3.53, 1e-12);
}

// 11. t 方向 (负侧拒绝): 实数据 Level + 手工三态
TEST(UecmBounds, T_Direction_NegativeSide) {
    // 实数据: t = −10.22 < I1(−3.22) → LevelRelation
    auto r = ardl::bounds_t_test(uecm_k1(2, {1}, "c"), ardl::PssCase::III);
    EXPECT_LT(r.stat, r.crit_vals[1][1]);
    EXPECT_EQ(r.conclusion, ardl::BoundsConclusion::LevelRelation);
    EXPECT_FALSE(r.has_p);
    EXPECT_TRUE(std::isnan(r.p_values[0]));

    // 手工: t=0 > I0(−2.86) → NoLevelRelation
    auto rn = ardl::bounds_t_test(manual_uecm(0.0, 0.0), ardl::PssCase::III);
    EXPECT_EQ(rn.conclusion, ardl::BoundsConclusion::NoLevelRelation);
    // 手工: t=−3.0 ∈ (−3.22, −2.86) → Inconclusive
    auto ri = ardl::bounds_t_test(manual_uecm(-3.0, 0.0), ardl::PssCase::III);
    EXPECT_EQ(ri.conclusion, ardl::BoundsConclusion::Inconclusive);
}

// 12. F 结论三态 (手工 UECM, case III k=1 5% CV = 4.9229/5.7217)
TEST(UecmBounds, F_Conclusion_ThreeStates) {
    // F = (π_y² + π_1²)/2 (vcov 单位阵, 约束 [π_y, π_1])
    auto rl = ardl::bounds_f_test(manual_uecm(3.0, 3.0), ardl::PssCase::III);
    EXPECT_NEAR(rl.stat, 9.0, 1e-12);   // > 5.7217 → Level
    EXPECT_EQ(rl.conclusion, ardl::BoundsConclusion::LevelRelation);

    auto rn = ardl::bounds_f_test(manual_uecm(1.0, 1.0), ardl::PssCase::III);
    EXPECT_NEAR(rn.stat, 1.0, 1e-12);   // < 4.9229 → NoLevel
    EXPECT_EQ(rn.conclusion, ardl::BoundsConclusion::NoLevelRelation);

    auto ri = ardl::bounds_f_test(manual_uecm(2.3, 2.3), ardl::PssCase::III);
    EXPECT_NEAR(ri.stat, 5.29, 1e-12);  // ∈ (4.9229, 5.7217) → Inconclusive
    EXPECT_EQ(ri.conclusion, ardl::BoundsConclusion::Inconclusive);
}

// 13. case-trend 相容 (R ARDL parse_case 镜像): 不匹配 → throw; ctt → throw
TEST(UecmBounds, CaseTrend_Compat_Throws) {
    const auto u_n = uecm_k1(2, {1}, "n");
    const auto u_c = uecm_k1(2, {1}, "c");
    const auto u_ct = uecm_k1(2, {1}, "ct");
    EXPECT_THROW(ardl::bounds_f_test(u_c, ardl::PssCase::I),
                 std::invalid_argument);   // case I 需 "n"
    EXPECT_THROW(ardl::bounds_f_test(u_n, ardl::PssCase::III),
                 std::invalid_argument);   // case III 需 "c"
    EXPECT_THROW(ardl::bounds_f_test(u_c, ardl::PssCase::IV),
                 std::invalid_argument);   // case IV 需 "ct"
    EXPECT_THROW(ardl::bounds_f_test(u_ct, ardl::PssCase::II),
                 std::invalid_argument);   // case II 需 "c"
    // ctt: ARDL 可拟合 (Step 2 扩展), 但无 case 映射
    ardl::ArdlSpec s;
    s.p = 1;
    s.q = {1};
    s.trend = "ctt";
    const auto u_ctt = ardl::uecm_from_ardl(
        ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s));
    EXPECT_THROW(ardl::bounds_f_test(u_ctt, ardl::PssCase::V),
                 std::invalid_argument);
}

// 14. alpha 档位: F 表 {10,5,1,0.1}% vs t 表 {10,5,2.5,1}% (PSS9 对偶面)
TEST(UecmBounds, Alpha_Grids_F_vs_T) {
    const auto u = uecm_k1(2, {1}, "c");
    // F: 四档全通过
    EXPECT_NO_THROW(ardl::bounds_f_test(u, ardl::PssCase::III, 0.10));
    EXPECT_NO_THROW(ardl::bounds_f_test(u, ardl::PssCase::III, 0.05));
    EXPECT_NO_THROW(ardl::bounds_f_test(u, ardl::PssCase::III, 0.01));
    EXPECT_NO_THROW(ardl::bounds_f_test(u, ardl::PssCase::III, 0.001));
    EXPECT_THROW(ardl::bounds_f_test(u, ardl::PssCase::III, 0.025),
                 std::invalid_argument);  // F 表无 2.5% 档
    EXPECT_THROW(ardl::bounds_f_test(u, ardl::PssCase::III, 0.5),
                 std::invalid_argument);
    // t: 2.5% 有, 0.1%/0.001% 无
    EXPECT_NO_THROW(ardl::bounds_t_test(u, ardl::PssCase::III, 0.025));
    EXPECT_THROW(ardl::bounds_t_test(u, ardl::PssCase::III, 0.001),
                 std::invalid_argument);
    EXPECT_THROW(ardl::bounds_t_test(u, ardl::PssCase::III, 0.5),
                 std::invalid_argument);
}

// 15. danish (PSS 2001 UK earnings 数据): F + 样本口径
TEST(UecmBounds, Danish_F_And_Sample) {
    ardl::ArdlSpec s;
    s.p = 2;
    s.q = {1, 1, 1};
    s.trend = "c";
    const auto u = ardl::uecm_from_ardl(ardl::ardl_fit(
        vec(ab::DANISH_LRM),
        {vec(ab::DANISH_LRY), vec(ab::DANISH_IBO), vec(ab::DANISH_IDE)}, s));
    EXPECT_EQ(u.n_obs_used, static_cast<Size>(ab::DANISH_BOUNDS[1]));
    auto r = ardl::bounds_f_test(u, ardl::PssCase::III);
    // danish 宏观数据 (未标准化) 条件数差 → 绝对容差 1e-9
    // (相对 ~4e-11, 仍远优于 R 交叉的 1e-8)
    EXPECT_NEAR(r.stat, ab::DANISH_BOUNDS[0], 1e-9);
    EXPECT_EQ(r.k, 3u);  // PSS6: lry/ibo/ide = 3 force 变量
    EXPECT_EQ(r.n_restrictions, 4u);
}

// 16. PSS10 回归断言: 本库 CV = 键 (3,3) ≠ statsmodels bounds_test 返回
//     (其内部用偏移键 (4,3)) — 防未来误镜像
TEST(UecmBounds, Danish_KOffset_Regression) {
    ardl::ArdlSpec s;
    s.p = 2;
    s.q = {1, 1, 1};
    s.trend = "c";
    const auto u = ardl::uecm_from_ardl(ardl::ardl_fit(
        vec(ab::DANISH_LRM),
        {vec(ab::DANISH_LRY), vec(ab::DANISH_IBO), vec(ab::DANISH_IDE)}, s));
    auto r = ardl::bounds_f_test(u, ardl::PssCase::III);
    for (Size i = 0; i < 4; ++i) {
        EXPECT_NEAR(r.crit_vals[0][i], ab::DANISH_CV_LIB[i], 1e-12);
        EXPECT_NEAR(r.crit_vals[1][i], ab::DANISH_CV_LIB[4 + i], 1e-12);
    }
    // 5% 档分离实证: lib [3.2290, 4.3223] vs statsmodels 返回 [2.8772, 4.0105]
    EXPECT_GT(std::abs(r.crit_vals[0][1] - ab::DANISH_CV_SM[1]), 0.1);
    EXPECT_GT(std::abs(r.crit_vals[1][1] - ab::DANISH_CV_SM[5]), 0.1);
}

// 17. PSS 原表交叉印证 (F 表 ~0.05 容差 — 不同模拟批次; t 表逐位 — 转录)
TEST(UecmBounds, PSS_OriginalTable_Anchor) {
    // F: statsmodels 3200 万次模拟 vs PSS 原文 4 万次 → 容差 0.05
    const auto lo = ardl::pss_crit_values(1, 3, false);
    const auto up = ardl::pss_crit_values(1, 3, true);
    EXPECT_NEAR(lo[1], 4.94, 0.05);   // PSS Table CI(iii) k=1 5%
    EXPECT_NEAR(up[1], 5.73, 0.05);
    // t: PSS 2001 Table CII(iii) 转录 (R ARDL), 逐位精确
    const auto u = uecm_k1(2, {1}, "c");
    const auto rt = ardl::bounds_t_test(u, ardl::PssCase::III);
    EXPECT_DOUBLE_EQ(rt.crit_vals[0][1], -2.86);
    EXPECT_DOUBLE_EQ(rt.crit_vals[1][1], -3.22);
}

// 18. 校验: ARDL(0,0) / k=0 bounds 拒绝 / ardl_fit hb 下限 (statsmodels
//     L586 镜像, Step 2 补校验回归)
TEST(UecmBounds, Validation_Throws) {
    // ARDL(0, ·) k=0 全退化: 无任何滞后 → 无 UECM
    ardl::ArdlSpec s00ar;
    s00ar.p = 0;
    s00ar.q = {};
    const auto r_00ar = ardl::ardl_fit(vec(ab::DATA_K1_Y), {}, s00ar);
    EXPECT_THROW(ardl::uecm_from_ardl(r_00ar), std::invalid_argument);

    // k=0: UECM 数学上可行 (p≥1), 但 bounds 无 k=0 键 (PSS6)
    ardl::ArdlSpec s1;
    s1.p = 1;
    s1.q = {};
    const auto u_ar = ardl::uecm_from_ardl(
        ardl::ardl_fit(vec(ab::DATA_K1_Y), {}, s1));
    EXPECT_EQ(u_ar.pi.size(), 1u);  // 仅 π_y
    EXPECT_THROW(ardl::bounds_f_test(u_ar, ardl::PssCase::III),
                 std::invalid_argument);
    EXPECT_THROW(ardl::bounds_t_test(u_ar, ardl::PssCase::III),
                 std::invalid_argument);

    // p=0 且 q=0: UECM 需差分至少 1 观测 → throw
    ardl::ArdlSpec s00;
    s00.p = 0;
    s00.q = {0};
    const auto r00 = ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s00);
    EXPECT_EQ(r00.hold_back_used, 0u);
    EXPECT_THROW(ardl::uecm_from_ardl(r00), std::invalid_argument);

    // ardl_fit 显式 hb < max(p, max q) → throw (负索引防护, statsmodels L586)
    ardl::ArdlSpec s_hb;
    s_hb.p = 2;
    s_hb.q = {1};
    s_hb.hold_back = 1;  // < max(2, 1) = 2
    EXPECT_THROW(ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s_hb),
                 std::invalid_argument);
    // hb=5 (≥ 2) 仍通过
    s_hb.hold_back = 5;
    EXPECT_NO_THROW(ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s_hb));
}
