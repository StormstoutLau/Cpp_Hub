// =============================================================================
// test_ardl_select.cpp - ARDL IC 网格阶数选择 (Phase 8 M1 Step 3, 10 用例)
//
// 基准: ardl_baseline.inc SEL_* (statsmodels 0.14.6 ardl_select_order,
// verify_ardl.py [S3] 生成): SEL_*_GRID = [maxlag, maxorder, hb_grid,
// n_act]; SEL_*_{IC} = [order..., ic_grid口径, ic_refit口径, nobs_refit]
// 幻觉点: A3 (同一样本: 全网格统一 hb = max(maxlag, maxorder)) /
//         网格 IC vs refit IC 双口径 (统一网格样本 T−hb_grid vs 最优模型
//         自有样本 T−max(p_best, max q_best) + A8 计账 n_acct = T−p_best) /
//         statsmodels 搜索域含变量排除候选 (maxorder+2)^k ⊃ 本库域
// =============================================================================

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/ardl/ardl_select.hpp"
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

}  // namespace

// 1. K1 BIC 全链路: order + ic_min (网格口径) + best_model.bic (refit 口径)
TEST(ArdlSelect, K1_Bic_VsStatsmodels) {
    auto r = ardl::ardl_select_order(vec(ab::DATA_K1_Y), exog_k1(), 3, 2,
                                      "bic");
    EXPECT_EQ(r.ic_used, "bic");
    EXPECT_EQ(r.p_opt, 1u);
    ASSERT_EQ(r.q_opt.size(), 1u);
    EXPECT_EQ(r.q_opt[0], 1u);
    EXPECT_NEAR(r.ic_min, ab::SEL_K1_BIC[2], 1e-10);          // 网格口径
    EXPECT_NEAR(r.best_model.bic, ab::SEL_K1_BIC[3], 1e-10);  // refit 口径
    EXPECT_EQ(r.n_obs_used, 247u);   // 统一 hb_grid = 3
    EXPECT_EQ(r.best_model.n_obs_used, 249u);  // refit hb = max(1,1) = 1
}

// 2. K1 AIC/HQIC 分支: 同数据不同 IC 入口
TEST(ArdlSelect, K1_AicHqic_VsStatsmodels) {
    auto ra = ardl::ardl_select_order(vec(ab::DATA_K1_Y), exog_k1(), 3, 2,
                                       "aic");
    EXPECT_EQ(ra.ic_used, "aic");
    EXPECT_EQ(ra.p_opt, 1u);
    EXPECT_EQ(ra.q_opt[0], 1u);
    EXPECT_NEAR(ra.ic_min, ab::SEL_K1_AIC[2], 1e-10);
    EXPECT_NEAR(ra.best_model.aic, ab::SEL_K1_AIC[3], 1e-10);

    auto rh = ardl::ardl_select_order(vec(ab::DATA_K1_Y), exog_k1(), 3, 2,
                                       "hqic");
    EXPECT_EQ(rh.ic_used, "hqic");
    EXPECT_EQ(rh.p_opt, 1u);
    EXPECT_EQ(rh.q_opt[0], 1u);
    EXPECT_NEAR(rh.ic_min, ab::SEL_K1_HQIC[2], 1e-10);
    EXPECT_NEAR(rh.best_model.hqic, ab::SEL_K1_HQIC[3], 1e-10);
}

// 3. K2 三 IC: p_best=2=hb_grid → 网格口径与 refit 口径重合 (基准印证)
TEST(ArdlSelect, K2_AllIC_VsStatsmodels) {
    for (const auto& [ic, sel] : std::vector<std::pair<const char*, const double*>>{
             {"aic", ab::SEL_K2_AIC}, {"bic", ab::SEL_K2_BIC},
             {"hqic", ab::SEL_K2_HQIC}}) {
        auto r = ardl::ardl_select_order(vec(ab::DATA_K2_Y), exog_k2(), 2, 1,
                                          ic);
        ASSERT_EQ(r.p_opt, static_cast<Size>(sel[0])) << ic;
        ASSERT_EQ(r.q_opt.size(), 2u) << ic;
        ASSERT_EQ(r.q_opt[0], static_cast<Size>(sel[1])) << ic;
        ASSERT_EQ(r.q_opt[1], static_cast<Size>(sel[2])) << ic;
        EXPECT_NEAR(r.ic_min, sel[3], 1e-10) << ic;
        // refit 口径 (基准 sel[4]; p_best=2=hb_grid → 与网格口径重合)
        if (std::string(ic) == "aic") {
            EXPECT_NEAR(r.best_model.aic, sel[4], 1e-10);
        } else if (std::string(ic) == "bic") {
            EXPECT_NEAR(r.best_model.bic, sel[4], 1e-10);
        } else {
            EXPECT_NEAR(r.best_model.hqic, sel[4], 1e-10);
        }
        EXPECT_EQ(r.n_obs_used, 248u);
        EXPECT_EQ(r.best_model.n_obs_used, 248u);
        // 网格口径 == refit 口径 (n 同 = 248, 仅浮点尾差)
        EXPECT_NEAR(r.ic_min, sel[4], 1e-9) << ic;
    }
}

// 4. 同一样本 (A3): 各候选项手动 ardl_fit(显式 hb_grid) → n_obs_used 恒等
TEST(ArdlSelect, SameSampleAcrossGrid) {
    const auto y = vec(ab::DATA_K1_Y);
    for (Size p = 0; p <= 3; ++p) {
        ardl::ArdlSpec s;
        s.p = p;
        s.q = {1};
        s.hold_back = 3;  // hb_grid = max(3, 2)
        auto r = ardl::ardl_fit(y, exog_k1(), s);
        ASSERT_EQ(r.n_obs_used, 247u) << "p=" << p;
        ASSERT_EQ(r.hold_back_used, 3u) << "p=" << p;
    }
    auto sel = ardl::ardl_select_order(y, exog_k1(), 3, 2, "bic");
    EXPECT_EQ(sel.n_obs_used, 247u);
}

// 5. 网格维度回显: k1 → 4×3=12; k2 → 3×2×2=12
TEST(ArdlSelect, GridDimensionsEcho) {
    auto r1 = ardl::ardl_select_order(vec(ab::DATA_K1_Y), exog_k1(), 3, 2);
    EXPECT_EQ(r1.n_p_grid, 4u);
    ASSERT_EQ(r1.q_grid_dims.size(), 1u);
    EXPECT_EQ(r1.q_grid_dims[0], 3u);
    EXPECT_EQ(r1.aic_grid.size(), 12u);
    EXPECT_EQ(r1.bic_grid.size(), 12u);
    EXPECT_EQ(r1.hqic_grid.size(), 12u);

    auto r2 = ardl::ardl_select_order(vec(ab::DATA_K2_Y), exog_k2(), 2, 1);
    EXPECT_EQ(r2.n_p_grid, 3u);
    ASSERT_EQ(r2.q_grid_dims.size(), 2u);
    EXPECT_EQ(r2.q_grid_dims[0], 2u);
    EXPECT_EQ(r2.q_grid_dims[1], 2u);
    EXPECT_EQ(r2.bic_grid.size(), 12u);
}

// 6. 扁平网格顺序: idx = p·4 + q0·2 + q1 (末变量最快) — 手动候选项对照
TEST(ArdlSelect, FlatGridOrdering) {
    auto r = ardl::ardl_select_order(vec(ab::DATA_K2_Y), exog_k2(), 2, 1,
                                      "bic");
    struct Pt { Size idx, p, q0, q1; };
    for (const Pt& pt : {Pt{0, 0, 0, 0}, Pt{1, 0, 0, 1}, Pt{5, 1, 0, 1},
                         Pt{11, 2, 1, 1}}) {
        ardl::ArdlSpec s;
        s.p = pt.p;
        s.q = {pt.q0, pt.q1};
        s.hold_back = 2;  // hb_grid
        auto m = ardl::ardl_fit(vec(ab::DATA_K2_Y), exog_k2(), s);
        EXPECT_NEAR(r.bic_grid[pt.idx], m.bic, 1e-10)
            << "idx=" << pt.idx;
    }
}

// 7. best_model refit 口径 A8 结构自证: n_acct = T − p_best (249),
//    llf/IC 公式重建 (不依赖基准数字 — SEL_K1_BIC[4]=249 为镜像印证)
TEST(ArdlSelect, BestModelRefitAccounting) {
    auto r = ardl::ardl_select_order(vec(ab::DATA_K1_Y), exog_k1(), 3, 2,
                                      "aic");
    const auto& bm = r.best_model;
    ASSERT_EQ(r.p_opt, 1u);
    // hb_eff = max(p, max q) = 1 → n_act = 249; A8: hb_early = p = 1,
    // n_acct = 249 (本用例两口径重合; ic_refit ≠ ic_grid 源于样本不同:
    // 网格 n=247 vs refit n=249)
    ASSERT_EQ(bm.n_obs_used, 249u);
    EXPECT_EQ(static_cast<Size>(ab::SEL_K1_AIC[4]), 249u);  // nobs_refit 镜像
    Real ssr = 0.0;
    for (Real e : bm.residuals) ssr += e * e;
    const Real n_a = 249.0;  // n_acct = T − p_best
    const Real llf = -n_a / 2.0 * (std::log(2.0 * std::acos(-1.0)) +
                                   std::log(ssr / n_a) + 1.0);
    EXPECT_NEAR(bm.loglik, llf, 1e-10);
    const Real Kp1 = static_cast<Real>(bm.params.size()) + 1.0;
    EXPECT_NEAR(bm.aic, -2.0 * llf + 2.0 * Kp1, 1e-9);
    // 网格口径 (n=247) 与 refit 口径 (n=249) 显著分离 (阻断性)
    EXPECT_GT(std::abs(r.ic_min - bm.aic), 1.0);
}

// 8. ic_min = 网格 argmin 自洽: 三 IC 逐一验证 (值 + 位置分解)
TEST(ArdlSelect, IcMinEqualsGridArgmin) {
    auto r = ardl::ardl_select_order(vec(ab::DATA_K1_Y), exog_k1(), 3, 2,
                                      "bic");
    // argmin 位置 (p 外层 4 × q 3, idx = p·3 + q)
    const auto& g = r.bic_grid;
    Size argmin = 0;
    for (Size i = 1; i < g.size(); ++i) {
        if (g[i] < g[argmin]) argmin = i;
    }
    EXPECT_EQ(argmin, r.p_opt * 3 + r.q_opt[0]);
    EXPECT_NEAR(r.ic_min, g[argmin], 0.0);
    EXPECT_NEAR(r.ic_min, *std::min_element(g.begin(), g.end()), 0.0);

    auto ra = ardl::ardl_select_order(vec(ab::DATA_K1_Y), exog_k1(), 3, 2,
                                       "aic");
    EXPECT_NEAR(ra.ic_min,
                *std::min_element(ra.aic_grid.begin(), ra.aic_grid.end()),
                0.0);
    auto rh = ardl::ardl_select_order(vec(ab::DATA_K1_Y), exog_k1(), 3, 2,
                                       "hqic");
    EXPECT_NEAR(rh.ic_min,
                *std::min_element(rh.hqic_grid.begin(), rh.hqic_grid.end()),
                0.0);
}

// 9. 校验: ic 非法 / trend 非法 / 不等长 / 样本过小 / k>10
TEST(ArdlSelect, ValidationThrows) {
    const auto y = vec(ab::DATA_K1_Y);
    const auto x = exog_k1();

    EXPECT_THROW(ardl::ardl_select_order(y, x, 3, 2, "xic"),
                 std::invalid_argument);
    EXPECT_THROW(ardl::ardl_select_order(y, x, 3, 2, "bic", "cttt"),
                 std::invalid_argument);

    auto xb = exog_k1();
    xb[0].pop_back();  // 不等长
    EXPECT_THROW(ardl::ardl_select_order(y, xb, 3, 2), std::invalid_argument);

    // 样本过小: hb_grid=200 → n_act=50 < 候选 npar (网格中途 throw)
    EXPECT_THROW(ardl::ardl_select_order(y, x, 200, 200),
                 std::invalid_argument);

    // k > 10 (ardl_fit 内校验, 首候选即 throw)
    std::vector<std::vector<Real>> x11(11, x[0]);
    EXPECT_THROW(ardl::ardl_select_order(y, x11, 1, 1),
                 std::invalid_argument);
}

// 10. 退化: maxlag=maxorder=0 → 网格 1×1 (p=0, q_j=0), hb_grid=0,
//     网格口径与 refit 口径重合 (n 均为 T; 构造式自证)
TEST(ArdlSelect, DegenerateMaxZero) {
    auto r = ardl::ardl_select_order(vec(ab::DATA_K1_Y), exog_k1(), 0, 0);
    EXPECT_EQ(r.p_opt, 0u);
    ASSERT_EQ(r.q_opt.size(), 1u);
    EXPECT_EQ(r.q_opt[0], 0u);
    EXPECT_EQ(r.n_p_grid, 1u);
    ASSERT_EQ(r.q_grid_dims.size(), 1u);
    EXPECT_EQ(r.q_grid_dims[0], 1u);
    EXPECT_EQ(r.n_obs_used, 250u);  // hb_grid = 0
    ASSERT_EQ(r.bic_grid.size(), 1u);
    // hb_grid=0 显式传 0 走默认分支 max(0,0)=0; hb_early=p=0 → n_acct=T
    // → 网格 IC 与 best_model IC 同口径同值
    EXPECT_EQ(r.best_model.n_obs_used, 250u);
    EXPECT_NEAR(r.bic_grid[0], r.best_model.bic, 1e-12);
    EXPECT_NEAR(r.aic_grid[0], r.best_model.aic, 1e-12);
    EXPECT_NEAR(r.ic_min, r.best_model.bic, 1e-12);
}
