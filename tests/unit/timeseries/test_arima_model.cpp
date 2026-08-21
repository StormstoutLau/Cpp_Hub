// =============================================================================
// test_arima_model.cpp - M1 ARIMA 估计测试 (27 用例, spec §1.2 测试矩阵
// + D-4: #16a 黄金锚 / #16b fast-vs-general 交叉验证 / #25 性能 T=1000)
//
// 基准: arima_baseline.inc
//   - R stats::arima CSS/CSS-ML (verify_arima.R 2026-08-18)
//   - forecast::Arima drift CSS-ML (漂移正解; loglik 与 statsmodels 一致)
//   - statsmodels innovations_mle (arima_statsmodels_baselines.json)
//
// 容差策略 (跨优化器对照, 7B GM 先例口径 — 先跑后按可达精度校准):
//   - CSS vs R method="CSS": 参数 ~2e-3 (SLSQP vs R optim, CSS 面较陡)
//   - CSS-ML vs R: 参数 ~5e-3 + loglik ~0.5 (ARMA(1,1) 似然面部分相消
//     平坦, R CSS-ML 落点 (0.458,−0.411) vs statsmodels (0.398,−0.358)
//     本身差 0.06 — 两成熟库的落点散布, 以 loglik 为主锚)
//   - AR2: n.cond = d + p (q 无关, arma12 p=1<q=2 实测 n.cond=1 定案)
//   - d=1: stats::arima 强制无均值 → 漂移伪根退化路径 (0.998) 也是
//     C++ include_drift=false 的合法对照; 漂移正解走 forecast
// =============================================================================
#include <gtest/gtest.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/rng.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/arima/arima_model.hpp"
#include "arima_baseline.inc"

namespace am = cpphub::v1::timeseries::arima;
namespace bl = cpphub::v1::timeseries::arima_baseline::v1;
using cpphub::Real;
using cpphub::Size;

static std::vector<Real> col(const double* a, Size n) {
    return std::vector<Real>(a, a + n);
}

// CSS 结果缓存 (每个 ~0.3s)
static const am::ArimaResult& css(Size i) {
    static std::vector<am::ArimaResult> cache(4, am::ArimaResult{});
    static bool done = false;
    if (!done) {
        const std::vector<Size> ps = {1, 2, 2, 1};
        const std::vector<Size> qs = {1, 2, 1, 2};
        const double* data[] = {bl::ARMA11, bl::ARMA22, bl::ARMA21,
                                bl::ARMA12};
        for (Size j = 0; j < 4; ++j) {
            cache[j] = am::arima_fit(col(data[j], bl::T),
                                     am::ArimaSpec{ps[j], 0, qs[j]},
                                     am::ArimaMethod::CSS);
        }
        done = true;
    }
    return cache[i];
}

static const am::ArimaResult& cssml(Size i) {
    static std::vector<am::ArimaResult> cache(4, am::ArimaResult{});
    static bool done = false;
    if (!done) {
        const std::vector<Size> ps = {1, 2, 2, 1};
        const std::vector<Size> qs = {1, 2, 1, 2};
        const double* data[] = {bl::ARMA11, bl::ARMA22, bl::ARMA21,
                                bl::ARMA12};
        for (Size j = 0; j < 4; ++j) {
            cache[j] = am::arima_fit(col(data[j], bl::T),
                                     am::ArimaSpec{ps[j], 0, qs[j]},
                                     am::ArimaMethod::CSS_ML);
        }
        done = true;
    }
    return cache[i];
}

// ---------------------------------------------------------------------------
// CSS vs R stats::arima method="CSS" (include.mean=FALSE)
// ---------------------------------------------------------------------------

// 1. arma11 CSS (φ, θ)
TEST(ArimaModel, CssArma11MatchesR) {
    const auto& r = css(0);
    ASSERT_EQ(r.params.phi.size(), 1u);
    ASSERT_EQ(r.params.theta.size(), 1u);
    EXPECT_NEAR(r.params.phi[0], bl::CSS_ARMA11[0], 2e-3);
    EXPECT_NEAR(r.params.theta[0], bl::CSS_ARMA11[1], 2e-3);
}

// 2. arma11 CSS sigma2 (SSR/(T−p) 口径)
TEST(ArimaModel, CssArma11Sigma2) {
    // R sigma2 = 0.94278 (SSR/T 与 SSR/(T−p) 差 ~0.3%, 双口径内校验)
    const auto& r = css(0);
    EXPECT_NEAR(r.params.sigma2, bl::CSS_ARMA11[2], 5e-3);
}

// 3. arma22 CSS (4 参数)
TEST(ArimaModel, CssArma22MatchesR) {
    const auto& r = css(1);
    ASSERT_EQ(r.params.phi.size(), 2u);
    ASSERT_EQ(r.params.theta.size(), 2u);
    EXPECT_NEAR(r.params.phi[0], bl::CSS_ARMA22[0], 3e-3);
    EXPECT_NEAR(r.params.phi[1], bl::CSS_ARMA22[1], 3e-3);
    EXPECT_NEAR(r.params.theta[0], bl::CSS_ARMA22[2], 3e-3);
    EXPECT_NEAR(r.params.theta[1], bl::CSS_ARMA22[3], 3e-3);
}

// 4. arma21 CSS (p≠q)
TEST(ArimaModel, CssArma21MatchesR) {
    const auto& r = css(2);
    EXPECT_NEAR(r.params.phi[0], bl::CSS_ARMA21[0], 2e-3);
    EXPECT_NEAR(r.params.phi[1], bl::CSS_ARMA21[1], 2e-3);
    EXPECT_NEAR(r.params.theta[0], bl::CSS_ARMA21[2], 2e-3);
}

// 5. AR2 定案: arma12 (p=1 < q=2) n_cond = d + p = 1 (与 q 无关)
TEST(ArimaModel, CssNCondEqualsDPlusP) {
    EXPECT_EQ(css(0).n_cond, 1u);  // arma11 p=1
    EXPECT_EQ(css(1).n_cond, 2u);  // arma22 p=2
    EXPECT_EQ(css(2).n_cond, 2u);  // arma21 p=2
    EXPECT_EQ(css(3).n_cond, 1u);  // arma12 p=1, q=2 → 1 (非 max(p,q))
}

// 6. arma12 CSS 参数 (n_cond 裁决的数据落点)
TEST(ArimaModel, CssArma12MatchesR) {
    const auto& r = css(3);
    EXPECT_NEAR(r.params.phi[0], bl::CSS_ARMA12[0], 2e-3);
    EXPECT_NEAR(r.params.theta[0], bl::CSS_ARMA12[1], 2e-3);
    EXPECT_NEAR(r.params.theta[1], bl::CSS_ARMA12[2], 2e-3);
}

// 7. CSS 残差: 前 p 个 = 0 (AR2 ε_{t<t0}=0, R 同: head 首值 0)
TEST(ArimaModel, CssResidualHeadZero) {
    const auto& r = css(0);
    ASSERT_EQ(r.residuals.size(), bl::T);
    EXPECT_NEAR(r.residuals[0], 0.0, 1e-15);
    // R resid[2] = 0.155796 (CSS 逐点对照, 1e-6 优化器落点层)
    // 注: R head 打印 0.155796409665093e+00 对应差分轴 t=1
}

// 8. CSS 收敛 + AIC 口径 (AR4: k = p+q+1+(σ²))
TEST(ArimaModel, CssConvergedAndAic) {
    const auto& r = css(0);
    EXPECT_TRUE(r.converged);
    const Real expect = -2.0 * r.loglik + 2.0 * 3.0;  // p+q=2, +σ²
    EXPECT_NEAR(r.aic, expect, 1e-9);
}

// ---------------------------------------------------------------------------
// CSS-ML vs R method="CSS-ML"
// ---------------------------------------------------------------------------

// 9. arma11 CSS-ML 参数 (似然面平坦, 落点容差放宽 + loglik 主锚)
TEST(ArimaModel, CssMlArma11MatchesR) {
    const auto& r = cssml(0);
    EXPECT_NEAR(r.params.phi[0], bl::CSSML_ARMA11[0], 5e-2);
    EXPECT_NEAR(r.params.theta[0], bl::CSSML_ARMA11[1], 5e-2);
}

// 10. arma11 CSS-ML loglik vs R (−417.73)
TEST(ArimaModel, CssMlArma11Loglik) {
    EXPECT_NEAR(cssml(0).loglik, bl::CSSML_ARMA11[3], 0.5);
}

// 11. arma22 CSS-ML
TEST(ArimaModel, CssMlArma22MatchesR) {
    const auto& r = cssml(1);
    EXPECT_NEAR(r.params.phi[0], bl::CSSML_ARMA22[0], 5e-2);
    EXPECT_NEAR(r.params.phi[1], bl::CSSML_ARMA22[1], 5e-2);
    EXPECT_NEAR(r.params.theta[0], bl::CSSML_ARMA22[2], 5e-2);
    EXPECT_NEAR(r.params.theta[1], bl::CSSML_ARMA22[3], 5e-2);
    EXPECT_NEAR(r.loglik, bl::CSSML_ARMA22[5], 0.5);
}

// 12. arma21 CSS-ML — φ/θ/σ²/loglik 全量对照 (D-4 修复后收紧)
//    v1.7 注记更正: 原注 "θ 谱等价类不唯一/似然面平坦脊" 实为
//    arma_acovf A 矩阵组装 bug (γ_bug = c(φ,θ)·γ_correct(θ*) 缩放
//    重参数化 — 集中化似然缩放不变 ⇒ loglik/φ̂ 恰好正确, θ̂/σ̂² 落在
//    重参数化像 (−0.805,0.822) vs R (−0.581,0.975)); 修复后 (θ,σ²)
//    恢复与 R 同落点, 断言从 "θ 仅方向" 收紧为数值层
TEST(ArimaModel, CssMlArma21MatchesR) {
    const auto& r = cssml(2);
    EXPECT_NEAR(r.params.phi[0], bl::CSSML_ARMA21[0], 1e-3);
    EXPECT_NEAR(r.params.phi[1], bl::CSSML_ARMA21[1], 1e-3);
    EXPECT_NEAR(r.params.theta[0], bl::CSSML_ARMA21[2], 5e-3);
    EXPECT_NEAR(r.params.sigma2, bl::CSSML_ARMA21[3], 5e-3);
    EXPECT_NEAR(r.loglik, bl::CSSML_ARMA21[4], 1e-6);
    EXPECT_LT(r.params.theta[0], 0.0);
}

// 13. arma12 CSS-ML
TEST(ArimaModel, CssMlArma12MatchesR) {
    const auto& r = cssml(3);
    EXPECT_NEAR(r.params.phi[0], bl::CSSML_ARMA12[0], 1e-2);
    EXPECT_NEAR(r.params.theta[0], bl::CSSML_ARMA12[1], 1e-2);
    EXPECT_NEAR(r.params.theta[1], bl::CSSML_ARMA12[2], 1e-2);
    EXPECT_NEAR(r.loglik, bl::CSSML_ARMA12[4], 0.5);
}

// 14. CSS-ML 精化单调 (ML 面上): CSS-ML ll ≥ CSS 解处的精确 ll
//     (AR3 注: CSS loglik 是统一高斯型口径, 与精确似然不可比 —
//      2026-08-18 实测修正: 单调性断言须在同一 ML 面上)
TEST(ArimaModel, CssMlRefinesCss) {
    const double* data[] = {bl::ARMA11, bl::ARMA22, bl::ARMA21, bl::ARMA12};
    const std::vector<Size> ps = {1, 2, 2, 1};
    const std::vector<Size> qs = {1, 2, 1, 2};
    for (Size i = 0; i < 4; ++i) {
        const std::vector<Real> z(data[i], data[i] + bl::T);
        const auto& rc = css(i);
        // CSS 解处的精确 loglik (与 CSS-ML 同面: z 无 demean, 无 drift)
        const auto pc = am::detail::arma_innovations_pieces(
            z, rc.params.phi, rc.params.theta);
        constexpr Real kTwoPi = 6.283185307179586476925286766559;
        const Real Tn = static_cast<Real>(bl::T);
        const Real s2 = pc.s_uv / Tn;
        const Real ll_css_exact = -0.5 * (Tn * std::log(kTwoPi * s2) + Tn
                                          + pc.sum_logv);
        EXPECT_GE(cssml(i).loglik, ll_css_exact - 1e-6);
    }
}

// ---------------------------------------------------------------------------
// Innovations vs statsmodels (arima_fit 路由)
// ---------------------------------------------------------------------------

// 15. arma11 Innovations vs statsmodels
TEST(ArimaModel, InnovArma11MatchesStatsmodels) {
    const auto r = am::arima_fit(col(bl::ARMA11, bl::T),
                                 am::ArimaSpec{1, 0, 1},
                                 am::ArimaMethod::Innovations);
    EXPECT_NEAR(r.params.phi[0], bl::SM_ARMA11[0], 5e-2);
    EXPECT_NEAR(r.params.theta[0], bl::SM_ARMA11[1], 5e-2);
    EXPECT_NEAR(r.loglik, bl::SM_ARMA11[3], 0.5);
}

// 16. Innovations arma21 — φ/θ/σ²/loglik 全量对照 (D-4 修复后收紧,
//     v1.7 "θ 谱等价" 误诊更正见 12 号注释)
TEST(ArimaModel, InnovArma21MatchesStatsmodels) {
    const auto r = am::arima_fit(col(bl::ARMA21, bl::T),
                                 am::ArimaSpec{2, 0, 1},
                                 am::ArimaMethod::Innovations);
    EXPECT_NEAR(r.params.phi[0], bl::SM_ARMA21[0], 1e-3);
    EXPECT_NEAR(r.params.phi[1], bl::SM_ARMA21[1], 1e-3);
    EXPECT_NEAR(r.params.theta[0], bl::SM_ARMA21[2], 5e-3);
    EXPECT_NEAR(r.params.sigma2, bl::SM_ARMA21[3], 5e-3);
    EXPECT_NEAR(r.loglik, bl::SM_ARMA21[4], 1e-6);
    EXPECT_LT(r.params.theta[0], 0.0);
}

// 16a. 通用 innovations 黄金锚 (B&D 2016 ITSF §5.2 MA(1) θ=−0.9,
// test_stattools.py L1315-1321): acov=[1.81,−0.9] → θ 行 [−0.4972,
// −0.6606, −0.7404], v=[1.81, 1.3625, 1.2155, 1.1436] — 通用版自 v1.7 起
// 为交叉验证参照 (D-4), 此锚锁死其正确性; fast 路径 v 须逐位同值
TEST(ArimaModel, InnovationsGoldenAnchorBdMa1) {
    const std::vector<Real> acov = {1.81, -0.9};
    std::vector<std::vector<Real>> th;
    std::vector<Real> v;
    am::detail::innovations_algo(acov, 4, th, v);
    ASSERT_EQ(v.size(), 4u);
    EXPECT_NEAR(v[0], 1.81, 5e-5);
    EXPECT_NEAR(v[1], 1.3625, 5e-5);
    EXPECT_NEAR(v[2], 1.2155, 5e-5);
    EXPECT_NEAR(v[3], 1.1436, 5e-5);
    ASSERT_EQ(th.size(), 4u);
    EXPECT_NEAR(th[1][0], -0.4972, 5e-5);
    EXPECT_NEAR(th[2][0], -0.6606, 5e-5);
    EXPECT_NEAR(th[3][0], -0.7404, 5e-5);
    // fast 路径 (MA(1) θ=−0.9, σ²=1 同模型): v 须与通用版逐位一致
    const std::vector<Real> z4 = {0.3, -1.2, 0.8, 0.1};
    const auto uv = am::detail::arma_innovations_uv(
        z4, std::vector<Real>{}, std::vector<Real>{-0.9});
    ASSERT_EQ(uv.second.size(), 4u);
    for (Size i = 0; i < 4; ++i) {
        EXPECT_NEAR(uv.second[i], v[i], 1e-12) << "i=" << i;
    }
}

// 16b. 快速带状 vs 通用 innovations 交叉验证 (D-4 防快但错)
// 两条路径数学等价 (带状精确性: Cov(W_t,y_s)=0 for s≤t−q−1 ⇒ θ 列≥q
// 恒零); 差异仅浮点消元次序 — v 逐位级 (探针 T=150 实测 ≤9e-16),
// u 经滤波反馈缓增 (≤6e-8), 按实测×100+ 安全余量定容差
TEST(ArimaModel, FastVsGeneralCrossCheck) {
    const Size T = 120;
    std::vector<Real> z(T);
    {
        cpphub::v1::Philox4x64 rng(2026, 23);
        for (Size t = 0; t < T; t += 2) {
            const auto w = cpphub::v1::box_muller(
                (rng() >> 11) * (1.0 / 9007199254740992.0),
                (rng() >> 11) * (1.0 / 9007199254740992.0));
            z[t] = w.first;
            if (t + 1 < T) z[t + 1] = w.second;
        }
    }
    // (p,q) 扫描: 白噪声/纯 AR/纯 MA/混合 (含 p>q 与 q>p 双向)
    const std::vector<std::vector<Real>> phis = {
        {}, {0.5}, {}, {0.5}, {0.5, -0.3}, {0.5}, {0.5, -0.3}};
    const std::vector<std::vector<Real>> thetas = {
        {}, {}, {0.4}, {0.4}, {0.4}, {0.4, 0.2}, {0.4, 0.2}};
    const std::vector<const char*> names = {
        "(0,0)", "(1,0)", "(0,1)", "(1,1)", "(2,1)", "(1,2)", "(2,2)"};
    for (Size c = 0; c < phis.size(); ++c) {
        // 通用参照 (O(T³), 全长 acov) — arma_acovf A 矩阵修复后
        // p≥2 组合恢复有效 ARMA γ, 通用算法不再抛 v 非正
        std::vector<Real> ar(1, 1.0);
        for (Real p : phis[c]) ar.push_back(-p);
        std::vector<Real> ma(1, 1.0);
        for (Real t : thetas[c]) ma.push_back(t);
        const auto gamma = am::detail::arma_acovf(ar, ma, 1.0, T);
        std::vector<std::vector<Real>> thg;
        std::vector<Real> vg;
        am::detail::innovations_algo(gamma, T, thg, vg);
        const auto ug = am::detail::innovations_filter(z, thg);
        // 快速路径 (O(T·(p+q)²))
        const auto uv = am::detail::arma_innovations_uv(z, phis[c], thetas[c]);
        const auto& uf = uv.first;
        const auto& vf = uv.second;
        Real max_du = 0.0, max_dv = 0.0, d_suv = 0.0, d_logv = 0.0;
        for (Size t = 0; t < T; ++t) {
            max_du = std::max(max_du, std::abs(uf[t] - ug[t]));
            max_dv = std::max(max_dv, std::abs(vf[t] - vg[t]));
            d_suv += std::abs(uf[t] * uf[t] / vf[t]
                              - ug[t] * ug[t] / vg[t]);
            d_logv += std::abs(std::log(vf[t]) - std::log(vg[t]));
        }
        EXPECT_LT(max_dv, 1e-12) << names[c] << " max_dv=" << max_dv;
        EXPECT_LT(max_du, 1e-5) << names[c] << " max_du=" << max_du;
        EXPECT_LT(d_suv, 1e-3) << names[c] << " d_suv=" << d_suv;
        EXPECT_LT(d_logv, 1e-10) << names[c] << " d_logv=" << d_logv;
    }
}

// ---------------------------------------------------------------------------
// d=1: 退化路径 + 漂移正解
// ---------------------------------------------------------------------------

// 17. d=1 无漂移: stats::arima 退化路径 (漂移未建模 → AR 伪根 0.998)
TEST(ArimaModel, D1NoDriftDegeneratePath) {
    const auto r = am::arima_fit(col(bl::ARIMA111D_LEVEL, bl::T),
                                 am::ArimaSpec{1, 1, 1, false},
                                 am::ArimaMethod::CSS);
    ASSERT_EQ(r.n_cond, 2u);  // d + p
    EXPECT_EQ(r.n_obs_used, bl::T - 1);
    // 退化路径: φ 接近 1 (漂移被 AR 吸收), 落点层容差
    EXPECT_NEAR(r.params.phi[0], bl::CSS_ARIMA111D[0], 2e-2);
    EXPECT_NEAR(r.params.theta[0], bl::CSS_ARIMA111D[1], 5e-2);
}

// 18. d=1 + drift: forecast::Arima CSS-ML 漂移正解 (loglik 同 statsmodels)
TEST(ArimaModel, D1DriftMatchesForecastArima) {
    const auto r = am::arima_fit(col(bl::ARIMA111D_LEVEL, bl::T),
                                 am::ArimaSpec{1, 1, 1, true},
                                 am::ArimaMethod::CSS_ML);
    EXPECT_NEAR(r.params.phi[0], bl::DRIFT_CSSML[0], 1e-2);
    EXPECT_NEAR(r.params.theta[0], bl::DRIFT_CSSML[1], 1e-2);
    EXPECT_NEAR(r.params.drift, bl::DRIFT_CSSML[2], 1e-2);
    EXPECT_NEAR(r.loglik, bl::DRIFT_CSSML[4], 0.5);
}

// 19. drift 恢复 DGP 真值 0.3 (夹具生成参数)
TEST(ArimaModel, D1DriftRecoversDgp) {
    const auto r = am::arima_fit(col(bl::ARIMA111D_LEVEL, bl::T),
                                 am::ArimaSpec{1, 1, 1, true},
                                 am::ArimaMethod::CSS_ML);
    // DGP: Δy = 0.3 + ARMA(1,1); drift 估计 0.357 (含小样本噪声, 真值附近)
    EXPECT_GT(r.params.drift, 0.25);
    EXPECT_LT(r.params.drift, 0.45);
}

// 20. AIC 口径 AR4 (drift 版): npar = p+q+1(drift)+1(σ²) = 4
//     (forecast::Arima 同口径: aic = 862.606 实测吻合)
TEST(ArimaModel, AicDriftParamCount) {
    const auto r = am::arima_fit(col(bl::ARIMA111D_LEVEL, bl::T),
                                 am::ArimaSpec{1, 1, 1, true},
                                 am::ArimaMethod::CSS_ML);
    const Real expect = -2.0 * r.loglik + 2.0 * 4.0;
    EXPECT_NEAR(r.aic, expect, 1e-9);
}

// ---------------------------------------------------------------------------
// 边界与配置 (AR7/AR8/AR5)
// ---------------------------------------------------------------------------

// 21. drift + d=0 拒绝 (AR5)
TEST(ArimaModel, DriftD0Rejected) {
    EXPECT_THROW(am::arima_fit(col(bl::ARMA11, bl::T),
                               am::ArimaSpec{1, 0, 1, true},
                               am::ArimaMethod::CSS),
                 std::invalid_argument);
}

// 22. NaN 输入拒绝
TEST(ArimaModel, NaNRejected) {
    auto d = col(bl::ARMA11, bl::T);
    d[10] = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(am::arima_fit(d, am::ArimaSpec{1, 0, 1}, am::ArimaMethod::CSS),
                 std::invalid_argument);
}

// 23. 样本过小拒绝
TEST(ArimaModel, SampleTooSmallRejected) {
    const std::vector<Real> d(6, 0.1);
    EXPECT_THROW(am::arima_fit(d, am::ArimaSpec{2, 0, 2}, am::ArimaMethod::CSS),
                 std::invalid_argument);
}

// 24. 多起始配置 (AR7): 关 HR 后仍收敛 (零起始可用)
TEST(ArimaModel, MultiStartConfig) {
    am::ArimaConfig cfg;
    cfg.use_hannan_rissanen = false;
    cfg.n_starting_points = 2;
    const auto r = am::arima_fit(col(bl::ARMA11, bl::T),
                                 am::ArimaSpec{1, 0, 1},
                                 am::ArimaMethod::CSS, cfg);
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.params.phi[0], bl::CSS_ARMA11[0], 2e-3);
}

// 25. 性能 P4 (checklist §15.4, D-2 裁决 2026-08-20) — T=1000 CSS-ML 多起始 < 10s
//
// 7B 先例模式 + §2.3 对策: n_starting_points=2 (成本近减半仍覆盖多起始
// 路径, #24 先例) + 收敛/参数落点断言防"快但错"。DGP: ARMA(2,1)
// φ=(0.5,−0.3), θ=0.4, σ=1, burn-in 200 (Philox 固定 seed)。
// 成本结构 (D-4 修正): 似然评估 O(T·(p+q)²) 带状快速算法 (探针证伪了
// v1.2 报告 "SLSQP 主导" 的预估 — 通用 O(T³) 算法实测 158.9s 击穿预算
// 16×, 移植后 0.04s); 预算 < 10s = 238× 安全余量 (CI 2-3× 慢径无忧)。
// ---------------------------------------------------------------------------
TEST(ArimaModel, PerformanceT1000CssMl) {
    const Size T = 1000, BURN = 200;
    const Real phi1 = 0.5, phi2 = -0.3, theta1 = 0.4;
    std::vector<Real> y;
    y.reserve(T + BURN);
    {
        cpphub::v1::Philox4x64 rng(2026, 22);
        Real y1 = 0.0, y2 = 0.0, e_prev = 0.0;
        for (Size t = 0; t < T + BURN; ++t) {
            const auto z = cpphub::v1::box_muller(
                (rng() >> 11) * (1.0 / 9007199254740992.0),
                (rng() >> 11) * (1.0 / 9007199254740992.0));
            const Real e = z.first;
            const Real yt = phi1 * y1 + phi2 * y2 + e + theta1 * e_prev;
            y2 = y1;
            y1 = yt;
            e_prev = e;
            y.push_back(yt);
        }
        y.erase(y.begin(), y.begin() + BURN);  // burn-in 丢弃
    }
    am::ArimaConfig cfg;
    cfg.use_hannan_rissanen = false;
    cfg.n_starting_points = 2;
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = am::arima_fit(y, am::ArimaSpec{2, 0, 1},
                                 am::ArimaMethod::CSS_ML, cfg);
    const auto elapsed = std::chrono::duration_cast<std::chrono::duration<Real>>(
        std::chrono::steady_clock::now() - t0);
    printf("[perf] ARIMA T=1000 CSS-ML (2 starts): %.3fs\n", elapsed.count());
    // 收敛 + 参数落点 (防"快但错", 跨优化器宽松容差)
    EXPECT_TRUE(r.converged) << r.message;
    ASSERT_EQ(r.params.phi.size(), 2u);
    ASSERT_EQ(r.params.theta.size(), 1u);
    EXPECT_NEAR(r.params.phi[0], phi1, 0.15);
    EXPECT_NEAR(r.params.phi[1], phi2, 0.15);
    EXPECT_NEAR(r.params.theta[0], theta1, 0.15);
    EXPECT_TRUE(std::isfinite(r.loglik));
    // 预算: < 10s (spec §15.4; 主控站实测见 printf, CI 2-3× 慢径余量)
    EXPECT_LT(elapsed.count(), 10.0) << "elapsed " << elapsed.count() << "s";
}
