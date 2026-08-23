// =============================================================================
// test_integration_phase8.cpp - Phase 8 端到端集成测试 (spec §5, 6 场景)
//
// 场景矩阵 (spec §5):
//   1. ARDL 完整链:    select → fit → UECM → bounds F+t (协整 DGP β 已知,
//                      k=2 双 force 变量 → LevelRelation)
//   2. ARIMA 三路径:   CSS-ML / innovations / statespace 同数据同口径
//                      (预去均值统一零均值) loglik 落点 + 参数一致
//   3. 单位根预检:     ADF 链 (7B 复用) 分类 I(0)/I(1)/I(2) → 混合
//                      I(0)/I(1) force 变量 bounds; I(2) 前置拒绝
//                      (文档契约: 调用方预检责任, 库不代查积分阶)
//   4. Kalman 诊断:    statespace v_t → LB/JB (7A residual_diagnostics);
//                      OLS 残差路径同接口打通 (k=0 纯 AR)
//   5. 双差分口径:     ARIMA(0,1,1) 嵌入 vs 先差分 — 两 loglik 数值差异
//                      显式记录 (SS5), 禁混比 (不做相等断言)
//   6. Lyapunov 守卫:  φ=1.01 → stationary_cov throw (非静默) → KF 谱半径
//                      守卫回退 diffuse (message 记录) + 回退与显式 diffuse
//                      递推逐位一致
//
// 数据: 场景 2/5/6 复用 ss_baseline.inc (与 statespace 单测同源夹具;
//   DATA = ARMA(1,1) T=120, DATA300 = ARMA(1,1) φ=0.6/θ=0.4 T=300,
//   gen_phase8_ss_inc.py seed=88); 场景 1/3/4 Philox 确定性生成 (7C
//   gen_noise 同款)
// 断言风格: 通路完整性 + 方向正确性 + 统计合理性 (集成不重复 1e-10 对照,
//   数值精度已由各模块单测锚定; 场景 2 三路径一致性例外 — spec §5 明示
//   "参数一致 1e-6", SLSQP ftol=1e-11 三方同配支撑)
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <vector>

#include "cpphub/core/rng.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/inference/residual_diagnostics.hpp"
#include "cpphub/timeseries/arima/arima_model.hpp"
#include "cpphub/timeseries/ardl/ardl_model.hpp"
#include "cpphub/timeseries/ardl/ardl_select.hpp"
#include "cpphub/timeseries/ardl/uecm_bounds.hpp"
#include "cpphub/timeseries/statespace/arima_statespace.hpp"
#include "cpphub/timeseries/statespace/ss_model.hpp"
#include "cpphub/timeseries/unit_root/adf_test.hpp"
// baseline 夹具 (与 statespace 单测同源)
#include "ss_baseline.inc"

namespace ad = cpphub::v1::timeseries::ardl;
namespace am = cpphub::v1::timeseries::arima;
namespace ss = cpphub::v1::timeseries::statespace;
namespace ur = cpphub::v1::timeseries::unit_root;
namespace em = cpphub::v1::econometrics;
namespace ssb = cpphub::v1::timeseries::ss_baseline;
using cpphub::Real;
using cpphub::Size;

namespace {

std::vector<Real> diff_series(const std::vector<Real>& y) {
    std::vector<Real> d(y.size() - 1);
    for (Size i = 0; i + 1 < y.size(); ++i) d[i] = y[i + 1] - y[i];
    return d;
}

/// N(0,1) 白噪声 (Philox + Box-Muller, 确定性; 7C 同款)
std::vector<Real> gen_noise(Size T, uint64_t seed) {
    cpphub::Philox4x64 rng(seed);
    auto u = [&rng]() {
        return (rng() >> 11) * (1.0 / 9007199254740992.0);
    };
    std::vector<Real> z(T);
    for (Size t = 0; t < T; ++t) {
        const auto [z1, z2] = cpphub::box_muller(u(), u());
        (void)z2;
        z[t] = z1;
    }
    return z;
}

}  // namespace

// ===========================================================================
// 场景 1: ARDL 完整链 select → fit → UECM → bounds F+t
// ===========================================================================
TEST(Phase8Integration, ArdlSelectToBoundsChain) {
    // DGP: x1/x2 独立 RW (I(1)); y = 0.8·x1 + 0.5·x2 + u, u = 0.4·u_{t-1}
    // + e (AR(1) 平稳误差) — 真协整 (β 已知), PSS bounds 标准设定
    const Size T = 250;
    const auto e1 = gen_noise(T, 101);
    const auto e2 = gen_noise(T, 102);
    const auto e3 = gen_noise(T, 103);
    std::vector<Real> x1(T), x2(T), y(T);
    x1[0] = e1[0];
    x2[0] = e2[0];
    Real u_prev = e3[0];
    y[0] = 0.8 * x1[0] + 0.5 * x2[0] + u_prev;
    for (Size t = 1; t < T; ++t) {
        x1[t] = x1[t - 1] + e1[t];
        x2[t] = x2[t - 1] + e2[t];
        const Real u = 0.4 * u_prev + e3[t];
        y[t] = 0.8 * x1[t] + 0.5 * x2[t] + u;
        u_prev = u;
    }

    // (a) IC 网格选阶 (k=2 双 force 变量; A3 同一 hold_back 样本)
    const auto sel = ad::ardl_select_order(y, {x1, x2}, 2, 2, "aic", "c");
    EXPECT_EQ(sel.n_p_grid, Size{3});               // maxlag+1
    ASSERT_EQ(sel.q_grid_dims.size(), Size{2});     // k=2
    EXPECT_EQ(sel.q_grid_dims[0], Size{3});         // maxorder+1 (uniform)
    EXPECT_EQ(sel.aic_grid.size(), Size{27});       // (maxlag+1)·(maxorder+1)^k = 3·3²
    EXPECT_EQ(sel.n_obs_used, T - 2);               // 统一 hb_grid = max(2,2)
    // UECM 前置: max(p, max q_j) ≥ 1 (差分至少损失 1 观测; p=q=0 退化)
    Size maxq = 0;
    for (Size qj : sel.q_opt) maxq = std::max(maxq, qj);
    ASSERT_GE(std::max(sel.p_opt, maxq), Size{1});
    // select → best_model refit 回显一致
    EXPECT_EQ(sel.best_model.spec.p, sel.p_opt);
    ASSERT_EQ(sel.best_model.spec.q.size(), sel.q_opt.size());
    for (Size j = 0; j < sel.q_opt.size(); ++j) {
        EXPECT_EQ(sel.best_model.spec.q[j], sel.q_opt[j]);
    }

    // (b) ARDL → UECM 重参数化 (不重估, 同一估计样本)
    const auto uecm = ad::uecm_from_ardl(sel.best_model);
    EXPECT_EQ(uecm.n_obs_used, sel.best_model.n_obs_used);
    ASSERT_EQ(uecm.pi.size(), Size{3});             // π_y + π_1 + π_2
    // 误差修正方向: 协整 → π_y < 0 (y 对偏离均衡负向调整)
    EXPECT_LT(uecm.pi[0], 0.0);

    // (c) bounds F (Case III, k=2): 真协整 → 拒绝 H0 (无水平关系)
    const auto fb = ad::bounds_f_test(uecm, ad::PssCase::III, 0.05);
    EXPECT_EQ(fb.k, Size{2});
    EXPECT_EQ(fb.n_restrictions, Size{3});  // Case III: π_y + π_1 + π_2
    EXPECT_FALSE(fb.is_t_test);
    EXPECT_TRUE(fb.has_p);
    // 界序结构: I(0) 下界 < I(1) 上界
    EXPECT_LT(fb.crit_vals[0][1], fb.crit_vals[1][1]);
    EXPECT_GT(fb.stat, fb.crit_vals[1][1]);  // F > I(1) 5% 上界
    EXPECT_EQ(fb.conclusion, ad::BoundsConclusion::LevelRelation);
    // p 值对偶面: I(1) null 下 p < 0.05 (响应面近似)
    EXPECT_LT(fb.p_values[1], 0.05);

    // (d) bounds t (Case III): π_y 显著为负 → 同结论 (方向与 F 相反)
    const auto tb = ad::bounds_t_test(uecm, ad::PssCase::III, 0.05);
    EXPECT_TRUE(tb.is_t_test);
    EXPECT_FALSE(tb.has_p);  // t 无响应面 p (A5)
    EXPECT_TRUE(std::isnan(tb.p_values[0]));
    EXPECT_TRUE(std::isnan(tb.p_values[1]));
    // t 界序结构: I(1) 上界更负 (两界均为负)
    EXPECT_GT(tb.crit_vals[0][1], tb.crit_vals[1][1]);
    EXPECT_LT(tb.stat, tb.crit_vals[1][1]);  // t < I(1) 5% 界 (负侧拒绝)
    EXPECT_EQ(tb.conclusion, ad::BoundsConclusion::LevelRelation);
}

// ===========================================================================
// 场景 2: ARIMA 精确似然三路径 (CSS-ML / innovations / statespace)
// ===========================================================================
TEST(Phase8Integration, ArimaThreePathExactLikelihood) {
    // DATA300 = 平稳 ARMA(1,1) (φ=0.6, θ=0.4, gen_phase8_ss_inc.py seed=88)
    // — 三路径共同输入。口径统一 (关键, AR5/AR6):
    //   预去均值 → A CSS-ML 零均值似然 (z 直入);
    //   B Innovations 内部 demean (对已去均值序列 = 数值恒等);
    //   C statespace include_constant=false 零均值 KF
    //   → 三者估计同一零均值 ARMA(1,1) 精确似然, loglik/参数可直接对拍
    const std::vector<Real> raw(ssb::DATA300, ssb::DATA300 + ssb::T300);
    Real mu = 0.0;
    for (Real v : raw) mu += v;
    mu /= static_cast<Real>(raw.size());
    std::vector<Real> z(raw.size());
    for (Size t = 0; t < raw.size(); ++t) z[t] = raw[t] - mu;

    const am::ArimaSpec spec{1, 0, 1};

    // (a) 路径 A: CSS 起始 + innovations 精确似然精化 (CSS-ML)
    const auto ra = am::arima_fit(z, spec, am::ArimaMethod::CSS_ML);
    ASSERT_TRUE(ra.converged);
    EXPECT_EQ(ra.n_obs_used, ssb::T300);  // d=0 全样本

    // (b) 路径 B: innovations 算法直接精确 MLE
    const auto rb = am::arima_fit(z, spec, am::ArimaMethod::Innovations);
    ASSERT_TRUE(rb.converged);

    // (c) 路径 C: statespace KF 精确似然 MLE (stationary init, burn=0)
    const auto rc = ss::arima_statespace_mle(z, spec);
    ASSERT_TRUE(rc.converged);
    EXPECT_EQ(rc.loglikelihood_burn, Size{0});  // d=0 且 AR 平稳

    // (d) 三路径落点一致: SLSQP ftol=1e-11 三方同配 (A/B/C 各自优化器
    //     收敛至同一最优), spec §5 "参数一致 1e-6"
    EXPECT_NEAR(ra.loglik, rb.loglik, 1e-6);
    EXPECT_NEAR(ra.loglik, rc.loglik, 1e-6);
    EXPECT_NEAR(rb.loglik, rc.loglik, 1e-6);
    EXPECT_NEAR(ra.params.phi[0], rb.params.phi[0], 1e-6);
    EXPECT_NEAR(ra.params.phi[0], rc.phi[0], 1e-6);
    EXPECT_NEAR(ra.params.theta[0], rb.params.theta[0], 1e-6);
    EXPECT_NEAR(ra.params.theta[0], rc.theta[0], 1e-6);
    // σ²/AIC 同口径 (k = p+q+1 含 σ²: AR4 与 K4 一致)
    EXPECT_NEAR(ra.params.sigma2, rb.params.sigma2, 1e-6);
    EXPECT_NEAR(ra.params.sigma2, rc.sigma2, 1e-6);
    EXPECT_NEAR(ra.aic, rc.aic, 1e-5);

    // (e) K1 常驻互检: C 内部 KF vs innovations gap (告警线 1e-5)
    EXPECT_LT(rc.innovations_gap, 1e-5);

    // (f) 统计合理性: DGP (φ, θ) = (0.6, 0.4) 邻域。φ 方向良态识别;
    //     θ 方向弱识别 — statsmodels 同 DGP 样本 (seed=88 族) T=120 的
    //     SARIMAX MLE θ̂ = 0.6199 (ss_baseline MLE_ARMA11_PARAMS[1]),
    //     本库 T=300 落 0.567 (向真值收敛) — 偏离系样本性质非实现偏差
    //     (实现正确性由 (d) 三路径 1e-6 对拍 + (e) K1 互检锚定)
    EXPECT_NEAR(rc.phi[0], 0.6, 0.15);
    EXPECT_NEAR(rc.theta[0], 0.4, 0.3);
}

// ===========================================================================
// 场景 3: 单位根预检 → 边限工作流 (混合 I(0)/I(1); I(2) 前置拒绝)
// ===========================================================================
TEST(Phase8Integration, UnitRootPrecheckBoundsWorkflow) {
    // DGP: x1 = AR(1) φ=0.5 (I(0)); x2 = RW (I(1)); x3 = cumsum(RW) (I(2));
    //      y = 0.6·x1 + 0.7·x2 + u (AR(1) 平稳误差) — 混合 I(0)/I(1) force
    //      变量 = PSS bounds 适用域 (边限检验对 force 变量积分阶混合稳健,
    //      但任何变量 I(2) 即越界 — 前置预检是调用方责任, 文档契约)
    const Size T = 250;
    const auto e1 = gen_noise(T, 201);
    const auto e2 = gen_noise(T, 202);
    const auto e3 = gen_noise(T, 203);
    const auto e4 = gen_noise(T, 204);
    std::vector<Real> x1(T), x2(T), x3(T), y(T);
    x1[0] = e1[0];
    x2[0] = e2[0];
    Real w = e3[0];  // x3 的增量载体 (自身 I(1))
    x3[0] = w;
    Real u_prev = e4[0];
    y[0] = 0.6 * x1[0] + 0.7 * x2[0] + u_prev;
    for (Size t = 1; t < T; ++t) {
        x1[t] = 0.5 * x1[t - 1] + e1[t];  // I(0)
        x2[t] = x2[t - 1] + e2[t];        // I(1)
        w += e3[t];                       // w: I(1)
        x3[t] = x3[t - 1] + w;            // x3 = cumsum(w): I(2)
        const Real u = 0.4 * u_prev + e4[t];
        y[t] = 0.6 * x1[t] + 0.7 * x2[t] + u;
        u_prev = u;
    }

    // (a) ADF 预检链 (7B 复用, auto lag + AIC): 逐变量积分阶分类
    const auto adf_x1 = ur::adf_test(x1, "c", 0, true);
    EXPECT_TRUE(adf_x1.reject_null);  // I(0)
    const auto adf_x2 = ur::adf_test(x2, "c", 0, true);
    EXPECT_FALSE(adf_x2.reject_null);  // 非 I(0) → 查差分
    const auto adf_dx2 = ur::adf_test(diff_series(x2), "c", 0, true);
    EXPECT_TRUE(adf_dx2.reject_null);  // Δx2 平稳 → x2 恰 I(1) ✓

    // (b) I(2) 检出 → 前置拒绝 (守门在调用方; 库无积分阶校验 — 契约演示)
    const auto adf_x3 = ur::adf_test(x3, "c", 0, true);
    ASSERT_FALSE(adf_x3.reject_null);  // 水平非平稳
    const auto d_x3 = diff_series(x3);
    const auto adf_dx3 = ur::adf_test(d_x3, "c", 0, true);
    ASSERT_FALSE(adf_dx3.reject_null);  // Δx3 仍单位根 → 至少 I(2)
    const auto adf_ddx3 = ur::adf_test(diff_series(d_x3), "c", 0, true);
    EXPECT_TRUE(adf_ddx3.reject_null);  // Δ²x3 平稳 → 恰 I(2) → 拒绝入列

    // (c) 有效混合集合 {x1, x2} → ARDL(1, 1, 1) → UECM → bounds (Case III)
    ad::ArdlSpec sp;
    sp.p = 1;
    sp.q = {1, 1};
    sp.trend = "c";
    const auto fit = ad::ardl_fit(y, {x1, x2}, sp);
    const auto uecm = ad::uecm_from_ardl(fit);
    EXPECT_LT(uecm.pi[0], 0.0);  // 误差修正方向
    const auto fb = ad::bounds_f_test(uecm, ad::PssCase::III, 0.05);
    EXPECT_EQ(fb.k, Size{2});
    EXPECT_GT(fb.stat, fb.crit_vals[1][1]);
    EXPECT_EQ(fb.conclusion, ad::BoundsConclusion::LevelRelation);
}

// ===========================================================================
// 场景 4: Kalman 诊断 — statespace v_t → LB/JB; OLS 残差路径打通
// ===========================================================================
TEST(Phase8Integration, KalmanInnovationsDiagnostics) {
    // DGP: 零均值 ARMA(1,1) φ=0.5, θ=0.3 (Philox 确定性)
    const Size T = 300;
    const auto e = gen_noise(T, 301);
    std::vector<Real> z(T);
    Real prev_e = 0.0;
    for (Size t = 0; t < T; ++t) {
        z[t] = (t > 0 ? 0.5 * z[t - 1] : 0.0) + e[t] + 0.3 * prev_e;
        prev_e = e[t];
    }

    // (a) statespace MLE → KF innovations v_t (全长)
    const am::ArimaSpec spec{1, 0, 1};
    const auto r = ss::arima_statespace_mle(z, spec);
    ASSERT_TRUE(r.converged);
    ASSERT_EQ(r.innovations.size(), T);

    // (b) v_t → 7A LB 检验: 正确设定 → 残差白噪声, 不拒绝
    const auto lb = em::ljung_box_test(r.innovations);
    EXPECT_EQ(lb.base.method_name, "Ljung-Box");
    EXPECT_GT(lb.base.p_value, 0.05);
    EXPECT_FALSE(lb.base.reject_null);

    // (c) v_t → 7A JB 正态性: 通路完整 + 字段有效 (高斯 DGP → 不强拒绝)
    const auto jb = em::jarque_bera_test(r.innovations);
    EXPECT_EQ(jb.base.method_name, "Jarque-Bera");
    EXPECT_TRUE(std::isfinite(jb.base.statistic));
    EXPECT_GE(jb.base.p_value, 0.0);
    EXPECT_LE(jb.base.p_value, 1.0);
    EXPECT_TRUE(std::isfinite(jb.skewness));
    EXPECT_TRUE(std::isfinite(jb.kurtosis));

    // (d) OLS 残差路径打通 (同一诊断接口): 纯 AR(3) OLS (ardl_fit k=0
    //     退化 = AR(p); ARMA(1,1) 的 AR(∞) 截断近似) 残差 → LB
    ad::ArdlSpec as;
    as.p = 3;
    as.trend = "c";
    const auto ar3 = ad::ardl_fit(z, {}, as);
    ASSERT_EQ(ar3.residuals.size(), T - 3);  // 估计样本 = T − p
    const auto lb_ols = em::ljung_box_test(ar3.residuals);
    EXPECT_GT(lb_ols.base.p_value, 0.05);
}

// ===========================================================================
// 场景 5: 双差分口径 — ARIMA(0,1,1) 嵌入 vs 先差分 (SS5 禁混比)
// ===========================================================================
TEST(Phase8Integration, DualDifferencingCaliber) {
    // data120 (与 statespace 单测同源); 嵌入 = levels 直入 (d=1 差分入
    // 状态, 混合初始化, burn=d=1); simple = 先差分样本 (T−1) + ARMA(0,1)
    // 形 (stationary init, burn=0)
    const std::vector<Real> y(ssb::DATA, ssb::DATA + ssb::T);
    const std::vector<Real> dy = diff_series(y);

    const auto me =
        ss::make_arima_statespace({}, {0.35}, 1, 0.9, false, false);
    const auto re = ss::kalman_filter(y, me, 1);
    ASSERT_TRUE(re.ok);

    const auto ms =
        ss::make_arima_statespace({}, {0.35}, 1, 0.9, false, true);
    const auto rs = ss::kalman_filter(dy, ms, 0);
    ASSERT_TRUE(rs.ok);

    // 样本口径回显: embedded 全长 T (burn 1 期不计入 llf), simple T−1
    EXPECT_EQ(re.n_obs, ssb::T);
    EXPECT_EQ(rs.n_obs, ssb::T - 1);

    // SS5: 两口径似然不同 (非恒等) 但时间平移等价 — gap = O(1/κ)
    // (κ=1e6 有限 diffuse 尾差, T=120 实测 ~1e-6 量级); 数值差异显式
    // 记录, 不做相等断言 (禁混比)
    const Real gap = std::fabs(re.loglik - rs.loglik);
    EXPECT_GT(gap, 0.0);
    EXPECT_LT(gap, 1e-5);
}

// ===========================================================================
// 场景 6: Lyapunov 守卫端到端 — φ=1.01 throw + KF diffuse 回退正确性
// ===========================================================================
TEST(Phase8Integration, LyapunovGuardEndToEnd) {
    // AR(1) φ=1.01 (单位根邻域): Lyapunov 解 P = σ²/(1−φ²) < 0 —
    // 静默负方差 = 幻觉温床, SS4 守卫必须非静默
    const std::vector<Real> y(ssb::DATA, ssb::DATA + ssb::T);
    const auto m = ss::make_arma_statespace({1.01}, {}, 1.0);

    // (a) 守卫非静默 1: stationary_cov 直接求解 → throw (非正定)
    EXPECT_THROW(ss::stationary_cov(m), std::runtime_error);

    // (b) 守卫非静默 2: KF 层 Stationary init → 谱半径守卫 (|λ| ≥ 1)
    //     自动回退 ApproximateDiffuse, message 记录 (非静默回退)
    const auto r = ss::kalman_filter(y, m, 0);
    EXPECT_TRUE(r.ok);
    EXPECT_NE(r.message.find("fell back"), std::string::npos);
    EXPECT_TRUE(std::isfinite(r.loglik));

    // (c) 回退正确: 显式 ApproximateDiffuse (同 κ=1e6) 与自动回退的
    //     递推初始化逐位相同 → loglik 恒等; 仅 message 差异 (显式路径
    //     无告警)
    auto m2 = m;
    m2.init.method = ss::SSInitMethod::ApproximateDiffuse;
    const auto r2 = ss::kalman_filter(y, m2, 0);
    EXPECT_TRUE(r2.ok);
    EXPECT_TRUE(r2.message.empty());
    EXPECT_NEAR(r.loglik, r2.loglik, 1e-12);
}
