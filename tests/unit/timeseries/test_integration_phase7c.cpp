// =============================================================================
// test_integration_phase7c.cpp - Phase 7C 端到端集成测试 (spec §7, 6 场景)
//
// 场景矩阵 (spec §7):
//   1. 单位根诊断全链: ADF/DF-GLS/NP 组合 → 差分决策 → ARIMA → 残差 LB;
//      断点→ZA→分段平稳 (M0+M1+7B 复用)
//   2. Granger 因果链:  I(1) 双序列 — 差分标准 F vs 水平 TY 增广 Wald 对比
//      (GR7: 共同 RW 因子 → 水平标准 F 虚假显著, 两条合法路径修正)
//   3. VAR→DY 溢出:    多资产 → IC 选阶 → 稳定性 → IRF/FEVD 双轨 → DY + 滚动
//      (M2 全链; V12 不稳定拦截)
//   4. 协整→VECM:      EG/Johansen rank → VECM → ECT 显著性 → β 投影空间
//      (M3 全链, CI8/CI9)
//   5. MIDAS 混频预测: 月度 y + 日度 x → MIDAS-DL vs U-MIDAS → MZ/DM (7B 复用);
//      MD3 期初起窗
//   6. GARCH-M 风险溢价: 收益率 → 三变体 → λ sandwich 显著性 (GM4) →
//      vs 无 M 模型 (AIC 对照)
//
// 数据: 全部复用各模块已锚定 baseline 夹具 (与单测同源, 集成只换链路视角);
//   断点/爆炸序列用 Philox 确定性生成
// 断言风格: 通路完整性 + 方向正确性 + 统计合理性 (集成不重复 1e-10 对照,
//   数值精度已由各模块单测锚定)
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "cpphub/core/rng.hpp"
#include "cpphub/core/types.hpp"
// 场景 1: 单位根 + ARIMA + LB
#include "cpphub/timeseries/arima/arima_model.hpp"
#include "cpphub/timeseries/unit_root/adf_test.hpp"
#include "cpphub/timeseries/unit_root/df_gls_test.hpp"
#include "cpphub/timeseries/unit_root/ng_perron_test.hpp"
#include "cpphub/timeseries/unit_root/zivot_andrews_test.hpp"
#include "cpphub/econometrics/inference/residual_diagnostics.hpp"
// 场景 2: Granger
#include "cpphub/timeseries/arima/granger_test.hpp"
// 场景 3: VAR/DY
#include "cpphub/timeseries/var/dy_spillover.hpp"
#include "cpphub/timeseries/var/fevd.hpp"
#include "cpphub/timeseries/var/irf.hpp"
#include "cpphub/timeseries/var/var_model.hpp"
#include "cpphub/timeseries/var/var_select.hpp"
// 场景 4: 协整/VECM
#include "cpphub/timeseries/cointegration/engle_granger.hpp"
#include "cpphub/timeseries/cointegration/johansen_test.hpp"
#include "cpphub/timeseries/cointegration/vecm_model.hpp"
// 场景 5: MIDAS + MZ/DM
#include "cpphub/timeseries/midas/midas_model.hpp"
#include "cpphub/timeseries/midas/mixed_freq_data.hpp"
#include "cpphub/econometrics/inference/specification_tests.hpp"
// 场景 6: GARCH-M vs GARCH
#include "cpphub/timeseries/garch/garch_m_model.hpp"
#include "cpphub/timeseries/garch/garch_model.hpp"
// baseline 夹具 (与各模块单测同源)
#include "coint_baseline.inc"
#include "gm_baseline.inc"
#include "granger_baseline.inc"
#include "midas_baseline.inc"
#include "unit_root_baseline.inc"
#include "var_baseline.inc"

namespace ur = cpphub::v1::timeseries::unit_root;
namespace am = cpphub::v1::timeseries::arima;
namespace gr = cpphub::v1::timeseries::granger;
namespace vt = cpphub::v1::timeseries::var;
namespace co = cpphub::v1::timeseries::cointegration;
namespace mm = cpphub::v1::timeseries::midas;
namespace ts = cpphub::v1::timeseries::garch;
namespace em = cpphub::v1::econometrics;
namespace urb = cpphub::v1::timeseries::unit_root::baseline;
namespace vb = cpphub::v1::timeseries::var_baseline;
namespace cb = cpphub::v1::timeseries::coint_baseline;
namespace gb = cpphub::v1::timeseries::granger_baseline;
namespace mb = cpphub::v1::timeseries::midas_baseline::v1;
namespace gmb = cpphub::v1::timeseries::garch::gm_baseline;
using cpphub::Real;
using cpphub::Size;

namespace {

std::vector<Real> diff_series(const std::vector<Real>& y) {
    std::vector<Real> d(y.size() - 1);
    for (Size i = 0; i + 1 < y.size(); ++i) d[i] = y[i + 1] - y[i];
    return d;
}

/// N(0,1) 白噪声 (Philox + Box-Muller, 确定性)
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
// 场景 1: 单位根诊断全链 (组合检验 → 差分决策 → ARIMA → 残差 LB; 断点→分段)
// ===========================================================================
TEST(Phase7CIntegration, UnitRootDiagnosticChain) {
    const std::vector<Real> rw(urb::Y_RW, urb::Y_RW + urb::T);

    // (a) M0 组合诊断: ADF + DF-GLS + Ng-Perron 三检验一致不拒绝 → I(1)
    //     (7B 场景 4 已加 KPSS+BH; 此处引入 M0 新回填的 NP 与 ZA 全链)
    const auto adf = ur::adf_test(rw, "c", 0, true);
    const auto dfgls = ur::df_gls_test(rw, "ct");
    const auto np = ur::ng_perron_test(rw, "ct");
    EXPECT_FALSE(adf.reject_null);         // ADF: 不拒绝单位根
    EXPECT_FALSE(dfgls.reject_null);       // DF-GLS: 不拒绝
    EXPECT_FALSE(np.reject_5pct[0]);       // NP MZα: 不拒绝
    EXPECT_FALSE(np.reject_5pct[1]);       // NP MZt: 不拒绝
    // NP 内部恒等式 (通路自检): MZt = MZα × MSB
    EXPECT_NEAR(np.mz_t, np.mz_alpha * np.msb, 1e-12);

    // (b) 差分决策 → 复检: Δy 平稳 → 判定 I(1), 建模转入差分域
    const auto dy = diff_series(rw);
    const auto adf_dy = ur::adf_test(dy, "c", 0, true);
    EXPECT_TRUE(adf_dy.reject_null);

    // (c) ARIMA(1,1,0) CSS-ML: 水平值直接拟合 (d=1 内部差分)
    const auto fit = am::arima_fit(rw, am::ArimaSpec{1, 1, 0},
                                   am::ArimaMethod::CSS_ML);
    ASSERT_TRUE(fit.converged);
    EXPECT_EQ(fit.n_obs_used, urb::T - 1);   // T − d
    // RW 差分后白噪声 → AR(1) 系数应近 0 (|φ| < 0.2)
    EXPECT_LT(std::fabs(fit.params.phi[0]), 0.2);

    // (d) 残差诊断收口: LB 检验 p > 0.05 → 残差白噪声, 模型充分
    const auto lb = em::ljung_box_test(fit.residuals);
    EXPECT_EQ(lb.base.method_name, "Ljung-Box");
    EXPECT_GT(lb.base.p_value, 0.05);
    EXPECT_FALSE(lb.base.reject_null);

    // (e) 断点→分段平稳 (Perron 1989 动机): 前 120 期 N(0,1), 后 120 期
    //     N(2,1) (均值断点 2σ, 序列本身分段平稳)。标准 ADF 检验力被断点
    //     吸收 → 不拒绝; ZA Model A (崩溃均值) 识别断点后拒绝单位根;
    //     断点前后分段各自 ADF 均拒绝。
    const Size Tb = 120, Ts = 240;
    std::vector<Real> broken(Ts);
    const auto noise = gen_noise(Ts, 2024u);
    for (Size t = 0; t < Ts; ++t) {
        broken[t] = noise[t] + (t >= Tb ? 2.0 : 0.0);
    }
    const auto adf_b = ur::adf_test(broken, "c", 0, true);
    EXPECT_FALSE(adf_b.reject_null);       // 均值断点扭曲标准 ADF (经典结果)
    const auto za = ur::zivot_andrews_test(broken, ur::ZAModel::A, 0, true);
    EXPECT_TRUE(za.reject_null);           // ZA 吸收断点后拒绝单位根
    // 断点定位: 0-based break_index 应落在真实断点附近 (±15 期)
    EXPECT_GE(za.break_index, Tb - 15);
    EXPECT_LE(za.break_index, Tb + 15);
    // ZA 统计量比标准 ADF 更负 (断点吸收 → 检验力恢复)
    EXPECT_LT(za.statistic, adf_b.statistic);
    // 分段平稳确认: 断点前后各自 ADF 拒绝单位根
    const std::vector<Real> seg1(broken.begin(), broken.begin() + Tb);
    const std::vector<Real> seg2(broken.begin() + Tb, broken.end());
    EXPECT_TRUE(ur::adf_test(seg1, "c", 0, true).reject_null);
    EXPECT_TRUE(ur::adf_test(seg2, "c", 0, true).reject_null);
}

// ===========================================================================
// 场景 2: Granger 因果链 — I(1) 双序列三条推断路径对比 (GR7)
//   DGP: y1 = common + 0.3e1, y2 = common + 0.3e2 (common = RW)
//   → 水平标准 F: 名义显著但参照分布失效 (I(1) 伪回归, GR7)
//   → 差分标准 F: 合法但检验的是差分域短期因果 — 差分消去共同趋势 → 无因果
//   → 水平 TY 增广 Wald: 合法且保留水平信息 — 协整系统 (场景 4 rank=1)
//     的水平因果真实存在 (共同 RW 历史经 y2 滞后传递) → 显著
//   两合法路径给出不同答案: 检验假设不同 (短期 vs 水平长短期), 均正确
// ===========================================================================
TEST(Phase7CIntegration, GrangerCausalityChainI1) {
    const std::vector<Real> y1(std::begin(gb::Y1), std::end(gb::Y1));
    const std::vector<Real> y2(std::begin(gb::Y2), std::end(gb::Y2));
    const Size T_i1 = y1.size();

    // (a) 诊断前提: 双序列均 I(1) (差分后平稳)
    EXPECT_FALSE(ur::adf_test(y1, "c", 0, true).reject_null);
    EXPECT_FALSE(ur::adf_test(y2, "c", 0, true).reject_null);
    EXPECT_TRUE(ur::adf_test(diff_series(y1), "c", 0, true).reject_null);
    EXPECT_TRUE(ur::adf_test(diff_series(y2), "c", 0, true).reject_null);

    // (b) 水平值标准 F (GR7 无效路径, 作对照): 名义显著
    const auto lvl = gr::granger_test(y1, y2, 1, 0, false);
    EXPECT_FALSE(lvl.has_ty);
    EXPECT_TRUE(std::isfinite(lvl.f_stat));
    // 共同 RW 因子 → y2 滞后对 y1 有"预测力" → 名义 p 显著;
    // 但 I(1) 水平值下 F/χ² 参照分布失效, 该显著性不可作推断
    EXPECT_LT(lvl.f_p, 0.05);

    // (c) 路径 1 — 差分标准 F (合法): 差分消去共同随机游走 → 短期无因果
    const auto dif = gr::granger_test(diff_series(y1), diff_series(y2),
                                      1, 0, false);
    EXPECT_TRUE(std::isfinite(dif.f_stat));
    EXPECT_GT(dif.f_p, 0.05);              // 差分域无 Granger 因果
    // df 语义: T−1 个差分观测, df2 = (T−1) − 1 − (2·1+1)
    EXPECT_EQ(dif.df2, T_i1 - 1 - 1 - 3);

    // (d) 路径 2 — 水平 TY 增广 Wald (合法, d_max=1): 水平因果显著 —
    //     与 Johansen rank=1 (场景 4) 自洽: 协整 ⇒ 至少单向 Granger 因果
    const auto ty = gr::granger_test(y1, y2, 1, 1, false);
    EXPECT_TRUE(ty.has_ty);
    EXPECT_TRUE(std::isfinite(ty.ty_wald_stat));
    EXPECT_LT(ty.ty_wald_p, 0.05);

    // (e) 三路径统计量互异 (推断路径分离的核心证据)
    EXPECT_NE(lvl.f_stat, dif.f_stat);
    EXPECT_NE(lvl.f_stat, ty.ty_wald_stat);
    EXPECT_NE(dif.f_stat, ty.ty_wald_stat);
    // 结论: 两合法路径按各自假设均正确 — 差分 F "短期无因果" (共同趋势
    // 已消去), TY "水平有因果" (协整信息保留); 水平 F 的显著结论数值上
    // 与 TY 同向, 但其参照分布无效 (GR7) — 合法推断必须走差分 F 或 TY
}

// ===========================================================================
// 场景 3: VAR→DY 溢出全链 (IC 选阶 → 稳定性 → IRF/FEVD 双轨 → DY + 滚动)
// ===========================================================================
TEST(Phase7CIntegration, VarToDySpilloverChain) {
    vt::MultivariateTSData d;
    d.columns = {std::vector<Real>(std::begin(vb::Y1), std::end(vb::Y1)),
                 std::vector<Real>(std::begin(vb::Y2), std::end(vb::Y2)),
                 std::vector<Real>(std::begin(vb::Y3), std::end(vb::Y3))};
    d.names = {"y1", "y2", "y3"};
    const Size T_var = d.T();

    // (a) IC 选阶: 四准则一致选 p=2 (DGP 为 VAR(2))
    for (const char* ic : {"aic", "bic", "hqic", "fpe"}) {
        const auto sel = vt::var_select_order(d, "c", 4, ic);
        EXPECT_EQ(sel.selected_lag, 2u) << ic;
    }

    // (b) VAR(2) 估计 + 稳定性 (V9): 平稳才可做 IRF/FEVD
    vt::VARSpec spec;
    spec.lag = 2;
    const auto fit = vt::var_fit(d, spec);
    EXPECT_EQ(fit.n_obs_used, T_var - 2);
    EXPECT_LT(fit.max_abs_eigenvalue, 1.0);
    EXPECT_TRUE(fit.is_strictly_stationary);

    // (c) IRF (V2/V3): Φ_0 = I; Ψ_0 = Cholesky P 下三角; 正交化方向非对称
    const auto irf = vt::var_irf(fit, 11);
    ASSERT_EQ(irf.theta.size(), 11u);
    EXPECT_DOUBLE_EQ(irf.phi[0](0, 0), 1.0);
    EXPECT_DOUBLE_EQ(irf.phi[0](1, 0), 0.0);
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            EXPECT_DOUBLE_EQ(irf.theta[0](i, j), 0.0);   // 上三角严格 0
    EXPECT_NE(irf.theta[1](0, 2), irf.theta[1](2, 0));   // 行=响应, 列=冲击

    // (d) FEVD 双轨 (V7/V8): Cholesky 与 DY 行和=1; PS 未归一
    for (auto fw : {vt::FevdFramework::Cholesky,
                    vt::FevdFramework::GeneralizedDY}) {
        const auto f = vt::var_fevd(fit, 10, fw);
        for (int i = 0; i < 3; ++i)
            EXPECT_NEAR(f.fevd.row(i).sum(), 1.0, 1e-12);
    }
    const auto f_ps = vt::var_fevd(fit, 10, vt::FevdFramework::GeneralizedPS);
    for (int i = 0; i < 3; ++i)
        EXPECT_NE(f_ps.fevd.row(i).sum(), 1.0);

    // (e) DY 静态溢出: TCI ∈ (0,100); 两个 DY2012 表恒等式:
    //     NET_j = TO_j − FROM_j; Σ_j TO_j = Σ_j FROM_j = TCI
    const auto sp = vt::dy_spillover_static(d, 10, "c", 2);
    EXPECT_GT(sp.tci, 0.0);
    EXPECT_LT(sp.tci, 100.0);
    Real sum_to = 0.0, sum_from = 0.0;
    for (int j = 0; j < 3; ++j) {
        EXPECT_NEAR(sp.net_spillover[j],
                    sp.to_spillover[j] - sp.from_spillover[j], 1e-12);
        sum_to += sp.to_spillover[j];
        sum_from += sp.from_spillover[j];
    }
    EXPECT_NEAR(sum_to, sp.tci, 1e-8);
    EXPECT_NEAR(sum_from, sp.tci, 1e-8);

    // (f) 滚动溢出: window=150 → 101 窗口; 末窗口 == 手动子样本静态
    const auto roll = vt::dy_spillover(d, 150, 10, "c", 2);
    ASSERT_EQ(roll.tci_path.size(), T_var - 150 + 1);
    EXPECT_NEAR(roll.tci, roll.tci_path.back(), 0.0);
    vt::MultivariateTSData sub;
    sub.columns.resize(3);
    for (int j = 0; j < 3; ++j) {
        sub.columns[j].assign(d.columns[j].begin() + (T_var - 150),
                              d.columns[j].end());
    }
    const auto tail = vt::dy_spillover_static(sub, 10, "c", 2);
    EXPECT_NEAR(roll.tci, tail.tci, 1e-10);
    // 滚动路径全有限 (波动率状态演化)
    for (Real v : roll.tci_path) EXPECT_TRUE(std::isfinite(v));

    // (g) V12 拦截: 爆炸 VAR (ρ=1.05, 双列独立噪声避免共线奇异)
    //     → max|eig|>1 → FEVD 拒绝 (溢出指数对不稳定系统无意义)
    const auto boom_n1 = gen_noise(100, 77u);
    const auto boom_n2 = gen_noise(100, 78u);
    vt::MultivariateTSData bad;
    bad.columns.resize(2);
    bad.columns[0].resize(100);
    bad.columns[1].resize(100);
    Real x = 0.0, yv = 0.0;
    for (Size t = 0; t < 100; ++t) {
        x = 1.05 * x + boom_n1[t];
        yv = 1.05 * yv + 0.4 * boom_n1[t] + boom_n2[t];
        bad.columns[0][t] = x;
        bad.columns[1][t] = yv;
    }
    vt::VARSpec spec_bad;
    spec_bad.lag = 1;
    const auto fit_bad = vt::var_fit(bad, spec_bad);
    EXPECT_FALSE(fit_bad.is_strictly_stationary);
    EXPECT_THROW(vt::var_fevd(fit_bad, 10, vt::FevdFramework::GeneralizedDY),
                 std::invalid_argument);
}

// ===========================================================================
// 场景 4: 协整→VECM 全链 (EG 筛选 → Johansen rank → VECM → ECT → β 空间)
// ===========================================================================
TEST(Phase7CIntegration, CointegrationToVecmChain) {
    vt::MultivariateTSData d;
    d.columns = {std::vector<Real>(std::begin(cb::Y1), std::end(cb::Y1)),
                 std::vector<Real>(std::begin(cb::Y2), std::end(cb::Y2)),
                 std::vector<Real>(std::begin(cb::Y3), std::end(cb::Y3))};
    d.names = {"y1", "y2", "y3"};
    const Size T_co = d.T();
    const std::vector<Real> y1(std::begin(cb::Y1), std::end(cb::Y1));
    const std::vector<Real> y2(std::begin(cb::Y2), std::end(cb::Y2));
    const std::vector<Real> y3(std::begin(cb::Y3), std::end(cb::Y3));

    // (a) EG 双变量筛选: (y1,y2) 拒绝无协整; (y1,y3) 对照不拒绝
    const auto eg12 = co::engle_granger(y1, y2, "c");
    const auto eg13 = co::engle_granger(y1, y3, "c");
    EXPECT_TRUE(eg12.reject_null);
    EXPECT_FALSE(eg13.reject_null);

    // (b) Johansen 多变量定秩: trace 5% → rank = 1 (与 EG 一致)
    const auto jo = co::coint_johansen(d, 0, 1);
    EXPECT_GT(jo.lr1(0), jo.cvt(0, 1));    // 拒绝 r=0
    EXPECT_LT(jo.lr1(1), jo.cvt(1, 1));    // 接受 r=1
    const Size r = co::select_coint_rank(d, 0, 1, "trace", 0.05);
    EXPECT_EQ(r, 1u);

    // (c) VECM(r=1, k=1, co): 夹具 DGP = 共同 RW + 独立噪声, 无拉回机制
    //     → α 载入小 (弱外生): 协整存在 (谱相关) 但不要求 ECT 显著;
    //     断言检验通路完整: 统计量与 EM2002 左尾临界值均有限
    const auto vecm = co::vecm_fit(d, r, 1, "co");
    EXPECT_EQ(vecm.rank, 1u);
    EXPECT_EQ(vecm.n_obs, T_co - 1 - 1);
    ASSERT_TRUE(vecm.has_ect_t);
    ASSERT_EQ(vecm.ect_t_stat.size(), 3u);
    for (int i = 0; i < 3; ++i) {
        EXPECT_TRUE(std::isfinite(vecm.ect_t_stat[i]));
        EXPECT_LT(vecm.ect_cv_5pct[i], 0.0);   // EM2002 左尾临界值
    }

    // (c') 对照: 真实拉回 DGP (VECM 生成, rank=1) → ECT 显著负
    //     Δy1 = −0.30·ect + e1;  Δy2 = −0.20·ect + e2  (ect = y1−y2)
    const auto ec_n1 = gen_noise(250, 5u);
    const auto ec_n2 = gen_noise(250, 6u);
    vt::MultivariateTSData ec;
    ec.columns.resize(2);
    ec.columns[0].resize(250);
    ec.columns[1].resize(250);
    Real a = 0.0, b = 0.0;
    for (Size t = 0; t < 250; ++t) {
        const Real ect = a - b;
        a += -0.30 * ect + ec_n1[t];
        b += -0.20 * ect + ec_n2[t];
        ec.columns[0][t] = a;
        ec.columns[1][t] = b;
    }
    const auto vecm_ec = co::vecm_fit(ec, 1, 0, "n");
    EXPECT_EQ(vecm_ec.n_obs, 250u - 0 - 1);
    ASSERT_TRUE(vecm_ec.has_ect_t);
    // 真实拉回被检出: 两方程 ECT t 均显著负 (t < EM2002 5% 左尾临界值)
    EXPECT_LT(vecm_ec.ect_t_stat[0], vecm_ec.ect_cv_5pct[0]);
    EXPECT_LT(vecm_ec.ect_t_stat[1], vecm_ec.ect_cv_5pct[1]);
    // α 符号与 DGP 一致 (负拉回)
    EXPECT_LT(vecm_ec.alpha(0, 0), 0.0);
    EXPECT_LT(vecm_ec.alpha(1, 0), 0.0);

    // (d) β 投影空间 (CI8): 默认归一 vs urca 归一 → P = β(β'β)⁻¹β' 一致
    const auto vecm_ur = co::vecm_fit(d, r, 1, "co", true);
    EXPECT_TRUE(vecm_ur.urca_normalization);
    EXPECT_DOUBLE_EQ(vecm_ur.beta(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(vecm.beta(0, 0), 1.0);    // r=1 时两归一数值相同
    auto projector = [](const Eigen::MatrixXd& b) {
        return b * (b.transpose() * b).inverse() * b.transpose();
    };
    const Eigen::MatrixXd p1 = projector(vecm.beta);
    const Eigen::MatrixXd p2 = projector(vecm_ur.beta);
    EXPECT_TRUE(p1.isApprox(p2, 1e-10));

    // (e) Π = αβ' 秩不变性 (CI9): SVD 奇异值恰 1 个非零 (rank=1)
    const Eigen::MatrixXd pi = vecm.alpha * vecm.beta.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(pi);
    EXPECT_GT(svd.singularValues()(0), 1e-8);
    EXPECT_LT(svd.singularValues()(1), 1e-8);
    EXPECT_LT(svd.singularValues()(2), 1e-8);
    // 归一不变: Π(默认) ≈ Π(urca)
    const Eigen::MatrixXd pi_ur = vecm_ur.alpha * vecm_ur.beta.transpose();
    EXPECT_TRUE(pi.isApprox(pi_ur, 1e-10));
}

// ===========================================================================
// 场景 5: MIDAS 混频预测链 (数据对齐 → DL vs U-MIDAS → MZ/DM 精度对比)
// ===========================================================================
TEST(Phase7CIntegration, MidasMixedFreqForecastChain) {
    mm::MixedFreqData data;
    data.y.assign(mb::Y, mb::Y + mb::N_LF);
    data.x.assign(mb::X, mb::X + mb::N_HF);
    data.m = mb::M;
    ASSERT_NO_THROW(data.validate());

    // (a) MD3 期初起窗: h_start=1 → 第一行 lag0 = x[m·t−1] (期初观测)
    const auto X_begin = mm::design_matrix(data, 2, 1);
    ASSERT_EQ(X_begin.size(), mm::aligned_y(data, 2, 1).size());
    // 期末对齐 (h_start=0): lag0 = x[m·t] (期末最新)
    const auto X_end = mm::design_matrix(data, 2, 0);
    EXPECT_NE(X_begin[0][0], X_end[0][0]);     // 两起窗取到不同高频观测

    // (b) 双模型拟合: MIDAS-DL (Nealmon 约束) vs U-MIDAS (无约束 OLS)
    const auto dl = mm::midas_fit(data, mm::MidasWeight::Nealmon,
                                  mm::MidasType::DL, 4);
    const auto um = mm::midas_fit(data, mm::MidasWeight::Nealmon,
                                  mm::MidasType::UMidas, 4);
    ASSERT_TRUE(dl.converged);
    ASSERT_TRUE(um.converged);
    EXPECT_EQ(um.lambda.empty(), true);        // U-MIDAS 无外层超参
    // 嵌套关系: 无约束 SSR ≤ 约束 NLS SSR (同一设计空间)
    EXPECT_LE(um.ssr, dl.ssr + 1e-9);

    // (c) MZ 精度回归 (7B 复用): actual = 对齐 y, forecast = fitted
    const auto actual = mm::aligned_y(data, 4);
    ASSERT_EQ(actual.size(), um.fitted.size());
    ASSERT_EQ(actual.size(), dl.fitted.size());
    const auto mz_um = em::mincer_zarnowitz_regression(actual, um.fitted);
    const auto mz_dl = em::mincer_zarnowitz_regression(actual, dl.fitted);
    // 两模型均有预测力 (R² > 0), U-MIDAS 无约束 ≥ DL 约束
    EXPECT_GT(mz_um.r_squared, 0.0);
    EXPECT_GT(mz_dl.r_squared, 0.0);
    EXPECT_GE(mz_um.r_squared, mz_dl.r_squared - 1e-9);

    // (d) DM 检验 (7B 复用): MSE 损失下两模型预测差异
    const auto dm = em::diebold_mariano_test(actual, um.fitted, dl.fitted,
                                             "mse", 1);
    EXPECT_EQ(dm.base.method_name, "Diebold-Mariano");
    EXPECT_TRUE(std::isfinite(dm.base.statistic));
    EXPECT_TRUE(std::isfinite(dm.base.p_value));
    // 嵌套关系 → U-MIDAS 损失不大于 DL (mean_loss_diff ≤ 0)
    EXPECT_LE(dm.mean_loss_diff, 1e-9);
}

// ===========================================================================
// 场景 6: GARCH-M 风险溢价链 (三变体 → λ sandwich 显著性 → vs 无 M 模型)
//   注: gm_baseline 夹具 (T=900, λ=0.5) 的 λ robust t ≈ 0.59 — GARCH-M
//   风险溢价低检验力是文献经典结果; 显著性链路用强信号模拟数据演示
//   (λ=1.2, T=3000, α(1+λ²)+β=0.915<1 平稳), 夹具只做点估计锚定
// ===========================================================================
TEST(Phase7CIntegration, GarchMRiskPremiumChain) {
    const std::vector<Real> y(gmb::Y_GM, gmb::Y_GM + gmb::T);

    // (a) 三变体估计: 均收敛, λ 为正 (DGP: vol 形 λ=0.50)
    for (auto form : {ts::GarchMForm::Volatility, ts::GarchMForm::Variance,
                      ts::GarchMForm::LogVariance}) {
        const auto r = ts::estimate_garch_m(y, form);
        ASSERT_TRUE(r.converged);
        EXPECT_GT(r.params.lambda, 0.0);
        EXPECT_LT(r.params.alpha + r.params.beta, 1.0);   // 方差平稳
    }
    // vol 形主对照: λ ≈ 0.5 (与 DGP 一致, 单测已 1e-4 锚定)
    const auto gm = ts::estimate_garch_m(y, ts::GarchMForm::Volatility);
    EXPECT_NEAR(gm.params.lambda, gmb::ARCH_VOL.lambda, 1e-4);

    // (b) 强信号模拟: r_t = μ + λ√h + √h·z, h_{t+1} = ω + α r² + β h
    //     (Philox 确定性, 预烧 200)
    const auto sim_garch_m = [](Size T, Real mu, Real lam, Real omega,
                                Real alpha, Real beta, uint64_t seed) {
        cpphub::Philox4x64 rng(seed);
        auto u = [&rng]() {
            return (rng() >> 11) * (1.0 / 9007199254740992.0);
        };
        const Size burn = 200;
        std::vector<Real> r(T);
        Real h = omega / (1.0 - beta - alpha * (1.0 + lam * lam));
        for (Size t = 0; t < T + burn; ++t) {
            const auto [z1, z2] = cpphub::box_muller(u(), u());
            (void)z2;
            const Real rt = mu + lam * std::sqrt(h) + std::sqrt(h) * z1;
            h = omega + alpha * rt * rt + beta * h;
            if (t >= burn) r[t - burn] = rt;
        }
        return r;
    };
    constexpr Size kTsim = 3000;
    constexpr Real kLam = 1.2, kOmega = 2e-6, kAlpha = 0.08, kBeta = 0.72;
    const auto ysim = sim_garch_m(kTsim, 0.02, kLam, kOmega, kAlpha, kBeta, 9u);

    // GM4 — λ sandwich 显著性: robust t = λ / se[1] > 1.96
    //     (BW QMLE 三明治, 强信号 + T=3000 → 风险溢价可检出)
    const auto gm_sim =
        ts::estimate_garch_m(ysim, ts::GarchMForm::Volatility);
    ASSERT_TRUE(gm_sim.converged);
    ASSERT_EQ(gm_sim.std_errors.size(), 5u);
    EXPECT_GT(gm_sim.std_errors[1], 0.0);
    const Real t_lambda = gm_sim.params.lambda / gm_sim.std_errors[1];
    EXPECT_GT(t_lambda, 1.96);                 // 5% 水平显著为正
    EXPECT_NEAR(gm_sim.params.lambda, kLam, 0.4);  // 点估计量级恢复
    // vcov 对称 + se = sqrt(diag) (三明治结构自检)
    for (Size a = 0; a < 5; ++a) {
        EXPECT_NEAR(gm_sim.std_errors[a],
                    std::sqrt(gm_sim.vcov[a][a]), 1e-15);
        for (Size b = 0; b < 5; ++b) {
            EXPECT_NEAR(gm_sim.vcov[a][b], gm_sim.vcov[b][a], 1e-10);
        }
    }

    // (c) vs 无 M 模型: 纯 GARCH(1,1) 忽略风险溢价 → 信息准则惩罚
    //     (强信号数据由 GARCH-M 生成, 含 M 模型更接近真 DGP → AIC 更低;
    //     参数数差 1, AIC = −2ℓ + 2k 可比)
    const auto plain_sim = ts::estimate_garch11(ysim);
    ASSERT_TRUE(plain_sim.converged);
    EXPECT_LT(gm_sim.aic, plain_sim.aic);
}
