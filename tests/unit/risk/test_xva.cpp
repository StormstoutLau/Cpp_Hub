// v1.2 Batch 3: xVA (CVA / DVA / FVA / BVA) 单元测试
// 覆盖: PDCurve / Exposure / compute_xva / compute_xva_single_period / compute_xva_mc
#include <gtest/gtest.h>
#include "cpphub/risk/xva.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"
#include <cmath>
#include <vector>
#include <functional>

using namespace cpphub;

// ============================================================
// PDCurve 测试
// ============================================================
TEST(PDCurveTest, FlatHazardRateSurvivalProb) {
    // h = 0.02 (200bp), Q(0, T) = exp(-0.02 * T)
    auto pd = PDCurve::flat(0.02);
    EXPECT_NEAR(pd.survival_prob(1.0), std::exp(-0.02 * 1.0), 1e-12);
    EXPECT_NEAR(pd.survival_prob(5.0), std::exp(-0.02 * 5.0), 1e-12);
    EXPECT_NEAR(pd.survival_prob(10.0), std::exp(-0.02 * 10.0), 1e-12);
    EXPECT_NEAR(pd.survival_prob(0.0), 1.0, 1e-15);
}

TEST(PDCurveTest, FlatHazardRateDefaultProb) {
    auto pd = PDCurve::flat(0.03);
    // PD(0, T) = 1 - exp(-h*T)
    EXPECT_NEAR(pd.default_prob(1.0), 1.0 - std::exp(-0.03), 1e-12);
    EXPECT_NEAR(pd.default_prob(10.0), 1.0 - std::exp(-0.30), 1e-12);
    EXPECT_NEAR(pd.default_prob(0.0), 0.0, 1e-15);
}

TEST(PDCurveTest, IntervalDefaultProb) {
    // PD(t1, t2) = Q(t1) - Q(t2) = exp(-h*t1) - exp(-h*t2)
    auto pd = PDCurve::flat(0.02);
    Real expected = std::exp(-0.02 * 1.0) - std::exp(-0.02 * 5.0);
    EXPECT_NEAR(pd.default_prob(1.0, 5.0), expected, 1e-12);
    // 退化: t1 = t2
    EXPECT_NEAR(pd.default_prob(3.0, 3.0), 0.0, 1e-15);
    // 退化: t1 > t2
    EXPECT_NEAR(pd.default_prob(5.0, 1.0), 0.0, 1e-15);
}

TEST(PDCurveTest, FromCDSSpreadSinglePointApprox) {
    // CDS spread s = 0.02, R = 0.40 → h = s / (1-R) = 0.02/0.6 = 0.0333...
    Real s = 0.02, R = 0.40;
    auto pd = PDCurve::from_cds_spread(s, R);
    Real h_expected = s / (1.0 - R);
    EXPECT_NEAR(pd.hazard_rate(5.0), h_expected, 1e-12);
    EXPECT_NEAR(pd.survival_prob(5.0), std::exp(-h_expected * 5.0), 1e-12);
}

TEST(PDCurveTest, PiecewiseConstantHazardRate) {
    // 分段常数: h=0.01 在 [0.01, 5], h=0.03 在 [5, 10]
    std::vector<Real> times = {1.0, 5.0, 10.0};
    std::vector<Real> hrs = {0.01, 0.03, 0.03};
    PDCurve pd(times, hrs);
    // T=3 在第一段: H(3) = 0.01 * 3
    EXPECT_NEAR(pd.survival_prob(3.0), std::exp(-0.01 * 3.0), 1e-12);
    // T=7 跨两段: H(7) = 0.01*5 + 0.03*(7-5) = 0.05 + 0.06 = 0.11
    EXPECT_NEAR(pd.survival_prob(7.0), std::exp(-0.11), 1e-12);
    // T=5 (节点): H(5) = 0.01*5 = 0.05
    EXPECT_NEAR(pd.survival_prob(5.0), std::exp(-0.05), 1e-12);
}

TEST(PDCurveTest, ZeroHazardRateGivesZeroPD) {
    auto pd = PDCurve::flat(0.0);
    EXPECT_NEAR(pd.default_prob(100.0), 0.0, 1e-15);
    EXPECT_NEAR(pd.survival_prob(100.0), 1.0, 1e-15);
}

TEST(PDCurveTest, InvalidParameters) {
    // 空向量
    EXPECT_THROW(PDCurve({}, {}), std::invalid_argument);
    // 负 hazard rate
    EXPECT_THROW(PDCurve({1.0}, {-0.01}), std::invalid_argument);
    // 非严格递增
    EXPECT_THROW(PDCurve({1.0, 1.0, 2.0}, {0.01, 0.02, 0.03}), std::invalid_argument);
    // size 不匹配
    EXPECT_THROW(PDCurve({1.0, 2.0}, {0.01}), std::invalid_argument);
    // from_cds_spread: R >= 1
    EXPECT_THROW(PDCurve::from_cds_spread(0.02, 1.0), std::invalid_argument);
}

// ============================================================
// Exposure (compute_exposure) 测试
// ============================================================
TEST(ExposureTest, AllPositiveValuesGivesZeroENE) {
    // 2 paths, 3 time points, all V > 0
    std::vector<std::vector<Real>> V = {
        {1.0, 2.0, 3.0},
        {0.5, 1.5, 2.5}
    };
    std::vector<Real> times = {0.0, 0.5, 1.0};
    auto curve = make_flat_curve(0.0);  // P(0,T) = 1
    auto prof = compute_exposure(V, times, curve);
    // EPE = mean of positives = (1.0+0.5)/2 = 0.75 at t=0
    EXPECT_NEAR(prof.epe[0], 0.75, 1e-12);
    EXPECT_NEAR(prof.ene[0], 0.0, 1e-15);
    EXPECT_NEAR(prof.epe[1], 1.75, 1e-12);
    EXPECT_NEAR(prof.ene[1], 0.0, 1e-15);
    EXPECT_NEAR(prof.epe[2], 2.75, 1e-12);
    EXPECT_NEAR(prof.ene[2], 0.0, 1e-15);
}

TEST(ExposureTest, AllNegativeValuesGivesZeroEPE) {
    std::vector<std::vector<Real>> V = {
        {-1.0, -2.0, -3.0},
        {-0.5, -1.5, -2.5}
    };
    std::vector<Real> times = {0.0, 0.5, 1.0};
    auto curve = make_flat_curve(0.0);
    auto prof = compute_exposure(V, times, curve);
    EXPECT_NEAR(prof.epe[0], 0.0, 1e-15);
    EXPECT_NEAR(prof.ene[0], 0.75, 1e-12);
    EXPECT_NEAR(prof.epe[1], 0.0, 1e-15);
    EXPECT_NEAR(prof.ene[1], 1.75, 1e-12);
}

TEST(ExposureTest, MixedValuesDecomposeToEPEandENE) {
    // 1 path, V(t=0) = +2, V(t=1) = -3
    std::vector<std::vector<Real>> V = {{2.0, -3.0}};
    std::vector<Real> times = {0.0, 1.0};
    auto curve = make_flat_curve(0.0);
    auto prof = compute_exposure(V, times, curve);
    EXPECT_NEAR(prof.epe[0], 2.0, 1e-12);
    EXPECT_NEAR(prof.ene[0], 0.0, 1e-15);
    EXPECT_NEAR(prof.epe[1], 0.0, 1e-15);
    EXPECT_NEAR(prof.ene[1], 3.0, 1e-12);
}

TEST(ExposureTest, DiscountFactorApplied) {
    // r = 0.05 → P(0, 1) = exp(-0.05)
    std::vector<std::vector<Real>> V = {{1.0, 1.0}};
    std::vector<Real> times = {0.0, 1.0};
    auto curve = make_flat_curve(0.05);
    auto prof = compute_exposure(V, times, curve);
    EXPECT_NEAR(prof.epe[0], 1.0, 1e-12);  // P(0,0)=1
    EXPECT_NEAR(prof.epe[1], std::exp(-0.05), 1e-12);  // 1 * exp(-0.05)
}

TEST(ExposureTest, SymmetricExposureEPEEqualsENE) {
    // 对称暴露: V = +a 和 -a 各 1 条
    std::vector<std::vector<Real>> V = {
        {2.0, 2.0},
        {-2.0, -2.0}
    };
    std::vector<Real> times = {0.0, 1.0};
    auto curve = make_flat_curve(0.0);
    auto prof = compute_exposure(V, times, curve);
    EXPECT_NEAR(prof.epe[0], prof.ene[0], 1e-12);
    EXPECT_NEAR(prof.epe[1], prof.ene[1], 1e-12);
    EXPECT_NEAR(prof.epe[0], 1.0, 1e-12);  // (2 + 0)/2 = 1
    EXPECT_NEAR(prof.ene[0], 1.0, 1e-12);  // (0 + 2)/2 = 1
}

TEST(ExposureTest, InvalidInputs) {
    auto curve = make_flat_curve(0.0);
    std::vector<std::vector<Real>> empty_V;
    EXPECT_THROW(compute_exposure(empty_V, {0.0, 1.0}, curve), std::invalid_argument);
    // 列数不匹配
    std::vector<std::vector<Real>> V = {{1.0, 2.0}, {1.0}};  // 第二条路径列数不对
    EXPECT_THROW(compute_exposure(V, {0.0, 1.0}, curve), std::invalid_argument);
}

// ============================================================
// compute_xva 单期测试 (trapezoidal 积分 → 解析公式匹配)
// ============================================================
TEST(XVATest, CVASignAndFormula) {
    // 单期 EPE=10, ENE=0, PD_c=0.05, R_c=0.40 → CVA = -0.6 * 10 * 0.05 = -0.30
    // 构造 exposure profile: 2 个时间点, t=0 和 t=T
    ExposureProfile prof;
    prof.times = {0.0, 1.0};
    prof.epe = {10.0, 10.0};  // 梯形 = 10
    prof.ene = {0.0, 0.0};
    auto pd_c = PDCurve::flat(0.05);   // h=0.05, T=1 → PD ≈ 1 - exp(-0.05) ≈ 0.04877
    auto pd_s = PDCurve::flat(0.0);    // 自身不违约 → DVA = 0
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.0;  // 关闭 FVA
    auto r = compute_xva(prof, pd_c, pd_s, cfg, 100.0);
    Real expected_pd_c = 1.0 - std::exp(-0.05);
    Real expected_cva = -(1.0 - 0.40) * 10.0 * expected_pd_c;
    EXPECT_NEAR(r.cva, expected_cva, 1e-10);
    EXPECT_LT(r.cva, 0.0);  // CVA ≤ 0
    EXPECT_NEAR(r.dva, 0.0, 1e-12);  // PD_self = 0
    EXPECT_NEAR(r.fva, 0.0, 1e-12);  // s_f = 0
    EXPECT_NEAR(r.bva, r.cva, 1e-12);
    EXPECT_NEAR(r.adjusted_price, 100.0 + r.cva, 1e-12);
}

TEST(XVATest, DVASignAndFormula) {
    // 对称场景: EPE=0, ENE=10, PD_self=0.05, R_self=0.40
    ExposureProfile prof;
    prof.times = {0.0, 1.0};
    prof.epe = {0.0, 0.0};
    prof.ene = {10.0, 10.0};
    auto pd_c = PDCurve::flat(0.0);    // 对手不违约 → CVA = 0
    auto pd_s = PDCurve::flat(0.05);   // 自身违约
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.0;
    auto r = compute_xva(prof, pd_c, pd_s, cfg);
    Real expected_pd_s = 1.0 - std::exp(-0.05);
    Real expected_dva = (1.0 - 0.40) * 10.0 * expected_pd_s;
    EXPECT_NEAR(r.dva, expected_dva, 1e-10);
    EXPECT_GT(r.dva, 0.0);  // DVA ≥ 0
    EXPECT_NEAR(r.cva, 0.0, 1e-12);
    EXPECT_NEAR(r.bva, r.dva, 1e-12);
}

TEST(XVATest, FVASymmetricExposure) {
    // EPE=10, ENE=0, s_f=0.01, T=1 → FVA = -0.01 * 10 * 1 = -0.10
    ExposureProfile prof;
    prof.times = {0.0, 1.0};
    prof.epe = {10.0, 10.0};
    prof.ene = {0.0, 0.0};
    auto pd_zero = PDCurve::flat(0.0);  // 关闭 CVA/DVA
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.01;
    auto r = compute_xva(prof, pd_zero, pd_zero, cfg);
    EXPECT_NEAR(r.cva, 0.0, 1e-12);
    EXPECT_NEAR(r.dva, 0.0, 1e-12);
    EXPECT_NEAR(r.fva, -0.10, 1e-10);
    EXPECT_LT(r.fva, 0.0);  // 正暴露 → funding cost → FVA < 0
}

TEST(XVATest, FVASymmetricExposureNegative) {
    // EPE=0, ENE=10 → net exposure 为负, FVA = -0.01 * (0 - 10) * 1 = +0.10
    ExposureProfile prof;
    prof.times = {0.0, 1.0};
    prof.epe = {0.0, 0.0};
    prof.ene = {10.0, 10.0};
    auto pd_zero = PDCurve::flat(0.0);
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.01;
    auto r = compute_xva(prof, pd_zero, pd_zero, cfg);
    EXPECT_NEAR(r.fva, 0.10, 1e-10);
    EXPECT_GT(r.fva, 0.0);  // 负暴露 → funding benefit → FVA > 0
}

TEST(XVATest, FVASymmetricExposureZeroWhenEPEEqualsENE) {
    // 对称暴露 EPE=ENE=5 → net = 0, FVA = 0
    ExposureProfile prof;
    prof.times = {0.0, 1.0};
    prof.epe = {5.0, 5.0};
    prof.ene = {5.0, 5.0};
    auto pd_zero = PDCurve::flat(0.0);
    XVAConfig cfg;
    cfg.funding_spread = 0.01;
    auto r = compute_xva(prof, pd_zero, pd_zero, cfg);
    EXPECT_NEAR(r.fva, 0.0, 1e-12);
}

TEST(XVATest, BVASumCVAplusDVAplusFVA) {
    // 混合场景: 全部非零
    ExposureProfile prof;
    prof.times = {0.0, 1.0};
    prof.epe = {10.0, 10.0};
    prof.ene = {5.0, 5.0};
    auto pd_c = PDCurve::flat(0.03);
    auto pd_s = PDCurve::flat(0.02);
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.005;
    auto r = compute_xva(prof, pd_c, pd_s, cfg, 50.0);
    EXPECT_NEAR(r.bva, r.cva + r.dva + r.fva, 1e-12);
    EXPECT_NEAR(r.adjusted_price, 50.0 + r.bva, 1e-12);
}

TEST(XVATest, SymmetryCVAEqualsDVAWhenSwapped) {
    // 当 PD_c = PD_s, R_c = R_s, EPE ↔ ENE 互换 → CVA ↔ DVA (符号相反, 绝对值相等)
    // 数学: CVA = -L * EPE * PD (≤0), DVA = +L * ENE * PD (≥0)
    // 所以 r1.cva = -r2.dva, r1.dva = -r2.cva
    ExposureProfile prof_epe;
    prof_epe.times = {0.0, 1.0};
    prof_epe.epe = {10.0, 10.0};
    prof_epe.ene = {0.0, 0.0};
    ExposureProfile prof_ene;
    prof_ene.times = {0.0, 1.0};
    prof_ene.epe = {0.0, 0.0};
    prof_ene.ene = {10.0, 10.0};
    auto pd = PDCurve::flat(0.05);
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.0;
    auto r1 = compute_xva(prof_epe, pd, pd, cfg);
    auto r2 = compute_xva(prof_ene, pd, pd, cfg);
    // r1: CVA<0, DVA=0; r2: CVA=0, DVA>0; r1.cva = -r2.dva
    EXPECT_NEAR(r1.cva, -r2.dva, 1e-12);
    EXPECT_NEAR(r1.dva, -r2.cva, 1e-12);
    EXPECT_NEAR(std::abs(r1.cva), std::abs(r2.dva), 1e-12);
}

// ============================================================
// compute_xva_single_period 与 compute_xva 一致性
// ============================================================
TEST(XVATest, SinglePeriodConsistentWithMultiPeriod) {
    // 单期 EPE=10, ENE=3, T=1, PD_c=0.05, PD_s=0.03, R=0.4, s_f=0.01
    // compute_xva_single_period:
    //   CVA = -0.6 * 10 * 0.05 = -0.30
    //   DVA = +0.6 * 3 * 0.03 = +0.054
    //   FVA = -0.01 * (10 - 3) * 1 = -0.07
    Real epe = 10.0, ene = 3.0, T = 1.0;
    Real pd_c_val = 0.05, pd_s_val = 0.03;
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.01;
    auto r_sp = compute_xva_single_period(epe, ene, T, pd_c_val, pd_s_val, cfg, 100.0);
    EXPECT_NEAR(r_sp.cva, -0.30, 1e-12);
    EXPECT_NEAR(r_sp.dva, 0.054, 1e-12);
    EXPECT_NEAR(r_sp.fva, -0.07, 1e-12);
    EXPECT_NEAR(r_sp.bva, -0.30 + 0.054 - 0.07, 1e-12);
    EXPECT_NEAR(r_sp.adjusted_price, 100.0 + r_sp.bva, 1e-12);
}

TEST(XVATest, SinglePeriodMatchesDiscreteOneStep) {
    // compute_xva with single step [0, T] vs compute_xva_single_period
    // 注意: compute_xva 使用 hazard rate PD, single_period 接收的是 PD 数值
    // 取 hazard = 0.05, T = 1 → PD(0,1) = 1 - exp(-0.05) ≈ 0.04877
    Real T = 1.0;
    Real h_c = 0.05, h_s = 0.03;
    Real pd_c_val = 1.0 - std::exp(-h_c * T);
    Real pd_s_val = 1.0 - std::exp(-h_s * T);
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.01;

    ExposureProfile prof;
    prof.times = {0.0, T};
    prof.epe = {10.0, 10.0};
    prof.ene = {3.0, 3.0};
    auto pd_c_curve = PDCurve::flat(h_c);
    auto pd_s_curve = PDCurve::flat(h_s);
    auto r_multi = compute_xva(prof, pd_c_curve, pd_s_curve, cfg, 100.0);
    auto r_sp = compute_xva_single_period(10.0, 3.0, T, pd_c_val, pd_s_val, cfg, 100.0);
    // 梯形公式对单步等价于端点平均 → EPE=ENE 恒定 → 与 single_period 精确一致
    EXPECT_NEAR(r_multi.cva, r_sp.cva, 1e-10);
    EXPECT_NEAR(r_multi.dva, r_sp.dva, 1e-10);
    // FVA: compute_xva 用 dt=1, single_period 用 T=1 → 一致
    EXPECT_NEAR(r_multi.fva, r_sp.fva, 1e-10);
}

// ============================================================
// 多期 trapezoidal 积分测试
// ============================================================
TEST(XVATest, MultiPeriodTrapezoidalIntegration) {
    // 4 个时间点 {0, 1, 2, 3}, EPE 线性增长
    // trapezoidal: 0.5*(epe[i]+epe[i+1]) * dt
    ExposureProfile prof;
    prof.times = {0.0, 1.0, 2.0, 3.0};
    prof.epe = {0.0, 1.0, 2.0, 3.0};  // 每段平均 = 0.5, 1.5, 2.5
    prof.ene = {0.0, 0.0, 0.0, 0.0};
    auto pd_c = PDCurve::flat(0.05);
    auto pd_s = PDCurve::flat(0.0);
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.0;
    auto r = compute_xva(prof, pd_c, pd_s, cfg);
    // CVA = -0.6 * Σ 0.5*(epe_i+epe_{i+1}) * ΔPD_i
    Real expected_cva = 0.0;
    for (Size i = 0; i + 1 < 4; ++i) {
        Real avg_epe = 0.5 * (prof.epe[i] + prof.epe[i + 1]);
        Real dPD = std::exp(-0.05 * prof.times[i]) - std::exp(-0.05 * prof.times[i + 1]);
        expected_cva += -0.6 * avg_epe * dPD;
    }
    EXPECT_NEAR(r.cva, expected_cva, 1e-10);
}

TEST(XVATest, ProfileTooFewPoints) {
    ExposureProfile prof;
    prof.times = {0.0};
    prof.epe = {1.0};
    prof.ene = {0.0};
    auto pd = PDCurve::flat(0.05);
    EXPECT_THROW(compute_xva(prof, pd, pd, XVAConfig{}), std::invalid_argument);
}

// ============================================================
// MC 端到端: 用 vanilla call option 验证
// ============================================================
TEST(XVAMCTest, VanillaCallCVAFromPathSim) {
    // 场景: bank 持有 call option (long call), V(t) = BS call price at time t
    // 用 MC 模拟 GBM 路径, 在每个 t_i 计算 call price, 然后计算 CVA
    // 验证: CVA < 0 (正暴露)
    Real S0 = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real T = 1.0;
    Size n_steps = 50;
    auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
    MultiAssetGBMPathGenerator gen(gen_cfg);

    // value_fn: 在 t 时刻 call 的 BS 价格
    // 剩余期限 tau = T - t, 用 BS 公式
    auto value_fn = [K, r, q, sigma, T](Real t, const std::vector<Real>& S_t) -> Real {
        Real tau = T - t;
        if (tau <= 0.0) {
            // 到期: payoff
            return std::max(S_t[0] - K, 0.0);
        }
        return AnalyticGreeksEngine::bsm_european(S_t[0], K, tau, r, q, sigma, true).price;
    };

    std::vector<Real> exposure_times = {0.0, 0.25, 0.5, 0.75, 1.0};
    auto pd_c = PDCurve::flat(0.05);   // 5% hazard
    auto pd_s = PDCurve::flat(0.0);    // bank 不违约
    auto disc = make_flat_curve(r);    // 折现曲线
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.0;

    Size n_paths = 5000;
    auto r_xva = compute_xva_mc(gen, value_fn, exposure_times,
                                  pd_c, pd_s, disc, cfg, n_paths, 42, 0.0);

    // CVA 应为负
    EXPECT_LT(r_xva.cva, 0.0);
    // DVA = 0 (PD_self = 0)
    EXPECT_NEAR(r_xva.dva, 0.0, 1e-12);
    // 验证量级: PD(0,1) ≈ 0.04877, LGD = 0.6
    // 平均 EPE_disc 应在 5-15 之间 (call price ~ 10)
    Real pd_c_val = 1.0 - std::exp(-0.05 * T);
    // 粗略估计: CVA ≈ -0.6 * avg_EPE * 0.04877
    // 反推 avg_EPE ≈ -CVA / (0.6 * 0.04877), 应在 5-15 范围
    Real implied_epe = -r_xva.cva / (0.6 * pd_c_val);
    EXPECT_GT(implied_epe, 3.0);
    EXPECT_LT(implied_epe, 15.0);
}

TEST(XVAMCTest, ShortCallPositionGivesDVA) {
    // bank short call (卖出 call): V(t) = -BS call price → EPE=0, ENE>0
    // 此时 CVA = 0, DVA > 0 (bank 自身违约时获益)
    Real S0 = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real T = 1.0;
    Size n_steps = 50;
    auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
    MultiAssetGBMPathGenerator gen(gen_cfg);

    // short call: value = -BS_price
    auto value_fn = [K, r, q, sigma, T](Real t, const std::vector<Real>& S_t) -> Real {
        Real tau = T - t;
        if (tau <= 0.0) return -std::max(S_t[0] - K, 0.0);
        return -AnalyticGreeksEngine::bsm_european(S_t[0], K, tau, r, q, sigma, true).price;
    };

    std::vector<Real> exposure_times = {0.0, 0.5, 1.0};
    auto pd_c = PDCurve::flat(0.0);    // 对手不违约
    auto pd_s = PDCurve::flat(0.04);   // bank 违约
    auto disc = make_flat_curve(r);
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.0;

    Size n_paths = 5000;
    auto r_xva = compute_xva_mc(gen, value_fn, exposure_times,
                                  pd_c, pd_s, disc, cfg, n_paths, 42);
    EXPECT_NEAR(r_xva.cva, 0.0, 1e-12);  // PD_c = 0
    EXPECT_GT(r_xva.dva, 0.0);  // ENE > 0, PD_s > 0
}

TEST(XVAMCTest, AdjustedPriceLowerThanRiskFreeForLongCall) {
    // long call + CVA → adjusted < risk_free
    Real S0 = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real T = 1.0;
    Size n_steps = 50;
    auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
    MultiAssetGBMPathGenerator gen(gen_cfg);

    auto value_fn = [K, r, q, sigma, T](Real t, const std::vector<Real>& S_t) -> Real {
        Real tau = T - t;
        if (tau <= 0.0) return std::max(S_t[0] - K, 0.0);
        return AnalyticGreeksEngine::bsm_european(S_t[0], K, tau, r, q, sigma, true).price;
    };

    std::vector<Real> exposure_times = {0.0, 0.5, 1.0};
    auto pd_c = PDCurve::flat(0.05);
    auto pd_s = PDCurve::flat(0.0);
    auto disc = make_flat_curve(r);
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.0;

    // risk-free price = BS price at t=0
    Real rf_price = AnalyticGreeksEngine::bsm_european(S0, K, T, r, q, sigma, true).price;
    Size n_paths = 5000;
    auto r_xva = compute_xva_mc(gen, value_fn, exposure_times,
                                  pd_c, pd_s, disc, cfg, n_paths, 42, rf_price);
    EXPECT_LT(r_xva.adjusted_price, rf_price);
    EXPECT_GT(r_xva.adjusted_price, 0.0);  // 仍为正
}

TEST(XVAMCTest, ExposureTimeOutOfRangeRejected) {
    Real S0 = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real T = 1.0;
    Size n_steps = 50;
    auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
    MultiAssetGBMPathGenerator gen(gen_cfg);
    auto value_fn = [](Real, const std::vector<Real>&) { return 0.0; };
    std::vector<Real> bad_times = {0.0, 2.0};  // 2.0 > T
    auto pd = PDCurve::flat(0.05);
    auto disc = make_flat_curve(r);
    EXPECT_THROW(compute_xva_mc(gen, value_fn, bad_times, pd, pd, disc, XVAConfig{}, 100, 42),
                  std::invalid_argument);
}

TEST(XVAMCTest, ZeroPathsRejected) {
    Real S0 = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real T = 1.0;
    Size n_steps = 50;
    auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
    MultiAssetGBMPathGenerator gen(gen_cfg);
    auto value_fn = [](Real, const std::vector<Real>&) { return 0.0; };
    auto pd = PDCurve::flat(0.05);
    auto disc = make_flat_curve(r);
    EXPECT_THROW(compute_xva_mc(gen, value_fn, {0.0, 1.0}, pd, pd, disc, XVAConfig{}, 0, 42),
                  std::invalid_argument);
}

// ============================================================
// FVA 双向符号测试 (CVA/DVA 关闭, 仅 FVA)
// ============================================================
TEST(XVATest, FVAOnlyPureFundingCost) {
    // EPE=10, ENE=0, s_f=0.02, T=2 → FVA = -0.02 * 10 * 2 = -0.40
    ExposureProfile prof;
    prof.times = {0.0, 1.0, 2.0};
    prof.epe = {10.0, 10.0, 10.0};
    prof.ene = {0.0, 0.0, 0.0};
    auto pd_zero = PDCurve::flat(0.0);
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.02;
    auto r = compute_xva(prof, pd_zero, pd_zero, cfg);
    EXPECT_NEAR(r.cva, 0.0, 1e-12);
    EXPECT_NEAR(r.dva, 0.0, 1e-12);
    // 两段 [0,1] 和 [1,2], 每段 avg_epe=10, dt=1, s_f=0.02
    // FVA = -0.02 * 10 * 1 + -0.02 * 10 * 1 = -0.40
    EXPECT_NEAR(r.fva, -0.40, 1e-10);
}

// ============================================================
// 综合: 跨产品场景 (long call + short put = long forward)
// 验证: forward 的暴露对称, EPE ≈ ENE
// ============================================================
TEST(XVAMCTest, LongForwardSymmetricExposure) {
    // long forward = long call (K) + short put (K), synthetic forward
    // V(t) = S(t) - K * exp(-r*(T-t))  (forward value at t)
    // 对 ATM forward, 暴露近似对称
    Real S0 = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real T = 1.0;
    Size n_steps = 50;
    auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
    MultiAssetGBMPathGenerator gen(gen_cfg);

    // forward value at t: V(t) = S(t) - K * exp(-r*(T-t))
    auto value_fn = [K, r, T](Real t, const std::vector<Real>& S_t) -> Real {
        Real tau = T - t;
        return S_t[0] - K * std::exp(-r * tau);
    };

    std::vector<Real> exposure_times = {0.0, 0.5, 1.0};
    auto pd_c = PDCurve::flat(0.05);
    auto pd_s = PDCurve::flat(0.05);
    auto disc = make_flat_curve(r);
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.005;

    Size n_paths = 10000;
    auto r_xva = compute_xva_mc(gen, value_fn, exposure_times,
                                  pd_c, pd_s, disc, cfg, n_paths, 42);
    // forward 暴露对称 → CVA ≈ -DVA (符号相反, 量级接近)
    // 注: GBM 下 forward 不完全对称 (lognormal 偏度), 但应接近
    EXPECT_LT(r_xva.cva, 0.0);
    EXPECT_GT(r_xva.dva, 0.0);
    // CVA + DVA 应相对小 (近似抵消), 主要残留来自 lognormal 偏度
    Real net_credit_adj = r_xva.cva + r_xva.dva;
    Real gross = std::abs(r_xva.cva) + r_xva.dva;
    // net / gross 应远小于 1 (强抵消)
    if (gross > 1e-6) {
        EXPECT_LT(std::abs(net_credit_adj) / gross, 0.5);
    }
}
