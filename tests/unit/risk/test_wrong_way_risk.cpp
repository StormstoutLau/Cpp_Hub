// v1.3 D: CVA with Wrong-Way Risk (WWR) 单元测试
// 覆盖: compute_cva_wwr_hw / compute_cva_wwr_copula / compute_cva_wwr_mc / estimate_exposure_volatility
#include <gtest/gtest.h>
#include "cpphub/risk/wrong_way_risk.hpp"
#include "cpphub/risk/xva.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include "cpphub/instruments/credit/credit_curve.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/core/math.hpp"
#include <cmath>
#include <vector>
#include <functional>
#include <limits>

using namespace cpphub;

namespace {

// 平坦暴露轮廓: EPE 常数, ENE=0, 时间网格覆盖 [0, T]
ExposureProfile flat_profile(Real T, Real epe, Size n = 5) {
    ExposureProfile p;
    for (Size i = 0; i < n; ++i) {
        p.times.push_back(T * static_cast<Real>(i) / static_cast<Real>(n - 1));
        p.epe.push_back(epe);
        p.ene.push_back(0.0);
    }
    return p;
}

XVAConfig default_cfg() {
    XVAConfig cfg;
    cfg.recovery_counterparty = 0.40;
    cfg.recovery_self = 0.40;
    cfg.funding_spread = 0.0;
    return cfg;
}

}  // namespace

// ============================================================
// Hull-White 近似: 定性行为 (测试要求 1-7)
// ============================================================
TEST(WrongWayRiskHWTest, ZeroRhoDegeneratesToIndependentCVA) {
    // 要求 1: rho_ww=0 时 WWR CVA == 独立 CVA (容差 1e-10)
    auto prof = flat_profile(1.0, 10.0);
    auto pd_c = PDCurve::flat(0.05);
    auto pd_s = PDCurve::flat(0.0);
    auto cfg = default_cfg();
    auto base = compute_xva(prof, pd_c, pd_s, cfg);
    WWRConfig wwr;
    wwr.rho_ww = 0.0;
    auto r = compute_cva_wwr_hw(prof, pd_c, pd_s, cfg, wwr, 0.3);
    EXPECT_NEAR(r.cva, base.cva, 1e-10);
    EXPECT_NEAR(r.dva, base.dva, 1e-12);
    EXPECT_NEAR(r.fva, base.fva, 1e-12);
    EXPECT_NEAR(r.bva, base.bva, 1e-10);
    EXPECT_NEAR(r.adjusted_price, base.adjusted_price, 1e-10);
}

TEST(WrongWayRiskHWTest, PositiveRhoAmplifiesCVA) {
    // 要求 2: rho_ww>0 → WWR 放大 (|CVA| > 独立 CVA)
    auto prof = flat_profile(1.0, 10.0);
    auto pd_c = PDCurve::flat(2.0);  // PD(0,1) = 1-e^{-2} = 0.865 > 0.5
    auto pd_s = PDCurve::flat(0.0);
    auto cfg = default_cfg();
    auto base = compute_xva(prof, pd_c, pd_s, cfg);
    WWRConfig wwr;
    wwr.rho_ww = 0.5;
    auto r = compute_cva_wwr_hw(prof, pd_c, pd_s, cfg, wwr, 0.6);
    EXPECT_LT(r.cva, 0.0);
    EXPECT_LT(base.cva, 0.0);
    EXPECT_LT(r.cva, base.cva);  // 更负
    EXPECT_GT(std::abs(r.cva), std::abs(base.cva));
}

TEST(WrongWayRiskHWTest, NegativeRhoReducesCVA) {
    // 要求 3: rho_ww<0 → RWR 缩小 (|CVA| < 独立 CVA)
    auto prof = flat_profile(1.0, 10.0);
    auto pd_c = PDCurve::flat(2.0);
    auto pd_s = PDCurve::flat(0.0);
    auto cfg = default_cfg();
    auto base = compute_xva(prof, pd_c, pd_s, cfg);
    WWRConfig wwr;
    wwr.rho_ww = -0.5;
    auto r = compute_cva_wwr_hw(prof, pd_c, pd_s, cfg, wwr, 0.6);
    EXPECT_LT(r.cva, 0.0);
    EXPECT_GT(r.cva, base.cva);  // 更接近 0
    EXPECT_LT(std::abs(r.cva), std::abs(base.cva));
}

TEST(WrongWayRiskHWTest, MonotonicInRho) {
    // 要求 4: |CVA| 随 rho_ww 单调递增
    auto prof = flat_profile(1.0, 10.0);
    auto pd_c = PDCurve::flat(2.0);  // PD_T=0.865, Φ^{-1}(PD_T)=1.10
    auto pd_s = PDCurve::flat(0.0);
    auto cfg = default_cfg();
    const std::vector<Real> rhos = {-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9};
    Real prev = 0.0;
    for (Real rho : rhos) {
        WWRConfig wwr;
        wwr.rho_ww = rho;
        auto r = compute_cva_wwr_hw(prof, pd_c, pd_s, cfg, wwr, 0.6);
        const Real mag = std::abs(r.cva);
        EXPECT_GE(mag, prev - 1e-12);
        prev = mag;
    }
}

TEST(WrongWayRiskHWTest, ExtremeRhoAmplifiesAtLeast2x) {
    // 要求 5: 极端 rho_ww→1 → CVA 显著放大 (≥ 2x 独立 CVA)
    auto prof = flat_profile(1.0, 10.0);
    auto pd_c = PDCurve::flat(2.0);
    auto pd_s = PDCurve::flat(0.0);
    auto cfg = default_cfg();
    auto base = compute_xva(prof, pd_c, pd_s, cfg);
    WWRConfig wwr;
    wwr.rho_ww = 0.99;
    auto r = compute_cva_wwr_hw(prof, pd_c, pd_s, cfg, wwr, 0.6);
    EXPECT_GE(std::abs(r.cva), 2.0 * std::abs(base.cva));
}

TEST(WrongWayRiskHWTest, ZeroPDGivesZeroCVA) {
    // 要求 6: PD=0 边界 → CVA = 0 (HW + Copula)
    auto prof = flat_profile(1.0, 10.0);
    auto pd_zero = PDCurve::flat(0.0);
    auto pd_s = PDCurve::flat(0.0);
    auto cfg = default_cfg();
    WWRConfig wwr;
    wwr.rho_ww = 0.7;
    auto hw = compute_cva_wwr_hw(prof, pd_zero, pd_s, cfg, wwr, 0.6);
    auto cop = compute_cva_wwr_copula(prof, pd_zero, pd_s, cfg, wwr, 0.6);
    EXPECT_NEAR(hw.cva, 0.0, 1e-12);
    EXPECT_NEAR(cop.cva, 0.0, 1e-12);
    EXPECT_NEAR(hw.bva, 0.0, 1e-12);
    EXPECT_NEAR(cop.bva, 0.0, 1e-12);
}

TEST(WrongWayRiskHWTest, FullRecoveryGivesZeroCVA) {
    // 要求 7: R=1 边界 → CVA = 0 (无损失, HW + Copula)
    auto prof = flat_profile(1.0, 10.0);
    auto pd_c = PDCurve::flat(0.05);
    auto pd_s = PDCurve::flat(0.0);
    XVAConfig cfg = default_cfg();
    cfg.recovery_counterparty = 1.0;  // 全回收 → LGD=0
    WWRConfig wwr;
    wwr.rho_ww = 0.7;
    auto hw = compute_cva_wwr_hw(prof, pd_c, pd_s, cfg, wwr, 0.6);
    auto cop = compute_cva_wwr_copula(prof, pd_c, pd_s, cfg, wwr, 0.6);
    EXPECT_NEAR(hw.cva, 0.0, 1e-12);
    EXPECT_NEAR(cop.cva, 0.0, 1e-12);
}

// ============================================================
// Copula vs HW 一致性 (要求 8)
// ============================================================
TEST(WrongWayRiskConsistencyTest, CopulaMatchesHWAtSmallRho) {
    // 小 rho (≤0.3) 时 Copula 与 HW 近似结果接近 (容差 10%)
    auto prof = flat_profile(1.0, 10.0);
    auto pd_c = PDCurve::flat(0.05);
    auto pd_s = PDCurve::flat(0.0);
    auto cfg = default_cfg();
    WWRConfig wwr;
    wwr.rho_ww = 0.2;
    const Real sigma_V = 0.5;
    auto hw = compute_cva_wwr_hw(prof, pd_c, pd_s, cfg, wwr, sigma_V);
    auto cop = compute_cva_wwr_copula(prof, pd_c, pd_s, cfg, wwr, sigma_V);
    EXPECT_LT(std::abs(hw.cva - cop.cva), 0.10 * std::abs(hw.cva));
}

TEST(WrongWayRiskConsistencyTest, CopulaZeroRhoDegeneratesToIndependent) {
    // 附加: Copula rho=0 → 退化为独立 CVA (容差 2%)
    auto prof = flat_profile(1.0, 10.0);
    auto pd_c = PDCurve::flat(0.05);
    auto pd_s = PDCurve::flat(0.0);
    auto cfg = default_cfg();
    auto base = compute_xva(prof, pd_c, pd_s, cfg);
    WWRConfig wwr;
    wwr.rho_ww = 0.0;
    auto cop = compute_cva_wwr_copula(prof, pd_c, pd_s, cfg, wwr, 0.5);
    EXPECT_NEAR(cop.cva, base.cva, 0.02 * std::abs(base.cva));
}

// ============================================================
// MC WWR (要求 9, 10)
// ============================================================
TEST(WrongWayRiskMCTest, MCWWRMatchesCopula) {
    // 要求 9: MC 与半解析 Copula 一致 (容差 5%, n_paths>=10000)
    const Real S0 = 100.0, sigma = 0.20, r = 0.05, q = 0.02, T = 1.0;
    const Size n_steps = 50;
    auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
    MultiAssetGBMPathGenerator gen(gen_cfg);

    // 暴露 = 持有 1 单位标的: V(t) = S(t) (恰好对数正态)
    auto value_fn = [](Real, const std::vector<Real>& S) -> Real { return S[0]; };
    std::vector<Real> exposure_times = {0.0, 0.25, 0.5, 0.75, 1.0};
    auto pd_c = PDCurve::flat(0.10);
    auto pd_s = PDCurve::flat(0.0);
    auto disc = make_flat_curve(r);
    auto cfg = default_cfg();
    WWRConfig wwr;
    wwr.rho_ww = 0.5;

    const Size n_paths = 20000;
    const uint64_t seed = 42;
    const Real mc_cva = compute_cva_wwr_mc(gen, value_fn, exposure_times, pd_c, disc, cfg, wwr, sigma, n_paths, seed);

    // 解析暴露轮廓: E[max(V,0)] = E[S(t)] = S0*e^{(r-q)t}, EPE_disc = E[S]*P_d
    ExposureProfile prof;
    prof.times = exposure_times;
    for (Real t : exposure_times) {
        prof.epe.push_back(S0 * std::exp((r - q) * t) * disc.discount_factor(t));
        prof.ene.push_back(0.0);
    }
    auto cop = compute_cva_wwr_copula(prof, pd_c, pd_s, cfg, wwr, sigma);
    EXPECT_LT(cop.cva, 0.0);
    EXPECT_LT(mc_cva, 0.0);
    EXPECT_NEAR(mc_cva, cop.cva, 0.05 * std::abs(cop.cva));
}

TEST(WrongWayRiskMCTest, MCZeroRhoDegeneratesToIndependent) {
    // 要求 10: MC WWR 在 rho_ww=0 时退化为独立 MC CVA
    const Real S0 = 100.0, sigma = 0.20, r = 0.05, q = 0.02, T = 1.0;
    const Size n_steps = 50;
    auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
    MultiAssetGBMPathGenerator gen(gen_cfg);
    auto value_fn = [](Real, const std::vector<Real>& S) -> Real { return S[0]; };
    std::vector<Real> exposure_times = {0.0, 0.25, 0.5, 0.75, 1.0};
    auto pd_c = PDCurve::flat(0.10);
    auto pd_s = PDCurve::flat(0.0);
    auto disc = make_flat_curve(r);
    auto cfg = default_cfg();
    WWRConfig wwr;
    wwr.rho_ww = 0.0;

    const Size n_paths = 20000;
    const uint64_t seed = 42;
    const Real mc_wwr = compute_cva_wwr_mc(gen, value_fn, exposure_times, pd_c, disc, cfg, wwr, sigma, n_paths, seed);
    auto r_xva = compute_xva_mc(gen, value_fn, exposure_times, pd_c, pd_s, disc, cfg, n_paths, seed);
    EXPECT_LT(r_xva.cva, 0.0);
    EXPECT_NEAR(mc_wwr, r_xva.cva, 0.10 * std::abs(r_xva.cva));
}

// ============================================================
// 暴露波动率估计 (要求 11)
// ============================================================
TEST(WrongWayRiskVolTest, ExposureVolatilityMatchesKnownGBM) {
    // V(t) = exp(sigma_true * sqrt(t) * Z) → 解析波动率 = sigma_true
    const Real sigma_true = 0.30;
    const Size n_paths = 50000;
    std::vector<Real> times = {0.5, 1.0, 2.0};
    std::vector<std::vector<Real>> samples(n_paths, std::vector<Real>(times.size()));
    for (Size p = 0; p < n_paths; ++p) {
        Philox4x64 rng(42, static_cast<uint64_t>(p));
        for (Size j = 0; j < times.size(); ++j) {
            samples[p][j] = std::exp(sigma_true * std::sqrt(times[j]) * next_normal(rng));
        }
    }
    const Real est = estimate_exposure_volatility(samples, times);
    EXPECT_NEAR(est, sigma_true, 0.02);
}

TEST(WrongWayRiskVolTest, ExposureVolatilityZeroForDeterministic) {
    // 确定性暴露 → 波动率 0
    const Size n_paths = 100;
    std::vector<Real> times = {1.0, 2.0};
    std::vector<std::vector<Real>> samples(n_paths, std::vector<Real>{5.0, 5.0});
    const Real est = estimate_exposure_volatility(samples, times);
    EXPECT_NEAR(est, 0.0, 1e-12);
}

// ============================================================
// 数值稳定性 (要求 12)
// ============================================================
TEST(WrongWayRiskStabilityTest, LargeRhoLongMaturityNoNaN) {
    // 大 rho (0.9) + 长期限 (10y) 不产生 NaN/Inf
    const Real T = 10.0;
    auto prof = flat_profile(T, 10.0, 11);
    auto pd_c = PDCurve::flat(0.05);
    auto pd_s = PDCurve::flat(0.0);
    auto cfg = default_cfg();
    WWRConfig wwr;
    wwr.rho_ww = 0.9;
    const Real sigma_V = 0.3;

    auto hw = compute_cva_wwr_hw(prof, pd_c, pd_s, cfg, wwr, sigma_V);
    auto cop = compute_cva_wwr_copula(prof, pd_c, pd_s, cfg, wwr, sigma_V);
    EXPECT_TRUE(std::isfinite(hw.cva));
    EXPECT_TRUE(std::isfinite(hw.bva));
    EXPECT_TRUE(std::isfinite(cop.cva));
    EXPECT_TRUE(std::isfinite(cop.bva));

    const Real S0 = 100.0, sigma = 0.20, r = 0.05, q = 0.0;
    const Size n_steps = 200;
    auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
    MultiAssetGBMPathGenerator gen(gen_cfg);
    auto value_fn = [](Real, const std::vector<Real>& S) -> Real { return S[0]; };
    std::vector<Real> exposure_times;
    for (Size i = 0; i <= 10; ++i) exposure_times.push_back(static_cast<Real>(i));
    auto disc = make_flat_curve(r);
    const Real mc = compute_cva_wwr_mc(gen, value_fn, exposure_times, pd_c, disc, cfg, wwr, sigma, 5000, 42);
    EXPECT_TRUE(std::isfinite(mc));
    EXPECT_LT(mc, 0.0);
}
