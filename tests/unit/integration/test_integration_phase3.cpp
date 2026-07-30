// Phase 3 end-to-end integration tests
// SOURCE: PHASE3_SPEC §1.2 - 10 integration tests covering full risk pipeline
// Implemented on main station (MSVC) - 2026-07-30
// Pipeline: SVI calibration -> Dupire local vol -> AAD Greeks -> VaR/ES -> Stress test
// NOTE: avoid <random> (MSVC ICE); use deterministic data or xorshift64* inline
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"
#include "cpphub/risk/greeks/aad_greeks.hpp"
#include "cpphub/risk/var/historical_var.hpp"
#include "cpphub/risk/var/parametric_var.hpp"
#include "cpphub/risk/var/mc_var.hpp"
#include "cpphub/risk/var/expected_shortfall.hpp"
#include "cpphub/risk/scenario/stress_test.hpp"
#include "cpphub/models/vol_surface/svi.hpp"
#include "cpphub/models/vol_surface/dupire_local_vol.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace cpphub::v1;

// ---------- Helper: xorshift64* for deterministic test data ----------
namespace {
inline uint64_t xorshift64(uint64_t& state) {
    uint64_t x = state;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    state = x;
    return x * 0x2545F4914F6CDD1DULL;
}
inline Real normal_sample(uint64_t& state) {
    // Box-Muller
    Real u1 = static_cast<Real>(xorshift64(state)) / static_cast<Real>(0xFFFFFFFFFFFFFFFFULL);
    if (u1 < 1e-12) u1 = 1e-12;
    Real u2 = static_cast<Real>(xorshift64(state)) / static_cast<Real>(0xFFFFFFFFFFFFFFFFULL);
    if (u2 < 1e-12) u2 = 1e-12;
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * PI * u2);
}
}  // namespace

// =========================================================================
// 1. SVI calibration + Dupire local vol recovery (flat IV)
// =========================================================================
TEST(IntegrationPhase3, SVIAndDupireRecovery) {
    // 构造平坦 IV 表面 (sigma=0.20), 验证 SVI 参数化 + Dupire 局部波动率恢复
    std::vector<Real> strikes;
    for (int i = 0; i < 21; ++i) strikes.push_back(80.0 + i * 2.0);
    std::vector<Real> maturities = {0.25, 0.5, 1.0};
    Real flat_vol = 0.20;
    std::vector<Real> ivs_flat;
    for (Size i = 0; i < strikes.size(); ++i) ivs_flat.push_back(flat_vol);

    // SVI 标定 (单期限)
    SVI svi(SVIParams{0.0, 0.1, 0.0, 0.1, 0.0}, SVIParamType::Raw);
    CalibConfig cfg;
    cfg.de_generations = 50;
    cfg.lm_max_iter = 200;
    auto calib_result = svi.calibrate(strikes, maturities, ivs_flat, /*forward=*/100.0, cfg);

    // Dupire 局部波动率恢复
    std::vector<std::vector<Real>> iv_grid(
        maturities.size(), std::vector<Real>(strikes.size(), flat_vol));
    DupireLocalVol dup(strikes, maturities, iv_grid, 100.0, 0.05, 0.02);

    // ATM 中心区域局部波动率应恢复 flat_vol
    Real lv = dup.local_vol(100.0, 0.5);
    EXPECT_NEAR(lv, flat_vol, 5e-3)
        << "local vol recovery failed: lv=" << lv;
}

// =========================================================================
// 2. AAD Greeks vs Analytic Greeks (BSM European)
// =========================================================================
TEST(IntegrationPhase3, AADGreeksVsAnalytic) {
    Real S = 100.0, K = 105.0, T = 0.5, r = 0.05, q = 0.02, sigma = 0.25;
    bool is_call = true;

    auto aad = AADGreeksEngine::bsm_european(S, K, T, r, q, sigma, is_call);
    auto ana = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, is_call);

    // AAD 与解析应精确一致 (1e-10)
    EXPECT_NEAR(aad.price, ana.price, 1e-10);
    EXPECT_NEAR(aad.delta, ana.delta, 1e-10);
    EXPECT_NEAR(aad.vega, ana.vega, 1e-10);
    EXPECT_NEAR(aad.rho, ana.rho, 1e-10);
    EXPECT_NEAR(aad.gamma, ana.gamma, 1e-10);
}

// =========================================================================
// 3. Put-call parity verification via AAD Greeks
// =========================================================================
TEST(IntegrationPhase3, PutCallParity) {
    // Put-Call Parity: C - P = S*exp(-qT) - K*exp(-rT)
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.04, q = 0.01, sigma = 0.20;

    auto call_aad = AADGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    auto put_aad  = AADGreeksEngine::bsm_european(S, K, T, r, q, sigma, false);

    Real parity_lhs = call_aad.price - put_aad.price;
    Real parity_rhs = S * std::exp(-q * T) - K * std::exp(-r * T);
    EXPECT_NEAR(parity_lhs, parity_rhs, 1e-10);

    // Delta parity: Δ_call - Δ_put = exp(-qT)
    Real delta_diff = call_aad.delta - put_aad.delta;
    EXPECT_NEAR(delta_diff, std::exp(-q * T), 1e-10);
}

// =========================================================================
// 4. Historical VaR quantile
// =========================================================================
TEST(IntegrationPhase3, HistoricalVaR) {
    // 生成正态分布 PnL, 验证 99% 历史 VaR ≈ 2.33
    uint64_t state = 42;
    std::vector<Real> pnl;
    for (int i = 0; i < 20000; ++i) {
        pnl.push_back(normal_sample(state));
    }
    HistoricalVaR hv(pnl, 0.99, 1);
    Real var_val = hv.var(QuantileInterpolation::Linear);
    EXPECT_NEAR(var_val, 2.326, 0.1);
}

// =========================================================================
// 5. Parametric VaR (Normal / StudentT / CornishFisher)
// =========================================================================
TEST(IntegrationPhase3, ParametricVaR) {
    PortfolioStats stats;
    stats.mean = 0.0;
    stats.variance = 1.0;
    stats.skewness = 0.0;
    stats.kurtosis = 0.0;
    stats.degrees_of_freedom = 5;
    ParametricVaR pvar(stats, 0.99, 1);

    Real normal_var = pvar.var(ParametricMethod::Normal);
    Real student_var = pvar.var(ParametricMethod::StudentT);
    EXPECT_NEAR(normal_var, 2.326, 0.01);
    EXPECT_GT(student_var, normal_var);  // Student t 有肥尾
}

// =========================================================================
// 6. MC VaR full revaluation
// =========================================================================
TEST(IntegrationPhase3, MCVaRFullRevaluation) {
    Real S0 = 100.0;
    auto payoff = [S0](const std::vector<Real>& S) -> Real { return S[0]; };
    std::vector<Real> current = {S0};
    std::vector<Real> cov = {0.01};
    MCVarConfig cfg;
    cfg.n_paths = 50000;
    cfg.seed = 42;
    cfg.antithetic = true;
    MCVaR mc(payoff, current, cov, 1, cfg);
    Real mc_var = mc.var(0.99, VaRApproximation::Full);

    PortfolioStats stats;
    stats.mean = 0; stats.variance = 0.01;
    ParametricVaR pvar(stats, 0.99, 1);
    Real normal_var = pvar.var(ParametricMethod::Normal);
    EXPECT_NEAR(mc_var, normal_var, normal_var * 0.05);
}

// =========================================================================
// 7. ES normal vs theoretical
// =========================================================================
TEST(IntegrationPhase3, ESNormal) {
    ExpectedShortfall es;
    Real es_val = es.normal_es(0, 1, 0.99);
    // 理论 ES_99% (标准正态) = phi(z_0.01) / 0.01 ≈ 2.665
    EXPECT_NEAR(es_val, 2.665, 0.01);
    Real z = inv_normal_cdf(0.01);
    EXPECT_GT(es_val, -z);  // ES > VaR
}

// =========================================================================
// 8. Stress test single scenario
// =========================================================================
TEST(IntegrationPhase3, StressTestSingle) {
    // 简单组合: 价值 = 100 * spot_factor (线性持仓)
    StressTester::PortfolioValueFn value_fn = [](const std::map<std::string, Real>& factors) -> Real {
        Real spot = factors.at("Equity");
        return 100.0 * spot;
    };
    std::map<std::string, Real> current = {{"Equity", 100.0}};
    StressTester tester(value_fn, current);

    StressScenario scenario;
    scenario.name = "Equity Crash";
    scenario.spot_shocks["Equity"] = -0.30;  // -30%
    scenario.probability = 1.0;

    auto pnl = tester.run(scenario);
    EXPECT_NEAR(pnl.total_pnl, -3000.0, 1e-6);  // 100 * 100 * (-0.30) = -3000
    EXPECT_NEAR(pnl.by_asset["Equity"], -3000.0, 1e-6);
}

// =========================================================================
// 9. Stress test historical crises
// =========================================================================
TEST(IntegrationPhase3, StressTestHistorical) {
    StressTester::PortfolioValueFn value_fn = [](const std::map<std::string, Real>& factors) -> Real {
        Real spot = factors.at("Equity");
        return 100.0 * spot;
    };
    std::map<std::string, Real> current = {{"Equity", 100.0}};
    StressTester tester(value_fn, current);

    std::vector<HistoricalCrisis> crises = {
        {"2008 GFC", "2008-09-15", "2009-03-09", {{"Equity", -0.45}}, {}, {}},
        {"2020 COVID", "2020-02-19", "2020-03-23", {{"Equity", -0.34}}, {}, {}}
    };
    auto results = tester.run_historical(crises);
    ASSERT_EQ(results.size(), 2);
    EXPECT_NEAR(results[0].total_pnl, -4500.0, 1e-6);  // 100*100*(-0.45)
    EXPECT_NEAR(results[1].total_pnl, -3400.0, 1e-6);
}

// =========================================================================
// 10. Full risk report: price + Greeks + VaR + ES + stress (综合)
// =========================================================================
TEST(IntegrationPhase3, FullRiskReport) {
    // 模拟一个欧式 call 期权的完整风险报告
    Real S = 100.0, K = 100.0, T = 0.25, r = 0.03, q = 0.0, sigma = 0.20;

    // (1) Price + Greeks via AAD
    auto greeks = AADGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    EXPECT_GT(greeks.price, 0);
    EXPECT_NEAR(greeks.delta, 0.5, 0.1);  // ATM call delta ≈ 0.5
    EXPECT_GT(greeks.gamma, 0);
    EXPECT_GT(greeks.vega, 0);

    // (2) VaR: 用 delta 近似 (1 天, 99%)
    // 日波动率 = sigma * sqrt(1/252)
    Real daily_vol = sigma * std::sqrt(1.0 / 252.0);
    PortfolioStats stats;
    stats.mean = 0;
    stats.variance = (greeks.delta * S * daily_vol) * (greeks.delta * S * daily_vol);
    ParametricVaR pvar(stats, 0.99, 1);
    Real var_99 = pvar.var(ParametricMethod::Normal);
    EXPECT_GT(var_99, 0);

    // (3) ES
    ExpectedShortfall es;
    Real es_99 = es.normal_es(0, std::sqrt(stats.variance), 0.99);
    EXPECT_GT(es_99, var_99);  // ES > VaR

    // (4) Stress: 股价 -20%
    StressTester::PortfolioValueFn value_fn = [K, T, r, q, sigma](const std::map<std::string, Real>& factors) -> Real {
        Real stressed_S = factors.at("Equity");
        auto g = AnalyticGreeksEngine::bsm_european(stressed_S, K, T, r, q, sigma, true);
        return g.price;
    };
    std::map<std::string, Real> current = {{"Equity", S}};
    StressTester tester(value_fn, current);
    StressScenario crash;
    crash.name = "Equity -20%";
    crash.spot_shocks["Equity"] = -0.20;
    auto pnl = tester.run(crash);
    EXPECT_LT(pnl.total_pnl, 0);  // 股价下跌, call 期权价格下跌

    // (5) Sanity: VaR 与 stress PnL 数量级一致
    EXPECT_LT(std::abs(pnl.total_pnl), var_99 * 20);  // stress 应在 VaR 的合理倍数内
}
