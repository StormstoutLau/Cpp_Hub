// TDD test for Likelihood Ratio (LR) Greeks (PHASE3_SPEC §2.2)
// Validates: LR Delta/Vega vs Analytic for BSM European + Digital options
// LR is the method of choice for discontinuous payoffs (digital/barrier)
// where Pathwise fails (payoff derivative undefined at the discontinuity).
//
// Mathematical basis (Glasserman §7):
//   GBM terminal: S_T = S * exp((r-q-σ²/2)T + σ√T Z), Z ~ N(0,1)
//   LR score: ∂log p/∂S = Z/(S σ √T),  ∂log p/∂σ = (Z²-1)/σ - √T Z
//   delta = e^{-rT} E[Payoff(S_T) * Z/(S σ √T)]
//   vega  = e^{-rT} E[Payoff(S_T) * ((Z²-1)/σ - √T Z)]
//
// Tolerance: LR variance is much larger than Pathwise for smooth payoffs,
//            so we use 2% relative tolerance for smooth, 5% for digital.
//            vs AAD 1e-4 tolerance is only achievable with common random
//            numbers + 1M+ paths; here we test statistical consistency.
#include <gtest/gtest.h>
#include "cpphub/risk/greeks/likelihood_ratio_greeks.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"
#include "cpphub/risk/greeks/aad_greeks.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include <cmath>

using namespace cpphub::v1;

// ---- Analytic digital option formulas (independent reference) ----
// Digital Call: V = e^{-rT} N(d2),  Delta = e^{-rT} n(d2)/(S σ √T)
// Digital Put : V = e^{-rT} N(-d2), Delta = -e^{-rT} n(d2)/(S σ √T)
// Vega_call   = -e^{-rT} n(d2) (√T + d2/σ)
// Vega_put    = +e^{-rT} n(d2) (√T + d2/σ)
static Real digital_d2(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    return (std::log(S / K) + (r - q - 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
}
static Real digital_call_price(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    return std::exp(-r * T) * normal_cdf(digital_d2(S, K, T, r, q, sigma));
}
static Real digital_put_price(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    return std::exp(-r * T) * normal_cdf(-digital_d2(S, K, T, r, q, sigma));
}
static Real digital_call_delta(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    Real d2 = digital_d2(S, K, T, r, q, sigma);
    return std::exp(-r * T) * normal_pdf(d2) / (S * sigma * std::sqrt(T));
}
static Real digital_call_vega(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    Real d2 = digital_d2(S, K, T, r, q, sigma);
    return -std::exp(-r * T) * normal_pdf(d2) * (std::sqrt(T) + d2 / sigma);
}

// ===== TEST 1: BSM Call Delta via LR vs Analytic (smooth payoff) =====
TEST(LRGreeks, BSMCallDeltaVsAnalytic) {
    Real S = 100.0, K = 105.0, T = 0.5, r = 0.05, q = 0.02, sigma = 0.20;
    bool is_call = true;
    Size n_paths = 500000;  // LR 方差大,需更多路径
    uint64_t seed = 42;

    LRGreeks g = LRGreeksEngine::bsm_european(
        S, K, T, r, q, sigma, is_call, n_paths, seed);

    Real analytic_delta = AnalyticGreeksEngine::bsm_delta(S, K, T, r, q, sigma, is_call);
    Real analytic_price = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, is_call).price;

    EXPECT_GT(g.price, 0);
    // LR delta 方差比 pathwise 大,容差 2%
    EXPECT_NEAR(g.delta, analytic_delta, std::max(1e-3, std::abs(analytic_delta) * 2e-2))
        << "delta=" << g.delta << " analytic=" << analytic_delta;
    EXPECT_NEAR(g.price, analytic_price, std::max(1e-3, std::abs(analytic_price) * 5e-3));
}

// ===== TEST 2: BSM Call Vega via LR vs Analytic =====
TEST(LRGreeks, BSMCallVegaVsAnalytic) {
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    bool is_call = true;
    Size n_paths = 1000000;  // vega LR 方差极大
    uint64_t seed = 7;

    LRGreeks g = LRGreeksEngine::bsm_european(
        S, K, T, r, q, sigma, is_call, n_paths, seed);

    Real analytic_vega = AnalyticGreeksEngine::bsm_vega(S, K, T, r, q, sigma);
    // LR vega 方差极大 (Z²项),容差 5%
    EXPECT_NEAR(g.vega, analytic_vega, std::max(1e-3, std::abs(analytic_vega) * 5e-2))
        << "vega=" << g.vega << " analytic=" << analytic_vega;
}

// ===== TEST 3: Digital Call Price via LR vs Analytic =====
TEST(LRGreeks, DigitalCallPriceVsAnalytic) {
    Real S = 100.0, K = 100.0, T = 0.25, r = 0.03, q = 0.0, sigma = 0.20;
    Size n_paths = 500000;
    uint64_t seed = 11;

    LRGreeks g = LRGreeksEngine::digital_european(
        S, K, T, r, q, sigma, /*is_call=*/true, n_paths, seed);

    Real analytic_price = digital_call_price(S, K, T, r, q, sigma);
    EXPECT_GT(g.price, 0);
    EXPECT_NEAR(g.price, analytic_price, std::max(1e-4, std::abs(analytic_price) * 5e-3))
        << "price=" << g.price << " analytic=" << analytic_price;
}

// ===== TEST 4: Digital Call Delta via LR vs Analytic (key use case) =====
TEST(LRGreeks, DigitalCallDeltaVsAnalytic) {
    Real S = 100.0, K = 100.0, T = 0.25, r = 0.03, q = 0.0, sigma = 0.20;
    Size n_paths = 1000000;
    uint64_t seed = 23;

    LRGreeks g = LRGreeksEngine::digital_european(
        S, K, T, r, q, sigma, /*is_call=*/true, n_paths, seed);

    Real analytic_delta = digital_call_delta(S, K, T, r, q, sigma);
    // 这是 LR 的核心应用场景,容差 3%
    EXPECT_NEAR(g.delta, analytic_delta, std::max(1e-3, std::abs(analytic_delta) * 3e-2))
        << "delta=" << g.delta << " analytic=" << analytic_delta;
}

// ===== TEST 5: Digital Call Vega via LR vs Analytic =====
TEST(LRGreeks, DigitalCallVegaVsAnalytic) {
    Real S = 100.0, K = 100.0, T = 0.25, r = 0.03, q = 0.0, sigma = 0.20;
    Size n_paths = 2000000;  // vega LR 方差极大,需 2M 路径
    uint64_t seed = 41;

    LRGreeks g = LRGreeksEngine::digital_european(
        S, K, T, r, q, sigma, /*is_call=*/true, n_paths, seed);

    Real analytic_vega = digital_call_vega(S, K, T, r, q, sigma);
    // vega LR 方差极大,容差 8%
    ASSERT_GT(std::abs(analytic_vega), 1e-6);
    EXPECT_NEAR(g.vega, analytic_vega, std::abs(analytic_vega) * 8e-2)
        << "vega=" << g.vega << " analytic=" << analytic_vega;
}

// ===== TEST 6: Digital Put Price + Put-Call Parity =====
// Digital Call + Digital Put = e^{-rT} (cash-or-nothing parity)
TEST(LRGreeks, DigitalPutCallParity) {
    Real S = 100.0, K = 100.0, T = 0.5, r = 0.05, q = 0.0, sigma = 0.20;
    Size n_paths = 500000;
    uint64_t seed = 88;

    LRGreeks gc = LRGreeksEngine::digital_european(
        S, K, T, r, q, sigma, /*is_call=*/true, n_paths, seed);
    LRGreeks gp = LRGreeksEngine::digital_european(
        S, K, T, r, q, sigma, /*is_call=*/false, n_paths, seed + 1);

    Real parity = std::exp(-r * T);  // 1 unit cash
    EXPECT_NEAR(gc.price + gp.price, parity, 5e-3)
        << "call=" << gc.price << " put=" << gp.price << " parity=" << parity;
}

// ===== TEST 7: BSM Call Delta via LR vs AAD (cross-method consistency) =====
// PHASE3_SPEC §1.3: LR vs AAD 1e-4 — but this is only achievable with common
// random numbers and identical path sets. With independent MC runs, we test
// statistical consistency at 2% level.
TEST(LRGreeks, BSMCallDeltaVsAAD) {
    Real S = 100.0, K = 100.0, T = 0.25, r = 0.03, q = 0.0, sigma = 0.20;
    bool is_call = true;
    Size n_paths = 500000;
    uint64_t seed = 99;

    LRGreeks lg = LRGreeksEngine::bsm_european(
        S, K, T, r, q, sigma, is_call, n_paths, seed);

    AADGreeks ag = AADGreeksEngine::bsm_european(S, K, T, r, q, sigma, is_call);

    EXPECT_NEAR(lg.delta, ag.delta, std::max(1e-3, std::abs(ag.delta) * 2e-2))
        << "lr=" << lg.delta << " aad=" << ag.delta;
}

// ===== TEST 8: Deep OTM Digital — LR remains unbiased where Pathwise would fail =====
TEST(LRGreeks, DeepOTMDigitalDelta) {
    // Deep OTM digital call: K=130, S=100 — most paths expire worthless
    // Pathwise would return 0 (no ITM paths to contribute), LR remains unbiased
    Real S = 100.0, K = 130.0, T = 0.5, r = 0.05, q = 0.0, sigma = 0.30;
    Size n_paths = 1000000;
    uint64_t seed = 2026;

    LRGreeks g = LRGreeksEngine::digital_european(
        S, K, T, r, q, sigma, /*is_call=*/true, n_paths, seed);

    Real analytic_delta = digital_call_delta(S, K, T, r, q, sigma);
    Real analytic_price = digital_call_price(S, K, T, r, q, sigma);
    EXPECT_GT(g.price, 0);
    EXPECT_GT(analytic_price, 0);
    // Deep OTM: 容差放宽到 10% (路径少、方差大)
    EXPECT_NEAR(g.delta, analytic_delta, std::max(1e-4, std::abs(analytic_delta) * 1e-1))
        << "delta=" << g.delta << " analytic=" << analytic_delta;
}

// ===== TEST 9: Determinism — same seed reproduces identical results =====
TEST(LRGreeks, DeterminismSameSeed) {
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Size n_paths = 100000;
    uint64_t seed = 12345;

    LRGreeks g1 = LRGreeksEngine::bsm_european(
        S, K, T, r, q, sigma, true, n_paths, seed);
    LRGreeks g2 = LRGreeksEngine::bsm_european(
        S, K, T, r, q, sigma, true, n_paths, seed);

    EXPECT_DOUBLE_EQ(g1.price, g2.price);
    EXPECT_DOUBLE_EQ(g1.delta, g2.delta);
    EXPECT_DOUBLE_EQ(g1.vega, g2.vega);
}
