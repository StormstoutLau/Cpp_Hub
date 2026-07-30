// TDD test for GreeksFactory (PHASE3_SPEC §2.2)
// Validates: Auto-dispatch routes vanilla→Analytic, digital→LR, explicit methods work
// Also validates: cross-method consistency (Analytic vs AAD vs Pathwise for vanilla)
#include <gtest/gtest.h>
#include "cpphub/risk/greeks/greeks_factory.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"
#include <cmath>

using namespace cpphub::v1;

// ===== TEST 1: Auto routes vanilla call to Analytic =====
TEST(GreeksFactory, AutoVanillaCallUsesAnalytic) {
    Real S = 100.0, K = 100.0, T = 0.5, r = 0.05, q = 0.02, sigma = 0.20;
    auto g = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::VanillaCall, GreeksMethod::Auto);
    EXPECT_EQ(g.method_used, GreeksMethod::Analytic);

    auto ref = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    EXPECT_NEAR(g.price, ref.price, 1e-12);
    EXPECT_NEAR(g.delta, ref.delta, 1e-12);
    EXPECT_NEAR(g.vega,  ref.vega,  1e-12);
    EXPECT_NEAR(g.gamma, ref.gamma, 1e-12);
}

// ===== TEST 2: Auto routes vanilla put to Analytic =====
TEST(GreeksFactory, AutoVanillaPutUsesAnalytic) {
    Real S = 100.0, K = 105.0, T = 1.0, r = 0.03, q = 0.01, sigma = 0.25;
    auto g = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::VanillaPut, GreeksMethod::Auto);
    EXPECT_EQ(g.method_used, GreeksMethod::Analytic);
    EXPECT_LT(g.delta, 0.0);  // put delta < 0
}

// ===== TEST 3: Auto routes digital call to LR (discontinuous payoff) =====
TEST(GreeksFactory, AutoDigitalCallUsesLR) {
    Real S = 100.0, K = 100.0, T = 0.25, r = 0.03, q = 0.0, sigma = 0.20;
    auto g = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::DigitalCall,
        GreeksMethod::Auto, /*n_paths=*/500000, /*seed=*/42);
    EXPECT_EQ(g.method_used, GreeksMethod::LR);

    // LR price 应接近解析解 e^{-rT} N(d2)
    Real d2 = (std::log(S / K) + (r - q - 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    Real analytic_price = std::exp(-r * T) * normal_cdf(d2);
    EXPECT_NEAR(g.price, analytic_price, std::max(1e-4, std::abs(analytic_price) * 5e-3));
}

// ===== TEST 4: Auto routes digital put to LR =====
TEST(GreeksFactory, AutoDigitalPutUsesLR) {
    Real S = 100.0, K = 100.0, T = 0.5, r = 0.05, q = 0.0, sigma = 0.20;
    auto g = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::DigitalPut,
        GreeksMethod::Auto, /*n_paths=*/500000, /*seed=*/7);
    EXPECT_EQ(g.method_used, GreeksMethod::LR);
    EXPECT_GT(g.price, 0.0);
}

// ===== TEST 5: Explicit Pathwise on vanilla call =====
TEST(GreeksFactory, ExplicitPathwiseVanillaCall) {
    Real S = 100.0, K = 100.0, T = 0.5, r = 0.05, q = 0.0, sigma = 0.20;
    auto g = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::VanillaCall,
        GreeksMethod::Pathwise, /*n_paths=*/200000, /*seed=*/11);
    EXPECT_EQ(g.method_used, GreeksMethod::Pathwise);

    auto ref = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    EXPECT_NEAR(g.delta, ref.delta, std::max(1e-3, std::abs(ref.delta) * 1e-2));
}

// ===== TEST 6: Explicit Pathwise on digital falls back to LR =====
TEST(GreeksFactory, ExplicitPathwiseDigitalFallsBackToLR) {
    Real S = 100.0, K = 100.0, T = 0.25, r = 0.03, q = 0.0, sigma = 0.20;
    auto g = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::DigitalCall,
        GreeksMethod::Pathwise, /*n_paths=*/300000, /*seed=*/23);
    EXPECT_EQ(g.method_used, GreeksMethod::LR);  // Fallback to LR
    EXPECT_FALSE(g.note.empty());  // Should have a note explaining the fallback
}

// ===== TEST 7: Explicit AAD on vanilla =====
TEST(GreeksFactory, ExplicitAADVanillaCall) {
    Real S = 100.0, K = 105.0, T = 0.5, r = 0.05, q = 0.02, sigma = 0.20;
    auto g = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::VanillaCall, GreeksMethod::AAD);
    EXPECT_EQ(g.method_used, GreeksMethod::AAD);

    auto ref = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    EXPECT_NEAR(g.delta, ref.delta, 1e-10);
    EXPECT_NEAR(g.vega,  ref.vega,  1e-10);
    EXPECT_NEAR(g.gamma, ref.gamma, 1e-10);
}

// ===== TEST 8: Explicit FD on vanilla =====
TEST(GreeksFactory, ExplicitFDVanillaCall) {
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    auto g = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::VanillaCall, GreeksMethod::FD);
    EXPECT_EQ(g.method_used, GreeksMethod::FD);

    auto ref = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    // FD 中心差分精度 ~1e-6
    EXPECT_NEAR(g.delta, ref.delta, 1e-4);
    EXPECT_NEAR(g.gamma, ref.gamma, 1e-4);
    EXPECT_NEAR(g.vega,  ref.vega,  1e-3);
}

// ===== TEST 9: Explicit FD on digital =====
TEST(GreeksFactory, ExplicitFDDigitalCall) {
    Real S = 100.0, K = 100.0, T = 0.25, r = 0.03, q = 0.0, sigma = 0.20;
    auto g = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::DigitalCall, GreeksMethod::FD);
    EXPECT_EQ(g.method_used, GreeksMethod::FD);

    Real d2 = (std::log(S / K) + (r - q - 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    Real analytic_price = std::exp(-r * T) * normal_cdf(d2);
    EXPECT_NEAR(g.price, analytic_price, 1e-6);
    // FD delta for digital will be noisy near strike; just check non-zero
    EXPECT_NE(g.delta, 0.0);
}

// ===== TEST 10: Cross-method consistency for vanilla call =====
// Analytic, AAD, Pathwise, FD should all agree on delta within MC noise
TEST(GreeksFactory, CrossMethodDeltaConsistency) {
    Real S = 100.0, K = 100.0, T = 0.5, r = 0.05, q = 0.0, sigma = 0.20;

    auto g_ana = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::VanillaCall, GreeksMethod::Analytic);
    auto g_aad = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::VanillaCall, GreeksMethod::AAD);
    auto g_pw  = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::VanillaCall,
        GreeksMethod::Pathwise, /*n_paths=*/500000, /*seed=*/1);
    auto g_fd  = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::VanillaCall, GreeksMethod::FD);

    // AAD vs Analytic: machine precision
    EXPECT_NEAR(g_aad.delta, g_ana.delta, 1e-10);
    // FD vs Analytic: 1e-4 (centered difference)
    EXPECT_NEAR(g_fd.delta,  g_ana.delta, 1e-4);
    // Pathwise vs Analytic: 1% (MC noise)
    EXPECT_NEAR(g_pw.delta,  g_ana.delta, std::max(1e-3, std::abs(g_ana.delta) * 1e-2));
}

// ===== TEST 11: Put-call parity via Factory =====
TEST(GreeksFactory, PutCallParity) {
    Real S = 100.0, K = 100.0, T = 0.5, r = 0.05, q = 0.0, sigma = 0.20;
    auto gc = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::VanillaCall, GreeksMethod::Analytic);
    auto gp = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, PayoffType::VanillaPut, GreeksMethod::Analytic);
    // C - P = S e^{-qT} - K e^{-rT}
    Real parity = S * std::exp(-q * T) - K * std::exp(-r * T);
    EXPECT_NEAR(gc.price - gp.price, parity, 1e-10);
}
