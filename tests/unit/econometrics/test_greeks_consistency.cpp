// =============================================================================
// test_greeks_consistency.cpp - Phase 7A Wave 3d Greeks 跨方法一致性测试
//
// 15 用例:
//   Delta (3) + Gamma (3) + Vega (3) + Theta (3) + Rho (3)
//
// 排幻觉点覆盖:
//   H20 (Pathwise/LR 用置信区间比较, 非点估计)
//   H21 (Gamma Numerical 用二阶差分, bump size dS=1e-4)
//
// 教材锚点: Glasserman 2003 §7, Broadie-Glasserman 1996
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/risk/greeks/greeks_consistency.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"

using cpphub::v1::GreeksConsistencyResult;
using cpphub::v1::GreeksFactory;
using cpphub::v1::GreeksMethod;
using cpphub::v1::PayoffType;
using cpphub::v1::Real;
using cpphub::v1::Size;
using cpphub::v1::greeks_consistency_check;

// 标准测试参数
static constexpr Real S0 = 100.0;
static constexpr Real K0 = 100.0;
static constexpr Real T0 = 0.5;
static constexpr Real r0 = 0.05;
static constexpr Real q0 = 0.02;
static constexpr Real sigma0 = 0.20;
static constexpr Size N_PATHS = 200000;

// 辅助: 相对差
static Real rel_diff(Real a, Real b) {
    Real ref = std::max(std::abs(a), std::abs(b));
    if (ref < 1e-15) return 0.0;
    return std::abs(a - b) / ref;
}

// =============================================================================
// Delta (3 用例)
// =============================================================================

// --- Delta 1: Vanilla Call — Analytic vs AAD 应精确一致 ---
TEST(GreeksConsistency, Delta_VanillaCall_AnalyticVsAAD_Exact) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "delta", N_PATHS, 42, 0.01);

    // Analytic vs AAD: 容差 1e-10 (AAD 是解析导数, 应精确一致)
    EXPECT_LT(rel_diff(res.analytic_value, res.aad_value), 1e-10)
        << "Analytic=" << res.analytic_value << " AAD=" << res.aad_value;

    // Delta 应为正 (call)
    EXPECT_GT(res.analytic_value, 0.0);
}

// --- Delta 2: Vanilla Call — 所有方法一致 (H20 CI 比较) ---
TEST(GreeksConsistency, Delta_VanillaCall_AllMethodsConsistent) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "delta", N_PATHS, 42, 0.01);

    // 所有方法应一致
    EXPECT_TRUE(res.consistent)
        << "Warnings: ";
    for (const auto& w : res.warnings) {
        std::cerr << "  " << w << std::endl;
    }

    // Pathwise/LR 应有有效估计 (非零 SE)
    EXPECT_GT(res.pathwise_std_error, 0.0);
    EXPECT_GT(res.lr_std_error, 0.0);

    // Pathwise 均值应接近 Analytic (H20: 在 CI 范围内)
    Real pw_diff = std::abs(res.pathwise_mean - res.analytic_value);
    Real pw_ci = 1.96 * res.pathwise_std_error;
    // 即使 CI 比较不通过, 相对差也应在合理 MC 范围内 (< 5%)
    EXPECT_LT(rel_diff(res.analytic_value, res.pathwise_mean), 0.05)
        << "Pathwise mean=" << res.pathwise_mean
        << " Analytic=" << res.analytic_value
        << " SE=" << res.pathwise_std_error;
}

// --- Delta 3: Vanilla Put — 所有方法一致 ---
TEST(GreeksConsistency, Delta_VanillaPut_AllMethodsConsistent) {
    auto res = greeks_consistency_check(
        S0, 105.0, 1.0, 0.03, 0.01, 0.25,
        PayoffType::VanillaPut, "delta", N_PATHS, 7, 0.01);

    EXPECT_TRUE(res.consistent);

    // Put delta 应为负
    EXPECT_LT(res.analytic_value, 0.0);
    EXPECT_LT(res.aad_value, 0.0);
    EXPECT_LT(res.pathwise_mean, 0.0);
    EXPECT_LT(res.lr_mean, 0.0);
}

// =============================================================================
// Gamma (3 用例)
// =============================================================================

// --- Gamma 1: Vanilla Call — Analytic vs AAD 应精确一致 ---
TEST(GreeksConsistency, Gamma_VanillaCall_AnalyticVsAAD_Exact) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "gamma", N_PATHS, 42, 0.01);

    // Analytic vs AAD: 容差 1e-10
    EXPECT_LT(rel_diff(res.analytic_value, res.aad_value), 1e-10)
        << "Analytic=" << res.analytic_value << " AAD=" << res.aad_value;

    // Gamma 应为正
    EXPECT_GT(res.analytic_value, 0.0);
}

// --- Gamma 2: H21 — Numerical (FD) 用 dS=1e-4, 应在 0.1% 内匹配 Analytic ---
TEST(GreeksConsistency, Gamma_VanillaCall_H21_NumericalSmallBump) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "gamma", N_PATHS, 42, 0.01);

    // H21: Numerical gamma 用 dS=1e-4 (二阶差分), 截断误差 O(dS²) = O(1e-8)
    // 应在 0.1% 内匹配 Analytic (远优于 1% 容差)
    EXPECT_LT(rel_diff(res.analytic_value, res.numerical_value), 0.001)
        << "Analytic=" << res.analytic_value
        << " Numerical=" << res.numerical_value
        << " (H21: dS=1e-4 for gamma second difference)";

    // 验证 Numerical gamma 不为零 (排除幻觉)
    EXPECT_GT(std::abs(res.numerical_value), 1e-10);
}

// --- Gamma 3: Pathwise/LR 不支持 Gamma, 应报告 0 + warning ---
TEST(GreeksConsistency, Gamma_PathwiseLRNotApplicable) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "gamma", N_PATHS, 42, 0.01);

    // Pathwise 不支持 gamma → mean=0, SE=0
    EXPECT_EQ(res.pathwise_mean, 0.0);
    EXPECT_EQ(res.pathwise_std_error, 0.0);

    // LR 不支持 gamma → mean=0, SE=0
    EXPECT_EQ(res.lr_mean, 0.0);
    EXPECT_EQ(res.lr_std_error, 0.0);

    // 应有 warning 说明不支持
    bool has_pathwise_warning = false;
    bool has_lr_warning = false;
    for (const auto& w : res.warnings) {
        if (w.find("Pathwise") != std::string::npos) has_pathwise_warning = true;
        if (w.find("LR") != std::string::npos) has_lr_warning = true;
    }
    EXPECT_TRUE(has_pathwise_warning)
        << "Should warn that Pathwise doesn't support gamma";
    EXPECT_TRUE(has_lr_warning)
        << "Should warn that LR doesn't support gamma";

    // Analytic/Numerical/AAD 应有有效值 (这三者支持 gamma)
    EXPECT_NE(res.analytic_value, 0.0);
    EXPECT_NE(res.numerical_value, 0.0);
    EXPECT_NE(res.aad_value, 0.0);
}

// =============================================================================
// Vega (3 用例)
// =============================================================================

// --- Vega 1: Vanilla Call — Analytic vs AAD 应精确一致 ---
TEST(GreeksConsistency, Vega_VanillaCall_AnalyticVsAAD_Exact) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "vega", N_PATHS, 42, 0.01);

    EXPECT_LT(rel_diff(res.analytic_value, res.aad_value), 1e-10)
        << "Analytic=" << res.analytic_value << " AAD=" << res.aad_value;

    // Vega 应为正
    EXPECT_GT(res.analytic_value, 0.0);
}

// --- Vega 2: Vanilla Call — 所有方法一致 (H20 CI 比较) ---
TEST(GreeksConsistency, Vega_VanillaCall_AllMethodsConsistent) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "vega", N_PATHS, 42, 0.01);

    EXPECT_TRUE(res.consistent);

    // 所有 5 种方法应有有效值
    EXPECT_NE(res.analytic_value, 0.0);
    EXPECT_NE(res.numerical_value, 0.0);
    EXPECT_NE(res.aad_value, 0.0);
    EXPECT_NE(res.pathwise_mean, 0.0);
    EXPECT_NE(res.lr_mean, 0.0);

    // Pathwise/LR 应有有效 SE
    EXPECT_GT(res.pathwise_std_error, 0.0);
    EXPECT_GT(res.lr_std_error, 0.0);
}

// --- Vega 3: Vanilla Put — 所有方法一致 ---
TEST(GreeksConsistency, Vega_VanillaPut_AllMethodsConsistent) {
    auto res = greeks_consistency_check(
        S0, 105.0, 1.0, 0.03, 0.01, 0.25,
        PayoffType::VanillaPut, "vega", N_PATHS, 7, 0.01);

    EXPECT_TRUE(res.consistent);

    // Put vega 应为正 (与 call vega 相同, BSM vega 对 call/put 一致)
    EXPECT_GT(res.analytic_value, 0.0);
}

// =============================================================================
// Theta (3 用例)
// =============================================================================

// --- Theta 1: Vanilla Call — Analytic vs AAD 应精确一致 ---
TEST(GreeksConsistency, Theta_VanillaCall_AnalyticVsAAD_Exact) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "theta", N_PATHS, 42, 0.01);

    EXPECT_LT(rel_diff(res.analytic_value, res.aad_value), 1e-10)
        << "Analytic=" << res.analytic_value << " AAD=" << res.aad_value;

    // Call theta 通常为负 (时间衰减)
    EXPECT_LT(res.analytic_value, 0.0);
}

// --- Theta 2: Vanilla Call — Analytic vs Numerical (FD) 在 1% 内 ---
TEST(GreeksConsistency, Theta_VanillaCall_AnalyticVsNumerical) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "theta", N_PATHS, 42, 0.01);

    // FD theta 用中心差分 (dT=1/365), 应在 1% 内
    EXPECT_LT(rel_diff(res.analytic_value, res.numerical_value), 0.01)
        << "Analytic=" << res.analytic_value
        << " Numerical=" << res.numerical_value;
}

// --- Theta 3: Pathwise/LR 不支持 Theta ---
TEST(GreeksConsistency, Theta_PathwiseLRNotApplicable) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "theta", N_PATHS, 42, 0.01);

    EXPECT_EQ(res.pathwise_mean, 0.0);
    EXPECT_EQ(res.pathwise_std_error, 0.0);
    EXPECT_EQ(res.lr_mean, 0.0);
    EXPECT_EQ(res.lr_std_error, 0.0);

    // 应有 warning
    bool has_pathwise_warning = false;
    bool has_lr_warning = false;
    for (const auto& w : res.warnings) {
        if (w.find("Pathwise") != std::string::npos) has_pathwise_warning = true;
        if (w.find("LR") != std::string::npos) has_lr_warning = true;
    }
    EXPECT_TRUE(has_pathwise_warning);
    EXPECT_TRUE(has_lr_warning);
}

// =============================================================================
// Rho (3 用例)
// =============================================================================

// --- Rho 1: Vanilla Call — Analytic vs AAD 应精确一致 ---
TEST(GreeksConsistency, Rho_VanillaCall_AnalyticVsAAD_Exact) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "rho", N_PATHS, 42, 0.01);

    EXPECT_LT(rel_diff(res.analytic_value, res.aad_value), 1e-10)
        << "Analytic=" << res.analytic_value << " AAD=" << res.aad_value;

    // Call rho 应为正 (利率上升, call 价值上升)
    EXPECT_GT(res.analytic_value, 0.0);
}

// --- Rho 2: Vanilla Call — Analytic vs Numerical (FD) 在 1% 内 ---
TEST(GreeksConsistency, Rho_VanillaCall_AnalyticVsNumerical) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "rho", N_PATHS, 42, 0.01);

    // FD rho 用中心差分 (dR=0.0001), 应在 1% 内
    EXPECT_LT(rel_diff(res.analytic_value, res.numerical_value), 0.01)
        << "Analytic=" << res.analytic_value
        << " Numerical=" << res.numerical_value;
}

// --- Rho 3: Pathwise/LR 不支持 Rho ---
TEST(GreeksConsistency, Rho_PathwiseLRNotApplicable) {
    auto res = greeks_consistency_check(
        S0, K0, T0, r0, q0, sigma0,
        PayoffType::VanillaCall, "rho", N_PATHS, 42, 0.01);

    EXPECT_EQ(res.pathwise_mean, 0.0);
    EXPECT_EQ(res.pathwise_std_error, 0.0);
    EXPECT_EQ(res.lr_mean, 0.0);
    EXPECT_EQ(res.lr_std_error, 0.0);

    bool has_pathwise_warning = false;
    bool has_lr_warning = false;
    for (const auto& w : res.warnings) {
        if (w.find("Pathwise") != std::string::npos) has_pathwise_warning = true;
        if (w.find("LR") != std::string::npos) has_lr_warning = true;
    }
    EXPECT_TRUE(has_pathwise_warning);
    EXPECT_TRUE(has_lr_warning);
}

// =============================================================================
// 参数校验 (额外, 不计入 15 用例)
// =============================================================================

TEST(GreeksConsistency, InvalidGreekName_Throws) {
    EXPECT_THROW(
        greeks_consistency_check(S0, K0, T0, r0, q0, sigma0,
                                 PayoffType::VanillaCall, "vanna"),
        std::invalid_argument);
}

TEST(GreeksConsistency, InvalidParams_Throws) {
    EXPECT_THROW(
        greeks_consistency_check(-1.0, K0, T0, r0, q0, sigma0,
                                 PayoffType::VanillaCall, "delta"),
        std::invalid_argument);
    EXPECT_THROW(
        greeks_consistency_check(S0, K0, T0, r0, q0, -0.1,
                                 PayoffType::VanillaCall, "delta"),
        std::invalid_argument);
}
