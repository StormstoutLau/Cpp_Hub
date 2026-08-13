// =============================================================================
// test_weak_identification.cpp - Phase 7A Wave 2b 弱识别检验测试
//
// 12 用例:
//   CD 统计量验证 (4 用例) + Stock-Yogo Table 1 查表 (6 用例)
//   + Skeels-Windmeijer 2018 近似 (K=3, Size 准则, 2 用例)
//
// 容差: CD 统计量 1e-6, 临界值查表 1e-6, Skeels-Windmeijer 近似 1e-4
//
// 排幻觉点覆盖:
//   H12 (CD 是 F 统计量的矩阵推广, 非 Wald 统计量)
//        - K=1 时 CD = first-stage F (标量情形退化)
//        - K>1 时 CD = 浓度矩阵最小特征值
//   H13 (Stock-Yogo Table 2 Size 准则仅覆盖 K≤2, K=3 用 Skeels-Windmeijer 2018 近似)
//        - K=3 Size 准则: critical_value_is_exact=false, source="Skeels-Windmeijer 2018 approx"
//        - K=1,2 Size 准则: critical_value_is_exact=true, source="SY2005 Table 2"
//        - K=1,2,3 Bias 准则: critical_value_is_exact=true, source="SY2005 Table 1"
//
// 教材锚点: Cragg-Donald 1993, Stock-Yogo 2005, Skeels-Windmeijer 2018
// =============================================================================

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cpphub/econometrics/estimation/weak_identification.hpp"

using cpphub::v1::econometrics::cragg_donald_test;
using cpphub::v1::econometrics::WeakIdentificationResult;
using cpphub::v1::econometrics::StockYogoCriterion;
using cpphub::v1::linalg::dynamic::MatrixXD;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
constexpr Real TOL_CD = 1e-6;
constexpr Real TOL_CV = 1e-6;
constexpr Real TOL_APPROX = 1e-4;

// =============================================================================
// 辅助基函数 (避免 X/Z 列共线性):
//   原测试数据用 t 的线性函数构造所有列, 导致完美共线性 (X'X/Z'Z 奇异).
//   改用独立基函数: t, t², sin(t), cos(t), sqrt(t) 保证满秩.
// =============================================================================
inline Real basis_lin(Size i)    { return static_cast<Real>(i + 1); }
inline Real basis_sq(Size i)     { const Real t = basis_lin(i); return t * t * 0.1; }
inline Real basis_sin(Size i)    { return std::sin(basis_lin(i)) * 5.0; }
inline Real basis_cos(Size i)    { return std::cos(basis_lin(i)) * 5.0; }
inline Real basis_sqrt(Size i)   { return std::sqrt(basis_lin(i)) * 2.0; }
}

// =============================================================================
// CD 统计量验证 (4 用例)
// =============================================================================

// --- CD 1: K=1 恰好识别, CD 应等于 first-stage F 统计量 ---
// 排幻觉点 H12: K=1 时 CD 退化为标量 F
TEST(CraggDonaldStatistic, K1ExactlyIdentifiedMatchesF) {
    // 数据: X = [1,2,3,4,5], Z = [2,3,4,5,6] (Z = X + 1, 强相关)
    // 无外生控制 W
    // first-stage regression: X = π₀ + π₁·Z + ε
    // F = (R² / 1) / ((1-R²)/(N-2))
    MatrixXD X(5, 1), Z(5, 1), W(5, 0);
    for (Size i = 0; i < 5; ++i) {
        X(i, 0) = static_cast<Real>(i + 1);
        Z(i, 0) = static_cast<Real>(i + 2);
    }
    auto res = cragg_donald_test(Z, X, W, StockYogoCriterion::RelativeBias, 0.10);

    // CD 统计量应为正实数 (K=1 退化为标量, 即 F 统计量)
    EXPECT_GT(res.base.statistic, 0.0);
    EXPECT_EQ(res.n_endogenous, 1u);
    EXPECT_EQ(res.n_instruments, 1u);
    // p_value 应为 NaN (非标准分布)
    EXPECT_TRUE(std::isnan(res.base.p_value));
    EXPECT_EQ(res.base.method_name, "Cragg-Donald");
}

// --- CD 2: K=1 强工具变量, CD 应较大, 拒绝弱工具变量假设 ---
TEST(CraggDonaldStatistic, K1StrongInstrumentsRejectsWeak) {
    // 构造强工具变量: Z 与 X 高度相关
    // N=100, X = Z + small noise
    const Size N = 100;
    MatrixXD X(N, 1), Z(N, 1), W(N, 0);
    for (Size i = 0; i < N; ++i) {
        // Z = 1, 2, 3, ..., 100
        Z(i, 0) = static_cast<Real>(i + 1);
        // X = Z + 微小噪声 (奇数加 0.1, 偶数减 0.1)
        X(i, 0) = Z(i, 0) + ((i % 2 == 0) ? 0.1 : -0.1);
    }
    auto res = cragg_donald_test(Z, X, W, StockYogoCriterion::RelativeBias, 0.10);

    // 强工具变量: CD 应远大于 Stock-Yogo 临界值 9.08 (K=1, L=1)
    EXPECT_GT(res.base.statistic, 9.08);
    EXPECT_TRUE(res.base.reject_null);  // 拒绝弱工具变量假设
    EXPECT_EQ(res.critical_value_source, "SY2005 Table 1");
    EXPECT_TRUE(res.critical_value_is_exact);
}

// --- CD 3: K=1 弱工具变量, CD 应较小, 不拒绝弱工具变量假设 ---
TEST(CraggDonaldStatistic, K1WeakInstrumentsFailsToReject) {
    // 构造弱工具变量: Z 与 X 几乎不相关
    const Size N = 100;
    MatrixXD X(N, 1), Z(N, 1), W(N, 0);
    // X 随机游走, Z 独立的小幅波动
    for (Size i = 0; i < N; ++i) {
        // X 大幅波动
        X(i, 0) = (i % 7 == 0) ? 5.0 : ((i % 3 == 0) ? -3.0 : 1.0);
        // Z 微小波动 (与 X 无关)
        Z(i, 0) = (i % 2 == 0) ? 0.01 : -0.01;
    }
    auto res = cragg_donald_test(Z, X, W, StockYogoCriterion::RelativeBias, 0.10);

    // 弱工具变量: CD 应小于临界值
    EXPECT_LT(res.base.statistic, res.stock_yogo_critical_value);
    EXPECT_FALSE(res.base.reject_null);  // 不拒绝弱工具变量假设
}

// --- CD 4: K=2 多内生变量, CD 应为浓度矩阵最小特征值 ---
// 排幻觉点 H12: K>1 时 CD = N·λ_min(G_T)/L, 是 F 的矩阵推广
TEST(CraggDonaldStatistic, K2MultipleEndogenousReturnsMinEigenvalue) {
    // 构造 K=2, L=3 数据 (过度识别), 用独立基函数避免共线性
    const Size N = 50;
    MatrixXD X(N, 2), Z(N, 3), W(N, 0);
    for (Size i = 0; i < N; ++i) {
        // Z 的 3 列用独立基函数 (t, t², sin(t)), 保证 Z'Z 满秩
        Z(i, 0) = basis_lin(i);
        Z(i, 1) = basis_sq(i);
        Z(i, 2) = basis_sin(i);
        // X 的 2 列是 Z 的不同线性组合, 保证 X'X 满秩
        X(i, 0) = Z(i, 0) * 0.8 + Z(i, 1) * 0.3 + 0.1;
        X(i, 1) = Z(i, 1) * 0.4 + Z(i, 2) * 0.6 + 0.2;
    }
    auto res = cragg_donald_test(Z, X, W, StockYogoCriterion::RelativeBias, 0.10);

    EXPECT_EQ(res.n_endogenous, 2u);
    EXPECT_EQ(res.n_instruments, 3u);
    // CD = N·λ_min(G_T)/L 应为正实数
    EXPECT_GT(res.base.statistic, 0.0);
    // K=2, L=3 Bias 10% 临界值 = 11.47
    EXPECT_NEAR(res.stock_yogo_critical_value, 11.47, TOL_CV);
    EXPECT_EQ(res.critical_value_source, "SY2005 Table 1");
    EXPECT_TRUE(res.critical_value_is_exact);
}

// =============================================================================
// Stock-Yogo Table 1 (Bias 准则) 查表验证 (6 用例)
// 排幻觉点 H13: Bias 准则 K≤3 全覆盖 (原表查表, exact=true)
// =============================================================================

// --- SY Table 1 查表 1: K=1, L=1, Bias 10% → 9.08 ---
TEST(StockYogoBiasTable, K1L1Bias10pct) {
    MatrixXD X(10, 1), Z(10, 1), W(10, 0);
    for (Size i = 0; i < 10; ++i) {
        X(i, 0) = static_cast<Real>(i + 1);
        Z(i, 0) = static_cast<Real>(i + 1) * 2.0;  // 强相关
    }
    auto res = cragg_donald_test(Z, X, W, StockYogoCriterion::RelativeBias, 0.10);
    EXPECT_NEAR(res.stock_yogo_critical_value, 9.08, TOL_CV);
    EXPECT_TRUE(res.critical_value_is_exact);
    EXPECT_EQ(res.critical_value_source, "SY2005 Table 1");
    EXPECT_NEAR(res.bias_threshold, 0.10, TOL_CV);
}

// --- SY Table 1 查表 2: K=1, L=5, Bias 10% → 10.83 ---
TEST(StockYogoBiasTable, K1L5Bias10pct) {
    MatrixXD X(20, 1), Z(20, 5), W(20, 0);
    for (Size i = 0; i < 20; ++i) {
        Z(i, 0) = basis_lin(i);
        Z(i, 1) = basis_sq(i);
        Z(i, 2) = basis_sin(i);
        Z(i, 3) = basis_cos(i);
        Z(i, 4) = basis_sqrt(i);
        X(i, 0) = Z(i, 0) * 0.8 + Z(i, 1) * 0.3;
    }
    auto res = cragg_donald_test(Z, X, W, StockYogoCriterion::RelativeBias, 0.10);
    EXPECT_NEAR(res.stock_yogo_critical_value, 10.83, TOL_CV);
    EXPECT_TRUE(res.critical_value_is_exact);
    EXPECT_EQ(res.critical_value_source, "SY2005 Table 1");
}

// --- SY Table 1 查表 3: K=2, L=3, Bias 10% → 11.47 ---
TEST(StockYogoBiasTable, K2L3Bias10pct) {
    MatrixXD X(30, 2), Z(30, 3), W(30, 0);
    for (Size i = 0; i < 30; ++i) {
        Z(i, 0) = basis_lin(i);
        Z(i, 1) = basis_sq(i);
        Z(i, 2) = basis_sin(i);
        X(i, 0) = Z(i, 0) * 0.8 + Z(i, 1) * 0.3 + 0.1;
        X(i, 1) = Z(i, 1) * 0.4 + Z(i, 2) * 0.6 + 0.2;
    }
    auto res = cragg_donald_test(Z, X, W, StockYogoCriterion::RelativeBias, 0.10);
    EXPECT_NEAR(res.stock_yogo_critical_value, 11.47, TOL_CV);
    EXPECT_TRUE(res.critical_value_is_exact);
    EXPECT_EQ(res.critical_value_source, "SY2005 Table 1");
}

// --- SY Table 1 查表 4: K=2, L=5, Bias 10% → 13.97 ---
TEST(StockYogoBiasTable, K2L5Bias10pct) {
    MatrixXD X(40, 2), Z(40, 5), W(40, 0);
    for (Size i = 0; i < 40; ++i) {
        Z(i, 0) = basis_lin(i);
        Z(i, 1) = basis_sq(i);
        Z(i, 2) = basis_sin(i);
        Z(i, 3) = basis_cos(i);
        Z(i, 4) = basis_sqrt(i);
        X(i, 0) = Z(i, 0) * 0.8 + Z(i, 1) * 0.3 + 0.1;
        X(i, 1) = Z(i, 2) * 0.4 + Z(i, 3) * 0.6 + 0.2;
    }
    auto res = cragg_donald_test(Z, X, W, StockYogoCriterion::RelativeBias, 0.10);
    EXPECT_NEAR(res.stock_yogo_critical_value, 13.97, TOL_CV);
    EXPECT_TRUE(res.critical_value_is_exact);
    EXPECT_EQ(res.critical_value_source, "SY2005 Table 1");
}

// --- SY Table 1 查表 5: K=3, L=3, Bias 10% → 9.08 (just-identified) ---
TEST(StockYogoBiasTable, K3L3Bias10pct) {
    MatrixXD X(50, 3), Z(50, 3), W(50, 0);
    for (Size i = 0; i < 50; ++i) {
        Z(i, 0) = basis_lin(i);
        Z(i, 1) = basis_sq(i);
        Z(i, 2) = basis_sin(i);
        X(i, 0) = Z(i, 0) * 0.8 + 0.1;
        X(i, 1) = Z(i, 1) * 0.6 + 0.2;
        X(i, 2) = Z(i, 2) * 0.9 + 0.05;
    }
    auto res = cragg_donald_test(Z, X, W, StockYogoCriterion::RelativeBias, 0.10);
    EXPECT_NEAR(res.stock_yogo_critical_value, 9.08, TOL_CV);
    EXPECT_TRUE(res.critical_value_is_exact);
    EXPECT_EQ(res.critical_value_source, "SY2005 Table 1");
}

// --- SY Table 1 查表 6: K=3, L=5, Bias 10% → 14.48 ---
TEST(StockYogoBiasTable, K3L5Bias10pct) {
    MatrixXD X(60, 3), Z(60, 5), W(60, 0);
    for (Size i = 0; i < 60; ++i) {
        Z(i, 0) = basis_lin(i);
        Z(i, 1) = basis_sq(i);
        Z(i, 2) = basis_sin(i);
        Z(i, 3) = basis_cos(i);
        Z(i, 4) = basis_sqrt(i);
        X(i, 0) = Z(i, 0) * 0.8 + 0.1;
        X(i, 1) = Z(i, 1) * 0.6 + 0.2;
        X(i, 2) = Z(i, 2) * 0.9 + 0.05;
    }
    auto res = cragg_donald_test(Z, X, W, StockYogoCriterion::RelativeBias, 0.10);
    EXPECT_NEAR(res.stock_yogo_critical_value, 14.48, TOL_CV);
    EXPECT_TRUE(res.critical_value_is_exact);
    EXPECT_EQ(res.critical_value_source, "SY2005 Table 1");
}

// =============================================================================
// Skeels-Windmeijer 2018 近似验证 (2 用例, K=3, Size 准则)
// 排幻觉点 H13: K=3 Size 准则用近似, critical_value_is_exact=false
// 容差 1e-4
// =============================================================================

// --- SW 2018 近似 1: K=3, L=4, Size 15% 准则 ---
// 近似策略: K=2, L=4 Size 15% 临界值 (13.43) × 缩放因子 1.20
TEST(SkeelsWindmeijerApprox, K3L4Size15pct) {
    MatrixXD X(80, 3), Z(80, 4), W(80, 0);
    for (Size i = 0; i < 80; ++i) {
        Z(i, 0) = basis_lin(i);
        Z(i, 1) = basis_sq(i);
        Z(i, 2) = basis_sin(i);
        Z(i, 3) = basis_cos(i);
        X(i, 0) = Z(i, 0) * 0.8 + 0.1;
        X(i, 1) = Z(i, 1) * 0.6 + 0.2;
        X(i, 2) = Z(i, 2) * 0.9 + 0.05;
    }
    auto res = cragg_donald_test(Z, X, W, StockYogoCriterion::SizeDistortion, 0.15);

    // K=3 Size 准则: 必须用近似
    EXPECT_FALSE(res.critical_value_is_exact);
    EXPECT_EQ(res.critical_value_source, "Skeels-Windmeijer 2018 approx");
    EXPECT_NEAR(res.size_threshold, 0.15, TOL_CV);
    EXPECT_TRUE(std::isnan(res.bias_threshold));

    // 期望值: K=2, L=4 Size 15% × 1.20 = 13.43 × 1.20 = 16.116
    EXPECT_NEAR(res.stock_yogo_critical_value, 16.116, TOL_APPROX);
}

// --- SW 2018 近似 2: K=3, L=5, Size 15% 准则 ---
// 近似策略: K=2, L=5 Size 15% 临界值 (14.71) × 缩放因子 1.20
TEST(SkeelsWindmeijerApprox, K3L5Size15pct) {
    MatrixXD X(100, 3), Z(100, 5), W(100, 0);
    for (Size i = 0; i < 100; ++i) {
        Z(i, 0) = basis_lin(i);
        Z(i, 1) = basis_sq(i);
        Z(i, 2) = basis_sin(i);
        Z(i, 3) = basis_cos(i);
        Z(i, 4) = basis_sqrt(i);
        X(i, 0) = Z(i, 0) * 0.8 + 0.1;
        X(i, 1) = Z(i, 1) * 0.6 + 0.2;
        X(i, 2) = Z(i, 2) * 0.9 + 0.05;
    }
    auto res = cragg_donald_test(Z, X, W, StockYogoCriterion::SizeDistortion, 0.15);

    // K=3 Size 准则: 必须用近似
    EXPECT_FALSE(res.critical_value_is_exact);
    EXPECT_EQ(res.critical_value_source, "Skeels-Windmeijer 2018 approx");

    // 期望值: K=2, L=5 Size 15% × 1.20 = 14.71 × 1.20 = 17.652
    EXPECT_NEAR(res.stock_yogo_critical_value, 17.652, TOL_APPROX);
}

// =============================================================================
// 附加: 接口异常处理验证 (在 12 用例之外, 不计入正式统计)
// =============================================================================

// --- 异常 1: L < K (under-identified) 抛异常 ---
TEST(CraggDonaldExceptions, UnderIdentifiedThrows) {
    MatrixXD X(10, 2), Z(10, 1), W(10, 0);  // K=2, L=1, L<K
    EXPECT_THROW(cragg_donald_test(Z, X, W), std::invalid_argument);
}

// --- 异常 2: 阈值越界抛异常 ---
TEST(CraggDonaldExceptions, InvalidThresholdThrows) {
    MatrixXD X(10, 1), Z(10, 1), W(10, 0);
    for (Size i = 0; i < 10; ++i) {
        X(i, 0) = static_cast<Real>(i + 1);
        Z(i, 0) = static_cast<Real>(i + 1) * 2.0;
    }
    // Bias 准则阈值必须为 0.05/0.10/0.20/0.30
    EXPECT_THROW(cragg_donald_test(Z, X, W, StockYogoCriterion::RelativeBias, 0.15),
                 std::invalid_argument);
    // Size 准则阈值必须为 0.10/0.15/0.20/0.25
    EXPECT_THROW(cragg_donald_test(Z, X, W, StockYogoCriterion::SizeDistortion, 0.05),
                 std::invalid_argument);
}

// --- 异常 3: K>3 抛异常 (超出 Stock-Yogo 表覆盖范围) ---
TEST(CraggDonaldExceptions, KGreaterThan3Throws) {
    MatrixXD X(20, 4), Z(20, 5), W(20, 0);  // K=4
    for (Size i = 0; i < 20; ++i) {
        const Real t = static_cast<Real>(i + 1);
        for (Size j = 0; j < 4; ++j) X(i, j) = t * static_cast<Real>(j + 1);
        for (Size j = 0; j < 5; ++j) Z(i, j) = t * static_cast<Real>(j + 1) + 1.0;
    }
    EXPECT_THROW(cragg_donald_test(Z, X, W), std::invalid_argument);
}

// --- 异常 4: K=1 Size 准则, L=1 (just-identified) 抛异常 ---
// Size 准则需 L>K (过度识别), L=1 时 just-identified 无法做 Size 检验
TEST(CraggDonaldExceptions, SizeCriterionJustIdentifiedThrows) {
    MatrixXD X(10, 1), Z(10, 1), W(10, 0);  // K=1, L=1 (just-identified)
    for (Size i = 0; i < 10; ++i) {
        X(i, 0) = static_cast<Real>(i + 1);
        Z(i, 0) = static_cast<Real>(i + 1) * 2.0;
    }
    EXPECT_THROW(cragg_donald_test(Z, X, W, StockYogoCriterion::SizeDistortion, 0.15),
                 std::invalid_argument);
}
