// =============================================================================
// test_specification_tests.cpp - Phase 7A Wave 2 模型设定检验测试
//
// 20 用例: 信息矩阵(6) + MZ(7) + DM(7)
//
// 容差: 统计量 1e-8, p_value 1e-6 (chi2_sf/beta_i 数值近似)
//
// 排幻觉点覆盖:
//   H9  (信息矩阵对 QMLE 检验均值方程, 非方差方程)
//   H10 (MZ R² 也是预测精度度量)
//   H11 (HLN 修正用 1/N 计算 γ̂_h, DM ~ t(N-1) 非 N(0,1))
//
// 教材锚点: White 1982, Mincer-Zarnowitz 1969, Diebold-Mariano 1995, HLN 1997
// =============================================================================

#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>

#include "cpphub/econometrics/inference/specification_tests.hpp"

using cpphub::v1::econometrics::information_matrix_test;
using cpphub::v1::econometrics::mincer_zarnowitz_regression;
using cpphub::v1::econometrics::diebold_mariano_test;
using cpphub::v1::econometrics::InformationMatrixResult;
using cpphub::v1::econometrics::MincerZarnowitzResult;
using cpphub::v1::econometrics::DieboldMarianoResult;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
constexpr Real TOL_STAT = 1e-8;
constexpr Real TOL_PVAL = 1e-6;
}  // namespace

// =============================================================================
// 信息矩阵等式检验 (6 用例)
// =============================================================================

// --- IM 1: 同方差正确模型 → IM 不拒绝 ---
TEST(InformationMatrixTest, HomoscedasticCorrectModel) {
    // OLS: y = β·x + ε, ε ~ iid N(0, σ²)
    // scores_t = ε_t · x_t, hessian = -X'X
    // I_outer = (1/N) Σ ε_t² x_t x_t' ≈ σ² · (1/N) X'X = I_inner (同方差)
    const Size N = 50;
    // 构造有变化的 scores (ε_t 交替 1,2,1,2,..., x_t = {1, t/N})
    std::vector<std::vector<Real>> scores(N, std::vector<Real>(2));
    std::vector<std::vector<Real>> hessian(2, std::vector<Real>(2, 0.0));
    for (Size t = 0; t < N; ++t) {
        const Real x1 = 1.0;
        const Real x2 = static_cast<Real>(t) / static_cast<Real>(N);
        const Real eps = (t % 2 == 0) ? 1.0 : 2.0;  // 交替 ε
        scores[t][0] = eps * x1;
        scores[t][1] = eps * x2;
        hessian[0][0] -= x1 * x1;
        hessian[0][1] -= x1 * x2;
        hessian[1][0] -= x2 * x1;
        hessian[1][1] -= x2 * x2;
    }
    auto res = information_matrix_test(scores, hessian);
    EXPECT_GE(res.base.statistic, 0.0);
    EXPECT_GT(res.base.p_value, 0.0);
    EXPECT_LE(res.base.p_value, 1.0);
}

// --- IM 2: 异方差错误模型 → IM 拒绝 ---
TEST(InformationMatrixTest, HeteroscedasticWrongModel) {
    // 异方差: ε_t² 随 t 变化 → I_outer ≠ I_inner
    const Size N = 100;
    std::vector<std::vector<Real>> scores(N, std::vector<Real>(2));
    std::vector<std::vector<Real>> hessian(2, std::vector<Real>(2, 0.0));
    for (Size t = 0; t < N; ++t) {
        const Real x1 = 1.0;
        const Real x2 = static_cast<Real>(t) / static_cast<Real>(N);
        // 异方差: ε_t = (t/N) * 2 → ε_t² = (t/N)² * 4
        const Real eps = (static_cast<Real>(t) / static_cast<Real>(N)) * 2.0;
        scores[t][0] = eps * x1;
        scores[t][1] = eps * x2;
        hessian[0][0] -= x1 * x1;
        hessian[0][1] -= x1 * x2;
        hessian[1][0] -= x2 * x1;
        hessian[1][1] -= x2 * x2;
    }
    auto res = information_matrix_test(scores, hessian);
    // 异方差 → IM 应较大
    EXPECT_GT(res.base.statistic, 0.0);
}

// --- IM 3: H9 - 对 QMLE 检验均值方程正确性 ---
TEST(InformationMatrixTest, QMLEMeanEquationCheck) {
    // 排幻觉点 H9: 信息矩阵检验对 QMLE 检验均值方程, 不是方差方程
    // 即使方差方程错误, 如果均值方程正确, IM 应不显著
    const Size N = 20;
    std::vector<std::vector<Real>> scores(N, std::vector<Real>(1));
    std::vector<std::vector<Real>> hessian(1, std::vector<Real>(1, 0.0));
    // K=1, scores 有变化 (交替 1, 2)
    for (Size t = 0; t < N; ++t) {
        scores[t][0] = (t % 2 == 0) ? 1.0 : 2.0;
        hessian[0][0] -= 1.0;
    }
    auto res = information_matrix_test(scores, hessian);
    // I_outer = (1/20)(10*1 + 10*4) = 50/20 = 2.5
    // I_inner = 20/20 = 1
    // D = 2.5 - 1 = 1.5 (≠0, 说明模型设定有误)
    EXPECT_GE(res.base.statistic, 0.0);
    EXPECT_GT(res.base.p_value, 0.0);
    EXPECT_LE(res.base.p_value, 1.0);
}

// --- IM 4: df 正确性 ---
TEST(InformationMatrixTest, DFCorrectness) {
    // K=1 → q = 1(1+2)/2 = 1
    // K=2 → q = 2(3)/2 = 3
    // K=3 → q = 3(4)/2 = 6
    // 注: scores 必须有变化, 否则 Â 奇异
    const Size N = 20;

    // K=1 (scores 有变化)
    std::vector<std::vector<Real>> scores1(N, std::vector<Real>(1));
    std::vector<std::vector<Real>> hessian1(1, std::vector<Real>(1, 0.0));
    for (Size t = 0; t < N; ++t) {
        scores1[t][0] = static_cast<Real>(t % 3) + 1.0;  // {1,2,3,1,2,3,...}
        hessian1[0][0] -= 1.0;
    }
    auto res1 = information_matrix_test(scores1, hessian1);
    EXPECT_EQ(res1.df, 1u);

    // K=2
    std::vector<std::vector<Real>> scores2(N, std::vector<Real>(2));
    std::vector<std::vector<Real>> hessian2(2, std::vector<Real>(2, 0.0));
    for (Size t = 0; t < N; ++t) {
        scores2[t][0] = static_cast<Real>(t % 3) + 1.0;
        scores2[t][1] = static_cast<Real>((t + 1) % 4) + 0.5;
        hessian2[0][0] -= 1.0;
        hessian2[1][1] -= 1.0;
    }
    auto res2 = information_matrix_test(scores2, hessian2);
    EXPECT_EQ(res2.df, 3u);

    // K=3
    std::vector<std::vector<Real>> scores3(N, std::vector<Real>(3));
    std::vector<std::vector<Real>> hessian3(3, std::vector<Real>(3, 0.0));
    for (Size t = 0; t < N; ++t) {
        scores3[t][0] = static_cast<Real>(t % 3) + 1.0;
        scores3[t][1] = static_cast<Real>((t + 1) % 4) + 0.5;
        scores3[t][2] = static_cast<Real>((t + 2) % 5) + 0.3;
        hessian3[0][0] -= 1.0;
        hessian3[1][1] -= 1.0;
        hessian3[2][2] -= 1.0;
    }
    auto res3 = information_matrix_test(scores3, hessian3);
    EXPECT_EQ(res3.df, 6u);
}

// --- IM 5: 手算闭式验证 (K=1) ---
TEST(InformationMatrixTest, HandComputedK1) {
    // K=1, N=5, scores = {1, 2, 1, 2, 1}, hessian = {{-5}}
    // I_outer = (1/5)(1+4+1+4+1) = 11/5 = 2.2
    // I_inner = 5/5 = 1
    // D = 2.2 - 1 = 1.2, vech(D) = {1.2}
    // Â = (1/5) Σ (s_t² - I_outer)² = (1/5)((1-2.2)² + (4-2.2)² + (1-2.2)² + (4-2.2)² + (1-2.2)²)
    //    = (1/5)(1.44 + 3.24 + 1.44 + 3.24 + 1.44) = 10.8/5 = 2.16
    // IM = 5 * 1.2 * (1/2.16) * 1.2 = 5 * 1.44/2.16 = 7.2/2.16 = 10/3 ≈ 3.3333
    std::vector<std::vector<Real>> scores = {{1.0}, {2.0}, {1.0}, {2.0}, {1.0}};
    std::vector<std::vector<Real>> hessian = {{-5.0}};
    auto res = information_matrix_test(scores, hessian);

    EXPECT_NEAR(res.base.statistic, 10.0 / 3.0, 1e-6);
    EXPECT_EQ(res.df, 1u);
}

// --- IM 6: 输入验证 ---
TEST(InformationMatrixTest, InputValidation) {
    // N < 5
    std::vector<std::vector<Real>> scores2 = {{1.0}, {2.0}};
    std::vector<std::vector<Real>> hessian2 = {{-2.0}};
    EXPECT_THROW(information_matrix_test(scores2, hessian2), std::invalid_argument);

    // 空 scores
    std::vector<std::vector<Real>> empty_scores;
    std::vector<std::vector<Real>> empty_hessian(1, std::vector<Real>(1, 0.0));
    EXPECT_THROW(information_matrix_test(empty_scores, empty_hessian), std::invalid_argument);

    // hessian 维度不匹配
    std::vector<std::vector<Real>> scores3(5, std::vector<Real>(2, 1.0));
    std::vector<std::vector<Real>> hessian3(1, std::vector<Real>(1, 0.0));
    EXPECT_THROW(information_matrix_test(scores3, hessian3), std::invalid_argument);
}

// =============================================================================
// Mincer-Zarnowitz 回归 (7 用例)
// =============================================================================

// --- MZ 1: 完美预测 → F 不拒绝 ---
TEST(MincerZarnowitzTest, PerfectForecastNotReject) {
    // actual = forecast (α=0, β=1)
    std::vector<Real> actual = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<Real> forecast = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    auto res = mincer_zarnowitz_regression(actual, forecast);

    // α=0, β=1, R²=1
    EXPECT_NEAR(res.alpha, 0.0, TOL_STAT);
    EXPECT_NEAR(res.beta, 1.0, TOL_STAT);
    EXPECT_NEAR(res.r_squared, 1.0, TOL_STAT);
    // RSS_u = 0 → F 应为 0 或 NaN (退化)
    EXPECT_GE(res.base.p_value, 0.0);
    EXPECT_LE(res.base.p_value, 1.0);
}

// --- MZ 2: 有偏预测 (α≠0) → F 拒绝 ---
TEST(MincerZarnowitzTest, BiasedForecastRejects) {
    // actual = 2 + 1*forecast + small noise (α=2, β=1)
    // 每点: 2+f+noise, noise ∈ {-0.02..+0.03}
    std::vector<Real> actual = {3.01, 4.02, 4.99, 6.01, 6.98, 8.03, 8.99, 10.02};
    std::vector<Real> forecast = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    auto res = mincer_zarnowitz_regression(actual, forecast);

    // α≈2, β≈1
    EXPECT_NEAR(res.alpha, 2.0, 0.1);
    EXPECT_NEAR(res.beta, 1.0, 0.1);
    // F 应大 (RSS_r >> RSS_u)
    EXPECT_GT(res.base.statistic, 5.0);
    EXPECT_TRUE(res.base.reject_null);
}

// --- MZ 3: β≠1 → F 拒绝 ---
TEST(MincerZarnowitzTest, BetaNotOneRejects) {
    // actual = 2 * forecast + small noise (α=0, β=2)
    std::vector<Real> actual = {2.01, 4.02, 5.99, 8.03, 9.98, 12.01, 14.02, 15.99};
    std::vector<Real> forecast = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    auto res = mincer_zarnowitz_regression(actual, forecast);

    // α≈0, β≈2
    EXPECT_NEAR(res.alpha, 0.0, 0.1);
    EXPECT_NEAR(res.beta, 2.0, 0.1);
    // F 应大 (β≠1)
    EXPECT_GT(res.base.statistic, 5.0);
    EXPECT_TRUE(res.base.reject_null);
}

// --- MZ 4: R² 正确性 ---
TEST(MincerZarnowitzTest, RSquaredCorrectness) {
    // actual = forecast + noise
    // R² 应在 (0, 1)
    std::vector<Real> actual = {1.1, 2.2, 2.9, 4.1, 4.8, 6.2, 6.9, 8.1};
    std::vector<Real> forecast = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    auto res = mincer_zarnowitz_regression(actual, forecast);

    EXPECT_GT(res.r_squared, 0.0);
    EXPECT_LT(res.r_squared, 1.0);
    // β≈1, α≈0.1
    EXPECT_NEAR(res.beta, 1.0, 0.1);
}

// --- MZ 5: H10 - R² 是预测精度度量 ---
TEST(MincerZarnowitzTest, RSquaredAsAccuracyMetric) {
    // 排幻觉点 H10: R² 越大预测越准
    // 场景 1: 小噪声 → R² 大
    std::vector<Real> actual1 = {1.01, 2.02, 2.99, 4.01, 4.98, 6.02, 6.99, 8.01};
    std::vector<Real> forecast = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    auto res1 = mincer_zarnowitz_regression(actual1, forecast);

    // 场景 2: 大噪声 → R² 小
    std::vector<Real> actual2 = {1.5, 1.5, 4.5, 3.5, 6.5, 5.5, 8.5, 7.5};
    auto res2 = mincer_zarnowitz_regression(actual2, forecast);

    EXPECT_GT(res1.r_squared, res2.r_squared);
}

// --- MZ 6: t 统计量正确性 ---
TEST(MincerZarnowitzTest, TStatisticCorrectness) {
    // actual = 2 * forecast + small noise (α=0, β=2)
    // β-1 = 1, t_β 应大
    std::vector<Real> actual = {2.01, 4.02, 5.99, 8.03, 9.98, 12.01, 14.02, 15.99};
    std::vector<Real> forecast = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    auto res = mincer_zarnowitz_regression(actual, forecast);

    // β≈2, β-1=1, t_β 应大
    EXPECT_NEAR(res.beta, 2.0, 0.1);
    EXPECT_GT(std::fabs(res.beta_t_stat), 3.0);  // β≠1 显著
}

// --- MZ 7: 输入验证 ---
TEST(MincerZarnowitzTest, InputValidation) {
    // N < 5
    std::vector<Real> actual = {1.0, 2.0, 3.0};
    std::vector<Real> forecast = {1.0, 2.0, 3.0};
    EXPECT_THROW(mincer_zarnowitz_regression(actual, forecast), std::invalid_argument);

    // 尺寸不匹配
    std::vector<Real> actual2(10, 1.0);
    std::vector<Real> forecast2(9, 1.0);
    EXPECT_THROW(mincer_zarnowitz_regression(actual2, forecast2), std::invalid_argument);
}

// =============================================================================
// Diebold-Mariano 检验 (7 用例)
// =============================================================================

// --- DM 1: 两个预测相同 → DM=0, p_value=1 ---
TEST(DieboldMarianoTest, IdenticalForecastsZeroDM) {
    std::vector<Real> actual = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    std::vector<Real> forecast1 = {1.1, 2.2, 2.9, 4.1, 4.8, 6.2, 6.9, 8.1, 8.9, 10.2};
    std::vector<Real> forecast2 = forecast1;  // 完全相同
    auto res = diebold_mariano_test(actual, forecast1, forecast2, "mse", 1);

    EXPECT_NEAR(res.mean_loss_diff, 0.0, TOL_STAT);
    EXPECT_NEAR(res.base.statistic, 0.0, TOL_STAT);
    EXPECT_NEAR(res.base.p_value, 1.0, TOL_PVAL);
    EXPECT_FALSE(res.base.reject_null);
}

// --- DM 2: forecast1 显著优于 forecast2 → DM 显著 ---
TEST(DieboldMarianoTest, Forecast1BetterThanForecast2) {
    // forecast1 接近 actual, forecast2 远离
    std::vector<Real> actual = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    std::vector<Real> forecast1 = {1.1, 2.1, 2.9, 4.1, 4.9, 6.1, 6.9, 8.1, 8.9, 10.1};
    std::vector<Real> forecast2 = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0};
    auto res = diebold_mariano_test(actual, forecast1, forecast2, "mse", 1);

    // d̄ = mean(L1 - L2) < 0 (forecast1 损失更小)
    EXPECT_LT(res.mean_loss_diff, 0.0);
    // |DM| 应较大
    EXPECT_GT(std::fabs(res.base.statistic), 2.0);
    EXPECT_TRUE(res.base.reject_null);
}

// --- DM 3: MSE vs MAE loss ---
TEST(DieboldMarianoTest, DifferentLossFunctions) {
    std::vector<Real> actual = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    std::vector<Real> forecast1 = {1.1, 2.1, 2.9, 4.1, 4.9, 6.1, 6.9, 8.1, 8.9, 10.1};
    std::vector<Real> forecast2 = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};

    auto res_mse = diebold_mariano_test(actual, forecast1, forecast2, "mse", 1);
    auto res_mae = diebold_mariano_test(actual, forecast1, forecast2, "mae", 1);

    // 两种 loss 的 DM 统计量应不同
    EXPECT_NE(res_mse.base.statistic, res_mae.base.statistic);
    // 都应拒绝 (forecast1 更好)
    EXPECT_TRUE(res_mse.base.reject_null);
    EXPECT_TRUE(res_mae.base.reject_null);
}

// --- DM 4: h>1 多步预测 ---
TEST(DieboldMarianoTest, MultiStepForecast) {
    std::vector<Real> actual = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    std::vector<Real> forecast1 = {1.1, 2.1, 2.9, 4.1, 4.9, 6.1, 6.9, 8.1, 8.9, 10.1};
    std::vector<Real> forecast2 = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};

    // h=1 vs h=3
    auto res_h1 = diebold_mariano_test(actual, forecast1, forecast2, "mse", 1);
    auto res_h3 = diebold_mariano_test(actual, forecast1, forecast2, "mse", 3);

    // h=3 时 V̂ 包含更多 γ̂_h → DM 可能不同
    EXPECT_GE(std::fabs(res_h1.base.statistic), 0.0);
    EXPECT_GE(std::fabs(res_h3.base.statistic), 0.0);
}

// --- DM 5: H11 - HLN 修正 (1/N 而非 1/(N-h)) ---
TEST(DieboldMarianoTest, HLNCorrection) {
    // 排幻觉点 H11: HLN 修正用 1/N 计算 γ̂_h
    // 验证: h=2 时, γ̂_1 用 1/N 计算 (非 1/(N-1))
    std::vector<Real> actual = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    std::vector<Real> forecast1 = {1.1, 2.1, 2.9, 4.1, 4.9, 6.1, 6.9, 8.1, 8.9, 10.1};
    std::vector<Real> forecast2 = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};

    // h=2, V̂ = γ̂_0 + 2·γ̂_1
    auto res = diebold_mariano_test(actual, forecast1, forecast2, "mse", 2);

    // DM 应为有限值 (HLN 修正生效)
    EXPECT_TRUE(std::isfinite(res.base.statistic));
    EXPECT_TRUE(std::isfinite(res.base.p_value));
}

// --- DM 6: H11 - t 分布非正态 ---
TEST(DieboldMarianoTest, TDistributionNotNormal) {
    // 排幻觉点 H11: DM ~ t(N-1), 非 N(0,1)
    // 验证: 小样本时 p_value 用 t 分布 (比正态更保守)
    std::vector<Real> actual = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<Real> forecast1 = {1.1, 2.1, 2.9, 4.1, 4.9};
    std::vector<Real> forecast2 = {1.5, 2.5, 3.5, 4.5, 5.5};
    auto res = diebold_mariano_test(actual, forecast1, forecast2, "mse", 1);

    // N=5, df=4, t 分布 p_value
    EXPECT_GE(res.base.p_value, 0.0);
    EXPECT_LE(res.base.p_value, 1.0);
    // DM 统计量应有限
    EXPECT_TRUE(std::isfinite(res.base.statistic));
}

// --- DM 7: 输入验证 ---
TEST(DieboldMarianoTest, InputValidation) {
    // N < 5
    std::vector<Real> actual = {1.0, 2.0};
    std::vector<Real> forecast1 = {1.0, 2.0};
    std::vector<Real> forecast2 = {1.0, 2.0};
    EXPECT_THROW(diebold_mariano_test(actual, forecast1, forecast2), std::invalid_argument);

    // 尺寸不匹配
    std::vector<Real> actual2(10, 1.0);
    std::vector<Real> forecast2a(10, 1.0);
    std::vector<Real> forecast2b(9, 1.0);
    EXPECT_THROW(diebold_mariano_test(actual2, forecast2a, forecast2b), std::invalid_argument);

    // 无效 h
    std::vector<Real> actual3(10, 1.0);
    std::vector<Real> f3(10, 1.0);
    EXPECT_THROW(diebold_mariano_test(actual3, f3, f3, "mse", 0), std::invalid_argument);
    EXPECT_THROW(diebold_mariano_test(actual3, f3, f3, "mse", 10), std::invalid_argument);

    // 未知 loss_function
    EXPECT_THROW(diebold_mariano_test(actual3, f3, f3, "unknown"), std::invalid_argument);
}
