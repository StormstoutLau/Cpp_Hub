// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.1 任务 2.1-2.2 - Gaussian MLE 测试
// 验证方法: 手算解析值 + 与 OLS 闭式解对齐 + 协方差三种形式渐近等价 (容差 1e-10)
//
// 排幻觉点:
//   G1: Gaussian MLE 系数 = OLS 闭式解 (identity link, 闭式解, 不走 Newton-Raphson)
//   G2: Gaussian MLE σ² = SSR/N (MLE 有偏估计), 非 SSR/(N-K) (OLS 无偏估计)
//   G3: Gaussian Hessian V = σ²·(X'X)^{-1} (含 σ² 因子), 非 (X'X)^{-1}
//        排幻觉: GLM bread = (X'WX)^{-1}, 但 Gaussian 的 W=1 不含 σ²,
//        需在协方差计算后乘 σ² (R sandwich::bread.glm 乘 dispersion)
//   G4: Gaussian Sandwich V = (X'X)^{-1} · X'diag(ε²)X · (X'X)^{-1} (= HC0 White 1980)
//        排幻觉: meat 用 raw residual ε=y-μ (非 Pearson residual), 无 1/σ² 因子
//   G5: 大样本下 Hessian ≈ Sandwich ≈ OPG (三者渐近等价, White 1982)
//
// 手算数据集 (复用 test_ols_hc.cpp N=4, K=2):
//   X = [1,1; 1,2; 1,3; 1,4],  y = [2; 3; 5; 7]
//   β = [0; 1.7]  (OLS 闭式解 = Gaussian MLE)
//   SSR = 0.30,  SST = 14.75,  ȳ = 4.25
//   σ²_MLE = SSR/N = 0.30/4 = 0.075
//   σ²_OLS = SSR/(N-K) = 0.30/2 = 0.15
//   (X'X)^{-1} = [3/2, -1/2; -1/2, 1/5]
//   Hessian V = σ²_MLE · (X'X)^{-1}:
//     V(0,0) = 0.075 · 1.5 = 0.1125
//     V(1,1) = 0.075 · 0.2 = 0.015
//     V(0,1) = 0.075 · (-0.5) = -0.0375
//   loglik (Gaussian MLE):
//     ℓ = -0.5·N·[log(2π) + log(σ²) + 1]
//       = -0.5·4·[1.8378770664 + (-2.5902671654) + 1]
//       = -2 · 0.2476099010
//       = -0.4952198020
//   残差: ε = y - Xβ = [2-1.7; 3-3.4; 5-5.1; 7-6.8] = [0.3; -0.4; -0.1; 0.2]
//   注意: 0.3²+0.4²+0.1²+0.2² = 0.09+0.16+0.01+0.04 = 0.30 = SSR ✓
//   Sandwich V (HC0):
//     X'diag(ε²)X = Σ ε_i² · x_i·x_i'
//     ε² = [0.09, 0.16, 0.01, 0.04]
//     x_1·x_1' = [1,1; 1,1],  x_2·x_2' = [1,2; 2,4]
//     x_3·x_3' = [1,3; 3,9],  x_4·x_4' = [1,4; 4,16]
//     X'diag(ε²)X = 0.09·[1,1;1,1] + 0.16·[1,2;2,4] + 0.01·[1,3;3,9] + 0.04·[1,4;4,16]
//               = [0.09+0.16+0.01+0.04, 0.09+0.32+0.03+0.16;
//                  0.09+0.32+0.03+0.16, 0.09+0.64+0.09+0.64]
//               = [0.30, 0.60; 0.60, 1.46]
//     (X'X)^{-1} = [1.5, -0.5; -0.5, 0.2]
//     V_Sandwich = (X'X)^{-1} · X'diag(ε²)X · (X'X)^{-1}
//               = [1.5, -0.5; -0.5, 0.2] · [0.30, 0.60; 0.60, 1.46] · [1.5, -0.5; -0.5, 0.2]
//     中间矩阵 A = (X'X)^{-1} · X'diag(ε²)X:
//       A(0,0) = 1.5·0.30 + (-0.5)·0.60 = 0.45 - 0.30 = 0.15
//       A(0,1) = 1.5·0.60 + (-0.5)·1.46 = 0.90 - 0.73 = 0.17
//       A(1,0) = -0.5·0.30 + 0.2·0.60 = -0.15 + 0.12 = -0.03
//       A(1,1) = -0.5·0.60 + 0.2·1.46 = -0.30 + 0.292 = -0.008
//     V = A · (X'X)^{-1}:
//       V(0,0) = 0.15·1.5 + 0.17·(-0.5) = 0.225 - 0.085 = 0.140
//       V(0,1) = 0.15·(-0.5) + 0.17·0.2 = -0.075 + 0.034 = -0.041
//       V(1,0) = -0.03·1.5 + (-0.008)·(-0.5) = -0.045 + 0.004 = -0.041
//       V(1,1) = -0.03·(-0.5) + (-0.008)·0.2 = 0.015 - 0.0016 = 0.0134
//     即 V_Sandwich = [0.140, -0.041; -0.041, 0.0134]  (HC0 形式)
//
// 注: Hessian V (G3) = 0.075·(X'X)^{-1} = [0.1125, -0.0375; -0.0375, 0.015]
//     与 Sandwich V [0.140, -0.041; -0.041, 0.0134] 不等 (因 N=4 小样本)
//     大样本下二者渐近等价 (因 E[ε²]=σ², X'diag(ε²)X ≈ σ²·X'X)

#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>

#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/covariance_type.hpp"
#include "cpphub/econometrics/core/data_types.hpp"
#include "cpphub/econometrics/core/estimation_result.hpp"
#include "cpphub/econometrics/core/estimator_base.hpp"
#include "cpphub/econometrics/estimation/mle.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;
using cpphub::v1::linalg::dynamic::MatrixXD;
using cpphub::v1::linalg::dynamic::VectorXD;

namespace {

// 手算数据集 N=4, K=2
MatrixXD make_X() {
    MatrixXD X(4, 2);
    X(0, 0) = 1.0; X(0, 1) = 1.0;
    X(1, 0) = 1.0; X(1, 1) = 2.0;
    X(2, 0) = 1.0; X(2, 1) = 3.0;
    X(3, 0) = 1.0; X(3, 1) = 4.0;
    return X;
}

VectorXD make_y() {
    VectorXD y(4);
    y(0) = 2.0; y(1) = 3.0; y(2) = 5.0; y(3) = 7.0;
    return y;
}

CrossSectionData make_data() {
    return make_cross_section(make_X(), make_y(), {"intercept", "x"}, "y");
}

}  // namespace

// =============================================================================
// §1 Gaussian MLE 系数 = OLS 闭式解 (排幻觉点 G1)
// =============================================================================
TEST(MLEGaussianTest, Coefficients_Equals_OLS_ClosedForm) {
    MLEEstimator mle(MLEFamily::Gaussian, CovarianceType::Hessian);
    EstimationResult r = mle.estimate(make_data());

    // β = [0; 1.7] (手算闭式解)
    ASSERT_EQ(r.coefficients.size(), 2u);
    EXPECT_NEAR(r.coefficients(0), 0.0, 1e-12);
    EXPECT_NEAR(r.coefficients(1), 1.7, 1e-12);

    // 收敛标志: Gaussian 走闭式解, 应标记为收敛, n_iter=1
    EXPECT_TRUE(mle.converged());
    EXPECT_EQ(mle.nIterations(), 1u);

    // 字段完整性
    EXPECT_EQ(r.n_obs, 4u);
    EXPECT_EQ(r.n_params, 2u);
    EXPECT_EQ(r.df_residual, 2u);
    EXPECT_EQ(r.cov_type, CovarianceType::Hessian);
}

// =============================================================================
// §2 Gaussian MLE 对数似然 (手算, 排幻觉点 G2: σ²=SSR/N)
// =============================================================================
TEST(MLEGaussianTest, LogLikelihood_HandComputed) {
    MLEEstimator mle(MLEFamily::Gaussian, CovarianceType::Hessian);
    EstimationResult r = mle.estimate(make_data());

    // σ²_MLE = 0.075, ℓ = -0.4952198020 (手算, 见文件头注释)
    EXPECT_NEAR(r.log_likelihood, -0.4952198020, 1e-10);
}

// =============================================================================
// §3 Gaussian MLE R² = 1 - SSR/SST (手算)
// =============================================================================
TEST(MLEGaussianTest, RSquared_HandComputed) {
    MLEEstimator mle(MLEFamily::Gaussian, CovarianceType::Hessian);
    EstimationResult r = mle.estimate(make_data());

    // SST = Σ(y-ȳ)² = 14.75, SSR = 0.30, R² = 1 - 0.30/14.75
    const Real expected_r2 = 1.0 - 0.30 / 14.75;
    EXPECT_NEAR(r.r_squared, expected_r2, 1e-12);

    // adj R² = 1 - (1-R²)·(N-1)/(N-K) = 1 - (1-R²)·3/2
    const Real expected_adj_r2 = 1.0 - (1.0 - expected_r2) * 3.0 / 2.0;
    EXPECT_NEAR(r.adj_r_squared, expected_adj_r2, 1e-12);
}

// =============================================================================
// §4 Gaussian MLE Hessian 协方差 = σ²·(X'X)^{-1} (排幻觉点 G3)
// =============================================================================
TEST(MLEGaussianTest, HessianVcov_Equals_Sigma2_XtX_inv) {
    MLEEstimator mle(MLEFamily::Gaussian, CovarianceType::Hessian);
    EstimationResult r = mle.estimate(make_data());

    // 排幻觉点 G3: Hessian V = σ²_MLE · (X'X)^{-1}
    //   σ²_MLE = 0.075
    //   (X'X)^{-1} = [1.5, -0.5; -0.5, 0.2]
    //   V = [0.1125, -0.0375; -0.0375, 0.015]
    ASSERT_EQ(r.vcov.rows(), 2u);
    ASSERT_EQ(r.vcov.cols(), 2u);
    EXPECT_NEAR(r.vcov(0, 0), 0.1125, 1e-12);
    EXPECT_NEAR(r.vcov(0, 1), -0.0375, 1e-12);
    EXPECT_NEAR(r.vcov(1, 0), -0.0375, 1e-12);
    EXPECT_NEAR(r.vcov(1, 1), 0.015, 1e-12);

    // 标准误 = sqrt(diag(V))
    ASSERT_EQ(r.std_errors.size(), 2u);
    EXPECT_NEAR(r.std_errors(0), std::sqrt(0.1125), 1e-12);
    EXPECT_NEAR(r.std_errors(1), std::sqrt(0.015), 1e-12);
}

// =============================================================================
// §5 Gaussian MLE Sandwich 协方差 = HC0 White 1980 (排幻觉点 G4)
// =============================================================================
TEST(MLEGaussianTest, SandwichVcov_Equals_HC0_White1980) {
    MLEEstimator mle(MLEFamily::Gaussian, CovarianceType::Sandwich);
    EstimationResult r = mle.estimate(make_data());

    // 排幻觉点 G4: V = (X'X)^{-1} · X'diag(ε²)X · (X'X)^{-1}
    // 手算值见文件头: V = [0.140, -0.041; -0.041, 0.0134]
    ASSERT_EQ(r.vcov.rows(), 2u);
    ASSERT_EQ(r.vcov.cols(), 2u);
    EXPECT_NEAR(r.vcov(0, 0), 0.140, 1e-10);
    EXPECT_NEAR(r.vcov(0, 1), -0.041, 1e-10);
    EXPECT_NEAR(r.vcov(1, 0), -0.041, 1e-10);
    EXPECT_NEAR(r.vcov(1, 1), 0.0134, 1e-10);
}

// =============================================================================
// §6 Gaussian MLE 通过 Estimator 基类多态调用 (API 一致性)
// =============================================================================
TEST(MLEGaussianTest, PolymorphicCallViaEstimatorBase) {
    std::unique_ptr<Estimator> est = std::make_unique<MLEEstimator>(
        MLEFamily::Gaussian, CovarianceType::Hessian);
    EXPECT_EQ(est->name(), "MLE(Gaussian)");
    EXPECT_EQ(est->estimatorClass(), EstimatorClass::Parametric);
    EXPECT_TRUE(est->isParametric());
    EXPECT_EQ(est->covarianceType(), CovarianceType::Hessian);

    EstimationResult r = est->estimate(make_data());
    EXPECT_NEAR(r.coefficients(1), 1.7, 1e-12);

    // clone 应返回独立的等价对象
    auto cloned = est->clone();
    EXPECT_EQ(cloned->name(), "MLE(Gaussian)");
    EstimationResult r2 = cloned->estimate(make_data());
    EXPECT_NEAR(r2.coefficients(1), 1.7, 1e-12);
}

// =============================================================================
// §7 Gaussian MLE 异常输入校验
// =============================================================================
TEST(MLEGaussianTest, ThrowsOnInsufficientObservations) {
    // N <= K: 不可识别
    MatrixXD X(2, 2);
    X(0, 0) = 1.0; X(0, 1) = 1.0;
    X(1, 0) = 1.0; X(1, 1) = 2.0;
    VectorXD y(2);
    y(0) = 1.0; y(1) = 2.0;
    auto data = make_cross_section(X, y, {"a", "b"}, "y");

    MLEEstimator mle(MLEFamily::Gaussian);
    EXPECT_THROW(mle.estimate(data), std::invalid_argument);
}

TEST(MLEGaussianTest, ThrowsOnUnsupportedCovarianceType) {
    // MLE 不支持 HC0/HC1/Cluster 等
    MLEEstimator mle(MLEFamily::Gaussian, CovarianceType::HC0);
    EXPECT_THROW(mle.estimate(make_data()), std::invalid_argument);
}

// =============================================================================
// §8 Gaussian MLE 不支持的分布族数据传入 PanelData 抛异常
// =============================================================================
TEST(MLEGaussianTest, ThrowsOnPanelData) {
    // 构造一个空 PanelData (MLEEstimator 仅支持 CrossSectionData)
    PanelData pd;
    pd.X = MatrixXD(2, 2);
    pd.y = VectorXD(2);

    MLEEstimator mle(MLEFamily::Gaussian);
    EconData data = pd;
    EXPECT_THROW(mle.estimate(data), std::invalid_argument);
}
