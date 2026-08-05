// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.1 任务 2.1-2.2 - Poisson MLE 测试
// 验证方法: link/mean 互逆 + IRLS weight + Newton-Raphson 收敛 + 手算解析解 (容差 1e-10)
//
// 排幻觉点:
//   P1: Poisson μ = exp(η) > 0 (始终非负), η = log(μ)
//   P2: IRLS weight w = μ (canonical log link, V(μ) = μ)
//   P3: ℓ_i = y·η - μ - log(y!) = y·η - exp(η) - lgamma(y+1)
//   P4: score = X'(y-μ), Hessian = -X'WX, W = diag(μ) (canonical)
//   P5: Poisson MLE 等价于 Poisson GLM (R glm(..., family=poisson)), IRLS = NR
//
// 手算数据集 A (无截距, N=2, K=1):
//   x = [1, 1],  y = [2, 4]
//   MLE: μ = ȳ = 3, β = log(3) ≈ 1.0986122886681098
//   W = [3, 3], X'WX = 3+3 = 6
//   Hessian V = (X'WX)^{-1} = 1/6 ≈ 0.16667
//   ε = [2-3, 4-3] = [-1, 1], ε² = [1, 1]
//   X'diag(ε²)X = 1+1 = 2
//   Sandwich V = (1/6)·2·(1/6) = 2/36 = 1/18 ≈ 0.05556
//   OPG V = (X'diag(ε²)X)^{-1} = 1/2 = 0.5
//   ℓ = 2·log(3) - 3 - log(2!) + 4·log(3) - 3 - log(4!)
//     = 6·log(3) - 6 - log(2) - log(24)
//     ≈ 6.591673732 - 6 - 0.693147181 - 3.178053830 ≈ -3.279527699
//
// 手算数据集 B (含截距, N=4, K=2):
//   X = [[1,0],[1,0],[1,1],[1,1]],  y = [1, 3, 2, 6]
//   组1 (x=0): ȳ=2, 组2 (x=1): ȳ=4
//   MLE: α = log(2), β = log(4) - log(2) = log(2)
//   β = (log2, log2) = (0.6931471805599453, 0.6931471805599453)
//   μ = [2, 2, 4, 4]
//   W = [2, 2, 4, 4]
//   X'WX = [[12, 8], [8, 8]]  (det=32)
//   Hessian V = (1/32)·[[8, -8], [-8, 12]] = [[0.25, -0.25], [-0.25, 0.375]]
//   ε = [-1, 1, -2, 2], ε² = [1, 1, 4, 4]
//   X'diag(ε²)X = [[10, 8], [8, 8]]
//   Sandwich V = A^{-1}·B·A^{-1} = [[0.125, -0.125], [-0.125, 0.25]]
//   OPG V = B^{-1} = (1/16)·[[8, -8], [-8, 10]] = [[0.5, -0.5], [-0.5, 0.625]]

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <random>
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
constexpr Real kLog2 = 0.69314718055994530941723212145817656807550013436025;
constexpr Real kLog3 = 1.0986122886681096913952452369225257046474905578227;

// 数据集 A: 无截距 N=2 K=1, β_MLE = log(3)
CrossSectionData make_data_A() {
    MatrixXD X(2, 1);
    X(0, 0) = 1.0;
    X(1, 0) = 1.0;
    VectorXD y(2);
    y(0) = 2.0;
    y(1) = 4.0;
    return make_cross_section(X, y, {"x"}, "y");
}

// 数据集 B: 含截距 N=4 K=2, β_MLE = (log2, log2)
CrossSectionData make_data_B() {
    MatrixXD X(4, 2);
    X(0, 0) = 1.0; X(0, 1) = 0.0;
    X(1, 0) = 1.0; X(1, 1) = 0.0;
    X(2, 0) = 1.0; X(2, 1) = 1.0;
    X(3, 0) = 1.0; X(3, 1) = 1.0;
    VectorXD y(4);
    y(0) = 1.0; y(1) = 3.0; y(2) = 2.0; y(3) = 6.0;
    return make_cross_section(X, y, {"const", "x"}, "y");
}

}  // namespace

// =============================================================================
// §1 Poisson link/mean 互逆 (排幻觉点 P1)
// =============================================================================
TEST(MLEPoissonTest, MeanLink_InverseRelation) {
    // 对一组 η 值, 验证 log(mean(η)) = η
    const Real eta_values[] = {-3.0, -1.0, -0.1, 0.0, 0.1, 1.0, 3.0, 5.0};
    for (Real eta : eta_values) {
        const Real mu = std::exp(eta);
        EXPECT_GT(mu, 0.0) << "Poisson mean must be positive";
        const Real eta_recovered = std::log(mu);
        EXPECT_NEAR(eta_recovered, eta, 1e-12)
            << "log(mean(eta=" << eta << ")) should be " << eta;
    }
}

// =============================================================================
// §2 Poisson IRLS weight = μ (排幻觉点 P2: w = μ for canonical log link)
// =============================================================================
TEST(MLEPoissonTest, IRLSWeight_Equals_Mean) {
    MLEEstimator mle(MLEFamily::Poisson);
    // β 使 μ 取已知值
    MatrixXD X(3, 1);
    X(0, 0) = 1.0; X(1, 0) = 1.0; X(2, 0) = 1.0;
    VectorXD beta(1);
    beta(0) = std::log(5.0);  // μ = 5 for all
    const VectorXD mu = mle.computeFittedMeans(X, beta);
    const VectorXD w = mle.computeIRLSWeights(X, beta);
    for (Size i = 0; i < 3; ++i) {
        EXPECT_NEAR(mu(i), 5.0, 1e-12);
        EXPECT_NEAR(w(i), mu(i), 1e-12) << "IRLS weight should equal mean (canonical)";
    }
}

// =============================================================================
// §3 Poisson loglik at β=0 (排幻觉点 P3: ℓ = -1 - log(y!))
// =============================================================================
TEST(MLEPoissonTest, LogLikelihood_AtBeta0) {
    // β=0 → η=0 → μ=1
    // ℓ_i = y·0 - 1 - log(y!) = -1 - lgamma(y+1)
    MLEEstimator mle(MLEFamily::Poisson);
    MatrixXD X(1, 1);
    X(0, 0) = 1.0;
    VectorXD beta(1);
    beta(0) = 0.0;

    // y=0: ℓ = -1 - log(1) = -1
    VectorXD y0(1); y0(0) = 0.0;
    EXPECT_NEAR(mle.computeLogLikelihood(X, y0, beta), -1.0, 1e-12);

    // y=1: ℓ = -1 - log(1) = -1
    VectorXD y1(1); y1(0) = 1.0;
    EXPECT_NEAR(mle.computeLogLikelihood(X, y1, beta), -1.0, 1e-12);

    // y=2: ℓ = -1 - log(2) ≈ -1.6931
    VectorXD y2(1); y2(0) = 2.0;
    EXPECT_NEAR(mle.computeLogLikelihood(X, y2, beta), -1.0 - kLog2, 1e-12);
}

// =============================================================================
// §4 Poisson MLE 无截距对称数据集 A: β = log(3)
// =============================================================================
TEST(MLEPoissonTest, SymmetricData_BetaEquals_LogMean) {
    MLEEstimator mle(MLEFamily::Poisson, CovarianceType::Hessian);
    mle.setMaxIter(100).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_A());

    ASSERT_EQ(r.coefficients.size(), 1u);
    EXPECT_NEAR(r.coefficients(0), kLog3, 1e-10);
    EXPECT_TRUE(mle.converged());
    EXPECT_LE(mle.nIterations(), 10u);
}

// =============================================================================
// §5 Poisson MLE loglik 手算 (数据集 A, ℓ ≈ -3.279527699)
// =============================================================================
TEST(MLEPoissonTest, LogLikelihood_HandComputed) {
    MLEEstimator mle(MLEFamily::Poisson, CovarianceType::Hessian);
    mle.setMaxIter(100).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_A());

    // ℓ = 6·log(3) - 6 - log(2) - log(24)
    const Real expected = 6.0 * kLog3 - 6.0 - kLog2 - std::lgamma(4.0 + 1.0);
    EXPECT_NEAR(r.log_likelihood, expected, 1e-10);
}

// =============================================================================
// §6 Poisson Hessian V = (X'WX)^{-1} (数据集 A, V = 1/6)
// =============================================================================
TEST(MLEPoissonTest, HessianVcov_HandComputed) {
    MLEEstimator mle(MLEFamily::Poisson, CovarianceType::Hessian);
    mle.setMaxIter(100).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_A());

    // X'WX = 6, V = 1/6
    ASSERT_EQ(r.vcov.rows(), 1u);
    EXPECT_NEAR(r.vcov(0, 0), 1.0 / 6.0, 1e-10);
    EXPECT_NEAR(r.std_errors(0), std::sqrt(1.0 / 6.0), 1e-10);
}

// =============================================================================
// §7 Poisson Sandwich V (数据集 A, V = 1/18)
// =============================================================================
TEST(MLEPoissonTest, SandwichVcov_HandComputed) {
    MLEEstimator mle(MLEFamily::Poisson, CovarianceType::Sandwich);
    mle.setMaxIter(100).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_A());

    // V = (1/6) · 2 · (1/6) = 1/18
    ASSERT_EQ(r.vcov.rows(), 1u);
    EXPECT_NEAR(r.vcov(0, 0), 1.0 / 18.0, 1e-10);
}

// =============================================================================
// §8 Poisson OPG V (数据集 A, V = 1/2)
// =============================================================================
TEST(MLEPoissonTest, OPGVcov_HandComputed) {
    MLEEstimator mle(MLEFamily::Poisson, CovarianceType::OPG);
    mle.setMaxIter(100).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_A());

    // V = (X'diag(ε²)X)^{-1} = 1/2
    ASSERT_EQ(r.vcov.rows(), 1u);
    EXPECT_NEAR(r.vcov(0, 0), 0.5, 1e-10);
}

// =============================================================================
// §9 Poisson MLE 含截距数据集 B: β = (log2, log2)
// =============================================================================
TEST(MLEPoissonTest, TwoGroup_BetaEquals_LogGroupMean) {
    MLEEstimator mle(MLEFamily::Poisson, CovarianceType::Hessian);
    mle.setMaxIter(100).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_B());

    ASSERT_EQ(r.coefficients.size(), 2u);
    EXPECT_NEAR(r.coefficients(0), kLog2, 1e-10);  // α = log(2)
    EXPECT_NEAR(r.coefficients(1), kLog2, 1e-10);  // β = log(2)
    EXPECT_TRUE(mle.converged());
}

// =============================================================================
// §10 Poisson Hessian V 含截距 (数据集 B, V = [[0.25,-0.25],[-0.25,0.375]])
// =============================================================================
TEST(MLEPoissonTest, HessianVcov_2D_HandComputed) {
    MLEEstimator mle(MLEFamily::Poisson, CovarianceType::Hessian);
    mle.setMaxIter(100).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_B());

    ASSERT_EQ(r.vcov.rows(), 2u);
    ASSERT_EQ(r.vcov.cols(), 2u);
    EXPECT_NEAR(r.vcov(0, 0), 0.25, 1e-10);
    EXPECT_NEAR(r.vcov(0, 1), -0.25, 1e-10);
    EXPECT_NEAR(r.vcov(1, 0), -0.25, 1e-10);
    EXPECT_NEAR(r.vcov(1, 1), 0.375, 1e-10);
}

// =============================================================================
// §11 Poisson Sandwich V 含截距 (数据集 B, V = [[0.125,-0.125],[-0.125,0.25]])
// =============================================================================
TEST(MLEPoissonTest, SandwichVcov_2D_HandComputed) {
    MLEEstimator mle(MLEFamily::Poisson, CovarianceType::Sandwich);
    mle.setMaxIter(100).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_B());

    ASSERT_EQ(r.vcov.rows(), 2u);
    EXPECT_NEAR(r.vcov(0, 0), 0.125, 1e-10);
    EXPECT_NEAR(r.vcov(0, 1), -0.125, 1e-10);
    EXPECT_NEAR(r.vcov(1, 0), -0.125, 1e-10);
    EXPECT_NEAR(r.vcov(1, 1), 0.25, 1e-10);
}

// =============================================================================
// §12 Poisson OPG V 含截距 (数据集 B, V = [[0.5,-0.5],[-0.5,0.625]])
// =============================================================================
TEST(MLEPoissonTest, OPGVcov_2D_HandComputed) {
    MLEEstimator mle(MLEFamily::Poisson, CovarianceType::OPG);
    mle.setMaxIter(100).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_B());

    ASSERT_EQ(r.vcov.rows(), 2u);
    EXPECT_NEAR(r.vcov(0, 0), 0.5, 1e-10);
    EXPECT_NEAR(r.vcov(0, 1), -0.5, 1e-10);
    EXPECT_NEAR(r.vcov(1, 0), -0.5, 1e-10);
    EXPECT_NEAR(r.vcov(1, 1), 0.625, 1e-10);
}

// =============================================================================
// §13 Poisson Score / Hessian 性质 (排幻觉点 P4: canonical link)
// =============================================================================
TEST(MLEPoissonTest, Score_Hessian_CanonicalProperties) {
    // 在 MLE 处, score = X'(y-μ) = 0
    MLEEstimator mle(MLEFamily::Poisson);
    MatrixXD X(2, 1);
    X(0, 0) = 1.0; X(1, 0) = 1.0;
    VectorXD y(2);
    y(0) = 2.0; y(1) = 4.0;

    VectorXD beta_mle(1);
    beta_mle(0) = kLog3;
    const VectorXD mu = mle.computeFittedMeans(X, beta_mle);

    // score = X'(y-μ) = 1·(2-3) + 1·(4-3) = 0
    const VectorXD score = VectorXD(X.eigen().transpose() * (y.eigen() - mu.eigen()));
    EXPECT_NEAR(score(0), 0.0, 1e-12);

    // 在 MLE 处, W = μ
    const VectorXD W = mle.computeIRLSWeights(X, beta_mle);
    EXPECT_NEAR(W(0), 3.0, 1e-12);
    EXPECT_NEAR(W(1), 3.0, 1e-12);
}

// =============================================================================
// §14 Poisson MLE 大样本: 三种协方差渐近等价
// =============================================================================
TEST(MLEPoissonTest, ThreeVcovTypes_ConvergeLargeSample) {
    // 构造大样本 Poisson 数据, 验证 Hessian/Sandwich/OPG 渐近等价
    const Size N = 2000;
    MatrixXD X(N, 2);
    VectorXD y(N);

    // y ~ Poisson(exp(0.5 + 0.3·x)), x ~ N(0,1)
    std::mt19937 rng(42);
    std::normal_distribution<Real> norm(0.0, 1.0);
    for (Size i = 0; i < N; ++i) {
        X(i, 0) = 1.0;
        X(i, 1) = norm(rng);
        const Real eta = 0.5 + 0.3 * X(i, 1);
        const Real mu = std::exp(eta);
        // Poisson 采样 via inverse CDF (Knuth)
        std::poisson_distribution<int> pois(static_cast<double>(mu));
        y(i) = static_cast<Real>(pois(rng));
    }
    auto data = make_cross_section(X, y, {"const", "x"}, "y");

    MLEEstimator mle_h(MLEFamily::Poisson, CovarianceType::Hessian);
    mle_h.setMaxIter(100).setTolerance(1e-10);
    EstimationResult rh = mle_h.estimate(data);

    MLEEstimator mle_s(MLEFamily::Poisson, CovarianceType::Sandwich);
    mle_s.setMaxIter(100).setTolerance(1e-10);
    EstimationResult rs = mle_s.estimate(data);

    MLEEstimator mle_o(MLEFamily::Poisson, CovarianceType::OPG);
    mle_o.setMaxIter(100).setTolerance(1e-10);
    EstimationResult ro = mle_o.estimate(data);

    ASSERT_TRUE(mle_h.converged());
    ASSERT_TRUE(mle_s.converged());
    ASSERT_TRUE(mle_o.converged());

    // 三种 SE 应渐近接近 (相对差异 < 15%)
    for (Size i = 0; i < 2; ++i) {
        const Real se_h = rh.std_errors(i);
        const Real se_s = rs.std_errors(i);
        const Real se_o = ro.std_errors(i);
        EXPECT_GT(se_h, 0.0);
        EXPECT_GT(se_s, 0.0);
        EXPECT_GT(se_o, 0.0);
        const Real rel_diff_hs = std::fabs(se_h - se_s) / se_h;
        const Real rel_diff_ho = std::fabs(se_h - se_o) / se_h;
        EXPECT_LT(rel_diff_hs, 0.15)
            << "Hessian vs Sandwich SE[" << i << "] diff " << rel_diff_hs;
        EXPECT_LT(rel_diff_ho, 0.20)
            << "Hessian vs OPG SE[" << i << "] diff " << rel_diff_ho;
    }
}

// =============================================================================
// §15 Poisson MLE 不支持错误协方差类型
// =============================================================================
TEST(MLEPoissonTest, ThrowsOnUnsupportedCovarianceType) {
    MLEEstimator mle(MLEFamily::Poisson, CovarianceType::Classical);
    EXPECT_THROW(mle.estimate(make_data_A()), std::invalid_argument);
}
