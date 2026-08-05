// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.1 任务 2.1-2.2 - Logistic MLE 测试
// 验证方法: link/mean 互逆 + Newton-Raphson 一步手算 + 对称数据 MLE 解析解 (容差 1e-10)
//
// 排幻觉点:
//   L1: Logistic mean(link(x))=x, link(mean(x))=x (logit/logistic 互逆)
//   L2: IRLS weight w=μ(1-μ) ∈ [0, 0.25] (μ(1-μ) ≤ 0.25, μ=0.5 时取最大)
//   L3: ℓ(y, η=0) = y·0 - log(1+e^0) = -log(2) (β=0 时所有观测对数似然相等)
//   E7: R glm 用 IRLS, C++ 用 Newton-Raphson + 解析 Hessian, canonical link 等价
//   E8: GLM bread = (X'WX)^{-1}, meat = X' diag(ε²) X (R sandwich::sandwich 实现)
//
// 手算数据集 A (对称, N=2, K=1, 无截距):
//   x = [1, 1],  y = [0, 1]
//   β₀=0: η=[0,0], μ=[0.5, 0.5]
//   score = 1·(0-0.5) + 1·(1-0.5) = 0  →  Δ=0, β₁=0 (一步收敛)
//   MLE: β = 0, μ = [0.5, 0.5]
//   ℓ = log(0.5) + log(0.5) = -2·log(2) ≈ -1.3862943611
//   W = [0.25, 0.25], X'WX = 0.25+0.25 = 0.5
//   Hessian V = (X'WX)^{-1} = 2
//   ε = [-0.5, 0.5], X'diag(ε²)X = 0.25+0.25 = 0.5
//   Sandwich V = (X'WX)^{-1} · X'diag(ε²)X · (X'WX)^{-1} = 2·0.5·2 = 2
//   (注: Hessian = Sandwich 因对称数据 ε²=μ(1-μ)=0.25)
//
// 手算数据集 B (N=4, K=1, 无截距, 均值 y=0.5):
//   x = [1, 1, 1, 1],  y = [1, 1, 0, 0]
//   MLE: μ = ȳ = 0.5 → β = 0 (因 x 全 1, μ=σ(logistic(β·1)))
//   score = Σ(y-μ) = 2-4·0.5 = 0 → β=0
//   W = 0.25 for all, X'WX = 4·0.25 = 1
//   Hessian V = 1
//   ε = [0.5, 0.5, -0.5, -0.5], X'diag(ε²)X = 4·0.25 = 1
//   Sandwich V = 1·1·1 = 1
//   ℓ = 4·(-log(2)) = -4·log(2) ≈ -2.7725887222
//
// 手算数据集 C (非对称, N=3, K=1, Newton-Raphson 一步):
//   x = [1, 2, 3],  y = [1, 0, 1]
//   β₀=0: μ=[0.5, 0.5, 0.5]
//   score = 1·0.5 + 2·(-0.5) + 3·0.5 = 0.5 - 1 + 1.5 = 1.0
//   X'WX = 0.25·(1+4+9) = 0.25·14 = 3.5
//   Δ = score / (X'WX) = 1/3.5 = 0.2857142857...
//   β₁ = 0 + 0.2857142857 ≈ 0.2857142857

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

// 数据集 A: 对称 N=2 K=1, β_MLE=0
CrossSectionData make_data_A() {
    MatrixXD X(2, 1);
    X(0, 0) = 1.0;
    X(1, 0) = 1.0;
    VectorXD y(2);
    y(0) = 0.0;
    y(1) = 1.0;
    return make_cross_section(X, y, {"x"}, "y");
}

// 数据集 B: 同质 N=4 K=1, β_MLE=0
CrossSectionData make_data_B() {
    MatrixXD X(4, 1);
    for (Size i = 0; i < 4; ++i) X(i, 0) = 1.0;
    VectorXD y(4);
    y(0) = 1.0; y(1) = 1.0; y(2) = 0.0; y(3) = 0.0;
    return make_cross_section(X, y, {"x"}, "y");
}

// 数据集 C: 非对称 N=3 K=1, Newton-Raphson 一步验证
CrossSectionData make_data_C() {
    MatrixXD X(3, 1);
    X(0, 0) = 1.0; X(1, 0) = 2.0; X(2, 0) = 3.0;
    VectorXD y(3);
    y(0) = 1.0; y(1) = 0.0; y(2) = 1.0;
    return make_cross_section(X, y, {"x"}, "y");
}

// 数据集 D: 完全分离 N=4 K=2 (含截距), MLE 不存在 (β→+∞)
//   X = [[1,1],[1,2],[1,3],[1,4]], y = [0,0,1,1]
//   决策边界 α+β·x = 0 位于 x=2.5, β→∞ 可完美分离
//   注: 无截距时不构成完全分离 (强制 μ(0)=0.5), 必须含截距
CrossSectionData make_data_D_separated() {
    MatrixXD X(4, 2);
    X(0, 0) = 1.0; X(0, 1) = 1.0;
    X(1, 0) = 1.0; X(1, 1) = 2.0;
    X(2, 0) = 1.0; X(2, 1) = 3.0;
    X(3, 0) = 1.0; X(3, 1) = 4.0;
    VectorXD y(4);
    y(0) = 0.0; y(1) = 0.0; y(2) = 1.0; y(3) = 1.0;
    return make_cross_section(X, y, {"const", "x"}, "y");
}

}  // namespace

// =============================================================================
// §1 Logistic link/mean 互逆 (排幻觉点 L1)
// =============================================================================
TEST(MLELogisticTest, MeanLink_InverseRelation) {
    // 对一组 η 值, 验证 link(mean(η)) = η
    const Real eta_values[] = {-5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0};
    for (Real eta : eta_values) {
        // mean(η) = 1/(1+e^{-η})
        const Real mu = 1.0 / (1.0 + std::exp(-eta));
        // link(μ) = log(μ/(1-μ))
        const Real eta_recovered = std::log(mu / (1.0 - mu));
        EXPECT_NEAR(eta_recovered, eta, 1e-10)
            << "link(mean(eta=" << eta << ")) should be " << eta;
    }
}

// =============================================================================
// §2 Logistic IRLS weight bounds (排幻觉点 L2: w=μ(1-μ) ≤ 0.25)
// =============================================================================
TEST(MLELogisticTest, IRLSWeight_BoundedByQuarter) {
    // 对一组 μ 值, 验证 w = μ(1-μ) ≤ 0.25
    const Real mu_values[] = {0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99};
    for (Real mu : mu_values) {
        const Real w = mu * (1.0 - mu);
        EXPECT_LE(w, 0.25 + 1e-15) << "w=" << w << " for mu=" << mu;
        EXPECT_GT(w, 0.0) << "w should be positive";
    }
    // μ=0.5 时 w 取最大 0.25
    EXPECT_NEAR(0.5 * 0.5, 0.25, 1e-15);
}

// =============================================================================
// §3 Logistic loglik at β=0 (排幻觉点 L3: ℓ=-log(2))
// =============================================================================
TEST(MLELogisticTest, LogLikelihood_AtBeta0) {
    // 单观测 ℓ(y, η=0) = -log(2) (因 ℓ = y·0 - log(1+e^0) = -log(2))
    MLEEstimator mle(MLEFamily::Logistic);
    MatrixXD X(1, 1);
    X(0, 0) = 1.0;
    VectorXD y(1);
    y(0) = 1.0;
    VectorXD beta(1);
    beta(0) = 0.0;

    const Real ll = mle.computeLogLikelihood(X, y, beta);
    EXPECT_NEAR(ll, -kLog2, 1e-15);

    // y=0 时同样 ℓ=-log(2)
    y(0) = 0.0;
    const Real ll0 = mle.computeLogLikelihood(X, y, beta);
    EXPECT_NEAR(ll0, -kLog2, 1e-15);
}

// =============================================================================
// §4 Logistic MLE 对称数据集 A: β=0 (手算解析解)
// =============================================================================
TEST(MLELogisticTest, SymmetricData_BetaIsZero) {
    MLEEstimator mle(MLEFamily::Logistic, CovarianceType::Hessian);
    mle.setMaxIter(50).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_A());

    // β = 0 (对称数据 MLE)
    ASSERT_EQ(r.coefficients.size(), 1u);
    EXPECT_NEAR(r.coefficients(0), 0.0, 1e-10);

    // loglik = -2·log(2)
    EXPECT_NEAR(r.log_likelihood, -2.0 * kLog2, 1e-10);

    // 收敛
    EXPECT_TRUE(mle.converged());
    EXPECT_LE(mle.nIterations(), 5u);  // 对称数据应一步收敛
}

// =============================================================================
// §5 Logistic MLE Hessian V = (X'WX)^{-1} (手算数据集 A, V=2)
// =============================================================================
TEST(MLELogisticTest, HessianVcov_HandComputed) {
    MLEEstimator mle(MLEFamily::Logistic, CovarianceType::Hessian);
    mle.setMaxIter(50).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_A());

    // 数据集 A: β=0, μ=[0.5,0.5], W=[0.25,0.25]
    // X'WX = 0.25·1 + 0.25·1 = 0.5
    // Hessian V = (X'WX)^{-1} = 2.0
    ASSERT_EQ(r.vcov.rows(), 1u);
    ASSERT_EQ(r.vcov.cols(), 1u);
    EXPECT_NEAR(r.vcov(0, 0), 2.0, 1e-10);

    // SE = sqrt(2)
    ASSERT_EQ(r.std_errors.size(), 1u);
    EXPECT_NEAR(r.std_errors(0), std::sqrt(2.0), 1e-10);
}

// =============================================================================
// §6 Logistic MLE Sandwich V (手算数据集 A, V=2)
// =============================================================================
TEST(MLELogisticTest, SandwichVcov_HandComputed) {
    MLEEstimator mle(MLEFamily::Logistic, CovarianceType::Sandwich);
    mle.setMaxIter(50).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_A());

    // 数据集 A: ε=[-0.5, 0.5], ε²=[0.25, 0.25]
    // X'diag(ε²)X = 0.25·1 + 0.25·1 = 0.5
    // (X'WX)^{-1} = 2
    // Sandwich V = 2·0.5·2 = 2
    ASSERT_EQ(r.vcov.rows(), 1u);
    EXPECT_NEAR(r.vcov(0, 0), 2.0, 1e-10);
}

// =============================================================================
// §7 Logistic MLE 同质数据集 B: β=0, ℓ=-4·log(2), V=1
// =============================================================================
TEST(MLELogisticTest, HomogeneousData_BetaZero) {
    MLEEstimator mle(MLEFamily::Logistic, CovarianceType::Hessian);
    mle.setMaxIter(50).setTolerance(1e-12);
    EstimationResult r = mle.estimate(make_data_B());

    // β = 0 (因 y 均值 = 0.5)
    ASSERT_EQ(r.coefficients.size(), 1u);
    EXPECT_NEAR(r.coefficients(0), 0.0, 1e-10);

    // ℓ = -4·log(2)
    EXPECT_NEAR(r.log_likelihood, -4.0 * kLog2, 1e-10);

    // Hessian V = 1 (X'WX = 4·0.25 = 1, V = 1)
    EXPECT_NEAR(r.vcov(0, 0), 1.0, 1e-10);
}

// =============================================================================
// §8 Logistic Newton-Raphson 一步手算 (数据集 C, β₀=0 → β₁=1/3.5)
// =============================================================================
TEST(MLELogisticTest, NewtonRaphson_SingleStep_HandComputed) {
    // 构造手动调用 Newton-Raphson 一步, 验证 β₁ = 1/3.5
    // mle.hpp 中 newtonRaphson 是 private, 通过 estimate 调用
    // 这里改用 IRLS weights / score 的间接验证
    MLEEstimator mle(MLEFamily::Logistic);
    MatrixXD X(3, 1);
    X(0, 0) = 1.0; X(1, 0) = 2.0; X(2, 0) = 3.0;
    VectorXD y(3);
    y(0) = 1.0; y(1) = 0.0; y(2) = 1.0;

    // β₀ = 0
    VectorXD beta0(1);
    beta0(0) = 0.0;

    // 在 β₀ = 0 处计算 μ, W
    const VectorXD mu = mle.computeFittedMeans(X, beta0);
    const VectorXD W = mle.computeIRLSWeights(X, beta0);

    // 验证 μ = [0.5, 0.5, 0.5]
    for (Size i = 0; i < 3; ++i) {
        EXPECT_NEAR(mu(i), 0.5, 1e-15);
        EXPECT_NEAR(W(i), 0.25, 1e-15);
    }

    // score = X'(y - μ)
    const VectorXD score = VectorXD(X.eigen().transpose() * (y.eigen() - mu.eigen()));
    EXPECT_NEAR(score(0), 1.0, 1e-15);  // 1·0.5 + 2·(-0.5) + 3·0.5 = 1.0

    // X'WX = 0.25·(1+4+9) = 3.5
    const Real XtWX = (W.eigen().array() * X.eigen().array().square()).sum();
    EXPECT_NEAR(XtWX, 3.5, 1e-15);

    // Δ = score / (X'WX) = 1/3.5 = 0.2857142857...
    const Real delta = score(0) / XtWX;
    EXPECT_NEAR(delta, 1.0 / 3.5, 1e-15);
    EXPECT_NEAR(delta, 0.2857142857142857, 1e-15);

    // β₁ = 0 + 0.2857142857
    const Real beta1 = beta0(0) + delta;
    EXPECT_NEAR(beta1, 0.2857142857142857, 1e-15);
}

// =============================================================================
// §9 Logistic MLE 完全分离 (排幻觉: 应抛异常 / 不收敛 / β 爆炸)
// =============================================================================
TEST(MLELogisticTest, CompleteSeparation_ThrowsOrDiverges) {
    MLEEstimator mle(MLEFamily::Logistic, CovarianceType::Hessian);
    mle.setMaxIter(100).setTolerance(1e-10);

    // 完全分离数据 y=[0,0,1,1], x=[1,2,3,4]
    // MLE β → +∞, 数值实现的三种合理行为:
    //   (a) 抛 runtime_error (X'WX 奇异, W=μ(1-μ)→0)
    //   (b) converged()=false (β 持续增大未达容差)
    //   (c) converged()=true 但 |β| > 10 (数值上"收敛"到发散解)
    // R glm() 同款行为: 返回大 β + 警告 "fitted probabilities numerically 0 or 1"
    bool throws = false;
    bool not_converged = false;
    bool beta_diverges = false;
    try {
        EstimationResult r = mle.estimate(make_data_D_separated());
        // 含截距完全分离: β_slope 或 |α| 应发散到较大值
        if (!mle.converged()) {
            not_converged = true;
        } else if (std::fabs(r.coefficients(0)) > 10.0
                   || std::fabs(r.coefficients(1)) > 10.0) {
            beta_diverges = true;
        }
    } catch (const std::runtime_error&) {
        throws = true;
    } catch (const std::exception&) {
        throws = true;
    }
    EXPECT_TRUE(throws || not_converged || beta_diverges)
        << "Complete separation should: throw, not converge, or return divergent beta";
}

// =============================================================================
// §10 Logistic Bernoulli 别名等价 (Bernoulli == Logistic)
// =============================================================================
TEST(MLELogisticTest, Bernoulli_AliasEquivalent) {
    MLEEstimator mle_logistic(MLEFamily::Logistic, CovarianceType::Hessian);
    MLEEstimator mle_bernoulli(MLEFamily::Bernoulli, CovarianceType::Hessian);

    EstimationResult r1 = mle_logistic.estimate(make_data_A());
    EstimationResult r2 = mle_bernoulli.estimate(make_data_A());

    EXPECT_NEAR(r1.coefficients(0), r2.coefficients(0), 1e-15);
    EXPECT_NEAR(r1.log_likelihood, r2.log_likelihood, 1e-15);
    EXPECT_NEAR(r1.vcov(0, 0), r2.vcov(0, 0), 1e-15);
}

// =============================================================================
// §11 Logistic 收敛性: 多起始值都收敛到同一解
// =============================================================================
TEST(MLELogisticTest, ConvergenceFromMultipleStarts) {
    // 数据集 C (非对称, MLE 存在), 不同起始值应收敛到同一 β
    MLEEstimator mle1(MLEFamily::Logistic, CovarianceType::Hessian);
    mle1.setMaxIter(200).setTolerance(1e-12);

    VectorXD start1(1); start1(0) = 0.0;
    VectorXD start2(1); start2(0) = 0.5;
    VectorXD start3(1); start3(0) = -0.5;

    mle1.setStartValues(start1);
    EstimationResult r1 = mle1.estimate(make_data_C());

    mle1.setStartValues(start2);
    EstimationResult r2 = mle1.estimate(make_data_C());

    mle1.setStartValues(start3);
    EstimationResult r3 = mle1.estimate(make_data_C());

    EXPECT_TRUE(mle1.converged());
    EXPECT_NEAR(r1.coefficients(0), r2.coefficients(0), 1e-8);
    EXPECT_NEAR(r1.coefficients(0), r3.coefficients(0), 1e-8);
    EXPECT_NEAR(r1.log_likelihood, r2.log_likelihood, 1e-8);
}

// =============================================================================
// §12 Logistic MLE Hessian vs Sandwich vs OPG 渐近等价性 (大样本)
// =============================================================================
TEST(MLELogisticTest, ThreeVcovTypes_ConvergeLargeSample) {
    // 构造大样本 (N=2000) 数据, 验证三种协方差矩阵对角元素渐近接近
    // 注: 小样本下三者有差异, 大样本下渐近等价 (White 1982)
    // 因 logistic MLE 无 σ² 修正问题, 三者直接对比 (不像 Gaussian)
    const Size N = 2000;
    MatrixXD X(N, 2);
    VectorXD y(N);

    // 构造: y = 1 if (1 + 0.5·x + N(0,1) > 0), x ~ N(0,1)
    std::mt19937 rng(42);
    std::normal_distribution<Real> norm(0.0, 1.0);
    for (Size i = 0; i < N; ++i) {
        X(i, 0) = 1.0;
        X(i, 1) = norm(rng);
        const Real eta_true = 1.0 + 0.5 * X(i, 1) + norm(rng);
        const Real mu_true = 1.0 / (1.0 + std::exp(-eta_true));
        // Bernoulli 采样
        const Real u = static_cast<Real>(rng()) / static_cast<Real>(rng.max());
        y(i) = (u < mu_true) ? 1.0 : 0.0;
    }
    auto data = make_cross_section(X, y, {"const", "x"}, "y");

    MLEEstimator mle_h(MLEFamily::Logistic, CovarianceType::Hessian);
    mle_h.setMaxIter(100).setTolerance(1e-10);
    EstimationResult rh = mle_h.estimate(data);

    MLEEstimator mle_s(MLEFamily::Logistic, CovarianceType::Sandwich);
    mle_s.setMaxIter(100).setTolerance(1e-10);
    EstimationResult rs = mle_s.estimate(data);

    MLEEstimator mle_o(MLEFamily::Logistic, CovarianceType::OPG);
    mle_o.setMaxIter(100).setTolerance(1e-10);
    EstimationResult ro = mle_o.estimate(data);

    ASSERT_TRUE(mle_h.converged());
    ASSERT_TRUE(mle_s.converged());
    ASSERT_TRUE(mle_o.converged());

    // 三种 SE 应渐近接近 (相对差异 < 10%)
    for (Size i = 0; i < 2; ++i) {
        const Real se_h = rh.std_errors(i);
        const Real se_s = rs.std_errors(i);
        const Real se_o = ro.std_errors(i);
        EXPECT_GT(se_h, 0.0);
        EXPECT_GT(se_s, 0.0);
        EXPECT_GT(se_o, 0.0);
        // 渐近等价: |se_h - se_s| / se_h < 0.1
        const Real rel_diff_hs = std::fabs(se_h - se_s) / se_h;
        const Real rel_diff_ho = std::fabs(se_h - se_o) / se_h;
        EXPECT_LT(rel_diff_hs, 0.10)
            << "Hessian vs Sandwich SE[" << i << "] diff " << rel_diff_hs;
        EXPECT_LT(rel_diff_ho, 0.15)
            << "Hessian vs OPG SE[" << i << "] diff " << rel_diff_ho;
    }
}
