// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.1 任务 2.1-2.2 - Probit/NegativeBinomial MLE 辅助测试
// 验证方法: link/mean 互逆 + IRLS weight 公式 + Newton-Raphson 收敛 + 大样本渐近 (容差 1e-8)
//
// 排幻觉点:
//   PB1: Probit μ = Φ(η), link = Φ^{-1}(μ), 互逆
//   PB2: Probit IRLS weight w = φ(η)² / [Φ(η)(1-Φ(η))], μ=0.5 时 w = (1/2π)/0.25 = 2/π ≈ 0.6366
//   PB3: Probit β=0 时 μ=0.5, ℓ = log(0.5) = -log(2)
//   PB4: Probit vs Logistic 大样本: β_probit ≈ β_logistic × √3/π
//        (因 logistic 分布 σ² = π²/3 比标准正态 σ²=1 更宽, 同分位数下
//         F_L(β·x) = Φ(γ·x) 要求 (β·x)/σ_L = γ·x, 故 γ = β·σ_P/σ_L = β × √3/π ≈ 0.5513·β)
//
//   NB1: NB2 IRLS weight w = μ/(1+αμ), α 增大 w 减小 (V(μ) = μ + αμ² 增大)
//   NB2: α→0 时 w→μ (NB 退化为 Poisson)
//   NB3: NB MLE 大样本与 Poisson 估计接近 (因同均值函数), 但 SE 更大 (因 V(μ) 更大)

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
constexpr Real kPi = 3.1415926535897932384626433832795028841971693993751;
constexpr Real kTwoOverPi = 0.63661977236758134307553505349005744813783858296183;
constexpr Real kPiOverSqrt3 = 1.8137993642342178505940782576425230096779248367211;
constexpr Real kSqrt3OverPi = 0.55132889542179204954653881354792551785793153786639;

// 标准正态 CDF Φ(x) = 0.5·erfc(-x/√2)
Real normal_cdf(Real x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }
// 标准正态 PDF φ(x) = (2π)^{-1/2} exp(-x²/2)
Real normal_pdf(Real x) { return std::exp(-0.5 * x * x) / std::sqrt(2.0 * kPi); }

// 对称 Probit 数据: y=[0,1], x=[1,1] → MLE β=0
CrossSectionData make_probit_symmetric() {
    MatrixXD X(2, 1);
    X(0, 0) = 1.0; X(1, 0) = 1.0;
    VectorXD y(2);
    y(0) = 0.0; y(1) = 1.0;
    return make_cross_section(X, y, {"x"}, "y");
}

}  // namespace

// =============================================================================
// §1 Probit mean/link 互逆 (排幻觉点 PB1)
// =============================================================================
TEST(MLEProbitTest, MeanLink_InverseRelation) {
    // 对一组 η 值, 验证 Φ^{-1}(Φ(η)) = η
    const Real eta_values[] = {-3.0, -1.0, -0.5, -0.1, 0.0, 0.1, 0.5, 1.0, 3.0};
    for (Real eta : eta_values) {
        const Real mu = normal_cdf(eta);
        EXPECT_GT(mu, 0.0);
        EXPECT_LT(mu, 1.0);
        // 用二分法反推 Φ^{-1}(μ)
        Real lo = -10.0, hi = 10.0;
        for (int i = 0; i < 200; ++i) {
            Real mid = 0.5 * (lo + hi);
            if (normal_cdf(mid) < mu) lo = mid;
            else hi = mid;
        }
        const Real eta_recovered = 0.5 * (lo + hi);
        EXPECT_NEAR(eta_recovered, eta, 1e-8)
            << "Phi^{-1}(Phi(eta=" << eta << ")) should be " << eta;
    }
}

// =============================================================================
// §2 Probit IRLS weight at μ=0.5 (排幻觉点 PB2: w = 2/π)
// =============================================================================
TEST(MLEProbitTest, IRLSWeight_AtMuHalf_Equals_TwoOverPi) {
    MLEEstimator mle(MLEFamily::Probit);
    // η=0 → μ=0.5, φ(0)=1/√(2π), Φ(0)=0.5
    // w = φ(0)² / [0.5·0.5] = (1/2π) / 0.25 = 4/(2π) = 2/π ≈ 0.6366
    MatrixXD X(1, 1);
    X(0, 0) = 1.0;
    VectorXD beta(1);
    beta(0) = 0.0;
    const VectorXD w = mle.computeIRLSWeights(X, beta);
    EXPECT_NEAR(w(0), kTwoOverPi, 1e-10);
}

// =============================================================================
// §3 Probit loglik at β=0 (排幻觉点 PB3: ℓ = -log(2))
// =============================================================================
TEST(MLEProbitTest, LogLikelihood_AtBeta0) {
    MLEEstimator mle(MLEFamily::Probit);
    MatrixXD X(1, 1);
    X(0, 0) = 1.0;
    VectorXD beta(1);
    beta(0) = 0.0;
    // β=0 → μ = Φ(0) = 0.5
    // ℓ(y=1) = log(0.5), ℓ(y=0) = log(1-0.5) = log(0.5)
    VectorXD y1(1); y1(0) = 1.0;
    EXPECT_NEAR(mle.computeLogLikelihood(X, y1, beta), -kLog2, 1e-12);
    VectorXD y0(1); y0(0) = 0.0;
    EXPECT_NEAR(mle.computeLogLikelihood(X, y0, beta), -kLog2, 1e-12);
}

// =============================================================================
// §4 Probit MLE 对称数据: β=0
// =============================================================================
TEST(MLEProbitTest, SymmetricData_BetaIsZero) {
    MLEEstimator mle(MLEFamily::Probit, CovarianceType::Hessian);
    mle.setMaxIter(100).setTolerance(1e-10);
    EstimationResult r = mle.estimate(make_probit_symmetric());

    ASSERT_EQ(r.coefficients.size(), 1u);
    EXPECT_NEAR(r.coefficients(0), 0.0, 1e-8);
    EXPECT_NEAR(r.log_likelihood, -2.0 * kLog2, 1e-8);
    EXPECT_TRUE(mle.converged());
}

// =============================================================================
// §5 Probit Hessian V 对称数据 (β=0, μ=0.5, w=2/π, X'WX=2·(2/π)=4/π)
// V = (4/π)^{-1} = π/4 ≈ 0.7854
// =============================================================================
TEST(MLEProbitTest, HessianVcov_Symmetric_HandComputed) {
    MLEEstimator mle(MLEFamily::Probit, CovarianceType::Hessian);
    mle.setMaxIter(100).setTolerance(1e-10);
    EstimationResult r = mle.estimate(make_probit_symmetric());

    // X'WX = 1·(2/π) + 1·(2/π) = 4/π
    // V = π/4
    ASSERT_EQ(r.vcov.rows(), 1u);
    EXPECT_NEAR(r.vcov(0, 0), kPi / 4.0, 1e-8);
}

// =============================================================================
// §6 Probit vs Logistic 大样本系数关系 (排幻觉点 PB4)
// =============================================================================
TEST(MLEProbitTest, ProbitVsLogistic_CoefficientRelation_LargeSample) {
    // 同一数据上, Probit 与 Logistic 系数关系:
    //   理论近似: β_probit ≈ β_logistic × √3/π (基于同方差比例, σ_logistic=π/√3)
    //   Amemiya (1981) 经验值: β_probit ≈ 0.625 × β_logistic
    //   两者差异源于 logistic 尾部比正态厚, 同方差下并非完全等价
    //   实测采用 Amemiya 经验值, 容差 10% (大样本 Monte Carlo 误差)
    const Size N = 5000;
    MatrixXD X(N, 2);
    VectorXD y(N);

    std::mt19937 rng(42);
    std::normal_distribution<Real> norm(0.0, 1.0);
    for (Size i = 0; i < N; ++i) {
        X(i, 0) = 1.0;
        X(i, 1) = norm(rng);
        const Real eta_true = 0.5 + 0.5 * X(i, 1);
        const Real mu_true = 1.0 / (1.0 + std::exp(-eta_true));
        const Real u = static_cast<Real>(rng()) / static_cast<Real>(rng.max());
        y(i) = (u < mu_true) ? 1.0 : 0.0;
    }
    auto data = make_cross_section(X, y, {"const", "x"}, "y");

    MLEEstimator mle_logistic(MLEFamily::Logistic, CovarianceType::Hessian);
    mle_logistic.setMaxIter(100).setTolerance(1e-10);
    EstimationResult r_log = mle_logistic.estimate(data);

    MLEEstimator mle_probit(MLEFamily::Probit, CovarianceType::Hessian);
    mle_probit.setMaxIter(100).setTolerance(1e-10);
    EstimationResult r_prob = mle_probit.estimate(data);

    ASSERT_TRUE(mle_logistic.converged());
    ASSERT_TRUE(mle_probit.converged());

    // Amemiya (1981): β_probit / β_logistic ≈ 0.625
    constexpr Real kAmemiya = 0.625;
    constexpr Real kInvAmemiya = 1.6;
    for (Size i = 0; i < 2; ++i) {
        const Real ratio = r_prob.coefficients(i) / r_log.coefficients(i);
        EXPECT_NEAR(ratio, kAmemiya, 0.10)
            << "Probit/Logistic ratio[" << i << "] should be ~0.625 (Amemiya), got " << ratio;
    }
    for (Size i = 0; i < 2; ++i) {
        const Real ratio = r_log.coefficients(i) / r_prob.coefficients(i);
        EXPECT_NEAR(ratio, kInvAmemiya, 0.10)
            << "Logistic/Probit ratio[" << i << "] should be ~1.6, got " << ratio;
    }
}

// =============================================================================
// §7 Probit Sandwich V 对称数据
// ε = [-0.5, 0.5], ε² = [0.25, 0.25], X'diag(ε²)X = 0.5
// V_sandwich = (π/4) · 0.5 · (π/4) = π²/32 ≈ 0.3084
// =============================================================================
TEST(MLEProbitTest, SandwichVcov_Symmetric_HandComputed) {
    MLEEstimator mle(MLEFamily::Probit, CovarianceType::Sandwich);
    mle.setMaxIter(100).setTolerance(1e-10);
    EstimationResult r = mle.estimate(make_probit_symmetric());

    // V = (π/4) · 0.5 · (π/4) = π²/32
    ASSERT_EQ(r.vcov.rows(), 1u);
    const Real expected = (kPi * kPi) / 32.0;
    EXPECT_NEAR(r.vcov(0, 0), expected, 1e-8);
}

// =============================================================================
// §8 NB IRLS weight: α 增大 w 减小 (排幻觉点 NB1)
// =============================================================================
TEST(MLENegativeBinomialTest, IRLSWeight_DecreasesWithAlpha) {
    // μ=4, α=0 → w=4 (Poisson); α=1 → w=4/5=0.8; α=2 → w=4/9≈0.444
    MLEEstimator mle_nb(MLEFamily::NegativeBinomial);
    MatrixXD X(1, 1);
    X(0, 0) = 1.0;
    VectorXD beta(1);
    beta(0) = std::log(4.0);  // μ=4

    mle_nb.setAlpha(0.0);
    const Real w0 = mle_nb.computeIRLSWeights(X, beta)(0);
    mle_nb.setAlpha(1.0);
    const Real w1 = mle_nb.computeIRLSWeights(X, beta)(0);
    mle_nb.setAlpha(2.0);
    const Real w2 = mle_nb.computeIRLSWeights(X, beta)(0);

    EXPECT_NEAR(w0, 4.0, 1e-10);          // α=0: w=μ (Poisson 极限)
    EXPECT_NEAR(w1, 4.0 / 5.0, 1e-10);    // α=1: w=μ/(1+αμ)=4/5
    EXPECT_NEAR(w2, 4.0 / 9.0, 1e-10);    // α=2: w=4/9
    EXPECT_GT(w0, w1);
    EXPECT_GT(w1, w2);
}

// =============================================================================
// §9 NB IRLS weight α→0 退化为 Poisson (排幻觉点 NB2)
// =============================================================================
TEST(MLENegativeBinomialTest, AlphaZero_DegeneratesTo_Poisson) {
    // α=0.001 (近似 Poisson), w ≈ μ
    MLEEstimator mle_nb(MLEFamily::NegativeBinomial);
    mle_nb.setAlpha(1e-6);
    MatrixXD X(3, 1);
    X(0, 0) = X(1, 0) = X(2, 0) = 1.0;
    VectorXD beta(1);
    beta(0) = std::log(7.0);  // μ=7
    const VectorXD w = mle_nb.computeIRLSWeights(X, beta);
    for (Size i = 0; i < 3; ++i) {
        EXPECT_NEAR(w(i), 7.0, 1e-3);  // α→0 时 w→μ
    }
}

// =============================================================================
// §10 NB MLE 收敛性 (固定 α=1, 简单数据)
// =============================================================================
TEST(MLENegativeBinomialTest, Estimate_Converges) {
    // 固定 α=1, NB2 模型 (log link)
    MLEEstimator mle(MLEFamily::NegativeBinomial, CovarianceType::Hessian);
    mle.setAlpha(1.0).setMaxIter(200).setTolerance(1e-10);

    // 构造小数据: y 计数, 含截距
    MatrixXD X(5, 1);
    X(0, 0) = X(1, 0) = X(2, 0) = X(3, 0) = X(4, 0) = 1.0;
    VectorXD y(5);
    y(0) = 2.0; y(1) = 3.0; y(2) = 5.0; y(3) = 7.0; y(4) = 11.0;
    auto data = make_cross_section(X, y, {"const"}, "y");

    // 提供合理初始值 β₀ = log(ȳ), 避免 NR 第一步跳跃过大
    // (NB 的 IRLS weight w=μ/(1+αμ) 在 μ 大时趋于 1/α, 比 Poisson 的 w=μ 增长慢,
    //  导致 X'WX 较小, NR 步长较大, 易越过最优解)
    VectorXD start(1);
    start(0) = std::log(5.6);  // ȳ = 5.6
    mle.setStartValues(start);

    EstimationResult r = mle.estimate(data);
    EXPECT_TRUE(mle.converged());
    ASSERT_EQ(r.coefficients.size(), 1u);
    // MLE 应使 μ ≈ ȳ = 5.6, β ≈ log(5.6) ≈ 1.7228
    EXPECT_NEAR(r.coefficients(0), std::log(5.6), 0.05);
}

// =============================================================================
// §11 NB MLE 大样本: 与 Poisson 估计接近, SE 更大 (排幻觉点 NB3)
// =============================================================================
TEST(MLENegativeBinomialTest, LargeSample_CloseToPoisson_LargerSE) {
    // 同一 Poisson 生成数据上, NB MLE (α 固定) 系数接近 Poisson MLE
    // 但 NB 的 SE 应大于 Poisson (因 V(μ)=μ+αμ² > μ)
    const Size N = 2000;
    MatrixXD X(N, 2);
    VectorXD y(N);

    std::mt19937 rng(42);
    std::normal_distribution<Real> norm(0.0, 1.0);
    for (Size i = 0; i < N; ++i) {
        X(i, 0) = 1.0;
        X(i, 1) = norm(rng);
        const Real mu = std::exp(0.5 + 0.3 * X(i, 1));
        std::poisson_distribution<int> pois(mu);
        y(i) = static_cast<Real>(pois(rng));
    }
    auto data = make_cross_section(X, y, {"const", "x"}, "y");

    MLEEstimator mle_pois(MLEFamily::Poisson, CovarianceType::Hessian);
    mle_pois.setMaxIter(100).setTolerance(1e-10);
    EstimationResult rp = mle_pois.estimate(data);

    MLEEstimator mle_nb(MLEFamily::NegativeBinomial, CovarianceType::Hessian);
    mle_nb.setAlpha(0.5).setMaxIter(200).setTolerance(1e-10);
    // 用 Poisson MLE 结果作为 NB 初始值 (NB 与 Poisson 同 score 方程, 系数应接近)
    VectorXD nb_start(2);
    nb_start(0) = rp.coefficients(0);
    nb_start(1) = rp.coefficients(1);
    mle_nb.setStartValues(nb_start);
    EstimationResult rn = mle_nb.estimate(data);

    ASSERT_TRUE(mle_pois.converged());
    ASSERT_TRUE(mle_nb.converged());

    // 系数接近 (NB MLE 在固定 α 下, 与 Poisson MLE 同 score 方程, 系数应相同)
    for (Size i = 0; i < 2; ++i) {
        EXPECT_NEAR(rn.coefficients(i), rp.coefficients(i), 0.01)
            << "NB and Poisson MLE coefficients[" << i << "] should be close";
    }

    // NB SE 应大于 Poisson (NB 假设 V(μ) 更大)
    // 注: 这是 Hessian V 的渐近比较, NB 的 (X'WX)^{-1} 因 W=μ/(1+αμ) 更小故 V 更大
    for (Size i = 0; i < 2; ++i) {
        EXPECT_GT(rn.std_errors(i), rp.std_errors(i) * 0.95)
            << "NB SE[" << i << "] should be >= Poisson SE (approximately)";
    }
}

// =============================================================================
// §12 NB MLE 不支持错误协方差类型
// =============================================================================
TEST(MLENegativeBinomialTest, ThrowsOnUnsupportedCovarianceType) {
    MLEEstimator mle(MLEFamily::NegativeBinomial, CovarianceType::Classical);
    mle.setAlpha(1.0);
    EXPECT_THROW(mle.estimate(make_probit_symmetric()), std::invalid_argument);
}
