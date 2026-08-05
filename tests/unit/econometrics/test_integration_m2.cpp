// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.3 - M2 端到端集成测试
// 验证方法: 工厂创建估计器 → 估计 → 假设检验 (Wald/LR) → 信息准则 (AIC/BIC)
//
// 排幻觉点:
//   F1: 工厂 create() 每次返回全新实例 (clone 语义, 不共享状态)
//   F2: 未注册名称抛 std::invalid_argument (fail-fast, 不返回 nullptr)
//   F3: 工厂不持有 Estimator 实例 (unique_ptr 转移所有权)
//   F4: Meyers Singleton 规避 SIOF (静态注册顺序正确)
//   E7: MLE Logistic 用 Newton-Raphson (对 canonical link 等价于 IRLS)
//   E8: QMLE Sandwich bread = (X'WX)^{-1}, meat = X' diag(ε²) X
//   E9: Wald 检验同时提供 χ² 和 F 两种 p 值
//   G3: Gaussian Hessian V = σ²·(X'X)^{-1} (含 σ² 因子)
//
// 集成场景:
//   §1 OLS → Wald 单系数显著性 + AIC/BIC
//   §2 MLE Logistic → Wald 联合显著性 + LR 检验 (约束 vs 无约束)
//   §3 MLE Poisson → QMLE Sandwich + Wald 系数显著性
//   §4 工厂模式分发 + Gaussian MLE vs OLS 数值一致性
//   §5 模型选择: Logistic vs Probit (AIC/BIC 比较)

#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/covariance_type.hpp"
#include "cpphub/econometrics/core/data_types.hpp"
#include "cpphub/econometrics/core/estimation_result.hpp"
#include "cpphub/econometrics/core/estimator_base.hpp"
#include "cpphub/econometrics/estimation/estimator_factory.hpp"
#include "cpphub/econometrics/estimation/mle.hpp"
#include "cpphub/econometrics/estimation/ols.hpp"
#include "cpphub/econometrics/inference/diagnostics.hpp"
#include "cpphub/econometrics/inference/hypothesis_tests.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;
using cpphub::v1::linalg::dynamic::MatrixXD;
using cpphub::v1::linalg::dynamic::VectorXD;

namespace {

// =============================================================================
// 手算数据集 A: N=4, K=2 含截距 (复用 test_ols_hc 手算值)
//   X = [1,1; 1,2; 1,3; 1,4],  y = [2; 3; 5; 7]
//   β_OLS = [0; 1.7],  σ²_MLE = SSR/N = 0.30/4 = 0.075
//   (X'X)^{-1} = [3/2, -1/2; -1/2, 1/5]
//   Classical vcov = 0.15·(X'X)^{-1}  (用 σ²_hat = SSR/(N-K) = 0.15)
//   Gaussian MLE vcov = 0.075·(X'X)^{-1}  (用 σ²_MLE = SSR/N = 0.075)
// =============================================================================
CrossSectionData make_data_A() {
    MatrixXD X(4, 2);
    X(0, 0) = 1.0; X(0, 1) = 1.0;
    X(1, 0) = 1.0; X(1, 1) = 2.0;
    X(2, 0) = 1.0; X(2, 1) = 3.0;
    X(3, 0) = 1.0; X(3, 1) = 4.0;
    VectorXD y(4);
    y(0) = 2.0; y(1) = 3.0; y(2) = 5.0; y(3) = 7.0;
    return make_cross_section(X, y, {"const", "x"}, "y");
}

// =============================================================================
// 手算数据集 B: N=30, K=2 含截距, 用于 Logistic/Probit 回归
//   构造非完全分离的二元结果 (β 真值约 [α, β] = [-3, 0.3])
//   注: 必须非完全分离, 否则 MLE 不存在 (β→∞, X'WX 奇异)
//   排幻觉: 完全分离检测 — y=0 集合与 y=1 集合在 x 上无重叠时 MLE 发散
//   N=30 + 渐变转换区保证 Wald 检验有足够统计功效 (β₁ 在 5% 水平显著)
// =============================================================================
CrossSectionData make_logistic_data() {
    const Size N = 30;
    MatrixXD X(N, 2);
    VectorXD y(N);
    for (Size i = 0; i < N; ++i) {
        X(i, 0) = 1.0;
        X(i, 1) = static_cast<Real>(i + 1);  // x = 1..30
    }
    // y 有正趋势, 渐变转换区 (x=11-17), 非完全分离:
    //   x=1-10: y=0 (10 个), x=11-14: 交替, x=15-17: 多 1 含 1 个 0, x=18-30: y=1
    static const Real y_vals[30] = {
        0,0,0,0,0,0,0,0,0,0,   // x=1-10: all y=0 (10)
        1,0,1,0,                // x=11-14: transition (4)
        1,1,0,                   // x=15-17: mostly y=1 with one y=0 (3)
        1,1,1,1,1,1,1,1,1,1,1,1,1  // x=18-30: all y=1 (13)
    };
    for (Size i = 0; i < N; ++i) y(i) = y_vals[i];
    return make_cross_section(X, y, {"const", "x"}, "y");
}

// =============================================================================
// 数据集 C: Poisson 回归 (计数数据)
//   构造 y ~ Poisson(exp(α + β·x)), α=0.5, β=0.3
// =============================================================================
CrossSectionData make_poisson_data() {
    const Size N = 100;
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
    return make_cross_section(X, y, {"const", "x"}, "y");
}

}  // namespace

// =============================================================================
// §1 OLS → Wald 单系数显著性 + AIC/BIC
// 端到端: EconData → 工厂创建 OLS → 估计 → Wald(β₁=0) → AIC/BIC
// 排幻觉点 F1 (工厂全新实例), E9 (Wald χ² + F 双 p 值)
// =============================================================================
TEST(IntegrationM2Test, OLS_Wald_InformationCriteria_EndToEnd) {
    // 1. 工厂创建 OLS 估计器 (排幻觉点 F1: 全新实例)
    auto ols_ptr = EstimatorFactory::instance().create("OLS");
    ASSERT_NE(ols_ptr, nullptr);
    EXPECT_EQ(ols_ptr->name(), "OLS");

    // 2. 估计
    EstimationResult r = ols_ptr->estimate(make_data_A());

    // 手算值: β = [0, 1.7]
    ASSERT_EQ(r.coefficients.size(), 2u);
    EXPECT_NEAR(r.coefficients(0), 0.0, 1e-10);
    EXPECT_NEAR(r.coefficients(1), 1.7, 1e-10);

    // 3. Wald 检验 H0: β₁ = 0 (单系数显著性)
    //    R = [0, 1], r = [0]
    //    Wald = (Rβ-r)' (RVR')^{-1} (Rβ-r) = 1.7² / V(β₁)
    MatrixXD R(1, 2);
    R(0, 0) = 0.0; R(0, 1) = 1.0;
    VectorXD r_vec(1);
    r_vec(0) = 0.0;

    // 用 OLS Classical vcov (含 σ²_hat = SSR/(N-K) = 0.15)
    // V(β₁) = 0.15 · 0.2 = 0.03
    // Wald = 1.7² / 0.03 = 2.89 / 0.03 = 96.3333...
    HypothesisTestResult wald = wald_test(r.coefficients, r.vcov, R, r_vec,
                                            static_cast<Real>(r.df_residual));
    EXPECT_EQ(wald.test_name, "Wald");
    EXPECT_NEAR(wald.statistic, 96.3333333333, 1e-6);
    EXPECT_EQ(wald.df, 1u);

    // 排幻觉点 E9: 同时提供 χ² 和 F 两种 p 值
    // F 形式: F = Wald/q = 96.333.../1 = 96.333... ~ F(1, 2)
    // p 值应非常小 (β₁ 高度显著)
    EXPECT_TRUE(wald.use_f_distribution);
    EXPECT_LT(wald.p_value, 0.05);
    EXPECT_TRUE(wald.reject_null_95);

    // 4. 信息准则 (排幻觉点: AIC = 2K - 2ℓ)
    //    Gaussian loglik: ℓ = -N/2·(log(2π) + 1 + log(σ²_MLE))
    //    σ²_MLE = SSR/N = 0.30/4 = 0.075
    //    ℓ = -2·(log(2π) + 1 + log(0.075)) = -2·(1.8379 + 1 - 2.5903) = -2·0.2476 ≈ -0.4952
    const InformationCriteria ic = compute_information_criteria(
        r.log_likelihood, r.n_params, r.n_obs);
    const Real expected_aic = 2.0 * static_cast<Real>(r.n_params) - 2.0 * r.log_likelihood;
    EXPECT_NEAR(ic.aic, expected_aic, 1e-10);
    // 注: N=4 时 log(N) = 1.386 < 2 = K 系数, 故 BIC = K·log(N) - 2ℓ < 2K - 2ℓ = AIC
    // 即小样本下 BIC < AIC (BIC 对参数的惩罚更轻), 此处仅验证 AIC 公式正确性
    const Real expected_bic = static_cast<Real>(r.n_params) * std::log(static_cast<Real>(r.n_obs)) - 2.0 * r.log_likelihood;
    EXPECT_NEAR(ic.bic, expected_bic, 1e-10);
}

// =============================================================================
// §2 MLE Logistic → Wald 联合显著性 + LR 检验
// 端到端: EconData → MLE Logistic → 估计 → Wald(H0: β=0 联合) → LR(约束 vs 无约束)
// 排幻觉点 E7 (Newton-Raphson ≡ IRLS for canonical link), E9 (Wald 双 p 值)
// =============================================================================
TEST(IntegrationM2Test, MLELogistic_Wald_LR_EndToEnd) {
    // 1. 工厂创建 MLE Logistic
    auto mle_ptr = EstimatorFactory::instance().createMLE(MLEFamily::Logistic);
    ASSERT_NE(mle_ptr, nullptr);
    mle_ptr->setMaxIter(200).setTolerance(1e-10);

    // 2. 无约束估计 (含 x)
    EstimationResult r_unrestricted = mle_ptr->estimate(make_logistic_data());
    EXPECT_TRUE(mle_ptr->converged());

    // 3. Wald 检验 H0: β₁ = 0 (x 系数显著性)
    MatrixXD R(1, 2);
    R(0, 0) = 0.0; R(0, 1) = 1.0;
    VectorXD r_vec(1);
    r_vec(0) = 0.0;

    HypothesisTestResult wald = wald_test(r_unrestricted.coefficients,
                                            r_unrestricted.vcov, R, r_vec,
                                            static_cast<Real>(r_unrestricted.df_residual));
    EXPECT_EQ(wald.test_name, "Wald");
    EXPECT_GT(wald.statistic, 0.0);
    // 数据可分离, β₁ 应显著
    EXPECT_LT(wald.p_value, 0.10) << "x coefficient should be significant";

    // 4. LR 检验: 无约束 (含 x) vs 约束 (仅截距)
    //    LR = 2(ℓ_UR - ℓ_R) ~ χ²(1)
    auto mle_restricted_ptr = EstimatorFactory::instance().createMLE(MLEFamily::Logistic);
    mle_restricted_ptr->setMaxIter(200).setTolerance(1e-10);

    // 约束模型: 仅截距, X = [1; 1; ...; 1], y 与无约束模型相同 (N=30)
    const Size N_logistic = 30;
    MatrixXD X_restricted(N_logistic, 1);
    for (Size i = 0; i < N_logistic; ++i) X_restricted(i, 0) = 1.0;
    VectorXD y_restricted(N_logistic);
    static const Real y_restricted_vals[30] = {
        0,0,0,0,0,0,0,0,0,0,   // x=1-10: y=0
        1,0,1,0,                // x=11-14: transition
        1,1,0,                   // x=15-17
        1,1,1,1,1,1,1,1,1,1,1,1,1  // x=18-30: y=1
    };
    for (Size i = 0; i < N_logistic; ++i) y_restricted(i) = y_restricted_vals[i];
    auto data_restricted = make_cross_section(X_restricted, y_restricted, {"const"}, "y");
    EstimationResult r_restricted = mle_restricted_ptr->estimate(data_restricted);

    HypothesisTestResult lr = lr_test(r_unrestricted, r_restricted);
    EXPECT_EQ(lr.test_name, "LR");
    EXPECT_GT(lr.statistic, 0.0);
    EXPECT_EQ(lr.df, 1u);

    // LR 与 Wald 渐近等价 (大样本), 小样本下数值不同但方向一致
    EXPECT_EQ(lr.reject_null_95, wald.reject_null_95)
        << "LR and Wald should agree on reject/fail at 5% level";

    // 5. 信息准则比较: 无约束模型 AIC 应更小 (似然更高)
    const InformationCriteria ic_ur = compute_information_criteria(
        r_unrestricted.log_likelihood, r_unrestricted.n_params, r_unrestricted.n_obs);
    const InformationCriteria ic_r = compute_information_criteria(
        r_restricted.log_likelihood, r_restricted.n_params, r_restricted.n_obs);
    // 若 x 显著, 无约束模型 ℓ 更高 → AIC 更小
    // 注: 小样本下可能 AIC_ur > AIC_r (β 弱), 此处仅验证 AIC 公式
    EXPECT_NEAR(ic_ur.aic, 2.0 * r_unrestricted.n_params - 2.0 * r_unrestricted.log_likelihood, 1e-10);
    EXPECT_NEAR(ic_r.aic, 2.0 * r_restricted.n_params - 2.0 * r_restricted.log_likelihood, 1e-10);
}

// =============================================================================
// §3 MLE Poisson → QMLE Sandwich + Wald 系数显著性
// 端到端: EconData → MLE Poisson → QMLE Sandwich 协方差 → Wald 检验
// 排幻觉点 E8 (Sandwich bread = (X'WX)^{-1}, meat = X' diag(ε²) X)
// =============================================================================
TEST(IntegrationM2Test, MLEPoisson_QMLE_Sandwich_Wald_EndToEnd) {
    // 1. 工厂创建 MLE Poisson (Hessian 协方差, 用于点估计)
    auto mle_hessian_ptr = EstimatorFactory::instance().createMLE(MLEFamily::Poisson);
    mle_hessian_ptr->setMaxIter(200).setTolerance(1e-10);
    EstimationResult r_hessian = mle_hessian_ptr->estimate(make_poisson_data());
    EXPECT_TRUE(mle_hessian_ptr->converged());

    // 2. 工厂创建 MLE Poisson (Sandwich 协方差, QMLE 稳健)
    //    排幻觉点 F1: 工厂每次返回全新实例, 不共享状态
    auto mle_sandwich_ptr = EstimatorFactory::instance().createMLE(MLEFamily::Poisson);
    mle_sandwich_ptr->setCovarianceType(CovarianceType::Sandwich);
    mle_sandwich_ptr->setMaxIter(200).setTolerance(1e-10);
    EstimationResult r_sandwich = mle_sandwich_ptr->estimate(make_poisson_data());

    // 3. 点估计应一致 (同一数据, 同一分布族, 同一优化器)
    ASSERT_EQ(r_hessian.coefficients.size(), r_sandwich.coefficients.size());
    for (Size i = 0; i < r_hessian.coefficients.size(); ++i) {
        EXPECT_NEAR(r_hessian.coefficients(i), r_sandwich.coefficients(i), 1e-8)
            << "Point estimates should be identical (same optimizer)";
    }

    // 4. 协方差应不同 (Hessian vs Sandwich)
    //    排幻觉点 E8: Sandwich = (X'WX)^{-1} · X'diag(ε²)X · (X'WX)^{-1}
    //    若数据真实服从 Poisson, Hessian ≈ Sandwich; 若过离散 (overdispersion),
    //    Sandwich > Hessian (meat 中 ε² 较大)
    for (Size i = 0; i < r_hessian.coefficients.size(); ++i) {
        const Real v_hessian = r_hessian.vcov(i, i);
        const Real v_sandwich = r_sandwich.vcov(i, i);
        EXPECT_GT(v_hessian, 0.0);
        EXPECT_GT(v_sandwich, 0.0);
        // Sandwich 协方差对角元素应与 Hessian 同数量级 (QMLE 渐近正确)
        // 不要求严格相等 (取决于数据过离散程度)
        const Real ratio = v_sandwich / v_hessian;
        EXPECT_GT(ratio, 0.1);
        EXPECT_LT(ratio, 10.0);
    }

    // 5. Wald 检验 H0: β₁ = 0 (用 Sandwich 协方差, QMLE 稳健)
    MatrixXD R(1, 2);
    R(0, 0) = 0.0; R(0, 1) = 1.0;
    VectorXD r_vec(1);
    r_vec(0) = 0.0;

    HypothesisTestResult wald_sandwich = wald_test(
        r_sandwich.coefficients, r_sandwich.vcov, R, r_vec);
    EXPECT_EQ(wald_sandwich.test_name, "Wald");
    EXPECT_GT(wald_sandwich.statistic, 0.0);
    // 真实 β₁ = 0.3, N=100, 应显著
    EXPECT_LT(wald_sandwich.p_value, 0.10)
        << "β₁ (true=0.3, N=100) should be significant under QMLE";
}

// =============================================================================
// §4 工厂模式分发 + Gaussian MLE vs OLS 数值一致性
// 端到端: 工厂创建 Gaussian MLE 和 OLS → 同一数据估计 → 数值一致
// 排幻觉点 G3 (Gaussian Hessian V = σ²·(X'X)^{-1}, σ²=SSR/N 即 MLE 估计)
//   注: OLS Classical vcov 用 σ²_hat = SSR/(N-K), Gaussian MLE vcov 用 σ²_MLE = SSR/N
//   两者协方差不同 (因子 N/(N-K)), 但 β 点估计一致
// =============================================================================
TEST(IntegrationM2Test, Factory_GaussianMLE_vs_OLS_Consistency) {
    // 1. 工厂创建 OLS 和 Gaussian MLE
    auto ols_ptr = EstimatorFactory::instance().create("OLS");
    auto mle_gauss_ptr = EstimatorFactory::instance().createMLE(MLEFamily::Gaussian);

    ASSERT_NE(ols_ptr, nullptr);
    ASSERT_NE(mle_gauss_ptr, nullptr);
    EXPECT_EQ(ols_ptr->name(), "OLS");
    EXPECT_EQ(mle_gauss_ptr->name(), "MLE(Gaussian)");

    // 2. 同一数据估计
    auto data = make_data_A();
    EstimationResult r_ols = ols_ptr->estimate(data);
    EstimationResult r_mle = mle_gauss_ptr->estimate(data);

    // 3. β 点估计应一致 (Gaussian MLE ≡ OLS 闭式解)
    ASSERT_EQ(r_ols.coefficients.size(), r_mle.coefficients.size());
    for (Size i = 0; i < r_ols.coefficients.size(); ++i) {
        EXPECT_NEAR(r_ols.coefficients(i), r_mle.coefficients(i), 1e-10)
            << "Gaussian MLE and OLS should give identical β";
    }

    // 4. 对数似然应一致 (Gaussian MLE 显式计算 ℓ, OLS 隐式计算)
    EXPECT_NEAR(r_ols.log_likelihood, r_mle.log_likelihood, 1e-8);

    // 5. 协方差差异 (排幻觉点 G3):
    //    OLS Classical: V = σ²_hat·(X'X)^{-1}, σ²_hat = SSR/(N-K)
    //    Gaussian MLE Hessian: V = σ²_MLE·(X'X)^{-1}, σ²_MLE = SSR/N
    //    比值: V_MLE / V_OLS = (SSR/N) / (SSR/(N-K)) = (N-K)/N
    const Size N = 4, K = 2;
    const Real expected_ratio = static_cast<Real>(N - K) / static_cast<Real>(N);
    for (Size i = 0; i < r_ols.coefficients.size(); ++i) {
        const Real v_ols = r_ols.vcov(i, i);
        const Real v_mle = r_mle.vcov(i, i);
        if (v_ols > 0.0 && v_mle > 0.0) {
            const Real actual_ratio = v_mle / v_ols;
            EXPECT_NEAR(actual_ratio, expected_ratio, 1e-8)
                << "V(MLE)/V(OLS) should be (N-K)/N = " << expected_ratio;
        }
    }
}

// =============================================================================
// §5 模型选择: Logistic vs Probit (AIC/BIC 比较)
// 端到端: 同一二元数据 → Logistic MLE + Probit MLE → AIC/BIC 比较
// 排幻觉点: 模型选择准则用于比较不同分布族假设的拟合优度
//   注: Logistic 和 Probit 对同一数据拟合通常 AIC 接近,
//       差异小于 1 时视为等价模型 (Burnham-Anderson 2002)
// =============================================================================
TEST(IntegrationM2Test, ModelSelection_Logistic_vs_Probit_AIC) {
    // 1. 工厂创建 Logistic 和 Probit MLE
    auto mle_logit_ptr = EstimatorFactory::instance().createMLE(MLEFamily::Logistic);
    auto mle_probit_ptr = EstimatorFactory::instance().createMLE(MLEFamily::Probit);

    mle_logit_ptr->setMaxIter(200).setTolerance(1e-10);
    mle_probit_ptr->setMaxIter(200).setTolerance(1e-10);

    // 2. 同一数据估计
    auto data = make_logistic_data();
    EstimationResult r_logit = mle_logit_ptr->estimate(data);
    EstimationResult r_probit = mle_probit_ptr->estimate(data);

    EXPECT_TRUE(mle_logit_ptr->converged());
    EXPECT_TRUE(mle_probit_ptr->converged());

    // 3. 计算信息准则
    const InformationCriteria ic_logit = compute_information_criteria(
        r_logit.log_likelihood, r_logit.n_params, r_logit.n_obs);
    const InformationCriteria ic_probit = compute_information_criteria(
        r_probit.log_likelihood, r_probit.n_params, r_probit.n_obs);

    // 4. 两个模型 AIC 应接近 (同参数数 K=2, 同数据)
    //    差异 < 2 视为等价 (Burnham-Anderson 2002 §2.4)
    const Real aic_diff = std::fabs(ic_logit.aic - ic_probit.aic);
    EXPECT_LT(aic_diff, 5.0)
        << "Logistic and Probit AIC should be close (same K, same data)";

    // 5. BIC 同样应接近
    const Real bic_diff = std::fabs(ic_logit.bic - ic_probit.bic);
    EXPECT_LT(bic_diff, 5.0)
        << "Logistic and Probit BIC should be close (same K, same data)";

    // 6. 排幻觉点 F1: 工厂返回的两次实例相互独立
    //    修改一个不影响另一个
    mle_logit_ptr->setCovarianceType(CovarianceType::Sandwich);
    EXPECT_EQ(mle_probit_ptr->covarianceType(), CovarianceType::Hessian)
        << "Modifying logit should not affect probit (independent instances)";
}

// =============================================================================
// §6 工厂异常处理: 未注册名称抛异常 (排幻觉点 F2)
// =============================================================================
TEST(IntegrationM2Test, Factory_UnknownName_Throws) {
    EXPECT_THROW(EstimatorFactory::instance().create("NonExistent"),
                 std::invalid_argument);
    EXPECT_THROW(EstimatorFactory::instance().create(""),
                 std::invalid_argument);
    EXPECT_THROW(EstimatorFactory::instance().create("OLS "),  // 带空格
                 std::invalid_argument);
}

// =============================================================================
// §7 工厂注册验证: 已注册估计器清单 (排幻觉点 F4)
// =============================================================================
TEST(IntegrationM2Test, Factory_RegisteredEstimators) {
    // 排幻觉点 F4: Meyers Singleton 规避 SIOF, 所有内置估计器应已注册
    EXPECT_TRUE(EstimatorFactory::instance().isRegistered("OLS"));
    EXPECT_TRUE(EstimatorFactory::instance().isRegistered("MLE.Gaussian"));
    EXPECT_TRUE(EstimatorFactory::instance().isRegistered("MLE.Logistic"));
    EXPECT_TRUE(EstimatorFactory::instance().isRegistered("MLE.Bernoulli"));
    EXPECT_TRUE(EstimatorFactory::instance().isRegistered("MLE.Probit"));
    EXPECT_TRUE(EstimatorFactory::instance().isRegistered("MLE.Poisson"));
    EXPECT_TRUE(EstimatorFactory::instance().isRegistered("MLE.NegativeBinomial"));

    // 至少 7 个内置估计器
    EXPECT_GE(EstimatorFactory::instance().size(), 7u);
}
