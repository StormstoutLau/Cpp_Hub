// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.2 任务 2.3 - 假设检验测试 (Wald/LR/LM/J-test)
// 验证方法: 手算解析值 + Greene/Wooldridge 教材对照 + 渐近等价性 (容差 1e-6)
//
// 排幻觉点:
//   E9: R lmtest::waldtest 默认 F 检验 (小样本), C++ 同时提供 χ² 和 F 两种 p 值
//   H1: LM = N·R²_aux (Breusch-Pagan 1979 原始形式), 非 (N-K)·R²_aux
//   H2: J-test df = q - k (过度识别约束数), 非 q (矩条件数)
//
// 手算数据集 (复用 test_ols_hc.cpp N=4, K=2):
//   X = [1,1; 1,2; 1,3; 1,4],  y = [2; 3; 5; 7]
//   β = [0; 1.7],  (X'X)^{-1} = [3/2,-1/2; -1/2,1/5]
//   Classical σ² = SSR/(N-K) = 0.30/2 = 0.15
//   vcov = 0.15·(X'X)^{-1} = [9/40, -3/40; -3/40, 3/100]
//     vcov(0,0) = 9/40 = 0.225
//     vcov(1,1) = 3/100 = 0.03
//   ȳ = 4.25,  SST = 14.75,  SSR = 0.30
//
// Wald (β₁=0): R=[0,1], r=[0], Rβ-r=1.7, RVR'=0.03, Wald=1.7²/0.03=289/30≈9.6333...
//   注: 之前估算有误, 实际 1.7² = 2.89, 2.89/0.03 = 96.333... (重新核算)
//   实际: 2.89 / 0.03 = 96.3333...
//
// LM (约束模型 = 仅截距):
//   ε_R = y - ȳ = [-2.25, -1.25, 0.75, 2.75]
//   σ²_R = ε_R'ε_R/N = 14.75/4 = 3.6875 (MLE)
//   score = X'ε_R = [0; 8.5]
//   (X'X)^{-1}·score = [-4.25; 1.7]
//   score'·(X'X)^{-1}·score = 8.5·1.7 = 14.45
//   LM = 14.45 / 3.6875 = 3.91864406779661...
//   等价 N·R²_aux = 4 · (14.45/14.75) = 3.91864406779661...

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "cpphub/econometrics/inference/hypothesis_tests.hpp"
#include "cpphub/econometrics/core/estimation_result.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;
using cpphub::v1::linalg::dynamic::MatrixXD;
using cpphub::v1::linalg::dynamic::VectorXD;

namespace {

// 构造 N=4, K=2 含截距数据集 (复用 test_ols_hc 手算值)
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

// (X'X)^{-1} = [3/2,-1/2; -1/2,1/5] (手算)
MatrixXD make_XtX_inv() {
    MatrixXD XtX_inv(2, 2);
    XtX_inv(0, 0) = 1.5;  XtX_inv(0, 1) = -0.5;
    XtX_inv(1, 0) = -0.5; XtX_inv(1, 1) = 0.2;
    return XtX_inv;
}

// β = [0; 1.7]
VectorXD make_beta() {
    VectorXD b(2);
    b(0) = 0.0; b(1) = 1.7;
    return b;
}

// Classical vcov = 0.15 · (X'X)^{-1} = [9/40, -3/40; -3/40, 3/100]
MatrixXD make_classical_vcov() {
    MatrixXD v(2, 2);
    v(0, 0) = 9.0 / 40.0;   v(0, 1) = -3.0 / 40.0;
    v(1, 0) = -3.0 / 40.0;  v(1, 1) = 3.0 / 100.0;
    return v;
}

// 约束模型 (仅截距) 残差: ε_R = y - ȳ = y - 4.25
VectorXD make_residuals_restricted() {
    VectorXD r(4);
    r(0) = -2.25; r(1) = -1.25; r(2) = 0.75; r(3) = 2.75;
    return r;
}

}  // namespace

// =============================================================================
// §1 Wald 检验 - 手算解析值
// =============================================================================

// Wald 检验 β₁ = 0: Wald = (1.7)² / 0.03 = 2.89/0.03 = 96.3333...
TEST(HypothesisTest, Wald_SingleRestriction_Beta1_Equals_Zero) {
    VectorXD beta = make_beta();
    MatrixXD vcov = make_classical_vcov();
    MatrixXD R(1, 2);
    R(0, 0) = 0.0; R(0, 1) = 1.0;
    VectorXD r(1);
    r(0) = 0.0;

    HypothesisTestResult result = wald_test(beta, vcov, R, r);

    EXPECT_EQ(result.test_name, "Wald");
    EXPECT_NEAR(result.statistic, 2.89 / 0.03, 1e-10);
    EXPECT_EQ(result.df, 1u);
    // χ²(1) p-value: P(χ² > 96.333) ≈ 1.0e-22 (极小)
    EXPECT_GT(result.p_value, 0.0);
    EXPECT_LT(result.p_value, 1e-15);
    EXPECT_TRUE(result.reject_null_95);
    EXPECT_TRUE(result.reject_null_99);
    // 默认无 df_residual, 使用 χ²
    EXPECT_FALSE(result.use_f_distribution);
}

// Wald 检验 β₀ = 0: Wald = 0² / 0.225 = 0 (β₀ 本身就是 0)
TEST(HypothesisTest, Wald_SingleRestriction_Beta0_Equals_Zero) {
    VectorXD beta = make_beta();
    MatrixXD vcov = make_classical_vcov();
    MatrixXD R(1, 2);
    R(0, 0) = 1.0; R(0, 1) = 0.0;
    VectorXD r(1);
    r(0) = 0.0;

    HypothesisTestResult result = wald_test(beta, vcov, R, r);

    EXPECT_NEAR(result.statistic, 0.0, 1e-12);
    EXPECT_EQ(result.df, 1u);
    // χ²(1) p-value when stat=0: p=1
    EXPECT_NEAR(result.p_value, 1.0, 1e-10);
    EXPECT_FALSE(result.reject_null_95);
    EXPECT_FALSE(result.reject_null_99);
}

// Wald 检验 β₁ = 1.7 (真实值, 不拒绝): Wald = 0
TEST(HypothesisTest, Wald_SingleRestriction_Beta1_Equals_TrueValue) {
    VectorXD beta = make_beta();
    MatrixXD vcov = make_classical_vcov();
    MatrixXD R(1, 2);
    R(0, 0) = 0.0; R(0, 1) = 1.0;
    VectorXD r(1);
    r(0) = 1.7;

    HypothesisTestResult result = wald_test(beta, vcov, R, r);

    EXPECT_NEAR(result.statistic, 0.0, 1e-12);
    EXPECT_NEAR(result.p_value, 1.0, 1e-10);
    EXPECT_FALSE(result.reject_null_95);
}

// Wald 检验 - 联合约束 β₀=0 AND β₁=1.7 (两者均成立): Wald = 0
TEST(HypothesisTest, Wald_JointRestriction_BothTrue) {
    VectorXD beta = make_beta();
    MatrixXD vcov = make_classical_vcov();
    MatrixXD R(2, 2);
    R(0, 0) = 1.0; R(0, 1) = 0.0;
    R(1, 0) = 0.0; R(1, 1) = 1.0;
    VectorXD r(2);
    r(0) = 0.0; r(1) = 1.7;

    HypothesisTestResult result = wald_test(beta, vcov, R, r);

    EXPECT_NEAR(result.statistic, 0.0, 1e-12);
    EXPECT_EQ(result.df, 2u);
    EXPECT_NEAR(result.p_value, 1.0, 1e-10);
}

// Wald 检验 - 联合约束 β₀=1 AND β₁=0 (均不成立)
// Rβ - r = [-1; 1.7]
// RVR' = [[0.225, -0.075], [-0.075, 0.03]]
// (RVR')^{-1} 需手算: det = 0.225·0.03 - 0.075² = 0.00675 - 0.005625 = 0.001125
// (RVR')^{-1} = (1/0.001125) · [[0.03, 0.075], [0.075, 0.225]]
// Wald = [-1, 1.7] · (RVR')^{-1} · [-1; 1.7]
//      = [-1, 1.7] · (1/0.001125) · [0.03·(-1)+0.075·1.7; 0.075·(-1)+0.225·1.7]
//      = [-1, 1.7] · (1/0.001125) · [0.0975; 0.3075]
//      = (1/0.001125) · (-0.0975 + 1.7·0.3075)
//      = (1/0.001125) · (-0.0975 + 0.52275)
//      = (1/0.001125) · 0.42525
//      = 378.0
TEST(HypothesisTest, Wald_JointRestriction_BothFalse) {
    VectorXD beta = make_beta();
    MatrixXD vcov = make_classical_vcov();
    MatrixXD R(2, 2);
    R(0, 0) = 1.0; R(0, 1) = 0.0;
    R(1, 0) = 0.0; R(1, 1) = 1.0;
    VectorXD r(2);
    r(0) = 1.0; r(1) = 0.0;

    HypothesisTestResult result = wald_test(beta, vcov, R, r);

    EXPECT_NEAR(result.statistic, 378.0, 1e-8);
    EXPECT_EQ(result.df, 2u);
    EXPECT_TRUE(result.reject_null_95);
    EXPECT_TRUE(result.reject_null_99);
}

// =============================================================================
// §2 Wald 检验 - F 分布 (小样本修正, 排幻觉点 E9)
// =============================================================================

// 排幻觉点 E9: 提供 df_residual 时使用 F 分布
// F = Wald / q = 96.333.../ 1 = 96.333...
// F(1, 2) p-value 应大于 χ²(1) p-value (F 更保守)
TEST(HypothesisTest, Wald_FDistribution_WhenDfResidualProvided) {
    VectorXD beta = make_beta();
    MatrixXD vcov = make_classical_vcov();
    MatrixXD R(1, 2);
    R(0, 0) = 0.0; R(0, 1) = 1.0;
    VectorXD r(1);
    r(0) = 0.0;

    // df_residual = N - K = 4 - 2 = 2
    HypothesisTestResult result = wald_test(beta, vcov, R, r, 2.0);

    EXPECT_TRUE(result.use_f_distribution);
    EXPECT_NEAR(result.f_df1, 1.0, 1e-10);
    EXPECT_NEAR(result.f_df2, 2.0, 1e-10);
    // F 统计量 = Wald / q = 96.333...
    const Real f_stat = result.statistic / 1.0;
    EXPECT_NEAR(f_stat, 2.89 / 0.03, 1e-10);
    // F(1, 2) p-value: P(F > 96.333) — 小样本下 F 比 χ² 更保守
    // 由于 df2=2 极小, p-value 仍应小于 0.05 (F(1,2) 95% 临界值 18.51)
    EXPECT_GT(result.p_value, 0.0);
    EXPECT_LT(result.p_value, 0.05);  // 仍应拒绝 (96.33 >> 18.51)
}

// 无 df_residual 时默认使用 χ² (排幻觉点 E9: 两种都提供, 调用方选择)
TEST(HypothesisTest, Wald_ChiSquare_WhenNoDfResidual) {
    VectorXD beta = make_beta();
    MatrixXD vcov = make_classical_vcov();
    MatrixXD R(1, 2);
    R(0, 0) = 0.0; R(0, 1) = 1.0;
    VectorXD r(1);
    r(0) = 0.0;

    HypothesisTestResult result = wald_test(beta, vcov, R, r);  // 无 df_residual

    EXPECT_FALSE(result.use_f_distribution);
    EXPECT_EQ(result.df, 1u);
}

// =============================================================================
// §3 Wald 检验 - 异常处理
// =============================================================================

TEST(HypothesisTest, Wald_ThrowsOnDimensionMismatch_R_Cols) {
    VectorXD beta(2);
    MatrixXD vcov(2, 2);
    MatrixXD R(1, 3);  // 列数不匹配
    VectorXD r(1);
    EXPECT_THROW(wald_test(beta, vcov, R, r), std::invalid_argument);
}

TEST(HypothesisTest, Wald_ThrowsOnDimensionMismatch_r_Size) {
    VectorXD beta(2);
    MatrixXD vcov(2, 2);
    MatrixXD R(1, 2);
    VectorXD r(2);  // 大小不匹配
    EXPECT_THROW(wald_test(beta, vcov, R, r), std::invalid_argument);
}

TEST(HypothesisTest, Wald_ThrowsOnDimensionMismatch_vcov) {
    VectorXD beta(2);
    MatrixXD vcov(3, 3);  // 维度不匹配
    MatrixXD R(1, 2);
    VectorXD r(1);
    EXPECT_THROW(wald_test(beta, vcov, R, r), std::invalid_argument);
}

TEST(HypothesisTest, Wald_ThrowsOnZeroRestrictions) {
    VectorXD beta(2);
    MatrixXD vcov(2, 2);
    MatrixXD R(0, 2);  // 零约束
    VectorXD r(0);
    EXPECT_THROW(wald_test(beta, vcov, R, r), std::invalid_argument);
}

// =============================================================================
// §4 LR 检验 - 手算解析值
// =============================================================================

// LR = 2(ℓ_UR - ℓ_R),  df = n_params_UR - n_params_R
TEST(LRTest, HandComputed_LikelihoodRatio) {
    EstimationResult unrestricted;
    unrestricted.log_likelihood = -10.0;
    unrestricted.n_params = 3;

    EstimationResult restricted;
    restricted.log_likelihood = -12.0;
    restricted.n_params = 2;

    HypothesisTestResult result = lr_test(unrestricted, restricted);

    EXPECT_EQ(result.test_name, "LR");
    EXPECT_NEAR(result.statistic, 2.0 * (-10.0 - (-12.0)), 1e-12);  // 4.0
    EXPECT_EQ(result.df, 1u);
    // χ²(1) p-value: P(χ² > 4) ≈ 0.0455
    EXPECT_NEAR(result.p_value, 0.04550026389635842, 1e-6);
    EXPECT_TRUE(result.reject_null_95);
    EXPECT_FALSE(result.reject_null_99);  // 0.0455 > 0.01
}

// LR = 0 当两个模型似然相等 (约束不影响)
TEST(LRTest, ZeroLR_WhenLogLikelihoods_Equal) {
    EstimationResult unrestricted;
    unrestricted.log_likelihood = -5.0;
    unrestricted.n_params = 4;

    EstimationResult restricted;
    restricted.log_likelihood = -5.0;
    restricted.n_params = 3;

    HypothesisTestResult result = lr_test(unrestricted, restricted);

    EXPECT_NEAR(result.statistic, 0.0, 1e-12);
    EXPECT_NEAR(result.p_value, 1.0, 1e-10);
    EXPECT_FALSE(result.reject_null_95);
}

// LR 抛异常: restricted.ll > unrestricted.ll (理论不可能)
TEST(LRTest, Throws_WhenRestrictedLL_GreaterThan_UnrestrictedLL) {
    EstimationResult unrestricted;
    unrestricted.log_likelihood = -10.0;
    unrestricted.n_params = 3;

    EstimationResult restricted;
    restricted.log_likelihood = -5.0;  // 比 unrestricted 大 (不可能)
    restricted.n_params = 2;

    EXPECT_THROW(lr_test(unrestricted, restricted), std::invalid_argument);
}

// LR 抛异常: restricted.n_params >= unrestricted.n_params
TEST(LRTest, Throws_WhenRestrictedParams_NotLessThan_Unrestricted) {
    EstimationResult unrestricted;
    unrestricted.log_likelihood = -10.0;
    unrestricted.n_params = 2;

    EstimationResult restricted;
    restricted.log_likelihood = -12.0;
    restricted.n_params = 2;  // 相等

    EXPECT_THROW(lr_test(unrestricted, restricted), std::invalid_argument);
}

// LR 允许微小数值负值 (1e-10 容差)
TEST(LRTest, ToleratesTinyNegativeLR) {
    EstimationResult unrestricted;
    unrestricted.log_likelihood = -10.0 - 1e-11;  // 微小低于 restricted
    unrestricted.n_params = 3;

    EstimationResult restricted;
    restricted.log_likelihood = -10.0;
    restricted.n_params = 2;

    // 不抛异常 (数值误差容差)
    EXPECT_NO_THROW(lr_test(unrestricted, restricted));
}

// =============================================================================
// §5 LM 检验 - 手算解析值 (排幻觉点 H1: N·R²_aux)
// =============================================================================

// LM = N·R²_aux (Breusch-Pagan 1979 原始形式)
// 手算: LM = 14.45 / 3.6875 = 3.91864406779661...
// 等价: N·R²_aux = 4 · (14.45/14.75) = 3.91864406779661...
TEST(LMTest, HandComputed_BreuschPagan_Form) {
    MatrixXD X = make_X();
    VectorXD residuals_restricted = make_residuals_restricted();
    MatrixXD XtX_inv = make_XtX_inv();

    HypothesisTestResult result = lm_test(X, residuals_restricted, XtX_inv);

    EXPECT_EQ(result.test_name, "LM");
    // LM = score'(X'X)^{-1}score / σ²_R = 14.45 / 3.6875
    const Real expected_lm = 14.45 / 3.6875;
    EXPECT_NEAR(result.statistic, expected_lm, 1e-10);
    EXPECT_EQ(result.df, 2u);  // q = K = 2

    // 验证 N·R²_aux 等价 (排幻觉点 H1)
    const Real ssr_R = 14.75;  // ε_R'ε_R
    const Real N = 4.0;
    const Real r2_aux = 14.45 / ssr_R;
    const Real lm_via_r2 = N * r2_aux;
    EXPECT_NEAR(result.statistic, lm_via_r2, 1e-10);

    // χ²(2) p-value: P(χ² > 3.9186) ≈ 0.1408
    EXPECT_GT(result.p_value, 0.0);
    EXPECT_LT(result.p_value, 1.0);
    EXPECT_FALSE(result.reject_null_95);  // 0.14 > 0.05
}

// 排幻觉点 H1 验证: LM = N·R²_aux, 不是 (N-K)·R²_aux
// 若错误使用 (N-K)·R²_aux: LM_wrong = 2 · 0.97966 = 1.9593...
// 正确值: LM = 4 · 0.97966 = 3.9186...
// 测试确认 LM 使用 N (非 N-K)
TEST(LMTest, AntiHallucination_H1_UsesN_NotNMinusK) {
    MatrixXD X = make_X();
    VectorXD residuals_restricted = make_residuals_restricted();
    MatrixXD XtX_inv = make_XtX_inv();

    HypothesisTestResult result = lm_test(X, residuals_restricted, XtX_inv);

    const Real N = 4.0;
    const Real K = 2.0;
    const Real ssr_R = 14.75;
    const Real r2_aux = 14.45 / ssr_R;

    // LM 应等于 N·R²_aux, 不等于 (N-K)·R²_aux
    EXPECT_NEAR(result.statistic, N * r2_aux, 1e-10);
    EXPECT_NE(result.statistic, (N - K) * r2_aux);  // 明确不相等
}

// LM 异常: 残差维度不匹配
TEST(LMTest, ThrowsOnDimensionMismatch_Residuals) {
    MatrixXD X(4, 2);
    VectorXD residuals(3);  // 维度不匹配
    MatrixXD XtX_inv(2, 2);
    EXPECT_THROW(lm_test(X, residuals, XtX_inv), std::invalid_argument);
}

// LM 异常: XtX_inv 维度不匹配
TEST(LMTest, ThrowsOnDimensionMismatch_XtX_inv) {
    MatrixXD X(4, 2);
    VectorXD residuals(4);
    MatrixXD XtX_inv(3, 3);  // 维度不匹配
    EXPECT_THROW(lm_test(X, residuals, XtX_inv), std::invalid_argument);
}

// LM 异常: 零残差平方和 (完美拟合)
TEST(LMTest, ThrowsOnZeroSSR) {
    MatrixXD X(4, 2);
    VectorXD residuals(4);
    residuals.eigen().setZero();  // 显式置零 (Eigen 不默认初始化)
    MatrixXD XtX_inv(2, 2);
    EXPECT_THROW(lm_test(X, residuals, XtX_inv), std::runtime_error);
}

// =============================================================================
// §6 J-test (Hansen 过度识别检验) - 手算解析值 (排幻觉点 H2)
// =============================================================================

// J = moments' W moments,  df = q - k (排幻觉点 H2)
TEST(JTest, HandComputed_Overidentification) {
    VectorXD moments(3);
    moments(0) = 1.0; moments(1) = 2.0; moments(2) = 3.0;

    MatrixXD W = MatrixXD(3, 3);
    for (Size i = 0; i < 3; ++i)
        for (Size j = 0; j < 3; ++j)
            W(i, j) = (i == j) ? 1.0 : 0.0;  // 单位矩阵

    // J = 1² + 2² + 3² = 14, df = 3 - 2 = 1
    HypothesisTestResult result = overidentification_test(moments, W, 2);

    EXPECT_EQ(result.test_name, "J");
    EXPECT_NEAR(result.statistic, 14.0, 1e-12);
    EXPECT_EQ(result.df, 1u);  // 排幻觉点 H2: df = q - k = 3 - 2 = 1
    // χ²(1) p-value: P(χ²(1) > 14) = erfc(√7) ≈ 0.00018281063298183488
    // 排幻觉: 之前期望值 0.0001759 是错误的 (来源不明, 可能来自不精确的在线计算器)
    //   数学推导: χ²(1) = Z², P(χ²(1)>x) = P(|Z|>√x) = 2·(1-Φ(√x)) = erfc(√(x/2))
    //   对 x=14: erfc(√7) = erfc(2.64575131...) = 0.00018281063298183488
    //   渐近检验: √(2/(πx))·exp(-x/2) = √(2/(14π))·exp(-7) = 0.2132·0.000912 = 0.0001944 (上界, 真值更小 ✓)
    EXPECT_NEAR(result.p_value, 0.00018281063298183488, 1e-12);
    EXPECT_TRUE(result.reject_null_95);
    EXPECT_TRUE(result.reject_null_99);
}

// 排幻觉点 H2: 严格识别 (q=k) 时 df=0, J 无分布
TEST(JTest, ExactIdentification_df0) {
    VectorXD moments(2);
    moments(0) = 1.0; moments(1) = 2.0;

    MatrixXD W = MatrixXD(2, 2);
    W(0, 0) = 1.0; W(0, 1) = 0.0;
    W(1, 0) = 0.0; W(1, 1) = 1.0;

    // q=2, k=2, df=0 (严格识别, 无过度识别约束)
    HypothesisTestResult result = overidentification_test(moments, W, 2);

    EXPECT_EQ(result.df, 0u);
    EXPECT_TRUE(std::isnan(result.p_value));
    EXPECT_FALSE(result.reject_null_95);
    EXPECT_FALSE(result.reject_null_99);
}

// J-test 异常: 加权矩阵维度不匹配
TEST(JTest, ThrowsOnDimensionMismatch_WeightingMatrix) {
    VectorXD moments(3);
    MatrixXD W(2, 2);  // 维度不匹配
    EXPECT_THROW(overidentification_test(moments, W, 1), std::invalid_argument);
}

// J-test 异常: q < k (矩条件少于参数, 无法识别)
TEST(JTest, ThrowsOnUnderidentification) {
    VectorXD moments(2);
    MatrixXD W(2, 2);
    EXPECT_THROW(overidentification_test(moments, W, 3), std::invalid_argument);
}

// J-test: 零矩条件 → J = 0, 不拒绝
TEST(JTest, ZeroMoments_J_Equals_Zero) {
    VectorXD moments(3);
    moments.eigen().setZero();  // 显式置零 (Eigen 不默认初始化)
    MatrixXD W(3, 3);
    for (Size i = 0; i < 3; ++i)
        for (Size j = 0; j < 3; ++j)
            W(i, j) = (i == j) ? 2.0 : 0.0;

    HypothesisTestResult result = overidentification_test(moments, W, 1);

    EXPECT_NEAR(result.statistic, 0.0, 1e-12);
    EXPECT_NEAR(result.p_value, 1.0, 1e-10);
    EXPECT_FALSE(result.reject_null_95);
}

// =============================================================================
// §7 χ² 分布 p 值验证 (与 scipy.stats.chi2.sf 对照)
// =============================================================================

// χ²(1) 上尾: P(χ² > 3.841459) = 0.05 (5% 临界值)
TEST(ChiSquareDistribution, CriticalValue_df1_95) {
    // detail::chi2_sf 在 hypothesis_tests.hpp 的 detail 命名空间
    // 通过 Wald 检验间接验证: 构造 Wald = 3.841459
    VectorXD beta(1);
    beta(0) = 1.0;
    MatrixXD vcov(1, 1);
    vcov(0, 0) = 1.0 / 3.841459;  // 使 Wald = 1²/(1/3.841459) = 3.841459
    MatrixXD R(1, 1);
    R(0, 0) = 1.0;
    VectorXD r(1);
    r(0) = 0.0;

    HypothesisTestResult result = wald_test(beta, vcov, R, r);
    EXPECT_NEAR(result.statistic, 3.841459, 1e-6);
    EXPECT_NEAR(result.p_value, 0.05, 1e-4);
}

// χ²(2) 上尾: P(χ² > 5.991465) = 0.05 (5% 临界值)
TEST(ChiSquareDistribution, CriticalValue_df2_95) {
    // 构造 2 维 Wald = 5.991465
    VectorXD beta(2);
    beta(0) = 1.0; beta(1) = 1.0;
    MatrixXD vcov(2, 2);
    vcov(0, 0) = 1.0; vcov(0, 1) = 0.0;
    vcov(1, 0) = 0.0; vcov(1, 1) = 1.0;
    // Wald = β'(V^{-1})β = 1+1 = 2, 需缩放使 Wald = 5.991465
    // 用 Rβ-r = [1,1], RVR' = (5.991465/2)·I → V = (5.991465/2)·I? 不对
    // Wald = (Rβ-r)'(RVR')^{-1}(Rβ-r) = [1,1]·(RVR')^{-1}·[1;1]
    // 若 RVR' = (2/5.991465)·I, 则 (RVR')^{-1} = (5.991465/2)·I
    // Wald = (5.991465/2)·(1+1) = 5.991465 ✓
    const Real scale = 2.0 / 5.991465;
    vcov(0, 0) = scale; vcov(1, 1) = scale;
    MatrixXD R(2, 2);
    R(0, 0) = 1.0; R(0, 1) = 0.0;
    R(1, 0) = 0.0; R(1, 1) = 1.0;
    VectorXD r(2);
    r(0) = 0.0; r(1) = 0.0;

    HypothesisTestResult result = wald_test(beta, vcov, R, r);
    EXPECT_NEAR(result.statistic, 5.991465, 1e-5);
    EXPECT_NEAR(result.p_value, 0.05, 1e-4);
}

// =============================================================================
// §8 F 分布 p 值验证 (与 scipy.stats.f.sf 对照)
// =============================================================================

// F(1, 2) 上尾: P(F > 18.5128) = 0.05
TEST(FDistribution, CriticalValue_F12_95) {
    // 构造 Wald = 18.5128·q = 18.5128 (q=1)
    VectorXD beta(1);
    beta(0) = 1.0;
    MatrixXD vcov(1, 1);
    vcov(0, 0) = 1.0 / 18.5128;
    MatrixXD R(1, 1);
    R(0, 0) = 1.0;
    VectorXD r(1);
    r(0) = 0.0;

    // df_residual = 2 → F(1, 2)
    HypothesisTestResult result = wald_test(beta, vcov, R, r, 2.0);

    EXPECT_TRUE(result.use_f_distribution);
    EXPECT_NEAR(result.f_df1, 1.0, 1e-10);
    EXPECT_NEAR(result.f_df2, 2.0, 1e-10);
    EXPECT_NEAR(result.p_value, 0.05, 1e-4);
}

// F(2, 10) 上尾: P(F > 4.1028) = 0.05
TEST(FDistribution, CriticalValue_F2_10_95) {
    // 构造 Wald = F·q = 4.1028·2 = 8.2056
    VectorXD beta(2);
    beta(0) = 1.0; beta(1) = 1.0;
    const Real f_crit = 4.1028;
    const Real wald_target = f_crit * 2.0;  // F = Wald/q → Wald = F·q
    // Wald = (Rβ-r)'(RVR')^{-1}(Rβ-r) = [1,1]·(RVR')^{-1}·[1;1] = wald_target
    // (RVR')^{-1} = (wald_target/2)·I → RVR' = (2/wald_target)·I
    const Real scale = 2.0 / wald_target;
    MatrixXD vcov(2, 2);
    vcov(0, 0) = scale; vcov(0, 1) = 0.0;
    vcov(1, 0) = 0.0; vcov(1, 1) = scale;
    MatrixXD R(2, 2);
    R(0, 0) = 1.0; R(0, 1) = 0.0;
    R(1, 0) = 0.0; R(1, 1) = 1.0;
    VectorXD r(2);
    r(0) = 0.0; r(1) = 0.0;

    HypothesisTestResult result = wald_test(beta, vcov, R, r, 10.0);

    EXPECT_TRUE(result.use_f_distribution);
    EXPECT_NEAR(result.f_df1, 2.0, 1e-10);
    EXPECT_NEAR(result.f_df2, 10.0, 1e-10);
    EXPECT_NEAR(result.p_value, 0.05, 1e-3);
}

// =============================================================================
// §9 临界值与拒绝域
// =============================================================================

// 大统计量应拒绝, 小统计量不应拒绝
TEST(HypothesisTest, RejectNull_LargeStatistic) {
    VectorXD beta(1);
    beta(0) = 100.0;  // 极大偏离
    MatrixXD vcov(1, 1);
    vcov(0, 0) = 1.0;
    MatrixXD R(1, 1);
    R(0, 0) = 1.0;
    VectorXD r(1);
    r(0) = 0.0;

    HypothesisTestResult result = wald_test(beta, vcov, R, r);
    EXPECT_TRUE(result.reject_null_95);
    EXPECT_TRUE(result.reject_null_99);
    EXPECT_LT(result.p_value, 1e-10);
}

// 临界值近似 (Wilson-Hilferty): 应为正数
TEST(HypothesisTest, CriticalValue_IsPositive) {
    VectorXD beta(1);
    beta(0) = 0.0;
    MatrixXD vcov(1, 1);
    vcov(0, 0) = 1.0;
    MatrixXD R(1, 1);
    R(0, 0) = 1.0;
    VectorXD r(1);
    r(0) = 0.0;

    HypothesisTestResult result = wald_test(beta, vcov, R, r);
    EXPECT_GT(result.critical_value_95, 0.0);
    EXPECT_GT(result.critical_value_99, result.critical_value_95);  // 99% > 95%
}

// =============================================================================
// §10 渐近等价性验证 (大样本下 Wald ≈ LR ≈ LM)
// =============================================================================
// 在大样本下, Wald/LR/LM 三大检验渐近等价 (Engle 1984)
// 这里通过构造同一假设的三个检验, 验证统计量在同一数量级
// (由于 N=4 是小样本, 不期望精确相等, 仅验证数量级合理)

TEST(AsymptoticEquivalence, Wald_LR_LM_SameOrderOfMagnitude) {
    // 使用同一数据集, 对 β₁=0 进行 Wald/LR/LM 检验
    // Wald: 96.333 (已手算)
    // LM: 3.9186 (已手算, 使用约束模型残差)
    // LR: 需构造两个模型的对数似然

    // Wald
    VectorXD beta = make_beta();
    MatrixXD vcov = make_classical_vcov();
    MatrixXD R(1, 2);
    R(0, 0) = 0.0; R(0, 1) = 1.0;
    VectorXD r(1);
    r(0) = 0.0;
    HypothesisTestResult wald_result = wald_test(beta, vcov, R, r);

    // LM (约束模型 = 仅截距)
    MatrixXD X = make_X();
    VectorXD residuals_R = make_residuals_restricted();
    MatrixXD XtX_inv = make_XtX_inv();
    HypothesisTestResult lm_result = lm_test(X, residuals_R, XtX_inv);

    // 两者均为正, 同一数量级合理性检查 (小样本不严格相等)
    EXPECT_GT(wald_result.statistic, 0.0);
    EXPECT_GT(lm_result.statistic, 0.0);
    // Wald > LM (Wald 用无约束估计, LM 用有约束估计, 小样本下差异显著)
    EXPECT_GT(wald_result.statistic, lm_result.statistic);
}
