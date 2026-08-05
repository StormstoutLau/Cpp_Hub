// SOURCE: PHASE6_IMPLEMENTATION_PLAN §5 - GMM 两步估计单元测试
// 验证方法: 手算解析值 + GMM 数学性质 (恰好识别=2SLS, 完全拟合 J=0, J≥0)
//
// 排幻觉点:
//   E10: Ŝ 用 moment matrix HAC (Z' diag(ε²) Z), 非 tangent matrix (数值导数)
//   E11: 工具变量构造 (Arellano-Bond 在 test_arellano_bond.cpp 验证)
//   H2:  J-test df = q - k (过度识别约束数), 非 q
//
// 手算数据集 1 (恰好识别, N=5, k=1, q=1):
//   X = [1, 2, 3, 4, 5]',  Z = [2, 3, 4, 5, 6]'  (Z = X+1)
//   y = [2.1, 3.9, 6.2, 7.8, 10.1]'  (y ≈ 2X + noise)
//   恰好识别: β̂_GMM = β̂_2SLS = Z'y / Z'X (k=q=1 标量情形)
//     Z'X = 2+6+12+20+30 = 70
//     Z'y = 2*2.1 + 3*3.9 + 4*6.2 + 5*7.8 + 6*10.1
//         = 4.2 + 11.7 + 24.8 + 39.0 + 60.6 = 140.3
//     β̂ = 140.3 / 70 = 2.0042857142857144...
//   J = 0 (df = q-k = 0, 恰好识别无过度识别约束)
//
// 手算数据集 2 (过度识别完全拟合, N=4, k=1, q=2):
//   X = [1, 2, 3, 4]',  Z = [[1,1],[2,4],[3,9],[4,16]]  (Z1=X, Z2=X²)
//   y = [2, 4, 6, 8]'  (y = 2X, 无误差)
//   β̂ = 2.0 (完全拟合), J = 0 (ε=0 → g=0)

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cpphub/econometrics/estimation/gmm.hpp"
#include "cpphub/econometrics/inference/hac_kernels.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;
using cpphub::v1::linalg::dynamic::MatrixXD;
using cpphub::v1::linalg::dynamic::VectorXD;

// =============================================================================
// 辅助: 构造恰好识别数据集 (N=5, k=1, q=1)
//   X = [1,2,3,4,5], Z = [2,3,4,5,6], y = [2.1, 3.9, 6.2, 7.8, 10.1]
//   β̂_2SLS = Z'y / Z'X = 140.3 / 70 = 2.0042857142857144
// =============================================================================
struct ExactlyIdentifiedData {
    MatrixXD X, Z;
    VectorXD y;
    Real expected_beta;
};

inline ExactlyIdentifiedData make_exactly_identified() {
    ExactlyIdentifiedData d;
    d.X = MatrixXD(5, 1);
    d.Z = MatrixXD(5, 1);
    d.y = VectorXD(5);
    // X = [1, 2, 3, 4, 5]
    for (Size i = 0; i < 5; ++i) d.X(i, 0) = static_cast<Real>(i + 1);
    // Z = [2, 3, 4, 5, 6]  (Z = X + 1)
    for (Size i = 0; i < 5; ++i) d.Z(i, 0) = static_cast<Real>(i + 2);
    // y = [2.1, 3.9, 6.2, 7.8, 10.1]
    d.y(0) = 2.1; d.y(1) = 3.9; d.y(2) = 6.2; d.y(3) = 7.8; d.y(4) = 10.1;
    // β̂ = Z'y / Z'X = 140.3 / 70
    d.expected_beta = 140.3 / 70.0;
    return d;
}

// =============================================================================
// 辅助: 构造过度识别完全拟合数据集 (N=4, k=1, q=2)
//   X = [1,2,3,4], Z = [[1,1],[2,4],[3,9],[4,16]], y = [2,4,6,8]
//   β̂ = 2.0 (y = 2X), J = 0 (ε=0)
// =============================================================================
struct OveridentifiedPerfectData {
    MatrixXD X, Z;
    VectorXD y;
    Real expected_beta;
};

inline OveridentifiedPerfectData make_overidentified_perfect() {
    OveridentifiedPerfectData d;
    d.X = MatrixXD(4, 1);
    d.Z = MatrixXD(4, 2);
    d.y = VectorXD(4);
    // X = [1, 2, 3, 4]
    for (Size i = 0; i < 4; ++i) d.X(i, 0) = static_cast<Real>(i + 1);
    // Z = [[1,1], [2,4], [3,9], [4,16]]  (Z1=X, Z2=X²)
    for (Size i = 0; i < 4; ++i) {
        d.Z(i, 0) = static_cast<Real>(i + 1);
        d.Z(i, 1) = static_cast<Real>((i + 1) * (i + 1));
    }
    // y = [2, 4, 6, 8]  (y = 2X)
    d.y(0) = 2.0; d.y(1) = 4.0; d.y(2) = 6.0; d.y(3) = 8.0;
    d.expected_beta = 2.0;
    return d;
}

// =============================================================================
// 辅助: 构造过度识别有噪声数据集 (N=6, k=1, q=2)
//   X = [1,2,3,4,5,6], Z = [[1,1],[2,4],[3,9],[4,16],[5,25],[6,36]]
//   y = [1.9, 4.2, 5.8, 8.3, 9.7, 12.4]  (y ≈ 2X + noise)
//   β̂ ≈ 2.0, J > 0 (有过度识别约束, 残差≠0)
// =============================================================================
struct OveridentifiedNoiseData {
    MatrixXD X, Z;
    VectorXD y;
};

inline OveridentifiedNoiseData make_overidentified_noise() {
    OveridentifiedNoiseData d;
    d.X = MatrixXD(6, 1);
    d.Z = MatrixXD(6, 2);
    d.y = VectorXD(6);
    for (Size i = 0; i < 6; ++i) {
        const Real x = static_cast<Real>(i + 1);
        d.X(i, 0) = x;
        d.Z(i, 0) = x;
        d.Z(i, 1) = x * x;
    }
    // y ≈ 2X + [-0.1, 0.2, -0.2, 0.3, -0.3, 0.4]
    d.y(0) = 1.9; d.y(1) = 4.2; d.y(2) = 5.8;
    d.y(3) = 8.3; d.y(4) = 9.7; d.y(5) = 12.4;
    return d;
}

// =============================================================================
// 测试 1: 恰好识别 - β̂ 等于 2SLS, J=0 (df=0)
// =============================================================================
TEST(GMMTwoStep, ExactlyIdentified_Equals2SLS) {
    const auto d = make_exactly_identified();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::TwoStep);

    // β̂ = 140.3/70 = 2.0042857142857144 (手算)
    EXPECT_NEAR(result.coefficients(0), d.expected_beta, 1e-10);

    // 恰好识别: df = q - k = 0, J 无意义 (NaN 或 0)
    EXPECT_EQ(result.j_df, 0u);
    // J 统计量: 恰好识别时数学上 J=0 (矩条件精确为 0)
    // 注: df=0 时 p 值为 NaN, 不检验 p_value
    EXPECT_NEAR(result.j_statistic, 0.0, 1e-8);
}

// =============================================================================
// 测试 2: 恰好识别 - 维度字段正确
// =============================================================================
TEST(GMMTwoStep, ExactlyIdentified_Dimensions) {
    const auto d = make_exactly_identified();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z);

    EXPECT_EQ(result.n_obs, 5u);
    EXPECT_EQ(result.n_params, 1u);
    EXPECT_EQ(result.n_moments, 1u);
    EXPECT_EQ(result.gmm_type, GMMType::TwoStep);
    EXPECT_TRUE(result.converged);
}

// =============================================================================
// 测试 3: 恰好识别 - 标准误为正, t 统计量 = β̂/SE
// =============================================================================
TEST(GMMTwoStep, ExactlyIdentified_StdErrors) {
    const auto d = make_exactly_identified();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z);

    const Real se = result.std_errors(0);
    EXPECT_GT(se, 0.0);
    EXPECT_NEAR(result.t_statistics(0), result.coefficients(0) / se, 1e-10);
    // 双侧 p 值: erfc(|t|/√2)
    const Real expected_p = std::erfc(std::fabs(result.t_statistics(0)) / std::sqrt(2.0));
    EXPECT_NEAR(result.p_values(0), expected_p, 1e-10);
    EXPECT_GE(result.p_values(0), 0.0);
    EXPECT_LE(result.p_values(0), 1.0);
}

// =============================================================================
// 测试 4: 过度识别完全拟合 - β̂=2.0, J=0 (ε=0)
// =============================================================================
TEST(GMMTwoStep, Overidentified_PerfectFit_JZero) {
    const auto d = make_overidentified_perfect();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::TwoStep);

    // β̂ = 2.0 (完全拟合 y=2X)
    EXPECT_NEAR(result.coefficients(0), d.expected_beta, 1e-10);

    // df = q - k = 2 - 1 = 1
    EXPECT_EQ(result.j_df, 1u);

    // J = 0 (残差=0 → 矩条件=0)
    EXPECT_NEAR(result.j_statistic, 0.0, 1e-8);
    // p 值 = 1.0 (J=0 → 不拒绝 H0)
    EXPECT_NEAR(result.j_pvalue, 1.0, 1e-8);
}

// =============================================================================
// 测试 5: 过度识别有噪声 - J > 0, df = 1
// =============================================================================
TEST(GMMTwoStep, Overidentified_Noise_JPositive) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::TwoStep);

    // β̂ ≈ 2.0 (真实值)
    EXPECT_NEAR(result.coefficients(0), 2.0, 0.1);

    // df = q - k = 2 - 1 = 1
    EXPECT_EQ(result.j_df, 1u);

    // J > 0 (有过度识别约束, 残差≠0)
    EXPECT_GT(result.j_statistic, 0.0);

    // p 值在 [0, 1] 区间
    EXPECT_GE(result.j_pvalue, 0.0);
    EXPECT_LE(result.j_pvalue, 1.0);
}

// =============================================================================
// 测试 6: J 统计量非负 (二次型性质)
// =============================================================================
TEST(GMMTwoStep, J_Statistic_NonNegative) {
    const auto d1 = make_exactly_identified();
    const auto d2 = make_overidentified_perfect();
    const auto d3 = make_overidentified_noise();

    const GMMResult r1 = gmm_linear_iv(d1.X, d1.y, d1.Z);
    const GMMResult r2 = gmm_linear_iv(d2.X, d2.y, d2.Z);
    const GMMResult r3 = gmm_linear_iv(d3.X, d3.y, d3.Z);

    EXPECT_GE(r1.j_statistic, 0.0);
    EXPECT_GE(r2.j_statistic, 0.0);
    EXPECT_GE(r3.j_statistic, 0.0);
}

// =============================================================================
// 测试 7: 欠识别 (q < k) 抛异常
// =============================================================================
TEST(GMMTwoStep, Underidentified_Throws) {
    // k=2, q=1 (q < k)
    MatrixXD X(5, 2);
    MatrixXD Z(5, 1);
    VectorXD y(5);
    for (Size i = 0; i < 5; ++i) {
        X(i, 0) = 1.0;
        X(i, 1) = static_cast<Real>(i + 1);
        Z(i, 0) = static_cast<Real>(i + 2);
        y(i) = static_cast<Real>(2 * (i + 1));
    }
    EXPECT_THROW(gmm_linear_iv(X, y, Z), std::invalid_argument);
}

// =============================================================================
// 测试 8: 维度不匹配抛异常
// =============================================================================
TEST(GMMTwoStep, DimensionMismatch_Throws) {
    MatrixXD X(5, 1);
    VectorXD y(4);  // y.size() != X.rows()
    MatrixXD Z(5, 1);
    EXPECT_THROW(gmm_linear_iv(X, y, Z), std::invalid_argument);

    VectorXD y2(5);
    MatrixXD Z2(4, 1);  // Z.rows() != X.rows()
    EXPECT_THROW(gmm_linear_iv(X, y2, Z2), std::invalid_argument);
}

// =============================================================================
// 测试 9: 协方差矩阵对称且对角线为正
// =============================================================================
TEST(GMMTwoStep, Vcov_Symmetric_PositiveDiagonal) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z);

    const Size k = result.n_params;
    for (Size i = 0; i < k; ++i) {
        EXPECT_GT(result.vcov(i, i), 0.0);  // 对角线为正
        for (Size j = i; j < k; ++j) {
            EXPECT_NEAR(result.vcov(i, j), result.vcov(j, i), 1e-12);  // 对称
        }
    }
}

// =============================================================================
// 测试 10: 目标函数值 = J 统计量
// =============================================================================
TEST(GMMTwoStep, ObjectiveEquals_JStatistic) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z);

    // GMM 目标函数 J = N · g' S^{-1} g = j_statistic
    EXPECT_NEAR(result.objective_value, result.j_statistic, 1e-10);
}

// =============================================================================
// 测试 11: 两步 GMM 第一步 = 2SLS (恰好识别时两步结果一致)
//   恰好识别 (q=k): 权重矩阵不影响估计量, 两步 GMM = 一步 GMM = 2SLS
// =============================================================================
TEST(GMMTwoStep, ExactlyIdentified_TwoStepEqualsOneStep) {
    const auto d = make_exactly_identified();

    // 两步 GMM
    const GMMResult r_twostep = gmm_linear_iv(d.X, d.y, d.Z, GMMType::TwoStep);
    // 迭代 GMM (恰好识别时也等于 2SLS)
    const GMMResult r_iterated = gmm_linear_iv(d.X, d.y, d.Z, GMMType::Iterated);

    // 恰好识别: 所有 GMM 类型给出相同 β̂
    EXPECT_NEAR(r_twostep.coefficients(0), r_iterated.coefficients(0), 1e-10);
    EXPECT_NEAR(r_twostep.coefficients(0), d.expected_beta, 1e-10);
}

// =============================================================================
// 测试 12: 迭代次数 = 2 (两步 GMM)
// =============================================================================
TEST(GMMTwoStep, TwoStep_NIter_Is2) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::TwoStep);
    EXPECT_EQ(result.n_iter, 2u);
}

// =============================================================================
// 测试 13: Hansen J 检验 df=1 时 p 值 = erfc(√(J/2))
//   排幻觉点: χ²(1) 的 SF = erfc(√(x/2)), 非 gammq 近似
// =============================================================================
TEST(GMMTwoStep, JTest_DF1_PValue_ErfcForm) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z);

    ASSERT_EQ(result.j_df, 1u);
    // df=1: P(χ²(1) > J) = erfc(√(J/2))
    const Real expected_p = std::erfc(std::sqrt(0.5 * result.j_statistic));
    EXPECT_NEAR(result.j_pvalue, expected_p, 1e-10);
}

// =============================================================================
// 测试 14: 大样本过度识别 - β̂ 收敛到真实值
//   构造 N=100 数据, β=3.0, 验证估计值接近 3.0
// =============================================================================
TEST(GMMTwoStep, LargeSample_ConvergesToTrue) {
    // 构造 N=100, k=1, q=2 的数据
    // y = 3*X + ε, Z = [X, X²] (X 外生, Z 与 ε 不相关)
    const Size N = 100;
    MatrixXD X(N, 1);
    MatrixXD Z(N, 2);
    VectorXD y(N);

    // 用固定模式生成确定性数据 (避免随机种子跨平台问题)
    for (Size i = 0; i < N; ++i) {
        const Real x = static_cast<Real>(i + 1) / 10.0;  // x ∈ [0.1, 10.0]
        X(i, 0) = x;
        Z(i, 0) = x;
        Z(i, 1) = x * x;
        // 确定性 "噪声": sin(x)*0.5 (均值≈0, 与 Z 弱相关)
        const Real noise = 0.5 * std::sin(x);
        y(i) = 3.0 * x + noise;
    }

    const GMMResult result = gmm_linear_iv(X, y, Z, GMMType::TwoStep);

    // β̂ ≈ 3.0 (噪声均值≈0, 大样本收敛)
    EXPECT_NEAR(result.coefficients(0), 3.0, 0.2);
    EXPECT_EQ(result.n_obs, N);
    EXPECT_EQ(result.n_moments, 2u);
    EXPECT_EQ(result.n_params, 1u);
}

// =============================================================================
// 测试 15: 弱工具变量检测 (Z'Z 奇异时抛异常)
// =============================================================================
TEST(GMMTwoStep, WeakInstruments_Throws) {
    // Z 全零列 → Z'Z 奇异
    MatrixXD X(5, 1);
    MatrixXD Z(5, 2);
    VectorXD y(5);
    for (Size i = 0; i < 5; ++i) {
        X(i, 0) = static_cast<Real>(i + 1);
        Z(i, 0) = 0.0;  // 全零
        Z(i, 1) = 0.0;  // 全零
        y(i) = static_cast<Real>(i + 1);
    }
    EXPECT_THROW(gmm_linear_iv(X, y, Z), std::runtime_error);
}
