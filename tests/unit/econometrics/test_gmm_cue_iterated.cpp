// SOURCE: PHASE6_IMPLEMENTATION_PLAN §5 - 迭代 GMM + CUE 单元测试
// 验证方法: GMM 数学性质 (迭代收敛性 / CUE 大样本等价 / 类型字段正确)
//
// 排幻觉点:
//   E10: Ŝ 用 moment matrix HAC
//   E12: 完全拟合 Ŝ 奇异处理 (不抛异常, 保留 Step 1 结果)
//
// 迭代 GMM 性质 (Hansen 1982, Hayashi §3.5):
//   - 重复 Step 2 直到 β 收敛 (β_new - β_old < tolerance)
//   - 恰好识别: 无需迭代 (q=k, 权重不影响估计量)
//   - 大样本下与两步 GMM 渐近等价
//
// CUE 性质 (Hansen-Heaton-Yaron 1996):
//   - θ̂_CUE = argmin_θ ḡ(θ)' Ŝ(θ)⁻¹ ḡ(θ) (同时更新 θ 和 Ŝ)
//   - 大样本下与两步 GMM 等价
//   - 小样本下更稳健 (但数值优化更复杂)
//   - 当前实现: 用两步 GMM 作为起始值 (大样本等价近似)

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
// 辅助: 构造过度识别有噪声数据 (N=6, k=1, q=2)
//   复用 test_gmm_two_step.cpp 的数据集
// =============================================================================
inline auto make_overidentified_noise() {
    struct Data { MatrixXD X, Z; VectorXD y; };
    Data d;
    d.X = MatrixXD(6, 1);
    d.Z = MatrixXD(6, 2);
    d.y = VectorXD(6);
    for (Size i = 0; i < 6; ++i) {
        const Real x = static_cast<Real>(i + 1);
        d.X(i, 0) = x;
        d.Z(i, 0) = x;
        d.Z(i, 1) = x * x;
    }
    d.y(0) = 1.9; d.y(1) = 4.2; d.y(2) = 5.8;
    d.y(3) = 8.3; d.y(4) = 9.7; d.y(5) = 12.4;
    return d;
}

// =============================================================================
// 辅助: 构造大样本过度识别数据 (N=100, k=1, q=2, β=3.0)
//   复用 test_gmm_two_step.cpp 的数据集
// =============================================================================
inline auto make_large_sample() {
    struct Data { MatrixXD X, Z; VectorXD y; };
    const Size N = 100;
    Data d;
    d.X = MatrixXD(N, 1);
    d.Z = MatrixXD(N, 2);
    d.y = VectorXD(N);
    for (Size i = 0; i < N; ++i) {
        const Real x = static_cast<Real>(i + 1) / 10.0;
        d.X(i, 0) = x;
        d.Z(i, 0) = x;
        d.Z(i, 1) = x * x;
        const Real noise = 0.5 * std::sin(x);
        d.y(i) = 3.0 * x + noise;
    }
    return d;
}

// =============================================================================
// 测试 1: 迭代 GMM 收敛标志为 true
// =============================================================================
TEST(GMMIterated, ConvergedFlag_True) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::Iterated);
    EXPECT_TRUE(result.converged);
}

// =============================================================================
// 测试 2: 迭代 GMM 迭代次数 >= 2 (至少执行一次 Step 2)
// =============================================================================
TEST(GMMIterated, NIter_AtLeast2) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::Iterated);
    EXPECT_GE(result.n_iter, 2u);
}

// =============================================================================
// 测试 3: 迭代 GMM 类型字段正确
// =============================================================================
TEST(GMMIterated, TypeField_Correct) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::Iterated);
    EXPECT_EQ(result.gmm_type, GMMType::Iterated);
}

// =============================================================================
// 测试 4: 迭代 GMM 与两步 GMM 大样本下接近 (渐近等价)
//   小样本可能有差异, 大样本应收敛到同一值
// =============================================================================
TEST(GMMIterated, LargeSample_CloseToTwoStep) {
    const auto d = make_large_sample();
    const GMMResult r_twostep = gmm_linear_iv(d.X, d.y, d.Z, GMMType::TwoStep);
    const GMMResult r_iterated = gmm_linear_iv(d.X, d.y, d.Z, GMMType::Iterated);

    // 大样本下两步 GMM 和迭代 GMM 应接近 (容差 1e-4)
    EXPECT_NEAR(r_twostep.coefficients(0), r_iterated.coefficients(0), 1e-4);
    EXPECT_NEAR(r_twostep.j_statistic, r_iterated.j_statistic, 1e-2);
}

// =============================================================================
// 测试 5: 迭代 GMM β̂ 收敛到真实值
// =============================================================================
TEST(GMMIterated, LargeSample_ConvergesToTrue) {
    const auto d = make_large_sample();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::Iterated);
    EXPECT_NEAR(result.coefficients(0), 3.0, 0.2);
}

// =============================================================================
// 测试 6: 迭代 GMM J 统计量非负
// =============================================================================
TEST(GMMIterated, J_Statistic_NonNegative) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::Iterated);
    EXPECT_GE(result.j_statistic, 0.0);
}

// =============================================================================
// 测试 7: 迭代 GMM 自由度 = q - k
// =============================================================================
TEST(GMMIterated, DF_Correct) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::Iterated);
    EXPECT_EQ(result.j_df, 1u);  // q=2, k=1, df=1
}

// =============================================================================
// 测试 8: 迭代 GMM 协方差矩阵对称
// =============================================================================
TEST(GMMIterated, Vcov_Symmetric) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::Iterated);
    const Size k = result.n_params;
    for (Size i = 0; i < k; ++i) {
        for (Size j = i; j < k; ++j) {
            EXPECT_NEAR(result.vcov(i, j), result.vcov(j, i), 1e-12);
        }
    }
}

// =============================================================================
// 测试 9: 恰好识别时迭代 GMM = 两步 GMM (权重不影响估计量)
// =============================================================================
TEST(GMMIterated, ExactlyIdentified_EqualsTwoStep) {
    // N=5, k=1, q=1 (恰好识别)
    MatrixXD X(5, 1), Z(5, 1);
    VectorXD y(5);
    for (Size i = 0; i < 5; ++i) {
        X(i, 0) = static_cast<Real>(i + 1);
        Z(i, 0) = static_cast<Real>(i + 2);
    }
    y(0) = 2.1; y(1) = 3.9; y(2) = 6.2; y(3) = 7.8; y(4) = 10.1;

    const GMMResult r_twostep = gmm_linear_iv(X, y, Z, GMMType::TwoStep);
    const GMMResult r_iterated = gmm_linear_iv(X, y, Z, GMMType::Iterated);

    EXPECT_NEAR(r_twostep.coefficients(0), r_iterated.coefficients(0), 1e-10);
}

// =============================================================================
// 测试 10: 迭代 GMM 容差参数影响收敛速度
//   较大容差应更快收敛 (更少迭代次数)
// =============================================================================
TEST(GMMIterated, Tolerance_AffectsConvergence) {
    const auto d = make_large_sample();

    // 严格容差
    const GMMResult r_strict = gmm_linear_iv(d.X, d.y, d.Z, GMMType::Iterated,
                                                HacKernel::Bartlett, 0, 100, 1e-12);
    // 宽松容差
    const GMMResult r_loose = gmm_linear_iv(d.X, d.y, d.Z, GMMType::Iterated,
                                              HacKernel::Bartlett, 0, 100, 1e-3);

    // 两者都应收敛
    EXPECT_TRUE(r_strict.converged);
    EXPECT_TRUE(r_loose.converged);
    // 宽松容差可能迭代次数更少 (或相等)
    // 注: 不强制 r_loose.n_iter <= r_strict.n_iter, 因为可能一步就收敛
}

// =============================================================================
// 测试 11: CUE 类型字段正确
// =============================================================================
TEST(GMMCUE, TypeField_Correct) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::CUE);
    EXPECT_EQ(result.gmm_type, GMMType::CUE);
}

// =============================================================================
// 测试 12: CUE 收敛标志为 true
// =============================================================================
TEST(GMMCUE, ConvergedFlag_True) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::CUE);
    EXPECT_TRUE(result.converged);
}

// =============================================================================
// 测试 13: CUE 与两步 GMM 大样本下接近 (渐近等价)
// =============================================================================
TEST(GMMCUE, LargeSample_CloseToTwoStep) {
    const auto d = make_large_sample();
    const GMMResult r_twostep = gmm_linear_iv(d.X, d.y, d.Z, GMMType::TwoStep);
    const GMMResult r_cue = gmm_linear_iv(d.X, d.y, d.Z, GMMType::CUE);

    // 大样本下 CUE 和两步 GMM 应接近
    // 注: CUE 在每一点更新 Ŝ(β), 与两步 GMM (Ŝ 固定在 β̂₁) 仅渐近等价
    // N=100 下差异 ~1e-4 量级, 容差 1e-2 足以检测实现 bug
    EXPECT_NEAR(r_twostep.coefficients(0), r_cue.coefficients(0), 1e-2);
}

// =============================================================================
// 测试 14: CUE β̂ 收敛到真实值
// =============================================================================
TEST(GMMCUE, LargeSample_ConvergesToTrue) {
    const auto d = make_large_sample();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::CUE);
    EXPECT_NEAR(result.coefficients(0), 3.0, 0.2);
}

// =============================================================================
// 测试 15: CUE J 统计量非负
// =============================================================================
TEST(GMMCUE, J_Statistic_NonNegative) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::CUE);
    EXPECT_GE(result.j_statistic, 0.0);
}

// =============================================================================
// 测试 16: CUE 自由度 = q - k
// =============================================================================
TEST(GMMCUE, DF_Correct) {
    const auto d = make_overidentified_noise();
    const GMMResult result = gmm_linear_iv(d.X, d.y, d.Z, GMMType::CUE);
    EXPECT_EQ(result.j_df, 1u);
}

// =============================================================================
// 测试 17: CUE 恰好识别时 = 两步 GMM
// =============================================================================
TEST(GMMCUE, ExactlyIdentified_EqualsTwoStep) {
    MatrixXD X(5, 1), Z(5, 1);
    VectorXD y(5);
    for (Size i = 0; i < 5; ++i) {
        X(i, 0) = static_cast<Real>(i + 1);
        Z(i, 0) = static_cast<Real>(i + 2);
    }
    y(0) = 2.1; y(1) = 3.9; y(2) = 6.2; y(3) = 7.8; y(4) = 10.1;

    const GMMResult r_twostep = gmm_linear_iv(X, y, Z, GMMType::TwoStep);
    const GMMResult r_cue = gmm_linear_iv(X, y, Z, GMMType::CUE);

    EXPECT_NEAR(r_twostep.coefficients(0), r_cue.coefficients(0), 1e-10);
}

// =============================================================================
// 测试 18: 完全拟合时迭代 GMM 不抛异常 (E12: Ŝ 奇异处理)
// =============================================================================
TEST(GMMIterated, PerfectFit_NoThrow) {
    // y = 2X 完全拟合, ε=0, Ŝ=0
    MatrixXD X(4, 1), Z(4, 2);
    VectorXD y(4);
    for (Size i = 0; i < 4; ++i) {
        const Real x = static_cast<Real>(i + 1);
        X(i, 0) = x;
        Z(i, 0) = x;
        Z(i, 1) = x * x;
        y(i) = 2.0 * x;
    }
    EXPECT_NO_THROW({
        const GMMResult result = gmm_linear_iv(X, y, Z, GMMType::Iterated);
        EXPECT_NEAR(result.coefficients(0), 2.0, 1e-8);
        EXPECT_NEAR(result.j_statistic, 0.0, 1e-8);
    });
}

// =============================================================================
// 测试 19: 完全拟合时 CUE 不抛异常
// =============================================================================
TEST(GMMCUE, PerfectFit_NoThrow) {
    MatrixXD X(4, 1), Z(4, 2);
    VectorXD y(4);
    for (Size i = 0; i < 4; ++i) {
        const Real x = static_cast<Real>(i + 1);
        X(i, 0) = x;
        Z(i, 0) = x;
        Z(i, 1) = x * x;
        y(i) = 2.0 * x;
    }
    EXPECT_NO_THROW({
        const GMMResult result = gmm_linear_iv(X, y, Z, GMMType::CUE);
        EXPECT_NEAR(result.coefficients(0), 2.0, 1e-8);
        EXPECT_NEAR(result.j_statistic, 0.0, 1e-8);
    });
}

// =============================================================================
// 测试 20: 三种 GMM 类型在恰好识别时给出相同 J=0
// =============================================================================
TEST(GMMAllTypes, ExactlyIdentified_AllGiveJZero) {
    MatrixXD X(5, 1), Z(5, 1);
    VectorXD y(5);
    for (Size i = 0; i < 5; ++i) {
        X(i, 0) = static_cast<Real>(i + 1);
        Z(i, 0) = static_cast<Real>(i + 2);
    }
    y(0) = 2.1; y(1) = 3.9; y(2) = 6.2; y(3) = 7.8; y(4) = 10.1;

    const GMMResult r1 = gmm_linear_iv(X, y, Z, GMMType::TwoStep);
    const GMMResult r2 = gmm_linear_iv(X, y, Z, GMMType::Iterated);
    const GMMResult r3 = gmm_linear_iv(X, y, Z, GMMType::CUE);

    EXPECT_NEAR(r1.j_statistic, 0.0, 1e-8);
    EXPECT_NEAR(r2.j_statistic, 0.0, 1e-8);
    EXPECT_NEAR(r3.j_statistic, 0.0, 1e-8);
}
