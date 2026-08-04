// SOURCE: PHASE6_IMPLEMENTATION_PLAN §2.3 - 双层 linalg 动态尺寸层 (ADR-013)
// 基准来源: Eigen 3.4.0 原生 API (https://eigen.tuxfamily.org/dox/group__DenseDecompositionBenchmark.html)
// 容差: 1e-14 (机器精度量级, 对照 Eigen 原生结果)
// 排幻觉点: 审计验收文档 §4.1 E1-E7 (Eigen3 版本/许可/CMake/隔离/LLT 稳定性)
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include "cpphub/core/linalg_dynamic.hpp"

using namespace cpphub::v1::linalg::dynamic;

namespace {
constexpr double TOL = 1e-14;

// 构造可逆对称正定矩阵 A = M' * M (保证 LLT 可用)
Eigen::MatrixXd make_spd_eigen(int n, unsigned seed = 42) {
    Eigen::MatrixXd M = Eigen::MatrixXd::Random(n, n);
    // Eigen::Random 范围 (-1, 1), 加 n*I 保证正定
    Eigen::MatrixXd A = M.transpose() * M + static_cast<double>(n) * Eigen::MatrixXd::Identity(n, n);
    return A;
}
}  // namespace

// =============================================================================
// §1 MatrixXD 构造与访问
// =============================================================================

TEST(LinalgDynamicTest, MatrixXD_Construct_Default) {
    MatrixXD m;
    EXPECT_EQ(m.rows(), 0);
    EXPECT_EQ(m.cols(), 0);
}

TEST(LinalgDynamicTest, MatrixXD_Construct_WithSize) {
    MatrixXD m(3, 4);
    EXPECT_EQ(m.rows(), 3);
    EXPECT_EQ(m.cols(), 4);
}

TEST(LinalgDynamicTest, MatrixXD_Access_And_Modify) {
    MatrixXD m(2, 2);
    m(0, 0) = 1.5;
    m(0, 1) = 2.5;
    m(1, 0) = 3.5;
    m(1, 1) = 4.5;
    EXPECT_DOUBLE_EQ(m(0, 0), 1.5);
    EXPECT_DOUBLE_EQ(m(0, 1), 2.5);
    EXPECT_DOUBLE_EQ(m(1, 0), 3.5);
    EXPECT_DOUBLE_EQ(m(1, 1), 4.5);
}

TEST(LinalgDynamicTest, MatrixXD_From_Eigen_Conversion) {
    Eigen::MatrixXd em = Eigen::MatrixXd::Constant(3, 3, 7.0);
    MatrixXD m(em);
    EXPECT_EQ(m.rows(), 3);
    EXPECT_EQ(m.cols(), 3);
    EXPECT_DOUBLE_EQ(m(0, 0), 7.0);
    EXPECT_DOUBLE_EQ(m(2, 2), 7.0);
}

TEST(LinalgDynamicTest, MatrixXD_To_Eigen_Conversion) {
    MatrixXD m(2, 2);
    m(0, 0) = 1.0;
    m(0, 1) = 2.0;
    m(1, 0) = 3.0;
    m(1, 1) = 4.0;
    const Eigen::MatrixXd& em = m.eigen();
    EXPECT_DOUBLE_EQ(em(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(em(1, 0), 3.0);
}

// =============================================================================
// §2 VectorXD 构造与访问
// =============================================================================

TEST(LinalgDynamicTest, VectorXD_Construct_And_Access) {
    VectorXD v(3);
    v(0) = 1.0;
    v(1) = 2.0;
    v(2) = 3.0;
    EXPECT_EQ(v.size(), 3);
    EXPECT_DOUBLE_EQ(v(0), 1.0);
    EXPECT_DOUBLE_EQ(v(2), 3.0);
}

TEST(LinalgDynamicTest, VectorXD_Dot_Product) {
    VectorXD v1(3), v2(3);
    v1(0) = 1.0; v1(1) = 2.0; v1(2) = 3.0;
    v2(0) = 4.0; v2(1) = 5.0; v2(2) = 6.0;
    // 1*4 + 2*5 + 3*6 = 32
    EXPECT_NEAR(v1.dot(v2), 32.0, TOL);
}

// =============================================================================
// §3 SVD 最小二乘求解 (OLS 核心)
// =============================================================================

TEST(LinalgDynamicTest, SVD_Solve_Overdetermined) {
    // 求解 Ax = b, A 为 4x3 (超定), 期望最小二乘解
    Eigen::MatrixXd emA(4, 3);
    emA << 1, 1, 1,
           1, 2, 4,
           1, 3, 9,
           1, 4, 16;
    MatrixXD A(emA);
    VectorXD b(4);
    b(0) = 2; b(1) = 5; b(2) = 10; b(3) = 17;
    // 真实模型 y = x^2 + 1 (拟合 x=1,2,3,4)
    // 系数 [const, x, x^2] = [1, 0, 1]
    VectorXD x = svd_solve(A, b);
    EXPECT_NEAR(x(0), 1.0, TOL);
    EXPECT_NEAR(x(1), 0.0, TOL);
    EXPECT_NEAR(x(2), 1.0, TOL);
}

TEST(LinalgDynamicTest, SVD_Solve_Square_Invertible) {
    Eigen::MatrixXd emA(3, 3);
    emA << 4, 1, 0,
           1, 3, 1,
           0, 1, 2;
    MatrixXD A(emA);
    VectorXD b(3);
    b(0) = 5; b(1) = 5; b(2) = 3;
    // A * [1;1;1] = [5;5;3]
    VectorXD x = svd_solve(A, b);
    EXPECT_NEAR(x(0), 1.0, TOL);
    EXPECT_NEAR(x(1), 1.0, TOL);
    EXPECT_NEAR(x(2), 1.0, TOL);
}

// =============================================================================
// §4 对称正定矩阵求逆 (X'X)^{-1} via LLT
// =============================================================================

TEST(LinalgDynamicTest, Inverse_Symmetric_SPD) {
    // A = [[2,1],[1,3]] SPD, A^{-1} = (1/5)[[3,-1],[-1,2]]
    Eigen::MatrixXd emA(2, 2);
    emA << 2, 1,
           1, 3;
    MatrixXD A(emA);
    MatrixXD inv = inverse_symmetric(A);
    EXPECT_NEAR(inv(0, 0), 3.0 / 5.0, TOL);
    EXPECT_NEAR(inv(0, 1), -1.0 / 5.0, TOL);
    EXPECT_NEAR(inv(1, 0), -1.0 / 5.0, TOL);
    EXPECT_NEAR(inv(1, 1), 2.0 / 5.0, TOL);
}

TEST(LinalgDynamicTest, Inverse_Symmetric_Larger) {
    // 5x5 SPD, A * A^{-1} = I
    int n = 5;
    Eigen::MatrixXd emA = make_spd_eigen(n);
    MatrixXD A(emA);
    MatrixXD inv = inverse_symmetric(A);
    // 验证 A * inv ≈ I
    Eigen::MatrixXd product = emA * inv.eigen();
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(n, n);
    EXPECT_NEAR((product - identity).cwiseAbs().maxCoeff(), 0.0, TOL);
}

// =============================================================================
// §5 Cholesky 分解
// =============================================================================

TEST(LinalgDynamicTest, Cholesky_Decompose) {
    Eigen::MatrixXd emA(3, 3);
    emA << 4, 2, 2,
           2, 5, 4,
           2, 4, 6;
    MatrixXD A(emA);
    MatrixXD L = cholesky_dynamic(A);
    // L 应为下三角, L * L' = A
    Eigen::MatrixXd emL = L.eigen();
    Eigen::MatrixXd reconstructed = emL * emL.transpose();
    EXPECT_NEAR((reconstructed - emA).cwiseAbs().maxCoeff(), 0.0, TOL);
    // 下三角验证
    EXPECT_NEAR(emL(0, 1), 0.0, TOL);
    EXPECT_NEAR(emL(0, 2), 0.0, TOL);
    EXPECT_NEAR(emL(1, 2), 0.0, TOL);
}

// =============================================================================
// §6 Ax=b 求解 (LU)
// =============================================================================

TEST(LinalgDynamicTest, Solve_Linear_System) {
    Eigen::MatrixXd emA(3, 3);
    emA << 2, 1, -1,
           -3, -1, 2,
           -2, 1, 2;
    MatrixXD A(emA);
    VectorXD b(3);
    b(0) = 8; b(1) = -11; b(2) = -3;
    // 解为 x = [2; 3; -1]
    VectorXD x = solve(A, b);
    EXPECT_NEAR(x(0), 2.0, TOL);
    EXPECT_NEAR(x(1), 3.0, TOL);
    EXPECT_NEAR(x(2), -1.0, TOL);
}

// =============================================================================
// §7 QR 分解
// =============================================================================

TEST(LinalgDynamicTest, QR_Decompose) {
    Eigen::MatrixXd emA(3, 3);
    emA << 12, -51, 4,
           6, 167, -68,
           -4, 24, -41;
    MatrixXD A(emA);
    MatrixXD Q = qr_decompose(A);
    // Q 应正交: Q * Q' = I
    Eigen::MatrixXd emQ = Q.eigen();
    Eigen::MatrixXd product = emQ * emQ.transpose();
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(3, 3);
    EXPECT_NEAR((product - identity).cwiseAbs().maxCoeff(), 0.0, TOL);
}

// =============================================================================
// §8 完整 SVD 分解
// =============================================================================

TEST(LinalgDynamicTest, SVD_Full_Decompose) {
    Eigen::MatrixXd emA(4, 2);
    emA << 1, 0,
           0, 1,
           1, 1,
           1, -1;
    MatrixXD A(emA);
    MatrixXD U;
    VectorXD S;
    MatrixXD V;
    svd_full(A, U, S, V);
    // A ≈ U * diag(S) * V'
    Eigen::MatrixXd reconstructed = U.eigen() * S.eigen().asDiagonal() * V.eigen().transpose();
    EXPECT_NEAR((reconstructed - emA).cwiseAbs().maxCoeff(), 0.0, TOL);
    // U 和 V 正交
    Eigen::MatrixXd UtU = U.eigen().transpose() * U.eigen();
    EXPECT_NEAR((UtU - Eigen::MatrixXd::Identity(U.rows(), U.rows())).cwiseAbs().maxCoeff(), 0.0, TOL);
}

// =============================================================================
// §9 矩阵运算 (加/乘/转置)
// =============================================================================

TEST(LinalgDynamicTest, Matrix_Multiply) {
    MatrixXD A(2, 3);
    A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
    A(1, 0) = 4; A(1, 1) = 5; A(1, 2) = 6;
    MatrixXD B(3, 2);
    B(0, 0) = 7; B(0, 1) = 8;
    B(1, 0) = 9; B(1, 1) = 10;
    B(2, 0) = 11; B(2, 1) = 12;
    MatrixXD C = A * B;
    // [1,2,3; 4,5,6] * [7,8; 9,10; 11,12] = [58, 64; 139, 154]
    EXPECT_DOUBLE_EQ(C(0, 0), 58.0);
    EXPECT_DOUBLE_EQ(C(0, 1), 64.0);
    EXPECT_DOUBLE_EQ(C(1, 0), 139.0);
    EXPECT_DOUBLE_EQ(C(1, 1), 154.0);
}

TEST(LinalgDynamicTest, Matrix_Transpose) {
    MatrixXD A(2, 3);
    A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
    A(1, 0) = 4; A(1, 1) = 5; A(1, 2) = 6;
    MatrixXD T = A.transpose();
    EXPECT_EQ(T.rows(), 3);
    EXPECT_EQ(T.cols(), 2);
    EXPECT_DOUBLE_EQ(T(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(T(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(T(2, 0), 3.0);
    EXPECT_DOUBLE_EQ(T(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(T(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(T(2, 1), 6.0);
}

// =============================================================================
// §10 排幻觉点验证: LLT 适用范围 (非 SPD 矩阵应抛异常)
// =============================================================================

TEST(LinalgDynamicTest, Inverse_Symmetric_NonSPD_Throws) {
    // 非正定矩阵, LLT 应失败
    Eigen::MatrixXd emA(2, 2);
    emA << 0, 1,
           1, 0;  // 特征值 ±1, 非正定
    MatrixXD A(emA);
    EXPECT_THROW(inverse_symmetric(A), std::runtime_error);
}
