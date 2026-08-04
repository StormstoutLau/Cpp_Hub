// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.3 任务 1.9 - Newey-West HAC 协方差 (排幻觉 E4/E5)
// TDD: 手算解析值 (精确分数验证), 容差 1e-12
//
// 手算数据集 1 (N=5, K=1, 无截距):
//   X=[1;2;3;4;5], y=[1;2;4;5;6], X'X=55, (X'X)^{-1}=1/55
//   β=67/55, u=[-12/55, -24/55, 19/55, 7/55, -5/55]
//   Ω_0 = 7106/3025, Ω_1 = -1264/3025, Ω_2 = -3453/3025, Ω_3 = 864/3025
//
// 手算数据集 2 (N=5, K=2, 含截距):
//   X=[[1,1],[1,2],[1,3],[1,4],[1,5]], y=[1,4,2,8,3]
//   (X'X)^{-1}=[[1.1,-0.3],[-0.3,0.1]], β=[1.2,0.8]
//   u=[-1.0, 1.2, -1.6, 3.6, -2.2]
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/econometrics/inference/hac_vcov.hpp"
#include "cpphub/econometrics/inference/hac_kernels.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;

// =============================================================================
// 辅助: 构造数据集 1 (K=1)
// =============================================================================
struct K1Dataset {
    MatrixXD X;
    VectorXD u;
    MatrixXD XtX_inv;
    K1Dataset() : X(5, 1), u(5), XtX_inv(1, 1) {
        X(0, 0) = 1.0; X(1, 0) = 2.0; X(2, 0) = 3.0; X(3, 0) = 4.0; X(4, 0) = 5.0;
        u(0) = -12.0 / 55.0; u(1) = -24.0 / 55.0; u(2) = 19.0 / 55.0;
        u(3) = 7.0 / 55.0;  u(4) = -5.0 / 55.0;
        XtX_inv(0, 0) = 1.0 / 55.0;
    }
};

// =============================================================================
// 辅助: 构造数据集 2 (K=2, 含截距)
// =============================================================================
struct K2Dataset {
    MatrixXD X;
    VectorXD u;
    MatrixXD XtX_inv;
    K2Dataset() : X(5, 2), u(5), XtX_inv(2, 2) {
        X(0, 0) = 1.0; X(0, 1) = 1.0;
        X(1, 0) = 1.0; X(1, 1) = 2.0;
        X(2, 0) = 1.0; X(2, 1) = 3.0;
        X(3, 0) = 1.0; X(3, 1) = 4.0;
        X(4, 0) = 1.0; X(4, 1) = 5.0;
        u(0) = -1.0; u(1) = 1.2; u(2) = -1.6; u(3) = 3.6; u(4) = -2.2;
        XtX_inv(0, 0) = 1.1;  XtX_inv(0, 1) = -0.3;
        XtX_inv(1, 0) = -0.3; XtX_inv(1, 1) = 0.1;
    }
};

// =============================================================================
// 测试 1: Bartlett K=1 L=1 精确值
// Ω = Ω_0 + Ω_1 = (7106-1264)/3025 = 5842/3025
// V = (1/55) * Ω * (1/55) = 5842/9150625
// =============================================================================
TEST(NeweyWestTest, Bartlett_K1_L1_ExactValue) {
    K1Dataset d;
    MatrixXD V = compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::Bartlett, 1);
    EXPECT_NEAR(V(0, 0), 5842.0 / 9150625.0, 1e-12);
}

// =============================================================================
// 测试 2: Bartlett K=1 L=2 精确值
// Ω = Ω_0 + (4/3)Ω_1 + (2/3)Ω_2 = 9356/9075
// V = 9356/27451875
// =============================================================================
TEST(NeweyWestTest, Bartlett_K1_L2_ExactValue) {
    K1Dataset d;
    MatrixXD V = compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::Bartlett, 2);
    EXPECT_NEAR(V(0, 0), 9356.0 / 27451875.0, 1e-12);
}

// =============================================================================
// 测试 3: Bartlett K=1 L=3 精确值 (排幻觉点 E5: w[L]=1/(L+1)=1/4 ≠ 0)
// Ω = Ω_0 + (3/2)Ω_1 + Ω_2 + (1/2)Ω_3 = 2189/3025
// V = 2189/9150625
// 若错误使用 w[l]=1-l/L, 则 w[3]=0, V 会等于 L=2 的值 (9356/27451875), 不等于此值
// =============================================================================
TEST(NeweyWestTest, Bartlett_K1_L3_ExactValue_WLastLagNonzero) {
    K1Dataset d;
    MatrixXD V = compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::Bartlett, 3);
    EXPECT_NEAR(V(0, 0), 2189.0 / 9150625.0, 1e-12);
    // 验证 w[3]≠0: 若 w[3]=0 (错误公式), 结果会等于 L=2 的值
    EXPECT_NE(V(0, 0), 0.0);  // w[3] 贡献使 V 为正
}

// =============================================================================
// 测试 4: max_lag=0 自动选择 NW 经验法则 (排幻觉点 E4)
// N=5: floor(4*(5/100)^(2/9)) = floor(2.0558) = 2
// 结果应与显式 max_lag=2 一致
// =============================================================================
TEST(NeweyWestTest, AutoMaxLag_NW_Rule) {
    K1Dataset d;
    MatrixXD V_auto = compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::Bartlett, 0);
    MatrixXD V_explicit = compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::Bartlett, 2);
    EXPECT_NEAR(V_auto(0, 0), V_explicit(0, 0), 1e-12);
    // 验证 NW 经验法则值
    EXPECT_NEAR(V_auto(0, 0), 9356.0 / 27451875.0, 1e-12);
}

// =============================================================================
// 测试 5: QS 内核 K=1 L=1 (通过已验证的 kernel_weight 计算期望值)
// Ω = Ω_0 + w_qs[1]*(Ω_1+Ω_1) = Ω_0 + 2*K_QS(0.5)*Ω_1
// =============================================================================
TEST(NeweyWestTest, QS_K1_L1_ViaKernelWeight) {
    K1Dataset d;
    MatrixXD V = compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::QuadraticSpectral, 1);
    // 独立计算期望值: kernel_weight 已在 test_hac_kernels.cpp 中独立验证
    Real w1 = kernel_weight(HacKernel::QuadraticSpectral, 0.5);
    Real omega_0 = 7106.0 / 3025.0;
    Real omega_1 = -1264.0 / 3025.0;
    Real expected_omega = omega_0 + w1 * 2.0 * omega_1;
    Real expected_v = expected_omega / (55.0 * 55.0);
    EXPECT_NEAR(V(0, 0), expected_v, 1e-12);
}

// =============================================================================
// 测试 6: Parzen 内核 K=1 L=1 精确值
// K_Parzen(0.5) = 1 - 6*0.25 + 6*0.125 = 0.25
// Ω = Ω_0 + 0.25*2*Ω_1 = Ω_0 + 0.5*Ω_1 = (7106-632)/3025 = 6474/3025
// V = 6474/9150625
// =============================================================================
TEST(NeweyWestTest, Parzen_K1_L1_ExactValue) {
    K1Dataset d;
    MatrixXD V = compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::Parzen, 1);
    EXPECT_NEAR(V(0, 0), 6474.0 / 9150625.0, 1e-12);
    // 验证 Parzen K(0.5)=0.25 (≠ Bartlett K(0.5)=0.5)
    EXPECT_NE(V(0, 0), 5842.0 / 9150625.0);
}

// =============================================================================
// 测试 7: Tukey-Hanning 内核 K=1 L=1 精确值
// K_Tukey(0.5) = (1+cos(π/2))/2 = 0.5 (与 Bartlett 相同)
// V = 5842/9150625 (与 Bartlett L=1 一致)
// =============================================================================
TEST(NeweyWestTest, TukeyHanning_K1_L1_EqualsBartlett) {
    K1Dataset d;
    MatrixXD V_tukey = compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::TukeyHanning, 1);
    MatrixXD V_bartlett = compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::Bartlett, 1);
    EXPECT_NEAR(V_tukey(0, 0), V_bartlett(0, 0), 1e-12);
    EXPECT_NEAR(V_tukey(0, 0), 5842.0 / 9150625.0, 1e-12);
}

// =============================================================================
// 测试 8: 所有内核 V_HAC 对称 (K=2)
// =============================================================================
TEST(NeweyWestTest, All_Kernels_Symmetric_Result) {
    K2Dataset d;
    const HacKernel kernels[] = {HacKernel::Bartlett, HacKernel::QuadraticSpectral,
                                  HacKernel::Parzen, HacKernel::TukeyHanning};
    for (HacKernel k : kernels) {
        MatrixXD V = compute_hac_vcov(d.X, d.u, d.XtX_inv, k, 1);
        EXPECT_NEAR(V(0, 1), V(1, 0), 1e-12)
            << "Symmetry failed for kernel " << to_string(k);
    }
}

// =============================================================================
// 测试 9: 所有内核 V_HAC 方差为正 (K=1)
// =============================================================================
TEST(NeweyWestTest, All_Kernels_Positive_Variance) {
    K1Dataset d;
    const HacKernel kernels[] = {HacKernel::Bartlett, HacKernel::QuadraticSpectral,
                                  HacKernel::Parzen, HacKernel::TukeyHanning};
    for (HacKernel k : kernels) {
        MatrixXD V = compute_hac_vcov(d.X, d.u, d.XtX_inv, k, 1);
        EXPECT_GT(V(0, 0), 0.0) << "Non-positive variance for kernel " << to_string(k);
    }
}

// =============================================================================
// 测试 10: Bartlett K=2 含截距 L=1 精确值
// Ω = [[6, 25.2], [25.2, 116.72]]
// V = [[708/625, -276/625], [-276/625, 122/625]] = [[1.1328, -0.4416], [-0.4416, 0.1952]]
// =============================================================================
TEST(NeweyWestTest, Bartlett_K2_WithIntercept_L1_ExactValue) {
    K2Dataset d;
    MatrixXD V = compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::Bartlett, 1);
    EXPECT_NEAR(V(0, 0), 708.0 / 625.0, 1e-12);   // 1.1328
    EXPECT_NEAR(V(0, 1), -276.0 / 625.0, 1e-12);  // -0.4416
    EXPECT_NEAR(V(1, 0), -276.0 / 625.0, 1e-12);  // -0.4416
    EXPECT_NEAR(V(1, 1), 122.0 / 625.0, 1e-12);   // 0.1952
}

// =============================================================================
// 测试 11: Bartlett K=2 含截距 L=2 精确值
// Ω_2 = [[9.44, 39.68], [20.8, 92.16]]
// Ω = [[502/75, 614/25], [614/25, 2442/25]]
// V = [[1276/1875, -142/625], [-142/625, 66/625]]
//   = [[0.680533..., -0.2272], [-0.2272, 0.1056]]
// =============================================================================
TEST(NeweyWestTest, Bartlett_K2_WithIntercept_L2_ExactValue) {
    K2Dataset d;
    MatrixXD V = compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::Bartlett, 2);
    EXPECT_NEAR(V(0, 0), 1276.0 / 1875.0, 1e-12);  // 0.680533...
    EXPECT_NEAR(V(0, 1), -142.0 / 625.0, 1e-12);   // -0.2272
    EXPECT_NEAR(V(1, 0), -142.0 / 625.0, 1e-12);   // -0.2272
    EXPECT_NEAR(V(1, 1), 66.0 / 625.0, 1e-12);     // 0.1056
}

// =============================================================================
// 测试 12: max_lag >= T 抛 std::invalid_argument
// =============================================================================
TEST(NeweyWestTest, MaxLag_Exceeds_Observations_Throws) {
    K1Dataset d;
    // max_lag = 5 = T (边界, 应抛异常)
    EXPECT_THROW(compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::Bartlett, 5),
                 std::invalid_argument);
    // max_lag = 6 > T (超出, 应抛异常)
    EXPECT_THROW(compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::Bartlett, 6),
                 std::invalid_argument);
}

// =============================================================================
// 测试 13: X 行数与 residuals 长度不匹配抛 std::runtime_error
// =============================================================================
TEST(NeweyWestTest, Dimension_Mismatch_Throws) {
    K1Dataset d;
    VectorXD short_u(4);
    short_u(0) = -12.0 / 55.0; short_u(1) = -24.0 / 55.0;
    short_u(2) = 19.0 / 55.0;  short_u(3) = 7.0 / 55.0;
    EXPECT_THROW(compute_hac_vcov(d.X, short_u, d.XtX_inv, HacKernel::Bartlett, 1),
                 std::runtime_error);
}

// =============================================================================
// 测试 14: prewhiten=true 抛 std::runtime_error (暂未实现)
// =============================================================================
TEST(NeweyWestTest, Prewhiten_True_Throws) {
    K1Dataset d;
    EXPECT_THROW(compute_hac_vcov(d.X, d.u, d.XtX_inv, HacKernel::Bartlett, 1, true),
                 std::runtime_error);
}

// =============================================================================
// 测试 15: 输出维度正确 (k × k)
// =============================================================================
TEST(NeweyWestTest, Result_Dimension_Correct) {
    K1Dataset d1;
    MatrixXD V1 = compute_hac_vcov(d1.X, d1.u, d1.XtX_inv, HacKernel::Bartlett, 1);
    EXPECT_EQ(V1.rows(), 1u);
    EXPECT_EQ(V1.cols(), 1u);

    K2Dataset d2;
    MatrixXD V2 = compute_hac_vcov(d2.X, d2.u, d2.XtX_inv, HacKernel::Bartlett, 1);
    EXPECT_EQ(V2.rows(), 2u);
    EXPECT_EQ(V2.cols(), 2u);
}
