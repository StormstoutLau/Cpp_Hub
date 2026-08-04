// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.3 任务 1.8 - HAC 内核 (排幻觉 E4/E5, Parzen 审计修复含 |u|)
#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <vector>

#include "cpphub/econometrics/inference/hac_kernels.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;

// =============================================================================
// Bartlett (Newey-West 1987)
// =============================================================================

TEST(HacKernelsTest, Bartlett_Kernel_Boundary_Values) {
    EXPECT_NEAR(kernel_weight(HacKernel::Bartlett, 0.0), 1.0, 1e-12);
    EXPECT_NEAR(kernel_weight(HacKernel::Bartlett, 0.5), 0.5, 1e-12);
    EXPECT_NEAR(kernel_weight(HacKernel::Bartlett, 1.0), 0.0, 1e-12);
    EXPECT_NEAR(kernel_weight(HacKernel::Bartlett, 1.5), 0.0, 1e-12);
}

TEST(HacKernelsTest, Bartlett_Kernel_Symmetry) {
    EXPECT_NEAR(kernel_weight(HacKernel::Bartlett, -0.3),
                kernel_weight(HacKernel::Bartlett, 0.3), 1e-12);
    EXPECT_NEAR(kernel_weight(HacKernel::Bartlett, -0.7),
                kernel_weight(HacKernel::Bartlett, 0.7), 1e-12);
}

TEST(HacKernelsTest, Bartlett_Weights_Vector) {
    // 排幻觉点 E5: 分母是 L+1=5, 非 L=4
    const auto w = kernel_weights(HacKernel::Bartlett, 4);
    ASSERT_EQ(w.size(), 5u);
    EXPECT_NEAR(w[0], 1.0, 1e-12);
    EXPECT_NEAR(w[1], 4.0 / 5.0, 1e-12);
    EXPECT_NEAR(w[2], 3.0 / 5.0, 1e-12);
    EXPECT_NEAR(w[3], 2.0 / 5.0, 1e-12);
    EXPECT_NEAR(w[4], 1.0 / 5.0, 1e-12);  // w[L] = 1/(L+1) = 0.2, 非 0
}

// =============================================================================
// Quadratic Spectral (Andrews 1991)
// =============================================================================

TEST(HacKernelsTest, QuadraticSpectral_Kernel_Boundary) {
    EXPECT_NEAR(kernel_weight(HacKernel::QuadraticSpectral, 0.0), 1.0, 1e-9);  // 极限
    EXPECT_GT(kernel_weight(HacKernel::QuadraticSpectral, 1.0), 0.0);
    EXPECT_NE(kernel_weight(HacKernel::QuadraticSpectral, 2.0), 0.0);  // 不截断 (有振荡小负瓣)
}

TEST(HacKernelsTest, QuadraticSpectral_Kernel_Value) {
    // K(1) = 25/(12*pi^2) * [sin(6pi/5)/(6pi/5) - cos(6pi/5)] ≈ 0.1379
    EXPECT_NEAR(kernel_weight(HacKernel::QuadraticSpectral, 1.0), 0.13786, 1e-3);
}

// =============================================================================
// Parzen (Gallant 1987), 审计修复含 |u|
// =============================================================================

TEST(HacKernelsTest, Parzen_Kernel_Boundary_Values) {
    EXPECT_NEAR(kernel_weight(HacKernel::Parzen, 0.0), 1.0, 1e-12);
    EXPECT_NEAR(kernel_weight(HacKernel::Parzen, 0.5), 0.25, 1e-12);  // 1-6*0.25+6*0.125
    EXPECT_NEAR(kernel_weight(HacKernel::Parzen, 1.0), 0.0, 1e-12);
    EXPECT_NEAR(kernel_weight(HacKernel::Parzen, 1.5), 0.0, 1e-12);  // 超出 support
}

TEST(HacKernelsTest, Parzen_Kernel_Symmetry) {
    EXPECT_NEAR(kernel_weight(HacKernel::Parzen, -0.3),
                kernel_weight(HacKernel::Parzen, 0.3), 1e-12);
    EXPECT_NEAR(kernel_weight(HacKernel::Parzen, -0.7),
                kernel_weight(HacKernel::Parzen, 0.7), 1e-12);  // 0.7>0.5 用 2(1-|u|)^3
}

TEST(HacKernelsTest, Parzen_Kernel_Piecewise) {
    EXPECT_NEAR(kernel_weight(HacKernel::Parzen, 0.3), 0.622, 1e-12);  // 1-6*0.09+6*0.027
    EXPECT_NEAR(kernel_weight(HacKernel::Parzen, 0.7), 0.054, 1e-12);  // 2*(0.3)^3
}

// =============================================================================
// Tukey-Hanning
// =============================================================================

TEST(HacKernelsTest, TukeyHanning_Kernel_Boundary) {
    EXPECT_NEAR(kernel_weight(HacKernel::TukeyHanning, 0.0), 1.0, 1e-12);
    EXPECT_NEAR(kernel_weight(HacKernel::TukeyHanning, 1.0), 0.0, 1e-12);
    EXPECT_NEAR(kernel_weight(HacKernel::TukeyHanning, 0.5), 0.5, 1e-12);  // (1+cos(pi/2))/2
}

TEST(HacKernelsTest, TukeyHanning_Symmetry) {
    EXPECT_NEAR(kernel_weight(HacKernel::TukeyHanning, -0.4),
                kernel_weight(HacKernel::TukeyHanning, 0.4), 1e-12);
}

// =============================================================================
// 批量权重向量 (所有内核)
// =============================================================================

TEST(HacKernelsTest, Kernel_Weights_All_Kernels) {
    const HacKernel kernels[] = {HacKernel::Bartlett, HacKernel::QuadraticSpectral,
                                 HacKernel::Parzen, HacKernel::TukeyHanning};
    for (HacKernel k : kernels) {
        const auto w = kernel_weights(k, 5);
        ASSERT_EQ(w.size(), 6u);  // max_lag + 1
        EXPECT_NEAR(w[0], 1.0, 1e-12);
        for (Real val : w) EXPECT_GE(val, 0.0);
    }
}

// =============================================================================
// select_max_lag (NW 1987 经验法则 / Andrews 1991 自动带宽)
// 排幻觉点 E4: andrews_optimal=false (默认) → NW, true → Andrews (需 ar1_coef)
// =============================================================================

TEST(HacKernelsTest, Select_Max_Lag_Bartlett_NW_Default) {
    // 排幻觉点 E4: andrews_optimal=false (默认) → NW 经验法则 floor(4*(T/100)^(2/9))
    // T=100: floor(4*1) = 4
    EXPECT_EQ(select_max_lag(100, HacKernel::Bartlett), 4u);
    EXPECT_EQ(select_max_lag(100, HacKernel::Bartlett, false, 0.0), 4u);
    // T=400: floor(4*(4)^(2/9)) = floor(4*1.3536) = floor(5.414) = 5
    EXPECT_EQ(select_max_lag(400, HacKernel::Bartlett, false, 0.0), 5u);
}

TEST(HacKernelsTest, Select_Max_Lag_Bartlett_Andrews_Optimal) {
    // 排幻觉点 E4: andrews_optimal=true, ar1_coef=0.5 → Andrews 1991 公式
    // alpha(1) = 4*0.25/(0.75)^2 = 1.7778
    // L* = floor(1.1447 * (1.7778 * 100)^(1/3)) = floor(1.1447 * 5.6202) = floor(6.4338) = 6
    const Size lag = select_max_lag(100, HacKernel::Bartlett, true, 0.5);
    EXPECT_GE(lag, 5u);
    EXPECT_LE(lag, 7u);  // 容差 ±1
}

TEST(HacKernelsTest, Select_Max_Lag_Bartlett_Andrews_Zero_Rho_Fallback) {
    // 排幻觉点 E4: andrews_optimal=true 但 ar1_coef=0 (无自相关) → 退化到 NW
    EXPECT_EQ(select_max_lag(100, HacKernel::Bartlett, true, 0.0), 4u);
}

TEST(HacKernelsTest, Select_Max_Lag_QS_Andrews) {
    // QS 用 b* = 1.3221 * (alpha(2)*T)^(1/5)
    // alpha(2) = 4*0.25/(0.75)^4 = 1.0/0.3164 = 3.1605
    // b* = floor(1.3221 * (3.1605 * 100)^(1/5)) = floor(1.3221 * 3.1698) = floor(4.191) = 4
    const Size bw = select_max_lag(100, HacKernel::QuadraticSpectral, true, 0.5);
    EXPECT_GT(bw, 0u);
    EXPECT_LE(bw, 6u);  // 容差 ±2
}

// =============================================================================
// to_string
// =============================================================================

TEST(HacKernelsTest, ToString_All_Kernels) {
    EXPECT_EQ(to_string(HacKernel::Bartlett), "Bartlett");
    EXPECT_EQ(to_string(HacKernel::QuadraticSpectral), "QuadraticSpectral");
    EXPECT_EQ(to_string(HacKernel::Parzen), "Parzen");
    EXPECT_EQ(to_string(HacKernel::TukeyHanning), "TukeyHanning");
}