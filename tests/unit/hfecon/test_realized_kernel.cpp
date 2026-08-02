// =============================================================================
// test_realized_kernel.cpp
// Phase 5 v1.4.1 - HFE Realized Kernel (BNS 2008 ECTA)
//
// 基准来源: R highfrequency 1.0.3 rKernelCov (generate_v141_baselines.R)
// 容差: 1e-15 (核函数解析值) / 1e-12 (R baseline 对标)
//
// SOURCE:
//   [BNS 2008] Barndorff-Nielsen, Hansen, Lunde, Shephard,
//              Econometrica 76(6), 1481-1536, doi:10.1111/j.1468-0262.2008.00837.x
//   [H-L 2006] Hansen & Lunde, JBES 24(2), 127-161, doi:10.1198/073500106000000072
//
// R 基准生成: tests/fixtures/hfe/generate_v141_baselines.R
// R 版本: 4.6.1 + highfrequency 1.0.3
// 生成时间: 2026-08-02
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/measures/kernels.hpp"
#include "cpphub/hfecon/measures/realized_kernel.hpp"
#include "cpphub/hfecon/noise/noise_variance.hpp"
#include "cpphub/hfecon/noise/bandwidth.hpp"
#include <vector>
#include <cmath>
#include <random>

using namespace cpphub::v1::hfecon;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

// ============================================================================
// 容差
// ============================================================================
constexpr Real TOL_KERNEL  = 1e-15;  // 核函数解析值
constexpr Real TOL_STRICT  = 1e-12;  // R baseline 对标 (无噪声合成数据)
constexpr Real TOL_STANDARD = 1e-10; // R baseline 对标 (默认)

// ============================================================================
// R baseline 数据 (generate_v141_baselines.R 2026-08-02 生成)
// R highfrequency 1.0.3, R 4.6.1
// ============================================================================

// B.1: 已知收益率序列 r = [0.01, -0.02, 0.03, -0.01, 0.02], n=5
//   gamma_0 (RV) = 0.0019
//   gamma_1 = -0.0013 (手算 + R 验证)
const std::vector<Real> RET_B1 = {0.01, -0.02, 0.03, -0.01, 0.02};
constexpr Real B1_RV      = 0.0019;
constexpr Real B1_GAMMA_1 = -0.0013;

// B.1 R baseline (rKernelCov 输出值, 2026-08-02)
// 所有核 H=1 DOF=T: -0.00135 (因为 k(0)=1 for all)
constexpr Real B1_H1_DOF_T  = -0.00135;
constexpr Real B1_H1_DOF_F  = -0.0007;

// 各核 H=2 DOF=T
constexpr Real B1_RECT_H2_DOF_T   = 0.002316666666666667;
constexpr Real B1_BART_H2_DOF_T   = 0.000483333333333333;
constexpr Real B1_SECOND_H2_DOF_T = 0.0014;
constexpr Real B1_PARZEN_H2_DOF_T = -0.000433333333333333;
constexpr Real B1_SEVENTH_H2_DOF_T = 0.0014859375;
constexpr Real B1_EIGHTH_H2_DOF_T  = 0.00178671875;

// 各核 H=3 DOF=T
constexpr Real B1_BART_H3_DOF_T   = 0.000261111111111111;
constexpr Real B1_PARZEN_H3_DOF_T = 0.000501851851851852;
constexpr Real B1_EPAN_H3_DOF_T   = 0.000520370370370370;
constexpr Real B1_MTH_H3_DOF_T    = 0.000089594116926347;

// B.2: GBM CASE3 returns (n=100, seed=42, sigma=0.01)
// R baseline (rKernelCov, DOF=TRUE)
constexpr Real B2_RV = 0.010558420207822800;
constexpr Real B2_BART_H1  = 0.011935982508821092;
constexpr Real B2_BART_H5  = 0.012321648318276536;
constexpr Real B2_PARZEN_H5  = 0.012058235293304572;
constexpr Real B2_PARZEN_H10 = 0.011293887673364924;
constexpr Real B2_TH_H10     = 0.009991940999737053;

// CASE3 价格序列 (R set.seed(42) GBM, 与 v1.4.0 共享)
std::vector<Real> r_case3_prices() {
    return std::vector<Real>{
        101.38039917595452, 100.80951930678135, 101.17625276642408, 101.81858984699473,
        102.23104330040793, 102.12260864846189, 103.67793931303507, 103.57984520742802,
        105.69176743555116, 105.62550457608084, 107.01281132027738, 109.48800645723944,
        107.97788160279309, 107.67727062775589, 107.5338095049555,  108.21985032207174,
        107.91266902871402, 105.08375782085513, 102.5502637483704,  103.91301863018465,
        103.59486924448341, 101.76586358015204, 101.59106069971737, 102.83258659625074,
        104.80004776280842, 104.34988550719638, 104.08177023796344, 102.26272240106542,
        102.73431454136686, 102.07891967627012, 102.54489858699878, 103.27022650783616,
        104.34473178668154, 103.71127977897947, 104.23629964143039, 102.46183081603857,
        101.66120415375512, 100.79983118852508, 98.395454060551174, 98.431003483923988,
        98.633978964810709, 98.278495919828828, 99.02643907679554,  98.309417634148389,
        96.973429388244014, 97.394057491444727, 96.60700611354865,  98.012231123074883,
        97.5902719916715,   98.232222712900949, 98.548966622169928, 97.779520992766194,
        99.332462747317521, 99.973127669179803, 100.06290448109721, 100.3400121863053,
        101.02393092902517, 101.11472441741427, 98.13311336912875,  98.413076475124512,
        98.052332359422962, 98.234123562975157, 98.807338936726623, 100.2001064478679,
        99.474002678242329, 100.7781691944837,  101.11719977608917, 102.17277970661081,
        103.11785780390004, 103.86389771145825, 102.78610482778274, 102.69344754230983,
        103.33576022608229, 102.35511239852755, 101.8010046441415,  102.39418642982434,
        103.18378568270319, 103.66342999507516, 102.74925863684794, 101.62543305377729,
        103.17441434013229, 103.44086674398767, 103.53239054954682, 103.40729910498969,
        102.17962169825519, 102.80687524182999, 102.58388274191515, 102.39657502742747,
        103.3567636733232,  104.20962124264086, 105.67048531652755, 105.16850611660557,
        105.85469687979499, 107.33754273540536, 106.15184672834231, 105.2420209750828,
        104.05767083080239, 102.55027156853281, 102.63232670461113, 103.30491982877388
    };
}

// 核函数解析值 (来自 R highfrequency 1.0.3 KK() 源码公式)
// k(0) = 1 for all kernels
// k(0.5):
constexpr Real K_RECT_05   = 1.0;
constexpr Real K_BART_05   = 0.5;
constexpr Real K_SECOND_05 = 0.75;       // 1 - 2*(0.5^3) = 1 - 0.25
constexpr Real K_EPAN_05   = 0.75;       // 1 - 0.5^2
constexpr Real K_CUBIC_05  = 0.5;        // 1 - 3*0.25 + 2*0.125 = 0.5
constexpr Real K_FIFTH_05  = 0.5;        // 1 - 10*0.125 + 15*0.0625 - 6*0.03125
constexpr Real K_SIXTH_05  = 0.65625;    // 1 - 15*0.0625 + 24*0.03125 - 10*0.015625
constexpr Real K_SEVENTH_05 = 0.7734375; // 1 - 21*0.03125 + 35*0.015625 - 15*0.0078125
constexpr Real K_EIGHTH_05  = 0.85546875;// 1 - 28*0.015625 + 48*0.0078125 - 21*0.00390625
constexpr Real K_PARZEN_05 = 0.25;       // 1 - 6*0.25 + 6*0.125 = 0.25 (x ≤ 0.5 branch)
constexpr Real K_TH_05     = 0.5;        // (1 + sin(π/2 - π/2))/2 = (1+0)/2
// ModifiedTukeyHanning k(0.5) = (1 - sin(π/2 - π*0.25))/2 = (1 - sin(π/4))/2 = (1 - √2/2)/2
constexpr Real K_MTH_05    = 0.14644660940672624;  // (1 - √2/2)/2

// k(1) (边界值, 部分 R 公式在此处非零)
constexpr Real K_SECOND_1  = -1.0;  // Discovery: R 1-2*1^3 = -1 (核函数为负!)

}  // namespace

// ============================================================================
// A. 核函数单测 (11 个, 对照 R highfrequency 1.0.3 KK() 源码公式)
// ============================================================================

TEST(HFE_Kernels, Rectangular) {
    // R KK() case 0: 恒为 1, 无支撑限制
    EXPECT_NEAR(kernel_value(KernelType::Rectangular, 0.0),  1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Rectangular, 0.5),  K_RECT_05, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Rectangular, 1.0),  1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Rectangular, 1.5),  1.0, TOL_KERNEL); // 不归零
}

TEST(HFE_Kernels, Bartlett) {
    // R KK() case 1: k(x) = 1 - x
    EXPECT_NEAR(kernel_value(KernelType::Bartlett, 0.0), 1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Bartlett, 0.5), K_BART_05, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Bartlett, 1.0), 0.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Bartlett, 1.5), 0.0, TOL_KERNEL); // 支撑外
}

TEST(HFE_Kernels, Second) {
    // R KK() case 2: k(x) = 1 - 2x³ (Discovery: R 偏离 BNS 2008 论文 1-x²)
    EXPECT_NEAR(kernel_value(KernelType::Second, 0.0), 1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Second, 0.5), K_SECOND_05, TOL_KERNEL);
    // k(1) = 1 - 2 = -1 (Discovery: 核函数为负, R 实现如此)
    EXPECT_NEAR(kernel_value(KernelType::Second, 1.0), K_SECOND_1, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Second, 1.5), 0.0, TOL_KERNEL);
}

TEST(HFE_Kernels, Epanechnikov) {
    // R KK() case 3: k(x) = 1 - x²
    EXPECT_NEAR(kernel_value(KernelType::Epanechnikov, 0.0), 1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Epanechnikov, 0.5), K_EPAN_05, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Epanechnikov, 1.0), 0.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Epanechnikov, 1.5), 0.0, TOL_KERNEL);
}

TEST(HFE_Kernels, Cubic) {
    // R KK() case 4: k(x) = 1 - 3x² + 2x³
    EXPECT_NEAR(kernel_value(KernelType::Cubic, 0.0), 1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Cubic, 0.5), K_CUBIC_05, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Cubic, 1.0), 0.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Cubic, 1.5), 0.0, TOL_KERNEL);
}

TEST(HFE_Kernels, Fifth) {
    // R KK() case 5: k(x) = 1 - 10x³ + 15x⁴ - 6x⁵
    EXPECT_NEAR(kernel_value(KernelType::Fifth, 0.0), 1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Fifth, 0.5), K_FIFTH_05, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Fifth, 1.0), 0.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Fifth, 1.5), 0.0, TOL_KERNEL);
}

TEST(HFE_Kernels, Sixth) {
    // R KK() case 6: k(x) = 1 - 15x⁴ + 24x⁵ - 10x⁶
    EXPECT_NEAR(kernel_value(KernelType::Sixth, 0.0), 1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Sixth, 0.5), K_SIXTH_05, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Sixth, 1.0), 0.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Sixth, 1.5), 0.0, TOL_KERNEL);
}

TEST(HFE_Kernels, Seventh) {
    // R KK() case 7: k(x) = 1 - 21x⁵ + 35x⁶ - 15x⁷
    EXPECT_NEAR(kernel_value(KernelType::Seventh, 0.0), 1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Seventh, 0.5), K_SEVENTH_05, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Seventh, 1.0), 0.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Seventh, 1.5), 0.0, TOL_KERNEL);
}

TEST(HFE_Kernels, Eighth) {
    // R KK() case 8: k(x) = 1 - 28x⁶ + 48x⁷ - 21x⁸
    EXPECT_NEAR(kernel_value(KernelType::Eighth, 0.0), 1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Eighth, 0.5), K_EIGHTH_05, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Eighth, 1.0), 0.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Eighth, 1.5), 0.0, TOL_KERNEL);
}

TEST(HFE_Kernels, Parzen) {
    // R KK() case 9: 分段
    //   x ≤ 0.5: 1 - 6x² + 6x³
    //   x > 0.5: 2(1-x)³
    EXPECT_NEAR(kernel_value(KernelType::Parzen, 0.0), 1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Parzen, 0.25), 0.71875, TOL_KERNEL); // 1-6*0.0625+6*0.015625
    EXPECT_NEAR(kernel_value(KernelType::Parzen, 0.5), K_PARZEN_05, TOL_KERNEL);
    // x=0.5 两分支交汇: 1-6*0.25+6*0.125 = 0.25; 2*(1-0.5)^3 = 0.25
    EXPECT_NEAR(kernel_value(KernelType::Parzen, 0.75), 0.03125, TOL_KERNEL); // 2*(0.25)^3
    EXPECT_NEAR(kernel_value(KernelType::Parzen, 1.0), 0.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::Parzen, 1.5), 0.0, TOL_KERNEL);
}

TEST(HFE_Kernels, TukeyHanningAndModified) {
    // R KK() case 10: TukeyHanning (1 + sin(π/2 - πx))/2 = (1 + cos(πx))/2
    EXPECT_NEAR(kernel_value(KernelType::TukeyHanning, 0.0), 1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::TukeyHanning, 0.5), K_TH_05, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::TukeyHanning, 1.0), 0.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::TukeyHanning, 1.5), 0.0, TOL_KERNEL);

    // R KK() case 11: ModifiedTukeyHanning (1 - sin(π/2 - π(1-x)²))/2
    EXPECT_NEAR(kernel_value(KernelType::ModifiedTukeyHanning, 0.0), 1.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::ModifiedTukeyHanning, 0.5), K_MTH_05, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::ModifiedTukeyHanning, 1.0), 0.0, TOL_KERNEL);
    EXPECT_NEAR(kernel_value(KernelType::ModifiedTukeyHanning, 1.5), 0.0, TOL_KERNEL);

    // parse_kernel_type 验证
    EXPECT_EQ(parse_kernel_type("rectangular"), KernelType::Rectangular);
    EXPECT_EQ(parse_kernel_type("Bartlett"), KernelType::Bartlett);
    EXPECT_EQ(parse_kernel_type("Parzen"), KernelType::Parzen);
    EXPECT_EQ(parse_kernel_type("TukeyHanning"), KernelType::TukeyHanning);
    EXPECT_EQ(parse_kernel_type("ModifiedTukeyHanning"), KernelType::ModifiedTukeyHanning);
    EXPECT_THROW(parse_kernel_type("unknown"), std::invalid_argument);
}

// ============================================================================
// B. R baseline 对标 (3 个, 硬编码 R rKernelCov 输出值)
// ============================================================================

TEST(HFE_RealizedKernel, RBaselineB1) {
    // R baseline: ret_b1 = [0.01, -0.02, 0.03, -0.01, 0.02], n=5
    // 严格对标 R rKernelCov, 容差 1e-12

    // RV (gamma_0) 和 gamma_1 验证
    auto res_h1 = RealizedKernel::estimate(RET_B1, KernelType::Rectangular, 1, false);
    EXPECT_NEAR(res_h1.rv, B1_RV, TOL_STRICT);
    EXPECT_NEAR(res_h1.gamma_1, B1_GAMMA_1, TOL_STRICT);
    EXPECT_EQ(res_h1.n_obs, 5u);
    EXPECT_EQ(res_h1.bandwidth, 1u);

    // H=1 所有核结果相同 (k(0)=1 for all), DOF=T
    for (KernelType kt : {KernelType::Rectangular, KernelType::Bartlett,
                          KernelType::Second, KernelType::Parzen, KernelType::TukeyHanning}) {
        auto res = RealizedKernel::estimate(RET_B1, kt, 1, true);
        EXPECT_NEAR(res.rk, B1_H1_DOF_T, TOL_STRICT)
            << "H=1 DOF=T should be same for all kernels";
    }
    // H=1 DOF=F
    auto res_h1_nodof = RealizedKernel::estimate(RET_B1, KernelType::Rectangular, 1, false);
    EXPECT_NEAR(res_h1_nodof.rk, B1_H1_DOF_F, TOL_STRICT);

    // H=2 DOF=T 各核
    EXPECT_NEAR(RealizedKernel::estimate(RET_B1, KernelType::Rectangular, 2, true).rk,
                B1_RECT_H2_DOF_T, TOL_STRICT);
    EXPECT_NEAR(RealizedKernel::estimate(RET_B1, KernelType::Bartlett, 2, true).rk,
                B1_BART_H2_DOF_T, TOL_STRICT);
    EXPECT_NEAR(RealizedKernel::estimate(RET_B1, KernelType::Second, 2, true).rk,
                B1_SECOND_H2_DOF_T, TOL_STRICT);
    EXPECT_NEAR(RealizedKernel::estimate(RET_B1, KernelType::Parzen, 2, true).rk,
                B1_PARZEN_H2_DOF_T, TOL_STRICT);
    EXPECT_NEAR(RealizedKernel::estimate(RET_B1, KernelType::Seventh, 2, true).rk,
                B1_SEVENTH_H2_DOF_T, TOL_STRICT);
    EXPECT_NEAR(RealizedKernel::estimate(RET_B1, KernelType::Eighth, 2, true).rk,
                B1_EIGHTH_H2_DOF_T, TOL_STRICT);

    // H=3 DOF=T 各核
    EXPECT_NEAR(RealizedKernel::estimate(RET_B1, KernelType::Bartlett, 3, true).rk,
                B1_BART_H3_DOF_T, TOL_STRICT);
    EXPECT_NEAR(RealizedKernel::estimate(RET_B1, KernelType::Parzen, 3, true).rk,
                B1_PARZEN_H3_DOF_T, TOL_STRICT);
    EXPECT_NEAR(RealizedKernel::estimate(RET_B1, KernelType::Epanechnikov, 3, true).rk,
                B1_EPAN_H3_DOF_T, TOL_STRICT);
    EXPECT_NEAR(RealizedKernel::estimate(RET_B1, KernelType::ModifiedTukeyHanning, 3, true).rk,
                B1_MTH_H3_DOF_T, TOL_STRICT);
}

TEST(HFE_RealizedKernel, RBaselineB2) {
    // R baseline: GBM CASE3 returns (n=100, seed=42, sigma=0.01)
    // 严格对标 R rKernelCov, DOF=TRUE, 容差 1e-12
    auto prices = r_case3_prices();
    auto returns = std::vector<Real>(prices.size(), 0.0);
    for (Size i = 1; i < prices.size(); ++i) {
        returns[i] = std::log(prices[i] / prices[i - 1]);
    }

    // RV
    auto res_h1 = RealizedKernel::estimate(returns, KernelType::Bartlett, 1, true);
    EXPECT_NEAR(res_h1.rv, B2_RV, TOL_STRICT);

    // Bartlett H=1/5
    EXPECT_NEAR(RealizedKernel::estimate(returns, KernelType::Bartlett, 1, true).rk,
                B2_BART_H1, TOL_STRICT);
    EXPECT_NEAR(RealizedKernel::estimate(returns, KernelType::Bartlett, 5, true).rk,
                B2_BART_H5, TOL_STRICT);

    // Parzen H=5/10
    EXPECT_NEAR(RealizedKernel::estimate(returns, KernelType::Parzen, 5, true).rk,
                B2_PARZEN_H5, TOL_STRICT);
    EXPECT_NEAR(RealizedKernel::estimate(returns, KernelType::Parzen, 10, true).rk,
                B2_PARZEN_H10, TOL_STRICT);

    // TukeyHanning H=10
    EXPECT_NEAR(RealizedKernel::estimate(returns, KernelType::TukeyHanning, 10, true).rk,
                B2_TH_H10, TOL_STRICT);

    // estimate_from_prices 一致性
    auto res_from_prices = RealizedKernel::estimate_from_prices(
        prices, KernelType::Bartlett, 5, true);
    EXPECT_NEAR(res_from_prices.rk, B2_BART_H5, TOL_STRICT);
}

TEST(HFE_RealizedKernel, NoiseRejectionAndGamma1) {
    // Part 1: gamma_1 验证 (使用 R baseline ret_b1)
    // RK(H=1, Rectangular, DOF=F) = gamma_0 + 2*k(0)*gamma_1 = RV + 2*gamma_1
    // => gamma_1 = (RK - RV) / 2
    auto res = RealizedKernel::estimate(RET_B1, KernelType::Rectangular, 1, false);
    Real gamma_1_extracted = (res.rk - res.rv) / 2.0;
    EXPECT_NEAR(gamma_1_extracted, B1_GAMMA_1, TOL_STRICT);
    EXPECT_NEAR(res.gamma_1, B1_GAMMA_1, TOL_STRICT);  // 直接字段验证

    // Part 2: 噪声稳健性 — MA(1) 收益率结构 (BNS 2008 标准噪声模型)
    // 关键: 价格观测噪声 p_obs = p_eff + ε 产生 MA(1) 收益率:
    //   r_obs_i = r_eff_i + ε_i - ε_{i-1}
    //   γ₁ = E[r_i * r_{i+1}] = -σ²_ε < 0  (MA(1) 负自协方差)
    //   γ_h = 0 for h ≥ 2
    //   RK = γ₀ + 2*k(0)*γ₁ < γ₀ = RV  (利用负自协方差修正噪声偏差)
    // 注意: 纯 i.i.d. 噪声 (直接加到收益率) 的 γ_h ≈ 0 (h>0), RK ≈ RV,
    //       无法体现 BNS 2008 噪声修正; 必须构造 MA(1) 结构
    std::mt19937_64 gen(12345);
    const Size n_obs = 500;
    const Real sigma_signal = 0.01;   // 有效收益率波动
    const Real sigma_noise  = 0.005;  // 微结构噪声 σ_ε
    std::normal_distribution<Real> sig_dist(0.0, sigma_signal);
    std::normal_distribution<Real> noi_dist(0.0, sigma_noise);

    // 直接构造 MA(1) 收益率: r_obs[i] = signal[i] + ε[i] - ε[i-1]
    // 理论: E[γ₀] = σ²_s + 2σ²_ε,  E[γ₁] = -σ²_ε,  E[γ_h] = 0 (h≥2)
    // RV ≈ n*(σ²_s + 2σ²_ε) = n*0.00015 = 0.075
    // RK(Bartlett,H≥1) ≈ n*σ²_s = n*0.0001 = 0.05  (噪声被修正)
    std::vector<Real> ma1_rets(n_obs);
    Real prev_eps = noi_dist(gen);  // ε_{-1}
    for (Size i = 0; i < n_obs; ++i) {
        Real sig = sig_dist(gen);
        Real eps = noi_dist(gen);
        ma1_rets[i] = sig + eps - prev_eps;
        prev_eps = eps;
    }

    Real rv_ma1 = 0.0;
    for (auto r : ma1_rets) rv_ma1 += r * r;

    // RK with Bartlett H=5 (利用 MA(1) 负自协方差修正噪声)
    auto rk_result = RealizedKernel::estimate(ma1_rets, KernelType::Bartlett, 5, true);
    // BNS 2008: MA(1) 噪声使 γ₁ < 0, RK = γ₀ + 2*k(0)*γ₁ + ... < γ₀ = RV
    EXPECT_LT(rk_result.rk, rv_ma1)
        << "RK should reduce noise bias with MA(1) noise: RK=" << rk_result.rk
        << " < RV=" << rv_ma1;
    // γ₁ 应为负 (MA(1) 噪声特征, 理论值 -σ²_ε = -2.5e-5)
    EXPECT_LT(rk_result.gamma_1, 0.0)
        << "gamma_1 should be negative for MA(1) noise structure";

    // Part 3: 异常处理
    EXPECT_THROW(RealizedKernel::estimate(RET_B1, KernelType::Bartlett, 0, true),
                 std::invalid_argument);  // H=0
    EXPECT_THROW(RealizedKernel::estimate({0.01}, KernelType::Bartlett, 2, true),
                 std::invalid_argument);  // n < H+1

    // Part 4: 噪声方差 + bandwidth (已知解析参数验证公式, 非 R 对标)
    // 使用已知 ω², IV, n 验证 BNS 2008 eq.51 公式正确性
    // 注意: NoiseVarianceEstimator 的 ω²=RV/(2n) 估计在高信号/噪声比时
    //       使 ω²/IV ≈ 1/(2n) 极小, 导致 H*→0; 这是估计器的已知局限,
    //       非 bug. 此处用已知参数直接验证公式实现.
    const Real omega2_known = 1e-4;   // ω² = 10⁻⁴
    const Real iv_known     = 1e-2;   // IV = 10⁻²
    const Size  n_known     = 500;
    // H* = c × ξ^(4/5) × (ω²/IV)^(2/5) × n^(3/5),  c=5.74 (Bartlett)
    // ξ² = IV + 2ω² = 0.0102,  ξ = 0.100995
    // H* = 5.74 × 0.100995^0.8 × (0.01)^0.4 × 500^0.6
    //    ≈ 5.74 × 0.1597 × 0.1585 × 41.63 ≈ 6.05
    const Real xi2_known = iv_known + 2.0 * omega2_known;
    const Real xi_known  = std::sqrt(xi2_known);
    const Real ratio_known = omega2_known / iv_known;
    const Real H_star_expected = 5.74 *
        std::pow(xi_known, 0.8) *
        std::pow(ratio_known, 0.4) *
        std::pow(static_cast<Real>(n_known), 0.6);

    Size H_opt = optimal_bandwidth(omega2_known, iv_known, n_known, KernelType::Bartlett);
    EXPECT_EQ(H_opt, static_cast<Size>(std::round(H_star_expected)));
    EXPECT_GT(H_opt, 3u) << "H* should be ~6 for these parameters";

    // 异常处理: 非正参数
    EXPECT_THROW(optimal_bandwidth(0.0, iv_known, n_known, KernelType::Bartlett),
                 std::invalid_argument);
    EXPECT_THROW(optimal_bandwidth(omega2_known, 0.0, n_known, KernelType::Bartlett),
                 std::invalid_argument);
    EXPECT_THROW(optimal_bandwidth(omega2_known, iv_known, 0, KernelType::Bartlett),
                 std::invalid_argument);

    // NoiseVarianceEstimator 基本验证 (纯噪声, ω² = RV/(2n))
    std::normal_distribution<Real> pure_noise_dist(0.0, 0.001);
    std::vector<Real> pure_noise(200);
    for (auto& r : pure_noise) r = pure_noise_dist(gen);
    auto nv = NoiseVarianceEstimator::estimate(pure_noise);
    EXPECT_GT(nv.omega2, 0.0);
    EXPECT_GT(nv.integrated_variance, 0.0);
    EXPECT_EQ(nv.n_obs, 200u);
    // ω² = RV/(2n) ≈ σ²/2 = 5e-7 (对 σ=0.001 纯噪声)
    EXPECT_NEAR(nv.omega2, 5e-7, 3e-7);  // 宽容区间 (随机波动)
}
