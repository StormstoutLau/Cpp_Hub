// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.2 任务 2.4 - 信息准则测试
// 验证方法: 手算解析值 + statsmodels 对照 (容差 1e-10)
// 排幻觉点 D1: AIC = 2K - 2ℓ (非 -2ℓ + 2K, 符号一致)
// 排幻觉点 D2: BIC = K·log(N) - 2ℓ (log 为自然对数, 非 log10)
// 排幻觉点 D3: HQ = 2K·log(log(N)) - 2ℓ (Hannan-Quinn 1979, 非 K·log(log(N)))
// 排幻觉点 D4: AICc = AIC + 2K(K+1)/(N-K-1) (Sugiura 1978, 当 N<=K+1 未定义)
#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>

#include "cpphub/econometrics/inference/diagnostics.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;

// =============================================================================
// 手算验证 (手算期望值)
// =============================================================================

// Case 1: ℓ = -10, K = 2, N = 100
// AIC = 2*2 - 2*(-10) = 4 + 20 = 24
// BIC = 2*log(100) - 2*(-10) = 2*4.605170186 + 20 = 29.210340372
// HQ  = 2*2*log(log(100)) - 2*(-10) = 4*log(4.605170186) + 20 = 4*1.527179817 + 20 = 26.108719268
// AICc = 24 + 2*2*3/(100-2-1) = 24 + 12/97 = 24.123711340
TEST(DiagnosticsTest, AIC_BIC_HQ_Case1_HandComputed) {
    const Real ll = -10.0;
    const Size K = 2;
    const Size N = 100;
    const InformationCriteria ic = compute_information_criteria(ll, K, N);

    EXPECT_NEAR(ic.aic, 24.0, 1e-10);
    EXPECT_NEAR(ic.bic, 2.0 * std::log(100.0) + 20.0, 1e-10);
    EXPECT_NEAR(ic.hq, 4.0 * std::log(std::log(100.0)) + 20.0, 1e-10);
    EXPECT_NEAR(ic.aicc, 24.0 + 12.0 / 97.0, 1e-10);
}

// Case 2: ℓ = -50.5, K = 3, N = 30
// AIC = 2*3 - 2*(-50.5) = 6 + 101 = 107
// BIC = 3*log(30) + 101 = 3*3.401197382 + 101 = 111.203592146
// HQ  = 2*3*log(log(30)) + 101 = 6*log(3.401197382) + 101 = 6*1.224191288 + 101 = 108.345147728
// AICc = 107 + 2*3*4/(30-3-1) = 107 + 24/26 = 107.923076923
TEST(DiagnosticsTest, AIC_BIC_HQ_Case2_HandComputed) {
    const Real ll = -50.5;
    const Size K = 3;
    const Size N = 30;
    const InformationCriteria ic = compute_information_criteria(ll, K, N);

    EXPECT_NEAR(ic.aic, 107.0, 1e-10);
    EXPECT_NEAR(ic.bic, 3.0 * std::log(30.0) + 101.0, 1e-10);
    EXPECT_NEAR(ic.hq, 6.0 * std::log(std::log(30.0)) + 101.0, 1e-10);
    EXPECT_NEAR(ic.aicc, 107.0 + 24.0 / 26.0, 1e-10);
}

// Case 3: ℓ = 0, K = 1, N = 10 (边界: 单参数)
// AIC = 2*1 - 0 = 2
// BIC = 1*log(10) - 0 = 2.302585093
// HQ  = 2*1*log(log(10)) - 0 = 2*log(2.302585093) = 2*0.834032445 = 1.668064890
// AICc = 2 + 2*1*2/(10-1-1) = 2 + 4/8 = 2.5
TEST(DiagnosticsTest, AIC_BIC_HQ_Case3_SingleParam) {
    const Real ll = 0.0;
    const Size K = 1;
    const Size N = 10;
    const InformationCriteria ic = compute_information_criteria(ll, K, N);

    EXPECT_NEAR(ic.aic, 2.0, 1e-10);
    EXPECT_NEAR(ic.bic, std::log(10.0), 1e-10);
    EXPECT_NEAR(ic.hq, 2.0 * std::log(std::log(10.0)), 1e-10);
    EXPECT_NEAR(ic.aicc, 2.5, 1e-10);
}

// =============================================================================
// 公式验证 (排幻觉点)
// =============================================================================

// 排幻觉点 D1: AIC = 2K - 2ℓ (非 -2ℓ + 2K, 符号一致但验证方向)
// 当 ℓ 增大 (似然更高), AIC 减小 (模型更好)
TEST(DiagnosticsTest, AIC_Decreases_When_LogLikelihood_Increases) {
    const Size K = 2;
    const Size N = 50;
    const InformationCriteria ic1 = compute_information_criteria(-20.0, K, N);
    const InformationCriteria ic2 = compute_information_criteria(-10.0, K, N);
    EXPECT_LT(ic2.aic, ic1.aic);
    EXPECT_LT(ic2.bic, ic1.bic);
    EXPECT_LT(ic2.hq, ic1.hq);
}

// 排幻觉点 D2: BIC 对 N 的惩罚比对 K 更重 (log(N) > 1 当 N > e)
// 大样本下 BIC > AIC (当 N 充分大)
TEST(DiagnosticsTest, BIC_GreaterThan_AIC_For_LargeN) {
    const Real ll = -10.0;
    const Size K = 2;
    const Size N = 100;  // log(100) ≈ 4.6 > 2 = K 系数
    const InformationCriteria ic = compute_information_criteria(ll, K, N);
    EXPECT_GT(ic.bic, ic.aic);
}

// 排幻觉点 D3: HQ 介于 AIC 和 BIC 之间 (Hannan-Quinn 1979 性质)
// 当 N > e^e ≈ 15.15 时, log(log(N)) > 1, HQ 惩罚介于 AIC (系数 2) 和 BIC (系数 log(N))
TEST(DiagnosticsTest, HQ_Between_AIC_and_BIC) {
    const Real ll = -15.0;
    const Size K = 3;
    const Size N = 100;
    const InformationCriteria ic = compute_information_criteria(ll, K, N);
    // AIC 系数 = 2K, HQ 系数 = 2K·log(log(N)), BIC 系数 = K·log(N)
    // 当 N=100: 2 < 2·log(log(100)) ≈ 3.05 < log(100) ≈ 4.6
    // 所以 AIC < HQ < BIC
    EXPECT_LT(ic.aic, ic.hq);
    EXPECT_LT(ic.hq, ic.bic);
}

// 排幻觉点 D4: AICc 小样本修正 - 当 N 很小时 AICc > AIC
TEST(DiagnosticsTest, AICc_Corrects_SmallSample) {
    const Real ll = -5.0;
    const Size K = 2;
    const Size N = 10;  // 小样本
    const InformationCriteria ic = compute_information_criteria(ll, K, N);
    EXPECT_GT(ic.aicc, ic.aic);
    // AICc = AIC + 2K(K+1)/(N-K-1) = AIC + 12/7
    EXPECT_NEAR(ic.aicc - ic.aic, 12.0 / 7.0, 1e-10);
}

// AICc 退化为 AIC 当 N → ∞
TEST(DiagnosticsTest, AICc_ConvergesTo_AIC_For_LargeN) {
    const Real ll = -10.0;
    const Size K = 2;
    const Size N = 1000000;  // 大样本
    const InformationCriteria ic = compute_information_criteria(ll, K, N);
    // 修正项 = 2K(K+1)/(N-K-1) = 12/999997 ≈ 1.2e-5, 趋近 0
    EXPECT_NEAR(ic.aicc, ic.aic, 1e-4);  // 修正项 < 1e-4
}

// =============================================================================
// 边界与异常
// =============================================================================

TEST(DiagnosticsTest, ThrowsOnZeroObs) {
    EXPECT_THROW(compute_information_criteria(0.0, 1, 0), std::invalid_argument);
}

TEST(DiagnosticsTest, ThrowsOnZeroParams) {
    EXPECT_THROW(compute_information_criteria(0.0, 0, 10), std::invalid_argument);
}

// AICc 当 N = K+1 (边界: 分母为零, 应返回 AIC 不抛异常)
TEST(DiagnosticsTest, AICc_Undefined_When_N_Equals_K_Plus_1) {
    const Real ll = -5.0;
    const Size K = 3;
    const Size N = 4;  // N = K + 1, 分母 N-K-1 = 0
    const InformationCriteria ic = compute_information_criteria(ll, K, N);
    // 不抛异常, 返回 AIC
    EXPECT_NEAR(ic.aicc, ic.aic, 1e-10);
}

// =============================================================================
// statsmodels 对照 (容差 1e-10)
// =============================================================================

// statsmodels OLS: ℓ = -10, K = 2, N = 100
// statsmodels.aic = -2*ll + 2*K = 20 + 4 = 24
// statsmodels.bic = -2*ll + K*log(N) = 20 + 2*4.605170186 = 29.210340372
TEST(DiagnosticsTest, MatchStatsmodels_Convention) {
    // statsmodels 使用 aic = -2*ll + 2*K (与 2*K - 2*ll 等价)
    const Real ll = -10.0;
    const Size K = 2;
    const Size N = 100;
    const InformationCriteria ic = compute_information_criteria(ll, K, N);

    const Real sm_aic = -2.0 * ll + 2.0 * static_cast<Real>(K);
    const Real sm_bic = -2.0 * ll + static_cast<Real>(K) * std::log(static_cast<Real>(N));

    EXPECT_NEAR(ic.aic, sm_aic, 1e-10);
    EXPECT_NEAR(ic.bic, sm_bic, 1e-10);
}
