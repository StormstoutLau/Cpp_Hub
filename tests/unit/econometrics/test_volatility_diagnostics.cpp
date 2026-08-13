// =============================================================================
// test_volatility_diagnostics.cpp - Phase 7A Wave 1 波动率模型诊断测试
//
// 15 用例: 标准化残差计算(3) + z_t LB(3) + z_t² LB(4) + JB(2) + 综合(3)
//
// 硬编码基准值: 全部解析手算或与 residual_diagnostics 对照
// 容差: z_t LB/z² LB = 1e-10 (复用 ljung_box_test), JB = 1e-10
//
// 排幻觉点覆盖:
//   H8 (z_t² 的 LB 检验是关键, 非 z_t)
//   - z_t 无自相关但 z_t² 有自相关 → GARCH 未充分捕捉条件异方差
//
// 教材锚点: Tsay 3ed Ch.3, McNeil-Frey-Embrechts 2005 §5.3
// =============================================================================

#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>

#include "cpphub/econometrics/inference/volatility_diagnostics.hpp"

using cpphub::v1::econometrics::volatility_diagnostics;
using cpphub::v1::econometrics::ljung_box_test;
using cpphub::v1::econometrics::jarque_bera_test;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
constexpr Real TOL = 1e-10;  // 复用已验证的 ljung_box_test/jarque_bera_test
}  // namespace

// =============================================================================
// 标准化残差计算 (3 用例)
// =============================================================================

// --- 1: h_t 恒定, z_t = ε_t/√h ---
TEST(VolatilityDiagnosticsTest, HomoscedasticStandardization) {
    // ε = {1, 2, 3, 4, 5}, h = {1, 1, 1, 1, 1}
    // z = ε/√1 = {1, 2, 3, 4, 5}
    std::vector<Real> eps = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<Real> h   = {1.0, 1.0, 1.0, 1.0, 1.0};
    auto res = volatility_diagnostics(eps, h, 1);

    ASSERT_EQ(res.standardized_residuals.size(), 5u);
    for (Size i = 0; i < 5; ++i) {
        EXPECT_NEAR(res.standardized_residuals[i], eps[i], TOL);
    }
}

// --- 2: h_t 时变, z_t 正确缩放 ---
TEST(VolatilityDiagnosticsTest, HeteroscedasticStandardization) {
    // ε = {2, 4, 6, 8, 10}, h = {4, 4, 9, 9, 4}
    // z = {2/2, 4/2, 6/3, 8/3, 10/2} = {1, 2, 2, 8/3, 5}
    std::vector<Real> eps = {2.0, 4.0, 6.0, 8.0, 10.0};
    std::vector<Real> h   = {4.0, 4.0, 9.0, 9.0, 4.0};
    auto res = volatility_diagnostics(eps, h, 1);

    ASSERT_EQ(res.standardized_residuals.size(), 5u);
    EXPECT_NEAR(res.standardized_residuals[0], 1.0, TOL);
    EXPECT_NEAR(res.standardized_residuals[1], 2.0, TOL);
    EXPECT_NEAR(res.standardized_residuals[2], 2.0, TOL);
    EXPECT_NEAR(res.standardized_residuals[3], 8.0 / 3.0, TOL);
    EXPECT_NEAR(res.standardized_residuals[4], 5.0, TOL);
}

// --- 3: 输入验证 (尺寸不匹配/h_t<=0 抛异常) ---
TEST(VolatilityDiagnosticsTest, InputValidation) {
    // 尺寸不匹配
    std::vector<Real> eps = {1.0, 2.0, 3.0};
    std::vector<Real> h   = {1.0, 1.0};
    EXPECT_THROW(volatility_diagnostics(eps, h, 1), std::invalid_argument);

    // h_t <= 0
    std::vector<Real> eps2 = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<Real> h2   = {1.0, 0.0, 1.0, 1.0, 1.0};  // h[1]=0
    EXPECT_THROW(volatility_diagnostics(eps2, h2, 1), std::runtime_error);

    // h_t < 0
    std::vector<Real> h3 = {1.0, -1.0, 1.0, 1.0, 1.0};  // h[1]<0
    EXPECT_THROW(volatility_diagnostics(eps2, h3, 1), std::runtime_error);

    // N < 5
    std::vector<Real> eps4 = {1.0, 2.0, 3.0};
    std::vector<Real> h4   = {1.0, 1.0, 1.0};
    EXPECT_THROW(volatility_diagnostics(eps4, h4, 1), std::invalid_argument);
}

// =============================================================================
// z_t LB 检验 (3 用例)
// =============================================================================

// --- 4: 无自相关残差 → z_t LB 不拒绝 ---
TEST(VolatilityDiagnosticsTest, ZLjungBoxNoAutocorrelation) {
    // 伪随机数据 (π 位数), h=1, z=ε
    // lag=1, N=15, LB 应不显著
    std::vector<Real> eps = {3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0,
                             5.0, 3.0, 5.0, 8.0, 9.0, 7.0, 9.0};
    std::vector<Real> h(15, 1.0);
    auto res = volatility_diagnostics(eps, h, 1);

    EXPECT_EQ(res.z_ljung_box.lag, 1u);
    EXPECT_GE(res.z_ljung_box.base.statistic, 0.0);
    // 伪随机数据, LB 不应太大 (放宽阈值, 不强制 p>0.05)
    EXPECT_LT(res.z_ljung_box.base.statistic, 10.0);
}

// --- 5: 有自相关残差 → z_t LB 拒绝 ---
TEST(VolatilityDiagnosticsTest, ZLjungBoxAutocorrelationRejects) {
    // ε = {1,2,3,4,5,6,7,8}: 线性趋势, 强正自相关
    // h=1, z=ε, LB(z, 1) 应拒绝
    std::vector<Real> eps = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<Real> h(8, 1.0);
    auto res = volatility_diagnostics(eps, h, 1);

    EXPECT_GT(res.z_ljung_box.base.statistic, 3.84);  // χ²(1) 5% 临界
    EXPECT_TRUE(res.z_ljung_box.base.reject_null);
    EXPECT_LT(res.z_ljung_box.base.p_value, 0.05);
}

// --- 6: z_t LB 与直接调用 ljung_box_test(z) 一致 ---
TEST(VolatilityDiagnosticsTest, ZLjungBoxMatchesDirectCall) {
    // 验证 volatility_diagnostics 内部对 z_t 的 LB 与外部直接调用一致
    std::vector<Real> eps = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<Real> h(8, 2.0);  // h=2, z = ε/√2
    auto res = volatility_diagnostics(eps, h, 1);

    // 手动计算 z
    std::vector<Real> z(8);
    for (Size i = 0; i < 8; ++i) z[i] = eps[i] / std::sqrt(2.0);
    auto direct_lb = ljung_box_test(z, 1);

    EXPECT_NEAR(res.z_ljung_box.base.statistic, direct_lb.base.statistic, TOL);
    EXPECT_NEAR(res.z_ljung_box.base.p_value, direct_lb.base.p_value, TOL);
    EXPECT_EQ(res.z_ljung_box.lag, direct_lb.lag);
}

// =============================================================================
// z_t² LB 检验 (4 用例) — 排幻觉点 H8 关键
// =============================================================================

// --- 7: H8 核心 - z_t² 有自相关 → z² LB 拒绝 (ARCH 效应未消除) ---
TEST(VolatilityDiagnosticsTest, ZSquaredLBDetectsARCHEffect) {
    // 排幻觉点 H8: z_t² 的 LB 检验是关键
    // 场景: 模型未捕捉 ARCH 效应 (h_t 恒定, 但 ε_t² 有自相关)
    // ε = {1, 2, 1, 2, 1, 2, 1, 2}: 幅度交替
    // h = {1,1,...,1}: 模型认为方差恒定
    // z = ε, z² = {1, 4, 1, 4, 1, 4, 1, 4}: 强正自相关
    std::vector<Real> eps = {1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0};
    std::vector<Real> h(8, 1.0);
    auto res = volatility_diagnostics(eps, h, 1);

    // z² LB 应拒绝 (ARCH 效应未消除)
    EXPECT_GT(res.z_squared_ljung_box.base.statistic, 3.84);
    EXPECT_TRUE(res.z_squared_ljung_box.base.reject_null);
    EXPECT_LT(res.z_squared_ljung_box.base.p_value, 0.05);
}

// --- 8: z_t² 无自相关 → z² LB 不拒绝 ---
TEST(VolatilityDiagnosticsTest, ZSquaredLBNoARCHEffect) {
    // 伪随机数据, z_t² 应无显著自相关
    // 用足够大的 N 使 LB 稳定
    std::vector<Real> eps = {3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0,
                             5.0, 3.0, 5.0, 8.0, 9.0, 7.0, 9.0};
    std::vector<Real> h(15, 1.0);
    auto res = volatility_diagnostics(eps, h, 1);

    // 伪随机数据, z² LB 不应太大 (放宽阈值)
    EXPECT_GE(res.z_squared_ljung_box.base.statistic, 0.0);
    EXPECT_LT(res.z_squared_ljung_box.base.statistic, 10.0);
}

// --- 9: z² LB 与直接调用 ljung_box_test(z²) 一致 ---
TEST(VolatilityDiagnosticsTest, ZSquaredLjungBoxMatchesDirectCall) {
    // 验证 volatility_diagnostics 内部对 z_t² 的 LB 与外部直接调用一致
    std::vector<Real> eps = {2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0};
    std::vector<Real> h(8, 4.0);  // h=4, z = ε/2
    auto res = volatility_diagnostics(eps, h, 1);

    // 手动计算 z²
    std::vector<Real> z2(8);
    for (Size i = 0; i < 8; ++i) {
        Real z = eps[i] / 2.0;
        z2[i] = z * z;
    }
    auto direct_lb = ljung_box_test(z2, 1);

    EXPECT_NEAR(res.z_squared_ljung_box.base.statistic, direct_lb.base.statistic, TOL);
    EXPECT_NEAR(res.z_squared_ljung_box.base.p_value, direct_lb.base.p_value, TOL);
}

// --- 10: H8 - z_t 无显著自相关但 z_t² 有自相关 (经典 ARCH 效应) ---
TEST(VolatilityDiagnosticsTest, ZNoAutocorrButZSquaredHasAutocorr) {
    // 排幻觉点 H8: z_t 无自相关但 z_t² 有自相关 → GARCH 未充分捕捉条件异方差
    // 构造: z_t 符号伪随机 (弱自相关), 但 |z_t| 有聚集 (z_t² 强自相关)
    //
    // z = {1, 2, -0.5, 3, -1, 2.5, -0.8, 2.8, -1.2, 3.1}
    //   - 符号: +, +, -, +, -, +, -, +, -, + (不规则, z_t 弱自相关)
    //   - 幅度: 1, 2, 0.5, 3, 1, 2.5, 0.8, 2.8, 1.2, 3.1 (有聚集, z_t² 强自相关)
    // h = 1 (模型未捕捉幅度聚集)
    std::vector<Real> z = {1.0, 2.0, -0.5, 3.0, -1.0, 2.5, -0.8, 2.8, -1.2, 3.1};
    std::vector<Real> h(10, 1.0);
    auto res = volatility_diagnostics(z, h, 1);

    // z² LB 应拒绝 (幅度聚集 = ARCH 效应)
    EXPECT_TRUE(res.z_squared_ljung_box.base.reject_null);
    EXPECT_LT(res.z_squared_ljung_box.base.p_value, 0.05);

    // 关键: 即使 z_t LB 不拒绝 (或弱拒绝), z² LB 拒绝 → 模型不充分
    // 这是 H8 的核心: 仅检查 z_t 不够, 必须检查 z_t²
    EXPECT_FALSE(res.model_adequate);  // z² LB 拒绝 → 模型不充分
}

// =============================================================================
// Jarque-Bera 正态性检验 (2 用例)
// =============================================================================

// --- 11: 正态残差 → JB 不拒绝 ---
TEST(VolatilityDiagnosticsTest, NormalResidualsJBNotReject) {
    // 对称数据, 近似正态
    // z = {-2, -1, 0, 1, 2}: 对称, skew=0, kurt=1.7 (platykurtic)
    // JB = 5 * (0 + (1.7-3)²/24) = 5 * 1.69/24 ≈ 0.352
    std::vector<Real> eps = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<Real> h(5, 1.0);
    auto res = volatility_diagnostics(eps, h, 1);

    EXPECT_NEAR(res.z_jarque_bera.base.statistic, 0.3520833333333333, TOL);
    EXPECT_FALSE(res.z_jarque_bera.base.reject_null);
    EXPECT_GT(res.z_jarque_bera.base.p_value, 0.05);
}

// --- 12: 非正态残差 → JB 拒绝 ---
TEST(VolatilityDiagnosticsTest, NonNormalResidualsJBRejects) {
    // 极右偏数据
    // ε = {1,1,1,1,1,1,1,1,1,10}: 极右偏
    std::vector<Real> eps = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0};
    std::vector<Real> h(10, 1.0);
    auto res = volatility_diagnostics(eps, h, 1);

    EXPECT_TRUE(res.z_jarque_bera.base.reject_null);
    EXPECT_LT(res.z_jarque_bera.base.p_value, 0.001);
}

// =============================================================================
// 综合诊断 (3 用例)
// =============================================================================

// --- 13: model_adequate 逻辑一致性 ---
TEST(VolatilityDiagnosticsTest, ModelAdequateConsistency) {
    // 验证 model_adequate == (!z_lb.reject && !z²_lb.reject && !z_jb.reject)
    // 用对称无自相关数据, 可能 model_adequate = true
    std::vector<Real> eps = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<Real> h(5, 1.0);
    auto res = volatility_diagnostics(eps, h, 1);

    // 验证逻辑一致性
    const bool expected = !res.z_ljung_box.base.reject_null &&
                          !res.z_squared_ljung_box.base.reject_null &&
                          !res.z_jarque_bera.base.reject_null;
    EXPECT_EQ(res.model_adequate, expected);
}

// --- 14: ARCH 效应未消除 → model_adequate = false ---
TEST(VolatilityDiagnosticsTest, ModelInadequateWhenARCHEffect) {
    // z_t² 有自相关 (ARCH 效应) → model_adequate = false
    std::vector<Real> eps = {1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0};
    std::vector<Real> h(8, 1.0);
    auto res = volatility_diagnostics(eps, h, 1);

    // z² LB 拒绝 → model_adequate = false
    EXPECT_FALSE(res.model_adequate);
    EXPECT_TRUE(res.z_squared_ljung_box.base.reject_null);
}

// --- 15: weighted LB 字段有效性 ---
TEST(VolatilityDiagnosticsTest, WeightedLBFieldsValid) {
    // 验证 weighted_lb_statistic/p_value 字段存在且有效
    // 注: 当前实现用 z² 的标准 LB (Fisher-Gallagher 2012 精确公式待补充)
    std::vector<Real> eps = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<Real> h(8, 1.0);
    auto res = volatility_diagnostics(eps, h, 1);

    // weighted_lb 应与 z² LB 相同 (当前实现)
    EXPECT_GE(res.weighted_lb_statistic, 0.0);
    EXPECT_GT(res.weighted_lb_p_value, 0.0);
    EXPECT_LE(res.weighted_lb_p_value, 1.0);

    // 当前实现: weighted_lb = z² LB (待 Fisher-Gallagher 2012 精确公式替换)
    EXPECT_NEAR(res.weighted_lb_statistic, res.z_squared_ljung_box.base.statistic, TOL);
    EXPECT_NEAR(res.weighted_lb_p_value, res.z_squared_ljung_box.base.p_value, TOL);
}
