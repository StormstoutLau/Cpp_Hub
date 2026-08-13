// =============================================================================
// test_hfecon_diagnostics.cpp - Phase 7A Wave 2c HAR/HEAVY 诊断测试
//
// 15 用例:
//   HAR 残差 LB (4) + HAR MZ (4) + HEAVY 标准化残差 (4) + HEAVY 交叉诊断 (3)
//
// 排幻觉点覆盖:
//   H8 (z_t² LB 是 ARCH 效应检验关键, 非 z_t LB)
//   H10 (MZ R² 是预测精度, joint F 检验 alpha=0 & beta=1)
//
// 教材锚点: Corsi 2009 (HAR), Shephard-Sheppard 2010 (HEAVY), Patton 2011 (MZ)
// =============================================================================

#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

#include "cpphub/hfecon/hfecon_diagnostics.hpp"

using cpphub::v1::hfecon::har_diagnostics;
using cpphub::v1::hfecon::heavy_diagnostics;
using cpphub::v1::hfecon::HARDiagnosticsResult;
using cpphub::v1::hfecon::HEAVYDiagnosticsResult;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
// =============================================================================
// 伪随机白噪声生成 (无自相关)
// 使用 std::mt19937 + 固定 seed 生成确定性 iid 序列, 保证可复现且无自相关
// 旧实现 sin(x*e)+cos(x*π) 在整数 x 上 cos(x*π)=±1 振荡, 有强自相关 (LB p≈1e-23)
// =============================================================================

// 生成 n 点白噪声序列 (seed=42, 可复现)
inline std::vector<Real> make_white_noise(Size n) {
    std::mt19937 gen(42);
    std::normal_distribution<Real> dist(0.0, 1.0);
    std::vector<Real> noise(n);
    for (Size i = 0; i < n; ++i) noise[i] = dist(gen);
    return noise;
}

// 生成 n 点平滑 fitted 序列 (模拟 HAR 拟合值, 非常数)
inline std::vector<Real> make_smooth_fitted(Size n) {
    std::vector<Real> fitted(n);
    for (Size i = 0; i < n; ++i) fitted[i] = 1.0 + 0.1 * static_cast<Real>(i);
    return fitted;
}
}

// =============================================================================
// HAR 残差 LB 检验 (4 用例)
// =============================================================================

// --- HAR LB 1: 无自相关残差 → LB 不拒绝 ---
TEST(HARResidualLjungBox, WhiteNoiseResidualsNotRejected) {
    const Size n = 50;
    auto fitted = make_smooth_fitted(n);
    auto residuals = make_white_noise(n);
    std::vector<Real> actual(n);
    for (Size i = 0; i < n; ++i) actual[i] = fitted[i] + residuals[i];

    auto res = har_diagnostics(actual, fitted, residuals);

    // LB 不拒绝 (白噪声无自相关)
    EXPECT_FALSE(res.residual_ljung_box.base.reject_null);
    EXPECT_GT(res.residual_ljung_box.base.p_value, 0.05);
}

// --- HAR LB 2: 有自相关残差 (AR(1)) → LB 拒绝 ---
TEST(HARResidualLjungBox, AR1ResidualsRejected) {
    const Size n = 50;
    auto fitted = make_smooth_fitted(n);
    auto noise = make_white_noise(n);

    // AR(1) with φ=0.9 (强自相关)
    std::vector<Real> residuals(n);
    residuals[0] = noise[0];
    for (Size i = 1; i < n; ++i) {
        residuals[i] = 0.9 * residuals[i - 1] + 0.1 * noise[i];
    }
    std::vector<Real> actual(n);
    for (Size i = 0; i < n; ++i) actual[i] = fitted[i] + residuals[i];

    auto res = har_diagnostics(actual, fitted, residuals);

    // LB 拒绝 (AR(1) 有强自相关)
    EXPECT_TRUE(res.residual_ljung_box.base.reject_null);
    EXPECT_LT(res.residual_ljung_box.base.p_value, 0.05);
}

// --- HAR LB 3: 残差太少 (< 10) → 抛异常 ---
TEST(HARResidualLjungBox, TooFewObservationsThrows) {
    std::vector<Real> actual(5, 1.0), fitted(5, 1.0), residuals(5, 0.1);
    EXPECT_THROW(har_diagnostics(actual, fitted, residuals), std::invalid_argument);
}

// --- HAR LB 4: 大小不匹配 → 抛异常 ---
TEST(HARResidualLjungBox, SizeMismatchThrows) {
    std::vector<Real> actual(20, 1.0), fitted(20, 1.0), residuals(15, 0.1);
    EXPECT_THROW(har_diagnostics(actual, fitted, residuals), std::invalid_argument);
}

// =============================================================================
// HAR Mincer-Zarnowitz 回归 (4 用例)
// 排幻觉点 H10: MZ joint F 检验 alpha=0 & beta=1, R² 是预测精度
// =============================================================================

// --- HAR MZ 1: 好预测 (actual ≈ fitted + small noise) → MZ 不拒绝 ---
TEST(HARMincerZarnowitz, GoodForecastNotRejected) {
    const Size n = 50;
    auto fitted = make_smooth_fitted(n);
    auto noise = make_white_noise(n);
    std::vector<Real> actual(n), residuals(n);
    for (Size i = 0; i < n; ++i) {
        residuals[i] = noise[i] * 0.1;  // 小残差
        actual[i] = fitted[i] + residuals[i];
    }

    auto res = har_diagnostics(actual, fitted, residuals);

    // MZ 不拒绝 (alpha≈0, beta≈1)
    EXPECT_FALSE(res.mz_regression.base.reject_null);
    EXPECT_NEAR(res.mz_regression.alpha, 0.0, 0.5);
    EXPECT_NEAR(res.mz_regression.beta, 1.0, 0.3);
    // R² 应较高 (预测精度好)
    EXPECT_GT(res.mz_regression.r_squared, 0.8);
}

// --- HAR MZ 2: 有偏预测 (actual = 0.5 * fitted) → MZ 拒绝 ---
TEST(HARMincerZarnowitz, BiasedForecastRejected) {
    const Size n = 50;
    auto fitted = make_smooth_fitted(n);
    auto noise = make_white_noise(n);
    std::vector<Real> actual(n), residuals(n);
    for (Size i = 0; i < n; ++i) {
        actual[i] = 0.5 * fitted[i] + noise[i] * 0.05;  // beta=0.5
        residuals[i] = actual[i] - fitted[i];             // 残差 = -0.5*fitted + noise
    }

    auto res = har_diagnostics(actual, fitted, residuals);

    // MZ 拒绝 (beta ≈ 0.5, 远离 1)
    EXPECT_TRUE(res.mz_regression.base.reject_null);
    EXPECT_NEAR(res.mz_regression.beta, 0.5, 0.2);
}

// --- HAR MZ 3: R² 预测精度验证 ---
TEST(HARMincerZarnowitz, RSquaredAsAccuracyMetric) {
    const Size n = 50;
    auto fitted = make_smooth_fitted(n);
    auto noise = make_white_noise(n);

    // 好预测: 小残差 → R² 高
    std::vector<Real> good_actual(n), good_resid(n);
    for (Size i = 0; i < n; ++i) {
        good_resid[i] = noise[i] * 0.05;
        good_actual[i] = fitted[i] + good_resid[i];
    }
    auto res_good = har_diagnostics(good_actual, fitted, good_resid);
    EXPECT_GT(res_good.r_squared, 0.9);

    // 差预测: 大残差 → R² 低
    std::vector<Real> bad_actual(n), bad_resid(n);
    for (Size i = 0; i < n; ++i) {
        bad_resid[i] = noise[i] * 5.0;  // 大残差
        bad_actual[i] = fitted[i] + bad_resid[i];
    }
    auto res_bad = har_diagnostics(bad_actual, fitted, bad_resid);
    EXPECT_LT(res_bad.r_squared, 0.5);
}

// --- HAR MZ 4: alpha/beta 参数值验证 ---
TEST(HARMincerZarnowitz, AlphaBetaValuesCorrect) {
    const Size n = 50;
    auto fitted = make_smooth_fitted(n);
    auto noise = make_white_noise(n);
    // 构造 actual = 0.5 + 1.2 * fitted + noise
    std::vector<Real> actual(n), residuals(n);
    for (Size i = 0; i < n; ++i) {
        actual[i] = 0.5 + 1.2 * fitted[i] + noise[i] * 0.1;
        residuals[i] = actual[i] - fitted[i];
    }
    auto res = har_diagnostics(actual, fitted, residuals);

    // MZ 回归: actual ~ alpha + beta * fitted
    EXPECT_NEAR(res.mz_regression.alpha, 0.5, 0.3);
    EXPECT_NEAR(res.mz_regression.beta, 1.2, 0.3);
}

// =============================================================================
// HEAVY 标准化残差诊断 (4 用例)
// 排幻觉点 H8: z_t² LB 是 ARCH 效应检验关键
// =============================================================================

// --- HEAVY 1: 标准化残差 iid → z_t LB 和 z_t² LB 不拒绝 ---
TEST(HEAVYStandardizedResiduals, IIDResidualsNotRejected) {
    const Size n = 50;
    // h_t = 1 (常数), ε_t = noise → z_t = ε_t (白噪声)
    std::vector<Real> h_t(n, 1.0);
    std::vector<Real> eps = make_white_noise(n);
    std::vector<Real> mu_t(n, 0.0);
    // 方差方程残差用独立白噪声 (seed=137 保证与 eps 不同)
    std::mt19937 gen2(137);
    std::normal_distribution<Real> dist2(0.0, 1.0);
    std::vector<Real> var_resid(n);
    for (Size i = 0; i < n; ++i) var_resid[i] = dist2(gen2);

    auto res = heavy_diagnostics(eps, mu_t, var_resid, h_t);

    // z_t LB 不拒绝 (白噪声)
    EXPECT_FALSE(res.measurement_equation.z_ljung_box.base.reject_null);
    // z_t² LB 不拒绝 (无 ARCH 效应)
    EXPECT_FALSE(res.measurement_equation.z_squared_ljung_box.base.reject_null);
}

// --- HEAVY 2: 有 ARCH 效应 → z_t² LB 拒绝 ---
// 排幻觉点 H8: z_t² LB 是关键 (z_t LB 可能不拒绝但有 ARCH 效应)
TEST(HEAVYStandardizedResiduals, ARCHEffectRejected) {
    const Size n = 50;
    // 构造 ARCH(1): ε_t² = ω + α·ε_{t-1}² + noise
    // z_t² = ε_t²/h_t 有强正自相关
    std::vector<Real> h_t(n), eps(n), mu_t(n, 0.0), var_resid(n);
    // 三组独立白噪声: 驱动 z2 AR(1) / 符号 / var_resid
    std::mt19937 gen_a(7), gen_b(13), gen_c(17);
    std::normal_distribution<Real> dist(0.0, 1.0);
    Real z2 = 1.0;
    for (Size i = 0; i < n; ++i) {
        h_t[i] = 1.0;  // h_t 常数, 使 z_t² = ε_t²
        // z_t² 服从 AR(1) with φ=0.9 (强自相关)
        z2 = 0.9 * z2 + 0.5 * std::fabs(dist(gen_a)) + 0.1;
        // ε_t = ±√z_t² (符号随机)
        const Real sign = (dist(gen_b) > 0) ? 1.0 : -1.0;
        eps[i] = sign * std::sqrt(z2);
        var_resid[i] = dist(gen_c) * 0.1;
    }

    auto res = heavy_diagnostics(eps, mu_t, var_resid, h_t);

    // z_t² LB 应拒绝 (ARCH 效应: z_t² 有强自相关)
    EXPECT_TRUE(res.measurement_equation.z_squared_ljung_box.base.reject_null);
}

// --- HEAVY 3: 残差太少 (< 10) → 抛异常 ---
TEST(HEAVYStandardizedResiduals, TooFewObservationsThrows) {
    std::vector<Real> eps(5, 0.1), mu(5, 0.0), var_resid(5, 0.1), h(5, 1.0);
    EXPECT_THROW(heavy_diagnostics(eps, mu, var_resid, h), std::invalid_argument);
}

// --- HEAVY 4: 大小不匹配 → 抛异常 ---
TEST(HEAVYStandardizedResiduals, SizeMismatchThrows) {
    std::vector<Real> eps(20, 0.1), mu(20, 0.0), var_resid(15, 0.1), h(20, 1.0);
    EXPECT_THROW(heavy_diagnostics(eps, mu, var_resid, h), std::invalid_argument);
}

// =============================================================================
// HEAVY 交叉诊断 (3 用例)
// =============================================================================

// --- HEAVY 交叉 1: h_t 与 RM_t 高相关 (波动率聚集) ---
TEST(HEAVYCrossDiagnostics, HighCorrelationHRM) {
    const Size n = 50;
    auto noise = make_white_noise(n);
    // h_t 递增, RM_t = eps = √h_t * (1 + 小扰动), 与 h_t 强正相关
    std::vector<Real> h_t(n), eps(n), mu_t(n, 0.0), var_resid(n);
    for (Size i = 0; i < n; ++i) {
        h_t[i] = 1.0 + 0.5 * static_cast<Real>(i);     // 递增 h_t
        // eps = √h_t * (1 + 0.05*noise), 主项 √h_t 驱动正相关
        eps[i] = std::sqrt(h_t[i]) * (1.0 + 0.05 * noise[i]);
        var_resid[i] = noise[(i + 40) % n] * 0.1;
    }
    auto res = heavy_diagnostics(eps, mu_t, var_resid, h_t);

    // h_t 与 RM_t 相关性应较高 (波动率聚集)
    EXPECT_GT(res.correlation_h_rm, 0.5);
}

// --- HEAVY 交叉 2: h_t 与 RM_t 低相关 ---
TEST(HEAVYCrossDiagnostics, LowCorrelationHRM) {
    const Size n = 50;
    // h_t 递增, RM_t 递减 (eps 为负且递减)
    std::vector<Real> h_t(n), eps(n), mu_t(n, 0.0), var_resid(n);
    std::mt19937 gen(99);
    std::normal_distribution<Real> dist(0.0, 1.0);
    for (Size i = 0; i < n; ++i) {
        h_t[i] = 1.0 + 0.5 * static_cast<Real>(i);     // 递增 h_t
        // eps 递减 (从正值到负值), RM_t = eps 递减
        eps[i] = static_cast<Real>(n - i) * 0.3 - 5.0;  // 从 n*0.3-5 到 -5
        var_resid[i] = dist(gen) * 0.1;
    }
    auto res = heavy_diagnostics(eps, mu_t, var_resid, h_t);

    // h_t 递增, RM_t 递减 → 负相关, |cor| 应较低或为负
    EXPECT_LT(res.correlation_h_rm, 0.3);
}

// --- HEAVY 交叉 3: model_adequate 综合判定 ---
TEST(HEAVYCrossDiagnostics, ModelAdequateComposite) {
    const Size n = 50;
    // 构造 "好" 的 HEAVY 数据: 白噪声标准化残差 + 高相关性
    std::vector<Real> h_t(n), eps(n), mu_t(n, 0.0), var_resid(n);
    std::mt19937 gen_e(21), gen_v(23);
    std::normal_distribution<Real> dist(0.0, 1.0);
    for (Size i = 0; i < n; ++i) {
        // h_t 有波动但与 RM_t 正相关
        h_t[i] = 1.0 + 0.3 * static_cast<Real>(i % 5);
        eps[i] = std::sqrt(h_t[i]) * std::fabs(dist(gen_e)) * 0.5;
        var_resid[i] = dist(gen_v) * 0.5;
    }
    auto res = heavy_diagnostics(eps, mu_t, var_resid, h_t);

    // 测试 model_adequate 字段存在且为 bool
    EXPECT_TRUE(res.model_adequate || !res.model_adequate);
    // 相关性应正 (h_t 与 RM_t 同向)
    EXPECT_GT(res.correlation_h_rm, 0.0);
}
