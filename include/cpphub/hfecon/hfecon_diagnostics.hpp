// =============================================================================
// hfecon_diagnostics.hpp - HAR/HEAVY 预测诊断 (P0-P1)
//
// Phase 7A Wave 2c: 高频计量模型诊断
//
// 教材锚点:
//   - Corsi 2009 (HAR 模型)
//   - Shephard-Sheppard 2010 (HEAVY 模型)
//   - Patton 2011 (RV 预测评估, MZ 回归)
//
// ADR-015 决策点: hfecon_diagnostics 调用 econometrics/inference/ 通用诊断
//   (residual_diagnostics / volatility_diagnostics / specification_tests),
//   通用诊断无 Eigen3 依赖, 不污染 hfecon 模块 (ADR-015 方案 B 保证)
//
// 包含 2 个诊断函数:
//   1. har_diagnostics: HAR 模型残差 LB + JB + MZ 回归 + R²
//   2. heavy_diagnostics: HEAVY 模型双方程标准化残差诊断 + h_t/RM_t 相关性
//
// 排幻觉点:
//   H8 (z_t² LB 是 ARCH 效应检验的关键, 非 z_t LB)
//   H10 (MZ R² 是预测精度指标, joint F 检验 alpha=0 & beta=1)
//
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训)
// =============================================================================
#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/inference/residual_diagnostics.hpp"  // LjungBoxResult, JarqueBeraResult
#include "cpphub/econometrics/inference/specification_tests.hpp"    // MincerZarnowitzResult
#include "cpphub/econometrics/inference/volatility_diagnostics.hpp" // VolatilityDiagnosticsResult

namespace cpphub {
inline namespace v1 {
namespace hfecon {

using econometrics::JarqueBeraResult;
using econometrics::LjungBoxResult;
using econometrics::MincerZarnowitzResult;
using econometrics::VolatilityDiagnosticsResult;

// =============================================================================
// HARDiagnosticsResult - HAR 模型诊断结果
//
// HAR (Heterogeneous Autoregressive) 模型假设:
//   1. 残差白噪声 (无自相关) → Ljung-Box 不拒绝
//   2. 残差正态 (JB 检验, 渐近推断前提) → Jarque-Bera 不拒绝
//   3. 预测无偏 (MZ: alpha=0, beta=1) → MZ joint F 不拒绝
//
// 排幻觉点 H10: MZ R² 是预测精度指标 (越高越好),
//               joint F 检验 alpha=0 & beta=1 是无偏性检验
// =============================================================================
struct HARDiagnosticsResult {
    // 残差诊断
    LjungBoxResult residual_ljung_box;  // 残差自相关 (应不拒绝)
    JarqueBeraResult residual_jb;       // 残差正态性 (HAR 假设, 应不拒绝)

    // 预测精度
    MincerZarnowitzResult mz_regression;  // 预测无偏性 (alpha=0 & beta=1, 应不拒绝)

    // 拟合优度
    Real r_squared;
    Real adjusted_r_squared;

    bool model_adequate;  // 所有检验均不拒绝 → true
};

// =============================================================================
// har_diagnostics - HAR 模型诊断
//
// 原理: HAR 残差应无自相关 (LB), 预测应无偏 (MZ: alpha=0, beta=1)
//
// @param actual_rv 实际 RV (观测值)
// @param fitted_rv HAR 拟合值 (预测值)
// @param residuals HAR 残差 (actual_rv - fitted_rv)
// @return HARDiagnosticsResult
// =============================================================================
inline HARDiagnosticsResult har_diagnostics(
    const std::vector<Real>& actual_rv,
    const std::vector<Real>& fitted_rv,
    const std::vector<Real>& residuals) {

    const Size n = residuals.size();
    if (n < 10) {
        throw std::invalid_argument("har_diagnostics: need at least 10 observations");
    }
    if (actual_rv.size() != n || fitted_rv.size() != n) {
        throw std::invalid_argument("har_diagnostics: size mismatch");
    }

    HARDiagnosticsResult result;

    // 1. 残差 Ljung-Box 检验 (lag 自动选择: min(10, n/5))
    const Size lb_lag = std::min(static_cast<Size>(10), n / 5);
    result.residual_ljung_box = econometrics::ljung_box_test(residuals, lb_lag);

    // 2. 残差 Jarque-Bera 检验
    result.residual_jb = econometrics::jarque_bera_test(residuals);

    // 3. Mincer-Zarnowitz 回归 (actual ~ alpha + beta * fitted)
    result.mz_regression = econometrics::mincer_zarnowitz_regression(actual_rv, fitted_rv);

    // 4. R² 和调整 R² (从 MZ 回归中获取, MZ R² = 预测精度)
    result.r_squared = result.mz_regression.r_squared;
    // 调整 R² 需要从残差计算 (MZ 回归是简单回归, k=2)
    // adj_r² = 1 - (1-R²)*(n-1)/(n-k), k=2 (截距 + 斜率)
    result.adjusted_r_squared = (n > 2)
        ? 1.0 - (1.0 - result.r_squared) * static_cast<Real>(n - 1) /
                    static_cast<Real>(n - 2)
        : 0.0;

    // 5. model_adequate: LB 不拒绝 AND JB 不拒绝 AND MZ 不拒绝
    result.model_adequate = !result.residual_ljung_box.base.reject_null &&
                            !result.residual_jb.base.reject_null &&
                            !result.mz_regression.base.reject_null;

    return result;
}

// =============================================================================
// HEAVYDiagnosticsResult - HEAVY 模型诊断结果
//
// HEAVY (High-frEquency-bAsed VolatilitY) 模型 (Shephard-Sheppard 2010):
//   1. 测量方程: RM_t = μ_t + ε_t, Var(ε_t|F_{t-1}) = h_t
//   2. 方差方程: h_t = ω + α·RM_{t-1}² + β·h_{t-1}
//
// 诊断:
//   - 测量方程标准化残差 z_t = ε_t/√h_t 应 iid N(0,1)
//   - 方差方程标准化残差应 iid N(0,1)
//   - h_t 与 RM_t 相关性应高 (接近 1, 波动率聚集效应)
//
// 排幻觉点 H8: z_t² LB 是 ARCH 效应检验关键 (非 z_t LB)
// =============================================================================
struct HEAVYDiagnosticsResult {
    VolatilityDiagnosticsResult variance_equation;     // h_t 方程标准化残差诊断
    VolatilityDiagnosticsResult measurement_equation;  // RM 方程标准化残差诊断

    // HEAVY 特有: 两方程交叉诊断
    Real correlation_h_rm;  // h_t 与 RM_t 的相关性 (应接近 1, 波动率聚集)
    bool model_adequate;
};

// =============================================================================
// heavy_diagnostics - HEAVY 模型诊断
//
// @param rm_residuals RM 方程残差 (ε_t = RM_t - μ_t)
// @param rm_conditional_means RM 条件均值 (μ_t, 用于重建 RM_t = μ_t + ε_t)
// @param variance_residuals 方差方程残差 (h_t - E[h_t|F_{t-1}])
// @param conditional_variances h_t (条件方差)
// @return HEAVYDiagnosticsResult
// =============================================================================
inline HEAVYDiagnosticsResult heavy_diagnostics(
    const std::vector<Real>& rm_residuals,
    const std::vector<Real>& rm_conditional_means,
    const std::vector<Real>& variance_residuals,
    const std::vector<Real>& conditional_variances) {

    const Size n = conditional_variances.size();
    if (n < 10) {
        throw std::invalid_argument("heavy_diagnostics: need at least 10 observations");
    }
    if (rm_residuals.size() != n || rm_conditional_means.size() != n ||
        variance_residuals.size() != n) {
        throw std::invalid_argument("heavy_diagnostics: size mismatch");
    }

    HEAVYDiagnosticsResult result;

    // 1. 测量方程诊断: 标准化残差 = ε_t / √h_t
    //    残差 = rm_residuals, 条件方差 = conditional_variances (h_t)
    result.measurement_equation =
        econometrics::volatility_diagnostics(rm_residuals, conditional_variances);

    // 2. 方差方程诊断: 标准化残差 = variance_residuals / √h_t
    //    (GARCH 型方差方程的扰动项方差与 h_t 成比例)
    result.variance_equation =
        econometrics::volatility_diagnostics(variance_residuals, conditional_variances);

    // 3. h_t 与 RM_t 的相关性 (RM_t = μ_t + ε_t)
    //    波动率聚集: h_t 高时 RM_t 也高, 相关性应接近 1
    Real mean_h = 0.0, mean_rm = 0.0;
    for (Size i = 0; i < n; ++i) {
        const Real rm_t = rm_conditional_means[i] + rm_residuals[i];
        mean_h += conditional_variances[i];
        mean_rm += rm_t;
    }
    mean_h /= static_cast<Real>(n);
    mean_rm /= static_cast<Real>(n);

    Real cov_h_rm = 0.0, var_h = 0.0, var_rm = 0.0;
    for (Size i = 0; i < n; ++i) {
        const Real rm_t = rm_conditional_means[i] + rm_residuals[i];
        const Real dh = conditional_variances[i] - mean_h;
        const Real drm = rm_t - mean_rm;
        cov_h_rm += dh * drm;
        var_h += dh * dh;
        var_rm += drm * drm;
    }

    const Real denom = std::sqrt(var_h * var_rm);
    result.correlation_h_rm = (denom > 1e-300) ? cov_h_rm / denom : 0.0;

    // 4. model_adequate: 双方程均 adequate AND 相关性 > 0.3
    //    (0.3 是宽松下限, 实践中 HEAVY 的 h_t/RM_t 相关性通常 > 0.5)
    result.model_adequate = result.variance_equation.model_adequate &&
                            result.measurement_equation.model_adequate &&
                            result.correlation_h_rm > 0.3;

    return result;
}

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
