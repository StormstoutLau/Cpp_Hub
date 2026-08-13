// =============================================================================
// volatility_diagnostics.hpp - 波动率模型标准化残差检验 (P0)
//
// Phase 7A Wave 1: 通用证伪统计量 (ADR-015 方案 B)
//
// 原理: 若波动率模型正确设定, z_t = ε_t/√h_t 应 ~ iid N(0,1)
//   - z_t 无自相关 → LB(z_t) 不拒绝
//   - z_t² 无自相关 → LB(z_t²) 不拒绝 (ARCH 效应已消除)
//   - z_t 正态 → JB(z_t) 不拒绝
//
// 排幻觉点 H8: z_t² 的 LB 检验是关键 (非仅 z_t)
//   z_t 无自相关但 z_t² 有自相关 → GARCH 未充分捕捉条件异方差
//
// ADR-015 方案 B: 仅依赖 core/, 不依赖 linalg_dynamic.hpp (Eigen3)
// 复用 residual_diagnostics.hpp 的 ljung_box_test / jarque_bera_test
//
// 教材锚点: Tsay 3ed Ch.3, McNeil-Frey-Embrechts 2005 §5.3
// =============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/inference/residual_diagnostics.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// 波动率模型标准化残差综合诊断结果 (复合诊断, 不组合 base — ADR-015 决策点 3)
struct VolatilityDiagnosticsResult {
    std::vector<Real> standardized_residuals;  // z_t = ε_t/√h_t
    LjungBoxResult z_ljung_box;                // z_t 自相关检验
    LjungBoxResult z_squared_ljung_box;        // z_t² 自相关检验 (ARCH 效应, 排幻觉点 H8 关键)
    JarqueBeraResult z_jarque_bera;            // z_t 正态性检验
    Real weighted_lb_statistic;                // 加权 LB 统计量 (Fisher-Gallagher 2012)
    Real weighted_lb_p_value;                  // 加权 LB p 值
    bool model_adequate;                       // 所有检验均不拒绝 → true
};

/// @brief 波动率模型标准化残差综合诊断
///
/// 检验 GARCH/EGARCH/GJR-GARCH/HEAVY 等波动率模型的标准化残差是否满足
/// iid N(0,1) 假设。同时执行 z_t 自相关、z_t² 自相关 (ARCH 效应)、
/// z_t 正态性三类检验。
///
/// @param residuals 原始残差 ε_t (长度 N)
/// @param conditional_variances 条件方差 h_t (长度 N, 必须 > 0)
/// @param lag LB 滞后阶数 (0 = 自动选择 min(10, N/5))
/// @return VolatilityDiagnosticsResult 综合诊断结果
inline VolatilityDiagnosticsResult volatility_diagnostics(
    const std::vector<Real>& residuals,
    const std::vector<Real>& conditional_variances,
    Size lag = 0) {

    const Size n = residuals.size();
    if (n < 5) {
        throw std::invalid_argument("volatility_diagnostics: need at least 5 observations");
    }
    if (conditional_variances.size() != n) {
        throw std::invalid_argument(
            "volatility_diagnostics: residuals and conditional_variances size mismatch");
    }

    // 计算标准化残差 z_t = ε_t / √h_t
    std::vector<Real> z(n);
    for (Size i = 0; i < n; ++i) {
        if (conditional_variances[i] <= 0.0) {
            throw std::runtime_error(
                "volatility_diagnostics: conditional variance must be positive");
        }
        z[i] = residuals[i] / std::sqrt(conditional_variances[i]);
    }

    // z_t 的 LB 检验 (自相关)
    LjungBoxResult z_lb = ljung_box_test(z, lag);

    // z_t² 的 LB 检验 (ARCH 效应, 排幻觉点 H8 关键)
    // 若 z_t² 仍有自相关, 说明 GARCH 未充分捕捉条件异方差
    std::vector<Real> z2(n);
    for (Size i = 0; i < n; ++i) {
        z2[i] = z[i] * z[i];
    }
    LjungBoxResult z2_lb = ljung_box_test(z2, z_lb.lag);  // 复用 z_lb 的 lag 保持一致

    // z_t 的 JB 检验 (正态性)
    JarqueBeraResult z_jb = jarque_bera_test(z);

    // 加权 LB (Fisher-Gallagher 2012)
    // 注: Fisher-Gallagher 2012 的精确加权公式涉及自相关系数的渐近方差修正,
    //     实现较复杂且 spec 未给出精确公式。保守处理: 暂用 z_t² 的标准 Ljung-Box,
    //     避免引入未经验证的加权公式幻觉。待 FG2012 精确公式补充后替换。
    const Real weighted_lb_stat = z2_lb.base.statistic;
    const Real weighted_lb_pval = z2_lb.base.p_value;

    // 模型充分性判断: 所有检验均不拒绝
    const bool model_ok =
        !z_lb.base.reject_null &&
        !z2_lb.base.reject_null &&
        !z_jb.base.reject_null;

    VolatilityDiagnosticsResult result;
    result.standardized_residuals = z;
    result.z_ljung_box = z_lb;
    result.z_squared_ljung_box = z2_lb;
    result.z_jarque_bera = z_jb;
    result.weighted_lb_statistic = weighted_lb_stat;
    result.weighted_lb_p_value = weighted_lb_pval;
    result.model_adequate = model_ok;
    return result;
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
