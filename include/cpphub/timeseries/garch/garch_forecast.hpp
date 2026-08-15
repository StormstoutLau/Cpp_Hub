// =============================================================================
// garch_forecast.hpp - 多步方差预测 (spec §2.0.5)
//
// Phase 7B v1.6 M1 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 教材锚点: Tsay 3ed Ch 5 / arch arch/univariate/base.py:991 (forecast 抽象)
//           + arch/univariate/volatility.py (GARCH.forecast 实现)
//
// 幻觉点防护 G13 (spec §2.0.5):
//   - 指数为 k-1 (非 k), 起点为 h_{T+1} (非 h_T)
//   - 错误公式 h_{T+k} = ω·(1-φ^k)/(1-φ) + φ^k·h_T 混淆起点与递归次数
//   - h_{T+1} = ω + α·ε²_T + β·h_T (已知 ε²_T, 非递归)
//   - h_{T+k} = ω·(1-φ^{k-1})/(1-φ) + φ^{k-1}·h_{T+1}, k ≥ 2
//   - 等价: h_{T+k} = σ̄² + φ^{k-1}·(h_{T+1} - σ̄²), σ̄² = ω/(1-α-β)
//   - k→∞: h_{T+k} → σ̄² (无条件方差)
//   实施: 直接递归 (arch 同款), 闭式公式由测试验证
// =============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/garch/garch_model.hpp"
#include "cpphub/timeseries/garch/egarch_model.hpp"
#include "cpphub/timeseries/garch/gjr_garch_model.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace garch {

/// @brief GARCH(1,1) 多步方差预测 (G13)
///
/// @param params      估计参数 (需 α,β ≥ 0, ω > 0)
/// @param last_variance  最后一期条件方差 h_T
/// @param last_residual  最后一期残差 ε_T
/// @param horizon     预测步数 k (0 → 空)
/// @return {h_{T+1}, ..., h_{T+k}} (直接递归)
inline std::vector<Real> forecast_garch11(const GarchParams& params,
                                          Real last_variance,
                                          Real last_residual,
                                          Size horizon) {
    if (!(last_variance > 0.0) || !std::isfinite(last_variance)) {
        throw std::invalid_argument("forecast_garch11: last_variance must be > 0");
    }
    if (params.omega <= 0.0 || params.alpha < 0.0 || params.beta < 0.0) {
        throw std::invalid_argument("forecast_garch11: invalid params");
    }
    std::vector<Real> fc(horizon);
    // Step 1: h_{T+1} = ω + α·ε²_T + β·h_T (已知 ε²_T)
    Real h_prev = params.omega + params.alpha * last_residual * last_residual
                  + params.beta * last_variance;
    for (Size k = 0; k < horizon; ++k) {
        if (k > 0) {
            // Step 2: E[ε²_{T+k-1}|F_T] = h_{T+k-1} (E[z²]=1)
            h_prev = params.omega
                     + (params.alpha + params.beta) * h_prev;
        }
        fc[k] = h_prev;
    }
    return fc;
}

/// @brief EGARCH(1,1) 多步方差预测 (Nelson 1991, spec §2.0.5)
///
/// 对 log(h) 递归再 exp:
///   Step 1: ln h_{T+1} = ω + β·ln h_T + α·z_T + γ·(|z_T| - E|z|) (已知 z_T)
///   Step 2 (k≥2): E[z_{T+k-1}|F_T]=0, E[|z_{T+k-1}||F_T]=E|z| → 非对称/
///           对称项期望为 0, ln h_{T+k} = ω + β·ln h_{T+k-1}
///   k→∞: ln h_{T+k} → ω/(1-β), h → exp(ω/(1-β)) (中位数式收敛)
///   注: 返回 exp(E[ln h]) (确定性递归), 非 E[h] (simulation 路径均值,
///   Jensen 差为正); arch analytic 仅支持 horizon=1, 多步用 simulation 对照
///
/// @param params            估计参数 (|β| < 1)
/// @param last_log_variance  ln h_T (最后期对数条件方差)
/// @param last_z            标准化残差 z_T = ε_T/√h_T
/// @param horizon           预测步数 k (0 → 空)
inline std::vector<Real> forecast_egarch(const EGarchParams& params,
                                         Real last_log_variance,
                                         Real last_z,
                                         Size horizon) {
    if (!std::isfinite(last_log_variance)) {
        throw std::invalid_argument("forecast_egarch: last_log_variance NaN");
    }
    if (!std::isfinite(last_z)) {
        throw std::invalid_argument("forecast_egarch: last_z NaN");
    }
    if (!(std::abs(params.beta) < 1.0)) {
        throw std::invalid_argument("forecast_egarch: |beta| >= 1");
    }
    std::vector<Real> fc(horizon);
    // Step 1: ln h_{T+1} 含已知 z_T (G6: 非对称/对称项用 z_T)
    Real ln_h_prev = params.omega + params.beta * last_log_variance
                     + params.alpha * last_z
                     + params.gamma * (std::abs(last_z) - EGARCH_E_ABS_Z);
    for (Size k = 0; k < horizon; ++k) {
        if (k > 0) {
            // Step 2: 新息项期望 0 (E[z]=0, E[|z|-E|z|]=0)
            ln_h_prev = params.omega + params.beta * ln_h_prev;
        }
        fc[k] = std::exp(ln_h_prev);
        if (!std::isfinite(fc[k])) {
            throw std::runtime_error("forecast_egarch: variance overflow");
        }
    }
    return fc;
}

/// @brief GJR-GARCH(1,1) 多步方差预测 (spec §2.0.5)
///
///   Step 1: h_{T+1} = ω + α·ε²_T + γ·I(ε_T<0)·ε²_T + β·h_T (已知 ε²_T)
///   Step 2 (k≥2): E[I(z<0)·ε²|F_T] = h/2 (对称分布, probe-2 实测) →
///           h_{T+k} = ω + φ·h_{T+k-1}, φ = α + γ/2 + β (G10 平稳性系数)
///   闭式: h_{T+k} = σ̄² + φ^{k-1}·(h_{T+1} - σ̄²), σ̄² = ω/(1-α-γ/2-β)
///   k→∞: h_{T+k} → σ̄²
inline std::vector<Real> forecast_gjr(const GjrGarchParams& params,
                                      Real last_variance,
                                      Real last_residual,
                                      Size horizon) {
    if (!(last_variance > 0.0) || !std::isfinite(last_variance)) {
        throw std::invalid_argument("forecast_gjr: last_variance must be > 0");
    }
    if (params.omega <= 0.0 || params.alpha < 0.0 || params.beta < 0.0) {
        throw std::invalid_argument("forecast_gjr: invalid params");
    }
    const Real phi = params.alpha + params.gamma / 2.0 + params.beta;
    if (!(phi < 1.0)) {
        throw std::invalid_argument("forecast_gjr: alpha+gamma/2+beta >= 1");
    }
    std::vector<Real> fc(horizon);
    // Step 1: h_{T+1} 含已知 ε²_T 与指示项 I(ε_T<0) (G7/G8)
    Real h_prev = params.omega
                  + params.alpha * last_residual * last_residual
                  + params.gamma * (last_residual < 0.0 ? 1.0 : 0.0)
                        * last_residual * last_residual
                  + params.beta * last_variance;
    for (Size k = 0; k < horizon; ++k) {
        if (k > 0) {
            // Step 2: 非对称项期望折半 (E[I(z<0)]=1/2)
            h_prev = params.omega + phi * h_prev;
        }
        fc[k] = h_prev;
    }
    return fc;
}

}  // namespace garch
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
