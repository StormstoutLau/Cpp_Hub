// =============================================================================
// ericsson_mackinnon_cv.hpp - Ericsson-MacKinnon 2002 ECT t 检验临界值查表
//
// Phase 7C v1.7 M3 (PHASE7C_SPEC.md §5.4; 决策 21, CI10)
//
// 用途: VECM 逐方程 ECT 系数 t 检验的非标准临界值 (rank=1 诊断场景)
//   — t 统计量本身用 ML 协方差 (与 statsmodels VECM summary 的 alpha t 一致);
//     但 H0 "该方程无误差修正" 下 ECT 项可能非平稳, t 分布非标准 → EM2002 查表
//
// 溯源: Fed IFDP 655 (1999 工作论文版, EM2002 同源) Table 2-5
//   (转录: em2002_ect_cv.inc; 三重验证: n=1 行 vs MacKinnon 2010 ≤0.0005,
//    表 7 K_ctt(3) T=51 有限样本 CV 精确复现, 144 行结构完整;
//    ⚠️ 正式发表版 (2002) 或与 WP 版微差 — 发表版付费墙, WP 版为可核验官方源)
//
// 公式: CV(T) = θ∞ + θ1/T + θ2/T² + θ3/T³  (T = ECM 回归有效样本量)
// =============================================================================

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace cointegration {

// 系数表 (自动生成, 勿手改; 布局 [case][n-1][水平][系数] 见 .inc 头注)
#include "cpphub/timeseries/cointegration/em2002_ect_cv.inc"

namespace detail {

inline Size em2002_case_index(const std::string& det) {
    if (det == "n") return 0;
    if (det == "c") return 1;
    if (det == "ct") return 2;
    if (det == "ctt") return 3;
    throw std::invalid_argument(
        "em2002_ect_cv: det must be n/c/ct/ctt (deterministic case)");
}

}  // namespace detail

// ---------------------------------------------------------------------------
// EM2002 ECT t 检验临界值
//
// @param n_vars ECM 系统变量总数 n (含 LHS; r=1 VECM 的 K) ∈ [1, 12]
// @param det 确定性情形 "n"(无确定项)/"c"(常数)/"ct"(常数+趋势)/"ctt"(+二次)
// @param level_pct 显著性水平百分数: 1.0 / 5.0 / 10.0
// @param T_eff ECM 回归有效样本量 (>0)
// @return 临界值 (左尾, 拒绝域 t < CV)
// ---------------------------------------------------------------------------
inline Real em2002_ect_critical_value(Size n_vars, const std::string& det,
                                      Real level_pct, Size T_eff) {
    if (n_vars < 1 || n_vars > 12) {
        throw std::invalid_argument(
            "em2002_ect_critical_value: n must be in [1, 12]");
    }
    Size lvl;
    if (level_pct == 1.0) {
        lvl = 0;
    } else if (level_pct == 5.0) {
        lvl = 1;
    } else if (level_pct == 10.0) {
        lvl = 2;
    } else {
        throw std::invalid_argument(
            "em2002_ect_critical_value: level_pct must be 1.0/5.0/10.0");
    }
    if (T_eff == 0) {
        throw std::invalid_argument(
            "em2002_ect_critical_value: T_eff must be > 0");
    }
    const Size c = detail::em2002_case_index(det);
    const Real* th = &EM2002_THETA[c][n_vars - 1][lvl][0];
    const Real inv_t = 1.0 / static_cast<Real>(T_eff);
    return th[0] + th[1] * inv_t + th[2] * inv_t * inv_t +
           th[3] * inv_t * inv_t * inv_t;
}

// 渐近临界值 (T → ∞: 纯 θ∞)
inline Real em2002_ect_critical_value_asymptotic(Size n_vars,
                                                 const std::string& det,
                                                 Real level_pct) {
    if (n_vars < 1 || n_vars > 12) {
        throw std::invalid_argument(
            "em2002_ect_critical_value_asymptotic: n must be in [1, 12]");
    }
    Size lvl;
    if (level_pct == 1.0) {
        lvl = 0;
    } else if (level_pct == 5.0) {
        lvl = 1;
    } else if (level_pct == 10.0) {
        lvl = 2;
    } else {
        throw std::invalid_argument(
            "em2002_ect_critical_value_asymptotic: level_pct must be 1/5/10");
    }
    const Size c = detail::em2002_case_index(det);
    return EM2002_THETA[c][n_vars - 1][lvl][0];
}

}  // namespace cointegration
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
