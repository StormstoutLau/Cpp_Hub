// =============================================================================
// osterwald_lenum_cv.hpp - Johansen 检验临界值查表: OL1992 + MHM96 双表
//
// Phase 7C v1.7 M3 (PHASE7C_SPEC.md §5.2; 决策 19)
//
// 双表并存 (JOHANSEN_DUAL_LIB_DIFF.md 冻结决策, 2026-08-19):
//   - OL1992 (urca 转录): ecdet 键控 {none, const, trend}, q=1..11, 10/5/1%
//     情形语义 (urca): none = 短回归无约束常数 (H1*); const = 常数限入协整
//     关系 (H1); trend = 趋势限入协整关系 (H*) + 无约束常数
//   - MHM96 (statsmodels 转录): det_order 键控 {-1, 0, 1}, n=1..12, 90/95/99%
//     JohansenResult.cvt/cvm 默认表源 (与 statsmodels coint_johansen 逐位一致,
//     select_coint_rank 复用同表 — B4)
//
// ⚠️ 双库情形映射警告 (CI5 + diff 报告 §2/§5):
//   - 仅 det_order=0 ↔ ecdet="none" 可映射 (计算已证 1e-10 一致)
//   - det_order=-1 (无确定项) 与 det_order=1 (预去势) 无 urca 对应;
//     urca const/trend 无 statsmodels 对应
//   - MC 裁决: statsmodels det_order∈{0,1} 的 MHM96 表与其自身统计量分布
//     不一致 (q=1 行呈 χ²(1) 型 — 受限情形特征; 其计算遵循无约束情形分布,
//     与 OL1992 none/trend 表吻合)。实际协整 rank 推断建议:
//     det_order=0 → ol1992_trace_cv("none", q); det_order=1 → 近似用
//     ol1992_*_cv("trend", q); det_order=-1 → MHM96 表自洽可用
//
// 排幻觉点 (spec §9.3): CI5 (statsmodels 一套表 + 真正的两套 = SM(MHM96)
//   vs urca(OL1992)); 表索引 q = N − r (非 r 本身)
// =============================================================================

#pragma once

#include <array>
#include <stdexcept>
#include <string>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace cointegration {

// 系数表 (自动生成, 勿手改; 布局与 static_assert 锚见各 .inc 头注)
#include "cpphub/timeseries/cointegration/ol1992_cv.inc"
#include "cpphub/timeseries/cointegration/mhm96_johansen_cv.inc"

namespace detail {

inline Size ecdet_index(const std::string& ecdet) {
    if (ecdet == "none") return 0;
    if (ecdet == "const") return 1;
    if (ecdet == "trend") return 2;
    throw std::invalid_argument(
        "ol1992: ecdet must be none/const/trend (urca semantics)");
}

inline Size det_order_index(int det_order) {
    if (det_order == -1) return 0;
    if (det_order == 0) return 1;
    if (det_order == 1) return 2;
    throw std::invalid_argument("mhm96: det_order must be -1/0/1");
}

}  // namespace detail

// ---------------------------------------------------------------------------
// OL1992 查表 (urca 语义, 10/5/1% 列序与 urca cval 一致)
// @param ecdet "none"/"const"/"trend" (urca ca.jo 语义, 直译不转译 §1.4-4)
// @param q 检验 r 时的自由方向数 = N − r ∈ [1, 11] (urca 仅 N<11 报告)
// ---------------------------------------------------------------------------
inline std::array<Real, 3> ol1992_trace_cv(const std::string& ecdet, Size q) {
    if (q < 1 || q > 11) {
        throw std::invalid_argument("ol1992_trace_cv: q = N-r must be in [1, 11]");
    }
    const Size e = detail::ecdet_index(ecdet);
    return {OL1992_TRACE[e][q - 1][0], OL1992_TRACE[e][q - 1][1],
            OL1992_TRACE[e][q - 1][2]};
}

inline std::array<Real, 3> ol1992_maxeig_cv(const std::string& ecdet, Size q) {
    if (q < 1 || q > 11) {
        throw std::invalid_argument("ol1992_maxeig_cv: q = N-r must be in [1, 11]");
    }
    const Size e = detail::ecdet_index(ecdet);
    return {OL1992_MAXEIG[e][q - 1][0], OL1992_MAXEIG[e][q - 1][1],
            OL1992_MAXEIG[e][q - 1][2]};
}

// ---------------------------------------------------------------------------
// MHM96 查表 (statsmodels 语义, 90/95/99% 列序与 c_sjt/c_sja 一致)
// @param det_order -1/0/1 (statsmodels coint_johansen 语义)
// @param n = N − r ∈ [1, 12]
// ---------------------------------------------------------------------------
inline std::array<Real, 3> mhm96_trace_cv(int det_order, Size n) {
    if (n < 1 || n > 12) {
        throw std::invalid_argument("mhm96_trace_cv: n = N-r must be in [1, 12]");
    }
    const Size d = detail::det_order_index(det_order);
    return {MHM96_TRACE[d][n - 1][0], MHM96_TRACE[d][n - 1][1],
            MHM96_TRACE[d][n - 1][2]};
}

inline std::array<Real, 3> mhm96_maxeig_cv(int det_order, Size n) {
    if (n < 1 || n > 12) {
        throw std::invalid_argument("mhm96_maxeig_cv: n = N-r must be in [1, 12]");
    }
    const Size d = detail::det_order_index(det_order);
    return {MHM96_MAXEIG[d][n - 1][0], MHM96_MAXEIG[d][n - 1][1],
            MHM96_MAXEIG[d][n - 1][2]};
}

}  // namespace cointegration
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
