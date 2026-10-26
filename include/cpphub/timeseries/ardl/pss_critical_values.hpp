// =============================================================================
// pss_critical_values.hpp - PSS 边限检验响应面临界值与 p 值
//
// Phase 8 v1.8 M1 (PHASE8_SPEC.md v1.1 §2.4)
//
// 理论锚点: Pesaran-Shin-Smith 2001 (J Appl Econ 16(3):289-326) Tables CI(i)-(v)
//   临界值来源: statsmodels 0.14.6 pss_critical_values.py 响应面系数机器转录
//   (每键 3200 万次模拟 vs PSS 原文 4 万次; 生成脚本 gen_phase8_pss_inc.py)
//
// 求值语义 (probe_pss_cv.py [B] 冻结, 与 statsmodels _pss_pvalue 数值自证 ≤1e-16):
//   临界值: PSS_CRIT_VALS[idx] 表常量直查 (非多项式求值)
//     四档 = {90%, 95%, 99%, 99.9%} (statsmodels crit_percentiles 口径, PSS9:
//     PSS 原表为 {10,5,2.5,1}%, 无 0.1% 档而有 2.5% 档 — 逐列对照 PSS 表会错读)
//   p 值: p = 1 − Φ(c₀ + c₁·x + c₂·x² + c₃·x³), x = log(stat)
//     系数切换: stat > stat_star ? small_p (2 阶) : large_p (3 阶)
//     (阈值比较用原始 stat 尺度, 非对数尺度 — statsmodels model.py L2453 实录)
//
// 排幻觉点:
//   PSS5  分位档错位: 本表四档 {10,5,1,0.1}%, 第 3 档是 1% (PSS 是 2.5%)
//   PSS6  k = force 变量数 (不含 y 自身滞后), 与 ARDL(p,q) 的 q 无关
//   PSS10 statsmodels bounds_test 调用层用 nvar = 1+k 偏移查表 (表键本身是 PSS
//         语义); 本接口 k 即 force 变量数, 不镜像该偏移 — danish 回归断言防误镜像
//   A5    仅渐近响应面; MC 模拟路径 (asymptotic=False) 与 Narayan 分档表均不做
//
// 依赖: 纯标量 (<cmath>/<array>/<stdexcept> + core/types.hpp), 无 Eigen —
//   M1 模块整体挂 cpphub_timeseries_mat (K5), 但本头文件零矩阵依赖
// =============================================================================
#pragma once

#include <array>
#include <cmath>
#include <stdexcept>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace ardl {

namespace detail {

// 系数表 (mackinnon_tables.inc 同款惯例: 于 namespace detail 内 include)
#include "cpphub/timeseries/ardl/pss_response_surface.inc"

/// 键 (k, case, I1) → 扁平索引 (constexpr, 供 static_assert 编译期求值)
/// k ∈ 1..10, case ∈ 1..5; 序 = k 外层, case 中层, I1 内层
constexpr Size pss_flat_index(Size k, Size pss_case, bool upper_i1) {
    return (k - 1) * 10 + (pss_case - 1) * 2 + (upper_i1 ? 1 : 0);
}

/// 响应面 p 值: Φ(−Σ cᵢ·log(stat)ⁱ) 上尾直算, 系数按 stat 与 stat_star 比较 (原始尺度)
inline Real pss_response_pvalue(Real stat, Size k, Size pss_case, bool upper_i1) {
    const Size i = pss_flat_index(k, pss_case, upper_i1);
    const Real x = std::log(stat);
    Real poly = 0.0;
    if (stat > PSS_STAT_STAR[i]) {
        // small_p: 2 阶 (3 系数), 升幂
        poly = PSS_SMALL_P[i][0] + x * (PSS_SMALL_P[i][1] + x * PSS_SMALL_P[i][2]);
    } else {
        // large_p: 3 阶 (4 系数), 升幂
        poly = PSS_LARGE_P[i][0] +
               x * (PSS_LARGE_P[i][1] + x * (PSS_LARGE_P[i][2] + x * PSS_LARGE_P[i][3]));
    }
    // 深尾上尾直算 (1−Φ 族第 4 例, RISK-036 同源): 旧式 1−Φ(poly) 在 poly
    // ≥8.25 灾难性抵消精确归零 (锚键 stat=100 真值 8.03e-26 实测报 0; 全
    // 100 键双支路均可达, RESEARCH_PSS_DEEPTAIL.md §2)。Φ(−x) = 0.5·erfc(x/√2)
    // 恒等式直算上尾: 中尾与旧式逐位一致或更准 (stat≤2 diff=0), 深尾保持
    // ≤1e-13 相对精度; stat≳1e5 真值 <1e-308 纯欠流归零, 非缺陷。
    return 0.5 * std::erfc(poly / std::sqrt(2.0));
}

/// 键域校验 (k ∈ 1..10, case ∈ 1..5; probe [A]: 表无 k=0 键)
inline void pss_check_key(Size k, Size pss_case) {
    if (k < 1 || k > 10) {
        throw std::invalid_argument(
            "PSS critical values: k (number of force variables) must be in [1, 10] "
            "(response surface coverage; statsmodels falls back to simulation beyond)");
    }
    if (pss_case < 1 || pss_case > 5) {
        throw std::invalid_argument("PSS critical values: case must be in [1, 5]");
    }
}

}  // namespace detail

// ---------------------------------------------------------------------------
// 编译期锚 (probe_pss_cv.py 冻结值; 四张表各至少一锚)
// ---------------------------------------------------------------------------
// 双锚 (spec §2.4): Case III k=1 5% 两界 — PSS 2001 Table CI(iii) 同位值
// 4.94/5.73 的 statsmodels 3200 万次模拟精化版 (Statalist Kripfganz 帖子旁证)
static_assert(detail::PSS_CRIT_VALS[detail::pss_flat_index(1, 3, false)][1] == 4.9228967,
              "PSS crit_vals anchor (1,3,I0) 5% = 4.9228967");
static_assert(detail::PSS_CRIT_VALS[detail::pss_flat_index(1, 3, true)][1] == 5.7217442,
              "PSS crit_vals anchor (1,3,I1) 5% = 5.7217442");
// 辅锚: stat_star / large_p / small_p 各一 (同键位, probe [B] 消费路径)
static_assert(detail::PSS_STAT_STAR[detail::pss_flat_index(1, 3, false)] ==
                  0.846390593107207,
              "PSS stat_star anchor (1,3,I0)");
static_assert(detail::PSS_LARGE_P[detail::pss_flat_index(1, 3, false)][0] == -0.62883,
              "PSS large_p anchor (1,3,I0) c0");
static_assert(detail::PSS_SMALL_P[detail::pss_flat_index(1, 3, false)][0] == -0.57887,
              "PSS small_p anchor (1,3,I0) c0");

// ---------------------------------------------------------------------------
// 公共接口
// ---------------------------------------------------------------------------

/// PSS 边限检验临界值 (渐近, 响应面表直查)
///
/// \param k         force 变量数 (PSS 语义, 不含 y 自身滞后; PSS6/PSS10)
/// \param pss_case  PSS Case 1..5 (编号跨 PSS/statsmodels/R ARDL/pssbounds 一致)
/// \param upper_i1  true = I(1) 上界 (x 含 I(1) 分量), false = I(0) 下界 (x 平稳)
/// \return 四档临界值 {10%, 5%, 1%, 0.1%} (statsmodels 口径, PSS9)
/// \throw std::invalid_argument k ∉ [1,10] 或 case ∉ [1,5]
std::array<Real, 4> pss_crit_values(Size k, Size pss_case, bool upper_i1) {
    detail::pss_check_key(k, pss_case);
    const auto& row = detail::PSS_CRIT_VALS[detail::pss_flat_index(k, pss_case, upper_i1)];
    return {row[0], row[1], row[2], row[3]};
}

/// PSS 边限检验 F 统计量 p 值界 (双 null 近似, 镜像 statsmodels p_values 语义)
///
/// 返回 {p_lower, p_upper}:
///   [0] = H0: 全 I(0) (下界 null) 下的 p 值 — 键 (k, case, I1=false)
///   [1] = H0: 全 I(1) (上界 null) 下的 p 值 — 键 (k, case, I1=true)
/// 两 p 值同时小于 α ⇒ 稳健拒绝; 同时大于 α ⇒ 稳健不拒绝; 介于 ⇒ 不确定区
/// (与三值结论 BoundsConclusion 的对偶判据, 由 uecm_bounds.hpp 消费)
///
/// 仅 F 统计量有响应面覆盖 (probe [D]: statsmodels 无 t p 值路径);
/// t 检验的 p 值不在本接口 (BoundsTestResult.has_p=false, A5 不做 MC)
///
/// \param stat  F 统计量 (Wald 二次型/约束数, PSS7; 必须 > 0)
/// \throw std::invalid_argument k/case 越界或 stat ≤ 0 / NaN
std::array<Real, 2> pss_pvalue_bounds(Real stat, Size k, Size pss_case) {
    detail::pss_check_key(k, pss_case);
    if (std::isnan(stat) || stat <= 0.0) {
        throw std::invalid_argument(
            "pss_pvalue_bounds: stat must be positive and non-NaN (log domain)");
    }
    return {detail::pss_response_pvalue(stat, k, pss_case, false),
            detail::pss_response_pvalue(stat, k, pss_case, true)};
}

}  // namespace ardl
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
