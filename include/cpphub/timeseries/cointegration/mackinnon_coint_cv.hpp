// =============================================================================
// mackinnon_coint_cv.hpp - MacKinnon 协整检验 p 值 (1994) 与临界值 (2010 响应面)
//
// Phase 7C v1.7 M3 (PHASE7C_SPEC.md §5.1; 决策 18)
//
// 消费方: engle_granger.hpp (EG 两步法的 Step 3/4)
//
// ⚠️ p 与 cv 不同源 (CI2, statsmodels issue #4138 官方确认):
//   - p 值: MacKinnon 1994 渐近近似 (mackinnonp, N=1..6)
//   - cv:   MacKinnon 2010 响应面 (mackinnoncrit, N=1..12) 含 nobs 小样本修正
//   小样本下可出现 "过 1% cv 而 p > 1%" — API 文档显式声明, 测试分列断言
//
// ⚠️ spec §5.1 原文 "1994 协整响应面 5 系数 4 次" 与实测不符 (diff 报告 §7
//   偏离记录): statsmodels coint() 的 cv 用 2010 表 4 系数 3 次;
//   1994 机制仅用于 p 值。本实现按实测转录。
//
// 公式 (与 7B mackinnon_cv.hpp 同约定, 系数升幂):
//   CV(T) = β∞ + β1/T + β2/T² + β3/T³     (T = nobs − 1, Stata egranger 约定)
//   p = Φ(d0 + d1·t + d2·t² [+ d3·t³])    (分段: ≤τ* smallp 二次 / >τ* largep 三次)
//   t > τ_max → p = 1;  t < τ_min → p = 0
// =============================================================================

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"  // detail::normal_cdf

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace cointegration {

// 系数表 (自动生成, 勿手改)
#include "cpphub/timeseries/cointegration/mackinnon1994_coint.inc"
#include "cpphub/timeseries/cointegration/mackinnon2010_coint_cv.inc"

namespace detail {

inline Size trend_index_1994(const std::string& trend) {
    if (trend == "n") return 0;
    if (trend == "c") return 1;
    if (trend == "ct") return 2;
    if (trend == "ctt") return 3;
    throw std::invalid_argument(
        "mackinnon_coint: trend must be n/c/ct/ctt");
}

// 多项式求值 (升幂): y = c0 + c1·x + c2·x² + ...
inline Real polyval_asc(const Real* coef, Size n, Real x) {
    Real y = 0.0;
    Real xk = 1.0;
    for (Size i = 0; i < n; ++i) {
        y += coef[i] * xk;
        xk *= x;
    }
    return y;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// MacKinnon 1994 协整检验 p 值 (渐近近似, 复刻 statsmodels mackinnonp)
//
// @param statistic EG 第二步的 t 统计量
// @param trend 第一步回归的确定性形式 "n"/"c"/"ct"/"ctt"
// @param n_vars 协整系统变量总数 N (EG 双变量 = 2; 表覆盖 N=1..6)
// ---------------------------------------------------------------------------
inline Real mackinnon_coint_p_value(Real statistic, const std::string& trend,
                                    Size n_vars) {
    if (n_vars < 1 || n_vars > 6) {
        throw std::invalid_argument(
            "mackinnon_coint_p_value: N must be in [1, 6] (1994 table range)");
    }
    const Size t = detail::trend_index_1994(trend);
    const Size i = n_vars - 1;
    const Real tmax = MK1994_TAUMAX[t][i];
    const Real big = 1e299;  // inf 编码 (N=1 n 情形)
    if (statistic > tmax && tmax >= big) return 1.0;
    if (statistic > tmax) return 1.0;
    if (statistic < MK1994_TAUMIN[t][i]) return 0.0;
    if (statistic <= MK1994_TAUSTAR[t][i]) {
        return unit_root::detail::normal_cdf(
            detail::polyval_asc(MK1994_SMALLP[t][i], 3, statistic));
    }
    return unit_root::detail::normal_cdf(
        detail::polyval_asc(MK1994_LARGEP[t][i], 4, statistic));
}

// ---------------------------------------------------------------------------
// MacKinnon 2010 协整检验临界值 (响应面, 复刻 statsmodels mackinnoncrit)
//
// @param n_vars 系统变量总数 N ∈ [1, 12] (EG 双变量 = 2)
// @param trend "c"/"ct"/"ctt"; "n" 无 2010 表 → {NaN, NaN, NaN}
//   (statsmodels coint(trend="n") 同样返回 NaN, CI2 注记)
// @param T_eff 有效样本 (statsmodels coint 传 nobs − 1; 0 → 渐近 β∞)
// @return {1%, 5%, 10%} 临界值
// ---------------------------------------------------------------------------
inline std::array<Real, 3> mackinnon_coint_critical_values(
    Size n_vars, const std::string& trend, Size T_eff) {
    if (n_vars < 1 || n_vars > 12) {
        throw std::invalid_argument(
            "mackinnon_coint_critical_values: N must be in [1, 12]");
    }
    if (trend == "n") {
        return {std::numeric_limits<Real>::quiet_NaN(),
                std::numeric_limits<Real>::quiet_NaN(),
                std::numeric_limits<Real>::quiet_NaN()};
    }
    if (T_eff == 0) {
        throw std::invalid_argument(
            "mackinnon_coint_critical_values: T_eff=0 not supported "
            "(pass asymptotic values via table)");
    }
    const Size i = n_vars - 1;
    const Real inv_t = 1.0 / static_cast<Real>(T_eff);
    const Real* tbl = nullptr;
    if (trend == "c") {
        tbl = &T2010_C_CV[i][0][0];
    } else if (trend == "ct") {
        tbl = &T2010_CT_CV[i][0][0];
    } else if (trend == "ctt") {
        tbl = &T2010_CTT_CV[i][0][0];
    } else {
        throw std::invalid_argument(
            "mackinnon_coint_critical_values: trend must be n/c/ct/ctt");
    }
    std::array<Real, 3> out{};
    for (Size lvl = 0; lvl < 3; ++lvl) {
        out[lvl] = detail::polyval_asc(&tbl[lvl * 4], 4, inv_t);
    }
    return out;
}

// 渐近临界值 (T → ∞: 纯 β∞ 列)
inline std::array<Real, 3> mackinnon_coint_critical_values_asymptotic(
    Size n_vars, const std::string& trend) {
    if (n_vars < 1 || n_vars > 12) {
        throw std::invalid_argument(
            "mackinnon_coint_critical_values_asymptotic: N in [1, 12]");
    }
    if (trend == "n") {
        return {std::numeric_limits<Real>::quiet_NaN(),
                std::numeric_limits<Real>::quiet_NaN(),
                std::numeric_limits<Real>::quiet_NaN()};
    }
    const Size i = n_vars - 1;
    if (trend == "c") {
        return {T2010_C_CV[i][0][0], T2010_C_CV[i][1][0], T2010_C_CV[i][2][0]};
    }
    if (trend == "ct") {
        return {T2010_CT_CV[i][0][0], T2010_CT_CV[i][1][0], T2010_CT_CV[i][2][0]};
    }
    if (trend == "ctt") {
        return {T2010_CTT_CV[i][0][0], T2010_CTT_CV[i][1][0],
                T2010_CTT_CV[i][2][0]};
    }
    throw std::invalid_argument(
        "mackinnon_coint_critical_values_asymptotic: trend must be n/c/ct/ctt");
}

}  // namespace cointegration
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
