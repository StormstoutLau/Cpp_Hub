// =============================================================================
// mackinnon_cv.hpp - MacKinnon 2010 response surface 临界值与 p 值 (spec §3.0.2)
//
// Phase 7B v1.6 M2 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 教材锚点: MacKinnon 2010 JBES / arch arch/unitroot/critical_values/
//
// 幻觉点防护 (spec §6.2):
//   U3: response surface 4 系数 3 次多项式 CV = c3/T³ + c2/T² + c1/T + c0
//       (非 5 系数 4 次 — 那是协整 N≥2 用的, ADF N=1 不用)
//   U13: ADF/PP p 值用 MacKinnon 2010 (非 1996 旧系数); PP 共享 ADF 表
//   U7: DF-GLS 临界值来自 arch 独立模拟 (非 ERS 1996 原表)
//
// 公式约定 (probe 实测 + 基准数值双重验证):
//   临界值: CV = c3/T³ + c2/T² + c1/T + c0   (系数升幂, T=inf → c0)
//   p 值:   z = Σ cᵢ·statⁱ (升幂!), p = Φ(z)
//     stat > tau_star → largep 多项式 (4 系数: 三次)
//     stat ≤ tau_star → smallp 多项式 (3 系数: 二次)
//     stat > tau_max → p=1; stat < tau_min → p=0
//   注: 升幂已用 MKP_ADF_C_T233_M1 手算验证 (降幂给 1-p, 方向错误)
// =============================================================================
#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"  // detail::normal_cdf 共享

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace unit_root {

// 系数表 (自动生成, 勿手改)
#include "cpphub/timeseries/unit_root/mackinnon_tables.inc"

namespace detail {

// 多项式求值 (升幂): y = c0 + c1·x + c2·x² + ...
inline Real polyval_ascending(const Real* coef, Size n_coef, Real x) {
    Real y = 0.0;
    Real xk = 1.0;
    for (Size i = 0; i < n_coef; ++i) {
        y += coef[i] * xk;
        xk *= x;
    }
    return y;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// MacKinnon 2010 response surface 临界值 (U-ADR3, U3/U13)
//
// @param test_type "adf" / "pp" (共享 DF 表) / "df_gls"
// @param trend_spec "nc"/"n", "c", "ct"
// @param T 样本量 (T→∞ 用 0 表示, 直接取 β_∞)
// @param n_params 仅支持 1 (N>1 为协整检验, 不在 v1.6 范围)
// @param p 显著性水平 (0.01 / 0.05 / 0.10)
// ---------------------------------------------------------------------------
inline Real mackinnon_critical_value(const std::string& test_type,
                                     const std::string& trend_spec,
                                     Size T, Size n_params, Real p) {
    if (n_params != 1) {
        throw std::invalid_argument(
            "mackinnon_critical_value: only n_params=1 supported "
            "(cointegration N>=2 not in v1.6 scope)");
    }

    // 表选择: adf/pp 共享, df_gls 独立
    bool is_dfgls;
    if (test_type == "adf" || test_type == "pp") {
        is_dfgls = false;
    } else if (test_type == "df_gls") {
        is_dfgls = true;
    } else {
        throw std::invalid_argument(
            "mackinnon_critical_value: test_type must be adf/pp/df_gls");
    }

    // 行索引: 1%/5%/10% (每 trend 3 行)
    Size row;
    if (p == 0.01) {
        row = 0;
    } else if (p == 0.05) {
        row = 1;
    } else if (p == 0.10) {
        row = 2;
    } else {
        throw std::invalid_argument(
            "mackinnon_critical_value: p must be 0.01, 0.05 or 0.10");
    }

    const Real* coef = nullptr;
    if (is_dfgls) {
        if (trend_spec == "c") {
            coef = DFGLS_CV_C[row];
        } else if (trend_spec == "ct") {
            coef = DFGLS_CV_CT[row];
        } else {
            throw std::invalid_argument(
                "mackinnon_critical_value: df_gls supports only c/ct");
        }
    } else {
        if (trend_spec == "nc" || trend_spec == "n") {
            coef = TAU2010_N[row];
        } else if (trend_spec == "c") {
            coef = TAU2010_C[row];
        } else if (trend_spec == "ct") {
            coef = TAU2010_CT[row];
        } else {
            throw std::invalid_argument(
                "mackinnon_critical_value: trend_spec must be nc/c/ct");
        }
    }

    if (T == 0) {
        return coef[0];  // T→∞: 纯 β_∞
    }
    const Real x = 1.0 / static_cast<Real>(T);
    // CV = c0 + c1/T + c2/T² + c3/T³ (升幂)
    return detail::polyval_ascending(coef, 4, x);
}

// ---------------------------------------------------------------------------
// MacKinnon 2010 p 值 (response surface)
// 分段: stat < tau_min → 0; stat > tau_max → 1;
//       stat > tau_star → largep 三次 (t 大 → p 大); 否则 smallp 二次
// 基准: MKP_ADF_C_T233_{M3p5,M2p5,M1} / MKP_DFGLS_C_{M3p5,M2} (1e-12)
// 排幻觉: t=-1 (> star) 手算 largep 升幂 z=0.6848 → Φ=0.7533 ✓;
//   若误走 smallp 得 z=0.7630 → Φ=0.7773 (实测失败值, 分支方向验证器) 对照)
// ---------------------------------------------------------------------------
inline Real mackinnon_p_value(Real statistic, const std::string& test_type,
                              const std::string& trend_spec,
                              Size T, Size n_params) {
    (void)T;  // arch mackinnonp 不依赖 nobs (表已按渐近校准)
    if (n_params != 1) {
        throw std::invalid_argument(
            "mackinnon_p_value: only n_params=1 supported");
    }

    bool is_dfgls;
    if (test_type == "adf" || test_type == "pp") {
        is_dfgls = false;
    } else if (test_type == "df_gls") {
        is_dfgls = true;
    } else {
        throw std::invalid_argument(
            "mackinnon_p_value: test_type must be adf/pp/df_gls");
    }

    const Real *tau_min, *tau_max, *tau_star, *smallp, *largep;
    Size n_small, n_large;
    if (is_dfgls) {
        if (trend_spec == "c") {
            tau_min = DFGLS_MIN_C;      n_small = 3;
            tau_max = DFGLS_MAX_C;
            tau_star = DFGLS_STAR_C;
            smallp = DFGLS_SMALLP_C;    n_large = 4;
            largep = DFGLS_LARGEP_C;
        } else if (trend_spec == "ct") {
            tau_min = DFGLS_MIN_CT;     n_small = 3;
            tau_max = DFGLS_MAX_CT;
            tau_star = DFGLS_STAR_CT;
            smallp = DFGLS_SMALLP_CT;   n_large = 4;
            largep = DFGLS_LARGEP_CT;
        } else {
            throw std::invalid_argument(
                "mackinnon_p_value: df_gls supports only c/ct");
        }
    } else {
        if (trend_spec == "nc" || trend_spec == "n") {
            tau_min = TAUMIN_N;         n_small = 3;
            tau_max = TAUMAX_N;
            tau_star = TAUSTAR_N;
            smallp = TAUSMALLP_N;       n_large = 4;
            largep = TAULARGEP_N;
        } else if (trend_spec == "c") {
            tau_min = TAUMIN_C;         n_small = 3;
            tau_max = TAUMAX_C;
            tau_star = TAUSTAR_C;
            smallp = TAUSMALLP_C;       n_large = 4;
            largep = TAULARGEP_C;
        } else if (trend_spec == "ct") {
            tau_min = TAUMIN_CT;        n_small = 3;
            tau_max = TAUMAX_CT;
            tau_star = TAUSTAR_CT;
            smallp = TAUSMALLP_CT;      n_large = 4;
            largep = TAULARGEP_CT;
        } else {
            throw std::invalid_argument(
                "mackinnon_p_value: trend_spec must be nc/c/ct");
        }
    }

    if (statistic > tau_max[0]) return 1.0;
    if (statistic < tau_min[0]) return 0.0;
    if (statistic > tau_star[0]) {
        return detail::normal_cdf(
            detail::polyval_ascending(largep, n_large, statistic));
    }
    return detail::normal_cdf(
        detail::polyval_ascending(smallp, n_small, statistic));
}

// ---------------------------------------------------------------------------
// KPSS p 值与临界值 (U6, arch kpss_crit 复刻)
//
// arch 算法 (arch/unitroot/unitroot.py::kpss_crit, probe 实测):
//   p = interp(stat, x, y) / 100   — x 为分位数 (升序), y 为百分位 (降序)
//   cv = interp([1, 5, 10], y[::-1], x[::-1]) — 反转后 y 升序
//   越界 clamp (np.interp 行为): stat > x_max → p = y_min/100 = 0.0001
//
// 表: KPSS 1992 100,000,000 次模拟分位数 (Hobijn et al. 1998 报告值)
// 基准: KPSSCRIT_C_P05/CV, KPSSCRIT_CT_P03/CV, KPSSCRIT_C_P10 (1e-12)
// ---------------------------------------------------------------------------

namespace detail {

// 升序 xp 上的线性插值, 越界 clamp (np.interp 兼容; 77 元素线性扫描足够)
inline Real interp_clamp(const Real* xp, const Real* fp, Size n, Real x) {
    if (x <= xp[0]) return fp[0];
    if (x >= xp[n - 1]) return fp[n - 1];
    for (Size i = 1; i < n; ++i) {
        if (xp[i] >= x) {
            const Real t = (x - xp[i - 1]) / (xp[i] - xp[i - 1]);
            return fp[i - 1] + t * (fp[i] - fp[i - 1]);
        }
    }
    return fp[n - 1];  // 不可达
}

}  // namespace detail

// KPSS p 值 (stat ≥ 0; trend 只支持 c/ct)
inline Real kpss_p_value(Real stat, const std::string& trend_spec) {
    const Real* x;
    const Real* y;
    Size n;
    if (trend_spec == "c") {
        y = KPSS_Y_C;
        x = KPSS_X_C;
        n = KPSS_N_C;
    } else if (trend_spec == "ct") {
        y = KPSS_Y_CT;
        x = KPSS_X_CT;
        n = KPSS_N_CT;
    } else {
        throw std::invalid_argument("kpss_p_value: trend_spec must be c/ct");
    }
    // x 升序可直接插; y 降序 (stat 越大 p 越小) 符合方向
    return detail::interp_clamp(x, y, n, stat) / 100.0;
}

// KPSS 临界值 {1%, 5%, 10%} — 在 (y 升序化, x) 上反向插值
inline std::vector<Real> kpss_critical_values(const std::string& trend_spec) {
    const Real* x;
    const Real* y;
    Size n;
    if (trend_spec == "c") {
        y = KPSS_Y_C;
        x = KPSS_X_C;
        n = KPSS_N_C;
    } else if (trend_spec == "ct") {
        y = KPSS_Y_CT;
        x = KPSS_X_CT;
        n = KPSS_N_CT;
    } else {
        throw std::invalid_argument(
            "kpss_critical_values: trend_spec must be c/ct");
    }
    // 反转: y 降序 → 升序 (x 同步反转保持配对)
    std::vector<Real> yr(n), xr(n);
    for (Size i = 0; i < n; ++i) {
        yr[i] = y[n - 1 - i];
        xr[i] = x[n - 1 - i];
    }
    std::vector<Real> out(3);
    const Real levels[3] = {1.0, 5.0, 10.0};
    for (Size j = 0; j < 3; ++j) {
        out[j] = detail::interp_clamp(yr.data(), xr.data(), n, levels[j]);
    }
    return out;
}

}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
