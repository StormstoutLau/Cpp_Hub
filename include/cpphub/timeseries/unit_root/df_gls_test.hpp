// =============================================================================
// df_gls_test.hpp - DF-GLS 单位根检验 (Elliott-Rothenberg-Stock 1996, spec §3.4)
//
// Phase 7B v1.6 M2 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 教材锚点: Elliott-Rothenberg-Stock 1996 / Perron-Qu 2007 (lag 选择)
// 对照库: Python arch 8.0.0 arch/unitroot/unitroot.py DFGLS
//
// 幻觉点防护 (spec §6.2):
//   U7: 临界值来自 arch 独立模拟表 (非 ERS 1996 原表)
//   U8: c̄ = -7.0 ("c" demean) / -13.5 ("ct" trend)
//   U9: GLS 变换用 ρ̄ = 1 + c̄/T (非 c̄ 本身); 第一项不变换 (arch:947-948)
//   U20: lag 选择 AIC (非 MAIC; arch:960-963), 用 OLS-detrended 序列 +
//        trend "n" (Perron-Qu 2007)
//   U21: p 值用 df_gls 独立 response surface (MacKinnon/arch 模拟)
//
// arch DFGLS._compute_statistic (unitroot.py:936-985 实测):
//   ct = c̄/nobs; dz[0]=z[0], dz[t]=z[t]-(1+ct)·z[t-1]; dy 同理
//   detrend_coef = lstsq(dz, dy); ỹ = y - z·coef
//   lags: _df_select_lags(y_ols_detrend, "n", max_lags, "aic")
//   Δỹ = γ·ỹ_{t-1} + Σδᵢ·Δỹ_{t-i} (trend "n", 无常数)
//   CV 用回归 nobs = (T-1) - lags
//
// 基准: DFGLS_CASES[4] (AIC lag 5/0) 1e-10
// =============================================================================
#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"
#include "cpphub/timeseries/unit_root/mackinnon_cv.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace unit_root {

struct DFGlsResult {
    Real statistic = 0.0;            ///< DF-GLS τ 统计量
    Real p_value = 0.0;              ///< p 值 (df_gls 独立 response surface, U21)
    Real critical_value_1pct = 0.0;
    Real critical_value_5pct = 0.0;
    Real critical_value_10pct = 0.0;
    Size n_lags = 0;                 ///< AIC 选择的 lag (U20, 非 MAIC)
    std::string trend_spec;          ///< "c" (demean) / "ct" (trend)
    Real c_bar = 0.0;                ///< c̄ = -7.0 ("c") / -13.5 ("ct") (U8)
    Real rho_bar = 0.0;              ///< ρ̄ = 1 + c̄/T (U9-修正)
    Size n_obs = 0;                  ///< 回归 nobs = (T-1) - n_lags
    bool reject_null = false;        ///< 是否拒绝 H0 (单位根)
    std::string summary;
};

// DF-GLS 检验
// H0: γ = 0 (单位根, 非平稳); H1: γ < 0 (平稳)
// @param data 原始序列 (水平值)
// @param trend_spec "c" (demean, c̄=-7.0) / "ct" (trend, c̄=-13.5)
// @param max_lag 最大滞后 (0 => AIC 自动, Schwert 候选上限 + arch cap)
inline DFGlsResult df_gls_test(const std::vector<Real>& data,
                               const std::string& trend_spec = "ct",
                               Size max_lag = 0) {
    const Size T = data.size();
    if (T < 5) {
        throw std::invalid_argument("df_gls_test: sample too small");
    }
    if (trend_spec != "c" && trend_spec != "ct") {
        throw std::invalid_argument("df_gls_test: trend_spec must be c/ct");
    }
    const Size tc = detail::trend_cols(trend_spec);

    // Step 1: GLS detrending (U9, ERS 1996)
    // 1.1 c̄ 与 ρ̄ = 1 + c̄/T (U8/U-ADR5)
    const Real c_bar = (trend_spec == "c") ? -7.0 : -13.5;
    const Real rho = 1.0 + c_bar / static_cast<Real>(T);

    // 1.2 z 趋势列 (c → [1], ct → [1, t]; ct 含常数 → t 平移不变)
    std::vector<std::vector<Real>> z(T);
    for (Size t = 0; t < T; ++t) {
        std::vector<Real> row;
        row.reserve(tc);
        row.push_back(1.0);
        if (tc >= 2) row.push_back(static_cast<Real>(t));
        z[t] = std::move(row);
    }
    // 1.3 GLS 变换 (第一项不变换, arch:945-948)
    std::vector<std::vector<Real>> dz(T);
    dz[0] = z[0];
    for (Size t = 1; t < T; ++t) {
        std::vector<Real> r(tc);
        for (Size j = 0; j < tc; ++j) r[j] = z[t][j] - rho * z[t - 1][j];
        dz[t] = std::move(r);
    }
    std::vector<Real> dy(T);
    dy[0] = data[0];
    for (Size t = 1; t < T; ++t) dy[t] = data[t] - rho * data[t - 1];

    // 1.4 OLS: dy ~ dz → δ̂ (arch pinv(lstsq), 列满秩时等价)
    const auto gls_fit = detail::ols_fit(dy, dz);

    // 1.5 GLS 退化: ỹ = y - z·δ̂
    std::vector<Real> yd(T);
    for (Size t = 0; t < T; ++t) {
        Real pred = 0.0;
        for (Size j = 0; j < tc; ++j) pred += gls_fit.beta[j] * z[t][j];
        yd[t] = data[t] - pred;
    }

    // Step 2: lag 选择 (U20: AIC, OLS-detrended 序列, trend "n")
    const auto ols_detrend = detail::ols_fit(data, z);
    std::vector<Real> yol(T);
    for (Size t = 0; t < T; ++t) {
        Real pred = 0.0;
        for (Size j = 0; j < tc; ++j) pred += ols_detrend.beta[j] * z[t][j];
        yol[t] = data[t] - pred;
    }
    const Size p = (max_lag == 0) ? select_lag_by_ic(yol, "n", 0, "aic")
                                  : max_lag;
    if (T - 1 <= p) {
        throw std::invalid_argument("df_gls_test: lag too large for sample");
    }

    // Step 3: DF-GLS 回归 (trend "n": Δỹₜ = γ·ỹₜ₋₁ + Σ δᵢ·Δỹₜ₋ᵢ + εₜ)
    const Size nobs = (T - 1) - p;
    std::vector<Real> lhs(nobs);
    std::vector<std::vector<Real>> X(nobs);
    for (Size i = 0; i < nobs; ++i) {
        lhs[i] = yd[p + i + 1] - yd[p + i];
        std::vector<Real> row;
        row.reserve(1 + p);
        row.push_back(yd[p + i]);  // ỹ_{t-1} → γ 位置 0 (trend "n")
        for (Size l = 1; l <= p; ++l) {
            row.push_back(yd[p + i + 1 - l] - yd[p + i - l]);
        }
        X[i] = std::move(row);
    }
    const auto ols = detail::ols_fit(lhs, X);
    const Real tau = ols.beta[0] / ols.bse[0];

    // Step 4: 临界值与 p 值 (U7/U21: df_gls 独立表, T = 回归 nobs)
    DFGlsResult res;
    res.statistic = tau;
    res.p_value = mackinnon_p_value(tau, "df_gls", trend_spec, nobs, 1);
    res.critical_value_1pct =
        mackinnon_critical_value("df_gls", trend_spec, nobs, 1, 0.01);
    res.critical_value_5pct =
        mackinnon_critical_value("df_gls", trend_spec, nobs, 1, 0.05);
    res.critical_value_10pct =
        mackinnon_critical_value("df_gls", trend_spec, nobs, 1, 0.10);
    res.n_lags = p;
    res.trend_spec = trend_spec;
    res.c_bar = c_bar;
    res.rho_bar = rho;
    res.n_obs = nobs;
    // H0: 单位根; 左尾拒绝
    res.reject_null = tau < res.critical_value_5pct;
    res.summary = "DF-GLS ( trend=" + trend_spec +
                  ", lags=" + std::to_string(p) + " ) H0: unit root; "
                  "H1: stationary";
    return res;
}

}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
