// =============================================================================
// adf_test.hpp - ADF 单位根检验 (spec §3.1)
//
// Phase 7B v1.6 M2 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 教材锚点: Dickey-Fuller 1979 / Said-Dickey 1984 / Hamilton 1994 Ch 17
// 对照库: Python arch 8.0.0 arch/unitroot/unitroot.py ADF (U-ADR9)
//
// 幻觉点防护 (spec §6.2):
//   U1: Schwert lag ceil(12·(T/100)^0.25) + arch 上限保护
//   U3: MacKinnon 2010 临界值, T 用 nobs (有效观测数, 非原始 T)
//   U4: τ 非标准分布 (p 值用 mackinnon_p_value, 非 Student-t)
//   AIC 模式: 评估样本固定截断 (T-1-max_lag), 最终回归用 (T-1-p) 样本
//     (ADF_AIC 基准 nobs=244/249 实证)
//
// 基准: ADF_FIXED[6] (lag=16, nobs=233) + ADF_AIC[4] (n_lags=5/0) 1e-12
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

struct ADFResult {
    Real statistic = 0.0;            ///< ADF τ 统计量 (U4: 非标准分布)
    Real p_value = 0.0;              ///< p 值 (MacKinnon 2010, U3)
    Real critical_value_1pct = 0.0;
    Real critical_value_5pct = 0.0;
    Real critical_value_10pct = 0.0;
    Size n_lags = 0;                 ///< 使用的滞后数
    std::string trend_spec;          ///< "nc"/"c"/"ct"
    std::vector<Real> coefficients;  ///< 回归系数 [trend..., γ, δ₁..δ_p]
    std::vector<Real> std_errors;    ///< 标准误
    Size n_obs = 0;                  ///< 有效观测数
    bool reject_null = false;        ///< 是否拒绝 H0 (单位根)
    std::string summary;             ///< 摘要 (H0/H1 方向)
};

// ADF 检验
// H0: γ = 0 (单位根, 非平稳); H1: γ < 0 (平稳)
// @param data 原始序列 (水平值)
// @param trend_spec "auto"/"nc"/"c"/"ct"
// @param max_lag 最大滞后 (0 => Schwert 自动 + arch 上限保护)
// @param use_aic_bic true => AIC 选择 lag (评估样本固定截断);
//                    false => 直接用 max_lag
inline ADFResult adf_test(const std::vector<Real>& data,
                          const std::string& trend_spec = "auto",
                          Size max_lag = 0,
                          bool use_aic_bic = false) {
    const Size T = data.size();
    if (T < 4) {
        throw std::invalid_argument("adf_test: sample too small");
    }

    // Step 1: 方程形式 (U2)
    std::string trend = trend_spec;
    if (trend == "auto") {
        trend = select_trend_spec(data);
    } else if (trend != "nc" && trend != "n" && trend != "c" &&
               trend != "ct") {
        throw std::invalid_argument(
            "adf_test: trend_spec must be auto/nc/c/ct");
    }
    const Size tc = detail::trend_cols(trend);

    // Step 2: lag 选择 (U1)
    if (max_lag == 0) {
        Size max_max = (T - 1) / 2 > 1 ? (T - 1) / 2 - 1 : 0;
        if (tc > max_max) max_max = 0; else max_max -= tc;
        max_lag = schwert_lag(T);
        if (max_lag > max_max) max_lag = max_max;
    }
    Size p = max_lag;
    if (use_aic_bic) {
        p = select_lag_by_ic(data, trend, max_lag, "aic");
    }
    if (T - 1 <= p) {
        throw std::invalid_argument("adf_test: lag too large for sample");
    }

    // Step 3: 最终 ADF 回归 (样本 T-1-p, 与 arch 一致)
    //   Δyₜ = [trend 列] + γ·yₜ₋₁ + Σ δᵢ·Δyₜ₋ᵢ + εₜ
    const Size nobs = (T - 1) - p;
    std::vector<Real> lhs(nobs);
    std::vector<std::vector<Real>> X(nobs);
    for (Size i = 0; i < nobs; ++i) {
        // 时刻 t = p + i + 1 (Δy_t), y_{t-1} = data[p+i]
        lhs[i] = data[p + i + 1] - data[p + i];
        std::vector<Real> row;
        row.reserve(tc + 1 + p);
        if (tc >= 1) row.push_back(1.0);
        if (tc >= 2) row.push_back(static_cast<Real>(i));
        row.push_back(data[p + i]);  // y_{t-1} → γ 位置 = tc
        for (Size l = 1; l <= p; ++l) {
            row.push_back(data[p + i + 1 - l] - data[p + i - l]);
        }
        X[i] = std::move(row);
    }
    const auto ols = detail::ols_fit(lhs, X);
    const Real gamma = ols.beta[tc];
    const Real gamma_se = ols.bse[tc];
    const Real tau = gamma / gamma_se;  // ADF τ (U4)

    // Step 4: MacKinnon 2010 (U3, T = nobs)
    ADFResult res;
    res.statistic = tau;
    res.p_value = mackinnon_p_value(tau, "adf", trend, nobs, 1);
    res.critical_value_1pct =
        mackinnon_critical_value("adf", trend, nobs, 1, 0.01);
    res.critical_value_5pct =
        mackinnon_critical_value("adf", trend, nobs, 1, 0.05);
    res.critical_value_10pct =
        mackinnon_critical_value("adf", trend, nobs, 1, 0.10);
    res.n_lags = p;
    res.trend_spec = trend;
    res.coefficients = ols.beta;
    res.std_errors = ols.bse;
    res.n_obs = nobs;
    // Step 5: H0: 单位根; 左尾拒绝
    res.reject_null = tau < res.critical_value_5pct;
    res.summary = std::string("ADF (") + trend +
                  (use_aic_bic ? ", AIC lag)" : ", fixed lag)") +
                  " H0: unit root; H1: stationary";
    return res;
}

}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
