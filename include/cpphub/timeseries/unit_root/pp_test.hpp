// =============================================================================
// pp_test.hpp - Phillips-Perron 单位根检验 (spec §3.2)
//
// Phase 7B v1.6 M2 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 教材锚点: Phillips-Perron 1988 / Hamilton 1994 Ch17 eq.17.6.17
// 对照库: Python arch 8.0.0 arch/unitroot/unitroot.py PhillipsPerron
//
// 幻觉点防护 (spec §6.2):
//   U5: PP 带宽用 Schwert 规则 ceil(12·(T/100)^0.25), 与 ADF 相同公式;
//       但 PP 无 ADF 的 min 上限保护 (arch unitroot.py:1135-1136 直接赋值)
//   U6: PP 用 Bartlett 核 (cov_nw), 非 QS 核
//   U3/U4: MacKinnon 2010 (PP 共享 ADF 表), τ 非标准分布
//   U5-sigma2 (三次修正, arch unitroot.py:1164-1166 实测):
//     第一项 ratio 用 gamma0 = SSR/n (无 df 修正);
//     第二项用 s = √(SSR/(n-k)) (df-corrected)
//     Z(tau) = sqrt(gamma0/lam2)·t_γ - 0.5·((lam2-gamma0)/lam)·(n·SE(γ)/s)
//     (第二项分母为 lam=√lam2 乘进系数, 非 2·σ²/2·γ̂/2·σ·σ_ε)
//
// 基准: PP_CASES[5] (250 样本 rw/ar × c/ct + 显式 lags=8) 1e-12
//       + PP_TINY 手算分量 (y=[1..10], trend n, lags=2) 1e-15
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

struct PPResult {
    Real statistic = 0.0;            ///< PP Z(tau) 统计量 (非标准分布)
    Real p_value = 0.0;              ///< p 值 (MacKinnon 2010, 共享 ADF 表)
    Real critical_value_1pct = 0.0;
    Real critical_value_5pct = 0.0;
    Real critical_value_10pct = 0.0;
    Size bandwidth = 0;              ///< NW 带宽 (U5/U-ADR8)
    std::string trend_spec;          ///< "nc"/"c"/"ct"
    Size n_obs = 0;                  ///< 有效观测数 (T-1, 无滞后损失)
    bool reject_null = false;        ///< 是否拒绝 H0 (单位根)
    std::string summary;             ///< 摘要 (H0/H1 方向)
};

// PP 检验
// H0: γ = 0 (单位根, 非平稳); H1: γ < 0 (平稳)
// @param data 原始序列 (水平值)
// @param trend_spec "auto"/"nc"/"c"/"ct"
// @param bandwidth NW 带宽 (0 => Schwert 自动, 无上限保护)
inline PPResult pp_test(const std::vector<Real>& data,
                        const std::string& trend_spec = "auto",
                        Size bandwidth = 0) {
    const Size T = data.size();
    if (T < 4) {
        throw std::invalid_argument("pp_test: sample too small");
    }

    // Step 1: 方程形式选择 (同 ADF)
    std::string trend = trend_spec;
    if (trend == "auto") {
        trend = select_trend_spec(data);
    } else if (trend != "nc" && trend != "n" && trend != "c" &&
               trend != "ct") {
        throw std::invalid_argument("pp_test: trend_spec must be auto/nc/c/ct");
    }
    const Size tc = detail::trend_cols(trend);

    // Step 2: 带宽选择 (U5, 无 ADF 式上限保护)
    if (bandwidth == 0) {
        bandwidth = schwert_lag(T);
    }

    // Step 3: PP 回归 (无滞后差分项, 仅 OLS)
    //   Δyₜ = [trend 列] + γ·yₜ₋₁ + εₜ, nobs = T-1
    const Size nobs = T - 1;
    if (bandwidth >= nobs) {
        throw std::invalid_argument("pp_test: bandwidth >= nobs");
    }
    std::vector<Real> lhs(nobs);
    std::vector<std::vector<Real>> X(nobs);
    for (Size i = 0; i < nobs; ++i) {
        // 时刻 t = i + 1 (Δy_t), y_{t-1} = data[i]
        lhs[i] = data[i + 1] - data[i];
        std::vector<Real> row;
        row.reserve(tc + 1);
        if (tc >= 1) row.push_back(1.0);
        if (tc >= 2) row.push_back(static_cast<Real>(i));
        row.push_back(data[i]);  // y_{t-1} → γ 位置 = tc
        X[i] = std::move(row);
    }
    const auto ols = detail::ols_fit(lhs, X);
    const Real gamma = ols.beta[tc];
    const Real sigma = ols.bse[tc];          // SE(γ̂)
    const Real t_gamma = gamma / sigma;      // t 统计量

    // Step 4: PP 修正 (U5/U6/U5-sigma2)
    const Real n = static_cast<Real>(nobs);
    const Real k = static_cast<Real>(ols.n_params);
    const Real gamma0 = ols.ssr / n;                          // SSR/n (无 df)
    const Real lam2 = long_run_variance(ols.resid, bandwidth);  // Bartlett
    const Real lam = std::sqrt(lam2);
    if (!(gamma0 > 0.0) || !(lam2 > 0.0)) {
        throw std::runtime_error(
            "pp_test: degenerate variance (zero residuals?)");
    }
    const Real s = std::sqrt(ols.ssr / (n - k));  // df-corrected
    const Real stat = std::sqrt(gamma0 / lam2) * t_gamma -
                      0.5 * ((lam2 - gamma0) / lam) * (n * sigma / s);

    // Step 5: 临界值与 p 值 (MacKinnon 2010, PP 共享 ADF 表, T = nobs)
    PPResult res;
    res.statistic = stat;
    res.p_value = mackinnon_p_value(stat, "pp", trend, nobs, 1);
    res.critical_value_1pct =
        mackinnon_critical_value("pp", trend, nobs, 1, 0.01);
    res.critical_value_5pct =
        mackinnon_critical_value("pp", trend, nobs, 1, 0.05);
    res.critical_value_10pct =
        mackinnon_critical_value("pp", trend, nobs, 1, 0.10);
    res.bandwidth = bandwidth;
    res.trend_spec = trend;
    res.n_obs = nobs;
    // H0: 单位根; 左尾拒绝
    res.reject_null = stat < res.critical_value_5pct;
    res.summary = "PP ( trend=" + trend +
                  ", bandwidth=" + std::to_string(bandwidth) +
                  " ) H0: unit root; H1: stationary";
    return res;
}

}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
