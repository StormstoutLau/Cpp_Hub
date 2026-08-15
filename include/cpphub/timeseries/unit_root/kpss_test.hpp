// =============================================================================
// kpss_test.hpp - KPSS 平稳性检验 (spec §3.3)
//
// Phase 7B v1.6 M2 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 教材锚点: Kwiatkowski-Phillips-Schmidt-Shin 1992 / Hobijn-Franses-Ooms 1998
// 对照库: Python arch 8.0.0 arch/unitroot/unitroot.py KPSS
//
// 幻觉点防护 (spec §6.2):
//   U10: KPSS 统计量非标准分布 (kpss_crit 100M MC 分位数插值, 非 χ²)
//   U11: 带宽默认 Hobijn et al. 1998 数据依赖法 (arch _autolag);
//        legacy 模式 = Schwert 规则 (显式传 bandwidth 即等价)
//   U12: H0: 平稳性 (与 ADF 相反!), 右尾拒绝
//   U11-kernel: Bartlett 核 (cov_nw), 非 QS 核
//
// arch _autolag 源码 (unitroot.py:1344-1375 实测):
//   covlags = int(n^(2/9))  (截断)
//   s0 = Σu²/n + Σ_{i=1..covlags} 2γ_i,  s1 = Σ_{i=1..covlags} i·2γ_i
//   ŝ = s1/s0,  γ̂ = 1.1447·(ŝ²)^(1/3)
//   lags = min(n, int(γ̂·n^(1/3)))  (截断)
//
// 基准: KPSS_CASES[5] (Hobijn 10/6 + legacy 16) 1e-10 + KPSS_TINY 手算 1e-15
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

// ---------------------------------------------------------------------------
// Hobijn et al. 1998 数据依赖带宽选择 (U11, U-ADR7; arch KPSS _autolag 复刻)
// 假设 Bartlett 核; 公开导出便于独立测试
// @param resids 已退化趋势的残差序列
// @return NW 带宽 (滞后阶数)
// ---------------------------------------------------------------------------
inline Size hobijn_bandwidth(const std::vector<Real>& resids) {
    const Size n = resids.size();
    if (n == 0) {
        throw std::invalid_argument("hobijn_bandwidth: empty residuals");
    }
    const Real nf = static_cast<Real>(n);
    const Size covlags = static_cast<Size>(std::pow(nf, 2.0 / 9.0));  // 截断

    Real s0 = 0.0;
    for (Real u : resids) s0 += u * u;
    s0 /= nf;
    Real s1 = 0.0;
    for (Size i = 1; i <= covlags; ++i) {
        Real cp = 0.0;
        for (Size t = i; t < n; ++t) cp += resids[t] * resids[t - i];
        cp /= nf / 2.0;  // arch: resids_prod /= nobs / 2 (即 2γ_i)
        s0 += cp;
        s1 += static_cast<Real>(i) * cp;
    }
    if (!(s0 > 0.0)) {
        throw std::runtime_error(
            "hobijn_bandwidth: s0 <= 0 (zero or degenerate residuals)");
    }
    const Real s_hat = s1 / s0;
    const Real pwr = 1.0 / 3.0;
    const Real gamma_hat = 1.1447 * std::pow(s_hat * s_hat, pwr);
    Real lags = gamma_hat * std::pow(nf, pwr);
    lags = std::floor(lags);              // int() 截断
    if (lags > nf) lags = nf;             // amin([nobs, ·])
    return static_cast<Size>(lags);
}

struct KPSSResult {
    Real statistic = 0.0;            ///< KPSS LM 统计量 (非标准分布, U10)
    Real p_value = 0.0;              ///< p 值 (arch kpss_crit MC 分位插值)
    Real critical_value_1pct = 0.0;
    Real critical_value_5pct = 0.0;
    Real critical_value_10pct = 0.0;
    Size bandwidth = 0;              ///< Hobijn 带宽 (U-ADR7, Bartlett 核)
    std::string trend_spec;          ///< "c" (level) / "ct" (trend)
    Size n_obs = 0;                  ///< 样本量 T
    bool reject_null = false;        ///< 是否拒绝 H0 (平稳性) — 注意 U12!
    std::string summary;             ///< 明确标注 H0/H1 方向
};

// KPSS 检验
// H0: 平稳性 (level/trend stationary) ← 与 ADF 相反! (U12)
// H1: 单位根 (非平稳)
// @param data 原始序列 (水平值)
// @param trend_spec "c" (level stationary) / "ct" (trend stationary)
// @param bandwidth NW 带宽 (0 => Hobijn et al. 1998 自动; 显式 Schwert 值
//        即 arch legacy 模式)
inline KPSSResult kpss_test(const std::vector<Real>& data,
                            const std::string& trend_spec = "c",
                            Size bandwidth = 0) {
    const Size T = data.size();
    if (T < 3) {
        throw std::invalid_argument("kpss_test: sample too small");
    }
    if (trend_spec != "c" && trend_spec != "ct") {
        throw std::invalid_argument("kpss_test: trend_spec must be c/ct");
    }

    // Step 1: 退化趋势 (OLS: c → 均值, ct → 线性趋势)
    const Size tc = detail::trend_cols(trend_spec);
    std::vector<std::vector<Real>> Z(T);
    for (Size t = 0; t < T; ++t) {
        std::vector<Real> row;
        row.reserve(tc);
        row.push_back(1.0);
        if (tc >= 2) row.push_back(static_cast<Real>(t));
        Z[t] = std::move(row);
    }
    const auto ols = detail::ols_fit(data, Z);
    const std::vector<Real>& u = ols.resid;

    // Step 2: 带宽选择 (U11: 默认 Hobijn; 显式值 = 用户/legacy 指定)
    if (bandwidth == 0) {
        bandwidth = hobijn_bandwidth(u);
    }
    if (bandwidth >= T) {
        throw std::invalid_argument("kpss_test: bandwidth >= nobs");
    }

    // Step 3: 长期方差 (Bartlett 核, U11-kernel) + Step 4: LM 统计量
    const Real lam = long_run_variance(u, bandwidth);
    if (!(lam > 0.0)) {
        throw std::runtime_error("kpss_test: degenerate long-run variance");
    }
    const Real nf = static_cast<Real>(T);
    Real ssum = 0.0;
    Real s = 0.0;
    for (Real e : u) {
        s += e;
        ssum += s * s;
    }
    const Real stat = ssum / (nf * nf) / lam;

    // Step 5: p 值与临界值 (U10: MC 分位数插值)
    KPSSResult res;
    res.statistic = stat;
    res.p_value = kpss_p_value(stat, trend_spec);
    const auto cvs = kpss_critical_values(trend_spec);
    res.critical_value_1pct = cvs[0];
    res.critical_value_5pct = cvs[1];
    res.critical_value_10pct = cvs[2];
    res.bandwidth = bandwidth;
    res.trend_spec = trend_spec;
    res.n_obs = T;
    // H0: 平稳性 — 右尾拒绝 (与 ADF 左尾相反, U12)
    res.reject_null = stat > res.critical_value_5pct;
    res.summary = "KPSS ( trend=" + trend_spec +
                  ", bandwidth=" + std::to_string(bandwidth) +
                  " ) H0: stationary; H1: unit root";
    return res;
}

}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
