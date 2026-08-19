// =============================================================================
// engle_granger.hpp - Engle-Granger 两步法协整检验 (AEG)
//
// Phase 7C v1.7 M3 (PHASE7C_SPEC.md §5.1)
//
// 教材锚点: Engle-Granger 1987 Econometrica 55(2):251-276;
//   MacKinnon 1994 JBES 12(2):167-176 (p 值) / MacKinnon 2010 (cv 响应面)
//
// 算法 (复刻 statsmodels coint, stattools.py L1702-1839; 容差 1e-10):
//   Step 1: OLS y0 ~ [y1, trend 列] ("n" 无 trend 列; 趋势列张成空间与
//           statsmodels add_trend 等价 — 残差不变)
//   Step 2: 残差 ADF (nc 形式, Schwert maxlag + AIC 自动滞后) → t 统计量
//           (复用 7B adf_test, 已对 statsmodels adfuller 1e-10 验证)
//   Step 3: cv = MacKinnon 2010 响应面, N=2, T_eff = nobs−1 (Stata egranger
//           约定); trend="n" → NaN (statsmodels 同)
//   Step 4: p = MacKinnon 1994 渐近近似 (N=2)
//   共线性保护: R² ≥ 1 − 100·√ε → 统计量 = −∞, p = 0 (statsmodels 同)
//
// 排幻觉点 (spec §9.3):
//   CI1: cv 按协整响应面 (N=2), 非 ADF N=1 表
//   CI2: p (1994 渐近) 与 cv (nobs−1 小样本修正) 不同源 — 小样本可现
//        "过 1% cv 而 p > 1%" (issue #4138), 测试分列断言
//   CI3: 方向依赖 — 本函数做 y0 ← y1; 反方向 = swap 两参重跑 (调用方显式选择;
//        方向无关替代 = phillips_ouliaris Pz)
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/cointegration/mackinnon_coint_cv.hpp"
#include "cpphub/timeseries/unit_root/adf_test.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace cointegration {

/// EG 两步法结果 (§5.1 v1.2 签名)
struct EGResult {
    Real statistic = 0.0;   ///< 第二步 ADF t 统计量 (共线时 = −∞)
    Real p_value = 0.0;     ///< MacKinnon 1994 渐近 p (CI2: 与 cv 不同源)
    Real cv_1pct = 0.0;     ///< MacKinnon 2010 响应面 (nobs−1 修正); "n" → NaN
    Real cv_5pct = 0.0;
    Real cv_10pct = 0.0;
    std::string trend = "c";   ///< "n"/"c"/"ct"/"ctt" (第一步回归形式)
    Size n_obs = 0;            ///< 总观测数
    bool reject_null = false;  ///< statistic < cv_5pct (H0: 无协整)
    std::string summary;
};

namespace detail {

// 第一步 OLS 残差: y0 ~ [y1, det 列]; 返回残差与中心化 R²
inline std::pair<std::vector<Real>, Real> eg_stage1_residuals(
    const std::vector<Real>& y0, const std::vector<Real>& y1,
    const std::string& trend) {
    const Size t = y0.size();
    const Size k_det = (trend == "n") ? 0 : (trend == "c") ? 1
                       : (trend == "ct") ? 2
                                         : 3;  // ctt
    const Size p = 1 + k_det;
    Eigen::MatrixXd X(t, p);
    X.col(0) = Eigen::VectorXd::Map(y1.data(), static_cast<Eigen::Index>(t));
    if (k_det >= 1) X.col(1).setOnes();
    for (Size j = 2; j <= k_det; ++j) {
        for (Size i = 0; i < t; ++i) {
            X(i, j) = std::pow(static_cast<Real>(i),
                               static_cast<Real>(j - 1));
        }
    }
    Eigen::VectorXd yv =
        Eigen::VectorXd::Map(y0.data(), static_cast<Eigen::Index>(t));
    Eigen::VectorXd beta = X.colPivHouseholderQr().solve(yv);
    Eigen::VectorXd resid = yv - X * beta;
    const Real ymean = yv.mean();
    const Real tss = (yv.array() - ymean).square().sum();
    const Real rss = resid.squaredNorm();
    const Real r2 = tss > 0.0 ? 1.0 - rss / tss : 1.0;
    std::vector<Real> out(resid.data(), resid.data() + t);
    return {out, r2};
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Engle-Granger 两步法协整检验
//
// @param y0 左侧变量 (被解释; CI3: 方向依赖 — 反方向请 swap 两参重跑)
// @param y1 右侧变量
// @param trend 第一步协整回归形式 "n"/"c"/"ct"/"ctt" (statsmodels 直译)
//   (跨库对照参数沿用字符串, §1.4-4)
// ---------------------------------------------------------------------------
inline EGResult engle_granger(const std::vector<Real>& y0,
                              const std::vector<Real>& y1,
                              const std::string& trend = "c") {
    if (trend != "n" && trend != "c" && trend != "ct" && trend != "ctt") {
        throw std::invalid_argument(
            "engle_granger: trend must be n/c/ct/ctt");
    }
    if (y0.size() != y1.size()) {
        throw std::invalid_argument("engle_granger: y0/y1 length mismatch");
    }
    const Size t = y0.size();
    if (t < 6) {
        throw std::invalid_argument("engle_granger: sample too small (T < 6)");
    }
    for (Size i = 0; i < t; ++i) {
        if (std::isnan(y0[i]) || std::isnan(y1[i])) {
            throw std::invalid_argument("engle_granger: NaN in input");
        }
    }

    EGResult res;
    res.trend = trend;
    res.n_obs = t;

    auto [resid, r2] = detail::eg_stage1_residuals(y0, y1, trend);

    // 共线性保护 (statsmodels: rsquared < 1 - 100*SQRTEPS 才做 ADF)
    const Real sqrteps =
        std::sqrt(std::numeric_limits<Real>::epsilon());
    if (r2 >= 1.0 - 100.0 * sqrteps) {
        res.statistic = -std::numeric_limits<Real>::infinity();
        res.p_value = 0.0;
    } else {
        // Step 2: 残差 ADF (nc, Schwert + AIC; 复用 7B 引擎)
        const auto adf = unit_root::adf_test(resid, "nc", 0, true);
        res.statistic = adf.statistic;
        res.p_value = mackinnon_coint_p_value(res.statistic, trend, 2);
    }

    // Step 3: 临界值 (N=2, T_eff = nobs−1; "n" → NaN)
    const auto cvs = mackinnon_coint_critical_values(2, trend, t - 1);
    res.cv_1pct = cvs[0];
    res.cv_5pct = cvs[1];
    res.cv_10pct = cvs[2];
    res.reject_null =
        !std::isnan(res.cv_5pct) && res.statistic < res.cv_5pct;

    res.summary = "Engle-Granger AEG test (H0: no cointegration; trend=" +
                  trend + ")";
    return res;
}

}  // namespace cointegration
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
