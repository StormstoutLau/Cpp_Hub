// =============================================================================
// var_select.hpp - VAR 滞后阶选择: IC 同样本 offset 搜索 (spec §4.1, V5)
//
// Phase 7C v1.7 M2 (PHASE7C_SPEC.md v1.2 §4.1)
//
// 一手源码语义落档 (statsmodels select_order L778-830):
//   - 候选 p = p_min..max_lag, p_min = 0 当 trend ≠ "n" (纯截距模型参与),
//     trend = "n" 时 p_min = 1
//   - 同样本: 对每个 p 用 offset = max_lag − p 截断 (V5 — 否则 IC 因
//     样本量不同而错序)
//   - max_lag 默认 round(12·(T/100)^{1/4}); max_estimable =
//     (T − K − ntrend)/(1 + K)
//   - selected = argmin(IC) + p_min
//   ⚠️ R vars::VARselect 轨迹从 p=1 起 (无 p=0), 其余逐位一致 (实测
//     verify_var.R: AIC 轨迹 p=1..4 与 statsmodels p=1..4 逐位同)
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/var/multivariate_data.hpp"
#include "cpphub/timeseries/var/var_model.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace var {

/// 滞后选择结果 (v1.2)
struct VARSelectResult {
    Size selected_lag = 0;
    std::string ic_used;
    std::vector<Real> aic, bic, hqic, fpe;  ///< 逐候选 p 全轨迹 (同样本, V5)
    Size n_obs_used = 0;
};

/// IC 滞后阶选择 (同样本 offset, V5). ic ∈ {"aic","bic","hqic","fpe"}
inline VARSelectResult var_select_order(const MultivariateTSData& data,
                                        const std::string& trend = "c",
                                        Size max_lag = 0,
                                        const std::string& ic = "aic") {
    data.validate();
    if (ic != "aic" && ic != "bic" && ic != "hqic" && ic != "fpe") {
        throw std::invalid_argument("var_select_order: ic must be aic/bic/hqic/fpe");
    }
    const Eigen::MatrixXd y = data.matrix();
    const Size T = y.rows(), K = y.cols();
    const Size kt = detail::trend_order(trend);

    if (max_lag == 0) {
        max_lag = static_cast<Size>(
            std::llround(12.0 * std::pow(static_cast<Real>(T) / 100.0, 0.25)));
    }
    const Size ntrend = (trend[0] == 'c') ? trend.size() : 0;
    const Size max_est = (T - K - ntrend) / (1 + K);
    if (max_lag > max_est) {
        throw std::invalid_argument(
            "var_select_order: max_lag too large for T and K");
    }
    if (max_lag < 1) throw std::invalid_argument("var_select_order: max_lag < 1");

    const Size p_min = (trend == "n") ? 1 : 0;
    VARSelectResult res;
    res.ic_used = ic;
    const Size n_cand = max_lag - p_min + 1;
    for (Size c = 0; c < n_cand; ++c) {
        const Size p = p_min + c;
        const Size off = max_lag - p;  // V5 同样本
        Eigen::MatrixXd Y, Z;
        const Eigen::MatrixXd ysub = y.bottomRows(T - off);
        detail::build_var_design(ysub, p, trend, Y, Z);
        const Size nobs = Y.rows();
        const Eigen::MatrixXd B = Z.householderQr().solve(Y);
        const Eigen::MatrixXd R = Y - Z * B;
        Eigen::MatrixXd smle = R.transpose() * R / static_cast<Real>(nobs);
        auto s = detail::info_criteria(smle, nobs, K, p, kt);
        res.aic.push_back(s.aic);
        res.bic.push_back(s.bic);
        res.hqic.push_back(s.hqic);
        res.fpe.push_back(s.fpe);
        res.n_obs_used = nobs;
    }

    const std::vector<Real>& traj = (ic == "aic")   ? res.aic
                                   : (ic == "bic")  ? res.bic
                                   : (ic == "hqic") ? res.hqic
                                                    : res.fpe;
    Size argmin = 0;
    Real best = std::numeric_limits<Real>::infinity();
    for (Size i = 0; i < traj.size(); ++i) {
        if (traj[i] < best) {
            best = traj[i];
            argmin = i;
        }
    }
    res.selected_lag = p_min + argmin;
    return res;
}

}  // namespace var
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
