// =============================================================================
// ardl_model.hpp - ARDL(p, q_1..q_k) 设计矩阵 + OLS + A8 计账镜像 (spec §2.1)
//
// Phase 8 v1.8 M1 Step 2 (PHASE8_SPEC.md v1.1 §2.1)
//
// 一手源码语义落档 (statsmodels 0.14.6 tsa/ardl/model.py 实录,
// verify_ardl.py [S1]/[S2] 逐值自证 1e-15; R ARDL 交叉待 Step 4):
//   - 模型形式 (causal=False, A6): 同期 x_{j,t} 在回归元内
//   - 设计矩阵列序 (冻结): [trend 列] + [y.L1..y.Lp] + [xj.L0..xj.Lq_j]
//   - trend 列值 = 全样本 1-based 再按 hold_back 切片 (t = hb+1 .. T) —
//     与 statsmodels VAR 的估计样本内重编号不同, 跨模块陷阱 (S1 实证)
//   - hold_back 默认 = max(p, max_j q_j) = _maxlag (L0 不贡献)
//   - OLS: Eigen HouseholderQR 求解; vcov = (X'X)⁻¹·s²_df,
//     s²_df = SSR/(n_act − K) (真实 df, OLS 结果固化 — 与 nobs 覆盖无关)
//   - ⚠️ 计账 quirk (幻觉点 A8, 精确机制): statsmodels 的 res.nobs/
//     sigma2/llf/aic/bic/hqic 用 n_acct = T − hb_early, 其中 hb_early 在
//     AutoReg.__init__ _check_lags 阶段固化 = 显式 hb 或 max(y_lags)=p
//     (仅 y 滞后, exog 滞后尚未并入), 之后不回写 → 当 max_j q_j > p 或
//     显式 hb > p 时 n_acct ≠ n_act = T − hold_back_used:
//       llf  = −n_acct/2·(ln 2π + ln(SSR/n_acct) + 1)
//       sigma2 = SSR/n_acct;  IC 惩罚项恒 (K+1) (σ² 计入):
//       aic  = −2·llf + 2·(K+1)
//       bic  = −2·llf + ln(n_acct)·(K+1)
//       hqic = −2·llf + 2·ln·ln(n_acct)·(K+1)
//     本实现逐式镜像 (Step 4: IC 公式差异属阻断性断言)
//   - "ctt" 档 statsmodels ARDL 不支持 (仅 n/c/t/ct) — 本库支持 (SPEC
//     裁决 #2), 基准 = 手工 numpy OLS (构造式直接保证); bounds 拒绝 ctt
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace ardl {

/// ARDL 估计设定 (spec §2.1)
struct ArdlSpec {
    Size p = 1;                    ///< y 自身滞后阶
    std::vector<Size> q;           ///< 各 force 变量滞后阶 (含零阶; k = q.size(), 每变量独立 q_j)
    std::string trend = "c";       ///< "n"/"c"/"ct"/"ctt" (statsmodels 直译 + ctt 扩展)
    bool causal = false;           ///< v1.8 恒 false (A6); 传 true → throw (v1.9 落地)
    Size hold_back = 0;            ///< 0 → 默认 max(p, max_j q_j) (statsmodels 同语义)
};

/// ARDL 估计结果 (spec §2.1 逐字段)
struct ArdlResult {
    std::vector<Real> params;      ///< OLS 系数 (列序 = 设计矩阵冻结序)
    std::vector<std::vector<Real>> vcov;  ///< (X'X)⁻¹σ̂²_df, 同列序 (bse = √diag)
    Real sigma2 = 0.0;             ///< SSR/n_acct (A8 计账镜像, = statsmodels res.sigma2)
    Real loglik = 0.0;             ///< A8 口径: T = n_acct (见头注释)
    Real aic = 0.0, bic = 0.0, hqic = 0.0;  ///< 镜像 statsmodels OLS 口径 (σ² 计入, 恒 K+1)
    std::vector<Real> fitted, residuals;    ///< 估计样本内 (t = hb .. T−1)
    Size n_obs_used = 0;           ///< n_act = T − hold_back_used (真实估计样本)
    Size hold_back_used = 0;
    ArdlSpec spec{};               ///< 回显 (含默认化后的 hold_back)
};

namespace detail {

/// trend → 确定性列数 ("ctt" 为本库扩展, statsmodels ARDL 无此档)
inline Size trend_order(const std::string& trend) {
    if (trend == "n") return 0;
    if (trend == "c") return 1;
    if (trend == "ct") return 2;
    if (trend == "ctt") return 3;
    throw std::invalid_argument("ardl_fit: trend must be one of n/c/ct/ctt");
}

}  // namespace detail

/// ARDL(p, q_1..q_k) OLS 估计 (causal=False 含同期项, A6)
///
/// x 为 k 个等长 force 变量; k=0 (纯 AR) 允许拟合, 但 bounds 检验拒绝
/// (PSS 响应面表无 k=0 键, k 域 = 1..10)
inline ArdlResult ardl_fit(const std::vector<Real>& y,
                           const std::vector<std::vector<Real>>& x,
                           const ArdlSpec& spec = ArdlSpec{}) {
    // ---- 校验 ----
    if (spec.causal) {
        throw std::invalid_argument(
            "ardl_fit: causal=true not supported in v1.8 (A6)");
    }
    const Size kt = detail::trend_order(spec.trend);
    const Size T = y.size();
    const Size k = x.size();
    if (T == 0) throw std::invalid_argument("ardl_fit: empty endog");
    if (k > 10) {
        throw std::invalid_argument("ardl_fit: k > 10 (beyond PSS surface)");
    }
    if (spec.q.size() != k) {
        throw std::invalid_argument("ardl_fit: q.size() must equal x.size()");
    }
    for (Size j = 0; j < k; ++j) {
        if (x[j].size() != T) {
            throw std::invalid_argument("ardl_fit: x[j] length != y length");
        }
        for (Real v : x[j]) {
            if (!std::isfinite(v)) {
                throw std::invalid_argument("ardl_fit: NaN/inf in exog");
            }
        }
    }
    for (Real v : y) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("ardl_fit: NaN/inf in endog");
        }
    }

    // ---- Step 1: hold_back_eff ----
    Size maxq = 0;
    for (Size qj : spec.q) maxq = std::max(maxq, qj);
    // 显式 hb < max(p, maxq) → 负索引风险; 镜像 statsmodels L586 校验
    if (spec.hold_back != 0 && spec.hold_back < std::max(spec.p, maxq)) {
        throw std::invalid_argument(
            "ardl_fit: hold_back must be >= max(p, max q_j)");
    }
    const Size hb_eff =
        (spec.hold_back != 0) ? spec.hold_back : std::max(spec.p, maxq);
    const Size n_act = T - hb_eff;

    // A8 计账: hb_early 在 statsmodels _check_lags 阶段固化 (仅 y 滞后)
    const Size hb_early = (spec.hold_back != 0) ? spec.hold_back : spec.p;
    const Size n_acct = T - hb_early;

    // ---- Step 2: 设计矩阵 (列序冻结) ----
    Size npar = kt + spec.p;
    for (Size qj : spec.q) npar += qj + 1;
    if (n_act <= npar) {
        throw std::invalid_argument(
            "ardl_fit: sample too small after hold_back");
    }

    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(n_act, npar);
    Eigen::VectorXd Y = Eigen::VectorXd::Zero(n_act);
    Size col = 0;
    if (spec.trend == "c" || spec.trend == "ct" || spec.trend == "ctt") {
        for (Size i = 0; i < n_act; ++i) X(i, col) = 1.0;
        ++col;
    }
    if (spec.trend == "ct" || spec.trend == "ctt") {
        // 全样本 1-based 再切片: t = hb_eff+1 .. T (S1 实证, ≠ VAR 约定)
        for (Size i = 0; i < n_act; ++i) {
            X(i, col) = static_cast<Real>(hb_eff + 1 + i);
        }
        ++col;
        if (spec.trend == "ctt") {
            for (Size i = 0; i < n_act; ++i) {
                const Real t = static_cast<Real>(hb_eff + 1 + i);
                X(i, col) = t * t;
            }
            ++col;
        }
    }
    for (Size l = 1; l <= spec.p; ++l) {  // y.L1 .. y.Lp
        for (Size i = 0; i < n_act; ++i) {
            X(i, col) = y[hb_eff + i - l];
        }
        ++col;
    }
    for (Size j = 0; j < k; ++j) {  // xj.L0 .. xj.Lq_j
        for (Size l = 0; l <= spec.q[j]; ++l) {
            for (Size i = 0; i < n_act; ++i) {
                X(i, col) = x[j][hb_eff + i - l];
            }
            ++col;
        }
    }
    for (Size i = 0; i < n_act; ++i) Y(i) = y[hb_eff + i];

    // ---- Step 3: OLS + vcov ----
    const Eigen::VectorXd beta = X.householderQr().solve(Y);
    const Eigen::VectorXd resid = Y - X * beta;
    const Real ssr = resid.squaredNorm();
    const Real s2_df = ssr / static_cast<Real>(n_act - npar);
    const Eigen::MatrixXd XtX_inv = (X.transpose() * X).lu().inverse();

    // ---- Step 4: A8 计账镜像 (llf/IC/sigma2 用 n_acct) ----
    const Real n_a = static_cast<Real>(n_acct);
    const Real Kp1 = static_cast<Real>(npar) + 1.0;
    const Real llf =
        -n_a / 2.0 * (std::log(2.0 * std::acos(-1.0)) + std::log(ssr / n_a) + 1.0);

    ArdlResult r;
    r.params.resize(npar);
    for (Size i = 0; i < npar; ++i) r.params[i] = beta(i);
    r.vcov.assign(npar, std::vector<Real>(npar, 0.0));
    for (Size i = 0; i < npar; ++i) {
        for (Size j = 0; j < npar; ++j) r.vcov[i][j] = s2_df * XtX_inv(i, j);
    }
    r.sigma2 = ssr / n_a;
    r.loglik = llf;
    r.aic = -2.0 * llf + 2.0 * Kp1;
    r.bic = -2.0 * llf + std::log(n_a) * Kp1;
    r.hqic = -2.0 * llf + 2.0 * std::log(std::log(n_a)) * Kp1;
    r.fitted.resize(n_act);
    r.residuals.resize(n_act);
    for (Size i = 0; i < n_act; ++i) {
        r.fitted[i] = (X.row(i) * beta)(0);
        r.residuals[i] = resid(i);
    }
    r.n_obs_used = n_act;
    r.hold_back_used = hb_eff;
    r.spec = spec;
    r.spec.hold_back = hb_eff;
    return r;
}

}  // namespace ardl
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
