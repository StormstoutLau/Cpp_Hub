// =============================================================================
// var_model.hpp - VAR(p) 逐方程 OLS 估计 + IC 五式 + 稳定性 (spec §4.1)
//
// Phase 7C v1.7 M2 (PHASE7C_SPEC.md v1.2 §4.1, 决策 9)
//
// 一手源码语义落档 (statsmodels 0.14.4 var_model.py 实录, 与 spec 冻结
// 0.14.6 L2281-2303 逐字一致; R vars 交叉一致):
//   - 设计矩阵 Z (get_var_endog): 列序 [trend(n/c/ct), y1.l1..yK.l1, y1.l2..]
//     lag-major — params 布局 K×(K·p + k_trend), 行=方程, 列=回归元;
//     与 R vars coef 顺序一致 (y1.l1, y2.l1, y3.l1, y1.l2, ...) 但 trend
//     位置: statsmodels trend 列在最前, vars const 在最后 — C++ 采用
//     statsmodels 布局 (主基准), 对照 vars 时仅平移列
//   - 逐方程 OLS: 同阶同回归元下与系统 GLS 数值等价 (决策 9;
//     statsmodels VAR.fit 的 method 参数未被引用, B3 双盲)
//   - Σ 分层 (V4): sigma_u = SSR/(T_eff − K·p − k_trend) (df 修正, IRF/
//     FEVD/GFEVD 用此, 对齐 statsmodels _chol_sigma_u 与 R covres);
//     sigma_u_mle = SSR/T_eff (IC/loglik 用此)
//   - IC 五式 (L2281-2303): fp = p·K² + K·k_trend; ld = logdet(Σ_mle);
//     aic = ld + 2/T·fp; bic = ld + ln(T)/T·fp; hqic = ld + 2·lnln(T)/T·fp;
//     fpe = ((T + df_model)/df_resid)^K · exp(ld), df_model = p·K + k_trend
//     (单方程, V6 指数 K)
//   - loglik (var_loglike L305-337): −(T·K/2)ln(2π) − (T/2)(ld + K)
//   - 稳定性 (V9): 伴随矩阵特征值; statsmodels roots = 1/eig (模 >1 稳定),
//     C++ 双输出 max_abs_eigenvalue (<1) + is_strictly_stationary (严格<1,
//     与 statsmodels is_stable 的 ≤ 容忍口径区分)
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/var/multivariate_data.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace var {

/// VAR 估计设定 (spec §4.1)
struct VARSpec {
    Size lag = 0;                       ///< 0 => 由 IC (aic) 选择
    std::string trend = "c";            ///< "n"/"c"/"ct" (statsmodels 直译)
    Size max_lag = 0;                   ///< IC 搜索上限 (0 => 12·(T/100)^{1/4})
    bool same_sample_ic = true;         ///< V5: offset 强制同一样本量
    Eigen::MatrixXd identification_P;   ///< P 注入 (决策 11): 0×0 => Cholesky
                                        ///< LLT 默认; 注入须下三角且 P·P'≈Σ
};

/// VAR 估计结果 (spec §4.1 逐字段)
struct VARResult {
    Eigen::MatrixXd coefficients;       ///< K×(K·p + k_trend) 行=方程
    Eigen::MatrixXd sigma_u_mle;        ///< ML 版 Σ = SSR/T (V4: ÷T)
    Real loglik = 0.0;
    Real aic = 0.0, bic = 0.0, hqic = 0.0, fpe = 0.0, logdet = 0.0;
    Real max_abs_eigenvalue = 0.0;      ///< V9: max|eig(伴随)|
    bool is_strictly_stationary = false;  ///< 严格 <1 判定
    Size n_obs_used = 0;
    Size lag = 0;                       ///< 实际滞后阶
    std::string trend = "c";
    Eigen::MatrixXd residuals;          ///< T_eff×K
    Eigen::MatrixXd sigma_u;            ///< df 修正版 Σ (IRF/FEVD 用, v1.2 附)
    std::vector<Eigen::MatrixXd> coeff_vcov;  ///< 逐方程 σ̂²(X'X)⁻¹
};

namespace detail {

/// trend → 确定性列数 (statsmodels get_trendorder 直译)
inline Size trend_order(const std::string& trend) {
    if (trend == "n") return 0;
    if (trend == "c") return 1;
    if (trend == "ct") return 2;
    throw std::invalid_argument("trend must be one of n/c/ct");
}

/// 构造 VAR 设计矩阵与响应 (0-based, 与 statsmodels get_var_endog 对齐):
///   Y = y[p .. T−1] (行序), Z = [det cols, y[p−1], y[p−2], ...]
/// det 列: "c" → 1; "ct" → (1, t) t 从 1 起按行 (statsmodels 时间列约定)
inline void build_var_design(const Eigen::MatrixXd& y, Size p,
                             const std::string& trend,
                             Eigen::MatrixXd& Y, Eigen::MatrixXd& Z) {
    const Size T = y.rows(), K = y.cols();
    const Size kt = trend_order(trend);
    const Size nobs = T - p;
    Y.resize(nobs, K);
    Z.resize(nobs, kt + K * p);
    for (Size i = 0; i < nobs; ++i) Y.row(i) = y.row(p + i);
    Size col = 0;
    if (trend == "c") {
        for (Size i = 0; i < nobs; ++i) Z(i, col) = 1.0;
        ++col;
    } else if (trend == "ct") {
        for (Size i = 0; i < nobs; ++i) {
            Z(i, col) = 1.0;
            Z(i, col + 1) = static_cast<Real>(i + 1);
        }
        col += 2;
    }
    for (Size l = 0; l < p; ++l) {
        for (Size j = 0; j < K; ++j) {
            for (Size i = 0; i < nobs; ++i) Z(i, col) = y(p + i - l - 1, j);
            ++col;
        }
    }
}

/// IC 五式 (V4/V6; statsmodels L2281-2303 逐字)
struct ICSet {
    Real aic, bic, hqic, fpe, logdet;
};

inline ICSet info_criteria(const Eigen::MatrixXd& sigma_mle, Size nobs, Size K,
                           Size p, Size k_trend) {
    // logdet(SPD): LLT 对角 log 和 (statsmodels logdet_symm 同法, 数值稳)
    const Eigen::LLT<Eigen::MatrixXd> llt(sigma_mle);
    if (llt.info() != Eigen::Success) {
        throw std::invalid_argument("info_criteria: sigma not PD");
    }
    const Real ld = 2.0 * llt.matrixL().toDenseMatrix().diagonal().array().log().sum();
    const Real fp = static_cast<Real>(p * K * K + K * k_trend);
    const Real n = static_cast<Real>(nobs);
    const Real df_model = static_cast<Real>(p * K + k_trend);
    const Real df_resid = n - df_model;
    ICSet ic;
    ic.logdet = ld;
    ic.aic = ld + 2.0 / n * fp;
    ic.bic = ld + std::log(n) / n * fp;
    ic.hqic = ld + 2.0 * std::log(std::log(n)) / n * fp;
    ic.fpe = (df_resid > 0)
                 ? std::pow((n + df_model) / df_resid, static_cast<Real>(K)) *
                       std::exp(ld)
                 : std::numeric_limits<Real>::infinity();
    return ic;
}

}  // namespace detail

/// VAR(p) 逐方程 OLS (决策 9). spec.lag = 0 时按 aic 选阶 (same_sample_ic).
inline VARResult var_fit(const MultivariateTSData& data, const VARSpec& spec = {}) {
    data.validate();
    if (spec.trend != "n" && spec.trend != "c" && spec.trend != "ct") {
        throw std::invalid_argument("var_fit: trend must be n/c/ct");
    }
    const Eigen::MatrixXd y = data.matrix();
    const Size T = y.rows(), K = y.cols();
    const Size kt = detail::trend_order(spec.trend);

    Size p = spec.lag;
    if (p == 0) {
        // IC 选阶: 借 var_select 的同样本机制 (V5), aic 准则
        // (避免循环依赖: 此处内联简版 — offset 同样本 aic 扫描)
        Size maxlag = spec.max_lag;
        if (maxlag == 0) {
            maxlag = static_cast<Size>(
                std::llround(12.0 * std::pow(static_cast<Real>(T) / 100.0, 0.25)));
        }
        const Size max_est = (T - K - kt) / (1 + K);
        if (maxlag > max_est) maxlag = max_est;
        if (maxlag < 1) throw std::invalid_argument("var_fit: no estimable lag");
        Real best = std::numeric_limits<Real>::infinity();
        Size best_p = 1;
        for (Size q = 1; q <= maxlag; ++q) {
            const Size off = maxlag - q;  // V5 同样本
            Eigen::MatrixXd Y, Z;
            const Eigen::MatrixXd ysub = y.bottomRows(T - off);
            detail::build_var_design(ysub, q, spec.trend, Y, Z);
            const Size nobs = Y.rows();
            const Eigen::MatrixXd B = Z.householderQr().solve(Y);
            const Eigen::MatrixXd R = Y - Z * B;
            Eigen::MatrixXd smle = R.transpose() * R / static_cast<Real>(nobs);
            auto ic = detail::info_criteria(smle, nobs, K, q, kt);
            if (ic.aic < best) {
                best = ic.aic;
                best_p = q;
            }
        }
        p = best_p;
    }

    if (p >= T - kt) throw std::invalid_argument("var_fit: lag too large for T");

    Eigen::MatrixXd Y, Z;
    detail::build_var_design(y, p, spec.trend, Y, Z);
    const Size nobs = Y.rows();
    const Size ncol = Z.cols();

    // 逐方程 OLS — 同回归元 ⇒ 矩阵形式一次求解等价 (决策 9)
    const Eigen::MatrixXd B = Z.householderQr().solve(Y);  // ncol×K
    const Eigen::MatrixXd R = Y - Z * B;                    // nobs×K 残差

    // Σ 分层 (V4)
    const Real df_resid_r = static_cast<Real>(nobs) - static_cast<Real>(K * p + kt);
    if (df_resid_r <= 0) throw std::invalid_argument("var_fit: df_resid <= 0");
    Eigen::MatrixXd SSR = R.transpose() * R;
    Eigen::MatrixXd sigma_mle = SSR / static_cast<Real>(nobs);
    Eigen::MatrixXd sigma_df = SSR / df_resid_r;

    auto ic = detail::info_criteria(sigma_mle, nobs, K, p, kt);

    // loglik (var_loglike 直译)
    const Real Kd = static_cast<Real>(K), nd = static_cast<Real>(nobs);
    const Real ll = -(nd * Kd / 2.0) * std::log(2.0 * std::acos(-1.0)) -
                    (nd / 2.0) * (ic.logdet + Kd);

    // 伴随矩阵特征值 (V9): A_comp = [[A1 A2 ... Ap],[I 0 ... 0],...]
    Eigen::MatrixXd comp = Eigen::MatrixXd::Zero(K * p, K * p);
    for (Size l = 0; l < p; ++l) {
        for (Size i = 0; i < K; ++i) {
            for (Size j = 0; j < K; ++j) {
                // B 行=回归元列, 列=方程; A_l[i,j] = 方程 i 系数 on y_j.l+1
                comp(i, l * K + j) = B(kt + l * K + j, i);
            }
        }
    }
    for (Size l = 1; l < p; ++l) {
        for (Size j = 0; j < K; ++j) comp(l * K + j, (l - 1) * K + j) = 1.0;
    }
    Real max_eig = 0.0;
    if (p > 0) {
        const Eigen::VectorXcd ev = comp.eigenvalues();
        for (Size i = 0; i < static_cast<Size>(ev.size()); ++i) {
            max_eig = std::max(max_eig, std::abs(ev[i]));
        }
    }

    // 逐方程系数协方差: σ̂²_i·(Z'Z)^{-1} (df 修正 σ, statsmodels cov_params)
    std::vector<Eigen::MatrixXd> vcov(K);
    {
        Eigen::MatrixXd ZtZ_inv =
            (Z.transpose() * Z).lu().inverse();
        for (Size i = 0; i < K; ++i) {
            vcov[i] = sigma_df(i, i) * ZtZ_inv;
        }
    }

    VARResult res;
    res.coefficients = B.transpose();  // K×(K·p+k_trend), 行=方程
    res.sigma_u_mle = sigma_mle;
    res.sigma_u = sigma_df;
    res.loglik = ll;
    res.aic = ic.aic;
    res.bic = ic.bic;
    res.hqic = ic.hqic;
    res.fpe = ic.fpe;
    res.logdet = ic.logdet;
    res.max_abs_eigenvalue = max_eig;
    res.is_strictly_stationary = (max_eig < 1.0);
    res.n_obs_used = nobs;
    res.lag = p;
    res.trend = spec.trend;
    res.residuals = R;
    res.coeff_vcov = std::move(vcov);
    return res;
}

}  // namespace var
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
