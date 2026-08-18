// =============================================================================
// hannan_rissanen.hpp - Hannan-Rissanen 初值估计 (spec §3.1 Step4 / AR7)
//
// Phase 7C v1.7 M1 (PHASE7C_SPEC.md v1.2 §3.1)
//
// 教材锚点: Hannan-Rissanen 1982; statsmodels
//   tsa/arima/estimators/hannan_rissanen.py (Step1/Step2 语义一手, 2026-08-18)
//
// 定位: 仅作多起始点的初值之一 (use_hannan_rissanen = true 时),
//   不承诺与 statsmodels HR 估计器逐位对照 (起始值不同 → 同一最优落点
//   差在优化器容差层; AR6 配对对照仅锁最终估计, 见 verify_arima.*)。
//   Step 3 bias correction 不实现 (起始值用途无必要, spec §10.5 边界)。
//
// 算法 (statsmodels hannan_rissanen L159-215 简化, 无 fixed params):
//   Step 1: 长 AR 阶 n_init = Schwert(T) (复用 unit_root::schwert_lag;
//           statsmodels 默认 max(10·log10(T), max_ar+max_ma+1), 起始值
//           用途下允许口径差) — Yule-Walker(mle) 系数 π, 残差 û_t
//   Step 2: OLS: z_t ~ z_{t−1..p} + û_{t−1..q}
//           样本对齐 (statsmodels L172/181): t 从 n_init + q 起 (0-based
//           z 索引), y = z[n_init+q:], AR 列 z[t−i], MA 列 û[t−j]
//   σ̂² = SSR/(n − k)
//
// 参数化: φ = AR 系数 (1 − φ₁B − …); θ = MA (1 + θ₁B + …) 正号 (AR1 同源)
// =============================================================================

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace arima {

/// HR 估计结果 (仅初值用途; 无 vcov/诊断)
struct HannanRissanenResult {
    std::vector<Real> phi;    ///< AR 初值 (1 − φ₁B − …)
    std::vector<Real> theta;  ///< MA 初值 (1 + θ₁B + …, 正号 AR1)
    Real sigma2 = 0.0;        ///< SSR/(n−k)
    std::vector<Real> resid;  ///< Step 2 OLS 残差
    Size n_init_ar = 0;       ///< Step 1 长 AR 阶 (回显)
};

namespace detail {

/// Yule-Walker AR(k) 系数 (Levinson-Durbin; statsmodels yule_walker
/// method='mle' 语义: acovf 用除以 T 的有偏自协方差)
/// 返回 {π₁..π_k}, γ 由序列直接计算 (FFT 不需要, T 适中)
inline std::vector<Real> yule_walker_coeffs(const std::vector<Real>& z, Size k) {
    const Size T = z.size();
    if (T <= k + 1) {
        throw std::invalid_argument("yule_walker: sample too small");
    }
    // 有偏自协方差 γ_j = (1/T)Σ_{t=0}^{T−j−1} z_t z_{t+j} (mle 口径)
    Real mean = 0.0;
    for (Real v : z) mean += v;
    mean /= static_cast<Real>(T);
    std::vector<Real> d(T);
    for (Size t = 0; t < T; ++t) d[t] = z[t] - mean;
    std::vector<Real> gamma(k + 1, 0.0);
    for (Size j = 0; j <= k; ++j) {
        Real s = 0.0;
        for (Size t = 0; t + j < T; ++t) s += d[t] * d[t + j];
        gamma[j] = s / static_cast<Real>(T);
    }
    // Levinson-Durbin 递归 (φ[0..m−1] 为 m 阶系数)
    std::vector<Real> phi(k, 0.0);
    Real sigma2 = gamma[0];
    if (sigma2 <= 0.0) {
        throw std::invalid_argument("yule_walker: zero variance");
    }
    for (Size m = 1; m <= k; ++m) {
        Real num = gamma[m];
        for (Size i = 1; i < m; ++i) {
            num -= phi[i - 1] * gamma[m - i];
        }
        const Real delta = num / sigma2;
        // φ_i ← φ_i − δ·φ_{m−i} (i = 1..m−1), 再置 φ_m = δ
        std::vector<Real> old_phi(phi.begin(), phi.begin() + (m - 1));
        for (Size i = 1; i <= m - 1; ++i) {
            phi[i - 1] = old_phi[i - 1] - delta * old_phi[m - 1 - i];
        }
        phi[m - 1] = delta;
        sigma2 = sigma2 * (1.0 - delta * delta);
        if (sigma2 <= 1e-300) {
            throw std::invalid_argument("yule_walker: degenerate (perfect fit)");
        }
    }
    return phi;
}

}  // namespace detail

/// @brief Hannan-Rissanen 两步初值估计
///
/// @param z 已差分/去均值序列 (demean 由调用方处理, statsmodels HR
///          innovations_mle 调用链 demean=False 因上游已减均值)
/// @param p AR 阶; @param q MA 阶
/// @param init_ar_order 0 => Schwert(T); >0 => 用户指定 Step1 长 AR 阶
/// @return {phi, theta, sigma2, resid}
/// @throws std::invalid_argument T 过小 / 零方差
inline HannanRissanenResult hannan_rissanen(const std::vector<Real>& z,
                                            Size p, Size q,
                                            Size init_ar_order = 0) {
    const Size T = z.size();
    const Size n_init = (init_ar_order > 0)
                              ? init_ar_order
                              : unit_root::schwert_lag(T);
    if (T < n_init + q + p + 3) {
        throw std::invalid_argument(
            "hannan_rissanen: sample too small for (p,q,init)");
    }

    // Step 1: 长 AR 残差 û (statsmodels: yule_walker(mle) → resid)
    std::vector<Real> pi;
    try {
        pi = detail::yule_walker_coeffs(z, n_init);
    } catch (const std::invalid_argument&) {
        // 退化序列: 零系数兜底 (起始值用途, 非终止路径)
        pi.assign(n_init, 0.0);
    }
    std::vector<Real> u(T, 0.0);  // u[t] 对 t < n_init 无定义, 置 0
    for (Size t = n_init; t < T; ++t) {
        Real m = 0.0;
        for (Size i = 0; i < n_init; ++i) {
            m += pi[i] * z[t - 1 - i];
        }
        u[t] = z[t] - m;
    }

    HannanRissanenResult res;
    res.n_init_ar = n_init;
    if (p == 0 && q == 0) {
        // statsmodels L127-129: 纯方差
        Real mean = 0.0;
        for (Real v : z) mean += v;
        mean /= static_cast<Real>(T);
        Real s2 = 0.0;
        for (Real v : z) s2 += (v - mean) * (v - mean);
        res.sigma2 = s2 / static_cast<Real>(T);
        return res;
    }

    // Step 2: OLS z_t ~ z_{t−1..p} + u_{t−1..q}, t = n_init+q .. T−1 (0-based)
    const Size t0 = n_init + q;
    const Size n = T - t0;
    const Size k = p + q;
    if (k == 0) {  // p=0,q=0 已在上方处理; 此处防御
        throw std::invalid_argument("hannan_rissanen: p=q=0 handled above");
    }
    std::vector<Real> y(n);
    std::vector<std::vector<Real>> X(n, std::vector<Real>(k, 0.0));
    for (Size i = 0; i < n; ++i) {
        const Size t = t0 + i;
        y[i] = z[t];
        for (Size a = 0; a < p; ++a) {
            X[i][a] = z[t - 1 - a];
        }
        for (Size b = 0; b < q; ++b) {
            X[i][p + b] = u[t - 1 - b];
        }
    }
    const auto fit = unit_root::detail::ols_fit(y, X);
    res.phi.assign(fit.beta.begin(), fit.beta.begin() + p);
    res.theta.assign(fit.beta.begin() + p, fit.beta.end());
    res.sigma2 = fit.ssr / static_cast<Real>(n - k);
    res.resid = fit.resid;
    return res;
}

}  // namespace arima
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
