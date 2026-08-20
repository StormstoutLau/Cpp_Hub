// =============================================================================
// arima_model.hpp - ARIMA(p,d,q) 估计: CSS / CSS-ML / Innovations (spec §3.1)
//
// Phase 7C v1.7 M1 (PHASE7C_SPEC.md v1.2 §3.1)
//
// 教材锚点: Box-Jenkins 5th ed Ch.7 (CSS) / B&D 2016 §5.2 (innovations,
//   innovations_mle.hpp) / forecast::Arima (drift 语义, AR5)
//
// 幻觉点防护 (spec §9.1):
//   AR1: MA 正号 (1+θB), 与 R/statsmodels 同号逐系数对照
//   AR2: n_cond = d + max(user_n.cond=0, p) — 与 q 无关 (spec 审计裁决,
//        arima.R L158-162; v1.0 的 max(p,q)+1 已修正); CSS 起始 ε_{t<t0}=0
//   AR3: loglik 口径 — CSS 报统一高斯型 −0.5·n·[log(2πσ̂²)+1]
//        (R 纯 CSS 只报 part log-likelihood, 不与 ML 混比; 文档声明);
//        CSS-ML/Innovations 报精确似然 (innovations 算法)
//   AR4: AIC 基于 T−d 个观测, k = p+q+1+(drift?1:0) (含 σ², R 口径)
//   AR5: drift = 差分截距 (仅 d≥1 有意义; d=0 无均值项 — 对齐 R
//        include.mean=FALSE 语义, 调用方需均值模型时自行预去均值并
//        用 intercept = μ/(1−Σφ) 换算, spec 决策 3)
//   AR6: method 配对对照 — CSS↔R method="CSS", CSS-ML↔R "CSS-ML",
//        Innovations↔statsmodels innovations_mle; 混配禁用
//   AR7: 多起始 {HR, CSS 解, CSS±扰动, 随机重启} 各精化取最优
//   AR8: Innovations 仅无缺失+无季节 (guard 在 innovations_mle)
//
// CSS 递归 (Step 2): z 轴起始 t0 = max(0, p) (AR2 的 n_cond − d);
//   ε_t = z_t − drift − Σφᵢz_{t−i} − Σθⱼε_{t−j}, ε_{t<t0} = 0 不累加;
//   SSR(β) = Σ_{t≥t0} ε_t², SLSQP 无约束 (bounds ±3 防护发散)
// =============================================================================

#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/calibration/optimizer.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/arima/hannan_rissanen.hpp"
#include "cpphub/timeseries/arima/innovations_mle.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace arima {

/// ARIMA 阶设定 (spec §3.1 冻结签名)
struct ArimaSpec {
    Size p = 0;
    Size d = 0;
    Size q = 0;
    bool include_drift = false;  ///< 仅 d≥1 有意义: 差分截距 (AR5)
};

/// ARIMA 参数 (MA 正号参数化, AR1)
struct ArimaParams {
    std::vector<Real> phi;   ///< AR: (1 − φ₁B − …)
    std::vector<Real> theta; ///< MA: (1 + θ₁B + …)
    Real drift = 0.0;
    Real sigma2 = 0.0;
};

enum class ArimaMethod { CSS, CSS_ML, Innovations };

/// ARIMA 估计结果 (spec §3.1 冻结签名 + v1.2 回显字段)
struct ArimaResult {
    ArimaParams params;
    ArimaMethod method{};
    Size n_cond = 0;         ///< CSS 起始条件数 (y 轴: d + max(0,p), AR2)
    Size n_obs_used = 0;     ///< T − d (差分后观测数)
    Real loglik = 0.0;
    Real aic = 0.0;
    std::vector<Real> residuals;              ///< ε_t (CSS) 或 u_t (ML/Innov)
    std::vector<std::vector<Real>> vcov;      ///< QMLE sandwich (可选, 空=未算)
    bool converged = false;
    Size n_iterations = 0;
    std::string message;
};

/// 估计配置 (spec §3.1 v1.2)
struct ArimaConfig {
    Size n_starting_points = 4;
    bool use_hannan_rissanen = true;
    bool compute_sandwich = false;   ///< 默认关; 纯 CSS 无意义
    Size seed = 42;
    SLSQP::Config optimizer_config = [] {
        SLSQP::Config c;
        c.max_iterations = 300;
        c.ftol = 1e-11;
        c.xtol = 1e-10;
        return c;
    }();
};

namespace detail {

/// CSS 递归残差 (AR2 口径): ε_{t<t0} 置 0, 返回全长 (T_z 个; 前 t0 个无效)
inline std::vector<Real> css_residuals(const std::vector<Real>& z,
                                       const std::vector<Real>& phi,
                                       const std::vector<Real>& theta,
                                       Real drift) {
    const Size T = z.size();
    const Size p = phi.size();
    const Size q = theta.size();
    const Size t0 = p;  // z 轴起始 = n_cond − d = max(0,p)
    std::vector<Real> eps(T, 0.0);
    for (Size t = t0; t < T; ++t) {
        Real e = z[t] - drift;
        for (Size i = 0; i < p; ++i) e -= phi[i] * z[t - 1 - i];
        for (Size j = 0; j < q; ++j) e -= theta[j] * eps[t - 1 - j];
        eps[t] = e;
    }
    return eps;
}

/// CSS 目标: SSR (x = {φ, θ, [drift]}); 惩罚防发散
inline Real css_ssr(const std::vector<Real>& x, const std::vector<Real>& z,
                    Size p, Size q, bool has_drift) {
    const std::vector<Real> phi(x.begin(), x.begin() + p);
    const std::vector<Real> theta(x.begin() + p, x.begin() + p + q);
    const Real drift = has_drift ? x[p + q] : 0.0;
    // 发散防护: 参数巨大时残差爆炸 → 检查有限性
    std::vector<Real> eps;
    try {
        eps = css_residuals(z, phi, theta, drift);
    } catch (...) {
        return 1e100;
    }
    Real ssr = 0.0;
    for (Size t = p; t < z.size(); ++t) {
        const Real e = eps[t];
        if (!std::isfinite(e)) return 1e100;
        ssr += e * e;
    }
    return std::isfinite(ssr) ? ssr : 1e100;
}

/// 差分 d 阶 (与 innovations_mle 同语义)
inline std::vector<Real> difference(const std::vector<Real>& y, Size d) {
    std::vector<Real> z = y;
    for (Size i = 0; i < d; ++i) {
        std::vector<Real> dz(z.size() - 1);
        for (Size t = 1; t < z.size(); ++t) dz[t - 1] = z[t] - z[t - 1];
        z = dz;
    }
    return z;
}

}  // namespace detail

/// @brief ARIMA(p,d,q) 估计 (三方法, spec §3.1)
///
/// @param data 原始序列 (无 NaN; Innovations 方法额外要求无缺失)
/// @param spec {p,d,q,include_drift}
/// @param method CSS (条件平方和) / CSS_ML (CSS 起始 + innovations 精确
///        似然精化) / Innovations (直接精确 MLE)
/// @param config 多起始/HR/sandwich/seed/SLSQP
/// @return ArimaResult (AR2/AR3/AR4 口径回显)
/// @throws std::invalid_argument T 过小 / NaN / drift 与 d 组合非法
///   (d=0 && include_drift: R forecast::Arima 语义下无意义 → 拒绝)
inline ArimaResult arima_fit(const std::vector<Real>& data,
                             const ArimaSpec& spec,
                             ArimaMethod method = ArimaMethod::CSS_ML,
                             const ArimaConfig& config = ArimaConfig{}) {
    const Size T = data.size();
    if (spec.d == 0 && spec.include_drift) {
        throw std::invalid_argument(
            "arima_fit: include_drift requires d >= 1 (AR5)");
    }
    for (Real v : data) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("arima_fit: NaN/Inf in data");
        }
    }
    if (T < spec.d + spec.p + spec.q + 5) {
        throw std::invalid_argument("arima_fit: sample too small");
    }

    // Step 1: 差分
    const std::vector<Real> z = detail::difference(data, spec.d);
    const Size Tz = z.size();  // T − d
    ArimaResult out;
    out.n_obs_used = Tz;
    out.n_cond = spec.d + std::max<Size>(0, spec.p);  // AR2 (y 轴)
    out.method = method;
    const bool has_drift = spec.include_drift;
    const Size k = spec.p + spec.q + (has_drift ? 1 : 0);

    if (method == ArimaMethod::Innovations) {
        // AR8: 委托 innovations_mle (guard 在彼)。
        // demean 语义 (2026-08-18 实测裁决): 无 drift 时 demean=true
        //   (= statsmodels innovations_mle 默认, AR6 配对; 实测 arma11
        //   逐位恢复 (0.398, −0.358, ll −416.932)); 有 drift 时 false
        //   (漂移作为显式参数走 CSS-ML 路径联合优化)
        const auto im = innovations_mle(z, spec.p, spec.q, 0, !has_drift,
                                        config.seed);
        out.params.phi = im.phi;
        out.params.theta = im.theta;
        out.params.sigma2 = im.sigma2;
        out.loglik = im.loglik;
        out.residuals = im.innovations;
        out.converged = im.converged;
        out.n_iterations = im.n_iterations;
        out.message = im.message + " [innovations]";
        if (has_drift) {
            // drift 未由 innovations_mle 估计 (其无截距参数):
            // 用 CSS-ML 路径重估 (语义完整), 此处实现为委托
            method = ArimaMethod::CSS_ML;  // fallthrough 到 CSS-ML
            out.message += " -> drift: CSS-ML path";
        } else {
            // AIC (AR4): k = p+q+1
            out.aic = -2.0 * out.loglik
                      + 2.0 * static_cast<Real>(spec.p + spec.q + 1);
            return out;
        }
    }

    // ---- CSS 阶段 (CSS 方法本体, 或 CSS-ML/Innovations+drift 的起始) ----
    const bool direct_css = (method == ArimaMethod::CSS);
    {
        std::vector<std::vector<Real>> starts;
        starts.push_back(std::vector<Real>(k, 0.0));
        if (config.use_hannan_rissanen && spec.p + spec.q > 0) {
            try {
                std::vector<Real> zc = z;
                if (has_drift) {  // HR 对去 drift 序列估计
                    Real m = 0.0;
                    for (Real v : z) m += v;
                    m /= static_cast<Real>(Tz);
                    for (Real& v : zc) v -= m;
                }
                const auto hr = hannan_rissanen(zc, spec.p, spec.q);
                std::vector<Real> x0 = detail::concat_params(hr.phi,
                                                             hr.theta);
                if (has_drift) {
                    Real m = 0.0;
                    for (Real v : z) m += v;
                    x0.push_back(m / static_cast<Real>(Tz));
                }
                starts.push_back(x0);
            } catch (...) {
                // HR 失败 → 仅零起始
            }
        }
        if (config.n_starting_points > starts.size()) {
            Philox4x64 rng(config.seed, 13);
            const Size extra = config.n_starting_points
                               - static_cast<Size>(starts.size());
            for (Size rep = 0; rep < extra; ++rep) {
                std::vector<Real> x0(k, 0.0);
                for (Size i = 0; i < k; ++i) {
                    const uint64_t r1 = rng();
                    const uint64_t r2 = rng();
                    const Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
                    const Real u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
                    const auto [a, b] = box_muller(u1, u2);
                    x0[i] = 0.25 * ((i % 2 == 0) ? a : b);
                }
                starts.push_back(x0);
            }
        }
        std::vector<Bounds> bounds(k, Bounds{-3.0, 3.0});
        if (has_drift) bounds[k - 1] = Bounds{-1e6, 1e6};  // drift 尺度自由
        Real best = std::numeric_limits<Real>::infinity();
        std::vector<Real> xbest;
        for (const auto& x0 : starts) {
            OptimizationResult r = SLSQP::minimize(
                [&](const std::vector<Real>& x) {
                    return detail::css_ssr(x, z, spec.p, spec.q, has_drift);
                },
                x0, bounds, {}, {}, config.optimizer_config);
            if (r.fx < best) {
                best = r.fx;
                xbest = r.x;
                out.converged = r.converged;
                out.n_iterations = r.n_iterations;
                out.message = r.message;
            }
        }
        out.params.phi.assign(xbest.begin(), xbest.begin() + spec.p);
        out.params.theta.assign(xbest.begin() + spec.p,
                                xbest.begin() + spec.p + spec.q);
        if (has_drift) out.params.drift = xbest[spec.p + spec.q];

        if (direct_css) {
            // CSS 输出 (AR3 统一高斯型 loglik)
            const Size n_used = Tz - std::max<Size>(0, spec.p);
            const Real n = static_cast<Real>(n_used);
            out.params.sigma2 = best / n;
            constexpr Real kTwoPi = 6.283185307179586476925286766559;
            out.loglik = -0.5 * n * (std::log(kTwoPi * out.params.sigma2)
                                     + 1.0);
            // AR4: k 参数含 σ²
            out.aic = -2.0 * out.loglik + 2.0 * static_cast<Real>(k + 1);
            out.residuals = detail::css_residuals(
                z, out.params.phi, out.params.theta, out.params.drift);
            out.message += " [CSS]";
            return out;
        }
    }

    // ---- CSS-ML 精化 (innovations 精确似然; CSS 解入起始集, AR7) ----
    {
        const std::vector<Real> z0 = z;
        // 目标: x = {φ, θ, [drift]} — drift 联合优化 (z − drift 后 pieces)
        auto nll_ml = [&](const std::vector<Real>& x) -> Real {
            const std::vector<Real> phi(x.begin(),
                                        x.begin() + spec.p);
            const std::vector<Real> theta(x.begin() + spec.p,
                                          x.begin() + spec.p + spec.q);
            const Real drift = has_drift ? x[spec.p + spec.q] : 0.0;
            std::vector<Real> zd(z0.size());
            for (Size t = 0; t < z0.size(); ++t) zd[t] = z0[t] - drift;
            return detail::concentrated_nll(
                detail::concat_params(phi, theta), zd, spec.p, spec.q);
        };
        // 起始集: {CSS 解, HR + 均值drift, 零, 扰动}
        std::vector<std::vector<Real>> starts;
        {
            std::vector<Real> xcss = detail::concat_params(
                out.params.phi, out.params.theta);
            if (has_drift) xcss.push_back(out.params.drift);
            starts.push_back(xcss);
        }
        if (config.use_hannan_rissanen && spec.p + spec.q > 0) {
            try {
                Real m = 0.0;
                for (Real v : z) m += v;
                m /= static_cast<Real>(Tz);
                std::vector<Real> zc = z;
                for (Real& v : zc) v -= m;
                const auto hr = hannan_rissanen(zc, spec.p, spec.q);
                std::vector<Real> x0 = detail::concat_params(hr.phi,
                                                             hr.theta);
                if (has_drift) x0.push_back(m);
                starts.push_back(x0);
            } catch (...) {
            }
        }
        starts.push_back(std::vector<Real>(k, 0.0));
        {
            Philox4x64 rng(config.seed, 17);
            std::vector<Real> x0(k, 0.0);
            for (Size i = 0; i < k; ++i) {
                const uint64_t r1 = rng();
                const Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
                const auto [a, b] = box_muller(u1, u1);
                x0[i] = 0.15 * a;
            }
            starts.push_back(x0);
        }
        std::vector<Bounds> bounds(k, Bounds{-0.999, 0.999});
        if (has_drift) bounds[k - 1] = Bounds{-1e6, 1e6};
        Real best = std::numeric_limits<Real>::infinity();
        std::vector<Real> xbest;
        bool conv = false;
        Size niter = 0;
        std::string msg;
        SLSQP::Config cfg_ml = config.optimizer_config;
        cfg_ml.max_iterations = 400;
        for (const auto& x0 : starts) {
            OptimizationResult r = SLSQP::minimize(nll_ml, x0, bounds, {},
                                                   {}, cfg_ml);
            if (r.fx < best) {
                best = r.fx;
                xbest = r.x;
                conv = r.converged;
                niter = r.n_iterations;
                msg = r.message;
            }
        }
        out.params.phi.assign(xbest.begin(), xbest.begin() + spec.p);
        out.params.theta.assign(xbest.begin() + spec.p,
                                xbest.begin() + spec.p + spec.q);
        if (has_drift) out.params.drift = xbest[spec.p + spec.q];

        // σ̂²/loglik: 精确 pieces (z − drift)
        const Real drift = out.params.drift;
        std::vector<Real> zd(Tz);
        for (Size t = 0; t < Tz; ++t) zd[t] = z[t] - drift;
        const auto pc = detail::arma_innovations_pieces(zd, out.params.phi,
                                                        out.params.theta);
        const Real Tn = static_cast<Real>(Tz);
        out.params.sigma2 = pc.s_uv / Tn;
        constexpr Real kTwoPi = 6.283185307179586476925286766559;
        out.loglik = -0.5 * (Tn * std::log(kTwoPi * out.params.sigma2)
                             + Tn + pc.sum_logv);
        out.aic = -2.0 * out.loglik + 2.0 * static_cast<Real>(k + 1);
        // 残差 = innovations u_t (z − drift 序列) — 快速路径 (D-4)
        out.residuals =
            detail::arma_innovations_uv(zd, out.params.phi,
                                        out.params.theta).first;
        out.converged = conv;
        out.n_iterations = niter;
        out.message = msg + " [CSS-ML]";
        return out;
    }
}

}  // namespace arima
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
