// =============================================================================
// garch_m_model.hpp - GARCH(1,1)-M 三变体 QMLE 估计 (spec §2.3)
//
// Phase 7C v1.7 M0 (PHASE7C_SPEC.md v1.2)
//
// 教材锚点: Engle-Lilien-Robins 1987 Econometrica 55(2):391-407
// 对照库: arch 8.0.0 ARCHInMean (主锚, form='vol'/'var'/'log', rescale=False,
//         1e-8~1e-10) + rugarch archm archpow=1|2 (solver 敏感, 1e-4, verify_gm.R)
//
// 幻觉点防护 (spec §9.5):
//   GM1: form 映射 — Variance g=h (arch 'var', rugarch archpow=2) /
//        Volatility g=√h (arch 'vol', archpow=1) / LogVariance g=log h
//        (arch 'log'; rugarch 无 log — 仅 arch 可锚)
//   GM2: ARCHInMean 为 arch 8.0 新增 (8.0.0 release notes 漏记, 以 tag 源码为准)
//   GM3: log 变体有 arch 数值基准 (升入 scope 依据)
//   GM4: λ 的 SE = QMLE sandwich V[1][1] (Bollerslev-Wooldridge, 非 OLS t 表)
//   GM5: 均值-方差耦合递归 (arch issue #269) — ε_{t-1} 含 λ·g(h_{t-1}),
//        λ 变化重写整条 h 路径; h_t 与 ε_t 逐点互依, 不可两步解耦
//
// arch 对齐语义 (本地 recursions_python.py L1094-1127 + mean.py L1590 一手):
//   - 递归: 每 t 先 update h_t (用 ε_{t-1}, h_{t-1}), 再 ε_t = y_t − μ − λ·g(h_t)
//   - h₁ = ω + (α+β)·σ²₀ (GARCHUpdater t=0 双 backcast, 与 7B G1 同源)
//   - backcast 用**无 in-mean 项**残差 (ARCHInMean.resids 丢末参 → 常数均值残差),
//     即 EWMA(y − ȳ)², 优化期冻结
//   - κ (λ) 起始值 0, 无约束 (spec: 原尺度边界 [−10, 10])
//   - form 幂次: 'log'→0 (专用 log 分支) / 'vol'→σ / 'var'→σ²
//
// 内部尺度变换 (7B 惯例延续): s = 1/sd(y−ȳ), y' = s·y, h' = s²·h, ε' = s·ε。
//   均值方程等变参数化 (逐 form, 由 s·[μ + λ·g(h)] = μ' + λ'·g(h') 解出):
//     vol:  μ' = s·μ,           λ' = λ      (g(h') = s·g(h))
//     var:  μ' = s·μ,           λ' = λ/s    (g(h') = s²·g(h))
//     log:  μ' = s·μ − 2sλ·lns, λ' = s·λ    (g(h') = g(h) + 2lns)
//   逆映射 + Jacobian (log 形 μ' 行含 λ' 交叉项 2lns/s) 见 sandwich 段。
// =============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/rng.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/calibration/optimizer.hpp"  // SLSQP (ADR-018)
#include "cpphub/timeseries/garch/garch_distribution.hpp"
#include "cpphub/timeseries/garch/garch_model.hpp"  // GarchConfig/backcast/detail

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace garch {

/// GARCH-M 均值方程风险项形式 (GM1)
enum class GarchMForm {
    Variance,      ///< g(h) = h    (arch form='var';  rugarch archpow=2)
    Volatility,    ///< g(h) = √h   (arch form='vol';  rugarch archpow=1)
    LogVariance    ///< g(h) = log h (arch form='log'; rugarch 无此形)
};

/// GARCH-M 参数 (θ = μ, λ, ω, α, β)
struct GarchMParams {
    Real mu = 0.0;       ///< 均值方程截距
    Real lambda = 0.0;   ///< 风险溢价系数 (g(h_t) 的系数)
    Real omega = 0.0;    ///< 方差方程常数 (ω > 0, 同 7B G4)
    Real alpha = 0.0;    ///< ARCH 系数 (α ≥ 0)
    Real beta = 0.0;     ///< GARCH 系数 (β ≥ 0)
};

/// GARCH-M 估计结果 (spec §2.3 v1.2 签名)
struct GarchMResult {
    GarchMParams params = {};                 ///< 含 λ; g(·) 形式由 form 决定 (GM1)
    std::vector<Real> conditional_variances;  ///< h_t (原尺度)
    std::vector<Real> residuals;              ///< ε_t = r_t − μ − λ·g(h_t) (GM5 耦合)
    std::vector<Real> std_residuals;          ///< z_t = ε_t/√h_t (G11 同源)
    Real log_likelihood = 0.0;                ///< 完整似然 (G3 同源, 含 −½log(2π))
    Real aic = 0.0;
    Real bic = 0.0;
    std::vector<std::vector<Real>> vcov;      ///< (5[+ν])² sandwich (μ,λ,ω,α,β) (GM4)
    std::vector<Real> std_errors;             ///< std_errors[1] = λ 的 robust SE (GM4)
    bool converged = false;
    Size n_iterations = 0;
    std::string message;
};

namespace detail {

/// g(h) 形式函数 (GM1; 公开以便测试)
inline Real gm_g(Real h, GarchMForm form) {
    switch (form) {
        case GarchMForm::Variance:
            return h;
        case GarchMForm::Volatility:
            return std::sqrt(h);
        case GarchMForm::LogVariance:
            return std::log(h);
    }
    throw std::invalid_argument("gm_g: unknown form");
}

/// GM5 耦合路径 (公开以便测试): h 与 ε 逐点互依
struct GmPath {
    std::vector<Real> h;    ///< 条件方差路径
    std::vector<Real> eps;  ///< 残差路径 (含 λ·g(h) 项)
};

/// GM5 耦合递归 (arch ARCHInMeanRecursion.recursion 逐行对齐):
///   h₁ = ω + (α+β)·σ²₀ (G1 双 backcast);  ε_t = y_t − μ − λ·g(h_t);
///   h_t = ω + α·ε²_{t−1} + β·h_{t−1}  ← ε_{t−1} 含 λ, λ 变化重写全路径
/// @throws std::runtime_error h ≤ 0 / 非有限 (优化层捕获 → 惩罚)
inline GmPath filter_garch_m(GarchMForm form, Real mu, Real lambda, Real omega,
                             Real alpha, Real beta, const std::vector<Real>& y,
                             Real sigma2_0) {
    const Size T = y.size();
    GmPath p;
    p.h.resize(T);
    p.eps.resize(T);
    for (Size t = 0; t < T; ++t) {
        const Real e2_prev = (t == 0) ? sigma2_0 : p.eps[t - 1] * p.eps[t - 1];
        const Real h_prev = (t == 0) ? sigma2_0 : p.h[t - 1];
        const Real ht = omega + alpha * e2_prev + beta * h_prev;
        if (!(ht > 0.0) || !std::isfinite(ht)) {
            throw std::runtime_error("filter_garch_m: invalid variance");
        }
        p.h[t] = ht;
        p.eps[t] = y[t] - mu - lambda * gm_g(ht, form);
        if (!std::isfinite(p.eps[t])) {
            throw std::runtime_error("filter_garch_m: non-finite residual");
        }
    }
    return p;
}

/// λ 原尺度 → 缩放尺度 (等变参数化; 公开以便测试)
inline Real gm_lambda_to_scaled(Real lambda, Real s, GarchMForm form) {
    switch (form) {
        case GarchMForm::Variance:
            return lambda / s;
        case GarchMForm::Volatility:
            return lambda;
        case GarchMForm::LogVariance:
            return lambda * s;
    }
    throw std::invalid_argument("gm_lambda_to_scaled: unknown form");
}

/// λ 缩放尺度 → 原尺度
inline Real gm_lambda_from_scaled(Real lambda_s, Real s, GarchMForm form) {
    switch (form) {
        case GarchMForm::Variance:
            return lambda_s * s;
        case GarchMForm::Volatility:
            return lambda_s;
        case GarchMForm::LogVariance:
            return lambda_s / s;
    }
    throw std::invalid_argument("gm_lambda_from_scaled: unknown form");
}

/// μ 原尺度 → 缩放尺度
inline Real gm_mu_to_scaled(Real mu, Real lambda, Real s, GarchMForm form) {
    const Real base = mu * s;
    if (form == GarchMForm::LogVariance) {
        return base - 2.0 * lambda * s * std::log(s);
    }
    return base;
}

/// μ 缩放尺度 → 原尺度
inline Real gm_mu_from_scaled(Real mu_s, Real lambda_s, Real s, GarchMForm form) {
    const Real base = mu_s / s;
    if (form == GarchMForm::LogVariance) {
        // μ' = s·μ − 2λ'·lns ⇒ μ = (μ' + 2λ'·lns)/s
        return base + 2.0 * lambda_s * std::log(s) / s;
    }
    return base;
}

/// ∂λ/∂λ' (Jacobian 用, 逐 form)
inline Real gm_lambda_jacobian(Real s, GarchMForm form) {
    switch (form) {
        case GarchMForm::Variance:
            return s;
        case GarchMForm::Volatility:
            return 1.0;
        case GarchMForm::LogVariance:
            return 1.0 / s;
    }
    throw std::invalid_argument("gm_lambda_jacobian: unknown form");
}

}  // namespace detail

/// @brief GARCH(1,1)-M QMLE 估计 (三变体)
///
/// 模型: r_t = μ + λ·g(h_t) + ε_t;  h_t = ω + α·ε²_{t−1} + β·h_{t−1} (GM5 耦合)
///
/// @param data 收益率序列 (水平值, 无 NaN, T ≥ 10)
/// @param form g(·) 形式 (GM1; 默认 Volatility = arch 默认 'vol')
/// @param config 复用 GarchConfig; initial_params = {μ,λ,ω,α,β[,ν]} (原尺度, 可选)
/// @return GarchMResult (vcov = (μ,λ,ω,α,β[,ν]) 全参数 QMLE sandwich, GM4)
/// @throws std::invalid_argument T<10 / NaN / 零方差
inline GarchMResult estimate_garch_m(const std::vector<Real>& data,
                                     GarchMForm form = GarchMForm::Volatility,
                                     const GarchConfig& config = GarchConfig{}) {
    const Size T = data.size();
    if (T < 10) {
        throw std::invalid_argument(
            "estimate_garch_m: need at least 10 observations");
    }
    Real mean = 0.0;
    for (Real v : data) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("estimate_garch_m: data contains NaN/Inf");
        }
        mean += v;
    }
    mean /= static_cast<Real>(T);

    Real var = 0.0;
    for (Size t = 0; t < T; ++t) {
        const Real e = data[t] - mean;
        var += e * e;
    }
    var /= static_cast<Real>(T);
    Real scale_abs = 0.0;
    for (Real v : data) scale_abs = std::max(scale_abs, std::abs(v));
    const Real rel_sd = (scale_abs > 0.0) ? std::sqrt(var) / scale_abs : 0.0;
    if (!(var > 0.0) || rel_sd < 1e-13) {
        throw std::invalid_argument("estimate_garch_m: zero variance data");
    }

    // 内部尺度变换: y' = s·y (var(y') = 1)
    const Real s = 1.0 / std::sqrt(var);
    const Real log_s = std::log(s);
    std::vector<Real> y_s(T);
    for (Size t = 0; t < T; ++t) y_s[t] = data[t] * s;
    const Real mean_s = mean * s;

    const bool has_nu = config.dist != GarchDist::Normal;
    const Size k = has_nu ? 6u : 5u;

    // backcast: 无 in-mean 项残差 (arch ARCHInMean.resids 丢末参 → EWMA(y−ȳ)²)
    std::vector<Real> eps0(T);
    for (Size t = 0; t < T; ++t) eps0[t] = y_s[t] - mean_s;
    const Real sigma2_0 = backcast_variance(eps0, config.backcast_lambda,
                                            config.backcast_window);

    // 逐观测似然 (θ' = μ',λ',ω',α,β[,ν]; GM5 耦合递归)
    const detail::PerObsLL pll = [&](const std::vector<Real>& theta) {
        const Real nu = has_nu ? theta[5] : 0.0;
        const detail::GmPath p = detail::filter_garch_m(
            form, theta[0], theta[1], theta[2], theta[3], theta[4], y_s, sigma2_0);
        std::vector<Real> ll(T);
        for (Size t = 0; t < T; ++t) {
            ll[t] = log_likelihood_term(p.eps[t], p.h[t], config.dist, nu);
        }
        return ll;
    };

    // 边界: λ' 等变换算 (原尺度 ±10, spec §2.3 Step3); ω'/α/β 同 7B G4
    const Real lam_hi = std::abs(detail::gm_lambda_to_scaled(10.0, s, form));
    std::vector<Bounds> bnds = {{mean_s - 10.0, mean_s + 10.0},
                                {-lam_hi, lam_hi},
                                {1e-8, 100.0},
                                {0.0, 1.0},
                                {0.0, 1.0}};
    if (has_nu) bnds.push_back({2.05, 500.0});
    const ConstraintFn ineq = [](const std::vector<Real>& x) {
        std::vector<Real> c = {x[2] - 1e-8, x[3], x[4],
                               1.0 - x[3] - x[4] - 1e-8};
        return c;
    };

    // 多起始 (scaled): κ 起始 0 (arch); 随机起点 λ' 加性扰动逃逸 0 陷阱
    const Real nu0 = 8.0;
    std::vector<std::vector<Real>> starts;
    auto push5 = [&](Real mu, Real lam, Real om, Real al, Real be) {
        std::vector<Real> st = {mu, lam, om, al, be};
        if (has_nu) st.push_back(nu0);
        starts.push_back(std::move(st));
    };
    push5(mean_s, 0.0, 0.10, 0.10, 0.85);  // HR
    push5(mean_s, 0.0, 0.05, 0.05, 0.90);  // ML
    if (config.initial_params.size() >= 5) {
        const auto& ip = config.initial_params;
        std::vector<Real> st = {
            detail::gm_mu_to_scaled(ip[0], ip[1], s, form),
            detail::gm_lambda_to_scaled(ip[1], s, form), ip[2] * s * s, ip[3],
            ip[4]};
        if (has_nu && ip.size() >= 6) {
            st.push_back(ip[5]);
        } else if (has_nu) {
            st.push_back(nu0);
        }
        starts.push_back(std::move(st));
    }
    {  // 随机: λ' 加性 ~N(0, 0.2²) 截断入界, ω' 乘性 10% (Philox 确定性)
        Philox4x64 rng(42, 11);
        const uint64_t r1 = rng();
        const uint64_t r2 = rng();
        const Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
        const Real u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
        const auto [z1, z2] = box_muller(u1, u2);
        Real lam_r = 0.2 * z1;
        lam_r = std::clamp(lam_r, -lam_hi, lam_hi);
        push5(mean_s, lam_r, 0.10 * (1.0 + 0.10 * z2), 0.10, 0.85);
    }
    if (!config.use_multistart) starts.resize(1);

    // SLSQP 多起始
    const detail::MultistartResult fit =
        detail::run_multistart(pll, starts, bnds, ineq, config.optimizer_config);

    // 最优参数处重组 (scaled)
    const Real mu_s = fit.x[0], lam_s = fit.x[1], om_s = fit.x[2], al_s = fit.x[3],
               be_s = fit.x[4];
    const Real nu_s = has_nu ? fit.x[5] : 0.0;
    const detail::GmPath path = detail::filter_garch_m(
        form, mu_s, lam_s, om_s, al_s, be_s, y_s, sigma2_0);
    Real ll_scaled = 0.0;
    for (Size t = 0; t < T; ++t) {
        ll_scaled +=
            log_likelihood_term(path.eps[t], path.h[t], config.dist, nu_s);
    }

    // 逆映射回原尺度 (等变参数化)
    GarchMResult result{};
    result.params.mu = detail::gm_mu_from_scaled(mu_s, lam_s, s, form);
    result.params.lambda = detail::gm_lambda_from_scaled(lam_s, s, form);
    result.params.omega = om_s / (s * s);
    result.params.alpha = al_s;
    result.params.beta = be_s;
    result.conditional_variances.resize(T);
    result.residuals.resize(T);
    result.std_residuals.resize(T);
    for (Size t = 0; t < T; ++t) {
        result.conditional_variances[t] = path.h[t] / (s * s);
        result.residuals[t] = path.eps[t] / s;
        result.std_residuals[t] =
            result.residuals[t] / std::sqrt(result.conditional_variances[t]);
    }
    result.log_likelihood = ll_scaled + static_cast<Real>(T) * log_s;
    result.aic = -2.0 * result.log_likelihood + 2.0 * static_cast<Real>(k);
    result.bic = -2.0 * result.log_likelihood
                 + static_cast<Real>(k) * std::log(static_cast<Real>(T));

    // GM4: sandwich (scaled 空间) + form 感知 Jacobian 映射
    // J = ∂(μ,λ,ω,α,β[,ν])/∂(μ',λ',ω',α,β[,ν]):
    //   vol: diag(1/s, 1, 1/s²);  var: diag(1/s, s, 1/s²)
    //   log: μ 行 = [1/s, 2·lns/s, 0,...], λ 行 = [0, 1/s, 0,...], ω 行 1/s²
    if (config.compute_sandwich) {
        try {
            const auto H = detail::numerical_hessian(pll, fit.x);
            const Size npar = fit.x.size();
            const auto G = detail::opg_scores(pll, fit.x);
            const Size Tn = G.empty() ? 0 : G[0].size();
            std::vector<Real> gbar(npar, 0.0);
            for (Size a = 0; a < npar; ++a) {
                Real acc = 0.0;
                for (Size t = 0; t < Tn; ++t) acc += G[a][t];
                gbar[a] = acc / static_cast<Real>(Tn);
            }
            std::vector<std::vector<Real>> S(npar, std::vector<Real>(npar, 0.0));
            const Real sc = static_cast<Real>(Tn) / static_cast<Real>(Tn - 1);
            for (Size t = 0; t < Tn; ++t) {
                for (Size a = 0; a < npar; ++a)
                    for (Size b = 0; b < npar; ++b) {
                        S[a][b] += (G[a][t] - gbar[a]) * (G[b][t] - gbar[b]);
                    }
            }
            for (Size a = 0; a < npar; ++a)
                for (Size b = 0; b < npar; ++b) S[a][b] *= sc;
            std::vector<std::vector<Real>> Hinv;
            if (!detail::invert_matrix(H, Hinv)) {
                result.message = "sandwich: Hessian inversion failed";
            } else {
                std::vector<std::vector<Real>> HS(npar,
                                                  std::vector<Real>(npar, 0.0));
                for (Size a = 0; a < npar; ++a)
                    for (Size b = 0; b < npar; ++b) {
                        Real acc = 0.0;
                        for (Size m = 0; m < npar; ++m) acc += Hinv[a][m] * S[m][b];
                        HS[a][b] = acc;
                    }
                std::vector<std::vector<Real>> V(npar, std::vector<Real>(npar, 0.0));
                for (Size a = 0; a < npar; ++a)
                    for (Size b = 0; b < npar; ++b) {
                        Real acc = 0.0;
                        for (Size m = 0; m < npar; ++m) acc += HS[a][m] * Hinv[m][b];
                        V[a][b] = acc;
                    }
                // Jacobian (scaled → original)
                std::vector<std::vector<Real>> J(npar,
                                                 std::vector<Real>(npar, 0.0));
                for (Size a = 0; a < npar; ++a) J[a][a] = 1.0;
                J[0][0] = 1.0 / s;                              // μ
                J[1][1] = detail::gm_lambda_jacobian(s, form);   // λ
                J[2][2] = 1.0 / (s * s);                        // ω
                if (form == GarchMForm::LogVariance) {
                    J[0][1] = 2.0 * log_s / s;                  // μ ← λ' 交叉
                }
                // V_orig = J·V·J'
                result.vcov.assign(npar, std::vector<Real>(npar, 0.0));
                for (Size a = 0; a < npar; ++a)
                    for (Size b = 0; b < npar; ++b) {
                        Real acc = 0.0;
                        for (Size i = 0; i < npar; ++i)
                            for (Size j = 0; j < npar; ++j) {
                                acc += J[a][i] * V[i][j] * J[b][j];
                            }
                        result.vcov[a][b] = acc;
                    }
                result.std_errors.resize(npar);
                for (Size a = 0; a < npar; ++a) {
                    result.std_errors[a] =
                        (result.vcov[a][a] > 0.0) ? std::sqrt(result.vcov[a][a])
                                                  : 0.0;
                }
            }
        } catch (const std::exception& e) {
            result.message = std::string("sandwich failed: ") + e.what();
        }
    }

    result.converged = fit.converged;
    result.n_iterations = fit.n_iterations;
    if (result.message.empty()) result.message = fit.message;
    return result;
}

}  // namespace garch
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
