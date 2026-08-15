// =============================================================================
// gjr_garch_model.hpp - GJR-GARCH(1,1) QMLE 估计 (spec §2.0.4)
//
// Phase 7B v1.6 M1 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 教材锚点: Glosten-Jagannathan-Runkle 1993 / Tsay 3ed Ch 6 /
//           arch arch_model(vol='GARCH', p=1, o=1, q=1)
// 优化器: SLSQP (ADR-018, calibration/optimizer.hpp)
//
// 幻觉点防护 (spec §6.1 + verify_gjr.py probe 实测 2026-08-15):
//   G7:  I(zₜ<0) ≡ I(εₜ<0) (σₜ>0 → 同号), arch 用 I(ε<0)
//   G8:  hₜ = ω + α·ε²ₜ₋₁ + γ·I(εₜ₋₁<0)·ε²ₜ₋₁ + β·hₜ₋₁;
//        γ>0 表杠杆效应 (负冲击放大波动)
//   G10: 平稳性条件 α + γ/2 + β < 1 (非 α+γ+β<1, E[I(z<0)]=1/2);
//        arch 约束矩阵对 o 项系数 -0.5
//   probe-1 (实测): h₁ = ω + (α+γ/2+β)·σ²₀ — arch compute_variance t=0
//        对 o 项用 0.5·backcast (probe_gjr_h1.py: 约定 A maxdiff 4.4e-16,
//        I(ε₀) 约定 3.6e-2)
//   G3/G9/G11/G16/G17: 同 GARCH(1,1) (garch_model.hpp)
//
// 内部尺度变换 (同 GARCH(1,1)): GJR 递归对 ε²/h 线性, MLE 尺度等变性
// 成立。内部 s = 1/sd(ε) 缩放, 估计后 ω = ω'/s² 映射回原尺度,
// V = J·V'·J', J = diag(1/s², 1, 1, 1, [1])。
//
// γ 边界 [-1, 1] (G-gamma-sign: γ 可为负, 反向杠杆; 仅平稳性约束限制上界)
// =============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <algorithm>

#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/calibration/optimizer.hpp"
#include "cpphub/timeseries/garch/garch_model.hpp"      // GarchConfig/detail
#include "cpphub/timeseries/garch/garch_distribution.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace garch {

// GJR-GARCH 参数 (Glosten-Jagannathan-Runkle 1993; 顺序 = arch: ω,α,γ,β)
struct GjrGarchParams {
    Real omega;   ///< 常数项 ω (> 0)
    Real alpha;   ///< ARCH 系数 α (≥ 0)
    Real gamma;   ///< 非对称系数 γ (G8 杠杆效应; ∈ [-1,1], 可为负)
    Real beta;    ///< GARCH 系数 β (≥ 0)
    Real nu;      ///< t/GED 自由度 (Normal 时忽略)
};

// GJR-GARCH 估计结果 (结构同 GarchResult, 参数替换为 GjrGarchParams)
struct GjrGarchResult {
    GjrGarchParams params;
    std::vector<Real> conditional_variances;
    std::vector<Real> residuals;
    std::vector<Real> std_residuals;
    Real log_likelihood;
    Real aic;
    Real bic;
    std::vector<std::vector<Real>> vcov;
    std::vector<Real> std_errors;
    bool converged;
    Size n_iterations;
    std::string message;
};

/// @brief GJR-GARCH(1,1) 方差递归 (公开以便测试, 与 arch compute_variance 一致)
///
/// h₁ = ω + (α + γ/2 + β)·σ²₀  (probe-1: arch t=0 对 o 项用 0.5·backcast,
///                              即 ε²₀ 与 h₀ 均取 backcast 且指示项取 1/2)
/// hₜ = ω + α·ε²ₜ₋₁ + γ·I(εₜ₋₁<0)·ε²ₜ₋₁ + β·hₜ₋₁  (t ≥ 2, G7/G8)
inline std::vector<Real> filter_gjr(const GjrGarchParams& params,
                                    const std::vector<Real>& residuals,
                                    Real sigma2_0) {
    const Size T = residuals.size();
    std::vector<Real> h(T);
    Real h_prev = sigma2_0;
    Real e2_prev = sigma2_0;
    Real ind_prev = 0.5;  // probe-1: t=0 指示项取 1/2 (期望值)
    for (Size t = 0; t < T; ++t) {
        const Real ht = params.omega + params.alpha * e2_prev
                        + params.gamma * ind_prev * e2_prev
                        + params.beta * h_prev;
        if (!(ht > 0.0) || !std::isfinite(ht)) {
            throw std::runtime_error("filter_gjr: invalid variance");
        }
        h[t] = ht;
        e2_prev = residuals[t] * residuals[t];
        h_prev = ht;
        ind_prev = (residuals[t] < 0.0) ? 1.0 : 0.0;  // G7: I(ε<0)
    }
    return h;
}

/// @brief GJR-GARCH(1,1) QMLE 估计 (spec §2.0.4)
///
/// 算法: 均值滤波 (G2) → EWMA backcast (G1) → SLSQP 多起始优化
/// (G4/G10 约束, G16 多起始) → sandwich 协方差 (G9) → 结果组装 (G17)。
///
/// @throws std::invalid_argument 若 T < 10, 含 NaN, 或零方差
inline GjrGarchResult estimate_gjr_garch(const std::vector<Real>& data,
                                         const GarchConfig& config = GarchConfig{}) {
    const Size T = data.size();
    if (T < 10) {
        throw std::invalid_argument("estimate_gjr_garch: need at least 10 observations");
    }
    Real mean = 0.0;
    for (Real v : data) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("estimate_gjr_garch: data contains NaN/Inf");
        }
        mean += v;
    }
    mean /= static_cast<Real>(T);

    // G2: 残差 = 数据 - 样本均值
    std::vector<Real> eps(T);
    for (Size t = 0; t < T; ++t) eps[t] = data[t] - mean;

    Real var = 0.0;
    for (Real e : eps) var += e * e;
    var /= static_cast<Real>(T);
    Real scale = 0.0;
    for (Real v : data) scale = std::max(scale, std::abs(v));
    const Real rel_sd = (scale > 0.0) ? std::sqrt(var) / scale : 0.0;
    if (!(var > 0.0) || rel_sd < 1e-13) {
        throw std::invalid_argument("estimate_gjr_garch: zero variance data");
    }

    // 内部尺度变换 (同 GARCH(1,1), 递归线性 → 等变性成立)
    const Real s = 1.0 / std::sqrt(var);
    const Real log_s = std::log(s);
    std::vector<Real> eps_s(T);
    for (Size t = 0; t < T; ++t) eps_s[t] = eps[t] * s;
    const Real var_s = var * s * s;  // = 1

    const bool has_nu = config.dist != GarchDist::Normal;
    const Size k = has_nu ? 5u : 4u;

    // G1: backcast (scaled 空间)
    const Real sigma2_0 =
        backcast_variance(eps_s, config.backcast_lambda, config.backcast_window);

    // 逐观测似然: probe-1 递归 h₁=ω+(α+γ/2+β)·σ²₀ (filter_gjr)
    detail::PerObsLL pll = [&](const std::vector<Real>& theta) {
        const Real nu = has_nu ? theta[4] : 0.0;
        const GjrGarchParams p_s{theta[0], theta[1], theta[2], theta[3], nu};
        const std::vector<Real> h = filter_gjr(p_s, eps_s, sigma2_0);
        std::vector<Real> ll(T);
        for (Size t = 0; t < T; ++t) {
            ll[t] = log_likelihood_term(eps_s[t], h[t], config.dist, nu);
        }
        return ll;
    };

    // G4/G10 约束 + 边界 (spec §2.0.4 Step 3; G-gamma-sign: γ ∈ [-1,1])
    std::vector<Bounds> bnds = {{1e-8, 100.0}, {0.0, 1.0}, {-1.0, 1.0}, {0.0, 1.0}};
    if (has_nu) bnds.push_back({2.05, 500.0});
    const ConstraintFn ineq = [](const std::vector<Real>& x) {
        // G10: α + γ/2 + β < 1 (o 项系数 -0.5, arch 约束矩阵同款)
        return std::vector<Real>{x[0] - 1e-8, x[1], x[3],
                                 1.0 - x[1] - x[2] / 2.0 - x[3] - 1e-8};
    };

    // G16 多起始 (scaled): HR/ML/用户/随机扰动
    const Real nu0 = 8.0;
    std::vector<std::vector<Real>> starts;
    auto push_start = [&](Real om, Real al, Real ga, Real be) {
        std::vector<Real> st = {om, al, ga, be};
        if (has_nu) st.push_back(nu0);
        starts.push_back(std::move(st));
    };
    push_start(0.10 * var_s, 0.10, 0.05, 0.85);  // HR
    push_start(0.05 * var_s, 0.05, 0.10, 0.90);  // ML
    if (config.initial_params.size() >= 4) {
        std::vector<Real> st = {config.initial_params[0] * s * s,
                                config.initial_params[1],
                                config.initial_params[2],
                                config.initial_params[3]};
        if (has_nu) {
            st.push_back(config.initial_params.size() >= 5
                             ? config.initial_params[4] : nu0);
        }
        starts.push_back(std::move(st));
    }
    {  // 随机扰动 (HR 基础 + 10% 高斯, Philox 确定性)
        Philox4x64 rng(44, 13);
        std::vector<Real> st = {0.10 * var_s, 0.10, 0.05, 0.85};
        if (has_nu) st.push_back(nu0);
        for (Size j = 0; j < st.size(); ++j) {
            const uint64_t r1 = rng();
            const uint64_t r2 = rng();
            const Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            const Real u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            const auto [z1, z2] = box_muller(u1, u2);
            (void)z2;
            st[j] *= (1.0 + 0.10 * z1);
        }
        starts.push_back(std::move(st));
    }
    if (!config.use_multistart) starts.resize(1);

    // SLSQP 多起始
    const detail::MultistartResult fit =
        detail::run_multistart(pll, starts, bnds, ineq, config.optimizer_config);

    // 最优参数处重组 (scaled)
    const Real om_s = fit.x[0], al_s = fit.x[1], ga_s = fit.x[2], be_s = fit.x[3];
    const Real nu_s = has_nu ? fit.x[4] : 0.0;
    const std::vector<Real> h_s =
        filter_gjr(GjrGarchParams{om_s, al_s, ga_s, be_s, nu_s}, eps_s, sigma2_0);
    Real ll_scaled = 0.0;
    for (Size t = 0; t < T; ++t) {
        ll_scaled += log_likelihood_term(eps_s[t], h_s[t], config.dist, nu_s);
    }

    // 映射回原尺度
    GjrGarchResult result{};
    result.params.omega = om_s / (s * s);
    result.params.alpha = al_s;
    result.params.gamma = ga_s;
    result.params.beta = be_s;
    result.params.nu = nu_s;
    result.conditional_variances.resize(T);
    result.residuals = eps;
    result.std_residuals.resize(T);
    for (Size t = 0; t < T; ++t) {
        result.conditional_variances[t] = h_s[t] / (s * s);
        result.std_residuals[t] = eps[t] / std::sqrt(result.conditional_variances[t]);
    }
    // ℓ = ℓ' + T·log s (h_s = h·s² 的 Jacobian 修正, 同 GARCH(1,1))
    result.log_likelihood = ll_scaled + static_cast<Real>(T) * log_s;
    // G17: AIC/BIC 用完整似然, k = 4 + 分布参数
    result.aic = -2.0 * result.log_likelihood + 2.0 * static_cast<Real>(k);
    result.bic = -2.0 * result.log_likelihood
                 + static_cast<Real>(k) * std::log(static_cast<Real>(T));

    // G9: sandwich V = H⁻¹·S·H⁻¹ (scaled 空间 + J 映射, 同 GARCH(1,1))
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
            const Real scl = static_cast<Real>(Tn) / static_cast<Real>(Tn - 1);
            for (Size t = 0; t < Tn; ++t)
                for (Size a = 0; a < npar; ++a)
                    for (Size b = 0; b < npar; ++b)
                        S[a][b] += (G[a][t] - gbar[a]) * (G[b][t] - gbar[b]);
            for (Size a = 0; a < npar; ++a)
                for (Size b = 0; b < npar; ++b) S[a][b] *= scl;
            std::vector<std::vector<Real>> Hinv;
            if (!detail::invert_matrix(H, Hinv)) {
                result.message = "sandwich: Hessian inversion failed";
            } else {
                std::vector<std::vector<Real>> HS(npar, std::vector<Real>(npar, 0.0));
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
                // 映射回原尺度: 仅 ω (下标 0) 行列乘 1/s²
                result.vcov.assign(npar, std::vector<Real>(npar, 0.0));
                for (Size a = 0; a < npar; ++a)
                    for (Size b = 0; b < npar; ++b) {
                        const Real ja = (a == 0) ? 1.0 / (s * s) : 1.0;
                        const Real jb = (b == 0) ? 1.0 / (s * s) : 1.0;
                        result.vcov[a][b] = ja * V[a][b] * jb;
                    }
                result.std_errors.resize(npar);
                for (Size a = 0; a < npar; ++a) {
                    result.std_errors[a] =
                        (result.vcov[a][a] > 0.0) ? std::sqrt(result.vcov[a][a]) : 0.0;
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
