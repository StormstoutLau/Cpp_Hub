// =============================================================================
// garch_model.hpp - GARCH(1,1) QMLE 估计 (spec §2.0.2)
//
// Phase 7B v1.6 M1 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 教材锚点: Bollerslev 1986 / Tsay 3ed Ch 5 / arch arch/univariate/mean.py
// 优化器: SLSQP (ADR-018, calibration/optimizer.hpp)
//
// 幻觉点防护 (spec §6.1):
//   G1: 方差初始值 EWMA backcast, λ=0.94, τ=min(75,T); 用前 τ 个残差;
//       归一化 wᵢ=λⁱ/Σλʲ (有限样本, 非无穷级数近似 (1-λ)Σλⁱ)
//       h₁ = ω + (α+β)·σ²₀ (arch 实测约定, probe_arch_convention.py 验证:
//       差异 1.4e-20 vs h₁=σ²₀ 的 2.3e-6) — 即 ε²₀ 与 h₀ 均取 backcast
//   G2: 残差 εₜ = rₜ - μ̂ (样本均值, 非 0)
//   G3: 似然含 -0.5·log(2π) (garch_distribution.hpp 承担)
//   G4: 约束 ω>0, α≥0, β≥0, α+β<1 (SLSQP bounds + 不等式)
//   G9: QMLE sandwich 协方差 V = H⁻¹·S·H⁻¹ (S = OPG 非中心化求和)
//   G11: 标准化残差 zₜ = εₜ/√hₜ
//   G13: 多步预测见 garch_forecast.hpp
//   G16: 多起始点避免局部最优 (HR/ML/user/random-perturb)
//   G17: AIC/BIC 用完整似然
//
// 内部尺度变换 (数值稳健性): 日频收益率方差 ~1e-4, ω 量级 ~1e-6 会使 SLSQP
// 数值差分退化。内部以 s = 1/sd(ε) 缩放 (MLE 尺度等变性), 估计后映射回
// 原尺度: ω = ω'/s², α/β/ν 不变, ℓ = ℓ' - T·log(s), V = J·V'·J',
// J = diag(1/s², 1, 1, [1])。最优解与 arch 未缩放结果一致 (等变性)。
// =============================================================================
#pragma once

#include <vector>
#include <functional>
#include <cmath>
#include <string>
#include <stdexcept>
#include <algorithm>

#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/calibration/optimizer.hpp"  // SLSQP (ADR-018)
#include "cpphub/timeseries/garch/garch_distribution.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace garch {

// GARCH(1,1) 参数
struct GarchParams {
    Real omega;    ///< 常数项 (ω > 0, G4)
    Real alpha;    ///< ARCH 系数 (α ≥ 0, G4)
    Real beta;     ///< GARCH 系数 (β ≥ 0, G4)
    Real nu;       ///< t/GED 自由度 (Normal 时忽略, G-ADR5)
};

// GARCH(1,1) 估计配置
struct GarchConfig {
    GarchDist dist = GarchDist::Normal;           ///< 分布族 (G-ADR5)
    Real backcast_lambda = 0.94;                  ///< EWMA backcast λ (G-ADR1, G1)
    Size backcast_window = 0;                     ///< backcast 窗口 τ (0 => min(75, T), arch 默认)
    bool use_multistart = true;                   ///< 多起始点 (G-ADR6, G16)
    Size n_multistart = 4;                        ///< 多起始点数 (HR/ML/user/random)
    bool compute_sandwich = true;                 ///< QMLE sandwich 协方差 (G-ADR3, G9)
    bool compute_diagnostics = true;              ///< 标准化残差诊断 (G-ADR4; 经
                                                  ///< garch_diagnostics.hpp 独立调用)
    Size bootstrap_reps = 1000;                   ///< JB 检验 Bootstrap 次数 (G-ADR4)
    std::vector<Real> initial_params;             ///< 用户起始点 {ω,α,β,[ν]} (可选, G16)
    SLSQP::Config optimizer_config = [] {
        SLSQP::Config c;
        c.max_iterations = 300;   // GARCH QMLE 需比默认 100 更深
        c.ftol = 1e-9;            // 1e-10 级参数精度 (G15)
        c.xtol = 1e-10;
        return c;
    }();
};

// GARCH(1,1) 估计结果
struct GarchResult {
    GarchParams params;                           ///< 参数估计
    std::vector<Real> conditional_variances;      ///< 条件方差 hₜ 序列
    std::vector<Real> residuals;                  ///< 残差 εₜ 序列
    std::vector<Real> std_residuals;              ///< 标准化残差 zₜ = εₜ/√hₜ (G11)
    Real log_likelihood;                          ///< 对数似然值 (G3: 含常数项)
    Real aic;                                     ///< AIC (G17: 用完整似然)
    Real bic;                                     ///< BIC (G17: 用完整似然)
    std::vector<std::vector<Real>> vcov;          ///< 参数协方差矩阵 (sandwich, G9)
    std::vector<Real> std_errors;                 ///< 标准误
    bool converged;                               ///< 收敛标志
    Size n_iterations;                            ///< 迭代次数
    std::string message;                          ///< 诊断消息
};

/// @brief EWMA backcast 方差 (G1, 公开以便测试)
///
/// arch `arch/univariate/volatility.py:1161-1168` 实测:
///   - 用**前 τ 个**残差 resids[:tau] (排幻觉点 G1-dir, 非末尾)
///   - 归一化 wᵢ = λⁱ/Σⱼλʲ (排幻觉点 G1-norm, sum(w)=1, 非无穷级数近似)
///   σ²₀ = Σᵢ₌₀^{τ-1} wᵢ·ε²ᵢ
inline Real backcast_variance(const std::vector<Real>& residuals,
                              Real lambda, Size window) {
    const Size T = residuals.size();
    if (T == 0) throw std::invalid_argument("backcast_variance: empty input");
    if (!(lambda > 0.0 && lambda < 1.0)) {
        throw std::invalid_argument("backcast_variance: lambda in (0,1)");
    }
    const Size tau = (window == 0) ? std::min<Size>(75, T) : std::min(window, T);
    Real wsum = 0.0, acc = 0.0, w = 1.0;
    for (Size i = 0; i < tau; ++i) {
        wsum += w;
        acc += w * residuals[i] * residuals[i];
        w *= lambda;
    }
    return acc / wsum;
}

/// @brief GARCH(1,1) 方差递归 (公开以便测试, 与 arch compute_variance 一致)
///
/// h₁ = ω + (α+β)·σ²₀ (G1: ε²₀ 与 h₀ 均取 backcast, arch 实测约定);
/// hₜ = ω + α·ε²ₜ₋₁ + β·hₜ₋₁ (t ≥ 2)
inline std::vector<Real> filter_garch11(const GarchParams& params,
                                        const std::vector<Real>& residuals,
                                        Real sigma2_0) {
    const Size T = residuals.size();
    std::vector<Real> h(T);
    Real h_prev = sigma2_0;
    Real e2_prev = sigma2_0;
    for (Size t = 0; t < T; ++t) {
        const Real ht = params.omega + params.alpha * e2_prev
                        + params.beta * h_prev;
        if (!(ht > 0.0) || !std::isfinite(ht)) {
            throw std::runtime_error("filter_garch11: invalid variance");
        }
        h[t] = ht;
        e2_prev = residuals[t] * residuals[t];
        h_prev = ht;
    }
    return h;
}

// ---------------------------------------------------------------------------
// 共享估计机制 (GARCH/EGARCH/GJR 三模型复用)
// ---------------------------------------------------------------------------
namespace detail {

// 逐观测对数似然函数: θ → {ℓ₁, ..., ℓ_T} (scaled 空间)
// 实现须保证: 给定 θ 重算 h 路径后逐观测求 ℓₜ
using PerObsLL = std::function<std::vector<Real>(const std::vector<Real>&)>;

// 目标函数 (负对数似然, 数值防护: 递归发散/异常 → 大惩罚)
inline ObjectiveFn make_neg_ll(const PerObsLL& pll) {
    return [pll](const std::vector<Real>& theta) -> Real {
        try {
            const auto ll = pll(theta);
            Real total = 0.0;
            for (Real v : ll) {
                if (!std::isfinite(v)) return 1e100;
                total += v;
            }
            if (!std::isfinite(total)) return 1e100;
            return -total;
        } catch (...) {
            return 1e100;
        }
    };
}

// 多起始 SLSQP 公共驱动 (G16): 每个起始点独立优化, 取目标最小者
struct MultistartResult {
    std::vector<Real> x;
    Real fx;
    bool converged;
    Size n_iterations;
    std::string message;
};

inline MultistartResult run_multistart(
    const PerObsLL& pll,
    const std::vector<std::vector<Real>>& starts,
    const std::vector<Bounds>& bounds,
    const ConstraintFn& ineq,  // c(x) >= 0 (向量)
    const SLSQP::Config& cfg) {
    const ObjectiveFn f = make_neg_ll(pll);
    MultistartResult best{};
    best.fx = std::numeric_limits<Real>::infinity();
    bool any = false;
    for (const auto& x0 : starts) {
        const OptimizationResult r =
            SLSQP::minimize(f, x0, bounds, {ineq}, {}, cfg);
        any = true;
        if (r.fx < best.fx) {
            best.x = r.x;
            best.fx = r.fx;
            best.converged = r.converged;
            best.n_iterations = r.n_iterations;
            best.message = r.message;
        }
    }
    if (!any) throw std::invalid_argument("run_multistart: no start points");
    return best;
}

// 矩阵求逆 (复用 calibration detail 的 Gauss-Jordan)
inline bool invert_matrix(const std::vector<std::vector<Real>>& A,
                          std::vector<std::vector<Real>>& inv) {
    const Size n = A.size();
    inv.assign(n, std::vector<Real>(n, 0.0));
    auto work = A;
    std::vector<std::vector<Real>> rhs(n, std::vector<Real>(n, 0.0));
    for (Size i = 0; i < n; ++i) rhs[i][i] = 1.0;
    for (Size col = 0; col < n; ++col) {
        auto Ac = work;
        auto b = rhs[col];
        if (!::cpphub::v1::detail::solve_linear_system(Ac, b, n)) return false;
        for (Size i = 0; i < n; ++i) inv[i][col] = b[i];
    }
    return true;
}

// 差分步长 (相对参数尺度; floor 0.1: 原尺度小参数如 EGARCH ω≈-0.003 需
// 绝对步长 ~1e-6 抑制 llf 舍入噪声, floor 1e-3 时 d=1e-8 噪声主导 vcov)
inline Real diff_step(Real x) { return 1e-5 * std::max(std::abs(x), 0.1); }

// 数值 Hessian (总似然 ℓ = Σℓₜ; 对角中央二阶差分 + 交叉四点差分)
inline std::vector<std::vector<Real>> numerical_hessian(const PerObsLL& pll,
                                                        const std::vector<Real>& x) {
    const Size n = x.size();
    auto total = [&](const std::vector<Real>& th) {
        const auto ll = pll(th);
        Real s = 0.0;
        for (Real v : ll) s += v;
        return s;
    };
    std::vector<std::vector<Real>> H(n, std::vector<Real>(n, 0.0));
    std::vector<Real> xp = x, xm = x;
    for (Size j = 0; j < n; ++j) {
        const Real d = diff_step(x[j]);
        xp[j] = x[j] + d; xm[j] = x[j] - d;
        H[j][j] = (total(xp) - 2.0 * total(x) + total(xm)) / (d * d);
        xp[j] = x[j]; xm[j] = x[j];
    }
    for (Size i = 0; i < n; ++i) {
        for (Size j = i + 1; j < n; ++j) {
            const Real di = diff_step(x[i]);
            const Real dj = diff_step(x[j]);
            std::vector<Real> a = x, b = x, c = x, d4 = x;
            a[i] += di; a[j] += dj;
            b[i] += di; b[j] -= dj;
            c[i] -= di; c[j] += dj;
            d4[i] -= di; d4[j] -= dj;
            const Real Hij = (total(a) - total(b) - total(c) + total(d4))
                             / (4.0 * di * dj);
            H[i][j] = Hij;
            H[j][i] = Hij;
        }
    }
    return H;
}

// OPG 得分矩阵 S = Σₜ gₜgₜ' (非中心化, G9)
// gₜⱼ = [ℓₜ(θ+δeⱼ) - ℓₜ(θ-δeⱼ)]/(2δ) (h 路径随 θ 重算 — arch 数值得分模式)
// 返回 G: G[j][t] = 第 j 参数方向上观测 t 的得分
inline std::vector<std::vector<Real>> opg_scores(const PerObsLL& pll,
                                                 const std::vector<Real>& x) {
    const Size n = x.size();
    const Size T = pll(x).size();
    std::vector<std::vector<Real>> G(n, std::vector<Real>(T));
    for (Size j = 0; j < n; ++j) {
        const Real d = diff_step(x[j]);
        std::vector<Real> xp = x, xm = x;
        xp[j] += d; xm[j] -= d;
        const auto lp = pll(xp);
        const auto lm = pll(xm);
        for (Size t = 0; t < T; ++t) G[j][t] = (lp[t] - lm[t]) / (2.0 * d);
    }
    return G;
}

}  // namespace detail

/// @brief GARCH(1,1) QMLE 估计
///
/// 算法 (spec §2.0.2): 均值滤波 (G2) → EWMA backcast (G1) → SLSQP 多起始
/// 优化 (G4 约束, G16 多起始) → sandwich 协方差 (G9) → 结果组装 (G17)。
///
/// @throws std::invalid_argument 若 T < 10, 含 NaN, 或零方差
inline GarchResult estimate_garch11(const std::vector<Real>& data,
                                    const GarchConfig& config = GarchConfig{}) {
    const Size T = data.size();
    if (T < 10) {
        throw std::invalid_argument("estimate_garch11: need at least 10 observations");
    }
    Real mean = 0.0;
    for (Real v : data) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("estimate_garch11: data contains NaN/Inf");
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
    // 近常数数据防护: 样本均值舍入可留 ~1e-16 级残差 (var ~1e-32),
    // 单纯 var>0 判不住; 用相对尺度 rel_sd < 1e-13 判常数列
    Real scale = 0.0;
    for (Real v : data) scale = std::max(scale, std::abs(v));
    const Real rel_sd = (scale > 0.0) ? std::sqrt(var) / scale : 0.0;
    if (!(var > 0.0) || rel_sd < 1e-13) {
        throw std::invalid_argument("estimate_garch11: zero variance data");
    }

    // 内部尺度变换 (头文件说明): s = 1/sd, 估计在单位方差尺度
    const Real s = 1.0 / std::sqrt(var);
    const Real log_s = std::log(s);
    std::vector<Real> eps_s(T);
    for (Size t = 0; t < T; ++t) eps_s[t] = eps[t] * s;
    const Real var_s = var * s * s;  // = 1

    const bool has_nu = config.dist != GarchDist::Normal;
    const Size k = has_nu ? 4u : 3u;

    // G1: backcast (scaled 空间)
    const Real sigma2_0 =
        backcast_variance(eps_s, config.backcast_lambda, config.backcast_window);

    // 逐观测似然: GARCH(1,1) 递归 h₁=ω+(α+β)·σ²₀ (G1, filter_garch11)
    detail::PerObsLL pll = [&](const std::vector<Real>& theta) {
        const Real nu = has_nu ? theta[3] : 0.0;
        const GarchParams p_s{theta[0], theta[1], theta[2], nu};
        const std::vector<Real> h = filter_garch11(p_s, eps_s, sigma2_0);
        std::vector<Real> ll(T);
        for (Size t = 0; t < T; ++t) {
            ll[t] = log_likelihood_term(eps_s[t], h[t], config.dist, nu);
        }
        return ll;
    };

    // G4 约束 + 边界 (arch StudentsT 边界: ν ∈ [2.05, 500])
    const std::vector<Bounds> bounds = [] {
        std::vector<Bounds> b = {{1e-8, 100.0}, {0.0, 1.0}, {0.0, 1.0}};
        return b;
    }();
    std::vector<Bounds> bnds = bounds;
    if (has_nu) bnds.push_back({2.05, 500.0});
    const ConstraintFn ineq = [has_nu](const std::vector<Real>& x) {
        std::vector<Real> c = {x[0] - 1e-8, x[1], x[2],
                               1.0 - x[1] - x[2] - 1e-8};
        return c;
    };

    // G16 多起始 (scaled): HR 高频 / ML / 用户 / 随机扰动
    std::vector<std::vector<Real>> starts;
    const Real nu0 = 8.0;  // arch StudentsT 默认起始
    auto push_start = [&](Real om, Real al, Real be) {
        std::vector<Real> st = {om, al, be};
        if (has_nu) st.push_back(nu0);
        starts.push_back(std::move(st));
    };
    push_start(0.10 * var_s, 0.10, 0.85);  // HR
    push_start(0.05 * var_s, 0.05, 0.90);  // ML
    if (config.initial_params.size() >= 3) {
        std::vector<Real> st = {config.initial_params[0] * s * s,
                                config.initial_params[1],
                                config.initial_params[2]};
        if (has_nu && config.initial_params.size() >= 4) {
            st.push_back(config.initial_params[3]);
        } else if (has_nu) {
            st.push_back(nu0);
        }
        starts.push_back(std::move(st));
    }
    {  // 随机扰动 (HR 基础 + 10% 高斯, Philox 确定性)
        Philox4x64 rng(42, 7);
        std::vector<Real> st = {0.10 * var_s, 0.10, 0.85};
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
    const Real om_s = fit.x[0], al_s = fit.x[1], be_s = fit.x[2];
    const Real nu_s = has_nu ? fit.x[3] : 0.0;
    const std::vector<Real> h_s =
        filter_garch11(GarchParams{om_s, al_s, be_s, nu_s}, eps_s, sigma2_0);
    Real ll_scaled = 0.0;
    for (Size t = 0; t < T; ++t) {
        ll_scaled += log_likelihood_term(eps_s[t], h_s[t], config.dist, nu_s);
    }

    // 映射回原尺度
    GarchResult result{};
    result.params.omega = om_s / (s * s);
    result.params.alpha = al_s;
    result.params.beta = be_s;
    result.params.nu = nu_s;
    result.conditional_variances.resize(T);
    result.residuals = eps;
    result.std_residuals.resize(T);
    for (Size t = 0; t < T; ++t) {
        result.conditional_variances[t] = h_s[t] / (s * s);
        result.std_residuals[t] = eps[t] / std::sqrt(result.conditional_variances[t]);
    }
    // ℓ'ₜ = ℓₜ - log s (h_s = h·s² → log 项移动 -log s), 故 ℓ = ℓ' + T·log s
    result.log_likelihood = ll_scaled + static_cast<Real>(T) * log_s;
    // G17: AIC/BIC 用完整似然, k = 3 + 分布参数
    result.aic = -2.0 * result.log_likelihood + 2.0 * static_cast<Real>(k);
    result.bic = -2.0 * result.log_likelihood
                 + static_cast<Real>(k) * std::log(static_cast<Real>(T));

    // G9: sandwich V = H⁻¹·S·H⁻¹ (scaled 空间), J = diag(1/s²,1,1,[1]) 映射
    // S 定义 (arch base.py:977-988 对齐, spec Step5 注授权调整):
    //   S = [n/(n-1)]·Σₜ(gₜ-ḡ)(gₜ-ḡ)' (中心化 OPG, 与 arch np.cov(ddof=1) 一致)
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
            const Real scale = static_cast<Real>(Tn) / static_cast<Real>(Tn - 1);
            for (Size t = 0; t < Tn; ++t) {
                for (Size a = 0; a < npar; ++a)
                    for (Size b = 0; b < npar; ++b) {
                        S[a][b] += (G[a][t] - gbar[a]) * (G[b][t] - gbar[b]);
                    }
            }
            for (Size a = 0; a < npar; ++a)
                for (Size b = 0; b < npar; ++b) S[a][b] *= scale;
            std::vector<std::vector<Real>> Hinv;
            if (!detail::invert_matrix(H, Hinv)) {
                result.message = "sandwich: Hessian inversion failed";
            } else {
                // V = H⁻¹·S·H⁻¹
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
                // 映射回原尺度: ω 行列乘 1/s² (α/β/ν 尺度不变)
                result.vcov.resize(npar, std::vector<Real>(npar, 0.0));
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
