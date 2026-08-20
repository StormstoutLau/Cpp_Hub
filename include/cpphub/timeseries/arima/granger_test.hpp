// =============================================================================
// granger_test.hpp - 线性 Granger 因果检验 (标准 F/χ²/LR + TY + HAC-Wald)
//
// Phase 7C v1.7 M1 (PHASE7C_SPEC.md v1.2 §3.2, 决策 5/6)
//
// 效应方程: effect_t = c + Σφᵢ·effect_{t−i} + Σβⱼ·cause_{t−j} + ε_t
//   H0: β₁ = … = β_p = 0 (cause 不 Granger-cause effect)
//
// 对照基准 (granger_statsmodels_baselines.txt, statsmodels 0.14.4):
//   4 统计量 vs grangercausalitytests 1e-10 (tsa/stattools.py L1473-1699
//   一手源码; lagmat2ds trim="both" 实测: 行 t = p..T−1, nobs = T − p)
//   TY Wald vs statsmodels OLS f_test (增广回归) 1e-8
//   HAC-Wald vs 显式 NW 三明治 (cov_hac_simple use_correction=False) 1e-8
//
// 幻觉点防护 (spec §9.1):
//   GR1 (df 公式): df_denom = nobs − (2p+1), nobs = T − p (非 T − 2p);
//        F = (SSR_r − SSR_u)/p · df_denom/SSR_u
//   GR2 (TY df=k): Wald χ² 自由度恒为 k (lag), 非 k + d_max
//   GR3 (d_max 外部): 不做模型内单位根自适应, 调用方预检验后传入
//   GR4 (增广阶不进约束): TY 估计 k+d_max 阶但约束矩阵仅前 k 阶 cause 滞后
//   GR5 (稳健须自建): statsmodels 四统计量均非异方差稳健; HAC-Wald 为
//        NW Bartlett 三明治 (w_l = 1 − l/(L+1)) 上的 Wald, χ²(k)
//   GR6 (方向): 显式 (cause, effect) 形参 — statsmodels 输入为
//        [effect, cause] (第一列被解释), 方向陷阱由此而来
//   GR7 (I(1) 失效): 水平值标准 F 分布推断失效 → 差分或 TY (集成场景 2)
//
// 数值路径: 正规方程 + Gauss-Jordan (probe_granger_prec.py 实测
//   夹具条件数 ≤ 43, 误差 ~1e-13 ≪ 容差 1e-10; QR 不必要)
// =============================================================================

#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/special_functions.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace granger {

/// Granger 因果检验结果 (spec §3.2 逐字段)
struct GrangerResult {
    // ---- 标准 4 统计量 (与 statsmodels grangercausalitytests 一一配对) ----
    Real f_stat = std::numeric_limits<Real>::quiet_NaN();  ///< ssr_ftest (主)
    Real f_p = std::numeric_limits<Real>::quiet_NaN();
    Real params_f_stat =
        std::numeric_limits<Real>::quiet_NaN();  ///< params_ftest (Wald F,
                                                 ///< 与 ssr_ftest 数学等价
                                                 ///< 的独立数值路径)
    Real params_f_p = std::numeric_limits<Real>::quiet_NaN();
    Real chi2_stat =
        std::numeric_limits<Real>::quiet_NaN();  ///< ssr_chi2test
    Real chi2_p = std::numeric_limits<Real>::quiet_NaN();
    Real lr_stat = std::numeric_limits<Real>::quiet_NaN();  ///< lrtest
    Real lr_p = std::numeric_limits<Real>::quiet_NaN();
    // ---- Toda-Yamamoto 增广 Wald (d_max > 0 时有效, GR2/GR3/GR4) ----
    Real ty_wald_stat = std::numeric_limits<Real>::quiet_NaN();
    Real ty_wald_p = std::numeric_limits<Real>::quiet_NaN();
    bool has_ty = false;
    // ---- HAC 稳健 Wald (with_hac = true 时有效, GR5) ----
    Real hac_wald_stat = std::numeric_limits<Real>::quiet_NaN();
    Real hac_wald_p = std::numeric_limits<Real>::quiet_NaN();
    bool has_hac = false;
    // ---- 自由度 (GR1) ----
    Size df1 = 0;  ///< df_num = p
    Size df2 = 0;  ///< df_denom = nobs − (2p+1), nobs = T − p
    std::string summary;  ///< 注明: 均非异方差稳健 (除 HAC 版)
};

namespace detail {

/// OLS 完整输出 (granger 需要 (X'X)⁻¹ 供 Wald/HAC, unit_root::detail::ols_fit
/// 不返回之; 求解/求逆复用同一 Gauss-Jordan 路径)
struct OlsJoint {
    std::vector<Real> beta;              ///< 系数 (k)
    std::vector<Real> resid;             ///< 残差 (n)
    std::vector<std::vector<Real>> xtx_inv;  ///< (X'X)⁻¹ (k×k)
    Real ssr = 0.0;
    Real sigma2 = 0.0;                   ///< SSR/(n−k) (df 修正)
    Size nobs = 0;
    Size k = 0;
};

inline OlsJoint ols_joint(const std::vector<Real>& y,
                          const std::vector<std::vector<Real>>& X) {
    const Size n = y.size();
    const Size k = X[0].size();
    if (n <= k) {
        throw std::invalid_argument("granger ols: insufficient observations");
    }

    // 正规方程
    std::vector<std::vector<Real>> XtX(k, std::vector<Real>(k, 0.0));
    std::vector<Real> Xty(k, 0.0);
    for (Size i = 0; i < k; ++i) {
        for (Size j = i; j < k; ++j) {
            Real s = 0.0;
            for (Size t = 0; t < n; ++t) s += X[t][i] * X[t][j];
            XtX[i][j] = s;
            XtX[j][i] = s;
        }
        Real s = 0.0;
        for (Size t = 0; t < n; ++t) s += X[t][i] * y[t];
        Xty[i] = s;
    }

    // Gauss-Jordan 消元解 beta (增广)
    std::vector<std::vector<Real>> A = XtX;
    std::vector<Real> b = Xty;
    for (Size col = 0; col < k; ++col) {
        Size piv = col;
        Real maxv = std::fabs(A[col][col]);
        for (Size r = col + 1; r < k; ++r) {
            if (std::fabs(A[r][col]) > maxv) {
                maxv = std::fabs(A[r][col]);
                piv = r;
            }
        }
        if (maxv < 1e-300) {
            throw std::runtime_error("granger ols: singular design matrix");
        }
        if (piv != col) {
            std::swap(A[piv], A[col]);
            std::swap(b[piv], b[col]);
        }
        const Real akk = A[col][col];
        for (Size r = col + 1; r < k; ++r) {
            const Real f = A[r][col] / akk;
            if (f == 0.0) continue;
            for (Size c = col; c < k; ++c) A[r][c] -= f * A[col][c];
            b[r] -= f * b[col];
        }
    }
    OlsJoint res;
    res.beta.resize(k);
    for (Size i = k; i-- > 0;) {
        Real s = b[i];
        for (Size j = i + 1; j < k; ++j) s -= A[i][j] * res.beta[j];
        res.beta[i] = s / A[i][i];
    }

    // 残差与 SSR
    res.resid.resize(n);
    res.ssr = 0.0;
    for (Size t = 0; t < n; ++t) {
        Real pred = 0.0;
        for (Size j = 0; j < k; ++j) pred += res.beta[j] * X[t][j];
        res.resid[t] = y[t] - pred;
        res.ssr += res.resid[t] * res.resid[t];
    }

    // (X'X)⁻¹ (XtX 已被消元破坏 → 用副本)
    unit_root::detail::invert_matrix(XtX);
    res.xtx_inv = std::move(XtX);
    res.nobs = n;
    res.k = k;
    res.sigma2 = res.ssr / static_cast<Real>(n - k);
    return res;
}

/// 线性 Wald 二次型: b'(V)⁻¹b, V 为 dim×dim (调用方保证对称正定)
inline Real wald_quadratic(const std::vector<Real>& bsub,
                           std::vector<std::vector<Real>> V) {
    const Size m = bsub.size();
    if (V.size() != m) {
        throw std::invalid_argument("wald_quadratic: dimension mismatch");
    }
    unit_root::detail::invert_matrix(V);
    Real w = 0.0;
    for (Size i = 0; i < m; ++i) {
        Real s = 0.0;
        for (Size j = 0; j < m; ++j) s += V[i][j] * bsub[j];
        w += bsub[i] * s;
    }
    return w;
}

/// NW 经验法则带宽 (hac_bandwidth=0 默认): floor(4·(T/100)^(2/9))
/// (与 econometrics select_max_lag NW 规则一致, E4; T = 回归有效样本)
inline Size nw_default_bandwidth(Size nobs) {
    const Real t = static_cast<Real>(nobs);
    return static_cast<Size>(
        std::floor(4.0 * std::pow(t / 100.0, 2.0 / 9.0)));
}

/// @brief NW-Bartlett 三明治 Wald on cause 滞后子块 (自包含, 供公开接口
///        与 White 退化 (L=0) 基准测试共用)
/// 设计 X = [own lags 1..p, cause lags 1..p, const], 行 t = p..T−1;
/// Ω = Σ_{l=0}^{L} w_l·(Γ_l + Γ_l'), w_l = 1 − l/(L+1),
///   Γ_l[i,j] = Σ_t X[t−l,i]·e_{t−l}·X[t,j]·e_t (statsmodels cov_hac_simple
///   use_correction=False 约定, 无小样本修正);
/// V = (X'X)⁻¹Ω(X'X)⁻¹; W = b_c'(V_c)⁻¹b_c, b_c = cause 系数块 [p, 2p)
/// @param L 显式带宽; **L=0 → White 退化** (Ω = Σ_t xu_t xu_t', 循环不执行;
///          公开 API 的 hac_bandwidth=0 语义是 NW 默认, White 仅此处可达)
/// @throws std::invalid_argument L ≥ nobs 或样本不足
inline Real hac_wald_statistic(const std::vector<Real>& cause,
                               const std::vector<Real>& effect,
                               Size p, Size L) {
    const Size T = effect.size();
    const Size nobs = T - p;
    const Size k = 2 * p + 1;
    if (static_cast<std::ptrdiff_t>(nobs) <=
        static_cast<std::ptrdiff_t>(k)) {
        throw std::invalid_argument("hac_wald_statistic: insufficient obs");
    }
    if (L >= nobs) {
        throw std::invalid_argument("hac_wald_statistic: L >= nobs");
    }
    std::vector<Real> lhs(nobs);
    std::vector<std::vector<Real>> X(nobs, std::vector<Real>(k));
    for (Size i = 0; i < nobs; ++i) {
        const Size t = p + i;
        lhs[i] = effect[t];
        for (Size j = 1; j <= p; ++j) X[i][j - 1] = effect[t - j];
        for (Size j = 1; j <= p; ++j) X[i][p + j - 1] = cause[t - j];
        X[i][2 * p] = 1.0;
    }
    const auto u = ols_joint(lhs, X);

    // Ω = Σ_t xu_t xu_t' + Σ_{l=1}^{L} w_l·(Γ_l + Γ_l')
    std::vector<std::vector<Real>> Om(k, std::vector<Real>(k, 0.0));
    for (Size i = 0; i < k; ++i) {
        for (Size j = i; j < k; ++j) {
            Real s = 0.0;
            for (Size t = 0; t < nobs; ++t) {
                s += X[t][i] * u.resid[t] * u.resid[t] * X[t][j];
            }
            for (Size l = 1; l <= L; ++l) {
                const Real w = 1.0 - static_cast<Real>(l) /
                                          static_cast<Real>(L + 1);
                Real c = 0.0;
                for (Size t = l; t < nobs; ++t) {
                    c += X[t - l][i] * u.resid[t - l] * u.resid[t] *
                         X[t][j];
                    c += X[t - l][j] * u.resid[t - l] * u.resid[t] *
                         X[t][i];
                }
                s += w * c;
            }
            Om[i][j] = s;
            Om[j][i] = s;
        }
    }
    // V = (X'X)⁻¹ Ω (X'X)⁻¹
    std::vector<std::vector<Real>> V(k, std::vector<Real>(k, 0.0));
    for (Size i = 0; i < k; ++i) {
        for (Size j = 0; j < k; ++j) {
            Real s = 0.0;
            for (Size a = 0; a < k; ++a) {
                for (Size b = 0; b < k; ++b) {
                    s += u.xtx_inv[i][a] * Om[a][b] * u.xtx_inv[b][j];
                }
            }
            V[i][j] = s;
        }
    }
    std::vector<Real> bsub(p);
    std::vector<std::vector<Real>> Vsub(p, std::vector<Real>(p, 0.0));
    for (Size i = 0; i < p; ++i) {
        bsub[i] = u.beta[p + i];
        for (Size j = 0; j < p; ++j) Vsub[i][j] = V[p + i][p + j];
    }
    return wald_quadratic(bsub, std::move(Vsub));
}

}  // namespace detail

/// @brief 线性 Granger 因果检验 (标准 F/χ²/LR + TY 增广 Wald + HAC-Wald)
/// @param cause  原因序列 (滞后项进入 effect 方程的候选解释变量)
/// @param effect 效应序列 (被解释变量; ⚠️ 参数顺序显式 (cause, effect),
///               消灭 statsmodels "第二列 cause 第一列" 方向陷阱 GR6)
/// @param lag    滞后阶 p (≥1; cause 与 effect 同阶, statsmodels 双变量约定)
/// @param d_max  TY 增广阶 (GR3: 外部单位根预检验给定, I(1) 惯例 1;
///               0 = 不做 TY → ty_wald_* = NaN)
/// @param with_hac 是否计算 HAC-Wald (false → hac_wald_* = NaN)
/// @param hac_bandwidth NW 带宽 L (0 = NW 经验法则默认; Bartlett w_l =
///               1 − l/(L+1), statsmodels cov_hac_simple 无小样本修正)
/// @return GrangerResult (§1.4-5 NaN 政策)
/// @throws std::invalid_argument 长度不齐/NaN/lag=0/样本不足/带宽越界
inline GrangerResult granger_test(const std::vector<Real>& cause,
                                  const std::vector<Real>& effect,
                                  Size lag = 1,
                                  Size d_max = 0,
                                  bool with_hac = true,
                                  Size hac_bandwidth = 0) {
    const Size T = effect.size();
    if (cause.size() != T) {
        throw std::invalid_argument("granger_test: length mismatch");
    }
    if (T == 0) {
        throw std::invalid_argument("granger_test: empty series");
    }
    if (lag == 0) {
        throw std::invalid_argument("granger_test: lag >= 1 required");
    }
    for (Size t = 0; t < T; ++t) {
        if (!std::isfinite(cause[t]) || !std::isfinite(effect[t])) {
            throw std::invalid_argument("granger_test: NaN/inf input");
        }
    }

    // ---- 标准回归 (statsmodels lagmat2ds 语义: 行 t = p..T−1, nobs = T−p) ----
    const Size p = lag;
    const Size nobs = T - p;                    // ≥ 2p+2 需校验
    const Size k_u = 2 * p + 1;                 // own + cause + const (GR1)
    if (static_cast<std::ptrdiff_t>(nobs) <=
        static_cast<std::ptrdiff_t>(k_u)) {
        throw std::invalid_argument(
            "granger_test: insufficient observations (need T > 3*lag + 1)");
    }

    std::vector<Real> lhs(nobs);
    std::vector<std::vector<Real>> Xr(nobs), Xu(nobs);
    for (Size i = 0; i < nobs; ++i) {
        const Size t = p + i;  // 0-based 当期
        lhs[i] = effect[t];
        Xr[i].reserve(p + 1);
        Xu[i].reserve(k_u);
        for (Size j = 1; j <= p; ++j) {
            const Real own = effect[t - j];
            Xr[i].push_back(own);
            Xu[i].push_back(own);
        }
        for (Size j = 1; j <= p; ++j) Xu[i].push_back(cause[t - j]);
        Xr[i].push_back(1.0);  // const 末列 (add_constant prepend=False)
        Xu[i].push_back(1.0);
    }

    const auto ur = detail::ols_joint(lhs, Xr);
    const auto uu = detail::ols_joint(lhs, Xu);

    GrangerResult res;
    res.df1 = p;
    res.df2 = nobs - k_u;
    const Real df_denom = static_cast<Real>(res.df2);

    // ssr_ftest (GR1): F = (SSR_r − SSR_u)/p · df_denom/SSR_u
    const Real ssr_diff = ur.ssr - uu.ssr;
    res.f_stat = ssr_diff / uu.ssr / static_cast<Real>(p) * df_denom;
    res.f_p = econometrics::detail::f_sf(static_cast<Real>(p), df_denom,
                                         res.f_stat);

    // ssr_chi2test: χ² = nobs·(SSR_r − SSR_u)/SSR_u
    res.chi2_stat = static_cast<Real>(nobs) * ssr_diff / uu.ssr;
    res.chi2_p =
        econometrics::detail::chi2_sf(static_cast<Real>(p), res.chi2_stat);

    // lrtest: LR = nobs·(ln SSR_r − ln SSR_u)
    res.lr_stat =
        static_cast<Real>(nobs) * (std::log(ur.ssr) - std::log(uu.ssr));
    res.lr_p =
        econometrics::detail::chi2_sf(static_cast<Real>(p), res.lr_stat);

    // params_ftest: cause 滞后子块 Wald F (V = σ̂²·(X'X)⁻¹, 独立数值路径)
    {
        std::vector<Real> bsub(p);
        std::vector<std::vector<Real>> Vsub(
            p, std::vector<Real>(p, 0.0));
        for (Size i = 0; i < p; ++i) {
            bsub[i] = uu.beta[p + i];
            for (Size j = 0; j < p; ++j) {
                Vsub[i][j] = uu.sigma2 * uu.xtx_inv[p + i][p + j];
            }
        }
        const Real w = detail::wald_quadratic(bsub, std::move(Vsub));
        res.params_f_stat = w / static_cast<Real>(p);
        res.params_f_p = econometrics::detail::f_sf(
            static_cast<Real>(p), df_denom, res.params_f_stat);
    }

    // ---- TY 增广 Wald (决策 5): 估计 p* = p + d_max 阶, 约束仅前 p 阶
    //      cause 滞后 (GR4), df = p (GR2) ----
    if (d_max > 0) {
        const Size ps = p + d_max;
        const Size n_ty = T - ps;
        const Size k_ty = 2 * ps + 1;
        if (static_cast<std::ptrdiff_t>(n_ty) <=
            static_cast<std::ptrdiff_t>(k_ty)) {
            throw std::invalid_argument(
                "granger_test: insufficient observations for TY "
                "(need T > 3*(lag+d_max) + 1)");
        }
        std::vector<Real> lhs_ty(n_ty);
        std::vector<std::vector<Real>> Xt(n_ty,
                                          std::vector<Real>(k_ty));
        for (Size i = 0; i < n_ty; ++i) {
            const Size t = ps + i;
            lhs_ty[i] = effect[t];
            for (Size j = 1; j <= ps; ++j) Xt[i][j - 1] = effect[t - j];
            for (Size j = 1; j <= ps; ++j) Xt[i][ps + j - 1] = cause[t - j];
            Xt[i][2 * ps] = 1.0;
        }
        const auto ut = detail::ols_joint(lhs_ty, Xt);
        std::vector<Real> bsub(p);
        std::vector<std::vector<Real>> Vsub(p,
                                            std::vector<Real>(p, 0.0));
        for (Size i = 0; i < p; ++i) {
            bsub[i] = ut.beta[ps + i];  // cause 滞后 1..p (GR4)
            for (Size j = 0; j < p; ++j) {
                Vsub[i][j] = ut.sigma2 * ut.xtx_inv[ps + i][ps + j];
            }
        }
        res.ty_wald_stat = detail::wald_quadratic(bsub, std::move(Vsub));
        res.ty_wald_p = econometrics::detail::chi2_sf(
            static_cast<Real>(p), res.ty_wald_stat);
        res.has_ty = true;
    }

    // ---- HAC-Wald (决策 6): NW Bartlett 三明治上的 Wald, χ²(p) (GR5) ----
    if (with_hac) {
        Size L = hac_bandwidth;
        if (L == 0) {
            L = detail::nw_default_bandwidth(nobs);
        }
        if (L >= nobs) {
            throw std::invalid_argument(
                "granger_test: hac_bandwidth >= nobs");
        }
        res.hac_wald_stat =
            detail::hac_wald_statistic(cause, effect, p, L);
        res.hac_wald_p = econometrics::detail::chi2_sf(
            static_cast<Real>(p), res.hac_wald_stat);
        res.has_hac = true;
    }

    res.summary = "Granger causality (lag=" + std::to_string(p) +
                  ", nobs=" + std::to_string(nobs) + "): standard "
                  "F/chi2/LR non-robust" +
                  (res.has_ty ? "; TY augmented Wald df=k" : "") +
                  (res.has_hac ? "; HAC-Wald NW robust" : "");
    return res;
}

}  // namespace granger
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
