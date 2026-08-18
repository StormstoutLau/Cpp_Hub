// =============================================================================
// midas_model.hpp - MIDAS-DL/AR(+AR*)/U-MIDAS 估计: 集中化 NLS (spec §6.3)
//
// Phase 7C v1.7 M4 (PHASE7C_SPEC.md v1.2 §6.3, 决策 22/23)
//
// 模型 (决策 22, 逐字复刻 midasr 公式语义):
//   MIDAS-DL:  y_t = μ + Σ_i θ_i·w_i(λ)·x_{tm−i+1} + ε_t
//              (w_1 ↔ h=0 期末最新 — W-dir 方向裁决 2026-08-18 定案:
//               probe_midas_form.R Form A 恢复 λ*=(5.03,−0.499) ≈ 真值)
//   MIDAS-AR:  + Σ_ℓ ρ_ℓ·y_{t−ℓ}   (ℓ = 1..max_low_lag, 普通无约束 AR 列)
//   MIDAS-AR*: ρ_ℓ = (1−φ)·φ^{ℓ−1} 单参数参数化 (外层 φ, spec §6.3)
//   U-MIDAS:   每个高频滞后无约束系数 (纯 OLS)
//
// 估计 — 集中化 NLS (决策 23, Ghysels-Qian 2019 profile 思想):
//   外层: SLSQP (ADR-018) 仅优化非线性形状超参 θ
//     Nealmon: θ = (λ₂..λ_K); NBeta: θ = (κ₁,κ₂,θ₀); AlmonP/HarStep: θ 空
//     (AlmonP/HarStep 权重关于参数线性 → 约束 OLS 解析等价, 无外层)
//   内层: 给定 θ → 线性列 Z(θ) → OLS(QR) 解 (μ, δ, ρ...) → 集中化 SSR(θ)
//   PolyStep: 断点为外部给定离散结构非连续参数 → midas_fit 不支持
//     (weights 逐点对照仍可用 polystep_weights)
//
// 对照基准 (verify_midas.R 2026-08-18, 新 DGP y/x 同数组):
//   W3 U-MIDAS:  coef=(2.0021768, 2.2389843, 1.4483865, 0.8124879, 0.5114662)
//                SSR=7.3324398, n_eff=99  — 1e-10 主锚
//   W4 NLS DL:   (μ,δ,λ₂)=(2.0044541, 5.0288963, −0.4945359), SSR=7.4968752
//                — 1e-6~1e-8 (midasr BFGS reltol=1e-12)
//   W5 随机起始: 同最优 (多起点收敛唯一)
//   W6 MIDAS-AR: (2.0353177, 5.0689495, −0.4797953, −0.0165488)
//
// 对齐 (MD3): design_matrix 行 j ↔ R mls na.omit 后期 j 逐行一致;
//   AR 列行 j 滞后 ℓ = y[j−1−ℓ] (0-based) ↔ R mls(y, ℓ, 1)
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
#include "cpphub/timeseries/midas/midas_weights.hpp"
#include "cpphub/timeseries/midas/mixed_freq_data.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace midas {

/// 模型类型 (spec §6.3)
enum class MidasType { DL, AR, ARStar, UMidas };

/// MIDAS 估计结果 (spec §6.3 逐字段)
struct MidasResult {
    Real intercept = 0.0;         ///< μ
    std::vector<Real> delta;      ///< 内层线性参数: (δ 尺度, AR 系数 ρ...);
                                  ///<   AlmonP/HarStep = 权重参数; U-MIDAS = X 系数
    std::vector<Real> lambda;     ///< 外层形状超参最终值 (Nealmon: λ₂..;
                                  ///<   NBeta: (κ₁,κ₂,θ₀); 线性族/U-MIDAS: 空)
    std::vector<Real> midas_coef; ///< 隐含高频滞后系数 (含尺度; U-MIDAS 同 delta_x)
    Real sigma2 = 0.0;            ///< SSR/(n−k) (df 修正, midas_u summary 同)
    Real loglik = 0.0;            ///< 高斯 MLE loglik (σ̂² = SSR/n)
    Real ssr = 0.0;
    bool converged = false;
    std::string message;
    Size n_obs = 0;
    MidasType type{};
    std::vector<Real> residuals, fitted;
    // hAh 权重检验 (K-Z 2012, MD6): 由 midas_diagnostics.hpp 填充, 三列全录
    Real hah_stat = 0.0;
    Real hah_p = std::numeric_limits<Real>::quiet_NaN();
    Size hah_df = 0;
};

namespace detail {

/// 权重基列 (给定外层 θ): 返回 n×Kb 线性列, 回归系数 = 权重参数
/// Nealmon/NBeta: 单列 = X·w̄(θ) (δ 归一 1, 内层系数 = δ)
/// AlmonP: Kb = K 列, 基 b_k[i] = (i+1)^k (midasr raw poly, i 从 1 起 MD1)
/// HarStep: Kb = 3 列 (B1=e₁, B2=(0,1/5×4,0…), B3=(1/20)×20, d=20)
inline std::vector<std::vector<Real>> weight_basis_columns(
    MidasWeight w, Size k_high, const std::vector<Real>& theta,
    const std::vector<std::vector<Real>>& X) {
    const Size n = X.size();
    std::vector<std::vector<Real>> cols;
    switch (w) {
    case MidasWeight::Nealmon: {
        // theta = (λ₂..λ_K); w̄ = nealmon(δ=1, θ)
        std::vector<Real> lam{1.0};
        lam.insert(lam.end(), theta.begin(), theta.end());
        const std::vector<Real> wbar = nealmon_weights(lam, k_high);
        cols.push_back(std::vector<Real>(n, 0.0));
        for (Size r = 0; r < n; ++r) {
            Real s = 0.0;
            for (Size i = 0; i < k_high; ++i) s += wbar[i] * X[r][i];
            cols[0][r] = s;
        }
        break;
    }
    case MidasWeight::NBeta: {
        if (theta.size() != 3) {
            throw std::invalid_argument("NBeta: theta = (k1, k2, th0)");
        }
        const std::vector<Real> wbar =
            nbeta_weights({1.0, theta[0], theta[1], theta[2]}, k_high);
        cols.push_back(std::vector<Real>(n, 0.0));
        for (Size r = 0; r < n; ++r) {
            Real s = 0.0;
            for (Size i = 0; i < k_high; ++i) s += wbar[i] * X[r][i];
            cols[0][r] = s;
        }
        break;
    }
    case MidasWeight::AlmonP: {
        // 基列数 K 由 theta 大小定 (默认 2: w = c + a·i); theta 空 → K=2
        const Size K = theta.empty() ? 2 : theta.size();
        for (Size k = 0; k < K; ++k) {
            cols.push_back(std::vector<Real>(n, 0.0));
            for (Size r = 0; r < n; ++r) {
                Real s = 0.0;
                for (Size i = 0; i < k_high; ++i) {
                    Real ip = 1.0;
                    for (Size p = 0; p < k; ++p) ip *= static_cast<Real>(i + 1);
                    s += ip * X[r][i];
                }
                cols[k][r] = s;
            }
        }
        break;
    }
    case MidasWeight::HarStep: {
        if (k_high != 20) {
            throw std::invalid_argument("HarStep: k_high = 20 (midasr)");
        }
        // B1 = e1; B2 = (0,1/5,1/5,1/5,1/5,0…); B3 = (1/20)×20
        for (Size k = 0; k < 3; ++k) cols.push_back(std::vector<Real>(n, 0.0));
        for (Size r = 0; r < n; ++r) {
            cols[0][r] = X[r][0];
            Real s2 = 0.0, s3 = 0.0;
            for (Size i = 1; i <= 4; ++i) s2 += X[r][i];
            for (Size i = 0; i < 20; ++i) s3 += X[r][i];
            cols[1][r] = s2 / 5.0;
            cols[2][r] = s3 / 20.0;
        }
        break;
    }
    case MidasWeight::PolyStep:
        throw std::invalid_argument(
            "PolyStep: step breakpoints are a discrete external structure, "
            "not continuously parameterizable — midas_fit unsupported "
            "(use polystep_weights for pointwise weights)");
    }
    return cols;
}

/// θ 维度: Nealmon = len(λ₂..) (默认 1); NBeta = 3; AlmonP/HarStep = 0
inline Size theta_dim(MidasWeight w, const std::vector<std::vector<Real>>& starts) {
    switch (w) {
    case MidasWeight::Nealmon:
        return starts.empty() ? 1 : starts[0].size();
    case MidasWeight::NBeta:
        return 3;
    case MidasWeight::AlmonP:
    case MidasWeight::HarStep:
    case MidasWeight::PolyStep:
        return 0;
    }
    return 0;
}

/// 默认多起点 (starts 空; MD8: 避开 λ=0 单一起点的平坦陷阱, 网格含三形状)
inline std::vector<std::vector<Real>> default_starts(
    MidasWeight w, const std::vector<std::vector<Real>>& user_starts) {
    if (!user_starts.empty()) return user_starts;
    switch (w) {
    case MidasWeight::Nealmon: {
        // 递减 (λ₂<0) / 均匀 (0) / 递增 三形状网格
        return {{-2.0}, {-1.0}, {-0.5}, {-0.1}, {0.0}, {0.1}, {0.5}};
    }
    case MidasWeight::NBeta: {
        return {{1.5, 1.5, 0.0}, {2.0, 2.0, 0.1}, {3.0, 2.0, 0.0},
                {2.0, 3.0, 0.05}, {1.0, 2.0, -0.1}};
    }
    default:
        return {{}};
    }
}

/// θ 的 SLSQP bounds
inline std::vector<Bounds> theta_bounds(MidasWeight w, Size dim) {
    switch (w) {
    case MidasWeight::Nealmon:
        return std::vector<Bounds>(dim, Bounds{-10.0, 10.0});
    case MidasWeight::NBeta:
        return {Bounds{0.05, 20.0}, Bounds{0.05, 20.0}, Bounds{-0.95, 0.95}};
    default:
        return {};
    }
}

}  // namespace detail

/// @brief MIDAS 估计主入口 (集中化 NLS / U-MIDAS 纯 OLS)
/// @param data 混频数据 (y n 期 + x n·m; design 行 j=⌊kmax/m⌋+1..n ↔ midasr
///        na.omit 后样本逐行一致)
/// @param weight 权重族 (U-MIDAS 忽略)
/// @param type DL / AR (max_low_lag 生效) / ARStar / UMidas
/// @param k_high 高频滞后数 d
/// @param starts 外层 θ 多起点 (空 ⇒ 默认网格, MD8)
/// @param seed 预留 (默认网格确定性, 无随机起始)
inline MidasResult midas_fit(const MixedFreqData& data, MidasWeight weight,
                             MidasType type, Size k_high = 1,
                             const std::vector<std::vector<Real>>& starts = {},
                             Size /*seed*/ = 42) {
    data.validate();
    if (k_high < 1) throw std::invalid_argument("midas_fit: k_high >= 1");
    if (weight == MidasWeight::PolyStep && type != MidasType::UMidas) {
        throw std::invalid_argument(
            "midas_fit: PolyStep unsupported (discrete breakpoints)");
    }
    const Size p_ar = (type == MidasType::DL || type == MidasType::UMidas)
                          ? 0
                          : data.max_low_lag;
    if (p_ar > 0 && p_ar + 2 >= data.y.size()) {
        throw std::invalid_argument("midas_fit: AR lags too large");
    }

    // ---- 设计矩阵与 y 对齐 (j0 = max(⌊kmax/m⌋+1, p_ar+1)) ----
    const Size n = data.y.size();
    const Size j0_x = (k_high - 1) / data.m + 1;
    const Size j0 = (p_ar > 0) ? std::max(j0_x, p_ar + 1) : j0_x;
    const Size n_eff = n - j0 + 1;
    if (n_eff < k_high + p_ar + 2) {
        throw std::invalid_argument("midas_fit: insufficient observations");
    }
    std::vector<std::vector<Real>> X(n_eff, std::vector<Real>(k_high));
    for (Size r = 0; r < n_eff; ++r) {
        const Size j = j0 + r;  // 1-based 低频期
        for (Size c = 0; c < k_high; ++c) {
            X[r][c] = data.x[data.m * j - c - 1];  // x[m·j − h], h=c (期末系)
        }
    }
    std::vector<Real> y_al(n_eff);
    for (Size r = 0; r < n_eff; ++r) y_al[r] = data.y[j0 + r - 1];
    // AR 列 (行 r, 滞后 ℓ): y[j−1−ℓ] 0-based = data.y[j0 + r − 1 − ℓ]
    std::vector<std::vector<Real>> ar_cols(
        p_ar, std::vector<Real>(n_eff, 0.0));
    for (Size r = 0; r < n_eff; ++r) {
        for (Size l = 1; l <= p_ar; ++l) {
            ar_cols[l - 1][r] = data.y[j0 + r - 1 - l];
        }
    }

    // ---- 内层: 给定 (θ, φ) 构造 [1, Z(θ), AR 部分] → OLS → SSR ----
    // AR 普通项: 每滞后一列 (内层系数); ARStar: 单列 c(φ) 系数固定 1
    // U-MIDAS 忽略 weight/starts (纯 OLS, dth 强制 0)
    const Size dth = (type == MidasType::UMidas)
                         ? 0
                         : detail::theta_dim(weight, starts);
    const bool has_phi = (type == MidasType::ARStar);
    const Size n_outer = dth + (has_phi ? 1 : 0);

    auto inner_ols = [&](const std::vector<Real>& theta, Real phi) {
        std::vector<std::vector<Real>> Z =
            (type == MidasType::UMidas)
                ? std::vector<std::vector<Real>>{}
                : detail::weight_basis_columns(weight, k_high, theta, X);
        // 组装设计阵 [1, Z..., X (UMidas), AR...]
        std::vector<std::vector<Real>> D;
        std::vector<Size> kind;  // 0=截距, 1=权重列(δ 系数), 2=UMidas 系数,
                                 // 3=AR 列(ρ 系数), 4=ARStar c(φ) 固定 1
        D.push_back(std::vector<Real>(n_eff, 1.0));
        kind.push_back(0);
        if (type == MidasType::UMidas) {
            for (Size c = 0; c < k_high; ++c) {
                D.push_back(std::vector<Real>(n_eff));
                for (Size r = 0; r < n_eff; ++r) D.back()[r] = X[r][c];
                kind.push_back(2);
            }
        } else {
            for (auto& zcol : Z) {
                D.push_back(std::move(zcol));
                kind.push_back(1);
            }
        }
        if (type == MidasType::AR) {
            for (Size l = 1; l <= p_ar; ++l) {
                D.push_back(ar_cols[l - 1]);
                kind.push_back(3);
            }
        } else if (type == MidasType::ARStar) {
            std::vector<Real> cph(n_eff, 0.0);
            for (Size r = 0; r < n_eff; ++r) {
                Real s = 0.0, ph = 1.0;
                for (Size l = 1; l <= p_ar; ++l) {
                    ph = (l == 1) ? 1.0 : ph * phi;
                    s += (1.0 - phi) * ph * ar_cols[l - 1][r];
                }
                cph[r] = s;
            }
            D.push_back(std::move(cph));
            kind.push_back(4);
        }
        // OLS (系数 + 残差), 固定列 (kind 4) 不进回归 → 从 y 中减去
        std::vector<Real> y_work = y_al;
        for (Size r = 0; r < n_eff; ++r) {
            for (Size k = 0; k < D.size(); ++k) {
                if (kind[k] == 4) y_work[r] -= D[k][r];
            }
        }
        std::vector<std::vector<Real>> Xreg;
        std::vector<Size> reg_kind;
        for (Size k = 0; k < D.size(); ++k) {
            if (kind[k] != 4) {
                Xreg.push_back(D[k]);
                reg_kind.push_back(kind[k]);
            }
        }
        // ols_fit 约定行主 X[obs][param] — D/Xreg 为列主, 转置
        std::vector<std::vector<Real>> Xrows(
            n_eff, std::vector<Real>(Xreg.size(), 0.0));
        for (Size r = 0; r < n_eff; ++r) {
            for (Size j = 0; j < Xreg.size(); ++j) Xrows[r][j] = Xreg[j][r];
        }
        const auto ols = unit_root::detail::ols_fit(y_work, Xrows);
        return std::make_pair(ols, reg_kind);
    };

    // 集中化 SSR(outer), outer = (θ..., φ?)
    auto concentrated_ssr = [&](const std::vector<Real>& outer) {
        const std::vector<Real> theta(outer.begin(),
                                      outer.begin() + static_cast<std::ptrdiff_t>(dth));
        const Real phi = has_phi ? outer.back() : 0.0;
        // ARStar 平凡陷阱: φ = 0 → ρ 全零 → c 列消失 → 退化为 DL
        auto pr = inner_ols(theta, phi);
        return pr.first.ssr;
    };

    // ---- U-MIDAS / 线性族: 单次 OLS, 无外层 ----
    MidasResult res;
    res.type = type;
    res.n_obs = n_eff;
    if (n_outer == 0) {
        auto pr = inner_ols({}, 0.0);
        const auto& ols = pr.first;
        const auto& reg_kind = pr.second;
        res.intercept = ols.beta[0];
        res.converged = true;
        res.message = "linear (concentrated OLS, no outer parameters)";
        Size ki = 1;
        if (type == MidasType::UMidas) {
            res.delta.assign(ols.beta.begin() + 1,
                             ols.beta.begin() + 1 + static_cast<std::ptrdiff_t>(k_high));
            res.midas_coef = res.delta;
            ki += k_high;
        } else {
            // AlmonP: delta = 权重参数 (K 个); HarStep: (p1,p2,p3)
            for (; ki < ols.beta.size(); ++ki) {
                if (reg_kind[ki] == 1) res.delta.push_back(ols.beta[ki]);
            }
            // 隐含权重: Almonp/HarStep 由参数重建
            if (weight == MidasWeight::AlmonP) {
                std::vector<Real> lam(res.delta.size(), 0.0);
                for (Size kk = 0; kk < res.delta.size(); ++kk) lam[kk] = res.delta[kk];
                res.midas_coef = almonp_weights(lam, k_high);
            } else if (weight == MidasWeight::HarStep) {
                res.midas_coef = harstep_weights(res.delta, k_high);
            }
        }
        // AR 系数 (kind 3 列, 权重列之后; 线性族分支)
        if (type == MidasType::AR) {
            for (Size k = 1; k < ols.beta.size(); ++k) {
                if (reg_kind[k] == 3) res.delta.push_back(ols.beta[k]);
            }
        }
        res.ssr = ols.ssr;
        res.residuals = ols.resid;
        const Size k_total = ols.n_params + (type == MidasType::ARStar ? 1 : 0);
        res.sigma2 = ols.ssr / static_cast<Real>(n_eff - k_total);
        constexpr Real kTwoPi = 6.283185307179586476925286766559;
        const Real s2mle = ols.ssr / static_cast<Real>(n_eff);
        res.loglik = -0.5 * static_cast<Real>(n_eff) *
                     (std::log(kTwoPi * s2mle) + 1.0);
        // fitted
        res.fitted.resize(n_eff);
        for (Size r = 0; r < n_eff; ++r) res.fitted[r] = y_al[r] - ols.resid[r];
        return res;
    }

    // ---- 非线性族: SLSQP 多起点 ----
    std::vector<std::vector<Real>> all_starts =
        detail::default_starts(weight, starts);
    // ARStar: 起点追加 φ ∈ {0, 0.5, −0.5}
    if (has_phi) {
        std::vector<std::vector<Real>> ext;
        for (const auto& s : all_starts) {
            for (Real phi : {0.0, 0.5, -0.5}) {
                std::vector<Real> v = s;
                v.push_back(phi);
                ext.push_back(v);
            }
        }
        all_starts = ext;
    }
    std::vector<Bounds> all_bnds = detail::theta_bounds(weight, dth);
    if (has_phi) all_bnds.push_back(Bounds{-0.98, 0.98});

    SLSQP::Config cfg;
    cfg.max_iterations = 500;
    cfg.ftol = 1e-12;
    cfg.xtol = 1e-10;

    Real best_f = std::numeric_limits<Real>::infinity();
    std::vector<Real> best_x;
    bool conv = false;
    std::string msg;
    for (const auto& x0 : all_starts) {
        if (x0.size() != n_outer) continue;
        OptimizationResult r = SLSQP::minimize(concentrated_ssr, x0, all_bnds,
                                               {}, {}, cfg);
        if (r.fx < best_f) {
            best_f = r.fx;
            best_x = r.x;
            conv = r.converged;
            msg = r.message;
        }
    }
    if (best_x.empty()) {
        throw std::runtime_error("midas_fit: no valid start (dim mismatch)");
    }

    auto pr = inner_ols(std::vector<Real>(best_x.begin(),
                                          best_x.begin() + static_cast<std::ptrdiff_t>(dth)),
                        has_phi ? best_x.back() : 0.0);
    const auto& ols = pr.first;
    const auto& reg_kind = pr.second;
    res.intercept = ols.beta[0];
    res.converged = conv;
    res.message = msg;
    // lambda = θ (Nealmon λ₂.. / NBeta (κ₁,κ₂,θ₀)) + (φ if ARStar)
    res.lambda.assign(best_x.begin(),
                      best_x.begin() + static_cast<std::ptrdiff_t>(dth));
    if (has_phi) res.lambda.push_back(best_x.back());
    // delta: δ (权重单列族, kind 1) + AR ρ (kind 3)
    for (Size k = 1; k < ols.beta.size(); ++k) {
        if (reg_kind[k] == 1 || reg_kind[k] == 3) {
            res.delta.push_back(ols.beta[k]);
        }
    }
    // midas_coef = δ·w̄(θ̂)
    {
        std::vector<Real> theta_hat(res.lambda.begin(),
                                    res.lambda.begin() + static_cast<std::ptrdiff_t>(dth));
        std::vector<Real> wbar;
        if (weight == MidasWeight::Nealmon) {
            std::vector<Real> lam{1.0};
            lam.insert(lam.end(), theta_hat.begin(), theta_hat.end());
            wbar = nealmon_weights(lam, k_high);
        } else {  // NBeta
            wbar = nbeta_weights({1.0, theta_hat[0], theta_hat[1], theta_hat[2]},
                                 k_high);
        }
        res.midas_coef.resize(k_high);
        for (Size i = 0; i < k_high; ++i) {
            res.midas_coef[i] = res.delta.empty() ? wbar[i] : res.delta[0] * wbar[i];
        }
    }
    res.ssr = ols.ssr;
    res.residuals = ols.resid;
    const Size k_total = ols.n_params + (has_phi ? 1 : 0);
    res.sigma2 = ols.ssr / static_cast<Real>(n_eff - k_total);
    constexpr Real kTwoPi = 6.283185307179586476925286766559;
    const Real s2mle = ols.ssr / static_cast<Real>(n_eff);
    res.loglik = -0.5 * static_cast<Real>(n_eff) * (std::log(kTwoPi * s2mle) + 1.0);
    res.fitted.resize(n_eff);
    for (Size r = 0; r < n_eff; ++r) res.fitted[r] = y_al[r] - ols.resid[r];
    return res;
}

}  // namespace midas
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
