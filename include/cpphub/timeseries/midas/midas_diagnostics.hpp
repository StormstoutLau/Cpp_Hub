// =============================================================================
// midas_diagnostics.hpp - MIDAS 诊断: hAh 权重检验 + 残差 LB/JB (spec §6.3)
//
// Phase 7C v1.7 M4 (PHASE7C_SPEC.md v1.2 §6.3, MD6)
//
// hAh restriction test (K-Z 2012; midasr hAh_test/prep_hAh 源码一手实录
// 2026-08-18):
//   H0: 高频滞后系数约束于权重参数化 (Nealmon/NBeta/AlmonP/HarStep)
//   Ha: 无约束 (U-MIDAS 系数)
//   X = [1, mls 列, AR 列] (n_eff × dk);  cfur = 无约束 OLS 系数
//   P = chol(X'X) (upper, P'P = X'X)
//   h.0 = P·(cfur − midas_coef_all)                  (dk×1)
//   Delta.0 = D0·(D0'X'X·D0)^{-1}·D0',  D0 = ∂midas_coef/∂params
//             (中心差分, 与 midasr 默认 numDeriv::jacobian 路径一致)
//   se2 = SSR_ur/(n_eff − dk)
//   stat = h.0'·[(I − P·Delta.0·P')/se2]·h.0
//   df = dk − nparam;  p = 1 − χ²cdf(stat, df)
//
// 对照基准 (verify_midas.R W7, 2026-08-18 新 DGP):
//   DL nealmon (k_high=4, 参数 3): stat=2.1080208, p=0.3485372, df=2
//
// 残差诊断: ljung_box_test / jarque_bera_test 复用 (Phase 7A, ADR-015)
// =============================================================================

#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cpphub/core/math.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/inference/residual_diagnostics.hpp"
#include "cpphub/timeseries/midas/midas_model.hpp"
#include "cpphub/timeseries/midas/midas_weights.hpp"
#include "cpphub/timeseries/midas/mixed_freq_data.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace midas {

/// hAh 检验结果 (三列全录, MD6)
struct HahResult {
    Real stat = 0.0;
    Real p = std::numeric_limits<Real>::quiet_NaN();
    Size df = 0;
};

namespace detail {

/// 上三角 Cholesky: A = U'U (A 对称正定)
inline std::vector<std::vector<Real>> chol_upper(
    const std::vector<std::vector<Real>>& A) {
    const Size n = A.size();
    std::vector<std::vector<Real>> U(n, std::vector<Real>(n, 0.0));
    for (Size i = 0; i < n; ++i) {
        for (Size j = i; j < n; ++j) {
            Real s = A[i][j];
            for (Size k = 0; k < i; ++k) s -= U[k][i] * U[k][j];
            if (i == j) {
                if (s <= 0.0) throw std::runtime_error("chol: not PD");
                U[i][i] = std::sqrt(s);
            } else {
                U[i][j] = s / U[i][i];
            }
        }
    }
    return U;
}

inline std::vector<std::vector<Real>> mat_mul(
    const std::vector<std::vector<Real>>& A,
    const std::vector<std::vector<Real>>& B) {
    const Size n = A.size(), m = B[0].size(), kc = B.size();
    std::vector<std::vector<Real>> C(n, std::vector<Real>(m, 0.0));
    for (Size i = 0; i < n; ++i)
        for (Size t = 0; t < kc; ++t) {
            const Real a = A[i][t];
            if (a == 0.0) continue;
            for (Size j = 0; j < m; ++j) C[i][j] += a * B[t][j];
        }
    return C;
}

inline std::vector<std::vector<Real>> mat_t(
    const std::vector<std::vector<Real>>& A) {
    const Size n = A.size(), m = A[0].size();
    std::vector<std::vector<Real>> T(m, std::vector<Real>(n));
    for (Size i = 0; i < n; ++i)
        for (Size j = 0; j < m; ++j) T[j][i] = A[i][j];
    return T;
}

inline std::vector<Real> mat_vec(const std::vector<std::vector<Real>>& A,
                                 const std::vector<Real>& v) {
    std::vector<Real> out(A.size(), 0.0);
    for (Size i = 0; i < A.size(); ++i) {
        Real s = 0.0;
        for (Size j = 0; j < v.size(); ++j) s += A[i][j] * v[j];
        out[i] = s;
    }
    return out;
}

}  // namespace detail

/// @brief hAh 权重约束检验 (K-Z 2012) — 对已估 MIDAS 模型
/// @param fit midas_fit 输出 (DL/AR/ARStar; U-MIDAS 本身即无约束备择, 拒绝)
/// @param data 原始混频数据 (与 fit 一致)
/// @param k_high 高频滞后数 (与 fit 一致)
/// @param weight 与 midas_fit 相同的权重族
/// @throws std::invalid_argument U-MIDAS / PolyStep / 维度不符
inline HahResult hAh_test(const MidasResult& fit, const MixedFreqData& data,
                          Size k_high, MidasWeight weight) {
    if (fit.type == MidasType::UMidas) {
        throw std::invalid_argument(
            "hAh_test: U-MIDAS is the unrestricted alternative itself");
    }
    if (weight == MidasWeight::PolyStep) {
        throw std::invalid_argument("hAh_test: PolyStep unsupported");
    }
    data.validate();
    const Size p_ar = (fit.type == MidasType::DL) ? 0 : data.max_low_lag;
    const Size n = data.y.size();
    const Size j0_x = (k_high - 1) / data.m + 1;
    const Size j0 = (p_ar > 0) ? std::max(j0_x, p_ar + 1) : j0_x;
    const Size n_eff = n - j0 + 1;

    // ---- 无约束设计阵 X = [1, mls 列, AR 列] (n_eff × dk) ----
    const Size dk = 1 + k_high + p_ar;
    std::vector<std::vector<Real>> X(n_eff, std::vector<Real>(dk, 0.0));
    for (Size r = 0; r < n_eff; ++r) {
        const Size j = j0 + r;
        X[r][0] = 1.0;
        for (Size c = 0; c < k_high; ++c) {
            X[r][1 + c] = data.x[data.m * j - c - 1];
        }
        for (Size l = 1; l <= p_ar; ++l) {
            X[r][1 + k_high + l - 1] = data.y[j0 + r - 1 - l];
        }
    }
    std::vector<Real> y_al(n_eff);
    for (Size r = 0; r < n_eff; ++r) y_al[r] = data.y[j0 + r - 1];
    const auto ur = unit_root::detail::ols_fit(y_al, X);

    // ---- 参数向量 prm = (μ, [权重参数块], θ..., ρ.../φ) ----
    // Nealmon/NBeta: 权重块 = (δ), θ = 外层形状; AlmonP/HarStep: 权重块 =
    //   (p₁..p_K) (线性), θ 空; ARStar 末位 φ, AR 追加 ρ
    const bool is_arstar = fit.type == MidasType::ARStar;
    const Size n_theta = fit.lambda.size() - (is_arstar ? 1 : 0);
    Size n_weight_params = 1;
    if (weight == MidasWeight::AlmonP) {
        n_weight_params = fit.delta.empty() ? 2 : fit.delta.size();
    } else if (weight == MidasWeight::HarStep) {
        n_weight_params = 3;
    }
    const Size nparam =
        1 + n_weight_params + n_theta + (is_arstar ? 1 : p_ar);

    std::vector<Real> prm0;
    prm0.push_back(fit.intercept);
    if (weight == MidasWeight::AlmonP || weight == MidasWeight::HarStep) {
        for (Size k = 0; k < n_weight_params; ++k) {
            prm0.push_back(k < fit.delta.size() ? fit.delta[k] : 0.0);
        }
    } else {
        prm0.push_back(fit.delta.empty() ? 0.0 : fit.delta[0]);
    }
    for (Size k = 0; k < n_theta; ++k) prm0.push_back(fit.lambda[k]);
    if (is_arstar) {
        prm0.push_back(fit.lambda.back());
    } else {
        for (Size l = 1; l <= p_ar; ++l) {
            const Size idx = n_weight_params + l - 1;
            prm0.push_back(idx < fit.delta.size() ? fit.delta[idx] : 0.0);
        }
    }

    // ---- all_coef(prm) → (μ, w₁..w_{k_high}, ρ...) 长度 dk ----
    auto all_coef = [&](const std::vector<Real>& prm) {
        std::vector<Real> coef(dk);
        coef[0] = prm[0];
        const Real dlt = prm[1];
        const std::vector<Real> theta_in(
            prm.begin() + 1 + static_cast<std::ptrdiff_t>(n_weight_params),
            prm.begin() + 1 + static_cast<std::ptrdiff_t>(n_weight_params + n_theta));
        const std::vector<Real> lin_prm(
            prm.begin() + 1,
            prm.begin() + 1 + static_cast<std::ptrdiff_t>(n_weight_params));
        switch (weight) {
        case MidasWeight::Nealmon: {
            std::vector<Real> lam{1.0};
            lam.insert(lam.end(), theta_in.begin(), theta_in.end());
            const auto wb = nealmon_weights(lam, k_high);
            for (Size i = 0; i < k_high; ++i) coef[1 + i] = dlt * wb[i];
            break;
        }
        case MidasWeight::NBeta: {
            const auto wb = nbeta_weights({1.0, theta_in[0], theta_in[1],
                                           theta_in[2]},
                                          k_high);
            for (Size i = 0; i < k_high; ++i) coef[1 + i] = dlt * wb[i];
            break;
        }
        case MidasWeight::AlmonP: {
            const auto w = almonp_weights(lin_prm, k_high);
            for (Size i = 0; i < k_high; ++i) coef[1 + i] = w[i];
            break;
        }
        case MidasWeight::HarStep: {
            const auto w = harstep_weights(lin_prm, k_high);
            for (Size i = 0; i < k_high; ++i) coef[1 + i] = w[i];
            break;
        }
        case MidasWeight::PolyStep:
            break;  // 前置拒绝
        }
        if (is_arstar) {
            const Real phi = prm.back();
            Real ph = 1.0;
            for (Size l = 1; l <= p_ar; ++l) {
                ph = (l == 1) ? 1.0 : ph * phi;
                coef[1 + k_high + l - 1] = (1.0 - phi) * ph;
            }
        } else {
            for (Size l = 1; l <= p_ar; ++l) {
                coef[1 + k_high + l - 1] =
                    prm[1 + n_weight_params + n_theta + l - 1];
            }
        }
        return coef;
    };

    // ---- D0 (dk × nparam): 中心差分 ----
    const std::vector<Real> c0 = all_coef(prm0);
    std::vector<std::vector<Real>> D0(dk, std::vector<Real>(nparam, 0.0));
    for (Size k = 0; k < nparam; ++k) {
        const Real hs = 1e-6 * std::max(1.0, std::fabs(prm0[k]));
        std::vector<Real> pp = prm0, pm = prm0;
        pp[k] += hs;
        pm[k] -= hs;
        const auto cp = all_coef(pp);
        const auto cm = all_coef(pm);
        for (Size i = 0; i < dk; ++i) D0[i][k] = (cp[i] - cm[i]) / (2.0 * hs);
    }

    // ---- XtX (对称), P = chol, h.0 ----
    std::vector<std::vector<Real>> XtX(dk, std::vector<Real>(dk, 0.0));
    for (Size i = 0; i < dk; ++i)
        for (Size jx = i; jx < dk; ++jx) {
            Real s = 0.0;
            for (Size r = 0; r < n_eff; ++r) s += X[r][i] * X[r][jx];
            XtX[i][jx] = s;
            XtX[jx][i] = s;
        }
    const auto P = detail::chol_upper(XtX);
    std::vector<Real> diff(dk);
    for (Size i = 0; i < dk; ++i) diff[i] = ur.beta[i] - c0[i];
    const std::vector<Real> h0 = detail::mat_vec(P, diff);

    // ---- Delta.0 = D0 (D0'XtX D0)^{-1} D0' ----
    const auto D0t = detail::mat_t(D0);
    const auto Mq = detail::mat_mul(detail::mat_mul(D0t, XtX), D0);
    auto Minv = Mq;
    unit_root::detail::invert_matrix(Minv);
    const auto Delta0 = detail::mat_mul(detail::mat_mul(D0, Minv), D0t);

    // ---- stat = h0'(I − P·Delta0·P')h0 / se2 ----
    const auto PDPt = detail::mat_mul(detail::mat_mul(P, Delta0),
                                      detail::mat_t(P));
    const Real se2 = ur.ssr / static_cast<Real>(n_eff - dk);
    const std::vector<Real> vh = detail::mat_vec(PDPt, h0);
    Real s2 = 0.0;
    for (Size i = 0; i < dk; ++i) s2 += h0[i] * (h0[i] - vh[i]);

    HahResult out;
    out.stat = s2 / se2;
    out.df = dk - nparam;
    if (out.df == 0) {
        out.p = std::numeric_limits<Real>::quiet_NaN();
    } else {
        out.p = 1.0 - regularized_lower_gamma(
                           static_cast<Real>(out.df) / 2.0, out.stat / 2.0);
    }
    return out;
}

/// @brief MIDAS 残差诊断 (LB 自相关 + JB 正态, Phase 7A 复用)
struct MidasResidualDiagnostics {
    econometrics::LjungBoxResult ljung_box;
    econometrics::JarqueBeraResult jarque_bera;
    bool residual_white = false;   ///< LB 不拒绝
    bool residual_normal = false;  ///< JB 不拒绝
};

inline MidasResidualDiagnostics midas_residual_diagnostics(
    const MidasResult& fit, Size lb_lag = 10) {
    MidasResidualDiagnostics d;
    d.ljung_box = econometrics::ljung_box_test(fit.residuals, lb_lag);
    d.jarque_bera = econometrics::jarque_bera_test(fit.residuals);
    d.residual_white = !d.ljung_box.base.reject_null;
    d.residual_normal = !d.jarque_bera.base.reject_null;
    return d;
}

}  // namespace midas
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
