// =============================================================================
// vecm_model.hpp - VECM 估计 (Johansen ML, 5 确定性情形 + β 双归一 + ECT t)
//
// Phase 7C v1.7 M3 (PHASE7C_SPEC.md §5.4; 决策 17/21)
//
// 教材锚点: Johansen 1988 JEDC / Lütkepohl 2005 Ch.7 (pp.286-299);
//   Ericsson-MacKinnon 2002 Econometrics J. 5(2):285-318 (ECT t 临界值)
//
// API 等价 statsmodels VECM(endog, k_ar_diff, coint_rank, deterministic).fit():
//   det ∈ {n, co, ci, lo, li} (CI4: 5 情形; 纯情形, 不支持 "colo" 等拼接)
//   alpha/beta/gamma/det_coef/det_coef_coint/sigma_u/llf 1e-10 对齐
//
// 算法 (复刻 statsmodels vecm.py L1000-1052):
//   p = k+1; y_lag1 = y[p−1:T−1] + [1 (ci)] + [arange(T)+p (li)]
//   ΔX = Δy 滞后 (行序 [Δy_{t−1} v0..vK, Δy_{t−2} v0..vK, ...]) +
//        [1 (co)] + [arange(T)+p+1 (lo)]
//   R0/R1 = Δy, y_lag1 对 ΔX 的 OLS 残差 (≡ SM 的 m 投影阵, 免 T×T 矩阵)
//   λ = eig(s11_·s10·s00⁻¹·s01·s11_) 降序; β̃ = s11_·V_r → β̃·inv(β̃[:r])
//     (前 r 行 = I_r, "Phillips 归一" 文献惯称; 列符号在归一中消去 → 确定性)
//   α = s01·β̃·inv(β̃'s11β̃);  Γ = (Δy−αβ̃'y_lag1)·ΔX'(ΔXΔX')⁻¹;
//   Σ = resid'resid/T;  llf (Lütkepohl 7.2.20)
//
// β 双归一 (决策 21, CI8):
//   默认: 前 r 行 = I_r (statsmodels); 开关 urca_normalization: 首变量 = 1
//   (β_new = β·diag(1/β[0,j]), α_new = α·diag(β[0,j]), Π = αβ' 不变)
//   ⚠️ r=1 时两种归一数值相同; 对照测试恒用投影矩阵 P=β(β'β)⁻¹β' (非逐元素)
//
// ECT t 检验 (CI10): t = α_j/SE, SE(α[j,i]) = sqrt(Σ_jj·mat1[i,i]),
//   mat1 = b_id·ω⁻¹·b_id' (statsmodels cov_params_default, Lütkepohl 7.2.21);
//   临界值用 EM2002 响应面 (非标准分布, ericsson_mackinnon_cv.hpp)。
//   仅 rank=1 计算 (逐方程单 ECT 情形; r>1 → NaN + has_ect_t=false, §1.4-5);
//   det→EM2002 case 映射: n→n, co/ci→c, lo/li→ct (近似, 文档声明)
//
// 排幻觉点 (spec §9.3): CI4 (5 情形) / CI8 (β 投影空间对照) /
//   CI9 (Π=αβ'=ΣAᵢ−I, αᵢ<0 拉回) / CI10 (EM2002 表非标准 t)
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/cointegration/ericsson_mackinnon_cv.hpp"
#include "cpphub/timeseries/var/multivariate_data.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace cointegration {

/// VECM 估计结果 (§5.4 v1.2 签名 + statsmodels 镜像附加字段)
struct VECMResult {
    Eigen::MatrixXd alpha;               ///< K×r 载入矩阵
    Eigen::MatrixXd beta;                ///< K×r; 默认前 r 行 = I_r (决策 21)
    std::vector<Eigen::MatrixXd> gamma;  ///< 逐滞后 Γᵢ (K×K, lag 1..k)
    Eigen::MatrixXd gamma_flat;   ///< K×(K·k) statsmodels gamma 原序 (附加)
    Eigen::MatrixXd det_coef;     ///< K×#det_out: co 常数列在前, lo 趋势列在后 (附加)
    Eigen::MatrixXd det_coef_coint;  ///< #det_coint×r: ci 常数行在前, li 趋势行在后 (附加)
    Eigen::MatrixXd sigma_u;      ///< K×K (ML: ÷T) (附加)
    Eigen::MatrixXd resid;        ///< K×T_eff (附加)
    Real loglik = 0.0;
    Size rank = 0;
    std::string det;             ///< 5 情形 {n, co, ci, lo, li} (CI4)
    Size k_ar_diff = 0;
    Size n_obs = 0;              ///< T_eff = T − k − 1
    std::vector<Real> ect_t_stat;   ///< K 个 (rank=1); 否则 NaN
    std::vector<Real> ect_cv_5pct;  ///< EM2002 5% 临界值 (rank=1); 否则 NaN
    bool has_ect_t = false;         ///< rank=1 才计算 (§1.4-5)
    bool urca_normalization = false;  ///< 回显
};

namespace detail {

// EM2002 case 映射: n→n; co/ci→c; lo/li→ct (近似)
inline std::string vecm_det_to_em2002(const std::string& det) {
    if (det == "n") return "n";
    if (det == "co" || det == "ci") return "c";
    return "ct";  // lo / li
}

// 对称正定矩阵的 -1/2 次幂 (SM _mat_sqrt 的逆; SVD≡特征分解, 对称情形等价)
inline Eigen::MatrixXd mat_inv_sqrt(const Eigen::MatrixXd& s) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        s, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::VectorXd d = svd.singularValues();
    for (Eigen::Index i = 0; i < d.size(); ++i) {
        if (d(i) <= 0.0) {
            throw std::runtime_error("vecm_fit: s11 not positive definite");
        }
        d(i) = std::sqrt(d(i));
    }
    // s11_ = inv(U·√Σ·V') — SVD 路径逐字复刻 SM _mat_sqrt+inv
    // (eigh 路径的 s11_ 与 SM 差 ~7e-11, 经 λ 放大破坏 llf 的 1e-10 对齐)
    return (svd.matrixU() * d.asDiagonal() * svd.matrixV().transpose())
        .inverse();
}

}  // namespace detail

// ---------------------------------------------------------------------------
// VECM 估计 (Johansen ML)
//
// @param data K 列水平数据 (K ≥ 2)
// @param rank 协整秩 r ∈ [1, K]
// @param k_ar_diff 滞后差分阶数 (≥ 0)
// @param det "n"/"co"/"ci"/"lo"/"li" (CI4: 5 情形; 常数+趋势组合 v1.8)
// @param urca_normalization false: β 前 r 行 = I_r (statsmodels);
//   true: 首变量 = 1 (urca 式; r=1 时两者相同)
// ---------------------------------------------------------------------------
inline VECMResult vecm_fit(const var::MultivariateTSData& data, Size rank,
                           Size k_ar_diff, const std::string& det = "n",
                           bool urca_normalization = false) {
    const bool det_ok = det == "n" || det == "co" || det == "ci" ||
                        det == "lo" || det == "li";
    if (!det_ok) {
        throw std::invalid_argument(
            "vecm_fit: det must be one of n/co/ci/lo/li");
    }
    data.validate();
    const Size k = data.K();
    if (rank < 1 || rank > k) {
        throw std::invalid_argument("vecm_fit: rank must be in [1, K]");
    }
    const Size t_tot = data.T();
    const Size p = k_ar_diff + 1;
    if (t_tot < p + 2) {
        throw std::invalid_argument(
            "vecm_fit: sample too small for k_ar_diff");
    }

    const Size t_eff = t_tot - p;
    const Eigen::MatrixXd y = data.matrix();  // T×K
    const Eigen::Index kk = static_cast<Eigen::Index>(k);
    const Eigen::Index tt = static_cast<Eigen::Index>(t_eff);

    // Δy (T−1)×K; dy1t (T_eff×K): Δy_t, t = p..T−1
    Eigen::MatrixXd dy(t_tot - 1, kk);
    for (Size i = 0; i < t_tot - 1; ++i) dy.row(i) = y.row(i + 1) - y.row(i);
    Eigen::MatrixXd dy1t = dy.bottomRows(tt);

    // y_lag1 ((K+det_c)×T_eff): y_{t−1} + [1 (ci)] + [arange(T_eff)+p (li)]
    const bool has_ci = (det == "ci");
    const bool has_li = (det == "li");
    const bool has_co = (det == "co");
    const bool has_lo = (det == "lo");
    const Size n_det_c = (has_ci ? 1u : 0u) + (has_li ? 1u : 0u);
    Eigen::MatrixXd y_lag1(kk + static_cast<Eigen::Index>(n_det_c), tt);
    for (Size i = 0; i < t_eff; ++i) {
        y_lag1.col(static_cast<Eigen::Index>(i)).head(kk) =
            y.row(p - 1 + i).transpose();
    }
    {
        Eigen::Index row = kk;
        if (has_ci) {
            for (Size i = 0; i < t_eff; ++i) y_lag1(row, i) = 1.0;
            ++row;
        }
        if (has_li) {
            for (Size i = 0; i < t_eff; ++i) {
                y_lag1(row, i) = static_cast<Real>(i) + static_cast<Real>(p);
            }
        }
    }

    // delta_x ((K·k + det_out)×T_eff): 行序 [Δy_{t−1} v0..vK, ...lag2...] +
    //   [1 (co)] + [arange(T_eff)+p+1 (lo)]
    const Size n_det_o = (has_co ? 1u : 0u) + (has_lo ? 1u : 0u);
    const Eigen::Index nx =
        static_cast<Eigen::Index>(k * k_ar_diff + n_det_o);
    Eigen::MatrixXd dx(nx, tt);
    for (Size i = 0; i < t_eff; ++i) {
        const Size t_reg = p + i;  // 回归时刻 t (0-based)
        for (Size j = 0; j < k_ar_diff; ++j) {
            for (Size v = 0; v < k; ++v) {
                // Δy_{t−1−j} = dy 行 (t−1−j)−1 = t_reg−2−j
                dx(static_cast<Eigen::Index>(j * k + v),
                   static_cast<Eigen::Index>(i)) = dy(t_reg - 2 - j, v);
            }
        }
    }
    {
        Eigen::Index row = static_cast<Eigen::Index>(k * k_ar_diff);
        if (has_co) {
            for (Size i = 0; i < t_eff; ++i) dx(row, i) = 1.0;
            ++row;
        }
        if (has_lo) {
            for (Size i = 0; i < t_eff; ++i) {
                dx(row, i) = static_cast<Real>(i) + static_cast<Real>(p) + 1.0;
            }
        }
    }

    // R0/R1: 对 ΔX 的 OLS 残差 (≡ SM m 投影阵: I − ΔX'(ΔXΔX')⁻¹ΔX;
    // K×T_eff / (K+det_c)×T_eff; 设计矩阵 ΔX' 为高矩阵, 正规方程与 SM 逐字一致)
    const Eigen::MatrixXd dxdx_inv = (dx * dx.transpose()).inverse();
    const Eigen::MatrixXd r0 =
        (dy1t - dx.transpose() * dxdx_inv * dx * dy1t).transpose();
    const Eigen::MatrixXd y_lag1_t = y_lag1.transpose();  // T_eff×(K+det_c)
    const Eigen::MatrixXd r1 =
        (y_lag1_t - dx.transpose() * dxdx_inv * dx * y_lag1_t).transpose();

    const Real tn = static_cast<Real>(t_eff);
    Eigen::MatrixXd s00 = r0 * r0.transpose() / tn;
    Eigen::MatrixXd s01 = r0 * r1.transpose() / tn;
    Eigen::MatrixXd s11 = r1 * r1.transpose() / tn;
    Eigen::MatrixXd s11_ = detail::mat_inv_sqrt(s11);

    // λ = eig(s11_·s10·s00⁻¹·s01·s11_) 降序; V 列与 λ 对应
    // 非对称 eig 路径逐字复刻 SM (vecm.py L451 np.linalg.eig);
    // meig 数学对称但 eigh 路径的 λ 与 SM 差 ~2e-11, 经 T/(2(1−λ₀))≈202×
    // 放大破坏 llf 的 1e-10 对齐 (虚部数值为 ~0, 取实部)
    Eigen::MatrixXd s01_s11_ = s01 * s11_;
    Eigen::MatrixXd meig = s01_s11_.transpose() * s00.inverse() * s01_s11_;
    Eigen::EigenSolver<Eigen::MatrixXd> es(meig);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("vecm_fit: eigen decomposition failed");
    }
    const Size m_dim = static_cast<Size>(meig.rows());
    std::vector<std::pair<Real, Eigen::Index>> order(m_dim);
    for (Size i = 0; i < m_dim; ++i) {
        order[i] = {es.eigenvalues()(static_cast<Eigen::Index>(i)).real(),
                    static_cast<Eigen::Index>(i)};
    }
    std::sort(order.begin(), order.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });
    Eigen::VectorXd lam_all(static_cast<Eigen::Index>(m_dim));
    Eigen::MatrixXd v_all(meig.rows(), meig.cols());
    for (Size i = 0; i < m_dim; ++i) {  // 降序 (SM: argsort[::-1])
        lam_all(static_cast<Eigen::Index>(i)) = order[i].first;
        v_all.col(static_cast<Eigen::Index>(i)) =
            es.eigenvectors().col(order[i].second).real();
    }

    const Eigen::Index rr = static_cast<Eigen::Index>(rank);
    // 未归一 β̃/α̃ (urca 归一基; β̃'s11β̃ = V'V = I ⇒ α̃ = s01·β̃)
    const Eigen::MatrixXd beta_raw = s11_ * v_all.leftCols(rr);  // (K+det_c)×r
    const Eigen::MatrixXd alpha_raw = s01 * beta_raw;
    Eigen::MatrixXd beta_t = beta_raw * beta_raw.topRows(rr).inverse();
    beta_t.topRows(rr) = Eigen::MatrixXd::Identity(rr, rr);  // 消舍入, 精确 I_r

    Eigen::MatrixXd alpha_t =
        s01 * beta_t * (beta_t.transpose() * s11 * beta_t).inverse();

    // Γ (K×(K·k+det_out)) = (Δy − αβ̃'y_lag1)·ΔX'(ΔXΔX')⁻¹
    const Eigen::MatrixXd dy_fit = dy1t.transpose() -
                                   alpha_t * beta_t.transpose() * y_lag1;
    Eigen::MatrixXd gamma_sm =
        dy_fit * dx.transpose() * (dx * dx.transpose()).inverse();

    const Eigen::MatrixXd resid_sm =
        dy_fit - gamma_sm * dx;  // K×T_eff
    Eigen::MatrixXd sigma_u = resid_sm * resid_sm.transpose() / tn;

    // llf (Lütkepohl 7.2.20, SM L1471-1485)
    // log|s00| 走 LU 路径 (= SM np.log(np.linalg.det)); Eigen 3×3 直接
    // determinant() 是余子式展开 (与 LAPACK LU 路径不同), 用 partialPivLu 对齐
    const Real logdet_s00 =
        std::log(s00.partialPivLu().determinant());
    Real sum_log = 0.0;
    for (Size i = 0; i < rank; ++i) {
        sum_log += std::log(1.0 - lam_all(static_cast<Eigen::Index>(i)));
    }
    const Real llf =
        -static_cast<Real>(k) * tn * std::log(2.0 * std::numbers::pi) / 2.0 -
        tn * (logdet_s00 + sum_log) / 2.0 -
        static_cast<Real>(k) * tn / 2.0;

    // ---- α 标准误 (statsmodels cov_params_default, Lütkepohl 7.2.21):
    // mat1 = b_id·ω⁻¹·b_id'; SE(α[j,i])² = Σ_jj·mat1[i,i]
    Eigen::MatrixXd b_y = beta_t.transpose() * y_lag1;  // r×T_eff
    const Eigen::Index nw = rr + nx;
    Eigen::MatrixXd omega(nw, nw);
    omega.block(0, 0, rr, rr) = b_y * b_y.transpose();
    omega.block(0, rr, rr, nx) = b_y * dx.transpose();
    omega.block(rr, 0, nx, rr) = dx * b_y.transpose();
    omega.block(rr, rr, nx, nx) = dx * dx.transpose();
    const Eigen::MatrixXd omega_inv = omega.inverse();

    const Eigen::Index nb = beta_t.rows();  // K + det_c
    Eigen::MatrixXd b_id = Eigen::MatrixXd::Zero(nb + nx, rr + nx);
    b_id.block(0, 0, nb, rr) = beta_t;
    b_id.block(nb, rr, nx, nx) = Eigen::MatrixXd::Identity(nx, nx);
    const Eigen::MatrixXd mat1 = b_id * omega_inv * b_id.transpose();

    // ---- 结果组装 ----
    VECMResult res;
    res.alpha = alpha_t;
    res.beta = beta_t.topRows(kk);
    res.det_coef_coint = beta_t.bottomRows(nb - kk);
    res.gamma_flat =
        gamma_sm.leftCols(static_cast<Eigen::Index>(k * k_ar_diff));
    res.det_coef = gamma_sm.rightCols(static_cast<Eigen::Index>(n_det_o));
    res.gamma.resize(k_ar_diff);
    for (Size j = 0; j < k_ar_diff; ++j) {
        res.gamma[j] = res.gamma_flat.middleCols(
            static_cast<Eigen::Index>(j * k), kk);
    }
    res.sigma_u = sigma_u;
    res.resid = resid_sm;
    res.loglik = llf;
    res.rank = rank;
    res.det = det;
    res.k_ar_diff = k_ar_diff;
    res.n_obs = t_eff;

    // ECT t (仅 rank=1; CI10; 于默认归一下计算)
    res.ect_t_stat.assign(k, std::numeric_limits<Real>::quiet_NaN());
    res.ect_cv_5pct.assign(k, std::numeric_limits<Real>::quiet_NaN());
    if (rank == 1) {
        const Real se_scale = std::sqrt(mat1(0, 0));
        const std::string em_case = detail::vecm_det_to_em2002(det);
        for (Size j = 0; j < k; ++j) {
            const Real se = std::sqrt(sigma_u(j, j)) * se_scale;
            res.ect_t_stat[j] = res.alpha(j, 0) / se;
            res.ect_cv_5pct[j] =
                em2002_ect_critical_value(k, em_case, 5.0, t_eff);
        }
        res.has_ect_t = true;
    }

    // urca 归一开关 (决策 21): 首变量 = 1 (r=1 时与默认相同)。
    // 基于未归一 β̃/α̃ 缩放 (默认归一的 β[0,j] = I_r 结构元, rank≥2 时
    // β[0,1] ≡ 0 不可作除数); Π = α̃β̃' 在两种归一下不变 (CI9)
    res.urca_normalization = urca_normalization;
    if (urca_normalization) {
        Eigen::MatrixXd beta_new = beta_raw;
        Eigen::MatrixXd alpha_new = alpha_raw;
        for (Size j = 0; j < rank; ++j) {
            const Eigen::Index jj = static_cast<Eigen::Index>(j);
            const Real b0 = beta_raw(0, jj);
            if (b0 == 0.0) {
                throw std::runtime_error(
                    "vecm_fit: urca normalization needs beta_raw[0,j] != 0");
            }
            beta_new.col(jj) = beta_raw.col(jj) / b0;
            alpha_new.col(jj) = alpha_raw.col(jj) * b0;
        }
        res.alpha = alpha_new;
        res.beta = beta_new.topRows(kk);
        res.det_coef_coint = beta_new.bottomRows(nb - kk);
    }
    return res;
}

}  // namespace cointegration
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
