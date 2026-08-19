// =============================================================================
// johansen_test.hpp - Johansen 协整秩检验 (迹 / 最大特征值, 3 det 情形)
//
// Phase 7C v1.7 M3 (PHASE7C_SPEC.md §5.2; 决策 19)
//
// 教材锚点: Johansen 1988 JEDC 12:231-254 (⚠️ 非 JASA, CI11);
//   Lütkepohl 2005 p.292 (RRR 算法); OL1992 OBES 54(3):461-472;
//   MHM96 (Queen's DP, johdist.f)
//
// API 等价 statsmodels coint_johansen(endog, det_order, k_ar_diff):
//   det_order ∈ {−1, 0, 1} (CI4: 3 情形; 5 情形归 VECM 类)
//   统计量算法复刻 LeSage/SM 管线 (vecm.py L603-737), 全 det × k 1e-10 对齐;
//   urca 交叉: det_order=0 ↔ ecdet="none", k_ar_diff = K−1 (diff 报告 §2)
//
// 算法 (diff 报告附录 B):
//   detrend(y, det_order): −1 原样 / 0 去均值 / 1 对 [1,t] 去势
//     (vander(linspace(-1,1,T)) 张成空间; 二轮去势阶 f = 0 (det≥0) 或 −1)
//   r0t = resid(Δy_t, Δy 滞后); rkt = resid(y_{t−k}, Δy 滞后)
//     (y_{t−k} ≡ y_{t−1}: 差被 Δy 滞后张成, diff 报告 §4 恒等式)
//   λ = eig(S₁₁⁻¹·S₁₀S₀₀⁻¹S₀₁) 降序 (CI6); T_eff = T−1−k
//   lr1[r] = −T·Σ_{i>r} ln(1−λᵢ); lr2[r] = −T·ln(1−λ_{r+1})
//   evec: β'S₁₁β = I 归一 (Cholesky 路径) + 首非零元符号约定
//     (CI8: 特征向量列符号依赖特征分解实现, 逐元素对照需列符号对齐)
//   cvt/cvm: MHM96 表 (c_sjt/c_sja 同源, 90/95/99); OL1992 独立查表见
//     osterwald_lenum_cv.hpp (双表 diff 报告 §5/§6 裁决)
//
// 排幻觉点 (spec §9.3): CI4 (3 vs 5 情形) / CI5 (双表源) / CI6 (迹公式+
//   λ 降序+有效 T) / CI7 (spec transitory/longrun 数学恒等, diff 报告 §3)
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/cointegration/osterwald_lenum_cv.hpp"
#include "cpphub/timeseries/var/multivariate_data.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace cointegration {

/// Johansen 检验结果 (§5.2 v1.2 签名; 属性名镜像 statsmodels)
struct JohansenResult {
    Eigen::VectorXd eig;    ///< λ̂ 降序 (CI6)
    Eigen::MatrixXd evec;   ///< β'S₁₁β = I 归一特征向量 (CI8: 列符号任意)
    Eigen::VectorXd lr1;    ///< 迹统计量, r = 0..N−1
    Eigen::VectorXd lr2;    ///< 最大特征值统计量, r = 0..N−1
    Eigen::MatrixXd cvt;    ///< N×3 (90/95/99), MHM96 (cv_source 回显)
    Eigen::MatrixXd cvm;    ///< N×3 (90/95/99), MHM96
    int det_order = 0;
    Size k_ar_diff = 0;
    Size n_obs = 0;         ///< 有效样本 T−1−k (CI6)
    std::string cv_source;  ///< "MHM96" (OL1992 独立查表 API 双对照)
};

namespace detail {

// statsmodels detrend: order=-1 原样; 0 去均值; 1 线性去势
// (vander(linspace(-1,1,T), order+1) 的列张成空间等价)
inline Eigen::MatrixXd johansen_detrend(const Eigen::MatrixXd& y, int order) {
    if (order == -1) return y;
    const Size t = static_cast<Size>(y.rows());
    const Size k = static_cast<Size>(y.cols());
    Eigen::MatrixXd x(t, order + 1);
    if (order == 0) {
        x.setOnes();
    } else {
        // linspace(-1, 1, T) 的 [t, 1] 列 (张成空间与 vander 一致)
        for (Size i = 0; i < t; ++i) {
            const Real u = -1.0 + 2.0 * static_cast<Real>(i) /
                                      static_cast<Real>(t - 1);
            x(i, 0) = u;
            x(i, 1) = 1.0;
        }
    }
    Eigen::MatrixXd beta = x.colPivHouseholderQr().solve(y);
    return y - x * beta;
}

// OLS 残差: y − X·lstsq(X, y) (statsmodels resid(y, x); pinv 与 QR 残差等价)
inline Eigen::MatrixXd ols_residual(const Eigen::MatrixXd& y,
                                    const Eigen::MatrixXd& x) {
    if (x.size() == 0) return y;
    Eigen::MatrixXd beta = x.colPivHouseholderQr().solve(y);
    return y - x * beta;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Johansen 协整秩检验 (API 等价 statsmodels coint_johansen)
//
// @param endog K 列水平数据 (K ≥ 2)
// @param det_order −1 (无确定项) / 0 (常数) / 1 (线性趋势) (CI4: 3 情形)
// @param k_ar_diff 滞后差分阶数 (≥ 0; urca 对照时 K_urca = k_ar_diff + 1)
// ---------------------------------------------------------------------------
inline JohansenResult coint_johansen(const var::MultivariateTSData& endog,
                                     int det_order, Size k_ar_diff) {
    if (det_order < -1 || det_order > 1) {
        throw std::invalid_argument(
            "coint_johansen: det_order must be -1/0/1");
    }
    endog.validate();
    const Size n = endog.K();
    if (n < 2) {
        throw std::invalid_argument("coint_johansen: K must be >= 2");
    }
    const Size t_tot = endog.T();
    if (t_tot < k_ar_diff + 3) {
        throw std::invalid_argument(
            "coint_johansen: sample too small for k_ar_diff");
    }

    Eigen::MatrixXd y = endog.matrix();  // T×K

    // f: 二轮去势阶 (SM: det_order > -1 → f=0 仅去均值; 否则 -1 不去势)
    const int f = (det_order > -1) ? 0 : -1;

    Eigen::MatrixXd e = detail::johansen_detrend(y, det_order);  // T×K

    // dx = diff(e): (T−1)×K
    const Size tdm = t_tot - 1;
    Eigen::MatrixXd dx(tdm, n);
    for (Size i = 0; i < tdm; ++i) dx.row(i) = e.row(i + 1) - e.row(i);

    // z = lagmat(dx, k)[k:]: 行 t (k..T−2) = [Δy_t, Δy_{t−1}, ..., Δy_{t−k+1}]
    // (对齐后回归元为 Δy_{t−1..t−k}; diff 报告附录 B)
    const Size t_eff = tdm - k_ar_diff;  // T−1−k
    Eigen::MatrixXd z(t_eff, n * k_ar_diff);
    if (k_ar_diff > 0) {
        for (Size i = 0; i < t_eff; ++i) {
            const Size t_reg = k_ar_diff + i;  // dx 行索引
            for (Size j = 0; j < k_ar_diff; ++j) {
                for (Size v = 0; v < n; ++v) {
                    z(i, j * n + v) = dx(t_reg - 1 - j, v);
                }
            }
        }
    }

    Eigen::MatrixXd dx_trim = dx.bottomRows(
        static_cast<Eigen::Index>(t_eff));
    dx_trim = detail::johansen_detrend(dx_trim, f);
    z = detail::johansen_detrend(z, f);

    Eigen::MatrixXd r0t = detail::ols_residual(dx_trim, z);  // T_eff×K

    // lx = e[1 : T−k] (水平, T_eff 行); ≡ y_{t−1} (span 恒等, diff 报告 §4)
    Eigen::MatrixXd lx(t_eff, n);
    for (Size i = 0; i < t_eff; ++i) lx.row(i) = e.row(i + 1);
    lx = detail::johansen_detrend(lx, f);
    Eigen::MatrixXd rkt = detail::ols_residual(lx, z);  // T_eff×K

    const Real tt = static_cast<Real>(t_eff);
    Eigen::MatrixXd skk = rkt.transpose() * rkt / tt;   // S₁₁
    Eigen::MatrixXd sk0 = rkt.transpose() * r0t / tt;   // S₁₀
    Eigen::MatrixXd s00 = r0t.transpose() * r0t / tt;   // S₀₀

    // A = S₁₁⁻¹·S₁₀S₀₀⁻¹S₀₁ (非对称; 特征值 = 平方典型相关)
    Eigen::MatrixXd sig = sk0 * s00.inverse() * sk0.transpose();
    Eigen::MatrixXd a = skk.inverse() * sig;

    Eigen::EigenSolver<Eigen::MatrixXd> es(a);
    const Eigen::Index nn = static_cast<Eigen::Index>(n);
    Eigen::VectorXd lam_real(nn);
    for (Eigen::Index i = 0; i < nn; ++i) lam_real(i) = es.eigenvalues()(i).real();
    Eigen::MatrixXd du = es.eigenvectors().real();  // 列 ↔ 特征值 (未排序)

    // 降序排序 (CI6)
    std::vector<Size> order(n);
    for (Size i = 0; i < n; ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](Size a_, Size b_) {
        return lam_real(static_cast<Eigen::Index>(a_)) >
               lam_real(static_cast<Eigen::Index>(b_));
    });
    Eigen::VectorXd a_sorted(nn);
    Eigen::MatrixXd du_sorted(nn, nn);
    for (Size i = 0; i < n; ++i) {
        a_sorted(static_cast<Eigen::Index>(i)) =
            lam_real(static_cast<Eigen::Index>(order[i]));
        du_sorted.col(static_cast<Eigen::Index>(i)) =
            du.col(static_cast<Eigen::Index>(order[i]));
    }

    // evec 归一: dt = du·inv(chol(du'S₁₁du)) → β'S₁₁β = I
    // (对已排序 du; 与 SM 先归一后排序等价 — 排序是列置换)
    Eigen::MatrixXd m = du_sorted.transpose() * skk * du_sorted;
    Eigen::LLT<Eigen::MatrixXd> llt(m);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("coint_johansen: S11 Cholesky failed");
    }
    Eigen::MatrixXd dt = du_sorted * llt.matrixL().solve(
        Eigen::MatrixXd::Identity(nn, nn));

    // 符号约定: 行主序首个非零元为正 (SM d.flat 扫描; 列符号仍任意 — CI8)
    Real first_nz = 0.0;
    for (Eigen::Index j = 0; j < nn && first_nz == 0.0; ++j) {
        for (Eigen::Index i = 0; i < nn; ++i) {
            if (dt(i, j) != 0.0) {
                first_nz = dt(i, j);
                break;
            }
        }
    }
    if (first_nz < 0.0) dt = -dt;

    // 统计量
    JohansenResult res;
    res.eig = a_sorted;
    res.evec = dt;
    res.lr1.resize(nn);
    res.lr2.resize(nn);
    for (Size r = 0; r < n; ++r) {
        Real acc = 0.0;
        for (Size i = r; i < n; ++i) {
            acc += std::log(1.0 - a_sorted(static_cast<Eigen::Index>(i)));
        }
        res.lr1(static_cast<Eigen::Index>(r)) = -tt * acc;
        res.lr2(static_cast<Eigen::Index>(r)) =
            -tt * std::log(1.0 - a_sorted(static_cast<Eigen::Index>(r)));
    }
    res.cvt.resize(nn, 3);
    res.cvm.resize(nn, 3);
    for (Size r = 0; r < n; ++r) {
        const auto cvt = mhm96_trace_cv(det_order, n - r);
        const auto cvm = mhm96_maxeig_cv(det_order, n - r);
        for (Size c = 0; c < 3; ++c) {
            res.cvt(static_cast<Eigen::Index>(r),
                    static_cast<Eigen::Index>(c)) = cvt[c];
            res.cvm(static_cast<Eigen::Index>(r),
                    static_cast<Eigen::Index>(c)) = cvm[c];
        }
    }
    res.det_order = det_order;
    res.k_ar_diff = k_ar_diff;
    res.n_obs = t_eff;
    res.cv_source = "MHM96";
    return res;
}

// ---------------------------------------------------------------------------
// 协整秩选择 (复刻 statsmodels select_coint_rank: 逐级检验至首次不拒绝)
//
// @param method "trace" / "maxeig"
// @param signif 0.1 / 0.05 / 0.01
// @return 建议秩 r (lr[r] < cv[r, signif] 即接受 r; 复用 cvt/cvm 同表, B4)
// ---------------------------------------------------------------------------
inline Size select_coint_rank(const var::MultivariateTSData& endog,
                              int det_order, Size k_ar_diff,
                              const std::string& method = "trace",
                              Real signif = 0.05) {
    if (method != "trace" && method != "maxeig") {
        throw std::invalid_argument(
            "select_coint_rank: method must be trace/maxeig");
    }
    Size col;
    if (signif == 0.1) {
        col = 0;
    } else if (signif == 0.05) {
        col = 1;
    } else if (signif == 0.01) {
        col = 2;
    } else {
        throw std::invalid_argument(
            "select_coint_rank: signif must be 0.1/0.05/0.01");
    }
    const JohansenResult jr = coint_johansen(endog, det_order, k_ar_diff);
    const Size n = endog.K();
    const Eigen::VectorXd& stat = (method == "trace") ? jr.lr1 : jr.lr2;
    const Eigen::MatrixXd& cv = (method == "trace") ? jr.cvt : jr.cvm;
    Size r = 0;
    while (r < n) {
        if (stat(static_cast<Eigen::Index>(r)) <
            cv(static_cast<Eigen::Index>(r),
               static_cast<Eigen::Index>(col))) {
            break;
        }
        ++r;
    }
    return r;
}

}  // namespace cointegration
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
