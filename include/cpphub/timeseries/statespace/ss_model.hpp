// =============================================================================
// ss_model.hpp - 状态空间模型 + 常规 Kalman 滤波 (DK Ch4 更新形) (spec §3.1)
//
// Phase 8 v1.8 M2 Step 2 (PHASE8_SPEC.md v1.1 §3.1)
//
// 教材锚点: Durbin-Koopman 2012 *Time Series Analysis by State Space Methods*
//   2nd ed Ch.4 (常规 KF, 更新形 P 递推) / Ch.5.6-5.7 (diffuse — 本期仅
//   approximate, SS1) / Harvey 1989 (ARIMA 状态空间, arima_statespace.hpp
//   消费本载体)
//
// probe_sarimax.py 冻结语义 (2026-08-22 全绿, tests/fixtures/timeseries/):
//   - 更新形递推 (probe manual_kf 与 statsmodels 0.14.6 数值等价实证 [D]:
//     ARMA(1,1) llf diff 2.81e-9 / ARIMA(1,1,1) 2.81e-8):
//       v_t = y_t − Z a_t − d;   F_t = Z P_t Z' + H
//       K_t = T P_t Z' F_t⁻¹;    a_{t+1} = T a_t + K_t v_t + c
//       P_{t+1} = T P_t (T − K_t Z)' + R Q R'
//     ⚠️ DK 更新形 — 与 Joseph 形/标准滤波形代数等价但增益编排不同, 不可
//     混搭 (SS2); K_t 含 T 前乘 (DK 预测增益, 非滤波增益)
//   - loglik 贡献 (t ≥ burn): −½(log F_t + v_t²/F_t + log 2π) (SS8 单变量
//     常数); burn = 贡献求和起点, 滤波全程运行 (innovations/pred_var
//     全长记录, 非丢观测 — probe [B 补充] 实证)
//   - 初始化 [E]: Stationary → a₁ = (I−T)⁻¹c, P₁ = Lyapunov 解;
//     ApproximateDiffuse → a₁ = 0, P₁ = κI (κ=1e6, statsmodels
//     approximate_diffuse_variance 默认口径); ARIMA 混合型 (κI_d ⊕
//     Lyapunov(ARMA 块)) 由 arima_statespace 经 Known 注入 — 本层不
//     硬编码块结构
//
// 幻觉点防护 (spec §7.2):
//   SS1: ApproximateDiffuse κ=1e6 为有限方差近似 (≠ exact diffuse κ→∞ 解析
//        极限, v1.9); FKV 2010: diffuse 似然不可普遍用于参数估计 — 本层
//        仅作初始化语义, 无似然目标函数角色
//   SS2: T 方向/更新形 (Harvey 形 AR 入第一列 + Z=(1,0,…,0) 由
//        arima_statespace 构造层负责)
//   SS3: R 为 m×r (1 ≤ r ≤ m) 非方阵, Q 为 r×r — R Q R' 秩 r 注入;
//        维度失配 → invalid_argument
//   SS4: Lyapunov 解后对称化 + LLT 正定守卫; 非正定 (单位根邻域 φ≈1.01
//        静默负方差) / 奇异 (单位根) → throw runtime_error
//   SS6: ARMA H=0 时 F 正定性全靠 P — F_t ≤ 0 (含 NaN) → throw
//   SS8: 单变量似然常数 = log 2π (多元才乘 p)
//
// 接口: header-only inline 自由函数 + 聚合体 (spec §1.4); Real/Size 出
//   core/types.hpp; Eigen3 载体 cpphub_timeseries_mat (K5, §6.2)
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace statespace {

/// 初始化方法 (spec §3.1; 无 ExactDiffuse 项 = SS1 边界, v1.9 评估)
enum class SSInitMethod { Known, Stationary, ApproximateDiffuse };

/// 初始化设定
struct SSInitSpec {
    SSInitMethod method = SSInitMethod::ApproximateDiffuse;
    Eigen::VectorXd a0;       ///< Known: E[α₁] (空 = 零向量)
    Eigen::MatrixXd P0;       ///< Known: Var[α₁] (必填, m×m)
    Real kappa = 1e6;         ///< ApproximateDiffuse 方差 (statsmodels 口径, K1)
};

/// 状态空间模型 (DK 记号, 单变量):
///   y_t = Z α_t + d + ε_t,          ε_t ~ N(0, H)      (H: 1×1)
///   α_{t+1} = T α_t + c + R η_t,    η_t ~ N(0, Q)      (R: m×r, Q: r×r, SS3)
struct StateSpaceModel {
    Eigen::MatrixXd Z;        ///< (1×m) 观测矩阵
    Eigen::MatrixXd T;        ///< (m×m) 状态转移
    Eigen::MatrixXd R;        ///< (m×r) 状态噪声载入 (1 ≤ r ≤ m)
    Eigen::MatrixXd Q;        ///< (r×r) 状态噪声方差
    Eigen::MatrixXd H;        ///< (1×1) 观测噪声方差 (ARMA: 零阵, SS6)
    Eigen::VectorXd d;        ///< (1) 观测截距 (空 = 0)
    Eigen::VectorXd c;        ///< (m) 状态截距 (空 = 0)
    SSInitSpec init{};
};

/// Kalman 滤波结果
struct KalmanResult {
    Real loglik = 0.0;                 ///< burn 后 Σ −½(log F_t + v_t²/F_t + log 2π)
    std::vector<Real> loglik_contribs; ///< 逐 t 贡献 (仅 t ≥ burn)
    std::vector<Real> innovations;     ///< v_t (burn 前后均记录)
    std::vector<Real> pred_var;        ///< F_t (全长)
    std::vector<Eigen::VectorXd> state_pred;      ///< a_{t+1|t} (t=0..T−1;
                                       ///< state_pred[0] = 初始先验 a₁;
                                       ///< innovations[t] = y[t] − Z·state_pred[t] − d)
    std::vector<Eigen::MatrixXd> state_cov_pred;  ///< P_{t+1|t}, 同上
    Size n_burn = 0;
    Size n_obs = 0;
    bool ok = false;
    std::string message;                ///< 非空 = Stationary→Diffuse 回退等记录
};

namespace detail {

/// 系统矩阵维度/有限性校验 (SS3; 失配 → invalid_argument)
inline void validate_model(const StateSpaceModel& m) {
    const Size ms = static_cast<Size>(m.T.rows());
    if (ms == 0 || m.T.cols() != m.T.rows())
        throw std::invalid_argument(
            "statespace: T must be a non-empty square matrix");
    const Size r = static_cast<Size>(m.R.cols());
    if (m.R.rows() != ms || r == 0 || r > ms)
        throw std::invalid_argument(
            "statespace: R must be (m x r) with 1 <= r <= m");
    if (m.Q.rows() != r || m.Q.cols() != r)
        throw std::invalid_argument("statespace: Q must be (r x r)");
    if (m.Z.rows() != 1 || m.Z.cols() != ms)
        throw std::invalid_argument("statespace: Z must be (1 x m)");
    if (m.H.rows() != 1 || m.H.cols() != 1)
        throw std::invalid_argument("statespace: H must be (1 x 1)");
    if (m.d.size() != 0 && m.d.size() != 1)
        throw std::invalid_argument("statespace: d must be empty or size 1");
    if (m.c.size() != 0 && m.c.size() != ms)
        throw std::invalid_argument("statespace: c must be empty or size m");
    if (!(m.Z.allFinite() && m.T.allFinite() && m.R.allFinite() &&
          m.Q.allFinite() && m.H.allFinite()))
        throw std::invalid_argument(
            "statespace: system matrices must be finite");
}

/// T 的谱半径 max|λ| (复特征值取模)
inline Real spectral_radius(const Eigen::MatrixXd& T) {
    Eigen::EigenSolver<Eigen::MatrixXd> es(T);
    Real rho = 0.0;
    for (Size i = 0; i < static_cast<Size>(es.eigenvalues().size()); ++i)
        rho = std::max(rho, std::abs(es.eigenvalues()(i)));
    return rho;
}

}  // namespace detail

/// 平稳初始化协方差: 解离散 Lyapunov 方程 P = T P T' + R Q R'
/// (直接法: vec(P) = (I − T⊗T)⁻¹ vec(RQR'), 列主序 vec)
/// 守卫 (SS4): fullPivLu 可逆性 + 解后对称化 + LLT 正定校验;
///   非正定 (单位根邻域 φ≈1.01 静默负方差) / 奇异 (单位根) → throw
///   runtime_error
inline Eigen::MatrixXd stationary_cov(const StateSpaceModel& m) {
    detail::validate_model(m);
    const Size ms = static_cast<Size>(m.T.rows());
    const Eigen::MatrixXd rqr = m.R * m.Q * m.R.transpose();

    // vec(AXB) = (B'⊗A)vec(X) → P = TPT' + RQR' 化为 (I − T⊗T)vec(P) = vec(RQR')
    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(ms * ms, ms * ms);
    for (Size i = 0; i < ms; ++i)
        for (Size j = 0; j < ms; ++j)
            for (Size k = 0; k < ms; ++k)
                for (Size l = 0; l < ms; ++l)
                    A(i * ms + k, j * ms + l) -= m.T(i, j) * m.T(k, l);
    Eigen::VectorXd rhs(ms * ms);
    for (Size i = 0; i < ms; ++i)
        for (Size j = 0; j < ms; ++j)
            rhs(i + j * ms) = rqr(i, j);

    Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
    if (!lu.isInvertible())
        throw std::runtime_error(
            "stationary_cov: singular Lyapunov system (unit-root T?)");
    const Eigen::VectorXd vec_p = lu.solve(rhs);

    Eigen::MatrixXd P(ms, ms);
    for (Size i = 0; i < ms; ++i)
        for (Size j = 0; j < ms; ++j)
            P(i, j) = vec_p(i + j * ms);
    P = 0.5 * (P + P.transpose());  // 对称化
    if (!P.allFinite())
        throw std::runtime_error(
            "stationary_cov: non-finite Lyapunov solution");
    Eigen::LLT<Eigen::MatrixXd> llt(P);
    if (llt.info() != Eigen::Success)
        throw std::runtime_error(
            "stationary_cov: Lyapunov solution not positive definite "
            "(spectral radius of T >= 1 neighbourhood, SS4)");
    return P;
}

/// 常规 Kalman 滤波 (DK Ch4 更新形; probe_sarimax.py manual_kf 镜像)
///
/// init.method=Stationary 且 T 谱半径 ≥ 1 → 内部转 ApproximateDiffuse 并
/// message 记录 (非静默, spec §3.1); burn = loglik 贡献求和起点 (滤波
/// 全程运行, 非丢观测)
inline KalmanResult kalman_filter(const std::vector<Real>& y,
                                  const StateSpaceModel& model,
                                  Size loglikelihood_burn = 0) {
    detail::validate_model(model);
    const Size n = static_cast<Size>(y.size());
    if (n == 0)
        throw std::invalid_argument("statespace: empty observation vector");
    for (Real v : y)
        if (!std::isfinite(v))
            throw std::invalid_argument(
                "statespace: y contains non-finite values");

    KalmanResult res;
    res.n_burn = loglikelihood_burn;
    res.n_obs = n;

    const Size ms = static_cast<Size>(model.T.rows());
    const Eigen::VectorXd zt = model.Z.transpose();  // (m×1)
    const Eigen::MatrixXd rqr = model.R * model.Q * model.R.transpose();
    const Real h = model.H(0, 0);
    const Real d0 = model.d.size() == 0 ? 0.0 : model.d(0);
    const Eigen::VectorXd cvec =
        model.c.size() == 0 ? Eigen::VectorXd::Zero(ms) : model.c;
    const Real log_2pi = std::log(2.0 * 3.14159265358979323846);

    // --- Step 1: 初始化 (三法) ---
    Eigen::VectorXd a(ms);
    Eigen::MatrixXd P(ms, ms);
    switch (model.init.method) {
        case SSInitMethod::Known: {
            if (model.init.P0.rows() != ms || model.init.P0.cols() != ms)
                throw std::invalid_argument(
                    "statespace: Known init requires P0 of size (m x m)");
            if (model.init.a0.size() != 0 &&
                model.init.a0.size() != ms)
                throw std::invalid_argument(
                    "statespace: a0 must be empty or size m");
            a = model.init.a0.size() == 0 ? Eigen::VectorXd::Zero(ms)
                                          : model.init.a0;
            P = model.init.P0;
            break;
        }
        case SSInitMethod::Stationary: {
            // 谱半径守卫: |λ| ≥ 1 (容差 1e-12 捕数值单位根) → 回退 diffuse
            if (detail::spectral_radius(model.T) >= 1.0 - 1e-12) {
                a.setZero(ms);
                P = model.init.kappa * Eigen::MatrixXd::Identity(ms, ms);
                res.message =
                    "statespace: stationary init requested but spectral "
                    "radius of T >= 1; fell back to ApproximateDiffuse "
                    "(kappa="
                    + std::to_string(model.init.kappa) + ")";
                break;
            }
            // a₁ = (I−T)⁻¹c (c 空则 0); P₁ = Lyapunov 解
            const Eigen::MatrixXd imt =
                Eigen::MatrixXd::Identity(ms, ms) - model.T;
            a = imt.partialPivLu().solve(cvec);
            P = stationary_cov(model);
            break;
        }
        case SSInitMethod::ApproximateDiffuse: {
            a.setZero(ms);
            P = model.init.kappa * Eigen::MatrixXd::Identity(ms, ms);
            break;
        }
    }

    // --- Step 2-4: 逐 t 递推 (DK Ch4 更新形) ---
    for (Size t = 0; t < n; ++t) {
        res.state_pred.push_back(a);
        res.state_cov_pred.push_back(P);

        const Real v = y[t] - (model.Z * a)(0, 0) - d0;
        const Real F = (model.Z * P * zt)(0, 0) + h;
        if (!(F > 0.0))
            throw std::runtime_error(
                "statespace: non-positive prediction variance F_t at t=" +
                std::to_string(t) + " (SS6)");
        res.innovations.push_back(v);
        res.pred_var.push_back(F);
        if (t >= loglikelihood_burn) {
            const Real contrib = -0.5 * (std::log(F) + v * v / F + log_2pi);
            res.loglik_contribs.push_back(contrib);
            res.loglik += contrib;
        }

        const Eigen::VectorXd K = model.T * P * zt / F;  // DK 预测增益
        const Eigen::VectorXd a_next = model.T * a + K * v + cvec;
        const Eigen::MatrixXd tmks = model.T - K * model.Z;  // (T − K Z)
        P = model.T * P * tmks.transpose() + rqr;  // ⚠️ 更新形 (非 Joseph)
        a = a_next;
    }
    res.ok = true;
    return res;
}

}  // namespace statespace
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
