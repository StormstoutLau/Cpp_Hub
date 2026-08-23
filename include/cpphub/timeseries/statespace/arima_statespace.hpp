// =============================================================================
// arima_statespace.hpp - ARMA Harvey 形 + ARIMA 嵌入差分 + 精确似然 MLE (spec §3.2)
//
// Phase 8 v1.8 M2 Step 3 (PHASE8_SPEC.md v1.1 §3.2)
//
// 教材锚点: Harvey 1989 *Forecasting, Structural Time Series Models and
//   the Kalman Filter* (ARMA 状态空间, Harvey 形) / Durbin-Koopman 2012
//   2nd ed Ch.4 (KF, ss_model.hpp 载体) / Monahan 1984 Biometrika
//   71(2):403-404 (AR 平稳化变换) / Box-Jenkins 5th ed Ch.7 (ARIMA)
//
// probe 冻结语义 (probe_sarimax.py + probe_d2_const.py + verify_drift_init.py
//   2026-08-22/23 全绿, tests/fixtures/timeseries/):
//   [布局] 状态 = [差分块 (d) | ARMA companion 块 (m₁ = max(p, q+1))],
//     m = d + m₁; Harvey 形 = AR 系数入 T 第一列 (上 Hessenberg companion)
//     + Z = [1]×d ⊕ [1, 0, …, 0] (hamilton_representation=False 默认, SS2)
//   [嵌入差分 T] 差分块上三角全 1 (triu) + 链接列 T[:d, d] = 1 + ARMA 块
//     companion; 黄金锚 ARIMA(1,1,1) (φ, θ):
//     T = [[1,1,0],[0,φ,1],[0,0,0]], R = [0,1,θ]', Z = [1,1,0]
//     (probe [A2] dT/dR = 0 逐元素实证; d=2 三组合同款, probe_d2_const [P2])
//   [常数注入] c[d] = 常数 (d=0: ARMA 方程截距 = state_intercept[0],
//     probe_d2_const [P1] llf diff 5.68e-14; d≥1: drift = 差分方程常数,
//     verify_drift_init 混合初始化 a₁^ARMA = (I−T_sub)⁻¹c_sub 实证) —
//     均走状态截距 c, 非观测截距 d_t
//   [初始化] d=0: Stationary (a₁=(I−T)⁻¹c, P₁=Lyapunov); d≥1: 混合 =
//     差分块 diffuse κ=1e6 ⊕ ARMA 块 stationary (Known 注入, probe [E]);
//     AR 非平稳 → 全 diffuse ApproximateDiffuse (statsmodels
//     enforce_stationarity=False 口径, probe [B])
//   [burn] 默认 = d (SPEC 开放问题 #3 裁决 = statsmodels 默认公式
//     k_states − _k_order = d, probe [B] 八组全验); AR 非平稳全 diffuse →
//     burn = d + m₁ (statsmodels k_states 口径); simple = 0
//   [Monahan] constrain(u): r[k]=u[k]/√(1+u[k]²), y[k][i]=y[k−1][i]+
//     r[k]·y[k−1][k−i−1], 输出 −y[n−1]·; MA = **−**constrain(u) (probe
//     [C] 负号实锤, 小值区 θ≈+u); σ² = u_σ² (transform_params 同构)
//   [SS5 双口径] embedded burn=d vs simple (Δ^d y 样本 ARMA, burn=0)
//     似然时间平移等价, gap = O(1/κ) ≈ 5e-7 (κ=1e6 有限 diffuse 尾差,
//     ss_baseline GOLD_LL_ARIMA111 vs _SIMPLE); burn=0 含 diffuse 伪贡献
//     (−½log κ ≈ −7.8/首期) 不可混比
//   [K4] σ² 非集中化全参数优化 (σ² = u_σ²); statsmodels fit 默认
//     concentrate_scale=True — 对拍必须同口径 (ss_baseline 全 False)
//
// 幻觉点防护 (spec §7.2):
//   SS1/SS9: approximate diffuse κ=1e6 仅作初始化语义 (FKV 2010: diffuse
//     似然不可普遍用于参数估计 — 非似然目标函数角色)
//   SS2: Harvey 形 (AR 入 T 第一列 + Z 前置); Hamilton 形 (AR 入第一行 +
//     Z 后置) 混用 = 模型级错误 — 本文件唯一 Harvey 形
//   SS5: 两差分口径似然不同 — 双跑差异显式回显, 禁混比断言
//   K1: innovations 常驻互检 (fast 带状 pieces; T=300 回归锚 2.6e-8,
//     告警线 1e-5; 含常数模型口径差异声明 — KF 参数化常数 vs 样本均值)
//
// v1.8 边界 (spec §8.5 不做清单): p=q=0 纯白噪声/纯随机游走不支持
//   (statsmodels _k_order=0 退化路径; throw invalid_argument);
//   exog / SARIMA 季节 / 缺失数据 / exact diffuse / 平滑器不做;
//   实施期签名增补: make_arima_statespace 追加第 7 defaulted 参数
//   constant_value (bool include_constant 无法携带数值, 混合初始化
//   a₀ = (I−T_sub)⁻¹c_sub 需要常数数值 — verify_drift_init 冻结语义;
//   纯 defaulted 追加非破坏, devlog 记录)
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/calibration/optimizer.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/arima/arima_model.hpp"
#include "cpphub/timeseries/arima/innovations_mle.hpp"
#include "cpphub/timeseries/statespace/ss_model.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace statespace {

/// ARIMA statespace MLE 配置 (spec §3.2 冻结签名)
struct ArimaSSConfig {
    bool simple_differencing = false;    ///< v1.8 默认 false = 嵌入差分 (K1 四件套之一); true 口径对照用例
    bool include_constant = false;       ///< d=0: ARMA 方程截距; d≥1: 差分方程常数 = levels drift (经 c[d] 注入)
    bool enforce_stationarity = true;    ///< AR 平稳化变换 (Monahan 1984, statsmodels _ar_trans_params 同构)
    bool enforce_invertibility = true;   ///< MA 可逆化变换 (_ma_trans_params 同构, 负号 probe [C])
    Size n_starting_points = 2;          ///< {innovations 解, CSS 解} 确定性初值 (7C 多起始惯例)
    Size seed = 42;                      ///< 预留 (本期无随机扰动)
    SLSQP::Config optimizer_config = [] {  ///< K3: 复用 ADR-018 (arima CSS-ML
        SLSQP::Config c;                   ///< 同款收紧 — SLSQP 裸默认 ftol=1e-6
        c.max_iterations = 300;            ///< 落点差 4e-5, 7C 惯例 1e-11 实测
        c.ftol = 1e-11;                    ///< llf 1e-6 级; MLE 对拍阈值依据)
        c.xtol = 1e-10;
        return c;
    }();
    bool cross_check_innovations = true; ///< K1: 常驻互检
};

/// ARIMA statespace MLE 结果 (spec §3.2 冻结签名)
struct ArimaSSResult {
    std::vector<Real> phi;         ///< AR 系数 (正号, 1 − φ₁B − …)
    std::vector<Real> theta;       ///< MA 系数 (正号, 1 + θ₁B + …, AR1)
    Real constant = 0.0;           ///< 实施期增补: 常数回显 (d=0 截距 / d≥1
                                   ///< drift; inc_const=false 时 0) — MLE
                                   ///< drift 用例参数对比需要 (make_arima_
                                   ///< statespace constant_value 同款增补)
    Real sigma2 = 0.0;             ///< 非集中化全参数优化 (K4)
    Real loglik = 0.0;             ///< KF loglik (burn 口径回显)
    Real aic = 0.0;                ///< −2·loglik + 2·(p+q+1+const) (statsmodels 对齐)
    Size loglikelihood_burn = 0;
    bool converged = false;
    Size n_iterations = 0;
    std::string message;
    Real innovations_loglik = 0.0; ///< 互检: innovations (T−d 差分样本) loglik
    Real innovations_gap = 0.0;    ///< |KF − innovations|; > 1e-5 → message 告警 (常驻协议)
    std::vector<Real> innovations; ///< v_t (KF 全长, 诊断)
};

namespace detail {

constexpr Real kPiSs = 3.14159265358979323846;

/// 值域裁剪 (Monahan 域防护 |u| ≤ 10 ⇔ |φ| ≤ 0.995)
inline Real clip_u(Real v) {
    if (!std::isfinite(v)) return 0.0;  // 非平稳解 → 零 (statsmodels L198-204 置零语义)
    return std::max(-10.0, std::min(10.0, v));
}

/// Monahan (1984) AR 平稳化正变换 — statsmodels constrain_stationary_
/// univariate 逐字: r[k] = u[k]/√(1+u[k]²); y[k][i] = y[k−1][i] +
/// r[k]·y[k−1][k−i−1] (i < k); y[k][k] = r[k]; 返回 −y[n−1]·
inline std::vector<Real> constrain_stationary(const std::vector<Real>& u) {
    const Size n = u.size();
    std::vector<std::vector<Real>> y(n, std::vector<Real>(n, 0.0));
    std::vector<Real> r(n, 0.0);
    for (Size k = 0; k < n; ++k)
        r[k] = u[k] / std::sqrt(1.0 + u[k] * u[k]);
    for (Size k = 0; k < n; ++k) {
        for (Size i = 0; i < k; ++i)
            y[k][i] = y[k - 1][i] + r[k] * y[k - 1][k - i - 1];
        y[k][k] = r[k];
    }
    std::vector<Real> out(n, 0.0);
    for (Size i = 0; i < n; ++i) out[i] = -y[n - 1][i];
    return out;
}

/// Monahan 逆变换 — unconstrain_stationary_univariate 逐字 (整组联合,
/// 不可逐系数); |c| → 1 时 u → ∞ (调用方 clip)
inline std::vector<Real> unconstrain_stationary(const std::vector<Real>& c) {
    const Size n = c.size();
    std::vector<std::vector<Real>> y(n, std::vector<Real>(n, 0.0));
    for (Size i = 0; i < n; ++i) y[n - 1][i] = -c[i];
    for (Size k = n - 1; k >= 1; --k) {
        for (Size i = 0; i < k; ++i)
            y[k - 1][i] =
                (y[k][i] - y[k][k] * y[k][k - i - 1]) / (1.0 - y[k][k] * y[k][k]);
    }
    std::vector<Real> x(n, 0.0);
    for (Size k = 0; k < n; ++k) {
        const Real rk = y[k][k];
        x[k] = rk / std::sqrt(1.0 - rk * rk);
    }
    return x;
}

/// d 阶差分
inline std::vector<Real> difference(const std::vector<Real>& y, Size d) {
    std::vector<Real> z = y;
    for (Size it = 0; it < d; ++it) {
        if (z.size() < 2)
            throw std::invalid_argument("difference: series too short");
        std::vector<Real> nx(z.size() - 1, 0.0);
        for (Size t = 1; t < z.size(); ++t) nx[t - 1] = z[t] - z[t - 1];
        z.swap(nx);
    }
    return z;
}

/// 样本均值
inline Real mean_of(const std::vector<Real>& v) {
    Real s = 0.0;
    for (Real x : v) s += x;
    return v.empty() ? 0.0 : s / static_cast<Real>(v.size());
}

/// AR companion 谱半径 (平稳性判定: ρ < 1; p=0 恒 0)
inline Real ar_spectral_radius(const std::vector<Real>& phi) {
    const Size p = phi.size();
    if (p == 0) return 0.0;
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero(p, p);
    for (Size i = 0; i < p; ++i) c(i, 0) = phi[i];
    for (Size i = 0; i + 1 < p; ++i) c(i + 1, i) = 1.0;
    return spectral_radius(c);
}

/// 固定参数 innovations 精确 loglik (K1 互检对手; fast 带状 pieces 与
/// arima::innovations_mle 同构; 零均值口径 — demean 由调用方决定):
///   ℓ = −0.5·[n·log(2πσ²) + Σ log v_t + (Σ u_t²/v_t)/σ²]
inline Real innovations_loglik_at(const std::vector<Real>& z,
                                  const std::vector<Real>& phi,
                                  const std::vector<Real>& theta,
                                  Real sigma2) {
    const arima::detail::LoglikePieces pc =
        arima::detail::arma_innovations_pieces(z, phi, theta);
    const Real n = static_cast<Real>(z.size());
    return -0.5 * (n * std::log(2.0 * kPiSs * sigma2) + pc.sum_logv +
                   pc.s_uv / sigma2);
}

/// (φ, θ) → 无约束域 u 前 p+q 项 (enforce 双开关; θ = −constrain 负号)
inline void params_to_u(const std::vector<Real>& phi,
                        const std::vector<Real>& theta,
                        bool enforce_ar, bool enforce_ma,
                        std::vector<Real>& u_out) {
    u_out.clear();
    u_out.reserve(phi.size() + theta.size());
    if (!phi.empty()) {
        if (enforce_ar) {
            const auto u = unconstrain_stationary(phi);
            for (Real v : u) u_out.push_back(clip_u(v));
        } else {
            for (Real v : phi) u_out.push_back(clip_u(v));
        }
    }
    if (!theta.empty()) {
        if (enforce_ma) {
            std::vector<Real> neg(theta.size(), 0.0);
            for (Size i = 0; i < theta.size(); ++i) neg[i] = -theta[i];
            const auto u = unconstrain_stationary(neg);
            for (Real v : u) u_out.push_back(clip_u(v));
        } else {
            for (Real v : theta) u_out.push_back(clip_u(v));
        }
    }
}

}  // namespace detail

/// ARMA(p,q) Harvey 形系统矩阵 (d=0): m₁ = max(p, q+1);
/// T 第一列 AR + 次对角单位 (companion), R = (1, θ₁, …)', Z = (1, 0, …, 0),
/// H = 0 (SS6), Q = [σ²]; init = Stationary (KF 内 a₁=(I−T)⁻¹c, P₁=Lyapunov)
inline StateSpaceModel make_arma_statespace(const std::vector<Real>& phi,
                                            const std::vector<Real>& theta,
                                            Real sigma2) {
    const Size p = phi.size();
    const Size q = theta.size();
    if (p + q == 0)
        throw std::invalid_argument(
            "make_arma_statespace: p=q=0 white noise not supported (v1.8)");
    if (!(sigma2 > 0.0) || !std::isfinite(sigma2))
        throw std::invalid_argument("make_arma_statespace: sigma2 must be > 0");
    for (Real v : phi)
        if (!std::isfinite(v))
            throw std::invalid_argument("make_arma_statespace: NaN in phi");
    for (Real v : theta)
        if (!std::isfinite(v))
            throw std::invalid_argument("make_arma_statespace: NaN in theta");

    const Size m1 = std::max(p, q + 1);
    StateSpaceModel m;
    m.Z = Eigen::MatrixXd::Zero(1, m1);
    m.Z(0, 0) = 1.0;
    m.T = Eigen::MatrixXd::Zero(m1, m1);
    for (Size i = 0; i < p; ++i) m.T(i, 0) = phi[i];
    for (Size i = 0; i + 1 < m1; ++i) m.T(i, i + 1) = 1.0;  // 超对角 (Harvey)
    m.R = Eigen::MatrixXd::Zero(m1, 1);
    m.R(0, 0) = 1.0;
    for (Size j = 0; j < q && j + 1 < m1; ++j) m.R(j + 1, 0) = theta[j];
    m.Q = Eigen::MatrixXd::Constant(1, 1, sigma2);
    m.H = Eigen::MatrixXd::Zero(1, 1);
    m.init.method = SSInitMethod::Stationary;
    return m;
}

/// ARIMA(p,d,q) 嵌入差分系统矩阵 (simple_differencing=true → 先差分后
/// ARMA 形, 数据由调用方差分; include_constant=true → c 布局开启,
/// c[d] = constant_value — d=0: ARMA 方程截距, d≥1: drift, probe 冻结)
///
/// 嵌入差分 (d ≥ 1): m = d + m₁; T = [triu(1) 差分块 | 链接列 | ARMA
/// companion]; R = (0…0, 1, θ…)'; Z = [1]×d ⊕ [1, 0…]; 初始化 = 混合
/// (差分块 diffuse κ ⊕ ARMA 块 stationary, Known 注入) 或全 diffuse
/// (AR 非平稳, statsmodels enforce=False 口径)
inline StateSpaceModel make_arima_statespace(
        const std::vector<Real>& phi, const std::vector<Real>& theta, Size d,
        Real sigma2, bool include_constant, bool simple_differencing = false,
        Real constant_value = 0.0) {
    const Size p = phi.size();
    const Size q = theta.size();
    if (p + q == 0)
        throw std::invalid_argument(
            "make_arima_statespace: p=q=0 (pure random walk) not supported "
            "(v1.8)");
    if (!(sigma2 > 0.0) || !std::isfinite(sigma2))
        throw std::invalid_argument("make_arima_statespace: sigma2 must be > 0");
    for (Real v : phi)
        if (!std::isfinite(v))
            throw std::invalid_argument("make_arima_statespace: NaN in phi");
    for (Real v : theta)
        if (!std::isfinite(v))
            throw std::invalid_argument("make_arima_statespace: NaN in theta");
    if (!std::isfinite(constant_value))
        throw std::invalid_argument("make_arima_statespace: NaN constant");

    // simple / d=0: ARMA 形 (无差分块); 常数在 c[0] (ARMA 方程截距)
    if (simple_differencing || d == 0) {
        StateSpaceModel m = make_arma_statespace(phi, theta, sigma2);
        if (include_constant) {
            const Size m1 = static_cast<Size>(m.T.rows());
            m.c = Eigen::VectorXd::Zero(m1);
            m.c(0) = constant_value;
        }
        return m;
    }

    // ---- d ≥ 1 嵌入差分 ----
    const Size m1 = std::max(p, q + 1);
    const Size m = d + m1;

    StateSpaceModel s;
    s.Z = Eigen::MatrixXd::Zero(1, m);
    for (Size i = 0; i <= d; ++i) s.Z(0, i) = 1.0;

    s.T = Eigen::MatrixXd::Zero(m, m);
    for (Size i = 0; i < d; ++i)              // 差分块上三角全 1
        for (Size j = i; j < d; ++j) s.T(i, j) = 1.0;
    for (Size i = 0; i < d; ++i) s.T(i, d) = 1.0;  // 链接列 T[:d, d]
    for (Size i = 0; i < p; ++i) s.T(d + i, d) = phi[i];  // AR 首列
    for (Size i = 0; i + 1 < m1; ++i) s.T(d + i, d + i + 1) = 1.0;  // companion 超对角

    s.R = Eigen::MatrixXd::Zero(m, 1);
    s.R(d, 0) = 1.0;
    for (Size j = 0; j < q && j + 1 < m1; ++j) s.R(d + j + 1, 0) = theta[j];

    s.Q = Eigen::MatrixXd::Constant(1, 1, sigma2);
    s.H = Eigen::MatrixXd::Zero(1, 1);

    if (include_constant) {
        s.c = Eigen::VectorXd::Zero(m);
        s.c(d) = constant_value;  // drift 注入 ARMA 块首行 (index d)
    }

    // ---- 初始化 (probe [E] + verify_drift_init 冻结): AR 块平稳 →
    // 混合 = 差分块 diffuse κ ⊕ ARMA 块 stationary (Known); 非平稳 →
    // 全 diffuse ApproximateDiffuse ----
    if (detail::ar_spectral_radius(phi) < 1.0) {
        const Eigen::MatrixXd t_sub = s.T.block(d, d, m1, m1);
        const Eigen::MatrixXd r_sub = s.R.block(d, 0, m1, 1);

        StateSpaceModel sub;  // ARMA 子模型 (Lyapunov 载体)
        sub.Z = Eigen::MatrixXd::Zero(1, m1);
        sub.Z(0, 0) = 1.0;
        sub.T = t_sub;
        sub.R = r_sub;
        sub.Q = s.Q;
        sub.H = s.H;
        const Eigen::MatrixXd p_sub = stationary_cov(sub);

        Eigen::VectorXd c_sub = Eigen::VectorXd::Zero(m1);
        if (include_constant) c_sub(0) = constant_value;
        const Eigen::VectorXd a_sub =
            (Eigen::MatrixXd::Identity(m1, m1) - t_sub).partialPivLu().solve(c_sub);

        s.init.method = SSInitMethod::Known;
        s.init.a0 = Eigen::VectorXd::Zero(m);
        s.init.a0.segment(d, m1) = a_sub;
        s.init.P0 = Eigen::MatrixXd::Zero(m, m);
        for (Size i = 0; i < d; ++i) s.init.P0(i, i) = s.init.kappa;
        s.init.P0.block(d, d, m1, m1) = p_sub;
    } else {
        s.init.method = SSInitMethod::ApproximateDiffuse;
    }
    return s;
}

namespace detail {

/// MLE 优化变量解包 (u = [u_φ(p), u_θ(q), u_σ, (u_c)])
struct ArimaParamsUnpacked {
    std::vector<Real> phi;
    std::vector<Real> theta;
    Real sigma2 = 0.0;
    Real constant = 0.0;
};

inline ArimaParamsUnpacked unpack_params(const std::vector<Real>& u, Size p,
                                         Size q, bool enforce_ar,
                                         bool enforce_ma, bool has_const) {
    ArimaParamsUnpacked up;
    if (p > 0) {
        if (enforce_ar) {
            up.phi = constrain_stationary(
                std::vector<Real>(u.begin(), u.begin() + p));
        } else {
            up.phi.assign(u.begin(), u.begin() + p);
        }
    }
    if (q > 0) {
        if (enforce_ma) {
            const auto c = constrain_stationary(
                std::vector<Real>(u.begin() + p, u.begin() + p + q));
            up.theta.resize(q, 0.0);
            for (Size i = 0; i < q; ++i) up.theta[i] = -c[i];  // 负号 probe [C]
        } else {
            up.theta.assign(u.begin() + p, u.begin() + p + q);
        }
    }
    up.sigma2 = u[p + q] * u[p + q];
    up.constant = has_const ? u[p + q + 1] : 0.0;
    return up;
}

}  // namespace detail

/// ARIMA 精确似然 MLE (Kalman 路径, spec §3.2)
///
/// 流程: 确定性多起始 {innovations_mle 解, arima CSS 解} → SLSQP
/// (K3, ADR-018) 无约束域优化 (Monahan 变换 + σ² = u²) → 最优解互检
/// (K1: innovations 同参数 loglik, gap > 1e-5 → message 告警)。
///
/// 初始化: AR 平稳 (|特征根|<1) → Stationary/混合; 否则 → ApproximateDiffuse
/// κ=1e6 (K1); burn 默认 = d (开放问题 #3 裁决 = statsmodels 默认公式),
/// AR 非平稳全 diffuse → burn = d + m₁, simple_differencing → burn = 0
///
/// @param data levels 序列 (差分由模型嵌入或 simple 模式内部处理)
/// @param spec p/d/q + include_drift (与 config.include_constant 取或)
/// @throws std::invalid_argument 空数据/NaN/p=q=0/T 过小
inline ArimaSSResult arima_statespace_mle(
        const std::vector<Real>& data, const arima::ArimaSpec& spec,
        const ArimaSSConfig& config = ArimaSSConfig{}) {
    const Size p = spec.p;
    const Size d = spec.d;
    const Size q = spec.q;
    const bool inc_const = config.include_constant || spec.include_drift;
    const Size m1 = std::max(p, q + 1);

    // ---- 校验 (§1.4-5) ----
    if (data.empty())
        throw std::invalid_argument("arima_statespace_mle: empty data");
    for (Real v : data)
        if (!std::isfinite(v))
            throw std::invalid_argument("arima_statespace_mle: NaN in data");
    if (p + q == 0)
        throw std::invalid_argument(
            "arima_statespace_mle: p=q=0 not supported (v1.8)");
    if (data.size() < d + m1 + 5)
        throw std::invalid_argument("arima_statespace_mle: T too small");

    const std::vector<Real> z = detail::difference(data, d);
    const std::vector<Real>& y_use =
        config.simple_differencing ? z : data;
    const Real mean_z = detail::mean_of(z);

    // ---- 优化变量布局: u = [u_φ(p), u_θ(q), u_σ, (u_c)] ----
    const Size n_u = p + q + 1 + (inc_const ? 1 : 0);
    const Size idx_sigma = p + q;
    const Size idx_const = p + q + 1;

    // ---- 似然评估 (neg_ll): 变换 → 系统矩阵 → KF (burn 口径) ----
    auto burn_of = [&](const std::vector<Real>& ph) -> Size {
        const bool nonstat = p > 0 && detail::ar_spectral_radius(ph) >= 1.0;
        if (config.simple_differencing) return nonstat ? m1 : 0;
        return nonstat ? d + m1 : d;
    };
    auto neg_ll = [&](const std::vector<Real>& u) -> Real {
        try {
            const auto up =
                detail::unpack_params(u, p, q, config.enforce_stationarity,
                                      config.enforce_invertibility, inc_const);
            if (!(up.sigma2 > 1e-12) || !std::isfinite(up.sigma2))
                return 1e100;
            const StateSpaceModel m = make_arima_statespace(
                up.phi, up.theta, d, up.sigma2, inc_const,
                config.simple_differencing, up.constant);
            const KalmanResult kr =
                kalman_filter(y_use, m, burn_of(up.phi));
            if (!kr.ok || !std::isfinite(kr.loglik)) return 1e100;
            return -kr.loglik;
        } catch (const std::exception&) {
            return 1e100;  // 守卫触发 (Lyapunov/F 正定) = 不可行参数
        }
    };

    // ---- 确定性多起始 (7C 惯例; 无随机扰动, seed 预留) ----
    std::vector<std::vector<Real>> starts;
    {
        // start 1: innovations_mle 解 (demean=true = statsmodels 默认口径;
        // 常数起始由均值推算 — innovations 无常数参数)
        try {
            const arima::InnovationsResult ir = arima::innovations_mle(
                data, p, q, d, true, config.seed);
            std::vector<Real> u0;
            detail::params_to_u(ir.phi, ir.theta,
                                config.enforce_stationarity,
                                config.enforce_invertibility, u0);
            const Real sig0 = std::sqrt(std::max(ir.sigma2, 1e-8));
            u0.push_back(sig0);
            if (inc_const) {
                // d≥1: drift ≈ mean(Δ^d y); d=0: 截距 ≈ μ(1−Σφ)
                Real c0 = mean_z;
                for (Real phv : ir.phi) c0 -= mean_z * phv;
                u0.push_back(std::isfinite(c0) ? c0 : 0.0);
            }
            if (u0.size() == n_u) starts.push_back(std::move(u0));
        } catch (const std::exception&) {
            // innovations 失败 (非平稳等) → 跳过该起始
        }
    }
    if (config.n_starting_points >= 2) {
        // start 2: arima CSS 解
        try {
            arima::ArimaSpec spec_css = spec;
            spec_css.include_drift = inc_const;
            const arima::ArimaResult ar = arima::arima_fit(
                data, spec_css, arima::ArimaMethod::CSS, arima::ArimaConfig{});
            std::vector<Real> u0;
            detail::params_to_u(ar.params.phi, ar.params.theta,
                                config.enforce_stationarity,
                                config.enforce_invertibility, u0);
            const Real sig0 = std::sqrt(std::max(ar.params.sigma2, 1e-8));
            u0.push_back(sig0);
            if (inc_const) {
                // CSS drift 直接可用; 退化时回退均值推导 c = μ(1−Σφ)
                Real c0 = ar.params.drift;
                if (!std::isfinite(c0)) {
                    c0 = mean_z;
                    for (Real phv : ar.params.phi) c0 -= mean_z * phv;
                }
                u0.push_back(std::isfinite(c0) ? c0 : 0.0);
            }
            if (u0.size() == n_u) starts.push_back(std::move(u0));
        } catch (const std::exception&) {
            // CSS 失败 → 跳过
        }
    }
    if (starts.empty()) {
        // 兜底: 零起始 (statsmodels 非平稳置零语义)
        starts.emplace_back(n_u, 0.0);
        starts[0][idx_sigma] = 1.0;
    }

    // ---- SLSQP 多起始优化 (K3: bounds 防护发散, 无约束项) ----
    SLSQP::Config opt = config.optimizer_config;
    std::vector<Bounds> bnds(n_u);
    for (Size i = 0; i < n_u; ++i) {
        if (i == idx_sigma)
            bnds[i] = {1e-4, 1e3};  // u_σ → σ² ∈ [1e-8, 1e6]
        else if (inc_const && i == idx_const)
            bnds[i] = {-1e3, 1e3};
        else
            bnds[i] = {-10.0, 10.0};  // Monahan 域防护
    }

    ArimaSSResult best;
    best.loglik = -std::numeric_limits<Real>::infinity();
    std::vector<Real> best_u;
    for (const auto& u0 : starts) {
        const OptimizationResult r =
            SLSQP::minimize(neg_ll, u0, bnds, {}, {}, opt);
        if (!std::isfinite(r.fx)) continue;
        if (r.fx < -best.loglik) {  // neg_ll 最小 ⇔ loglik 最大
            best.loglik = -r.fx;
            best_u = r.x;
            best.converged = r.converged;
            best.n_iterations = r.n_iterations;
        }
    }

    if (best_u.empty()) {
        best.loglik = 0.0;
        best.message = "arima_statespace_mle: all starts failed";
        return best;
    }

    // ---- 最优解回显 (重跑 KF: loglik + v_t + burn 口径) ----
    const auto up = detail::unpack_params(
        best_u, p, q, config.enforce_stationarity,
        config.enforce_invertibility, inc_const);
    const StateSpaceModel m_best = make_arima_statespace(
        up.phi, up.theta, d, up.sigma2, inc_const,
        config.simple_differencing, up.constant);
    const KalmanResult kr =
        kalman_filter(y_use, m_best, burn_of(up.phi));
    if (!kr.ok || !std::isfinite(kr.loglik)) {
        // 优化器返回解与守卫冲突 (边缘参数) — 保留最优 fx 回显
        best.message = "arima_statespace_mle: final KF guard tripped";
        return best;
    }

    best.phi = up.phi;
    best.theta = up.theta;
    best.constant = up.constant;
    best.sigma2 = up.sigma2;
    best.loglik = kr.loglik;
    best.loglikelihood_burn = burn_of(up.phi);
    best.innovations = kr.innovations;  // v_t (KF 全长, 诊断)

    // AIC (statsmodels 对齐): −2·llf + 2·k, k = p+q+1+const (σ² 计入, K4)
    const Real k_par = static_cast<Real>(n_u);
    best.aic = -2.0 * best.loglik + 2.0 * k_par;

    // ---- K1 常驻互检: innovations 同参数精确似然 ----
    // 口径 (verify_drift_init [K1] + gen_phase8_ss_inc [I] 冻结):
    //   中心化 = 模型隐含均值 μ = c/(1−Σφ) (d≥1: Δy − drift 同款),
    //   非样本均值 — 同参数严格等价:
    //   d=0/simple burn=0: 全样本恒等 (gap ~1e-8 浮点尾差);
    //   d≥1 embedded burn=d: 时间平移等价 (gap O(1/κ) ≈ 1e-6);
    //   AR 非平稳 (全 diffuse burn=d+m₁) 或 μ 不可识别 (1−Σφ ≈ 0):
    //   样本数不同, SS5 禁混比 → 跳过 + message 声明
    const bool nonstat =
        p > 0 && detail::ar_spectral_radius(up.phi) >= 1.0;
    Real sphi = 0.0;
    for (Real phv : up.phi) sphi += phv;
    const bool mu_ok = (1.0 - sphi) > 1e-12;
    if (config.cross_check_innovations && !nonstat &&
        (!inc_const || mu_ok)) {
        std::vector<Real> z_cc = z;  // Δ^d y (simple 模式 = y_use)
        if (inc_const) {
            const Real mu = up.constant / (1.0 - sphi);
            for (Real& v : z_cc) v -= mu;
        }
        try {
            best.innovations_loglik = detail::innovations_loglik_at(
                z_cc, up.phi, up.theta, up.sigma2);
            best.innovations_gap =
                std::abs(best.loglik - best.innovations_loglik);
            if (best.innovations_gap > 1e-5)
                best.message = "cross-check WARNING: innovations gap " +
                               std::to_string(best.innovations_gap) +
                               " > 1e-5";
        } catch (const std::exception&) {
            best.message = "cross-check skipped: innovations pieces failed";
        }
    } else if (config.cross_check_innovations) {
        best.message = nonstat
            ? "cross-check skipped: nonstationary (full diffuse, "
              "burn=d+m1 vs innovations full sample — SS5)"
            : "cross-check skipped: mean unidentified (1-sum(phi)~0)";
    }

    return best;
}

}  // namespace statespace
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub