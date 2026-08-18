// =============================================================================
// irf.hpp - VAR 脉冲响应: Φ 递推 + Ψ = Φ·P 正交化 + bootstrap 带 (spec §4.2)
//
// Phase 7C v1.7 M2 (PHASE7C_SPEC.md v1.2 §4.2, V3/V13)
//
// 公式 (Lütkepohl 2005 Ch.2):
//   Φ_0 = I;  Φ_l = Σ_{i=1..min(l,p)} Φ_{l−i}·A_i   (伴随递推)
//   Θ_h[i,j]: 行 = 响应变量 i, 列 = 冲击变量 j (V3: 转置即全错)
//   正交化: Ψ_h = Φ_h·P, P = 下三角 Cholesky (V2)
//
// 一手源码语义:
//   - statsmodels _chol_sigma_u = np.linalg.cholesky(sigma_u) ← 用
//     df 修正版 Σ (÷(T−Kp−k)), 非 MLE 版 (V2 对齐细节)
//   - R/MATLAB chol 返回上三角: t(chol(Σ)) == numpy chol (V2 等价断言)
//   - P 注入 (决策 11): VARSpec.identification_P 非空时优先, 仅做
//     下三角性校验 (P·P'≈Σ 责任在调用方)
//   - V13: bootstrap 带与 δ 法带不同源, 容差断言仅适用点估计;
//     带实现 = 中心化残差移动块 bootstrap (决策 14 思想), 块长默认
//     ⌈(T_eff/100)^{2/9}·8⌉ Politis-White 风格
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cpphub/core/rng.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/var/multivariate_data.hpp"
#include "cpphub/timeseries/var/var_model.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace var {

/// IRF 结果 (spec §4.2)
struct IRFResult {
    std::vector<Eigen::MatrixXd> theta;  ///< Ψ_0..Ψ_{H−1} (正交化, V3 行=响应)
    bool has_bands = false;              ///< bootstrap 带是否计算 (V13)
    std::vector<Eigen::MatrixXd> irf_lower, irf_upper;  ///< 2.5%/97.5% 分位
    std::vector<Eigen::MatrixXd> phi;    ///< 未正交化 Φ_0..Φ_{H−1} (附)
    Eigen::MatrixXd P;                   ///< 正交化矩阵回显 (V2)
};

namespace detail {

/// 伴随系数 A_1..A_p 从 VARResult 提取 (K×K each)
inline std::vector<Eigen::MatrixXd> var_companion_blocks(const VARResult& fit) {
    const Size K = fit.coefficients.rows();
    const Size p = fit.lag;
    const Size kt = fit.coefficients.cols() - K * p;
    std::vector<Eigen::MatrixXd> A(p, Eigen::MatrixXd::Zero(K, K));
    for (Size l = 0; l < p; ++l) {
        for (Size i = 0; i < K; ++i) {
            for (Size j = 0; j < K; ++j) A[l](i, j) = fit.coefficients(i, kt + l * K + j);
        }
    }
    return A;
}

/// Φ_0..Φ_{H−1} 递推
inline std::vector<Eigen::MatrixXd> phi_recursion(
    const std::vector<Eigen::MatrixXd>& A, Size H) {
    const Size K = A.empty() ? 0 : A[0].rows();
    const Size p = A.size();
    std::vector<Eigen::MatrixXd> phi(H);
    for (Size h = 0; h < H; ++h) {
        phi[h] = Eigen::MatrixXd::Zero(K, K);
        if (h == 0) {
            phi[0].setIdentity();
            continue;
        }
        phi[h].setZero();
        for (Size i = 1; i <= std::min<Size>(h, p); ++i) {
            phi[h] += phi[h - i] * A[i - 1];
        }
    }
    return phi;
}

/// P 矩阵 (V2): 注入优先 (校验下三角), 否则 LLT(df 修正 Σ).matrixL()
inline Eigen::MatrixXd orthogonal_P(const VARResult& fit,
                                    const Eigen::MatrixXd& injected) {
    if (injected.size() > 0) {
        if (injected.rows() != injected.cols() ||
            injected.rows() != fit.sigma_u.rows()) {
            throw std::invalid_argument("var_irf: identification_P size mismatch");
        }
        const Size n = injected.rows();
        for (Size i = 0; i < n; ++i) {
            for (Size j = i + 1; j < n; ++j) {
                if (injected(i, j) != 0.0) {
                    throw std::invalid_argument(
                        "var_irf: identification_P must be lower triangular");
                }
            }
        }
        return injected;
    }
    Eigen::LLT<Eigen::MatrixXd> llt(fit.sigma_u);
    if (llt.info() != Eigen::Success) {
        throw std::invalid_argument("var_irf: sigma_u not PD for Cholesky");
    }
    return llt.matrixL();
}

}  // namespace detail

/// VAR 脉冲响应 (V3). bootstrap=true 时计算 95% 移动块 bootstrap 带 (V13).
inline IRFResult var_irf(const VARResult& fit, Size horizon = 10,
                         bool bootstrap = false, Size n_boot = 1000,
                         Size seed = 42) {
    if (horizon == 0) throw std::invalid_argument("var_irf: horizon must be >= 1");
    if (fit.coefficients.size() == 0) {
        throw std::invalid_argument("var_irf: empty fit");
    }
    const auto A = detail::var_companion_blocks(fit);
    const auto phi = detail::phi_recursion(A, horizon);

    const Eigen::MatrixXd P =
        detail::orthogonal_P(fit, Eigen::MatrixXd{});  // irf 层默认 Cholesky
    std::vector<Eigen::MatrixXd> psi(horizon);
    for (Size h = 0; h < horizon; ++h) psi[h] = phi[h] * P;

    IRFResult res;
    res.theta = std::move(psi);
    res.phi = phi;
    res.P = P;

    if (bootstrap) {
        if (n_boot == 0) throw std::invalid_argument("var_irf: n_boot must be >= 1");
        const Size K = fit.coefficients.rows();
        const Size p = fit.lag;
        const Size T_eff = fit.n_obs_used;
        const Size kt = fit.coefficients.cols() - K * p;
        const Eigen::MatrixXd R = fit.residuals;  // T_eff×K
        // 中心化残差
        Eigen::MatrixXd Rc = R.rowwise() - R.colwise().mean();
        // 块长 (Politis-White 风格简化): ⌈8·(T/100)^{2/9}⌉
        const Size block_len = std::max<Size>(
            1, static_cast<Size>(std::ceil(
                   8.0 * std::pow(static_cast<Real>(T_eff) / 100.0, 2.0 / 9.0))));

        // 重构原始 y 需要: fit 未保留原始数据 — 用残差 + 拟合值重建不可行
        // (无原数据); 改为在 (Y, Z) 空间重抽配对块: 对回归观测行做块重抽,
        // 重建 Y* = Z·B + e*, 重估 B* 后递推 IRF (与重抽 y* 渐近等价,
        // 固定 Z 的 wild-paired 变体; V13: 带仅描述不确定性, 不进容差断言)
        // → 需要原始数据: 由 residuals 与 coefficients 不可逆推 Z — 因此
        // bootstrap 需要原始 MultivariateTSData, 见 var_irf_bootstrap(带数据重载)

        std::vector<Eigen::MatrixXd> lower(horizon), upper(horizon);
        for (Size h = 0; h < horizon; ++h) {
            lower[h] = Eigen::MatrixXd::Constant(K, K,
                             std::numeric_limits<double>::quiet_NaN());
            upper[h] = lower[h];
        }
        res.has_bands = false;
        res.irf_lower = std::move(lower);
        res.irf_upper = std::move(upper);
        (void)seed; (void)block_len; (void)Rc;
        // 说明: 无原数据的 bootstrap 由 var_irf_bootstrap(data, ...) 提供
    }
    return res;
}

/// 带 bootstrap 置信带的 IRF (重载: 提供原始数据).
/// 中心化残差移动块重抽 → y* 递推 (初始 p 期取真实观测) → 重估 → IRF
/// 95% percentile 带 (V13: 仅不确定性描述, 容差断言不适用).
inline IRFResult var_irf_bootstrap(const MultivariateTSData& data,
                                   const VARResult& fit, Size horizon = 10,
                                   Size n_boot = 1000, Size seed = 42) {
    IRFResult base = var_irf(fit, horizon);
    const Size K = fit.coefficients.rows();
    const Size p = fit.lag;
    const Eigen::MatrixXd y = data.matrix();
    const Size T = y.rows();
    const Eigen::MatrixXd R = fit.residuals;
    const Size T_eff = R.rows();
    Eigen::MatrixXd Rc = R.rowwise() - R.colwise().mean();
    const Size block_len = std::max<Size>(
        1, static_cast<Size>(std::ceil(
               8.0 * std::pow(static_cast<Real>(T_eff) / 100.0, 2.0 / 9.0))));

    Philox4x64 rng(static_cast<uint64_t>(seed), 0x49524642ULL);
    std::vector<std::vector<Eigen::MatrixXd>> samples(
        horizon, std::vector<Eigen::MatrixXd>(n_boot));
    Eigen::MatrixXd ystar = Eigen::MatrixXd::Zero(T, K);
    Eigen::MatrixXd eall = Eigen::MatrixXd::Zero(T_eff, K);
    for (Size b = 0; b < n_boot; ++b) {
        // 块重抽 (环形跨块允许)
        Size filled = 0;
        while (filled < T_eff) {
            const Size start = static_cast<Size>(rng() % T_eff);
            for (Size l = 0; l < block_len && filled < T_eff; ++l) {
                eall.row(filled++) = Rc.row((start + l) % T_eff);
            }
        }
        ystar.topRows(p) = y.topRows(p);
        // 递推 y*: 使用估计系数 (固定 B, residual bootstrap);
        // 确定性列 [c: 截距, ct: 截距+趋势 (时间值 = i+1, i 为回归行号)]
        const auto A = detail::var_companion_blocks(fit);
        const Size kt = fit.coefficients.cols() - K * p;
        for (Size t = p; t < T; ++t) {
            Eigen::VectorXd acc = Eigen::VectorXd::Zero(K);
            if (kt >= 1) acc += fit.coefficients.col(0);
            if (kt >= 2) acc += fit.coefficients.col(1) * static_cast<Real>(t - p + 1);
            for (Size l = 0; l < p; ++l) acc += A[l] * ystar.row(t - l - 1).transpose();
            ystar.row(t) = acc + eall.row(t - p);
        }
        // 重估
        MultivariateTSData dstar;
        dstar.columns.resize(K);
        for (Size j = 0; j < K; ++j) {
            dstar.columns[j].resize(T);
            for (Size t = 0; t < T; ++t) dstar.columns[j][t] = ystar(t, j);
        }
        VARSpec sp;
        sp.lag = p;
        sp.trend = fit.trend;
        VARResult fstar = var_fit(dstar, sp);
        const auto Ab = detail::var_companion_blocks(fstar);
        const auto phib = detail::phi_recursion(Ab, horizon);
        Eigen::LLT<Eigen::MatrixXd> llt(fstar.sigma_u);
        if (llt.info() != Eigen::Success) continue;
        const Eigen::MatrixXd Pb = llt.matrixL();
        for (Size h = 0; h < horizon; ++h) samples[h][b] = phib[h] * Pb;
    }
    // percentile 2.5/97.5 (逐元素)
    base.has_bands = true;
    base.irf_lower.resize(horizon);
    base.irf_upper.resize(horizon);
    std::vector<Real> vals(n_boot);
    for (Size h = 0; h < horizon; ++h) {
        base.irf_lower[h] = Eigen::MatrixXd::Zero(K, K);
        base.irf_upper[h] = Eigen::MatrixXd::Zero(K, K);
        for (Size i = 0; i < K; ++i) {
            for (Size j = 0; j < K; ++j) {
                for (Size b = 0; b < n_boot; ++b) vals[b] = samples[h][b](i, j);
                std::sort(vals.begin(), vals.end());
                const Real q025 = vals[static_cast<Size>(0.025 * (n_boot - 1))];
                const Real q975 = vals[static_cast<Size>(0.975 * (n_boot - 1))];
                base.irf_lower[h](i, j) = q025;
                base.irf_upper[h](i, j) = q975;
            }
        }
    }
    return base;
}

}  // namespace var
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
