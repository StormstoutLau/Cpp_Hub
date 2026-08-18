// =============================================================================
// fevd.hpp - FEVD 双轨: Cholesky 正交化 + GFEVD (DY/PS 框架) (spec §4.2)
//
// Phase 7C v1.7 M2 (PHASE7C_SPEC.md v1.2 §4.2, 决策 12/15, V7/V8/V12)
//
// 公式:
//   Cholesky 轨 (行和精确 = 1, V7):
//     FEVD_{ij}(h) = Σ_{l<h} Ψ_l[i,j]² / Σ_{j'} Σ_{l<h} Ψ_l[i,j']²
//     (Ψ = Φ·P, P = 下三角 Cholesky; 与 statsmodels FEVD L2370-2399
//      的 cumsum(Ψ²)/mse_diag 恒等 — 正交化情形行和 = 响应方差)
//   GFEVD DY 2012 框架 (默认, 溢出指数专用; R Spillover g.fevd 主基准):
//     num[i,j] = σ_jj⁻¹ · Σ_{l<h} [(Φ_l Σ)[i,j]]²     ← 分子系数 σ_jj⁻¹
//     den[i]   = Σ_{l<h} (Φ_l Σ Φ_l')[i,i]
//     θ̃_ij = (num/den) 行归一化 (V7 行和 = 1)
//     ⚠️ Σ 用 df 修正版 (SSR/(T−Kp−k)), 对齐 g.fevd 的
//        summary(x)$covres — 非 MLE 版 (一手源码 L54 实录)
//   GFEVD PS 1998 框架 (可选输出, 文献对照):
//     分子系数 = σ_ii⁻¹ (响应变量方差), 配未归一分母, 行和 ≠ 1 (V8:
//     两框架归一化后数值不同, 除非各 σ 相等; API 注明框架且不混用)
//   V12: 不稳定 VAR (max|eig 伴随| ≥ 1) 前置拦截 FEVD
//   V11: GIRF 隐含假设 (Kim 2013) — Σ 直接进入分子, 无需正交化
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/var/irf.hpp"
#include "cpphub/timeseries/var/multivariate_data.hpp"
#include "cpphub/timeseries/var/var_model.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace var {

/// FEVD 框架 (spec §4.2)
enum class FevdFramework { Cholesky, GeneralizedDY, GeneralizedPS };

/// FEVD 结果 (spec §4.2)
struct FEVDResult {
    Eigen::MatrixXd fevd;       ///< K×K (Cholesky/DY: 行和=1; PS: 行和≠1)
    FevdFramework framework{};
    Size horizon = 0;
};

/// FEVD 计算 (V7/V8/V12)
inline FEVDResult var_fevd(const VARResult& fit, Size horizon = 10,
                           FevdFramework fw = FevdFramework::GeneralizedDY) {
    if (horizon == 0) throw std::invalid_argument("var_fevd: horizon must be >= 1");
    if (fit.coefficients.size() == 0) {
        throw std::invalid_argument("var_fevd: empty fit");
    }
    // V12: 不稳定 VAR 拦截 (FEVD 无意义)
    if (!fit.is_strictly_stationary) {
        throw std::invalid_argument(
            "var_fevd: VAR not strictly stationary (companion max |eig| >= 1); "
            "FEVD undefined (V12)");
    }
    const Size K = fit.coefficients.rows();

    FEVDResult res;
    res.framework = fw;
    res.horizon = horizon;

    if (fw == FevdFramework::Cholesky) {
        const Eigen::MatrixXd P =
            detail::orthogonal_P(fit, Eigen::MatrixXd{});
        Eigen::MatrixXd num = Eigen::MatrixXd::Zero(K, K);
        const auto A = detail::var_companion_blocks(fit);
        const auto phi = detail::phi_recursion(A, horizon);
        for (Size l = 0; l < horizon; ++l) {
            const Eigen::MatrixXd psi = phi[l] * P;
            num += psi.cwiseProduct(psi);
        }
        // 行归一化 (行和=1, V7)
        for (Size i = 0; i < K; ++i) {
            const Real rs = num.row(i).sum();
            if (rs <= 0) throw std::invalid_argument("var_fevd: zero row sum");
            num.row(i) /= rs;
        }
        res.fevd = num;
        return res;
    }

    // GFEVD (V11: Σ 直接进入, 无正交化; Σ = df 修正版, 对齐 g.fevd covres)
    const Eigen::MatrixXd& S = fit.sigma_u;
    const auto A = detail::var_companion_blocks(fit);
    const auto phi = detail::phi_recursion(A, horizon);
    Eigen::MatrixXd num = Eigen::MatrixXd::Zero(K, K);
    Eigen::VectorXd den = Eigen::VectorXd::Zero(K);
    for (Size l = 0; l < horizon; ++l) {
        const Eigen::MatrixXd PS = phi[l] * S;
        num += PS.cwiseProduct(PS);
        den += (PS * phi[l].transpose()).diagonal();
    }
    Eigen::MatrixXd theta = Eigen::MatrixXd::Zero(K, K);
    for (Size i = 0; i < K; ++i) {
        if (den(i) <= 0) throw std::invalid_argument("var_fevd: zero denominator");
        for (Size j = 0; j < K; ++j) {
            if (fw == FevdFramework::GeneralizedDY) {
                // DY: σ_jj⁻¹ 分子系数
                theta(i, j) = num(i, j) / S(j, j) / den(i);
            } else {
                // PS: σ_ii⁻¹ 分子系数, 未归一 (行和 ≠ 1)
                theta(i, j) = num(i, j) / S(i, i) / den(i);
            }
        }
    }
    if (fw == FevdFramework::GeneralizedDY) {
        for (Size i = 0; i < K; ++i) {
            const Real rs = theta.row(i).sum();
            if (rs <= 0) throw std::invalid_argument("var_fevd: zero row sum");
            theta.row(i) /= rs;
        }
    }
    res.fevd = theta;
    return res;
}

}  // namespace var
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
