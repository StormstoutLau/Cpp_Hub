// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.3 任务 1.9 - HAC 协方差 (Newey-West / Andrews)
// 排幻觉点:
//   E4 (Andrews 自动带宽 vs NW 经验法则): max_lag=0 默认走 NW 经验法则
//      floor(4*(T/100)^(2/9)), 非 Andrews 自动带宽。调用方需要 Andrews 时应先拟合
//      AR(1) 再通过 select_max_lag(..., andrews_optimal=true, rho) 显式传 max_lag。
//   E5 (Bartlett w[l] = 1-l/(L+1), 非 1-l/L): 权重分母为 L+1, 复用已验证的
//      kernel_weights() (R sandwich::kweights 实测一致)。
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训)
#pragma once

#include <stdexcept>
#include <vector>

#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/covariance_type.hpp"
#include "cpphub/econometrics/inference/hac_kernels.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

using linalg::dynamic::MatrixXD;
using linalg::dynamic::VectorXD;

/// @brief 计算 HAC (Heteroskedasticity and Autocorrelation Consistent) 协方差矩阵
/// @param X 设计矩阵 (T × k)
/// @param residuals 残差向量 (T)
/// @param XtX_inv (X'X)^{-1} (k × k)
/// @param kernel HAC 内核类型 (Bartlett/QS/Parzen/TukeyHanning)
/// @param max_lag 最大滞后 L (若为 0, 自动选择: NW 经验法则 floor(4*(T/100)^(2/9)))
/// @param prewhiten 是否预白化 (Andrews-Monahan 1992, 默认 false, 暂未实现)
/// @param ar1_coef AR(1) 系数 rho (预留给 prewhiten 使用, 当前版本未使用)
/// @return HAC 协方差矩阵 (k × k)
/// @throws std::runtime_error 若 prewhiten=true (暂未实现) 或 X 维度与 residuals 不匹配
/// @throws std::invalid_argument 若 max_lag >= T 或 XtX_inv 维度不匹配
///
/// 公式: V_HAC = (X'X)^{-1} * Ω * (X'X)^{-1}
///   Ω = w[0] * Ω_0 + Σ_{l=1}^{L} w[l] * (Ω_l + Ω_l')
///   Ω_0[i][j] = Σ_t X[t][i] * u[t]² * X[t][j]              (contemporaneous)
///   Ω_l[i][j] = Σ_{t=l}^{T-1} X[t-l][i] * u[t-l] * u[t] * X[t][j]   (lag l > 0)
///
/// 排幻觉点 E5: 权重 w[l] = K(l/(L+1)), 由 kernel_weights() 计算 (非 K(l/L))
/// 排幻觉点 E4: max_lag=0 默认走 NW 经验法则, 非 Andrews 自动带宽
inline MatrixXD compute_hac_vcov(const MatrixXD& X, const VectorXD& residuals,
                                  const MatrixXD& XtX_inv, HacKernel kernel,
                                  Size max_lag = 0, bool prewhiten = false,
                                  Real ar1_coef = 0.0) {
    // ar1_coef 预留给 prewhiten 实现, 当前版本未使用 (排幻觉点 E4: 不走 Andrews 自动带宽)
    (void)ar1_coef;

    // prewhiten (Andrews-Monahan 1992) 暂未实现
    if (prewhiten) {
        throw std::runtime_error(
            "compute_hac_vcov: prewhiten (Andrews-Monahan 1992) not yet implemented");
    }

    const Size T = X.rows();
    const Size k = X.cols();

    // 维度校验
    if (T != residuals.size()) {
        throw std::runtime_error(
            "compute_hac_vcov: dimension mismatch (X.rows() != residuals.size())");
    }
    if (XtX_inv.rows() != k || XtX_inv.cols() != k) {
        throw std::invalid_argument(
            "compute_hac_vcov: XtX_inv dimension mismatch (expected k×k)");
    }
    if (T == 0 || k == 0) {
        throw std::invalid_argument("compute_hac_vcov: empty X matrix");
    }

    // 排幻觉点 E4: max_lag=0 默认走 NW 经验法则 (非 Andrews 自动带宽)
    if (max_lag == 0) {
        max_lag = select_max_lag(T, kernel, false, 0.0);
    }

    if (max_lag >= T) {
        throw std::invalid_argument(
            "compute_hac_vcov: max_lag must be < T (number of observations)");
    }

    // 排幻觉点 E5: 权重 w[l] = K(l/(L+1)), 由已验证的 kernel_weights() 计算
    const std::vector<Real> w = kernel_weights(kernel, max_lag);

    // 计算 Ω = w[0]*Ω_0 + Σ_{l=1}^{L} w[l]*(Ω_l + Ω_l')
    // 利用 Ω 的对称性 (Ω_0 对称, Ω_l+Ω_l' 对称), 只算上三角再镜像
    MatrixXD Omega(k, k);
    for (Size i = 0; i < k; ++i) {
        for (Size j = i; j < k; ++j) {
            Real val = 0.0;

            // Ω_0[i][j] = Σ_t X[t][i] * u[t]² * X[t][j]  (w[0] = K(0) = 1 对所有内核)
            for (Size t = 0; t < T; ++t) {
                val += w[0] * X(t, i) * residuals(t) * residuals(t) * X(t, j);
            }

            // Ω_l[i][j] + Ω_l[j][i] = Σ_t u[t-l]*u[t] * (X[t-l][i]*X[t][j] + X[t-l][j]*X[t][i])
            for (Size l = 1; l <= max_lag; ++l) {
                Real cross = 0.0;
                for (Size t = l; t < T; ++t) {
                    const Real uu = residuals(t - l) * residuals(t);
                    cross += uu * (X(t - l, i) * X(t, j) + X(t - l, j) * X(t, i));
                }
                val += w[l] * cross;
            }

            Omega(i, j) = val;
            Omega(j, i) = val;  // Ω 对称
        }
    }

    // V_HAC = (X'X)^{-1} * Ω * (X'X)^{-1}
    return XtX_inv * Omega * XtX_inv;
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
