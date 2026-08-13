// =============================================================================
// ols_simple.hpp - 轻量级 OLS 估计 (std::vector + Gauss-Jordan)
//
// Phase 7A Wave 0: detail/ 公共基础设施 (ADR-015 方案 B)
//
// 设计依据: har_model.hpp::ols_estimate (L186-307, 已验证) 的 Gauss-Jordan 实现模式
// 与 ols_estimate 的有意差异 (ADR-015 决策点 5):
//   - 不自动添加常数列 (由调用方决定)
//   - 不计算 adj_r_squared / llh (诊断辅助回归不需要)
//   - 仅返回系数 + fitted + residuals + r_squared
//
// 规模适用范围: N=百级到千级, K<10 (ADR015 §7.4 H8 修正)
// 生产级 OLS 估计 (大 N, 性能敏感) 用 econometrics/estimation/ols.hpp (Eigen3)
//
// ADR-015 方案 B: 仅依赖 core/, 不依赖 linalg_dynamic.hpp (Eigen3)
// =============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {
namespace detail {

// 轻量级 OLS 估计: y = X beta + epsilon (含常数项由调用方决定)
// 实现方式: std::vector + Gauss-Jordan 消元 (partial pivoting)
// 参考: har_model.hpp::ols_estimate (L186-307)
//
// @param y 因变量 (长度 N)
// @param X 设计矩阵 (N×K, 不含常数列; 调用方若需常数项请自行添加一列 1.0)
// @param fitted_values 输出拟合值
// @param residuals 输出残差
// @param r_squared 输出 R²
// @return OLS 系数 (长度 K)
inline std::vector<Real> ols_simple(
    const std::vector<Real>& y,
    const std::vector<std::vector<Real>>& X,
    std::vector<Real>& fitted_values,
    std::vector<Real>& residuals,
    Real& r_squared) {

    const Size n = y.size();
    if (n == 0) {
        throw std::invalid_argument("ols_simple: empty y");
    }
    const Size k = X.empty() ? 0 : X[0].size();
    const Size p = k;  // 不自动加常数项 (与 ols_estimate 的有意差异)

    if (p == 0) {
        throw std::invalid_argument("ols_simple: empty design matrix");
    }
    if (n < p) {
        throw std::invalid_argument("ols_simple: insufficient observations");
    }

    // 正规方程: (X^T X) beta = X^T y
    std::vector<std::vector<Real>> XtX(p, std::vector<Real>(p, 0.0));
    std::vector<Real> Xty(p, 0.0);
    for (Size i = 0; i < p; ++i) {
        for (Size j = 0; j < p; ++j) {
            Real s = 0.0;
            for (Size t = 0; t < n; ++t) {
                s += X[t][i] * X[t][j];
            }
            XtX[i][j] = s;
        }
        Real s = 0.0;
        for (Size t = 0; t < n; ++t) {
            s += X[t][i] * y[t];
        }
        Xty[i] = s;
    }

    // Gauss-Jordan 消元求解 (partial pivoting)
    std::vector<std::vector<Real>> A = XtX;
    std::vector<Real> b = Xty;
    for (Size col = 0; col < p; ++col) {
        // Partial pivot
        Size piv = col;
        Real maxv = std::fabs(A[col][col]);
        for (Size r = col + 1; r < p; ++r) {
            Real v = std::fabs(A[r][col]);
            if (v > maxv) { maxv = v; piv = r; }
        }
        if (maxv < 1e-15) {
            throw std::runtime_error("ols_simple: singular design matrix");
        }
        if (piv != col) {
            std::swap(A[piv], A[col]);
            std::swap(b[piv], b[col]);
        }
        Real akk = A[col][col];
        for (Size r = col + 1; r < p; ++r) {
            Real f = A[r][col] / akk;
            if (f == 0.0) continue;
            for (Size c = col; c < p; ++c) A[r][c] -= f * A[col][c];
            b[r] -= f * b[col];
        }
    }
    // Back substitution
    std::vector<Real> beta(p);
    for (Size i = p; i-- > 0;) {
        Real s = b[i];
        for (Size j = i + 1; j < p; ++j) s -= A[i][j] * beta[j];
        beta[i] = s / A[i][i];
    }

    // 计算 fitted values 和 residuals
    fitted_values.resize(n);
    residuals.resize(n);
    for (Size t = 0; t < n; ++t) {
        Real pred = 0.0;
        for (Size j = 0; j < k; ++j) {
            pred += beta[j] * X[t][j];
        }
        fitted_values[t] = pred;
        residuals[t] = y[t] - pred;
    }

    // R^2
    Real y_mean = 0.0;
    for (Size t = 0; t < n; ++t) y_mean += y[t];
    y_mean /= static_cast<Real>(n);
    Real ss_tot = 0.0, ss_res = 0.0;
    for (Size t = 0; t < n; ++t) {
        Real dy = y[t] - y_mean;
        ss_tot += dy * dy;
        ss_res += residuals[t] * residuals[t];
    }
    r_squared = (ss_tot > 1e-300) ? (1.0 - ss_res / ss_tot) : 0.0;

    return beta;
}

}  // namespace detail
}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
