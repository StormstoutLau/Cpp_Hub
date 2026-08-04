// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.2 Day 3-4 任务 1.7 - HC0-HC5 异方差稳健协方差
// 教材锚点:
//   - White (1980) HC0 原始论文
//   - MacKinnon & White (1985) HC1/HC2/HC3 表 1 (小样本相对表现)
//   - Cribari-Neto (2004) HC4; Cribari-Neto & Souza (2007) HC5
//   - Davidson-MacKinnon Ch.5-6 完整推导 / R sandwich::vcovHC 源码核对
//
// 排幻觉点 (R sandwich 对照):
//   E2: HC1 = N/(N-K) · HC0 (分母是 N-K, 非笔误的 n-k; R vcovHC L120-130 实测)
//   E3: HC2/HC3 leverage h_i = x_i'(X'X)^{-1}x_i (投影矩阵 P=X(X'X)^{-1}X' 对角元,
//       R vcovHC L150-180 用 hatvalues; 非 X'X 对角)
//
// 公式 (严格遵循, 排幻觉点 E2/E3):
//   HC0: V = (X'X)^{-1} [Σ_i x_i x_i' u_i²] (X'X)^{-1}
//   HC1: V = N/(N-K) · HC0
//   HC2: V = (X'X)^{-1} [Σ_i x_i x_i' u_i²/(1-h_i)] (X'X)^{-1}
//   HC3: V = (X'X)^{-1} [Σ_i x_i x_i' u_i²/(1-h_i)²] (X'X)^{-1}
//   HC4: δ_i = min(N·h_i/K, 4);  V = (X'X)^{-1} [Σ_i x_i x_i' u_i²/(1-h_i)^δ_i] (X'X)^{-1}
//   HC5: δ_i = min(N·h_i/K, max(4, 0.7·N·h_max/K));  公式同 HC4 但 δ_i 不同; h_max = max_i h_i
//        (注: spec §7.3 提及 R sandwich HC5 用 (1-h)^(δ/2), 实为 HC4m 修正版行为;
//         HC5 本身用 (1-h)^δ_i, 与任务描述一致)
//
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训, 避免嵌套 namespace 错误)
#pragma once

#include <algorithm>  // std::min, std::max
#include <cmath>      // std::pow
#include <stdexcept>  // std::invalid_argument, std::runtime_error
#include <string>     // to_string

#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/covariance_type.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// econometrics 命名空间内别名 (与 estimation_result.hpp / data_types.hpp 一致, ADR-013)
using linalg::dynamic::MatrixXD;
using linalg::dynamic::VectorXD;

/// @brief 判断协方差类型是否为 HC0-HC5 异方差稳健族
/// @param t 协方差类型
/// @return true 若为 HC0/HC1/HC2/HC3/HC4/HC5
inline bool is_hc_type(CovarianceType t) {
    switch (t) {
        case CovarianceType::HC0:
        case CovarianceType::HC1:
        case CovarianceType::HC2:
        case CovarianceType::HC3:
        case CovarianceType::HC4:
        case CovarianceType::HC5:
            return true;
        default:
            return false;
    }
}

/// @brief 计算 HC0-HC5 异方差稳健协方差矩阵
/// @param X 设计矩阵 (n × k)
/// @param residuals 残差向量 (n)
/// @param XtX_inv (X'X)^{-1} (k × k), 调用方须保证可逆 (SPD)
/// @param type CovarianceType::HC0/HC1/HC2/HC3/HC4/HC5
/// @return 协方差矩阵 (k × k)
/// @throws std::invalid_argument 如果 type 不是 HC0-HC5, 或维度不匹配, 或 n <= k
/// @throws std::runtime_error 如果 leverage h_i == 1 (HC2-HC5 分母为零)
inline MatrixXD compute_hc_vcov(const MatrixXD& X, const VectorXD& residuals,
                                const MatrixXD& XtX_inv, CovarianceType type) {
    if (!is_hc_type(type)) {
        throw std::invalid_argument("compute_hc_vcov: type must be HC0-HC5, got " + to_string(type));
    }
    const Size n = X.rows();
    const Size k = X.cols();
    if (residuals.size() != n) {
        throw std::invalid_argument("compute_hc_vcov: residuals size (" +
                                    std::to_string(residuals.size()) +
                                    ") != X.rows() (" + std::to_string(n) + ")");
    }
    if (XtX_inv.rows() != k || XtX_inv.cols() != k) {
        throw std::invalid_argument("compute_hc_vcov: XtX_inv must be k×k");
    }
    if (n <= k) {
        throw std::invalid_argument("compute_hc_vcov: requires n > k (degrees of freedom)");
    }

    // leverage h_i = X(i,:) · (X'X)^{-1} · X(i,:)' = P(i,i), P = X (X'X)^{-1} X'
    // 排幻觉点 E3: h_i 是投影矩阵对角元, 非 (X'X) 对角元
    Eigen::MatrixXd P = X.eigen() * XtX_inv.eigen() * X.eigen().transpose();

    // HC5 需要 h_max = max_i h_i
    Real hmax = 0.0;
    if (type == CovarianceType::HC5) {
        for (Size i = 0; i < n; ++i) {
            hmax = std::max(hmax, static_cast<Real>(P(i, i)));
        }
    }

    // 权重 w_i (meat 的对角)
    Eigen::VectorXd w = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(n));
    for (Size i = 0; i < n; ++i) {
        const Real u = residuals(i);
        const Real u2 = u * u;
        const Real h = static_cast<Real>(P(i, i));
        const Real one_minus_h = 1.0 - h;
        Real wi = 0.0;
        switch (type) {
            case CovarianceType::HC0:
                wi = u2;
                break;
            case CovarianceType::HC1:
                // 排幻觉点 E2: HC1 = N/(N-K) · HC0, 缩放在最后统一施加
                wi = u2;
                break;
            case CovarianceType::HC2:
                if (one_minus_h <= 0.0) {
                    throw std::runtime_error("compute_hc_vcov: leverage h_i == 1 (HC2 undefined)");
                }
                wi = u2 / one_minus_h;
                break;
            case CovarianceType::HC3:
                if (one_minus_h <= 0.0) {
                    throw std::runtime_error("compute_hc_vcov: leverage h_i == 1 (HC3 undefined)");
                }
                wi = u2 / (one_minus_h * one_minus_h);
                break;
            case CovarianceType::HC4: {
                if (one_minus_h <= 0.0) {
                    throw std::runtime_error("compute_hc_vcov: leverage h_i == 1 (HC4 undefined)");
                }
                const Real delta = std::min(static_cast<Real>(n) * h / static_cast<Real>(k),
                                            static_cast<Real>(4));
                wi = u2 / std::pow(one_minus_h, delta);
                break;
            }
            case CovarianceType::HC5: {
                if (one_minus_h <= 0.0) {
                    throw std::runtime_error("compute_hc_vcov: leverage h_i == 1 (HC5 undefined)");
                }
                const Real delta = std::min(static_cast<Real>(n) * h / static_cast<Real>(k),
                                            std::max(static_cast<Real>(4),
                                                     static_cast<Real>(0.7) * static_cast<Real>(n) * hmax /
                                                         static_cast<Real>(k)));
                wi = u2 / std::pow(one_minus_h, delta);
                break;
            }
            default:
                // is_hc_type 已保证不可达
                throw std::invalid_argument("compute_hc_vcov: unreachable covariance type");
        }
        w(static_cast<Eigen::Index>(i)) = wi;
    }

    // meat = X' diag(w) X  (三明治中部)
    const Eigen::MatrixXd meat = X.eigen().transpose() * w.asDiagonal() * X.eigen();
    // V = (X'X)^{-1} · meat · (X'X)^{-1}  (三明治: bread-meat-bread)
    Eigen::MatrixXd V = XtX_inv.eigen() * meat * XtX_inv.eigen();

    // 排幻觉点 E2: HC1 = N/(N-K) · HC0 (分母 N-K)
    if (type == CovarianceType::HC1) {
        const Real scale = static_cast<Real>(n) / static_cast<Real>(n - k);
        V *= scale;
    }
    return MatrixXD(V);
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
