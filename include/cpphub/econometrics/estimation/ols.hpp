// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.2 Day 3-4 任务 1.6 - OLSEstimator (排幻觉点 E1)
// 经典普通最小二乘 (OLS) 估计器 + Classical/HC0-HC5 协方差支持
//
// 教材锚点: Davidson-MacKinnon Ch.2-3, Greene Ch.3-4, R lm() + sandwich::vcovHC
//
// 排幻觉点 E1 (对照 R lm()):
//   R lm() 默认在 X 中隐式添加截距列; C++ OLSEstimator 不自动添加截距,
//   由用户在 X 中显式构造 (若需截距, 第一列应全为 1). 若 X 第一列全为 1 视为含截距.
//
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训, 避免嵌套 namespace 错误)
#pragma once

#include <cmath>       // std::log, std::sqrt, std::lgamma, std::pow, std::abs, std::exp
#include <limits>      // std::numeric_limits
#include <memory>      // std::unique_ptr, std::make_unique
#include <stdexcept>   // std::invalid_argument, std::runtime_error
#include <string>      // std::string, std::to_string
#include <type_traits> // std::decay_t, std::is_same_v

#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/covariance_type.hpp"
#include "cpphub/econometrics/core/data_types.hpp"
#include "cpphub/econometrics/core/estimation_result.hpp"
#include "cpphub/econometrics/core/estimator_base.hpp"
#include "cpphub/econometrics/core/special_functions.hpp"
#include "cpphub/econometrics/inference/hc_standard_errors.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// econometrics 命名空间内别名 (与 data_types.hpp / estimation_result.hpp 一致, ADR-013)
using linalg::dynamic::MatrixXD;
using linalg::dynamic::VectorXD;

// =============================================================================
// detail: t 分布双侧 p 值
//   依赖 special_functions.hpp 中的 kTwoPi / betacf / beta_i (共享定义点)
//   t 双侧 p 值: p = I_x(df/2, 1/2), x = df/(df + t²)
// =============================================================================
namespace detail {

/// @brief t 分布双侧 p 值: P(|T| > |t|)
/// @param t 观测 t 统计量
/// @param df 自由度 (残差自由度 n-k)
/// @return 双侧 p 值 ∈ [0,1]; df <= 0 返回 NaN
inline Real t_two_sided_pvalue(Real t, Real df) {
    if (df <= 0.0) return std::numeric_limits<Real>::quiet_NaN();
    const Real x = df / (df + t * t);
    return beta_i(0.5 * df, 0.5, x);
}

}  // namespace detail

// =============================================================================
// OLSEstimator - 普通最小二乘估计器
//
// 设计:
//   - 不自动添加截距 (排幻觉点 E1), 用户在 X 中显式构造
//   - 支持 CovarianceType::Classical 与 HC0-HC5; 其他类型抛 invalid_argument
//   - (X'X)^{-1} 用 LLT (inverse_symmetric), 奇异 X 抛 runtime_error
//   - 仅支持 CrossSectionData; PanelData/TimeSeriesData 抛 invalid_argument
// =============================================================================
class OLSEstimator : public Estimator {
public:
    OLSEstimator() = default;

    /// @brief 指定协方差类型的构造器
    /// @param cov_type 协方差类型 (Classical 或 HC0-HC5)
    explicit OLSEstimator(CovarianceType cov_type) { cov_type_ = cov_type; }

    /// @brief 执行 OLS 估计
    /// @param data 横截面数据 (仅 CrossSectionData)
    /// @return 估计结果 (系数/标准误/t/p/协方差/R²/对数似然等)
    EstimationResult estimate(const EconData& data) override;

    std::string name() const override { return "OLS"; }
    EstimatorClass estimatorClass() const override { return EstimatorClass::Parametric; }
    std::unique_ptr<Estimator> clone() const override;

    /// @brief 计算拟合值 ŷ = X·β
    VectorXD computeFittedValues(const MatrixXD& X, const VectorXD& beta) const;

    /// @brief 计算残差 u = y - X·β
    VectorXD computeResiduals(const MatrixXD& X, const VectorXD& y, const VectorXD& beta) const;

    /// @brief 计算投影矩阵 P = X·(X'X)^{-1}·X' (n×n, 对称幂等, 对角 = leverage h_i)
    MatrixXD computeProjectionMatrix(const MatrixXD& X, const MatrixXD& XtX_inv) const;

    /// @brief 计算 R² = 1 - SSR/SST
    Real computeRSquared(const VectorXD& y, const VectorXD& fitted, const VectorXD& residuals) const;

    /// @brief 计算调整后 R² = 1 - (1-R²)·(n-1)/(n-k)
    Real computeAdjustedRSquared(Real r_squared, Size n, Size k) const;

    /// @brief 计算联合显著性 F 统计量 = (R²/k) / ((1-R²)/(n-k))
    Real computeFStatistic(Real r_squared, Size n, Size k) const;

private:
    /// @brief 从 EconData variant 提取 X, y (仅 CrossSectionData)
    /// @throws std::invalid_argument 若为 PanelData/TimeSeriesData
    void extractXY(const EconData& data, MatrixXD& X, VectorXD& y) const;
};

// =============================================================================
// 内联实现 (header-only, inline 避免 ODR 冲突)
// =============================================================================

inline std::unique_ptr<Estimator> OLSEstimator::clone() const {
    return std::make_unique<OLSEstimator>(*this);
}

inline void OLSEstimator::extractXY(const EconData& data, MatrixXD& X, VectorXD& y) const {
    std::visit([&](auto&& d) {
        using T = std::decay_t<decltype(d)>;
        if constexpr (std::is_same_v<T, CrossSectionData>) {
            X = d.X;
            y = d.y;
        } else {
            // PanelData / TimeSeriesData: OLS 横截面估计器不支持
            throw std::invalid_argument(
                "OLSEstimator::extractXY: only CrossSectionData is supported "
                "(PanelData/TimeSeriesData require panel/time-series estimators)");
        }
    }, data);
}

inline VectorXD OLSEstimator::computeFittedValues(const MatrixXD& X, const VectorXD& beta) const {
    return VectorXD(X.eigen() * beta.eigen());
}

inline VectorXD OLSEstimator::computeResiduals(const MatrixXD& X, const VectorXD& y,
                                               const VectorXD& beta) const {
    return VectorXD(y.eigen() - X.eigen() * beta.eigen());
}

inline MatrixXD OLSEstimator::computeProjectionMatrix(const MatrixXD& X, const MatrixXD& XtX_inv) const {
    // P = X · (X'X)^{-1} · X' (n×n)
    return MatrixXD(X.eigen() * XtX_inv.eigen() * X.eigen().transpose());
}

inline Real OLSEstimator::computeRSquared(const VectorXD& y, const VectorXD& /*fitted*/,
                                          const VectorXD& residuals) const {
    const Real ssr = residuals.eigen().squaredNorm();
    const Real ybar = y.eigen().mean();
    const Real sst = (y.eigen().array() - ybar).matrix().squaredNorm();
    if (sst <= 0.0) {
        throw std::runtime_error("OLSEstimator::computeRSquared: zero total sum of squares (SST=0)");
    }
    return 1.0 - ssr / sst;
}

inline Real OLSEstimator::computeAdjustedRSquared(Real r_squared, Size n, Size k) const {
    if (n <= k) {
        throw std::invalid_argument("OLSEstimator::computeAdjustedRSquared: requires n > k");
    }
    return 1.0 - (1.0 - r_squared) * static_cast<Real>(n - 1) / static_cast<Real>(n - k);
}

inline Real OLSEstimator::computeFStatistic(Real r_squared, Size n, Size k) const {
    if (n <= k) {
        throw std::invalid_argument("OLSEstimator::computeFStatistic: requires n > k");
    }
    if (r_squared >= 1.0) {
        throw std::runtime_error("OLSEstimator::computeFStatistic: R² >= 1 (perfect fit, F undefined)");
    }
    return (r_squared / static_cast<Real>(k)) /
           ((1.0 - r_squared) / static_cast<Real>(n - k));
}

inline EstimationResult OLSEstimator::estimate(const EconData& data) {
    // 1. 提取 X, y (仅 CrossSectionData)
    MatrixXD X;
    VectorXD y;
    extractXY(data, X, y);

    const Size n = X.rows();
    const Size k = X.cols();
    if (k == 0) {
        throw std::invalid_argument("OLSEstimator::estimate: X has zero columns");
    }
    if (y.size() != n) {
        throw std::invalid_argument("OLSEstimator::estimate: X.rows() != y.size()");
    }
    // n > k (自由度要求); n < k 报错, n == k 自由度为零无法推断
    if (n <= k) {
        throw std::invalid_argument(
            "OLSEstimator::estimate: requires n > k (got n=" + std::to_string(n) +
            ", k=" + std::to_string(k) + ")");
    }

    // 2. (X'X)^{-1} via LLT (SPD 检测); 奇异 X 抛 runtime_error
    //    排幻觉点: OLS 用正规方程 β = (X'X)^{-1} X'y, 不自动添加截距 (E1)
    const MatrixXD XtX = MatrixXD(X.eigen().transpose() * X.eigen());
    const MatrixXD XtX_inv = inverse_symmetric(XtX);  // LLT, 非 SPD 抛 runtime_error
    const VectorXD Xty = VectorXD(X.eigen().transpose() * y.eigen());
    const VectorXD beta = VectorXD(XtX_inv.eigen() * Xty.eigen());

    // 3. 拟合值 / 残差 / R²
    const VectorXD fitted = computeFittedValues(X, beta);
    const VectorXD residuals = computeResiduals(X, y, beta);
    const Real r2 = computeRSquared(y, fitted, residuals);
    const Real adj_r2 = computeAdjustedRSquared(r2, n, k);
    const Size df_resid = n - k;

    // 4. 协方差矩阵 (Classical 或 HC0-HC5)
    MatrixXD vcov;
    const CovarianceType ct = cov_type_;
    if (ct == CovarianceType::Classical) {
        const Real sigma2 = residuals.eigen().squaredNorm() / static_cast<Real>(df_resid);
        vcov = MatrixXD(XtX_inv.eigen() * sigma2);
    } else if (is_hc_type(ct)) {
        vcov = compute_hc_vcov(X, residuals, XtX_inv, ct);
    } else {
        throw std::invalid_argument(
            "OLSEstimator::estimate: unsupported covariance type '" + to_string(ct) +
            "' (OLS supports Classical and HC0-HC5)");
    }

    // 5. 标准误 / t 统计量 / 双侧 p 值
    VectorXD se(k);
    VectorXD tstat(k);
    VectorXD pval(k);
    for (Size i = 0; i < k; ++i) {
        Real v = vcov(i, i);
        if (v < 0.0) v = 0.0;  // 数值保护: 协方差对角不应为负
        const Real s = std::sqrt(v);
        se(i) = s;
        if (s > 0.0) {
            const Real t = beta(i) / s;
            tstat(i) = t;
            pval(i) = detail::t_two_sided_pvalue(t, static_cast<Real>(df_resid));
        } else {
            tstat(i) = std::numeric_limits<Real>::infinity();
            pval(i) = 0.0;
        }
    }

    // 6. 高斯对数似然 (MLE): σ²_MLE = SSR/n
    //    LL = -n/2 · (log(2π) + 1 + log(σ²_MLE))
    const Real ssr = residuals.eigen().squaredNorm();
    const Real sigma2_mle = ssr / static_cast<Real>(n);
    Real loglik = 0.0;
    if (sigma2_mle > 0.0) {
        loglik = -0.5 * static_cast<Real>(n) *
                 (std::log(detail::kTwoPi) + 1.0 + std::log(sigma2_mle));
    } else {
        loglik = std::numeric_limits<Real>::infinity();  // 完美拟合
    }

    // 7. 填充结果
    EstimationResult result;
    result.coefficients = beta;
    result.std_errors = se;
    result.t_statistics = tstat;
    result.p_values = pval;
    result.vcov = vcov;
    result.log_likelihood = loglik;
    result.r_squared = r2;
    result.adj_r_squared = adj_r2;
    result.n_obs = n;
    result.n_params = k;
    result.df_residual = df_resid;
    result.cov_type = ct;
    return result;
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
