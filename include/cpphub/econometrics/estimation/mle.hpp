// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.1 任务 2.1-2.2 - MLE/QMLE 估计器 (排幻觉点 E7/E8)
// 经典极大似然估计 (MLE) + 拟极大似然 (QMLE) Sandwich 协方差
//
// 教材锚点:
//   - Greene 8ed Ch.17 (Logit/Probit/Poisson 等离散选择模型)
//   - McCullagh-Nelder 1989 GLM (canonical link + IRLS)
//   - White 1982 (QMLE Sandwich: A^{-1} B A^{-1})
//   - Wooldridge 2010 §12.3 (GLM QMLE 渐近理论)
//
// 算法: Newton-Raphson + 解析梯度/Hessian (对 canonical link 等价于 IRLS)
//   迭代: β_{k+1} = β_k + (X'WX)^{-1} X'(y - μ)
//   其中 μ = g^{-1}(Xβ) 为均值函数, W = diag(w_i) 为 IRLS 权重
//
// 排幻觉点 E7 (R glm IRLS vs C++ Newton-Raphson):
//   R glm() 对 canonical link 用 IRLS (Iteratively Reweighted Least Squares);
//   C++ 用 Newton-Raphson + 解析 Hessian.
//   对 canonical link (Logit/Log/Identity), 两者数学等价:
//     IRLS: β_{k+1} = (X'WX)^{-1} X'Wz, z = Xβ + (y-μ)/w·w = Xβ + (y-μ)·W^{-1}·W → X'(y-μ) + X'WXβ
//     NR:   β_{k+1} = β + (-H)^{-1}·score = β + (X'WX)^{-1}·X'(y-μ)
//     两者化简后完全一致.
//   对 non-canonical link (Probit), NR 与 IRLS 等价但 W 的形式不同.
//
// 排幻觉点 E8 (GLM bread = (X'WX)^{-1}, meat = X' diag(ε²) X):
//   R sandwich::sandwich() 对 GLM 的实现:
//     bread = (X'WX)^{-1}  (Hessian 逆, W 为 IRLS 权重)
//     meat  = X' diag(ε²) X  (OPG 外积, ε = y - μ 为原始残差)
//   注: meat 用的是原始残差 (y - μ), 不是 score 的外积 (但对 canonical link 两者一致,
//       因 score_i = x_i(y_i - μ_i), score_i score_i' = x_i x_i' (y_i-μ_i)²).
//
// 协方差类型 (CovarianceType):
//   - Hessian: V = -H^{-1} = (X'WX)^{-1}  (经典 MLE, 假设模型正确)
//   - OPG (Outer Product of Gradients): V = (Σ g_i g_i')^{-1}  (Berndt-Hall-Hall-Hausman)
//   - Sandwich (QMLE): V = A^{-1} B A^{-1} = (X'WX)^{-1} X' diag(ε²) X (X'WX)^{-1}
//   注: 对正确指定的 MLE, 三者渐近等价; QMLE 仅需均值正确 (White 1982)
//
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训)
#pragma once

#include <cmath>       // std::log, std::exp, std::sqrt, std::fabs, std::abs
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

namespace cpphub {
inline namespace v1 {
namespace econometrics {

using linalg::dynamic::MatrixXD;
using linalg::dynamic::VectorXD;

// =============================================================================
// MLEFamily - MLE 分布族 (canonical link 标注)
//
// 注: Bernoulli 与 Logistic 等价 (Bernoulli y∈{0,1}, logit link);
//     保留 Bernoulli 名称仅为 API 兼容, 实现复用 Logistic.
// =============================================================================
enum class MLEFamily {
    Gaussian,           ///< 正态 (identity link, 闭式解 = OLS)
    Logistic,           ///< Logistic 回归 (logit link, canonical)
    Bernoulli,          ///< Bernoulli (同 Logistic, 别名)
    Probit,             ///< Probit 回归 (probit link, non-canonical)
    Poisson,            ///< Poisson 回归 (log link, canonical)
    NegativeBinomial    ///< 负二项 (log link, 含离散参数 α)
};

/// @brief MLEFamily 枚举转字符串
inline std::string to_string(MLEFamily f) {
    switch (f) {
        case MLEFamily::Gaussian:          return "Gaussian";
        case MLEFamily::Logistic:          return "Logistic";
        case MLEFamily::Bernoulli:         return "Bernoulli";
        case MLEFamily::Probit:            return "Probit";
        case MLEFamily::Poisson:           return "Poisson";
        case MLEFamily::NegativeBinomial:  return "NegativeBinomial";
    }
    return "Unknown";
}

// =============================================================================
// detail: MLE 内部辅助函数 (link/逆 link/权重/对数似然)
//
// 设计: 每个分布族提供 4 个函数:
//   - mean_fn(eta): 逆 link, μ = g^{-1}(Xβ) = E[y|x]
//   - link_fn(mu):  link, η = g(μ) = Xβ
//   - irls_weight(mu): IRLS 权重 w = (dμ/dη)² / V(μ) (for canonical link, w = V(μ))
//   - loglik_term(y, mu, alpha): 单观测对数似然 ℓ_i
//
// 数值稳定性:
//   - Logistic: 用 1/(1+exp(-η)) 计算, 大 |η| 时用 log1p 避免溢出
//   - Poisson: μ = exp(η), 大 η 时截断避免 inf
//   - Probit: 用标准正态 CDF (erfc 近似, Abramowitz-Stegun 7.1.26)
// =============================================================================
namespace detail {

// kTwoPi / kLogTwoPi 已移至 core/special_functions.hpp (共享定义点, 消除跨头文件重复)

// -----------------------------------------------------------------------------
// 标准正态 PDF φ(x) = (2π)^{-1/2} exp(-x²/2)
// -----------------------------------------------------------------------------
inline Real normal_pdf(Real x) {
    return std::exp(-0.5 * x * x) / std::sqrt(kTwoPi);
}

// -----------------------------------------------------------------------------
// 标准正态 CDF Φ(x) = 0.5 * erfc(-x/√2)
//   Abramowitz-Stegun 7.1.26: erfc 近似精度 ~1e-7 (足够用于 Probit MLE)
//   注: C++11 std::erfc 已达机器精度, 直接使用
// -----------------------------------------------------------------------------
inline Real normal_cdf(Real x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// =============================================================================
// Gaussian (identity link): μ = η, w = 1, V(μ) = σ²
// =============================================================================
struct GaussianFamily {
    static Real mean_fn(Real eta) { return eta; }
    static Real link_fn(Real mu) { return mu; }
    static Real irls_weight(Real /*mu*/, Real /*alpha*/) { return 1.0; }
    // ℓ_i = -0.5 * [log(2π) + log(σ²) + (y-μ)²/σ²]
    // 注: σ² 由样本估计 (MLE: σ² = SSR/N), 不在 β 迭代中
    static Real loglik_term(Real y, Real mu, Real sigma2) {
        if (sigma2 <= 0.0) return -std::numeric_limits<Real>::infinity();
        return -0.5 * (kLogTwoPi + std::log(sigma2) + (y - mu) * (y - mu) / sigma2);
    }
};

// =============================================================================
// Logistic / Bernoulli (logit link, canonical): μ = 1/(1+e^{-η}), w = μ(1-μ)
// =============================================================================
struct LogisticFamily {
    // μ = 1 / (1 + exp(-η)), 数值稳定版
    static Real mean_fn(Real eta) {
        if (eta >= 0.0) {
            const Real e = std::exp(-eta);
            return 1.0 / (1.0 + e);
        }
        const Real e = std::exp(eta);
        return e / (1.0 + e);
    }
    // η = log(μ / (1-μ))
    static Real link_fn(Real mu) {
        if (mu <= 0.0) return -std::numeric_limits<Real>::infinity();
        if (mu >= 1.0) return std::numeric_limits<Real>::infinity();
        return std::log(mu / (1.0 - mu));
    }
    // canonical link: w = μ(1-μ)
    static Real irls_weight(Real mu, Real /*alpha*/) {
        return mu * (1.0 - mu);
    }
    // ℓ_i = y·log(μ) + (1-y)·log(1-μ)
    // 数值稳定: 用 η 直接计算, 避免 log(0)
    static Real loglik_term(Real y, Real eta, Real /*alpha*/) {
        // ℓ = y·η - log(1 + e^η)  (canonical form)
        // 数值稳定 log(1+e^η): η<0 时 = log(1+e^η), η>=0 时 = η + log(1+e^{-η})
        Real log1pexp;
        if (eta >= 0.0) {
            log1pexp = eta + std::log1p(std::exp(-eta));
        } else {
            log1pexp = std::log1p(std::exp(eta));
        }
        return y * eta - log1pexp;
    }
};

// =============================================================================
// Probit (probit link, non-canonical): μ = Φ(η), w = φ(η)² / [Φ(η)(1-Φ(η))]
// =============================================================================
struct ProbitFamily {
    // μ = Φ(η)
    static Real mean_fn(Real eta) {
        return normal_cdf(eta);
    }
    // η = Φ^{-1}(μ)
    static Real link_fn(Real mu) {
        if (mu <= 0.0) return -std::numeric_limits<Real>::infinity();
        if (mu >= 1.0) return std::numeric_limits<Real>::infinity();
        // Newton-Raphson for inverse normal CDF (Acklam 近似 + 精炼)
        // 简化: 用二分法 (足够用于 link_fn, 通常不调用)
        Real lo = -10.0, hi = 10.0;
        for (int i = 0; i < 100; ++i) {
            Real mid = 0.5 * (lo + hi);
            if (normal_cdf(mid) < mu) lo = mid;
            else hi = mid;
        }
        return 0.5 * (lo + hi);
    }
    // non-canonical: w = φ(η)² / [Φ(η)(1-Φ(η))]
    static Real irls_weight(Real mu, Real /*alpha*/) {
        // eta = Φ^{-1}(μ), 但这里传 μ, 需反推 eta
        // 实际迭代中直接用 eta 计算 w, 此函数仅在需要时用 mu 近似
        if (mu <= 1e-12 || mu >= 1.0 - 1e-12) return 1e-12;
        const Real eta = link_fn(mu);
        const Real pdf = normal_pdf(eta);
        return pdf * pdf / (mu * (1.0 - mu));
    }
    // ℓ_i = y·log(Φ(η)) + (1-y)·log(1-Φ(η))
    static Real loglik_term(Real y, Real eta, Real /*alpha*/) {
        const Real p = normal_cdf(eta);
        const Real p_clipped = std::max(std::min(p, 1.0 - 1e-15), 1e-15);
        if (y > 0.5) return std::log(p_clipped);
        return std::log(1.0 - p_clipped);
    }
};

// =============================================================================
// Poisson (log link, canonical): μ = e^η, w = μ
// =============================================================================
struct PoissonFamily {
    // μ = exp(η), 截断避免溢出
    static Real mean_fn(Real eta) {
        if (eta > 50.0) return std::exp(50.0);  // 截断
        if (eta < -50.0) return std::exp(-50.0);
        return std::exp(eta);
    }
    // η = log(μ)
    static Real link_fn(Real mu) {
        if (mu <= 0.0) return -50.0;
        return std::log(mu);
    }
    // canonical: w = μ
    static Real irls_weight(Real mu, Real /*alpha*/) {
        return mu;
    }
    // ℓ_i = y·log(μ) - μ - log(y!)
    //   log(y!) 用 std::lgamma(y+1) (Stirling 精确版)
    static Real loglik_term(Real y, Real eta, Real /*alpha*/) {
        const Real mu = mean_fn(eta);
        return y * eta - mu - std::lgamma(y + 1.0);
    }
};

// =============================================================================
// Negative Binomial (log link, canonical for fixed α):
//   μ = e^η, V(μ) = μ + α·μ², w = μ²/(μ + α·μ²) = μ/(1+α·μ)
//   离散参数 α 固定时为 canonical; α 估计时为两步法
// =============================================================================
struct NegativeBinomialFamily {
    static Real mean_fn(Real eta) {
        if (eta > 50.0) return std::exp(50.0);
        if (eta < -50.0) return std::exp(-50.0);
        return std::exp(eta);
    }
    static Real link_fn(Real mu) {
        if (mu <= 0.0) return -50.0;
        return std::log(mu);
    }
    // w = μ / (1 + α·μ)
    static Real irls_weight(Real mu, Real alpha) {
        if (mu <= 0.0) return 1e-12;
        return mu / (1.0 + alpha * mu);
    }
    // ℓ_i = log Γ(y+1/α) - log Γ(y+1) - log Γ(1/α)
    //       + y·log(α·μ/(1+α·μ)) - (1/α)·log(1+α·μ)
    static Real loglik_term(Real y, Real eta, Real alpha) {
        const Real mu = mean_fn(eta);
        const Real inv_alpha = 1.0 / alpha;
        const Real log_alpha_mu = std::log(alpha * mu + 1e-300);
        const Real log_1_plus_alpha_mu = std::log1p(alpha * mu);
        return std::lgamma(y + inv_alpha) - std::lgamma(y + 1.0) - std::lgamma(inv_alpha)
               + y * (log_alpha_mu - log_1_plus_alpha_mu)
               - inv_alpha * log_1_plus_alpha_mu;
    }
};

}  // namespace detail

// =============================================================================
// MLEEstimator - 极大似然估计器
//
// 支持 6 种分布族 (Gaussian/Logistic/Bernoulli/Probit/Poisson/NegativeBinomial)
// 协方差: Hessian (默认) / OPG / Sandwich (QMLE)
// 优化: Newton-Raphson (IRLS for canonical link), 解析梯度/Hessian
//
// 用法:
//   MLEEstimator mle(MLEFamily::Logistic, CovarianceType::Hessian);
//   mle.setMaxIter(100).setTolerance(1e-10);
//   EstimationResult r = mle.estimate(data);
//
// 注: Gaussian family 退化为 OLS, 直接用闭式解 (不走 Newton-Raphson)
// =============================================================================
class MLEEstimator : public Estimator {
public:
    /// @brief 默认构造 (Gaussian + Hessian)
    MLEEstimator()
        : family_(MLEFamily::Gaussian)
        , max_iter_(100)
        , tol_(1e-10)
        , alpha_(0.0) {}

    /// @brief 指定分布族与协方差类型
    /// @param family MLE 分布族
    /// @param cov_type 协方差类型 (Hessian/OPG/Sandwich)
    explicit MLEEstimator(MLEFamily family, CovarianceType cov_type = CovarianceType::Hessian)
        : family_(family)
        , max_iter_(100)
        , tol_(1e-10)
        , alpha_(0.0) {
        cov_type_ = cov_type;
    }

    /// @brief 设置负二项离散参数 α (仅 NegativeBinomial 有效)
    MLEEstimator& setAlpha(Real alpha) { alpha_ = alpha; return *this; }

    /// @brief 设置最大迭代次数
    MLEEstimator& setMaxIter(Size n) { max_iter_ = n; return *this; }

    /// @brief 设置收敛容差 (β 相对变化)
    MLEEstimator& setTolerance(Real t) { tol_ = t; return *this; }

    /// @brief 设置初始值 (空 = 零初始化)
    MLEEstimator& setStartValues(const VectorXD& start) { start_ = start; return *this; }

    EstimationResult estimate(const EconData& data) override;
    std::string name() const override;
    EstimatorClass estimatorClass() const override { return EstimatorClass::Parametric; }
    std::unique_ptr<Estimator> clone() const override;

    /// @brief 获取最后一次迭代的收敛信息
    Size nIterations() const { return n_iter_; }
    bool converged() const { return converged_; }

    /// @brief 计算给定 β 下的对数似然值
    /// @param X 设计矩阵
    /// @param y 响应向量
    /// @param beta 系数
    /// @return 对数似然 ℓ(β) = Σ ℓ_i
    Real computeLogLikelihood(const MatrixXD& X, const VectorXD& y, const VectorXD& beta) const;

    /// @brief 计算 IRLS 权重对角向量 w_i
    /// @param X 设计矩阵
    /// @param beta 系数
    /// @return w 向量 (n)
    VectorXD computeIRLSWeights(const MatrixXD& X, const VectorXD& beta) const;

    /// @brief 计算拟合均值 μ_i = g^{-1}(Xβ)
    VectorXD computeFittedMeans(const MatrixXD& X, const VectorXD& beta) const;

    /// @brief 计算线性预测 η = Xβ
    VectorXD computeLinearPredictor(const MatrixXD& X, const VectorXD& beta) const;

private:
    MLEFamily family_;
    Size max_iter_;
    Real tol_;
    Real alpha_;          ///< NB 离散参数 (仅 NegativeBinomial)
    VectorXD start_;      ///< 初始值 (空 = 零)
    Size n_iter_ = 0;     ///< 实际迭代次数
    bool converged_ = false;

    /// @brief 从 EconData variant 提取 X, y (仅 CrossSectionData)
    void extractXY(const EconData& data, MatrixXD& X, VectorXD& y) const;

    /// @brief Newton-Raphson 迭代核心
    /// @param X, y, beta_init 初始值
    /// @return 收敛后的 β, 同时更新 n_iter_/converged_
    VectorXD newtonRaphson(const MatrixXD& X, const VectorXD& y, const VectorXD& beta_init);

    /// @brief 计算 Hessian 协方差 V = (X'WX)^{-1} (canonical link: Hessian 逆)
    MatrixXD computeHessianVcov(const MatrixXD& X, const VectorXD& W_diag) const;

    /// @brief 计算 OPG 协方差 V = (Σ g_i g_i')^{-1} (Berndt-Hall-Hall-Hausman)
    /// @param X 设计矩阵, residuals = y - μ
    MatrixXD computeOPGVcov(const MatrixXD& X, const VectorXD& residuals) const;

    /// @brief 计算 Sandwich 协acobian V = A^{-1} B A^{-1} (White 1982 QMLE)
    /// @param A_inv = (X'WX)^{-1}, B = X' diag(ε²) X
    MatrixXD computeSandwichVcov(const MatrixXD& A_inv, const MatrixXD& X,
                                  const VectorXD& residuals) const;
};

// =============================================================================
// 内联实现
// =============================================================================

inline std::string MLEEstimator::name() const {
    return "MLE(" + to_string(family_) + ")";
}

inline std::unique_ptr<Estimator> MLEEstimator::clone() const {
    auto cloned = std::make_unique<MLEEstimator>(family_, cov_type_);
    cloned->max_iter_ = max_iter_;
    cloned->tol_ = tol_;
    cloned->alpha_ = alpha_;
    cloned->start_ = start_;
    return cloned;
}

inline void MLEEstimator::extractXY(const EconData& data, MatrixXD& X, VectorXD& y) const {
    std::visit([&](auto&& d) {
        using T = std::decay_t<decltype(d)>;
        if constexpr (std::is_same_v<T, CrossSectionData>) {
            X = d.X;
            y = d.y;
        } else {
            throw std::invalid_argument(
                "MLEEstimator::extractXY: only CrossSectionData is supported");
        }
    }, data);
}

inline VectorXD MLEEstimator::computeLinearPredictor(const MatrixXD& X,
                                                       const VectorXD& beta) const {
    return VectorXD(X.eigen() * beta.eigen());
}

inline VectorXD MLEEstimator::computeFittedMeans(const MatrixXD& X,
                                                   const VectorXD& beta) const {
    const VectorXD eta = computeLinearPredictor(X, beta);
    const Size n = eta.size();
    VectorXD mu(n);
    for (Size i = 0; i < n; ++i) {
        switch (family_) {
            case MLEFamily::Gaussian:
                mu(i) = detail::GaussianFamily::mean_fn(eta(i));
                break;
            case MLEFamily::Logistic:
            case MLEFamily::Bernoulli:
                mu(i) = detail::LogisticFamily::mean_fn(eta(i));
                break;
            case MLEFamily::Probit:
                mu(i) = detail::ProbitFamily::mean_fn(eta(i));
                break;
            case MLEFamily::Poisson:
                mu(i) = detail::PoissonFamily::mean_fn(eta(i));
                break;
            case MLEFamily::NegativeBinomial:
                mu(i) = detail::NegativeBinomialFamily::mean_fn(eta(i));
                break;
        }
    }
    return mu;
}

inline VectorXD MLEEstimator::computeIRLSWeights(const MatrixXD& X,
                                                   const VectorXD& beta) const {
    const VectorXD mu = computeFittedMeans(X, beta);
    const Size n = mu.size();
    VectorXD w(n);
    for (Size i = 0; i < n; ++i) {
        switch (family_) {
            case MLEFamily::Gaussian:
                w(i) = detail::GaussianFamily::irls_weight(mu(i), alpha_);
                break;
            case MLEFamily::Logistic:
            case MLEFamily::Bernoulli:
                w(i) = detail::LogisticFamily::irls_weight(mu(i), alpha_);
                break;
            case MLEFamily::Probit: {
                // Probit: w = φ(η)² / [Φ(η)(1-Φ(η))], 需 η
                const Real eta = (X.eigen().row(i) * beta.eigen())(0);
                const Real pdf = detail::normal_pdf(eta);
                const Real p = mu(i);
                const Real denom = p * (1.0 - p);
                w(i) = (denom > 1e-15) ? (pdf * pdf / denom) : 1e-15;
                break;
            }
            case MLEFamily::Poisson:
                w(i) = detail::PoissonFamily::irls_weight(mu(i), alpha_);
                break;
            case MLEFamily::NegativeBinomial:
                w(i) = detail::NegativeBinomialFamily::irls_weight(mu(i), alpha_);
                break;
        }
    }
    return w;
}

inline Real MLEEstimator::computeLogLikelihood(const MatrixXD& X,
                                                  const VectorXD& y,
                                                  const VectorXD& beta) const {
    const VectorXD eta = computeLinearPredictor(X, beta);
    const Size n = y.size();
    Real ll = 0.0;

    if (family_ == MLEFamily::Gaussian) {
        // Gaussian: 需 σ², 用 MLE 估计 σ² = SSR/N
        const VectorXD mu = computeFittedMeans(X, beta);
        const Real ssr = (y.eigen() - mu.eigen()).squaredNorm();
        const Real sigma2 = ssr / static_cast<Real>(n);
        for (Size i = 0; i < n; ++i) {
            ll += detail::GaussianFamily::loglik_term(y(i), mu(i), sigma2);
        }
        return ll;
    }

    for (Size i = 0; i < n; ++i) {
        switch (family_) {
            case MLEFamily::Logistic:
            case MLEFamily::Bernoulli:
                ll += detail::LogisticFamily::loglik_term(y(i), eta(i), alpha_);
                break;
            case MLEFamily::Probit:
                ll += detail::ProbitFamily::loglik_term(y(i), eta(i), alpha_);
                break;
            case MLEFamily::Poisson:
                ll += detail::PoissonFamily::loglik_term(y(i), eta(i), alpha_);
                break;
            case MLEFamily::NegativeBinomial:
                ll += detail::NegativeBinomialFamily::loglik_term(y(i), eta(i), alpha_);
                break;
            case MLEFamily::Gaussian:
                break;  // 已在上方处理
        }
    }
    return ll;
}

inline VectorXD MLEEstimator::newtonRaphson(const MatrixXD& X,
                                              const VectorXD& y,
                                              const VectorXD& beta_init) {
    const Size n = X.rows();
    const Size k = X.cols();

    VectorXD beta = beta_init;
    converged_ = false;
    n_iter_ = 0;

    for (Size iter = 0; iter < max_iter_; ++iter) {
        n_iter_ = iter + 1;

        // 1. 计算 μ, W, score
        const VectorXD eta = computeLinearPredictor(X, beta);
        const VectorXD mu = computeFittedMeans(X, beta);
        const VectorXD W = computeIRLSWeights(X, beta);

        // score = X'(y - μ)
        const Eigen::VectorXd residual = y.eigen() - mu.eigen();
        const Eigen::VectorXd score = X.eigen().transpose() * residual;

        // 2. Hessian: H = -X'WX (canonical), NR 步长 = (X'WX)^{-1} X'(y-μ)
        //    Probit: H = -X'WX (non-canonical, W 同定义)
        Eigen::MatrixXd XtWX = Eigen::MatrixXd::Zero(k, k);
        for (Size i = 0; i < n; ++i) {
            XtWX += W(i) * (X.eigen().row(i).transpose() * X.eigen().row(i));
        }

        // 3. 解 (X'WX)·delta = score  via LLT (SPD)
        Eigen::LLT<Eigen::MatrixXd> llt(XtWX);
        if (llt.info() != Eigen::Success) {
            // X'WX 奇异: 完全分离 (Logistic) 或 μ=0 (Poisson)
            throw std::runtime_error(
                "MLEEstimator::newtonRaphson: X'WX is singular (possible complete separation)");
        }
        const Eigen::VectorXd delta = llt.solve(score);

        // 4. 更新 β
        const Eigen::VectorXd beta_new = beta.eigen() + delta;
        const Real rel_change = delta.norm() / (beta_new.norm() + 1e-300);
        beta = VectorXD(beta_new);

        // 5. 收敛检查: 相对变化 < tol
        if (rel_change < tol_) {
            converged_ = true;
            break;
        }
    }

    return beta;
}

inline MatrixXD MLEEstimator::computeHessianVcov(const MatrixXD& X,
                                                   const VectorXD& W_diag) const {
    // V = (X'WX)^{-1} via LLT
    const Size n = X.rows();
    const Size k = X.cols();
    Eigen::MatrixXd XtWX = Eigen::MatrixXd::Zero(k, k);
    for (Size i = 0; i < n; ++i) {
        XtWX += W_diag(i) * (X.eigen().row(i).transpose() * X.eigen().row(i));
    }
    Eigen::LLT<Eigen::MatrixXd> llt(XtWX);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("MLEEstimator: X'WX singular for Hessian vcov");
    }
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(k, k);
    return MatrixXD(llt.solve(I));
}

inline MatrixXD MLEEstimator::computeOPGVcov(const MatrixXD& X,
                                               const VectorXD& residuals) const {
    // OPG: V = (Σ g_i g_i')^{-1}, g_i = x_i · (y_i - μ_i)
    // 排幻觉点 E8: meat = X' diag(ε²) X, OPG 逆 = (X' diag(ε²) X)^{-1}
    const Size n = X.rows();
    const Size k = X.cols();
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(k, k);
    for (Size i = 0; i < n; ++i) {
        const Real eps2 = residuals(i) * residuals(i);
        G += eps2 * (X.eigen().row(i).transpose() * X.eigen().row(i));
    }
    Eigen::LLT<Eigen::MatrixXd> llt(G);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("MLEEstimator: OPG matrix singular");
    }
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(k, k);
    return MatrixXD(llt.solve(I));
}

inline MatrixXD MLEEstimator::computeSandwichVcov(const MatrixXD& A_inv,
                                                    const MatrixXD& X,
                                                    const VectorXD& residuals) const {
    // 排幻觉点 E8: Sandwich = A^{-1} · B · A^{-1}
    //   A = X'WX (Hessian), A^{-1} = (X'WX)^{-1} = bread
    //   B = X' diag(ε²) X = meat (ε = y - μ 原始残差)
    const Size n = X.rows();
    const Size k = X.cols();
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(k, k);
    for (Size i = 0; i < n; ++i) {
        const Real eps2 = residuals(i) * residuals(i);
        B += eps2 * (X.eigen().row(i).transpose() * X.eigen().row(i));
    }
    // V = A^{-1} · B · A^{-1}
    const Eigen::MatrixXd result = A_inv.eigen() * B * A_inv.eigen();
    return MatrixXD(result);
}

inline EstimationResult MLEEstimator::estimate(const EconData& data) {
    // 1. 提取 X, y
    MatrixXD X;
    VectorXD y;
    extractXY(data, X, y);

    const Size n = X.rows();
    const Size k = X.cols();
    if (k == 0) {
        throw std::invalid_argument("MLEEstimator::estimate: X has zero columns");
    }
    if (y.size() != n) {
        throw std::invalid_argument("MLEEstimator::estimate: X.rows() != y.size()");
    }
    if (n <= k) {
        throw std::invalid_argument("MLEEstimator::estimate: requires n > k");
    }

    // 协方差类型检查
    const CovarianceType ct = cov_type_;
    if (ct != CovarianceType::Hessian && ct != CovarianceType::OPG
        && ct != CovarianceType::Sandwich) {
        throw std::invalid_argument(
            "MLEEstimator::estimate: unsupported covariance type '" + to_string(ct) +
            "' (MLE supports Hessian/OPG/Sandwich)");
    }

    // 2. 初始值
    VectorXD beta0;
    if (start_.size() == k) {
        beta0 = start_;
    } else {
        beta0 = VectorXD(k);
        beta0.eigen().setZero();
    }

    // 3. 估计 β
    VectorXD beta;
    if (family_ == MLEFamily::Gaussian) {
        // Gaussian MLE = OLS 闭式解: β = (X'X)^{-1} X'y
        const MatrixXD XtX = MatrixXD(X.eigen().transpose() * X.eigen());
        const MatrixXD XtX_inv = inverse_symmetric(XtX);
        const VectorXD Xty = VectorXD(X.eigen().transpose() * y.eigen());
        beta = VectorXD(XtX_inv.eigen() * Xty.eigen());
        converged_ = true;
        n_iter_ = 1;
    } else {
        beta = newtonRaphson(X, y, beta0);
    }

    if (!converged_) {
        // 不抛异常, 但在结果中标记 (调用方可检查 nIterations/converged)
        // 注: 严格场景应抛异常, 但 MLE 数值优化可能因容差设置不同而"未收敛"
        //     结果仍可能可用, 由调用方判断
    }

    // 4. 计算 μ, residuals, W
    const VectorXD eta = computeLinearPredictor(X, beta);
    const VectorXD mu = computeFittedMeans(X, beta);
    const VectorXD residuals = VectorXD(y.eigen() - mu.eigen());
    const VectorXD W = computeIRLSWeights(X, beta);

    // 5. 协方差矩阵
    // 排幻觉点 G3 (Gaussian Hessian V = σ²·(X'X)^{-1}, 含 σ² 因子):
    //   GLM bread = (X'WX)^{-1}, 对 Gaussian 的 W=1 不含 σ², 需在 Hessian 协方差
    //   上乘 σ² (R sandwich::bread.glm 乘 dispersion 同款做法)
    //   OPG/Sandwich 用 raw residual ε=y-μ, meat 自带 σ² 信息, 无需修正
    //   非 Gaussian 族 (Logistic/Poisson/...) 的 irls_weight = V(μ) 已含正确尺度, 无需修正
    MatrixXD vcov;
    const MatrixXD A_inv = computeHessianVcov(X, W);  // (X'WX)^{-1} = bread

    if (ct == CovarianceType::Hessian) {
        if (family_ == MLEFamily::Gaussian) {
            // Gaussian: V = σ²_MLE · (X'X)^{-1}
            const Real sigma2_mle = residuals.eigen().squaredNorm() / static_cast<Real>(n);
            vcov = MatrixXD(A_inv.eigen() * sigma2_mle);
        } else {
            vcov = A_inv;
        }
    } else if (ct == CovarianceType::OPG) {
        if (family_ == MLEFamily::Gaussian) {
            // Gaussian: V = σ⁴ · (X' diag(ε²) X)^{-1}
            //   原因: g_i = x_i·ε_i/σ², OPG = (Σ g_i g_i')^{-1} = σ⁴ · (X' diag(ε²) X)^{-1}
            //   mle.hpp 默认 OPG 用 g_i = x_i·ε_i (无 1/σ²), 需乘 σ⁴ 校正
            const Real sigma2_mle = residuals.eigen().squaredNorm() / static_cast<Real>(n);
            const Real sigma4_mle = sigma2_mle * sigma2_mle;
            const MatrixXD opg_raw = computeOPGVcov(X, residuals);
            vcov = MatrixXD(opg_raw.eigen() * sigma4_mle);
        } else {
            vcov = computeOPGVcov(X, residuals);
        }
    } else {  // Sandwich
        // Sandwich V = (X'WX)^{-1} · X' diag(ε²) X · (X'WX)^{-1}
        // 对 Gaussian: W=1, 退化为 HC0 White 1980, 无需 σ² 修正
        // 对非 Gaussian: W=V(μ), bread 含正确尺度, meat 含 ε² 自带尺度
        vcov = computeSandwichVcov(A_inv, X, residuals);
    }

    // 6. 标准误 / z 统计量 / p 值 (大样本用正态近似, 非 t)
    //    排幻觉点: MLE 用 z (标准正态), 不用 t (t 是 OLS 小样本)
    VectorXD se(k);
    VectorXD zstat(k);
    VectorXD pval(k);
    for (Size i = 0; i < k; ++i) {
        Real v = vcov(i, i);
        if (v < 0.0) v = 0.0;  // 数值保护
        const Real s = std::sqrt(v);
        se(i) = s;
        if (s > 0.0) {
            const Real z = beta(i) / s;
            zstat(i) = z;
            // 双侧 p 值: p = 2·(1 - Φ(|z|)) = erfc(|z|/√2)
            pval(i) = std::erfc(std::fabs(z) / std::sqrt(2.0));
        } else {
            zstat(i) = std::numeric_limits<Real>::infinity();
            pval(i) = 0.0;
        }
    }

    // 7. 对数似然
    const Real loglik = computeLogLikelihood(X, y, beta);

    // 8. R² (仅 Gaussian 有意义; 其他用 McFadden pseudo-R²)
    Real r2 = 0.0;
    Real adj_r2 = 0.0;
    if (family_ == MLEFamily::Gaussian) {
        const Real ssr = residuals.eigen().squaredNorm();
        const Real ybar = y.eigen().mean();
        const Real sst = (y.eigen().array() - ybar).matrix().squaredNorm();
        r2 = (sst > 0.0) ? (1.0 - ssr / sst) : 0.0;
        adj_r2 = (n > k) ? (1.0 - (1.0 - r2) * static_cast<Real>(n - 1) / static_cast<Real>(n - k)) : 0.0;
    } else {
        // McFadden pseudo-R² = 1 - ℓ/ℓ_0, ℓ_0 = 仅截距模型的对数似然
        // 简化: 用 0 (pseudo-R² 不在所有分布族有定义, 由调用方按需计算)
        r2 = 0.0;
        adj_r2 = 0.0;
    }

    // 9. 填充结果
    EstimationResult result;
    result.coefficients = beta;
    result.std_errors = se;
    result.t_statistics = zstat;  // MLE 中实为 z, 复用字段
    result.p_values = pval;
    result.vcov = vcov;
    result.log_likelihood = loglik;
    result.r_squared = r2;
    result.adj_r_squared = adj_r2;
    result.n_obs = n;
    result.n_params = k;
    result.df_residual = n - k;
    result.cov_type = ct;
    return result;
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
