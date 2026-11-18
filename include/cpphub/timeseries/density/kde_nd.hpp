// =============================================================================
// kde_nd.hpp - 多维 Gaussian KDE (全协方差带宽, P2-3)
//
// 语义: 逐步镜像 scipy `gaussian_kde` (probe_multivariate_kde.py):
//   factor = n^(-1/(d+4))                  (kScott, scipy 默认)
//          = (n(d+2)/4)^(-1/(d+4))          (kSilverman)
//          = 用户给定                       (kUserFactor)
//   H      = factor² · S                    (S = 样本协方差, ddof=1)
//   pdf(x) = (2π)^(-d/2) |H|^(-1/2) / n Σ_i exp(−½ (x−x_i)ᵀ H⁻¹ (x−x_i))
//   log_pdf = max-shift log-sum-exp (深尾返回 −inf 而非 NaN)
//
// 注意: 与 1D `KernelDensity1D` 的带宽语义不同 — 1D 用 statsmodels/R 的绝对 h
// (1.059·A·n^(-1/5)), 本类用 scipy 的无量纲因子 × 协方差矩阵; 两者在 d=1 时
// 数值不同 (scott: n^(-1/5)·sd vs 1.059·A·n^(-1/5)), 不互相替代.
//
// 带宽矩阵 H 通过 Cholesky 分解求值 (header-only, 不依赖 Eigen3).
// =============================================================================

#pragma once

#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace density {

/// 多维 KDE 带宽规则 (scipy gaussian_kde 语义)
enum class MultivariateBandwidthRule {
    kScott,       ///< factor = n^(-1/(d+4))
    kSilverman,   ///< factor = (n(d+2)/4)^(-1/(d+4))
    kUserFactor   ///< factor = config.factor (> 0)
};

struct MultivariateKdeConfig {
    MultivariateBandwidthRule rule = MultivariateBandwidthRule::kScott;
    /// kUserFactor 时使用; H = factor²·S (无量纲因子)
    Real factor = 0.0;
};

/// 多维 Gaussian 核密度估计 (全协方差带宽, O(N d²) 单点求值)
class KernelDensityND {
public:
    /// @param points 行主序 n×dim 样本 (每行一个观测)
    /// @throws std::invalid_argument n < 2 / dim < 1 / 非有限值 /
    ///         factor 非正 / 协方差矩阵非正定
    KernelDensityND(std::span<const Real> points, Size n, Size dim,
                    MultivariateKdeConfig config = {})
        : n_(n), dim_(dim) {
        if (dim_ < 1) {
            throw std::invalid_argument("KDE: dimension must be >= 1");
        }
        if (n_ < 2) {
            throw std::invalid_argument("KDE: sample size must be >= 2");
        }
        if (points.size() != n_ * dim_) {
            throw std::invalid_argument("KDE: points size must be n * dim");
        }
        for (Real v : points) {
            if (!std::isfinite(v)) {
                throw std::invalid_argument("KDE: sample contains non-finite value");
            }
        }
        points_.assign(points.begin(), points.end());

        factor_ = resolve_factor(config);
        compute_covariance();
        h_.assign(dim_ * dim_, 0.0);
        const Real f2 = factor_ * factor_;
        for (Size i = 0; i < dim_ * dim_; ++i) {
            h_[i] = f2 * cov_[i];
        }
        cholesky_factor();
    }

    /// 密度 pdf(x) (x 长度 = dim)
    Real pdf(std::span<const Real> x) const {
        require_point(x);
        const Real* p = points_.data();
        std::vector<Real> z(dim_);
        Real s0 = 0.0;
        Real s1 = 0.0;
        Size i = 0;
        for (; i + 2 <= n_; i += 2) {
            s0 += std::exp(-0.5 * quadratic(x.data(), p + i * dim_, z.data()));
            s1 += std::exp(-0.5 * quadratic(x.data(), p + (i + 1) * dim_, z.data()));
        }
        if (i < n_) {
            s0 += std::exp(-0.5 * quadratic(x.data(), p + i * dim_, z.data()));
        }
        return (s0 + s1) * std::exp(log_norm_);
    }

    /// 对数密度 (max-shift log-sum-exp)
    Real log_pdf(std::span<const Real> x) const {
        require_point(x);
        std::vector<Real> z(dim_);
        Real m = -std::numeric_limits<Real>::infinity();
        for (Size i = 0; i < n_; ++i) {
            m = std::max(m, -0.5 * quadratic(x.data(), &points_[i * dim_], z.data()));
        }
        if (!std::isfinite(m)) {
            return -std::numeric_limits<Real>::infinity();
        }
        Real s = 0.0;
        for (Size i = 0; i < n_; ++i) {
            s += std::exp(-0.5 * quadratic(x.data(), &points_[i * dim_], z.data()) - m);
        }
        return log_norm_ + m + std::log(s);
    }

    /// 带宽矩阵元素 H(i,j) = factor²·S(i,j)
    Real bandwidth(Size i, Size j) const {
        if (i >= dim_ || j >= dim_) {
            throw std::out_of_range("KDE: bandwidth index out of range");
        }
        return h_[i * dim_ + j];
    }
    Real factor() const noexcept { return factor_; }
    Size dimension() const noexcept { return dim_; }
    Size size() const noexcept { return n_; }

private:
    Real resolve_factor(const MultivariateKdeConfig& config) const {
        const Real n_real = static_cast<Real>(n_);
        const Real d4 = static_cast<Real>(dim_ + 4);
        Real f = 0.0;
        switch (config.rule) {
            case MultivariateBandwidthRule::kScott:
                f = std::pow(n_real, -1.0 / d4);
                break;
            case MultivariateBandwidthRule::kSilverman:
                f = std::pow(n_real * (static_cast<Real>(dim_) + 2.0) / 4.0,
                             -1.0 / d4);
                break;
            case MultivariateBandwidthRule::kUserFactor:
                f = config.factor;
                break;
        }
        if (!(f > 0.0) || !std::isfinite(f)) {
            throw std::invalid_argument("KDE: bandwidth factor must be positive and finite");
        }
        return f;
    }

    void compute_covariance() {
        std::vector<Real> mean(dim_, 0.0);
        for (Size i = 0; i < n_; ++i) {
            for (Size j = 0; j < dim_; ++j) {
                mean[j] += points_[i * dim_ + j];
            }
        }
        for (Size j = 0; j < dim_; ++j) {
            mean[j] /= static_cast<Real>(n_);
        }
        cov_.assign(dim_ * dim_, 0.0);
        for (Size i = 0; i < n_; ++i) {
            for (Size a = 0; a < dim_; ++a) {
                const Real da = points_[i * dim_ + a] - mean[a];
                for (Size b = 0; b <= a; ++b) {
                    cov_[a * dim_ + b] += da * (points_[i * dim_ + b] - mean[b]);
                }
            }
        }
        const Real inv = 1.0 / static_cast<Real>(n_ - 1);
        for (Size a = 0; a < dim_; ++a) {
            for (Size b = 0; b <= a; ++b) {
                const Real v = cov_[a * dim_ + b] * inv;
                cov_[a * dim_ + b] = v;
                cov_[b * dim_ + a] = v;
            }
        }
    }

    /// H = L Lᵀ (下三角, 行主序)
    void cholesky_factor() {
        chol_.assign(dim_ * dim_, 0.0);
        Real log_det = 0.0;
        for (Size j = 0; j < dim_; ++j) {
            Real s = h_[j * dim_ + j];
            for (Size k = 0; k < j; ++k) {
                const Real ljk = chol_[j * dim_ + k];
                s -= ljk * ljk;
            }
            if (!(s > 0.0) || !std::isfinite(s)) {
                throw std::invalid_argument(
                    "KDE: bandwidth matrix is not positive definite");
            }
            const Real ljj = std::sqrt(s);
            chol_[j * dim_ + j] = ljj;
            log_det += std::log(ljj);
            for (Size i = j + 1; i < dim_; ++i) {
                Real t = h_[i * dim_ + j];
                for (Size k = 0; k < j; ++k) {
                    t -= chol_[i * dim_ + k] * chol_[j * dim_ + k];
                }
                chol_[i * dim_ + j] = t / ljj;
            }
        }
        log_norm_ = -0.5 * static_cast<Real>(dim_) * std::log(2.0 * kPi) -
                    log_det - std::log(static_cast<Real>(n_));
    }

    /// (x−xi)ᵀ H⁻¹ (x−xi) = ‖L⁻¹(x−xi)‖², 前代 O(d²); z 为每次求值复用的暂存
    Real quadratic(const Real* x, const Real* xi, Real* z) const {
        for (Size j = 0; j < dim_; ++j) {
            Real t = x[j] - xi[j];
            for (Size k = 0; k < j; ++k) {
                t -= chol_[j * dim_ + k] * z[k];
            }
            z[j] = t / chol_[j * dim_ + j];
        }
        Real q = 0.0;
        for (Size j = 0; j < dim_; ++j) {
            q += z[j] * z[j];
        }
        return q;
    }

    void require_point(std::span<const Real> x) const {
        if (x.size() != dim_) {
            throw std::invalid_argument("KDE: evaluation point dimension mismatch");
        }
        for (Real v : x) {
            if (!std::isfinite(v)) {
                throw std::invalid_argument("KDE: evaluation point contains non-finite value");
            }
        }
    }

    static constexpr Real kPi = 3.141592653589793238462643383279502884;

    std::vector<Real> points_;
    std::vector<Real> cov_;
    std::vector<Real> h_;
    std::vector<Real> chol_;
    Size n_ = 0;
    Size dim_ = 0;
    Real factor_ = 0.0;
    Real log_norm_ = 0.0;
};

}  // namespace density
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
