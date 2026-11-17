// =============================================================================
// kde.hpp - 1D Gaussian KDE + 带宽敏感性网格 (Phase 9 M1, PHASE9_SPEC v1.1 §3.2)
//
// pdf(x)     = 1/(n·h·√(2π)) · Σ_i exp(−½((x−x_i)/h)²)
// log_pdf(x) = log_norm + logsumexp(−½ z_i²)  (max-shift, 深尾不产 NaN)
//
// 与 scipy gaussian_kde(bw_method=factor) 同核同 h 时 pdf 一致 (probe P6:
// 7.494e-16); 带宽规则语义见 bandwidth.hpp.
//
// 范围 (D9-1/D9-8): 仅 1D Gaussian 核. P2-1 已落地 0 处反射边界修正
// (reflection=true, 支撑 [0,∞)); P2-5 已落地权重 (混合权重 w_i/W,
// 带宽规则仍按未加权样本解析 — statsmodels/sklearn 同语义).
// =============================================================================

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include "cpphub/core/constants.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/density/bandwidth.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace density {

/// KDE 配置
struct KdeConfig {
    Bandwidth bandwidth{};
    /// true: 0 处反射边界修正 (支撑 [0,∞), 样本须严格为正)
    bool reflection = false;
};

/// 一维 Gaussian 核密度估计 (header-only, O(N) 单点求值)
///
/// 无权:   pdf(x) = 1/(n h √(2π)) Σ exp(−½ z_i²)
/// 加权:   pdf(x) = 1/(h √(2π)) Σ (w_i/W) exp(−½ z_i²)     (P2-5)
///   - 权重仅改混合权重; 带宽规则按**未加权**样本解析 (statsmodels
///     `kdensity` / sklearn `sample_weight` 同语义; scipy gaussian_kde
///     用加权协方差 + neff, 本库不采用 — 见 SPEC §3.2)
///   - 权重仅比例有意义; 允许 0 权重; neff=(Σw)²/Σw² 仅作诊断
/// 反射:   pdf_R(x) = 1/(h √(2π)) Σ p_i [φ(z_i) + φ((x+x_i)/h)], x ≥ 0  (P2-1)
///   - ∫_0^∞ pdf_R = 1 (每个观测的两半核质量相加为 1); x < 0 返回 0 / −inf
class KernelDensity1D {
public:
    /// @throws std::invalid_argument n < 2 / 非有限样本 / h 非正 /
    ///         reflection 且样本非严格正
    explicit KernelDensity1D(std::span<const Real> sample, KdeConfig config = {})
        : KernelDensity1D(sample, std::span<const Real>(), config) {}

    /// @param weights 非负有限权重 (可未归一; 仅比例有意义); 空 span = 等权
    /// @throws std::invalid_argument 权重长度不符 / 含负值或非有限值 / 和为零
    KernelDensity1D(std::span<const Real> sample, std::span<const Real> weights,
                    KdeConfig config = {})
        : reflection_(config.reflection), rule_(config.bandwidth.rule) {
        sample_.assign(sample.begin(), sample.end());
        n_ = sample_.size();
        if (reflection_) {
            for (Real xi : sample_) {
                if (!(xi > 0.0)) {
                    throw std::invalid_argument(
                        "KDE: reflection requires strictly positive sample");
                }
            }
        }
        h_ = resolve_bandwidth(sample, config.bandwidth);
        inv_h_ = 1.0 / h_;
        if (weights.empty()) {
            weighted_ = false;
            neff_ = static_cast<Real>(n_);
            log_norm_ = -std::log(static_cast<Real>(n_)) - std::log(h_) -
                        std::log(SQRT_2PI);
        } else {
            init_weights(weights);
        }
    }

    /// 密度 (双累加器解 FP 加法依赖链; UCRT exp2 实测比 exp 慢 8×, 禁用)
    Real pdf(Real x) const {
        if (reflection_ && x < 0.0) {
            return 0.0;
        }
        const Real* p = sample_.data();
        Real s0 = 0.0;
        Real s1 = 0.0;
        Size i = 0;
        if (weighted_) {
            const Real* w = weights_.data();
            for (; i + 2 <= n_; i += 2) {
                s0 += w[i] * kernel_at(x, p[i]);
                s1 += w[i + 1] * kernel_at(x, p[i + 1]);
            }
            if (i < n_) {
                s0 += w[i] * kernel_at(x, p[i]);
            }
        } else {
            for (; i + 2 <= n_; i += 2) {
                s0 += kernel_at(x, p[i]);
                s1 += kernel_at(x, p[i + 1]);
            }
            if (i < n_) {
                s0 += kernel_at(x, p[i]);
            }
        }
        return (s0 + s1) * std::exp(log_norm_);
    }

    /// 对数密度: max-shift log-sum-exp, 深尾返回 -inf 而非 NaN
    Real log_pdf(Real x) const {
        if (reflection_ && x < 0.0) {
            return -std::numeric_limits<Real>::infinity();
        }
        Real m = -std::numeric_limits<Real>::infinity();
        if (weighted_) {
            for (Size i = 0; i < n_; ++i) {
                if (!std::isfinite(log_weights_[i])) {
                    continue;  // 0 权重
                }
                m = std::max(m, log_weights_[i] + log_kernel_at(x, sample_[i]));
            }
        } else {
            for (Size i = 0; i < n_; ++i) {
                m = std::max(m, log_kernel_at(x, sample_[i]));
            }
        }
        if (!std::isfinite(m)) {
            return -std::numeric_limits<Real>::infinity();
        }
        Real s = 0.0;
        if (weighted_) {
            for (Size i = 0; i < n_; ++i) {
                if (!std::isfinite(log_weights_[i])) {
                    continue;
                }
                s += std::exp(log_weights_[i] + log_kernel_at(x, sample_[i]) - m);
            }
        } else {
            for (Size i = 0; i < n_; ++i) {
                s += std::exp(log_kernel_at(x, sample_[i]) - m);
            }
        }
        return log_norm_ + m + std::log(s);
    }

    /// 网格批量求值 (等价于逐点 pdf)
    std::vector<Real> pdf(std::span<const Real> grid) const {
        std::vector<Real> out;
        out.reserve(grid.size());
        for (Real x : grid) {
            out.push_back(pdf(x));
        }
        return out;
    }

    /// 解析后的绝对带宽 h (数据单位; 权重不影响 h)
    Real bandwidth() const noexcept { return h_; }
    BandwidthRule rule() const noexcept { return rule_; }
    std::size_t size() const noexcept { return n_; }
    /// 是否使用权重 (P2-5)
    bool weighted() const noexcept { return weighted_; }
    /// 有效样本量 (Σw)²/Σw²; 无权时 = n (仅诊断, 不参与带宽)
    Real effective_sample_size() const noexcept { return neff_; }

private:
    /// 单点未加权核值; reflection 时加另一半核 (x+x_i)/h
    Real kernel_at(Real x, Real xi) const {
        const Real z = (x - xi) * inv_h_;
        Real k = std::exp(-0.5 * z * z);
        if (reflection_) {
            const Real zr = (x + xi) * inv_h_;
            k += std::exp(-0.5 * zr * zr);
        }
        return k;
    }

    /// log(核值之和): reflection 时两项取 logaddexp 的稳定形式
    Real log_kernel_at(Real x, Real xi) const {
        const Real z = (x - xi) * inv_h_;
        const Real a = -0.5 * z * z;
        if (!reflection_) {
            return a;
        }
        const Real zr = (x + xi) * inv_h_;
        const Real b = -0.5 * zr * zr;
        const Real hi = std::max(a, b);
        const Real lo = std::min(a, b);
        return hi + std::log1p(std::exp(lo - hi));
    }

    void init_weights(std::span<const Real> weights) {
        if (weights.size() != n_) {
            throw std::invalid_argument("KDE: weights size must match sample size");
        }
        Real sum = 0.0;
        Real sum_sq = 0.0;
        for (Real w : weights) {
            if (!std::isfinite(w) || w < 0.0) {
                throw std::invalid_argument(
                    "KDE: weights must be finite and non-negative");
            }
            sum += w;
            sum_sq += w * w;
        }
        if (!(sum > 0.0) || !std::isfinite(sum) || !(sum_sq > 0.0)) {
            throw std::invalid_argument("KDE: weights must have positive finite sum");
        }
        weighted_ = true;
        neff_ = sum * sum / sum_sq;
        weights_.resize(n_);
        log_weights_.resize(n_);
        for (Size i = 0; i < n_; ++i) {
            const Real w = weights[i] / sum;
            weights_[i] = w;
            log_weights_[i] =
                (w > 0.0) ? std::log(w) : -std::numeric_limits<Real>::infinity();
        }
        log_norm_ = -std::log(h_) - std::log(SQRT_2PI);
    }

    std::vector<Real> sample_;
    std::vector<Real> weights_;      ///< 归一化 w_i/W (加权时)
    std::vector<Real> log_weights_;  ///< log(w_i/W) (0 权重为 −inf)
    Size n_ = 0;
    Real h_ = 0.0;
    Real inv_h_ = 0.0;     ///< 1/h (内层乘法)
    Real log_norm_ = 0.0;  ///< 无权: −log(n h √(2π)); 加权: −log(h √(2π))
    Real neff_ = 0.0;
    bool weighted_ = false;
    bool reflection_ = false;
    BandwidthRule rule_ = BandwidthRule::kScott;
};

/// 多带宽敏感性诊断点
struct BandwidthSensitivityPoint {
    Real bandwidth = 0.0;      ///< 绝对 h
    Real log_likelihood = 0.0; ///< Σ log pdf(x_i) (不含留一修正, P2)
    Real pdf_at_mode = 0.0;    ///< 评估网格上的 pdf 最大值
};

/// @brief 对每个绝对 h 重建 KDE, 返回 log-likelihood 与网格峰值
/// @throws std::invalid_argument 样本非法 / h 非正 / eval_grid 为空
inline std::vector<BandwidthSensitivityPoint> bandwidth_sensitivity(
    std::span<const Real> sample, std::span<const Real> bandwidths,
    std::span<const Real> eval_grid) {
    if (eval_grid.empty()) {
        throw std::invalid_argument("KDE: eval_grid must be non-empty");
    }
    std::vector<BandwidthSensitivityPoint> out;
    out.reserve(bandwidths.size());
    for (Real h : bandwidths) {
        KdeConfig cfg;
        cfg.bandwidth.rule = BandwidthRule::kUser;
        cfg.bandwidth.value = h;
        cfg.bandwidth.is_factor = false;
        const KernelDensity1D kde(sample, cfg);

        BandwidthSensitivityPoint p;
        p.bandwidth = kde.bandwidth();
        p.log_likelihood = 0.0;
        for (Real xi : sample) {
            p.log_likelihood += kde.log_pdf(xi);
        }
        p.pdf_at_mode = kde.pdf(eval_grid[0]);
        for (Size i = 1; i < eval_grid.size(); ++i) {
            p.pdf_at_mode = std::max(p.pdf_at_mode, kde.pdf(eval_grid[i]));
        }
        out.push_back(p);
    }
    return out;
}

}  // namespace density
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
