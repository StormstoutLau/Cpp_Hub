// =============================================================================
// kde.hpp - 1D Gaussian KDE + 带宽敏感性网格 (Phase 9 M1, PHASE9_SPEC v1.1 §3.2)
//
// pdf(x)     = 1/(n·h·√(2π)) · Σ_i exp(−½((x−x_i)/h)²)
// log_pdf(x) = log_norm + logsumexp(−½ z_i²)  (max-shift, 深尾不产 NaN)
//
// 与 scipy gaussian_kde(bw_method=factor) 同核同 h 时 pdf 一致 (probe P6:
// 7.494e-16); 带宽规则语义见 bandwidth.hpp.
//
// 范围 (D9-1/D9-8): 仅 1D Gaussian 核、正支撑数据; reflection 边界修正
// 为 P2 预留, v1.9 显式抛 not-implemented.
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
    /// P2 预留; v1.9 必须 false, true 抛 std::logic_error
    bool reflection = false;
};

/// 一维 Gaussian 核密度估计 (header-only, O(N) 单点求值)
class KernelDensity1D {
public:
    /// @throws std::logic_error reflection=true (P2 未实现)
    /// @throws std::invalid_argument n < 2 / 非有限样本 / h 非正
    explicit KernelDensity1D(std::span<const Real> sample, KdeConfig config = {}) {
        if (config.reflection) {
            throw std::logic_error(
                "KDE: reflection boundary correction is not implemented (P2)");
        }
        h_ = resolve_bandwidth(sample, config.bandwidth);
        rule_ = config.bandwidth.rule;
        sample_.assign(sample.begin(), sample.end());
        log_norm_ = -std::log(static_cast<Real>(sample_.size())) - std::log(h_) -
                    std::log(SQRT_2PI);
        inv_h_ = 1.0 / h_;  // pdf 内层用乘法 (perf 预算 D9-10)
    }

    /// 密度: 1/(n h √(2π)) Σ exp(−½ z²)
    /// 双累加器解 FP 加法依赖链 (D9-10 perf; UCRT exp2 实测比 exp 慢 8×,
    /// 禁用 exp2 变换 — 2026-11-11 probe 记录)
    Real pdf(Real x) const {
        const Real* p = sample_.data();
        const Size n = sample_.size();
        Real s0 = 0.0;
        Real s1 = 0.0;
        Size i = 0;
        for (; i + 2 <= n; i += 2) {
            const Real z0 = (x - p[i]) * inv_h_;
            const Real z1 = (x - p[i + 1]) * inv_h_;
            s0 += std::exp(-0.5 * z0 * z0);
            s1 += std::exp(-0.5 * z1 * z1);
        }
        if (i < n) {
            const Real z = (x - p[i]) * inv_h_;
            s0 += std::exp(-0.5 * z * z);
        }
        return (s0 + s1) * std::exp(log_norm_);
    }

    /// 对数密度: max-shift log-sum-exp, 深尾返回 -inf 而非 NaN
    Real log_pdf(Real x) const {
        Real m = -std::numeric_limits<Real>::infinity();
        for (Real xi : sample_) {
            const Real z = (x - xi) * inv_h_;
            m = std::max(m, -0.5 * z * z);
        }
        if (!std::isfinite(m)) {
            return -std::numeric_limits<Real>::infinity();
        }
        Real s = 0.0;
        for (Real xi : sample_) {
            const Real z = (x - xi) * inv_h_;
            s += std::exp(-0.5 * z * z - m);
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

    /// 解析后的绝对带宽 h (数据单位)
    Real bandwidth() const noexcept { return h_; }
    BandwidthRule rule() const noexcept { return rule_; }
    std::size_t size() const noexcept { return sample_.size(); }

private:
    std::vector<Real> sample_;
    Real h_ = 0.0;
    Real inv_h_ = 0.0;     ///< 1/h (内层乘法)
    Real log_norm_ = 0.0;  ///< −log(n h √(2π))
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
