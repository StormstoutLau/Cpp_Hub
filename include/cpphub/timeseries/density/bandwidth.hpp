// =============================================================================
// bandwidth.hpp - KDE 带宽类型与规则解析 (Phase 9 M1, PHASE9_SPEC v1.1 §3.1)
//
// 同名不同物的带宽语义在此显式区分 (调研 A1-A25):
//   kScott           1.059 · A · n^(-1/5)   (statsmodels bw_scott / R bw.nrd)
//   kSilverman       0.9   · A · n^(-1/5)   (statsmodels bw_silverman / R bw.nrd0)
//   kNormalReference C     · A · n^(-1/5)   (statsmodels 默认; C 由 probe P4 反解)
//   kUser + is_factor=false  绝对 h          (statsmodels / R 语义)
//   kUser + is_factor=true   h = value · sd(ddof=1) (scipy gaussian_kde bw_method)
// 其中 A = min(sd(ddof=1), IQR/1.349), IQR = q75 − q25 (线性插值分位数).
//
// 不提供 sklearn KernelDensity(bandwidth='scott') 的 "n^(-1/(d+4)) 裸因子当 h"
// 语义 (调研 A1-A3): 该因子无量纲, 当成绝对 h 会让密度随数据单位变化.
// =============================================================================

#pragma once

#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace density {

/// 带宽规则 (枚举数值仅差 ~0.02%, 保留独立枚举以便未来换核)
enum class BandwidthRule { kUser, kScott, kSilverman, kNormalReference };

/// 带宽规格
struct Bandwidth {
    BandwidthRule rule = BandwidthRule::kScott;
    /// kUser: 绝对 h (is_factor=false) 或无量纲因子 (is_factor=true)
    Real value = 0.0;
    /// 仅 kUser 有效: true 时 h = value · sd(ddof=1) (scipy 语义)
    bool is_factor = false;
};

namespace detail {

/// statsmodels normal_reference 常数 (probe P4 反解:
/// 0.337107228 / (0.918303663 · 200^(-1/5))) — 非 1.059 (审计 A2)
inline constexpr Real kNormalReferenceConstant = 1.059223840025653;

/// 线性插值分位数 (numpy.percentile 默认 "linear"): pos = q·(n−1)
inline Real linear_quantile(std::span<const Real> sorted, Real q) {
    const Size n = sorted.size();
    const Real pos = q * static_cast<Real>(n - 1);
    const Size lo = static_cast<Size>(pos);
    if (lo >= n - 1) {
        return sorted[n - 1];
    }
    const Real frac = pos - static_cast<Real>(lo);
    return sorted[lo] + frac * (sorted[lo + 1] - sorted[lo]);
}

/// 样本标准差 (ddof=1)
inline Real sample_sd(std::span<const Real> sample) {
    const Real n = static_cast<Real>(sample.size());
    Real mean = 0.0;
    for (Real v : sample) {
        mean += v;
    }
    mean /= n;
    Real ss = 0.0;
    for (Real v : sample) {
        const Real d = v - mean;
        ss += d * d;
    }
    return std::sqrt(ss / (n - 1.0));
}

}  // namespace detail

/// @brief 规则 → 绝对带宽 h (数据单位)
/// @throws std::invalid_argument 样本 n < 2 / 含非有限值 / 解析出的 h 非正或非有限
inline Real resolve_bandwidth(std::span<const Real> sample, const Bandwidth& bw) {
    if (sample.size() < 2) {
        throw std::invalid_argument("KDE: sample size must be >= 2");
    }
    for (Real v : sample) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("KDE: sample contains non-finite value");
        }
    }

    if (bw.rule == BandwidthRule::kUser) {
        Real h = bw.value;
        if (bw.is_factor) {
            h *= detail::sample_sd(sample);
        }
        if (!std::isfinite(h) || h <= 0.0) {
            throw std::invalid_argument("KDE: bandwidth must be positive and finite");
        }
        return h;
    }

    const Real sd = detail::sample_sd(sample);
    std::vector<Real> sorted(sample.begin(), sample.end());
    std::sort(sorted.begin(), sorted.end());
    const Real iqr = detail::linear_quantile(sorted, 0.75) -
                     detail::linear_quantile(sorted, 0.25);
    const Real a = std::min(sd, iqr / 1.349);
    const Real n_pow = std::pow(static_cast<Real>(sample.size()), -0.2);

    Real h = 0.0;
    switch (bw.rule) {
        case BandwidthRule::kScott:
            h = 1.059 * a * n_pow;
            break;
        case BandwidthRule::kSilverman:
            h = 0.9 * a * n_pow;
            break;
        case BandwidthRule::kNormalReference:
            h = detail::kNormalReferenceConstant * a * n_pow;
            break;
        case BandwidthRule::kUser:
            break;  // 上面已返回
    }
    if (!std::isfinite(h) || h <= 0.0) {
        throw std::invalid_argument("KDE: bandwidth must be positive and finite");
    }
    return h;
}

}  // namespace density
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
