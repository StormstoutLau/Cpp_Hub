// =============================================================================
// bandwidth.hpp - KDE 带宽类型与规则解析 (Phase 9 M1, PHASE9_SPEC v1.1 §3.1)
//
// 同名不同物的带宽语义在此显式区分 (调研 A1-A25):
//   kScott           1.059 · A · n^(-1/5)   (statsmodels bw_scott / R bw.nrd)
//   kSilverman       0.9   · A · n^(-1/5)   (statsmodels bw_silverman / R bw.nrd0)
//   kNormalReference C     · A · n^(-1/5)   (statsmodels 默认; C 由 probe P4 反解)
//   kSj              Sheather-Jones (1991) solve-the-equation plug-in
//                     (R bw.SJ(method="ste") 同算法: 分箱 pair counts + phi4/phi6)
//   kIsj             Botev-Grotowski-Kroese (2010) improved SJ (diffusion)
//                     (KDEpy improved_sheather_jones 同算法: 2^10 网格 + DCT-II)
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
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace density {

/// 带宽规则 (Scott/Silverman/NR 数值仅差 ~0.02%, 保留独立枚举以便未来换核;
/// kSj/kIsj 为数据驱动 plug-in 选择器, 计算成本 O(n²)/O(1024²))
enum class BandwidthRule {
    kUser,
    kScott,
    kSilverman,
    kNormalReference,
    kSj,
    kIsj
};

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

inline constexpr Real kPi = 3.141592653589793238462643383279502884;

/// Sheather-Jones (1991) solve-the-equation plug-in 带宽
/// 逐步镜像 R `bw.SJ(method="ste")`:
///   - nb=1000 分箱 pair counts (R bw_pair_cnts/bw_den_binned; O(n + nb²))
///   - φ4(h)=Σ 对距离核泛函, φ6(h)=六阶导泛函 (Scott 1992 (6.67)/(6.69))
///   - scale=min(sd, IQR/1.349); a=1.24·scale·n^(-1/7); b=1.23·scale·n^(-1/9)
///   - TD=-φ6(b); alph2=1.357·(φ4(a)/TD)^(1/7)
///   - 解 h = (c1/φ4(alph2·h^(5/7)))^(1/5), c1=1/(2√π·n)
/// 与 R 的差异仅在求根精度 (R 默认 tol=0.1·lower 较宽松; 本实现求到机器精度).
inline Real sj_bandwidth(std::span<const Real> sample, Real scale) {
    const Size n = sample.size();
    constexpr Size kNb = 1000;  // R bw.SJ 默认 nb
    const Real n_real = static_cast<Real>(n);

    Real xmin = sample[0];
    Real xmax = sample[0];
    for (Real v : sample) {
        xmin = std::min(xmin, v);
        xmax = std::max(xmax, v);
    }
    const Real dd = (xmax - xmin) * 1.01 / static_cast<Real>(kNb);
    if (!(dd > 0.0) || !std::isfinite(dd)) {
        throw std::invalid_argument("KDE: data are constant in bandwidth calculation");
    }

    // R bw_pair_cnts 的分箱: trunc(|x|/dd)·sign(x) 后平移, tabulate 到 nb 箱
    std::vector<long long> bins(n);
    long long min_bin = std::numeric_limits<long long>::max();
    for (Size i = 0; i < n; ++i) {
        bins[i] = static_cast<long long>(std::trunc(sample[i] / dd));
        min_bin = std::min(min_bin, bins[i]);
    }
    std::vector<Real> w(kNb, 0.0);
    for (Size i = 0; i < n; ++i) {
        const long long idx = bins[i] - min_bin;  // 0-based (R 的 1-based 减 1)
        if (idx >= 0 && idx < static_cast<long long>(kNb)) {
            w[static_cast<Size>(idx)] += 1.0;
        }
    }
    std::vector<Real> cnt(kNb, 0.0);
    for (Size ii = 0; ii < kNb; ++ii) {
        const Real wi = w[ii];
        if (wi == 0.0) {
            continue;
        }
        cnt[0] += wi * (wi - 1.0);
        for (Size jj = 0; jj < ii; ++jj) {
            if (w[jj] != 0.0) {
                cnt[ii - jj] += wi * w[jj];
            }
        }
    }
    cnt[0] *= 0.5;

    // φ_{2r}(h): 高斯核 r 阶导数泛函 (R bw_phi4/bw_phi6)
    const auto phi = [&](Real h, int order) {
        Real sum = 0.0;
        for (Size k = 0; k < cnt.size(); ++k) {
            const Real t = static_cast<Real>(k) * dd / h;
            const Real d2 = t * t;
            if (d2 >= 1000.0) {
                break;
            }
            Real term = 0.0;
            if (order == 4) {
                term = std::exp(-d2 / 2.0) * (d2 * d2 - 6.0 * d2 + 3.0);
            } else {
                term = std::exp(-d2 / 2.0) *
                       (d2 * d2 * d2 - 15.0 * d2 * d2 + 45.0 * d2 - 15.0);
            }
            sum += term * cnt[k];
        }
        if (order == 4) {
            sum = 2.0 * sum + 3.0 * n_real;
        } else {
            sum = 2.0 * sum - 15.0 * n_real;
        }
        return sum / (n_real * (n_real - 1.0) * std::pow(h, order + 1)) /
               std::sqrt(2.0 * kPi);
    };

    const Real a = 1.24 * scale * std::pow(n_real, -1.0 / 7.0);
    const Real b = 1.23 * scale * std::pow(n_real, -1.0 / 9.0);
    const Real c1 = 1.0 / (2.0 * std::sqrt(kPi) * n_real);
    const Real td = -phi(b, 6);
    if (!(td > 0.0) || !std::isfinite(td)) {
        throw std::invalid_argument("KDE: sample is too sparse to find TD");
    }
    const Real alph2 = 1.357 * std::pow(phi(a, 4) / td, 1.0 / 7.0);
    if (!std::isfinite(alph2)) {
        throw std::invalid_argument("KDE: sample is too sparse to find alph2");
    }

    const Real hmax = 1.144 * scale * std::pow(n_real, -0.2);
    Real lower = 0.1 * hmax;
    Real upper = hmax;
    const auto fsd = [&](Real h) {
        return std::pow(c1 / phi(alph2 * std::pow(h, 5.0 / 7.0), 4), 0.2) - h;
    };
    if (!(std::isfinite(fsd(lower)) && std::isfinite(fsd(upper)))) {
        throw std::invalid_argument("KDE: SJ bandwidth functional is not finite");
    }
    int itry = 1;
    while (fsd(lower) * fsd(upper) > 0.0) {
        if (itry > 99) {
            throw std::invalid_argument("KDE: no solution in the specified range");
        }
        if (itry % 2) {
            upper *= 1.2;
        } else {
            lower /= 1.2;
        }
        ++itry;
    }

    // 二分求根 (机器精度; 与 R uniroot 的解一致到其容差内)
    Real lo = lower;
    Real hi = upper;
    Real f_lo = fsd(lo);
    for (int iter = 0; iter < 200; ++iter) {
        const Real mid = 0.5 * (lo + hi);
        const Real f_mid = fsd(mid);
        if (f_mid == 0.0 || (hi - lo) <= 1e-15 * std::max(1.0, std::fabs(mid))) {
            lo = mid;
            hi = mid;
            break;
        }
        if ((f_lo < 0.0) == (f_mid < 0.0)) {
            lo = mid;
            f_lo = f_mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

/// Botev-Grotowski-Kroese (2010) improved Sheather-Jones (diffusion) 带宽
/// 逐步镜像 KDEpy `improved_sheather_jones`:
///   - 2^10 网格 (autogrid: min−max(0.5R,6) … max+max(0.5R,6)) + 线性分箱
///   - DCT-II 系数 a, a2=a[1:]², I_sq=(1..N−1)²; N_unique clamp [50,1050]
///   - 固定点迭代 t=ξγ⁵(t) (ell=7, s=6..2), brentq → t*; h=√t*·R
/// DCT 直接 O(M²) 求值 (M=1024, ~1M 次 cos; 选择器仅在构造时调用一次).
inline Real isj_bandwidth(std::span<const Real> sample) {
    const Size n = sample.size();
    constexpr Size kGrid = 1024;
    Real xmin = sample[0];
    Real xmax = sample[0];
    for (Real v : sample) {
        xmin = std::min(xmin, v);
        xmax = std::max(xmax, v);
    }
    const Real range = xmax - xmin;
    if (!(range > 0.0)) {
        throw std::invalid_argument("KDE: data are constant in bandwidth calculation");
    }
    const Real border = std::max(0.5 * range, 6.0);
    const Real lo = xmin - border;
    const Real hi = xmax + border;
    const Real dx = (hi - lo) / static_cast<Real>(kGrid - 1);
    const Real w = 1.0 / static_cast<Real>(n);
    std::vector<Real> hist(kGrid, 0.0);
    for (Real v : sample) {
        const Real t = (v - lo) / dx;
        const long long i0 = static_cast<long long>(std::floor(t));
        if (i0 < 0 || i0 >= static_cast<long long>(kGrid - 1)) {
            continue;  // 数据在网格内, 仅防御
        }
        const Real frac = t - static_cast<Real>(i0);
        hist[static_cast<Size>(i0)] += (1.0 - frac) * w;
        hist[static_cast<Size>(i0) + 1] += frac * w;
    }

    // 唯一值数: KDEpy 的 _root 只用 clamped N 计算容差, 而 _fixed_point 通过
    // args 收到**原始唯一值数** (N 在 _root 内钳制但 args 已构造) — 重离散
    // 数据下两者差异显著, 此处忠实镜像 (ar200: h 0.0052 vs 0.890)
    std::vector<Real> sorted(sample.begin(), sample.end());
    std::sort(sorted.begin(), sorted.end());
    const Size uniq = static_cast<Size>(std::unique(sorted.begin(), sorted.end()) -
                                        sorted.begin());
    const Real n_fp = static_cast<Real>(uniq);  // 固定点内部使用
    const Real nu = static_cast<Real>(std::min<Size>(std::max<Size>(uniq, 50), 1050));

    // DCT-II: a_k = 2 Σ_j hist_j cos(π k (2j+1) / (2M))
    std::vector<Real> coef(kGrid, 0.0);
    for (Size k = 0; k < kGrid; ++k) {
        Real sum = 0.0;
        const Real ck = kPi * static_cast<Real>(k) / static_cast<Real>(2 * kGrid);
        for (Size j = 0; j < kGrid; ++j) {
            sum += hist[j] * std::cos(ck * static_cast<Real>(2 * j + 1));
        }
        coef[k] = 2.0 * sum;
    }
    std::vector<Real> a2(kGrid - 1);
    for (Size k = 1; k < kGrid; ++k) {
        a2[k - 1] = coef[k] * coef[k];
    }

    const auto fsum = [&](Real t, int ell) {
        Real sum = 0.0;
        for (Size i = 0; i + 1 < kGrid; ++i) {
            const Real isq = static_cast<Real>((i + 1) * (i + 1));
            sum += std::pow(isq, ell) * a2[i] *
                   std::exp(-isq * kPi * kPi * t);
        }
        return 0.5 * std::pow(kPi, 2.0 * ell) * sum;
    };

    const auto fixed_point = [&](Real t) {
        Real f = fsum(t, 7);
        if (!(f > 0.0)) {
            return -1.0;
        }
        for (int s = 6; s >= 2; --s) {
            Real odd = 1.0;
            for (int q = 1; q <= 2 * s; q += 2) {
                odd *= static_cast<Real>(q);
            }
            const Real k0 = odd / std::sqrt(2.0 * kPi);
            const Real cst = (1.0 + std::pow(0.5, s + 0.5)) / 3.0;
            const Real time =
                std::pow(2.0 * cst * k0 / (n_fp * f), 2.0 / (3.0 + 2.0 * s));
            f = fsum(time, s);
        }
        const Real t_opt = std::pow(2.0 * n_fp * std::sqrt(kPi) * f, -0.4);
        return t - t_opt;
    };

    Real tol = 10e-12 + 0.01 * (nu - 50.0) / 1000.0;
    while (fixed_point(tol) <= 0.0) {
        tol *= 2.0;
        if (tol >= 1.0) {
            throw std::invalid_argument(
                "KDE: ISJ bandwidth selection did not converge");
        }
    }
    Real lo_t = 0.0;
    Real hi_t = tol;
    for (int iter = 0; iter < 200; ++iter) {
        const Real mid = 0.5 * (lo_t + hi_t);
        if (fixed_point(mid) <= 0.0) {
            lo_t = mid;
        } else {
            hi_t = mid;
        }
        if ((hi_t - lo_t) <= 1e-15 * std::max(1.0, std::fabs(mid))) {
            break;
        }
    }
    const Real t_star = 0.5 * (lo_t + hi_t);
    return std::sqrt(t_star) * range;
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
        case BandwidthRule::kSj:
            h = detail::sj_bandwidth(sample, a);
            break;
        case BandwidthRule::kIsj:
            h = detail::isj_bandwidth(sample);
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
