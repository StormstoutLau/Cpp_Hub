// =============================================================================
// transfer_entropy.hpp - KSG-CMI 迁移熵 + effective TE + 块置换 p 值
//                        (Phase 9 M3, PHASE9_SPEC v1.1 §5)
//
// 定义 (Schreiber 2000; Barnett 2009 Eq. 10):
//   TE(X→Y) = CMI(Y_t ; X_{t−delay}, …, X_{t−k·delay} | Y_{t−delay}, …, Y_{t−l·delay})
//           = I(Y_t ; X_past | Y_past)
// 内部固定调用顺序 CMI(target_t, source_past, target_past) (spec §5.1),
// 避免 source/target 写反; 输出 nats.
//
// 字段语义:
//   te            KSG 条件 MI 估计 (nats), 不截断负值
//   te_effective  te − mean(TE(source 随机置换)) (Marschinski-Kantz 2002)
//                 仅诊断量, 不作推断依据 (审计 B1: 零假设下符号随种子波动)
//   te_bits       te / ln2
//   p_value       块置换零分布 (Dimpfl-Peter 2013, circular block shuffle
//                 只打乱 source, 保两序列各自自相关), 与 n_shuffles 独立抽样
//                 (审计 B3): p = (1 + #{null >= te}) / (1 + n_block_boot)
//   n_effective   T − max(k,l)·delay
//
// 不提供 "因果" 字段: TE 是条件 MI 的估计, 不是因果发现工具
// (James-Barnett-Crutchfield 2016; Smirnov 2013).
// 复杂度 O(N² · (k+l)) 每次 TE 评估 (暴力 KSG, D9-3).
// =============================================================================

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

#include "cpphub/core/constants.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/information/ksg.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace information {

/// 块置换的块长规则 (Q11)
enum class BlockLengthRule {
    kSqrt,         ///< ceil(√T) (v1.9 默认口径)
    kPolitisWhite  ///< Politis-White (2004) 自动块长
};

/// 迁移熵配置
struct TransferEntropyConfig {
    std::size_t source_history = 1;  ///< 源嵌入阶 k
    std::size_t target_history = 1;  ///< 目标嵌入阶 l
    std::size_t delay = 1;           ///< 嵌入间隔
    std::size_t ksg_neighbors = 3;   ///< KSG 近邻数
    bool standardize = true;         ///< 原始序列逐条 z-score (ddof=1)
    std::size_t n_shuffles = 100;    ///< effective TE 的随机置换次数
    std::size_t n_block_boot = 100;  ///< p 值的块置换次数 (与 n_shuffles 独立)
    std::size_t block_length = 0;    ///< >0: 显式块长 (优先级最高); 0: 按 block_rule
    /// Q11: block_length=0 时的块长规则
    ///   kSqrt          ceil(√T) (v1.9 默认, 保持既有 p 值锚不变)
    ///   kPolitisWhite  Politis-White (2004) 自动块长 (强持久序列推荐)
    BlockLengthRule block_rule = BlockLengthRule::kSqrt;
    std::uint64_t seed = 0;
};

/// 迁移熵结果
struct TransferEntropyResult {
    Real te = 0.0;            ///< nats
    Real te_effective = 0.0;  ///< te − mean(shuffled te)
    Real te_bits = 0.0;       ///< te / ln2
    Real p_value = 1.0;       ///< 块置换零分布
    std::size_t n_effective = 0;
    std::size_t block_length_used = 0;  ///< p 值实际使用的块长 (Q11 诊断)
};

namespace detail {

/// 嵌入后的连续坐标块; 行布局 [X_past(k), Y_t, Y_past(l)] 使 CMI 的
/// 四个空间 (联合 / x-z / y-z / z) 都对应连续坐标区间
struct EmbeddedSeries {
    std::vector<Real> points;  ///< n × dim
    Size n = 0;
    Size k = 0;
    Size l = 0;
    Size dim = 0;
};

inline EmbeddedSeries build_embedding(std::span<const Real> source,
                                      std::span<const Real> target, Size k, Size l,
                                      Size delay, Size start, Size n_eff) {
    EmbeddedSeries e;
    e.n = n_eff;
    e.k = k;
    e.l = l;
    e.dim = k + l + 1;
    e.points.resize(n_eff * e.dim);
    for (Size i = 0; i < n_eff; ++i) {
        const Size t = start + i;
        Real* row = &e.points[i * e.dim];
        for (Size j = 0; j < k; ++j) {
            row[j] = source[t - (j + 1) * delay];
        }
        row[k] = target[t];
        for (Size j = 0; j < l; ++j) {
            row[k + 1 + j] = target[t - (j + 1) * delay];
        }
    }
    return e;
}

/// 坐标区间 [c0, c1) 上的 Chebyshev 距离
inline Real cheb_range(const Real* a, const Real* b, Size c0, Size c1) {
    Real d = 0.0;
    for (Size c = c0; c < c1; ++c) {
        d = std::max(d, std::fabs(a[c] - b[c]));
    }
    return d;
}

/// Frenzel-Pompe CMI (多维嵌入): ψ(k) − ⟨ψ(nxz+1)+ψ(nyz+1)−ψ(nz+1)⟩
/// x = Y_t, y = X_past, z = Y_past (spec §5.1 固定调用顺序)
inline Real cmi_embedded(const EmbeddedSeries& e, Size k_nn) {
    const Size n = e.n;
    const Size dim = e.dim;
    const Size k = e.k;
    std::vector<Real> buf(k_nn, std::numeric_limits<Real>::infinity());
    Real sum = 0.0;
    for (Size i = 0; i < n; ++i) {
        const Real* pi = &e.points[i * dim];

        // 联合空间第 k 近邻距离 (排除自身)
        std::fill(buf.begin(), buf.end(), std::numeric_limits<Real>::infinity());
        for (Size j = 0; j < n; ++j) {
            if (j == i) {
                continue;
            }
            const Real d = cheb_range(pi, &e.points[j * dim], 0, dim);
            if (d < buf[k_nn - 1]) {
                Size pos = k_nn - 1;
                while (pos > 0 && buf[pos - 1] > d) {
                    buf[pos] = buf[pos - 1];
                    --pos;
                }
                buf[pos] = d;
            }
        }
        const Real radius = std::nextafter(buf[k_nn - 1], 0.0);

        // 三个边际计数共用同一半径; 计数含自身, 故 ψ(n·) ≡ ψ((n−1)+1)
        Size nxz = 0;
        Size nyz = 0;
        Size nz = 0;
        for (Size j = 0; j < n; ++j) {
            const Real* pj = &e.points[j * dim];
            const Real dz = cheb_range(pi, pj, k + 1, dim);
            const Real dy = std::fabs(pi[k] - pj[k]);
            const Real dx = cheb_range(pi, pj, 0, k);
            if (std::max(dy, dz) <= radius) {
                ++nxz;
            }
            if (std::max(dx, dz) <= radius) {
                ++nyz;
            }
            if (dz <= radius) {
                ++nz;
            }
        }
        sum += digamma(static_cast<Real>(nxz)) +
               digamma(static_cast<Real>(nyz)) -
               digamma(static_cast<Real>(nz));
    }
    return digamma(static_cast<Real>(k_nn)) - sum / static_cast<Real>(n);
}

/// Circular block shuffle: 从随机起点按块长环绕采样, 保块内自相关
inline void circular_block_shuffle(std::span<const Real> src, std::vector<Real>& out,
                                   Size block_len, std::mt19937_64& rng) {
    const Size n = src.size();
    out.resize(n);
    std::uniform_int_distribution<Size> pick(0, n - 1);
    Size pos = 0;
    while (pos < n) {
        const Size s = pick(rng);
        for (Size j = 0; j < block_len && pos < n; ++j, ++pos) {
            out[pos] = src[(s + j) % n];
        }
    }
}

/// Politis-White (2004) + Patton-Politis-White (2009) 修正的 circular
/// bootstrap 最优块长 (flat-top 核; 与 arch 8.0 `optimal_block_length`
/// circular 分支同算法, 逐步镜像其 autocovariance 截断与核权重)。
/// 返回未取整的实数块长; 调用方 ceil 并 clamp 到 [1, n]。
inline Real politis_white_block_length_raw(std::span<const Real> x) {
    const Size nobs = x.size();
    Real mean = 0.0;
    for (Real v : x) {
        mean += v;
    }
    mean /= static_cast<Real>(nobs);

    std::vector<Real> eps(nobs);
    for (Size i = 0; i < nobs; ++i) {
        eps[i] = x[i] - mean;
    }

    const Real n = static_cast<Real>(nobs);
    const Real b_max = std::ceil(std::min(3.0 * std::sqrt(n), n / 3.0));
    const Size kn = std::max<Size>(5, static_cast<Size>(std::log10(n)));
    const Size m_max = static_cast<Size>(std::ceil(std::sqrt(n))) + kn;
    const Real cv = 2.0 * std::sqrt(std::log10(n) / n);

    std::vector<Real> acv(m_max + 1, 0.0);
    std::vector<Real> abs_acorr(m_max + 1, 0.0);
    Size opt_m = 0;
    bool has_opt = false;
    for (Size i = 0; i <= m_max; ++i) {
        const Size i1 = i + 1;
        Real v1 = 0.0;
        Real v2 = 0.0;
        Real cross = 0.0;
        for (Size j = i1; j < nobs; ++j) {
            v1 += eps[j] * eps[j];
        }
        for (Size j = 0; j + i1 < nobs; ++j) {
            v2 += eps[j] * eps[j];
        }
        for (Size j = i; j < nobs; ++j) {
            cross += eps[j] * eps[j - i];
        }
        acv[i] = cross / n;
        const Real denom = std::sqrt(v1 * v2);
        abs_acorr[i] = (denom > 0.0) ? std::fabs(cross) / denom : 0.0;
        if (!has_opt && i >= kn) {
            bool all_small = true;
            for (Size j = i - kn; j < i; ++j) {
                if (!(abs_acorr[j] < cv)) {
                    all_small = false;
                    break;
                }
            }
            if (all_small) {
                opt_m = i - kn;
                has_opt = true;
            }
        }
    }

    Size m = has_opt ? 2 * std::max<Size>(opt_m, 1) : m_max;
    m = std::min(m, m_max);

    Real g = 0.0;
    Real lr_acv = acv[0];
    const Real m_real = static_cast<Real>(m);
    for (Size k = 1; k <= m; ++k) {
        const Real ratio = static_cast<Real>(k) / m_real;
        const Real lam = (ratio <= 0.5) ? 1.0 : 2.0 * (1.0 - ratio);
        g += 2.0 * lam * static_cast<Real>(k) * acv[k];
        lr_acv += 2.0 * lam * acv[k];
    }
    const Real d_cb = (4.0 / 3.0) * lr_acv * lr_acv;
    if (!(d_cb > 0.0)) {
        return 1.0;  // 常数序列 → 无自相关
    }
    const Real b = std::cbrt((2.0 * g * g) / d_cb) * std::cbrt(n);
    return std::min(b, b_max);
}

}  // namespace detail

/// @brief Politis-White (2004) 最优 circular block bootstrap 块长 (整数, ≥ 1)
///
/// 与 arch 8.0 `optimal_block_length(..., "circular")` 同算法; 差异仅在
/// 末步 ceil + 下限 1 (arch 对 i.i.d. 可返回 0.0 浮点值)。
/// 强持久序列下该块长显著大于 ceil(√T), 使块置换零分布保留与原序列
/// 相同量级的自相关 (Q11; 见 TransferEntropyConfig::block_rule)。
///
/// @throws std::invalid_argument n < 20 (自相关估计不可靠) / 非有限值
inline Size politis_white_block_length(std::span<const Real> x) {
    if (x.size() < 20) {
        throw std::invalid_argument(
            "PW block length: requires n >= 20 observations");
    }
    detail::require_finite(x);
    const Real b = detail::politis_white_block_length_raw(x);
    return std::max<Size>(1, static_cast<Size>(std::ceil(b)));
}

/// @brief KSG-CMI 迁移熵 TE(X→Y), nats
/// @throws std::invalid_argument 长度不等 / 空序列 / 参数非法 /
///         非有限值 / 有效样本不足 (T − max(k,l)·delay < ksg_neighbors+1)
inline TransferEntropyResult transfer_entropy(std::span<const Real> source,
                                              std::span<const Real> target,
                                              const TransferEntropyConfig& config = {}) {
    if (source.size() != target.size()) {
        throw std::invalid_argument("TE: source and target must have equal length");
    }
    if (config.source_history < 1 || config.target_history < 1 || config.delay < 1 ||
        config.ksg_neighbors < 1) {
        throw std::invalid_argument(
            "TE: source_history, target_history, delay and ksg_neighbors must be >= 1");
    }
    detail::require_finite(source);
    detail::require_finite(target);

    const Size T = source.size();
    const Size start =
        std::max(config.source_history, config.target_history) * config.delay;
    if (T <= start) {
        throw std::invalid_argument("TE: insufficient observations for embedding");
    }
    const Size n_eff = T - start;
    if (n_eff < config.ksg_neighbors + 1) {
        throw std::invalid_argument("TE: insufficient observations for ksg_neighbors");
    }

    std::vector<Real> xs(source.begin(), source.end());
    std::vector<Real> ys(target.begin(), target.end());
    if (config.standardize) {
        detail::standardize_inplace(xs);
        detail::standardize_inplace(ys);
    }

    const auto build = [&](std::span<const Real> src) {
        return detail::build_embedding(src, std::span<const Real>(ys),
                                       config.source_history, config.target_history,
                                       config.delay, start, n_eff);
    };

    TransferEntropyResult result;
    result.n_effective = n_eff;
    result.te = detail::cmi_embedded(build(std::span<const Real>(xs)), config.ksg_neighbors);
    result.te_bits = result.te / LN_2;

    std::mt19937_64 rng(config.seed);

    // effective TE: 随机置换 source, 破坏全部依赖 (Marschinski-Kantz)
    if (config.n_shuffles > 0) {
        std::vector<Real> permuted = xs;
        Real sum = 0.0;
        for (Size s = 0; s < config.n_shuffles; ++s) {
            std::shuffle(permuted.begin(), permuted.end(), rng);
            sum += detail::cmi_embedded(build(std::span<const Real>(permuted)),
                                        config.ksg_neighbors);
        }
        result.te_effective = result.te - sum / static_cast<Real>(config.n_shuffles);
    } else {
        result.te_effective = result.te;
    }

    // p 值: 块置换零分布 (只打乱 source, 与 n_shuffles 独立抽样)
    if (config.n_block_boot == 0) {
        result.p_value = 1.0;
    } else {
        Size block = config.block_length;
        if (block == 0) {
            if (config.block_rule == BlockLengthRule::kPolitisWhite && T >= 20) {
                block = politis_white_block_length(std::span<const Real>(xs));
            } else {
                block = static_cast<Size>(
                    std::ceil(std::sqrt(static_cast<Real>(T))));
            }
        }
        block = std::min(std::max<Size>(block, 1), T);
        result.block_length_used = block;

        std::vector<Real> surrogate;
        Size exceed = 0;
        for (Size b = 0; b < config.n_block_boot; ++b) {
            detail::circular_block_shuffle(std::span<const Real>(xs), surrogate, block,
                                           rng);
            const Real null_te = detail::cmi_embedded(
                build(std::span<const Real>(surrogate)), config.ksg_neighbors);
            if (null_te >= result.te) {
                ++exceed;
            }
        }
        result.p_value =
            (1.0 + static_cast<Real>(exceed)) / (1.0 + static_cast<Real>(config.n_block_boot));
    }

    return result;
}

}  // namespace information
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
