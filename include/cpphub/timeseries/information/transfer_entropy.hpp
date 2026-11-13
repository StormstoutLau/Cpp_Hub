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

/// 迁移熵配置
struct TransferEntropyConfig {
    std::size_t source_history = 1;  ///< 源嵌入阶 k
    std::size_t target_history = 1;  ///< 目标嵌入阶 l
    std::size_t delay = 1;           ///< 嵌入间隔
    std::size_t ksg_neighbors = 3;   ///< KSG 近邻数
    bool standardize = true;         ///< 原始序列逐条 z-score (ddof=1)
    std::size_t n_shuffles = 100;    ///< effective TE 的随机置换次数
    std::size_t n_block_boot = 100;  ///< p 值的块置换次数 (与 n_shuffles 独立)
    std::size_t block_length = 0;    ///< 0 => ceil(sqrt(T))
    std::uint64_t seed = 0;
};

/// 迁移熵结果
struct TransferEntropyResult {
    Real te = 0.0;            ///< nats
    Real te_effective = 0.0;  ///< te − mean(shuffled te)
    Real te_bits = 0.0;       ///< te / ln2
    Real p_value = 1.0;       ///< 块置换零分布
    std::size_t n_effective = 0;
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

}  // namespace detail

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
            block = static_cast<Size>(
                std::ceil(std::sqrt(static_cast<Real>(T))));
        }
        block = std::min(std::max<Size>(block, 1), T);

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
