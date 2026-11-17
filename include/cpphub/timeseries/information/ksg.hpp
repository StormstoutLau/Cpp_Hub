// =============================================================================
// ksg.hpp - KSG Algorithm 1 互信息 + Frenzel-Pompe 条件互信息 (Phase 9 M2,
//           PHASE9_SPEC v1.1 §4.2)
//
// Kraskov-Stoegbauer-Grassberger 2004:
//   MI   = ψ(k) + ψ(N) − ⟨ψ(n_x+1) + ψ(n_y+1)⟩                    (Eq. 8)
//   CMI  = ψ(k) − ⟨ψ(n_xz+1) + ψ(n_yz+1) − ψ(n_z+1)⟩        (Frenzel-Pompe)
//
// 约定 (冻结, 对齐 sklearn _compute_mi_cc + probe_ksg.py K0 逐位一致):
//   - 联合空间 Chebyshev (max norm): (x,y) 或 (x,y,z)
//   - ε(i) = 到第 k 近邻的联合距离 (排除自身, 对应 sklearn kneighbors()
//     在 X=None 时不把查询点算作自身近邻)
//   - 边际计数用 nextafter(ε, 0) 上的 <= 比较 (严格小于 ε 的等价实现;
//     ε=0 的重复点场景与 sklearn query_radius 语义一致), 计数减 1 (自身)
//   - 返回 nats; 不截断负值 (sklearn 的 max(0,·) 语义不在本库)
//   - standardize: 逐维 z-score (ddof=1); jitter: 相对 mean|x| 的 N(0,1)
//     噪声, 单条 RNG 流按 x→y→z 顺序注入 (seed 固定, 用于 ties)
//   - P2-4 (2026-11-17): 联合空间 kNN 与边际半径计数改由 Chebyshev k-d tree
//     加速 (density/kdtree.hpp), 典型 O(N log N); 所有距离/计数与暴力参考
//     逐位一致 (树只改访问路径); 高维退化时最坏仍 O(N²)
//   - 文档化缩放敏感性: 混合单位数据务必 standardize=true (调研风险 ②)
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

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/density/kdtree.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace information {

/// KSG 估计器配置
struct KsgConfig {
    std::size_t k = 3;             ///< 近邻数 (k >= 1)
    bool standardize = true;       ///< 逐维 z-score (ddof=1)
    double jitter_scale = 1e-10;   ///< 相对 mean|x| 的噪声幅度; 0 = 关闭
    std::uint64_t seed = 0;        ///< jitter RNG 种子
};

namespace detail {

/// digamma ψ(x) (x > 0): 递推至 x >= 16 后渐近级数到 x^-10
/// (与 scipy.special.digamma 实测最大绝对误差 1.2e-15)
inline Real digamma(Real x) {
    if (!(x > 0.0)) {
        throw std::invalid_argument("KSG: digamma requires x > 0");
    }
    Real result = 0.0;
    while (x < 16.0) {
        result -= 1.0 / x;
        x += 1.0;
    }
    const Real inv = 1.0 / x;
    const Real inv2 = inv * inv;
    result += std::log(x) - 0.5 * inv -
              inv2 * (1.0 / 12.0 -
                      inv2 * (1.0 / 120.0 -
                              inv2 * (1.0 / 252.0 -
                                      inv2 * (1.0 / 240.0 - inv2 * (1.0 / 132.0)))));
    return result;
}

inline void require_finite(std::span<const Real> v) {
    for (Real u : v) {
        if (!std::isfinite(u)) {
            throw std::invalid_argument("KSG: input contains non-finite value");
        }
    }
}

/// 逐维 z-score (ddof=1); 常量输入 (sd=0) 保持原样 (距离全 0 的退化 ties)
inline void standardize_inplace(std::vector<Real>& v) {
    const Real n = static_cast<Real>(v.size());
    Real mean = 0.0;
    for (Real u : v) {
        mean += u;
    }
    mean /= n;
    Real ss = 0.0;
    for (Real u : v) {
        const Real d = u - mean;
        ss += d * d;
    }
    const Real sd = std::sqrt(ss / (n - 1.0));
    if (!(sd > 0.0) || !std::isfinite(sd)) {
        return;
    }
    for (Real& u : v) {
        u = (u - mean) / sd;
    }
}

/// 相对 mean|v| 的 N(0,1) 噪声 (共用一条 RNG 流, 保证各维独立且可复现)
inline void add_jitter_inplace(std::vector<Real>& v, Real scale, std::mt19937_64& rng) {
    const Real n = static_cast<Real>(v.size());
    Real mabs = 0.0;
    for (Real u : v) {
        mabs += std::fabs(u);
    }
    mabs /= n;
    const Real amp = scale * mabs;
    if (!(amp > 0.0) || !std::isfinite(amp)) {
        return;
    }
    std::normal_distribution<Real> nd(0.0, 1.0);
    for (Real& u : v) {
        u += amp * nd(rng);
    }
}

/// 升序样本上的精确 1D 半径计数 (含查询点自身): |v − center| <= radius
/// 两侧分别用同一减法表达式二分 (而非 center ± radius 的舍入边界), 与暴力
/// `fabs(v − center) <= radius` 逐位一致; O(log n)
inline Size count_within_sorted_1d(const std::vector<Real>& sorted, Real center,
                                   Real radius) {
    // 左侧 v <= center: 谓词 (center − v) <= radius 随 v 递增由 false 变 true
    const Size hi = static_cast<Size>(
        std::upper_bound(sorted.begin(), sorted.end(), center) - sorted.begin());
    Size lo = 0;
    Size r = hi;
    while (lo < r) {
        const Size m = lo + (r - lo) / 2;
        if (center - sorted[m] <= radius) {
            r = m;
        } else {
            lo = m + 1;
        }
    }
    const Size left_count = hi - lo;
    // 右侧 v > center: 谓词 (v − center) <= radius 随 v 递增由 true 变 false
    Size l2 = hi;
    Size r2 = sorted.size();
    while (l2 < r2) {
        const Size m = l2 + (r2 - l2) / 2;
        if (sorted[m] - center <= radius) {
            l2 = m + 1;
        } else {
            r2 = m;
        }
    }
    return left_count + (l2 - hi);
}

/// 公共校验: 长度一致 / k 合法 / N >= k+1 / 值有限
inline void validate_knn_inputs(std::span<const Real> x, std::span<const Real> y,
                                const KsgConfig& config) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("KSG: x and y must have equal length");
    }
    if (config.k < 1) {
        throw std::invalid_argument("KSG: k must be >= 1");
    }
    if (x.size() < config.k + 1) {
        throw std::invalid_argument("KSG: sample size must be >= k + 1");
    }
    require_finite(x);
    require_finite(y);
}

}  // namespace detail

/// @brief KSG Algorithm 1 互信息 (nats)
/// @throws std::invalid_argument 长度不等 / k < 1 / N < k+1 / 非有限值
inline Real ksg_mutual_info(std::span<const Real> x, std::span<const Real> y,
                            const KsgConfig& config = {}) {
    detail::validate_knn_inputs(x, y, config);
    const Size n = x.size();

    std::vector<Real> xs(x.begin(), x.end());
    std::vector<Real> ys(y.begin(), y.end());
    if (config.standardize) {
        detail::standardize_inplace(xs);
        detail::standardize_inplace(ys);
    }
    if (config.jitter_scale > 0.0) {
        std::mt19937_64 rng(config.seed);
        detail::add_jitter_inplace(xs, config.jitter_scale, rng);
        detail::add_jitter_inplace(ys, config.jitter_scale, rng);
    }

    // P2-4: 联合 (x,y) k-d tree — k+1 含自身 ≡ 排除自身的第 k 近邻
    std::vector<Real> joint(n * 2);
    for (Size i = 0; i < n; ++i) {
        joint[2 * i] = xs[i];
        joint[2 * i + 1] = ys[i];
    }
    static constexpr Size kJointDims[2] = {0, 1};
    const density::KdTree joint_tree(std::span<const Real>(joint), n, 2,
                                     std::span<const Size>(kJointDims, 2));
    std::vector<Real> sx = xs;
    std::vector<Real> sy = ys;
    std::sort(sx.begin(), sx.end());
    std::sort(sy.begin(), sy.end());

    Real sum_x = 0.0;
    Real sum_y = 0.0;
    for (Size i = 0; i < n; ++i) {
        const Real eps = joint_tree.kth_distance_including_self(i, config.k + 1);
        const Real radius = std::nextafter(eps, 0.0);
        const Size nx = detail::count_within_sorted_1d(sx, xs[i], radius) - 1;
        const Size ny = detail::count_within_sorted_1d(sy, ys[i], radius) - 1;
        sum_x += detail::digamma(static_cast<Real>(nx) + 1.0);
        sum_y += detail::digamma(static_cast<Real>(ny) + 1.0);
    }

    const Real n_real = static_cast<Real>(n);
    return detail::digamma(n_real) + detail::digamma(static_cast<Real>(config.k)) -
           sum_x / n_real - sum_y / n_real;
}

/// @brief Frenzel-Pompe 条件互信息 I(X;Y|Z) (nats)
/// @throws std::invalid_argument 长度不等 / k < 1 / N < k+1 / 非有限值
inline Real ksg_conditional_mutual_info(std::span<const Real> x, std::span<const Real> y,
                                        std::span<const Real> z,
                                        const KsgConfig& config = {}) {
    detail::validate_knn_inputs(x, y, config);
    if (x.size() != z.size()) {
        throw std::invalid_argument("KSG: x, y and z must have equal length");
    }
    detail::require_finite(z);
    const Size n = x.size();

    std::vector<Real> xs(x.begin(), x.end());
    std::vector<Real> ys(y.begin(), y.end());
    std::vector<Real> zs(z.begin(), z.end());
    if (config.standardize) {
        detail::standardize_inplace(xs);
        detail::standardize_inplace(ys);
        detail::standardize_inplace(zs);
    }
    if (config.jitter_scale > 0.0) {
        std::mt19937_64 rng(config.seed);
        detail::add_jitter_inplace(xs, config.jitter_scale, rng);
        detail::add_jitter_inplace(ys, config.jitter_scale, rng);
        detail::add_jitter_inplace(zs, config.jitter_scale, rng);
    }

    // P2-4: 联合 (x,y,z) + 边际 (x,z)/(y,z) k-d tree; z 边际走有序数组二分
    std::vector<Real> joint(n * 3);
    for (Size i = 0; i < n; ++i) {
        joint[3 * i] = xs[i];
        joint[3 * i + 1] = ys[i];
        joint[3 * i + 2] = zs[i];
    }
    static constexpr Size kDimsXYZ[3] = {0, 1, 2};
    static constexpr Size kDimsXZ[2] = {0, 2};
    static constexpr Size kDimsYZ[2] = {1, 2};
    const std::span<const Real> pts(joint);
    const density::KdTree tree_xyz(pts, n, 3, std::span<const Size>(kDimsXYZ, 3));
    const density::KdTree tree_xz(pts, n, 3, std::span<const Size>(kDimsXZ, 2));
    const density::KdTree tree_yz(pts, n, 3, std::span<const Size>(kDimsYZ, 2));
    std::vector<Real> sz = zs;
    std::sort(sz.begin(), sz.end());

    Real sum = 0.0;
    for (Size i = 0; i < n; ++i) {
        const Real eps = tree_xyz.kth_distance_including_self(i, config.k + 1);
        const Real radius = std::nextafter(eps, 0.0);
        const Size nxz = tree_xz.count_within(i, radius) - 1;
        const Size nyz = tree_yz.count_within(i, radius) - 1;
        const Size nz = detail::count_within_sorted_1d(sz, zs[i], radius) - 1;
        sum += detail::digamma(static_cast<Real>(nxz) + 1.0) +
               detail::digamma(static_cast<Real>(nyz) + 1.0) -
               detail::digamma(static_cast<Real>(nz) + 1.0);
    }

    return detail::digamma(static_cast<Real>(config.k)) - sum / static_cast<Real>(n);
}

}  // namespace information
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
