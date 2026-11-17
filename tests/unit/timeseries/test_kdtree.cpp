// =============================================================================
// test_kdtree.cpp - Phase 9 P2-4: Chebyshev k-d tree (TDD 前置)
//
// 状态: P2-4 已实现 (include/cpphub/timeseries/density/kdtree.hpp),
//       被 ksg.hpp / transfer_entropy.hpp 用于加速 kNN 与半径计数.
//
// 判据:
//   - kth_distance_including_self / count_within 与独立暴力 oracle 逐位一致
//     (含重复点/并列距离; 1-5 维 + 列子集)
//   - 全重合数据: 距离恒 0, radius=0 计数 = n
//   - 边界: 空样本 / 空列集 / 列越界 / 非有限值 / 查询越界 / k 越界 /
//     radius<0 / 未构建
// =============================================================================

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cpphub/timeseries/density/kdtree.hpp"

namespace den = cpphub::v1::timeseries::density;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

// ---------------- 独立暴力 oracle (与 kdtree.hpp 无共享代码) ----------------

Real brute_cheb(const std::vector<Real>& pts, Size stride, const std::vector<Size>& dims,
                Size a, Size b) {
    Real d = 0.0;
    for (Size c : dims) {
        d = std::max(d, std::fabs(pts[a * stride + c] - pts[b * stride + c]));
    }
    return d;
}

Real brute_kth(const std::vector<Real>& pts, Size n, Size stride,
               const std::vector<Size>& dims, Size query, Size k) {
    std::vector<Real> dist;
    dist.reserve(n);
    for (Size j = 0; j < n; ++j) {
        dist.push_back(brute_cheb(pts, stride, dims, query, j));
    }
    std::nth_element(dist.begin(), dist.begin() + static_cast<std::ptrdiff_t>(k - 1),
                     dist.end());
    return dist[k - 1];
}

Size brute_count(const std::vector<Real>& pts, Size n, Size stride,
                 const std::vector<Size>& dims, Size query, Real radius) {
    Size c = 0;
    for (Size j = 0; j < n; ++j) {
        if (brute_cheb(pts, stride, dims, query, j) <= radius) {
            ++c;
        }
    }
    return c;
}

// 重复点友好数据: 整数 LCG → 101 个网格值 (n > 101 时必然出现重复坐标)
std::vector<Real> grid_points(Size n, Size dim, std::uint64_t seed) {
    std::vector<Real> pts(n * dim);
    std::uint64_t s = seed;
    for (Size i = 0; i < n * dim; ++i) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        pts[i] = static_cast<Real>((s >> 33) % 101) / 10.0;
    }
    return pts;
}

std::vector<Size> all_dims(Size dim) {
    std::vector<Size> v(dim);
    for (Size c = 0; c < dim; ++c) {
        v[c] = c;
    }
    return v;
}

}  // namespace

TEST(KdTree, KthDistanceMatchesBruteForce) {
    constexpr Size n = 257;
    for (Size dim = 1; dim <= 5; ++dim) {
        const std::vector<Real> pts = grid_points(n, dim, 11 + dim);
        const std::vector<Size> dims = all_dims(dim);
        const den::KdTree tree(pts, n, dim, dims);
        for (Size k : {Size{1}, Size{2}, Size{3}, Size{7}, n}) {
            for (Size q = 0; q < n; ++q) {
                EXPECT_DOUBLE_EQ(tree.kth_distance_including_self(q, k),
                                 brute_kth(pts, n, dim, dims, q, k))
                    << "dim=" << dim << " k=" << k << " q=" << q;
            }
        }
    }
}

TEST(KdTree, CountWithinMatchesBruteForce) {
    constexpr Size n = 257;
    const Real inf = std::numeric_limits<Real>::infinity();
    for (Size dim = 1; dim <= 5; ++dim) {
        const std::vector<Real> pts = grid_points(n, dim, 31 + dim);
        const std::vector<Size> dims = all_dims(dim);
        const den::KdTree tree(pts, n, dim, dims);
        for (Real radius : {0.0, 0.05, 0.1, 0.25, 0.7, 2.0, inf}) {
            for (Size q = 0; q < n; ++q) {
                EXPECT_EQ(tree.count_within(q, radius),
                          brute_count(pts, n, dim, dims, q, radius))
                    << "dim=" << dim << " r=" << radius << " q=" << q;
            }
        }
    }
}

TEST(KdTree, SubsetDimsMatchesBruteForce) {
    // stride=4 但仅 {0,2} 参与距离: 验证列子集语义 (对应 KSG 的 (x,z) 边际)
    constexpr Size n = 200;
    constexpr Size stride = 4;
    const std::vector<Real> pts = grid_points(n, stride, 77);
    const std::vector<Size> dims = {0, 2};
    const den::KdTree tree(pts, n, stride, dims);
    for (Size q = 0; q < n; ++q) {
        EXPECT_DOUBLE_EQ(tree.kth_distance_including_self(q, 4),
                         brute_kth(pts, n, stride, dims, q, 4));
        EXPECT_EQ(tree.count_within(q, 0.3), brute_count(pts, n, stride, dims, q, 0.3));
    }
}

TEST(KdTree, AllIdenticalPoints) {
    constexpr Size n = 37;
    constexpr Size dim = 3;
    std::vector<Real> pts(n * dim, 1.5);
    const std::vector<Size> dims = all_dims(dim);
    const den::KdTree tree(pts, n, dim, dims);
    for (Size k : {Size{1}, Size{5}, n}) {
        EXPECT_DOUBLE_EQ(tree.kth_distance_including_self(0, k), 0.0);
    }
    EXPECT_EQ(tree.count_within(0, 0.0), n);
    EXPECT_EQ(tree.count_within(n - 1, 1e-12), n);
    EXPECT_EQ(tree.size(), n);
    EXPECT_EQ(tree.dims(), dim);
}

TEST(KdTree, MatchesBruteForceForAllKthOnSmallSample) {
    constexpr Size n = 9;
    constexpr Size dim = 2;
    const std::vector<Real> pts = grid_points(n, dim, 5);
    const std::vector<Size> dims = all_dims(dim);
    const den::KdTree tree(pts, n, dim, dims);
    for (Size k = 1; k <= n; ++k) {
        for (Size q = 0; q < n; ++q) {
            EXPECT_DOUBLE_EQ(tree.kth_distance_including_self(q, k),
                             brute_kth(pts, n, dim, dims, q, k));
        }
    }
}

TEST(KdTree, Validation) {
    const std::vector<Real> pts = {0.0, 1.0, 2.0, 3.0};
    const std::vector<Size> dims2 = {0, 1};
    const std::vector<Size> empty_dims;
    const std::vector<Size> bad_dim = {0, 4};  // stride=2, 列 4 越界
    EXPECT_THROW(den::KdTree(pts, 0, 2, dims2), std::invalid_argument);
    EXPECT_THROW(den::KdTree(pts, 2, 2, empty_dims), std::invalid_argument);
    EXPECT_THROW(den::KdTree(pts, 2, 2, bad_dim), std::invalid_argument);
    EXPECT_THROW(den::KdTree(pts, 3, 2, dims2), std::invalid_argument);  // 缓冲区不足
    const std::vector<Real> nan_pts = {0.0, std::numeric_limits<Real>::quiet_NaN()};
    EXPECT_THROW(den::KdTree(nan_pts, 1, 2, dims2), std::invalid_argument);

    const den::KdTree tree(pts, 2, 2, dims2);
    EXPECT_THROW(tree.kth_distance_including_self(2, 1), std::invalid_argument);
    EXPECT_THROW(tree.kth_distance_including_self(0, 0), std::invalid_argument);
    EXPECT_THROW(tree.kth_distance_including_self(0, 3), std::invalid_argument);
    EXPECT_THROW(tree.count_within(0, -0.1), std::invalid_argument);
    EXPECT_THROW(tree.count_within(0, std::numeric_limits<Real>::quiet_NaN()),
                 std::invalid_argument);

    const den::KdTree unbuilt;
    EXPECT_THROW(unbuilt.kth_distance_including_self(0, 1), std::logic_error);
    EXPECT_THROW(unbuilt.count_within(0, 0.0), std::logic_error);
}
