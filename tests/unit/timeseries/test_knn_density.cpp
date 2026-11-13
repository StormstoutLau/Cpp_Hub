// =============================================================================
// test_knn_density.cpp - Phase 9 M2: 暴力 Chebyshev kNN 密度 (TDD 前置)
//
// 状态: M2 已实现 (include/cpphub/timeseries/density/knn_density.hpp),
//       已接入 tests/CMakeLists.txt.
//
// 定义: Loftsgaarden-Quesenberry 1D 密度
//   pdf(x) = k / (n * 2 * r_k(x)), r_k = 到第 k 近邻的绝对距离
//   边界外按 r_k 截断, 不做归一化 (spec §4.1).
// =============================================================================

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include "cpphub/timeseries/density/knn_density.hpp"

namespace de = cpphub::v1::timeseries::density;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

const std::vector<Real> kSample = {0.0, 1.0, 2.0, 3.0, 4.0};

de::KnnConfig knn_config(Size k) {
    de::KnnConfig cfg;
    cfg.k = k;
    return cfg;
}

}  // namespace

TEST(KnnDensity, GoldenValuesK2) {
    de::KnnDensity1D knn(std::span<const Real>(kSample), knn_config(2));
    EXPECT_EQ(knn.k(), 2u);
    EXPECT_NEAR(knn.pdf(1.5), 0.4, 1e-15);
    EXPECT_NEAR(knn.pdf(0.0), 0.2, 1e-15);
    EXPECT_NEAR(knn.pdf(2.0), 0.2, 1e-15);
    EXPECT_NEAR(knn.pdf(10.0), 0.02857142857142857, 1e-15);
}

TEST(KnnDensity, GoldenValuesK4) {
    de::KnnDensity1D knn(std::span<const Real>(kSample), knn_config(4));
    EXPECT_NEAR(knn.pdf(1.5), 0.26666666666666666, 1e-15);
    EXPECT_NEAR(knn.pdf(0.0), 0.13333333333333333, 1e-15);
    EXPECT_NEAR(knn.pdf(2.0), 0.2, 1e-15);
    EXPECT_NEAR(knn.pdf(10.0), 0.044444444444444446, 1e-15);
}

TEST(KnnDensity, GridMatchesScalar) {
    de::KnnDensity1D knn(std::span<const Real>(kSample), knn_config(2));
    const std::vector<Real> grid = {-2.0, 0.0, 1.5, 2.0, 6.0};
    const auto out = knn.pdf(std::span<const Real>(grid));
    ASSERT_EQ(out.size(), grid.size());
    for (Size i = 0; i < grid.size(); ++i) {
        EXPECT_NEAR(out[i], knn.pdf(grid[i]), 1e-15);
    }
}

TEST(KnnDensity, OutsideRangeIsFinitePositive) {
    de::KnnDensity1D knn(std::span<const Real>(kSample), knn_config(2));
    for (Real x : {-100.0, -10.0, 10.0, 100.0}) {
        const Real v = knn.pdf(x);
        EXPECT_TRUE(std::isfinite(v));
        EXPECT_GT(v, 0.0);
    }
}

TEST(KnnDensity, KEqualsNMinusOne) {
    de::KnnDensity1D knn(std::span<const Real>(kSample), knn_config(4));
    EXPECT_NEAR(knn.pdf(2.0), 0.2, 1e-15);
}

TEST(KnnDensity, Validation) {
    std::vector<Real> empty;
    std::vector<Real> nan_sample = {0.0, std::numeric_limits<Real>::quiet_NaN()};

    EXPECT_THROW(de::KnnDensity1D(std::span<const Real>(empty), knn_config(1)),
                 std::invalid_argument);
    EXPECT_THROW(de::KnnDensity1D(std::span<const Real>(nan_sample), knn_config(1)),
                 std::invalid_argument);
    EXPECT_THROW(de::KnnDensity1D(std::span<const Real>(kSample), knn_config(0)),
                 std::invalid_argument);
    EXPECT_THROW(de::KnnDensity1D(std::span<const Real>(kSample), knn_config(5)),
                 std::invalid_argument);
}
