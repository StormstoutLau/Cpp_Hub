// =============================================================================
// test_preaveraged_cov.cpp
// Phase 5 v1.4.2 Wave B - Pre-averaged (Subsampled) 协方差估计测试
//
// 对标: R highfrequency 1.0.3 rAVGCov
// 容差: 1e-12 (合成数据, 严格)
//
// SOURCE: PHASE5_HFE_SPEC §5.1, §5.3 D2, §5.5
//   R highfrequency 1.0.3 R/realizedMeasures.R rAVGCov (L790-951)
//
// 关键幻觉排除 (spec §5.3 D2):
//   单资产加 (m+1)/m 系数 (m = alignPeriod/k = scalingFraction)
//   多资产不加该校正系数
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/measures/preaveraged_cov.hpp"
#include "cpphub/hfecon/measures/cov_utils.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace cpphub::v1::hfecon;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
constexpr Real TOL_STRICT = 1e-12;
constexpr Real TOL_LOOSE  = 1e-10;
}  // namespace

// =============================================================================
// TEST 1: 单资产常数价格 — rAVGCov 应为 0
// =============================================================================
TEST(PreaveragedCovTest, SingleAssetConstantPrice) {
    std::vector<Real> prices(100, 100.0);
    auto result = PreaveragedCov::estimate({prices}, 5, 1, false, false);
    EXPECT_EQ(result.n_assets, 1u);
    EXPECT_EQ(result.scaling_fraction, 5u);
    ASSERT_EQ(result.cov.size(), 1u);
    EXPECT_NEAR(result.cov[0], 0.0, TOL_STRICT);
}

// =============================================================================
// TEST 2: 单资产 Brownian motion — rAVGCov 应接近 RV
//   构造 GBM 路径 (固定种子), sigma=0.001, n=500
//   s=1 时 rAVGCov 退化为标准 RV
//   s>1 时 rAVGCov 应与 RV 数量级一致
// =============================================================================
TEST(PreaveragedCovTest, SingleAssetBrownianMotion) {
    // 简单 LCG 伪随机数 (固定种子, 跨平台一致)
    std::vector<Real> prices(500);
    prices[0] = 100.0;
    uint64_t seed = 12345;
    for (Size i = 1; i < 500; ++i) {
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g = (u - 0.5) * 2.0;  // [-1, 1] 均匀
        prices[i] = prices[i-1] * std::exp(0.001 * g);  // sigma=0.001
    }

    // s=1: 退化为标准 RV
    auto r1 = PreaveragedCov::estimate({prices}, 1, 1, false, false);
    EXPECT_EQ(r1.scaling_fraction, 1u);

    // 计算标准 RV (单 tick log returns)
    Real rv = 0.0;
    for (Size i = 1; i < 500; ++i) {
        Real r = std::log(prices[i]) - std::log(prices[i-1]);
        rv += r * r;
    }
    EXPECT_NEAR(r1.cov[0], rv, TOL_STRICT);  // s=1 应严格等于 RV

    // s=5: rAVGCov 应与 RV 数量级一致 (含 (m+1)/m 校正)
    auto r5 = PreaveragedCov::estimate({prices}, 5, 1, false, false);
    EXPECT_GT(r5.cov[0], 0.0);
    Real ratio = r5.cov[0] / rv;
    EXPECT_GT(ratio, 0.5);
    EXPECT_LT(ratio, 2.0);
}

// =============================================================================
// TEST 3: 双资产完全相关 — 非对角线应等于对角线
//   两个相同价格序列, 协方差矩阵应对称且非对角线为正
// =============================================================================
TEST(PreaveragedCovTest, DualAssetPerfectCorrelation) {
    std::vector<Real> prices(200);
    prices[0] = 100.0;
    uint64_t seed = 77;
    for (Size i = 1; i < 200; ++i) {
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g = (u - 0.5) * 2.0;
        prices[i] = prices[i-1] * std::exp(0.001 * g);
    }

    // 两个相同资产
    auto result = PreaveragedCov::estimate({prices, prices}, 5, 1, false, false);
    EXPECT_EQ(result.n_assets, 2u);
    ASSERT_EQ(result.cov.size(), 4u);

    // 对角线应相等
    EXPECT_NEAR(result.cov[0], result.cov[3], TOL_STRICT);

    // 完全相关: 非对角线应等于对角线 (相同资产)
    // 注意: 多资产不应用 (m+1)/m 校正, 但对角线应用了
    // 完全相关时 cov[1] = cov[0] (因为双资产 cross-product 退化)
    EXPECT_GT(result.cov[1], 0.0);
    EXPECT_NEAR(result.cov[1], result.cov[2], TOL_STRICT);  // 对称性

    // 验证非对角线与对角线的关系
    // 双资产 cross-product = sum(r1_i * r2_i) / s (无校正)
    // 单资产 RV with correction = (m+1)/m * sum(r_i^2) / s
    // 完全相关时 r1 = r2, 所以 cross-product = sum(r_i^2) / s
    // 而 diagonal = (m+1)/m * sum(r_i^2) / s
    // 所以 cov[1] / cov[0] ≈ m / (m+1), 其中 m = (N-s)/s
    // 注意: 由于 bucket 边界处理, 实际 ratio 与理论值有小差异, 用宽松容差
    Real ratio = result.cov[1] / result.cov[0];
    Real N_ret = static_cast<Real>(prices.size() - 1);  // returns 数量
    Real s_real = 5.0;
    Real m = (N_ret - s_real) / s_real;  // m = (N-s)/s
    Real expected_ratio = m / (m + 1.0);
    // 宽松容差: bucket 边界效应导致 ~0.5% 偏差
    EXPECT_NEAR(ratio, expected_ratio, 0.01);
}

// =============================================================================
// TEST 4: k 参数 — k>1 时 scalingFraction = alignPeriod/k
// =============================================================================
TEST(PreaveragedCovTest, KParameter) {
    std::vector<Real> prices(200);
    prices[0] = 100.0;
    for (Size i = 1; i < 200; ++i) {
        prices[i] = prices[i-1] * 1.001;
    }

    // alignPeriod=6, k=2: s = 3
    auto r1 = PreaveragedCov::estimate({prices}, 6, 2, false, false);
    EXPECT_EQ(r1.scaling_fraction, 3u);

    // alignPeriod=6, k=3: s = 2
    auto r2 = PreaveragedCov::estimate({prices}, 6, 3, false, false);
    EXPECT_EQ(r2.scaling_fraction, 2u);

    // alignPeriod=4, k=2: s = 2 (与 alignPeriod=2, k=1 等价)
    auto r3 = PreaveragedCov::estimate({prices}, 4, 2, false, false);
    auto r4 = PreaveragedCov::estimate({prices}, 2, 1, false, false);
    EXPECT_EQ(r3.scaling_fraction, r4.scaling_fraction);
    // 由于 (m+1)/m 校正因子不同 (m=2 vs m=2), 应严格相等
    EXPECT_NEAR(r3.cov[0], r4.cov[0], TOL_STRICT);
}

// =============================================================================
// TEST 5: 异常处理 — 空输入 / 长度不一致 / k=0 / alignPeriod%k!=0
// =============================================================================
TEST(PreaveragedCovTest, ExceptionHandling) {
    // 空矩阵
    std::vector<std::vector<Real>> empty;
    EXPECT_THROW(PreaveragedCov::estimate(empty), std::invalid_argument);

    // 长度不一致
    std::vector<std::vector<Real>> mismatched = {
        {100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0},
        {100.0, 101.0, 102.0}
    };
    EXPECT_THROW(PreaveragedCov::estimate(mismatched, 1, 1), std::invalid_argument);

    // k=0
    std::vector<std::vector<Real>> valid = {
        {100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0}
    };
    EXPECT_THROW(PreaveragedCov::estimate(valid, 5, 0), std::invalid_argument);

    // alignPeriod % k != 0
    EXPECT_THROW(PreaveragedCov::estimate(valid, 5, 2), std::invalid_argument);

    // N < s+1
    EXPECT_THROW(PreaveragedCov::estimate(valid, 10, 1), std::invalid_argument);
}

// =============================================================================
// TEST 6: 相关矩阵转换 — cor=true 时对角线应为 1
// =============================================================================
TEST(PreaveragedCovTest, CorrelationMatrix) {
    std::vector<Real> p1(200), p2(200);
    p1[0] = 100.0;
    p2[0] = 100.0;
    uint64_t seed1 = 11, seed2 = 22;
    for (Size i = 1; i < 200; ++i) {
        seed1 = (6364136223846793005ULL * seed1 + 1442695040888963407ULL);
        seed2 = (6364136223846793005ULL * seed2 + 1442695040888963407ULL);
        const Real u1 = static_cast<Real>((seed1 >> 11) & 0xFFFFFF) / 16777216.0;
        const Real u2 = static_cast<Real>((seed2 >> 11) & 0xFFFFFF) / 16777216.0;
        p1[i] = p1[i-1] * std::exp(0.001 * (u1 - 0.5) * 2.0);
        p2[i] = p2[i-1] * std::exp(0.001 * (u2 - 0.5) * 2.0);
    }

    auto result = PreaveragedCov::estimate({p1, p2}, 5, 1, true, false);
    EXPECT_EQ(result.n_assets, 2u);
    ASSERT_EQ(result.cov.size(), 4u);

    // 对角线应为 1
    EXPECT_NEAR(result.cov[0], 1.0, TOL_LOOSE);
    EXPECT_NEAR(result.cov[3], 1.0, TOL_LOOSE);

    // 非对角线应在 [-1, 1]
    EXPECT_GE(result.cov[1], -1.0 - TOL_LOOSE);
    EXPECT_LE(result.cov[1], 1.0 + TOL_LOOSE);
    EXPECT_NEAR(result.cov[1], result.cov[2], TOL_STRICT);  // 对称性
}

// =============================================================================
// TEST 7: 单资产便捷接口 — estimate_univariate
// =============================================================================
TEST(PreaveragedCovTest, UnivariateInterface) {
    std::vector<Real> prices(100);
    prices[0] = 100.0;
    for (Size i = 1; i < 100; ++i) {
        prices[i] = prices[i-1] * 1.0005;
    }

    Real v1 = PreaveragedCov::estimate_univariate(prices, 5, 1);
    auto r1 = PreaveragedCov::estimate({prices}, 5, 1, false, false);
    EXPECT_NEAR(v1, r1.cov[0], TOL_STRICT);

    // k=0 异常
    EXPECT_THROW(PreaveragedCov::estimate_univariate(prices, 5, 0),
                 std::invalid_argument);
    // alignPeriod % k != 0 异常
    EXPECT_THROW(PreaveragedCov::estimate_univariate(prices, 5, 2),
                 std::invalid_argument);
}
