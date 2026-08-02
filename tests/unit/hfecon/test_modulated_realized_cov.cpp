// =============================================================================
// test_modulated_realized_cov.cpp
// Phase 5 v1.4.2 Wave B - Modulated Realized Covariance 测试
//
// 对标: R highfrequency 1.0.3 rMRCov + crv + preavbi + hatreturn + gfunction
// 容差: 1e-12 (合成数据, 严格)
//
// SOURCE: PHASE5_HFE_SPEC §5.1, §5.3 D3, §5.5
//   R highfrequency 1.0.3 R/internalPreaveringEstimators.R
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/measures/modulated_realized_cov.hpp"
#include "cpphub/hfecon/measures/cov_utils.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace cpphub::v1::hfecon;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
constexpr Real TOL_STRICT = 1e-12;
constexpr Real TOL_LOOSE = 1e-8;
}  // namespace

// =============================================================================
// TEST 1: 单资产 — 常数价格序列 (RV=0)
//   prices = [100, 100, ..., 100]
//   ret = [0, 0, ..., 0]
//   hatreturn = [0, 0, ..., 0]
//   crv = 0 - 0 = 0
// =============================================================================
TEST(ModulatedRealizedCovTest, SingleAssetConstant) {
    std::vector<Real> prices(100, 100.0);

    auto result = ModulatedRealizedCov::estimate({prices});
    EXPECT_EQ(result.n_assets, 1u);
    EXPECT_EQ(result.n_obs, 100u);
    ASSERT_EQ(result.cov.size(), 1u);
    EXPECT_NEAR(result.cov[0], 0.0, TOL_STRICT);
}

// =============================================================================
// TEST 2: 单资产 Brownian motion — crv 应接近 RV (无噪声)
//   BM 无微结构噪声, crv 应与 RV 数量级一致
// =============================================================================
TEST(ModulatedRealizedCovTest, SingleAssetBrownian) {
    // 固定种子的简单 LCG (跨平台一致)
    std::vector<Real> prices(1000);
    prices[0] = 100.0;
    uint64_t seed = 42;
    for (Size i = 1; i < 1000; ++i) {
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g = (u - 0.5) * 2.0;  // [-1, 1] 均匀
        prices[i] = prices[i - 1] * std::exp(0.001 * g);  // sigma=0.001
    }

    auto result = ModulatedRealizedCov::estimate({prices});
    EXPECT_EQ(result.n_assets, 1u);
    EXPECT_EQ(result.n_obs, 1000u);
    ASSERT_EQ(result.cov.size(), 1u);

    // 计算 RV 对比
    Real rv = 0.0;
    for (Size i = 1; i < 1000; ++i) {
        const Real r = std::log(prices[i]) - std::log(prices[i - 1]);
        rv += r * r;
    }

    // 无噪声时, crv 应与 RV 数量级一致 (容差: 因子 0.1 ~ 10)
    EXPECT_GT(result.cov[0], 0.0);
    EXPECT_GT(rv, 0.0);
    const Real ratio = result.cov[0] / rv;
    EXPECT_GT(ratio, 0.1);
    EXPECT_LT(ratio, 10.0);
}

// =============================================================================
// TEST 3: 双资产完全相关 — 协方差 = 方差
//   prices1 = prices2 (完全同步, 完全相关)
//   preavbi(prices1, prices2) ≈ crv(prices1) ?
//
//   注意: crv 用 1/(sqrt(N)*theta*psi2kn) * sum(r1^2)
//         preavbi 用 N/(N-kn+2) * 1/(psi2*kn) * sum(r1*r2)
//         公式不同, 但数量级一致
// =============================================================================
TEST(ModulatedRealizedCovTest, DualAssetPerfectlyCorrelated) {
    std::vector<Real> prices(400);
    prices[0] = 100.0;
    uint64_t seed = 7;
    for (Size i = 1; i < 400; ++i) {
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g = (u - 0.5) * 2.0;
        prices[i] = prices[i - 1] * std::exp(0.001 * g);
    }

    // 双资产完全相关
    auto result = ModulatedRealizedCov::estimate({prices, prices});
    EXPECT_EQ(result.n_assets, 2u);
    EXPECT_EQ(result.n_obs, 400u);
    ASSERT_EQ(result.cov.size(), 4u);

    // 完全相关时, cov[off-diag] 应与 cov[diag] 同号且接近 (但公式不同)
    EXPECT_GT(result.cov[0], 0.0);  // var1
    EXPECT_GT(result.cov[3], 0.0);  // var2
    // var1 == var2 (因为 prices1 == prices2)
    EXPECT_NEAR(result.cov[0], result.cov[3], TOL_LOOSE);
    // cov12 同号
    EXPECT_GT(result.cov[1] * result.cov[0], 0.0);
    // 矩阵对称
    EXPECT_NEAR(result.cov[1], result.cov[2], TOL_STRICT);
}

// =============================================================================
// TEST 4: theta 参数敏感性 — 不同 theta 应产生不同结果
// =============================================================================
TEST(ModulatedRealizedCovTest, ThetaSensitivity) {
    std::vector<Real> prices(400);
    prices[0] = 100.0;
    uint64_t seed = 13;
    for (Size i = 1; i < 400; ++i) {
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g = (u - 0.5) * 2.0;
        prices[i] = prices[i - 1] * std::exp(0.001 * g);
    }

    auto r1 = ModulatedRealizedCov::estimate({prices}, 0.5);
    auto r2 = ModulatedRealizedCov::estimate({prices}, 0.8);
    auto r3 = ModulatedRealizedCov::estimate({prices}, 1.0);

    // 不同 theta 应产生不同的 crv 值
    EXPECT_NE(r1.cov[0], r2.cov[0]);
    EXPECT_NE(r2.cov[0], r3.cov[0]);
    EXPECT_NE(r1.cov[0], r3.cov[0]);

    // kn 应随 theta 增大而增大
    EXPECT_EQ(r1.kn, static_cast<Size>(std::floor(0.5 * std::sqrt(400.0))));
    EXPECT_EQ(r2.kn, static_cast<Size>(std::floor(0.8 * std::sqrt(400.0))));
    EXPECT_EQ(r3.kn, static_cast<Size>(std::floor(1.0 * std::sqrt(400.0))));
}

// =============================================================================
// TEST 5: 双资产 pairwise 模式 — 对角线用 crv, 非对角线用 preavbi
// =============================================================================
TEST(ModulatedRealizedCovTest, PairwiseMode) {
    std::vector<Real> prices1(400), prices2(400);
    prices1[0] = 100.0;
    prices2[0] = 100.0;
    uint64_t seed = 17;
    for (Size i = 1; i < 400; ++i) {
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g1 = (u - 0.5) * 2.0;
        prices1[i] = prices1[i - 1] * std::exp(0.001 * g1);
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u2 = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g2 = (u2 - 0.5) * 2.0;
        prices2[i] = prices2[i - 1] * std::exp(0.001 * g2);
    }

    auto r_matrix = ModulatedRealizedCov::estimate({prices1, prices2}, 0.8, false, false);
    auto r_pairwise = ModulatedRealizedCov::estimate({prices1, prices2}, 0.8, true, false);

    EXPECT_EQ(r_pairwise.n_assets, 2u);
    EXPECT_EQ(r_pairwise.n_obs, 400u);

    // Pairwise 的对角线应与单资产 crv 一致
    auto r_single1 = ModulatedRealizedCov::estimate({prices1}, 0.8);
    auto r_single2 = ModulatedRealizedCov::estimate({prices2}, 0.8);
    EXPECT_NEAR(r_pairwise.cov[0], r_single1.cov[0], TOL_STRICT);
    EXPECT_NEAR(r_pairwise.cov[3], r_single2.cov[0], TOL_STRICT);

    // 矩阵模式和 pairwise 模式的对角线可能不同 (矩阵版用 N/(N-kn+2) 系数, crv 用 1/sqrt(N) 系数)
    // 但应同号
    EXPECT_GT(r_matrix.cov[0] * r_pairwise.cov[0], 0.0);

    // 矩阵对称性
    EXPECT_NEAR(r_pairwise.cov[1], r_pairwise.cov[2], TOL_STRICT);
}

// =============================================================================
// TEST 6: 异常处理 — 空输入, 长度不一致, kn 过小
// =============================================================================
TEST(ModulatedRealizedCovTest, ExceptionHandling) {
    // 空输入
    EXPECT_THROW(ModulatedRealizedCov::estimate({}), std::invalid_argument);

    // 长度不一致
    std::vector<Real> p1(100, 100.0);
    std::vector<Real> p2(200, 100.0);
    EXPECT_THROW(ModulatedRealizedCov::estimate({p1, p2}), std::invalid_argument);

    // kn 过小 (N=4, theta=0.8 -> kn = floor(0.8*2) = 1 < 2)
    std::vector<Real> short_prices(4, 100.0);
    EXPECT_THROW(ModulatedRealizedCov::estimate({short_prices}, 0.8),
                 std::invalid_argument);
}

// =============================================================================
// TEST 7: makePsd 投影 — 非对角线噪声校正可能导致非 PSD, makePsd 应恢复 PSD
// =============================================================================
TEST(ModulatedRealizedCovTest, MakePsdProjection) {
    std::vector<Real> prices1(400), prices2(400);
    prices1[0] = 100.0;
    prices2[0] = 100.0;
    uint64_t seed = 23;
    for (Size i = 1; i < 400; ++i) {
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g1 = (u - 0.5) * 2.0;
        prices1[i] = prices1[i - 1] * std::exp(0.001 * g1);
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u2 = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g2 = (u2 - 0.5) * 2.0;
        prices2[i] = prices2[i - 1] * std::exp(0.001 * g2);
    }

    // 不做 PSD 投影
    auto r1 = ModulatedRealizedCov::estimate({prices1, prices2}, 0.8, false, false, false);
    // 做 PSD 投影
    auto r2 = ModulatedRealizedCov::estimate({prices1, prices2}, 0.8, false, true, false);

    EXPECT_EQ(r1.n_assets, 2u);
    EXPECT_EQ(r2.n_assets, 2u);

    // r2 应为 PSD: 所有特征值 >= -tolerance
    // 2x2 矩阵 [[a, b], [b, d]] 的特征值 = (a+d)/2 ± sqrt(((a-d)/2)^2 + b^2)
    // 最小特征值 = (a+d)/2 - sqrt(((a-d)/2)^2 + b^2) >= 0
    const Real a = r2.cov[0];
    const Real d = r2.cov[3];
    const Real b = r2.cov[1];
    const Real trace = a + d;
    const Real disc = std::sqrt((a - d) * (a - d) / 4.0 + b * b);
    const Real min_eig = trace / 2.0 - disc;
    EXPECT_GE(min_eig, -TOL_LOOSE);

    // 对称性
    EXPECT_NEAR(r2.cov[1], r2.cov[2], TOL_STRICT);
}
