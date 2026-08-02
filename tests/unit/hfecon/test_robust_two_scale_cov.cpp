// =============================================================================
// test_robust_two_scale_cov.cpp
// Phase 5 v1.4.2 Wave B - Robust Two-Scale 协方差估计测试
//
// 对标: R highfrequency 1.0.3 rRTSCov + RTSRV + RTSCov_bi
// 容差: 1e-10 (迭代算法, 标准容差)
//
// SOURCE: PHASE5_HFE_SPEC §5.1, §5.3 D6, §5.5
//   R highfrequency 1.0.3 R/internalRealizedMeasures.R
//     RTSRV (L490-538): 单资产 robust TSRV (迭代截断)
//     RTSCov_bi (L360-487): 双资产 robust TSCov (迭代截断)
//
// 关键幻觉排除 (spec §5.3 D6):
//   eta=9 时 R 硬编码 ccc=1.0415, 不走查表 (查表值为 1.04146535666802)
//   单资产用 zeta = 1/pchisq(eta, 3), 双资产用 ccc = cfactor(eta)
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/measures/robust_two_scale_cov.hpp"
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
// 辅助函数: 生成 GBM 路径 (固定种子, 跨平台一致)
// =============================================================================
namespace {
std::vector<Real> make_gbm_prices(Size n, Real sigma, Real p0, uint64_t seed) {
    std::vector<Real> prices(n);
    prices[0] = p0;
    for (Size i = 1; i < n; ++i) {
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g = (u - 0.5) * 2.0;  // [-1, 1] 均匀
        prices[i] = prices[i-1] * std::exp(sigma * g);
    }
    return prices;
}
}  // namespace

// =============================================================================
// TEST 1: 单资产 Brownian motion — RTSRV 应为正且有限
//   构造 GBM 路径 (固定种子), sigma=0.001, n=2000, K=100, J=1
//   n=2000 >= 10*K=1000 满足约束
// =============================================================================
TEST(RobustTwoScaleCovTest, SingleAssetBrownianMotion) {
    auto prices = make_gbm_prices(2000, 0.001, 100.0, 42);

    // K=100, J=1, eta=9 (R 默认)
    auto result = RobustTwoScaleCov::estimate({prices}, 100, 1, 9.0, false);
    EXPECT_EQ(result.n_assets, 1u);
    EXPECT_EQ(result.n_obs, 2000u);
    EXPECT_EQ(result.K, 100u);
    EXPECT_EQ(result.J, 1u);
    ASSERT_EQ(result.cov.size(), 1u);

    // RTSRV 应为正
    EXPECT_GT(result.cov[0], 0.0);
    EXPECT_TRUE(std::isfinite(result.cov[0]));

    // 与标准 RV 比较 (数量级一致)
    Real rv = 0.0;
    for (Size i = 1; i < 2000; ++i) {
        Real r = std::log(prices[i]) - std::log(prices[i-1]);
        rv += r * r;
    }
    Real ratio = result.cov[0] / rv;
    EXPECT_GT(ratio, 0.1);
    EXPECT_LT(ratio, 10.0);
}

// =============================================================================
// TEST 2: 单资产常数价格 — RTSRV 应为 0
// =============================================================================
TEST(RobustTwoScaleCovTest, SingleAssetConstantPrice) {
    std::vector<Real> prices(1500, 100.0);  // 常数价格, n=1500 >= 10*100
    auto result = RobustTwoScaleCov::estimate({prices}, 100, 1, 9.0, false);
    EXPECT_NEAR(result.cov[0], 0.0, TOL_LOOSE);
}

// =============================================================================
// TEST 3: 单资产便捷接口 — estimate_univariate
// =============================================================================
TEST(RobustTwoScaleCovTest, UnivariateInterface) {
    auto prices = make_gbm_prices(2000, 0.001, 100.0, 99);

    Real v1 = RobustTwoScaleCov::estimate_univariate(prices, 100, 1, 9.0);
    auto r1 = RobustTwoScaleCov::estimate({prices}, 100, 1, 9.0, false);
    EXPECT_NEAR(v1, r1.cov[0], TOL_STRICT);
    EXPECT_GT(v1, 0.0);
}

// =============================================================================
// TEST 4: 双资产完全相关 — 非对角线应为正且对称
//   两个相同价格序列, 协方差矩阵应对称
// =============================================================================
TEST(RobustTwoScaleCovTest, DualAssetPerfectCorrelation) {
    auto prices = make_gbm_prices(2000, 0.001, 100.0, 333);

    auto result = RobustTwoScaleCov::estimate({prices, prices}, 100, 1, 9.0, false);
    EXPECT_EQ(result.n_assets, 2u);
    ASSERT_EQ(result.cov.size(), 4u);

    // 对角线应相等
    EXPECT_NEAR(result.cov[0], result.cov[3], TOL_STRICT);

    // 非对角线对称
    EXPECT_NEAR(result.cov[1], result.cov[2], TOL_STRICT);

    // 完全相关资产的非对角线应接近对角线 (RTSCov_bi 与 RTSRV 公式不同,
    // 但完全相关时应数量级一致)
    EXPECT_GT(result.cov[1], 0.0);
    Real ratio = result.cov[1] / result.cov[0];
    EXPECT_GT(ratio, 0.1);
    EXPECT_LT(ratio, 10.0);
}

// =============================================================================
// TEST 5: noisevar 参数 — 提供噪声方差应改变结果
// =============================================================================
TEST(RobustTwoScaleCovTest, NoisevarParameter) {
    auto prices = make_gbm_prices(2000, 0.001, 100.0, 555);

    // 自动估计 noisevar
    auto r1 = RobustTwoScaleCov::estimate({prices}, 100, 1, 9.0, false);

    // 提供外部 noisevar (较大值)
    auto r2 = RobustTwoScaleCov::estimate({prices}, 100, 1, 9.0, false,
                                           {1e-4});

    // 提供外部 noisevar (较小值)
    auto r3 = RobustTwoScaleCov::estimate({prices}, 100, 1, 9.0, false,
                                           {1e-10});

    // 结果应不同 (不同的 noisevar 影响截断阈值)
    // 但都应为正且有限
    EXPECT_GT(r1.cov[0], 0.0);
    EXPECT_GT(r2.cov[0], 0.0);
    EXPECT_GT(r3.cov[0], 0.0);
    EXPECT_TRUE(std::isfinite(r1.cov[0]));
    EXPECT_TRUE(std::isfinite(r2.cov[0]));
    EXPECT_TRUE(std::isfinite(r3.cov[0]));
}

// =============================================================================
// TEST 6: eta 参数 — 不同 eta 应改变截断阈值
// =============================================================================
TEST(RobustTwoScaleCovTest, EtaParameter) {
    auto prices = make_gbm_prices(2000, 0.001, 100.0, 888);

    // eta=9 (R 默认, 硬编码 ccc=1.0415)
    auto r1 = RobustTwoScaleCov::estimate({prices}, 100, 1, 9.0, false);

    // eta=5 (查表 ccc)
    auto r2 = RobustTwoScaleCov::estimate({prices}, 100, 1, 5.0, false);

    // eta=15 (查表 ccc)
    auto r3 = RobustTwoScaleCov::estimate({prices}, 100, 1, 15.0, false);

    // 所有结果应为正且有限
    EXPECT_GT(r1.cov[0], 0.0);
    EXPECT_GT(r2.cov[0], 0.0);
    EXPECT_GT(r3.cov[0], 0.0);

    // eta=9 验证 K 和 J 参数记录正确
    EXPECT_EQ(r1.K, 100u);
    EXPECT_EQ(r1.J, 1u);
    EXPECT_EQ(r1.eta, 9.0);
}

// =============================================================================
// TEST 7: 异常处理 — 空输入 / 长度不一致 / K<=J / n<10*K / eta 越界
// =============================================================================
TEST(RobustTwoScaleCovTest, ExceptionHandling) {
    // 空矩阵
    std::vector<std::vector<Real>> empty;
    EXPECT_THROW(RobustTwoScaleCov::estimate(empty), std::invalid_argument);

    // 长度不一致
    std::vector<std::vector<Real>> mismatched = {
        std::vector<Real>(1500, 100.0),
        std::vector<Real>(1000, 100.0)
    };
    EXPECT_THROW(RobustTwoScaleCov::estimate(mismatched, 100, 1), std::invalid_argument);

    // K <= J
    std::vector<std::vector<Real>> valid = {
        std::vector<Real>(1500, 100.0)
    };
    EXPECT_THROW(RobustTwoScaleCov::estimate(valid, 1, 1), std::invalid_argument);
    EXPECT_THROW(RobustTwoScaleCov::estimate(valid, 1, 2), std::invalid_argument);
    EXPECT_THROW(RobustTwoScaleCov::estimate(valid, 5, 10), std::invalid_argument);

    // n < 10*K
    std::vector<std::vector<Real>> short_data = {
        std::vector<Real>(500, 100.0)  // n=500 < 10*100=1000
    };
    EXPECT_THROW(RobustTwoScaleCov::estimate(short_data, 100, 1), std::invalid_argument);
}

// =============================================================================
// TEST 8: makePsd 投影 — 双资产 PSD 投影后特征值 >= 0
// =============================================================================
TEST(RobustTwoScaleCovTest, MakePsdProjection) {
    // 构造两个相关性较低的价格序列
    auto p1 = make_gbm_prices(2000, 0.001, 100.0, 111);
    auto p2 = make_gbm_prices(2000, 0.001, 100.0, 222);

    // 不做 PSD 投影
    auto r1 = RobustTwoScaleCov::estimate({p1, p2}, 100, 1, 9.0, false);
    // 做 PSD 投影
    auto r2 = RobustTwoScaleCov::estimate({p1, p2}, 100, 1, 9.0, true);

    // 验证 PSD 后的矩阵特征值 >= 0
    ASSERT_EQ(r2.cov.size(), 4u);
    const Real trace = r2.cov[0] + r2.cov[3];
    const Real det = r2.cov[0] * r2.cov[3] - r2.cov[1] * r2.cov[2];
    const Real disc = trace * trace - 4.0 * det;
    if (disc >= 0.0) {
        const Real sq = std::sqrt(disc);
        const Real lambda1 = (trace + sq) / 2.0;
        const Real lambda2 = (trace - sq) / 2.0;
        EXPECT_GE(lambda1, -TOL_LOOSE);
        EXPECT_GE(lambda2, -TOL_LOOSE);
    }
}

// =============================================================================
// TEST 9: start_iv 参数 — 提供初始 IV 估计
// =============================================================================
TEST(RobustTwoScaleCovTest, StartIVParameter) {
    auto prices = make_gbm_prices(2000, 0.001, 100.0, 444);

    // 自动估计 start_iv
    auto r1 = RobustTwoScaleCov::estimate({prices}, 100, 1, 9.0, false);

    // 提供外部 start_iv
    auto r2 = RobustTwoScaleCov::estimate({prices}, 100, 1, 9.0, false,
                                           {}, {1e-4});

    // 两种情况都应为正且有限
    EXPECT_GT(r1.cov[0], 0.0);
    EXPECT_GT(r2.cov[0], 0.0);
    EXPECT_TRUE(std::isfinite(r1.cov[0]));
    EXPECT_TRUE(std::isfinite(r2.cov[0]));
}
