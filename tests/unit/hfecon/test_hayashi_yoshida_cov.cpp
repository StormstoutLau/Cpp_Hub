// =============================================================================
// test_hayashi_yoshida_cov.cpp
// Phase 5 v1.4.2 Wave A - Hayashi-Yoshida 协方差估计测试
//
// 对标: R highfrequency 1.0.3 rHYCov + pcovcc + makePsd
// 容差: 1e-12 (合成数据, 严格)
//
// SOURCE: PHASE5_HFE_SPEC §5.1, §5.3 D4/D5, §5.5
//   R highfrequency 1.0.3 rHYCov (v142_source_dump.txt L895-1035)
//   R highfrequency 1.0.3 src/realizedMeasures.cpp pcovcc
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/measures/hayashi_yoshida_cov.hpp"
#include "cpphub/hfecon/measures/cov_utils.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace cpphub::v1::hfecon;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
constexpr Real TOL_STRICT = 1e-12;
}  // namespace

// =============================================================================
// TEST 1: 单资产常数序列 — RV = 0, cov = [[0]]
// =============================================================================
TEST(HayashiYoshidaCovTest, SingleAssetConstant) {
    // 常数收益率 -> RV = 0
    std::vector<std::vector<Real>> returns = {
        {0.0, 0.0, 0.0, 0.0, 0.0}
    };

    auto result = HayashiYoshidaCov::estimate(returns);
    EXPECT_EQ(result.n_assets, 1u);
    EXPECT_EQ(result.n_obs, 5u);
    ASSERT_EQ(result.cov.size(), 1u);
    EXPECT_NEAR(result.cov[0], 0.0, TOL_STRICT);
}

// =============================================================================
// TEST 2: 双资产 period=1 (同步) — 退化为标准 realized covariance
//   r1 = [0.01, 0.02, -0.01, 0.03, -0.02]
//   r2 = [0.02, 0.01, -0.02, 0.01, -0.01]
//   var1 = 0.0001+0.0004+0.0001+0.0009+0.0004 = 0.0019
//   var2 = 0.0004+0.0001+0.0004+0.0001+0.0001 = 0.0011
//   cov12 = 0.0002+0.0002+0.0002+0.0003+0.0002 = 0.0011
// =============================================================================
TEST(HayashiYoshidaCovTest, DualAssetPeriod1) {
    std::vector<std::vector<Real>> returns = {
        {0.01, 0.02, -0.01, 0.03, -0.02},
        {0.02, 0.01, -0.02, 0.01, -0.01}
    };

    auto result = HayashiYoshidaCov::estimate(returns, 1, false, false);
    EXPECT_EQ(result.n_assets, 2u);
    EXPECT_EQ(result.n_obs, 5u);
    ASSERT_EQ(result.cov.size(), 4u);

    // 期望值 (手算)
    constexpr Real EXP_VAR1 = 0.0019;
    constexpr Real EXP_VAR2 = 0.0011;
    constexpr Real EXP_COV12 = 0.0011;

    EXPECT_NEAR(result.cov[0], EXP_VAR1, TOL_STRICT);   // (0,0)
    EXPECT_NEAR(result.cov[3], EXP_VAR2, TOL_STRICT);   // (1,1)
    EXPECT_NEAR(result.cov[1], EXP_COV12, TOL_STRICT);  // (0,1)
    EXPECT_NEAR(result.cov[2], EXP_COV12, TOL_STRICT);  // (1,0)
}

// =============================================================================
// TEST 3: 双资产 period=2 — HY 聚合, 与 period=1 不同
//   r1 = [0.01, 0.02, 0.03, 0.04]
//   r2 = [0.01, 0.01, 0.01, 0.01]
//
//   period=1 (同步): cov12 = (1+2+3+4)*0.0001 = 0.001
//   period=2 (HY 聚合):
//     ap[0] = a[0]+a[1] = 0.03, atp[0] = 2
//     ap[1] = a[2]+a[3] = 0.07, atp[1] = 4
//     i=0: tmp_ret = b[0]+b[1] = 0.02 (bt[1]=2 == atp[0]=2), ans = 0.03*0.02 = 0.0006
//     i=1: tmp_ret = b[2]+b[3] = 0.02 (bt[3]=4 == atp[1]=4), ans = 0.07*0.02 = 0.0014
//     cov12 = 0.0006 + 0.0014 = 0.002
// =============================================================================
TEST(HayashiYoshidaCovTest, DualAssetPeriod2) {
    std::vector<std::vector<Real>> returns = {
        {0.01, 0.02, 0.03, 0.04},
        {0.01, 0.01, 0.01, 0.01}
    };

    // period=1: 标准同步协方差
    auto r1 = HayashiYoshidaCov::estimate(returns, 1, false, false);
    constexpr Real EXP_COV12_P1 = 0.001;  // (1+2+3+4)*0.0001
    EXPECT_NEAR(r1.cov[1], EXP_COV12_P1, TOL_STRICT);

    // period=2: HY 聚合
    auto r2 = HayashiYoshidaCov::estimate(returns, 2, false, false);
    constexpr Real EXP_COV12_P2 = 0.002;  // 0.0006 + 0.0014
    EXPECT_NEAR(r2.cov[1], EXP_COV12_P2, TOL_STRICT);

    // 对角线不受 period 影响 (始终为 RV)
    EXPECT_NEAR(r2.cov[0], r1.cov[0], TOL_STRICT);
    EXPECT_NEAR(r2.cov[3], r1.cov[3], TOL_STRICT);

    // period=1 和 period=2 的非对角线必须不同
    EXPECT_NE(r1.cov[1], r2.cov[1]);
}

// =============================================================================
// TEST 4: make_psd 直接测试 — 非 PSD 矩阵投影
//   输入: [[1, 2], [2, 1]] — 特征值 {3, -1} (非 PSD)
//   期望: [[1.5, 1.5], [1.5, 1.5]] — 特征值 {3, 0}
// =============================================================================
TEST(HayashiYoshidaCovTest, MakePsdProjection) {
    // 非 PSD 对称矩阵
    std::vector<Real> mat = {1.0, 2.0,
                             2.0, 1.0};

    auto psd = make_psd(mat, 2);
    ASSERT_EQ(psd.size(), 4u);

    // 期望值: Q*diag(3,0)*Q^T = [[1.5, 1.5], [1.5, 1.5]]
    EXPECT_NEAR(psd[0], 1.5, TOL_STRICT);
    EXPECT_NEAR(psd[1], 1.5, TOL_STRICT);
    EXPECT_NEAR(psd[2], 1.5, TOL_STRICT);
    EXPECT_NEAR(psd[3], 1.5, TOL_STRICT);

    // 验证 PSD: 所有特征值 >= 0
    // 2x2: lambda = (trace ± sqrt(trace^2 - 4*det)) / 2
    const Real trace = psd[0] + psd[3];
    const Real det = psd[0] * psd[3] - psd[1] * psd[2];
    const Real disc = trace * trace - 4.0 * det;
    ASSERT_GE(disc, -TOL_STRICT);  // 判别式 >= 0
    const Real sq = std::sqrt(std::max(disc, 0.0));
    const Real lambda1 = (trace + sq) / 2.0;
    const Real lambda2 = (trace - sq) / 2.0;
    EXPECT_GE(lambda1, -TOL_STRICT);
    EXPECT_GE(lambda2, -TOL_STRICT);
}

// =============================================================================
// TEST 5: 异常处理 — 空输入 / 长度不一致 / period=0
// =============================================================================
TEST(HayashiYoshidaCovTest, ExceptionHandling) {
    // 空矩阵
    std::vector<std::vector<Real>> empty;
    EXPECT_THROW(HayashiYoshidaCov::estimate(empty), std::invalid_argument);

    // 长度不一致
    std::vector<std::vector<Real>> mismatched = {
        {0.01, 0.02, 0.03},
        {0.01, 0.02}
    };
    EXPECT_THROW(HayashiYoshidaCov::estimate(mismatched), std::invalid_argument);

    // 空收益率
    std::vector<std::vector<Real>> empty_returns = {
        {},
        {}
    };
    EXPECT_THROW(HayashiYoshidaCov::estimate(empty_returns), std::invalid_argument);

    // period=0
    std::vector<std::vector<Real>> valid = {
        {0.01, 0.02, 0.03},
        {0.01, 0.02, 0.03}
    };
    EXPECT_THROW(HayashiYoshidaCov::estimate(valid, 0), std::invalid_argument);
}

// =============================================================================
// TEST 6: cor=TRUE 相关矩阵转换
//   使用 TEST 2 的数据, 验证相关系数
//   cor = cov12 / sqrt(var1 * var2) = 0.0011 / sqrt(0.0019 * 0.0011)
// =============================================================================
TEST(HayashiYoshidaCovTest, CorrelationMatrix) {
    std::vector<std::vector<Real>> returns = {
        {0.01, 0.02, -0.01, 0.03, -0.02},
        {0.02, 0.01, -0.02, 0.01, -0.01}
    };

    auto result = HayashiYoshidaCov::estimate(returns, 1, false, true);
    ASSERT_EQ(result.cov.size(), 4u);

    // 对角线 = 1
    EXPECT_NEAR(result.cov[0], 1.0, TOL_STRICT);
    EXPECT_NEAR(result.cov[3], 1.0, TOL_STRICT);

    // 相关系数 = 0.0011 / sqrt(0.0019 * 0.0011)
    const Real exp_cor = 0.0011 / std::sqrt(0.0019 * 0.0011);
    EXPECT_NEAR(result.cov[1], exp_cor, TOL_STRICT);
    EXPECT_NEAR(result.cov[2], exp_cor, TOL_STRICT);
}
