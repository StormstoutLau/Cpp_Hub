// =============================================================================
// test_kde.cpp - Phase 9 M1: 1D Gaussian KDE + 带宽语义 (TDD 前置)
//
// 状态: M1 已实现 (include/cpphub/timeseries/density/{bandwidth,kde}.hpp),
//       已接入 tests/CMakeLists.txt, 10/10 PASS.
//
// 基准: probe_sklearn_kde.py (8/8 PASS) + scipy 黄金值 (n=5 手算样本)
//   - 带宽公式: Scott 1.059 / Silverman 0.9 / normal_reference C≈1.05922
//   - 同核同 h: scipy gaussian_kde pdf 1e-14
//   - 语义: AbsoluteH vs DimensionlessFactor 显式区分
// =============================================================================

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

#include "cpphub/timeseries/density/kde.hpp"

namespace de = cpphub::v1::timeseries::density;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

// 5 点样本: q25/q75 恰为数据点, 无插值歧义
const std::vector<Real> kSample = {-1.0, -0.5, 0.0, 0.5, 1.0};
constexpr Real kSd = 0.7905694150420949;        // ddof=1
constexpr Real kIqr = 1.0;
constexpr Real kA = 0.7412898443291327;        // min(sd, IQR/1.349)
constexpr Real kScott = 0.5689708405001331;    // 1.059 * A * n^-0.2
constexpr Real kSilverman = 0.48354462365450407;  // 0.9 * A * n^-0.2
constexpr Real kNormalRef = 0.5690911034345367;   // C≈1.059223840025653 * A * n^-0.2

const std::vector<Real> kGrid = {-1.0, -0.5, 0.0, 0.5, 1.0};
const std::vector<Real> kPdf = {
    0.2797878600285868, 0.3765226177459382, 0.3963462649864382,
    0.3765226177459382, 0.2797878600285868};
const std::vector<Real> kLogPdf = {
    -1.2737236057242571, -0.9767771598521500, -0.9254670432669307,
    -0.9767771598521500, -1.2737236057242571};

de::Bandwidth rule_bandwidth(de::BandwidthRule rule) {
    de::Bandwidth bw;
    bw.rule = rule;
    return bw;
}

de::Bandwidth absolute_bandwidth(Real h) {
    de::Bandwidth bw;
    bw.rule = de::BandwidthRule::kUser;
    bw.value = h;
    bw.is_factor = false;
    return bw;
}

de::Bandwidth factor_bandwidth(Real factor) {
    de::Bandwidth bw;
    bw.rule = de::BandwidthRule::kUser;
    bw.value = factor;
    bw.is_factor = true;
    return bw;
}

de::KdeConfig config_with(const de::Bandwidth& bw) {
    de::KdeConfig cfg;
    cfg.bandwidth = bw;
    return cfg;
}

}  // namespace

TEST(KdeBandwidth, ScottMatchesFormula) {
    EXPECT_NEAR(de::resolve_bandwidth(kSample, rule_bandwidth(de::BandwidthRule::kScott)),
                kScott, 1e-15);
}

TEST(KdeBandwidth, SilvermanMatchesFormula) {
    EXPECT_NEAR(de::resolve_bandwidth(kSample, rule_bandwidth(de::BandwidthRule::kSilverman)),
                kSilverman, 1e-15);
}

TEST(KdeBandwidth, NormalReferenceMatchesFrozenConstant) {
    // C 由 probe_sklearn_kde.py P4 从 statsmodels 实测反解冻结 (≈1.05922, 非 1.059)
    EXPECT_NEAR(de::resolve_bandwidth(kSample, rule_bandwidth(de::BandwidthRule::kNormalReference)),
                kNormalRef, 1e-15);
}

TEST(KdeBandwidth, AbsoluteAndFactorSemantics) {
    EXPECT_NEAR(de::resolve_bandwidth(kSample, absolute_bandwidth(0.5)), 0.5, 1e-15);
    EXPECT_NEAR(de::resolve_bandwidth(kSample, factor_bandwidth(0.5 / kSd)), 0.5, 1e-15);
}

TEST(KdePdf, MatchesScipyGolden) {
    de::KernelDensity1D kde(std::span<const Real>(kSample), config_with(absolute_bandwidth(0.5)));
    ASSERT_EQ(kde.size(), kSample.size());
    EXPECT_NEAR(kde.bandwidth(), 0.5, 1e-15);
    EXPECT_EQ(kde.rule(), de::BandwidthRule::kUser);
    for (Size i = 0; i < kGrid.size(); ++i) {
        EXPECT_NEAR(kde.pdf(kGrid[i]), kPdf[i], 1e-14) << "grid[" << i << "]";
        EXPECT_NEAR(kde.log_pdf(kGrid[i]), kLogPdf[i], 1e-14) << "grid[" << i << "]";
    }
}

TEST(KdePdf, LogPdfMatchesLogOfPdf) {
    de::KernelDensity1D kde(std::span<const Real>(kSample), config_with(absolute_bandwidth(0.5)));
    for (Real x : kGrid) {
        EXPECT_NEAR(kde.log_pdf(x), std::log(kde.pdf(x)), 1e-14);
    }
}

TEST(KdePdf, GridMatchesScalar) {
    de::KernelDensity1D kde(std::span<const Real>(kSample), config_with(absolute_bandwidth(0.5)));
    const auto out = kde.pdf(std::span<const Real>(kGrid));
    ASSERT_EQ(out.size(), kGrid.size());
    for (Size i = 0; i < kGrid.size(); ++i) {
        EXPECT_NEAR(out[i], kde.pdf(kGrid[i]), 1e-15);
    }
}

TEST(KdeSensitivity, ReturnsFinitePoints) {
    const std::vector<Real> bandwidths = {0.2, 0.5, 1.0};
    const std::vector<Real> eval_grid = {-1.0, 0.0, 1.0};
    const auto pts = de::bandwidth_sensitivity(
        std::span<const Real>(kSample),
        std::span<const Real>(bandwidths),
        std::span<const Real>(eval_grid));
    ASSERT_EQ(pts.size(), bandwidths.size());
    for (Size i = 0; i < pts.size(); ++i) {
        EXPECT_NEAR(pts[i].bandwidth, bandwidths[i], 1e-15);
        EXPECT_TRUE(std::isfinite(pts[i].log_likelihood));
        EXPECT_TRUE(std::isfinite(pts[i].pdf_at_mode));
        EXPECT_GT(pts[i].pdf_at_mode, 0.0);
    }
}

TEST(KdeValidation, RejectsInvalidInputs) {
    std::vector<Real> empty;
    std::vector<Real> one = {0.0};
    std::vector<Real> nan_sample = {0.0, std::numeric_limits<Real>::quiet_NaN()};
    const de::KdeConfig cfg = config_with(absolute_bandwidth(0.5));

    EXPECT_THROW(de::KernelDensity1D(std::span<const Real>(empty), cfg), std::invalid_argument);
    EXPECT_THROW(de::KernelDensity1D(std::span<const Real>(one), cfg), std::invalid_argument);
    EXPECT_THROW(de::KernelDensity1D(std::span<const Real>(nan_sample), cfg), std::invalid_argument);

    EXPECT_THROW(de::KernelDensity1D(std::span<const Real>(kSample),
                                    config_with(absolute_bandwidth(0.0))),
                 std::invalid_argument);
    EXPECT_THROW(de::KernelDensity1D(std::span<const Real>(kSample),
                                    config_with(absolute_bandwidth(-0.1))),
                 std::invalid_argument);
}

TEST(KdeValidation, ReflectionNotImplemented) {
    de::KdeConfig cfg = config_with(absolute_bandwidth(0.5));
    cfg.reflection = true;
    EXPECT_THROW(de::KernelDensity1D(std::span<const Real>(kSample), cfg), std::logic_error);
}

// spec D9-10 性能预算 (2026-11-11 G 门禁修订 50→150ms): 实测 52-58ms,
// UCRT 标量 exp ~5.5ns/次为硬下限, 暴力算法已最优 (inv_h 乘法/双累加器
// 无增益; exp2 变换实测慢 8× 禁用); 150ms 留 CI runner 1.5-2× 慢速余量
TEST(KdePerformance, N10kM1kUnderBudget) {
    std::mt19937_64 rng(3000);
    std::normal_distribution<Real> nd(0.0, 1.0);
    std::vector<Real> sample(10000);
    for (Real& v : sample) {
        v = nd(rng);
    }
    const de::KernelDensity1D kde(std::span<const Real>(sample),
                                  config_with(absolute_bandwidth(0.1)));

    const auto t0 = std::chrono::steady_clock::now();
    Real acc = 0.0;
    for (Size i = 0; i < 1000; ++i) {
        acc += kde.pdf(-3.0 + 6.0 * static_cast<Real>(i) / 999.0);
    }
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - t0)
                        .count();
    EXPECT_TRUE(std::isfinite(acc));  // 防优化器消除求值
    EXPECT_LT(static_cast<double>(ms), 150.0)
        << "KDE N=10k/M=1k took " << ms << " ms (budget 150 ms)";
}
