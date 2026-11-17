// =============================================================================
// test_kde.cpp - Phase 9 M1: 1D Gaussian KDE + 带宽语义 (TDD 前置)
//
// 状态: M1 已实现 (include/cpphub/timeseries/density/{bandwidth,kde}.hpp),
//       已接入 tests/CMakeLists.txt; M1 11 例 + P2-2 SJ/ISJ 6 例 = 17.
//
// 基准: probe_sklearn_kde.py (8/8 PASS) + scipy 黄金值 (n=5 手算样本)
//       + probe_sj_isj.py (P2-2: SJ=R bw.SJ ste / ISJ=KDEpy 翻译, 三组夹具)
//   - 带宽公式: Scott 1.059 / Silverman 0.9 / normal_reference C≈1.05922
//   - 数据驱动: kSj / kIsj (锚 1e-10 / 1e-8; SJ 严格尺度同变)
//   - 同核同 h: scipy gaussian_kde pdf 1e-14
//   - 语义: AbsoluteH vs DimensionlessFactor 显式区分
// =============================================================================

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <cstdint>
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

// ---- P2-2 SJ/ISJ 锚 (probe_sj_isj.py 参考实现, 确定性夹具) ----
// grid150: x_i = ((i*7919 + 12345) % 1000)/100
// lcg200 : x_i = ((s_i >> 11) % 1000003)/100000 (uint64 LCG)
// ar200  : x_i = floor(3*x_{i-1}/4) + u_i, u_i=(s_i>>33)%4
constexpr Real kSjGrid150 = 0.950442595810742;
constexpr Real kIsjGrid150 = 0.549966080404110;
constexpr Real kSjLcg200 = 0.640167818657877;
constexpr Real kIsjLcg200 = 0.319826651877444;
constexpr Real kSjAr200 = 0.081279027344518;
constexpr Real kIsjAr200 = 0.005215310828978;

std::vector<Real> grid150() {
    std::vector<Real> x(150);
    for (Size i = 0; i < x.size(); ++i) {
        x[i] = static_cast<Real>((i * 7919 + 12345) % 1000) / 100.0;
    }
    return x;
}

std::vector<Real> lcg200() {
    std::vector<Real> x(200);
    std::uint64_t s = 12345;
    for (Size i = 0; i < x.size(); ++i) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        x[i] = static_cast<Real>((s >> 11) % 1000003ULL) / 100000.0;
    }
    return x;
}

std::vector<Real> ar200() {
    std::vector<Real> x(200);
    std::uint64_t s = 12345;
    for (Size i = 0; i < x.size(); ++i) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        const Real u = static_cast<Real>((s >> 33) % 4ULL);
        x[i] = (i == 0) ? u : std::floor(3.0 * x[i - 1] / 4.0) + u;
    }
    return x;
}

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

// P2-2: SJ (Sheather-Jones 1991, R bw.SJ ste 同算法) / ISJ (Botev 2010, KDEpy 同算法)
TEST(KdeBandwidth, SjMatchesReference) {
    const auto g = grid150();
    const auto l = lcg200();
    const auto a = ar200();
    EXPECT_NEAR(de::resolve_bandwidth(std::span<const Real>(g),
                                      rule_bandwidth(de::BandwidthRule::kSj)),
                kSjGrid150, 1e-10);
    EXPECT_NEAR(de::resolve_bandwidth(std::span<const Real>(l),
                                      rule_bandwidth(de::BandwidthRule::kSj)),
                kSjLcg200, 1e-10);
    EXPECT_NEAR(de::resolve_bandwidth(std::span<const Real>(a),
                                      rule_bandwidth(de::BandwidthRule::kSj)),
                kSjAr200, 1e-10);
}

TEST(KdeBandwidth, IsjMatchesReference) {
    const auto g = grid150();
    const auto l = lcg200();
    const auto a = ar200();
    EXPECT_NEAR(de::resolve_bandwidth(std::span<const Real>(g),
                                      rule_bandwidth(de::BandwidthRule::kIsj)),
                kIsjGrid150, 1e-8);
    EXPECT_NEAR(de::resolve_bandwidth(std::span<const Real>(l),
                                      rule_bandwidth(de::BandwidthRule::kIsj)),
                kIsjLcg200, 1e-8);
    EXPECT_NEAR(de::resolve_bandwidth(std::span<const Real>(a),
                                      rule_bandwidth(de::BandwidthRule::kIsj)),
                kIsjAr200, 1e-8);
}

TEST(KdeBandwidth, SjIsjDifferFromRuleOfThumb) {
    const auto g = grid150();
    const Real sj = de::resolve_bandwidth(std::span<const Real>(g),
                                          rule_bandwidth(de::BandwidthRule::kSj));
    const Real isj = de::resolve_bandwidth(std::span<const Real>(g),
                                           rule_bandwidth(de::BandwidthRule::kIsj));
    const Real nr = de::resolve_bandwidth(std::span<const Real>(g),
                                          rule_bandwidth(de::BandwidthRule::kNormalReference));
    EXPECT_GT(std::fabs(sj - nr), 0.05 * nr) << "sj=" << sj << " nr=" << nr;
    EXPECT_GT(std::fabs(isj - nr), 0.05 * nr) << "isj=" << isj << " nr=" << nr;
}

TEST(KdeBandwidth, SjScaleEquivariant) {
    const auto g = grid150();
    std::vector<Real> scaled = g;
    for (Real& v : scaled) {
        v *= 10.0;
    }
    const Real base = de::resolve_bandwidth(std::span<const Real>(g),
                                            rule_bandwidth(de::BandwidthRule::kSj));
    const Real factor = de::resolve_bandwidth(std::span<const Real>(scaled),
                                              rule_bandwidth(de::BandwidthRule::kSj));
    EXPECT_NEAR(factor, 10.0 * base, 1e-9 * base)
        << "SJ 应随数据尺度线性缩放 (base=" << base << ")";
}

TEST(KdeValidation, SjIsjRejectConstantSample) {
    std::vector<Real> constant(50, 1.5);
    EXPECT_THROW(de::resolve_bandwidth(std::span<const Real>(constant),
                                       rule_bandwidth(de::BandwidthRule::kSj)),
                 std::invalid_argument);
    EXPECT_THROW(de::resolve_bandwidth(std::span<const Real>(constant),
                                       rule_bandwidth(de::BandwidthRule::kIsj)),
                 std::invalid_argument);
}

TEST(KdePdf, SjIsjConstructKde) {
    const auto g = grid150();
    for (de::BandwidthRule rule : {de::BandwidthRule::kSj, de::BandwidthRule::kIsj}) {
        de::KernelDensity1D kde(std::span<const Real>(g), config_with(rule_bandwidth(rule)));
        EXPECT_NEAR(kde.bandwidth(),
                    de::resolve_bandwidth(std::span<const Real>(g), rule_bandwidth(rule)),
                    1e-12);
        for (Real x : {-1.0, 0.0, 3.0, 6.0, 9.0}) {
            EXPECT_TRUE(std::isfinite(kde.pdf(x)));
            EXPECT_TRUE(std::isfinite(kde.log_pdf(x)));
        }
    }
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
// 覆盖率 CI 为 Debug -O0 + --coverage 插桩构建, 预算口径不适用 → skip
TEST(KdePerformance, N10kM1kUnderBudget) {
#if defined(CPPHUB_COVERAGE_BUILD)
    GTEST_SKIP() << "perf budget (PHASE9_SPEC D9-10) asserted on optimized "
                    "builds; coverage job is Debug -O0 + --coverage";
#endif
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
