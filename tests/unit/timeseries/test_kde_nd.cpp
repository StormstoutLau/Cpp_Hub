// =============================================================================
// test_kde_nd.cpp - P2-3: 多维 Gaussian KDE (全协方差, scipy 语义)
//
// 基准: probe_multivariate_kde.py (scipy gaussian_kde 1e-14 vs Cholesky 参考)
//   - 2D fixture (n=120) / 3D fixture (n=150), 确定性 uint64 LCG
//   - scott / silverman / user factor 0.37, pdf 与 log_pdf 锚
//   - d=2 时 scott 与 silverman 因子恒等 ((d+2)/4 = 1), d=3 可区分
// =============================================================================

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include "cpphub/timeseries/density/kde_nd.hpp"

namespace de = cpphub::v1::timeseries::density;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

constexpr Real kInv2p53 = 1.0 / 9007199254740992.0;  // 2^53

Real next_unit(std::uint64_t& s) {
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<Real>(s >> 11) * kInv2p53;
}

std::vector<Real> points_2d(Size n = 120) {
    std::uint64_t su = 12345;
    std::uint64_t sv = 54321;
    std::vector<Real> pts(n * 2);
    for (Size i = 0; i < n; ++i) {
        const Real x = 10.0 * next_unit(su);
        const Real y = 0.6 * x + 2.0 * next_unit(sv);
        pts[2 * i] = x;
        pts[2 * i + 1] = y;
    }
    return pts;
}

std::vector<Real> points_3d(Size n = 150) {
    std::uint64_t su = 12345;
    std::uint64_t sv = 54321;
    std::uint64_t sw = 99999;
    std::vector<Real> pts(n * 3);
    for (Size i = 0; i < n; ++i) {
        const Real x = 10.0 * next_unit(su);
        const Real y = 0.6 * x + 2.0 * next_unit(sv);
        const Real z = -0.3 * x + 0.4 * y + 1.5 * next_unit(sw);
        pts[3 * i] = x;
        pts[3 * i + 1] = y;
        pts[3 * i + 2] = z;
    }
    return pts;
}

const std::vector<Real> kGrid2 = {0.0, 0.0, 1.0, 1.0, 2.5, 0.5,
                                  -1.0, 2.0, 8.0, 6.0};
const std::vector<Real> kGrid3 = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.5, 0.5,
                                  1.0, -1.0, 2.0, 0.0, 8.0, 6.0, 4.0};

// probe_multivariate_kde.py 冻结锚 (scipy 1.15.3 / statsmodels 0.14.4 环境)
const std::vector<Real> kPdf2Scott = {
    0.01467814638810232, 0.04751010408479, 9.640541453933522e-07,
    6.457892553210157e-05, 0.048926795518021};
const std::vector<Real> kLogPdf2Scott = {
    -4.221395531608666, -3.046812872988193, -13.852118376492706,
    -9.64762243055358, -3.0174300669759373};
const std::vector<Real> kPdf2User = {
    0.01448976296307656, 0.05209555068452817, 1.680984486689797e-08,
    9.396738660393664e-06, 0.04965612924066889};
const std::vector<Real> kPdf3Scott = {
    0.004205309048917656, 0.032082258900030754, 2.68839752340275e-06,
    1.6527145795626954e-10, 1.7413291565322156e-30};
const std::vector<Real> kLogPdf3Scott = {
    -5.471407492677005, -3.439452083873329, -12.826565258067085,
    -22.523431794151577, -68.52290408514067};
const std::vector<Real> kPdf3Silverman = {
    0.004288477633910604, 0.03304405438242568, 1.765377185412385e-06,
    5.862620689461026e-11, 2.6731203903349294e-32};
const std::vector<Real> kLogPdf3Silverman = {
    -5.451823472911608, -3.409913626525286, -13.247146187626987,
    -23.55983930268556, -72.6994765002599};

std::span<const Real> point_at(const std::vector<Real>& grid, Size dim, Size i) {
    return std::span<const Real>(grid.data() + i * dim, dim);
}

}  // namespace

TEST(KdeNd, Scott2dMatchesReference) {
    const auto pts = points_2d();
    de::KernelDensityND kde(std::span<const Real>(pts), 120, 2);
    EXPECT_NEAR(kde.factor(), std::pow(120.0, -1.0 / 6.0), 1e-15);
    for (Size i = 0; i < kPdf2Scott.size(); ++i) {
        EXPECT_NEAR(kde.pdf(point_at(kGrid2, 2, i)), kPdf2Scott[i], 1e-12)
            << "pdf[" << i << "]";
        EXPECT_NEAR(kde.log_pdf(point_at(kGrid2, 2, i)), kLogPdf2Scott[i], 1e-10)
            << "log_pdf[" << i << "]";
    }
}

TEST(KdeNd, UserFactorMatchesReference) {
    de::MultivariateKdeConfig cfg;
    cfg.rule = de::MultivariateBandwidthRule::kUserFactor;
    cfg.factor = 0.37;
    const auto pts = points_2d();
    de::KernelDensityND kde(std::span<const Real>(pts), 120, 2, cfg);
    EXPECT_NEAR(kde.factor(), 0.37, 1e-15);
    for (Size i = 0; i < kPdf2User.size(); ++i) {
        EXPECT_NEAR(kde.pdf(point_at(kGrid2, 2, i)), kPdf2User[i], 1e-12)
            << "pdf[" << i << "]";
    }
}

TEST(KdeNd, Scott3dMatchesReference) {
    const auto pts = points_3d();
    de::KernelDensityND kde(std::span<const Real>(pts), 150, 3);
    EXPECT_NEAR(kde.factor(), std::pow(150.0, -1.0 / 7.0), 1e-15);
    for (Size i = 0; i < kPdf3Scott.size(); ++i) {
        EXPECT_NEAR(kde.pdf(point_at(kGrid3, 3, i)), kPdf3Scott[i], 1e-12)
            << "pdf[" << i << "]";
        EXPECT_NEAR(kde.log_pdf(point_at(kGrid3, 3, i)), kLogPdf3Scott[i], 1e-10)
            << "log_pdf[" << i << "]";
    }
}

TEST(KdeNd, Silverman3dMatchesReference) {
    de::MultivariateKdeConfig cfg;
    cfg.rule = de::MultivariateBandwidthRule::kSilverman;
    const auto pts = points_3d();
    de::KernelDensityND kde(std::span<const Real>(pts), 150, 3, cfg);
    EXPECT_NEAR(kde.factor(), std::pow(150.0 * 5.0 / 4.0, -1.0 / 7.0), 1e-15);
    for (Size i = 0; i < kPdf3Silverman.size(); ++i) {
        EXPECT_NEAR(kde.pdf(point_at(kGrid3, 3, i)), kPdf3Silverman[i], 1e-12)
            << "pdf[" << i << "]";
        EXPECT_NEAR(kde.log_pdf(point_at(kGrid3, 3, i)), kLogPdf3Silverman[i], 1e-10)
            << "log_pdf[" << i << "]";
    }
}

TEST(KdeNd, TwoDimScottEqualsSilvermanFactor) {
    // (d+2)/4 = 1 → d=2 时两规则因子恒等 (scipy 同款性质)
    const auto pts = points_2d();
    de::MultivariateKdeConfig s;
    s.rule = de::MultivariateBandwidthRule::kSilverman;
    de::KernelDensityND a(std::span<const Real>(pts), 120, 2);
    de::KernelDensityND b(std::span<const Real>(pts), 120, 2, s);
    EXPECT_DOUBLE_EQ(a.factor(), b.factor());
}

TEST(KdeNd, LogPdfMatchesLogOfPdf) {
    const auto pts = points_3d();
    de::KernelDensityND kde(std::span<const Real>(pts), 150, 3);
    for (Size i = 0; i < kGrid3.size() / 3; ++i) {
        const auto x = point_at(kGrid3, 3, i);
        EXPECT_NEAR(kde.log_pdf(x), std::log(kde.pdf(x)), 1e-12) << "grid[" << i << "]";
    }
}

TEST(KdeNd, BandwidthMatrixIsFactorSquaredCovariance) {
    const auto pts = points_2d();
    de::KernelDensityND kde(std::span<const Real>(pts), 120, 2);
    const Real f2 = kde.factor() * kde.factor();
    Real mean0 = 0.0;
    Real mean1 = 0.0;
    for (Size i = 0; i < 120; ++i) {
        mean0 += pts[2 * i];
        mean1 += pts[2 * i + 1];
    }
    mean0 /= 120.0;
    mean1 /= 120.0;
    Real c00 = 0.0;
    Real c01 = 0.0;
    Real c11 = 0.0;
    for (Size i = 0; i < 120; ++i) {
        const Real d0 = pts[2 * i] - mean0;
        const Real d1 = pts[2 * i + 1] - mean1;
        c00 += d0 * d0;
        c01 += d0 * d1;
        c11 += d1 * d1;
    }
    c00 /= 119.0;
    c01 /= 119.0;
    c11 /= 119.0;
    EXPECT_NEAR(kde.bandwidth(0, 0), f2 * c00, 1e-12);
    EXPECT_NEAR(kde.bandwidth(0, 1), f2 * c01, 1e-12);
    EXPECT_NEAR(kde.bandwidth(1, 0), f2 * c01, 1e-12);
    EXPECT_NEAR(kde.bandwidth(1, 1), f2 * c11, 1e-12);
    EXPECT_THROW(kde.bandwidth(2, 0), std::out_of_range);
}

TEST(KdeNdValidation, RejectsInvalidInputs) {
    const auto pts = points_2d(10);
    EXPECT_THROW(de::KernelDensityND(std::span<const Real>(), 0, 2), std::invalid_argument);
    EXPECT_THROW(de::KernelDensityND(std::span<const Real>(pts), 1, 2),
                 std::invalid_argument);
    EXPECT_THROW(de::KernelDensityND(std::span<const Real>(pts), 10, 0),
                 std::invalid_argument);
    EXPECT_THROW(de::KernelDensityND(std::span<const Real>(pts), 10, 3),
                 std::invalid_argument);

    std::vector<Real> bad = pts;
    bad[3] = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(de::KernelDensityND(std::span<const Real>(bad), 10, 2),
                 std::invalid_argument);

    de::MultivariateKdeConfig user;
    user.rule = de::MultivariateBandwidthRule::kUserFactor;
    user.factor = 0.0;
    EXPECT_THROW(de::KernelDensityND(std::span<const Real>(pts), 10, 2, user),
                 std::invalid_argument);
    user.factor = -0.1;
    EXPECT_THROW(de::KernelDensityND(std::span<const Real>(pts), 10, 2, user),
                 std::invalid_argument);
}

TEST(KdeNdValidation, RejectsSingularCovariance) {
    std::vector<Real> collinear;
    for (Size i = 0; i < 20; ++i) {
        const Real v = static_cast<Real>(i);
        collinear.push_back(v);
        collinear.push_back(2.0 * v + 1.0);
    }
    EXPECT_THROW(de::KernelDensityND(std::span<const Real>(collinear), 20, 2),
                 std::invalid_argument);
}

TEST(KdeNdValidation, RejectsBadEvaluationPoint) {
    const auto pts = points_2d(10);
    de::KernelDensityND kde(std::span<const Real>(pts), 10, 2);
    const std::vector<Real> wrong_dim = {0.0, 1.0, 2.0};
    EXPECT_THROW(kde.pdf(std::span<const Real>(wrong_dim)), std::invalid_argument);
    const std::vector<Real> nan_pt = {0.0, std::numeric_limits<Real>::quiet_NaN()};
    EXPECT_THROW(kde.log_pdf(std::span<const Real>(nan_pt)), std::invalid_argument);
}
