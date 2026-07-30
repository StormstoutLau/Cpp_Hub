#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

using namespace cpphub;

TEST(CoreMath, ErfKnownValues) {
    EXPECT_NEAR(std::erf(0.0), 0.0, 1e-15);
    EXPECT_NEAR(std::erf(1.0), 0.8427007929497149, 1e-14);
    EXPECT_NEAR(std::erf(-1.0), -0.8427007929497149, 1e-14);
    EXPECT_NEAR(std::erf(2.0), 0.9953222650189527, 1e-14);
}

TEST(CoreMath, ErfcKnownValues) {
    EXPECT_NEAR(std::erfc(0.0), 1.0, 1e-15);
    EXPECT_NEAR(std::erfc(1.0), 0.1572992070502851, 1e-14);
}

TEST(CoreMath, NormalPdf) {
    EXPECT_NEAR(normal_pdf(0.0), INV_SQRT_2PI, 1e-15);
    EXPECT_NEAR(normal_pdf(1.0), 0.24197072451914337, 1e-14);
}

TEST(CoreMath, NormalCdfKnownValues) {
    EXPECT_NEAR(normal_cdf(0.0), 0.5, 1e-15);
    EXPECT_NEAR(normal_cdf(1.0), 0.8413447460685429, 1e-14);
    EXPECT_NEAR(normal_cdf(-1.0), 0.15865525393145707, 1e-14);
    EXPECT_NEAR(normal_cdf(3.0), 0.9986501019683699, 1e-14);
}

TEST(CoreMath, InvNormalCdfKnownValues) {
    EXPECT_NEAR(inv_normal_cdf(0.5), 0.0, 1e-10);
    EXPECT_NEAR(inv_normal_cdf(0.95), 1.6448536269514722, 1e-10);
    EXPECT_NEAR(inv_normal_cdf(0.975), 1.959963984540054, 1e-10);
    EXPECT_NEAR(inv_normal_cdf(0.99), 2.3263478740408408, 1e-10);
}

TEST(CoreMath, InvNormalCdfRoundTrip) {
    for (double p : {0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99}) {
        double x = inv_normal_cdf(p);
        double p_back = normal_cdf(x);
        EXPECT_NEAR(p, p_back, 1e-10);
    }
}

TEST(CoreMath, InvNormalCdfBoundary) {
    EXPECT_TRUE(std::isinf(inv_normal_cdf(0.0)));
    EXPECT_TRUE(std::isinf(inv_normal_cdf(1.0)));
}

TEST(CoreMath, BesselI0) {
    EXPECT_NEAR(bessel_i0(0.0), 1.0, 1e-14);
    EXPECT_NEAR(bessel_i0(1.0), 1.2660658777520082, 1e-12);
    EXPECT_NEAR(bessel_i0(5.0), 27.23987182360445, 1e-10);
}

TEST(CoreMath, BesselI1) {
    EXPECT_NEAR(bessel_i1(0.0), 0.0, 1e-14);
    EXPECT_NEAR(bessel_i1(1.0), 0.5651591039924850, 1e-12);
    EXPECT_NEAR(bessel_i1(5.0), 24.33564214245054, 1e-10);
}

TEST(CoreMath, NormalPdfSymmetry) {
    EXPECT_NEAR(normal_pdf(-1.0), normal_pdf(1.0), 1e-15);
    EXPECT_NEAR(normal_pdf(-2.0), normal_pdf(2.0), 1e-15);
}

TEST(CoreMath, NormalCdfSymmetry) {
    for (double x : {0.5, 1.0, 2.0, 3.0}) {
        EXPECT_NEAR(normal_cdf(-x), 1.0 - normal_cdf(x), 1e-14);
    }
}

TEST(CoreMath, ErfWrapperConsistency) {
    for (double x : {0.0, 0.5, 1.0, 2.0}) {
        EXPECT_NEAR(std::erf(x), cpphub::erf(x), 1e-15);
    }
    EXPECT_NEAR(std::erf(0.0) + std::erfc(0.0), 1.0, 1e-15);
}

TEST(CoreMath, InvNormalCdfMonotonic) {
    double prev = -1e100;
    for (double p = 0.01; p < 1.0; p += 0.1) {
        double x = inv_normal_cdf(p);
        EXPECT_GT(x, prev);
        prev = x;
    }
}

TEST(CoreMath, InvNormalCdfSmallP) {
    double x = inv_normal_cdf(0.001);
    EXPECT_NEAR(normal_cdf(x), 0.001, 1e-10);
}
