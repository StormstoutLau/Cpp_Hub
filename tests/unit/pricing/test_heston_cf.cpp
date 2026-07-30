#include <gtest/gtest.h>
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include <cmath>

using namespace cpphub;

TEST(HestonCFStandalone, AtZeroReturnsOne) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, 0.0};
    Complex u(0, 0);
    Real tau = 1.0;
    Complex phi = heston_characteristic_function(u, tau, 100.0, p);
    EXPECT_NEAR(std::real(phi), 1.0, 1e-12);
    EXPECT_NEAR(std::imag(phi), 0.0, 1e-12);
}

TEST(HestonCFStandalone, UnitModulusForRealU) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, 0.0};
    Real tau = 1.0;
    for (Real u_val = 0.1; u_val <= 5.0; u_val += 0.1) {
        Complex u(u_val, 0);
        Complex phi = heston_characteristic_function(u, tau, 100.0, p);
        EXPECT_LE(std::abs(phi), 1.0 + 1e-12);
    }
}

TEST(HestonCFStandalone, DecaysAtInfinity) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, 0.0};
    Real tau = 1.0;
    Complex u(50.0, 0);
    Complex phi = heston_characteristic_function(u, tau, 100.0, p);
    EXPECT_LT(std::abs(phi), 0.01);
}

TEST(HestonCFStandalone, SchoutensTableU0_5) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, 0.0};
    Complex u(0.5, 0.0);
    Real tau = 1.0;
    Complex phi = heston_characteristic_function(u, tau, 100.0, p);
    EXPECT_NEAR(std::real(phi), -0.6573742159702644, 1e-10);
    EXPECT_NEAR(std::imag(phi), 0.7466039816575294, 1e-10);
}

TEST(HestonCFStandalone, SchoutensTableU1_0) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, 0.0};
    Complex u(1.0, 0.0);
    Real tau = 1.0;
    Complex phi = heston_characteristic_function(u, tau, 100.0, p);
    EXPECT_NEAR(std::real(phi), -0.1231706804396436, 1e-10);
    EXPECT_NEAR(std::imag(phi), -0.9715285746424274, 1e-10);
}

TEST(HestonCFStandalone, SchoutensTableU2_0) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, 0.0};
    Complex u(2.0, 0.0);
    Real tau = 1.0;
    Complex phi = heston_characteristic_function(u, tau, 100.0, p);
    EXPECT_NEAR(std::real(phi), -0.8932404438080451, 1e-10);
    EXPECT_NEAR(std::imag(phi), 0.2240614476780658, 1e-10);
}

TEST(HestonCFStandalone, SchoutensTableUComplex) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, 0.0};
    Complex u(0.5, 0.5);
    Real tau = 1.0;
    Complex phi = heston_characteristic_function(u, tau, 100.0, p);
    EXPECT_LT(std::abs(phi), 1.0);
    EXPECT_TRUE(std::isfinite(std::real(phi)));
    EXPECT_TRUE(std::isfinite(std::imag(phi)));
}

TEST(HestonCFStandalone, NoBranchCutDiscontinuity) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, 0.0};
    Real tau = 1.0;
    const Real eps = 1e-8;
    for (Real u_val = 0.1; u_val <= 10.0; u_val += 0.1) {
        Complex phi_fwd = heston_characteristic_function(Complex(u_val + eps, 0), tau, 100.0, p);
        Complex phi_bwd = heston_characteristic_function(Complex(u_val - eps, 0), tau, 100.0, p);
        Complex phi_mid = heston_characteristic_function(Complex(u_val, 0), tau, 100.0, p);
        Real smoothness = std::abs(phi_fwd + phi_bwd - Real(2) * phi_mid);
        EXPECT_LT(smoothness, 1e-8);
    }
}
