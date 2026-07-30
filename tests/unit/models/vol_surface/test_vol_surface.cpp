// VolSurface unit tests — bilinear / cubic-spline interpolation, SVI fit,
// calendar & butterfly arbitrage checks.
// Implemented on main station (MSVC) - 2026-07-31
#include <gtest/gtest.h>
#include "cpphub/models/vol_surface/vol_surface.hpp"
#include "cpphub/calibration/calibrator.hpp"  // bsm_implied_vol
#include "cpphub/monte_carlo/control_variate.hpp"  // bsm_call_price
#include <vector>
#include <cmath>

using namespace cpphub::v1;

namespace {
// Build a flat smile surface: constant IV = 0.2 across strikes & maturities.
VolSurface make_flat_surface() {
    std::vector<Real> strikes = {80.0, 90.0, 100.0, 110.0, 120.0};
    std::vector<Real> maturities = {0.25, 0.5, 1.0};
    std::vector<std::vector<Real>> ivs(maturities.size(),
                                        std::vector<Real>(strikes.size(), 0.20));
    VolSurface s(strikes, maturities, ivs);
    s.set_market(100.0, 0.0, 0.0);
    return s;
}

// Build a smile surface with mild skew using SVI total variance at T=1.
// w(k) = 0.04 + 0.4 * (rho*k + sqrt(k^2 + sigma^2))
VolSurface make_smile_surface() {
    std::vector<Real> strikes;
    std::vector<Real> maturities = {0.25, 0.5, 1.0};
    Real F = 100.0;
    Real a = 0.04, b = 0.4, rho = -0.2, sigma = 0.1, m = 0.0;
    for (int i = -4; i <= 4; ++i) strikes.push_back(F * std::exp(i * 0.1));
    Size nK = strikes.size();
    std::vector<std::vector<Real>> ivs(maturities.size(), std::vector<Real>(nK));
    for (Size j = 0; j < maturities.size(); ++j) {
        Real T = maturities[j];
        for (Size i = 0; i < nK; ++i) {
            Real k = std::log(strikes[i] / F);
            Real w = a + b * (rho * (k - m) + std::sqrt((k - m) * (k - m) + sigma * sigma));
            ivs[j][i] = std::sqrt(w / T);
        }
    }
    VolSurface s(strikes, maturities, ivs);
    s.set_market(F, 0.0, 0.0);
    return s;
}
}  // namespace

TEST(VolSurface, ConstructionAndValidation) {
    EXPECT_NO_THROW(make_flat_surface());

    // Strikes must be strictly ascending
    std::vector<Real> bad_strikes = {100.0, 90.0, 110.0};
    std::vector<Real> T = {1.0};
    std::vector<std::vector<Real>> ivs(1, std::vector<Real>(3, 0.2));
    EXPECT_THROW(VolSurface(bad_strikes, T, ivs), std::invalid_argument);

    // Dimension mismatch
    std::vector<Real> ok_strikes = {90.0, 100.0, 110.0};
    std::vector<std::vector<Real>> bad_ivs(1, std::vector<Real>(2, 0.2));
    EXPECT_THROW(VolSurface(ok_strikes, T, bad_ivs), std::invalid_argument);
}

TEST(VolSurface, BilinearInterpolationFlatSurface) {
    auto s = make_flat_surface();
    // All interpolations should return 0.2 exactly on flat surface
    EXPECT_NEAR(s.implied_vol(100.0, 0.5, InterpolationMethod::Bilinear), 0.20, 1e-12);
    EXPECT_NEAR(s.implied_vol(95.0, 0.375, InterpolationMethod::Bilinear), 0.20, 1e-12);
    EXPECT_NEAR(s.implied_vol(105.0, 0.75, InterpolationMethod::Bilinear), 0.20, 1e-12);
}

TEST(VolSurface, BilinearInterpolationGridPoints) {
    auto s = make_smile_surface();
    const auto& strikes = s.strikes();
    const auto& maturities = s.maturities();
    const auto& vols = s.vols();
    // At grid nodes, bilinear must reproduce input exactly
    for (Size j = 0; j < maturities.size(); ++j) {
        for (Size i = 0; i < strikes.size(); ++i) {
            Real iv = s.implied_vol(strikes[i], maturities[j],
                                     InterpolationMethod::Bilinear);
            EXPECT_NEAR(iv, vols[j][i], 1e-12)
                << "j=" << j << " i=" << i;
        }
    }
}

TEST(VolSurface, CubicSplineInterpolationGridPoints) {
    auto s = make_smile_surface();
    const auto& strikes = s.strikes();
    const auto& maturities = s.maturities();
    const auto& vols = s.vols();
    // Cubic spline must also reproduce grid points exactly
    for (Size j = 0; j < maturities.size(); ++j) {
        for (Size i = 0; i < strikes.size(); ++i) {
            Real iv = s.implied_vol(strikes[i], maturities[j],
                                     InterpolationMethod::CubicSpline);
            EXPECT_NEAR(iv, vols[j][i], 1e-10)
                << "j=" << j << " i=" << i;
        }
    }
}

TEST(VolSurface, InterpolationBetweenMethodsConsistent) {
    auto s = make_flat_surface();
    // On a flat surface both methods must agree to machine precision
    for (Real K = 85.0; K <= 115.0; K += 5.0) {
        for (Real T = 0.3; T <= 0.9; T += 0.1) {
            Real vb = s.implied_vol(K, T, InterpolationMethod::Bilinear);
            Real vc = s.implied_vol(K, T, InterpolationMethod::CubicSpline);
            EXPECT_NEAR(vb, vc, 1e-12) << "K=" << K << " T=" << T;
        }
    }
}

TEST(VolSurface, TotalVarianceAndLogMoneyness) {
    auto s = make_flat_surface();
    s.set_market(100.0, 0.0, 0.0);
    Real K = 100.0, T = 1.0;
    EXPECT_NEAR(s.total_variance(K, T), 0.04, 1e-12);  // 0.2^2 * 1
    EXPECT_NEAR(s.log_moneyness(K, T), 0.0, 1e-12);
    EXPECT_NEAR(s.log_moneyness(110.0, T), std::log(1.1), 1e-12);
}

TEST(VolSurface, CallPutPriceOnFlatSurface) {
    auto s = make_flat_surface();
    s.set_market(100.0, 0.05, 0.0);
    Real K = 100.0, T = 1.0, S = 100.0, r = 0.05, q = 0.0;
    Real C = s.call_price(K, T, S, r, q);
    Real C_ref = bsm_call_price(S, K, T, r, q, 0.20);
    EXPECT_NEAR(C, C_ref, 1e-10);
}

TEST(VolSurface, CalendarArbitrageFlatOk) {
    auto s = make_flat_surface();
    // Flat 0.2 across T — w(T) = 0.04*T strictly increasing, no arbitrage
    EXPECT_TRUE(s.check_calendar_arbitrage());
}

TEST(VolSurface, CalendarArbitrageViolationDetected) {
    // Construct a surface where shorter maturity has higher IV → w decreases in T
    std::vector<Real> strikes = {90.0, 100.0, 110.0};
    std::vector<Real> maturities = {0.25, 1.0};
    std::vector<std::vector<Real>> ivs = {
        {0.40, 0.40, 0.40},  // T=0.25, w = 0.16 * 0.25 = 0.04
        {0.10, 0.10, 0.10}   // T=1.0,  w = 0.01 * 1.0  = 0.01  -> violation
    };
    VolSurface s(strikes, maturities, ivs);
    EXPECT_FALSE(s.check_calendar_arbitrage());
}

TEST(VolSurface, ButterflyArbitrageFlatOk) {
    auto s = make_flat_surface();
    s.set_market(100.0, 0.0, 0.0);
    // Flat IV → BSM call is convex in K → no butterfly arbitrage
    EXPECT_TRUE(s.check_butterfly_arbitrage());
}

TEST(VolSurface, SVIFitRoundTrip) {
    // Generate IVs from a known SVI, fit, and check IV reconstruction at grid points.
    auto s = make_smile_surface();
    s.set_market(100.0, 0.0, 0.0);
    // Fit at T=1.0 (last maturity in the surface)
    Real T_target = 1.0;
    EXPECT_NO_THROW(s.fit_svi(T_target));
    const SVI* slice = s.svi_slice(T_target);
    ASSERT_NE(slice, nullptr);
    // At grid strikes for the closest maturity, fitted IV should match input within tolerance.
    // fit_svi uses a lightweight config (DE 30 pop × 100 gen + LM 100 iter) for speed;
    // expect ~5% accuracy which is sufficient for surface construction round-trip.
    const auto& strikes = s.strikes();
    const auto& vols = s.vols();
    Size j_closest = 0;
    Real best_dist = std::abs(s.maturities()[0] - T_target);
    for (Size j = 1; j < s.maturities().size(); ++j) {
        Real d = std::abs(s.maturities()[j] - T_target);
        if (d < best_dist) { best_dist = d; j_closest = j; }
    }
    for (Size i = 0; i < strikes.size(); ++i) {
        Real K = strikes[i];
        Real k = std::log(K / 100.0);
        Real iv_fit = slice->implied_vol(k, s.maturities()[j_closest]);
        EXPECT_NEAR(iv_fit, vols[j_closest][i], 5e-2)
            << "i=" << i << " K=" << K;
    }
}
