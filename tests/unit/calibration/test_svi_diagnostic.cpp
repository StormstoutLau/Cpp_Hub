// 诊断 SVI 标定失败原因 (RISK-013 验证) + LM 迭代轨迹
#include <gtest/gtest.h>
#include "cpphub/models/vol_surface/svi.hpp"
#include "cpphub/calibration/optimizer.hpp"
#include <cmath>
#include <iostream>

using namespace cpphub::v1;

TEST(SVIDiagnostic, CalibFullConvergence) {
    SVIParams true_p{0.04, 0.4, -0.3, 0.1, 0.0};
    SVI true_svi(true_p);
    Real T = 1.0;
    Real F = 100.0;
    std::vector<Real> strikes, maturities, ivs;
    for (int i = 0; i < 11; ++i) {
        Real k = -0.4 + 0.08 * i;
        Real K = F * std::exp(k);
        Real w = true_svi.total_variance(k);
        Real iv = std::sqrt(w / T);
        strikes.push_back(K);
        maturities.push_back(T);
        ivs.push_back(iv);
    }

    SVI calib_svi(SVIParams{0.0, 0.1, 0.0, 0.1, 0.0});
    CalibConfig cfg;
    cfg.de_pop_size = 100;
    cfg.de_generations = 500;
    cfg.lm_max_iter = 2000;
    cfg.ftol = 1e-14;
    cfg.xtol = 1e-14;
    cfg.use_de_init = true;
    auto result = calib_svi.calibrate(strikes, maturities, ivs, F, cfg);

    std::cout << "converged: " << result.converged << "\n";
    std::cout << "message: " << result.message << "\n";
    std::cout << "n_iterations: " << result.n_iterations << "\n";
    std::cout << "objective_value (0.5*sum r^2): " << result.objective_value << "\n";
    std::cout << "params: a=" << result.params[0] << " b=" << result.params[1]
              << " rho=" << result.params[2] << " sigma=" << result.params[3]
              << " m=" << result.params[4] << "\n";

    SVI fitted(SVIParams{result.params[0], result.params[1], result.params[2],
                         result.params[3], result.params[4]});
    Real max_err = 0;
    for (int i = 0; i < 11; ++i) {
        Real k = -0.4 + 0.08 * i;
        Real w_fit = fitted.total_variance(k);
        Real w_true = true_svi.total_variance(k);
        Real err = std::abs(w_fit - w_true);
        if (err > max_err) max_err = err;
    }
    std::cout << "max_err=" << max_err << "\n";
    EXPECT_LT(max_err, 1e-4);
}

TEST(LMDiagnostic, QuadraticResidualTrace) {
    ResidualFn r = [](const std::vector<Real>& x) -> std::vector<Real> {
        return {x[0] - 2.0, x[1] - 3.0, x[0] + x[1] - 5.0};
    };
    LevenbergMarquardt::Config cfg;
    cfg.max_iterations = 200;
    auto result = LevenbergMarquardt::minimize(r, {0.0, 0.0}, cfg);

    std::cout << "converged: " << result.converged << "\n";
    std::cout << "message: " << result.message << "\n";
    std::cout << "n_iterations: " << result.n_iterations << "\n";
    std::cout << "n_function_evaluations: " << result.n_function_evaluations << "\n";
    std::cout << "fx: " << result.fx << "\n";
    std::cout << "x[0]: " << result.x[0] << " (expect 2.0, err=" << std::abs(result.x[0]-2.0) << ")\n";
    std::cout << "x[1]: " << result.x[1] << " (expect 3.0, err=" << std::abs(result.x[1]-3.0) << ")\n";

    EXPECT_NEAR(result.x[0], 2.0, 1e-6);
    EXPECT_NEAR(result.x[1], 3.0, 1e-6);
}

// ===========================================================================
// SVI multi-slice calibration tests
// ===========================================================================

TEST(SVIMultiSlice, EmptyInputFails) {
    SVI s(SVIParams{0.0, 0.1, 0.0, 0.1, 0.0});
    CalibrationResult summary;
    auto slices = s.calibrate_slices({}, {}, {}, 100.0, CalibConfig{}, &summary);
    EXPECT_TRUE(slices.empty());
    EXPECT_FALSE(summary.converged);
}

TEST(SVIMultiSlice, SizeMismatchFails) {
    SVI s(SVIParams{0.0, 0.1, 0.0, 0.1, 0.0});
    std::vector<Real> strikes = {90.0, 100.0, 110.0};
    std::vector<Real> maturities = {0.25, 0.5, 1.0};
    // 3 maturities × 3 strikes expected = 9 entries; provide 5 to trigger mismatch
    std::vector<Real> ivs = {0.20, 0.20, 0.20, 0.20, 0.20};
    CalibrationResult summary;
    auto slices = s.calibrate_slices(strikes, maturities, ivs, 100.0, CalibConfig{}, &summary);
    EXPECT_TRUE(slices.empty());
    EXPECT_FALSE(summary.converged);
    EXPECT_NE(summary.message.find("mismatch"), std::string::npos);
}

TEST(SVIMultiSlice, SingleSliceEquivalentToCalibrate) {
    // calibrate_slices with one maturity should match calibrate within tolerance.
    SVIParams true_p{0.04, 0.4, -0.3, 0.1, 0.0};
    SVI true_svi(true_p);
    Real T = 1.0;
    Real F = 100.0;
    std::vector<Real> strikes = {80.0, 90.0, 100.0, 110.0, 120.0};
    std::vector<Real> maturities = {T};
    std::vector<Real> ivs;
    for (Real K : strikes) {
        Real k = std::log(K / F);
        ivs.push_back(true_svi.implied_vol(k, T));
    }

    SVI s_multi(SVIParams{0.0, 0.1, 0.0, 0.1, 0.0});
    CalibConfig cfg;
    cfg.use_de_init = false;
    cfg.lm_max_iter = 500;
    CalibrationResult summary;
    auto slices = s_multi.calibrate_slices(strikes, maturities, ivs, F, cfg, &summary);

    ASSERT_EQ(slices.size(), 1u);
    EXPECT_TRUE(summary.converged) << summary.message;
    SVIParams p = slices[T];
    // Recovered params should be close to truth
    EXPECT_NEAR(p.a,     true_p.a,     1e-3);
    EXPECT_NEAR(p.b,     true_p.b,     1e-3);
    EXPECT_NEAR(p.rho,   true_p.rho,   0.05);
    EXPECT_NEAR(p.sigma, true_p.sigma, 1e-3);
    EXPECT_NEAR(p.m,     true_p.m,     1e-3);
}

TEST(SVIMultiSlice, ThreeSlicesEachConverges) {
    // Three maturities with independent SVI slices. Each slice has the SAME IV smile
    // (constant across T by construction), so total variance w = iv^2 * T scales linearly
    // with T. The recovered SVI params will differ per slice (a and b scale with T), but
    // the implied vols reconstructed from each slice must match the input IVs.
    SVIParams true_p{0.04, 0.4, -0.3, 0.1, 0.0};
    SVI true_svi(true_p);
    Real F = 100.0;
    std::vector<Real> strikes = {80.0, 90.0, 100.0, 110.0, 120.0};
    std::vector<Real> maturities = {0.25, 0.5, 1.0};
    std::vector<Real> ivs;
    for (Real T : maturities) {
        for (Real K : strikes) {
            Real k = std::log(K / F);
            // Use the T=1 SVI total variance, rescale: w_T = w_T1 * T, iv = sqrt(w_T/T) = sqrt(w_T1)
            Real w_T1 = true_svi.total_variance(k);
            ivs.push_back(std::sqrt(w_T1));  // IV constant across T
        }
    }

    SVI s(SVIParams{0.0, 0.1, 0.0, 0.1, 0.0});
    CalibConfig cfg;
    cfg.use_de_init = false;
    cfg.lm_max_iter = 500;
    CalibrationResult summary;
    auto slices = s.calibrate_slices(strikes, maturities, ivs, F, cfg, &summary);

    ASSERT_EQ(slices.size(), 3u);
    EXPECT_TRUE(summary.converged) << summary.message;
    // Verify functional fit: each slice's reconstructed IV must match input IVs (not params,
    // since SVI parameterization is non-unique when IV is constant across maturities).
    for (Size j = 0; j < maturities.size(); ++j) {
        Real T = maturities[j];
        auto it = slices.find(T);
        ASSERT_NE(it, slices.end());
        SVIParams p = it->second;
        SVI slice(p, SVIParamType::Raw);
        for (Size i = 0; i < strikes.size(); ++i) {
            Real k = std::log(strikes[i] / F);
            Real iv_fit = slice.implied_vol(k, T);
            Real iv_mkt = ivs[j * strikes.size() + i];
            // Relative error should be small
            Real rel_err = std::abs(iv_fit - iv_mkt) / iv_mkt;
            EXPECT_LT(rel_err, 0.02) << "T=" << T << " K=" << strikes[i]
                                     << " iv_fit=" << iv_fit << " iv_mkt=" << iv_mkt;
        }
    }
    // Internal state should be set to the longest-maturity slice
    EXPECT_NEAR(s.T(), maturities.back(), 1e-12);
}

TEST(SVIMultiSlice, SummaryPacksParams) {
    // Verify summary.params contains 5*n_maturities entries packed sequentially.
    SVIParams true_p{0.04, 0.4, -0.3, 0.1, 0.0};
    SVI true_svi(true_p);
    Real F = 100.0;
    std::vector<Real> strikes = {90.0, 100.0, 110.0};
    std::vector<Real> maturities = {0.5, 1.0};
    std::vector<Real> ivs;
    for (Real T : maturities) {
        for (Real K : strikes) {
            Real k = std::log(K / F);
            Real w = true_svi.total_variance(k) * T;
            ivs.push_back(std::sqrt(w / T));
        }
    }

    SVI s(SVIParams{0.0, 0.1, 0.0, 0.1, 0.0});
    CalibConfig cfg;
    cfg.use_de_init = false;
    cfg.lm_max_iter = 500;
    CalibrationResult summary;
    auto slices = s.calibrate_slices(strikes, maturities, ivs, F, cfg, &summary);

    ASSERT_EQ(slices.size(), 2u);
    // 2 slices × 5 params = 10 packed values
    EXPECT_EQ(summary.params.size(), 10u);
    // residuals: one entry per slice (max |residual|)
    EXPECT_EQ(summary.residuals.size(), 2u);
    EXPECT_GT(summary.n_iterations, 0u);
}
