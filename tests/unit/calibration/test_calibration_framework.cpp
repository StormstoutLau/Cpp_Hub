#include <gtest/gtest.h>
#include "cpphub/calibration/optimizer.hpp"
#include "cpphub/models/vol_surface/svi.hpp"
#include <cmath>

using namespace cpphub::v1;

TEST(M3CompileCheck, LMQuadraticResidual) {
    // Quadratic residuals: r = [x0-2, x1-3, x0+x1-5], min at (2,3)
    ResidualFn r = [](const std::vector<Real>& x) -> std::vector<Real> {
        return {x[0] - 2.0, x[1] - 3.0, x[0] + x[1] - 5.0};
    };
    LevenbergMarquardt::Config cfg;
    cfg.max_iterations = 200;
    auto result = LevenbergMarquardt::minimize(r, {0.0, 0.0}, cfg);
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(result.x[0], 2.0, 1e-6);
    EXPECT_NEAR(result.x[1], 3.0, 1e-6);
}

TEST(M3CompileCheck, NelderMeadQuadratic) {
    ObjectiveFn f = [](const std::vector<Real>& x) -> Real {
        return (x[0] - 3.0)*(x[0] - 3.0) + (x[1] + 2.0)*(x[1] + 2.0);
    };
    NelderMead::Config cfg;
    cfg.ftol = 1e-12;
    cfg.xtol = 1e-12;
    auto result = NelderMead::minimize(f, {0.0, 0.0}, cfg);
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(result.x[0], 3.0, 1e-6);
    EXPECT_NEAR(result.x[1], -2.0, 1e-6);
}

TEST(M3CompileCheck, DEGlobalSearch) {
    ObjectiveFn f = [](const std::vector<Real>& x) -> Real {
        return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    };
    std::vector<Bounds> bounds = {{-5, 5}, {-5, 5}, {-5, 5}};
    DifferentialEvolution::Config cfg;
    cfg.max_generations = 100;
    cfg.tol = 1e-10;
    auto result = DifferentialEvolution::minimize(f, bounds, cfg);
    EXPECT_LT(std::abs(result.x[0]), 1e-3);
    EXPECT_LT(std::abs(result.x[1]), 1e-3);
    EXPECT_LT(std::abs(result.x[2]), 1e-3);
}

TEST(M3CompileCheck, SVITotalVariance) {
    SVIParams p{0.04, 0.4, 0.0, 0.1, 0.0};
    SVI s(p);
    EXPECT_NEAR(s.total_variance(0.0), 0.08, 1e-12);
    EXPECT_NEAR(s.implied_vol(0.0, 1.0), std::sqrt(0.08), 1e-12);
}

TEST(M3CompileCheck, SVINoArbitrage) {
    SVIParams p{0.04, 0.4, 0.0, 0.1, 0.0};
    SVI s(p);
    EXPECT_TRUE(s.check_butterfly_arbitrage());
    EXPECT_TRUE(s.check_butterfly_at(0.0));
}

TEST(M3CompileCheck, SVICalibration) {
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
    // RISK-013 修复: DE 100 代 × 50 种群对 5 维 SVI 不足,需增加到 500 代 × 100 种群
    // 才能可靠找到全局最优 (DE 收敛后 LM 1 次迭代即满足 gtol)
    cfg.de_pop_size = 100;
    cfg.de_generations = 500;
    cfg.lm_max_iter = 2000;
    cfg.ftol = 1e-14;
    cfg.xtol = 1e-14;
    auto result = calib_svi.calibrate(strikes, maturities, ivs, F, cfg);
    EXPECT_TRUE(result.converged);
    // Check fitted total variance matches at sample points (functional fit)
    SVI fitted(SVIParams{result.params[0], result.params[1], result.params[2],
                         result.params[3], result.params[4]});
    for (int i = 0; i < 11; ++i) {
        Real k = -0.4 + 0.08 * i;
        Real w_fit = fitted.total_variance(k);
        Real w_true = true_svi.total_variance(k);
        EXPECT_NEAR(w_fit, w_true, 1e-4) << "k=" << k;
    }
}
