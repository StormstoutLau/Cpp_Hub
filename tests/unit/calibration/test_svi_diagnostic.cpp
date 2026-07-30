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
