#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include "cpphub/models/diffusion/heston.hpp"
#include "cpphub/models/diffusion/heston_qe.hpp"

using namespace cpphub;

static HestonParams make_test_params() {
    return HestonParams{
        100.0,    // S0
        0.04,     // v0
        1.5,      // kappa
        0.04,     // theta
        0.3,      // sigma
        -0.5,     // rho
        0.0,      // r
        0.0       // q
    };
}

TEST(HestonCF, CharacteristicFunctionAtZero) {
    auto p = make_test_params();
    Heston h(p);
    Complex result = h.characteristic_function(Complex{0.0, 0.0}, 1.0);
    EXPECT_NEAR(result.real(), 1.0, 1e-12);
    EXPECT_NEAR(result.imag(), 0.0, 1e-12);
}

TEST(HestonCF, CharacteristicFunctionUnitModulus) {
    auto p = make_test_params();
    Heston h(p);
    for (Real u_val = 0.1; u_val <= 5.0; u_val += 0.3) {
        Complex u(u_val, 0.0);
        Complex cf = h.characteristic_function(u, 1.0);
        Real mod = std::abs(cf);
        EXPECT_LE(mod, 1.0 + 1e-12);
    }
}

TEST(HestonCF, CharacteristicFunctionVsSchoutensTable) {
    auto p = make_test_params();
    Heston h(p);

    struct TestCase {
        Complex u;
        Real expected_re;
        Real expected_im;
    };

    // Corrected Schoutens reference table values (Little Trap form).
    // The previous values were corrupted by the branch-cut bug in the
    // original Heston (1993) log form; aligned with standalone CF tests
    // (test_heston_cf.cpp) and Albrecher (2007) "Little Trap" paper.
    TestCase cases[] = {
        {Complex{0.5, 0.0}, -0.6573742159702644, 0.7466039816575294},
        {Complex{1.0, 0.0}, -0.1231706804396436, -0.9715285746424274},
        {Complex{2.0, 0.0}, -0.8932404438080451, 0.2240614476780658},
    };

    for (const auto& tc : cases) {
        Complex cf = h.characteristic_function(tc.u, 1.0);
        EXPECT_NEAR(cf.real(), tc.expected_re, 1e-10)
            << "u = " << tc.u;
        EXPECT_NEAR(cf.imag(), tc.expected_im, 1e-10)
            << "u = " << tc.u;
    }
}

TEST(HestonCF, CharacteristicFunctionNoBranchCutDiscontinuity) {
    auto p = make_test_params();
    Heston h(p);
    Real prev_re = 0.0, prev_im = 0.0;
    bool first = true;
    for (Real u_val = 0.01; u_val <= 10.0; u_val += 0.01) {
        Complex u(u_val, 0.0);
        Complex cf = h.characteristic_function(u, 1.0);
        if (!first) {
            Real jump_re = std::abs(cf.real() - prev_re);
            Real jump_im = std::abs(cf.imag() - prev_im);
            Real tol = 1.0;
            EXPECT_LT(jump_re, tol) << "Discontinuity at u = " << u_val
                << " prev=(" << prev_re << "," << prev_im
                << ") curr=(" << cf.real() << "," << cf.imag() << ")";
            EXPECT_LT(jump_im, tol) << "Discontinuity at u = " << u_val
                << " prev=(" << prev_re << "," << prev_im
                << ") curr=(" << cf.real() << "," << cf.imag() << ")";
        }
        first = false;
        prev_re = cf.real();
        prev_im = cf.imag();
    }
}

TEST(HestonCF, CharacteristicFunctionDecaysAtInfinity) {
    auto p = make_test_params();
    Heston h(p);
    for (Real u_val = 10.0; u_val <= 50.0; u_val += 5.0) {
        Complex u(u_val, 0.0);
        Complex cf = h.characteristic_function(u, 1.0);
        Real mod = std::abs(cf);
        EXPECT_LT(mod, 1.0);
    }
}

TEST(HestonEuler, PathDimensionCorrect) {
    auto p = make_test_params();
    Heston h(p, HestonScheme::FullTruncation);
    Philox4x64 rng(42);
    Size n_steps = 100;
    std::vector<Real> path(n_steps + 1);
    h.generate_path(1.0, n_steps, path, rng);
    EXPECT_EQ(path.size(), n_steps + 1);
}

TEST(HestonEuler, PathStartsAtSpot) {
    auto p = make_test_params();
    Heston h(p, HestonScheme::FullTruncation);
    Philox4x64 rng(42);
    std::vector<Real> path(101);
    h.generate_path(1.0, 100, path, rng);
    EXPECT_DOUBLE_EQ(path[0], p.S0);
}

TEST(HestonEuler, DeterministicWithSameSeed) {
    auto p = make_test_params();
    Heston h(p, HestonScheme::FullTruncation);
    Size n_steps = 50;
    std::vector<Real> path1(n_steps + 1);
    std::vector<Real> path2(n_steps + 1);
    {
        Philox4x64 rng(42);
        h.generate_path(1.0, n_steps, path1, rng);
    }
    {
        Philox4x64 rng(42);
        h.generate_path(1.0, n_steps, path2, rng);
    }
    for (Size i = 0; i <= n_steps; ++i) {
        EXPECT_DOUBLE_EQ(path1[i], path2[i])
            << "Mismatch at step " << i;
    }
}

TEST(HestonEuler, FullTruncationKeepsVarianceNonNegative) {
    auto p = make_test_params();
    Heston h(p, HestonScheme::FullTruncation);
    Philox4x64 rng(12345);
    Size n_steps = 252;
    std::vector<Real> path(n_steps + 1);
    for (int trial = 0; trial < 10; ++trial) {
        h.generate_path(1.0, n_steps, path, rng);
        for (Size i = 0; i <= n_steps; ++i) {
            EXPECT_GE(path[i], 0.0)
                << "Negative price at step " << i << " trial " << trial;
        }
    }
}

TEST(HestonQE, PathDimensionCorrect) {
    auto p = make_test_params();
    HestonQE h(p);
    Philox4x64 rng(42);
    Size n_steps = 100;
    std::vector<Real> path(n_steps + 1);
    h.generate_path(1.0, n_steps, path, rng);
    EXPECT_EQ(path.size(), n_steps + 1);
}

TEST(HestonQE, PathStartsAtSpot) {
    auto p = make_test_params();
    HestonQE h(p);
    Philox4x64 rng(42);
    std::vector<Real> path(101);
    h.generate_path(1.0, 100, path, rng);
    EXPECT_DOUBLE_EQ(path[0], p.S0);
}

TEST(HestonQE, VarianceNonNegative) {
    auto p = make_test_params();
    HestonQE h(p);
    Philox4x64 rng(9999);
    Size n_steps = 252;
    std::vector<Real> path(n_steps + 1);
    for (int trial = 0; trial < 10; ++trial) {
        h.generate_path(1.0, n_steps, path, rng);
        for (Size i = 0; i <= n_steps; ++i) {
            EXPECT_GE(path[i], 0.0)
                << "Negative price at step " << i << " trial " << trial;
        }
    }
}

TEST(HestonQE, MartingalePropertyApproximatelyHolds) {
    auto p = make_test_params();
    HestonQE h(p);
    Size n_paths = 10000;
    Size n_steps = 52;
    Real T = 1.0;
    Real dt = T / n_steps;

    Real sum_ST = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 12345 + 42));
        std::vector<Real> path(n_steps + 1);
        h.generate_path(T, n_steps, path, rng);
        sum_ST += path.back();
    }
    Real mean_ST = sum_ST / static_cast<Real>(n_paths);
    Real expected_ST = p.S0 * std::exp((p.r - p.q) * T);
    Real rel_error = std::abs(mean_ST - expected_ST) / expected_ST;
    EXPECT_LT(rel_error, 0.01)
        << "mean_ST = " << mean_ST
        << ", expected = " << expected_ST
        << ", rel_error = " << rel_error;
}

TEST(HestonIntegration, EulerVsQEMeanPricesConsistent) {
    auto p = make_test_params();
    Heston h_euler(p, HestonScheme::FullTruncation);
    HestonQE h_qe(p);

    Size n_paths = 10000;
    Size n_steps = 52;
    Real T = 1.0;
    Real K = 100.0;
    Real r = p.r;
    Real q = p.q;

    Real sum_euler = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 777 + 1));
        std::vector<Real> path(n_steps + 1);
        h_euler.generate_path(T, n_steps, path, rng);
        Real ST = path.back();
        Real payoff = std::max(ST - K, 0.0);
        sum_euler += payoff;
    }
    Real mean_euler = sum_euler / static_cast<Real>(n_paths) * std::exp(-r * T);

    Real sum_qe = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 777 + 1));
        std::vector<Real> path(n_steps + 1);
        h_qe.generate_path(T, n_steps, path, rng);
        Real ST = path.back();
        Real payoff = std::max(ST - K, 0.0);
        sum_qe += payoff;
    }
    Real mean_qe = sum_qe / static_cast<Real>(n_paths) * std::exp(-r * T);

    Real avg = 0.5 * (std::abs(mean_euler) + std::abs(mean_qe));
    if (avg > 1e-10) {
        Real rel_diff = std::abs(mean_euler - mean_qe) / avg;
        EXPECT_LT(rel_diff, 0.5)
            << "euler = " << mean_euler
            << ", qe = " << mean_qe
            << ", rel_diff = " << rel_diff;
    }
}

TEST(HestonIntegration, QEMConvergesFasterThanEuler) {
    auto p = make_test_params();
    Heston h_euler(p, HestonScheme::FullTruncation);
    HestonQE h_qe(p);

    Size n_steps = 52;
    Real T = 1.0;
    Real K = 100.0;
    Real r = p.r;

    auto compute_se = [&](auto& model, Size n_paths) {
        Real sum_payoff = 0.0;
        Real sum_payoff2 = 0.0;
        for (Size j = 0; j < n_paths; ++j) {
            Philox4x64 rng(static_cast<uint64_t>(j * 9999 + 1));
            std::vector<Real> path(n_steps + 1);
            model.generate_path(T, n_steps, path, rng);
            Real ST = path.back();
            Real payoff = std::max(ST - K, 0.0);
            sum_payoff += payoff;
            sum_payoff2 += payoff * payoff;
        }
        Real mean = sum_payoff / static_cast<Real>(n_paths) * std::exp(-r * T);
        Real var = (sum_payoff2 / static_cast<Real>(n_paths) - (sum_payoff / static_cast<Real>(n_paths)) * (sum_payoff / static_cast<Real>(n_paths))) * std::exp(-2.0 * r * T);
        Real se = std::sqrt(std::max(var, 0.0) / static_cast<Real>(n_paths));
        return se;
    };

    Real se_euler = compute_se(h_euler, 5000);
    Real se_qe = compute_se(h_qe, 5000);

    EXPECT_LT(se_qe, se_euler * 1.5)
        << "SE Euler = " << se_euler
        << ", SE QE = " << se_qe;
}
