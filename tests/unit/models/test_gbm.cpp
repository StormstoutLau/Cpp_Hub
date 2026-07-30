#include <gtest/gtest.h>
#include "cpphub/models/diffusion/gbm.hpp"
#include <cmath>
#include <vector>

using namespace cpphub;

TEST(GBM, PathStartsAtSpot) {
    GBMParams p{100.0, 0.05, 0.2};
    GBM gbm(p);
    Philox4x64 rng(12345);
    std::vector<Real> path(11);
    gbm.generate_path(1.0, 10, path, rng);
    EXPECT_DOUBLE_EQ(path[0], 100.0);
}

TEST(GBM, PathDimensionCorrect) {
    GBMParams p{100.0, 0.05, 0.2};
    GBM gbm(p);
    Philox4x64 rng(12345);
    Size n_steps = 10;
    std::vector<Real> path(n_steps + 1);
    gbm.generate_path(1.0, n_steps, path, rng);
    EXPECT_EQ(path.size(), n_steps + 1);
}

TEST(GBM, DeterministicWithSameSeed) {
    GBMParams p{100.0, 0.05, 0.2};
    GBM gbm(p);
    Philox4x64 rng1(42);
    Philox4x64 rng2(42);
    Size n_steps = 100;
    std::vector<Real> path1(n_steps + 1);
    std::vector<Real> path2(n_steps + 1);
    gbm.generate_path(1.0, n_steps, path1, rng1);
    gbm.generate_path(1.0, n_steps, path2, rng2);
    for (Size i = 0; i <= n_steps; ++i) {
        EXPECT_DOUBLE_EQ(path1[i], path2[i]);
    }
}

TEST(GBM, ExactSolutionMatchesExpectation) {
    GBMParams p{100.0, 0.05, 0.3};
    GBM gbm(p);
    Real T = 1.0;
    Real dt = T / 1;
    Size n_paths = 10000;
    Real sum = 0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(j);
        uint64_t r1 = rng();
        uint64_t r2 = rng();
        double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
        double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
        auto [z1, z2] = box_muller(u1, u2);
        Real ST = gbm.evolve(p.S0, dt, z1);
        sum += ST;
    }
    Real avg = sum / static_cast<Real>(n_paths);
    Real expected = p.S0 * std::exp(p.mu * T);
    Real rel_error = std::abs(avg - expected) / expected;
    EXPECT_LT(rel_error, 0.01);
}

TEST(GBM, EvolveFunctionCorrect) {
    GBMParams p{100.0, 0.05, 0.2};
    GBM gbm(p);
    Real S = 120.0;
    Real dt = 0.5;
    Real Z = 0.0;
    Real result = gbm.evolve(S, dt, Z);
    Real expected = S * std::exp((p.mu - 0.5 * p.sigma * p.sigma) * dt);
    EXPECT_DOUBLE_EQ(result, expected);
}
