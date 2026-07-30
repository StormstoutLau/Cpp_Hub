#include <gtest/gtest.h>
#include "cpphub/monte_carlo/sobol.hpp"
#include "cpphub/monte_carlo/brownian_bridge.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/core/math.hpp"
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>

using namespace cpphub::v1;

// ========== Sobol 序列测试 (6 用例) ==========

TEST(SobolSequence, FirstPointIsZero)
{
    SobolSequence seq(5);
    auto p0 = seq(0);
    for (Size i = 0; i < p0.size(); ++i) {
        EXPECT_EQ(p0[i], 0.0);
    }
}

TEST(SobolSequence, SecondPointIsHalf)
{
    SobolSequence seq(5);
    auto p1 = seq(1);
    for (Size i = 0; i < p1.size(); ++i) {
        EXPECT_DOUBLE_EQ(p1[i], 0.5);
    }
}

TEST(SobolSequence, Dimension1)
{
    SobolSequence seq(1);
    for (uint64_t n = 0; n < 64; ++n) {
        auto p = seq(n);
        // van der Corput in base 2: radical inverse of n
        double vdc = 0.0;
        double denom = 0.5;
        uint64_t m = n;
        while (m > 0) {
            if (m & 1) vdc += denom;
            m >>= 1;
            denom *= 0.5;
        }
        EXPECT_NEAR(p[0], vdc, 1e-15);
    }
}

TEST(SobolSequence, Uniformity)
{
    SobolSequence seq(3);
    const int N = 1000;
    std::vector<double> sums(3, 0.0);
    for (int i = 0; i < N; ++i) {
        auto p = seq(i);
        for (int d = 0; d < 3; ++d) {
            sums[d] += p[d];
        }
    }
    for (int d = 0; d < 3; ++d) {
        double mean = sums[d] / N;
        EXPECT_NEAR(mean, 0.5, 0.01);
    }
}

TEST(SobolSequence, Discrepancy)
{
    SobolSequence seq(2);
    const int N = 1000;
    std::vector<std::vector<double>> pts(N, std::vector<double>(2));
    for (int i = 0; i < N; ++i) {
        auto p = seq(i);
        pts[i][0] = p[0];
        pts[i][1] = p[1];
    }
    // Compute star discrepancy over grid
    const int grid = 20;
    double max_dev = 0.0;
    for (int i = 0; i <= grid; ++i) {
        for (int j = 0; j <= grid; ++j) {
            double a = static_cast<double>(i) / grid;
            double b = static_cast<double>(j) / grid;
            if (a == 0.0 || b == 0.0) continue;
            int count = 0;
            for (int k = 0; k < N; ++k) {
                if (pts[k][0] < a && pts[k][1] < b) ++count;
            }
            double frac = static_cast<double>(count) / N;
            double vol = a * b;
            double dev = std::abs(frac - vol);
            if (dev > max_dev) max_dev = dev;
        }
    }
    EXPECT_LT(max_dev, 0.02);
}

TEST(SobolSequence, Reproducibility)
{
    SobolSequence seq1(4, 42);
    SobolSequence seq2(4, 42);
    for (uint64_t n = 0; n < 100; ++n) {
        auto p1 = seq1(n);
        auto p2 = seq2(n);
        for (Size d = 0; d < p1.size(); ++d) {
            EXPECT_DOUBLE_EQ(p1[d], p2[d]);
        }
    }
}

// ========== Brownian Bridge 测试 (3 用例) ==========

TEST(BrownianBridge, Endpoints)
{
    const Size n_steps = 32;
    const Real T = 1.0;
    const int N = 10000;
    BrownianBridge bb(n_steps, T);
    SobolSequence seq(n_steps);
    double sum_WT = 0.0, sum_WT2 = 0.0;
    for (int i = 1; i <= N; ++i) {
        auto uniforms = seq(i);
        auto path = bb.generate_path(uniforms);
        Real WT = path.back();
        sum_WT += WT;
        sum_WT2 += WT * WT;
    }
    double mean_WT = sum_WT / N;
    double var_WT = sum_WT2 / N - mean_WT * mean_WT;
    EXPECT_NEAR(mean_WT, 0.0, 0.05);
    EXPECT_NEAR(var_WT, T, 0.1);
}

TEST(BrownianBridge, Increments)
{
    const Size n_steps = 16;
    const Real T = 1.0;
    const int N = 500;
    BrownianBridge bb(n_steps, T);
    SobolSequence seq(n_steps);
    for (int i = 1; i <= N; ++i) {
        auto uniforms = seq(i);
        auto incs = bb.generate(uniforms);
        auto path = bb.generate_path(uniforms);
        Real sum_incs = 0.0;
        for (auto& d : incs) sum_incs += d;
        EXPECT_NEAR(sum_incs, path.back(), 1e-12);
    }
}

TEST(BrownianBridge, VarianceReduction)
{
    const Size n_steps = 32;
    const Real T = 1.0;
    const int N = 4000;
    // Method 1: Sobol WITHOUT Brownian bridge (plain uniform-to-normal per dim)
    SobolSequence seq1(n_steps);
    double sum1 = 0.0;
    for (int i = 1; i <= N; ++i) {
        auto u = seq1(i);
        double WT = 0.0;
        for (Size j = 0; j < n_steps; ++j) {
            WT += std::sqrt(T / n_steps) * inv_normal_cdf(u[j]);
        }
        sum1 += std::exp(WT);
    }
    double est1 = sum1 / N;
    // Method 2: Sobol WITH Brownian bridge
    SobolSequence seq2(n_steps);
    BrownianBridge bb(n_steps, T);
    double sum2 = 0.0;
    for (int i = 1; i <= N; ++i) {
        auto uniforms = seq2(i);
        auto incs = bb.generate(uniforms);
        Real WT = 0.0;
        for (auto& d : incs) WT += d;
        sum2 += std::exp(WT);
    }
    double est2 = sum2 / N;
    double exact = std::exp(0.5 * T);
    double err1 = std::abs(est1 - exact);
    double err2 = std::abs(est2 - exact);
    // Variance reduction means BB gives much smaller error
    EXPECT_GT(err1, err2 * 5.0);
}

// ========== QMC 集成测试 (3 用例) ==========

TEST(QMC, IntegrationE)
{
    const int N = 65536;
    SobolSequence seq(1);
    BrownianBridge bb(1, 1.0);
    double sum = 0.0;
    for (int i = 1; i <= N; ++i) {
        auto uniforms = seq(i);
        auto path = bb.generate_path(uniforms);
        sum += std::exp(path.back());
    }
    double estimate = sum / N;
    double exact = std::exp(0.5);
    EXPECT_NEAR(estimate, exact, 1e-3);
}

TEST(QMC, AsianOption)
{
    const Size n_steps = 12;
    const Real T = 1.0;
    const Real S0 = 100.0;
    const Real K = 100.0;
    const Real r = 0.05;
    const Real sigma = 0.2;
    const int N = 1048576;
    SobolSequence seq(n_steps);
    BrownianBridge bb(n_steps, T);
    double sum_payoff = 0.0;
    for (int i = 1; i <= N; ++i) {
        auto uniforms = seq(i);
        auto path = bb.generate_path(uniforms);
        double prod = 1.0;
        for (Size j = 0; j < n_steps; ++j) {
            double t = static_cast<double>(j + 1) / n_steps * T;
            double St = S0 * std::exp((r - 0.5 * sigma * sigma) * t + sigma * path[j]);
            prod *= St;
        }
        double geom_avg = std::pow(prod, 1.0 / n_steps);
        double payoff = std::exp(-r * T) * std::max(geom_avg - K, 0.0);
        sum_payoff += payoff;
    }
    double mc_price = sum_payoff / N;
    // Analytic price: geometric Asian option in GBM
    // G = S0 * exp((r - sigma^2/2) * mean_t + sigma * mean_W)
    // ln(G) ~ N(ln(S0) + (r - sigma^2/2)*T_adj, sigma^2 * T_adj)
    double dt = T / n_steps;
    double sum_t = 0.0;
    for (Size j = 0; j < n_steps; ++j) {
        sum_t += static_cast<double>(j + 1) * dt;
    }
    double T_avg = sum_t / n_steps;
    // Var(mean_W) = Var(1/n * sum(W(t_i)))
    // For BB construction, W(t_i) = Brownian motion at t_i
    // Cov(W(t_i), W(t_j)) = min(t_i, t_j)
    double var_meanW = 0.0;
    for (Size i = 0; i < n_steps; ++i) {
        for (Size j = 0; j < n_steps; ++j) {
            double ti = static_cast<double>(i + 1) * dt;
            double tj = static_cast<double>(j + 1) * dt;
            var_meanW += std::min(ti, tj);
        }
    }
    var_meanW /= (n_steps * n_steps);
    double mu_adj = (r - 0.5 * sigma * sigma) * T_avg;
    double sigma_adj = sigma * std::sqrt(var_meanW);
    double d1 = (std::log(S0 / K) + mu_adj + sigma_adj * sigma_adj) / sigma_adj;
    double d2 = (std::log(S0 / K) + mu_adj) / sigma_adj;
    double analytic = std::exp(-r * T) * (S0 * std::exp(mu_adj + 0.5 * sigma_adj * sigma_adj) * normal_cdf(d1) - K * normal_cdf(d2));
    EXPECT_NEAR(mc_price, analytic, 1e-4);
}

TEST(QMC, VarianceReductionVsPseudoRandom)
{
    const int N = 1024;
    const Real T = 1.0;
    double exact = std::exp(0.5 * T);
    // PRNG estimator variance: Var[exp(Z)] = exp(2) - exp(1) for Z~N(0,1)
    // Theoretical: Var = exp(2*1) - exp(1)^2 = e^2 - e
    double theo_var = 4.67077427047161;  // exp(2) - exp(1)
    // QMC MSE (squared error from true value)
    SobolSequence seq(1);
    BrownianBridge bb(1, T);
    double sum_qmc = 0.0;
    for (int i = 1; i <= N; ++i) {
        auto uniforms = seq(i);
        auto path = bb.generate_path(uniforms);
        sum_qmc += std::exp(path.back());
    }
    double est_qmc = sum_qmc / N;
    double mse_qmc = (est_qmc - exact) * (est_qmc - exact);
    EXPECT_GT(theo_var / mse_qmc / N, 10.0);
}
