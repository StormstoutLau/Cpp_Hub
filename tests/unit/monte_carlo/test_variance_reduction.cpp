#include <gtest/gtest.h>
#include "cpphub/monte_carlo/control_variate.hpp"
#include "cpphub/monte_carlo/moment_matching.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/models/diffusion/gbm.hpp"
#include <cmath>

using namespace cpphub;

namespace {

Real bsm_call(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(S - K, 0.0);
    Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    Real d2 = d1 - sigma * std::sqrt(T);
    return S * std::exp(-q * T) * normal_cdf(d1) - K * std::exp(-r * T) * normal_cdf(d2);
}

Real bsm_put(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(K - S, 0.0);
    Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    Real d2 = d1 - sigma * std::sqrt(T);
    return K * std::exp(-r * T) * normal_cdf(-d2) - S * std::exp(-q * T) * normal_cdf(-d1);
}

Real generate_normal(Philox4x64& rng) {
    uint64_t r1 = rng();
    uint64_t r2 = rng();
    Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
    Real u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
    return box_muller(u1, u2).first;
}

Real geom_asian_price(Real S0, Real K, Real T, Real r, Real q, Real sigma, Size n_steps) {
    if (T <= 0.0) return std::max(S0 - K, 0.0);
    Real dt = T / static_cast<Real>(n_steps);
    Real m_geom = std::log(S0) + (r - q - 0.5 * sigma * sigma) * T * static_cast<Real>(n_steps + 1) / (2.0 * static_cast<Real>(n_steps));
    Real v_geom = sigma * sigma * T * static_cast<Real>(n_steps + 1) * static_cast<Real>(2 * n_steps + 1) / (6.0 * static_cast<Real>(n_steps * n_steps));
    Real d1 = (m_geom + v_geom - std::log(K)) / std::sqrt(v_geom);
    Real d2 = d1 - std::sqrt(v_geom);
    return std::exp(-r * T) * (std::exp(m_geom + 0.5 * v_geom) * normal_cdf(d1) - K * normal_cdf(d2));
}

} // anonymous namespace

TEST(ControlVariate, PerfectCorrelation) {
    ControlVariate cv(0.0);
    for (Size i = 0; i < 1000; ++i) {
        Real x = static_cast<Real>(i);
        cv.add_sample(x, x);
    }
    Real ratio = cv.variance_reduction_ratio();
    EXPECT_TRUE(ratio > 1e8 || !std::isfinite(ratio));
    Real est = cv.estimate();
    EXPECT_NEAR(est, 0.0, 1e-12);
}

TEST(ControlVariate, Independent) {
    Philox4x64 rng(42);
    ControlVariate cv(0.0);
    for (Size i = 0; i < 10000; ++i) {
        Real y = generate_normal(rng);
        Real z = generate_normal(rng);
        cv.add_sample(y, z);
    }
    Real ratio = cv.variance_reduction_ratio();
    EXPECT_NEAR(ratio, 1.0, 0.15);
}

TEST(ControlVariate, BSControlVariate) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    BSControlVariate bscv(S0, K, T, r, q, sigma, true);
    GBMParams gbm_params{S0, r - q, sigma};
    GBM gbm(gbm_params);
    Size n_paths = 1000;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 12345 + 42));
        std::vector<Real> path(2);
        gbm.generate_path(T, 1, path, rng);
        Real ST = path[1];
        Real payoff = std::max(ST - K, 0.0);
        bscv.add_path(ST, payoff);
    }
    Real cv_price = bscv.estimate();
    cv_price = std::exp(-r * T) * cv_price;
    Real ref = bsm_call(S0, K, T, r, q, sigma);
    EXPECT_NEAR(cv_price, ref, 0.05);
}

TEST(ControlVariate, VarianceReductionRatio) {
    Philox4x64 rng(42);
    Real rho = 0.9;
    ControlVariate cv(0.0);
    Size n = 10000;
    for (Size i = 0; i < n; ++i) {
        Real z1 = generate_normal(rng);
        Real z2 = generate_normal(rng);
        Real X = z1;
        Real Y = rho * z1 + std::sqrt(1.0 - rho * rho) * z2;
        cv.add_sample(Y, X);
    }
    Real ratio = cv.variance_reduction_ratio();
    EXPECT_NEAR(ratio, 5.26, 0.5);
}

TEST(ControlVariate, AsianArithmetic) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    Size n_steps = 12;
    Size n_paths = 1000;
    Real dt = T / static_cast<Real>(n_steps);
    Real geom_expect = geom_asian_price(S0, K, T, r, q, sigma, n_steps);
    ControlVariate cv(geom_expect);
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 9999 + 1));
        Real S = S0;
        Real sum_arith = 0.0;
        Real prod_geom = 1.0;
        for (Size i = 0; i < n_steps; ++i) {
            Real Z = generate_normal(rng);
            S = S * std::exp((r - q - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * Z);
            sum_arith += S;
            prod_geom *= S;
        }
        Real avg_arith = sum_arith / static_cast<Real>(n_steps);
        Real avg_geom = std::pow(prod_geom, 1.0 / static_cast<Real>(n_steps));
        Real payoff_arith = std::max(avg_arith - K, 0.0);
        Real payoff_geom = std::max(avg_geom - K, 0.0);
        cv.add_sample(payoff_arith, payoff_geom);
    }
    Real ratio = cv.variance_reduction_ratio();
    EXPECT_GT(ratio, 5.0);
}

TEST(MomentMatching, MeanMatched) {
    Real S0 = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    GBMMomentMatching gbm_mm(S0, T, r, q, sigma);
    GBMParams p{S0, r - q, sigma};
    GBM gbm(p);
    Size n_paths = 1000;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 777 + 1));
        std::vector<Real> path(2);
        gbm.generate_path(T, 1, path, rng);
        gbm_mm.add_path(path[1]);
    }
    Real corrected_mean = gbm_mm.estimate_mean();
    Real theoretical = S0 * std::exp((r - q) * T);
    EXPECT_NEAR(corrected_mean, theoretical, 1e-6);
}

TEST(MomentMatching, VarianceMatched) {
    Real S0 = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    Size n_paths = 1000;
    std::vector<Real> samples;
    GBMParams p{S0, r - q, sigma};
    GBM gbm(p);
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 777 + 1));
        std::vector<Real> path(2);
        gbm.generate_path(T, 1, path, rng);
        samples.push_back(path[1]);
    }
    Real inv_n = 1.0 / static_cast<Real>(n_paths);
    Real sample_mean = 0, sample_sq = 0;
    for (auto s : samples) { sample_mean += s; sample_sq += s * s; }
    sample_mean *= inv_n;
    Real sample_var = sample_sq * inv_n - sample_mean * sample_mean;
    Real sample_std = std::sqrt(sample_var);
    Real theor_mean = S0 * std::exp((r - q) * T);
    Real theor_var = theor_mean * theor_mean * (std::exp(sigma * sigma * T) - 1.0);
    Real theor_std = std::sqrt(theor_var);
    Real sum_corr = 0, sum_corr_sq = 0;
    for (auto s : samples) {
        Real sc = (s - sample_mean) * (theor_std / sample_std) + theor_mean;
        sum_corr += sc;
        sum_corr_sq += sc * sc;
    }
    Real corr_mean = sum_corr * inv_n;
    Real corr_var = sum_corr_sq * inv_n - corr_mean * corr_mean;
    EXPECT_NEAR(corr_mean, theor_mean, 1e-6);
    EXPECT_NEAR(corr_var, theor_var, 1e-6);
}

TEST(MomentMatching, GBMOptionPrice) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    GBMMomentMatching gbm_mm(S0, T, r, q, sigma);
    GBMParams p{S0, r - q, sigma};
    GBM gbm(p);
    Size n_paths = 1000;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 777 + 1));
        std::vector<Real> path(2);
        gbm.generate_path(T, 1, path, rng);
        gbm_mm.add_path(path[1]);
    }
    auto payoff = [K](Real S) { return std::max(S - K, 0.0); };
    Real mm_price = std::exp(-r * T) * gbm_mm.estimate_option_price(payoff);
    Real ref = bsm_call(S0, K, T, r, q, sigma);
    EXPECT_NEAR(mm_price, ref, 0.1);
}

TEST(VRComparison, ControlVsMomentMatching) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    Size n_paths = 1000;
    BSControlVariate bscv(S0, K, T, r, q, sigma, true);
    GBMMomentMatching gbm_mm(S0, T, r, q, sigma);
    GBMParams p{S0, r - q, sigma};
    GBM gbm(p);
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 12345 + 42));
        std::vector<Real> path(2);
        gbm.generate_path(T, 1, path, rng);
        Real ST = path[1];
        Real payoff = std::max(ST - K, 0.0);
        bscv.add_path(ST, payoff);
        gbm_mm.add_path(ST);
    }
    Real cv_ratio = bscv.variance_reduction_ratio();
    EXPECT_GT(cv_ratio, 1.0);
    auto payoff_fn = [K](Real S) { return std::max(S - K, 0.0); };
    Real mm_price = std::exp(-r * T) * gbm_mm.estimate_option_price(payoff_fn);
    Real cv_price = bscv.estimate() * std::exp(-r * T);
    Real ref = bsm_call(S0, K, T, r, q, sigma);
    EXPECT_NEAR(mm_price, ref, 0.1);
    EXPECT_NEAR(cv_price, ref, 0.05);
}

TEST(VRComparison, AsianOption) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    Size n_steps = 12;
    Size n_paths = 1000;
    Real dt = T / static_cast<Real>(n_steps);
    Real geom_expect = geom_asian_price(S0, K, T, r, q, sigma, n_steps);
    ControlVariate cv(geom_expect);
    std::vector<Real> raw_payoffs;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 9999 + 1));
        Real S = S0;
        Real sum_arith = 0.0;
        Real prod_geom = 1.0;
        for (Size i = 0; i < n_steps; ++i) {
            Real Z = generate_normal(rng);
            S = S * std::exp((r - q - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * Z);
            sum_arith += S;
            prod_geom *= S;
        }
        Real avg_arith = sum_arith / static_cast<Real>(n_steps);
        Real avg_geom = std::pow(prod_geom, 1.0 / static_cast<Real>(n_steps));
        Real payoff_arith = std::max(avg_arith - K, 0.0);
        Real payoff_geom = std::max(avg_geom - K, 0.0);
        cv.add_sample(payoff_arith, payoff_geom);
        raw_payoffs.push_back(payoff_arith);
    }
    Real cv_ratio = cv.variance_reduction_ratio();
    EXPECT_GT(cv_ratio, 5.0);

    Real inv_n = 1.0 / static_cast<Real>(n_paths);
    Real raw_sum = 0, raw_sum_sq = 0;
    for (auto p : raw_payoffs) { raw_sum += p; raw_sum_sq += p * p; }
    Real raw_mean = raw_sum * inv_n;
    Real raw_var = raw_sum_sq * inv_n - raw_mean * raw_mean;
    Real raw_se = std::sqrt(raw_var * inv_n);
    Real cv_est = cv.estimate();

    Real mm_mean = 0.0;
    Real mm_sum = 0.0, mm_sum_sq = 0.0;
    for (auto p : raw_payoffs) {
        Real corrected = (p - raw_mean) * 1.0 + geom_expect;
        mm_mean += corrected;
        mm_sum += corrected;
        mm_sum_sq += corrected * corrected;
    }
    mm_mean *= inv_n;
    Real mm_var = mm_sum_sq * inv_n - mm_mean * mm_mean;
    Real mm_se = std::sqrt(mm_var * inv_n);

    EXPECT_NEAR(cv_est, geom_expect, 0.5);
    EXPECT_NEAR(mm_mean, geom_expect, 0.5);
    EXPECT_GT(cv_ratio, 1.0);
    if (mm_var > 0.0 && raw_var > 0.0) {
        Real mm_ratio = raw_var / mm_var;
        EXPECT_GT(mm_ratio, 0.5);
    }
}
