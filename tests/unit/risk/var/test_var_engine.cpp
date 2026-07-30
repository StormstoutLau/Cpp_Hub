#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/risk/var/historical_var.hpp"
#include "cpphub/risk/var/parametric_var.hpp"
#include "cpphub/risk/var/mc_var.hpp"
#include "cpphub/risk/var/expected_shortfall.hpp"
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>

using namespace cpphub::v1;

// ========== HistoricalVaR Tests (6) ==========

TEST(HistoricalVaR, BasicQuantile) {
    std::mt19937_64 gen(42);
    std::normal_distribution<Real> dist(0, 1);
    std::vector<Real> pnl;
    for (int i = 0; i < 10000; ++i) {
        pnl.push_back(dist(gen));
    }
    HistoricalVaR hv(pnl, 0.99, 1);
    Real var_val = hv.var(QuantileInterpolation::Linear);
    EXPECT_NEAR(var_val, 2.33, 0.3);
}

TEST(HistoricalVaR, LinearInterpolation) {
    std::vector<Real> pnl;
    for (int i = 1; i <= 100; ++i) {
        pnl.push_back(static_cast<Real>(i));
    }
    HistoricalVaR hv(pnl, 0.95, 1);
    Real var_val = hv.var(QuantileInterpolation::Linear);
    EXPECT_NEAR(var_val, 95.05, 0.01);
}

TEST(HistoricalVaR, ConservativeInterpolation) {
    std::vector<Real> pnl;
    for (int i = 1; i <= 100; ++i) {
        pnl.push_back(static_cast<Real>(i));
    }
    HistoricalVaR hv(pnl, 0.95, 1);
    Real var_val = hv.var(QuantileInterpolation::Conservative);
    EXPECT_NEAR(var_val, 96.0, 0.01);
}

TEST(HistoricalVaR, RollingWindow) {
    std::vector<Real> pnl(200, 0);
    for (Size i = 0; i < 200; ++i) {
        pnl[i] = static_cast<Real>(i % 50);
    }
    auto rolling = HistoricalVaR::rolling_var(pnl, 50, 0.95,
                                              QuantileInterpolation::Linear);
    EXPECT_EQ(rolling.size(), 150);
    EXPECT_GT(rolling[0], 0);
}

TEST(HistoricalVaR, WeightedVaR) {
    std::mt19937_64 gen(42);
    std::normal_distribution<Real> dist(0, 1);
    std::vector<Real> pnl;
    for (int i = 0; i < 500; ++i) {
        pnl.push_back(dist(gen));
    }
    HistoricalVaR hv(pnl, 0.99, 1);
    Real w_var = hv.weighted_var(0.99, QuantileInterpolation::Linear);
    Real plain_var = hv.var(QuantileInterpolation::Linear);
    EXPECT_NE(w_var, plain_var);
}

TEST(HistoricalVaR, BootstrapCI) {
    std::mt19937_64 gen(42);
    std::normal_distribution<Real> dist(0, 1);
    std::vector<Real> pnl;
    for (int i = 0; i < 2000; ++i) {
        pnl.push_back(dist(gen));
    }
    HistoricalVaR hv(pnl, 0.99, 1);
    auto ci = hv.bootstrap_ci(500, 0.95);
    EXPECT_LT(ci.first, ci.second);
    EXPECT_GT(ci.first, 0);
    EXPECT_GT(ci.second, 0);
}

// ========== ParametricVaR Tests (8) ==========

TEST(ParametricVaR, NormalVaR) {
    PortfolioStats stats;
    stats.mean = 0;
    stats.variance = 1;
    ParametricVaR pvar(stats, 0.99, 1);
    Real var_val = pvar.var(ParametricMethod::Normal);
    EXPECT_NEAR(var_val, 2.326, 0.01);
}

TEST(ParametricVaR, NormalVaRWithDrift) {
    PortfolioStats stats;
    stats.mean = 0.05;
    stats.variance = 1;
    ParametricVaR pvar(stats, 0.99, 1);
    Real var_val = pvar.var(ParametricMethod::Normal);
    EXPECT_NEAR(var_val, -(0.05 + inv_normal_cdf(0.01)), 0.01);
}

TEST(ParametricVaR, StudentTVaR) {
    PortfolioStats stats;
    stats.mean = 0;
    stats.variance = 1;
    stats.degrees_of_freedom = 5;
    ParametricVaR pvar(stats, 0.99, 1);
    Real normal_var = pvar.var(ParametricMethod::Normal);
    Real t_var = pvar.var(ParametricMethod::StudentT);
    EXPECT_GT(t_var, normal_var);
}

TEST(ParametricVaR, StudentTQuantile) {
    PortfolioStats stats;
    stats.mean = 0;
    stats.variance = 1;
    stats.degrees_of_freedom = 5;
    ParametricVaR pvar(stats, 0.99, 1);
    Real var_val = pvar.var(ParametricMethod::StudentT);
    EXPECT_NEAR(var_val, 2.607, 0.05);
}

TEST(ParametricVaR, CornishFisherPositiveSkew) {
    PortfolioStats stats;
    stats.mean = 0;
    stats.variance = 1;
    stats.skewness = 0.5;
    stats.kurtosis = 0;
    ParametricVaR pvar(stats, 0.99, 1);
    Real cf_var = pvar.var(ParametricMethod::CornishFisher);
    Real normal_var = pvar.var(ParametricMethod::Normal);
    EXPECT_LT(cf_var, normal_var);
}

TEST(ParametricVaR, CornishFisherFatTail) {
    PortfolioStats stats;
    stats.mean = 0;
    stats.variance = 1;
    stats.skewness = 0;
    stats.kurtosis = 2;
    ParametricVaR pvar(stats, 0.99, 1);
    Real cf_var = pvar.var(ParametricMethod::CornishFisher);
    Real normal_var = pvar.var(ParametricMethod::Normal);
    EXPECT_GT(cf_var, normal_var);
}

TEST(ParametricVaR, EstimateStats) {
    std::mt19937_64 gen(42);
    std::normal_distribution<Real> dist(0.01, 0.02);
    std::vector<Real> pnl;
    for (int i = 0; i < 10000; ++i) {
        pnl.push_back(dist(gen));
    }
    PortfolioStats stats = ParametricVaR::estimate_stats(pnl);
    EXPECT_NEAR(stats.mean, 0.01, 0.005);
    EXPECT_NEAR(stats.variance, 0.0004, 0.0002);
    EXPECT_NEAR(stats.skewness, 0, 0.5);
    EXPECT_NEAR(stats.kurtosis, 0, 1.0);
}

TEST(ParametricVaR, FromReturns) {
    Size n_assets = 2;
    Size n_obs = 1000;
    std::mt19937_64 gen(42);
    std::normal_distribution<Real> dist1(0.0001, 0.01);
    std::normal_distribution<Real> dist2(0.0002, 0.015);
    std::vector<std::vector<Real>> returns(n_assets, std::vector<Real>(n_obs));
    for (Size i = 0; i < n_obs; ++i) {
        returns[0][i] = dist1(gen);
        returns[1][i] = dist2(gen);
    }
    std::vector<Real> weights = {0.6, 0.4};
    PortfolioStats stats = ParametricVaR::from_returns(returns, weights);
    EXPECT_NEAR(stats.mean, 0.6 * 0.0001 + 0.4 * 0.0002, 0.001);
    EXPECT_GT(stats.variance, 0);
}

// ========== MCVaR Tests (7) ==========

TEST(MCVaR, FullRevaluation) {
    Real S0 = 100;
    auto payoff = [S0](const std::vector<Real>& S) -> Real {
        return S[0];
    };
    std::vector<Real> current = {S0};
    std::vector<Real> cov = {0.01};
    Size n_factors = 1;
    MCVarConfig cfg;
    cfg.n_paths = 100000;
    cfg.seed = 1234;
    cfg.antithetic = true;
    MCVaR mc(payoff, current, cov, n_factors, cfg);
    Real mc_var = mc.var(0.99, VaRApproximation::Full);
    PortfolioStats stats;
    stats.mean = 0;
    stats.variance = 0.01;
    ParametricVaR pvar(stats, 0.99, 1);
    Real normal_var = pvar.var(ParametricMethod::Normal);
    EXPECT_NEAR(mc_var, normal_var, normal_var * 0.05);
}

TEST(MCVaR, DeltaApprox) {
    Real S0 = 100;
    auto payoff = [S0](const std::vector<Real>& S) -> Real {
        return S[0];
    };
    std::vector<Real> current = {S0};
    std::vector<Real> cov = {0.01};
    Size n_factors = 1;
    MCVarConfig cfg;
    cfg.n_paths = 100000;
    cfg.seed = 1234;
    cfg.antithetic = true;
    MCVaR mc(payoff, current, cov, n_factors, cfg);
    std::vector<Real> delta = {1.0};
    Real mc_var = mc.var(0.99, VaRApproximation::Delta, delta);
    PortfolioStats stats;
    stats.mean = 0;
    stats.variance = 0.01;
    ParametricVaR pvar(stats, 0.99, 1);
    Real normal_var = pvar.var(ParametricMethod::Normal);
    EXPECT_NEAR(mc_var, normal_var, normal_var * 0.02);
}

TEST(MCVaR, DeltaGammaApprox) {
    Real K = 100;
    Real S0 = 100;
    auto payoff = [K](const std::vector<Real>& S) -> Real {
        Real val = S[0] - K;
        return val > 0 ? val : 0;
    };
    std::vector<Real> current = {S0};
    std::vector<Real> cov = {0.01};
    Size n_factors = 1;
    MCVarConfig cfg;
    cfg.n_paths = 100000;
    cfg.seed = 1234;
    cfg.antithetic = false;
    MCVaR mc(payoff, current, cov, n_factors, cfg);
    std::vector<Real> delta = {0.5};
    std::vector<Real> gamma = {0.02};
    Real full_var = mc.var(0.99, VaRApproximation::Full);
    Real delta_var = mc.var(0.99, VaRApproximation::Delta, delta);
    Real dg_var = mc.var(0.99, VaRApproximation::DeltaGamma, delta, gamma);
    Real full_minus_dg = std::abs(full_var - dg_var);
    Real full_minus_delta = std::abs(full_var - delta_var);
    EXPECT_LE(full_minus_dg, full_minus_delta);
}

TEST(MCVaR, AntitheticVarianceReduction) {
    Real S0 = 100;
    auto payoff = [S0](const std::vector<Real>& S) -> Real {
        return S[0];
    };
    std::vector<Real> current = {S0};
    std::vector<Real> cov = {0.01};
    Size n_factors = 1;
    MCVarConfig cfg_plain;
    cfg_plain.n_paths = 50000;
    cfg_plain.seed = 42;
    cfg_plain.antithetic = false;
    MCVaR mc_plain(payoff, current, cov, n_factors, cfg_plain);
    auto pnl_plain = mc_plain.simulate_pnl_full();
    Real mean_plain = std::accumulate(pnl_plain.begin(), pnl_plain.end(), Real(0)) / pnl_plain.size();
    Real var_plain = 0;
    for (auto v : pnl_plain) {
        Real d = v - mean_plain;
        var_plain += d * d;
    }
    var_plain /= pnl_plain.size();

    MCVarConfig cfg_at;
    cfg_at.n_paths = 50000;
    cfg_at.seed = 42;
    cfg_at.antithetic = true;
    MCVaR mc_at(payoff, current, cov, n_factors, cfg_at);
    auto pnl_at = mc_at.simulate_pnl_full();
    Real mean_at = std::accumulate(pnl_at.begin(), pnl_at.end(), Real(0)) / pnl_at.size();
    Real var_at = 0;
    for (auto v : pnl_at) {
        Real d = v - mean_at;
        var_at += d * d;
    }
    var_at /= pnl_at.size();

    EXPECT_LT(var_at, var_plain);
}

TEST(MCVaR, MultiVariateNormal) {
    Size nf = 2;
    std::vector<Real> cov = {1.0, 0.5, 0.5, 1.0};
    std::vector<Real> current = {0, 0};
    auto payoff = [](const std::vector<Real>& S) -> Real {
        return S[0] + S[1];
    };
    MCVarConfig cfg;
    cfg.n_paths = 100000;
    cfg.seed = 42;
    cfg.antithetic = false;
    MCVaR mc(payoff, current, cov, nf, cfg);
    auto pnl = mc.simulate_pnl_full();
    Real mean = std::accumulate(pnl.begin(), pnl.end(), Real(0)) / pnl.size();
    EXPECT_NEAR(mean, 0, 0.1);
}

TEST(MCVaR, StandardError) {
    Real S0 = 100;
    auto payoff = [S0](const std::vector<Real>& S) -> Real {
        return S[0];
    };
    std::vector<Real> current = {S0};
    std::vector<Real> cov = {0.01};
    Size n_factors = 1;
    MCVarConfig cfg;
    cfg.n_paths = 100000;
    cfg.seed = 42;
    cfg.antithetic = true;
    MCVaR mc(payoff, current, cov, n_factors, cfg);
    Real mc_var = mc.var(0.99, VaRApproximation::Full);
    Real se = mc.standard_error(0.99, 200);
    EXPECT_LT(se, mc_var * 0.01);
}

TEST(MCVaR, Reproducibility) {
    Real S0 = 100;
    auto payoff = [S0](const std::vector<Real>& S) -> Real {
        return S[0];
    };
    std::vector<Real> current = {S0};
    std::vector<Real> cov = {0.01};
    Size n_factors = 1;
    MCVarConfig cfg;
    cfg.n_paths = 50000;
    cfg.seed = 42;
    cfg.antithetic = true;
    MCVaR mc1(payoff, current, cov, n_factors, cfg);
    MCVaR mc2(payoff, current, cov, n_factors, cfg);
    Real v1 = mc1.var(0.99, VaRApproximation::Full);
    Real v2 = mc2.var(0.99, VaRApproximation::Full);
    EXPECT_EQ(v1, v2);
}

// ========== ExpectedShortfall Tests (4) ==========

TEST(ES, NormalES) {
    ExpectedShortfall es;
    Real es_val = es.normal_es(0, 1, 0.99);
    Real z = inv_normal_cdf(0.01);
    Real phi_z = normal_pdf(z);
    Real expected = phi_z / 0.01;
    EXPECT_NEAR(es_val, expected, 0.001);
    EXPECT_GT(es_val, 2.326);
    EXPECT_NEAR(es_val, 2.665, 0.01);
}

TEST(ES, StudentTES) {
    ExpectedShortfall es;
    Real es_normal = es.normal_es(0, 1, 0.99);
    Real es_student = es.student_t_es(0, 1, 5, 0.99);
    EXPECT_GT(es_student, es_normal);
}

TEST(ES, FromLosses) {
    std::mt19937_64 gen(42);
    std::normal_distribution<Real> dist(0, 1);
    std::vector<Real> losses;
    for (int i = 0; i < 50000; ++i) {
        losses.push_back(-dist(gen));
    }
    ExpectedShortfall es;
    Real sample_es = es.from_losses(losses, 0.99);
    Real z = inv_normal_cdf(0.01);
    Real phi_z = normal_pdf(z);
    Real theoretical_es = phi_z / 0.01;
    EXPECT_NEAR(sample_es, theoretical_es, 0.3);
}

TEST(ES, TailAverage) {
    std::mt19937_64 gen(42);
    std::normal_distribution<Real> dist(0, 1);
    std::vector<Real> losses;
    for (int i = 0; i < 10000; ++i) {
        losses.push_back(-dist(gen));
    }
    ExpectedShortfall es;
    Real es_val = es.tail_average(losses, 0.99, 100);
    Real z = inv_normal_cdf(0.01);
    Real phi_z = normal_pdf(z);
    Real theoretical = phi_z / 0.01;
    EXPECT_NEAR(es_val, theoretical, 0.5);
}
