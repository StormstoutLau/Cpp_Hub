#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/core/math.hpp"
#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>

namespace cpphub {
inline namespace v1 {

enum class VaRApproximation { Full, DeltaGamma, Delta };

struct MCVarConfig {
    Size n_paths = 100000;
    uint64_t seed = 42;
    bool antithetic = true;
    bool use_control_variate = false;
};

class MCVaR {
public:
    MCVaR(std::function<Real(const std::vector<Real>&)> portfolio_value_fn,
          std::vector<Real> current_risk_factors,
          std::vector<Real> covariance_matrix,
          Size n_factors,
          MCVarConfig config = MCVarConfig{})
        : portfolio_value_fn_(std::move(portfolio_value_fn))
        , current_factors_(std::move(current_risk_factors))
        , covariance_(std::move(covariance_matrix))
        , n_factors_(n_factors)
        , config_(config)
    {
        portfolio_value0_ = portfolio_value_fn_(current_factors_);
    }

    std::vector<Real> simulate_pnl_full() const {
        Size n = config_.n_paths;
        Size n_draws = config_.antithetic ? (n + 1) / 2 : n;
        Philox4x64 rng(config_.seed);
        std::vector<Real> pnl;
        pnl.reserve(config_.antithetic ? n_draws * 2 : n_draws);
        std::vector<Real> L = cholesky_decomp(covariance_, n_factors_);
        for (Size i = 0; i < n_draws; ++i) {
            std::vector<Real> dR = simulate_multivariate_normal(L, rng);
            std::vector<Real> factors(n_factors_);
            for (Size j = 0; j < n_factors_; ++j) {
                factors[j] = current_factors_[j] + dR[j];
            }
            Real v = portfolio_value_fn_(factors);
            pnl.push_back(v - portfolio_value0_);
            if (config_.antithetic) {
                for (Size j = 0; j < n_factors_; ++j) {
                    factors[j] = current_factors_[j] - dR[j];
                }
                v = portfolio_value_fn_(factors);
                pnl.push_back(v - portfolio_value0_);
            }
        }
        return pnl;
    }

    std::vector<Real> simulate_pnl_delta_gamma(const std::vector<Real>& delta,
                                                const std::vector<Real>& gamma) const {
        Size n = config_.n_paths;
        Size n_draws = config_.antithetic ? (n + 1) / 2 : n;
        Philox4x64 rng(config_.seed);
        std::vector<Real> pnl;
        pnl.reserve(config_.antithetic ? n_draws * 2 : n_draws);
        std::vector<Real> L = cholesky_decomp(covariance_, n_factors_);
        for (Size i = 0; i < n_draws; ++i) {
            std::vector<Real> dR = simulate_multivariate_normal(L, rng);
            Real delta_term = 0;
            for (Size j = 0; j < n_factors_; ++j) {
                delta_term += delta[j] * dR[j];
            }
            Real gamma_term = 0;
            for (Size j = 0; j < n_factors_; ++j) {
                for (Size k = 0; k < n_factors_; ++k) {
                    gamma_term += gamma[j * n_factors_ + k] * dR[j] * dR[k];
                }
            }
            pnl.push_back(delta_term + 0.5 * gamma_term);
            if (config_.antithetic) {
                delta_term = 0;
                for (Size j = 0; j < n_factors_; ++j) {
                    delta_term += delta[j] * (-dR[j]);
                }
                gamma_term = 0;
                for (Size j = 0; j < n_factors_; ++j) {
                    for (Size k = 0; k < n_factors_; ++k) {
                        gamma_term += gamma[j * n_factors_ + k] * (-dR[j]) * (-dR[k]);
                    }
                }
                pnl.push_back(delta_term + 0.5 * gamma_term);
            }
        }
        return pnl;
    }

    std::vector<Real> simulate_pnl_delta(const std::vector<Real>& delta) const {
        Size n = config_.n_paths;
        Size n_draws = config_.antithetic ? (n + 1) / 2 : n;
        Philox4x64 rng(config_.seed);
        std::vector<Real> pnl;
        pnl.reserve(config_.antithetic ? n_draws * 2 : n_draws);
        std::vector<Real> L = cholesky_decomp(covariance_, n_factors_);
        for (Size i = 0; i < n_draws; ++i) {
            std::vector<Real> dR = simulate_multivariate_normal(L, rng);
            Real delta_val = 0;
            for (Size j = 0; j < n_factors_; ++j) {
                delta_val += delta[j] * dR[j];
            }
            pnl.push_back(delta_val);
            if (config_.antithetic) {
                delta_val = 0;
                for (Size j = 0; j < n_factors_; ++j) {
                    delta_val += delta[j] * (-dR[j]);
                }
                pnl.push_back(delta_val);
            }
        }
        return pnl;
    }

    Real var(Real confidence = 0.99, VaRApproximation approx = VaRApproximation::Full,
             const std::vector<Real>& delta = {}, const std::vector<Real>& gamma = {}) {
        std::vector<Real> pnl;
        switch (approx) {
            case VaRApproximation::Full:
                pnl = simulate_pnl_full();
                break;
            case VaRApproximation::DeltaGamma:
                pnl = simulate_pnl_delta_gamma(delta, gamma);
                break;
            case VaRApproximation::Delta:
                pnl = simulate_pnl_delta(delta);
                break;
        }
        std::sort(pnl.begin(), pnl.end());
        Size idx = static_cast<Size>((1.0 - confidence) * pnl.size());
        if (idx >= pnl.size()) idx = pnl.size() - 1;
        return -pnl[idx];
    }

    Real standard_error(Real confidence, Size n_bootstrap = 1000) const {
        std::vector<Real> pnl = simulate_pnl_full();
        Size n = pnl.size();
        std::mt19937_64 gen(9999);
        std::uniform_int_distribution<Size> dist(0, n - 1);
        std::vector<Real> bootstrap_vars;
        bootstrap_vars.reserve(n_bootstrap);
        for (Size b = 0; b < n_bootstrap; ++b) {
            std::vector<Real> sample;
            sample.reserve(n);
            for (Size i = 0; i < n; ++i) {
                sample.push_back(pnl[dist(gen)]);
            }
            std::sort(sample.begin(), sample.end());
            Size idx = static_cast<Size>((1.0 - confidence) * sample.size());
            if (idx >= sample.size()) idx = sample.size() - 1;
            bootstrap_vars.push_back(-sample[idx]);
        }
        Real mean = std::accumulate(bootstrap_vars.begin(), bootstrap_vars.end(), Real(0))
                    / static_cast<Real>(n_bootstrap);
        Real var_sum = 0;
        for (auto v : bootstrap_vars) {
            Real d = v - mean;
            var_sum += d * d;
        }
        return std::sqrt(var_sum / static_cast<Real>(n_bootstrap));
    }

private:
    std::function<Real(const std::vector<Real>&)> portfolio_value_fn_;
    std::vector<Real> current_factors_;
    std::vector<Real> covariance_;
    Size n_factors_;
    MCVarConfig config_;
    Real portfolio_value0_;

    static std::vector<Real> cholesky_decomp(const std::vector<Real>& cov, Size n) {
        std::vector<Real> L(n * n, 0);
        Real jitter = 1e-10;
        for (Size j = 0; j < n; ++j) {
            Real s = 0;
            for (Size k = 0; k < j; ++k) {
                s += L[j * n + k] * L[j * n + k];
            }
            Real val = cov[j * n + j] - s;
            if (val <= 1e-15) val = jitter;
            L[j * n + j] = std::sqrt(val);
            for (Size i = j + 1; i < n; ++i) {
                Real s2 = 0;
                for (Size k = 0; k < j; ++k) {
                    s2 += L[i * n + k] * L[j * n + k];
                }
                L[i * n + j] = (cov[i * n + j] - s2) / L[j * n + j];
            }
        }
        return L;
    }

    std::vector<Real> simulate_multivariate_normal(const std::vector<Real>& L, Philox4x64& rng) const {
        std::vector<Real> Z(n_factors_);
        for (Size i = 0; i < n_factors_; i += 2) {
            uint64_t r1 = rng();
            uint64_t r2 = rng();
            double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            if (u1 == 0) u1 = 1e-16;
            auto [z1, z2] = box_muller(u1, u2);
            Z[i] = z1;
            if (i + 1 < n_factors_) Z[i + 1] = z2;
        }
        std::vector<Real> dR(n_factors_, 0);
        for (Size i = 0; i < n_factors_; ++i) {
            for (Size j = 0; j <= i; ++j) {
                dR[i] += L[i * n_factors_ + j] * Z[j];
            }
        }
        return dR;
    }
};

}  // namespace v1
}  // namespace cpphub
