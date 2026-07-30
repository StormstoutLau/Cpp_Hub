#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace cpphub {
inline namespace v1 {

enum class ParametricMethod { Normal, StudentT, CornishFisher };

struct PortfolioStats {
    Real mean{};
    Real variance{};
    Real skewness{};
    Real kurtosis{};
    Real degrees_of_freedom{};
};

class ParametricVaR {
public:
    ParametricVaR(PortfolioStats stats, Real confidence = 0.99, Size horizon_days = 1)
        : stats_(stats), confidence_(confidence), horizon_days_(horizon_days) {}

    Real var(ParametricMethod method) const {
        switch (method) {
            case ParametricMethod::Normal: return normal_var();
            case ParametricMethod::StudentT: return student_t_var();
            case ParametricMethod::CornishFisher: return cornish_fisher_var();
        }
        return 0;
    }

    static PortfolioStats estimate_stats(const std::vector<Real>& pnl_history) {
        PortfolioStats stats;
        Size n = pnl_history.size();
        if (n == 0) return stats;
        Real mean = std::accumulate(pnl_history.begin(), pnl_history.end(), Real(0)) / n;
        Real var = 0, skew = 0, kurt = 0;
        for (auto x : pnl_history) {
            Real d = x - mean;
            var += d * d;
            skew += d * d * d;
            kurt += d * d * d * d;
        }
        var /= n;
        Real stddev = std::sqrt(var);
        skew = skew / n / (stddev * stddev * stddev);
        kurt = kurt / n / (var * var) - 3.0;
        stats.mean = mean;
        stats.variance = var;
        stats.skewness = skew;
        stats.kurtosis = kurt;
        stats.degrees_of_freedom = estimate_t_dof(kurt);
        return stats;
    }

    static Real estimate_t_dof(Real kurtosis) {
        if (kurtosis <= 0) return 30;
        Real dof = 4.0 / kurtosis + 4.0;
        if (dof < 2.1) dof = 2.1;
        return dof;
    }

    static PortfolioStats from_returns(const std::vector<std::vector<Real>>& asset_returns,
                                        const std::vector<Real>& weights) {
        PortfolioStats stats;
        Size n_assets = weights.size();
        Size n_obs = asset_returns.empty() ? 0 : asset_returns[0].size();
        if (n_assets == 0 || n_obs == 0) return stats;
        std::vector<Real> portfolio_returns(n_obs, 0);
        for (Size i = 0; i < n_obs; ++i) {
            for (Size j = 0; j < n_assets; ++j) {
                portfolio_returns[i] += weights[j] * asset_returns[j][i];
            }
        }
        return estimate_stats(portfolio_returns);
    }

private:
    PortfolioStats stats_;
    Real confidence_;
    Size horizon_days_;

    Real normal_var() const {
        Real z = inv_normal_cdf(1.0 - confidence_);
        return -(stats_.mean + z * std::sqrt(stats_.variance * static_cast<Real>(horizon_days_)));
    }

    Real student_t_var() const {
        Real dof = stats_.degrees_of_freedom;
        if (dof <= 2) dof = 2.1;
        Real t = student_t_quantile(1.0 - confidence_, dof);
        Real scale = std::sqrt(stats_.variance * (dof - 2.0) / dof * static_cast<Real>(horizon_days_));
        return -(stats_.mean + t * scale);
    }

    Real cornish_fisher_var() const {
        Real z = inv_normal_cdf(1.0 - confidence_);
        Real z2 = z * z;
        Real z3 = z2 * z;
        Real skew = stats_.skewness;
        Real kurt = stats_.kurtosis;
        Real z_cf = z + (z2 - 1.0) * skew / 6.0
                     + (z3 - 3.0 * z) * kurt / 24.0
                     - (2.0 * z3 - 5.0 * z) * skew * skew / 36.0;
        return -(stats_.mean + z_cf * std::sqrt(stats_.variance * static_cast<Real>(horizon_days_)));
    }

    static Real betainc(Real x, Real a, Real b) {
        if (x <= 0) return 0;
        if (x >= 1) return 1;
        if (x > (a + 1.0) / (a + b + 2.0)) {
            return 1.0 - betainc(1.0 - x, b, a);
        }
        Real lbeta = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
        Real front = std::exp(a * std::log(x) + b * std::log(1.0 - x) - lbeta) / a;
        Real f = 1.0;
        Real C = 1.0;
        Real D = 0.0;
        Real tiny = 1e-30;
        for (int m = 1; m <= 200; ++m) {
            Real a_m;
            if (m % 2 == 0) {
                int k = m / 2;
                a_m = static_cast<Real>(k) * (b - static_cast<Real>(k)) * x
                      / ((a + static_cast<Real>(m) - 1.0) * (a + static_cast<Real>(m)));
            } else {
                int k = (m - 1) / 2;
                a_m = -(a + static_cast<Real>(k)) * (a + b + static_cast<Real>(k)) * x
                      / ((a + static_cast<Real>(m) - 1.0) * (a + static_cast<Real>(m)));
            }
            D = 1.0 + a_m * D;
            if (D == 0) D = tiny;
            C = 1.0 + a_m / C;
            if (C == 0) C = tiny;
            D = 1.0 / D;
            Real delta = C * D;
            f *= delta;
            if (std::abs(delta - 1.0) < 1e-14) break;
        }
        return front / f;
    }

    Real student_t_quantile(Real p, Real dof) const {
        Real x = inv_normal_cdf(p);
        for (int iter = 0; iter < 100; ++iter) {
            Real t = x;
            bool neg = (t < 0);
            Real abs_t = std::abs(t);
            Real cdf_val;
            if (abs_t < 1e-15) {
                cdf_val = 0.5;
            } else {
                Real x_beta = dof / (dof + t * t);
                Real ib = betainc(x_beta, dof / 2.0, 0.5);
                cdf_val = neg ? 0.5 * ib : 1.0 - 0.5 * ib;
            }
            Real f = cdf_val - p;
            if (std::abs(f) < 1e-15) break;
            Real pdf_val = std::exp(std::lgamma((dof + 1.0) / 2.0) - std::lgamma(dof / 2.0))
                           * std::pow(1.0 + t * t / dof, -(dof + 1.0) / 2.0)
                           / std::sqrt(dof * PI);
            if (std::abs(pdf_val) < 1e-30) break;
            x -= f / pdf_val;
        }
        return x;
    }

    Real student_t_pdf(Real x, Real dof) const {
        Real t1 = std::exp(std::lgamma((dof + 1.0) / 2.0) - std::lgamma(dof / 2.0));
        Real t2 = std::sqrt(dof * PI);
        return t1 / t2 * std::pow(1.0 + x * x / dof, -(dof + 1.0) / 2.0);
    }
};

}  // namespace v1
}  // namespace cpphub
