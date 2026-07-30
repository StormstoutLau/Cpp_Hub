#pragma once
#include "cpphub/core/types.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>

namespace cpphub {
inline namespace v1 {

enum class QuantileInterpolation { Linear, Conservative, Empirical };

class HistoricalVaR {
public:
    HistoricalVaR(const std::vector<Real>& pnl_history, Real confidence = 0.99,
                  Size horizon_days = 1)
        : confidence_(confidence), horizon_days_(horizon_days)
    {
        losses_.reserve(pnl_history.size());
        for (auto pnl : pnl_history) {
            losses_.push_back(-pnl);
        }
    }

    Real var(QuantileInterpolation interp = QuantileInterpolation::Linear) const {
        if (losses_.empty()) return 0;
        std::vector<Real> sorted = losses_;
        std::sort(sorted.begin(), sorted.end());
        Real q = 1.0 - confidence_;
        Real val = quantile_interpolate(sorted, q, interp);
        return -val * std::sqrt(static_cast<Real>(horizon_days_));
    }

    static std::vector<Real> rolling_var(const std::vector<Real>& pnl_history,
                                          Size window_size, Real confidence,
                                          QuantileInterpolation interp = QuantileInterpolation::Linear) {
        std::vector<Real> result;
        if (pnl_history.size() <= window_size) return result;
        result.reserve(pnl_history.size() - window_size);
        for (Size t = window_size; t < pnl_history.size(); ++t) {
            HistoricalVaR hv(std::vector<Real>(pnl_history.begin() + t - window_size, pnl_history.begin() + t),
                           confidence, 1);
            result.push_back(hv.var(interp));
        }
        return result;
    }

    Real weighted_var(Real decay = 0.99,
                      QuantileInterpolation interp = QuantileInterpolation::Linear) const {
        if (losses_.empty()) return 0;
        Size n = losses_.size();
        std::vector<std::pair<Real, Real>> weighted;
        weighted.reserve(n);
        Real weight_sum = 0;
        for (Size i = 0; i < n; ++i) {
            Real w = (1.0 - decay) * std::pow(decay, static_cast<Real>(n - 1 - i));
            weighted.push_back({losses_[i], w});
            weight_sum += w;
        }
        std::sort(weighted.begin(), weighted.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        Real cum_weight = 0;
        Real q = 1.0 - confidence_;
        Real target = q * weight_sum;
        for (Size i = 0; i < n - 1; ++i) {
            Real next_cum = cum_weight + weighted[i].second;
            if (next_cum >= target) {
                if (interp == QuantileInterpolation::Conservative) {
                    return -weighted[i].first * std::sqrt(static_cast<Real>(horizon_days_));
                } else if (interp == QuantileInterpolation::Linear) {
                    Real frac = (target - cum_weight) / weighted[i].second;
                    return -(weighted[i].first + frac * (weighted[i + 1].first - weighted[i].first))
                           * std::sqrt(static_cast<Real>(horizon_days_));
                } else {
                    return -weighted[i].first * std::sqrt(static_cast<Real>(horizon_days_));
                }
            }
            cum_weight = next_cum;
        }
        return -weighted.back().first * std::sqrt(static_cast<Real>(horizon_days_));
    }

    std::pair<Real, Real> bootstrap_ci(Size n_bootstrap = 10000, Real ci_level = 0.95) const {
        if (losses_.empty()) return {0, 0};
        std::vector<Real> bootstrap_vars;
        bootstrap_vars.reserve(n_bootstrap);
        Size n = losses_.size();
        std::mt19937_64 gen(12345);
        std::uniform_int_distribution<Size> dist(0, n - 1);
        for (Size b = 0; b < n_bootstrap; ++b) {
            std::vector<Real> sample;
            sample.reserve(n);
            for (Size i = 0; i < n; ++i) {
                sample.push_back(losses_[dist(gen)]);
            }
            HistoricalVaR hv(sample, confidence_, 1);
            bootstrap_vars.push_back(hv.var(QuantileInterpolation::Linear));
        }
        std::sort(bootstrap_vars.begin(), bootstrap_vars.end());
        Real lower_q = (1.0 - ci_level) / 2.0;
        Real upper_q = 1.0 - lower_q;
        Size lower_idx = static_cast<Size>(lower_q * n_bootstrap);
        Size upper_idx = static_cast<Size>(upper_q * n_bootstrap);
        if (lower_idx >= n_bootstrap) lower_idx = n_bootstrap - 1;
        if (upper_idx >= n_bootstrap) upper_idx = n_bootstrap - 1;
        return {bootstrap_vars[lower_idx], bootstrap_vars[upper_idx]};
    }

private:
    std::vector<Real> losses_;
    Real confidence_;
    Size horizon_days_;

    Real quantile_interpolate(const std::vector<Real>& sorted_losses, Real q,
                              QuantileInterpolation interp) const {
        Size n = sorted_losses.size();
        if (n == 0) return 0;
        if (n == 1) return sorted_losses[0];
        if (q <= 0) return sorted_losses[0];
        if (q >= 1) return sorted_losses[n - 1];

        if (interp == QuantileInterpolation::Conservative) {
            Real index = q * (n - 1);
            Size lo = static_cast<Size>(std::floor(index));
            Size hi = static_cast<Size>(std::ceil(index));
            if (lo >= n) lo = n - 1;
            if (hi >= n) hi = n - 1;
            if (lo == hi) return sorted_losses[lo];
            if (std::abs(sorted_losses[lo]) >= std::abs(sorted_losses[hi])) {
                return sorted_losses[lo];
            } else {
                return sorted_losses[hi];
            }
        } else if (interp == QuantileInterpolation::Linear) {
            Real index = q * (n - 1);
            Size lo = static_cast<Size>(std::floor(index));
            Size hi = static_cast<Size>(std::ceil(index));
            if (lo == hi) return sorted_losses[lo];
            Real frac = index - lo;
            return sorted_losses[lo] + frac * (sorted_losses[hi] - sorted_losses[lo]);
        } else {
            Size idx = static_cast<Size>(q * n);
            if (idx >= n) idx = n - 1;
            return sorted_losses[idx];
        }
    }
};

}  // namespace v1
}  // namespace cpphub
