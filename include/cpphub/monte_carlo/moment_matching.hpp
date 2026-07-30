#pragma once
#include "cpphub/core/types.hpp"
#include <vector>
#include <functional>
#include <cmath>

namespace cpphub {
inline namespace v1 {

class MomentMatchingVR {
public:
    MomentMatchingVR(Real theoretical_mean, Real theoretical_variance)
        : mu_(theoretical_mean), var_(theoretical_variance) {}

    void add_sample(Real sample) {
        ++n_;
        sum_ += sample;
        sum_sq_ += sample * sample;
    }

    Real estimate() const {
        return mu_;
    }

    Real variance_reduction_ratio() const {
        if (n_ < 2 || var_ <= 0.0) return 1.0;
        Real inv_n = 1.0 / static_cast<Real>(n_);
        Real sample_mean = sum_ * inv_n;
        Real sample_var = sum_sq_ * inv_n - sample_mean * sample_mean;
        if (sample_var <= 0.0) return 1.0;
        return sample_var / var_;
    }

    void reset() {
        n_ = 0;
        sum_ = sum_sq_ = 0.0;
    }

    Size sample_count() const noexcept { return n_; }

private:
    Real mu_, var_;
    Size n_{0};
    Real sum_{0}, sum_sq_{0};
};

class GBMMomentMatching {
public:
    GBMMomentMatching(Real S0, Real T, Real r, Real q, Real sigma)
        : S0_(S0), T_(T), r_(r), q_(q), sigma_(sigma) {}

    void add_path(Real terminal_spot) {
        samples_.push_back(terminal_spot);
        mm_.add_sample(terminal_spot);
    }

    Real estimate_mean() const {
        Real mean = S0_ * std::exp((r_ - q_) * T_);
        return mean;
    }

    Real estimate_option_price(const std::function<Real(Real)>& payoff) const {
        if (samples_.empty()) return 0.0;
        Size n = samples_.size();
        Real inv_n = 1.0 / static_cast<Real>(n);

        Real sample_sum = 0.0, sample_sq = 0.0;
        for (auto s : samples_) {
            sample_sum += s;
            sample_sq += s * s;
        }
        Real sample_mean = sample_sum * inv_n;
        Real sample_var = sample_sq * inv_n - sample_mean * sample_mean;
        if (sample_var <= 0.0) {
            Real sum_p = 0.0;
            for (auto s : samples_) sum_p += payoff(s);
            return sum_p * inv_n;
        }
        Real sample_std = std::sqrt(sample_var);

        Real theor_mean = S0_ * std::exp((r_ - q_) * T_);
        Real theor_var = theor_mean * theor_mean * (std::exp(sigma_ * sigma_ * T_) - 1.0);
        if (theor_var <= 0.0) {
            Real sum_p = 0.0;
            for (auto s : samples_) sum_p += payoff(s);
            return sum_p * inv_n;
        }
        Real theor_std = std::sqrt(theor_var);

        Real sum_payoff = 0.0;
        for (auto s : samples_) {
            Real s_corrected = (s - sample_mean) * (theor_std / sample_std) + theor_mean;
            sum_payoff += payoff(s_corrected);
        }
        return sum_payoff * inv_n;
    }

    Size sample_count() const noexcept { return samples_.size(); }

private:
    MomentMatchingVR mm_{0.0, 0.0};
    std::vector<Real> samples_;
    Real S0_, T_, r_, q_, sigma_;
};

}  // namespace v1
}  // namespace cpphub
