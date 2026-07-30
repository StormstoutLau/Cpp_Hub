#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include <cmath>
#include <limits>
#include <vector>
#include <functional>

namespace cpphub {
inline namespace v1 {

inline Real bsm_call_price(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(S - K, 0.0);
    Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    Real d2 = d1 - sigma * std::sqrt(T);
    return S * std::exp(-q * T) * normal_cdf(d1) - K * std::exp(-r * T) * normal_cdf(d2);
}

inline Real bsm_put_price(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(K - S, 0.0);
    Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    Real d2 = d1 - sigma * std::sqrt(T);
    return K * std::exp(-r * T) * normal_cdf(-d2) - S * std::exp(-q * T) * normal_cdf(-d1);
}

class ControlVariate {
public:
    explicit ControlVariate(Real control_expectation)
        : control_expectation_(control_expectation) {}

    void add_sample(Real option_payoff, Real control_value) {
        ++n_;
        sum_y_ += option_payoff;
        sum_x_ += control_value;
        sum_yy_ += option_payoff * option_payoff;
        sum_xx_ += control_value * control_value;
        sum_xy_ += option_payoff * control_value;
    }

    Real estimate() const {
        if (n_ == 0) return 0.0;
        Real inv_n = 1.0 / static_cast<Real>(n_);
        Real mean_y = sum_y_ * inv_n;
        Real mean_x = sum_x_ * inv_n;
        if (n_ < 2) return mean_y;
        Real cov = sum_xy_ * inv_n - mean_x * mean_y;
        Real var_x = sum_xx_ * inv_n - mean_x * mean_x;
        if (var_x <= 0.0) return mean_y;
        Real beta = cov / var_x;
        return mean_y - beta * (mean_x - control_expectation_);
    }

    Real variance_reduction_ratio() const {
        if (n_ < 2) return 1.0;
        Real inv_n = 1.0 / static_cast<Real>(n_);
        Real mean_y = sum_y_ * inv_n;
        Real mean_x = sum_x_ * inv_n;
        Real var_y = sum_yy_ * inv_n - mean_y * mean_y;
        Real var_x = sum_xx_ * inv_n - mean_x * mean_x;
        if (var_y <= 0.0 || var_x <= 0.0) return 1.0;
        Real cov = sum_xy_ * inv_n - mean_x * mean_y;
        Real rho = cov / std::sqrt(var_y * var_x);
        Real rho2 = rho * rho;
        if (rho2 >= 1.0) return std::numeric_limits<Real>::max();
        return 1.0 / (1.0 - rho2);
    }

    void reset() {
        n_ = 0;
        sum_y_ = sum_x_ = 0.0;
        sum_yy_ = sum_xx_ = sum_xy_ = 0.0;
    }

    Size sample_count() const noexcept { return n_; }

private:
    Real control_expectation_;
    Size n_{0};
    Real sum_y_{0}, sum_x_{0};
    Real sum_yy_{0}, sum_xx_{0}, sum_xy_{0};
};

class BSControlVariate {
public:
    BSControlVariate(Real S0, Real K, Real T, Real r, Real q, Real sigma, bool is_call)
        : bs_price_(is_call ? bsm_call_price(S0, K, T, r, q, sigma)
                            : bsm_put_price(S0, K, T, r, q, sigma))
        , K_(K), is_call_(is_call)
        , cv_(bs_price_ * std::exp(r * T)) {}

    void add_path(Real terminal_spot, Real option_payoff) {
        Real control_value = is_call_
            ? std::max(terminal_spot - K_, 0.0)
            : std::max(K_ - terminal_spot, 0.0);
        cv_.add_sample(option_payoff, control_value);
    }

    Real estimate() const { return cv_.estimate(); }

    Real variance_reduction_ratio() const { return cv_.variance_reduction_ratio(); }

    Size sample_count() const noexcept { return cv_.sample_count(); }

private:
    Real bs_price_;
    Real K_;
    bool is_call_;
    ControlVariate cv_;
};

}  // namespace v1
}  // namespace cpphub
