#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

namespace detail {

inline Real chi2_cdf(Real x, Size k) {
    if (x <= 0) return 0;
    Real kf = static_cast<Real>(k);
    Real t = std::pow(x / kf, Real(1) / Real(3));
    Real mean = Real(1) - Real(2) / (Real(9) * kf);
    Real var = Real(2) / (Real(9) * kf);
    Real z = (t - mean) / std::sqrt(var);
    return normal_cdf(z);
}

} // namespace detail

struct BacktestResult {
    Size n_violations;
    Size n_observations;
    Real violation_rate;
    Real expected_rate;
    Real p_value;
    bool reject_null;
};

class KupiecPOF {
public:
    static BacktestResult test(Size n_violations, Size n_observations, Real confidence) {
        Real p = Real(1) - confidence;
        Real x = static_cast<Real>(n_violations);
        Real N = static_cast<Real>(n_observations);
        Real violation_rate = x / N;

        if (n_violations == 0) {
            Real L0 = std::pow(Real(1) - p, N);
            Real L1 = Real(1);
            Real LR = Real(-2) * std::log(L0 / L1);
            Real p_value = Real(1) - detail::chi2_cdf(LR, Size(1));
            return {n_violations, n_observations, violation_rate, p, p_value, p_value < Real(0.05)};
        }
        if (n_violations == n_observations) {
            Real L0 = std::pow(p, N);
            Real L1 = Real(1);
            Real LR = Real(-2) * std::log(L0 / L1);
            Real p_value = Real(1) - detail::chi2_cdf(LR, Size(1));
            return {n_violations, n_observations, violation_rate, p, p_value, p_value < Real(0.05)};
        }

        Real L0 = std::pow(Real(1) - p, N - x) * std::pow(p, x);
        Real L1 = std::pow(Real(1) - violation_rate, N - x) * std::pow(violation_rate, x);
        Real LR = Real(-2) * std::log(L0 / L1);
        Real p_value = Real(1) - detail::chi2_cdf(LR, Size(1));

        return {n_violations, n_observations, violation_rate, p, p_value, p_value < Real(0.05)};
    }

    static BacktestResult test(const std::vector<Real>& var_series,
                                const std::vector<Real>& realized_losses,
                                Real confidence) {
        Size n_violations = 0;
        Size n = std::min(var_series.size(), realized_losses.size());
        for (Size i = 0; i < n; ++i) {
            if (realized_losses[i] > var_series[i]) {
                ++n_violations;
            }
        }
        return test(n_violations, n, confidence);
    }
};

class ChristoffersenIID {
public:
    static BacktestResult test(const std::vector<Real>& var_series,
                                const std::vector<Real>& realized_losses) {
        Size n = std::min(var_series.size(), realized_losses.size());
        std::vector<int> breaches(n, 0);
        for (Size i = 0; i < n; ++i) {
            if (realized_losses[i] > var_series[i]) {
                breaches[i] = 1;
            }
        }

        Size n00 = 0, n01 = 0, n10 = 0, n11 = 0;
        for (Size i = 1; i < n; ++i) {
            if (breaches[i-1] == 0 && breaches[i] == 0) ++n00;
            if (breaches[i-1] == 0 && breaches[i] == 1) ++n01;
            if (breaches[i-1] == 1 && breaches[i] == 0) ++n10;
            if (breaches[i-1] == 1 && breaches[i] == 1) ++n11;
        }

        Size n0 = n00 + n01;
        Size n1 = n10 + n11;

        Real pi = static_cast<Real>(n01 + n11) / static_cast<Real>(n - 1);
        Real pi0 = (n0 > 0) ? static_cast<Real>(n01) / static_cast<Real>(n0) : Real(0);
        Real pi1 = (n1 > 0) ? static_cast<Real>(n11) / static_cast<Real>(n1) : Real(0);

        Real L_IID = std::pow(Real(1) - pi, static_cast<Real>(n - 1 - n01 - n11)) * std::pow(pi, static_cast<Real>(n01 + n11));
        Real L_nonIID = std::pow(Real(1) - pi0, static_cast<Real>(n00)) * std::pow(pi0, static_cast<Real>(n01))
                      * std::pow(Real(1) - pi1, static_cast<Real>(n10)) * std::pow(pi1, static_cast<Real>(n11));

        Real LR_ind = Real(0);
        if (L_IID > 0 && L_nonIID > 0) {
            LR_ind = Real(-2) * std::log(L_IID / L_nonIID);
        }

        Real p_value = Real(1) - detail::chi2_cdf(LR_ind, Size(1));
        Size v = n01 + n11;

        return {v, n, static_cast<Real>(v) / static_cast<Real>(n), pi, p_value, p_value < Real(0.05)};
    }

    static BacktestResult joint_test(const std::vector<Real>& var_series,
                                      const std::vector<Real>& realized_losses,
                                      Real confidence) {
        Size n = std::min(var_series.size(), realized_losses.size());
        Real p = Real(1) - confidence;
        Size v = 0;
        for (Size i = 0; i < n; ++i) {
            if (realized_losses[i] > var_series[i]) ++v;
        }

        // LR_cc = LR_pof + LR_ind  (joint test)
        auto pof_result = KupiecPOF::test(v, n, confidence);
        auto ind_result = test(var_series, realized_losses);

        Real LR_cc = Real(-2) * std::log(
            std::pow(Real(1) - p, static_cast<Real>(n - v)) * std::pow(p, static_cast<Real>(v))
        );

        if (v > 0 && v < n) {
            Real vrate = static_cast<Real>(v) / static_cast<Real>(n);
            LR_cc = Real(-2) * std::log(
                std::pow(Real(1) - p, static_cast<Real>(n - v)) * std::pow(p, static_cast<Real>(v))
                / (std::pow(Real(1) - vrate, static_cast<Real>(n - v)) * std::pow(vrate, static_cast<Real>(v)))
            );
        }

        Real ind_LR = Real(-2) * std::log(
            std::pow(Real(1) - ind_result.violation_rate, static_cast<Real>(n - 1 - (v - (n > 0 && realized_losses[0] > var_series[0] ? 1 : 0)))) 
            * std::pow(ind_result.violation_rate, static_cast<Real>(v - (n > 0 && realized_losses[0] > var_series[0] ? 1 : 0)))
        );

        Real LR_cc_total = LR_cc;

        Real p_value = Real(1) - detail::chi2_cdf(LR_cc_total, Size(2));

        return {v, n, static_cast<Real>(v) / static_cast<Real>(n), p, p_value, p_value < Real(0.05)};
    }
};

enum class BaselZone { Green, Yellow, Red };

struct BaselTrafficLightResult {
    BaselZone zone;
    Size n_violations;
    Size n_observations;
    Real capital_multiplier;
    std::string description;
};

class BaselTrafficLight {
public:
    static BaselTrafficLightResult assess(Size n_violations, Size n_observations = 250) {
        BaselZone zone;
        Real multiplier;
        std::string desc;

        if (n_violations <= 4) {
            zone = BaselZone::Green;
            multiplier = Real(3.0);
            desc = "Green zone: model OK";
        } else if (n_violations <= 9) {
            zone = BaselZone::Yellow;
            static const Real multipliers[] = {Real(3.4), Real(3.5), Real(3.65), Real(3.75), Real(3.85)};
            multiplier = multipliers[n_violations - 5];
            desc = "Yellow zone: model may be inaccurate";
        } else {
            zone = BaselZone::Red;
            multiplier = Real(4.0);
            desc = "Red zone: model is inaccurate";
        }

        return {zone, n_violations, n_observations, multiplier, desc};
    }
};

} // namespace v1
} // namespace cpphub
