#pragma once
#include "cpphub/core/types.hpp"
#include <vector>
#include <map>
#include <string>
#include <functional>
#include <cmath>
#include <random>
#include <algorithm>
#include <limits>

namespace cpphub {
inline namespace v1 {

struct SensitivityResult {
    std::string factor_name;
    Real shock;
    Real pnl;
    Real pct_change;
};

class SensitivityAnalysis {
public:
    using ValueFn = std::function<Real(const std::map<std::string, Real>&)>;

    SensitivityAnalysis(ValueFn value_fn, std::map<std::string, Real> current_factors)
        : value_fn_(std::move(value_fn)), current_factors_(std::move(current_factors)) {
        base_value_ = value_fn_(current_factors_);
    }

    std::vector<SensitivityResult> single_factor(
        const std::vector<std::string>& factors,
        Real shock_pct = Real(0.01)) {
        std::vector<SensitivityResult> results;
        results.reserve(factors.size());

        for (const auto& factor : factors) {
            auto shocked = current_factors_;
            if (shocked.find(factor) != shocked.end()) {
                shocked[factor] = current_factors_.at(factor) * (Real(1) + shock_pct);
            } else {
                shocked[factor] = shock_pct;
            }
            Real new_value = value_fn_(shocked);
            Real pnl = new_value - base_value_;
            Real pct_change = (std::abs(base_value_) > std::numeric_limits<Real>::epsilon())
                ? pnl / base_value_
                : Real(0);

            results.push_back({factor, shock_pct, pnl, pct_change});
        }

        return results;
    }

    std::vector<std::vector<SensitivityResult>> multi_factor(
        const std::vector<std::string>& factors,
        const std::vector<Real>& shock_levels) {
        std::vector<std::vector<SensitivityResult>> results;
        results.resize(factors.size());

        for (Size i = 0; i < factors.size(); ++i) {
            results[i].reserve(shock_levels.size());
            for (Real level : shock_levels) {
                auto shocked = current_factors_;
                if (shocked.find(factors[i]) != shocked.end()) {
                    shocked[factors[i]] = current_factors_.at(factors[i]) * (Real(1) + level);
                } else {
                    shocked[factors[i]] = level;
                }
                Real new_value = value_fn_(shocked);
                Real pnl = new_value - base_value_;
                Real pct_change = (std::abs(base_value_) > std::numeric_limits<Real>::epsilon())
                    ? pnl / base_value_
                    : Real(0);

                results[i].push_back({factors[i], level, pnl, pct_change});
            }
        }

        return results;
    }

    std::vector<SensitivityResult> directional(
        const std::map<std::string, Real>& direction,
        const std::vector<Real>& magnitudes) {
        std::vector<SensitivityResult> results;
        results.reserve(magnitudes.size());

        for (Real mag : magnitudes) {
            auto shocked = current_factors_;
            Real total_pnl = Real(0);
            for (const auto& [factor, dir] : direction) {
                Real shock = dir * mag;
                if (shocked.find(factor) != shocked.end()) {
                    shocked[factor] = current_factors_.at(factor) * (Real(1) + shock);
                } else {
                    shocked[factor] = shock;
                }
            }
            Real new_value = value_fn_(shocked);
            Real pnl = new_value - base_value_;
            Real pct_change = (std::abs(base_value_) > std::numeric_limits<Real>::epsilon())
                ? pnl / base_value_
                : Real(0);

            std::string desc;
            for (const auto& [f, d] : direction) {
                if (!desc.empty()) desc += "+";
                desc += f;
            }
            results.push_back({desc, mag, pnl, pct_change});
        }

        return results;
    }

    Real worst_case(const std::map<std::string, std::pair<Real, Real>>& factor_ranges,
                    Size n_samples = 10000) {
        std::mt19937_64 rng(42);
        Real worst_pnl = Real(0);

        std::vector<std::string> factor_names;
        std::vector<std::uniform_real_distribution<Real>> dists;
        for (const auto& [name, range] : factor_ranges) {
            factor_names.push_back(name);
            dists.emplace_back(range.first, range.second);
        }

        Size nf = factor_names.size();

        for (Size s = 0; s < n_samples; ++s) {
            auto shocked = current_factors_;
            for (Size i = 0; i < nf; ++i) {
                Real shock_val = dists[i](rng);
                if (shocked.find(factor_names[i]) != shocked.end()) {
                    shocked[factor_names[i]] = current_factors_.at(factor_names[i]) * (Real(1) + shock_val);
                } else {
                    shocked[factor_names[i]] = shock_val;
                }
            }
            Real new_value = value_fn_(shocked);
            Real pnl = new_value - base_value_;
            if (pnl < worst_pnl) {
                worst_pnl = pnl;
            }
        }

        return worst_pnl;
    }

private:
    ValueFn value_fn_;
    std::map<std::string, Real> current_factors_;
    Real base_value_;
};

} // namespace v1
} // namespace cpphub
