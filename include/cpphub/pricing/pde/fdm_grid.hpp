#pragma once
#include <vector>
#include <cmath>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {

class FDMGrid {
public:
    FDMGrid(Size n_points, Real S0, Real K, Real sigma, Real T, Real alpha = 0.2)
        : S_(n_points) {
        Real c = alpha;
        Real s_min = S0 * std::exp(-c * sigma * std::sqrt(T) * 5.0);
        Real s_max = S0 * std::exp( c * sigma * std::sqrt(T) * 5.0);
        build_sinh_grid(n_points, s_min, s_max, S0, c);
    }

    Size size() const noexcept { return S_.size(); }
    Real s(Size i) const noexcept { return S_[i]; }
    Real ds(Size i) const noexcept { return S_[i + 1] - S_[i]; }
    Real ds_avg(Size i) const noexcept { return Real(0.5) * (S_[i + 1] - S_[i - 1]); }

    Real s_min() const noexcept { return S_.front(); }
    Real s_max() const noexcept { return S_.back(); }

    const std::vector<Real>& points() const noexcept { return S_; }

private:
    std::vector<Real> S_;

    void build_sinh_grid(Size n, Real S_min, Real S_max, Real S0, Real alpha) {
        Real center = (S_min + S_max) / 2.0;
        Real range = (S_max - S_min) / 2.0;
        Real inv_sinh_alpha = 1.0 / std::sinh(alpha);
        for (Size i = 0; i < n; ++i) {
            Real xi = -1.0 + 2.0 * static_cast<Real>(i) / static_cast<Real>(n - 1);
            S_[i] = center + range * std::sinh(alpha * xi) * inv_sinh_alpha;
        }
    }
};

class TimeGrid {
public:
    TimeGrid(Size n_steps, Real T)
        : t_(n_steps + 1), dt_(T / static_cast<Real>(n_steps)) {
        for (Size i = 0; i <= n_steps; ++i) {
            t_[i] = static_cast<Real>(i) * dt_;
        }
    }

    Size size() const noexcept { return t_.size(); }
    Real t(Size i) const noexcept { return t_[i]; }
    Real dt() const noexcept { return dt_; }
    Real T() const noexcept { return t_.back(); }

private:
    std::vector<Real> t_;
    Real dt_;
};

}  // namespace v1
}  // namespace cpphub
