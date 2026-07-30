#pragma once
// SOURCE: PHASE3_SPEC §4.2 - Volatility surface (interpolation + SVI parameterization)
// Implemented on main station (MSVC) - 2026-07-31
// 2D implied vol grid with bilinear / cubic spline interpolation,
// SVI slice fitting, calendar/butterfly arbitrage checks.
// For Dupire local volatility, use DupireLocalVol (dupire_local_vol.hpp)
//   which has a self-contained IV grid path and does not depend on this class.
#include "cpphub/core/types.hpp"
#include "cpphub/models/vol_surface/svi.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"  // bsm_call_price / bsm_put_price
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

enum class InterpolationMethod {
    Bilinear,
    Bicubic,
    CubicSpline,
    SVI
};

class VolSurface {
public:
    VolSurface();
    VolSurface(const std::vector<Real>& strikes,
                const std::vector<Real>& maturities,
                const std::vector<std::vector<Real>>& implied_vols);

    // Set market context (forward, rate, dividend) for BSM inversion
    void set_market(Real forward, Real r, Real q) {
        forward_ = forward; r_ = r; q_ = q;
    }

    Real implied_vol(Real K, Real T,
                     InterpolationMethod method = InterpolationMethod::Bilinear) const;
    Real total_variance(Real K, Real T) const;
    Real log_moneyness(Real K, Real T) const;

    // Derivatives (for Dupire): finite-difference on the IV grid
    Real dC_dT(Real K, Real T) const;
    Real dC_dK(Real K, Real T) const;
    Real d2C_dK2(Real K, Real T) const;

    // BSM prices using this surface's IV
    Real call_price(Real K, Real T, Real S, Real r, Real q) const;
    Real put_price(Real K, Real T, Real S, Real r, Real q) const;

    // Data access
    const std::vector<Real>& strikes() const { return strikes_; }
    const std::vector<Real>& maturities() const { return maturities_; }
    const std::vector<std::vector<Real>>& vols() const { return implied_vols_; }

    // SVI fit at a given maturity slice
    void fit_svi(Real T);
    const SVI* svi_slice(Real T) const;

    // Arbitrage checks
    bool check_calendar_arbitrage() const;
    bool check_butterfly_arbitrage() const;

private:
    std::vector<Real> strikes_;
    std::vector<Real> maturities_;
    std::vector<std::vector<Real>> implied_vols_;  // [T_idx][K_idx]
    Real forward_ = 100.0;
    Real r_ = 0.0;
    Real q_ = 0.0;

    // Cached SVI slices keyed by maturity
    mutable std::map<Real, std::unique_ptr<SVI>> svi_slices_;
    mutable bool splines_built_ = false;
    // Cubic spline second derivatives per maturity row
    mutable std::vector<std::vector<Real>> spline_d2y_;

    Real bilinear_interp(Real K, Real T) const;
    Real cubic_spline_interp(Real K, Real T) const;
    void build_splines() const;
    // Find index i such that strikes_[i] <= K < strikes_[i+1]
    Size find_strike_index(Real K) const;
    Size find_maturity_index(Real T) const;
};

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

inline VolSurface::VolSurface() = default;

inline VolSurface::VolSurface(
        const std::vector<Real>& strikes,
        const std::vector<Real>& maturities,
        const std::vector<std::vector<Real>>& implied_vols)
    : strikes_(strikes), maturities_(maturities), implied_vols_(implied_vols) {
    if (strikes_.size() < 2) {
        throw std::invalid_argument("VolSurface: need at least 2 strikes");
    }
    if (maturities_.size() < 1) {
        throw std::invalid_argument("VolSurface: need at least 1 maturity");
    }
    if (implied_vols_.size() != maturities_.size()) {
        throw std::invalid_argument("VolSurface: implied_vols rows != maturities");
    }
    for (Size j = 0; j < maturities_.size(); ++j) {
        if (implied_vols_[j].size() != strikes_.size()) {
            throw std::invalid_argument("VolSurface: implied_vols cols != strikes");
        }
    }
    // Verify strikes sorted ascending
    for (Size i = 1; i < strikes_.size(); ++i) {
        if (strikes_[i] <= strikes_[i - 1]) {
            throw std::invalid_argument("VolSurface: strikes must be strictly ascending");
        }
    }
    for (Size j = 1; j < maturities_.size(); ++j) {
        if (maturities_[j] <= maturities_[j - 1]) {
            throw std::invalid_argument("VolSurface: maturities must be strictly ascending");
        }
    }
}

inline Size VolSurface::find_strike_index(Real K) const {
    if (K <= strikes_.front()) return 0;
    if (K >= strikes_.back()) return strikes_.size() - 2;
    Size lo = 0, hi = strikes_.size() - 1;
    while (hi - lo > 1) {
        Size mid = (lo + hi) / 2;
        if (strikes_[mid] <= K) lo = mid; else hi = mid;
    }
    return lo;
}

inline Size VolSurface::find_maturity_index(Real T) const {
    if (T <= maturities_.front()) return 0;
    if (T >= maturities_.back()) return maturities_.size() - 2;
    Size lo = 0, hi = maturities_.size() - 1;
    while (hi - lo > 1) {
        Size mid = (lo + hi) / 2;
        if (maturities_[mid] <= T) lo = mid; else hi = mid;
    }
    return lo;
}

inline Real VolSurface::bilinear_interp(Real K, Real T) const {
    // Clamp T to grid
    Size j = find_maturity_index(T);
    Size j1 = std::min(j + 1, maturities_.size() - 1);
    Real t0 = maturities_[j], t1 = maturities_[j1];
    Real wt = (t1 > t0) ? (T - t0) / (t1 - t0) : 0.0;
    wt = std::max(0.0, std::min(1.0, wt));

    // Clamp K to grid
    Size i = find_strike_index(K);
    Size i1 = std::min(i + 1, strikes_.size() - 1);
    Real k0 = strikes_[i], k1 = strikes_[i1];
    Real wk = (k1 > k0) ? (K - k0) / (k1 - k0) : 0.0;
    wk = std::max(0.0, std::min(1.0, wk));

    Real v00 = implied_vols_[j][i];
    Real v01 = implied_vols_[j][i1];
    Real v10 = implied_vols_[j1][i];
    Real v11 = implied_vols_[j1][i1];
    Real v0 = v00 * (1.0 - wk) + v01 * wk;
    Real v1 = v10 * (1.0 - wk) + v11 * wk;
    return v0 * (1.0 - wt) + v1 * wt;
}

// Natural cubic spline (Press et al. Numerical Recipes)
inline void VolSurface::build_splines() const {
    if (splines_built_) return;
    Size n = strikes_.size();
    Size nT = maturities_.size();
    spline_d2y_.resize(nT, std::vector<Real>(n, 0.0));
    for (Size j = 0; j < nT; ++j) {
        std::vector<Real>& y2 = spline_d2y_[j];
        const std::vector<Real>& y = implied_vols_[j];
        std::vector<Real> u(n, 0.0);
        y2[0] = 0.0;  // natural spline
        u[0] = 0.0;
        for (Size i = 1; i < n - 1; ++i) {
            Real h0 = strikes_[i] - strikes_[i - 1];
            Real h1 = strikes_[i + 1] - strikes_[i];
            Real sig = h0 / (h0 + h1);
            Real p = sig * y2[i - 1] + 2.0;
            y2[i] = (sig - 1.0) / p;
            u[i] = (6.0 * ((y[i + 1] - y[i]) / h1 - (y[i] - y[i - 1]) / h0) /
                    (h0 + h1) - sig * u[i - 1]) / p;
        }
        y2[n - 1] = 0.0;  // natural spline
        for (Size i = n - 1; i-- > 0;) {
            y2[i] = y2[i] * y2[i + 1] + u[i];
        }
    }
    splines_built_ = true;
}

inline Real VolSurface::cubic_spline_interp(Real K, Real T) const {
    build_splines();
    // Interpolate in K using cubic spline at the two surrounding maturities,
    // then linear in T (bilinear in T for stability).
    Size j = find_maturity_index(T);
    Size j1 = std::min(j + 1, maturities_.size() - 1);
    Real t0 = maturities_[j], t1 = maturities_[j1];
    Real wt = (t1 > t0) ? (T - t0) / (t1 - t0) : 0.0;
    wt = std::max(0.0, std::min(1.0, wt));

    auto spline_at = [&](Size jj, Real kk) -> Real {
        Size i = find_strike_index(kk);
        Size i1 = std::min(i + 1, strikes_.size() - 1);
        Real h = strikes_[i1] - strikes_[i];
        if (h <= 0.0) return implied_vols_[jj][i];
        Real a = (strikes_[i1] - kk) / h;
        Real b = (kk - strikes_[i]) / h;
        const std::vector<Real>& y = implied_vols_[jj];
        const std::vector<Real>& y2 = spline_d2y_[jj];
        return a * y[i] + b * y[i1] +
               ((a * a * a - a) * y2[i] + (b * b * b - b) * y2[i1]) * h * h / 6.0;
    };

    Real v0 = spline_at(j, K);
    Real v1 = spline_at(j1, K);
    return v0 * (1.0 - wt) + v1 * wt;
}

inline Real VolSurface::implied_vol(Real K, Real T, InterpolationMethod method) const {
    if (T <= maturities_.front()) {
        // Extrapolate flatly in T at shortest maturity
        T = maturities_.front();
    }
    if (T >= maturities_.back()) {
        T = maturities_.back();
    }
    switch (method) {
        case InterpolationMethod::Bilinear:
        case InterpolationMethod::Bicubic:
            return bilinear_interp(K, T);
        case InterpolationMethod::CubicSpline:
            return cubic_spline_interp(K, T);
        case InterpolationMethod::SVI: {
            // Find nearest cached SVI slice
            auto it = svi_slices_.lower_bound(T);
            if (it == svi_slices_.end() || it->first > T) {
                if (it != svi_slices_.begin()) --it;
            }
            if (it == svi_slices_.end()) {
                // No SVI fit; fall back to cubic spline
                return cubic_spline_interp(K, T);
            }
            Real k = std::log(K / forward_);
            return it->second->implied_vol(k, T);
        }
    }
    return bilinear_interp(K, T);
}

inline Real VolSurface::total_variance(Real K, Real T) const {
    Real iv = implied_vol(K, T);
    return iv * iv * T;
}

inline Real VolSurface::log_moneyness(Real K, Real T) const {
    (void)T;
    return std::log(K / forward_);
}

inline Real VolSurface::call_price(Real K, Real T, Real S, Real r, Real q) const {
    Real sigma = implied_vol(K, T);
    return bsm_call_price(S, K, T, r, q, sigma);
}

inline Real VolSurface::put_price(Real K, Real T, Real S, Real r, Real q) const {
    Real sigma = implied_vol(K, T);
    return bsm_put_price(S, K, T, r, q, sigma);
}

// Finite-difference derivatives on the IV grid (central differences)
inline Real VolSurface::dC_dT(Real K, Real T) const {
    // Use central difference in T with BSM call prices
    Real dT = std::max(1e-4, 1e-4 * T);
    Real C_plus = call_price(K, T + dT, forward_, r_, q_);
    Real C_minus = call_price(K, T - dT, forward_, r_, q_);
    return (C_plus - C_minus) / (2.0 * dT);
}

inline Real VolSurface::dC_dK(Real K, Real T) const {
    Real dK = std::max(1e-4, 1e-4 * K);
    Real C_plus = call_price(K + dK, T, forward_, r_, q_);
    Real C_minus = call_price(K - dK, T, forward_, r_, q_);
    return (C_plus - C_minus) / (2.0 * dK);
}

inline Real VolSurface::d2C_dK2(Real K, Real T) const {
    Real dK = std::max(1e-3, 1e-3 * K);
    Real C_plus = call_price(K + dK, T, forward_, r_, q_);
    Real C_mid = call_price(K, T, forward_, r_, q_);
    Real C_minus = call_price(K - dK, T, forward_, r_, q_);
    return (C_plus - 2.0 * C_mid + C_minus) / (dK * dK);
}

inline void VolSurface::fit_svi(Real T) {
    // Find the maturity row closest to T
    Size j = 0;
    Real best_dist = std::abs(maturities_[0] - T);
    for (Size jj = 1; jj < maturities_.size(); ++jj) {
        Real d = std::abs(maturities_[jj] - T);
        if (d < best_dist) { best_dist = d; j = jj; }
    }
    Real T_actual = maturities_[j];

    SVIParams initial{0.04, 0.4, -0.2, 0.2, 0.0};  // reasonable defaults
    auto svi = std::make_unique<SVI>(initial, SVIParamType::Raw);
    CalibConfig cfg;
    cfg.de_pop_size = 30;
    cfg.de_generations = 100;
    cfg.lm_max_iter = 100;
    svi->calibrate(strikes_, {T_actual}, implied_vols_[j], forward_, cfg);
    svi_slices_[T_actual] = std::move(svi);
}

inline const SVI* VolSurface::svi_slice(Real T) const {
    auto it = svi_slices_.find(T);
    if (it != svi_slices_.end()) return it->second.get();
    // Try nearest
    it = svi_slices_.lower_bound(T);
    if (it != svi_slices_.end()) return it->second.get();
    if (it != svi_slices_.begin()) { --it; return it->second.get(); }
    return nullptr;
}

// Calendar arbitrage: total variance w(k, T) must be non-decreasing in T
inline bool VolSurface::check_calendar_arbitrage() const {
    if (maturities_.size() < 2) return true;
    for (Size j = 0; j < maturities_.size() - 1; ++j) {
        Real T1 = maturities_[j];
        Real T2 = maturities_[j + 1];
        for (Size i = 0; i < strikes_.size(); ++i) {
            Real K = strikes_[i];
            Real w1 = implied_vols_[j][i] * implied_vols_[j][i] * T1;
            Real w2 = implied_vols_[j + 1][i] * implied_vols_[j + 1][i] * T2;
            if (w2 < w1 - 1e-10) return false;
        }
    }
    return true;
}

// Butterfly arbitrage: BSM call price must be convex in K
// C(K-dK) - 2C(K) + C(K+dK) >= 0 (positive second derivative)
inline bool VolSurface::check_butterfly_arbitrage() const {
    for (Size j = 0; j < maturities_.size(); ++j) {
        Real T = maturities_[j];
        for (Size i = 1; i < strikes_.size() - 1; ++i) {
            Real K = strikes_[i];
            Real dK = std::min(K - strikes_[i - 1], strikes_[i + 1] - K);
            Real C_left = call_price(K - dK, T, forward_, r_, q_);
            Real C_mid = call_price(K, T, forward_, r_, q_);
            Real C_right = call_price(K + dK, T, forward_, r_, q_);
            Real second_deriv = (C_left - 2.0 * C_mid + C_right) / (dK * dK);
            if (second_deriv < -1e-8) return false;
        }
    }
    return true;
}

}  // namespace v1
}  // namespace cpphub
