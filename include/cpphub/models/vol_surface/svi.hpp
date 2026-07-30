#pragma once
// SOURCE: PHASE3_SPEC §4.2 - SVI volatility surface parameterization
// Implemented on main station (MSVC) - 2026-07-30
// SVI raw: w(k) = a + b * (rho*(k-m) + sqrt((k-m)^2 + sigma^2))
// No-arbitrage: butterfly (g(k) >= 0) + calendar (w_T1 <= w_T2)
// Ref: Gatheral & Jacquier (2014), "Arbitrage-free SVI volatility surfaces"
#include "cpphub/core/types.hpp"
#include "cpphub/calibration/optimizer.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace cpphub {
inline namespace v1 {

// SVI 原始参数化: w(k) = a + b * (ρ*(k-m) + sqrt((k-m)² + σ²))
// 其中 w = 总方差 (σ²_impl * T), k = log-moneyness (log(K/F))
struct SVIParams {
    Real a;       // 总方差水平 (a + b*|ρ|*σ >= 0 required for w >= 0)
    Real b;       // 斜率系数 (b >= 0)
    Real rho;     // 相关性 (-1 < ρ < 1)
    Real sigma;   // 弯曲度 (σ > 0)
    Real m;       // 中心偏移
};

enum class SVIParamType { Raw, Natural, JumpWings };

class SVI {
public:
    SVI(SVIParams params, SVIParamType type = SVIParamType::Raw);

    // Total variance w(k) = sigma_impl^2 * T
    Real total_variance(Real k) const;
    // Implied vol = sqrt(w(k) / T)
    Real implied_vol(Real k, Real T) const;
    // dw/dk
    Real first_derivative(Real k) const;
    // d²w/dk²
    Real second_derivative(Real k) const;

    // Butterfly arbitrage: g(k) >= 0 for all k (Gatheral-Jacquier)
    bool check_butterfly_arbitrage() const;
    bool check_butterfly_at(Real k) const;
    // Returns k values where g(k) < 0 (empty if no violation)
    std::vector<Real> find_arbitrage_violations(const std::vector<Real>& k_grid) const;

    // Parameter conversions
    static SVIParams raw_to_natural(const SVIParams& raw, Real T);
    static SVIParams natural_to_raw(const SVIParams& nat, Real T);
    static SVIParams raw_to_jump_wings(const SVIParams& raw, Real T);
    static SVIParams jump_wings_to_raw(const SVIParams& jw, Real T);

    // Calibrate to market quotes (strikes, maturities, implied vols) at fixed maturity T
    CalibrationResult calibrate(
        const std::vector<Real>& strikes,
        const std::vector<Real>& maturities,
        const std::vector<Real>& implied_vols,
        Real forward,
        const CalibConfig& cfg = CalibConfig{});

    SVIParams params() const { return params_; }
    SVIParamType type() const { return type_; }
    Real T() const { return T_; }

private:
    SVIParams params_;
    SVIParamType type_;
    Real T_ = 1.0;  // maturity (used for unit conversion)

    // Internal: work in raw params
    Real raw_total_variance(Real k) const;
    Real raw_dvar_dk(Real k) const;
    Real raw_d2var_dk2(Real k) const;
    // g(k) = (1 - k*w'/2w)^2 - (w'/2)^2 * (w + 0.25) + w''/2  ( Gatheral eq.)
    Real butterfly_g_function(Real k) const;
};

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

inline SVI::SVI(SVIParams params, SVIParamType type)
    : params_(std::move(params)), type_(type) {
    if (type_ != SVIParamType::Raw) {
        // Convert to raw for internal computation
        if (type_ == SVIParamType::Natural) {
            params_ = natural_to_raw(params_, T_);
        } else {  // JumpWings
            params_ = jump_wings_to_raw(params_, T_);
        }
        type_ = SVIParamType::Raw;
    }
}

inline Real SVI::raw_total_variance(Real k) const {
    Real km = k - params_.m;
    Real s = std::sqrt(km * km + params_.sigma * params_.sigma);
    return params_.a + params_.b * (params_.rho * km + s);
}

inline Real SVI::total_variance(Real k) const {
    return raw_total_variance(k);
}

inline Real SVI::implied_vol(Real k, Real T) const {
    if (T <= 0.0) {
        throw std::invalid_argument("SVI::implied_vol: T must be positive");
    }
    Real w = raw_total_variance(k);
    if (w < 0.0) return 0.0;  // violation, clamp
    return std::sqrt(w / T);
}

inline Real SVI::raw_dvar_dk(Real k) const {
    Real km = k - params_.m;
    Real s = std::sqrt(km * km + params_.sigma * params_.sigma);
    // dw/dk = b * (rho + km / s)
    return params_.b * (params_.rho + km / s);
}

inline Real SVI::first_derivative(Real k) const {
    return raw_dvar_dk(k);
}

inline Real SVI::raw_d2var_dk2(Real k) const {
    Real km = k - params_.m;
    Real s = std::sqrt(km * km + params_.sigma * params_.sigma);
    // d²w/dk² = b * (1/s - km²/s³) = b * σ² / s³
    Real s3 = s * s * s;
    return params_.b * params_.sigma * params_.sigma / s3;
}

inline Real SVI::second_derivative(Real k) const {
    return raw_d2var_dk2(k);
}

inline Real SVI::butterfly_g_function(Real k) const {
    Real w = raw_total_variance(k);
    Real wp = raw_dvar_dk(k);
    Real wpp = raw_d2var_dk2(k);
    if (w <= 0.0) return -1.0;  // w must be positive
    Real t = 1.0 - k * wp / (2.0 * w);
    Real g = t * t - wp * wp * 0.25 * (w + 0.25) + wpp * 0.5;
    return g;
}

inline bool SVI::check_butterfly_at(Real k) const {
    return butterfly_g_function(k) >= 0.0;
}

inline bool SVI::check_butterfly_arbitrage() const {
    // Check g(k) >= 0 on a dense grid. Theoretical min of g occurs near k=m.
    // Sample k from m-5*sigma to m+5*sigma with 1001 points
    Real k_min = params_.m - 5.0 * params_.sigma;
    Real k_max = params_.m + 5.0 * params_.sigma;
    Size n = 1001;
    for (Size i = 0; i < n; ++i) {
        Real k = k_min + (k_max - k_min) * static_cast<Real>(i) / static_cast<Real>(n - 1);
        if (butterfly_g_function(k) < 0.0) return false;
    }
    // Also check w >= 0 at extremes
    if (raw_total_variance(k_min) < 0.0) return false;
    if (raw_total_variance(k_max) < 0.0) return false;
    return true;
}

inline std::vector<Real> SVI::find_arbitrage_violations(const std::vector<Real>& k_grid) const {
    std::vector<Real> viol;
    for (Real k : k_grid) {
        if (butterfly_g_function(k) < 0.0) viol.push_back(k);
    }
    return viol;
}

// Natural parameters: (nu, psi, p, c, v_tilde) where
//   nu = a + b*sigma*sqrt(1-rho^2)  (minimum variance)
//   psi = b*rho/(2*sqrt(nu))        (ATM skew)
//   p = rho                          (correlation, same)
//   c = b/(2*nu)                     (curvature)
//   v_tilde = a                      (level offset)
inline SVIParams SVI::raw_to_natural(const SVIParams& raw, Real /*T*/) {
    Real nu = raw.a + raw.b * raw.sigma * std::sqrt(1.0 - raw.rho * raw.rho);
    Real psi = raw.b * raw.rho / (2.0 * std::sqrt(std::max(nu, 1e-12)));
    // Pack into SVIParams: {nu, psi, rho(p), c, m}
    // Note: natural params don't fit cleanly; we store as {nu, psi, p, c, m}
    SVIParams nat;
    nat.a = nu;
    nat.b = psi;
    nat.rho = raw.rho;  // p
    nat.sigma = raw.b / (2.0 * std::max(nu, 1e-12));  // c
    nat.m = raw.m;
    return nat;
}

inline SVIParams SVI::natural_to_raw(const SVIParams& nat, Real /*T*/) {
    // nat = {nu, psi, p, c, m}
    Real nu = nat.a;
    Real psi = nat.b;
    Real p = nat.rho;
    Real c = nat.sigma;
    SVIParams raw;
    raw.m = nat.m;
    raw.rho = p;
    // b = 2 * c * nu
    raw.b = 2.0 * c * nu;
    // sigma = sqrt(nu) * sqrt(1-p^2) / c   (from nu = a + b*sigma*sqrt(1-p^2) and b=2*c*nu)
    // But a = nu - b*sigma*sqrt(1-p^2) => need a. For pure natural->raw we need a.
    // This conversion is underdetermined without a; assume v_tilde = a = 0 (common choice)
    raw.a = 0.0;
    Real sqrt1p2 = std::sqrt(std::max(1.0 - p * p, 0.0));
    raw.sigma = (c > 0.0) ? (sqrt1p2 / c) : 1.0;
    // Recompute a from nu
    raw.a = nu - raw.b * raw.sigma * sqrt1p2;
    return raw;
}

// Jump-Wings (JW) parameters: (v_tilde, psi, p, c, v_tilde_atm)
//   v_tilde = (a + b*sigma*sqrt(1-rho^2)) / T   (ATM variance / T)
//   psi = rho/(2*sqrt(sigma*(1-rho^2))) * ... (ATM skew)
//   p = rho
//   c = b*T / (2*nu)   (curvature, time-scaled)
//   v_tilde_atm = (a + b*rho*m + b*sqrt(m^2+sigma^2)) / T  (ATM total variance)
// For simplicity we implement the standard Gatheral-Jacquier (2014) JW form.
inline SVIParams SVI::raw_to_jump_wings(const SVIParams& raw, Real T) {
    if (T <= 0.0) throw std::invalid_argument("raw_to_jump_wings: T must be positive");
    Real sqrt1p2 = std::sqrt(std::max(1.0 - raw.rho * raw.rho, 0.0));
    Real nu = raw.a + raw.b * raw.sigma * sqrt1p2;  // min variance
    Real atm_var = raw.a + raw.b * (raw.rho * (-raw.m) + std::sqrt(raw.m * raw.m + raw.sigma * raw.sigma));
    // JW: {v_tilde_min/T, psi, p, c, v_tilde_atm/T}
    SVIParams jw;
    jw.a = nu / T;                                              // v_tilde_min
    jw.b = (nu > 0.0) ? (raw.b * raw.rho / (2.0 * std::sqrt(nu))) : 0.0;  // psi (ATM skew)
    jw.rho = raw.rho;                                            // p
    jw.sigma = (nu > 0.0) ? (raw.b / (2.0 * nu)) : 0.0;         // c (curvature)
    jw.m = atm_var / T;                                          // v_tilde_atm
    return jw;
}

inline SVIParams SVI::jump_wings_to_raw(const SVIParams& jw, Real T) {
    if (T <= 0.0) throw std::invalid_argument("jump_wings_to_raw: T must be positive");
    // jw = {v_tilde_min, psi, p, c, v_tilde_atm}
    Real v_min = jw.a * T;
    Real psi = jw.b;
    Real p = jw.rho;
    Real c = jw.sigma;
    Real atm_var = jw.m * T;
    SVIParams raw;
    raw.rho = p;
    raw.b = 2.0 * c * v_min;
    Real sqrt1p2 = std::sqrt(std::max(1.0 - p * p, 0.0));
    raw.sigma = (c > 0.0) ? (sqrt1p2 / c) : 1.0;
    raw.a = v_min - raw.b * raw.sigma * sqrt1p2;
    // Solve for m: atm_var = a + b*(rho*(-m) + sqrt(m^2 + sigma^2))
    // This is transcendental; use m=0 as approximation (valid for symmetric smiles)
    raw.m = 0.0;
    return raw;
}

inline CalibrationResult SVI::calibrate(
        const std::vector<Real>& strikes,
        const std::vector<Real>& maturities,
        const std::vector<Real>& implied_vols,
        Real forward,
        const CalibConfig& cfg) {
    CalibrationResult result;
    if (strikes.empty() || maturities.empty()) {
        result.converged = false;
        result.message = "empty input data";
        return result;
    }

    // For a single-maturity SVI slice, pick the first maturity (or the one with most quotes)
    // Here we calibrate to the first maturity for simplicity.
    Real T = maturities[0];
    T_ = T;

    // Build target: total variance w = sigma_impl^2 * T at log-moneyness k = log(K/F)
    std::vector<Real> k_targets;
    std::vector<Real> w_targets;
    Size n = strikes.size();
    // Assume strikes are for the first maturity (1D calibration)
    for (Size i = 0; i < n; ++i) {
        Real k = std::log(strikes[i] / forward);
        Real w = implied_vols[i] * implied_vols[i] * T;
        k_targets.push_back(k);
        w_targets.push_back(w);
    }

    // Residual function: r_i = w_model(k_i) - w_target_i
    ResidualFn residual = [&](const std::vector<Real>& x) -> std::vector<Real> {
        SVIParams p{x[0], x[1], x[2], x[3], x[4]};
        // Clamp to valid region
        if (p.b < 0.0) p.b = 0.0;
        if (p.sigma <= 0.0) p.sigma = 1e-4;
        if (p.rho <= -1.0) p.rho = -0.999;
        if (p.rho >= 1.0) p.rho = 0.999;
        SVI s(p, SVIParamType::Raw);
        std::vector<Real> r(n);
        for (Size i = 0; i < n; ++i) {
            r[i] = s.total_variance(k_targets[i]) - w_targets[i];
        }
        return r;
    };

    // Initial guess: a=min(w), b=0.1, rho=0, sigma=0.1, m=0
    Real w_min = *std::min_element(w_targets.begin(), w_targets.end());
    std::vector<Real> x0 = {w_min * 0.9, 0.1, 0.0, 0.1, 0.0};

    // DE global search first (if enabled), then LM refine
    std::vector<Real> x_init = x0;
    if (cfg.use_de_init) {
        std::vector<Bounds> bounds = {
            {-1.0, 5.0},     // a
            {0.0, 5.0},      // b
            {-0.99, 0.99},   // rho
            {1e-4, 5.0},     // sigma
            {-2.0, 2.0}      // m
        };
        ObjectiveFn obj = [&](const std::vector<Real>& xx) -> Real {
            auto r = residual(xx);
            Real s = 0.0;
            for (Real v : r) s += v * v;
            return s;
        };
        DifferentialEvolution::Config de_cfg;
        de_cfg.population_size = cfg.de_pop_size;
        de_cfg.max_generations = cfg.de_generations;
        de_cfg.seed = cfg.seed;
        auto de_result = DifferentialEvolution::minimize(obj, bounds, de_cfg);
        x_init = de_result.x;
    }

    // LM refine
    LevenbergMarquardt::Config lm_cfg;
    lm_cfg.max_iterations = cfg.lm_max_iter;
    lm_cfg.ftol = cfg.ftol;
    lm_cfg.xtol = cfg.xtol;
    auto lm_result = LevenbergMarquardt::minimize(residual, x_init, lm_cfg);

    result.params = lm_result.x;
    result.objective_value = lm_result.fx;
    result.n_iterations = lm_result.n_iterations;
    result.converged = lm_result.converged;
    result.message = lm_result.message;
    result.residuals = residual(lm_result.x);

    // Update internal params
    params_ = SVIParams{lm_result.x[0], lm_result.x[1], lm_result.x[2], lm_result.x[3], lm_result.x[4]};
    type_ = SVIParamType::Raw;
    return result;
}

}  // inline namespace v1
}  // namespace cpphub
