#pragma once
// SOURCE: PHASE3_SPEC §4.1 - Model calibrators (Heston/SABR/SVI)
// Implemented on main station (MSVC) - 2026-07-31
// HestonCalibrator: CF integral pricing (Carr-Madan) + IV calibration (DE + LM)
// SABRCalibrator: Hagan 2002 explicit IV formula + IV calibration (DE + LM)
// SVI calibration is available via SVI::calibrate() in svi.hpp
// NOTE: CalibrationResult and CalibConfig are defined in optimizer.hpp
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/calibration/optimizer.hpp"
#include "cpphub/calibration/objective.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"  // bsm_call_price / bsm_put_price
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <complex>

namespace cpphub {
inline namespace v1 {

class Calibrator {
public:
    virtual ~Calibrator() = default;
    virtual CalibrationResult calibrate(
        const std::vector<MarketQuote>& quotes,
        const CalibConfig& cfg = CalibConfig{}) = 0;
    virtual std::string name() const = 0;
};

struct HestonParams {
    Real v0, kappa, theta, sigma_v, rho;
};

struct SABRParams {
    Real alpha, beta, nu, rho;
};

// ---------------------------------------------------------------------------
// BSM implied volatility inversion (Newton-Raphson with vega)
// Needed by both Heston and SABR calibrators when market quotes are prices
// but calibration target is IV (common practice for equity/FX vol surfaces)
// ---------------------------------------------------------------------------
inline Real bsm_implied_vol(Real C_market, Real S, Real K, Real T,
                             Real r, Real q, bool is_call,
                             Real tol = 1e-10, int max_iter = 50) {
    if (T <= 0.0 || C_market <= 0.0) return 0.0;
    // Initial guess: Brenner-Subrahmanyam (1988) ATM approximation
    Real sigma = std::sqrt(2.0 * M_PI / T) * C_market / S;
    if (sigma <= 0.001) sigma = 0.20;
    sigma = std::max(0.001, std::min(5.0, sigma));

    for (int iter = 0; iter < max_iter; ++iter) {
        Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
        Real d2 = d1 - sigma * std::sqrt(T);
        Real price = is_call
            ? (S * std::exp(-q * T) * normal_cdf(d1) - K * std::exp(-r * T) * normal_cdf(d2))
            : (K * std::exp(-r * T) * normal_cdf(-d2) - S * std::exp(-q * T) * normal_cdf(-d1));
        Real vega = S * std::exp(-q * T) * normal_pdf(d1) * std::sqrt(T);
        if (vega < 1e-12) break;
        Real diff = price - C_market;
        if (std::abs(diff) < tol) return sigma;
        Real step = diff / vega;
        sigma -= step;
        if (sigma <= 0.001) sigma = 0.0015;
        if (sigma >= 5.0) sigma = 4.99;
    }
    return sigma;
}

// ---------------------------------------------------------------------------
// Heston call price via Carr-Madan Fourier inversion (direct integral)
// C(k) = S e^{-qT} - K e^{-rT} (1 - F1) / 2 + ... (simplified, using P1+P2)
// We use the standard Gil-Pelaez inversion:
//   P_j = 0.5 + (1/π) ∫_0^∞ Im[ e^{-iu ln K} φ_j(u) / (iu) ] du
//   C = S e^{-qT} P1 - K e^{-rT} P2
// where φ_1(u) = φ(u - i) * e^{...} and φ_2(u) = φ(u)
// ---------------------------------------------------------------------------
namespace detail {

inline Real heston_call_price_cf(Real S, Real K, Real T, Real r, Real q,
                                  const HestonParams& hp, Size n_quad = 4096) {
    if (T <= 0.0) return std::max(S - K, 0.0);

    HestonCFParams p;
    p.v0 = hp.v0;
    p.kappa = hp.kappa;
    p.theta = hp.theta;
    p.sigma = hp.sigma_v;
    p.rho = hp.rho;
    p.r = r;
    p.q = q;

    const Complex I(0, 1);
    Real log_K = std::log(K);

    // Trapezoidal quadrature on [0, umax]
    Real umax = 200.0;  // sufficient for typical params
    Real du = umax / static_cast<Real>(n_quad);

    // P2 = Pr(S_T >= K) under risk-neutral measure
    // Gil-Pelaez: P2 = 0.5 + (1/π) ∫_0^∞ Re[ e^{-iu log K} φ(u) / (iu) ] du
    // Note: Re[f/(iu)] = Im[f]/u  (the standard form). Using .imag() here was a bug
    // since Im[f/(iu)] = -Re[f]/u, which produces large negative values near u=0.
    Real integral2 = 0.0;
    for (Size i = 0; i < n_quad; ++i) {
        Real u = (static_cast<Real>(i) + 0.5) * du;
        Complex iu = I * u;
        Complex phi = heston_characteristic_function(u, T, S, p);
        Complex integrand = std::exp(-iu * log_K) * phi / iu;
        integral2 += integrand.real() * du;
    }
    Real P2 = 0.5 + integral2 / M_PI;

    // P1 = Pr(S_T >= K) under stock measure (use φ(u - i) / S e^{(r-q)T})
    // φ_1(u) = φ(u - i) * e^{-i u (r-q)T + (r-q)T} / ... actually:
    // φ_1(u) = e^{i u (r-q)T} * φ(u - i) * S / (e^{(r-q)T} S) = φ(u-i) * e^{...} ...
    // Standard form: φ_1(u) = φ(u - i, S, p) where φ uses shifted characteristic fn.
    // Simpler: use the identity P1 = S e^{(r-q)T} / K * (something), but easier to
    // directly compute φ_1 by calling heston_characteristic_function with u -> u - I
    // and multiplying by the discount factor.
    // The standard Carr-Madan formulation: φ_1(u) = φ(u - i) / φ(-i)
    // where φ(-i) = E[S_T/S] = e^{(r-q)T} (martingale condition).
    Real integral1 = 0.0;
    Complex phi_neg_i = heston_characteristic_function(Complex(0, -1), T, S, p);
    for (Size i = 0; i < n_quad; ++i) {
        Real u = (static_cast<Real>(i) + 0.5) * du;
        Complex u_shifted = u - I;  // u - i
        Complex phi1 = heston_characteristic_function(u_shifted, T, S, p) / phi_neg_i;
        Complex iu = I * u;
        Complex integrand = std::exp(-iu * log_K) * phi1 / iu;
        integral1 += integrand.real() * du;
    }
    Real P1 = 0.5 + integral1 / M_PI;

    Real call = S * std::exp(-q * T) * P1 - K * std::exp(-r * T) * P2;
    return std::max(call, 0.0);
}

// SABR Hagan 2002 implied vol formula (normal vol for lognormal model)
// σ_B(K, F) = α / (F^(1-β) * (z χ(z) + ... ))  (truncated expansion)
// Standard Hagan formula for lognormal implied vol:
inline Real sabr_implied_vol_hagan(Real F, Real K, Real T,
                                    const SABRParams& sp) {
    if (T <= 0.0 || F <= 0.0 || K <= 0.0) return 0.0;
    Real alpha = sp.alpha;
    Real beta = sp.beta;
    Real nu = sp.nu;
    Real rho = sp.rho;

    Real logFK = std::log(F / K);
    Real FK_beta = std::pow(F * K, (1.0 - beta) * 0.5);
    Real one_minus_beta = 1.0 - beta;

    // z = (nu/alpha) * (F*K)^((1-beta)/2) * log(F/K)
    Real z = (nu / alpha) * FK_beta * logFK;

    // χ(z) = ln[(sqrt(1-2ρz+z²) + z - ρ) / (1 - ρ)] / z
    Real chi_z;
    if (std::abs(z) < 1e-8) {
        // L'Hôpital: χ(z) → 1 as z → 0
        chi_z = 1.0;
    } else {
        Real sqrt_term = std::sqrt(1.0 - 2.0 * rho * z + z * z);
        Real num = std::log((sqrt_term + z - rho) / (1.0 - rho));
        chi_z = num / z;
    }

    // ATM case (K == F): z = 0, χ = 1
    // σ_atm = α / F^(1-β) * [1 + ((1-β)² α²/(24 F^(2-2β)) + ρβα/(4F^(1-β)) + (2-3ρ²)ν²/24) T]
    if (std::abs(logFK) < 1e-10) {
        Real F_beta = std::pow(F, one_minus_beta);
        Real term1 = one_minus_beta * one_minus_beta * alpha * alpha /
                     (24.0 * F_beta * F_beta);
        Real term2 = rho * beta * nu * alpha / (4.0 * F_beta);
        Real term3 = (2.0 - 3.0 * rho * rho) * nu * nu / 24.0;
        return (alpha / F_beta) * (1.0 + (term1 + term2 + term3) * T);
    }

    // General case
    Real F_beta = std::pow(F, one_minus_beta);
    Real term1 = one_minus_beta * one_minus_beta * alpha * alpha /
                 (24.0 * FK_beta * FK_beta);
    Real term2 = rho * beta * nu * alpha / (4.0 * FK_beta);
    Real term3 = (2.0 - 3.0 * rho * rho) * nu * nu / 24.0;

    Real numerator = alpha * (1.0 + (term1 + term2 + term3) * T);
    Real denominator = FK_beta * (1.0 + one_minus_beta * one_minus_beta *
                                  logFK * logFK / 24.0 +
                                  one_minus_beta * one_minus_beta *
                                  one_minus_beta * one_minus_beta *
                                  logFK * logFK * logFK * logFK / 1920.0);
    return numerator / denominator * (z / chi_z);
}

}  // namespace detail

// ---------------------------------------------------------------------------
// HestonCalibrator
// Parameter vector x = [v0, kappa, theta, sigma_v, rho]
// Objective: sum_i w_i (IV_model_i - IV_market_i)^2
// Constraints: v0>0, kappa>0, theta>0, sigma_v>0, |rho|<1, Feller 2κθ>σ²
// ---------------------------------------------------------------------------
class HestonCalibrator : public Calibrator {
public:
    CalibrationResult calibrate(
        const std::vector<MarketQuote>& quotes,
        const CalibConfig& cfg = CalibConfig{}) override;
    std::string name() const override { return "HestonCalibrator"; }

    static std::vector<Bounds> default_bounds() {
        return {{1e-4, 1.0},      // v0
                {1e-4, 10.0},     // kappa
                {1e-4, 1.0},       // theta
                {1e-4, 5.0},       // sigma_v
                {-0.99, 0.99}};    // rho
    }
    static bool check_feller(const HestonParams& p) {
        return 2.0 * p.kappa * p.theta > p.sigma_v * p.sigma_v;
    }
    HestonParams extract_params(const std::vector<Real>& x) const {
        return HestonParams{x[0], x[1], x[2], x[3], x[4]};
    }

    // Set market context (spot, rate, dividend) — required before calibrate
    void set_market(Real S, Real r, Real q) { S_ = S; r_ = r; q_ = q; }

private:
    Real S_ = 100.0;
    Real r_ = 0.0;
    Real q_ = 0.0;
};

inline CalibrationResult HestonCalibrator::calibrate(
        const std::vector<MarketQuote>& quotes, const CalibConfig& cfg) {
    CalibrationResult result;
    if (quotes.empty()) {
        result.converged = false;
        result.message = "empty quotes";
        return result;
    }

    // IV model function: given Heston params (v0, kappa, theta, sigma, rho),
    // return model IV for each quote (K, T).
    auto iv_fn = [this](const std::vector<Real>& x, Real K, Real T) -> Real {
        HestonParams hp{x[0], x[1], x[2], x[3], x[4]};
        // Clamp to valid region
        if (hp.v0 <= 0.0) hp.v0 = 1e-4;
        if (hp.kappa <= 0.0) hp.kappa = 1e-4;
        if (hp.theta <= 0.0) hp.theta = 1e-4;
        if (hp.sigma_v <= 0.0) hp.sigma_v = 1e-4;
        if (hp.rho <= -1.0) hp.rho = -0.999;
        if (hp.rho >= 1.0) hp.rho = 0.999;
        Real price = detail::heston_call_price_cf(S_, K, T, r_, q_, hp);
        // Invert BSM to get IV
        return bsm_implied_vol(price, S_, K, T, r_, q_, true);
    };

    auto obj = ObjectiveFunction::make_iv_objective(iv_fn, quotes,
                 WeightingScheme::RelativeError);

    // DE global search
    std::vector<Real> x_init(5, 0.0);
    if (cfg.use_de_init) {
        auto bounds = default_bounds();
        DifferentialEvolution::Config de_cfg;
        de_cfg.population_size = cfg.de_pop_size;
        de_cfg.max_generations = cfg.de_generations;
        de_cfg.seed = cfg.seed;
        auto de_result = DifferentialEvolution::minimize(obj.to_objective_fn(), bounds, de_cfg);
        x_init = de_result.x;
    } else {
        // Reasonable defaults
        x_init = {0.04, 1.0, 0.04, 0.3, -0.5};
    }

    // LM refine
    LevenbergMarquardt::Config lm_cfg;
    lm_cfg.max_iterations = cfg.lm_max_iter;
    lm_cfg.ftol = cfg.ftol;
    lm_cfg.xtol = cfg.xtol;
    auto lm_result = LevenbergMarquardt::minimize(obj.to_residual_fn(), x_init, lm_cfg);

    result.params = lm_result.x;
    result.objective_value = lm_result.fx;
    result.n_iterations = lm_result.n_iterations;
    result.converged = lm_result.converged;
    result.message = lm_result.message;
    result.residuals = obj.residuals(lm_result.x);

    // Feller check
    HestonParams hp = extract_params(result.params);
    if (!check_feller(hp)) {
        result.message += " [WARNING: Feller condition violated]";
    }
    return result;
}

// ---------------------------------------------------------------------------
// SABRCalibrator
// Parameter vector x = [alpha, beta, nu, rho]
// Objective: sum_i w_i (IV_model_i - IV_market_i)^2
// beta is usually fixed (0.5 equity, 1.0 FX, 0.0 rates); here we calibrate all 4
// Constraints: alpha>0, 0<=beta<=1, nu>0, |rho|<1
// ---------------------------------------------------------------------------
class SABRCalibrator : public Calibrator {
public:
    CalibrationResult calibrate(
        const std::vector<MarketQuote>& quotes,
        const CalibConfig& cfg = CalibConfig{}) override;
    std::string name() const override { return "SABRCalibrator"; }

    static std::vector<Bounds> default_bounds() {
        return {{1e-4, 5.0},      // alpha
                {0.0, 1.0},        // beta
                {1e-4, 5.0},       // nu
                {-0.99, 0.99}};    // rho
    }
    SABRParams extract_params(const std::vector<Real>& x) const {
        return SABRParams{x[0], x[1], x[2], x[3]};
    }

    void set_market(Real F, Real r, Real q) {
        F_ = F; r_ = r; q_ = q;
        (void)q_;  // unused for SABR (uses forward F = S e^{(r-q)T})
    }

private:
    Real F_ = 100.0;
    Real r_ = 0.0;
    Real q_ = 0.0;
};

inline CalibrationResult SABRCalibrator::calibrate(
        const std::vector<MarketQuote>& quotes, const CalibConfig& cfg) {
    CalibrationResult result;
    if (quotes.empty()) {
        result.converged = false;
        result.message = "empty quotes";
        return result;
    }

    // IV model function via Hagan 2002
    auto iv_fn = [this](const std::vector<Real>& x, Real K, Real T) -> Real {
        SABRParams sp{x[0], x[1], x[2], x[3]};
        if (sp.alpha <= 0.0) sp.alpha = 1e-4;
        if (sp.beta < 0.0) sp.beta = 0.0;
        if (sp.beta > 1.0) sp.beta = 1.0;
        if (sp.nu <= 0.0) sp.nu = 1e-4;
        if (sp.rho <= -1.0) sp.rho = -0.999;
        if (sp.rho >= 1.0) sp.rho = 0.999;
        // Forward at maturity: F = S * e^{(r-q)T}
        Real F = F_ * std::exp((r_ - q_) * T);
        return detail::sabr_implied_vol_hagan(F, K, T, sp);
    };

    auto obj = ObjectiveFunction::make_iv_objective(iv_fn, quotes,
                 WeightingScheme::RelativeError);

    std::vector<Real> x_init(4, 0.0);
    if (cfg.use_de_init) {
        auto bounds = default_bounds();
        DifferentialEvolution::Config de_cfg;
        de_cfg.population_size = cfg.de_pop_size;
        de_cfg.max_generations = cfg.de_generations;
        de_cfg.seed = cfg.seed;
        auto de_result = DifferentialEvolution::minimize(obj.to_objective_fn(), bounds, de_cfg);
        x_init = de_result.x;
    } else {
        x_init = {0.2, 0.5, 0.3, -0.2};  // Reasonable defaults
    }

    LevenbergMarquardt::Config lm_cfg;
    lm_cfg.max_iterations = cfg.lm_max_iter;
    lm_cfg.ftol = cfg.ftol;
    lm_cfg.xtol = cfg.xtol;
    auto lm_result = LevenbergMarquardt::minimize(obj.to_residual_fn(), x_init, lm_cfg);

    result.params = lm_result.x;
    result.objective_value = lm_result.fx;
    result.n_iterations = lm_result.n_iterations;
    result.converged = lm_result.converged;
    result.message = lm_result.message;
    result.residuals = obj.residuals(lm_result.x);
    return result;
}

}  // inline namespace v1
}  // namespace cpphub
