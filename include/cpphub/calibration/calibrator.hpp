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
#include "cpphub/pricing/analytic/bates_cf.hpp"             // make_bates_cf (Bates CF)
#include "cpphub/pricing/analytic/cev_analytic.hpp"         // cev_call_price (CEV analytic)
#include "cpphub/pricing/fourier/characteristic_functions.hpp"  // make_vg_cf
#include "cpphub/pricing/analytic/vg_analytic.hpp"               // VGParams
#include "cpphub/pricing/fourier/cos_method.hpp"            // COSEngine (Fang-Oosterlee 2009)
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
//
// 数值加固: 纯 Newton 在短到期/深度 OTM 时会发散 (sigma 冲出 [0, 5] 或振荡),
// 导致 IV 崩到下界。改为 Newton + 二分混合:
//   - 先建立 [lo, hi] 价格包围区间
//   - Newton 步留在区间内则接受, 越界则退化为二分
//   - 每次更新后收缩包围区间, 保证单调收敛
// ---------------------------------------------------------------------------
inline Real bsm_implied_vol(Real C_market, Real S, Real K, Real T,
                             Real r, Real q, bool is_call,
                             Real tol = 1e-10, int max_iter = 50) {
    if (T <= 0.0 || C_market <= 0.0) return 0.0;

    const Real S_disc = S * std::exp(-q * T);
    const Real K_disc = K * std::exp(-r * T);
    const Real sqT = std::sqrt(T);

    // 期权价格函数 (sigma > 0)
    auto bsm_price = [&](Real sigma) -> Real {
        Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqT);
        Real d2 = d1 - sigma * sqT;
        if (is_call) {
            return S_disc * normal_cdf(d1) - K_disc * normal_cdf(d2);
        } else {
            return K_disc * normal_cdf(-d2) - S_disc * normal_cdf(-d1);
        }
    };

    // 无套利价格边界: 内在价值下限 与 资产上界
    Real intrinsic = is_call ? (S_disc - K_disc) : (K_disc - S_disc);
    Real min_price = std::max(intrinsic, 0.0);
    Real max_price = is_call ? S_disc : K_disc;
    if (C_market <= min_price * (1.0 + 1e-12) + 1e-14) return 0.001;
    if (C_market >= max_price * (1.0 - 1e-12) - 1e-14) return 5.0;

    // 初始包围区间 [lo, hi]: 保证 price(lo) <= C <= price(hi)
    Real lo = 0.0001, hi = 5.0;
    while (bsm_price(hi) < C_market && hi < 100.0) hi *= 2.0;
    if (bsm_price(hi) < C_market) return 5.0;
    while (bsm_price(lo) > C_market) lo *= 0.5;

    // Initial guess: Brenner-Subrahmanyam (1988) ATM approximation
    Real sigma = std::sqrt(2.0 * M_PI / T) * C_market / S;
    sigma = std::max(lo, std::min(hi, sigma));

    for (int iter = 0; iter < max_iter; ++iter) {
        Real p = bsm_price(sigma);
        if (std::abs(p - C_market) < tol) return sigma;

        Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqT);
        Real vega = S_disc * normal_pdf(d1) * sqT;

        // 候选 Newton 步
        Real sigma_new = (vega > 1e-14)
            ? sigma - (p - C_market) / vega
            : 0.5 * (lo + hi);
        // 越界或步进异常 → 二分
        if (sigma_new <= lo || sigma_new >= hi ||
            !(sigma_new > 0.0) || std::isnan(sigma_new)) {
            sigma_new = 0.5 * (lo + hi);
        }

        // 根据符号收缩包围区间
        if (p > C_market) hi = sigma; else lo = sigma;
        sigma = sigma_new;

        if (hi - lo < 1e-10 * std::max(1.0, hi)) {
            return 0.5 * (lo + hi);
        }
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
        de_cfg.lambda_reg = cfg.lambda_reg;
        de_cfg.params_prior = cfg.params_prior;
        de_cfg.early_stop_rmse = cfg.early_stop_rmse;
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
    lm_cfg.lambda_reg = cfg.lambda_reg;
    lm_cfg.params_prior = cfg.params_prior;
    lm_cfg.early_stop_rmse = cfg.early_stop_rmse;
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

    // Full 4-parameter bounds (alpha, beta, nu, rho)
    static std::vector<Bounds> default_bounds() {
        return {{1e-4, 5.0},      // alpha
                {0.0, 1.0},        // beta
                {1e-4, 5.0},       // nu
                {-0.99, 0.99}};    // rho
    }
    // Reduced 3-parameter bounds when beta is fixed (alpha, nu, rho)
    static std::vector<Bounds> default_bounds_fixed_beta() {
        return {{1e-4, 5.0},      // alpha
                {1e-4, 5.0},       // nu
                {-0.99, 0.99}};    // rho
    }
    SABRParams extract_params(const std::vector<Real>& x) const {
        if (has_fixed_beta_) {
            return SABRParams{x[0], fixed_beta_, x[1], x[2]};
        }
        return SABRParams{x[0], x[1], x[2], x[3]};
    }

    void set_market(Real F, Real r, Real q) {
        F_ = F; r_ = r; q_ = q;
        (void)q_;  // unused for SABR (uses forward F = S e^{(r-q)T})
    }

    // Fix beta to a constant value (common practice: equity=0.5, FX=0.0 or 1.0, rates=0.5).
    // After calling this, calibrate() operates on a 3-parameter (alpha, nu, rho) problem.
    void set_fixed_beta(Real beta) {
        if (beta < 0.0 || beta > 1.0) {
            throw std::invalid_argument("SABRCalibrator::set_fixed_beta: beta must be in [0,1]");
        }
        fixed_beta_ = beta;
        has_fixed_beta_ = true;
    }
    void clear_fixed_beta() { has_fixed_beta_ = false; }
    bool has_fixed_beta() const { return has_fixed_beta_; }
    Real fixed_beta() const { return fixed_beta_; }

private:
    Real F_ = 100.0;
    Real r_ = 0.0;
    Real q_ = 0.0;
    Real fixed_beta_ = 0.5;       // valid only when has_fixed_beta_ = true
    bool has_fixed_beta_ = false;
};

inline CalibrationResult SABRCalibrator::calibrate(
        const std::vector<MarketQuote>& quotes, const CalibConfig& cfg) {
    CalibrationResult result;
    if (quotes.empty()) {
        result.converged = false;
        result.message = "empty quotes";
        return result;
    }

    // Branch on fixed-beta mode:
    //   - has_fixed_beta_ = true  → 3-param (alpha, nu, rho), beta = fixed_beta_
    //   - has_fixed_beta_ = false → 4-param (alpha, beta, nu, rho)
    if (has_fixed_beta_) {
        // IV model function via Hagan 2002 (3-param form, beta fixed)
        auto iv_fn = [this](const std::vector<Real>& x, Real K, Real T) -> Real {
            SABRParams sp{x[0], fixed_beta_, x[1], x[2]};
            if (sp.alpha <= 0.0) sp.alpha = 1e-4;
            if (sp.nu <= 0.0) sp.nu = 1e-4;
            if (sp.rho <= -1.0) sp.rho = -0.999;
            if (sp.rho >= 1.0) sp.rho = 0.999;
            Real F = F_ * std::exp((r_ - q_) * T);
            return detail::sabr_implied_vol_hagan(F, K, T, sp);
        };

        auto obj = ObjectiveFunction::make_iv_objective(iv_fn, quotes,
                     WeightingScheme::RelativeError);

        std::vector<Real> x_init(3, 0.0);
        if (cfg.use_de_init) {
            auto bounds = default_bounds_fixed_beta();
            DifferentialEvolution::Config de_cfg;
            de_cfg.population_size = cfg.de_pop_size;
            de_cfg.max_generations = cfg.de_generations;
            de_cfg.seed = cfg.seed;
            de_cfg.lambda_reg = cfg.lambda_reg;
            de_cfg.params_prior = cfg.params_prior;
            de_cfg.early_stop_rmse = cfg.early_stop_rmse;
            auto de_result = DifferentialEvolution::minimize(obj.to_objective_fn(), bounds, de_cfg);
            x_init = de_result.x;
        } else {
            x_init = {0.2, 0.3, -0.2};  // (alpha, nu, rho) defaults
        }

        LevenbergMarquardt::Config lm_cfg;
        lm_cfg.max_iterations = cfg.lm_max_iter;
        lm_cfg.ftol = cfg.ftol;
        lm_cfg.xtol = cfg.xtol;
        lm_cfg.lambda_reg = cfg.lambda_reg;
        lm_cfg.params_prior = cfg.params_prior;
        lm_cfg.early_stop_rmse = cfg.early_stop_rmse;
        auto lm_result = LevenbergMarquardt::minimize(obj.to_residual_fn(), x_init, lm_cfg);

        result.params = lm_result.x;
        result.objective_value = lm_result.fx;
        result.n_iterations = lm_result.n_iterations;
        result.converged = lm_result.converged;
        result.message = lm_result.message;
        result.residuals = obj.residuals(lm_result.x);
        return result;
    }

    // Full 4-parameter calibration (original path)
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
        de_cfg.lambda_reg = cfg.lambda_reg;
        de_cfg.params_prior = cfg.params_prior;
        de_cfg.early_stop_rmse = cfg.early_stop_rmse;
        auto de_result = DifferentialEvolution::minimize(obj.to_objective_fn(), bounds, de_cfg);
        x_init = de_result.x;
    } else {
        x_init = {0.2, 0.5, 0.3, -0.2};  // Reasonable defaults
    }

    LevenbergMarquardt::Config lm_cfg;
    lm_cfg.max_iterations = cfg.lm_max_iter;
    lm_cfg.ftol = cfg.ftol;
    lm_cfg.xtol = cfg.xtol;
    lm_cfg.lambda_reg = cfg.lambda_reg;
    lm_cfg.params_prior = cfg.params_prior;
    lm_cfg.early_stop_rmse = cfg.early_stop_rmse;
    auto lm_result = LevenbergMarquardt::minimize(obj.to_residual_fn(), x_init, lm_cfg);

    result.params = lm_result.x;
    result.objective_value = lm_result.fx;
    result.n_iterations = lm_result.n_iterations;
    result.converged = lm_result.converged;
    result.message = lm_result.message;
    result.residuals = obj.residuals(lm_result.x);
    return result;
}

// ---------------------------------------------------------------------------
// BatesCalibrator
// Parameter vector x = [v0, kappa, theta, sigma_v, rho, lambda, mu_J, sigma_J]
// 8-parameter calibration (Heston 5 + Merton jump 3)
// Objective: sum_i w_i (IV_model_i - IV_market_i)^2
// Pricing: Bates CF (make_bates_cf) + COS method (Fang-Oosterlee 2009),
//          IV recovered via bsm_implied_vol
// Constraints: v0>0, kappa>0, theta>0, sigma_v>0, |rho|<1,
//              lambda>0, mu_J∈R, sigma_J>0
// ---------------------------------------------------------------------------
struct BatesParams {
    Real v0, kappa, theta, sigma_v, rho;  // Heston 部分
    Real lambda, mu_J, sigma_J;            // Merton 跳跃部分
};

class BatesCalibrator : public Calibrator {
public:
    CalibrationResult calibrate(
        const std::vector<MarketQuote>& quotes,
        const CalibConfig& cfg = CalibConfig{}) override;
    std::string name() const override { return "BatesCalibrator"; }

    static std::vector<Bounds> default_bounds() {
        return {{1e-4, 1.0},      // v0
                {1e-4, 10.0},     // kappa
                {1e-4, 1.0},      // theta
                {1e-4, 5.0},      // sigma_v
                {-0.99, 0.99},    // rho
                {1e-4, 5.0},      // lambda (跳跃强度)
                {-0.5, 0.5},      // mu_J (跳跃均值)
                {1e-4, 1.0}};     // sigma_J (跳跃波动率)
    }
    static bool check_feller(const BatesParams& p) {
        return 2.0 * p.kappa * p.theta > p.sigma_v * p.sigma_v;
    }
    BatesParams extract_params(const std::vector<Real>& x) const {
        return BatesParams{x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]};
    }
    void set_market(Real S, Real r, Real q) { S_ = S; r_ = r; q_ = q; }

private:
    Real S_ = 100.0;
    Real r_ = 0.0;
    Real q_ = 0.0;
};

inline CalibrationResult BatesCalibrator::calibrate(
        const std::vector<MarketQuote>& quotes, const CalibConfig& cfg) {
    CalibrationResult result;
    if (quotes.empty()) {
        result.converged = false;
        result.message = "empty quotes";
        return result;
    }

    auto bounds = default_bounds();
    auto iv_fn = [this, &bounds](const std::vector<Real>& x, Real K, Real T) -> Real {
        BatesCFParams p;
        p.v0      = std::clamp(x[0], bounds[0].lower, bounds[0].upper);
        p.kappa   = std::clamp(x[1], bounds[1].lower, bounds[1].upper);
        p.theta   = std::clamp(x[2], bounds[2].lower, bounds[2].upper);
        p.sigma   = std::clamp(x[3], bounds[3].lower, bounds[3].upper);
        p.rho     = std::clamp(x[4], bounds[4].lower, bounds[4].upper);
        p.lambda  = std::clamp(x[5], bounds[5].lower, bounds[5].upper);
        p.mu_J    = std::clamp(x[6], bounds[6].lower, bounds[6].upper);
        p.sigma_J = std::clamp(x[7], bounds[7].lower, bounds[7].upper);
        p.r = r_;
        p.q = q_;
        auto phi = make_bates_cf(S_, r_, q_, p, T);
        COSEngine::Config cos_cfg;
        cos_cfg.n_terms = 256;
        COSEngine engine(phi, S_, r_, q_, T, cos_cfg);
        Real price = engine.price_call(K);
        return bsm_implied_vol(price, S_, K, T, r_, q_, true);
    };

    auto obj = ObjectiveFunction::make_iv_objective(iv_fn, quotes,
                 WeightingScheme::RelativeError);

    // DE global search
    std::vector<Real> x_init(8, 0.0);
    if (cfg.use_de_init) {
        DifferentialEvolution::Config de_cfg;
        de_cfg.population_size = cfg.de_pop_size;
        de_cfg.max_generations = cfg.de_generations;
        de_cfg.seed = cfg.seed;
        de_cfg.lambda_reg = cfg.lambda_reg;
        de_cfg.params_prior = cfg.params_prior;
        de_cfg.early_stop_rmse = cfg.early_stop_rmse;
        auto de_result = DifferentialEvolution::minimize(obj.to_objective_fn(), bounds, de_cfg);
        x_init = de_result.x;
    } else {
        x_init = {0.04, 1.0, 0.04, 0.3, -0.5, 0.2, 0.0, 0.1};
    }

    // LM refine
    LevenbergMarquardt::Config lm_cfg;
    lm_cfg.max_iterations = cfg.lm_max_iter;
    lm_cfg.ftol = cfg.ftol;
    lm_cfg.xtol = cfg.xtol;
    lm_cfg.lambda_reg = cfg.lambda_reg;
    lm_cfg.params_prior = cfg.params_prior;
    lm_cfg.early_stop_rmse = cfg.early_stop_rmse;
    auto lm_result = LevenbergMarquardt::minimize(obj.to_residual_fn(), x_init, lm_cfg);

    result.params = lm_result.x;
    result.objective_value = lm_result.fx;
    result.n_iterations = lm_result.n_iterations;
    result.converged = lm_result.converged;
    result.message = lm_result.message;
    result.residuals = obj.residuals(lm_result.x);

    BatesParams bp = extract_params(result.params);
    if (!check_feller(bp)) {
        result.message += " [WARNING: Feller condition violated]";
    }
    return result;
}

// ---------------------------------------------------------------------------
// VGCalibrator
// Parameter vector x = [sigma, nu, theta]
// 3-parameter calibration (VG pure-Lévy process, no jump intensity)
// Objective: sum_i w_i (IV_model_i - IV_market_i)^2
// Pricing: VG CF (make_vg_cf) + COS method, IV via bsm_implied_vol
// Feller 条件: 1 - theta*nu - sigma^2*nu/2 > 0 (标定中 clamp nu)
// 参数结构 VGParams 复用 cpphub/pricing/analytic/vg_analytic.hpp 中的定义
// ---------------------------------------------------------------------------

// Clamp nu to satisfy the VG Feller condition 1 - theta*nu - sigma^2*nu/2 > 0
inline Real vg_clamp_feller(Real sigma, Real nu, Real theta) {
    // 1 - theta*nu - sigma^2*nu/2 > 0  ⇔  nu*(theta + sigma^2/2) < 1
    Real denom = theta + 0.5 * sigma * sigma;
    if (denom > 0.0) {
        nu = std::min(nu, 0.95 / denom);  // 0.95 留出安全余量, 保证严格 > 0
    }
    return std::max(nu, 1e-6);
}

class VGCalibrator : public Calibrator {
public:
    CalibrationResult calibrate(
        const std::vector<MarketQuote>& quotes,
        const CalibConfig& cfg = CalibConfig{}) override;
    std::string name() const override { return "VGCalibrator"; }

    static std::vector<Bounds> default_bounds() {
        return {{1e-4, 2.0},      // sigma
                {1e-4, 5.0},      // nu
                {-1.0, 1.0}};     // theta
    }
    static bool check_feller(const VGParams& p) {
        return 1.0 - p.theta * p.nu - 0.5 * p.sigma * p.sigma * p.nu > 0.0;
    }
    VGParams extract_params(const std::vector<Real>& x) const {
        return VGParams{x[0], x[1], x[2]};
    }
    void set_market(Real S, Real r, Real q) { S_ = S; r_ = r; q_ = q; }

private:
    Real S_ = 100.0;
    Real r_ = 0.0;
    Real q_ = 0.0;
};

inline CalibrationResult VGCalibrator::calibrate(
        const std::vector<MarketQuote>& quotes, const CalibConfig& cfg) {
    CalibrationResult result;
    if (quotes.empty()) {
        result.converged = false;
        result.message = "empty quotes";
        return result;
    }

    auto bounds = default_bounds();
    auto iv_fn = [this, &bounds](const std::vector<Real>& x, Real K, Real T) -> Real {
        Real sigma = std::clamp(x[0], bounds[0].lower, bounds[0].upper);
        Real nu    = std::clamp(x[1], bounds[1].lower, bounds[1].upper);
        Real theta = std::clamp(x[2], bounds[2].lower, bounds[2].upper);
        // VG Feller 条件: 1 - theta*nu - sigma^2*nu/2 > 0, 违反时 clamp nu
        nu = vg_clamp_feller(sigma, nu, theta);
        auto phi = make_vg_cf(S_, r_, q_, sigma, nu, theta, T);
        COSEngine::Config cos_cfg;
        cos_cfg.n_terms = 256;
        COSEngine engine(phi, S_, r_, q_, T, cos_cfg);
        Real price = engine.price_call(K);
        return bsm_implied_vol(price, S_, K, T, r_, q_, true);
    };

    auto obj = ObjectiveFunction::make_iv_objective(iv_fn, quotes,
                 WeightingScheme::RelativeError);

    std::vector<Real> x_init(3, 0.0);
    if (cfg.use_de_init) {
        DifferentialEvolution::Config de_cfg;
        de_cfg.population_size = cfg.de_pop_size;
        de_cfg.max_generations = cfg.de_generations;
        de_cfg.seed = cfg.seed;
        de_cfg.lambda_reg = cfg.lambda_reg;
        de_cfg.params_prior = cfg.params_prior;
        de_cfg.early_stop_rmse = cfg.early_stop_rmse;
        auto de_result = DifferentialEvolution::minimize(obj.to_objective_fn(), bounds, de_cfg);
        x_init = de_result.x;
    } else {
        x_init = {0.2, 0.3, -0.1};
    }

    LevenbergMarquardt::Config lm_cfg;
    lm_cfg.max_iterations = cfg.lm_max_iter;
    lm_cfg.ftol = cfg.ftol;
    lm_cfg.xtol = cfg.xtol;
    lm_cfg.lambda_reg = cfg.lambda_reg;
    lm_cfg.params_prior = cfg.params_prior;
    lm_cfg.early_stop_rmse = cfg.early_stop_rmse;
    auto lm_result = LevenbergMarquardt::minimize(obj.to_residual_fn(), x_init, lm_cfg);

    result.params = lm_result.x;
    result.objective_value = lm_result.fx;
    result.n_iterations = lm_result.n_iterations;
    result.converged = lm_result.converged;
    result.message = lm_result.message;
    result.residuals = obj.residuals(lm_result.x);

    VGParams vp = extract_params(result.params);
    if (!check_feller(vp)) {
        result.message += " [WARNING: VG Feller condition violated]";
    }
    return result;
}

// ---------------------------------------------------------------------------
// CEVCalibrator
// Parameter vector x = [sigma, beta]
// 2-parameter calibration (CEV elasticity + volatility scale)
// Objective: sum_i w_i (IV_model_i - IV_market_i)^2
// Pricing: cev_analytic.hpp (noncentral chi2 analytic), IV via bsm_implied_vol
// Constraints: sigma>0, 0<beta<1 (0=normal, 1=lognormal, middle = CEV)
// ---------------------------------------------------------------------------
struct CEVCalibParams {
    Real sigma, beta;
};

class CEVCalibrator : public Calibrator {
public:
    CalibrationResult calibrate(
        const std::vector<MarketQuote>& quotes,
        const CalibConfig& cfg = CalibConfig{}) override;
    std::string name() const override { return "CEVCalibrator"; }

    static std::vector<Bounds> default_bounds() {
        return {{1e-4, 2.0},      // sigma
                {0.01, 0.99}};    // beta (0=正常, 1=对数正态, 中间为 CEV)
    }
    CEVCalibParams extract_params(const std::vector<Real>& x) const {
        return CEVCalibParams{x[0], x[1]};
    }
    void set_market(Real S, Real r, Real q) { S_ = S; r_ = r; q_ = q; }

private:
    Real S_ = 100.0;
    Real r_ = 0.0;
    Real q_ = 0.0;
};

inline CalibrationResult CEVCalibrator::calibrate(
        const std::vector<MarketQuote>& quotes, const CalibConfig& cfg) {
    CalibrationResult result;
    if (quotes.empty()) {
        result.converged = false;
        result.message = "empty quotes";
        return result;
    }

    auto bounds = default_bounds();
    auto iv_fn = [this, &bounds](const std::vector<Real>& x, Real K, Real T) -> Real {
        CEVParams p;
        p.sigma = std::clamp(x[0], bounds[0].lower, bounds[0].upper);
        // cev_analytic 要求 beta < 1; clamp 到安全上界避免数值退化
        p.beta = std::clamp(x[1], bounds[1].lower, 0.999);
        Real price = cev_call_price(S_, K, T, r_, q_, p);
        return bsm_implied_vol(price, S_, K, T, r_, q_, true);
    };

    auto obj = ObjectiveFunction::make_iv_objective(iv_fn, quotes,
                 WeightingScheme::RelativeError);

    std::vector<Real> x_init(2, 0.0);
    if (cfg.use_de_init) {
        DifferentialEvolution::Config de_cfg;
        de_cfg.population_size = cfg.de_pop_size;
        de_cfg.max_generations = cfg.de_generations;
        de_cfg.seed = cfg.seed;
        de_cfg.lambda_reg = cfg.lambda_reg;
        de_cfg.params_prior = cfg.params_prior;
        de_cfg.early_stop_rmse = cfg.early_stop_rmse;
        auto de_result = DifferentialEvolution::minimize(obj.to_objective_fn(), bounds, de_cfg);
        x_init = de_result.x;
    } else {
        x_init = {0.25, 0.5};
    }

    LevenbergMarquardt::Config lm_cfg;
    lm_cfg.max_iterations = cfg.lm_max_iter;
    lm_cfg.ftol = cfg.ftol;
    lm_cfg.xtol = cfg.xtol;
    lm_cfg.lambda_reg = cfg.lambda_reg;
    lm_cfg.params_prior = cfg.params_prior;
    lm_cfg.early_stop_rmse = cfg.early_stop_rmse;
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
