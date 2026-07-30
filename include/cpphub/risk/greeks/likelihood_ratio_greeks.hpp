#pragma once
// SOURCE: PHASE3_SPEC §2.2 - Likelihood Ratio (LR) Greeks
// Algorithm: d Price / d θ = E[discount * Payoff(S_T) * (∂ log p / ∂θ)]
//   For GBM: S_T = S * exp((r-q-σ²/2)T + σ√T Z), Z ~ N(0,1)
//   Score functions (Glasserman §7, Eq 7.34-7.35):
//     ∂ log p / ∂S = Z / (S σ √T)
//     ∂ log p / ∂σ = (Z² - 1)/σ - √T Z
//   delta = e^{-rT} * E[Payoff(S_T) * Z / (S σ √T)]
//   vega  = e^{-rT} * E[Payoff(S_T) * ((Z²-1)/σ - √T Z)]
//
// Applicability:
//   - Discontinuous payoffs (digital, barrier, lookback with kink) — PRIMARY use case
//   - Smooth payoffs also work, but variance is much larger than Pathwise
//   - Gamma not directly supported (would need 2nd-order score, O(1/T) variance)
//
// vs Pathwise: LR has higher variance but works for ANY payoff (no smoothness
// requirement). LR is unbiased where Pathwise returns 0 (deep OTM digital).
// vs AAD: AAD is exact but requires differentiable payoff; LR works for
//   indicator payoffs where AAD's max() kink causes biased gradient.
//
// Ref: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" §7.2-7.3
// TDD: validated vs AnalyticGreeksEngine for BSM, vs closed-form for digital
// Implementation: inline xorshift64* RNG to avoid MSVC <random> ICE (RISK-014 sibling)
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include <cmath>

namespace cpphub {
inline namespace v1 {

struct LRGreeks {
    Real price;
    Real delta;     // d Price / d S
    Real vega;      // d Price / d σ
    Real gamma;     // not directly supported by LR (set to 0)
};

class LRGreeksEngine {
public:
    // BSM European option (smooth payoff) via LR estimator
    // Higher variance than Pathwise; use this for cross-validation only.
    static LRGreeks bsm_european(
        Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call,
        Size n_paths = 100000, uint64_t seed = 42) {

        auto rng = make_rng(seed);
        Real drift = (r - q - 0.5 * sigma * sigma) * T;
        Real vol_sqrt_T = sigma * std::sqrt(T);
        Real disc = std::exp(-r * T);

        // Pre-compute score function denominators
        Real inv_S_sigma_sqrtT = 1.0 / (S * sigma * std::sqrt(T));
        Real inv_sigma = 1.0 / sigma;

        Real sum_price = 0.0;
        Real sum_delta = 0.0;
        Real sum_vega  = 0.0;

        for (Size i = 0; i < n_paths; ++i) {
            Real Z = next_normal(rng);
            Real ST = S * std::exp(drift + vol_sqrt_T * Z);

            // Smooth payoff: max(ST-K, 0) or max(K-ST, 0)
            Real payoff = is_call ? (ST - K) : (K - ST);
            if (payoff < 0.0) payoff = 0.0;

            // LR scores
            Real score_S = Z * inv_S_sigma_sqrtT;
            Real score_v = (Z * Z - 1.0) * inv_sigma - std::sqrt(T) * Z;

            sum_price += payoff;
            sum_delta += payoff * score_S;
            sum_vega  += payoff * score_v;
        }

        Real inv_n = 1.0 / static_cast<Real>(n_paths);
        return LRGreeks{
            .price = disc * sum_price * inv_n,
            .delta = disc * sum_delta * inv_n,
            .vega  = disc * sum_vega  * inv_n,
            .gamma = 0.0
        };
    }

    // Digital (binary cash-or-nothing) European option via LR estimator
    // This is the PRIMARY use case: Pathwise fails (payoff indicator has
    // zero derivative a.e.), AAD biased (max(·,0) kink), LR unbiased.
    // Payoff = 1{ST > K} for call, 1{ST < K} for put.
    static LRGreeks digital_european(
        Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call,
        Size n_paths = 100000, uint64_t seed = 42) {

        auto rng = make_rng(seed);
        Real drift = (r - q - 0.5 * sigma * sigma) * T;
        Real vol_sqrt_T = sigma * std::sqrt(T);
        Real disc = std::exp(-r * T);

        Real inv_S_sigma_sqrtT = 1.0 / (S * sigma * std::sqrt(T));
        Real inv_sigma = 1.0 / sigma;

        Real sum_price = 0.0;
        Real sum_delta = 0.0;
        Real sum_vega  = 0.0;

        for (Size i = 0; i < n_paths; ++i) {
            Real Z = next_normal(rng);
            Real ST = S * std::exp(drift + vol_sqrt_T * Z);

            // Indicator payoff: 1 unit cash if ITM at expiry
            Real payoff = is_call ? (ST > K ? 1.0 : 0.0)
                                  : (ST < K ? 1.0 : 0.0);

            // LR scores (identical to smooth case — score depends on density, not payoff)
            Real score_S = Z * inv_S_sigma_sqrtT;
            Real score_v = (Z * Z - 1.0) * inv_sigma - std::sqrt(T) * Z;

            sum_price += payoff;
            sum_delta += payoff * score_S;
            sum_vega  += payoff * score_v;
        }

        Real inv_n = 1.0 / static_cast<Real>(n_paths);
        return LRGreeks{
            .price = disc * sum_price * inv_n,
            .delta = disc * sum_delta * inv_n,
            .vega  = disc * sum_vega  * inv_n,
            .gamma = 0.0
        };
    }

private:
    // xorshift64* state container
    struct RngState { uint64_t s; };

    static RngState make_rng(uint64_t seed) {
        return RngState{seed ? seed : 0x9E3779B97F4A7C15ULL};
    }

    static Real next_uniform(RngState& rng) {
        uint64_t x = rng.s;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        rng.s = x;
        uint64_t r = x * 0x2545F4914F6CDD1DULL;
        return static_cast<Real>(r >> 11) * (1.0 / 9007199254740992.0);
    }

    static Real next_normal(RngState& rng) {
        // Box-Muller (uses one uniform, discards second — simpler, slightly wasteful)
        Real u1 = next_uniform(rng);
        Real u2 = next_uniform(rng);
        if (u1 < 1e-300) u1 = 1e-300;
        Real mag = std::sqrt(-2.0 * std::log(u1));
        Real ang = 2.0 * PI * u2;
        return mag * std::cos(ang);
    }
};

}  // inline namespace v1
}  // namespace cpphub
