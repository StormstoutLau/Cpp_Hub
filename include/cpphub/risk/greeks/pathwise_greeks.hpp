#pragma once
// SOURCE: PHASE3_SPEC §2.2 - Pathwise Greeks for European options
// Algorithm: d Price / d θ = E[discount * d Payoff / d θ]
//   For GBM: S_T = S * exp((r-q-σ²/2)T + σ√T Z)
//   ∂S_T/∂S = S_T / S
//   ∂S_T/∂σ = S_T * (√T Z - σT) = S_T * (W_T - σT)
// Applicability: smooth payoff (European call/put, arithmetic Asian).
// For discontinuous payoff (digital/barrier), use LRGreeks.
// Ref: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" §7
// TDD: validated vs AnalyticGreeksEngine (1e-3 relative tolerance, 200k paths)
// Implementation: inline xorshift64* RNG to avoid MSVC <random> ICE (RISK-014 sibling)
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include <cmath>

namespace cpphub {
inline namespace v1 {

struct PathwiseGreeks {
    Real price;
    Real delta;     // d Price / d S
    Real vega;      // d Price / d σ
    Real gamma;     // d² Price / d S² (pathwise 二阶,需要 S_T > K 区域的二阶展开;此处用 0 表示不支持)
};

class PathwiseGreeksEngine {
public:
    // BSM European option via Pathwise estimator
    // Single-step GBM (terminal-only simulation, no path dependence)
    static PathwiseGreeks bsm_european(
        Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call,
        Size n_paths = 100000, uint64_t seed = 42) {

        // Inline xorshift64* RNG (避免 <random> MSVC ICE)
        uint64_t state = seed ? seed : 0x9E3779B97F4A7C15ULL;
        auto next_uniform = [&state]() -> Real {
            uint64_t x = state;
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            state = x;
            uint64_t r = x * 0x2545F4914F6CDD1DULL;
            // 取高 53 位作为 [0,1) double
            return static_cast<Real>(r >> 11) * (1.0 / 9007199254740992.0);
        };
        auto next_normal = [&]() -> Real {
            // Box-Muller
            Real u1 = next_uniform();
            Real u2 = next_uniform();
            if (u1 < 1e-300) u1 = 1e-300;
            Real mag = std::sqrt(-2.0 * std::log(u1));
            Real ang = 2.0 * PI * u2;
            // 只用第一个,第二个丢弃 (简化,可优化用掉)
            return mag * std::cos(ang);
        };

        Real drift = (r - q - 0.5 * sigma * sigma) * T;
        Real vol_sqrt_T = sigma * std::sqrt(T);
        Real disc = std::exp(-r * T);

        Real sum_price = 0.0;
        Real sum_delta = 0.0;
        Real sum_vega  = 0.0;

        Real sqrt_T = std::sqrt(T);

        for (Size i = 0; i < n_paths; ++i) {
            Real Z = next_normal();
            Real W_T = vol_sqrt_T * Z;  // Brownian at T = σ√T Z
            Real ST = S * std::exp(drift + W_T);

            // Payoff
            Real payoff = is_call ? (ST - K) : (K - ST);
            if (payoff < 0.0) payoff = 0.0;

            // Pathwise derivatives (in-the-money 区域)
            bool itm = is_call ? (ST > K) : (ST < K);
            if (itm) {
                // d Payoff / d S = sign * dS_T/dS = sign * S_T / S
                Real sign = is_call ? 1.0 : -1.0;
                Real dPayoff_dS = sign * ST / S;
                // d Payoff / d σ = sign * dS_T/dσ
                //   S_T = S * exp((r-q-σ²/2)T + σ√T Z)
                //   dS_T/dσ = S_T * (√T Z - σT)   [NOT σ√T Z - σT]
                Real dPayoff_dsigma = sign * ST * (sqrt_T * Z - sigma * T);

                sum_delta += dPayoff_dS;
                sum_vega  += dPayoff_dsigma;
            }
            sum_price += payoff;
        }

        Real inv_n = 1.0 / static_cast<Real>(n_paths);
        return PathwiseGreeks{
            .price = disc * sum_price * inv_n,
            .delta = disc * sum_delta * inv_n,
            .vega  = disc * sum_vega  * inv_n,
            .gamma = 0.0  // pathwise 二阶不直接支持
        };
    }
};

}  // inline namespace v1
}  // namespace cpphub
