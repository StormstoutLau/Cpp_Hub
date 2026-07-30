#pragma once
#include "cpphub/core/types.hpp"
#include <functional>

namespace cpphub {
inline namespace v1 {

struct NumericalGreeks {
    Real delta, gamma, vega, theta, rho;
};

class NumericalGreeksEngine {
public:
    using PriceFn = std::function<Real(Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call)>;

    static NumericalGreeks bsm_european(
        Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call,
        const PriceFn& price_fn,
        Real dS = 0.01, Real dSigma = 0.0001, Real dT = 1.0 / 365.0, Real dR = 0.0001)
    {
        NumericalGreeks g;

        g.delta = central_difference(price_fn, S, K, T, r, q, sigma, is_call, dS);

        g.gamma = second_difference(price_fn, S, K, T, r, q, sigma, is_call, dS);

        {
            Real v_up = price_fn(S, K, T, r, q, sigma + dSigma, is_call);
            Real v_down = price_fn(S, K, T, r, q, sigma - dSigma, is_call);
            g.vega = (v_up - v_down) / (Real(2.0) * dSigma);
        }

        {
            Real t_up = price_fn(S, K, T + dT, r, q, sigma, is_call);
            Real t_down = price_fn(S, K, T - dT, r, q, sigma, is_call);
            g.theta = -(t_up - t_down) / (Real(2.0) * dT);
        }

        {
            Real r_up = price_fn(S, K, T, r + dR, q, sigma, is_call);
            Real r_down = price_fn(S, K, T, r - dR, q, sigma, is_call);
            g.rho = (r_up - r_down) / (Real(2.0) * dR);
        }

        return g;
    }

    static Real forward_difference(const PriceFn& fn,
        Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call, Real dS)
    {
        Real base = fn(S, K, T, r, q, sigma, is_call);
        Real bumped = fn(S + dS, K, T, r, q, sigma, is_call);
        return (bumped - base) / dS;
    }

    static Real central_difference(const PriceFn& fn,
        Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call, Real dS)
    {
        Real up = fn(S + dS, K, T, r, q, sigma, is_call);
        Real down = fn(S - dS, K, T, r, q, sigma, is_call);
        return (up - down) / (Real(2.0) * dS);
    }

    static Real second_difference(const PriceFn& fn,
        Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call, Real dS)
    {
        Real up = fn(S + dS, K, T, r, q, sigma, is_call);
        Real down = fn(S - dS, K, T, r, q, sigma, is_call);
        Real base = fn(S, K, T, r, q, sigma, is_call);
        return (up - Real(2.0) * base + down) / (dS * dS);
    }
};

}  // namespace v1
}  // namespace cpphub
