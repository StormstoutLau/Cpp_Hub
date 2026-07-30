#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include <cmath>

namespace cpphub {
inline namespace v1 {

struct AnalyticGreeks {
    Real price, delta, gamma, vega, theta, rho, vanna, vomma;
};

class AnalyticGreeksEngine {
public:
    static AnalyticGreeks bsm_european(Real S, Real K, Real T, Real r, Real q,
                                        Real sigma, bool is_call)
    {
        Real d1 = (std::log(S / K) + (r - q + Real(0.5) * sigma * sigma) * T) / (sigma * std::sqrt(T));
        Real d2 = d1 - sigma * std::sqrt(T);

        Real n_d1 = normal_pdf(d1);
        Real N_d1 = normal_cdf(d1);
        Real N_d2 = normal_cdf(d2);

        Real exp_neg_qT = std::exp(-q * T);
        Real exp_neg_rT = std::exp(-r * T);

        Real price;
        Real delta, gamma, vega, theta, rho, vanna, vomma;

        Real price_call = S * exp_neg_qT * N_d1 - K * exp_neg_rT * N_d2;
        Real delta_call = exp_neg_qT * N_d1;
        Real gamma_val = n_d1 * exp_neg_qT / (S * sigma * std::sqrt(T));
        Real vega_val = S * exp_neg_qT * n_d1 * std::sqrt(T);
        Real theta_call = -S * n_d1 * sigma * exp_neg_qT / (Real(2.0) * std::sqrt(T))
                          - r * K * exp_neg_rT * N_d2
                          + q * S * exp_neg_qT * N_d1;
        Real rho_call = K * T * exp_neg_rT * N_d2;
        Real vanna_val = -exp_neg_qT * n_d1 * d2 / sigma;
        Real vomma_val = S * exp_neg_qT * n_d1 * std::sqrt(T) * d1 * d2 / sigma;

        if (is_call) {
            price = price_call;
            delta = delta_call;
            theta = theta_call;
            rho = rho_call;
        } else {
            price = price_call - S * exp_neg_qT + K * exp_neg_rT;
            delta = exp_neg_qT * (N_d1 - Real(1.0));
            theta = -S * n_d1 * sigma * exp_neg_qT / (Real(2.0) * std::sqrt(T))
                    + r * K * exp_neg_rT * normal_cdf(-d2)
                    - q * S * exp_neg_qT * normal_cdf(-d1);
            rho = -K * T * exp_neg_rT * normal_cdf(-d2);
        }

        gamma = gamma_val;
        vega = vega_val;
        vanna = vanna_val;
        vomma = vomma_val;

        return {price, delta, gamma, vega, theta, rho, vanna, vomma};
    }

    static Real bsm_delta(Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call)
    {
        Real d1 = (std::log(S / K) + (r - q + Real(0.5) * sigma * sigma) * T) / (sigma * std::sqrt(T));
        Real N_d1 = normal_cdf(d1);
        Real exp_neg_qT = std::exp(-q * T);
        if (is_call) return exp_neg_qT * N_d1;
        else         return exp_neg_qT * (N_d1 - Real(1.0));
    }

    static Real bsm_gamma(Real S, Real K, Real T, Real r, Real q, Real sigma)
    {
        Real d1 = (std::log(S / K) + (r - q + Real(0.5) * sigma * sigma) * T) / (sigma * std::sqrt(T));
        return normal_pdf(d1) * std::exp(-q * T) / (S * sigma * std::sqrt(T));
    }

    static Real bsm_vega(Real S, Real K, Real T, Real r, Real q, Real sigma)
    {
        Real d1 = (std::log(S / K) + (r - q + Real(0.5) * sigma * sigma) * T) / (sigma * std::sqrt(T));
        return S * std::exp(-q * T) * normal_pdf(d1) * std::sqrt(T);
    }

    static Real bsm_theta(Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call)
    {
        Real d1 = (std::log(S / K) + (r - q + Real(0.5) * sigma * sigma) * T) / (sigma * std::sqrt(T));
        Real d2 = d1 - sigma * std::sqrt(T);
        Real n_d1 = normal_pdf(d1);
        Real exp_neg_qT = std::exp(-q * T);
        Real exp_neg_rT = std::exp(-r * T);
        Real common = -S * n_d1 * sigma * exp_neg_qT / (Real(2.0) * std::sqrt(T));
        if (is_call) {
            return common - r * K * exp_neg_rT * normal_cdf(d2) + q * S * exp_neg_qT * normal_cdf(d1);
        } else {
            return common + r * K * exp_neg_rT * normal_cdf(-d2) - q * S * exp_neg_qT * normal_cdf(-d1);
        }
    }

    static Real bsm_rho(Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call)
    {
        Real d1 = (std::log(S / K) + (r - q + Real(0.5) * sigma * sigma) * T) / (sigma * std::sqrt(T));
        Real d2 = d1 - sigma * std::sqrt(T);
        Real exp_neg_rT = std::exp(-r * T);
        if (is_call) return K * T * exp_neg_rT * normal_cdf(d2);
        else         return -K * T * exp_neg_rT * normal_cdf(-d2);
    }

    static Real bsm_vanna(Real S, Real K, Real T, Real r, Real q, Real sigma)
    {
        Real d1 = (std::log(S / K) + (r - q + Real(0.5) * sigma * sigma) * T) / (sigma * std::sqrt(T));
        Real d2 = d1 - sigma * std::sqrt(T);
        return -std::exp(-q * T) * normal_pdf(d1) * d2 / sigma;
    }

    static Real bsm_vomma(Real S, Real K, Real T, Real r, Real q, Real sigma)
    {
        Real d1 = (std::log(S / K) + (r - q + Real(0.5) * sigma * sigma) * T) / (sigma * std::sqrt(T));
        Real d2 = d1 - sigma * std::sqrt(T);
        return S * std::exp(-q * T) * normal_pdf(d1) * std::sqrt(T) * d1 * d2 / sigma;
    }
};

}  // namespace v1
}  // namespace cpphub
