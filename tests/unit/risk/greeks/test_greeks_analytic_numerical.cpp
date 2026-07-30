#include <gtest/gtest.h>
#include <cmath>
#include <functional>

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"
#include "cpphub/risk/greeks/greeks_numerical.hpp"
#include "cpphub/risk/greeks/ad_dual.hpp"

namespace cpphub {
inline namespace v1 {

namespace {

using cpphub::v1::Real;
using cpphub::v1::normal_cdf;
using cpphub::v1::normal_pdf;

constexpr Real S = 100.0;
constexpr Real K = 100.0;
constexpr Real T = 1.0;
constexpr Real r = 0.05;
constexpr Real q = 0.0;
constexpr Real sigma = 0.2;

Real bsm_call_price(Real S_, Real K_, Real T_, Real r_, Real q_, Real sigma_)
{
    Real d1 = (std::log(S_ / K_) + (r_ - q_ + Real(0.5) * sigma_ * sigma_) * T_) / (sigma_ * std::sqrt(T_));
    Real d2 = d1 - sigma_ * std::sqrt(T_);
    return S_ * std::exp(-q_ * T_) * normal_cdf(d1) - K_ * std::exp(-r_ * T_) * normal_cdf(d2);
}

Real bsm_put_price(Real S_, Real K_, Real T_, Real r_, Real q_, Real sigma_)
{
    return bsm_call_price(S_, K_, T_, r_, q_, sigma_) - S_ * std::exp(-q_ * T_) + K_ * std::exp(-r_ * T_);
}

Real bsm_price(Real S_, Real K_, Real T_, Real r_, Real q_, Real sigma_, bool is_call)
{
    if (is_call) return bsm_call_price(S_, K_, T_, r_, q_, sigma_);
    else         return bsm_put_price(S_, K_, T_, r_, q_, sigma_);
}

template<typename T>
T bsm_price_dual(T S_, T K_, T T_, T r_, T q_, T sigma_, bool is_call)
{
    T d1 = (log(S_ / K_) + (r_ - q_ + sigma_ * sigma_ / T(2.0)) * T_) / (sigma_ * sqrt(T_));
    T d2_ = d1 - sigma_ * sqrt(T_);
    T call = S_ * exp(-q_ * T_) * normal_cdf_dual(d1) - K_ * exp(-r_ * T_) * normal_cdf_dual(d2_);
    if (is_call) return call;
    return call - S_ * exp(-q_ * T_) + K_ * exp(-r_ * T_);
}

template<typename T>
T normal_cdf_dual_generic(T x) {
    return T(0.5) * (T(1.0) + erf(x / sqrt(T(2.0))));
}

} // anonymous namespace

TEST(AnalyticGreeks, BSMCallPrice)
{
    auto g = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    EXPECT_NEAR(g.price, 10.4506, 1e-4);
}

TEST(AnalyticGreeks, BSMCallDelta)
{
    auto g = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    EXPECT_NEAR(g.delta, 0.6368, 1e-4);
}

TEST(AnalyticGreeks, BSMCallGamma)
{
    auto g = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    EXPECT_NEAR(g.gamma, 0.0188, 1e-4);
}

TEST(AnalyticGreeks, BSMCallVega)
{
    auto g = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    EXPECT_NEAR(g.vega, 37.5214, 5e-3);
}

TEST(AnalyticGreeks, BSMCallTheta)
{
    auto g = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    EXPECT_NEAR(g.theta, -6.4141, 1e-4);
}

TEST(AnalyticGreeks, BSMCallRho)
{
    auto g = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    EXPECT_NEAR(g.rho, 53.2324, 1e-4);
}

TEST(AnalyticGreeks, BSMPutGreeks)
{
    auto call = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    auto put = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, false);

    Real exp_neg_qT = std::exp(-q * T);
    EXPECT_NEAR(call.delta - put.delta, exp_neg_qT, 1e-10);
    EXPECT_NEAR(call.gamma, put.gamma, 1e-10);
    EXPECT_NEAR(call.vega, put.vega, 1e-10);
    EXPECT_NEAR(call.vanna, put.vanna, 1e-10);
    EXPECT_NEAR(call.vomma, put.vomma, 1e-10);

    Real exp_neg_rT = std::exp(-r * T);
    Real call_price = bsm_call_price(S, K, T, r, q, sigma);
    Real put_price = bsm_put_price(S, K, T, r, q, sigma);
    EXPECT_NEAR(call_price - put_price, S * exp_neg_qT - K * exp_neg_rT, 1e-10);
}

TEST(AnalyticGreeks, HigherOrderGreeks)
{
    auto g = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);

    EXPECT_NEAR(g.vanna, -0.2814, 1e-3);
    EXPECT_NEAR(g.vomma, 9.8501, 1e-3);
    EXPECT_LT(g.vanna, 0);
    EXPECT_GT(g.vomma, 0);
}

TEST(NumericalGreeks, CentralDelta)
{
    auto analytic = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    auto numerical = NumericalGreeksEngine::bsm_european(S, K, T, r, q, sigma, true, bsm_price);
    EXPECT_NEAR(numerical.delta, analytic.delta, 1e-6);
}

TEST(NumericalGreeks, CentralGamma)
{
    auto analytic = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    auto numerical = NumericalGreeksEngine::bsm_european(S, K, T, r, q, sigma, true, bsm_price);
    EXPECT_NEAR(numerical.gamma, analytic.gamma, 5e-4);
}

TEST(NumericalGreeks, CentralVega)
{
    auto analytic = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    auto numerical = NumericalGreeksEngine::bsm_european(S, K, T, r, q, sigma, true, bsm_price);
    EXPECT_NEAR(numerical.vega, analytic.vega, 5e-4);
}

TEST(NumericalGreeks, CentralTheta)
{
    auto analytic = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    auto numerical = NumericalGreeksEngine::bsm_european(S, K, T, r, q, sigma, true, bsm_price);
    EXPECT_NEAR(numerical.theta, analytic.theta, 1e-3);
}

TEST(NumericalGreeks, CentralRho)
{
    auto analytic = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);
    auto numerical = NumericalGreeksEngine::bsm_european(S, K, T, r, q, sigma, true, bsm_price);
    EXPECT_NEAR(numerical.rho, analytic.rho, 5e-4);
}

TEST(NumericalGreeks, BumpSizeSensitivity)
{
    auto analytic = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);

    Real bumps[] = {0.001, 0.01, 0.1};
    Real delta_errors[3];
    Real gamma_errors[3];

    for (int i = 0; i < 3; ++i) {
        auto num = NumericalGreeksEngine::bsm_european(S, K, T, r, q, sigma, true, bsm_price,
            bumps[i], 0.0001, 1.0/365.0, 0.0001);
        delta_errors[i] = std::abs(num.delta - analytic.delta);
        gamma_errors[i] = std::abs(num.gamma - analytic.gamma);
    }

    EXPECT_LT(delta_errors[0], delta_errors[2]);
    EXPECT_LT(gamma_errors[0], gamma_errors[2]);
    EXPECT_LT(delta_errors[1], 1e-6);
    EXPECT_LT(gamma_errors[1], 1e-4);
}

TEST(ADDual, BSMCallDelta)
{
    auto analytic = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);

    dual S_val = S;
    auto price_fn = [&](dual S_dual) -> dual {
        return bsm_price_dual(S_dual, dual(K), dual(T), dual(r), dual(q), dual(sigma), true);
    };
    double delta_ad = derivative(price_fn, wrt(S_val), at(S_val));

    EXPECT_NEAR(delta_ad, analytic.delta, 1e-10);
}

TEST(ADDual, BSMCallVega)
{
    auto analytic = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);

    dual sigma_val = sigma;
    auto price_fn = [&](dual sigma_dual) -> dual {
        return bsm_price_dual(dual(S), dual(K), dual(T), dual(r), dual(q), sigma_dual, true);
    };
    double vega_ad = derivative(price_fn, wrt(sigma_val), at(sigma_val));

    EXPECT_NEAR(vega_ad, analytic.vega, 1e-10);
}

TEST(ADDual, BSMCallRho)
{
    auto analytic = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);

    dual r_val = r;
    auto price_fn = [&](dual r_dual) -> dual {
        return bsm_price_dual(dual(S), dual(K), dual(T), r_dual, dual(q), dual(sigma), true);
    };
    double rho_ad = derivative(price_fn, wrt(r_val), at(r_val));

    EXPECT_NEAR(rho_ad, analytic.rho, 1e-10);
}

TEST(ADDual, BSMPutDelta)
{
    auto analytic = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, false);

    dual S_val = S;
    auto price_fn = [&](dual S_dual) -> dual {
        return bsm_price_dual(S_dual, dual(K), dual(T), dual(r), dual(q), dual(sigma), false);
    };
    double delta_ad = derivative(price_fn, wrt(S_val), at(S_val));

    EXPECT_NEAR(delta_ad, analytic.delta, 1e-10);
}

TEST(ADDual, BSMGammaViaDual)
{
    auto analytic = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, true);

    dual2nd S_val = S;
    auto price_fn = [&](dual2nd S_dual) -> dual2nd {
        dual2nd sigma_dual = sigma;
        dual2nd K_dual = K;
        dual2nd T_dual = T;
        dual2nd r_dual = r;
        dual2nd q_dual = q;
        dual2nd d1 = (log(S_dual / K_dual) + (r_dual - q_dual + sigma_dual * sigma_dual / dual2nd(2.0)) * T_dual) / (sigma_dual * sqrt(T_dual));
        dual2nd d2 = d1 - sigma_dual * sqrt(T_dual);
        dual2nd call = S_dual * exp(-q_dual * T_dual) * normal_cdf_dual_generic(d1)
                     - K_dual * exp(-r_dual * T_dual) * normal_cdf_dual_generic(d2);
        return call;
    };
    double gamma_ad = derivative<2>(price_fn, wrt(S_val), at(S_val));

    EXPECT_NEAR(gamma_ad, analytic.gamma, 1e-8);
}

TEST(ADDual, CrosValidationWithNumerical)
{
    auto numerical = NumericalGreeksEngine::bsm_european(S, K, T, r, q, sigma, true, bsm_price);

    dual S_val = S;
    auto delta_fn = [&](dual S_dual) -> dual {
        return bsm_price_dual(S_dual, dual(K), dual(T), dual(r), dual(q), dual(sigma), true);
    };
    double delta_ad = derivative(delta_fn, wrt(S_val), at(S_val));
    EXPECT_NEAR(delta_ad, numerical.delta, 1e-6);

    dual sigma_val = sigma;
    auto vega_fn = [&](dual sigma_dual) -> dual {
        return bsm_price_dual(dual(S), dual(K), dual(T), dual(r), dual(q), sigma_dual, true);
    };
    double vega_ad = derivative(vega_fn, wrt(sigma_val), at(sigma_val));
    EXPECT_NEAR(vega_ad, numerical.vega, 1e-6);

    dual r_val = r;
    auto rho_fn = [&](dual r_dual) -> dual {
        return bsm_price_dual(dual(S), dual(K), dual(T), r_dual, dual(q), dual(sigma), true);
    };
    double rho_ad = derivative(rho_fn, wrt(r_val), at(r_val));
    EXPECT_NEAR(rho_ad, numerical.rho, 1e-6);
}

} // namespace v1
} // namespace cpphub
