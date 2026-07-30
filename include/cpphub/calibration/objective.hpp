#pragma once
// SOURCE: PHASE3_SPEC §4.1 - Objective functions for calibration
// Implemented on main station (MSVC) - 2026-07-31
// Weighting schemes: Price / Vega / RelativeError / Mixed
// residuals[i] = w_i * (model_i - market_i), where w_i depends on scheme
// model_i is computed by per-quote model function (params, K, T) -> Real
#include "cpphub/core/types.hpp"
#include "cpphub/calibration/optimizer.hpp"
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>

namespace cpphub {
inline namespace v1 {

enum class WeightingScheme {
    PriceWeighted,    // w_i = 1 / |price_i|^0.5  (equalize price magnitude)
    VegaWeighted,     // w_i = 1 / max(|vega_i|, eps)  (downweight liquid ATM)
    RelativeError,    // w_i = 1 / |market_i|  (relative error)
    Mixed             // lambda * price + (1-lambda) * iv
};

struct MarketQuote {
    Real strike;
    Real maturity;
    Real market_price;
    Real implied_vol;
    Real vega;  // optional, for vega-weighted
};

// Per-quote model function: given params and a quote, return model prediction
using PerQuoteModelFn = std::function<Real(const std::vector<Real>&, const MarketQuote&)>;

class ObjectiveFunction {
public:
    // Construct from a per-quote model function comparing against market_price.
    ObjectiveFunction(
        PerQuoteModelFn model_fn,
        const std::vector<MarketQuote>& quotes,
        WeightingScheme scheme = WeightingScheme::VegaWeighted,
        Real lambda_price = 1.0);

    // Compatibility constructor: single-value model (ignores per-quote structure)
    ObjectiveFunction(
        std::function<Real(const std::vector<Real>&)> model_fn,
        const std::vector<MarketQuote>& quotes,
        WeightingScheme scheme = WeightingScheme::VegaWeighted,
        Real lambda_price = 1.0);

    Real evaluate(const std::vector<Real>& params) const;
    std::vector<Real> residuals(const std::vector<Real>& params) const;
    Size n_quotes() const { return quotes_.size(); }
    WeightingScheme scheme() const { return scheme_; }

    // Build an IV objective: iv_fn(params, K, T) -> model implied vol
    // Comparison field: quotes_[i].implied_vol
    static ObjectiveFunction make_iv_objective(
        std::function<Real(const std::vector<Real>&, Real, Real)> iv_fn,
        const std::vector<MarketQuote>& quotes,
        WeightingScheme scheme = WeightingScheme::VegaWeighted);

    // Build a price objective: price_fn(params, K, T) -> model price
    // Comparison field: quotes_[i].market_price
    static ObjectiveFunction make_price_objective(
        std::function<Real(const std::vector<Real>&, Real, Real)> price_fn,
        const std::vector<MarketQuote>& quotes,
        WeightingScheme scheme = WeightingScheme::VegaWeighted);

    // Convert to a ResidualFn suitable for LevenbergMarquardt
    ResidualFn to_residual_fn() const;

    // Convert to an ObjectiveFn suitable for NelderMead / DifferentialEvolution
    ObjectiveFn to_objective_fn() const;

private:
    PerQuoteModelFn model_fn_;
    std::vector<MarketQuote> quotes_;
    WeightingScheme scheme_;
    Real lambda_price_;
    std::vector<Real> weights_;
    bool compare_iv_ = false;  // if true, compare model_fn_ output to implied_vol

    void compute_weights();
    Real market_value(const MarketQuote& q) const {
        return compare_iv_ ? q.implied_vol : q.market_price;
    }
};

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

inline void ObjectiveFunction::compute_weights() {
    Size n = quotes_.size();
    weights_.resize(n, 1.0);
    const Real eps = 1e-12;
    for (Size i = 0; i < n; ++i) {
        Real w = 1.0;
        switch (scheme_) {
            case WeightingScheme::PriceWeighted: {
                Real p = std::abs(quotes_[i].market_price);
                w = (p > eps) ? 1.0 / std::sqrt(p) : 1.0 / std::sqrt(eps);
                break;
            }
            case WeightingScheme::VegaWeighted: {
                Real v = std::abs(quotes_[i].vega);
                w = (v > eps) ? 1.0 / v : 1.0 / eps;
                break;
            }
            case WeightingScheme::RelativeError: {
                Real m = std::abs(market_value(quotes_[i]));
                w = (m > eps) ? 1.0 / m : 1.0 / eps;
                break;
            }
            case WeightingScheme::Mixed:
                w = 1.0;
                break;
        }
        weights_[i] = w;
    }
}

inline ObjectiveFunction::ObjectiveFunction(
        PerQuoteModelFn model_fn,
        const std::vector<MarketQuote>& quotes,
        WeightingScheme scheme,
        Real lambda_price)
    : model_fn_(std::move(model_fn)),
      quotes_(quotes),
      scheme_(scheme),
      lambda_price_(lambda_price) {
    if (quotes_.empty()) {
        throw std::invalid_argument("ObjectiveFunction: quotes cannot be empty");
    }
    compute_weights();
}

inline ObjectiveFunction::ObjectiveFunction(
        std::function<Real(const std::vector<Real>&)> model_fn,
        const std::vector<MarketQuote>& quotes,
        WeightingScheme scheme,
        Real lambda_price)
    : ObjectiveFunction(
          [model_fn](const std::vector<Real>& p, const MarketQuote&) -> Real {
              return model_fn(p);
          },
          quotes, scheme, lambda_price) {}

inline Real ObjectiveFunction::evaluate(const std::vector<Real>& params) const {
    Real sse = 0.0;
    for (Size i = 0; i < quotes_.size(); ++i) {
        Real model = model_fn_(params, quotes_[i]);
        Real market = market_value(quotes_[i]);
        Real diff = model - market;
        sse += weights_[i] * diff * diff;
    }
    return 0.5 * sse;
}

inline std::vector<Real> ObjectiveFunction::residuals(const std::vector<Real>& params) const {
    std::vector<Real> r(quotes_.size());
    for (Size i = 0; i < quotes_.size(); ++i) {
        Real model = model_fn_(params, quotes_[i]);
        Real market = market_value(quotes_[i]);
        r[i] = weights_[i] * (model - market);
    }
    return r;
}

inline ObjectiveFunction ObjectiveFunction::make_iv_objective(
        std::function<Real(const std::vector<Real>&, Real, Real)> iv_fn,
        const std::vector<MarketQuote>& quotes,
        WeightingScheme scheme) {
    PerQuoteModelFn wrapped = [iv_fn](const std::vector<Real>& p, const MarketQuote& q) -> Real {
        return iv_fn(p, q.strike, q.maturity);
    };
    ObjectiveFunction obj(wrapped, quotes, scheme, 1.0);
    obj.compare_iv_ = true;
    obj.compute_weights();  // recompute with compare_iv_ flag set
    return obj;
}

inline ObjectiveFunction ObjectiveFunction::make_price_objective(
        std::function<Real(const std::vector<Real>&, Real, Real)> price_fn,
        const std::vector<MarketQuote>& quotes,
        WeightingScheme scheme) {
    PerQuoteModelFn wrapped = [price_fn](const std::vector<Real>& p, const MarketQuote& q) -> Real {
        return price_fn(p, q.strike, q.maturity);
    };
    return ObjectiveFunction(wrapped, quotes, scheme, 1.0);
}

inline ResidualFn ObjectiveFunction::to_residual_fn() const {
    return [this](const std::vector<Real>& params) -> std::vector<Real> {
        return this->residuals(params);
    };
}

inline ObjectiveFn ObjectiveFunction::to_objective_fn() const {
    return [this](const std::vector<Real>& params) -> Real {
        return this->evaluate(params);
    };
}

}  // inline namespace v1
}  // namespace cpphub
