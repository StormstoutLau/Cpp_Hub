#pragma once
// SOURCE: PHASE3_SPEC §4.1 - Objective functions for calibration
// TODO: Implement by A station agent
#include "cpphub/core/types.hpp"
#include "cpphub/calibration/optimizer.hpp"
#include <vector>
#include <functional>

namespace cpphub {
inline namespace v1 {

enum class WeightingScheme {
    PriceWeighted,
    VegaWeighted,
    RelativeError,
    Mixed
};

struct MarketQuote {
    Real strike;
    Real maturity;
    Real market_price;
    Real implied_vol;
    Real vega;  // optional, for vega-weighted
};

class ObjectiveFunction {
public:
    ObjectiveFunction(
        std::function<Real(const std::vector<Real>&)> model_fn,
        const std::vector<MarketQuote>& quotes,
        WeightingScheme scheme = WeightingScheme::VegaWeighted,
        Real lambda_price = 1.0);

    Real evaluate(const std::vector<Real>& params) const;
    std::vector<Real> residuals(const std::vector<Real>& params) const;
    Size n_quotes() const { return quotes_.size(); }

    static ObjectiveFunction make_iv_objective(
        std::function<Real(const std::vector<Real>&, Real, Real)> iv_fn,
        const std::vector<MarketQuote>& quotes,
        WeightingScheme scheme = WeightingScheme::VegaWeighted);

private:
    std::function<Real(const std::vector<Real>&)> model_fn_;
    std::vector<MarketQuote> quotes_;
    WeightingScheme scheme_;
    Real lambda_price_;
};

}  // inline namespace v1
}  // namespace cpphub
