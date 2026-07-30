// SOURCE: BUILD_PLAN §2.2 - PricingEngine Strategy 基类
#pragma once
#include <string>

namespace cpphub {
inline namespace v1 {

struct Greeks {
    double delta{0};
    double gamma{0};
    double vega{0};
    double theta{0};
    double rho{0};
};

class PricingEngine {
public:
    virtual ~PricingEngine() = default;
    virtual double price(/* VanillaOption, 待后续定义 */) const = 0;
    virtual Greeks greeks(/* VanillaOption */) const = 0;
    virtual std::string name() const = 0;
};

}  // namespace v1
}  // namespace cpphub
