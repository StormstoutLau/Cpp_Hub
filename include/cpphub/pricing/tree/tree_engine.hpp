#pragma once
#include <memory>
#include <variant>
#include <string>
#include <vector>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/instruments/payoff/payoff.hpp"
#include "cpphub/pricing/tree/binomial.hpp"
#include "cpphub/pricing/tree/trinomial.hpp"

namespace cpphub {
inline namespace v1 {

enum class OptionType {
    European,
    American,
    Bermudan
};

class TreeEngine {
public:
    explicit TreeEngine(BinomialTreeEngine engine) : engine_(std::move(engine)) {}
    explicit TreeEngine(TrinomialTreeEngine engine) : engine_(std::move(engine)) {}

    Real price(const PayOff& payoff, OptionType option_type = OptionType::European,
                const std::vector<Size>& bermudan_steps = {}) const {
        return std::visit([&](const auto& eng) -> Real {
            using T = std::decay_t<decltype(eng)>;
            if constexpr (std::is_same_v<T, BinomialTreeEngine>) {
                switch (option_type) {
                    case OptionType::European: return eng.price_european(payoff);
                    case OptionType::American: return eng.price_american(payoff);
                    case OptionType::Bermudan: return eng.price_bermudan(payoff, bermudan_steps);
                }
            } else {
                switch (option_type) {
                    case OptionType::European: return eng.price_european(payoff);
                    case OptionType::American: return eng.price_american(payoff);
                    case OptionType::Bermudan: return eng.price_bermudan(payoff, bermudan_steps);
                }
            }
            return 0.0;
        }, engine_);
    }

    Real delta(const PayOff& payoff, OptionType option_type = OptionType::European) const {
        return std::visit([&](const auto& eng) -> Real {
            using T = std::decay_t<decltype(eng)>;
            if constexpr (std::is_same_v<T, BinomialTreeEngine>) {
                auto g = eng.greeks(payoff, option_type == OptionType::American);
                return g.delta;
            } else {
                auto g = eng.greeks(payoff, option_type == OptionType::American);
                return g.delta;
            }
        }, engine_);
    }

    Real gamma(const PayOff& payoff, OptionType option_type = OptionType::European) const {
        return std::visit([&](const auto& eng) -> Real {
            using T = std::decay_t<decltype(eng)>;
            if constexpr (std::is_same_v<T, BinomialTreeEngine>) {
                auto g = eng.greeks(payoff, option_type == OptionType::American);
                return g.gamma;
            } else {
                auto g = eng.greeks(payoff, option_type == OptionType::American);
                return g.gamma;
            }
        }, engine_);
    }

    Real theta(const PayOff& payoff, OptionType option_type = OptionType::European) const {
        return std::visit([&](const auto& eng) -> Real {
            using T = std::decay_t<decltype(eng)>;
            if constexpr (std::is_same_v<T, BinomialTreeEngine>) {
                auto g = eng.greeks(payoff, option_type == OptionType::American);
                return g.theta;
            } else {
                auto g = eng.greeks(payoff, option_type == OptionType::American);
                return g.theta;
            }
        }, engine_);
    }

    std::string name() const {
        return std::visit([](const auto& eng) { return eng.name(); }, engine_);
    }

private:
    std::variant<BinomialTreeEngine, TrinomialTreeEngine> engine_;
};

}  // namespace v1
}  // namespace cpphub
