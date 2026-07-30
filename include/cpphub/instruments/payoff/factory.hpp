#pragma once
#include "payoff.hpp"
#include "cpphub/core/error.hpp"
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace cpphub {
inline namespace v1 {

class PayOffFactory {
public:
    using Creator = std::function<std::unique_ptr<PayOff>(const std::map<std::string, double>&)>;

private:
    std::map<std::string, Creator> registry_;
    PayOffFactory() = default;
public:
    PayOffFactory(const PayOffFactory&) = delete;
    PayOffFactory& operator=(const PayOffFactory&) = delete;

    static PayOffFactory& instance() {
        static PayOffFactory factory;
        return factory;
    }

    void register_payoff(const std::string& name, Creator creator) {
        registry_[name] = std::move(creator);
    }

    std::unique_ptr<PayOff> create(const std::string& name,
                                    const std::map<std::string, double>& params = {}) const {
        auto it = registry_.find(name);
        if (it == registry_.end()) {
            throw CppHubException("Unknown PayOff type: " + name,
                                   ErrorCode::InvalidArgument);
        }
        return it->second(params);
    }

    std::vector<std::string> list() const {
        std::vector<std::string> names;
        for (const auto& [name, _] : registry_) names.push_back(name);
        return names;
    }

    bool is_registered(const std::string& name) const {
        return registry_.count(name) > 0;
    }
};

template <typename T>
class PayOffHelper {
public:
    PayOffHelper(const std::string& name) {
        PayOffFactory::instance().register_payoff(name,
            [](const std::map<std::string, double>& params) -> std::unique_ptr<PayOff> {
                double strike = 0.0;
                auto it = params.find("strike");
                if (it != params.end()) strike = it->second;
                return std::make_unique<T>(strike);
            });
    }
};

#define REGISTER_PAYOFF(name, classname) \
    static ::cpphub::v1::PayOffHelper<classname> \
    classname##_reg(#name)

}  // namespace v1
}  // namespace cpphub
