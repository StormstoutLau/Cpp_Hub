#pragma once
#include "payoff.hpp"
#include <memory>
#include <stdexcept>

namespace cpphub {
inline namespace v1 {

class PayOffBridge {
    std::unique_ptr<PayOff> ptr_;
public:
    PayOffBridge() = default;

    explicit PayOffBridge(std::unique_ptr<PayOff> ptr) : ptr_(std::move(ptr)) {}

    PayOffBridge(const PayOff& p) : ptr_(p.clone()) {}

    PayOffBridge(const PayOffBridge& other)
        : ptr_(other.ptr_ ? other.ptr_->clone() : nullptr) {}

    PayOffBridge(PayOffBridge&& other) noexcept = default;

    PayOffBridge& operator=(const PayOffBridge& other) {
        if (this != &other) {
            ptr_ = other.ptr_ ? other.ptr_->clone() : nullptr;
        }
        return *this;
    }

    PayOffBridge& operator=(PayOffBridge&& other) noexcept = default;

    ~PayOffBridge() = default;

    double operator()(double spot) const {
        if (!ptr_) {
            throw std::runtime_error("PayOffBridge: empty bridge");
        }
        return (*ptr_)(spot);
    }

    const PayOff& get() const { return *ptr_; }
    PayOff& get() { return *ptr_; }

    explicit operator bool() const noexcept { return static_cast<bool>(ptr_); }
};

}  // namespace v1
}  // namespace cpphub
