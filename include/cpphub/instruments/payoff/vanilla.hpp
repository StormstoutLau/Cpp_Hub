#pragma once
#include "payoff.hpp"
#include <algorithm>
#include <memory>
#include <string>

namespace cpphub {
inline namespace v1 {

class CallPayOff : public PayOff {
    double strike_;
public:
    explicit CallPayOff(double strike) : strike_(strike) {}
    double operator()(double spot) const override {
        return std::max(spot - strike_, 0.0);
    }
    std::unique_ptr<PayOff> clone() const override {
        return std::make_unique<CallPayOff>(*this);
    }
    std::string name() const override { return "Call"; }
};

class PutPayOff : public PayOff {
    double strike_;
public:
    explicit PutPayOff(double strike) : strike_(strike) {}
    double operator()(double spot) const override {
        return std::max(strike_ - spot, 0.0);
    }
    std::unique_ptr<PayOff> clone() const override {
        return std::make_unique<PutPayOff>(*this);
    }
    std::string name() const override { return "Put"; }
};

class DigitalCallPayOff : public PayOff {
    double strike_;
    double payment_;
public:
    explicit DigitalCallPayOff(double strike, double payment = 1.0)
        : strike_(strike), payment_(payment) {}
    double operator()(double spot) const override {
        return spot > strike_ ? payment_ : 0.0;
    }
    std::unique_ptr<PayOff> clone() const override {
        return std::make_unique<DigitalCallPayOff>(*this);
    }
    std::string name() const override { return "DigitalCall"; }
};

class DigitalPutPayOff : public PayOff {
    double strike_;
    double payment_;
public:
    explicit DigitalPutPayOff(double strike, double payment = 1.0)
        : strike_(strike), payment_(payment) {}
    double operator()(double spot) const override {
        return spot < strike_ ? payment_ : 0.0;
    }
    std::unique_ptr<PayOff> clone() const override {
        return std::make_unique<DigitalPutPayOff>(*this);
    }
    std::string name() const override { return "DigitalPut"; }
};

class DoubleDigitalPayOff : public PayOff {
    double lower_level_;
    double upper_level_;
    double payment_;
public:
    DoubleDigitalPayOff(double lower, double upper, double payment = 1.0)
        : lower_level_(lower), upper_level_(upper), payment_(payment) {}
    double operator()(double spot) const override {
        return (spot >= lower_level_ && spot <= upper_level_) ? payment_ : 0.0;
    }
    std::unique_ptr<PayOff> clone() const override {
        return std::make_unique<DoubleDigitalPayOff>(*this);
    }
    std::string name() const override { return "DoubleDigital"; }
};

}  // namespace v1
}  // namespace cpphub
