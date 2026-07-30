#pragma once
#include "payoff.hpp"
#include <algorithm>
#include <memory>
#include <string>

namespace cpphub {
inline namespace v1 {

enum class AsianType { Arithmetic, Geometric };
enum class AsianPayoffType { Call, Put };

class AsianPayOff : public PayOff {
    double strike_;
    AsianType avg_type_;
    AsianPayoffType payoff_type_;
public:
    AsianPayOff(double strike, AsianType avg_type = AsianType::Arithmetic,
                AsianPayoffType payoff_type = AsianPayoffType::Call)
        : strike_(strike), avg_type_(avg_type), payoff_type_(payoff_type) {}

    double operator()(double spot) const override {
        if (payoff_type_ == AsianPayoffType::Call) {
            return std::max(spot - strike_, 0.0);
        } else {
            return std::max(strike_ - spot, 0.0);
        }
    }
    std::unique_ptr<PayOff> clone() const override {
        return std::make_unique<AsianPayOff>(*this);
    }
    std::string name() const override {
        std::string s = "Asian";
        s += (avg_type_ == AsianType::Arithmetic ? "Arith" : "Geom");
        s += (payoff_type_ == AsianPayoffType::Call ? "Call" : "Put");
        return s;
    }
    AsianType average_type() const { return avg_type_; }
};

enum class BarrierType {
    UpIn,
    DownIn,
    UpOut,
    DownOut
};

class BarrierPayOff : public PayOff {
    double strike_;
    double barrier_;
    BarrierType type_;
    std::unique_ptr<PayOff> inner_;

    bool is_triggered(double spot) const {
        switch (type_) {
            case BarrierType::UpIn:
            case BarrierType::UpOut:
                return spot >= barrier_;
            case BarrierType::DownIn:
            case BarrierType::DownOut:
                return spot <= barrier_;
        }
        return false;
    }
public:
    BarrierPayOff(double strike, double barrier, BarrierType type,
                  std::unique_ptr<PayOff> inner)
        : strike_(strike), barrier_(barrier), type_(type),
          inner_(std::move(inner)) {}

    double operator()(double spot) const override {
        bool triggered = is_triggered(spot);
        switch (type_) {
            case BarrierType::UpIn:
            case BarrierType::DownIn:
                return triggered ? (*inner_)(spot) : 0.0;
            case BarrierType::UpOut:
            case BarrierType::DownOut:
                return triggered ? 0.0 : (*inner_)(spot);
        }
        return 0.0;
    }
    std::unique_ptr<PayOff> clone() const override {
        return std::make_unique<BarrierPayOff>(strike_, barrier_, type_,
                                                inner_->clone());
    }
    std::string name() const override {
        std::string s = "Barrier";
        switch (type_) {
            case BarrierType::UpIn: s += "UpIn"; break;
            case BarrierType::DownIn: s += "DownIn"; break;
            case BarrierType::UpOut: s += "UpOut"; break;
            case BarrierType::DownOut: s += "DownOut"; break;
        }
        return s;
    }
    double barrier() const { return barrier_; }
    BarrierType barrier_type() const { return type_; }
};

enum class LookbackType { CallOnMin, CallOnMax, PutOnMin, PutOnMax };

class LookbackPayOff : public PayOff {
    double strike_;
    LookbackType type_;
public:
    LookbackPayOff(double strike, LookbackType type)
        : strike_(strike), type_(type) {}

    double operator()(double spot) const override {
        switch (type_) {
            case LookbackType::CallOnMin:
                return std::max(spot - strike_, 0.0);
            case LookbackType::CallOnMax:
                return std::max(spot - strike_, 0.0);
            case LookbackType::PutOnMin:
                return std::max(strike_ - spot, 0.0);
            case LookbackType::PutOnMax:
                return std::max(strike_ - spot, 0.0);
        }
        return 0.0;
    }
    std::unique_ptr<PayOff> clone() const override {
        return std::make_unique<LookbackPayOff>(*this);
    }
    std::string name() const override {
        std::string s = "Lookback";
        switch (type_) {
            case LookbackType::CallOnMin: s += "CallOnMin"; break;
            case LookbackType::CallOnMax: s += "CallOnMax"; break;
            case LookbackType::PutOnMin: s += "PutOnMin"; break;
            case LookbackType::PutOnMax: s += "PutOnMax"; break;
        }
        return s;
    }
};

}  // namespace v1
}  // namespace cpphub
