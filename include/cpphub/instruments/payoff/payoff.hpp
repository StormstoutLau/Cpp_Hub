// SOURCE: Joshi Ch.3 - PayOff 抽象基类 + 虚拟构造函数
#pragma once
#include <memory>
#include <string>

namespace cpphub {
inline namespace v1 {

class PayOff {
public:
    virtual ~PayOff() = default;
    virtual double operator()(double spot) const = 0;
    virtual std::unique_ptr<PayOff> clone() const = 0;
    virtual std::string name() const = 0;
};

}  // namespace v1
}  // namespace cpphub
