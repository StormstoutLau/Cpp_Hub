// SOURCE: Joshi Ch.3 - PayOff 抽象基类 + 虚拟构造函数
#pragma once
#include <memory>
#include <string>

namespace cpphub {
inline namespace v1 {

// RISK-031 (ADR-021 P3): payoff 类型识别枚举 — 替代 name() 字符串分派
// (QuantLib QuEP 2 先例: 边界策略与 payoff 解耦; name() 可被自定义派生类
// 冒名劫持, enum 不能)。Unknown 为基类默认 → 自定义 payoff fail-fast。
enum class PayOffKind {
    Call,
    Put,
    DigitalCall,
    DigitalPut,
    DoubleDigital,
    Unknown
};

class PayOff {
public:
    virtual ~PayOff() = default;
    virtual double operator()(double spot) const = 0;
    virtual std::unique_ptr<PayOff> clone() const = 0;
    virtual std::string name() const = 0;
    // 非纯虚 (默认 Unknown): 现有/用户派生类源码兼容; 内置 vanilla 派生 override。
    // PDE 引擎据此选择边界条件 (price_barrier 对 Unknown fail-fast throw)。
    virtual PayOffKind kind() const { return PayOffKind::Unknown; }
};

}  // namespace v1
}  // namespace cpphub
