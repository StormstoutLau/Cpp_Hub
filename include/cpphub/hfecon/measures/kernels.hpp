// SOURCE: PHASE5_HFE_SPEC §4.2
//   [BNS 2008] Barndorff-Nielsen, Hansen, Lunde, Shephard,
//              Econometrica 76(6), 1481-1536, doi:10.1111/j.1468-0262.2008.00837.x
// R 对照: highfrequency 1.0.3 KK() — realizedMeasures.cpp L16-74 (CRAN 源码实测)
#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// KernelType: 12 种核函数 (R listAvailableKernels() 实测 2026-08-02)
// =============================================================================
// R 源码 KK(double x, int type) — realizedMeasures.cpp L16-74
// 注意: R 实现与 BNS 2008 Table 1 有差异 (Second/Seventh/Eighth 系数不同)
// 决策: C++ 严格对标 R 源码, 保证 R baseline 数值一致
// =============================================================================

enum class KernelType {
    Rectangular,          // R "rectangular",  type=0: k(x) = 1 (全部, 不归零)
    Bartlett,             // R "Bartlett",     type=1: k(x) = 1 - x
    Second,               // R "Second",       type=2: k(x) = 1 - 2x³
    Epanechnikov,         // R "Epanechnikov", type=3: k(x) = 1 - x²
    Cubic,                // R "Cubic",        type=4: k(x) = 1 - 3x² + 2x³
    Fifth,                // R "Fifth",        type=5: k(x) = 1 - 10x³ + 15x⁴ - 6x⁵
    Sixth,                // R "Sixth",        type=6: k(x) = 1 - 15x⁴ + 24x⁵ - 10x⁶
    Seventh,              // R "Seventh",      type=7: k(x) = 1 - 21x⁵ + 35x⁶ - 15x⁷
    Eighth,               // R "Eighth",       type=8: k(x) = 1 - 28x⁶ + 48x⁷ - 21x⁸
    Parzen,               // R "Parzen",       type=9: 分段 (Bartlett-Parzen)
    TukeyHanning,         // R "TukeyHanning", type=10: (1 + sin(π/2 - πx))/2
    ModifiedTukeyHanning  // R "ModifiedTukeyHanning", type=11: (1 - sin(π/2 - π(1-x)²))/2
};

// =============================================================================
// kernel_value: 核函数值 k(x)
// =============================================================================
// 行为 (严格对标 R highfrequency 1.0.3 KK() 源码):
//   - Rectangular: 恒返回 1.0 (R KK() case 0, 无支撑限制)
//   - 其他核: |x| > 1 时返回 0 (支撑外归零, R KK() 不检查但 estimator
//             中 x=(h-1)/H ∈ [0,1) 永不超 1, C++ 增加支撑检查用于数学正确性)
//   - |x| ≤ 1: 返回 R KK() 公式值 (即使为负, 如 Second k(1)=-1)
//
// Discovery (R 实现偏离 BNS 2008 论文):
//   1. Second: R = 1-2x³, BNS 2008 Table 1 = 1-x²
//   2. Seventh/Eighth: R 多项式阶数和系数与 BNS 2008 不同
//   3. Second k(1) = -1 < 0 (核函数为负, 数学上不合理但 R 实现如此)
//   4. TukeyHanning: R sin(π/2-πx) = BNS (1+cos(πx))/2 (等价)
//
// 异常: 未知 KernelType 抛 invalid_argument
// =============================================================================
inline Real kernel_value(KernelType type, Real x) {
    // Rectangular: 恒为 1 (R KK() case 0, 无支撑限制)
    if (type == KernelType::Rectangular) {
        return 1.0;
    }

    // 其他核: |x| > 1 时归零 (支撑外)
    const Real ax = std::fabs(x);
    if (ax > 1.0) {
        return 0.0;
    }

    // |x| ≤ 1: 计算 R KK() 公式
    switch (type) {
        case KernelType::Bartlett:
            return 1.0 - ax;
        case KernelType::Second:
            return 1.0 - 2.0 * ax * ax * ax;  // 1 - 2x³
        case KernelType::Epanechnikov:
            return 1.0 - ax * ax;              // 1 - x²
        case KernelType::Cubic:
            return 1.0 - 3.0 * ax * ax + 2.0 * ax * ax * ax;  // 1 - 3x² + 2x³
        case KernelType::Fifth: {
            const Real x2 = ax * ax;
            const Real x3 = x2 * ax;
            const Real x4 = x3 * ax;
            const Real x5 = x4 * ax;
            return 1.0 - 10.0 * x3 + 15.0 * x4 - 6.0 * x5;
        }
        case KernelType::Sixth: {
            const Real x2 = ax * ax;
            const Real x3 = x2 * ax;
            const Real x4 = x3 * ax;
            const Real x5 = x4 * ax;
            const Real x6 = x5 * ax;
            return 1.0 - 15.0 * x4 + 24.0 * x5 - 10.0 * x6;
        }
        case KernelType::Seventh: {
            const Real x2 = ax * ax;
            const Real x3 = x2 * ax;
            const Real x4 = x3 * ax;
            const Real x5 = x4 * ax;
            const Real x6 = x5 * ax;
            const Real x7 = x6 * ax;
            return 1.0 - 21.0 * x5 + 35.0 * x6 - 15.0 * x7;
        }
        case KernelType::Eighth: {
            const Real x2 = ax * ax;
            const Real x3 = x2 * ax;
            const Real x4 = x3 * ax;
            const Real x5 = x4 * ax;
            const Real x6 = x5 * ax;
            const Real x7 = x6 * ax;
            const Real x8 = x7 * ax;
            return 1.0 - 28.0 * x6 + 48.0 * x7 - 21.0 * x8;
        }
        case KernelType::Parzen:
            // R KK() case 9: 分段
            // x > 0.5: 2(1-x)³
            // x ≤ 0.5: 1 - 6x² + 6x³
            if (ax > 0.5) {
                const Real t = 1.0 - ax;
                return 2.0 * t * t * t;
            } else {
                const Real x2 = ax * ax;
                const Real x3 = x2 * ax;
                return 1.0 - 6.0 * x2 + 6.0 * x3;
            }
        case KernelType::TukeyHanning:
            // R KK() case 10: (1 + sin(π/2 - πx))/2
            // 等价于 BNS 2008: (1 + cos(πx))/2
            {
                constexpr Real PI = 3.14159265358979323846;
                return (1.0 + std::sin(PI / 2.0 - PI * ax)) / 2.0;
            }
        case KernelType::ModifiedTukeyHanning:
            // R KK() case 11: (1 - sin(π/2 - π(1-x)²))/2
            {
                constexpr Real PI = 3.14159265358979323846;
                const Real t = 1.0 - ax;
                const Real t2 = t * t;
                return (1.0 - std::sin(PI / 2.0 - PI * t2)) / 2.0;
            }
        default:
            throw std::invalid_argument("kernel_value: unknown KernelType");
    }
}

// =============================================================================
// parse_kernel_type: 从 R 字符串解析 KernelType
// =============================================================================
// 大小写敏感, 与 R listAvailableKernels() 返回值一致
// 异常: 未知字符串抛 invalid_argument
// =============================================================================
inline KernelType parse_kernel_type(const std::string& name) {
    if (name == "rectangular")           return KernelType::Rectangular;
    if (name == "Bartlett")              return KernelType::Bartlett;
    if (name == "Second")                return KernelType::Second;
    if (name == "Epanechnikov")          return KernelType::Epanechnikov;
    if (name == "Cubic")                 return KernelType::Cubic;
    if (name == "Fifth")                 return KernelType::Fifth;
    if (name == "Sixth")                 return KernelType::Sixth;
    if (name == "Seventh")               return KernelType::Seventh;
    if (name == "Eighth")                return KernelType::Eighth;
    if (name == "Parzen")                return KernelType::Parzen;
    if (name == "TukeyHanning")          return KernelType::TukeyHanning;
    if (name == "ModifiedTukeyHanning")  return KernelType::ModifiedTukeyHanning;
    throw std::invalid_argument(
        "parse_kernel_type: unknown kernel name '" + name + "'");
}

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
