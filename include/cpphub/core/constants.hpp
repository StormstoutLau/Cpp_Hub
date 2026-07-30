// SOURCE: PHASE1_SPEC §2.4 - 数学常数
#pragma once
#include <numbers>
#include <cmath>
#include "cpphub/core/types.hpp"

// 避免与 <cfloat> 中的 DBL_EPSILON 宏冲突
#ifdef DBL_EPSILON
#undef DBL_EPSILON
#endif

namespace cpphub {
inline namespace v1 {

// C++20 std::numbers 提供编译期常数
inline constexpr Real PI = std::numbers::pi_v<Real>;
inline constexpr Real SQRT_2 = std::numbers::sqrt2_v<Real>;
inline constexpr Real INV_SQRT_2 = Real(1) / SQRT_2;
// std::sqrt 非 constexpr, 用 2/sqrt(2pi) 的等价形式: sqrt(2)*INV_SQRT_2PI
// 或直接使用数值常量 (2.5066282746310005024157652848110...)
inline constexpr Real SQRT_2PI = Real(2.5066282746310005024157652848110);
inline constexpr Real INV_SQRT_2PI = Real(1) / SQRT_2PI;
inline constexpr Real LN_2 = std::numbers::ln2_v<Real>;
inline constexpr Real CPPHUB_DBL_EPSILON = Real(2.220446049250313e-16);
// 兼容旧代码: 保留 DBL_EPSILON 名字但放在命名空间内, 避免宏冲突
inline constexpr Real EPSILON = CPPHUB_DBL_EPSILON;

}  // namespace v1
}  // namespace cpphub
