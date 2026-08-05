// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.2 - 共享特殊函数 (消除跨头文件重复定义)
// 用途: 提供 econometrics 命名空间内多个头文件共享的数学常数与特殊函数
//
// 背景: ols.hpp / mle.hpp / hypothesis_tests.hpp 原各自独立定义 kTwoPi / betacf / beta_i,
//   当同一 TU 同时包含多个头文件时 (如 test_integration_m2.cpp) 触发 C2374/C2084 重定义错误.
//   本头文件作为唯一定义点, 消除重复.
//
// 算法锚点: Numerical Recipes (Press et al.) betai/betacf (Lentz 连分式)
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训)
#pragma once

#include <cmath>       // std::log, std::sqrt, std::lgamma, std::exp, std::abs
#include <limits>      // std::numeric_limits

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {
namespace detail {

// =============================================================================
// 数学常数
// =============================================================================
constexpr Real kTwoPi = 6.2831853071795864769252867665590057683943387987502;
constexpr Real kLogTwoPi = 1.8378770664093454835606594728112352797227949472756;

// =============================================================================
// 正则化不完全贝塔函数 I_x(a,b) = B_x(a,b) / B(a,b)
//   算法: 连分式 (Lentz 方法, Numerical Recipes betacf), 不依赖外部统计库
//   用途:
//     - t 分布双侧 p 值: p = I_x(df/2, 1/2), x = df/(df + t²)
//     - F 分布上尾概率: P(F > f) = I_x(df2/2, df1/2), x = df2/(df2 + df1·f)
// =============================================================================

/// @brief 连分式展开 (Lentz 方法), 用于 I_x(a,b)
inline Real betacf(Real a, Real b, Real x) {
    const int MAXIT = 300;
    const Real EPS = 3e-16;
    const Real FPMIN = 1e-300;
    const Real qab = a + b;
    const Real qap = a + 1.0;
    const Real qam = a - 1.0;
    Real c = 1.0;
    Real d = 1.0 - qab * x / qap;
    if (std::abs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    Real h = d;
    for (int m = 1; m <= MAXIT; ++m) {
        const Real mf = static_cast<Real>(m);
        const int m2 = 2 * m;
        // 偶步
        Real aa = mf * (b - mf) * x / ((qam + static_cast<Real>(m2)) * (a + static_cast<Real>(m2)));
        d = 1.0 + aa * d;
        if (std::abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        // 奇步
        aa = -(a + mf) * (qab + mf) * x / ((a + static_cast<Real>(m2)) * (qap + static_cast<Real>(m2)));
        d = 1.0 + aa * d;
        if (std::abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        const Real del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < EPS) break;
    }
    return h;
}

/// @brief 正则化不完全贝塔 I_x(a,b) = B_x(a,b)/B(a,b)
///   注: 直接使用 std::lgamma (C++标准, 机器精度), 不经 gammln 包装
inline Real beta_i(Real a, Real b, Real x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    const Real lbeta = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b);
    const Real bt = std::exp(lbeta + a * std::log(x) + b * std::log(1.0 - x));
    // 选择收敛更快的方向
    if (x < (a + 1.0) / (a + b + 2.0)) {
        return bt * betacf(a, b, x) / a;
    }
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

}  // namespace detail
}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
