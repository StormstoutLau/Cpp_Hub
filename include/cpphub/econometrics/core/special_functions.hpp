// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.2 - 共享特殊函数 (消除跨头文件重复定义)
// 用途: 提供 econometrics 命名空间内多个头文件共享的数学常数与特殊函数
//
// 背景: ols.hpp / mle.hpp / hypothesis_tests.hpp 原各自独立定义 kTwoPi / betacf / beta_i,
//   当同一 TU 同时包含多个头文件时 (如 test_integration_m2.cpp) 触发 C2374/C2084 重定义错误.
//   本头文件作为唯一定义点, 消除重复.
//
// Phase 7A 扩展: chi2_sf/gammp/gammq 等分布函数从 hypothesis_tests.hpp 移入,
//   供 residual_diagnostics.hpp 等 ADR-015 方案 B 头文件复用 (不引入 Eigen3 依赖).
//
// 算法锚点: Numerical Recipes (Press et al.) betai/betacf/gammp/gammq
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训)
#pragma once

#include <cmath>       // std::log, std::sqrt, std::lgamma, std::exp, std::abs, std::erfc
#include <limits>      // std::numeric_limits
#include <stdexcept>   // std::domain_error

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

// =============================================================================
// 不完全 Gamma 函数 (Numerical Recipes gammp/gammq)
//   用于 χ² 分布 p 值计算: P(χ² > x) = gammq(df/2, x/2)
//
// 注: gammln 改用 std::lgamma (C++标准, 精度更高), 替代 Lanczos 近似
//   原 Lanczos 近似对小参数 (如 0.5) 有 1e-4 级误差, std::lgamma 达机器精度
// =============================================================================

constexpr Real kGammaEuler = 0.57721566490153286060651209008240243104215933593992;

/// @brief 对数 Gamma 函数 (委托 std::lgamma, C++标准保证精度)
inline Real gammln(Real xx) {
    return std::lgamma(xx);
}

// 不完全 Gamma 函数 P(a,x) = γ(a,x)/Γ(a) (级数展开, 适用于 x < a+1)
inline Real gser(Real a, Real x) {
    const Real EPS = 3e-16;
    const int ITMAX = 300;
    if (x <= 0.0) return 0.0;
    Real gln = gammln(a);
    Real ap = a;
    Real sum = 1.0 / a;
    Real del = sum;
    for (int n = 0; n < ITMAX; ++n) {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (std::abs(del) < std::abs(sum) * EPS) break;
    }
    return sum * std::exp(-x + a * std::log(x) - gln);
}

// 不完全 Gamma 函数 Q(a,x) = 1 - P(a,x) (连分式展开, 适用于 x >= a+1)
inline Real gcf(Real a, Real x) {
    const Real EPS = 3e-16;
    const Real FPMIN = 1e-300;
    const int ITMAX = 300;
    Real gln = gammln(a);
    Real b = x + 1.0 - a;
    Real c = 1.0 / FPMIN;
    Real d = 1.0 / b;
    Real h = d;
    for (int i = 1; i <= ITMAX; ++i) {
        const Real an = -static_cast<Real>(i) * (static_cast<Real>(i) - a);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        const Real del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < EPS) break;
    }
    return std::exp(-x + a * std::log(x) - gln) * h;
}

/// @brief 不完全 Gamma 函数 P(a,x) = γ(a,x)/Γ(a)
inline Real gammp(Real a, Real x) {
    if (x < 0.0 || a <= 0.0) {
        throw std::domain_error("gammp: invalid arguments (x<0 or a<=0)");
    }
    if (x < a + 1.0) {
        return gser(a, x);
    }
    return 1.0 - gcf(a, x);
}

/// @brief 不完全 Gamma 函数 Q(a,x) = 1 - P(a,x) (上尾)
inline Real gammq(Real a, Real x) {
    if (x < 0.0 || a <= 0.0) {
        throw std::domain_error("gammq: invalid arguments (x<0 or a<=0)");
    }
    if (x < a + 1.0) {
        return 1.0 - gser(a, x);
    }
    return gcf(a, x);
}

/// @brief χ² 分布上尾概率 P(χ²_df > x)
///   通用形式: gammq(df/2, x/2)
///   df=1 精确形式: P(χ²(1) > x) = erfc(√(x/2))
///   df=2 精确形式: P(χ²(2) > x) = exp(-x/2)
inline Real chi2_sf(Real df, Real x) {
    if (df <= 0.0) return std::numeric_limits<Real>::quiet_NaN();
    if (x <= 0.0) return 1.0;
    if (df == 1.0) {
        return std::erfc(std::sqrt(0.5 * x));
    }
    if (df == 2.0) {
        return std::exp(-0.5 * x);
    }
    return gammq(0.5 * df, 0.5 * x);
}

/// @brief χ² 分布下尾概率 P(χ²_df ≤ x) = 1 - chi2_sf(df, x)
inline Real chi2_cdf(Real df, Real x) {
    if (df <= 0.0) return std::numeric_limits<Real>::quiet_NaN();
    if (x <= 0.0) return 0.0;
    return 1.0 - chi2_sf(df, x);
}

/// @brief χ² 分布分位函数 (临界值, Wilson-Hilferty 近似)
///   仅用于临界值显示, 不用于精确推断
inline Real chi2_ppf_approx(Real df, Real p) {
    (void)p;
    const Real z = (p > 0.975) ? 2.3263 : 1.6449;
    const Real t = 1.0 - 2.0 / (9.0 * df) + z * std::sqrt(2.0 / (9.0 * df));
    return df * t * t * t;
}

/// @brief F 分布上尾概率 P(F_{df1,df2} > f)
///   排幻觉点: 参数顺序为 (df2/2, df1/2), 不是 (df1/2, df2/2)
inline Real f_sf(Real df1, Real df2, Real f) {
    if (df1 <= 0.0 || df2 <= 0.0) return std::numeric_limits<Real>::quiet_NaN();
    if (f <= 0.0) return 1.0;
    const Real x = df2 / (df2 + df1 * f);
    return beta_i(0.5 * df2, 0.5 * df1, x);
}

}  // namespace detail
}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
