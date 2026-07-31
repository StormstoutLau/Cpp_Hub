#pragma once
// SOURCE: Li, Ding (2015) "Solving nonlinear fractional differential equations
//         by Adams-Bashforth-Moulton method" (Numerical Algorithms)
// SOURCE: Diethelm, Ford, Freed (2002) "A Predictor-Corrector Approach for the
//         Numerical Solution of Fractional Differential Equations"
//
// 分数阶 ODE 求解器: 求解 Caputo 分数阶初值问题
//   D^α y(t) = f(t, y(t)),  y(0) = y0,  α ∈ (0, 1]
//
// Caputo 分数阶导数 (α ∈ (0,1)):
//   D^α f(t) = (1/Γ(1-α)) ∫₀ᵗ f'(s) / (t-s)^α ds
//
// 等价 Volterra 积分方程:
//   y(t) = y0 + (1/Γ(α)) ∫₀ᵗ (t-s)^{α-1} f(s, y(s)) ds
//
// 在均匀网格 t_j = j·Δt 上使用 Adams-Bashforth-Moulton 预测校正 (产品积分法):
//
// 预测步 (Adams-Bashforth, 左矩形核积分):
//   y^P_{n+1} = y0 + Σⱼ₌₀ⁿ bⱼ f(tⱼ, yⱼ),  bⱼ = (Δt^α/Γ(α+1))·[(n+1-j)^α - (n-j)^α]
//
// 校正步 (Adams-Moulton, 梯形核积分):
//   y_{n+1} = y0 + (Δt^α/Γ(α+2))·[ f(t_{n+1}, y^P_{n+1})
//              + a₀ f(t₀,y₀) + Σⱼ₌₁ⁿ aⱼ f(tⱼ,yⱼ) ]
//   a₀ = n^{α+1} - (n-α)(n+1)^α
//   aⱼ = (n-j+2)^{α+1} + (n-j)^{α+1} - 2(n-j+1)^{α+1}   (1 ≤ j ≤ n)
//
// 该方法的阶数为 min(2, 1+α) (对 t^α 类奇异初值), 对常数右端精确.
// 当 α = 1 时精确退化为标准 Adams-Bashforth-Moulton (梯形法则, 2 阶).
#include "cpphub/core/types.hpp"
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace cpphub {
inline namespace v1 {

using ComplexODEFn = std::function<Complex(Real, Complex)>;

namespace detail {

// Γ(x), x > 0
inline Real gamma_function(Real x) { return std::tgamma(x); }

// 预测核权重: b_j = (n+1-j)^α - (n-j)^α
inline Real adams_predictor_weight(Real alpha, Size n, Size j) {
    Real upper = static_cast<Real>(n + 1 - j);
    Real lower = static_cast<Real>(n - j);
    return std::pow(upper, alpha) - std::pow(lower, alpha);
}

// 校正核权重 (内部点, 1 ≤ j ≤ n):
// a_j = (n-j+2)^{α+1} + (n-j)^{α+1} - 2(n-j+1)^{α+1}
inline Real adams_corrector_weight(Real alpha, Size n, Size j) {
    Real a = static_cast<Real>(n - j + 2);
    Real b = static_cast<Real>(n - j);
    Real c = static_cast<Real>(n - j + 1);
    Real ap = alpha + Real(1);
    return std::pow(a, ap) + std::pow(b, ap) - Real(2) * std::pow(c, ap);
}

}  // namespace detail

// 求解复值分数阶 ODE, 返回 [y(t_0), y(t_1), ..., y(t_N)]
// t_grid: 均匀网格 t_j = j*Δt, j=0,...,N (含端点)
// y0: 初始值 y(0)
// f: 右端函数 f(t, y)
// alpha: 分数阶 α ∈ (0, 1] (α=1 时退化为标准 ODE)
// n_corrector: 校正迭代次数 (默认 1, 增加可提高精度)
std::vector<Complex> solve_fractional_ode(
    const std::vector<Real>& t_grid,
    Complex y0,
    const ComplexODEFn& f,
    Real alpha,
    Size n_corrector = 1);

// 分数阶 Riemann-Liouville 积分: I^α g(t) = (1/Γ(α)) ∫₀^t (t-s)^{α-1} g(s) ds
// 在均匀网格上用产品梯形法则 (与校正步核权重一致, 对常数/线性被积函数精确)
Real fractional_integral(
    const std::vector<Real>& t_grid,
    const std::vector<Real>& g,
    Real alpha,
    Real T_idx);  // 积分上限索引

Complex fractional_integral_complex(
    const std::vector<Real>& t_grid,
    const std::vector<Complex>& g,
    Real alpha,
    Real T_idx);

inline std::vector<Complex> solve_fractional_ode(
    const std::vector<Real>& t_grid,
    Complex y0,
    const ComplexODEFn& f,
    Real alpha,
    Size n_corrector) {
    if (t_grid.size() < 2) {
        throw std::invalid_argument(
            "solve_fractional_ode: t_grid must contain at least 2 points");
    }
    if (!(alpha > Real(0)) || !(alpha <= Real(1))) {
        throw std::invalid_argument("solve_fractional_ode: alpha must be in (0, 1]");
    }
    if (n_corrector == 0) n_corrector = 1;

    const Size N = t_grid.size() - 1;
    const Real dt = t_grid[1] - t_grid[0];
    if (!(dt > Real(0))) {
        throw std::invalid_argument("solve_fractional_ode: t_grid must be increasing");
    }
    for (Size k = 2; k < t_grid.size(); ++k) {
        Real diff = std::abs((t_grid[k] - t_grid[k - 1]) - dt);
        if (diff > Real(1e-10) * std::max(Real(1), std::abs(dt))) {
            throw std::invalid_argument("solve_fractional_ode: t_grid must be uniform");
        }
    }

    const Real gamma_a1 = detail::gamma_function(alpha + Real(1));
    const Real gamma_a2 = detail::gamma_function(alpha + Real(2));
    const Real dt_a = std::pow(dt, alpha);
    const Real c_pred = dt_a / gamma_a1;   // Δt^α / Γ(α+1)
    const Real c_corr = dt_a / gamma_a2;   // Δt^α / Γ(α+2)

    std::vector<Complex> y(N + 1);
    std::vector<Complex> fv(N + 1);
    y[0] = y0;
    fv[0] = f(t_grid[0], y0);

    for (Size n = 0; n < N; ++n) {
        // ---- 预测步 (Adams-Bashforth) ----
        Complex sum_p(0, 0);
        for (Size j = 0; j <= n; ++j) {
            sum_p += fv[j] * detail::adams_predictor_weight(alpha, n, j);
        }
        Complex y_pred = y0 + c_pred * sum_p;

        // ---- 校正步 (Adams-Moulton): 先累加已知点固定部分 ----
        Real a0 = std::pow(static_cast<Real>(n), alpha + Real(1))
                - (static_cast<Real>(n) - alpha) * std::pow(static_cast<Real>(n + 1), alpha);
        Complex fixed = a0 * fv[0];
        for (Size j = 1; j <= n; ++j) {
            fixed += detail::adams_corrector_weight(alpha, n, j) * fv[j];
        }
        Complex base = y0 + c_corr * fixed;

        // ---- 校正迭代: 仅新点 f(t_{n+1}, y) 随迭代变化 ----
        Complex y_cur = y_pred;
        for (Size it = 0; it < n_corrector; ++it) {
            y_cur = base + c_corr * f(t_grid[n + 1], y_cur);
        }

        y[n + 1] = y_cur;
        fv[n + 1] = f(t_grid[n + 1], y_cur);
    }

    return y;
}

namespace detail {

// 通用分数阶积分核: 在网格索引 m (>=1) 处计算 I^α g(t_m)
// 使用与 ODE 校正步一致的产品梯形权重:
//   I^α g(t_m) = (Δt^α/Γ(α+2)) · [ g_m + a₀^{(N)} g₀ + Σⱼ₌₁ᴺ aⱼ^{(N)} gⱼ ]
//   N = m-1, a₀^{(N)} = N^{α+1} - (N-α)(N+1)^α
//           aⱼ^{(N)} = (N-j+2)^{α+1} + (N-j)^{α+1} - 2(N-j+1)^{α+1}
template <typename T>
T fractional_integral_impl(const std::vector<Real>& t_grid, const std::vector<T>& g,
                           Real alpha, Real T_idx) {
    if (t_grid.size() < 2) {
        throw std::invalid_argument("fractional_integral: t_grid too small");
    }
    if (!(alpha > Real(0)) || !(alpha <= Real(1))) {
        throw std::invalid_argument("fractional_integral: alpha must be in (0, 1]");
    }
    if (g.size() != t_grid.size()) {
        throw std::invalid_argument("fractional_integral: g size must match t_grid");
    }
    if (T_idx <= Real(0)) return T(0);
    Size m = static_cast<Size>(std::llround(T_idx));
    if (m >= t_grid.size()) m = t_grid.size() - 1;

    const Size N = m - 1;
    const Real dt = t_grid[1] - t_grid[0];
    const Real gamma_a2 = gamma_function(alpha + Real(2));
    const Real c_corr = std::pow(dt, alpha) / gamma_a2;

    T sum = g[m];
    Real a0 = std::pow(static_cast<Real>(N), alpha + Real(1))
            - (static_cast<Real>(N) - alpha) * std::pow(static_cast<Real>(N + 1), alpha);
    sum += static_cast<T>(a0) * g[0];
    for (Size j = 1; j <= N; ++j) {
        sum += static_cast<T>(adams_corrector_weight(alpha, N, j)) * g[j];
    }
    return static_cast<T>(c_corr) * sum;
}

}  // namespace detail

inline Real fractional_integral(const std::vector<Real>& t_grid, const std::vector<Real>& g,
                                Real alpha, Real T_idx) {
    return detail::fractional_integral_impl(t_grid, g, alpha, T_idx);
}

inline Complex fractional_integral_complex(const std::vector<Real>& t_grid,
                                           const std::vector<Complex>& g,
                                           Real alpha, Real T_idx) {
    return detail::fractional_integral_impl(t_grid, g, alpha, T_idx);
}

}  // namespace v1
}  // namespace cpphub
