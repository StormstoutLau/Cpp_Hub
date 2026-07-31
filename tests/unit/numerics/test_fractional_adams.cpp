// 分数阶 ODE 求解器单元测试 (fractional_adams.hpp)
// 覆盖: 常函数 / 线性函数精确解 / Mittag-Leffler 指数 / α=1 退化 RK4 /
//       复值守恒 / 分数阶积分 / 复值 Mittag-Leffler
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <complex>
#include <functional>
#include <algorithm>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/numerics/fractional_adams.hpp"

using namespace cpphub;

namespace {

std::vector<Real> make_uniform_grid(Real T, Size n) {
    std::vector<Real> g(n + 1);
    for (Size k = 0; k <= n; ++k) {
        g[k] = T * static_cast<Real>(k) / static_cast<Real>(n);
    }
    return g;
}

// Mittag-Leffler 函数 E_α(z) = Σ_{k=0}^∞ z^k / Γ(αk + 1) (级数截断)
Complex mittag_leffler(Real alpha, Complex z) {
    Complex sum(0, 0);
    for (Size k = 0; k < 400; ++k) {
        Complex term = std::pow(z, static_cast<Real>(k))
                     / std::tgamma(alpha * static_cast<Real>(k) + 1.0);
        sum += term;
        if (k > 0 && std::abs(term) < 1e-15 * std::max(Real(1), std::abs(sum))) break;
    }
    return sum;
}

// 经典 4 阶 Runge-Kutta (用于 α=1 退化对照)
Real rk4_integrate(const std::function<Real(Real, Real)>& f, Real y0, Real T, Size n) {
    const Real dt = T / static_cast<Real>(n);
    Real t = 0.0, y = y0;
    for (Size i = 0; i < n; ++i) {
        Real k1 = f(t, y);
        Real k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1);
        Real k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2);
        Real k4 = f(t + dt, y + dt * k3);
        y += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t += dt;
    }
    return y;
}

}  // namespace

// 1. 常函数: f(t,y)=0, y0=c → y(t)=c (Caputo 导数零)
TEST(FractionalAdamsTest, ConstantFunction) {
    auto grid = make_uniform_grid(1.0, 64);
    auto y = solve_fractional_ode(grid, Complex(2.5, 0),
                                  [](Real, Complex) { return Complex(0, 0); }, 0.5, 1);
    ASSERT_EQ(y.size(), grid.size());
    for (Size i = 0; i < y.size(); ++i) {
        EXPECT_NEAR(std::real(y[i]), 2.5, 1e-12) << "idx=" << i;
        EXPECT_NEAR(std::imag(y[i]), 0.0, 1e-12) << "idx=" << i;
    }
}

// 2. 线性函数: f(t,y)=1, y0=0 → y(t)=t^α/Γ(α+1) (产品梯形对常数被积函数精确)
TEST(FractionalAdamsTest, LinearFunctionExact) {
    const Real alpha = 0.5;
    const Size N = 200;
    auto grid = make_uniform_grid(1.0, N);
    auto y = solve_fractional_ode(grid, Complex(0, 0),
                                  [](Real, Complex) { return Complex(1, 0); }, alpha, 1);
    for (Size i = 0; i <= N; i += 25) {
        Real exact = std::pow(grid[i], alpha) / std::tgamma(alpha + 1.0);
        EXPECT_NEAR(std::real(y[i]), exact, 1e-8) << "t=" << grid[i];
    }
}

// 3. 指数函数: f(t,y)=λy, y0=1 → y(t)=E_α(λt^α) (Mittag-Leffler)
TEST(FractionalAdamsTest, MittagLefflerExponential) {
    const Real alpha = 0.5;
    const Real lambda = -1.0;
    const Real T = 1.0;
    const Size N = 1000;
    auto grid = make_uniform_grid(T, N);
    auto y = solve_fractional_ode(grid, Complex(1, 0),
                                  [lambda](Real, Complex h) { return lambda * h; }, alpha, 1);
    Complex exact = mittag_leffler(alpha, lambda * std::pow(T, alpha));
    EXPECT_NEAR(std::real(y[N]), std::real(exact), 1e-3);
    EXPECT_NEAR(std::imag(y[N]), std::imag(exact), 1e-3);
}

// 4. α=1 退化: 退化为标准 ODE y'=y-t²+1, 与 RK4 / 解析解对比 (容差 1e-6)
TEST(FractionalAdamsTest, AlphaOneMatchesRK4) {
    const Real T = 1.0;
    const Size N = 2000;
    auto grid = make_uniform_grid(T, N);
    auto f = [](Real t, Complex y) {
        return y - Complex(t * t, 0) + Complex(1, 0);
    };
    auto y = solve_fractional_ode(grid, Complex(0.5, 0), f, 1.0, 1);

    Real exact = (T + 1.0) * (T + 1.0) - 0.5 * std::exp(T);
    EXPECT_NEAR(std::real(y[N]), exact, 1e-6) << "vs analytic";

    Real rk4 = rk4_integrate([](Real t, Real yr) { return yr - t * t + 1.0; },
                             0.5, T, N);
    EXPECT_NEAR(std::real(y[N]), rk4, 1e-6) << "vs RK4";
}

// 5. 复值 ODE: f(t,y)=iy, y0=1 → |y(t)|=1 (α=1 退化为梯形法则, 保模良好)
TEST(FractionalAdamsTest, ComplexConservation) {
    const Real T = 2.0 * PI;
    const Size N = 2000;
    auto grid = make_uniform_grid(T, N);
    auto y = solve_fractional_ode(grid, Complex(1, 0),
                                  [](Real, Complex h) { return Complex(0, 1) * h; }, 1.0, 1);
    EXPECT_NEAR(std::abs(y[N]), 1.0, 1e-6);
    EXPECT_NEAR(std::real(y[N]), std::cos(T), 1e-4);
    EXPECT_NEAR(std::imag(y[N]), std::sin(T), 1e-4);
}

// 6. 分数阶积分: I^α[1] = t^α/Γ(α+1), I^α[t] = t^{α+1}/Γ(α+2)
TEST(FractionalAdamsTest, FractionalIntegral) {
    const Real alpha = 0.5;
    const Size N = 200;
    auto grid = make_uniform_grid(1.0, N);
    std::vector<Real> ones(N + 1, 1.0);
    std::vector<Real> tvals(N + 1);
    for (Size i = 0; i <= N; ++i) tvals[i] = grid[i];

    for (Size i : {N / 2, N}) {
        Real I1 = fractional_integral(grid, ones, alpha, static_cast<Real>(i));
        Real exact1 = std::pow(grid[i], alpha) / std::tgamma(alpha + 1.0);
        EXPECT_NEAR(I1, exact1, 1e-9) << "idx=" << i;

        Real It = fractional_integral(grid, tvals, alpha, static_cast<Real>(i));
        Real exact_t = std::pow(grid[i], alpha + 1.0) / std::tgamma(alpha + 2.0);
        EXPECT_NEAR(It, exact_t, 1e-9) << "idx=" << i;
    }

    // 复值版本与实值一致
    std::vector<Complex> cvals(tvals.begin(), tvals.end());
    Complex Ic = fractional_integral_complex(grid, cvals, alpha, static_cast<Real>(N));
    EXPECT_NEAR(std::real(Ic), std::pow(1.0, alpha + 1.0) / std::tgamma(alpha + 2.0), 1e-9);
    EXPECT_NEAR(std::imag(Ic), 0.0, 1e-9);
}

// 7. 复值分数阶 ODE: f(t,y)=iy, y0=1, α=0.5 → y(T)=E_0.5(i·T^α)
TEST(FractionalAdamsTest, ComplexMittagLeffler) {
    const Real alpha = 0.5;
    const Real T = 1.0;
    const Size N = 2000;
    auto grid = make_uniform_grid(T, N);
    auto y = solve_fractional_ode(grid, Complex(1, 0),
                                  [](Real, Complex h) { return Complex(0, 1) * h; }, alpha, 1);
    Complex exact = mittag_leffler(alpha, Complex(0, 1) * std::pow(T, alpha));
    EXPECT_NEAR(std::real(y[N]), std::real(exact), 1e-3);
    EXPECT_NEAR(std::imag(y[N]), std::imag(exact), 1e-3);
}

// 8. 参数校验: 非法 alpha / 网格 抛异常
TEST(FractionalAdamsTest, InvalidArguments) {
    auto grid = make_uniform_grid(1.0, 32);
    ComplexODEFn f = [](Real, Complex) { return Complex(0, 0); };
    EXPECT_THROW(solve_fractional_ode(grid, Complex(0, 0), f, 0.0, 1), std::invalid_argument);
    EXPECT_THROW(solve_fractional_ode(grid, Complex(0, 0), f, 1.5, 1), std::invalid_argument);
    EXPECT_THROW(solve_fractional_ode({}, Complex(0, 0), f, 0.5, 1), std::invalid_argument);
    std::vector<Real> nonuniform{0.0, 0.3, 1.0};
    EXPECT_THROW(solve_fractional_ode(nonuniform, Complex(0, 0), f, 0.5, 1),
                 std::invalid_argument);
    EXPECT_THROW(fractional_integral(grid, std::vector<Real>(grid.size() - 1, 1.0),
                                     0.5, 10.0), std::invalid_argument);
}
