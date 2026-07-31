#pragma once
// SOURCE: El Euch, Rosenbaum (2018) "The characteristic function of rough Heston
//         model" (Mathematical Finance), eq. 3.4
// SOURCE: El Euch, Fukasawa, Rosenbaum (2018) "The microstructural foundations of
//         leverage effect and volatility feedback"
//
// Rough Heston 特征函数 (El Euch-Rosenbaum 2018):
//   方差过程: D^α v(t) = κ(θ - v(t)) + σ√v(t) dW(t),  α = H + 1/2 ∈ (0.5, 1)
//   特征函数 φ(u) = E[exp(iu·ln S_T)] 由分数阶 Riccati ODE 给出:
//     D^α h(u,t) = -(1/2)(u² + iu) + (ρσiu - κ)·h(u,t) + (σ²/2)·h(u,t)²,  h(u,0) = 0
//     φ(u) = exp(iu·ln S₀ + iu(r-q)T + v₀·h(u,T) + κθ·I^α h(u,T))
//   其中 I^α 为 Riemann-Liouville 分数阶积分 (见 fractional_adams.hpp).
//
// 符号约定: 本库 heston_characteristic_function (Albrecher log-of-ratio 形式) 的
// Riccati 右端为 -½(u²+iu) + (ρσiu-κ)·h + ½σ²h² (二次项为正, 与 El Euch-Rosenbaum
// 论文的 ½(u²+iu) - ... - ½σ²h² 相差整体负号, 对应共轭形式). 为使 H=0.5 (α=1)
// 时精确退化为标准 Heston CF (容差 1e-4), 此处采用与本库 Heston 一致的符号.
//
// 数值方法: h 用 Adams-Bashforth-Moulton 预测校正求解 (fractional_adams.hpp),
// I^α h 用产品梯形法则在同一网格上求值.
//
// 退化验证: α = 1 (H = 0.5) 时 D^1 h = dh/dt, I^1 h = ∫h ds,
// 退化为标准 Heston Riccati ODE, 与 heston_characteristic_function 一致
// (容差 1e-4, 用 heston_cf.hpp 作为基准).
#include "cpphub/core/types.hpp"
#include "cpphub/numerics/fractional_adams.hpp"
#include "cpphub/pricing/fourier/characteristic_functions.hpp"
#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace cpphub {
inline namespace v1 {

struct RoughHestonCFParams {
    Real H;        // Hurst 参数 ∈ (0, 0.5], α = H + 0.5 (H=0.5 允许以支持 α=1 退化验证)
    Real kappa;    // 均值回归速度
    Real theta;    // 长期方差
    Real sigma;    // vol of vol
    Real rho;      // 相关性 ∈ (-1, 1)
    Real v0;       // 初始方差
    Real S0;       // 标的现价
    Real r;        // 利率
    Real q;        // 股息率
    Real T;        // 到期时间
    Size n_steps;  // Adams 求解器网格步数 (默认 200, 精度要求高时 500)
};

// 验证参数有效性
inline void validate_rough_heston_cf_params(const RoughHestonCFParams& p) {
    if (!(p.H > Real(0)) || !(p.H <= Real(0.5))) {
        throw std::invalid_argument("rough_heston_cf: H must be in (0, 0.5]");
    }
    if (p.kappa < Real(0)) {
        throw std::invalid_argument("rough_heston_cf: kappa must be non-negative");
    }
    if (p.theta < Real(0)) {
        throw std::invalid_argument("rough_heston_cf: theta must be non-negative");
    }
    if (!(p.sigma > Real(0))) {
        throw std::invalid_argument("rough_heston_cf: sigma must be positive");
    }
    if (!(p.rho > Real(-1)) || !(p.rho < Real(1))) {
        throw std::invalid_argument("rough_heston_cf: rho must be in (-1, 1)");
    }
    if (p.v0 < Real(0)) {
        throw std::invalid_argument("rough_heston_cf: v0 must be non-negative");
    }
    if (!(p.S0 > Real(0))) {
        throw std::invalid_argument("rough_heston_cf: S0 must be positive");
    }
    if (!(p.T > Real(0))) {
        throw std::invalid_argument("rough_heston_cf: T must be positive");
    }
    if (p.n_steps < 8) {
        throw std::invalid_argument("rough_heston_cf: n_steps too small (>= 8 required)");
    }
}

// α = H + 0.5
inline Real rough_heston_alpha(Real H) { return H + Real(0.5); }

// 求解分数阶 Riccati ODE, 返回 h(u, t_grid) (复值序列)
// 用于调试和测试, 不直接对外暴露但便于验证
inline std::vector<Complex> solve_rough_heston_riccati(
    Complex u,
    const RoughHestonCFParams& p,
    const std::vector<Real>& t_grid) {
    validate_rough_heston_cf_params(p);
    const Complex i(0, 1);
    const Real alpha = rough_heston_alpha(p.H);
    // D^α h = a + b·h + c·h²  (符号与本库 heston_characteristic_function 一致,
    // 使 H=0.5 精确退化为标准 Heston Riccati)
    Complex a = -Real(0.5) * (u * u + i * u);
    Complex b = p.rho * p.sigma * i * u - p.kappa;
    Complex c = Real(0.5) * p.sigma * p.sigma;
    ComplexODEFn f = [a, b, c](Real, Complex h) -> Complex {
        return a + b * h + c * h * h;
    };
    return solve_fractional_ode(t_grid, Complex(0, 0), f, alpha, 1);
}

namespace detail {

// 均匀网格 t_j = j·T/n_steps, j=0..n_steps
inline std::vector<Real> make_rough_heston_grid(Real T, Size n_steps) {
    std::vector<Real> t_grid(n_steps + 1);
    for (Size k = 0; k <= n_steps; ++k) {
        t_grid[k] = T * static_cast<Real>(k) / static_cast<Real>(n_steps);
    }
    return t_grid;
}

}  // namespace detail

// Rough Heston 特征函数 φ(u) = E[exp(iu·ln S_T)]
inline Complex rough_heston_characteristic_function(Complex u, const RoughHestonCFParams& p) {
    validate_rough_heston_cf_params(p);
    std::vector<Real> t_grid = detail::make_rough_heston_grid(p.T, p.n_steps);
    std::vector<Complex> h = solve_rough_heston_riccati(u, p, t_grid);

    const Real alpha = rough_heston_alpha(p.H);
    Complex Ih = fractional_integral_complex(t_grid, h, alpha,
                                             static_cast<Real>(p.n_steps));

    const Complex i(0, 1);
    Complex phase = i * u * (std::log(p.S0) + (p.r - p.q) * p.T);
    Complex exponent = phase + Complex(p.v0, 0) * h.back()
                     + Complex(p.kappa * p.theta, 0) * Ih;
    return std::exp(exponent);
}

// CharFn 工厂 (COSEngine 兼容)
// 闭包预计算 t_grid, 每次调用 φ(u) 只解一次 Riccati
inline CharFn make_rough_heston_cf(const RoughHestonCFParams& p) {
    validate_rough_heston_cf_params(p);
    RoughHestonCFParams params = p;
    std::vector<Real> t_grid = detail::make_rough_heston_grid(params.T, params.n_steps);

    return [params, t_grid](Complex u) -> Complex {
        std::vector<Complex> h = solve_rough_heston_riccati(u, params, t_grid);

        const Real alpha = rough_heston_alpha(params.H);
        Complex Ih = fractional_integral_complex(t_grid, h, alpha,
                                                 static_cast<Real>(params.n_steps));

        const Complex i(0, 1);
        Complex phase = i * u * (std::log(params.S0) + (params.r - params.q) * params.T);
        Complex exponent = phase + Complex(params.v0, 0) * h.back()
                         + Complex(params.kappa * params.theta, 0) * Ih;
        return std::exp(exponent);
    };
}

}  // namespace v1
}  // namespace cpphub
