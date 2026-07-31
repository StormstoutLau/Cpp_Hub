# G. Rough Heston 解析层 (分数阶 Riccati CF + Adams 求解器) — 任务规范

## 执行站: A 站 (scott-lau-NEX.local)
## 模型: opencode/deepseek-v4-flash-free --auto

## 背景
Rough Heston (El Euch-Rosenbaum 2018) 将经典 Heston 方差过程的驱动噪声替换为 fractional Brownian motion,
方差过程服从 fractional CIR: D^α v(t) = κ(θ - v(t)) + σ√v(t) dW(t), α = H + 1/2 ∈ (0.5, 1).
特征函数无闭式解, 需数值求解分数阶 Riccati ODE.

## 文献
- El Euch, Rosenbaum (2018) "The characteristic function of rough Heston model" (Mathematical Finance)
- El Euch, Fukasawa, Rosenbaum (2018) "The microstructural foundations of leverage effect and volatility feedback"
- Li, Ding (2015) "Solving nonlinear fractional differential equations by Adams-Bashforth-Moulton method"
- Jaisson, Rosenbaum (2016) "Rough fractional volatility: Markovian approximation"

## 复用资源
- `pricing/analytic/heston_cf.hpp`: HestonCFParams + heston_characteristic_function (退化验证基准)
- `pricing/fourier/characteristic_functions.hpp`: CharFn 类型别名 + make_heston_cf 工厂模式参考
- `pricing/fourier/cos_method.hpp`: COSEngine (定价验证)
- `core/types.hpp`: Real, Complex, Size

## 需新增文件
1. `include/cpphub/numerics/fractional_adams.hpp` — 分数阶 ODE 求解器 (Caputo 导数 + Adams 预测校正)
2. `include/cpphub/pricing/analytic/rough_heston_cf.hpp` — Rough Heston 特征函数
3. `tests/unit/numerics/test_fractional_adams.cpp` — 求解器单元测试
4. `tests/unit/pricing/test_rough_heston_cf.cpp` — CF 单元测试
5. `tests/CMakeLists.txt` — 注册测试 (用 cpphub_add_test 宏)

## 数学规范

### 1. Caputo 分数阶导数 (α ∈ (0,1))
D^α f(t) = (1/Γ(1-α)) ∫₀^t f'(s) / (t-s)^α ds

### 2. 分数阶 Riccati ODE (El Euch-Rosenbaum 2018, eq. 3.4)
D^α h(u,t) = (1/2)(u² + iu) + (ρσiu - κ)·h(u,t) - (σ²/2)·h(u,t)²
h(u, 0) = 0

其中 u 为复数 (特征函数自变量), h 为复值函数.

### 3. 特征函数
φ(u) = exp(iu·ln S₀ + iu(r-q)T + v₀·h(u,T) + κθ·I^α h(u,T))

其中 I^α h(u,T) = (1/Γ(α)) ∫₀^T (T-s)^{α-1} h(u,s) ds  (Riemann-Liouville 分数阶积分)

### 4. H→0.5 (α→1) 退化
α=1 时 D^1 h = dh/dt, I^1 h = ∫h ds, 退化为标准 Heston Riccati ODE:
dh/dt = (1/2)(u²+iu) + (ρσiu-κ)h - (σ²/2)h²
解析解即经典 Heston CF (用 heston_characteristic_function 验证, 容差 1e-4).

## 接口规范

### fractional_adams.hpp
```cpp
#pragma once
#include "cpphub/core/types.hpp"
#include <vector>
#include <functional>
#include <cmath>

namespace cpphub {
inline namespace v1 {

// ============ Caputo 分数阶导数数值实现 ============
// 求解 D^α y(t) = f(t, y(t)), y(0) = y0, α ∈ (0,1)
// 使用 Adams-Bashforth-Moulton 预测校正 (Li-Ding 2015)
//
// 预测步 (Adams-Bashforth):
//   y_{n+1}^P = Σⱼ₌₀ⁿ bⱼ · f(tⱼ, yⱼ)
// 校正步 (Adams-Moulton):
//   y_{n+1} = Σⱼ₌₀ⁿ bⱼ · f(tⱼ, yⱼ) + cₙ · f(t_{n+1}, y_{n+1}^P)
//
// 系数: bⱼ = (Δt^α / Γ(α+1)) · [binom(α, n-j) - binom(α, n-j-1)]  (j < n)
//        bₙ = (Δt^α / Γ(α+1)) · binom(α, 0)
//        cₙ = (Δt^α / Γ(α+2))
// 其中 binom(α, k) = α(α-1)...(α-k+1)/k!

using ComplexODEFn = std::function<Complex(Real, Complex)>;

// 求解复值分数阶 ODE, 返回 [y(t_0), y(t_1), ..., y(t_N)]
// t_grid: 均匀网格 t_j = j*Δt, j=0,...,N (含端点)
// y0: 初始值 y(0)
// f: 右端函数 f(t, y)
// alpha: 分数阶 α ∈ (0,1)
// n_corrector: 校正迭代次数 (默认 1, 增加可提高精度)
std::vector<Complex> solve_fractional_ode(
    const std::vector<Real>& t_grid,
    Complex y0,
    const ComplexODEFn& f,
    Real alpha,
    Size n_corrector = 1);

// 分数阶 Riemann-Liouville 积分: I^α g(t) = (1/Γ(α)) ∫₀^t (t-s)^{α-1} g(s) ds
// 在均匀网格上用梯形/矩形积分
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

}  // namespace v1
}  // namespace cpphub
```

### rough_heston_cf.hpp
```cpp
#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/pricing/fourier/characteristic_functions.hpp"  // CharFn
#include <vector>

namespace cpphub {
inline namespace v1 {

struct RoughHestonCFParams {
    Real H;        // Hurst 参数 ∈ (0, 0.5), α = H + 0.5
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
void validate_rough_heston_cf_params(const RoughHestonCFParams& p);

// α = H + 0.5
inline Real rough_heston_alpha(Real H) { return H + 0.5; }

// 求解分数阶 Riccati ODE, 返回 h(u, t_grid) (复值序列)
// 用于调试和测试, 不直接对外暴露但便于验证
std::vector<Complex> solve_rough_heston_riccati(
    Complex u,
    const RoughHestonCFParams& p,
    const std::vector<Real>& t_grid);

// Rough Heston 特征函数 φ(u) = E[exp(iu·ln S_T)]
Complex rough_heston_characteristic_function(Complex u, const RoughHestonCFParams& p);

// CharFn 工厂 (COSEngine 兼容)
CharFn make_rough_heston_cf(const RoughHestonCFParams& p);

}  // namespace v1
}  // namespace cpphub
```

## 测试要求

### test_fractional_adams.cpp (至少 6 测试)
1. **常函数**: f(t,y)=0, y0=c → y(t)=c (Caputo 导数零)
2. **线性函数**: f(t,y)=1, y0=0 → y(t)=t^α/Γ(α+1) (Caputo 导数精确解)
3. **指数函数**: f(t,y)=λy, y0=1 → y(t)=E_α(λt^α) (Mittag-Leffler, 数值匹配)
4. **α=1 退化**: α=1 时退化为标准 ODE y'=f(t,y), 与 RK4 对比 (容差 1e-6)
5. **复值 ODE**: f(t,y)=iy, y0=1 → |y(t)|=1 (守恒)
6. **分数阶积分**: I^α [1] = t^α/Γ(α+1), 数值匹配解析

### test_rough_heston_cf.cpp (至少 12 测试)
1. **u=0 归一**: φ(0) = 1
2. **模长 ≤ 1**: |φ(u)| ≤ 1 对多个 u
3. **共轭对称**: φ(-u) = conj(φ(u)) (实值分布)
4. **H→0.5 退化为 Heston**: H=0.499 (α≈1) 时 φ 与 heston_characteristic_function 一致 (容差 1e-3, 因 α≠1 有小偏差)
5. **H=0.5 精确退化**: H=0.5 (α=1) 时退化为标准 Heston, 容差 1e-4
6. **COS Call-Put 平价**: COS 定价满足 call - put = S₀e^{-qT} - Ke^{-rT} (容差 1e-4)
7. **COS vs Heston 一致性**: H=0.5 时 Rough Heston COS 价格 == Heston COS 价格 (容差 1e-4)
8. **ATM smile 倾斜**: H<0.5 时 IV skew 比 Heston 更陡 (rough volatility 特征)
9. **rho 单调性**: rho→-1 时 IV skew 加深
10. **参数校验**: H≥0.5 / H≤0 / sigma≤0 抛异常
11. **数值稳定性**: 大 |u| (|u|=100) 不产生 NaN/Inf
12. **网格收敛**: n_steps 从 100→500, CF 值收敛 (差值 < 1e-4)

## 实现要点

### Adams-Bashforth-Moulton 预测校正
1. **网格**: 均匀 t_j = j*Δt, Δt = T/n_steps
2. **预测** (显式): y_{n+1}^P = Σⱼ bⱼ f(tⱼ, yⱼ)
3. **校正** (隐式): y_{n+1} = Σⱼ bⱼ f(tⱼ, yⱼ) + c·f(t_{n+1}, y_{n+1}^P)
4. **系数预计算**: bⱼ 只依赖 α 和索引差, 可缓存; 用 std::lgamma 计算 Γ 函数
5. **复值支持**: Complex 运算, 注意 f(t,y) 对 h 的二次项 -0.5σ²h²

### Rough Heston CF 计算
1. **预解 Riccati**: 对每个 u, 在 t_grid 上求解 h(u, t), 得 h(u, T) 和 I^α h(u, T)
2. **特征函数**: φ(u) = exp(iu·lnS₀ + iu(r-q)T + v₀·h(u,T) + κθ·I^α h(u,T))
3. **缓存策略**: make_rough_heston_cf 工厂返回闭包, 内部预计算 t_grid 和 Adams 系数, 每次调用 φ(u) 只解一次 Riccati
4. **性能**: 单次 CF 求解约 O(N²) (Adams 求和), N=200 时 ~0.1ms, COS 256 项约 25ms

### H→0.5 退化验证
- α=1 时 D^1 h = dh/dt, 分数阶积分 I^1 h = ∫h ds
- 此时 Riccati ODE 退化为标准 Heston, 有半解析解
- 用 heston_characteristic_function (Albrecher trap) 作为基准, 容差 1e-4

## 验证标准
- 编译: g++ -std=c++17 -O2 (A 站 GCC 13.3.0)
- 测试: 全部通过
- 跨平台: 主站 MSVC 2022 编译通过
- 精度: H=0.5 退化容差 1e-4, COS 平价容差 1e-4
