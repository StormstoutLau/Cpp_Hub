# Cpp_Hub 使用说明

> 版本: 1.0
> 更新: 2026-08-01
> 范围: 当前代码仓库全部模块（含 Phase 1-3 + Phase 4 LITE + 维护迭代内容）
> 定位: 工程级 C++20 量化定价库，header-first，附带 Python 绑定与可选 CUDA MC

本文档面向**库使用者**与**二次开发者**，提供：
1. 仓库目录结构与模块组织
2. 构建与安装流程
3. 各模块公共接口与调用方法
4. 每个算法的**收敛性 / 精度 / 复杂度**理论分析
5. 典型使用示例

权威性声明：本文档与 `docs/architecture/ARCHITECTURE.md`、`docs/phases/phase{1-4}/PHASEx_SPEC.md` 互补，前者描述"为什么这样设计"，本文档描述"如何使用"。

---

## 1. 仓库目录结构

```
Cpp_Hub/
├── include/cpphub/                # Header-only 主库
│   ├── c_api/                     # C ABI（条件编译 CPPHUB_ENABLE_C_API=ON）
│   │   └── cpphub_c_api.h
│   ├── core/                      # 基础层：types/config/constants/error/math/rng/linalg/simd
│   ├── numerics/                  # 数值方法：fractional_adams（分数阶 ODE）
│   ├── instruments/               # 金融工具
│   │   ├── payoff/                # PayOff 抽象 + vanilla/exotic + factory/bridge
│   │   ├── ir/                    # 利率衍生品：OIS/FRA/IRS/BasisSwap/CapFloor/Swaption + 曲线
│   │   └── credit/                # 信用衍生品：CDS/BasketCDS/CDO/Copula/TRS/CreditSpreadOption
│   ├── models/                    # 随机过程模型
│   │   ├── diffusion/             # 扩散模型：GBM/Heston/HestonQE/Bates/CEV/SABR/VG/RoughBergomi/RoughHeston
│   │   ├── ir/                    # 利率模型：ShortRate(Vasicek/CIR/HW/G2++)/HJM/LMM
│   │   └── vol_surface/           # 波动率曲面：SVI/SSVI/VolSurface/DupireLocalVol
│   ├── pricing/                   # 定价引擎
│   │   ├── analytic/              # 解析解：BSM/Heston CF/Bates CF/CEV/VG/SABR-Hagan/RoughBergomi CF/RoughHeston CF/路径依赖解析
│   │   ├── fourier/               # 傅里叶方法：COS(Fang-Oosterlee)/FFT(Carr-Madan)/CF 工厂
│   │   ├── monte_carlo/           # MC 引擎：MCEngine/QMC/Sobol/BrownianBridge/LSMC/路径生成器/方差缩减
│   │   ├── pde/                   # PDE 方法：1D(CN+PSOR)/2D(ADI Craig-Sneyd/Hundsdorfer-Verwer)/Thomas/FDMGrid
│   │   ├── tree/                  # 树形方法：Binomial(CRR/JR/Tian/Leisen-Reimer)/Trinomial(Explicit/Implicit/Hybrid)
│   │   └── engine.hpp             # 统一定价入口
│   ├── monte_carlo/               # 通用方差缩减：ConditionalMC/ControlVariate/ImportanceSampling/MomentMatching/StratifiedSampling
│   ├── calibration/               # 校准：Calibrator(Heston/SABR/Bates/VG/CEV)/Objective/Optimizer(LM/DE/NelderMead)
│   ├── risk/                      # 风险管理
│   │   ├── greeks/                # Greeks：Analytic/AAD/Pathwise/LR/Numerical/Factory + autodiff 集成
│   │   ├── var/                   # VaR：Historical/Parametric(Normal/StudentT/CF)/MC/ES/Backtesting(Kupiec/Christoffersen/Basel)
│   │   ├── scenario/              # 情景分析：Sensitivity/StressTest(Basel FRTB + 历史危机)
│   │   ├── xva.hpp                # XVA：CVA/DVA/FVA/BVA
│   │   ├── wrong_way_risk.hpp     # WWR：Hull-White/Copula(Pykhtin-Zhu)/MC
│   │   └── pfe_sa_ccr.hpp         # PFE + SA-CCR(Basel III)
│   └── performance/gpu/           # GPU MC（CUDA 可选，CPU stub fallback）
├── src/                           # 非 header-only 源文件（CUDA kernel + CPU stub）
├── tests/                         # GoogleTest 单元 + 集成测试
├── benchmarks/                    # 性能基准
├── python/                        # Python 绑定（nanobind）+ pytest
├── third_party/autodiff/          # vendored autodiff 头文件（MIT）
├── docs/                          # 架构 / Phase 规格 / ADR / 审计 / 任务
└── CMakeLists.txt
```

### 模块依赖层次（自底向上）

```
Layer 0: core/ (types, config, constants, error, math, rng, linalg, simd)
   ↑
Layer 1: numerics/ (fractional_adams)  ← 依赖 core
   ↑
Layer 2: instruments/payoff/ + models/diffusion/  ← 依赖 core
   ↑
Layer 3: pricing/ (analytic, fourier, monte_carlo, pde, tree)  ← 依赖 Layer 0-2
   ↑
Layer 4: calibration/ + models/vol_surface/ + models/ir/ + instruments/ir/credit/  ← 依赖 Layer 0-3
   ↑
Layer 5: risk/ (greeks, var, scenario, xva, wwr, pfe_sa_ccr)  ← 依赖 Layer 0-4
   ↑
Layer 6: monte_carlo/ (通用方差缩减) + performance/gpu/  ← 依赖 Layer 0-3
```

---

## 2. 构建与安装

### 2.1 依赖要求

- C++20 编译器：MSVC 19.x / GCC 13+ / Clang 14+
- CMake ≥ 3.25
- Python ≥ 3.9（可选，用于 Python 绑定）
- CUDA ≥ 12（可选，用于 GPU MC）
- 第三方：GoogleTest（测试）、nanobind（Python 绑定）、autodiff（vendored）

### 2.2 标准构建

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build --build-config Release --output-on-failure
```

### 2.3 可选功能

```bash
# Python 绑定
pip install -e python/
pytest python/tests/

# GPU MC（主控站 RTX 4060 可选）
cmake -S . -B build_cuda -DCMAKE_BUILD_TYPE=Release -DCPPHUB_ENABLE_CUDA=ON
cmake --build build_cuda --config Release

# C ABI（版本化符号 cpphub_v1_*）
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCPPHUB_ENABLE_C_API=ON
```

### 2.4 跨平台浮点确定性

由编译选项保证 IEEE-754 严格模式：
- GCC: `-ffp-contract=off`
- MSVC: `/fp:precise`

跨平台验收基准：MSVC (Win10) + GCC (Ubuntu NEX/GTR-Pro) 三平台 **286/286 位精确一致**。

---

## 3. 各模块详细说明

### 3.1 核心层 `core/`

#### 3.1.1 类型与配置

**文件**: `core/types.hpp`, `core/config.hpp`, `core/constants.hpp`, `core/error.hpp`

| 类型/宏 | 定义 | 说明 |
|---|---|---|
| `Real` | `double` | IEEE 754 双精度 |
| `Complex` | `std::complex<double>` | 复数 |
| `Size` | `std::size_t` | 无符号尺寸 |
| `Index` | `std::ptrdiff_t` | 带符号索引 |
| `PI`, `SQRT_2`, `INV_SQRT_2`, `SQRT_2PI`, `INV_SQRT_2PI`, `LN_2` | `constexpr Real` | 数学常数（取自 `std::numbers`） |
| `CPPHUB_DBL_EPSILON` | `2.220446049250313e-16` | 机器精度 |
| `CPPHUB_VERSION_MAJOR/MINOR/PATCH` | `1.0.0` | 版本号 |
| `CPPHUB_HAS_AVX2/AVX512/NEON` | 编译期检测 | SIMD 能力 |
| `ErrorCode` | enum class | Success/InvalidArgument/OutOfRange/ConvergenceFailure/NotConverged/NumericalError/NotImplemented |
| `CppHubException` | `: std::runtime_error` | 统一异常 |

调用示例：

```cpp
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
using namespace cpphub::v1;

Real x = PI * 2.0;              // 6.283185307179586
Real eps = CPPHUB_DBL_EPSILON;  // 2.220446049250313e-16
```

#### 3.1.2 数学函数 `core/math.hpp`

| 函数 | 算法理论依据 | 精度 | 复杂度 |
|---|---|---|---|
| `erf(x)` / `erfc(x)` | 包装 `std::erf` / `std::erfc` | 机器精度 | O(1) |
| `normal_pdf(x)` | `(1/√(2π)) exp(-x²/2)` | 机器精度 | O(1) |
| `normal_cdf(x)` | `0.5 * erfc(-x/√2)` | 机器精度 | O(1) |
| `inv_normal_cdf(p)` | **Hastings 有理近似 + 3 次 Halley 类迭代** | 接近机器精度 | O(1) |
| `bessel_i0(x)` / `bessel_i1(x)` | 幂级数展开（\|x\|<3.75 时 kmax=15，否则 kmax=50） | 迭代至收敛 | O(50) |
| `regularized_lower_gamma(a,x)` | **Numerical Recipes §6.2**：x<a+1 级数展开，x≥a+1 Lentz 连分式 | ~1e-14 | O(N) 收敛 |
| `regularized_upper_gamma(a,x)` | `1 - regularized_lower_gamma(a,x)` | ~1e-14 | O(N) |
| `noncentral_chi2_cdf(x,k,λ)` | Poisson 加权 Gamma CDF 级数；小 λ≤50 直接级数，大 λ>50 对数空间递推 | 小 λ ~1e-12，大 λ ~1e-6 | O(λ/2 ± 5σ) |

**注意**：`inv_normal_cdf` 用的是 **Hastings 有理近似（Abramowitz-Stegun 26.2.23 风格）+ Halley 迭代**，**非 Acklam 算法**。Hastings 初值精度约 1e-3~1e-4，经 3 次 Halley 迭代后逼近机器精度。

调用示例：

```cpp
#include "cpphub/core/math.hpp"
using namespace cpphub::v1;

Real p = normal_cdf(1.96);        // ≈ 0.975
Real z = inv_normal_cdf(0.975);   // ≈ 1.96
Real nc = noncentral_chi2_cdf(3.0, 2.0, 4.0);  // 非中心卡方 CDF
```

#### 3.1.3 随机数生成 `core/rng.hpp`

**算法**: **Philox 4x64-10**（Salmon et al. 2011, *Parallel Random Numbers: As Easy as 1, 2, 3*, SC）

| 属性 | 值 |
|---|---|
| 算法 | Philox 4x64-10 计数式 ARNG（10 轮 Feistel-like 混合） |
| 周期 | **2^128**（128 位计数器空间） |
| 状态 | `counter_[2]` + `key_[2]`（共 256 位） |
| 输出 | 64 位无符号整数，均匀→[0,1) 用 `u>>11 · 1/2^53`（53 位精度） |
| 正态变换 | Box-Muller：`z1=√(-2 ln u1)·cos(2π u2)`, `z2=√(-2 ln u1)·sin(2π u2)` |
| 相关化 | `generate_correlated(rng, L)`：先独立正态 z，再 `L·z`（L 为 Cholesky 因子） |

**并行性**: counter-based 设计天然支持并行 — 不同 (counter, key) 独立无重叠。GPU 实现与 CPU 位精确一致（相同 seed + counter 编号 → 相同 Z）。

调用示例：

```cpp
#include "cpphub/core/rng.hpp"
using namespace cpphub::v1;

Philox4x64 rng(42, 0);           // seed=42, stream=0
uint64_t u = rng();              // 均匀 [0, 2^64)
auto [z1, z2] = rng.box_muller(u1, u2);  // 标准正态
auto z4 = rng.normal_simd(rng);  // SIMD 4 元素正态
```

#### 3.1.4 线性代数 `core/linalg.hpp`

| 函数 | 算法 | 复杂度 | 备注 |
|---|---|---|---|
| `cholesky<N>(A)` | Cholesky-Banachiewicz 按列分解 `A=L·Lᵀ` | **O(N³)** | 要求正定，否则抛 `NumericalError` |
| `thomas_algorithm<N>(a,b,c,d)` | 追赶法（前消元 + 回代） | **O(N)** | 三对角系统，无主元 |
| `Matrix::determinant()` | ≥4×4 用**带部分主元高斯消元**（LU 形式） | O(N³) | 主元阈值 1e-30 |

定长模板 `Vector<N>` / `Matrix<R,C>` 适合小规模（N≤10）线性代数。

#### 3.1.5 SIMD `core/simd.hpp`

`f64x4`（AVX2 256 位 / 4×double）与 `f64x8`（2×256 位拼接）。AVX2 分支用 `_mm256_*_pd` 内建；非 AVX2 标量回退（逐元素循环）。**注意**：`exp`/`log` 即便在 AVX2 分支也是标量实现（store→逐元素→load 回）。

### 3.2 数值方法层 `numerics/`

#### `numerics/fractional_adams.hpp` — 分数阶 ODE 求解器

**用途**: 求解 Caputo 分数阶初值问题 `D^α y(t) = f(t,y)`, α∈(0,1]，用于 Rough Heston CF 的分数阶 Riccati ODE。

**算法**: **Diethelm-Ford-Freed (2002) 预测校正** + **Li, Ding (2015) Adams-Bashforth-Moulton 产品积分法**

| 属性 | 值 |
|---|---|
| 方法阶数 | **min(2, 1+α)**（对 t^α 类奇异初值） |
| α=1 退化 | 精确退化为标准 Adams-Bashforth-Moulton（梯形法则，2 阶） |
| 复杂度 | **O(N²)**（每步核求和 O(N)） |
| n_corrector | 默认 1，可增大以提高精度 |
| 网格要求 | 必须均匀（容差 1e-10·max(1,\|Δt\|)） |

调用示例：

```cpp
#include "cpphub/numerics/fractional_adams.hpp"
using namespace cpphub::v1;

std::vector<Real> t_grid(101);
for (Size i = 0; i <= 100; ++i) t_grid[i] = 0.01 * i;
ComplexODEFn f = [](Real t, Complex y) { return -0.5 * y; };
auto y = solve_fractional_ode(t_grid, Complex(1.0, 0.0), f, /*alpha=*/0.75);
```

### 3.3 金融工具层 `instruments/`

#### 3.3.1 PayOff 抽象 `instruments/payoff/`

- `PayOff`（抽象基类）：`operator()(Spot)`, `clone()`
- `PayOffCall` / `PayOffPut` / `PayOffDigitalCall` / `PayOffDigitalPut` / `PayOffDoubleDigital`
- `PathDependentAsian` / `PathDependentBarrier` / `PathDependentLookback`（路径依赖）
- `PayOffFactory` + `PayOffBridge`（工厂 + 桥接模式）

调用示例：

```cpp
#include "cpphub/instruments/payoff/vanilla.hpp"
using namespace cpphub::v1;

PayOffCall call(100.0);
Real payoff = call(105.0);  // 5.0
```

#### 3.3.2 利率衍生品 `instruments/ir/`

| 类 | 定价方法 | 算法理论 | 复杂度 |
|---|---|---|---|
| `ZeroCurve` | 分段线性/对数线性插值 + OIS 拔靴 | 求解各期限零息利率使贴现因子匹配市场价 | 插值 O(log n)；拔靴 O(N²) |
| `OIS` | 单曲线 NPV | 固定端 vs 浮动端贴现求和 | O(N) |
| `FRA` | 多曲线 OIS 贴现 | 远期曲线算远期 + OIS 曲线贴现 | O(1) |
| `IRS` | 单/多曲线 NPV + par swap rate | 固定端 = 固定利率 × 年份分数 × 贴现因子之和 | O(N) |
| `BasisSwap` | 多曲线 NPV + par spread | 两腿浮动，不同基准曲线 | O(N) |
| `CapFloor` | Black 76 / 短利率模型 | 每个 caplet 当远期利率欧式期权，Black 76 闭式 | O(N) |
| `Swaption` | Hull-White + **Jamshidian 分解** | 求 r* 使 IRS(r*)=K，分解为零息债期权组合 | O(N) + 牛顿法 3-5 次 |

#### 3.3.3 信用衍生品 `instruments/credit/`

| 类 | 方法 | 算法理论 | 复杂度 |
|---|---|---|---|
| `CreditCurve` | hazard rate 模型 | P(τ>t)=exp(-∫h(s)ds) | O(N) |
| `CDS` | NPV + par spread + 拔靴 | 保费端 + 赔付端贴现 | 拔靴 O(N²) |
| `Copula` | Gaussian / One-Factor / t-Copula | Cholesky 分解 + CDF 逆映射 | Cholesky O(N³) 一次性 + O(N²) 每路径；One-Factor O(N) 每路径 |
| `BasketCDS` | MC / 半解析(Li 2000) | One-Factor 条件独立性，对 M 积分 | MC O(M·N²)；半解析 O(K·N·log N) |
| `CDO` | LHP 近似 + MC | Vasicek 单因子闭式（无限同质资产假设） | LHP 闭式 O(K)；MC O(M·N²) |
| `CreditSpreadOption` | Black 76 + 生存调整 | 利差对数正态，乘以 P(τ>T_option) | O(1) |
| `TRS` | 两腿 NPV | 资产收益端 vs 融资端，信用风险折现 | O(N) |

### 3.4 模型层 `models/`

#### 3.4.1 扩散模型 `models/diffusion/`

所有模型继承 `StochasticProcess` 抽象基类：`generate_path(T, n_steps, path, rng)`, `characteristic_function(u, tau)`, `spot()`。

| 模型 | 离散化方案 | 强收敛阶 | 弱收敛阶 | 复杂度 | 算法理论依据 |
|---|---|---|---|---|---|
| `GBM` | **精确解**（闭式） | ∞（无离散误差） | ∞ | O(N) | `S(t+Δt)=S(t)·exp((μ-σ²/2)Δt+σ√Δt·Z)` |
| `Heston` (基础) | Euler / **Full Truncation (Lord 2010)** | 0.5 | 1.0 | O(N) | 漂移用 v_n，扩散用 max(v_n,0)，保证非负无偏 |
| `HestonQE` | **Quadratic Exponential (Andersen 2008)** + 鞅修正 | **1.0** | 1.0 | O(N) | ψ≤1.5 用非中心卡方二次近似；ψ>1.5 用指数混合；Feller 违反时仍高精度 |
| `Bates` | Euler（方差） + 精确跳跃 | 0.5（方差） | 1.0 | O(N+λT) | Heston + Merton 对数正态跳跃，CF 乘积律 |
| `CEV` | Euler / Log-Euler | 0.5 | 1.0 | O(N) | `dS=μSdt+σS^β dW`，β=1 退化 GBM |
| `SABR` | Euler / Log-Euler | 0.5 | 1.0 | O(N) | `dF=αF^β dW1`, `dα=να dW2`，Cholesky 相关化 |
| `VarianceGamma` | **精确模拟**（Gamma subordinator） | ∞ | ∞ | O(N) | `X(t)=θG(t)+σW(G(t))`，G~Gamma(t/ν,ν) |
| `RoughBergomi` | RL-fBm **Cholesky 精确** + log-Euler 价格 | 0.5（价格） | 1.0 | Cholesky O(N³) 一次性 + O(N²) 每路径 | Bayer-Friz-Gatheral 2016；W̃^H=√(2H+1)∫(t-s)^{H-1/2}dW |
| `RoughBergomi` (Hybrid) | **Hybrid Scheme (BLP 2017)** | 0.5 | 1.0 | **O(N·b)** 每路径 | 近端精确核 + 远端平顶近似，b=1 时近似 O(N) |
| `RoughHeston` | **分数阶 Euler (Smith 2017)** + Volterra 核 | 弱 1.0 | 1.0 | 核预计算 O(N²) + O(N²) 每路径 | El Euch-Rosenbaum 2018；α=H+1/2∈(0.5,1) |

调用示例：

```cpp
#include "cpphub/models/diffusion/heston_qe.hpp"
using namespace cpphub::v1;

HestonParams p{0.04, 2.0, 0.04, 0.3, -0.7, 100.0, 0.05, 0.0};
HestonQE heston(p);
std::vector<Real> path(252);
Philox4x64 rng(42);
heston.generate_path(1.0, 251, path, rng);
```

#### 3.4.2 利率模型 `models/ir/`

| 模型 | 闭式解 | 模拟方案 | 算法理论 |
|---|---|---|---|
| `Vasicek` / `CIR` | 零息债 P=A·exp(-B·r0)；债期权解析 | Euler | 仿射短利率模型 |
| `HullWhite` | 债期权闭式（含 θ(t) 拟合）；Swaption 用 **Jamshidian 分解** | Euler | Vasicek 时变 θ(t) 版本，精确拟合初始期限结构 |
| `G2PlusPlus` | 零息债 + 债期权闭式（Brigo-Mercurio） | Euler | 两因子 Hull-White，r=x+y+φ |
| `HJM` | 无套利漂移 α(t,T)=Σσ_i(t,T)·∫σ_i(t,u)du | Euler（远期曲线演化） | 直接建模 f(t,T)，非马尔可夫 |
| `LMM` | caplet Black 76 闭式；swaption **Rebonato 近似** | log-Euler | 离散期限 LIBOR，spot/terminal measure |

### 3.5 定价引擎层 `pricing/`

#### 3.5.1 解析定价 `pricing/analytic/`

| 文件 | 模型 | 算法理论 | 精度 | 复杂度 |
|---|---|---|---|---|
| `heston_cf.hpp` | Heston CF | **Albrecher 2007 "Little Trap"**（log-of-ratio 形式，避免分支切割不连续） | O(1) | O(1) |
| `bates_cf.hpp` | Bates CF | Heston CF × Merton 跳跃 CF，r̃=r-λm drift 调整 | O(1) | O(1) |
| `cev_analytic.hpp` | CEV 欧式 | **Schroder 1989 非中心卡方分布**；β=1 回退 BS | 依赖 `noncentral_chi2_cdf` | O(1) |
| `path_dependent_analytic.hpp` | 几何亚式 / 障碍 / 回望 | **Kemna-Vorst 1990** / **Reiner-Rubinstein 1991** / **Goldman-Sosin-Gatto 1979 + Conze-Viswanathan 1991** | 闭式 | O(1) |
| `rough_bergomi_cf.hpp` | rBergomi 近似 CF | **累积量展开（Edgeworth/Escher 型）**，只保留 c1/c2/c3 舍 c4 保证 \|φ\|≤1 | Edgeworth 前三阶 | O(1) |
| `rough_heston_cf.hpp` | Rough Heston CF | **El Euch-Rosenbaum 2018 分数阶 Riccati ODE**，Adams-Bashforth-Moulton 求解 | O(n_steps²) | O(n_steps²) |
| `sabr_hagan.hpp` | SABR 隐含波动率 | **Hagan 2002 渐近 O(σ²)**，**Oblój 2008 修正 x(z) 符号** | 渐近 O(σ²) | O(1) |
| `vg_analytic.hpp` | VG CF + 累积量 | Madan-Carr-Chang 1998 闭式 CF | 闭式 | O(1) |

#### 3.5.2 傅里叶方法 `pricing/fourier/`

##### COS 方法 `cos_method.hpp`

**算法**: **Fang-Oosterlee (2009) COS (Fourier-cosine series) 方法**

核心公式：
```
v(x₀,t₀) ≈ e^{-rT}·Σ'_{k=0}^{N-1} Re[φ(kπ/(b-a))·e^{ikπ(x₀-a)/(b-a)}]·V_k
```
- φ 为 ln S_T 的特征函数
- V_k 为 payoff 系数（Call: χ_k - K·ψ_k；Put: K·ψ_k - χ_k）
- Σ' 表示首项权重减半

| 属性 | 值 |
|---|---|
| **收敛阶** | **指数收敛 O(e^{-N})**（对光滑密度）；非光滑 O(1/N) |
| 截断区间 | 自动通过 CF 一阶/二阶矩估计，默认 L=10 标准差 |
| 默认 N | 256 项 |
| 复杂度 | 单 K 定价 O(N)；批量 K 个行权价 O(K·N) |
| 典型 N | GBM 128-512，Heston 512 |

调用示例：

```cpp
#include "cpphub/pricing/fourier/cos_method.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
using namespace cpphub::v1;

HestonCFParams p{0.04, 2.0, 0.04, 0.3, -0.7, 0.05, 0.0};
auto phi = make_heston_cf(100.0, 0.05, 0.0, p, 1.0);
COSEngine::Config cfg;
cfg.n_terms = 512;
COSEngine cos(phi, 100.0, 0.05, 0.0, 1.0, cfg);
Real call = cos.price_call(100.0);
```

##### FFT 方法 `fft_method.hpp`

**算法**: **Carr-Madan (1999) FFT 方法**

| 属性 | 值 |
|---|---|
| Damping factor α | **默认 1.5**（范围 1.5-5，短期期权需较大值） |
| FFT 点数 N | 默认 4096（必须 2 的幂） |
| η（v 步长） | 默认 0.25 |
| λ（strike 网格分辨率） | `λ = 2π/(N·η)` ≈ 0.00614 |
| 积分加权 | **Simpson 法则** |
| FFT 内核 | 自实现 **Radix-2 Cooley-Tukey**（位反转 + 蝶形） |
| **复杂度** | **O(N log N)**（单次调用计算整个 strike 网格 N 个价格） |

#### 3.5.3 蒙特卡洛 `pricing/monte_carlo/`

##### MCEngine `mc_engine.hpp`

通用 MC 引擎，支持 Antithetic + Control Variate 方差缩减。

| 属性 | 值 |
|---|---|
| 方差缩减 | Antithetic（`X_av=(X(Z)+X(-Z))/2`）+ Control Variate（`β*=Cov(X,Y)/Var(Y)`） |
| **收敛阶** | **O(1/√N)** |
| 95% CI | price ± 1.96·std_error |
| 组合 | Antithetic 后应用 CV（Glasserman 4.3） |

##### QMC 引擎 `qmc_engine.hpp`

| 属性 | 值 |
|---|---|
| 算法 | Sobol 低偏差序列 + Brownian Bridge + **Randomized QMC (RQMC)** |
| **收敛阶** | **O((log N)^d / N)**（Koksma-Hlawka 不等式），vs 标准 MC O(1/√N) |
| RQMC | R 个独立 scrambled Sobol 序列（Owen 随机化），`SE=std(θ̂_r)/√R` |
| 默认 n_paths | 4096（建议 2^k） |
| 默认 n_replicates | 16 |

##### Sobol 序列 `sobol.hpp`

| 属性 | 值 |
|---|---|
| 方向数集 | **Joe & Kuo (2008) Table S1** |
| 维度上限 | dim 1 特殊（`v_i=1<<(63-i)`）+ dim 2~20 表内；**dim 21+ 用 dim 20 fallback**（质量略降） |
| 精度 | 64-bit |
| Scrambling | Owen hash-based digital scrambling（简化实现） |

##### Brownian Bridge `brownian_bridge.hpp`

| 类 | 算法 | 复杂度 |
|---|---|---|
| `BrownianBridge` | 中点递归构造：W(T)=√T·Z_0，递归填充中点 | O(N) |
| `MultiAssetBrownianBridge` | Cholesky 相关化多资产 | O(N·d²) |
| `BrownianBridgePCA` | **理论最优 PCA 构造**（Jacobi 特征分解） | O((nd)³) |

**方差缩减配置**: BB 将最重要维度（Z_0）控制 W(T)，有效维度从 n 降到 ~log2(n)，与 Sobol 配合显著提升收敛速度。

##### LSMC 引擎 `lsmc_engine.hpp`

Longstaff-Schwartz 美式/百慕大期权定价。

| 属性 | 值 |
|---|---|
| Basis | Laguerre（默认，物理学家版本）/ Hermite（概率学家 He_n）/ Monomial / Chebyshev（第一类 T_n） |
| 默认阶数 | 3 |
| 回归方法 | **OLS 正规方程 (XᵀX+λI)β=XᵀY**，可选 Ridge 正则化 |
| 解线性方程组 | Gauss-Jordan 消元 |
| 默认 n_paths | 10000 |
| 默认 n_steps | 50 |
| 归一化 | x = S/K（Longstaff-Schwartz 2001） |
| Antithetic | 默认开启 |

##### 路径生成器 `path_generator.hpp`, `multi_asset_path_generator.hpp`

| 方案 | 公式 | 强收敛阶 |
|---|---|---|
| Exact（GBM 闭式解） | `S_{t+dt}=S_t·exp((r-q-0.5σ²)dt+σ√dt·Z)` | ∞ |
| Euler-Maruyama | `S_{t+dt}=S_t·(1+(r-q)dt+σ√dt·Z)` | 0.5 |
| Milstein | Euler + `0.5σ²·Z(Z-1)·dt` | 1.0 |

多资产 GBM 用 Cholesky 分解相关矩阵 R=LL^T，`dW=L·dZ`。

##### 方差缩减 `variance_reduction.hpp`, `antithetic.hpp`

`VREngine` 统一 Decorator，组合 Antithetic + Control Variate。

| 方法 | 原理 | 方差缩减比 |
|---|---|---|
| Antithetic | `X_av=(X(Z)+X(-Z))/2`，payoff 单调时 Cov(X(Z),X(-Z))<0 | ~2x |
| Control Variate | `X_cv=X-β(Y-μ_Y)`, `β*=Cov(X,Y)/Var(Y)` | `1/(1-ρ²)`，ρ→1 时趋无穷 |

#### 3.5.4 PDE 方法 `pricing/pde/`

##### 1D PDE `pde_engine.hpp` + `fdm_scheme.hpp` + `fdm_grid.hpp` + `thomas_solver.hpp`

| 格式 | 算法 | 收敛阶 | 稳定性 |
|---|---|---|---|
| ExplicitEuler | 显式 Euler | 一阶 | 有稳定性限制 (Δt/Δx² ≤ 0.5) |
| ImplicitEuler | 隐式 Euler | 一阶 | **L 稳定**（damping 高频） |
| **CrankNicolson** | theta scheme（θ=0.5） | **二阶时间精度** | 无条件稳定 |
| **RannacherSmoothing** | 前 n_warmup 步 ImplicitEuler + 后续 CrankNicolson | 二阶 | 处理非光滑 payoff（行权价 kink、数字期权）的 CN 振荡 |

**美式定价**: Crank-Nicolson + **PSOR**（Projected Successive Over-Relaxation）求解 LCP。
- PSOR 参数：`omega=1.5`, `tol=1e-8`, `max_iter=2000`
- **收敛条件**: ω∈(1,2)
- 迭代：`V[i]=max(payoff[i], old_val+ω·(y-old_val))`

**网格**: sinh 变换使网格在 S0 附近密集、两端稀疏。

**Thomas 算法**: 三对角系统求解，**O(n)**。

##### 2D PDE (Heston) `pde_engine_2d.hpp` + `fdm_2d.hpp`

**算法**: **ADI（Alternating Direction Implicit）** for Heston PDE

Heston PDE（log-S 变换）：
```
dV/dt = 0.5v·d²V/dx² + (r-q-0.5v)·dV/dx + 0.5ξ²v·d²V/dv² + κ(θ-v)·dV/dv + ρξv·d²V/(dxdv) - rV
```

算子分裂：L = L_x + L_v + L_xv + L_0
- L_x / L_v: 隐式（Thomas 算法求解三对角）
- L_xv（混合导数）: 显式 5-point stencil `(V[i+1,j+1]-V[i-1,j+1]-V[i+1,j-1]+V[i-1,j-1])/(4·dx·dv)`

| 格式 | 算法 | 收敛阶 |
|---|---|---|
| **CraigSneyd** | theta=0.5，两个校正阶段显式加入混合导数 | **二阶时间精度** |
| **HundsdorferVerwer** | theta=0.5（最优 0.5+1/√12），混合导数折叠进完整算子 | **二阶时间精度** |

**边界条件**:
- x_min（S→0）: call V=0，put V=K·exp(-rτ)
- x_max（S→∞）: call V=S_max·exp(-qτ)-K·exp(-rτ)，put V=0
- v_min（v=0）: 退化 PDE（BSM with σ=0）
- v_max: 线性外推（Neumann-like）

#### 3.5.5 树形方法 `pricing/tree/`

##### 二叉树 `binomial.hpp`

| 格式 | 算法 | 收敛阶 |
|---|---|---|
| CRR (Cox-Ross-Rubinstein 1979) | `u=e^{σ√dt}`, `d=1/u`, `p=(e^{(r-q)dt}-d)/(u-d)` | **O(1/n)**（锯齿振荡） |
| Jarrow-Rudd | `u=e^{(r-q-σ²/2)dt+σ√dt}`, `d=e^{(r-q-σ²/2)dt-σ√dt}`, `p=0.5` | O(1/n) |
| Tian (1993) | 三阶矩匹配 | O(1/n) |
| **Leisen-Reimer (1986)** | **Peizer-Pratt 二项分布反演**，强制 n 为奇数 | **O(1/n²)**（二阶收敛，消除偶数步振荡） |

##### 三叉树 `trinomial.hpp`

| 格式 | 算法 | 收敛阶 |
|---|---|---|
| Explicit | Kamrad-Ritchken (1991) 风格，`dx=σ√dt` | O(1/n) |
| Implicit | `dx=σ√(3·dt)`（更宽支撑） | O(1/n) |
| Hybrid | 当前实现等价于 Explicit | O(1/n) |

**复杂度**: 空间 O(n)，时间 O(n²)。

### 3.6 通用方差缩减层 `monte_carlo/`

| 文件 | 方法 | 原理 | 方差缩减比 |
|---|---|---|---|
| `conditional_mc.hpp` | Conditional MC | `Var(X)=Var(E[X\|Y])+E[Var(X\|Y)]`，用 E[X\|Y] 消除后一项 | 恒 ≤ 1（不增） |
| `control_variate.hpp` | Control Variate | `X_cv=X-β(Y-μ_Y)`, `β*=Cov(X,Y)/Var(Y)` | `1/(1-ρ²)` |
| `importance_sampling.hpp` | Importance Sampling | **Girsanov 测度变换**，`L=dQ/dP=exp(θ·Z+0.5·θ²)`，最优 θ*=-d2 | 深度 OTM 数十倍 |
| `moment_matching.hpp` | Moment Matching | 调整样本使样本矩匹配理论矩 | sample_var/theor_var |
| `stratified_sampling.hpp` | Stratified Sampling | [0,1] 分 K 层，**Neyman 最优分配** `n_k∝p_k·σ_k` | 恒 ≤ 1（不增） |

### 3.7 校准层 `calibration/`

#### 优化器 `optimizer.hpp`

| 优化器 | 算法 | 收敛性 | 复杂度 |
|---|---|---|---|
| **LevenbergMarquardt** | `(JᵀJ+λI)dx=-Jᵀr`，J 由中心差分（O(h²)），λ 自适应 | **近极小点二次收敛**（Gauss-Newton 极限），总体超线性 | 每次 iter O(n²) Jac + O(n³) 解方程 |
| **NelderMead** | 标准单纯形（reflect α=1, expand γ=2, contract ρ=0.5, shrink σ=0.5） | 线性（慢），对非光滑/噪声鲁棒 | O(iter×n) |
| **DifferentialEvolution** | **rand/1/bin**：变异 `v=x_r1+F(x_r2-x_r3)`，F=0.8，CR=0.9 | **全局概率 1 收敛**，无梯度 | O(pop×gen×f_eval) |

#### 目标函数 `objective.hpp`

加权方案：
- PriceWeighted: `w=1/√|price|`
- VegaWeighted: `w=1/max(|vega|,ε)`（downweight ATM）
- RelativeError: `w=1/|market|`（默认）
- Mixed: λ·price + (1-λ)·iv

#### 标定器 `calibrator.hpp`

| 标定器 | 参数 | 定价方法 | 约束 |
|---|---|---|---|
| `HestonCalibrator` | 5 参数 [v0,κ,θ,σ,ρ] | Carr-Madan + Gil-Pelaez（4096 梯形）→ IV | Feller: `2κθ>σ²` |
| `SABRCalibrator` | 4 参数 [α,β,ν,ρ] | Hagan 2002 显式 IV | β∈[0,1] |
| `BatesCalibrator` | 8 参数 [v0,κ,θ,σ,ρ,λ,μ_J,σ_J] | Bates CF + **COS 256 项** → IV | Feller + λ>0 |
| `VGCalibrator` | 3 参数 [σ,ν,θ] | VG CF + COS 256 项 → IV | Feller: `1-θν-σ²ν/2>0` |
| `CEVCalibrator` | 2 参数 [σ,β] | CEV 解析定价（非中心卡方） → IV | β∈(0,1) |

**BSM 隐含波动率反演** `bsm_implied_vol`: Newton + 二分混合，保单调收敛，tol 1e-10，50 次。

**统一流程**: DE 全局寻优 → LM 精修，目标 RelativeError 加权 IV 残差。

### 3.8 波动率曲面层 `models/vol_surface/`

#### SVI `svi.hpp`

**参数化**: **Gatheral 2004 raw 参数** `w(k) = a + b(ρ(k-m) + √((k-m)²+σ²))`
- a: 总方差水平，b: 斜率（≥0），ρ: 相关性（|ρ|<1），σ: 弯曲度（>0），m: 中心偏移

**无套利检查**: 蝶式套利 `g(k) = (1 - k·w'/2w)² - (w'/2)²·(w+0.25) + w''/2 ≥ 0`（Gatheral-Jacquier 2014），1001 点网格。

**参数转换**: Raw ↔ Natural ↔ Jump-Wings。

#### SSVI `ssvi.hpp`

**参数化**: **Gatheral-Jacquier 2014** `w(k,θ) = θ/2·(1 + φ(θ)·ρ·k + √((φ(θ)·k+ρ)² + 1-ρ²))`
- θ: ATM 总方差，φ(θ): 期限结构函数

**无套利条件**:
1. φ(θ) > 0
2. ∂_θ(θ·φ(θ)) > 0（Calendar 套利）
3. |ρ| < 1
4. φ(θ)·θ·(1+|ρ|) < 4（充分 Butterfly 条件，Theorem 4.2）
5. φ(θ)·(θ·(1+|ρ|) + 2·ψ·√(1-ρ²)) < 4（严格 Butterfly，Theorem 4.4），ψ=∂_θ(θ·φ(θ))/φ(θ)

**预设参数化**:
- Heston-like: `φ(θ)=η·θ^{-λ}`, 0<λ<1
- Power-law: `φ(θ)=η·θ^{-γ}`, 0<γ<1/2（保证严格无蝶式套利）

#### Dupire 局部波动率 `dupire_local_vol.hpp`

**公式**: **Dupire (1994)** `σ²_loc(K,T) = (∂C/∂T + qC + (r-q)K·∂C/∂K) / (0.5·K²·∂²C/∂K²)`

- C(K,T) 由 IV 网格双线性插值 + BSM 解析 call 价
- 偏导用中心差分：`h_K=1e-3·K`, `h_T=1e-4·max(T,0.1)`
- **精度**: O(h²) 截断，但 IV 插值误差会放大（分母 ∂²C/∂K² 数值敏感，极 OTM 退化返回 0）
- **MC 验证**: 50 步 Euler，xorshift64* + Box-Muller

### 3.9 风险管理层 `risk/`

#### 3.9.1 Greeks `risk/greeks/`

| 方法 | 算法 | 精度 | 复杂度 | 适用场景 |
|---|---|---|---|---|
| `AnalyticGreeksEngine` | BSM 闭式公式 | 机器精度 | O(1) | Vanilla BSM |
| `AADGreeksEngine` | **autodiff reverse mode**（外部 autodiff 库 `autodiff::var`） | 机器精度（一阶） | 前向 O(N)，反向 O(N)，内存 O(N) | 路径依赖/篮子 |
| `PathwiseGreeksEngine` | **Glasserman 2003 §7**：`dPrice/dθ=E[disc·dPayoff/dθ]` | MC 误差 O(1/√N) | O(N_paths) | 光滑 payoff |
| `LRGreeksEngine` | **Glasserman 2003 §7.2-7.3**：`dPrice/dθ=E[disc·Payoff·(∂log p/∂θ)]` | MC 误差 O(1/√N) | O(N_paths) | **不连续 payoff**（digital/barrier） |
| `NumericalGreeksEngine` | 有限差分（中心差分 O(h²)） | O(h²) 截断 + ε/h 舍入 | O(1) per Greek | 兜底 |
| `GreeksFactory` | Auto 派发 | — | 零开销 | 统一入口 |

**Auto 派发策略**（PHASE3_SPEC §2.2）:
- Vanilla + BSM → **Analytic**
- 不连续 payoff → **LR**（Pathwise 在 digital 上 delta=0 a.e.，AAD 在 max kink 有偏）
- 光滑无闭式 → **Pathwise**
- 路径依赖/篮子 → **AAD**
- 兜底 → **FD**

#### 3.9.2 VaR `risk/var/`

| 方法 | 算法 | 理论 | 精度/复杂度 |
|---|---|---|---|
| `HistoricalVaR` | 历史模拟 + 分位插值 + BRW 衰减 + bootstrap | √T 时间缩放 | 排序 O(N log N)；bootstrap O(n_boot×N log N) |
| `ParametricVaR` | Normal / StudentT / Cornish-Fisher | Normal: `VaR=-(μ+z_{1-c}σ√h)`；CF: `z_cf=z+(z²-1)S/6+(z³-3z)K/24-(2z³-5z)S²/36` | Normal 闭式 O(1)；StudentT Newton 100 次 + betainc |
| `MCVaR` | Full / Delta-Gamma / Delta + antithetic | Cholesky 分解协方差，`ΔR=LZ`，`Π≈δᵀΔR+½ΔRᵀΓΔR` | Cholesky O(n³)；MC 误差 O(1/√n_paths) |
| `ExpectedShortfall` | 经验 / 正态 / StudentT / CF | 正态: `ES=-μ+σφ(z)/(1-c)` | 正态闭式；StudentT Newton + betainc |
| `Backtesting` | **Kupiec POF** / **Christoffersen** / Basel 交通灯 | Kupiec LR~χ²(1)；Christoffersen LR_ind~χ²(1)，LR_cc~χ²(2)；χ² 用 Wilson-Hilferty 近似 | O(N) |

**Basel 交通灯**（250 日 99%）: Green ≤4 (×3.0)，Yellow 5-9 (×3.4-3.85)，Red ≥10 (×4.0)。

#### 3.9.3 情景分析 `risk/scenario/`

- `SensitivityAnalysis`: 单因子/多因子/方向性/最坏情形（MC 搜索）
- `StressTester`: Basel FRTB 7 场景 + 历史危机（2008 GFC / 2020 COVID）+ 情景加权 ES

#### 3.9.4 XVA `risk/xva.hpp`

**算法**（假设 τ 与 V(t) 独立）:
- **CVA** = `-(1-R_c)·∫EPE_disc(t)·dPD_c(t)`（梯形离散）
- **DVA** = `+(1-R_self)·∫ENE_disc(t)·dPD_self(t)`（IFRS 13）
- **FVA** (Hull-White symmetric) = `-∫(EPE-ENE)·s_f(t)·P_d(t)·dt`
- **BVA** = CVA + DVA + FVA
- `adjusted_price = risk_free_price + BVA`

文献: Pykhtin & Zhu (2007), Gregory (2015), Hull & White (2012/2014), Brigo Pallavicini Papatheodorou (2013)。

#### 3.9.5 Wrong-Way Risk `risk/wrong_way_risk.hpp`

| 方法 | 算法 | 理论 | 复杂度 |
|---|---|---|---|
| Hull-White 近似 | 闭式 | `CVA_wwr≈CVA_ind·exp(ρσ_V√T·Φ⁻¹(PD_T)+½(ρσ_V√T)²)` (Hull-White 2012) | O(1) |
| **Pykhtin-Zhu Copula** | 半解析 | 单因子 Gaussian Copula，**20 节点 Gauss-Hermite 求积** | O(20) 指数收敛 |
| MC 模拟 | 全路径 | 每路径相关 [M_V, M_D]，M_D 决定 τ | O(n_paths) |

#### 3.9.6 PFE + SA-CCR `risk/pfe_sa_ccr.hpp`

- **PFE**: MC 路径算 EE/ENE/EPE/EEPE/PFE（经验分位）
- **SA-CCR (BCBS d279, 2019 修订)**:
  - `RC = max(V-C, 0)`
  - `PFE_addon = mul_addon × aggregate_addon`，`mul_addon=1.4`
  - `aggregate_addon = Σ RF × max(|L-S|, 0.4×(L+S))`（supervisory floor）
  - `EAD = 1.4 × (RC + PFE_addon)`
  - MF: <1Y 用 `√(T/1Y)`

### 3.10 GPU MC `performance/gpu/`

**算法**: 每线程 1 条路径，Philox counter RNG + Box-Muller，GBM 精确解一步到期。

| 属性 | 值 |
|---|---|
| RNG | Philox4x64-10（与 CPU 同算法，位精确一致） |
| 定价 | `price = E[exp(-rT)·payoff(S_T)]`，`S_T=S·exp((r-σ²/2)T+σ√T·Z)` |
| MC 误差 | O(1/√N) |
| CPU/GPU 一致性 | 价格相对差 < 1e-12（求和顺序不同） |
| 单测容差 | 价格 1e-6，Greeks 1e-4 |
| 目标硬件 | RTX 4060（Ada Lovelace, SM 8.9, 8 GB VRAM） |
| block_size | 256（warp 32 的倍数） |

---

## 4. 算法理论分析汇总

### 4.1 收敛阶汇总

| 算法类别 | 方法 | 收敛阶 | 理论依据 |
|---|---|---|---|
| **蒙特卡洛** | 标准 MC | O(1/√N) | 中心极限定理 |
| | Sobol QMC | O((log N)^d / N) | Koksma-Hlawka 不等式 |
| | RQMC | 保持 QMC + 无偏误差估计 | Owen 随机化 |
| **SDE 离散化** | Euler-Maruyama | 强 0.5，弱 1.0 | Kloeden-Platen |
| | Milstein | 强 1.0，弱 1.0 | Kloeden-Platen |
| | Heston Full Truncation | 强 0.5，弱 1.0 | Lord et al. 2010 |
| | Heston QE | **强 1.0** | Andersen 2008 |
| | GBM Exact | ∞（无离散误差） | 闭式解 |
| | VG Exact | ∞ | Gamma subordinator |
| **PDE 时间** | Explicit Euler | 一阶 | 有稳定性限制 |
| | Implicit Euler | 一阶 | L 稳定 |
| | Crank-Nicolson | **二阶**（θ=0.5） | 无条件稳定 |
| | Rannacher Smoothing | 二阶 | 处理非光滑 payoff |
| | ADI Craig-Sneyd | **二阶**（θ=0.5） | in 't Hout & Foulon 2010 |
| | ADI Hundsdorfer-Verwer | **二阶**（最优 θ=0.5+1/√12） | Hundsdorfer & Verwer 2003 |
| | PSOR | 收敛 ω∈(1,2) | LCP 求解 |
| **树形** | CRR / Jarrow-Rudd / Tian | O(1/n)（锯齿振荡） | Cox-Ross-Rubinstein 1979 |
| | **Leisen-Reimer** | **O(1/n²)** | Peizer-Pratt 反演 |
| | Trinomial | O(1/n) | Kamrad-Ritchken 1991 |
| **傅里叶** | COS 方法 | **指数收敛 O(e^{-N})**（光滑密度） | Fang-Oosterlee 2009 |
| | FFT (Carr-Madan) | O(N log N) | Carr-Madan 1999 |
| **分数阶 ODE** | Adams-Bashforth-Moulton | min(2, 1+α) | Diethelm-Ford-Freed 2002 |
| **优化** | Levenberg-Marquardt | 近极小点二次 | Gauss-Newton 极限 |
| | Differential Evolution | 全局概率 1 收敛 | Storn-Price 1997 |
| | Nelder-Mead | 线性（慢） | 单纯形 |

### 4.2 精度基准（实测）

| 模块 | 基准 | 实测精度 |
|---|---|---|
| BSM 解析 vs benchmark | 1e-12 | ✅ |
| Heston CF vs Schoutens table | 1e-8 | ✅（Albrecher Little Trap 分支切割修正后） |
| MC 收敛阶 | -0.5 ± 0.05 | ✅ |
| Sobol QMC 方差缩减 | ≥ 10x | ✅ |
| AAD Greeks vs 解析 | 1e-10 | ✅ |
| Greeks 四法一致（Analytic/Pathwise/LR/AAD） | 1e-6 | ✅ |
| VaR 回测 p-value | > 0.05 | ✅ |
| 标定目标函数 | < 1e-6 | ✅ |
| 跨平台浮点确定性 | 位精确 | ✅（MSVC + 2× GCC） |
| CPU/GPU MC 一致性 | 价格相对差 < 1e-12 | ✅ |
| Heston CF H=0.5 退化标准 Heston | 1e-4 | ✅（Rough Heston） |
| rbergomi Hybrid vs Cholesky | < 1e-12 | ✅（b=n_steps 时退化） |

### 4.3 复杂度汇总

| 操作 | 复杂度 | 备注 |
|---|---|---|
| Cholesky 分解 | O(N³) | 一次性 |
| Thomas 算法 | O(N) | 三对角 |
| 矩阵行列式（≥4×4） | O(N³) | 部分主元高斯消元 |
| Philox RNG | O(1) | 单次 |
| Box-Muller | O(1) | 单次 |
| 相关正态生成 | O(N²) | 含矩阵-向量乘 |
| 分数阶 Adams 求解 | O(N²) | 每步核求和 |
| AAD 反向传播 | O(N) | N=算子节点数 |
| MC 定价 | O(N_paths × path_cost) | |
| Sobol 序列生成 | O(dim) | 单次 |
| Brownian Bridge | O(N) | 单资产 |
| PCA 路径构造 | O((nd)³) | Jacobi 特征分解 |
| LSMC 回归 | O(N_paths × basis²) | OLS 正规方程 |
| COS 定价 | O(N) 单 K，O(K·N) 批量 | N=n_terms |
| FFT 定价 | O(N log N) | 整个 strike 网格 |
| 1D PDE 单步 | O(N) | Thomas |
| 2D PDE ADI 单步 | O(N_x × N_v) | 两个方向 Thomas |
| PSOR 迭代 | O(iter × N) | iter≤2000 |
| 二叉树定价 | O(n²) | n=步数 |
| 三叉树定价 | O(n²) | |
| DE 优化 | O(pop × gen × f_eval) | |
| LM 优化 | O(iter × n²) Jac + O(n³) 解方程 | |
| Heston CF 积分 | O(n_quad=4096) | 梯形 |
| Rough Heston CF | O(n_steps²) | 分数阶 Adams |
| rBergomi Cholesky | O(N³) 一次性 + O(N²) 每路径 | |
| rBergomi Hybrid | O(N·b) 每路径 | b=1 时近似 O(N) |

---

## 5. 使用示例

### 5.1 BSM 欧式期权定价 + Greeks

```cpp
#include "cpphub/risk/greeks/greeks_analytic.hpp"
using namespace cpphub::v1;

Real S=100, K=100, T=1.0, r=0.05, q=0.0, sigma=0.2;
auto greeks = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, /*is_call=*/true);
// greeks.price, .delta, .gamma, .vega, .theta, .rho, .vanna, .vomma
```

### 5.2 Heston COS 定价

```cpp
#include "cpphub/pricing/fourier/cos_method.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
using namespace cpphub::v1;

HestonCFParams p{0.04, 2.0, 0.04, 0.3, -0.7, 0.05, 0.0};
auto phi = make_heston_cf(100.0, 0.05, 0.0, p, 1.0);
COSEngine cos(phi, 100.0, 0.05, 0.0, 1.0);
Real call = cos.price_call(100.0);  // 指数收敛
```

### 5.3 Heston QE 蒙特卡洛

```cpp
#include "cpphub/models/diffusion/heston_qe.hpp"
#include "cpphub/pricing/monte_carlo/mc_engine.hpp"
using namespace cpphub::v1;

HestonParams p{0.04, 2.0, 0.04, 0.3, -0.7, 100.0, 0.05, 0.0};
HestonQE heston(p);
Philox4x64 rng(42);
std::vector<Real> payoffs;
for (Size i = 0; i < 100000; ++i) {
    std::vector<Real> path(2);
    heston.generate_path(1.0, 1, path, rng);
    payoffs.push_back(std::exp(-0.05*1.0) * std::max(path[1] - 100.0, 0.0));
}
auto result = MCEngine::make_result(payoffs, /*df=*/1.0);
// result.price, .std_error, .ci_lower, .ci_upper
```

### 5.4 Sobol QMC + Brownian Bridge 路径依赖定价

```cpp
#include "cpphub/pricing/monte_carlo/qmc_engine.hpp"
#include "cpphub/pricing/monte_carlo/path_dependent_payoffs.hpp"
using namespace cpphub::v1;

QMCConfig cfg;
cfg.n_paths = 4096;
cfg.n_replicates = 16;
cfg.use_brownian_bridge = true;

auto payoff = make_asian_payoff(100.0, OptionType::Call, AsianAverageType::Arithmetic);
auto result = price_path_dependent_qmc_single(100.0, 0.2, 0.05, 0.0, 1.0, 50, payoff, cfg);
// result.price, .std_error, .variance_reduction
```

### 5.5 PDE 美式期权定价

```cpp
#include "cpphub/pricing/pde/pde_engine.hpp"
#include "cpphub/instruments/payoff/vanilla.hpp"
using namespace cpphub::v1;

PDEEngineConfig cfg;
cfg.scheme = FDMSchemeType::CrankNicolson;
cfg.rannacher_warmup = 4;  // 处理行权价 kink
PDEEngine engine(cfg);

PayOffPut put(100.0);
Real price = engine.price_american(put, 100.0, 100.0, 1.0, 0.05, 0.0, 0.2);
auto greeks = engine.greeks(put, 100.0, 100.0, 1.0, 0.05, 0.0, 0.2, /*american=*/true);
```

### 5.6 二叉树 Leisen-Reimer 定价

```cpp
#include "cpphub/pricing/tree/binomial.hpp"
#include "cpphub/instruments/payoff/vanilla.hpp"
using namespace cpphub::v1;

BinomialParams params{100.0, 100.0, 1.0, 0.05, 0.0, 0.2, 1001};
BinomialTreeEngine engine(params, BinomialType::LeisenReimer);  // O(1/n²) 收敛
PayOffCall call(100.0);
Real price = engine.price_american(call);
```

### 5.7 Heston 模型校准

```cpp
#include "cpphub/calibration/calibrator.hpp"
using namespace cpphub::v1;

std::vector<MarketQuote> quotes = /* 市场数据 */;
CalibConfig cfg;
cfg.de_generations = 200;
cfg.lm_max_iter = 200;

HestonCalibrator calibrator;
calibrator.set_market(100.0, 0.05, 0.0);
auto result = calibrator.calibrate(quotes, cfg);
// result.params = [v0, kappa, theta, sigma_v, rho]
// result.objective_value, .converged, .covariance
bool feller = HestonCalibrator::check_feller(calibrator.extract_params(result.params));
```

### 5.8 SVI 波动率曲面校准 + 无套利检查

```cpp
#include "cpphub/models/vol_surface/svi.hpp"
using namespace cpphub::v1;

SVI svi(SVIParams{0.0, 0.1, 0.0, 0.1, 0.0}, SVIParamType::Raw);
CalibConfig cfg;
auto result = svi.calibrate(strikes, maturities, implied_vols, /*forward=*/100.0, cfg);
bool no_butterfly = svi.check_butterfly_arbitrage();
```

### 5.9 VaR + ES + 回测

```cpp
#include "cpphub/risk/var/historical_var.hpp"
#include "cpphub/risk/var/parametric_var.hpp"
#include "cpphub/risk/var/expected_shortfall.hpp"
#include "cpphub/risk/var/backtesting.hpp"
using namespace cpphub::v1;

// 历史 VaR
HistoricalVaR hvar(pnl_history, 0.99, 1);
Real var = hvar.var(QuantileInterpolation::Linear);

// 参数化 VaR (Cornish-Fisher)
PortfolioStats stats = ParametricVaR::estimate_stats(pnl_history);
ParametricVaR pvar(stats, 0.99, 1);
Real cf_var = pvar.var(ParametricMethod::CornishFisher);

// ES
ExpectedShortfall es;
Real es_99 = es.normal_es(0, 1, 0.99);  // ≈ 2.665

// 回测
KupiecPOF kupiec;
auto bt = kupiec.test(var_series, realized_losses, 0.99);
// bt.p_value > 0.05 → 不拒绝模型
```

### 5.10 CVA + Wrong-Way Risk

```cpp
#include "cpphub/risk/xva.hpp"
#include "cpphub/risk/wrong_way_risk.hpp"
using namespace cpphub::v1;

XVAConfig xva_cfg;
xva_cfg.recovery_counterparty = 0.40;
WWRConfig wwr_cfg;
wwr_cfg.rho_ww = 0.3;  // WWR 相关性

// Hull-White 近似
auto cva_wwr = compute_cva_wwr_hwu(profile, pd_c, pd_self, xva_cfg, wwr_cfg, sigma_V);

// Pykhtin-Zhu Copula (20 节点 Gauss-Hermite)
auto cva_copula = compute_cva_wwr_copula(profile, pd_c, pd_self, xva_cfg, wwr_cfg, sigma_V);
```

### 5.11 Python 绑定

```python
import cpphub

# BSM 定价
price = cpphub.bsm_call_price(100.0, 100.0, 1.0, 0.05, 0.0, 0.2)

# Heston MC
result = cpphub.heston_mc_price(
    S=100.0, K=100.0, T=1.0, r=0.05, q=0.0,
    v0=0.04, kappa=2.0, theta=0.04, sigma=0.3, rho=-0.7,
    n_paths=100000, seed=42
)
# result.price, result.std_error

# Greeks
greeks = cpphub.bsm_greeks(100.0, 100.0, 1.0, 0.05, 0.0, 0.2, is_call=True)
# greeks.delta, greeks.gamma, greeks.vega, ...

# VaR
var = cpphub.historical_var(pnl_list, confidence=0.99)
```

---

## 6. 文档关联

| 文档 | 内容 | 与本文档关系 |
|---|---|---|
| `docs/architecture/ARCHITECTURE.md` | 分层架构、模块契约、数据流、ADR 摘要 | 设计 rationale |
| `docs/phases/phase{1-4}/PHASEx_SPEC.md` | 阶段执行规格书 | 验收标准 |
| `docs/decisions/ADR_INDEX.md` | 14 项架构决策记录 | 决策依据 |
| `docs/audit/AUDIT_CHECKLIST.md` | 4 阶段加权审计 | 质量门禁 |
| `BUILD_PLAN.md` | 12 周/4 阶段路线图 | 里程碑 |
| `TRACEABILITY_REPORT.md` | 48 项技术声明溯源 | 零幻觉验证 |
| `docs/tasks/CALIBRATION_IMPROVEMENT_TODO.md` | 校准模块 5 项改进方向 | v1.1+ 维护 |

---

**文档版本**: 1.0
**覆盖范围**: 全部模块（core / numerics / instruments / models / pricing / monte_carlo / calibration / risk / performance）
**理论分析**: 收敛阶 / 精度 / 复杂度 三维度全覆盖
**示例**: 11 个典型用例（BSM / Heston COS / Heston QE MC / Sobol QMC / PDE 美式 / Leisen-Reimer / 校准 / SVI / VaR / CVA-WWR / Python）
