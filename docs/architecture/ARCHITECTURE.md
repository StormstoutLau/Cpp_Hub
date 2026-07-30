# Cpp_Hub 总体架构文档

> 版本: 1.0  
> 日期: 2026-07-29  
> 规范: SPEC v1.0 (Architecture Decision Records + Phase-gated Delivery)

---

## 1. 系统概览

### 1.1 定位与边界
```
Cpp_Hub = 轻量级定价内核 + Python绑定 + 可扩展插件架构
```
- **In-scope**: 欧式/美式/障碍/亚式期权定价、Heston/SABR/利率模型、MC/PDE/Tree/解析解引擎、Greeks(AAD/Pathwise)、模型标定
- **Out-of-scope**: 完整术语结构/日历/货币体系(用QuantLib)、行情接入/Bloomberg API、Web UI/可视化、分布式计算集群

### 1.2 核心指标 (SLA)
| 指标 | 目标 | 验收方式 |
|------|------|----------|
| 单欧式期权定价延迟 | < 1 μs (批量 < 50 ns/期权) | Google Benchmark CI 守护 |
| MC 百万路径/秒 | > 50M paths/s (AVX2) | `benchmarks/mc_benchmark.cpp` |
| 编译时间 (增量) | < 5s | CMake + Unity Build |
| Python 调用开销 | < 500 ns | nanobind 零拷贝 |
| 二进制体积 | < 10 MB | `strip --strip-all` |
| 数值精度 | 1e-12 vs 解析解 | `tests/validation/` 对标 Haug/Fang-Oosterlee |

---

## 2. 分层架构

```
┌─────────────────────────────────────────────────────────────────┐
│                      Python Binding Layer (nanobind)            │
│  cpphub._core  |  cpphub.batch  |  cpphub.calib  |  cpphub.io   │
├─────────────────────────────────────────────────────────────────┤
│                        C ABI Stable Interface                    │
│  cpphub_price_batch(), cpphub_calibrate(), cpphub_greeks_batch() │
├─────────────────────────────────────────────────────────────────┤
│                      Performance Layer (Header-only)             │
│  simd.hpp  |  rng.hpp  |  parallel.hpp  |  batch.hpp  | gpu.hpp  │
├─────────────────────────────────────────────────────────────────┤
│                      Pricing Engines (Strategy)                  │
│  AnalyticEngine  |  MCEngine  |  PDEEngine  |  TreeEngine       │
│  COSEngine  |  FourierEngine  |  LSMCEngine  |  QMCEngine       │
├─────────────────────────────────────────────────────────────────┤
│                      Models (Template Method)                    │
│  StochasticProcess  |  DiffusionModels(GBM/Heston/SABR/Bates)   │
│  ShortRateModels(Vasicek/CIR/HW/G2++)  |  LMM/HJM               │
├─────────────────────────────────────────────────────────────────┤
│                      Instruments (Bridge + Factory)              │
│  PayOff/PayOffBridge  |  VanillaOption/ExoticOption             │
│  PayOffFactory(Registry)  |  Exercise(European/American/Bermudan)│
├─────────────────────────────────────────────────────────────────┤
│                      Core Infrastructure (Header-only)           │
│  math.hpp  |  linalg.hpp  |  datetime.hpp  |  error.hpp         │
│  constants.hpp  |  types.hpp  |  memory.hpp  |  config.hpp      │
└─────────────────────────────────────────────────────────────────┘
```

### 2.1 关键设计决策 (ADR 摘要)

| ADR | 决策 | 理由 | 替代方案 |
|-----|------|------|----------|
| ADR-001 | Header-only core + 单一共享库 | 编译快、部署简单、Python轮子小 | 全静态库 / 多共享库 |
| ADR-002 | Bridge + Virtual Constructor (Joshi Ch.4) | 开闭原则、无切片拷贝、值语义 | `std::variant` / CRTP |
| ADR-003 | Factory + 静态注册模板 (Joshi Ch.10) | 运行时扩展无需重编译、单例工厂 | 手工 switch / 插件系统 |
| ADR-004 | 计数器 RNG (Philox/Threefry) + 确定性分块 | 位精确复现、SIMD友好、无锁并行 | `std::mt19937` / PCG |
| ADR-005 | 可移植 SIMD 抽象层 (simde + 内置) | 单代码路径支持 AVX2/AVX-512/NEON | 手写汇编 / Vc / Highway |
| ADR-006 | nanobind 绑定 + 批量 NumPy API | 零拷贝、编译快、二进制小 | pybind11 / ctypes / C API |
| ADR-007 | AAD (伴随模式) 统一 Greeks | 一次扫描全希腊值、无数值误差 | 有限差分 / Pathwise / LR |
| ADR-008 | COS/FFT 谱方法作为解析引擎补充 | 指数收敛、Heston/Levy 通用 | 仅解析解 / 仅 MC |

---

## 3. 模块契约

### 3.1 核心层 (`include/cpphub/core/`)

| 文件 | 职责 | 关键类型/函数 | 依赖 |
|------|------|---------------|------|
| `config.hpp` | 编译期特性检测、版本号 | `CPPHUB_VERSION`, `CPPHUB_HAS_AVX2`, `CPPHUB_HAS_OPENMP` | 无 |
| `types.hpp` | 基础类型别名 | `Real=double`, `Complex=std::complex<Real>`, `Size=size_t` | `<cstddef>` |
| `constants.hpp` | 数学常数 | `PI`, `SQRT_2`, `INV_SQRT_2PI`, `DBL_EPS` | 无 |
| `math.hpp` | 特殊函数 | `erf`, `erfc`, `cdf_normal`, `cdf_normal_inv`, `bessel_i`, `gamma` | `<cmath>` |
| `linalg.hpp` | 表达式模板线性代数 (固定尺寸, 定价专用) | `Matrix<R,C>`, `Vector<N>`, `cholesky`, `thomas_algorithm` | `simd.hpp` |

> **linalg.hpp 范围说明**: 当前仅覆盖定价所需的固定尺寸矩阵运算（编译期 `R×C`，栈分配）。未来扩展计量模块（Newey-West、Fama-MacBeth、GARCH、因子回归）时，需新增 `linalg_dynamic.hpp` 封装 Eigen3 处理动态尺寸大矩阵（运行时 `N×K`，堆分配 + BLAS 加速 + SVD/QR/LU）。详见 ADR-013 和 `docs/research/ECONOMETRICS_LANDSCAPE.md`。
| `datetime.hpp` | 日期/日历/年基 | `Date`, `Calendar`, `DayCount`, `Schedule`, `Frequency` | 无 |
| `simd.hpp` | 可移植 SIMD | `f64x4`, `f64x8`, `exp`, `log`, `sqrt`, `hsum`, `blend` | `simde` (内置) |
| `rng.hpp` | 无状态计数器 RNG | `Philox4x64`, `Threefry4x64`, `normal_simd` | `Random123` (头文件) |
| `parallel.hpp` | 确定性并行 | `ThreadPool`, `parallel_for`, `blocked_range` | OpenMP / TBB |
| `batch.hpp` | 批量计算入口 | `bsm_price_batch`, `bsm_greeks_batch`, `heston_cos_batch` | 以上所有 |
| `memory.hpp` | 内存工具 | `aligned_allocator`, `object_pool`, `stack_allocator` | `<memory>` |
| `error.hpp` | 错误码/异常 | `Result<T, ErrorCode>`, `CppHubException` | `<system_error>` |

### 3.2 定价引擎层 (`include/cpphub/pricing/`)

```cpp
// 统一引擎接口 (Strategy Pattern)
class PricingEngine {
public:
    virtual ~PricingEngine() = default;
    virtual double price(const VanillaOption&) const = 0;
    virtual Greeks greeks(const VanillaOption&) const = 0;
    virtual std::string name() const = 0;
};

// 具体引擎
class AnalyticBSEngine : public PricingEngine { ... };
class MCEngine : public PricingEngine { MCConfig cfg_; ... };
class PDEEngine : public PricingEngine { FDMScheme scheme_; ... };
class COSEngine : public PricingEngine { COSConfig cfg_; ... };
class TreeEngine : public PricingEngine { TreeType type_; ... };
```

**PayOff 接口设计（双接口方案，支持路径相关期权）**:

原设计 `operator()(double spot)` 无法表达 Asian/Barrier/Lookback 等路径相关期权。
采用双接口方案，路径相关期权使用 `PathDependentPayOff` 接口（参考 Joshi Ch.7-8）：

```cpp
// 终值 PayOff: 仅依赖到期 spot（Call/Put/Digital/DoubleDigital）
class PayOff {
public:
    virtual ~PayOff() = default;
    virtual double operator()(double spot) const = 0;
    virtual std::unique_ptr<PayOff> clone() const = 0;
    virtual std::string name() const = 0;
};

// 路径相关 PayOff: 依赖整条路径（Asian/Barrier/Lookback/Cliquet）
// 参考 Joshi Ch.7 PathDependent, Ch.8 ExoticEngine
class PathDependentPayOff {
public:
    virtual ~PathDependentPayOff() = default;
    // 传入完整路径 (times, spots), 返回现金流向量
    virtual std::vector<CashFlow> cashFlows(const Path& path) const = 0;
    // 路径相关期权可能需要多个观察时间点
    virtual std::vector<double> lookbackTimes() const = 0;
    virtual std::unique_ptr<PathDependentPayOff> clone() const = 0;
    virtual std::string name() const = 0;
};

// VanillaOption 持有 PayOffBridge (终值)
// ExoticOption 持有 PathDependentPayOffBridge (路径相关)
// MC 引擎通过 concept/if constexpr 分发到对应接口
```

**引擎选择策略**:
| 期权类型 | 推荐引擎 | 备选 |
|----------|----------|------|
| 欧式 BS | `AnalyticBSEngine` | `COSEngine` (Heston) |
| 欧式 Heston | `COSEngine(N=256)` | `AnalyticHestonCFEngine` |
| 美式/百慕大 | `PDEEngine(CN+PSOR)` | `TreeEngine(LeisenReimer)` |
| 亚式/障碍/回望 | `MCEngine(QMC+BB)` | `PDEEngine` (低维) |
| 篮子/最优/最差 | `MCEngine(MultiAsset)` | `COSEngine` (约数) |

### 3.3 模型层 (`include/cpphub/models/`)

```cpp
// Template Method: 固定算法骨架，派生类填空
class StochasticProcess {
public:
    virtual size_t dimension() const = 0;
    virtual void generate_path(Real T, Size steps, Span<Real> path, RNG& rng) const = 0;
    // 可选：特征函数 (用于 Fourier 方法)
    virtual std::complex<Real> characteristic_function(std::complex<Real> u, Real tau) const { return {}; }
};

// 具体模型
class GBM : public StochasticProcess { ... };           // 1因子
class Heston : public StochasticProcess { ... };        // 2因子
class Bates : public StochasticProcess { ... };         // Heston + Jump
class SABR : public StochasticProcess { ... };          // 随机波动率
class VarianceGamma : public StochasticProcess { ... }; // Levy
```

### 3.4 标定层 (`include/cpphub/calibration/`)

```cpp
struct CalibrationProblem {
    std::vector<MarketQuote> quotes;  // (strike, expiry, price/iv, type)
    ModelParameters initial_guess;
    WeightingScheme weights;          // vega-weighted / price-weighted
    Bounds bounds;                    // 参数边界约束
};

class Calibrator {
public:
    virtual ModelParameters calibrate(const CalibrationProblem&) = 0;
    virtual ~Calibrator() = default;
};

class LevenbergMarquardtCalibrator : public Calibrator { ... };
class DifferentialEvolutionCalibrator : public Calibrator { ... }; // 全局搜索
```

---

## 4. 数据流与并发模型

### 4.1 单期权定价流
```
Option + Model + Engine
       │
       ▼
┌──────────────────┐
│ PricingEngine::price() │
└────────┬─────────┘
         │
    ┌────┴────┐
    ▼         ▼
 Analytic   Numeric(MC/PDE/Tree)
    │         │
    ▼         ▼
 Result   StatisticsMC
            │
            ▼
       ┌────────┐
       │ Greeks │ (AAD/Pathwise/FD)
       └────────┘
```

### 4.2 批量定价流 (SIMD + 并行)
```
批量输入: spots[N], strikes[N], vols[N], rates[N], expiries[N]
       │
       ▼
┌─────────────────────────┐
│ blocked_range(0, N, 64) │  ← 固定块大小，位精确复现
└───────────┬─────────────┘
            │
    ┌───────┼───────┐
    ▼       ▼       ▼
 Thread1  Thread2  ThreadK  (每线程独立 Philox seed = base + block_id)
    │       │       │
    ▼       ▼       ▼
 SIMD循环  SIMD循环  SIMD循环  (f64x4/f64x8 向量化)
    │       │       │
    └───────┼───────┘
            ▼
       写入输出数组
```

### 4.3 确定性并行保证
```cpp
// 关键不变量：相同输入 → 位精确相同输出，与线程数无关
class DeterministicMC {
    static constexpr size_t N_BLOCKS = 64;  // 固定块数
    Philox4x64 rng_;  // seed = base_seed + block_id
    
    void parallel_price(...) {
        #pragma omp parallel for schedule(static)  // 静态调度！
        for (size_t b = 0; b < N_BLOCKS; ++b) {
            Philox4x64 block_rng(base_seed + b);
            process_block(block_rng, ...);
        }
    }
};
```

---

## 5. 部署与分发

### 5.1 构建产物
| 产物 | 目标 | 说明 |
|------|------|------|
| `libcpphub.{so,dll,dylib}` | C++ 链接 | 稳定 C ABI，版本号 `cpphub_abi_v1` |
| `cpphub_python.{so,pyd}` | Python 轮子 | `pip install cpphub` |
| `cpphub_benchmarks` | 性能回归 | CI 运行，结果上传为 artifact |

### 5.2 版本策略
- **语义化版本**: `MAJOR.MINOR.PATCH`
- **ABI 稳定**: MAJOR 变更才破坏 C ABI
- **Python 轮子**: `cpphub-X.Y.Z-cp311-cp311-manylinux_2_17_x86_64.whl`

### 5.3 支持矩阵
| 平台 | 编译器 | SIMD | Python | 状态 |
|------|--------|------|--------|------|
| Linux x86_64 | GCC 11+/Clang 14+ | AVX2/AVX-512 | 3.10-3.12 | ✅ Tier 1 |
| Windows x64 | MSVC 19.35+ | AVX2 | 3.10-3.12 | ✅ Tier 1 |
| macOS arm64 | Apple Clang 14+ | NEON | 3.10-3.12 | ✅ Tier 1 |
| Linux aarch64 | GCC 11+ | NEON | 3.10+ | 🟡 Tier 2 |

---

## 6. 质量门禁

| 门禁 | 工具 | 阈值 | 失败动作 |
|------|------|------|----------|
| 编译警告 | `-Wall -Wextra -Wpedantic -Werror` | 0 | 阻塞合并 |
| 编译选项 | 数值路径禁用 `-ffast-math`，用 `-ffp-contract=off`；CI 用 `-march=x86-64-v3` 固定目标 | 位精确复现测试通过 | 阻塞合并 |
| 静态分析 | Clang-Tidy (modernize, performance, bugprone) | 0 | 阻塞合并 |
| 单元测试 | Catch2 | 覆盖率 > 90% | 阻塞合并 |
| 数值验证 | `tests/validation/` | 相对误差 < 1e-8 (解析) / 1e-3 (MC) | 阻塞合并 |
| 性能基准 | Google Benchmark | 不得回退 > 5% | 阻塞合并 (需人工复核) |
| 内存安全 | ASan/UBSan/MSan (CI) | 0 报错 | 阻塞合并 |
| 线程安全 | TSan (CI) | 0 数据竞争 | 阻塞合并 |
| Python 绑定 | pytest + hypothesis | 100% API 覆盖 | 阻塞合并 |

> **编译选项与位精确复现的关系**：`-ffast-math` 启用 `-funsafe-math-optimizations`，允许编译器重排浮点运算（破坏结合律）、用 FMA 替换 mul+add，导致不同线程数的归约顺序产生不同结果。因此数值路径（MC/Greeks/AAD）必须禁用 `-ffast-math`。仅批量纯函数（BS batch）可单独启用并从位精确门禁豁免。`-march=native` 在 CI runner 上因 CPU 型号不固定会导致 `SIGILL` 或结果不一致，改用 `-march=x86-64-v3`（AVX2）固定目标。

---

## 7. 扩展点

| 扩展点 | 接口 | 示例 |
|--------|------|------|
| 新 PayOff | `PayOffFactory::registerPayOff(name, creator)` | `LookbackPayOff` |
| 新随机过程 | 继承 `StochasticProcess` + `ModelFactory::register` | `RoughBergomi` |
| 新定价引擎 | 继承 `PricingEngine` + `EngineFactory::register` | `DeepBSDEEngine` |
| 新标定算法 | 继承 `Calibrator` + `CalibratorFactory::register` | `BayesianCalibrator` |
| 新方差缩减 | 继承 `VarianceReduction` (Decorator) | `ImportanceSamplingVR` |
| GPU 后端 | 实现 `GpuBackend` 抽象接口 | `CudaBackend`, `HipBackend` |

---

## 8. 附录：关键文件映射表

| 功能区 | 核心文件 | 行数估计 | 来源溯源 |
|--------|----------|----------|----------|
| PayOff Bridge | `instruments/payoff/payoff_bridge.hpp` | ~80 | Joshi Ch.4 p.72-75 |
| PayOff Factory | `instruments/payoff/factory.hpp` | ~120 | Joshi Ch.10 p.175-184 |
| StatisticsMC | `pricing/monte_carlo/statistics.hpp` | ~100 | Joshi Ch.5 p.84-90 |
| Sobol QMC | `pricing/monte_carlo/quasi_monte_carlo.hpp` | ~200 | Duffy Ch.2.4 p.78-84 |
| Heston COS | `pricing/analytic/heston_cos.hpp` | ~150 | Fang & Oosterlee (2009) + QuantKernel |
| Thomas 算法 | `pricing/pde/thomas_algorithm.hpp` | ~80 | Duffy Ch.18 p.323-326 |
| AAD Greeks | `risk/greeks/aad_greeks.hpp` | ~200 | pricer/adjoint.hpp |
| Philox RNG | `core/rng.hpp` | ~150 | pricer/rng.hpp + Random123 |
| SIMD 层 | `core/simd.hpp` | ~300 | pricer/simd.hpp + simde |
| nanobind 绑定 | `python/bindings.cpp` | ~300 | OptiCore 模式 |

---

**文档维护**: 每个 Phase 结束必须更新本文档对应章节，并记录 ADR 变更