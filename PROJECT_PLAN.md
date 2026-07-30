> ⚠️ **STATUS: DEPRECATED (2026-07-29)**
>
> 本文档为早期规划草稿，与权威文档存在以下冲突：
> - 命名空间：本文档用 `cpp_hub`，权威文档用 `cpphub`
> - Python 绑定：本文档用 `pybind11`，权威文档用 `nanobind`
> - Phase 划分：本文档 6 Phase，权威文档 4 Phase
> - 目录结构：本文档平铺式，权威文档 `include/cpphub/` + `src/` 分层
> - API 风格：本文档链式 Builder，权威文档简单构造
>
> **权威文档以以下为准**：
> - `BUILD_PLAN.md` - 总体构建计划
> - `docs/architecture/ARCHITECTURE.md` - 架构设计
> - `docs/phases/phase{1-4}/PHASEx_SPEC.md` - 阶段执行规格书
> - `docs/decisions/ADR_INDEX.md` - 架构决策记录
>
> 本文档仅保留作历史参考，不再维护。所有新内容写入权威文档。

# Cpp_Hub: 本地高性能金融计算C++脚本库 - 构建计划文档（已弃用）

> 基于三本经典金融工程C++教材分析制定：
> 1. **《C++ Design Patterns and Derivatives Pricing》** (Joshi) - 设计模式与衍生品定价
> 2. **《Modeling Derivatives in C++》** (Duffy) - 衍生品建模与数值方法
> 3. **《金融工程师的C++》** (Duffy) - C++金融工程基础与最佳实践

---

## 1. 项目愿景与目标

### 1.1 核心目标
构建一个**高性能、模块化、可扩展**的本地C++金融计算库，支持：
- **衍生品定价**：欧式/美式/亚式/障碍/期权组合等
- **数值方法**：Monte Carlo、PDE/有限差分、树形模型、准蒙特卡洛
- **风险管理**：Greeks计算、VaR、ES、情景分析
- **利率模型**：短利率模型、HJM、LMM、SABR
- **随机过程**：GBM、Heston、SABR、跳跃扩散、CEV
- **高性能计算**：SIMD、多线程、内存池、表达式模板

### 1.2 设计原则（综合三本书精华）
| 原则 | 来源 | 实施方式 |
|------|------|----------|
| **开闭原则** | Joshi Ch.3,10 | Bridge模式+Factory模式，PayOff/Parameter开放扩展封闭修改 |
| **桥接模式** | Joshi Ch.4 | `PayOffBridge`/`ParameterBridge` 分离接口与实现 |
| **虚拟构造函数** | Joshi Ch.4 | `clone()` 模式实现多态拷贝 |
| **三法则/五法则** | Joshi Ch.4, Duffy Ch.4 | 显式实现拷贝/移动/赋值/析构 |
| **const正确性** | Duffy Ch.4 | 全接口const限定，引用传递避免拷贝 |
| **表达式模板** | Duffy Ch.5 | 矩阵/向量运算零开销抽象 |
| **策略模式** | Joshi全书 | 随机过程/数值方法/求解器可插拔 |
| **RAII + 智能指针** | 现代C++ | `unique_ptr`/`shared_ptr`管理生命周期 |

---

## 2. 库架构设计

### 2.1 模块层级结构
```
Cpp_Hub/
├── core/                    # 核心基础设施
│   ├── memory/              # 内存池、对象池、SIMD分配器
│   ├── math/                # 数学工具：特殊函数、插值、积分
│   ├── random/              # 随机数引擎、分布、准随机序列
│   ├── linalg/              # 线性代数：矩阵、向量、Cholesky、SVD
│   └── utils/               # 日期、日历、工具函数
│
├── instruments/             # 金融工具定义
│   ├── payoff/              # PayOff层级（Bridge模式）
│   │   ├── PayOff.hpp
│   │   ├── PayOffBridge.hpp
│   │   ├── VanillaPayOff.hpp (Call/Put/Digital/DoubleDigital)
│   │   ├── ExoticPayOff.hpp (Asian/Barrier/Lookback/Cliquet)
│   │   └── PayOffFactory.hpp (注册表模式)
│   ├── option/              # 期权合约
│   │   ├── VanillaOption.hpp
│   │   ├── ExoticOption.hpp
│   │   └── OptionFactory.hpp
│   ├── bond/                # 固收产品
│   ├── swap/                # 互换
│   └── structured/          # 结构性产品
│
├── models/                  # 随机过程与模型
│   ├── process/             # 过程基类与实现
│   │   ├── StochasticProcess.hpp
│   │   ├── GBM.hpp
│   │   ├── Heston.hpp
│   │   ├── SABR.hpp
│   │   ├── Bates.hpp (Jump-Diffusion)
│   │   ├── CEV.hpp
│   │   └── MultiFactorProcess.hpp
│   ├── ir/                  # 利率模型
│   │   ├── ShortRateModel.hpp (Vasicek/CIR/Hull-White/G2++)
│   │   ├── HJM.hpp
│   │   └── LMM.hpp (Libor Market Model)
│   └── volatility/          # 波动率曲面
│       ├── VolSurface.hpp
│       ├── SABRCalibration.hpp
│       └── LocalVol.hpp (Dupire)
│
├── pricing/                 # 定价引擎
│   ├── engine/              # 定价引擎基类
│   │   ├── PricingEngine.hpp
│   │   ├── AnalyticEngine.hpp (Black-Scholes等闭式解)
│   │   ├── MonteCarloEngine.hpp
│   │   ├── PDEngine.hpp (显式/隐式/Crank-Nicolson)
│   │   ├── TreeEngine.hpp (CRR/JR/Tian/Leisen-Reimer)
│   │   └── FourierEngine.hpp (COS/CONV/FFT)
│   ├── mc/                  # Monte Carlo专用
│   │   ├── PathGenerator.hpp
│   │   ├── VarianceReduction.hpp (Control Variate/Antithetic/Importance)
│   │   ├── QuasiRandom.hpp (Sobol/Halton/Faure)
│   │   ├── BrownianBridge.hpp
│   │   └── RNG.hpp (Mersenne Twister/PCG/Philox)
│   ├── pde/                 # PDE求解器
│   │   ├── FDMGrid.hpp
│   │   ├── FDMScheme.hpp (Explicit/Implicit/CN/Douglas)
│   │   └── BoundaryCondition.hpp
│   └── analytic/            # 解析解公式库
│       ├── BlackScholes.hpp
│       ├── Black76.hpp
│       ├── Bachelier.hpp
│       ├── HestonAnalytic.hpp
│       └── SABRAnalytic.hpp
│
├── risk/                    # 风险管理
│   ├── greeks/              # Greeks计算
│   │   ├── AnalyticGreeks.hpp
│   │   ├── FiniteDifferenceGreeks.hpp
│   │   ├── PathwiseGreeks.hpp (自动微分/路径法)
│   │   └── LikelihoodRatioGreeks.hpp
│   ├── var/                 # VaR/ES
│   │   ├── HistoricalVaR.hpp
│   │   ├── ParametricVaR.hpp
│   │   ├── MonteCarloVaR.hpp
│   │   └── ExpectedShortfall.hpp
│   └── scenario/            # 情景分析/压力测试
│
├── calibration/             # 模型标定
│   ├── Optimizer.hpp (Levenberg-Marquardt/Nelder-Mead/DE)
│   ├── CalibrationHelper.hpp
│   └── ObjectiveFunction.hpp
│
├── performance/             # 高性能优化层
│   ├── simd/                # AVX2/AVX-512内联函数
│   ├── parallel/            # 线程池、任务图、OpenMP/TBB
│   ├── memory/              # 对齐分配器、对象池
│   └── expression/          # 表达式模板(Eigen风格)
│
├── io/                      # I/O与集成
│   ├── serialization/       # JSON/Protobuf/MessagePack
│   ├── excel/               # XLW/XLL接口
│   ├── python/              # pybind11绑定
│   └── market_data/         // 行情数据接口
│
├── test/                    # 测试体系
│   ├── unit/                # 单元测试
│   ├── integration/         # 集成测试
│   ├── benchmark/           # 性能基准
│   └── validation/          # 数值验证(对比已知解)
│
├── examples/                # 示例程序
│   ├── pricing_examples.cpp
│   ├── calibration_examples.cpp
│   ├── risk_examples.cpp
│   └── performance_demo.cpp
│
└── benchmarks/              # 基准测试数据
```

### 2.2 核心设计模式应用矩阵

| 组件 | 设计模式 | 关键类 | 书籍参考 |
|------|----------|--------|----------|
| PayOff层级 | Bridge + Virtual Constructor | `PayOff`, `PayOffBridge`, `PayOffFactory` | Joshi Ch.3,4,10 |
| 参数项 | Bridge | `Parameter`, `ParameterConstant`, `ParameterPiecewise` | Joshi Ch.4 |
| 随机过程 | Strategy + Template Method | `StochasticProcess`, `GBM`, `Heston` | Duffy Ch.2, Joshi |
| 定价引擎 | Strategy + Bridge | `PricingEngine`, `MCEngine`, `PDEEngine` | Joshi全书 |
| 路径生成器 | Builder | `PathGenerator`, `MultiPathGenerator` | Duffy Ch.2 |
| 方差缩减 | Decorator | `ControlVariate`, `Antithetic`, `ImportanceSampling` | Duffy Ch.2.5 |
| 准随机数 | Strategy | `SobolSequence`, `HaltonSequence` | Duffy Ch.2.4 |
| 矩阵运算 | Expression Template | `Matrix`, `Vector`, `Expression` | Duffy Ch.5 |
| 标定器 | Template Method | `Calibrator`, `ObjectiveFunction` | Duffy + 现代优化 |

---

## 3. 关键技术实现细节

### 3.1 PayOff Bridge模式（Joshi Ch.4核心）

```cpp
// core/payoff/PayOff.hpp
class PayOff {
public:
    virtual ~PayOff() = default;
    virtual double operator()(double spot) const = 0;
    virtual std::unique_ptr<PayOff> clone() const = 0;  // 虚拟构造函数
    virtual std::string name() const = 0;
};

// core/payoff/PayOffBridge.hpp
class PayOffBridge {
    std::unique_ptr<PayOff> ptr_;
public:
    PayOffBridge() = default;
    PayOffBridge(const PayOff& payoff) : ptr_(payoff.clone()) {}
    PayOffBridge(const PayOffBridge& other) : ptr_(other.ptr_->clone()) {}
    PayOffBridge(PayOffBridge&&) = default;
    PayOffBridge& operator=(const PayOffBridge& other) {
        if (this != &other) ptr_ = other.ptr_->clone();
        return *this;
    }
    PayOffBridge& operator=(PayOffBridge&&) = default;
    double operator()(double spot) const { return (*ptr_)(spot); }
    const PayOff& get() const { return *ptr_; }
};

// core/payoff/PayOffFactory.hpp - 注册表模式
class PayOffFactory {
    using Creator = std::function<std::unique_ptr<PayOff>(const nlohmann::json&)>;
    std::unordered_map<std::string, Creator> registry_;
public:
    static PayOffFactory& instance();
    void registerPayOff(const std::string& name, Creator creator);
    std::unique_ptr<PayOff> create(const std::string& name, const nlohmann::json& params);
    PayOffBridge createBridge(const std::string& name, const nlohmann::json& params);
};
```

### 3.2 Monte Carlo引擎（Duffy Ch.2 + Joshi优化）

```cpp
// pricing/mc/PathGenerator.hpp
template<typename Process>
class PathGenerator {
    Process process_;
    size_t nSteps_, nPaths_;
    std::unique_ptr<RNG> rng_;
    bool useBrownianBridge_;
public:
    // 生成单条路径
    std::vector<double> generatePath(double spot, double T) const;
    // 批量生成（SIMD友好）
    void generatePaths(double spot, double T, double* output) const;
    // 多因子相关路径
    void generateCorrelatedPaths(const std::vector<double>& spots, double T,
                                  const Matrix& chol, double* output) const;
};

// pricing/mc/MonteCarloEngine.hpp
template<typename PayoffType, typename ProcessType>
class MonteCarloEngine : public PricingEngine {
    ProcessType process_;
    size_t nPaths_, nSteps_;
    std::unique_ptr<VarianceReduction> vr_;
    bool useQuasiRandom_;
public:
    double price(const VanillaOption& option) const override;
    std::vector<double> pricePathDependent(const ExoticOption& option) const;
    Greeks greeks(const Option& option, GreekFlags flags) const;
};

// pricing/mc/VarianceReduction.hpp - Decorator模式
class VarianceReduction {
public:
    virtual ~VarianceReduction() = default;
    virtual void apply(std::vector<double>& payoffs, const Path& path) const = 0;
};

class AntitheticVariates : public VarianceReduction { ... };
class ControlVariate : public VarianceReduction { ... };
class ImportanceSampling : public VarianceReduction { ... };
class MomentMatching : public VarianceReduction { ... };
```

### 3.3 Heston模型与解析特征函数定价

```cpp
// models/process/Heston.hpp
struct HestonParams {
    double v0, kappa, theta, sigma, rho;
    // Feller条件检查: 2*kappa*theta > sigma*sigma
    bool checkFeller() const { return 2*kappa*theta > sigma*sigma; }
};

class HestonProcess : public StochasticProcess {
    HestonParams params_;
public:
    // 特征函数（用于COS/FFT定价）
    std::complex<double> characteristicFunction(std::complex<double> u, double tau) const;
    // 精确模拟
    void simulateExact(double v0, double dt, double& v_out, double& int_v_out) const;
    // QE方案
    void simulateQE(double v0, double dt, double& v_out) const;
};

// pricing/analytic/HestonAnalytic.hpp
class HestonAnalyticEngine : public AnalyticEngine {
    // COS方法
    double priceCOS(const VanillaOption& opt, const HestonParams& params, 
                    size_t N = 4096, double L = 12) const;
    // FFT方法
    double priceFFT(const VanillaOption& opt, const HestonParams& params) const;
    // 积分法（Gauss-Laguerre）
    double priceIntegration(const VanillaOption& opt, const HestonParams& params) const;
};
```

### 3.4 PDE引擎（Duffy Ch.12+现代改进）

```cpp
// pricing/pde/FDMGrid.hpp
class FDMGrid {
    std::vector<double> x_;  // 空间网格（非均匀，集中在strike附近）
    std::vector<double> t_;  // 时间网格
    Matrix u_;               // 解网格 (space x time)
public:
    // 非均匀网格生成：sinh变换集中在at-the-money
    static std::vector<double> generateNonUniformGrid(double x_min, double x_max, 
                                                      size_t N, double concentration = 1.0);
};

// pricing/pde/FDMScheme.hpp
class FDMScheme {
public:
    virtual void step(const Matrix& u_n, Matrix& u_n1, double dt) const = 0;
    virtual ~FDMScheme() = default;
};

class ExplicitEuler : public FDMScheme { ... };
class ImplicitEuler : public FDMScheme { ... };
class CrankNicolson : public FDMScheme { ... };  // 二阶精度，无条件稳定
class DouglasScheme : public FDMScheme { ... };  // 多因子ADI

// pricing/pde/PDEEngine.hpp
template<typename Scheme = CrankNicolson>
class PDEEngine : public PricingEngine {
    std::unique_ptr<Scheme> scheme_;
    BoundaryCondition bc_lower_, bc_upper_;
public:
    double price(const VanillaOption& opt, const StochasticProcess& process) const;
    Greeks computeGreeks(const VanillaOption& opt, const StochasticProcess& process) const;
    // 早期行使处理（美式期权）
    double priceAmerican(const AmericanOption& opt, const StochasticProcess& process) const;
};
```

### 3.5 高性能SIMD与并行化

```cpp
// performance/simd/MathSIMD.hpp
#if defined(__AVX2__)
#include <immintrin.h>

inline __m256d exp_avx2(__m256d x) {
    // 近似指数函数，最大相对误差 < 1e-6
    const __m256d one = _mm256_set1_pd(1.0);
    const __m256d ln2 = _mm256_set1_pd(0.6931471805599453);
    // ... 多项式近似实现
}

inline void black_scholes_simd(const double* spots, const double* strikes,
                               double vol, double rate, double time,
                               double* calls, double* puts, size_t n) {
    for (size_t i = 0; i + 3 < n; i += 4) {
        __m256d S = _mm256_loadu_pd(spots + i);
        __m256d K = _mm256_loadu_pd(strikes + i);
        // 向量化BS公式计算
        // ...
        _mm256_storeu_pd(calls + i, call_result);
        _mm256_storeu_pd(puts + i, put_result);
    }
}
#endif

// performance/parallel/ThreadPool.hpp
class ThreadPool {
    std::vector<std::thread> workers_;
    moodycamel::ConcurrentQueue<std::function<void()>> tasks_;
    std::atomic<bool> stop_{false};
public:
    template<typename F, typename... Args>
    auto enqueue(F&& f, Args&&... args) 
        -> std::future<std::invoke_result_t<F, Args...>>;
    
    // 并行for循环
    template<typename F>
    void parallel_for(size_t begin, size_t end, F&& f);
    
    // 任务图执行
    void execute(TaskGraph&& graph);
};
```

---

## 4. 开发里程碑与交付物

### Phase 1: 基础设施 (Week 1-2)
| 任务 | 交付物 | 验收标准 |
|------|--------|----------|
| CMake + Conan/vcpkg配置 | `CMakeLists.txt`, `conanfile.txt` | 单命令构建，支持Release/Debug/RelWithDebInfo |
| 核心数学库 | `core/math/`, `core/linalg/` | 单测覆盖>90%，对比Eigen/Boost验证 |
| 随机数引擎 | `core/random/` | 通过TestU01/BigCrush测试，Sobol序列验证 |
| 日期/日历系统 | `core/utils/Date.hpp`, `Calendar.hpp` | 支持ISDA日历，节假日数据可配置 |
| 内存池/SIMD工具 | `performance/memory/`, `performance/simd/` | 基准测试显示>2x加速 |

### Phase 2: 核心定价框架 (Week 3-4)
| 任务 | 交付物 | 验收标准 |
|------|--------|----------|
| PayOff体系 | `instruments/payoff/` | Factory注册新PayOff无需修改现有代码 |
| Parameter Bridge | `models/parameter/` | 支持常数/分段常数/函数对象参数 |
| 随机过程基类 | `models/process/StochasticProcess.hpp` | GBM/Heston/CEV通过同一接口驱动MC |
| Monte Carlo引擎 | `pricing/mc/MonteCarloEngine.hpp` | 欧式期权定价误差<0.1%，Greeks通过有限差分验证 |
| 解析引擎 | `pricing/analytic/BlackScholes.hpp` | 对比Haug(1998)基准数据全通过 |

### Phase 3: 进阶模型与数值方法 (Week 5-7)
| 任务 | 交付物 | 验收标准 |
|------|--------|----------|
| Heston模型 | `models/process/Heston.hpp` | 特征函数验证，COS/FFT定价对比文献值 |
| SABR模型 | `models/volatility/SABR.hpp` | 隐含波动率拟合误差<1bp |
| PDE引擎 | `pricing/pde/` | 欧式期权收敛阶2，美式期权早期行使边界正确 |
| 树形模型 | `pricing/tree/` | CRR/Leisen-Reimer收敛验证 |
| 准蒙特卡洛 | `pricing/mc/QuasiRandom.hpp` | Sobol序列方差缩减>10x vs伪随机 |

### Phase 4: 风险管理与标定 (Week 8-9)
| 任务 | 交付物 | 验收标准 |
|------|--------|----------|
| Greeks体系 | `risk/greeks/` | 路径法/似然比法/解析法一致性验证 |
| VaR/ES引擎 | `risk/var/` | 历史回测覆盖率在置信区间内 |
| 模型标定器 | `calibration/` | Heston/SABR标定在真实市场数据收敛 |
| 情景分析 | `risk/scenario/` | 支持历史/假设/相关性冲击场景 |

### Phase 5: 高性能优化与集成 (Week 10-11)
| 任务 | 交付物 | 验收标准 |
|------|--------|----------|
| SIMD向量化 | `performance/simd/` | 批量BS定价>5x标量版本 |
| 多线程并行 | `performance/parallel/` | 100万路径MC线性加速比>0.85 |
| Python绑定 | `io/python/` | `pip install cpp_hub` 可用，API符合Python习惯 |
| Excel加载项 | `io/excel/` | XLL加载，UDF函数在Excel中正常调用 |
| 序列化 | `io/serialization/` | JSON/Protobuf往返无损 |

### Phase 6: 测试、文档、发布 (Week 12)
| 任务 | 交付物 |
|------|--------|
| 完整测试套件 | 单测/集成/基准/验证全覆盖，CI/CD流水线 |
| API文档 | Doxygen + Sphinx生成，含数学公式推导 |
| 用户指南 | 入门教程、进阶用法、性能调优指南 |
| 基准报告 | 对比QuantLib/OpenGamma/自研基线 |
| v1.0 Release | GitHub Release + Conan Center Index提交 |

---

## 5. 质量保证体系

### 5.1 测试金字塔
```
                    ┌─────────────────┐
                    │  Validation     │  ← 对比已知解析解、文献基准、商业软件
                    │  (Accuracy)     │
            ┌───────┴───────┴───────┐
            │   Benchmark           │  ← 性能回归测试，基准数据版本化
            │  (Performance)        │
    ┌───────┴───────┬───────┴───────┐
    │ Integration   │  Integration  │  ← 端到端定价/标定/风险流程
    │  (Pricing)    │   (Risk)      │
┌───┴───┐   ┌───────┴───────┐   ┌───┴───┐
│ Unit  │   │     Unit      │   │ Unit  │  ← 单元测试覆盖>90%
│ Math  │   │   Models      │   │ Utils │
└───────┘   └───────────────┘   └───────┘
```

### 5.2 数值验证基准
| 类别 | 参考来源 | 容差 |
|------|----------|------|
| 欧式期权BS | Haug (1998) "Complete Guide to Option Pricing Formulas" | 1e-10 |
| 美式期权 | Broadie & Detemple (1996) 基准值 | 1e-4 |
| Heston | Heston (1993) 原论文表格值 / Kahl & Jäckel (2005) | 1e-6 |
| SABR隐含波动率 | Hagan et al. (2002) 展开式 | 1bp |
| 亚式期权 | Kemna & Vorst (1990) 几何亚式解析解 | 1e-6 |
| 障碍期权 | Hui (1996) / Geman & Yor (1996) | 1e-5 |

### 5.3 持续集成配置
```yaml
# .github/workflows/ci.yml
jobs:
  build_test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        config: [Debug, Release, RelWithDebInfo]
        compiler: [gcc-13, clang-17]
    steps:
      - uses: actions/checkout@v4
      - name: Install deps
        run: sudo apt-get install -y libboost-all-dev libeigen3-dev
      - name: Configure
        run: cmake -B build -DCMAKE_BUILD_TYPE=${{matrix.config}} -DCPP_HUB_BUILD_TESTS=ON
      - name: Build
        run: cmake --build build --parallel
      - name: Unit Tests
        run: ctest --test-dir build -j4 --output-on-failure
      - name: Benchmark
        if: matrix.config == 'Release'
        run: ./build/benchmarks/cpp_hub_benchmarks --benchmark_format=json > bench.json
      - name: Upload Benchmark
        uses: actions/upload-artifact@v4
        with:
          name: benchmark-${{matrix.config}}-${{matrix.compiler}}
          path: bench.json
```

---

## 6. 依赖管理与构建

### 6.1 核心依赖（最小化原则）
| 依赖 | 用途 | 版本要求 | 备选 |
|------|------|----------|------|
| **Eigen3** | 线性代数、矩阵运算 | >= 3.4 | 可选：自实现表达式模板 |
| **Boost** | 随机数、数学特殊函数、日期时间、序列化 | >= 1.83 | 部分可用标准库替代 |
| **fmt** | 格式化输出、日志 | >= 10.0 | std::format (C++20) |
| **nlohmann/json** | JSON序列化、配置 | >= 3.11 |  |
| **pybind11** | Python绑定 | >= 2.11 | nanobind |
| **CLI11** | 命令行界面 | >= 2.4 |  |
| **spdlog** | 高性能日志 | >= 1.13 |  |

### 6.2 可选依赖（性能/功能增强）
| 依赖 | 用途 | 启用条件 |
|------|------|----------|
| **Intel TBB** | 任务并行、并发容器 | `CPP_HUB_USE_TBB=ON` |
| **OpenMP** | SIMD循环并行化 | `CPP_HUB_USE_OPENMP=ON` |
| **Intel MKL** | 加速线性代数/FFT/RNG | `CPP_HUB_USE_MKL=ON` |
| **CUDA** | GPU加速Monte Carlo | `CPP_HUB_USE_CUDA=ON` |
| **QuantLib** | 互操作性验证 | `CPP_HUB_USE_QL=ON` |

### 6.3 CMake预设配置
```cmake
# CMakePresets.json
{
  "version": 6,
  "configurePresets": [
    {
      "name": "default",
      "displayName": "Default Release Build",
      "description": "Optimized build with all features",
      "generator": "Ninja",
      "binaryDir": "${sourceDir}/build/release",
      "cacheVariables": {
        "CMAKE_BUILD_TYPE": "Release",
        "CPP_HUB_BUILD_TESTS": "ON",
        "CPP_HUB_BUILD_BENCHMARKS": "ON",
        "CPP_HUB_BUILD_PYTHON": "ON",
        "CPP_HUB_USE_TBB": "ON",
        "CPP_HUB_USE_OPENMP": "ON"
      }
    },
    {
      "name": "debug",
      "inherits": "default",
      "displayName": "Debug Build",
      "cacheVariables": { "CMAKE_BUILD_TYPE": "Debug" }
    },
    {
      "name": "minimal",
      "inherits": "default",
      "displayName": "Minimal Dependencies",
      "cacheVariables": {
        "CPP_HUB_USE_TBB": "OFF",
        "CPP_HUB_USE_OPENMP": "OFF",
        "CPP_HUB_BUILD_PYTHON": "OFF"
      }
    }
  ]
}
```

---

## 7. API设计示例

### 7.1 统一的定价接口
```cpp
// 统一入口
#include <cpp_hub/pricing/EngineFactory.hpp>
#include <cpp_hub/instruments/OptionFactory.hpp>
#include <cpp_hub/models/ModelFactory.hpp>

using namespace cpp_hub;

int main() {
    // 1. 创建期权
    auto option = OptionFactory::createEuropean(
        OptionType::Call, 100.0, 100.0, 1.0  // K, S, T
    );
    
    // 2. 创建模型（从配置/标定结果）
    auto model = ModelFactory::createHeston(
        HestonParams{.v0=0.04, .kappa=2.0, .theta=0.04, .sigma=0.3, .rho=-0.7}
    );
    
    // 3. 选择定价引擎
    auto engine = EngineFactory::create<MonteCarloEngine>(model)
        .withPaths(1'000'000)
        .withSteps(100)
        .withVarianceReduction({AntitheticVariates{}, ControlVariate{}})
        .withQuasiRandom(SobolSequence{})
        .build();
    
    // 4. 定价
    double price = engine->price(option);
    Greeks greeks = engine->greeks(option, GreekFlags::All);
    
    std::cout << fmt::format("Price: {:.6f}, Delta: {:.6f}, Vega: {:.6f}\n",
                             price, greeks.delta, greeks.vega);
    
    // 5. 批量定价（SIMD加速）
    std::vector<OptionPtr> portfolio = loadPortfolio("portfolio.json");
    auto results = engine->priceBatch(portfolio);
}
```

### 7.2 Python绑定接口
```python
# pip install cpp_hub
import cpp_hub as ch

# 期权定价
opt = ch.EuropeanOption(type='call', strike=100, spot=100, maturity=1.0)
model = ch.HestonModel(v0=0.04, kappa=2.0, theta=0.04, sigma=0.3, rho=-0.7)

engine = ch.MonteCarloEngine(model, n_paths=1_000_000, n_steps=100)
engine.enable_antithetic()
engine.enable_control_variate()
engine.use_sobol()

price = engine.price(opt)
greeks = engine.greeks(opt)
print(f"Price: {price:.6f}, Delta: {greeks.delta:.6f}")

# 批量定价
portfolio = ch.load_portfolio("portfolio.json")
results = engine.price_batch(portfolio)

# 标定
calibrator = ch.HestonCalibrator(market_data)
params = calibrator.calibrate(market_vols, strikes, maturities)
```

---

## 8. 风险与对策

| 风险 | 等级 | 缓解措施 |
|------|------|----------|
| 数值稳定性（Heston特征函数积分发散） | 高 | 实现Kahl-Jäckel旋转轮廓积分，提供回退方案 |
| 模板膨胀导致编译慢/二进制大 | 中 | 显式实例化常用类型，隐藏实现细节在.cpp |
| 依赖地狱 | 低 | Conan包管理，提供header-only核心子集 |
| 性能回归 | 中 | CI集成基准测试，阈值回归自动失败 |
| 标定不收敛 | 高 | 多起始点+全局优化(DE)+局部优化(LM)混合策略 |
| 内存泄漏/UB | 高 | 启用ASan/TSan/MSan，智能指针全覆盖，Clang-Tidy全规则 |

---

## 9. 参考书籍内容映射表

| 功能模块 | Joshi书章节 | Duffy《Modeling》章节 | Duffy《C++金融工程师》章节 |
|----------|-------------|----------------------|---------------------------|
| PayOff Bridge | Ch.3, 4 | - | Ch.3 |
| Factory模式 | Ch.10 | - | - |
| Parameter Bridge | Ch.4 | - | - |
| Monte Carlo框架 | Ch.5,6 | Ch.2 | - |
| 方差缩减 | - | Ch.2.5, 2.7 | - |
| 准随机序列 | - | Ch.2.4 | - |
| Brownian Bridge | - | Ch.2.9 | - |
| 跳跃扩散/CEV | - | Ch.2.10 | - |
| Heston解析解 | - | Ch.9 | - |
| SABR/波动率曲面 | - | Ch.8, 11 | - |
| 利率模型(Vasicek/CIR/HW) | - | Ch.10 | - |
| PDE/有限差分 | - | Ch.12 | - |
| 树形模型 | - | Ch.13 | - |
| C++基础/类设计 | - | - | Ch.3,4 |
| 操作符重载/表达式模板 | - | - | Ch.5 |
| 内存管理/RAII | - | - | Ch.6 |
| 继承/多态/模板 | - | - | Ch.7,8,9 |
| 设计模式应用 | 全书 | - | Ch.11+ |

---

## 10. 立即行动项

1. **初始化仓库结构**
   ```bash
   cd F:\Cpp_Hub
   mkdir -p core/{math,random,linalg,memory,utils} instruments/{payoff,option,bond,swap} models/{process,ir,volatility} pricing/{engine,mc,pde,analytic,tree} risk/{greeks,var,scenario} calibration performance/{simd,parallel,memory,expression} io/{serialization,excel,python,market_data} test/{unit,integration,benchmark,validation} examples benchmarks
   ```

2. **配置CMake + Conan + CI**
   - 创建`CMakeLists.txt`根文件
   - 创建`conanfile.txt`依赖声明
   - 创建`.github/workflows/ci.yml`

3. **实现核心数学模块**（Week 1重点）
   - `core/math/SpecialFunctions.hpp` - 正态分布、误差函数、贝塞尔函数
   - `core/linalg/Matrix.hpp` - 表达式模板矩阵类
   - `core/random/RNG.hpp` - 统一随机数接口

4. **实现PayOff Bridge + Factory**（Week 2重点）
   - 基类、Bridge、具体PayOff、注册表工厂
   - 单元测试覆盖所有PayOff类型

---

**文档版本**: 1.0  
**创建日期**: 2026-07-29  
**基于教材**: Joshi《C++ Design Patterns and Derivatives Pricing 2nd Ed》、Duffy《Modeling Derivatives in C++》、Duffy《Financial Engineer's C++》  
**目标交付**: 12周内v1.0 Release