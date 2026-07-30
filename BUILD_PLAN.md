# Cpp_Hub: 本地高性能金融计算C++脚本库 - 构建计划文档

> **数据源**：三本经典教材 + 10+ 现代开源量化库深度分析
> - 《C++ Design Patterns and Derivatives Pricing》 (Joshi) - 设计模式/架构核心
> - 《Modeling Derivatives in C++》 (Duffy) - 数值方法/工程实践
> - 《金融工程师的C++》 (Duffy) - C++金融工程基础/最佳实践
> - 现代开源库：QuantKernel, OptiCore, pricer, mape, qflib, QuantPricer, HPC-Pricing-Kernel, QuantSphere

---

## 1. 项目定位与核心目标

### 1.1 定位
**Cpp_Hub** 是一个**模块化、高性能、可嵌入**的本地C++金融计算库，定位为：
- ❌ 不是 QuantLib 的替代品（太重、太复杂、编译慢）
- ✅ 是**轻量级定价内核** + **Python绑定** + **可扩展插件架构**
- ✅ 面向：量化研究原型验证、风险实时计算、结构化产品定价、教学演示

### 1.2 核心指标（参考现代库基准）

| 指标 | 目标 | 参考基准 |
|------|------|----------|
| 单欧式期权定价延迟 | < 1 μs | OptiCore: <1ms/10k ≈ 0.1μs |
| Monte Carlo 百万路径/秒 | > 50M paths/s | QuantKernel: ~20M/s (batch) |
| 编译时间 (增量) | < 5s | Header-only + 模块化 |
| Python调用开销 | < 500ns | nanobind 零拷贝 |
| 二进制体积 | < 10 MB | 单一共享库 |

### 1.3 设计原则（融合三本书精华 + 现代实践）

| 原则 | 来源 | 实施 |
|------|------|------|
| **Bridge + Virtual Constructor** | Joshi Ch.3-4 | `PayOffBridge`/`ParameterBridge` 分离接口实现 |
| **Strategy + Template Method** | Joshi 全书 | 定价引擎/随机过程/求解器可插拔 |
| **Factory + Registry** | Joshi Ch.10 | `PayOffFactory`/`ModelFactory` 运行时注册 |
| **Expression Templates** | Duffy Ch.5 | 矩阵/向量运算零开销抽象 |
| **const-correctness + RAII** | Duffy Ch.4 | 全接口 const、智能指针管理生命周期 |
| **SIMD + Counter-based RNG** | pricer/QuantKernel | 可移植 SIMD 层 + Philox/Threefry 无状态 RNG |
| **Deterministic Parallelism** | pricer/mape | 固定分块 + 确定性种子 → 位精确复现 |
| **Header-only Core + C ABI** | QuantKernel/mape | 核心仅头文件、稳定 C API 边界、Python 绑定层分离 |

---

## 2. 库架构设计

### 2.1 目录结构

```
Cpp_Hub/
├── cmake/                    # CMake 模块、工具链、依赖管理
├── include/
│   └── cpphub/
│       ├── core/             # 核心基础设施
│       │   ├── config.hpp            # 编译配置、SIMD检测、版本
│       │   ├── types.hpp             # 基础类型别名、Real/Complex
│       │   ├── constants.hpp         # 数学常数
│       │   ├── memory.hpp            # 对齐分配器、对象池、栈分配器
│       │   ├── simd.hpp              # 可移植 SIMD 抽象层 (AVX2/AVX-512/NEON)
│       │   ├── rng.hpp               # 无状态计数器 RNG (Philox/Threefry)
│       │   ├── math.hpp              # 特殊函数、插值、积分、求根
│       │   ├── linalg.hpp            # 表达式模板矩阵/向量 (Eigen-lite)
│       │   ├── datetime.hpp          # 日期、日历、频率、DCF
│       │   └── error.hpp             # 错误码、异常层级
│       │
│       ├── instruments/      # 金融工具定义
│       │   ├── payoff/
│       │   │   ├── payoff.hpp            # PayOff 基类 + clone() 虚拟构造
│       │   │   ├── payoff_bridge.hpp     # Bridge 模式封装
│       │   │   ├── vanilla.hpp           # Call/Put/Digital/DoubleDigital
│       │   │   ├── exotic.hpp            # Asian/Barrier/Lookback/Cliquet
│       │   │   ├── composite.hpp         # 组合/篮子/最优/最差/价差
│       │   │   └── factory.hpp           # 注册表模式工厂
│       │   ├── option/
│       │   │   ├── vanilla_option.hpp
│       │   │   ├── exotic_option.hpp
│       │   │   └── exercise.hpp          # European/American/Bermudan
│       │   ├── fixed_income/
│       │   │   ├── bond.hpp
│       │   │   ├── swap.hpp
│       │   │   └── cashflow.hpp
│       │   └── structured/
│       │
│       ├── models/             # 随机过程与模型
│       │   ├── process.hpp             # StochasticProcess 基类 (Template Method)
│       │   ├── diffusion/
│       │   │   ├── gbm.hpp
│       │   │   ├── heston.hpp
│       │   │   ├── sabr.hpp
│       │   │   ├── bates.hpp           # Jump-diffusion
│       │   │   ├── cev.hpp
│       │   │   └── variance_gamma.hpp
│       │   ├── ir/
│       │   │   ├── short_rate.hpp      # Vasicek/CIR/HullWhite/G2++
│       │   │   ├── hjm.hpp
│       │   │   └── lmm.hpp
│       │   ├── vol_surface/
│       │   │   ├── vol_surface.hpp
│       │   │   ├── sabr_calibration.hpp
│       │   │   ├── svi.hpp
│       │   │   └── dupire_local_vol.hpp
│       │   └── factory.hpp             # ModelFactory 注册表
│       │
│       ├── pricing/            # 定价引擎 (Strategy 模式)
│       │   ├── engine.hpp              # PricingEngine 基类
│       │   ├── analytic/
│       │   │   ├── black_scholes.hpp
│       │   │   ├── black76.hpp
│       │   │   ├── bachelier.hpp
│       │   │   ├── heston_cf.hpp       # 特征函数 + COS/FFT/CONV
│       │   │   ├── sabr_hagan.hpp
│       │   │   └── variance_gamma_cf.hpp
│       │   ├── monte_carlo/
│       │   │   ├── mc_engine.hpp
│       │   │   ├── path_generator.hpp
│       │   │   ├── rng_factory.hpp
│       │   │   ├── variance_reduction.hpp  # Antithetic/ControlVariate/Importance
│       │   │   ├── quasi_monte_carlo.hpp   # Sobol/Halton/Faure
│       │   │   ├── brownian_bridge.hpp
│       │   │   ├── longstaff_schwartz.hpp  # LSMC 美式/百慕大
│       │   │   └── greeks_mc.hpp           # Pathwise/LR/AAD
│       │   ├── pde/
│       │   │   ├── fdm_grid.hpp
│       │   │   ├── fdm_scheme.hpp        # Explicit/Implicit/CN/Douglas/ADI
│       │   │   ├── boundary.hpp
│       │   │   ├── thomas_algorithm.hpp
│       │   │   └── psor_solver.hpp       # 美式早期行使
│       │   ├── tree/
│       │   │   ├── binomial.hpp          # CRR/JR/Tian/Leisen-Reimer
│       │   │   └── trinomial.hpp
│       │   └── fourier/
│       │       ├── cos_method.hpp
│       │       ├── fft_convolution.hpp
│       │       └── hilbert_transform.hpp
│       │
│       ├── risk/               # 风险管理
│       │   ├── greeks/
│       │   │   ├── analytic_greeks.hpp
│       │   │   ├── fd_greeks.hpp
│       │   │   ├── pathwise_greeks.hpp
│       │   │   ├── likelihood_ratio_greeks.hpp
│       │   │   └── aad_greeks.hpp        # 伴随自动微分
│       │   ├── var/
│       │   │   ├── historical_var.hpp
│       │   │   ├── parametric_var.hpp
│       │   │   ├── mc_var.hpp
│       │   │   └── expected_shortfall.hpp
│       │   └── scenario/
│       │       ├── stress_test.hpp
│       │       └── sensitivity.hpp
│       │
│       ├── calibration/        # 模型标定
│       │   ├── optimizer.hpp           # Levenberg-Marquardt/Nelder-Mead/DE
│       │   ├── objective.hpp
│       │   └── helpers.hpp
│       │
│       ├── performance/        # 高性能优化层
│       │   ├── parallel.hpp            # 线程池、任务图、确定性并行
│       │   ├── batch.hpp               # 批量计算 API (SIMD 友好)
│       │   └── gpu.hpp                 # CUDA/HIP 可选后端 (预留)
│       │
│       └── python/             # Python 绑定 (nanobind)
│           ├── bindings.cpp
│           ├── numpy_interop.hpp
│           └── pandas_interop.hpp
│
├── src/
│   ├── core/                 # 核心实现 (仅需编译的部分)
│   ├── models/
│   ├── pricing/
│   └── python/               # nanobind 模块入口
│
├── tests/
│   ├── unit/                 # 单元测试
│   ├── integration/          # 集成测试
│   ├── validation/           # 数值验证 (对比解析解/基准值)
│   └── benchmark/            # 性能基准
│
├── benchmarks/               # 独立基准程序
├── examples/                 # 示例程序
├── python/                   # Python 包装包
│   ├── cpphub/
│   ├── pyproject.toml
│   └── setup.py
│
├── third_party/              # 子模块依赖
│   ├── nanobind/
│   ├── catch2/
│   ├── fmt/
│   ├── simde/                # 可移植 SIMD
│   └── random123/            # Philox/Threefry RNG
│
├── CMakeLists.txt
├── vcpkg.json
├── conanfile.py
├── BUILD_PLAN.md
└── README.md
```

### 2.2 核心类图（关键模式）

```cpp
// ============ Bridge + Virtual Constructor (Joshi Ch.4) ============
class PayOff {
public:
    virtual ~PayOff() = default;
    virtual double operator()(double spot) const = 0;
    virtual std::unique_ptr<PayOff> clone() const = 0;  // Virtual Constructor
    virtual std::string name() const = 0;
};

class PayOffBridge {
    std::unique_ptr<PayOff> ptr_;
public:
    PayOffBridge() = default;
    PayOffBridge(const PayOff& p) : ptr_(p.clone()) {}
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

// ============ Strategy: PricingEngine ============
class PricingEngine {
public:
    virtual ~PricingEngine() = default;
    virtual double price(const VanillaOption& opt) const = 0;
    virtual Greeks greeks(const VanillaOption& opt) const = 0;
    virtual std::string name() const = 0;
};

// 具体引擎
class AnalyticBSEngine : public PricingEngine { ... };
class MCEngine : public PricingEngine { ... };
class PDEEngine : public PricingEngine { ... };
class COSEngine : public PricingEngine { ... };

// ============ Template Method: StochasticProcess ============
class StochasticProcess {
public:
    virtual ~StochasticProcess() = default;
    virtual size_t dimension() const = 0;
    virtual void generatePath(double T, size_t steps, 
                              Span<double> path, RNG& rng) const = 0;
    // 可选：特征函数、转移密度等
    virtual std::complex<double> characteristicFunction(
        std::complex<double> u, double tau) const { return {}; }
};

class GBM : public StochasticProcess { ... };
class Heston : public StochasticProcess { ... };
```

---

## 3. 关键模块实现规格

### 3.1 核心基础设施

#### 3.1.1 可移植 SIMD 层 (`core/simd.hpp`)
```cpp
// 基于 simde + 编译器内置，零依赖
namespace cpphub::simd {
    using f64x4 = simde::native<double, 4>;  // AVX2: 4x double
    using f64x8 = simde::native<double, 8>;  // AVX-512: 8x double
    
    // 基础算术
    CPPHUB_FORCE_INLINE f64x4 add(f64x4 a, f64x4 b) { return a + b; }
    CPPHUB_FORCE_INLINE f64x4 mul(f64x4 a, f64x4 b) { return a * b; }
    CPPHUB_FORCE_INLINE f64x4 exp(f64x4 x);    // 向量化 exp
    CPPHUB_FORCE_INLINE f64x4 log(f64x4 x);    // 向量化 log
    CPPHUB_FORCE_INLINE f64x4 sqrt(f64x4 x);   // 向量化 sqrt
    
    // 归约
    CPPHUB_FORCE_INLINE double hsum(f64x4 x);  // 水平求和
    
    // 掩码/混合
    CPPHUB_FORCE_INLINE f64x4 blend(f64x4 a, f64x4 b, f64x4 mask);
}
```
> **参考**：`pricer/simd.hpp`, `QuantKernel` batch kernels

#### 3.1.2 无状态计数器 RNG (`core/rng.hpp`)
```cpp
// Philox-4x64 / Threefry-4x64 (Random123)
// 关键特性：draw(i) = f(seed, i) 纯函数 → 完全可并行、确定性、SIMD 友好
namespace cpphub::rng {
    struct Philox4x64 {
        using state_type = std::array<uint64_t, 2>;  // key
        using counter_type = std::array<uint64_t, 4>;
        
        explicit Philox4x64(uint64_t seed = 0) : key_{seed, seed ^ 0x9E3779B97F4A7C15} {}
        
        // 生成 4 个 u64 (可直接转 double)
        std::array<uint64_t, 4> operator()(counter_type ctr) const noexcept;
        
        // SIMD 批量生成
        void generate(f64x4* out, size_t count, uint64_t base_ctr) const noexcept;
    };
    
    // Box-Muller 正态分布 (SIMD 版)
    void normal(Philox4x64 rng, f64x4* out, size_t count) noexcept;
}
```
> **参考**：`pricer/rng.hpp`, `QuantKernel` batch RNG

#### 3.1.3 表达式模板线性代数 (`core/linalg.hpp`)
```cpp
// Eigen-lite: 仅实现金融需要的固定/小尺寸矩阵运算
namespace cpphub::linalg {
    template<size_t R, size_t C> class Matrix;
    template<size_t N> using Vector = Matrix<N, 1>;
    
    // 表达式模板：M = A * B + C 无临时对象
    template<typename E> class MatrixExpr;
    
    // Cholesky 分解 (相关矩阵)
    template<size_t N> void cholesky(Matrix<N,N>& L, const Matrix<N,N>& A);
    
    // Thomas 算法 (三对角) - PDE 必需
    template<size_t N> void thomasAlgorithm(
        const Vector<N>& a, const Vector<N>& b, const Vector<N>& c,
        Vector<N>& d, Vector<N>& x);
}
```
> **参考**：Duffy Ch.5, `pricer/simd_mc.hpp`

---

### 3.2 PayOff 体系 (Bridge + Factory)

```cpp
// instruments/payoff/vanilla.hpp
class CallPayOff : public PayOff {
    double strike_;
public:
    explicit CallPayOff(double K) : strike_(K) {}
    double operator()(double S) const override { return std::max(S - strike_, 0.0); }
    std::unique_ptr<PayOff> clone() const override { return std::make_unique<CallPayOff>(*this); }
    std::string name() const override { return "Call"; }
};

// instruments/payoff/exotic.hpp
// 路径相关期权使用 PathDependentPayOff 接口 (非 PayOff)
// 参考 Joshi Ch.7 PathDependent, Ch.8 ExoticEngine
class AsianPayOff : public PathDependentPayOff {
    double strike_; bool is_geometric_;  // 几何平均有解析解
    std::vector<double> lookbackTimes_;  // 观察时间点
public:
    AsianPayOff(double K, bool geom, std::vector<double> times)
        : strike_(K), is_geometric_(geom), lookbackTimes_(std::move(times)) {}
    // 返回现金流: 路径均值 vs strike
    std::vector<CashFlow> cashFlows(const Path& path) const override;
    std::vector<double> lookbackTimes() const override { return lookbackTimes_; }
    std::unique_ptr<PathDependentPayOff> clone() const override {
        return std::make_unique<AsianPayOff>(*this);
    }
};

class BarrierPayOff : public PathDependentPayOff {
    double strike_, barrier_; int barrier_type_;  // 8种类型
    std::vector<double> lookbackTimes_;
public:
    // 需要路径最大/最小值判断是否触障
    std::vector<CashFlow> cashFlows(const Path& path) const override;
    // ...
};

// instruments/payoff/factory.hpp - 注册表模式 (Joshi Ch.10)
class PayOffFactory {
    using Creator = std::function<std::unique_ptr<PayOff>(const nlohmann::json&)>;
    std::unordered_map<std::string, Creator> registry_;
public:
    static PayOffFactory& instance();
    void registerPayOff(const std::string& name, Creator creator);
    std::unique_ptr<PayOff> create(const std::string& name, const nlohmann::json& params) const;
    PayOffBridge createBridge(const std::string& name, const nlohmann::json& params) const;
};

// 静态注册宏
#define REGISTER_PAYOFF(name, type) \
    namespace { \
        struct name##_registrar { \
            name##_registrar() { \
                PayOffFactory::instance().registerPayOff(#name, \
                    [](const json& p) { return std::make_unique<type>(p); }); \
            } \
        } name##_reg; \
    }

// 内置注册
REGISTER_PAYOFF(Call, CallPayOff)
REGISTER_PAYOFF(Put, PutPayOff)
REGISTER_PAYOFF(Digital, DigitalPayOff)
REGISTER_PAYOFF(Asian, AsianPayOff)
```

---

### 3.3 Monte Carlo 引擎 (高性能核心)

```cpp
// pricing/monte_carlo/mc_engine.hpp
struct MCConfig {
    size_t n_paths = 100'000;
    size_t n_steps = 100;
    bool antithetic = true;
    bool control_variate = true;
    bool quasi_random = false;  // Sobol
    size_t seed = 0;
    int n_threads = 0;  // 0 = auto
};

class MCEngine : public PricingEngine {
    MCConfig config_;
    std::unique_ptr<VarianceReduction> vr_;
    
public:
    explicit MCEngine(MCConfig cfg = {}) : config_(cfg) {
        if (config_.antithetic) vr_ = std::make_unique<AntitheticVR>();
        if (config_.control_variate) vr_ = std::make_unique<ControlVariateVR>(vr_.release());
    }
    
    double price(const VanillaOption& opt) const override {
        const auto& process = opt.process();
        const auto& payoff = opt.payoff();
        double spot = opt.spot();
        double T = opt.expiry();
        
        // 确定性分块并行
        const size_t n_blocks = 64;  // 固定块数，位精确复现
        const size_t paths_per_block = (config_.n_paths + n_blocks - 1) / n_blocks;
        
        double sum = 0.0, sum_sq = 0.0;
        
        #pragma omp parallel for reduction(+:sum,sum_sq) schedule(static)
        for (size_t block = 0; block < n_blocks; ++block) {
            Philox4x64 rng(config_.seed + block);  // 每块独立种子
            double block_sum = 0.0, block_sq = 0.0;
            
            for (size_t p = 0; p < paths_per_block; ++p) {
                double path_sum = 0.0;
                // 路径生成 (SIMD 向量化)
                // ...
                double payoff_val = payoff(path_sum);
                if (vr_) payoff_val = vr_->adjust(payoff_val, path_sum);
                block_sum += payoff_val;
                block_sq += payoff_val * payoff_val;
            }
            sum += block_sum;
            sum_sq += block_sq;
        }
        
        double mean = sum / config_.n_paths;
        double std_err = std::sqrt((sum_sq / config_.n_paths - mean * mean) / config_.n_paths);
        return mean * std::exp(-opt.rate() * T);  // 贴现
    }
    
    Greeks greeks(const VanillaOption& opt) const override {
        // Pathwise / LR / AAD 实现
        // 复用同一路径，仅调整 payoff 求导
    }
};
```

> **关键优化参考**：
> - `pricer/parallel_simd.hpp` - 确定性分块 + SIMD
> - `QuantKernel` batch API - 零开销批量入口
> - `mape` - 每线程独立 RNG 流

---

### 3.4 解析/谱方法引擎 (高精度、低延迟)

```cpp
// pricing/analytic/heston_cf.hpp - Heston 特征函数 + COS 方法
class HestonCFEngine : public PricingEngine {
    struct Config {
        size_t N = 256;       // COS 级数项数
        double L = 12.0;      // 积分截断范围
        bool use_fft = false; // FFT 加速
    } cfg_;
    
    // Heston 特征函数 (复数版本，支持 SIMD)
    std::complex<double> charFn(std::complex<double> u, double tau,
                                double v0, double kappa, double theta,
                                double sigma, double rho) const;
    
public:
    double price(const VanillaOption& opt) const override {
        // COS 方法: 价格 = e^{-rT} * Re[ sum_{k=0}^{N-1} ' Re[ phi(k) * exp(-ikx0) ] ]
        // 指数收敛，N=256 即达机器精度
        // 参考: Fang & Oosterlee (2009), QuantKernel COS batch
    }
};

// pricing/fourier/cos_method.hpp - 通用 COS 引擎
template<typename CharFn>
class COSEngine : public PricingEngine {
    CharFn char_fn_;
    size_t N_;
    double L_;
    
public:
    double price(const VanillaOption& opt) const override {
        // 通用 COS 实现，接受任意特征函数
        // 支持: Heston, Bates, VG, CGMY, SABR (近似), 黑斯
    }
};
```

> **参考**：`QuantKernel` (heston_price_cf_batch, cos_method_fang_oosterlee_price_batch), `hpc-stochastic-pricing-kernel` (COS engine)

---

### 3.5 PDE/FDM 引擎 (美式/障碍/局部波动率)

```cpp
// pricing/pde/fdm_scheme.hpp
enum class FDMScheme { Explicit, Implicit, CrankNicolson, Douglas, ADI_CraigSneyd };

class FDMEngine : public PricingEngine {
    struct GridConfig {
        size_t n_spot = 400;
        size_t n_time = 2000;
        double spot_max_mult = 4.0;  // S_max = strike * mult
        FDMScheme scheme = FDMScheme::CrankNicolson;
    } grid_;
    
    // 非均匀网格：sinh 变换集中在 ATM
    std::vector<double> buildNonUniformGrid(double S0, double K, double sigma, double T) const;
    
    // Thomas 算法 O(N) 求解三对角
    void solveTridiagonal(const std::vector<double>& a,
                          const std::vector<double>& b,
                          const std::vector<double>& c,
                          std::vector<double>& d,
                          std::vector<double>& x) const;
    
    // PSOR 投影法 - 美式早期行使
    void psorStep(std::vector<double>& V, const std::vector<double>& payoff,
                  double omega, double tol, int max_iter) const;
};
```

> **参考**：Duffy Ch.12, `QuantPricer` FDM Crank-Nicolson + Thomas, `QuantSphere` PDE PSOR

---

### 3.6 自动微分 Greeks (AAD)

```cpp
// risk/greeks/aad_greeks.hpp
// 使用逆模式 AD (伴随法) - 一次扫描得所有 Greek
// 参考: pricer/adjoint.hpp, mape (dual numbers)

namespace cpphub::ad {
    // 采用 autodiff 库 (MIT, header-only), 不自研 Tape
    // ADR-007 Revised: std::function 类型擦除阻止内联, 性能比成熟库慢 10-50x
    #include <autodiff/forward/dual.hpp>
    #include <autodiff/reverse/var.hpp>
    using namespace autodiff;
    // 前向模式: dual (单个 Greek)
    // 反向模式: var (全 Greeks 一次扫描)
    // 自定义接口适配层: GreeksEngine 内部使用, 不跨 Python 边界
}

// 使用示例：Black-Scholes AAD Greeks
Greeks computeBSGreeksAAD(const VanillaOption& opt) {
    using autodiff::var;
    var S = opt.spot();
    var K = opt.strike();
    var r = opt.rate();
    var sigma = opt.vol();
    var T = opt.expiry();

    // Black-Scholes 公式用 var 计算
    var d1 = (log(S/K) + (r + sigma*sigma*0.5)*T) / (sigma*sqrt(T));
    var d2 = d1 - sigma*sqrt(T);
    var price = S * cdf(d1) - K * exp(-r*T) * cdf(d2);

    autodiff::derive(price, autodiff::wrt(S, sigma, r, T));  // 一次扫描全 Greeks

    return { .delta = S.derivative(), .vega = sigma.derivative(), ... };
}
```

---

## 4. Python 绑定架构 (nanobind)

```cpp
// python/bindings.cpp
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/optional.h>

namespace nb = nanobind;
using namespace cpphub;

NB_MODULE(cpphub, m) {
    m.doc() = "CppHub: High-performance quantitative finance library";
    
    // ===== 核心类型 =====
    nb::class_<PayOffBridge>(m, "PayOff")
        .def(nb::init<const std::string&, const nb::dict&>(), "name"_a, "params"_a)
        .def("__call__", &PayOffBridge::operator(), "spot"_a);
    
    nb::class_<VanillaOption>(m, "VanillaOption")
        .def(nb::init<PayOffBridge, double, double, double, double>(),
             "payoff"_a, "spot"_a, "strike"_a, "rate"_a, "vol"_a, "expiry"_a)
        .def_rw("spot", &VanillaOption::spot_)
        .def_rw("rate", &VanillaOption::rate_)
        .def_rw("vol", &VanillaOption::vol_);
    
    // ===== 定价引擎 =====
    nb::class_<PricingEngine>(m, "PricingEngine")
        .def("price", &PricingEngine::price, "option"_a)
        .def("greeks", &PricingEngine::greeks, "option"_a);
    
    nb::class_<AnalyticBSEngine, PricingEngine>(m, "AnalyticBSEngine")
        .def(nb::init<>());
    
    nb::class_<MCEngine, PricingEngine>(m, "MCEngine")
        .def(nb::init<MCConfig>(), "config"_a = MCConfig{});
    
    nb::class_<COSEngine, PricingEngine>(m, "COSEngine")
        .def(nb::init<HestonParams, COSEngine::Config>());
    
    // ===== 批量 API (NumPy 零拷贝) =====
    m.def("bsm_price_batch", [](nb::ndarray<nb::numpy, double> spots,
                                 nb::ndarray<nb::numpy, double> strikes,
                                 nb::ndarray<nb::numpy, double> rates,
                                 nb::ndarray<nb::numpy, double> vols,
                                 nb::ndarray<nb::numpy, double> expiries,
                                 const char* opt_type) {
        // 直接指针访问，SIMD 向量化
        size_t n = spots.shape(0);
        auto out = nb::ndarray<nb::numpy, double>::create(n);
        bsm_price_batch_impl(spots.data(), strikes.data(), rates.data(),
                             vols.data(), expiries.data(), out.mutable_data(), n, opt_type);
        return out;
    }, "spots"_a, "strikes"_a, "rates"_a, "vols"_a, "expiries"_a, "opt_type"_a);
    
    // 批量 Greeks
    m.def("bsm_greeks_batch", &bsm_greeks_batch_impl, ...);
    
    // ===== 模型标定 =====
    nb::class_<HestonCalibrator>(m, "HestonCalibrator")
        .def(nb::init<>())
        .def("calibrate", &HestonCalibrator::calibrate, "market_data"_a, "config"_a = CalibConfig{});
}
```

> **参考**：`OptiCore` (nanobind, 4× faster compile), `QuantKernel` (C ABI + ctypes), `qflib` (Python C API)

---

## 5. 构建系统与依赖管理

### 5.1 CMake 配置要点

```cmake
# CMakeLists.txt 关键片段
cmake_minimum_required(VERSION 3.25)
project(CppHub VERSION 1.0.0 LANGUAGES CXX)

# C++20 模块支持 (MSVC 19.35+, GCC 13+, Clang 16+)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 编译优化
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(
        $<$<CXX_COMPILER_ID:MSVC>:/O2 /arch:AVX2 /fp:precise /GL>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3 -march=x86-64-v3 -ffp-contract=off -flto>
    )
    # 注意: 禁用 -ffast-math 以保证位精确复现(ADR-004).
    # -march=x86-64-v3 固定 AVX2 目标, 避免 -march=native 在 CI runner 上的不一致.
    # 批量纯函数可单独 target_compile_options(... -ffast-math) 豁免.
    # OpenMP SIMD (不需要 threading runtime)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/openmp:experimental> 
                            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-fopenmp-simd>)
    endif()
endif()

# 头文件库目标
add_library(cpphub_core INTERFACE)
target_include_directories(cpphub_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(cpphub_core INTERFACE cxx_std_20)

# 编译单元 (仅需编译的部分：Python 绑定、SIMD 运行时检测等)
add_library(cpphub SHARED src/core/simd_dispatch.cpp src/python/bindings.cpp)
target_link_libraries(cpphub PUBLIC cpphub_core)
target_link_libraries(cpphub PRIVATE nanobind::nanobind ${OpenMP_CXX_FLAGS})

# Python 包
find_package(Python3 REQUIRED COMPONENTS Development)
nanobind_add_module(cpphub_python src/python/bindings.cpp)
target_link_libraries(cpphub_python PRIVATE cpphub)

# 依赖管理 (vcpkg/Conan)
# vcpkg.json: nanobind, catch2, fmt, simde, random123
```

### 5.2 依赖矩阵

| 依赖 | 用途 | 版本 | 备选 |
|------|------|------|------|
| nanobind | Python 绑定 | ≥2.0 | pybind11 (较慢) |
| Catch2 | 测试 | ≥3.5 | GoogleTest |
| fmt | 格式化/日志 | ≥10.0 | std::format (C++20) |
| simde | 可移植 SIMD | ≥0.8 | 手写内联汇编 |
| Random123 | Philox/Threefry RNG | ≥1.15 | 自实现 |
| nlohmann/json | 配置/序列化 | ≥3.11 | - |
| OpenMP | 线程并行/SIMD | 系统提供 | TBB/HPX |
| CUDA/HIP | GPU 后端 (可选) | ≥12.0 | - |

---

## 6. 开发路线图 (里程碑)

### 版本与 Phase 对应关系

| 版本 | Phase | 内容 | 预计工期 |
|------|-------|------|----------|
| **v1.0** | Phase 1 + Phase 2 (削减后) | BS/Heston + Analytic/MC/PDE/Tree + Sobol QMC + 基础 Greeks | 6-8 周 |
| **v1.1** | Phase 3 (+ Phase 2 推迟内容) | 进阶模型 (SABR/IR/Levy) + COS/FFT/LSMC + AAD Greeks + VaR/ES + 标定 + 波动率曲面 | 3-4 周 |
| **v2.0** | Phase 4 | GPU/分布式/Excel/云原生/多平台发布 | 3-4 周 |

### Phase 1: 核心内核 (Week 1-3) ✅ **v1.0 核心 - 最小可行产品**
- [ ] `core/` 基础设施：SIMD、RNG、数学、线性代数、日期时间
- [ ] `instruments/payoff/`：Bridge 模式、Vanilla/Exotic、Factory 注册表
- [ ] `models/process/`：GBM、Heston、基类模板方法
- [ ] `pricing/analytic/`：Black-Scholes、Bachelier、Black76 (解析解 + Greeks)
- [ ] `pricing/monte_carlo/`：基础 MC 引擎 + Antithetic + 确定性并行
- [ ] Python 绑定：`VanillaOption`、`AnalyticBSEngine`、`MCEngine`、批量 API
- [ ] 单元测试 + 数值验证 (对比已知解)
- [ ] CI: Windows/Linux/macOS + MSVC/GCC/Clang

### Phase 2: 进阶数值方法 (Week 4-6) ✅ **v1.0 核心 - 削减后**
- [ ] **Heston 模型完善** (特征函数 + QE/Exact 模拟)
- [ ] **PDE/FDM**：Crank-Nicolson + Thomas + PSOR (美式/障碍)
- [ ] **树形模型**：CRR/JR/Tian + Leisen-Reimer (高阶收敛)
- [ ] **准蒙特卡洛**：Sobol 序列 + Brownian Bridge
- [ ] **基础方差缩减**：Control Variate、Moment Matching
- ⏸️ **v1.1 推迟**: SABR/Bates/VG/CEV 模型、利率模型、COS/FFT/CONV、LSMC、标定框架、波动率曲面

### Phase 3: 风险与标定 (Week 7-10) 📋 **v1.1**
- [ ] **Phase 2 推迟内容**: SABR/Bates/VG/CEV/利率模型 + COS/FFT/LSMC + 高级方差缩减
- [ ] **Greeks 体系**：解析/有限差分/Pathwise/LR/AAD 统一接口 (AAD 采用 autodiff 库)
- [ ] **VaR/ES**：历史模拟/参数/MC
- [ ] **标定框架**：Levenberg-Marquardt + DE + 多目标 (价格 + Greeks)
- [ ] **波动率曲面**：SVI/SSVI 无套利参数化 + Dupire 局部波动率

### Phase 4: 生产级特性 (Week 11-14) 🚀 **v2.0**
- [ ] **GPU 后端** (CUDA)：MC 路径生成、COS/FFT、PDE 求解器
- [ ] **分布式计算**：MPI/ZeroMQ 多节点 MC
- [ ] **持久化**：Arrow/Parquet 行情数据、模型序列化
- [ ] **Excel 加载项** (XLL)：xlwings / xlOil 集成
- [ ] **云原生**：gRPC 微服务、Prometheus 监控、OpenTelemetry 追踪
- [ ] **文档与示例**：Sphinx + Jupyter Notebooks
- [ ] **多平台发布**：Wheel + Conda Forge + vcpkg + Conan + Docker
- [ ] **性能回归测试**：基准守护 (PR > 5% 回归失败)

### v1.0 Scope 削减方案 (2026-07-29 评审后新增)

> **背景**: 12 周单人开发无法覆盖原 Phase 1-4 全部 scope（QuantLib 团队 20+ 年累积）。
> 按减法原则，将原 scope 分为"v1.0 保留 / v1.1 推迟 / v2.0 砍掉"三档。

| 模块 | v1.0 保留 (Phase 1+2 核心) | v1.1 推迟 (Phase 3) | v2.0 砍掉 (Phase 4) |
|------|------|------|------|
| **随机过程** | GBM, Heston | SABR, Bates, VG, CEV | HJM, LMM |
| **定价引擎** | Analytic, MC, PDE, Tree | COS, FFT, LSMC | CONV |
| **Greeks** | Analytic, FD | Pathwise, LR, AAD (autodiff) | - |
| **VaR** | 历史, 参数 | MC VaR, ES | 情景分析 |
| **标定** | Levenberg-Marquardt | DE | SVI/SSVI/Dupire |
| **Python** | nanobind + Wheel | vcpkg 或 Conan | Conda Forge |
| **SIMD** | CPU AVX2 + OpenMP | - | GPU/CUDA |
| **分布式** | - | - | MPI, ZeroMQ |
| **Excel** | - | - | XLL (xlOil) |
| **部署** | - | - | gRPC, K8s, Prometheus |

**v1.0 交付目标**: 6-8 周可交付一个真实可用的 C++ 定价内核 + Python 绑定，覆盖 BS/Heston + Analytic/MC/PDE/Tree + 基础 Greeks + 历史/参数 VaR + LM 标定。

**v1.1 交付目标**: +3-4 周，覆盖进阶引擎 (COS/FFT/LSMC) + AAD Greeks + MC VaR/ES。

**v2.0**: GPU/分布式/Excel/部署，视需求和资源再议。

---

## 7. 验收标准与测试策略

### 7.1 数值验证矩阵

| 模型/方法 | 基准 | 容差 | 测试用例来源 |
|-----------|------|------|--------------|
| BS 欧式 | 解析解 | 1e-12 | Haug (2007) |
| BS Greeks | 解析解 | 1e-10 | 自动微分验证 |
| Heston COS | Fang & Oosterlee 表1 | 1e-8 | 论文基准值 |
| Heston MC | COS 价格 | 1e-3 (相对) | 同一参数集 |
| 美式 Put (树) | Broadie-Detemple | 1e-4 | 基准文献 |
| 美式 Put (PDE) | 树形/Leisen-Reimer | 1e-4 | 交叉验证 |
| 亚式期权 (MC) | 几何平均解析解 | 1e-3 | 控制变量验证 |
| 障碍期权 | Haug 解析公式 | 1e-4 | 8 种障碍类型 |
| SABR IV | Hagan 公式 | 1e-8 | 标定往返测试 |

### 7.2 性能基准 (CI 守护)

```cpp
// benchmarks/pricing_benchmark.cpp
BENCHMARK(BS_Analytic_Scalar) {
    for (auto _ : state) {
        double p = bsm_price(spot, strike, rate, vol, expiry, 'C');
        benchmark::DoNotOptimize(p);
    }
}

BENCHMARK(BS_Analytic_Batch_AVX2) {
    alignas(64) double spots[1000], strikes[1000], rates[1000], vols[1000], expiries[1000];
    // 初始化...
    for (auto _ : state) {
        bsm_price_batch_avx2(spots, strikes, rates, vols, expiries, out, 1000, 'C');
    }
}

BENCHMARK(MC_GBM_1M_Paths) {
    MCEngine engine({.n_paths=1'000'000, .n_steps=1, .antithetic=true});
    VanillaOption opt(...);
    for (auto _ : state) {
        double p = engine.price(opt);
        benchmark::DoNotOptimize(p);
    }
}

// CI 阈值 (GitHub Actions ubuntu-latest, 2 vCPU)
// BS_Analytic_Batch_AVX2: > 50M prices/s
// MC_GBM_1M_Paths: > 20M paths/s
```

### 7.3 对抗性测试

```cpp
// tests/adversarial/test_edge_cases.cpp
TEST_CASE("Extreme parameters") {
    // 零利率、零波动率、极短/极长期限
    CHECK(bsm_price(100, 100, 0.0, 0.0, 1e-6, 'C') == 0.0);
    CHECK(bsm_price(100, 100, 0.05, 0.0, 1.0, 'C') == 100*(1-exp(-0.05)));
    
    // 深度实值/虚值
    CHECK(bsm_price(1000, 100, 0.05, 0.2, 1.0, 'C') == doctest::Approx(900).epsilon(1e-6));
    CHECK(bsm_price(10, 100, 0.05, 0.2, 1.0, 'P') == doctest::Approx(90).epsilon(1e-6));
    
    // 高波动率
    CHECK(bsm_price(100, 100, 0.05, 5.0, 1.0, 'C') > 0);
    
    // 负利率 (Bachelier)
    CHECK(bachelier_price(100, 100, -0.01, 0.01, 1.0, 'C') >= 0);
}

TEST_CASE("Idempotency - same input same output") {
    MCEngine engine({.n_paths=100000, .seed=42});
    VanillaOption opt(...);
    double p1 = engine.price(opt);
    double p2 = engine.price(opt);
    CHECK(p1 == p2);  // 位精确复现
}

TEST_CASE("Thread count independence") {
    // 1, 2, 4, 8, 16 线程结果位精确相同
    for (int t : {1,2,4,8,16}) {
        MCEngine engine({.n_paths=1000000, .n_threads=t, .seed=123});
        CHECK(engine.price(opt) == reference_price);
    }
}
```

---

## 8. 与现有库的集成策略

| 场景 | 推荐方案 |
|------|----------|
| **需要完整术语结构/日历/货币** | 集成 QuantLib (仅链接所需模块) |
| **需要 Bloomberg/Refinitiv 数据** | 使用 `blpapi`/`ema` C++ SDK，CppHub 仅做计算 |
| **需要 Web 可视化** | Python 层用 `FastAPI` + `Plotly Dash`，CppHub 通过 nanobind 调用 |
| **需要数据库存储** | `SQLite` + `sqlite_orm` 或 `DuckDB` 分析查询 |
| **需要消息队列** | `ZeroMQ`/`Redis Streams` 分发定价任务 |
| **容器化部署** | 多阶段 Docker: build (C++) → runtime (Python + .so) |

---

## 9. 许可证与分发

- **核心库**: Apache-2.0 (允许商业闭源嵌入)
- **Python 包**: PyPI 发布 `pip install cpphub`
- **Conda**: `conda install -c conda-forge cpphub`
- **Docker**: `ghcr.io/yourname/cpphub:latest` (包含 Python 环境)
- **系统包**: vcpkg/Conan recipe 提供

---

## 10. 立即行动清单 (Day 1-3)

```bash
# 1. 初始化仓库
mkdir -p Cpp_Hub && cd Cpp_Hub
git init

# 2. 创建目录结构
mkdir -p include/cpphub/{core,instruments/{payoff,option},models/{process,diffusion,ir},pricing/{analytic,monte_carlo,pde,tree,fourier},risk/{greeks,var},calibration,performance,python}
mkdir -p src/{core,models,pricing,python} tests/{unit,integration,validation,benchmark} benchmarks examples python/cpphub third_party

# 3. 添加子模块依赖
git submodule add https://github.com/nanobind/nanobind third_party/nanobind
git submodule add https://github.com/catchorg/Catch2 third_party/catch2
git submodule add https://github.com/fmtlib/fmt third_party/fmt
git submodule add https://github.com/simd-everywhere/simde third_party/simde
git submodule add https://github.com/DES-Distributed-Environment/Random123 third_party/random123

# 4. 编写 CMakeLists.txt + vcpkg.json + conanfile.py
# 5. 实现 core/simd.hpp + core/rng.hpp (最高优先级，后续所有性能依赖此二)
# 6. 实现 instruments/payoff/ (Bridge + Factory)
# 7. 实现 pricing/analytic/black_scholes.hpp (标量 + SIMD 批量)
# 8. 编写 Python 绑定骨架 + 第一个可运行示例
# 9. 配置 GitHub Actions CI (Windows/Linux/macOS, MSVC/GCC/Clang)
# 10. 跑通 `pip install -e python/` 并验证 `import cpphub; cpphub.bsm_price_batch(...)`
```

---

## 附录：关键参考代码映射表

| 功能 | 参考库 | 关键文件 |
|------|--------|----------|
| Bridge PayOff | Joshi 书 + QuantKernel | `payoff_bridge.hpp`, `QuantKernel/include/payoff.hpp` |
| SIMD 抽象层 | pricer | `pricer/simd.hpp` |
| Counter-based RNG | pricer + Random123 | `pricer/rng.hpp`, `Random123/Philox.h` |
| 确定性并行 MC | pricer/mape | `pricer/parallel_simd.hpp`, `mape/src/engine/monte_carlo.cpp` |
| COS 方法 | QuantKernel + hpc-pricing | `QuantKernel/src/cos_method.cpp`, `hpc-stochastic-pricing-kernel/src/engines/cos_engine.hpp` |
| Heston 特征函数 | QuantKernel | `QuantKernel/src/models/heston.cpp` |
| AAD Greeks | pricer + mape | `pricer/adjoint.hpp`, `mape/src/engine/autodiff.cpp` |
| LSMC | QuantPricer + hpc-pricing | `QuantPricer/src/mc/lsmc_engine.cpp`, `hpc-stochastic-pricing-kernel/src/engines/lsmc_engine.hpp` |
| nanobind 绑定 | OptiCore | `OptiCore/src/bindings.cpp` |
| 表达式模板线性代数 | Duffy Ch.5 + Eigen | `Duffy/linalg/`, `Eigen/Core` |
| SVI 标定 | QuantLib + pricer | `pricer/svi.hpp`, `QuantLib/ql/termstructures/volatility/equityfx/svi.hpp` |

---

**文档版本**: 1.0  
**创建日期**: 2026-07-29  
**维护者**: Cpp_Hub Core Team  
**下一步**: 执行「立即行动清单」第 1-10 项，建立可编译、可测试、可打包的最小骨架