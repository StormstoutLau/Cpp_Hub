# 架构决策记录 (ADR) 索引

> 格式遵循 MADR (Markdown Architectural Decision Records)  
> 状态：`Proposed` → `Accepted` → `Superseded` / `Deprecated`  
> 所有 ADR 必须在实施前创建，Accepted 后不可变更，仅能 Superseded

---

## ADR 列表

| 编号 | 标题 | 状态 | 日期 | 关联 Phase |
|------|------|------|------|------------|
| ADR-001 | Header-only Core + 单一共享库 | Accepted | 2026-07-29 | Phase 1 |
| ADR-002 | Bridge + Virtual Constructor (PayOff) | Accepted | 2026-07-29 | Phase 1 |
| ADR-003 | Factory + 静态注册模板 | Accepted | 2026-07-29 | Phase 1 |
| ADR-004 | 计数器 RNG (Philox/Threefry) + 确定性分块并行 | Accepted | 2026-07-29 | Phase 1 |
| ADR-005 | 可移植 SIMD 抽象层 | Accepted | 2026-07-29 | Phase 1 |
| ADR-006 | nanobind Python 绑定 + 批量 NumPy API | Accepted | 2026-07-29 | Phase 1 |
| ADR-007 | AAD (伴随模式) 统一 Greeks | Accepted | 2026-07-29 | Phase 3 |
| ADR-008 | COS/FFT 谱方法作为解析引擎补充 | Accepted | 2026-07-29 | Phase 2 |
| ADR-009 | C ABI 稳定边界 + 版本化 | Accepted | 2026-07-31 | Phase 4 LITE |
| ADR-010 | GPU 后端可选编译 + 运行时回退 | Accepted | 2026-07-31 | Phase 4 LITE |
| ADR-011 | 分布式计算: MPI Master-Worker + 确定性聚合 | Proposed | 2026-07-29 | Phase 4 |
| ADR-012 | xlOil XLL 加载项 (异步 UDF + 缓存) | Proposed | 2026-07-29 | Phase 4 |
| ADR-013 | 双层线性代数架构 (固定尺寸 + 动态尺寸) | Accepted | 2026-08-04 | Phase 6 |
| ADR-014 | 标定 (Calibration) vs 估计 (Estimation) 的分离 | Accepted | 2026-07-29 | Phase 6 |
| ADR-015 | 证伪统计量模块边界 (通用 vs 模块特定) | Accepted | 2026-08-12 | Phase 7A |
| ADR-016 | 金融时间序列实施边界 (18 项) | Accepted | 2026-08-15 | Phase 7B |
| ADR-017 | 时序模块命名空间 (`cpphub::v1::timeseries`) | Accepted | 2026-08-15 | Phase 7B |

---

## ADR-001: Header-only Core + 单一共享库

**状态**: Accepted  
**日期**: 2026-07-29  
**决策者**: 架构组  
**关联 Phase**: 1

### 背景
需要平衡编译速度、部署简便性、Python 轮子大小、C++ 链接灵活性。

### 决策
- **核心层** (`include/cpphub/core/`, `instruments/`, `models/`, `pricing/`, `risk/`) 设计为 **header-only**，仅依赖标准库和头文件库 (simde, Random123, fmt)
- **仅编译单元**: `src/core/simd_dispatch.cpp` (运行时 SIMD 分发), `src/python/bindings.cpp` (nanobind 模块)
- 产出单一共享库 `libcpphub.{so,dll,dylib}` + Python 模块 `cpphub_python.{so,pyd}`
- CMake `INTERFACE` 库 `cpphub_core` 导出包含目录与编译特性

### 理由
1. **编译速度**: 头文件库无需链接，增量编译仅重新编译变更的 TU
2. **部署简单**: Python wheel 仅包含 1 个 `.so` + 头文件 (可选)
3. **零开销抽象**: 模板/内联函数编译期展开，LTO 跨 TU 优化
4. **ABI 稳定**: 仅导出 C 风格符号 + nanobind 模块入口，C++ 类不跨 ABI 边界

### 后果
- 模板代码膨胀风险 → 显式实例化常用类型 (`Matrix<4,4>`, `Vector<3>`)
- 编译期依赖传播 → 模块化设计，最小化头文件包含
- 调试符号大 → Release 构建 `strip --strip-all`

### 替代方案评估
| 方案 | 优点 | 缺点 | 结论 |
|------|------|------|------|
| 全静态库 | 无 ABI 问题 | 部署繁琐，Python wheel 需静态链接 | ❌ |
| 多共享库 (core/pricing/risk) | 模块化 | 版本兼容性地狱，加载慢 | ❌ |
| 纯头文件 (含 Python 绑定) | 无编译 | pybind11/nanobind 需编译 | ❌ |

---

## ADR-002: Bridge + Virtual Constructor (PayOff)

**状态**: Accepted  
**日期**: 2026-07-29  
**来源**: Joshi "C++ Design Patterns" Ch.4  
**关联 Phase**: 1

### 背景
`PayOff` 为抽象基类，含纯虚 `operator()(double)`。`VanillaOption` 需存储任意派生 PayOff，支持拷贝/赋值/销毁，且不切片。

### 决策
实现 **Bridge 模式** + **虚拟构造函数**:
```cpp
class PayOff {
public:
    virtual ~PayOff() = default;
    virtual double operator()(double spot) const = 0;
    virtual std::unique_ptr<PayOff> clone() const = 0;  // Virtual Constructor
};

class PayOffBridge {
    std::unique_ptr<PayOff> ptr_;
public:
    PayOffBridge() = default;
    PayOffBridge(const PayOff& p) : ptr_(p.clone()) {}
    PayOffBridge(const PayOffBridge& o) : ptr_(o.ptr_->clone()) {}
    PayOffBridge(PayOffBridge&&) = default;
    PayOffBridge& operator=(const PayOffBridge& o) { 
        if (this != &o) ptr_ = o.ptr_->clone(); return *this; 
    }
    double operator()(double spot) const { return (*ptr_)(spot); }
};
```
派生类实现 `clone()` 返回 `std::make_unique<Derived>(*this)`。

### 理由
1. **值语义**: `VanillaOption` 含 `PayOffBridge` 成员，表现如值类型
2. **开闭原则**: 新增 PayOff 无需修改 `VanillaOption` 或 `PayOffBridge`
3. **无切片**: `clone()` 保证派生类型完整复制
4. **Rule of 5**: Bridge 显式实现拷贝/移动/赋值/析构

### 后果
- 每次拷贝分配堆内存 (`new`) → Phase 4 引入对象池优化
- `unique_ptr` 单一所有权 → 共享所有权需改 `shared_ptr` (当前不需要)
- **路径相关期权限制**: `PayOff::operator()(double spot)` 仅适用于终值期权（Call/Put/Digital）。Asian/Barrier/Lookback 等路径相关期权需使用独立的 `PathDependentPayOff` 接口（参考 Joshi Ch.7-8），该接口接收完整路径并返回现金流向量。`VanillaOption` 持有 `PayOffBridge`，`ExoticOption` 持有 `PathDependentPayOffBridge`。详见 ARCHITECTURE.md §3.2。

---

## ADR-003: Factory + 静态注册模板

**状态**: Accepted  
**日期**: 2026-07-29  
**来源**: Joshi Ch.10  
**关联 Phase**: 1

### 背景
运行时根据字符串名称创建 PayOff/Model/Engine，新增类型零修改现有代码。

### 决策
**单例工厂 + 模板辅助类自动注册**:
```cpp
class PayOffFactory {
    using Creator = std::function<std::unique_ptr<PayOff>(const nlohmann::json&)>;
    std::unordered_map<std::string, Creator> registry_;
public:
    static PayOffFactory& instance() { static PayOffFactory f; return f; }
    void registerPayOff(const std::string& name, Creator c) { registry_[name] = std::move(c); }
    std::unique_ptr<PayOff> create(const std::string& name, const nlohmann::json& params) const;
    PayOffBridge createBridge(const std::string& name, const nlohmann::json& params) const;
};

template<class T>
struct PayOffRegistrar {
    PayOffRegistrar(const std::string& name) {
        PayOffFactory::instance().registerPayOff(name, [](const json& p) {
            return std::make_unique<T>(p);
        });
    }
};

// 使用: 每个 .cpp 文件静态实例化
static PayOffRegistrar<CallPayOff> reg_call("Call");
static PayOffRegistrar<PutPayOff> reg_put("Put");
```

### 理由
1. **零配置**: 只需链接包含注册代码的 `.cpp` 即可自动注册
2. **线程安全**: C++11 静态局部变量初始化保证
3. **参数化**: JSON 参数支持任意复杂构造 (DoubleDigital 需两行权价)
4. **可测试**: 单测可注册 Mock PayOff 验证工厂逻辑

### 后果
- 动态库加载顺序影响注册 → 显式链接所有注册单元
- 注册表无锁读 → 仅启动期写入，运行期只读

---

## ADR-004: 计数器 RNG + 确定性分块并行

**状态**: Accepted  
**日期**: 2026-07-29  
**来源**: pricer/rng.hpp, Random123, mape  
**关联 Phase**: 1

### 背景
Monte Carlo 需高质量、可并行、可复现的随机数。`std::mt19937` 有状态、难并行、种子管理复杂。

### 决策
采用 **Philox4x64 / Threefry4x64** (Random123) 计数器 RNG:
- **纯函数**: `draw(i) = f(key, counter=i)` 无内部状态
- **SIMD 友好**: 一次生成 4/8 个 64-bit 整数
- **并行安全**: 线程/块/进程只需不同 `key` 或 `counter` 偏移
- **确定性分块**: 固定 64 个块，每块独立种子 = `base_seed + block_id`

```cpp
struct Philox4x64 {
    using key_type = std::array<uint64_t, 2>;
    using ctr_type = std::array<uint64_t, 4>;
    key_type key_;
    explicit Philox4x64(uint64_t seed=0) : key_{seed, seed ^ 0x9E3779B97F4A7C15} {}
    ctr_type operator()(ctr_type ctr) const noexcept;  // 10 轮 Feistel
    // SIMD 版本: 一次生成 4 个 double
    void generate(double* out, size_t n, uint64_t base_ctr) const noexcept;
};
```

### 理由
1. **位精确复现**: 相同种子、相同路径数、不同线程数 → 位精确相同结果
2. **零同步开销**: 无锁、无原子操作、无线程局部存储竞争
3. **数学质量**: 通过 TestU01 BigCrush，周期 2^256
4. **可移植**: 纯 C++ 模板实现，无汇编，支持 CPU/GPU

### 后果
- 需维护 Random123 子模块版本锁定
- Box-Muller 正态变换需 SIMD 实现 (Phase 1 实现)

---

## ADR-005: 可移植 SIMD 抽象层

**状态**: Accepted  
**日期**: 2026-07-29  
**来源**: pricer/simd.hpp, simde, QuantKernel  
**关联 Phase**: 1

### 背景
AVX2/AVX-512/NEON 内在函数不兼容，手写多版本维护成本高。

### 决策
使用 **simde** (GitHub: simd-everywhere/simde) 作为可移植抽象层，仅在 dispatch 层写编译器特定代码:

```cpp
// core/simd.hpp
#if defined(__AVX2__)
    #include <immintrin.h>
    using f64x4 = __m256d;
    constexpr size_t SIMD_WIDTH = 4;
#elif defined(__AVX512F__)
    #include <immintrin.h>
    using f64x4 = __m512d;  // 实际 8 个 double
    constexpr size_t SIMD_WIDTH = 8;
#elif defined(__ARM_NEON)
    #include <arm_neon.h>
    using f64x4 = float64x2_t;  // NEON 128-bit = 2 double
    constexpr size_t SIMD_WIDTH = 2;
#else
    // 标量回退
    struct f64x4 { double v[4]; };
    constexpr size_t SIMD_WIDTH = 1;
#endif

// 统一接口
inline f64x4 simd_add(f64x4 a, f64x4 b);
inline f64x4 simd_mul(f64x4 a, f64x4 b);
inline f64x4 simd_exp(f64x4 x);  // 多项式近似
inline f64x4 simd_log(f64x4 x);
inline f64x4 simd_sqrt(f64x4 x);
inline double simd_hsum(f64x4 x);  // 水平求和
```

运行时检测 (`simd_dispatch.cpp`):
```cpp
enum class SimdLevel { Scalar, NEON, AVX2, AVX512 };
SimdLevel detect_simd_level() {
    #if defined(__AVX512F__) return SimdLevel::AVX512;
    #elif defined(__AVX2__) return SimdLevel::AVX2;
    #elif defined(__ARM_NEON) return SimdLevel::NEON;
    #else return SimdLevel::Scalar;
}
```

### 理由
1. **单一代码库**: 业务逻辑用统一 `simd_*` 接口
2. **零开销**: 内联展开为原生指令
3. **渐进增强**: 新指令集只需添加 dispatch 分支
4. **测试简单**: 标量回退作为参考实现验证正确性

---

## ADR-006: nanobind Python 绑定 + 批量 NumPy API

**状态**: Accepted  
**日期**: 2026-07-29  
**来源**: OptiCore, QuantKernel  
**关联 Phase**: 1

### 背景
Python 绑定需零拷贝、编译快、二进制小、支持 NumPy 批量计算。

### 决策
**nanobind** (而非 pybind11) + **批量 API 设计**:

```cpp
// 标量 API (便利性)
m.def("bsm_price", &bsm_price, "spot"_a, "strike"_a, "rate"_a, "vol"_a, "expiry"_a, "type"_a='C');

// 批量 API (性能) - 零拷贝 NumPy 指针直达
m.def("bsm_price_batch", [](nb::ndarray<nb::numpy, double> spots,
                              nb::ndarray<nb::numpy, double> strikes,
                              nb::ndarray<nb::numpy, double> rates,
                              nb::ndarray<nb::numpy, double> vols,
                              nb::ndarray<nb::numpy, double> expiries,
                              char type) {
    size_t n = spots.shape(0);
    auto out = nb::ndarray<nb::numpy, double>::create(n);
    bsm_price_batch_impl(spots.data(), strikes.data(), rates.data(),
                         vols.data(), expiries.data(), out.mutable_data(), n, type);
    return out;
}, "spots"_a, "strikes"_a, "rates"_a, "vols"_a, "expiries"_a, "type"_a='C');
```

**关键约束**:
- 所有批量函数 `noexcept`，输入数组要求 C-contiguous、aligned(64)
- 异常通过 `Result<T,ErrorCode>` 转换，绑定层统一 `throw` Python 异常
- GIL 释放: `nb::gil_scoped_release` 包裹耗时计算

### 理由
| 维度 | nanobind | pybind11 | ctypes/C API |
|------|----------|----------|--------------|
| 编译速度 | 4x 快 | 基准 | N/A |
| 二进制大小 | 5x 小 | 基准 | 最小 |
| NumPy 互操作 | 原生 ndarray | 需转换 | 手动 |
| C++20 支持 | 完整 | 部分 | 无 |

### 后果
- MSVC 需 `/std:c++20` 且 nanobind 2.3+
- 批量 API 文档需明确内存布局要求

---

## ADR-007: AAD (伴随自动微分) 统一 Greeks

**状态**: Accepted (Revised 2026-07-29)  
**日期**: 2026-07-29  
**来源**: autodiff 库 (https://autodiff.github.io), pricer/adjoint.hpp  
**关联 Phase**: 3

### 背景
Greeks 计算方法多 (解析/FD/Pathwise/LR/AAD)，需统一接口、自动选择最优方法、支持高阶。
原设计自研 `Tape` + `std::function<void()>` 反向传播，但 `std::function` 类型擦除阻止内联，性能比成熟库慢 10-50x，且内存估算偏低（实际 64-96 字节/节点 vs 声称 32 字节）。

### 决策
**采用 autodiff 库** (MIT 许可, header-only) 实现 AAD，不自研 Tape:

```cpp
#include <autodiff/forward/dual.hpp>
#include <autodiff/reverse/var.hpp>

using namespace autodiff;

// 前向模式: 单个 Greek
dual price_dual(const dual& S, const dual& K, const dual& r, 
                const dual& sigma, const dual& T);
// 自动获取 Delta/Vega 等

// 反向模式: 全 Greeks 一次扫描
var price_var(const var& S, const var& K, const var& r,
              const var& sigma, const var& T);
// price_var.backward() 后获取所有偏导
```

**统一 Greeks 接口** (不变):
```cpp
enum class GreeksMethod { Auto, Analytic, Pathwise, LR, FD, AAD };

struct GreeksEngine {
    Greeks compute(const Option& opt, const Model& model, GreeksMethod m = Auto) {
        if (m == Auto) {
            if (opt.is_european() && model.has_analytic_greeks()) return analytic_greeks(opt, model);
            if (opt.payoff_is_smooth() && model.supports_pathwise()) return pathwise_greeks(opt, model);
            if (!opt.payoff_is_smooth() && model.supports_lr()) return lr_greeks(opt, model);
            if (model.supports_aad()) return aad_greeks(opt, model);  // 内部用 autodiff
            return fd_greeks(opt, model);
        }
        // 显式分发...
    }
};
```

### 理由
1. **不自研**: autodiff 库经过充分测试，性能优化（表达式模板、内存池），自研 Tape 无法达到同等质量
2. **一次扫描全 Greeks**: 反向模式 O(1) 获取所有阶 Greeks
3. **精度**: 机器精度级，无数值微分误差
4. **统一性**: 复杂模型 (Heston/篮子/路径相关) 统一用 AAD
5. **可验证**: 其他方法作为回归测试基准

### 后果
- 新增依赖: autodiff (header-only, MIT), 通过 vcpkg 或 submodule 引入
- 编译时间略增: 模板代码膨胀，但远小于自研 Tape 的维护成本
- 与 nanobind 集成: autodiff 的 `dual`/`var` 类型不跨 Python 边界，绑定层用 `double` 转换

---

## ADR-008: COS/FFT 谱方法作为解析引擎补充

**状态**: Accepted  
**日期**: 2026-07-29  
**来源**: Fang & Oosterlee (2009), QuantKernel, hpc-stochastic-pricing-kernel  
**关联 Phase**: 2

### 背景
Heston/Levy 模型有闭式特征函数但无简单解析解。MC 收敛慢，PDE 高维困难。

### 决策
实现 **COS 方法** (Fang-Oosterlee) 作为主力谱引擎，FFT/CONV 作为备选:

```cpp
template<typename CharFn>
class COSEngine : public PricingEngine {
    CharFn char_fn_;
    size_t N_ = 256;  // 级数项数
    double L_ = 12.0; // 积分截断 [c1 - L*sqrt(c2), c1 + L*sqrt(c2)]
    
    double price(const VanillaOption& opt) const override {
        // 1. 计算 cumulants c1, c2, c4 (特征函数导数)
        // 2. 积分区间 [a,b] = [c1 - L*sqrt(c2), c1 + L*sqrt(c2)]
        // 3. 计算 payoff 系数 V_k (Call/Put 解析积分)
        // 4. 级数求和: sum_{k=0}^{N-1} ' Re[ phi(k*pi/(b-a)) * exp(-i*k*pi*(x-a)/(b-a)) ] * V_k
        // 5. 折现: exp(-r*T) * result
    }
};
```

**支持模型**: 任何提供 `characteristic_function(u, tau)` 的模型 (Heston, Bates, VG, CGMY, BS, SABR 近似)

### 理由
| 方法 | 收敛阶 | 适用模型 | 实现复杂度 |
|------|--------|----------|------------|
| COS | 指数 | 所有仿射/Levy | 低 |
| FFT | 指数 | 所有仿射 | 中 (需阻尼/插值) |
| CONV | 指数 | 波动率曲面校准 | 高 |

COS: N=256 达机器精度，单期权 < 1ms，无插值误差。

### 后果
- 需实现 Kahl-Jäckel 旋转轮廓避免分支切割 (复数 log/sqrt)
- cumulants 数值稳定计算 (高阶导数用复步微分)

---

## ADR-009: C ABI 稳定边界 + 版本化

**状态**: Accepted (2026-07-31, Phase 4 LITE 已实施)
**日期**: 2026-07-29 (Proposed) / 2026-07-31 (Accepted)
**关联 Phase**: 4 LITE

### 背景
C++ ABI 不稳定，需为多语言绑定 (Rust, C#, Julia, Excel XLL) 提供稳定接口。

### 决策
导出纯 C 接口 (`extern "C"`), 版本化符号 `cpphub_v1_*`:
```c
// include/cpphub/c_api/cpphub_c_api.h
#ifdef __cplusplus
extern "C" {
#endif

#define CPPHUB_ABI_VERSION 1

typedef struct { double delta, gamma, vega, theta, rho; } cpphub_greeks_t;
typedef struct { double price; double std_err; } cpphub_mc_result_t;

int cpphub_v1_bsm_price_batch(const double* spots, const double* strikes,
                              const double* rates, const double* vols,
                              const double* expiries, double* prices,
                              size_t n, char opt_type);

int cpphub_v1_heston_price(double S, double K, double T, double r, double q,
                           double v0, double kappa, double theta,
                           double sigma_v, double rho, double* price);

int cpphub_v1_bsm_implied_vol(double C_market, double S, double K, double T,
                              double r, double q, int is_call, double* iv);

int cpphub_v1_mc_price(double S, double K, double T, double r, double q,
                       double sigma, int is_call, size_t n_paths,
                       unsigned long long seed, cpphub_mc_result_t* result);

const char* cpphub_v1_get_last_error();
int cpphub_v1_get_abi_version();

#ifdef __cplusplus
}
#endif
```
实现文件 `src/c_api/c_api.cpp` 仅依赖 C++ 核心内部符号，不暴露任何 C++ 类型。

### 实施细节 (Phase 4 LITE)
- **构建**: `src/c_api/CMakeLists.txt` 编译为静态库 `cpphub_c_api`, 条件编译 `CPPHUB_ENABLE_C_API=ON` (默认 OFF, 不影响 header-only 核心)
- **错误处理**: 所有 C++ 异常在 C 边界捕获, 转为错误码 (0=success, -1=invalid_argument, -2=runtime_error, -3=unknown)
- **线程安全**: `cpphub_v1_get_last_error()` 使用 `thread_local` 存储, 每个 OS 线程独立错误消息
- **已实现函数**:
  - `cpphub_v1_bsm_price_batch`: 批量 BSM 定价 (委托 `bsm_call_price`/`bsm_put_price`)
  - `cpphub_v1_heston_price`: Heston call 定价 (委托 `detail::heston_call_price_cf`, Carr-Madan 傅里叶反演)
  - `cpphub_v1_bsm_implied_vol`: Newton-Raphson IV 反推 (委托 `bsm_implied_vol`)
  - `cpphub_v1_mc_price`: GBM 蒙特卡洛定价 (Philox RNG + Box-Muller)
  - `cpphub_v1_get_abi_version`: 返回 `CPPHUB_ABI_VERSION` (=1)
  - `cpphub_v1_get_last_error`: 返回 `thread_local` 错误字符串
- **未实现 (留待 v2.0+)**: `cpphub_v1_mc_price` 的 JSON 配置版本 (需 JSON 解析依赖), Excel XLL 直接调用接口 (ADR-012)

### 理由
1. **语言无关**: 任何支持 C FFI 的语言可直接调用 (Rust `extern "C"`, C# `DllImport`, Python `ctypes`)
2. **ABI 稳定**: 仅 POD 类型、错误码、不变符号
3. **版本共存**: `v1` `v2` 同进程加载互不干扰
4. **Excel XLL**: 直接加载 DLL 调用 `cpphub_v1_*` (ADR-012)
5. **零开销**: 静态链接时 LTO 内联, 动态链接时仅一次指针间接

### 后果
- C++ 异常需在 C 边界捕获转错误码 (已实现 `try/catch` 包装)
- 生命周期管理: 当前仅无状态函数, 不需 opaque handle; v2.0+ 引入 `cpphub_context_t*` 时再加引用计数
- 版本兼容: `CPPHUB_ABI_VERSION` 编译期常量, 不兼容变更需升 v2

---

## ADR-010: GPU 后端可选编译 + 运行时回退

**状态**: Accepted (2026-07-31, Phase 4 LITE 已实施)
**日期**: 2026-07-29 (Proposed) / 2026-07-31 (Accepted)
**关联 Phase**: 4 LITE

### 背景
GPU 加速显著但非所有部署环境有 CUDA/ROCm。

### 提议决策
- **CMake 选项**: `CPPHUB_ENABLE_CUDA=ON/OFF` (默认 OFF)
- **编译时**: 仅启用时编译 `.cu`，链接 `cudart`
- **运行时**: 无 CUDA 时自动编译 CPU stub (`gpu_mc_cpu_stub.cpp`),15 个 GPU MC 测试仍可在 CPU 上运行
- **RNG 一致性**: GPU Philox4x64-10 与 CPU `cpphub::core::rng` 算法完全一致,同 seed+counter → 同 Z (位精确)

```cmake
option(CPPHUB_ENABLE_CUDA "Build CUDA GPU MC kernel (requires NVIDIA GPU + nvcc)" OFF)
if(CPPHUB_ENABLE_CUDA)
    enable_language(CUDA)
    # ... gpu_mc.cu compiled with nvcc
else()
    # ... gpu_mc_cpu_stub.cpp compiled as fallback
endif()
```

```cpp
enum class Backend { Auto, CPU, CUDA, HIP };

class GpuDispatcher {
    static bool cuda_available_;
public:
    static void initialize() { cuda_available_ = (cudaGetDeviceCount(nullptr) == cudaSuccess); }
    static bool is_available(Backend b) {
        if (b == Backend::CUDA) return cuda_available_;
        return true;  // CPU 总可用
    }
};
```

### 理由
1. **零成本抽象**: 未启用 CUDA 时无任何编译/链接依赖
2. **部署灵活**: 同一 wheel 包含 CPU/GPU 代码，运行时自动选择
3. **CI 友好**: 无 GPU 的 CI 跑 CPU 路径，有 GPU 跑对比测试

---

## ADR-011: 分布式计算: MPI Master-Worker + 确定性聚合

**状态**: Proposed  
**日期**: 2026-07-29  
**关联 Phase**: 4

### 背景
超大规模 MC (亿级路径) 单机内存/算力不足。

### 提议决策
- **MPI 通信**: 标准 MPI-3 非阻塞 `MPI_Isend/Irecv` + `MPI_Reduce`
- **任务分发**: Master 静态分配固定块 (64块/rank)，种子 = `base_seed + global_block_id`
- **聚合**: `MPI_Reduce(MPI_SUM)` 求和 `sum` 和 `sum_sq`，Master 计算最终价格/标准误
- **容错**: 定期 `MPI_Bcast` checkpoint (每 1000 路径)，Rank 失效时重分配其块

```cpp
class DistributedMCEngine {
    MPI_Comm comm_;
    int rank_, size_;
    MCConfig cfg_;
    
    double price(const Option& opt) {
        // 1. 广播参数 + 基础种子
        // 2. 计算本 rank 负责块范围 [rank*blocks_per_rank, (rank+1)*blocks_per_rank)
        // 3. 本地计算 local_sum, local_sum_sq
        // 4. MPI_Reduce 到 rank 0
        // 5. rank 0 计算最终结果并广播 (可选)
    }
};
```

### 理由
1. **位精确复现**: 固定分块 + 确定性种子 → 与单机完全一致
2. **线性加速**: 无通信瓶颈 (仅两次 Reduce)
3. **标准化**: MPI 可移植、成熟、HPC 环境原生支持
3. **最小依赖**: 仅需 MPI 库，无额外框架

---

## ADR-012: xlOil XLL 加载项 (异步 UDF + 缓存)

**状态**: Proposed  
**日期**: 2026-07-29  
**关联 Phase**: 4

### 背景
Excel 是量化团队主要工作台，需高性能、异步、类型安全的 UDF。

### 提议决策
使用 **xlOil** (现代 C++17 XLL 框架，无 COM/ATL):

```cpp
// xll_addin.cpp
#include <xloil/ExcelObj.h>
#include <xloil/FunctionRegistry.h>
using namespace xloil;

// 同步 UDF (快速解析解)
XLOIL_EXPORT int CPPHUB_BSM_PRICE(double spot, double strike, double rate, 
                                  double vol, double expiry, const char* type,
                                  double* result) {
    *result = cpphub::bsm_price(spot, strike, rate, vol, expiry, *type);
    return 0;
}

// 异步 UDF (慢速 MC) - xlOil 原生支持
XLOIL_EXPORT XLOPER12* CPPHUB_MC_PRICE_ASYNC(double spot, double strike, 
                                             double rate, double vol, double expiry,
                                             const char* type, int n_paths,
                                             XLOPER12* handle) {
    // 返回异步句柄，Excel 显示 #GETTING_DATA... 直到完成
    auto task = std::make_shared<AsyncTask>([=]{
        return cpphub::mc_price(spot, strike, rate, vol, expiry, *type, n_paths);
    });
    return AsyncTask::register_task(task, handle);
}

// 缓存 UDF (依赖追踪)
XLOIL_EXPORT int CPPHUB_CACHED_PRICE(double spot, double strike, ...) {
    static LRUCache<std::tuple<...>, double> cache(10000);
    auto key = std::make_tuple(spot, strike, ...);
    if (auto it = cache.find(key); it != cache.end()) return it->second;
    double val = cpphub::bsm_price(...);
    cache.insert(key, val);
    return val;
}
```

**类型映射**: `double`/`int`/`string`/`Range`/`Array` → C++ 原生类型，自动转换

### 理由
1. **零样板代码**: 模板元编程自动生成 XLL 注册
2. **异步原生**: Excel 12+ 异步 UDF 标准，无阻塞 UI
3. **类型安全**: 编译期检查参数类型/数量
4. **性能**: 直接指针访问 Excel 内存，无 Variant 开销
5. **现代 C++**: 无 ATL/COM/MFC，头文件库，单 DLL

### 后果
- 需 xlOil 子模块锁定版本 (当前 1.8.x)
- Windows 专用，跨平台需单独方案 (PyXLL/Office.js)

---

## ADR-013: 双层线性代数架构 (固定尺寸 + 动态尺寸)

**状态**: Accepted (2026-08-04, Phase 6 实施前提升)
**版本归属**: v1.5 (Phase 6 经典参数计量模块引入时实施)

### 背景

当前 `core/linalg.hpp` 采用编译期固定尺寸的 Eigen-lite 设计（`template<size_t R, size_t C> class Matrix`），仅覆盖定价模块所需的 Cholesky 分解和 Thomas 算法。

计量经济学模块（未来 v1.1+）的矩阵运算需求与定价完全不同：

| 维度 | 定价需求 | 计量需求 |
|---|---|---|
| 矩阵尺寸 | 小且固定（4×4 相关矩阵、PDE 网格 <1000） | 大且动态（回归 N×K，N=万级，K=百级） |
| 核心操作 | Cholesky、Thomas | OLS (X'X)⁻¹、SVD、QR、稀疏矩阵 |
| 内存布局 | 栈分配 | 堆分配 + 可能需分块 |
| 并行 | OpenMP SIMD | BLAS Level 3 |
| 稀疏性 | 不需要 | 面板数据需要 |

如果强行扩展固定尺寸 `Matrix<R,C>` 到动态尺寸，会破坏定价模块的栈分配和全内联优化优势。

详细调研见 `docs/research/ECONOMETRICS_LANDSCAPE.md` §10。

### 决策

采用双层 linalg 架构：

1. **保留 `core/linalg.hpp`** (Eigen-lite, 固定尺寸)
   - 命名空间: `cpphub::linalg`
   - 用于: 定价模块（相关矩阵、PDE 网格）
   - 特性: 编译期 `R×C`，栈分配，表达式模板，全内联
   - 依赖: 无（header-only，不引入 Eigen3）

2. **新增 `core/linalg_dynamic.hpp`** (封装 Eigen3, 动态尺寸)
   - 命名空间: `cpphub::linalg::dynamic`
   - 用于: 计量模块（回归、SVD、QR、因子分解）
   - 特性: 运行时 `N×K`，堆分配，BLAS 加速，稀疏矩阵
   - 依赖: Eigen3 (header-only，可选 MKL/OpenBLAS 后端)

```cpp
// core/linalg_dynamic.hpp (v1.1+ 新增)
namespace cpphub::linalg::dynamic {
    class MatrixXD { Eigen::MatrixXd data_; /* ... */ };
    class VectorXD { Eigen::VectorXd data_; /* ... */ };
    class SparseMatrix { Eigen::SparseMatrix<double> data_; /* ... */ };

    // 计量核心操作 (委托 Eigen3 + LAPACK)
    MatrixXD svd(const MatrixXD& A);
    MatrixXD qr(const MatrixXD& A);
    VectorXD lstsq(const MatrixXD& A, const VectorXD& b);
    MatrixXD cholesky_dynamic(const MatrixXD& A);
}
```

### 理由

1. **关注点分离**: 定价和计量的矩阵运算特征完全不同，强行统一会牺牲性能
2. **编译时间**: Eigen3 编译慢（头文件 ~1.5MB），定价模块不应承受此成本
3. **依赖最小化**: v1.0/v1.1 定价模块保持 header-only，不引入新依赖
4. **渐进式引入**: 计量模块在 v1.1+ 才需要，此时引入 Eigen3 不影响已稳定的定价内核
5. **性能不妥协**: 定价用栈分配 + SIMD，计量用堆分配 + BLAS，各自最优

### 替代方案评估

| 方案 | 优势 | 劣势 | 结论 |
|---|---|---|---|
| A. 扩展固定尺寸 Matrix 支持动态 | 单一接口 | 破坏栈分配优化，模板膨胀 | ❌ 拒绝 |
| B. 直接用 Eigen3 替换 Eigen-lite | 代码少 | 定价模块编译时间增加，依赖膨胀 | ❌ 拒绝 |
| C. 双层 linalg（本方案） | 各自最优 | 两个命名空间 | ✅ 采纳 |
| D. 自研动态矩阵 | 无外部依赖 | 重复造轮子，SVD/QR 实现复杂 | ❌ 拒绝 |

### 后果

- 两个 linalg 命名空间: `linalg::` (固定) 和 `linalg::dynamic::` (动态)
- v1.1+ 引入 Eigen3 依赖（通过 vcpkg 或 submodule）
- 计量模块（`cpphub/econometrics/`）依赖 `linalg_dynamic.hpp`
- 定价模块不依赖 Eigen3，保持 header-only
- 未来可选: v2.0+ 引入 Intel MKL 或 OpenBLAS 作为 Eigen3 后端

---

## ADR-014: 标定 (Calibration) vs 估计 (Estimation) 的分离

**状态**: Accepted (2026-07-29)
**版本**: v1.1+ 规划

### 背景

Cpp_Hub 的"标定"模块 (Phase 3 / v1.1, [PHASE3_SPEC §4](../phases/phase3/PHASE3_SPEC.md)) 是**衍生品标定**: 用期权市场价格反推模型参数 (如 Heston 参数), 目标是让模型价格匹配市场 IV 面。这是工程导向的, 共享 optimizer.hpp (LM/DE), 通常不计算标准误差。

但用户 Research OS 的"因子失效诊断"方向需要的是**统计估计**: 用历史时间序列拟合统计模型参数 (如 GARCH), 目标是推断。这是推断导向的, 必须计算标准误差 (Newey-West / sandwich / Hessian), 必须做模型检验 (LR / Wald / LM)。

两者性质不同, 不应混在同一模块。

### 决策

将"标定"和"估计"分离为两个独立模块:

```
cpphub/calibration/           # 衍生品标定 (工程, 匹配市场) - Phase 3 / v1.1
├── optimizer.hpp              # LM + DE + Nelder-Mead (共享)
├── calibrator.hpp             # HestonCalibrator / SABRCalibrator / SVICalibrator
└── objective.hpp              # VEGA 加权 / 价格加权

cpphub/econometrics/estimation/  # 统计估计 (推断, 拟合历史) - v1.2+
├── mle.hpp                       # MLE / QMLE
├── gmm.hpp                       # GMM
├── bootstrap.hpp                 # Block / Wild bootstrap
├── standard_errors.hpp           # Hessian / Sandwich / Newey-West
└── hypothesis_tests.hpp          # Wald / LR / LM 检验
```

**共享**: optimizer.hpp (LM/DE 可复用于 MLE 优化), linalg_dynamic.hpp (ADR-013)
**不共享**: 目标函数, 标准误差计算, 模型检验

### 理由

| 维度 | 标定 (calibration/) | 估计 (estimation/) |
|---|---|---|
| **目标** | 匹配当前市场 | 推断历史规律 |
| **数据** | 截面 (IV 面) | 时间序列 |
| **标准误差** | 通常不计算 | **必须计算** |
| **模型检验** | 残差分析 | LR/Wald/LM + AIC/BIC |
| **输出** | 参数点估计 | 参数估计 + 置信区间 + p 值 |
| **典型场景** | Heston 参数 → 匹配 IV | GARCH 参数 → 拟合收益 |
| **C++ 现状** | QuantLib 有 | **完全没有** |

### 方案评估

| 方案 | 优点 | 缺点 | 结论 |
|---|---|---|---|
| A. 合并到 calibration/ | 单一模块 | 混淆两种性质不同的操作 | ❌ 拒绝 |
| B. 标定在 calibration/, 估计在 econometrics/ (本方案) | 性质分离, 共享 optimizer | 两个目录 | ✅ 采纳 |
| C. 完全独立, 不共享 optimizer | 解耦最大化 | 重复实现 LM/DE | ❌ 拒绝 |

### 后果

- v1.1 完成衍生品标定 (calibration/), 不含统计估计
- v1.2 (或 v1.3) 新增统计估计 (econometrics/estimation/), 复用 optimizer.hpp
- 统计估计模块依赖 linalg_dynamic.hpp (ADR-013), 需要先引入 Eigen3
- 用户 Research OS 的"因子失效诊断"方向依赖此模块 (GARCH MLE + Newey-West + Bootstrap)
- 详见 [PHASE3_SPEC §5](../phases/phase3/PHASE3_SPEC.md) 和 [ECONOMETRICS_LANDSCAPE.md §11](../research/ECONOMETRICS_LANDSCAPE.md)

---

## ADR-015: 证伪统计量模块边界 (通用 vs 模块特定)

**状态**: Accepted (2026-08-12)
**版本归属**: v1.6 (Phase 7A)
**关联 Phase**: 7A
**调研依据**: [ADR015_FALSIFICATION_MODULE_BOUNDARY_RESEARCH.md](../research/ADR015_FALSIFICATION_MODULE_BOUNDARY_RESEARCH.md) v1.2 (经三轮审计排幻觉, 共修正 10 个幻觉点)
**执行规格**: [PHASE7A_FALSIFICATION_SPEC.md](../phases/phase7/PHASE7A_FALSIFICATION_SPEC.md)

### 背景

Phase 7A 规划为 v1.0-v1.5 全部已实现模块 (risk/pricing/hfecon/econometrics) 系统性补齐证伪统计量, 共 11 个新增头文件 + 172 个测试用例。

11 个头文件采用"5 通用 + 6 模块特定"混合归属, 引出 5 个待决策的边界问题:

1. **Eigen3 隔离边界**: 通用诊断放 `econometrics/inference/` (链接 `cpphub_econometrics`, 含 Eigen3), hfecon/risk/pricing 调用时需链接 `cpphub_econometrics`, 违反 ADR-013 (cpphub_core 不污染 Eigen3)
2. **weak_identification 归属**: spec 标注"GMM 专用", 但物理位置在"通用" inference/ 目录
3. **TestResult 通用基类**: 是否引入统一结果结构, 参考 R `htest` S3 类
4. **命名空间不对称**: risk/pricing 无子命名空间 (`cpphub::v1`), econometrics/hfecon 有 (`cpphub::v1::econometrics`/`cpphub::v1::hfecon`)
5. **OLS 重复实现**: 通用诊断的辅助回归若不依赖 Eigen3, 需用 std::vector 重新实现 OLS, 与 `econometrics/estimation/ols.hpp` (Eigen3) 存在重复

**代码库事实** (全部已验证, 详见调研报告 §2):
- CMake `cpphub_core` 不链接 Eigen3 (CMakeLists.txt L41), `cpphub_econometrics` 链接 `eigen3_interface` (L43)
- hfecon/risk/pricing 当前零依赖 econometrics (Grep 验证)
- `har_model.hpp::ols_estimate` (L186-307) 已证明 std::vector + Gauss-Jordan OLS 可行, 不依赖 Eigen3

### 决策

采用 **方案 B: 通用诊断不依赖 Eigen3**, 5 个决策点如下:

#### 决策点 1: 方案选择 — B (通用诊断不依赖 Eigen3)

通用诊断头文件 (`econometrics/inference/` 下) **仅依赖 `core/`**, 不 `#include "cpphub/core/linalg_dynamic.hpp"`, 接口参数不用 `linalg::dynamic::MatrixXD`:
- 纯序列检验 (JB/LB/Box-Pierce/DM): 直接用 `vector<Real>` 实现
- 回归检验 (BG/BP/White/MZ/CUSUM/Andrews): 用 `detail/ols_simple.hpp` 实现辅助回归 (std::vector + Gauss-Jordan)

#### 决策点 2: weak_identification 移到 estimation/

`weak_identification.hpp` 从 `inference/` 移到 `estimation/`, 因 Cragg-Donald 浓度矩阵 `G_T = (X̃'X̃)^{-1/2} X̃'Z̃ (Z̃'Z̃)^{-1} Z̃'X̃ (X̃'X̃)^{-1/2}` 需矩阵特征值分解 (Eigen3 SVD), GMM/IV 专用。`estimation/` 已链接 `cpphub_econometrics` (含 Eigen3)。

#### 决策点 3: TestResultBase 组合方式 (含复合诊断例外)

引入 `struct TestResultBase { Real statistic; Real p_value; std::string method_name; bool reject_null; }`, 各 Result 结构体通过**组合**方式复用。依据: R `htest` S3 类有统一字段结构 (statistic/parameters/p.value/estimate/null.value/alternative/method/data.name)。

**复合诊断例外**: 复合诊断结构 (如 `VolatilityDiagnosticsResult`/`HARDiagnosticsResult`/`HEAVYDiagnosticsResult`) 聚合多个子检验, 不对应单一 htest, **不组合 `TestResultBase`**, 而是通过子结构 (如 `LjungBoxResult.base`) 间接获得统一接口。

**判定准则**:
- 单一检验统计量 (JB/LB/BG/BP/White/IM/MZ/DM/CD/CUSUM) → 组合 `TestResultBase`
- 聚合多个子检验的复合诊断 (VolatilityDiagnostics/HARDiagnostics/HEAVYDiagnostics) → 无 base 字段

#### 决策点 4: 命名空间维持现状 (Phase 7A scope)

risk/pricing 维持现有命名空间 `cpphub::v1` (无子命名空间), 新增诊断头文件 (risk_diagnostics/greeks_consistency/pricing_diagnostics) 落在 `cpphub::v1`。

**理由**: 新增子命名空间的影响范围涉及 risk/ 和 pricing/ 下所有现有头文件, 超出 Phase 7A scope (Phase 7A 是"新增证伪统计量", 不应触发现有代码的大规模命名空间调整)。可通过 `using` 声明过渡, 应在独立 ADR 中决策。

#### 决策点 5: OLS 重复实现可接受

通用诊断用 `detail/ols_simple.hpp` (std::vector + Gauss-Jordan) 实现辅助回归, ~50-80 行重复代码, 换取 Eigen3 隔离。

**ols_simple 与 ols_estimate 的有意差异**:
- `ols_simple` **不自动添加常数列**: 由调用方决定 (BG/BP/White 辅助回归设计矩阵各不相同)
- `ols_simple` **不计算 adj_r_squared / llh**: 诊断辅助回归仅需 R² 用于 LM 统计量
- `ols_estimate` 的"先例"含义: 指 Gauss-Jordan + partial pivoting 的**实现模式**先例, 非 API 完全复刻

### 11 个头文件归属

| # | 文件 | 归属 | 依赖 Eigen3 | Wave | 理由 |
|---|------|------|------------|------|------|
| 1 | residual_diagnostics.hpp | econometrics/inference/ | 否 | 1 | JB/LB/BG/BP/White, 用 std::vector OLS |
| 2 | volatility_diagnostics.hpp | econometrics/inference/ | 否 | 1 | 标准化残差检验, 仅需 vector<Real> |
| 3 | specification_tests.hpp | econometrics/inference/ | 否 | 2 | IM/MZ/DM, MZ 用 std::vector OLS |
| 4 | structural_break.hpp | econometrics/inference/ | 否 | 3 | CUSUM/Andrews, 用 std::vector OLS |
| 5 | conduction_metrics.hpp | econometrics/inference/ | 否 | - | v2.0+ 预留 |
| 6 | **weak_identification.hpp** | **econometrics/estimation/** | **是** | 2 | Cragg-Donald 需特征值分解, GMM 专用 |
| 7 | risk_diagnostics.hpp | risk/var/ | 否 | 1 | DQ/Berkowitz/ES 后验 |
| 8 | greeks_consistency.hpp | risk/greeks/ | 否 | 3 | Greeks 跨方法一致性 |
| 9 | pricing_diagnostics.hpp | pricing/ | 否 | 3 | IV 拟合优度 |
| 10 | hfecon_diagnostics.hpp | hfecon/models/ | 否 | **2** | HAR 残差, 跨 Wave 依赖: Wave 1 (residual/volatility) + Wave 2 (specification_tests MZ) |
| 11 | jump_test_diagnostics.hpp | hfecon/tests/ | 否 | 3 | 跳跃检验多重修正 |

**新增公共基础设施** (`econometrics/inference/detail/`):
- `test_result_base.hpp`: TestResultBase 通用结果基结构
- `ols_simple.hpp`: 轻量级 OLS (std::vector + Gauss-Jordan, 不依赖 Eigen3)

### 归属判定准则

| 准则 | 通用层 (econometrics/inference/) | 模块特定 (各自模块下) |
|------|--------------------------------|---------------------|
| 1. 适用范围 | >=2 个模块会调用 | 仅 1 个模块使用 |
| 2. 输入依赖 | 只需残差/矩阵等通用量 | 需模块特定对象 |
| 3. Eigen3 依赖 | **不依赖 linalg_dynamic** (用 std::vector) | 可依赖 |
| 4. 数学基础 | 经典统计文献 | 领域文献 |
| 5. 复用潜力 | 未来新模块会复用 | 不太可能跨模块复用 |

**特例**: weak_identification 虽满足准则 1 (GMM + 未来 LIML), 但不满足准则 3 (需 Eigen3), 移到 estimation/。

### 理由

1. **ADR-013 兼容**: 通用诊断不引入 Eigen3, hfecon/risk/pricing 调用时无需链接 cpphub_econometrics
2. **数学可行**: 除 Cragg-Donald 外, 所有通用诊断可用 std::vector 实现 (调研报告 §3.1 已分析)
3. **工程先例**: `har_model.hpp::ols_estimate` (L186-307) 已证明 std::vector OLS 可行
4. **业界对照**: 类似 R `sandwich` 包的"模块化构建块"设计, 通用工具不依赖具体模型; R `htest` S3 类提供统一结果结构
5. **依赖方向清晰**: hfecon→通用诊断 是单向, 无循环风险
6. **Wave 归属修正**: hfecon_diagnostics 横跨 Wave 1 (residual/volatility) 和 Wave 2 (specification_tests), 归入 Wave 2 以满足依赖闭包

### 替代方案评估

| 方案 | 优势 | 劣势 | 结论 |
|---|---|---|---|
| A. 维持现状 (5 通用 + 6 模块特定, weak_identification 在 inference/) | 无迁移成本 | M1 (Eigen3 隔离) 未解决, 违反 ADR-013 | ❌ 拒绝 |
| **B. 通用诊断不依赖 Eigen3** (本方案) | ADR-013 兼容, 依赖方向清晰 | OLS 重复实现 ~50-80 行 | ✅ 采纳 |
| C. 新建 diagnostics/ 顶层目录 | 通用诊断独立于 econometrics | 引入新顶层目录, 破坏现有架构简洁性 | ❌ 拒绝 |

### 后果

- **新增依赖**: hfecon→econometrics (首次引入, hfecon_diagnostics 调用通用 volatility_diagnostics/specification_tests)
- **OLS 重复**: `detail/ols_simple.hpp` (std::vector) 与 `estimation/ols.hpp` (Eigen3) 存在 ~50-80 行重复, 工程上可接受
- **命名空间不对称**: risk/pricing 诊断头文件落在 `cpphub::v1`, econometrics/hfecon 诊断在各自子命名空间, 未来需独立 ADR 统一
- **Wave 2 扩展**: hfecon_diagnostics 从 Wave 1 移到 Wave 2, 因依赖 specification_tests 的 MincerZarnowitzResult
- **TestResultBase 例外**: 复合诊断 (VolatilityDiagnosticsResult 等) 不组合 base 字段, 通过子结构间接获得统一接口
- **v2.0+ 预留**: `conduction_metrics.hpp` 为信息论事前度量预留空文件 (见 INFORMATION_THEORY_METRICS_RESEARCH.md)

### 关联

- **调研报告**: [ADR015_FALSIFICATION_MODULE_BOUNDARY_RESEARCH.md](../research/ADR015_FALSIFICATION_MODULE_BOUNDARY_RESEARCH.md) v1.2
- **执行规格**: [PHASE7A_FALSIFICATION_SPEC.md](../phases/phase7/PHASE7A_FALSIFICATION_SPEC.md)
- **前置 ADR**: [ADR-013](#adr-013-双层线性代数架构-固定尺寸--动态尺寸) (双层 linalg, Eigen3 隔离边界)
- **关联 ADR**: [ADR-014](#adr-014-标定-calibration-vs-估计-estimation-的分离) (calibration vs estimation 分离, estimation/ 目录已存在)

---

## ADR-016: 金融时间序列实施边界 (18 项)

**状态**: Accepted
**日期**: 2026-08-15
**版本归属**: v1.6 (Phase 7B M1/M2)
**关联 Phase**: 7B
**完整文档**: [ADR-016_FINANCIAL_TIMESERIES_BOUNDARY.md](ADR-016_FINANCIAL_TIMESERIES_BOUNDARY.md)
**调研依据**: [FINANCIAL_TIMESERIES_RESEARCH.md](../research/FINANCIAL_TIMESERIES_RESEARCH.md) v3.2 §13-14

### 背景

v1.6 GARCH 族与单位根检验实施中存在 18 个边界决策点 (GARCH 7 项 + 单位根 11 项), 需在实施前明确方案。

**决策原则**: (1) 优先与 `arch`/`rugarch`/`urca` 开源库对齐; (2) 理论正确性优先于实现便利; (3) 预留扩展接口 (v1.7 多元 GARCH / 协整)。

### 决策

#### GARCH 族边界决策 (7 项)

| # | 决策点 | 决策 | 依据 |
|---|--------|------|------|
| G-ADR1 | 方差初始值策略 | **backcast** (EWMA λ=0.94) | `arch` 默认, 小样本稳健 (G1) |
| G-ADR2 | 参数约束方法 | **SLSQP** | 支持不等式约束; `arch` 用 SLSQP (G15/G22); v1.5 optimizer 需扩展 |
| G-ADR3 | QMLE 协方差估计 | **sandwich** H⁻¹SH⁻¹ | Bollerslev-Wooldridge 1992; 非正态下 H⁻¹ 低估 SE (G9) |
| G-ADR4 | 标准化残差 JB 检验临界值 | **Bootstrap** | 小样本下 hₜ 有估计误差, 渐近临界值过度拒绝; 复用 v1.5 Bootstrap (G11) |
| G-ADR5 | t 分布自由度估计 | **联合 QMLE** | Bollerslev 1987; 两步法效率低且偏差大 (G14) |
| G-ADR6 | 优化器起始点策略 | **多起始** | GARCH 似然非凸, 多起始避免局部最优 (G16) |
| G-ADR7 | GARCH-M 参数化 | **c·σₜ** | `arch` 默认 ('vol' 标准差形式); Engle-Lilien-Robbins 1987 原始参数化 (G18) |

#### 单位根检验边界决策 (11 项)

| # | 决策点 | 决策 | 依据 |
|---|--------|------|------|
| U-ADR1 | ADF lag 选择规则 | **Schwert** floor(12(T/100)^0.25) | `arch` 默认; Schwert 1989 (U1) |
| U-ADR2 | ADF 方程形式选择 | **自动** (趋势检验后选择) | `arch` 默认自动; 减少用户误选 (U2) |
| U-ADR3 | 临界值来源 | **MacKinnon 2010** response surface | `arch` 已用 2010; 4 系数, 精度最高 (U3/U13) |
| U-ADR4 | DF-GLS demean 临界值 | **arch 独立系数表** | ERS 原表仅含 trend; demean 需 `arch` 扩展 (U7) |
| U-ADR5 | DF-GLS c̄值 | **ERS 固定值** c̄₁=-7.0, c̄₂=-13.5 | ERS 1996; detrending 回归必须含 c̄ (U8) |
| U-ADR6 | KPSS H0 方向标注 | **平稳性** (H₀=平稳) | KPSS 1992; API 必须明确标注, 避免与 ADF 混淆 (U12) |
| U-ADR7 | KPSS 长期方差带宽 | **Andrews 自动** | KPSS 1992; `arch` 默认 (U11) |
| U-ADR8 | PP 检验带宽 | **NW automatic** | Phillips-Perron 1988; `arch` 默认; PP 用 Bartlett 核 (U5/U6) |
| U-ADR9 | 基准对照软件 | **Python arch** | `arch` 用 MacKinnon 2010 (最新); `urca` 可能用旧版 |
| U-ADR10 | 方差比检验变体 | **全部** (Z₁/Z₂ + Chow-Denning + Chen-Deo + debiased) | 覆盖所有变体; 多重检验复用 Phase 7A (U14-U17) |
| U-ADR11 | 结构断点检验 | **推迟 v1.7** (Zivot-Andrews) | v1.6 scope 已满; Phase 7A 已有 CUSUM/Andrews 通用框架 (U18) |

### 后果

- **C1 优化器扩展** (前置): v1.5 `optimizer.hpp` 需扩展 SLSQP (G-ADR2)
- **C2 Estimator 接口**: 已支持 `TimeSeriesData`, MLE/OLS 派生类 visit 逻辑需补齐
- **C4 命名空间**: 时序模型统一落 `cpphub::v1::timeseries` (见 ADR-017)
- **C6 数据载体**: `EconData` 已含 `TimeSeriesData`, 仅缺 `MultivariateTSData` (阶段 5 前置)
- **测试基准**: Python `arch` + R `rugarch`, 容差 1e-10 至 1e-12
- **幻觉点核查**: 60 个幻觉点 (G1-G23, U1-U22) 实施时逐点核查

### 关联

- **完整文档**: [ADR-016_FINANCIAL_TIMESERIES_BOUNDARY.md](ADR-016_FINANCIAL_TIMESERIES_BOUNDARY.md)
- **调研报告**: [FINANCIAL_TIMESERIES_RESEARCH.md](../research/FINANCIAL_TIMESERIES_RESEARCH.md) v3.2 §13-14
- **前置 ADR**: [ADR-013](#adr-013-双层线性代数架构-固定尺寸--动态尺寸), [ADR-014](#adr-014-标定-calibration-vs-估计-estimation-的分离)
- **关联 ADR**: [ADR-017](#adr-017-时序模块命名空间-cpphubv1timeseries) (时序模块命名空间, 约束 C4)

---

## ADR-017: 时序模块命名空间 (`cpphub::v1::timeseries`)

**状态**: Accepted
**日期**: 2026-08-15
**版本归属**: v1.6 (Phase 7B)
**关联 Phase**: 7B
**完整文档**: [ADR-017_TIMESERIES_NAMESPACE.md](ADR-017_TIMESERIES_NAMESPACE.md)
**调研依据**: [FINANCIAL_TIMESERIES_RESEARCH.md](../research/FINANCIAL_TIMESERIES_RESEARCH.md) v3.2 §22 (约束 C4)
**前置 ADR**: [ADR-015](#adr-015-证伪统计量模块边界-通用-vs-模块特定) (命名空间不对称问题, 决策点 4)

### 背景

Cpp_Hub v1.0-v1.5 命名空间不对称: risk/pricing 落 `cpphub::v1` (无子命名空间), econometrics/hfecon 落子命名空间。Phase 7B 新增金融时间序列模块 (GARCH/单位根/ARIMA/MIDAS/VAR/DCC/因果检验) 需明确命名空间归属。

### 决策

**时序模型统一落 `cpphub::v1::timeseries` 子命名空间**, 与 `econometrics`/`hfecon` 平级。

```cpp
namespace cpphub {
inline namespace v1 {
namespace timeseries {

namespace garch { class GarchModel; /* ... */ }           // v1.6 M1
namespace unit_root { class ADFTest; /* ... */ }          // v1.6 M2
namespace arima { class ARIMAModel; /* ... */ }           // v1.6 M3
namespace midas { class MIDASModel; /* ... */ }           // v1.6 M4
namespace var { class VARModel; /* ... */ }               // v1.7
namespace multivariate_vol { class DCCModel; /* ... */ }  // v1.7
namespace causality { class GrangerTest; /* ... */ }      // v1.7+

}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
```

**物理目录**: `include/cpphub/timeseries/{garch,unit_root,arima,midas,var,multivariate_vol,causality}/`

### 理由

1. **与 econometrics/hfecon 平级**: 时序模块作为独立方法论体系, 应享有同等待遇
2. **方法论独立性**: 金融时间序列与横截面计量有本质差异 (波动率是研究对象、时间维度是核心、动态系统建模、信息论度量)
3. **避免顶层污染**: 防止 `cpphub::v1` 符号爆炸 (GARCH 7+ 类 + 单位根 6+ 类 + VAR 5+ 类 + 因果 10+ 类), 避免 `VARModel` 与 risk 的 `VaR` 混淆
4. **渐进式引入**: 不影响现有 risk/pricing/econometrics/hfecon 代码, 现有时序功能 (har_model/volatility_diagnostics) 保持原位
5. **为 v1.7+ 预留空间**: VAR/DCC/协整/长记忆/非线性/MS-AR 等未来扩展有清晰归属

### 替代方案评估

| 方案 | 优点 | 缺点 | 结论 |
|------|------|------|------|
| A. 落 `cpphub::v1` (如 risk/pricing) | 与 risk/pricing 一致 | 顶层符号爆炸, 与 VaR 混淆 | ❌ 拒绝 |
| **B. 落 `cpphub::v1::timeseries`** | 与 econometrics/hfecon 平级 | 新增子命名空间 | ✅ 采纳 |
| C. 落 `cpphub::v1::econometrics::timeseries` | 复用 econometrics | 时序 ≠ 计量, 层级过深 | ❌ 拒绝 |
| D. 统一所有模块到子命名空间 | 完全对称 | 大规模迁移, 超出 v1.6 scope | ❌ 拒绝 (留待 v2.0) |

### 后果

- **新模块归属**: v1.6+ 所有时序头文件落 `include/cpphub/timeseries/`, 命名空间 `cpphub::v1::timeseries`
- **现有模块不迁移**: hfecon 的 `har_model.hpp`、econometrics 的 `volatility_diagnostics.hpp` 保持原位
- **CMake 目标**: 新增 `cpphub_timeseries` INTERFACE 库, 或合并到 `cpphub_econometrics` (若依赖 Eigen3)
- **Eigen3 隔离**: 按 ADR-013 原则, 需 Eigen3 的时序方法 (GARCH 估计/VAR) 放 `cpphub_econometrics`, 纯序列检验 (方差比) 放 `cpphub_core`
- **命名空间不对称**: risk/pricing 无子命名空间问题留待 v2.0 独立 ADR (ADR-018) 统一决策
- **C ABI 版本化**: 时序模块用 `cpphub_v1_6_*` / `cpphub_v1_7_*` 版本前缀 (约束 C5)

### 关联

- **完整文档**: [ADR-017_TIMESERIES_NAMESPACE.md](ADR-017_TIMESERIES_NAMESPACE.md)
- **调研报告**: [FINANCIAL_TIMESERIES_RESEARCH.md](../research/FINANCIAL_TIMESERIES_RESEARCH.md) v3.2 §22
- **前置 ADR**: [ADR-013](#adr-013-双层线性代数架构-固定尺寸--动态尺寸) (Eigen3 隔离), [ADR-015](#adr-015-证伪统计量模块边界-通用-vs-模块特定) (命名空间不对称, 决策点 4)
- **关联 ADR**: [ADR-016](#adr-016-金融时间序列实施边界-18-项) (金融时间序列实施边界, 18 项决策)
- **后续 ADR**: risk/pricing 子命名空间统一留待 v2.0 独立 ADR (编号顺延; ADR-018 已用于 SLSQP 优化器边界)

---

## ADR-018: SLSQP 优化器实现边界

**状态**: Accepted
**日期**: 2026-08-15
**版本归属**: v1.6 (Phase 7B 阶段 1 前置)
**关联 Phase**: 7B
**完整文档**: [ADR-018_SLSQP_BOUNDARY.md](ADR-018_SLSQP_BOUNDARY.md)
**调研依据**: [SLSQP_EXTENSION_RESEARCH.md](../research/SLSQP_EXTENSION_RESEARCH.md) v1.0
**前置 ADR**: [ADR-016](#adr-016-金融时间序列实施边界-18-项) (G-ADR2: GARCH 参数约束方法选 SLSQP)

### 背景

Phase 7B 约束 C1: 现有 optimizer 仅 LM/NelderMead/DE 无约束优化器, 无法满足 GARCH 参数约束 (ω>0, α≥0, β≥0, α+β<1)。ADR-016 G-ADR2 决策采用 SLSQP (对齐 arch 默认 solnp 路径的数值对照需求)。

### 决策 (7 项)

1. QP 求解: **active-set** (KKT 系统求解等式约束 QP, 与 scipy 一致)
2. BFGS 更新: **Damped** (曲率条件不满足时 damped 版本保证 B 正定, Nocedal-Wright Eq 18.15)
3. 约束表达: **c_i(x) ≥ 0** (与 scipy `f_ieqcons` / arch `a·x−b` 一致)
4. 数值差分: **中心差分** (O(h²), 步长自适应)
5. 边界处理: **QP 内直接处理** (避免 2n 冗余约束)
6. Eigen3: **不引入** (ADR-013 隔离; n<10 用 std::vector + Gauss-Jordan)
7. 等式约束: **支持** (未来 ARIMA/MIDAS 固定持久性等需求)

### 后果

- GARCH/EGARCH/GJR QMLE 的约束优化路径打通; ADR-019 决策 23 复用 (MIDAS 集中化 NLS 外层) 与 C5 (GM QMLE)
- 接口: `SLSQP` 类 + `ConstraintFn`/`ConstraintJacobianFn` 回调

### 关联

- **完整文档**: [ADR-018_SLSQP_BOUNDARY.md](ADR-018_SLSQP_BOUNDARY.md)
- **后续 ADR**: [ADR-019](#adr-019-v17-多变量时序与混频实施边界-263-项) (决策 23/C5 复用本优化器)

---

## ADR-019: v1.7 多变量时序与混频实施边界 (26+3 项)

**状态**: Accepted
**日期**: 2026-08-17
**版本归属**: v1.7 (Phase 7C)
**关联 Phase**: 7C
**完整文档**: [ADR-019_V17_TIMESERIES_BOUNDARY.md](ADR-019_V17_TIMESERIES_BOUNDARY.md)
**调研依据**: [PHASE7C_RESEARCH.md](../research/PHASE7C_RESEARCH.md) v1.1 (3 agent 126 条全量审计)
**复核记录**: [ADR019_REVIEW_PILOT.md](../research/ADR019_REVIEW_PILOT.md) v1.1 (**R1-R4 调研门禁首轮 pilot, 全通过**)
**前置 ADR**: [ADR-016](#adr-016-金融时间序列实施边界-18-项) / [ADR-017](#adr-017-时序模块命名空间-cpphubv1timeseries) / [ADR-018](#adr-018-slsqp-优化器实现边界)

### 背景

v1.7 扩展多变量时序 (VAR/IRF/FEVD/DY/协整三件套) 与混频回归 (MIDAS), 回填三项 (Ng-Perron / Zivot-Andrews / GARCH-M)。首个运行于 DEVELOPMENT_WORKFLOW v1.1 阶段 0 (断言分级证据框架) 的调研→复核→冻结闭环: 调研报告经 126 条全量审计 (8 处实质错误修正), 决策集经机械探针 + 双盲重推导 (5/5 TRUE) + R2 抽检 (1 处引文失准当场修正), R1-R4 门禁全通过后冻结。

### 决策 (26+3 项, 摘要)

- **ARIMA/Granger (1-8)**: CSS+CSS-ML + innovations 精确 MLE; (1+θB) 正号参数化; d>0 无均值 + drift 选项; 多起始点; 不做 SARIMA/wild bootstrap
- **VAR/DY (9-16)**: 逐方程 OLS; IC 用 Lütkepohl 约定 + 同样本 offset; Cholesky 统一 Eigen LLT 下三角 + P 注入; FEVD 双轨 (Cholesky + generalized), DY 指数基于 GFEVD (自实现, 主基准 R Spillover); 不做 SVAR/BVAR/TVP-VAR
- **协整 (17-21)**: EG + Johansen 等价 API; EG 用 MacKinnon 1994 响应面; Johansen 临界值主录 Osterwald-Lenum 1992 + 双库 diff 前置冻结; PO 纳入 (Pz 优先); VECM β 默认 Phillips 归一
- **MIDAS (22-26)**: DL/U-MIDAS/MIDAS-AR 逐字复刻 midasr; 集中化 NLS (外层 SLSQP λ + 内层解析 δ); midasr 0.9 唯一主基准 (夹具固化 CSV); log-sum-exp 防溢出
- **回填 (NP/ZA/GM)**: NP 复用 ERS GLS 去势 + 原文公式钉死基准; ZA 三模型 + trim 参数化; GM 三变体双锚 (rugarch archpow ↔ arch 8.0 form) + fix() 互验三步法

### 后果

- **门禁制度实证**: R2 抓到复核报告自身 1 处形态 II 转述引文并当场修正 — "生成端引证, 审计端不采信"配置的首轮拦截收益; 双盲产出强于原链的证据 (B3: statsmodels `method` 参数实现层面未被使用)
- **正向**: 决策集全部阻断性断言双源, 可直接冻结进 PHASE7C_SPEC
- **负向**: MIDAS 夹具固化成本; Johansen 双库 diff 前置工作; GFEVD 自实现维护
- **开放问题**: SARIMA / SVAR 族 / ARDL-PSS / v1.8+ MIDAS 扩展 / 假设区 H1-H2 ([待定] 进 spec 开放问题节)

### 关联

- **完整文档**: [ADR-019_V17_TIMESERIES_BOUNDARY.md](ADR-019_V17_TIMESERIES_BOUNDARY.md)
- **调研报告**: [PHASE7C_RESEARCH.md](../research/PHASE7C_RESEARCH.md) v1.1
- **复核 pilot**: [ADR019_REVIEW_PILOT.md](../research/ADR019_REVIEW_PILOT.md) v1.1
- **门禁规范**: [ASSERTION_EVIDENCE_FRAMEWORK.md](../ASSERTION_EVIDENCE_FRAMEWORK.md) v1.1
- **后续**: PHASE7C_SPEC 编写 (R 清零 → spec 冻结 → 实施 → G 清零 → 合并)