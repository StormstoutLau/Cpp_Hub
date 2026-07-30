# Phase 1 执行规格书 - 核心内核 MVP

> **版本归属**: **v1.0 核心** (v1.0 = Phase 1 + Phase 2 削减后核心)
> **目标**: 3 周内交付可编译、可测试、可 Python 调用的最小可行定价内核
> **范围**: BS 欧式期权解析解 + MC 引擎(反变量+Sobol) + PayOff Bridge/Factory + 核心基础设施
> **里程碑**: M1(Week 1) 基础设施 → M2(Week 2) PayOff/MC/解析解 → M3(Week 3) Python 绑定/CI/文档

---

## 1. 交付物清单

### 1.1 必须编译通过的目标
| 目标 | 类型 | 入口点 |
|------|------|--------|
| `cpphub_core` | Interface Library | 仅头文件 |
| `cpphub` | Shared Library | `src/core/simd_dispatch.cpp`, `src/python/bindings.cpp` |
| `cpphub_python` | nanobind Module | `python/cpphub/__init__.py` |
| `cpphub_tests` | Executable | `tests/unit/*.cpp` |
| `cpphub_benchmarks` | Executable | `benchmarks/*.cpp` |

### 1.2 必须通过的测试
| 测试套件 | 用例数 | 覆盖模块 |
|----------|--------|----------|
| `test_core_math` | 15 | `math.hpp`, `constants.hpp` |
| `test_linalg` | 12 | `linalg.hpp` (Cholesky, Thomas) |
| `test_datetime` | 10 | `datetime.hpp` |
| `test_simd` | 8 | `simd.hpp` (AVX2/NEON 路径) |
| `test_rng` | 10 | `rng.hpp` (Philox 确定性、正态分布) |
| `test_payoff_bridge` | 12 | `payoff_bridge.hpp` (拷贝/移动/赋值/clone) |
| `test_payoff_factory` | 8 | `factory.hpp` (注册/创建/参数化) |
| `test_bsm_analytic` | 20 | `black_scholes.hpp` (价格/Greeks/边界) |
| `test_mc_engine` | 15 | `mc_engine.hpp` (反变量/控制变量/QMC/确定性) |
| `test_integration` | 5 | 端到端：Option → Engine → Price/Greeks |
| `test_python_bindings` | 12 | `bindings.cpp` (标量/批量/异常传播) |

### 1.3 必须达到的数值基准
| 基准 | 容差 | 参考值来源 |
|------|------|------------|
| BS Call/Put 价格 | 1e-12 | Haug (2007) Table 1-1 |
| BS Greeks (Delta/Gamma/Vega/Theta/Rho) | 1e-10 | 自动微分验证 |
| MC 欧式 Call (1M paths, antithetic) | 相对误差 < 0.5% | 解析解 |
| MC 标准误差收敛阶 | O(1/√N) | log-log 回归斜率 -0.5 ± 0.05 |
| Sobol QMC 方差缩减 | > 10x vs 伪随机 | 同路径数对比 |
| 确定性复现 (同种子/不同线程数) | 位精确相同 | `CHECK(price_1t == price_8t)` |

### 1.4 必须达到的性能基准 (GitHub Actions ubuntu-latest)
| 基准 | 相对目标 (推荐) | 绝对参考值 | 失败阈值 |
|------|------|------|----------|
| `BSM_Analytic_Batch_AVX2/10000` | > 5x `QuantLib::AnalyticEuropeanEngine` | ~30M prices/s | < 2x QL 或 < 10M/s |
| `MC_GBM_1M_Paths_Antithetic` | > 3x `QuantLib::MCEuropeanEngine` | ~25M paths/s | < 1.5x QL 或 < 8M/s |
| `MC_GBM_1M_Paths_Sobol` | > 3x `QuantLib::MCEuropeanEngine` (伪随机) | ~20M paths/s | < 1.5x QL 或 < 6M/s |
| `Cholesky_10x10` | > 5x NumPy `np.linalg.cholesky` | ~5M ops/s | < 2x NumPy 或 < 1M/s |

> **注意**: 绝对值依赖硬件 (CI runner CPU 型号不固定)，作为参考值而非硬性门禁。相对目标 (vs QuantLib/NumPy) 是主要门禁标准。绝对值失败时需人工复核硬件差异。

---

## 2. Week 1: 基础设施 (M1)

### 2.1 目录骨架创建 (Day 1)
```bash
# 执行脚本: scripts/bootstrap_phase1.sh
mkdir -p include/cpphub/{core,instruments/payoff,pricing/{analytic,monte_carlo},python}
mkdir -p src/{core,python}
mkdir -p tests/{unit,validation}
mkdir -p benchmarks
mkdir -p third_party/{nanobind,catch2,fmt,simde,random123}
mkdir -p python/cpphub
mkdir -p cmake/modules
```

### 2.2 依赖锁定 (Day 1)
| 依赖 | 版本 | 获取方式 | 说明 |
|------|------|----------|------|
| nanobind | 2.3.0 | git submodule | Python 绑定 |
| Catch2 | 3.5.3 | git submodule | 单元测试 |
| simde | 0.8.1 | git submodule | 可移植 SIMD |
| Random123 | 1.15 | git submodule | Philox/Threefry |
| fmt | 10.2.1 | vcpkg | 日志/格式化 (可选，用 std::format 替代) |

**统一依赖管理策略**: vcpkg 管理成熟库 (fmt)，submodule 管理需精确控制的头文件库 (nanobind/simde/Random123/Catch2)。
v1.0 **不使用 Conan**，避免双包管理器复杂度。vcpkg.json 仅包含 fmt。

### 2.3 CMake 配置 (Day 1-2)
```cmake
# CMakeLists.txt 关键片段
cmake_minimum_required(VERSION 3.25)
project(CppHub VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 编译优化
# 注意: 数值路径(MC/Greeks)禁用 -ffast-math, 因为它允许编译器重排浮点运算,
# 破坏结合律, 与"位精确复现"门禁(PHASE1_SPEC §1.3)矛盾.
# 仅批量纯函数(BS batch)可单独启用 -ffast-math, 但需从位精确门禁豁免.
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(
        $<$<CXX_COMPILER_ID:MSVC>:/O2 /arch:AVX2 /fp:precise /GL>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3 -march=x86-64-v3 -ffp-contract=off -flto>
    )
    # OpenMP SIMD only (no threading runtime)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/openmp:experimental>
                            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-fopenmp-simd>)
    endif()
endif()
# 批量纯函数目标可单独启用 fast-math:
# target_compile_options(cpphub_batch PRIVATE $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-ffast-math>)

# 头文件库
add_library(cpphub_core INTERFACE)
target_include_directories(cpphub_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(cpphub_core INTERFACE cxx_std_20)

# 共享库 (仅编译绑定+SIMD dispatch)
add_library(cpphub SHARED src/core/simd_dispatch.cpp src/python/bindings.cpp)
target_link_libraries(cpphub PUBLIC cpphub_core)
target_link_libraries(cpphub PRIVATE nanobind::nanobind ${OpenMP_CXX_FLAGS})

# Python 模块
find_package(Python3 REQUIRED COMPONENTS Development)
nanobind_add_module(cpphub_python src/python/bindings.cpp)
target_link_libraries(cpphub_python PRIVATE cpphub)

# 测试
enable_testing()
add_subdirectory(tests)
add_subdirectory(benchmarks)
```

### 2.4 核心基础设施实现 (Day 2-5)

| 文件 | 关键内容 | 行数 | 验收标准 |
|------|----------|------|----------|
| `core/config.hpp` | 特性宏、版本、编译器检测 | ~80 | `CPPHUB_HAS_AVX2` 在 MSVC/GCC/Clang 下正确 |
| `core/types.hpp` | `Real`, `Complex`, `Size`, `Index` | ~30 | 类型别名编译通过 |
| `core/constants.hpp` | `PI`, `SQRT_2`, `INV_SQRT_2PI`, `DBL_EPS` | ~40 | 数值与标准库一致 |
| `core/math.hpp` | `erf`, `erfc`, `cdf_normal`, `cdf_normal_inv`, `bessel_i` | ~120 | 对比 `std::erf`/`boost::math` 相对误差 < 1e-14 |
| `core/linalg.hpp` | `Matrix<R,C>`, `Vector<N>`, `cholesky`, `thomas_algorithm` | ~250 | Cholesky 分解对称正定矩阵；Thomas 算法三对角求解 |
| `core/datetime.hpp` | `Date`, `Calendar`(TARGET/NULL), `DayCount`(ACT/360, 30/360), `Schedule` | ~300 | ISDA 日期滚动规则、节假日处理 |
| `core/simd.hpp` | `f64x4`, `f64x8`, `exp`, `log`, `sqrt`, `hsum`, `blend`, `load/store` | ~300 | AVX2 路径生成正确汇编；NEON 回退工作 |
| `core/rng.hpp` | `Philox4x64`, `normal_simd`, `generate_correlated` (Cholesky) | ~200 | 确定性：相同种子相同输出；正态分布 KS 检验 p>0.05 |
| `core/parallel.hpp` | `ThreadPool`, `parallel_for`, `blocked_range` | ~150 | 固定块数并行、位精确复现 |
| `core/batch.hpp` | `bsm_price_batch`, `bsm_greeks_batch` | ~150 | SIMD 向量化、对齐内存访问 |
| `core/memory.hpp` | `aligned_allocator<64>`, `object_pool`, `stack_allocator` | ~100 | 64 字节对齐、无内存泄漏 |
| `core/error.hpp` | `Result<T,ErrorCode>`, `CppHubException`, 错误码枚举 | ~80 | 异常安全、错误码可转字符串 |

**每日必须**: `cmake --build build --target cpphub_tests && ctest --output-on-failure`

---

## 3. Week 2: PayOff/定价引擎/解析解 (M2)

### 3.1 PayOff 体系 (Day 8-10)

| 文件 | 关键类/函数 | 行数 | 溯源 |
|------|-------------|------|------|
| `instruments/payoff/payoff.hpp` | `PayOff` (纯虚 `operator()`, `clone()`, `name()`) | ~60 | Joshi Ch.3 |
| `instruments/payoff/payoff_bridge.hpp` | `PayOffBridge` (Bridge + Rule of 5 + `clone()`) | ~80 | Joshi Ch.4 p.72-75 |
| `instruments/payoff/vanilla.hpp` | `CallPayOff`, `PutPayOff`, `DigitalPayOff`, `DoubleDigitalPayOff` | ~100 | Joshi Ch.3, Ch.3.6 |
| `instruments/payoff/exotic.hpp` | `AsianPayOff`(算术/几何), `BarrierPayOff`(8类型), `LookbackPayOff` | ~150 | Duffy Ch.2.8 |
| `instruments/payoff/factory.hpp` | `PayOffFactory` (单例+map), `PayOffHelper<T>` 自动注册宏 | ~120 | Joshi Ch.10 p.175-184 |

**注册宏使用示例**:
```cpp
// instruments/payoff/vanilla.cpp
REGISTER_PAYOFF(Call, CallPayOff)
REGISTER_PAYOFF(Put, PutPayOff)
REGISTER_PAYOFF(Digital, DigitalPayOff)
REGISTER_PAYOFF(DoubleDigital, DoubleDigitalPayOff)
```

### 3.2 Monte Carlo 引擎 (Day 10-12)

| 文件 | 关键类/函数 | 行数 | 关键算法 |
|------|-------------|------|----------|
| `pricing/monte_carlo/mc_config.hpp` | `MCConfig` (n_paths, n_steps, antithetic, control_variate, qmc, seed, n_threads) | ~60 | - |
| `pricing/monte_carlo/statistics.hpp` | `StatisticsMC` (DumpOneResult, GetResultsSoFar, clone) | ~80 | Joshi Ch.5 |
| `pricing/monte_carlo/variance_reduction.hpp` | `AntitheticVR`, `ControlVariateVR` (Decorator 模式) | ~120 | Duffy Ch.2.5 |
| `pricing/monte_carlo/quasi_monte_carlo.hpp` | `SobolSequence` (Direction numbers, Gray code), `SobolPathGenerator` | ~200 | Duffy Ch.2.4 p.78-84 |
| `pricing/monte_carlo/path_generator.hpp` | `PathGenerator<Process>` (Euler/Milstein, Brownian Bridge) | ~150 | Duffy Ch.2.2, Ch.2.9 |
| `pricing/monte_carlo/mc_engine.hpp` | `MCEngine : PricingEngine` (price, greeks, batch_price) | ~200 | 确定性分块并行 |

**MCEngine::price 伪码 (确定性并行)**:
```cpp
double MCEngine::price(const VanillaOption& opt) const {
    const auto& process = opt.process();
    const auto& payoff = opt.payoff();
    const size_t n_blocks = 64;  // 固定块数
    const size_t paths_per_block = (cfg_.n_paths + n_blocks - 1) / n_blocks;
    
    double sum = 0, sum_sq = 0;
    #pragma omp parallel for schedule(static) reduction(+:sum,sum_sq)
    for (size_t b = 0; b < n_blocks; ++b) {
        Philox4x64 rng(cfg_.seed + b);
        double block_sum = 0, block_sq = 0;
        for (size_t p = 0; p < paths_per_block; ++p) {
            auto path = generate_path(process, opt.expiry(), cfg_.n_steps, rng);
            double payoff_val = payoff(path);
            if (cfg_.antithetic) {
                auto path_anti = generate_antithetic_path(...);
                payoff_val = 0.5 * (payoff_val + payoff(path_anti));
            }
            if (cfg_.control_variate) payoff_val = apply_control_variate(payoff_val, path);
            block_sum += payoff_val;
            block_sq += payoff_val * payoff_val;
        }
        sum += block_sum;
        sum_sq += block_sq;
    }
    double mean = sum / cfg_.n_paths;
    double std_err = std::sqrt((sum_sq / cfg_.n_paths - mean * mean) / cfg_.n_paths);
    return mean * std::exp(-opt.rate() * opt.expiry());  // 贴现
}
```

### 3.3 解析解引擎 (Day 12-14)

| 文件 | 关键函数 | 行数 | 验收 |
|------|----------|------|------|
| `pricing/analytic/black_scholes.hpp` | `bsm_price`, `bsm_greeks`, `bsm_implied_vol` (Newton+Bisect) | ~150 | 1e-12 vs Haug |
| `pricing/analytic/bachelier.hpp` | `bachelier_price`, `bachelier_greeks` (负利率) | ~80 | 1e-12 |
| `pricing/analytic/black76.hpp` | `black76_price` (期货期权) | ~60 | 1e-12 |
| `pricing/analytic/analytic_engine.hpp` | `AnalyticBSEngine : PricingEngine` | ~60 | 统一接口 |

**BSM Greeks 实现要求**:
```cpp
struct Greeks { double delta, gamma, vega, theta, rho; };

Greeks bsm_greeks(Real S, Real K, Real r, Real sigma, Real T, char type) {
    // 使用统一 d1/d2 计算，避免重复
    // Delta: type=='C' ? N(d1) : N(d1)-1
    // Gamma: N'(d1) / (S*sigma*sqrt(T))
    // Vega: S * N'(d1) * sqrt(T)
    // Theta: -(S*N'(d1)*sigma)/(2*sqrt(T)) - r*K*exp(-rT)*N(type=='C'?d2:-d2)
    // Rho: type=='C' ? K*T*exp(-rT)*N(d2) : -K*T*exp(-rT)*N(-d2)
}
```

### 3.4 期权/行权/工厂 (Day 14)

| 文件 | 关键类 | 行数 |
|------|--------|------|
| `instruments/option/vanilla_option.hpp` | `VanillaOption` (payoff_bridge, spot, rate, vol, expiry, process) | ~80 |
| `instruments/option/exercise.hpp` | `Exercise` (European, American, Bermudan 日期向量) | ~60 |
| `instruments/option/factory.hpp` | `OptionFactory::createEuropean/Vanilla/Barrier/Asian` | ~100 |

---

## 4. Week 3: Python 绑定 / CI / 文档 / 发布 (M3)

### 4.1 nanobind 绑定 (Day 15-17)

```cpp
// src/python/bindings.cpp
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/optional.h>

#include "cpphub/pricing/analytic/black_scholes.hpp"
#include "cpphub/pricing/monte_carlo/mc_engine.hpp"
#include "cpphub/instruments/payoff/factory.hpp"
#include "cpphub/instruments/option/vanilla_option.hpp"
#include "cpphub/core/batch.hpp"

namespace nb = nanobind;
using namespace cpphub;

NB_MODULE(cpphub, m) {
    m.doc() = "CppHub: High-performance quantitative finance library";

    // 标量 API
    m.def("bsm_price", &bsm_price, "spot"_a, "strike"_a, "rate"_a, "vol"_a, "expiry"_a, "type"_a='C');
    m.def("bsm_greeks", &bsm_greeks, "spot"_a, "strike"_a, "rate"_a, "vol"_a, "expiry"_a, "type"_a='C');
    m.def("bsm_implied_vol", &bsm_implied_vol, "price"_a, "spot"_a, "strike"_a, "rate"_a, "expiry"_a, "type"_a='C');

    // 批量 API (零拷贝 NumPy)
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

    m.def("bsm_greeks_batch", &bsm_greeks_batch_impl, ...);

    // PayOff
    nb::class_<PayOffBridge>(m, "PayOff")
        .def(nb::init<const std::string&, const nb::dict&>(), "name"_a, "params"_a)
        .def("__call__", &PayOffBridge::operator(), "spot"_a);

    // VanillaOption
    nb::class_<VanillaOption>(m, "VanillaOption")
        .def(nb::init<PayOffBridge, double, double, double, double>(),
             "payoff"_a, "spot"_a, "rate"_a, "vol"_a, "expiry"_a)
        .def_rw("spot", &VanillaOption::spot_)
        .def_rw("rate", &VanillaOption::rate_)
        .def_rw("vol", &VanillaOption::vol_);

    // Engines
    nb::class_<PricingEngine>(m, "PricingEngine")
        .def("price", &PricingEngine::price, "option"_a)
        .def("greeks", &PricingEngine::greeks, "option"_a);

    nb::class_<AnalyticBSEngine, PricingEngine>(m, "AnalyticBSEngine")
        .def(nb::init<>());

    nb::class_<MCEngine, PricingEngine>(m, "MCEngine")
        .def(nb::init<MCConfig>(), "config"_a = MCConfig{})
        .def("price", &MCEngine::price, "option"_a)
        .def("greeks", &MCEngine::greeks, "option"_a);

    // MCConfig
    nb::class_<MCConfig>(m, "MCConfig")
        .def(nb::init<>())
        .def_rw("n_paths", &MCConfig::n_paths)
        .def_rw("n_steps", &MCConfig::n_steps)
        .def_rw("antithetic", &MCConfig::antithetic)
        .def_rw("control_variate", &MCConfig::control_variate)
        .def_rw("quasi_random", &MCConfig::quasi_random)
        .def_rw("seed", &MCConfig::seed)
        .def_rw("n_threads", &MCConfig::n_threads);
}
```

### 4.2 Python 包装包 (Day 17)

```
python/
├── cpphub/
│   ├── __init__.py          # from ._cpphub import *; 版本号
│   ├── _version.py          # 生成自 CMake
│   ├── pricing.py           # 高级封装: price_european(), calibrate_heston()
│   ├── instruments.py       # Option(), PayOff() 便捷构造
│   └── utils.py             # date_utils, calendar_utils
├── pyproject.toml           # [build-system] = scikit-build-core + nanobind
├── setup.py                 # 最小配置
└── tests/
    ├── test_bindings.py     # pytest + hypothesis
    └── test_numerics.py     # 对比 QuantLib/解析解
```

**pyproject.toml 关键片段**:
```toml
[build-system]
requires = ["scikit-build-core>=0.10", "nanobind>=2.3", "wheel"]
build-backend = "scikit_build_core.build"

[project]
name = "cpphub"
version = "1.0.0"
description = "High-performance quantitative finance library"
requires-python = ">=3.10"
dependencies = ["numpy>=1.24", "pandas>=2.0"]

[tool.scikit-build]
cmake.minimum-version = "3.25"
cmake.args = ["-DCMAKE_BUILD_TYPE=Release"]
wheel.packages = ["python/cpphub"]
```

### 4.3 CI/CD 配置 (Day 18)

**.github/workflows/ci.yml**:
```yaml
name: CI
on: [push, pull_request]
jobs:
  build-test:
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, windows-latest, macos-latest]
        build_type: [Release, Debug]
        compiler: 
          - {os: ubuntu-latest, cc: gcc-13, cxx: g++-13}
          - {os: ubuntu-latest, cc: clang-17, cxx: clang++-17}
          - {os: windows-latest, cc: cl, cxx: cl}
          - {os: macos-latest, cc: clang, cxx: clang++}
    steps:
      - uses: actions/checkout@v4
        with: {submodules: recursive}
      - name: Install deps (Linux)
        if: runner.os == 'Linux'
        run: |
          sudo apt-get update && sudo apt-get install -y ${{ matrix.cc }} ${{ matrix.cxx }} libomp-dev
      - name: Configure
        run: cmake -B build -DCMAKE_BUILD_TYPE=${{ matrix.build_type }} -DCMAKE_C_COMPILER=${{ matrix.cc }} -DCMAKE_CXX_COMPILER=${{ matrix.cxx }}
      - name: Build
        run: cmake --build build --parallel
      - name: Test
        run: ctest --test-dir build --output-on-failure -j2
      - name: Benchmark (Release only)
        if: matrix.build_type == 'Release'
        run: ./build/benchmarks/cpphub_benchmarks --benchmark_format=json > bench.json
      - name: Upload benchmark
        uses: actions/upload-artifact@v4
        with: {name: bench-${{ matrix.os }}-${{ matrix.cc }}, path: bench.json}

  python-wheel:
    needs: build-test
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with: {python-version: '3.11'}
      - run: pip install scikit-build-core nanobind wheel
      - run: pip install -e ./python -v
      - run: pytest python/cpphub/tests/ -v
      - run: pip wheel ./python -w dist/
      - uses: actions/upload-artifact@v4
        with: {name: wheels, path: dist/}

  performance-guard:
    needs: build-test
    runs-on: ubuntu-latest
    if: github.event_name == 'pull_request'
    steps:
      - uses: actions/download-artifact@v4
        with: {name: bench-ubuntu-latest-gcc, path: bench}
      - run: |
          python scripts/perf_guard.py bench/bench.json --threshold 0.05
          # 如果任何基准回退 > 5% 则失败
```

**scripts/perf_guard.py**:
```python
import json, sys
def main():
    with open(sys.argv[1]) as f:
        curr = {b['name']: b['real_time'] for b in json.load(f)['benchmarks']}
    # 从 main 分支下载基准基线 (GitHub Pages 或 artifact)
    # 简化：内嵌基线
    baseline = {...}  # 从 main 分支 CI 固化
    for name, time in curr.items():
        if name in baseline:
            regress = (time - baseline[name]) / baseline[name]
            if regress > float(sys.argv[3].split('=')[1]):
                print(f"REGRESSION: {name} +{regress*100:.1f}%")
                sys.exit(1)
    print("All benchmarks within threshold")
if __name__ == '__main__': main()
```

### 4.4 文档生成 (Day 19)

```bash
# docs/Doxyfile 关键设置
PROJECT_NAME = "CppHub"
INPUT = ../include ../src
RECURSIVE = YES
GENERATE_HTML = YES
GENERATE_XML = YES
USE_MDFILE_AS_MAINPAGE = ../README.md
ENABLE_PREPROCESSING = YES
MACRO_EXPANSION = YES
EXPAND_ONLY_PREDEF = YES
PREDEFINED = CPPHUB_API=
```

**Sphinx Python 文档**:
```bash
# docs/conf.py
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.viewcode',
              'sphinx.ext.intersphinx', 'nbsphinx']
autodoc_mock_imports = ['numpy', 'pandas']
```

### 4.5 发布清单 (Day 20-21)

| 项目 | 命令 | 验证 |
|------|------|------|
| 版本标签 | `git tag -a v1.0.0 -m "Phase 1 MVP"` | `git describe --tags` |
| GitHub Release | `gh release create v1.0.0 --notes-file RELEASE.md` | 页面可访问 |
| PyPI 上传 | `twine upload dist/*` | `pip install cpphub==1.0.0` 成功 |
| Conda Forge | 提交 PR 到 `conda-forge/cpphub-feedstock` | `conda install -c conda-forge cpphub` |
| vcpkg Port | 提交 PR 到 `microsoft/vcpkg/ports/cpphub` | `vcpkg install cpphub` |
| Conan Center | 提交 PR 到 `conan-io/conan-center-index` | `conan install cpphub/1.0.0@` |

---

## 5. 验收检查表 (Phase 1 完成标志)

| 类别 | 检查项 | 通过标准 | 签名 |
|------|--------|----------|------|
| **编译** | 三平台三编译器零警告 | `cmake --build build --target cpphub --parallel` 无警告 |  |
| **单测** | 所有单测通过 | `ctest --output-on-failure` 全绿 |  |
| **数值** | 验收矩阵全达标 | `tests/validation/run_all.py` 生成报告 |  |
| **性能** | 基准达标且无回退 | `scripts/perf_guard.py` 通过 |  |
| **内存** | ASan/UBSan/TSan 0 报错 | CI 全绿 |  |
| **Python** | `pip install -e python` + `pytest` 全绿 | `import cpphub; cpphub.bsm_price_batch(...)` 工作 |  |
| **文档** | `make -C docs html` 无警告 | `docs/build/html/index.html` 存在 |  |
| **打包** | Wheel/Sdist/Conan/vcpkg 均可安装 | 干净环境 `pip install cpphub` 成功 |  |
| **溯源** | 每个 .hpp 头部有 SOURCE 标注 | `grep -r "SOURCE:" include/` 无遗漏 |  |

---

## 6. 风险与缓解

| 风险 | 影响 | 可能性 | 缓解措施 |
|------|------|--------|----------|
| SIMD 内在函数跨编译器差异 | 编译失败/错误数值 | 高 | 使用 `simde` 抽象层，仅在 dispatch 层写编译器特定代码 |
| nanobind 与 MSVC 兼容性 | Python 轮子构建失败 | 中 | CI 包含 MSVC + nanobind 测试，必要时回退 pybind11 |
| Philox 种子分块逻辑 Bug | 并行结果不确定性 | 中 | 单测 `test_rng_determinism` 强制 1/2/4/8 线程位精确对比 |
| Sobol 方向数生成错误 | QMC 收敛阶异常 | 低 | 对比 Duffy 书中 Table 2.1 前 10 个方向数 |
| 批量 API 内存对齐崩溃 | Segfault | 中 | `aligned_alloc(64)` + `assume_aligned`，单测覆盖非对齐输入 |
| 异常跨 C++/Python 边界崩溃 | 解释器崩溃 | 低 | 所有对外函数 `noexcept` + `Result<T,ErrorCode>`，绑定层统一转异常 |

---

**Phase 1 负责人**: _______________  
**审核人**: _______________  
**开始日期**: _______________  
**预计结束**: _______________