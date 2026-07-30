# Phase 4 执行规格书 - 高性能优化、生产化与发布

> **⚠️ 范围调整记录 (2026-07-31 评审)**:
> 原始 v2.0 scope (GPU + MPI + Excel + gRPC + 云原生 + 多平台发布) 为商业产品发布规格,
> 与单人研究 OS 的实际需求严重错配 (A/B 站无 CUDA 无法跨平台验证; Excel/gRPC/K8s 非研究必需;
> PhD 申请在即需聚焦研究产出)。经评审决定执行 **Phase 4 LITE (研究优先)** 范围:
>
> **Phase 4 LITE 实际交付** (1-2 周):
> - **M1**: Phase 3 审计整改 (E1/E2 性能基准 + D2 SSVI 跨期限 + G2 Python 交叉验证 + G4 基准索引)
> - **M2**: nanobind Python 绑定 (核心模块: BSM/Heston/MC/Greeks/VaR)
> - **M3**: GPU MC (主控站 RTX 4060 单独实现, CMake `option(CPPHUB_ENABLE_CUDA)` 条件编译, 无 CUDA 时 CPU 回退, 不纳入 A/B 站跨平台回归) + 数学附录文档
>
> **明确推迟到 v2.0+ 的项** (留待 PhD 入学后再评估):
> - MPI 分布式计算 (A/B 站无法验证)
> - Excel XLL 加载项 (非研究必需)
> - gRPC 微服务 / Prometheus / OpenTelemetry (生产化留待产品化阶段)
> - Docker / K8s / Conda Forge / vcpkg / Conan 多平台发布 (无发布需求)
> - Arrow/Parquet 持久化 (现有 DuckDB 已满足需求)
>
> **原始 v2.0 scope 保留在下方作为参考**, 实际执行以 LITE scope 为准。
>
> ---

> **版本归属**: **v2.0** (v2.0 = Phase 4 完整 scope，含 GPU/分布式/Excel/云原生)
> **目标**: 3 周内交付 GPU 后端、分布式计算、持久化、Excel/云集成、完整文档、v2.0 正式发布
> **前置**: Phase 1-3 (v1.0 + v1.1) 全部通过
> **里程碑**: M1(Week 1) GPU/分布式/持久化 → M2(Week 2) Excel/云/监控 → M3(Week 3) 文档/发布/验收
>
> **v2.0 Scope 声明** (2026-07-29 评审):
> - 本 Phase 为 v2.0，在 v1.0 + v1.1 基础上添加高性能计算与生产化能力
> - **v2.0 核心交付**:
>   - GPU 后端 (CUDA): MC 路径生成、PDE 求解器
>   - 分布式计算 (MPI/ZeroMQ): 多节点 MC
>   - 持久化 (Arrow/Parquet): 行情数据、模型序列化
>   - Excel XLL 加载项 (xlOil)
>   - 云原生: gRPC 微服务、Prometheus 监控、OpenTelemetry 追踪
>   - 多平台发布: Wheel + Conda Forge + vcpkg + Conan + Docker
> - **前置依赖**: v1.0 (Phase 1+2 核心) + v1.1 (Phase 3 风险/标定) 必须先完成
> - **发布版本号**: v2.0.0 (非 v1.0.0，因为 v1.0 已用于 Phase 1+2 的初始发布)

---

## 1. 交付物清单

### 1.1 新增编译目标
| 目标 | 类型 | 关键源文件 |
|------|------|------------|
| `cpphub_cuda` | Shared Library (可选) | `src/performance/gpu/*.cu`, `src/performance/gpu/*.cpp` |
| `cpphub_mpi` | Static Library (可选) | `src/performance/distributed/*.cpp` |
| `cpphub_excel` | XLL Add-in | `src/io/excel/*.cpp` (xlOil/xlwings) |
| `cpphub_python` (增量) | nanobind Module | `src/python/bindings_v2.cpp` (流式/异步 API) |

### 1.2 必须通过的新增测试
| 测试套件 | 用例数 | 覆盖模块 |
|----------|--------|----------|
| `test_gpu_mc` | 15 | `gpu_mc.cu` (路径生成/定价/显存管理) |
| `test_gpu_pde` | 10 | `gpu_pde.cu` (Thomas 算法/ADI GPU 版) |
| `test_mpi_mc` | 10 | `mpi_mc.cpp` (块分发/聚合/容错) |
| `test_persistence` | 15 | `arrow_io.cpp`, `parquet_io.cpp`, `model_serialization.cpp` |
| `test_excel_xll` | 20 | `xll_addin.cpp` (UDF/异步/缓存) |
| `test_streaming_api` | 10 | `async_engine.hpp` (Future/Callback/流式批量) |
| `test_monitoring` | 8 | `metrics.hpp`, `tracing.hpp`, `health_check.cpp` |
| `test_integration_phase4` | 12 | 端到端：GPU 加速标定→分布式风险→Excel 报表→监控告警 |

### 1.3 必须达到的性能基准
| 基准 | 目标 | 验收方式 |
|------|------|----------|
| GPU MC (Heston, 1M paths) | > 500M paths/s (RTX 4090) | `benchmarks/gpu_mc_benchmark.cu` |
| GPU PDE (400x2000 网格) | < 10ms (vs CPU 100ms) | `benchmarks/gpu_pde_benchmark.cu` |
| 分布式 MC (8 节点) | 线性加速比 > 0.9 | `benchmarks/mpi_scaling_benchmark.cpp` |
| 模型序列化/反序列化 | > 100 MB/s | `benchmarks/serialization_benchmark.cpp` |
| Excel UDF 调用开销 | < 5 μs/调用 | `benchmarks/excel_udf_benchmark.cpp` |
| 流式批量 API 吞吐 | > 1M prices/s 持续 | `benchmarks/streaming_benchmark.cpp` |

---

## 2. Week 1: GPU/分布式/持久化 (M1)

### 2.1 GPU 后端 (`include/cpphub/performance/gpu/`)

| 文件 | 关键类/方法 | 行数 | 核心内核 |
|------|-------------|------|----------|
| `gpu_config.hpp` | `GpuConfig` (device_id, stream_count, memory_pool_size) | ~80 | 设备管理 |
| `gpu_rng.cu` | `PhiloxGPU` (设备端计数器 RNG, 并行生成) | ~200 | `curand` 替代 |
| `gpu_mc.cu` | `GpuMCEngine` (路径生成/定价/归约 全 GPU) | ~300 | **核心**: 单线程单路径、共享内存归约 |
| `gpu_pde.cu` | `GpuPDEEngine` (Thomas 算法并行/ADI GPU 版) | ~250 | 循环展开、寄存器压力优化 |
| `gpu_memory.hpp` | `GpuMemoryPool` (cudaMallocAsync/内存池/统一内存) | ~150 | 零拷贝主机-设备 |
| `gpu_dispatch.hpp` | `GpuDispatcher` (CPU/GPU 自动回退/混合执行) | ~150 | 运行时决策 |

**GPU MC 内核设计规范**:
```cpp
// kernel: 每个线程处理 1 条路径 (或 warp 处理 1 条路径用于多因子)
__global__ void mc_kernel_gbm(
    const Real* __restrict__ spots,
    const Real* __restrict__ strikes,
    const Real* __restrict__ rates,
    const Real* __restrict__ vols,
    const Real* __restrict__ expiries,
    Real* __restrict__ prices,
    const Real* __restrict__ greeks,  // 可选输出
    Size n, char opt_type,
    PhiloxState rng_state) {
    Size idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    // Philox 设备端生成
    Philox4x64 rng(rng_state, idx);
    Real S = spots[idx];
    Real K = strikes[idx];
    Real r = rates[idx];
    Real v = vols[idx];
    Real T = expiries[idx];
    
    // 精确解一步到期 (欧式)
    Real drift = (r - 0.5*v*v) * T;
    Real vol_sqrt_t = v * sqrt(T);
    Real Z = rng.normal();  // Box-Muller 设备端
    Real ST = S * exp(drift + vol_sqrt_t * Z);
    
    Real payoff = (opt_type == 'C') ? fmaxf(ST - K, 0) : fmaxf(K - ST, 0);
    prices[idx] = payoff * exp(-r * T);
    
    // Greeks (可选): Pathwise 同路径计算
    if (greeks) { ... }
}

// Host 启动
void launch_gpu_mc(Size n, int block_size=256) {
    dim3 grid((n + block_size - 1) / block_size);
    mc_kernel_gbm<<<grid, block_size>>>(...);
    cudaDeviceSynchronize();
}
```

### 2.2 分布式计算 (`include/cpphub/performance/distributed/`)

| 文件 | 关键类/方法 | 行数 | 协议/模式 |
|------|-------------|------|----------|
| `mpi_context.hpp` | `MpiContext` (RAII 初始化/Finalize, 通信器管理) | ~100 | MPI-3 |
| `distributed_mc.hpp` | `DistributedMCEngine` (块分发/种子同步/结果聚合) | ~250 | Master-Worker |
| `distributed_pde.hpp` | `DistributedPDEEngine` (域分解/重叠通信/迭代同步) | ~200 | Domain Decomposition |
| `fault_tolerance.hpp` | `CheckpointManager` (周期性检查点/故障恢复/任务重分配) | ~200 | Checkpoint/Restart |
| `task_queue.hpp` | `DistributedTaskQueue` (Redis/RabbitMQ/ZeroMQ 后端) | ~150 | 异步任务分发 |

**分布式 MC 确定性聚合**:
```cpp
class DistributedMCEngine {
    MpiContext mpi_;
    MCConfig cfg_;
    
    double price(const Option& opt) {
        // 1. Master 广播参数 + 基础种子
        // 2. 每个 rank 负责固定块 (rank * blocks_per_rank ... (rank+1)*blocks_per_rank)
        //    种子 = base_seed + block_id (全局唯一)
        // 3. 本地计算 sum, sum_sq
        // 4. MPI_Reduce (MPI_SUM) 到 Master
        // 5. Master 计算最终价格
        // 关键: 静态块分配 + 确定性种子 → 位精确复现
    }
};
```

### 2.3 持久化与序列化 (`include/cpphub/io/serialization/`)

| 文件 | 关键类/方法 | 行数 | 格式/用途 |
|------|-------------|------|----------|
| `arrow_io.hpp` | `ArrowSerializer` (Schema 定义/RecordBatch 读写/零拷贝) | ~200 | 内存/网络传输 |
| `parquet_io.hpp` | `ParquetWriter/Reader` (分区/压缩/谓词下推) | ~200 | 磁盘存储/查询 |
| `model_serialization.hpp` | `ModelSerializer` (JSON/Protobuf/MessagePack 模型参数/校准结果) | ~200 | 模型版本管理 |
| `market_data_io.hpp` | `MarketDataStore` (时序/曲面/期权链/增量更新) | ~200 | 行情数据仓库 |

**模型序列化 Schema 示例 (Protobuf)**:
```protobuf
// model.proto
message HestonModel {
    string model_id = 1;
    string version = 2;
    int64 calibrated_at = 3;
    double v0 = 4;
    double kappa = 5;
    double theta = 6;
    double sigma = 7;
    double rho = 8;
    CalibrationInfo calib_info = 9;
}

message CalibrationInfo {
    double objective_value = 1;
    int iterations = 2;
    bool converged = 3;
    map<string, double> parameter_errors = 4;
    map<string, double> parameter_correlations = 5;
}
```

---

## 3. Week 2: Excel 集成/云原生/监控 (M2)

### 3.1 Excel XLL 加载项 (`include/cpphub/io/excel/`)

| 文件 | 关键类/方法 | 行数 | 技术栈 |
|------|-------------|------|--------|
| `xll_addin.cpp` | `XLL_Addin` (xlOil: C++17 头文件库, 无 COM) | ~300 | xlOil |
| `xll_udf.hpp` | `UDF_Registry` (函数注册/类型映射/异步/缓存) | ~200 | 模板元编程 |
| `xll_async.hpp` | `AsyncUDF` (Excel 12+ 异步回调/进度条/取消) | ~150 | XLL SDK |
| `xll_cache.hpp` | `UDFCache` (LRU/TTL/依赖追踪/失效通知) | ~150 | 缓存一致性 |
| `xll_types.hpp` | `ExcelTypeMapper` (Range/Array/Variant ↔ C++ 类型) | ~100 | 类型安全 |

**UDF 示例 (xlOil 风格)**:
```cpp
// 同步 UDF
XL_REGISTER_FUNCTION(CPPHUB_BSM_PRICE, "CPPHUB.BSM.PRICE",
    "Black-Scholes European option price",
    "Spot", "Strike", "Rate", "Vol", "Expiry", "Type",
    [](double S, double K, double r, double v, double T, const char* type) -> double {
        return cpphub::bsm_price(S, K, r, v, T, *type);
    });

// 异步 UDF (大批量)
XL_REGISTER_ASYNC_FUNCTION(CPPHUB_BSM_PRICE_BATCH, "CPPHUB.BSM.PRICE.BATCH",
    "Batch Black-Scholes pricing",
    "Spots", "Strikes", "Rates", "Vols", "Expiries", "Type",
    [](const xl::Array<double>& spots, const xl::Array<double>& strikes,
       const xl::Array<double>& rates, const xl::Array<double>& vols,
       const xl::Array<double>& expiries, const char* type) -> xl::AsyncResult<xl::Array<double>> {
        return cpphub::bsm_price_batch_async(spots, strikes, rates, vols, expiries, *type);
    });
```

### 3.2 云原生与流式 API (`include/cpphub/performance/` + `python/`)

| 文件 | 关键类/方法 | 行数 | 特性 |
|------|-------------|------|------|
| `async_engine.hpp` | `AsyncPricingEngine` (Future/Promise/回调/进度流) | ~200 | C++20 协程 |
| `streaming_batch.hpp` | `StreamingBatchProcessor` (背压/窗口/检查点/恰一次) | ~250 | 响应式流 |
| `service_mesh.hpp` | `GrpcService` (gRPC 定价/标定/风险微服务) | ~300 | Protobuf/gRPC |
| `metrics.hpp` | `MetricsCollector` (Prometheus 指标/Histogram/Summary/Gauge) | ~150 | 可观测性 |
| `tracing.hpp` | `DistributedTracing` (OpenTelemetry/W3C TraceContext) | ~150 | 分布式追踪 |

**Python 异步/流式 API**:
```python
# cpphub/async_pricing.py
import asyncio
from cpphub import AsyncPricingEngine, StreamBatch

async def price_stream(options: AsyncIterable[Option]) -> AsyncIterable[PriceResult]:
    engine = AsyncPricingEngine(mc_config)
    async for batch in StreamBatch(options, batch_size=10000):
        results = await engine.price_batch_async(batch)
        for r in results:
            yield r

# 使用
async def main():
    async for result in price_stream(option_generator()):
        print(f"{result.option_id}: {result.price:.6f} ± {result.std_err:.6f}")

# gRPC 服务定义 (pricing.proto)
service PricingService {
    rpc Price(PriceRequest) returns (PriceResponse);
    rpc PriceBatch(stream PriceRequest) returns (stream PriceResponse);
    rpc Calibrate(CalibrationRequest) returns (CalibrationResponse);
    rpc HealthCheck(HealthRequest) returns (HealthResponse);
}
```

### 3.3 监控与可观测性

```cpp
// metrics.hpp - Prometheus 兼容
class MetricsCollector {
    // Counter: 定价请求总数、错误数
    // Histogram: 定价延迟 (buckets: 1μs, 10μs, 100μs, 1ms, 10ms)
    // Gauge: 当前活跃计算任务数、GPU 显存使用
    // Summary: MC 标准误差分布、标定收敛步数
    
    static MetricsCollector& instance();
    void record_price_latency(double us);
    void record_mc_paths_per_sec(double pps);
    void record_calibration_iterations(int iter);
    void increment_errors(const char* error_type);
};

// tracing.hpp - OpenTelemetry
class Tracer {
    // 定价请求 Span: option_id, model, engine, path_count
    // 标定 Span: model, method, iterations, objective
    // 导出: OTLP/Jaeger/Zipkin
};
```

---

## 4. Week 3: 文档完善、发布工程、v1.0 验收 (M3)

### 4.1 文档体系完善

| 文档 | 位置 | 关键内容 | 受众 |
|------|------|----------|------|
| **用户指南** | `docs/user_guide/` | 快速开始/核心概念/常见工作流/性能调优 | 量化研究员 |
| **API 参考** | `docs/api/` | Doxygen 生成 C++ API / Sphinx 生成 Python API | 开发者 |
| **数学附录** | `docs/math/` | 所有模型推导/数值方法证明/参考文献 | 模型验证 |
| **部署指南** | `docs/deployment/` | Docker/K8s/系统依赖/GPU/分布式配置 | 运维/DevOps |
| **迁移指南** | `docs/migration/` | 从 QuantLib/自研库迁移/接口对照表 | 架构师 |
| **最佳实践** | `docs/best_practices/` | 数值稳定性/标定陷阱/性能分析/测试策略 | 全员 |

**文档构建管道**:
```yaml
# .github/workflows/docs.yml
- name: Build C++ API (Doxygen)
  run: doxygen docs/Doxyfile
- name: Build Python API (Sphinx)
  run: cd docs && make html SPHINXOPTS="-W"
- name: Build Math Appendix (LaTeX)
  run: cd docs/math && latexmk -pdf appendix.tex
- name: Deploy to GitHub Pages
  uses: peaceiris/actions-gh-pages@v3
  with: {github_token: ${{ secrets.GITHUB_TOKEN }}, publish_dir: ./docs/build/html}
```

### 4.2 发布工程清单

| 步骤 | 命令/动作 | 验证标准 | 负责人 |
|------|-----------|----------|--------|
| 版本确定 | `git tag -a v2.0.0 -m "Release v2.0.0"` | 语义化版本、CHANGELOG 更新 | Release Manager |
| 多平台 Wheel | `cibuildwheel --platform linux/macos/windows` | `pip install cpphub==2.0.0` 全平台通过 | CI/CD |
| Conda Forge | PR 到 `conda-forge/cpphub-feedstock` | `conda install -c conda-forge cpphub` | 维护者 |
| vcpkg Port | PR 到 `microsoft/vcpkg/ports/cpphub` | `vcpkg install cpphub:x64-windows` | 维护者 |
| Conan Center | PR 到 `conan-io/conan-center-index` | `conan install cpphub/2.0.0@` | 维护者 |
| Docker 镜像 | `docker build -t ghcr.io/org/cpphub:2.0.0 .` | `docker run --rm ghcr.io/org/cpphub:2.0.0 python -c "import cpphub"` | CI/CD |
| GitHub Release | `gh release create v2.0.0 --notes-file RELEASE.md --generate-notes` | 页面包含二进制、Wheel、源码、校验和 | Release Manager |
| 文档部署 | 合并 `docs` 分支触发 Pages 部署 | `https://org.github.io/cpphub/` 可访问 | CI/CD |
| 宣传通知 | 发送邮件/Slack/社区论坛 | 利益相关者确认收到 | PM |

### 4.3 v2.0 最终验收矩阵

| 维度 | 指标 | v2.0 门槛 | 实际值 | 通过 |
|------|------|-----------|--------|------|
| **功能完整性** | 核心模型覆盖 | BS/Heston/SABR/Bates/VG/IR/LocalVol | | ☐ |
| | 定价引擎 | Analytic/MC/PDE/Tree/Fourier/LSMC | | ☐ |
| | Greeks 体系 | Analytic/Pathwise/LR/FD/AAD 全覆盖 | | ☐ |
| | 风险度量 | VaR/ES/压力测试/场景分析 | | ☐ |
| | 标定框架 | LM/DE/NM + 多目标 + 诊断 | | ☐ |
| | 高性能 | GPU MC > 500M paths/s, 分布式线性加速 > 0.9 | | ☐ |
| | 生产化 | Excel XLL, gRPC, Prometheus, Docker, K8s | | ☐ |
| **数值正确性** | 解析解对比 | 1e-12 (BS), 1e-8 (Heston COS) | | ☐ |
| | MC 收敛阶 | O(1/√N) ± 0.05 | | ☐ |
| | PDE 收敛阶 | CN O(Δt², Δx²) | | ☐ |
| | 标定重现性 | 同种子同结果、目标函数 < 1e-6 | | ☐ |
| | CPU/GPU 一致性 | 双精度位精确或相对差 < 1e-12 | | ☐ |
| **性能** | 单欧式期权延迟 | < 1 μs (批量 < 50 ns) | | ☐ |
| | MC 吞吐 | > 50M paths/s (CPU), > 500M (GPU) | | ☐ |
| | 标定时间 | Heston < 10s, SABR < 1s | | ☐ |
| | 内存占用 | < 500MB (万期权组合全 Greeks) | | ☐ |
| **工程质量** | 编译警告 | 0 (三平台三编译器) | | ☐ |
| | 单测覆盖率 | > 90% (核心模块 > 95%) | | ☐ |
| | 静态分析 | Clang-Tidy 0 报错 | | ☐ |
| | Sanitizers | ASan/UBSan/TSan 0 报错 | | ☐ |
| | 文档完整性 | 所有公开 API 有文档、数学附录齐全 | | ☐ |
| **生产就绪** | Python Wheel | manylinux2014 / macOS arm64/x64 / Windows | | ☐ |
| | C++ ABI 稳定 | v2.0.x 补丁版本不破坏 ABI | | ☐ |
| | 监控指标 | Prometheus/Grafana 开箱即用 | | ☐ |
| | 部署文档 | Docker/K8s/裸金属/GPU 全覆盖 | | ☐ |

---

## 5. 风险与缓解 (Phase 4)

| 风险 | 影响 | 可能性 | 缓解措施 |
|------|------|--------|----------|
| GPU 内核数值差异 (单精度/双精度/排序) | 结果不一致 | 高 | 强制双精度 `Real`; 使用 Kahan 求和; 单测强制 CPU/GPU 位精确对比 |
| CUDA 版本兼容性 (11.x/12.x) | 构建失败 | 中 | CMake `find_package(CUDA)` 版本范围; 提供预编译 PTX/SASS |
| MPI 死锁/通信开销大 | 分布式不加速 | 中 | 非阻塞通信 `MPI_Isend/Irecv`; 重叠计算通信; 任务粒度自适应 |
| xlOil/XLL 加载崩溃 | Excel 不可用 | 中 | 独立进程隔离 (RTD/RTD Server); 详细错误日志; 兼容性测试矩阵 |
| gRPC/Protobuf 版本冲突 | 微服务通信失败 | 低 | 固定 Protobuf 版本; `protoc` 版本锁定; 契约测试 |
| 发布流程手动步骤遗漏 | 版本不一致/遗漏平台 | 中 | 完全自动化 `release.yml`; 必须通过所有检查才能发布 |

---

**Phase 4 负责人**: _______________  
**审核人**: _______________  
**开始日期**: _______________  
**预计结束**: _______________