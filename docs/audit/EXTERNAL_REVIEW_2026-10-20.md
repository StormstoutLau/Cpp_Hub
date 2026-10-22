# 外部代码评审核验与 Issue 清单

> **文档定位**: 对一份外部（第三方 LLM）代码评审意见的**逐条核验**结果与**可跟踪 Issue 清单**。
> **评审来源**: 外部评审快照（HIGH 4 / MEDIUM 5 / LOW 1 + 2 项安全移交），快照时点介于 S2（2026-10-18）与 S5（2026-10-20）之间——依据：评审可见 `cpphub_tool_registry.hpp`（1160 行，S2 定型）但称"CI 无 Python job"（S5 已修）。
> **核验方法**: 本机读码（E3，当前工作树）+ 一次源比对（E1，Random123 `philox.h`）+ 数值手算；行号抽查约 20 处。
> **核验结论**: **评审 12 条实质全部成立，无幻觉、无捏造**；1 条因本会话变更而过时（CI Python job）；评审未覆盖 1 条同族问题（本库自补）。
> **编号约定**: 延续 `docs/DEVELOPMENT_LOG.md` 的 RISK 序列（现存最大 RISK-022，2026-07-31），本清单登记 **RISK-023 ~ RISK-035**。

---

## 一、核验结论摘要

| 项 | 数量 | 说明 |
|---|---|---|
| 外部评审条目（含安全移交） | 12 | HIGH 4 + MEDIUM 5 + LOW 1 + 安全 2 |
| 裁定 **成立** | 11 | 逐条给出文件:行证据锚点（§三） |
| 裁定 **部分成立** | 1 | RISK-032：文档/版本子项成立；CI Python job 子项**已由本会话 S5 修复** |
| 裁定 **不成立/幻觉** | 0 | 行号引用绝大多数逐字命中，仅 3 处轻微偏移（±1~±2 行） |
| 本库**增补** | 2 | RISK-035（S5 冻结 golden 时发现）、RISK-036（本次修复时发现，均评审未覆盖） |
| **合计登记** | **14** | 见 §二 主表 |

**证据分级**（沿用 `ASSERTION_EVIDENCE_FRAMEWORK.md`）：

- **E1** 一次源逐字比对（Random123 `philox.h`；版本号/行数计数命令）
- **E3** 静态读码（绑定当前工作树 commit，含未提交改动）
- **手算** 数值量级验算（下溢阈值、卡方近似偏差）

---

## 二、Issue 主表

| 编号 | 模块 | 问题摘要 | 评审级 | 裁定 | 优先级 | 状态 |
|---|---|---|---|---|---|---|
| RISK-023 | risk/var | Christoffersen 联合检验只含覆盖率分量却按 χ²(2) 取值 | HIGH | 成立 | **P0** | **已修复**（2026-10-21，见 §三） |
| RISK-024 | risk/var | 回测似然连乘下溢 → 静默失效（NaN / LR=0） | HIGH | 成立 | **P0** | **已修复**（2026-10-21，见 §三） |
| RISK-025 | performance/gpu | PIMPL 裸指针 + 默认析构 → Impl 与 CUDA stream 永久泄漏 | HIGH | 成立 | **P0** | **已修复**（2026-10-21，见 §三） |
| RISK-026 | performance/gpu | CUDA 调用/kernel launch 零错误检查；block_size 未校验 | HIGH | 成立 | **P0** | **已修复**（2026-10-21，见 §三） |
| RISK-027 | core/rng | "Philox4x64-10" 名实不符；无 KAT；GPU/CPU 流不一致 | MEDIUM | 成立 | P1 | 调研完成，待裁决（[spec/p1-fixes/RESEARCH.md §3.1/§6.4 U-1](../../spec/p1-fixes/RESEARCH.md)） |
| RISK-028 | risk/var | χ² CDF 用 Wilson-Hilferty 近似，k=1 临界带内翻转 | MEDIUM | 成立 | **P0** | **已修复**（2026-10-21，见 §三） |
| RISK-029 | 构建 | Release 写死 `-march=x86-64-v3` / `/arch:AVX2` | MEDIUM | 成立 | P2 | 待处理 |
| RISK-030 | pricing/pde | const 接口写 `mutable` 计数器 → 共享引擎数据竞争 | MEDIUM | 成立 | P2 | 待处理 |
| RISK-031 | pricing/pde | 边界按 `payoff.name()` 字符串分派；Neumann 配置被静默忽略 | MEDIUM | 成立 | P2 | 待处理 |
| RISK-032 | docs/版本/CI | README GPU 测试数不实；版本号三处不一致；无覆盖率 | LOW | **部分成立** | P1 | 调研完成，待裁决（RESEARCH §3.2 / U-2、U-3）；CI 子项已关闭（S5） |
| RISK-033 | 构建/供应链 | GoogleTest `FetchContent` 无 `URL_HASH` | 安全 | 成立 | P1 | 调研完成（RESEARCH §3.3，已取得 commit SHA 与归档 SHA256） |
| RISK-034 | c_api | `json_min` 解析外部输入（总量无上限）→ 移交安全复核 | 安全 | 成立（待复核） | P1 | 调研完成（RESEARCH §3.4 / U-4），加固待裁决 |
| RISK-035 | c_api（MCP） | `chow_denning_joint` 把联合 p 值填进 `statistic` 字段 | 自补 | 成立 | P2 | 待裁定 |
| RISK-036 | risk/var | χ² 上尾按 `1-P` 计算 → p ≲ 1e-16 时归零（判定不受影响） | 自补 | 成立 | P3 | 已登记（待裁定：是否改精确上尾） |

---

## 三、逐条证据与处置

### RISK-023 — Christoffersen 联合检验分量缺失（HIGH）

**事实**（`include/cpphub/risk/var/backtesting.hpp:122-158`）：

- `pof_result`（:133）计算后从未使用；
- `ind_LR`（:148-151）计算后于 :153 被丢弃；
- `LR_cc_total = LR_cc`（:153）**只含覆盖率分量**，却在 :155 按 `chi2_cdf(..., 2)` 取值。

**核验补充（对评审建议的修正）**：被丢弃的 `ind_LR` 表达式**本身也不是合法的 LR_ind**——它用无条件违约率 `ind_result.violation_rate` 冒充转移概率。正确来源是 `ChristoffersenIID::test` 内 :107-114 已算好的、基于 n00/n01/n10/n11 转移表的 `LR_ind`，应将其暴露后相加：

```
LR_cc = LR_pof + LR_ind  ~ χ²(2)
```

**现有用例为何没抓到**：`tests/unit/risk/var/test_backtesting.cpp:72-83` 的 `JointTest` 仅断言 `ind.p_value <= joint.p_value + 1e-6`——在 bug 下依然通过（弱断言）。

**触发场景**：违约率恰等于 `1-confidence` 但违约集中于连续区间 → `LR_pof≈0` → χ²(2) 的 p≈1 → 不拒绝（正确检验应拒绝）。

**处置**：(1) 按上式修 `joint_test`，删无用变量；(2) 新增**失败用例**：等率聚集序列（`joint.reject_null == true`）；(3) 同步收紧现有 `JointTest` 断言。

**实施（2026-10-21）**：

- `joint_test` 重写：`LR_cc = LR_pof + LR_ind`，其中 `LR_pof` 由 `(n, v)` 的对数似然比给出、`LR_ind` 由**转移计数**（n00/n01/n10/n11）给出，两者共用与 `test()` 相同的 `independence_lr()`（单一实现，避免再次分叉）；被丢弃的 `pof_result`/`ind_LR` 伪变量（含用 `violation_rate` 冒充转移概率的表达式）全部删除，`LR_cc_total` 中间变量取消；
- 附带修正（同族口径问题）：`ChristoffersenIID::test` 返回的 `n_violations` 由 `n01+n11`（漏计首期）改为**全样本直接计数**，与 `KupiecPOF`/`joint_test` 口径统一 —— 首期为违约时旧值少计 1（`violation_rate` 分母为 n 而分子来自 n-1 样本，属同一类静默错数）。

**验证**（新增 3 用例 + 收紧 1 用例）：

| 用例 | 断言 | 修复前 | 修复后 |
|---|---|---|---|
| `JointTestSumsLRComponents` | 温和聚簇（LR_cc≈3.89）→ `p == exp(-LR_cc/2)`（χ²(2) 上尾闭式）±1e-12 | **失败**：LR_cc 仅含 LR_pof=0.19 → p≈0.67（期望 0.1426） | 通过 |
| `JointTestClusteredAtExpectedRate` | 5% 违约率 + 5 个 5 连簇 → `reject_null == true`、p<1e-6 | **失败**：LR_cc≈0 → p≈1 → 静默不拒绝 | 通过 |
| `JointTestIndependentCorrectCoverage` | 孤立违约 + 覆盖正确 → `p == exp(-LR_cc/2)` 且不拒绝 | **失败**：p≈0.951 vs 期望 0.9504（±1e-12） | 通过 |
| `JointTest`（收紧） | 聚簇序列 `ind/joint.reject_null == true` 且 `joint.p ≤ ind.p + 1e-12` | 失败（旧 p≈1） | 通过 |
| `ViolationCountsConsistentAcrossTests` | 首期违约序列：Kupiec/ind/joint 三者 `n_violations` 均为 3 | **失败**：`test()` 返回 2 | 通过 |

---

### RISK-024 — 回测似然连乘下溢（HIGH）

**事实**：`backtesting.hpp:58-60`（Kupiec）、`:107-114`（Christoffersen）直接计算 `L = Π p^x (1-p)^(N-x)`。

**数值核验（手算）**：

| 场景 | 量级 | 后果 |
|---|---|---|
| Christoffersen：p=0.05, n=3750, x≈187 | L ≈ 10^(-79.3-243.4) ≈ 10^-322.7 | 低于 `DBL_TRUE_MIN`(4.94e-324) 边缘 → **次正规、精度尽失**（非精确 0）；n≳3850 时为 0 |
| Kupiec：N=5000, x=2500 | L0≈10^-3308=0，L1=2^-5000≈10^-1505=0 | `log(0/0)` → **p=NaN**，`reject_null=false`（50% 违约率都不拒绝） |

**核验补充（评审未列的第二条静默通道）**：完美聚集时 `pi0=0`、`pi1=1` → `(1-pi1)^n10 = 0` → `L_nonIID = 0` → 命中 :112 的 `L_IID>0 && L_nonIID>0` 守卫 → `LR_ind` 直接置 0、p=1。**这条与下溢无关，改成对数似然不会自动修掉**，守卫语义需一并重写（区分"数值退化"与"模型完美拟合"）。

**处置**：(1) 全部改对数似然 `x·log(p) + (N-x)·log1p(-p)`，约定 `0·log0 = 0`；(2) 重写 :112 守卫；(3) 补大 N（≥3750）与聚集两类回归用例。

**实施（2026-10-21）**：

- `detail::xlogy` / `xlog1m` / `log_lik_binomial`（0·log0 = 0 约定）落地，三处似然（Kupiec POF、独立性 IID/非 IID）全部改为对数空间；
- **守卫重写 = 守卫删除**：原 `L_IID>0 && L_nonIID>0` 守卫无法区分"数值退化"与"完美拟合"。对数空间中"完美聚集"（pi0=0 / pi1=1）由约定自然给出 `logL_non = 0`（该模型以概率 1 拟合观测），`LR_ind = -2·logL_IID > 0` —— 证据保留，不再被误置 0；
- `KupiecPOF::test` 的 `v=0` / `v=N` 两个特判分支由同一约定覆盖（`xlogy(0, 0) = 0`、`xlogy(N, 1) = 0`），已删除；
- 输入边界（原实现静默产出 NaN/∞）：`n_observations == 0`、`n_violations > n_observations`、`confidence ∉ (0,1)` → `std::invalid_argument`；`ChristoffersenIID` 观测数 < 2（转移计数需 n-1 ≥ 1）→ `std::invalid_argument`。

**验证**（新增 4 用例）：

| 用例 | 断言 | 修复前 | 修复后 |
|---|---|---|---|
| `KupiecPOF.LargeSampleNoUnderflow` | v=2500/N=5000：p 有限、`reject_null`、p<1e-10 | **失败**：L0=L1=0 → `log(0/0)` = NaN → p=NaN、`reject_null` 恒 false | 通过 |
| `ChristoffersenIID.LargeSampleNoUnderflow` | n=4000 的 200 连簇：p 有限、`reject_null`、p<1e-10 | **失败**：L_IID≈1e-343 下溢 → 守卫置 LR_ind=0 → p=1 | 通过 |
| `KupiecPOF.RejectsInvalidInputs` | v>N / N=0 / confidence=1.0 各抛 `invalid_argument` | **失败**（不抛） | 通过 |
| `ChristoffersenIID.RejectsTooFewObservations` | n=1 与 n=0 各抛 `invalid_argument` | **失败**（pi = v/0 = NaN） | 通过 |

---

### RISK-025 — GPU PIMPL 泄漏（HIGH）

**事实**：

- `include/cpphub/performance/gpu/gpu_mc.hpp:110`：`Impl* impl_ = nullptr;`（裸指针）
- `src/performance/gpu/gpu_mc.cu:296`：`~GpuMCEngine() = default;` → :282-290 的 `~Impl`（含 `cudaStreamDestroy`）**永不执行**
- :298-299：默认移动仅拷贝指针 → 移出对象与移入对象共享同一 `Impl`
- `src/performance/gpu/gpu_mc_cpu_stub.cpp:41`：同一问题（stub 版本泄漏 `Impl` 对象）
- :293-294 另存 `config_` 与 `Impl::config` 重复

**核验补充**：当前析构为空，故失败模式是**泄漏**而非悬垂/双重释放；一旦改用 `unique_ptr` 而不同步修正移动语义，会从"泄漏"变成"双重释放"——修复须一并处理。

**处置**：`std::unique_ptr<Impl>` + 析构在 .cu/.cpp 内离线定义（`= default` 不能用于不完整类型）+ 删除重复 `config_`；补"循环构造 N 个引擎"的回归（可用 CUDA 内存查询或 stream 计数间接验证）。

**实施（2026-10-21）**：

- `gpu_mc.hpp`：`Impl* impl_` → `std::unique_ptr<Impl>`（+ `#include <memory>`）；删除与 `Impl::config` 重复的 `config_` 成员；构造期异常语义写入 `@throws` 契约；
- `gpu_mc.cu` / `gpu_mc_cpu_stub.cpp`：析构与移动保持离线 `= default`——`unique_ptr` 同时解决"释放"与"所有权转移"两个问题（避免评审 §核验补充指出的"改指针不改移动语义 → 双重释放"陷阱）；`Impl` 在**全部资源就绪后**才 `++` 计数，构造抛异常不产生计数偏差；
- 诊断钩子 `live_impl_count()`（`gpu_mc.hpp:87-89`）标注为"测试专用、非 API 契约"。

**验证**（新增 2 用例，`test_gpu_mc.cpp:320-357`）：

| 用例 | 断言 | 修复前（实机 RTX 4060） | 修复后 |
|---|---|---|---|
| `EngineDestructorReleasesImpl` | 作用域内 count==1 → 出作用域==0；循环 64 次后==0 | **失败**：泄漏只增不减，实测累计 72 个 Impl（含 72 个未销毁 stream） | 通过 |
| `EngineMoveTransfersOwnership` | 移动构造后==1；移动赋值释放 dst2 旧 Impl 后==1；全部析构后==0 | **失败**：count 单调递增（73→74），无任何释放 | 通过 |

---

### RISK-026 — CUDA 错误检查缺失与 block_size 未校验（HIGH）

**事实**（`src/performance/gpu/gpu_mc.cu`）：

- :268-280 `cudaSetDevice`/`cudaStreamCreate` 返回码未检查；
- :328-340 `cudaMalloc`/`cudaMemcpy`、:355 kernel launch、:359 `cudaStreamSynchronize`、:366-369 D2H 回拷**全部未检查**；
- :318 `result.used_gpu = true` 在 launch **之前**设置 → 失败链下调用方拿到 `price=0, SE=0, used_gpu=true` 的**伪结果**；
- :163-168 / :187-192 `__shfl_down_sync(0xFFFFFFFF, ...)` 在全掩码下要求 32 lane 全部活跃 → `block_size` 非 32 倍数（如 100）为**未定义行为**；
- `block_size=2048` 同时超过主流 GPU `maxThreadsPerBlock=1024` 与默认 48KB smem（4×2048×8B=64KB）→ launch 失败（静默）。

**处置**：(1) `CUDA_CHECK` 宏（失败抛 `std::runtime_error`，复用既有 stub 的异常语义）；(2) 构造期校验 `block_size % 32 == 0` 且 `≤ maxThreadsPerBlock`；(3) 设备内存 RAII（`unique_ptr` + 自定义 deleter）；(4) 失败时**重置** `used_gpu=false`；(5) warp 归约掩码改为 `__activemask()` 或按实际活跃 lane 构造。

**实施（2026-10-21）**（逐项对照，含一处对处置项的修正）：

| 处置项 | 落地 | 覆盖范围 |
|---|---|---|
| (1) 错误检查宏 | `CPPHUB_CUDA_CHECK`（`gpu_mc.cu:314-315`）→ `cuda_check()` 抛 `std::runtime_error`，附 `cudaGetErrorString` + 调用点文件:行 | 14 处检查点：`cudaGetDeviceCount` / `cudaSetDevice` / `cudaStreamCreate` / `cudaMalloc` / 4×H2D + 4×D2H `cudaMemcpy` / `cudaGetLastError`（launch 后）/ `cudaStreamSynchronize` |
| (2) 构造期校验 | `Impl` 构造：`block_size` 非 32 的正整数倍 → `std::invalid_argument`；`device_id` 越界（按 `cudaGetDeviceCount`）→ `std::runtime_error`；`block_size > device.max_threads_per_block` → `std::invalid_argument` | 校验全部**前置**，失败不创建任何资源（含计数器）；顺带关闭 48KB smem 隐患：`block_size ≤ 1024` ⇒ smem 最大 4×1024×8B = 32KB < 48KB 默认上限 |
| (3) 设备内存 RAII | `DeviceBuffer`（`gpu_mc.cu:280-307`）：构造 `cudaMalloc`（失败抛异常）、析构 `cudaFree`；删除原 4 处显式 `cudaFree` | 覆盖 `cudaMalloc` 半途失败与后续任何早退路径 |
| (4) `used_gpu` 语义 | 赋值移至全部 CUDA 调用成功之后（`gpu_mc.cu:469`） | 失败链不再返回 `price=0 + SE=0 + used_gpu=true` 的伪结果 |
| (5) warp 归约掩码 | **改为由校验保证全掩码安全**，未改用 `__activemask()` | 依据：`block_size` 被强制为 32 的正整数倍 ⇒ 归约点所有 lane 必然活跃，全掩码无 UB；`__activemask()` 在收敛后仅反映"当前活跃"，语义上更易引入新错误。已在 `gpu_mc.cu:165-167` 就地注明该不变量 |

**附带修正**（同族问题，一并处理）：析构中冗余且未检查返回值的 `cudaDeviceSynchronize` 移除（`cudaStreamDestroy` 本身保证流上残留工作完成后自动回收）；`Impl` 构造中 `stream_count` 双分支雷同的 `cudaStreamCreate` 合并为单次调用；设备属性查询失败（`gpu_query_device` 返回空 info）由隐性降级改为显式 `std::runtime_error`。

**验证**（新增 2 用例 + 全量回归）：

| 用例 | 覆盖 | 修复前 | 修复后 |
|---|---|---|---|
| `EngineRejectsInvalidBlockSize` | block_size = 100 / 0 / -32 / 4096 各抛 `std::invalid_argument` | **失败**（不抛，静默接受） | 通过 |
| `EngineRejectsUnknownDevice` | device_id = 9999 抛 `std::runtime_error` | **失败**（静默伪结果） | 通过 |

平台与结果：RTX 4060（SM 8.9）+ nvcc 13.2 + MSVC 19.50；`test_gpu_mc` **19/19 通过**（含既有 15 用例，定价/确定性/SE 收敛/性能全绿，证明 RAII 重构未改变数值路径）；CPU stub 路径（`CPPHUB_ENABLE_CUDA=OFF`）`cpphub_cuda` 编译通过。

---

### RISK-027 — RNG 名实不符（MEDIUM，待裁定）

**一次源比对**（Random123 `philox.h`，master，2026-10-21 取得）：

| 项 | Random123 原文 | Cpp_Hub `core/rng.hpp` | 结论 |
|---|---|---|---|
| 4x64 乘数 0 | `PHILOX_M4x64_0 = 0xD2E7470EE14C6C93` | `PHILOX_M4x64_0 = 0xD2B74407B1CE6E93`（:18） | 该值实为 **`PHILOX_M2x64_0`** → 命名错误 |
| 4x64 乘数 1 | `PHILOX_M4x64_1 = 0xCA5A826395121157` | `0x9CCD9A8C5AE13B41`（:19） | **philox.h 全文不存在**该常量 |
| key 增量 | 64 位引擎用 `W64_*`（`0x9E3779B97F4A7C15`/`0xBB67AE8584CAA73B`） | `W32_*`（:20-21，32 位 Weyl） | 用于 64 位引擎，非标准 |
| 轮函数 | 2x64：`{hi(M2x64_0·c0)^k0^c1, lo}`（单乘法、单 key 字）；4x64：4 字状态 | `c0'=hi1^lo0^k0, c1'=hi0^lo1^k1`（:82-100，2 字状态） | **两者都不是**；可逆性未见论证 |
| KAT | 官方 KAT 向量 | `tests/unit/core/test_rng.cpp` 9 用例全为自洽/统计型 | 无已知答案测试 |
| CPU/GPU 流一致性 | — | GPU 取 counter=path 的 (r0,r1)（`gpu_mc.cu:132-136`）；CPU `normal_simd` 取 4 个连续 counter 的 r[0]（`rng.hpp:113-117`） | 同 seed 下**正态数流不同**；`README.md:40` 的 "bit-exact match with CPU" 不成立；`test_gpu_mc.cpp` 无 CPU/GPU 逐样本比对用例（grep 零命中）|

**裁定**：这是**自定义 2 字计数器式 RNG**，与 Random123 Philox 系列无对应关系。

**处置（二选一，需裁决）**：

- (a) 换标准 Random123 Philox + 官方 KAT + CPU/GPU 逐样本比对 → **代价**：所有依赖该 RNG 的 MC/GPU 数值基线需整体重算；
- (b) 保留现状，改名为自定义 RNG + 修正 `README.md:40` 与 `gpu_mc.hpp:10-13` 的位精确表述 + 补统计检验基线。

---

### RISK-028 — χ² CDF 近似（MEDIUM）

**事实**：`backtesting.hpp:14-22` 用 Wilson-Hilferty 近似；`include/cpphub/core/math.hpp:100` 已有精确 `regularized_lower_gamma(a, x)`（Numerical Recipes §6.2，精度 ~1e-14），未被使用。

**数值核验（k=1 手算）**：

| LR | 代码 p | 真值 p | 后果 |
|---|---|---|---|
| 2.706 | 0.0958 | 0.1000 | 系统偏低 ~0.004 |
| 3.800 | 0.0486 | 0.0513 | **翻转**（代码拒绝、真值不拒绝） |
| 3.841 | 0.0473 | 0.0500 | 临界点附近偏离 0.0027 |

**处置**：`chi2_cdf(x, k) = regularized_lower_gamma(k/2, x/2)`（一行替换）；补临界值断言用例（LR=3.841 → p=0.05±1e-6）。

**实施（2026-10-21）**：按上述一行替换落地（`backtesting.hpp` detail 层），`x ≤ 0` 仍返回 0（保证数值上略负的 LR 不触发 `regularized_lower_gamma` 的 `a>0, x≥0` 前置抛错）。

**验证**（新增 3 用例，全部覆盖旧近似的两处失效面）：

| 用例 | 断言 | 修复前 | 修复后 |
|---|---|---|---|
| `Chi2Cdf.MatchesExactGammaAndClosedForm` | k=1 闭式 `CDF(x)=erf(√(x/2))` ±1e-14；标准临界值（3.841459→0.95、2.705543→0.90、6.634897→0.99、5.991465(k=2)→0.95、9.487729(k=4)→0.95）±1e-6 | **失败**（Wilson-Hilferty 偏差 ~2.8e-3） | 通过 |
| `KupiecPOF.PValueMatchesExactChi2Tail` | 6/250@99% → `p == erfc(√(LR/2))` ±1e-10 | **失败**：旧 p≈0.0569 vs 精确 0.0593（偏差 2.4e-3） | 通过 |
| `ChristoffersenIID.PValueMatchesExactChi2Tail` | 温和聚簇（LR_ind≈3.71）→ `p == erfc(√(LR_ind/2))` ±1e-10 | **失败**：旧 p≈0.0513 vs 精确 0.0542 | 通过 |

**同族核查（E3，确认本库其余 χ² 路径无需处理）**：`econometrics/core/special_functions.hpp:176-193` 的 `chi2_sf`/`chi2_cdf` 已是**精确实现**（df=1 用 `erfc(√(x/2))`、df=2 用 `exp(-x/2)`、一般 df 用 `gammq` 连分式），且该模块的 `chi2_ppf_approx` 虽为 Wilson-Hilferty，但已显式限定用途（"仅用于临界值显示，不用于精确推断"）→ 静默错风险不存在，本次不改。**副产物**：上述 df=1/df=2 闭式与本清单测试所用作独立参考（`erfc`/`exp`）逐式一致，构成跨来源佐证。

**范围说明（本次未做，非评审项）**：上尾仍按 `1 - P` 计算 ⇒ p ≲ 1e-16 时归零（**不影响**拒绝/不拒绝判定，因 LR 极大时两种取值都判拒绝）。精确上尾需 `Q(a,x)` 直算（在 `core/math.hpp` 暴露精确上不完全 Gamma），属 P0 之外，已登记 **RISK-036** 跟踪。

---

### RISK-029 — 架构指令集硬编码（MEDIUM）

**事实**：`CMakeLists.txt:15-21` 在 `CMAKE_BUILD_TYPE == "Release"` 时无条件追加 `/arch:AVX2`（MSVC）/ `-march=x86-64-v3`（其余）。**评审措辞"无条件"不准确——有 Release 条件**，但**未按架构条件化**这一定性成立：aarch64/Apple Silicon 上 `-march=x86-64-v3` 直接编译失败；生成物在不支持 AVX2 的 CPU 上 SIGILL；`README.md:78` 声明支持 Clang 14+；`include/cpphub/core/cpu_features.hpp:26-63` 有运行时检测但与此处编译选项未联动。

**范围界定**：`python/CMakeLists.txt`（wheel 构建路径）不含这些选项 → **wheel 分发不受影响**；受影响的是根工程 Release 构建。

**处置**：`option(CPPHUB_NATIVE_ARCH "Enable x86-64-v3/AVX2 codegen" OFF)` 控制，按 `CMAKE_SYSTEM_PROCESSOR` 判平台；如需性能，走运行时派发（`cpu_features.hpp` 已在位）。

---

### RISK-030 — PDE `mutable` 计数器数据竞争（MEDIUM）

**事实**：`include/cpphub/pricing/pde/pde_engine.hpp:259` `mutable Size last_total_iter_count_`；写点 :91、:125（均在 `const` 方法 `price_american` 内）；读点 :255 `last_total_iterations()`。两个线程对同一 `const PDEEngine` 并发调用 `price_american` 即数据竞争（UB），读者拿到不确定值。

**处置**：迭代次数并入返回值（如 `PricingResult{price, iterations}`），或最少改为 `std::atomic<Size>`。属 API 变更，走 ADR 记录。

---

### RISK-031 — PDE 边界字符串分派与 Neumann 失效（MEDIUM）

**事实**：

- `pde_engine.hpp:190-203`：`price_barrier` 按 `payoff.name() == "Call"` 二选一，否则一律 `bsm_put` 远场边界；可见 payoff 名称为 `"Call"`/`"Put"`/`"DigitalCall"`/`"DigitalPut"`/`"DoubleDigital"`（`include/cpphub/instruments/payoff/vanilla.hpp:20/33/48/63/79`）→ **除 Call/Put 外全部落到 put 远场**（数字期权、自定义 payoff 结果有偏）；
- :143-144 注释声明"Dirichlet 或线性外推 (Neumann)"，但 `price_barrier`（:181-208 循环体）**从不调用 `finalize_boundary`**；Neumann 仅在 :285 与 :308-312 生效 → `BoundaryType::Neumann` 在障碍路径被**静默忽略**。

**处置**：渐变式边界改由 `PayOff` 提供（虚函数）或显式传入边界策略；Neumann 要么实现，要么在障碍路径显式报错（拒绝静默降级）。

---

### RISK-032 — 文档/版本/CI 一致性（LOW，部分成立）

| 子项 | 事实 | 状态 |
|---|---|---|
| README GPU 测试数 | `README.md:107` "15 GPU MC tests still run on CPU" vs `src/performance/gpu/CMakeLists.txt:93,113-114`（无 CUDA 时 `test_gpu_mc` **不构建**）+ `gpu_mc_cpu_stub.cpp:58,67`（stub 直接抛异常） | 成立，待修 |
| 版本号三处不一致 | `CMakeLists.txt:2`=`1.0.0`；`python/pyproject.toml:24`=`1.4.0`；`python/src/cpphub_python.cpp:2043`=`1.4.0`；git tag 已到 `v1.8.0` | 成立，待修 |
| CI 无 Python job | **已由本会话 S5 修复**：`.github/workflows/ci.yml:132-177` 新增 `python` job（windows-2022 / ubuntu-latest × py3.11 / py3.12） | **已关闭** |
| 无覆盖率 | CI 无 gcov/llvm-cov、无阈值 | 成立，待办 |
| 文件过长 | `include/cpphub/calibration/optimizer.hpp` **1244 行**、`include/cpphub/c_api/cpphub_tool_registry.hpp` **1160 行**（逐字命中） | 成立，待决策 |

**处置**：版本号单源化（`project(VERSION)` → `configure_file` 注入 pyproject/`__version__`）；修 README；CI 加覆盖率；大文件拆分与"拆子项目"建议**不采纳物理拆分**（见 §四）。

---

### RISK-033 — GoogleTest FetchContent 无哈希（安全）

**事实**：`CMakeLists.txt:63-66` 以 tag tarball URL 拉取 googletest v1.14.0，**无 `URL_HASH`**。GitHub tag 可被重指，供应链面非空。

**处置**：加 `URL_HASH SHA256=<...>`（或改 `GIT_TAG <完整 40 位 SHA>`）；两者取其一，勿只依赖 tag 名。

---

### RISK-034 — `json_min` 外部输入面（安全，待复核）

**事实**：`include/cpphub/c_api/detail/json_min.hpp` 已有**深度上限**（`:27 kMaxParseDepth = 64`）与**尾随垃圾拒绝**（`:166`）；但 ABI 入口 `cpphub_v1_tool_call` 对 `args_json` **无总长上限**——长度校验在解析之后按语义字段进行（`kMaxSeriesLength=10000` / `kMaxColumns=64` 在 `cpphub_tool_registry.hpp:50-52`）。

**风险定性**：本地 stdio 场景，现实风险低；风险面是"超大输入的内存放大"。**移交 security-review 复核**（与既有 `TRAE-security-review` 流程），本清单只登记不结论。

---

### RISK-035 — （本库增补）`chow_denning_joint` 的 `statistic` 语义漂移（P2，待裁定）

**事实**：`include/cpphub/c_api/cpphub_tool_registry.hpp:540-544` 把 `r.chow_denning_p_value` 同时传给 `statistic` 与 `p_value` 两个字段；库内 `VarianceRatioResult` 只有 `chow_denning_p_value`（无联合统计量）。

**状态**：已在 `spec/mcp-server/CHECKLIST.md` §1.5 登记为 F9；`python/tests/fixtures/mcp_golden.json` 的 `meta.known_issues` 已标注"冻结值反映现状而非期望语义"。

**待裁定**：改为 `statistic=null` + note 注明"仅 p 值"（工具契约变更 + golden 重生成），或维持现状并在描述中说明。

---

### RISK-036 — （本库增补）χ² 上尾精度下限（P3，已登记）

**事实**：`backtesting.hpp` 的 p 值一律按 `1 - chi2_cdf(LR, k)` 计算。CDF 侧已是精确值（`regularized_lower_gamma`），但上尾由减法得到 ⇒ 当 `LR` 大（p ≲ 1e-16，约 `LR ≳ 74` @k=1）时 `1 - P` 在双精度下舍入为 0。

**影响界定**：**不影响任何判定**——该区间内无论 p 取 0 还是 1e-45，`reject_null` 均为 true；影响面仅限于"p 值量级不可读"（报告里显示 0 而非 1e-45）。本次 RISK-023/024/028 的回归用例已按此界定设计（深度尾部只断言判定与量级上界，不断言精确值）。

**处置候选**（待裁定）：在 `core/math.hpp` 暴露精确上不完全 Gamma `Q(a,x) = regularized_upper_gamma`（现实现为 `1 - P`，NR §6.2 的连分式分支本就直接算出 Q，可无成本返回），随后 `p = Q(k/2, LR/2)`。**代价**：触碰 core 头文件（全库重编译），收益仅在 p < 1e-16 的可读性。

---

## 四、依赖关系纠正与建议采纳情况

**一处上游建议被证伪**：曾有意见将 RISK-023/024/028 归为"MCP 直接相关、S6 前必修"（理由：诊断工具与回测共用 `backtesting.hpp`）。**读码核对不成立**——`backtesting.hpp` 的消费方仅 `python/src/cpphub_python.cpp:66`（Python 绑定）与两处测试；MCP 链路（`cpphub_tool_registry.hpp` → `garch_diagnostics.hpp`）**不含**它。故这三项是**库质量/风险 API 问题**，不是 S6 阻塞项；是否在 S6 前修属优先级选择而非依赖约束。

**评审建议采纳情况**：

| 建议 | 处置 |
|---|---|
| 先修 `backtesting.hpp` 三件套 + 补回归 | **采纳**（P0） |
| 重构 GPU（unique_ptr / CUDA_CHECK / block_size） | **采纳**（P0） |
| RNG 换标准 Philox 并加 KAT | **采纳但需裁决**（代价 = MC 基线重算；见 RISK-027） |
| 构建去 `-march` + CI 增 Python/覆盖率/ARM 矩阵 | **部分采纳**：Python job 已在 S5 完成；其余见 RISK-029/032 |
| 统一版本号、修 README | **采纳**（P1） |
| 拆分 >800 行文件 | **部分采纳**：`optimizer.hpp` 可拆；`cpphub_tool_registry.hpp` 是"单一登记点"设计的产物，拆分会破坏 7 工具集中登记（ADR-015/设计 §3.3），倾向保留并加区间注释 |
| 把计量/时序/MCP 拆成独立子项目或子库 | **不采纳物理拆分**：ADR-013 的 INTERFACE 分层（`cpphub_core` / `cpphub_econometrics` / `cpphub_timeseries_mat`）已完成依赖隔离；物理拆分与 header-only 单树分发、单一 wheel 冲突。改为"CI 分层 job" |

---

## 五、优先级排序

| 优先级 | 编号 | 动作 | 理由 |
|---|---|---|---|
| ~~P0~~ **已完成** | RISK-023 / 024 / 028 | `backtesting.hpp` 三件套（对数似然 + 精确 χ² + LR_cc=LR_pof+LR_ind + 守卫删除 + 输入边界），2026-10-21 落地，`test_backtesting` 21/21 + 相邻三套件全绿 | 静默错——统计检验给出**错误且不报错**的结论 |
| ~~P0~~ **已完成** | RISK-025 / 026 | GPU PIMPL + CUDA_CHECK + block_size 校验（2026-10-21 落地，`test_gpu_mc` 19/19 通过） | 资源泄漏 + 伪结果（`used_gpu=true` 的 0 值） |
| P1 | RISK-027 | RNG 裁决（换实现 or 更名） | 影响 README 承诺与 MC 基线；换实现需冻结一次基线重算 |
| P1 | RISK-032 | 版本号单源化 + README 修正 + CI 覆盖率 | 一致性/可信度面 |
| P1 | RISK-033 / 034 | 加 `URL_HASH`；`json_min` 移交安全复核 | 供应链与输入面 |
| P2 | RISK-029 / 030 / 031 | 架构开关化；PDE 迭代计数与边界策略重构 | 平台/API 变更，走 ADR |
| P2 | RISK-035 | F9 语义裁决（`statistic=null` + golden 重生成，或维持现状） | MCP 工具契约 |
| P3 | RISK-036 | χ² 上尾精确化（`Q(a,x)` 直算） | 仅影响 p<1e-16 的可读性, 判定不受影响 |

---

## 六、证据边界与未验证项（诚实声明）

1. **RISK-025/026 已补实机复现**（2026-10-21，RTX 4060）：泄漏链以 `live_impl_count()` 计数复现（修复前：64 次循环构造 + 前序用例累计 72 个 Impl 未释放），配置校验链以 4 组非法 `block_size` + 越界 `device_id` 用例复现。**仍未实测**的分支：launch 期错误（如显存耗尽）与 `cudaMalloc` 失败路径——需构造真实失败场景，当前仅由代码路径覆盖保证，未做故障注入。
2. **RISK-027 的"轮函数不可逆"**：仅确认"**没有可逆性论证**"（评审原话即如此），未给出可逆性证明或反例。
3. **行号偏移记录**（不影响裁定）：`gpu_mc.cu` 析构在 :296（评审写 297）、移动在 :298-299（评审写 299-300）；`pyproject.toml` 版本在 :24（评审写 22）。其余抽查项逐字命中。
4. **RISK-024 的次正规判定**：n=3750 时是"精度尽失的次正规数"而非精确 0（阈值 `DBL_TRUE_MIN=4.94e-324`）；评审"下溢为 0"的方向正确、幅度略乐观；n≳3850 时为精确 0。
5. **RISK-034 未做深度安全分析**：仅登记输入面与移交，结论以 security-review 结果为准。
6. **RISK-023/024/028 的验证边界**（2026-10-21）：用例全部为**构造序列**（确定性人工违约模式）+ 闭式参考（χ² 尾用 `erfc` / `exp` 独立计算），未接入真实回测数据；LR_cc = LR_pof + LR_ind 的分解已数值验证到 1e-12（经 χ²(2) 闭式间接观测，因为该 API 不暴露统计量本身），**未与 R 包**（`rugarch::vaRtest` 等）**逐一对照**。另：Python 绑定为静态编译产物，wheel 需重建后才包含本次修复（CI 的 python job 覆盖该路径，本地未重建）。

---

## 七、计数与维护约定

**计数（RULE-2 逐项口径）**：外部评审 12 条（HIGH 4 / MEDIUM 5 / LOW 1 / 安全 2）+ 本库增补 2 条（RISK-035 / 036）= **登记 14**；其中 成立 12 / 部分成立 1（RISK-032）/ 不成立 0；**已关闭子项 1**（CI Python job，S5）；**已修复 5**（RISK-023 / 024 / 025 / 026 / 028，2026-10-21）。

**修改历史**：

| 日期 | 变更 |
|---|---|
| 2026-10-21 | 建立核验清单（RISK-023~035 登记） |
| 2026-10-21 | RISK-025 / RISK-026 修复落地并实机验证，状态转"已修复"（§三 各附实施/验证记录） |
| 2026-10-21 | RISK-023 / RISK-024 / RISK-028 修复落地（`test_backtesting` 21/21 + 相邻 3 套件全绿），状态转"已修复"；新增 RISK-036（χ² 上尾精度下限，P3）；消费方声明更新（`USAGE_GUIDE` §3.9.2 去除 Wilson-Hilferty 描述） |
| 2026-10-21 | P1 四项（027/032/033/034）调研完成 → [spec/p1-fixes/RESEARCH.md](../../spec/p1-fixes/RESEARCH.md)；主表状态转"调研完成，待裁决"，未决项 U-1~U-4 见该文件 §6.4 |
| 2026-10-21 | P1 设计完成（Step 3）→ [spec/p1-fixes/DESIGN.md](../../spec/p1-fixes/DESIGN.md)；U-1~U-5 全部裁定（A2 vendor Random123 / 版本两轴 / 覆盖率仅报告 / 长度+节点预算） |
| 2026-10-21 | P1 设计 Review（Step 4）完成：DESIGN 全部 `文件:行` 逐条重放 + 外部 E1 复核（GitHub API tag→commit `f8d7d77c…8571`、raw philox.h L175/L207-L247 全真），8 处行号/计数级勘误回写（无裁决性错误）→ DESIGN §12 |
| 2026-10-21 | P1 实施文档（Step 5）完成 → [spec/p1-fixes/IMPLEMENTATION.md](../../spec/p1-fixes/IMPLEMENTATION.md)（三构建目录矩阵 / U-6/U-8/U-10 探针与退路 / 四项 diff 级设计 / TDD 红→绿矩阵 / 4 提交序列 / A1-A10 映射）；两处设计勘误连带回写 DESIGN（E-1 KAT 夹具 .hpp 化→§11 D-4 / E-2 P3 定性为回归锁），新增 GPU 测试探针 rng_device_probe 细化 |
| 2026-10-21 | P1 实施 Review（Step 6）完成：IMPLEMENTATION 全部 `文件:行` 重读 + 代码草稿静态推演（编译/越界/进制/边界算术）+ 与 DESIGN/RESEARCH 对齐；**13 处勘误回写**（3 处真编译错误：重复 namespace alias、探针先用后定义/命名空间错层、test_gpu_mc 缺 include；1 处越界写；1 处 py3.9 兼容；3 处 CI 可执行性：缺 `--check` 门禁、`python3`、artifact glob；2 处文档自洽；1 处幻觉 URL 清除；其余口径），无裁决性错误 → IMPLEMENTATION §10 |
| 2026-10-21 | P1 验收清单（Step 7）完成 → [spec/p1-fixes/CHECKLIST.md](../../spec/p1-fixes/CHECKLIST.md)：G 6 项（前置门禁 U-6/7/8/10、基线、离线 gtest）+ C 78 项（R-033 6 / R-032 19 / R-034 12 / R-027 28 / 回归性能 6 / 提交回滚 7）= 84 项；含硬约束（禁放宽容差、禁改夹具、禁退回 .txt、禁恢复自写 Philox）；A1-A10 全覆盖映射 |
| 2026-10-22 | P1 文档一致性对齐（Step 8）完成：RESEARCH U-1~U-5 ↔ DESIGN 5/5、A1-A10 三文档 10/10/10、测试名 9/9 归一化核对；**7 处不一致闭环**（1×P2：test_rng 既有用例 8→**9**，RESEARCH E1-6 早已记 9；6×P3：DESIGN §6.2 行区间 :25-29→:25-28、DESIGN §4.4 与 IMPL 的 `python3`/`coverage*.html` 对齐并补 `--check` 门禁、D-04/D-05 口径归一化、C-79 补全测试名、§1.1 URL 以 RESEARCH §7.1 填实）→ CHECKLIST §1 |
| 2026-10-22 | P1 实施（Step 9）G 系列门禁全部通过：U-7 裁定 **1.9.0**；U-8 vendor Random123 commit `9545ff64…1150`（codeload 直连，28 文件 + LICENSE + UPSTREAM.md）；U-6 KAT 列序确认（实际在 `tests/kat_vectors`，10 轮 3 行）+ 夹具生成（第 1 行答案与 device 输出交叉验证）；**U-10 探针发现 Windows+nvcc 下 `MULHILO64_MSVC_INTRIN` 遮蔽 `CUDA_INTRIN`（philox.h L199-202）→ device 编译失败**，修复 = device pass 宏守卫（DESIGN §11 D-6，sm_89 验证输出 = 官方 KAT）；基线三目录全绿（2605 / 2587 / 2548）；R-033 形态经裁决改为 **URL@commit**（本机 git 协议不通，D-5）；实施期勘误 E-6~E-10 见 IMPLEMENTATION §9 |

**复算命令**（可重放）：

```powershell
# 主表行数（RISK-023~036，应为 14）
(Get-Content docs\audit\EXTERNAL_REVIEW_2026-10-20.md | Select-String -Pattern '^\| RISK-0\d\d ').Count
# 行数引用复核
(Get-Content include\cpphub\calibration\optimizer.hpp).Count
(Get-Content include\cpphub\c_api\cpphub_tool_registry.hpp).Count
```

**维护约定**：

- 每条 issue 修复时：在本文件 §二 主表更新状态，并在 `docs/DEVELOPMENT_LOG.md` 追加修复条目（沿用 RISK-NNN 编号）；
- 关闭时须附**验证方式与结果**（回归用例名 + 平台），与 RISK-017~022 的既有格式一致；
- RISK-032 的 CI 子项与 RISK-035 已分别回写 `spec/mcp-server/CHECKLIST.md` §7 与 §1.5。

---

**核验人**: Scott（鹏） · **核验日期**: 2026-10-21 · **评审快照**: 约 2026-10-19~20 · **最近修复**: 2026-10-21（RISK-023/024/025/026/028）
