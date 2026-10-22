# 调研报告：P1 整改（RISK-027 / 032 / 033 / 034）

> **Feature**: p1-fixes（外部评审 P1 四项整改）
> **Spec 步骤**: Step 1–2（调研 → Review 形成调研文档）
> **规范依据**: Spec_Workflow 10 步流程 + [ASSERTION_EVIDENCE_FRAMEWORK.md](../../docs/ASSERTION_EVIDENCE_FRAMEWORK.md)
> **上游**: [EXTERNAL_REVIEW_2026-10-20.md](../../docs/audit/EXTERNAL_REVIEW_2026-10-20.md) §五 优先级表（P1 行）
> **状态**: **调研完成（Step 1–2）** — U-1~U-5 **全部已裁定**（2026-10-21），可进入 DESIGN（Step 3）
> **v2 变更（2026-10-21）**: 补外部锚点（Random123 / TensorFlow / JAX / NumPy / NVIDIA / cuRAND / simdjson / Jackson CVE）+ RISK-027 与 RISK-034 的**多维度对比**；新增 CPU/GPU 位精确实测（E1-9）
> **调研日期**: 2026-10-21

---

## §0 断言统计表

| 证据等级 | 条数 | 说明 |
|---|---|---|
| **E1** 一次源（页面级） | 1 | Random123 `philox.h` 全文抓取，逐行核对常量与轮次 |
| **E1** 可复算命令 / 本地实测 | 9 | 归档哈希、tag→commit、include 扇出、用例计数、全量 ctest、**CPU/GPU 位精确探针** |
| **E3** 静态读码 | 14 | 实现结构、文档承诺、版本六源、CI 配置、输入面保护、许可证与 third_party |
| **检索级**（厂商/社区文档） | 5 | Random123 README、NVIDIA 浮点白皮书、simdjson 限制、Jackson CVE 三例 |
| **估算（未实测，不得进入结论）** | 1 | `json_min::Value` 内存放大倍数 |
| **未核实（明示边界）** | 4 | 见 §5.2 |

**可复算命令**（E1，供审计端独立重数）：

```powershell
# E1-1 googletest v1.14.0 归档哈希（867764 字节）
Invoke-WebRequest -Uri https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz -OutFile $env:TEMP\g114.tar.gz
(Get-FileHash $env:TEMP\g114.tar.gz -Algorithm SHA256).Hash.ToLower()

# E1-2 tag → commit（单行输出 = lightweight tag，无 ^{} peeled 行）
git ls-remote https://github.com/google/googletest refs/tags/v1.14.0

# E1-3 版本/标签现状
git tag --sort=-v:refname | Select-Object -First 3
git describe --tags --always

# E1-4 rng.hpp 直接 include 扇出（应 76 = 42 头 + 32 测试 + 2 实现）
(Select-String -Path (Get-ChildItem -Recurse -Include *.hpp,*.cpp,*.cu -Path include,src,python/src,tests,benchmarks) -Pattern 'include\s+"cpphub/core/rng\.hpp"').Count

# E1-5 GPU 测试用例数（应 19）
(Select-String -Path tests\unit\performance\test_gpu_mc.cpp -Pattern '^TEST\(').Count

# E1-6 RNG 用例数（应 9）
(Select-String -Path tests\unit\core\test_rng.cpp -Pattern '^TEST\(').Count

# E1-7 全量回归观测（build/ 目录，C API=OFF）
ctest --test-dir build -C Release

# E1-8 一次源常量（Random123 philox.h，行号为该文件）
#   L207 PHILOX_M2x64_0 = 0xD2B74407B1CE6E93
#   L210 PHILOX_M4x64_0 = 0xD2E7470EE14C6C93
#   L213 PHILOX_M4x64_1 = 0xCA5A826395121157
#   L225 PHILOX_W64_0   = 0x9E3779B97F4A7C15
#   L231 PHILOX_W32_0   = 0x9E3779B9
#   L247 PHILOX4x64_DEFAULT_ROUNDS = 10

# E1-9 CPU/GPU 位精确实测（临时探针 tmp/probe_rng_bit.cu，未入库；nvcc 13.2 / sm_89）
#   n=20000, seed=42, counter=path: uniform-word mismatches = 0
#                                   z1(Box-Muller) bit mismatches = 2870 (14.35%)
```

---

## 1. 调研目标

对 P1 四项给出**可裁决的事实基础**，不做实现：

| 项 | 问题 |
|---|---|
| RISK-027 | 现有 RNG 到底是何物？换标准 Philox 的真实代价是多少？改名代价多大？两者谁更破坏契约？ |
| RISK-032 | 版本号到底有几个真值源、各是多少？README 哪些数字已失真？CI 有没有覆盖率/CUDA/ARM？ |
| RISK-033 | 唯一的外部依赖是什么？如何钉死（URL_HASH / commit SHA / vendor）？ |
| RISK-034 | 外部 JSON 输入的现有保护边界在哪？放大风险是否真实、威胁模型是什么？ |

---

## 2. 调研方法

- **E1**：`git ls-remote` 解析 tag → commit；本地下载归档后 `Get-FileHash`；`Select-String` 计数；**一次源抓取**（Random123 `philox.h`）；**本地实测探针**（CPU/GPU 逐位比对，`tmp/probe_rng_bit.cu`，用后删除）。
- **E3**：直接读码（`include/cpphub/core/rng.hpp`、`src/performance/gpu/gpu_mc.cu`、`include/cpphub/c_api/detail/json_min.hpp`、`src/c_api/c_api.cpp`、`CMakeLists.txt`、`python/*`、`.github/workflows/ci.yml`、`README.md`、`docs/USAGE_GUIDE.md`、`docs/architecture/ARCHITECTURE.md`、`LICENSE`、`third_party/`）。
- **运行观测**：`ctest --test-dir build -C Release`（2557 注册 / 2554 构建并通过，3 个 `NOT_BUILT` 为 `CPPHUB_ENABLE_C_API=OFF` 下的占位，非回归）。
- **检索级**：厂商/社区文档（Random123、TensorFlow、JAX、NumPy、NVIDIA、cuRAND、simdjson、NVD CVE），仅用于外部锚点，**不单独支撑任何结论**（见 §5.2）。

---

## 3. 调研发现

### 3.1 RISK-027 — RNG 名实不符

#### 3.1.1 实现事实（E3）

| 事实 | 证据 |
|---|---|
| 类名 `Philox4x64`，但状态是 **2 字计数器 + 2 字密钥** | [rng.hpp:24-25](../../include/cpphub/core/rng.hpp#L24-L25) `counter_[2]` / `key_[2]` |
| `compute()` 产出 **2 字**，`operator()` 只取 `r[0]` | [rng.hpp:83-101](../../include/cpphub/core/rng.hpp#L83-L101)、[:36-40](../../include/cpphub/core/rng.hpp#L36-L40) |
| 常量 0 取的是 **Philox2x64 的乘数**；常量 1 在 Random123 `philox.h` 全文不存在 | [rng.hpp:18-19](../../include/cpphub/core/rng.hpp#L18-L19)（一次源比对已在 [EXTERNAL_REVIEW §三 RISK-027](../../docs/audit/EXTERNAL_REVIEW_2026-10-20.md) 完成） |
| 密钥增量用 **32 位 Weyl 常量**（`0x9E3779B9` / `0xBB67AE85`）而非 64 位引擎的 `W64_*` | [rng.hpp:20-21](../../include/cpphub/core/rng.hpp#L20-L21)、[:96-97](../../include/cpphub/core/rng.hpp#L96-L97) |
| 轮函数为 **双乘数交叉混合**（`c0' = hi(c1·M1) ^ lo(c0·M0) ^ k0`），既非 Random123 2x64 轮（单乘法 + `^c1`），也非 4x64 轮（4 字状态） | [rng.hpp:89-98](../../include/cpphub/core/rng.hpp#L89-L98) |
| `next4()` = **连续推进 4 次计数器**、各取首字 | [rng.hpp:48-58](../../include/cpphub/core/rng.hpp#L48-L58) |
| GPU 端同算法的独立副本，counter = `[path, 0]`，只用 `z1` | [gpu_mc.cu:40-75](../../src/performance/gpu/gpu_mc.cu#L40-L75)、[:131-144](../../src/performance/gpu/gpu_mc.cu#L131-L144) |

**结论**：这是一个**自定义 2 字计数器 RNG**，与 Random123 Philox 系列无对应关系（可逆性亦无论证）。

#### 3.1.2 消费面（E1，76 个文件）

| 类别 | 数量 | 代表 |
|---|---|---|
| `include/` 头文件 | **42** | `models/diffusion/*`（8）、`pricing/monte_carlo/*`（9）、`timeseries/{garch,arima,var,midas}/*`（9）、`econometrics/resampling/*`（5）、`instruments/credit/*`（3）、`models/ir/*`（3）、`risk/{xva,pfe_sa_ccr,mc_var}`（3）、`monte_carlo/*`（2） |
| 实现文件 | **2** | `src/c_api/c_api.cpp:11`、`python/src/cpphub_python.cpp:29` |
| 测试/基准 | **32** | `tests/unit/core/test_rng.cpp`、`tests/unit/performance/benchmark_philox.cpp` 等 |

Python 公共接口已暴露该类：[cpphub_python.cpp:311](../../python/src/cpphub_python.cpp#L311) `nb::class_<Philox4x64>(m, "Philox4x64", ...)`。

#### 3.1.3 数值依赖分类（决定"换实现要重算多少"）

| 类别 | 结论 | 证据 |
|---|---|---|
| 官方 KAT（known-answer test） | **不存在**。`test_rng.cpp` 9 个用例全为自洽/统计型 | [test_rng.cpp](../../tests/unit/core/test_rng.cpp)（E1-6 = 9） |
| RNG 输出黄金向量 | **不存在**。含 `rng.hpp` 的测试中 `EXPECT_DOUBLE_EQ` 共 12 个文件，抽样核对全部锚定**解析常数/确定性变换/同 seed 复现**，而非 RNG 样本统计 | `test_zivot_andrews.cpp:173` `critical_1pct=-5.34`；`test_ng_perron.cpp:91-93` 临界值；`test_garch_m_model.cpp:94-96` 确定性映射；`test_sabr.cpp:491` 同 seed 路径逐元素相等 |
| 基准夹具（`tests/fixtures/`） | **不受影响**：由 Python/R 脚本（`gen_*.py` / `verify_*.py`）生成，与 C++ RNG 无关 | `tests/fixtures/{econometrics,hfe,timeseries}` |
| 统计容差型断言 | **需全量复跑**；固定 seed 的样本窗口（如违规率 ∈ (0.004, 0.020)）理论上可能越界，属"需重选 seed"而非实现错误 | `test_integration_phase7b.cpp:122-125` |
| 跨平台位精确承诺 | 与 RNG 无关（同一实现、同一编译选项下自洽），换 RNG 不改变"位精确可复现"性质 | `README.md:121`、`ARCHITECTURE.md:248` |

**关键判定**：换实现**不需要**重算任何外部数值基准（R/Python 夹具），只需**全量复跑**并处理少量种子敏感窗口。这与"换 RNG 会让所有 MC 基线失效"的初始担忧**不符**——该担忧缺乏证据支持。

#### 3.1.4 CPU/GPU 样本级不一致（README 承诺的直接反例）

| 侧 | 消费方式 | 证据 |
|---|---|---|
| CPU | `normal_simd` 调 `next4()` → **4 个连续计数器**，各取首字做 2 组 Box-Muller | [rng.hpp:112-121](../../include/cpphub/core/rng.hpp#L112-L121) |
| GPU | 每条路径 **1 个计数器** `[path, 0]`，取 `(r0, r1)` 做 1 组 Box-Muller，仅用 `z1` | [gpu_mc.cu:131-144](../../src/performance/gpu/gpu_mc.cu#L131-L144) |

即使两侧算法逐位相同，**消费映射不同** ⇒ 同 seed 下正态流不同 ⇒ `README.md:40` 的 "bit-exact match with CPU" 不成立。

#### 3.1.4b 本地实测（E1-9，决定性）

为避免仅凭推断下结论，用临时探针（`tmp/probe_rng_bit.cu`，未入库）在同一台 RTX 4060 上做逐位比对：两侧使用**完全相同的轮函数**、`counter=path`、`seed=42`、n=20000。

| 比对对象 | 不一致数 | 结论 |
|---|---|---|
| RNG 原始 64 位输出（uniform word） | **0 / 20000** | 整数核 **位精确** |
| `z1 = √(-2 ln u1)·cos(2π u2)`（Box-Muller 正态） | **2870 / 20000（14.35%）** | **正态样本末位不同**（相对差 ~1e-16） |

**判定**：CPU/GPU 的差异**不在 RNG**，而在 `log`/`cos` 的实现——CUDA 设备端数学库与主机 libm 对同一输入不保证逐位相同（NVIDIA 浮点白皮书明确说明"不同数学库不能期望对同一输入算出完全相同的结果"，见 §3.1.7）。因此 `README.md:40` 的承诺**只能在 uniform 层级成立**，正态层级的位精确需要额外的共享变换实现，或把承诺降级为"uniform 流位精确 + 正态容差"。

#### 3.1.5 文档承诺清单（E3，需随裁决一并修正）

| 位置 | 原文要点 |
|---|---|
| [README.md:40](../../README.md#L40) | "CUDA Philox4x64-10 RNG, bit-exact match with CPU `cpphub::core::rng`" |
| [USAGE_GUIDE.md:185,189](../../docs/USAGE_GUIDE.md#L185) | "**Philox 4x64-10**（Salmon et al. 2011）"、"10 轮 Feistel-like 混合" |
| [USAGE_GUIDE.md:196](../../docs/USAGE_GUIDE.md#L196) | "GPU 实现与 CPU 位精确一致（相同 seed + counter 编号 → 相同 Z）" |
| [USAGE_GUIDE.md:694](../../docs/USAGE_GUIDE.md#L694) | GPU MC 表："Philox4x64-10（与 CPU 同算法，位精确一致）" |
| [ARCHITECTURE.md:68](../../docs/architecture/ARCHITECTURE.md#L68) | ADR-004 "计数器 RNG (Philox/Threefry)"，参照 `Random123` |
| [ARCHITECTURE.md:91,331](../../docs/architecture/ARCHITECTURE.md#L91) | `rng.hpp` 提供 `Philox4x64`, `Threefry4x64`（**后者实际不存在**） |

#### 3.1.6 外部锚点（检索级 + E1）

| # | 来源 | 关键事实（原文摘录） | 对本 feature 的含义 |
|---|---|---|---|
| X1 | Random123 `philox.h`（**E1 页面级**） | L207 `PHILOX_M2x64_0=0xD2B74407B1CE6E93`；L210 `PHILOX_M4x64_0=0xD2E7470EE14C6C93`；L213 `PHILOX_M4x64_1=0xCA5A826395121157`；L225 `PHILOX_W64_0=0x9E3779B97F4A7C15`；L247 `PHILOX4x64_DEFAULT_ROUNDS 10` | 现仓库的"M4x64_0"实为 **2x64 乘数**；4x64 的两个乘数与 64 位 Weyl 常量均有唯一权威值 |
| X2 | Random123 README（检索级） | "for CPUs (C and C++) and GPUs (CUDA and OpenCL)"；"implemented entirely in header files"；"passed all of the rigorous SmallCrush, Crush and BigCrush tests"；"**not** suitable for use in cryptography" | 同一套头文件可同时覆盖 CPU/GPU，天然消除"两份实现漂移" |
| X3 | Random123 文档 `ExamplesREADME`（检索级） | `examples/kat_vectors` 含"a few thousand Known Answer Test vectors（method, counter, key, answer）"，`katc.c` **在 serial C/C++/CUDA/OpenCL 间共用** | 可同时验证 CPU 与 GPU 实现 → 两份实现均钉在权威向量上 |
| X4 | TensorFlow 文档（检索级） | `tf.random.Generator.from_seed(123, alg="philox")`；支持 `"philox"` 与 `"threefry"` | 工业级框架把 Philox 作为一等 RNG 选项 |
| X5 | JAX 文档（检索级） | 默认 `threefry2x32`，可选 `philox2x32` | counter-based RNG 是 JAX 可分裂 PRNG 的基础 |
| X6 | NumPy 兼容性政策（检索级） | `Generator`"does not provide a version compatibility guarantee… as better algorithms evolve the bit stream may change"；`RandomState` 则冻结以保位级兼容 | **换流是可接受的工程实践**，前提是：明确声明 + 提供冻结路径 |
| X7 | NVIDIA 浮点白皮书（检索级，§5.1） | "different math libraries cannot be expected to compute exactly the same result for a given input. This applies to GPU programming as well… GPU 用 CUDA math library，CPU 用 host C library" | 正态层级位精确不可能仅靠"同 RNG"达成（已被 E1-9 实测证实） |
| X8 | cuRAND 文档（检索级） | 内置 `CURAND_RNG_PSEUDO_PHILOX4_32_10`；"Each thread computes 4 random numbers… most efficient… multiple of 4 times number of threads" | Philox 的"每线程/每块 4 字"映射是厂商推荐用法，支撑 U-1 的块对齐方案 |
| X9 | Salmon et al. 2011（SC'11，Random123 README 引用） | 提出 counter-based RNG（Philox/Threefry），面向并行与 GPU | 学术锚点：不需要全局状态，counter 决定样本 |

#### 3.1.7 四方案多维度对比

> 维度：**成本**（工作量/回归面）、**收益**（正确性/可验证性）、**兼容**（API/数值/文档）、**维护**（上游/依赖）、**升级潜力**（counter 空间/GPU/扩展）。

| 维度 | **A. 自实现真 Philox4x64-10** | **A2. vendor Random123 + 薄封装** | **C. 自实现真 Philox2x64-10** | **B. 保留算法（改名或仅改文档）** |
|---|---|---|---|---|
| **成本** | 中：`rng.hpp` 轮函数/常量/状态重写（~60 行）+ GPU 副本（~40 行）+ 文档 6 处 + KAT/对齐用例；全量复跑 | 中低：`third_party/random123`（BSD-3，header-only）+ 适配层（CUDA 共用）；文档 + 用例；全量复跑 | 中：改动面小于 A（保留 2 字状态与单字输出语义），但**必须改名** | 低（仅文档）/ 中（改名波及 76 文件） |
| **收益** | 名称/实现/文档一致；可用 Kat 向量；uniform 层位精确可保证 | 同 A，且**零自研轮函数**（权威实现直接可用，含 CUDA 宏） | 同 A（但名需改） | 零正确性收益（改名后仍需自证）；仅"消除文档虚假" |
| **兼容·Python API** | ✅ 不变（`cpphub.Philox4x64` 及方法签名保留） | ✅ 不变（封装同名类） | ❌ 破坏（4x64 名不实 → 需改名为 2x64） | ❌ 改名破坏；仅改文档则保留 |
| **兼容·42 头文件** | ✅ 零改动 | ✅ 零改动（接口层不变） | ✅ 零改动（若保留类名则名实不符） | ⚠️ 改名则全改 |
| **兼容·数值** | 全部模拟样本变化 → 全量复跑 + 少量 seed 重选 | 同 A | 同 A | 零变化 |
| **维护·上游** | 自维护（轮函数固定，风险低但需自证可逆/质量） | **上游活跃**（D. E. Shaw Research），版本可钉 commit；许可证 BSD-3 与项目 Apache-2.0 兼容（需保留版权声明） | 自维护 | 自维护（且无外部参考实现） |
| **维护·缺陷面** | 自研 2 份（CPU/GPU）可能再次漂移 | 单一权威实现 + 共用 KAT → 漂移面最小 | 自研 2 份 | 自研 2 份（已存在漂移风险） |
| **升级潜力·counter** | 2^256（4 字），远超任何 MC 用量 | 同 A | 2^128（2 字），仍足够 | 2^128 |
| **升级潜力·GPU** | 每路径 1 块产出 4 字 → 与 cuRAND 推荐映射一致（X8） | 同 A，且 cuRAND/OpenCL 生态同构 | 每计数器 2 字 | 每计数器 2 字（现为 1 次计算取 2 字） |
| **外部锚点** | X1/X2/X3/X4/X5/X9 | X1/X2/X3/X7/X8 + 许可证 | X1/X3 | 无（自定义算法无任何外部参考） |
| **结论** | ✅ 可选 | ✅ **推荐**（若接受第三方） | ⚠️ 仅在"必须零依赖"时才考虑 | ❌ 不推荐（改名破坏 API，改文档不解决名实问题） |

**推荐（二选一，取决于第三方依赖政策）**：

- 若允许新增 header-only 依赖 → **A2**（vendor Random123）：一致性最强、自研面最小、KAT 现成、CUDA 头文件级支持。
- 若坚持零新增依赖 → **A**：自实现，但**必须**用 Random123 `kat_vectors` 做 KAT（CPU 与 GPU 各跑一遍），以此抵消"自研轮函数"的不确定性。

**无论 A 还是 A2，都必须同时裁定 U-1**（CPU/GPU 样本映射）与 **文档承诺层级**（uniform 位精确 vs 正态位精确），否则 README:40 仍不成立。

---

### 3.2 RISK-032 — 版本号 / README / 覆盖率

#### 3.2.1 版本号真值源（E3，六处）

| # | 位置 | 值 |
|---|---|---|
| 1 | [CMakeLists.txt:2](../../CMakeLists.txt#L2) `project(CppHub VERSION ...)` | **1.0.0** |
| 2 | [config.hpp:5-8](../../include/cpphub/core/config.hpp#L5-L8) `CPPHUB_VERSION_*` 宏 | **1.0.0** |
| 3 | [python/CMakeLists.txt:11](../../python/CMakeLists.txt#L11) | **1.4.0** |
| 4 | [pyproject.toml:24](../../python/pyproject.toml#L24) | **1.4.0** |
| 5 | [cpphub_python.cpp:2043](../../python/src/cpphub_python.cpp#L2043) `m.attr("__version__")` | **"1.4.0"** |
| 6 | git tag（E1-3） | 最新 **v1.8.0**；HEAD = `v1.8.0-16-g0223c92` |

另有 README 叙事停在 v1.4.3（[README.md:7](../../README.md#L7)），C ABI 版本 `CPPHUB_ABI_VERSION 1`（[cpphub_c_api.h:25](../../include/cpphub/c_api/cpphub_c_api.h#L25)，与上述版本**无冲突**，属独立 ABI 契约）。

**问题本质不是"三处不一致"，而是缺少版本语义定义**：`1.0.0`（库 API 版本）与 `1.8.0`（仓库发布 tag）可能**刻意不同**。设计阶段必须先裁定语义（U-2），再选注入机制，否则"单源化"会把两种不同含义强行合并。

#### 3.2.2 README 失真项

| 声称 | 事实 | 证据 |
|---|---|---|
| [README.md:107](../../README.md#L107) "Without CUDA, a CPU stub is compiled automatically and **15 GPU MC tests still run on CPU**" | **双重失真**：(a) 用例数已是 **19**（E1-5）；(b) 无 CUDA 时 `test_gpu_mc` **根本不构建**（`test_gpu_mc: SKIPPED (no CUDA)`），且 stub 的 `price_*` 直接抛异常 | [src/performance/gpu/CMakeLists.txt:93-114](../../src/performance/gpu/CMakeLists.txt#L93-L114)、[gpu_mc_cpu_stub.cpp:66-81](../../src/performance/gpu/gpu_mc_cpu_stub.cpp#L66-L81) |
| [README.md:21](../../README.md#L21) "1412 / 1412 位精确一致" | 当前注册 **2557** 用例；本机构建通过 **2554**（3 个 `NOT_BUILT` 属 C API 关闭下的占位） | E1-7 |
| [README.md:40](../../README.md#L40) bit-exact match | **分层判定**（§3.1.4 / E1-9）：uniform 词层级**成立**（0/20000 不一致）；正态样本层级**不成立**（14.35% 末位不同） | E1-9 |

#### 3.2.3 CI 现状（E3，[ci.yml](../../.github/workflows/ci.yml)）

| job | 矩阵 | 覆盖 |
|---|---|---|
| `build-test` | windows-2022 / ubuntu-latest | 核心库 + 全量 ctest |
| `c-api` | windows-2022 / ubuntu-latest（`-DCPPHUB_ENABLE_C_API=ON`） | `ctest -R CApi` |
| `python` | 2 OS × py3.11/3.12 | `pip install ./python[test,mcp]` + pytest |

**缺口**：无 gcov/lcov/llvm-cov/Codecov/OpenCppCoverage；无覆盖率阈值；**无 CUDA job**（GPU 19 用例在 CI 中零覆盖）；无 macOS/ARM 矩阵（`-march=x86-64-v3` 的 RISK-029 因此也不会被 CI 暴露）。

---

### 3.3 RISK-033 — GoogleTest 供应链

| 项 | 事实 |
|---|---|
| 唯一 FetchContent 依赖 | `FetchContent_Declare(googletest URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz)`，**无 `URL_HASH`**（[CMakeLists.txt:63-66](../../CMakeLists.txt#L63-L66)） |
| 归档 SHA256（E1-1） | `8ad598c73ad796e0d8280b082cebd82a630d73e73cd3c70057938a6501bba5d7`（867764 字节，2026-10-21 下载） |
| tag → commit（E1-2） | `refs/tags/v1.14.0` = `f8d7d77c06936315286eb55f8de22cd23c188571`（`ls-remote` 无 `^{}` 行 ⇒ lightweight tag，该值即 commit） |
| 现有缓解 | CI 用 `actions/cache` 缓存 `build/_deps`（key `gtest-1.14.0-*`），但**首次拉取仍无完整性校验** |
| 本地现状 | `build/_deps/googletest-src` 已缓存；本会话另用 `-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=...` 离线复用 |

**三种钉法对比**：

| 方案 | 强度 | 代价 | 备注 |
|---|---|---|---|
| `URL_HASH SHA256=<E1-1>` | 内容级完整 | 1 行 | **GitHub 自动归档不承诺长期字节不变**（历史上有归档重生成先例）——哈希可能在未来某日失效 |
| `GIT_REPOSITORY + GIT_TAG f8d7d77c…`（40 位 commit） | 身份级稳定 | 改 1 段 + 需 git | 不受归档重生成影响；tag 被重指亦无效 |
| vendor 到 `third_party/` | 最强（离线） | 仓库体积 + 第三方代码入库 | 与 project 现有 `third_party/eigen`（submodule）策略需统一 |

**调研建议**：`GIT_TAG <40 位 commit>` 为主；如坚持 URL 方式则必须同时加 `URL_HASH`。二者取其一即可满足评审要求，但**不要**只钉 tag 名。

---

### 3.4 RISK-034 — `json_min` 外部输入面

#### 3.4.1 现有保护（E3）

| 保护 | 状态 | 证据 |
|---|---|---|
| 递归深度上限 | ✅ `kMaxParseDepth = 64`，对象/数组各设守卫 | [json_min.hpp:27](../../include/cpphub/c_api/detail/json_min.hpp#L27)、[:233](../../include/cpphub/c_api/detail/json_min.hpp#L233)、[:270](../../include/cpphub/c_api/detail/json_min.hpp#L270) |
| 尾随垃圾拒绝 | ✅ `pos_ != s_.size()` → fail | [json_min.hpp:166](../../include/cpphub/c_api/detail/json_min.hpp#L166) |
| 数字语法校验 + 正确舍入 | ✅ 手写语法 + `from_chars` 一次转换 | [json_min.hpp:353-390](../../include/cpphub/c_api/detail/json_min.hpp#L353-L390) |
| UTF-8/代理对校验 | ✅ lone surrogate 拒绝 | [json_min.hpp:326-338](../../include/cpphub/c_api/detail/json_min.hpp#L326-L338) |
| **输入总长度上限** | ❌ 无 | [c_api.cpp:292-331](../../src/c_api/c_api.cpp#L292-L331)：仅判空指针，随后 `parse(args_json, ...)`（隐式 `string_view` ⇒ `strlen`） |
| **数组/对象元素数上限** | ❌ 无（`arr.push_back` / `obj.emplace_back` 无界） | [json_min.hpp:253](../../include/cpphub/c_api/detail/json_min.hpp#L253)、[:283](../../include/cpphub/c_api/detail/json_min.hpp#L283) |
| **节点总数预算** | ❌ 无 | — |
| 语义上限 | ⚠️ 在**解析之后**按字段校验（`kMaxSeriesLength=10000`、`kMaxColumns=64`） | [cpphub_tool_registry.hpp:50-52](../../include/cpphub/c_api/cpphub_tool_registry.hpp#L50-L52) |

#### 3.4.2 放大结构

`Value` 同时内嵌 `std::string`、`std::vector<Value>`、`std::vector<std::pair<std::string, Value>>`，即**每个节点**至少含 3 个容器头 + 1 个 double + tag。

> **估算（未实测，禁止作为结论）**：单节点开销约 10² 字节量级；对 `[0,0,0,…]` 这类最小 token（2 字节/元素）输入，理论上放大 ~几十倍。**已标注为未验证**，若进入设计需实测 `sizeof(Value)`（见 §5.2）。

#### 3.4.3 威胁模型

- 入口为 **C ABI 本地调用**；`cpphub_v1_tool_call` 的调用方是同进程的 Python 绑定 / MCP stdio 服务端，**无网络暴露路径**（仓库内无 socket/HTTP server）。
- 因此风险定性：**本地资源放大（DoS-at-self）**，非远程可利用漏洞；与评审"移交 security-review 复核"的定性一致。
- 但仍应设上限：MCP 客户端（Agent）可产生任意大 JSON，误触发大输入不该让进程 OOM。

**加固候选**（供设计裁决）：(i) 入口 `strlen` 上限（如 8 MiB）+ 明确错误码；(ii) 解析器节点计数预算（如 1e6）；(iii) 数组/对象元素数上限。建议 **(i)+(ii)**：成本低、覆盖"任意放大"，且不改变合法输入语义。

#### 3.4.4 外部锚点（检索级）

| # | 来源 | 关键事实 | 对本 feature 的含义 |
|---|---|---|---|
| Y1 | simdjson 文档 | 单文档上限 **4 GiB**；默认嵌套深度上限 **1024**（`DEFAULT_MAX_DEPTH`） | 主流高性能解析器采用"**尺寸上限 + 深度上限**"双约束 |
| Y2 | Jackson CVE-2026-29062 | `StreamReadConstraints.maxNestingDepth` 默认 **500**；绕过导致 `StackOverflowError` → DoS；3.1.0 修复 | 深度限制是**被真实 CVE 验证过**的必要控制 |
| Y3 | Jackson CVE-2026-50193 | 深层嵌套 JSON 读取 + 序列化 → 资源耗尽 DoS | 放大不仅在解析，也在后续处理；入口应收紧 |
| Y4 | Jackson CVE-2026-18401 | `maxNumberLength`（默认 **1000 字符**）在异步解析器中未生效 → 任意长数字 token 造成内存分配/CPU 耗尽 | **单 token 长度**同样是被 CVE 验证的控制点 |
| Y5 | JSON 安全实践（检索级） | 建议显式深度上限"32 或 64"，高安全场景用迭代式解析器 | 本库 `kMaxParseDepth=64` **已达标**，缺的是尺寸/节点约束 |

#### 3.4.5 三方案多维度对比

| 维度 | **1. 仅入口长度上限** | **2. 长度 + 节点预算**（推荐） | **3. 长度 + 节点 + 元素数** | **4. 不加固（仅登记）** |
|---|---|---|---|---|
| **成本** | 最低（1 常量 + 1 分支 + 1 用例） | 低（长度上限 + 解析器计数 + 2 用例） | 中（多 2 常量/2 错误语义/更多用例） | 零 |
| **收益·防放大** | 覆盖"大输入"（线性放大） | 覆盖**任意**放大（含小体积高扇出） | 同 2，且错误信息可按容器区分 | 无 |
| **收益·对齐外部实践** | 部分（Y1 的尺寸项） | ✅ 完整对齐 Y1 双约束 + Y2 的深度已具备 | 超额 | ❌ |
| **兼容（合法输入）** | 仅在 >8 MiB 时拒绝（合法序列远小于此） | 节点预算 1e6 ≈ 单序列 10k × 64 列仍在界内 | 同 2（元素上限需按工具语义校准） | 完全不变 |
| **维护** | 1 常量 | 2 常量 + 解析器内 1 个计数器 | 4 常量 + 4 错误文案 | 无 |
| **升级潜力** | 未来流式解析时天然保留 | 节点预算可迁移为"最大文档复杂度"指标 | 元素上限与未来 schema 校验重叠 | 技术债累积 |
| **外部锚点** | Y1 | Y1/Y2/Y3/Y5 | Y1~Y5 | 无 |
| **结论** | ⚠️ 防护有缺口（`[0,0,…]` 放大仍在界内但节点数不受控） | ✅ **推荐** | 可接受但冗余（节点预算已覆盖放大） | ❌ |

**补充判定（Y4 对照）**：本库数字解析在 `from_chars` 阶段对超范围数字返回失败（[json_min.hpp:384-387](../../include/cpphub/c_api/detail/json_min.hpp#L384-L387)），**不会**产生 Jackson 式"任意长数字入树"；其扫描成本为 O(输入长度)，已被方案 1/2 的长度上限覆盖。因此无需单独的 `maxNumberLength` 常量。

---

## 4. 综合分析

### 4.1 四项成本/风险矩阵

| 项 | 改动半径 | 数值风险 | 契约风险 | 建议顺序 |
|---|---|---|---|---|
| RISK-033 | 1 行（+CI cache key） | 无 | 无 | **1** |
| RISK-032 | README + 版本注入机制 + CI | 无 | 版本语义未定（U-2） | **2** |
| RISK-034 | 2 文件（c_api.cpp / json_min.hpp）+ 用例 | 无 | 合法输入语义不变（需回归） | **3** |
| RISK-027 | `rng.hpp`（或第三方 vendor）+ GPU 副本 + 6 处文档 + 新用例 | 全量复跑；少量 seed 敏感窗口 | Python API 名保留（A/A2） | **4** |

> **RISK-027 的关键新证据**（E1-9）：现状 RNG 的 uniform 流在 CPU/GPU 间**已经**位精确（0/20000），差异只在 Box-Muller 的正态（14.35% 末位不同）。因此"换 RNG"本身不会修复 README:40 的正态承诺——**必须同时裁定 U-1 的承诺层级**。

### 4.2 交叉影响

- **RISK-032 × RISK-027**：`README.md:40` 与 `USAGE_GUIDE §3.1.3`（RNG；Step 4 勘误：原误标 §3.8）的 RNG 承诺应在同一批改动中订正，避免两次改同一段。
- **RISK-032 × RISK-033**：CI cache key 含 `gtest-1.14.0`；若改为 commit 钉法，cache key 应同步更新。
- **RISK-034 × MCP（S6）**：输入上限属 MCP 对外契约的一部分，宜在 S6 真实客户端联调前落地。
- **RISK-032 覆盖率 × RISK-029**：无 ARM 矩阵，故 `-march=x86-64-v3` 问题在 CI 中不可见；覆盖率/矩阵改造可合并规划（但 RISK-029 属 P2，不在本 feature 范围）。

### 4.3 与既有基线的冲突检查

- 换 RNG **不**使 `tests/fixtures/` 的任何文件失效（§3.1.3）。
- 换 RNG **会**改变所有 C++ 模拟样本，与 `README.md:21` 的 "1412/1412 位精确一致" 无冲突（该数字本就过时且指跨平台一致性，非跨 RNG）。

---

## 5. 幻觉抑制审查（Step 2 Review）

### 5.1 已核实（可重放）

见 §0 的 9 条 E1 命令（含 Random123 一次源行号与 CPU/GPU 位精确探针）+ §3 各表下的 `文件:行`。核心结论均有代码锚点或实测输出，无"凭印象"陈述。

### 5.2 未核实 / 证据边界（诚实声明）

| # | 未核实项 | 影响 | 处置 |
|---|---|---|---|
| 1 | `sizeof(json_min::Value)` 未实测 | §3.4.2 的放大倍数是估算 | DESIGN 前用编译器实测并替换为 E1 |
| 2 | 换 RNG 后具体哪些统计窗口会红 | 只能给出"需复跑"的判断，无法给出数量 | 实施方案中列为 Step 1（先跑基线，记录红项） |
| 3 | `README.md:21` 的 1412 是否为历史真实值 | 本文件只断言"当前不等于 2557"，不断言当时错误 | 保持该措辞 |
| 4 | GitHub 自动归档的长期字节稳定性 | 只是选型考量，未取一次源 | 已在 §3.3 标注"历史上"为弱表述；建议本身不依赖该事实（commit 钉法优先） |
| 5 | NVIDIA 浮点白皮书 §5.1 的**页面级**命中 | 抓取在 §5.1 前被截断，X7 仅为**检索级**摘录 | 结论已由 **E1-9 实测**独立支撑（14.35% 末位不同），不依赖 X7 单独成立 |
| 6 | Random123 `kat_vectors` 的具体内容 | 仅确认其存在与跨 C/C++/CUDA/OpenCL 共用（检索级 X3），未逐条下载 | 实施阶段取一次源并入库为测试夹具 |

### 5.3 待修正项

- 无（本节调研未修改任何源码/文档事实性内容）。

---

## 6. 对设计的输入

### 6.1 关键约束

1. **Python API 不得破坏**（`cpphub.Philox4x64`、`__call__`、`next4`、`discard` 语义保留）。
2. **42 个头文件不得改动**（决定 RNG 变更必须收敛在 `rng.hpp` 内部）。
3. 版本号注入必须**不引入第三方构建依赖**（现状仅 googletest 一个 FetchContent）。
4. `json_min` 保持 header-only、无第三方 JSON 库（spec/mcp-server U-6 定案）。
5. 所有新增保护需附 TDD 红用例（沿用本仓库纪律）。

### 6.2 方案候选（供 DESIGN 选择）

| 项 | 候选 | 备注 |
|---|---|---|
| RISK-027 | A 自实现真 Philox4x64-10 / **A2 vendor Random123 + 薄封装（已裁定）** / C 真 2x64 / B 保留算法 | §3.1.7 多维度对比；A2 需配 U-1 的块对齐与分层承诺 |
| RISK-032 版本注入 | (a) `CMake project(VERSION)` → `configure_file` 生成 `config.hpp` + CI 校验 pyproject/`__version__`；(b) 根目录 `VERSION` 文件单一源，三处读它；(c) 仅 CI 一致性门禁（最小改动） | **U-2 已裁定**：(b) 的"库版本与 tag 分离 + VERSION 文件" |
| RISK-032 覆盖率 | (a) Linux GCC `--coverage` + gcovr + 报告 artifact；(b) 加阈值门禁；(c) Windows OpenCppCoverage | **U-3 已裁定**：(a) 新增独立 job、仅出报告 |
| RISK-033 | GIT_TAG commit（推荐）/ URL_HASH / vendor | §3.3 |
| RISK-034 | 1 仅长度上限 / **2 长度 + 节点预算（已裁定）** / 3 长度 + 节点 + 元素数 / 4 不加固 | §3.4.5 多维度对比 |

### 6.3 风险

| 风险 | 触发 | 缓解 |
|---|---|---|
| 换 RNG 后统计窗口偶发红 | 固定 seed 样本恰好越界 | 允许重选 seed，但需在提交说明中记录"seed 变更，非阈值放宽" |
| 版本语义误判 | 把库版本与仓库 tag 强行合并 | U-2 先行裁决 |
| 覆盖率接入拖慢 CI | 全量 `--coverage` 使编译时间翻倍 | 仅新增独立 job，不与主 job 合并 |
| 输入上限过严 | 合法 10k 序列 + 多列可能接近上限 | 上限值按 `kMaxSeriesLength × kMaxColumns` 上界估算后取 10× 余量 |

### 6.4 未决项登记（四元块）

> **裁决状态（2026-10-21）**：**U-1~U-5 全部已裁定**。以下保留候选与理由，供 DESIGN 引用。

**U-1｜CPU/GPU 随机流的对齐映射与承诺层级 —— ✅ 已裁定**

- 问题：4x64 下每条路径用一个 4 字块（4 uniforms → 2 normals）还是沿用"每输出一计数器"？`README:40` 承诺到哪一层？
- 裁定（Scott，2026-10-21）：**(a) 块对齐 + 分层承诺**——CPU/GPU 每步消耗一个 4 字块；承诺改为"**uniform 流位精确；正态样本容差 1e-15**"。
- 理由：E1-9 已证明"同 RNG 也不保证正态位精确"（14.35% 末位不同），共享超越函数的收益不抵成本。
- DESIGN 须落实：`next4()` 语义变更说明、CPU `normal_simd` 与 GPU 内核的块映射一致、README/USAGE_GUIDE 措辞。

**U-2｜版本语义与单源机制 —— ✅ 已裁定**

- 裁定（Scott，2026-10-21）：**库版本与仓库 tag 分离 + 根目录 `VERSION` 文件为唯一源**。
- 落地要点：`configure_file` 注入 `config.hpp`；CI 校验 `pyproject.toml` / `__version__` 与 `VERSION` 一致；README 说明两者关系。

**U-3｜覆盖率门禁形式 —— ✅ 已裁定**

- 裁定（Scott，2026-10-21）：**新增独立 CI job，仅生成并上传报告，不设阈值**。阈值评估留到基线建立后的下一发布。

**U-4｜`json_min` 加固强度 —— ✅ 已裁定**

- 裁定（Scott，2026-10-21）：**方案 2：入口长度上限 + 解析器节点预算**（§3.4.5）。
- 落地要点：长度上限在 `cpphub_v1_tool_call` 入口（`strlen` 前判）；节点预算在 `Parser` 内计数；错误经现有 `CPPHUB_ERR_INVALID_ARG` 下发；附 TDD 红用例。

**U-5｜RNG 实现路线 —— ✅ 已裁定**

- 裁定（Scott，2026-10-21）：**A2 — vendor Random123 + 薄封装**。
- 落地要点：`third_party/random123`（BSD-3，保留版权声明，与项目 Apache-2.0 兼容）；`rng.hpp` 以 `Philox4x64` 公共类封装 Random123 的 `philox4x64_10`（保持 Python API 与 42 头文件不变）；GPU 侧直接使用同一套 Random123 头（其含 `R123_CUDA_DEVICE` 支持）；KAT 以 `examples/kat_vectors` 为夹具（CPU 与 GPU 各跑一遍）；按 U-1 做块对齐与分层承诺。

---

## 7. 参考文献

**一次源（E1）**

1. Random123 `philox.h` — https://raw.githubusercontent.com/DEShawResearch/random123/main/include/Random123/philox.h （E1-8 行号锚点）
2. Google Test v1.14.0 归档与 tag（E1-1 / E1-2）

**检索级（厂商/社区/学术）**

3. Random123 README — https://github.com/DEShawResearch/random123 （CUDA/OpenCL 支持、header-only、BigCrush、非密码学用途）
4. Random123 *Examples, Tests and Benchmarks* — https://thesalmons.org/john/random123/releases/latest/docs/ExamplesREADME.html （`examples/kat_vectors`、跨语言共用 KAT）
5. Salmon, Moraes, Dror, Shaw (2011), *Parallel Random Numbers: As Easy as 1, 2, 3*, SC'11 — http://dl.acm.org/citation.cfm?doid=2063405
6. TensorFlow `tf.random.Generator` — https://www.tensorflow.org/api_docs/python/tf/random/Generator （`alg="philox"` / `"threefry"`）
7. JAX `jax.random` — https://docs.jax.dev/en/latest/jax.random.html （默认 threefry2x32，可选 philox2x32）
8. NumPy Random Generator 兼容性 — https://numpy.org/doc/stable/reference/random/compatibility.html
9. NVIDIA *Floating Point and IEEE 754*（§5.1 Mathematical Function Accuracy）— https://docs.nvidia.com/cuda/archive/9.1/floating-point/index.html
10. cuRAND 文档 — https://docs.nvidia.com/cuda/curand/group__HOST.html （`CURAND_RNG_PSEUDO_PHILOX4_32_10`）
11. simdjson *The Basics* — https://simdjson.org/api/0.8.0/md_doc_basics.html （4 GiB / depth 1024）
12. NVD CVE-2026-29062 / CVE-2026-50193 / CVE-2026-18401（Jackson 深度与数字长度约束）

**本仓库**

13. [ARCHITECTURE.md](../../docs/architecture/ARCHITECTURE.md) ADR-004；[spec/mcp-server/RESEARCH.md](../mcp-server/RESEARCH.md) U-6（json_min 范围定案）。

---

## 附录 A：证据登记块

```text
E1-1 googletest v1.14.0 tarball SHA256 = 8ad598c73ad796e0d8280b082cebd82a630d73e73cd3c70057938a6501bba5d7 (867764 B)
E1-2 refs/tags/v1.14.0 -> f8d7d77c06936315286eb55f8de22cd23c188571 (lightweight)
E1-3 tags: v1.8.0 / v1.7.1 / v1.7.0 ; HEAD = v1.8.0-16-g0223c92 (2026-10-20)
E1-4 rng.hpp direct include = 76 (42 include/ + 32 tests/ + 2 src/)
E1-5 test_gpu_mc TEST count = 19
E1-6 test_rng TEST count = 9
E1-7 ctest build/ = 2557 registered, 2554 built & passed, 3 NOT_BUILT (C_API=OFF)
E1-8 Random123 philox.h: L207 M2x64_0=0xD2B74407B1CE6E93 ; L210 M4x64_0=0xD2E7470EE14C6C93 ;
      L213 M4x64_1=0xCA5A826395121157 ; L225 W64_0=0x9E3779B97F4A7C15 ;
      L231 W32_0=0x9E3779B9 ; L247 PHILOX4x64_DEFAULT_ROUNDS=10
E1-9 probe(tmp/probe_rng_bit.cu, nvcc 13.2, sm_89, n=20000, seed=42, counter=path):
      uniform-word mismatches = 0 ; Box-Muller z1 bit mismatches = 2870 (14.35%)
```

**调研人**: Scott（鹏） · **复核方式**: 本文件所有 `文件:行` 均已在当前工作树（HEAD `0223c92` + 未提交改动）读取确认；检索级条目均附 URL，且结论不单独依赖检索级证据。
