# 设计文档：P1 整改（RISK-027 / 032 / 033 / 034）

> **Feature**: p1-fixes
> **Spec 步骤**: Step 3（DESIGN）— 上游 [RESEARCH.md](./RESEARCH.md)（Step 1–2，U-1~U-5 已裁定）
> **状态**: **Step 8 文档一致性对齐完成**（2026-10-22，7 处不一致已闭环：1×P2 计数 + 6×P3，见 [CHECKLIST.md](./CHECKLIST.md) §1.2）— 可进入 Step 9（TDD 实施）
> **设计日期**: 2026-10-21

---

## §1 范围

### 1.1 本次范围

| 编号 | 目标 | 裁定依据 |
|---|---|---|
| RISK-033 | GoogleTest 依赖钉到不可变身份（commit SHA） | RESEARCH §3.3 |
| RISK-032 | 版本两轴模型 + README 事实修正 + CI 覆盖率 job + C API 测试过滤修复 | U-2 / U-3 |
| RISK-034 | `json_min` 入口长度上限 + 节点预算 | U-4（方案 2） |
| RISK-027 | vendor Random123 + `Philox4x64` 薄封装 + 块对齐 + KAT | U-5（A2）/ U-1 |

### 1.2 明确不做（避免范围蔓延）

- 不实现 Threefry（`ARCHITECTURE.md:91` 该条为文档错误，按 §5.6 删除）。
- 不换 PCG/xoshiro/mt19937；不引入非 Random123 的 RNG 依赖。
- 不实现"共享超越函数"以追求正态位精确（U-1 否决）。
- 不设覆盖率阈值、不接入 Codecov（U-3）。
- 不改 `-march`/`/arch:AVX2`（RISK-029，P2）；不改 PDE（RISK-030/031，P2）。
- 不重写 README 的 changelog 段落，只改**事实性错误**与新增版本模型说明。

---

## §2 总体方案与依赖顺序

```
R-033 (供应链钉死)          ── 独立，最先（1 行 + CI cache key）
   ↓
R-032 (版本/README/CI)      ── 依赖 R-033 的 cache key 命名；覆盖率 job 独立
   ↓
R-034 (json_min 加固)       ── 独立于 R-027；测试受 C_API 开关约束
   ↓
R-027 (RNG A2)              ── 最大改动；最后执行，需全量复跑与基线记录
```

**关键交叉约束**：R-027 与 R-032 共享 `README.md:40` / `USAGE_GUIDE §3.1.3`（RNG，:182-196）与 `§3.10`（GPU MC，:688-697）的文档段，改动必须合并为一次提交，避免二次修改同一段。（Step 4 勘误 F3：原文误标 §3.8。）

---

## §3 RISK-033 — 供应链钉死

### 3.1 变更

[CMakeLists.txt:61-68](../../CMakeLists.txt#L61-L68) 由 tag 归档 URL 改为 **commit-addressed 归档 URL**：

```cmake
# 供应链 (RISK-033): 以 commit 为不可变身份钉住 googletest v1.14.0
# 不用 tag 名 / tag 归档 URL: tag 可被重指, tag 归档字节不承诺稳定。
# 不用 GIT_REPOSITORY: Step 9 G-05 实测本机 git 协议到 github.com:443 不通,
#   而 codeload HTTPS 可达 (URL 形态同时保证本地与 CI 都能全新 configure)。
FetchContent_Declare(
    googletest
    URL https://codeload.github.com/google/googletest/tar.gz/f8d7d77c06936315286eb55f8de22cd23c188571
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
```

**依据（E1）**：tag→commit 解析 `refs/tags/v1.14.0` = `f8d7d77c…8571`（lightweight tag，无 peeled 行 ⇒ 该值即 commit）；URL 以该 40 位 SHA 为路径，身份不可变。（Step 9 裁决：原 `GIT_REPOSITORY` 形态因本机 git 不通改为 URL，见 §11 D-5。）

### 3.2 CI cache key

`gtest-1.14.0-*` → `gtest-f8d7d77c-*`。key 定义实际共 **2 处**（[ci.yml:59](../../.github/workflows/ci.yml#L59) `build-test`、[ci.yml:110](../../.github/workflows/ci.yml#L110) `c-api`，经 matrix 展开运行期共 4 个 key 实例）；`python` job 只做 pip 缓存、无 `_deps` 缓存，无需改。（Step 4 勘误 F1：原文"三个 job 共 4 处"计数错误。）

### 3.3 保留的逃生口

`-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=<本地路径>` 继续可用（本机离线构建验证已用过）。

### 3.4 备选（不采纳，记录理由）

- **`GIT_REPOSITORY` + `GIT_TAG <sha>`**（原设计）：身份同样不可变，但本机 git 协议到 github.com:443 不通（G-05 实测 exit 128），每次全新 configure 都需 `-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=<本地>` 才能工作；CI 虽可用，但本地开发主路径被阻断。故改为 URL 形态（D-5）。
- **`URL_HASH SHA256=8ad598c7…a5d7`**：哈希绑定"当前归档字节"；若 GitHub 重生成归档，构建会在未来某日失败。commit 身份不受此影响。

---

## §4 RISK-032 — 版本、README、CI

### 4.1 版本模型（两轴，各自单一源）

| 轴 | 单一源 | 注入路径 | 消费者 |
|---|---|---|---|
| **C++ 库 API 版本** | 仓库根 `VERSION`（新文件） | `project(CppHub VERSION …)` 读取；`tools/sync_version.py` 生成 `include/cpphub/core/version.hpp` | `CPPHUB_VERSION{,_MAJOR,_MINOR,_PATCH}` 宏 |
| **Python 发行版本** | `python/pyproject.toml` | scikit-build-core → CMake 变量 `SKBUILD_PROJECT_VERSION` → 编译宏 `CPPHUB_PY_VERSION` | `cpphub.__version__` |
| 仓库发布 tag | git tag | 发布流程校验 `tag == pyproject.version` | README / GitHub Release |

**初始值**：`VERSION` = `1.0.0`（库 API 版本维持现状，符合 U-2 的"分离"语义；本 feature 不改库 API 契约）。`pyproject.version` = **1.9.0**（下一发布号，随 tag 同步；发布号可在实施时由 Scott 调整，机制不变）。

> **为什么不用 `configure_file` 直接写进源码树**：源码树 header-only 使用方式（`-I include`）要求 `config.hpp` 可独立编译；把版本头生成到构建目录会破坏该用法。故采用"`VERSION` 为唯一真值 + 脚本生成 + CI `--check` 拒绝漂移"，等价于单源强制，同时保留源码树可编译性。此为实现层面对 RESEARCH U-2 措辞的收敛（记录于 §11 D-1）。

### 4.2 文件级设计

**A. 新增 `VERSION`（仓库根）**

```
1.0.0
```

**B. 新增 `tools/sync_version.py`**

- `python tools/sync_version.py`：读 `VERSION` → 写 `include/cpphub/core/version.hpp`。
- `python tools/sync_version.py --check`：只校验，不一致退出码 1（CI 用）。
- `version.hpp` **随源码树提交入库**（脚本生成后提交）：源码树 header-only 用法（clone 后直接 `-I include`）没有构建步骤可先生成它，`config.hpp` 的 include 必须在 clone 后即可解析；单一真值仍为 `VERSION`，CI `--check` 保证两者不漂移。（Step 4 补全 F8：入库策略原文未言明。）
- 生成内容（固定格式，含"请勿手改"头注释）：

```cpp
// 本文件由 tools/sync_version.py 从仓库根 VERSION 生成，请勿手改。
#pragma once
#define CPPHUB_VERSION_MAJOR 1
#define CPPHUB_VERSION_MINOR 0
#define CPPHUB_VERSION_PATCH 0
#define CPPHUB_VERSION "1.0.0"
```

**C. 修改 [config.hpp:4-8](../../include/cpphub/core/config.hpp#L4-L8)**

删除手写版本宏，改为：

```cpp
// 版本（单一源：仓库根 VERSION → tools/sync_version.py 生成）
#include "cpphub/core/version.hpp"
```

（保留宏名不变 ⇒ 现有包含方零改动。）

**D. 根 [CMakeLists.txt:2](../../CMakeLists.txt#L2)**

```cmake
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" CPPHUB_LIB_VERSION LIMIT_COUNT 1)
project(CppHub VERSION ${CPPHUB_LIB_VERSION} LANGUAGES CXX)
```

**E. [python/CMakeLists.txt:10-11](../../python/CMakeLists.txt#L10-L11)**

```cmake
# 版本单一源: wheel 构建时由 scikit-build-core 注入 SKBUILD_PROJECT_VERSION;
# 直接构建 python/ 时回退到仓库根 VERSION (库版本), 保证无硬编码。
if(DEFINED SKBUILD_PROJECT_VERSION)
    set(CPPHUB_PY_VERSION "${SKBUILD_PROJECT_VERSION}")
else()
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../VERSION" CPPHUB_PY_VERSION LIMIT_COUNT 1)
endif()
project(cpphub_python VERSION ${CPPHUB_PY_VERSION} LANGUAGES CXX)
```

并新增：

```cmake
target_compile_definitions(_core PRIVATE CPPHUB_PY_VERSION="${CPPHUB_PY_VERSION}")
```

**F. [cpphub_python.cpp:2043](../../python/src/cpphub_python.cpp#L2043)**

```cpp
m.attr("__version__") = CPPHUB_PY_VERSION;  // 由构建系统注入, 不再硬编码
```

**G. [python/pyproject.toml:24](../../python/pyproject.toml#L24)**

`version = "1.9.0"`（与发布 tag 同步；见 §4.2 初始值说明）。

### 4.3 README 事实修正

| 位置 | 现状 | 改为 |
|---|---|---|
| [README.md:40](../../README.md#L40) | "CUDA Philox4x64-10 RNG, bit-exact match with CPU" | "Philox4x64-10 (Random123): **uniform stream bit-exact** with CPU; normal samples agree to 1e-15 (CUDA/host math libraries differ in last bits)" |
| [README.md:107](../../README.md#L107) | "Without CUDA, a CPU stub is compiled automatically and 15 GPU MC tests still run on CPU." | "Without CUDA, a CPU stub is compiled and `test_gpu_mc` is **skipped** (`test_gpu_mc: SKIPPED (no CUDA)`); GPU entry points throw `std::runtime_error`. GPU tests require `-DCPPHUB_ENABLE_CUDA=ON`." |
| [README.md:21](../../README.md#L21) | "1412 / 1412 位精确一致" | 删除冻结数字，改为"全量测试通过（见 CI）；跨平台位精确指**同编译器族内**可复现，非跨 RNG/跨数学库" |
| [README.md:121](../../README.md#L121) | "Cross-platform float determinism \| bit-exact" | 加注"（同 RNG 版本 + 同编译选项）" |
| README 顶部（§版本） | 无版本模型说明 | 新增 3 行：C++ API 版本见 `VERSION`；Python 发行版本见 `pyproject.toml`；仓库 tag 独立 |

### 4.4 覆盖率 CI job（U-3：仅报告，无阈值）

新增独立 job（不改现有三个 job 的语义）：

```yaml
  coverage:
    name: Coverage (ubuntu-latest, GCC)
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
        with: { submodules: recursive }
      - name: Install toolchain
        run: sudo apt-get update && sudo apt-get install -y ninja-build gcc g++ && python3 -m pip install gcovr
      - name: Configure (coverage)
        run: |
          cmake -S . -B build_cov -G Ninja \
            -DCMAKE_BUILD_TYPE=Debug \
            -DCMAKE_CXX_FLAGS="--coverage -O0 -g" \
            -DCMAKE_EXE_LINKER_FLAGS="--coverage"
      - name: Build
        run: cmake --build build_cov --parallel 4
      - name: Test
        working-directory: build_cov
        run: ctest --output-on-failure --parallel 4
      - name: gcovr report
        run: |
          gcovr --root . --filter 'include/cpphub/' --exclude 'third_party' \
                --xml coverage.xml --html-details coverage.html
      - uses: actions/upload-artifact@v4
        with:
          name: coverage-report
          path: |
            coverage.xml
            coverage*.html
```

**理由与边界**：
- 另在 `build-test` job（Linux only）加版本漂移门禁步骤 `python3 tools/sync_version.py --check`（Step 8 对齐：§4.2 B 的要求在 CI 中的落地位置；实现细节见 [IMPLEMENTATION.md](./IMPLEMENTATION.md) §3.6.2）。
- `Debug` ⇒ 根 CMake 的 Release 架构硬编码（`-march`）不生效，避免 RISK-029 干扰。
- 不设阈值 ⇒ 不会因历史模块低覆盖而误伤合入（U-3）。
- 报告只过滤 `include/cpphub/`，排除测试与第三方。

### 4.5 修复 C API 测试过滤（既有缺口）

**事实（E3）**：[ci.yml:125](../../.github/workflows/ci.yml#L125) 用 `ctest -R "CApi"`。该正则只匹配 `CApi.*` 与 `CApiTools.*`；`JsonMin*`（14 用例）与 `ToolRegistry.*`（14 用例）**在 CI 中从未运行**。（Step 4 勘误 F2：ToolRegistry 实测 14，原文误记 15。）

**改为**（在 C API 配置下运行该配置构建出的全部测试）：

```yaml
      - name: Run C API Tests
        working-directory: build
        run: ctest --build-config Release --output-on-failure --parallel 4
```

（C API 配置构建的测试集 = 全量 + c_api 增量，与 `build-test` 重叠但可接受；如需最小化，用 `-R "CApi|JsonMin|ToolRegistry"`。设计采用**不带过滤**，避免下次再漏。）

---

## §5 RISK-034 — `json_min` 加固

### 5.1 常量与位置

| 常量 | 值 | 位置 | 理由 |
|---|---|---|---|
| `kMaxArgsJsonBytes` | `8 * 1024 * 1024`（8 MiB） | [cpphub_tool_registry.hpp](../../include/cpphub/c_api/cpphub_tool_registry.hpp#L49-L52) 现有语义上限旁 | 属**工具面输入策略**，与 `kMaxSeriesLength` / `kMaxColumns` 同处 |
| `kMaxParseNodes` | `1'000'000` | [json_min.hpp:27](../../include/cpphub/c_api/detail/json_min.hpp#L27) `kMaxParseDepth` 旁 | 属**解析器实现约束**，与深度上限同层 |

**上界校验**：合法最大输入 ≈ `kMaxSeriesLength(10k) × kMaxColumns(64)` ≈ 640k 个数值节点；`1e6` 留 ~1.5× 余量。若未来放宽序列上限，两者需同步。

### 5.2 入口长度检查（[c_api.cpp:292-331](../../src/c_api/c_api.cpp#L292-L331)）

在 `tool_name`/`args_json` 判空之后、`find_tool` 之前插入：

```cpp
const std::size_t args_len = std::strlen(args_json);
if (args_len > cpphub::v1::c_api::kMaxArgsJsonBytes) {
    set_error("args_json exceeds limit (" + std::to_string(args_len) +
              " > " + std::to_string(cpphub::v1::c_api::kMaxArgsJsonBytes) + " bytes)");
    return CPPHUB_ERR_INVALID_ARG;
}
```

（`<cstring>` 已在 [c_api.cpp:18](../../src/c_api/c_api.cpp#L18) 包含，无需新增；`strlen` 前提是 C ABI 契约的 NUL 结尾字符串。）

### 5.3 节点预算（`json_min::detail::Parser`）

```cpp
// 成员
std::size_t nodes_ = 0;

// parse_value(out, depth) 入口
if (++nodes_ > kMaxParseNodes) {
    return fail("maximum node count exceeded");
}
```

**计数语义**：每调用一次 `parse_value` = 1 个节点（null/bool/number/string/array/object 各计 1）；对象键（key）不计入（其规模已被长度上限约束，且键不是 JSON 值）。该语义写入常量注释，测试按此断言。

### 5.4 错误映射

两类失败都走既有路径：`json_min::parse` 返回 false → `cpphub_v1_tool_call` 返回 `CPPHUB_ERR_INVALID_ARG`，错误文案经 `cpphub_v1_get_last_error` 下发（[c_api.cpp:321-326](../../src/c_api/c_api.cpp#L321-L326)）。**不新增错误码**（保持 C ABI 版本 1 的纯增/不变约定）。

### 5.5 测试

| 用例 | 文件 | 断言 |
|---|---|---|
| `JsonMinBudget.RejectsWideArrayOverBudget` | `test_json_min.cpp` | `[0,0,…]`（`kMaxParseNodes+2` 个元素）→ `parse == false`，error 含 "node" |
| `JsonMinBudget.AcceptsAtBoundary` | 同上 | 恰 `kMaxParseNodes` 个节点 → `parse == true` |
| `JsonMinBudget.DeepNestingStillLimitedByDepth` | 同上 | 回归：深度超限仍报 "nesting depth"（两类限制互不遮蔽） |
| `CApiTools.RejectsOversizedArgsJson` | `test_c_api_tools.cpp` | 构造 `kMaxArgsJsonBytes+1` 字节字符串 → `CPPHUB_ERR_INVALID_ARG`，`get_last_error` 含 "exceeds limit" |
| `CApiTools.OversizedCheckPrecedesParse` | 同上 | 超长但语法非法 → 报 "exceeds limit" 而非 JSON 语法错误（证明检查前置） |

**注意**：这些测试位于 `if(CPPHUB_ENABLE_C_API)` 块内，本地需在 `build_capi` 验证；CI 由 §4.5 修复后的 `c-api` job 覆盖。

---

## §6 RISK-027 — RNG（A2：vendor Random123）

### 6.1 vendor 布局

```
third_party/random123/
├── include/Random123/        # 上游头文件原样拷贝 (philox.h, array.h, features/*)
├── LICENSE                   # BSD-3-Clause 原文
└── UPSTREAM.md               # 来源 URL + 固定 commit + 拷贝日期
```

- 拷贝整个上游 `include/Random123/`（约十余个头文件），**不做修改**；`UPSTREAM.md` 记录固定 commit（实施时取 `main` 的具体 SHA 并写入）。
- 许可证：BSD-3-Clause，与项目 Apache-2.0 兼容；保留版权声明（`LICENSE` + 头文件自带 notice）。
- 不使用 git submodule（避免下游 wheel 构建追加 `submodules: recursive` 负担；本项目 Eigen 已是 submodule，不再增加）。

### 6.2 包含路径

| 目标 | 新增 |
|---|---|
| 根 [CMakeLists.txt:25-28](../../CMakeLists.txt#L25-L28) `cpphub_core` | `${CMAKE_CURRENT_SOURCE_DIR}/third_party/random123/include` |
| [src/performance/gpu/CMakeLists.txt:36-39](../../src/performance/gpu/CMakeLists.txt#L36-L39) `cpphub_cuda` | `${GPU_THIRD_PARTY_DIR}/random123/include` |
| [python/CMakeLists.txt:27-31](../../python/CMakeLists.txt#L27-L31) `cpphub_headers` | `${CMAKE_CURRENT_SOURCE_DIR}/../third_party/random123/include` |

### 6.3 `rng.hpp` 公共契约（类名与方法签名不变）

```cpp
class Philox4x64 {
    r123::Philox4x64::ctr_type ctr_{};    // 4×u64 计数器
    r123::Philox4x64::key_type key_{};    // 2×u64 密钥
    r123::Philox4x64::ctr_type block_{};  // 当前块 4 字
    int idx_ = 4;                          // 4 = 块已耗尽
public:
    explicit Philox4x64(uint64_t seed = 0, uint64_t stream = 0);
    uint64_t operator()();
    void discard(uint64_t n);
    std::array<uint64_t, 4> next4();
    // min()/max()/result_type 保留
};
```

| 方法 | 新语义 | 与旧实现差异（行为变更，须写文档） |
|---|---|---|
| 构造 | `ctr={0,0,0,0}`, `key={seed, stream}` | 同 |
| `operator()` | 返回**流顺序**的下一字；块耗尽时 `ctr++` 并重算整块 | 旧：每次 `ctr++` 后取首字（等价于"每块只取 1 字"） |
| `next4()` | 返回流顺序的 4 个字（实现为 4× `operator()`） | 旧：4 个计数器的首字；**输出值改变** |
| `discard(n)` | 丢弃 n 个**输出字**：`ctr += n/4`，再消费 `n%4` | 旧：`ctr += n`（当时 1 counter = 1 输出，等价） |

- `normal_simd` / `generate_correlated` **代码不变**（它们只调用 `next4()` / `operator()`）。
- **块对齐不变量**：`normal_simd` 每步恰消耗 4 字 = 1 块 ⇒ 若 `normal_simd` 从块边界开始，则每步保持块边界；GPU 每路径恰用 1 块 ⇒ **CPU 第 p 步与 GPU 第 p 条路径使用同一 counter 块与同一 4 字**。

### 6.4 GPU 映射（`gpu_mc.cu`）

```cpp
#if defined(__CUDA_ARCH__)
// Windows+nvcc: 关闭 MSVC 内建分派, 让 device 走 __umul64hi (D-6 / Step 9 探针)
#define R123_USE_MULHILO64_MSVC_INTRIN 0
#endif
#include <Random123/philox.h>

// 内核内（每路径 1 块）
r123::Philox4x64 rng;
r123::Philox4x64::ctr_type ctr{};
ctr.v[0] = static_cast<uint64_t>(path);           // 其余 3 字为 0
r123::Philox4x64::key_type key{};
key.v[0] = seed; key.v[1] = 0;
const auto r = rng(ctr, key);                      // 4 字
double u1 = uint64_to_unit(r.v[0]);
double u2 = uint64_to_unit(r.v[1]);
double z1, z2; gpu_box_muller(u1, u2, z1, z2);     // 仍只用 z1（保持单路径单正态）
```

删除设备端自写的 `d_PHILOX_M4x64_*` 常量（[gpu_mc.cu:42-45](../../src/performance/gpu/gpu_mc.cu#L42-L45)）与 `philox4x64_10()`（[:59-77](../../src/performance/gpu/gpu_mc.cu#L59-L77)）。

### 6.5 KAT（known-answer test）

| 项 | 设计 |
|---|---|
| 夹具 | `tests/fixtures/rng/kat_philox4x64_10.hpp`：**生成头文件**（Step 5 勘误 E-1 / §11 D-4，弃 .txt+运行时解析），从上游 `examples/kat_vectors` 抽取 `philox4x64`、10 轮的行（每行 = counter[4] key[2] answer[4]）；沿用 grunfeld 先例（生成头 + `target_include_directories`），规避 ctest 工作目录差异 |
| CPU 用例 | `RngKat.Philox4x64MatchesKnownAnswers`：逐行构造 ctr/key，直接调 `r123::Philox4x64`，比对 4 字 |
| GPU 用例 | `GpuMC.PhiloxMatchesKnownAnswers`：把向量拷到设备，内核逐行计算，比对 4 字（证明设备端编译出同一实现） |
| CPU/GPU 对齐用例 | `GpuMC.UniformStreamMatchesCpu`：对 p=0..9999，比较 CPU `next4()` 与 GPU 块输出**逐字相同**（把 E1-9 固化为回归） |
| 正态层级用例 | `GpuMC.NormalSamplesWithinTolerance`：同一 u1/u2 下 CPU/GPU `z1` 相对差 ≤ 1e-15（固化 U-1 的分层承诺） |

> **夹具格式**：上游 `kat_vectors` 的列序在本阶段未能页面级确认（下载被拒），列为实施期确认项 U-6；**列序解析移入一次性抽取脚本**（tmp/，不入库），夹具头文件只含编译期常量，不预写死。

### 6.6 文档修正（与 §4.3 合并提交）

| 位置 | 改为 |
|---|---|
| [USAGE_GUIDE.md:185-196](../../docs/USAGE_GUIDE.md#L185-L196) | 算法：Random123 `philox4x64_10`（Salmon et al. 2011）；counter 2^256；输出 4×u64/块；`next4()` = 4 个字；分层承诺（uniform 位精确 / 正态 1e-15） |
| [USAGE_GUIDE.md:693](../../docs/USAGE_GUIDE.md#L693) | GPU 行同步为同一表述 |
| [ARCHITECTURE.md:91](../../docs/architecture/ARCHITECTURE.md#L91) | 删除不存在的 `Threefry4x64`；写 `Philox4x64`（Random123 vendor） |
| [ARCHITECTURE.md:68](../../docs/architecture/ARCHITECTURE.md#L68) | ADR-004 备注 `Random123` 已 vendor 至 `third_party/random123` |

### 6.7 数值影响与复跑策略

1. **夹具不受影响**：`tests/fixtures/**` 由 Python/R 生成（RESEARCH §3.1.3），不重算。
2. **全量复跑**：`ctest` 全量，记录失败清单。
3. **失败处置**：仅允许两类修复——(a) 若断言是"种子上偶然越界"的统计窗口，**重选 seed** 并在提交说明记录"seed 变更，非阈值放宽"；(b) 若是实现错误，修实现。**禁止**放宽统计容差。
4. **性能记录**：`benchmark_philox` 改造前后各跑一次，把吞吐写入 DEVELOPMENT_LOG（用户要求性能改动附前后数据）。不设硬门禁，但若下降超过 2× 需在日志中说明可接受性。

### 6.8 兼容性

| 面 | 影响 |
|---|---|
| C++ 公共 API | 无（类名、构造函数、方法签名不变） |
| 42 个 include 方 | 无（仅 `rng.hpp` 内部实现变化） |
| Python API | 无签名变化；`next4()` **数值与语义变化**（4 字 = 1 块）；docstring 更新 |
| 数值流 | 全部改变（可接受，RESEARCH X6：NumPy 同样"算法演进改变比特流"） |
| 轮次/常量 | 由 `PHILOX4x64_DEFAULT_ROUNDS=10` 与上游常量接管 |

---

## §7 测试计划（TDD 红→绿顺序）

| 阶段 | 红用例 | 预期红 | 实现 |
|---|---|---|---|
| P1 | `CApiTools.RejectsOversizedArgsJson` / `OversizedCheckPrecedesParse`、`JsonMinBudget.*` | 全部失败（无上限） | §5.2 + §5.3 |
| P2 | `RngKat.Philox4x64MatchesKnownAnswers`、`GpuMC.UniformStreamMatchesCpu`、`GpuMC.PhiloxMatchesKnownAnswers` | 失败（KAT 无、uniform 对齐未验证） | §6.3 + §6.4 |
| P3 | `GpuMC.NormalSamplesWithinTolerance` | **非红用例（Step 5 勘误 E-2）**：实测末位差 ≈1e-16 相对差，低于 1e-15 容差；且实施前 GPU 探针不存在、无法链接。本用例是**回归锁**而非红驱动 | 无代码改动；实现块对齐后以 1e-15 容差通过（U-1 契约） |
| P4 | 版本一致性：`tools/sync_version.py --check`（CI 脚本，非 gtest） | 初始失败（无 VERSION） | §4.2 |
| P5 | 全量回归 | 记录并处置 | §6.7 |

**新增/修改文件预估**：

| 类型 | 文件 |
|---|---|
| 新增 | `VERSION`、`tools/sync_version.py`、`include/cpphub/core/version.hpp`、`third_party/random123/**`、`tests/fixtures/rng/kat_philox4x64_10.hpp`（生成头，D-4）、`include/cpphub/performance/gpu/rng_device_probe.hpp`（GPU 测试探针声明，Step 5 细化） |
| 修改 | `CMakeLists.txt`、`python/CMakeLists.txt`、`python/pyproject.toml`、`python/src/cpphub_python.cpp`、`include/cpphub/core/config.hpp`、`include/cpphub/core/rng.hpp`、`src/performance/gpu/gpu_mc.cu`、`src/performance/gpu/CMakeLists.txt`、`src/c_api/c_api.cpp`、`include/cpphub/c_api/cpphub_tool_registry.hpp`、`include/cpphub/c_api/detail/json_min.hpp`、`.github/workflows/ci.yml`、`README.md`、`docs/USAGE_GUIDE.md`、`docs/architecture/ARCHITECTURE.md`、`tests/unit/core/test_rng.cpp`、`tests/unit/performance/test_gpu_mc.cpp`、`tests/unit/c_api/test_json_min.cpp`、`tests/unit/c_api/test_c_api_tools.cpp` |

---

## §8 风险与回滚

| 风险 | 概率 | 影响 | 缓解 / 回滚 |
|---|---|---|---|
| Random123 头在 nvcc 下编译失败（宏/特性检测） | 中 | GPU 构建断 | 先做 1 个最小 .cu 探针验证；失败则退回方案 A（自实现），不影响其余 P1 项 |
| KAT 夹具格式与预期不符 | 中 | 测试无法解析 | U-6：实施期先下载查看首行再写解析器；必要时手工抽取 20 行 |
| 换 RNG 导致部分统计窗口红 | 高 | 需重选 seed | §6.7 处置规则（禁止放宽容差） |
| `__version__` 在非 scikit-build 路径缺失宏 | 低 | 编译失败 | python/CMakeLists 的回退分支提供 `CPPHUB_PY_VERSION` |
| 覆盖率 job 时长过长 | 中 | CI 排队 | 独立 job；`--parallel 4`；只跑一次全量 |
| 节点预算误伤合法输入 | 低 | 工具调用被拒 | 上界 1e6 vs 合法最大 ~6.4e5；测试覆盖边界 |

**回滚单元**：四项各自独立提交；R-027 的 vendor 目录与 `rng.hpp` 在同一提交，可整体 revert。

---

## §9 验收标准（供 CHECKLIST 使用）

| # | 标准 | 验证方式 |
|---|---|---|
| A1 | 无 `URL` tag 依赖；`FETCHCONTENT` 钉 commit | 读 `CMakeLists.txt` + 清缓存 configure |
| A2 | `VERSION`、`version.hpp`、`__version__`、`pyproject` 三者关系符合两轴模型 | `tools/sync_version.py --check` + pytest `__version__` 断言 |
| A3 | README 三处事实错误已修，且无冻结的测试总数 | grep `15 GPU MC` / `1412` 零命中 |
| A4 | 覆盖率 job 产出 artifact | CI 运行后下载 `coverage-report` |
| A5 | C API 测试在 CI 中实际运行 | `c-api` job 日志含 `JsonMin*` / `ToolRegistry.*` |
| A6 | 8 MiB+1 输入被拒；节点预算边界正确 | P1 用例全绿 |
| A7 | KAT（CPU/GPU）全绿；uniform 流 CPU/GPU 逐字相同 | P2 用例全绿 |
| A8 | 正态层级 CPU/GPU 相对差 ≤ 1e-15 | P3 用例全绿 |
| A9 | 全量回归通过（允许 §6.7 的 seed 重选） | `ctest` 全绿 + 失败处置记录 |
| A10 | 性能前后数据已记录 | `benchmark_philox` 两次输出 |

---

## §10 未决/实施期确认项

| 编号 | 项 | 处置 |
|---|---|---|
| U-6 | 上游 `kat_vectors` 的确切列序 | 实施第一步下载并查看首行；不预写死 |
| U-7 | `pyproject.version` 用 1.9.0 还是其它发布号 | 由 Scott 在实施时确认；机制不受影响 |
| U-8 | vendor Random123 的固定 commit | 实施时取 `main` 当前 SHA 写入 `UPSTREAM.md` |
| U-9 | 覆盖率是否覆盖测试代码本身 | 本设计不覆盖（只测 `include/cpphub/`）；如需再议 |
| U-10 | Random123 的调用形态与 64 位开关 | 设计按 `rng(ctr, key)` + `ctr_type::v[]` 书写；"上游 `philox.h` 第 175 行为 `#if R123_USE_PHILOX_64BIT`" **已于 Step 4 页面级复核为真**（2026-10-21，raw philox.h）。实施第一步仍需：用最小探针验证调用形态，并确认 `R123_USE_PHILOX_64BIT` 在 MSVC/nvcc 下默认开启。若形态不同，仅调整封装内部，公共契约不变 |

---

## §11 与调研裁决的偏差记录

| # | 调研表述 | 设计收敛 | 理由 |
|---|---|---|---|
| D-1 | U-2 落地要点写"`configure_file` 注入 `config.hpp`" | 改为 `VERSION` → `tools/sync_version.py` → **源码树内** `version.hpp` + CI `--check` | 纯 `configure_file` 到构建目录会破坏源码树 header-only 用法；语义（单一真值 + 禁止漂移）不变 |
| D-2 | A2 表述"薄封装" | 明确为"保持类名/方法签名的适配层；`next4()` 语义变更需文档化" | 4x64 块语义与旧 2 字实现不可避免存在输出语义差 |
| D-3 | 覆盖率"新增 job 仅出报告" | 额外修复 `ctest -R "CApi"` 过滤缺口 | 该缺口使 JsonMin/ToolRegistry 测试从未在 CI 运行，属同一"CI 可信度"目标 |
| D-4 | KAT 夹具 ".txt + 运行时解析器"（Step 5 勘误 E-1） | 改**生成头文件** `kat_philox4x64_10.hpp`（grunfeld 先例：生成头 + `target_include_directories`），列序解析（U-6）移入一次性抽取脚本（tmp/，不入库） | ctest 工作目录差异已在此前用生成头规避过一次（grunfeld）；GPU 用例（test_gpu_mc）同样直接 include，运行时解析器需在两侧各实现一遍 |
| D-5 | R-033 "`GIT_REPOSITORY` + `GIT_TAG <sha>`"（Step 9 裁决，Scott 2026-10-22） | 改 **`URL https://codeload.github.com/google/googletest/tar.gz/<sha>`**（+ `DOWNLOAD_EXTRACT_TIMESTAMP TRUE`） | 本机 git 协议到 github.com:443 不通（实测 exit 128），codeload HTTPS 可达；URL 仍以 commit SHA 为身份，不可变性语义不变，且本地/CI 均可全新 configure |
| D-6 | R-027 GPU include（Step 9 探针发现） | `gpu_mc.cu` 在 `#include <Random123/philox.h>` **之前**加 device pass 宏守卫：`#if defined(__CUDA_ARCH__)` → `#define R123_USE_MULHILO64_MSVC_INTRIN 0` → `#endif` | Windows+nvcc 下 `MULHILO64_MSVC_INTRIN` 优先于 `CUDA_INTRIN`（philox.h L199-202 分派顺序），device 代码解析到 host 内建 `_umul128` → 编译失败；守卫后 device 走 `__umul64hi`，已在 sm_89 上验证输出与官方 KAT 逐字一致 |

---

## §12 Step 4 设计 Review 记录（2026-10-21）

**方法**：DESIGN 全部 `文件:行` 引用对当前工作树逐条重读；外部事实按 E1 重放（GitHub API tag→commit、raw philox.h 行号锚点）；与 RESEARCH §6.4 裁决（U-1~U-5）逐项对齐。

**核对通过**（关键项）：

- §3.1 SHA：GitHub API 页面级 `refs/tags/v1.14.0` → `f8d7d77c06936315286eb55f8de22cd23c188571`，`type: "commit"`（lightweight，即 commit 本体）。`git ls-remote` 本机被网络阻断（exit 128），以 API 复核替代。
- E1-8 六个常量行号（L207/210/213/225/231/247）与 **L175 `#if R123_USE_PHILOX_64BIT`** 全部复核为真（raw philox.h）。
- §4.2 各行号：config.hpp:4-8 / CMakeLists.txt:2、:15-21（Release `-march` 守卫）、:25-29 / python/CMakeLists.txt:10-11、:27-31 / cpphub_python.cpp:2043 / pyproject.toml:24。
- §4.3 README:21/:40/:107/:121 现状引用逐字属实；"SKIPPED (no CUDA)" 与 gpu CMakeLists:114 实际文案一致。
- §4.5：ci.yml:125 过滤缺口属实；c_api 4 个测试目标确在 `if(CPPHUB_ENABLE_C_API)` 内（tests/CMakeLists.txt:203-222）。
- §5 全部锚点：registry:49-52 / json_min.hpp:27、:192（`parse_value(Value&, int)` 签名）/ c_api.cpp:292-331、:321-326 / `CPPHUB_ERR_INVALID_ARG`（cpphub_c_api.h:30）。
- §6.2 三处 include 路径行号；GPU stub 分支不编 `gpu_mc.cu`，无需 random123 路径。
- §6.3 旧三方法语义描述与 rng.hpp 实际实现一致；§6.4 `gpu_box_muller` 函数名存在；§6.6 ARCHITECTURE:68/:91（:91 确含 `Threefry4x64` 假条目）；§6.7 `benchmark_philox.cpp` 存在；Eigen 为 submodule（.gitmodules）。

**勘误（已回写本文档）**：

| # | 位置 | 原文 | 事实 | 证据 |
|---|---|---|---|---|
| F1 | §3.2 | "三个 job 共 4 处：build-test 2、c-api…" | key 定义 **2 处**（build-test/c-api 各 1，运行期 4 实例）；python 无 gtest 缓存 | [ci.yml:59](../../.github/workflows/ci.yml#L59) / [:110](../../.github/workflows/ci.yml#L110) |
| F2 | §4.5 | "ToolRegistry.*（15 用例）" | 实为 **14** | `^(TEST\|TEST_F)\(ToolRegistry` 全仓计数 = 14 |
| F3 | §2 | "USAGE_GUIDE §3.8" | RNG 段为 **§3.1.3**（:183）；GPU 表为 §3.10（:688）；§3.8 是波动率曲面 | USAGE_GUIDE 标题行索引 |
| F4 | §6.6 | "USAGE_GUIDE.md:694" | RNG 行为 **:693**（694 是定价行） | USAGE_GUIDE :691-697 |
| F5 | §5.2 | "需新增 `#include <cstring>`" | c_api.cpp:18 **已包含** | [c_api.cpp:18](../../src/c_api/c_api.cpp#L18) |
| F6 | §6.4 | "gpu_mc.cu:40-75" | 常量 :42-45、`philox4x64_10()` :59-77 | [gpu_mc.cu:42-77](../../src/performance/gpu/gpu_mc.cu#L42-L77) |
| F7 | §10 U-10 | "philox.h 第 175 行"原未核实 | **已核实为真** | raw philox.h L175（本次页面级） |
| F8 | §4.2 B | version.hpp 入库策略未言明 | 补充"生成物随源码树提交 + CI `--check` 防漂移"（设计完整性补全，非事实错误） | 源码树 header-only 用法推演 |

**连带勘误**：RESEARCH §4.2 同源的 "USAGE_GUIDE §3.8" 标注一并修正为 §3.1.3。

**结论**：无方向性/裁决性错误；8 处勘误均为行号、计数、文案级别，不改变任何设计决策；U-1~U-5 裁决映射完整（§1.1）。**可进入 Step 5（IMPLEMENTATION）**。

---

**设计人**: Scott（鹏） · **复核方式**: 所有 `文件:行` 均对当前工作树读取确认（Step 4 逐条重放，见 §12）；外部依赖版本以 E1 证据（RESEARCH §0 + Step 4 复核）为准。
