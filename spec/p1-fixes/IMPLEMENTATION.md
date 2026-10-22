# 实施文档：P1 整改（RISK-027 / 032 / 033 / 034）

> **Feature**: p1-fixes
> **Spec 步骤**: Step 5（IMPLEMENTATION）— 上游 [DESIGN.md](./DESIGN.md)（Step 4 review 完成，§12）+ [RESEARCH.md](./RESEARCH.md)（U-1~U-5 已裁定）
> **状态**: **Step 6 review 完成**（2026-10-21，13 处勘误已回写，见 §10）— 可进入 Step 7（CHECKLIST）
> **实施环境**: Windows 11 + MSVC 2022 + RTX 4060（SM 8.9）+ CUDA 12.x + Python 3.12

---

## §0 实施环境与基线

### 0.1 三构建目录矩阵

| 目录 | 配置 | 用途 |
|---|---|---|
| `build/` | 默认（无 C_API / 无 CUDA） | 核心回归（R-027 主战场） |
| `build_capi/` | `-DCPPHUB_ENABLE_C_API=ON` | R-034 测试（test_json_min / test_c_api_tools 均在 `if(CPPHUB_ENABLE_C_API)` 块内，[tests/CMakeLists.txt:203-222](../../tests/CMakeLists.txt#L203-L222)） |
| `build_cuda/` | `-DCPPHUB_ENABLE_CUDA=ON`（本机自动 `CMAKE_CUDA_ARCHITECTURES=89`，[CMakeLists.txt:88-90](../../CMakeLists.txt#L88-L90)） | R-027 GPU 侧（gpu_mc.cu / test_gpu_mc / 探针） |

### 0.2 离线 gtest 复用

本机网络对 GitHub 受限（Step 4 实测 `git ls-remote` exit 128）。`FetchContent` 改 commit 后首次 configure 需拉取，若网络不通，沿用 DESIGN §3.3 逃生口：

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=<本地 googletest 解包目录>
```

（目录内容须为 f8d7d77c 对应源码；CI 侧不受影响。）

### 0.3 基线记录（实施前）

```powershell
ctest --test-dir build --output-on-failure        # 记录通过数（基线 N_core）
ctest --test-dir build_capi --output-on-failure   # 基线 N_capi
ctest --test-dir build_cuda --output-on-failure   # 基线 N_cuda（含 gpu_mc）
```

三项全绿为前提；任何先红 = 停止实施，先修复。

### 0.4 与 DESIGN 的章节对照（Step 8 对齐补入）

| 本文档 | DESIGN | 关系 |
|---|---|---|
| §0 环境与基线 | §4.1（构建/版本）、§6.7（复跑） | 执行前提 |
| §1 前置门禁 | §10（U-6/U-7/U-8/U-10）、§7 KAT 夹具、§8 风险退路 | 阻塞项落地 |
| §2 | §3（RISK-033） | 1:1 |
| §3 | §4（RISK-032） | 1:1 |
| §4 | §5（RISK-034） | 1:1 |
| §5 | §6（RISK-027）+ §4.3 文档行 | 1:1 + 合并提交约束 |
| §6 测试矩阵 | §7（TDD 红→绿顺序） | 展开 |
| §7 提交序列 | §8（回滚单元） | 展开 |
| §8 验收映射 | §9（A1-A10） | 1:1 |
| §9 实施期偏差 | §11（D-1~D-4） | E-1~E-5 承接 |
| §10 Step 6 review | 无（流程产物） | — |

---

## §1 实施第一步：阻塞确认与探针

按顺序执行，全部通过才进入 §2。任何一步失败 → 按"退路"列处置并回写 DESIGN。

| # | 项 | 动作 | 退路 |
|---|---|---|---|
| 1 | **U-7**（发布号） | ~~向 Scott 确认~~ **已裁定（2026-10-22）：1.9.0**（G-01 ✅） | — |
| 2 | **U-8**（vendor 源） | ~~下载~~ **已执行（2026-10-22）**：GitHub API 取 main SHA = `9545ff64…1150`；codeload 直连下载归档 430 644 B；拷贝 `include/Random123/**`（28 文件）+ LICENSE + UPSTREAM.md（G-02 ✅） | 网络阻断则经 ghproxy 镜像 / 手动下载 |
| 3 | **U-6**（KAT 列序） | ~~下载上游 `examples/kat_vectors`~~ **已执行**：实际路径为 **`tests/kat_vectors`**（非 examples，勘误 E-7）；列序 = `nameNxW R CTR[4] KEY[2] EXPECTED[4]`；philox4x64 **10 轮仅 3 行**（勘误 E-8）；夹具已生成（G-03 ✅） | 人工抽取直接写进头文件 |
| 4 | **U-10**（调用形态探针） | 最小 .cu 探针（不入库）验证三件事 | 见 1.2 |
| 5 | KAT 夹具生成 | ~~一次性脚本~~ **已执行（2026-10-22）**：[tmp/extract_kat.py](../../tmp/extract_kat.py) → `tests/fixtures/rng/kat_philox4x64_10.hpp`，`kKatCount = 3`（G-03 ✅） | 见 §1.3 |

### 1.1 vendor 布局（U-8 落地物）

```
third_party/random123/
├── include/Random123/        # 上游头文件原样拷贝（philox.h, array.h, features/* 等十余个）
├── LICENSE                   # 上游 LICENSE 原文（BSD-3-Clause）
└── UPSTREAM.md               # 来源 URL + 固定 commit + 拷贝日期
```

`UPSTREAM.md` 模板：

```markdown
# Random123（vendored）

- 上游: https://github.com/DEShawResearch/random123（RESEARCH §7.1 E1 锚点；raw philox.h 即取自该仓）
- 固定 commit: <U-8 实测 SHA>（main 分支，<拷贝日期>）
- 拷贝范围: include/Random123/** 原样，无任何修改
- 许可证: BSD-3-Clause（见 LICENSE），与本项目 Apache-2.0 兼容
- 依据: RISK-027 / DESIGN §6.1（不用 submodule：避免下游 wheel 构建追加 recursive 负担）
```

### 1.2 U-10 探针（tmp/probe_r123.cu，不入库）—— **已执行（2026-10-22），结论见下**

**执行结果**：host（MSVC）侧正常（`PHILOX_64BIT=1`，输出 = 官方 KAT）；**device（nvcc）侧首轮失败**：`philox.h(356): calling a __host__ function("mulhilo64") from a __device__ function`。根因：Windows 下 `R123_USE_MULHILO64_MSVC_INTRIN=1` 与 `R123_USE_MULHILO64_CUDA_INTRIN=1` 同时成立，而 philox.h L199-202 的分派顺序让 MSVC 内建（host-only `_umul128`）优先。
**修复（已验证）**：device pass 关闭 MSVC 分派 ——

```cpp
#if defined(__CUDA_ARCH__)
#define R123_USE_MULHILO64_MSVC_INTRIN 0
#endif
#include <Random123/philox.h>
```

修复后 `nvcc -arch=sm_89` 编译通过，device 输出 `16554d9eca36314c db20fe9d672d0fdc d7e772cee186176b 7e68b68aec7ba23b` = 官方 KAT 第 1 行，逐字一致。此守卫**必须**进入 `gpu_mc.cu`（见 §5.4 第 2' 条）与探针段。记录为 DESIGN §11 D-6。

原始探针代码（形态验证用）：

```cuda
// U-10 探针: 验证 Random123 调用形态 + 64 位开关 (RISK-027 前置)
#include <Random123/philox.h>
#include <cstdio>
#include <cstdint>

__global__ void probe_kernel(uint64_t* out) {
    r123::Philox4x64::ctr_type c{};  c.v[0] = 0;   // 期望形态: 实例(ctr, key) → ctr_type
    r123::Philox4x64::key_type k{};  k.v[0] = 0;
    r123::Philox4x64::ctr_type r = r123::Philox4x64{}(c, k);
    for (int i = 0; i < 4; ++i) out[i] = r.v[i];
}

int main() {
#if defined(R123_USE_PHILOX_64BIT)
    std::printf("R123_USE_PHILOX_64BIT = %d\n", R123_USE_PHILOX_64BIT);
#else
    std::printf("R123_USE_PHILOX_64BIT undefined!\n");
#endif
    // host 端同式计算 4 字, cudaMemcpy 比对 → 全等则形态正确
    ...  // 省略: cudaMalloc/copy/比较, 与 host 逐字相等才继续
}
```

```powershell
nvcc -I third_party/random123/include tmp\probe_r123.cu -o tmp\probe_r123.exe
tmp\probe_r123.exe
```

**判定**：(a) nvcc 编译通过；(b) `R123_USE_PHILOX_64BIT` 非 0（上游 [philox.h:175](https://github.com/AMitherd/random123/blob/main/include/Random123/philox.h#L175) `#if` 分支生效，Step 4 已页面级核实该行存在）；(c) 设备端 4 字与 host 端逐字相同。
**退路**：若特性宏在 nvcc/MSVC 下探测失败 → 只在 `UPSTREAM.md` 记录最小补丁（如 `-DR123_USE_PHILOX_64BIT=1`），并回写 DESIGN §10 U-10；若完全无法编译 → 停止 R-027，回 DESIGN 重议（其余三项不受影响）。

### 1.3 KAT 夹具形态（D-4：.hpp 生成头文件，弃 .txt+运行时解析器）

**决策**：夹具用生成头文件 `tests/fixtures/rng/kat_philox4x64_10.hpp`，沿用 grunfeld 先例（[tests/CMakeLists.txt:336-337](../../tests/CMakeLists.txt#L336-L337) `target_include_directories(... fixtures/econometrics)` + 生成头 [grunfeld_data.hpp](../../tests/fixtures/econometrics/grunfeld_data.hpp)）。理由：

1. ctest 工作目录差异 → 运行时解析 .txt 需要 CWD 敏感的路径逻辑（既有仓库已用生成头规避过一次）；
2. 列序解析（U-6）收进一次性抽取脚本，产物是零解析成本的编译期常量；
3. GPU 侧（test_gpu_mc）同样直接 include 该头。

**抽取脚本**（`tmp/extract_kat.py`，一次性，不入库）：读上游 `kat_vectors` → 过滤 `philox4x64` 且 10 轮的行 → 按实际列序（U-6）产出：

```cpp
// 生成头: 由 tmp/extract_kat.py 从上游 kat_vectors 抽取 (commit <U-8>), 请勿手改。
#pragma once
#include <cstdint>
#include <cstddef>

namespace cpphub::v1::rng_kat {

struct KatVector {
    std::uint64_t ctr[4];
    std::uint64_t key[2];
    std::uint64_t answer[4];  // philox4x64_10 (Random123 默认轮数)
};

// 上游 philox4x64 10 轮**全部** 3 行 (0/0, fff../fff.., pi-digits)
inline constexpr KatVector kKatVectors[] = {
    {{0,0,0,0}, {0,0}, {0x16554d9eca36314cULL, 0xdb20fe9d672d0fdcULL,
                      0xd7e772cee186176bULL, 0x7e68b68aec7ba23bULL}},
    // ... 共 3 行（实际生成头见 tests/fixtures/rng/kat_philox4x64_10.hpp）
};
inline constexpr std::size_t kKatCount = sizeof(kKatVectors) / sizeof(kKatVectors[0]);

}  // namespace cpphub::v1::rng_kat
```

**实际抽取结果（Step 9 已落盘）**：3 行向量（上游 philox4x64 10 轮的**全部**行，非抽样）；第 1 行答案与我方 nvcc 探针 device 输出逐字一致（独立交叉验证）。

---

## §2 RISK-033 — googletest 钉 commit（提交 1）

TDD 不适用（构建面变更，无 gtest 红用例）；验证 = 重放。

### 2.1 [CMakeLists.txt:61-68](../../CMakeLists.txt#L61-L68)

原（tag 归档 URL）→ 改为（DESIGN §3.1 / §11 D-5 裁决：URL@commit 形态）：

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

### 2.2 cache key（2 处定义）

- [ci.yml:59](../../.github/workflows/ci.yml#L59)（build-test）：`key: gtest-1.14.0-${{ runner.os }}-${{ matrix.compiler }}` → `gtest-f8d7d77c-…`
- [ci.yml:110](../../.github/workflows/ci.yml#L110)（c-api）：`key: gtest-1.14.0-${{ runner.os }}-capi` → `gtest-f8d7d77c-…`

python job 无 `_deps` 缓存（仅 :151 pip cache），不动。

### 2.3 验证（E1 可重放）

```powershell
# 清缓存全新 configure（URL@commit 形态：本机 codeload 可达，无需覆盖）
cmake -S . -B build_verify -G "Visual Studio 17 2022"
# 归档解包目录内容确认（URL 模式无 git 元数据，改用内容指纹）
Get-Content build_verify\_deps\googletest-src\CMakeLists.txt | Select-String 'GOOGLETEST_VERSION'  # 期望 1.14.0
ctest --test-dir build_verify -C Release --output-on-failure               # 全绿后删除 build_verify
```

（Step 9 勘误 E-6：原 C-13 用 `git rev-parse` 验证 commit，但 URL 解包无 `.git`，改以版本内容 + 归档 URL 的 SHA 路径为身份证据。）

---

## §3 RISK-032 — 版本两轴 + README 事实 + CI（提交 2）

### 3.1 TDD 序（本项是唯一有真红用例的流程面）

1. **红**：先改 [python/pyproject.toml:24](../../python/pyproject.toml#L24) `version = "1.4.0"` → `"1.9.0"`（U-7 确认后；若 Scott 改号，替换之）。
2. **红**：新增 [python/tests/test_version.py](../../python/tests/test_version.py)：

```python
"""RISK-032: __version__ 单一源一致性（A2 验收）。

直接正则抽取 pyproject 的 version, 避免 tomllib 的 3.11+ 门槛
（pyproject requires-python >= 3.9, 本测试随包可在 3.9 上运行）。
"""
import pathlib
import re

import cpphub


def test_version_matches_pyproject() -> None:
    pyproject = pathlib.Path(__file__).resolve().parents[1] / "pyproject.toml"
    text = pyproject.read_text(encoding="utf-8")
    m = re.search(r'^version\s*=\s*"([^"]+)"', text, re.M)
    assert m is not None, "pyproject.toml 缺少 project.version"
    assert cpphub.__version__ == m.group(1)
```

   （此刻 `__version__` 硬编码 `"1.4.0"` ≠ `1.9.0` → 红。本地跑法：先 `pip install "./python[test]"` 旧码装好，改 pyproject 后重跑 pytest。）
3. **绿**：完成 3.2 步骤 E/F（构建注入）→ 重装 wheel → pytest 绿。

### 3.2 文件级步骤（顺序执行）

| # | 文件 | 改动 |
|---|---|---|
| A | `VERSION`（新增） | 内容一行 `1.0.0`（+换行） |
| B | `tools/sync_version.py`（新增） | 完整脚本见 3.3 |
| C | `include/cpphub/core/version.hpp` | 运行 `python tools/sync_version.py` 生成 → **随源码树提交**（DESIGN §4.2 B / F8） |
| D | [config.hpp:4-8](../../include/cpphub/core/config.hpp#L4-L8) | 删 4 行手写宏，改为下方片段 |
| E | [python/CMakeLists.txt:10-11](../../python/CMakeLists.txt#L10-L11) + :69 后 | 两处，见 3.4 |
| F | [cpphub_python.cpp:2043](../../python/src/cpphub_python.cpp#L2043) | `m.attr("__version__") = "1.4.0";` → 宏读取，见 3.5 |
| G | [ci.yml:125](../../.github/workflows/ci.yml#L125) | 去 `-R "CApi"` 过滤（DESIGN §4.5），见 3.6 |
| H | [ci.yml:127](../../.github/workflows/ci.yml#L127) 前 | 插入 coverage job（DESIGN §4.4 YAML，已按 R-5/R-6 回写 DESIGN） |
| I | README.md :21 / :107 / :121 / 顶部 | DESIGN §4.3 表 4 行（**:40 与 USAGE_GUIDE 段归 R-027 提交，见 §5.6，避免二次改同段**） |

config.hpp（保留宏名 ⇒ 包含方零改动）：

```cpp
// 版本（单一源: 仓库根 VERSION → tools/sync_version.py 生成, RISK-032）
#include "cpphub/core/version.hpp"
```

### 3.3 tools/sync_version.py（完整）

```python
#!/usr/bin/env python3
"""RISK-032: 版本单一源同步。

- 无参数: 读仓库根 VERSION → 写 include/cpphub/core/version.hpp
- --check: 只校验, 不一致退出码 1（CI 用, 防生成物与真值漂移）

单一真值 = VERSION。version.hpp 为生成物且随源码树提交
（源码树 header-only 用法 clone 后即可编译, 无生成步骤可先行）。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION_FILE = ROOT / "VERSION"
HEADER = ROOT / "include" / "cpphub" / "core" / "version.hpp"

TEMPLATE = """\
// 本文件由 tools/sync_version.py 从仓库根 VERSION 生成, 请勿手改。
// 单一真值: VERSION（RISK-032 两轴模型: 本文件 = C++ 库 API 版本;
// Python 发行版本见 python/pyproject.toml）。
#pragma once
#define CPPHUB_VERSION_MAJOR {major}
#define CPPHUB_VERSION_MINOR {minor}
#define CPPHUB_VERSION_PATCH {patch}
#define CPPHUB_VERSION "{version}"
"""


def read_version() -> str:
    text = VERSION_FILE.read_text(encoding="utf-8").strip()
    if not re.fullmatch(r"\d+\.\d+\.\d+", text):
        sys.exit(f"VERSION 文件格式非法: {text!r} (期望 MAJOR.MINOR.PATCH)")
    return text


def expected_content(version: str) -> str:
    major, minor, patch = version.split(".")
    return TEMPLATE.format(major=major, minor=minor, patch=patch, version=version)


def main() -> int:
    version = read_version()
    content = expected_content(version)
    if "--check" in sys.argv[1:]:
        if HEADER.exists() and HEADER.read_text(encoding="utf-8") == content:
            print(f"version.hpp 与 VERSION 一致: {version}")
            return 0
        sys.exit(f"version.hpp 与 VERSION ({version}) 不一致 — 运行: python tools/sync_version.py")
    # 用 Path.open 而非 write_text(newline=...): 后者需 Python 3.10,
    # 而本项目 requires-python >= 3.9 (Step 6 勘误 R-3)。
    with HEADER.open("w", encoding="utf-8", newline="\n") as f:
        f.write(content)
    print(f"已写入 {HEADER.relative_to(ROOT)} ({version})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

（`sys.exit(str)` 即退出码 1 + stderr 输出；`newline="\n"` 保证跨平台行尾稳定，否则 `--check` 会在 Windows 上因 CRLF 误报。）

### 3.4 python/CMakeLists.txt（两处）

第一处，替换 [:10-11](../../python/CMakeLists.txt#L10-L11)（`cmake_minimum_required` 保留）：

```cmake
# RISK-032: Python 发行版本单一源 = pyproject.toml。
# wheel 构建时 scikit-build-core 注入 SKBUILD_PROJECT_VERSION;
# 直接构建 python/ 目录时回退到仓库根 VERSION (库版本), 保证无硬编码。
if(DEFINED SKBUILD_PROJECT_VERSION)
    set(CPPHUB_PY_VERSION "${SKBUILD_PROJECT_VERSION}")
else()
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../VERSION" CPPHUB_PY_VERSION LIMIT_COUNT 1)
endif()
project(cpphub_python VERSION ${CPPHUB_PY_VERSION} LANGUAGES CXX)
```

第二处，[:69](../../python/CMakeLists.txt#L69) `target_link_libraries(_core PRIVATE cpphub_headers)` 之后追加：

```cmake
# RISK-032: __version__ 由构建系统注入 (cpphub_python.cpp 消费 CPPHUB_PY_VERSION 宏)
target_compile_definitions(_core PRIVATE CPPHUB_PY_VERSION="${CPPHUB_PY_VERSION}")
```

### 3.5 cpphub_python.cpp:2043

```cpp
#ifndef CPPHUB_PY_VERSION
#define CPPHUB_PY_VERSION "unknown"  // 构建系统未注入时的兜底 (RISK-032)
#endif
m.attr("__version__") = CPPHUB_PY_VERSION;  // 由构建系统注入, 不再硬编码
```

（`#ifndef` 紧邻赋值处，兜底值为 `"unknown"` —— 不再写任何真实版本号，漂移可被 test_version.py 捕获。）

### 3.6 ci.yml（三处改动）

**3.6.1 过滤修复（[ci.yml:125](../../.github/workflows/ci.yml#L125)）**

```yaml
      - name: Run C API Tests
        working-directory: build
        run: ctest --build-config ${{ env.CMAKE_BUILD_TYPE }} --output-on-failure --parallel 4
```

（保持 `${{ env.CMAKE_BUILD_TYPE }}` 而非 DESIGN 片段中字面 `Release` —— 与 env 定义同源；语义等价。此差异记入 §9 E-3。）

**3.6.2 版本漂移门禁（新增步骤 —— Step 6 勘误 R-4）**

DESIGN §4.2 B 要求"CI `--check` 防漂移"，但 §4.4 YAML 未含该步骤。在 `build-test` job 的 `Configure CMake` 之后加：

```yaml
      - name: Check version header (RISK-032)
        if: runner.os == 'Linux'
        shell: bash
        run: python3 tools/sync_version.py --check
```

（`ubuntu-latest` 自带 `python3`；脚本仅用标准库。限定 Linux 避免 Windows python 别名差异。）

**3.6.3 coverage job（DESIGN §4.4 YAML，两处实施修正）**

按 DESIGN §4.4 原文插入到 c-api job 之后、python job 注释块（[:127-131](../../.github/workflows/ci.yml#L127-L131)）之前，但修正两点：

1. `python -m pip install gcovr` → **`python3 -m pip install gcovr`**（ubuntu-latest 无 `python` 别名；Step 6 勘误 R-5）；
2. `upload-artifact` 的 `path` 中 `coverage.html` → **`coverage*.html`**（`--html-details` 生成主文件 + 每模块明细，只传主文件会使下钻链接失效；Step 6 勘误 R-6）。

其余参数（Debug 避开 `-march` 守卫、`--coverage`、`--filter 'include/cpphub/'`、无阈值）不变。

### 3.7 根 CMakeLists.txt:2

```cmake
# RISK-032: 版本单一源 = 仓库根 VERSION
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" CPPHUB_LIB_VERSION LIMIT_COUNT 1)
project(CppHub VERSION ${CPPHUB_LIB_VERSION} LANGUAGES CXX)
```

### 3.8 验证

```powershell
python tools/sync_version.py --check          # 退出码 0
cmake -S . -B build                            # 重配, 确认 project 版本 = 1.0.0
ctest --test-dir build --output-on-failure     # 全绿
pip install "./python[test]" ; python -m pytest python/tests/test_version.py -v   # 绿, __version__=1.9.0
# A3: README 冻结数字清零
# grep "1412" README.md / grep "15 GPU MC" README.md  → 零命中
```

---

## §4 RISK-034 — json_min 输入预算（提交 3）

### 4.1 常量（2 处新增）

[json_min.hpp:27](../../include/cpphub/c_api/detail/json_min.hpp#L27) `kMaxParseDepth` 之后：

```cpp
// 节点预算上限（RISK-034）：每次 parse_value 调用计 1 个节点
// （null/bool/number/string/array/object 各计 1），对象键不计入——
// 键的规模已被调用方入口长度上限约束，且键不是 JSON 值。
// 合法上界 ≈ kMaxSeriesLength(1e4) × kMaxColumns(64) ≈ 6.4e5 个数值节点，
// 1e6 留 ~1.5× 余量；未来放宽序列上限时两者需同步。
inline constexpr std::size_t kMaxParseNodes = 1'000'000;
```

[cpphub_tool_registry.hpp:52](../../include/cpphub/c_api/cpphub_tool_registry.hpp#L52) `kMaxColumns` 之后：

```cpp
// args_json 字节上限（RISK-034）：属工具面输入策略（与序列/列数上限同层），
// 由 c_api.cpp tool_call 入口在 strlen 后前置检查，防超大输入在解析前占用内存。
inline constexpr std::size_t kMaxArgsJsonBytes = 8u * 1024u * 1024u;  // 8 MiB
```

### 4.2 解析器节点计数（Parser）

- [Parser](../../include/cpphub/c_api/detail/json_min.hpp#L157) 私有成员区（:174 `err_` 旁）加 `std::size_t nodes_ = 0;`
- [parse_value](../../include/cpphub/c_api/detail/json_min.hpp#L192) 入口（switch 之前）加：

```cpp
        if (++nodes_ > kMaxParseNodes) {
            return fail("maximum node count exceeded");
        }
```

**计数语义**（写入常量注释并被测试锚定）：每调用一次 `parse_value` = 1 节点；对象键（`parse_string` 直调）不经过 `parse_value` ⇒ 天然不计入。

### 4.3 入口长度检查（[c_api.cpp](../../src/c_api/c_api.cpp#L292-L331)）

在 `args_json` 判空（:303-306）之后、`find_tool`（:308）之前插入：

```cpp
        // RISK-034: 入口长度上限, 前置于解析与工具查找
        const std::size_t args_len = std::strlen(args_json);
        if (args_len > cpphub::v1::c_api::kMaxArgsJsonBytes) {
            set_error("args_json exceeds limit (" + std::to_string(args_len) +
                      " > " + std::to_string(cpphub::v1::c_api::kMaxArgsJsonBytes) +
                      " bytes)");
            return CPPHUB_ERR_INVALID_ARG;
        }
```

（`<cstring>` 已在 [c_api.cpp:18](../../src/c_api/c_api.cpp#L18)；错误沿用 `CPPHUB_ERR_INVALID_ARG`，**不新增错误码**——C ABI v1 纯增约定。）

### 4.4 测试（5 用例，DESIGN §5.5 全集）

test_json_min.cpp 追加（既有 `RejectsDepthOverflow` :321 同风格）：

```cpp
// ---------------------------------------------------------------------------
// RISK-034: 节点预算（kMaxParseNodes）
// 注意: 文件顶部 (:18) 已有 `namespace jm = cpphub::v1::json_min;`,
// 同一作用域重复声明 namespace alias = 编译错误 (Step 6 勘误 R-1), 此处不得再写。
// ---------------------------------------------------------------------------

TEST(JsonMinBudget, AcceptsAtBoundary) {
    // N 元素数组 = 1(数组) + N(元素) 节点 → kMaxParseNodes-1 个元素恰为 1e6 节点
    std::string text = "[";
    for (std::size_t i = 0; i < jm::kMaxParseNodes - 1; ++i) {
        if (i) text += ',';
        text += '0';
    }
    text += ']';
    std::string err;
    jm::Value v;
    ASSERT_TRUE(jm::parse(text, v, err)) << err;
    ASSERT_TRUE(v.is_array());
    EXPECT_EQ(v.as_array().size(), jm::kMaxParseNodes - 1);
}

TEST(JsonMinBudget, RejectsWideArrayOverBudget) {
    // kMaxParseNodes+2 个元素 → 1'000'003 节点, 超限拒绝且错误含 "node"
    std::string text = "[";
    for (std::size_t i = 0; i < jm::kMaxParseNodes + 2; ++i) {
        if (i) text += ',';
        text += '0';
    }
    text += ']';
    std::string err;
    jm::Value v;
    EXPECT_FALSE(jm::parse(text, v, err));
    EXPECT_NE(std::string(err).find("node"), std::string::npos) << err;
}

TEST(JsonMinBudget, DeepNestingStillLimitedByDepth) {
    // 回归锁: 节点预算不遮蔽深度限制（小节点数 + 深嵌套仍报 nesting depth）
    std::string deep;
    for (int i = 0; i < 200; ++i) deep.push_back('[');
    deep += "1";
    for (int i = 0; i < 200; ++i) deep.push_back(']');
    std::string err;
    jm::Value v;
    EXPECT_FALSE(jm::parse(deep, v, err));
    EXPECT_NE(std::string(err).find("nesting depth"), std::string::npos) << err;
}
```

test_c_api_tools.cpp 追加（文件顶部补 `#include "cpphub/c_api/cpphub_tool_registry.hpp"`，为取 `kMaxArgsJsonBytes`）：

```cpp
// ---------------------------------------------------------------------------
// RISK-034: 入口长度上限（kMaxArgsJsonBytes, 8 MiB）
// ---------------------------------------------------------------------------
TEST(CApiTools, RejectsOversizedArgsJson) {
    const std::size_t limit = cpphub::v1::c_api::kMaxArgsJsonBytes;
    std::string args = "{\"series\":[" + std::string(limit, '0') + "]}";
    ASSERT_GT(args.size(), limit);
    size_t need = 0;
    EXPECT_EQ(cpphub_v1_tool_call("cpphub_unit_root_diagnose", args.c_str(),
                                  nullptr, &need),
              CPPHUB_ERR_INVALID_ARG);
    EXPECT_NE(std::string(cpphub_v1_get_last_error()).find("exceeds limit"),
              std::string::npos);
}

TEST(CApiTools, OversizedCheckPrecedesParse) {
    // 超长 + 语法非法 → 报 "exceeds limit" 而非 JSON 语法错误（证明检查前置）
    std::string args(cpphub::v1::c_api::kMaxArgsJsonBytes + 1, '!');
    size_t need = 0;
    EXPECT_EQ(cpphub_v1_tool_call("cpphub_unit_root_diagnose", args.c_str(),
                                  nullptr, &need),
              CPPHUB_ERR_INVALID_ARG);
    const std::string err = cpphub_v1_get_last_error();
    EXPECT_NE(err.find("exceeds limit"), std::string::npos) << err;
    EXPECT_EQ(err.find("invalid args_json"), std::string::npos) << err;
}
```

**TDD 红→绿性质**：`RejectsWideArrayOverBudget` / `RejectsOversizedArgsJson` / `OversizedCheckPrecedesParse` = **红**（现无上限：超宽解析成功 / 超长报 `invalid args_json` 语法错）；`AcceptsAtBoundary` / `DeepNestingStillLimitedByDepth` = **回归锁**（实现前后均应绿，防回归）。
**资源**：1e6 元素用例峰值 ≈ 1e6 × sizeof(Value)（Value 含 2 个 `std::vector` + `std::string`，约 96–112 B ⇒ ~100 MB，叠加 vector 扩容峰值 ~150–200 MB）+ ~2 MB 文本；CI runner（7 GB+）可承受。8 MiB 用例为一次 strlen + 一次字符串构造，可忽略。

### 4.5 验证（build_capi）

```powershell
cmake -S . -B build_capi -DCPPHUB_ENABLE_C_API=ON
cmake --build build_capi --config Release --parallel
ctest --test-dir build_capi --build-config Release -R "JsonMinBudget|CApiTools" --output-on-failure
ctest --test-dir build_capi --build-config Release --output-on-failure   # 全量回归
```

（性能影响：入口 strlen 一次 O(n) + 每值一次计数比较，相对解析本身可忽略，不单列 benchmark。）

---

## §5 RISK-027 — vendor Random123 + 块对齐 + KAT（提交 4，最大改动）

前置：§1 全部完成（U-8 vendor 就位、U-10 探针通过、KAT 夹具已生成）。

### 5.1 include 路径（三处，DESIGN §6.2）

| 目标 | 位置 | 新增行 |
|---|---|---|
| `cpphub_core` | [CMakeLists.txt:25-28](../../CMakeLists.txt#L25-L28) `target_include_directories` 列表 | `${CMAKE_CURRENT_SOURCE_DIR}/third_party/random123/include` |
| `cpphub_cuda`（CUDA 分支） | [src/performance/gpu/CMakeLists.txt:36-39](../../src/performance/gpu/CMakeLists.txt#L36-L39) | `${GPU_THIRD_PARTY_DIR}/random123/include` |
| `cpphub_headers` | [python/CMakeLists.txt:27-31](../../python/CMakeLists.txt#L27-L31) | `${CMAKE_CURRENT_SOURCE_DIR}/../third_party/random123/include` |

stub 分支（:72-75）不编 gpu_mc.cu，**不需要** random123 路径（DESIGN §12 已确认）。

### 5.2 测试目标 CMake（R-027 相关）

- [tests/CMakeLists.txt:22](../../tests/CMakeLists.txt#L22) `cpphub_add_test(test_rng ...)` 之后：

```cmake
# RISK-027: KAT 夹具为生成头 (grunfeld 模式), test_rng 直接 include
target_include_directories(test_rng PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/rng)
```

- [src/performance/gpu/CMakeLists.txt](../../src/performance/gpu/CMakeLists.txt#L93-L114) `test_gpu_mc` 块内（:102 后）：

```cmake
    # RISK-027: KAT 夹具 + CPU/GPU 对齐用例
    target_include_directories(test_gpu_mc PRIVATE ${GPU_TESTS_DIR}/fixtures/rng)
```

- [src/performance/gpu/CMakeLists.txt:27-34](../../src/performance/gpu/CMakeLists.txt#L27-L34) `FILE_SET cuda_headers` FILES 追加一行（保持头文件登记惯例）：

```cmake
            ${GPU_INCLUDE_DIR}/cpphub/performance/gpu/rng_device_probe.hpp
```

### 5.3 rng.hpp 全量重写（公共契约不变，DESIGN §6.3）

**边界**：删 [L10-13](../../include/cpphub/core/rng.hpp#L10-L13)（MSVC intrin，上游接管）、[L18-21](../../include/cpphub/core/rng.hpp#L18-L21)（错误常量：`0xD2B74407B1CE6E93` 实为 Random123 的 **M2x64_0**，`0x9CCD9A8C5AE13B41` 非任何 Philox4x64 常量，`PHILOX_W32_*` 是 4x32 的 W——名实不符根因）、[L23-102](../../include/cpphub/core/rng.hpp#L23-L102)（手写类）。**保留** [L104-136](../../include/cpphub/core/rng.hpp#L104-L136)（box_muller / normal_simd / generate_correlated 逐字不动）与新 include 区。

新实现（全量替换后文件头部）：

```cpp
#pragma once
#include <cstdint>
#include <array>
#include <utility>
#include <limits>
#include <Random123/philox.h>   // RISK-027: vendor Random123 (third_party/random123)
#include "cpphub/core/constants.hpp"
#include "cpphub/core/simd.hpp"
#include "cpphub/core/linalg.hpp"

namespace cpphub {
inline namespace v1 {

// Philox4x64-10 (Random123): counter 4xu64, key 2xu64, 每次调用产出 4 字。
// 本类为流式适配层: 逐字输出 + 块缓冲, 对外语义见 DESIGN §6.3
// （operator() = 流顺序下一字; next4() = 恰 1 个完整块 = 与 GPU 每路径 1 块对齐）。
class Philox4x64 {
    r123::Philox4x64::ctr_type ctr_{};    // 4xu64 计数器（下一块的 counter）
    r123::Philox4x64::key_type key_{};    // 2xu64 密钥 = {seed, stream}
    r123::Philox4x64::ctr_type block_{};  // 当前缓冲块的 4 字
    int idx_ = 4;                          // 4 = 块已耗尽

public:
    using result_type = uint64_t;
    static constexpr uint64_t min() { return 0; }
    static constexpr uint64_t max() { return UINT64_MAX; }

    explicit Philox4x64(uint64_t seed = 0, uint64_t stream = 0) {
        key_.v[0] = seed;
        key_.v[1] = stream;
    }

    // 流顺序的下一字; 块耗尽时以当前 counter 计算整块并前移
    uint64_t operator()() {
        if (idx_ == 4) refill();
        return block_.v[idx_++];
    }

    // 恰 1 个完整块（4 字）。调用方若从块边界出发则保持块对齐
    // （normal_simd 即此用法 → CPU 第 p 步 = GPU 第 p 路径同一 counter 同一 4 字）。
    std::array<uint64_t, 4> next4() {
        return {operator()(), operator()(), operator()(), operator()()};
    }

    // 丢弃 n 个输出字（DESIGN §6.3: ctr += n/4, 再消费 n%4）
    void discard(uint64_t n) {
        if (n == 0) return;
        const uint64_t avail = (idx_ == 4) ? 0 : static_cast<uint64_t>(4 - idx_);
        if (n < avail) {
            idx_ += static_cast<int>(n);
            return;
        }
        n -= avail;
        idx_ = 4;
        const uint64_t blocks = n / 4;
        if (blocks > 0) add_to_ctr(blocks);
        const uint64_t rem = n % 4;
        if (rem > 0) {
            refill();
            idx_ = static_cast<int>(rem);
        }
        // rem == 0: 保持 idx_=4, 下次 operator() 触发 refill
    }

private:
    void refill() {
        block_ = r123::Philox4x64{}(ctr_, key_);
        bump_ctr();
        idx_ = 0;
    }
    void bump_ctr() {
        for (int i = 0; i < 4; ++i) {
            if (++ctr_.v[i] != 0) break;  // 256 位小端进位
        }
    }
    void add_to_ctr(uint64_t k) {
        for (int i = 0; i < 4 && k != 0; ++i) {
            const uint64_t old = ctr_.v[i];
            ctr_.v[i] = old + k;
            k = (ctr_.v[i] < old) ? 1 : 0;
        }
    }
};
```

（后接原 L104-136 三个函数，原样保留。轮数/常量由上游 `Philox4x64 = Philox4x64_R<PHILOX4x64_DEFAULT_ROUNDS=10>` 接管。）

**discard 语义自检**（对应既有 [test_rng.cpp:93-101](../../tests/unit/core/test_rng.cpp#L93-L101) `PhiloxDiscard`）：`discard(50)` → 12 完整块 + 2 字内消费，随后 `operator()` 与逐字推进的 rng1 输出逐字相同 ✓。

### 5.4 gpu_mc.cu diff

1. 头注释（:4-8）"同样常量 (PHILOX_M4x64_0/1, PHILOX_W32_0/1)" 句改为："同一 Random123 Philox4x64-10 实现（third_party/random123），常量/轮数由上游接管"。
2. include 区（:19-31）加 `#include <Random123/philox.h>`。
2'. **include 前置宏守卫（Step 9 探针结论，DESIGN §11 D-6，必需）**：

```cpp
// Windows+nvcc: 关闭 MSVC 内建分派, 让 device 代码走 __umul64hi (Risk-027 探针实测)
// 否则 philox.h L199-202 让 host-only 的 _umul128 被 device 代码调用 → 编译失败
#if defined(__CUDA_ARCH__)
#define R123_USE_MULHILO64_MSVC_INTRIN 0
#endif
#include <Random123/philox.h>
```

（探针段（§5.5）同样需要该守卫，因为它也在 .cu 中 include philox.h。）
3. 删 :42-45（4 个 `__device__ __constant__`）、:47-55（gpu_mulhi/mullo）、:57-77（philox4x64_10）；保留匿名命名空间内 [gpu_box_muller:81-87](../../src/performance/gpu/gpu_mc.cu#L81-L87) 与 [uint64_to_unit:90-92](../../src/performance/gpu/gpu_mc.cu#L90-L92)。
4. 内核调用点 [（:133-141）](../../src/performance/gpu/gpu_mc.cu#L133-L141) 替换为（DESIGN §6.4）：

```cpp
        // Philox: Random123 vendor (RISK-027), key=[seed, 0], counter=[path, 0, 0, 0]
        r123::Philox4x64::ctr_type ctr{};
        ctr.v[0] = static_cast<uint64_t>(path);
        r123::Philox4x64::key_type key{};
        key.v[0] = seed;
        const r123::Philox4x64::ctr_type r = r123::Philox4x64{}(ctr, key);

        double u1 = uint64_to_unit(r.v[0]);
        double u2 = uint64_to_unit(r.v[1]);
        double z1, z2;
        gpu_box_muller(u1, u2, z1, z2);
```

（u1/u2 选取 `r.v[0]/r.v[1]` 与 CPU `normal_simd` 第 1 块的 `u[0]/u[1]` 同位 → uniform 层位精确；仍只用 z1。）

### 5.5 GPU 测试探针（新文件 + gpu_mc.cu 末尾追加）

**新头** `include/cpphub/performance/gpu/rng_device_probe.hpp`：

```cpp
// GPU RNG 测试探针（RISK-027）: 把设备端 Philox 块/KAT/Box-Muller 结果导出到 host,
// 供 test_gpu_mc.cpp（.cpp, 不能含 __global__）做 CPU/GPU 对齐断言。
// 仅在 CPPHUB_CUDA_AVAILABLE 构建下有定义（gpu_mc.cu 末尾实现）;
// CPU stub 构建下未链接 —— 调用方（test_gpu_mc）本身也仅在 CUDA 分支构建。
#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>

namespace cpphub::v1::gpu_probe {

struct KatVector {
    std::uint64_t ctr[4];
    std::uint64_t key[2];
};

// 第 p 块 = Philox4x64(ctr={p,0,0,0}, key={seed,0}) 的 4 字;
// 返回 4*n_blocks 个字, 块优先顺序 (block0[0..3], block1[0..3], ...)。
std::vector<std::uint64_t> philox_blocks(std::uint64_t seed, std::size_t n_blocks);

// 逐 KAT 向量在设备端计算, 返回 4*vectors.size() 个字（每向量 4 字, 按输入序）。
std::vector<std::uint64_t> philox_kat(const std::vector<KatVector>& vectors);

// 与 MC 内核同一 uint64_to_unit + gpu_box_muller 路径的 z1 批量结果（第 p 路径）。
std::vector<double> box_muller_z1_batch(std::uint64_t seed, std::size_t n);

}  // namespace cpphub::v1::gpu_probe
```

**gpu_mc.cu 文件末尾**（三个 `}  // namespace` 收尾之前）。Step 6 勘误 R-2 修正了三处：结构体改为直接复用探针头的 `gpu_probe::KatVector`（原草稿 `KatVectorIn` 先用后定义）；**不得**再写 `namespace cpphub { inline namespace v1 {`（会错层嵌套成 `cpphub::v1::cpphub`），只写 `namespace gpu_probe {`；三个 kernel 均补越界守卫与 `cudaGetLastError()` 检查。

```cpp
// 头文件区 (:19-31) 追加:
#include "cpphub/performance/gpu/rng_device_probe.hpp"

namespace {

constexpr int kProbeBlock = 128;

__global__ void probe_philox_block_kernel(uint64_t seed, uint64_t* out, uint64_t n) {
    const uint64_t path = blockIdx.x * blockDim.x + threadIdx.x;
    if (path >= n) return;  // grid 向上取整 → 必须守卫, 否则越界写
    r123::Philox4x64::ctr_type ctr{};
    ctr.v[0] = path;
    r123::Philox4x64::key_type key{};
    key.v[0] = seed;
    const r123::Philox4x64::ctr_type r = r123::Philox4x64{}(ctr, key);
    for (int i = 0; i < 4; ++i) out[path * 4 + i] = r.v[i];
}

__global__ void probe_philox_kat_kernel(const gpu_probe::KatVector* vecs,
                                        uint64_t* out, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    r123::Philox4x64::ctr_type ctr{};
    for (int j = 0; j < 4; ++j) ctr.v[j] = vecs[i].ctr[j];
    r123::Philox4x64::key_type key{};
    key.v[0] = vecs[i].key[0];
    key.v[1] = vecs[i].key[1];
    const r123::Philox4x64::ctr_type r = r123::Philox4x64{}(ctr, key);
    for (int j = 0; j < 4; ++j) out[i * 4 + j] = r.v[j];
}

__global__ void probe_box_muller_z1_kernel(uint64_t seed, double* out, uint64_t n) {
    const uint64_t path = blockIdx.x * blockDim.x + threadIdx.x;
    if (path >= n) return;
    r123::Philox4x64::ctr_type ctr{};
    ctr.v[0] = path;
    r123::Philox4x64::key_type key{};
    key.v[0] = seed;
    const r123::Philox4x64::ctr_type r = r123::Philox4x64{}(ctr, key);
    double z1, z2;
    gpu_box_muller(uint64_to_unit(r.v[0]), uint64_to_unit(r.v[1]), z1, z2);
    out[path] = z1;
}

unsigned probe_grid(std::size_t n) {
    return static_cast<unsigned>((n + kProbeBlock - 1) / kProbeBlock);
}

}  // namespace

namespace gpu_probe {   // 已在 cpphub::v1 内, 直接开子命名空间

std::vector<std::uint64_t> philox_blocks(std::uint64_t seed, std::size_t n_blocks) {
    if (n_blocks == 0) return {};
    std::vector<std::uint64_t> host(n_blocks * 4);
    DeviceBuffer dev(host.size() * sizeof(std::uint64_t));
    probe_philox_block_kernel<<<probe_grid(n_blocks), kProbeBlock>>>(
        seed, reinterpret_cast<std::uint64_t*>(dev.get()),
        static_cast<std::uint64_t>(n_blocks));
    CPPHUB_CUDA_CHECK(cudaGetLastError());
    CPPHUB_CUDA_CHECK(cudaMemcpy(host.data(), dev.get(),
                                 host.size() * sizeof(std::uint64_t),
                                 cudaMemcpyDeviceToHost));
    return host;
}

std::vector<std::uint64_t> philox_kat(const std::vector<KatVector>& vectors) {
    if (vectors.empty()) return {};
    std::vector<std::uint64_t> out(vectors.size() * 4);
    DeviceBuffer d_in(vectors.size() * sizeof(KatVector));
    DeviceBuffer d_out(out.size() * sizeof(std::uint64_t));
    CPPHUB_CUDA_CHECK(cudaMemcpy(d_in.get(), vectors.data(),
                                 vectors.size() * sizeof(KatVector),
                                 cudaMemcpyHostToDevice));
    probe_philox_kat_kernel<<<probe_grid(vectors.size()), kProbeBlock>>>(
        reinterpret_cast<const KatVector*>(d_in.get()),
        reinterpret_cast<std::uint64_t*>(d_out.get()),
        static_cast<int>(vectors.size()));
    CPPHUB_CUDA_CHECK(cudaGetLastError());
    CPPHUB_CUDA_CHECK(cudaMemcpy(out.data(), d_out.get(),
                                 out.size() * sizeof(std::uint64_t),
                                 cudaMemcpyDeviceToHost));
    return out;
}

std::vector<double> box_muller_z1_batch(std::uint64_t seed, std::size_t n) {
    if (n == 0) return {};
    std::vector<double> out(n);
    DeviceBuffer dev(n * sizeof(double));
    probe_box_muller_z1_kernel<<<probe_grid(n), kProbeBlock>>>(
        seed, dev.get(), static_cast<std::uint64_t>(n));
    CPPHUB_CUDA_CHECK(cudaGetLastError());
    CPPHUB_CUDA_CHECK(cudaMemcpy(out.data(), dev.get(), n * sizeof(double),
                                 cudaMemcpyDeviceToHost));
    return out;
}

}  // namespace gpu_probe
```

（`DeviceBuffer` 以 `Real*` 承载任意字节——u64/KatVector/double 数组均为整块分配回收；`reinterpret_cast` 仅用于 `cudaMalloc` 的 `void*` 语义等价转换。kernel 启动后统一 `CPPHUB_CUDA_CHECK(cudaGetLastError())`，沿本文件 RISK-026 惯例。探针函数只进 CUDA 分支静态库，stub 分支无需实现。）

### 5.6 KAT 与对齐测试（4 用例）

**test_rng.cpp**（既有 **9** 用例全为语义/统计断言、无固值，预期直接绿——实施时验证并记录）：`#include "kat_philox4x64_10.hpp"` 加到[文件头部 include 区](../../tests/unit/core/test_rng.cpp#L1-L6)，用例追加到文件尾部：

```cpp
// RISK-027: 官方 KAT（上游 kat_vectors 生成头, tmp/extract_kat.py 抽取）
TEST(RngKat, Philox4x64MatchesKnownAnswers) {
    for (std::size_t i = 0; i < rng_kat::kKatCount; ++i) {
        const auto& kv = rng_kat::kKatVectors[i];
        r123::Philox4x64::ctr_type ctr{};
        for (int j = 0; j < 4; ++j) ctr.v[j] = kv.ctr[j];
        r123::Philox4x64::key_type key{};
        key.v[0] = kv.key[0];
        key.v[1] = kv.key[1];
        const auto r = r123::Philox4x64{}(ctr, key);
        for (int j = 0; j < 4; ++j) {
            ASSERT_EQ(r.v[j], kv.answer[j]) << "vector " << i << " word " << j;
        }
    }
}
```

（需 `#include <Random123/philox.h>`——rng.hpp 已带出。）

**test_gpu_mc.cpp 追加**（suite 沿用既有 `GpuMC`）：

```cpp
#include "cpphub/core/rng.hpp"                       // CPU 侧 Philox4x64 / box_muller (R-8)
#include "cpphub/performance/gpu/rng_device_probe.hpp"
#include "kat_philox4x64_10.hpp"

// P2a: 设备端 Random123 与官方 KAT 一致（证明 nvcc 编译出同一实现）
TEST(GpuMC, PhiloxMatchesKnownAnswers) {
    std::vector<gpu_probe::KatVector> vecs(rng_kat::kKatCount);
    for (std::size_t i = 0; i < rng_kat::kKatCount; ++i) {
        for (int j = 0; j < 4; ++j) vecs[i].ctr[j] = rng_kat::kKatVectors[i].ctr[j];
        vecs[i].key[0] = rng_kat::kKatVectors[i].key[0];
        vecs[i].key[1] = rng_kat::kKatVectors[i].key[1];
    }
    const auto got = gpu_probe::philox_kat(vecs);
    for (std::size_t i = 0; i < rng_kat::kKatCount; ++i) {
        for (int j = 0; j < 4; ++j) {
            ASSERT_EQ(got[i * 4 + j], rng_kat::kKatVectors[i].answer[j])
                << "vector " << i << " word " << j;
        }
    }
}

// P2b: uniform 流 CPU/GPU 逐字相同（把 E1-9 固化为回归; U-1 分层承诺第 1 层）
TEST(GpuMC, UniformStreamMatchesCpu) {
    constexpr std::size_t kN = 10000;
    Philox4x64 cpu_rng(12345);
    const auto gpu = gpu_probe::philox_blocks(12345, kN);
    for (std::size_t p = 0; p < kN; ++p) {
        const auto blk = cpu_rng.next4();
        for (int j = 0; j < 4; ++j) {
            ASSERT_EQ(blk[j], gpu[p * 4 + j]) << "block " << p << " word " << j;
        }
    }
}

// P3: 正态层级 1e-15 相对差（U-1 分层承诺第 2 层; 回归锁, 非红驱动）
TEST(GpuMC, NormalSamplesWithinTolerance) {
    constexpr std::size_t kN = 10000;
    Philox4x64 cpu_rng(12345);
    const auto gpu = gpu_probe::box_muller_z1_batch(12345, kN);
    for (std::size_t p = 0; p < kN; ++p) {
        auto u = cpu_rng.next4();
        double u1 = (u[0] >> 11) * (1.0 / 9007199254740992.0);
        double u2 = (u[1] >> 11) * (1.0 / 9007199254740992.0);
        const double zc1 = box_muller(u1, u2).first;  // 只用 z1, 与 GPU 一致 (R-9)
        const double diff = std::fabs(zc1 - gpu[p]);
        const double denom = std::fabs(zc1);
        ASSERT_TRUE(denom > 1.0 ? diff / denom <= 1e-15 : diff <= 1e-15)
            << "path " << p << " cpu=" << zc1 << " gpu=" << gpu[p];
    }
}
```

（P3 容差形式 = DESIGN §6.5「相对差 ≤ 1e-15，|z|<1 时绝对兜底」；**Step 5 勘误**：该用例非红驱动——实测末位差 ≈1e-16 相对差本已低于 1e-15，且实施前探针不存在无法链接；性质为回归锁。已回写 DESIGN §7。）

### 5.7 文档修正（并入本提交，DESIGN §6.6 / §2 合并约束）

| 位置 | 改为 |
|---|---|
| [README.md:40](../../README.md#L40) | DESIGN §4.3 行 1 文案（uniform bit-exact + normal 1e-15） |
| [USAGE_GUIDE.md:183-196](../../docs/USAGE_GUIDE.md#L183-L196)（§3.1.3，**逐行改，R-7**） | :185 算法→Random123 `philox4x64_10`（Salmon et al. 2011）；**:190 周期 `2^128` → `2^256`（counter 4×u64）**；**:191 状态 `counter_[2]+key_[2]（256 位）` → `counter[4]+key[2]（384 位）`**；:192 输出→「每块 4×u64；`operator()` 逐字、`next4()` 一整块」；**:196 并行性**→分层承诺（uniform 位精确 / 正态 1e-15）+ 源码树使用者需 `-I third_party/random123/include` |
| [USAGE_GUIDE.md:693](../../docs/USAGE_GUIDE.md#L693) | GPU 行同步为同一表述（§3.10 表格仅此一行涉及 RNG；:696 的价格 1e-12 是求和顺序差异，与非 RNG 位级无关，不改） |
| [ARCHITECTURE.md:91](../../docs/architecture/ARCHITECTURE.md#L91) | 删 `Threefry4x64` 假条目，写 `Philox4x64`（Random123 vendor） |
| [ARCHITECTURE.md:68](../../docs/architecture/ARCHITECTURE.md#L68) | ADR-004 备注 vendor 路径 |

### 5.8 复跑处置与性能记录（DESIGN §6.7）

1. 夹具不受影响（tests/fixtures/** 由 Python/R 生成，不重算）。
2. 全量复跑三目录（§0.1），记录失败清单；处置仅两类：(a) 统计窗口红 → **重选 seed** 并在提交说明记录；(b) 实现错误 → 修实现。**禁止放宽容差**。
3. 性能前后：`benchmark_philox`（[tests/unit/performance/benchmark_philox.cpp](../../tests/unit/performance/benchmark_philox.cpp)，**路径在 tests/ 下而非 benchmarks/**；信息性、不入 ctest；`N=1e8` 逐字 `rng()` 调用）改造前后各跑一次，吞吐写入 [DEVELOPMENT_LOG.md](../../docs/DEVELOPMENT_LOG.md) 新增小节「R-027 前后吞吐」。旧值与新值**不可直接比回归**（算法实现更换），记录为"实现更换基线"。
   理论量级：旧实现每输出 1 字 = 1 次 10 轮 compute（产 2 字用 1 字）；新实现每 4 字 = 1 次块调用 ⇒ 每字轮数降到 1/4，**吞吐上限约 4×**；缓冲分支开销以实测为准，若下降 >2× 须在日志说明可接受性。
   注：benchmark 输出的平台标签（MSVC `_umul128` / GCC `__uint128_t`，[:15-25](../../tests/unit/performance/benchmark_philox.cpp#L15-L25)）描述的是旧实现的乘法路径，换 vendor 后仅剩"编译平台"含义；不改代码，记录时注明语境。
4. Python docstring：`next4()` 语义变化（4 字 = 1 块）在绑定层如有暴露需同步（实施时 grep `next4` 于 python/src）。

---

## §6 测试执行矩阵（TDD 红→绿总序）

| 阶段 | 用例 | 红性质 | 绿条件 | 命令（目录） |
|---|---|---|---|---|
| P1 | §4.4 三红 + 两回归锁 | 入口无上限 / 无节点预算 | §4.2+§4.3 落地 | `ctest -R "JsonMinBudget\|CApiTools"`（build_capi） |
| P2 | `RngKat.Philox4x64MatchesKnownAnswers`、`GpuMC.PhiloxMatchesKnownAnswers`、`GpuMC.UniformStreamMatchesCpu` | 实施前 KAT 不存在/探针无法链接（编译期红） | §5.3+§5.4+§5.5 落地 | test_rng（build）、test_gpu_mc（build_cuda） |
| P3 | `GpuMC.NormalSamplesWithinTolerance` | **非红（Step 5 勘误）**：≈1e-16 < 1e-15，回归锁 | 随 P2 | test_gpu_mc（build_cuda） |
| P4 | `tools/sync_version.py --check` + `test_version.py` | 无 version.hpp / `__version__` 硬编码 | §3.2 C/E/F 落地 | 脚本 + `pytest python/tests/test_version.py` |
| P5 | 全量回归（三目录） | 记录并按 §5.8 处置 | 全绿 | §0.3 命令 |

**既有测试影响预判**：test_rng **9** 用例（语义/统计断言，无固值）预期直接绿；`PhiloxDiscard`/`PhiloxNext4Consistency` 已按新语义自检（§5.3）。统计类（NormalSimdStats / GenerateCorrelated，容差 0.05/0.1）换流后按 §5.8 规则处置。

---

## §7 提交序列（4 commits，各含验证）

| # | 信息草案 | 内容 | 提交前动作 |
|---|---|---|---|
| 1 | `build(RISK-033): pin googletest to immutable commit f8d7d77c (v1.14.0)` | §2 全部 | §2.3 重放记录 |
| 2 | `feat(RISK-032): single-source version model + README fact fixes + CI coverage & C API test filter` | §3 全部（README :21/:107/:121 + 顶部版本模型；不含 :40） | sync --check + pytest + grep A3 |
| 3 | `feat(RISK-034): json_min input budgets (8 MiB entry cap + 1e6 node budget)` | §4 全部 | build_capi 全绿 |
| 4 | `feat(RISK-027): vendor Random123, rewrite Philox4x64 wrapper, GPU/CPU block alignment + KAT` | §5 全部 + README:40 + USAGE_GUIDE + ARCHITECTURE | 三目录全绿 + 性能前后记录 |

回滚单元：四项各自独立；R-027 的 vendor 目录与 rng.hpp 同提交，可整体 revert（DESIGN §8）。
用户规约：**每个提交前 review 排除幻觉**（diff 逐行对 IMPLEMENTATION 与 DESIGN 意图），再提交。

---

## §8 验收映射（A1-A10 → 实施证据）

| # | 标准（DESIGN §9） | 实施节 | 证据命令 |
|---|---|---|---|
| A1 | 钉 commit | §2 | `git -C build\_deps\googletest-src rev-parse HEAD` |
| A2 | 两轴模型 | §3.1/§3.3 | `sync_version.py --check` + test_version.py |
| A3 | README 事实修正 | §3.2 I + §5.7 | grep `1412` / `15 GPU MC` 零命中 |
| A4 | 覆盖率 artifact | §3.2 H | CI 下载 `coverage-report` |
| A5 | C API 测试实际运行 | §3.2 G | c-api job 日志含 `JsonMin*` / `ToolRegistry.*` |
| A6 | 8 MiB+1 拒 + 节点边界 | §4 | P1 用例全绿 |
| A7 | KAT + uniform 位精确 | §5.6 | P2 用例全绿 |
| A8 | 正态 1e-15 | §5.6 | P3 用例全绿 |
| A9 | 全量回归（允许重选 seed） | §5.8 | 三目录 ctest 全绿 + 处置记录 |
| A10 | 性能前后数据 | §5.8 | benchmark_philox 两次输出 |

---

## §9 实施期勘误与设计偏差记录

| # | 项 | 说明 |
|---|---|---|
| E-1 | KAT 夹具 .txt → .hpp | DESIGN §6.5/§7 原为 .txt+运行时解析；改生成头（grunfeld 先例），列序解析移入一次性脚本。已回写 DESIGN（§11 D-4） |
| E-2 | P3 红性质 | DESIGN §7 原记"失败（14.35% 末位不同）"；实测 ≈1e-16 相对差 < 1e-15 容差，且实施前探针无法链接 → 定性为回归锁。已回写 DESIGN §7 |
| E-3 | ci.yml :125 修复保留 `env.CMAKE_BUILD_TYPE` | DESIGN §4.5 片段字面 `Release`；语义等价，与 env 同源 |
| E-4 | benchmark_philox 路径 | DESIGN §6.7 未写路径；实测在 `tests/unit/performance/`（非 benchmarks/ 子目录） |
| E-5 | GPU 探针接口 | DESIGN §6.5 仅述用例；具体导出（rng_device_probe.hpp 三函数）为实施细化，公共契约不变 |
| E-6 | C-13 验证方式 | URL 解包无 `.git`（`git rev-parse` 不可用，实测落回父仓库 HEAD）→ 改为"归档 URL 以 commit SHA 为路径 + 解包内容版本指纹（GOOGLETEST_VERSION 1.14.0）" |
| E-7 | KAT 夹具来源路径 | DESIGN §6.5/IMPLEMENTATION 原写上游 `examples/kat_vectors`；**实际在 `tests/kat_vectors`**（examples 下无该文件，404 实测） |
| E-8 | philox4x64 10 轮向量条数 | 设计预期"抽取 20~40 行"；上游**仅 3 行**（已验证覆盖 0/0、全 f、pi-digits 三类边界），`kKatCount = 3`，无需人工补充 |
| E-9 | `build/` 实际配置 | IMPLEMENTATION §0.1 假设 `build/` 为默认（无 C_API）；实测 **C_API=ON**（与 `build_capi` 重复）。基线按实际配置记录（2553 项）；不改变门禁语义 |
| E-10 | Random123 nvcc 宏守卫 | 见 DESIGN §11 D-6：`gpu_mc.cu` 与探针段必须在 include 前定义 device pass 守卫，否则 device 代码调用 host-only `_umul128` 编译失败 |
| E-11 | CMake 4.4 `gtest_discover_tests` JSON 临时文件竞争 | 全新目录高并行构建（`--parallel 8`）时，`GoogleTestAddTests` 的 `cmake_test_discovery_<hash>.json` 文件名 hash 不含 target 名 → 多 target 争抢同一文件，POST_BUILD 读到编译日志而 JSON 解析失败（`MSB3073` + `ParseTestList.cmake:96`）。**既有基建缺陷，与本次改动无关**；串行构建（`--parallel 1`）或低并行度可规避。CI 用 `--parallel 4`，风险存在但非本次范围 |
| E-12 | CTest 注册表代际污染（build/ vs 全新目录） | 事实：① exe 层面两目录**完全一致**（逐 exe `--gtest_list_tests`：各 157 exe / **2590 用例**）；② `build/` 注册 2573 项，其中**248 组同名重复**（多代发现文件共存：3032 个 `*.cmake` vs 干净目录 1406 个）→ 唯一注册仅 2325，反漏 265 个真实用例（如 `CGMYCharactericFunction.*` 9 个全缺）；③ 干净目录 `tmp/verify_capi`（C_API=ON）经"删除 exe + 发现文件 → 串行重链"后注册 **2590**，与 exe 用例数逐一相等、零重复。判定：**R-032 回归证据以干净目录 2590/2590 全绿为准**（`tmp/ctest_r032_final.log`）；`build/` 的 2573/2573 亦全绿，但既跑重复又漏用例，属既有目录污染（同 E-9/E-11 族），**后续 C-91 三目录回归必须用干净目录**，否则数字不可信 |

---

## §10 Step 6 实施 Review 记录（2026-10-21）

**方法**：本文档全部 `文件:行` 引用对当前工作树逐条重读（`Read`/`Grep -n` 双口径）；全部代码草稿做静态正确性推演（编译可行性、越界/进制/边界算术、生命周期）；与 DESIGN §1-§9、RESEARCH U-1~U-5 逐项对齐。

**核对通过（V 清单）**

| # | 结论 | 证据 |
|---|---|---|
| V-1 | 全部行号锚点属实：README :21/:40/:107/:121（grep 复核，非偏移）、USAGE_GUIDE :183-196/:693、ARCHITECTURE :68/:91、config.hpp :4-8、CMakeLists :2/:25-28/:61-68/:88-90/:110、ci.yml :59/:110/:125/:127-131、python/CMakeLists :10-11/:27-31/:69、cpphub_python.cpp :2043、pyproject :24、rng.hpp :10-13/:18-21/:23-102/:104-136、gpu_mc.cu :4-8/:19-31/:42-45/:47-55/:57-77/:81-87/:90-92/:133-141、json_min :27/:157/:174/:192/:473-481、registry :50/:52、tests/CMakeLists :22/:203-222/:336-337、gpu/CMakeLists :27-34/:36-39/:72-75/:93-114/:102、test_rng :93-101 | 逐条重读 |
| V-2 | `PHILOX_M4x64_*` / `PHILOX_W32_*` / `d_PHILOX_*` **无库外消费者**（仅 rng.hpp 内部与 gpu_mc.cu 自用）→ 删除安全 | grep 35 行，无第三方引用 |
| V-3 | `gpu_mulhi`/`gpu_mullo` 仅被 `philox4x64_10` 调用 → 随其删除 | gpu_mc.cu 内 grep |
| V-4 | `json_min::parse` 每次调用新建 `detail::Parser`（:473-481）→ `nodes_` 成员无需跨调用重置 | 读码 |
| V-5 | test_c_api_tools 仅链 `cpphub_c_api`，而该库 PUBLIC→`cpphub_timeseries_mat`（Eigen 路径 + `EIGEN_MPL2_ONLY` 传递）→ 测试内 include registry 合法，无额外 CMake 改动 | src/c_api/CMakeLists.txt:9-13 |
| V-6 | `.gitignore`：`build*/` 覆盖 build_capi/build_cuda/build_verify；`tmp/` 已忽略；`third_party/random123/`、`tests/fixtures/rng/`、`VERSION` 不被忽略（可入库）；`*.zip` 忽略（U-8 下载物不入库） | .gitignore 逐段 |
| V-7 | rng.hpp 消费方**恰好 42 个头文件**（grep 76 文件 = 42 headers + 34 tests/sources）→ DESIGN §6.8 计数准确，非 42 幻觉 | grep count 后目视分类 |
| V-8 | test_gpu_mc 现为 **19** 个 TEST（:39-358），README:107 的 "15 GPU MC tests" 双重失真（数量 + CPU 可跑）→ DESIGN 替换文案回避数字，正确 | test_gpu_mc.cpp 通读 |
| V-9 | 根 CMake :110 **无条件** `add_subdirectory(src/performance/gpu)` → CUDA off 时 stub 分支生效、`test_gpu_mc` 不构建（:113-114 `SKIPPED (no CUDA)`）→ README:107 替换文案属实 | 根 CMake :104-110 + gpu CMake :93/:113 |
| V-10 | `discard(n)` 新实现推演正确：`discard(50)` → 12 块 + 2 字内消费，与逐字推进 50 次的流一致；块内部分 discard（avail 分支）与 256 位进位（`add_to_ctr`/`bump_ctr`）均正确 | 手算 + 既有 `PhiloxDiscard` 用例 |
| V-11 | 节点预算边界算术：N 元素数组 = 1+N 节点；`kMaxParseNodes-1` 元素 = 恰 1e6（不触发 `>`）；`kMaxParseNodes+2` 元素 = 1'000'003（在 1'000'001 次调用失败） | 手算 |
| V-12 | 探针 kernel 的 grid 向上取整 ⇒ 原草稿缺越界守卫（已修）；`KatVectorIn` 先用后定义、命名空间重复嵌套均为真编译错误（已修） | 静态推演 |

**勘误（已回写本文档）**

| # | 位置 | 原稿 | 问题 | 修正 |
|---|---|---|---|---|
| R-1 | §4.4 | 追加块内再写 `namespace jm = cpphub::v1::json_min;` | 与 [test_json_min.cpp:18](../../tests/unit/c_api/test_json_min.cpp#L18) 重复声明 namespace alias → **编译错误** | 删除该行，改注释提醒 |
| R-2 | §5.5 | `KatVectorIn` 在 kernel 之后定义；文件末尾重开 `namespace cpphub { inline namespace v1 {`；kernel 无越界守卫 | 三处均为真编译/运行错误（错层嵌套、先用后定义、grid 取整越界写） | 复用 `gpu_probe::KatVector`；只开 `namespace gpu_probe`；补 `path >= n` 守卫 + `cudaGetLastError()` |
| R-3 | §3.3 | `HEADER.write_text(..., newline="\n")` | `write_text` 的 `newline` 参数需 **Python 3.10**，项目 `requires-python >= 3.9` | 改 `HEADER.open("w", ..., newline="\n")` |
| R-4 | §3.6 | 无 CI `--check` 步骤 | DESIGN §4.2 B 要求"CI `--check` 防漂移"，§4.4 YAML 未包含 → 门禁缺失 | 新增 build-test（Linux）步骤 `python3 tools/sync_version.py --check` |
| R-5 | §3.2 H | 照抄 DESIGN `python -m pip install gcovr` | ubuntu-latest 无 `python` 别名 → coverage job 挂 | 改 `python3 -m pip` |
| R-6 | §3.2 H | artifact 仅 `coverage.html` | `--html-details` 产出多文件，明细链接失效 | 改 `coverage*.html` |
| R-7 | §5.7 | "counter 2^256 / 4×u64 / 分层承诺" 概括 | USAGE_GUIDE 现文 :190 周期 `2^128`、:191 状态 `counter_[2]+key_[2]`、:192 输出、:196 "位精确一致" 会与新实现**自相矛盾** | 列为逐行改（:190/:191/:192/:196） |
| R-8 | §5.6 | test_gpu_mc 未加 `#include "cpphub/core/rng.hpp"` | 新用例用 CPU `Philox4x64`/`box_muller`，test_gpu_mc 现无此包含 → 编译失败 | 补齐 include |
| R-9 | §5.6 | `auto [zc1, zc2] = box_muller(...)`，`zc2` 未用 | MSVC /W3 噪音（且无意义） | 改 `box_muller(...).first` |
| R-10 | §1.1 | `https://github.com/AMitherd/random123` | **无 E1 依据的仓库 URL（幻觉风险）** | 改为待 U-8 页面级确认，候选仅作提示 |
| R-11 | §5.8 | "写入 DEVELOPMENT_LOG" 未给路径；未说明 benchmark 平台标签语境 | 落地歧义；标签（`_umul128`/`__uint128_t`）在 vendor 后语义变化 | 明确 [DEVELOPMENT_LOG.md](../../docs/DEVELOPMENT_LOG.md) + 标签语境说明 + 理论 ~4× 量级 |
| R-12 | §3.2 H | "[:127] 前" 表述含糊 | 插入点在 c-api job 结束与 python job 注释块之间 | §3.6.3 精确化（:127-131 前） |
| R-13 | §4.4 | "~120 MB 量级" | 估算口径不明（未含 vector 扩容） | 改为 sizeof(Value) 96–112 B ⇒ 峰值 150–200 MB |

**非勘误（复核后确认原稿正确，避免过度修正）**：README :107/:121 行号经 `grep -n` 复核属实（Read 偏移读数一度疑似 106/120，以 grep 为准）；DESIGN §4.3 文案与 gpu CMake `SKIPPED (no CUDA)` 逐字一致；DESIGN §6.8 "42 个 include 方" 精确；`.gitignore`、`parse()` 生命周期、`cpphub_c_api` 的 Eigen 传递均无需改动。

**结论**：无方向性/裁决性错误；13 处勘误为编译错误（R-1/R-2/R-8）、运行错误（R-2 越界）、工具链兼容（R-3）、CI 可执行性（R-4/R-5/R-6）、文档自洽（R-7/R-11）、幻觉清除（R-10）与估算口径（R-13），不改变任何设计裁决与验收口径。**可进入 Step 7（CHECKLIST）**。

---

**实施人**: Scott（鹏） · **依据**: DESIGN（Step 4 review 后）+ 本文档 §1 探针结论 · **规约**: TDD 红→绿、禁放宽容差、性能改动附前后数据、提交前 review
