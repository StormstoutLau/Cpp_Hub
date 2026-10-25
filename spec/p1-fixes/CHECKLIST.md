# 审查验收 Checklist：P1 整改（p1-fixes）

---

id: P1-FIXES-CHECKLIST
type: acceptance
version: 1.0
status: audited
date: 2026-10-21
depends: [P1-FIXES-IMPLEMENTATION, P1-FIXES-DESIGN, P1-FIXES-RESEARCH]
upstream: EXTERNAL_REVIEW_2026-10-20

---

> **Feature**: p1-fixes（RISK-027 / 032 / 033 / 034）
> **状态**: **验收通过**（§2-§8 全部 ✅，2026-10-22 签收；§1 一致性核对已于 2026-10-22 完成；Step 10 ADD 审计修复已回写，见 [AUDIT.md](./AUDIT.md)）
> **Spec 步骤**: Step 7（本清单）→ Step 8（§1 文档一致性）→ Step 9（TDD 实施）→ Step 10（ADD 审计）
> **基于实施**: [IMPLEMENTATION.md](./IMPLEMENTATION.md)（Step 6 review 后）· **基于设计**: [DESIGN.md](./DESIGN.md)（Step 4 review 后）· **调研**: [RESEARCH.md](./RESEARCH.md)
> **计数基准**：G 系列 **6**（§2）+ C 系列 **78**（§3 6 / §4 19 / §5 12 / §6 28 / §7 6 / §8 7）= **84** 项；口径 = 逐项编号计数（`grep -c '^| [GC]-'`，非事后推算）

---

## 0. 使用方式与判定口径

1. **状态标记**：`☐` 未执行 / `✅` 通过 / `❌` 失败 / `⏭` 不适用（须在证据列写理由）。
2. **证据必须可重放**：命令 + 关键输出（行数、计数、exit code）写入"证据"列，不接受"已验证"式空口记录。
3. **禁止项（硬约束，违反即 ❌）**：
   - 禁止放宽容差换取绿；统计窗口红只允许**重选 seed**，并在提交说明写明"seed 变更，非阈值放宽"。
   - 禁止修改 `tests/fixtures/**` 既有数据（本 feature 不重算基线）。
   - 禁止把 KAT 夹具退回运行时解析 `.txt`（D-4）。
   - 禁止在 `gpu_mc.cu` 恢复自写 Philox（`d_PHILOX_*` / `philox4x64_10`）。
4. **三构建目录**（IMPLEMENTATION §0.1）：`build/`（核心）、`build_capi/`（C ABI）、`build_cuda/`（CUDA）。Windows 多配置生成器一律带 `-C Release`。

---

## 1. 文档一致性（Step 8，2026-10-22 完成）

> **执行记录**：同一 agent 体系单视角自查（非独立复审）；机械命令 = `grep -c` 与逐条 `Read`。发现 **7 处不一致（1×P2 + 6×P3）**，全部在本轮闭环（见 §1.2）。

### 1.1 核对结果

| # | 核对维度 | 方法 | 状态 | 证据 |
|---|---|---|---|---|
| D-01 | RESEARCH U-1~U-5 裁决在 DESIGN 落地 | 逐条对照 [RESEARCH §6.4](./RESEARCH.md) 与 DESIGN | ✅ | U-1→§6.3 块对齐 + §6.5 分层承诺 + §4.3 README；U-2→§4.1/§4.2 + §11 D-1（`configure_file` 收敛）；U-3→§4.4（无阈值）；U-4→§5.1-§5.4（入口 strlen + 节点预算 + INVALID_ARG）；U-5→§6.1-§6.4（vendor/封装/GPU 共用/KAT）。**5/5** |
| D-02 | DESIGN §9 A1-A10 ↔ 本清单 §9 一一对应 | 双侧 `grep -c '^\| A\d+ \|'` | ✅ | DESIGN 10 / IMPLEMENTATION 10 / CHECKLIST 10（30 行总命中） |
| D-03 | IMPLEMENTATION 每节可追溯到 DESIGN | 章节对照表 | ✅ | 本轮新增 [IMPLEMENTATION §0.4](./IMPLEMENTATION.md) 对照表（§0-§9 全覆盖；§10 为流程产物，无 DESIGN 对应） |
| D-04 | 四类测试名三文档一致 | `grep` 归一化比对 | ✅ | 9 个新测试（RNG 4 + JSON/CApi 5）三文档均命中。**归一化规则**：ctest 名 = `Suite.Name`；IMPLEMENTATION 代码草稿的 `TEST(Suite, Name)` 视为同一测试（原"逐字一致"口径不成立，见 S8-4） |
| D-05 | 无失效**生效用法**残留 | `grep` | ✅ | `write_text(content` 0；`kat_philox4x64_10.txt` 在 IMPLEMENTATION 0；`KatVectorIn` 仅 3 处（§5.5 禁止说明 + §10 V-12/R-2 勘误）；探针段重开 `namespace cpphub` 仅 1 处禁止说明。勘误/禁止语境引用不计（原"零命中"口径误伤，见 S8-5） |

### 1.2 Step 8 发现与闭环

| # | 严重性 | 发现 | 证据 | 处置 |
|---|---|---|---|---|
| S8-1 | P3 | DESIGN §6.2 把根 CMake include 区间写成 `:25-29`，而 :29 实为 `target_compile_features` | 读码 :25-28 为 include 列表 | 改 `:25-28`（DESIGN §6.2 已回写） |
| S8-2 | P3 | DESIGN §4.4 coverage job 与 IMPLEMENTATION §3.6.3 不一致（`python` vs `python3`；`coverage.html` vs `coverage*.html`） | 两侧对照 | 回写 DESIGN §4.4（`python3` + glob）；并补一行指向 §3.6.2 的 `--check` 门禁 |
| S8-3 | **P2** | **计数错误**：IMPLEMENTATION §5.6/§6 与 CHECKLIST C-71 写"test_rng 既有 8 用例"，实测 **9**（RESEARCH §3.1 E1-6 早已记为 9） | `grep -c '^TEST\(' tests/unit/core/test_rng.cpp` = 9 | 三处改 9（IMPLEMENTATION ×2、CHECKLIST ×1） |
| S8-4 | P3 | CHECKLIST D-04 口径"逐字一致"不可执行：IMPLEMENTATION 用 gtest 形 `TEST(Suite, Name)` | 机械 grep 显示 IMPLEMENTATION 对 5 个 JSON 测试名为 0 命中 | D-04 改为归一化规则（§1.1） |
| S8-5 | P3 | CHECKLIST D-05 口径"零命中"会误伤勘误说明，且拼写为 `WriteText` | grep 命中 4 行（均在禁止/勘误语境） | D-05 改为"零生效用法"+ 修正大小写 |
| S8-6 | P3 | CHECKLIST C-79 只写 `RngKat`，缺完整测试名，D-04 无法通过 | 对照 DESIGN §6.5 全名 | C-79 补 `RngKat.Philox4x64MatchesKnownAnswers` |
| S8-7 | P3 | IMPLEMENTATION §1.1 上游 URL 为占位"待确认"，而 RESEARCH §7.1 已有 E1 锚点 `DEShawResearch/random123` | RESEARCH §7 参考文献 1/3 | 以 RESEARCH 为准填入 URL，删除占位 |

> 历史记录不入勘误：DESIGN §12（Step 4）仍保留当时的 `:25-29` 表述，属已签核的审计快照；以 §6.2 为准。

---

## 2. 实施前置门禁（G 系列）

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| G-01 | U-7 发布号已确认 | 向 Scott 确认；写入 `python/pyproject.toml:24` | 取值确定（默认 1.9.0）且与提交说明一致 | **Scott 裁定 1.9.0（2026-10-22）** | ✅ |
| G-02 | U-8 vendor 身份可审计 | 读 `third_party/random123/UPSTREAM.md` | 含**页面级核实**的仓库 URL + 40 位 commit SHA + 拷贝日期 + 范围 + 许可证 | commit `9545ff6413f258be2f04c1d319d99aaef7521150`（GitHub API 直连确认）；codeload 归档 430 644 B；**27 文件** + LICENSE 落地；UPSTREAM.md 已写 | ✅ |
| G-03 | U-6 KAT 列序已确认 | 打开上游 `kat_vectors` 首行，与抽取脚本输出比对 | 脚本产出的 `ctr/key/answer` 三组数值与上游行**逐字一致** | 路径 `tests/kat_vectors`（非 examples，E-7）；列序 `nameNxW R CTR[4] KEY[2] EXPECTED[4]`；philox4x64 10 轮 3 行（E-8）；夹具 `tests/fixtures/rng/kat_philox4x64_10.hpp` 已生成（kKatCount=3），第 1 行答案与 nvcc device 输出逐字一致 | ✅ |
| G-04 | U-10 探针通过（Random123 在 nvcc 下可用） | 运行 `tmp/probe_r123.cu`（IMPLEMENTATION §1.2） | ① nvcc 编译通过；② `R123_USE_PHILOX_64BIT` 非 0；③ 设备 4 字 == host 4 字 | host：`PHILOX_64BIT=1` + 输出=KAT ✓；device **首轮失败**（`mulhilo64` host/device 错配）→ 加 device pass 守卫后 `nvcc -arch=sm_89` 通过，device 输出 `16554d9e…7ba23b` = KAT 第 1 行逐字一致（D-6 / E-10） | ✅ |
| G-05 | 基线全绿（改动前） | `ctest --test-dir build -C Release -j 8 --output-on-failure`；同 build_capi / build_cuda | 三目录 0 失败 | build **2605/2605**（重配清 4 项 NOT_BUILT 残留；残留源 = 旧 configure 注册，E-9）；build_capi **2587/2587**；build_cuda **2548/2548**（含 `gpu_mc` 通过） | ✅ |
| G-06 | 离线 gtest 逃生口可用（网络受限时） | `-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=<本地目录>` 配置成功 | configure 无下载错误 | 本地 `build/_deps/googletest-src`（GOOGLETEST_VERSION=1.14.0）可用为覆盖源；**且已裁决改用 URL@commit（codeload 直连可达 200），主路径不再依赖覆盖**（DESIGN §11 D-5） | ✅ |

> G-01~G-06 任一未通过 → 停止实施，不得进入 §3。

---

## 3. RISK-033 — 供应链钉死

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-11 | 无 tag 归档 URL 依赖 | `Select-String -Path CMakeLists.txt -Pattern 'archive/refs/tags'` | 0 命中 | **0 命中**（2026-10-22 实测） | ✅ |
| C-12 | 以 commit 为不可变身份 | 读 `CMakeLists.txt` googletest 块 | URL 形如 `https://codeload.github.com/google/googletest/tar.gz/<40 位 SHA>`（URL@commit 形态，D-5）；无 tag 名、无 `GIT_REPOSITORY` | 内容锚（行号随提交序列漂移, 2026-10-22 复测 URL 位于 :72）：URL = `.../tar.gz/f8d7d77c06936315286eb55f8de22cd23c188571` + `DOWNLOAD_EXTRACT_TIMESTAMP TRUE`。**归档 SHA-256 留存（ADD 审计 OBS-1）**: `7FF5DB23DE232A39CBB5C9F5143C355885E30AC596161A6B9FC50C4538BFBF01`（868 904 B, codeload 实测 2026-10-22） | ✅ |
| C-13 | 全新 configure 解析到同一份源码 | 全新目录 `cmake -S . -B build_verify`；读 `build_verify\_deps\googletest-src\CMakeLists.txt` | configure 无下载错误；`GOOGLETEST_VERSION 1.14.0`；且 URL 路径 SHA = `f8d7d77c06936315286eb55f8de22cd23c188571`（URL 解包无 `.git`，身份证据 = SHA 路径 + 版本指纹，E-6） | configure exit 0；`urlinfo.txt`: `url(s)=https://codeload.github.com/google/googletest/tar.gz/f8d7d77c…8571`；`set(GOOGLETEST_VERSION 1.14.0)` | ✅ |
| C-14 | cache key 同步且计数正确 | `Select-String -Path .github/workflows/ci.yml -Pattern 'gtest-1.14.0'`；再查 `gtest-f8d7d77c` | 旧 key 0 命中；新 key **2** 命中（build-test 与 c-api 两处, 2026-10-22 复测为 :59/:115） | 旧 **0** / 新 **2**（2026-10-22 实测） | ✅ |
| C-15 | python job 未被误改 | 读 ci.yml :132-177 | 仍为 pip cache，无 `_deps` | python job 段 `_deps` 命中 **0** | ✅ |
| C-16 | 清缓存后测试全绿 | `ctest --test-dir build_verify -C Release -j 8` | 0 失败；验完删除 `build_verify` | **2536/2536**（2026-10-22；含 2 项 NOT_BUILT 残留重建后全绿）；`build_verify` 已删除 | ✅ |

---

## 4. RISK-032 — 版本两轴 + 文档事实 + CI

### 4.1 版本单一源

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-21 | `VERSION` 存在且格式合法 | `Get-Content VERSION` | 单行 `1.0.0`（LF 结尾） | 内容 | ✅ 6 字节 `31 2E 30 2E 30 0A`，仅 LF 无 CR |
| C-22 | `version.hpp` 为生成物且已入库 | `git ls-files include/cpphub/core/version.hpp`；读文件头 | 已被跟踪；含"请勿手改"头注释 + 4 个宏 | `git ls-files` 输出 | ✅ 由 `tools/sync_version.py` 生成（356 B，CR 计数 0）；入库随提交 2 |
| C-23 | `--check` 门禁生效 | `python tools/sync_version.py --check; echo $LASTEXITCODE` | exit 0；篡改 version.hpp 后重跑 exit 1 且 stderr 指明不一致 | 两次 exit code | ✅ 原态 exit 0（"一致: 1.0.0"）；篡改 `"1.0.0"`→`"0.9.9"` 后 exit 1 + stderr "不一致 — 运行:"；恢复后 exit 0 |
| C-24 | config.hpp 不再手写版本 | `Select-String -Path include/cpphub/core/config.hpp -Pattern 'define CPPHUB_VERSION'` | 仅 `#include "cpphub/core/version.hpp"`，宏定义 0 命中 | 计数 0 + include 行 | ✅ 宏定义 0 命中；:4-5 为注释 + `#include "cpphub/core/version.hpp"` |
| C-25 | 根 CMake 从 VERSION 取版本 | `cmake -S . -B build` 后读 `build/CMakeCache.txt` 或 `cmake -LA` | `CppHub` 项目版本 = 1.0.0（无第二处硬编码） | configure 输出 | ✅ configure exit 0；`CMAKE_PROJECT_VERSION:STATIC=1.0.0`（MAJOR=1 / MINOR=0 / PATCH=0） |
| C-26 | Python 轴注入正确 | 读 `python/CMakeLists.txt` | 有 `SKBUILD_PROJECT_VERSION` 分支 + `../VERSION` 回退；`_core` 带 `CPPHUB_PY_VERSION` 编译定义 | 两段代码 | ✅ :12-20 分支；:80-81 `target_compile_definitions(_core PRIVATE CPPHUB_PY_VERSION=...)`。回退分支实测：standalone configure exit 0 且 `CMAKE_PROJECT_VERSION=1.0.0`（来自 `../VERSION`） |
| C-27 | `__version__` 与 pyproject 一致 | `pip install "./python[test]"`；`python -m pytest python/tests/test_version.py -v` | 1 passed；`__version__` == pyproject `version` | pytest 输出 | ✅ wheel 重装成功（1.9.0）；`__version__ = 1.9.0`；`test_version_matches_pyproject PASSED`。python 全量 `111 passed, 1 skipped` |
| C-28 | C-27 为真红→绿 | 先在旧 wheel（`__version__="1.4.0"`）上跑该测试 | 改 pyproject 后、改代码前：**FAILED**；改 3.2 E/F 后重装：**PASSED** | 两次输出 | ✅ 红：`AssertionError: assert '1.4.0' == '1.9.0'`（旧 wheel）；绿：重装 1.9.0 后 `PASSED` |
| C-29 | 无第三处版本硬编码 | `Select-String -Path python/src/cpphub_python.cpp -Pattern '1\.4\.0'`；全仓 `src/` 查 `1\.4\.0` | 0 命中（`docs/` 历史记录除外） | 计数 0 | ✅ `python/src/cpphub_python.cpp` 0 命中；`src/` 全树 0 命中 |

### 4.2 README 事实修正

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-31 | 冻结总数清零 | `Select-String -Path README.md -Pattern '1412'` | 0 命中 | 计数 0 | ✅ 0 命中；:27 改为"全量测试通过（见 CI）；位精确复现限**同编译器族内 + 同 RNG 版本 + 同编译选项**" |
| C-32 | GPU 测试数失真清除 | `Select-String -Path README.md -Pattern '15 GPU MC'` | 0 命中；:107 改为"stub 编译 + `test_gpu_mc` skipped"表述 | 计数 0 + :107 全文 | ✅ 0 命中；:113 "a CPU stub is compiled and `test_gpu_mc` is **skipped** (`test_gpu_mc: SKIPPED (no CUDA)`); GPU entry points throw `std::runtime_error`" |
| C-33 | RNG 承诺分层（:46） | 读 README :46 | uniform bit-exact + normal 1e-15 双承诺 | ✅ README:46 "Philox4x64-10 (Random123): **uniform stream bit-exact** with CPU; normal samples agree to 1e-15"（原计划 :40，实际行漂移至 :46） |
| C-34 | 跨平台限定（:121） | 读 README :121 | 注明"同 RNG 版本 + 同编译选项" | 行全文 | ✅ :127 "Cross-platform float determinism \| bit-exact \| ✅ (MSVC + 2× GCC；同 RNG 版本 + 同编译选项)" |
| C-35 | 顶部版本模型说明 | 读 README 顶部 | 3 行：C++ API 版本见 `VERSION` / Python 发行版本见 `pyproject.toml` / tag 独立 | 段落 | ✅ :9-13 "版本模型（两轴，彼此独立）"三行齐全 |

### 4.3 CI 可执行性

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-36 | C API 过滤缺口已修 | `Select-String -Path .github/workflows/ci.yml -Pattern '-R "CApi"'` | 0 命中；Run C API Tests 步骤无 `-R` | 计数 0 + 步骤全文 | ✅ 0 命中；:128-130 `ctest --build-config ${{ env.CMAKE_BUILD_TYPE }} --output-on-failure --parallel 4`（无 `-R`）。干净目录实测：`CApi*` 26 个 + `JsonMin\|ToolRegistry` 28 个（合计 54）在无过滤下全部执行 |
| C-37 | 版本漂移门禁进 CI | 查 ci.yml build-test job | 含 `python3 tools/sync_version.py --check`（`if: runner.os == 'Linux'`） | 步骤全文 | ✅ :74-77 build-test job 内 "Check version header (RISK-032)"，`if: runner.os == 'Linux'`，`python3 tools/sync_version.py --check` |
| C-38 | coverage job 存在且可跑 | 查 ci.yml | Debug + `--coverage` + `--filter 'include/cpphub/'` + **`python3 -m pip install gcovr`** + artifact `coverage*.html` / `coverage.xml` | job 全文 | ✅ coverage job 现位于 ci.yml :136-164（2026-10-22 复测）：Debug、`--coverage -O0 -g`、`gcovr --filter 'include/cpphub/' --exclude 'third_party'`、`python3 -m pip install gcovr`、artifact `coverage.xml` + `coverage*.html` |
| C-39 | YAML 语法有效 | `python -c "import yaml; yaml.safe_load(open('.github/workflows/ci.yml',encoding='utf-8')); print('ok')"`（无 PyYAML 时以推送后 Actions 解析为准） | ok / Actions 不报 syntax error | 输出 | ✅ `yaml ok`（PyYAML，Python 3.12） |
| C-40 | A5 实证：JsonMin/ToolRegistry 在 CI 真跑 | 看 c-api job 日志 | 日志出现 `JsonMinBudget.*`、`ToolRegistry.*`、`CApiTools.*` | 日志片段 | ✅ run 39755779788 C ABI job（ubuntu #119747299808 / windows #119747299759）：`(JsonMinBudget\|ToolRegistry\|CApiTools)\..*Passed` = 30 行/平台，双平台 2600/2600 全绿；本审计新增 3 用例（RejectsFirstNodeOverBudget #1341 / ObjectKeysNotCountedInNodeBudget #1342 / AcceptsArgsJsonAtExactLimit #1366）均在列且 Passed。注：同 run 其余 job 失败（guard 自匹配等 4 簇）另行修复，见 DEVELOPMENT_LOG 2026-10-22 CI 条目 |

---

## 5. RISK-034 — `json_min` 输入预算

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-41 | 节点预算常量与语义注释 | 读 `json_min.hpp`（`kMaxParseDepth` 之后） | `kMaxParseNodes = 1'000'000`；注释写明"每次 parse_value 计 1，键不计入" | 代码片段 | ✅ :29-34 常量 + 5 行注释（含"键不计入"与合法上界推导） |
| C-42 | 入口长度常量 | 读 `cpphub_tool_registry.hpp`（`kMaxColumns` 之后） | `kMaxArgsJsonBytes = 8u*1024u*1024u` + 策略层注释 | 代码片段 | ✅ :53-55 `inline constexpr std::size_t kMaxArgsJsonBytes = 8u * 1024u * 1024u;` + 策略层注释 |
| C-43 | 计数实现位置正确 | 读 `Parser::parse_value` 入口 | `if (++nodes_ > kMaxParseNodes) return fail("maximum node count exceeded");` | 代码片段 | ✅ :200-204 `parse_value` 入口（switch 前）；:182 成员 `std::size_t nodes_ = 0;` |
| C-44 | 入口检查**前置于**工具查找与解析 | 读 `cpphub_v1_tool_call` | strlen 检查位于 args 判空之后、`find_tool` 之前；返回 `CPPHUB_ERR_INVALID_ARG` | 代码片段 | ✅ c_api.cpp :308-315；`std::strlen` + `kMaxArgsJsonBytes` 比较 + `set_error("args_json exceeds limit (...)")` + 返回 `CPPHUB_ERR_INVALID_ARG` |
| C-45 | 边界算术正确（接受侧） | 运行 `JsonMinBudget.AcceptsAtBoundary` | 999'999 元素（=1e6 节点）parse 成功 | 用例通过 | ✅ PASSED（142 ms）；`as_array().size() == kMaxParseNodes - 1` |
| C-46 | 边界算术正确（拒绝侧） | 运行 `JsonMinBudget.RejectsWideArrayOverBudget` | 1'000'002 元素拒绝，error 含 `node` | 用例通过 + 错误串 | ✅ PASSED（144 ms）；红态曾 `Actual: true`（无上限时解析成功），绿态拒绝且 err 含 "node" |
| C-47 | 两类限制互不遮蔽 | 运行 `JsonMinBudget.DeepNestingStillLimitedByDepth` | 深嵌套仍报 `nesting depth` | 用例通过 | ✅ PASSED（0 ms）；err 含 "nesting depth" |
| C-48 | 8 MiB+1 被拒且文案正确 | 运行 `CApiTools.RejectsOversizedArgsJson` | `-1`；`get_last_error` 含 `exceeds limit` | 用例通过 | ✅ PASSED（4 ms）；返回 `CPPHUB_ERR_INVALID_ARG`，err 含 "exceeds limit" |
| C-49 | 检查前置可证 | 运行 `CApiTools.OversizedCheckPrecedesParse` | 超长且非法 → 报 `exceeds limit`，**不含** `invalid args_json` | 用例通过 | ✅ PASSED（4 ms）；红态曾报 `invalid args_json: invalid number`（证明未前置），绿态仅 `exceeds limit` |
| C-50 | 红→绿性质成立 | 在改动前跑 §5 用例 | C-46/C-48/C-49 为 **FAILED**；C-45/C-47 为 PASSED（回归锁） | 改动前后两次 ctest 输出 | ✅ 仅加常量未加逻辑时实测：`JsonMinBudget` 2 PASSED / 1 FAILED（AcceptsAtBoundary、DeepNesting PASSED）；`CApiTools` 0 PASSED / 2 FAILED；实施后 `test_json_min` 17/17、`test_c_api_tools` 10/10 |
| C-51 | 无新增错误码 / ABI 不破 | 读 `cpphub_c_api.h`；`cpphub_v1_get_abi_version()` | 错误码集合不变；ABI 版本仍 = 1 | 头文件 diff + 用例 | ✅ `cpphub_c_api.h` 未改动（`git diff` 0 命中）；复用 `CPPHUB_ERR_INVALID_ARG` |
| C-52 | 资源开销可接受 | 记录 `ctest -R "JsonMinBudget"` 用时 | 两个 1e6 节点用例合计 < 10 s（本机实测记录） | 耗时数字 | ✅ 三个用例合计 287 ms（142 + 144 + 0），远低于 10 s |

---

## 6. RISK-027 — vendor Random123 + 块对齐 + KAT

### 6.1 供应链与构建

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-61 | vendor 目录完整 | `Get-ChildItem third_party/random123 -Recurse` | `include/Random123/**` + `LICENSE` + `UPSTREAM.md` | `include/Random123/**` **27** 文件 + LICENSE + UPSTREAM.md（现树实测） | ✅ |
| C-62 | 上游头文件未被修改 | 与 U-8 下载源逐文件比对（hash 或 `diff -r`） | 差异仅限新增 `UPSTREAM.md`/目录布局 | 上游 zip（`tmp/r123_dl/random123.zip`，SHA-256 `1E0C81EE…E01F0`）解包后逐文件比较：**27/27 SHA-256 全等，0 diff**；LICENSE 逐字节相同 | ✅ |
| C-63 | 三处 include 路径补齐 | 读根 CMake :25-28、gpu CMake :36-39、python CMake :27-31 | 各含 `random123/include` | 根 `CMakeLists.txt:30`、`src/performance/gpu/CMakeLists.txt:40`、`python/CMakeLists.txt:40` 均含 | ✅ |
| C-64 | stub 分支未被污染 | 读 gpu CMake :65-88 | 无 random123 路径、无 `.cu` | :67-90 仅 `gpu_mc_cpu_stub.cpp` + `gpu_config.cpp`，无 random123 路径 | ✅ |
| C-65 | 库外无人引用旧常量 | `Get-ChildItem -Recurse -Include *.cpp,*.hpp,*.cu \| Select-String -Pattern 'PHILOX_M4x64\|PHILOX_W32\|d_PHILOX'` | 除 `spec/`、`docs/` 外 **0** 命中 | 0 命中（仅两处 `#include "kat_philox4x64_10.hpp"` 命中文件名子串，非常量） | ✅ |

### 6.2 `rng.hpp` 契约

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-66 | 公共签名不变 | 读 rng.hpp | 类名 `Philox4x64`；`result_type`、`min/max`、`Philox4x64(uint64_t=0,uint64_t=0)`、`operator()`、`discard(uint64_t)`、`next4()->std::array<uint64_t,4>` 全保留 | `rng.hpp:17-83`：类名 :17、result_type/min/max :24-26、构造 :28、operator() :34、next4 :41、discard :46 | ✅ |
| C-67 | 旧实现残留清零 | `Select-String -Path include/cpphub/core/rng.hpp -Pattern 'mulhi\|mullo\|PHILOX_'` | 0 命中；含 `#include <Random123/philox.h>` | 全文 0 命中；`rng.hpp:6` 为 include 行 | ✅ |
| C-68 | `box_muller`/`normal_simd`/`generate_correlated` 未改 | `git diff -U0 include/cpphub/core/rng.hpp` | 三函数区间零 diff | diff 末 hunk 止于 `:81`（类私有段），:85-120 三函数零 hunk | ✅ |
| C-69 | 块对齐语义 | 读 `next4()`；跑 `Rng.PhiloxNext4Consistency` | `next4() == 4×operator()`；用例通过 | `rng.hpp:41-43` 逐字等于 4 次 `operator()`；test_rng 10/10 通过 | ✅ |
| C-70 | discard 语义 | 跑 `Rng.PhiloxDiscard` | 通过（50 → 12 块 + 2 字内） | test_rng 10/10 通过 | ✅ |
| C-71 | 既有 **9** 用例无需改动即绿 | `git diff tests/unit/core/test_rng.cpp`（除新增 KAT 段） | 既有用例零修改；全部通过 | `+17` 行（1 include + `RngKat` 15 行），既有 9 用例零改动；10/10 通过 | ✅ |

### 6.3 GPU 侧

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-72 | 设备端自写 Philox 已删 | `Select-String -Path src/performance/gpu/gpu_mc.cu -Pattern 'd_PHILOX|philox4x64_10|gpu_mulhi|gpu_mullo'` | 0 命中 | 0 命中（gpu_mc.cu 全文） | ✅ |
| C-73 | 内核走 Random123 | 读内核调用点 | `ctr.v[0]=path`、`key.v[0]=seed`、`r.v[0]/r.v[1]` 进 `uint64_to_unit` | `gpu_mc.cu:104-111` 符合 | ✅ |
| C-74 | 探针头存在且仅 CUDA 分支提供实现 | 读 `rng_device_probe.hpp`；`Select-String -Path src/performance/gpu/gpu_mc_cpu_stub.cpp -Pattern 'gpu_probe'` | 头文件 3 个函数声明；stub 0 命中（不做未定义引用） | 头 :20/:23/:26 三声明；stub 0 命中 | ✅ |
| C-75 | 探针 kernel 安全 | 读 gpu_mc.cu 探针段 | 三个 kernel 均有 `path >= n` 守卫；launch 后 `cudaGetLastError()`；空输入早退 | `gpu_mc.cu:503/515/527` 守卫；:553/572/585 `CPPHUB_CUDA_CHECK(cudaGetLastError())`；:547/561/580 空输入早退（2026-10-22 复测） | ✅ |
| C-76 | 探针头登记进 FILE_SET | 读 gpu CMake :27-34 | FILES 含 `rng_device_probe.hpp` | `src/performance/gpu/CMakeLists.txt:34` | ✅ |

### 6.4 KAT 与对齐测试

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-77 | 夹具为生成头（非 .txt） | `Test-Path tests/fixtures/rng/kat_philox4x64_10.hpp`；`Test-Path ...\.txt` | hpp = True；txt = False | hpp=True / txt=False | ✅ |
| C-78 | 夹具内容来自上游 | 抽查 ≥3 行与上游 `tests/kat_vectors` 逐字比对 | 完全一致；`kKatCount = 3`（上游 10 轮全部行） | 上游 `kat_vectors` 3 行逐字一致（`16554d9e…` / `87b092c3…` / `a528f454…`）；`kKatCount=3` | ✅ |
| C-79 | CPU KAT：`RngKat.Philox4x64MatchesKnownAnswers` | `ctest --test-dir build -C Release -R "RngKat"` | 1/1 通过 | build 与 `tmp/verify_capi` 均 Passed | ✅ |
| C-80 | GPU KAT | `ctest --test-dir build_cuda -C Release -R "gpu_mc"`，读 `GpuMC.PhiloxMatchesKnownAnswers` 结果 | 通过（设备端与上游一致） | build_cuda `gpu_mc` 22/22（2026-10-22，RTX 4060） | ✅ |
| C-81 | uniform 流逐字对齐 | 同上，`GpuMC.UniformStreamMatchesCpu` | p=0..9999 共 40'000 字全等 | 22/22 中的 `GpuMC.UniformStreamMatchesCpu` Passed（10000 块×4 字逐字相等） | ✅ |
| C-82 | 正态层级承诺 | 同上，`GpuMC.NormalSamplesWithinTolerance` | 相对差 ≤1e-15（\|z\|<1 时绝对 ≤1e-15） | 22/22 中的 `GpuMC.NormalSamplesWithinTolerance` Passed | ✅ |
| C-83 | P3 定性记录正确 | 读 DESIGN §7 与 IMPLEMENTATION §5.6 | 均写明"回归锁、非红驱动" | DESIGN §11 D-4 / IMPLEMENTATION §5.6 写明 KAT 为回归锚，非红驱动 | ✅ |
| C-84 | 测试 include 路径齐全 | 读 tests/CMakeLists :22 附近、gpu CMake test_gpu_mc 块 | 两处 `target_include_directories(... fixtures/rng)`；test_gpu_mc.cpp 含 `rng.hpp` | `tests/CMakeLists.txt:24`、`src/performance/gpu/CMakeLists.txt:106`；test_gpu_mc.cpp 含 `cpphub/core/rng.hpp` + 探针头 | ✅ |

### 6.5 文档同步

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-85 | USAGE_GUIDE §3.1.3 自洽 | 读 :183-196 | 算法=Random123；周期=2^256；状态=counter[4]+key[2]（384 位）；输出=4×u64/块；承诺=分层 | :185-196 五项全符合 | ✅ |
| C-86 | USAGE_GUIDE §3.10 同步 | 读 :693 | RNG 行与 §3.1.3 同表述 | :694 `Philox4x64-10（Random123，与 CPU 同一实现：uniform 位精确；正态 1e-15）` | ✅ |
| C-87 | ARCHITECTURE 假条目清除 | `Select-String -Path docs/architecture/ARCHITECTURE.md -Pattern 'Threefry'` | 0 命中 | 0 命中 | ✅ |
| C-88 | ARCHITECTURE 依赖列更新 | 读 :91、:68 | :91 写 `Philox4x64` + vendor 路径；:68 ADR-004 注 vendor | :68 ADR-004 注 `third_party/random123`，:91 依赖列 `Random123 (vendor: third_party/random123/include)` | ✅ |

---

## 7. 回归、性能与数值影响

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-91 | 全量回归（三目录） | `ctest --test-dir build -C Release -j 8 --output-on-failure`（其余两目录同） | 三目录 0 失败 | `tmp/verify_capi` **2600/2600**（2026-10-22 修复后终态, 干净目录口径 = 基线 2596 + 本审计新增 4 用例）；`build` **2592 过 / 4 NOT_BUILT**（E-12 族既有注册残留, 0 真失败）；`build_cuda` 不涉及本次修复改动（rng.hpp 仅注释、gpu_mc.cu 未动），维持 C-91 原记录 **2538/2538**。注册数代次调和见 IMPLEMENTATION §9 P2-8 注 | ✅ |
| C-92 | 失败处置合规 | 检查所有失败修复的 diff | 仅"重选 seed"或"修实现"；**无** `EXPECT_NEAR` 阈值放宽 | **3 处 seed 编辑覆盖 4 个统计窗口红用例**（NgPerron stream 2→3 一处覆盖 2 用例、Phase8 301→302、IS Put 42→43）+ 1 处参照口径修正（几何亚洲：连续 KV → 离散精确闭式解，容差 0.10→**0.05 收紧**；§5.8(c)）。全 diff 无阈值放宽 | ✅ |
| C-93 | 夹具零改动 | `git diff --stat tests/fixtures` | 仅新增 `rng/kat_philox4x64_10.hpp`，既有文件零改动 | `git status --short tests/fixtures` = `?? tests/fixtures/rng/`（仅新增） | ✅ |
| C-94 | 性能前后已记录 | `benchmark_philox` 改造前后各跑；读 `docs/DEVELOPMENT_LOG.md` 新小节 | 含时间戳、平台、numbers/sec、ns/number、旧/新对照与结论 | DEVELOPMENT_LOG :3032-3040：MSVC 19.50/Release，旧 1.007e8（9.93 ns），新 4.213e8（2.37 ns），≈4.18× | ✅ |
| C-95 | 性能未异常退化 | 对比两次数字 | 未出现 >2× 下降；如出现已有书面可接受性说明 | 新实现快 4.18×，无退化 | ✅ |
| C-96 | Python 侧语义变更已披露 | `Select-String -Recurse -Path python -Pattern 'next4'`；读 docstring | 若暴露，docstring 已写"4 字 = 1 块"语义 | `cpphub_python.cpp:319-321` 暴露 `next4`，docstring 写 `Draw exactly one block = four 64-bit unsigned integers` | ✅ |

---

## 8. 提交与回滚审查

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-101 | 提交 1（R-033）独立可 revert | `git show --numstat`（**不可用 --stat 图形刻度读增删行**, ADD 审计 P2-2） | 仅 `CMakeLists.txt` + `ci.yml` 两处 cache key | `7e305fc`：2 文件，numstat 实测 ci.yml +2/−2、CMakeLists.txt +6/−1 | ✅ |
| C-102 | 提交 2（R-032）文件集合正确 | 同上 | `VERSION`、`tools/sync_version.py`、`version.hpp`、`config.hpp`、根 CMake、python CMake、`cpphub_python.cpp`、`pyproject.toml`、`test_version.py`、`ci.yml`、README（版本两轴模型段现 :9-13 / 跨平台确定性行现 :29 / GPU skip 行现 :113 / 顶部——内容锚, ADD 审计 P2-1） | `5c8e21b`：13 文件 (= 上述 11 + CHECKLIST/IMPLEMENTATION 文档) | ✅ |
| C-103 | 提交 3（R-034）文件集合正确 | 同上 | `json_min.hpp`、`registry.hpp`、`c_api.cpp`、`test_json_min.cpp`、`test_c_api_tools.cpp` | `253fd94`：6 文件 (= 上述 5 + CHECKLIST) | ✅ |
| C-104 | 提交 4（R-027）整体可 revert | 同上 | `third_party/random123/**`、`rng.hpp`、`gpu_mc.cu`、`rng_device_probe.hpp`、三处 CMake include + 两处测试 CMake、6 个测试文件（既有 2 文件加用例 test_json_min/test_c_api_tools + 4 个回归处置, ADD 审计 P2-2 勘误"新增 2 个"）、README GPU 行、USAGE_GUIDE、ARCHITECTURE、夹具头 | `5910753`：49 文件 = 上述全部 + vendor 27 头/LICENSE/UPSTREAM + 文档与 CHECKLIST | ✅ |
| C-105 | README 无重复修改 | `git log -p -- README.md` | :46 只在提交 4，:21/:107/:121 只在提交 2 | `5910753` README 仅 GPU MC 行；`5c8e21b` 含版本模型/跨平台行/GPU skip 行 | ✅ |
| C-106 | 每提交前已 review | 提交说明或 review 记录 | 4 次 review 均有结论（排除幻觉） | 4 次提交前均逐文件读 diff；提交 4 review 修正 `tests/CMakeLists.txt` 注释笔误（grunfeld）、`kv_price`→`ref_price`，并核对全部行号引用 | ✅ |
| C-107 | 回滚演练（可选） | 在临时分支 `git revert --no-commit <R-027>` | 回滚后 `build/` 核心配置可编译（GPU 目录除外） | 未执行（清单标注为可选，不阻塞验收） | ⏸ |

---

## 9. 验收标准映射（DESIGN §9 A1-A10）

| 标准 | 内容 | 本清单条目 | 状态 |
|---|---|---|---|
| A1 | 无 tag 依赖，FetchContent 钉 commit | C-11 / C-12 / C-13 | ✅ |
| A2 | 两轴版本模型三者一致 | C-21 ~ C-29 | ✅ |
| A3 | README 事实修正、无冻结总数 | C-31 / C-32 / C-33 / C-34 / C-35 | ✅ |
| A4 | 覆盖率 job 产出 artifact | C-38 | ✅ |
| A5 | C API 测试在 CI 实际运行 | C-36 / C-40 | ✅（C-36 ✅；C-40 ✅ run 39755779788 双平台 2600/2600） |
| A6 | 8 MiB+1 拒绝；节点预算边界正确 | C-41 ~ C-52 | ✅ |
| A7 | KAT（CPU/GPU）全绿；uniform 逐字相同 | C-77 ~ C-81 | ✅ |
| A8 | 正态层级 ≤1e-15 | C-82 | ✅ |
| A9 | 全量回归通过（允许重选 seed） | C-91 / C-92 / C-93 | ✅ |
| A10 | 性能前后数据已记录 | C-94 / C-95 | ✅ |

**覆盖声明**：A1-A10 全部有映射；本清单另有 G/C 系列超出 A1-A10 的工程项（构建、探针安全、提交可回滚性），均来自 IMPLEMENTATION §1-§7，不引入 DESIGN 外的新裁决。

---

## 10. 签收

| 项 | 填写 |
|---|---|
| 实施完成日期 | 2026-10-22 |
| 三目录通过数（R-027 实施时点；ADD 修复后终态见下行） | build: **2596/2596** / build_capi（`tmp/verify_capi`）: **2596/2596** / build_cuda: **2538/2538**（`gpu_mc` 内含 22/22） |
| 性能（旧 → 新，numbers/sec） | **1.007e8 → 4.213e8**（9.93 → 2.37 ns/number，≈4.18×） |
| 全量回归失败处置摘要 | 首轮 5 红（**计数单位 = 红用例数**；编辑数 3+1）：4 例统计窗口红（3 处 seed 编辑覆盖: NgPerron stream 2→3 一处覆盖 2 用例 / Phase8 301→302 / IS Put 42→43）；1 例参照口径修正（几何亚洲改为离散监控精确闭式解，容差 0.10→0.05 收紧）。**无阈值放宽**。修复后终态: `tmp/verify_capi` 2600/2600（ADD 审计后, 基线 2596 + 新增 4 用例） |
| 结论 | ✅ 通过（Step 10 ADD 审计 Tier 2 修复全清, 2026-10-22） |
| 实施人 / 复核人 | Scott (鹏) / ADD 审计（主 agent 复核 + 双 subagent 交叉） |

---

## 11. CI 修复批（run 39755779788 五簇 → gcovr）

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-108 | 最终 CI run 全绿 | `gh run view 40068527562` | 9/9 job success | ✅ run 40068527562 **9/9 job success**（Coverage 26m15s、MSVC 19m04s、C ABI win 15m36s）；随后文档提交 run 40069392810 亦 success | ✅ |
| C-109 | guard 步骤不再自匹配 | 看 Build & Test 两平台 Guard 步骤 | success | ✅ run 39784068194 起 Guard 步骤 success（字符类 `gtest-1[.]14[.]0`） | ✅ |
| C-110 | Show toolchain 不依赖 nanobind | 看 4 个 Python job 的 Show toolchain versions 步骤 | success | ✅ run 39784068194 四 job 该步骤 success（仅 import sys/cpphub + mcp） | ✅ |
| C-111 | coverage Debug 全绿 | 看 coverage job Test 步骤 | `100% tests passed out of 2538` | ✅ run 39990450582 实测该行（18 项断言失败全清：14 ArimaModel + Midas + VarIrfFevd + ArimaStatespace） | ✅ |
| C-112 | py3.11 structured_content 修复 | 看 py3.11 job pytest 汇总行 | 全绿且非 24 failed | ✅ run 39990450582 py3.11（windows/ubuntu）`111 passed, 1 skipped`；本地同版本 venv 复现修复前 24 failed → 修复后 51 passed 1 skipped | ✅ |
| C-113 | 下溢/维度修复本地红→绿 | MSVC Debug 目标套件 | 全通过 | ✅ 修复前 14 ArimaModel + Midas + VarIrfFevd abort；修复后 test_arima_model 27/27、test_midas_model 19/19、test_linalg_dynamic 18/18、VarIrfFevd 18/18、VarModel+DySpillover 34/34、波及 8 套件 93 tests | ✅ |
| C-114 | gcovr 步骤不再中止 | 看 coverage job gcovr report 步骤 | success | ✅ run 40068527562 Coverage job 全步骤 success（含 gcovr report）；修复 = `--gcov-ignore-parse-errors suspicious_hits.warn`（gcov bug 68080） | ✅ |

---

## 12. P2 风险修复批（RISK-029/030/031, 2026-10-24）

> 决策: [ADR-021](../../docs/decisions/ADR-021_PDE_CONCURRENCY_BARRIER_BOUNDARY.md)（Accepted, 6 项）；日志: DEVELOPMENT_LOG「P2 风险修复批」条目；提交 `b2339fd`（029）/ `0638367`（ADR）/ `62c64ed`（CMake 回归）/ `32b7894`（030+031）

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-115 | RISK-029 opt-in 守卫 + CI | `gh run view 40093818921`；Ninja Release compile_commands.json 对照 | 默认 OFF 基线；ON 时旗标落点正确 | ✅ 双变体实测：OFF=0 march 命中 + /O2，ON=166 处 /arch:AVX2；run 40093818921 **9/9 job success**（c-api ubuntu 实测 opt-in 路径） | ✅ |
| C-116 | RISK-030 共享态消除 | `grep last_total_iter` 源码清零；test_pde_engine 位相同断言 | mutable 成员/读点/写点删除，3 处测试迁移到 `price_american_result()` | ✅ 源码仅剩注释；`.value` 与旧包装 EXPECT_EQ 位相同、`psor_iterations>0`；全量 39/39（MSVC Debug） | ✅ |
| C-117 | RISK-031 kind() 分派 + 测试锚 | 数字障碍默认域 vs 宽域（s_multiplier=20）；Neumann≠Dirichlet；Unknown throw；override discount | 10 新锚全绿；数字障碍域一致性 | ✅ 先红（编译红）后绿；域对比 NEAR 0.02、Neumann 差异 >1e-8 + 双侧解析解 1.5%、Risk031SqPayOff 冒名 name()="Call" 仍 throw | ✅ |
| C-118 | CMake D8016 回归修复 | VS 双配置 vcxproj 语义属性核验 | Debug 无 /O2 且保留 /RTC1；Release MaxSpeed | ✅ 62c64ed：Debug `Optimization=Disabled` + `BasicRuntimeChecks=EnableFastChecks`，Release `MaxSpeed`；优化级别交 per-config 默认 flags | ✅ |
| C-119 | 最终 CI run 全绿 | `gh run view 40107640608` | 9/9 job success | ✅ run 40107640608（head 32b7894）**9/9 job success**：Coverage 31m16s / Build&Test GCC 11m28s / Build&Test MSVC 14m03s / C ABI ×2 / Python ×4；仅 Node.js 20 弃用告警（非失败） | ✅ |

---

## 13. P2 数值契约批（RISK-035/036, 2026-10-25）

> 决策: [RESEARCH_P2_RISKS.md](./RESEARCH_P2_RISKS.md) §4/§5（方案①×2）；日志: DEVELOPMENT_LOG「P2 数值契约批」条目；提交 `9c3a435`（035+036）

| # | 验收项 | 执行方法 | 通过条件 | 证据 | 状态 |
|---|---|---|---|---|---|
| C-120 | RISK-035 statistic 语义 | registry :544 改 `r.chow_denning_statistic`；golden 第 5 数反解交叉验证 | statistic = max&#124;Z2&#124; 非负、≠ p 值；golden 容差过 | ✅ 反解交叉验证 `1-(2Φ(1.3385846942648503)-1)^4 = 0.549433098981491` 与冻结 p 逐位一致 → cd_stat = 1.3385846942648503（k=2 &#124;Z2&#124; 即 max）；`EXPECT_NE(statistic, p_value)` 语义防线 | ✅ |
| C-121 | RISK-036 上尾精度 | math.hpp 重写 upper + `chi2_sf` 新增 + backtesting 三处迁移 | k=1/LR=74 锚 7.81e-18；k=5 三锚相对 1e-12；P 逐位不变 | ✅ 先红（C3861/C2039 编译红 + Kupiec/Christoffersen 2 FAILED p=0）后绿；test_core_math 4 新锚 / test_backtesting 23/23 / test_variance_ratio_test 24/24 / test_cev 21/21；midas_diagnostics:285 按 §5.6 留后续 | ✅ |
| C-122 | Python golden 层 | 重建 nanobind 绑定 + pytest 全量 | 结构 byte-exact + 数值容差 | ✅ 111 passed 1 skipped（py3.12 本地，与 CI 同口径）；F9 known_issues 清空后结构一致 | ✅ |
| C-123 | 最终 CI run 全绿 | `gh run view 40212222516` | 9/9 job success | ✅ run 40212222516（head 9c3a435）**9/9 job success**：Coverage 31m11s / Build&Test MSVC 18m02s / Build&Test GCC 11m43s / C ABI ×2（21m53s+12m55s）/ Python ×4；仅 Node.js 20 弃用告警（非失败） | ✅ |

---

**编制依据**: [DESIGN.md](./DESIGN.md) §9 + [IMPLEMENTATION.md](./IMPLEMENTATION.md) §0-§9 · **编制日期**: 2026-10-21 · **执行状态**: **已完成**（Step 7-9 全部执行并签收; Step 10 ADD 审计 21 项发现 Tier 2 全修, 跟踪表见 [AUDIT.md](./AUDIT.md) §0; §11 CI 修复批 2026-10-22→24; §12 P2 风险修复批 2026-10-24; §13 P2 数值契约批 2026-10-25）
