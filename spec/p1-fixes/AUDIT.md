# ADD 审计报告：P1 整改（p1-fixes）代码-设计-文档三方对齐

---

id: P1-FIXES-AUDIT
type: audit
version: 1.0
status: closed
date: 2026-10-22
depends: [P1-FIXES-CHECKLIST, P1-FIXES-IMPLEMENTATION, P1-FIXES-DESIGN, P1-FIXES-RESEARCH]
baseline: 5741143

---

> **Feature**: p1-fixes（RISK-027 / 032 / 033 / 034）
> **Spec 步骤**: Step 10（10 步流程末环）
> **审计框架**: audit-driven-development v3.5+（Detection → Repair；本文为 Detection 产出，Repair 人机协作、永不自动修复）
> **评分口径**: `score = max(0, 100 − (P0×20 + P1×8 + P2×3 + P3×1))`

---

## 0. 修复跟踪表

| 编号 | 严重度 | 模块 | 描述 | 状态 | 提交 | 测试验证 |
|---|---|---|---|---|---|---|
| P2-1 | P2 | 跨模块 | spec/docs 行号锚点漂移集群（C-12/C-14/C-38/C-75/C-102/README 引用/IMPLEMENTATION §2.1） | ✅ fixed | 5ab43d2 | 内容锚化 + PowerShell 复测行号（:72/:59/:115/:136/:503 等） |
| P2-2 | P2 | 跨模块 | CHECKLIST §8 提交证据数字失实（C-101 stat、C-104「新增 2 个测试文件」） | ✅ fixed | 5ab43d2 | `git show --numstat` 实测 +2/−2、+6/−1 回写 |
| P2-3 | P2 | R-032 | python/CMakeLists.txt standalone 回退分支 `__version__` 落 C++ 轴 1.0.0，与 test_version 断言 pyproject 1.9.0 互斥 | ✅ fixed | 5ab43d2 | standalone configure exit 0 + `CMAKE_PROJECT_VERSION=1.9.0`；pytest test_version PASSED |
| P2-4 | P2 | R-034 | 「对象键不计入节点预算」语义承诺测试锚定未兑现（DESIGN §5.3 / IMPLEMENTATION §4.2） | ✅ fixed | 5ab43d2 | `JsonMinBudget.ObjectKeysNotCountedInNodeBudget` PASSED |
| P2-5 | P2 | R-027 | ADR_INDEX.md:16/:180 仍宣称 Philox4x64 / Threefry4x64，与 ARCHITECTURE 已删假条目口径矛盾 | ✅ fixed | 5ab43d2 | 两处勘误回写，diff review 确认 |
| P2-6 | P2 | 跨模块 | CHECKLIST 头部 `status: pending` / :16「§2-§8 全部 ☐」/ :267「待 Step 9」与正文全 ✅ 矛盾；§10 结论行与签名人空置 | ✅ fixed | 5ab43d2 | 头部 audited + §10 结论 ✅ 通过 + 签收回填 |
| P2-7 | P2 | 跨模块 | 回归处置计数单位不统一：DEVELOPMENT_LOG「3 例 seed」（编辑数）vs CHECKLIST C-92「4 例重选种子」（单位错标）vs §10「4 例统计窗口」（红用例数）；事实 = 3 处 seed 编辑覆盖 4 个红用例 + 1 例参照修正 | ✅ fixed | 5ab43d2 | 三处统一为「3 处编辑覆盖 4 红用例 + 1 参照修正（5 红）」口径 |
| P2-8 | P2 | 跨模块 | 三目录历史测试计数（E-9=2553 / G-05=2605 / E-12=2573,2590 / C-91=2596；build_cuda 2548→2538）未补一行调和说明 | ✅ fixed | 5ab43d2 | IMPLEMENTATION §9 末补代次调和注（登记事实, 不作无证据归因） |
| P3-1 | P3 | R-033 | cache key 与 CMake URL SHA 同步无 CI 防回退断言 | ✅ fixed | 5ab43d2 | guard 步骤入 ci.yml；PyYAML 校验通过；变量拼接防 grep 自匹配 |
| P3-2 | P3 | R-032 | version.hpp 生成头 2 行两轴注释超出 DESIGN §4.2 B 模板（仅 1 行） | ✅ fixed | 5ab43d2 | DESIGN §4.2 补两轴注 |
| P3-3 | P3 | R-032 | README:7 保留 v1.4.0–v1.4.3 冻结叙事，两轴模型下易误读为当前版本 | ✅ fixed | 5ab43d2 | 改无版本号功能描述 + 指向版本历史表 |
| P3-4 | P3 | R-034 | 节点预算缺最小拒绝量用例（恰 +1 节点） | ✅ fixed | 5ab43d2 | `JsonMinBudget.RejectsFirstNodeOverBudget` PASSED |
| P3-5 | P3 | R-034 | args_json 缺「恰等 8 MiB 应接受」用例，`>` 语义无回归锁 | ✅ fixed | 5ab43d2 | `CApiTools.AcceptsArgsJsonAtExactLimit` PASSED |
| P3-6 | P3 | R-027 | rng.hpp:15 注释「对外语义不变」失准（API 不变但输出流已变） | ✅ fixed | 5ab43d2 | 注释指向 USAGE_GUIDE §3.1.3；2600/2600 无回归 |
| P3-7 | P3 | R-027 | PHASE1_SPEC.md:78 仍写 Random123「git submodule \| 1.15」 | ✅ fixed | 5ab43d2 | 单元格补 in-tree vendor 弃用注 |
| P3-8 | P3 | R-027 | IMPLEMENTATION §5.8 处置枚举缺第 (c) 类「测试参照口径错误 → 修参照」 | ✅ fixed | 5ab43d2 | §5.8 补 (c) 类 + 实例（非放宽容差） |
| P3-9 | P3 | R-027 | `discard(n)` 整块倍数路径（rem==0 延迟 refill）无用例锚定 | ✅ fixed | 5ab43d2 | `Rng.PhiloxDiscardExactBlocks` PASSED |
| P3-10 | P3 | 跨模块 | DESIGN §6.5 仍写「examples/kat_vectors」，E-7 勘误（tests/）未回写 | ✅ fixed | 5ab43d2 | §6.5 勘误标注回写 |
| P3-11 | P3 | 跨模块 | DESIGN/IMPLEMENTATION 头部状态行停更（写「可进入 Step 9/Step 7」，正文已含 Step 9 执行证据） | ✅ fixed | 5ab43d2 | 两文档头部更新为 Step 9 完成 + Step 10 回写 |
| P3-12 | P3 | 跨模块 | A5/C-40 ⏸（待 CI 日志）无跟踪条款约定何时补证 | ✅ fixed | 5ab43d2 | C-40 补跟踪条款（推送后取证回填） |
| OBS-1 | 观察 | R-033 | URL_HASH 缺位属留痕的有意取舍（DESIGN §3.4 裁决），不计分；建议留存一次 codeload@commit 归档 SHA256 备查 | ✅ done | 5ab43d2 | 归档 SHA-256 留存于 CHECKLIST C-12 证据列 |

> **修复后终态（2026-10-22）**：`tmp/verify_capi`（干净目录口径）**2600/2600** 全绿（= 修复前基线 2596 + 本审计新增 4 用例）；`pytest test_version.py` PASSED；`ci.yml` PyYAML 解析通过。21 项全部闭合，预期评分 87 → **A**（P2 全清）；文档项与代码项均已回写（P3 全清, 观察项 done）。

---

## 1. 审计范围与方法

### 1.1 Phase 0 质量门

spec 四件套按 5 维评分（可测试约束 / 模块映射 / 接口契约 / 修正项 / 跨模块契约）= **5/5，Tier A** → 全量审计（Phase 1-4），无评分上限限制。

### 1.2 审计维度与执行方式

基线 = `5741143`（工作树干净）。5 个审计 subagent 并行（4 个风险项模块审计 + 1 个跨模块契约审计），全部只读；主 agent 对 7 项关键 P2 做源码级复核（Phase 3.5 幻觉防护），全部属实。

| 维度 | 覆盖对象 | 模块评分 |
|---|---|---|
| R-033 供应链 | CMakeLists.txt:64-76、ci.yml 缓存段、7e305fc diff | A- |
| R-032 版本两轴 | VERSION/version.hpp/sync_version.py/python 构建链/README/ci.yml、5c8e21b diff | A- |
| R-034 JSON 预算 | json_min.hpp/registry/c_api 入口/新测试、253fd94 diff | A |
| R-027 Random123 | vendor/rng.hpp/KAT 夹具/gpu_mc.cu/回归处置 4 文件/docs、5910753 diff | A |
| 跨模块契约 | A1-A10 映射链、提交边界、CHECKLIST 证据抽检 11 项、spec 元数据 | B+ |

### 1.3 汇总评分

**P0 = 0，P1 = 0，P2 = 8，P3 = 12**（subagent 原始 25 条发现去重合并，映射见附录 A）。

```
score = 100 − (0×20 + 0×8 + 8×3 + 12×1) = 64 → C+（55-64）
```

**解读**：数值分被跨模块文档证据漂移单一根因拖累（P2-1/P2-2/P2-6/P2-7/P2-8 同源：后续提交移动行号 + 计数单位混用）；四个模块级评分 A/A-/A-/A 表明**实现与设计本体高度对齐，0 P0/0 P1**——核心不变式（供应链不可变身份、版本单一源、预算硬上限、vendor 哈希全等、KAT 逐字一致、无容差放宽）全部独立复验成立。

---

## 2. P2 发现详情（均已主 agent 复核，除 P2-8 采信 subagent 中置信度）

### P2-3（最实质，代码相邻）standalone 回退分支两轴混轴
- **evidence**: python/CMakeLists.txt:15-19（回退 `../VERSION` → CPPHUB_PY_VERSION=1.0.0）· python/tests/test_version.py:17（断言 `__version__ == pyproject version`，即 1.9.0）
- **spec_ref**: DESIGN §4.2 E + 初始值节
- **claim**: 非 wheel 路径（直接 `cmake -S python`）下 `__version__` = 1.0.0，与 test_version 断言互斥；回退分支仅验到 configure exit 0，未测 pytest。
- **fix_cost**: 低 · **regression_risk**: low · **fix_suggestion**: 二选一——回退分支改为解析 pyproject version；或 test_version 对 standalone 构建跳过并注释（wheel 主路径不受影响）。

### P2-4「键不计入」语义承诺未兑现测试
- **evidence**: include/cpphub/c_api/detail/json_min.hpp:29-34（注释承诺）与 :257-258（键直调 parse_string，行为正确）· tests/unit/c_api/test_json_min.cpp:350-389（无键语义用例）
- **spec_ref**: DESIGN §5.3 / IMPLEMENTATION §4.2（「计数语义……被测试锚定」）
- **fix_cost**: 低 · **fix_suggestion**: 补一个「大量键少量值」用例断言成功解析，反向锚定键不计入。

### P2-5 ADR_INDEX 与 ARCHITECTURE 口径矛盾
- **evidence**: docs/decisions/ADR_INDEX.md:16、:180（「Philox/Threefry」「Philox4x64 / Threefry4x64」）vs docs/architecture/ARCHITECTURE.md:68/:91（已删 Threefry 假条目）；DEVELOPMENT_LOG:3041 自述「ARCHITECTURE 删 Threefry4x64 假条目、ADR-004 注明 vendor」但 ADR_INDEX 未同步。
- **spec_ref**: DESIGN §6.6（修正表只列 ARCHITECTURE 两处，属 spec 覆盖缺口）
- **fix_cost**: 低（改 2 行）· **fix_suggestion**: ADR_INDEX.md:16/:180 补「Threefry 未实现（2026-10-22 R-027 勘误）」备注。

### P2-6 CHECKLIST 头部状态与正文矛盾 + §10 签收空置
- **evidence**: CHECKLIST.md:8（status: pending）、:16（「§2-§8 全部 ☐」）、:267（「执行状态: 待 Step 9」）vs :70-248 全部 ✅；:262 结论行 ☐ 未勾、:263 实施人/复核人空。
- **fix_cost**: 低 · **fix_suggestion**: 头部改「Step 9 已完成（2026-10-22），Step 10 审计中」，§10 勾选结论（按 P2 修复决策定「通过」或「有条件通过」），填签名人。

### P2-7 回归处置计数单位不统一
- **evidence**: docs/DEVELOPMENT_LOG.md:3043-3044（「重选 seed 3 例」= 编辑数 +「参照口径修正 1 例」）· CHECKLIST.md:213 C-92（「4 例重选种子……覆盖 4 用例」——单位错标，实为 3 处编辑覆盖 4 用例）· CHECKLIST.md:261 §10（「4 例统计窗口」= 红用例数，单位未声明）
- **spec_ref**: IMPLEMENTATION §5.8
- **fix_cost**: 低 · **fix_suggestion**: 三处统一为「3 处 seed 变更覆盖 4 个统计窗口红用例 + 1 处参照口径修正（1 个红用例），合计 5 红 / 4 文件」。

### P2-8 三目录历史计数未调和（confidence: medium，subagent 复核）
- **evidence**: IMPLEMENTATION.md:1045/1048（E-9=2553、G-05 引 E-9 但写 2605、E-12=2573/2590）· CHECKLIST.md:74/212（G-05=2605、C-91=2596）· build_cuda 2548→2538（G-05 vs C-91）无解释
- **fix_cost**: 低 · **fix_suggestion**: 在 E-9/E-12 或 C-91 补一行：各数字对应的 configure 代次与残留清理事件（E-9 残留源 / E-12 注册表污染 / 重配清 NOT_BUILT）。

### P2-1 行号锚点漂移集群
- **evidence**: CMakeLists.txt:72（URL 行；C-12 锚 :63-71）· ci.yml:115（c-api cache key；C-14 锚 :110）· ci.yml:136（coverage job；C-38 锚 :132）· README 实际 :27/:113/:127（C-102 锚 :21/:107/:121）· gpu_mc.cu 守卫实际 :503/:515/:527（C-75 锚 −1）· IMPLEMENTATION §2.1 锚 :61-68（实际块 :64-76）
- **成因**: 7e305fc 之后 5c8e21b 等提交头部增行（CMakeLists +2、ci.yml build-test +41）。
- **fix_cost**: 低 · **fix_suggestion**: 集中刷新一遍锚点，或改用内容锚（URL 字面量 / `gtest-f8d7d77c`）。

### P2-2 CHECKLIST §8 提交证据数字失实
- **evidence**: git show --numstat 7e305fc = ci.yml +2/−2、CMakeLists.txt +6/−1（C-101 写 +4/−2、+7/−1，疑将 --stat 图形刻度误读为计数）· C-104「6 个测试文件（新增 2 个…）」失实——test_rng.cpp/test_gpu_mc.cpp 均为既有文件加用例，无新增测试 .cpp
- **fix_cost**: 低 · **fix_suggestion**: 改为实测数字；「新增 2 个」改「既有 2 文件加用例」。

---

## 3. P3 发现（紧凑表）

| # | evidence | spec_ref | claim | fix_cost |
|---|---|---|---|---|
| P3-1 | ci.yml 全文无 URL SHA / cache key 断言步骤 | DESIGN §3.2 | 无防回退门禁，同步靠人工 | <10-line |
| P3-2 | version.hpp:1-3 vs DESIGN §4.2 B 模板 | DESIGN §4.2 B | 生成头多 2 行注释（脚本以 TEMPLATE 为准，行为无差） | 1-line |
| P3-3 | README.md:7 | RESEARCH §3.2.1:211 | v1.4.0–v1.4.3 冻结叙事遗留 | 低 |
| P3-4 | test_json_min.cpp:365-377 | DESIGN §5.5 | 缺恰 +1 节点最小拒绝量用例（spec 本身只要求 +2，属 spec 级盲区） | 低 |
| P3-5 | test_c_api_tools.cpp:334-355 | DESIGN §5.2 | 缺恰等 8 MiB 接受侧用例，`>` 语义无锁 | 低 |
| P3-6 | rng.hpp:15 | DESIGN §11 D-2、§6.8 | 「对外语义不变」应改「对外 API 不变，输出流语义变更见 USAGE_GUIDE §3.1.3」 | 1-line |
| P3-7 | PHASE1_SPEC.md:78 vs UPSTREAM.md:3-8 | DESIGN §6.1 | submodule/1.15 遗留，与 in-tree vendor 矛盾 | 低 |
| P3-8 | test_multi_asset_mc.cpp:147-166 | IMPLEMENTATION §5.8 | 处置枚举缺第 (c) 类「参照口径错误→修参照」（本例方向合规且已透明记录） | 低 |
| P3-9 | rng.hpp:57-62 vs test_rng.cpp:94-102 | DESIGN §6.3 | discard rem==0 延迟 refill 路径无用例 | 低 |
| P3-10 | DESIGN.md:370 vs E-7 勘误 | IMPLEMENTATION §9 | examples/kat_vectors 未回写 | 1-line |
| P3-11 | IMPLEMENTATION.md:5,91 · DESIGN.md:5,477 | 两文档头部 | 状态行停更于 Step 9 前 | 1-line |
| P3-12 | CHECKLIST.md:128,243 | §4.3/§9 | A5/C-40 ⏸ 无补证跟踪条款 | 1-line |

**观察项 OBS-1**（不计分）：URL_HASH 缺位为 DESIGN §3.4 留痕裁决；可选加固 = 留存一次 codeload@commit 归档 SHA256 备查。

---

## 4. 修复优先级矩阵

- **Tier 1（P0，立即）**：空。
- **Tier 2（本轮）**：P2-3（唯一代码相邻项，1 处分支/测试决策）→ P2-4/P3-4/P3-5（测试补强，与 P2-3 可并入一次 TDD 微循环）→ P2-5/P2-6/P2-7/P2-8/P2-1/P2-2（文档一致性批量回写）。
- **Tier 3（可选）**：P3-1~P3-12 + OBS-1。

排序依据：代码语义 > 测试锚定 > 文档证据；1 行修复优先。

## 5. 测试盲区汇总

1. R-034：节点预算恰 +1 拒绝、8 MiB 恰等接受、键不计入语义（P2-4/P3-4/P3-5）。
2. R-032：standalone 回退路径、`CPPHUB_PY_VERSION="unknown"` 兜底分支（cpphub_python.cpp:2043-2045）零覆盖（P2-3）。
3. R-027：discard 整块倍数路径（P3-9）；KAT 仅覆盖 0/全 f/pi 三类边界，中间 counter/key 域无官方锚（上游仅有 3 行，覆盖度声明诚实）。
4. R-033：cache key↔URL SHA 同步无 CI 断言（P3-1）。

## 6. CHECKLIST 静态证据抽检表（跨模块 subagent，11 项）

| 抽检项 | 声明 | 重验结果 | 一致？ |
|---|---|---|---|
| C-11 | archive/refs/tags 0 命中 | 0 命中 | ✅ |
| C-14 | 旧 key 0 / 新 key 2 | 0 / 2（:59,:115；锚行漂移→P2-1） | ✅ |
| C-29 | 1.4.0 清零 | src/ 与 python/src 0 命中 | ✅ |
| C-31 | README 1412 清零 | 0 命中（含「15 GPU MC」） | ✅ |
| C-21 | VERSION 单行 LF | `1.0.0`、lf（逐字节 hex 未验） | ✅ |
| S8-3 | test_rng 既有 9 用例 | 253fd94 处 =9；现树 10 = 9+RngKat | ✅ |
| G-02 | UPSTREAM.md 关键字段 | URL+40 位 SHA+日期+范围+许可证齐全 | ✅ |
| C-61 | vendor 27 文件 | git ls-files = 29 = 27 头 + LICENSE + UPSTREAM | ✅ |
| §10 vs 正文 | 2596/2596/2538、1.007e8→4.213e8 | 与 C-91/C-94 一致 | ✅ |
| D-02 | A 表 30 行 | 10+10+10 实测 | ✅ |
| C-92 vs 日志 | seed 计数 | 3 编辑 vs「4 例」——单位错标 | ❌（P2-7） |

**独立复验成立的核心不变式**（主 agent / subagent 双重验证）：自写 Philox 残留 0；vendor 27 头 SHA-256 全等 + LICENSE 逐字节相同；KAT 3 向量与上游逐字一致；D-6 宏守卫位于首次 include 之前且前置头无传递包含；回归处置零容差放宽；预算常量/入口检查/错误文案与 spec 逐字一致；README 分段提交边界（C-105）零串门。

## 7. 修复后预期评分

- P2 全清（8 项）：64 → **87（A）**
- P2 + P3 全清（20 项 + OBS 落实）：→ **100（A+）**

## 8. 下一步建议

- A.（推荐）Tier 2 一次闭环：P2-3 走 TDD 微循环（红→绿），P2-4/P3-4/P3-5 三个边界用例并入同轮；其余 P2 文档批量回写后统一提交。
- B. Tier 3 顺手项并入 A 的提交；P3-1（CI 断言）可推迟到 C-40 CI 日志取证同一推送。
- C. A5/C-40：推送后以 GitHub Actions 日志补证，与 P3-12 跟踪条款一并闭环。

## 附录 A：发现去重映射

| 报告 ID | subagent 原始 ID |
|---|---|
| P2-1 | R033-1 + XC-06 + R032-4 + R027-5 |
| P2-2 | R033-2 ≡ XC-04 + XC-05 |
| P2-3 | R032-1 |
| P2-4 | R034-3 |
| P2-5 | R027-1 |
| P2-6 | XC-01（含 §10 结论/签收空置） |
| P2-7 | XC-02 |
| P2-8 | XC-03 |
| P3-1~P3-12 | R033-4；R032-2、R032-3；R034-1、R034-2；R027-2、R027-3、R027-4、R027-6；XC-07、XC-08、XC-09 |
| OBS-1 | R033-3 |

## 附录 B：验证局限（审计诚实声明）

- runtime 证据（三目录 ctest 计数、pytest 111 passed、nvcc 输出、性能 4.18×）未重跑，仅核对文档内部与代码一致性。
- P2-8 采信 subagent 中置信度；VERSION 逐字节 hex、version.hpp 字节数未独立验证（git show / ls-files --eol 替代覆盖）。
- 上游 kat_vectors 原始文件比对依赖工作树 tmp/r123_dl 解包树（vendor 27 头与之一致已独立复验）；codeload 可达性与 tag→commit 映射未在线复核。
- 工具降级记录：subagent 阶段 Glob/Read 各有超时，均以 Grep/git 只读命令替代完成，证据行号均为实读。
