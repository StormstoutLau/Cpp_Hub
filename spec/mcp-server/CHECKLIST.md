# 审查验收 Checklist：Cpp_Hub MCP Server（mcp-server）

---
id: MCP-SERVER-CHECKLIST
type: design
version: 1.0
status: pending
date: 2026-10-18
depends: [MCP-SERVER-IMPLEMENTATION, MCP-SERVER-DESIGN]
upstream: null
---

> **Feature**: mcp-server（Cpp_Hub 对外 MCP 接口层）
> **状态**: 待验收（全部验收项 ☐；实现尚未开始）
> **Spec 步骤**: Step 7–8, 10
> **基于实施**: [IMPLEMENTATION.md](./IMPLEMENTATION.md) · **基于设计**: [DESIGN.md](./DESIGN.md)
> **计数基准**：本节各表逐项登记；验收项总数见 §10.1（口径 = 逐项计数，非事后推算，RULE-2）

---

## 1. 文档一致性验收（Step 8）

> **执行记录（2026-10-18，独立 pass）**：按 RULE-4 标注 `` `[自查·单视角]` ``（同一 agent 体系执行，非同源复审）；机械核对命令见 §1.5。Step 8 发现 **1 × P3（F4）**；Step 9 实施期**追加 1 项核对维度**并发现 **1 × P2（F5）**——两项均已闭环，见 §1.5。

### 1.1 RESEARCH.md ↔ DESIGN.md 对齐

| 检查项 | 状态 | 说明 |
|---|---|---|
| DESIGN §2.1 七行决策均可追溯到 RESEARCH 锚点 | ✅ `[自查·单视角]` | 七行锚点（§3.1.1/§3.1.2/§3.6/K3/A19/A21 + MATRIX E-*）grep 全命中 |
| RESEARCH 关键发现（K1–K4 约束）被 DESIGN 采纳 | ✅ `[自查·单视角]` | K1→§8.1 无状态；K2→§4.1 命名；K3→§7 超长；K4→§10.1/§10.2 安全 |
| 无文档间矛盾（含已失效声明未被引用） | ✅ `[自查·单视角]` | H-2「25k tokens」在 DESIGN/IMPL 零出现（已复核） |

### 1.2 DESIGN.md ↔ IMPLEMENTATION.md 对齐

| 检查项 | 状态 | 说明 |
|---|---|---|
| DESIGN §3.2 六模块在 IMPL §3 均有对应实施 | ✅ `[自查·单视角]` | 6/6 同名同职责（IMPL 顺序不同，集合相同） |
| IMPL 接口签名与 DESIGN §4 一致 | ✅ `[自查·单视角]` | 内容无缺失；**小节编号非 1:1（F4）已按建议 (a) 处置**——IMPL §4 首部增「↔ DESIGN 映射表」（Step 8 回写），DESIGN 未改动 |
| DESIGN §4.1「出参要点」粒度 ↔ §6.1/§6.2 出参结构粒度一致 | ✅ `[自查·单视角·Step 9 追加]` | **本项为 Step 9 实施期补入的核对维度**（原 §1.2 无此维度，故 Step 8 漏检 F5）；7/7 工具为多统计量 → 由 §6.2 `ToolResult.results[]` 承接，见 §1.5 与 [DECISION_MATRIX.md](./DECISION_MATRIX.md) §9–§12 |
| DESIGN §8 六条不变式在 IMPL §7 均有实施位置 | ✅ `[自查·单视角]` | 6/6 逐条对应 |
| 无设计未覆盖的实施 | ✅ `[自查·单视角]` | IMPL 未新增 DESIGN 外模块；B2 明确列终态不实施 |

### 1.3 IMPLEMENTATION.md ↔ CHECKLIST.md 对齐

| 检查项 | 状态 | 说明 |
|---|---|---|
| IMPL §8 测试策略的功能点在本 checklist 有验收项 | ✅ `[自查·单视角]` | IMPL §8.1 四模块 ↔ CHECKLIST §2.1–2.5 |
| 本 checklist 验收项均可追溯到 IMPL 章节 | ✅ `[自查·单视角]` | 每节标注（IMPL §3.x / §4 / §6 / §7） |

### 1.4 术语一致性

| 术语 | RESEARCH.md | DESIGN.md | IMPLEMENTATION.md | 一致 |
|---|---|---|---|---|
| 诊断信封 | DiagnosticEnvelope 提及 | §6.1 `DiagnosticEnvelope` | §4.1 同名 dataclass | ✅ |
| 工具面 | ≤10 / 诊断类组合工具 | §4.1 同 | §3.3 注册表规模断言 | ✅ |
| 骨干 | JSON-through-C-ABI | §3.1 同 | §3.1 `c_api_ext` | ✅ |
| 只读 | K4 安全约束 | §8 不变式 4 | §3.4 `readOnlyHint=True` | ✅ |

### 1.5 Step 8 发现与机械核对记录

**核对命令**（可重放）：模块名逐一比对（DESIGN §3.2 表 vs IMPL `### 3.x`）；不变式 `^[0-9]+\. \*\*` vs `^\| [0-9]+ ` 计数；接口 `^### 4\.\d` 双侧列举。

| 严重性 | 发现 | 证据 | 处置建议 |
|---|---|---|---|
| **P3（F4）** | 四文档**接口小节编号非 1:1**：DESIGN §4.1（工具契约）由 IMPL §4.3（工具声明）承接；IMPL §4.1（DiagnosticEnvelope）实为 DESIGN §6.1 内容 | 双侧 `### 4.x` 列举（§1.5 命令） | 二选一：(a) 在 IMPL §4 增「↔ DESIGN 映射表」；(b) 统一编号。**已按 (a) 处置**（IMPL §4 首部映射表，改动最小，不动已定稿 DESIGN）；**状态：已闭环** |
| 注 | DESIGN §4.3 与 IMPL §3.4 的 SDK 形态已在 Step 6 双向对账订正（F3） | DESIGN §4.3 订正记录 | 已闭环 |
| **P2（F5）核对维度缺失** | Step 8 核对表**未含**「出参粒度语义」维度 → §4.1 的一揽子（7/7 多统计量）与 §6.1 的单方法信封的**粒度冲突漏检**；Step 9 编写 registry handler 时触发发现 | DESIGN §4.1 出参要点列 vs §6.1 字段（7/7 复核，复算命令见 MATRIX §9） | **补入核对项**（见 §1.2 第 3 行）；粒度冲突经 Type-1 七维矩阵裁定为**方案 (a) 容器层**，并回写 DESIGN §6.2 / DESIGN §4.1 / IMPL §4.1 → **已闭环**（[DECISION_MATRIX.md](./DECISION_MATRIX.md) §9–§12） |
| **P3（F6）Conclusion 方向语义** | `Conclusion.n_reject` 跨方法计数会**错误**（KPSS 的 H0 是平稳性，与 ADF 族相反，U12） | 读码 E3：`kpss_test.hpp` L92-95、`unit_root` U12 陷阱 | **实现前**细化为方向显式（`primary_method`/`reject_null`/`n_total`），回写 DESIGN §6.2 → **已闭环** |
| **P3（F7）入参词汇未核验** | DESIGN §4.1 初版入参为拟写值，与库内实测签名不符（7/7 工具皆有差异，如 `spec.model:"ctt"` 不存在、variance-ratio 无 `bootstrap`、vecm 缺必填 `rank`） | 逐签名 E3 读码（绑 commit `4f6c74f`） | 入参列**整体重写**为实测词汇 + 加注 F7 对照表，IMPL §4.3 声明同步更正 → **已闭环** |

## 2. 功能验收

### 2.1 `json_serializer`（IMPL §3.2）

| 验收项 | 测试方法 | 通过条件 | 状态 | 证据 |
|---|---|---|---|---|
| double 往返位精度 | `test_json_min.cpp` 往返用例 | `memcmp` 位级相同 | ✅ | `test_json_min.exe` 14/14 通过（2026-10-18）；位比较用 `std::bit_cast<uint64_t>` 等价位级比对；另有 Python `repr()` 外部 oracle 8 组逐字符一致 |
| 嵌套对象/数组 | 同上 | 结构与值均一致 | ✅ | 同上；含 `dump→parse→dump` 规范化稳定性 |
| 字符串转义/UTF-8 | 同上 | 转义与多字节正确 | ✅ | 同上；含 `\uXXXX`、代理对 U+1F600、非 ASCII 直通 |
| null/布尔/空数组边界 | 同上 | 不崩溃且语义正确 | ✅ | 同上 |

### 2.2 `tool_registry`（IMPL §3.3）

| 验收项 | 测试方法 | 通过条件 | 状态 | 证据 |
|---|---|---|---|---|
| 7 工具注册完整 | `test_tool_registry.cpp` | 名称集合 = DESIGN §4.1 七个 | ✅ | `test_tool_registry.exe` 14/14 通过（2026-10-18）；顺序与 DESIGN §4.1 表一致；`tool_count()==7` |
| 名称唯一 + 字符集合规 | 同上 | 无重复；正则 `[A-Za-z0-9_.-]{1,128}` | ✅ | 同上；`std::regex_match` 7/7 通过 + `set` 去重断言 + `tool_count() <= 10`（不变式 3） |
| schema 哈希快照锁定 | 同上 | 快照与登记一致（防投毒式变更） | ✅ | 同上；FNV-1a64(name‖description‖inputSchema) 七值冻结在测试内（`e14c201b5ae55741` … `9c2f7816784078c5`） |

### 2.3 `c_api_ext`（IMPL §3.1）

| 验收项 | 测试方法 | 通过条件 | 状态 | 证据 |
|---|---|---|---|---|
| `tool_list` 返回完整 schema JSON | `test_c_api_tools.cpp` | 可解析且含 7 工具 | ☐ | |
| `tool_call` 闭环（逐工具） | 同上 | 7 工具各至少 1 正常调用返回信封 | ☐ | |
| 输出与既有 6 符号共存不冲突 | 同上 + 链接测试 | 符号解析无冲突 | ☐ | |

### 2.4 `mcp_server` 壳（IMPL §3.4）

| 验收项 | 测试方法 | 通过条件 | 状态 | 证据 |
|---|---|---|---|---|
| 壳可启动（stdio 子进程） | `pytest test_mcp_tools.py` | 进程启动且响应 `tools/list` | ☐ | |
| 7 工具经壳暴露 | 同上 | 列表名称与注册表一致 | ☐ | |
| 类型化返回生成 `structured_content` | 同上 | 信封字段完整出现 | ☐ | |
| 壳内无业务逻辑 | 代码审查 | 转发层无数值计算 | ☐ | |

### 2.5 `input_validator`（IMPL §3.5）

| 验收项 | 测试方法 | 通过条件 | 状态 | 证据 |
|---|---|---|---|---|
| 数组长度超限拒绝 | 边界用例 | 长度 = `kMaxSeriesLength+1` 被拒 | ✅ | `test_tool_registry.cpp::RejectsInvalidInput`；>10000 元素返回错误含 "max length"（C++ 侧校验在 `cpphub_tool_registry.hpp::detail` 内联实现——IMPL §2.3 未列 `input_validator` 文件，实际落地于注册表参数校验） |
| 非有限值（NaN/Inf）拒绝 | 同上 | 返回 INVALID_ARG | ✅ | 同上；`is_finite` 检查 → 错误信息含 "non-finite"（Python 侧 INVALID_ARG 由 S3 ABI 层映射） |
| 必填缺失/类型错误 | Pydantic 用例 | `ValidationError` | ☐ | C++ 侧已拒绝（缺必填/类型不符 → `ToolOutcome.error`）；**Pydantic 侧待 S4 壳** |

## 3. 接口验收

### 3.1 `cpphub_v1_tool_list`（DESIGN §4.2）

| 验收项 | 通过条件 | 状态 |
|---|---|---|
| 签名与 DESIGN §4.2 一致 | 参数/返回类型完全一致 | ☐ |
| 正常调用返回可用 JSON | 可被解析且非空 | ☐ |
| 缓冲不足返回 `CPPHUB_ERR_BUFFER_TOO_SMALL` | 返回 -4 且回填所需长度 | ☐ |
| 二次分配后调用成功 | 大缓冲返回 `CPPHUB_OK` | ☐ |

### 3.2 `cpphub_v1_tool_call`（DESIGN §4.2）

| 验收项 | 通过条件 | 状态 |
|---|---|---|
| 签名与 DESIGN §4.2 一致 | 参数/返回类型完全一致 | ☐ |
| 正常输入返回正确信封 | 字段与诊断结果一致 | ☐ |
| 未知工具名返回 INVALID_ARG | 错误码正确 | ☐ |
| 缓冲不足返回 -4 且可重试 | 语义与 `tool_list` 一致 | ☐ |

### 3.3 MCP 工具声明（DESIGN §4.1）

| 验收项 | 通过条件 | 状态 |
|---|---|---|
| 7 工具名与入参 schema 与 DESIGN §4.1 一致 | 逐项比对 | ☐ |
| 描述含文献陷阱注释 | GR6/CI1 等注释出现在 description | ☐ |
| 只读注解生效 | `readOnlyHint=True` | ☐ |
| 工具名默认取函数名之覆盖生效 | 名称与预期一致 | ☐ |

## 4. 不变式验收

| 不变式（DESIGN §8） | 验证方法 | 状态 |
|---|---|---|
| 1 无状态 | 同一调用连续两次结果逐位相同 | ✅ |
| 2 工具名唯一 + 字符集 | 注册表快照断言 | ✅ |
| 3 工具面 ≤10 | 单测断言 `size() <= 10` | ✅ |
| 4 只读 | 注解快照 + 无副作用审查 | ☐ |
| 5 版本协商可达 | 联调（U-4）modern/legacy 双路径 | ☐ |
| 6 出参 schema 只增不改 | 信封字段快照锁定 | ☐ |

## 5. 错误处理验收

| 错误场景（IMPL §6） | 触发方式 | 预期行为 | 状态 |
|---|---|---|---|
| 工具名未知 | 传不存在名称 | `CPPHUB_ERR_INVALID_ARG` | ☐ |
| 参数类型/范围非法 | 传错类型 | `-32602` / `ValidationError` | ☐ |
| 缓冲不足 | 传 0 容量 | `-4` + 回填长度 | ☐ |
| 序列超长 | >10k 数组 | INVALID_ARG + 可操作提示 | ☐ |
| 协议版本不匹配 | 请求未支持版本 | `UnsupportedProtocolVersionError` | ☐ |
| 内核运行时异常 | 注入异常 | `CPPHUB_ERR_RUNTIME` + last_error | ☐ |

## 6. 性能验收

| 性能指标 | 目标 | 实测 | 状态 |
|---|---|---|---|
| `tool_list` 首次响应 | < 5 ms（含懒构建） | | ☐ |
| `tool_call` 透传开销（1k 点） | < 1 ms（不含诊断本体） | | ☐ |
| JSON 往返 10k 点数组 | < 20 ms | | ☐ |

## 7. 兼容性验收

| 环境 | 版本 | 状态 | 说明 |
|---|---|---|---|
| Python | 3.11 | ☐ | CI 目标下限（U-7 取严） |
| Python | 3.12 | ☐ | |
| OS | Windows 2022 (MSVC) | ☐ | 既有 CI 矩阵 |
| OS | Ubuntu (GCC) | ☐ | 既有 CI 矩阵 |

## 8. ADD 审计（Step 10）

### 8.1 Spec 质量门（Phase 0）

| 维度 | 得分（0-1） | 说明 |
|---|---|---|
| 可测试约束 | | 验收项均含可测通过条件 |
| 模块映射 | | IMPL §3 六模块 ↔ DESIGN §3.2 |
| 接口契约 | | DESIGN §4 ↔ IMPL §4 逐项 |
| 修正项 | | F1/F2/F3 均已修正并回写 |
| 跨模块契约 | | ABI ↔ 壳 ↔ 注册表 三层接口 |
| **总分** | /5 | **档位**: |

### 8.2 ADD 审计发现

| 严重性 | 发现 | 证据 | 修复建议 |
|---|---|---|---|
| P1 | （实现后填写） | | |
| P2 | （实现后填写） | | |
| P3 | （实现后填写） | | |

### 8.3 ADD Iron Law 检查

| 检查项 | 状态 |
|---|---|
| 断言恒真式（无恒真断言冒充验证） | ☐ |
| 单文件检查盲区（跨文件契约已覆盖） | ☐ |
| 设计文档独有约束有测试（K1–K4 / 六不变式） | ☐ |
| 修正阻断性项有测试（F1 新增错误码 -4 / F3 SDK 形态） | ☐ |

## 9. 文档完整性

| 文档 | 存在 | 与实现一致 |
|---|---|---|
| RESEARCH.md | ✅ | ☐ |
| DESIGN.md | ✅ | ☐ |
| IMPLEMENTATION.md | ✅ | ☐ |
| CHECKLIST.md（本文件） | ✅ | ☐ |
| DECISION_MATRIX.md | ✅ | ☐ |
| 开发日志（DEVELOPMENT_LOG.md） | ☐ | ☐ |

## 10. 验收结论

### 10.1 验收统计

| 类别 | 总数 | 通过 | 失败 | 待办 |
|---|---|---|---|---|
| 文档一致性 | 14 | 14 | 0 | 0 |
| 功能 | 17 | 9 | 0 | 8 |
| 接口 | 12 | 0 | 0 | 12 |
| 不变式 | 6 | 3 | 0 | 3 |
| 错误处理 | 6 | 0 | 0 | 6 |
| 性能 | 3 | 0 | 0 | 3 |
| 兼容性 | 4 | 0 | 0 | 4 |
| ADD 审计 | 4 | 0 | 0 | 4 |
| **总计** | **66** | **26** | **0** | **40** |

> **统计口径（RULE-2）**：上表来自本文件**逐项核对**（每项含独立章节与状态单元格 `☐`/`✅`/`⚠️`），非事后汇总推算。
> **计数可重放**（E1）：状态单元格逐节计数 = 文档一致性 §1(14) / 功能 §2(17) / 接口 §3(12) / 不变式 §4(6) / 错误处理 §5(6) / 性能 §6(3) / 兼容性 §7(4) / ADD 审计 §8.3(4)，合计 **66**；§9(12) 文档完整性 / §10.2(3) 验收决定 / §11(6) 后续行动 **不计入**本表八类。
> **复算命令**：`(Get-Content spec\mcp-server\CHECKLIST.md | Select-String -Pattern '☐|✅|⚠️' | Measure-Object).Count`（按 `^## (\d+)\.` 分段归集；§1 已由 Step 8 独立 pass 执行并闭环 F4，Step 9 追加 1 项核对维度后闭环 F5）。
> **规范形态**：章节号固定 `### 10.1 验收统计`；数值单元格为纯十进制整数；逐行 `总数 = 通过 + 失败 + 待办`（14 = 14+0+0；17 = 9+0+8；6 = 3+0+3；…），总计行按列求和（14+17+12+6+6+3+4+4 = 66；通过列 = 26；待办列 = 40）。
> **进行中说明**：Step 9 TDD 逐 S 步回写——**S1**（`json_serializer`）§2.1 四行 ✅；**S2**（`tool_registry`）§2.2 三行 + §2.5 两行 + §4 不变式 1–3 行 ✅；F5 容器层/F6 Conclusion/F7 入参 三处订正已回写四文档。S3–S6 待续。

### 10.2 验收决定

- [ ] **验收通过**：所有 P1 项通过，无阻塞性问题
- [ ] **有条件通过**：条件 = U-4 客户端联调完成 + U-7 Python 下限 CI 确认
- [ ] **验收失败**：原因 =

### 10.3 签字

| 角色 | 签字 | 日期 |
|---|---|---|
| 实施者 | | |
| 审查者 | | |

## 11. 后续行动

| 行动 | 责任人 | 期限 | 状态 |
|---|---|---|---|
| Step 8 四文档一致性独立 pass | — | — | ✅ |
| Step 9 TDD 实现（S1–S5） | | | ☐ |
| U-4 真实客户端联调 | | | ☐ |
| U-7 Python 下限 CI 确认 | | | ☐ |
| 更新 docs/README.md 索引（spec/ 四件套状态） | | | ☐ |
| 记录 DEVELOPMENT_LOG.md | | | ☐ |

---

**Review 签字**: _________ 日期: _________
