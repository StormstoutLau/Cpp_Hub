# 决策矩阵：MCP Server 未决项裁定（U-1 / U-2 / U-3，追加 F5）

---
id: MCP-SERVER-DECISION-MATRIX
type: design
version: 1.1
status: draft
date: 2026-10-18
depends: [MCP-SERVER-RESEARCH, FWK-ASSERTION, ADR-009, ADR-010]
upstream: null
---

> **用途**：对 [RESEARCH.md](./RESEARCH.md) §6.4 登记的三项未决项做「选项 × 维度」逐格取证裁定（§0–§8）；并追加 Step 9 实施期新发现的 **F5（工具出参形状）** 裁定（§9–§12）。
> **规范依据**：Spec_Workflow `spec/templates/DECISION_MATRIX_TEMPLATE.md` + `SPEC_PROCESS.md`「迭代收敛与未决项登记」节；证据分级沿用 `docs/ASSERTION_EVIDENCE_FRAMEWORK.md`。
> **取证日期**：§0–§8 = 2026-10-17；§9–§12 = 2026-10-18（全部页面级复核；读数日期同此）

## 0. 前置判据与 Type 判定

| 判据 | U-1 版本策略 | U-2 形态终态 | U-3 工具面边界 |
|---|---|---|---|
| 可逆性 | **Type-1**：声明的受支持协议版本 = 与客户端的契约 | **Type-2**：B1 是明确过渡点，替换为 B2 不改变工具契约 | **Type-1**：工具面是 agent 依赖的公开契约，移除工具即破坏兼容 |
| 候选数 | 2（硬绑 2026-07-28 / dual-era） | 2（B1 Python 壳 / B2 C++ 原生） | 2（诊断类限定 / 诊断+定价全暴露） |
| 必要性 | 选错 → legacy 客户端直接失败（返工） | 无返工点（先做再改） | 选错 → 工具面膨胀后收缩 = 破坏兼容 |

**本轮判定**：U-1、U-3 为 Type-1 → 全维度取证；U-2 为 Type-2 → 矩阵仅作依据记录（执行「先做再改 + 记录决策」）。

---

## 1. 证据锚点（页面级复核，A 类）

| ID | 断言 | URL | 原文引文（≤3 行） | 级别 |
|----|------|-----|------------------|------|
| **E-1** | 早期修订用 `initialize` 握手建立会话；契约由"双方协商单版本"定义 | `modelcontextprotocol.io/specification/latest/basic/lifecycle` → 重定向至 `.../2026-07-28/basic/versioning` | `**Legacy** : protocol versions that establish a session with an initialize handshake (2025-11-25 and earlier).` | 页面级 |
| **E-2** | 2026-07-28 及以后为"现代"修订：**无握手**，版本逐请求声明 | 同上 | `**Modern** : protocol versions that convey version, identity, and capabilities as per-request metadata (revision 2026-07-28 and later).` / `There is no negotiation handshake. Every request carries its protocol version...` | 页面级 |
| **E-3** | 版本不匹配时服务端**必须**回 `UnsupportedProtocolVersionError` 并列出受支持版本，客户端据此重试 | 同上 | `it **MUST** respond with an UnsupportedProtocolVersionError listing the versions it does support` / `"supported": ["2026-07-28", "2025-11-25"]` | 页面级 |
| **E-4** | **Legacy 客户端 × Modern 服务端 = 失败**（无 fall-forward 机制） | 同上（Compatibility Matrix） | `Legacy \| Modern \| Fails. stdio: the server rejects initialize with a JSON-RPC error...` / `Legacy clients have no fall-forward mechanism.` | 页面级 |
| **E-5** | 官方文档 Versioning 页文字标称 current = 2025-11-25（与 `/specification/latest/` 解析结果不一致） | `modelcontextprotocol.io/specification/versioning` | `The **current** protocol version is [**2025-11-25**]` | 页面级（**冲突登记**） |
| **E-6** | 2026-07-28 已正式发布（含无状态核心 / MRTR / 可缓存列表 / 授权加固 / 扩展框架 / Tier 1 SDK 更新） | `blog.modelcontextprotocol.io/posts/2026-07-28/` | `The 2026-07-28 Model Context Protocol specification is out, bringing a stateless protocol core, Multi Round-Trip Requests, header-based routing, cacheable list results, authorization hardening, a formal extensions framework, and updated Tier 1 SDKs.` | 页面级 |
| **E-7** | 官方 SDK 仅 10 个语言，Tier 1 = TypeScript / Python / C# / Go；**无 C++ 官方 SDK** | `modelcontextprotocol.io/docs/sdk` | `\| Python \| modelcontextprotocol/python-sdk \| Tier 1 \|`（表中无 C++ 行） | 页面级 |
| **E-8** | stdio 传输 = 本地同机子进程通信，**无网络开销**，性能最优 | `modelcontextprotocol.io/docs/2026-07-28/learn/architecture` | `**Stdio transport** : Uses standard input/output streams for direct process communication between local processes on the same machine, providing optimal performance with no network overhead.` | 页面级 |
| **E-9** | 存在零依赖原生 C/C++ 社区 SDK（纯 stdio，无 Python/Node 运行时） | `ibm.github.io/tsar-mcp/` | `This SDK is a zero-dependency C/C++ framework for building Model Context Protocol (MCP) servers as native executables — without Python, Node.js, package managers, or background runtimes.` | 页面级 |
| **E-10** | 「包装既有 API」是官方点名的工具设计**反模式** | `anthropic.com/engineering/writing-tools-for-agents` | `More tools don't always lead to better outcomes. A common error we've observed is tools that merely wrap existing software functionality or API endpoints—whether or not the tools are appropriate for agents.` | 页面级 |
| **E-11** | 工具过多/重叠会分散 agent；工具应**合并多步操作** | 同上 | `Too many tools or overlapping tools can also distract agents from pursuing efficient strategies.` / `Tools can consolidate functionality, handling potentially multiple discrete operations (or API calls) under the hood.` | 页面级 |
| **E-12** | 社区范式 = 官方 reference server 均为**窄域工作流**（7 个：Everything/Fetch/Filesystem/Git/Memory/Sequential Thinking/Time） | `modelcontextprotocol.io/examples` | `[**Git**](...) - Tools to read, search, and manipulate Git repositories`（其余同构，共 7 个） | 页面级 |
| **E-13** | 社区量化类 MCP server 至少 5 个（分析/决策辅助型工具面） | 检索级（HPSILab / Quant Research Assistant / tradingview-mcp / MonteWalk / FactorHub） | 未做页面级复核 | ⚠️ 检索级 |

---

## 2. U-1 矩阵：规范版本策略（Type-1）

| 维度 | 选项 A：硬绑 2026-07-28（modern-only） | 选项 B：dual-era（2026-07-28 modern + 2025-11-25 legacy） |
|------|--------------------------------------|--------------------------------------------------------|
| **成本** | 低——单一代码路径；用官方 SDK 默认行为（E1/E7；2026-10-17） | 中——需处理 `initialize` legacy 分支与 `server/discover` 探测（E-1/E-2；2026-10-17） |
| **收益** | 中——实现最简、贴合最新原语 | 高——覆盖存量 legacy 客户端，避免 E-4 的"直接失败"（页面级；2026-10-17） |
| **可维护性** | 中——单路径易维护，但客户端生态变化需再适配 | 中——两条路径，但官方已定义分代语义与探测规则（E-1/E-2；2026-10-17） |
| **兼容** | **低**——E-4 明确：Legacy client × Modern server = Fails（页面级） | **高**——官方 Compatibility Matrix 全档位可达 Work（E-4；2026-10-17） |
| **升级潜力** | 高——天然对齐后续 modern 修订 | 高——dual-era 服务端可"同时服务两代"（E-4 末段） |
| **短期** | 中——可快速出可用 server | 中——多一层探测逻辑，首轮略慢 |
| **长期** | 中——须回头补 legacy 支持（返工） | 高——一次到位，无回退点压力 |

**格内复算说明**：本表无计数/占比类声明；全部为定性评估 + 证据等级，符合 R7「无法复算者不得写入格内」。

---

## 3. U-2 矩阵：形态终态（Type-2，依据记录）

| 维度 | 选项 A：B1 Python 薄壳（官方 Tier 1 SDK） | 选项 B：B2 C++ 原生 stdio 单 exe（社区 SDK） |
|------|------------------------------------------|---------------------------------------------|
| **成本** | 低——复用现有 nanobind；官方 SDK 开箱（E-7；2026-10-17） | 中-高——无官方 C++ SDK，须引社区实现（E-7/E-9；2026-10-17） |
| **收益** | 高——最快验证工具设计；tool schema 可迭代 | 高——零 Python 分发、无网络端口（E-9；2026-10-17） |
| **可维护性** | 高——Tier 1 官方维护承诺（E-7） | 中——社区 SDK 维护承诺弱于 Tier 1（E-7 vs E-9） |
| **兼容** | 高——Tier 1 SDK 协议合规 + 类型安全（E-7） | 中——协议合规由社区实现保证 |
| **升级潜力** | 高——随官方 SDK 升级（E-6「updated Tier 1 SDKs」） | 中——需跟随规范自行升级 |
| **短期** | **高**——立即可用 | 低——需先做协议层/工具注册验证 |
| **长期** | 中——与"header-only C++"品牌叙事有张力 | **高**——单二进制分发，契合 C++ 库定位 |

**结论**：Type-2 → **先做 B1，B2 为终态**；切换触发条件见 §5。

---

## 4. U-3 矩阵：工具面边界（Type-1）

| 维度 | 选项 A：仅诊断类组合工具（≤10） | 选项 B：诊断类 + 定价类全暴露 |
|------|-------------------------------|----------------------------|
| **成本** | 低——工具数少、schema 小 | 中——须为众多定价函数设计 schema |
| **收益** | 高——命中差异化资产（`TestResult` + D1–D23 陷阱注释） | 中——覆盖面广但同质化（BS/Heston 定价 1998 年即有） |
| **可维护性** | 高——面窄，契约稳定 | 低——定价 API 演进即需改工具面 |
| **兼容** | 高——≤10 工具，agent 选择负担小（E-11） | 低——工具重叠/过多分散 agent（E-11 页面级） |
| **升级潜力** | 高——诊断类可增量扩展 | 中——工具面膨胀后收缩 = 破坏兼容 |
| **短期** | 高——先交付 7 个诊断工具 | 中——首轮工作量大 |
| **长期** | 高——契合 E-10（避免"包装既有 API"反模式）与 E-12（窄域范式） | 低——落入 E-10 点名的反模式 |

---

## 5. 收敛出口（本轮结论）

**本轮结论（明确建议）**：

1. **U-1 → 选 B（dual-era），目标修订 = `2026-07-28`（modern），并接受 `2025-11-25` 及更早（legacy `initialize` 握手）**。一句话理由：页面级 Compatibility Matrix 明确 `Legacy client × Modern server = Fails` 且 legacy 客户端**无 fall-forward 机制**，硬绑单版本会造成存量客户端直接不可用；dual-era 是官方定义的合法且完备的形态。
   - **首轮落地**：直接用官方 **Python Tier 1 SDK**，由其承载协商（E-7）；Cpp_Hub 侧只需声明受支持版本列表。
2. **U-2 → 选 A 先行（B1 Python 薄壳），B 为终态（B2 C++ 原生 stdio）**。一句话理由：官方**无 C++ SDK**（E-7），原生路径必须依赖社区 SDK（E-9，维护承诺弱于 Tier 1），故先以 Tier 1 路径低成本验证工具契约，再在 schema 冻结后切换。
3. **U-3 → 选 A（仅诊断类组合工具，≤10）**。一句话理由：官方将"包装既有 API"点名为反模式（E-10）、工具过多会分散 agent（E-11）、官方 reference servers 全为窄域工作流（E-12）；定价类是同质化能力，诊断类是本库差异化资产。

**未闭合的未决项**（四元块）：

- **U-4 客户端实际支持矩阵**
  - **判据**：≥1 个真实客户端（Trae / Claude Desktop / Cursor）完成 `tools/list` + `tools/call` 闭环，且协商路径被日志证实为 modern 或 legacy
  - **触发条件**：首次端到端联调（可机械识别：客户端日志出现 `UnsupportedProtocolVersionError` 或 `initialize`）
  - **当前倾向**：以官方 SDK 默认 dual-era 行为承载，先按 legacy 兼容设计
  - **证据等级**：待核（无本机客户端实测）

- **U-5 官方文档 Versioning 页与 `/specification/` 修订状态不一致的最终定性**
  - **判据**：官方站点两者一致，或 SEP/官方回复界定"current"与"latest"语义差异
  - **触发条件**：下次规范迭代（版本号变更）时复查（可机械识别：versioning 页文字变化）
  - **当前倾向**：以 `/specification/latest/` 解析结果与 spec 内 Modern/Legacy 定义为实施依据（E-1/E-2/E-6）
  - **证据等级**：E1（页面级，冲突已登记）

---

## 6. 幻觉排除与订正记录（本轮实证）

| # | 先前声明 | 复核结果 | 处置 |
|---|---------|---------|------|
| **H-1** | "2026-07-28 为现行正式规范"（原 A9，**检索级**） | 页面级复核发现官方 *Versioning* 页文字写 `current ... 2025-11-25`（E-5），与 `/specification/latest/` 解析到 `2026-07-28`（E-1/E-2）**不一致** | **冲突登记**；改以 spec 的 Modern/Legacy 定义 + `/latest/` 解析为实施依据；差异列入 U-5 待官方澄清 |
| **H-2** | "Claude Code 工具响应默认上限 25,000 tokens"（原 A14，**检索级**） | 页面级 WebFetch + grep **未命中**该数值 | **降级为未验证**，**不写入**建议依据（从论证中剔除） |
| **H-3** | "官方 SDK 已更新 Tier 1"（表述含混） | 页面级 SDK 表核实 Tier 1 = TS/Python/C#/Go（E-7） | 精确化；并**新增事实**：官方 SDK **无 C++** |
| **H-4** | "C++ 原生 MCP 可行"（原 B3，**推断**） | 页面级证据支持"可行"（E-9），但补充"无官方 C++ SDK、须依赖社区实现" | 精确化：可行性成立，**维护风险上修** |
| **H-5** | 社区量化 MCP server 工具面"以诊断/分析为主"（原 B2/B4） | 仍为**检索级**，未做页面级复核 | 如实标注 ⚠️ 检索级；**不单独支撑** U-3 结论（U-3 由 E-10/E-11/E-12 页面级证据支撑） |

> **纪律声明**：按 `ASSERTION_EVIDENCE_FRAMEWORK` §3「版本号/章节号/数值常量三类各有链接」要求，H-1/H-2 两条正是"数值/版本凭印象"高发位，本轮各捕获 1 例——**页面级复核是唯一可靠拦截层**。

---

## 7. 社区案例（可复用模式）

| 案例 | 模式 | 对本决策的作用 |
|------|------|---------------|
| 官方 reference servers（7 个，E-12） | 每 server 只覆盖一个窄域工作流 | 支持 U-3 选 A（诊断域 = 一个窄域） |
| 官方 SDK 分层（Tier 1–3，E-7） | 官方维护承诺分级；Python 为 Tier 1 | 支持 U-2 选 A（Tier 1 低风险路径） |
| IBM TSAR-MCP（E-9） | 零依赖原生 C/C++ + 纯 stdio + ssh 远端 | 证明 U-2 的 B2 终态技术可行 |
| Streamable HTTP 传输（E-8 同页） | 远端传输需 OAuth 等鉴权 | 支持 RESEARCH §C4「远程 HTTP 延后」 |
| OWASP MCP 安全清单（检索级） | 工具响应视为不可信数据 | 支撑 K4 安全约束（待补页面级） |

---

## 8. 反模式与诚实边界

**反模式（显式避免）**：
- 分析瘫痪：本轮 **关闭 U-1/U-2/U-3 三项**，仅新开 U-4/U-5 **且均带触发条件** → 构成「迭代净进展」。
- criteria gaming：矩阵格只作依据，**结论未被指标优化所替代**——三项结论均以"要不要返工 / 会不会破坏兼容"为唯一判据。

**诚实边界**：
- 本矩阵为**过程性形态**：无门禁、无校验器；可机械校验的仅"格是否齐全 / 证据等级是否非空 / 触发条件是否存在 / 引用是否完整"，**不核结论对错**——结论正确性留人工裁决（Scott）。
- E-13 与 OWASP 条目仍为**检索级**，若进入 DESIGN 阻断性引用须先补页面级复核。

---

## 9. F5 问题陈述（Step 9 实施期发现，P2）

**发现**：DESIGN §4.1 的「出参要点」列与 DESIGN §6.1 的 `DiagnosticEnvelope` **粒度互斥**——前者 7/7 工具皆为**多统计量一揽子**，后者为**单方法单标量**信封。

| 出处 | 逐字内容（读数 2026-10-18） |
|---|---|
| DESIGN §4.1 出参要点 | `cpphub_unit_root_diagnose` → `ADF/PP/KPSS/DF-GLS/NP/ZA 一揽子 + 综合结论`；`cpphub_cointegration_scan` → `EG + Johansen 迹/λmax + rank 判定`；`cpphub_vecm_fit` → `α/β/ECT t + EM2002 临界值 + Π=αβ′`；`cpphub_dy_spillover` → `TCI/TO/FROM/NET 矩阵 + 滚动选项`（其余 3 个同构） |
| DESIGN §6.1 | `class DiagnosticEnvelope: tool / method / statistic / p_value / cv_5pct / reject_null / notes`（**单** method、**单** statistic） |
| IMPL §4.1 | 声明 `**签名一致性**：与 DESIGN §6.1 一致 ✅` —— 该声明在**粒度**意义上失实 |

**7/7 复核命令（E1，可重放）**：`Select-String -Path spec\mcp-server\DESIGN.md -Pattern '^\| .cpphub_.*\|' | Measure-Object` → 7 行；逐行比对"出参要点"列是否含 `+`/`一揽子`/`矩阵` 等并列语义 → **7/7 为多统计量**。

**Step 8 未捕获的原因（我的核对缺项，已登记）**：CHECKLIST §1.2 只核"小节编号映射是否 1:1"，**未核"出参粒度语义"**——属核对表维度缺失，非他人文档缺陷。

**可逆性判定**：工具出参形状 = **agent 依赖的公开契约**（发布后变更即破坏兼容，且客户端 SHOULD 校验 schema，见 E-15）→ **Type-1**，须全维度取证。

---

## 10. F5 补充证据锚点（页面级复核，2026-10-18）

| ID | 断言 | URL | 原文引文（≤3 行） | 级别 |
|----|------|-----|------------------|------|
| **E-14** | `structuredContent` 可为**任意 JSON 值**（对象/数组/字符串/数字/布尔/null），须符合 `outputSchema`；官方并给出**数组根** outputSchema 示例（`list_users`） | `modelcontextprotocol.io/specification/2026-07-28/server/tools` | `**Structured** content is returned as a JSON value in the structuredContent field of a result. This can be any JSON value (object, array, string, number, boolean, or null)...` / `Example tool with array output schema: {"type": "array", "items": {...}}` | 页面级（官方，A） |
| **E-15** | 声明 `outputSchema` 时：服务端 **MUST** 输出符合它；客户端 **SHOULD** 校验 | 同上（Output Schema 节） | `Servers **MUST** provide structured results that conform to this schema.` / `Clients **SHOULD** validate structured results against this schema.` | 页面级（官方，A） |
| **E-16** | 组合关键字（`anyOf`/`oneOf`/`allOf`/`if`/`then`/`else`）**"校验代价高"**，实现 SHOULD 施加边界（深度/子 schema 数/时间预算）以防校验器 DoS | `modelcontextprotocol.io/specification/2026-07-28/basic#composition-keyword-resource-use` | `Composition keywords (anyOf, oneOf, allOf, if/then/else) and $defs enable expressive schemas but can be expensive to validate. Implementations SHOULD apply reasonable bounds... to prevent a malicious schema from acting as a Denial-of-Service vector against the validator.` | 页面级（官方，A） |
| **E-17** | 协议自身对"多形态结果"的官方模式 = **统一容器 + 显式判别串**，而非 schema 层条件分支 | `modelcontextprotocol.io/specification/2026-07-28/basic`（ResultType 节） | `MCP supports polymorphic result types, allowing servers to return different structures based on the outcome of the request. The resultType field is a string that clients can use to determine how to parse and handle the result object.` | 页面级（官方，A） |
| **E-18** | 第三方 MCP 框架（ActionMCP）明确：**inputSchema 须扁平**，而 **outputSchema 支持无限制嵌套（含对象数组）**；并警示"SDK 严格度在上升，仅靠 LLM 客户端测试会得到假信心，须用 MCP Inspector 验证" | `raw.githubusercontent.com/seuros/action_mcp/master/TOOLS.MD` | `\| Nesting Level \| Flat only (single level) \| Unlimited nesting depth \|` / `**SDK strictness is increasing** - What works today with "vibed" clients may break tomorrow with strict validation` | 页面级（第三方，B） |
| **E-19** | 社区惯例：工具结果为 `status`/`data`/`metadata` 的**统一信封**（metadata 与 data 同层） | `theneuralbase.com/mcp/learn/intermediate/query-result-formatting/` | `The article's JSON structure (status/data/metadata) is the recommended MCP standard.` | ⚠️ 检索级（第三方，C） |
| **E-20** | 社区工程指南：**inputSchema** 应"尽可能扁平"（嵌套 → token 增多、LLM 认知负荷、解析失败率上升） | `kansei-link.com/en/insights/mcp-tool-schema-design-guide-2026.html` | `Principle 3: keep inputSchema as flat as possible ... The deeper the nesting: (1) more tokens, (2) more LLM cognitive load, (3) higher parse-failure rate.` | ⚠️ 检索级（第三方，C） |
| **E-21** | SEP-2106：`outputSchema` 由"仅 object"放宽为**任意合法 JSON Schema**；`structuredContent` 同步放宽 | `modelcontextprotocol.io/seps/2106-json-schema-2020-12` | `Update outputSchema types to allow any valid JSON Schema (remove type: "object" constraint)` | ⚠️ 检索级（官方 SEP，B；未做页面级复核） |

**幻觉排除（本轮）**：

| # | 先前/候选声明 | 复核结果 | 处置 |
|---|---|---|---|
| **H-6** | "协议限定 outputSchema 根必须为 object" | 页面级 E-14 给出 **array 根**官方示例；E-21 显示该约束已被移除 | **证伪**，从论证剔除 |
| **H-7** | "方案 (b) 会因 oneOf 而被客户端拒绝" | E-16 只说组合关键字"校验代价高、SHOULD 加边界"（DoS 防护），**并未**说校验器会拒绝 | **修正为**：`oneOf` 引入**校验成本与实现分歧**，而非"被拒"；且 (b) 若把 `methods` 设为可选字段可完全避开 `oneOf` → **该论据不成立**，从 (b) 的减分项中剔除 |
| **H-8** | "E-20 的扁平化建议适用于出参" | E-20 原文明确针对 **inputSchema**；(E-18) 更是把 input（扁平）与 output（可深嵌套）**分列** | **修正为**：E-20 **不构成**出参形状的依据，仅用于确认"输入扁平 / 输出可嵌套"的分野 |

> **纪律声明**：H-7 与 H-8 两条都是"**凭直觉把结论当证据**"的典型——若不做页面级复核，会分别夸大 (b) 的技术缺陷、误用输入侧指南去论证出参。**本轮 (a)/(b) 的真实区分度因此缩小**，裁定改为完全依据 §11 的七维打分。

---

## 11. F5 矩阵：工具出参形状（Type-1）

| 维度 | 选项 (a)：容器层 —— `{tool, results:[envelope…], conclusion}` | 选项 (b)：扁平扩展 —— 信封 + 可选 `methods:[…]` |
|------|-----------------------------------------------------------|------------------------------------------------|
| **成本** | 低——容器是纯增量；信封本体零改写；C++ 侧仅多一层 `json::Value::array` 装配（E-14/E-17；2026-10-18） | 低-中——需**重新定义** §6.1 七个字段在"一揽子"场景下的语义（顶层 `statistic` 何时是检验统计量、何时是"综合结论"），并回写 DESIGN §6.1（E-15；2026-10-18） |
| **收益** | 高——逐行语义**唯一确定**且与 `TestResultBase`（ADR-015）一一对齐；单工具亦走同一路径（`results` 长度 1） | 中——少一层嵌套、多数场景 token 略省；但顶层字段语义随工具漂移 |
| **可维护性** | 高——新增统计量 = 向 `results` 追加一行，**零 schema 变更** | 中——每新增一类"顶层语义"都要扩例，模型与校验逻辑随之分叉 |
| **兼容** | 高——不变式 6「出参 schema 只增不改」天然满足（容器层为新增字段，行内字段不变，E-14/E-15） | **中-低**——顶层字段**语义**变更属"改"而非"增"，违反不变式 6；客户端 SHOULD 校验（E-15）时语义不稳更难严格校验 |
| **升级潜力** | 高——与协议自身判别模式同构（**统一容器 + 显式判别字段**，E-17）；未来加工具只加 `results` 枚举值 | 中——需为每个新工具设计"综合结论"的顶层字段含义 |
| **短期** | 中——首轮多写一层容器与 `conclusion` 固定形状 | **高**——实现最少，一次装配即返回 |
| **长期** | **高**——契约稳定、无需返工 | 低——多统计量语义持续漂移，最终仍需重构为容器 → **返工** |
| **社区案例佐证** | E-18（output 可无限嵌套/对象数组）+ E-19（status/data/metadata 统一信封）→ **主流形态** | 未见同类"变体扁平信封"的社区先例（E-19 的统一信封恰是反例方向） |

**格内复算说明（R7 同构）**：本表无计数/占比类声明，全部为定性评估 + 证据等级。§9 的 `7/7 多统计量` 已附复算命令；`results 长度 1` 为设计约定非计数。

**反模式（显式避免）**：
- **criteria gaming**：不以"写起来少"给 (b) 加分——(b) 的短期优势合法，但其**长期返工**与**不变式 6 破坏**是实质代价；反之 (a) 的短期成本合法，不得因"看起来更规范"而免检。
- **分析瘫痪**：本轮**关闭 F5**，不新开无触发条件的未决项。

**潜在论据的排除记录**（防"为结论找证据"）：H-7（`oneOf` 被拒）与 H-8（扁平化指南适用于出参）**均被证伪并从打分中剔除**；若保留，会人为放大 (a) 的优势——**本轮结论不依赖这两条**。

---

## 12. 收敛出口（本轮结论）

**本轮结论（明确建议）**：

> **F5 → 选 (a) 容器层**：`tool_call` 出参 = `{ "tool": str, "results": [DiagnosticEnvelope…], "conclusion": {…} }`；**`DiagnosticEnvelope` 本体严格保持 DESIGN §6.1 原样**（7 字段、单方法），`results` 长度 ≥1（单统计量工具长度 1），`conclusion` 固定为 `{reject_null: bool, n_reject: int, n_total: int}`。

**一句话理由**：(a) 对 §6.1 **零改写**、逐行与 `TestResultBase` 一一对齐、容器层为纯增量（满足不变式 6），且与协议自身"统一容器 + 显式判别字段"的多形态结果模式（E-17）同构；而 (b) 需要**重新定义**顶层字段在"一揽子"场景的语义（属"改"而非"增"），并已在 H-7/H-8 两条论据被证伪后仍无法回避**长期返工**。

**回写清单（裁定后执行，属 Step 8 追加轮）**：
1. DESIGN §6.1 —— 保留信封原样，**新增** `ToolResult` 容器结构定义。
2. DESIGN §4.1 —— 出参要点列补注"由 §6.1 信封的 `results[]` 数组承载"。
3. IMPL §4.1 —— 更正 `与 DESIGN §6.1 一致 ✅` → `信封与 DESIGN §6.1 一致；容器层为 F5 新增`。
4. CHECKLIST §1.2 —— **新增**核对项「出参粒度语义一致（信封 vs 一揽子）」，并把本次未捕获记为核对维度缺失。

**未闭合的未决项**：**（无）**——F5 本轮闭合；U-4/U-5/U-7 承前述登记，触发条件不变。

**Type 判定依据**：出参形状一经发布即为 agent 依赖的公开契约（客户端 SHOULD 校验 schema，E-15），无无痛回退点 → **Type-1**，故启用全维度矩阵而非"先做再改"。

---

**Review 签字**: _________ 日期: _________

