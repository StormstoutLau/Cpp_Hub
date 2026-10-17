# 调研报告：Cpp_Hub MCP Server（供 Agent 编程工具消费）

---
id: MCP-SERVER-RESEARCH
type: design
version: 1.0
status: draft
date: 2026-10-16
depends: [FWK-ASSERTION, ADR-009, ADR-010, EXTERNAL-INTERFACE-RESEARCH]
upstream: null
---

> **Feature**: mcp-server（Cpp_Hub 对外 MCP 接口层）
> **Spec 步骤**: Step 1–2（调研 → Review 形成调研文档）
> **规范依据**: Spec_Workflow 10 步流程 + [ASSERTION_EVIDENCE_FRAMEWORK.md](../../docs/ASSERTION_EVIDENCE_FRAMEWORK.md)（权威源 = github.com/StormstoutLau/Spec_Workflow）
> **上游**: 前置调研《对外接口扩展》（本地留存） v1.0（§4 JSON-through-ABI 骨干 / §5 P0–P3 路径 / §6 工具设计）
> **状态**: draft — 待 Scott 裁决 §6.4 三项未决项后进入 DESIGN（Step 3）

---

## §0 断言统计表

| 级别 | 计数 | 说明 |
|------|------|------|
| **A 事实类** | 22 | 每条附 URL + 原文引文（≤3 行）；其中「页面级复核」7 条、「检索级」11 条、「代码取证 E1/E3」4 条（标注见 §5.1/§5.2） |
| **B 推断类** | 5 | 附逐步推理链 + 各依赖源 A 级证据；存在性断言 B1/B2 附搜索库名清单 |
| **C 判断类** | 4 | rationale + 假设声明，标 [判断] |
| **假设区** | 3 | 无 URL+引文，禁止进入正文结论（附录 C） |
| 阻断性断言 | 7 | A1/A2/A3/A5/A6/A7/A10 影响 spec 冻结；其中 6 条页面级复核，A3（stdio）为检索级待二核 |

> **计数可重放命令**（E1，供审计端独立重数）——**真值源 = 附录 B 机读登记块**（唯一枚举），正文标签为其引用：
> ```
> (Select-String -Path docs\research\MCP_SERVER_RESEARCH.md -Pattern '"id":"A\d+"' | Measure-Object).Count   # 22
> (Select-String -Path docs\research\MCP_SERVER_RESEARCH.md -Pattern '"id":"B\d+"' | Measure-Object).Count   # 5
> (Select-String -Path docs\research\MCP_SERVER_RESEARCH.md -Pattern '"id":"C\d+"' | Measure-Object).Count   # 4
> ```
> 本表由人工计数，脚本化生成待 P-003 类工具就绪；声明 = 重数纪律。

---

## 1. 调研目标

**核心问题**：
1. **形态选择**：本仓库（header-only C++20 库 + 现有 C ABI/Python 绑定）应如何落地 MCP server？Python 壳、C++ 原生、还是 JSON-through-C-ABI 骨干？
2. **社区可用模式**：MCP 生态（官方规范 / 官方 SDK / 社区 C++ SDK / 量化类 server）提供了哪些可复用的设计模式与反模式？
3. **本仓库适配**：结合本库差异化资产（诊断类 `TestResult` 结构、位精确确定性、23 个 R 源码排幻觉点 D1–D23），MCP 工具面应暴露什么、如何命名、如何控 token？
4. **安全约束**：面向公开仓库 + Agent 消费，MCP 有哪些强制安全要求（prompt injection / tool poisoning / 只读优先）？

---

## 2. 调研方法

### 2.1 使用的工具

| 工具 | 用途 | 查询/目标 |
|------|------|-----------|
| WebSearch | 通用检索 | MCP 2026 规范/传输/安全/工具设计/量化案例（8 组查询） |
| WebFetch | 页面级取证 | modelcontextprotocol.io 规范页 ×2、IBM TSAR-MCP ×1 |
| Read/Glob/Grep | 本仓库代码取证 | `src/c_api/c_api.cpp`、`python/src/cpphub_python.cpp`、`CMakeLists.txt` |
| Read | 规范取证 | `F:\Spec_Workflow\{SPEC_PROCESS.md, AGENTS.md, README.md, spec/templates/*}` |

### 2.2 调研范围

- **时间范围**：以 2026 年现行规范为基准（MCP 规范版本节奏快，2025-03-26 → 2025-06-18 → 2025-11-25 → **2026-07-28**）
- **领域**：MCP 协议规范、MCP SDK 生态、开发者工具类 server、量化金融类 server、Agent 工具设计方法论
- **排除**：非公开/需授权的 SaaS 专有 server 实现细节；远程 HTTP 部署的 OAuth 全流程（本轮仅做可行性判断，见 C4）

---

## 3. 调研发现

### 3.1 MCP 协议现状（权威规范）

#### 3.1.1 无状态核心

- **[A1]** MCP 是**无状态协议**。 [页面级复核]
  - URL: `https://modelcontextprotocol.io/specification/2026-07-28/basic#statelessness`
  - 引文: `The Model Context Protocol (MCP) is a **stateless protocol** : all the information needed to process a request is contained in the request itself.`
- **[A2]** 服务器**禁止**依赖同连接的前序请求建立上下文。 [页面级复核]
  - URL: 同上
  - 引文: `Servers **MUST NOT** rely on prior requests over the same connection to establish context (e.g., capabilities, protocol version, client identity). Every request supplies this metadata in its _meta field.`

> **对本 feature 的直接含义**：MCP server 不应把"大数组状态"寄存在服务端会话中；本仓库诊断函数是**纯函数**（输入序列 → 输出 `TestResult`），与无状态核心天然契合（见 B5）。

#### 3.1.2 传输

- **[A3]** stdio 传输 = 客户端把 server 作为**子进程**启动，经标准流收发 JSON-RPC。 [检索级]
  - URL: `https://modelcontextprotocol.io/specification/2026-07-28/basic/transports/stdio`
  - 引文: `In the stdio transport, the client launches the MCP server as a subprocess.` / `The server reads JSON-RPC messages from stdin and writes JSON-RPC messages to stdout.`
- **[A4]** Streamable HTTP 传输用于远程 server，支持标准 HTTP 鉴权（推荐 OAuth）。 [检索级]
  - URL: `https://modelcontextprotocol.io/docs/2026-07-28/learn/architecture`
  - 引文: `Streamable HTTP transport: Uses HTTP POST for client-to-server messages with optional Server-Sent Events for streaming capabilities.`

#### 3.1.3 工具（Tools）原语

- **[A5]** 工具发现 = `tools/list` 请求。 [页面级复核]
  - URL: `https://modelcontextprotocol.io/specification/2026-07-28/server/tools`
  - 引文: `To discover available tools, clients send a tools/list request.`
- **[A6]** 工具调用 = `tools/call` 请求，参数 `name` + `arguments`。 [页面级复核]
  - URL: 同上
  - 引文: `To invoke a tool, clients send a tools/call request:` / `"params": { "name": "get_weather", "arguments": { "location": "New York" } }`
- **[A7]** 工具名约束：长度 1–128、大小写敏感、仅 `A-Za-z0-9_-` 与 `.`、server 内唯一。 [页面级复核]
  - URL: `https://modelcontextprotocol.io/specification/2026-07-28/server/tools#tool-names`
  - 引文: `Tool names **SHOULD** be between 1 and 128 characters in length (inclusive).` / `Tool names **SHOULD** be unique within a server.`
- **[A8]** 2026-07-28 引入 **input-required 工具结果**（工具可在执行中请求用户补充输入，客户端重试需换 JSON-RPC `id`）。 [页面级复核]
  - URL: `https://modelcontextprotocol.io/specification/2026-07-28/server/tools#input-required-tool-results`
  - 引文: `Note that the JSON-RPC id **MUST** be different between the initial request and the retry.`
- **[A9]** 2026-07-28 规范为**正式发布版**，含无状态核心 / Multi Round-Trip Requests / 可缓存列表结果 / 幻展框架 / 正式弃用策略。 [检索级]
  - URL: `https://blog.modelcontextprotocol.io/posts/2026-07-28/`
  - 引文: `The 2026-07-28 Model Context Protocol specification is out, bringing a stateless protocol core, Multi Round-Trip Requests, header-based routing, cacheable list results, authorization hardening, a formal extensions framework, and updated Tier 1 SDKs.`

### 3.2 SDK 与实现生态

- **[A10]** **IBM TSAR-MCP**：零依赖原生 C/C++ MCP SDK，纯 `stdio`，无需 Python/Node/运行时。 [页面级复核]
  - URL: `https://ibm.github.io/tsar-mcp/`
  - 引文: `This SDK is a zero-dependency C/C++ framework for building Model Context Protocol (MCP) servers as native executables — without Python, Node.js, package managers, or background runtimes.`
  - 引文: `TSAR-MCP takes the opposite approach: pure stdio, paired with standard ssh for remote access. No open network ports.`
- **[A11]** 官方 Python SDK 包名 `mcp`；**FastMCP 已并入官方 SDK**（`pip install "mcp[cli]"` / `pip install fastmcp`）。 [检索级]
  - URL: `https://pypi.org/project/mcp/1.9.4/` / `https://pypi.org/project/fastmcp/1.0/`
  - 引文: `pip install "mcp[cli]"` / `🎉 FastMCP has been added to the official MCP SDK!`
- **[A12]** 官方 **reference servers** 含 Filesystem / Git / Fetch / Memory / Sequential Thinking / Time。 [检索级]
  - URL: `https://modelcontextprotocol.io/examples`
  - 引文: `Filesystem - Secure file operations with configurable access controls` / `Git - Tools to read, search, and manipulate Git repositories`
- **[A13]** 存在 **官方 MCP Registry**（服务发现入口）。 [检索级]
  - URL: `https://registry.modelcontextprotocol.io/`
  - 引文: `Official MCP Registry — Discover Model Context Protocol servers`

### 3.3 工具设计方法论（Agent 消费侧）

- **[A14]** Anthropic 官方工具编写指南：**限制响应 token 量**（Claude Code 默认上限 25,000 tokens）。 [检索级]
  - URL: `https://www.anthropic.com/engineering/writing-tools-for-agents`
  - 引文: `For Claude Code, we restrict tool responses to 25,000 tokens by default.`
- **[A15]** 社区工具设计共识：**动词+宾语命名**、**同族工具共享命名空间前缀**、禁止泛化名（如 `get_data`）。 [检索级]
  - URL: `https://www.gingerlabs.ai/blog/mcp-tool-schema-design`
  - 引文: `Tool name uses verb + object pattern (create_invoice, not invoice)` / `Related tools share a namespace prefix`
- **[A16]** 客户端最佳实践：工具定义应可缓存，避免重复 `tools/list` 往返。 [检索级]
  - URL: `https://modelcontextprotocol.io/docs/2025-11-25/develop/clients/client-best-practices`
  - 引文: `Cache tool definitions — Once fetched from a server, memoize the definition host-side so re-injecting it later doesn't need another tools/list round trip.`

### 3.4 安全约束（公开库 + Agent 消费）

- **[A17]** OWASP MCP Security Cheat Sheet：**把每个工具响应视为不可信数据**（含已批准 server 的响应）。 [检索级]
  - URL: `https://cheatsheetseries.owasp.org/cheatsheets/MCP_Security_Cheat_Sheet.html`
  - 引文: `Treat every tool response as untrusted data, including responses from approved servers.`
- **[A18]** MCP 规范 Security Best Practices：server **禁止**接受非本 server 签发的 token；会话 ID 须安全非确定性。 [检索级]
  - URL: `https://modelcontextprotocol.io/specification/draft/basic/security_best_practices`
  - 引文: `MCP servers MUST NOT accept any tokens that were not explicitly issued for the MCP server.`

### 3.5 社区案例（量化 / 开发者工具）

- **[B1]** **社区存在 ≥4 个独立 C++ MCP SDK 实现**（存在性断言）。
  - 推理链：① 检索 `MCP server C++ implementation SDK`；② 命中并逐一登记：IBM **tsar-mcp**（官方 GitHub `IBM/tsar-mcp`）、**cppmcp**（LobeHub 收录）、**cxxmcp**（`caomengxuan666/cxxmcp`，ServerBuilder 风格 API）、**cpp-mcp**（`hkr04/cpp-mcp`，HTTP 配置式）、以及社区 plugin 架构实现（`larkascending503/mcp_server`，CMake 构建）。
  - 依赖源 A 级证据：见 A10（TSAR-MCP 页面级）；其余为检索级（LobeHub / skillget / project-awesome 收录页）。
  - 结论：C++ 侧 MCP 实现路径**已有多家成熟实现**，非孤例。
- **[B2]** **社区存在多个量化/交易类 MCP server**（存在性断言）。
  - 推理链：① 检索 `quantitative finance MCP server`；② 命中并登记：**Quant Research Assistant**（HMM regime / 回测 / 期权 / 蒙特卡洛 / 配对交易 / 因子模型 / 组合优化）、**HPSILab Quant Finance MCP**（8-tool，期权/隐含波动率/蒙特卡洛/回测）、**tradingview-mcp**（40+ tools，回测/走查验证）、**MonteWalk**（风险指标/蒙特卡洛/回测）、**FactorHub MCP**（因子/回测）。
  - 说明：均为**分析/决策辅助型**工具面，非原始数值计算库包装。

### 3.6 本仓库现状（代码取证，E3 静态读码 + E1 命令）

- **[A19]** C ABI 当前仅 **6 个** `extern "C"` 符号：`cpphub_v1_get_abi_version` / `cpphub_v1_get_last_error` / `cpphub_v1_bsm_price_batch` / `cpphub_v1_heston_price` / `cpphub_v1_bsm_implied_vol` / `cpphub_v1_mc_price`。 [E3]
  - 证据：`src/c_api/c_api.cpp` 逐行读码；头文件 `include/cpphub/c_api/cpphub_c_api.h`
  - 复核命令：`Select-String -Path src\c_api\c_api.cpp -Pattern 'cpphub_v1_' | Measure-Object`
- **[A20]** Python 绑定现有 **58 个** `m.def(`（nanobind），停留在定价层。 [E1]
  - 证据：`python/src/cpphub_python.cpp`
  - 复核命令：`(Select-String -Path python\src\cpphub_python.cpp -Pattern '\bm\.def\(' | Measure-Object).Count`
- **[A21]** 构建开关：`CPPHUB_ENABLE_C_API`（默认 OFF）/ `CPPHUB_ENABLE_PYTHON`（默认 OFF）/ `CPPHUB_ENABLE_CUDA`（默认 OFF）。 [E3]
  - 证据：`CMakeLists.txt` L115 / L126 / L80
- **[A22]** 测试资产：`tests/` 下 **155 个** `test_*.cpp`，`tests/validation/` 目录存在（基准对齐层）。 [E1]
  - 复核命令：`(Get-ChildItem tests -Recurse -Filter "test_*.cpp" | Measure-Object).Count`

---

## 4. 综合分析

### 4.1 关键发现总结

1. **MCP 已是"无状态 + 双传输 + 工具原语"的成熟协议，规范版本 2026-07-28**。[置信度: ★★★★★]
2. **C++ 原生 MCP 已有零依赖先例（IBM TSAR-MCP，纯 stdio）**，因此"零 Python 依赖的单 exe"路径**无需等到后期**。[置信度: ★★★★☆]
3. **官方 SDK 与社区生态均指向"薄壳 + 有语义的工具面"**，而非把库 API 逐函数映射。[置信度: ★★★★☆]
4. **量化类 MCP server 的工具面高度集中在"诊断/分析/决策辅助"**，与 Cpp_Hub 的"验证方法论"定位一致。[置信度: ★★★★☆]
5. **安全是硬约束**：工具响应视为不可信数据、只读优先、窄工具面。[置信度: ★★★★★]

> **B/C 类断言索引**（正文引用锚点，完整定义见附录 B）：
> [B1] C++ MCP SDK 多家实现 · [B2] 量化类 MCP server 多家 · [B3] 原生单 exe 零依赖可行 · [B4] 工具面以诊断/分析为主 · [B5] 无状态核心与本库纯函数契合；
> [C1] 形态路线（B1→B2 + JSON-through-ABI）· [C2] 工具面边界（诊断优先、≤10）· [C3] 安全默认（只读/标记/校验）· [C4] 发布策略（本地 stdio 优先、Registry 登记、远程 HTTP 延后）。

### 4.2 技术 landscape

| 形态 | 代表 | 优点 | 缺点 | 对本库适用性 |
|------|------|------|------|-------------|
| Python 薄壳（B1） | 官方 `mcp` / FastMCP | 落地快、调试易；可直接复用现有 nanobind | 需 Python 环境；与"header-only C++"品牌叙事有张力 | ✅ 首轮验证工具面 |
| C++ 原生 stdio（B2） | IBM TSAR-MCP、cpp-mcp | 零 Python、单 exe、无网络端口 | 需自实现 JSON-RPC/schema 或引第三方；工具面锁定后改动成本高 | ✅ 终态 |
| JSON-through-C-ABI 骨干 | 本仓库 前置调研《对外接口扩展》（本地留存） §4 提案 | 一次定义全语言受益；新功能只改 C++ 注册 | 需先扩 C ABI（现仅 6 函数） | ✅ 推荐底座 |
| 远程 Streamable HTTP | 厂商型 server | 远程共享、可托管 | OAuth/鉴权/托管成本；公开库无服务端 | ⏸ 本轮不做（C4） |

**反模式（调研中明确）**：把库的每个函数平铺成同质工具（工具面膨胀 → 选择困难 + token 膨胀 + 注入面扩大）；工具描述写"是什么"而非"何时用/返回什么/副作用"。

### 4.3 研究空白（本 feature 可填补）

- 现有量化类 MCP server 多为"数据+指标查询"型，**没有**把"统计检验的结论语义（统计量/p 值/临界值/`reject_null`）+ 文献级陷阱注释"作为一等输出——Cpp_Hub 的 `TestResult` 结构 + 23 点排幻觉注释（D1–D23）恰好补此空白。
- 现有 C++ MCP SDK 均未针对"位精确确定性"这一本库独有卖点设计（无"跨平台一致性"暴露面）。

---

## 5. 幻觉抑制审查（Step 2 Review）

### 5.1 文献/来源验证

| 引用 | URL | 验证方式 | 级别 |
|------|-----|---------|------|
| [A1][A2] | spec/2026-07-28/basic | WebFetch 页面 + 逐字比对（L181 / L183） | ✅ 页面级 |
| [A5][A6][A7][A8] | spec/2026-07-28/server/tools | WebFetch 页面 + 逐字比对（L90 / L102 / L177-185 / L124） | ✅ 页面级 |
| [A10] | ibm.github.io/tsar-mcp | WebFetch 页面 + 逐字比对 | ✅ 页面级 |
| [A3][A4] | spec …/transports/stdio · …/learn/architecture | WebSearch 返回片段（未做页面级 grep） | ⚠️ 检索级 |
| [A9] | blog.modelcontextprotocol.io/posts/2026-07-28/ | WebSearch 返回片段 | ⚠️ 检索级 |
| [A11][A13] | pypi.org · registry.modelcontextprotocol.io | WebSearch 返回片段 | ⚠️ 检索级 |
| [A12] | modelcontextprotocol.io/examples | WebSearch 返回片段 | ⚠️ 检索级 |
| [A14][A15][A16] | anthropic.com · gingerlabs.ai · mcp docs | WebSearch 返回片段 | ⚠️ 检索级 |
| [A17][A18] | owasp.org · mcp spec security | WebSearch 返回片段 | ⚠️ 检索级 |

> **纪律说明**（沿用本生态既有教训「检索级 ≠ 页面级复核」）：11 条检索级 A 类未做页面级 grep 复核，若进入 spec 阻断性引用，须先补页面级复核或第二独立源。

### 5.2 技术声明验证

| 声明 | 来源 | 验证状态 |
|------|------|---------|
| MCP 2026-07-28 为现行正式规范 | A9 | ✅（官方 blog，检索级） |
| stdio 为本地进程传输、无需网络端口 | A3 + A10 | ✅ 双源 |
| 官方 Python SDK 即 `mcp` 包、FastMCP 已并入 | A11 | ✅（PyPI 两条目一致） |
| 工具名约束 1–128 / 唯一 | A7 | ✅ 页面级 |
| C ABI 仅 6 函数 | A19 | ✅ E3 读码 + 命令复核 |
| Python 绑定 58 个 def | A20 | ✅ E1 命令 |

### 5.3 待修正项

- [ ] 对 11 条「检索级」A 类补页面级复核（进入 DESIGN 前，阻断性项优先）
- [ ] 核实本仓库 MCP 客户端（Trae / Claude Desktop / Cursor）对 2026-07-28 规范的支持矩阵（假设区 H1）

---

## 6. 对设计的输入

### 6.1 可用的技术方案

1. **底座**：扩展 C ABI 为 JSON-through-ABI（`cpphub_v1_tool_list` + `cpphub_v1_tool_call` + JSON 序列化层），与现有 6 个标量符号并存（前置调研《对外接口扩展》（本地留存） §4）。
2. **形态**：`B1` Python 薄壳先行（复用官方 `mcp` SD​K 与现有 nanobind）→ `B2` C++ 原生 stdio 单 exe（参照 TSAR-MCP 零依赖模式）。
3. **工具面**：**诊断组合优先**——`unit_root_diagnose`（ADF/PP/KPSS/DF-GLS/NP/ZA 一揽子）、`cointegration_scan`、`vecm_fit`、`garch_diagnose`、`granger_causality`、`dy_spillover`、`variance_ratio_test`（详见 前置调研《对外接口扩展》（本地留存） §6）。
4. **输出语义**：每个工具返回结构化结论（`statistic` / `p_value` / `cv_5pct` / `reject_null` / `method`），并携带**文献陷阱注释**（如 GR6 列序方向、CI1 p 与 cv 分列）——description 中随工具下发。

### 6.2 关键约束

- **[K1] 无状态约束**（源自 A1/A2）：不得依赖会话状态；每个请求自包含。
- **[K2] 工具名约束**（源自 A7）：`verb_object` + 命名空间前缀（如 `cpphub_` / `ts_`），1–128 字符，`[A-Za-z0-9_-.]`，server 内唯一。
- **[K3] token 预算约束**（源自 A14）：单个工具响应须设默认上限，长数组须分页/截断。
- **[K4] 安全约束**（源自 A17/A18）：只读优先；工具描述在安装时校验；禁止接受非本 server 签发的 token。

### 6.3 风险

| 风险 | 说明 | 缓解 |
|------|------|------|
| R-1 规范漂移 | MCP 规范迭代快（约 3–5 个月一版） | 锁定实现到 2026-07-28；用官方 SDK 而非手写协议层 |
| R-2 工具面膨胀 | 平铺 155 测试/58 绑定函数 → agent 选择困难 + token 爆炸 | 工具数量硬上限（建议 ≤10）；组合优于单点（C2 判断） |
| R-3 注入面 | 公开库 server 被不可信输入驱动 | 只读 + 输出标记为数据 + 输入 schema 强校验 |
| R-4 双轨漂移 | Python 壳与 C++ 原生两套实现漂移 | JSON-through-ABI 单一骨干，前端仅协议薄壳 |

### 6.4 未决项登记（四元块）

- **U-1 规范版本锁定与客户端兼容**
  - **判据**：目标客户端（≥1 个真实 MCP 客户端）完成 `tools/list` + `tools/call` 闭环且无版本协商失败
  - **触发条件**：首次端到端联调；或客户端侧报 `protocolVersion` 不匹配
  - **当前倾向**：锁定 **2026-07-28**，用官方 SDK 承载版本协商
  - **证据等级**：E3（规范页读码）/ 待核（客户端支持矩阵）

- **U-2 形态终态（B1 → B2 切换时机）**
  - **判据**：B1 完成 7 个诊断工具的端到端调用，且工具 schema 稳定 ≥1 个迭代
  - **触发条件**：工具面冻结 + 需要"零 Python 分发"时
  - **当前倾向**：**先 B1 后 B2**（低成本验证工具设计），底座同为 JSON-through-ABI
  - **证据等级**：C（判断，基于 §4.2 对比表）

- **U-3 是否暴露定价类工具**
  - **判据**：真实 agent 调用日志中出现定价类需求，且诊断类工具未覆盖
  - **触发条件**：≥3 次由 agent 侧提出的定价调用诉求
  - **当前倾向**：**暂不暴露**——优先诊断类（差异化资产），定价类按需增量
  - **证据等级**：B（B2 社区案例归纳）/ C（判断）

---

## 7. 参考文献

- Model Context Protocol. *Overview / Statelessness*（2026-07-28）. https://modelcontextprotocol.io/specification/2026-07-28/basic
- Model Context Protocol. *Tools*（2026-07-28）. https://modelcontextprotocol.io/specification/2026-07-28/server/tools
- Model Context Protocol. *Transports / stdio*. https://modelcontextprotocol.io/specification/2026-07-28/basic/transports/stdio
- Model Context Protocol. *Architecture overview*. https://modelcontextprotocol.io/docs/2026-07-28/learn/architecture
- Model Context Protocol Blog. *The 2026-07-28 Specification*. https://blog.modelcontextprotocol.io/posts/2026-07-28/
- Model Context Protocol. *Security Best Practices*. https://modelcontextprotocol.io/specification/draft/basic/security_best_practices
- Model Context Protocol. *Example Servers*. https://modelcontextprotocol.io/examples
- IBM Research. *TSAR-MCP: Zero-Dependency C/C++ SDK*. https://ibm.github.io/tsar-mcp/
- Anthropic. *Writing effective tools for AI agents*. https://www.anthropic.com/engineering/writing-tools-for-agents
- OWASP. *MCP Security Cheat Sheet*. https://cheatsheetseries.owasp.org/cheatsheets/MCP_Security_Cheat_Sheet.html
- PyPI. *mcp* / *fastmcp*. https://pypi.org/project/mcp/ · https://pypi.org/project/fastmcp/
- 本仓库：对外接口调研（本地留存） v1.0；`ADR-009`（C ABI）；`ADR-010`（Python nanobind）

---

## 附录 B：断言机读登记块

```assertions
[
  {"id":"A1","level":"A","claim":"MCP 是无状态协议","source":"https://modelcontextprotocol.io/specification/2026-07-28/basic#statelessness","quote":"The Model Context Protocol (MCP) is a stateless protocol : all the information needed to process a request is contained in the request itself.","verified":"page-level"},
  {"id":"A2","level":"A","claim":"服务器禁止依赖前序请求建立上下文","source":"https://modelcontextprotocol.io/specification/2026-07-28/basic#statelessness","quote":"Servers MUST NOT rely on prior requests over the same connection to establish context","verified":"page-level"},
  {"id":"A3","level":"A","claim":"stdio 传输=客户端启动子进程, stdin/stdout 收发 JSON-RPC","source":"https://modelcontextprotocol.io/specification/2026-07-28/basic/transports/stdio","quote":"The server reads JSON-RPC messages from stdin and writes JSON-RPC messages to stdout.","verified":"search-level"},
  {"id":"A4","level":"A","claim":"Streamable HTTP 用于远程, 推荐 OAuth","source":"https://modelcontextprotocol.io/docs/2026-07-28/learn/architecture","quote":"Streamable HTTP transport: Uses HTTP POST for client-to-server messages","verified":"search-level"},
  {"id":"A5","level":"A","claim":"工具发现=tools/list","source":"https://modelcontextprotocol.io/specification/2026-07-28/server/tools","quote":"To discover available tools, clients send a tools/list request.","verified":"page-level"},
  {"id":"A6","level":"A","claim":"工具调用=tools/call(name,arguments)","source":"https://modelcontextprotocol.io/specification/2026-07-28/server/tools","quote":"To invoke a tool, clients send a tools/call request:","verified":"page-level"},
  {"id":"A7","level":"A","claim":"工具名 1-128 字符, 大小写敏感, server 内唯一","source":"https://modelcontextprotocol.io/specification/2026-07-28/server/tools#tool-names","quote":"Tool names SHOULD be between 1 and 128 characters in length (inclusive).","verified":"page-level"},
  {"id":"A8","level":"A","claim":"input-required 工具结果; 重试须换 JSON-RPC id","source":"https://modelcontextprotocol.io/specification/2026-07-28/server/tools#input-required-tool-results","quote":"Note that the JSON-RPC id MUST be different between the initial request and the retry.","verified":"page-level"},
  {"id":"A9","level":"A","claim":"2026-07-28 为正式规范版本","source":"https://blog.modelcontextprotocol.io/posts/2026-07-28/","quote":"The 2026-07-28 Model Context Protocol specification is out","verified":"search-level"},
  {"id":"A10","level":"A","claim":"IBM TSAR-MCP 零依赖原生 C/C++ SDK, 纯 stdio","source":"https://ibm.github.io/tsar-mcp/","quote":"This SDK is a zero-dependency C/C++ framework for building Model Context Protocol (MCP) servers as native executables","verified":"page-level"},
  {"id":"A11","level":"A","claim":"官方 Python SDK=mcp 包; FastMCP 已并入","source":"https://pypi.org/project/fastmcp/1.0/","quote":"FastMCP has been added to the official MCP SDK!","verified":"search-level"},
  {"id":"A12","level":"A","claim":"官方 reference servers 含 Filesystem/Git/Fetch/Memory/Time","source":"https://modelcontextprotocol.io/examples","quote":"Filesystem - Secure file operations with configurable access controls","verified":"search-level"},
  {"id":"A13","level":"A","claim":"存在官方 MCP Registry","source":"https://registry.modelcontextprotocol.io/","quote":"Official MCP Registry","verified":"search-level"},
  {"id":"A14","level":"A","claim":"Claude Code 工具响应默认上限 25000 tokens","source":"https://www.anthropic.com/engineering/writing-tools-for-agents","quote":"For Claude Code, we restrict tool responses to 25,000 tokens by default.","verified":"search-level"},
  {"id":"A15","level":"A","claim":"工具名用 动词+宾语 + 命名空间前缀","source":"https://www.gingerlabs.ai/blog/mcp-tool-schema-design","quote":"Tool name uses verb + object pattern (create_invoice, not invoice)","verified":"search-level"},
  {"id":"A16","level":"A","claim":"客户端应缓存工具定义","source":"https://modelcontextprotocol.io/docs/2025-11-25/develop/clients/client-best-practices","quote":"Cache tool definitions","verified":"search-level"},
  {"id":"A17","level":"A","claim":"OWASP: 工具响应视为不可信数据","source":"https://cheatsheetseries.owasp.org/cheatsheets/MCP_Security_Cheat_Sheet.html","quote":"Treat every tool response as untrusted data, including responses from approved servers.","verified":"search-level"},
  {"id":"A18","level":"A","claim":"server 禁止接受非本 server 签发的 token","source":"https://modelcontextprotocol.io/specification/draft/basic/security_best_practices","quote":"MCP servers MUST NOT accept any tokens that were not explicitly issued for the MCP server.","verified":"search-level"},
  {"id":"A19","level":"A","claim":"C ABI 仅 6 个 cpphub_v1_* 符号","source":"repo:src/c_api/c_api.cpp","quote":"cpphub_v1_get_abi_version / cpphub_v1_get_last_error / cpphub_v1_bsm_price_batch / cpphub_v1_heston_price / cpphub_v1_bsm_implied_vol / cpphub_v1_mc_price","verified":"E3-code-read"},
  {"id":"A20","level":"A","claim":"Python 绑定 58 个 m.def","source":"repo:python/src/cpphub_python.cpp","quote":"(Select-String -Path python\\src\\cpphub_python.cpp -Pattern '\\bm\\.def\\(' | Measure-Object).Count = 58","verified":"E1-command"},
  {"id":"A21","level":"A","claim":"CMake 开关 C_API/PYTHON/CUDA 默认 OFF","source":"repo:CMakeLists.txt","quote":"option(CPPHUB_ENABLE_C_API ...) OFF","verified":"E3-code-read"},
  {"id":"A22","level":"A","claim":"155 个 test_*.cpp; tests/validation 存在","source":"repo:tests/","quote":"(Get-ChildItem tests -Recurse -Filter test_*.cpp | Measure-Object).Count = 155","verified":"E1-command"},
  {"id":"B1","level":"B","claim":"社区存在 >=4 个独立 C++ MCP SDK","chain":["检索 MCP server C++ SDK","命中 tsar-mcp/cppmcp/cxxmcp/cpp-mcp/mcp_server(plugin)","逐一登记仓库名"],"deps":["A10"]},
  {"id":"B2","level":"B","claim":"社区存在多个量化/交易类 MCP server","chain":["检索 quantitative finance MCP server","命中 Quant Research Assistant/HPSILab/tradingview-mcp/MonteWalk/FactorHub"],"deps":[]},
  {"id":"B3","level":"B","claim":"原生 C++ MCP 单 exe 零运行时依赖可行","chain":["stdio 传输=子进程 stdin/stdout (A3)","TSAR-MCP 已以纯 stdio 零依赖实现 (A10)","结论: 本库可同路径产出单 exe"],"deps":["A3","A10"]},
  {"id":"B4","level":"B","claim":"量化类 MCP 工具面以分析/诊断为主","chain":["B2 案例工具清单观察","均为指标/回测/诊断型, 非原始数值库逐函数映射"],"deps":["B2"]},
  {"id":"B5","level":"B","claim":"无状态核心与本库纯函数诊断契合","chain":["MCP 无状态 (A1/A2)","本库诊断函数输入序列->输出 TestResult, 无会话态","结论: 无需服务端状态管理"],"deps":["A1","A2","A19"]},
  {"id":"C1","level":"C","claim":"形态建议: B1 Python 壳先行 -> B2 原生, 底座 JSON-through-ABI","rationale":"低成本验证工具设计; 社区 SW 证据显示 C++ 原生可行但工具面锁定成本高"},
  {"id":"C2","level":"C","claim":"工具面优先诊断类, 组合优于单点, 数量<=10","rationale":"社区 token 预算约束(A14) + 差异化资产定位"},
  {"id":"C3","level":"C","claim":"安全默认: 只读 + 响应标记为数据 + 描述校验 + 窄面","rationale":"A17/A18 硬约束 + 公开库暴露面控制"},
  {"id":"C4","level":"C","claim":"发布: 本地 stdio 优先, Registry 登记, 远程 HTTP 延后","rationale":"公开库无服务端; 避免 OAuth/托管成本"}
]
```

---

## 附录 C：假设区（无 URL+引文，禁止进入正文结论）

- **H1** 本仓库目标 MCP 客户端（Trae / Claude Desktop / Cursor）对 `protocolVersion: 2026-07-28` 的支持矩阵未知——未做实测。
- **H2** 官方 MCP Registry 的服务提交与审核流程细节未核验（A13 仅证明 Registry 存在）。
- **H3** 假设 agent 侧对"诊断类组合工具"的接受度高于"单点工具"——依据为 B4 的案例归纳，非实测调用日志。

---

**Review 签字**: _________ 日期: _________
