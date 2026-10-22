# 设计文档：Cpp_Hub MCP Server（mcp-server）

---
id: MCP-SERVER-DESIGN
type: design
version: 1.2
status: draft
date: 2026-10-18
depends: [MCP-SERVER-RESEARCH, MCP-SERVER-DECISION-MATRIX, ADR-009, ADR-010, ADR-015]
upstream: null
---

> **Feature**: mcp-server（Cpp_Hub 对外 MCP 接口层）
> **状态**: draft — Step 4 独立 pass ✅（v1.1）；Step 6 触发 §4.3 双向对账订正（v1.2）
> **Spec 步骤**: Step 3–4
> **基于调研**: [RESEARCH.md](./RESEARCH.md) · 裁定依据: [DECISION_MATRIX.md](./DECISION_MATRIX.md)

---

## 1. 设计目标

让 **Agent 编程工具**（Trae / Claude Desktop / Cursor 等 MCP 客户端）能以**最窄的工具面**调用 Cpp_Hub 的**诊断能力**，并获得**机器可判的结构化结论**（统计量 / p 值 / 临界值 / `reject_null`）与**文献陷阱注释**——即把本库的差异化资产（`TestResult` + D1–D23 排幻觉标注）交付给 agent，而非把库 API 逐函数平铺。

## 2. 设计依据

### 2.1 调研结论 → 设计决策

| 调研发现（证据锚点） | 设计决策 | 引用 |
|---|---|---|
| MCP 为**无状态**协议；服务器 MUST NOT 依赖前序请求（A1/A2，页面级） | 工具实现为**纯函数**；不引入服务端会话态 | RESEARCH §3.1.1 |
| 无握手原语为 legacy（`2025-11-25` 及更早）；`Legacy client × Modern server = Fails`（E-1/E-4，页面级） | **U-1 → dual-era**：目标 `2026-07-28`，兼容 legacy | DECISION_MATRIX §2/§5 |
| 官方 SDK **无 C++**；Tier 1 = TS/Python/C#/Go（E-7，页面级） | **U-2 → B1（Python Tier 1 SDK 薄壳）先行**，B2 为终态 | DECISION_MATRIX §3/§5 |
| 「包装既有 API」为官方点名反模式；工具过多分散 agent（E-10/E-11，页面级） | **U-3 → 仅诊断类组合工具（≤10）**，不暴露定价类 | DECISION_MATRIX §4/§5 |
| stdio = 本地子进程、无网络开销（E-8，页面级） | 首选 **stdio 传输**；Streamable HTTP 延后 | RESEARCH §3.1.2 / C4 |
| 官方 reference servers 均为窄域工作流（E-12，页面级） | 本 server 限定「计量诊断」单一窄域 | DECISION_MATRIX §7 |
| C ABI 现仅 6 个标量符号（A19，E3 读码） | 新增 **JSON-through-ABI 骨干**（2 个动态符号） | RESEARCH §3.6 / §6.1 |

### 2.2 相关 ADR

| ADR | 决策摘要 | 对本设计的影响 |
|---|---|---|
| ADR-009（C ABI 稳定边界 + 版本化） | C ABI 为跨语言稳定边界，符号版本化（`cpphub_v1_*`） | 扩展符号沿用 `cpphub_v1_` 前缀与既有错误码机制 |
| ADR-010（nanobind Python 绑定） | Python 经 nanobind 暴露 | B1 薄壳复用现有绑定；热路径直绑、冷路径 JSON 透传 |
| ADR-015（证伪统计量模块边界：通用 vs 模块特定） | 诊断结果结构统一 | 工具输出信封以 `TestResultBase` 字段为基准 |

### 2.3 职责边界

- **职责内**：MCP 协议壳、工具注册表与 schema、JSON (反)序列化层、输入校验与错误映射、版本协商声明。
- **职责外**：数值算法本身（属 C++ 内核）、定价数学、任何研究 IP 内容、远端托管/鉴权服务。

## 3. 架构设计

### 3.1 整体架构

```
Agent / MCP Client
      │  stdio (JSON-RPC 2.0, 逐请求携带 _meta.protocolVersion)
      ▼
┌───────────────────────────────────────────────┐
│ MCP Server 壳  (B1: 官方 Python mcp SDK)       │  ← 无业务逻辑, 仅协议与转发
│  · tools/list  ──► 工具注册表 (schema + 描述)   │
│  · tools/call  ──► 统一入口 tool_call           │
└───────────────────────────────────────────────┘
      │  JSON in / JSON out
      ▼
┌───────────────────────────────────────────────┐
│ JSON-through-C-ABI 骨干                        │
│  cpphub_v1_tool_list() / cpphub_v1_tool_call() │  ← 新增 2 动态符号, 与现有 6 标量符号并存
│  + JSON 序列化层 (TestResult → JSON)            │
└───────────────────────────────────────────────┘
      │
      ▼
┌───────────────────────────────────────────────┐
│ C++ 内核 (header-only, 不动)                    │
│  include/cpphub/econometrics/**  →  TestResult │
└───────────────────────────────────────────────┘
```

**B2 终态**：仅将最上层「MCP Server 壳」由 Python 换为 C++ 原生 stdio 单 exe（参照零依赖社区 SDK 形态，E-9）；骨干与其下完全复用。

### 3.2 模块划分

| 模块 | 职责 | 输入 | 输出 | 依赖 |
|---|---|---|---|---|
| `mcp_server`（壳） | 协议实现、工具注册、转发 | JSON-RPC 请求 | JSON-RPC 响应 | 官方 MCP SDK |
| `tool_registry` | 工具名/schema/描述（含陷阱注释） | 注册调用 | `tools/list` 负载 | 无 |
| `json_serializer` | `TestResult` ↔ JSON 互转 | C++ 结构 / JSON 文本 | 反之 | 无（见 U-6） |
| `c_api_ext` | `tool_list` / `tool_call` 动态符号 | JSON 字符串 | JSON 字符串 + 错误码 | `json_serializer` |
| `diagnostics` | 诊断函数（已在库内） | 序列 / 规格 | `TestResult` | C++ 内核（ADR-015） |
| `input_validator` | schema 校验、长度/类型边界 | 原始参数 | 校验结果或错误 | 无 |

### 3.3 数据流

`agent 请求 → 壳解析 → 参数校验（失败即错误响应）→ 序列化为 JSON → cpphub_v1_tool_call → C++ 诊断 → TestResult → JSON → 壳封装 → agent`

### 3.4 控制流

1. 客户端启动 server（stdio 子进程）→ 2. `server/discover` 或 `tools/list` → 3. 客户端调用 `tools/call` → 4. 壳校验并转发 → 5. ABI 调度到诊断函数 → 6. 结果回填 → 7. 版本不匹配时改走 `UnsupportedProtocolVersionError` 协商路径（E-3）。

## 4. 接口定义

### 4.1 MCP 工具契约（`tools/list` 的 `inputSchema`）

> 命名遵循 `cpphub_` 命名空间前缀 + `verb_object`（E-10/E-11 页面级指引）；字符集 `[A-Za-z0-9_-.]`、≤128 字符（A7 页面级）。

| 工具名 | 入参（JSON Schema 摘要） | 出参要点 |
|---|---|---|
| `cpphub_unit_root_diagnose` | `{series: number[], trend?: "auto"\|"nc"\|"c"\|"ct" (默认 auto), tests?: ("adf"\|"pp"\|"kpss"\|"df_gls"\|"ng_perron"\|"za")[] (默认全部), za_model?: "A"\|"B"\|"C" (默认 C), max_lag?: integer (0=自动)}` | ADF/PP/KPSS/DF-GLS/NP/ZA 一揽子 + 综合结论 |
| `cpphub_variance_ratio_test` | `{series: number[], k?: integer[] (默认 [2,5,10,20]), use_debiased?: bool (默认 true), trend?: "c"\|"n" (默认 c)}` | Lo-MacKinlay Z1/Z2 + CLM debiased + Chow-Denning SMM 联合 |
| `cpphub_cointegration_scan` | `{series: number[][], eg_trend?: "n"\|"c"\|"ct"\|"ctt" (默认 c), det_order?: -1\|0\|1 (默认 0), k_ar_diff?: integer (默认 0)}` | EG（K=2 时）+ Johansen 迹/λmax + rank 判定 |
| `cpphub_vecm_fit` | `{series: number[][], rank: integer, k_ar_diff?: integer (默认 0), det?: "n"\|"co"\|"ci"\|"lo"\|"li" (默认 n), urca_normalization?: bool (默认 false)}` | α/β/ECT t + EM2002 临界值 + Π=αβ′ |
| `cpphub_garch_diagnose` | `{series: number[], model: "GARCH"\|"EGARCH"\|"GJR"\|"GARCH-M", dist?: "normal"\|"t"\|"ged" (默认 normal), form?: "var"\|"vol"\|"log" (仅 GARCH-M, 默认 vol)}` | 参数 + sandwich SE + LB(z)/LB(z²)/JB 残差诊断 + 稳定性 |
| `cpphub_granger_causality` | `{cause: number[], effect: number[], lag?: integer (默认 1), d_max?: integer (默认 0=不做 TY), with_hac?: bool (默认 true), hac_bandwidth?: integer (默认 0)}` | F/χ²/LR + TY 增广 Wald + HAC-Wald |
| `cpphub_dy_spillover` | `{series: number[][], horizon?: integer (默认 10), lag?: integer (默认 0=IC 自动), trend?: "c"\|"n" (默认 c), rolling?: integer (滚动窗口, 须 > 2K; 缺省=全样本静态)}` | TCI/TO/FROM/NET 矩阵 + 滚动 tci_path/net_path |

> **⚠️ F7 订正（2026-10-18，Step 9 实施期 API 核验触发，P3）**：上表**入参列已于 Step 9 按库内实际签名逐项重写**——Step 3 初版入参为**未核验的拟写值**，与实现期实测不符，逐项对照如下（IMPLEMENTATION.md §4.3 曾声明"与 DESIGN §4.1 的 7 工具名/入参一致 ✅"，该声明在**入参**意义上未经 API 核验，一并更正）：
>
> | 工具 | 初版拟写（未核验） | 实测不符处 |
> |---|---|---|
> | unit_root | `spec:{model:"c"\|"ct"\|"ctt"}` | 单位根族 trend 词汇为 `auto/nc/c/ct`（KPSS/DF-GLS/NP 仅 `c/ct`），**无 `ctt`**；`model` 是 ZA 专有量（A/B/C） |
> | variance_ratio | `q`, `bootstrap?` | 库内为 `k_list`（多 horizon）+ `use_debiased`/`trend`；**该 API 无 wild bootstrap 参数**，出参要点中的 "wild bootstrap" 同步更正 |
> | cointegration | `spec:{trend, coint_trend}` | 库内为 `eg_trend`（EG 两步法）+ `det_order`/`k_ar_diff`（Johansen） |
> | vecm | `lag`, `determ` | 库内为 `rank`（必填）+ `k_ar_diff` + `det`（`n/co/ci/lo/li` 五情形） |
> | garch | `model` | 需另加 `dist`（normal/t/ged）与 GARCH-M 的 `form`（var/vol/log） |
> | granger | `lag`, `d_max` | 需另加 `with_hac`/`hac_bandwidth`；参数顺序 **显式 `(cause, effect)`**（GR6 方向陷阱） |
> | dy_spillover | `lag`, `horizon`, `rolling?` | 需另加 `trend`；`rolling` 即 `window`，须 `> 2·K`（否则库内抛错） |
>
> **证据等级**：E3（静态读码，绑 commit `4f6c74f` 的 `include/cpphub/timeseries/**` 逐签名）；**处置**：上表已改为实测词汇，§4.3 与 IMPL §4.3 同步更正。

> **出参承载（F5 订正，2026-10-18）**：本表"出参要点"列所列的多统计量内容**全部**由 §6.2 `ToolResult.results[]` 承载（每项即 §6.1 信封），顶层另带 `conclusion`；**不采用**"变体扁平信封"（裁定依据见 [DECISION_MATRIX.md](./DECISION_MATRIX.md) §11–§12）。

**工具数上限 ≤10**（含后续增量），description 中内嵌文献陷阱注释（如 GR6 列序方向 / CI1 p 与 cv 分列），使 agent 消费时自带免疫。

### 4.2 C ABI 扩展符号（`include/cpphub/c_api/`）

```c
/* 枚举注册工具: name + JSON schema + description */
int cpphub_v1_tool_list(char* out_json, size_t* out_len);
/* 统一调用入口: 参数与结果均 JSON 字符串 */
int cpphub_v1_tool_call(const char* tool_name,
                        const char* args_json,
                        char* out_json, size_t* out_len);
```

- 缓冲语义为 JSON 扩展**新增约定**（独立复核 E1 探针证实：现有 6 符号均为调用方分配的标量/批量 double 数组接口，无 buffer-query 先例）：`out_len` 入参 = 缓冲容量、出参 = 实际长度；不足时返回 `CPPHUB_ERR_BUFFER_TOO_SMALL` 并回填所需长度。
- 错误码：**复用**既有 `CPPHUB_OK` / `CPPHUB_ERR_INVALID_ARG` / `CPPHUB_ERR_RUNTIME` / `CPPHUB_ERR_UNKNOWN`，并**新增** `CPPHUB_ERR_BUFFER_TOO_SMALL` (-4)——独立复核证实该码现有头文件中不存在，须随本 feature 定义；细节经 `cpphub_v1_get_last_error()`。

### 4.3 Python 薄壳（B1）

```python
# python/src/cpphub/mcp_server.py  (SDK v2 实测 API; 见 §4.3 注)
from mcp.server import MCPServer              # SDK 2.3.0 入口 (页面级取证)
from mcp.types import ToolAnnotations

mcp = MCPServer("cpphub")

@mcp.tool(name="cpphub_unit_root_diagnose",
          description=CPPHUB_TOOL_DOCS["unit_root_diagnose"],
          annotations=ToolAnnotations(readOnlyHint=True))
def unit_root_diagnose(series: list[float], spec: dict | None = None, tests: list[str] | None = None) -> dict:
    """执行单位根六检验一揽子诊断, 返回结构化 TestResult 信封。"""
    return cpphub_abi.tool_call("unit_root_diagnose", {"series": series, "spec": spec, "tests": tests})
```

> **§4.3 订正记录（2026-10-18，Step 6 核验触发）**：本小节初版写作 `from mcp.server.fastmcp import FastMCP` / `FastMCP("cpphub")`，属**官方 SDK v1 形态**。Step 6 依 PyPI（`mcp` = 2.3.0）与官方 SDK 文档（页面级）核验，实测 v2 入口为 `from mcp.server import MCPServer`，工具以 `@mcp.tool()` 声明、默认取函数名并可经装饰器覆盖，返回类型化对象自动产出 `structured_content`。已按实测订正；**精确 kwargs 形式仍待联调轮逐字核验**（不视为已验证契约）。

## 5. 替代方案

### 5.1 方案 A：JSON-through-C-ABI 骨干 + B1 Python 壳（**选定**）

- 描述：一次定义「工具注册 + 统一 JSON 入口」于 C ABI，协议壳退化为转发层。
- 优点：新功能只改 C++ 侧注册，所有前端自动获得；ABI 稳定性由既有版本化机制保证；B1 可用 Tier 1 SDK 快速验证工具设计。
- 缺点：需先扩 C ABI 与序列化层（前置成本）。
- **选择理由**：唯一同时满足「最窄壳 / 多语言复用 / 无官方 C++ SDK 约束」的形态（E-7/E-9）。

### 5.2 方案 B：直接上 C++ 原生 stdio（跳过 B1）（**否决**）

- 优点：一步到位、零 Python 分发。
- 否决理由：官方无 C++ SDK（E-7 页面级），须依赖社区实现，维护承诺弱于 Tier 1；工具 schema 尚未冻结即锁死实现路径，返工风险高（U-2 判为 **Type-2**，应先做再改）。

### 5.3 方案 C：逐函数 nanobind 直绑为 MCP 工具（**否决**）

- 优点：无 ABI 改动、复用现有绑定。
- 否决理由：**落入官方点名反模式**（E-10：工具仅包装既有 API）；且把 58 个 `m.def` 平铺为工具面会突破 ≤10 上限、放大注入面（E-11）。

## 6. 数据结构

### 6.1 诊断结果信封（`DiagnosticEnvelope`，**逐方法**）

```python
@dataclass
class DiagnosticEnvelope:
    tool: str            # 工具名
    method: str          # 检验方法名 (对应 TestResultBase.method_name)
    statistic: float
    p_value: float | None
    cv_5pct: float | None
    reject_null: bool
    notes: list[str]     # 文献陷阱注释 (随工具下发, 如 GR6 列序方向)
```

### 6.2 工具出参容器（`ToolResult`，**所有工具的统一出参**）

> **⚠️ F5 订正（2026-10-18，Step 9 实施期发现；裁定见 [DECISION_MATRIX.md](./DECISION_MATRIX.md) §9–§12）**：§6.1 初版的小节标题曾写作"所有工具的统一出参"，但 §4.1 的 7 个工具**全部**为多统计量一揽子（7/7 复核，复算命令见 MATRIX §9），单方法信封无法承载。经七维矩阵裁定为**方案 (a) 容器层**：§6.1 信封本体**原样保留**，另加容器承载 `results[]`。

```python
@dataclass
class ToolResult:
    tool: str                            # 工具名 (与 §4.1 一致)
    results: list[DiagnosticEnvelope]    # 逐方法信封, 长度 >= 1 (单统计量工具长度为 1)
    conclusion: Conclusion               # 综合结论 (形状固定, 避免引入变体 schema)

@dataclass
class Conclusion:
    primary_method: str   # 主检验方法名 (= results 中某项的 method)
    reject_null: bool     # 恒等于主检验行的 reject_null (H0 方向由主检验定义)
    n_total: int          # results 的方法数 (纯计数, 无方向语义)
```

> **⚠️ F6 细化（2026-10-18，P3，**实现前**冻结）**：初版 `Conclusion` 曾设想 `n_reject`（"一揽子中 `reject_null == True` 的方法数"），但各方法的 **H0 方向不一致**——KPSS 的 H0 是**平稳性**，与 ADF/PP/DF-GLS/NP/ZA 的 H0（单位根）**相反**（unit_root U12 陷阱，读码 E3）。跨方法直接计数会产出**错误的**综合结论，这正是 §11 裁定中否掉方案 (b) 的"语义漂移"同类风险。故改为**方向显式**：`reject_null` **只取主检验行**，`n_total` 为纯计数；跨方法的同/异向一致性由 agent 读 `results[]` 各行 `notes` 自行判断。

**与不变式 6 的关系**：容器层为**纯新增字段**（§6.1 的 7 个字段一字未改），属"只增不改"；`conclusion` 形状**先固定**，否则会重新引入变体 schema（MATRIX §12 回写清单第 1 项）。

## 7. 错误处理

| 错误场景 | 处理方式 | 用户可见信息 |
|---|---|---|
| 工具名未知 | ABI 返回 `CPPHUB_ERR_INVALID_ARG` | JSON-RPC 错误 + 可用工具名列表 |
| 参数类型/范围非法 | `input_validator` 拦截，不进入内核 | JSON-RPC `-32602` + 字段级说明 |
| 序列超长（>10k） | 拒绝并提示分页/截断策略 | 明确的可操作提示（K3 token 预算） |
| 协议版本不匹配 | 回 `UnsupportedProtocolVersionError`（`-32022`），列 `supported` | 客户端据此重试（E-3） |
| 内核运行时异常 | `CPPHUB_ERR_RUNTIME` + `get_last_error()` | 摘要信息，不含内部路径 |

## 8. 不变式（Invariants）

1. **无状态**：任一工具调用仅依赖本次请求内容，不依赖前序请求（A1/A2）。
2. **工具名唯一**且符合 `[A-Za-z0-9_-.]`、长度 ≤128（A7）。
3. **工具面 ≤10**；新增工具须同步更新注册表与文档，禁止"临时直通"。
4. **只读**：所有工具无副作用（无写文件、无网络外发）。
5. **版本协商可达**：modern 与 legacy 客户端均可完成 `tools/list` + `tools/call`。
6. **出参 schema 稳定**：新增字段只追加，不移除/不改名既有字段（工具面为 Type-1 契约）。

## 9. 幻觉抑制审查（Step 4 Review）

> ⚠️ **RULE-1 标注**：§9.1–9.3 原始勾选与正文同一次生成（`[自查·并发]`）；**2026-10-18 已由独立 pass 复核**（异基座会话，不采信自查声明、全部 E1 探针重放），通过项升级为 `[已复核·独立pass]`，发现项经修正后升级——全过程见 §9.4。

### 9.1 设计基于已验证的调研结论

- [x] `[已复核·独立pass]` 所有设计决策可追溯到 RESEARCH / DECISION_MATRIX（§2.1 映射表 7 行锚点全部解析命中：RESEARCH §3.1.1/§3.1.2/§3.6/K3/A19/A21，MATRIX E-1/E-3/E-7/E-8/E-9/E-10/E-11/E-12）
- [x] `[已复核·独立pass]` 无未经验证的假设（FastMCP 签名已标 ⚠️ 待核；遗留待核项已入 §11 U-4/U-5/U-6 四元块）
- [x] `[已复核·独立pass]` 无论证驱动的归因扭曲（已失效声明 H-2「25k tokens」未被引用；**复核发现 2 处仓库事实归因失实 F1/F2，已于 v1.1 修正**）

### 9.2 替代方案审查

- [x] `[已复核·独立pass]` 列出 2 个替代方案（§5.2 C++ 原生 / §5.3 逐函数直绑）
- [x] `[已复核·独立pass]` 每个替代方案有明确否决理由（均绑定页面级锚点 E-7/E-9/E-10/E-11）

### 9.3 职责边界审查

- [x] `[已复核·独立pass]` 职责边界清晰（§2.3：协议壳/注册表/序列化/校验 vs 数值算法；§10.1 「不改动 C++ 内核数值路径」可执行）
- [x] `[已复核·独立pass]` 不越界吞并其他研究范式（不承载研究 IP、不做远端鉴权服务）

### 9.4 独立 pass 复核记录（2026-10-18）

> 复核者：异基座独立会话（与起草会话不同模型）。输入 = 本文档 v1.0 + 仓库代码 + E1 探针重放；**不采信**起草侧自查声明。

**E1 探针与结果**（第三方可原样重放）：

| # | 探针（命令） | 目标声明 | 结果 |
|---|---|---|---|
| P1 | `Get-Content include\cpphub\c_api\cpphub_c_api.h \| Select-String 'CPPHUB_ERR\|CPPHUB_OK'` | §4.2 错误码「复用既有」 | ❌ **F1**：`CPPHUB_ERR_BUFFER_TOO_SMALL` 不存在（实际仅 OK/INVALID_ARG/RUNTIME/UNKNOWN 4 码）→ 已修正为「新增 (-4)」 |
| P2 | 探针1 全符号清单 `cpphub_v1_*`（6 个）+ 签名审查 | §4.2「缓冲语义沿用现有约定」 | ❌ **F2**：现有 6 符号均为调用方分配 double 数组接口，无 out_json/size_t* 先例 → 已修正为「新增约定」 |
| P3 | 7 工具底层头文件逐一探测（adf_test / pp_test / kpss_test / df_gls_test / ng_perron_test / zivot_andrews_test / variance_ratio_test / johansen_test / engle_granger / vecm_model / garch_* / granger_test / dy_spillover） | §4.1 七工具可实现性 | ✅ 全部存在（含 `pp_test.hpp` + `test_pp_test.cpp`；首轮探针因 `-First 3` 截断产生 PP 缺失**假阳性**，精确复探排除——截断探针教训落档） |
| P4 | `test_result_base.hpp` L25-29 读码 | §6.1/§10.1 以 `TestResultBase` 为基准 | ✅ 4 字段 (statistic/p_value/method_name/reject_null) 与信封映射一致；`cv_5pct` 库内 5 处真实使用 |
| P5 | `git grep -n -E "GR6\|CI1" -- "*.md"` | §4.1 陷阱注释引用 | ✅ GR6（显式 cause/effect 方向）、CI1（MacKinnon 响应面 1e-12）均在受控文档落码 |
| P6 | §2.1 引用锚点逐一在 RESEARCH/MATRIX grep | 可追溯性 | ✅ 全部解析命中 |
| P7 | §2.1「C ABI 现仅 6 个标量符号」（A19 重放） | 符号计数 | ✅ 6 个精确 |

**发现汇总**：2 × P2（F1/F2，文档级失实，已修正并回写本节）/ 0 × P1 / 其余通过。
**修正生效**：v1.0 → v1.1（§4.2 两处措辞 + 本节新增）。
**遗留**：复核者与起草者同为本仓库 agent 体系（异基座但同人配置），完全异源复审（RULE-4 升级路径）待 U-4 端到端联调轮一并执行；Review 签字栏留待 Scott 确认。

## 10. 对实施的输入

### 10.1 关键工程约束

- C ABI 为**可选构建**（默认 OFF）：须 `-DCPPHUB_ENABLE_C_API=ON`；B1 壳另需 `-DCPPHUB_ENABLE_PYTHON=ON`（A21）。
- **不改动 C++ 内核数值路径**；本 feature 只新增 ABI/序列化/壳三层。
- 工具输出以 `TestResultBase` 字段为基准（ADR-015），不得另造并行结果结构。
- 工具描述文本视为**面向模型的输入**，须同时接受安全扫描（工具投毒防护）。

### 10.2 风险与缓解

| 风险 | 缓解措施 |
|---|---|
| R-1 规范漂移（MCP 迭代快） | 一律经官方 SDK 承载协议层；版本清单集中声明，便于随规范扩展 |
| R-2 工具面膨胀 | 硬上限 ≤10 + 代码评审门禁；组合工具优先 |
| R-3 注入面 | 只读 + 参数 schema 强校验 + 工具描述哈希钉（工具投毒防护） |
| R-4 双轨漂移（B1/B2 两套实现） | 业务逻辑只在 ABI 层一处；前端仅协议壳 |
| R-5 依赖版本虚报 | Step 6 强制核验 SDK 实际安装版本与 API 签名（§4.3 ⚠️ 义务） |

## 11. 未决项登记（四元块）与决策矩阵

**11.1 未决项登记**：

- **U-4 客户端实际支持矩阵**
  - **判据**：≥1 个真实客户端完成 `tools/list` + `tools/call` 闭环，且协商路径被日志证实
  - **触发条件**：首次端到端联调（可机械识别：出现 `UnsupportedProtocolVersionError` 或 `initialize`）
  - **当前倾向**：以官方 SDK 默认 dual-era 行为承载，先按 legacy 兼容设计
  - **证据等级**：待核

- **U-5 官方文档 Versioning 页与 `/specification/latest/` 不一致的定性**
  - **判据**：官方站点两者一致，或 SEP/官方说明界定 "current" 与 "latest" 语义
  - **触发条件**：下次规范迭代（版本号变更）时复查
  - **当前倾向**：以 `/latest/` 解析 + spec 内 Modern/Legacy 定义为实施依据
  - **证据等级**：E1（页面级，冲突已于 DECISION_MATRIX §6 H-1 登记）

- **U-6 JSON 序列化实现路线（Type-2，不启动矩阵）**
  - **判据**：序列化层可独立单测（往返一致 + 数值位精度保持）
  - **触发条件**：Step 5 编写 IMPLEMENTATION 时确定（可机械识别：是否新增第三方依赖）
  - **当前倾向**：**自实现最小 JSON 读写**（避免新增第三方依赖，与"零依赖"倾向一致）；若工作量超预期再评估引入成熟库
  - **证据等级**：C（判断）

**11.2 决策矩阵**：U-1 / U-2 / U-3 三项 Type-1 决策已按 [DECISION_MATRIX.md](./DECISION_MATRIX.md) 全维度取证并收敛（§5 收敛出口）；U-6 为 Type-2，**不启动矩阵**（先做再改）。

---

**Review 签字**: _________ 日期: _________
