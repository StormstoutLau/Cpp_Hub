# 实施文档：Cpp_Hub MCP Server（mcp-server）

---
id: MCP-SERVER-IMPLEMENTATION
type: design
version: 1.0
status: draft
date: 2026-10-18
depends: [MCP-SERVER-DESIGN, MCP-SERVER-RESEARCH, MCP-SERVER-DECISION-MATRIX]
upstream: null
---

> **Feature**: mcp-server（Cpp_Hub 对外 MCP 接口层）
> **状态**: draft — Step 6 依赖核验已执行（本文件 §2.2/§9.1 记录核验结果与 1 处设计偏差）
> **Spec 步骤**: Step 5–6
> **基于设计**: [DESIGN.md](./DESIGN.md) v1.1 · **基于调研**: [RESEARCH.md](./RESEARCH.md)

---

## 1. 实施概述

分三层落地，**不改动 C++ 内核数值路径**：(1) **C ABI 扩展层**——新增 JSON-through-ABI 的 2 个动态符号与最小 JSON 读写；(2) **B1 协议壳**——用官方 MCP Python SDK 暴露 7 个诊断工具，把调用透传到 ABI；(3) **构建与测试接入**——CMake 开关 + 契约测试 + 工具 schema 快照测试。B2（C++ 原生 stdio）为终态，本轮不实施。

## 2. 工程细节

### 2.1 技术栈

| 组件 | 技术 | 版本 | 验证状态 |
|---|---|---|---|
| 协议 SDK | 官方 MCP Python SDK（PyPI `mcp`） | **>= 2.3.0** | ✅ 已查 PyPI（页面级，2026-10-18） |
| 运行语言 | Python | **>= 3.10**（PyPI 元数据）/ 3.11+（SDK 文档） | ⚠️ **双源冲突**，见 §2.2 |
| Schema/校验 | Pydantic（SDK 传递依赖） | 随 SDK | ✅ SDK 文档明确以 Pydantic 建模（页面级） |
| 内核 | C++20 header-only | 既有 | ✅ 不变 |
| C ABI | 既有 `include/cpphub/c_api/` | `CPPHUB_ABI_VERSION` | ✅ 读码 E3 |
| 构建 | CMake ≥ 3.25 | 既有 | ✅ |
| JSON（C++ 侧） | **自实现最小读写** | 无第三方 | ✅ 由 U-6 定案，见 §3.3 |
| 测试 | GoogleTest（C++）+ pytest（Python） | 既有 | ✅ |

### 2.2 依赖版本验证（Step 6 门禁）

| 依赖 | 声明 | 上游实测 | 兼容性 | 验证方式 |
|---|---|---|---|---|
| `mcp` | >= 2.3.0 | **2.3.0**（PyPI 页标题 "mcp 2.3.0"、"Release files for mcp 2.3.0"） | ✅ 存在且为当前版 | WebFetch PyPI（页面级） |
| Python 下限 | — | PyPI 元数据 `Requires: Python >=3.10`；SDK 文档页写 `Python 3.11+` | ⚠️ **冲突** | WebFetch 双源（页面级） |
| Pydantic | 传递依赖 | SDK 文档示例 `from pydantic import BaseModel, Field` | ✅ | WebFetch SDK 文档（页面级） |
| `ToolAnnotations` | 传递依赖 | SDK 文档 `from mcp.types import ToolAnnotations` | ✅ | 同上 |

> **冲突处置（未掩盖）**：Python 下限双源不一致。**实施取 `>=3.11`（更严者）** 作为 CI 目标下限，并在 `pyproject.toml` 声明 `requires-python = ">=3.10"` 与 PyPI 元数据一致；冲突登记为 **U-7**（四元块见 §10）。
>
> **无虚构库/函数**：上表所有标识均来自官方页面取证；本仓库不引入任何未验证第三方（`third_party/` 现仅 `autodiff`/`eigen`，E1 实测）。

### 2.3 文件结构（新增/改动）

```
cpphub/                                  # 仓库根
├── include/cpphub/c_api/
│   ├── cpphub_c_api.h                   # [改] +2 符号声明 + CPPHUB_ERR_BUFFER_TOO_SMALL(-4)
│   └── cpphub_tool_registry.hpp         # [新] 工具注册表 (name/schema/description/handler)
├── include/cpphub/c_api/detail/
│   └── json_min.hpp                     # [新] 最小 JSON 读写 (无第三方依赖)
├── src/c_api/
│   └── c_api.cpp                        # [改] 实现 cpphub_v1_tool_list / cpphub_v1_tool_call
├── python/src/
│   ├── cpphub_python.cpp                # [不变] 既有 nanobind 绑定
│   └── mcp_server.py                    # [新] B1 协议壳 (官方 SDK)
├── python/tests/
│   └── test_mcp_tools.py                # [新] 工具契约 + schema 快照测试
├── tests/c_api/
│   └── test_c_api_tools.cpp             # [新] JSON 往返 + 缓冲协商 + 错误码契约
└── spec/mcp-server/                     # 本 feature 文档 (四件套)
```

## 3. 模块实施

### 3.1 `c_api_ext`（C ABI 扩展层）

**职责**（DESIGN §3.2）：暴露 `tool_list` / `tool_call` 动态符号。

**接口签名**（与 DESIGN §4.2 一致）：

```c
int cpphub_v1_tool_list(char* out_json, size_t* out_len);
int cpphub_v1_tool_call(const char* tool_name, const char* args_json,
                        char* out_json, size_t* out_len);
```

**实施要点**：
- `out_len` 为**新增约定**（DESIGN §4.2 独立 pass F2 已订正）：入参 = 容量、出参 = 所需长度；`*out_len` 不足时返回新增错误码 `CPPHUB_ERR_BUFFER_TOO_SMALL` (-4) 并回填所需长度，供调用方二次分配。
- 复用既有线程局部错误字符串：失败时经 `cpphub_v1_get_last_error()` 取细节（读码 E3：头文件已有该符号）。
- 符号命名沿用既有 `cpphub_v1_` 前缀（ADR-009 版本化纪律）；ABI 版本号随符号变更递增。

**低效操作排除**：

| 潜在低效 | 排除措施 |
|---|---|
| 每次 `tool_list` 重建 schema 字符串 | 注册表在首次调用时惰性构建并缓存（只读，无状态语义不冲突：缓存是常量，非请求态） |
| 大数组 JSON 逐字符拼接导致 O(n²) | 预估容量 + 单次预留缓冲；数组走整数 → 十进制快速路径 |
| 同一参数多次解析 | `tool_call` 内 JSON 解析一次，传值而非传文本到 handler |

### 3.2 `json_serializer`（C++ 侧最小 JSON）

**职责**：`TestResult` ↔ JSON 互转。

**`U-6` 定案（Type-2，本轮确定）**：**C++ 侧自实现最小 JSON 读写**（`json_min.hpp`）。

- 论据（E1）：`third_party/` 实测仅 `autodiff`/`eigen`，**无任何 JSON 库**；header-only 库引入 nlohmann/json 会显著增加编译期成本与依赖面。
- 范围限定：仅需支持本 feature 所需子集（对象 / 数组 / 字符串 / 数字 / 布尔 / null），**不做通用 JSON 库**；`double` 往返须保证 **IEEE-754 位精度**（沿用最短往返表示）。
- Python 侧**不**自实现：schema 与校验交给 SDK/Pydantic（§3.4）。

**低效操作排除**：解析器为单遍、无回溯；字符串解码一次性分配。

### 3.3 `tool_registry`（工具注册表）

**职责**：集中登记 7 个工具的 `name` / `inputSchema` / `description` / handler。

**实施要点**：
- 单一登记点：新增工具只改此处 → 所有前端（B1/B2/未来语言）自动获得（DESIGN §5.1 收益）。
- `description` 内嵌文献陷阱注释（GR6 列序方向 / CI1 p 与 cv 分列），使 agent 消费自带免疫。
- **描述文本哈希钉**：注册时对 (name, description, inputSchema) 三元组取哈希写测试快照 → 防"工具投毒式"静默变更（DESIGN §10.2 R-3）。

### 3.4 `mcp_server`（B1 协议壳）

**职责**：协议实现、工具注册、转发。

**接口签名**（**已按 SDK v2 实测 API 订正**，见 §9.1 F3）：

```python
# python/src/mcp_server.py
from mcp.server import MCPServer              # SDK v2 入口 (2.3.0 实测)
from mcp.types import ToolAnnotations

mcp = MCPServer("cpphub")

@mcp.tool(
    name="cpphub_unit_root_diagnose",         # 覆盖式命名 (SDK 默认取函数名)
    description=CPPHUB_TOOL_DOCS["unit_root_diagnose"],   # 含陷阱注释
    annotations=ToolAnnotations(readOnlyHint=True),        # 映射 DESIGN §8 不变式 4 (只读)
)
def unit_root_diagnose(
    series: Annotated[list[float], Field(description="价格/收益序列")],
    spec: Literal["c", "ct", "ctt"] = "c",
) -> DiagnosticEnvelope:
    """单位根六检验一揽子诊断。"""
    return cpphub_abi.tool_call("unit_root_diagnose", {"series": series, "spec": spec})
```

**实施要点**：
- **返回类型化对象** → SDK 自动同时产出 `content`（模型读）与 `structured_content`（客户端读）——正是 DESIGN §6.1 信封的落地形态。
- 参数用 `Annotated[..., Field(...)]` / `Literal` / Pydantic 模型描述，schema 由 SDK 生成，无需手写 JSON Schema。
- 透传层为唯一业务入口：壳内**不**做数值逻辑（DESIGN §2.3 职责边界）。

**待核项（联调轮）**：`@mcp.tool()` 的 `name` / `description` / `annotations` 精确 kwargs 形式以实际安装版本为准（SDK 文档 §"Names, titles, and annotations" 已确认可覆盖，但未逐字取证 kwargs 名）。

### 3.5 `input_validator`

**职责**：schema 校验、长度/类型边界。

**实施要点**：C++ 侧做**硬边界**（数组长度 ≤ `kMaxSeriesLength`、类型/必填、非有限值拒绝）；Python 侧由 Pydantic 做**类型与取值域**校验。两侧都不信任输入。

### 3.6 `diagnostics`（既有，不改）

**职责**：诊断函数本体（DESIGN §2.3 职责外）。

**实施要点**：**零改动**；仅在注册表内以函数指针/`std::function` 适配为 handler。底层 7 工具对应头文件已逐一实测存在（`adf_test.hpp` / `pp_test.hpp` / `kpss_test.hpp` / `df_gls_test.hpp` / `ng_perron_test.hpp` / `zivot_andrews_test.hpp` / `variance_ratio_test.hpp` / `engle_granger.hpp` / `johansen_test.hpp` / `vecm_model.hpp` / `garch_*` / `granger_test.hpp` / `dy_spillover.hpp`）。

## 4. 接口实施

> **编号映射（Step 8 回写，F4 处置 (a)）**：本文章节号与 DESIGN 非 1:1（信封属设计的数据结构节），映射如下：

| IMPL 小节 | 承接自 DESIGN | 说明 |
|---|---|---|
| §4.1 `DiagnosticEnvelope` | §6.1 数据结构 | 出参信封在 DESIGN 归入数据结构节 |
| §4.2 C ABI 2 符号 | §4.2 接口定义 | 编号一致 |
| §4.3 `mcp` 工具声明 | §4.1 工具契约 + §4.3 Python 薄壳 | 工具名/入参 schema 本体在 §4.1；SDK 声明形态在 §4.3 |

### 4.1 `DiagnosticEnvelope`（统一出参）

```python
@dataclass
class DiagnosticEnvelope:
    tool: str
    method: str
    statistic: float
    p_value: float | None
    cv_5pct: float | None
    reject_null: bool
    notes: list[str]
```

**字段基准**：与 `TestResultBase`（`statistic` / `p_value` / `method_name` / `reject_null`，读码 E3 实测 4 字段）对齐；`cv_5pct` 为库内既有用法（5 处）。**签名一致性**：与 DESIGN §6.1 一致 ✅

### 4.2 C ABI 2 符号

**签名一致性**：与 DESIGN §4.2 一致 ✅（含 F1/F2 订正结果）

### 4.3 `mcp` 工具声明

**签名一致性**：与 DESIGN §4.1 的 7 工具名/入参一致 ✅；**SDK 调用形态与 DESIGN §4.3 示意不一致**——DESIGN 写的是 v1 `FastMCP`，实测 v2 为 `MCPServer`（§9.1 F3）。**处置：本文档为准，DESIGN §4.3 已同步订正。**

## 5. 兼容性

### 5.1 Python 版本兼容

| Python 版本 | 支持 | 说明 |
|---|---|---|
| 3.11+ | ✅ | CI 目标下限（取双源更严者） |
| 3.10 | ⚠️ | PyPI 元数据允许；未在 CI 覆盖，属 U-7 |

### 5.2 依赖兼容

| 依赖组合 | 兼容性 | 验证 |
|---|---|---|
| `mcp 2.3.0` + Pydantic | ✅ | SDK 文档示例即用 Pydantic |
| `mcp` + 既有 nanobind 绑定 | ⚠️ 待验证 | 二者共处同一 Python 环境，须在 Step 9 验证 import 无冲突 |
| C++20 (MSVC/GCC) + `json_min.hpp` | ✅ | 纯头文件、无第三方 |

### 5.3 向后兼容

- **C ABI**：仅**新增**符号，既有 6 符号签名与语义不变（ADR-009 版本化）。
- **协议**：dual-era（DESIGN §8 不变式 5）；现代与 legacy 客户端均可完成 `tools/list` + `tools/call`。
- **工具面**：出参字段只增不改（DESIGN §8 不变式 6）。

## 6. 错误处理实施

| 错误场景（DESIGN §7） | 异常/错误码 | 处理位置 | 测试位置 |
|---|---|---|---|
| 工具名未知 | `CPPHUB_ERR_INVALID_ARG` | `tool_call` 注册表查找 | `test_c_api_tools.cpp` |
| 参数类型/范围非法 | `CPPHUB_ERR_INVALID_ARG` / Pydantic `ValidationError` | `input_validator` | 两侧单测 |
| 缓冲不足 | `CPPHUB_ERR_BUFFER_TOO_SMALL` (-4) | `tool_list` / `tool_call` | 缓冲协商用例 |
| 序列超长（>10k） | `CPPHUB_ERR_INVALID_ARG` + 可操作提示 | `input_validator` | 边界用例 |
| 协议版本不匹配 | `UnsupportedProtocolVersionError` (`-32022`) | SDK 承载 | 联调（U-4） |
| 内核运行时异常 | `CPPHUB_ERR_RUNTIME` + `get_last_error()` | `tool_call` catch | 异常注入用例 |

## 7. 不变式实施

| 不变式（DESIGN §8） | 实施位置 | 验证方式 |
|---|---|---|
| 1 无状态 | 注册表缓存为常量；无请求态 | 连续调用同一工具两次，结果逐位相同 |
| 2 工具名唯一 + 字符集 | 注册表登记时断言 | 注册表快照测试（唯一性 + 正则） |
| 3 工具面 ≤10 | 注册表规模断言 | 单测断言 `registry.size() <= 10` |
| 4 只读 | `ToolAnnotations(readOnlyHint=True)` + 无副作用实现 | 注解快照 + 代码审查 |
| 5 版本协商可达 | dual-era 声明（SDK 承载） | 联调（U-4） |
| 6 出参 schema 稳定 | 信封字段只追加 | 快照测试锁字段集 |

## 8. 测试策略

### 8.1 单元测试

| 模块 | 测试文件 | 目标用例数 | 覆盖点 |
|---|---|---|---|
| `json_min.hpp` | `tests/c_api/test_json_min.cpp` | ≥15 | 往返位精度 / 转义 / 嵌套 / 边界 |
| `tool_registry` | `tests/c_api/test_tool_registry.cpp` | ≥10 | 唯一性 / 规模 ≤10 / schema 快照 |
| `c_api_ext` | `tests/c_api/test_c_api_tools.cpp` | ≥12 | 缓冲协商 / 错误码 / JSON 往返 |
| B1 壳 | `python/tests/test_mcp_tools.py` | ≥9 | 7 工具冒烟 + 信封字段 + 只读注解 |

### 8.2 集成测试

- **ABI 契约测试**（扩既有 `tests/c_api/`）：`tool_list` → 解析 → 逐工具 `tool_call` 闭环。
- **三平台一致性**：沿用既有 CI 矩阵（Win10 MSVC / Ubuntu GCC×2），JSON 输出须**逐字节一致**。

### 8.3 属性测试

- 不变式 1（无状态）与 6（schema 稳定）以快照 + 重复调用差分实现；不引入新测试框架。

## 9. 幻觉抑制审查（Step 6 Review）

> ⚠️ **RULE-1**：本节与正文**同一次生成**，勾选标 `[自查·并发]`，**不计独立复核**；独立 pass 待执行。

### 9.1 依赖版本验证

- [ ] `[自查·并发]` 所有依赖版本已查询 PyPI 确认存在（`mcp` = 2.3.0，页面级取证）
- [ ] `[自查·并发]` **无虚构的库或函数**——**发现 F3：DESIGN §4.3 的 `FastMCP` 为 SDK v1 形态，实测 v2 为 `MCPServer`**（已订正）
- [ ] `[自查·并发]` 版本兼容性已验证（Python 下限双源冲突 → 登记 U-7，取更严者）

### 9.2 接口签名验证

- [ ] `[自查·并发]` 所有接口签名与 DESIGN §4 一致（**例外：§4.3 SDK 形态 → 已按实测订正并回写 DESIGN**）
- [ ] `[自查·并发]` 所有签名可实现（`MCPServer` + `@mcp.tool()` + Pydantic 均为官方文档实证）

### 9.3 实施与设计对齐

- [ ] `[自查·并发]` 6 个模块逐一对应 DESIGN §3.2 六行
- [ ] `[自查·并发]` 每个接口对应 DESIGN §4 相应小节
- [ ] `[自查·并发]` 无设计未覆盖的实施（B2 明确列为终态、本轮不做）

### 9.4 低效操作排除

- [ ] `[自查·并发]` 无 O(n²) 可优化项（JSON 容量预估；解析单遍）
- [ ] `[自查·并发]` 无多余 I/O（无临时文件）
- [ ] `[自查·并发]` 无重复计算（schema 缓存 + 单次解析）

**Step 6 发现汇总**：1 × P2（**F3** SDK API 版本偏差，已订正设计）+ 1 × 待裁（Python 下限双源冲突 → U-7）。

## 10. 实施步骤与未决项

### 10.1 实施步骤

1. **S1** `json_min.hpp` + 单测 → 2. **S2** `cpphub_tool_registry.hpp` + 2 符号声明（含新增错误码）→ 3. **S3** `c_api.cpp` 实现 + ABI 契约测试 → 4. **S4** B1 壳 `mcp_server.py`（7 工具）→ 5. **S5** pytest + 三平台 JSON 逐字节一致验证 → 6. **S6** 与真实 MCP 客户端联调（U-4）。

### 10.2 未决项登记（四元块）

- **U-4 客户端实际支持矩阵**（承 DESIGN §11）— 判据/触发/倾向/证据等级同 DESIGN。
- **U-5 官方 Versioning 页不一致定性**（承 DESIGN §11）。
- **U-7 Python 下限双源冲突**
  - **判据**：CI 在 3.10 与 3.11 均跑通，或上游统一声明
  - **触发条件**：Step 9 首次 CI 运行（可机械识别：3.10 job 结果）
  - **当前倾向**：取 **>=3.11** 为 CI 目标；`pyproject.toml` 声明 `>=3.10` 以匹配 PyPI 元数据
  - **证据等级**：E1（页面级双源，冲突已登记）
- **U-6 JSON 序列化路线** — **本轮已闭合**：C++ 侧自实现最小 JSON；Python 侧交 SDK/Pydantic（见 §3.2/§3.3）。

---

**Review 签字**: _________ 日期: _________
