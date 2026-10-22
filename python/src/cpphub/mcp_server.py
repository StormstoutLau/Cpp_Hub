"""Cpp_Hub MCP server — B1 薄壳 (spec/mcp-server DESIGN §4.3)。

职责边界 (DESIGN §2.3)：本模块只做协议声明与转发。
- 工具语义 / inputSchema 硬校验 / 数值计算全部在 C++ 注册表 (`cpphub/c_api/`)；
- 本层不实现任何数值逻辑，也不复制工具描述（描述取自 `mcp_tool_list`）；
- 出参统一为 `ToolResult` 容器 (DESIGN §6.2)，由 `_core.mcp_tool_call` 原样返回 JSON。

运行 (stdio):
    python -m cpphub.mcp_server

依赖 (可选 extra):
    pip install "cpphub[mcp]"
"""

from __future__ import annotations

import json
from typing import Annotated, Any, Literal, Optional, TypedDict

from mcp.server import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import ToolAnnotations
from pydantic import Field

from . import _core

__all__ = [
    "DiagnosticEnvelope",
    "Conclusion",
    "ToolResult",
    "mcp",
    "main",
]


# ---------------------------------------------------------------------------
# 出参结构 (DESIGN §6.1 / §6.2)。TypedDict 让 SDK 直接以此生成 outputSchema
# 与 structured_content，不再包一层 {"result": ...}。
# ---------------------------------------------------------------------------
class DiagnosticEnvelope(TypedDict):
    """逐方法诊断信封 (对应 TestResultBase 的 statistic/p_value/method_name/reject_null)。"""

    tool: str
    method: str
    statistic: Optional[float]  # C++ 侧 NaN -> JSON null
    p_value: Optional[float]
    cv_5pct: Optional[float]
    reject_null: bool
    notes: list[str]  # 文献陷阱注释 (随工具下发)


class Conclusion(TypedDict):
    """综合结论：只取主检验行，避免跨方法 H0 方向不一致 (F6)。"""

    primary_method: str
    reject_null: bool
    n_total: int


class ToolResult(TypedDict):
    """所有工具的统一出参容器 (F5)。"""

    tool: str
    results: list[DiagnosticEnvelope]
    conclusion: Conclusion


_READ_ONLY = ToolAnnotations(readOnlyHint=True)

mcp = MCPServer("cpphub")

# 工具名/描述的唯一真值源在 C++ 注册表；import 时即暴露不一致 (快速失败)
_SPECS = {spec["name"]: spec for spec in json.loads(_core.mcp_tool_list())}


def _describe(name: str) -> str:
    return _SPECS[name]["description"]


def _call(tool: str, **kwargs: Any) -> ToolResult:
    """转发到 C ABI 统一入口；参数非法/内核失败转成 ToolError (模型可读)。"""
    args = {key: value for key, value in kwargs.items() if value is not None}
    try:
        payload = json.dumps(args, allow_nan=False, separators=(",", ":"))
    except ValueError as exc:
        raise ToolError(f"{tool}: non-finite number in arguments ({exc})") from exc
    try:
        raw = _core.mcp_tool_call(tool, payload)
    except (ValueError, RuntimeError) as exc:
        raise ToolError(f"{tool}: {exc}") from exc
    return json.loads(raw)


@mcp.tool(
    name="cpphub_unit_root_diagnose",
    description=_describe("cpphub_unit_root_diagnose"),
    annotations=_READ_ONLY,
)
def unit_root_diagnose(
    series: Annotated[list[float], Field(description="价格/收益水平序列")],
    trend: Literal["auto", "nc", "c", "ct"] = "auto",
    tests: Optional[
        list[Literal["adf", "pp", "kpss", "df_gls", "ng_perron", "za"]]
    ] = None,
    za_model: Literal["A", "B", "C"] = "C",
    max_lag: Annotated[int, Field(ge=0, description="0=自动")] = 0,
) -> ToolResult:
    return _call(
        "cpphub_unit_root_diagnose",
        series=series,
        trend=trend,
        tests=tests,
        za_model=za_model,
        max_lag=max_lag,
    )


@mcp.tool(
    name="cpphub_variance_ratio_test",
    description=_describe("cpphub_variance_ratio_test"),
    annotations=_READ_ONLY,
)
def variance_ratio_test(
    series: Annotated[list[float], Field(description="价格/收益水平序列")],
    k: Annotated[
        Optional[list[int]], Field(description="horizon 列表, 逐项 >= 2; 缺省 [2,5,10,20]")
    ] = None,
    use_debiased: bool = True,
    trend: Literal["c", "n"] = "c",
) -> ToolResult:
    return _call(
        "cpphub_variance_ratio_test",
        series=series,
        k=k,
        use_debiased=use_debiased,
        trend=trend,
    )


@mcp.tool(
    name="cpphub_cointegration_scan",
    description=_describe("cpphub_cointegration_scan"),
    annotations=_READ_ONLY,
)
def cointegration_scan(
    series: Annotated[
        list[list[float]], Field(description="K 列等长序列 (K>=2), 外层为变量列")
    ],
    eg_trend: Literal["n", "c", "ct", "ctt"] = "c",
    det_order: Literal[-1, 0, 1] = 0,
    k_ar_diff: Annotated[int, Field(ge=0)] = 0,
) -> ToolResult:
    return _call(
        "cpphub_cointegration_scan",
        series=series,
        eg_trend=eg_trend,
        det_order=det_order,
        k_ar_diff=k_ar_diff,
    )


@mcp.tool(
    name="cpphub_vecm_fit",
    description=_describe("cpphub_vecm_fit"),
    annotations=_READ_ONLY,
)
def vecm_fit(
    series: Annotated[
        list[list[float]], Field(description="K 列等长水平序列 (K>=2)")
    ],
    rank: Annotated[int, Field(ge=1, description="协整秩 r, 须在 [1,K]")],
    k_ar_diff: Annotated[int, Field(ge=0)] = 0,
    det: Literal["n", "co", "ci", "lo", "li"] = "n",
    urca_normalization: bool = False,
) -> ToolResult:
    return _call(
        "cpphub_vecm_fit",
        series=series,
        rank=rank,
        k_ar_diff=k_ar_diff,
        det=det,
        urca_normalization=urca_normalization,
    )


@mcp.tool(
    name="cpphub_garch_diagnose",
    description=_describe("cpphub_garch_diagnose"),
    annotations=_READ_ONLY,
)
def garch_diagnose(
    series: Annotated[list[float], Field(description="收益序列")],
    model: Literal["GARCH", "EGARCH", "GJR", "GARCH-M"],
    dist: Literal["normal", "t", "ged"] = "normal",
    form: Annotated[
        Literal["var", "vol", "log"], Field(description="仅 GARCH-M 有效")
    ] = "vol",
) -> ToolResult:
    return _call(
        "cpphub_garch_diagnose",
        series=series,
        model=model,
        dist=dist,
        form=form,
    )


@mcp.tool(
    name="cpphub_granger_causality",
    description=_describe("cpphub_granger_causality"),
    annotations=_READ_ONLY,
)
def granger_causality(
    cause: Annotated[
        list[float], Field(description="原因序列 (进入 effect 方程的滞后项)")
    ],
    effect: Annotated[list[float], Field(description="效应序列 (被解释变量)")],
    lag: Annotated[int, Field(ge=1)] = 1,
    d_max: Annotated[int, Field(ge=0, description="TY 增广阶; 0=不做 TY")] = 0,
    with_hac: bool = True,
    hac_bandwidth: Annotated[int, Field(ge=0)] = 0,
) -> ToolResult:
    return _call(
        "cpphub_granger_causality",
        cause=cause,
        effect=effect,
        lag=lag,
        d_max=d_max,
        with_hac=with_hac,
        hac_bandwidth=hac_bandwidth,
    )


@mcp.tool(
    name="cpphub_dy_spillover",
    description=_describe("cpphub_dy_spillover"),
    annotations=_READ_ONLY,
)
def dy_spillover(
    series: Annotated[
        list[list[float]], Field(description="K 列等长序列 (K>=2)")
    ],
    horizon: Annotated[int, Field(ge=1)] = 10,
    lag: Annotated[int, Field(ge=0, description="0=每窗口 IC 自动选阶")] = 0,
    trend: Literal["c", "n"] = "c",
    rolling: Annotated[
        int, Field(ge=0, description="滚动窗口大小; 0=全样本静态; 须 > 2*K")
    ] = 0,
) -> ToolResult:
    return _call(
        "cpphub_dy_spillover",
        series=series,
        horizon=horizon,
        lag=lag,
        trend=trend,
        rolling=rolling,
    )


def main() -> None:
    """以 stdio 传输启动 MCP 服务器 (本地子进程, 无网络开销)。"""
    mcp.run()


if __name__ == "__main__":
    main()
