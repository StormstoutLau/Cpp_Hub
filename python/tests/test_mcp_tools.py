"""MCP 壳契约测试 (spec/mcp-server S5)。

四层：
1. 注册表 ↔ SDK 声明：工具名集合、入参属性集合（防两处 schema 漂移）；
2. 进程内驱动 MCPServer：`tools/list` 与 `tools/call`（无需真实客户端，U-4 留给 S6）；
3. 出参结构：`structured_content` 必须是 DESIGN §6.2 的 ToolResult 容器；
4. 跨平台一致性：结构化 byte-exact + 数值容差（golden 快照）。

一致性口径（裁定 2026-10-20）：
  * **结构 byte-exact**：键序 / 键集 / 类型 / 字符串 / 布尔 / 整数 / null 逐字节一致；
    浮点数在结构中替换为 `<num>` 占位，**notes 内的浮点字面量**同样替换（整数保留原值）。
  * **数值容差**：JSON 数值按遍历序逐项比较，默认 `rtol=1e-9, atol=1e-12`；
    迭代估计类工具（GARCH QMLE）放宽到 `rtol=1e-3`，见 `TOLERANCES`。
  * 输入序列由**整数 LCG** 生成（不调用 libm），保证各平台输入逐位相同。

重新生成 golden（仅在有意变更契约时执行）：
    CPPHUB_MCP_REGEN_GOLDEN=1 pytest python/tests/test_mcp_tools.py -k write_golden

运行: pytest python/tests/test_mcp_tools.py -v
"""

import asyncio
import json
import math
import os
import re
from pathlib import Path

import pytest

mcp = pytest.importorskip("mcp", reason="需要 pip install 'cpphub[mcp]'")
pytest.importorskip("mcp.server", reason="需要 mcp >= 2.3.0 (MCPServer)")

from cpphub import _core  # noqa: E402
from cpphub import mcp_server  # noqa: E402
from mcp.server.mcpserver.exceptions import ToolError  # noqa: E402

TOOL_NAMES = [
    "cpphub_unit_root_diagnose",
    "cpphub_variance_ratio_test",
    "cpphub_cointegration_scan",
    "cpphub_vecm_fit",
    "cpphub_garch_diagnose",
    "cpphub_granger_causality",
    "cpphub_dy_spillover",
]

GOLDEN_PATH = Path(__file__).parent / "fixtures" / "mcp_golden.json"
REGEN = os.environ.get("CPPHUB_MCP_REGEN_GOLDEN", "") == "1"

DEFAULT_TOL = (1e-9, 1e-12)
TOLERANCES = {
    # QMLE 迭代终止条件对舍入路径敏感, 参数层面放宽; 结构仍 byte-exact
    "cpphub_garch_diagnose": (1e-3, 1e-9),
}

# notes 内嵌的浮点字面量 (整数保留, 以便 n_lags=16 之类仍逐字节比对)
_FLOAT_IN_TEXT = re.compile(r"-?\d+\.\d+(?:[eE][+-]?\d+)?")


# ---------------------------------------------------------------------------
# 确定性输入 (整数 LCG, 不依赖 libm → 跨平台逐位一致)
# ---------------------------------------------------------------------------
def _uniforms(seed: int, n: int) -> list[float]:
    state = seed
    out = []
    for _ in range(n):
        state = (1103515245 * state + 12345) % 2147483648
        out.append(state / 2147483648.0)  # 精确二进制小数
    return out


def _noise(seed: int, n: int) -> list[float]:
    """Irwin-Hall(12) 近似 N(0,1)：12 个均匀数之和减 6，全为精确运算。"""
    u = _uniforms(seed, n * 12)
    return [sum(u[i * 12:(i + 1) * 12]) - 6.0 for i in range(n)]


def _walk(seed: int, n: int, step: float = 0.5, start: float = 100.0) -> list[float]:
    out = []
    level = start
    for z in _noise(seed, n):
        level += step * z
        out.append(level)
    return out


def _trend_series(n: int = 300) -> list[float]:
    return _walk(2026, n)


def _coint_cols(n: int = 300) -> list[list[float]]:
    x = _walk(11, n)
    e = _noise(12, n)
    return [x, [0.8 * v + 3.0 + 0.3 * w for v, w in zip(x, e)]]


def _three_cols(n: int = 300) -> list[list[float]]:
    a = _walk(21, n)
    b = [v + 0.6 * w for v, w in zip(a, _noise(22, n))]
    c = [0.5 * v + 1.5 * w for v, w in zip(a, _noise(23, n))]
    return [a, b, c]


def _returns(n: int = 400) -> list[float]:
    """确定性 GARCH(1,1) 路径：sigma2_t = 1e-5 + 0.08 r^2 + 0.90 sigma2。"""
    z = _noise(31, n)
    var = 1e-4
    out = []
    for v in z:
        var = 1e-5 + 0.08 * (out[-1] ** 2 if out else 0.0) + 0.90 * var
        out.append(math.sqrt(var) * v)
    return out


CASES = {
    "cpphub_unit_root_diagnose": {
        "series": _trend_series(),
        "tests": ["adf", "pp", "kpss", "df_gls", "ng_perron", "za"],
    },
    "cpphub_variance_ratio_test": {"series": _trend_series(), "k": [2, 5, 10, 20]},
    "cpphub_cointegration_scan": {"series": _coint_cols(), "k_ar_diff": 1},
    "cpphub_vecm_fit": {"series": _coint_cols(), "rank": 1, "k_ar_diff": 1},
    "cpphub_garch_diagnose": {"series": _returns(), "model": "GARCH"},
    "cpphub_granger_causality": {
        "cause": _trend_series(),
        "effect": [v + w for v, w in zip(_trend_series(), _noise(41, 300))],
        "lag": 2,
    },
    "cpphub_dy_spillover": {"series": _three_cols(), "lag": 1},
}


def _list_tools():
    return asyncio.run(mcp_server.mcp.list_tools())


def _call(name: str, arguments: dict):
    return asyncio.run(mcp_server.mcp.call_tool(name, arguments))


# ---------------------------------------------------------------------------
# 一致性比较器
# ---------------------------------------------------------------------------
def _normalize(value):
    """结构骨架：浮点 → <num>；notes 内浮点字面量 → <num>；整数/文本/布尔/null 原样。"""
    if isinstance(value, dict):
        return {key: _normalize(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_normalize(item) for item in value]
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, float):
        return "<num>"
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        return _FLOAT_IN_TEXT.sub("<num>", value)
    raise TypeError(f"unexpected JSON value type: {type(value)!r}")


def _skeleton(payload) -> str:
    return json.dumps(_normalize(payload), ensure_ascii=False, separators=(",", ":"))


def _numbers(value, out: list[float]) -> None:
    if isinstance(value, dict):
        for item in value.values():
            _numbers(item, out)
    elif isinstance(value, list):
        for item in value:
            _numbers(item, out)
    elif isinstance(value, float):
        out.append(value)


def _snapshot(payload) -> dict:
    nums: list[float] = []
    _numbers(payload, nums)
    return {"skeleton": _skeleton(payload), "numbers": nums}


def _load_golden() -> dict:
    if not GOLDEN_PATH.exists():
        pytest.skip(f"缺少 golden 快照: {GOLDEN_PATH}")
    return json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# 1. 注册表 ↔ SDK 声明
# ---------------------------------------------------------------------------
def test_registry_has_seven_tools():
    specs = json.loads(_core.mcp_tool_list())
    assert [s["name"] for s in specs] == TOOL_NAMES
    assert len(specs) <= 10  # DESIGN §8 不变式 3


def test_sdk_exposes_same_tool_names():
    assert [t.name for t in _list_tools()] == TOOL_NAMES


@pytest.mark.parametrize("name", TOOL_NAMES)
def test_input_schema_properties_match_registry(name: str):
    specs = {s["name"]: s for s in json.loads(_core.mcp_tool_list())}
    cpp_schema = specs[name]["inputSchema"]
    sdk_schema = next(t for t in _list_tools() if t.name == name).input_schema
    assert set(sdk_schema["properties"]) == set(cpp_schema["properties"])
    assert set(sdk_schema.get("required", [])) == set(cpp_schema.get("required", []))


@pytest.mark.parametrize("name", TOOL_NAMES)
def test_description_comes_from_registry(name: str):
    specs = {s["name"]: s for s in json.loads(_core.mcp_tool_list())}
    tool = next(t for t in _list_tools() if t.name == name)
    assert tool.description == specs[name]["description"]


def test_read_only_annotation():
    for tool in _list_tools():
        assert tool.annotations is not None
        assert tool.annotations.read_only_hint is True


# ---------------------------------------------------------------------------
# 2. tools/call 进程内驱动
# ---------------------------------------------------------------------------
def test_call_variance_ratio_returns_container():
    result = _call("cpphub_variance_ratio_test", {"series": _trend_series()})
    assert not result.is_error
    payload = result.structured_content
    assert payload["tool"] == "cpphub_variance_ratio_test"
    assert len(payload["results"]) >= 1
    assert payload["conclusion"]["primary_method"] in {
        row["method"] for row in payload["results"]
    }
    assert payload["conclusion"]["n_total"] == len(payload["results"])
    # content 同时存在 (agent 读文本, 客户端读结构)
    assert result.content


def test_call_unit_root_returns_container():
    result = _call(
        "cpphub_unit_root_diagnose",
        {"series": _trend_series(), "tests": ["adf", "kpss"]},
    )
    assert not result.is_error
    payload = result.structured_content
    methods = {row["method"] for row in payload["results"]}
    assert {"adf", "kpss"} <= methods
    # F6: KPSS 的 H0 方向与 ADF 相反, 综合结论只取主检验行
    kpss = next(r for r in payload["results"] if r["method"] == "kpss")
    assert any("H0" in note for note in kpss["notes"])


def test_envelope_fields_are_nullable_numbers():
    payload = _call("cpphub_variance_ratio_test", {"series": _trend_series()}).structured_content
    for row in payload["results"]:
        assert set(row) == {
            "tool",
            "method",
            "statistic",
            "p_value",
            "cv_5pct",
            "reject_null",
            "notes",
        }
        for key in ("statistic", "p_value", "cv_5pct"):
            assert row[key] is None or isinstance(row[key], float)
        assert isinstance(row["notes"], list)


# ---------------------------------------------------------------------------
# 3. 错误路径 (统一经 ToolError 暴露可读消息)
# ---------------------------------------------------------------------------
def test_unknown_tool_is_tool_error():
    with pytest.raises(ToolError) as excinfo:
        _call("cpphub_not_a_tool", {})
    assert "cpphub_not_a_tool" in str(excinfo.value)


def test_invalid_argument_message_surfaces():
    with pytest.raises(ToolError) as excinfo:
        _call("cpphub_variance_ratio_test", {"series": _trend_series(), "k": [1]})
    assert "k" in str(excinfo.value)


def test_non_finite_series_rejected():
    series = _trend_series()
    series[10] = float("nan")
    with pytest.raises(ToolError) as excinfo:
        _call("cpphub_variance_ratio_test", {"series": series})
    assert "finite" in str(excinfo.value).lower()


# ---------------------------------------------------------------------------
# 4. 一致性
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", TOOL_NAMES)
def test_shell_is_pure_passthrough(name: str):
    """壳不得改写 ABI 出参: structured_content 与 _core 原始 JSON 完全相等。"""
    args = CASES[name]
    via_shell = _call(name, args).structured_content
    via_abi = json.loads(
        _core.mcp_tool_call(name, json.dumps(args, separators=(",", ":")))
    )
    assert via_shell == via_abi


@pytest.mark.parametrize("name", TOOL_NAMES)
def test_repeated_calls_are_identical(name: str):
    """DESIGN §8 不变式 1: 无状态 → 同一输入两次调用逐位相同。"""
    first = _snapshot(_call(name, CASES[name]).structured_content)
    second = _snapshot(_call(name, CASES[name]).structured_content)
    assert first == second


@pytest.mark.skipif(not REGEN, reason="仅在 CPPHUB_MCP_REGEN_GOLDEN=1 时写快照")
def test_write_golden():
    tools = {
        name: _snapshot(_call(name, args).structured_content)
        for name, args in CASES.items()
    }
    golden = {
        "meta": {
            "note": "spec/mcp-server 跨平台一致性快照; 结构 byte-exact + 数值容差",
            "tolerance": {"default": list(DEFAULT_TOL), **{
                k: list(v) for k, v in TOLERANCES.items()
            }},
            "rules": "浮点-><num>; notes 内浮点字面量-><num>; 整数/文本/布尔/null 原样",
            "known_issues": [],
        },
        "tools": tools,
    }
    GOLDEN_PATH.parent.mkdir(parents=True, exist_ok=True)
    GOLDEN_PATH.write_text(
        json.dumps(golden, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


@pytest.mark.skipif(REGEN, reason="正在重新生成 golden")
@pytest.mark.parametrize("name", TOOL_NAMES)
def test_golden_structure_is_byte_exact(name: str):
    expected = _load_golden()["tools"][name]
    got = _snapshot(_call(name, CASES[name]).structured_content)
    assert got["skeleton"] == expected["skeleton"]


@pytest.mark.skipif(REGEN, reason="正在重新生成 golden")
@pytest.mark.parametrize("name", TOOL_NAMES)
def test_golden_numbers_within_tolerance(name: str):
    expected = _load_golden()["tools"][name]["numbers"]
    got = _snapshot(_call(name, CASES[name]).structured_content)["numbers"]
    assert len(got) == len(expected), f"{name}: 数值个数不一致"
    rtol, atol = TOLERANCES.get(name, DEFAULT_TOL)
    for index, (value, ref) in enumerate(zip(got, expected)):
        assert math.isclose(value, ref, rel_tol=rtol, abs_tol=atol), (
            f"{name}: 第 {index} 个数值 {value!r} 偏离 golden {ref!r} "
            f"(rtol={rtol}, atol={atol})"
        )
