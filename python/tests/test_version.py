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