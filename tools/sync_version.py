#!/usr/bin/env python3
"""RISK-032: 版本单一源同步。

- 无参数: 读仓库根 VERSION → 写 include/cpphub/core/version.hpp
- --check: 只校验, 不一致退出码 1（CI 用, 防生成物与真值漂移）

单一真值 = VERSION。version.hpp 为生成物且随源码树提交
（源码树 header-only 用法 clone 后即可编译, 无生成步骤可先行）。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION_FILE = ROOT / "VERSION"
HEADER = ROOT / "include" / "cpphub" / "core" / "version.hpp"

TEMPLATE = """\
// 本文件由 tools/sync_version.py 从仓库根 VERSION 生成, 请勿手改。
// 单一真值: VERSION（RISK-032 两轴模型: 本文件 = C++ 库 API 版本;
// Python 发行版本见 python/pyproject.toml）。
#pragma once
#define CPPHUB_VERSION_MAJOR {major}
#define CPPHUB_VERSION_MINOR {minor}
#define CPPHUB_VERSION_PATCH {patch}
#define CPPHUB_VERSION "{version}"
"""


def read_version() -> str:
    text = VERSION_FILE.read_text(encoding="utf-8").strip()
    if not re.fullmatch(r"\d+\.\d+\.\d+", text):
        sys.exit(f"VERSION 文件格式非法: {text!r} (期望 MAJOR.MINOR.PATCH)")
    return text


def expected_content(version: str) -> str:
    major, minor, patch = version.split(".")
    return TEMPLATE.format(major=major, minor=minor, patch=patch, version=version)


def main() -> int:
    version = read_version()
    content = expected_content(version)
    if "--check" in sys.argv[1:]:
        if HEADER.exists() and HEADER.read_text(encoding="utf-8") == content:
            print(f"version.hpp 与 VERSION 一致: {version}")
            return 0
        sys.exit(f"version.hpp 与 VERSION ({version}) 不一致 — 运行: python tools/sync_version.py")
    # 用 Path.open 而非 write_text(newline=...): 后者需 Python 3.10,
    # 而本项目 requires-python >= 3.9 (Step 6 勘误 R-3)。
    with HEADER.open("w", encoding="utf-8", newline="\n") as f:
        f.write(content)
    print(f"已写入 {HEADER.relative_to(ROOT)} ({version})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
