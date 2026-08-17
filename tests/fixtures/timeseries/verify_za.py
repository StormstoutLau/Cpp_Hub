# -*- coding: utf-8 -*-
# verify_za.py - statsmodels 0.14.6 Zivot-Andrews 对照基准生成 (Phase 7C M0)
#
# 用途: 生成固定合成数据 + statsmodels (Baum 近似) 三模型基准值,
#       供 test_zivot_andrews.cpp 硬编码断言 (spec §3.2.1, 容差 1e-10)
# 环境: conda open-webui (statsmodels 0.14.6): C:/Users/Peng/.conda/envs/open-webui/python.exe
#
# 已知对照边界 (2026-08-17 实测归因, 见 zivot_andrews_test.hpp 头注):
#   - "c"/"ct" (Model A/C): stat/p/lag/bpidx 全对 ~1e-14
#   - "t" (Model B): stat/p/lag 对, 但 bpidx 差 1 — statsmodels "t" 模型
#     DT 用 cutoff−1 边界 + trend[0:] (与 urca 形不同); B 主对照 = urca (verify_za.R)
#
# 运行: python verify_za.py  (打印基准 + 落盘 za_statsmodels_baselines.json
#       + 生成 tests/unit/timeseries/za_baseline.inc)
import json
import os

import numpy as np
from statsmodels.tsa.stattools import zivot_andrews

HERE = os.path.dirname(os.path.abspath(__file__))
INC_PATH = os.path.normpath(os.path.join(HERE, "..", "..", "unit", "timeseries",
                                         "za_baseline.inc"))

# R urca 1.3-4 ur.za(y, model=m, lag=1) 全精度转录 (verify_za.R 2026-08-17 输出,
# 12 位小数; C++ 实现与之逐位一致, 容差 1e-8)。bpoint 为 1-based (C++ break_index+1)。
# 更新流程: 改数据/模型后重跑 verify_za.R, 将输出转录至此 (12 位冗余 >> 1e-8 容差)。
URCA = {
    "intercept": (-3.724924811426, 68),
    "trend": (-2.795429231053, 53),
    "both": (-3.760689783846, 68),
}


def make_za_data() -> np.ndarray:
    """固定合成数据 (RandomState(42), T=120, 真断点 0-based 70, 崩溃均值 + 趋势)"""
    rng = np.random.RandomState(42)
    T = 120
    e = rng.standard_normal(T + 50)[50:]  # 与 C++/R 侧同一构造
    y = np.cumsum(e)
    y[70:] += 4.0
    y = y + 0.05 * np.arange(T)
    return y


def gen_inc(y: np.ndarray, out: dict) -> None:
    """生成 tests/unit/timeseries/za_baseline.inc (数据数组 + 双库基线)"""
    vals = [repr(float(v)) for v in y]
    rows = ["    " + ", ".join(vals[i:i + 4]) + ","
            for i in range(0, len(vals), 4)]

    def case(reg, name):
        c = out[reg]
        return (f"constexpr ZaCase {name}{{{repr(c['stat'])}, {repr(c['p'])}, "
                f"{c['bpidx']}, {c['lag']}}};  // regression=\"{reg}\"")

    ua, ub, uc = URCA["intercept"], URCA["trend"], URCA["both"]
    inc = "\n".join([
        "// 自动生成: tests/fixtures/timeseries/verify_za.py (statsmodels 0.14.6)",
        "//   数据: RandomState(42), T=120, 真断点 0-based 70 (崩溃均值 + 趋势)",
        "//   sm 基线: zivot_andrews(trim=0.15, autolag=\"AIC\") — Baum 近似模式",
        "//   urca 基线: verify_za.R 输出转录 (lag=1, 12 位小数, 2026-08-17 逐位一致)",
        "// 勿手改 — 重新生成请运行该脚本",
        "#pragma once",
        '#include "cpphub/core/types.hpp"',
        "namespace cpphub { inline namespace v1 { namespace timeseries {",
        "namespace unit_root { namespace za_baseline {",
        "",
        "constexpr Size T = 120;",
        "constexpr Size TRUE_BREAK = 70;      // 0-based 构造断点 (y[70:] += 4.0)",
        "constexpr Size GRID_015 = 84;        // trim=0.15 候选数 (b in [19,102])",
        "",
        "constexpr Real Y_ZA[T] = {",
        *rows,
        "};",
        "",
        "// statsmodels 基线 (bpidx 0-based, statsmodels 约定)",
        "struct ZaCase { Real stat, p; Size bpidx, lag; };",
        case("c", "SM_A"),
        "// Model B bpidx quirk: C++(urca 形) = SM_B.bpidx - 1 (DT cutoff 边界差异,",
        "//     stat/p 仍 ~1e-14 一致; 见 zivot_andrews_test.hpp 头注)",
        case("t", "SM_B"),
        case("ct", "SM_C"),
        "",
        "// urca 1.3-4 ur.za(lag=1, 无 trim) 基线 — 固定 lag 主模式对照 (容差 1e-8)",
        "// bpoint 为 1-based: C++ break_index = bpoint - 1",
        f"constexpr Real URCA_A_STAT = {repr(ua[0])};   // intercept",
        f"constexpr Real URCA_B_STAT = {repr(ub[0])};   // trend",
        f"constexpr Real URCA_C_STAT = {repr(uc[0])};   // both",
        f"constexpr Size URCA_A_BPOINT = {ua[1]};  // 1-based",
        f"constexpr Size URCA_B_BPOINT = {ub[1]};",
        f"constexpr Size URCA_C_BPOINT = {uc[1]};",
        "",
        "}}}}}  // namespace cpphub::v1::timeseries::unit_root::za_baseline",
        "",
    ])
    with open(INC_PATH, "w", encoding="utf-8", newline="\n") as f:
        f.write(inc)
    print("inc ->", INC_PATH)


def main() -> None:
    y = make_za_data()
    np.savetxt(os.path.join(HERE, "za_smoke_data.csv"), y, delimiter=",",
               header="y", comments="")
    out = {}
    for reg in ("c", "t", "ct"):
        stat, p, cvd, baselag, bpidx = zivot_andrews(
            y, trim=0.15, regression=reg, autolag="AIC")
        out[reg] = {"stat": float(stat), "p": float(p), "lag": int(baselag),
                    "bpidx": int(bpidx),
                    "cv": {k: float(v) for k, v in cvd.items()}}
    path = os.path.join(HERE, "za_statsmodels_baselines.json")
    json.dump(out, open(path, "w"), indent=1)
    print(json.dumps(out, indent=1))
    print("baselines ->", path)
    gen_inc(y, out)


if __name__ == "__main__":
    main()
