# -*- coding: utf-8 -*-
# verify_np_stata.py - Stata dfgls 逐 k MAIC/σ̂² 对照基准生成 (Phase 7C M0 → C-3 v1.1)
#
# 状态: Stata 18 MP 已装机 (C:\Program Files\Stata18\StataMP-64.exe, 2026-08-20 探针
#   验证: 手册例 Min MAIC = -6.136692@lag1 逐字复现, r(results) 列名/行序一手确认) —
#   本脚本从"占位"升级为实际管线, C6 预批降级条件解除
#
# 双场景:
#   A (smoke):     np_smoke_data.csv (RandomState(42) T=200 随机游走+轻趋势),
#                  maxlag(14) notrend → np_stata_baselines.csv (逐 k 14 行)
#   B (manual):    lutkepohl2 ln_inv (T=92, Stata 手册 Example 1 数据),
#                  maxlag(11) trend → np_stata_manual_example.csv (逐 k 11 行)
#                  + lutkepohl2_ln_inv.csv (原始序列, 供 C++ 侧读入)
#
# 关键口径 (v1.1 审计一手取证, 禁用默认):
#   - r(results) 列名 = k MAIC SIC RMSE DFGLS (dfgls.ado v1.1.1); 无 sd 列
#   - 行序降序 k=kmax..1, 行名全为 r1 (无意义) — 对照必须按 k 列值匹配
#   - RMSE = σ̂ (非 σ̂²); 平方后对照 C++ sigma2_k
#   - Stata Schwert 默认: 手册 floor[12((T+1)/100)^.25] 与 ado int(12(T/100)^.25)
#     文档-源码矛盾 (T=200 均为 14, C++ ceil(12(T/100)^.25)=15) — 两侧显式 maxlag 对齐
#
# 用法:
#   1. python verify_np_stata.py --emit-dofile   # 生成 np_dfgls.do + np_smoke_data.csv
#   2. "C:\Program Files\Stata18\StataMP-64.exe" /e do np_dfgls.do   # 批处理
#   3. python verify_np_stata.py --parse          # 解析 → 对照表
import argparse
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

STATA_EXE = r"C:\Program Files\Stata18\StataMP-64.exe"

# 显式 maxlag: Stata 手册/ado 两公式 T=200 均得 14; C++ schwert_lag(200)=ceil=15
# 两侧禁用默认, 一律 14 (C++ 侧 max_lag=14)
SMOKE_KMAX = 14
# 手册例: T=92, Stata 默认 (两公式均 11) — 显式 11
MANUAL_KMAX = 11


def make_np_data() -> np.ndarray:
    """固定合成数据 (RandomState(42), T=200, 随机游走 + 轻趋势)

    与 C++ test_ng_perron.cpp gen_rw 不同源 (Philox vs RandomState) — Stata 对照
    只锚 MAIC/σ̂² 逐 k 序列 (同一 CSV 数据进 Stata 与 C++ 即可, DGP 无须同构)
    """
    rng = np.random.RandomState(42)
    T = 200
    e = rng.standard_normal(T)
    y = np.cumsum(e) + 0.02 * np.arange(T)
    return y


DOFILE = """\
* np_dfgls.do — Stata dfgls 逐 k MAIC/RMSE 基准生成 (Stata 18 MP, 2026-08-20)
* r(results) 列名 k MAIC SIC RMSE DFGLS; 行序降序 k=kmax..1 (按 k 列值匹配); RMSE=σ̂
version 18.0
set more off

* --- 场景 A: smoke T=200, notrend, maxlag(14) 显式 ---
import delimited using "np_smoke_data.csv", clear case(preserve) stringcols(_all)
destring y, replace force
* dfgls 需 tsset (无时间列 → 用行号; 手册例 lutkepohl2 自带 quarterly tsset)
generate t = _n
tsset t
dfgls y, maxlag({kmax_smoke}) notrend
matrix R = r(results)
clear
svmat double R, names(col)
export delimited using "np_stata_baselines.csv", replace

* --- 场景 B: 手册例 lutkepohl2 ln_inv, trend, maxlag(11) 显式 ---
use "https://www.stata-press.com/data/r18/lutkepohl2.dta", clear
* 保留时间变量 qtr (tsset 依赖; dfgls 需 tsset)
keep ln_inv qtr
* lutkepohl2 的 ln_inv 为 float 存储 (~8 位) — recast double (无损) 后全精度
* 导出, 保证 C++ 读入值与 Stata 内部计算值逐位一致 (float→double 精确提升)
recast double ln_inv
format ln_inv %21.0g
export delimited using "lutkepohl2_ln_inv.csv", replace
dfgls ln_inv, maxlag({kmax_manual})
matrix R2 = r(results)
clear
svmat double R2, names(col)
export delimited using "np_stata_manual_example.csv", replace
exit, clear
"""


def emit_inc() -> None:
    """从管线产物生成 np_stata_baseline.inc (数值字符串原样透传, 无损)"""
    import csv as _csv

    def read_col(path, col):
        with open(path, newline="") as f:
            return [row[col] for row in _csv.DictReader(f)]

    smoke_y = read_col(os.path.join(HERE, "np_smoke_data.csv"), "y")
    ln_inv = read_col(os.path.join(HERE, "lutkepohl2_ln_inv.csv"), "ln_inv")
    smoke = sorted(
        _csv.DictReader(open(os.path.join(HERE, "np_stata_baselines.csv"))),
        key=lambda r: int(float(r["k"])))
    manual = sorted(
        _csv.DictReader(open(os.path.join(HERE, "np_stata_manual_example.csv"))),
        key=lambda r: int(float(r["k"])))
    assert len(smoke_y) == 200 and len(ln_inv) == 92
    assert len(smoke) == SMOKE_KMAX and len(manual) == MANUAL_KMAX
    assert [int(float(r["k"])) for r in smoke] == list(range(1, SMOKE_KMAX + 1))
    assert [int(float(r["k"])) for r in manual] == list(range(1, MANUAL_KMAX + 1))

    def fixnum(s):
        # Stata 导出前导小数点风格 → C++ 字面量 (-.xxx → -0.xxx / .xxx → 0.xxx)
        s = s.strip()
        if s.startswith("-."):
            return "-0." + s[2:]
        if s.startswith("."):
            return "0." + s[1:]
        return s

    def arr(vals, per=3, indent="    "):
        vals = [fixnum(v) for v in vals]
        lines = []
        for i in range(0, len(vals), per):
            lines.append(indent + ", ".join(vals[i:i + per]) + ",")
        return "\n".join(lines)

    inc = f"""\
// np_stata_baseline.inc - NP 模块 Stata 18 dfgls 对照基准
// (自动生成: python verify_np_stata.py --emit-inc; 数值字符串原样透传, 无损)
//
// 数据: np_smoke_data.csv (RandomState(42) T=200 随机游走+轻趋势) +
//       lutkepohl2 ln_inv (T=92, Stata 手册 Example 1; float 存储 →
//       recast double 全精度导出, 与 Stata 内部计算值逐位一致)
// 基准: Stata 18.0 MP dfgls 装机实测 (2026-08-20, r(results) 全精度导出):
//   - smoke: maxlag({SMOKE_KMAX}) notrend → 逐 k=1..{SMOKE_KMAX} MAIC/RMSE
//   - manual: maxlag({MANUAL_KMAX}) trend → 逐 k=1..{MANUAL_KMAX} MAIC/RMSE
//     (手册例 Min MAIC = -6.136692 at lag 1 实机复现)
// 口径: r(results) 行序降序 (本文件按 k 升序重排); RMSE = σ̂, σ̂² = RMSE²;
//   对照 = ng_perron_test(同数据, "c"/"ct", 同 maxlag) — σ̂² 分母修正
//   (ng_perron_test.hpp, SSR/n_r) 后实测差: MAIC ≤ 8e-14, RMSE ≤ 3e-16
#pragma once
namespace cpphub {{ inline namespace v1 {{ namespace timeseries {{
namespace np_stata_baseline {{
inline namespace v1 {{

inline constexpr Size T_SMOKE = 200;
inline constexpr double SMOKE_Y[] = {{
{arr(smoke_y)}
}};

inline constexpr Size T_MANUAL = 92;
inline constexpr double LN_INV[] = {{
{arr(ln_inv)}
}};

inline constexpr Size KMAX_SMOKE = {SMOKE_KMAX};
inline constexpr double SMOKE_MAIC[] = {{  // k=1..{SMOKE_KMAX} 升序
{arr([r["MAIC"] for r in smoke])}
}};
inline constexpr double SMOKE_RMSE[] = {{  // σ̂, k=1..{SMOKE_KMAX} 升序
{arr([r["RMSE"] for r in smoke])}
}};

inline constexpr Size KMAX_MANUAL = {MANUAL_KMAX};
inline constexpr double MANUAL_MAIC[] = {{  // k=1..{MANUAL_KMAX} 升序
{arr([r["MAIC"] for r in manual])}
}};
inline constexpr double MANUAL_RMSE[] = {{  // σ̂, k=1..{MANUAL_KMAX} 升序
{arr([r["RMSE"] for r in manual])}
}};

}}  // inline namespace v1
}}  // namespace np_stata_baseline
}}  // namespace timeseries
}}  // inline namespace v1
}}  // namespace cpphub
"""
    out = os.path.normpath(os.path.join(
        HERE, "..", "..", "unit", "timeseries", "np_stata_baseline.inc"))
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write(inc)
    print("inc ->", out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--emit-dofile", action="store_true")
    ap.add_argument("--parse", action="store_true")
    ap.add_argument("--emit-inc", action="store_true")
    args = ap.parse_args()

    if args.emit_dofile:
        y = make_np_data()
        csv = os.path.join(HERE, "np_smoke_data.csv")
        np.savetxt(csv, y, delimiter=",", header="y", comments="")
        do = os.path.join(HERE, "np_dfgls.do")
        with open(do, "w", encoding="utf-8", newline="\n") as f:
            f.write(
                DOFILE.format(kmax_smoke=SMOKE_KMAX, kmax_manual=MANUAL_KMAX)
            )
        print("do-file ->", do)
        print("数据    ->", csv)
        print("下一步  ->", f'"{STATA_EXE}" /e do np_dfgls.do (cwd=fixtures/timeseries)')
    elif args.emit_inc:
        emit_inc()
    elif args.parse:
        import csv as _csv

        for name, path in [
            ("smoke (T=200, notrend, kmax=14)",
             os.path.join(HERE, "np_stata_baselines.csv")),
            ("manual (lutkepohl2 ln_inv, trend, kmax=11)",
             os.path.join(HERE, "np_stata_manual_example.csv")),
        ]:
            rows = list(_csv.DictReader(open(path)))
            rows.sort(key=lambda r: int(float(r["k"])))  # 按 k 升序
            print(f"\n=== {name}: {len(rows)} rows ===")
            print("k, MAIC, RMSE (σ̂) — C++ 对照: maic[k] / sqrt(sigma2_k[k])")
            for r in rows:
                print(f"{int(float(r['k']))}, {float(r['MAIC']):.12f}, "
                      f"{float(r['RMSE']):.12f}")
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
