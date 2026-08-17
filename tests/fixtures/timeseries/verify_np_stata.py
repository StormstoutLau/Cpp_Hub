# -*- coding: utf-8 -*-
# verify_np_stata.py - Stata dfgls 逐 k MAIC/σ̂² 对照基准生成 (Phase 7C M0)
#
# 状态: 占位 (readiness C6: Stata 可编程访问未验证, 降级路径预批 2026-08-17)
#   - 当前 C++ 侧基准 = 原文公式钉死 + 恒等式自检 (1e-12) + np_tables 精确 +
#     模拟方向断言 (test_ng_perron.cpp 18 用例)
#   - Stata 批处理可用后: 运行本脚本生成 np_stata_baselines.csv, 补
#     "逐 k MAIC/σ̂² vs Stata dfgls r(results) 1e-10" 硬编码断言 (spec §1.3)
#
# 用法 (Stata 装机后):
#   1. python verify_np_stata.py --emit-dofile   # 生成 np_dfgls.do
#   2. stata-mp -b do np_dfgls.do                 # 批处理 (或 StataMP-Noinst 版)
#   3. python verify_np_stata.py --parse          # 解析日志 → CSV 基准
#
# 数据: 与 C++ test_ng_perron.cpp 同构造 — Philox 替代为 numpy RandomState
#       (跨语言固定); 生成后以 %.17g 落盘 np_smoke_data.csv 供 Stata import
import argparse
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def make_np_data() -> np.ndarray:
    """固定合成数据 (RandomState(42), T=200, 随机游走 + 轻趋势)

    与 C++ 侧 gen_rw 不同源 (Philox vs RandomState) — Stata 对照只锚
    MAIC/σ̂² 逐 k 序列 (同一 CSV 数据进 Stata 与 C++ 即可, DGP 无须同构)
    """
    rng = np.random.RandomState(42)
    T = 200
    e = rng.standard_normal(T)
    y = np.cumsum(e) + 0.02 * np.arange(T)
    return y


DOFILE = """\
* np_dfgls.do - 逐 k MAIC/σ̂² 导出 (Stata 15+; dfgls 需 Stata 12+)
import delimited using np_smoke_data.csv, clear case(preserve) stringcols(_all)
destring y, replace force
* dfgls y, maxlag(K) notrend: r(results) 逐 k 行 [k, maic, sd] (sd = σ̂²(k))
dfgls y, maxlag({kmax}) notrend
matrix R = r(results)
clear
svmat double R, names(col)
rename k k_
rename maic maic_
rename sd sigma2_k_
export delimited using np_stata_baselines.csv, replace
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--emit-dofile", action="store_true")
    ap.add_argument("--parse", action="store_true")
    args = ap.parse_args()

    y = make_np_data()
    csv = os.path.join(HERE, "np_smoke_data.csv")
    if args.emit_dofile:
        np.savetxt(csv, y, delimiter=",", header="y", comments="")
        do = os.path.join(HERE, "np_dfgls.do")
        with open(do, "w", encoding="utf-8", newline="\n") as f:
            f.write(DOFILE.format(kmax=13))  # Schwert(200)=13
        print("do-file ->", do, "(数据 ->", csv, ")")
    elif args.parse:
        import csv as _csv
        path = os.path.join(HERE, "np_stata_baselines.csv")
        rows = list(_csv.DictReader(open(path)))
        print(f"k*, MAIC, sigma2_k  ({len(rows)} rows) — 对照 C++ maic/sigma2_k 1e-10")
        for r in rows[:5]:
            print(r)
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
