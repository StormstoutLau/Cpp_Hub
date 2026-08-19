# -*- coding: utf-8 -*-
# verify_eg.py - Engle-Granger 基准生成: statsmodels coint (Phase 7C M3)
#
# 对照点 (spec §5.1, 容差表):
#   统计量/p 值/临界值 vs statsmodels coint 1e-10 (CI1/CI2)
#   p (1994 渐近) 与 cv (nobs-1 小样本修正) 分列断言 (CI2, issue #4138)
#   双方向输出 (CI3: LHS 选择改变统计量)
#
# 语义 (一手源码 stattools.py L1702-1839):
#   Step1: OLS y0 ~ [y1, trend 列] (add_trend prepend=False; "n" 无趋势列)
#   Step2: adfuller(resid, maxlag=None, autolag="aic", regression="n")
#   cv = mackinnoncrit(N=k_vars, regression=trend, nobs=nobs-1)  ← 2010 响应面
#       (trend="n" 时 cv = NaN ×3)
#   p  = mackinnonp(t, regression=trend, N=k_vars)               ← 1994 渐近
#
# 输出: eg_statsmodels_baselines.txt (tag value 格式, C++ 测试读取)

import os

import numpy as np
from statsmodels.tsa.stattools import coint

HERE = os.path.dirname(os.path.abspath(__file__))


def emit(fh, tag, obj):
    a = np.asarray(obj, dtype=float).ravel()
    fh.write("%s\t%s\n" % (tag, " ".join("%.17g" % v for v in a)))


def main():
    data = np.loadtxt(os.path.join(HERE, "coint_smoke_data.csv"),
                      delimiter=",", skiprows=1)
    y1, y2, y3 = data[:, 0], data[:, 1], data[:, 2]

    # (tag, y0, y1): 协整对双方向 + 非协整对
    pairs = [
        ("y1_y2", y1, y2),   # 协整 (y1 ← y2)
        ("y2_y1", y2, y1),   # 协整反方向 (CI3: 统计量不同)
        ("y1_y3", y1, y3),   # 无协整
    ]
    out = os.path.join(HERE, "eg_statsmodels_baselines.txt")
    with open(out, "w") as fh:
        for tag, y0, y1 in pairs:
            for trend in ("n", "c", "ct", "ctt"):
                t, p, crit = coint(y0, y1, trend=trend)
                base = "eg_%s_%s" % (tag, trend)
                emit(fh, base + "_t", t)
                emit(fh, base + "_p", p)
                emit(fh, base + "_cv", crit)   # [1%, 5%, 10%]; "n" → NaN×3
    print("wrote", out)


if __name__ == "__main__":
    main()
