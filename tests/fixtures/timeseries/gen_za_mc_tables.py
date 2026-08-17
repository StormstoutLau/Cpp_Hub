# -*- coding: utf-8 -*-
# gen_za_mc_tables.py - statsmodels 0.14.6 vs arch 8.0.0 ZA MC 表交叉比对 + C++ 字面量生成
# 零手抄: 两库均从本地安装程序化读取 (2026-08-17, M0 za_mc_cv.inc 生成器, 沿 gen_mackinnon_tables.py 惯例)
# 复验: 库升版后重跑 — 144 值双库比对 + 锚索引输出; 输出块经拼装写入 critical_values/za_mc_cv.inc
import ast
import re
import sys

import numpy as np

SM = r"C:\Users\Peng\.conda\envs\open-webui\Lib\site-packages\statsmodels\tsa\stattools.py"
src = open(SM, encoding="utf-8").read()

# 解析 ZivotAndrewsUnitRoot.__init__ 内三块 self._x = ( (pct, val), ... )
def extract(name):
    m = re.search(r"self\._" + name + r"\s*=\s*(\(.*?\))\s*\n\s*self\._za_critical_values",
                  src, re.S)
    if m is None:
        # 末块 (_ct) 后面跟的不是 _za_critical_values 赋值行, 放宽
        m = re.search(r"self\._" + name + r"\s*=\s*(\(.*?\))\s*\)", src, re.S)
    return ast.literal_eval(m.group(1))

sm_tables = {k: extract(k) for k in ("c", "t", "ct")}

# arch 侧
from arch.unitroot.critical_values.zivot_andrews import za_critical_values  # noqa: E402

n_diff = 0
for k in ("c", "t", "ct"):
    sm = np.asarray(sm_tables[k], dtype=float)
    ar = np.asarray(za_critical_values[k], dtype=float)
    assert sm.shape == ar.shape, (k, sm.shape, ar.shape)
    for i in range(sm.shape[0]):
        if sm[i, 0] != ar[i, 0] or sm[i, 1] != ar[i, 1]:
            n_diff += 1
            print(f"DIFF {k}[{i}]: sm={tuple(sm[i])} arch={tuple(ar[i])}")
print(f"cross-check: {3*sm.shape[0]} 值, 差异 {n_diff}, 每模型 {sm.shape[0]} 点")

# 锚索引 (按 pct 查)
for k in ("c", "t", "ct"):
    sm = sm_tables[k]
    for pct in (0.001, 0.1, 1.0, 5.0, 10.0, 99.9):
        idx = [i for i, p in enumerate(sm) if p[0] == pct]
        print(f"{k} pct={pct}: idx={idx} val={[sm[i][1] for i in idx]}")

# 生成 C++ 字面量
out = open(r"F:\Cpp_Hub\build\za_mc_gen_out.txt", "w", encoding="utf-8")
name_map = {"c": "ZA_MC_C_GOLDEN", "t": "ZA_MC_T_GOLDEN", "ct": "ZA_MC_CT_GOLDEN"}
for k in ("c", "t", "ct"):
    sm = sm_tables[k]
    out.write(f"constexpr std::array<ZaMcPoint, {len(sm)}> {name_map[k]}{{{{\n")
    for p, v in sm:
        out.write(f"    {{{p:.3f}, {v}}},\n")
    out.write("}};\n\n")
out.close()
print("C++ block -> F:/Cpp_Hub/build/za_mc_gen_out.txt")
