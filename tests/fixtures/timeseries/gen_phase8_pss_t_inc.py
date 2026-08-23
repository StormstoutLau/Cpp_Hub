# -*- coding: utf-8 -*-
# gen_phase8_pss_t_inc.py - PSS 2001 Table CII (t 界) 转录生成 (Phase 8 M1 Step 4)
#
# 来源: R ARDL 0.2.5 crit_val_bounds_pss2001$t (export_pss2001_t_bounds.R 导出,
#   pss2001_t_bounds.csv) — PSS 2001 Table CII(i)/(iii)/(v) 转录, 与原文逐位
#   一致 (case iii k=1 5% = -2.86/-3.22, 原文锚)
#
# ⚠️ 与 statsmodels F 表 (pss_response_surface.inc) 的差异 (PSS9 对偶面):
#   - 档位: t 表 {10, 5, 2.5, 1}% vs F 表 {10, 5, 1, 0.1}%
#     (F 表第 3 档是 1% 无 2.5%; t 表第 4 档是 1% 无 0.1% — 两表第 3/4 档语义不同)
#   - k 域: t 表原表 k=0..10, 本转录仅 k=1..10 (bounds 检验语义 k>=1,
#     与 F 表键域一致; k=0 行不消费)
#   - 精度: 原表两位小数 (非响应面), 转录逐位保留
#
# 用法: & C:\Users\Peng\.conda\envs\open-webui\python.exe gen_phase8_pss_t_inc.py
# 产出: ../../include/cpphub/timeseries/ardl/pss_t_bounds.inc
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "..", "include", "cpphub",
                   "timeseries", "ardl", "pss_t_bounds.inc")

CASES = {"i": 0, "iii": 1, "v": 2}          # 表内 case 索引
ALPHAS = [0.100, 0.050, 0.025, 0.010]       # 档序 (PSS 原表)

vals = {}  # (case_idx, k, i1, alpha_idx) -> value
with open(os.path.join(HERE, "pss2001_t_bounds.csv")) as f:
    for row in csv.DictReader(f):
        k = int(row["k"])
        if k < 1 or k > 10:
            continue
        ci = CASES[row["case"]]
        ai = ALPHAS.index(float(row["alpha"]))
        vals[(ci, k, 0, ai)] = float(row["I0"])
        vals[(ci, k, 1, ai)] = float(row["I1"])

assert len(vals) == 3 * 10 * 2 * 4, f"expected 240 values, got {len(vals)}"

lines = [
    "// ============================================================================",
    "// pss_t_bounds.inc - PSS 2001 Table CII t 检验临界值 (自动生成, 勿手改)",
    "//",
    "// 来源: R ARDL 0.2.5 crit_val_bounds_pss2001$t (PSS 2001 Table CII(i)/(iii)/(v)",
    "// 转录, export_pss2001_t_bounds.R 导出) — statsmodels 无 t 表 (bounds_t_test",
    "// 不存在的库, probe [D]); R ARDL 是唯一软件对照源, 表值与 PSS 原文逐位一致",
    "// (case iii k=1 5%: -2.86/-3.22 = 原文锚)",
    "//",
    "// 布局: PSS_T_CRIT[case_idx][k-1][i1][alpha_idx]",
    "//   case_idx: 0=I, 1=III, 2=V (t 检验仅此三 case, A4)",
    "//   i1: 0=I(0) 下界, 1=I(1) 上界 (均为负值, I1 更负)",
    "//   alpha_idx: 0={10%}, 1={5%}, 2={2.5%}, 3={1%}",
    "// ⚠️ 档位 {10,5,2.5,1}% ≠ F 表 {10,5,1,0.1}% (PSS9 对偶面: 两表第 3/4 档",
    "// 语义不同 — F 表第 3 档 1% 无 2.5%; t 表第 4 档 1% 无 0.1%)",
    "// k 域 1..10 (原表含 k=0 行, bounds 语义不消费 — 与 F 表键域一致)",
    "// ============================================================================",
    "// (namespace 由消费方 uecm_bounds.hpp 的 detail 包装)",
    "",
    "inline constexpr double PSS_T_CRIT[3][10][2][4] = {",
]
for ci in range(3):
    cname = ["I", "III", "V"][ci]
    lines.append(f"    {{  // Case {cname}")
    for k in range(1, 11):
        lines.append(f"        {{  // k = {k}")
        for i1 in range(2):
            iname = ["I(0)", "I(1)"][i1]
            row = ", ".join(repr(vals[(ci, k, i1, ai)]) for ai in range(4))
            lines.append(f"            {{ {row} }},  // {iname}")
        lines.append("        },")
    lines.append("    },")
lines += [
    "};",
    "",
]

with open(OUT, "w", newline="\n") as f:
    f.write("\n".join(lines))
print(f"written: {os.path.normpath(OUT)}  values = {len(vals)}")
# 往返校验
chk = {(1, 1, 0, 1): -2.86, (1, 1, 1, 1): -3.22,   # iii k=1 5%
       (0, 1, 0, 1): -1.95, (2, 10, 1, 3): None}   # i k=1 5% / v k=10 1%
assert vals[(1, 1, 0, 1)] == -2.86 and vals[(1, 1, 1, 1)] == -3.22
assert vals[(0, 1, 0, 1)] == -1.95
print("anchor check ✓ (iii k=1 5% = -2.86/-3.22; i k=1 5% I0 = -1.95)")
