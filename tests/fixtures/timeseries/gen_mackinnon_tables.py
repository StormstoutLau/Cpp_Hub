# gen_mackinnon_tables.py - 从 arch 8.0.0 提取 MacKinnon 2010 系数表 → C++ .inc
# 表来源 (MIT 许可, 数值属公开科学常数):
#   arch/unitroot/critical_values/dickey_fuller.py: tau_2010 / tau_min / tau_max /
#     tau_star / tau_small_p / tau_large_p
#   arch/unitroot/critical_values/dfgls.py: dfgls_cv_approx / dfgls_tau_min /
#     dfgls_tau_max / dfgls_tau_star / dfgls_small_p / dfgls_large_p
#   arch/unitroot/critical_values/kpss.py: kpss_critical_values (插值表)
import math

import numpy as np
from arch.unitroot.critical_values.dickey_fuller import (
    tau_2010, tau_min, tau_max, tau_star, tau_small_p, tau_large_p,
)
from arch.unitroot.critical_values.dfgls import (
    dfgls_cv_approx, dfgls_tau_min, dfgls_tau_max, dfgls_tau_star,
    dfgls_small_p, dfgls_large_p,
)
from arch.unitroot.critical_values.kpss import kpss_critical_values

OUT = r"f:\Cpp_Hub\include\cpphub\timeseries\unit_root\mackinnon_tables.inc"


def fmt(x):
    x = float(x)
    if math.isinf(x):
        return "std::numeric_limits<Real>::infinity()"
    if math.isnan(x):
        return "std::numeric_limits<Real>::quiet_NaN()"
    return repr(x)


def dump_vec(v, name):
    v = np.asarray(v, dtype=float).ravel()
    body = ", ".join(fmt(x) for x in v)
    return f"inline constexpr Real {name}[] = {{{body}}};\n"


def dump_mat(m, name):
    m = np.atleast_2d(np.asarray(m, dtype=float))
    rows = []
    for r in m:
        rows.append("    {" + ", ".join(fmt(x) for x in r) + "}")
    return (f"inline constexpr Real {name}[][4] = {{\n" +
            ",\n".join(rows) + "\n};\n")


chunks = []
chunks.append(
"""// =============================================================================
// mackinnon_tables.inc - MacKinnon 2010 response surface 系数表 (自动生成)
//
// 生成脚本: tests/fixtures/timeseries/gen_mackinnon_tables.py (arch 8.0.0 提取)
// 来源: arch/unitroot/critical_values/{dickey_fuller,dfgls,kpss}.py (MIT)
// 重新生成: python tests/fixtures/timeseries/gen_mackinnon_tables.py
//
// 约定 (与 arch polyval 一致):
//   CV(p) 行向量 [c0, c1, c2, c3] (升幂): CV = c3/T³ + c2/T² + c1/T + c0
//   p 值多项式系数升幂: p = Φ(Σ cᵢ·statⁱ)
//   本文件被 mackinnon_cv.hpp 直接 #include (位于 namespace 内)
// =============================================================================
#pragma once

"""
)

chunks.append("// --- ADF (tau, MacKinnon 2010) — 每 trend 3 行 (1%/5%/10%), 4 系数 ---\n")
for reg in ["n", "c", "ct", "ctt"]:
    t = np.asarray(tau_2010[reg][0], dtype=float)  # (3, 4)
    chunks.append(dump_mat(t, f"TAU2010_{reg.upper()}"))
for reg in ["n", "c", "ct", "ctt"]:
    chunks.append(dump_vec(tau_min[reg][0], f"TAUMIN_{reg.upper()}"))
    chunks.append(dump_vec(tau_max[reg][0], f"TAUMAX_{reg.upper()}"))
    chunks.append(dump_vec(tau_star[reg][0], f"TAUSTAR_{reg.upper()}"))
    chunks.append(dump_vec(tau_small_p[reg][0], f"TAUSMALLP_{reg.upper()}"))
    chunks.append(dump_vec(tau_large_p[reg][0], f"TAULARGEP_{reg.upper()}"))

chunks.append("\n// --- DF-GLS (arch 独立模拟, 非 ERS 1996 原表, U7) ---\n")
for reg in ["c", "ct"]:
    chunks.append(dump_mat(dfgls_cv_approx[reg], f"DFGLS_CV_{reg.upper()}"))
    chunks.append(dump_vec(dfgls_tau_min[reg], f"DFGLS_MIN_{reg.upper()}"))
    chunks.append(dump_vec(dfgls_tau_max[reg], f"DFGLS_MAX_{reg.upper()}"))
    chunks.append(dump_vec(dfgls_tau_star[reg], f"DFGLS_STAR_{reg.upper()}"))
    chunks.append(dump_vec(dfgls_small_p[reg], f"DFGLS_SMALLP_{reg.upper()}"))
    chunks.append(dump_vec(dfgls_large_p[reg], f"DFGLS_LARGEP_{reg.upper()}"))

chunks.append("\n// --- KPSS 1992 插值表 (y=pct, x=quantile; 100,000,000 次模拟) ---\n")
for reg in ["c", "ct"]:
    tbl = np.asarray(kpss_critical_values[reg], dtype=float)  # (n, 2)
    chunks.append(f"inline constexpr Size KPSS_N_{reg.upper()} = {tbl.shape[0]};\n")
    chunks.append(dump_vec(tbl[:, 0], f"KPSS_Y_{reg.upper()}"))
    chunks.append(dump_vec(tbl[:, 1], f"KPSS_X_{reg.upper()}"))

with open(OUT, "w", encoding="utf-8", newline="\n") as f:
    f.write("".join(chunks))

# 自检: 打印关键值
print("tau_2010 c 5% T=100:", np.polyval(np.asarray(tau_2010["c"][0][1])[::-1], 0.01))
print("dfgls ct 5% T=100:", np.polyval(np.asarray(dfgls_cv_approx["ct"][1])[::-1], 0.01))
print("kpss c table rows:", np.asarray(kpss_critical_values["c"]).shape)
print("written:", OUT)
