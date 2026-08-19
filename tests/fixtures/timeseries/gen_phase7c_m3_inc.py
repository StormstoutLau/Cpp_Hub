# -*- coding: utf-8 -*-
# gen_phase7c_m3_inc.py - M3 协整临界值表 .inc 生成 (Phase 7C)
#
# 五件套 (输出到 include/cpphub/timeseries/cointegration/, 头文件消费):
#   1. ol1992_cv.inc          — OL1992 Johansen 表 (urca ca.jo 源码常量转录, 3 ecdet × 2 统计量 × q=1..11 × 10/5/1%)
#   2. mhm96_johansen_cv.inc  — MHM96 Johansen 表 (statsmodels coint_tables 转录, 3 det × 2 统计量 × n=1..12 × 90/95/99%)
#   3. mackinnon1994_coint.inc— MacKinnon 1994 p 值近似表 (statsmodels adfvalues 转录, 4 trend × N=1..6)
#   4. mackinnon2010_coint_cv.inc — MacKinnon 2010 协整 CV 响应面 (statsmodels adfvalues 转录, c/ct/ctt × N=1..12 × 3 水平 × 4 系数)
#   5. em2002_ect_cv.inc      — Ericsson-MacKinnon 2002 ECT t 响应面 (IFDP 655 PDF 转录, 4 case × n=1..12 × 3 水平 × 4 系数)
#
# 溯源:
#   OL1992  ← tests/fixtures/timeseries/ca_jo_source.txt (urca 1.3-4 deparse, 2026-08-19)
#   MHM96   ← statsmodels 0.14.4 coint_tables.py (ejcp0-2/tjcp0-2, LeSage johansen.m 移植)
#   1994    ← statsmodels 0.14.4 adfvalues.py (tau_star/min/max + smallp/largep, 已乘 scaling)
#   2010    ← statsmodels 0.14.4 adfvalues.py (tau_{c,ct,ctt}_2010, N=1..12)
#   EM2002  ← tests/fixtures/timeseries/em2002_ect_cv.csv (IFDP 655 密码解码转录, 三重验证)
#
# 双库 diff 冻结决策 (docs/phases/phase7/JOHANSEN_DUAL_LIB_DIFF.md §6):
#   统计量主对照 = statsmodels; cvt/cvm 默认 = MHM96 (SM 行为等价); OL1992 独立查表 API
#   ⚠️ SM det_order∈{0,1} 的 MHM96 表与其自身统计量分布不一致 (MC 裁决), 文档警示

import csv
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join("..", "..", "..", "include", "cpphub", "timeseries",
                      "cointegration")
OUTDIR = os.path.abspath(os.path.join(HERE, OUTDIR))

GUARD = ("// 本文件由 tests/fixtures/timeseries/gen_phase7c_m3_inc.py 自动生成, 勿手改\n"
         "// (生成日期 2026-08-19; 再生成: python gen_phase7c_m3_inc.py)\n")


def fmt(v):
    return repr(float(v))


def emit_array(fh, name, shape_comment, arr):
    """arr: numpy ndarray, 任意维; 输出嵌套花括号初始化的 constexpr 数组"""
    flat = arr.ravel()
    dims = arr.shape

    def build(dim_offsets, depth):
        if depth == len(dims):
            return fmt(flat[sum(o * s for o, s in zip(dim_offsets, strides))])
        parts = []
        for i in range(dims[depth]):
            parts.append(build(dim_offsets + [i], depth + 1))
        return "{\n" + ", ".join(parts) + "\n}"

    strides = []
    s = 1
    for d in reversed(dims):
        strides.insert(0, s)
        s *= d
    body = build([], 0)
    dim_str = "".join("[%d]" % d for d in dims)
    fh.write("// %s\n" % shape_comment)
    fh.write("inline constexpr double %s%s = %s;\n\n" % (name, dim_str, body))


# ---------------------------------------------------------------------------
# 1. OL1992 (urca ca.jo 源码转录)
# ---------------------------------------------------------------------------
def parse_urca_arrays():
    src = open(os.path.join(HERE, "ca_jo_source.txt")).read()
    tables = {}
    for name in ("cv.none", "cv.const", "cv.trend"):
        m = re.search(re.escape(name) + r"\s*<-\s*array\(c\(([^)]+)\)",
                      src, re.S)
        vals = [float(v) for v in m.group(1).replace("\n", " ").split(",")]
        assert len(vals) == 66, (name, len(vals))
        # R column-major, dim c(11,3,2): [,,1]=maxeig, [,,2]=trace
        a = np.array(vals).reshape((2, 3, 11)).transpose(0, 2, 1)
        tables[name] = {"maxeig": a[0], "trace": a[1]}  # 各 (11 q × 3 pct)
    return tables


# ---------------------------------------------------------------------------
# 2. MHM96 (statsmodels coint_tables)
# ---------------------------------------------------------------------------
def load_mhm96():
    from statsmodels.tsa import coint_tables as ct
    return {
        "trace": {0: ct.tjcp0.copy(), 1: ct.tjcp1.copy(), 2: ct.tjcp2.copy()},
        "maxeig": {0: ct.ejcp0.copy(), 1: ct.ejcp1.copy(), 2: ct.ejcp2.copy()},
    }


# ---------------------------------------------------------------------------
# 3. MacKinnon 1994 p 值表 (statsmodels adfvalues, 已含 scaling)
# ---------------------------------------------------------------------------
def load_mk1994():
    from statsmodels.tsa import adfvalues as av
    trends = ("n", "c", "ct", "ctt")
    out = {
        "star": np.array([av._tau_stars[t] for t in trends]),   # 4×6
        "min": np.array([av._tau_mins[t] for t in trends]),     # 4×6
        "max": np.array([av._tau_maxs[t] for t in trends]),     # 4×6
        "smallp": np.array([av._tau_smallps[t] for t in trends]),  # 4×6×3
        "largep": np.array([av._tau_largeps[t] for t in trends]),  # 4×6×4
    }
    return out


# ---------------------------------------------------------------------------
# 4. MacKinnon 2010 协整 CV 响应面 (statsmodels adfvalues)
# ---------------------------------------------------------------------------
def load_mk2010():
    from statsmodels.tsa import adfvalues as av
    return {"c": av.tau_c_2010.copy(), "ct": av.tau_ct_2010.copy(),
            "ctt": av.tau_ctt_2010.copy()}  # 各 12×3×4 (N × 水平 × [β∞,β1,β2,β3])


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # ---- 1. ol1992_cv.inc ----
    ur = parse_urca_arrays()
    ecdets = ("none", "const", "trend")
    trace = np.stack([ur["cv." + e]["trace"] for e in ecdets])   # 3×11×3
    maxeig = np.stack([ur["cv." + e]["maxeig"] for e in ecdets])
    with open(os.path.join(OUTDIR, "ol1992_cv.inc"), "w", encoding="utf-8") as fh:
        fh.write("// =============================================================================\n")
        fh.write("// ol1992_cv.inc - Osterwald-Lenum 1992 Johansen 检验临界值表 (urca 转录)\n")
        fh.write("//\n")
        fh.write("// Phase 7C v1.7 M3 (PHASE7C_SPEC.md §5.2; 消费方: osterwald_lenum_cv.hpp)\n")
        fh.write("// 溯源: urca 1.3-4 ca.jo 源码常量 deparse (F:/R/win-library/4.6,\n")
        fh.write("//       ca_jo_source.txt, 2026-08-19); 论文锚: OL1992 OBES 54(3):461-472\n")
        fh.write(GUARD)
        fh.write("//\n")
        fh.write("// 布局: [ecdet][q-1][pct], ecdet 0=none/1=const/2=trend (urca 语义:\n")
        fh.write("//   none=短回归无约束常数 H1*, const=常数限入协整关系 H1, trend=趋势限入\n")
        fh.write("//   协整关系 H* + 无约束常数); q = N - r (检验 r 时的自由方向数);\n")
        fh.write("//   pct: 0=10%, 1=5%, 2=1% (urca cval 列序 10pct/5pct/1pct)\n")
        fh.write("// ⚠️ 双库对照注意 (JOHANSEN_DUAL_LIB_DIFF.md): urca ecdet 情形与\n")
        fh.write("//   statsmodels det_order 情形仅 det0↔none 可映射, 余不对应\n")
        fh.write("// =============================================================================\n\n")
        emit_array(fh, "OL1992_TRACE",
                   "迹统计量临界值 [ecdet 3][q 11][pct 3] (10/5/1%)", trace)
        emit_array(fh, "OL1992_MAXEIG",
                   "最大特征值统计量临界值 [ecdet 3][q 11][pct 3] (10/5/1%)", maxeig)
        # static_assert 锚 (转录保真; 置于数组声明后): urca cv.none trace q=3 = 28.71/31.52/37.22
        fh.write("static_assert(OL1992_TRACE[0][2][0] == 28.71, \"OL1992 anchor: none trace q3 10pct\");\n")
        fh.write("static_assert(OL1992_TRACE[0][2][1] == 31.52, \"OL1992 anchor: none trace q3 5pct\");\n")
        fh.write("static_assert(OL1992_TRACE[0][2][2] == 37.22, \"OL1992 anchor: none trace q3 1pct\");\n")
        fh.write("static_assert(OL1992_TRACE[1][0][1] == 9.24, \"OL1992 anchor: const trace q1 5pct\");\n")
        fh.write("static_assert(OL1992_TRACE[2][1][1] == 25.32, \"OL1992 anchor: trend trace q2 5pct\");\n")
        fh.write("static_assert(OL1992_MAXEIG[2][1][1] == 18.96, \"OL1992 anchor: trend maxeig q2 5pct\");\n")

    # ---- 2. mhm96_johansen_cv.inc ----
    mhm = load_mhm96()
    mtrace = np.stack([mhm["trace"][d] for d in (0, 1, 2)])   # 3×12×3
    mmaxeig = np.stack([mhm["maxeig"][d] for d in (0, 1, 2)])
    with open(os.path.join(OUTDIR, "mhm96_johansen_cv.inc"), "w", encoding="utf-8") as fh:
        fh.write("// =============================================================================\n")
        fh.write("// mhm96_johansen_cv.inc - MacKinnon-Haug-Michelis 1996 Johansen 检验临界值表\n")
        fh.write("//\n")
        fh.write("// Phase 7C v1.7 M3 (PHASE7C_SPEC.md §5.2; 消费方: johansen_test.hpp —\n")
        fh.write("//   JohansenResult.cvt/cvm 默认表源, 与 statsmodels coint_johansen 逐位一致)\n")
        fh.write("// 溯源: statsmodels 0.14.4 coint_tables.py (LeSage johansen.m 移植,\n")
        fh.write("//       MacKinnon 1996 johdist.f 计算); 论文锚: MHM96 (Queen's DP)\n")
        fh.write(GUARD)
        fh.write("//\n")
        fh.write("// 布局: [det_order][n-1][pct], det_order 0=-1/1=0/2=1 (statsmodels 语义);\n")
        fh.write("//   n = N - r; pct: 0=90%, 1=95%, 2=99% (statsmodels c_sjt/c_sja 列序)\n")
        fh.write("// ⚠️ MC 裁决 (JOHANSEN_DUAL_LIB_DIFF.md §5): det_order∈{0,1} 的本表与其\n")
        fh.write("//   自身统计量分布不一致 (χ²(1) 型 q=1 行为受限情形特征); 实际推断建议\n")
        fh.write("//   用 OL1992 查表 (det0→none 表; det1→trend 表近似) 或 det_order=-1\n")
        fh.write("// =============================================================================\n\n")
        emit_array(fh, "MHM96_TRACE",
                   "迹统计量临界值 [det_order 3][n 12][pct 3] (90/95/99%)", mtrace)
        emit_array(fh, "MHM96_MAXEIG",
                   "最大特征值统计量临界值 [det_order 3][n 12][pct 3] (90/95/99%)", mmaxeig)
        # static_assert 锚 (转录保真; 置于数组声明后): tjcp1[2] = 27.0669/29.7961/35.4628
        fh.write("static_assert(MHM96_TRACE[1][2][0] == 27.0669, \"MHM96 anchor: det0 trace n3 90pct\");\n")
        fh.write("static_assert(MHM96_TRACE[1][2][1] == 29.7961, \"MHM96 anchor: det0 trace n3 95pct\");\n")
        fh.write("static_assert(MHM96_TRACE[1][2][2] == 35.4628, \"MHM96 anchor: det0 trace n3 99pct\");\n")
        fh.write("static_assert(MHM96_TRACE[0][1][1] == 12.3212, \"MHM96 anchor: det-1 trace n2 95pct\");\n")
        fh.write("static_assert(MHM96_MAXEIG[2][2][0] == 21.8731, \"MHM96 anchor: det1 maxeig n3 90pct\");\n")

    # ---- 3. mackinnon1994_coint.inc ----
    mk = load_mk1994()
    with open(os.path.join(OUTDIR, "mackinnon1994_coint.inc"), "w", encoding="utf-8") as fh:
        fh.write("// =============================================================================\n")
        fh.write("// mackinnon1994_coint.inc - MacKinnon 1994 p 值近似表 (N≥1, 协整检验适用)\n")
        fh.write("//\n")
        fh.write("// Phase 7C v1.7 M3 (PHASE7C_SPEC.md §5.1; 消费方: mackinnon_coint_cv.hpp)\n")
        fh.write("// 溯源: statsmodels 0.14.4 adfvalues.py (mackinnonp; smallp/largep 已乘\n")
        fh.write("//       scaling [1,1,1e-2]/[1,1e-1,1e-1,1e-2]); 论文锚: MacKinnon 1994\n")
        fh.write("//       JBES 12(2):167-176 Tables 1-4\n")
        fh.write(GUARD)
        fh.write("//\n")
        fh.write("// 布局: [trend][N-1], trend 0=n/1=c/2=ct/3=ctt, N=1..6 (表覆盖域;\n")
        fh.write("//   EG 协整 k_vars=N; N>6 statsmodels 越界, 本实现同判)\n")
        fh.write("// p 值分段: stat > max → 1; stat < min → 0; stat ≤ star → smallp 二次;\n")
        fh.write("//   否则 largep 三次; p = Φ(polyval) — 系数升幂 (与 7B mackinnon_cv.hpp 同约定)\n")
        fh.write("// =============================================================================\n\n")
        emit_array(fh, "MK1994_TAUSTAR",
                   "τ* 分段界 [trend 4][N 6]", mk["star"])
        emit_array(fh, "MK1994_TAUMIN",
                   "τ_min 左尾截断 [trend 4][N 6]", mk["min"])
        emit_array(fh, "MK1994_TAUMAX",
                   "τ_max 右尾截断 [trend 4][N 6] (N=1 n 情形为 inf → 以 1e300 编码)",
                   np.where(np.isinf(mk["max"]), 1e300, mk["max"]))
        emit_array(fh, "MK1994_SMALLP",
                   "small p 多项式系数 (升幂, 3 系数二次) [trend 4][N 6][3]", mk["smallp"])
        emit_array(fh, "MK1994_LARGEP",
                   "large p 多项式系数 (升幂, 4 系数三次) [trend 4][N 6][4]", mk["largep"])
        # static_assert 锚 (转录保真; 置于数组声明后)
        fh.write("static_assert(MK1994_TAUSTAR[1][1] == -2.62, \"MK1994 anchor: tau_star_c N=2\");\n")
        fh.write("static_assert(MK1994_TAUMAX[1][1] == 0.92, \"MK1994 anchor: tau_max_c N=2\");\n")
        fh.write("static_assert(MK1994_SMALLP[1][1][0] == 2.92, \"MK1994 anchor: c smallp N=2 d0\");\n")
        fh.write("static_assert(MK1994_SMALLP[1][1][1] == 1.5012, \"MK1994 anchor: c smallp N=2 d1\");\n")
        fh.write("static_assert(MK1994_SMALLP[1][1][2] == 0.039796, \"MK1994 anchor: c smallp N=2 d2 (scaling 1e-2)\");\n")
        fh.write("static_assert(MK1994_LARGEP[2][1][1] == 0.5272, \"MK1994 anchor: ct largep N=2 d1 (scaling 1e-1)\");\n")

    # ---- 4. mackinnon2010_coint_cv.inc ----
    mk10 = load_mk2010()
    with open(os.path.join(OUTDIR, "mackinnon2010_coint_cv.inc"), "w", encoding="utf-8") as fh:
        fh.write("// =============================================================================\n")
        fh.write("// mackinnon2010_coint_cv.inc - MacKinnon 2010 协整检验 CV 响应面 (N=1..12)\n")
        fh.write("//\n")
        fh.write("// Phase 7C v1.7 M3 (PHASE7C_SPEC.md §5.1 Step3; 消费方: mackinnon_coint_cv.hpp)\n")
        fh.write("// 溯源: statsmodels 0.14.4 adfvalues.py tau_{c,ct,ctt}_2010; 论文锚:\n")
        fh.write("//       MacKinnon 2010 \"Critical Values for Cointegration Tests\"\n")
        fh.write("//       (Queen's WP 1227)\n")
        fh.write(GUARD)
        fh.write("//\n")
        fh.write("// ⚠️ spec §5.1 原文称 \"1994 协整响应面 5 系数 4 次\" — 实测 statsmodels\n")
        fh.write("//   coint() 的 cv 用 2010 表 4 系数 3 次 (polyval(val.T, 1/nobs)); 本表按\n")
        fh.write("//   实测转录 (JOHANSEN_DUAL_LIB_DIFF.md §7 偏离记录; 1994 机制仅用于 p 值)\n")
        fh.write("//\n")
        fh.write("// 布局: [N-1][水平][系数], 水平 0=1%/1=5%/2=10%; 系数升幂 [β∞, β1, β2, β3]\n")
        fh.write("//   CV(T) = β∞ + β1/T + β2/T² + β3/T³ (T = nobs-1, Stata egranger 约定)\n")
        fh.write("//   trend=\"n\" 无表 (statsmodels 返回 NaN; EG trend=n 时 cv 为 NaN)\n")
        fh.write("//   N=1 行与 7B mackinnon_tables.inc 的 ADF 表同源 (适用域分离, 决策 18)\n")
        fh.write("// =============================================================================\n\n")
        emit_array(fh, "T2010_C_CV",
                   "trend=c 响应面 [N 12][水平 3][系数 4]", mk10["c"])
        emit_array(fh, "T2010_CT_CV",
                   "trend=ct 响应面 [N 12][水平 3][系数 4]", mk10["ct"])
        emit_array(fh, "T2010_CTT_CV",
                   "trend=ctt 响应面 [N 12][水平 3][系数 4]", mk10["ctt"])
        # static_assert 锚 (转录保真; 置于数组声明后)
        fh.write("static_assert(T2010_C_CV[1][0][0] == -3.89644, \"M2010 anchor: c N=2 1pct b_inf\");\n")
        fh.write("static_assert(T2010_C_CV[1][1][1] == -6.1101, \"M2010 anchor: c N=2 5pct b1\");\n")
        fh.write("static_assert(T2010_CT_CV[2][0][3] == 104.244, \"M2010 anchor: ct N=3 1pct b3\");\n")
        fh.write("static_assert(T2010_CTT_CV[2][1][0] == -4.45311, \"M2010 anchor: ctt N=3 5pct b_inf\");\n")

    # ---- 5. em2002_ect_cv.inc ----
    rows = []
    with open(os.path.join(HERE, "em2002_ect_cv.csv")) as fh:
        for r in csv.DictReader(fh):
            rows.append(r)
    cases = {"n": 0, "c": 1, "ct": 2, "ctt": 3}
    sizes = {"1%": 0, "5%": 1, "10%": 2}
    em = np.full((4, 12, 3, 4), np.nan)
    for r in rows:
        em[cases[r["case"]], int(r["n"]) - 1, sizes[r["size"]]] = [
            float(r["theta_inf"]), float(r["theta1"]), float(r["theta2"]),
            float(r["theta3"])]
    with open(os.path.join(OUTDIR, "em2002_ect_cv.inc"), "w", encoding="utf-8") as fh:
        fh.write("// =============================================================================\n")
        fh.write("// em2002_ect_cv.inc - Ericsson-MacKinnon 2002 ECT t 检验响应面 (n=1..12)\n")
        fh.write("//\n")
        fh.write("// Phase 7C v1.7 M3 (PHASE7C_SPEC.md §5.4; 消费方: ericsson_mackinnon_cv.hpp)\n")
        fh.write("// 溯源: Fed IFDP 655 (1999 工作论文版, EM2002 同源) Table 2-5 密码解码转录\n")
        fh.write("//       (em2002_ifdp655.pdf, cmr10092 字体 cipher=ASCII+3; 生成器\n")
        fh.write("//       gen_em2002_tables.py; 三重验证: n=1 行 vs MacKinnon 2010 ≤0.0005,\n")
        fh.write("//       表 7 K_ctt(3) T=51 有限样本 CV 精确复现, 144 行结构完整)\n")
        fh.write("//       论文锚: Ericsson-MacKinnon 2002 Econometrics J. 5(2):285-318\n")
        fh.write("//       ⚠️ 正式发表版 (2002) 表与 WP 版或微差; 发表版付费墙, 本转录为可核验的\n")
        fh.write("//       官方免费版本\n")
        fh.write(GUARD)
        fh.write("//\n")
        fh.write("// 布局: [case][n-1][水平][系数], case 0=n(无确定项)/1=c(常数)/2=ct(常数+趋势)\n")
        fh.write("//   /3=ctt(常数+趋势+二次趋势); 水平 0=1%/1=5%/2=10%; 系数升幂 [θ∞,θ1,θ2,θ3]\n")
        fh.write("//   CV(T) = θ∞ + θ1/T + θ2/T² + θ3/T³  (T = ECM 回归有效样本量)\n")
        fh.write("//   n = ECM 系统变量总数 (含 LHS); θ3 原表以整数+尾随句点印刷 (−66. = −66)\n")
        fh.write("// =============================================================================\n\n")
        emit_array(fh, "EM2002_THETA",
                   "EM2002 ECT t 响应面 [case 4][n 12][水平 3][系数 4]", em)
        # static_assert 锚 (转录保真; 置于数组声明后)
        fh.write("static_assert(EM2002_THETA[1][0][0][0] == -3.4307, \"EM2002 anchor: c n=1 1pct t_inf\");\n")
        fh.write("static_assert(EM2002_THETA[1][0][0][1] == -6.52, \"EM2002 anchor: c n=1 1pct t1\");\n")
        fh.write("static_assert(EM2002_THETA[3][2][0][0] == -4.8399, \"EM2002 anchor: ctt n=3 1pct t_inf\");\n")
        fh.write("static_assert(EM2002_THETA[3][2][0][3] == -136.0, \"EM2002 anchor: ctt n=3 1pct t3\");\n")
        fh.write("static_assert(EM2002_THETA[0][0][1][1] == -0.35, \"EM2002 anchor: n n=1 5pct t1\");\n")

    print("wrote 5 .inc files to", OUTDIR)


if __name__ == "__main__":
    sys.exit(main())
