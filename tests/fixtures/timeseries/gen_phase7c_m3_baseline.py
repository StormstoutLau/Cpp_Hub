# -*- coding: utf-8 -*-
# gen_phase7c_m3_baseline.py - M3 协整测试基准 .inc 生成 (Phase 7C)
#
# 数据源 (verify_*.py/.R 生成, 固化):
#   eg_statsmodels_baselines.txt      (statsmodels coint 0.14.4)
#   johansen_statsmodels_baselines.txt (statsmodels coint_johansen 0.14.4)
#   vecm_statsmodels_baselines.txt    (statsmodels VECM 0.14.4)
#   po_urca_baselines.txt             (urca 1.3-4 ca.po)
#
# 输出: tests/unit/timeseries/coint_baseline.inc
#   命名: EG_<PAIR>_<TREND>_{T,P,CV} / JO_DET<d>_K<k>_{EIG,LR1,LR2,CVT,CVM,EVEC,NOBS,RANK_*}
#         / VECM_<DET>_R<r>_K<k>_{ALPHA,BETA,GAMMA,DET_COEF,DET_COEF_COINT,SIGMA_U,LLF,STDERR_ALPHA,...}
#         / PO_<PAIR>_<TYPE>_<DEMEAN>_<LAG>{,_CVAL,_LMAX}

import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(HERE, "..", "..", "unit", "timeseries"))


def fmt(v):
    v = float(v)
    if v != v:  # NaN → 合法 constexpr 表达式
        return "NQ"
    return repr(v)


def sanitize(name):
    # 标识符合法化: 负号 → M (jo_det-1 → JO_DETM1), 小数点 → _
    return name.replace("-", "M").replace(".", "_")


def array(name, vals):
    body = ",\n    ".join(fmt(v) for v in vals)
    return "inline constexpr double %s[] = {\n    %s,\n};\n" % (name, body)


def load(path):
    d = {}
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 2 and parts[1].strip():
            d[parts[0]] = [float(v) for v in parts[1].split()]
    return d


def main():
    eg = load(os.path.join(HERE, "eg_statsmodels_baselines.txt"))
    jo = load(os.path.join(HERE, "johansen_statsmodels_baselines.txt"))
    ve = load(os.path.join(HERE, "vecm_statsmodels_baselines.txt"))
    po = load(os.path.join(HERE, "po_urca_baselines.txt"))

    # 夹具数据 (coint_smoke_data.csv, T=250)
    import csv
    cols = [[], [], []]
    with open(os.path.join(HERE, "coint_smoke_data.csv")) as fh:
        for row in csv.reader(fh):
            if row and row[0] != "y1":
                for j in range(3):
                    cols[j].append(float(row[j]))

    out = []
    out.append("// =============================================================================\n")
    out.append("// coint_baseline.inc - M3 协整模块测试基准 (Phase 7C v1.7)\n")
    out.append("//\n")
    out.append("// 生成: tests/fixtures/timeseries/gen_phase7c_m3_baseline.py (勿手改)\n")
    out.append("// 源: eg/johansen/vecm = statsmodels 0.14.4; po = urca 1.3-4 (R 4.6.1)\n")
    out.append("// 数据: coint_smoke_data.csv (T=250; y1/y2 协整, y3 独立 RW)\n")
    out.append("// 容差: statsmodels 1e-10 / urca 1e-8 (spec §1.3)\n")
    out.append("// =============================================================================\n\n")
    out.append("#include <limits>\n\n")
    out.append("namespace cpphub {\ninline namespace v1 {\nnamespace timeseries {\n")
    out.append("namespace coint_baseline {\n\n")
    out.append("// NaN 编码 (EG trend=n 的 cv; §1.4-5)\n")
    out.append("inline constexpr double NQ =\n")
    out.append("    std::numeric_limits<double>::quiet_NaN();\n\n")

    out.append(array("Y1", cols[0]))
    out.append(array("Y2", cols[1]))
    out.append(array("Y3", cols[2]))

    # ---- EG ----
    for tag, vals in eg.items():
        name = sanitize("EG_" + tag[3:].upper())  # eg_y1_y2_c_t → EG_Y1_Y2_C_T
        out.append(array(name, vals))

    # ---- Johansen ----
    for tag, vals in jo.items():
        name = sanitize("JO_" + tag[3:].upper())  # jo_det0_k1_eig → JO_DET0_K1_EIG
        out.append(array(name, vals))

    # ---- VECM ----
    for tag, vals in ve.items():
        name = sanitize("VECM_" + tag[5:].upper())  # vecm_n_r1_k1_alpha → VECM_N_R1_K1_ALPHA
        out.append(array(name, vals))

    # ---- PO ----
    for tag, vals in po.items():
        name = sanitize("PO_" + tag[3:].upper())  # po_y1y2_Pu_none_short → PO_Y1Y2_PU_NONE_SHORT
        out.append(array(name, vals))

    # ---- urca 交叉 (JOHANSEN_DUAL_LIB_DIFF.md §2: det_order=0 ↔ ecdet="none",
    #      k_ar_diff = K−1; 1e-8 容差) ----
    ur = load(os.path.join(HERE, "johansen_urca_grid.txt"))
    for tag, vals in ur.items():
        if "_cvt" in tag or "_cvm" in tag:
            continue  # 表值另由 OL1992 static_assert 锚定
        name = sanitize("UR_" + tag[3:].upper())  # ur_K2_none_eig → UR_K2_NONE_EIG
        out.append(array(name, vals))

    out.append("}  // namespace coint_baseline\n")
    out.append("}  // namespace timeseries\n")
    out.append("}  // inline namespace v1\n")
    out.append("}  // namespace cpphub\n")

    path = os.path.join(OUT, "coint_baseline.inc")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("".join(out))
    print("wrote", path, "(%d arrays)" % (len(eg) + len(jo) + len(ve) + len(po)))


if __name__ == "__main__":
    main()

