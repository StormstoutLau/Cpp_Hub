# -*- coding: utf-8 -*-
# gen_phase7c_m2_inc.py - VAR/DY 基准 .inc 生成 (Phase 7C M2)
#
# 输入: var_smoke_data.csv + var_statsmodels_baselines.json (主基准)
#       + var_r_values.txt (vars 交叉 + Spillover GFEVD/DY 主基准)
# 输出: tests/unit/timeseries/var_baseline.inc
#   namespace cpphub::v1::timeseries::var_baseline (嵌套先例 arima_baseline)
#
# 锚点容差分层 (spec §5):
#   statsmodels 主: 系数/IC/IRF/FEVD 1e-10
#   vars 交叉: IC 1e-8 (VARselect 轨迹与 statsmodels 逐位一致, 双保险)
#   Spillover 主: GFEVD/DY 指数 1e-8
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "unit", "timeseries", "var_baseline.inc")


def fmt(v):
    return repr(float(v))


def arr(name, vals):
    return (f"inline constexpr std::array<Real, {len(vals)}> {name} = {{\n"
            "    " + ", ".join(fmt(v) for v in vals) + "\n};\n")


def main():
    with open(os.path.join(HERE, "var_statsmodels_baselines.json")) as f:
        sm = json.load(f)
    rvals = {}
    with open(os.path.join(HERE, "var_r_values.txt")) as f:
        for line in f:
            if not line.strip():
                continue
            key, _, rest = line.partition("\t")
            rvals[key] = [float(x) for x in rest.split()]

    ycols = [[], [], []]
    with open(os.path.join(HERE, "var_smoke_data.csv")) as f:
        for row in csv.reader(f):
            if row and row[0].startswith("y"):
                continue
            for j in range(3):
                ycols[j].append(float(row[j]))

    fx = sm["fixtures"]
    v2 = fx["var2_c"]

    parts = []
    parts.append(
        "// var_baseline.inc - Phase 7C M2 VAR/DY 黄金基准 (自动生成, 勿手改)\n"
        "// 源: verify_var.py (statsmodels 0.14.4 主基准) + dump_var_r_values.R\n"
        "//     (vars 交叉 + Spillover 0.1.1 裁剪装载 GFEVD/DY 主基准)\n"
        "// DGP: 稳定 VAR(2) K=3 T=250 seed=42 (verify_var.py gen_var2)\n"
        "// 锚点: SM 系数/IC/IRF/FEVD 1e-10; vars IC 1e-8; Spillover 1e-8\n"
        "#pragma once\n"
        "#include <array>\n\n"
        "#include \"cpphub/core/types.hpp\"\n\n"
        "namespace cpphub {\n"
        "inline namespace v1 {\n"
        "namespace timeseries {\n"
        "namespace var_baseline {\n\n")

    parts.append("// --- 夹具数据 250×3 (行=时间) ---\n")
    for j in range(3):
        parts.append(arr(f"Y{j + 1}", ycols[j]))

    parts.append("\n// --- statsmodels 主基准: VAR(2) trend=c ---\n")
    parts.append(arr("SM_PARAMS", [x for row in v2["params"] for x in row]))
    parts.append(arr("SM_SIGMA_U", [x for row in v2["sigma_u"] for x in row]))
    parts.append(arr("SM_SIGMA_U_MLE", [x for row in v2["sigma_u_mle"] for x in row]))
    parts.append(arr("SM_SCALARS", [v2["loglik"], v2["aic"], v2["bic"],
                                    v2["hqic"], v2["fpe"], v2["detomega"]]))
    parts.append("// statsmodels roots = 1/|eig(comp)| 升序; max_eig = 1/roots[0]\n")
    parts.append(arr("SM_ROOTS_ASC", v2["roots"]))
    parts.append(arr("SM_MAX_EIG", [1.0 / v2["roots"][0]]))
    parts.append(arr("SM_CHOL_LOWER", [x for row in fx["chol_sigma_u"] for x in row]))

    parts.append("\n// --- select_order maxlags=4 轨迹 (p=0..4, V5 同样本) ---\n")
    so = fx["select_order_max4"]["ics"]
    parts.append(arr("SM_SEL_AIC", so["aic"]))
    parts.append(arr("SM_SEL_BIC", so["bic"]))
    parts.append(arr("SM_SEL_HQIC", so["hqic"]))
    parts.append(arr("SM_SEL_FPE", so["fpe"]))

    parts.append("\n// --- IRF 正交化 (Ψ_h, h=0/1/2/10, V3 行=响应) ---\n")
    irfs = fx["orth_irf_h10"]  # 11×3×3
    for h in (0, 1, 2, 10):
        parts.append(arr(f"SM_IRF_H{h}", [x for row in irfs[h] for x in row]))
    parts.append(arr("SM_PHI_H1", [x for row in fx["phi_ma_rep_h10"][1] for x in row]))

    parts.append("\n// --- FEVD Cholesky 轨 (V7 行和=1; decomp 布局 eq×period×comp) ---\n")
    dec = fx["fevd_orth_periods10"]  # 3×10×3
    parts.append(arr("SM_FEVD_H10", [x for eq in dec for x in eq[9]]))
    parts.append(arr("SM_FEVD_H1", [x for eq in dec for x in eq[0]]))
    decr = fx["reorder_201"]["fevd"]
    parts.append(arr("SM_REORDER_FEVD", [x for eq in decr for x in eq[9]]))
    parts.append(arr("SM_REORDER_PARAMS",
                     [x for row in fx["reorder_201"]["params"] for x in row]))

    parts.append("\n// --- trend=n / ct 单点 + p=0 ---\n")
    parts.append(arr("SM_N_PARAMS", [x for row in fx["var2_n"]["params"] for x in row]))
    parts.append(arr("SM_N_SCALARS", [fx["var2_n"]["loglik"], fx["var2_n"]["aic"],
                                      fx["var2_n"]["bic"], fx["var2_n"]["hqic"],
                                      fx["var2_n"]["fpe"]]))
    parts.append(arr("SM_CT_SCALARS", [fx["var2_ct"]["loglik"], fx["var2_ct"]["aic"],
                                       fx["var2_ct"]["bic"], fx["var2_ct"]["hqic"],
                                       fx["var2_ct"]["fpe"]]))
    parts.append(arr("SM_P0_SCALARS", [fx["var0_c"]["aic"], fx["var0_c"]["bic"],
                                       fx["var0_c"]["hqic"], fx["var0_c"]["fpe"]]))
    parts.append(arr("SM_N_SIGMA_MLE",
                     [x for row in fx["var2_n"]["sigma_u_mle"] for x in row]))

    parts.append("\n// --- R vars 交叉 ---\n")
    parts.append("// vars_coef 行序 = 回归元 (y1.l1 y2.l1 y3.l1 y1.l2 ... const)\n"
                 "// ⚠️ 与 statsmodels params 布局不同: SM trend 列在最前\n")
    parts.append(arr("R_VARS_COEF", rvals["vars_coef"]))
    parts.append(arr("R_VARS_LOGLIK", rvals["vars_loglik"]))
    parts.append(arr("R_VARS_CHOL", rvals["vars_chol_lower"]))
    parts.append(arr("R_VARS_FEVD_H10", rvals["vars_fevd_orth_h10"]))
    parts.append("// VARselect 轨迹 p=1..4 (无 p=0; 与 SM_SEL_*[1..4] 对照)\n")
    parts.append(arr("R_SEL_AIC", rvals["vars_varselect_aic"]))
    parts.append(arr("R_SEL_SC", rvals["vars_varselect_sc"]))
    parts.append(arr("R_SEL_FPE", rvals["vars_varselect_fpe"]))

    parts.append("\n// --- R Spillover 0.1.1 主基准 (GFEVD/DY) ---\n")
    parts.append(arr("SP_GFEVD_H10", rvals["gfevd_dy_h10"]))
    parts.append(arr("SP_GFEVD_H10_RAW", rvals["gfevd_dy_h10_raw"]))
    parts.append(arr("SP_TABLE", rvals["gspillover_table_h10"]))
    parts.append(arr("SP_TO", rvals["gspillover_to"]))
    parts.append(arr("SP_FROM", rvals["gspillover_from"]))
    parts.append(arr("SP_TCI", rvals["gspillover_tci"]))
    parts.append(arr("SP_NET", rvals["net"]))
    parts.append(arr("SP_GFEVD_H50", rvals["gfevd_dy_h50"]))
    parts.append(arr("SP_TCI_H50", rvals["gspillover_tci_h50"]))
    parts.append("// roll.spillover w=150 step=1 逐窗口 TCI (101 个)\n")
    parts.append(arr("SP_ROLL_TCI", rvals["roll_tci_all"]))

    parts.append("\n}  // namespace var_baseline\n"
                 "}  // namespace timeseries\n"
                 "}  // inline namespace v1\n"
                 "}  // namespace cpphub\n")

    with open(OUT, "w") as f:
        f.write("".join(parts))
    print("written:", os.path.normpath(OUT))


if __name__ == "__main__":
    main()
