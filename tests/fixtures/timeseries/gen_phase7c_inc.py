# -*- coding: utf-8 -*-
# gen_phase7c_inc.py - M1/M4 baseline .inc 生成 (Phase 7C)
#
# 数据源: arima_smoke_data.csv / midas_smoke_data.csv (verify_arima.py /
#         verify_midas.R 生成, 固化)
# 基准值: 硬编码转录自
#   - arima_statsmodels_baselines.json (statsmodels innovations_mle)
#   - verify_arima.R 输出 2026-08-18 (R stats::arima CSS/CSS-ML + forecast)
#   - verify_midas.R 输出 2026-08-18 (midasr midas_u/midas_r/hAh, reltol=1e-12)
#
# 输出: tests/unit/timeseries/{arima_baseline.inc, midas_baseline.inc}
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "unit", "timeseries")


def fmt(v):
    return repr(float(v))


def read_smoke(name):
    with open(os.path.join(HERE, name)) as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    keys = None
    for ln in lines:  # #columns 可能在首或末 (两遍扫描)
        if ln.startswith("#columns="):
            keys = ln[len("#columns="):].split(",")
    cols = {k: [] for k in keys}
    for ln in lines:
        if ln.startswith("#"):
            continue
        for k, v in zip(keys, ln.split(",")):
            cols[k].append(float(v))
    return cols


def array(name, vals):
    body = ",\n    ".join(fmt(v) for v in vals)
    return f"inline constexpr double {name}[] = {{\n    {body},\n}};\n"


# ---------------------------------------------------------------------------
# arima_baseline.inc
# ---------------------------------------------------------------------------
A = read_smoke("arima_smoke_data.csv")

R_CSS = {
    # arma11: ar1, ma1, sigma2, n_cond
    "arma11": [0.507151481290735, -0.466362070410424, 0.942778544428568, 1],
    "arma22": [0.102005174085234, -0.655121846388472, -0.264709005166464,
               0.727013657335629, 0.987977426909865, 2],
    "arma21": [0.516042187287411, -0.289624762794409, -0.584213647784576,
               0.979580994346941, 2],
    "arma12": [0.259465951504607, -0.416783025587859, 0.303971376849179,
               1.00905066181426, 1],
}
R_CSSML = {
    "arma11": [0.457955151898394, -0.410772283239053, 0.948363373492533,
               -417.730619789911, 841.461239579821],
    "arma22": [0.0535081324641637, -0.846471559767537, -0.172678207733824,
               0.880970241416009, 0.99029917210364, -424.407990689178],
    "arma21": [0.514403490857483, -0.288902360743633, -0.581174028304493,
               0.975184796286224, -422.067419916139],
    "arma12": [0.259372168617419, -0.416205445250941, 0.302206727300177,
               1.00561612169720, -426.614393480063],
}
# d=1 退化路径 (stats::arima 对 d≥1 强制无均值, 漂移被 AR 近单位根吸收)
R_CSS_111D = [0.996315302474128, -0.953887362433824, 1.05118054207756, 2]
R_CSSML_111D = [0.998466178588142, -0.978150735056852, 1.03449142604403,
                -430.039497959073]
# forecast::Arima drift CSS-ML (漂移正解; loglik 与 statsmodels 一致)
R_DRIFT = [-0.00162289764516166, 0.0998774744538597, 0.357341429939636,
           1.03085532012423, -427.30303213856]

SM = {  # statsmodels innovations_mle: phi..., theta..., sigma2, loglik
    "arma11": [0.398025655212007, -0.357566899092779, 0.943329157923882,
               -416.931675303541],
    "arma22": [0.0529313131189783, -0.846000482922429, -0.1726026045894231,
               0.8803834556007019, 0.9891789183651445, -424.2381225582666],
    "arma21": [0.516702718758647, -0.288699155517682, -0.584423395693958,
               0.974215960453365, -421.919449727613],
    "arma12": [0.254826308370333, -0.412692316706425, 0.301259136110341,
               1.004784699118182, -426.489786734614],
    # 差分后 demean (nobs=299): phi, theta, sigma2, loglik
    "arima111d": [-0.00156870107191634, 0.0998201832766277,
                  1.020512361733972, -427.303045095959],
}

with open(os.path.join(OUT, "arima_baseline.inc"), "w") as f:
    f.write("// arima_baseline.inc - M1 基准 (自动生成 gen_phase7c_inc.py)\n")
    f.write("// 数据: arima_smoke_data.csv (verify_arima.py, %.17g)\n")
    f.write("// 基准: R stats::arima CSS/CSS-ML + forecast::Arima (verify_arima.R\n")
    f.write("//       2026-08-18) + statsmodels innovations_mle (.json)\n")
    f.write("// 语义: n.cond = d + p (q 无关, arma12 实测定案);\n")
    f.write("//       d>=1 stats::arima 强制无均值 → 漂移伪根路径 (退化对照);\n")
    f.write("//       漂移正解 = forecast::Arima include.drift (loglik 同 statsmodels)\n")
    f.write("#pragma once\n")
    f.write("namespace cpphub { inline namespace v1 { namespace timeseries {\n")
    f.write("namespace arima_baseline {\ninline namespace v1 {\n")
    f.write(f"inline constexpr Size T = {len(A['arma11'])};\n")
    f.write(f"inline constexpr Size T_DIFF = {len(A['arma11']) - 1};\n")
    for k in ["arma11", "arma22", "arma21", "arma12", "arima111d_level"]:
        f.write(array(k.upper(), A[k]))
    f.write(array("CSS_ARMA11", R_CSS["arma11"]))
    f.write(array("CSS_ARMA22", R_CSS["arma22"]))
    f.write(array("CSS_ARMA21", R_CSS["arma21"]))
    f.write(array("CSS_ARMA12", R_CSS["arma12"]))
    f.write(array("CSSML_ARMA11", R_CSSML["arma11"]))
    f.write(array("CSSML_ARMA22", R_CSSML["arma22"]))
    f.write(array("CSSML_ARMA21", R_CSSML["arma21"]))
    f.write(array("CSSML_ARMA12", R_CSSML["arma12"]))
    f.write(array("CSS_ARIMA111D", R_CSS_111D))
    f.write(array("CSSML_ARIMA111D", R_CSSML_111D))
    f.write(array("DRIFT_CSSML", R_DRIFT))
    f.write(array("SM_ARMA11", SM["arma11"]))
    f.write(array("SM_ARMA22", SM["arma22"]))
    f.write(array("SM_ARMA21", SM["arma21"]))
    f.write(array("SM_ARMA12", SM["arma12"]))
    f.write(array("SM_ARIMA111D", SM["arima111d"]))
    f.write("}  // inline namespace v1\n}  // namespace arima_baseline\n")
    f.write("}  // namespace timeseries\n}  // inline namespace v1\n}  // namespace cpphub\n")
print("written arima_baseline.inc")

# ---------------------------------------------------------------------------
# midas_baseline.inc
# ---------------------------------------------------------------------------
M = read_smoke("midas_smoke_data.csv")  # y, x1..x3 各 100
y = M["y"]
x = []
for i in range(len(y)):  # 低频宽表逐行展开 = 高频原序
    x.extend([M["x1"][i], M["x2"][i], M["x3"][i]])

W = {  # W1 权重逐点 (midasr 源函数直调, 1e-12)
    "NEALMON1": [0.45505423392341127, 0.27600434470659363,
                 0.16740509727844333, 0.10153632409155181],
    "NEALMON2": [0.1373751707747993, 0.14298155784893038,
                 0.14298155784893038, 0.1373751707747993,
                 0.12681326572987286, 0.11247327702266779],
    "NBETA1": [7.1054273576009934e-16, 0.44999999999999968,
               0.39999999999999969, 0.14999999999999988,
               1.5777218104420222e-31],
    "NBETA2": [0.029411766708790883, 0.079504500933496161,
               0.092774620989598502, 0.096618222736228873,
               0.092774620989598516, 0.079504500933496161,
               0.029411766708790883],
    "ALMONP": [0.14000000000000001, 0.16, 0.16000000000000003,
               0.14000000000000001, 0.10000000000000001],
    "POLYSTEP": [0.5, 0.5, 0.20000000000000001, 0.20000000000000001,
                 0.20000000000000001, 0.10000000000000001,
                 0.10000000000000001, 0.10000000000000001],
    "HARSTEP": [0.66499999999999992, 0.065000000000000002,
                0.0050000000000000001, 0.0050000000000000001],  # w1,w2,w6,w20
}
UMIDAS = {  # W3 midas_u: coef(intercept,w1..w4), sigma2, SSR, n_eff
    "COEF": [2.0021767720635286, 2.238984315040462, 1.4483864990838673,
             0.81248790537010396, 0.51146622724423074],
    "SIGMA2": 0.078004678227699617,
    "SSR": 7.3324397534037642,
}
NLS = {  # W4 midas_r nealmon: coef(mu,delta,l2), SSR, midas.coef(mu,w1..w4)
    "COEF": [2.0044540668461845, 5.0288962801347097, -0.49453592787709455],
    "SSR": 7.4968752407615744,
    "MIDAS_COEF": [2.0044540668461845, 2.2769678129701965, 1.3886176048512693,
                   0.84685380334276827, 0.51645705897047522],
}
NLS2 = [2.0044540625241631, 5.0288962672449653, -0.49453588858886466]  # W5
AR = [2.0353177125387552, 5.0689494875477887, -0.47979526023042079,
      -0.016548798378619139]  # W6: mu, delta, l2, rho1
HAH = [2.1080208397230265, 0.34853716215988106, 2]  # W7: stat, p, df
W_TRUE = [2.275271169617056, 1.3800217235329681, 0.83702548639221652,
          0.507681620457759]

with open(os.path.join(OUT, "midas_baseline.inc"), "w") as f:
    f.write("// midas_baseline.inc - M4 基准 (自动生成 gen_phase7c_inc.py)\n")
    f.write("// 数据: midas_smoke_data.csv (verify_midas.R, m=3, 低频宽表)\n")
    f.write("//       y 100 期 + x 300 (逐行展开) — CSV 与 C++ 同数组\n")
    f.write("// 基准: midasr 0.9 (verify_midas.R 2026-08-18, reltol=1e-12)\n")
    f.write("// 方向: w_1 ↔ h=0 期末最新 (W-dir Form A 恢复 λ* 定案)\n")
    f.write("#pragma once\n")
    f.write("namespace cpphub { inline namespace v1 { namespace timeseries {\n")
    f.write("namespace midas_baseline {\ninline namespace v1 {\n")
    f.write("inline constexpr Size N_LF = 100;\n")
    f.write("inline constexpr Size N_HF = 300;\n")
    f.write("inline constexpr Size M = 3;\n")
    f.write(array("Y", y))
    f.write(array("X", x))
    for k, v in W.items():
        f.write(array(k, v))
    f.write(array("UMIDAS_COEF", UMIDAS["COEF"]))
    f.write(f"inline constexpr double UMIDAS_SIGMA2 = {fmt(UMIDAS['SIGMA2'])};\n")
    f.write(f"inline constexpr double UMIDAS_SSR = {fmt(UMIDAS['SSR'])};\n")
    f.write(array("NLS_COEF", NLS["COEF"]))
    f.write(f"inline constexpr double NLS_SSR = {fmt(NLS['SSR'])};\n")
    f.write(array("NLS_MIDAS_COEF", NLS["MIDAS_COEF"]))
    f.write(array("NLS2_COEF", NLS2))
    f.write(array("AR_COEF", AR))
    f.write(array("HAH", HAH))
    f.write(array("W_TRUE", W_TRUE))
    f.write("}  // inline namespace v1\n}  // namespace midas_baseline\n")
    f.write("}  // namespace timeseries\n}  // inline namespace v1\n}  // namespace cpphub\n")
print("written midas_baseline.inc")
