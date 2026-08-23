# -*- coding: utf-8 -*-
# verify_ardl.py - ARDL 基准生成 + 语义探明 (Phase 8 M1 Step 2)
#
# 对照点 (PHASE8_SPEC §1.3):
#   ARDL 系数/σ̂²/IC vs statsmodels ARDL.fit 1e-10
#   ardl_select_order 最优阶 + IC 值 vs statsmodels 1e-10 (Step 3 消费, 一次生成)
#
# 语义探明 [S] (statsmodels 0.14.6 tsa/ardl, 一手源码 + 手工重建自证):
#   S1 列序: 手工设计矩阵 (冻结序, spec §2.1 Step 2) vs statsmodels params/bse
#      逐值 diff — 同时钉死列序 / trend 列值 (1-based 于估计样本) / hold_back 默认
#   S2 IC 公式: 候选式数值识别 (回归口径 −2llf + penalty·k vs AutoReg 逐观测口径)
#
# 用法: & C:\Users\Peng\.conda\envs\open-webui\python.exe verify_ardl.py
#   产出: ardl_smoke_data.csv (R 交叉用) + ardl_statsmodels_baselines.json
#         + ../../unit/timeseries/ardl_baseline.inc (数据+基准全嵌, 测试自包含)
import itertools
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "unit", "timeseries")


# ---------------------------------------------------------------------------
# 夹具 (RandomState 跨语言固定, verify_arima.py 同模式)
# ---------------------------------------------------------------------------
def gen_k1(T=250, seed=52):
    """k=1: x1 = AR(1) 0.6; y = 1.2 + 0.55·y_{t-1} + 0.5·x1_t + 0.25·x1_{t-1} + ε"""
    rng = np.random.RandomState(seed)
    burn = 100
    n = T + burn
    e_x = rng.standard_normal(n)
    e_y = rng.standard_normal(n)
    x1 = np.zeros(n)
    for t in range(1, n):
        x1[t] = 0.6 * x1[t - 1] + e_x[t]
    y = np.zeros(n)
    for t in range(1, n):
        y[t] = 1.2 + 0.55 * y[t - 1] + 0.5 * x1[t] + 0.25 * x1[t - 1] + e_y[t]
    return y[burn:], x1[burn:]


def gen_k2(T=250, seed=53):
    """k=2: x1/x2 = AR(1) 0.5/−0.4; y 双 force + AR(2)"""
    rng = np.random.RandomState(seed)
    burn = 100
    n = T + burn
    e = rng.standard_normal((n, 3))
    x1 = np.zeros(n)
    x2 = np.zeros(n)
    y = np.zeros(n)
    for t in range(1, n):
        x1[t] = 0.5 * x1[t - 1] + e[t, 0]
        x2[t] = -0.4 * x2[t - 1] + e[t, 1]
        y[t] = (0.8 + 0.4 * y[t - 1] - 0.3 * y[t - 2] + 0.6 * x1[t]
                + 0.3 * x2[t] + e[t, 2])
    return y[burn:], x1[burn:], x2[burn:]


# ---------------------------------------------------------------------------
# 手工 ARDL OLS (spec §2.1 Step 2 冻结序; 与 statsmodels 自证的参照实现)
#
# statsmodels 0.14.6 ARDL 语义实录 (本脚本 [S] 段逐条自证):
#   列序      = [trend 列 (const/trend/trend², prepend)] + [y.L1..y.Lp]
#               + [xj.L0..xj.Lq, j=1..k]  (源码 _construct_regressors)
#   trend 值  = 全样本 1-based 再按 hold_back 切片: t = hb+1 .. T
#               (非估计样本内重编号 1..n — 与 statsmodels VAR 约定不同!)
#   hold_back = 默认 max(p, max_j q_j) = _maxlag (L0 不贡献)
#   bse/vcov  = s²·(X'X)⁻¹, s² = SSR/(n_act − K), n_act = T − hold_back (真实 df)
#   ⚠️ 计账 quirk (幻觉点 A8, 精确机制源码实锤):
#      AutoReg.__init__ L262 self.nobs = _y.shape[0] 用 ARDL placeholder
#      _y = (T − hb_early, 0); hb_early 在 _check_lags 阶段 = 显式 hb 或
#      max(y_lags)=p (仅 y 滞后, exog 滞后尚未并入!). 之后
#      _construct_regressors 才把 exog 滞后并入 _maxlag 定稿 hold_back,
#      但 nobs 不回写 → res.nobs = n_acct = T − hb_early ≠ n_act (当
#      max_j q_j > p 或显式 hb 时暴露). res.sigma2/llf/aic/bic/hqic 全部
#      走 n_acct 分母; bse/vcov 走真实 df (OLS 结果固化, 不受覆盖影响)
#      → llf = −n_acct/2·(ln2π + ln(SSR/n_acct) + 1); sigma2 = SSR/n_acct;
#        aic = −2llf + 2(K+1); bic = −2llf + ln(n_acct)·(K+1);
#        hqic = −2llf + 2·lnln(n_acct)·(K+1)  (σ² 计入惩罚, 恒 K+1)
#   ⚠️ order 语义: int = range(causal, q+1) 全展开; list = 同一显式滞后列表
#      套所有变量 (非逐变量!); 逐变量阶必须 dict (源码 _format_order)
#   ⚠️ statsmodels ARDL 无 "ctt" 档 (仅 n/c/t/ct); "t" 档本库不采
# ---------------------------------------------------------------------------
def manual_ardl(y, xs, p, qs, trend, hold_back=None):
    T = len(y)
    k = len(xs)
    qs = list(qs) if isinstance(qs, (list, tuple)) else [qs] * k if k else []
    hb = hold_back if hold_back is not None else max([p] + qs)
    nobs = T - hb                      # n_act: 真实估计样本
    # A8: nobs 计账在 AutoReg.__init__ 固化 = T − hb_early, hb_early =
    # 显式 hb 或 max(y_lags)=p (仅 y 滞后) — exog 滞后不参与
    hb_early = hold_back if hold_back is not None else p
    n_acct = T - hb_early              # statsmodels 计账样本 (A8 quirk)
    cols = []
    names = []
    if trend in ("c", "ct", "ctt"):
        cols.append(np.ones(nobs))
        names.append("const")
    if trend in ("ct", "ctt"):
        tvec = np.arange(hb + 1, T + 1, dtype=float)  # 全样本 1-based 再切片
        cols.append(tvec)
        names.append("trend")
        if trend == "ctt":
            cols.append(tvec ** 2)
            names.append("trend_squared")
    for i in range(1, p + 1):  # y.L1 .. y.Lp
        cols.append(y[hb - i: T - i])
        names.append(f"y.L{i}")
    for j, (x, q) in enumerate(zip(xs, qs)):  # xj.L0 .. xj.Lq
        for l in range(0, q + 1):
            cols.append(x[hb - l: T - l])
            names.append(f"x{j + 1}.L{l}")
    X = np.column_stack(cols)
    Y = y[hb:]
    XtX = X.T @ X
    beta = np.linalg.solve(XtX, X.T @ Y)
    resid = Y - X @ beta
    ssr = float(resid @ resid)
    K = X.shape[1]
    s2_df = ssr / (nobs - K)           # bse/vcov 用 (真实 df)
    s2_acct = ssr / n_acct             # res.sigma2 用 (A8)
    XtX_inv = np.linalg.inv(XtX)
    vcov = XtX_inv * s2_df
    bse = np.sqrt(np.diag(vcov))
    llf = -n_acct / 2.0 * (np.log(2 * np.pi) + np.log(ssr / n_acct) + 1.0)
    aic = -2 * llf + 2.0 * (K + 1)
    bic = -2 * llf + np.log(n_acct) * (K + 1)
    hqic = -2 * llf + 2.0 * np.log(np.log(n_acct)) * (K + 1)
    return {
        "beta": beta, "bse": bse, "names": names, "nobs": nobs,
        "hold_back": hb, "K": K, "ssr": ssr, "s2_df": s2_df,
        "s2_acct": s2_acct, "llf": llf, "aic": aic, "bic": bic,
        "hqic": hqic, "n_acct": n_acct,
    }


def sm_fit(y, xs, p, qs, trend, hold_back=None):
    """statsmodels ARDL(...).fit() — 主基准 (pandas 输入 → params 带 .index 名)"""
    import pandas as pd
    from statsmodels.tsa.ardl import ARDL
    if isinstance(qs, (list, tuple)):
        # ⚠️ statsmodels order 语义 (源码 _format_order 实录): int = 全展开
        # range(causal, q+1); list = 同一显式滞后列表套所有变量 (非逐变量!);
        # 逐变量最大阶必须用 dict {列名: q_j} — 语义陷阱已在 [S1] 探明
        order = {f"x{j + 1}": int(q) for j, q in enumerate(qs)}
    else:
        order = qs
    endog = pd.Series(np.asarray(y), name="y")
    exog = pd.DataFrame(np.column_stack(xs), columns=[f"x{j + 1}" for j in range(len(xs))]) if xs else None
    m = ARDL(endog, p, exog, order if xs else None, trend=trend,
             causal=False, hold_back=hold_back)
    return m.fit()


def manual_uecm_from_ardl(ardl_params, ardl_vcov, p, qs, trend):
    """ARDL → UECM 系数 B 变换 (SPEC §2.3 Step 1, 不重估; C++ uecm_from_ardl 镜像)

    列序 (statsmodels UECM._construct_regressors 实录, blocks 序):
      UECM = [det(kt), π_y, π_1..π_k, ψ_1..ψ_{p-1}, γ_{1,0}..γ_{1,q_1-1}, ..., γ_{k,·}]
      ARDL = [det(kt), φ_1..φ_p, θ_{1,0}..θ_{1,q_1}, ..., θ_{k,·}]
    代数 (恒等重参数化, OLS 下与重估等价):
      π_y = Σφ_i − 1 (仿射, b0 = −1); π_j = Σ_{l=0..q_j} θ_{j,l}
      ψ_i = −Σ_{r=i+1..p} φ_r (i = 1..p−1)
      γ_{j,0} = θ_{j,0}; γ_{j,l} = −Σ_{r=l+1..q_j} θ_{j,r} (l ≥ 1)
    """
    kt = {"n": 0, "c": 1, "ct": 2, "ctt": 3}[trend]
    k = len(qs)
    n_a = len(ardl_params)
    n_u = kt + 1 + k + (p - 1 if p >= 1 else 0) + sum(qs)
    B = np.zeros((n_u, n_a))
    b0 = np.zeros(n_u)
    off = []
    cur = kt + p
    for q in qs:
        off.append(cur)
        cur += q + 1
    for i in range(kt):  # det 恒等
        B[i, i] = 1.0
    row = kt  # π_y
    for i in range(p):
        B[row, kt + i] = 1.0
    b0[row] = -1.0
    for j in range(k):  # π_j
        for l in range(qs[j] + 1):
            B[kt + 1 + j, off[j] + l] = 1.0
    for i in range(1, p):  # ψ_i
        for r in range(i + 1, p + 1):
            B[kt + 1 + k + (i - 1), kt + r - 1] = -1.0
    row0 = kt + 1 + k + (p - 1 if p >= 1 else 0)  # γ_{j,l}
    for j in range(k):
        for l in range(qs[j]):
            row = row0 + sum(qs[:j]) + l
            if l == 0:
                B[row, off[j]] = 1.0
            else:
                for r in range(l + 1, qs[j] + 1):
                    B[row, off[j] + r] = -1.0
    return B @ ardl_params + b0, B @ ardl_vcov @ B.T


# ---------------------------------------------------------------------------
# [S] 语义探明 (首次运行审阅 → 结论冻结于本文件头)
# ---------------------------------------------------------------------------
FIXTURES = {}
SPECS = [
    # name: (fx, p, qs, trend, hold_back)
    ("c_p2q1",  "k1", 2, [1], "c", None),
    ("n_p2q1",  "k1", 2, [1], "n", None),
    ("ct_p1q2", "k1", 1, [2], "ct", None),
    ("ctt_p1q1", "k1", 1, [1], "ctt", None),
    ("c_k2",    "k2", 2, [1, 2], "c", None),
    ("c_p0q1",  "k1", 0, [1], "c", None),
    ("c_p2q0",  "k1", 2, [0], "c", None),
    ("c_hb5",   "k1", 2, [1], "c", 5),
    ("ar2",     "k1", 2, [], "c", None),  # k=0 纯 AR (无 exog)
]

# UECM/bounds 对拍矩阵 (Step 4): trend 由 case 映射 (statsmodels bounds_test
# 内部按 case 重构 UECM 的 trend; 对拍前提 = ARDL trend 与之一致)
CASE_TREND = {1: "n", 2: "c", 3: "c", 4: "ct", 5: "ct"}
UECM_SPECS = [
    # (name, fx, p, qs, case)
    ("uecm_c1_k1", "k1", 2, [1], 1),
    ("uecm_c2_k1", "k1", 2, [1], 2),
    ("uecm_c3_k1", "k1", 2, [1], 3),
    ("uecm_c4_k1", "k1", 2, [1], 4),
    ("uecm_c5_k1", "k1", 2, [1], 5),
    ("uecm_c3_k2", "k2", 2, [1, 2], 3),
    ("uecm_c5_k2", "k2", 2, [1, 2], 5),
]


def section_probe():
    print("=" * 72)
    print("[S1] 手工设计矩阵 (冻结序) vs statsmodels params/bse 逐值 diff")
    print("     (ctt 无 statsmodels 对照 — ARDL 无该档, 跳过)")
    for name, fx, p, qs, trend, hb in SPECS:
        y, xs = FIXTURES[fx]
        man = manual_ardl(y, xs, p, qs, trend, hb)
        if trend == "ctt":
            print(f"  {name:8s} trend=ctt p={p} q={qs}: n_act={man['nobs']} "
                  f"n_acct={man['n_acct']} K={man['K']} (manual only)")
            continue
        res = sm_fit(y, xs, p, qs, trend, hb)
        sm_names = list(res.params.index)
        dp = float(np.max(np.abs(np.asarray(res.params) - man["beta"])))
        db = float(np.max(np.abs(np.asarray(res.bse) - man["bse"])))
        ok_names = sm_names == man["names"]
        print(f"  {name:8s} trend={trend:3s} p={p} q={qs}: "
              f"nobs(sm)={res.nobs} n_act={man['nobs']} n_acct={man['n_acct']} "
              f"dparam={dp:.2e} dbse={db:.2e} names_match={ok_names}")
        if not ok_names:
            print(f"    sm : {sm_names}")
            print(f"    man: {man['names']}")

    print("=" * 72)
    print("[S2] 镜像公式自证 (n_acct = T − hb_early 计账: llf/sigma2/aic/bic/hqic)")
    for name, fx, p, qs, trend, hb in SPECS:
        if trend == "ctt":
            continue
        y, xs = FIXTURES[fx]
        man = manual_ardl(y, xs, p, qs, trend, hb)
        res = sm_fit(y, xs, p, qs, trend, hb)
        checks = {
            "llf": abs(res.llf - man["llf"]),
            "aic": abs(res.aic - man["aic"]),
            "bic": abs(res.bic - man["bic"]),
            "hqic": abs(res.hqic - man["hqic"]),
            "sigma2": abs(res.sigma2 - man["s2_acct"]),
            "nobs": abs(res.nobs - man["n_acct"]),
        }
        worst = max(checks.values())
        print(f"  {name:8s}: " + " ".join(f"{k}={v:.1e}" for k, v in checks.items())
              + ("  ✓" if worst < 1e-10 else "  ✗✗✗"))

    print("=" * 72)
    print("[S3] hold_back 默认 = max(p, max_j q_j) 抽查 (显式 hb=5 用例在 S1 c_hb5)")


# ---------------------------------------------------------------------------
# [V] 基准生成 (JSON + .inc)
# ---------------------------------------------------------------------------
def section_verify():
    from statsmodels.tsa.ardl import ardl_select_order

    baselines = {}
    for name, fx, p, qs, trend, hb in SPECS:
        y, xs = FIXTURES[fx]
        man = manual_ardl(y, xs, p, qs, trend, hb)
        if trend == "ctt":
            # 无 statsmodels 对照: 基准 = 手工 numpy OLS (C++ 侧只验证
            # 数值管线与列序; ctt 的语义正确性由构造式 [1,t,t²] 直接保证)
            entry = {
                "params": [float(v) for v in man["beta"]],
                "bse": [float(v) for v in man["bse"]],
                "llf": float(man["llf"]),
                "aic": float(man["aic"]),
                "bic": float(man["bic"]),
                "hqic": float(man["hqic"]),
                "nobs": int(man["n_acct"]),
                "hold_back": int(man["hold_back"]),
                "K": int(man["K"]),
                "sigma2_df": float(man["s2_df"]),
                "sigma2_acct": float(man["s2_acct"]),
                "manual_max_param_diff": 0.0,
                "no_sm_counterpart": True,
                "spec": {"p": p, "q": qs, "trend": trend, "hold_back": hb},
            }
            baselines[name] = entry
            print(f"{name}: (manual only, no statsmodels ctt) "
                  f"llf={entry['llf']:.6f}")
            continue
        res = sm_fit(y, xs, p, qs, trend, hb)
        entry = {
            "params": [float(v) for v in res.params],
            "bse": [float(v) for v in res.bse],
            "llf": float(res.llf),
            "aic": float(res.aic),
            "bic": float(res.bic),
            "hqic": float(res.hqic),
            "nobs": int(res.nobs),
            "hold_back": int(man["hold_back"]),
            "K": int(man["K"]),
            "sigma2_df": float(man["s2_df"]),
            "sigma2_acct": float(man["s2_acct"]),
            "manual_max_param_diff": float(np.max(np.abs(
                np.asarray(res.params) - man["beta"]))),
            "spec": {"p": p, "q": qs, "trend": trend, "hold_back": hb},
        }
        baselines[name] = entry
        print(f"{name}: dparam={entry['manual_max_param_diff']:.2e} "
              f"llf={entry['llf']:.6f} aic={entry['aic']:.6f}")

    # ardl_select_order 基准 (Step 3 消费)
    #
    # ⚠️ 选阶 IC 计账 (幻觉点 A8 连带, 源码实锤): 网格内 compute_ics 用
    #    n_act = T − hb_grid (统一 hb = max(maxlag, maxorder), 可比性),
    #    lstsq OLS + Frisch-Waugh 对 always(deterministic) 正交化;
    #    而 sel.model refit 后 fit().aic 走 A8 计账 n_acct = T − p_best.
    #    当 maxorder > p_best 时两口径 IC 值不同 → 双值入基准供 Step 3 断言
    sel_specs = [
        ("sel_k1", "k1", 3, 2),
        ("sel_k2", "k2", 2, 1),
    ]
    for name, fx, maxlag, maxorder in sel_specs:
        y, xs = FIXTURES[fx]
        sel = ardl_select_order(y, maxlag, np.column_stack(xs), maxorder,
                                trend="c")
        hb_grid = max(maxlag, maxorder)
        n_act_grid = len(y) - hb_grid
        entry = {}
        for ic in ("aic", "bic", "hqic"):
            s = getattr(sel, ic)  # pd.Series: index=IC 值升序, value=(p, {var: q})
            p_best, dl = s.iloc[0]
            order = [int(p_best)] + [int(v) if v is not None else 0
                                     for v in dl.values()]
            assert tuple(order) == tuple(int(v) for v in sel.model.ardl_order)
            res_best = sel.model.fit()  # refit: A8 计账 nobs = T − hb_early
            entry[ic] = {
                "order": order,
                "ic_value": float(s.index[0]),       # 网格口径 (n_act 统一 hb)
                "ic_value_refit": float(getattr(res_best, ic)),  # A8 口径
                "nobs_refit": int(res_best.nobs),
            }
        baselines[name] = {"grid": {"maxlag": maxlag, "maxorder": maxorder,
                                    "hold_back": hb_grid, "n_act": n_act_grid},
                           "sel": entry}
        print(f"{name}: bic→{entry['bic']['order']} "
              f"ic={entry['bic']['ic_value']:.6f} "
              f"(refit {entry['bic']['ic_value_refit']:.6f}) "
              f"n_act={n_act_grid}")

    with open(os.path.join(HERE, "ardl_statsmodels_baselines.json"), "w") as f:
        json.dump(baselines, f, indent=1)
    print("written: ardl_statsmodels_baselines.json")
    return baselines


# ---------------------------------------------------------------------------
# [S4] UECM 重参数化 + bounds F/t 基准 (Step 4 消费)
#
# 语义实录 (statsmodels 0.14.6 + R ARDL 0.2.5, probe_r_ardl.R/gen_r_ardl_bounds.R):
#   - UECM.from_ardl: 重构模型实例再 fit (数学上 = B 变换, OLS 等价 — [S4]
#     手工 B 变换 vs UECM.fit() params/bse 逐值 diff 自证)
#   - bounds_test: 内部按 case 重构 UECM (trend = case 映射, 无显式 hb);
#     Wald = (Rβ)'(RVR')⁻¹(Rβ)/q, 约束行 = statsmodels rest 索引实录
#   - PSS10: statsmodels 查表键 k = nvar = 1+k_PSS (偏移); 本库 CV/p 用
#     表键直查 (k_PSS) — statsmodels bounds_test 返回的 CV/p 不对照
#   - t 检验: statsmodels 无; 构造式 t = π_y/bse + R ARDL 交叉 (A4)
#   - R ARDL: case-trend 不相容 throw (parse_case 实录, 本库镜像);
#     F/t 对 trend 列线性变换不变 (Wald 不变性 → time() 尺度无关)
# ---------------------------------------------------------------------------
def section_uecm():
    from statsmodels.tsa.ardl import UECM
    from statsmodels.tsa.ardl.model import _pss_pvalue
    import statsmodels.tsa.ardl.pss_critical_values as pcv

    r_bounds = {}
    with open(os.path.join(HERE, "r_ardl_bounds.csv")) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            nm, fv, tv = ln.split(",")
            r_bounds[nm] = (float(fv), None if tv == "NA" else float(tv))

    print("=" * 72)
    print("[S4] UECM B 变换 vs statsmodels UECM.from_ardl().fit() + bounds F/t")
    baselines = {}
    for name, fx, p, qs, case in UECM_SPECS:
        trend = CASE_TREND[case]
        y, xs = FIXTURES[fx]
        res_ardl = sm_fit(y, xs, p, qs, trend)
        par_u, vcov_u = manual_uecm_from_ardl(
            np.asarray(res_ardl.params), np.asarray(res_ardl.cov_params()),
            p, qs, trend)
        bse_u = np.sqrt(np.diag(vcov_u))
        res_uecm = UECM.from_ardl(res_ardl.model).fit()
        dp = float(np.max(np.abs(np.asarray(res_uecm.params) - par_u)))
        db = float(np.max(np.abs(np.asarray(res_uecm.bse) - bse_u)))

        # F: statsmodels bounds_test + 手工 Wald 镜像 (C++ bounds_f_test 的
        # 直接参照; rest 索引 = statsmodels bounds_test L2399-2408 实录)
        bt = res_uecm.bounds_test(case=case)
        cov_u = np.asarray(res_uecm.cov_params())
        params_u = np.asarray(res_uecm.params)
        nvar = 1 + len(qs)
        rest = {1: np.arange(nvar), 2: np.arange(nvar + 1),
                3: np.arange(1, nvar + 1), 4: np.arange(1, nvar + 2),
                5: np.arange(2, nvar + 2)}[case]
        R = np.zeros((len(rest), cov_u.shape[1]))
        for i, loc in enumerate(rest):
            R[i, loc] = 1.0
        coef = R @ params_u
        wald = float(coef.T @ np.linalg.solve(R @ cov_u @ R.T, coef) / len(rest))

        k_pss = len(qs)
        cv_lo = [float(v) for v in pcv.crit_vals[(k_pss, case, False)]]
        cv_up = [float(v) for v in pcv.crit_vals[(k_pss, case, True)]]
        p_lo = float(_pss_pvalue(bt.stat, k_pss, case, False))
        p_up = float(_pss_pvalue(bt.stat, k_pss, case, True))

        kt = {"n": 0, "c": 1, "ct": 2}[trend]
        t_sm = float(params_u[kt] / np.asarray(res_uecm.bse)[kt])
        F_R, t_R = r_bounds[name]

        entry = {
            "params": [float(v) for v in par_u],
            "bse": [float(v) for v in bse_u],
            "nobs": int(len(res_uecm.resid)), "K": int(len(par_u)),
            "F": float(bt.stat), "F_manual": wald, "F_R": F_R,
            "t_sm": t_sm, "t_R": t_R,
            "p_lo": p_lo, "p_up": p_up,
            "cv_lo": cv_lo, "cv_up": cv_up,
            "dparam": dp, "dbse": db,
            "spec": {"p": p, "q": qs, "trend": trend, "case": case},
        }
        baselines[name] = entry
        dFR = abs(entry["F"] - F_R) if F_R is not None else float("nan")
        dFRm = abs(entry["F"] - wald)
        print(f"  {name}: dparam={dp:.2e} dbse={db:.2e} "
              f"F={entry['F']:.6f} (manual {wald:.6f}, R {F_R:.6f}) "
              f"dF_sm_R={dFR:.2e} dF_sm_manual={dFRm:.2e}")
        assert dp < 1e-10 and db < 1e-10, "B 变换与 statsmodels UECM 不一致"
        assert dFRm < 1e-10, "手工 Wald 与 bounds_test 不一致"
        if F_R is not None:
            assert dFR < 1e-8, f"F 与 R 偏离 {dFR:.2e}"

    # danish 案例 (PSS10 回归断言数据): k_PSS=3, p=2, q=1, trend="c", case 3
    from statsmodels.datasets import danish_data
    dd = danish_data.load_pandas().data
    res_d = sm_fit(dd.lrm.values, [dd.lry.values, dd.ibo.values,
                                   dd.ide.values], 2, [1, 1, 1], "c")
    res_ud = UECM.from_ardl(res_d.model).fit()
    bt_d = res_ud.bounds_test(case=3)
    cv_lib = [float(v) for v in pcv.crit_vals[(3, 3, False)]] + \
             [float(v) for v in pcv.crit_vals[(3, 3, True)]]
    cv_sm = [float(v) for v in bt_d.crit_vals["lower"].values] + \
            [float(v) for v in bt_d.crit_vals["upper"].values]
    assert cv_lib != cv_sm, "danish: 键 (3,3) 与 (4,3) 意外重合 (PSS10 失效?)"
    danish = {
        "lrm": [float(v) for v in dd.lrm.values],
        "lry": [float(v) for v in dd.lry.values],
        "ibo": [float(v) for v in dd.ibo.values],
        "ide": [float(v) for v in dd.ide.values],
        "F": float(bt_d.stat), "nobs": int(len(res_ud.resid)),
        "cv_lib": cv_lib, "cv_sm": cv_sm,
    }
    baselines["danish"] = danish
    print(f"  danish: F={danish['F']:.6f} "
          f"lib_cv_5%=[{cv_lib[1]:.4f},{cv_lib[5]:.4f}] "
          f"sm_cv_5%=[{cv_sm[1]:.4f},{cv_sm[5]:.4f}] (PSS10 偏移 ✓)")

    print("written: (uecm/danish 并入 baselines, main 统一落 JSON)")
    return baselines


def fmt(v):
    return repr(float(v))


def arr(name, vals):
    body = ",\n    ".join(fmt(v) for v in vals)
    return f"inline constexpr double {name}[] = {{\n    {body},\n}};\n"


def write_inc(baselines):
    y1, (x1,) = FIXTURES["k1"]
    y2, (x21, x22) = FIXTURES["k2"]
    parts = [
        "// ============================================================================",
        "// ardl_baseline.inc - ARDL 黄金基准 (自动生成 by verify_ardl.py, 勿手改)",
        "//",
        "// 数据: RandomState 夹具 (跨语言固定, repr(float) 全精度); 基准:",
        "// statsmodels 0.14.6 ARDL.fit (params/bse/llf/IC); SCALARS = [llf,",
        "// aic, bic, hqic, sigma2_df, sigma2_acct, hold_back, nobs, K]",
        "// ============================================================================",
        "namespace cpphub::v1::ardl_baseline {",
        "",
        arr("DATA_K1_Y", y1),
        arr("DATA_K1_X1", x1),
        arr("DATA_K2_Y", y2),
        arr("DATA_K2_X1", x21),
        arr("DATA_K2_X2", x22),
    ]
    for name, fx, p, qs, trend, hb in SPECS:
        b = baselines[name]
        parts.append(arr(f"PARAMS_{name.upper()}", b["params"]))
        parts.append(arr(f"BSE_{name.upper()}", b["bse"]))
        parts.append(arr(f"SCALARS_{name.upper()}",
                         [b["llf"], b["aic"], b["bic"], b["hqic"],
                          b["sigma2_df"], b["sigma2_acct"],
                          b["hold_back"], b["nobs"], b["K"]]))
    for name, fx, maxlag, maxorder in [("sel_k1", "k1", 3, 2),
                                       ("sel_k2", "k2", 2, 1)]:
        b = baselines[name]
        s = b["sel"]
        # [maxlag, maxorder, hb_grid, n_act_grid] + 每 IC: order + [ic, ic_refit, nobs_refit]
        parts.append(arr(f"{name.upper()}_GRID",
                         [maxlag, maxorder, b["grid"]["hold_back"],
                          b["grid"]["n_act"]]))
        for ic in ("aic", "bic", "hqic"):
            e = s[ic]
            parts.append(arr(f"{name.upper()}_{ic.upper()}",
                             e["order"] + [e["ic_value"],
                                           e["ic_value_refit"],
                                           e["nobs_refit"]]))
    # UECM/bounds 基准 (Step 4):
    #   UECM_<N>_PARAMS/BSE   = B 变换系数/标准误 (vs statsmodels UECM 1e-10)
    #   BOUNDS_<N>_STAT       = [F, F_R, p_lo, p_up, nobs, K]
    #   BOUNDS_<N>_CV         = [I0 四档, I1 四档] (本库 k_PSS 键直查)
    #   BOUNDS_<N>_T          = [t_R, t_sm] (仅 case 1/3/5; statsmodels 无 t)
    #   DANISH_*              = PSS10 k 偏移回归断言 (k_PSS=3, case 3;
    #                           cv_lib = 键 (3,3) ≠ cv_sm = bounds_test 偏移键 (4,3))
    for name, fx, p, qs, case in UECM_SPECS:
        b = baselines[name]
        n = name.upper().replace("UECM_", "", 1)  # uecm_c1_k1 → C1_K1
        parts.append(arr(f"UECM_{n}_PARAMS", b["params"]))
        parts.append(arr(f"UECM_{n}_BSE", b["bse"]))
        parts.append(arr(f"BOUNDS_{n}_STAT",
                         [b["F"], b["F_R"], b["p_lo"], b["p_up"],
                          b["nobs"], b["K"]]))
        parts.append(arr(f"BOUNDS_{n}_CV", b["cv_lo"] + b["cv_up"]))
        if b["t_R"] is not None:
            parts.append(arr(f"BOUNDS_{n}_T", [b["t_R"], b["t_sm"]]))
    d = baselines["danish"]
    parts.append(arr("DANISH_LRM", d["lrm"]))
    parts.append(arr("DANISH_LRY", d["lry"]))
    parts.append(arr("DANISH_IBO", d["ibo"]))
    parts.append(arr("DANISH_IDE", d["ide"]))
    parts.append(arr("DANISH_BOUNDS", [d["F"], d["nobs"]]))
    parts.append(arr("DANISH_CV_LIB", d["cv_lib"]))
    parts.append(arr("DANISH_CV_SM", d["cv_sm"]))
    parts.append("}  // namespace cpphub::v1::ardl_baseline")
    path = os.path.join(OUT, "ardl_baseline.inc")
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(parts) + "\n")
    print(f"written: {path}")


def write_csv():
    y1, (x1,) = FIXTURES["k1"]
    y2, (x21, x22) = FIXTURES["k2"]
    cols = {"k1_y": y1, "k1_x1": x1, "k2_y": y2, "k2_x1": x21, "k2_x2": x22}
    keys = list(cols)
    lens = {len(v) for v in cols.values()}
    assert len(lens) == 1, f"column length mismatch: {lens}"
    with open(os.path.join(HERE, "ardl_smoke_data.csv"), "w", newline="\n") as f:
        for i in range(lens.pop()):
            f.write(",".join("%.17g" % cols[k][i] for k in keys) + "\n")
        f.write("#columns=" + ",".join(keys) + "\n")
    print("written: ardl_smoke_data.csv")


def main():
    global FIXTURES
    y1, x1 = gen_k1()
    y2, x21, x22 = gen_k2()
    FIXTURES = {"k1": (y1, [x1]), "k2": (y2, [x21, x22])}
    section_probe()
    baselines = section_verify()
    baselines.update(section_uecm())
    with open(os.path.join(HERE, "ardl_statsmodels_baselines.json"), "w") as f:
        json.dump(baselines, f, indent=1)
    write_inc(baselines)
    write_csv()


if __name__ == "__main__":
    main()
