# -*- coding: utf-8 -*-
# verify_var.py - VAR 基准生成: statsmodels VAR (Phase 7C M2)
#
# 对照点 (spec §5.1/§5.2):
#   系数/IC 五式/IRF/FEVD vs statsmodels var_model 1e-10
#   select_order 同样本 offset 轨迹 (V5)
#   Cholesky 下三角 P = np.linalg.cholesky(sigma_u) (V2)
#
# 一手源码语义落档 (statsmodels 0.14.4, var_model.py, 与 spec 冻结 0.14.6 的
# L2281-2303 逐字一致; 版本差异记录):
#   - info_criteria (L2281-2303):
#       free_params = k_ar·K² + K·k_exog;  ld = logdet(Σ_mle)
#       aic = ld + 2/nobs·fp;  bic = ld + ln(nobs)/nobs·fp
#       hqic = ld + 2·lnln(nobs)/nobs·fp
#       fpe = ((nobs + df_model)/df_resid)^K · exp(ld)   (V6: 指数 K,
#         df_model = k_ar·K + k_exog 单方程, n* 含确定性项)
#   - sigma_u_mle = sigma_u · df_resid/nobs = SSR/T  (V4: ÷T)
#     ⚠️ sigma_u 本体 = SSR/(nobs − K·k_ar − k_trend) (df 修正)
#   - _chol_sigma_u = np.linalg.cholesky(sigma_u)  ← 用 df 修正版 Σ!
#     (V2: numpy cholesky = 下三角; IRF/FEVD 正交化基准用此矩阵)
#   - select_order (L778-830): offset = maxlags − p (同样本, V5),
#     p_min = 0 (trend ≠ "n" 时纯截距模型参与), maxlags 默认
#     round(12·(T/100)^{1/4}); selected = argmin + p_min
#   - var_loglike (L305-337):
#       ll = −(nobs·K/2)·ln(2π) − (nobs/2)·(ld + K)
#   - FEVD (L2370-2399): fevd[i] = cumsum(Ψ²)[i] / mse_diag[i] 逐行;
#     orth_irfs = Φ_h·P (P=下三角 chol); 与行和归一化恒等 (正交化情形)
#
# 用法: python verify_var.py  (生成 var_smoke_data.csv + var_statsmodels_baselines.json)
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# DGP: 稳定 VAR(2), K=3, T=250 (非对角交互 + 非对角协方差, 保证
# Granger/溢出/Cholesky 次序敏感非平凡)
A1_TRUE = np.array([[0.5, 0.1, 0.0],
                    [0.0, 0.4, 0.2],
                    [0.1, 0.0, 0.3]])
A2_TRUE = np.array([[-0.2, 0.0, 0.1],
                    [0.05, -0.1, 0.0],
                    [0.0, 0.1, -0.05]])
C_TRUE = np.array([0.2, -0.1, 0.05])
SIGMA_TRUE = np.array([[1.0, 0.4, 0.2],
                       [0.4, 1.5, 0.3],
                       [0.2, 0.3, 0.8]])


def gen_var2(T=250, seed=42, burn=200):
    """稳定 VAR(2) 模拟 (RandomState 跨固定)"""
    rng = np.random.RandomState(seed)
    P = np.linalg.cholesky(SIGMA_TRUE)
    n = T + burn
    y = np.zeros((n, 3))
    eps = rng.standard_normal((n, 3)) @ P.T
    for t in range(2, n):
        y[t] = C_TRUE + A1_TRUE @ y[t - 1] + A2_TRUE @ y[t - 2] + eps[t]
    return y[burn:]


def main():
    from statsmodels.tsa.api import VAR

    y = gen_var2()
    np.savetxt(os.path.join(HERE, "var_smoke_data.csv"), y,
               delimiter=",", fmt="%.17g",
               header="y1,y2,y3", comments="")

    out = {"statsmodels_version": None, "fixtures": {}}
    import statsmodels
    out["statsmodels_version"] = statsmodels.__version__

    # --- V1: 主拟合 p=2 trend=c ---
    model = VAR(y)
    fit = model.fit(2, trend="c")
    f1 = {
        "params": fit.params.tolist(),            # K×(K·p+k_trend)
        "sigma_u": fit.sigma_u.tolist(),          # df 修正版
        "sigma_u_mle": fit.sigma_u_mle.tolist(),  # SSR/T (V4)
        "nobs": int(fit.nobs),
        "df_model": int(fit.df_model),
        "df_resid": int(fit.df_resid),
        "loglik": float(fit.llf) if hasattr(fit, "llf") else None,
        "aic": float(fit.aic), "bic": float(fit.bic),
        "hqic": float(fit.hqic), "fpe": float(fit.fpe),
        "detomega": float(fit.detomega),
        "roots": np.sort(np.abs(fit.roots)).tolist(),  # 伴随矩阵特征值模
        "is_stable": bool(fit.is_stable(verbose=False)),
    }
    # llf: statsmodels 无直接 llf 属性? var_loglike 手算
    from statsmodels.tsa.vector_ar.var_model import var_loglike
    f1["loglik"] = float(var_loglike(fit.resid, fit.sigma_u_mle, fit.nobs))
    out["fixtures"]["var2_c"] = f1

    # --- V2: Cholesky 下三角 (V2, 用 df 修正版 Σ) ---
    Pchol = np.linalg.cholesky(fit.sigma_u)
    out["fixtures"]["chol_sigma_u"] = Pchol.tolist()

    # --- V3: trend=n 与 trend=ct 单点 ---
    for trend, key in (("n", "var2_n"), ("ct", "var2_ct")):
        ft = model.fit(2, trend=trend)
        kt = ft.k_trend
        out["fixtures"][key] = {
            "params": ft.params.tolist(),
            "sigma_u_mle": ft.sigma_u_mle.tolist(),
            "aic": float(ft.aic), "bic": float(ft.bic),
            "hqic": float(ft.hqic), "fpe": float(ft.fpe),
            "k_trend": int(kt),
            "loglik": float(var_loglike(ft.resid, ft.sigma_u_mle, ft.nobs)),
            "nobs": int(ft.nobs),
        }

    # --- V4: select_order 同样本轨迹 (maxlags=4, V5) ---
    sel = model.select_order(4, trend="c")
    out["fixtures"]["select_order_max4"] = {
        "ics": {k: [float(x) for x in v] for k, v in sel.ics.items()},
        "selected": {k: int(v) for k, v in sel.selected_orders.items()},
    }

    # --- V5: IRF 正交化 (V3, horizon 10) ---
    irfs = fit.orth_ma_rep(maxn=10)   # Ψ_0..Ψ_10 (11 个)
    out["fixtures"]["orth_irf_h10"] = irfs.tolist()
    phi = fit.ma_rep(maxn=10)         # 未正交化 Φ_0..Φ_10
    out["fixtures"]["phi_ma_rep_h10"] = phi.tolist()

    # --- V6: FEVD Cholesky 轨 (V7 行和=1) ---
    fevd = fit.fevd(10)
    out["fixtures"]["fevd_orth_periods10"] = fevd.decomp.tolist()

    # --- V7: 变量重排敏感性 (V2 注入/reorder) ---
    fit_r = fit.reorder([2, 0, 1])
    out["fixtures"]["reorder_201"] = {
        "params": fit_r.params.tolist(),
        "fevd": fit_r.fevd(10).decomp.tolist(),
    }

    # --- V8: p=0 纯截距 (select_order p_min=0 成员) ---
    f0 = model.fit(0, trend="c")
    out["fixtures"]["var0_c"] = {
        "params": f0.params.tolist(),
        "aic": float(f0.aic), "bic": float(f0.bic),
        "hqic": float(f0.hqic), "fpe": float(f0.fpe),
        "sigma_u_mle": f0.sigma_u_mle.tolist(),
    }

    path = os.path.join(HERE, "var_statsmodels_baselines.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print("written:", path)
    print("select_order selected:", sel.selected_orders)
    print("aic trace:", [round(x, 4) for x in sel.ics["aic"]])


if __name__ == "__main__":
    main()
