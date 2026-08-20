# -*- coding: utf-8 -*-
# verify_granger.py - Granger 因果检验基准生成: statsmodels grangercausalitytests
#                    (Phase 7C M1)
#
# 对照点 (spec §3.2, 容差表):
#   4 统计量 (ssr_ftest / ssr_chi2test / lrtest / params_ftest) vs statsmodels
#   1e-10 (GR1 df 公式 / GR6 方向 — 显式 (cause, effect) 形参复现 SM 方向)
#   TY Wald vs statsmodels OLS f_test (增广回归) 1e-8 (GR2 df=k / GR4 增广阶不进约束)
#   HAC-Wald vs 显式 NW 三明治 (cov_hac_simple 无小样本修正) 1e-8 (GR5)
#
# statsmodels 语义 (一手源码 tsa/stattools.py L1473-1699, 0.14.4):
#   输入 x = [effect, cause] (第二列 cause 第一列 — 方向陷阱 GR6 的来源)
#   dta = lagmat2ds: [effect_t, effect_{t−1..t−p}, cause_{t−1..t−p}]
#     (lagmat trim="both" 实测 2026-08-19: 行 t = p+1..T, nobs = T − p)
#   受约束:  effect ~ [own lags, const]           (SSR_r)
#   无约束:  effect ~ [own lags, cause lags, const] (SSR_u)
#   ssr_ftest:   F = (SSR_r − SSR_u)/p·df_u/SSR_u,  df_u = nobs − (2p+1)
#   ssr_chi2:    χ² = nobs·(SSR_r − SSR_u)/SSR_u
#   lrtest:      LR = nobs·(log SSR_r − log SSR_u)
#   params_ftest: Wald F on cause-lag 子块 (R V R', V = σ̂²(X'X)⁻¹) — 与
#                ssr_ftest 数学等价, 数值独立路径 (两列分录断言)
#
# TY (Toda-Yamamoto, 决策 5): 增广回归 p* = k + d_max 阶 (GR4: 全部估计),
#   statsmodels f_test 约束 cause 前 k 阶滞后 (增广阶不进约束矩阵),
#   Wald χ² = k·F (σ̂² = SSR/(n−k_full) df 修正路径, 探测实测等价 8e-14),
#   p 值 χ²(k) (GR2: df = k, 非 k+d_max)
#
# HAC-Wald (决策 6): NW Bartlett 三明治 (statsmodels cov_hac_simple,
#   use_correction=False — get_robustcov_results('HAC') 默认约定, 探测实测
#   fvalue = W/p 严格一致 2.5e-14): W = b'(R V R')⁻¹b, V = (X'X)⁻¹Ω(X'X)⁻¹,
#   Ω = Σ_{l=0}^{L} w_l(Ω_l + Ω_l'), w_l = 1 − l/(L+1), p 值 χ²(k) (GR5)
#
# 输出: granger_statsmodels_baselines.txt
#   gr_<tag>_p<p>   10 值: [f, f_p, params_f, params_f_p, chi2, chi2_p,
#                           lr, lr_p, df_denom, df_num]
#   grty_<tag>_p<k>_d<d>  6 值: [wald, chi2_p, sm_f, sm_f_p, df_num, df_denom]
#   grhac_<tag>_p<p>_l<L> 6 值: [wald, chi2_p, sm_f, sm_f_p, df_num, df_denom]

import os

import numpy as np
from scipy import stats as sps
from statsmodels.regression.linear_model import OLS
from statsmodels.tsa.stattools import grangercausalitytests

HERE = os.path.dirname(os.path.abspath(__file__))


def stats_row(effect, cause, p):
    """调用 statsmodels (输入 [effect, cause]: 第二列 cause 第一列)。"""
    res = grangercausalitytests(np.column_stack([effect, cause]), [p],
                                verbose=False)
    r = res[p][0]
    f, f_p, df_denom, df_num = r["ssr_ftest"]
    pf, pf_p, pf_dfd, pf_dfn = r["params_ftest"]
    c2, c2_p, c2_df = r["ssr_chi2test"]
    lr, lr_p, lr_df = r["lrtest"]
    return [float(f), float(f_p), float(pf), float(pf_p),
            float(c2), float(c2_p), float(lr), float(lr_p),
            float(df_denom), float(df_num)]


def design(effect, cause, p):
    """lagmat2ds 语义: 行 t = p..T−1 (0-based), nobs = T − p。

    X = [own lags 1..p, cause lags 1..p, const] (const 末列, prepend=False)
    """
    T = len(effect)
    t = np.arange(p, T)
    lhs = effect[t]
    own = np.column_stack([effect[t - j] for j in range(1, p + 1)])
    ext = np.column_stack([cause[t - j] for j in range(1, p + 1)])
    return lhs, np.hstack([own, ext, np.ones((len(t), 1))])


def ty_row(effect, cause, k, d_max):
    """TY 增广 Wald: 回归 p* = k+d 阶, 约束 cause 前 k 阶 (GR4), df = k (GR2)。"""
    lhs, X = design(effect, cause, k + d_max)
    n, kc = X.shape
    res = OLS(lhs, X).fit()
    R = np.zeros((k, kc))
    for i in range(k):
        R[i, (k + d_max) + i] = 1.0
    ft = res.f_test(R)
    f = float(ft.fvalue)
    wald = k * f  # Wald χ² (df 修正 σ² 路径, 探测实测 = k·F 精确)
    chi2_p = float(sps.chi2.sf(wald, k))
    return [wald, chi2_p, f, float(ft.pvalue),
            k, n - kc]


def hac_wald(effect, cause, p, L):
    """NW-Bartlett 三明治 Wald (cov_hac_simple, use_correction=False)。"""
    lhs, X = design(effect, cause, p)
    n, kc = X.shape
    res = OLS(lhs, X).fit()
    b = res.params
    XtXi = np.linalg.inv(X.T @ X)
    xu = X * res.resid[:, None]
    Om = xu.T @ xu
    for lag in range(1, L + 1):
        s = xu[lag:].T @ xu[:-lag]
        Om += (1.0 - lag / (L + 1.0)) * (s + s.T)
    V = XtXi @ Om @ XtXi
    lo, hi = p, 2 * p  # cause 滞后系数块
    W = float(b[lo:hi] @ np.linalg.inv(V[lo:hi, lo:hi]) @ b[lo:hi])
    # statsmodels 交叉 (use_t=True → F 约定): fvalue = W/p, p 值 F(p, n−kc)
    rh = res.get_robustcov_results(cov_type="HAC", maxlags=L)
    R = np.zeros((p, kc))
    for i in range(p):
        R[i, p + i] = 1.0
    ft = rh.f_test(R)
    return W, float(ft.fvalue), float(ft.pvalue)


def main():
    d = np.loadtxt(os.path.join(HERE, "granger_smoke_data.csv"),
                   delimiter=",", skiprows=1)
    x, y, z = d[:, 0], d[:, 1], d[:, 2]
    ci = np.loadtxt(os.path.join(HERE, "coint_smoke_data.csv"),
                    delimiter=",", skiprows=1)
    y1, y2 = ci[:, 0], ci[:, 1]

    # (tag, effect, cause, lags): 真因果/反向/零假设/I(1) 水平 (GR7 场景)
    cases = [
        ("fwd", y, x, [1, 2, 4]),   # x → y 真因果 (DGP 一阶)
        ("rev", x, y, [1, 2]),      # y ↛ x 反向 (GR6 方向)
        ("null", y, z, [1, 2]),     # z ↛ y 零假设
        ("i1", y1, y2, [1, 2]),     # I(1) 协整对水平值 (GR7, coint 夹具)
    ]

    # TY: k 阶约束 + d_max 增广 (I(1) 惯例 d=1; fwd 含 d=2 覆盖多步增广)
    ty_cases = [
        ("fwd", y, x, [(1, 1), (2, 1), (2, 2)]),
        ("rev", x, y, [(1, 1)]),
        ("i1", y1, y2, [(1, 1), (2, 1)]),
    ]

    # HAC: 显式带宽 L (3 = 惯例; 0 = White 退化; 4 = NW 规则默认
    #   floor(4·(249/100)^(2/9)) = 4, 锁定 hac_bandwidth=0 的默认行为)
    hac_cases = [
        ("fwd", y, x, [(1, 3), (2, 3), (1, 0), (1, 4)]),
        ("rev", x, y, [(1, 3)]),
        ("null", y, z, [(1, 3)]),
        ("i1", y1, y2, [(1, 3)]),
    ]

    out = os.path.join(HERE, "granger_statsmodels_baselines.txt")
    with open(out, "w") as fh:
        for tag, eff, cau, lags in cases:
            for p in lags:
                row = stats_row(eff, cau, p)
                fh.write("gr_%s_p%d\t%s\n"
                         % (tag, p, " ".join("%.17g" % v for v in row)))
        for tag, eff, cau, kds in ty_cases:
            for k, dm in kds:
                row = ty_row(eff, cau, k, dm)
                fh.write("grty_%s_p%d_d%d\t%s\n"
                         % (tag, k, dm, " ".join("%.17g" % v for v in row)))
        for tag, eff, cau, pls in hac_cases:
            for p, L in pls:
                W, smf, smp = hac_wald(eff, cau, p, L)
                row = [W, float(sps.chi2.sf(W, p)), smf, smp,
                       p, len(eff) - p - (2 * p + 1)]
                fh.write("grhac_%s_p%d_l%d\t%s\n"
                         % (tag, p, L, " ".join("%.17g" % v for v in row)))
    print("wrote", out)


if __name__ == "__main__":
    main()
