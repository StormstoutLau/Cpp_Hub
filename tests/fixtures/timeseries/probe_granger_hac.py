# -*- coding: utf-8 -*-
# probe_granger_hac.py - statsmodels HAC f_test 语义探测 (一次性)
# 问题: get_robustcov_results('HAC').f_test 的 fvalue/pvalue 用什么分布约定?
import os

import numpy as np
from scipy import stats
from statsmodels.regression.linear_model import OLS

HERE = os.path.dirname(os.path.abspath(__file__))

d = np.loadtxt(os.path.join(HERE, "granger_smoke_data.csv"), delimiter=",", skiprows=1)
y, x = d[:, 1], d[:, 0]  # effect=y, cause=x

p = 2
T = len(y)
t = np.arange(p, T)
lhs = y[t]
own = np.column_stack([y[t - j] for j in range(1, p + 1)])
ext = np.column_stack([x[t - j] for j in range(1, p + 1)])
Xu = np.hstack([own, ext, np.ones((len(t), 1))])  # [own, ext, const] (SM prepend=False)
Xr = np.hstack([own, np.ones((len(t), 1))])

res_u = OLS(lhs, Xu).fit()
L = 3
res_hac = res_u.get_robustcov_results(cov_type="HAC", maxlags=L)
R = np.column_stack((np.zeros((p, p)), np.eye(p), np.zeros((p, 1))))
ft = res_hac.f_test(R)
print("f_test: fvalue=%r pvalue=%r df_num=%r df_denom=%r"
      % (float(ft.fvalue), float(ft.pvalue), ft.df_num, ft.df_denom))
print("use_t:", res_hac.use_t)

# 显式 NW 三明治 (无小样本修正, cov_hac_simple use_correction=False)
b = res_u.params
n, k = Xu.shape
XtXi = np.linalg.inv(Xu.T @ Xu)
u = res_u.resid
xu = Xu * u[:, None]
S = xu.T @ xu
for lag in range(1, L + 1):
    s = xu[lag:].T @ xu[:-lag]
    S += (1.0 - lag / (L + 1.0)) * (s + s.T)
V = XtXi @ S @ XtXi
bR = b[p:2 * p]
VR = V[p:2 * p, p:2 * p]
W = float(bR @ np.linalg.inv(VR) @ bR)
print("explicit W(chi2) =", W)
print("chi2.sf(W, p)   =", stats.chi2.sf(W, p))
print("fvalue == W/p?  ", abs(float(ft.fvalue) - W / p))
print("fvalue == W?    ", abs(float(ft.fvalue) - W))
print("pvalue==chi2?   ", abs(float(ft.pvalue) - stats.chi2.sf(W, p)))
print("pvalue==F?      ", abs(float(ft.pvalue) - stats.f.sf(W / p, p, n - k)))

# TY: 增广回归 (k+d lags), 约束前 k 个 cause 滞后
k_lag, dmax = 2, 1
pp = k_lag + dmax
t2 = np.arange(pp, T)
lhs2 = y[t2]
own2 = np.column_stack([y[t2 - j] for j in range(1, pp + 1)])
ext2 = np.column_stack([x[t2 - j] for j in range(1, pp + 1)])
Xt = np.hstack([own2, ext2, np.ones((len(t2), 1))])
rt = OLS(lhs2, Xt).fit()
Rty = np.zeros((k_lag, Xt.shape[1]))
for i in range(k_lag):
    Rty[i, pp + i] = 1.0  # cause lags 1..k (GR4: 增广阶不进约束)
fty = rt.f_test(Rty)
print("\nTY f_test: F=%r p=%r df=(%r,%r)" % (float(fty.fvalue), float(fty.pvalue),
                                             fty.df_num, fty.df_denom))
Wty = float(fty.fvalue) * k_lag
print("TY Wald chi2 = F*k =", Wty, " chi2.sf:", stats.chi2.sf(Wty, k_lag))
# 显式 Wald (df 修正 sigma2)
bt = rt.params
n2, k2 = Xt.shape
s2 = rt.ssr / (n2 - k2)
XtXi2 = np.linalg.inv(Xt.T @ Xt)
bR2 = bt[pp:pp + k_lag]
VR2 = s2 * XtXi2[pp:pp + k_lag, pp:pp + k_lag]
W2 = float(bR2 @ np.linalg.inv(VR2) @ bR2)
print("explicit TY W =", W2, " diff:", abs(W2 - Wty))
