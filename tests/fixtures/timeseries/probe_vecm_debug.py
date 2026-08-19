# -*- coding: utf-8 -*-
# probe_vecm_debug.py - VECM C++ 管线逐步复现 vs statsmodels 内部量 (定位分歧)

import numpy as np
from numpy.linalg import inv
from statsmodels.tsa.vector_ar.vecm import VECM, _endog_matrices, _sij

y = np.loadtxt(
    r"F:\Cpp_Hub\tests\fixtures\timeseries\coint_smoke_data.csv",
    delimiter=",", skiprows=1)
K, T = 3, y.shape[0]
k_ar_diff, rank, det = 1, 1, "n"

y_1_T, delta_y_1_T, y_lag1, delta_x = _endog_matrices(
    y.T, None, None, k_ar_diff, det)
print("shapes:", y_1_T.shape, delta_y_1_T.shape, y_lag1.shape, delta_x.shape)

s00, s01, s10, s11, s11_, lambd, v = _sij(delta_x, delta_y_1_T, y_lag1)
print("SM λ:", np.round(lambd, 8))

# ---- C++ 管线复现 (det=n, p=2, t_eff=248) ----
p = k_ar_diff + 1
t_eff = T - p
dy = np.diff(y, axis=0)                       # (T−1)×K
dy1t = dy[p - 1:, :]                            # (248)×K
y_lag1_cpp = y[p - 1:T - 1, :].T               # K×248
dx_cpp = np.zeros((K * k_ar_diff, t_eff))
for i in range(t_eff):
    t_reg = p + i
    for j in range(k_ar_diff):
        for vv in range(K):
            dx_cpp[j * K + vv, i] = dy[t_reg - 2 - j, vv]

# R0/R1 = OLS 残差
B0 = dy1t.T @ dx_cpp.T @ inv(dx_cpp @ dx_cpp.T)
r0 = dy1t.T - B0 @ dx_cpp
B1 = y_lag1_cpp @ dx_cpp.T @ inv(dx_cpp @ dx_cpp.T)
r1 = y_lag1_cpp - B1 @ dx_cpp
s00_c = r0 @ r0.T / t_eff
s01_c = r0 @ r1.T / t_eff
s11_c = r1 @ r1.T / t_eff

print("s00 diff:", np.max(np.abs(s00 - s00_c)))
print("s01 diff:", np.max(np.abs(s01 - s01_c)))
print("s11 diff:", np.max(np.abs(s11 - s11_c)))

# s11_ 与 λ
ev = np.linalg.eigh(s11_c)
s11_c_ = ev[1] @ np.diag(1 / np.sqrt(ev[0])) @ ev[1].T
print("s11_ diff:", np.max(np.abs(s11_ - s11_c_)))

M = s01_c @ s11_c_
meig = M.T @ inv(s00_c) @ M
lam_c, v_c = np.linalg.eigh(meig)
lam_c = lam_c[::-1]
v_c = v_c[:, ::-1]
print("C++ λ:", np.round(lam_c, 8))
print("λ diff:", np.max(np.abs(lambd - lam_c)))

beta_c = s11_c_ @ v_c[:, :rank]
beta_c = beta_c @ inv(beta_c[:rank])
alpha_c = s01_c @ beta_c @ inv(beta_c.T @ s11_c @ beta_c)
print("C++ beta:", beta_c.ravel())
print("C++ alpha:", alpha_c.ravel())

# statsmodels 结果
res = VECM(y, k_ar_diff=k_ar_diff, coint_rank=rank, deterministic=det).fit()
print("SM beta:  ", res.beta.ravel())
print("SM alpha: ", res.alpha.ravel())
