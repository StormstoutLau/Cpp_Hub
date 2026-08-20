# -*- coding: utf-8 -*-
# probe_granger_prec.py - 数值路径探测: 正规方程 vs QR (一次性)
# 目的: 决定 granger_test.hpp 的 OLS 求解路径 (spec 容差 1e-10)
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def build(eff, cau, p):
    T = len(eff)
    t = np.arange(p, T)
    lhs = eff[t]
    own = np.column_stack([eff[t - j] for j in range(1, p + 1)])
    ext = np.column_stack([cau[t - j] for j in range(1, p + 1)])
    one = np.ones((len(t), 1))
    return lhs, np.hstack([own, one]), np.hstack([own, ext, one])


def ols_ne(y, X):
    XtX = X.T @ X
    b = np.linalg.solve(XtX, X.T @ y)
    r = y - X @ b
    return b, r @ r, np.linalg.inv(XtX)


def ols_qr(y, X):
    Q, R = np.linalg.qr(X)
    b = np.linalg.solve(R, Q.T @ y)
    r = y - X @ b
    Rinv = np.linalg.inv(R)
    return b, r @ r, Rinv @ Rinv.T


d = np.loadtxt(os.path.join(HERE, "granger_smoke_data.csv"), delimiter=",", skiprows=1)
ci = np.loadtxt(os.path.join(HERE, "coint_smoke_data.csv"), delimiter=",", skiprows=1)
cases = [
    ("fwd_p1", d[:, 1], d[:, 0], 1),
    ("fwd_p2", d[:, 1], d[:, 0], 2),
    ("i1_p1", ci[:, 0], ci[:, 1], 1),
    ("i1_p2", ci[:, 0], ci[:, 1], 2),
]

ref = {}
for tag, eff, cau, p in cases:
    y, Xr, Xu = build(eff, cau, p)
    _, ssr_r, _ = ols_qr(y, Xr)
    _, ssr_u, _ = ols_qr(y, Xu)
    n = len(y)
    df = n - (2 * p + 1)
    f = (ssr_r - ssr_u) / ssr_u / p * df
    ref[tag] = (ssr_r, ssr_u, f)

for tag, eff, cau, p in cases:
    y, Xr, Xu = build(eff, cau, p)
    for name, fn in [("NE", ols_ne), ("QR", ols_qr)]:
        _, ssr_r, _ = fn(y, Xr)
        _, ssr_u, _ = fn(y, Xu)
        n = len(y)
        df = n - (2 * p + 1)
        f = (ssr_r - ssr_u) / ssr_u / p * df
        s0, u0, f0 = ref[tag]
        print("%-7s %-3s dSSR_r=%9.2e dSSR_u=%9.2e dF=%9.2e  cond(Xu)=%8.2e"
              % (tag, name, abs(ssr_r - s0), abs(ssr_u - u0), abs(f - f0),
                 np.linalg.cond(Xu)))
