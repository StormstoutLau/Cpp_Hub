# -*- coding: utf-8 -*-
# verify_drift_init.py — drift 混合初始化语义验证 (arima_statespace.hpp 前置)
#
# 待实证: ARIMA(0,1,1)+drift 的 statsmodels llf 是否等于手工 KF:
#   a1 = (0, drift, 0)  [差分块 0 + ARMA 块平稳均值 (I−T_sub)⁻¹c_sub]
#   P1 = diag(κ=1e6) ⊕ Lyapunov(T_sub, R_sub Q R_sub')
#   c_vec = (0, drift, 0)  [state_intercept, 注入 index d=1]
#
# 依据源码: initialization.py L718-721 (stationary mean = solve(I−T, c))
#   + L681-686 (block 递归用 sub-block 矩阵) — 本脚本数值闭环
#
# 用法: & C:\Users\Peng\.conda\envs\open-webui\python.exe verify_drift_init.py
import numpy as np
from statsmodels.tsa.statespace.sarimax import SARIMAX


def gen_arma11(T=120, seed=88):
    rng = np.random.RandomState(seed)
    burn = 200
    n = T + burn
    e = rng.standard_normal(n)
    y = np.zeros(n)
    for t in range(1, n):
        y[t] = 0.6 * y[t - 1] + e[t] + 0.4 * e[t - 1]
    return y[burn:]


def manual_kf(y, Z, T, R, Q, H, a1, P1, c=None, burn=0):
    T_n = len(y)
    m = T.shape[0]
    if c is None:
        c = np.zeros(m)
    a = a1.astype(float).copy()
    P = P1.astype(float).copy()
    Zf = np.asarray(Z).flatten()
    Rm = np.asarray(R).reshape(m, -1)
    Qm = np.asarray(Q).reshape(Rm.shape[1], Rm.shape[1])
    h = float(np.asarray(H).flatten()[0])
    ll = 0.0
    for t in range(T_n):
        v = y[t] - Zf @ a
        F = float(Zf @ P @ Zf + h)
        if t >= burn:
            ll += -0.5 * (np.log(F) + v * v / F + np.log(2 * np.pi))
        K = T @ P @ Zf / F
        a = T @ a + K * v + c
        P = T @ P @ (T - np.outer(K, Zf)).T + Rm @ Qm @ Rm.T
    return ll


def lyapunov(T, RQR):
    P = RQR.copy()
    for _ in range(500):
        P = T @ P @ T.T + RQR
    return P


def main():
    y = gen_arma11(T=120)
    drift, theta, s2 = 0.2, 0.5, 0.9

    mod = SARIMAX(y, order=(0, 1, 1), trend="c",
                  concentrate_scale=False, enforce_stationarity=True,
                  enforce_invertibility=True)
    params = np.array([drift, theta, s2])
    ll_sm = mod.loglike(params)
    print(f"statsmodels llf (burn={mod.loglikelihood_burn}) = {ll_sm:.12f}")

    # 手工系统 (probe [A3] 布局: m = d + max(p, q+1) = 1 + 2 = 3)
    Z = np.array([1.0, 1.0, 0.0])
    T = np.array([[1.0, 1.0, 0.0],
                  [0.0, 0.0, 1.0],
                  [0.0, 0.0, 0.0]])
    R = np.array([[0.0], [1.0], [theta]])
    Q = np.array([[s2]])
    H = np.zeros((1, 1))
    c = np.array([0.0, drift, 0.0])  # state_intercept[1] = drift (index d)

    # 混合初始化: 差分块 diffuse κ; ARMA 块 stationary
    kappa = 1e6
    Tsub = T[1:, 1:]                      # [[0,1],[0,0]]
    Rsub = R[1:, :]                        # [1, θ]'
    RQR_sub = Rsub @ Q @ Rsub.T
    Psub = lyapunov(Tsub, RQR_sub)
    # 闭式: σ²[[1+θ², θ],[θ, θ²]]
    Psub_closed = s2 * np.array([[1 + theta**2, theta], [theta, theta**2]])
    print(f"Lyapunov 闭式 diff = {np.max(np.abs(Psub - Psub_closed)):.2e}")

    a1_block = np.linalg.solve(np.eye(2) - Tsub, c[1:])
    print(f"ARMA 块平稳均值 (I−T)⁻¹c = {a1_block} (期望 = [drift, 0] = "
          f"[{drift}, 0.0])")

    a1 = np.zeros(3)
    a1[1:] = a1_block
    P1 = np.zeros((3, 3))
    P1[0, 0] = kappa
    P1[1:, 1:] = Psub

    ll_man = manual_kf(y, Z, T, R, Q, H, a1, P1, c=c, burn=1)
    print(f"手工 KF llf (混合初始化, burn=1) = {ll_man:.12f}")
    print(f"diff = {abs(ll_sm - ll_man):.3e}")
    assert abs(ll_sm - ll_man) < 1e-8, "混合初始化语义与 statsmodels 不符!"

    # 对照: a1 = 0 (无漂移均值) 的错误初始化 → diff 应显著
    ll_wrong = manual_kf(y, Z, T, R, Q, H, np.zeros(3), P1, c=c, burn=1)
    print(f"[对照] a1=0 错误初始化 diff = {abs(ll_sm - ll_wrong):.3e} "
          f"(应 > 1e-8, 证明平稳均值必要)")

    # ---- innovations 互检 (K1): Δy − drift 零均值 MA(1) 精确似然 ----
    # 手工 Durbin-Levinson (gen_phase8_ss_inc.py 同款独立实现)
    def manual_innovations_llf(z, phi, theta, sigma2):
        n = len(z)
        p, q = len(phi), len(theta)
        psi = np.zeros(n)
        psi[0] = 1.0
        for j in range(1, n):
            s = 0.0
            for i in range(1, min(j, p) + 1):
                s += phi[i - 1] * psi[j - i]
            if j <= q:
                s += theta[j - 1]
            psi[j] = s
        gamma = np.array([sigma2 * np.sum(psi[:n - l] * psi[l:])
                          for l in range(n)])
        v = np.zeros(n)
        v[0] = gamma[0]
        phi_dl = np.zeros((n, n))
        ll = -0.5 * (np.log(2 * np.pi * v[0]) + z[0] ** 2 / v[0])
        for t in range(1, n):
            if t == 1:
                kappa_ = gamma[1] / gamma[0]
                phi_dl[1, 0] = kappa_
            else:
                s = gamma[t]
                for j in range(1, t):
                    s -= phi_dl[t - 1, j - 1] * gamma[t - j]
                kappa_ = s / v[t - 1]
                for j in range(1, t):
                    phi_dl[t, j - 1] = (phi_dl[t - 1, j - 1]
                                        - kappa_ * phi_dl[t - 1, t - 1 - j])
                phi_dl[t, t - 1] = kappa_
            v[t] = v[t - 1] * (1 - kappa_ * kappa_)
            e = z[t]
            for j in range(1, t + 1):
                e -= phi_dl[t, j - 1] * z[t - j]
            ll += -0.5 * (np.log(2 * np.pi * v[t]) + e * e / v[t])
        return ll

    dy = np.diff(y) - drift
    ll_dl = manual_innovations_llf(dy, [], [theta], s2)
    print(f"\n[K1] innovations (Δy − drift) llf = {ll_dl:.12f}")
    print(f"[K1] gap |KF − innovations| = {abs(ll_sm - ll_dl):.3e}")
    print("     (d=1 期望 O(1/κ) ≈ 1e-6 量级)")

    print("\n✓ drift 混合初始化语义实证冻结")


if __name__ == "__main__":
    main()
