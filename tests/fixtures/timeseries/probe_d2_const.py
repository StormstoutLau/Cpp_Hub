# -*- coding: utf-8 -*-
# probe_d2_const.py - arima_statespace.hpp 实现前最后两个语义钉子
#
# [P1] d=0 + trend='c': 常数 c 直入 state_intercept[0] (AR 方程截距,
#      非均值 μ)? 用手工 KF (probe_sarimax.py manual_kf 同款) 对照 llf
# [P2] d=2 嵌入差分 T 布局: triu(2×2 全 1) + T[:2, 2]=1 链接列 +
#      companion ARMA 块 — 逐元素对照 statsmodels + llf 全量对齐
#
# 运行: & C:\Users\Peng\.conda\envs\open-webui\python.exe probe_d2_const.py
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
    """手工 Kalman 滤波 (DK Ch4 更新形) — ss_model.hpp 镜像"""
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
    m = T.shape[0]
    A = np.eye(m * m)
    for i in range(m):
        for j in range(m):
            for k in range(m):
                for l in range(m):
                    A[i * m + k, j * m + l] -= T[i, j] * T[k, l]
    rhs = np.zeros(m * m)
    for i in range(m):
        for j in range(m):
            rhs[i + j * m] = RQR[i, j]
    vp = np.linalg.solve(A, rhs)
    P = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            P[i, j] = vp[i + j * m]
    return 0.5 * (P + P.T)


def main():
    y = gen_arma11(T=120)

    # ================= [P1] d=0 + trend='c' =================
    print("=" * 74)
    print("[P1] ARMA(1,1) + trend='c': 常数直入 state_intercept[0]?")
    c0, phi, theta, s2 = 0.3, 0.5, 0.4, 1.0
    mod = SARIMAX(y, order=(1, 0, 1), trend="c", concentrate_scale=False,
                  enforce_stationarity=True, enforce_invertibility=True)
    ll_sm = mod.loglike(np.array([c0, phi, theta, s2]))
    # 手工: c_vec = (c0, 0)' 直入 state_intercept (截距语义)
    Z = np.array([1.0, 0.0])
    T = np.array([[phi, 1.0], [0.0, 0.0]])
    R = np.array([[1.0], [theta]])
    Q = np.array([[s2]])
    H = np.zeros((1, 1))
    c = np.array([c0, 0.0])
    P1_m = lyapunov(T, R @ Q @ R.T)
    a1 = np.linalg.solve(np.eye(2) - T, c)
    ll_man = manual_kf(y, Z, T, R, Q, H, a1, P1_m, c=c, burn=0)
    print(f"  statsmodels llf = {ll_sm:.12f}")
    print(f"  手工 KF (c 直入 state_intercept[0]) = {ll_man:.12f}")
    print(f"  diff = {abs(ll_sm - ll_man):.3e}  (<1e-9 → 截距语义实证)")
    assert abs(ll_sm - ll_man) < 1e-9
    # 对照: 若 c 是均值 μ → 截距应填 μ(1−φ); 用错误语义 diff 应显著
    ll_wrong = manual_kf(y, Z, T, R, Q, H,
                         np.linalg.solve(np.eye(2) - T,
                                         np.array([c0 * (1 - phi), 0.0])),
                         P1_m, c=np.array([c0 * (1 - phi), 0.0]), burn=0)
    print(f"  [对照] 均值语义 (c·(1−φ)) diff = {abs(ll_sm - ll_wrong):.3e} "
          f"(应 >> 1e-9)")

    # ================= [P2] d=2 嵌入差分布局 =================
    print()
    print("=" * 74)
    print("[P2] ARIMA(p,2,q) 嵌入差分: T 逐元素 + llf (burn=2)")
    for (p, q, params) in [
        (0, 1, np.array([0.4, 1.0])),          # ARIMA(0,2,1)
        (1, 1, np.array([0.5, 0.3, 1.0])),     # ARIMA(1,2,1)
        (2, 1, np.array([0.4, -0.2, 0.3, 1.0])),  # ARIMA(2,2,1)
    ]:
        mod = SARIMAX(y, order=(p, 2, q), concentrate_scale=False,
                      enforce_stationarity=True, enforce_invertibility=True)
        m1 = max(p, q + 1)
        m = 2 + m1
        ll_sm = mod.loglike(params)
        # statsmodels 实际矩阵
        mod.update(params)
        T_sm = np.asarray(mod.ssm.transition).reshape(m, m)
        R_sm = np.asarray(mod.ssm.selection).reshape(m, 1)
        Z_sm = np.asarray(mod.ssm.design).reshape(1, m)
        # 手工构造: triu(2×2) + 链接列 T[:2,2]=1 + companion(m1) + AR 首列
        T_man = np.zeros((m, m))
        T_man[0, 0] = T_man[0, 1] = T_man[1, 1] = 1.0   # triu 全 1
        T_man[0, 2] = T_man[1, 2] = 1.0                  # 链接列
        for i in range(m1 - 1):
            T_man[2 + i, 3 + i] = 1.0                    # companion 超对角
        for i in range(p):
            T_man[2 + i, 2] = params[i]                  # AR 首列
        R_man = np.zeros((m, 1))
        R_man[2, 0] = 1.0
        for j in range(q):
            R_man[3 + j, 0] = params[p + j]
        Z_man = np.zeros(m)
        Z_man[:3] = 1.0                                  # [1]*(d+1) + [0]*
        dT = np.max(np.abs(T_sm - T_man))
        dR = np.max(np.abs(R_sm - R_man))
        dZ = np.max(np.abs(Z_sm.flatten() - Z_man))
        # 手工 llf: 混合初始化 (diff 块 κ / ARMA 块 stationary) + burn=2
        kappa = 1e6
        Tsub = T_man[2:, 2:]
        Rsub = R_man[2:]
        Q = np.array([[params[-1]]])
        Psub = lyapunov(Tsub, Rsub @ Q @ Rsub.T)
        a1 = np.zeros(m)
        a1[2:] = np.linalg.solve(np.eye(m1) - Tsub, np.zeros(m1))
        P1_m = np.zeros((m, m))
        P1_m[0, 0] = P1_m[1, 1] = kappa
        P1_m[2:, 2:] = Psub
        ll_man = manual_kf(y, Z_man, T_man, R_man, Q, np.zeros((1, 1)),
                           a1, P1_m, c=None, burn=2)
        print(f"  ({p},2,{q}): dT={dT:.2e} dR={dR:.2e} dZ={dZ:.2e} | "
              f"llf diff = {abs(ll_sm - ll_man):.3e}")
        assert dT < 1e-14 and dR < 1e-14 and dZ < 1e-14
        # d=2 双 diffuse 态: 有限 κ=1e6 尾差 O(1/κ) (~1e-7 量级, 非布局错误)
        assert abs(ll_sm - ll_man) < 1e-6
    print("  ✓ d=2 布局 (triu + 链接列 + companion) 逐元素实证 + llf 对齐")


if __name__ == "__main__":
    main()
