# -*- coding: utf-8 -*-
# gen_phase8_ss_inc.py — Phase 8 M2 SARIMAX/statespace 基准 .inc 生成器
#
# 用途: test_arima_statespace.cpp 的 statsmodels 0.14.6 对照基准
#   (黄金锚 llf 固定参数 + MLE fit 结果 + Monahan 变换锚 + 互检 gap 锚)
#
# 语义口径 (probe_sarimax.py 冻结 + 本脚本 [I] 段互检实验 2026-08-23):
#   - 全部 concentrate_scale=False (K4 非集中化, params 含 σ² 末位)
#   - llf 均为 burn=d 口径 (statsmodels 默认; burn=贡献求和起点)
#   - embedded burn=d vs simple/dy 精确似然: 逐项时间平移等价
#     (t=1..T−1 ↔ t=0..T−2, v/F 逐项相同), gap = O(1/κ)=1e-6 量级
#     (κ=1e6 有限近似尾差) — K1 互检告警线 1e-5 的依据
#   - embedded burn=0 含 diffuse 伪贡献 (t=0 的 −½log(κ) ≈ −7.8)
#     — SS5: 与差分样本口径不可混比
#   - trend='c' 常数注入 state_intercept[ARMA 块首行] (d=0: c[0] 均值
#     截距; d≥1: c[d] drift) — 本轮探针实证, SPEC "d=0 观测截距" 描述
#     修正为 state_intercept (statsmodels 对齐)
#
# 重新生成: & C:\Users\Peng\.conda\envs\open-webui\python.exe gen_phase8_ss_inc.py
import pathlib

import numpy as np
from statsmodels.tsa.statespace.sarimax import SARIMAX
from statsmodels.tsa.statespace.tools import (
    constrain_stationary_univariate as constrain,
    unconstrain_stationary_univariate as unconstrain,
)

OUT = pathlib.Path(__file__).resolve().parents[3] / "tests" / "unit" / \
    "timeseries" / "ss_baseline.inc"


def fmt(v) -> str:
    return repr(float(v))


def manual_innovations_llf(z, phi, theta, sigma2):
    """手工 Durbin-Levinson innovations 精确似然 (B&D 5.1.3; 独立双盲实现,
    与 C++ innovations_mle 同构) — 零均值 ARMA(p,q) 全样本 loglik"""
    n = len(z)
    p, q = len(phi), len(theta)
    # ψ 权重 (截断 n 项)
    psi = np.zeros(n)
    psi[0] = 1.0
    for j in range(1, n):
        s = 0.0
        for i in range(1, min(j, p) + 1):
            s += phi[i - 1] * psi[j - i]
        if j <= q:
            s += theta[j - 1]
        psi[j] = s
    # γ_l = σ² Σ_{k≤n−1−l} ψ_k ψ_{k+l}
    gamma = np.array([sigma2 * np.sum(psi[:n - l] * psi[l:])
                      for l in range(n)])
    # Durbin-Levinson: φ_{t,·} 预测系数, v_t 新息方差
    v = np.zeros(n)
    v[0] = gamma[0]
    phi_dl = np.zeros((n, n))  # phi_dl[t, j-1] = φ_{t,j}
    ll = -0.5 * (np.log(2 * np.pi * v[0]) + z[0] ** 2 / v[0])
    for t in range(1, n):
        if t == 1:
            kappa = gamma[1] / gamma[0]
            phi_dl[1, 0] = kappa
        else:
            s = gamma[t]
            for j in range(1, t):
                s -= phi_dl[t - 1, j - 1] * gamma[t - j]
            kappa = s / v[t - 1]
            for j in range(1, t):
                phi_dl[t, j - 1] = (phi_dl[t - 1, j - 1]
                                    - kappa * phi_dl[t - 1, t - 1 - j])
            phi_dl[t, t - 1] = kappa
        v[t] = v[t - 1] * (1 - kappa * kappa)
        e = z[t]
        for j in range(1, t + 1):
            e -= phi_dl[t, j - 1] * z[t - j]
        ll += -0.5 * (np.log(2 * np.pi * v[t]) + e * e / v[t])
    return ll


def gen_arma11(T=120, seed=88):
    """与 probe_sarimax.py 同款生成器 (黄金锚数据可复算)"""
    rng = np.random.RandomState(seed)
    burn = 200
    n = T + burn
    e = rng.standard_normal(n)
    y = np.zeros(n)
    for t in range(1, n):
        y[t] = 0.6 * y[t - 1] + e[t] + 0.4 * e[t - 1]
    return y[burn:]


def main() -> None:
    y = gen_arma11(T=120)
    y300 = gen_arma11(T=300, seed=88)

    gold = {}
    # ---- [G] 黄金锚 llf (固定参数, burn=d 默认) ----
    mod = SARIMAX(y, order=(1, 0, 1), concentrate_scale=False,
                  enforce_stationarity=True, enforce_invertibility=True)
    gold["GOLD_LL_ARMA11"] = mod.loglike(np.array([0.55, 0.35, 0.9]))

    mod11 = SARIMAX(y, order=(1, 1, 1), concentrate_scale=False,
                    enforce_stationarity=True, enforce_invertibility=True)
    gold["GOLD_LL_ARIMA111"] = mod11.loglike(np.array([0.55, 0.35, 0.9]))

    modd = SARIMAX(y, order=(0, 1, 1), trend="c", concentrate_scale=False,
                   enforce_stationarity=True, enforce_invertibility=True)
    gold["GOLD_LL_ARIMA011_DRIFT"] = modd.loglike(np.array([0.2, 0.5, 0.9]))

    # SS5 双口径: simple differencing (Δy 样本 ARMA(1,1) burn=0)
    dy = np.diff(y)
    modsd = SARIMAX(dy, order=(1, 0, 1), concentrate_scale=False,
                    enforce_stationarity=True, enforce_invertibility=True)
    gold["GOLD_LL_ARIMA111_SIMPLE"] = modsd.loglike(np.array([0.55, 0.35, 0.9]))
    # embedded burn=0 (含 diffuse 伪贡献, 差异显式)
    modb0 = SARIMAX(y, order=(1, 1, 1), concentrate_scale=False,
                    enforce_stationarity=True, enforce_invertibility=True,
                    loglikelihood_burn=0)
    gold["GOLD_LL_ARIMA111_BURN0"] = modb0.loglike(np.array([0.55, 0.35, 0.9]))

    # ---- [M] MLE 基准 (fit, 非集中化) ----
    mle = {}
    fits = [
        ("ARMA11", (1, 0, 1), False),
        ("ARIMA111", (1, 1, 1), False),
        ("ARIMA011D", (0, 1, 1), True),      # drift (trend='c')
        ("ARMA22", (2, 0, 2), False),        # 高阶 Monahan
    ]
    for name, order, trend_c in fits:
        m = SARIMAX(y, order=order,
                    trend="c" if trend_c else None,
                    concentrate_scale=False,
                    enforce_stationarity=True, enforce_invertibility=True)
        r = m.fit(disp=False)
        mle[f"MLE_{name}_PARAMS"] = list(r.params)
        mle[f"MLE_{name}_LLF"] = float(r.llf)
        mle[f"MLE_{name}_AIC"] = float(r.aic)
        mle[f"MLE_{name}_K"] = len(r.params)
        print(f"{name}: params={np.round(r.params, 6)} llf={r.llf:.8f} "
              f"aic={r.aic:.6f} k={len(r.params)}")

    # ---- [T] Monahan 变换锚 (constrain 输出 = statsmodels AR 系数口径) ----
    mon_in = [np.array([0.3, -0.2, 0.5, -0.1]),
              np.array([1.5, 0.8, -2.0, 0.4, 1.2]),
              np.array([0.0, 0.0]),
              np.array([-3.0])]
    mon_out = [constrain(u) for u in mon_in]
    # 往返校验
    for u, c in zip(mon_in, mon_out):
        back = unconstrain(c)
        assert np.max(np.abs(back - u)) < 1e-12, f"Monahan 往返失败 {u}"

    # ---- [I] 互检 gap 锚 (T=300, K1 协议) ----
    # d=0: KF(steady) llf vs 手工 Durbin-Levinson 精确似然 (同一数据
    #      y300) — 数学恒等, 预期 1e-8 级 (浮点尾差)
    # d=1: embedded(burn=1) vs simple(dy 样本) — 时间平移等价,
    #      gap = O(1/κ) ≈ 1e-6 (κ=1e6 有限 diffuse 尾差)
    pe = np.array([0.6, 0.4, 1.0])
    me0 = SARIMAX(y300, order=(1, 0, 1), concentrate_scale=False,
                  enforce_stationarity=True, enforce_invertibility=True)
    ll_d0_kf = me0.loglike(pe)
    ll_d0_dl = manual_innovations_llf(y300, [0.6], [0.4], 1.0)
    gap_d0 = abs(ll_d0_kf - ll_d0_dl)
    me = SARIMAX(y300, order=(1, 1, 1), concentrate_scale=False,
                 enforce_stationarity=True, enforce_invertibility=True)
    ll_emb = me.loglike(pe)
    ms = SARIMAX(np.diff(y300), order=(1, 0, 1), concentrate_scale=False,
                 enforce_stationarity=True, enforce_invertibility=True)
    ll_sd = ms.loglike(pe)
    gap_d1 = abs(ll_emb - ll_sd)
    print(f"互检 gap: d=0 {gap_d0:.3e}  d=1 {gap_d1:.3e}")
    gold["GAP_T300_D0"] = gap_d0
    gold["GAP_T300_D1"] = gap_d1

    # ---- 写 .inc ----
    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("// ss_baseline.inc - M2 statespace 基准 (自动生成\n")
        f.write("//   gen_phase8_ss_inc.py, statsmodels 0.14.6, 2026-08-23)\n")
        f.write("// 数据: gen_arma11(T=120/300, seed=88) 与 probe_sarimax.py 同款\n")
        f.write("// 口径: concentrate_scale=False (K4 非集中化); llf burn=d;\n")
        f.write("//   embedded(burn=d) vs simple(dy) 时间平移等价 gap=O(1/κ);\n")
        f.write("//   burn=0 含 diffuse 伪贡献 (SS5 不可混比)\n")
        f.write("#pragma once\n")
        f.write("#include \"cpphub/core/types.hpp\"\n")
        f.write("namespace cpphub { inline namespace v1 { namespace timeseries {\n")
        f.write("namespace ss_baseline {\ninline namespace v1 {\n\n")

        f.write("inline constexpr Size T = 120;\n")
        f.write("inline constexpr Size T300 = 300;\n")
        f.write(f"inline constexpr double DATA[{len(y)}] = {{\n")
        for i in range(0, len(y), 4):
            f.write("    " + ", ".join(fmt(v) for v in y[i:i + 4]) + ",\n")
        f.write("};\n\n")
        f.write(f"inline constexpr double DATA300[{len(y300)}] = {{\n")
        for i in range(0, len(y300), 4):
            f.write("    " + ", ".join(fmt(v) for v in y300[i:i + 4]) + ",\n")
        f.write("};\n\n")

        for k, v in gold.items():
            f.write(f"inline constexpr double {k} = {fmt(v)};\n")
        f.write("\n")

        for k, v in mle.items():
            if isinstance(v, list):
                f.write(f"inline constexpr double {k}[] = {{\n")
                for i in range(0, len(v), 4):
                    f.write("    " + ", ".join(fmt(x) for x in v[i:i + 4])
                            + ",\n")
                f.write("};\n")
            elif isinstance(v, float):
                f.write(f"inline constexpr double {k} = {fmt(v)};\n")
            else:
                f.write(f"inline constexpr Size {k} = {v};\n")
        f.write("\n")

        f.write("// Monahan 变换锚: constrain(u) → statsmodels AR 系数口径\n")
        f.write("// (MA 系数 = −constrain(u), probe [C]); 往返 ≤1e-12\n")
        f.write("inline constexpr double MONAHAN_IN[][5] = {\n")
        for u in mon_in:
            row = list(u) + [0.0] * (5 - len(u))
            f.write("    {" + ", ".join(fmt(v) for v in row) + "},\n")
        f.write("};\n")
        f.write("inline constexpr Size MONAHAN_LEN[] = {"
                + ", ".join(str(len(u)) for u in mon_in) + "};\n")
        f.write("inline constexpr double MONAHAN_OUT[][5] = {\n")
        for c in mon_out:
            row = list(c) + [0.0] * (5 - len(c))
            f.write("    {" + ", ".join(fmt(v) for v in row) + "},\n")
        f.write("};\n\n")

        f.write("}  // inline namespace v1\n}  // namespace ss_baseline\n")
        f.write("}  // namespace timeseries\n}  // namespace v1\n}  // namespace cpphub\n")
    print(f"OK 写入 {OUT}")


if __name__ == "__main__":
    main()
