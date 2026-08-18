# -*- coding: utf-8 -*-
# verify_arima.py - ARIMA 基准生成: statsmodels innovations_mle (Phase 7C M1)
#
# 对照点 (spec §1.3):
#   ARIMA innovations vs statsmodels innovations_mle (无缺失非季节) 1e-10
#
# 语义对齐 (estimators/innovations.py 一手, 2026-08-18):
#   - statsmodels: 差分 d → demean=True → HR 起始 → scipy minimize
#     全参数 (φ,θ,σ²) on −arma_loglike
#   - C++: 差分 d → demean → {HR,0,扰动} SLSQP on 集中化 nll (σ² 解析)
#   - 同一似然面; 落点差 = 优化器容差层 (实测落档, 7B GM 先例口径)
#   ⚠️ d=0 夹具数据非零均值 (测 demean 语义); C++ innovations_mle(demean=True)
#
# 用法: python verify_arima.py   (生成 arima_smoke_data.csv + 基准 JSON)
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def gen_arma11(T=300, seed=42):
    """ARMA(1,1): φ=0.5, θ=−0.4, σ=1, burn-in 100 (RandomState 跨语言固定)"""
    rng = np.random.RandomState(seed)
    phi, theta, e = 0.5, -0.4, rng.standard_normal(T + 100)
    z = np.zeros(T + 100)
    for t in range(1, T + 100):
        z[t] = phi * z[t - 1] + e[t] + theta * e[t - 1]
    return z[100:]


def gen_arma22(T=300, seed=43):
    """ARMA(2,2): φ=(0.3,−0.2), θ=(−0.5,0.3) — 多峰逃逸夹具 (AR7)"""
    rng = np.random.RandomState(seed)
    phi, theta = [0.3, -0.2], [-0.5, 0.3]
    e = rng.standard_normal(T + 100)
    z = np.zeros(T + 100)
    for t in range(2, T + 100):
        z[t] = (phi[0] * z[t - 1] + phi[1] * z[t - 2] + e[t]
                + theta[0] * e[t - 1] + theta[1] * e[t - 2])
    return z[100:]


def gen_arma21(T=300, seed=45):
    """ARMA(2,1): p≠q — AR2 裁决夹具 (R CSS n.cond = d+max(p,q) vs d+p)"""
    rng = np.random.RandomState(seed)
    phi, theta = [0.4, -0.25], [-0.6]
    e = rng.standard_normal(T + 100)
    z = np.zeros(T + 100)
    for t in range(2, T + 100):
        z[t] = (phi[0] * z[t - 1] + phi[1] * z[t - 2] + e[t]
                + theta[0] * e[t - 1])
    return z[100:]


def gen_arma12(T=300, seed=46):
    """ARMA(1,2): p<q — AR2 定案夹具 (max(p,q)=2 vs p=1 → n.cond 1 or 2)"""
    rng = np.random.RandomState(seed)
    phi, theta = [0.35], [-0.5, 0.25]
    e = rng.standard_normal(T + 100)
    z = np.zeros(T + 100)
    for t in range(2, T + 100):
        z[t] = (phi[0] * z[t - 1] + e[t]
                + theta[0] * e[t - 1] + theta[1] * e[t - 2])
    return z[100:]


def gen_arima111_drift(T=300, seed=44):
    """ARIMA(1,1,1)+drift: Δy = 0.3 + ARMA(1,1)"""
    rng = np.random.RandomState(seed)
    d = gen_arma11(T, seed) + 0.3
    return np.cumsum(d)


def main():
    from statsmodels.tsa.arima.estimators.innovations import innovations_mle

    fixtures = {
        "arma11": {"data": gen_arma11(), "order": (1, 0, 1), "demean": True},
        "arma22": {"data": gen_arma22(), "order": (2, 0, 2), "demean": True},
        "arma21": {"data": gen_arma21(), "order": (2, 0, 1), "demean": True},
        "arma12": {"data": gen_arma12(), "order": (1, 0, 2), "demean": True},
        # innovations+drift: statsmodels 无 drift 参数 → demean=True 吸收
        # (drift ≈ 差分序列均值; C++ 侧 innovations_mle(demean=True) 对齐)
        "arima111d": {"data": np.diff(gen_arima111_drift()),
                      "order": (1, 0, 1), "demean": True},
    }

    out = {}
    for name, fx in fixtures.items():
        p, other = innovations_mle(fx["data"], order=fx["order"],
                                   demean=fx["demean"])
        llf = None
        try:
            from statsmodels.tsa.innovations.arma_innovations import (
                arma_loglike,
            )
            z = fx["data"] - fx["data"].mean()
            llf = float(arma_loglike(
                z, ar_params=-np.asarray(p.reduced_ar_poly.coef[1:]),
                ma_params=np.asarray(p.reduced_ma_poly.coef[1:]),
                sigma2=float(p.sigma2)))
        except Exception as exc:  # noqa: BLE001
            llf = None
            print(f"{name}: llf eval failed: {exc}")
        out[name] = {
            "phi": [float(v) for v in p.ar_params],
            "theta": [float(v) for v in p.ma_params],
            "sigma2": float(p.sigma2),
            "loglik": llf,
            "demean": fx["demean"],
            "order": list(fx["order"]),
            "nobs": int(len(fx["data"])),
            "converged": bool(other["minimize_results"].success),
        }
        print(name, json.dumps(out[name], indent=1))

    # 数据 CSV (%.17g 全精度, C++/R 同源消费; 各列同长 T=300;
    # arima111d 存 level, 差分由消费方自行 diff — 列长一致原则)
    with open(os.path.join(HERE, "arima_smoke_data.csv"), "w") as f:
        cols = {"arma11": fixtures["arma11"]["data"],
                "arma22": fixtures["arma22"]["data"],
                "arma21": fixtures["arma21"]["data"],
                "arma12": fixtures["arma12"]["data"],
                "arima111d_level": gen_arima111_drift()}
        keys = list(cols)
        lens = {len(v) for v in cols.values()}
        assert len(lens) == 1, f"column length mismatch: {lens}"
        n = lens.pop()
        for i in range(n):
            f.write(",".join("%.17g" % cols[k][i] for k in keys) + "\n")
        f.write("#columns=" + ",".join(keys) + "\n")

    with open(os.path.join(HERE, "arima_statsmodels_baselines.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("written: arima_statsmodels_baselines.json")


if __name__ == "__main__":
    main()
