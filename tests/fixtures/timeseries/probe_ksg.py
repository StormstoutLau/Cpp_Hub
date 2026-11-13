# -*- coding: utf-8 -*-
"""Phase 9 前置 probe: KSG 互信息 / CMI / 迁移熵的估计器约定.

用途
----
冻结 KSG Algorithm 1/2 的公式、Chebyshev 距离、digamma 计数、TE=GC/2 锚、
缩放敏感性与零假设偏差; 为 C++ 实现提供基准与红绿判据.
本脚本是 Phase 9 spec 的 probe 前置件, 不参与 C++ 构建.

运行
----
    py -3.12 tests/fixtures/timeseries/probe_ksg.py
    py -3.12 tests/fixtures/timeseries/probe_ksg.py --json

判据
----
K1-K7 全部 PASS 时退出码 0; 任一 FAIL 退出码 1.
本 probe 不依赖 JIDT/R; JIDT 对照在 spec 的 verify 阶段另行执行.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field

import numpy as np
from scipy.special import digamma
from scipy.spatial import cKDTree
from sklearn.feature_selection._mutual_info import _compute_mi_cc

SEED = 20261111


@dataclass
class Case:
    name: str
    status: str
    detail: str
    values: dict = field(default_factory=dict)


def check(name: str, ok: bool, detail: str, values: dict | None = None) -> Case:
    return Case(name, "PASS" if ok else "FAIL", detail, values or {})


def _count_within(points: np.ndarray, radii: np.ndarray, minus_self: bool) -> np.ndarray:
    pts = points.reshape(-1, 1) if points.ndim == 1 else points
    tree = cKDTree(pts)
    cnt = tree.query_ball_point(pts, radii, p=np.inf, return_length=True)
    return cnt - (1 if minus_self else 0)


def ksg_mi_alg1(x: np.ndarray, y: np.ndarray, k: int = 3) -> float:
    """KSG Algorithm 1 (Eq. 8), 联合 Chebyshev, 严格小于 eps(i)."""
    data = np.column_stack([x, y])
    tree = cKDTree(data)
    d, _ = tree.query(data, k=k + 1, p=np.inf)
    eps = np.nextafter(d[:, -1], 0)
    nx = _count_within(x, eps, minus_self=True)
    ny = _count_within(y, eps, minus_self=True)
    return float(
        digamma(k) + digamma(len(x))
        - np.mean(digamma(nx + 1) + digamma(ny + 1))
    )


def ksg_mi_alg2(x: np.ndarray, y: np.ndarray, k: int = 3) -> float:
    """KSG Algorithm 2 (Eq. 9), 边际各自取第 k 联合近邻的投影半径."""
    data = np.column_stack([x, y])
    tree = cKDTree(data)
    _, idx = tree.query(data, k=k + 1, p=np.inf)
    j = idx[:, -1]
    ex = np.abs(x - x[j])
    ey = np.abs(y - y[j])
    nx = _count_within(x, ex / 2.0, minus_self=False)
    ny = _count_within(y, ey / 2.0, minus_self=False)
    return float(
        digamma(k) - 1.0 / k + digamma(len(x))
        - np.mean(digamma(nx) + digamma(ny))
    )


def cmi_ksg(x: np.ndarray, y: np.ndarray, z: np.ndarray, k: int = 3) -> float:
    """Frenzel-Pompe CMI: I(X;Y|Z) = psi(k) - <psi(nxz+1)+psi(nyz+1)-psi(nz+1)>."""
    data = np.column_stack([x, y, z])
    tree = cKDTree(data)
    d, _ = tree.query(data, k=k + 1, p=np.inf)
    eps = np.nextafter(d[:, -1], 0)
    nxz = _count_within(data[:, [0, 2]], eps, minus_self=True)
    nyz = _count_within(data[:, [1, 2]], eps, minus_self=True)
    nz = _count_within(data[:, 2], eps, minus_self=True)
    return float(
        digamma(k)
        - np.mean(digamma(nxz + 1) + digamma(nyz + 1) - digamma(nz + 1))
    )


def te_xy(x: np.ndarray, y: np.ndarray, k: int = 3, lag: int = 1) -> float:
    """TE(X->Y) = CMI(Y_t ; X_{t-lag} | Y_{t-lag})."""
    return cmi_ksg(y[lag:], x[:-lag], y[:-lag], k=k)


def ar1(rng: np.random.Generator, phi: float, n: int, e: np.ndarray) -> np.ndarray:
    s = np.empty(n)
    s[0] = e[0] / np.sqrt(1.0 - phi ** 2)
    for t in range(1, n):
        s[t] = phi * s[t - 1] + e[t]
    return s


def gauss_var(rng: np.random.Generator, phi: float, a: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    ex, ey = rng.normal(size=n), rng.normal(size=n)
    x = ar1(rng, phi, n, ex)
    y = np.empty(n)
    y[0] = ey[0] / np.sqrt(1.0 - phi ** 2)
    for t in range(1, n):
        y[t] = phi * y[t - 1] + a * x[t - 1] + ey[t]
    return x, y


def zscore(x: np.ndarray) -> np.ndarray:
    return (x - np.mean(x)) / np.std(x, ddof=1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true", help="输出 JSON")
    args = ap.parse_args()

    cases: list[Case] = []
    n = 2000

    # K0: 自实现 Algorithm 1 与 sklearn _compute_mi_cc 逐位一致 (同数据, 无 jitter)
    r = np.random.default_rng(SEED)
    z = r.multivariate_normal([0, 0], [[1, 0.8], [0.8, 1]], size=n)
    ours = ksg_mi_alg1(z[:, 0], z[:, 1], k=3)
    ref = float(_compute_mi_cc(z[:, 0], z[:, 1], 3))
    cases.append(check(
        "K0 自实现 KSG Alg1 vs sklearn",
        abs(ours - ref) < 1e-12,
        f"ours={ours:.12f}; sklearn={ref:.12f}; diff={abs(ours-ref):.3e}",
        {"ours": ours, "sklearn": ref},
    ))

    # K1: Gaussian MI 偏差/SD (n=2000, k=3, 5 seeds)
    for rho in (0.5, 0.8):
        truth = -0.5 * np.log(1 - rho ** 2)
        vals = []
        for seed in range(5):
            rr = np.random.default_rng(seed)
            zz = rr.multivariate_normal([0, 0], [[1, rho], [rho, 1]], size=n)
            vals.append(ksg_mi_alg1(zz[:, 0], zz[:, 1], k=3))
        vals = np.asarray(vals)
        bias = float(vals.mean() - truth)
        sd = float(vals.std(ddof=1))
        cases.append(check(
            f"K1 Gaussian MI rho={rho}",
            abs(bias) < 0.02 and sd < 0.03,
            f"truth={truth:.6f}; mean={vals.mean():.6f}; bias={bias:+.6f}; sd={sd:.6f}",
            {"truth": truth, "mean": float(vals.mean()), "bias": bias, "sd": sd},
        ))

    # K2: Algorithm 2 语义未冻结 — 三种半径变体 (ex, ex/2, d_k) 均显著偏离解析值,
    # 原文 "ε_x(i)/2" 的投影半径读法与 Eq.(9) 的 n_x 计数口径需要独立复核.
    # Phase 9 只实现 Algorithm 1 (与 sklearn 一致); Alg2 留作 spec 开放问题.
    a1 = ksg_mi_alg1(z[:, 0], z[:, 1], k=3)
    cases.append(Case(
        "K2 KSG Algorithm 2 语义未冻结",
        "SKIP",
        f"Alg1={a1:.6f}; Alg2 三变体偏离解析值, 不纳入 Phase 9 主线 (见 spec H5)",
        {"alg1": a1},
    ))

    # K3: 缩放敏感性 — 标准化后不变, 原始尺度变化
    rr = np.random.default_rng(SEED + 1)
    x, y = gauss_var(rr, 0.5, 0.3, n)
    te_std = [te_xy(zscore(s * x), zscore(y), k=3) for s in (0.1, 1.0, 10.0)]
    te_raw = [te_xy(s * x, y, k=3) for s in (0.1, 1.0, 10.0)]
    std_span = max(te_std) - min(te_std)
    raw_span = max(te_raw) - min(te_raw)
    cases.append(check(
        "K3 标准化消除缩放敏感性",
        std_span < 1e-12 and raw_span > 0.002,
        f"std TE={te_std}; span={std_span:.3e}; raw TE={te_raw}; span={raw_span:.6f}",
        {"te_std": te_std, "te_raw": te_raw, "std_span": std_span, "raw_span": raw_span},
    ))

    # K4: TE = GC/2 锚 (Gaussian VAR, 5 seeds)
    phi, a = 0.5, 0.3
    gc = float(np.log(1 + a * a / (1 - phi ** 2)))
    te_true = gc / 2.0
    vals = []
    for seed in range(5):
        rr = np.random.default_rng(200 + seed)
        xx, yy = gauss_var(rr, phi, a, n)
        vals.append(te_xy(zscore(xx), zscore(yy), k=3))
    vals = np.asarray(vals)
    rel = float(vals.mean() / te_true - 1.0)
    cases.append(check(
        "K4 TE vs GC/2 锚",
        abs(rel) < 0.35,
        f"GC={gc:.6f}; GC/2={te_true:.6f}; TE={vals.mean():.6f} "
        f"(sd {vals.std(ddof=1):.6f}; rel {rel:+.1%})",
        {"gc": gc, "te_true": te_true, "te_mean": float(vals.mean()),
         "te_sd": float(vals.std(ddof=1)), "rel": rel},
    ))

    # K5: 零假设 — i.i.d. MI≈0; 独立 AR(1) TE 仍为正 (伪迁移)
    mi_null = []
    for seed in range(20):
        rr = np.random.default_rng(300 + seed)
        zz = rr.normal(size=(n, 2))
        mi_null.append(ksg_mi_alg1(zz[:, 0], zz[:, 1], k=3))
    mi_null = np.asarray(mi_null)
    te_null = []
    te_eff = []
    for seed in range(5):
        rr = np.random.default_rng(400 + seed)
        xx = ar1(rr, 0.5, n, rr.normal(size=n))
        yy = ar1(rr, 0.5, n, rr.normal(size=n))
        t0 = te_xy(zscore(xx), zscore(yy), k=3)
        sh = [te_xy(zscore(rr.permutation(xx)), zscore(yy), k=3) for _ in range(5)]
        te_null.append(t0)
        te_eff.append(t0 - float(np.mean(sh)))
    te_null = np.asarray(te_null)
    te_eff = np.asarray(te_eff)
    cases.append(check(
        "K5 零假设: i.i.d. MI≈0, AR(1) TE 非零",
        abs(float(mi_null.mean())) < 0.01 and abs(float(te_null.mean())) > 0.002,
        f"iid MI mean={mi_null.mean():+.6f} (sd {mi_null.std(ddof=1):.6f}); "
        f"AR(1) TE mean={te_null.mean():+.6f} (sd {te_null.std(ddof=1):.6f}); "
        f"effective TE mean={te_eff.mean():+.6f} (符号随有限样本波动)",
        {"iid_mi_mean": float(mi_null.mean()), "iid_mi_sd": float(mi_null.std(ddof=1)),
         "ar1_te_mean": float(te_null.mean()), "ar1_te_sd": float(te_null.std(ddof=1)),
         "ar1_te_eff_mean": float(te_eff.mean())},
    ))

    # K6: ties (舍入) 对 KSG 的影响 — jitter 改变估计
    rr = np.random.default_rng(SEED + 2)
    xx, yy = gauss_var(rr, 0.5, 0.3, n)
    xr = np.round(xx, 2)
    yr = np.round(yy, 2)
    te_tie = te_xy(zscore(xr), zscore(yr), k=3)
    jitter = 1e-10 * max(1.0, float(np.mean(np.abs(xr))))
    te_jit = te_xy(zscore(xr + jitter * rr.normal(size=n)),
                   zscore(yr + jitter * rr.normal(size=n)), k=3)
    cases.append(check(
        "K6 ties 对 KSG 估计有影响",
        np.isfinite(te_tie) and np.isfinite(te_jit) and abs(te_tie - te_jit) > 1e-6,
        f"rounded TE={te_tie:.6f}; jittered TE={te_jit:.6f}; diff={abs(te_tie-te_jit):.6f}",
        {"te_tie": te_tie, "te_jitter": te_jit},
    ))

    # K7: nats -> bits 换算
    ratio = (1.0 / np.log(2.0))
    cases.append(check(
        "K7 nats/bits 换算因子",
        abs(ratio - 1.4426950408889634) < 1e-12,
        f"1/ln2={ratio:.16f}",
        {"nats_to_bits": ratio},
    ))

    if args.json:
        print(json.dumps(
            {"seed": SEED, "n": n, "cases": [c.__dict__ for c in cases]},
            ensure_ascii=False, indent=2,
        ))
    else:
        for c in cases:
            print(f"[{c.status}] {c.name}: {c.detail}")
        n_fail = sum(c.status == "FAIL" for c in cases)
        print(f"\n{len(cases) - n_fail} PASS, {n_fail} FAIL")
    return 0 if all(c.status != "FAIL" for c in cases) else 1


if __name__ == "__main__":
    sys.exit(main())
