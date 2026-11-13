# -*- coding: utf-8 -*-
"""Phase 9 前置 probe: KDE 带宽语义与跨库数值对照.

用途
----
冻结 sklearn / scipy / statsmodels / R 的带宽语义差异, 为 C++ KDE 实现提供
可复现的基准与红绿判据. 本脚本是 Phase 9 spec 的 probe 前置件, 不参与 C++ 构建.

运行
----
    py -3.12 tests/fixtures/timeseries/probe_sklearn_kde.py
    py -3.12 tests/fixtures/timeseries/probe_sklearn_kde.py --json

判据
----
P1-P8 全部 PASS 时退出码 0; 任一 FAIL 退出码 1.
R 未安装时 P7 标 SKIP (不阻断), JIDT/R 包不在本 probe 范围.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field

import numpy as np
from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity
from statsmodels.nonparametric.bandwidths import (
    bw_normal_reference,
    bw_scott,
    bw_silverman,
)
from statsmodels.nonparametric.kde import KDEUnivariate

SEED = 20261111
N = 200


@dataclass
class Case:
    name: str
    status: str
    detail: str
    values: dict = field(default_factory=dict)


def check(name: str, ok: bool, detail: str, values: dict | None = None) -> Case:
    return Case(name, "PASS" if ok else "FAIL", detail, values or {})


def r_bandwidths(data: np.ndarray) -> tuple[float, float] | None:
    """调用 base R 的 bw.nrd0/bw.nrd; 无 Rscript 时返回 None."""
    rscript = shutil.which("Rscript")
    if rscript is None:
        cand = r"C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
        rscript = cand if os.path.exists(cand) else None
    if rscript is None:
        return None
    fd, path = tempfile.mkstemp(suffix=".csv")
    os.close(fd)
    try:
        np.savetxt(path, data, fmt="%.17g")
        rpath = path.replace("\\", "/")
        code = (
            f'x <- scan("{rpath}", quiet=TRUE); '
            'cat(sprintf("%.17g %.17g", bw.nrd0(x), bw.nrd(x)))'
        )
        out = subprocess.run(
            [rscript, "-e", code], capture_output=True, text=True, timeout=120
        )
        if out.returncode != 0:
            return None
        a, b = out.stdout.strip().split()
        return float(a), float(b)
    finally:
        os.unlink(path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true", help="输出 JSON")
    args = ap.parse_args()

    rng = np.random.default_rng(SEED)
    data = rng.normal(size=N)
    X = data.reshape(-1, 1)
    sd = float(np.std(data, ddof=1))
    n, d = N, 1
    cases: list[Case] = []

    # P1: sklearn 字符串带宽 = 裸因子 (不乘数据尺度)
    sk_scott = float(KernelDensity(bandwidth="scott").fit(X).bandwidth_)
    sk_silver = float(KernelDensity(bandwidth="silverman").fit(X).bandwidth_)
    exp_scott = n ** (-1 / (d + 4))
    exp_silver = (n * (d + 2) / 4) ** (-1 / (d + 4))
    cases.append(check(
        "P1 sklearn 字符串带宽为裸因子",
        abs(sk_scott - exp_scott) < 1e-12 and abs(sk_silver - exp_silver) < 1e-12,
        f"scott={sk_scott:.9f} (exp {exp_scott:.9f}); "
        f"silverman={sk_silver:.9f} (exp {exp_silver:.9f})",
        {"sk_scott": sk_scott, "sk_silverman": sk_silver},
    ))

    # P2: sklearn 带宽不随数据尺度变化 (尺度语义陷阱)
    X10 = (10 * data).reshape(-1, 1)
    sk_scott_10 = float(KernelDensity(bandwidth="scott").fit(X10).bandwidth_)
    cases.append(check(
        "P2 sklearn 带宽对尺度不变 (即不乘 sd)",
        abs(sk_scott_10 - sk_scott) < 1e-15,
        f"sigma=1 -> {sk_scott:.9f}; sigma=10 -> {sk_scott_10:.9f}",
        {"sk_scott_sigma1": sk_scott, "sk_scott_sigma10": sk_scott_10},
    ))

    # P3: scipy 因子无量纲, 实际 h = factor * sd(data)
    sc = gaussian_kde(data)
    sc_factor = float(sc.factor)
    sc_h = sc_factor * sd
    cases.append(check(
        "P3 scipy factor 为无量纲乘子",
        abs(sc_factor - exp_scott) < 1e-12 and abs(sc_h - exp_scott * sd) < 1e-12,
        f"factor={sc_factor:.9f}; sd={sd:.9f}; h={sc_h:.9f}",
        {"scipy_factor": sc_factor, "scipy_h": sc_h, "sd": sd},
    ))

    # P4: statsmodels 三个选择器的绝对 h 公式
    A = min(sd, float(np.subtract(*np.percentile(data, [75, 25]))) / 1.349)
    sm_scott = float(bw_scott(data))
    sm_silver = float(bw_silverman(data))
    sm_nr = float(bw_normal_reference(data))
    cases.append(check(
        "P4 statsmodels 选择器为绝对 h (1.059/0.9/C)",
        abs(sm_scott - 1.059 * A * n ** -0.2) < 1e-12
        and abs(sm_silver - 0.9 * A * n ** -0.2) < 1e-12
        and abs(sm_nr - sm_scott) < 0.02 * sm_scott,
        f"scott={sm_scott:.9f}; silverman={sm_silver:.9f}; normal_ref={sm_nr:.9f}; A={A:.9f}",
        {"sm_scott": sm_scott, "sm_silverman": sm_silver, "sm_normal_ref": sm_nr, "A": A},
    ))

    # P5: statsmodels 文档列出 cv_ml/cv_ls 但实现不支持 (fft 两条路径)
    cv_ok = True
    cv_msgs = []
    for fft in (True, False):
        for bw in ("cv_ml", "cv_ls"):
            try:
                KDEUnivariate(data).fit(bw=bw, fft=fft)
            except ValueError as exc:
                cv_msgs.append(f"{bw}/fft={fft}: {exc}")
            else:
                cv_ok = False
    cases.append(check(
        "P5 statsmodels cv_ml/cv_ls 不可用",
        cv_ok and len(cv_msgs) == 4,
        "; ".join(cv_msgs),
        {"messages": cv_msgs},
    ))

    # P6: 同核同 h 时 sklearn vs scipy 密度一致 (1e-10)
    h = 0.3
    grid = np.array([-1.0, -0.3, 0.0, 0.3, 1.0])
    sk_density = np.exp(
        KernelDensity(kernel="gaussian", bandwidth=h, atol=0.0, rtol=0.0)
        .fit(X)
        .score_samples(grid.reshape(-1, 1))
    )
    sc_density = np.asarray(gaussian_kde(data, bw_method=h / sd).pdf(grid))
    max_diff = float(np.max(np.abs(sk_density - sc_density)))
    rel = float(np.max(np.abs(sk_density - sc_density) / np.maximum(sc_density, 1e-300)))
    cases.append(check(
        "P6 sklearn/scipy 同 h 密度一致",
        max_diff < 1e-10,
        f"h={h}; max|diff|={max_diff:.3e}; max rel={rel:.3e}",
        {"h": h, "max_abs_diff": max_diff, "max_rel_diff": rel,
         "sk_density": sk_density.tolist(), "sc_density": sc_density.tolist()},
    ))

    # P7: R base bw.nrd0/bw.nrd 公式 (无 R 则 SKIP)
    rb = r_bandwidths(data)
    if rb is None:
        cases.append(Case("P7 R bw.nrd0/bw.nrd", "SKIP", "Rscript 不可用"))
    else:
        exp_nrd0 = 0.9 * min(sd, float(np.subtract(*np.percentile(data, [75, 25]))) / 1.34) * n ** -0.2
        exp_nrd = 1.06 * min(sd, float(np.subtract(*np.percentile(data, [75, 25]))) / 1.34) * n ** -0.2
        cases.append(check(
            "P7 R bw.nrd0/bw.nrd 公式",
            abs(rb[0] - exp_nrd0) < 1e-10 and abs(rb[1] - exp_nrd) < 1e-10,
            f"nrd0={rb[0]:.9f} (exp {exp_nrd0:.9f}); nrd={rb[1]:.9f} (exp {exp_nrd:.9f})",
            {"r_nrd0": rb[0], "r_nrd": rb[1], "exp_nrd0": exp_nrd0, "exp_nrd": exp_nrd},
        ))

    # P8: 加权时 scipy 用 neff, sklearn 用名义 n (带宽不随 weight 变化)
    w = rng.random(N)
    w /= w.sum()
    sc_w = gaussian_kde(data, weights=w)
    kd_w = KernelDensity(bandwidth="scott").fit(X, sample_weight=w)
    neff = float(np.sum(w) ** 2 / np.sum(w ** 2))
    cases.append(check(
        "P8 加权带宽: scipy 用 neff, sklearn 用 n",
        abs(float(sc_w.neff) - neff) < 1e-10
        and abs(float(kd_w.bandwidth_) - sk_scott) < 1e-15,
        f"scipy neff={sc_w.neff:.6f} (exp {neff:.6f}); "
        f"sklearn bw={float(kd_w.bandwidth_):.9f} (unweighted {sk_scott:.9f})",
        {"scipy_neff": float(sc_w.neff), "sklearn_bw_weighted": float(kd_w.bandwidth_)},
    ))

    if args.json:
        print(json.dumps(
            {"seed": SEED, "n": N, "cases": [c.__dict__ for c in cases]},
            ensure_ascii=False, indent=2,
        ))
    else:
        for c in cases:
            print(f"[{c.status}] {c.name}: {c.detail}")
        n_fail = sum(c.status == "FAIL" for c in cases)
        n_skip = sum(c.status == "SKIP" for c in cases)
        print(f"\n{len(cases) - n_fail - n_skip} PASS, {n_fail} FAIL, {n_skip} SKIP")
    return 0 if all(c.status != "FAIL" for c in cases) else 1


if __name__ == "__main__":
    sys.exit(main())
