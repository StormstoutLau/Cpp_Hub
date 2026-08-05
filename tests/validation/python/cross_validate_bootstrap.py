# -*- coding: utf-8 -*-
"""
Cpp_Hub Phase 6 v1.5 M4 - Bootstrap Python 交叉验证

对比 Cpp_Hub C++ 实现与 Python (arch / numpy) 的 Bootstrap 估计:
  1. Pearson 相关系数 r (法学院数据)
  2. 配对 Bootstrap CI (B=999)
  3. Wild Bootstrap (Rademacher)
  4. Block Bootstrap (arch.bootstrap.CircularBlockBootstrap, 可选)
  5. BCa CI (arch.bootstrap.IIDBootstrap.conf_int, 可选)

对照 C++ test_paired/wild/block/cluster_bootstrap。

输出 JSON 格式基准值, 供 C++ 测试加载验证。

依赖:
  - numpy (必需)
  - arch (可选, 提供 Block / BCa; 不可用时用纯 numpy 实现配对/Wild)

用法:
    python cross_validate_bootstrap.py > benchmarks_bootstrap.json

参考:
  - Cpp_Hub: tests/unit/econometrics (test_paired/wild/block/cluster_bootstrap)
  - Efron-Tibshirani (1993) "An Introduction to the Bootstrap"
  - Liu (1988) Wild Bootstrap
  - Politis-Romano (1994) Circular Block Bootstrap
  - 容差说明: 1e-6 (确定性种子配对/Wild), 1e-2 (Block/BCa 优化器)
"""
from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy import stats as st

try:
    from arch.bootstrap import (CircularBlockBootstrap, IIDBootstrap,
                                 StationaryBootstrap)
    HAS_ARCH = True
except ImportError as e:
    HAS_ARCH = False
    print(f"# WARNING: arch not available ({e}); "
          f"using pure numpy for paired/wild bootstrap", file=sys.stderr)


# ---------------------------------------------------------------------------
# 法学院数据 (Efron-Tibshirani 1993, 硬编码)
# ---------------------------------------------------------------------------

LSAT = np.array([576, 635, 558, 578, 666, 580, 555, 661, 651, 605,
                  653, 575, 545, 572, 594], dtype=float)
GPA = np.array([3.39, 3.30, 2.81, 3.03, 3.44, 3.07, 3.00, 3.43, 3.36,
                 3.13, 3.12, 2.74, 2.76, 2.88, 2.96], dtype=float)
N_LAW = 15


def pearson_r(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson 相关系数。"""
    return float(st.pearsonr(x, y)[0])


# ---------------------------------------------------------------------------
# 配对 Bootstrap (纯 numpy)
# ---------------------------------------------------------------------------

def paired_bootstrap(x: np.ndarray, y: np.ndarray, B: int = 999,
                      seed: int = 42) -> Dict:
    """配对 Bootstrap CI (镜像 C++ test_paired_bootstrap)。

    用法: np.random.choice(15, 15, replace=True) 生成索引, 在 (x, y) 配对样本
    上重采样, 计算 Pearson r 的百分位 CI。
    """
    rng = np.random.default_rng(seed)
    n = len(x)
    r_stars = np.empty(B)
    for b in range(B):
        idx = rng.choice(n, n, replace=True)
        r_stars[b] = st.pearsonr(x[idx], y[idx])[0]
    lo, hi = np.percentile(r_stars, [2.5, 97.5])
    return {
        "paired_bootstrap_ci": [float(lo), float(hi)],
        "paired_bootstrap_mean": float(np.mean(r_stars)),
        "paired_bootstrap_std": float(np.std(r_stars, ddof=1)),
        "paired_bootstrap_B": int(B),
        "paired_bootstrap_seed": int(seed),
        "paired_bootstrap_first10_indices": rng.choice(n, n, replace=True)[:10].tolist(),
    }


# ---------------------------------------------------------------------------
# Wild Bootstrap (Rademacher, 纯 numpy)
# ---------------------------------------------------------------------------

def wild_bootstrap(x: np.ndarray, y: np.ndarray, B: int = 999,
                    seed: int = 42) -> Dict:
    """Wild Bootstrap (Rademacher) - 镜像 C++ test_wild_bootstrap。

    模型: y = β0 + β1*x + ε
    y* = fitted + v * resid, v ∈ {-1, +1} (Rademacher)
    在每次重采样中重新拟合 OLS, 计算 β1* 或 r*。
    这里采用 r* (相关系数) 作为统计量, 与配对 Bootstrap 对照。
    """
    rng = np.random.default_rng(seed)
    n = len(x)
    # 原始 OLS 拟合
    X = np.column_stack([np.ones(n), x])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    fitted = X @ beta
    resid = y - fitted

    r_stars = np.empty(B)
    beta1_stars = np.empty(B)
    for b in range(B):
        v = rng.choice([-1.0, 1.0], size=n, replace=True)
        y_star = fitted + v * resid
        r_stars[b] = st.pearsonr(x, y_star)[0]
        beta1_stars[b] = np.linalg.lstsq(X, y_star, rcond=None)[0][1]
    lo_r, hi_r = np.percentile(r_stars, [2.5, 97.5])
    lo_b, hi_b = np.percentile(beta1_stars, [2.5, 97.5])
    return {
        "wild_bootstrap_ci": [float(lo_r), float(hi_r)],
        "wild_bootstrap_ci_beta1": [float(lo_b), float(hi_b)],
        "wild_bootstrap_mean": float(np.mean(r_stars)),
        "wild_bootstrap_std": float(np.std(r_stars, ddof=1)),
        "wild_bootstrap_B": int(B),
        "wild_bootstrap_seed": int(seed),
        "wild_bootstrap_distribution": "Rademacher",
    }


# ---------------------------------------------------------------------------
# Block Bootstrap (arch, 可选)
# ---------------------------------------------------------------------------

def block_bootstrap_arch(x: np.ndarray, y: np.ndarray, B: int = 999,
                          block_size: int = 3, seed: int = 42) -> Dict:
    """Block Bootstrap (arch.bootstrap.CircularBlockBootstrap)。

    镜像 C++ test_block_bootstrap。
    """
    if not HAS_ARCH:
        return {
            "block_bootstrap_available": False,
            "block_bootstrap_note": "arch package not available; skipped",
        }

    def stat_func(idx):
        # arch 传递的是索引
        return st.pearsonr(x[idx], y[idx])[0]

    bs = CircularBlockBootstrap(block_size, x, y, seed=seed)
    results = bs.apply(stat_func, B)
    r_stars = results.flatten()
    lo, hi = np.percentile(r_stars, [2.5, 97.5])
    return {
        "block_bootstrap_ci": [float(lo), float(hi)],
        "block_bootstrap_mean": float(np.mean(r_stars)),
        "block_bootstrap_std": float(np.std(r_stars, ddof=1)),
        "block_bootstrap_B": int(B),
        "block_bootstrap_block_size": int(block_size),
        "block_bootstrap_method": "arch.bootstrap.CircularBlockBootstrap",
        "block_bootstrap_available": True,
    }


# ---------------------------------------------------------------------------
# BCa CI (arch, 可选)
# ---------------------------------------------------------------------------

def bca_bootstrap_arch(x: np.ndarray, y: np.ndarray, B: int = 999,
                        seed: int = 42) -> Dict:
    """BCa CI (arch.bootstrap.IIDBootstrap.conf_int)。

    镜像 C++ 实现 (若有)。BCa = Bias-Corrected and Accelerated。
    """
    if not HAS_ARCH:
        return {
            "bca_ci_available": False,
            "bca_ci_note": "arch package not available; skipped",
        }

    def stat_func(*args):
        # arch 传参风格: IIDBootstrap 把每个数组作为位置参数
        # 新版 arch: stat_func(x, y) 接收多个位置参数
        if len(args) == 2:
            xx, yy = args
        else:
            xx, yy = args[0]
        return st.pearsonr(xx, yy)[0]

    bs = IIDBootstrap(x, y, seed=seed)
    try:
        ci = bs.conf_int(stat_func, 1000, method="bca",
                          size=0.05)
        # ci 是 2 x ? 矩阵: [lower, upper]
        ci_arr = np.asarray(ci)
        if ci_arr.ndim == 2:
            lo, hi = float(ci_arr[0, 0]), float(ci_arr[1, 0])
        else:
            lo, hi = float(ci_arr[0]), float(ci_arr[1])
        return {
            "bca_ci": [lo, hi],
            "bca_ci_method": "arch.bootstrap.IIDBootstrap.conf_int (bca)",
            "bca_ci_available": True,
            "bca_ci_B": 1000,
        }
    except Exception as e:
        print(f"# WARNING: arch BCa failed: {e}", file=sys.stderr)
        return {
            "bca_ci_available": False,
            "bca_ci_error": str(e),
        }


# ---------------------------------------------------------------------------
# numpy 兜底 Block Bootstrap (近似)
# ---------------------------------------------------------------------------

def block_bootstrap_numpy(x: np.ndarray, y: np.ndarray, B: int = 999,
                            block_size: int = 3, seed: int = 42) -> Dict:
    """纯 numpy Circular Block Bootstrap (arch 不可用时的兜底)。"""
    rng = np.random.default_rng(seed)
    n = len(x)
    n_blocks = int(np.ceil(n / block_size))
    r_stars = np.empty(B)
    for b in range(B):
        idx = np.empty(n, dtype=int)
        pos = 0
        for _ in range(n_blocks):
            start = rng.integers(0, n)
            for j in range(block_size):
                if pos < n:
                    idx[pos] = (start + j) % n
                    pos += 1
        r_stars[b] = st.pearsonr(x[idx], y[idx])[0]
    lo, hi = np.percentile(r_stars, [2.5, 97.5])
    return {
        "block_bootstrap_ci": [float(lo), float(hi)],
        "block_bootstrap_mean": float(np.mean(r_stars)),
        "block_bootstrap_std": float(np.std(r_stars, ddof=1)),
        "block_bootstrap_B": int(B),
        "block_bootstrap_block_size": int(block_size),
        "block_bootstrap_method": "numpy circular block (fallback)",
        "block_bootstrap_available": True,
    }


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------

def generate_benchmarks(B: int = 999, seed: int = 42,
                          block_size: int = 3) -> Dict:
    """生成所有 Bootstrap 基准值。"""
    r_hat = pearson_r(LSAT, GPA)

    benchmarks = {
        "metadata": {
            "generator": "cross_validate_bootstrap.py",
            "version": "1.0",
            "phase": "Phase 6 v1.5 M4",
            "cpp_reference": [
                "test_paired_bootstrap",
                "test_wild_bootstrap",
                "test_block_bootstrap",
                "test_cluster_bootstrap",
            ],
            "libraries": ["numpy", "arch"],
            "arch_available": HAS_ARCH,
            "dataset": "Law School (Efron-Tibshirani 1993, hardcoded)",
            "n_samples": int(N_LAW),
            "B": int(B),
            "seed": int(seed),
            "description": "Bootstrap Python 交叉验证",
            "tolerance_note": "1e-6 (确定性种子配对/Wild), "
                              "1e-2 (Block/BCa 优化器)",
            "references": [
                "Efron-Tibshirani (1993) An Introduction to the Bootstrap",
                "Liu (1988) Wild Bootstrap",
                "Politis-Romano (1994) Circular Block Bootstrap",
                "Efron (1987) BCa CI",
            ],
        },
        "law_school_data": {
            "LSAT": LSAT.tolist(),
            "GPA": GPA.tolist(),
        },
        "law_school_pearson_r": float(r_hat),
    }

    # 配对 Bootstrap
    try:
        benchmarks.update(paired_bootstrap(LSAT, GPA, B=B, seed=seed))
    except Exception as e:
        benchmarks["paired_bootstrap_error"] = str(e)

    # Wild Bootstrap
    try:
        benchmarks.update(wild_bootstrap(LSAT, GPA, B=B, seed=seed))
    except Exception as e:
        benchmarks["wild_bootstrap_error"] = str(e)

    # Block Bootstrap (优先 arch, 否则 numpy 兜底)
    if HAS_ARCH:
        try:
            benchmarks.update(block_bootstrap_arch(LSAT, GPA, B=B,
                                                     block_size=block_size,
                                                     seed=seed))
        except Exception as e:
            print(f"# WARNING: arch block bootstrap failed ({e}); "
                  f"falling back to numpy", file=sys.stderr)
            benchmarks.update(block_bootstrap_numpy(LSAT, GPA, B=B,
                                                     block_size=block_size,
                                                     seed=seed))
    else:
        benchmarks.update(block_bootstrap_numpy(LSAT, GPA, B=B,
                                                  block_size=block_size,
                                                  seed=seed))

    # BCa CI (仅 arch 可用时)
    if HAS_ARCH:
        benchmarks.update(bca_bootstrap_arch(LSAT, GPA, B=B, seed=seed))
    else:
        benchmarks["bca_ci_available"] = False
        benchmarks["bca_ci_note"] = ("arch package not available; "
                                       "BCa not computed (numpy fallback "
                                       "for BCa is non-trivial)")

    benchmarks["tolerances"] = {
        "pearson_r": 1e-12,                # 解析公式
        "paired_bootstrap_ci": 1e-6,       # 确定性种子
        "wild_bootstrap_ci": 1e-6,         # 确定性种子
        "block_bootstrap_ci": 1e-2,        # 不同 RNG / 循环顺序
        "bca_ci": 1e-2,                    # BCa 涉及 jackknife
    }
    return benchmarks


# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Cpp_Hub Bootstrap 交叉验证基准生成")
    parser.add_argument("--B", type=int, default=999,
                        help="Bootstrap 重复次数 (default: 999)")
    parser.add_argument("--seed", type=int, default=42,
                        help="随机种子 (default: 42)")
    parser.add_argument("--block-size", type=int, default=3,
                        help="Block Bootstrap 块大小 (default: 3)")
    parser.add_argument("--output", "-o", default=None,
                        help="输出文件 (default: stdout)")
    args = parser.parse_args()

    benchmarks = generate_benchmarks(B=args.B, seed=args.seed,
                                       block_size=args.block_size)
    output = json.dumps(benchmarks, indent=2, default=str)

    if args.output:
        with open(args.output, "w") as f:
            f.write(output)
        print(f"# Written to {args.output}", file=sys.stderr)
    else:
        print(output)


if __name__ == "__main__":
    main()
