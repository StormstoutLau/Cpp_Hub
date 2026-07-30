#!/usr/bin/env python3
"""
Cpp_Hub Phase 4 LITE - G2 整改项: VaR 模块 Python 交叉验证

对比 Cpp_Hub C++ 实现与 Python (numpy/scipy/statsmodels) 的 VaR 计算:
  1. Historical VaR (Linear/Conservative/Empirical 插值)
  2. Parametric VaR (Normal 分布假设)
  3. Weighted (BRM decay) VaR
  4. MC VaR (确定性种子)

输出 JSON 格式基准值,供 C++ 测试 test_python_cross_validation.cpp 加载验证。

依赖: numpy, scipy, statsmodels (可选, fallback 到 numpy)

用法:
    python cross_validate_var.py > benchmarks_var.json
    python cross_validate_var.py --seed 42 --n-samples 5000

参考:
  - Cpp_Hub: include/cpphub/risk/var/historical_var.hpp
  - Cpp_Hub: include/cpphub/risk/var/parametric_var.hpp
  - Cpp_Hub: include/cpphub/risk/var/mc_var.hpp
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy import stats as st

try:
    import statsmodels.api as sm
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("# WARNING: statsmodels not available, using numpy/scipy only",
          file=sys.stderr)


# ---------------------------------------------------------------------------
# VaR 计算函数 (镜像 C++ HistoricalVaR 实现)
# ---------------------------------------------------------------------------

def historical_var_linear(pnl: List[float], confidence: float = 0.99) -> float:
    """历史 VaR (Linear 插值) - 镜像 C++ QuantileInterpolation::Linear。

    C++ 实现: index = q * (n-1), lo=floor, hi=ceil, frac=index-lo,
              sorted[lo] + frac*(sorted[hi]-sorted[lo])
    注意 C++ 用 losses = -pnl, 然后 var = -quantile * sqrt(horizon)
    """
    losses = np.sort(-np.array(pnl))
    n = len(losses)
    q = 1.0 - confidence
    index = q * (n - 1)
    lo = int(np.floor(index))
    hi = int(np.ceil(index))
    if lo == hi:
        val = losses[lo]
    else:
        frac = index - lo
        val = losses[lo] + frac * (losses[hi] - losses[lo])
    return -val  # horizon=1, sqrt(1)=1


def historical_var_conservative(pnl: List[float], confidence: float = 0.99) -> float:
    """历史 VaR (Conservative 插值) - 镜像 C++ QuantileInterpolation::Conservative。

    C++ 实现: 在 lo/hi 之间选择 |值| 更大的 (即更保守的损失)
    """
    losses = np.sort(-np.array(pnl))
    n = len(losses)
    q = 1.0 - confidence
    index = q * (n - 1)
    lo = int(np.floor(index))
    hi = int(np.ceil(index))
    if lo >= n:
        lo = n - 1
    if hi >= n:
        hi = n - 1
    if lo == hi:
        return -losses[lo]
    if abs(losses[lo]) >= abs(losses[hi]):
        return -losses[lo]
    else:
        return -losses[hi]


def historical_var_empirical(pnl: List[float], confidence: float = 0.99) -> float:
    """历史 VaR (Empirical 插值) - 镜像 C++ QuantileInterpolation::Empirical。

    C++ 实现: idx = int(q * n), sorted[idx]
    """
    losses = np.sort(-np.array(pnl))
    n = len(losses)
    q = 1.0 - confidence
    idx = int(q * n)
    if idx >= n:
        idx = n - 1
    return -losses[idx]


def parametric_var_normal(mean: float, std: float,
                           confidence: float = 0.99) -> float:
    """参数化 VaR (Normal 分布假设) - 镜像 C++ ParametricVaR。

    VaR = -(mean + std * z_alpha), 其中 z_alpha = Φ⁻¹(1-confidence)
    (losses = -pnl, 所以 VaR 是正数表示损失)
    """
    z = st.norm.ppf(1.0 - confidence)
    return -(mean + std * z)


def weighted_var_brm(pnl: List[float], decay: float = 0.99,
                      confidence: float = 0.99) -> float:
    """BRM 衰减加权 VaR - 镜像 C++ HistoricalVaR::weighted_var。

    权重: w_i = (1-decay) * decay^(n-1-i)
    """
    losses = -np.array(pnl)
    n = len(losses)
    weights = np.array([(1.0 - decay) * decay ** (n - 1 - i)
                         for i in range(n)])
    # 按损失排序
    order = np.argsort(losses)
    sorted_losses = losses[order]
    sorted_weights = weights[order]

    cum_weight = 0.0
    total_weight = np.sum(weights)
    q = 1.0 - confidence
    target = q * total_weight

    for i in range(n - 1):
        next_cum = cum_weight + sorted_weights[i]
        if next_cum >= target:
            frac = (target - cum_weight) / sorted_weights[i]
            val = (sorted_losses[i] +
                   frac * (sorted_losses[i + 1] - sorted_losses[i]))
            return -val
        cum_weight = next_cum
    return -sorted_losses[-1]


def mc_var_deterministic(n_paths: int, seed: int = 42,
                          confidence: float = 0.99) -> float:
    """MC VaR (确定性种子) - 使用 numpy MT19937 镜像 C++ MCVaR 行为。

    注意: C++ 用 Philox4x64 RNG,numpy 用 MT19937,数值不会位精确相等,
    但统计量 (均值/方差/VaR) 应在 MC 误差范围内一致 (容差 1e-2)。
    这里仅作为统计一致性参考,不作为位精确基准。
    """
    rng = np.random.default_rng(seed)
    # 模拟 GBM 终值: S_T = S0 * exp((r-σ²/2)T + σ√T Z)
    S0, K, T, r, sigma = 100.0, 100.0, 1.0, 0.05, 0.20
    Z = rng.standard_normal(n_paths)
    ST = S0 * np.exp((r - 0.5 * sigma ** 2) * T + sigma * np.sqrt(T) * Z)
    payoffs = np.maximum(ST - K, 0.0) * np.exp(-r * T)
    pnl = payoffs - np.mean(payoffs)
    losses = -pnl
    losses.sort()
    q = 1.0 - confidence
    index = q * (n_paths - 1)
    lo = int(np.floor(index))
    hi = int(np.ceil(index))
    frac = index - lo
    val = losses[lo] + frac * (losses[hi] - losses[lo])
    return -val


# ---------------------------------------------------------------------------
# 基准数据生成 (确定性,固定种子)
# ---------------------------------------------------------------------------

def generate_pnl_history(n_samples: int = 5000, seed: int = 42) -> List[float]:
    """生成确定性 PnL 历史 (正态分布 + 少量厚尾)"""
    rng = np.random.default_rng(seed)
    # 主体: N(0.001, 0.02) 的日收益
    normal_part = rng.normal(0.001, 0.02, n_samples)
    # 5% 概率出现厚尾 (混合分布)
    tail_mask = rng.random(n_samples) < 0.05
    tail_part = rng.normal(-0.005, 0.06, n_samples)
    pnl = np.where(tail_mask, tail_part, normal_part)
    return pnl.tolist()


def generate_benchmarks(seed: int = 42, n_samples: int = 5000) -> Dict:
    """生成所有 VaR 基准值"""
    pnl = generate_pnl_history(n_samples, seed)

    # 基本统计量
    mean_pnl = float(np.mean(pnl))
    std_pnl = float(np.std(pnl, ddof=0))  # 总体标准差 (与 C++ 一致)

    benchmarks = {
        "metadata": {
            "generator": "cross_validate_var.py",
            "version": "1.0",
            "seed": seed,
            "n_samples": n_samples,
            "statsmodels_available": HAS_STATSMODELS,
            "description": "VaR 交叉验证基准 (Python numpy/scipy/statsmodels)",
            "cpp_reference": [
                "include/cpphub/risk/var/historical_var.hpp",
                "include/cpphub/risk/var/parametric_var.hpp",
                "include/cpphub/risk/var/mc_var.hpp",
            ],
        },
        "inputs": {
            "pnl_history_first10": pnl[:10],  # 供 C++ 测试 sanity check
            "pnl_mean": mean_pnl,
            "pnl_std": std_pnl,
            "pnl_min": float(np.min(pnl)),
            "pnl_max": float(np.max(pnl)),
        },
        "historical_var": {
            "confidence_99": {
                "linear": historical_var_linear(pnl, 0.99),
                "conservative": historical_var_conservative(pnl, 0.99),
                "empirical": historical_var_empirical(pnl, 0.99),
            },
            "confidence_95": {
                "linear": historical_var_linear(pnl, 0.95),
                "conservative": historical_var_conservative(pnl, 0.95),
                "empirical": historical_var_empirical(pnl, 0.95),
            },
        },
        "parametric_var": {
            "confidence_99": parametric_var_normal(mean_pnl, std_pnl, 0.99),
            "confidence_95": parametric_var_normal(mean_pnl, std_pnl, 0.95),
            "z_alpha_99": float(st.norm.ppf(0.01)),
            "z_alpha_95": float(st.norm.ppf(0.05)),
        },
        "weighted_var": {
            "decay_0_99": {
                "confidence_99": weighted_var_brm(pnl, 0.99, 0.99),
                "confidence_95": weighted_var_brm(pnl, 0.99, 0.95),
            },
            "decay_0_95": {
                "confidence_99": weighted_var_brm(pnl, 0.95, 0.99),
            },
        },
        "mc_var": {
            "n_paths_10000": {
                "var_99": mc_var_deterministic(10000, seed, 0.99),
                "var_95": mc_var_deterministic(10000, seed, 0.95),
                "tolerance": 2e-2,  # MC 统计容差
                "note": "MC VaR 统计一致性参考 (不同 RNG,容差 2e-2)",
            },
        },
        "tolerances": {
            "historical_linear": 1e-12,      # 算法相同,应位精确
            "historical_conservative": 1e-12,
            "historical_empirical": 1e-12,
            "parametric": 1e-12,             # scipy.stats.norm.ppf 精确
            "weighted": 1e-10,               # 累加顺序可能微差
            "mc": 2e-2,                      # 不同 RNG,统计一致
        },
    }

    # 如果有 statsmodels,增加 statsmodels 交叉验证
    if HAS_STATSMODELS:
        try:
            # statsmodels 的 quantile regression 或 ETS 可用于 VaR
            # 这里用 statsmodels.stats.descriptivestats 的分位数
            from statsmodels.stats.descriptivestats import describe
            desc = describe(np.array(pnl))
            benchmarks["statsmodels_validation"] = {
                "available": True,
                "description": "statsmodels describe() 交叉验证",
                "mean": float(desc["mean"]),
                "std": float(np.sqrt(desc["var"])),
                "nobs": int(desc["nobs"]),
            }
        except Exception as e:
            benchmarks["statsmodels_validation"] = {
                "available": False,
                "error": str(e),
            }
    else:
        benchmarks["statsmodels_validation"] = {
            "available": False,
            "note": "statsmodels not installed, using numpy/scipy only",
        }

    return benchmarks


# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Cpp_Hub VaR 交叉验证基准生成")
    parser.add_argument("--seed", type=int, default=42,
                        help="随机种子 (default: 42)")
    parser.add_argument("--n-samples", type=int, default=5000,
                        help="PnL 样本数 (default: 5000)")
    parser.add_argument("--output", "-o", default=None,
                        help="输出文件 (default: stdout)")
    args = parser.parse_args()

    benchmarks = generate_benchmarks(args.seed, args.n_samples)

    output = json.dumps(benchmarks, indent=2, default=str)
    if args.output:
        with open(args.output, "w") as f:
            f.write(output)
        print(f"# Written to {args.output}", file=sys.stderr)
    else:
        print(output)


if __name__ == "__main__":
    main()
