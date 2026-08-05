# -*- coding: utf-8 -*-
"""
Cpp_Hub Phase 6 v1.5 M1 - OLS + HC0-3 异方差稳健协方差 Python 交叉验证

对比 Cpp_Hub C++ 实现与 Python (numpy + statsmodels) 的 OLS 估计:
  1. OLS 系数估计 (Longley 数据集)
  2. HC0-HC3 异方差稳健协方差矩阵
  3. R-squared, adjusted R-squared
  4. F-statistic

对照 C++ test_ols_hc。

输出 JSON 格式基准值, 供 C++ 测试加载验证。

依赖: numpy, statsmodels

用法:
    python cross_validate_ols.py > benchmarks_ols.json

参考:
  - Cpp_Hub: tests/unit/econometrics (test_ols_hc)
  - Longley (1967) "An Appraisal of Least-Squares Programs..."
  - White (1980) HC0; MacKinnon-White (1985) HC1/HC2/HC3
  - 容差说明: 1e-8 (数值稳定场景)
"""
from __future__ import annotations

import argparse
import json
import sys
from typing import Dict

import numpy as np

try:
    import statsmodels.api as sm
    from statsmodels.datasets import longley
    HAS_STATSMODELS = True
except ImportError as e:
    HAS_STATSMODELS = False
    print(f"# WARNING: statsmodels not available ({e}), "
          f"using numpy fallback", file=sys.stderr)


# ---------------------------------------------------------------------------
# numpy fallback OLS (近似, 无 HC1/HC2/HC3 完整实现)
# ---------------------------------------------------------------------------

def ols_numpy_fallback(X: np.ndarray, y: np.ndarray) -> Dict:
    """纯 numpy OLS (供 statsmodels 不可用时降级)。"""
    n, k = X.shape
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    rss = float(resid @ resid)
    y_mean = float(np.mean(y))
    tss = float(np.sum((y - y_mean) ** 2))
    r2 = 1.0 - rss / tss
    adj_r2 = 1.0 - (1.0 - r2) * (n - 1) / (n - k)
    # HC0: X^T diag(r^2) X 逆 * (X^T X)^-1
    meat = X.T @ (X * (resid ** 2)[:, None])
    hc0 = XtX_inv @ meat @ XtX_inv
    # HC1: HC0 * n/(n-k)
    hc1 = hc0 * n / (n - k)
    f_stat = (r2 / (k - 1)) / ((1.0 - r2) / (n - k))
    return {
        "longley_coefficients": beta.tolist(),
        "longley_hc0_vcov": hc0.tolist(),
        "longley_hc1_vcov": hc1.tolist(),
        "longley_hc2_vcov": hc0.tolist(),  # 近似
        "longley_hc3_vcov": hc0.tolist(),  # 近似
        "longley_r_squared": float(r2),
        "longley_adj_r_squared": float(adj_r2),
        "longley_f_stat": float(f_stat),
        "longley_resid_ss": rss,
        "_fallback": True,
    }


# ---------------------------------------------------------------------------
# statsmodels 主路径
# ---------------------------------------------------------------------------

def compute_longley_ols() -> Dict:
    """Longley 数据集 OLS + HC0-HC3。"""
    data = longley.load_pandas()
    # 模型: Employed ~ GNP.deflator + GNP + Unemployed + Armed.Forces + Population + Year
    y = data.endog.values.astype(float)
    X_df = data.exog  # statsmodels Longley 自带常数列 (const)
    X = X_df.values.astype(float)

    # 拟合 OLS (默认 OLS 协方差)
    model = sm.OLS(y, X)
    res_ols = model.fit()

    beta = res_ols.params
    r2 = res_ols.rsquared
    adj_r2 = res_ols.rsquared_adj
    f_stat = res_ols.fvalue

    # HC0-HC3 协方差
    hc_results = {}
    for hc in ["HC0", "HC1", "HC2", "HC3"]:
        try:
            res_hc = model.fit(cov_type=hc)
            hc_results[hc] = res_hc.cov_params()
        except Exception as e:
            print(f"# WARNING: {hc} fit failed: {e}", file=sys.stderr)
            hc_results[hc] = np.full((X.shape[1], X.shape[1]), np.nan)

    return {
        "longley_coefficients": beta.tolist(),
        "longley_hc0_vcov": hc_results["HC0"].tolist(),
        "longley_hc1_vcov": hc_results["HC1"].tolist(),
        "longley_hc2_vcov": hc_results["HC2"].tolist(),
        "longley_hc3_vcov": hc_results["HC3"].tolist(),
        "longley_r_squared": float(r2),
        "longley_adj_r_squared": float(adj_r2),
        "longley_f_stat": float(f_stat),
        "longley_resid_ss": float(res_ols.ssr),
        "longley_n_obs": int(res_ols.nobs),
        "longley_n_params": int(X.shape[1]),
        "longley_exog_names": list(X_df.columns),
        "_fallback": False,
    }


def generate_benchmarks() -> Dict:
    """生成所有 OLS 基准值。"""
    benchmarks = {
        "metadata": {
            "generator": "cross_validate_ols.py",
            "version": "1.0",
            "phase": "Phase 6 v1.5 M1",
            "cpp_reference": "test_ols_hc",
            "libraries": ["numpy", "statsmodels"],
            "dataset": "Longley (statsmodels.datasets.longley)",
            "model": "Employed ~ GNP.deflator + GNP + Unemployed + Armed.Forces + Population + Year",
            "statsmodels_available": HAS_STATSMODELS,
            "description": "OLS + HC0-3 异方差稳健协方差 Python 交叉验证",
            "tolerance_note": "1e-8 (数值稳定场景)",
            "references": [
                "Longley (1967) An Appraisal of Least-Squares Programs...",
                "White (1980) Heteroskedasticity-consistent covariance (HC0)",
                "MacKinnon-White (1985) HC1/HC2/HC3",
            ],
        },
    }

    if HAS_STATSMODELS:
        try:
            ols_result = compute_longley_ols()
            benchmarks.update(ols_result)
        except Exception as e:
            print(f"# WARNING: statsmodels Longley OLS failed: {e}",
                  file=sys.stderr)
            # 降级到 numpy
            data = longley.load_pandas() if "longley" in dir() else None
            if data is not None:
                y = data.endog.values.astype(float)
                X = data.exog.values.astype(float)
                benchmarks.update(ols_numpy_fallback(X, y))
            benchmarks["error"] = str(e)
    else:
        # numpy 降级路径 (硬编码 Longley 数据,因 statsmodels 不可用无法加载)
        # 这里仅填充占位,真实降级需要 C++ 端提供数据
        benchmarks["longley_coefficients"] = []
        benchmarks["note"] = ("statsmodels not available; "
                              "Longley data not loaded in fallback mode.")
        benchmarks["_fallback"] = True

    benchmarks["tolerances"] = {
        "coefficients": 1e-8,        # OLS 解析解,数值稳定
        "hc_vcov": 1e-8,             # HC0-HC3 公式直接,数值稳定
        "r_squared": 1e-10,          # 确定性公式
        "f_statistic": 1e-8,
    }
    return benchmarks


# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Cpp_Hub OLS + HC0-3 交叉验证基准生成")
    parser.add_argument("--output", "-o", default=None,
                        help="输出文件 (default: stdout)")
    args = parser.parse_args()

    benchmarks = generate_benchmarks()
    output = json.dumps(benchmarks, indent=2, default=str)

    if args.output:
        with open(args.output, "w") as f:
            f.write(output)
        print(f"# Written to {args.output}", file=sys.stderr)
    else:
        print(output)


if __name__ == "__main__":
    main()
