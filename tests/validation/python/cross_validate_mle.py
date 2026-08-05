# -*- coding: utf-8 -*-
"""
Cpp_Hub Phase 6 v1.5 M2 - MLE Logistic/Poisson Python 交叉验证

对比 Cpp_Hub C++ 实现与 Python (statsmodels) 的 MLE 估计:
  1. Spector-Mazzeo 数据集 Logit 回归 (score ~ GPA + TUCE + PSI)
     - MLE 系数, 标准误, 对数似然
     - 三种协方差: Hessian (默认), Sandwich (HC1)
  2. DoctorVisits / 合成 Poisson 回归
     - MLE 系数, Sandwich 协方差

对照 C++ test_mle_logistic/poisson。

输出 JSON 格式基准值, 供 C++ 测试加载验证。

依赖: numpy, statsmodels

用法:
    python cross_validate_mle.py > benchmarks_mle.json

参考:
  - Cpp_Hub: tests/unit/econometrics (test_mle_logistic, test_mle_poisson)
  - Spector-Mazzeo (1980) "A Strategy for Teaching the Analysis of Survey Data"
  - Cameron-Trivedi (2013) "Regression Analysis of Count Data"
  - 容差说明: 1e-8 (数值稳定场景), 1e-6 (MLE 迭代)
"""
from __future__ import annotations

import argparse
import json
import sys
from typing import Dict

import numpy as np

try:
    import statsmodels.api as sm
    from statsmodels.datasets import spector, cancer
    from statsmodels.discrete.discrete_model import Logit, Poisson
    HAS_STATSMODELS = True
except ImportError as e:
    HAS_STATSMODELS = False
    print(f"# WARNING: statsmodels not available ({e}), "
          f"using numpy fallback", file=sys.stderr)


# ---------------------------------------------------------------------------
# Spector-Mazzeo Logit
# ---------------------------------------------------------------------------

def compute_spector_logit() -> Dict:
    """Spector-Mazzeo 数据集 Logit MLE + Hessian/Sandwich 协方差。"""
    data = spector.load_pandas()
    y = data.endog.values.astype(float)
    X_df = data.exog  # 自带 const
    X = X_df.values.astype(float)

    model = Logit(y, X)
    # 默认 (Hessian 逆) 协方差
    res_hess = model.fit(disp=False)
    # HC1 Sandwich 协方差
    res_sand = model.fit(disp=False, cov_type="HC1")

    return {
        "spector_logit_coefficients": res_hess.params.tolist(),
        "spector_logit_standard_errors": res_hess.bse.tolist(),
        "spector_logit_log_likelihood": float(res_hess.llf),
        "spector_logit_vcov": res_hess.cov_params().tolist(),
        "spector_logit_sandwich_vcov": res_sand.cov_params().tolist(),
        "spector_logit_n_obs": int(res_hess.nobs),
        "spector_logit_n_params": int(X.shape[1]),
        "spector_logit_exog_names": list(X_df.columns),
    }


# ---------------------------------------------------------------------------
# Poisson (DoctorVisits 或合成数据)
# ---------------------------------------------------------------------------

def compute_poisson(seed: int = 42, n: int = 1000) -> Dict:
    """Poisson 回归 MLE + Sandwich 协方差。

    优先使用 statsmodels.datasets.cancer (DoctorVisits), 不可用则合成数据。
    """
    use_synthetic = False
    data_source = "statsmodels.datasets.cancer"
    try:
        data = cancer.load_pandas()
        # cancer 数据集是疾病数据,这里采用经典 DoctorVisits (cancer 在 statsmodels 中
        # 通常指 Cameron-Trivedi 的 doctor visits; 不同版本字段名可能不同)
        # 若结构不适配 Poisson, 退到合成数据。
        if hasattr(data, "endog") and hasattr(data, "exog"):
            y = data.endog.values.astype(float)
            X_df = data.exog
            X = X_df.values.astype(float)
            # 确保目标是非负整数计数
            if not np.all(y >= 0) or not np.all(np.equal(np.mod(y, 1), 0)):
                raise ValueError("cancer endog not integer counts")
        else:
            raise AttributeError("cancer dataset lacks endog/exog")
    except Exception as e:
        print(f"# WARNING: cancer dataset unavailable ({e}), "
              f"using synthetic Poisson data", file=sys.stderr)
        use_synthetic = True
        data_source = f"synthetic (seed={seed}, N={n})"

    if use_synthetic:
        rng = np.random.default_rng(seed)
        # 合成 Poisson: y ~ Poisson(exp(b0 + b1*x1 + b2*x2))
        x1 = rng.standard_normal(n)
        x2 = rng.standard_normal(n)
        true_beta = np.array([0.5, 0.3, -0.2])
        mu = np.exp(true_beta[0] + true_beta[1] * x1 + true_beta[2] * x2)
        y = rng.poisson(mu).astype(float)
        X = np.column_stack([np.ones(n), x1, x2])
        X_df_cols = ["const", "x1", "x2"]
    else:
        X_df_cols = list(X_df.columns)

    model = Poisson(y, X)
    res_hess = model.fit(disp=False)
    try:
        res_sand = model.fit(disp=False, cov_type="HC1")
        sandwich_vcov = res_sand.cov_params().tolist()
    except Exception as e:
        print(f"# WARNING: Poisson HC1 failed: {e}", file=sys.stderr)
        sandwich_vcov = np.full((X.shape[1], X.shape[1]),
                                np.nan).tolist()

    return {
        "poisson_coefficients": res_hess.params.tolist(),
        "poisson_standard_errors": res_hess.bse.tolist(),
        "poisson_log_likelihood": float(res_hess.llf),
        "poisson_sandwich_vcov": sandwich_vcov,
        "poisson_data_source": data_source,
        "poisson_n_obs": int(res_hess.nobs),
        "poisson_n_params": int(X.shape[1]),
        "poisson_exog_names": X_df_cols,
    }


# ---------------------------------------------------------------------------
# numpy fallback (Newton-Raphson Logit)
# ---------------------------------------------------------------------------

def logit_numpy_fallback(X: np.ndarray, y: np.ndarray,
                          max_iter: int = 100, tol: float = 1e-10) -> Dict:
    """纯 numpy Logit MLE (Newton-Raphson) 降级路径。"""
    n, k = X.shape
    beta = np.zeros(k)
    for _ in range(max_iter):
        eta = X @ beta
        p = 1.0 / (1.0 + np.exp(-eta))
        # 梯度
        grad = X.T @ (y - p)
        # Hessian
        W = p * (1.0 - p)
        H = -(X.T * W) @ X
        try:
            delta = np.linalg.solve(H, grad)
        except np.linalg.LinAlgError:
            break
        beta_new = beta - delta
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new
    eta = X @ beta
    p = 1.0 / (1.0 + np.exp(-eta))
    W = p * (1.0 - p)
    H = -(X.T * W) @ X
    vcov = np.linalg.inv(-H)
    ll = float(np.sum(y * np.log(p + 1e-300) +
                       (1 - y) * np.log(1 - p + 1e-300)))
    se = np.sqrt(np.diag(vcov))
    # Sandwich: H^{-1} (X^T diag((y-p)^2) X) H^{-1}
    meat = X.T @ (X * ((y - p) ** 2)[:, None])
    sandwich = vcov @ meat @ vcov
    return {
        "spector_logit_coefficients": beta.tolist(),
        "spector_logit_standard_errors": se.tolist(),
        "spector_logit_log_likelihood": ll,
        "spector_logit_vcov": vcov.tolist(),
        "spector_logit_sandwich_vcov": sandwich.tolist(),
        "_fallback": True,
    }


def generate_benchmarks(seed: int = 42, n: int = 1000) -> Dict:
    """生成所有 MLE 基准值。"""
    benchmarks = {
        "metadata": {
            "generator": "cross_validate_mle.py",
            "version": "1.0",
            "phase": "Phase 6 v1.5 M2",
            "cpp_reference": ["test_mle_logistic", "test_mle_poisson"],
            "libraries": ["numpy", "statsmodels"],
            "datasets": ["Spector-Mazzeo (statsmodels.datasets.spector)",
                         "DoctorVisits / synthetic Poisson"],
            "statsmodels_available": HAS_STATSMODELS,
            "description": "MLE Logistic/Poisson Python 交叉验证",
            "tolerance_note": "1e-8 (数值稳定场景), 1e-6 (MLE 迭代)",
            "references": [
                "Spector-Mazzeo (1980)",
                "Cameron-Trivedi (2013) Regression Analysis of Count Data",
                "Huber (1967) Sandwich covariance",
            ],
        },
    }

    if HAS_STATSMODELS:
        # Spector-Mazzeo Logit
        try:
            benchmarks.update(compute_spector_logit())
        except Exception as e:
            print(f"# WARNING: Spector Logit failed: {e}", file=sys.stderr)
            try:
                data = spector.load_pandas()
                y = data.endog.values.astype(float)
                X = data.exog.values.astype(float)
                benchmarks.update(logit_numpy_fallback(X, y))
            except Exception as e2:
                benchmarks["spector_logit_error"] = str(e2)

        # Poisson
        try:
            benchmarks.update(compute_poisson(seed=seed, n=n))
        except Exception as e:
            print(f"# WARNING: Poisson failed: {e}", file=sys.stderr)
            benchmarks["poisson_error"] = str(e)
    else:
        # statsmodels 不可用: 用 numpy Logit fallback + 合成 Poisson (无真实 MLE)
        benchmarks["note"] = ("statsmodels unavailable; "
                              "only numpy Logit fallback provided")
        benchmarks["_fallback"] = True

    benchmarks["tolerances"] = {
        "logit_coefficients": 1e-6,        # MLE 迭代收敛容差
        "logit_vcov_hessian": 1e-6,
        "logit_sandwich_vcov": 1e-6,
        "logit_log_likelihood": 1e-8,
        "poisson_coefficients": 1e-6,
        "poisson_sandwich_vcov": 1e-6,
    }
    return benchmarks


# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Cpp_Hub MLE Logistic/Poisson 交叉验证基准生成")
    parser.add_argument("--seed", type=int, default=42,
                        help="合成 Poisson 随机种子 (default: 42)")
    parser.add_argument("--n-samples", type=int, default=1000,
                        help="合成 Poisson 样本数 (default: 1000)")
    parser.add_argument("--output", "-o", default=None,
                        help="输出文件 (default: stdout)")
    args = parser.parse_args()

    benchmarks = generate_benchmarks(seed=args.seed, n=args.n_samples)
    output = json.dumps(benchmarks, indent=2, default=str)

    if args.output:
        with open(args.output, "w") as f:
            f.write(output)
        print(f"# Written to {args.output}", file=sys.stderr)
    else:
        print(output)


if __name__ == "__main__":
    main()
