# -*- coding: utf-8 -*-
"""
Cpp_Hub Phase 6 v1.5 M3 - GMM + CUE Python 交叉验证

对比 Cpp_Hub C++ 实现与 Python (linearmodels / statsmodels / numpy) 的 GMM 估计:
  1. 两步 GMM (Hansen 1982)
  2. CUE (Continuously Updating Estimator, Hansen-Heaton-Yaron 1996)
  3. J 检验统计量与 p 值

对照 C++ test_gmm_two_step/cue_iterated。

输出 JSON 格式基准值, 供 C++ 测试加载验证。

依赖 (按优先级):
  1. linearmodels (IVGMM, IVGMMCUE) - 首选
  2. statsmodels.sandbox.regression.gmm (IVGMM) - 备选
  3. numpy 手动实现两步 GMM - 兜底

用法:
    python cross_validate_gmm.py > benchmarks_gmm.json

参考:
  - Cpp_Hub: tests/unit/econometrics (test_gmm_two_step, test_gmm_cue_iterated)
  - Hansen (1982) "Large Sample Properties of GMM Estimators"
  - Hansen-Heaton-Yaron (1996) "Finite-Sample Properties of CUE"
  - 容差说明: 1e-6 (GMM 迭代收敛), 1e-8 (J 统计量公式)
"""
from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, Tuple

import numpy as np
from scipy import stats as st


# ---------------------------------------------------------------------------
# 依赖探测 (linearmodels > statsmodels.gmm > numpy 手动)
# ---------------------------------------------------------------------------

HAS_LINEARMODELS = False
HAS_SM_GMM = False

try:
    from linearmodels.iv import IVGMM, IVGMMCUE
    HAS_LINEARMODELS = True
except ImportError as e:
    print(f"# WARNING: linearmodels not available ({e})", file=sys.stderr)

if not HAS_LINEARMODELS:
    try:
        from statsmodels.sandbox.regression.gmm import IVGMM as SM_IVGMM
        HAS_SM_GMM = True
    except ImportError as e:
        print(f"# WARNING: statsmodels.gmm not available ({e}); "
              f"will use numpy manual GMM", file=sys.stderr)


# ---------------------------------------------------------------------------
# 合成 IV 数据 (set.seed(42), N=500)
# ---------------------------------------------------------------------------

def generate_iv_data(seed: int = 42, n: int = 500) -> Dict:
    """合成 IV 数据 (镜像 C++ test_gmm_two_step 数据生成)。

    结构:
      内生: x = 0.5*z + 0.3*error_x + noise
      y = 1.0 + 2.0*x + error_y
      工具: Z = [1, z, z2] (z2 是额外工具)
    其中 (error_x, error_y) 服从联合正态, 使 x 内生 (与 y 的扰动相关)。
    """
    rng = np.random.default_rng(seed)
    # 联合误差: error_x 和 error_y 相关
    mean = [0.0, 0.0]
    cov = [[1.0, 0.5], [0.5, 1.0]]
    errors = rng.multivariate_normal(mean, cov, size=n)
    error_x = errors[:, 0]
    error_y = errors[:, 1]

    # 工具变量 z 和额外工具 z2
    z = rng.standard_normal(n)
    z2 = rng.standard_normal(n)
    # 内生变量 x
    noise = rng.standard_normal(n) * 0.1
    x = 0.5 * z + 0.3 * error_x + noise
    # 结果 y
    y = 1.0 + 2.0 * x + error_y

    return {
        "y": y, "x": x, "z": z, "z2": z2,
        "n": n, "seed": seed,
    }


# ---------------------------------------------------------------------------
# linearmodels 路径
# ---------------------------------------------------------------------------

def _to_linearmodels_df(data: Dict):
    """构造 linearmodels 所需的 DataFrame。"""
    import pandas as pd
    df = pd.DataFrame({
        "y": data["y"],
        "x": data["x"],
        "z": data["z"],
        "z2": data["z2"],
        "const": 1.0,
    })
    return df


def compute_gmm_linearmodels(data: Dict) -> Dict:
    """使用 linearmodels.IVGMM / IVGMMCUE。"""
    df = _to_linearmodels_df(data)

    # 公式: y ~ [x] / [const, z, z2]  (1 个内生, 3 个工具, 过度识别 2)
    formula = "y ~ 1 + [x ~ 1 + z + z2]"
    # 注意: linearmodels 公式语法 [endog ~ instruments]
    # 实际正确语法: y ~ [x ~ z + z2] (const 自动加)
    formula = "y ~ 1 + [x ~ z + z2]"

    result = {}

    # 两步 GMM
    try:
        from linearmodels.iv import IVGMM
        mod = IVGMM.from_formula(formula, df)
        res = mod.fit(cov_type="robust", iter_limit=2, display=False)
        result["gmm_twostep_coefficients"] = res.params.tolist()
        result["gmm_twostep_method"] = "linearmodels.IVGMM"
        # J 统计量
        j_stat = float(res.j_stat.stat)
        j_pvalue = float(res.j_stat.pval)
        result["gmm_j_statistic"] = j_stat
        result["gmm_j_pvalue"] = j_pvalue
    except Exception as e:
        print(f"# WARNING: linearmodels IVGMM failed: {e}", file=sys.stderr)
        result["gmm_twostep_error"] = str(e)

    # CUE
    try:
        from linearmodels.iv import IVGMMCUE
        mod_cue = IVGMMCUE.from_formula(formula, df)
        res_cue = mod_cue.fit(cov_type="robust", display=False)
        result["gmm_cue_coefficients"] = res_cue.params.tolist()
        result["gmm_cue_method"] = "linearmodels.IVGMMCUE"
        # 若两步 GMM 失败, 用 CUE 的 J 统计量
        if "gmm_j_statistic" not in result:
            result["gmm_j_statistic"] = float(res_cue.j_stat.stat)
            result["gmm_j_pvalue"] = float(res_cue.j_stat.pval)
    except Exception as e:
        print(f"# WARNING: linearmodels IVGMMCUE failed (skipping CUE): {e}",
              file=sys.stderr)
        result["gmm_cue_error"] = str(e)
        result["gmm_cue_skipped"] = True

    return result


# ---------------------------------------------------------------------------
# statsmodels.sandbox.regression.gmm 路径
# ---------------------------------------------------------------------------

def compute_gmm_statsmodels(data: Dict) -> Dict:
    """使用 statsmodels.sandbox.regression.gmm.IVGMM。"""
    y = data["y"]
    x = data["x"]
    z = data["z"]
    z2 = data["z2"]
    n = data["n"]

    # 内生回归矩阵 X = [1, x], 工具矩阵 Z = [1, z, z2]
    X = np.column_stack([np.ones(n), x])
    Z = np.column_stack([np.ones(n), z, z2])

    result = {}
    try:
        # statsmodels IVGMM 需要外生 + 内生 + 工具
        # 接口: IVGMM(endog, exog, instrument) 其中 exog 包含外生 + 内生
        mod = SM_IVGMM(y, X, Z)
        res = mod.fit(maxiter=2, method="bfgs", disp=0)
        result["gmm_twostep_coefficients"] = res.params.tolist()
        result["gmm_twostep_method"] = "statsmodels.sandbox.regression.gmm.IVGMM"
        # J 检验
        try:
            j_stat = float(res.j_test().stat)  # 接口可能因版本不同
            j_pvalue = float(res.j_test().pval)
            result["gmm_j_statistic"] = j_stat
            result["gmm_j_pvalue"] = j_pvalue
        except Exception as e:
            print(f"# WARNING: statsmodels J-test failed: {e}",
                  file=sys.stderr)
    except Exception as e:
        print(f"# WARNING: statsmodels IVGMM failed: {e}", file=sys.stderr)
        result["gmm_twostep_error"] = str(e)

    # statsmodels 无原生 CUE, 标记跳过
    result["gmm_cue_skipped"] = True
    result["gmm_cue_method"] = "skipped (statsmodels.gmm has no CUE)"
    return result


# ---------------------------------------------------------------------------
# numpy 手动两步 GMM (兜底)
# ---------------------------------------------------------------------------

def two_step_gmm_numpy(y: np.ndarray, X: np.ndarray, Z: np.ndarray,
                        W0: np.ndarray = None) -> Tuple[np.ndarray,
                                                          np.ndarray,
                                                          float, float]:
    """手动两步 GMM。

    矩条件: E[Z_i (y_i - X_i β)] = 0
    样本矩: g(β) = (1/n) Z^T (y - X β)
    目标: J(β) = n * g(β)^T W g(β)

    Step 1: W0 = (Z^T Z)^{-1} → β_1 = (X^T Z W0 Z^T X)^{-1} X^T Z W0 Z^T y
    Step 2: S = (1/n) Σ Z_i Z_i^T ε_i^2 (White), W1 = S^{-1}
            β_2 = (X^T Z W1 Z^T X)^{-1} X^T Z W1 Z^T y
    J = n * g(β_2)^T W1 g(β_2) ~ χ²_{L-K} (L=工具数, K=参数数)
    """
    n, k = X.shape
    L = Z.shape[1]
    if W0 is None:
        W0 = np.linalg.inv(Z.T @ Z)

    # Step 1
    ZX = Z.T @ X
    Zy = Z.T @ y
    beta1 = np.linalg.solve(ZX.T @ W0 @ ZX, ZX.T @ W0 @ Zy)

    # 计算 Step 1 残差 → 估计 S
    resid1 = y - X @ beta1
    # White 估计量 S = (1/n) Σ z_i z_i^T ε_i^2
    S = (Z * resid1[:, None]).T @ (Z * resid1[:, None]) / n
    W1 = np.linalg.inv(S)

    # Step 2
    beta2 = np.linalg.solve(ZX.T @ W1 @ ZX, ZX.T @ W1 @ Zy)

    # J 统计量
    resid2 = y - X @ beta2
    g = Z.T @ resid2 / n
    J_stat = float(n * g.T @ W1 @ g)
    df = L - k
    J_pvalue = float(1.0 - st.chi2.cdf(J_stat, df))

    # 协方差 (稳健): (X^T Z W1 Z^T X)^{-1} / n
    vcov = np.linalg.inv(ZX.T @ W1 @ ZX) / n
    return beta2, vcov, J_stat, J_pvalue


def cue_numpy(y: np.ndarray, X: np.ndarray, Z: np.ndarray,
              x0: np.ndarray = None) -> Tuple[np.ndarray, float]:
    """手动 CUE (Continuously Updating Estimator)。

    CUE 同时优化 β 和 W(β):
      β_CUE = argmin_β  n * g(β)^T S(β)^{-1} g(β)
    其中 S(β) = (1/n) Σ z_i z_i^T (y_i - x_i β)^2
    """
    from scipy.optimize import minimize
    n, k = X.shape
    L = Z.shape[1]

    def objective(beta):
        resid = y - X @ beta
        g = Z.T @ resid / n
        S = (Z * resid[:, None]).T @ (Z * resid[:, None]) / n
        try:
            Sinv = np.linalg.inv(S)
        except np.linalg.LinAlgError:
            return 1e20
        return float(n * g.T @ Sinv @ g)

    if x0 is None:
        # 用 2SLS 作起点
        W0 = np.linalg.inv(Z.T @ Z)
        ZX = Z.T @ X
        Zy = Z.T @ y
        x0 = np.linalg.solve(ZX.T @ W0 @ ZX, ZX.T @ W0 @ Zy)

    res = minimize(objective, x0, method="Nelder-Mead",
                    options={"xatol": 1e-10, "fatol": 1e-12,
                              "maxiter": 5000, "maxfev": 5000})
    return res.x, float(res.fun)


def compute_gmm_numpy(data: Dict) -> Dict:
    """纯 numpy 手动两步 GMM + CUE。"""
    y = data["y"]
    x = data["x"]
    z = data["z"]
    z2 = data["z2"]
    n = data["n"]

    X = np.column_stack([np.ones(n), x])
    Z = np.column_stack([np.ones(n), z, z2])

    beta2, vcov, J_stat, J_pvalue = two_step_gmm_numpy(y, X, Z)
    result = {
        "gmm_twostep_coefficients": beta2.tolist(),
        "gmm_twostep_vcov": vcov.tolist(),
        "gmm_j_statistic": float(J_stat),
        "gmm_j_pvalue": float(J_pvalue),
        "gmm_twostep_method": "numpy manual 2-step GMM",
    }

    try:
        beta_cue, J_cue = cue_numpy(y, X, Z, x0=beta2)
        result["gmm_cue_coefficients"] = beta_cue.tolist()
        result["gmm_cue_method"] = "numpy manual CUE (Nelder-Mead)"
        # CUE 的 J 统计量 (使用 CUE 自身 β)
        result["gmm_cue_j_statistic"] = float(J_cue)
    except Exception as e:
        print(f"# WARNING: numpy CUE failed: {e}", file=sys.stderr)
        result["gmm_cue_skipped"] = True
        result["gmm_cue_error"] = str(e)

    return result


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------

def generate_benchmarks(seed: int = 42, n: int = 500) -> Dict:
    """生成所有 GMM 基准值。"""
    data = generate_iv_data(seed=seed, n=n)

    benchmarks = {
        "metadata": {
            "generator": "cross_validate_gmm.py",
            "version": "1.0",
            "phase": "Phase 6 v1.5 M3",
            "cpp_reference": ["test_gmm_two_step", "test_gmm_cue_iterated"],
            "libraries": ["linearmodels", "statsmodels", "numpy"],
            "linearmodels_available": HAS_LINEARMODELS,
            "statsmodels_gmm_available": HAS_SM_GMM,
            "dataset": f"synthetic IV (seed={seed}, N={n})",
            "data_structure": {
                "endogenous": "x = 0.5*z + 0.3*error_x + noise",
                "outcome": "y = 1.0 + 2.0*x + error_y",
                "instruments": "Z = [1, z, z2]",
            },
            "true_params": {"const": 1.0, "x": 2.0},
            "description": "GMM + CUE Python 交叉验证",
            "tolerance_note": "1e-6 (GMM 迭代收敛), 1e-8 (J 统计量公式)",
            "references": [
                "Hansen (1982) Large Sample Properties of GMM Estimators",
                "Hansen-Heaton-Yaron (1996) CUE",
            ],
        },
        "data_summary": {
            "n": int(data["n"]),
            "y_mean": float(np.mean(data["y"])),
            "y_std": float(np.std(data["y"])),
            "x_mean": float(np.mean(data["x"])),
            "x_std": float(np.std(data["x"])),
            "instruments_count": 3,  # [1, z, z2]
            "overidentification_df": 3 - 2,  # L - K = 1
        },
    }

    # 优先 linearmodels, 然后 statsmodels.gmm, 最后 numpy 手动
    if HAS_LINEARMODELS:
        try:
            gmm_result = compute_gmm_linearmodels(data)
            benchmarks.update(gmm_result)
            benchmarks["gmm_backend"] = "linearmodels"
        except Exception as e:
            print(f"# WARNING: linearmodels path failed: {e}", file=sys.stderr)
            # 降级到 statsmodels
            if HAS_SM_GMM:
                benchmarks.update(compute_gmm_statsmodels(data))
                benchmarks["gmm_backend"] = "statsmodels.gmm"
            else:
                benchmarks.update(compute_gmm_numpy(data))
                benchmarks["gmm_backend"] = "numpy"
    elif HAS_SM_GMM:
        benchmarks.update(compute_gmm_statsmodels(data))
        benchmarks["gmm_backend"] = "statsmodels.gmm"
    else:
        benchmarks.update(compute_gmm_numpy(data))
        benchmarks["gmm_backend"] = "numpy"

    benchmarks["tolerances"] = {
        "gmm_twostep_coefficients": 1e-6,
        "gmm_cue_coefficients": 1e-5,        # CUE 优化器更敏感
        "gmm_j_statistic": 1e-6,
        "gmm_j_pvalue": 1e-6,
    }
    return benchmarks


# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Cpp_Hub GMM + CUE 交叉验证基准生成")
    parser.add_argument("--seed", type=int, default=42,
                        help="随机种子 (default: 42)")
    parser.add_argument("--n-samples", type=int, default=500,
                        help="样本数 (default: 500)")
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
