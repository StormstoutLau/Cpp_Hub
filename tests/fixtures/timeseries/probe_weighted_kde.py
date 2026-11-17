"""P2-1/P2-5 参考探针: 反射边界修正与加权 KDE 语义.

用途
----
冻结 statsmodels / scipy 的权重语义差异, 并为 C++ `KernelDensity1D`
的权重路径 (P2-5) 与 0 处反射边界修正 (P2-1) 提供可复现锚点.

关键结论
--------
- statsmodels `KDEUnivariate.fit(fft=False, weights=w)`:
    density(x) = Σ w_i K((x-x_i)/h) / (Σw · h)  (权重仅比例有意义)
  带宽规则 (`bw='normal_reference'` 等) **完全忽略权重** (用未加权样本).
- sklearn `KernelDensity(sample_weight=w)`: 字符串带宽同样忽略权重.
- scipy `gaussian_kde(weights=w)`: 另一套语义 — 加权协方差 + neff^(-1/5);
  本库不采用, 仅登记为同名不同物.
- 反射法 (正支撑 [0,∞)): f_R(x) = Σ p_i [φ((x-x_i)/h)+φ((x+x_i)/h)]/h,
  ∫_0^∞ f_R = 1 且 f_R(0) = 2 f_plain(0) (精确恒等式).

运行
----
    python tests/fixtures/timeseries/probe_weighted_kde.py
退出码: W1-W5 全 PASS 为 0, 任一 FAIL 为 1.
"""

import sys

import numpy as np
from scipy.integrate import quad
from scipy.stats import gaussian_kde, norm
from statsmodels.nonparametric.kde import KDEUnivariate

# 5 点夹具 (加权; 与 test_kde.cpp 的 kSample 一致)
SAMPLE5 = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
WEIGHTS5 = np.array([1.0, 2.0, 3.0, 2.0, 1.0])
GRID5 = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
H5 = 0.5

# 正支撑夹具 (反射; 与 test_kde.cpp 一致)
POS_SAMPLE = np.array([0.2, 0.5, 1.0, 2.0])
POS_WEIGHTS = np.array([1.0, 2.0, 3.0, 4.0])
POS_GRID = np.array([0.0, 0.1, 0.5, 1.0, 2.0, 5.0])
H_POS = 0.4


def direct_weighted(x, w, h, grid):
    """Σ (w_i/Σw) φ((g-x_i)/h) / h (独立参考实现)."""
    x = np.asarray(x, dtype=float)
    wn = np.asarray(w, dtype=float)
    wn = wn / wn.sum()
    z = (np.asarray(grid, dtype=float)[:, None] - x) / h
    return norm.pdf(z).dot(wn) / h


def direct_reflection(x, w, h, grid):
    """Σ (w_i/Σw) [φ((g-x_i)/h) + φ((g+x_i)/h)] / h."""
    x = np.asarray(x, dtype=float)
    wn = np.asarray(w, dtype=float)
    wn = wn / wn.sum()
    g = np.asarray(grid, dtype=float)[:, None]
    return (norm.pdf((g - x) / h) + norm.pdf((g + x) / h)).dot(wn) / h


def statsmodels_density(x, w, h, grid):
    """statsmodels 自身 kernel 的加权密度 (fft=False 路径)."""
    k = KDEUnivariate(x)
    k.fit(bw=h, fft=False, weights=w)
    wn = np.asarray(w, dtype=float) / np.sum(w)
    return k.kernel((np.asarray(grid)[:, None] - x) / h).dot(wn) / h


def check(name, ok, detail):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    return ok


def main():
    ok = True

    # W1: statsmodels 加权密度 == 直接公式; 输出 C++ 锚点
    sm = statsmodels_density(SAMPLE5, WEIGHTS5, H5, GRID5)
    direct = direct_weighted(SAMPLE5, WEIGHTS5, H5, GRID5)
    max_diff = float(np.max(np.abs(sm - direct)))
    ok &= check(
        "W1 statsmodels 加权密度 == Σwφ/(Σw·h)",
        max_diff < 1e-15,
        f"max|diff|={max_diff:.3e}; anchors={[float(v) for v in sm]}",
    )

    # W2: statsmodels 带宽规则忽略权重
    k_w = KDEUnivariate(SAMPLE5)
    k_w.fit(fft=False, weights=WEIGHTS5)
    k_0 = KDEUnivariate(SAMPLE5)
    k_0.fit(fft=False)
    ok &= check(
        "W2 statsmodels 带宽忽略权重",
        float(k_w.bw) == float(k_0.bw),
        f"weighted={float(k_w.bw):.17g}; unweighted={float(k_0.bw):.17g}",
    )

    # W3: 权重仅比例有意义 (×7 不变)
    sm7 = statsmodels_density(SAMPLE5, 7.0 * WEIGHTS5, H5, GRID5)
    d7 = float(np.max(np.abs(sm7 - sm)))
    ok &= check(
        "W3 权重尺度不变",
        d7 < 1e-15,
        f"max|diff(7w, w)|={d7:.3e}",
    )

    # W4: scipy 加权语义 (登记为不同实现)
    sc = gaussian_kde(SAMPLE5, weights=WEIGHTS5)
    neff = float(WEIGHTS5.sum() ** 2 / np.sum(WEIGHTS5**2))
    ok &= check(
        "W4 scipy 加权: factor=neff^(-1/5), 加权协方差",
        abs(float(sc.neff) - neff) < 1e-12
        and abs(float(sc.factor) - neff ** -0.2) < 1e-12,
        f"neff={float(sc.neff):.15g}; factor={float(sc.factor):.15g}; "
        f"weighted sd={float(np.sqrt(sc.covariance[0, 0])):.15g}",
    )

    # W5: 反射恒等式 (质量 1; f_R(0)=2 f_plain(0)); 输出 C++ 锚点
    mass = quad(lambda t: float(direct_reflection(POS_SAMPLE, np.ones(4), H_POS, [t])[0]),
                0.0, np.inf, limit=400)[0]
    f_r0 = float(direct_reflection(POS_SAMPLE, np.ones(4), H_POS, [0.0])[0])
    f_p0 = float(direct_weighted(POS_SAMPLE, np.ones(4), H_POS, [0.0])[0])
    refl = direct_reflection(POS_SAMPLE, np.ones(4), H_POS, POS_GRID)
    wrefl = direct_reflection(POS_SAMPLE, POS_WEIGHTS, H_POS, POS_GRID)
    ok &= check(
        "W5 反射: ∫_0^∞=1 且 f_R(0)=2f(0)",
        abs(mass - 1.0) < 1e-12 and abs(f_r0 - 2.0 * f_p0) < 1e-15,
        f"mass={mass:.15g}; f_R(0)={f_r0:.17g}; 2f(0)={2.0 * f_p0:.17g}; "
        f"anchors={[float(v) for v in refl]}; "
        f"weighted anchors={[float(v) for v in wrefl]}",
    )

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
