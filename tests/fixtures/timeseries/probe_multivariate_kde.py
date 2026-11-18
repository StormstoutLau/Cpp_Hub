"""P2-3 参考探针: 多变量 Gaussian KDE 语义与数值锚.

用途
----
冻结 scipy `gaussian_kde` 的全协方差带宽语义, 为 C++ `KernelDensityND`
(P2-3) 提供可复现锚点:

    factor  = n^(-1/(d+4))                 ('scott', 默认)
            = (n(d+2)/4)^(-1/(d+4))        ('silverman')
            = 用户给定                       (bw_method=float)
    H       = factor^2 * S                 (S = 样本协方差, ddof=1)
    pdf(x)  = (2π)^(-d/2) |H|^(-1/2) / n Σ_i exp(-½ (x-x_i)^T H^-1 (x-x_i))

数据为确定性 uint64 LCG 夹具 (跨语言可复现), 无 RNG 版本依赖.

运行
----
    python tests/fixtures/timeseries/probe_multivariate_kde.py
退出码: W1-W3 全 PASS 为 0, 任一 FAIL 为 1.
"""

import sys

import numpy as np
from scipy.stats import gaussian_kde

N2 = 120
N3 = 150


def _lcg(n, seed=12345):
    s = np.uint64(seed)
    out = np.empty(n)
    with np.errstate(over="ignore"):
        for i in range(n):
            s = np.uint64(s) * np.uint64(6364136223846793005) + np.uint64(
                1442695040888963407
            )
            out[i] = (int(s) >> 11) / 2.0**53
    return out


def points_2d(n=N2):
    u = _lcg(n, 12345)
    v = _lcg(n, 54321)
    x = 10.0 * u
    y = 0.6 * x + 2.0 * v
    return np.column_stack([x, y])


def points_3d(n=N3):
    u = _lcg(n, 12345)
    v = _lcg(n, 54321)
    w = _lcg(n, 99999)
    x = 10.0 * u
    y = 0.6 * x + 2.0 * v
    z = -0.3 * x + 0.4 * y + 1.5 * w
    return np.column_stack([x, y, z])


GRID2 = np.array([[0.0, 0.0], [1.0, 1.0], [2.5, 0.5], [-1.0, 2.0], [8.0, 6.0]])
GRID3 = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.5, 0.5, 1.0],
                  [-1.0, 2.0, 0.0], [8.0, 6.0, 4.0]])


def direct_density(points, cov, factor, grid):
    """Cholesky 直接参考实现 (与 scipy 参数化一致)."""
    h = factor**2 * cov
    chol = np.linalg.cholesky(h)
    d = points.shape[1]
    log_norm = -0.5 * d * np.log(2.0 * np.pi) - np.sum(np.log(np.diag(chol))) \
        - np.log(points.shape[0])
    out = np.empty(len(grid))
    for i, x in enumerate(grid):
        z = np.linalg.solve(chol, x[:, None] - points.T)
        quad = np.sum(z * z, axis=0)
        out[i] = np.exp(log_norm) * np.sum(np.exp(-0.5 * quad))
    return out


def check(name, ok, detail):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    return ok


def main():
    ok = True
    for label, pts, grid in (("2d", points_2d(), GRID2), ("3d", points_3d(), GRID3)):
        d = pts.shape[1]
        cov = np.cov(pts, rowvar=False, ddof=1)
        for method in ("scott", "silverman"):
            kde = gaussian_kde(pts.T, bw_method=method)
            ref = direct_density(pts, cov, float(kde.factor), grid)
            got = np.asarray(kde.pdf(grid.T))
            diff = float(np.max(np.abs(got - ref)))
            ok &= check(
                f"{label}/{method} scipy == Cholesky 参考",
                diff < 1e-14,
                f"factor={float(kde.factor):.17g}; max|diff|={diff:.3e}",
            )
            print(f"    anchors({label},{method}) = "
                  f"{[float(v) for v in got]}")
            print(f"    logpdf({label},{method}) = "
                  f"{[float(v) for v in kde.logpdf(grid.T)]}")
        # 用户因子语义: H = factor^2 * S
        kde_u = gaussian_kde(pts.T, bw_method=0.37)
        ref_u = direct_density(pts, cov, 0.37, grid)
        got_u = np.asarray(kde_u.pdf(grid.T))
        ok &= check(
            f"{label}/user factor H=f^2 S",
            float(np.max(np.abs(got_u - ref_u))) < 1e-14,
            f"max|diff|={float(np.max(np.abs(got_u - ref_u))):.3e}; "
            f"anchors={[float(v) for v in got_u]}",
        )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
