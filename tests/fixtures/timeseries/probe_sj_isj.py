"""P2-2 参考探针: SJ (Sheather-Jones 1991) 与 ISJ (Botev et al. 2010) 带宽.

- SJ: 逐步翻译 R `stats::bw.SJ(x, method="ste")` 的实现
  (R 源码 src/library/stats/R/bandwidths.R + src/bandwidths.c:
   bw_pair_cnts 分箱, bw_phi4/bw_phi6, scale=min(sd,IQR/1.349),
   a=1.24*scale*n^(-1/7), b=1.23*scale*n^(-1/9), alph2=1.357*(phi4(a)/TD)^(1/7),
   uniroot 求 fSD(h)=(c1/phi4(alph2*h^(5/7)))^(1/5)-h).
  本探针用高精度二分求根, 与 R 默认宽松容差 (tol=0.1*lower) 的差异在 R 自身
  容差内; 若安装 R, 可用 bw.SJ(x, tol=1e-14) 做逐位对照.

- ISJ: 逐步翻译 KDEpy `KDEpy.bw_selection.improved_sheather_jones`
  (n=2^10 网格, autogrid boundary_abs=6/boundary_rel=0.5, 线性分箱,
   DCT-II, _fixed_point + brentq), 返回 sqrt(t*)*R.

数据为确定性整数 LCG 夹具 (跨语言可复现, 无 libm 依赖).
"""

import numpy as np
from scipy import fftpack
from scipy.optimize import brentq


def grid_series(n: int = 150) -> np.ndarray:
    """LCG 网格序列: x_i = ((i*7919 + 12345) % 1000)/100 (确定性)."""
    return np.array([((i * 7919 + 12345) % 1000) / 100.0 for i in range(n)])


def ar_series(n: int = 200, seed: int = 12345) -> np.ndarray:
    """整数 AR(0.75): x_i = 3*x_{i-1}//4 + u_i, u_i ∈ {0..3} (uint64 LCG)."""
    s = np.uint64(seed)
    x = np.zeros(n)
    for i in range(n):
        s = np.uint64(s) * np.uint64(6364136223846793005) + np.uint64(1442695040888963407)
        u = float((int(s) >> 33) % 4)
        x[i] = u if i == 0 else (3.0 * x[i - 1]) // 4.0 + u
    return x


def lcg_series(n: int = 200, seed: int = 12345) -> np.ndarray:
    """多唯一值连续型夹具: x_i = ((s_i >> 11) % 1000003)/100000 (确定性)."""
    s = np.uint64(seed)
    x = np.zeros(n)
    with np.errstate(over="ignore"):
        for i in range(n):
            s = np.uint64(s) * np.uint64(6364136223846793005) + np.uint64(
                1442695040888963407
            )
            x[i] = ((int(s) >> 11) % 1000003) / 100000.0
    return x


def _pair_counts(x: np.ndarray, nb: int = 1000):
    """R bw_pair_cnts (binned 路径): 返回 (dd, cnt[0..nb])."""
    xmin, xmax = x.min(), x.max()
    dd = (xmax - xmin) * 1.01 / nb
    if dd == 0.0:
        raise ValueError("data are constant")
    bins = np.trunc(np.abs(x) / dd) * np.sign(x)
    bins = bins - bins.min() + 1
    w = np.bincount(bins.astype(np.int64), minlength=nb + 2)
    cnt = np.zeros(nb + 1)
    for ii in range(1, len(w)):
        wi = w[ii]
        if wi == 0:
            continue
        cnt[0] += wi * (wi - 1.0)
        for jj in range(1, ii):
            if w[jj]:
                cnt[ii - jj] += wi * w[jj]
    cnt[0] *= 0.5
    return dd, cnt


def sj_bandwidth(x: np.ndarray, nb: int = 1000) -> float:
    """R bw.SJ(method='ste') 的等价实现 (精确求根)."""
    x = np.asarray(x, dtype=float)
    n = x.size
    dd, cnt = _pair_counts(x, nb)
    k = np.arange(cnt.size)
    delta = (k * dd) ** 2

    def phi_poly(h, power, poly_deg):
        d = delta / (h * h)
        m = d < 1000.0
        d = d[m]
        c = cnt[m]
        if poly_deg == 4:
            term = np.exp(-d / 2.0) * (d * d - 6.0 * d + 3.0)
            total = 2.0 * np.sum(term * c) + 3.0 * n
            return total / (n * (n - 1.0) * h**5) / np.sqrt(2.0 * np.pi)
        else:
            term = np.exp(-d / 2.0) * (
                d * d * d - 15.0 * d * d + 45.0 * d - 15.0
            )
            total = 2.0 * np.sum(term * c) - 15.0 * n
            return total / (n * (n - 1.0) * h**7) / np.sqrt(2.0 * np.pi)

    def phi4(h):
        return phi_poly(h, 5, 4)

    def phi6(h):
        return phi_poly(h, 7, 6)

    sd = np.std(x, ddof=1)
    q25, q75 = np.percentile(x, [25, 75])
    scale = min(sd, (q75 - q25) / 1.349)
    a = 1.24 * scale * n ** (-1.0 / 7.0)
    b = 1.23 * scale * n ** (-1.0 / 9.0)
    c1 = 1.0 / (2.0 * np.sqrt(np.pi) * n)
    td = -phi6(b)
    if not np.isfinite(td) or td <= 0:
        raise ValueError("sample is too sparse to find TD")
    hmax = 1.144 * scale * n ** (-1.0 / 5.0)
    lower, upper = 0.1 * hmax, hmax
    alph2 = 1.357 * (phi4(a) / td) ** (1.0 / 7.0)

    def f(h):
        return (c1 / phi4(alph2 * h ** (5.0 / 7.0))) ** 0.2 - h

    # 区间扩展与 R 一致 (交替 upper*1.2 / lower/1.2)
    itry = 1
    while f(lower) * f(upper) > 0:
        if itry > 99:
            raise ValueError("no solution in the specified range")
        if itry % 2:
            upper *= 1.2
        else:
            lower /= 1.2
        itry += 1
    return brentq(f, lower, upper, xtol=1e-15, rtol=1e-15, maxiter=500)


def isj_bandwidth(x: np.ndarray) -> float:
    """KDEpy improved_sheather_jones 的等价实现."""
    x = np.asarray(x, dtype=float)
    n = 2**10
    rng = x.max() - x.min()
    if rng == 0.0:
        raise ValueError("data are constant")
    border = max(0.5 * rng, 6.0)
    lo, hi = x.min() - border, x.max() + border
    grid = np.linspace(lo, hi, num=n)
    dx = (hi - lo) / (n - 1)

    # 线性分箱 (KDEpy linbin_numpy, 权重 1/N)
    t = (x - lo) / dx
    frac, integ = np.modf(t)
    integ = integ.astype(np.int64)
    hist = np.zeros(n)
    w = np.full(x.size, 1.0 / x.size)
    np.add.at(hist, integ, (1.0 - frac) * w)
    np.add.at(hist, integ + 1, frac * w)

    coef = fftpack.dct(hist)  # DCT-II, norm=None
    i_sq = np.arange(1, n, dtype=float) ** 2
    a2 = coef[1:] ** 2
    # KDEpy 细节: _root 只用 clamped N 计算容差, 而 _fixed_point 通过 args 收到
    # **原始唯一值数** (N 在 _root 内被钳制但 args 已构造) — 重离散数据差异显著
    n_fp = float(len(np.unique(x)))
    nu = min(max(len(np.unique(x)), 50), 1050)

    def fsum(t, ell):
        return 0.5 * np.pi ** (2 * ell) * np.sum(
            i_sq**ell * a2 * np.exp(-i_sq * np.pi**2 * t)
        )

    def fixed_point(t):
        f = fsum(t, 7)
        if f <= 0:
            return -1.0
        for s in range(6, 1, -1):  # KDEpy: reversed(range(2, ell)), ell=7
            odd = np.prod(np.arange(1, 2 * s + 1, 2, dtype=float))
            k0 = odd / np.sqrt(2.0 * np.pi)
            const = (1.0 + 0.5 ** (s + 0.5)) / 3.0
            time = (2.0 * const * k0 / (n_fp * f)) ** (2.0 / (3.0 + 2.0 * s))
            f = fsum(time, s)
        t_opt = (2.0 * n_fp * np.sqrt(np.pi) * f) ** (-2.0 / 5.0)
        return t - t_opt

    tol = 10e-12 + 0.01 * (nu - 50) / 1000
    while fixed_point(tol) <= 0:
        tol *= 2.0
        if tol >= 1.0:
            raise ValueError("root finding did not converge")
    t_star = brentq(fixed_point, 0.0, tol, xtol=1e-15, rtol=1e-15, maxiter=500)
    return float(np.sqrt(t_star) * rng)


if __name__ == "__main__":
    # KDEpy _fixed_point doctest 自检 (保证翻译正确)
    def _fp_check(t, i_sq, a2, nu):
        def fsum(t, ell):
            return 0.5 * np.pi ** (2 * ell) * np.sum(
                i_sq**ell * a2 * np.exp(-i_sq * np.pi**2 * t)
            )

        f = fsum(t, 7)
        for s in range(6, 1, -1):
            odd = np.prod(np.arange(1, 2 * s + 1, 2, dtype=float))
            k0 = odd / np.sqrt(2.0 * np.pi)
            const = (1.0 + 0.5 ** (s + 0.5)) / 3.0
            time = (2.0 * const * k0 / (nu * f)) ** (2.0 / (3.0 + 2.0 * s))
            f = fsum(time, s)
        t_opt = (2.0 * nu * np.sqrt(np.pi) * f) ** (-2.0 / 5.0)
        return t - t_opt

    print("doctest1:", _fp_check(0.01, np.arange(1, 51), np.arange(1, 51), 50))
    print("doctest2:", _fp_check(0.07, np.arange(1, 11), np.arange(1, 11), 25))

    for name, data in [
        ("grid150", grid_series()),
        ("ar200", ar_series()),
        ("lcg200", lcg_series()),
    ]:
        print(
            f"{name}: sj={sj_bandwidth(data):.15f} isj={isj_bandwidth(data):.15f} "
            f"n={data.size} uniq={len(np.unique(data))} "
            f"sd={np.std(data, ddof=1):.15f} range={data.max()-data.min():.15f}"
        )
