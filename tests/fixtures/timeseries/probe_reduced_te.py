"""H4 探针: Kirkley (2025) reduced TE (Eq. 12/13) 在离散化 Gaussian VAR 上的偏差

背景 (spec §10.1 H4): "Kirkley reduced TE 是否能把 +25% 偏差压到 ±10%"。
关键前提: reduced TE 定义在**有限基数符号序列**上 (Eq. 8 的组合计数);
连续 KSG-CMI 的 +25% 偏置属于另一估计器现象, 不能由 reduced TE 直接消除。
本脚本:
  1. 实现 Eq. 10 (组合条件熵 H_C) 与 Eq. 12/13 (reduced TE), 与 plug-in 离散 TE
     (Eq. 3) 对照;
  2. 在 Gaussian VAR(φ=0.5, a=0.3) 分位分箱 (C=4) 上扫描 T, 以 T_ref=200000 的
     plug-in 值为参照, 比较两种估计的有限样本偏差;
  3. 零假设 (两条独立 AR(1)) 下 reduced TE 的负值比例 (MDL 显著性语义)。

用法: python tests/fixtures/timeseries/probe_reduced_te.py
"""
import numpy as np

PHI = 0.5
A = 0.3
C = 4          # 分箱数 (符号基数)
K = 1          # x 历史阶
L = 1          # y 历史阶
TS = [200, 500, 2000, 20000]
T_REF = 200000


def gen_var(n: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    x = np.empty(n)
    y = np.empty(n)
    x[0] = rng.normal() / np.sqrt(1.0 - PHI * PHI)
    y[0] = rng.normal() / np.sqrt(1.0 - PHI * PHI)
    ez = rng.normal(size=(n - 1, 2))
    for t in range(1, n):
        x[t] = PHI * x[t - 1] + ez[t - 1, 0]
        y[t] = PHI * y[t - 1] + A * x[t - 1] + ez[t - 1, 1]
    return x, y


def quantize(x: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """按参照样本的分位点分箱 (out-of-sample 口径, 避免小样本分位漂移)"""
    edges = np.quantile(ref, np.linspace(0.0, 1.0, C + 1)[1:-1])
    return np.clip(np.searchsorted(edges, x, side="right"), 0, C - 1)


def _log_choose(a: np.ndarray, b: int) -> np.ndarray:
    """log C(a + b - 1, b - 1), 组合数用 lgamma 稳定计算"""
    from scipy.special import gammaln
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    return gammaln(a + b) - gammaln(a + 1.0) - gammaln(b)


def _log_factorial(counts: np.ndarray) -> float:
    from scipy.special import gammaln
    counts = np.asarray(counts, dtype=np.float64)
    counts = counts[counts > 0]
    return float(np.sum(gammaln(counts + 1.0)))


def _entropy_from_counts(counts: np.ndarray) -> float:
    p = counts.astype(np.float64)
    p = p[p > 0]
    p /= p.sum()
    return float(-np.sum(p * np.log2(p)))


def tables(yf: np.ndarray, yp: np.ndarray, xp: np.ndarray) -> dict:
    n_qr = np.bincount(yf * C + yp, minlength=C * C).astype(np.float64)
    n_rs = np.bincount(yp * C + xp, minlength=C * C).astype(np.float64)
    n_qrs = np.bincount((yf * C + yp) * C + xp,
                        minlength=C * C * C).astype(np.float64)
    n_r = np.bincount(yp, minlength=C).astype(np.float64)
    n_q = np.bincount(yf, minlength=C).astype(np.float64)
    return {"n_qr": n_qr, "n_rs": n_rs, "n_qrs": n_qrs, "n_r": n_r, "n_q": n_q}


def plugin_te(t: dict, n: int) -> float:
    """Eq. 3/5: H_S(y+|Ypast) - H_S(y+|Ypast,Xpast) (bits)"""
    # H(y+|Ypast) = H(y+,Ypast) - H(Ypast); counts 为联合表
    h_qr = _entropy_from_counts(t["n_qr"])
    h_r = _entropy_from_counts(t["n_r"])
    h_qrs = _entropy_from_counts(t["n_qrs"])
    h_rs = _entropy_from_counts(t["n_rs"])
    return (h_qr - h_r) - (h_qrs - h_rs)


def reduced_te(t: dict, n: int) -> float:
    """Eq. 12/13: H_C(y+|Ypast) - H_C(y+|Ypast,Xpast) (bits)"""
    log2 = np.log(2.0)
    delta = (_log_choose(t["n_r"], C).sum() -
             _log_choose(t["n_rs"], C).sum()) / n / log2
    ratio = (_log_factorial(t["n_qrs"]) + _log_factorial(t["n_r"]) -
             _log_factorial(t["n_qr"]) - _log_factorial(t["n_rs"])) / n / log2
    return delta + ratio


def estimate(yf: np.ndarray, yp: np.ndarray, xp: np.ndarray) -> tuple[float, float]:
    n = len(yf)
    t = tables(yf, yp, xp)
    return plugin_te(t, n), reduced_te(t, n)


def main() -> None:
    # 参照: 大样本 plug-in (偏差 O(1/T_ref))
    xr, yr = gen_var(T_REF, 20261117)
    xq = quantize(xr, xr)
    yq = quantize(yr, yr)
    yf_r, yp_r, xp_r = yq[L:], yq[L - 1:-1], xq[K:-1]
    # y_t 与 x_t 对齐: 预测 y_{t+1}, 历史 y_t/x_t
    yf_r = yq[1:]
    yp_r = yq[:-1]
    xp_r = xq[:-1]
    ref_plugin, ref_reduced = estimate(yf_r, yp_r, xp_r)
    print(f"reference T={T_REF}: plugin={ref_plugin:.6f} bits "
          f"reduced={ref_reduced:.6f} bits")

    print(f"\nC={C} K={K} L={L}; 偏差相对 T_ref 的参照值")
    print(f"{'T':>7} {'plugin':>10} {'bias%':>8} {'reduced':>10} {'bias%':>8}")
    for T in TS:
        x, y = gen_var(T, 20261117)
        # 分位点用参照样本 (out-of-sample 口径)
        xq = quantize(x, xr)
        yq = quantize(y, yr)
        yf, yp, xp = yq[1:], yq[:-1], xq[:-1]
        pl, rd = estimate(yf, yp, xp)
        b_pl = 100.0 * (pl - ref_plugin) / ref_plugin
        b_rd = 100.0 * (rd - ref_reduced) / ref_reduced
        print(f"{T:7d} {pl:10.6f} {b_pl:8.1f} {rd:10.6f} {b_rd:8.1f}")

    # 零假设: 独立 AR(1), 统计 reduced TE 的符号 (MDL 显著性语义)
    neg = 0
    runs = 40
    for s in range(runs):
        rng = np.random.default_rng(7000 + s)
        n = 500
        xs = np.empty(n)
        ys = np.empty(n)
        xs[0] = rng.normal()
        ys[0] = rng.normal()
        e = rng.normal(size=(n - 1, 2))
        for t in range(1, n):
            xs[t] = 0.8 * xs[t - 1] + e[t - 1, 0]
            ys[t] = 0.8 * ys[t - 1] + e[t - 1, 1]
        xs_q = np.digitize(xs, np.quantile(xs, np.linspace(0, 1, C + 1)[1:-1]))
        ys_q = np.digitize(ys, np.quantile(ys, np.linspace(0, 1, C + 1)[1:-1]))
        _, rd = estimate(ys_q[1:], ys_q[:-1], xs_q[:-1])
        if rd < 0:
            neg += 1
    print(f"\nnull (independent AR(1), n=500, C={C}): "
          f"reduced TE < 0 in {neg}/{runs} runs")


if __name__ == "__main__":
    main()
