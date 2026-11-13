# =============================================================================
# gen_phase9_baselines.py - Phase 9 C++ 测试黄金值生成器 (spec §1.1 交付物)
#
# 基准来源:
#   - scipy.stats.gaussian_kde     同核同 h pdf (probe P6 实测互差 7.494e-16)
#   - sklearn.feature_selection._mutual_info._compute_mi_cc (KSG Alg1 私有
#     API, 仅生成期使用, 审计 C4: 不进运行时测试)
#   - probe_sklearn_kde.py / probe_ksg.py (tests/fixtures/timeseries/) 的
#     语义冻结输出 (带宽公式常数 C≈1.05922 / GC=ln(1+a²/(1-φ²)) / TE=GC/2)
#
# 用途: 复跑输出 test_kde.cpp / test_knn_density.cpp / test_ksg.cpp 内嵌
#   黄金值; 与测试文件内常量逐位对照即为验收 (G3 容差分层 §1.3)。
#   测试文件不 include 本脚本输出 (.inc), 黄金值以内嵌常量 + 头注溯源
#   交付 (Phase 7C/8 同款模式, 避免 CI 依赖 Python)。
#
# 复现: python gen_phase9_baselines.py  (numpy/scipy/sklearn 装机环境)
# 版本: numpy 2.x / scipy 1.x / scikit-learn 1.x (probe 报告 §1 冻结版本)
# =============================================================================

import numpy as np
from scipy.stats import gaussian_kde
from sklearn.feature_selection._mutual_info import _compute_mi_cc


def kde_section() -> None:
    """test_kde.cpp: 5 点样本 {-1,-0.5,0,0.5,1} 的带宽常数 + pdf/log_pdf 黄金值."""
    x = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    n = x.size
    sd = float(np.std(x, ddof=1))
    q25, q75 = np.percentile(x, [25, 75])
    a = min(sd, (q75 - q25) / 1.349)
    n_pow = float(n ** -0.2)
    print("[test_kde] sample = {-1, -0.5, 0, 0.5, 1}")
    print(f"  kSd        = {sd!r}          # ddof=1")
    print(f"  kIqr       = {float(q75 - q25)!r}                      # q75-q25")
    print(f"  kA         = {float(a)!r}   # min(sd, IQR/1.349)")
    print(f"  kScott     = {float(1.059 * a * n_pow)!r}  # 1.059*A*n^-0.2")
    print(f"  kSilverman = {float(0.9 * a * n_pow)!r}  # 0.9*A*n^-0.2")
    # statsmodels normal_reference 常数 C: probe P4 反解冻结 (audit A2)
    c_nr = 1.059223840025653
    print(f"  kNormalRef = {float(c_nr * a * n_pow)!r}  # C*A*n^-0.2, C={c_nr}")

    h = 0.5
    grid = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    kde = gaussian_kde(x, bw_method=h / sd)  # factor 语义: h = value*sd
    pdf = kde.pdf(grid)
    print(f"  h={h}; kGrid = {{-1, -0.5, 0, 0.5, 1}}")
    print("  kPdf    = {" + ", ".join(f"{v:.17g}" for v in pdf) + "}")
    print("  kLogPdf = {" + ", ".join(f"{v:.17g}" for v in np.log(pdf)) + "}")


def knn_section() -> None:
    """test_knn_density.cpp: 样本 {0..4} 的 Loftsgaarden-Quesenberry 密度."""
    s = np.array([0.0, 1.0, 2.0, 3.0, 4.0])

    def knn_pdf(x0: float, k: int) -> float:
        d = np.sort(np.abs(s - x0))
        r = d[k - 1]
        return k / (len(s) * 2 * r)

    print("[test_knn_density] sample = {0,1,2,3,4}")
    for x0 in (1.5, 0.0, 2.0, 10.0):
        for k in (2, 4):
            print(f"  pdf({x0:4.1f}, k={k}) = {float(knn_pdf(x0, k))!r}")


def ksg_section() -> None:
    """test_ksg.cpp: 40 点二元 Gaussian (rho=0.7, seed 12345) sklearn 黄金值."""
    rng = np.random.default_rng(12345)
    rho = 0.7
    z = rng.multivariate_normal([0, 0], [[1, rho], [rho, 1]], size=40)
    xs, ys = z[:, 0], z[:, 1]
    print("[test_ksg] 40-pt bivariate Gaussian rho=0.7, numpy seed 12345")
    print("  kX = {" + ", ".join(f"{v:.17g}" for v in xs) + "}")
    print("  kY = {" + ", ".join(f"{v:.17g}" for v in ys) + "}")
    print(f"  kMiGoldenK3 = {float(_compute_mi_cc(xs, ys, 3))!r}  # raw, no jitter")
    print(f"  kMiGoldenK5 = {float(_compute_mi_cc(xs, ys, 5))!r}  # raw, no jitter")
    print(f"  analytic(rho=0.8, integration anchor) = {float(-0.5 * np.log(1 - 0.8 ** 2))!r}")


def te_section() -> None:
    """test_transfer_entropy.cpp / test_integration_phase9.cpp 解析锚."""
    phi, a = 0.5, 0.3
    gc = float(np.log(1.0 + a * a / (1.0 - phi * phi)))
    print("[test_te] Gaussian VAR(phi=0.5, a=0.3)")
    print(f"  kGc = {gc!r}  # ln(1+a^2/(1-phi^2)) (Barnett 2009)")
    print(f"  kTe = {gc / 2.0!r}  # GC/2 (nats)")


if __name__ == "__main__":
    kde_section()
    knn_section()
    ksg_section()
    te_section()
