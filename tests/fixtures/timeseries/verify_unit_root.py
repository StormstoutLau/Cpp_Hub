# verify_unit_root.py - M2 单位根/方差比 arch 8.0.0 对照基准生成
# 输出: tests/unit/timeseries/unit_root_baseline.inc (HFE 硬编码策略)
# 用法: python tests/fixtures/timeseries/verify_unit_root.py
import numpy as np
from arch.unitroot import ADF, PhillipsPerron, KPSS, DFGLS, VarianceRatio
from arch.unitroot.unitroot import mackinnonp, mackinnoncrit, kpss_crit
from arch.utility.cov import cov_nw

OUT = r"f:\Cpp_Hub\tests\unit\timeseries\unit_root_baseline.inc"

# ---------------------------------------------------------------- 模拟数据
rng = np.random.default_rng(20260815)
T = 250
e = rng.standard_normal(T + 200)
y_rw = np.cumsum(e[100:100 + T])            # 随机游走 (单位根)
y_ar = np.empty(T)                          # 平稳 AR(1) phi=0.5
y_ar[0] = e[0] / np.sqrt(1 - 0.25)
for t in range(1, T):
    y_ar[t] = 0.5 * y_ar[t - 1] + e[t]
y_tr = 0.05 * np.arange(1, T + 1) + y_ar    # 趋势 + AR 噪声 (自动选 ct)

TP = 501
rp = rng.standard_normal(TP + 200)
r_rw = 0.0005 + 0.01 * rp[100:100 + TP]     # 随机游走价格收益率
p_rw = 100 * np.exp(np.cumsum(r_rw))
r_ar = np.empty(TP)                         # AR(1) 收益率 (VR 应拒绝)
r_ar[0] = 0.3 * 0.0 + rp[0] / np.sqrt(1 - 0.09)
for t in range(1, TP):
    r_ar[t] = 0.3 * r_ar[t - 1] + 0.01 * rp[t]
p_ar = 100 * np.exp(np.cumsum(0.0005 + r_ar))

tiny_y = np.arange(1.0, 11.0)               # [1..10] 手算小样本
# VR 小样本: 上升趋势 + 交替噪声 (完美线性会使方差=0 → VR nan)
tiny_p = np.array([1.0, 1.3, 1.1, 1.6, 1.2, 1.9, 1.5, 2.0, 1.7, 2.2, 1.8])

SCHWERT_250 = int(np.ceil(12.0 * (T / 100.0) ** 0.25))  # = 16

L = []  # 输出行


def w(line=""):
    L.append(line)


def fmt(x):
    """numpy 标量 → C++ double 字面量 (numpy 2.x repr 会输出 np.float64(...))"""
    return repr(float(x))


def tag(x):
    """数值 → C++ 标识符片段 (3.5 → '3p5')"""
    return f"{abs(x):g}".replace(".", "p")


def vec(name, arr, per=4):
    arr = np.asarray(arr, dtype=float)
    w(f"constexpr Real {name}[{len(arr)}] = {{")
    for i in range(0, len(arr), per):
        w("    " + ", ".join(repr(float(x)) for x in arr[i:i + per]) + ",")
    w("};")
    w()


w("// 自动生成: tests/fixtures/timeseries/verify_unit_root.py (arch 8.0.0)")
w("// 勿手改 — 重新生成请运行该脚本")
w("#pragma once")
w("namespace cpphub { inline namespace v1 { namespace timeseries {")
w("namespace unit_root { namespace baseline {")
w()
w(f"constexpr Size T = {T};")
w(f"constexpr Size TP = {TP};")
w(f"constexpr Size TINY_N = 10;")
w(f"constexpr Size TINY_P_N = 11;")
w(f"constexpr Real SCHWERT_250 = {SCHWERT_250};")  # int
w()
vec("Y_RW", y_rw)
vec("Y_AR", y_ar)
vec("Y_TR", y_tr)
vec("P_RW", p_rw)
vec("P_AR", p_ar)
vec("TINY_Y", tiny_y)
vec("TINY_P", tiny_p)

# ---------------------------------------------------------------- ADF (固定 lag=16)
w("// --- ADF 固定 lag=16 (Schwert(250)) ---")
w("struct AdfCase { Real stat, p, cv1, cv5, cv10; Size nobs, n_lags; };")
w(f"constexpr AdfCase ADF_FIXED[6] = {{")
for series, sname in [(y_rw, "rw"), (y_ar, "ar")]:
    for tr in ["n", "c", "ct"]:
        r = ADF(series, lags=SCHWERT_250, trend=tr)
        cv = r.critical_values
        w(f"    {{{fmt(r.stat)}, {fmt(r.pvalue)}, {fmt(cv['1%'])}, "
          f"{fmt(cv['5%'])}, {fmt(cv['10%'])}, {int(r.nobs)}, {int(r.lags)}}},"
          f"  // {sname}/{tr}")
w("};")
w()

# ---------------------------------------------------------------- ADF AIC 自动 lag
w("// --- ADF AIC 自动 lag (arch lags=None, method='aic') ---")
w("constexpr AdfCase ADF_AIC[4] = {")
for series, sname in [(y_rw, "rw"), (y_ar, "ar")]:
    for tr in ["n", "c"]:
        r = ADF(series, trend=tr)
        cv = r.critical_values
        w(f"    {{{fmt(r.stat)}, {fmt(r.pvalue)}, {fmt(cv['1%'])}, "
          f"{fmt(cv['5%'])}, {fmt(cv['10%'])}, {int(r.nobs)}, {int(r.lags)}}},"
          f"  // {sname}/{tr}")
w("};")
w()

# ---------------------------------------------------------------- PP
w("// --- PP (default bandwidth=Schwert) + 显式 lags=8 + tiny 手算 ---")
w("struct PpCase { Real stat, p, cv1, cv5, cv10; Size nobs, lags; };")
w("constexpr PpCase PP_CASES[5] = {")
for series, sname in [(y_rw, "rw"), (y_ar, "ar")]:
    for tr in ["c", "ct"]:
        r = PhillipsPerron(series, trend=tr)
        cv = r.critical_values
        w(f"    {{{fmt(r.stat)}, {fmt(r.pvalue)}, {fmt(cv['1%'])}, "
          f"{fmt(cv['5%'])}, {fmt(cv['10%'])}, {int(r.nobs)}, {int(r.lags)}}},"
          f"  // {sname}/{tr}")
r = PhillipsPerron(y_rw, lags=8, trend="c")
cv = r.critical_values
w(f"    {{{fmt(r.stat)}, {fmt(r.pvalue)}, {fmt(cv['1%'])}, "
  f"{fmt(cv['5%'])}, {fmt(cv['10%'])}, {int(r.nobs)}, {int(r.lags)}}},  # rw/c lags=8"
  .replace("#", "//"))
w("};")
w()

# PP tiny 手算分量: y=[1..10], trend n, lags=2
from statsmodels.regression.linear_model import OLS
rhs = tiny_y[:-1][:, None]
res = OLS(tiny_y[1:], rhs).fit()
u = res.resid
n, k = res.nobs, 1
lam2 = cov_nw(u, 2, demean=False)
s2 = u @ u / (n - k)
gamma0 = s2 * (n - k) / n
sigma = res.bse[0]
rho = res.params[0]
stat_tau = np.sqrt(gamma0 / lam2) * ((rho - 1) / sigma) - 0.5 * (
    (lam2 - gamma0) / np.sqrt(lam2)) * (n * sigma / np.sqrt(s2))
w("// PP tiny 手算分量 (y=[1..10], trend n, lags=2)")
w(f"constexpr Real PP_TINY_RHO = {fmt(rho)};")
w(f"constexpr Real PP_TINY_BSE = {fmt(sigma)};")
w(f"constexpr Real PP_TINY_LAM2 = {fmt(lam2)};")
w(f"constexpr Real PP_TINY_S2 = {fmt(s2)};")
w(f"constexpr Real PP_TINY_GAMMA0 = {fmt(gamma0)};")
w(f"constexpr Real PP_TINY_STAT = {fmt(stat_tau)};")
r = PhillipsPerron(tiny_y, lags=2, trend="n")
w(f"constexpr Real PP_TINY_ARCHSTAT = {fmt(r.stat)};  // 交叉验证")
w()

# ---------------------------------------------------------------- KPSS
w("// --- KPSS (Hobijn autolag + legacy Schwert) ---")
w("struct KpssCase { Real stat, p, cv1, cv5, cv10; Size lags; };")
w("constexpr KpssCase KPSS_CASES[5] = {")
for series, sname in [(y_rw, "rw"), (y_ar, "ar")]:
    for tr in ["c", "ct"]:
        r = KPSS(series, trend=tr)
        cv = r.critical_values
        w(f"    {{{fmt(r.stat)}, {fmt(r.pvalue)}, {fmt(cv['1%'])}, "
          f"{fmt(cv['5%'])}, {fmt(cv['10%'])}, {int(r.lags)}}},  // {sname}/{tr}")
r = KPSS(y_rw, lags=-1, trend="c")  # legacy Schwert
cv = r.critical_values
w(f"    {{{fmt(r.stat)}, {fmt(r.pvalue)}, {fmt(cv['1%'])}, "
  f"{fmt(cv['5%'])}, {fmt(cv['10%'])}, {int(r.lags)}}},  // rw/c legacy")
w("};")
w()

# KPSS tiny 手算: y=[1..10], trend c, lags=2
z = np.ones((10, 1))
resk = OLS(tiny_y, z).fit()
uk = resk.resid
lam_k = cov_nw(uk, 2, demean=False)
s_cum = np.cumsum(uk)
kpss_tiny = 1 / (10 ** 2) * (s_cum ** 2).sum() / lam_k
w("// KPSS tiny 手算 (y=[1..10], trend c, lags=2)")
w(f"constexpr Real KPSS_TINY_LAM = {fmt(lam_k)};")
w(f"constexpr Real KPSS_TINY_STAT = {fmt(kpss_tiny)};")
r = KPSS(tiny_y, lags=2, trend="c")
w(f"constexpr Real KPSS_TINY_ARCHSTAT = {fmt(r.stat)};  // 交叉验证")
w()

# cov_nw 手算小例 (LRV 单元测试): u=[1..5], lags=2, demean=False
u5 = np.arange(1.0, 6.0)
w(f"constexpr Real COVNW_U5_L2 = {fmt(cov_nw(u5, 2, demean=False))};")
w(f"constexpr Real COVNW_U5_L0 = {fmt(cov_nw(u5, 0, demean=False))};")
w()

# ---------------------------------------------------------------- DF-GLS
w("// --- DF-GLS (AIC 自动 lag) ---")
w("struct DfglsCase { Real stat, p, cv1, cv5, cv10; Size nobs, lags; };")
w("constexpr DfglsCase DFGLS_CASES[4] = {")
for series, sname in [(y_rw, "rw"), (y_ar, "ar")]:
    for tr in ["c", "ct"]:
        r = DFGLS(series, trend=tr)
        cv = r.critical_values
        w(f"    {{{fmt(r.stat)}, {fmt(r.pvalue)}, {fmt(cv['1%'])}, "
          f"{fmt(cv['5%'])}, {fmt(cv['10%'])}, {int(r.nobs)}, {int(r.lags)}}},"
          f"  // {sname}/{tr}")
w("};")
w()

# ---------------------------------------------------------------- VarianceRatio
w("// --- VarianceRatio (robust/debiased 组合, trend=c) ---")
w("struct VrCase { Real vr, stat, p; };")
combos = [(True, True), (True, False), (False, True)]
w(f"constexpr VrCase VR_CASES[{2 * 2 * 3}] = {{")
for series, sname in [(p_rw, "rw"), (p_ar, "ar")]:
    for k in [2, 5]:
        for (rob, deb) in combos:
            r = VarianceRatio(series, lags=k, trend="c",
                              debiased=deb, robust=rob)
            w(f"    {{{fmt(r.vr)}, {fmt(r.stat)}, {fmt(r.pvalue)}}},"
              f"  // {sname}/k={k}/rob={rob}/deb={deb}")
w("};")
# trend="n" 变体
w("constexpr VrCase VR_N_CASES[2] = {")
for k in [2, 5]:
    r = VarianceRatio(p_rw, lags=k, trend="n", debiased=True, robust=True)
    w(f"    {{{fmt(r.vr)}, {fmt(r.stat)}, {fmt(r.pvalue)}}},  // rw/k={k}/n")
w("};")
w()

# VR tiny 手算分量: P=[1..11], k=2, trend=c, robust+debiased
dy = np.diff(tiny_p)
nq = len(dy)
mu = (tiny_p[-1] - tiny_p[0]) / (len(tiny_p) - 1)
sig1 = np.sum((dy - mu) ** 2) / nq
dyq = tiny_p[2:] - tiny_p[:-2]
sigq = np.sum((dyq - 2 * mu) ** 2) / (nq * 2)
sig1_d = sig1 * nq / (nq - 1)
m = 2 * (nq - 2 + 1) * (1 - 2 / nq)
sigq_d = sigq * (nq * 2) / m
z2 = (dy - mu) ** 2
scale = np.sum(z2) ** 2
theta = 0.0
for kk in range(1, 2):
    delta = nq * z2[kk:] @ z2[:-kk] / scale
    theta += 4 * (1 - kk / 2) ** 2 * delta
vr_d = sigq_d / sig1_d
stat_d = np.sqrt(nq) * (vr_d - 1) / np.sqrt(theta)
vr_h = sigq / sig1
var_h = (2 * (2 * 2 - 1) * (2 - 1)) / (3 * 2)
stat_h = np.sqrt(nq) * (vr_h - 1) / np.sqrt(var_h)
w("// VR tiny 手算分量 (P=[1..11], k=2, trend=c)")
w(f"constexpr Real VR_TINY_MU = {fmt(mu)};")
w(f"constexpr Real VR_TINY_SIG1 = {fmt(sig1)};")
w(f"constexpr Real VR_TINY_SIGQ = {fmt(sigq)};")
w(f"constexpr Real VR_TINY_VR = {fmt(vr_h)};")
w(f"constexpr Real VR_TINY_VARH = {fmt(var_h)};")
w(f"constexpr Real VR_TINY_STATH = {fmt(stat_h)};")
w(f"constexpr Real VR_TINY_THETA = {fmt(theta)};")
w(f"constexpr Real VR_TINY_VRD = {fmt(vr_d)};")
w(f"constexpr Real VR_TINY_STATD = {fmt(stat_d)};")
r = VarianceRatio(tiny_p, lags=2, trend="c", debiased=True, robust=True)
w(f"constexpr Real VR_TINY_ARCHVR = {fmt(r.vr)};  // 交叉验证")
w(f"constexpr Real VR_TINY_ARCHSTAT = {fmt(r.stat)};  // 交叉验证")
r2 = VarianceRatio(tiny_p, lags=2, trend="c", debiased=False, robust=False)
w(f"constexpr Real VR_TINY_ARCHVR_H = {fmt(r2.vr)};  // 交叉验证 (homosk)")
w(f"constexpr Real VR_TINY_ARCHSTAT_H = {fmt(r2.stat)};  // 交叉验证 (homosk)")
w()

# ---------------------------------------------------------------- MacKinnon CV
w("// --- MacKinnon CV/P 参考值 ---")
for tr in ["n", "c", "ct"]:
    for nn in [233, 100]:
        cva = mackinnoncrit(regression=tr, nobs=nn)
        w(f"constexpr Real MKCV_ADF_{tr.upper()}_{nn}[] = "
          f"{{{fmt(cva[0])}, {fmt(cva[1])}, {fmt(cva[2])}}};")
    cva = mackinnoncrit(regression=tr, nobs=np.inf)
    w(f"constexpr Real MKCV_ADF_{tr.upper()}_INF[] = "
      f"{{{fmt(cva[0])}, {fmt(cva[1])}, {fmt(cva[2])}}};")
for tr in ["c", "ct"]:
    for nn in [233, 100]:
        cva = mackinnoncrit(regression=tr, nobs=nn, dist_type="dfgls")
        w(f"constexpr Real MKCV_DFGLS_{tr.upper()}_{nn}[] = "
          f"{{{fmt(cva[0])}, {fmt(cva[1])}, {fmt(cva[2])}}};")
    cva = mackinnoncrit(regression=tr, nobs=np.inf, dist_type="dfgls")
    w(f"constexpr Real MKCV_DFGLS_{tr.upper()}_INF[] = "
      f"{{{fmt(cva[0])}, {fmt(cva[1])}, {fmt(cva[2])}}};")
w()
for st in [-3.5, -2.5, -1.0]:
    w(f"constexpr Real MKP_ADF_C_T233_M{tag(st)} = "
      f"{fmt(mackinnonp(st, regression='c', num_unit_roots=1))};")
for st in [-3.5, -2.0]:
    w(f"constexpr Real MKP_DFGLS_C_M{tag(st)} = "
      f"{fmt(mackinnonp(st, regression='c', dist_type='dfgls'))};")
w()

# ---------------------------------------------------------------- KPSS crit
w("// --- KPSS p 值插值 ---")
pv, cvv = kpss_crit(0.5, "c")
w(f"constexpr Real KPSSCRIT_C_P05 = {fmt(pv)};")
w(f"constexpr Real KPSSCRIT_C_CV[] = "
  f"{{{fmt(cvv[0])}, {fmt(cvv[1])}, {fmt(cvv[2])}}};")
pv, cvv = kpss_crit(0.3, "ct")
w(f"constexpr Real KPSSCRIT_CT_P03 = {fmt(pv)};")
w(f"constexpr Real KPSSCRIT_CT_CV[] = "
  f"{{{fmt(cvv[0])}, {fmt(cvv[1])}, {fmt(cvv[2])}}};")
pv, cvv = kpss_crit(10.0, "c")   # 超界 → p=0.01 (表下限)
w(f"constexpr Real KPSSCRIT_C_P10 = {fmt(pv)};  // 超上界截断")
w()

w("}  // namespace baseline")
w("}  // namespace unit_root")
w("}  // namespace timeseries")
w("}  // namespace v1")
w("}  // namespace cpphub")

with open(OUT, "w", encoding="utf-8", newline="\n") as f:
    f.write("\n".join(L) + "\n")

print("written:", OUT)
print("Schwert(250) =", SCHWERT_250)
print("sanity: ADF rw/c fixed:", ADF(y_rw, lags=16, trend="c").stat,
      "p=", ADF(y_rw, lags=16, trend="c").pvalue)
print("sanity: ADF ar/c fixed:", ADF(y_ar, lags=16, trend="c").stat)
print("sanity: PP tiny cross-check:", stat_tau, "vs arch", r.stat if False else
      PhillipsPerron(tiny_y, lags=2, trend="n").stat)
print("sanity: KPSS tiny:", kpss_tiny, "vs arch", KPSS(tiny_y, lags=2, trend="c").stat)
print("sanity: VR tiny deb:", stat_d, "vs arch",
      VarianceRatio(tiny_p, lags=2).stat)
