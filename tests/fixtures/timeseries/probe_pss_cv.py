# probe_pss_cv.py — M1 前置任务 (PHASE8_SPEC §2.4 Step 1 / §10)
# 目的: 探明 statsmodels 0.14.6 pss_critical_values.py 的
#   (a) 系数文件布局 (键空间 / 多项式阶 / 求值输入)
#   (b) crit_vals 与 p_values 的生成路径 (表常量 vs 多项式求值)
#   (c) bounds_test 内部 k 语义 (nvar = len(ardl_order) 是否 = k_PSS)
#   (d) t 检验存在性 (statsmodels 是否有 bounds_t_test)
# 用法: & C:\Users\Peng\.conda\envs\open-webui\python.exe probe_pss_cv.py
import numpy as np
from statsmodels.tsa.ardl import ARDL, UECM
from statsmodels.tsa.ardl import pss_critical_values as pcv
from statsmodels.datasets import danish_data

print("=" * 70)
print("[A] 键空间布局")
keys = sorted(pcv.large_p.keys())
ks = sorted({k for k, c, i in keys})
cases = sorted({c for k, c, i in keys})
i1s = sorted({i for k, c, i in keys})
print(f"large_p: {len(pcv.large_p)} 键; k ∈ {ks}; case ∈ {cases}; I1 ∈ {i1s}")
print(f"small_p: {len(pcv.small_p)} 键; stat_star: {len(pcv.stat_star)} 键; crit_vals: {len(pcv.crit_vals)} 键")
print(f"crit_percentiles = {pcv.crit_percentiles}")
print(f"large_p 系数长度: {sorted({len(v) for v in pcv.large_p.values()})}")
print(f"small_p 系数长度: {sorted({len(v) for v in pcv.small_p.values()})}")
print(f"crit_vals 值长度: {sorted({len(v) for v in pcv.crit_vals.values()})}")
n_digits = sum(len(v) for v in pcv.large_p.values()) + sum(len(v) for v in pcv.small_p.values()) \
    + len(pcv.stat_star) + sum(len(v) for v in pcv.crit_vals.values())
print(f"转录总数字量: {n_digits}")

print("=" * 70)
print("[B] 求值语义 (源码 model.py L2447-2455 _pss_pvalue 实录)")
print("  key = (k, case, i1); log_stat = log(stat)")
print("  p = small_p if stat > stat_star else large_p  (注意: 阈值比较用原始 stat)")
print("  x = [log_stat**i for i in range(len(p))]  →  p_val = 1 - norm.cdf(x @ p)")
# 数值自证: 手工重算 _pss_pvalue
from statsmodels.tsa.ardl.model import _pss_pvalue
for (k, case, i1), expect in [((1, 3, False), None), ((1, 3, True), None)]:
    for stat in [2.0, 4.9228967, 5.7217442, 10.0]:
        ref = _pss_pvalue(stat, k, case, i1)
        key = (k, case, i1)
        thr = pcv.stat_star[key]
        p = pcv.small_p[key] if stat > thr else pcv.large_p[key]
        x = np.log(stat)
        manual = 1 - 0.5 * (1 + __import__("math").erf((np.polyval(p[::-1], x)) / np.sqrt(2)))
        print(f"  stat={stat:9.4f} key={key}: statsmodels={ref:.10f} 手工正态CDF={manual:.10f} diff={abs(ref-manual):.2e}")

print("=" * 70)
print("[C] bounds_test 内部 k 语义 (danish_data, k_PSS=3)")
data = danish_data.load_pandas().data
lrm = data.lrm
exog = data[["lry", "ibo", "ide"]]
ardl = ARDL(lrm, 2, exog, 1, trend="c", causal=False).fit()
uecm_res = UECM.from_ardl(ardl.model).fit()
bt = uecm_res.bounds_test(case=3)
print(f"stat = {bt.stat:.6f}")
print(f"crit_vals:\n{bt.crit_vals}")
print(f"p_values:\n{bt.p_values}")
res = uecm_res
print(f"uecm ardl_order = {res.model.ardl_order}  (len={len(res.model.ardl_order)})")
print(f"nvar = k_key = {len(res.model.ardl_order)}  vs k_PSS = 3 (lry/ibo/ide)")
print(f"uecm param names = {list(res.params.index)}")
# 对照: 手工用 k=3 与 k=4 键查 CV, 看哪个与 bounds_test 返回一致
for k_try in [3, 4]:
    lo = pcv.crit_vals[(k_try, 3, False)]
    up = pcv.crit_vals[(k_try, 3, True)]
    print(f"  手工键 (k={k_try},3): 5% lower={lo[1]:.7f} upper={up[1]:.7f}")

print("=" * 70)
print("[D] t 检验存在性")
print("statsmodels ARDL/UECM bounds_test: 仅 F (Wald/q); 无 bounds_t_test —",
      not hasattr(ardl, "bounds_t_test"))
import statsmodels.tsa.ardl as ardl_mod
print("ardl 模块导出:", [n for n in dir(ardl_mod) if "bound" in n.lower() or "t_test" in n.lower()])

print("=" * 70)
print("[E] k 键语义判定 (核心发现)")
print("Case III 5% CV 全 k 序列 (statsmodels 表, 键 (k,3,I1)):")
for k in range(1, 11):
    lo = pcv.crit_vals[(k, 3, False)][1]
    up = pcv.crit_vals[(k, 3, True)][1]
    print(f"  k={k:2d}: I(0)={lo:.7f}  I(1)={up:.7f}")
print("已知锚: PSS 原文 Table CI(iii) k=1 5% = 4.94/5.73 (Statalist Kripfganz 帖子旁证)")
print("→ statsmodels 表键 (1,3)=4.9229/5.7217 ≈ PSS k=1 ⇒ 表键语义 = PSS k (force 变量数)")
print("→ 但 bounds_test 内部 k = nvar = len(ardl_order) = 1 + k_PSS (含 y 的水平项总数)")
print(f"→ danish 实测: k_PSS=3 (lry/ibo/ide) → bounds_test 用键 (4,3) [CV 2.8772/4.0105]")
print(f"   而 PSS 论文 danish 案例 (Table 7) 用 k=3 [statsmodels 键 (3,3): 3.2290/4.3223]")
print("⇒ 【发现】statsmodels bounds_test 查表 k 与其表键语义 (PSS k) 偏移 +1;")
print("   R ARDL 0.2.5 bounds_f_test 的 k = number of independent variables (PSS 语义,")
print("   rdrr.io 文档 null.value 字段实录) — 三方中 statsmodels 调用层是异类")
print("⇒ F 统计量计算不依赖 k 键 (Wald 只用系数+vcov+约束矩阵), stat 对照不受影响;")
print("   仅 CV/p_values 查表口径需裁决: PSS 语义 (k_PSS) vs 镜像 statsmodels (1+k_PSS)")

print("=" * 70)
print("[F] 裁决落档 (2026-08-21, Scott 确认推荐项)")
print("k 查表口径 = PSS 语义 (k = force 变量数):")
print("  1. bounds_f/t_test 内部查表 k = k_PSS (与 PSS 原文/R ARDL 一致)")
print("  2. pss_crit_values(k, case, i1) 接口的 k 即 force 变量数 (表键语义直译)")
print("  3. 对照矩阵: stat vs statsmodels 1e-10 (不受 k 影响) / CV+p vs R ARDL")
print("     (1e-8) + vs statsmodels 表键直查 (精确) — statsmodels bounds_test")
print("     返回的 CV/p 不逐值对照, 文档声明其 k 偏移 +1 (PSS10 新幻觉点)")
print("  4. 依据: ADR-016 决策原则 (2) 理论正确性优先 + PSS6 定义 + R ARDL 实证")

print("=" * 70)
print("[G] .inc 全量 diff (生成后回读, 1200 数字独立往返校验)")
import re
from pathlib import Path

inc = Path(__file__).resolve().parents[3] / "include" / "cpphub" / "timeseries" / "ardl" / "pss_response_surface.inc"
text = inc.read_text(encoding="utf-8")

def parse_block(name: str, ncols: int):
    m = re.search(rf"PSS_{name}\[?\]?[^=]*=\s*\{{(.*?)\n\}};", text, re.S)
    rows = []
    for row_m in re.finditer(r"\{([^{}]*)\}|^    ([\d.eE+-]+,?\s*)$",
                             m.group(1), re.S | re.M):
        s = row_m.group(1) or row_m.group(2)
        vals = [float(v) for v in s.replace(",", " ").split()]
        if vals:
            rows.append(vals)
    return rows

def parse_star():
    m = re.search(r"PSS_STAT_STAR\[\]\s*=\s*\{(.*?)\n\};", text, re.S)
    return [float(v) for v in m.group(1).replace(",", " ").split()]

inc_large = parse_block("LARGE_P", 4)
inc_small = parse_block("SMALL_P", 3)
inc_star = parse_star()
inc_crit = parse_block("CRIT_VALS", 4)
print(f"  回读行数: large={len(inc_large)} small={len(inc_small)} "
      f"star={len(inc_star)} crit={len(inc_crit)} (各预期 100)")

def idx(k, case, i1):
    return (k - 1) * 10 + (case - 1) * 2 + (1 if i1 else 0)

diffs = 0
for (k, case, i1), vals in pcv.large_p.items():
    if [float(v) for v in inc_large[idx(k, case, i1)]] != [float(v) for v in vals]:
        diffs += 1
for (k, case, i1), vals in pcv.small_p.items():
    if [float(v) for v in inc_small[idx(k, case, i1)]] != [float(v) for v in vals]:
        diffs += 1
for (k, case, i1), v in pcv.stat_star.items():
    if float(inc_star[idx(k, case, i1)]) != float(v):
        diffs += 1
for (k, case, i1), vals in pcv.crit_vals.items():
    if [float(v) for v in inc_crit[idx(k, case, i1)]] != [float(v) for v in vals]:
        diffs += 1
n_total = len(pcv.large_p) * 4 + len(pcv.small_p) * 3 + len(pcv.stat_star) + len(pcv.crit_vals) * 4
print(f"  逐值 diff: {diffs} 处不一致 / {n_total} 数字 (预期 0)")
if diffs == 0 and n_total == 1200:
    print("  ✅ 全量 diff 通过 — .inc 与 statsmodels 逐位一致 (风险 #1 钉死)")
else:
    print("  ❌ DIFF FAILED — 检查生成脚本")
    raise SystemExit(1)
