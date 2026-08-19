# -*- coding: utf-8 -*-
# gen_em2002_tables.py - EM2002 (IFDP 655) 响应面表精确转录
#
# 源: Ericsson-MacKinnon 2002 IFDP 655 (1999 WP 版, Fed 官方免费 PDF)
#   Table 2 (p31): 无确定项; Table 3 (p32): 常数; Table 4 (p33): 常数+趋势;
#   Table 5 (p34): 常数+趋势+二次趋势
# PDF 字体 cmr10092 自定义编码: cipher = chr(plain_ord + 3), '~' = minus
# 公式 (经表 7 有限样本 CV 在 T=51 精确复现验证): CV(T) = θ∞ + θ1/T + θ2/T² + θ3/T³
#
# 三重验证:
#   V1: n=1 行 θ∞ vs MacKinnon 2010 (tau_2010 表 N=1) ≤ 0.001
#   V2: 行结构 12 n × 3 sizes × 4 系数 + s.e. + S.E. 完整
#   V3: 表 7 K_ctt(3) 有限 CV (T=51) 复现 (内嵌断言)
#
# 输出: em2002_ect_cv.csv (case, n, size, theta_inf, se, theta1, theta2, theta3, reg_se)

import csv

import fitz

PDF = r"F:\Cpp_Hub\tests\fixtures\timeseries\em2002_ifdp655.pdf"
OUT = r"F:\Cpp_Hub\tests\fixtures\timeseries\em2002_ect_cv.csv"

CIPHER = {}
for plain_ord in range(33, 127):
    CIPHER[chr(plain_ord + 3)] = chr(plain_ord)


def dec_word(w):
    out = []
    for ch in w:
        if ch == '~':
            out.append('-')
        else:
            out.append(CIPHER.get(ch, ch))
    return "".join(out)


def is_num(s):
    if not s:
        return False
    body = s.lstrip('-').rstrip('.')
    return body.replace('.', '', 1).isdigit() and body.count('.') <= 1


PAGES = {  # PDF 页号 (0-based) → case
    30: "n",
    31: "c",
    32: "ct",
    33: "ctt",
}

rows = []
doc = fitz.open(PDF)
for pno, case in PAGES.items():
    page = doc[pno]
    words = page.get_text("words")  # (x0, y0, x1, y1, word, block, line, wordno)
    # 按行分组 (y0 圆整), 行内按 x 排序
    lines = {}
    for w in words:
        ykey = round(w[1] / 3.0)  # 3pt 容差
        lines.setdefault(ykey, []).append(w)
    n_cur = None
    for ykey in sorted(lines):
        ws = sorted(lines[ykey], key=lambda w: w[0])
        toks = [dec_word(w[4]) for w in ws]
        toks = [t for t in toks if t.strip()]
        if not toks:
            continue
        # 数据行: [n] size θ∞ (s.e.) θ1 θ2 θ3 S.E.  或  size θ∞ (s.e.) θ1 θ2 θ3 S.E.
        if len(toks) >= 7 and toks[0].isdigit() and toks[1] in ("1%", "5%", "10%"):
            n_cur = int(toks[0])
            rest = toks[1:]
        elif len(toks) >= 6 and toks[0] in ("1%", "5%", "10%"):
            rest = toks
            if n_cur is None:
                continue
        else:
            continue
        size = rest[0]
        nums = rest[1:]
        # θ∞ (s.e.) θ1 θ2 θ3 S.E. → 6 个数值 token (s.e. 带括号)
        if len(nums) != 6:
            continue
        the = nums[0]
        se = nums[1].strip("()")
        t1, t2, t3, rse = nums[2], nums[3], nums[4], nums[5]
        vals = [the, se, t1, t2, t3, rse]
        if not all(is_num(v) for v in [the, t1, t2, t3]):
            continue
        rows.append([case, n_cur, size] + [float(v) for v in vals])
doc.close()

# ---- V2 结构验证 ----
expect = {"n": [], "c": [], "ct": [], "ctt": []}
for r in rows:
    expect[r[0]].append((r[1], r[2]))
for case, keys in expect.items():
    assert len(keys) == 36, (case, len(keys))
    want = [(n, s) for n in range(1, 13) for s in ("1%", "5%", "10%")]
    assert keys == want, (case, keys[:6])

# ---- V1 n=1 行 vs MacKinnon 2010 ----
M2010_N1 = {  # trend → {size: θ∞}
    "n": {"1%": -2.56574, "5%": -1.94100, "10%": -1.61682},
    "c": {"1%": -3.43035, "5%": -2.86154, "10%": -2.56677},
    "ct": {"1%": -3.95877, "5%": -3.41049, "10%": -3.12705},
    "ctt": {"1%": -4.37113, "5%": -3.83239, "10%": -3.55326},
}
print("V1: n=1 行 θ∞ vs MacKinnon 2010 (|Δ| 应 ≤ 0.001):")
for r in rows:
    if r[0] in M2010_N1 and r[1] == 1:
        ref = M2010_N1[r[0]][r[2]]
        diff = abs(r[3] - ref)
        assert diff < 0.001, (r, ref, diff)
        print("  %-4s %-3s θ∞=%9.4f  M2010=%9.5f  |Δ|=%.4f" %
              (r[0], r[2], r[3], ref, diff))

# ---- V3 表 7 K_ctt(3) T=51 有限样本 CV 复现 ----
ctt3 = {r[2]: r for r in rows if r[0] == "ctt" and r[1] == 3}
T = 51.0
finites = {"1%": -5.09, "5%": -4.38, "10%": -4.03}
print("\nV3: 表 7 K_ctt(3) 有限样本 CV (T=51) 复现:")
for size, want in finites.items():
    r = ctt3[size]
    cv = (r[3] + r[5] / T + r[6] / T**2 + r[7] / T**3)
    print("  %-3s 计算=%8.4f 论文表7=%6.2f  (舍入一致: %s)" %
          (size, cv, want, abs(round(cv, 2) - want) < 0.005))
    assert abs(round(cv, 2) - want) < 0.005, (size, cv, want)

with open(OUT, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["case", "n", "size", "theta_inf", "se_theta_inf",
                "theta1", "theta2", "theta3", "reg_se"])
    w.writerows(rows)
print("\nwrote", OUT, "(%d rows)" % len(rows))
