# -*- coding: utf-8 -*-
# probe_em2002_fonts.py - θ3 列尾部字形字体检查 (vs θ∞ 小数点)

import fitz

PDF = r"F:\Cpp_Hub\tests\fixtures\timeseries\em2002_ifdp655.pdf"

doc = fitz.open(PDF)
page = doc[33]  # PDF 页 34 = Table 5 (ctt)
raw = page.get_text("rawdict")
for block in raw["blocks"]:
    if block["type"] != 0:
        continue
    for line in block["lines"]:
        for span in line["spans"]:
            chars = "".join(c["c"] for c in span["chars"])
            # θ∞ n=1 1% = ~716:47; θ3 n=1 1% = ~991; S.E. = 3133;7<
            if chars.strip() in ("~716:47", "~991", "3133;7<", ":17", "~4418:"):
                fonts = {}
                for c in span["chars"]:
                    fonts.setdefault(c["c"], []).append(
                        (span["font"], round(c["origin"][0], 1)))
                print(repr(chars), "font:", span["font"])
                for ch, occ in fonts.items():
                    print("   char %r x=%s" % (ch, occ))
doc.close()
