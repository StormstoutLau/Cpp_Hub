# -*- coding: utf-8 -*-
# probe_em2002_formula.py - EM2002 全文解码定位公式段

import fitz

PDF = r"F:\Cpp_Hub\tests\fixtures\timeseries\em2002_ifdp655.pdf"

CIPHER = {}
for plain_ord in range(33, 127):
    CIPHER[chr(plain_ord + 3)] = chr(plain_ord)


def dec(s):
    out = []
    for ch in s:
        if ch == '~':
            out.append('-')
        else:
            out.append(CIPHER.get(ch, ch))
    return "".join(out)


doc = fitz.open(PDF)
for pno in range(len(doc)):
    d = dec(doc[pno].get_text())
    if "esponse" in d or "urfaces" in d:
        # 打印含关键词的行及上下文
        lines = d.split("\n")
        for i, ln in enumerate(lines):
            if "esponse" in ln or "urfaces" in ln:
                lo = max(0, i - 12)
                hi = min(len(lines), i + 12)
                print("===== PDF PAGE %d, 行 %d 附近 =====" % (pno + 1, i))
                print("\n".join(lines[lo:hi]))
                print()
doc.close()
