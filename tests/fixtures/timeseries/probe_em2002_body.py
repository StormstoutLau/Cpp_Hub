# -*- coding: utf-8 -*-
# probe_em2002_body.py - EM2002 正文 (p15-20) 解码找响应面公式

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
for pno in range(14, 20):  # PDF 页 15-20
    d = dec(doc[pno].get_text())
    print("===== PDF PAGE %d =====" % (pno + 1))
    print(d)
doc.close()
