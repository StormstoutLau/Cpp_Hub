# -*- coding: utf-8 -*-
# render_em2002.py - EM2002 IFDP655 表页: PyMuPDF 文本提取 (Table 2-5 响应面)

import fitz

PDF = r"F:\Cpp_Hub\tests\fixtures\timeseries\em2002_ifdp655.pdf"

doc = fitz.open(PDF)
for pno in (30, 31, 32, 33):  # 0-based: PDF 页 31-34
    page = doc[pno]
    print("\n===== PDF PAGE %d =====" % (pno + 1))
    print(page.get_text())
doc.close()
