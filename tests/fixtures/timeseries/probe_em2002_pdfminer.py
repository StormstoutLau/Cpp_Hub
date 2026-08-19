# -*- coding: utf-8 -*-
# probe_em2002_pdfminer.py - pdfminer 提取 EM2002 表页 (对照 cipher 解码)

from pdfminer.high_level import extract_text

PDF = r"F:\Cpp_Hub\tests\fixtures\timeseries\em2002_ifdp655.pdf"
text = extract_text(PDF, page_numbers=[33])  # PDF 页 34 (Table 5, ctt)
print(text[:3000])
