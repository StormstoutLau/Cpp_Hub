# -*- coding: utf-8 -*-
# probe_em2002_pdf.py - EM2002 (IFDP 655) 全页文本转储 (表定位用)

PDF = r"F:\Cpp_Hub\tests\fixtures\timeseries\em2002_ifdp655.pdf"
OUT = r"F:\Cpp_Hub\tests\fixtures\timeseries\em2002_pages_raw.txt"

try:
    from pypdf import PdfReader
except ImportError:
    from PyPDF2 import PdfReader

reader = PdfReader(PDF)
print("pages:", len(reader.pages))

with open(OUT, "w", encoding="utf-8") as fh:
    for i, page in enumerate(reader.pages):
        text = page.extract_text() or ""
        fh.write("\n===== PAGE %d =====\n" % (i + 1))
        fh.write(text)
print("wrote", OUT)
