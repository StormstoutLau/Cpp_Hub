# -*- coding: utf-8 -*-
# probe_granger_layout.py - lagmat/lagmat2ds 布局探测 (一次性, 不留基准)
import numpy as np
from statsmodels.tsa.tsatools import lagmat, lagmat2ds

x = np.arange(1, 8, dtype=float)  # T=7
m = lagmat(x, 2, trim="both", original="in")
print("lagmat(x,2,both,in) shape:", m.shape)
print(m)
X = np.column_stack([x, x[::-1]])
d = lagmat2ds(X, 2, trim="both", dropex=1)
print("lagmat2ds shape:", d.shape)
print(d)
