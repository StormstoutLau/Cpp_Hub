# -*- coding: utf-8 -*-
# probe_sm_innovations.py — statsmodels innovations loglik 复杂度探针 (C-4/ARIMA, 2026-08-20)
# 目的: 对照 C++ innovations_algo O(T^3) 退化嫌疑 — 实测 statsmodels Cython
#       arma_innovations.arma_loglike (innovations_mle 的目标函数, 同一路径)
#       在 T=300/600/1000 的单次调用耗时, 确认基准侧复杂度
#       (O(T*(p+q)^2) vs O(T^3))
import time

import numpy as np
from statsmodels.tsa.innovations import arma_innovations

print("module:", arma_innovations.__file__)
rng = np.random.RandomState(42)

prev = None
for T in (300, 600, 1000):
    y = rng.standard_normal(T)
    ar = np.array([-0.5, 0.3])   # ARMA(2,1)
    ma = np.array([0.4])
    # 预热 (首次含 import/编译开销)
    arma_innovations.arma_loglike(y, ar_params=ar, ma_params=ma, sigma2=1.0)
    n = 200
    t0 = time.perf_counter()
    for _ in range(n):
        arma_innovations.arma_loglike(y, ar_params=ar, ma_params=ma, sigma2=1.0)
    per = (time.perf_counter() - t0) / n * 1000
    ratio = f" | x{per/prev:.2f} vs prev" if prev else ""
    print(f"T={T:5d}: arma_loglike {per:9.4f} ms/call{ratio}")
    prev = per
