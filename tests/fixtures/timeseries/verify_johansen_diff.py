# -*- coding: utf-8 -*-
# verify_johansen_diff.py - statsmodels 侧 Johansen 双库 diff 网格 (M3 前置, §6.2.1)
#
# 冻结结论见 docs/phases/phase7/JOHANSEN_DUAL_LIB_DIFF.md:
#   - 参数映射: urca ca.jo(K) ≙ coint_johansen(k_ar_diff = K-1)
#   - 情形映射: SM det_order=0 ≡ urca ecdet="none" (短回归无约束常数); det -1/1 无 urca 对应
#   - urca spec 参数数学恒等 (y_{t-1}-y_{t-K} ∈ span(Z1))
#
# 网格: det_order ∈ {-1, 0, 1} × k_ar_diff ∈ {1, 2} → johansen_sm_grid.txt
# 配对分析: probe_johansen_pair.py (读双 grid, 输出匹配对)

import os
import numpy as np
from statsmodels.tsa.vector_ar.vecm import coint_johansen

HERE = os.path.dirname(os.path.abspath(__file__))


def emit(fh, tag, obj):
    a = np.asarray(obj, dtype=float).ravel()
    fh.write("%s\t%s\n" % (tag, " ".join("%.17g" % v for v in a)))


def main():
    y = np.loadtxt(os.path.join(HERE, "coint_smoke_data.csv"),
                   delimiter=",", skiprows=1)
    out = os.path.join(HERE, "johansen_sm_grid.txt")
    with open(out, "w") as fh:
        for det in (-1, 0, 1):
            for k in (1, 2):
                tag = "sm_det%d_k%d" % (det, k)
                r = coint_johansen(y, det, k)
                emit(fh, tag + "_eig", r.eig)
                emit(fh, tag + "_lr1", r.lr1)   # trace, r=0..N-1
                emit(fh, tag + "_lr2", r.lr2)   # maxeig, r=0..N-1
                emit(fh, tag + "_cvt", r.cvt)   # N×3 行序 r=0..N-1 (行=rank, 列=90/95/99)
                emit(fh, tag + "_cvm", r.cvm)
                emit(fh, tag + "_evec", r.evec)
                fh.write("%s_nobs\t%d\n" % (tag, r.rkt.shape[0]))
    print("wrote", out)


if __name__ == "__main__":
    main()
