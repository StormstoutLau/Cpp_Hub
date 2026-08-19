# -*- coding: utf-8 -*-
# verify_johansen.py - Johansen 基准生成: statsmodels coint_johansen (Phase 7C M3)
#
# 对照点 (spec §5.2 + JOHANSEN_DUAL_LIB_DIFF.md 冻结决策):
#   eig/lr1/lr2/cvt/cvm/evec vs statsmodels 1e-10, 全 det_order ∈ {-1,0,1} × k ∈ {1,2}
#   urca 交叉 (det_order=0 ↔ ecdet="none", k_ar_diff = K-1): 由 verify_johansen_diff.R 承担
#   select_coint_rank (trace/maxeig × signif) 逐级检验复刻
#
# 语义 (一手源码 vecm.py L603-737, 详见 diff 报告附录 B):
#   detrend(endog, det_order) → dx/lagmat/trim → detrend(·, f) → resid(dx,z)/resid(lx,z)
#   eig(inv(skk)·sk0·inv(s00)·sk0') + Cholesky β'S11β=I 归一 + 首非零元符号
#   lr1[r] = -t·Σ_{i>r} ln(1-λi); lr2[r] = -t·ln(1-λ_{r+1}); t = T-1-k
#   cvt[r,:] = c_sjt(N-r, det_order)  (MHM96, 90/95/99)
#
# 输出: johansen_statsmodels_baselines.txt

import os

import numpy as np
from statsmodels.tsa.vector_ar.vecm import coint_johansen, select_coint_rank

HERE = os.path.dirname(os.path.abspath(__file__))


def emit(fh, tag, obj):
    a = np.asarray(obj, dtype=float).ravel()
    fh.write("%s\t%s\n" % (tag, " ".join("%.17g" % v for v in a)))


def main():
    y = np.loadtxt(os.path.join(HERE, "coint_smoke_data.csv"),
                   delimiter=",", skiprows=1)
    out = os.path.join(HERE, "johansen_statsmodels_baselines.txt")
    with open(out, "w") as fh:
        for det in (-1, 0, 1):
            for k in (1, 2):
                base = "jo_det%d_k%d" % (det, k)
                r = coint_johansen(y, det, k)
                emit(fh, base + "_eig", r.eig)              # λ 降序
                emit(fh, base + "_lr1", r.lr1)              # trace r=0..N-1
                emit(fh, base + "_lr2", r.lr2)              # maxeig r=0..N-1
                emit(fh, base + "_cvt", r.cvt)              # N×3 (90/95/99), MHM96
                emit(fh, base + "_cvm", r.cvm)
                emit(fh, base + "_evec", r.evec)            # β'S11β=I 归一 (列符号任意)
                fh.write("%s_nobs\t%d\n" % (base, r.rkt.shape[0]))
                for method in ("trace", "maxeig"):
                    for signif in (0.1, 0.05, 0.01):
                        rk = select_coint_rank(y, det, k, method=method,
                                               signif=signif)
                        fh.write("%s_rank_%s_%g\t%d\n"
                                 % (base, method, signif, rk.rank))
    print("wrote", out)


if __name__ == "__main__":
    main()
