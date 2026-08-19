# -*- coding: utf-8 -*-
# verify_vecm.py - VECM 基准生成: statsmodels VECM (Phase 7C M3)
#
# 对照点 (spec §5.4):
#   alpha/beta/gamma/det_coef/det_coef_coint/sigma_u/llf vs statsmodels 1e-10
#   β 对照恒用投影矩阵 P = β(β'β)⁻¹β' (CI8; 本文件同时输出归一后 β 逐元素值
#   —— statsmodels β 前 r 行 = I_r 归一, r=1 时逐元素可直接对照)
#
# 语义 (一手源码 vecm.py):
#   _endog_matrices (L259-372): p = k_ar_diff+1;
#     y_lag1 = y[p-1:-1] + [1 (ci)] + [arange(T)+p (li, coint)]
#     delta_x = Δy 滞后 + [1 (co)] + [arange(T)+p+1 (lo, 外部)]
#   _sij (L415-458): R0/R1 (m = I - ΔX'(ΔXΔX')⁻¹ΔX); λ = eig(s11_·s10·s00⁻¹·s01·s11_)
#   β̃ = s11_·V_r → β̃·inv(β̃[:r]) (前 r 行 = I_r, L1017)
#   α = s01·β̃;  Γ = (Δy-αβ̃'y_lag1)·ΔX'(ΔXΔX')⁻¹;  Σ = resid'resid/T
#   llf (L1471): -KT·ln(2π)/2 - T·(ln|s00| + Σ_{i<r} ln(1-λi))/2 - KT/2
#   β 拆分 (L1418): [beta(K×r); det_coef_coint(#det_coint×r)]
#   Γ 拆分 (L1419): [gamma(K×K·k), det_coef(K×#det_out)]
#
# 输出: vecm_statsmodels_baselines.txt

import os

import numpy as np
from statsmodels.tsa.vector_ar.vecm import VECM

HERE = os.path.dirname(os.path.abspath(__file__))


def emit(fh, tag, obj):
    a = np.asarray(obj, dtype=float).ravel()
    fh.write("%s\t%s\n" % (tag, " ".join("%.17g" % v for v in a)))


def main():
    y = np.loadtxt(os.path.join(HERE, "coint_smoke_data.csv"),
                   delimiter=",", skiprows=1)
    out = os.path.join(HERE, "vecm_statsmodels_baselines.txt")
    with open(out, "w") as fh:
        for det in ("n", "co", "ci", "lo", "li"):
            for rank in (1, 2):
                for k in (1, 2):
                    base = "vecm_%s_r%d_k%d" % (det, rank, k)
                    res = VECM(y, k_ar_diff=k, coint_rank=rank,
                               deterministic=det).fit()
                    emit(fh, base + "_alpha", res.alpha)
                    emit(fh, base + "_stderr_alpha", res.stderr_alpha)  # ECT t 基准 (CI10)
                    emit(fh, base + "_beta", res.beta)              # K×r, 前 r 行 = I_r
                    emit(fh, base + "_gamma", res.gamma)
                    emit(fh, base + "_det_coef", res.det_coef)      # K×#det_out (co/lo)
                    emit(fh, base + "_det_coef_coint",
                         res.det_coef_coint)                        # #det_coint×r (ci/li)
                    emit(fh, base + "_sigma_u", res.sigma_u)
                    emit(fh, base + "_llf", res.llf)
                    # resid (T_eff×K) 前 4 个时点, 行主序: [t0v0 t0v1 t0v2 t1v0 ...]
                    emit(fh, base + "_resid_head", res.resid[:4])
                    fh.write("%s_nobs\t%d\n" % (base, res.nobs))
    print("wrote", out)


if __name__ == "__main__":
    main()
