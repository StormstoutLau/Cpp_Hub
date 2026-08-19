# -*- coding: utf-8 -*-
# probe_johansen_mc2.py - q=1 行的修正 MC (真 rank=1 DGP, 使 r<=1 检验处于原假设)
#
# probe_johansen_mc.py 的 q=1 行缺陷: DGP 真 rank=0 时测 r<=1, 非表分布设定。
# 此处 DGP: y2 = RW, y1 = y2 + e (真 rank=1, 协整向量 [1,-1]), 取 lr1[1]。
# 检验: SM det0 计算的 q=1 分布 ≈ urca none q=1 [6.5, 8.18, 11.65] 还是
#       SM tjcp1 q=1 = χ²(1) [2.7055, 3.8415, 6.6349]?

import numpy as np
from statsmodels.tsa.vector_ar.vecm import coint_johansen

NREP = 4000
T = 250


def mc_q1(det_order):
    rng = np.random.RandomState(42)
    s1 = np.empty(NREP)
    for i in range(NREP):
        e = rng.standard_normal((T, 2))
        y2 = np.cumsum(e[:, 1])
        y1 = y2 + e[:, 0]          # 真 rank=1
        y = np.column_stack([y1, y2])
        r = coint_johansen(y, det_order, 1)
        s1[i] = r.lr1[1]           # r=1 检验 → q=1
    return np.quantile(s1, [0.90, 0.95, 0.99])


def main():
    print("MC: y2=RW, y1=y2+e (真 rank=1), T=%d, %d reps, k=1, 取 lr1[1] (q=1)"
          % (T, NREP))
    print("经验分位数 [90%% 95%% 99%%]:")
    for det in (-1, 0, 1):
        q = mc_q1(det)
        print("  SM det=%+d: [%7.4f %7.4f %7.4f]" % (det, q[0], q[1], q[2]))
    print()
    print("对照表值 q=1:")
    print("  SM tjcp0 (det=-1 表): [ 2.9762  4.1296  6.9406]")
    print("  SM tjcp1 (det= 0 表): [ 2.7055  3.8415  6.6349]  (= χ²(1))")
    print("  SM tjcp2 (det= 1 表): [ 2.7055  3.8415  6.6349]  (= χ²(1))")
    print("  urca none  q=1:       [ 6.5000  8.1800 11.6500]")


if __name__ == "__main__":
    main()
