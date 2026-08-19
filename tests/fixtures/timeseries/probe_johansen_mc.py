# -*- coding: utf-8 -*-
# probe_johansen_mc.py - Monte Carlo 裁决: SM 三套 MHM96 表 vs urca OL1992 表
# 哪套表对应 coint_johansen 统计量的真实分布
#
# 背景: SM det0 计算已证 ≡ urca ecdet="none" 计算 (1e-10), 但两库临界值表不同:
#   SM tjcp0[0] = [2.9762, 4.1296, 6.9406]  (det=-1 表 q=1)
#   SM tjcp1[0] = [2.7055, 3.8415, 6.6349]  (det= 0 表 q=1, = χ²(1) 分位数)
#   SM tjcp2[0] = [2.7055, 3.8415, 6.6349]  (det= 1 表 q=1, 同 χ²(1))
#   urca none q=1 = [6.5, 8.18, 11.65]; const q=1 = [7.52, 9.24, 12.97]
# DGP: N=2 独立随机游走 (真 rank=0, 所有原假设为真), T=250, 4000 reps
# N=1+det>=0 不可用 (SM 形状 bug, pinv 收到 1-dim), 用 N=2 的 lr1[1] 同为 q=1

import numpy as np
from statsmodels.tsa.vector_ar.vecm import coint_johansen

NREP = 4000
T = 250


def mc(det_order):
    rng = np.random.RandomState(42)
    s0 = np.empty(NREP)  # r=0 → q=2
    s1 = np.empty(NREP)  # r=1 → q=1
    for i in range(NREP):
        y = np.cumsum(rng.standard_normal((T, 2)), axis=0)
        r = coint_johansen(y, det_order, 1)
        s0[i], s1[i] = r.lr1[0], r.lr1[1]
    return (np.quantile(s0, [0.90, 0.95, 0.99]),
            np.quantile(s1, [0.90, 0.95, 0.99]))


def main():
    print("MC: N=2 独立 RW, T=%d, %d reps, k_ar_diff=1" % (T, NREP))
    print("经验分位数 [90%% 95%% 99%%] (有限样本 T=250 vs 渐近表, 允许 ~±0.3):")
    tables = {
        -1: ("SM tjcp0 (det=-1 表)",
             np.array([2.9762, 4.1296, 6.9406]),
             np.array([10.4741, 12.3212, 16.3640])),
        0: ("SM tjcp1 (det= 0 表)",
            np.array([2.7055, 3.8415, 6.6349]),
            np.array([13.4294, 15.4943, 19.9349])),
        1: ("SM tjcp2 (det= 1 表)",
            np.array([2.7055, 3.8415, 6.6349]),
            np.array([16.1619, 18.3985, 23.1485])),
    }
    urca_none = ("urca none (OL1992)",
                 np.array([6.5, 8.18, 11.65]),
                 np.array([15.66, 17.95, 23.52]))
    for det in (-1, 0, 1):
        q2, q1 = mc(det)
        name, t_q1, t_q2 = tables[det]
        print("\nSM det=%+d  q=2 (r=0): [%7.4f %7.4f %7.4f]" %
              (det, q2[0], q2[1], q2[2]))
        print("  vs %-22s [%7.4f %7.4f %7.4f]  maxdiff=%.3f" %
              (name, t_q2[0], t_q2[1], t_q2[2],
               np.max(np.abs(q2 - t_q2))))
        print("  vs %-22s [%7.4f %7.4f %7.4f]  maxdiff=%.3f" %
              (urca_none[0], urca_none[2][0], urca_none[2][1],
               urca_none[2][2], np.max(np.abs(q2 - urca_none[2]))))
        print("SM det=%+d  q=1 (r=1): [%7.4f %7.4f %7.4f]" %
              (det, q1[0], q1[1], q1[2]))
        print("  vs %-22s [%7.4f %7.4f %7.4f]  maxdiff=%.3f" %
              (name, t_q1[0], t_q1[1], t_q1[2],
               np.max(np.abs(q1 - t_q1))))
        print("  vs %-22s [%7.4f %7.4f %7.4f]  maxdiff=%.3f" %
              (urca_none[0], urca_none[1][0], urca_none[1][1],
               urca_none[1][2], np.max(np.abs(q1 - urca_none[1]))))


if __name__ == "__main__":
    main()
