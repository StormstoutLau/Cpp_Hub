# -*- coding: utf-8 -*-
# probe_johansen_pair.py - SM × urca 网格自动配对 (M3 双库映射裁决)
#
# 对每对 (sm_tag, urca_tag): 比较 eig 前 3 个 + lr1 + lr2 的最大绝对差
# 输出: 全部配对差异矩阵 + ≤1e-8 的匹配结论

import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(path):
    d = {}
    for line in open(path):
        parts = line.strip().split("\t")
        if len(parts) == 2 and parts[1].strip():
            d[parts[0]] = np.array([float(v) for v in parts[1].split()])
    return d


def main():
    sm = load(os.path.join(HERE, "johansen_sm_grid.txt"))
    ur = load(os.path.join(HERE, "johansen_urca_grid.txt"))

    sm_tags = sorted({k.rsplit("_", 1)[0] for k in sm
                      if k.endswith("_eig")})
    ur_tags = sorted({k.rsplit("_", 1)[0] for k in ur
                      if k.endswith("_eig")})

    print("=== 配对差异矩阵 (max |Δ| over eig[0:3] + lr1 + lr2) ===")
    print("%-16s %-28s %12s %12s %12s" %
          ("SM", "urca", "eig_diff", "lr1_diff", "lr2_diff"))
    matches = []
    for st in sm_tags:
        for ut in ur_tags:
            de = np.max(np.abs(sm[st + "_eig"][:3] - ur[ut + "_eig"][:3]))
            d1 = np.max(np.abs(sm[st + "_lr1"] - ur[ut + "_lr1"]))
            d2 = np.max(np.abs(sm[st + "_lr2"] - ur[ut + "_lr2"]))
            dmax = max(de, d1, d2)
            if dmax < 1e-6:
                matches.append((st, ut, de, d1, d2))
            elif dmax < 0.5:
                print("%-16s %-28s %12.3e %12.3e %12.3e" %
                      (st, ut, de, d1, d2))
    print()
    print("=== ≤1e-6 精确匹配对 ===")
    for st, ut, de, d1, d2 in matches:
        print("%-16s %-28s eig=%.3e lr1=%.3e lr2=%.3e" % (st, ut, de, d1, d2))
    print()
    print("=== 临界值表对照 (匹配对的 cvt, 行序均翻转后 r=0..N-1) ===")
    for st, ut, _, _, _ in matches:
        dcvt = np.max(np.abs(sm[st + "_cvt"] - ur[ut + "_cvt"]))
        dcvm = np.max(np.abs(sm[st + "_cvm"] - ur[ut + "_cvm"]))
        print("%-16s %-28s cvt_maxdiff=%.4f cvm_maxdiff=%.4f" %
              (st, ut, dcvt, dcvm))


if __name__ == "__main__":
    main()
