# -*- coding: utf-8 -*-
# gen_coint_fixture.py - 协整夹具 DGP (Phase 7C M3, 全 M3 脚本共用)
#
# DGP: rank=1 I(1) 系统, K=3, T=250, seed=42 (RandomState 跨版本固定, M2 惯例)
#   common_t = Σ e0 (随机游走共同因子)
#   y1 = common + 0.3·e1,  y2 = common + 0.3·e2   ⇒ y1−y2 平稳, β=(1,−1,0)
#   y3 = Σ e3 (独立随机游走)
#   → Johansen trace: r=0 拒绝 (≥1 协整), r=1 大概率不拒绝 (rank=1 非平凡)
#   → EG (y1,y2): 协整拒绝; (y1,y3): 无协整 (对照) — 双方向/双情形都非平凡
#
# 用法: python gen_coint_fixture.py  (生成 coint_smoke_data.csv)
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def gen_coint(T=250, seed=42):
    rng = np.random.RandomState(seed)
    e = rng.standard_normal((T, 4))
    common = np.cumsum(e[:, 0])
    y1 = common + 0.3 * e[:, 1]
    y2 = common + 0.3 * e[:, 2]
    y3 = np.cumsum(e[:, 3])
    return np.column_stack([y1, y2, y3])


def main():
    y = gen_coint()
    np.savetxt(os.path.join(HERE, "coint_smoke_data.csv"),
               y, fmt="%.17g", delimiter=",",
               header="y1,y2,y3", comments="")
    print("coint_smoke_data.csv written:", y.shape)
    # 自检: y1-y2 平稳性远强于 y1 本身 (ADF 粗查)
    from statsmodels.tsa.stattools import adfuller
    r_c = adfuller(y[:, 0] - y[:, 1], regression="n", autolag=None, maxlag=0)
    r_l = adfuller(y[:, 0], regression="n", autolag=None, maxlag=0)
    print("ADF(y1-y2) nc stat=%.4f  ADF(y1) nc stat=%.4f" % (r_c[0], r_l[0]))


if __name__ == "__main__":
    main()
