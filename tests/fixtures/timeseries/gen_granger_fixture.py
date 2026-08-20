# -*- coding: utf-8 -*-
# gen_granger_fixture.py - Granger 因果夹具 DGP (Phase 7C M1, verify_granger.py 共用)
#
# DGP: 三列平稳序列, T=250, seed=42 (RandomState 跨版本固定, M2/M3 惯例)
#   x (cause):  x_t = 0.6·x_{t−1} + u_t                  (AR(1) 驱动源)
#   y (effect): y_t = 0.3 + 0.5·y_{t−1} + 0.4·x_{t−1} + e_t
#               → x Granger-cause y (真因果, 一阶滞后), y ↛ x (反向无因果)
#   z (null):   z_t = 0.5·z_{t−1} + w_t                  (与 y/x 独立, 零假设对照)
#
# 用法: python gen_granger_fixture.py  (生成 granger_smoke_data.csv)
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def gen_granger(T=250, seed=42):
    rng = np.random.RandomState(seed)
    u = rng.standard_normal(T)      # x 新息
    e = rng.standard_normal(T)      # y 新息
    w = rng.standard_normal(T)      # z 新息
    x = np.zeros(T)
    y = np.zeros(T)
    z = np.zeros(T)
    for t in range(1, T):
        x[t] = 0.6 * x[t - 1] + u[t]
        y[t] = 0.3 + 0.5 * y[t - 1] + 0.4 * x[t - 1] + e[t]
        z[t] = 0.5 * z[t - 1] + w[t]
    return np.column_stack([x, y, z])


def main():
    d = gen_granger()
    np.savetxt(os.path.join(HERE, "granger_smoke_data.csv"),
               d, fmt="%.17g", delimiter=",",
               header="x,y,z", comments="")
    print("granger_smoke_data.csv written:", d.shape)


if __name__ == "__main__":
    main()
