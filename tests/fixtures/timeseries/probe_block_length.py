"""Q11 块长对照探针: Politis-White (2004) 自动块长 vs ceil(sqrt(T))

参考实现: arch 8.0 `arch.bootstrap.optimal_block_length` (Patton-Politis-White
2009 修正, circular 分支)。C++ 侧 `politis_white_block_length()` 逐步镜像该
算法; 本脚本生成**整数序列** (LCG + 整数 AR(1)) 作为跨语言精确可复现的夹具,
因为整数运算无 libm 参与, 两边逐位一致。

用法: python tests/fixtures/timeseries/probe_block_length.py
输出: 各序列的 arch b_sb / b_cb (浮点) + ceil 后整数块长。
"""
import numpy as np
from arch.bootstrap import optimal_block_length

N = 500
S0 = 12345
MASK = (1 << 64) - 1
MUL = 6364136223846793005
INC = 1442695040888963407


def lcg_series(ar: bool) -> np.ndarray:
    """整数序列 (精确可复现): u_i ∈ {0,1,2,3}

    ar=True : x_i = 3*x_{i-1}//4 + u_i   (φ≈0.75 持久序列, x 恒非负)
    ar=False: x_i = u_i                  (i.i.d.)
    """
    s = S0
    x = np.zeros(N, dtype=np.int64)
    prev = 0
    for i in range(N):
        s = (s * MUL + INC) & MASK
        u = (s >> 33) % 4
        if ar:
            prev = (3 * prev) // 4 + int(u)
            x[i] = prev
        else:
            x[i] = u
    return x


def report(name: str, x: np.ndarray) -> None:
    df = optimal_block_length(x.astype(np.float64))
    b_sb = float(df["stationary"].iloc[0])
    b_cb = float(df["circular"].iloc[0])
    print(f"{name}: n={len(x)} stationary={b_sb:.12f} circular={b_cb:.12f} "
          f"ceil_circular={int(np.ceil(b_cb))}")
    print(f"  first10={list(x[:10])} last5={list(x[-5:])}")


def main() -> None:
    iid = lcg_series(ar=False)
    ar = lcg_series(ar=True)
    report("iid", iid)
    report("ar075", ar)
    print(f"sqrtT={int(np.ceil(np.sqrt(N)))}")


if __name__ == "__main__":
    main()
