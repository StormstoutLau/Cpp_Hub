# Johansen 双库 diff 报告 (冻结) — statsmodels vs urca

> **状态**: FROZEN (2026-08-19, M3 前置任务, spec §5.2/§6.2.1/决策 19)
> **证据脚本**: `tests/fixtures/timeseries/verify_johansen_diff.py` + `verify_johansen_diff.R`
> (网格转储) + `probe_johansen_pair.py` (自动配对) + `probe_johansen_mc.py`/`probe_johansen_mc2.py`
> (Monte Carlo 分布裁决) + `dump_ca_jo.R` (urca 源码导出)
> **数据**: `coint_smoke_data.csv` (T=250, N=3, y1/y2 协整 + y3 独立 RW)
> **版本**: statsmodels 0.14.4 / R 4.6.1 + urca (F:\R\win-library\4.6)

---

## 1. 参数映射 (上轮会话"双库差异显著"的根因裁决)

**statsmodels `k_ar_diff` = urca `K` − 1。**

- SM `coint_johansen(endog, det_order, k_ar_diff)`: k_ar_diff = 滞后差分阶数
- urca `ca.jo(x, K, ...)`: K = 水平 VAR 阶数 (含 K−1 个滞后差分)
- 上轮会话用 `coint_johansen(k_ar_diff=2)` 对照 `ca.jo(K=2)` — 模型阶数错位 (VAR(3) vs VAR(2)), 差异由此而来, **并非算法分歧**

## 2. 情形映射 (数值验证, probe_johansen_pair.py)

| SM det_order | SM 计算 (源码语义) | urca ecdet | 数值匹配 |
|---|---|---|---|
| −1 | 无任何确定性项 (H₂) | — | **无 urca 对应** (urca 三情形均含常数) |
| 0 | 水平去均值 + 辅回归隐含常数 → 短回归无约束常数 (H₁*) | "none" (Z1 显式含常数, `Z1 <- cbind(1, Z1)`) | ✓ **精确匹配** (eig 1.8e-13 / lr 7.4e-11) |
| 1 | 水平预去势 [1,t] + 去均值 (LeSage 管线) | — | **无 urca 对应** (urca "trend" 为受限趋势, 计算不同) |
| — | — | "const" (常数限入协整关系, H₁) | 无 SM 对应 |
| — | — | "trend" (趋势限入 CE + 无约束常数, H*) | 无 SM 对应 |

匹配对 (自动配对, ≤1e-6):

```
sm_det0_k1  ≡ ur_K2_none   eig=1.829e-13 lr1=7.367e-11 lr2=7.383e-11
sm_det0_k2  ≡ ur_K3_none   eig=3.116e-13 lr1=1.072e-10 lr2=1.074e-10
```

## 3. urca `spec` 参数恒等性 (数学恒等, 非数值巧合)

`spec="transitory"` 用 y_{t−1}, `"longrun"` 用 y_{t−K}, 但

```
y_{t−1} − y_{t−K} = Σ_{i=1}^{K−1} Δy_{t−i} ∈ span(Z1)
```

故两者对 Z1 回归的残差 RK 完全相同 → 特征值/统计量/特征向量张成空间恒等。
**urca spec 参数真空** (对任意 ecdet/K 成立)。

## 4. SM `lx = endog[1 : nobs−k]` 索引疑点裁决 (非 bug)

SM 水平回归用 y_{t−k} (而非 y_{t−1}), 但差分回归元含 Δy_{t−1}..Δy_{t−k}, 而

```
y_{t−1} − y_{t−k} = Σ_{i=1}^{k−1} Δy_{t−i} ∈ span(z)
```

→ rkt 与"标准 transitory 对齐"的残差恒等, 对任意 k_ar_diff 成立。
**数值证据**: sm_det0_k2 ≡ ur_K3_none (transitory) 至 1e-10。SM 实现数学等价, 无需复刻异常。

## 5. 临界值表: 两库不同源且 SM det∈{0,1} 内部不一致 (MC 裁决)

- SM: MHM96 量化值 (LeSage johansen.m 转录, `coint_tables.py`, n≤12, 三套 p∈{−1,0,1} × trace/maxeig)
- urca: OL1992 表 (源码常量转录, n≤11, 三套 ecdet × trace/maxeig)
- **名义对应情形 (det0 ↔ none) 的两套表数值显著不同**, 且 MC (N=2, T=250, 4000 reps) 裁决:

| 统计量来源 | MC 经验分位数 [90/95/99%] | 匹配的表 | 不匹配的表 |
|---|---|---|---|
| det=−1, q=2 | [10.56, 12.29, 15.81] | tjcp0 [10.47, 12.32, 16.36] ✓ | urca none ✗ |
| det=−1, q=1 (真 rank=1 DGP) | [2.93, 4.07, 7.10] | tjcp0 [2.98, 4.13, 6.94] ✓ | — |
| det=0, q=2 | [15.91, 18.39, 22.78] | **urca none [15.66, 17.95, 23.52] ✓** | tjcp1 [13.43, 15.49, 19.93] ✗ (差 2.9) |
| det=0, q=1 (真 rank=1 DGP) | [6.59, 8.18, 11.70] | **urca none [6.50, 8.18, 11.65] ✓ (三位吻合)** | tjcp1 = χ²(1) [2.71, 3.84, 6.63] ✗ |
| det=1, q=2 | [21.61, 24.15, 29.60] | ≈ urca trend [22.76, 25.32, 30.45] (差 ~1.2, 有限样本可解释) | tjcp2 [16.16, 18.40, 23.15] ✗ (差 6.5) |
| det=1, q=1 (真 rank=1 DGP) | [9.80, 11.61, 15.88] | ≈ urca trend [10.49, 12.25, 16.26] ✓ | tjcp2 = χ²(1) ✗ |

**结论**: tjcp1/tjcp2 在 q=1 处均为 χ²(1) 分位数 (受限确定性情形的特征), 而 SM
det0/det1 的**计算**遵循无约束情形分布 (与 urca none/trend 表吻合)。即 LeSage 管线的
p 索引把受限情形的表挂到了无约束情形的计算上 — **statsmodels det_order∈{0,1} 的
cvt/cvm 与其自身统计量分布不一致** (det_order=−1 自洽)。urca 三情形计算↔表自洽。

## 6. 冻结决策 (C++ 实现口径)

1. **统计量主对照 = statsmodels** (算法复刻 LeSage/SM 管线: 去势→lagmat→残差化→RRR 特征值):
   - `coint_johansen` eig/lr1/lr2/evec vs SM **1e-10**, 全部 det_order ∈ {−1,0,1} × k_ar_diff ≥ 1
   - urca 交叉验证 (仅可映射情形): det_order=0 ↔ ecdet="none", k_ar_diff = K−1, 容差 1e-8 (实测 1e-10)
   - `select_coint_rank` 复刻 SM 逐级检验 (lr[r] < cv[r, signif_idx] 即接受 r)
2. **临界值双表并存**:
   - `JohansenResult.cvt/cvm` 默认 = MHM96 表 (SM 转录, det_order 键控) — 与 SM
     `cvt/cvm` 及 `select_coint_rank` 行为逐位一致, `cv_source = "MHM96"` 回显。
     原因: OL1992 (urca 转录) 按 ecdet 键控, 仅 none↔det0 可映射, det_order=−1/1
     **无 OL1992 表可用**, 无法作为 cvt/cvm 的统一默认源 (spec 决策 19 "OL1992 主录"
     在此不可行, 由本报告裁决改录, 回溯记录于 §7)
   - `osterwald_lenum_cv.hpp` 提供 OL1992 独立查表 (urca 源码常量转录, ecdet 键控,
     constexpr + static_assert, 精确相等) — 双对照 + 用户可用正确表源做推断
   - 头文件文档显式警告: SM det_order∈{0,1} 的表与其统计量分布不一致 (§5); 做实际
     协整 rank 推断建议用 OL1992 查表 (det0→none 表; det1→trend 表近似) 或 det_order=−1
3. **evec 归一**: 复刻 SM (β'·S₁₁·β = I + 首非零元符号约定); 与 urca 对照需首行归一后比较
   (urca V 为 V[1,j]=1 归一); C++ 测试与 SM 对照按列符号对齐后逐元素 (特征向量列符号
   依赖 LAPACK 特征分解次序, 逐元素裸比不稳健)

## 7. 与 spec 的偏离记录 (回溯, 非静默)

| spec 原文 | 冻结修改 | 依据 |
|---|---|---|
| §5.2 "cvt/cvm: OL1992 主录" | 默认表源改 MHM96, OL1992 降为独立查表 API + 双对照 | OL1992 无 det_order=−1/1 对应表 (§5); MHM96 为唯一全覆盖表源且保 `select_coint_rank` ≡ SM |
| §5.2 "Johansen vs urca ca.jo 1e-8" (暗示全情形) | 限定为 det_order=0 ↔ ecdet="none" 映射情形 | 双库确定性情形集不合 (§2), 仅此情形可映射 |
| §9.3 CI7 "与 statsmodels 对照必须用 transitory" | spec 恒等性裁决 (§3), transitory/longrun 无差别 | y_{t−1}−y_{t−K} ∈ span(Z1) 恒等式 |
| 文件清单 `ol1992_cv.inc` (仅 Johansen 一表) | 追加 MHM96 表转录 (`ol1992_cv.inc` 同目录伴生) | 决策 19 "双对照" 需要 SM 表源落地 |

---

## 附录 A: 复现命令

```
python tests/fixtures/timeseries/verify_johansen_diff.py     # SM 网格 → johansen_sm_grid.txt
Rscript tests/fixtures/timeseries/verify_johansen_diff.R     # urca 网格 → johansen_urca_grid.txt
python tests/fixtures/timeseries/probe_johansen_pair.py      # 自动配对 → §2 匹配对
python tests/fixtures/timeseries/probe_johansen_mc.py        # MC 裁决 q=2 行
python tests/fixtures/timeseries/probe_johansen_mc2.py       # MC 裁决 q=1 行 (真 rank=1 DGP)
Rscript tests/fixtures/timeseries/dump_ca_jo.R               # urca ca.jo 源码导出 (OL1992 常量转录源)
```

## 附录 B: SM coint_johansen 算法语义 (一手源码, vecm.py L603-737)

```
detrend(y, order): order=-1 原样; 0 去均值; 1 对 [1,t] (vander(linspace(-1,1,T),2)) OLS 残差
f = 0 若 det_order > -1 否则 -1     ← 二轮去势阶数: det>=0 时仅去均值
endog ← detrend(endog, det_order)
dx = diff(endog);  z = lagmat(dx, k)[k:];  dx = dx[k:]
z, dx ← detrend(·, f);  r0t = resid(dx, z)                    # Δy_t ~ Δy_{t-1..t-k}
lx = endog[1 : T−k];  lx ← detrend(lx, f);  rkt = resid(lx, z) # y_{t-k} ~ Δy lags (≡ y_{t-1}, §4)
skk = rkt'rkt/t;  sk0 = rkt'r0t/t;  s00 = r0t'r0t/t
au, du = eig(inv(skk) · sk0·inv(s00)·sk0')
dt = du · inv(chol(du' skk du))     # β'S₁₁β = I 归一
降序排序; 全矩阵乘 sign(首非零元)    # 符号约定 (列符号仍任意, 依赖 LAPACK)
lr1[r] = −t·Σ_{i>r} ln(1−λᵢ);  lr2[r] = −t·ln(1−λ_{r+1});  t = rkt 行数 = T−1−k
cvt[r,:] = c_sjt(N−r, det_order);  cvm[r,:] = c_sja(N−r, det_order)   # MHM96
```
