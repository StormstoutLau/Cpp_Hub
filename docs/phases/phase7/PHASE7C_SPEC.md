# Phase 7C 执行规格书 - 多变量时序与混频模块 (v1.7 M0-M4)

> **版本归属**: **v1.7** (Phase 7C)
> **目标**: 实现 v1.7 全部时序扩展 — 回填三项 (M0) + ARIMA/Granger (M1) + VAR/IRF/FEVD/DY (M2) + 协整三件套 (M3) + MIDAS (M4), 补齐 Research OS 因子诊断链的多变量与混频基础设施
> **覆盖范围**:
>   - M0: Ng-Perron M 检验族 + Zivot-Andrews 结构断点单位根 + GARCH(1,1)-M 三变体
>   - M1: ARIMA(p,d,q) CSS/CSS-ML/innovations MLE + 线性 Granger (标准 F/χ²/LR + Toda-Yamamoto + HAC-Wald)
>   - M2: VAR(p) 逐方程 OLS + IC 五式 + IRF + FEVD 双轨 (Cholesky/generalized) + Diebold-Yilmaz 溢出指数
>   - M3: Engle-Granger + Johansen + VECM + Phillips-Ouliaris
>   - M4: MIDAS-DL/MIDAS-AR/U-MIDAS + nealmon/nbeta/almonp/polystep/harstep 权重族 + 集中化 NLS
> **前置**:
>   - Phase 7B (v1.6) 已验收: 2207 测试三平台 + CI 全绿 (commit 0e62ce7)
>   - 可复用: GARCH QMLE 框架 + SLSQP (ADR-018) + ADF/DF-GLS 引擎 + MacKinnon 2010 表管线 + HAC/OLS/QR + block bootstrap + multiple_test_correction
> **里程碑顺序**: M0 (框架已热, 先行) → M1 ∥ M4 (并行) → M2 (依赖 M1) → M3 (依赖 M2)
>
> **Scope 声明** (ADR-019 冻结, 26+3 项决策):
> - **ARIMA**: 仅 CSS + CSS-ML + innovations 精确 MLE (限无缺失/无季节); Kalman 留 v1.8; 不做 SARIMA / wild bootstrap
> - **VAR**: 仅逐方程 OLS; FEVD 双轨; GFEVD 自实现 (statsmodels 无); 不做 SVAR/BVAR/TVP-VAR
> - **协整**: EG + Johansen (3 det 情形 API) + VECM (5 情形) + PO (Pu/Pz); ARDL/PSS 出 scope
> - **MIDAS**: DL/AR/U-MIDAS + 5 权重函数; 不做 midas_nlpr/sp/qr/imidas_r/amweights (v1.8+)
> - **回填**: NP 四统计量 + MAIC/MBIC/seq-t; ZA 三模型 A/B/C; GM 三变体 (h/√h/log h)
> - **推迟 v1.8**: 多元 GARCH (CCC/DCC) / Kalman / 长记忆族 (APARCH/FIGARCH/IGARCH/ARFIMA/HYGARCH) / KDE/KNN/TE
> - **对照库版本冻结**: Python arch 8.0.0 (GM 主锚) / statsmodels 0.14.6 / R midasr 0.9 (MIDAS 唯一主基准) / R urca (Johansen/ZA/PO 交叉) / R rugarch (GM archpow) / R vars (VAR 交叉) / R Spillover (GFEVD/DY 主基准)
>
> **关联**:
> - [ADR-019: v1.7 多变量时序与混频实施边界 (26+3 项)](../../decisions/ADR-019_V17_TIMESERIES_BOUNDARY.md) (Accepted 2026-08-17)
> - [PHASE7C_RESEARCH.md](../../research/PHASE7C_RESEARCH.md) v1.1 (调研报告, 3 agent 126 条全量审计, 8 处实质错误已修正)
> - [ADR019_REVIEW_PILOT.md](../../research/ADR019_REVIEW_PILOT.md) v1.1 (R1-R4 门禁首轮 pilot, 5/5 双盲 TRUE)
> - [ADR-016](../../decisions/ADR-016_FINANCIAL_TIMESERIES_BOUNDARY.md) / [ADR-017](../../decisions/ADR-017_TIMESERIES_NAMESPACE.md) / [ADR-018](../../decisions/ADR-018_SLSQP_BOUNDARY.md)
> - [PHASE7B_FINANCIAL_TS_SPEC.md](./PHASE7B_FINANCIAL_TS_SPEC.md) (v1.6 前置)

---

## 0. 门禁记录 (新工作流: R 清零 → spec 冻结)

| 门禁 | 状态 | 记录 |
|------|------|------|
| R1 断言统计表 | ✅ | 调研报告 v1.1 + 复核 pilot §0 |
| R2 A 类证据可核验 | ✅ (1 处修正) | 复核 pilot §3: B1 引文转述失准当场修正 |
| R3 B 类登记+审计 | ✅ | `assertion_audit.py` 探针 B1-B4 SURVIVED + 双盲 5/5 TRUE |
| R4 阻断性清零 | ✅ | 无 FALSIFIED/CONFLICT; 阻断性断言全部双源; H1/H2 入本 spec §13 [待定] |

**本 spec 引用的全部阻断性断言均经双源验证** (调研报告 v1.1 审计 + 复核 pilot 双盲), 实施期发现冲突时回溯 ADR-019 修订流程, 不得静默改公式。

---

## 1. 交付物清单

### 1.1 新增文件结构 (按 ADR-017 命名空间组织)

```
# M0 回填: NP/ZA 落 unit_root, GM 落 garch (复用 Phase 7B 框架)
include/cpphub/timeseries/unit_root/
├── ng_perron_test.hpp          # 新增 M0: MZα/MZt/MSB/MPT + MAIC/MBIC/seq-t (复用 GLS 去势)
├── zivot_andrews_test.hpp      # 新增 M0: 三模型 A/B/C + 断点搜索 + trim
└── np_tables.hpp               # 新增 M0: NP 2001 Table 1 临界值 (constexpr + static_assert)
include/cpphub/timeseries/garch/
└── garch_m_model.hpp           # 新增 M0: GARCH(1,1)-M 三变体 g∈{h,√h,log h}

# M1: ARIMA + Granger (纯标量, 无 Eigen3, C1)
include/cpphub/timeseries/arima/
├── arima_model.hpp             # 新增 M1: CSS + CSS-ML (SLSQP/BFGS 精化)
├── innovations_mle.hpp         # 新增 M1: innovations 精确 MLE (B&D 2016 ITSF §5.2)
├── hannan_rissanen.hpp         # 新增 M1: HR 起始值 (多起始点策略之一)
└── granger_test.hpp            # 新增 M1: 标准 F/χ²/LR + TY 增广 Wald + HAC-Wald

# M2: VAR/DY (需 Eigen3 → cpphub_timeseries_mat, 见 §8.2)
include/cpphub/timeseries/var/
├── var_model.hpp               # 新增 M2: 逐方程 OLS + is_stable 双输出 + P 注入/变量重排
├── var_select.hpp              # 新增 M2: IC 五式 + select_order 同样本 offset
├── irf.hpp                     # 新增 M2: 正交化 IRF + block bootstrap 置信带
├── fevd.hpp                    # 新增 M2: Cholesky FEVD (行和=1) + GFEVD (行归一化, 双系数框架)
└── dy_spillover.hpp            # 新增 M2: DY 溢出指数 (TO/FROM/NET/TCI + 滚动)

# M3: 协整 (需 Eigen3)
include/cpphub/timeseries/cointegration/
├── engle_granger.hpp           # 新增 M3: EG 两步法 + MacKinnon 1994 响应面
├── johansen_test.hpp           # 新增 M3: 迹/最大特征值 + 3 det 情形
├── phillips_ouliaris.hpp       # 新增 M3: Pu/Pz 双实现
├── vecm_model.hpp              # 新增 M3: VECM 5 情形 + β 双归一 + ECT t 检验
├── osterwald_lenum_cv.hpp      # 新增 M3: OL1992 表 (urca 转录 constexpr + static_assert)
├── mackinnon_coint_cv.hpp      # 新增 M3: MacKinnon 1994 协整响应面 (与 2010 ADF 表分文件)
└── ericsson_mackinnon_cv.hpp   # 新增 M3: EM2002 ECT t 临界值查表

# M4: MIDAS (纯标量, 无 Eigen3)
include/cpphub/timeseries/midas/
├── mixed_freq_data.hpp         # 新增 M4: mls 对齐 (lag0=期末) + MixedFreqData
├── midas_weights.hpp           # 新增 M4: nealmon/nbeta/almonp/polystep/harstep + log-sum-exp
├── midas_model.hpp             # 新增 M4: MIDAS-DL/AR(+AR*)/U-MIDAS + 集中化 NLS
└── midas_diagnostics.hpp       # 新增 M4: 残差诊断 (复用 Phase 7A) + hAh 权重检验

# 数据载体 (C4)
include/cpphub/timeseries/var/multivariate_data.hpp    # 新增 M2: MultivariateTSData (K 列等长)
include/cpphub/timeseries/midas/mixed_freq_data.hpp    # 新增 M4: 低频 y + 高频 x + m

# verify 脚本 (不入版本控制, 沿用 C6 模式: probe_*.py 源码核查 + verify_*.py/.R 基准 + .inc 黄金值)
tests/fixtures/timeseries/
├── verify_np_stata.py          # Stata dfgls r(results) 逐 k MAIC → CSV (唯一确定性免费对照)
├── verify_np_semi.py           # EViews/Julia 抄录值半基准 (1e-6~1e-8, Julia 降权)
├── verify_za.py                # statsmodels zivot_andrews (Baum 近似模式)
├── verify_za.R                 # urca ur.za (固定 lag, trim 放开)
├── verify_gm.py                # arch 8.0 ARCHInMean (form='vol'/'var'/'log')
├── verify_gm.R                 # rugarch archpow=1/2 (含 fix() 三步法)
├── verify_arima.R              # R stats::arima CSS/CSS-ML (method 配对)
├── verify_arima.py             # statsmodels ARIMA innovations_mle (无缺失非季节)
├── verify_granger.py           # grangercausalitytests 4 统计量
├── verify_var.py               # statsmodels VAR (系数/IC/IRF/FEVD)
├── verify_var.R                # vars::VAR + VARselect 交叉
├── verify_gfevd.R              # R Spillover g.fevd/G.spillover (GFEVD 主基准)
├── verify_eg.py                # statsmodels coint
├── verify_johansen.py          # statsmodels coint_johansen
├── verify_johansen_diff.R      # urca ca.jo vs statsmodels 双库 diff (前置冻结主对照)
├── verify_vecm.py              # statsmodels VECM (β 张成空间对照)
├── verify_po.R                 # urca ca.po
├── verify_midas.R              # midasr 0.9 夹具生成 (收紧 optim control + set.seed)
├── verify_midas_u.R            # midas_u 纯 OLS 锚 (1e-10)
└── critical_values/
    ├── np_table1.inc           # NP 2001 Table 1 (4 统计量 × 常数/趋势 × 1/5/10%)
    ├── za1992_cv.inc           # ZA 1992 论文表 (主)
    ├── za_mc_cv.inc            # statsmodels/arch 共用 MC 表 (可选对照)
    ├── ol1992_cv.inc           # OL1992 (urca 转录, constexpr + static_assert)
    ├── mackinnon1994_coint.inc # MacKinnon 1994 协整响应面 (N≥2)
    └── em2002_ect_cv.inc       # Ericsson-MacKinnon 2002 ECT t 临界值
```

### 1.2 新增测试套件

| 测试套件 | 用例数 | 覆盖模块 | 里程碑 |
|----------|--------|----------|--------|
| `test_ng_perron` | 18 | 四统计量 + MAIC/MBIC/seq-t + 恒等式自检 + 临界值 | M0 |
| `test_zivot_andrews` | 15 | 三模型 + DU/DT 构造 + trim + 断点搜索 | M0 |
| `test_garch_m_model` | 16 | 三变体 + λ sandwich + 均值方差耦合 | M0 |
| `test_arima_model` | 24 | CSS/CSS-ML + n.cond + drift + 多起始 + AIC 口径 | M1 |
| `test_innovations_mle` | 12 | innovations 精确 MLE + 无缺失前提 + GLS trend | M1 |
| `test_granger_causality` | 16 | 4 统计量 + 方向 + TY df + HAC-Wald | M1 |
| `test_var_model` | 20 | 逐方程 OLS + IC 五式 + 同样本 + 稳定性 + P 注入 | M2 |
| `test_var_irf_fevd` | 18 | IRF 方向 + FEVD 行和 + GFEVD 双框架 | M2 |
| `test_dy_spillover` | 12 | TO/FROM/NET/TCI + H 可配 + 滚动窗口 | M2 |
| `test_engle_granger` | 14 | 两步法 + 1994 响应面 + p/cv 不同源声明 | M3 |
| `test_johansen_test` | 18 | 迹/最大特征值 + 3 情形 + λ 降序 + OL1992/双表 diff | M3 |
| `test_vecm_model` | 16 | 5 情形 + β 双归一 (投影空间对照) + ECT t 查表 | M3 |
| `test_phillips_ouliaris` | 10 | Pu 方向依赖 / Pz 方向无关 | M3 |
| `test_midas_weights` | 16 | 5 权重族 + 索引起点 + log-sum-exp + δ 尺度 | M4 |
| `test_midas_model` | 18 | DL/AR/U-MIDAS + 集中化 NLS + 多起点 | M4 |
| `test_integration_phase7c` | 8 | 端到端 6 场景 (§7) | 全部 |

**新增测试总数**: ~263 (累计 ~2470)

### 1.3 必须达到的数值基准 (容差分层, C7)

| 基准 | 容差 | 验收方式 |
|------|------|----------|
| NP 逐 k MAIC/σ̂² vs Stata dfgls r(results) | **1e-10** | `test_ng_perron.cpp` |
| NP 四统计量 vs EViews 抄录 (Julia 降权旁证) | 1e-6~1e-8 | `verify_np_semi.py` |
| NP 恒等式 MZt ≡ MZα×MSB (内部自检) | 1e-12 | `test_ng_perron.cpp` |
| NP 渐近 5% 临界值 vs NP 2001 Table 1 | **精确相等** | `np_table1.inc` |
| ZA 统计量/断点 vs statsmodels (Baum 模式, lag 冻结) | 1e-10 | `test_zivot_andrews.cpp` |
| ZA vs urca ur.za (固定 lag, trim 放开) | 1e-8 | `verify_za.R` |
| ZA 临界值 vs ZA1992 论文表 (主) / MC 表 (对照) | 精确相等 | `za1992_cv.inc` |
| GM vs arch 8.0 ARCHInMean (form='vol'/'var'/'log') | 1e-8~1e-10 | `test_garch_m_model.cpp` |
| GM vs rugarch archpow=1/2 (solver 敏感) | 1e-4 + fix() 三步法 | `verify_gm.R` |
| ARIMA CSS vs R stats::arima method="CSS" (method 配对) | 1e-10 | `verify_arima.R` |
| ARIMA CSS-ML vs R method="CSS-ML" | 1e-8 | `verify_arima.R` |
| ARIMA innovations vs statsmodels innovations_mle (无缺失非季节) | 1e-10 | `verify_arima.py` |
| Granger 4 统计量 vs statsmodels grangercausalitytests | 1e-10 | `verify_granger.py` |
| TY Wald vs 文献数值例 (Zapata-Gil 1999 / 参考实现) | 1e-8 | `test_granger_causality.cpp` |
| VAR 系数/IC/IRF/FEVD vs statsmodels var_model | 1e-10 | `verify_var.py` |
| VAR IC vs R vars::VARselect (交叉) | 1e-8 | `verify_var.R` |
| GFEVD/DY vs R Spillover g.fevd/G.spillover | 1e-8 | `verify_gfevd.R` |
| EG 统计量/临界值 vs statsmodels coint | 1e-10 | `verify_eg.py` |
| Johansen eig/lr1/lr2 vs statsmodels coint_johansen (transitory) | 1e-10 | `verify_johansen.py` |
| Johansen vs urca ca.jo (OL1992 表 diff 报告, 冻结主对照) | 1e-8 | `verify_johansen_diff.R` |
| OL1992 表 vs urca 源码常量 | **精确相等** | `ol1992_cv.inc` + static_assert |
| VECM 系数 vs statsmodels VECM; β 比较投影矩阵 P=β(β'β)⁻¹β' | 1e-10 | `verify_vecm.py` |
| PO Pu/Pz vs urca ca.po | 1e-8 | `verify_po.R` |
| U-MIDAS vs midasr midas_u (纯 OLS) | **1e-10** | `verify_midas_u.R` |
| MIDAS NLS vs midasr midas_r 夹具 (收紧 reltol=1e-12, maxit=10000, set.seed) | 1e-6~1e-8 | `verify_midas.R` |
| nealmon/nbeta 权重 vs midasr 逐点 | 1e-12 | `test_midas_weights.cpp` |
| log-sum-exp 与裸公式差 (非溢出区间) | <1e-14 | `test_midas_weights.cpp` |

---

## 2. M0: 回填三项 (先行, 框架已热)

> **定位**: Phase 7B 推迟项回填, 复用已验收框架 (DF-GLS 去势 / ADF 引擎 / GARCH QMLE)
> **命名空间**: NP/ZA 落 `cpphub::v1::timeseries::unit_root`, GM 落 `...::garch` (C2)

### 2.1 Ng-Perron M 检验族 (`ng_perron_test.hpp`)

**教材锚点**: Ng-Perron 2001 *Econometrica* 69(6):1519-1554 (DOI 10.1111/1468-0262.00256; ⚠️ 非 J.Econometrics, NP1)
**对照生态** (调研域 E 关键发现): 无任何成熟开源库输出 M 族 (statsmodels/urca 均无; Stata dfgls 仅输出 DF-GLS τ + 逐 k MAIC, 不输出 M 统计量, NP6) → 基准 = **原文公式钉死 + 文献 Table + Stata 逐 k MAIC (确定性) + EViews/Julia 半基准**

**接口签名**:

```cpp
namespace cpphub::v1::timeseries::unit_root {

struct NgPerronResult {
    // 四统计量 (NP 2001 §4)
    Real mz_alpha;      ///< MZα = (T⁻¹y_T² − s²) / (2·T⁻²Σỹ²)
    Real mz_t;          ///< MZt = MZα × MSB  (恒等式, 自检用)
    Real msb;           ///< MSB = √(T⁻²Σỹ² / s²)
    Real mpt;           ///< MPT, 分情形 (NP4): 常数 c̄=−7 / 趋势 c̄=−13.5
    // 逐项临界值 (1/5/10%) 与结论
    Real cv_1pct[4];    ///< {MZα, MZt, MSB, MPT}
    Real cv_5pct[4];
    Real cv_10pct[4];
    bool reject_5pct[4];
    // 滞后选择
    Size selected_lag;          ///< MAIC 最优 k
    std::vector<Real> maic;     ///< 逐 k MAIC(k) 轨迹 (与 Stata dfgls r(results) 对照)
    std::vector<Real> sigma2_k; ///< 逐 k σ̂²(k) (Stata 对照点, 1e-10)
    std::string trend_spec;     ///< "c" (c̄=−7) / "ct" (c̄=−13.5)
    std::string summary;
};

NgPerronResult ng_perron_test(const std::vector<Real>& data,
                               const std::string& trend_spec = "ct",
                               Size max_lag = 0);   // 0 => Schwert 上限

}  // namespace
```

**算法步骤**:

```
Step 1: GLS 去势 (复用 df_gls_test.hpp 的 ERS 变换)
  1.1 c̄: "c" → −7.0; "ct" → −13.5; ρ̄ = 1 + c̄/T
  1.2 ỹ = GLS 去势后序列 (与 DF-GLS 完全同一变换, 不重复实现)

Step 2: ADF 型辅助回归 (对 Δỹ, 逐 k)
  2.1 Δỹ_t = β₀·ỹ_{t−1} + Σ_{i=1}^{k} δᵢ·Δỹ_{t−i} + ε_t   (无常数, 已去势)
  2.2 σ̂²(k) = SSR(k)/(T − k_max − 1)   (固定样本 T − k_max 口径, τ_T 求和同口径)
  2.3 τ_T(k) = β̂₀² · Σ_{t=k_max+1}^{T} ỹ²_{t−1} / σ̂²(k)
      ⚠️ NP2 已 4 源裁决 (BC wp369 + AU ng_perron00 双工作稿 eq.(12) + Stata 手册 + Zivot 讲义):
      λ̂−λ̃ 差形式为 v1.0 幻觉, 零命中
  2.4 MAIC(k) = ln σ̂²(k) + C_T·(τ_T(k) + k)/(T − k_max),  C_T=2
      MBIC(k): C_T = ln(T − k_max);  seq-t: 逐 k t 统计量终止规则 (三准则全存)
  2.5 k* = argmin MAIC(k), k = 0..max_lag

Step 3: 长期方差 s² 与 AR 谱校正 (NP3)
  3.1 AR 谱密度对 **Δỹ** (去势序列的差分) 拟合, 非 ỹ 水平值:
      s²_AR = σ̂²(k*) / (1 − Σφ̂)²   (对水平拟合则分母→0 爆炸, NP3)
      φ̂ 为 Δỹ 的 AR(k*) 系数 (Step 2 回归的 δᵢ)

Step 4: 四统计量 (k = k*, Σỹ² = Σ ỹ²_{t−1} 固定样本口径)
  4.1 MZα = (T⁻¹·ỹ_T² − s²) / (2·T⁻²·Σỹ²)
  4.2 MSB = √(T⁻²·Σỹ² / s²)
  4.3 MZt = MZα × MSB   (恒等式, 断言 1e-12 自检)
  4.4 MPT 分情形 (NP4, v1.1 审计修正):
      常数 (c̄=−7):    MPT = (c̄²·T⁻²Σỹ² − c̄·T⁻¹·ỹ_T²) / s²
                       = (49·T⁻²Σỹ² + 7·T⁻¹·ỹ_T²) / s²
      趋势 (c̄=−13.5): MPT = (c̄²·T⁻²Σỹ² + (1−c̄)·T⁻¹·ỹ_T²) / s²
                       末项系数 **+14.5** (非 −c̄=13.5, H1 Julia 疑偏即此处)

Step 5: 临界值 (np_table1.inc, NP 2001 Table 1 转录)
  渐近 5% 锚: 常数 {−8.10, −1.98, 0.233, 3.17}; 趋势 {−17.30, −2.91, 0.168, 5.48}
  拒绝方向: MZα/MZt 越负越拒绝; MSB/MPT 越小越拒绝 (NP5)
```

**幻觉点映射**: NP1 (Econometrica) / NP2 (τ_T 4 源裁决) / NP3 (AR 谱对差分) / NP4 (MPT 分情形+14.5) / NP5 (方向+来源) / NP6 (Stata 非基准, 仅 MAIC 对照)

**测试矩阵要点**: MZt≡MZα×MSB 恒等式 1e-12; 逐 k MAIC/σ̂² vs Stata dfgls 1e-10; MPT 趋势末项 +14.5 vs (−c̄) 区分断言; 常数/趋势临界值精确相等; 模拟 AR(1) φ=0.95 → 四统计量正确拒绝方向。

---

### 2.2 Zivot-Andrews 结构断点单位根 (`zivot_andrews_test.hpp`)

**教材锚点**: Zivot-Andrews 1992 *JBES* 10(3):251-270 (⚠️ statsmodels 文档误写刊名 "J. Business & Economic Studies", 勿照抄)
**对照库**: statsmodels `zivot_andrews` (0.11.0 引入, PR #6014; Baum 2004/2015 近似 — 滞后一次性预选) / R urca `ur.za` (固定 lag, 无 trim) / arch (与 statsmodels 共用 MC 表)

**接口签名**:

```cpp
namespace cpphub::v1::timeseries::unit_root {

enum class ZAModel { A, B, C };   ///< A: 崩溃均值 / B: 断裂趋势 / C: 两者 (ZA 1992 三模型)

struct ZAResult {
    Real statistic;           ///< min over Tb of t(α̂)  (ZA5: 最负, 非 max)
    Size break_index;         ///< 最优断点位置
    Real critical_1pct, critical_5pct, critical_10pct;  ///< ZA1992 论文表 (主)
    Real p_value_baum;        ///< Baum 近似 p 值 (对照 statsmodels 用)
    ZAModel model;
    Size n_lags;
    Real trim;                ///< 默认 0.15 (上限 0.333)
    std::vector<Real> t_stats_path;  ///< 逐断点 t(α̂) 轨迹 (诊断用)
    bool reject_null;         ///< H0: 带单断点单位根
};

ZAResult zivot_andrews_test(const std::vector<Real>& data,
                            ZAModel model = ZAModel::C,
                            Size fixed_lag = 0,        // 0 => Baum 式预选 (对照 statsmodels)
                            Real trim = 0.15);

}  // namespace
```

**算法步骤**:

```
Step 1: 断点网格: Tb ∈ [trim·T, (1−trim)·T] (默认 0.15; urca 复现须 trim 放开, ZA2)
Step 2: 每个候选 Tb, 估计增广 DF 回归 (滞后策略按模式):
  Model A: Δy_t = c + α·y_{t−1} + β·DU_t + Σ δᵢ Δy_{t−i} + ε_t
  Model B: Δy_t = c + β·t + γ·DT_t + α·y_{t−1} + Σ δᵢ Δy_{t−i} + ε_t
  Model C: Δy_t = c + β·t + γ·DT_t + θ·DU_t + α·y_{t−1} + Σ δᵢ Δy_{t−i} + ε_t
  ⚠️ ZA3: DU_t = 1(t > Tb); DT_t = (t − Tb)·1(t > Tb)  (断点后重新计时, 非 DU·t 全局)
Step 3: 滞后策略双模式:
  - 主模式: 固定 lag (与 urca 对齐, 参数 fixed_lag)
  - 对照模式: Baum 式一次性预选 (固定 T 全样本选 lag 后搜索断点, 对齐 statsmodels)
Step 4: 统计量 = min_{Tb} t(α̂) (ZA5); 记录 argmin 断点
Step 5: 临界值: 主用 ZA1992 论文表 (za1992_cv.inc, 对齐 urca); MC 表 (za_mc_cv.inc,
  c 1% = −5.27644, ⚠️ 非 −5.83 (0.1% 分位值); 5% −4.81067; 10% −4.56618; t 1% −5.03421;
  ct 1% −5.57556) 作为可选对照, API 注明采用哪套
```

**幻觉点映射**: ZA1 (三库算法差异 → lag 冻结) / ZA2 (trim 参数化) / ZA3 (DT 构造) / ZA4 (双临界值表, MC c 1%=−5.27644) / ZA5 (取最负)

---

### 2.3 GARCH(1,1)-M (`garch_m_model.hpp`)

**教材锚点**: Engle-Lilien-Robins 1987 Econometrica 55(2):391-407 / Glosten-Jagannathan-Runkle 1993
**对照库**: arch 8.0.0 `ARCHInMean` (form='vol'/'var'/'log'/float; ⚠️ 8.0.0 release notes 漏记, 以 v8.0.0 tag 源码为准) + rugarch `archm=TRUE, archpow=1|2`

**接口签名**:

```cpp
namespace cpphub::v1::timeseries::garch {

enum class GarchMForm { Variance, Volatility, LogVariance };
///< rugarch 对偶: archpow=2↔Variance (g=h), archpow=1↔Volatility (g=√h);
/// arch 8.0 对偶: form='var'↔Variance, 'vol'↔Volatility, 'log'↔LogVariance
/// (GM1: archpow 1=σ_t, 2=σ_t²; GM3: log 变体有 arch 数值基准, 升入 scope)

struct GarchMParams {
    Real mu;      ///< 均值方程截距
    Real lambda;  ///< 风险溢价系数 (g(h_t) 的系数)
    Real omega, alpha, beta;   ///< GARCH(1,1) 方差参数 (同 Phase 7B)
};

GarchMResult estimate_garch_m(const std::vector<Real>& data,
                              GarchMForm form = GarchMForm::Volatility,
                              const GarchConfig& config = GarchConfig{});
///< 均值方程: r_t = μ + λ·g(h_t) + ε_t;  ε_t = r_t − μ − λ·g(h_t)
///< 结果含 λ 的 QMLE sandwich (Bollerslev-Wooldridge) SE (GM4: 非 OLS t 表)
```

**算法步骤**:

```
Step 1: 联合参数 θ = (μ, λ, ω, α, β); 多起始点沿用 GARCH 框架 (C5)
Step 2: 递归 (⚠️ GM5 均值-方差耦合, arch issue #269):
  2.1 h_1 = backcast (复用 Phase 7B G1 EWMA)
  2.2 t ≥ 2: ε_{t−1} = r_{t−1} − μ − λ·g(h_{t−1})   ← λ 变化重写整条 h 路径
      h_t = ω + α·ε²_{t−1} + β·h_{t−1}
  2.3 ε_t = r_t − μ − λ·g(h_t);  似然同 Phase 7B G3 完整形式
Step 3: SLSQP 约束同 GARCH(1,1) G4 + λ 无约束 (边界 [−10, 10])
Step 4: sandwich V = H⁻¹SH⁻¹ (复用); λ 的 SE 取 V[1][1] (GM4)
Step 5: 对照三步法 (verify_gm.R, 风险 #6):
  5.1 C++ vs arch 8.0 (主锚, form 三值逐一对偶, 1e-8~1e-10)
  5.2 C++ vs rugarch (archpow=1/2, 容差 1e-4 — solver 精度差异)
  5.3 fix() 互验: rugarch 估计值代入 ARCHInMean.fix() 求 llf, 隔离似然差 vs 优化器落点差
  ⚠️ 两包默认都重标定数据 (llf 差 n·log(s²)): arch rescale=False / rugarch 同关, 或 Jacobian 校正
```

**幻觉点映射**: GM1 (archpow 映射) / GM2 (arch 8.0 才有) / GM3 (log 变体有基准) / GM4 (sandwich SE) / GM5 (耦合递归)

---

## 3. M1: ARIMA + 线性 Granger (与 M4 并行)

> **命名空间**: `cpphub::v1::timeseries::arima` (纯标量, 无 Eigen3)
> **教材锚点**: Box-Jenkins 5th ed Ch.7 (CSS) / Brockwell-Davis 2016 ITSF §5.2 (innovations) / Hyndman fpp2 §8.1 (差分/drift) / Granger 1969 / Toda-Yamamoto 1995 / Lütkepohl 2005 **§3.6.1** p.102 (Granger Wald; ⚠️ 非 §7.2.2)

### 3.1 ARIMA 模型 (`arima_model.hpp`)

**对照库**: R `stats::arima` (CSS/CSS-ML) + statsmodels `tsa.arima.model.ARIMA` (innovations_mle)

**接口签名**:

```cpp
namespace cpphub::v1::timeseries::arima {

struct ArimaSpec {
    Size p, d, q;
    bool include_drift = false;   ///< 仅 d=1 有意义: 差分截距 (forecast::Arima 语义, AR5)
    ///< d>0 默认无均值项 (R stats::arima 忽略 include.mean; 决策 3)
};

struct ArimaParams {
    std::vector<Real> phi;   ///< AR: (1 − φ₁B − … − φ_p B^p)
    std::vector<Real> theta; ///< MA: (1 + θ₁B + … + θ_q B^q)  ⚠️ 正号参数化 (AR1)
    Real drift = 0.0;
    Real sigma2;
};

enum class ArimaMethod { CSS, CSS_ML, Innovations };

struct ArimaResult {
    ArimaParams params;
    Real loglik;            ///< 差分后数据似然: −(T−d)/2·[log(2πσ²)+1] (AR3/AR4 口径)
    Real aic;               ///< 基于 T−d 个观测 (AR4)
    std::vector<Real> residuals;
    std::vector<std::vector<Real>> vcov;   ///< sandwich (QMLE) 可选
    bool converged;
};

ArimaResult arima_fit(const std::vector<Real>& data, const ArimaSpec& spec,
                      ArimaMethod method = ArimaMethod::CSS_ML,
                      const ArimaConfig& config = ArimaConfig{});
///< config: n_starting_points (HR/CSS 预优化 + 随机重启, AR7/G-ADR6 同源)
///< Innovations 仅允许: 无缺失 + 无季节 + q ≥ 0 (AR8 限定, 对标 innovations_mle)
```

**算法步骤**:

```
Step 1: 差分 d 阶 → z_t (T−d 个观测); drift = 差分截距 (d=1 时可选)
Step 2: CSS (决策 1):
  2.1 n.cond = d + max(user_n.cond, p)   ⚠️ AR2: 与 q 无关 (arima.R L158-162;
      v1.0 的 max(p,q)+1 会导致与 R 起始点不可对照 — 已审计修正)
  2.2 ε_t = 0 (t < n.cond 前 innovations 视为 0)
  2.3 CSS(θ) = Σ_{t≥n.cond} ε_t²(θ), 递归: ε_t = z_t − Σφᵢ z_{t−i} − Σθⱼ ε_{t−j} (−drift)
Step 3: CSS-ML: CSS 解为起始 → 高斯似然 BFGS/SLSQP 精化 (C5);
  loglik = −(T−d)/2·[log(2πσ̂²) + 1], σ̂² = SSR/(T−d)  (AR3: 对照统一高斯形式;
  R 对纯 CSS 只报 "part log-likelihood", 不与 ML 混比)
Step 4: 多起始点 (AR7): {HR 初值, CSS 解, CSS 解±扰动, 随机重启} 各精化取最优 loglik
  (似然多局部最优是已知问题, 否则 1e-10 对照假失败, 风险 #3)
Step 5: Innovations MLE (`innovations_mle.hpp`): B&D 2016 ITSF §5.2 算法
  (Levinson-Durbin 型革新递归), 无缺失时与 Kalman 精确 MLE 同似然;
  exog/trend 经 GLS 另行处理 (非同一似然, 不互比)
Step 6: 残差诊断复用 Phase 7A (LB/JB)
```

**幻觉点映射**: AR1 ((1+θB) 与 R/statsmodels 同号) / AR2 (n.cond 与 q 无关) / AR3 (loglik 口径) / AR4 (AIC 基于 T−d) / AR5 (drift 语义) / AR6 (method 配对对照, CSS 与精确 MLE 不互照 1e-10) / AR7 (多起始) / AR8 (innovations 限定)

**测试要点**: CSS vs R method="CSS" 1e-10 (method 配对); CSS-ML vs R "CSS-ML" 1e-8; innovations vs statsmodels 1e-10 (无缺失非季节); MA(1) θ₁ 符号与 R 逐系数一致 (AR1); d=1 drift vs forecast::Arima; T−d 口径 AIC; 局部最优逃逸 (模拟 ARMA(2,2) 已知双峰案例)。

### 3.2 Granger 因果检验 (`granger_test.hpp`)

**对照库**: statsmodels `grangercausalitytests` (4 统计量) + R vars::causality (系统 F 交叉)

**接口签名**:

```cpp
namespace cpphub::v1::timeseries::granger {

struct GrangerResult {
    Real f_stat, f_p;              ///< ssr_ftest (主)
    Real params_f_stat, params_f_p;///< params_ftest (第二 F 构造)
    Real chi2_stat, chi2_p;        ///< ssr_chi2test
    Real lr_stat, lr_p;            ///< lrtest
    // Toda-Yamamoto
    Real ty_wald_stat, ty_wald_p;  ///< df = k (GR2), 非 k+d_max
    Size df1, df2;
    // HAC 稳健
    Real hac_wald_stat, hac_wald_p;///< NW vcov 上的 Wald (GR5)
    std::string summary;           ///< 注明: 均非异方差稳健 (除 HAC 版)
};

GrangerResult granger_test(const std::vector<Real>& cause,
                           const std::vector<Real>& effect,
                           Size lag = 1,
                           Size d_max = 0,       ///< TY 用: 外部单位根预检验给定 (GR3, 默认 0=不做 TY)
                           bool with_hac = true, Size hac_bandwidth = 0);
///< ⚠️ 显式 (cause, effect) 形参, 消灭 statsmodels "第二列 cause 第一列" 方向陷阱 (GR6)
///< effect 方程: effect_t = c + Σφᵢ effect_{t−i} + Σβⱼ cause_{t−j} + ε_t
```

**算法步骤**:

```
Step 1: 标准 F (GR1): F = ((RSS_R − RSS_U)/m) / (RSS_U/(T − k_u))
  双变量同阶 p: m = p, k_u = 2p + 1 (受约束方程 = 仅自身滞后 + 截距)
Step 2: χ² = T·(RSS_R/RSS_U − 1) ≈ F·m; LR = T·(log RSS_R − log RSS_U)
  (与 statsmodels 四统计量一一配对, 容差 1e-10)
Step 3: TY 增广 Wald (决策 5): 估计 VAR(k + d_max) 水平值,
  Wald 约束仅前 k 阶 cause 滞后 (GR4: 增广阶全部估计但不进约束矩阵), df = k (GR2)
  d_max 由调用方单位根预检验给定 (I(1) 惯例 d_max=1; GR3 不做模型内自适应)
Step 4: HAC-Wald (决策 6): NW 协方差 (复用 v1.5) 上的 Wald, χ²(k) (GR5)
```

**幻觉点映射**: GR1 (df 公式) / GR2 (TY df=k) / GR3 (d_max 外部) / GR4 (不丢样本) / GR5 (稳健须自建) / GR6 (方向) / GR7 (I(1) 水平值标准 F 失效 → 差分/TY)

---

## 4. M2: VAR + IRF/FEVD + Diebold-Yilmaz (依赖 M1)

> **命名空间**: `cpphub::v1::timeseries::var` (需 Eigen3, §8.2)
> **教材锚点**: Lütkepohl 2005 Ch.2 (IRF/FEVD) Ch.3 (估计) Ch.4 pp.146-150 (IC) / Sims 1980 / Koop-Pesaran-Potter 1996 / Pesaran-Shin 1998 *Economics Letters* 58:17-29 (⚠️ 非 J.Econometrics, V1) / Diebold-Yilmaz 2012 IJF
> **对照库**: statsmodels var_model (主) + R vars (交叉) + R Spillover (GFEVD/DY 主基准)

### 4.1 VAR 估计与滞后选择 (`var_model.hpp` / `var_select.hpp`)

```cpp
namespace cpphub::v1::timeseries::var {

struct VARSpec {
    Size lag = 0;                  ///< 0 => 由 IC 选择
    std::string trend = "c";       ///< "n"/"c"/"ct"
    Size max_lag = 0;              ///< IC 搜索上限 (0 => Schwert 型)
    bool same_sample_ic = true;    ///< V5: offset 机制强制同一样本量
};

struct VARResult {
    Eigen::MatrixXd coefficients;  ///< K×(K·p + k_trend) 逐方程同回归元 (决策 9)
    Eigen::MatrixXd sigma_u_mle;   ///< ML 版 Σ_u = SSR/T (V4: ÷T, 非 ÷(T−k))
    Real loglik, aic, bic, hqic, fpe, logdet;
    Real max_abs_eigenvalue;       ///< V9: max|eig(伴随矩阵)|
    bool is_strictly_stationary;   ///< 严格 <1 判定 (与 statsmodels ≤1 区分, 双输出)
    Size n_obs_used;
    // IRF/FEVD 结果 (§4.2)
};

VARResult var_fit(const MultivariateTSData& data, const VARSpec& spec = {});
///< 估计: 逐方程 OLS — 同阶同回归元下与系统 GLS 数值等价 (决策 9; B3 双盲:
///< statsmodels VAR.fit 的 method 参数在实现中未被引用, 文档+实现双重仅 'ols')

// IC 五式 (V4/V6, statsmodels var_model.py L2282-2303 已本地直读确认, 风险 #7 关闭):
//   fp = p·K² + K·k_trend;  ld = logdet(Σ_mle)
//   aic = ld + 2/T·fp;  bic = ld + ln(T)/T·fp;  hqic = ld + 2·lnln(T)/T·fp
//   fpe = ((T + df_model)/df_resid)^K · exp(ld)
//   df_model = p·K + k_trend (单方程); df_resid = T − df_model
// select_order: 对每个候选 p 用 offset 截断到同一样本 (V5), 否则 IC 因样本量错序

// Cholesky 与识别 (决策 11, V2):
//   P = Eigen::LLT(Σ).matrixL()   ← 下三角 (numpy/Eigen 约定)
//   ⚠️ R/MATLAB chol 返回上三角 — 误用使递归识别假设反转
//   API 暴露: set_identification_P(MatrixXd) 注入自定义 P + reorder_variables(indices)
```

**幻觉点映射**: V2 (三角方向) / V4 (ML 协方差 ÷T) / V5 (同样本) / V6 (FPE 指数 K, n* 含确定性项) / V9 (稳定性双输出) / V12 (不稳定 VAR 前置拦截 FEVD)

### 4.2 IRF 与 FEVD 双轨 (`irf.hpp` / `fevd.hpp`)

```
IRF (V3):
  Φ_0 = I;  Φ_l = Σ_{i=1}^{min(l,p)} Φ_{l−i}·A_i   (伴随递推)
  Θ_h[i,j]: 行 = 响应变量 i, 列 = 冲击变量 j (转置即全错)
  正交化: Ψ_h = Φ_h·P (P=下三角 Cholesky)
  置信带: block bootstrap 复用 v1.5 (决策 14); δ 法渐近 SE 二期
  ⚠️ V13: bootstrap 带与 δ 法带不同源, 容差断言仅适用点估计

FEVD Cholesky 轨 (行和精确 = 1):
  FEVD_{ij}(h) = Σ_{l<h} Ψ_l[i,j]² / Σ_{j'} Σ_{l<h} Ψ_l[i,j']²

GFEVD 轨 (决策 12/15, 自实现 — statsmodels 无, B1/B2 双盲确认):
  广义 (未正交化) 脉冲: 分子系数按目标框架二选一 (V8, v1.1 审计关键修正):
    DY 2012 框架 (默认, 溢出指数专用):
      ω_ij(h) = σ_jj⁻¹ · Σ_{l<h} (e_i' Φ_l Σ e_j)² / Σ_{l<h} e_i' Φ_l Σ Φ_l' e_i
      行归一化: θ̃_ij = ω_ij / Σ_{j'} ω_ij j'   (行和 = 1, V7)
    PS 1998 框架 (可选输出, 文献对照):
      原式 (12) 分子系数 = σ_ii⁻¹ (响应变量方差), 配未归一分母, 行和 ≠ 1
  ⚠️ 两框架行归一化后数值不同 (除非各 σ 相等), API 必须注明框架且不混用
```

### 4.3 DY 溢出指数 (`dy_spillover.hpp`)

```
输入: θ̃ (GFEVD DY 框架, H 步, 默认 H=10 可配置 — V10: 论文日度惯例非数学常数)
  TCI   = 100 · Σ_{i≠j} θ̃_ij / (N·H_effective)   (总溢出, H_effective=H)
  TO_j  = 100 · Σ_{i≠j} θ̃_ij / N   (j 发出)
  FROM_j= 100 · Σ_{j'≠j} θ̃_{j j'} / N   (j 接收)
  NET_j = TO_j − FROM_j
  滚动窗口: window 参数 (默认 120, 开放问题 d); 每窗口完整重估 VAR+GFEVD
主基准: R Spillover 包 g.fevd/G.spillover (1e-8)
```

**测试要点**: 系数/IC/IRF/FEVD vs statsmodels 1e-10; IC vs vars::VARselect 1e-8; Cholesky vs `np.linalg.cholesky` (下三角) 与 R `t(chol(Σ))` 等价断言; GFEVD vs R Spillover 1e-8; PS vs DY 框架差异断言 (σ 不等时归一化后数值不同); FEVD 行和=1 (1e-12); 不稳定 VAR 抛异常 (V12); DY H=10 vs H=50 敏感性。

---

## 5. M3: 协整三件套 (依赖 M2)

> **命名空间**: `cpphub::v1::timeseries::cointegration` (需 Eigen3)
> **教材锚点**: Engle-Granger 1987 / Johansen 1988 *JEDC* 12:231-254 (⚠️ 非 JASA, CI11) / MacKinnon 1994 JBES 12(2):167-176 / Osterwald-Lenum 1992 OBES 54(3):461-472 / Ericsson-MacKinnon 2002 Econometrics J. 5(2):285-318 / Lütkepohl 2005 p.292

### 5.1 Engle-Granger (`engle_granger.hpp`)

```
Step 1: OLS y₁ ~ y₂ (+trend), 残差 û
Step 2: û 的 ADF 回归 (nc 形式), 统计量 t(α̂)
Step 3: 临界值 = MacKinnon 1994 协整响应面 (mackinnon_coint_cv.hpp):
  CV = β_∞ + β₀/T + β₁/T² + β₂/T³ + β₃/T⁴  (N≥2 用 5 系数 4 次; ⚠️ 与 ADF N=1 的
  2010 表 4 系数分文件存储, 决策 18)
  按 N = 系统变量总数 与 step1 trend 索引 (CI1)
Step 4: p 值 = 1994 渐近近似 (复刻 statsmodels); ⚠️ CI2: statsmodels 的 p (渐近) 与
  cv (含 nobs 小样本修正) 不同源 (issue #4138 官方确认), 可现 "过 1% cv 而 p>1%",
  API 文档显式声明两者差异, 测试断言分别对照
Step 5: 方向依赖 (CI3): LHS 选择改变统计量 → 提供双方向输出; 方向无关替代 = PO-Pz
```

### 5.2 Johansen (`johansen_test.hpp`)

```
API 等价 coint_johansen(endog, det_order, k_ar_diff): det_order ∈ {−1, 0, 1} (3 情形, CI4;
  5 情形归 VECM 类)
算法 (Lütkepohl p.292):
  1. 辅助回归 Δy_t ~ Δy 滞后 (R₀t), y_{t−1} ~ Δy 滞后 (R₁t)   ← transitory 形式
     ⚠️ CI7: urca spec="longrun" 是 ΠX_{t−k} 形式, 与 statsmodels 对照必须用 transitory
  2. S₀₀=ΣR₀R₀'/T, S₀₁, S₁₁; 广义特征值问题 |λS₁₁ − S₀₁S₁₀S₀₀⁻¹... | 的 λ̂ 降序 (CI6)
  3. λ_trace(r) = −T·Σ_{i=r+1}^N ln(1−λ̂ᵢ);  λ_max(r) = −T·ln(1−λ̂_{r+1})
     T = 剔除预样本后的有效样本 (CI6)
  4. rank 从 r=0 逐级检验至首次不拒绝
临界值 (决策 19):
  主录 OL1992 (ol1992_cv.inc: urca 源码常量转录 constexpr + static_assert, 仅 N<11)
  statsmodels 内嵌表 (MHM96, c_sjt/c_sja — B1 双盲: 全库仅此一套) 双对照
  ⚠️ 实施第一步: 双库数值 diff 报告 (verify_johansen_diff.R) → 决定主对照并冻结
```

### 5.3 Phillips-Ouliaris (`phillips_ouliaris.hpp`)

```
Pu (残差基, 方向依赖) / Pz (对协整向量归一化不变, 方向无关) 双实现, Pz 优先 (决策 20, CI12)
长期方差: Bartlett 核 (复用 long_run_variance)
对照: urca ca.po (demean/lag/type 参数配对, 1e-8)
```

### 5.4 VECM (`vecm_model.hpp`)

```
5 情形 {n, co, ci, lo, li} (决策 17/CI4)
Δy_t = α(β'y_{t−1} + 协整常数/趋势) + Σ Γᵢ Δy_{t−i} + 外生 + ε_t
β 归一化 (决策 21, CI8):
  默认: 前 r×r 块 = I_r ("Phillips 归一" 文献惯称, statsmodels vecm.py L1002 实现)
  开关: urca 式首变量 = 1
  ⚠️ coint_johansen 函数本身不归一 (特征向量符号任意, issue #5517)
  ⚠️ 对照测试比较 β 张成空间: P = β(β'β)⁻¹β' 投影矩阵 (非逐元素)
Π = αβ' = ΣAᵢ − I (CI9 符号: αᵢ < 0 = 正向偏离被拉回)
ECT 系数 t 检验: Ericsson-MacKinnon 2002 响应面查表 (决策 21, CI10; 非标准 t)
```

**M3 测试要点**: EG vs statsmodels coint 1e-10 (统计量 + cv 分列断言); Johansen eig/lr1/lr2 vs coint_johansen 1e-10 (transitory); OL1992 表 static_assert; VECM 系数 vs statsmodels VECM 1e-10 + β 投影空间; ECT t vs EM2002 表; PO vs urca 1e-8; 阶数/rank 逐级检验序列正确性。

---

## 6. M4: MIDAS 混频回归 (与 M1 并行)

> **命名空间**: `cpphub::v1::timeseries::midas` (纯标量, 无 Eigen3)
> **教材锚点**: Ghysels-Santa-Clara-Valkanov 2006 J.Econometrics 131:59-95 / 综述 Ghysels-Sinko-Valkanov **2007** Econometric Reviews 26(1):53-90 (⚠️ "2006 综述"不存在, MD5) / Foroni-Marcellino-Schumacher 2011 DP 35 → 2015 JRSS-A 178(1) (U-MIDAS, MD4) / Clements-Galvao 2008 JBES (MIDAS-AR) / Ghysels-Qian 2019 E&S 9:1-16 (分离估计) / Kvedaras-Zemlys **2012** EL 116(2):250-254 (权重检验, MD6)
> **对照库**: R midasr 0.9 **唯一主数值基准** (B5 双盲: Python 生态四角度验证无维护良好 MIDAS 回归实现; arch 仅 MIDASHyperbolic 波动率过程; PyPI MIDASpy 为缺失值插补包, 名称撞车)

### 6.1 混频对齐 (`mixed_freq_data.hpp`)

```
mls(x, m, h): 列 = x[m·t − h], t = 1..n (低频索引)
  ⚠️ MD3: lag0 = x_{tm} = 低频期内**最后一个**高频观测 (期末对齐);
  "期初可得信息" 预测必须从 k ≥ m 起窗
  m ∤ length(x) → 报错 (midasr 同)
MixedFreqData { 低频 y (n), 高频 x (N = n·m + 预烧), m, 最大低频滞后 }
```

### 6.2 权重族 (`midas_weights.hpp`)

```cpp
enum class MidasWeight { Nealmon, NBeta, AlmonP, PolyStep, HarStep };
// ⚠️ 两套网格并存 (决策 22, MD1):
//   nealmon: i = 1..d  (从 1 起!)  w_i = λ₁·exp(Σ_k λ_{k+1}·i^k) / Σ_{i=1}^d exp(·)
//            λ₁ = δ 为自由总尺度 (Σw = δ, 独立线性参数, MD2)
//   nbetaMT: xi = (0..d−1)/(d−1)  (从 0 起, 与 nealmon 不同!) + κ₁κ₂ → β pdf → 归一化
// log-sum-exp 防溢出 (决策 25, MD7): 裸 exp 在 λ₂·d² ≳ 709 溢出;
//   logΣexp(η) = max(η) + log Σ exp(η − max(η)); 非溢出区间与裸公式差 < 1e-14 (断言)
// λ/a/b 全自由实数 (无非负约束 — exp/β 核归一化天然 ≥ 0, MD7)
```

### 6.3 模型与估计 (`midas_model.hpp`)

```
模型 (决策 22):
  MIDAS-DL:  y_t = μ + Σ_{j,i} θ_{j,i}·w_i(λ)·x_{(t−j)m−i} + ε_t   (逐字复刻 midas_r 公式)
  MIDAS-AR:  + Σ ρᵢ y_{t−i};  AR*: 约束 ρ = φ^(i−1)·(1−φ) 参数化
  U-MIDAS:   每个高频滞后无约束系数 (纯 OLS)

估计 — 集中化 NLS (决策 23, Ghysels-Qian 2019):
  外层: SLSQP (ADR-018) 仅优化 λ (权重超参)
  内层: 给定 λ, 权重固定 → X(λ) 线性于 (δ, 截距, AR 系数) → OLS(QR) 解析解
        → 集中化 SSR(λ)
  多起点: λ 网格 × {递减, 驼峰, 均匀} 三形状 (MD8: λ=0 起点平坦陷阱)
备选: 联合 NLS (全参数一起, 对照 midas_r Ofunction 语义)
诊断 (midas_diagnostics.hpp): 残差 LB/JB 复用; hAh 权重检验 (K-Z 2012) 三列全录
```

**夹具策略 (风险 #4, 决策 24)**:
```
verify_midas.R 生成规范 (全部固化 CSV):
  - midasr 0.9 锁版本; set.seed(42); optim control = list(reltol=1e-12, maxit=10000)
  - 每夹具记录: 收敛码 + SSR 双判据 (convergence && SSR 一致才采信)
  - midas_u 纯 OLS 锚: 1e-10 (无优化器敏感性)
  - midas_r NLS 锚: 1e-6~1e-8 分层 (BFGS 收敛容差所限)
  - hAh_test 输出三列全录 (统计量/p/df)
```

**M4 测试要点**: 权重逐点 vs midasr 1e-12 (含 i=1 起点与 xi=0 起点两套); log-sum-exp 差 <1e-14; mls 期末对齐 (lag0 = x_tm) 断言; U-MIDAS vs midas_u 1e-10; MIDAS-DL/AR vs 夹具 1e-6~1e-8; 集中化 NLS vs 联合 NLS 一致性 (同一最优); 多起点逃逸。

---

## 7. 端到端集成测试 (`test_integration_phase7c.cpp`)

| # | 场景 | 链路 | 验证点 |
|---|------|------|--------|
| 1 | 单位根诊断全链 | ADF/DF-GLS/NP/ZA 组合 → 差分决策 → ARIMA → 残差 LB | M0+M1+7B 复用; 断点→分段平稳 |
| 2 | Granger 因果链 | I(1) 双序列: 差分标准 F vs 水平 TY 增广 Wald 对比 | GR7 (标准 F 失效场景) |
| 3 | VAR→DY 溢出 | 多资产收益率 → IC 选阶 → 稳定性 → IRF/FEVD 双轨 → DY 矩阵 + 滚动 | M2 全链; V12 拦截 |
| 4 | 协整→VECM | EG/Johansen rank → VECM → ECT 显著性 → β 投影空间 | M3 全链 |
| 5 | MIDAS 混频预测 | 月度 y + 日度 x → MIDAS-DL vs U-MIDAS → 预测精度 (MZ/DM 复用 7B) | M4; MD3 期初起窗 |
| 6 | GARCH-M 风险溢价 | 收益率 → 三变体 → λ sandwich 显著性 → vs 无 M 模型 | M0; GM4/GM5 |

---

## 8. 依赖关系与复用

### 8.1 模块依赖图

```
M0 回填
├── ng_perron_test ← df_gls (GLS 去势复用) + np_table1.inc + Stata MAIC 对照
├── zivot_andrews_test ← ADF 引擎 (滞后复用) + za1992_cv.inc
└── garch_m_model ← garch_model (backcast/似然/SLSQP/sandwich 全复用, GM5 耦合递归改)
M1 ARIMA/Granger
├── arima_model ← SLSQP (C5) + hannan_rissanen (新, 起始值)
├── innovations_mle ← 纯新 (B&D 递归)
└── granger_test ← OLS + HAC/NW (v1.5 复用) + Wald
M2 VAR/DY
├── var_model/var_select ← OLS + Eigen3 (LLT/特征值)
├── irf ← block bootstrap (v1.5)
├── fevd ← var_model (Φ_l 递推)
└── dy_spillover ← fevd (GFEVD)
M3 协整
├── engle_granger ← OLS + ADF 引擎 + mackinnon1994_coint.inc
├── johansen_test ← Eigen3 自伴广义特征值 + ol1992_cv.inc
├── phillips_ouliaris ← long_run_variance (7B 复用)
└── vecm_model ← johansen (rank) + OLS + EM2002 表
M4 MIDAS
├── mixed_freq_data ← 纯新
├── midas_weights ← 纯新 (log-sum-exp)
└── midas_model ← SLSQP (C5 外层) + OLS/QR (v1.5 内层) + midas_diagnostics ← Phase 7A
```

### 8.2 Eigen3 隔离 (C1/ADR-013)

| 模块 | Eigen3 | 载体 |
|------|--------|------|
| M0 (NP/ZA/GM) | 否 | `cpphub_timeseries` (现有 INTERFACE) |
| M1 (ARIMA/Granger) | 否 | `cpphub_timeseries` |
| M2 (VAR/IRF/FEVD/DY) | **是** | 新 INTERFACE target `cpphub_timeseries_mat` (链接 `Eigen3::Eigen`) |
| M3 (协整/VECM) | **是** | `cpphub_timeseries_mat` |
| M4 (MIDAS) | 否 | `cpphub_timeseries` |

单变量头文件**不得** include Eigen; var/cointegration 头文件不得被单变量模块反向依赖。物理目录均在 `include/cpphub/timeseries/` (ADR-017), 以 CMake target 划隔离边界 (ADR-013 精神)。

### 8.3 C ABI (C3)

新增 C 导出 (若 Python 绑定/CDLL 需要) 统一 `cpphub_v1_7_*` 前缀; M0-M4 默认仅头文件消费, 不强制新增 C ABI 面。

---

## 9. 幻觉点核查清单 (实施时逐点验证, verify 脚本逐条落盘)

> 核查模式沿用 Phase 7B: probe_*.py 打印对照库源码行 + 小样本手算 + .inc 黄金值断言。
> 编号总计: AR 8 + GR 7 + V 13 + CI 12 + MD 8 + NP 6 + ZA 5 + GM 5 = **64 编号**
> (调研报告头部 "56 项" 为汇总口径统计, 与分表编号数不一致, 以本表为准)。

### 9.1 ARIMA/Granger (AR1-AR8, GR1-GR7)

| ID | 影响级 | 核查方式 | 容差 |
|----|--------|---------|------|
| AR1 (MA 同号 (1+θB)) | 极高 | verify_arima.R: R stats::arima MA 系数逐个对照 | 1e-10 |
| AR2 (n.cond=d+max(user,p), 与 q 无关) | 极高 | arima.R L158-162 probe + CSS SSR 逐位对照 | 1e-10 |
| AR3 (loglik 口径) | 高 | R "part log-likelihood" vs 本实现完整高斯型, 文档声明 | 1e-10 |
| AR4 (AIC 基于 T−d) | 高 | R AIC 数值对照 (SARIMAX simple_differencing 同理) | 1e-10 |
| AR5 (drift=差分截距) | 中 | forecast::Arima 语义对照 | 1e-8 |
| AR6 (method 配对) | 极高 | CSS↔CSS, ML↔ML 矩阵; 混配禁用 | - |
| AR7 (多起始) | 高 | 合成双峰似然案例逃逸断言 | - |
| AR8 (innovations 限定) | 高 | 无缺失/无季节 guard + statsmodels innovations_mle 对照 | 1e-10 |
| GR1 (F df 公式) | 高 | 手算小样本 + statsmodels ssr_ftest | 1e-10 |
| GR2 (TY df=k) | 极高 | Zapata-Gil 1999 数值例 | 1e-8 |
| GR3 (d_max 外部) | 中 | API 契约 (无自适应) | - |
| GR4 (增广阶不进约束) | 高 | 约束矩阵单元测试 | - |
| GR5 (非稳健默认) | 高 | HAC 版 vs 标准版差异断言 | 1e-8 |
| GR6 (方向) | 极高 | 显式 (cause,effect) 形参 + statsmodels 方向复现断言 | 1e-10 |
| GR7 (I(1) 失效) | 中 | 集成场景 2 | - |

### 9.2 VAR/DY (V1-V13)

| ID | 影响级 | 核查方式 | 容差 |
|----|--------|---------|------|
| V1 (PS 1998 = EL) | 低 | 文献引用正确性 (spec 内已冻结) | - |
| V2 (下三角) | 极高 | Eigen LLT vs np.linalg.cholesky vs R t(chol) 三方 | 1e-12 |
| V3 (Θ[i,j] 方向) | 极高 | 单位冲击脉冲响应方向断言 | 1e-12 |
| V4 (Σ ÷T) | 极高 | info_criteria 五式 vs var_model.py L2282-2303 | 1e-10 |
| V5 (同样本 offset) | 高 | select_order vs statsmodels | 1e-10 |
| V6 (FPE 形状) | 高 | ((T+df)/df')^K·exp(ld) 逐项 | 1e-10 |
| V7 (GFEVD 行归一) | 极高 | 行和=1 断言 + R Spillover | 1e-8 |
| V8 (双系数框架) | 极高 | σ_jj⁻¹ (DY) vs σ_ii⁻¹ (PS) 两输出差异断言 | 1e-8 |
| V9 (稳定性双输出) | 高 | max eig + 严格 <1 | 1e-12 |
| V10 (H 可配) | 中 | H=10/50 敏感性测试 | - |
| V11 (GIRF 隐含假设) | 低 | 文档注明 (Kim 2013) | - |
| V12 (不稳定拦截) | 高 | 异常路径测试 | - |
| V13 (bootstrap≠δ法) | 中 | 容差断言仅点估计 | - |

### 9.3 协整 (CI1-CI12)

| ID | 影响级 | 核查方式 | 容差 |
|----|--------|---------|------|
| CI1 (1994 响应面按 N 索引) | 极高 | mackinnon1994_coint.inc vs statsmodels c_sjt 系数 | 精确/1e-12 |
| CI2 (p/cv 不同源) | 高 | issue #4138 复现: 分列断言 + API 文档 | - |
| CI3 (EG 方向依赖) | 中 | 双方向输出差异断言 | 1e-10 |
| CI4 (3 vs 5 情形) | 高 | API 边界 (coint_johansen 3 / VECM 5) | - |
| CI5 (statsmodels 一套表) | 高 | 探针已确认 (B1); ol1992 vs MHM96 双表 diff 报告 | 精确 |
| CI6 (迹公式+λ降序+有效 T) | 极高 | coint_johansen lr1/lr2 对照 | 1e-10 |
| CI7 (transitory 对照) | 极高 | urca ca.jo spec="transitory" 对照 | 1e-8 |
| CI8 (β 投影空间对照) | 极高 | P=β(β'β)⁻¹β' 断言 (非逐元素) | 1e-10 |
| CI9 (Π 符号) | 中 | αᵢ<0 拉回方向断言 | - |
| CI10 (EM2002 表) | 高 | ericsson_mackinnon_cv.inc 查表 | 精确 |
| CI11 (文献出处) | 低 | spec 已冻结 (JEDC/JAE) | - |
| CI12 (Pu/Pz 双实现) | 高 | urca ca.po 双 type 对照 | 1e-8 |

### 9.4 MIDAS (MD1-MD8)

| ID | 影响级 | 核查方式 | 容差 |
|----|--------|---------|------|
| MD1 (nealmon i=1 起) | 极高 | midasr nealmon 逐点权重 | 1e-12 |
| MD2 (δ 独立线性参数) | 高 | 集中化内层 OLS 解析消去断言 | 1e-12 |
| MD3 (lag0=期末) | 极高 | mls(x, m, 0) == x[m·t] 断言 | 精确 |
| MD4 (U-MIDAS 出处) | 低 | 文献已冻结 (F-M-S 2011/2015) | - |
| MD5 (综述 2007) | 低 | 同上 | - |
| MD6 (K-Z 2012) | 低 | 同上 | - |
| MD7 (log-sum-exp) | 极高 | 非溢出区间裸公式差 + 溢出区间有限值断言 | <1e-14 |
| MD8 (多起点) | 高 | λ=0 平坦陷阱逃逸测试 | - |

### 9.5 回填 (NP1-NP6, ZA1-ZA5, GM1-GM5)

| ID | 影响级 | 核查方式 | 容差 |
|----|--------|---------|------|
| NP1 (Econometrica) | 低 | 文献冻结 | - |
| NP2 (τ_T 公式) | 极高 | 4 源已裁决; Stata 逐 k MAIC 对照 | 1e-10 |
| NP3 (AR 谱对差分) | 极高 | 对水平拟合爆炸反例断言 | - |
| NP4 (MPT 分情形 +14.5) | 极高 | 趋势情形末项系数断言 (+14.5 ≠ 13.5) | 1e-12 |
| NP5 (方向+Table 1) | 极高 | np_table1.inc static_assert | 精确 |
| NP6 (Stata 非基准) | 中 | 仅 MAIC/σ̂² 列对照 | 1e-10 |
| ZA1 (三库 lag 策略) | 高 | 双模式 (固定/Baum) 分列对照 | 1e-10/1e-8 |
| ZA2 (trim 参数化) | 高 | urca 复现 trim 放开 | 1e-8 |
| ZA3 (DT 断点后重新计时) | 极高 | DT_t=(t−Tb)·1(t>Tb) 单元断言 | 精确 |
| ZA4 (双临界值表) | 高 | za1992 主 + MC (c 1%=−5.27644 ≠ −5.83) | 精确 |
| ZA5 (取最负) | 极高 | 单调断点路径 min 断言 | 1e-12 |
| GM1 (archpow 1=σ/2=σ²) | 极高 | rugarch 手册映射断言 | - |
| GM2 (arch 8.0 才有) | 中 | 版本 guard (文档) | - |
| GM3 (log 变体有基准) | 高 | arch form='log' 对照 | 1e-8 |
| GM4 (sandwich SE) | 高 | BW robust vs arch | 1e-10 |
| GM5 (耦合递归) | 极高 | λ 扰动 → h 路径全变断言 (issue #269) | 1e-10 |

---

## 10. 验收标准 (G 门禁)

### 10.1 编译与测试
- [ ] MSVC Release: 全量通过 (新增 ~263 + 现有 2207 ≈ 2470)
- [ ] A 站 GCC: fresh clone + rebuild + ctest (submodules: recursive — CI 教训)
- [ ] B 站 GCC: fresh clone + rebuild + ctest
- [ ] 三平台无数值偏差 (容差分层表 §1.3)

### 10.2 数值基准 (§1.3 全表)
- [ ] M0: Stata MAIC 1e-10 / NP 恒等式 1e-12 / ZA 双模式 / GM 双锚 + fix 三步法
- [ ] M1: R CSS method 配对 1e-10 / innovations 1e-10 / Granger 4 统计量 1e-10
- [ ] M2: statsmodels VAR 1e-10 / vars 交叉 1e-8 / GFEVD vs Spillover 1e-8
- [ ] M3: coint/coint_johansen 1e-10 / OL1992 精确 / urca 交叉 1e-8 / β 投影空间
- [ ] M4: midas_u 1e-10 / midas_r 夹具 1e-6~1e-8 / 权重逐点 1e-12 / LSE <1e-14

### 10.3 幻觉点覆盖
- [ ] §9 五域 64 编号全部核查, probe/verify 脚本逐条落盘
- [ ] 前置工程: Johansen 双库 diff 报告冻结主对照 (决策 19)
- [ ] 临界值 .inc 六件套全部 static_assert

### 10.4 文档与流程对齐
- [ ] ADR-019 (26+3 项) 全部实施且无静默偏离 (冲突回溯 ADR 修订)
- [ ] ADR-017 命名空间 + ADR-013 Eigen 隔离 (§8.2 target 划界) 遵守
- [ ] C1-C7 兼容性约束逐项满足
- [ ] DEVELOPMENT_LOG.md 逐里程碑更新; PHASE7C_FINAL_ACCEPTANCE.md + checklist

### 10.5 Scope 边界 (不做清单)
- [ ] 未实现: SARIMA / wild bootstrap / SVAR / BVAR / TVP-VAR / ARDL-PSS / midas_nlpr/sp/qr/imidas_r / amweights / DCC / Kalman / 长记忆族
- [ ] 单变量模块零 Eigen 依赖 (grep 断言)

---

## 11. 风险与缓解 (调研风险表全量继承)

| # | 风险 | 等级 | 缓解 |
|---|------|------|------|
| 1 | Johansen OL1992 表录入 (无机器可读版) | 高 | urca 源码常量转录 constexpr + static_assert; **双库 diff 前置冻结** (verify_johansen_diff.R 为 M3 第一任务) |
| 2 | NP 无权威开源数值基准 | 高 | 原文公式 > Stata 逐 k MAIC (1e-10) > EViews 抄录 (Julia 降权); MZt≡MZα×MSB 恒等式自检 |
| 3 | ARIMA 似然多局部极值 → 1e-10 假失败 | 高 | 多起始 (HR/CSS/随机) + method 配对对照矩阵 (AR6) |
| 4 | midasr 唯一主基准 + BFGS 容差 | 中 | 夹具收紧 reltol=1e-12/maxit=10000/set.seed/锁版本 + midas_u 1e-10 锚 + convergence/SSR 双判据 |
| 5 | ZA 三库滞后策略不可逐位对齐 | 中 | 双模式 API (固定 lag 主 / Baum 对照); trim 参数化 |
| 6 | GM vs rugarch solver 系数差 | 中 | 1e-4 + fix() 互验三步法 + rescale 统一关闭 (llf 差 n·log(s²)) |
| 9 | PQ 2007 MAIC 与 NP 2001 微差 | 低 | 以 NP 2001 为准 (H2 [待定]) |

---

## 12. 实施顺序

```
M0 回填 (先行, 框架已热)
  ├── Week 1: NP (GLS 去势复用 + MAIC + np_table1.inc + Stata 对照)
  ├── Week 1-2: ZA (ADF 引擎复用 + za1992/za_mc 双表 + 双模式)
  └── Week 2: GM (GARCH 复用 + arch 8.0/rugarch 双锚 + fix 三步法)
M1 ∥ M4 (并行)
  ├── M1 Week 2-4: CSS → CSS-ML → HR 多起始 → innovations → Granger 4+TY+HAC
  └── M4 Week 2-4: mls → 权重族 (LSE) → U-MIDAS (midas_u 锚) → DL/AR 集中化 NLS → 夹具
M2 (依赖 M1)
  └── Week 4-6: cpphub_timeseries_mat target → var_fit/IC → IRF/FEVD 双轨 → GFEVD/DY (Spillover 对照)
M3 (依赖 M2)
  └── Week 6-8: 双库 diff 前置 → EG+1994 表 → Johansen+OL1992 → PO → VECM+EM2002
收尾 Week 8-9: 集成 6 场景 + 三平台验证 + PHASE7C_FINAL_ACCEPTANCE.md
```

---

## 13. 开放问题 ([待定] 项, 沿用 Phase 7B 推迟惯例)

| # | 问题 | 状态 | 处置 |
|---|------|------|------|
| a | DY 滚动窗口宽默认值 (120? 数据频率自适应?) | [待定] | M2 实施时以月度 120 / 日度 200 提案, ADR-020 或 spec v1.1 定稿 |
| b | H1: Julia MacroEconometricModels.jl NP MPT 常数情形末项系数疑偏 (v0.8.0 年轻包) | [待定-假设区] | 查证路径: NP 2001 出版版 p.1531 Table 1 脚注 vs ngperron.jl; 不阻塞 (Julia 仅半基准) |
| c | H2: Perron-Qu 2007 修正版 MAIC 与 NP 2001 小样本微差 | [待定-假设区] | v1.8 Gretl `adf --gls --perron-qu` 实测对照 |
| d | SARIMA / wild bootstrap (Granger) | scope 外 | v1.8+ |
| e | SVAR/BVAR/TVP-VAR | scope 外 | v1.8+ (DSGE 场景) |
| f | ARDL/PSS 边界检验 | scope 外 | v1.8+ |
| g | midas_nlpr/sp/qr/imidas_r, amweights, gompertz/nakagami/lcauchy/genexp | scope 外 | v1.8+ |
| h | 多元 GARCH (CCC/DCC) / Kalman | scope 外 | v1.8+ |

---

## 14. 关联文档

- [ADR-019: v1.7 实施边界 (26+3 项)](../../decisions/ADR-019_V17_TIMESERIES_BOUNDARY.md)
- [ADR_INDEX.md](../../decisions/ADR_INDEX.md) (ADR-016/017/018/019)
- [PHASE7C_RESEARCH.md](../../research/PHASE7C_RESEARCH.md) v1.1 (调研全文 + 审计记录 + 幻觉点详表)
- [ADR019_REVIEW_PILOT.md](../../research/ADR019_REVIEW_PILOT.md) v1.1 (R1-R4 pilot)
- [ASSERTION_EVIDENCE_FRAMEWORK.md](../../ASSERTION_EVIDENCE_FRAMEWORK.md) v1.1 (调研/审计规范)
- [DEVELOPMENT_WORKFLOW.md](../../DEVELOPMENT_WORKFLOW.md) v1.1 §2.0/§4.2 (阶段 0 与 R/G 双门禁)
- [PHASE7B_FINANCIAL_TS_SPEC.md](./PHASE7B_FINANCIAL_TS_SPEC.md) (v1.6 复用基础)
