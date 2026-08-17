# Phase 7C 审计验收 Checklist

> **版本归属**: v1.7 (Phase 7C M0-M4)
> **关联 Spec**: [PHASE7C_SPEC.md](./PHASE7C_SPEC.md) v1.1
> **关联 ADR**: ADR-013 / ADR-016 / ADR-017 / ADR-018 / ADR-019
> **编写日期**: 2026-08-17 (实施前冻结验收框架, 全部 `[ ]` 待实施后逐项审计签署)
> **验收人**: Scott (self-review, 符合项目 solo developer 工作流)

> **§0 编写期 spec 对齐审计 (2026-08-17)**: 本 checklist 由 spec v1.1 逐条映射生成, 编写过程执行文档对齐清点, 发现并修正 2 处计数偏差 — ① 新增头文件实列 **25 个** (devlog 曾记 "23": unit_root 3 + garch 1 + arima 4 + var 6 (含 multivariate_data) + cointegration 7 + midas 4, 以 spec §1.1 文件树逐行清点为准); ② verify 脚本实列 **19 个** (devlog 曾记 "20"), 另有临界值 `.inc` 6 件。spec 正文无总数声明无需改; devlog 为历史记录不回改, 以本表为准。其余对齐项: 测试 16 套 263 用例 ✓ / 幻觉点 64 编号 ✓ / 数值基准 28 行 ✓ / 集成 6 场景 ✓ / 风险 7 项 (#7/#8 已关闭不列) ✓ / 开放问题 8 项 (a-c 已裁决入 §10/§14) ✓。

---

## 使用说明

- 每项标注 `[ ]` 待验收，`[x]` 已通过，`[!]` 未通过（需附 issue 描述, 附录 A 模板）
- 幻觉点前缀: **AR** (ARIMA) / **GR** (Granger) / **V** (VAR/DY) / **CI** (协整) / **MD** (MIDAS) / **NP** / **ZA** / **GM** — 共 64 编号 (spec §9)
- 容差标注: 1e-10 = 相对误差 ≤ 1e-10; "精确" = 数学定义/表值完全一致
- 基准库 (版本冻结): Python `arch` 8.0.0 (GM 主锚) / `statsmodels` 0.14.6 / R `midasr` 0.9 (MIDAS 唯一主基准) / R `urca` (Johansen/ZA/PO) / R `rugarch` (GM archpow) / R `vars` (VAR 交叉) / R `Spillover` (GFEVD/DY 主基准) / Stata `dfgls` (NP 逐 k MAIC)
- **R 门禁回溯条款**: 实施期若发现实现与 spec 冻结公式冲突, 须回溯 ADR-019 修订流程并记录于 §14, **禁止静默改公式** (spec §0)
- **对照禁令** (开放问题裁决生效): Julia 常数情形 MPT 禁作对照 (H1, 仅趋势可用); v1.8 前禁以 arch DFGLS 的 k 选择充当 MAIC-PQ (H2)

---

## 1. 交付物完整性

### 1.1 新增头文件 (25 个, 按 ADR-017 命名空间)

**M0 回填 (4 个, unit_root + garch)**

| # | 文件 | 状态 | 备注 |
|---|------|------|------|
| 1.1.1 | `include/cpphub/timeseries/unit_root/ng_perron_test.hpp` | [ ] | MZα/MZt/MSB/MPT + MAIC/MBIC/seq-t |
| 1.1.2 | `include/cpphub/timeseries/unit_root/zivot_andrews_test.hpp` | [ ] | 三模型 A/B/C + trim + 断点搜索 |
| 1.1.3 | `include/cpphub/timeseries/unit_root/np_tables.hpp` | [ ] | NP 2001 Table 1 (constexpr + static_assert) |
| 1.1.4 | `include/cpphub/timeseries/garch/garch_m_model.hpp` | [ ] | GARCH(1,1)-M 三变体 |

**M1 ARIMA/Granger (4 个, 纯标量无 Eigen)**

| # | 文件 | 状态 | 备注 |
|---|------|------|------|
| 1.1.5 | `include/cpphub/timeseries/arima/arima_model.hpp` | [ ] | CSS + CSS-ML |
| 1.1.6 | `include/cpphub/timeseries/arima/innovations_mle.hpp` | [ ] | B&D 2016 §5.2 精确 MLE |
| 1.1.7 | `include/cpphub/timeseries/arima/hannan_rissanen.hpp` | [ ] | HR 起始值 |
| 1.1.8 | `include/cpphub/timeseries/arima/granger_test.hpp` | [ ] | F/χ²/LR + TY + HAC-Wald |

**M2 VAR/DY (6 个, 需 Eigen3 → cpphub_timeseries_mat)**

| # | 文件 | 状态 | 备注 |
|---|------|------|------|
| 1.1.9 | `include/cpphub/timeseries/var/multivariate_data.hpp` | [ ] | MultivariateTSData (C4) |
| 1.1.10 | `include/cpphub/timeseries/var/var_model.hpp` | [ ] | 逐方程 OLS + P 注入/重排 |
| 1.1.11 | `include/cpphub/timeseries/var/var_select.hpp` | [ ] | IC 五式 + 同样本 offset |
| 1.1.12 | `include/cpphub/timeseries/var/irf.hpp` | [ ] | 正交化 IRF + bootstrap 带 |
| 1.1.13 | `include/cpphub/timeseries/var/fevd.hpp` | [ ] | Cholesky + GFEVD 双轨 |
| 1.1.14 | `include/cpphub/timeseries/var/dy_spillover.hpp` | [ ] | TO/FROM/NET/TCI + 滚动 |

**M3 协整 (7 个, 需 Eigen3)**

| # | 文件 | 状态 | 备注 |
|---|------|------|------|
| 1.1.15 | `include/cpphub/timeseries/cointegration/engle_granger.hpp` | [ ] | EG 两步法 |
| 1.1.16 | `include/cpphub/timeseries/cointegration/johansen_test.hpp` | [ ] | 迹/最大特征值, 3 det 情形 |
| 1.1.17 | `include/cpphub/timeseries/cointegration/phillips_ouliaris.hpp` | [ ] | Pu/Pz 双实现 |
| 1.1.18 | `include/cpphub/timeseries/cointegration/vecm_model.hpp` | [ ] | 5 情形 + β 双归一 + ECT t |
| 1.1.19 | `include/cpphub/timeseries/cointegration/osterwald_lenum_cv.hpp` | [ ] | OL1992 表 (urca 转录) |
| 1.1.20 | `include/cpphub/timeseries/cointegration/mackinnon_coint_cv.hpp` | [ ] | MacKinnon 1994 协整响应面 |
| 1.1.21 | `include/cpphub/timeseries/cointegration/ericsson_mackinnon_cv.hpp` | [ ] | EM2002 ECT t 临界值 |

**M4 MIDAS (4 个, 纯标量无 Eigen)**

| # | 文件 | 状态 | 备注 |
|---|------|------|------|
| 1.1.22 | `include/cpphub/timeseries/midas/mixed_freq_data.hpp` | [ ] | mls 期末对齐 (C4) |
| 1.1.23 | `include/cpphub/timeseries/midas/midas_weights.hpp` | [ ] | 5 权重族 + log-sum-exp |
| 1.1.24 | `include/cpphub/timeseries/midas/midas_model.hpp` | [ ] | DL/AR/U-MIDAS + 集中化 NLS |
| 1.1.25 | `include/cpphub/timeseries/midas/midas_diagnostics.hpp` | [ ] | 残差诊断 + hAh 检验 |

### 1.2 基准验证脚本 (19 个) + 临界值表 (6 件, 沿用 7B 全提交惯例)

| # | 脚本 | 状态 | 备注 |
|---|------|------|------|
| 1.2.1 | `tests/fixtures/timeseries/verify_np_stata.py` | [ ] | Stata dfgls r(results) 逐 k MAIC → CSV |
| 1.2.2 | `tests/fixtures/timeseries/verify_np_semi.py` | [ ] | EViews/Julia 抄录半基准 (Julia 常数 MPT 禁用) |
| 1.2.3 | `tests/fixtures/timeseries/verify_za.py` | [ ] | statsmodels (Baum 模式) |
| 1.2.4 | `tests/fixtures/timeseries/verify_za.R` | [ ] | urca ur.za (固定 lag, trim 放开) |
| 1.2.5 | `tests/fixtures/timeseries/verify_gm.py` | [ ] | arch 8.0 ARCHInMean form 三值 |
| 1.2.6 | `tests/fixtures/timeseries/verify_gm.R` | [ ] | rugarch archpow=1/2 + fix() 三步法 |
| 1.2.7 | `tests/fixtures/timeseries/verify_arima.R` | [ ] | R stats::arima CSS/CSS-ML (method 配对) |
| 1.2.8 | `tests/fixtures/timeseries/verify_arima.py` | [ ] | statsmodels innovations_mle |
| 1.2.9 | `tests/fixtures/timeseries/verify_granger.py` | [ ] | grangercausalitytests 4 统计量 |
| 1.2.10 | `tests/fixtures/timeseries/verify_var.py` | [ ] | statsmodels VAR 系数/IC/IRF/FEVD |
| 1.2.11 | `tests/fixtures/timeseries/verify_var.R` | [ ] | vars::VAR + VARselect |
| 1.2.12 | `tests/fixtures/timeseries/verify_gfevd.R` | [ ] | R Spillover g.fevd/G.spillover |
| 1.2.13 | `tests/fixtures/timeseries/verify_eg.py` | [ ] | statsmodels coint |
| 1.2.14 | `tests/fixtures/timeseries/verify_johansen.py` | [ ] | statsmodels coint_johansen |
| 1.2.15 | `tests/fixtures/timeseries/verify_johansen_diff.R` | [ ] | **M3 前置任务**: 双库 diff 冻结主对照 |
| 1.2.16 | `tests/fixtures/timeseries/verify_vecm.py` | [ ] | statsmodels VECM (β 投影空间) |
| 1.2.17 | `tests/fixtures/timeseries/verify_po.R` | [ ] | urca ca.po |
| 1.2.18 | `tests/fixtures/timeseries/verify_midas.R` | [ ] | midasr 0.9 夹具 (收紧 control + seed) |
| 1.2.19 | `tests/fixtures/timeseries/verify_midas_u.R` | [ ] | midas_u 纯 OLS 锚 (1e-10) |
| 1.2.20 | `critical_values/np_table1.inc` | [ ] | NP 2001 Table 1 |
| 1.2.21 | `critical_values/za1992_cv.inc` | [ ] | ZA 论文表 (主) |
| 1.2.22 | `critical_values/za_mc_cv.inc` | [ ] | MC 表 (c 1%=−5.27644) |
| 1.2.23 | `critical_values/ol1992_cv.inc` | [ ] | OL1992 (urca 转录) |
| 1.2.24 | `critical_values/mackinnon1994_coint.inc` | [ ] | MacKinnon 1994 响应面 (N≥2) |
| 1.2.25 | `critical_values/em2002_ect_cv.inc` | [ ] | EM2002 ECT t |

### 1.3 测试套件 (16 套, ~263 用例)

| # | 测试套件 | 用例数 | 状态 | 里程碑 |
|---|----------|--------|------|--------|
| 1.3.1 | `test_ng_perron` | 18 | [ ] | M0 |
| 1.3.2 | `test_zivot_andrews` | 15 | [ ] | M0 |
| 1.3.3 | `test_garch_m_model` | 16 | [ ] | M0 |
| 1.3.4 | `test_arima_model` | 24 | [ ] | M1 |
| 1.3.5 | `test_innovations_mle` | 12 | [ ] | M1 |
| 1.3.6 | `test_granger_causality` | 16 | [ ] | M1 |
| 1.3.7 | `test_var_model` | 20 | [ ] | M2 |
| 1.3.8 | `test_var_irf_fevd` | 18 | [ ] | M2 |
| 1.3.9 | `test_dy_spillover` | 12 | [ ] | M2 |
| 1.3.10 | `test_engle_granger` | 14 | [ ] | M3 |
| 1.3.11 | `test_johansen_test` | 18 | [ ] | M3 |
| 1.3.12 | `test_vecm_model` | 16 | [ ] | M3 |
| 1.3.13 | `test_phillips_ouliaris` | 10 | [ ] | M3 |
| 1.3.14 | `test_midas_weights` | 16 | [ ] | M4 |
| 1.3.15 | `test_midas_model` | 18 | [ ] | M4 |
| 1.3.16 | `test_integration_phase7c` | 8 | [ ] | 端到端 6 场景 (§9) |

---

## 2. 编译与跨平台测试

### 2.1 主控站 (Windows MSVC Release)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 2.1.1 | CMake 配置成功 (含新 target `cpphub_timeseries_mat`) | [ ] | §8.2 Eigen 隔离 |
| 2.1.2 | MSVC Release 编译零警告零错误 (`/utf-8` 沿用) | [ ] | |
| 2.1.3 | 全量 ctest 通过 (新增 ~263 + 现有 2207 ≈ 2470) | [ ] | |
| 2.1.4 | 无现有测试退化 (Phase 1-7B 全部仍通过) | [ ] | |
| 2.1.5 | SLSQP 12/12 仍通过 (ADR-018 无退化) | [ ] | M1/M4 外层复用 |

### 2.2 A 站 (Ubuntu GCC)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 2.2.1 | fresh clone + rebuild (submodules: recursive) | [ ] | CI 教训 |
| 2.2.2 | GCC 编译零警告零错误 | [ ] | |
| 2.2.3 | ctest 全量通过 | [ ] | |
| 2.2.4 | 与主控站数值一致 (容差分层 §3-§7) | [ ] | |

### 2.3 B 站 (Ubuntu GCC)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 2.3.1 | fresh clone + rebuild | [ ] | |
| 2.3.2 | GCC 编译零警告零错误 | [ ] | |
| 2.3.3 | ctest 全量通过 | [ ] | |
| 2.3.4 | 与主控站数值一致 | [ ] | |

### 2.4 三平台一致性与 CI

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 2.4.1 | M0-M4 全部测试三平台无数值偏差 | [ ] | |
| 2.4.2 | GitHub Actions 全绿 (Windows bash shell + submodules 修复沿用) | [ ] | |

---

## 3. M0 数值基准对照 (回填三项)

### 3.1 Ng-Perron — vs Stata dfgls (确定性) + EViews/Julia 半基准

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 3.1.1 | 逐 k MAIC(k)/σ̂²(k) vs Stata dfgls r(results) | 1e-10 | [ ] | 唯一确定性免费对照 |
| 3.1.2 | τ_T(k) = β̂₀²·Σỹ²/σ̂²(k) (固定样本 T−k_max) | 1e-10 | [ ] | NP2 4 源裁决公式 |
| 3.1.3 | MZt ≡ MZα×MSB 恒等式自检 | 1e-12 | [ ] | |
| 3.1.4 | MPT 常数末项 −c̄=+7 / 趋势末项 +14.5 (≠13.5 断言) | 1e-12 | [ ] | NP4; H1 裁决后 Julia 常数 MPT 禁对照 |
| 3.1.5 | 四统计量 vs EViews 抄录 (Julia 趋势情形旁证) | 1e-6~1e-8 | [ ] | |
| 3.1.6 | 渐近 5% 临界值 vs NP Table 1 (常数 −8.10/−1.98/0.233/3.17; 趋势 −17.30/−2.91/0.168/5.48) | 精确 | [ ] | np_table1.inc |
| 3.1.7 | 拒绝方向: MZα/MZt 越负, MSB/MPT 越小 | 方向性 | [ ] | NP5 |
| 3.1.8 | AR 谱密度对 Δỹ (非水平) 拟合 | 1e-10 | [ ] | NP3 (水平拟合爆炸反例断言) |

### 3.2 Zivot-Andrews — vs statsmodels (Baum 模式) + urca (固定 lag)

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 3.2.1 | 统计量/断点 vs statsmodels (lag 冻结 Baum 预选) | 1e-10 | [ ] | |
| 3.2.2 | vs urca ur.za (固定 lag, trim 放开) | 1e-8 | [ ] | ZA1/ZA2 |
| 3.2.3 | DT_t = (t−Tb)·1(t>Tb) 断点后重新计时 | 精确 | [ ] | ZA3 |
| 3.2.4 | 统计量 = min t(α̂) (最负, 非 max) | 1e-12 | [ ] | ZA5 |
| 3.2.5 | 临界值: ZA1992 论文表 (主) | 精确 | [ ] | za1992_cv.inc |
| 3.2.6 | MC 表对照 (c 1%=−5.27644 ≠ −5.83) | 精确 | [ ] | ZA4 |
| 3.2.7 | trim 参数化 (默认 0.15, 上限 0.333) | 精确 | [ ] | |

### 3.3 GARCH-M — vs arch 8.0 (主锚) + rugarch (1e-4 + 三步法)

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 3.3.1 | form='vol' (g=√h) vs arch ↔ rugarch archpow=1 | 1e-8~1e-10 / 1e-4 | [ ] | 双锚 |
| 3.3.2 | form='var' (g=h) ↔ archpow=2 | 同上 | [ ] | GM1 |
| 3.3.3 | form='log' vs arch (log 变体数值基准) | 1e-8 | [ ] | GM3 |
| 3.3.4 | λ 的 QMLE sandwich SE (BW robust) | 1e-10 | [ ] | GM4 |
| 3.3.5 | 均值-方差耦合: λ 扰动重写 h 路径断言 | 1e-10 | [ ] | GM5 (issue #269) |
| 3.3.6 | fix() 互验三步法 (隔离似然差 vs 落点差) | 流程 | [ ] | rescale 统一关闭 |
| 3.3.7 | loglik 含完整常数项 | 1e-15 | [ ] | 沿用 G3 |

---

## 4. M1 数值基准对照 (ARIMA + Granger)

### 4.1 ARIMA — vs R stats::arima + statsmodels (method 配对)

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 4.1.1 | CSS vs R method="CSS" (method 配对, 禁混配) | 1e-10 | [ ] | AR6 |
| 4.1.2 | n.cond = d + max(user, p), 与 q 无关 | 1e-10 | [ ] | AR2 (arima.R L158-162) |
| 4.1.3 | CSS-ML vs R method="CSS-ML" | 1e-8 | [ ] | |
| 4.1.4 | MA 系数 (1+θB) 正号与 R 逐系数一致 | 1e-10 | [ ] | AR1 |
| 4.1.5 | loglik/AIC 基于 T−d 差分观测 | 1e-10 | [ ] | AR3/AR4 |
| 4.1.6 | d=1 drift = 差分截距 (forecast::Arima 语义) | 1e-8 | [ ] | AR5 |
| 4.1.7 | innovations vs statsmodels innovations_mle (无缺失非季节) | 1e-10 | [ ] | AR8 guard |
| 4.1.8 | 多起始点逃逸 (合成双峰似然案例) | 流程 | [ ] | AR7 |

### 4.2 Granger — vs statsmodels + 文献数值例

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 4.2.1 | 4 统计量 (ssr_ftest/params_ftest/ssr_chi2test/lrtest) 全对 | 1e-10 | [ ] | |
| 4.2.2 | 显式 (cause, effect) 方向 (statsmodels 第二列 cause 第一列复现断言) | 1e-10 | [ ] | GR6 |
| 4.2.3 | F df 公式: m=p, k_u=2p+1 手算对照 | 1e-10 | [ ] | GR1 |
| 4.2.4 | TY Wald df=k (增广阶不进约束矩阵) | 1e-8 | [ ] | GR2/GR4, Zapata-Gil 1999 |
| 4.2.5 | HAC-Wald (NW vcov) vs 标准版差异断言 | 1e-8 | [ ] | GR5 |
| 4.2.6 | I(1) 水平标准 F 失效场景 (集成场景 2) | 方向性 | [ ] | GR7 |

---

## 5. M2 数值基准对照 (VAR/IRF/FEVD/DY)

### 5.1 VAR 估计与 IC — vs statsmodels + vars

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 5.1.1 | 系数矩阵 (逐方程 OLS) vs statsmodels | 1e-10 | [ ] | 决策 9 |
| 5.1.2 | IC 五式 (aic/bic/hqic/fpe/logdet) vs var_model.py | 1e-10 | [ ] | V4/V6 |
| 5.1.3 | Σ_mle = SSR/T (÷T 非 ÷(T−k)) | 1e-10 | [ ] | V4 |
| 5.1.4 | select_order 同样本 offset | 1e-10 | [ ] | V5 |
| 5.1.5 | IC vs R vars::VARselect | 1e-8 | [ ] | |
| 5.1.6 | Cholesky 下三角 (Eigen LLT = np.linalg.cholesky = R t(chol)) | 1e-12 | [ ] | V2 |
| 5.1.7 | is_stable 双输出 (max|eig| + 严格 <1) | 1e-12 | [ ] | V9 |

### 5.2 IRF/FEVD 双轨 — vs statsmodels + R Spillover

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 5.2.1 | IRF Θ_h[i,j] 方向 (行=响应, 列=冲击) | 1e-12 | [ ] | V3 |
| 5.2.2 | Cholesky FEVD 行和精确=1 | 1e-12 | [ ] | |
| 5.2.3 | GFEVD (DY 框架 σ_jj⁻¹ + 行归一化) vs R Spillover | 1e-8 | [ ] | 决策 15 |
| 5.2.4 | PS 框架 (σ_ii⁻¹) 可选输出, 与 DY 归一化后数值差异断言 | 1e-8 | [ ] | V8 |
| 5.2.5 | 不稳定 VAR 拦截 FEVD (异常) | 流程 | [ ] | V12 |
| 5.2.6 | bootstrap 置信带仅点估计容差 | 流程 | [ ] | V13 |

### 5.3 DY 溢出指数 — vs R Spillover

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 5.3.1 | TCI/TO/FROM/NET vs G.spillover | 1e-8 | [ ] | |
| 5.3.2 | window 必填无默认 + 强制 window>2·N | 精确 | [ ] | §13-a 裁决 |
| 5.3.3 | 频率默认表 {日 200/周 200/月 60} × H {10/10 周/12 月} | 精确 | [ ] | 月 60 非 120 |
| 5.3.4 | H 可配置敏感性 (H=10 vs 50) | 方向性 | [ ] | V10 |
| 5.3.5 | 滚动窗口 step=1 每窗口全重估 | 流程 | [ ] | |

---

## 6. M3 数值基准对照 (协整)

### 6.1 Engle-Granger — vs statsmodels coint

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 6.1.1 | 统计量 vs statsmodels coint | 1e-10 | [ ] | |
| 6.1.2 | MacKinnon 1994 响应面 (5 系数 4 次, N≥2, 与 2010 ADF 表分文件) | 1e-12 | [ ] | CI1, 决策 18 |
| 6.1.3 | p 值 (渐近) 与 cv (小样本修正) 分列断言 + API 文档声明不同源 | 流程 | [ ] | CI2 (#4138) |
| 6.1.4 | 双方向输出差异 | 1e-10 | [ ] | CI3 |

### 6.2 Johansen — vs statsmodels + urca (diff 前置)

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 6.2.1 | **前置: verify_johansen_diff.R 双库 diff 报告冻结主对照** | 流程 | [ ] | 决策 19 第一任务 |
| 6.2.2 | eig/lr1/lr2 vs coint_johansen (transitory 形式) | 1e-10 | [ ] | CI7 |
| 6.2.3 | λ̂ 降序 + 有效 T 口径 | 1e-10 | [ ] | CI6 |
| 6.2.4 | OL1992 表 vs urca 源码常量 (static_assert) | 精确 | [ ] | ol1992_cv.inc |
| 6.2.5 | statsmodels 内嵌表 (MHM96) 双对照 | 精确 | [ ] | CI5/B1 |
| 6.2.6 | 3 det 情形 API 边界 (5 情形归 VECM) | 流程 | [ ] | CI4 |

### 6.3 VECM + PO — vs statsmodels + urca

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 6.3.1 | VECM 系数 vs statsmodels (5 情形) | 1e-10 | [ ] | |
| 6.3.2 | β 对照用投影矩阵 P=β(β'β)⁻¹β' (非逐元素) | 1e-10 | [ ] | CI8 |
| 6.3.3 | β 双归一 (前 r×r=I_r 默认 / urca 首变量开关) | 1e-10 | [ ] | 决策 21 |
| 6.3.4 | ECT t 检验 EM2002 查表 | 精确 | [ ] | CI10 |
| 6.3.5 | Π=αβ' 符号方向 (αᵢ<0 拉回) | 方向性 | [ ] | CI9 |
| 6.3.6 | Pu (方向依赖) / Pz (方向无关) vs urca ca.po | 1e-8 | [ ] | CI12, 决策 20 |

---

## 7. M4 数值基准对照 (MIDAS)

### 7.1 权重族与对齐 — vs midasr 0.9

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 7.1.1 | nealmon 权重逐点 (i=1..d 从 1 起) | 1e-12 | [ ] | MD1 |
| 7.1.2 | nbetaMT (xi 从 0 起, 两套网格并存) | 1e-12 | [ ] | |
| 7.1.3 | δ 独立线性参数 (Σw=δ, 内层解析消去) | 1e-12 | [ ] | MD2 |
| 7.1.4 | mls lag0 = 期末 x_{tm} 断言 | 精确 | [ ] | MD3 |
| 7.1.5 | log-sum-exp 与裸公式差 (非溢出区间) | <1e-14 | [ ] | 决策 25/MD7 |
| 7.1.6 | 溢出区间 (λ₂·d²≳709) 有限值断言 | 流程 | [ ] | |

### 7.2 模型估计 — vs midasr 夹具

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 7.2.1 | U-MIDAS vs midas_u (纯 OLS 锚) | 1e-10 | [ ] | 风险 #4 |
| 7.2.2 | MIDAS-DL vs midas_r 夹具 (reltol=1e-12, maxit=10000, set.seed) | 1e-6~1e-8 | [ ] | 决策 24 |
| 7.2.3 | MIDAS-AR (+AR*) vs 夹具 | 1e-6~1e-8 | [ ] | |
| 7.2.4 | 集中化 NLS ≡ 联合 NLS 同一最优 | 1e-6 | [ ] | 决策 23 |
| 7.2.5 | 多起点 λ 网格×{递减,驼峰,均匀} 逃逸 (λ=0 平坦陷阱) | 流程 | [ ] | MD8 |
| 7.2.6 | 夹具 convergence + SSR 双判据 + hAh 三列全录 | 流程 | [ ] | |

---

## 8. 幻觉点逐项核查 (64 编号, probe/verify 脚本逐条落盘)

### 8.1 ARIMA/Granger (AR1-AR8, GR1-GR7)

| ID | 影响级 | 核查内容 | 容差 | 状态 | 脚本 |
|----|--------|----------|------|------|------|
| AR1 | 极高 | MA 两库同号 (1+θB) | 1e-10 | [ ] | verify_arima.R |
| AR2 | 极高 | n.cond=d+max(user,p) 与 q 无关 | 1e-10 | [ ] | |
| AR3 | 高 | loglik 口径统一高斯形式 | 1e-10 | [ ] | |
| AR4 | 高 | AIC 基于 T−d | 1e-10 | [ ] | |
| AR5 | 中 | drift=差分截距 | 1e-8 | [ ] | |
| AR6 | 极高 | method 配对对照 (禁混配) | - | [ ] | |
| AR7 | 高 | 多起始 (HR/CSS/随机) | - | [ ] | |
| AR8 | 高 | innovations 限定 (无缺失/无季节) | 1e-10 | [ ] | verify_arima.py |
| GR1 | 高 | F df 公式 (m=p, k_u=2p+1) | 1e-10 | [ ] | verify_granger.py |
| GR2 | 极高 | TY df=k | 1e-8 | [ ] | |
| GR3 | 中 | d_max 外部给定 | - | [ ] | |
| GR4 | 高 | 增广阶不进约束矩阵 | - | [ ] | |
| GR5 | 高 | 默认非稳健 (HAC 自建) | 1e-8 | [ ] | |
| GR6 | 极高 | 方向: 显式 (cause,effect) | 1e-10 | [ ] | |
| GR7 | 中 | I(1) 水平标准 F 失效 | - | [ ] | 集成 2 |

### 8.2 VAR/DY (V1-V13)

| ID | 影响级 | 核查内容 | 容差 | 状态 | 脚本 |
|----|--------|----------|------|------|------|
| V1 | 低 | PS 1998 = Econ. Letters | - | [ ] | 文献冻结 |
| V2 | 极高 | Cholesky 下三角 | 1e-12 | [ ] | verify_var.py |
| V3 | 极高 | Θ[i,j] 行响应列冲击 | 1e-12 | [ ] | |
| V4 | 极高 | IC 用 ML Σ (÷T) | 1e-10 | [ ] | |
| V5 | 高 | select_order 同样本 offset | 1e-10 | [ ] | |
| V6 | 高 | FPE 指数 K, n* 含 det | 1e-10 | [ ] | |
| V7 | 极高 | GFEVD 行归一化 (行和=1) | 1e-8 | [ ] | verify_gfevd.R |
| V8 | 极高 | 双系数框架 (DY σ_jj⁻¹ / PS σ_ii⁻¹) | 1e-8 | [ ] | |
| V9 | 高 | 稳定性双输出 (严格 <1) | 1e-12 | [ ] | |
| V10 | 中 | H 可配置 (非数学常数) | - | [ ] | |
| V11 | 低 | GIRF 隐含假设注明 (Kim 2013) | - | [ ] | 文档 |
| V12 | 高 | 不稳定 VAR 拦截 FEVD | - | [ ] | |
| V13 | 中 | bootstrap≠δ法, 容差仅点估计 | - | [ ] | |

### 8.3 协整 (CI1-CI12)

| ID | 影响级 | 核查内容 | 容差 | 状态 | 脚本 |
|----|--------|----------|------|------|------|
| CI1 | 极高 | 1994 响应面按 N 索引 | 精确/1e-12 | [ ] | .inc + verify_eg.py |
| CI2 | 高 | p/cv 不同源 (#4138) | - | [ ] | |
| CI3 | 中 | EG 方向依赖 (双方向输出) | 1e-10 | [ ] | |
| CI4 | 高 | 3 vs 5 情形 API 边界 | - | [ ] | |
| CI5 | 高 | statsmodels 一套表 (MHM96) | 精确 | [ ] | verify_johansen_diff.R |
| CI6 | 极高 | 迹公式+λ降序+有效T | 1e-10 | [ ] | verify_johansen.py |
| CI7 | 极高 | transitory 对照 urca | 1e-8 | [ ] | |
| CI8 | 极高 | β 投影空间对照 | 1e-10 | [ ] | verify_vecm.py |
| CI9 | 中 | Π 符号 (αᵢ<0 拉回) | - | [ ] | |
| CI10 | 高 | EM2002 查表 | 精确 | [ ] | .inc |
| CI11 | 低 | 文献出处 (JEDC/JAE) | - | [ ] | 文献冻结 |
| CI12 | 高 | Pu/Pz 双实现 | 1e-8 | [ ] | verify_po.R |

### 8.4 MIDAS (MD1-MD8)

| ID | 影响级 | 核查内容 | 容差 | 状态 | 脚本 |
|----|--------|----------|------|------|------|
| MD1 | 极高 | nealmon i=1 起 | 1e-12 | [ ] | verify_midas.R |
| MD2 | 高 | δ 独立线性参数 | 1e-12 | [ ] | |
| MD3 | 极高 | lag0=期末对齐 | 精确 | [ ] | |
| MD4 | 低 | U-MIDAS 出处 (F-M-S 2011/2015) | - | [ ] | 文献冻结 |
| MD5 | 低 | 综述 GSV 2007 | - | [ ] | |
| MD6 | 低 | K-Z 2012 | - | [ ] | |
| MD7 | 极高 | log-sum-exp 防溢出 | <1e-14 | [ ] | test_midas_weights |
| MD8 | 高 | 多起点逃逸 | - | [ ] | |

### 8.5 回填 (NP1-NP6, ZA1-ZA5, GM1-GM5)

| ID | 影响级 | 核查内容 | 容差 | 状态 | 脚本 |
|----|--------|----------|------|------|------|
| NP1 | 低 | NP 2001 = Econometrica | - | [ ] | 文献冻结 |
| NP2 | 极高 | τ_T 4 源公式 | 1e-10 | [ ] | verify_np_stata.py |
| NP3 | 极高 | AR 谱对差分拟合 | - | [ ] | |
| NP4 | 极高 | MPT 分情形 (+7/+14.5) | 1e-12 | [ ] | |
| NP5 | 极高 | 方向+Table 1 | 精确 | [ ] | np_table1.inc |
| NP6 | 中 | Stata 非基准 (仅 MAIC 列) | 1e-10 | [ ] | |
| ZA1 | 高 | 三库 lag 策略 (双模式) | 1e-10/1e-8 | [ ] | verify_za.py/.R |
| ZA2 | 高 | trim 参数化 | 1e-8 | [ ] | |
| ZA3 | 极高 | DT 断点后重新计时 | 精确 | [ ] | |
| ZA4 | 高 | 双临界值表 (MC c 1%=−5.27644) | 精确 | [ ] | .inc |
| ZA5 | 极高 | 取最负 min | 1e-12 | [ ] | |
| GM1 | 极高 | archpow 1=σ/2=σ² | - | [ ] | verify_gm.R |
| GM2 | 中 | arch 8.0 才有 ARCHInMean | - | [ ] | 版本 guard |
| GM3 | 高 | log 变体有基准 (form='log') | 1e-8 | [ ] | verify_gm.py |
| GM4 | 高 | λ sandwich SE (BW) | 1e-10 | [ ] | |
| GM5 | 极高 | 均值方差耦合递归 | 1e-10 | [ ] | |

**scope 外推迟项** (保持 [ ], 不计入统计): SARIMA / wild bootstrap / SVAR / BVAR / TVP-VAR / ARDL-PSS / midas 扩展族 / DCC / Kalman / 长记忆族 (开放问题 d-h, v1.8+)

---

## 9. 端到端集成测试 (spec §7, 6 场景)

| # | 场景 | 验证点 | 状态 |
|---|------|--------|------|
| 9.1 | 单位根诊断全链 (ADF/DF-GLS/NP/ZA → 差分 → ARIMA → LB) | M0+M1+7B 复用; 断点→分段平稳 | [ ] |
| 9.2 | Granger 因果链 (差分 F vs 水平 TY 对比) | GR7 失效场景 | [ ] |
| 9.3 | VAR→DY 溢出全链 (IC→稳定→IRF/FEVD→DY+滚动) | M2 全链; V12 拦截 | [ ] |
| 9.4 | 协整→VECM (rank→ECT 显著性→β 投影) | M3 全链 | [ ] |
| 9.5 | MIDAS 混频预测 (DL vs U-MIDAS, MZ/DM 复用 7B) | M4; MD3 期初起窗 | [ ] |
| 9.6 | GARCH-M 风险溢价 (三变体→λ 显著性→vs 无 M) | M0; GM4/GM5 | [ ] |

---

## 10. ADR 对齐验证

### 10.1 ADR-019 (26+3 项决策)

| # | 决策 | 实施一致性 | 状态 |
|---|------|-----------|------|
| 10.1.1 | D1 ARIMA 仅 CSS+CSS-ML+innovations (Kalman v1.8) | [ ] | |
| 10.1.2 | D2 (1+θB) 正号参数化 | [ ] | |
| 10.1.3 | D3 d>0 无均值 + drift 选项 | [ ] | |
| 10.1.4 | D4 Granger 主基准 statsmodels 4 统计量 | [ ] | |
| 10.1.5 | D5 TY 仅增广 Wald, d_max 用户传入 | [ ] | |
| 10.1.6 | D6 HAC-Wald 复用 NW | [ ] | |
| 10.1.7 | D7 ARIMA 多起始点 | [ ] | |
| 10.1.8 | D8 不做 SARIMA/wild bootstrap | [ ] | |
| 10.1.9 | D9 VAR 逐方程 OLS (method='ols') | [ ] | |
| 10.1.10 | D10 IC Lütkepohl 约定 + offset | [ ] | |
| 10.1.11 | D11 Cholesky LLT 下三角 + P 注入/重排 | [ ] | |
| 10.1.12 | D12 FEVD 双轨, DY 基于后者 | [ ] | |
| 10.1.13 | D13 is_stable 双输出 | [ ] | |
| 10.1.14 | D14 IRF 带 block bootstrap | [ ] | |
| 10.1.15 | D15 GFEVD 自实现 (R Spillover 主基准) | [ ] | |
| 10.1.16 | D16 不做 SVAR/BVAR/TVP-VAR | [ ] | |
| 10.1.17 | D17 EG + coint_johansen 等价 API (3 情形) | [ ] | |
| 10.1.18 | D18 EG 临界值 1994 响应面 (分文件) | [ ] | |
| 10.1.19 | D19 Johansen 主录 OL1992 + 双库 diff 前置 | [ ] | |
| 10.1.20 | D20 PO 纳入 (Pu/Pz, Pz 优先) | [ ] | |
| 10.1.21 | D21 VECM β 双归一 + EM2002 | [ ] | |
| 10.1.22 | D22 MIDAS scope (DL/AR/U-MIDAS + 5 权重) | [ ] | |
| 10.1.23 | D23 集中化 NLS (SLSQP 外层 + 内层解析) | [ ] | |
| 10.1.24 | D24 midasr 0.9 唯一主基准 (夹具固化) | [ ] | |
| 10.1.25 | D25 log-sum-exp 防溢出 | [ ] | |
| 10.1.26 | D26 MIDAS 扩展 scope 外 | [ ] | |
| 10.1.27 | 回填 NP (ERS 去势复用 + 四统计量 + NP Table 1) | [ ] | |
| 10.1.28 | 回填 ZA (三模型 + 双模式 + trim) | [ ] | |
| 10.1.29 | 回填 GM (三变体双锚 + fix 三步法) | [ ] | |

### 10.2 ADR-017 / ADR-013 / ADR-018 / R 门禁流程

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 10.2.1 | 命名空间: arima/var/cointegration/midas 落 `cpphub::v1::timeseries::*` | [ ] | ADR-017 |
| 10.2.2 | 头文件 #include 在 namespace 外 | [ ] | project_memory 教训 |
| 10.2.3 | 新 CMake target `cpphub_timeseries_mat` 链接 Eigen3 (仅 var/cointegration) | [ ] | ADR-013 |
| 10.2.4 | 单变量模块 (M0/M1/M4) 零 Eigen include (grep 断言) | [ ] | C1 |
| 10.2.5 | 单变量模块无 var/cointegration 反向依赖 | [ ] | §8.2 |
| 10.2.6 | SLSQP 复用不修改 (12/12 无退化) | [ ] | ADR-018 |
| 10.2.7 | C ABI 新导出 `cpphub_v1_7_*` 前缀 (如启用) | [ ] | C3 |
| 10.2.8 | 实施期公式冲突回溯 ADR-019 修订 (无静默改公式) | [ ] | R 门禁条款 |
| 10.2.9 | 对照禁令执行 (Julia 常数 MPT / arch≠MAIC-PQ) | [ ] | H1/H2 裁决 |

---

## 11. Scope 边界验证

### 11.1 已实现 (scope 内)

| # | 检查项 | 状态 |
|---|--------|------|
| 11.1.1 | NP 四统计量 + MAIC/MBIC/seq-t | [ ] |
| 11.1.2 | ZA 三模型 + 双模式 | [ ] |
| 11.1.3 | GARCH-M 三变体 | [ ] |
| 11.1.4 | ARIMA CSS/CSS-ML/innovations | [ ] |
| 11.1.5 | Granger 4 统计量 + TY + HAC | [ ] |
| 11.1.6 | VAR + IC + IRF + FEVD 双轨 + DY | [ ] |
| 11.1.7 | EG + Johansen + VECM + PO | [ ] |
| 11.1.8 | MIDAS-DL/AR/U-MIDAS + 5 权重族 | [ ] |

### 11.2 未实现 (scope 外, 确认未越界)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 11.2.1 | SARIMA 未实现 | [ ] | v1.8+ |
| 11.2.2 | Granger wild bootstrap 未实现 | [ ] | v1.8+ |
| 11.2.3 | SVAR/BVAR/TVP-VAR 未实现 | [ ] | v1.8+ |
| 11.2.4 | ARDL/PSS 未实现 | [ ] | v1.8+ |
| 11.2.5 | midas_nlpr/sp/qr/imidas_r/amweights 未实现 | [ ] | v1.8+ |
| 11.2.6 | gompertzp/nakagamip/lcauchyp/genexp 未实现 | [ ] | v1.8+ (实名带 p 后缀) |
| 11.2.7 | DCC/CCC 未实现 | [ ] | v1.8+/v1.9 |
| 11.2.8 | Kalman/statespace 未实现 | [ ] | v1.8+ |
| 11.2.9 | 长记忆族 (APARCH/FIGARCH/IGARCH/ARFIMA/HYGARCH) 未实现 | [ ] | v1.8+ |

---

## 12. 代码规范验证

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 12.1 | 头文件 #include 在 namespace 外 | [ ] | C2065 教训 |
| 12.2 | 无参数名遮蔽函数名 | [ ] | v1.4.2 教训 |
| 12.3 | 接口签名与 spec §2-§6 一致 | [ ] | |
| 12.4 | Result 结构体含 spec 全部字段 | [ ] | |
| 12.5 | 异常处理 (invalid_argument/runtime_error) | [ ] | |
| 12.6 | 无 magic number (constexpr + static_assert) | [ ] | .inc 六件套 |
| 12.7 | 排幻觉点注释标注 (如 `// NP4: MPT 分情形 +14.5`) | [ ] | |
| 12.8 | 双框架/双模式 API 注明适用域 (GFEVD 框架/ZA 模式/PQ vs NP) | [ ] | 7C 特有 |
| 12.9 | 源码 UTF-8 (MSVC /utf-8 沿用) | [ ] | C4819 教训 |

---

## 13. 复用验证 (不重复实现)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 13.1 | NP 复用 DF-GLS 去势 (同变换不重写) | [ ] | spec §2.1 Step1 |
| 13.2 | ZA 复用 ADF 引擎滞后框架 | [ ] | |
| 13.3 | GM 复用 GARCH backcast/似然/SLSQP/sandwich | [ ] | 仅耦合递归新增 |
| 13.4 | ARIMA 外层复用 SLSQP (ADR-018) | [ ] | C5 |
| 13.5 | Granger HAC 复用 v1.5 hac/NW | [ ] | 决策 6 |
| 13.6 | IRF 置信带复用 block bootstrap (v1.5) | [ ] | 决策 14 |
| 13.7 | PO 长期方差复用 long_run_variance (7B) | [ ] | |
| 13.8 | MIDAS 内层复用 OLS/QR (v1.5) | [ ] | 决策 23 |
| 13.9 | 诊断复用 Phase 7A (JB/LB/multiple_test_correction) | [ ] | |
| 13.10 | 复用偏差 (如 7B Issue #2 型) 逐项归因记录 | [ ] | 豁免须理由 |

---

## 14. 文档对齐验证

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 14.1 | ADR-019 26+3 项全部实施且无静默偏离 | [ ] | 冲突回溯记录 |
| 14.2 | ADR-017 命名空间 + ADR-013 Eigen 隔离遵守 | [ ] | |
| 14.3 | C1-C7 兼容性约束逐项满足 | [ ] | 调研第四部分 |
| 14.4 | DEVELOPMENT_LOG.md 逐里程碑更新 | [ ] | |
| 14.5 | PHASE7C_FINAL_ACCEPTANCE.md 编写 | [ ] | 收尾 |
| 14.6 | spec 幻觉点注释保留到代码 | [ ] | |
| 14.7 | verify 脚本可追溯 (源码打印+手算+baseline) | [ ] | 附录 B 模式 |
| 14.8 | R 脚本含 `.libPaths()` 显式加载 | [ ] | Rscript 教训 |
| 14.9 | 开放问题 [待定] 项处置记录 (a-c 已裁决生效, H1/H2 假设区) | [ ] | spec §13 v1.1 |

---

## 15. 性能验证

| # | 检查项 | 基准 | 状态 | 备注 |
|---|--------|------|------|------|
| 15.1 | NP T=1000 (逐 k MAIC 搜索) | < 1 sec | [ ] | |
| 15.2 | ZA T=500 断点网格搜索 | < 5 sec | [ ] | |
| 15.3 | GM T=5000 三变体 | < 10 sec | [ ] | |
| 15.4 | ARIMA T=1000 CSS-ML 多起始 | < 10 sec | [ ] | |
| 15.5 | VAR K=5, T=500 + IC 扫描 | < 5 sec | [ ] | |
| 15.6 | Johansen N=5, T=500 | < 2 sec | [ ] | |
| 15.7 | MIDAS NLS T=250, m=22 | < 10 sec | [ ] | |
| 15.8 | 全量 ctest (~2470) | < 45 min | [ ] | 7B 基线 614s + 增量 |

> **溯源修复 (v1.2)**: 本节 8 项预算已冻结于 spec §8.8 性能预算 (此前为 checklist 编写期自行推导, 无 spec 源头 — 现恢复 "验收点 ← spec" 单向追溯); 另 spec §1.4 接口总则 + §8.4-8.7 工程约束 (构建矩阵/零破坏/精度/线程) 同步生效, 对应本 checklist §2/§10/§12 验收项。

---

## 16. 风险项验收 (spec §11 全量)

| # | 风险 | 缓解落实 | 状态 |
|---|------|---------|------|
| 16.1 | Johansen OL1992 录入 (无机器可读版) | urca 转录 static_assert + 双库 diff 前置冻结 | [ ] |
| 16.2 | NP 无权威开源基准 | 原文>Stata MAIC>EViews; MZt≡MZα×MSB 自检 | [ ] |
| 16.3 | ARIMA 多局部极值假失败 | 多起始 + method 配对矩阵 | [ ] |
| 16.4 | midasr 唯一主基准 BFGS 容差 | 夹具收紧 + midas_u 1e-10 锚 + 双判据 | [ ] |
| 16.5 | ZA 三库不可逐位对齐 | 双模式 API + trim 参数化 | [ ] |
| 16.6 | GM vs rugarch solver 差 | 1e-4 + fix() 三步法 + rescale 统一关 | [ ] |
| 16.7 | PQ 2007 vs NP 2001 微差 | NP 2001 为准 (H2 裁决, PQ v1.8 Gretl 对照) | [ ] |

---

## 17. 最终签字

### 17.1 验收统计 (编写期冻结验收点框架, 实施后填写)

| 维度 | 总项数 | 通过 | 未通过 | 通过率 |
|------|--------|------|--------|--------|
| 1. 交付物完整性 | 66 | | | |
| 2. 编译与跨平台 | 15 | | | |
| 3. M0 数值基准 | 22 | | | |
| 4. M1 数值基准 | 14 | | | |
| 5. M2 数值基准 | 18 | | | |
| 6. M3 数值基准 | 16 | | | |
| 7. M4 数值基准 | 12 | | | |
| 8. 幻觉点核查 | 64 (+10 scope 外不计) | | | |
| 9. 端到端集成 | 6 | | | |
| 10. ADR 对齐 | 38 | | | |
| 11. Scope 边界 | 17 | | | |
| 12. 代码规范 | 9 | | | |
| 13. 复用验证 | 10 | | | |
| 14. 文档对齐 | 9 | | | |
| 15. 性能验证 | 8 | | | |
| 16. 风险项 | 7 | | | |
| **总计** | **331** | | | |

### 17.2 验收结论

- [ ] **全部通过** — Phase 7C 验收完成, 可进入 v1.7 发布
- [ ] **有条件通过** — 存在未通过项, 附 issue 列表, 限期修复
- [ ] **未通过** — 存在极高严重性未通过项, 需返工

### 17.3 签字

| 角色 | 姓名 | 日期 | 签字 |
|------|------|------|------|
| 实施者 | Scott (w/ Claude GLM-5.3) | 待填 | |
| 审计者 | Scott (self-review, 326 项逐项审计) | 待填 | |
| 架构组 | Scott (solo developer 工作流) | 待填 | |

---

## 附录 A: 未通过项 Issue 模板

```
Issue #___
关联检查项: §__._
问题描述:
复现步骤:
预期结果:
实际结果:
严重性: [极高/高/中/低]
修复方案:
修复状态: [ ] 待修复 [ ] 已修复 [ ] 已验证
```

## 附录 B: verify 脚本核查模式 (沿用 7B, 参照 v1.4.0/1.4.1)

每个 verify 脚本必须包含:
1. **对照库源码打印**: `inspect.getsource()` (Python) / `print()` 关键行 (R: 直接引用 tarball 源码行号)
2. **小样本手算**: 构造 T=5~10 简单序列手动计算预期值
3. **基准数值**: 调用基准库 API 生成 (R/midasr/urca/rugarch/vars/Spillover; Stata dfgls r(results) 抄录 CSV)
4. **C++ 输出对照**: 逐数值对比
5. **容差判定**: 相对误差 ≤ spec 标注容差 (分层: 1e-10/1e-8/1e-6/1e-4/精确)
6. **夹具固化**: midasr 输出写 CSV (含 convergence 码 + SSR + set.seed 信息)

## 附录 C: 里程碑验收节奏 (建议)

| 里程碑 | 验收节点 | 覆盖节 |
|--------|---------|--------|
| M0 | 逐项审计 §3 + §8.5 + §1.1 (1.1.1-1.1.4) | 先行 |
| M1 | §4 + §8.1 | 与 M4 并行 |
| M4 | §7 + §8.4 | |
| M2 | §5 + §8.2 + DY 窗口裁决项 (§5.3.2-5.3.3) | |
| M3 | §6 + §8.3 + 双库 diff 前置 (§6.2.1) | |
| 收尾 | §2/§9-§17 全量 + FINAL_ACCEPTANCE | 三平台 |
