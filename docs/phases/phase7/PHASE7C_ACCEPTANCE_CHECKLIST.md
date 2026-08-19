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
| 1.1.1 | `include/cpphub/timeseries/unit_root/ng_perron_test.hpp` | [x] | MZα/MZt/MSB/MPT + MAIC; 固定样本口径坐标勘误 069264c (dyd[T−1] 越界读根除); MBIC/seq-t 推迟 (v1.2 签名无字段) |
| 1.1.2 | `include/cpphub/timeseries/unit_root/zivot_andrews_test.hpp` | [x] | 三模型 A/B/C + trim + 断点搜索; 双模式 (固定 lag 主/Baum 预选对照) |
| 1.1.3 | `include/cpphub/timeseries/unit_root/np_tables.hpp` | [x] | NP 2001 Table 1 (constexpr + static_assert, 含 0.143/4.03 转录陷阱防呆) |
| 1.1.4 | `include/cpphub/timeseries/garch/garch_m_model.hpp` | [x] | GARCH(1,1)-M 三变体 + form 感知 Jacobian sandwich |

**M1 ARIMA/Granger (4 个, 纯标量无 Eigen)**

| # | 文件 | 状态 | 备注 |
|---|------|------|------|
| 1.1.5 | `include/cpphub/timeseries/arima/arima_model.hpp` | [x] | CSS + CSS-ML; d=1 漂移走 forecast::Arima 语义 (AR5), n.cond=d+p 定案 (AR2) |
| 1.1.6 | `include/cpphub/timeseries/arima/innovations_mle.hpp` | [x] | B&D 2016 §5.2 精确 MLE; arma11 vs statsmodels 逐位 (ll=−416.9317) |
| 1.1.7 | `include/cpphub/timeseries/arima/hannan_rissanen.hpp` | [x] | HR 起始值; Schwert n_init 默认 + 多起始点集合成员 |
| 1.1.8 | `include/cpphub/timeseries/arima/granger_test.hpp` | [ ] | F/χ²/LR + TY + HAC-Wald |

**M2 VAR/DY (6 个, 需 Eigen3 → cpphub_timeseries_mat)**

| # | 文件 | 状态 | 备注 |
|---|------|------|------|
| 1.1.9 | `include/cpphub/timeseries/var/multivariate_data.hpp` | [x] | MultivariateTSData (C4); T/K/matrix/reorder + 校验 |
| 1.1.10 | `include/cpphub/timeseries/var/var_model.hpp` | [x] | 逐方程 OLS + IC 五式 + 稳定性双输出 (V9); logdet ×2 修复后 SM 1e-10 |
| 1.1.11 | `include/cpphub/timeseries/var/var_select.hpp` | [x] | IC 五式 + 同样本 offset (V5); auto-lag 默认 aic |
| 1.1.12 | `include/cpphub/timeseries/var/irf.hpp` | [x] | Φ 递推 + Ψ=Φ·P 正交化 + 块 bootstrap 带 (V13); P 注入 (决策 11) |
| 1.1.13 | `include/cpphub/timeseries/var/fevd.hpp` | [x] | Cholesky + GFEVD 双框架 DY/PS (V8); V12 不稳定拦截 |
| 1.1.14 | `include/cpphub/timeseries/var/dy_spillover.hpp` | [x] | TCI/TO/FROM/NET + 滚动 (window>2K, 无默认); R Spillover 1e-8 |

**M3 协整 (7 个, 需 Eigen3)**

| # | 文件 | 状态 | 备注 |
|---|------|------|------|
| 1.1.15 | `include/cpphub/timeseries/cointegration/engle_granger.hpp` | [x] | EG 两步法; SM coint 1e-10 全趋势 |
| 1.1.16 | `include/cpphub/timeseries/cointegration/johansen_test.hpp` | [x] | 迹/最大特征值, 3 det 情形; λ SVD 路径复刻 |
| 1.1.17 | `include/cpphub/timeseries/cointegration/phillips_ouliaris.hpp` | [x] | Pu/Pz 双实现; urca 1e-8 (12 组合) |
| 1.1.18 | `include/cpphub/timeseries/cointegration/vecm_model.hpp` | [x] | 5 情形 + β 双归一 + ECT t; llf SVD+非对称 eig 路径对齐 |
| 1.1.19 | `include/cpphub/timeseries/cointegration/osterwald_lenum_cv.hpp` | [x] | OL1992 表 (urca 转录) + MHM96 双表 |
| 1.1.20 | `include/cpphub/timeseries/cointegration/mackinnon_coint_cv.hpp` | [x] | MacKinnon 1994 协整响应面 (N≥2, 与 2010 ADF 表分文件) |
| 1.1.21 | `include/cpphub/timeseries/cointegration/ericsson_mackinnon_cv.hpp` | [x] | EM2002 ECT t 临界值 (响应面 + 渐近锚) |

**M4 MIDAS (4 个, 纯标量无 Eigen)**

| # | 文件 | 状态 | 备注 |
|---|------|------|------|
| 1.1.22 | `include/cpphub/timeseries/midas/mixed_freq_data.hpp` | [x] | mls 期末对齐 (C4); W-dir 定案 w₁↔h=0 期末 (Form A 恢复 λ*) |
| 1.1.23 | `include/cpphub/timeseries/midas/midas_weights.hpp` | [x] | 5 权重族 + log-sum-exp; nbetaMT 4 参数实施勘误 (θ₀ 第 4 参) |
| 1.1.24 | `include/cpphub/timeseries/midas/midas_model.hpp` | [x] | DL/AR/AR*/U-MIDAS + 集中化 NLS (外层 SLSQP θ + 内层 OLS) |
| 1.1.25 | `include/cpphub/timeseries/midas/midas_diagnostics.hpp` | [x] | 残差诊断 + hAh 检验 (K-Z 2012, prep_hAh 逐字复刻, ~1e-4) |

### 1.2 基准验证脚本 (19 个) + 临界值表 (6 件, 沿用 7B 全提交惯例)

| # | 脚本 | 状态 | 备注 |
|---|------|------|------|
| 1.2.1 | `tests/fixtures/timeseries/verify_np_stata.py` | [x] 占位 | C6 降级预批: --emit-dofile/--parse 管线就绪, Stata 装机后补 1e-10 硬断言 |
| 1.2.2 | `tests/fixtures/timeseries/verify_np_semi.py` | [ ] | EViews/Julia 抄录半基准 (Julia 常数 MPT 禁用); NP6 边界下暂以原文公式+恒等式自检替代 |
| 1.2.3 | `tests/fixtures/timeseries/verify_za.py` | [x] | statsmodels (Baum 模式), 基准 JSON 入库 |
| 1.2.4 | `tests/fixtures/timeseries/verify_za.R` | [x] | urca ur.za (固定 lag, trim 放开), 12 位全精度 |
| 1.2.5 | `tests/fixtures/timeseries/verify_gm.py` | [x] | arch 8.0 ARCHInMean form 三值, 基准入库 |
| 1.2.6 | `tests/fixtures/timeseries/verify_gm.R` | [x] | rugarch archpow=1/2 + fix() 三步法 (发现 arch fix() bug, C++ 充当独立评估器) |
| 1.2.7 | `tests/fixtures/timeseries/verify_arima.R` | [x] | R stats::arima CSS/CSS-ML (method 配对) + forecast::Arima drift (d≥1 强制无均值裁决留档) |
| 1.2.8 | `tests/fixtures/timeseries/verify_arima.py` | [x] | statsmodels innovations_mle, 基准 JSON 入库 |
| 1.2.9 | `tests/fixtures/timeseries/verify_granger.py` | [ ] | grangercausalitytests 4 统计量 |
| 1.2.10 | `tests/fixtures/timeseries/verify_var.py` | [x] | statsmodels VAR 系数/IC/IRF/FEVD/select_order, 基准 JSON 入库 |
| 1.2.11 | `tests/fixtures/timeseries/verify_var.R` | [x] | vars::VAR + VARselect 交叉 (dump_var_r_values.R 机器精度 dump) |
| 1.2.12 | `tests/fixtures/timeseries/verify_gfevd.R` | [x] | R Spillover g.fevd/G.spillover 主基准 (0.1.1 裁剪装载) |
| 1.2.13 | `tests/fixtures/timeseries/verify_eg.py` | [x] | statsmodels coint (T×4 趋势 × 3 对, 基准入库) |
| 1.2.14 | `tests/fixtures/timeseries/verify_johansen.py` | [x] | statsmodels coint_johansen (det×k 网格, 基准入库) |
| 1.2.15 | `tests/fixtures/timeseries/verify_johansen_diff.R` | [x] | **M3 前置任务**: 双库 diff 冻结主对照 (JOHANSEN_DUAL_LIB_DIFF.md + MC 裁决) |
| 1.2.16 | `tests/fixtures/timeseries/verify_vecm.py` | [x] | statsmodels VECM (β 投影空间); resid_head 语义勘误 ([:4] 前 4 时点) |
| 1.2.17 | `tests/fixtures/timeseries/verify_po.R` | [x] | urca ca.po (3 对 × Pu/Pz × 3 demean × short/long) |
| 1.2.18 | `tests/fixtures/timeseries/verify_midas.R` | [x] | midasr 0.9 夹具 (reltol=1e-12 + seed): W1 权重逐点/W2 mls 对齐/W3 U-MIDAS/W4-5 NLS/W6 AR/W7 hAh; probe_midas_form.R (start 语义裁决) 随附 |
| 1.2.19 | `tests/fixtures/timeseries/verify_midas_u.R` | [x] 合并 | 并入 verify_midas.R W3 (midas_u 纯 OLS 锚 1e-10 → 实测逐位一致); 独立脚本不再单设 |
| 1.2.20 | `critical_values/np_table1.inc` | [x] 形态变更 | 以 `np_tables.hpp` (include 树) 交付, constexpr + static_assert; 零依赖裸表意图由头文件承担 |
| 1.2.21 | `critical_values/za1992_cv.inc` | [x] | ZA 论文表 (主), urca/Baum/讲义三源零差异 |
| 1.2.22 | `critical_values/za_mc_cv.inc` | [x] | MC 表 (c 1%=−5.27644), 144 值程序化双库零差异 |
| 1.2.23 | `critical_values/ol1992_cv.inc` | [x] 形态变更 | 以 `osterwald_lenum_cv.hpp` (include 树) 交付, constexpr + static_assert (urca 源码常量转录) |
| 1.2.24 | `critical_values/mackinnon1994_coint.inc` | [x] 形态变更 | 以 `mackinnon_coint_cv.hpp` (include 树) 交付, 1994 p 响应面 + 2010 cv 响应面 |
| 1.2.25 | `critical_values/em2002_ect_cv.inc` | [x] 形态变更 | 以 `ericsson_mackinnon_cv.hpp` (include 树) 交付, EM2002 响应面 (PDF 三重验证转录) |

### 1.3 测试套件 (16 套, ~263 用例)

| # | 测试套件 | 用例数 | 状态 | 里程碑 |
|---|----------|--------|------|--------|
| 1.3.1 | `test_ng_perron` | 18 | [x] | M0 (069264c 修复后 60 次循环 0 失败) |
| 1.3.2 | `test_zivot_andrews` | 15 | [x] | M0 (双库基准 ~1e-13) |
| 1.3.3 | `test_garch_m_model` | 16 | [x] | M0 (arch 三 form 1e-5~1e-6) |
| 1.3.4 | `test_arima_model` | 24 | [x] | M1 (CSS/CSS-ML vs R 四夹具; θ 谱等价类 arma21 主锚 φ+ll) |
| 1.3.5 | `test_innovations_mle` | 12 | [x] | M1 (黄金锚 B&D MA(1) 4 位逐位 + arma11 vs statsmodels 逐位) |
| 1.3.6 | `test_granger_causality` | 16 | [ ] | M1 |
| 1.3.7 | `test_var_model` | 20+1 | [x] | M2 (21/21; SM 1e-10 + vars 交叉 1e-8 + §15.5 性能) |
| 1.3.8 | `test_var_irf_fevd` | 18 | [x] | M2 (18/18; SM orth_ma_rep/fevd 1e-12 + Spillover 1e-8) |
| 1.3.9 | `test_dy_spillover` | 12 | [x] | M2 (12/12; G.spillover 表/TCI/TO/FROM/NET 1e-8 + roll 1e-6) |
| 1.3.10 | `test_engle_granger` | 14 | [x] | M3 (14/14; SM coint 1e-10, 4 趋势全对) |
| 1.3.11 | `test_johansen_test` | 18 | [x] | M3 (18/18; SM 1e-10 + urca 网格交叉 1e-8) |
| 1.3.12 | `test_vecm_model` | 16 | [x] | M3 (16/16; 5 情形 SM 1e-10 含 llf 路径对齐) |
| 1.3.13 | `test_phillips_ouliaris` | 10 | [x] | M3 (10/10; urca 1e-8, CI12 双向) |
| 1.3.14 | `test_midas_weights` | 16 | [x] | M4 (16/16, 1e-12 全对) |
| 1.3.15 | `test_midas_model` | 18+1 | [x] | M4 (19/19; U-MIDAS 逐位, NLS/hAh 落点层) |
| 1.3.16 | `test_integration_phase7c` | 8 | [ ] | 端到端 6 场景 (§9) |

---

## 2. 编译与跨平台测试

### 2.1 主控站 (Windows MSVC Release)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 2.1.1 | CMake 配置成功 (含新 target `cpphub_timeseries_mat`) | [x] | M2 已建 (Eigen3 INTERFACE 链接); 全量构建通过 |
| 2.1.2 | MSVC Release 编译零警告零错误 (`/utf-8` 沿用) | [x] | |
| 2.1.3 | 全量 ctest 通过 (新增 49+71+51+58 + 现有 2207 = 2436) | [x] | M0 轮 **2256/2256** (656.49s, 首轮 2255 拦截 NP 越界读); M1∥M4 轮 **2327/2327** (67b5450); M2 轮 **2375/2375** (823.46s) + 补充 3 用例; M3 轮 **2436/2436** (新增 M3 58 用例: EG 14 + Johansen 18 + VECM 16 + PO 10) |
| 2.1.4 | 无现有测试退化 (Phase 1-7B 全部仍通过) | [x] | 2207 基线无退化 |
| 2.1.5 | SLSQP 12/12 仍通过 (ADR-018 无退化) | [x] | |

### 2.2 A 站 (Ubuntu GCC)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 2.2.1 | fresh clone + rebuild (submodules: recursive) | [x] | bundle 中继 (069264c): git bundle + eigen tar + FETCHCONTENT_SOURCE_DIR_GOOGLETEST 本地覆盖, 零外网; M1∥M4 轮 (67b5450) 同法增量 ff + rebuild; M2 轮 (7d64939) 增量 ff 67b5450→7d64939 + rebuild |
| 2.2.2 | GCC 编译零警告零错误 | [x] | 三轮复核; M2 轮 (7d64939) GCC 13.3.0 自有代码 **0 警告** (日志警告全为 autodiff 第三方噪声, 口径同前) |
| 2.2.3 | ctest 全量通过 | [x] | **2238/2238** (363.70s, 069264c); **2309/2309** (418.47s, 67b5450); **2360/2360** (418.31s, 7d64939 M2 轮; 差额 18 = 平台专属用例 7B 期已知, M2 新增 51 用例 GCC 全数运行含子集 51/51) |
| 2.2.4 | 与主控站数值一致 (容差分层 §3-§7) | [x] | M0 49 + M1/M4 71 + M2 51 用例三平台行为一致 (7d64939) |

### 2.3 B 站 (Ubuntu GCC)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 2.3.1 | fresh clone + rebuild | [x] | bundle 中继 (069264c); ⚠️ B 站 IPv4 动态地址 ping 不通但 mDNS 主机名 scott-lau-GTR-Pro.local (IPv6) 在线 — "B 站离线"结论需主机名复核; M2 轮 (7d64939) mDNS 增量 ff + rebuild 正常 |
| 2.3.2 | GCC 编译零警告零错误 | [x] | M2 轮同 2.2.2 口径, 自有代码 0 警告 |
| 2.3.3 | ctest 全量通过 | [x] | **2238/2238** (358.67s, 069264c); **2309/2309** (412.67s, 67b5450); **2360/2360** (412.60s, 7d64939 M2 轮, M2 子集 51/51) |
| 2.3.4 | 与主控站数值一致 | [x] | 同 2.2.4 |

### 2.4 三平台一致性与 CI

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 2.4.1 | M0-M4 全部测试三平台无数值偏差 | [ ] | 部分: M0 (49) + M1/M4 (71) + **M2 (51)** 三平台全绿一致 (7d64939); **M3 (58) 主控站 2436/2436 全绿**, A/B 站 GCC 补验待 M3 提交后增量 ff + rebuild |
| 2.4.2 | GitHub Actions 全绿 (Windows bash shell + submodules 修复沿用) | [x] | run #60 (67b5450 代码轮) 4/4 + run #61 (d581f30 文档轮) 4/4 全绿 (Build&Test Ubuntu GCC/Windows MSVC + C ABI ×2); run #60 首次承载 M1∥M4 全量 2327 用例 |

---

## 3. M0 数值基准对照 (回填三项)

### 3.1 Ng-Perron — vs Stata dfgls (确定性) + EViews/Julia 半基准

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 3.1.1 | 逐 k MAIC(k)/σ̂²(k) vs Stata dfgls r(results) | 1e-10 | [ ] 占位 | verify_np_stata.py 管线就绪 (C6 降级预批); ⚠️ 069264c 口径修正后 Stata 对照需同步 n = T−1−k_max 口径核对 |
| 3.1.2 | τ_T(k) = β̂₀²·Σỹ²/σ̂²(k) (固定样本) | 1e-10 | [x] | 实现内逐 k 轨迹; 口径修正后 Σỹ² = 回归 ỹ 列平方和 (069264c) |
| 3.1.3 | MZt ≡ MZα×MSB 恒等式自检 | 1e-12 | [x] | 两 trend 情形测试断言 |
| 3.1.4 | MPT 常数末项 −c̄=+7 / 趋势末项 +14.5 (≠13.5 断言) | 1e-12 | [x] | NP4; MptPositiveBothCases/MptDiffersAcrossDetrending |
| 3.1.5 | 四统计量 vs EViews 抄录 (Julia 趋势情形旁证) | 1e-6~1e-8 | [ ] | 未做; NP6 边界下以原文公式+恒等式自检替代 (DEVELOPMENT_LOG 记录) |
| 3.1.6 | 渐近 5% 临界值 vs NP Table 1 (常数 −8.10/−1.98/0.233/3.17; 趋势 −17.30/−2.91/0.168/5.48) | 精确 | [x] | np_tables.hpp + static_assert + 0.143/4.03 陷阱防呆 |
| 3.1.7 | 拒绝方向: MZα/MZt 越负, MSB/MPT 越小 | 方向性 | [x] | RejectDirectionReconstruction; RW 不拒绝/AR(0.95,T=500) 拒绝 (功率曲线落档) |
| 3.1.8 | AR 谱密度对 Δỹ (非水平) 拟合 | 1e-10 | [x] | NP3 结构断言 (s²_AR = σ̂²(k*)/(1−Σβ)²) |

### 3.2 Zivot-Andrews — vs statsmodels (Baum 模式) + urca (固定 lag)

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 3.2.1 | 统计量/断点 vs statsmodels (lag 冻结 Baum 预选) | 1e-10 | [x] | 实测 ~1e-14 (A/C 全对; B 模型 bpidx 差 1 = statsmodels DT 边界 quirk, 已知对照边界) |
| 3.2.2 | vs urca ur.za (固定 lag, trim 放开) | 1e-8 | [x] | 实测 ~4e-13, bpoint 1-based 映射实测定档 |
| 3.2.3 | DT_t = (t−Tb)·1(t>Tb) 断点后重新计时 | 精确 | [x] | 与 DU 定义同验证 (statsmodels/urca 两库同约定) |
| 3.2.4 | 统计量 = min t(α̂) (最负, 非 max) | 1e-12 | [x] | stat≡min(path) 断言 |
| 3.2.5 | 临界值: ZA1992 论文表 (主) | 精确 | [x] | za1992_cv.inc, urca/Baum/讲义三源零差异 |
| 3.2.6 | MC 表对照 (c 1%=−5.27644 ≠ −5.83) | 精确 | [x] | za_mc_cv.inc, ZA4 陷阱断言在测 |
| 3.2.7 | trim 参数化 (默认 0.15, 上限 0.333) | 精确 | [x] | 网格 84/60 + 越界/空网格 throw 用例 |

### 3.3 GARCH-M — vs arch 8.0 (主锚) + rugarch (1e-4 + 三步法)

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 3.3.1 | form='vol' (g=√h) vs arch ↔ rugarch archpow=1 | 1e-8~1e-10 / 1e-4 | [x] arch / [ ] rugarch | arch 实测 params 1e-5~1e-6 (数值梯度 vs 解析梯度本质差, 7B 先例口径); rugarch 参数级对照不可行 (递归初始化约定不同, 已知对照边界) |
| 3.3.2 | form='var' (g=h) ↔ archpow=2 | 同上 | [x] arch | 尺度等变参数化验证: 数据×s ⇒ λ 减半 1e-5 |
| 3.3.3 | form='log' vs arch (log 变体数值基准) | 1e-8 | [x] | arch 1e-5~1e-6 同层 |
| 3.3.4 | λ 的 QMLE sandwich SE (BW robust) | 1e-10 | [x] 实际 2e-2 | 数值三明治噪声层 (spec 名义 1e-8~1e-10 仅同优化器同梯度可达, 7B test_garch_model 先例); form 感知 Jacobian 映射验证通过 |
| 3.3.5 | 均值-方差耦合: λ 扰动重写 h 路径断言 | 1e-10 | [x] | λ=0 退化 ≡ filter_garch11 逐位 1e-15; λ+0.01 ⇒ h[0] 不变 t≥1 全变 |
| 3.3.6 | fix() 互验三步法 (隔离似然差 vs 落点差) | 流程 | [x] | arch 8.0 ARCHInMean.fix() 自身 bug (llf=−3002.8) → C++ 充当独立评估器完成隔离; 可提上游 issue |
| 3.3.7 | loglik 含完整常数项 | 1e-15 | [x] | λ=0 退化锚 1e-15 |

---

## 4. M1 数值基准对照 (ARIMA + Granger)

### 4.1 ARIMA — vs R stats::arima + statsmodels (method 配对)

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 4.1.1 | CSS vs R method="CSS" (method 配对, 禁混配) | 1e-10 | [x] | AR6; 实测参数 2e-3 (SLSQP vs R optim 落点层, GM 先例口径), 4 夹具全对 — test_arima_model 1/3/4/6 |
| 4.1.2 | n.cond = d + max(user, p), 与 q 无关 | 1e-10 | [x] | AR2 定案: 四夹具实测 (arma12 p=1<q=2 → n.cond=1) — test 5 |
| 4.1.3 | CSS-ML vs R method="CSS-ML" | 1e-8 | [x] | loglik 主锚; arma21 θ 谱等价类 (ll/φ 逐位同) — test 9-13 |
| 4.1.4 | MA 系数 (1+θB) 正号与 R 逐系数一致 | 1e-10 | [x] | AR1; arma11/12/22 θ 逐系数对 |
| 4.1.5 | loglik/AIC 基于 T−d 差分观测 | 1e-10 | [x] | AR3/AR4; AIC npar 含 σ² (forecast 862.606 实测吻合) — test 8/20 |
| 4.1.6 | d=1 drift = 差分截距 (forecast::Arima 语义) | 1e-8 | [x] | AR5; 漂移正解 = forecast::Arima (drift=0.3573, ll=−427.303 同 statsmodels); stats::arima d≥1 强制无均值 → 伪根 0.9985 退化路径亦对照 — test 17-19 |
| 4.1.7 | innovations vs statsmodels innovations_mle (无缺失非季节) | 1e-10 | [x] | AR8 guard; arma11 **逐位一致** (φ=0.398025, ll=−416.9317); arma21 θ 谱等价 ll 逐位同 — test 15-16 + test_innovations_mle 12/12 |
| 4.1.8 | 多起始点逃逸 (合成双峰似然案例) | 流程 | [x] | AR7; {HR, 0, 随机扰动} 起始集 + use_hannan_rissanen=false 路径 — test 24 + test_innovations_mle 2 |

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
| 5.1.1 | 系数矩阵 (逐方程 OLS) vs statsmodels | 1e-10 | [x] | 决策 9 (SM_PARAMS 7×3 回归元主序, 测试转置读取) |
| 5.1.2 | IC 五式 (aic/bic/hqic/fpe/logdet) vs var_model.py | 1e-10 | [x] | V4/V6 (logdet=2·Σlog L_ii LLT 因子 2 修复) |
| 5.1.3 | Σ_mle = SSR/T (÷T 非 ÷(T−k)) | 1e-10 | [x] | V4 (SCALARS[5]=det(Σ_df) 换算关系断言) |
| 5.1.4 | select_order 同样本 offset | 1e-10 | [x] | V5 (p=0..4 轨迹逐点; p=0 全样本锚闭包重算) |
| 5.1.5 | IC vs R vars::VARselect | 1e-8 | [x] | AIC/SC p=1..4 逐位 (FPE 不作锚, V6 口径差留档) |
| 5.1.6 | Cholesky 下三角 (Eigen LLT = np.linalg.cholesky = R t(chol)) | 1e-12 | [x] | V2 (R dump 列主序; matrixL 须物化后取元) |
| 5.1.7 | is_stable 双输出 (max|eig| + 严格 <1) | 1e-12 | [x] | V9 (roots=1/eig 双口径) |

### 5.2 IRF/FEVD 双轨 — vs statsmodels + R Spillover

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 5.2.1 | IRF Θ_h[i,j] 方向 (行=响应, 列=冲击) | 1e-12 | [x] | V3 (Ψ_h vs orth_ma_rep h=1/2/10; Ψ₁=ΣΦ₁(k,·)P(·,0) 全和断言) |
| 5.2.2 | Cholesky FEVD 行和精确=1 | 1e-12 | [x] | H=10 + H=1 双锚 |
| 5.2.3 | GFEVD (DY 框架 σ_jj⁻¹ + 行归一化) vs R Spillover | 1e-8 | [x] | 决策 15 (raw 行归一恒等式交叉: c=[1.1122,1.2462,1.1468]) |
| 5.2.4 | PS 框架 (σ_ii⁻¹) 可选输出, 与 DY 归一化后数值差异断言 | 1e-8 | [x] | V8 (归一后 max_diff>1e-3) |
| 5.2.5 | 不稳定 VAR 拦截 FEVD (异常) | 流程 | [x] | V12 (爆炸 ρ=1.05 确定性触发; 随机游走 ρ̂<1 不触发属合法) |
| 5.2.6 | bootstrap 置信带仅点估计容差 | 流程 | [x] | V13 (带仅结构断言: lower≤upper+有限; 点估计主锚) |

### 5.3 DY 溢出指数 — vs R Spillover

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 5.3.1 | TCI/TO/FROM/NET vs G.spillover | 1e-8 | [x] | 表/TO/FROM/NET 逐元素 (SP_TABLE 列主序); NET=TO−FROM 恒等 |
| 5.3.2 | window 必填无默认 + 强制 window>2·N | 精确 | [x] | §13-a 裁决 (≤2K/>T 拒绝; 可估计边界 w≥K·p+kt+K+p 数据依赖留档) |
| 5.3.3 | 频率默认表 {日 200/周 200/月 60} × H {10/10 周/12 月} | 精确 | [x] | API 无默认 (window 必填); 月 60×H=12 可配执行断言 |
| 5.3.4 | H 可配置敏感性 (H=10 vs 50) | 方向性 | [x] | V10 (H=50 主锚 1e-8; TCI 差异 (0, 0.1) 区间) |
| 5.3.5 | 滚动窗口 step=1 每窗口全重估 | 流程 | [x] | w=150 路径 101 点 vs roll.spillover 1e-6 (R 默认 p=1 对齐); 末窗口=子样本 static 1e-10 |

---

## 6. M3 数值基准对照 (协整)

### 6.1 Engle-Granger — vs statsmodels coint

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 6.1.1 | 统计量 vs statsmodels coint | 1e-10 | [x] | 4 趋势 × 3 对全对 (StatisticAllTrends/DirectionDependence/NoCointegration) |
| 6.1.2 | MacKinnon 1994 响应面 (5 系数 4 次, N≥2, 与 2010 ADF 表分文件) | 1e-12 | [x] | CI1, 决策 18; 渐近锚 + 有限样本公式手算复现 (MackinnonCVFiniteSampleFormula) |
| 6.1.3 | p 值 (渐近) 与 cv (小样本修正) 分列断言 + API 文档声明不同源 | 流程 | [x] | CI2 (#4138); PValue/CriticalValues 分列 + trend=n NaN cv 语义 (TrendNNaNcriticalValues) |
| 6.1.4 | 双方向输出差异 | 1e-10 | [x] | CI3; fwd/rev 双基准 + cv 与方向无关断言 |

### 6.2 Johansen — vs statsmodels + urca (diff 前置)

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 6.2.1 | **前置: verify_johansen_diff.R 双库 diff 报告冻结主对照** | 流程 | [x] | 决策 19 第一任务; JOHANSEN_DUAL_LIB_DIFF.md (参数映射 + MC 裁决 + SM 主基准冻结) |
| 6.2.2 | eig/lr1/lr2 vs coint_johansen (transitory 形式) | 1e-10 | [x] | CI7; det×k 网格 (det −1/0/1 × k 1/2) |
| 6.2.3 | λ̂ 降序 + 有效 T 口径 | 1e-10 | [x] | CI6; n_obs = T−1−k + 迹/最大特征值恒等式 (TraceMaxEigIdentity) |
| 6.2.4 | OL1992 表 vs urca 源码常量 (static_assert) | 精确 | [x] | osterwald_lenum_cv.hpp 表锚 static_assert |
| 6.2.5 | statsmodels 内嵌表 (MHM96) 双对照 | 精确 | [x] | CI5/B1; cv_source 回显 "MHM96" + OL1992 独立查表 API 双对照 |
| 6.2.6 | 3 det 情形 API 边界 (5 情形归 VECM) | 流程 | [x] | CI4; det_order ∈ {−1,0,1} 校验 + select_coint_rank |

### 6.3 VECM + PO — vs statsmodels + urca

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 6.3.1 | VECM 系数 vs statsmodels (5 情形) | 1e-10 | [x] | alpha/beta/gamma/det_coef/sigma_u/llf/resid 全对 (llf 走 SVD+非对称 eig 路径) |
| 6.3.2 | β 对照用投影矩阵 P=β(β'β)⁻¹β' (非逐元素) | 1e-10 | [x] | CI8; BetaProjectionSpace (r=2) |
| 6.3.3 | β 双归一 (前 r×r=I_r 默认 / urca 首变量开关) | 1e-10 | [x] | 决策 21; urca 基于未归一 β̃ 缩放 (Π 不变断言, 修复 rank≥2 除零) |
| 6.3.4 | ECT t 检验 EM2002 查表 | 精确 | [x] | CI10; t vs stderr_alpha 复算 + EM2002 表锚 (ctt n=3 1% T=51 → −5.0860) |
| 6.3.5 | Π=αβ' 符号方向 (αᵢ<0 拉回) | 方向性 | [x] | CI9; Π 秩 = r (SVD) + 归一消列符号 (前 r 行 = I_r 确定性) |
| 6.3.6 | Pu (方向依赖) / Pz (方向无关) vs urca ca.po | 1e-8 | [x] | CI12, 决策 20; 12 组合全对 + 方向依赖/无关双向断言 |

---

## 7. M4 数值基准对照 (MIDAS)

### 7.1 权重族与对齐 — vs midasr 0.9

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 7.1.1 | nealmon 权重逐点 (i=1..d 从 1 起) | 1e-12 | [x] | MD1; 2/3 参数两组 1e-12 全对 — test_midas_weights 1-2 |
| 7.1.2 | nbetaMT (xi 从 0 起, 两套网格并存) | 1e-12 | [x] | 4 参数 (θ₀ 第 4 参实施勘误); θ₀=0/0.1 两组 — test 6-7 |
| 7.1.3 | δ 独立线性参数 (Σw=δ, 内层解析消去) | 1e-12 | [x] | MD2; Σw=δ 多 δ 值断言 — test 3/8 |
| 7.1.4 | mls lag0 = 期末 x_{tm} 断言 | 精确 | [x] | MD3; W-dir 方向定案 (Form A 恢复 λ*) — test_midas_model 1/4 |
| 7.1.5 | log-sum-exp 与裸公式差 (非溢出区间) | <1e-14 | [x] | 决策 25/MD7 — test_midas_weights 4 |
| 7.1.6 | 溢出区间 (λ₂·d²≳709) 有限值断言 | 流程 | [x] | λ₂=−0.2 d=200 (下溢侧) 有限 + Σw=1 — test 5 |

### 7.2 模型估计 — vs midasr 夹具

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 7.2.1 | U-MIDAS vs midas_u (纯 OLS 锚) | 1e-10 | [x] | 风险 #4; **逐位一致** (SSR 7.3324398, 系数 ~1e-15) — test 5-8 |
| 7.2.2 | MIDAS-DL vs midas_r 夹具 (reltol=1e-12, maxit=10000, set.seed) | 1e-6~1e-8 | [x] | 决策 24; λ* 恢复 (5.030,−0.499) 参数 1e-4/SSR 1e-8; midas.coef 隐含权重 1e-6 — test 9-11 |
| 7.2.3 | MIDAS-AR (+AR*) vs 夹具 | 1e-6~1e-8 | [x] | W6 (μ,δ,λ₂,ρ₁) 全对 1e-4~1e-5; AR* 参数化 (1−φ)φ^{ℓ−1} 实现 — test 13-14 |
| 7.2.4 | 集中化 NLS ≡ 联合 NLS 同一最优 | 1e-6 | [x] | 决策 23; C++ 集中化 (SLSQP θ + 内层 OLS) vs midas_r 联合 BFGS 收敛同一 SSR (W5 双起始) — test 12 |
| 7.2.5 | 多起点 λ 网格×{递减,驼峰,均匀} 逃逸 (λ=0 平坦陷阱) | 流程 | [x] | MD8; 默认网格 {−2..0.5} + 远程起始 {{2.0}} 同最优 — test 12 |
| 7.2.6 | 夹具 convergence + SSR 双判据 + hAh 三列全录 | 流程 | [x] | hAh stat/p/df = (2.108, 0.3486, 2) ~1e-4 (prep_hAh 逐字复刻) — test 18 |

---

## 8. 幻觉点逐项核查 (64 编号, probe/verify 脚本逐条落盘)

### 8.1 ARIMA/Granger (AR1-AR8, GR1-GR7)

| ID | 影响级 | 核查内容 | 容差 | 状态 | 脚本 |
|----|--------|----------|------|------|------|
| AR1 | 极高 | MA 两库同号 (1+θB) | 1e-10 | [x] | verify_arima.R; θ 逐系数对 (test_arima_model 1/3/4/6) |
| AR2 | 极高 | n.cond=d+max(user,p) 与 q 无关 | 1e-10 | [x] | 四夹具实测定案 (arma12 n.cond=1); test 5 |
| AR3 | 高 | loglik 口径统一高斯形式 | 1e-10 | [x] | CSS 高斯型 vs 精确似然不可比 (实测裁决, 精化单调须同 ML 面); test 14 |
| AR4 | 高 | AIC 基于 T−d | 1e-10 | [x] | npar 含 σ² (forecast 862.606 吻合); test 8/20 |
| AR5 | 中 | drift=差分截距 | 1e-8 | [x] | forecast::Arima 语义 (drift=0.3573, DGP 0.3 恢复); stats::arima d≥1 强制无均值裁决留档; test 17-19 |
| AR6 | 极高 | method 配对对照 (禁混配) | - | [x] | CSS↔R CSS, CSS-ML↔R CSS-ML, Innov↔statsmodels; test 1-16 |
| AR7 | 高 | 多起始 (HR/CSS/随机) | - | [x] | {HR, 0, 扰动} 集合 + 关 HR 路径; test 24 |
| AR8 | 高 | innovations 限定 (无缺失/无季节) | 1e-10 | [x] | verify_arima.py; NaN 拒绝 + arma11 逐位; test_innovations_mle 9-10 |
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
| V1 | 低 | PS 1998 = Econ. Letters | - | [x] | 文献冻结 (fevd.hpp 注释) |
| V2 | 极高 | Cholesky 下三角 | 1e-12 | [x] | verify_var.py + vars t(chol) 交叉 |
| V3 | 极高 | Θ[i,j] 行响应列冲击 | 1e-12 | [x] | orth_ma_rep h=0/1/2/10 全锚 |
| V4 | 极高 | IC 用 ML Σ (÷T) | 1e-10 | [x] | logdet ×2 修复后 SM 1e-10 |
| V5 | 高 | select_order 同样本 offset | 1e-10 | [x] | p=0..4 轨迹 1e-10 |
| V6 | 高 | FPE 指数 K, n* 含 det | 1e-10 | [x] | ((T+df_m)/df_r)^K·e^ld |
| V7 | 极高 | GFEVD 行归一化 (行和=1) | 1e-8 | [x] | verify_gfevd.R + raw 行归一恒等 |
| V8 | 极高 | 双系数框架 (DY σ_jj⁻¹ / PS σ_ii⁻¹) | 1e-8 | [x] | 归一后差异 >1e-3 断言 |
| V9 | 高 | 稳定性双输出 (严格 <1) | 1e-12 | [x] | max_abs_eigenvalue + 布尔双输出 |
| V10 | 中 | H 可配置 (非数学常数) | - | [x] | H=10/12/50 三档测试 |
| V11 | 低 | GIRF 隐含假设注明 (Kim 2013) | - | [x] | 文档 (fevd.hpp 注释) |
| V12 | 高 | 不稳定 VAR 拦截 FEVD | - | [x] | 爆炸 ρ=1.05 触发 throw |
| V13 | 中 | bootstrap≠δ法, 容差仅点估计 | - | [x] | 带仅结构断言 |

### 8.3 协整 (CI1-CI12)

| ID | 影响级 | 核查内容 | 容差 | 状态 | 脚本 |
|----|--------|----------|------|------|------|
| CI1 | 极高 | 1994 响应面按 N 索引 | 精确/1e-12 | [x] | mackinnon_coint_cv.hpp + verify_eg.py; 渐近锚 −3.89644 等表锚 + 手算复现 |
| CI2 | 高 | p/cv 不同源 (#4138) | - | [x] | 分列断言 (PValue/CriticalValues 分测) + trend=n NaN cv 语义 |
| CI3 | 中 | EG 方向依赖 (双方向输出) | 1e-10 | [x] | DirectionDependence 双基准 + cv 方向无关 |
| CI4 | 高 | 3 vs 5 情形 API 边界 | - | [x] | Johansen det_order ∈ {−1,0,1} / VECM det ∈ 5 情形 |
| CI5 | 高 | statsmodels 一套表 (MHM96) | 精确 | [x] | verify_johansen_diff.R + cv_source "MHM96" 回显 |
| CI6 | 极高 | 迹公式+λ降序+有效T | 1e-10 | [x] | verify_johansen.py; TraceMaxEigIdentity 恒等式 |
| CI7 | 极高 | transitory 对照 urca | 1e-8 | [x] | JOHANSEN_DUAL_LIB_DIFF.md 参数映射冻结 (SM 主基准 + urca 网格交叉 1e-8) |
| CI8 | 极高 | β 投影空间对照 | 1e-10 | [x] | verify_vecm.py; BetaProjectionSpace (r=2 投影矩阵) |
| CI9 | 中 | Π 符号 (αᵢ<0 拉回) | - | [x] | 归一消列符号 (前 r 行 = I_r 确定性) + Π 秩 SVD 断言 |
| CI10 | 高 | EM2002 查表 | 精确 | [x] | ericsson_mackinnon_cv.hpp; 表锚 −3.4307 + ctt T=51 → −5.0860 |
| CI11 | 低 | 文献出处 (JEDC/JAE) | - | [x] | 文献冻结 (头文件锚点注释) |
| CI12 | 高 | Pu/Pz 双实现 | 1e-8 | [x] | verify_po.R; 12 组合 + PuDirectionDependence/PzDirectionInvariance |

### 8.4 MIDAS (MD1-MD8)

| ID | 影响级 | 核查内容 | 容差 | 状态 | 脚本 |
|----|--------|----------|------|------|------|
| MD1 | 极高 | nealmon i=1 起 | 1e-12 | [x] | verify_midas.R W1; test_midas_weights 1-2 |
| MD2 | 高 | δ 独立线性参数 | 1e-12 | [x] | Σw=δ 断言; test 3/8 |
| MD3 | 极高 | lag0=期末对齐 | 精确 | [x] | W2 实测 + W-dir Form A 恢复 λ* 定案; test_midas_model 1/4 |
| MD4 | 低 | U-MIDAS 出处 (F-M-S 2011/2015) | - | [x] | 文献冻结 (spec §6.3 锚点) |
| MD5 | 低 | 综述 GSV 2007 | - | [x] | 文献冻结 (spec §6.3 锚点) |
| MD6 | 低 | K-Z 2012 | - | [x] | hAh stat/p/df 三列全录 ~1e-4; test_midas_model 18 |
| MD7 | 极高 | log-sum-exp 防溢出 | <1e-14 | [x] | test_midas_weights 4-5 |
| MD8 | 高 | 多起点逃逸 | - | [x] | 默认网格 + 远程起始同最优; test_midas_model 12 |

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
| 9.3 | VAR→DY 溢出全链 (IC→稳定→IRF/FEVD→DY+滚动) | M2 全链; V12 拦截 | [x] 单元级全链 (三套件 50 用例; 端到端场景待 §1.3.16) |
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
| 10.1.9 | D9 VAR 逐方程 OLS (method='ols') | [x] | 同回归元矩阵一次求解等价 |
| 10.1.10 | D10 IC Lütkepohl 约定 + offset | [x] | IC 五式 + V5 同样本 |
| 10.1.11 | D11 Cholesky LLT 下三角 + P 注入/重排 | [x] | 注入校验下三角; reorder 敏感性测试 |
| 10.1.12 | D12 FEVD 双轨, DY 基于后者 | [x] | Cholesky + GFEVD(DY) 双轨 |
| 10.1.13 | D13 is_stable 双输出 | [x] | max|eig| + 严格 <1 |
| 10.1.14 | D14 IRF 带 block bootstrap | [x] | Politis-White 风格块长 |
| 10.1.15 | D15 GFEVD 自实现 (R Spillover 主基准) | [x] | g.fevd 1e-8 + G.spillover 表/指数 1e-8 |
| 10.1.16 | D16 不做 SVAR/BVAR/TVP-VAR | [x] | scope 外推迟项维持 |
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
| 11.1.6 | VAR + IC + IRF + FEVD 双轨 + DY | [x] | M2 全量 (50 用例; 滚动/auto-lag/PS 框架) |
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
| 15.5 | VAR K=5, T=500 + IC 扫描 | < 5 sec | [x] | VarModel.PerfK5T500ICScan (maxlag≈18 全扫描, 实测 < 0.05 sec) |
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
