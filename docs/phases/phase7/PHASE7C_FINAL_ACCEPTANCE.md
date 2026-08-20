# Phase 7C v1.7 最终验收报告 — 多变量时序与混频模块

> **文档类型**: 最终验收报告 (Final Acceptance Report)
> **审计对象**: Phase 7C v1.7 完整 scope (M0 回填 + M1 ARIMA/Granger + M2 VAR/DY + M3 协整 + M4 MIDAS + §7 端到端集成)
> **关联文档**:
>   - [PHASE7C_SPEC.md](./PHASE7C_SPEC.md) v1.2 (执行规格)
>   - [PHASE7C_ACCEPTANCE_CHECKLIST.md](./PHASE7C_ACCEPTANCE_CHECKLIST.md) (逐项审计清单, 331 项, 64 幻觉点逐点状态)
>   - [ADR-019 多变量时序实施边界](../../decisions/ADR-019_MULTIVARIATE_TS_BOUNDARY.md) / ADR-016 / ADR-017 / ADR-018
>   - [DEVELOPMENT_LOG.md](../../DEVELOPMENT_LOG.md) Phase 7C 五轮实施记录
> **验收日期**: 2026-08-19 (初验) / 2026-08-20 (C-1/C-2 闭环收口)
> **验收员**: Scott (self-review, 符合项目 solo developer 工作流)
> **验收标准**: 全量回归三平台通过 + 多基准库数值对照 (statsmodels/urca/midasr/vars/Spillover 分层容差) + 64 幻觉点全覆盖 + Eigen3 隔离 (单变量模块零 Eigen) + ADR-019 26+3 决策全对齐

---

## 1. 验收基本信息

| 项 | 值 |
|----|-----|
| 项目名称 | Cpp_Hub |
| 审计阶段 | Phase 7C v1.7 (多变量时序与混频 — M0 回填 NP/ZA/GM + M1 ARIMA/Granger + M2 VAR/IRF/FEVD/DY + M3 协整三件套 + M4 MIDAS + 端到端集成) |
| 前置版本 | v1.6 Phase 7B 金融时间序列 (2207/2207 主控 / 2189×2 GCC) |
| 后置版本 | v1.7 多变量时序与混频 (**2458/2458** 主控站 / **2440/2440** ×2 A/B 站, 三平台全绿) |
| Phase 7C 新增测试 | **251 个** (M0 回填 49 + M1∥M4 71 + M2 51 + M3 58 + 收尾轮 Granger 16 + 集成 6) |
| 累计测试 | 2458 (主控站 MSVC) / 2440 ×2 (A/B 站 GCC, 实测全绿) |
| 验收方法 | TDD 实现 + 多基准库分层对照 (1e-10/1e-8/1e-6 分层) + 三平台 bundle 中继验证 + ADR 逐项对齐审计 |
| 里程碑轮次 | M0 (069264c) → M1∥M4 (67b5450) → M2 (7d64939) → M3 (c008f46) → 收尾轮 (Granger+集成) |

> 主控站与 GCC 18 个测试差额为平台专属测试 (Windows/MSVC-only 用例), 7B 期已知口径, 非功能差异。

---

## 2. 全量回归测试结果

### 2.1 三平台矩阵 (2026-08-20 全轮闭环实测)

| 平台 | 编译器 | 测试总数 | 通过数 | 失败数 | 测试耗时 | 验证方式 | 状态 |
|------|--------|---------|--------|--------|---------|---------|------|
| 主控站 (Win10) | MSVC 19.50 (VS2026) | 2458 | 2458 | 0 | 141.63 sec (-j8) | 增量构建 + 全量 ctest | ✅ |
| A 工作站 (`scott-lau-NEX.local`) | GCC 13.3 | 2440 | 2440 | 0 | 428.73 sec | bundle 增量 ff c008f46→006cec6 + rebuild + ctest | ✅ |
| B 工作站 (`scott-lau-GTR-Pro.local`) | GCC 13.3 | 2440 | 2440 | 0 | 445.17 sec | 同法 (mDNS 通道) | ✅ |

> 五轮全部三平台全绿闭环 (M0-M3 四轮 229 用例 + 收尾轮 22 用例 Granger 16 + 集成 6); A/B 站自有代码 0 警告 (10 条 = autodiff 第三方 9 + GCC13 系统头误报 1, 口径与前四轮一致)。

### 2.2 GitHub Actions CI (run #60/#61/#62, 67b5450/d581f30/006cec6)

| Job | 结果 |
|-----|------|
| Build & Test (windows-2022, MSVC) | ✅ success |
| Build & Test (ubuntu-latest, GCC) | ✅ success |
| C ABI (windows-2022) | ✅ success |
| C ABI (ubuntu-latest) | ✅ success |

> run #60 首次承载 M1∥M4 全量 2327 用例双平台绿; **run #62 (收尾轮 006cec6) completed success — 首次承载全量 2458 用例双平台绿** (2026-08-20 确认, 条件 C-2 闭环)。

### 2.3 测试数演进 (五轮)

| 轮次 | 内容 | 新增测试 | 主控站累计 | A/B 站累计 | Commit |
|------|------|---------|-----------|-----------|--------|
| M0 | NP (18) + ZA (15) + GM (16) 回填 | 49 | 2256 | 2238 | 069264c |
| M1∥M4 | ARIMA (24+12) + MIDAS (16+19) | 71 | 2327 | 2309 | 67b5450 |
| M2 | VAR (21) + IRF/FEVD (18) + DY (12) | 51 | 2378 | 2360 | 7d64939 |
| M3 | EG (14) + Johansen (18) + VECM (16) + PO (10) | 58 | 2436 | 2418 | c008f46 |
| 收尾 | Granger (16) + 端到端集成 (6) | 22 | **2458** | **2440** (实测双绿) | fc35652+006cec6 |

### 2.4 跨平台验证执行细节

- A/B 站 github.com 直连受阻, 全程沿用 7B 期 **bundle 中继工作流**: 主控站 `git bundle create` + scp 至站内, `git clone/fetch <bundle>` + cmake `-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=<本地源>` 零外网依赖验证; 五轮全部增量 ff (无 fresh clone 重复成本)。
- B 站 IPv4 动态地址 ping 不通但 mDNS 主机名 (IPv6) 通道正常 — "B 站离线"结论已作废, 以主机名复核为准。
- 差额 18 (主控 2436 − A/B 2418, M3 轮口径) = 平台专属用例, 与 7B 期一致; **M0-M3 新增 229 用例 GCC 双站全数运行, 三平台零数值偏差**。

---

## 3. 实施范围汇总

### 3.1 新增头文件 (25 个, 按 ADR-017 命名空间)

| 里程碑 | 文件 | 用途 |
|--------|------|------|
| M0 回填 (4) | `unit_root/ng_perron_test.hpp` / `unit_root/zivot_andrews_test.hpp` / `unit_root/np_tables.hpp` / `garch/garch_m_model.hpp` | MZα/MZt/MSB/MPT + MAIC; ZA 三模型 + trim; NP Table 1 constexpr; GARCH-M 三变体 |
| M1 (4) | `arima/arima_model.hpp` / `arima/innovations_mle.hpp` / `arima/hannan_rissanen.hpp` / `arima/granger_test.hpp` | CSS + CSS-ML; B&D 精确 MLE; HR 起始值; F/χ²/LR + TY + HAC-Wald |
| M2 (6) | `var/multivariate_data.hpp` / `var_model.hpp` / `var_select.hpp` / `irf.hpp` / `fevd.hpp` / `dy_spillover.hpp` | 数据载体; 逐方程 OLS + IC 五式; auto-lag; Φ 递推 + 块 bootstrap; Cholesky + GFEVD 双轨; TCI/TO/FROM/NET + 滚动 |
| M3 (7) | `cointegration/engle_granger.hpp` / `johansen_test.hpp` / `phillips_ouliaris.hpp` / `vecm_model.hpp` / `osterwald_lenum_cv.hpp` / `mackinnon_coint_cv.hpp` / `ericsson_mackinnon_cv.hpp` | EG 两步法; 迹/最大特征值; Pu/Pz; VECM 5 情形 + β 双归一; OL1992/MHM96; MacKinnon 1994; EM2002 ECT t |
| M4 (4) | `midas/mixed_freq_data.hpp` / `midas_weights.hpp` / `midas_model.hpp` / `midas_diagnostics.hpp` | mls 期末对齐; 5 权重族 + log-sum-exp; DL/AR/AR*/U-MIDAS 集中化 NLS; 残差诊断 + hAh |

**Eigen3 隔离 (ADR-013/017)**: 单变量模块 (unit_root/garch/arima/midas = 12 头文件) **零 Eigen include** (grep 断言); Eigen include 恰好 9 处全落 var (5) + cointegration (4); 新 CMake target `cpphub_timeseries_mat` 链接 Eigen3 INTERFACE — 依赖单向: 单变量 ← 多变量, 无反向。

### 3.2 新增测试套件 (16 套, 251 用例)

见 §2.3 演进表。spec 名义 263 用例, 实施定案 251 (集成场景 8→6: 每场景一条完整链路 TEST, 场景内多断言; 总数核算 2207+251=2458 主控站吻合)。基准锚定沿用硬编码 `constexpr` baseline 策略 (`coint_baseline.inc` 398 数组 / `var_baseline.inc` 等), 测试零运行时外部依赖。

### 3.3 基准验证脚本 (19 个) + 临界值表 (6 件)

| 类别 | 数量 | 清单 (tests/fixtures/timeseries/) |
|------|------|------|
| verify Python (statsmodels/arch) | 8 | np_stata (占位管线) / za / gm / arima / granger / var / eg / johansen / vecm |
| verify R (urca/vars/Spillover/midasr/rugarch) | 8 | za / gm / arima / var / gfevd / johansen_diff / po / midas |
| probe/gen 辅助 | 3+ | probe_midas_form / dump_var_r_values / gen_phase7c_{m2,m3}_baseline |
| 临界值表 (.inc → 4 件以 include 树头文件交付) | 6 | np_table1 / za1992_cv / za_mc_cv / ol1992_cv / mackinnon1994_coint / em2002_ect_cv — 后 4 件 constexpr + static_assert 形态 (np_tables/osterwald_lenum_cv/mackinnon_coint_cv/ericsson_mackinnon_cv) |

> 基准库版本冻结: statsmodels 0.14.x (M1/M2/M3 主锚) / R urca (PO 主 + Johansen/ZA 交叉) / R midasr 0.9 (MIDAS 唯一主基准) / R vars (VAR 交叉) / R Spillover (GFEVD/DY 主基准) / R rugarch (GM archpow) / Stata dfgls (NP MAIC, 占位) / arch 8.0 (7B 主锚沿用)。

---

## 4. 数值基准与幻觉点核查

### 4.1 基准策略 (分层容差)

- **1e-10 逐位层**: 同库同路径对照 — statsmodels (EG/Johansen/VECM/VAR/Granger 四统计量/innovations MLE/U-MIDAS)、midasr midas_u (实测逐位)
- **1e-8 交叉层**: R 库交叉 — urca (PO 12 组合/Johansen 网格/ZA 固定 lag)、vars (IC 逐位)、Spillover (GFEVD 表/DY 指数/滚动)
- **1e-4~1e-6 落点层**: 不同优化器落点差 — GM vs arch (数值梯度 vs 解析)、ARIMA CSS vs R optim、midasr NLS (reltol=1e-12 + seed 夹具收紧)
- **精确层**: 临界值表 (OL1992/ZA1992/MC 表/NP Table 1/MacKinnon 1994/EM2002) constexpr + static_assert 锚

### 4.2 幻觉点覆盖 (64 编号全核查, checklist §8 逐点状态)

| 族 | 编号 | 状态 |
|----|------|------|
| ARIMA (AR1-AR8) | 8 | ✅ 全落码 (n.cond 定案/drift 语义/θ 谱等价类/多起始) |
| Granger (GR1-GR7) | 7 | ✅ 全落码 (df 公式/TY df=k/HAC 三态/列序陷阱/方向复现/I(1) 失效) |
| VAR/DY (V1-V13) | 13 | ✅ 全落码 (logdet 因子 2/Σ_mle ÷T/布局三态/LLT 陷阱/V12 拦截) |
| 协整 (CI1-CI12) | 12 | ✅ 全落码 (p 与 cv 分列/llf 数值路径/β 双归一/EM2002 三重验证/Pu-Pz 方向) |
| MIDAS (MD1-MD8) | 8 | ✅ 全落码 (mls 期末/nbetaMT θ₀/log-sum-exp/期初起窗) |
| 回填 NP/ZA/GM (NP1-6, ZA1-5, GM1-5) | 16 | ✅ 全落码 (MAIC 口径/MC 表 ≠ 论文表/archpow 映射/fix() bug) |

### 4.3 实施期关键幻觉点实测记录 (spec/文档 vs 源码差异)

1. **V6 logdet 漏因子 2** (核心 bug): LLT `Σ log(L_ii)` 少乘 2 (det=ΠL_ii²) ⇒ IC 全家偏差 — 手算 det 反推定位, 一行修复解 12 个失败
2. **CI6 llf 数值路径**: eigh 对称路径与 statsmodels 差 ~7e-11, 经 λ→log(1−λ₀)→T/(2(1−λ₀)) 放大 **202×** 破坏 llf 1e-10 — SVD + 非对称 eig + LU det 逐字对齐后 diff=1.5e-11
3. **AR5 d≥1 均值语义**: stats::arima d≥1 强制无均值 (伪根 0.9985 退化), 正解 = forecast::Arima drift (差分截距) — 双路径对照留档
4. **GR6 statsmodels 列序**: 二维输入第二列 cause 第一列 — 方向陷阱注释落码 + 复现断言
5. **NP 固定样本口径坐标错位**: dyd[T−1] 堆越界读 (069264c 勘误, 60 次循环 0 失败根除)
6. **MD nbetaMT 4 参数**: midasr 实施与文献 θ₀ 第 4 参勘误
7. **GM arch fix() bug**: arch 8.0 ARCHInMean.fix() llf=−3002.8 自身缺陷 — C++ 充当独立评估器完成隔离验证 (可提上游 issue)
8. **Johansen 双库 diff**: SM 与 urca 辅助回归 level 项索引不同 (y_{t−k} vs y_{t−1}) — JOHANSEN_DUAL_LIB_DIFF.md 冻结 + MC 裁决渐近等价小样本不同 → SM 主基准

---

## 5. 端到端集成测试 (spec §7, 6/6 通过)

| # | 场景 | 验证链路 | 关键断言 | 状态 |
|---|------|---------|---------|------|
| 1 | 单位根诊断全链 | RW 水平 (ADF/DF-GLS/NP/ZA) → 差分 → ARIMA → LB | 不拒绝→拒绝→残差白噪声; 断点→分段平稳 | ✅ |
| 2 | Granger 因果链 | I(1) 双序列三路对照 | 水平 F 伪显著 / 差分 F 真不显著 / 水平 TY 真显著 (GR7) | ✅ |
| 3 | VAR→DY 溢出全链 | IC→稳定→IRF/FEVD 双轨→DY 静态+滚动 | ΣTO=ΣFROM=TCI 恒等式 + V12 爆炸 VAR (ρ=1.05) 拦截 | ✅ |
| 4 | 协整→VECM 全链 | EG 筛选→Johansen rank→VECM→ECT→β 空间 | 夹具链 + 真实拉回 DGP: ECT t<cv 双方程 + α 双负 | ✅ |
| 5 | MIDAS 混频预测 | 月-日对齐 (mls 期末) + DL/U-MIDAS + 预测 | MD3 期初起窗 + 收敛性 | ✅ |
| 6 | GARCH-M 风险溢价 | 三变体→λ 显著性→vs 无 M | 强信号模拟 (T=3000, λ=1.2): sandwich t>1.96 + λ 恢复 ±0.4 | ✅ |

**实施期修正原则**: 夹具 DGP 检验力不足时**新增强信号对照** (场景 4 拉回 DGP / 场景 6 λ=1.2) 而非放宽断言; 恒等式/方向断言错误**修正断言**而非改实现 (场景 2 TY 方向 / 场景 3 DY 恒等式); 数据生成缺陷**修数据** (场景 3 独立双噪声源防共线)。

---

## 6. 验收标准核查 (spec §8 对照)

| 标准 | 要求 | 实际 | 状态 |
|------|------|------|------|
| 测试通过率 | 全量通过 | 主控 2458/2458; A/B 站 M0-M3 四轮 2418/2418×2, 收尾轮待执行 | ✅/⏳ |
| 数值基准 | 分层容差 (§4.1) | 1e-10 主锚 + 1e-8 交叉 + 精确表锚全覆盖 | ✅ |
| 幻觉点覆盖 | 64 编号 | 64/64 逐点核查落码 (+10 scope 外不计) | ✅ |
| Eigen3 隔离 | 单变量零引入 | 12 头文件 0 Eigen; 9 include 全落 var/cointegration | ✅ |
| 跨平台一致性 | 三平台 ctest | M0-M3 229 用例三平台全绿; 收尾轮 22 用例 C-1 | ✅/⏳ |
| 命名空间 | ADR-017 | `cpphub::v1::timeseries::{arima,var,cointegration,midas}` | ✅ |
| ADR-019 | 26+3 决策 | 29/29 逐项对齐, 无静默偏离 | ✅ |
| Scope 边界 | §11 不越界 | SARIMA/wild bootstrap/SVAR/BVAR/TVP-VAR/ARDL/MIDAS 扩展/DCC/Kalman/长记忆族 grep 零匹配 | ✅ |

---

## 7. 复用验证 (不重复实现, checklist §13 十项全勾)

| 复用项 | 来源 |
|--------|------|
| NP ← DF-GLS 去势 | 7B df_gls_test (ERS 变换复用) |
| ZA ← ADF 引擎滞后框架 | 7B adf_test |
| GM ← GARCH backcast/似然/SLSQP/sandwich | 7B garch_model (λ=0 退化 ≡ filter_garch11 逐位证复用) |
| ARIMA/MIDAS 外层 ← SLSQP | ADR-018 optimizer (12/12 无退化 ×5 轮) |
| Granger HAC ← NW vcov | v1.5 hac (Schwert 规则同源) |
| IRF 置信带 ← block bootstrap | v1.5 |
| PO 长期方差 ← long_run_variance | 7B unit_root_common |
| MIDAS 内层 ← OLS/QR | v1.5 |
| 诊断 ← JB/LB/multiple_test_correction | Phase 7A |
| EG 第二步 ← adf_test (nc + Schwert + AIC) | 7B |

---

## 8. 实施期问题与修复记录

### 8.1 测试设计与断言问题 (收尾轮 6 处 + 前期各轮留档)

| # | 问题 | 修复 |
|---|------|------|
| 1 | 场景 2 TY 断言方向与 DGP 预期相反 | 修正为 p<0.05, 三路对照叙事理顺 |
| 2 | 场景 3 DY 恒等式混淆 (TO+FROM vs ΣTO=ΣFROM=TCI) | 修正恒等式 + NET 逐项 1e-12 |
| 3 | 场景 3 爆炸 VAR 双列同噪声共线 (sigma not PD) | 独立双噪声源 (seed 77/78) |
| 4 | 场景 4 夹具 DGP 无拉回 → ECT 断言必败 | 夹具断有限性 + 新增真实拉回 DGP 对照 |
| 5 | 场景 6 λ 检验力不足 (t=0.59) | 强信号模拟 T=3000 λ=1.2 |
| 6 | gtest_discover_tests 缓存过期 (6 场景仅注册 1) | 删 exe 强制重链接 + 发现 |

### 8.2 基准事实修正 (spec 幻觉)

见 §4.3 (8 项)。全部通过 probe 脚本 (源码打印 + 小样本手算 + 双库 diff) 验证后修正; R 门禁回溯条款执行: 无静默改公式, 4 处实施勘误均属对照边界或勘误类 (§10.2.8)。

### 8.3 环境/基础设施问题

| # | 问题 | 处置 |
|---|------|------|
| 1 | 收尾期命令执行通道 (RunCommand/MCP) 间歇性失效 | A/B 站收尾轮顺延 (C-1), 主控站验证先行闭环 |
| 2 | B 站 IPv4 ping 不通误判离线 | mDNS 主机名复核通道正常, 结论作废 |
| 3 | 同文件并行 Edit 竞态 (M2 轮再犯 v1.2 教训) | 串行 Edit 铁律重申 |

---

## 9. 遗留事项与后续行动

### 9.1 v1.7 发布前条件 (checklist §17.2; C-1/C-2 已闭环)

| # | 事项 | 优先级 | 状态 |
|---|------|--------|------|
| C-1 | A/B 站收尾轮验证 | 高 | ✅ **闭环 (2026-08-20)**: A 站 2440/2440 (428.73s) / B 站 2440/2440 (445.17s) |
| C-2 | 收尾轮 commit 推送 + CI run 确认 | 高 | ✅ **闭环 (2026-08-20)**: fc35652+006cec6 推送, CI run #62 success (全量 2458 双平台) |
| C-3 | Stata 装机后补 NP MAIC 1e-10 硬断言 (管线就绪) | 中 (v1.8 前) | 待办 |
| C-4 | 性能独立用例补设 (ZA/GM/ARIMA 规模用例) 或 N/A-级放行留档 | 低 (v1.8 评估) | 待办 |

> **v1.7.0 发布条件全部满足** — 可执行 tag 发布。

### 9.2 后续版本 scope (v1.8+, spec §9.2 推迟项)

- SARIMA / Granger wild bootstrap / Kalman statespace / ARDL-PSS
- SVAR / BVAR / TVP-VAR / DCC / CCC 多元波动率
- MIDAS 扩展族 (nlpr/sp/qr/imidas_r/amweights/带 p 后缀权重)
- 长记忆族 (APARCH/FIGARCH/IGARCH/ARFIMA/HYGARCH)

---

## 10. 验收结论

| 维度 | 权重 | 状态 | 备注 |
|------|------|------|------|
| 测试通过率 | 25% | ✅ | 三平台全绿: 主控 2458/2458 + A/B 站 2440/2440 ×2 (五轮零退化) |
| 数值基准 | 25% | ✅ | 五库分层对照 (1e-10 主锚/1e-8 交叉/精确表锚) |
| 幻觉点覆盖 | 15% | ✅ | 64/64 逐点核查落码 |
| ADR 对齐 | 10% | ✅ | ADR-019 26+3 全对齐 + Eigen 隔离 + scope 边界 |
| 集成测试 | 10% | ✅ | 6 场景端到端全通过 (三平台一致) |
| 复用验证 | 10% | ✅ | 十项复用零重复实现 (退化锚逐位证明) |
| CI/跨平台 | 5% | ✅ | run #60/#61/#62 三连绿; #62 首载全量 2458 双平台 |
| **总计** | **100%** | **✅ 有条件通过 → 发布条件已满足** | 331 项 324 过 (97.9%), 未过 7 项全为预批占位/N-A 级放行; C-1/C-2 已闭环 (2026-08-20), C-3/C-4 为 v1.8 前非阻塞项 |

**验收结论**: ✅ **有条件通过 → v1.7.0 发布条件全部满足** — Phase 7C v1.7 在功能、数值、幻觉点、ADR、scope、规范、复用、风险、CI/跨平台九个维度全部达标; 三平台全量回归零退化 (主控 2458/2458 + A/B 站 2440/2440 ×2, 2026-08-20 实测); CI run #62 首次承载全量 2458 用例双平台绿。仅余 C-3 (Stata NP 硬断言) / C-4 (性能独立用例) 两项 v1.8 前非阻塞项。

**关键成果**:

1. **251 个新测试** (16 套件), 累计 2458 (主控), M0-M3 229 用例三平台零数值偏差
2. **25 头文件全 header-only**: 单变量 12 个零 Eigen + 多变量 13 个 Eigen 隔离 (新 target cpphub_timeseries_mat)
3. **64 幻觉点逐点核查**: 含 logdet 因子 2 / llf 数值路径 202× 放大 / arch fix() bug / Johansen 双库 diff 四项 spec 级发现
4. **6 库基准矩阵**: statsmodels/urca/midasr/vars/Spillover/rugarch 分层容差全对齐 + 6 临界值表 constexpr 锚
5. **6 场景端到端集成**: 单位根→ARIMA / Granger 三路 / VAR→DY / 协整→VECM / MIDAS 混频 / GARCH-M 风险溢价全链路打通
6. **v1.7 时序版图闭环**: 单变量 (7B GARCH/单位根/方差比 + 7C NP/ZA/GM/ARIMA/Granger/MIDAS) + 多变量 (VAR/IRF/FEVD/DY/协整/VECM) — C++ 侧与 statsmodels/urca 生态互补的完整金融时间序列基础设施

**Reviewer**: Scott (self-review, 2026-08-19)
