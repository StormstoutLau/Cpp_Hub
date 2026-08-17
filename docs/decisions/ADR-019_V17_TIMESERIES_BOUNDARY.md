# ADR-019: v1.7 多变量时序与混频模块实施边界 (26+3 项)

**状态**: Accepted
**日期**: 2026-08-17
**版本归属**: v1.7 (Phase 7C)
**关联 Phase**: 7C
**决策者**: 架构组
**调研依据**: [PHASE7C_RESEARCH.md](../research/PHASE7C_RESEARCH.md) v1.1 (已经 3 agent 126 条全量审计)
**复核记录**: [ADR019_REVIEW_PILOT.md](../research/ADR019_REVIEW_PILOT.md) v1.1 — **新 R1-R4 调研门禁首轮 pilot, 全部通过后冻结**
**前置 ADR**: [ADR-016](ADR_INDEX.md#adr-016-金融时间序列实施边界-18-项) (金融时间序列实施边界) / [ADR-017](ADR_INDEX.md#adr-017-时序模块命名空间-cpphubv1timeseries) (命名空间) / [ADR-018](ADR_INDEX.md#adr-018-slsqp-优化器实现边界) (SLSQP, 约束 C5 复用)

---

## 背景

Phase 7B (v1.6) 交付 GARCH 族 + 单位根 + 方差比模块后, v1.7 扩展至多变量时序 (VAR/协整) 与混频回归 (MIDAS), 并回填三项遗留 (Ng-Perron / Zivot-Andrews / GARCH-M)。

调研阶段首次运行于 DEVELOPMENT_WORKFLOW.md v1.1 阶段 0 (断言分级证据框架): 调研报告经 3 agent 126 条全量审计 (修正 8 处实质错误), 决策集再经 ADR-019 复核 pilot (Phase 1 机械探针 + Phase 2 双盲重推导 5/5 TRUE + R2 引文抽检), R1-R4 门禁全部通过。

## 决策

### ARIMA/Granger (8 项)

| # | 决策 |
|---|------|
| 1 | ARIMA v1 只做 CSS + CSS-ML (CSS 起始→高斯似然 BFGS 精化); innovations 算法实现精确 MLE (限无缺失/无季节, 对标 `innovations_mle`); Kalman 留 v1.8 |
| 2 | 采用 (1+θB) 正号参数化, 与 R/statsmodels 直接对照 |
| 3 | d>0 默认无均值项, 提供显式 `drift` 选项 (forecast::Arima 语义) |
| 4 | Granger 主基准 statsmodels `grangercausalitytests` (4 统计量全对), 次基准 vars::causality (系统 F) |
| 5 | TY 只做增广 Wald, d_max 用户传入 (复用 ADF/PP/KPSS 预检验) |
| 6 | HAC-Wald 稳健版复用 Newey-West, 不新增估计器 |
| 7 | ARIMA 多起始点策略 (HR/CSS 预优化 + 随机重启), 否则 1e-10 对照假失败 |
| 8 | 不做: 季节 SARIMA / wild bootstrap p 值 (开放问题) |

### VAR/DY (8 项)

| # | 决策 |
|---|------|
| 9 | 估计只做逐方程 OLS (与系统 GLS 数值等价), 对齐 `method='ols'` [复核 B3: method 参数在 statsmodels 实现中未被使用, 文档与实现双重成立] |
| 10 | IC 用 Lütkepohl Ch.4/statsmodels 约定 (ML logdet + 系统参数量), select_order 强制同样本 (offset 机制) |
| 11 | Cholesky 统一 Eigen `LLT.matrixL()` 下三角; 暴露 P 注入接口 + 变量重排接口 (排序敏感性检验) |
| 12 | FEVD 双轨: Cholesky (行和=1) + generalized (行归一化); DY 指数仅基于后者, H 默认 10 可配置 |
| 13 | `is_stable` 返回 max\|eig\| + 严格 <1 双输出 |
| 14 | IRF 置信带复用 block bootstrap; δ 法渐近 SE 二期 |
| 15 | GFEVD 自实现 (statsmodels 无 [复核 B1/B2 双盲确认], 主基准 R Spillover 包 `g.fevd`/`G.spillover`) |
| 16 | 不做: SVAR/BVAR/TVP-VAR (v1.8+ DSGE 场景) |

### 协整 (5 项)

| # | 决策 |
|---|------|
| 17 | 先 EG (AEG) + coint_johansen 等价 API (3 det 情形); VECM 类再开 5 情形 |
| 18 | EG 临界值新增 MacKinnon 1994 协整响应面 (与 ADF 2010 表分文件, 注明 N=1 vs N≥2 适用域); p 值复刻 statsmodels 渐近版并文档声明与 cv 的差异 |
| 19 | Johansen 临界值主录 Osterwald-Lenum 1992 (urca 源码常量转录 constexpr + static_assert), statsmodels 内嵌表双对照 [复核 B1: statsmodels 仅一套 Johansen 表 (c_sjt/c_sja), 全库穷举确认]; **先跑双库数值 diff 决定主对照并冻结** |
| 20 | PO 检验纳入 (Pu/Pz 双实现, Pz 优先); ARDL/PSS 2001 出 scope |
| 21 | VECM β 默认 Phillips 归一 (前 r 行=I_r), urca 式首变量归一做开关; ECT t 检验仅给 Ericsson-MacKinnon 2002 查表 |

### MIDAS (5 项)

| # | 决策 |
|---|------|
| 22 | scope 内: MIDAS-DL (nealmon/nbeta/almonp/polystep/harstep 权重族, 两套网格并存) + U-MIDAS + MIDAS-AR (含 AR*); 逐字复刻 midasr 公式 |
| 23 | 估计: 集中化 NLS — 外层 SLSQP 只优化 λ, 内层 OLS(QR) 解析解 δ/截距/AR (Ghysels-Qian 2019); 联合 NLS 备选; 多起点 λ 网格×{递减/驼峰/均匀} |
| 24 | midasr 0.9 唯一主基准 [复核 B5 双盲: Python 生态四角度验证无维护良好 MIDAS 回归实现]: 夹具用收紧容差的 midasr 输出固化 CSV |
| 25 | exp 内部 log-sum-exp 防溢出 (验证与裸公式差 <1e-14) |
| 26 | scope 外 (v1.8+): midas_nlpr/sp/qr/imidas_r, amweights, gompertz/nakagami/lcauchy/genexp |

### 回填 (追加 3 项)

| # | 决策 |
|---|------|
| NP | 复用 ERS GLS 去势 (c̄=−7/−13.5) + 四统计量 + MAIC/MBIC/seq-t; 基准 = NP 原文公式钉死 + 文献 Table + EViews/Julia 输出对齐 (1e-6~1e-8); 临界值照 mackinnon_tables.inc 管线录 NP Table 1 |
| ZA | A/B/C 三模型, 固定 lag 主模式 + Baum 式预选 lag 对照模式; trim 默认 0.15 参数化; 临界值主用 ZA1992 论文表 (对齐 urca), MC 表可选 |
| GM | GARCH(1,1)-M + g∈{h,√h,log h} 三变体; 双锚 rugarch archpow=1/2 ↔ arch 8.0 `form='vol'/'var'/'log'` 精确对偶 (log 变体有数值对照不再限自洽测试); 系数容差 1e-4 + fix() 互验三步法; EGARCH-M/GJR-M 后续 |

## 依据与门禁记录

- 阻断性断言 → 双源状态映射: 见 [ADR019_REVIEW_PILOT.md](../research/ADR019_REVIEW_PILOT.md) §2
- Phase 1 机械探针: B1-B4 PROBE_SURVIVED (本地 statsmodels 0.14.6)
- Phase 2 双盲重推导: 5/5 TRUE (独立 auditor 一手取证)
- R1-R4 门禁: 全部通过 (2026-08-17), 含 1 处 R2 引文修正 (形态 II, 已同步登记块)
- 假设区 H1 (Julia NP MPT 系数疑偏) / H2 (PQ 2007 小样本微差) 不入正文, 以 [待定] 进 spec 开放问题节

## 后果

- **正向**: 决策集全部断言经双源验证, 可直接冻结进 PHASE7C_SPEC; GARCH-M log 变体获得数值对照路径; GFEVD/协整边界清晰
- **负向/成本**: MIDAS 无 Python 基准 → 夹具固化流程成本; Johansen 双库 diff 前置工作; GFEVD 自实现维护负担
- **开放问题**: SARIMA/wild bootstrap (决策 8), SVAR/BVAR/TVP-VAR (决策 16), ARDL/PSS (决策 20), v1.8+ MIDAS 扩展 (决策 26), H1/H2 假设区

## 关联

- **完整调研**: [PHASE7C_RESEARCH.md](../research/PHASE7C_RESEARCH.md) v1.1 (域调研 + 审计记录)
- **复核 pilot**: [ADR019_REVIEW_PILOT.md](../research/ADR019_REVIEW_PILOT.md) v1.1 (R1-R4 首轮验证)
- **门禁规范**: [ASSERTION_EVIDENCE_FRAMEWORK.md](../ASSERTION_EVIDENCE_FRAMEWORK.md) v1.1 / [DEVELOPMENT_WORKFLOW.md](../DEVELOPMENT_WORKFLOW.md) v1.1 §2.0/§4.2
- **兼容性约束**: C1-C7 见调研报告第四部分 (Eigen3 隔离 / 命名空间 / C ABI v1_7 前缀 / SLSQP 复用 / verify 脚本模式 / 容差分层)
- **后续**: PHASE7C_SPEC 编写 (R4 清零 → spec 冻结 → 实施 → G 系列门禁)
