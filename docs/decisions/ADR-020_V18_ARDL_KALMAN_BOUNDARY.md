# ADR-020: v1.8 ARDL/PSS 边限检验与 Kalman 状态空间实施边界 (12 项)

**状态**: Accepted
**日期**: 2026-08-21
**版本归属**: v1.8 (Phase 8)
**关联 Phase**: 8 (M1/M2)
**决策者**: 架构组 (Scott 2026-08-21 裁决: 12 项全部按调研建议直接定稿)
**调研依据**: [PHASE8_RESEARCH.md](../research/PHASE8_RESEARCH.md) v1.1 (双 agent 并行域调研 + 审计二轮通过, 3 处修正)
**Scope 依据**: [PHASE8_SCOPE_RESEARCH.md](../research/PHASE8_SCOPE_RESEARCH.md) v1.3 (D-1 已裁决: 双免费午餐)
**前置 ADR**: [ADR-019](ADR_INDEX.md#adr-019-v17-多变量时序与混频实施边界-263-项) 决策 1 ("Kalman 留 v1.8") / 决策 20 (ARDL/PSS 出 v1.7 scope) / [ADR-017](ADR_INDEX.md#adr-017-时序模块命名空间-cpphubv1timeseries) (命名空间) / [ADR-018](ADR_INDEX.md#adr-018-slsqp-优化器实现边界) (SLSQP, K3 复用)

---

## 背景

Phase 7C (v1.7) 闭环后 (v1.7.0/v1.7.1 双 tag 发布), Phase 8 按 D-1 裁决执行"双免费午餐"结构: M1 ARDL/PSS 边限检验 + M2 Kalman/statespace (兑现 ADR-019 决策 1 遗留的 "Kalman 留 v1.8") + M3 战略选项 (D-2 另行裁决, 不阻塞)。

"免费午餐"依据 (scope 调研预判 + 域调研确认): 双端基准成熟 (statsmodels ≥0.13 ARDL 全家 + statespace 完整开源 / R ARDL 0.2.5 经 JAE 论文验证) / 复用率高 (UECM = OLS + Wald; Kalman = Eigen 载体 + innovations 交叉验证) / 成本周级。

调研报告 v1.0 经审计二轮复验 (v1.1): 修正 3 处 (Chandrasekhar 版本归属 v0.12.0 三方一手实证 / statsmodels nsim 文档 bug 实锤 / Narayan k≤7 边界降级 B-), Statalist Kripfganz 帖子独立印证 PSS Case III k=1 锚点。12 项决策建议中 10 项依据一手来源 (PSS 2001 原文逐表提取 / statsmodels 0.14.6 源码 / R ARDL 0.2.5 源码), 无阻断性问题。

## 决策

### M1 ARDL/PSS 边限检验 (7 项)

| # | 决策点 | 决策 |
|---|--------|------|
| A1 | CV 实现路线 | **响应面多项式** (statsmodels `pss_critical_values.py` 路线: ~500 行系数 + 3 阶多项式求值, 键 (k, case, I1)), PSS 2001 原表值作 static_assert 锚。拒绝硬编码 PSS 表 (pssbounds 模式, 无法给 p 值) 与 MC 现场模拟 (引入随机性) |
| A2 | 模型范围 | ARDL(p,q) 单方程, k 个 force 变量 (每变量独立 q_j); **不做** panel/quantile/Fourier/NARDL (ardlverse 系全部排除, CV 错误实锤) |
| A3 | 阶数选择 | 网格搜索 argmin IC (aic/bic/hqic), 对齐 `ardl_select_order` 语义 (同一 hold_back 样本); **不做** glob=True 指数搜索 (v1.9 评估) |
| A4 | t 检验覆盖 | **仅 Case I/III/V** (PSS Table CII 边界, 原文仅此三 case 有 t 界); Case II/IV 请求 t → 报错并说明 (镜像 pssbounds 行为) |
| A5 | 小样本 CV | 响应面天然按 T 插值 (KS 2020 路线); **不做** Narayan 分档表 (k≤7 边界两轮未获一手来源, 降级 B-) + **不做** MC 模拟 (asymptotic=False 语义); 小样本警告由文档提示 |
| A6 | 同期项 | UECM 含同期 Δx_j (causal=False 语义); causal 模式 v1.9 |
| A7 | 稳健推断 | v1.8 仅 cov_type='nonrobust' (表分布仅同方差成立); HC/HAC Wald v1.9 评估 |

### M2 Kalman/statespace (5 项)

| # | 决策点 | 决策 |
|---|--------|------|
| K1 | scope | SSModel 基类 (DK 记号) + 常规 Kalman (DK Ch4) + 平稳初始化 (Lyapunov + 对称/正定守卫) + approximate diffuse (κ=1e6) + ARMA(p,q) Harvey 形 + 嵌入差分 ARIMA 精确似然 MLE + innovations 互检; 对齐 statsmodels SARIMAX 默认口径四件套 (simple_differencing=False / κ=1e6 / 参数变换 / MA 正号) |
| K2 | 推迟 v1.9+ | exact diffuse (de Jong) / 平滑器 RTS / 缺失数据 / 多元 (VARMAX/DynamicFactor/UCM) / Chandrasekhar / sqrt filter / univariate filter (多元才需要); 每项有明确代价证据, 平滑非似然必需 |
| K3 | MLE 优化器 | 复用 ADR-018 SLSQP/BFGS (ARIMA CSS-ML 同款) + transform_params 平稳化变换 (镜像 statsmodels enforce_stationarity); 零新增基础设施 |
| K4 | σ² 处理 | **非集中化** (镜像 statsmodels SARIMAX 全参数优化), 与 innovations_mle 的集中化形成口径差异声明; 交叉验证以 loglik 为主锚 (7C 先例: 落点差在优化器容差层) |
| K5 | Eigen 载体 | 新增 `include/cpphub/timeseries/statespace/` 挂 `cpphub_timeseries_mat` (INTERFACE) — 镜像 var/cointegration 隔离模式 (CMakeLists L52-58 现成模式) |

## 依据与审计记录

- 调研: PHASE8_RESEARCH.md v1.0 双 agent 并行域调研 (WebSearch/WebFetch + PSS 2001 原文 PDF 逐表提取 + statsmodels/R/Stata 源码逐行比对) + 本地代码库复用映射自查
- 审计二轮 (v1.1, 2026-08-21): 全部断言独立双盲复验, 3 处修正 —
  1. Chandrasekhar 版本归属 v0.12.0 (原误记 v0.14; 三方一手实证: 官方 release notes + GitHub tag v0.12.0 2020-08-27 + Chad Fulton 作者博客/PR #6411)
  2. statsmodels 稳定版文档 nsim 说明 bug 实锤 (与 asymptotic/seed 语义自相矛盾, 源码逻辑为准)
  3. Narayan 2005 k≤7 边界降级 B- (Case II-V 覆盖与 T=30-80 分档核实, k 上限未获一手来源) → 直接支撑 A5 不采 Narayan 分档表
- 独立印证: Statalist Kripfganz 帖子印证 PSS Case III k=1 5% 锚点 (F 4.94/5.73, t −2.86/−3.22), 与 statsmodels `crit_vals` 数值吻合
- 幻觉点核查清单 19 条 (域 A 9 条 + 域 B 10 条) 已就绪, 供 PHASE8_SPEC §幻觉点核查直接引用

## 后果

- **正向**: M1 (OLS+Wald 纯复用, 风险低) 先行热身, M2 (Eigen statespace 载体) 兑现 ADR-019 遗留; innovations 仲裁协议 (T=300, 2.6e-8/1.2e-6) 从一次性升级为常驻互检; 两模块完全正交可并行 (7C M1∥M4 模式)
- **负向/成本**: statsmodels 响应面系数转录 (~500 数字) 需逐键 static_assert 锚 + verify 脚本全量 diff; Case II/IV 受限确定性项的 UECM 约束构造需 R/statsmodels 双端数值对照; Lyapunov 负方差静默失效是最高风险项 (守卫强制)
- **开放问题** (spec 阶段裁决): bounds_test 返回结构是否含 p 值 (建议含, statsmodels 口径) / ARDL trend 枚举语义细节 / loglikelihood_burn 默认值
- **不阻塞**: D-2 (M3 战略选项) 排期最晚, 可与 M1/M2 开发并行决策

## 关联

- **完整调研**: [PHASE8_RESEARCH.md](../research/PHASE8_RESEARCH.md) v1.1 (域调研 + 审计二轮记录)
- **Scope 调研**: [PHASE8_SCOPE_RESEARCH.md](../research/PHASE8_SCOPE_RESEARCH.md) v1.3 (D-1 裁决 + 双免费午餐支撑材料 §5.1)
- **前置 ADR**: [ADR-017](ADR_INDEX.md#adr-017-时序模块命名空间-cpphubv1timeseries) (timeseries 命名空间) / [ADR-018](ADR_INDEX.md#adr-018-slsqp-优化器实现边界) (K3 复用) / [ADR-019](ADR_INDEX.md#adr-019-v17-多变量时序与混频实施边界-263-项) (决策 1/20 本 ADR 兑现)
- **兼容性约束**: C1-C5 见调研报告第四部分 (零破坏 v1.7.1 API 冻结 / Eigen 隔离 / 对照库版本冻结 statsmodels 0.14.6 + R ARDL 0.2.5 / C ABI 不进 / 构建矩阵 cxx_std_20)
- **后续**: PHASE8_SPEC 编写 (M0) → M1 ARDL/PSS → M2 Kalman → M3 (D-2 待裁决)
