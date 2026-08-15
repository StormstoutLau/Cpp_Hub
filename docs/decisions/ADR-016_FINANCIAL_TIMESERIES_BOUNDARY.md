# ADR-016: 金融时间序列实施边界 (18 项)

**状态**: Accepted
**日期**: 2026-08-15
**版本归属**: v1.6 (Phase 7B M1/M2)
**关联 Phase**: 7B
**决策者**: 架构组
**调研依据**: [FINANCIAL_TIMESERIES_RESEARCH.md](../research/FINANCIAL_TIMESERIES_RESEARCH.md) v3.2 §13-14 (经 4 子 agent review 综合分析, 60 个幻觉点全部验证, G-ADR1/G-ADR7 决策依据已修正; 2026-08-15 二次审计修正 U-ADR1/5/7/8/10 决策依据与 spec 对齐)

---

## 背景

v1.6 GARCH 族与单位根检验实施中存在 18 个边界决策点 (GARCH 7 项 + 单位根 11 项), 需在实施前明确方案, 避免:
1. 数值不一致 (与 R `rugarch`/`arch`/`urca` 基准对照失败)
2. 理论偏差 (如 QMLE 协方差用 Hessian 逆而非 sandwich, 非正态下 SE 被低估)
3. 局部最优 (GARCH 似然非凸, 单一起始点易陷入次优解)
4. 临界值过时 (MacKinnon 1996 vs 2010 response surface 精度差异)

**决策原则**:
1. **开源库对齐**: 优先与 `arch`/`rugarch`/`urca` 对齐, 保证 R/Python 基准数值一致 (容差 1e-10 至 1e-12)
2. **理论正确性优先**: 理论正确性优先于实现便利 (如 sandwich 估计器比 Hessian 逆复杂但理论正确)
3. **预留扩展接口**: 为 v1.7 多元 GARCH (CCC/DCC) 与协整 (Johansen) 预留接口

---

## 决策

### A. GARCH 族边界决策 (7 项)

| # | 决策点 | 候选方案 | 决策 | 依据 |
|---|--------|---------|------|------|
| **G-ADR1** | 方差初始值策略 | (a) ε²₁ (b) 样本方差 (c) backcast (EWMA λ=0.94) | **(c) backcast** | `arch` 默认 (EWMA λ=0.94, τ=min(75,T)); 小样本稳健 (G1); 注: `rugarch` 默认样本方差, EWMA 为可选 |
| **G-ADR2** | 参数约束方法 | (a) Nelder-Mead 无约束 (b) SLSQP 约束优化 (c) log 重参数化 | **(b) SLSQP** | 支持不等式约束 (ω>0, α>0, β>0, α+β<1); `arch` 用 SLSQP, `rugarch` 用 solnp (G15/G22); v1.5 optimizer 需扩展 (约束 C1) |
| **G-ADR3** | QMLE 协方差估计 | (a) Hessian 逆 (b) sandwich H⁻¹SH⁻¹ (c) OPG | **(b) sandwich** | Bollerslev-Wooldridge 1992; 非正态下 H⁻¹ 低估 SE (G9); sandwich 仅需均值正确, 稳健性最强 |
| **G-ADR4** | 标准化残差 JB 检验临界值 | (a) 渐近 χ² (b) Bootstrap (c) 有限样本修正 | **(b) Bootstrap** | 小样本下 hₜ 有估计误差, 渐近临界值过度拒绝; 复用 v1.5 Bootstrap 实现 (G11) |
| **G-ADR5** | t 分布自由度估计 | (a) 两步矩估计 (b) 联合 QMLE | **(b) 联合 QMLE** | Bollerslev 1987; 两步法效率低且偏差大 (G14); 联合估计与 (ω,α,β) 同时优化 |
| **G-ADR6** | 优化器起始点策略 | (a) 单一起始 (b) 多起始取最优 | **(b) 多起始** | GARCH 似然非凸, 多起始避免局部最优 (G16); 默认 4 个起始点 (HR 高频/ML 极大似然/用户指定/随机扰动) |
| **G-ADR7** | GARCH-M 参数化 | (a) c·σₜ (b) c·σ²ₜ (c) c·log(σ²ₜ) | **(a) c·σₜ** | `arch` 默认 ('vol' 标准差形式, `form="vol"`); Engle-Lilien-Robbins 1987 原始参数化 (G18) |

### B. 单位根检验边界决策 (11 项)

| # | 决策点 | 候选方案 | 决策 | 依据 |
|---|--------|---------|------|------|
| **U-ADR1** | ADF lag 选择规则 | (a) 立方根 floor(T^1/3) (b) Schwert ceil(12(T/100)^0.25) (c) AIC/BIC | **(b) Schwert** | `arch` 默认; Schwert 1989; ceil (向上取整, 非 floor); 可选 AIC/BIC 作为备选 (U1) |
| **U-ADR2** | ADF 方程形式选择 | (a) 手动指定 nc/c/ct (b) 自动 (趋势检验后选择) | **(b) 自动** | `arch` 默认自动; 基于趋势显著性选择 c/ct; 减少用户误选 (U2) |
| **U-ADR3** | 临界值来源 | (a) MacKinnon 1996 (b) MacKinnon 2010 response surface (c) 查表 | **(b) MacKinnon 2010** | `arch` 已用 2010; 4 系数 response surface, 精度最高 (U3/U13) |
| **U-ADR4** | DF-GLS demean 临界值 | (a) ERS 1996 原表 (b) `arch` 独立系数表 | **(b) arch 独立表** | ERS 原表仅含 trend; demean 需 `arch` 扩展系数 (U7) |
| **U-ADR5** | DF-GLS c̄值 | (a) 0 (b) ERS 固定值 c̄=-7.0 (demean) / c̄=-13.5 (trend) | **(b) ERS 固定值** | ERS 1996; trend_spec=="c" => c̄=-7.0 (demean, 弱 detrending); trend_spec=="ct" => c̄=-13.5 (trend, 强 detrending); detrending 回归必须含 c̄ (U8); 与 arch `dfgls.py` 实现一致 |
| **U-ADR6** | KPSS H0 方向标注 | (a) 单位根 (错误) (b) 平稳性 | **(b) 平稳性** | KPSS 1992; API 必须明确标注 H0/H1 方向, 避免与 ADF 混淆 (U12) |
| **U-ADR7** | KPSS 长期方差带宽 | (a) 固定 (b) Hobijn et al. 1998 数据依赖法 (c) Schwert | **(b) Hobijn et al. 1998** | KPSS 1992; `arch` KPSS `_autolag` 默认 (基于 Bartlett 核假设); 带宽影响长期方差估计 (U11); legacy 模式 (lags=-1) 用 Schwert 规则; 核函数为 Bartlett (非 QS, arch `unitroot.py:1336,1352-1357` 实测) |
| **U-ADR8** | PP 检验带宽 | (a) 固定 (b) Schwert 规则 ceil(12(T/100)^0.25) (c) Andrews | **(b) Schwert 规则** | Phillips-Perron 1988; `arch` 默认 (与 ADF lag 相同公式, `unitroot.py:1135-1136` 实测); PP 用 Bartlett 核 (通过 cov_nw, 与 KPSS 相同; U5/U6) |
| **U-ADR9** | 基准对照软件 | (a) R urca (b) Python arch (c) R tseries | **(b) Python arch** | `arch` 用 MacKinnon 2010 (最新); `urca` 可能用旧版; 统一用 `arch` 作为 C++ 基准 |
| **U-ADR10** | 方差比检验变体 | (a) 仅 Z₁/Z₂ (b) Z₁/Z₂ + Chow-Denning (c) 全部含 debiased | **(c) 全部** | 覆盖 Lo-MacKinlay + Chow-Denning + CLM debiased (Campbell-Lo-MacKinlay 1997, 非 Chen-Deo 2006); 多重检验复用 Phase 7A `multiple_test_correction` (U14-U17); Chow-Denning 用 Z₂ (异方差稳健) + SMM(m,∞) 联合分布 |
| **U-ADR11** | 结构断点检验 | (a) v1.6 纳入 Zivot-Andrews (b) 推迟到 v1.7 | **(b) 推迟 v1.7** | v1.6 scope 已满; ZA 检验需独立临界值表; Phase 7A 已有 CUSUM/Andrews 通用框架 (U18) |

---

## 理由

### GARCH 族决策理由

1. **G-ADR1 backcast**: EWMA backcast 比样本方差更稳健, 对早期观测赋予递减权重, 避免 GARCH 启动期方差跳变; `arch` 默认策略, 保证基准一致
2. **G-ADR2 SLSQP**: GARCH 参数约束 (ω>0, α>0, β>0, α+β<1) 是必要条件 (违反则方差可能为负或不平稳), SLSQP 是处理不等式约束的标准方法; log 重参数化虽可避免约束, 但 Jacobian 复杂且数值不稳定
3. **G-ADR3 sandwich**: QMLE 在非正态分布下仍一致但协方差需用 sandwich 估计器 (H⁻¹SH⁻¹), Hessian 逆仅在正态分布下正确; 金融收益非正态是常态, 必须用 sandwich
4. **G-ADR4 Bootstrap**: 标准化残差 zₜ = εₜ/√hₜ 含估计误差, 渐近 χ² 临界值在小样本下过度拒绝 H₀; Bootstrap 通过重采样 zₜ 修正此偏差
5. **G-ADR5 联合 QMLE**: t 分布自由度 ν 与 (ω,α,β) 存在耦合, 两步法 (先估 GARCH 再估 ν) 忽略此耦合导致效率损失; 联合 QMLE 一次优化所有参数
6. **G-ADR6 多起始**: GARCH 对数似然函数非凸 (多局部最优), 单一起始点可能陷入次优解; 多起始策略提高全局最优命中率
7. **G-ADR7 c·σₜ**: Engle-Lilien-Robbins 1987 原始 GARCH-M 用标准差形式 (风险溢价线性于 σₜ), `arch` 默认; σ²ₜ 和 log(σ²ₜ) 是替代参数化, 但与原始理论不符

### 单位根检验决策理由

1. **U-ADR1 Schwert**: Schwert 1989 基于 ARMA 谱窗的 lag 选择规则 `ceil(12·(T/100)^0.25)` (向上取整, 非 floor; arch `unitroot.py:378-380, 1136, 1326` 三处实测一致), 对多种数据生成过程稳健; `arch` 默认, 保证基准一致
2. **U-ADR2 自动**: 用户可能误选方程形式 (如对明显有趋势的数据选 nc), 自动选择基于趋势显著性检验, 减少人为错误
3. **U-ADR3 MacKinnon 2010**: 2010 版 response surface 方法用 4 系数多项式逼近临界值, 精度远高于 1996 版查表; `arch` 已采用
4. **U-ADR4 arch 独立表**: ERS 1996 原表仅含 trend case, demean case 需 `arch` 扩展的独立系数; 用 ERS 原表会导致 demean 检验临界值错误
5. **U-ADR5 ERS 固定值**: c̄=-7.0 (demean, trend_spec=="c"), c̄=-13.5 (trend, trend_spec=="ct") 是 ERS 1996 推导的 detrending 回归常数, 省略会导致 GLS detrending 错误; 与 arch `dfgls.py` 实现一致
6. **U-ADR6 平稳性**: KPSS 的 H₀ 是平稳性 (与 ADF 相反), API 必须明确标注, 避免用户混淆
7. **U-ADR7 Hobijn et al. 1998 / U-ADR8 Schwert 规则**: KPSS 用 Hobijn et al. 1998 数据依赖法 (arch KPSS `_autolag` 默认, 基于 Bartlett 核假设); PP 用 Schwert 规则 `ceil(12·(T/100)^0.25)` (与 ADF lag 相同公式, arch `unitroot.py:1135-1136` 实测); 两者均用 Bartlett 核 (非 QS 核), 带宽影响长期方差估计
8. **U-ADR9 Python arch**: `arch` 用最新 MacKinnon 2010, `urca` 可能用旧版, 统一用 `arch` 保证基准最新
9. **U-ADR10 全部变体**: 方差比检验有多个变体 (Z₁/Z₂/Chow-Denning/CLM debiased), 全部实现覆盖所有用例; CLM debiased 是 Campbell-Lo-MacKinlay 1997 重叠块偏差修正 (arch `debiased=True` 实现, 非 Chen-Deo 2006); Chow-Denning 联合检验用 Z₂ (异方差稳健) + SMM(m,∞) 分布; 多重检验复用 Phase 7A 的 `multiple_test_correction`
10. **U-ADR11 推迟 ZA**: Zivot-Andrews 结构断点检验需独立临界值表, v1.6 scope 已满 (GARCH + 单位根 + 方差比); Phase 7A 已有 CUSUM/Andrews 通用结构断点框架, ZA 可在 v1.7 补充

---

## 替代方案评估

### G-ADR2 参数约束方法

| 方案 | 优点 | 缺点 | 结论 |
|------|------|------|------|
| (a) Nelder-Mead 无约束 | 实现简单 | 无法处理约束, 参数可能违反平稳性 | ❌ 拒绝 |
| **(b) SLSQP** | 标准约束优化, `arch` 对齐 | 需扩展 v1.5 optimizer | ✅ 采纳 |
| (c) log 重参数化 | 无需约束优化器 | Jacobian 复杂, 数值不稳定 | ❌ 拒绝 |

### U-ADR3 临界值来源

| 方案 | 优点 | 缺点 | 结论 |
|------|------|------|------|
| (a) MacKinnon 1996 | 经典 | 精度低, 查表插值误差 | ❌ 拒绝 |
| **(b) MacKinnon 2010** | 精度最高, `arch` 对齐 | 需实现 response surface | ✅ 采纳 |
| (c) 查表 | 实现简单 | 精度低, 维护成本高 | ❌ 拒绝 |

---

## 后果

### 实施约束

1. **C1 优化器扩展** (前置): v1.5 `optimizer.hpp` 需扩展 SLSQP 约束优化 (G-ADR2), 阶段 1 实施前完成
2. **C2 Estimator 接口**: `Estimator` 基类已支持 `TimeSeriesData` (variant 第三分支), 仅 MLE/OLS 派生类 visit 逻辑需补齐
3. **C4 命名空间**: 时序模型统一落 `cpphub::v1::timeseries` (见 ADR-017)
4. **C6 数据载体**: `EconData` 已含 `TimeSeriesData` (单变量), 仅缺 `MultivariateTSData` (多变量, 阶段 5 前置)

### 测试基准

- **基准软件**: Python `arch` (U-ADR9), R `rugarch` (GARCH 族交叉验证)
- **容差**: 1e-10 至 1e-12 (与 v1.4 HFE 模块一致)
- **核查模式**: 参照 v1.4.0/1.4.1 R `highfrequency` 核查模式, 编写 `verify_X.R`/`verify_X.py` 脚本, 打印源码 + 小样本手算

### 幻觉点核查

18 项决策点对应的幻觉点 (G1-G23, U1-U22) 必须在实施时逐点核查:
- **GARCH 族**: G1 (backcast), G7 (EGARCH 非对称项), G9 (QMLE sandwich), G11 (JB Bootstrap), G14 (t 分布联合估计), G15/G22 (SLSQP vs solnp), G16 (多起始), G18 (GARCH-M 参数化)
- **单位根**: U1 (Schwert lag, ceil 非 floor), U2 (自动方程形式), U3/U13 (MacKinnon 2010, 4 系数 3 次多项式), U5 (PP 带宽用 Schwert 规则, 非 NW automatic), U6 (PP 用 Bartlett 核, 与 KPSS 相同), U7 (DF-GLS demean), U8 (c̄值, c→-7.0/ct→-13.5), U9 (GLS detrending 用 ρ̄=1+c̄/T), U11 (KPSS 带宽用 Hobijn et al. 1998, 非 Andrews; 核用 Bartlett, 非 QS), U12 (KPSS H0 方向), U14-vr (VR(k) 不再除 k), U15 (Z₂ 三重: 因子 4/4 阶矩 δⱼ/√T), U16 (Chow-Denning 用 Z₂ + SMM 分布), U17 (CLM debiased, 非 Chen-Deo 2006), U18 (ZA 推迟)

---

## 关联

- **调研报告**: [FINANCIAL_TIMESERIES_RESEARCH.md](../research/FINANCIAL_TIMESERIES_RESEARCH.md) v3.2 §13-14
- **前置 ADR**: [ADR-013](ADR_INDEX.md#adr-013-双层线性代数架构-固定尺寸--动态尺寸) (Eigen3 隔离边界), [ADR-014](ADR_INDEX.md#adr-014-标定-calibration-vs-估计-estimation-的分离) (calibration vs estimation 分离)
- **关联 ADR**: [ADR-017](ADR-017_TIMESERIES_NAMESPACE.md) (时序模块命名空间, 约束 C4)
- **后续 ADR**: 待编写 ADR-018 (SLSQP optimizer 扩展, 约束 C1), 若需要
