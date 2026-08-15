# Phase 7B 审计验收 Checklist

> **版本归属**: v1.6 (Phase 7B M1/M2)
> **关联 Spec**: [PHASE7B_FINANCIAL_TS_SPEC.md](./PHASE7B_FINANCIAL_TS_SPEC.md)
> **关联 ADR**: ADR-013 / ADR-016 / ADR-017 / ADR-018
> **验收日期**: ___________
> **验收人**: ___________

---

## 使用说明

- 每项标注 `[ ]` 待验收，`[x]` 已通过，`[!]` 未通过（需附 issue 描述）
- **G (Garch)** = M1 GARCH 族幻觉点 (G1-G23)
- **U (Unit root)** = M2 单位根与方差比幻觉点 (U1-U22)
- 容差标注: 1e-10 = 相对误差 ≤ 1e-10; "精确" = 数学定义完全一致
- 基准库: Python `arch` (主基准, U-ADR9) + R `rugarch` (GARCH 交叉验证) + R `vrtest` (Chow-Denning)

---

## 1. 交付物完整性

### 1.1 M1 GARCH 族头文件 (6 个)

| # | 文件 | 状态 | 备注 |
|---|------|------|------|
| 1.1.1 | `include/cpphub/timeseries/garch/garch_model.hpp` | [ ] | GARCH(1,1) QMLE + backcast + 递归方差 |
| 1.1.2 | `include/cpphub/timeseries/garch/egarch_model.hpp` | [ ] | EGARCH QMLE + log 方差 + 杠杆效应 |
| 1.1.3 | `include/cpphub/timeseries/garch/gjr_garch_model.hpp` | [ ] | GJR-GARCH QMLE + 非对称项 |
| 1.1.4 | `include/cpphub/timeseries/garch/garch_distribution.hpp` | [ ] | Normal/t/GED 似然函数 |
| 1.1.5 | `include/cpphub/timeseries/garch/garch_forecast.hpp` | [ ] | 多步方差预测 |
| 1.1.6 | `include/cpphub/timeseries/garch/garch_diagnostics.hpp` | [ ] | 标准化残差诊断 (复用 Phase 7A) |

### 1.2 M2 单位根与方差比头文件 (7 个)

| # | 文件 | 状态 | 备注 |
|---|------|------|------|
| 1.2.1 | `include/cpphub/timeseries/unit_root/adf_test.hpp` | [ ] | ADF 检验 |
| 1.2.2 | `include/cpphub/timeseries/unit_root/pp_test.hpp` | [ ] | PP 检验 |
| 1.2.3 | `include/cpphub/timeseries/unit_root/kpss_test.hpp` | [ ] | KPSS 检验 |
| 1.2.4 | `include/cpphub/timeseries/unit_root/df_gls_test.hpp` | [ ] | DF-GLS 检验 |
| 1.2.5 | `include/cpphub/timeseries/unit_root/variance_ratio_test.hpp` | [ ] | 方差比检验 (4 变体) |
| 1.2.6 | `include/cpphub/timeseries/unit_root/mackinnon_cv.hpp` | [ ] | MacKinnon 2010 临界值 |
| 1.2.7 | `include/cpphub/timeseries/unit_root/unit_root_common.hpp` | [ ] | 共享工具 (Schwert/长期方差) |

### 1.3 基准验证脚本 (17 个, 已提交版本控制)

> **修正说明 (2026-08-15)**: 原计划 8 个 verify 脚本不入版本控制; 实施时按 v1.4.1 可追溯惯例改为**全部提交** (`tests/fixtures/timeseries/`, commit `1441fbb`)。`.gitignore` 排除的仅限含版权的第三方源码副本。单位根 5 检验 (ADF/PP/KPSS/DF-GLS/VR) 合并为单一 `verify_unit_root.py`; 另提交 12 个 probe 脚本 (arch 源码逐点核查过程可追溯)。

| # | 脚本 | 状态 | 备注 |
|---|------|------|------|
| 1.3.1 | `tests/fixtures/timeseries/verify_garch.py` | [ ] | arch GARCH(1,1) 对照 |
| 1.3.2 | `tests/fixtures/timeseries/verify_egarch.py` | [ ] | arch EGARCH 对照 |
| 1.3.3 | `tests/fixtures/timeseries/verify_gjr.py` | [ ] | arch GJR 对照 |
| 1.3.4 | `tests/fixtures/timeseries/verify_unit_root.py` | [ ] | arch 单位根 5 检验合并对照 (ADF/PP/KPSS/DF-GLS/VR) |
| 1.3.5 | `tests/fixtures/timeseries/gen_mackinnon_tables.py` | [ ] | MacKinnon 2010 系数表 arch 源码 → C++ 生成器 |
| 1.3.6 | `probe_arch_unitroot(1-6).py` 等 12 个 probe | [ ] | arch 源码逐点核查 (convention/cov/egarch/unitroot×6/egarch_sim/egarch_src/gjr_h1) |

### 1.4 测试套件 (13 套, 203 用例)

> **修正说明 (2026-08-15)**: 原列 12 套 191 例; 实际交付 13 套 203 例 — 补 `test_unit_root_common` (12 例), 且 garch_distribution 13 / variance_ratio 22 / mackinnon_cv 12 / integration 5 与原估计有出入。M1 = 84, M2 = 114, 集成 = 5。

| # | 测试套件 | 用例数 | 状态 | 备注 |
|---|----------|--------|------|------|
| 1.4.1 | `test_garch_model` | 20 | [ ] | M1 |
| 1.4.2 | `test_egarch_model` | 18 | [ ] | M1 |
| 1.4.3 | `test_gjr_garch_model` | 18 | [ ] | M1 |
| 1.4.4 | `test_garch_distribution` | 13 | [ ] | M1 |
| 1.4.5 | `test_garch_diagnostics` | 15 | [ ] | M1 |
| 1.4.6 | `test_adf_test` | 20 | [ ] | M2 |
| 1.4.7 | `test_pp_test` | 15 | [ ] | M2 |
| 1.4.8 | `test_kpss_test` | 15 | [ ] | M2 |
| 1.4.9 | `test_df_gls_test` | 18 | [ ] | M2 |
| 1.4.10 | `test_variance_ratio_test` | 22 | [ ] | M2 |
| 1.4.11 | `test_mackinnon_cv` | 12 | [ ] | M2 |
| 1.4.12 | `test_unit_root_common` | 12 | [ ] | M2 (原清单遗漏) |
| 1.4.13 | `test_integration_phase7b` | 5 | [ ] | 端到端集成 (spec §4 即 5 场景, 原列 10 有误) |

---

## 2. 编译与跨平台测试

> **实测数据 (2026-08-15, 已完成)**: 主控 MSVC **2207/2207** (614.21 sec) / A 站 GCC **2189/2189** (364.47 sec) / B 站 GCC **2189/2189** (358.56 sec), 零失败; GitHub Actions CI run #47 (commit `a1b7215`) 4/4 job 全绿 (仓库首个)。MSVC 与 GCC 18 个测试差额为平台专属用例, 非功能差异。A/B 站当日 github 阻断, 采用 bundle 中继 (git bundle + googletest/eigen 源 tar scp + FETCHCONTENT_SOURCE_DIR 本地覆盖)。状态列留白待正式验收签署。

### 2.1 主控站 (Windows MSVC Release)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 2.1.1 | CMake 配置成功 (`cmake -DCPPHUB_ENABLE_C_API=ON`) | [ ] | 实测通过 (2026-08-15) |
| 2.1.2 | MSVC Release 编译零警告零错误 | [ ] | 实测通过 (MSVC 19.50) |
| 2.1.3 | 全量 ctest 通过 (新增 ~201 + 现有 = ~2205) | [ ] | 实测 **2207/2207** (614.21 sec) |
| 2.1.4 | 无现有测试退化 (Phase 1-7A 全部仍通过) | [ ] | 实测零回归 |
| 2.1.5 | SLSQP 12/12 测试仍通过 (ADR-018 无退化) | [ ] | 实测通过 (Calib 套件 12/12) |

### 2.2 A 站 (Ubuntu GCC)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 2.2.1 | fresh clone + rebuild 成功 | [ ] | 实测通过 (bundle 中继 + eigen tar, GCC 13.3) |
| 2.2.2 | GCC 编译零警告零错误 | [ ] | 实测通过 |
| 2.2.3 | ctest 全量通过 | [ ] | 实测 **2189/2189** (364.47 sec) |
| 2.2.4 | 与主控站数值一致 (容差 1e-10) | [ ] | 实测一致 (硬编码基准三平台同源) |

### 2.3 B 站 (Ubuntu GCC)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 2.3.1 | fresh clone + rebuild 成功 | [ ] | 实测通过 (bundle 中继 + eigen tar, GCC 13.3) |
| 2.3.2 | GCC 编译零警告零错误 | [ ] | 实测通过 |
| 2.3.3 | ctest 全量通过 | [ ] | 实测 **2189/2189** (358.56 sec) |
| 2.3.4 | 与主控站数值一致 (容差 1e-10) | [ ] | 实测一致 (硬编码基准三平台同源) |

### 2.4 三平台一致性

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 2.4.1 | SLSQP + GARCH + 单位根测试三平台无数值偏差 | [ ] | 实测一致 (含 CI run #47 双平台全绿) |
| 2.4.2 | 无平台相关编译宏差异 (如 `timegm`/`_mkgmtime`) | [ ] | 实测无新增平台宏 (7B 全 header-only, 纯 std) |

---

## 3. M1 GARCH 族数值基准对照

### 3.1 GARCH(1,1) — vs Python `arch`

| # | 验证点 | 容差 | 状态 | arch 配置 |
|---|--------|------|------|-----------|
| 3.1.1 | ω/α/β 参数估计 | 1e-10 | [ ] | `arch_model(o='Constant', vol='GARCH')` |
| 3.1.2 | backcast 方差 σ²₀ (λ=0.94) | 1e-12 | [ ] | 前 τ 个残差, 有限样本归一化 |
| 3.1.3 | 条件方差序列 hₜ | 1e-10 | [ ] | |
| 3.1.4 | sandwich 协方差 V (robust) | 1e-10 | [ ] | `cov_type='robust'` |
| 3.1.5 | 标准误 SE = √diag(V) | 1e-10 | [ ] | |
| 3.1.6 | 多步预测 h_{T+k} (horizon=10) | 1e-10 | [ ] | |
| 3.1.7 | 预测收敛到 σ̄²=ω/(1-α-β) | 1e-6 | [ ] | k→∞ |
| 3.1.8 | 约束违反检测 (ω<0 拒绝) | - | [ ] | |
| 3.1.9 | 含 NaN 数据拒绝 | - | [ ] | |
| 3.1.10 | T<10 拒绝 | - | [ ] | |
| 3.1.11 | 大样本 T=5000 性能 < 5 sec | - | [ ] | |

### 3.2 GARCH(1,1) — vs R `rugarch`

| # | 验证点 | 容差 | 状态 | rugarch 配置 |
|---|--------|------|------|-------------|
| 3.2.1 | ω/α/β 参数 | 1e-8 | [ ] | `ugarchfit(model="sGARCH")` |
| 3.2.2 | solver 差异 (SLSQP vs solnp) | 1e-8 | [ ] | G22 跨库核查 |

### 3.3 t-GARCH(1,1)

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 3.3.1 | ν 联合估计 | 1e-8 | [ ] | `dist='t'`, G14 联合 QMLE |
| 3.3.2 | ν 边界 [2.05, 500.0] | 精确 | [ ] | arch StudentsT 默认 |
| 3.3.3 | ν→1000 退化到 Normal | 1e-6 | [ ] | 渐近行为 |

### 3.4 GED-GARCH(1,1)

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 3.4.1 | GED ν=2 退化为正态 | 1e-10 | [ ] | |
| 3.4.2 | GED 对数似然公式 | 1e-15 | [ ] | G-GED1: 指数 ν, 缩放 c 含 2^(-1/ν) |
| 3.4.3 | 批量似然求和 = 逐项相加 | 1e-15 | [ ] | |

### 3.5 GARCH 似然与 AIC/BIC

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 3.5.1 | Normal 似然含 -0.5·log(2π) 常数项 | 1e-15 | [ ] | G3 |
| 3.5.2 | AIC/BIC 用完整似然 (含常数项) | 1e-10 | [ ] | G17 |
| 3.5.3 | AIC/BIC vs arch | 1e-10 | [ ] | |

### 3.6 EGARCH — vs Python `arch`

| # | 验证点 | 容差 | 状态 | arch 配置 |
|---|--------|------|------|-----------|
| 3.6.1 | ω/α/β/γ 参数 | 1e-10 | [ ] | `vol='EGARCH'` |
| 3.6.2 | E\|z\| = √(2/π) ≈ 0.7979 | 1e-15 | [ ] | G5 |
| 3.6.3 | 非对称项用 zₜ₋₁ (非 εₜ₋₁) | 1e-10 | [ ] | G6 |
| 3.6.4 | 条件方差序列 hₜ | 1e-10 | [ ] | |
| 3.6.5 | sandwich 协方差 V | 1e-10 | [ ] | `cov_type='robust'` |
| 3.6.6 | 参数映射 (spec α=arch γ, spec γ=arch α) | 1e-10 | [ ] | G23 符号约定 |
| 3.6.7 | t-EGARCH ν 联合估计 | 1e-8 | [ ] | |
| 3.6.8 | 平稳性 \|β\|<1 检测 | - | [ ] | |
| 3.6.9 | 多步预测 (simulation 对照) | 1e-6 | [ ] | arch 仅 horizon=1 analytic |
| 3.6.10 | 杠杆效应 (γ<0 负冲击放大) | 方向性 | [ ] | |

### 3.7 GJR-GARCH — vs Python `arch`

| # | 验证点 | 容差 | 状态 | arch 配置 |
|---|--------|------|------|-----------|
| 3.7.1 | ω/α/γ/β 参数 | 1e-10 | [ ] | `o=1` |
| 3.7.2 | 平稳性 α+γ/2+β<1 | 1e-12 | [ ] | G10, arch 约束矩阵 -0.5 |
| 3.7.3 | γ 边界 [-1, 1] (允许负值) | 精确 | [ ] | G-gamma-sign |
| 3.7.4 | I(z<0) ≡ I(ε<0) | 1e-15 | [ ] | G7 |
| 3.7.5 | 条件方差序列 hₜ | 1e-10 | [ ] | |
| 3.7.6 | sandwich 协方差 V | 1e-10 | [ ] | |
| 3.7.7 | t-GJR ν 联合估计 | 1e-8 | [ ] | |
| 3.7.8 | 杠杆效应 (γ>0 负冲击放大) | 方向性 | [ ] | G8 |

### 3.8 GARCH 标准化残差诊断

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 3.8.1 | zₜ = εₜ/√hₜ (非 εₜ/σₜ) | 1e-15 | [ ] | G11 |
| 3.8.2 | JB 检验 (Bootstrap 临界值) | 1e-6 | [ ] | G-ADR4 |
| 3.8.3 | LB zₜ vs R `Box.test(z, lag=10)` | 1e-8 | [ ] | |
| 3.8.4 | LB z²ₜ vs R `Box.test(z², lag=10)` | 1e-8 | [ ] | G12 ARCH 效应 |
| 3.8.5 | LB 滞后自动选择 (floor(log(T)) 或 Schwert) | - | [ ] | G12 |

### 3.9 GARCH(1,1) 多起始点

| # | 验证点 | 状态 | 备注 |
|---|--------|------|------|
| 3.9.1 | 4 起始点策略 (HR/ML/user/random) | [ ] | G16, G-ADR6 |
| 3.9.2 | 取对数似然最大者 | [ ] | |
| 3.9.3 | 避免局部最优 (非凸似然) | [ ] | |

---

## 4. M2 单位根检验数值基准对照

### 4.1 ADF 检验 — vs Python `arch.unitroot.ADF`

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 4.1.1 | ADF τ 统计量 | 1e-10 | [ ] | |
| 4.1.2 | Schwert lag `ceil(12·(T/100)^0.25)` | 1e-10 | [ ] | U1, 向上取整 |
| 4.1.3 | 自动方程形式 (nc/c/ct) | 1e-10 | [ ] | U2 |
| 4.1.4 | ADF nc 方程 (无常数无趋势) | 1e-10 | [ ] | |
| 4.1.5 | ADF c 方程 (常数) | 1e-10 | [ ] | |
| 4.1.6 | ADF ct 方程 (常数+趋势) | 1e-10 | [ ] | |
| 4.1.7 | AIC/BIC lag 选择 (备选) | 1e-10 | [ ] | |
| 4.1.8 | τ 非标准分布 (非 Student-t) | - | [ ] | U4 |
| 4.1.9 | p 值 (MacKinnon 2010) | 1e-12 | [ ] | U3/U13 |

### 4.2 PP 检验 — vs Python `arch.unitroot.PhillipsPerron`

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 4.2.1 | PP Z(tau) 统计量 | 1e-10 | [ ] | |
| 4.2.2 | 带宽 = Schwert `ceil(12·(T/100)^0.25)` | 精确 | [ ] | U5, 非 NW automatic |
| 4.2.3 | Bartlett 核 (通过 cov_nw) | - | [ ] | U6, 与 KPSS 相同 |
| 4.2.4 | 短期方差 s² = SSR/(n-k) (df-corrected) | 1e-10 | [ ] | U5-sigma2 |
| 4.2.5 | Z(tau) 分母 = 2·σ·s (非 2·σ², 非 2·σ·σ_ε) | 精确 | [ ] | 三次修正 |
| 4.2.6 | 共享 MacKinnon 2010 临界值 | 1e-12 | [ ] | 与 ADF 相同 |

### 4.3 KPSS 检验 — vs Python `arch.unitroot.KPSS`

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 4.3.1 | KPSS LM 统计量 | 1e-10 | [ ] | |
| 4.3.2 | H0: 平稳性 (与 ADF 相反) | - | [ ] | U12, API 明确标注 |
| 4.3.3 | Bartlett 核 (非 QS 核) | - | [ ] | U11-kernel |
| 4.3.4 | Hobijn et al. 1998 带宽 (默认) | 1e-10 | [ ] | U11, 非 Andrews |
| 4.3.5 | legacy Schwert 带宽 (备选) | 1e-10 | [ ] | lags=-1 |
| 4.3.6 | "c" (level stationary) | 1e-10 | [ ] | |
| 4.3.7 | "ct" (trend stationary) | 1e-10 | [ ] | |
| 4.3.8 | 临界值 (KPSS 1992 原表/MC) | 1e-6 | [ ] | U10 非标准分布 |

### 4.4 DF-GLS 检验 — vs Python `arch.unitroot.DFGLS`

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 4.4.1 | DF-GLS τ 统计量 | 1e-10 | [ ] | |
| 4.4.2 | c̄=-7.0 (trend_spec=="c", demean) | 精确 | [ ] | U8 |
| 4.4.3 | c̄=-13.5 (trend_spec=="ct", trend) | 精确 | [ ] | U8 |
| 4.4.4 | ρ̄ = 1 + c̄/T (非 c̄ 本身) | 1e-10 | [ ] | U9-rho |
| 4.4.5 | GLS 变换: y* = [y₁, y₂-ρ̄·y₁, ...] | 1e-10 | [ ] | U9, 第一项不变换 |
| 4.4.6 | lag 选择 = AIC (非 MAIC) | 1e-10 | [ ] | U20, arch 默认 method="aic" |
| 4.4.7 | 临界值 (arch 独立模拟, 非 ERS 原表) | 1e-10 | [ ] | U7 |
| 4.4.8 | p 值 (MacKinnon 1994/2010) | 1e-12 | [ ] | U21 |
| 4.4.9 | "c" (demean) | 1e-10 | [ ] | |
| 4.4.10 | "ct" (trend) | 1e-10 | [ ] | |

### 4.5 MacKinnon 2010 临界值

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 4.5.1 | 4 系数 3 次多项式 CV(p)=β_∞+β_0/T+β_1/T²+β_2/T³ | 1e-12 | [ ] | U3, 非 5 系数 4 次 |
| 4.5.2 | ADF nc/c/ct 三套独立系数 | 1e-12 | [ ] | 不能混用 |
| 4.5.3 | PP 共享 ADF 系数 | 1e-12 | [ ] | |
| 4.5.4 | DF-GLS c/ct 独立系数 | 1e-12 | [ ] | U7 |
| 4.5.5 | p 值 response surface | 1e-12 | [ ] | |

---

## 5. M2 方差比检验数值基准对照

### 5.1 Lo-MacKinlay Z₁/Z₂ — vs Python `arch.unitroot.VarianceRatio`

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 5.1.1 | VR(k) = σ̂²_k / σ̂²_1 (不再除 k) | 1e-10 | [ ] | U14-vr |
| 5.1.2 | σ̂²_k = (1/(Tk))·Σ(Rₜ^(k)-k·μ̂)² (含 1/k, 减 k·μ̂) | 1e-10 | [ ] | U17-demean |
| 5.1.3 | 默认 trend="c" 强制去均值 | - | [ ] | U17-demean |
| 5.1.4 | μ̂ = 样本平均收益 (trend="n" 时 μ̂=0) | 1e-10 | [ ] | |
| 5.1.5 | Z₁ = [VR(k)-1] / √[2(2k-1)(k-1)/(3kT)] | 1e-10 | [ ] | U14, 1/T 放入方差 |
| 5.1.6 | Z₁ p 值 = 2·(1-Φ(\|Z₁\|)) | 1e-10 | [ ] | |
| 5.1.7 | Z₂ δⱼ = T·Σ(rₜ-μ̂)²·(rₜ₋ⱼ-μ̂)² / [Σ(rₜ-μ̂)²]² (4 阶矩) | 1e-10 | [ ] | U15(b) |
| 5.1.8 | Z₂ Var = 4·Σ((k-j)/k)²·δⱼ (因子 4, 权重有平方) | 1e-10 | [ ] | U15(a) |
| 5.1.9 | Z₂ = √T·[VR(k)-1] / √θ (含 √T 因子) | 1e-10 | [ ] | U15(c) |
| 5.1.10 | Z₂ p 值 = 2·(1-Φ(\|Z₂\|)) | 1e-10 | [ ] | |
| 5.1.11 | δⱼ 求和范围 j=1 到 k-1 | 1e-10 | [ ] | |

### 5.2 CLM debiased — vs Python `arch` (debiased=True)

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 5.2.1 | σ̂²_1 ← σ̂²_1·T/(T-1) (Bessel 修正) | 1e-10 | [ ] | U17 |
| 5.2.2 | σ̂²_k ← σ̂²_k·T²/[(T-k+1)(T-k)] (重叠块修正) | 1e-10 | [ ] | U17 |
| 5.2.3 | VR_debiased = VR·T(T-1)/[(T-k+1)(T-k)] | 1e-10 | [ ] | |
| 5.2.4 | 归属: CLM 1997 (非 Chen-Deo 2006) | - | [ ] | U17 |

### 5.3 Chow-Denning 联合检验 — vs R `vrtest::Chow.Denning`

| # | 验证点 | 容差 | 状态 | 备注 |
|---|--------|------|------|------|
| 5.3.1 | CD = max_k \|Z₂(k)\| (用 Z₂, 非 Z₁) | 1e-8 | [ ] | U16 |
| 5.3.2 | 联合 p 值 = 1-[2·Φ(\|CD\|)-1]^m (SMM(m,∞) 分布) | 1e-8 | [ ] | U16-smm, m=\|k_list\| |
| 5.3.3 | 多 horizon {2, 5, 10, 20} | 1e-8 | [ ] | |
| 5.3.4 | 非 Bonferroni 修正 | - | [ ] | |

---

## 6. 幻觉点逐项核查 (verify 脚本)

### 6.1 GARCH 族幻觉点 (G1-G23)

| ID | 影响级 | 核查内容 | 容差 | 状态 | verify 脚本 |
|----|--------|----------|------|------|-------------|
| G1 | 高 | backcast 用前 τ 个残差, 有限样本归一化 | 1e-12 | [ ] | verify_garch.py |
| G1-dir | 高 | arch 用前 τ 个残差 (非末尾) | 1e-12 | [ ] | |
| G1-norm | 高 | 权重 wᵢ=λⁱ/Σλʲ (非无穷级数近似) | 1e-12 | [ ] | |
| G2 | 中 | ε₀ = μ̂ (非 0) | 1e-10 | [ ] | |
| G3 | 中 | 似然含 -0.5·log(2π) | 1e-15 | [ ] | |
| G4 | 高 | 参数约束 ω>0, α≥0, β≥0, α+β<1 | 1e-12 | [ ] | |
| G5 | 高 | E\|z\| = √(2/π) | 1e-15 | [ ] | verify_egarch.py |
| G6 | 高 | 非对称项用 zₜ₋₁ (非 εₜ₋₁) | 1e-10 | [ ] | |
| G7 | 低 | I(z<0) ≡ I(ε<0) | 1e-15 | [ ] | |
| G8 | 高 | GJR 符号 (γ>0 杠杆) | 方向性 | [ ] | verify_gjr.py |
| G9 | 高 | sandwich V = H⁻¹SH⁻¹ (非 Hessian 逆) | 1e-10 | [ ] | |
| G10 | 高 | 平稳性 α+γ/2+β<1 (非 α+γ+β<1) | 1e-12 | [ ] | |
| G11 | 中 | zₜ = εₜ/√hₜ | 1e-15 | [ ] | |
| G12 | 中 | z² LB 滞后 = floor(log(T)) 或 Schwert | 1e-8 | [ ] | |
| G13 | 中 | 预测指数 k-1, 起点 h_{T+1} | 1e-6 | [ ] | |
| G14 | 高 | t 分布 ν 联合 QMLE | 1e-8 | [ ] | |
| G15 | 高 | SLSQP 优化器 (ADR-018) | 1e-10 | [ ] | |
| G16 | 中 | 多起始点 (4 个) | - | [ ] | |
| G17 | 中 | AIC/BIC 用完整似然 | 1e-10 | [ ] | |
| G-gamma-sign | 中 | GJR γ 边界 [-1,1] (允许负值) | 精确 | [ ] | |
| G-nu-bound | 中 | ν 边界 [2.05, 500.0] (三模型一致) | 精确 | [ ] | |
| G-GED1 | 高 | GED 指数 ν, 缩放 c 含 2^(-1/ν) | 1e-15 | [ ] | |
| G23 | 高 | EGARCH 参数映射 (spec α=arch γ) | 1e-10 | [ ] | |
| G22 | 中 | 跨库 solver 差异 (arch SLSQP vs rugarch solnp) | 1e-8 | [ ] | |

**推迟项** (不在本次验收范围):

| ID | 状态 | 备注 |
|----|------|------|
| G18 | [ ] | GARCH-M (推迟 v1.6+) |
| G19 | [ ] | APARCH δ (推迟 v1.6+) |
| G20 | [ ] | FIGARCH 截断 (推迟 v1.6+) |
| G21 | [ ] | IGARCH 非平稳 (推迟 v1.6+) |

### 6.2 单位根与方差比幻觉点 (U1-U22)

| ID | 影响级 | 核查内容 | 容差 | 状态 | verify 脚本 |
|----|--------|----------|------|------|-------------|
| U1 | 高 | Schwert `ceil(12·(T/100)^0.25)` (非 floor) | 1e-10 | [ ] | verify_adf.py |
| U1-round | 高 | arch 3 处实现 (ADF/PP/KPSS legacy) 全部 ceil | 精确 | [ ] | |
| U2 | 高 | 自动方程形式 (趋势显著性) | 1e-10 | [ ] | |
| U3 | 极高 | MacKinnon 2010: 4 系数 3 次多项式 (非 5 系数 4 次) | 1e-12 | [ ] | verify_adf.py |
| U3-coef | 极高 | arch tau_2010 每行 4 元素 | 精确 | [ ] | |
| U4 | 高 | ADF τ 非标准分布 (非 Student-t) | - | [ ] | |
| U5 | 高 | PP 带宽 = Schwert (非 floor(4·(T/100)^(2/9))) | 精确 | [ ] | verify_pp.py |
| U5-sigma2 | 中 | PP 短期方差 s²=SSR/(n-k) (df-corrected, 非 SSR/n) | 1e-10 | [ ] | |
| U6 | 高 | PP 用 Bartlett 核 (与 KPSS 相同) | - | [ ] | |
| U7 | 极高 | DF-GLS 临界值 arch 独立模拟 (非 ERS 原表) | 1e-10 | [ ] | verify_df_gls.py |
| U8 | 高 | c̄=-7.0 (c), c̄=-13.5 (ct) | 精确 | [ ] | |
| U9 | 高 | GLS detrending (非 OLS) | 1e-10 | [ ] | |
| U9-rho | 高 | ρ̄ = 1+c̄/T (非 c̄ 本身) | 1e-10 | [ ] | |
| U10 | 高 | KPSS 非标准分布 (非 χ²) | 1e-6 | [ ] | verify_kpss.py |
| U11 | 高 | KPSS 带宽 Hobijn et al. 1998 (非 Andrews) | 1e-10 | [ ] | |
| U11-kernel | 高 | KPSS 用 Bartlett 核 (非 QS 核) | - | [ ] | |
| U12 | 极高 | KPSS H0: 平稳性 (与 ADF 相反) | - | [ ] | |
| U13 | 中 | MacKinnon 2010 (非 1996 旧版) | 1e-12 | [ ] | |
| U14 | 高 | Z₁ = [VR-1]/√[2(2k-1)(k-1)/(3kT)] | 1e-10 | [ ] | verify_variance_ratio.py |
| U14-vr | 极高 | VR(k) = σ̂²_k/σ̂²_1 (不再除 k) | 1e-10 | [ ] | |
| U15 | 极高 | Z₂ 三重: 因子 4 / 4 阶矩 δⱼ / √T | 1e-10 | [ ] | |
| U16 | 高 | Chow-Denning 用 Z₂ + SMM(m,∞) 分布 | 1e-8 | [ ] | |
| U16-smm | 高 | p = 1-[2·Φ(\|CD\|)-1]^m (非 Bonferroni) | 1e-8 | [ ] | |
| U17 | 高 | CLM debiased (非 Chen-Deo 2006) | 1e-8 | [ ] | |
| U17-demean | 高 | 默认 trend="c" 强制去均值, σ̂²_k 减 k·μ̂ | 1e-10 | [ ] | |
| U18 | 中 | 多检验多重修正 (复用 Phase 7A) | 1e-10 | [ ] | |
| U20 | 高 | DF-GLS lag 用 AIC (非 MAIC) | 1e-10 | [ ] | |
| U20-method | 高 | arch 默认 method="aic"; MAIC 仅用于 M 检验族 | 精确 | [ ] | |
| U21 | 中 | MacKinnon 1994 vs 2010 版本区分 | - | [ ] | |

**推迟项** (不在本次验收范围):

| ID | 状态 | 备注 |
|----|------|------|
| U19 | [ ] | Ng-Perron M 检验族 (推迟 v1.7) |
| U22 | [ ] | Zivot-Andrews 结构断点 (推迟 v1.7) |

---

## 7. 端到端集成测试

> **修正说明 (2026-08-15)**: 实际交付为 **5 个 TEST 场景** (与 spec §4 表格 1:1, 即 7.1-7.5); 7.6-7.10 为复用/方向性验证点, 内嵌于 5 场景断言或由对应单测套件覆盖 (Kupiec/MZ/DM/BH 在场景 1/3/4 内, t 分布残差与 EGARCH 杠杆在 test_garch_diagnostics / test_egarch_model 内)。

| # | 场景 | 验证点 | 状态 | 备注 |
|---|------|--------|------|------|
| 7.1 | GARCH→VaR 集成 | GARCH 波动率 → VaR → Backtest 通路完整 | [ ] | |
| 7.2 | ADF→伪回归诊断 | PE/PB 因子值 ADF, 单位根则差分 | [ ] | |
| 7.3 | GARCH vs HAR 对比 | 同一收益率, GARCH vs HAR 预测精度 | [ ] | MZ/DM 检验复用 |
| 7.4 | 多检验多重修正 | ADF+PP+KPSS+DF-GLS + BH 修正 | [ ] | U18 |
| 7.5 | GARCH 标准化残差诊断 | JB+LB+z²LB 全流程 | [ ] | G11/G12 |
| 7.6 | GARCH-VaR Backtest 复用 Phase 7A | risk 模块 VaR 计算 | [ ] | |
| 7.7 | Mincer-Zarnowitz 回归复用 Phase 7A | specification_tests | [ ] | |
| 7.8 | Diebold-Mariano 检验复用 Phase 7A | specification_tests | [ ] | |
| 7.9 | t-GARCH 标准化残差 t 分布检验 | volatility_diagnostics | [ ] | |
| 7.10 | EGARCH 杠杆效应检测 | α<0 负冲击放大波动 (spec 命名) | [ ] | G23 |

---

## 8. ADR 对齐验证

### 8.1 ADR-016 (18 项决策)

| # | 决策点 | 实施一致性 | 状态 | 备注 |
|---|--------|-----------|------|------|
| 8.1.1 | G-ADR1 backcast (EWMA λ=0.94) | [ ] | | |
| 8.1.2 | G-ADR2 SLSQP 约束优化 | [ ] | | |
| 8.1.3 | G-ADR3 sandwich 协方差 | [ ] | | |
| 8.1.4 | G-ADR4 JB Bootstrap 临界值 | [ ] | | |
| 8.1.5 | G-ADR5 t 分布联合 QMLE | [ ] | | |
| 8.1.6 | G-ADR6 多起始点 (4 个) | [ ] | | |
| 8.1.7 | G-ADR7 GARCH-M c·σₜ (推迟) | [ ] | | |
| 8.1.8 | U-ADR1 Schwert ceil | [ ] | | |
| 8.1.9 | U-ADR2 自动方程形式 | [ ] | | |
| 8.1.10 | U-ADR3 MacKinnon 2010 (4 系数) | [ ] | | |
| 8.1.11 | U-ADR4 DF-GLS arch 独立表 | [ ] | | |
| 8.1.12 | U-ADR5 c̄=-7.0(c)/-13.5(ct) | [ ] | | |
| 8.1.13 | U-ADR6 KPSS H0: 平稳性 | [ ] | | |
| 8.1.14 | U-ADR7 KPSS Hobijn 带宽 + Bartlett 核 | [ ] | | |
| 8.1.15 | U-ADR8 PP Schwert 带宽 + Bartlett 核 | [ ] | | |
| 8.1.16 | U-ADR9 Python arch 基准 | [ ] | | |
| 8.1.17 | U-ADR10 方差比全部变体 (含 CLM debiased) | [ ] | | |
| 8.1.18 | U-ADR11 ZA 推迟 v1.7 | [ ] | | |

### 8.2 ADR-017 (命名空间)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 8.2.1 | GARCH 族落 `cpphub::v1::timeseries::garch` | [ ] | |
| 8.2.2 | 单位根落 `cpphub::v1::timeseries::unit_root` | [ ] | |
| 8.2.3 | 头文件 #include 在 namespace 外 | [ ] | project_memory 教训 |

### 8.3 ADR-018 (SLSQP)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 8.3.1 | SLSQP 复用 (不修改 ADR-018 实现) | [ ] | |
| 8.3.2 | SLSQP 12/12 测试无退化 | [ ] | |
| 8.3.3 | GARCH 约束通过 SLSQP 不等式约束实现 | [ ] | |

### 8.4 ADR-013 (Eigen3 隔离)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 8.4.1 | M1 GARCH 不引入 Eigen3 | [ ] | 用 std::vector + SLSQP |
| 8.4.2 | M2 单位根不引入 Eigen3 | [ ] | 用 detail::ols_simple (Gauss-Jordan) |
| 8.4.3 | M2 方差比不引入 Eigen3 | [ ] | 纯序列运算 |
| 8.4.4 | cpphub_timeseries 为 header-only INTERFACE 库 | [ ] | 不链接 Eigen3 |

---

## 9. Scope 边界验证

### 9.1 已实现 (在 scope 内)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 9.1.1 | GARCH(1,1) 已实现 | [ ] | |
| 9.1.2 | EGARCH 已实现 | [ ] | |
| 9.1.3 | GJR-GARCH 已实现 | [ ] | |
| 9.1.4 | ADF 已实现 | [ ] | |
| 9.1.5 | PP 已实现 | [ ] | |
| 9.1.6 | KPSS 已实现 | [ ] | |
| 9.1.7 | DF-GLS 已实现 | [ ] | |
| 9.1.8 | Lo-MacKinlay Z₁/Z₂ 已实现 | [ ] | |
| 9.1.9 | Chow-Denning 已实现 | [ ] | |
| 9.1.10 | CLM debiased 已实现 | [ ] | |

### 9.2 未实现 (在 scope 外, 确认未越界)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 9.2.1 | APARCH 未实现 | [ ] | 推迟 v1.6+ |
| 9.2.2 | FIGARCH 未实现 | [ ] | 推迟 v1.6+ |
| 9.2.3 | IGARCH 未实现 | [ ] | 推迟 v1.6+ |
| 9.2.4 | GARCH-M 未实现 | [ ] | 推迟 v1.6+ |
| 9.2.5 | Ng-Perron M 检验族未实现 | [ ] | 推迟 v1.7 |
| 9.2.6 | Zivot-Andrews 未实现 | [ ] | 推迟 v1.7 |
| 9.2.7 | ARIMA 未实现 | [ ] | M3 |
| 9.2.8 | MIDAS 未实现 | [ ] | M4 |
| 9.2.9 | VAR/DCC 未实现 | [ ] | v1.7 |

---

## 10. 代码规范验证

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 10.1 | 头文件 `#include` 在 namespace 外 | [ ] | project_memory 教训 (C2065) |
| 10.2 | 无参数名遮蔽函数名 (如 `make_psd` vs `bool make_psd`) | [ ] | v1.4.2 教训 |
| 10.3 | `extern` 声明不在 namespace 内 | [ ] | timegm 教训 |
| 10.4 | 接口签名与 spec 一致 (参数顺序/类型/命名) | [ ] | |
| 10.5 | Result 结构体含所有 spec 定义字段 | [ ] | |
| 10.6 | 异常处理 (std::invalid_argument / std::runtime_error) | [ ] | |
| 10.7 | 无 magic number (常数用 constexpr) | [ ] | |
| 10.8 | 排幻觉点注释标注 (如 `// 排幻觉点 G13: ...`) | [ ] | |

---

## 11. 复用验证 (不重复实现)

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 11.1 | GARCH 诊断复用 Phase 7A `volatility_diagnostics` | [ ] | JB/LB 直接调用 |
| 11.2 | JB Bootstrap 复用 v1.5 `block_bootstrap` | [ ] | |
| 11.3 | HAC 协方差复用 v1.5 `hac_vcov` | [ ] | PP/KPSS 长期方差 |
| 11.4 | Bartlett 核权重复用 `hac_kernels` | [ ] | |
| 11.5 | OLS 辅助回归复用 `detail::ols_simple` | [ ] | ADF/DF-GLS |
| 11.6 | 多重检验修正复用 Phase 7A `multiple_test_correction` | [ ] | U18 |
| 11.7 | SLSQP 复用 ADR-018 `optimizer.hpp` | [ ] | 不修改 |
| 11.8 | Mincer-Zarnowitz 复用 Phase 7A `specification_tests` | [ ] | 集成测试 |

---

## 12. 文档对齐验证

| # | 检查项 | 状态 | 备注 |
|---|--------|------|------|
| 12.1 | ADR-016 18 项决策全部实施 | [ ] | |
| 12.2 | ADR-017 命名空间遵守 | [ ] | |
| 12.3 | ADR-018 SLSQP 复用 | [ ] | |
| 12.4 | DEVELOPMENT_LOG.md 更新 Phase 7B 进度 | [ ] | |
| 12.5 | spec 中所有幻觉点注释保留到代码 | [ ] | |
| 12.6 | verify 脚本可追溯 (打印 arch 源码 + 手算) | [ ] | |
| 12.7 | R 基准脚本含 `.libPaths()` 显式加载 | [ ] | Rscript 教训 |

---

## 13. 性能验证

| # | 检查项 | 基准 | 状态 | 备注 |
|---|--------|------|------|------|
| 13.1 | GARCH(1,1) T=5000 估计 | < 5 sec | [ ] | |
| 13.2 | EGARCH T=5000 估计 | < 8 sec | [ ] | |
| 13.3 | GJR-GARCH T=5000 估计 | < 8 sec | [ ] | |
| 13.4 | ADF T=1000 | < 1 sec | [ ] | |
| 13.5 | 方差比 T=1000, k_list={2,5,10,20} | < 1 sec | [ ] | |
| 13.6 | 全量 ctest (~2205 用例) | < 30 min | [ ] | 三平台 |

---

## 14. 风险项验收

| # | 风险 | 缓解措施落实 | 状态 | 备注 |
|----|------|-------------|------|------|
| 14.1 | MacKinnon 2010 系数表录入错误 | arch 源码逐系数核对 | [ ] | |
| 14.2 | DF-GLS c/ct 临界值混用 | 独立系数表, c/ct 分别存储 | [ ] | U7 |
| 14.3 | GARCH 局部最优 | 多起始点 (G-ADR6, 4 起始) | [ ] | |
| 14.4 | backcast 实现错误 | arch `volatility.py:1161-1168` 源码核对 | [ ] | |
| 14.5 | EGARCH 符号约定不一致 | arch vs rugarch 源码核对 | [ ] | G23 |
| 14.6 | Chen-Deo 归属错误 | 已修正为 CLM 1997 | [ ] | U17 |
| 14.7 | 跨库 solver 差异 | arch SLSQP 主基准, rugarch 放宽 1e-8 | [ ] | G22 |
| 14.8 | PP σ²_ε df-corrected | 已修正为 s²=SSR/(n-k) | [ ] | U5-sigma2 |
| 14.9 | GJR γ 约束过严 | 已改为 [-1,1] 允许负值 | [ ] | G-gamma-sign |
| 14.10 | ν 边界不一致 | 已统一 [2.05, 500.0] | [ ] | G-nu-bound |

---

## 15. 最终签字

### 15.1 验收统计

| 维度 | 总项数 | 通过 | 未通过 | 通过率 |
|------|--------|------|--------|--------|
| 1. 交付物完整性 | 43 | | | |
| 2. 编译与跨平台 | 15 | | | |
| 3. M1 GARCH 数值基准 | 48 | | | |
| 4. M2 单位根数值基准 | 38 | | | |
| 5. M2 方差比数值基准 | 19 | | | |
| 6. 幻觉点核查 | 53 | | | |
| 7. 端到端集成 | 10 | | | |
| 8. ADR 对齐 | 28 | | | |
| 9. Scope 边界 | 19 | | | |
| 10. 代码规范 | 8 | | | |
| 11. 复用验证 | 8 | | | |
| 12. 文档对齐 | 7 | | | |
| 13. 性能验证 | 6 | | | |
| 14. 风险项 | 10 | | | |
| **总计** | **312** | | | |

> 注: §1 项数 33→43 (§1.3 脚本 8→17 项, §1.4 套件 12→13 套), 总计 302→312 (2026-08-15 修正)。

### 15.2 验收结论

- [ ] **全部通过** — Phase 7B 验收完成, 可进入 v1.6 发布
- [ ] **有条件通过** — 存在未通过项, 附 issue 列表, 限期修复
- [ ] **未通过** — 存在极高严重性未通过项, 需返工

### 15.3 签字

| 角色 | 姓名 | 日期 | 签字 |
|------|------|------|------|
| 实施者 | | | |
| 审计者 | | | |
| 架构组 | | | |

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

## 附录 B: verify 脚本核查模式 (参照 v1.4.0/1.4.1)

每个 verify 脚本必须包含:
1. **arch 源码打印**: `print(getsource(arch.unitroot.XXX))` 或 `inspect.getsource()`
2. **小样本手算**: 构造 T=5~10 的简单序列, 手动计算预期值
3. **arch baseline**: 调用 arch API 生成基准数值
4. **C++ 输出对照**: 打印 C++ 实现结果, 与 arch baseline 逐数值对比
5. **容差判定**: 相对误差 ≤ spec 标注容差
