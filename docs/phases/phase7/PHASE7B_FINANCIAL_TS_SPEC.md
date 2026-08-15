# Phase 7B 执行规格书 - 金融时间序列模块 (v1.6 M1/M2)

> **版本归属**: **v1.6** (Phase 7B, 与 Phase 7A 证伪统计量并行)
> **目标**: 实现金融时间序列核心模块 — GARCH 族 (M1) + 单位根与方差比检验 (M2), 为 Research OS 因子失效诊断提供波动率建模与伪回归检测基础设施
> **覆盖范围**:
>   - M1: GARCH(1,1) / EGARCH / GJR-GARCH QMLE + sandwich 协方差 + 标准化残差诊断
>   - M2: ADF / PP / KPSS / DF-GLS 单位根检验 + Lo-MacKinlay / Chow-Denning / CLM debiased 方差比检验
> **前置**:
>   - Phase 6 (v1.5) 三平台 1767/1767 测试通过
>   - Phase 7A (v1.6) 证伪统计量补齐 1962/1962 测试通过 (volatility_diagnostics/multiple_test_correction 等 13 头文件可复用)
>   - **SLSQP 优化器扩展已完成** (ADR-018, commit `e2f3d5c`, 12/12 测试三平台通过)
> **里程碑**: M1 GARCH 族 → M2 单位根与方差比 (M1/M2 可并行, 无相互依赖)
>
> **Scope 声明**:
> - **严格聚焦金融时间序列** (波动率建模 + 单位根检验 + 方差比检验), 不含 ARIMA (M3)、MIDAS (M4)、VAR/DCC (v1.7)
> - **GARCH 族仅实现 3 个核心模型** (GARCH(1,1)/EGARCH/GJR-GARCH), APARCH/FIGARCH/IGARCH/GARCH-M 推迟到 v1.6+ 或 v1.7
> - **单位根仅实现 4 个核心检验** (ADF/PP/KPSS/DF-GLS), Ng-Perron M 检验族/Zivot-Andrews 推迟到 v1.7
> - **方差比实现 4 变体** (Lo-MacKinlay Z₁/Z₂ + Chow-Denning + CLM debiased)
> - **复用 v1.5/v1.6 基础设施**: MLEEstimator 框架、HAC (Newey-West)、Bootstrap (4 类)、Phase 7A volatility_diagnostics、SLSQP 优化器
> - **教材锚点**: Tsay 3ed (2010) Ch 3-7 / Hamilton 1994 Ch 17 / Bollerslev 1986 / Nelson 1991 / GJR 1993 / KPSS 1992 / ERS 1996 / Lo-MacKinlay 1988 / MacKinnon 2010
> - **对照库**: Python `arch` (主基准, U-ADR9) + R `rugarch` (GARCH 交叉验证), 容差 1e-10 至 1e-12
>
> **关联**:
> - [ADR-016: 金融时间序列实施边界 (18 项)](../../decisions/ADR-016_FINANCIAL_TIMESERIES_BOUNDARY.md) (Accepted 2026-08-15)
> - [ADR-017: 时序模块命名空间 `cpphub::v1::timeseries`](../../decisions/ADR-017_TIMESERIES_NAMESPACE.md) (Accepted 2026-08-15)
> - [ADR-018: SLSQP 优化器实现边界](../../decisions/ADR-018_SLSQP_BOUNDARY.md) (Accepted 2026-08-15, 已实现)
> - [FINANCIAL_TIMESERIES_RESEARCH.md](../../research/FINANCIAL_TIMESERIES_RESEARCH.md) v3.2 (调研报告, 4 子 agent review, 60 幻觉点验证)
> - [PHASE7A_FALSIFICATION_SPEC.md](./PHASE7A_FALSIFICATION_SPEC.md) (证伪统计量 spec, M1 复用 volatility_diagnostics)
> - [PHASE6_ECONOMETRICS_SPEC.md](../phase6/PHASE6_ECONOMETRICS_SPEC.md) (v1.5 基础设施: MLE/OLS/HAC/Bootstrap)

---

## 1. 交付物清单

### 1.1 新增文件结构 (按 ADR-017 命名空间组织)

```
# M1: GARCH 族 (cpphub::v1::timeseries::garch)
include/cpphub/timeseries/garch/
├── garch_model.hpp           # 新增 P0: GARCH(1,1) QMLE + 递归方差 + backcast
├── egarch_model.hpp          # 新增 P0: EGARCH QMLE + 对数方差 + 杠杆效应
├── gjr_garch_model.hpp       # 新增 P0: GJR-GARCH QMLE + 非对称项
├── garch_distribution.hpp    # 新增 P0: 分布族 (Normal / t / GED) + 似然函数
├── garch_forecast.hpp        # 新增 P0: 多步方差预测 (收敛到无条件方差)
└── garch_diagnostics.hpp     # 新增 P0: 标准化残差诊断 (调用 Phase 7A volatility_diagnostics)

# M2: 单位根与方差比 (cpphub::v1::timeseries::unit_root)
include/cpphub/timeseries/unit_root/
├── adf_test.hpp              # 新增 P0: ADF 检验 (Schwert lag + 自动方程形式 + MacKinnon 2010)
├── pp_test.hpp               # 新增 P0: PP 检验 (Bartlett 核 + Schwert 带宽)
├── kpss_test.hpp             # 新增 P0: KPSS 检验 (H0: 平稳性 + Hobijn 带宽 + Bartlett 核)
├── df_gls_test.hpp           # 新增 P0: DF-GLS 检验 (GLS detrending + AIC lag + arch 独立临界值)
├── variance_ratio_test.hpp   # 新增 P0: Lo-MacKinlay Z₁/Z₂ + Chow-Denning + CLM debiased
├── mackinnon_cv.hpp          # 新增 P0: MacKinnon 2010 response surface 临界值 (ADF/PP 共享)
└── unit_root_common.hpp      # 新增 P0: 共享工具 (Schwert lag / 自动方程形式 / 长期方差)

# 测试基准脚本 (不入版本控制, .gitignore 排除)
tests/fixtures/timeseries/
├── verify_garch.py           # Python arch 对照基准生成
├── verify_egarch.py          # Python arch EGARCH 对照
├── verify_gjr.py             # Python arch GJR 对照
├── verify_adf.py             # Python arch ADF 对照
├── verify_pp.py              # Python arch PP 对照
├── verify_kpss.py            # Python arch KPSS 对照
├── verify_df_gls.py          # Python arch DF-GLS 对照
└── verify_variance_ratio.py  # Python arch 方差比对照
```

### 1.2 新增测试套件

| 测试套件 | 用例数 | 覆盖模块 | 优先级 |
|----------|--------|----------|--------|
| `test_garch_model` | 20 | GARCH(1,1) QMLE + backcast + 预测 + sandwich | P0 |
| `test_egarch_model` | 18 | EGARCH QMLE + 杠杆效应 + E\|z\| 常数 | P0 |
| `test_gjr_garch_model` | 18 | GJR-GARCH QMLE + 非对称项 + 平稳性条件 | P0 |
| `test_garch_distribution` | 12 | Normal/t/GED 似然 + t 分布联合估计 | P0 |
| `test_garch_diagnostics` | 15 | 标准化残差 JB/LB/z² LB (复用 Phase 7A) | P0 |
| `test_adf_test` | 20 | Schwert lag + nc/c/ct + MacKinnon 2010 | P0 |
| `test_pp_test` | 15 | Bartlett 核 + Schwert 带宽 | P0 |
| `test_kpss_test` | 15 | H0: 平稳性 + Hobijn 带宽 + Bartlett 核 | P0 |
| `test_df_gls_test` | 18 | GLS detrending + AIC + arch 独立临界值 | P0 |
| `test_variance_ratio_test` | 20 | Z₁/Z₂ + Chow-Denning + CLM debiased | P0 |
| `test_mackinnon_cv` | 10 | MacKinnon 2010 response surface 系数 | P0 |
| `test_integration_phase7b` | 10 | 端到端: GARCH→VaR + ADF→伪回归诊断 | - |

**新增测试总数**: ~201

### 1.3 必须达到的数值基准

| 基准 | 容差 | 验收方式 |
|------|------|----------|
| GARCH(1,1) 参数 vs Python `arch.arch_model(o='Constant', vol='GARCH')` | 1e-10 | `test_garch_model.cpp` |
| GARCH(1,1) 参数 vs R `rugarch::ugarchfit(model="sGARCH")` | 1e-8 | `verify_garch.py` 生成基准 |
| GARCH sandwich 协方差 vs Python `arch` (cov_type='robust') | 1e-10 | `test_garch_model.cpp` |
| GARCH backcast 方差 vs Python `arch` (backcast=EWMA λ=0.94) | 1e-12 | `test_garch_model.cpp` |
| EGARCH 参数 vs Python `arch.arch_model(vol='EGARCH')` | 1e-10 | `test_egarch_model.cpp` |
| EGARCH E\|z\| = √(2/π) vs arch 源码 `np.sqrt(2/np.pi)` | 1e-15 | `test_egarch_model.cpp` |
| GJR-GARCH 参数 vs Python `arch.arch_model(vol='GARCH', o=1)` | 1e-10 | `test_gjr_garch_model.cpp` |
| GJR 平稳性 α+γ/2+β<1 vs arch 约束矩阵 | 1e-12 | `test_gjr_garch_model.cpp` |
| t-GARCH 自由度 ν vs Python `arch.arch_model(dist='t')` | 1e-8 | `test_garch_distribution.cpp` |
| ADF 统计量 vs Python `arch.unitroot.ADF` | 1e-10 | `test_adf_test.cpp` |
| ADF 临界值 vs MacKinnon 2010 response surface | 1e-12 | `test_mackinnon_cv.cpp` |
| PP 统计量 vs Python `arch.unitroot.PhillipsPerron` | 1e-10 | `test_pp_test.cpp` |
| KPSS 统计量 vs Python `arch.unitroot.KPSS` | 1e-10 | `test_kpss_test.cpp` |
| DF-GLS 统计量 vs Python `arch.unitroot.DFGLS` | 1e-10 | `test_df_gls_test.cpp` |
| DF-GLS c̄=-7.0 (c) / c̄=-13.5 (ct) vs ERS 1996 | 精确 | `test_df_gls_test.cpp` |
| Lo-MacKinlay Z₁/Z₂ vs Python `arch.unitroot.VarianceRatio` | 1e-10 | `test_variance_ratio_test.cpp` |
| CLM debiased VR vs Python `arch` (debiased=True) | 1e-10 | `test_variance_ratio_test.cpp` |
| Chow-Denning 联合检验 vs R `vrtest::Chow.Denning` | 1e-8 | `test_variance_ratio_test.cpp` |

---

## 2. M1: GARCH 族 (Week 1-3)

> **定位**: 金融时间序列核心模块 — 波动率建模
> **教材锚点**: Tsay 3ed Ch 5-6 / Bollerslev 1986 / Nelson 1991 / GJR 1993 / Bollerslev-Wooldridge 1992
> **对照库**: Python `arch.arch_model` (主基准) + R `rugarch::ugarchfit` (交叉验证)
> **命名空间**: `cpphub::v1::timeseries::garch` (ADR-017)

### 2.0 公共基础设施

#### 2.0.1 `garch_distribution.hpp` - 分布族与似然函数

**教材锚点**: Bollerslev 1986 (Normal) / Bollerslev 1987 RES (t) / Nelson 1991 (GED)

**接口签名**:

```cpp
#pragma once
#include <vector>
#include <functional>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace garch {

// GARCH 分布族 (G-ADR5: t 分布与 (ω,α,β) 联合 QMLE)
enum class GarchDist {
    Normal,   ///< 标准正态: ℓ = -0.5·[log(2π) + log(hₜ) + ε²ₜ/hₜ]
    StudentT, ///< Student-t (ν 自由度): 含 Γ((ν+1)/2)/Γ(ν/2) 项
    GED       ///< 广义误差分布 (Nelson 1991): 含 Γ(1/ν)/Γ(3/ν) 项
};

// 单观测对数似然 (G3: 必须含 -0.5·log(2π) 常数项)
// @param residual 残差 εₜ
// @param variance 条件方差 hₜ
// @param dist 分布族
// @param nu t/GED 自由度/形状参数 (Normal 时忽略)
// @return 单观测对数似然 ℓₜ
// @throws std::invalid_argument 若 variance <= 0 或 nu <= 2.05 (t 分布需 ν>2 保证方差存在, arch 下界 2.05)
Real log_likelihood_term(Real residual, Real variance,
                          GarchDist dist, Real nu = 0.0);

// 批量对数似然 (求和)
Real log_likelihood(const std::vector<Real>& residuals,
                    const std::vector<Real>& variances,
                    GarchDist dist, Real nu = 0.0);

// 解析梯度 (可选, 未提供时 SLSQP 用中心差分)
// @return {dℓ/dω, dℓ/dα, dℓ/dβ, dℓ/dν} (Normal 时第 4 项为 0)
std::vector<Real> log_likelihood_gradient(const std::vector<Real>& residuals,
                                           const std::vector<Real>& variances,
                                           GarchDist dist, Real nu = 0.0);

}  // namespace garch
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
```

**幻觉点映射**:
- **G3**: 对数似然必须含 `-0.5·log(2π)` 常数项 (arch `arch/univariate/distribution.py` Normal 类 `lls = -0.5*(log(2*pi)+log(sigma2)+resids**2/sigma2)` 已验证)
- **G14**: t 分布自由度 ν 与 (ω,α,β) **联合 QMLE** (非两步矩估计), Bollerslev 1987 RES 69(3):542
- **G17**: AIC/BIC 跨模型比较必须用完整似然 (含常数项), 否则 t-GARCH vs Normal-GARCH 比较失真

**算法步骤**:

1. Normal: `ℓₜ = -0.5·[log(2π) + log(hₜ) + ε²ₜ/hₜ]`
2. StudentT (Bollerslev 1987):
   ```
   ℓₜ = log(Γ((ν+1)/2)) - log(Γ(ν/2)) - 0.5·log(ν-2) - 0.5·log(π)
        - 0.5·log(hₜ) - ((ν+1)/2)·log(1 + ε²ₜ/((ν-2)·hₜ))
   ```
   注: ν>2 保证方差存在; ν→∞ 时退化为 Normal
3. GED (Nelson 1991, arch `arch/univariate/distribution.py:1137-1142` 实测):
   ```
   c = 2^(-1/ν) · sqrt(Γ(1/ν) / Γ(3/ν))  (缩放常数, 含 2^(-1/ν) 因子)
   ℓₜ = log(ν) - log(c) - log(Γ(1/ν)) - (1+1/ν)·log(2) - 0.5·log(hₜ)
        - 0.5·|εₜ/(sqrt(hₜ)·c)|^ν
   ```
   **排幻觉点 G-GED1** (arch 源码实测): 指数为 **ν** (非 2ν); 缩放常数 c 含 **2^(-1/ν)** 因子 (非裸 λ);
   常数项为 `log(ν) - log(c) - log(Γ(1/ν)) - (1+1/ν)·log(2)` (非 `log(ν/(λ·2^(1+1/ν)·Γ(1/ν)))`)

**测试矩阵**:

| 用例 | 验证点 | 容差 |
|------|--------|------|
| Normal 单观测似然 vs 手算 (ε=0.1, h=0.04) | G3 常数项 | 1e-15 |
| t 分布 ν=5 vs arch `dist='t'` | G14 联合估计 | 1e-10 |
| t 分布 ν→1000 退化到 Normal | 渐近行为 | 1e-6 |
| GED ν=2 (正态特殊情况) | GED 退化 | 1e-10 |
| 批量似然求和 vs 逐项相加 | 一致性 | 1e-15 |

---

#### 2.0.2 `garch_model.hpp` - GARCH(1,1) QMLE

**教材锚点**: Bollerslev 1986 / Tsay 3ed Ch 5 / arch 源码 `arch/univariate/mean.py` (`arch_model()` 函数)

**接口签名**:

```cpp
#pragma once
#include <vector>
#include "cpphub/core/types.hpp"
#include "cpphub/calibration/optimizer.hpp"  // SLSQP (ADR-018)
#include "cpphub/timeseries/garch/garch_distribution.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace garch {

// GARCH(1,1) 参数
struct GarchParams {
    Real omega;    ///< 常数项 (ω > 0, G4)
    Real alpha;    ///< ARCH 系数 (α ≥ 0, G4)
    Real beta;     ///< GARCH 系数 (β ≥ 0, G4)
    Real nu;       ///< t/GED 自由度 (Normal 时忽略, G-ADR5)
};

// GARCH(1,1) 估计配置
struct GarchConfig {
    GarchDist dist = GarchDist::Normal;           ///< 分布族 (G-ADR5)
    Real backcast_lambda = 0.94;                  ///< EWMA backcast λ (G-ADR1, G1)
    Size backcast_window = 0;                     ///< backcast 窗口 τ (0 => min(75, T), arch 默认)
    bool use_multistart = true;                   ///< 多起始点 (G-ADR6, G16)
    Size n_multistart = 4;                        ///< 多起始点数 (HR/ML/user/random)
    bool compute_sandwich = true;                 ///< QMLE sandwich 协方差 (G-ADR3, G9)
    bool compute_diagnostics = true;              ///< 标准化残差诊断 (G-ADR4, G11)
    Size bootstrap_reps = 1000;                   ///< JB 检验 Bootstrap 次数 (G-ADR4)
    SLSQP::Config optimizer_config;               ///< SLSQP 配置 (ADR-018)
};

// GARCH(1,1) 估计结果
struct GarchResult {
    GarchParams params;                           ///< 参数估计
    std::vector<Real> conditional_variances;      ///< 条件方差 hₜ 序列
    std::vector<Real> residuals;                  ///< 残差 εₜ 序列
    std::vector<Real> std_residuals;              ///< 标准化残差 zₜ = εₜ/√hₜ (G11)
    Real log_likelihood;                          ///< 对数似然值 (G3: 含常数项)
    Real aic;                                     ///< AIC (G17: 用完整似然)
    Real bic;                                     ///< BIC (G17: 用完整似然)
    std::vector<std::vector<Real>> vcov;          ///< 参数协方差矩阵 (sandwich, G9)
    std::vector<Real> std_errors;                 ///< 标准误
    bool converged;                               ///< 收敛标志
    Size n_iterations;                            ///< 迭代次数
    std::string message;                          ///< 诊断消息
};

// GARCH(1,1) QMLE 估计
// @param data 收益率序列 (rₜ, 长度 T)
// @param config 估计配置
// @return 估计结果 (含条件方差/残差/标准化残差/协方差)
// @throws std::invalid_argument 若 T < 10 或 data 含 NaN
GarchResult estimate_garch11(const std::vector<Real>& data,
                              const GarchConfig& config = GarchConfig{});

// GARCH(1,1) 多步方差预测 (G13)
// @param params 估计参数
// @param last_variance 最后一期条件方差 h_T
// @param last_residual 最后一期残差 ε_T
// @param horizon 预测步数 k
// @return k 步方差预测 {h_{T+1}, ..., h_{T+k}}
// 注: h_{T+1} = ω + α·ε²_T + β·h_T  (已知 ε²_T)
//     h_{T+k} = ω·(1-φ^{k-1})/(1-φ) + φ^{k-1}·h_{T+1}  (k≥1, φ=α+β)
//     等价: h_{T+k} = σ̄² + φ^{k-1}·(h_{T+1} - σ̄²), σ̄²=ω/(1-α-β)
//     k→∞: h_{T+k} → σ̄² (无条件方差)
// 排幻觉点 G13: 指数为 k-1 (非 k), 起点为 h_{T+1} (非 h_T)
//   错误公式 h_{T+k} = ω·(1-φ^k)/(1-φ) + φ^k·h_T 混淆了起点与递归次数
std::vector<Real> forecast_garch11(const GarchParams& params,
                                    Real last_variance, Real last_residual,
                                    Size horizon);

}  // namespace garch
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
```

**幻觉点映射**:
- **G1**: 方差初始值用 EWMA backcast, λ=0.94, τ=min(75,T) (arch 默认)
  - **排幻觉点 G1-dir** (arch `arch/univariate/volatility.py:1161-1168` 实测): arch 用**前 τ 个**残差 `resids[:tau]`, 非**末尾** τ 个
  - **排幻觉点 G1-norm** (arch 同上): 归一化为 `wᵢ = λⁱ/Σⱼλʲ` (有限样本归一化, sum(w)=1), 非 `(1-λ)·Σλⁱ` (无穷级数近似)
  - 正确公式: `σ²₀ = Σᵢ₌₀^{τ-1} (λⁱ/Σⱼ₌₀^{τ-1}λʲ)·ε²ᵢ` (用前 τ 个残差, 归一化权重)
- **G2**: 残差初始值 ε₀ 用样本均值 μ̂ (非 0)
- **G3**: 对数似然含 `-0.5·log(2π)` 常数项
- **G4**: 参数约束 `ω>0, α≥0, β≥0, α+β<1` (SLSQP 不等式约束, ADR-018)
- **G9**: QMLE 协方差用 sandwich `V = H⁻¹·S·H⁻¹` (非 Hessian 逆, G-ADR3)
- **G11**: 标准化残差 `zₜ = εₜ/√hₜ` (非 εₜ/σₜ)
- **G13**: 多步预测收敛到无条件方差 `σ̄²=ω/(1-α-β)`
- **G16**: 多起始点避免局部最优 (G-ADR6)
- **G17**: AIC/BIC 用完整似然 (G3 常数项)

**算法步骤** (GARCH(1,1) QMLE):

```
输入: 收益率序列 r = {r₁, ..., r_T}, 配置 config

Step 1: 预处理
  1.1 计算样本均值 μ̂ = mean(r)
  1.2 计算残差 εₜ = rₜ - μ̂ (G2)
  1.3 计算 backcast 方差 (G1, arch `arch/univariate/volatility.py:1161-1168`):
      τ = config.backcast_window == 0 ? min(75, T) : config.backcast_window
      wᵢ = λⁱ / Σⱼ₌₀^{τ-1} λʲ   (有限样本归一化, sum(w)=1)
      σ²₀ = Σᵢ₌₀^{τ-1} wᵢ · ε²ᵢ   (用前 τ 个残差, 非末尾)

Step 2: 定义目标函数 (负对数似然)
  2.1 给定参数 θ = (ω, α, β, [ν]), 递归计算条件方差:
      h₁ = σ²₀ (backcast)
      hₜ = ω + α·ε²ₜ₋₁ + β·hₜ₋₁  (t = 2, ..., T)
  2.2 计算对数似然:
      ℓ(θ) = Σₜ log_likelihood_term(εₜ, hₜ, dist, ν)  (G3)
  2.3 返回 -ℓ(θ) (SLSQP 最小化)

Step 3: 定义约束 (G4, ADR-018 SLSQP)
  3.1 不等式约束 (c_i(x) >= 0, ADR-018 约定):
      c₁ = ω       (ω > 0, 实际 ω >= 1e-8 避免 0)
      c₂ = α       (α >= 0)
      c₃ = β       (β >= 0)
      c₄ = 1 - α - β - 1e-8  (α + β < 1, 平稳性)
  3.2 边界约束 (Bounds, arch `arch/univariate/distribution.py:505,508` 实测 StudentsT 类边界):
      ω ∈ [1e-8, 100], α ∈ [0, 1], β ∈ [0, 1], ν ∈ [2.05, 500.0]  (arch 默认, 非 [2.1, 100])

Step 4: 多起始点优化 (G-ADR6, G16)
  注 (arch `arch/univariate/volatility.py:1224-1244` 实测): arch 用 64 候选 (alphas×gammas×abg = 4×4×4) 计算
  似然, 选最优起始点后做 1 次 SLSQP; 本 spec 用 4 起始点各做 SLSQP 取最优, 策略不同
  但结果应一致 (全局最优唯一时). 实施时可选 arch 策略 (64 候选选优) 或 spec 策略.
  4.1 起始点 1 (HR 高频): ω=0.1·var, α=0.1, β=0.85
  4.2 起始点 2 (ML 极大似然): ω=0.05·var, α=0.05, β=0.9
  4.3 起始点 3 (用户指定, 可选): config 提供或跳过
  4.4 起始点 4 (随机扰动): 在起始点 1 基础上加 10% 高斯扰动
  4.5 对每个起始点调用 SLSQP::minimize, 取对数似然最大者

Step 5: 计算 QMLE 协方差 (G-ADR3, G9)
  注 (arch `arch/univariate/base.py:977-988` 实测): arch 的 sandwich 实现为
    H /= nobs; inv_H = inv(H); score_cov = np.cov(scores.T)  # 中心化, 除 (nobs-1)
    V = inv_H · score_cov · inv_H / nobs
  本 spec 用非中心化 OPG 求和 S = Σ gₜ gₜ', 两者尺度差 nobs/(nobs-1) 因子 (大样本可忽略).
  verify 脚本需核对 arch_result.param_cov vs spec 的 V, 必要时调整 S 定义.
  5.1 计算 Hessian H (数值差分, ADR-018 numerical_hessian)
  5.2 计算外积得分 S = Σₜ (∂℗ₜ/∂θ)·(∂℗ₜ/∂θ)'  (OPG, 非中心化)
  5.3 sandwich: V = H⁻¹·S·H⁻¹  (G9)
  5.4 标准误 SE = sqrt(diag(V))

Step 6: 标准化残差诊断 (G-ADR4, G11)
  6.1 zₜ = εₜ / √hₜ
  6.2 调用 Phase 7A volatility_diagnostics:
      - JB 正态性 (Bootstrap 临界值, G-ADR4)
      - LB 自相关 (zₜ)
      - LB ARCH 效应 (z²ₜ, G12: 滞后用 floor(log(T)) 或 Schwert)

Step 7: 返回结果
  7.1 填充 GarchResult (参数/方差/残差/标准化残差/似然/AIC/BIC/协方差/SE)
```

**测试矩阵**:

| 用例 | 验证点 | 容差 |
|------|--------|------|
| GARCH(1,1) Normal 参数 vs arch | ω/α/β 数值一致 | 1e-10 |
| GARCH(1,1) backcast vs arch (λ=0.94) | σ²₀ 一致 | 1e-12 |
| GARCH(1,1) 条件方差序列 vs arch | hₜ 序列一致 | 1e-10 |
| GARCH(1,1) sandwich 协方差 vs arch (robust) | V 矩阵一致 | 1e-10 |
| GARCH(1,1) 标准误 vs arch | SE 一致 | 1e-10 |
| GARCH(1,1) 多步预测 vs arch.forecast(horizon=10) | h_{T+k} 一致 | 1e-10 |
| GARCH(1,1) 多步预测收敛到 σ̄²=ω/(1-α-β) | k→∞ 行为 | 1e-6 |
| GARCH(1,1) 多起始点避免局部最优 | 4 起始点取最优 | - |
| GARCH(1,1) 约束违反检测 (ω<0 时拒绝) | G4 约束生效 | - |
| t-GARCH(1,1) ν 联合估计 vs arch (dist='t') | ν 数值一致 | 1e-8 |
| 标准化残差 JB 检验 (Bootstrap) | G-ADR4 | 1e-6 |
| 标准化残差 LB 检验 vs R `Box.test(z, lag=10)` | G11 | 1e-8 |
| z²ₜ LB 检验 (ARCH 效应) vs R `Box.test(z², lag=10)` | G12 | 1e-8 |
| AIC/BIC vs arch (完整似然) | G17 | 1e-10 |
| 收敛失败检测 (T<10 或全 0 数据) | 异常处理 | - |
| 与 rugarch 交叉验证 (sGARCH) | R 对照 | 1e-8 |
| 大样本 T=5000 性能 | < 5 sec | - |
| 含 NaN 数据拒绝 | 异常处理 | - |

---

#### 2.0.3 `egarch_model.hpp` - EGARCH QMLE

**教材锚点**: Nelson 1991 / Tsay 3ed Ch 6 / arch 源码 `arch/univariate/mean.py` (`arch_model()` vol='EGARCH')

**接口签名**:

```cpp
// EGARCH 参数 (Nelson 1991)
// 排幻觉点 G23 (arch `arch/univariate/recursions_python.py:319-348` 实测):
//   arch 参数顺序: alpha[p] = symmetric |z|-E|z| 项, gamma[o] = asymmetric z 项
//   本 spec 参数顺序与 Nelson 1991 一致, 与 arch 相反:
//     spec alpha = Nelson θ = asymmetric z 项 = arch gamma[o]
//     spec gamma = Nelson γ = symmetric |z|-E|z| 项 = arch alpha[p]
//   verify_egarch.py 需做参数映射: arch_alpha = spec_gamma, arch_gamma = spec_alpha
struct EGarchParams {
    Real omega;   ///< 常数项 ω
    Real alpha;   ///< asymmetric 系数 (θ in Nelson 1991, 乘以 zₜ₋₁; 对应 arch gamma[o])
    Real beta;    ///< GARCH 系数 β (持续性)
    Real gamma;   ///< symmetric size 系数 (γ in Nelson 1991, 乘以 (|z|-E|z|); 对应 arch alpha[p])
    Real nu;      ///< t/GED 自由度 (Normal 时忽略)
};

// EGARCH 估计结果 (结构同 GarchResult, 参数替换为 EGarchParams)
struct EGarchResult {
    EGarchParams params;
    std::vector<Real> conditional_variances;
    std::vector<Real> residuals;
    std::vector<Real> std_residuals;
    Real log_likelihood;
    Real aic;
    Real bic;
    std::vector<std::vector<Real>> vcov;
    std::vector<Real> std_errors;
    bool converged;
    Size n_iterations;
    std::string message;
};

// EGARCH QMLE 估计
EGarchResult estimate_egarch(const std::vector<Real>& data,
                              const GarchConfig& config = GarchConfig{});
```

**幻觉点映射**:
- **G5**: `E|z| = √(2/π) ≈ 0.7979` (非 `2/π`), arch 源码 `arch/univariate/volatility.py` EGARCH 类 `norm_const = np.sqrt(2/np.pi)` 已验证
- **G6**: 非对称项用标准化残差 `zₜ₋₁ = εₜ₋₁/√hₜ₋₁` (非未标准化 εₜ₋₁)
- **G23**: EGARCH θ 符号约定 (α 在本 spec 中对应 Nelson 1991 的 θ), 实施时核查 arch vs rugarch 一致性
- **G3/G9/G11/G16/G17**: 同 GARCH(1,1)

**算法步骤** (EGARCH QMLE):

```
输入: 收益率序列 r = {r₁, ..., r_T}

Step 1: 预处理 (同 GARCH(1,1) Step 1)

Step 2: 定义目标函数
  2.1 给定参数 θ = (ω, α, β, γ, [ν]), 递归计算对数方差:
      E|z| = sqrt(2/π)  (G5, arch 源码验证)
      log(h₁) = log(σ²₀)  (backcast)
      zₜ = εₜ/√hₜ  (标准化残差, G6)
      log(hₜ) = ω + β·log(hₜ₋₁) + α·zₜ₋₁ + γ·(|zₜ₋₁| - E|z|)  (Nelson 1991, G6)
      hₜ = exp(log(hₜ))  (指数变换保证正定性)
  2.2 计算对数似然 ℓ(θ) (同 GARCH(1,1) Step 2.2)
  2.3 返回 -ℓ(θ)

Step 3: 定义约束
  3.1 EGARCH 无方差非负约束 (log 变换自动保证 hₜ > 0)
  3.2 平稳性约束: |β| < 1 (EGARCH 平稳条件不同于 GARCH)
  3.3 边界: ω ∈ [-10, 10], α ∈ [-1, 1], β ∈ [0, 1), γ ∈ [-1, 1], ν ∈ [2.05, 500.0]
  注 (排幻觉点 G-nu-bound): ν 边界统一用 [2.05, 500.0] (arch StudentsT 默认, `distribution.py:505,508` 实测),
    非 [2.1, 100] (旧版); GARCH(1,1)/EGARCH/GJR-GARCH 三模型一致

Step 4-7: 同 GARCH(1,1) (多起始/sandwich/诊断/返回)
```

**测试矩阵**:

| 用例 | 验证点 | 容差 |
|------|--------|------|
| EGARCH Normal 参数 vs arch (vol='EGARCH') | ω/α/β/γ 一致 | 1e-10 |
| E\|z\| = √(2/π) vs arch 源码 | G5 常数 | 1e-15 |
| EGARCH 非对称项用 zₜ₋₁ (非 εₜ₋₁) | G6 | 1e-10 |
| EGARCH 杠杆效应 (γ<0 时负冲击放大波动) | 方向性 | - |
| EGARCH 条件方差序列 vs arch | hₜ 一致 | 1e-10 |
| EGARCH sandwich 协方差 vs arch (robust) | V 一致 | 1e-10 |
| EGARCH 符号约定 (θ vs -θ, G23) | arch/rugarch 一致 | 1e-10 |
| t-EGARCH ν 联合估计 vs arch | ν 一致 | 1e-8 |
| EGARCH 平稳性 \|β\|<1 检测 | 约束生效 | - |
| EGARCH 多步预测 | 收敛行为 | 1e-6 |

---

#### 2.0.4 `gjr_garch_model.hpp` - GJR-GARCH QMLE

**教材锚点**: GJR 1993 / Tsay 3ed Ch 6 / arch 源码 `arch/univariate/mean.py` (`arch_model()` o=1)

**接口签名**:

```cpp
// GJR-GARCH 参数 (Glosten-Jagannathan-Runkle 1993)
struct GjrGarchParams {
    Real omega;   ///< 常数项 ω
    Real alpha;   ///< ARCH 系数 α
    Real gamma;   ///< 非对称系数 γ (杠杆效应, G8)
    Real beta;    ///< GARCH 系数 β
    Real nu;      ///< t/GED 自由度
};

// GJR-GARCH 估计结果 (结构同 GarchResult)
struct GjrGarchResult {
    GjrGarchParams params;
    // ... (同 GarchResult 其他字段)
};

// GJR-GARCH QMLE 估计
GjrGarchResult estimate_gjr_garch(const std::vector<Real>& data,
                                   const GarchConfig& config = GarchConfig{});
```

**幻觉点映射**:
- **G7**: 非对称指示函数 `I(zₜ<0)` 与 `I(εₜ<0)` 数学等价 (σ>0), arch 用 `I(ε<0)`, 影响级别"低"
- **G8**: `hₜ = ω + α·ε²ₜ₋₁ + γ·I(zₜ₋₁<0)·ε²ₜ₋₁ + β·hₜ₋₁`; γ>0 表杠杆效应 (负冲击放大波动)
- **G10**: 平稳性条件 `α + γ/2 + β < 1` (非 `α+γ+β<1`, 因 E[I(z<0)]=1/2), arch 约束矩阵 `a[k+1,p+1:p+o+1]=-0.5` 已验证
- **G3/G9/G11/G16/G17**: 同 GARCH(1,1)

**算法步骤** (GJR-GARCH QMLE):

```
Step 1: 预处理 (同 GARCH(1,1))

Step 2: 定义目标函数
  2.1 递归计算条件方差:
      h₁ = σ²₀ (backcast)
      zₜ = εₜ/√hₜ
      Iₜ = (zₜ < 0) ? 1 : 0  (G7: 等价于 εₜ<0)
      hₜ = ω + α·ε²ₜ₋₁ + γ·Iₜ₋₁·ε²ₜ₋₁ + β·hₜ₋₁  (G8)
  2.2 计算对数似然 (同 GARCH(1,1))
  2.3 返回 -ℓ(θ)

Step 3: 定义约束 (G4, G10)
  3.1 不等式约束:
      c₁ = ω       (ω > 0)
      c₂ = α       (α >= 0)
      c₃ = β       (β >= 0)
      c₄ = 1 - α - γ/2 - β - 1e-8  (G10: α + γ/2 + β < 1)
  注 (排幻觉点 G-gamma-sign): γ 无 ≥ 0 约束! GJR 1993 原文允许 γ 为负 (反向杠杆效应);
    arch 实践中 γ 边界为 [-1, 1], 不强制方向. 仅平稳性约束 α+γ/2+β<1 限制 γ 上界.
  3.2 边界: ω ∈ [1e-8, 100], α ∈ [0, 1], γ ∈ [-1, 1], β ∈ [0, 1], ν ∈ [2.05, 500.0]
  注 (排幻觉点 G-nu-bound): ν 边界统一用 [2.05, 500.0] (arch StudentsT 默认), 非 [2.1, 100]

Step 4-7: 同 GARCH(1,1)
```

**测试矩阵**:

| 用例 | 验证点 | 容差 |
|------|--------|------|
| GJR-GARCH 参数 vs arch (o=1) | ω/α/γ/β 一致 | 1e-10 |
| GJR 平稳性 α+γ/2+β<1 vs arch 约束矩阵 | G10 | 1e-12 |
| GJR 杠杆效应 (γ>0 时负冲击放大波动) | G8 方向性 | - |
| I(z<0) ≡ I(ε<0) (G7, σ>0) | 等价性 | 1e-15 |
| GJR 条件方差 vs arch | hₜ 一致 | 1e-10 |
| GJR sandwich 协方差 vs arch | V 一致 | 1e-10 |
| t-GJR ν 联合估计 vs arch | ν 一致 | 1e-8 |

---

#### 2.0.5 `garch_forecast.hpp` - 多步方差预测

**教材锚点**: Tsay 3ed Ch 5 / arch `arch/univariate/base.py:991` (forecast 抽象方法) + `arch/univariate/volatility.py` (具体实现)

**接口签名**:

```cpp
// GARCH(1,1) 多步方差预测 (G13)
std::vector<Real> forecast_garch11(const GarchParams& params,
                                    Real last_variance, Real last_residual,
                                    Size horizon);

// EGARCH 多步方差预测 (Nelson 1991)
// 注 (arch `arch/univariate/volatility.py:2755-2774` 实测): arch EGARCH 仅支持 horizon=1 解析预测,
// 多步必须用 simulation 方法; 本 spec 提供闭式递归 (对 log(h) 递归再 exp),
// 需在 verify_egarch.py 中用 arch simulation 方法对照 (非 analytic)
std::vector<Real> forecast_egarch(const EGarchParams& params,
                                   Real last_log_variance, Real last_z,
                                   Size horizon);

// GJR-GARCH 多步方差预测
// 注: E[I(z<0)] = 1/2, 预测时非对称项折半
std::vector<Real> forecast_gjr(const GjrGarchParams& params,
                                Real last_variance, Real last_residual,
                                Size horizon);
```

**幻觉点映射**:
- **G13**: GARCH(1,1) k 步预测 `h_{T+k} = ω·(1-φ^{k-1})/(1-φ) + φ^{k-1}·h_{T+1}` (φ=α+β, k≥1), 收敛到 `σ̄²=ω/(1-α-β)`
  - **排幻觉点**: 指数为 `k-1` (非 `k`), 起点为 `h_{T+1}` (非 `h_T`)
  - 错误公式 `ω·(1-φ^k)/(1-φ) + φ^k·h_T` 混淆了起点与递归次数: h_{T+1} 含 ε²_T 信息, 不能用 φ·h_T 替代
  - 等价形式: `h_{T+k} = σ̄² + φ^{k-1}·(h_{T+1} - σ̄²)` (指数衰减到无条件方差)

**算法步骤**:

```
GARCH(1,1) 预测 (G13):
  Step 1: h_{T+1} = ω + α·ε²_T + β·h_T  (已知 ε²_T, 非递归)
  Step 2: k≥2 时, E[ε²_{T+k-1}|F_T] = h_{T+k-1} (zₜ~N(0,1), E[z²]=1)
          h_{T+k} = ω + (α+β)·h_{T+k-1}  (一阶递归, φ=α+β)
  Step 3: 闭式解 (以 h_{T+1} 为起点, 递归 k-1 次):
          h_{T+k} = ω·(1-φ^{k-1})/(1-φ) + φ^{k-1}·h_{T+1}
          等价: h_{T+k} = σ̄² + φ^{k-1}·(h_{T+1} - σ̄²)
  Step 4: k→∞: h_{T+k} → σ̄² = ω/(1-α-β)  (无条件方差)

  实施建议: 直接递归计算 (arch `arch/univariate/volatility.py` forecast 方法实现), 闭式公式仅用于验证收敛性
```

---

#### 2.0.6 `garch_diagnostics.hpp` - 标准化残差诊断

**教材锚点**: Tsay 3ed Ch 5 / Phase 7A `volatility_diagnostics.hpp`

**接口签名**:

```cpp
#pragma once
#include <vector>
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/inference/volatility_diagnostics.hpp"  // Phase 7A 复用

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace garch {

// GARCH 标准化残差诊断结果 (聚合多个子检验)
struct GarchDiagnosticsResult {
    // L1 模型内诊断
    econometrics::JarqueBeraResult jb_test;           ///< 正态性 (Bootstrap 临界值, G-ADR4)
    econometrics::LjungBoxResult lb_z;                ///< zₜ 自相关 (G11)
    econometrics::LjungBoxResult lb_z_squared;        ///< z²ₜ ARCH 效应 (G12)
    
    // 诊断结论
    bool passes_normality;                            ///< JB 检验是否通过 (p > 0.05)
    bool passes_no_autocorr;                          ///< LB zₜ 检验是否通过
    bool passes_no_arch_effect;                       ///< LB z²ₜ 检验是否通过
    std::string summary;                              ///< 诊断摘要
};

// 执行 GARCH 标准化残差诊断
// @param std_residuals 标准化残差 zₜ = εₜ/√hₜ
// @param n_bootstrap Bootstrap 次数 (G-ADR4, 0 表示用渐近临界值)
// @param lb_lag LB 滞后数 (0 表示自动: floor(log(T)) 或 Schwert, G12)
// @return 诊断结果
GarchDiagnosticsResult diagnose_garch_residuals(
    const std::vector<Real>& std_residuals,
    Size n_bootstrap = 1000,
    Size lb_lag = 0);

}  // namespace garch
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
```

**复用说明**: 直接调用 Phase 7A `volatility_diagnostics.hpp` 的 `jarque_bera_test()` / `ljung_box_test()` 函数, 不重复实现。Bootstrap 临界值复用 v1.5 `block_bootstrap.hpp`。

**幻觉点映射**:
- **G11**: `zₜ = εₜ/√hₜ` (非 εₜ/σₜ, 标准差符号陷阱)
- **G12**: z²ₜ LB 滞后用 `floor(log(T))` 或 Schwert 自适应 (非固定 10)
- **G-ADR4**: JB 检验用 Bootstrap 临界值 (非渐近 χ², 小样本过度拒绝)

---

## 3. M2: 单位根与方差比检验 (Week 2-4, 与 M1 并行)

> **定位**: 伪回归检测基础设施 — 因子值 (PE/PB) 单位根检验
> **教材锚点**: Hamilton 1994 Ch 17 / Dickey-Fuller 1979 / Phillips-Perron 1988 / KPSS 1992 / ERS 1996 / Lo-MacKinlay 1988
> **对照库**: Python `arch.unitroot` (主基准, U-ADR9)
> **命名空间**: `cpphub::v1::timeseries::unit_root` (ADR-017)

### 3.0 公共基础设施

#### 3.0.1 `unit_root_common.hpp` - 共享工具

**接口签名**:

```cpp
#pragma once
#include <vector>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace unit_root {

// Schwert lag 选择规则 (U-ADR1, U1)
// ceil(12·(T/100)^0.25), arch 默认 (向上取整, 非 floor)
// 排幻觉点 U1-round (arch `unitroot.py:378-380, 1136, 1326` 实测):
//   arch 3 处实现 (ADF/PP/KPSS legacy) 全部用 int(ceil(12.0 * (nobs/100)**0.25)),
//   非 floor; ADF/DFGLS 有上限 min(lag, max((T-1)//2 - 1, 0)), PP/KPSS legacy 无上限
// @param T 样本量
// @return Schwert lag 数
Size schwert_lag(Size T);

// AIC/BIC lag 选择 (备选, U-ADR1)
// @param residuals 残差序列
// @param max_lag 最大滞后
// @param criterion "AIC" 或 "BIC"
// @return 最优 lag
Size select_lag_by_ic(const std::vector<Real>& residuals,
                       Size max_lag, const std::string& criterion);

// 自动方程形式选择 (U-ADR2, U2)
// 基于趋势显著性检验选择 nc/c/ct
// @param data 原始序列
// @return "nc" (无常数无趋势) / "c" (常数) / "ct" (常数+趋势)
std::string select_trend_spec(const std::vector<Real>& data);

// 长期方差估计 (KPSS/PP 共用, U-ADR7/U-ADR8)
// @param residuals 残差序列
// @param bandwidth 带宽 (0 表示自动选择)
// @param kernel 核函数 ("Bartlett" / "QS" / "Parzen")
// @return 长期方差 σ²
Real long_run_variance(const std::vector<Real>& residuals,
                        Size bandwidth = 0,
                        const std::string& kernel = "Bartlett");

}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
```

**幻觉点映射**:
- **U1**: Schwert 规则 `ceil(12·(T/100)^0.25)` (向上取整, 非 floor; arch `unitroot.py:378-380, 1136, 1326` 三处实测一致), arch 默认 (非立方根 `floor(T^(1/3))`)
- **U2**: 自动方程形式基于趋势显著性 (非手动指定)
- **U11/U5**: 长期方差带宽自动选择 (KPSS 用 Hobijn et al. 1998 数据依赖法, PP 用 Schwert 规则; 均用 Bartlett 核)

---

#### 3.0.2 `mackinnon_cv.hpp` - MacKinnon 2010 临界值

**教材锚点**: MacKinnon 2010 *JBES* / arch `arch/unitroot/critical_values/dickey_fuller.py`

**接口签名**:

```cpp
#pragma once
#include <vector>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace unit_root {

// MacKinnon 2010 response surface 临界值 (U-ADR3, U3/U13)
// 公式 (arch `unitroot.py:1919-1935` + `arch/unitroot/critical_values/dickey_fuller.py` 实测):
//   CV(p) = β_∞ + β_0·T^(-1) + β_1·T^(-2) + β_2·T^(-3)
//   4 系数 (含 β_∞), 3 次多项式, 非 5 系数 4 次
//   排幻觉点 U3-coef: arch tau_2010 数组每行 4 元素 (如 [-3.43, -6.54, -16.79, -79.43]),
//     polyval([c3,c2,c1,c0], 1/T) = c3/T³ + c2/T² + c1/T + c0
//   注: MacKinnon 1996/2010 的 5 系数 4 次多项式仅用于协整检验 (N≥2), ADF N=1 不用
// @param test_type 检验类型 ("adf" / "pp" / "df_gls")
// @param trend_spec 方程形式 ("nc" / "c" / "ct")
// @param T 样本量
// @param n_params 模型参数数 (通常 1, ADF 的 α)
// @param p 显著性水平 (0.01 / 0.05 / 0.10)
// @return 临界值
// @throws std::invalid_argument 若 test_type/trend_spec/p 不支持
Real mackinnon_critical_value(const std::string& test_type,
                               const std::string& trend_spec,
                               Size T, Size n_params, Real p);

// MacKinnon 2010 p 值 (response surface)
// @param statistic 检验统计量
// @param test_type / trend_spec / T / n_params 同上
// @return p 值
Real mackinnon_p_value(Real statistic, const std::string& test_type,
                        const std::string& trend_spec,
                        Size T, Size n_params);

}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
```

**幻觉点映射**:
- **U3**: MacKinnon 2010 response surface: **4 系数 (β_∞ + 3 个 T 负幂项), 3 次多项式** (非 5 系数 4 次, arch `arch/unitroot/critical_values/dickey_fuller.py` 实测每行 4 元素; 非 1996 旧版查表)
- **U13**: ADF/PP p 值用 MacKinnon 2010 (非 1996 旧系数)
- **U7**: DF-GLS demean (c) 和 trend (ct) 临界值来自 arch 独立模拟 (非 ERS 1996 原表)
- **U21**: MacKinnon 1994 vs 2010 版本区分 (arch `arch/unitroot/unitroot.py:832` DFGLS 类 docstring 明确引用)

**算法步骤**:

```
MacKinnon 2010 response surface (arch `unitroot.py:1919-1935` 实测):
  CV(p) = β_∞ + β_0/T + β_1/T² + β_2/T³    ← 4 系数, 3 次多项式 (非 5 系数 4 次)

  系数表 (arch `arch/unitroot/critical_values/dickey_fuller.py`, 按 test_type × trend_spec × p 索引):
  - ADF nc/c/ct 三套独立系数 (U3: 不能混用), 每套每行 4 元素
  - PP 共享 ADF 系数 (同为 DF 分布)
  - DF-GLS c/ct 独立系数 (U7: arch `arch/unitroot/critical_values/dfgls.py` 独立模拟, 非 ERS 原表), 每行 4 元素
  - 注: 5 系数 4 次多项式仅用于协整检验 (N≥2), ADF N=1 不用

p 值计算:
  1. 构造 MacKinnon 2010 的 p 值 response surface (额外系数)
  2. 给定 statistic, 反向插值得到 p
```

---

### 3.1 ADF 检验 (`adf_test.hpp`)

**教材锚点**: Dickey-Fuller 1979 / Said-Dickey 1984 / Hamilton 1994 Ch 17

**接口签名**:

```cpp
struct ADFResult {
    Real statistic;                ///< ADF τ 统计量 (U4: 非标准分布, 非 Student-t)
    Real p_value;                  ///< p 值 (MacKinnon 2010, U3)
    Real critical_value_1pct;      ///< 1% 临界值
    Real critical_value_5pct;      ///< 5% 临界值
    Real critical_value_10pct;     ///< 10% 临界值
    Size n_lags;                   ///< 使用的滞后数 (Schwert, U1)
    std::string trend_spec;        ///< 方程形式 "nc"/"c"/"ct" (U2)
    std::vector<Real> coefficients;  ///< 回归系数
    std::vector<Real> std_errors;    ///< 标准误
    Size n_obs;                    ///< 有效观测数
    bool reject_null;              ///< 是否拒绝 H0 (单位根)
    std::string summary;           ///< 摘要
};

// ADF 检验 (U-ADR1/U-ADR2/U-ADR3)
// @param data 原始序列 (价格/因子值)
// @param trend_spec 方程形式 ("auto" / "nc" / "c" / "ct", U-ADR2)
// @param max_lag 最大滞后 (0 表示 Schwert 自动, U-ADR1)
// @param use_aic_bic 是否用 AIC/BIC 选择 lag (false => Schwert)
// @return ADF 检验结果
ADFResult adf_test(const std::vector<Real>& data,
                    const std::string& trend_spec = "auto",
                    Size max_lag = 0,
                    bool use_aic_bic = false);
```

**幻觉点映射**:
- **U1**: Schwert 规则 `ceil(12·(T/100)^0.25)` (向上取整, 非 floor; arch `unitroot.py:378-380` 实测)
- **U2**: 自动方程形式基于趋势显著性 (arch 默认)
- **U3**: MacKinnon 2010 临界值 (非 1996 旧版)
- **U4**: ADF τ 统计量非标准分布 (非 Student-t)
- **U18**: 多个单位根检验需多重检验修正 (复用 Phase 7A `multiple_test_correction`)

**算法步骤**:

```
Step 1: 方程形式选择 (U2, U-ADR2)
  if trend_spec == "auto":
    1.1 对 data 回归 Δyₜ = α + β·t + γ·yₜ₋₁ + Σ δᵢ·Δyₜ₋ᵢ + εₜ (ct 形式)
    1.2 检验 β 显著性 (t 检验)
    1.3 若 β 显著 (p < 0.05), 选 "ct"; 否则检验 α 显著性
    1.4 若 α 显著, 选 "c"; 否则选 "nc"
  else: 使用指定 trend_spec

Step 2: Lag 选择 (U1, U-ADR1)
  if max_lag == 0:
    2.1 max_lag = schwert_lag(T) = ceil(12·(T/100)^0.25)  (向上取整)
    2.2 max_lag = min(max_lag, max((T-1)//2 - 1, 0))  (arch 上限保护)
  if use_aic_bic:
    2.3 对 lag = 0, 1, ..., max_lag 估计 ADF 回归, 选 AIC/BIC 最小者
  else:
    2.4 直接使用 max_lag

Step 3: ADF 回归
  3.1 构造回归变量:
      Δyₜ = α + β·t + γ·yₜ₋₁ + Σ_{i=1}^{p} δᵢ·Δyₜ₋ᵢ + εₜ
      (nc: α=β=0; c: β=0; ct: 全保留)
  3.2 OLS 估计, 得到 γ̂ 及其 t 统计量 (ADF τ)
  3.3 τ = γ̂ / SE(γ̂)

Step 4: 临界值与 p 值 (U3, U-ADR3)
  4.1 调用 mackinnon_critical_value("adf", trend_spec, T, 1, p)
  4.2 调用 mackinnon_p_value(τ, "adf", trend_spec, T, 1)

Step 5: 判断
  5.1 H0: γ = 0 (单位根, 非平稳)
  5.2 H1: γ < 0 (平稳)
  5.3 reject_null = (τ < critical_value_5pct)
```

---

### 3.2 PP 检验 (`pp_test.hpp`)

**教材锚点**: Phillips-Perron 1988

**接口签名**:

```cpp
struct PPResult {
    Real statistic;                ///< PP τ 统计量
    Real p_value;                  ///< p 值 (MacKinnon 2010)
    Real critical_value_1pct;
    Real critical_value_5pct;
    Real critical_value_10pct;
    Size bandwidth;                ///< NW 带宽 (U-ADR8)
    std::string trend_spec;        ///< "nc"/"c"/"ct"
    bool reject_null;
    std::string summary;
};

// PP 检验 (U-ADR8, U5/U6)
PPResult pp_test(const std::vector<Real>& data,
                  const std::string& trend_spec = "auto",
                  Size bandwidth = 0);  // 0 => NW automatic
```

**幻觉点映射**:
- **U5**: PP 带宽用 Schwert 规则 `ceil(12·(T/100)^0.25)` (与 ADF 相同, 非 `floor(4·(T/100)^(2/9))`; arch `unitroot.py:1135-1136` 实测)
- **U6**: PP 用 Bartlett 核 (通过 `cov_nw`, 与 KPSS 相同; arch `unitroot.py:1156` 实测)
- **U3/U4**: 同 ADF (MacKinnon 2010, 非标准分布)

**算法步骤**:

```
Step 1: 方程形式选择 (同 ADF Step 1)

Step 2: 带宽选择 (U5, U-ADR8)
  if bandwidth == 0:
    2.1 Schwert 规则: bandwidth = ceil(12·(T/100)^0.25)  (与 ADF lag 相同公式, arch 默认)
    2.2 注: PP 无 ADF 的 min 上限保护 (arch `unitroot.py:1135-1136` 直接赋值)

Step 3: PP 回归 (无滞后差分项, 仅 OLS)
  3.1 Δyₜ = α + β·t + γ·yₜ₋₁ + εₜ (ct 形式)
  3.2 OLS 估计, 得到 γ̂ 及其 t 统计量 t_γ
  3.3 计算残差 ε̂ₜ

Step 4: PP 修正 (Phillips-Perron 1988, arch `unitroot.py:1164-1166` 实测)
  4.1 长期方差: σ² = long_run_variance(ε̂, bandwidth, "Bartlett")  (U6, 通过 cov_nw)
  4.2 短期方差 (排幻觉点 U5-sigma2, df-corrected):
      s² = SSR / (n - k)  (df-corrected, n=样本量, k=回归参数数; arch `unitroot.py` 实测)
      注: 非 mean(ε̂²) = SSR/n (无 df 修正); 两者差 n/(n-k) 因子, 小样本下影响 1e-10 容差
  4.3 PP Z(tau) 统计量 (Phillips-Perron 1988 eq.11, Hamilton 1994 Ch17 eq.17.6.17):
      Z(tau) = (s/σ) · t_γ - (T·SE(γ̂)·(σ² - s²)) / (2·σ·s)
      排幻觉点 (三次修正): 第二项分母为 **2·σ·s** (长期 std × 短期 df-corrected std),
        非 2·σ² (前次修正仍错), 非 2·γ̂ (原始错误), 非 2·σ·σ_ε (σ_ε=√(SSR/n) 无 df 修正)
      arch 源码 (`unitroot.py:1164-1166`):
        `sqt(gamma0/lam2)*((rho-1)/sigma) - 0.5*((lam2-gamma0)/lam)*(n*sigma/s)`
        其中 gamma0=σ²_ε=SSR/n (无 df 修正, 用于第一项 ratio),
        lam2=σ² (长期方差), lam=σ (长期 std), s=df-corrected σ_ε=√(SSR/(n-k)) (用于第二项),
        sigma=SE(γ̂)
      注: arch 第一项用 gamma0 (SSR/n), 第二项用 s (SSR/(n-k)); 本 spec 统一用 s (df-corrected)
        以避免两套定义, 大样本下 n/(n-k)→1 数值一致; verify 脚本需核对 arch 两项差异

Step 5: 临界值与 p 值 (同 ADF Step 4, 共享 MacKinnon 系数)
```

---

### 3.3 KPSS 检验 (`kpss_test.hpp`)

**教材锚点**: KPSS 1992

**接口签名**:

```cpp
struct KPSSResult {
    Real statistic;                ///< KPSS LM 统计量
    Real critical_value_1pct;
    Real critical_value_5pct;
    Real critical_value_10pct;
    Size bandwidth;                ///< Hobijn et al. 1998 带宽 (U-ADR7, 基于 Bartlett 核)
    std::string trend_spec;        ///< "c" (level stationary) / "ct" (trend stationary)
    bool reject_null;              ///< 是否拒绝 H0 (平稳性)
    std::string summary;           ///< 明确标注 H0/H1 方向 (U12)
};

// KPSS 检验 (U-ADR6/U-ADR7, U11/U12)
// H0: 平稳性 (level 或 trend stationary)  ← 注意: 与 ADF 相反!
// H1: 单位根 (非平稳)
// 注: 长期方差用 Bartlett 核 (cov_nw), 与 PP 相同 (非 QS 核; arch `unitroot.py:1336` 实测)
KPSSResult kpss_test(const std::vector<Real>& data,
                      const std::string& trend_spec = "c",  // "c" or "ct"
                      Size bandwidth = 0);  // 0 => Hobijn et al. 1998 自动 (Bartlett 核)
```

**幻觉点映射**:
- **U10**: KPSS 统计量非标准分布 (需查 KPSS 1992 原表或 MC 模拟, 非 χ² 近似)
- **U11**: 带宽用 Hobijn et al. 1998 数据依赖法 (arch KPSS `_autolag` 默认, 基于 Bartlett 核假设; legacy `lags=-1` 用 Schwert 规则)
- **U12**: **H0: 平稳性** (与 ADF 相反), API 必须明确标注
- **U11-kernel (排幻觉点)**: KPSS 用 **Bartlett 核** (通过 `cov_nw`), **非 QS 核**; 与 PP 完全相同 (arch `unitroot.py:1336, 1352-1357` 实测 docstring 明确 "Assumes Bartlett / Newey-West kernel")

**算法步骤**:

```
Step 1: 退化趋势 (detrending)
  1.1 若 trend_spec == "c": yₜ = α + uₜ (OLS 退化均值)
  1.2 若 trend_spec == "ct": yₜ = α + β·t + uₜ (OLS 退化趋势)
  1.3 得到残差 ûₜ

Step 2: 累积残差
  2.1 Sₜ = Σᵢ₌₁ᵗ ûᵢ  (部分和)

Step 3: 长期方差估计 (U11, U-ADR7)
  3.1 if bandwidth == 0: Hobijn et al. 1998 数据依赖法自动选择 (arch KPSS `_autolag` 默认)
      legacy 模式 (lags=-1): Schwert 规则 ceil(12·(T/100)^0.25)
  3.2 σ² = long_run_variance(û, bandwidth, "Bartlett")  (KPSS 用 Bartlett 核, 与 PP 相同; 非 QS 核)

Step 4: KPSS 统计量
  4.1 LM = (1/T²) · Σₜ Sₜ² / σ²  (KPSS 1992 公式)

Step 5: 临界值 (U10)
  5.1 KPSS 1992 原表 (或 arch 内置 MC 模拟表)
  5.2 H0: 平稳性, reject_null = (LM > critical_value_5pct)
```

---

### 3.4 DF-GLS 检验 (`df_gls_test.hpp`)

**教材锚点**: Elliott-Rothenberg-Stock 1996 / arch `arch/unitroot/unitroot.py:832` (DFGLS 类)

**接口签名**:

```cpp
struct DFGlsResult {
    Real statistic;                ///< DF-GLS τ 统计量
    Real p_value;                  ///< p 值 (MacKinnon 1994/2010, U21)
    Real critical_value_1pct;
    Real critical_value_5pct;
    Real critical_value_10pct;
    Size n_lags;                   ///< AIC 选择的 lag (arch 默认 method="aic", 非 MAIC, U20-修正)
    std::string trend_spec;        ///< "c" (demean) / "ct" (trend)
    Real c_bar;                    ///< c̄ 值 (-7.0 for "c" demean, -13.5 for "ct" trend, U8/U-ADR5)
    Real rho_bar;                  ///< ρ̄ = 1 + c̄/T (GLS 变换 AR(1) 系数, U9-修正)
    bool reject_null;
    std::string summary;
};

// DF-GLS 检验 (U-ADR4/U-ADR5, U7/U8/U9/U20)
// 注: c̄ 值 (ERS 1996): trend_spec=="c" => c̄=-7.0 (demean, 弱 detrending);
//     trend_spec=="ct" => c̄=-13.5 (trend, 强 detrending). 调研报告 §12 U8 为正确标注.
// 注: GLS 变换用 ρ̄=1+c̄/T (非 c̄ 直接, U9-修正, arch `unitroot.py:947-959` 实测)
// 注: lag 选择 arch 默认 AIC (非 MAIC, U20-修正, arch `unitroot.py:911-923` 实测)
DFGlsResult df_gls_test(const std::vector<Real>& data,
                         const std::string& trend_spec = "ct",  // "c" or "ct"
                         Size max_lag = 0);  // 0 => AIC 自动 (arch 默认, 非 MAIC)
```

**幻觉点映射**:
- **U7**: DF-GLS demean (c) 和 trend (ct) 临界值来自 arch 独立模拟 (非 ERS 1996 原表)
- **U8**: c̄ 值 (ERS 1996): trend_spec=="c" => c̄=-7.0 (demean, 弱 detrending); trend_spec=="ct" => c̄=-13.5 (trend, 强 detrending). arch `dfgls.py` 实现一致
- **U9**: 必须用 GLS detrending (非 OLS 代替)
- **U20**: DF-GLS lag 用 AIC (arch 默认 method="aic", 非 MAIC; Ng-Perron 2001 MAIC 仅用于 M 检验族 MZα/MZt/MSB/MPT), 非 Schwert
- **U21**: MacKinnon 1994 vs 2010 版本区分

**算法步骤**:

```
Step 1: GLS detrending (U9, ERS 1996, arch `unitroot.py:947-959` 实测)
  1.1 确定 c̄ (U8, U-ADR5, 与 arch 一致):
      if trend_spec == "c":  c̄ = -7.0   (demean, 弱 detrending)
      if trend_spec == "ct": c̄ = -13.5  (trend, 强 detrending)
  1.2 计算 AR(1) 系数 ρ̄ (排幻觉点 U9-rho, ERS 1996 §2):
      ρ̄ = 1 + c̄/T   (local-to-unity 参数, 非 c̄ 本身!)
      注: c̄=-7.0, T=100 时 ρ̄=0.93 (非 -7.0), c̄=-13.5, T=100 时 ρ̄=0.865 (非 -13.5)
  1.3 构造 GLS 变换 (arch delta_y/delta_z):
      y* = [y₁, y₂ - ρ̄·y₁, y₃ - ρ̄·y₂, ..., y_T - ρ̄·y_{T-1}]
      z* = [z₁, z₂ - ρ̄·z₁, ..., z_T - ρ̄·z_{T-1}]  (z 为趋势列: c→[1,...,1], ct→[1,t])
      注: 第一项 y₁/z₁ 不变换 (arch delta_y[0]=y[0] 保持原值)
  1.4 OLS 回归 y* on z*, 得到趋势系数 δ̂
  1.5 GLS 退化: ỹₜ = yₜ - δ̂'·zₜ

Step 2: Lag 选择 (arch `unitroot.py:911-923, 962-971` 实测)
  注 (排幻觉点 U20-method): arch DFGLS 默认 method="aic" (非 MAIC!), 调用与 ADF 完全相同的
    _df_select_lags; MAIC 用于 Ng-Perron M 检验族 (MZα/MZt/MSB/MPT), 不用于 ERS DF-GLS.
    arch 用 OLS detrending (非 GLS) 做 lag 选择, ref: Perron-Qu 2007 (Economics Letters 94:12-19).
  if max_lag == 0:
    2.1 对 OLS-detrended 序列 ỹ_ols 做 ADF lag 选择 (AIC 默认, 或 BIC/t-stat)
    2.2 max_lag = schwert_lag(T), 对 lag=0..max_lag 用 _df_select_lags 选最优
  else: 使用 max_lag

Step 3: DF-GLS 回归
  3.1 Δỹₜ = α + γ·ỹₜ₋₁ + Σ δᵢ·Δỹₜ₋ᵢ + εₜ
  3.2 OLS 估计, τ = γ̂ / SE(γ̂)

Step 4: 临界值与 p 值 (U7, U21)
  4.1 调用 mackinnon_critical_value("df_gls", trend_spec, T, 1, p)
      (arch 独立模拟表, 非 ERS 原表)
  4.2 H0: γ = 0 (单位根), reject_null = (τ < critical_value_5pct)
```

---

### 3.5 方差比检验 (`variance_ratio_test.hpp`)

**教材锚点**: Lo-MacKinlay 1988 / Chow-Denning 1993 / Campbell-Lo-MacKinlay 1997 (CLM debiased)
**对照库**: Python `arch.unitroot.VarianceRatio` (主基准, `arch/unitroot/unitroot.py:1603-1775` VarianceRatio 类) + R `vrtest::Chow.Denning` (Chow-Denning 联合检验, arch 未实现)

**接口签名**:

```cpp
struct VarianceRatioResult {
    Real vr_statistic;             ///< VR(k) = σ̂²_k / σ̂²_1 (σ̂²_k 已含 1/k 因子, U14-vr)
    Real z1_statistic;             ///< Lo-MacKinlay Z₁ (同方差, U14)
    Real z1_p_value;
    Real z2_statistic;             ///< Lo-MacKinlay Z₂ (异方差稳健, U15)
    Real z2_p_value;
    std::vector<Real> chow_denning_stats;  ///< 多个 k 的 Z₂ 统计量 (U16, 异方差稳健)
    Real chow_denning_p_value;     ///< 联合检验 p 值 (SMM 分布, U16-smm)
    Real clm_debiased_statistic;   ///< CLM debiased Z 统计量 (U17, Campbell-Lo-MacKinlay 1997)
    Real clm_debiased_p_value;
    Size k;                        ///< 检验 horizon
    bool reject_null;              ///< H0: 随机游走
    std::string summary;
};

// 方差比检验 (U-ADR10, U14-U17)
// @param data 价格序列 (非收益率!)
// @param k 检验 horizon (2, 5, 10, 20 等)
// @param use_debiased 是否用 CLM 小样本 debiased 修正 (U17, 默认 true)
// @param trend "c" (去均值, arch 默认) / "n" (零漂移, 不去均值)
// 注: arch 默认 trend="c" 强制去均值 (μ̂ = 样本平均收益); trend="n" 时 μ̂=0
VarianceRatioResult variance_ratio_test(const std::vector<Real>& data,
                                         Size k = 2,
                                         bool use_debiased = true,
                                         const std::string& trend = "c");

// 多 horizon 方差比 + Chow-Denning 联合检验
// @param k_list horizon 列表 (如 {2, 5, 10, 20})
// 注: Chow-Denning 用 Z₂ (异方差稳健) 而非 Z₁; arch 未实现 CD, 需参考 R vrtest
VarianceRatioResult variance_ratio_test_multi(
    const std::vector<Real>& data,
    const std::vector<Size>& k_list = {2, 5, 10, 20},
    bool use_debiased = true,
    const std::string& trend = "c");
```

**幻觉点映射**:
- **U14-vr (排幻觉点, 高严重性)**: VR(k) = σ̂²_k / σ̂²_1 (σ̂²_k 已含 1/k 因子, **不再除 k**)
  - arch `unitroot.py:1745` (`sigma2_q = sum(...) / (nq*q)`) + `unitroot.py:1768` (`self._vr = sigma2_q / sigma2_1`)
  - 错误公式 `VR = σ̂²_k / (k·σ̂²_1)` 在 σ̂²_k 已含 1/k 时多除一次 k, 导致 VR 偏小 k 倍
- **U14**: Z₁ = [VR(k)-1] / √[2(2k-1)(k-1)/(3kT)] (同方差, 1/T 放入方差) — 经核对正确
- **U15 (排幻觉点三重错, 高严重性)**: Z₂ 方差公式三处错误:
  - (a) 前置因子: spec 写 2, arch 实际为 **4** (`unitroot.py:1766`: `4 * (1 - k/q)**2.0`)
  - (b) δⱼ 定义: spec 写 2 阶矩自协方差, arch 实际为 **4 阶矩** (`unitroot.py:1764`):
    `δⱼ = T·Σₜ(rₜ-μ̂)²·(rₜ₋ⱼ-μ̂)² / [Σₜ(rₜ-μ̂)²]²` (CLM 2.4.43, Lo-MacKinlay 1988 原文)
  - (c) Z₂ 缺 √T 因子: spec 写 `Z₂ = [VR-1]/√θ`, arch 实际为 `Z₂ = √T·[VR-1]/√θ` (`unitroot.py:1771`)
  - 注: robust 方差 θ 是 O(1) (不含 1/T), 必须乘 √T 标准化; Z₁ 的 1/T 放入方差是代数等价, 但 Z₂ 不能照搬
- **U16 (排幻觉点, 中严重性)**: Chow-Denning 联合检验:
  - arch **未实现** CD (仅单 horizon), 需参考 R `vrtest::Chow.Denning`
  - 应使用 **Z₂ (异方差稳健)** 而非 Z₁ (金融数据标准实践, 与 spec U15 强调 Z₂ 一致)
  - 联合 p 值用 **SMM(m,∞) 分布** (Studentized Maximum Modulus), 非 Bonferroni:
    `p = 1 - [2·Φ(|CD|) - 1]^m` (m = |k_list|)
- **U17 (排幻觉点归属错误, 中-高严重性)**: arch 的 `debiased` 是 **CLM (Campbell-Lo-MacKinlay 1997) 重叠块偏差修正**, **非 Chen-Deo 2006**
  - arch `unitroot.py:1748-1752` 实测: `σ̂²_1 *= T/(T-1)`, `σ̂²_k *= T²/[(T-k+1)(T-k)]`
  - 等价于: `VR_debiased = VR · T(T-1)/[(T-k+1)(T-k)]` (乘法修正因子)
  - arch docstring 仅引用 CLM, 未提 Chen-Deo; Chen-Deo 2006 实为幂变换联合检验 (R `vrtest::Chen.Deo`), 与此无关
- **U17-demean (排幻觉点, 中严重性)**: 去均值是 arch 默认强制行为 (trend="c"), 非可选
  - arch `unitroot.py:1730-1733`: 默认 `μ̂ = (y[-1]-y[0])/(nobs-1)` = 样本平均收益
  - σ̂²_k 公式须减 `k·μ̂` (spec 原漏减): `σ̂²_k = (1/(Tk))·Σ(Rₜ^(k) - k·μ̂)²`
  - 仅 trend="n" (零漂移) 时 μ̂=0
- **U18**: 多检验需多重修正 (复用 Phase 7A `multiple_test_correction`)

**算法步骤**:

```
Step 1: 计算收益率与漂移
  1.1 rₜ = log(Pₜ) - log(Pₜ₋₁)  (对数收益率)
  1.2 T = len(r)
  1.3 漂移估计 (arch 默认 trend="c" 强制去均值):
      if trend == "c":  μ̂ = (1/T)·Σₜ rₜ  (样本平均收益, arch 用 (y[-1]-y[0])/(nobs-1) 等价)
      if trend == "n":  μ̂ = 0  (零漂移, 不去均值)

Step 2: 计算 VR(k) (U14-vr, arch `unitroot.py:1737, 1745, 1768` 实测)
  2.1 σ̂²_1 = (1/T)·Σₜ (rₜ - μ̂)²                    (1 期方差, 去均值)
  2.2 σ̂²_k = (1/(T·k))·Σₜ (Rₜ^(k) - k·μ̂)²           (k 期方差, 已含 1/k 因子, 减 k·μ̂)
      其中 Rₜ^(k) = rₜ + rₜ₋₁ + ... + rₜ₋ₖ₊₁  (k 期累积收益)
  2.3 VR(k) = σ̂²_k / σ̂²_1                            (不再除 k! σ̂²_k 已含 1/k)

Step 3: Lo-MacKinlay Z₁ (同方差, U14, 经核对正确)
  3.1 Var(VR(k)) = 2(2k-1)(k-1) / (3kT)              (1/T 放入方差)
  3.2 Z₁ = [VR(k) - 1] / sqrt(Var(VR(k)))
  3.3 p_value = 2·(1 - Φ(|Z₁|))                      (arch `unitroot.py:1774` 一致)

Step 4: Lo-MacKinlay Z₂ (异方差稳健, U15 三重修正, arch `unitroot.py:1759-1771` 实测)
  4.1 计算 δⱼ (4 阶矩, CLM 2.4.43, arch `unitroot.py:1764`):
      z2ₜ = (rₜ - μ̂)²
      scale = (Σₜ z2ₜ)²
      δⱼ = T · Σₜ z2ₜ·z2ₜ₋ⱼ / scale    (j = 1, ..., k-1)
      注: 非 2 阶矩自协方差 γ̂(j); 这是 4 阶矩, 捕捉异方差结构
  4.2 Var(VR(k))_robust = θ = 4·Σ_{j=1}^{k-1} ((k-j)/k)² · δⱼ
      (前置因子 4, 非 2; 权重 ((k-j)/k)² 有平方; θ 是 O(1), 不含 1/T)
  4.3 Z₂ = √T · [VR(k) - 1] / sqrt(θ)    (含 √T 因子! arch `unitroot.py:1771`)
  4.4 p_value = 2·(1 - Φ(|Z₂|))

Step 5: CLM debiased 小样本修正 (U17, Campbell-Lo-MacKinlay 1997, arch `unitroot.py:1748-1752` 实测)
  注: 非 Chen-Deo 2006! arch debiased 是 CLM 重叠块偏差修正.
  5.1 σ̂²_1 ← σ̂²_1 · T/(T-1)                           (Bessel 修正)
  5.2 σ̂²_k ← σ̂²_k · T²/[(T-k+1)(T-k)]                 (重叠块修正)
      等价于: VR_debiased = VR · T(T-1)/[(T-k+1)(T-k)]
  5.3 用 VR_debiased 重新计算 Z₁/Z₂ (Step 3/4)

Step 6: Chow-Denning 联合检验 (U16, 多 horizon, arch 未实现, 参考 R vrtest)
  6.1 对每个 k ∈ k_list 计算 Z₂(k)    (异方差稳健, 非 Z₁)
  6.2 CD = max_k |Z₂(k)|
  6.3 联合 p 值 = 1 - [2·Φ(|CD|) - 1]^m    (SMM(m,∞) 分布, m=|k_list|, 非 Bonferroni)
```

---

## 4. 端到端集成测试 (`test_integration_phase7b.cpp`)

**测试场景**:

| 用例 | 场景 | 验证点 |
|------|------|--------|
| GARCH→VaR 集成 | GARCH 波动率 → VaR 计算 → Backtest | GARCH-VaR 通路完整 |
| ADF→伪回归诊断 | 对 PE/PB 因子值做 ADF, 若单位根则差分 | 伪回归检测流程 |
| GARCH vs HAR 对比 | 同一收益率序列, GARCH 与 HAR 预测精度对比 | MZ/DM 检验复用 |
| 多检验多重修正 | ADF+PP+KPSS+DF-GLS 联合检验 + BH 修正 | U18 多重检验 |
| GARCH 标准化残差诊断 | GARCH 估计后 JB+LB+z²LB 全流程 | G11/G12 |

---

## 5. 依赖关系与复用

### 5.1 模块依赖图

```
M1 GARCH 族
├── garch_distribution.hpp (独立, 仅依赖 core/)
├── garch_model.hpp
│   ├── ← SLSQP (ADR-018, 已实现)
│   ├── ← garch_distribution.hpp
│   ├── ← Phase 7A volatility_diagnostics (JB/LB)
│   └── ← v1.5 block_bootstrap (JB Bootstrap 临界值)
├── egarch_model.hpp ← garch_distribution.hpp + SLSQP
├── gjr_garch_model.hpp ← garch_distribution.hpp + SLSQP
├── garch_forecast.hpp (独立)
└── garch_diagnostics.hpp ← Phase 7A volatility_diagnostics

M2 单位根与方差比
├── unit_root_common.hpp (独立, 仅依赖 core/)
├── mackinnon_cv.hpp (独立, 系数表)
├── adf_test.hpp ← unit_root_common + mackinnon_cv + detail::ols_simple
├── pp_test.hpp ← unit_root_common + mackinnon_cv (Bartlett 核, cov_nw)
├── kpss_test.hpp ← unit_root_common (Bartlett 核, 与 PP 相同; Hobijn 带宽; 临界值表)
├── df_gls_test.hpp ← unit_root_common + mackinnon_cv (AIC, 非 MAIC)
└── variance_ratio_test.hpp ← unit_root_common + Phase 7A multiple_test_correction
```

### 5.2 Eigen3 隔离分析

| 模块 | 是否需 Eigen3 | 理由 |
|------|--------------|------|
| GARCH 族 (M1) | **否** | 参数少 (3-5 个), 用 `std::vector` + SLSQP (ADR-018 已用 std::vector) |
| 单位根 (M2) | **否** | OLS 辅助回归用 `detail::ols_simple` (Gauss-Jordan, ADR-015 方案 B) |
| 方差比 (M2) | **否** | 纯序列运算, 无矩阵 |

**结论**: M1/M2 均不需 Eigen3, 可放 `cpphub_timeseries` INTERFACE 库 (header-only), 不链接 Eigen3。Eigen3 隔离边界 (ADR-013) 不受影响。

---

## 6. 幻觉点核查清单 (实施时逐点验证)

### 6.1 GARCH 族 (G1-G23, 23 项)

> **核查模式**: 参照 v1.4.0/1.4.1 R highfrequency 核查模式, 编写 `verify_X.py` 脚本, 打印 arch 源码 + 小样本手算

| ID | 影响级 | 核查方式 | 验收容差 |
|----|--------|---------|---------|
| G1 (backcast) | 高 | arch `arch/univariate/volatility.py:1161-1168` backcast 参数 | 1e-12 |
| G2 (ε₀=μ̂) | 中 | arch `arch/univariate/mean.py` residuals 初始化 | 1e-10 |
| G3 (log(2π)) | 中 | arch `arch/univariate/distribution.py` `lls = -0.5*(log(2*pi)+...)` | 1e-15 |
| G4 (参数约束) | 高 | arch 约束矩阵 `a[i,i]=1.0` | 1e-12 |
| G5 (E\|z\|=√(2/π)) | 高 | arch `arch/univariate/volatility.py` EGARCH `np.sqrt(2/np.pi)` | 1e-15 |
| G6 (z 标准化) | 高 | arch `arch/univariate/volatility.py` `errors=rng(nobs+burn)` | 1e-10 |
| G7 (I(ε<0)≡I(z<0)) | 低 | 数学等价 (σ>0) | 1e-15 |
| G8 (GJR 符号) | 高 | GJR 1993 摘要 | 方向性 |
| G9 (sandwich) | 高 | arch `cov_type='robust'` | 1e-10 |
| G10 (α+γ/2+β<1) | 高 | arch 约束矩阵 `-0.5` | 1e-12 |
| G11 (z=ε/√h) | 中 | 数学定义 | 1e-15 |
| G12 (z² LB 滞后) | 中 | floor(log(T)) 或 Schwert | 1e-8 |
| G13 (预测收敛) | 中 | σ̄²=ω/(1-α-β) | 1e-6 |
| G14 (t 联合估计) | 高 | Bollerslev 1987 RES 69(3):542 | 1e-8 |
| G15 (SLSQP) | 高 | arch 默认 SLSQP (ADR-018 已实现) | 1e-10 |
| G16 (多起始) | 中 | 4 起始点策略 | - |
| G17 (AIC/BIC 完整似然) | 中 | 含 log(2π) | 1e-10 |
| G18 (GARCH-M c·σ) | 中 | arch `form="vol"` | 1e-10 (推迟) |
| G19 (APARCH δ) | 高 | 推迟到 v1.6+ | - |
| G20 (FIGARCH 截断) | 高 | 推迟到 v1.6+ | - |
| G21 (IGARCH 非平稳) | 中 | 推迟到 v1.6+ | - |
| G22 (跨库 solver) | 中 | arch SLSQP vs rugarch solnp | 1e-8 |
| G23 (EGARCH 符号) | 高 | arch vs rugarch θ 符号核查 | 1e-10 |

**本 spec 覆盖**: G1-G17 (17 项), G18 推迟 (GARCH-M), G19-G21 推迟 (APARCH/FIGARCH/IGARCH), G22/G23 实施时核查。

### 6.2 单位根与方差比 (U1-U22, 22 项)

| ID | 影响级 | 核查方式 | 验收容差 |
|----|--------|---------|---------|
| U1 (Schwert) | 高 | arch `unitroot.py:378-380` ceil(12·(T/100)^0.25) | 1e-10 |
| U2 (自动方程) | 高 | arch `trend="auto"` | 1e-10 |
| U3 (MacKinnon 2010) | 极高 | arch `critical_values/dickey_fuller.py` 4 系数 | 1e-12 |
| U4 (非 Student-t) | 高 | DF 分布非标准 | - |
| U5 (PP Schwert 带宽) | 高 | arch `unitroot.py:1135-1136` ceil(12·(T/100)^0.25) (非 floor(4·(T/100)^(2/9))) | 1e-10 |
| U6 (PP Bartlett) | 高 | arch `unitroot.py:1156` cov_nw (Bartlett, 与 KPSS 相同) | - |
| U7 (DF-GLS arch 表) | 极高 | arch `critical_values/dfgls.py` 独立模拟 | 1e-10 |
| U8 (c̄ 值) | 高 | ERS 1996 c̄=-7.0 (c), c̄=-13.5 (ct) | 精确 |
| U9 (GLS detrending) | 高 | ERS 1996 方法, ρ̄=1+c̄/T | 1e-10 |
| U10 (KPSS 非标准) | 高 | KPSS 1992 原表/MC | 1e-6 |
| U11 (KPSS Hobijn 带宽) | 高 | arch `unitroot.py:1352-1357` Bartlett 核 (非 QS, 非 Andrews) | 1e-10 |
| U12 (KPSS H0 方向) | 极高 | H0: 平稳性 (与 ADF 相反) | - |
| U13 (MacKinnon 2010) | 中 | 非 1996 旧版 | 1e-12 |
| U14 (Z₁ 公式) | 高 | Lo-MacKinlay 1988, 1/T 放入方差 (经核对正确) | 1e-10 |
| U14-vr (VR(k) 定义) | 极高 | arch `unitroot.py:1745,1768` VR=σ̂²_k/σ̂²_1 (不再除k) | 1e-10 |
| U15 (Z₂ 异方差) | 极高 | arch `unitroot.py:1759-1771` 三重: 因子4/4阶矩δⱼ/√T | 1e-10 |
| U16 (Chow-Denning) | 高 | arch 未实现; R vrtest; 用Z₂非Z₁; SMM(m,∞)分布 p=1-[2Φ(\|CD\|)-1]^m | 1e-8 |
| U17 (CLM debiased) | 高 | arch `unitroot.py:1748-1752` CLM修正(非Chen-Deo); VR·T(T-1)/[(T-k+1)(T-k)] | 1e-8 |
| U17-demean (去均值) | 高 | arch `unitroot.py:1730-1733` 默认trend="c"强制去均值; σ̂²_k减k·μ̂ | 1e-10 |
| U18 (多重修正) | 中 | 复用 Phase 7A | 1e-10 |
| U19 (Ng-Perron M) | 高 | 推迟到 v1.7 | - |
| U20 (MAIC→AIC) | 高 | DF-GLS 用 AIC (非 MAIC; arch `unitroot.py:909` 默认 method="aic") | 1e-10 |
| U21 (MacKinnon 版本) | 中 | 1994 vs 2010 区分 | - |
| U22 (ZA 推迟) | 中 | 推迟到 v1.7 | - |

**本 spec 覆盖**: U1-U18 (18 项, 含 U14-vr/U17-demean 新增子项), U19/U22 推迟 (Ng-Perron/ZA), U20/U21 已核查 (AIC/MacKinnon 版本)。

---

## 7. 验收标准

### 7.1 编译与测试

- [ ] MSVC Release: 全量测试通过 (新增 ~201 + 现有 2004 = ~2205)
- [ ] A 站 GCC: fresh clone + rebuild + ctest 通过
- [ ] B 站 GCC: fresh clone + rebuild + ctest 通过
- [ ] 三平台 SLSQP + GARCH + 单位根测试无数值偏差 (容差 1e-10)

### 7.2 数值基准对照

- [ ] GARCH(1,1) vs Python `arch`: 参数 1e-10, 协方差 1e-10
- [ ] EGARCH vs Python `arch`: 参数 1e-10, E\|z\| 1e-15
- [ ] GJR-GARCH vs Python `arch`: 参数 1e-10, 平稳性 1e-12
- [ ] ADF/PP/KPSS/DF-GLS vs Python `arch`: 统计量 1e-10
- [ ] MacKinnon 2010 临界值: 1e-12
- [ ] 方差比 Z₁/Z₂ vs Python `arch`: 1e-10
- [ ] GARCH vs R `rugarch` 交叉验证: 1e-8

### 7.3 幻觉点覆盖

- [ ] G1-G17 (17 项) 全部核查, 编写 verify 脚本
- [ ] U1-U18 (18 项, 含 U14-vr/U17-demean) 全部核查, 编写 verify 脚本
- [ ] G22/G23 (跨库差异/符号) 实施时核查
- [ ] U20/U21 (AIC/MacKinnon 版本) 已核查, 实施时验证

### 7.4 文档对齐

- [ ] ADR-016 (18 项决策) 全部实施
- [ ] ADR-017 (命名空间 `cpphub::v1::timeseries`) 遵守
- [ ] ADR-018 (SLSQP) 复用, 不修改
- [ ] DEVELOPMENT_LOG.md 更新 Phase 7B 进度

### 7.5 Scope 边界

- [ ] 仅实现 GARCH(1,1)/EGARCH/GJR-GARCH (3 模型)
- [ ] 仅实现 ADF/PP/KPSS/DF-GLS (4 检验) + 方差比 (4 变体)
- [ ] 未引入 Eigen3 (M1/M2 均用 std::vector + Gauss-Jordan)
- [ ] 未实现 APARCH/FIGARCH/IGARCH/GARCH-M (推迟 v1.6+)
- [ ] 未实现 Ng-Perron M/Zivot-Andrews (推迟 v1.7)
- [ ] 未实现 ARIMA/MIDAS/VAR (M3/M4/v1.7)

---

## 8. 风险与缓解

| 风险 | 影响 | 缓解 |
|------|------|------|
| MacKinnon 2010 系数表录入错误 | 临界值错误, 检验结论反转 | arch `arch/unitroot/critical_values/dickey_fuller.py` 源码逐系数核对 |
| DF-GLS c/ct 临界值混用 | U7 极高影响 | 独立系数表, c/ct 分别存储 |
| GARCH 局部最优 | 参数估计次优 | 多起始点 (G-ADR6, 4 起始) |
| backcast 实现错误 | 方差初始化偏差 | arch `arch/univariate/volatility.py:1161-1168` 源码核对 |
| EGARCH 符号约定不一致 | G23 系数符号反转 | arch vs rugarch 源码核对 |
| Chen-Deo 2006 归属错误 | U17 修正因子引用错误 | arch debiased 实为 CLM 1997 修正, 非 Chen-Deo 2006; Chen-Deo 是幂变换联合检验 (R vrtest) |
| 跨库 solver 差异 | G22 参数微小差异 | arch SLSQP 为主基准, rugarch 容差放宽到 1e-8 |

---

## 9. 实施顺序建议

```
Week 1: M1 公共基础设施
  ├── garch_distribution.hpp (Normal/t/GED 似然)
  └── garch_model.hpp 骨架 (backcast + 递归方差 + SLSQP 调用)

Week 2: M1 模型实现
  ├── garch_model.hpp 完整 (sandwich + 诊断)
  ├── egarch_model.hpp (E|z| + log 方差)
  └── gjr_garch_model.hpp (非对称项 + 平稳性)

Week 3: M1 测试 + M2 并行启动
  ├── GARCH 测试 vs arch (20+18+18+12+15 用例)
  ├── mackinnon_cv.hpp (系数表录入)
  └── unit_root_common.hpp (Schwert/长期方差)

Week 4: M2 实现 + 集成
  ├── adf_test.hpp / pp_test.hpp / kpss_test.hpp / df_gls_test.hpp
  ├── variance_ratio_test.hpp (Z₁/Z₂/Chow-Denning/CLM debiased)
  └── test_integration_phase7b.cpp (端到端)

Week 5: 跨平台验证 + 验收
  ├── A/B 站 GCC fresh clone + rebuild + ctest
  └── PHASE7B_FINAL_ACCEPTANCE.md
```

---

## 10. 关联文档

- [ADR-016: 金融时间序列实施边界 (18 项)](../../decisions/ADR-016_FINANCIAL_TIMESERIES_BOUNDARY.md)
- [ADR-017: 时序模块命名空间](../../decisions/ADR-017_TIMESERIES_NAMESPACE.md)
- [ADR-018: SLSQP 优化器实现边界](../../decisions/ADR-018_SLSQP_BOUNDARY.md)
- [FINANCIAL_TIMESERIES_RESEARCH.md](../../research/FINANCIAL_TIMESERIES_RESEARCH.md) v3.2
- [PHASE7A_FALSIFICATION_SPEC.md](./PHASE7A_FALSIFICATION_SPEC.md)
- [PHASE6_ECONOMETRICS_SPEC.md](../phase6/PHASE6_ECONOMETRICS_SPEC.md)
