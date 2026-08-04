# Phase 6 执行规格书 - 经典参数计量模块

> **版本归属**: **v1.5** (v1.5 = Phase 6 完整 scope)
> **目标**: 交付研究级 C++ 经典参数计量基础设施,填补 C++ 计量生态空白
> **前置**: Phase 1-5 (v1.0-v1.4.3) 全部通过, Eigen3 依赖就绪 (ADR-013)
> **里程碑**: M1(Week 1-2) OLS + HC/HAC/Cluster → M2(Week 3-4) MLE/QMLE + 检验 → M3(Week 5-6) GMM → M4(Week 7-8) Bootstrap + 验收
>
> **Scope 声明** (2026-08-04 评审):
> - **严格聚焦参数方法** (v1.3 规划的 4 波 12 项), 不含半参数/非参数实现
> - **架构预留扩展接口**: `Estimator` 基类设计支持未来半参数/非参数派生 (v1.6+/v2.0+)
> - **与 calibration/ 分离** (ADR-014): 工程标定 vs 统计估计性质不同, 不混模块
> - **依赖 Eigen3** (ADR-013): 动态矩阵 `linalg::dynamic::MatrixXD`, 不污染定价模块
> - **教材锚点**: Greene 8ed / Wooldridge CS 2ed / Hayashi / Davidson-MacKinnon 四主教材
> - **对照库**: statsmodels / R sandwich / arch / linearmodels (跨语言数值验证 1e-8)
> - **Research OS 价值**: 为"因子失效诊断"方向提供统计推断基础设施 (GARCH MLE + Newey-West + Bootstrap)

---

## 1. 交付物清单

### 1.1 新增目录结构

```
include/cpphub/econometrics/         # 经典参数计量模块 (header-only)
├── core/                            # 模块内部基础设施
│   ├── data_types.hpp               # EconData / PanelData / TimeSeries 数据载体
│   ├── estimation_result.hpp        # EstimationResult 统一输出结构
│   ├── covariance_type.hpp          # CovarianceType 枚举 (HC0-3/HAC/Cluster/OPG/Hessian/Sandwich)
│   └── estimator_base.hpp           # Estimator 抽象基类 (扩展接口预留)
├── estimation/                      # 估计器 (ADR-014)
│   ├── ols.hpp                      # OLS + HC0-HC3 标准误差
│   ├── mle.hpp                      # MLE / QMLE (Logistic/Poisson/NB/Gaussian GLM)
│   ├── gmm.hpp                      # GMM (两步/迭代/CUE) + Arellano-Bond
│   └── estimator_factory.hpp        # 工厂 + 静态注册 (ADR-003 风格)
├── inference/                       # 推断
│   ├── standard_errors.hpp          # Hessian / Sandwich / Newey-West / Cluster
│   ├── hac_kernels.hpp              # Bartlett / Quadratic Spectral / Parzen 内核
│   ├── hypothesis_tests.hpp         # Wald / LR / LM 三大检验
│   └── diagnostics.hpp              # AIC / BIC / HQ / Overidentification J-test
├── resampling/                      # Bootstrap (避免与 monte_carlo/ 混淆)
│   ├── bootstrap_base.hpp           # BootstrapEngine 基类
│   ├── paired_bootstrap.hpp         # 配对 Bootstrap (i.i.d.)
│   ├── wild_bootstrap.hpp           # Wild Bootstrap (异方差稳健, Mammen/Rademacher)
│   ├── block_bootstrap.hpp          # Block Bootstrap (Politis-Romano stationary)
│   └── cluster_bootstrap.hpp        # Cluster Bootstrap (聚类)
└── data/                            # 数据处理
    ├── panel_data.hpp               # 面板数据 (长表/宽表转换, 平衡/非平衡)
    └── formula.hpp                  # Formula 接口 (R 风格 y ~ x1 + x2, 可选)
```

### 1.2 新增编译目标

| 目标 | 类型 | 关键源文件 |
|------|------|------------|
| `cpphub` (增量) | Shared Library | header-only, 无新增 .cpp (除可能的 Eigen3 实例化) |
| `cpphub_econometrics` | Interface Library | CMake INTERFACE, 导出 Eigen3 依赖 |

### 1.3 必须通过的新增测试

| 测试套件 | 用例数 | 覆盖模块 |
|----------|--------|----------|
| `test_ols_hc` | 20 | `ols.hpp` (OLS + HC0-HC3, Longley/Nerlove 数据) |
| `test_hac_kernels` | 15 | `hac_kernels.hpp` (Bartlett/QS/Parzen 边界与权重) |
| `test_newey_west` | 15 | `standard_errors.hpp` (Newey-West 1987 表对照) |
| `test_cluster_se` | 10 | `standard_errors.hpp` (Liang-Zeger 聚类) |
| `test_mle_logistic` | 15 | `mle.hpp` (Logistic, Greene 表 17.x) |
| `test_mle_poisson` | 10 | `mle.hpp` (Poisson, Cameron-Trivedi 数据) |
| `test_qmle_sandwich` | 10 | `mle.hpp` (QMLE + Sandwich, White 1982) |
| `test_hypothesis_tests` | 20 | `hypothesis_tests.hpp` (Wald/LR/LM, Greene 表 5.x) |
| `test_gmm_two_step` | 15 | `gmm.hpp` (两步 GMM, Hayashi 表 3.x) |
| `test_gmm_cue` | 10 | `gmm.hpp` (CUE, Hansen-Heaton-Yaron 1996) |
| `test_arellano_bond` | 10 | `gmm.hpp` (动态面板, Arellano-Bond 1991) |
| `test_paired_bootstrap` | 10 | `paired_bootstrap.hpp` (Efron-Tibshirani 1993 法学院数据) |
| `test_wild_bootstrap` | 10 | `wild_bootstrap.hpp` (Mammen/Rademacher) |
| `test_block_bootstrap` | 10 | `block_bootstrap.hpp` (Politis-Romano 1994) |
| `test_cluster_bootstrap` | 5 | `cluster_bootstrap.hpp` (Cameron-Gelbach-Miller 2008) |
| `test_integration_phase6` | 10 | 端到端: OLS→HC/HAC→检验→Bootstrap |

### 1.4 必须达到的数值基准

| 基准 | 容差 | 验收方式 |
|------|------|----------|
| OLS 系数 vs Greene 表 3.x (Longley, 16 obs × 6 regressors) | 1e-10 | `test_ols_hc.cpp` |
| HC0-3 标准误差 vs R `sandwich::vcovHC` (Nerlove 145 obs) | 1e-8 | `test_ols_hc.cpp` |
| Newey-West HAC vs R `sandwich::NeweyWest` (Longley + AR(1) 合成噪声) | 1e-8 | `test_newey_west.cpp` |
| 聚类 SE vs R `sandwich::vcovCL` (Grunfeld, Greene/Baltagi 10 firms 子集) | 1e-8 | `test_cluster_se.cpp` |
| MLE Logistic vs Greene 表 17.x (Spector-Mazzeo 32 obs) | 1e-8 | `test_mle_logistic.cpp` |
| MLE Poisson vs statsmodels `Poisson` (DoctorVisits 5,190 obs) | 1e-8 | `test_mle_poisson.cpp` |
| QMLE Sandwich vs R `sandwich::sandwich` | 1e-8 | `test_qmle_sandwich.cpp` |
| Wald/LR/LM vs Greene 表 5.x | 1e-6 | `test_hypothesis_tests.cpp` |
| 两步 GMM vs Hayashi 表 3.x (Hansen-Singleton 1982 CCAPM, Bartlett/QS 内核) | 1e-6 | `test_gmm_two_step.cpp` |
| CUE vs linearmodels `IVGMMCUE` | 1e-6 | `test_gmm_cue.cpp` |
| Arellano-Bond vs Stata `xtabond` (1991 论文 abdata, 140 UK firms 1976-1984) | 1e-6 | `test_arellano_bond.cpp` |
| Bootstrap CI 覆盖率 vs Efron-Tibshirani (1993) 法学院 15 obs | 1e-3 (随机性) | `test_paired_bootstrap.cpp` |
| 跨语言 (statsmodels/R sandwich/arch) | 1e-8 | CI matrix |

---

## 2. 架构决策与扩展接口

### 2.1 模块边界 (ADR-014 实施)

```
cpphub/calibration/              # 衍生品标定 (工程, Phase 3, 已完成)
├── optimizer.hpp                # LM + DE + Nelder-Mead (共享)
├── calibrator.hpp               # HestonCalibrator / SABRCalibrator
└── objective.hpp                # VEGA 加权 / 价格加权

cpphub/econometrics/             # 统计估计 (推断, Phase 6, 本规格)
├── core/                        # 基类与数据类型
├── estimation/                  # 估计器 (OLS/MLE/GMM)
├── inference/                   # 推断 (标准误差/检验/诊断)
├── resampling/                  # Bootstrap
└── data/                        # 面板数据 / Formula

cpphub/hfecon/                   # 高频计量 (Phase 5, 已完成 v1.4.3)
└── ...
```

**共享**: `calibration/optimizer.hpp` (LM/DE 复用于 MLE 优化), `core/linalg_dynamic.hpp` (ADR-013)
**不共享**: 目标函数, 标准误差计算, 模型检验

### 2.2 双层 linalg 依赖 (ADR-013 实施)

```cpp
// core/linalg_dynamic.hpp (Phase 6 引入, 封装 Eigen3)
namespace cpphub::linalg::dynamic {
    class MatrixXD { Eigen::MatrixXd data_; /* ... */ };
    class VectorXD { Eigen::VectorXd data_; /* ... */ };
    class SparseMatrix { Eigen::SparseMatrix<double> data_; /* ... */ };

    // 计量核心操作 (委托 Eigen3 + LAPACK)
    MatrixXD svd(const MatrixXD& A);
    MatrixXD qr(const MatrixXD& A);
    VectorXD lstsq(const MatrixXD& A, const VectorXD& b);     // OLS 核心
    MatrixXD cholesky_dynamic(const MatrixXD& A);
    MatrixXD inverse_symmetric(const MatrixXD& A);            // (X'X)^{-1}
    VectorXD solve(const MatrixXD& A, const VectorXD& b);     // Ax = b
}
```

**依赖策略**:
- Eigen3 通过 vcpkg 或 submodule 引入 (header-only, ~1.5MB)
- 可选 MKL/OpenBLAS 后端 (v2.0+ 评估)
- 定价模块 (`pricing/`, `models/`, `risk/`) **不依赖** Eigen3, 保持 header-only

### 2.3 Estimator 抽象基类 (扩展接口预留)

```cpp
// core/estimator_base.hpp
#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/data_types.hpp"
#include "cpphub/econometrics/core/estimation_result.hpp"
#include "cpphub/econometrics/core/covariance_type.hpp"
#include <memory>
#include <string>

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// =============================================================================
// 估计器类别枚举 (扩展接口预留)
// =============================================================================
enum class EstimatorClass {
    Parametric,        // 参数方法 (v1.5 实现: OLS/MLE/GMM)
    Semiparametric,    // 半参数方法 (v1.6+ 预留: Yatchew/Series 估计)
    Nonparametric,     // 非参数方法 (v2.0+ 预留: Li-Racine 核估计)
    MachineLearning    // 机器学习方法 (v2.0+ 预留: DML/Lasso)
};

// =============================================================================
// Estimator 抽象基类
//
// 设计目标:
//   1. 统一所有估计器的接口 (OLS/MLE/GMM/未来半参数/非参数)
//   2. 支持协方差类型切换 (HC/HAC/Cluster/Sandwich)
//   3. 预留 isParametric()/isSemiparametric()/isNonparametric() 扩展点
//   4. 支持 clone() 虚拟构造 (ADR-002 风格, 用于工厂与多态容器)
//
// 扩展示例 (v1.6+ 半参数):
//   class SemiparametricEstimator : public Estimator {
//       bool isSemiparametric() const override { return true; }
//       EstimationResult estimate(const EconData& data) override;
//   };
// =============================================================================
class Estimator {
public:
    virtual ~Estimator() = default;

    // 核心估计接口 (所有估计器必须实现)
    virtual EstimationResult estimate(const EconData& data) = 0;

    // 估计器标识
    virtual std::string name() const = 0;
    virtual EstimatorClass estimatorClass() const { return EstimatorClass::Parametric; }

    // 扩展点: 类型查询 (默认参数方法, 半参数/非参数派生类 override)
    virtual bool isParametric() const { return true; }
    virtual bool isSemiparametric() const { return false; }
    virtual bool isNonparametric() const { return false; }

    // 协方差类型设置 (支持运行时切换标准误差计算方式)
    virtual void setCovarianceType(CovarianceType type) { cov_type_ = type; }
    CovarianceType covarianceType() const { return cov_type_; }

    // 虚拟构造 (ADR-002 风格, 用于工厂注册与多态容器)
    virtual std::unique_ptr<Estimator> clone() const = 0;

protected:
    CovarianceType cov_type_ = CovarianceType::Classical;
};

// =============================================================================
// InferenceResult 推断结果 (标准误差 + 检验)
// =============================================================================
struct InferenceResult {
    VectorXD coefficients;
    VectorXD std_errors;
    VectorXD t_statistics;
    VectorXD p_values;
    MatrixXD vcov;                    // 完整协方差矩阵
    CovarianceType cov_type;
    Real log_likelihood;
    Real r_squared;
    Real adj_r_squared;
    Size n_obs;
    Size n_params;
    Size df_residual;
};

// =============================================================================
// BootstrapResult 重采样结果
// =============================================================================
struct BootstrapResult {
    VectorXD coef_mean;
    VectorXD coef_std;
    MatrixXD coef_vcov;
    std::vector<VectorXD> bootstrap_samples;
    Real lower_ci;                    // 百分位置信区间下界
    Real upper_ci;                    // 百分位置信区间上界
    Size n_replicates;
    Size n_failed;
};

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
```

### 2.4 数据载体设计

```cpp
// core/data_types.hpp
#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/linalg_dynamic.hpp"  // ADR-013
#include <vector>
#include <string>
#include <variant>

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// =============================================================================
// EconData 统一数据载体
//
// 设计: variant 承载三种数据形态, 估计器通过 visit 模式分发
// =============================================================================
struct CrossSectionData {
    MatrixXD X;                       // N × K 设计矩阵
    VectorXD y;                       // N × 1 响应向量
    std::vector<std::string> x_names;
    std::string y_name;
};

struct PanelData {
    MatrixXD X;
    VectorXD y;
    std::vector<Index> entity_id;     // 个体标识 (1..N)
    std::vector<Index> time_id;       // 时间标识 (1..T)
    std::vector<std::string> x_names;
    std::string y_name;
    bool balanced;
};

struct TimeSeriesData {
    VectorXD y;
    MatrixXD X;
    std::vector<Real> timestamps;     // 可选时间戳 (用于 HAC lag 选择)
    std::vector<std::string> x_names;
    std::string y_name;
};

using EconData = std::variant<CrossSectionData, PanelData, TimeSeriesData>;

// 便利构造函数
CrossSectionData make_cross_section(const MatrixXD& X, const VectorXD& y);
PanelData make_panel(const MatrixXD& X, const VectorXD& y,
                     const std::vector<Index>& entity,
                     const std::vector<Index>& time);
TimeSeriesData make_time_series(const VectorXD& y, const MatrixXD& X);

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
```

### 2.5 协方差类型枚举

```cpp
// core/covariance_type.hpp
#pragma once
#include <string>

namespace cpphub {
inline namespace v1 {
namespace econometrics {

enum class CovarianceType {
    // OLS 经典
    Classical,                        // σ² (X'X)^{-1}
    HC0,                              // White (1980)
    HC1,                              // MacKinnon-White (1985) 小样本调整
    HC2,                              // MacKinnon-White (1985) leverage 调整
    HC3,                              // MacKinnon-White (1985) Jackknife
    HC4,                              // Cribari-Neto (2004) 高杠杆调整
    HC5,                              // Cribari-Neto-Souza (2007)
    // 时间序列
    HAC_Bartlett,                     // Newey-West (1987)
    HAC_QuadraticSpectral,            // Andrews (1991)
    HAC_Parzen,                       // Gallant (1987)
    // 聚类
    Cluster,                          // Liang-Zeger (1986)
    Cluster_TwoWay,                   // Cameron-Gelbach-Miller (2011) 双向聚类
    // MLE
    OPG,                              // Outer Product of Gradients
    Hessian,                          // 观测信息矩阵
    Sandwich,                         // A^{-1} B A^{-1} (White 1982 QMLE)
    // 扩展预留 (v1.6+)
    Bootstrap,                        // Bootstrap 协方差 (v1.5 第四波)
    Custom                            // 用户自定义 (扩展点)
};

std::string to_string(CovarianceType type);
CovarianceType parse_covariance_type(const std::string& s);

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
```

### 2.6 扩展接口设计说明 (半参数/非参数预留)

**预留原则**: v1.5 不实现半参数/非参数, 但接口设计支持未来派生:

| 扩展场景 | 实现方式 | 预计版本 |
|---------|---------|---------|
| **半参数回归** (Yatchew 2003) | `SemiparametricEstimator : public Estimator` override `isSemiparametric()` | v1.6+ |
| **核估计** (Li-Racine 2007) | `KernelEstimator : public Estimator` override `isNonparametric()` | v2.0+ |
| **局部线性回归** | `LocalLinearEstimator : public Estimator` | v2.0+ |
| **Series 估计** (Newey 1997) | `SeriesEstimator : public Estimator` | v1.6+ |
| **双重机器学习 DML** (Chernozhukov 2018) | `DMLEstimator : public Estimator` override `estimatorClass() → MachineLearning` | v2.0+ |
| **Lasso/Post-Double Selection** | `LassoEstimator : public Estimator` | v2.0+ |

**关键设计约束**:
- `EstimationResult` 结构通用, 不限定参数维度 (支持非参数的函数估计结果, 通过 `fitted_values` 字段)
- `CovarianceType::Custom` 允许用户注入自定义协方差计算 (Bootstrap / 半参数)
- `Estimator::clone()` 支持多态容器与工厂注册 (ADR-003 风格)

---

## 3. Week 1-2: OLS + HC/HAC/Cluster (M1)

### 3.1 OLS 估计器 (`estimation/ols.hpp`)

```cpp
#pragma once
#include "cpphub/econometrics/core/estimator_base.hpp"
#include "cpphub/econometrics/inference/standard_errors.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// =============================================================================
// OLS 估计器
//
// 教材锚点:
//   - Greene 8ed Ch.3-4 (OLS 几何与渐近)
//   - Davidson-MacKinnon Ch.2-3 (Frisch-Waugh-Lovell 定理在 Ch.2 §2.4-2.5, OLS 统计性质在 Ch.3)
//
// 算法 (Davidson-MacKinnon):
//   1. β = (X'X)^{-1} X'y
//   2. ε = y - Xβ
//   3. σ² = ε'ε / (N - K)
//   4. V_classical = σ² (X'X)^{-1}
//   5. 根据 cov_type_ 计算 HC/HAC/Cluster 协方差
// =============================================================================
class OLSEstimator : public Estimator {
public:
    EstimationResult estimate(const EconData& data) override;
    std::string name() const override { return "OLS"; }
    std::unique_ptr<Estimator> clone() const override { return std::make_unique<OLSEstimator>(*this); }

    // OLS 专属诊断
    Real computeRSquared() const;
    Real computeAdjustedRSquared() const;
    Real computeFStatistic() const;          // 联合显著性
    VectorXD computeFittedValues() const;
    VectorXD computeResiduals() const;
    MatrixXD computeProjectionMatrix() const;  // P = X(X'X)^{-1}X' (用于 HC2/HC3 leverage)

private:
    MatrixXD X_;
    VectorXD y_;
    VectorXD beta_;
    VectorXD residuals_;
    MatrixXD XtX_inv_;
};

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
```

### 3.2 HC 标准误差 (`inference/standard_errors.hpp`)

```cpp
// =============================================================================
// HC0-HC5 异方差稳健标准误差
//
// 教材锚点:
//   - White (1980) HC0 原始论文
//   - MacKinnon-White (1985) HC1/HC2/HC3 表 1 (小样本相对表现)
//   - Davidson-MacKinnon Ch.5-6 完整推导
//
// 公式 (Davidson-MacKinnon / R sandwich::vcovHC 源码核对):
//   V_HC0 = (X'X)^{-1} [Σ_i x_i x_i' ε_i²] (X'X)^{-1}
//   V_HC1 = (N/(N-K)) · V_HC0
//   V_HC2 = (X'X)^{-1} [Σ_i x_i x_i' ε_i² / (1 - h_i)] (X'X)^{-1}
//           h_i = x_i'(X'X)^{-1}x_i  (投影矩阵 P = X(X'X)^{-1}X' 的对角元素, leverage)
//   V_HC3 = (X'X)^{-1} [Σ_i x_i x_i' ε_i² / (1 - h_i)²] (X'X)^{-1}   Jackknife
//   V_HC4 = (X'X)^{-1} [Σ_i x_i x_i' ε_i² / (1 - h_i)^δ_i] (X'X)^{-1}
//           δ_i = min(N·h_i/K, 4)   (Cribari-Neto 2004)
//   V_HC5 = (X'X)^{-1} [Σ_i x_i x_i' ε_i² / (1 - h_i)^(δ_i/2)] (X'X)^{-1}
//           δ_i = min(N·h_i/K, max(4, 0.7·N·h_max/K))   (Cribari-Neto-Souza 2007, 经 Ng-Wilcox 2009 PDF 确认)
//           h_max = max(h_1, ..., h_N), R sandwich::vcovHC HC5 实现用 sqrt((1-h)^δ) 等价于 (1-h)^(δ/2)
// =============================================================================
MatrixXD compute_hc_vcov(const MatrixXD& X, const VectorXD& residuals,
                          const MatrixXD& XtX_inv, CovarianceType type);

// =============================================================================
// Newey-West HAC (Bartlett / QS / Parzen 内核)
//
// 教材锚点:
//   - Newey-West (1987) RCSTAT 表 (Bartlett)
//   - Andrews (1991) Econometrica (QS 最优)
//
// 公式 (Newey-West):
//   V_HAC = (X'X)^{-1} [Ω_0 + Σ_{l=1}^{L} K(l/L) (Ω_l + Ω_l')] (X'X)^{-1}
//   Ω_0 = Σ_t x_t x_t' ε_t²
//   Ω_l = Σ_{t=l+1}^{T} x_t x_{t-l}' ε_t ε_{t-l}
//
// 内核:
//   Bartlett (Newey-West): K(u) = 1 - |u|     for |u| ≤ 1
//   Quadratic Spectral:    K(u) = 25/(12π²u²) [sin(6πu/(5))/(6πu/(5)) - cos(6πu/5)]  for u ≠ 0
//                          K(0) = 1
//   Parzen (Gallant 1987):  K(u) = 1 - 6u² + 6|u|³   for 0 ≤ |u| ≤ 0.5
//                                   2(1-|u|)³       for 0.5 < |u| ≤ 1
//                          (对称内核, 公式含 |u| 而非 u, 否则负侧失真)
//
// Lag 选择:
//   Newey-West (1987) 论文经验法则: L = floor(4·(T/100)^{2/9})
//   R sandwich::NeweyWest 默认: 调用 bwNeweyWest() 自动带宽
//       → 基于 AR(1) 拟合的 Andrews (1991) 最优公式: L = 1.1447·[α(1)·T]^{1/3}
//       → sandwich 1.1-0 起 prewhite=FALSE 时直接使用此公式, 不再回退到 NW 经验法则
//   C++ 实现: 默认 max_lag=0 触发 bwNeweyWest 等价公式, 文档化两者差异
// =============================================================================
MatrixXD compute_hac_vcov(const MatrixXD& X, const VectorXD& residuals,
                          const MatrixXD& XtX_inv, CovarianceType type,
                          Size max_lag = 0,                    // 0 = 自动选择
                          bool prewhiten = false);             // Andrews-Monahan (1992) 预白化

// =============================================================================
// 聚类标准误差 (Liang-Zeger 1986)
//
// 教材锚点:
//   - Liang-Zeger (1986) Biometrika (原始)
//   - Cameron-Miller (2015) "A Practitioner's Guide to Cluster-Robust Inference"
//
// 公式 (Liang-Zeger):
//   V_cluster = (X'X)^{-1} [Σ_g X_g' ε_g ε_g' X_g] (X'X)^{-1}
//   X_g, ε_g 为第 g 个聚类的子矩阵/子向量
//   小样本调整: · G/(G-1) · (N-1)/(N-K)
//
// 双向聚类 (Cameron-Gelbach-Miller 2011):
//   V_twoway = V_cluster(g1) + V_cluster(g2) - V_cluster(g1∩g2)
// =============================================================================
MatrixXD compute_cluster_vcov(const MatrixXD& X, const VectorXD& residuals,
                               const MatrixXD& XtX_inv,
                               const std::vector<Index>& cluster_id,
                               bool twoway = false,
                               const std::vector<Index>& cluster_id2 = {});
```

### 3.3 HAC 内核 (`inference/hac_kernels.hpp`)

```cpp
// =============================================================================
// HAC 内核函数 (Andrews 1991 统一框架)
//
// 排幻觉点 (R sandwich 对照):
//   - R sandwich::kweights Bartlett 实测: w[l] = 1 - l/(L+1)  (非 1 - l/L)
//   - R sandwich::kweights QS 实测: 与 Andrews (1991) 公式一致
//   - lag 索引差异: R 用 0-based, 论文用 1-based, C++ 实现需明确
// =============================================================================
enum class HacKernel { Bartlett, QuadraticSpectral, Parzen, TukeyHanning };

Real kernel_weight(HacKernel kernel, Real lag_ratio);   // u = l/(L+1)
Size  select_max_lag(Size n_obs, HacKernel kernel,
                     bool andrews_optimal = false,
                     Real ar1_coef = 0.0);
```

---

## 4. Week 3-4: MLE/QMLE + 假设检验 (M2)

### 4.1 MLE 估计器 (`estimation/mle.hpp`)

```cpp
// =============================================================================
// MLE / QMLE 估计器
//
// 教材锚点:
//   - Greene 8ed Ch.14-17 (MLE + QMLE)
//   - Wooldridge CS 2ed Ch.12-13 (QMLE 渐近正态性完整证明)
//   - White (1982) QMLE 信息矩阵等式
//
// 算法 (Greene):
//   1. 指定 log-likelihood: ℓ(θ) = Σ_i log f(y_i | x_i, θ)
//   2. 数值最大化: 复用 calibration/optimizer.hpp 的 LM/BFGS
//   3. 协方差 (三种):
//      OPG:      V = (Σ g_i g_i')^{-1}              g_i = ∂log f_i/∂θ
//      Hessian:  V = -H^{-1}                        H = ∂²ℓ/∂θ∂θ'
//      Sandwich: V = A^{-1} B A^{-1}                A = -H, B = Σ g_i g_i'
//   4. QMLE: 即使分布误设, Sandwich 协方差仍渐近正确 (White 1982)
// =============================================================================
enum class MLEFamily {
    Gaussian,      // OLS 等价 (用于检验)
    Logistic,      // 二元 Logit
    Probit,        // 二元 Probit
    Poisson,       // 计数数据
    NegativeBinomial,
    Bernoulli,
    Exponential,
    Gamma,
    Custom         // 用户自定义 log-likelihood (扩展点)
};

class MLEEstimator : public Estimator {
public:
    explicit MLEEstimator(MLEFamily family);

    EstimationResult estimate(const EconData& data) override;
    std::string name() const override;
    std::unique_ptr<Estimator> clone() const override;

    // MLE 专属
    void setStartValues(const VectorXD& theta0);
    void setOptimizationConfig(const OptimizationConfig& cfg);  // 复用 calibration/optimizer.hpp
    Real logLikelihood() const;
    Real computeAIC() const;   // 2K - 2ℓ
    Real computeBIC() const;   // K·log(N) - 2ℓ
    Real computeHQ() const;    // 2K·log(log(N)) - 2ℓ  (Hannan-Quinn)

    // 扩展点: 用户自定义 log-likelihood
    using LogLikelihoodFn = std::function<Real(const VectorXD& theta, const MatrixXD& X, const VectorXD& y)>;
    using ScoreFn = std::function<VectorXD(const VectorXD& theta, const MatrixXD& X, const VectorXD& y)>;
    using HessianFn = std::function<MatrixXD(const VectorXD& theta, const MatrixXD& X, const VectorXD& y)>;
    void setCustomLogLikelihood(LogLikelihoodFn ll, ScoreFn score = nullptr, HessianFn hessian = nullptr);

private:
    MLEFamily family_;
    VectorXD theta_;
    Real log_lik_;
    OptimizationConfig opt_cfg_;
    LogLikelihoodFn custom_ll_;
    ScoreFn custom_score_;
    HessianFn custom_hessian_;

    // 内置 log-likelihood 实现
    Real gaussianLogLik(const VectorXD& theta, const MatrixXD& X, const VectorXD& y) const;
    Real logisticLogLik(const VectorXD& theta, const MatrixXD& X, const VectorXD& y) const;
    Real poissonLogLik(const VectorXD& theta, const MatrixXD& X, const VectorXD& y) const;
    // ... 其他族
};
```

### 4.2 假设检验 (`inference/hypothesis_tests.hpp`)

```cpp
// =============================================================================
// Wald / LR / LM 三大检验
//
// 教材锚点:
//   - Wooldridge CS 2ed Ch.12 (渐近等价性证明)
//   - Greene 8ed Ch.5 (显式公式)
//
// 三检验关系 (Engle 1984):
//   Wald:  无约束模型 + 线性约束 Rβ = r
//   LR:    有约束 vs 无约束 log-likelihood 差
//   LM:    有约束模型 + score 检验
//   渐近 χ²(q), q = 约束个数
//
// 公式 (Greene):
//   Wald = (Rβ̂ - r)' [R V R']^{-1} (Rβ̂ - r)
//   LR = 2 (ℓ_UR - ℓ_R)
//   LM = ε_R' X (X'X)^{-1} X' ε_R / σ²_R    ( score · I^{-1} · score' )
// =============================================================================
struct HypothesisTestResult {
    std::string test_name;             // "Wald" / "LR" / "LM"
    Real statistic;
    Real p_value;
    Size df;                           // 自由度 (约束个数 q)
    Real critical_value_95;
    Real critical_value_99;
    bool reject_null_95;
    bool reject_null_99;
};

// 线性约束 Rβ = r
HypothesisTestResult wald_test(const VectorXD& beta, const MatrixXD& vcov,
                                const MatrixXD& R, const VectorXD& r);

// LR 检验 (需要两个估计结果)
HypothesisTestResult lr_test(const EstimationResult& unrestricted,
                              const EstimationResult& restricted);

// LM 检验 (score test, 需要有约束估计的残差)
HypothesisTestResult lm_test(const MatrixXD& X, const VectorXD& residuals_restricted,
                              const MatrixXD& XtX_inv);

// 过度识别检验 (GMM J-test, Hansen 1982)
HypothesisTestResult overidentification_test(const VectorXD& moments,
                                              const MatrixXD& weighting_matrix,
                                              Size n_params);

// 信息准则
struct InformationCriteria {
    Real aic;    // 2K - 2ℓ
    Real bic;    // K·log(N) - 2ℓ
    Real hq;     // 2K·log(log(N)) - 2ℓ
};
InformationCriteria compute_information_criteria(Real log_likelihood, Size n_params, Size n_obs);
```

---

## 5. Week 5-6: GMM (M3)

### 5.1 GMM 估计器 (`estimation/gmm.hpp`)

```cpp
// =============================================================================
// GMM 估计器 (Hansen 1982)
//
// 教材锚点:
//   - Hayashi Ch.3-4 (GMM 主线, Hansen 1982 定理 3.1-3.4)
//   - Greene 8ed Ch.13 (GMM 综述)
//   - Cochrane Ch.10-11 (资产定价 GMM, SDF 框架)
//
// 算法 (Hayashi):
//   1. 定义矩条件: g_i(θ) = z_i · ε_i(θ), E[g_i] = 0
//   2. 目标函数: J(θ) = n · ḡ(θ)' W ḡ(θ), ḡ = (1/n) Σ g_i
//   3. 两步 GMM:
//      Step 1: W₁ = I, 估计 θ̂₁
//      Step 2: W₂ = Ŝ^{-1}(θ̂₁), 估计 θ̂₂   (Ŝ = 长期方差, Newey-West HAC)
//   4. 迭代 GMM: 重复 Step 2 直到收敛
//   5. CUE (Continuously Updating): θ̂_CUE = argmin_θ ḡ(θ)' Ŝ(θ)^{-1} ḡ(θ)
//      (Hansen-Heaton-Yaron 1996; 注意 Ŝ 依赖 θ, 与两步 GMM 的 Ŝ(θ̂₁) 不同)
//
// 过度识别检验:
//   J = n · ḡ(θ̂)' Ŝ^{-1} ḡ(θ̂) ~ χ²(q - k)   q = 矩条件数, k = 参数数
// =============================================================================
enum class GMMType {
    TwoStep,       // Hansen (1982) 两步
    Iterated,      // 迭代直到收敛
    CUE            // Continuously Updating (Hansen-Heaton-Yaron 1996)
};

class GMMEstimator : public Estimator {
public:
    GMMEstimator();

    EstimationResult estimate(const EconData& data) override;
    std::string name() const override { return "GMM"; }
    std::unique_ptr<Estimator> clone() const override { return std::make_unique<GMMEstimator>(*this); }

    // GMM 专属配置
    void setGMMType(GMMType type) { gmm_type_ = type; }
    void setMaxIterations(Size n) { max_iter_ = n; }
    void setConvergenceTolerance(Real tol) { tol_ = tol; }
    void setHACConfig(HacKernel kernel, Size max_lag, bool prewhiten);

    // 矩条件接口 (必须由用户指定)
    using MomentFn = std::function<VectorXD(
        const VectorXD& theta,
        const MatrixXD& instruments,
        const VectorXD& endogenous,
        const MatrixXD& exogenous)>;
    void setMomentFunction(MomentFn moment_fn);

    // 工具变量便捷接口
    void setIVSpecification(const MatrixXD& instruments, const MatrixXD& endogenous,
                             const MatrixXD& exogenous);

    // Arellano-Bond 动态面板
    EstimationResult estimateArellanoBond(const PanelData& panel,
                                           Size max_lags_instruments,
                                           bool twostep = true);

    // J-test (过度识别检验)
    HypothesisTestResult overidentificationTest() const;

private:
    GMMType gmm_type_ = GMMType::TwoStep;
    Size max_iter_ = 50;
    Real tol_ = 1e-8;
    HacKernel hac_kernel_ = HacKernel::Bartlett;
    Size hac_max_lag_ = 0;
    bool hac_prewhiten_ = false;
    MomentFn moment_fn_;

    VectorXD theta_;
    MatrixXD weighting_matrix_;
    MatrixXD vcov_;
    Real j_statistic_;
    Size n_moments_;
    Size n_params_;

    // 内部: Ŝ 长期方差 (Newey-West HAC)
    MatrixXD computeLongRunCovariance(const MatrixXD& moment_matrix) const;
};
```

### 5.2 GMM 算法详细步骤 (Hayashi §3.5)

```cpp
// =============================================================================
// 两步 GMM 实现伪代码
//
// 输入: 矩函数 g_i(θ), 数据 (Z, Y, X), 起始值 θ₀
//
// Step 1 (W₁ = I):
//   θ̂₁ = argmin_θ [ Σ_i g_i(θ) ]' [ Σ_i g_i(θ) ]
//   使用 calibration/optimizer.hpp 的 LM/BFGS
//
// Step 2 (W₂ = Ŝ^{-1}):
//   1. 计算 moment_matrix[i] = g_i(θ̂₁)   (N × q 矩阵)
//   2. 计算 Ŝ = Newey-West HAC(moment_matrix, kernel, max_lag)
//   3. W₂ = Ŝ^{-1}
//   4. θ̂₂ = argmin_θ [ Σ_i g_i(θ) ]' W₂ [ Σ_i g_i(θ) ]
//
// 协方差:
//   V = (G' W₂ G)^{-1}   G = ∂ḡ/∂θ (q × k 矩阵)
//   或更稳健: V = (G' W₂ G)^{-1} G' W₂ Ŝ W₂ G (G' W₂ G)^{-1}
//   当 W₂ = Ŝ^{-1} 时两者相等
//
// 迭代 GMM:
//   重复 Step 2 直到 ||θ̂_{k+1} - θ̂_k|| < tol
//
// CUE:
//   θ̂ = argmin_θ [ Σ_i g_i(θ) ]' Ŝ(θ)^{-1} [ Σ_i g_i(θ) ]
//   注意 Ŝ 依赖 θ, 数值优化更复杂
// =============================================================================
```

---

## 6. Week 7-8: Bootstrap (M4)

### 6.1 Bootstrap 基类 (`resampling/bootstrap_base.hpp`)

```cpp
// =============================================================================
// Bootstrap 引擎基类
//
// 教材锚点:
//   - Davison-Hinkley (1997) "Bootstrap Methods and Their Application"
//   - Cameron-Trivedi Ch.11 (配对/非参数/残差/Wild 完整方法)
//   - Politis-Romano (1994) Stationary Bootstrap
//
// 扩展接口: 支持 v1.6+ Wild Cluster Bootstrap (Cameron-Gelbach-Miller 2008)
// =============================================================================
class BootstrapEngine {
public:
    virtual ~BootstrapEngine() = default;

    // 核心重采样接口
    virtual BootstrapResult resample(const Estimator& estimator,
                                      const EconData& data,
                                      Size n_replicates = 999,
                                      std::uint64_t seed = 42) = 0;

    virtual std::string name() const = 0;

    // 配置
    void setConfidenceLevel(Real alpha) { alpha_ = alpha; }   // 默认 0.05
    void setRNG(std::uint64_t seed) { seed_ = seed; }

protected:
    Real alpha_ = 0.05;
    std::uint64_t seed_ = 42;

    // 百分位置信区间 (Efron's percentile)
    std::pair<Real, Real> percentileCI(const std::vector<Real>& samples, Real alpha) const;

    // BCa 区间 (bias-corrected and accelerated, Davison-Hinkley §5.3)
    std::pair<Real, Real> bcaCI(const std::vector<Real>& samples, Real alpha,
                                 const Estimator& estimator, const EconData& data) const;
};
```

### 6.2 四种 Bootstrap 实现

```cpp
// =============================================================================
// 配对 Bootstrap (i.i.d. 场景, Efron 1979)
//
// 算法:
//   for b in 1..B:
//     1. 从 {(x_i, y_i)} 有放回采样 N 对 → {(x*_i, y*_i)}
//     2. 在重采样数据上估计 θ̂*_b
//   3. 协方差: V_boot = (1/(B-1)) Σ (θ̂*_b - θ̄*)(θ̂*_b - θ̄*)'
// =============================================================================
class PairedBootstrap : public BootstrapEngine {
public:
    BootstrapResult resample(const Estimator& estimator, const EconData& data,
                              Size n_replicates, std::uint64_t seed) override;
    std::string name() const override { return "PairedBootstrap"; }
};

// =============================================================================
// Wild Bootstrap (异方差稳健, Wu 1986 / Liu 1988 / Mammen 1993)
//
// 算法:
//   for b in 1..B:
//     1. 生成权重 v_i (Rademacher: P(v=1)=0.5, P(v=-1)=0.5; 或 Mammen/Webb6)
//     2. y*_i = x_i'β̂ + v_i · ε̂_i   (固定 X, 残差乘权重)
//     3. 在 (X, y*) 上估计 β̂*_b
//
// 三种权重分布:
//   Rademacher: v = ±1 w.p. 0.5 (最保守, 推荐小样本, 默认分布; Cameron-Gelbach-Miller 2008 推荐)
//   Mammen (1993): v = (1-√5)/2 w.p. (5+√5)/10, v = (1+√5)/2 w.p. (5-√5)/10
//                  (满足 E[v]=0, E[v²]=1, 二阶矩匹配, 历史"默认"但已被 Rademacher 取代)
//   Webb6 (Webb 2018): 6 点分布 {-√(3/2), -1, -√(1/2), √(1/2), 1, √(3/2)} 等概率
//                  (改进小样本稳健性, 适用于少数聚类 wild bootstrap; 源: MacKinnon-Webb 2018 Econometrics Journal)
// =============================================================================
class WildBootstrap : public BootstrapEngine {
public:
    enum class WeightDistribution { Mammen, Rademacher, Webb6 };  // Webb6: 6 点分布
    explicit WildBootstrap(WeightDistribution dist = WeightDistribution::Rademacher);

    BootstrapResult resample(const Estimator& estimator, const EconData& data,
                              Size n_replicates, std::uint64_t seed) override;
    std::string name() const override { return "WildBootstrap"; }

private:
    WeightDistribution dist_;
};

// =============================================================================
// Block Bootstrap (时间序列, Politis-Romano 1994 stationary)
//
// 算法 (Politis-Romano):
//   for b in 1..B:
//     1. 随机块长度 L ~ 几何分布(参数 p = 1/expected_block_length)
//     2. 从原序列随机选起始点, 取长度 L 的块
//     3. 拼接块直到长度 N
//     4. 在重采样序列上估计 θ̂*_b
//
// 块长度选择:
//   Politis-White (2004) 自动选择: 基于自相关衰减
//   经验法则: L = N^{1/3} (Politis-Romano 1994)
// =============================================================================
class BlockBootstrap : public BootstrapEngine {
public:
    enum class BlockType { Stationary, Circular };  // Politis-Romano vs Politis-White

    explicit BlockBootstrap(BlockType type = BlockType::Stationary,
                             Size expected_block_length = 0);  // 0 = 自动选择

    BootstrapResult resample(const Estimator& estimator, const EconData& data,
                              Size n_replicates, std::uint64_t seed) override;
    std::string name() const override { return "BlockBootstrap"; }

    // Politis-White (2004) 自动块长度选择
    Size selectBlockLength(const VectorXD& series) const;

private:
    BlockType type_;
    Size block_length_;
};

// =============================================================================
// Cluster Bootstrap (聚类数据, Cameron-Gelbach-Miller 2008)
//
// 算法:
//   for b in 1..B:
//     1. 从 G 个聚类中有放回采样 G 个聚类
//     2. 拼接为重采样数据 (保留聚类内相关性)
//     3. 在重采样数据上估计 θ̂*_b
//
// 注意: 当聚类数 G < 20 时, 标准 cluster bootstrap 不可靠
//       推荐使用 Wild Cluster Bootstrap (Cameron-Gelbach-Miller 2008)
// =============================================================================
class ClusterBootstrap : public BootstrapEngine {
public:
    explicit ClusterBootstrap(const std::vector<Index>& cluster_id);

    BootstrapResult resample(const Estimator& estimator, const EconData& data,
                              Size n_replicates, std::uint64_t seed) override;
    std::string name() const override { return "ClusterBootstrap"; }

private:
    std::vector<Index> cluster_id_;
};
```

---

## 7. 测试矩阵与基准对照

### 7.1 教材数值基准

| 测试类别 | 基准来源 | 容差 | 数据集 |
|---------|---------|------|--------|
| OLS 系数 | Greene 表 3.x | 1e-10 | Longley (16 obs, 6 vars) |
| HC0-3 标准误差 | MacKinnon-White (1985) 表 1 | 1e-6 | Nerlove 电力 (145 obs) |
| HAC Bartlett | Newey-West (1987) 表 | 1e-6 | Longley + 添加 AR(1) 噪声 |
| 聚类 SE | R `sandwich::vcovCL` | 1e-8 | Grunfeld (Greene/Baltagi 10 firms × 20 years 子集) |
| MLE Logistic | Greene 表 17.x | 1e-8 | Spector-Mazzeo (32 obs) |
| MLE Poisson | Cameron-Trivedi (2013) | 1e-8 | Doctor visits (5,190 obs, 1977 Australian Health Survey) |
| QMLE Sandwich | R `sandwich::sandwich` | 1e-8 | 同上 |
| Wald/LR/LM | Greene 表 5.x | 1e-6 | Longley |
| 两步 GMM | Hayashi 表 3.x (Bartlett/QS 内核) | 1e-6 | Hansen-Singleton (1982) CCAPM |
| CUE | linearmodels `IVGMMCUE` | 1e-6 | 同上 |
| Arellano-Bond | Stata `xtabond` 文档 | 1e-6 | Arellano-Bond (1991) 论文 abdata, 140 UK firms 1976-1984 |
| Bootstrap CI | Efron-Tibshirani (1993) 法学院数据 | 1e-3 (随机性) | Law school (15 obs) |

### 7.2 跨语言对照矩阵

| v1.5 模块 | Python 对照 | R 对照 | Stata 对照 |
|----------|-------------|--------|-----------|
| OLS + HC0-3 | statsmodels `OLS` `cov_type='HC0'/'HC1'/'HC2'/'HC3'` | `sandwich::vcovHC` `type='HC0'..'HC4'` | `reg, vce(robust)` / `vce(hc3)` |
| Newey-West HAC | statsmodels `cov_type='HAC'` `cov_kwds={'maxlags':L}` | `sandwich::NeweyWest` `prewhite=FALSE` | `newey y x, lag(L)` |
| 聚类 SE | statsmodels `cov_type='cluster'` `cov_kwds={'groups':g}` | `sandwich::vcovCL` `cluster=gg` | `reg, vce(cluster id)` |
| 双向聚类 | statsmodels `cov_type='cluster'` 双 groups (需手动组合) | `sandwich::vcovCL` `cluster=list(g1,g2)` | `cluster2` (user-written, Baum) 或 `vce(cluster ...)` 手动 |
| MLE Logistic | statsmodels `Logit` | `glm(..., family=binomial)` | `logit` |
| MLE Poisson | statsmodels `Poisson` | `glm(..., family=poisson)` | `poisson` |
| QMLE Sandwich | statsmodels GLM `cov_type='HC0'` (GLM 默认 Sandwich) | `sandwich::sandwich` | `glm, vce(robust)` |
| Wald | statsmodels `result.wald_test()` 或 `result.f_test()` | `lmtest::waldtest` / `aod::wald.test` | `test` / `testparm` |
| LR | statsmodels `result.compare_lr_test(restricted)` | `lmtest::lrtest` | `lrtest` |
| LM | statsmodels `statsmodels.stats.diagnostic.het_breuschpagan` (BP) 或 `compare_lm_test` | `lmtest::bptest` (BP) | `estat hettest` |
| 两步 GMM | linearmodels `IVGMM` | `gmm::gmm` | `ivreg2 gmm` |
| CUE | linearmodels `IVGMMCUE` | `gmm::gmm` `type='cue'` | `ivreg2 cue` |
| Arellano-Bond | linearmodels `IVGMM` (动态面板) 或 `AbsorbingLS` | `plm::pgmm` | `xtabond` / `xtabond2` |
| Paired Bootstrap | arch `bootstrap.IIDBootstrap` (或自实现) | `boot::boot` `sim='ordinary'` | `bootstrap` |
| Wild Bootstrap | 自实现 (arch 无 WildBootstrap 类, 通过 `IIDBootstrap` + `extra_kwargs` 传递权重) | `multiwayvcov::cluster.boot` `boot_type='wild'` `wild_type='rademacher'/'mammen'/'norm'` | `boottest` (Roodman) |
| Block Bootstrap | arch `bootstrap.StationaryBootstrap` / `MovingBlockBootstrap` | `boot::tsboot` `sim='fixed'/'geom'` | `bootstrap, cluster` |
| Cluster Bootstrap | 自实现 (arch 无 ClusterBootstrap 类) 或 `multiwayvcov::cluster.boot` `boot_type='xy'` | `multiwayvcov::cluster.boot` (pairs/wild cluster bootstrap) | `boottest, cluster` |

### 7.3 排幻觉点预研 (R 源码对照清单)

| 编号 | 模块 | 幻觉点 | R 源码核对位置 | C++ 实现差异 |
|------|------|--------|--------------|-------------|
| E1 | OLS | R `lm()` 默认含截距, X 矩阵首列为 1 | `lm.fit` 源码 | 显式构造截距列, 用户控制 |
| E2 | HC1 | R `sandwich::vcovHC` HC1 = N/(N-K) · HC0, 非经典 n/(n-k) | `vcovHC` L120-130 | 明确小样本调整公式 |
| E3 | HC2/HC3 | R 实测 leverage h_i = x_i'(X'X)^{-1}x_i (对角元素) | `vcovHC` L150-180 | 计算 projection matrix 对角 |
| E4 | Newey-West | R `sandwich::NeweyWest` 默认调用 `bwNeweyWest()` 自动带宽 (基于 AR(1) 拟合的 Andrews 1991 公式 `L = 1.1447·[α(1)·T]^{1/3}`), **非** Newey-West 1987 论文经验法则 `floor(4·(T/100)^{2/9})`. 用户若需 NW 经验法则须显式 `lag = floor(4*(N/100)^(2/9))` | `NeweyWest.R` `bwNeweyWest` | C++ 默认实现 Andrews 自动带宽, 显式参数可选 NW 经验法则 |
| E5 | HAC Bartlett | R `sandwich::kweights` w[l] = 1 - l/(L+1), 非 1 - l/L | `kweights` L20-40 | 实测 R 源码后标注 |
| E6 | Cluster | R `sandwich::vcovCL` 小样本调整 G/(G-1)·(N-1)/(N-K) | `vcovCL` L80-100 | 明确小样本调整 |
| E7 | MLE Logistic | R `glm` 默认用 IRLS, 非通用 BFGS | `glm.fit` 源码 | IRLS 与 BFGS 数值对比 |
| E8 | QMLE Sandwich | R `sandwich::sandwich` bread = (X'WX)^{-1}, meat = X' diag(ε²) X | `sandwich` L30-60 | 明确 GLM bread 公式 |
| E9 | Wald | R `lmtest::waldtest` 默认 F 检验 (小样本), 非渐近 χ² | `waldtest` 源码 | 提供两种统计量 |
| E10 | GMM Ŝ | R `gmm::gmm` Ŝ 用 tangent matrix, 非 moment matrix HAC | `gmm` L200-250 | Hayashi 教材公式 vs R 实现 |
| E11 | Arellano-Bond | R `plm::pgmm` 工具变量矩阵构造 (GMM-style instruments) | `pgmm` 源码 | 严格按 Arellano-Bond 1991 |
| E12 | Wild Bootstrap | R `multiwayvcov::cluster.boot(boot_type='wild')` 默认 `wild_type='rademacher'` (±1 等概率), Mammen 通过 `wild_type='mammen'` 显式指定 (v=(1±√5)/2, 概率 (5±√5)/10), 另有 `wild_type='norm'` (标准正态). 注意: `multiwayvcov` 包**无** `cluster.wild` 函数, wild bootstrap 通过 `cluster.boot` + `boot_type='wild'` 实现; 默认 Rademacher 遵循 Cameron-Gelbach-Miller 2008 (非 Davidson-Flachaire 2008) | `cluster.boot` 源码 (multiwayvcov 1.2.3) | 三种分布都实现, 默认 Rademacher (CGM 2008) |

### 7.4 测试矩阵设计原则

1. **三层数值验证**:
   - 第一层: 教材表基准 (Greene/Hayashi/MacKinnon-White 表)
   - 第二层: 跨语言对照 (statsmodels/R sandwich/arch, 容差 1e-8)
   - 第三层: 渐近性质验证 (大样本下 Wald ≈ LR ≈ LM)

2. **边界条件测试**:
   - 完美共线性 (X'X 奇异) → SVD fallback 或报错
   - 单聚类 (G=1) → 退化为 robust SE
   - 极小样本 (N < K) → 报错或正则化
   - 空值/缺失值处理策略

3. **随机性测试** (Bootstrap):
   - 固定种子保证可复现
   - 100 次重复的均值方差验证
   - 与理论协方差对比 (蒙特卡洛覆盖概率)

---

## 8. 实施波次与依赖关系

### 8.1 4 波 12 项实施路径

```
M1 (Week 1-2): 第一波 - OLS + HC/HAC/Cluster
  ├─ 1.1 Eigen3 引入 + linalg_dynamic.hpp (ADR-013)
  ├─ 1.2 Estimator 基类 + 数据载体
  ├─ 1.3 OLS + HC0-HC3 (Davidson-MacKinnon)
  ├─ 1.4 Newey-West HAC (Bartlett/QS/Parzen)
  └─ 1.5 聚类 SE (Liang-Zeger + 双向)
       ↓ 依赖 (Estimator 基类)

M2 (Week 3-4): 第二波 - MLE/QMLE + 检验
  ├─ 2.1 MLE 框架 + 6 种族 (Gaussian/Logit/Probit/Poisson/NB/Bernoulli)
  ├─ 2.2 QMLE + Sandwich (White 1982)
  ├─ 2.3 Wald/LR/LM 三大检验
  └─ 2.4 AIC/BIC/HQ 信息准则 + J-test
       ↓ 依赖 (InferenceResult 结构)

M3 (Week 5-6): 第三波 - GMM
  ├─ 3.1 两步 GMM (Hansen 1982)
  ├─ 3.2 迭代 GMM
  ├─ 3.3 CUE (Hansen-Heaton-Yaron 1996)
  └─ 3.4 Arellano-Bond 动态面板
       ↓ 依赖 (HAC 内核复用 M1)

M4 (Week 7-8): 第四波 - Bootstrap
  ├─ 4.1 Bootstrap 基类 + 百分位/BCa CI
  ├─ 4.2 配对 Bootstrap + Wild Bootstrap (Mammen/Rademacher)
  ├─ 4.3 Block Bootstrap (Politis-Romano stationary)
  ├─ 4.4 Cluster Bootstrap
  └─ 4.5 端到端集成测试
```

### 8.2 与已有模块的依赖关系

```
v1.5 依赖 (v1.4.3 已具备):
  ├─ core/types.hpp                 (Real, Size, Index)
  ├─ core/linalg_dynamic.hpp        (ADR-013, v1.5 新增)
  ├─ calibration/optimizer.hpp      (LM/DE 复用于 MLE/GMM, ADR-014)
  └─ core/rng.hpp                   (Philox, 用于 Bootstrap)

v1.5 不依赖:
  ├─ pricing/ (定价模块)
  ├─ models/ (模型模块)
  ├─ risk/ (风险模块)
  └─ hfecon/ (高频计量模块, 独立)

v1.5 为 v1.6+ 预留:
  ├─ Estimator 基类 (半参数/非参数派生)
  ├─ CovarianceType::Custom (用户自定义)
  ├─ MLEFamily::Custom (用户自定义 log-likelihood)
  └─ BootstrapEngine 基类 (Wild Cluster Bootstrap 扩展)
```

---

## 9. 风险与缓解

### 9.1 技术风险

| 风险 | 概率 | 影响 | 缓解措施 |
|------|------|------|---------|
| Eigen3 编译时间膨胀 | 中 | 中 | 定价模块隔离, header-only + 编译单元分离 |
| R 源码幻觉 (公式差异) | 高 | 高 | 严格按 HFE 模块排幻觉流程, E1-E12 预研清单 |
| GMM 数值优化不收敛 | 中 | 高 | 多起始值 + CUE fallback + IRLS 对比 |
| Bootstrap 蒙特卡洛覆盖概率偏差 | 中 | 中 | 固定种子 + 大样本 (B ≥ 999) + BCa 区间 |
| 跨平台数值差异 (MSVC vs GCC) | 低 | 低 | 测试容差 1e-8 (非 1e-12), 容差矩阵分级 |
| Arellano-Bond 工具变量矩阵膨胀 | 中 | 中 | 稀疏矩阵 (Eigen::SparseMatrix) + 分块计算 |

### 9.2 范围风险

| 风险 | 缓解措施 |
|------|---------|
| 范围蔓延 (半参数/非参数内容混入) | 严格遵循"参数方法 4 波 12 项", 半参数推至 v1.6+ |
| 教材公式与开源实现差异 | E1-E12 排幻觉点清单, 每个差异单独标注 |
| Stata 基准不可获取 | 以 R/Python 为主要基准, Stata 仅文档对照 |

### 9.3 性能风险

性能优化**推迟至 Phase 7 独立优化波次**, v1.5 聚焦正确性与可读性:

| 场景 | v1.5 目标 | v2.0+ 优化方向 |
|------|----------|---------------|
| OLS (N=10000, K=100) | < 100ms | BLAS Level 3 + 多线程 |
| GMM 两步 (N=1000, q=50) | < 500ms | 并行 HAC + SVD 分块 |
| Bootstrap (B=999, N=1000) | < 30s | OpenMP 并行 + Philox 分块 (ADR-004) |
| Arellano-Bond (N=1000, T=10) | < 1s | 稀疏矩阵 + 分块 GMM |

---

## 10. 验收检查表 (Phase 6)

| 类别 | 检查项 | 通过标准 |
|------|--------|----------|
| **编译** | 增量编译 < 60s | `cmake --build build --target cpphub` |
| **依赖** | Eigen3 引入不影响定价模块 | 定价测试全绿 (v1.4.3 回归 1412/1412) |
| **单测** | 新增 175+ 测试全绿 | `ctest -R "phase6"` |
| **OLS** | vs Greene 表 3.x 1e-10 | `test_ols_hc.cpp` |
| **HC0-3** | vs MacKinnon-White (1985) 1e-6 | `test_ols_hc.cpp` |
| **HAC** | vs Newey-West (1987) 1e-6 | `test_newey_west.cpp` |
| **Cluster SE** | vs R `sandwich::vcovCL` 1e-8 | `test_cluster_se.cpp` |
| **MLE** | vs Greene 表 17.x 1e-8 | `test_mle_*.cpp` |
| **QMLE** | vs R `sandwich::sandwich` 1e-8 | `test_qmle_sandwich.cpp` |
| **Wald/LR/LM** | vs Greene 表 5.x 1e-6 | `test_hypothesis_tests.cpp` |
| **两步 GMM** | vs Hayashi 表 3.x 1e-6 | `test_gmm_two_step.cpp` |
| **CUE** | vs linearmodels `IVGMMCUE` 1e-6 | `test_gmm_cue.cpp` |
| **Arellano-Bond** | vs Stata `xtabond` 1e-6 | `test_arellano_bond.cpp` |
| **Bootstrap** | vs Efron-Tibshirani (1993) 1e-3 | `test_*_bootstrap.cpp` |
| **跨语言** | statsmodels/R sandwich/arch 1e-8 | CI matrix |
| **跨平台** | MSVC + GCC 一致 | A/B 站验证 |
| **扩展接口** | Estimator 基类支持派生 | 手动 Mock 半参数派生类测试 |
| **Python** | nanobind 绑定 (可选) | `cpphub.OLSEstimator`, `cpphub.GMMEstimator` |

---

## 11. 与 Research OS 的战略对接

### 11.1 "因子失效诊断"方向支撑

| Research OS 组件 | v1.5 提供的基础设施 |
|-----------------|-------------------|
| 因子回归 + Newey-West HAC | M1 第三波 Newey-West HAC |
| 因子组合 t-检验 (多重检验前) | M2 Wald/LR/LM 检验 |
| GMM 资产定价 (SDF 框架) | M3 GMM + Hansen-Jagannathan 距离 (v1.6+) |
| Romano-Wolf 多重检验 | M4 Bootstrap 基类 (复用于 StepM) |
| 因子 FDR (Efron 2010) | M2 MLE 框架 + Empirical Bayes (v1.6+) |

### 11.2 与 Phase 5 HFE 模块的协同

| HFE 输出 | v1.5 消费场景 |
|---------|-------------|
| HAR 模型预测 | HAR 系数检验 (Wald test on β_daily/β_weekly/β_monthly) |
| HEAVY 模型参数 | MLE 标准误差 + QMLE Sandwich 诊断 |
| Realized Measures 序列 | HAR/HEAVY 回归的 Newey-West HAC |

### 11.3 v1.6+ 扩展路线图

| 版本 | 扩展内容 | 依赖 |
|------|---------|------|
| v1.6 | 半参数回归 (Yatchew / Series) | Estimator 基类 (v1.5) |
| v1.6 | 贝叶斯估计 (MCMC) | optimizer.hpp + linalg_dynamic.hpp |
| v1.7 | 状态空间 + Kalman 滤波 | linalg_dynamic.hpp + 时间序列数据载体 |
| v1.8 | 因子诊断专用 (FDR/SPA/MCS/StepM) | Bootstrap 基类 + MLE 框架 |
| v2.0 | 机器学习推断 (DML/Lasso) | Estimator 基类 (MachineLearning 派生) |
| v2.0 | 非参数方法 (Li-Racine 核估计) | Estimator 基类 (Nonparametric 派生) |

---

**Phase 6 负责人**: _______________
**审核人**: _______________
**开始日期**: _______________
**预计结束**: _______________

---

## 附录 A: 教材锚点速查

| v1.5 模块 | 主教材 | 专题教材 | 对照开源库 | 测试基准 |
|----------|--------|---------|-----------|---------|
| OLS + HC/HAC/Cluster | Greene 8ed Ch.3-5 | Davidson-MacKinnon Ch.2-6 (FWL Ch.2, OLS Ch.3, HC Ch.5, HAC Ch.6) | R `sandwich` / statsmodels | MacKinnon-White (1985) 表 1 |
| MLE/QMLE | Greene 8ed Ch.14-17 | Wooldridge CS 2ed Ch.12-13 | statsmodels `discrete`/`glm` | Greene 表 17.x |
| GMM | Hayashi Ch.3-4 | Greene 8ed Ch.13 / Cochrane Ch.10-11 | linearmodels `IVGMM`/`IVGMMCUE` | Hayashi 表 3.x |
| Bootstrap | Cameron-Trivedi Ch.11 | Efron-Tibshirani (1993) / Davison-Hinkley (1997) | arch `bootstrap` / R `boot` | Efron-Tibshirani 法学院数据 |
| 假设检验 | Wooldridge CS 2ed Ch.12-15 | Greene 8ed Ch.5-6 | statsmodels `wald_test`/`compare_lr_test` | Greene 表 5.x |

## 附录 B: ADR 关联

| ADR | 标题 | v1.5 实施关联 |
|-----|------|-------------|
| ADR-002 | Bridge + Virtual Constructor | Estimator::clone() 虚拟构造 |
| ADR-003 | Factory + 静态注册模板 | estimator_factory.hpp |
| ADR-004 | 计数器 RNG | Bootstrap 用 Philox |
| ADR-006 | nanobind Python 绑定 | v1.5 Python 绑定 (可选) |
| ADR-013 | 双层 linalg | linalg_dynamic.hpp (Eigen3) |
| ADR-014 | 标定 vs 估计分离 | econometrics/ vs calibration/ |

## 附录 C: 排幻觉点清单 (E1-E12)

详见 §7.3 测试矩阵, 每个 v1.5 子模块实施前必须完成对应 E 编号的 R 源码核对, 在头文件注释中标注差异 (参考 HFE 模块 D1-D23 风格)。
