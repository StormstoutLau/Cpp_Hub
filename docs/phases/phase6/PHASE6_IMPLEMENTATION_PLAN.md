# Phase 6 实施方案 - 经典参数计量模块 (v1.5)

> **版本**: v1.5 (Phase 6 完整 scope)
> **关联规格**: [PHASE6_ECONOMETRICS_SPEC.md](./PHASE6_ECONOMETRICS_SPEC.md) (审计通过 2026-08-04)
> **关联审计**: [PHASE6_AUDIT_REPORT.md](./PHASE6_AUDIT_REPORT.md) (spec 审计, 24 处幻觉点已修复)
> **关联审计**: [PHASE6_IMPLEMENTATION_AUDIT_REPORT.md](./PHASE6_IMPLEMENTATION_AUDIT_REPORT.md) (实施方案审计, 14 处幻觉点已修复)
> **关联 ADR**: ADR-013 (Accepted, 双层 linalg), ADR-014 (Accepted, calibration/estimation 分离)
> **前置**: v1.4.3 全量回归 1412/1412 通过, Eigen3 引入路径已确认
> **目标**: 交付研究级 C++ 经典参数计量基础设施, 175+ 测试, 跨语言数值对照 1e-8

---

## 1. 实施总览

### 1.1 版本与波次

| 波次 | 周期 | 内容 | 测试数 | 依赖 |
|------|------|------|--------|------|
| **M1** | Week 1-2 | OLS + HC/HAC/Cluster | 60 | Eigen3 + Estimator 基类 |
| **M2** | Week 3-4 | MLE/QMLE + 假设检验 | 55 | M1 InferenceResult |
| **M3** | Week 5-6 | GMM (两步/迭代/CUE) + Arellano-Bond | 35 | M1 HAC 内核复用 |
| **M4** | Week 7-8 | Bootstrap (4 种) + 端到端 | 35 | M1-M3 估计器 |
| **合计** | 8 周 | 4 波 12 项 | **185** | - |

### 1.2 依赖矩阵

```
外部依赖 (v1.5 新增):
  └─ Eigen 3.4.0 (header-only, MPL2.0, ~1.5MB; 注: 5.0.1 已发布但 3.4.0 为 C++14 兼容稳定版)
     ├─ 引入方式: git submodule → third_party/eigen
     ├─ CMake: find_package(Eigen3 CONFIG) 或 INTERFACE include
     └─ 可选后端: MKL/OpenBLAS (v2.0+ 评估, v1.5 不启用)

内部依赖 (v1.4.3 已具备):
  ├─ core/types.hpp           (Real=double, Size, Index)
  ├─ core/rng.hpp             (Philox4x64, Bootstrap 用)
  ├─ calibration/optimizer.hpp (LM/BFGS, MLE/GMM 复用, ADR-014)
  └─ core/linalg.hpp          (固定尺寸, 定价模块用, 不变)

不依赖 (隔离保证):
  ├─ pricing/   (定价模块, header-only 保持)
  ├─ models/    (模型模块)
  ├─ risk/      (风险模块)
  └─ hfecon/    (高频计量, 独立模块)
```

### 1.3 工程约束

| 约束 | 说明 |
|------|------|
| **Header-only** | `include/cpphub/econometrics/` 全部头文件, 无 .cpp (除可能的 Eigen3 显式实例化) |
| **CMake INTERFACE** | `cpphub_econometrics` 导出 Eigen3 依赖, 不污染 `cpphub_core` |
| **C++20** | 与项目一致, 不引入 C++23 特性 |
| **命名空间** | `cpphub::inline namespace v1::econometrics` |
| **异常安全** | 数学错误抛 `std::invalid_argument`/`std::runtime_error`, 不用异常做控制流 |
| **编码规范** | 头文件 `#include` 在 namespace 外 (project_memory 教训), 8 空格缩进 |

---

## 2. 工程基础设施 (M1 第 1-3 天)

### 2.1 Eigen3 引入 (ADR-013 实施)

**策略**: git submodule + CMake INTERFACE, 不用 vcpkg (避免系统依赖)

```bash
# 1. 添加 Eigen3 作为 submodule (锁定 3.4.0 tag)
cd f:/Cpp_Hub
git submodule add https://gitlab.com/libeigen/eigen.git third_party/eigen
cd third_party/eigen
git checkout 3.4.0
```

**CMake 配置** (`CMakeLists.txt` 增量修改):

```cmake
# === Phase 6 v1.5: Eigen3 + 计量模块 (ADR-013) ===
# Eigen3 header-only, 不编译 Eigen 自身
set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen)
add_library(eigen3_interface INTERFACE)
target_include_directories(eigen3_interface INTERFACE ${EIGEN3_INCLUDE_DIR})
# 禁用 Eigen 对 Boost 等的可选依赖
target_compile_definitions(eigen3_interface INTERFACE EIGEN_MPL2_ONLY)

# 计量模块 INTERFACE 库 (ADR-013: 双层 linalg, 与 cpphub_core 隔离)
add_library(cpphub_econometrics INTERFACE)
target_link_libraries(cpphub_econometrics INTERFACE cpphub_core eigen3_interface)
target_include_directories(cpphub_econometrics INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
# MSVC: 抑制 Eigen 内部警告 (C4505 未引用局部函数等)
if(MSVC)
    target_compile_options(cpphub_econometrics INTERFACE /wd4505 /wd4714)
endif()
```

**隔离验证**: `cpphub_core` 不链接 `eigen3_interface`, 定价模块编译不受影响。

### 2.2 目录结构创建

```
include/cpphub/econometrics/         # 新建 (M1 第一天)
├── core/
│   ├── data_types.hpp               # EconData variant (CrossSection/Panel/TimeSeries)
│   ├── estimation_result.hpp        # EstimationResult + InferenceResult + BootstrapResult
│   ├── covariance_type.hpp          # CovarianceType 枚举 + to_string
│   └── estimator_base.hpp           # Estimator 抽象基类 + EstimatorClass 枚举
├── estimation/
│   ├── ols.hpp                      # M1
│   ├── mle.hpp                      # M2
│   ├── gmm.hpp                      # M3
│   └── estimator_factory.hpp        # M2 (ADR-003 风格)
├── inference/
│   ├── standard_errors.hpp          # M1 (HC/HAC/Cluster)
│   ├── hac_kernels.hpp              # M1 (Bartlett/QS/Parzen)
│   ├── hypothesis_tests.hpp         # M2 (Wald/LR/LM)
│   └── diagnostics.hpp              # M2 (AIC/BIC/HQ/J-test)
├── resampling/
│   ├── bootstrap_base.hpp           # M4
│   ├── paired_bootstrap.hpp         # M4
│   ├── wild_bootstrap.hpp           # M4
│   ├── block_bootstrap.hpp          # M4
│   └── cluster_bootstrap.hpp        # M4
└── data/
    ├── panel_data.hpp               # M1 (面板数据载体)
    └── formula.hpp                  # M2 (可选, R 风格 formula)
```

### 2.3 双层 linalg 实施 (ADR-013)

**`include/cpphub/core/linalg_dynamic.hpp`** (新建, 封装 Eigen3):

```cpp
#pragma once
// ADR-013: 双层 linalg 架构 - 动态尺寸层 (封装 Eigen3)
// 用于: econometrics/ 模块 (回归/SVD/QR/稀疏)
// 不用于: pricing/models/risk (保持 core/linalg.hpp 固定尺寸)
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace linalg::dynamic {

class MatrixXD {
public:
    MatrixXD() = default;
    MatrixXD(Size rows, Size cols) : data_(rows, cols) {}
    // ... 构造/访问器/运算符 (委托 Eigen, 不暴露 Eigen 类型)
    Real& operator()(Size i, Size j) { return data_(i, j); }
    Real operator()(Size i, Size j) const { return data_(i, j); }
    Size rows() const { return data_.rows(); }
    Size cols() const { return data_.cols(); }
    Real* data() { return data_.data(); }
    const Real* data() const { return data_.data(); }
    // 显式转换 (边界明确, 不隐式泄漏 Eigen)
    const Eigen::MatrixXd& eigen() const { return data_; }
    Eigen::MatrixXd& eigen() { return data_; }
private:
    Eigen::MatrixXd data_;
};

class VectorXD { /* 类似, 封装 Eigen::VectorXd */ };
class SparseMatrix { /* 封装 Eigen::SparseMatrix<double> */ };

// 计量核心操作 (委托 Eigen + LAPACK)
MatrixXD svd_solve(const MatrixXD& A, const VectorXD& b);   // 最小二乘
MatrixXD inverse_symmetric(const MatrixXD& A);              // (X'X)^{-1} via LLT
MatrixXD cholesky_dynamic(const MatrixXD& A);
VectorXD solve(const MatrixXD& A, const VectorXD& b);       // Ax=b via LU
MatrixXD qr_decompose(const MatrixXD& A);
void svd_full(const MatrixXD& A, MatrixXD& U, VectorXD& S, MatrixXD& V);

}  // namespace linalg::dynamic
}  // namespace v1
}  // namespace cpphub
```

**设计要点**:
- `MatrixXD` 封装 Eigen3, 不暴露 `Eigen::MatrixXd` 到公共 API (除 `eigen()` 显式转换)
- 估计器内部使用 `MatrixXD`, 不直接 `#include <Eigen/Dense>` (降低耦合)
- `inverse_symmetric` 用 `Eigen::LLT` (对称正定), 比 `Eigen::PartialPivLU` 快约 2x (理论复杂度 n³/3 vs 2n³/3, 适用于 SPD 矩阵无需选主元; 注: LLT 适用范围窄于 LU, 非 SPD 矩阵仍需 LU)

---

## 3. M1 详细实施 (Week 1-2: OLS + HC/HAC/Cluster)

### 3.1 Day 1-2: 基础设施 + Estimator 基类

**任务 1.1**: 创建 `core/linalg_dynamic.hpp` (§2.3)
- 实现 `MatrixXD`/`VectorXD`/`SparseMatrix` 封装
- 实现 `svd_solve`/`inverse_symmetric`/`cholesky_dynamic`/`solve`/`qr_decompose`/`svd_full`
- **测试**: `test_linalg_dynamic.cpp` (10 用例, 对照 Eigen 原生 API, 容差 1e-14)

**任务 1.2**: 创建 `core/covariance_type.hpp`
- `CovarianceType` 枚举 (Classical/HC0-5/HAC_Bartlett/HAC_QS/HAC_Parzen/Cluster/Cluster_TwoWay/OPG/Hessian/Sandwich/Bootstrap/Custom)
- `to_string(CovarianceType)` 函数

**任务 1.3**: 创建 `core/data_types.hpp`
- `CrossSectionData`/`PanelData`/`TimeSeriesData` 结构
- `EconData = std::variant<...>` 别名
- `make_cross_section`/`make_panel`/`make_time_series` 便利函数

**任务 1.4**: 创建 `core/estimation_result.hpp`
- `EstimationResult` (coefficients/std_errors/t_statistics/p_values/vcov/log_likelihood/r_squared)
- `InferenceResult` (扩展 EstimationResult, 加 Wald/LR/LM 结果)
- `BootstrapResult` (coef_mean/coef_std/coef_vcov/bootstrap_samples/lower_ci/upper_ci)

**任务 1.5**: 创建 `core/estimator_base.hpp`
- `EstimatorClass` 枚举 (Parametric/Semiparametric/Nonparametric/MachineLearning)
- `Estimator` 抽象基类 (estimate/name/estimatorClass/isParametric/isSemiparametric/isNonparametric/setCovarianceType/clone)
- **关键**: `clone()` 纯虚函数 (ADR-002 风格, 用于工厂注册)

**Day 1-2 测试**: `test_econometrics_core.cpp` (15 用例)
- MatrixXD 构造/访问/运算
- VectorXD 操作
- EconData variant 构造与访问
- EstimationResult 序列化
- Estimator 派生类 Mock (验证 clone/estimate 多态)

### 3.2 Day 3-4: OLS + HC0-HC3

**任务 1.6**: 创建 `estimation/ols.hpp`
- `OLSEstimator : public Estimator`
- `estimate(const EconData&)`: β = (X'X)^{-1} X'y via `linalg::dynamic::inverse_symmetric`
- OLS 专属: `computeRSquared`/`computeAdjustedRSquared`/`computeFStatistic`/`computeFittedValues`/`computeResiduals`/`computeProjectionMatrix`
- **排幻觉点 E1**: 显式构造截距列 (用户控制, 非自动), 对照 R `lm()` 默认含截距

**任务 1.7**: 创建 `inference/standard_errors.hpp` (HC 部分)
- `compute_hc_vcov(X, residuals, XtX_inv, type)`: HC0-HC5
- **排幻觉点 E2**: HC1 = N/(N-K) · HC0 (非 n/(n-k))
- **排幻觉点 E3**: HC2/HC3 leverage h_i = x_i'(X'X)^{-1}x_i (投影矩阵对角)
- HC4: δ_i = min(N·h_i/K, 4) (Cribari-Neto 2004)
- HC5: δ_i = min(N·h_i/K, max(4, 0.7·N·h_max/K)) (Cribari-Neto-Souza 2007, 审计修复)

**Day 3-4 测试**: `test_ols_hc.cpp` (20 用例)
- OLS 系数 vs Greene 表 3.x (Longley 16 obs × 6 regressors, 容差 1e-10)
- HC0-3 标准误差 vs R `sandwich::vcovHC` (Nerlove 145 obs, 容差 1e-8)
- HC4/HC5 vs R `sandwich::vcovHC` (容差 1e-8)
- 边界: 完美共线性 (X'X 奇异) → SVD fallback 或报错
- 边界: N < K → 报错

### 3.3 Day 5-6: HAC 内核 + Newey-West

**任务 1.8**: 创建 `inference/hac_kernels.hpp`
- `HacKernel` 枚举 (Bartlett/QuadraticSpectral/Parzen/TukeyHanning)
- `kernel_weight(kernel, lag_ratio)`: u = l/(L+1)
- `select_max_lag(n_obs, kernel, andrews_optimal, ar1_coef)`
- **排幻觉点 E5**: Bartlett w[l] = 1 - l/(L+1), 非 1 - l/L (R `sandwich::kweights` 实测)
- Parzen: K(u) = 1 - 6u² + 6|u|³ (审计修复, 含 |u|)
- QS: K(u) = 25/(12π²u²) [sin(6πu/5)/(6πu/5) - cos(6πu/5)]

**任务 1.9**: 扩展 `inference/standard_errors.hpp` (HAC 部分)
- `compute_hac_vcov(X, residuals, XtX_inv, type, max_lag, prewhiten)`
- **排幻觉点 E4**: 默认 `max_lag=0` 触发 Andrews 自动带宽 `L = 1.1447·[α(1)·T]^{1/3}` (R `bwNeweyWest` 等价), 非 NW 1987 经验法则 `floor(4·(T/100)^{2/9})`
- Andrews-Monahan (1992) 预白化选项
- 公式: V_HAC = (X'X)^{-1} [Ω_0 + Σ_{l=1}^{L} K(l/L) (Ω_l + Ω_l')] (X'X)^{-1}

**Day 5-6 测试**: `test_hac_kernels.cpp` (15 用例) + `test_newey_west.cpp` (15 用例)
- 内核边界: K(0)=1, K(1)=0 (Bartlett), K(1)=0 (Parzen)
- 内核对称性: K(-u) = K(u)
- Newey-West HAC vs R `sandwich::NeweyWest` (Longley + AR(1) 合成噪声, 容差 1e-8)
- Andrews 自动带宽 vs NW 经验法则 (文档化差异, 两个都测)
- 预白化开关对比

### 3.4 Day 7-8: 聚类 SE + M1 集成

**任务 1.10**: 扩展 `inference/standard_errors.hpp` (Cluster 部分)
- `compute_cluster_vcov(X, residuals, XtX_inv, cluster_id, twoway, cluster_id2)`
- **排幻觉点 E6**: 小样本调整 G/(G-1)·(N-1)/(N-K) (R `sandwich::vcovCL`)
- 双向聚类: V_twoway = V(g1) + V(g2) - V(g1∩g2) (Cameron-Gelbach-Miller 2011)

**任务 1.11**: 创建 `data/panel_data.hpp`
- `PanelData` 长表/宽表转换
- 平衡/非平衡面板检测
- 聚类标识提取

**Day 7-8 测试**: `test_cluster_se.cpp` (10 用例)
- 聚类 SE vs R `sandwich::vcovCL` (Grunfeld 10 firms × 20 years, 容差 1e-8)
- 双向聚类 vs R `sandwich::vcovCL` (cluster=list(g1,g2))
- 边界: G=1 退化为 robust SE
- 边界: 聚类内观测数不等

**M1 集成测试**: `test_integration_m1.cpp` (5 用例)
- 端到端: EconData → OLS → HC/HAC/Cluster → EstimationResult
- 与 HFE HAR 模型对照 (HAR 系数 + Newey-West HAC)

### 3.5 M1 验收检查

| 检查项 | 通过标准 |
|--------|----------|
| Eigen3 编译 | `cmake --build build` 无错误, 增量编译 < 60s |
| 定价模块隔离 | v1.4.3 回归 1412/1412 全绿 (Eigen3 不污染) |
| OLS | vs Greene 表 3.x 1e-10 (Longley) |
| HC0-5 | vs R `sandwich::vcovHC` 1e-8 (Nerlove) |
| HAC | vs R `sandwich::NeweyWest` 1e-8 (Longley + AR(1)) |
| Cluster | vs R `sandwich::vcovCL` 1e-8 (Grunfeld) |
| 测试数 | 60+ 新增测试全绿 |
| 跨平台 | MSVC + GCC 一致 (A/B 站验证) |

---

## 4. M2 详细实施 (Week 3-4: MLE/QMLE + 检验)

### 4.1 Day 9-11: MLE 框架 + 6 种族

**任务 2.1**: 创建 `estimation/mle.hpp`
- `MLEFamily` 枚举 (Gaussian/Logistic/Probit/Poisson/NegativeBinomial/Bernoulli/Exponential/Gamma/Custom)
- `MLEEstimator : public Estimator`
- 复用 `calibration/optimizer.hpp` 的 LM/BFGS (ADR-014)
- 内置 log-likelihood: `gaussianLogLik`/`logisticLogLik`/`poissonLogLik`/`probitLogLik`/`nbLogLik`/`bernoulliLogLik`
- 扩展点: `setCustomLogLikelihood(ll, score, hessian)` (用户自定义)
- **排幻觉点 E7**: R `glm` 用 IRLS, C++ 用 BFGS, 文档化数值差异

**任务 2.2**: 协方差计算 (三种)
- OPG: V = (Σ g_i g_i')^{-1}
- Hessian: V = -H^{-1} (数值 Hessian via 中心差分)
- Sandwich: V = A^{-1} B A^{-1} (White 1982 QMLE)
- **排幻觉点 E8**: GLM bread = (X'WX)^{-1}, meat = X' diag(ε²) X (R `sandwich::sandwich`)

**Day 9-11 测试**: `test_mle_logistic.cpp` (15 用例) + `test_mle_poisson.cpp` (10 用例)
- MLE Logistic vs Greene 表 17.x (Spector-Mazzeo 32 obs, 容差 1e-8)
- MLE Poisson vs statsmodels `Poisson` (DoctorVisits 5190 obs, 容差 1e-8)
- QMLE Sandwich vs R `sandwich::sandwich` (容差 1e-8)
- 收敛性: 多起始值 + IRLS 对比
- 边界: 完全分离 (Logistic) → 报错或正则化

### 4.2 Day 12-13: 假设检验

**任务 2.3**: 创建 `inference/hypothesis_tests.hpp`
- `wald_test(beta, vcov, R, r)`: (Rβ̂-r)'[RVR']^{-1}(Rβ̂-r) ~ χ²(q)
- `lr_test(unrestricted, restricted)`: 2(ℓ_UR - ℓ_R) ~ χ²(q)
- `lm_test(X, residuals_restricted, XtX_inv)`: score test
- `overidentification_test(moments, weighting_matrix, n_params)`: Hansen J-test
- **排幻觉点 E9**: R `lmtest::waldtest` 默认 F 检验 (小样本), C++ 提供两种 (χ² + F)

**任务 2.4**: 创建 `inference/diagnostics.hpp`
- `compute_information_criteria(log_likelihood, n_params, n_obs)`: AIC/BIC/HQ
- AIC = 2K - 2ℓ, BIC = K·log(N) - 2ℓ, HQ = 2K·log(log(N)) - 2ℓ

**Day 12-13 测试**: `test_hypothesis_tests.cpp` (20 用例) + `test_diagnostics.cpp` (10 用例)
- Wald/LR/LM vs Greene 表 5.x (Longley, 容差 1e-6)
- 渐近等价性: 大样本下 Wald ≈ LR ≈ LM
- AIC/BIC/HQ vs statsmodels (容差 1e-10)
- J-test vs Hayashi 表 (容差 1e-6)

### 4.3 Day 14: 工厂 + M2 集成

**任务 2.5**: 创建 `estimation/estimator_factory.hpp`
- `EstimatorFactory` 单例 (ADR-003 风格)
- `EstimatorRegistrar<T>` 模板自动注册
- 注册: OLS/MLE/Gaussian/Logistic/Poisson/...

**M2 集成测试**: `test_integration_m2.cpp` (5 用例)
- 端到端: EconData → MLE → Wald/LR/LM → AIC/BIC
- 因子回归 + Wald 检验 (Research OS 场景)

### 4.4 M2 验收检查

| 检查项 | 通过标准 |
|--------|----------|
| MLE Logistic | vs Greene 表 17.x 1e-8 (Spector-Mazzeo) |
| MLE Poisson | vs statsmodels `Poisson` 1e-8 (DoctorVisits 5190 obs) |
| QMLE Sandwich | vs R `sandwich::sandwich` 1e-8 |
| Wald/LR/LM | vs Greene 表 5.x 1e-6 |
| AIC/BIC/HQ | vs statsmodels 1e-10 |
| 测试数 | 55+ 新增测试全绿 |
| 扩展接口 | MLEFamily::Custom 自定义似然测试 |

---

## 5. M3 详细实施 (Week 5-6: GMM + Arellano-Bond)

### 5.1 Day 15-17: 两步 GMM

**任务 3.1**: 创建 `estimation/gmm.hpp`
- `GMMType` 枚举 (TwoStep/Iterated/CUE)
- `GMMEstimator : public Estimator`
- `MomentFn` 函数类型 (用户指定矩条件)
- `setIVSpecification(instruments, endogenous, exogenous)`: IV 便捷接口
- 两步 GMM 算法 (Hayashi §3.5):
  - Step 1: W₁=I, 估计 θ̂₁ via BFGS
  - Step 2: W₂=Ŝ^{-1}(θ̂₁), 估计 θ̂₂
  - Ŝ = Newey-West HAC (复用 M1 `compute_hac_vcov`)
- **排幻觉点 E10**: R `gmm::gmm` Ŝ 用 tangent matrix, C++ 按 Hayashi 教材用 moment matrix HAC (文档化差异)

**Day 15-17 测试**: `test_gmm_two_step.cpp` (15 用例)
- 两步 GMM vs Hayashi 表 3.x (Hansen-Singleton 1982 CCAPM, Bartlett/QS 内核, 容差 1e-6)
- vs linearmodels `IVGMM` (Python, 容差 1e-8)
- vs R `gmm::gmm` (容差 1e-8)
- J-test 过度识别检验
- 边界: 严格识别 (q=k) → J=0
- 边界: 弱工具变量 → 警告

### 5.2 Day 18-19: 迭代 GMM + CUE

**任务 3.2**: 扩展 `gmm.hpp` (Iterated + CUE)
- 迭代 GMM: 重复 Step 2 直到 ||θ̂_{k+1} - θ̂_k|| < tol
- CUE: θ̂_CUE = argmin_θ ḡ(θ)' Ŝ(θ)^{-1} ḡ(θ) (审计修复符号)
- **关键**: CUE 中 Ŝ 依赖 θ, 数值优化更复杂 (需 trust-region 或 grid search 起始值)
- 收敛诊断: 多起始值 + CUE fallback

**Day 18-19 测试**: `test_gmm_cue.cpp` (10 用例)
- CUE vs linearmodels `IVGMMCUE` (Python, 容差 1e-6)
- CUE vs R `gmm::gmm type='cue'` (容差 1e-6)
- 迭代 GMM 收敛性 (vs 两步 GMM 对比)
- CUE 与两步 GMM 大样本等价性

### 5.3 Day 20-21: Arellano-Bond

**任务 3.3**: 扩展 `gmm.hpp` (Arellano-Bond)
- `estimateArellanoBond(panel, max_lags_instruments, twostep)`
- GMM-style instruments 构造 (Arellano-Bond 1991)
- 工具变量矩阵: 使用 `SparseMatrix` (避免矩阵膨胀)
- **排幻觉点 E11**: R `plm::pgmm` 工具变量矩阵构造 (GMM-style instruments), C++ 严格按 Arellano-Bond 1991

**Day 20-21 测试**: `test_arellano_bond.cpp` (10 用例)
- Arellano-Bond vs Stata `xtabond` 文档 (abdata 140 UK firms 1976-1984, 容差 1e-6)
- vs R `plm::pgmm` (容差 1e-6)
- AR(1)/AR(2) 序列相关检验
- Sargan/Hansen 过度识别检验

### 5.4 M3 验收检查

| 检查项 | 通过标准 |
|--------|----------|
| 两步 GMM | vs Hayashi 表 3.x 1e-6 (Hansen-Singleton CCAPM) |
| CUE | vs linearmodels `IVGMMCUE` 1e-6 |
| Arellano-Bond | vs Stata `xtabond` 1e-6 (abdata) |
| 测试数 | 35+ 新增测试全绿 |
| 稀疏矩阵 | Arellano-Bond 工具变量矩阵用 SparseMatrix |

---

## 6. M4 详细实施 (Week 7-8: Bootstrap + 端到端)

### 6.1 Day 22-24: Bootstrap 基类 + 配对/Wild

**任务 4.1**: 创建 `resampling/bootstrap_base.hpp`
- `BootstrapEngine` 抽象基类
- `resample(estimator, data, n_replicates, seed)` 纯虚
- `percentileCI(samples, alpha)`: 百分位置信区间
- `bcaCI(samples, alpha, estimator, data)`: BCa 区间 (Davison-Hinkley §5.3)
- RNG: 复用 `core/rng.hpp` Philox4x64 (ADR-004, 确定性分块并行)

**任务 4.2**: 创建 `resampling/paired_bootstrap.hpp`
- `PairedBootstrap : public BootstrapEngine`
- 算法: 有放回采样 N 对 → 估计 → 协方差 V_boot = (1/(B-1)) Σ (θ̂*_b - θ̄*)(θ̂*_b - θ̄*)'

**任务 4.3**: 创建 `resampling/wild_bootstrap.hpp`
- `WildBootstrap : public BootstrapEngine`
- `WeightDistribution` 枚举 (Mammen/Rademacher/Webb6)
- 算法: y*_i = x_i'β̂ + v_i · ε̂_i (固定 X, 残差乘权重)
- Rademacher: v = ±1 w.p. 0.5 (默认, Cameron-Gelbach-Miller 2008 推荐; 源: multiwayvcov 文档 "The default is the Rademacher distribution, following CGM (2008)")
- Mammen: v = (1±√5)/2 w.p. (5±√5)/10 (Mammen 1993)
- Webb6: 6 点分布 {-√(3/2), -1, -√(1/2), √(1/2), 1, √(3/2)} 等概率 (Webb 2018, MacKinnon-Webb 2018 Econometrics Journal; 注: R multiwayvcov 1.2.3 未内置 Webb6, 需通过 wild_type=function() 自定义)

**Day 22-24 测试**: `test_paired_bootstrap.cpp` (10 用例) + `test_wild_bootstrap.cpp` (10 用例)
- 配对 Bootstrap CI vs Efron-Tibshirani (1993) 法学院 15 obs (容差 1e-3)
- Wild Bootstrap vs R `multiwayvcov::cluster.boot(boot_type='wild')` (容差 1e-3)
- 三种权重分布对比 (Rademacher/Mammen/Webb6)
- 固定种子可复现性
- 蒙特卡洛覆盖概率 (大样本, 95% CI 覆盖率 0.93-0.97)

### 6.2 Day 25-26: Block + Cluster Bootstrap

**任务 4.4**: 创建 `resampling/block_bootstrap.hpp`
- `BlockBootstrap : public BootstrapEngine`
- `BlockType` 枚举 (Stationary/Circular)
- Politis-Romano (1994) stationary: 块长度 L ~ 几何分布
- `selectBlockLength(series)`: Politis-White (2004) 自动选择
- 经验法则: L = N^{1/3}

**任务 4.5**: 创建 `resampling/cluster_bootstrap.hpp`
- `ClusterBootstrap : public BootstrapEngine`
- 算法: 从 G 个聚类中有放回采样 G 个聚类
- 边界: G < 20 时警告 (推荐 Wild Cluster Bootstrap)

**Day 25-26 测试**: `test_block_bootstrap.cpp` (10 用例) + `test_cluster_bootstrap.cpp` (5 用例)
- Block Bootstrap vs R `boot::tsboot` (容差 1e-3)
- Cluster Bootstrap vs R `multiwayvcov::cluster.boot(boot_type='xy')` (容差 1e-3)
- 块长度自动选择 vs Politis-White
- 时间序列自相关保留验证

### 6.3 Day 27-28: 端到端集成 + M4 验收

**任务 4.6**: `test_integration_phase6.cpp` (10 用例)
- 端到端: OLS → HC/HAC → Wald → Bootstrap CI
- 端到端: MLE → QMLE Sandwich → LR → Bootstrap
- 端到端: GMM → J-test → Bootstrap
- 与 HFE 模块协同: HAR 系数 Wald 检验 + Wild Bootstrap

### 6.4 M4 验收检查

| 检查项 | 通过标准 |
|--------|----------|
| 配对 Bootstrap | vs Efron-Tibshirani 1e-3 (法学院 15 obs) |
| Wild Bootstrap | vs R `multiwayvcov::cluster.boot` 1e-3 |
| Block Bootstrap | vs R `boot::tsboot` 1e-3 |
| Cluster Bootstrap | vs R `multiwayvcov::cluster.boot` 1e-3 |
| 测试数 | 35+ 新增测试全绿 |
| 端到端 | 10 个集成测试通过 |
| 总测试数 | 185+ (v1.4.3 1412 + v1.5 185 = 1597+) |

---

## 7. 兼容性保障

### 7.1 回归测试策略

**每波次结束前必须通过**:
```bash
# 1. v1.4.3 全量回归 (确保不破坏)
cd build && ctest --output-on-failure
# 预期: 1412/1412 通过 (v1.4.3 基线)

# 2. v1.5 新增测试
ctest -R "phase6|econometrics|ols|hc|hac|cluster|mle|gmm|bootstrap"
# 预期: 185+ 通过
```

**隔离保障**:
- `cpphub_core` 不链接 `eigen3_interface` → 定价模块编译零影响
- `econometrics/` 头文件不 `#include "cpphub/pricing/*"` 或 `"cpphub/models/*"`
- 仅依赖: `core/types.hpp`, `core/rng.hpp`, `calibration/optimizer.hpp`, `core/linalg_dynamic.hpp`

### 7.2 跨平台验证 (G4 gate 硬要求)

| 平台 | 编译器 | 验证内容 |
|------|--------|---------|
| 主控站 (Win10) | MSVC 19.3x | 全量 1597+ 测试 |
| A 站 (Ubuntu) | GCC 13.x | fresh clone + rebuild + ctest |
| B 站 (Ubuntu) | GCC 13.x | fresh clone + rebuild + ctest |

**跨平台数值差异**: MSVC vs GCC 在 Eigen3 SVD/QR 上可能有 1e-12 级差异, 测试容差 1e-8 足够。

### 7.3 编译时间控制

| 风险 | 缓解 |
|------|------|
| Eigen3 模板膨胀 | 定价模块隔离, 仅 econometrics/ 编译时包含 Eigen |
| header-only 依赖传播 | `linalg_dynamic.hpp` 仅在 econometrics/ 内部使用 |
| MSVC C4505 警告 | `/wd4505` 抑制 Eigen 内部警告 |
| 增量编译 | Eigen3 头文件不变, 仅 econometrics/ 变更时重编译 |

**目标**: 增量编译 (单文件修改) < 60s, 全量编译 < 5min (主控站)。

---

## 8. 性能考量

### 8.1 v1.5 性能目标 (正确性优先, 性能次要)

| 场景 | v1.5 目标 | 优化方式 |
|------|----------|---------|
| OLS (N=10000, K=100) | < 100ms | Eigen3 BLAS Level 3 (自动) |
| HC0-3 (N=10000, K=100) | < 50ms | 向量化 Σ x_i x_i' ε_i² |
| Newey-West HAC (N=1000, L=10) | < 20ms | 内核权重预计算 |
| MLE Logistic (N=1000, K=10) | < 200ms | BFGS + 数值 Hessian |
| GMM 两步 (N=1000, q=50) | < 500ms | HAC 复用 + SVD 分块 |
| Bootstrap (B=999, N=1000) | < 30s | OpenMP 并行 (Philox 分块) |
| Arellano-Bond (N=1000, T=10) | < 1s | SparseMatrix + 分块 GMM |

### 8.2 Bootstrap 并行化 (ADR-004 应用)

```cpp
// Bootstrap 并行策略 (复用 ADR-004 Philox 计数器 RNG)
// 每个 replicate 独立种子: seed_b = base_seed + b (确定性分块)
// OpenMP 并行 (v1.5 可选, v2.0+ 默认)
#pragma omp parallel for schedule(dynamic) if(B > 100)
for (Size b = 0; b < B; ++b) {
    Philox4x64 rng(base_seed_, b);  // 确定性分块
    auto sample = resampleOnce(estimator, data, rng);
    #pragma omp critical
    results.push_back(sample);
}
```

**关键**: 即使不启用 OpenMP, Philox 分块保证可复现性 (相同种子/相同 B → 位精确相同结果)。

### 8.3 性能优化推迟到 Phase 7

v1.5 聚焦**正确性与可读性**, 不做以下优化:
- BLAS Level 3 手动调优 (Eigen3 自带)
- SIMD 向量化 (Eigen3 自带)
- 内存池 (Bootstrap 样本复用)
- GPU 加速 (Phase 7 评估)

---

## 9. 扩展性设计

### 9.1 Estimator 基类扩展点

```cpp
// v1.6+ 半参数扩展示例 (无需修改 v1.5 代码)
class SemiparametricEstimator : public Estimator {
    bool isSemiparametric() const override { return true; }
    EstimatorClass estimatorClass() const override { return EstimatorClass::Semiparametric; }
    EstimationResult estimate(const EconData& data) override;
    std::unique_ptr<Estimator> clone() const override {
        return std::make_unique<SemiparametricEstimator>(*this);
    }
};
```

### 9.2 扩展点清单

| 扩展场景 | 接口 | 预计版本 |
|---------|------|---------|
| 半参数回归 (Yatchew) | `SemiparametricEstimator : public Estimator` | v1.6+ |
| 核估计 (Li-Racine) | `KernelEstimator : public Estimator` | v2.0+ |
| DML (Chernozhukov 2018) | `DMLEstimator : public Estimator` | v2.0+ |
| Lasso/Post-Double Selection | `LassoEstimator : public Estimator` | v2.0+ |
| 自定义协方差 | `CovarianceType::Custom` + 用户函数 | v1.5 (已预留) |
| 自定义似然 | `MLEFamily::Custom` + `setCustomLogLikelihood` | v1.5 (已预留) |
| Wild Cluster Bootstrap | `WildClusterBootstrap : public BootstrapEngine` | v1.6+ |
| Hansen-Jagannathan 距离 | `HJDistanceEstimator : public Estimator` | v1.6+ |

### 9.3 工厂注册扩展

```cpp
// v1.6+ 新增估计器自动注册 (ADR-003 风格)
static EstimatorRegistrar<SemiparametricEstimator> reg_semi("Semiparametric");
// 无需修改 EstimatorFactory 代码, 链接即注册
```

---

## 10. 风险控制

### 10.1 技术风险

| 风险 | 概率 | 影响 | 缓解 |
|------|------|------|------|
| Eigen3 编译时间膨胀 | 中 | 中 | 定价模块隔离 + INTERFACE 库 |
| R 源码幻觉 (公式差异) | 高 | 高 | E1-E12 排幻觉点预研清单 (审计已修复) |
| GMM 数值优化不收敛 | 中 | 高 | 多起始值 + CUE fallback + IRLS 对比 |
| Bootstrap 覆盖概率偏差 | 中 | 中 | 固定种子 + B ≥ 999 + BCa 区间 |
| Arellano-Bond 工具变量膨胀 | 中 | 中 | SparseMatrix + 分块计算 |
| 跨平台数值差异 | 低 | 低 | 容差 1e-8 (非 1e-12) |

### 10.2 范围风险

| 风险 | 缓解 |
|------|------|
| 范围蔓延 (半参数混入) | 严格遵循"参数方法 4 波 12 项", 半参数推至 v1.6+ |
| 教材公式与开源实现差异 | E1-E12 排幻觉点清单 (每个差异单独标注) |
| Stata 基准不可获取 | 以 R/Python 为主要基准, Stata 仅文档对照 |

### 10.3 实施风险

| 风险 | 缓解 |
|------|------|
| Eigen3 版本不兼容 | 锁定 3.4.0 tag, submodule 固定 |
| MSVC Eigen 警告 | `/wd4505 /wd4714` 抑制 |
| A/B 站 GCC 编译 | Eigen3 纯 header-only, 无平台依赖 |

---

## 11. 交付物清单

### 11.1 代码交付

| 类型 | 数量 | 位置 |
|------|------|------|
| 头文件 | 18 | `include/cpphub/econometrics/` |
| 测试文件 | 16 | `tests/unit/econometrics/` |
| R 基准生成脚本 | 8 | `tests/fixtures/econometrics/` |
| Python 对照脚本 | 4 | `tests/validation/python/` |

### 11.2 文档交付

| 文档 | 状态 |
|------|------|
| [PHASE6_ECONOMETRICS_SPEC.md](./PHASE6_ECONOMETRICS_SPEC.md) | ✅ 审计通过 |
| [PHASE6_AUDIT_REPORT.md](./PHASE6_AUDIT_REPORT.md) | ✅ 24 处幻觉点已修复 |
| [PHASE6_IMPLEMENTATION_PLAN.md](./PHASE6_IMPLEMENTATION_PLAN.md) | ✅ 本文档 |
| ADR-013 (双层 linalg) | ✅ Accepted |
| ADR-014 (calibration/estimation 分离) | ✅ Accepted |

### 11.3 测试交付

| 波次 | 测试数 | 累计 |
|------|--------|------|
| M1 | 60 | 60 |
| M2 | 55 | 115 |
| M3 | 35 | 150 |
| M4 | 35 | 185 |
| **合计** | **185** | (v1.4.3 1412 + v1.5 185 = **1597**) |

---

## 12. 实施时间表

| 周次 | 波次 | 任务 | 交付 |
|------|------|------|------|
| W1 | M1 | Eigen3 + linalg_dynamic + Estimator 基类 + OLS | 25 测试 |
| W2 | M1 | HC/HAC/Cluster + 集成 | 35 测试 |
| W3 | M2 | MLE 框架 + 6 种族 | 25 测试 |
| W4 | M2 | QMLE + 检验 + 诊断 + 工厂 | 30 测试 |
| W5 | M3 | 两步 GMM + 迭代 GMM | 20 测试 |
| W6 | M3 | CUE + Arellano-Bond | 15 测试 |
| W7 | M4 | Bootstrap 基类 + 配对/Wild | 20 测试 |
| W8 | M4 | Block/Cluster + 端到端 + 验收 | 15 测试 |
| **合计** | - | **4 波 12 项** | **185 测试** |

---

## 13. 与 Research OS 的战略对接

### 13.1 因子失效诊断支撑

| Research OS 组件 | v1.5 提供的基础设施 | 使用波次 |
|-----------------|-------------------|---------|
| 因子回归 + Newey-West HAC | M1 HAC (Bartlett/QS/Parzen) | M1 |
| 因子组合 t-检验 | M2 Wald 检验 | M2 |
| GMM 资产定价 (SDF) | M3 GMM + J-test | M3 |
| Romano-Wolf 多重检验 | M4 Bootstrap 基类 (复用于 StepM) | M4 |
| 因子 FDR (Efron 2010) | M2 MLE 框架 + Empirical Bayes (v1.6+) | M2 |

### 13.2 与 HFE 模块协同

| HFE 输出 | v1.5 消费场景 |
|---------|-------------|
| HAR 模型预测 | HAR 系数检验 (Wald on β_daily/β_weekly/β_monthly) |
| HEAVY 模型参数 | MLE 标准误差 + QMLE Sandwich 诊断 |
| 已实现度量序列 | Newey-West HAC (RV 序列自相关) |

---

## 14. 验收检查表 (Phase 6 最终)

| 类别 | 检查项 | 通过标准 |
|------|--------|----------|
| **编译** | 增量编译 < 60s | `cmake --build build --target cpphub` |
| **依赖** | Eigen3 不影响定价模块 | v1.4.3 回归 1412/1412 全绿 |
| **单测** | 185+ 新增测试全绿 | `ctest -R "phase6\|econometrics"` |
| **OLS** | vs Greene 表 3.x 1e-10 | `test_ols_hc.cpp` |
| **HC0-5** | vs MacKinnon-White (1985) 1e-6 | `test_ols_hc.cpp` |
| **HAC** | vs R `sandwich::NeweyWest` 1e-8 | `test_newey_west.cpp` |
| **Cluster** | vs R `sandwich::vcovCL` 1e-8 | `test_cluster_se.cpp` |
| **MLE** | vs Greene 表 17.x 1e-8 | `test_mle_*.cpp` |
| **QMLE** | vs R `sandwich::sandwich` 1e-8 | `test_qmle_sandwich.cpp` |
| **Wald/LR/LM** | vs Greene 表 5.x 1e-6 | `test_hypothesis_tests.cpp` |
| **GMM** | vs Hayashi 表 3.x 1e-6 | `test_gmm_two_step.cpp` |
| **CUE** | vs linearmodels `IVGMMCUE` 1e-6 | `test_gmm_cue.cpp` |
| **Arellano-Bond** | vs Stata `xtabond` 1e-6 | `test_arellano_bond.cpp` |
| **Bootstrap** | vs Efron-Tibshirani 1e-3 | `test_*_bootstrap.cpp` |
| **跨语言** | statsmodels/R sandwich/arch 1e-8 | CI matrix |
| **跨平台** | MSVC + GCC 一致 | A/B 站验证 |
| **扩展接口** | Estimator 基类支持派生 | Mock 半参数派生类测试 |
| **总测试** | 1597+ 全绿 | `ctest` 全量 |

---

## 附录 A: CMake 增量修改清单

```cmake
# === 1. 根 CMakeLists.txt 新增 (Phase 6 v1.5) ===

# Eigen3 (header-only, ADR-013)
set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen)
add_library(eigen3_interface INTERFACE)
target_include_directories(eigen3_interface INTERFACE ${EIGEN3_INCLUDE_DIR})
target_compile_definitions(eigen3_interface INTERFACE EIGEN_MPL2_ONLY)

# 计量模块 INTERFACE 库
add_library(cpphub_econometrics INTERFACE)
target_link_libraries(cpphub_econometrics INTERFACE cpphub_core eigen3_interface)
if(MSVC)
    target_compile_options(cpphub_econometrics INTERFACE /wd4505 /wd4714)
endif()

# === 2. tests/CMakeLists.txt 新增 (Phase 6 v1.5) ===

# Phase 6 v1.5: 经典参数计量模块测试
function(cpphub_add_econ_test name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE cpphub_econometrics GTest::gtest_main)
    gtest_discover_tests(${name})
endfunction()

# M1: OLS + HC/HAC/Cluster
cpphub_add_econ_test(test_linalg_dynamic unit/econometrics/test_linalg_dynamic.cpp)
cpphub_add_econ_test(test_econometrics_core unit/econometrics/test_econometrics_core.cpp)
cpphub_add_econ_test(test_ols_hc unit/econometrics/test_ols_hc.cpp)
cpphub_add_econ_test(test_hac_kernels unit/econometrics/test_hac_kernels.cpp)
cpphub_add_econ_test(test_newey_west unit/econometrics/test_newey_west.cpp)
cpphub_add_econ_test(test_cluster_se unit/econometrics/test_cluster_se.cpp)
cpphub_add_econ_test(test_integration_m1 unit/econometrics/test_integration_m1.cpp)

# M2: MLE/QMLE + 检验
cpphub_add_econ_test(test_mle_logistic unit/econometrics/test_mle_logistic.cpp)
cpphub_add_econ_test(test_mle_poisson unit/econometrics/test_mle_poisson.cpp)
cpphub_add_econ_test(test_qmle_sandwich unit/econometrics/test_qmle_sandwich.cpp)
cpphub_add_econ_test(test_hypothesis_tests unit/econometrics/test_hypothesis_tests.cpp)
cpphub_add_econ_test(test_diagnostics unit/econometrics/test_diagnostics.cpp)
cpphub_add_econ_test(test_integration_m2 unit/econometrics/test_integration_m2.cpp)

# M3: GMM
cpphub_add_econ_test(test_gmm_two_step unit/econometrics/test_gmm_two_step.cpp)
cpphub_add_econ_test(test_gmm_cue unit/econometrics/test_gmm_cue.cpp)
cpphub_add_econ_test(test_arellano_bond unit/econometrics/test_arellano_bond.cpp)

# M4: Bootstrap
cpphub_add_econ_test(test_paired_bootstrap unit/econometrics/test_paired_bootstrap.cpp)
cpphub_add_econ_test(test_wild_bootstrap unit/econometrics/test_wild_bootstrap.cpp)
cpphub_add_econ_test(test_block_bootstrap unit/econometrics/test_block_bootstrap.cpp)
cpphub_add_econ_test(test_cluster_bootstrap unit/econometrics/test_cluster_bootstrap.cpp)
cpphub_add_econ_test(test_integration_phase6 unit/econometrics/test_integration_phase6.cpp)
```

---

## 附录 B: R 基准生成脚本清单

| 脚本 | 用途 | 输出 |
|------|------|------|
| `generate_ols_baselines.R` | OLS + HC0-5 (Longley/Nerlove) | `ols_baselines.json` |
| `generate_hac_baselines.R` | Newey-West HAC (Longley + AR(1)) | `hac_baselines.json` |
| `generate_cluster_baselines.R` | 聚类 SE (Grunfeld) | `cluster_baselines.json` |
| `generate_mle_baselines.R` | MLE Logistic/Poisson (Spector-Mazzeo/DoctorVisits) | `mle_baselines.json` |
| `generate_gmm_baselines.R` | GMM 两步/CUE (Hansen-Singleton) | `gmm_baselines.json` |
| `generate_arellano_bond_baselines.R` | Arellano-Bond (abdata) | `arellano_bond_baselines.json` |
| `generate_bootstrap_baselines.R` | Bootstrap CI (法学院数据) | `bootstrap_baselines.json` |
| `verify_econometrics.R` | R 源码核对 (排幻觉点 E1-E12) | `verify_log.txt` |

**策略**: 沿用 HFE v1.4.0/v1.4.1 硬编码 baseline 策略 (避免 CI 依赖 R 环境), R 脚本仅用于 baseline 生成与版本控制可追溯。

---

## 附录 C: 跨语言对照脚本清单

| 脚本 | Python 库 | 对照模块 |
|------|-----------|---------|
| `cross_validate_ols.py` | statsmodels `OLS` | OLS + HC0-3 |
| `cross_validate_mle.py` | statsmodels `Logit`/`Poisson` | MLE |
| `cross_validate_gmm.py` | linearmodels `IVGMM`/`IVGMMCUE` | GMM + CUE |
| `cross_validate_bootstrap.py` | arch `bootstrap` | Bootstrap |

**容差**: 跨语言 1e-8 (数值稳定场景), 1e-6 (迭代算法), 1e-3 (Bootstrap 随机性)。

---

**实施方案状态**: ✅ 完成, 审计通过 (14 处幻觉点已修复), 可进入 M1 实施
**前置确认**: ADR-013 Accepted, ADR-014 Accepted, spec 审计通过, 实施方案审计通过
**下一步**: 执行 §2.1 Eigen3 引入 + §2.3 linalg_dynamic.hpp 实施
