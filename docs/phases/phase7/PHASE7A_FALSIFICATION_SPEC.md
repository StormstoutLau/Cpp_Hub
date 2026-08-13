# Phase 7A 执行规格书 - 全模块证伪统计量补齐

> **版本归属**: **v1.6** (Phase 7A, 与 Phase 7B 金融时间序列并行)
> **目标**: 系统性补齐 v1.0-v1.5 **全部已实现模块**缺失的证伪统计量, 为 Research OS 因子失效诊断提供事后诊断基础设施
> **覆盖范围**:
>   - Phase 3 (v1.0): 风险管理模块 (`risk/var/`, `risk/greeks/`, `risk/xva.hpp`, `risk/pfe_sa_ccr.hpp`)
>   - Phase 4 (v1.1-v1.2): 定价与 Greeks (`pricing/`, `models/`)
>   - Phase 5 (v1.4): 高频计量模块 (`hfecon/`)
>   - Phase 6 (v1.5): 经典参数计量模块 (`econometrics/`)
> **前置**: Phase 6 (v1.5) 三平台 1767/1767 测试全部通过
> **里程碑**: W1 P0 基础诊断 → W2 P1 模型设定检验 → W3 P2 风险/定价/Greeks 诊断 → W4 集成验收
>
> **Scope 声明**:
> - **严格聚焦事后证伪统计量** (post-estimation falsification statistics), 不含信息论事前度量 (v2.0+ scope, 见 INFORMATION_THEORY_METRICS_RESEARCH.md)
> - **通用统计量** (JB/LB/BG/BP/White/MZ/DM) 放在 `econometrics/inference/` 下, 各模块复用
> - **模块特定诊断** 放在各自模块下 (risk_diagnostics / greeks_consistency / pricing_diagnostics / hfecon_diagnostics)
> - **与 Phase 7B 协同**: 波动率模型诊断 (GARCH/HEAVY 标准化残差) 依赖 Phase 7B GARCH 估计器
> - **教材锚点**: Greene 8ed / Wooldridge CS 2ed / Davidson-MacKinnon 2004 / Tsay 3ed / McNeil-Frey-Embrechts 2005 / Engle-Manganelli 2004
> - **对照库**: statsmodels / R `lmtest` `tseries` `strucchange` `forecast` `rugarch` `PerformanceAnalytics`
>
> **关联**:
> - [ADR-015: 证伪统计量模块边界](../../decisions/ADR_INDEX.md#adr-015-证伪统计量模块边界-通用-vs-模块特定) (Accepted 2026-08-12, 方案 B: 通用诊断不依赖 Eigen3)
> - [ADR015_FALSIFICATION_MODULE_BOUNDARY_RESEARCH.md](../../research/ADR015_FALSIFICATION_MODULE_BOUNDARY_RESEARCH.md) v1.2 (调研报告, 三轮审计排幻觉)
> - [FINANCIAL_TIMESERIES_RESEARCH.md](../../research/FINANCIAL_TIMESERIES_RESEARCH.md) §6-7 (证伪统计量矩阵)
> - [INFORMATION_THEORY_METRICS_RESEARCH.md](../../research/INFORMATION_THEORY_METRICS_RESEARCH.md) (信息论事前度量, v2.0+ scope)
> - [PHASE6_ECONOMETRICS_SPEC.md](../phase6/PHASE6_ECONOMETRICS_SPEC.md) (v1.5 基础设施)

---

## 1. 交付物清单

### 1.1 新增文件结构 (按模块群组织)

```
# 模块群 4: 经典计量 — 通用证伪统计量 (各模块复用, ADR-015 方案 B: 不依赖 Eigen3)
include/cpphub/econometrics/inference/
├── hypothesis_tests.hpp              # 已有 v1.5: Wald/LR/LM
├── diagnostics.hpp                   # 已有 v1.5: AIC/BIC/HQ
├── hac_kernels.hpp                   # 已有 v1.5: Bartlett/QS/Parzen
├── hc_standard_errors.hpp            # 已有 v1.5: HC0-HC5
├── hac_vcov.hpp                      # 已有 v1.5: Newey-West
├── cluster_vcov.hpp                  # 已有 v1.5: Cluster SE
├── residual_diagnostics.hpp          # 新增 P0: JB/LB/BG/BP/White (通用, 仅依赖 core/)
├── volatility_diagnostics.hpp        # 新增 P0: 标准化残差检验 (GARCH/HEAVY 通用, 仅依赖 core/)
├── specification_tests.hpp           # 新增 P1: 信息矩阵/MZ/DM (通用, 仅依赖 core/)
├── structural_break.hpp              # 新增 P2: CUSUM/Andrews (通用, 仅依赖 core/)
├── conduction_metrics.hpp            # 预留 v2.0+: 信息论度量空文件
└── detail/                           # 新增: 通用诊断内部基础设施 (ADR-015 方案 B)
    ├── test_result_base.hpp          # 新增: TestResultBase 通用结果基结构 (组合方式)
    └── ols_simple.hpp                # 新增: 轻量级 OLS (std::vector + Gauss-Jordan, 不依赖 Eigen3)

# 模块群 4: 经典计量 — GMM 弱识别 (需 Eigen3, ADR-015 决策点 2)
include/cpphub/econometrics/estimation/
├── ols.hpp                           # 已有 v1.5 (Eigen3)
├── mle.hpp                           # 已有 v1.5 (Eigen3)
├── gmm.hpp                           # 已有 v1.5 (Eigen3)
└── weak_identification.hpp           # 新增 P1: Cragg-Donald/Stock-Yogo (GMM 专用, 依赖 Eigen3)

# 模块群 1: 风险管理 — 风险模型回测补齐
include/cpphub/risk/var/
├── backtesting.hpp                   # 已有 v1.0: Kupiec/Christoffersen/Traffic Light
└── risk_diagnostics.hpp              # 新增 P0-P2: McNeil-Frey DQ/Berkowitz/MC收敛/ES后验

# 模块群 1: 风险管理 — Greeks 一致性
include/cpphub/risk/greeks/
├── greeks_factory.hpp                # 已有 v1.1
└── greeks_consistency.hpp            # 新增 P2: Analytic vs Numerical vs AAD vs Pathwise vs LR

# 模块群 2: 定价 — 模型拟合优度
include/cpphub/pricing/
├── engine.hpp                        # 已有 v1.2
└── pricing_diagnostics.hpp           # 新增 P2: IV 微笑拟合优度/价格残差诊断

# 模块群 3: 高频计量 — HAR/HEAVY/跳跃诊断
include/cpphub/hfecon/models/
├── har_model.hpp                     # 已有 v1.4.2
├── heavy_model.hpp                   # 已有 v1.4.2
└── hfecon_diagnostics.hpp            # 新增 P0-P1: HAR 残差 LB/MZ (调用通用 volatility_diagnostics)

include/cpphub/hfecon/tests/
├── bns_jump_test.hpp                 # 已有 v1.4
└── jump_test_diagnostics.hpp         # 新增 P2: 跳跃检验多重修正 (Bonferroni/BH)
```

> **ADR-015 方案 B 约束** (Accepted 2026-08-12, 见 [ADR_INDEX.md §ADR-015](../../decisions/ADR_INDEX.md#adr-015-证伪统计量模块边界-通用-vs-模块特定)):
> - 通用诊断头文件 (`inference/` 下) **仅依赖 `core/`**, 不 `#include "cpphub/core/linalg_dynamic.hpp"`, 接口参数不用 `linalg::dynamic::MatrixXD`
> - 回归检验 (BG/BP/White/MZ/CUSUM/Andrews) 用 `detail/ols_simple.hpp` 实现辅助回归 (参考 `har_model.hpp` ols_estimate 的 Gauss-Jordan 实现模式, 决策点 5)
> - `weak_identification.hpp` 移到 `estimation/` (需 Eigen3 特征值分解, Cragg-Donald 浓度矩阵 `G_T`, 决策点 2)
> - `risk/` 和 `pricing/` 维持现有命名空间 `cpphub::v1` (无子命名空间, 决策点 4)
> - TestResultBase 组合方式, 复合诊断 (VolatilityDiagnosticsResult 等) 例外不组合 (决策点 3)
> - 详见 [ADR015 调研报告 v1.2](../../research/ADR015_FALSIFICATION_MODULE_BOUNDARY_RESEARCH.md) (三轮审计排幻觉, 共修正 10 个幻觉点)

### 1.2 新增测试套件

| 测试套件 | 用例数 | 覆盖模块 | 优先级 |
|----------|--------|----------|--------|
| `test_residual_diagnostics` | 25 | JB/LB/BG/BP/White (通用) | P0 |
| `test_volatility_diagnostics` | 15 | 标准化残差 + z² LB (GARCH/HEAVY 通用) | P0 |
| `test_specification_tests` | 20 | 信息矩阵/MZ/DM (通用) | P1 |
| `test_weak_identification` | 12 | Cragg-Donald/Stock-Yogo (GMM) | P1 |
| `test_structural_break` | 15 | CUSUM/Andrews (通用) | P2 |
| `test_risk_diagnostics` | 20 | McNeil-Frey DQ/Berkowitz/MC收敛/ES后验 | P0-P2 |
| `test_greeks_consistency` | 15 | Greeks 跨方法一致性 | P2 |
| `test_pricing_diagnostics` | 15 | IV 拟合优度/价格残差 | P2 |
| `test_hfecon_diagnostics` | 15 | HAR 残差 LB/MZ | P0-P1 |
| `test_jump_test_diagnostics` | 10 | 跳跃检验多重修正 | P2 |
| `test_integration_phase7a` | 10 | 端到端: 全模块诊断流程 | - |

**新增测试总数**: 172

### 1.3 必须达到的数值基准

| 基准 | 容差 | 验收方式 |
|------|------|----------|
| Jarque-Bera vs R `tseries::jarque.bera.test` | 1e-10 | `test_residual_diagnostics.cpp` |
| Ljung-Box vs R `Box.test(type="Ljung")` | 1e-10 | `test_residual_diagnostics.cpp` |
| Breusch-Godfrey vs R `lmtest::bgtest` | 1e-8 | `test_residual_diagnostics.cpp` |
| Breusch-Pagan vs R `lmtest::bptest` | 1e-8 | `test_residual_diagnostics.cpp` |
| White vs R `lmtest::bptest(~fitted+fitted²)` | 1e-8 | `test_residual_diagnostics.cpp` |
| GARCH z² LB vs R `rugarch` residuals + `Box.test` | 1e-8 | `test_volatility_diagnostics.cpp` |
| 信息矩阵 vs statsmodels | 1e-6 | `test_specification_tests.cpp` |
| Mincer-Zarnowitz vs R `lm(actual~fitted)` | 1e-10 | `test_specification_tests.cpp` |
| Diebold-Mariano vs R `forecast::dm.test` | 1e-8 | `test_specification_tests.cpp` |
| Cragg-Donald vs Stata `ivreg2` first-stage | 1e-6 | `test_weak_identification.cpp` |
| CUSUM vs R `strucchange::efp` | 1e-6 | `test_structural_break.cpp` |
| Andrews vs R `strucchange::Fstats` | 1e-6 | `test_structural_break.cpp` |
| McNeil-Frey DQ vs R `vmoinput::DynamicQuantileTest` | 1e-6 | `test_risk_diagnostics.cpp` |
| Berkowitz vs R `VaR::Berkowitz` | 1e-6 | `test_risk_diagnostics.cpp` |
| Greeks 一致性 (Analytic vs AAD) | 1e-6 | `test_greeks_consistency.cpp` |
| IV 拟合优度 χ² | 1e-6 | `test_pricing_diagnostics.cpp` |
| HAR 残差 LB vs R `highfrequency` | 1e-8 | `test_hfecon_diagnostics.cpp` |
| Bonferroni/BH 修正 vs 手算 | 1e-10 | `test_jump_test_diagnostics.cpp` |

---

## 2. 模块群 4: 经典计量 — 通用证伪统计量 (Week 1-2)

> **定位**: 通用统计量工具, 被 risk/hfecon/econometrics 三个模块群复用

### 2.0 `detail/` 公共基础设施 (P0, 前置)

**教材锚点**: R `htest` S3 类 (EnvStats 文档已核实, ADR015 §4.1), `har_model.hpp::ols_estimate` (L186-307 已验证)

**ADR-015 方案 B 决策依据**: 通用诊断头文件仅依赖 `core/`, 不依赖 `linalg_dynamic.hpp` (Eigen3)。回归检验需要 OLS 辅助回归, 但可参考 `har_model.hpp` 用 `std::vector` + Gauss-Jordan 实现, 不必引入 Eigen3。各 `Result` 结构体需统一基字段 (statistic/p_value/method_name/reject_null), 参考 R `htest` 类的设计。**例外**: 复合诊断 (如 `VolatilityDiagnosticsResult`) 聚合多个子检验, 本身不对应单一 htest, 无需 `base` 字段。

**位置**: `include/cpphub/econometrics/inference/detail/`

#### 2.0.1 `test_result_base.hpp` - 通用结果基结构

**设计依据**: R `htest` S3 类有统一字段 (statistic/parameters/p.value/estimate/null.value/alternative/method/data.name, ADR015 §4.1 已核实), 便于跨检验统一打印/序列化/比较。C++ 中采用**组合方式** (非继承) 复用, 避免 vtable 开销, 同时保持 POD 友好。

**接口签名**:

```cpp
#pragma once
#include <string>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {
namespace detail {

// 通用检验结果基结构 (组合方式复用)
// 对应 R htest S3 类的 statistic/p.value/method 三个核心字段
// reject_null 为 Cpp_Hub 扩展字段, 便于业务层快速判断
struct TestResultBase {
    Real statistic;           // 检验统计量 (JB/LB/CD/DM/...)
    Real p_value;             // p 值 (无 p 值时为 NaN, 如 Cragg-Donald 非标准分布)
    std::string method_name;  // 方法名 ("Jarque-Bera"/"Ljung-Box"/...)
    bool reject_null;         // 在 significance_level 下是否拒绝 H0
};

}  // namespace detail
}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
```

**复用示例**:

```cpp
struct JarqueBeraResult {
    detail::TestResultBase base;  // statistic=JB, p_value, method="Jarque-Bera", reject_null
    Real skewness;
    Real kurtosis;
};
```

#### 2.0.2 `ols_simple.hpp` - 轻量级 OLS (std::vector + Gauss-Jordan)

**设计依据**: `har_model.hpp::ols_estimate` (L186-307, 已验证) 用 `std::vector` + Gauss-Jordan 消元实现 OLS, 不依赖 Eigen3。通用诊断的辅助回归 (BG/BP/White/MZ/CUSUM/Andrews) 复用同一 Gauss-Jordan 实现模式, 抽取为公共函数避免每个诊断文件各自重写。与 `ols_estimate` 的有意差异: `ols_simple` 不自动添加常数列 (由调用方决定)、不计算 `adj_r_squared`/`llh` (诊断辅助回归不需要), 仅返回系数 + fitted + residuals + r_squared。

**规模适用范围**: N=百级到千级, K<10 (ADR015 §7.4 H8 修正)。CUSUM/Andrews 用递归/滚动 OLS, N 可能是全样本 (千级), 仍在 std::vector + Gauss-Jordan 可接受范围内 (har_model.hpp 已有先例)。生产级 OLS 估计 (大 N, 性能敏感) 仍用 `econometrics/estimation/ols.hpp` (Eigen3)。

**接口签名**:

```cpp
#pragma once
#include <vector>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {
namespace detail {

// 轻量级 OLS 估计: y = X beta + epsilon (含常数项由调用方决定)
// 实现方式: std::vector + Gauss-Jordan 消元 (partial pivoting)
// 参考: har_model.hpp::ols_estimate (L186-307)
//
// @param y 因变量 (长度 N)
// @param X 设计矩阵 (N×K, 不含常数列; 调用方若需常数项请自行添加一列 1.0)
// @param fitted_values 输出拟合值
// @param residuals 输出残差
// @param r_squared 输出 R²
// @return OLS 系数 (长度 K)
std::vector<Real> ols_simple(
    const std::vector<Real>& y,
    const std::vector<std::vector<Real>>& X,
    std::vector<Real>& fitted_values,
    std::vector<Real>& residuals,
    Real& r_squared);

}  // namespace detail
}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
```

**与 `estimation/ols.hpp` 的关系** (ADR015 §7.4):

| 维度 | `estimation/ols.hpp` | `detail/ols_simple.hpp` |
|------|---------------------|------------------------|
| 依赖 | Eigen3 (`linalg_dynamic`) | `core/` (std::vector) |
| 用途 | 生产级 OLS 估计 (大 N, 性能敏感) | 诊断辅助回归 (N≤千级, K<10) |
| 实现 | Eigen3 SVD/QR | Gauss-Jordan 消元 |
| 重复度 | — | ~50-80 行 (参考 har_model.hpp Gauss-Jordan 实现模式) |

**测试矩阵**: 复用 `har_model.hpp::ols_estimate` 的现有测试 (已验证 N=百级场景), 不单独新增测试套件。

### 2.1 `residual_diagnostics.hpp` - 残差诊断 (P0)

**教材锚点**: Greene 8ed §4.8/§9.5/§13.7, Tsay 3ed §2

**ADR-015 方案 B**: 仅依赖 `core/`, 不依赖 `linalg_dynamic.hpp` (Eigen3)。BG/BP/White 辅助回归用 `detail/ols_simple.hpp`。

**接口签名**:

```cpp
#pragma once
#include <vector>
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/inference/detail/test_result_base.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// Jarque-Bera 正态性检验 (Jarque-Bera 1987)
// JB = N * [S²/6 + (K-3)²/24] ~ χ²(2)
// 排幻觉点 H1: 分母用样本标准差 s (有偏), 非 σ
// 排幻觉点 H2: 公式用峰度 K, K-3 才是超额峰度
struct JarqueBeraResult {
    detail::TestResultBase base;  // statistic=JB, p_value, method="Jarque-Bera", reject_null
    Real skewness;
    Real kurtosis;        // 非超额峰度
};

JarqueBeraResult jarque_bera_test(const std::vector<Real>& residuals);

// Ljung-Box 残差自相关检验 (Ljung-Box 1978)
// LB = N(N+2) Σ ρ_h²/(N-h) ~ χ²(m)
// 排幻觉点 H3: Ljung-Box 用 N(N+2)/(N-h) 加权, 非 Box-Pierce 的 N
// 排幻觉点 H4: lag 自动选择 = min(10, N/5) (Hyndman 推荐)
struct LjungBoxResult {
    detail::TestResultBase base;  // statistic=LB, p_value, method="Ljung-Box", reject_null
    Size lag;
    std::vector<Real> autocorrelations;
};

LjungBoxResult ljung_box_test(const std::vector<Real>& residuals, Size lag = 0);

// Breusch-Godfrey LM 自相关检验 (Breusch 1978, Godfrey 1978)
// 辅助回归 e = Xγ + Σ δ_h e_{t-h} + u, LM = N·R²_aux ~ χ²(p)
// 排幻觉点 H5: 辅助回归必须包含原 X
struct BreuschGodfreyResult {
    detail::TestResultBase base;  // statistic=LM, p_value, method="Breusch-Godfrey", reject_null
    Size lag;
};

// ADR-015 方案 B: X 用 std::vector<std::vector<Real>>, 不依赖 Eigen3
// 辅助回归内部用 detail::ols_simple()
BreuschGodfreyResult breusch_godfrey_test(
    const std::vector<std::vector<Real>>& X,  // N×K 解释变量矩阵 (不含常数列)
    const std::vector<Real>& residuals,
    Size lag);

// Breusch-Pagan 异方差检验 (Breusch-Pagan 1979, Koenker 1981 修正)
// 排幻觉点 H6: 默认 Koenker 修正 (用 e², 非 e²/σ²), R `bptest` 默认也是 studentized
struct BreuschPaganResult {
    detail::TestResultBase base;  // statistic=LM, p_value, method="Breusch-Pagan", reject_null
};

BreuschPaganResult breusch_pagan_test(
    const std::vector<std::vector<Real>>& X,  // ADR-015 方案 B: 不依赖 Eigen3
    const std::vector<Real>& residuals);

// White 异方差检验 (White 1980)
// 辅助回归 e² = Zγ + u, Z = [X, X²交叉项, X²]
// 排幻觉点 H7: 高维 q = K(K+1)/2, N > q 强制检查
struct WhiteResult {
    detail::TestResultBase base;  // statistic=LM, p_value, method="White", reject_null
};

WhiteResult white_test(
    const std::vector<std::vector<Real>>& X,  // ADR-015 方案 B: 不依赖 Eigen3
    const std::vector<Real>& residuals,
    bool include_cross_terms = true);

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
```

**测试矩阵** (25 用例): JB(5) + LB(5) + BG(5) + BP(5) + White(5), 全部 R `lmtest`/`tseries` 对照

### 2.2 `volatility_diagnostics.hpp` - 波动率模型标准化残差检验 (P0)

**教材锚点**: Tsay 3ed Ch.3, McNeil-Frey-Embrechts 2005 §5.3

**适用范围**: GARCH/EGARCH/GJR-GARCH (Phase 7B) + HEAVY (Phase 5)

**ADR-015 方案 B**: 仅依赖 `core/`, 不依赖 `linalg_dynamic.hpp` (Eigen3)。

**接口签名**:

```cpp
#pragma once
#include <vector>
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/inference/residual_diagnostics.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// 波动率模型标准化残差检验
// 原理: 若模型正确设定, z_t = ε_t/√h_t 应 ~ iid N(0,1)
// 排幻觉点 H8: z_t² 的 LB 检验是关键 (非仅 z_t)
//   z_t 无自相关但 z_t² 有自相关 → GARCH 未充分捕捉条件异方差
struct VolatilityDiagnosticsResult {
    std::vector<Real> standardized_residuals;  // z_t = ε_t/√h_t
    LjungBoxResult z_ljung_box;                // z_t 自相关
    LjungBoxResult z_squared_ljung_box;        // z_t² 自相关 (ARCH 效应, 关键)
    JarqueBeraResult z_jarque_bera;            // z_t 正态性
    Real weighted_lb_statistic;                // 加权 LB (Fisher-Gallagher 2012)
    Real weighted_lb_p_value;
    bool model_adequate;                       // 所有检验均不拒绝 → true
};

/// @param residuals 原始残差 ε_t
/// @param conditional_variances 条件方差 h_t
/// @param lag LB 滞后阶数 (0 = 自动)
VolatilityDiagnosticsResult volatility_diagnostics(
    const std::vector<Real>& residuals,
    const std::vector<Real>& conditional_variances,
    Size lag = 0);

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
```

**测试矩阵** (15 用例): 标准化残差计算(3) + z_t LB(3) + z_t² LB(4) + JB(2) + 综合(3)

### 2.3 `specification_tests.hpp` - 模型设定与预测检验 (P1)

**教材锚点**: White 1982, Mincer-Zarnowitz 1969, Diebold-Mariano 1995, Harvey-Leybourne-Newbold 1997

**ADR-015 方案 B**: 仅依赖 `core/`, 不依赖 `linalg_dynamic.hpp` (Eigen3)。MZ 回归用 `detail/ols_simple.hpp`。

**接口签名**:

```cpp
#pragma once
#include <vector>
#include <functional>
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/inference/detail/test_result_base.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// 信息矩阵等式检验 (White 1982)
// IM = N · vec(I_outer - I_inner)' [Â]^{-1} vec(I_outer - I_inner) ~ χ²(q)
// 排幻觉点 H9: 对 QMLE 检验均值方程正确性, 不是方差方程
struct InformationMatrixResult {
    detail::TestResultBase base;  // statistic=IM, p_value, method="Information Matrix", reject_null
    Size df;  // q = K(K+1)/2
};

// ADR-015 方案 B: scores/hessian 用 std::vector<std::vector<Real>>, 不依赖 Eigen3
InformationMatrixResult information_matrix_test(
    const std::vector<std::vector<Real>>& scores,   // N×K 得分矩阵
    const std::vector<std::vector<Real>>& hessian); // K×K Hessian 矩阵

// Mincer-Zarnowitz 回归 (Mincer-Zarnowitz 1969)
// y_t = α + β·ŷ_t + ε_t, H0: α=0 且 β=1
// 排幻觉点 H10: R² 也是预测精度度量
struct MincerZarnowitzResult {
    detail::TestResultBase base;  // statistic=joint F, p_value, method="Mincer-Zarnowitz", reject_null
    Real alpha;
    Real beta;
    Real alpha_t_stat;
    Real beta_t_stat;
    Real r_squared;
};

MincerZarnowitzResult mincer_zarnowitz_regression(
    const std::vector<Real>& actual,
    const std::vector<Real>& forecast);

// Diebold-Mariano 检验 (Diebold-Mariano 1995, HLN 1997 修正)
// DM = d̄ / √(γ̂_0 + 2Σ γ̂_h/N) ~ t(N-1)
// 排幻觉点 H11: HLN 修正用 1/N 计算 γ̂_h, DM ~ t(N-1) 非 N(0,1)
struct DieboldMarianoResult {
    detail::TestResultBase base;  // statistic=DM, p_value, method="Diebold-Mariano", reject_null
    Real mean_loss_diff;
};

DieboldMarianoResult diebold_mariano_test(
    const std::vector<Real>& actual,
    const std::vector<Real>& forecast1,
    const std::vector<Real>& forecast2,
    const std::string& loss_function = "mse",  // "mse"/"mae"/"hmahe"
    Size h = 1);

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
```

**测试矩阵** (20 用例): 信息矩阵(6) + MZ(7) + DM(7)

### 2.4 `weak_identification.hpp` - GMM 弱识别检验 (P1)

**教材锚点**: Cragg-Donald 1993, Stock-Yogo 2005

**ADR-015 决策点 2**: 位于 `econometrics/estimation/` (非 `inference/`), 依赖 Eigen3 (Cragg-Donald 浓度矩阵需特征值分解)。

**接口签名**:

```cpp
#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/econometrics/inference/detail/test_result_base.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// Cragg-Donald 弱识别统计量 (Cragg-Donald 1993)
// CD = min eigenvalue(G_T) 其中 G_T = (X̃'X̃)^{-1/2} X̃'Z̃ (Z̃'Z̃)^{-1} Z̃'X̃ (X̃'X̃)^{-1/2}
//   X̃, Z̃ 为 partialling out 外生控制 W 后的残差矩阵: X̃ = M_W·X, Z̃ = M_W·Z
//   注: X 在外层 (两侧 (X̃'X̃)^{-1/2}), Z 在内层 (见 ADR015 §3.3 核实)
// 排幻觉点 H12: CD 是 F 统计量的矩阵推广, 非 Wald 统计量
// 排幻觉点 H13: Stock-Yogo 2005 临界值表覆盖范围见 §8.3
//   - Table 1 (Bias): K≤3 全覆盖 (L=1..30)
//   - Table 2 (Size): 仅 K≤2 覆盖, K=3 用 Skeels-Windmeijer 2018 解析近似
enum class StockYogoCriterion {
    RelativeBias,    // SY2005 TSLS Bias 表, K≤3 全覆盖 (L≤30)
    SizeDistortion   // SY2005 TSLS Size 表, K≤2 原表, K=3 用 Skeels-Windmeijer 2018 近似
};

struct WeakIdentificationResult {
    detail::TestResultBase base;       // statistic=CD, p_value=N/A (非标准分布), method="Cragg-Donald", reject_null
    Size n_instruments;               // L
    Size n_endogenous;                // K
    StockYogoCriterion criterion;
    Real bias_threshold;              // 0.05/0.10/0.20/0.30 (RelativeBias 时有效)
    Real size_threshold;              // 0.10/0.15/0.20/0.25 (SizeDistortion 时有效)
    Real stock_yogo_critical_value;
    bool critical_value_is_exact;     // true=原表查表, false=解析近似
    std::string critical_value_source;// "SY2005 Table 1/2" 或 "Skeels-Windmeijer 2018 approx"
    std::string interpretation;
};

WeakIdentificationResult cragg_donald_test(
    const linalg::dynamic::MatrixXD& Z,  // 工具变量矩阵 (N × L)
    const linalg::dynamic::MatrixXD& X_endog,  // 内生变量矩阵 (N × K)
    const linalg::dynamic::MatrixXD& X_exog,   // 外生变量矩阵 (N × K_exog)
    StockYogoCriterion criterion = StockYogoCriterion::RelativeBias,
    Real threshold = 0.10);  // bias: 0.05/0.10/0.20/0.30, size: 0.10/0.15/0.20/0.25

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
```

**测试矩阵** (12 用例): CD 统计量(4) + Stock-Yogo Table 1 查表 (K=1,2,3 × L=1..5, 6 用例, 硬编码原表值) + Skeels-Windmeijer 2018 近似 (K=3, Size 准则, 2 用例, 容差 1e-4)

### 2.5 `structural_break.hpp` - 结构断点检验 (P2)

**教材锚点**: Brown-Durbin-Evans 1975, Andrews 1993, Hansen 1997

**ADR-015 方案 B**: 仅依赖 `core/`, 不依赖 `linalg_dynamic.hpp` (Eigen3)。递归/滚动 OLS 用 `detail/ols_simple.hpp`。

**接口签名**:

```cpp
#pragma once
#include <vector>
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/inference/detail/test_result_base.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// CUSUM 结构稳定性检验 (Brown-Durbin-Evans 1975)
// 排幻觉点 H14: 用递归残差, 非普通残差 (Greene 8ed §9.6)
struct CusumResult {
    detail::TestResultBase base;  // statistic=max|CUSUM|, p_value, method="CUSUM", reject_null
    std::vector<Real> cusum_path;
    std::vector<Real> confidence_band;
    Size breakpoint_estimate;
};

// ADR-015 方案 B: X 用 std::vector<std::vector<Real>>, 不依赖 Eigen3
CusumResult cusum_test(
    const std::vector<std::vector<Real>>& X,
    const std::vector<Real>& y,
    Real significance_level = 0.05);

// Andrews 未知断点检验 (Andrews 1993)
// supLR = max_π LR(π), Hansen 1997 p 值
// 排幻觉点 H15: p 值用 Hansen 1997 非标准分布, 非 χ²/F
struct AndrewsBreakpointResult {
    detail::TestResultBase base;  // statistic=supLR, p_value, method="Andrews", reject_null
    Size breakpoint_estimate;
    Real breakpoint_fraction;
};

AndrewsBreakpointResult andrews_breakpoint_test(
    const std::vector<std::vector<Real>>& X,  // ADR-015 方案 B: 不依赖 Eigen3
    const std::vector<Real>& y,
    Real trim_fraction = 0.15);

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
```

**测试矩阵** (15 用例): CUSUM(7) + Andrews(8)

---

## 3. 模块群 1: 风险管理 — 风险模型回测补齐 (Week 1, 3)

### 3.1 `risk_diagnostics.hpp` - VaR/ES 后验诊断 (P0-P2)

**教材锚点**: Engle-Manganelli 2004 (DQ), Berkowitz 2001, McNeil-Frey-Embrechts 2005

**依赖**: 已有 `risk/var/backtesting.hpp` (Kupiec/Christoffersen)

**ADR-015 决策点 4**: 命名空间维持 `cpphub::v1` (无 `risk` 子命名空间, 与现有 backtesting.hpp 一致)。

**接口签名**:

```cpp
#pragma once
#include <vector>
#include <functional>
#include "cpphub/core/types.hpp"
#include "cpphub/risk/var/backtesting.hpp"

namespace cpphub {
inline namespace v1 {

// McNeil-Frey 动态量化检验 (Engle-Manganelli 2004)
// 辅助回归 Hit_t = X_t'γ + u_t, DQ = (Hit-π)'X(X'X)^{-1}X'(Hit-π)/(π(1-π)) ~ χ²(q)
// 排幻觉点 H16: X_t 包含 [1, Hit_{t-1}, VaR_t, ...], 非仅 Hit_{t-1}
struct DynamicQuantileResult {
    Real dq_statistic;
    Size df;
    Real p_value;
    bool reject_correct_coverage;
};

DynamicQuantileResult dynamic_quantile_test(
    const std::vector<Real>& returns,
    const std::vector<Real>& var_forecasts,
    Real confidence_level,
    const std::vector<Size>& hit_lags = {1, 2, 3, 4});

// Berkowitz 尾部检验 (Berkowitz 2001)
// z = Φ⁻¹(F_model(r)), 检验 z ~ N(0,1) 且无自相关
// 排幻觉点 H17: 需模型 CDF, 非经验 CDF (检验的是"模型正确性"非"分布形状")
struct BerkowitzResult {
    Real lr_statistic;    // ~ χ²(3)
    Real lr_mean;
    Real lr_variance;
    Real lr_autocorr;
    Real p_value;
    bool reject_correct_distribution;
};

BerkowitzResult berkowitz_test(
    const std::vector<Real>& returns,
    std::function<Real(Real)> model_cdf,
    Size lag = 1);

// MC 收敛性诊断 (MCVar 专用)
// 原理: MC 标准误差应随 √N 衰减, 偏离则提示路径数不足或方差过大
// 排幻觉点 H18: 标准误差用批次均值法 (batch means), 非单次估计
struct MCConvergenceResult {
    Real estimated_std_error;    // 估计的 MC 标准误差
    Size n_paths;
    Size n_batches;              // 批次数
    std::vector<Real> batch_means;
    Real convergence_rate;       // 应接近 1/√2 ≈ 0.707 (路径数翻倍时 SE 减半)
    bool converged;              // SE < tolerance
};

MCConvergenceResult mc_convergence_diagnosis(
    const std::vector<Real>& mc_estimates,  // 每条路径的估计值
    Size n_batches = 20,
    Real tolerance = 1e-4);

// ES 后验检验 (Expected Shortfall backtesting)
// 原理: 超越条件均值应等于 ES 预测
// 排幻觉点 H19: ES 后验需在 VaR 超越条件下, 非全样本
//   ES_realized = mean(r_t | r_t < -VaR_t), 应接近 ES_forecast
struct ESBacktestResult {
    Real es_forecast_mean;
    Real es_realized_mean;      // 超越条件下的实际均值
    Real bias;                  // es_realized - es_forecast
    Size n_violations;
    Real t_stat;                // (es_realized - es_forecast) / SE
    Real p_value;
    bool reject_correct_es;
};

ESBacktestResult es_backtest(
    const std::vector<Real>& returns,
    const std::vector<Real>& var_forecasts,
    const std::vector<Real>& es_forecasts);

}  // namespace v1
}  // namespace cpphub
```

**测试矩阵** (20 用例): DQ(8) + Berkowitz(7) + MC收敛(3) + ES后验(2)

### 3.2 `greeks_consistency.hpp` - Greeks 跨方法一致性检验 (P2)

**教材锚点**: Glasserman 2003 §7 (Monte Carlo Greeks), Broadie-Glasserman 1996

**依赖**: 已有 `risk/greeks/greeks_factory.hpp` (Analytic/Numerical/AAD/Pathwise/LR)

**接口签名**:

```cpp
#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/risk/greeks/greeks_factory.hpp"

namespace cpphub {
inline namespace v1 {

// Greeks 跨方法一致性检验
// 原理: 同一 Greeks 用不同方法计算应数值一致 (容差内)
//   - Analytic vs Numerical: FD 精度受 bump size 影响, 应 < 1%
//   - Analytic vs AAD: 应精确一致 (AAD 是解析导数)
//   - Pathwise vs LR: 方差不同但均值应一致
//   - Numerical vs AAD: 应 < 0.1%
//
// 排幻觉点 H20: Pathwise/LR 有随机性, 用置信区间比较而非点估计
//   |mean_pathwise - analytic| < z_{0.975} * SE_pathwise
// 排幻觉点 H21: Gamma 的 Numerical 用二阶差分, bump size 需更小 (1e-4 而非 1e-2)
struct GreeksConsistencyResult {
    Real analytic_value;
    Real numerical_value;
    Real aad_value;
    Real pathwise_mean;
    Real pathwise_std_error;
    Real lr_mean;
    Real lr_std_error;
    Real max_discrepancy;        // 最大相对偏差
    bool consistent;             // 所有方法在容差内一致
    std::vector<std::string> warnings;  // 不一致的方法对
};

GreeksConsistencyResult greeks_consistency_check(
    Real S, Real K, Real T, Real r, Real q, Real sigma,
    PayoffType payoff,
    const std::string& greek_name,  // "delta"/"gamma"/"vega"/"theta"/"rho"
    Size n_paths = 100000,
    uint64_t seed = 42,
    Real relative_tolerance = 0.01);

}  // namespace v1
}  // namespace cpphub
```

**测试矩阵** (15 用例): Delta(3) + Gamma(3) + Vega(3) + Theta(3) + Rho(3)

---

## 4. 模块群 2: 定价 — 模型拟合优度诊断 (Week 3)

### 4.1 `pricing_diagnostics.hpp` - IV 微笑拟合与价格残差 (P2)

**教材锚点**: Gatheral 2006 (SVI), Fengler 2009 (IV surface)

**接口签名**:

```cpp
#pragma once
#include <vector>
#include "cpphub/core/types.hpp"
#include "cpphub/core/linalg_dynamic.hpp"

namespace cpphub {
inline namespace v1 {

// 隐含波动率微笑拟合优度检验
// 原理: 模型 IV vs 市场 IV 的残差应无系统性偏差
//   χ² = Σ (IV_model - IV_market)² / σ_market² ~ χ²(N-1)
// 排幻觉点 H22: 权重用市场 IV 的 Bid-Ask 宽度, 非简单方差
struct IVFitGoodnessResult {
    Real chi_squared;
    Size degrees_of_freedom;
    Real p_value;
    Real rmse;                    // 均方根误差
    Real max_abs_residual;        // 最大绝对残差
    std::vector<Real> residuals;
    bool reject_good_fit;
};

IVFitGoodnessResult iv_fit_goodness_test(
    const std::vector<Real>& strikes,
    const std::vector<Real>& maturities,
    const std::vector<Real>& iv_market,
    const std::vector<Real>& iv_model,
    const std::vector<Real>& iv_bid_ask_spread);  // 权重 = 1/spread

// 模型价格 vs 市场价格残差诊断
// 原理: 检验残差无自相关/无偏/同方差
struct PriceResidualDiagnostics {
    Real mean_residual;
    Real std_residual;
    Real t_stat_bias;             // H0: mean = 0
    Real p_value_bias;
    bool has_bias;
    std::vector<Real> residuals;
};

PriceResidualDiagnostics price_residual_analysis(
    const std::vector<Real>& market_prices,
    const std::vector<Real>& model_prices);

}  // namespace v1
}  // namespace cpphub
```

**测试矩阵** (15 用例): IV 拟合优度(8) + 价格残差(7)

---

## 5. 模块群 3: 高频计量 — HAR/HEAVY/跳跃诊断 (Week 1-2)

### 5.1 `hfecon_diagnostics.hpp` - HAR/HEAVY 预测诊断 (P0-P1)

**教材锚点**: Corsi 2009, Shephard-Sheppard 2010, Patton 2011 (RV 预测评估)

**依赖**: 已有 `hfecon/models/har_model.hpp`, `hfecon/models/heavy_model.hpp`

**接口签名**:

```cpp
#pragma once
#include <vector>
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/inference/residual_diagnostics.hpp"
#include "cpphub/econometrics/inference/specification_tests.hpp"
#include "cpphub/econometrics/inference/volatility_diagnostics.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// HAR 模型诊断
// 原理: HAR 残差应无自相关 (LB), 预测应无偏 (MZ)
struct HARDiagnosticsResult {
    // 残差诊断
    econometrics::LjungBoxResult residual_ljung_box;  // 残差自相关
    econometrics::JarqueBeraResult residual_jb;        // 残差正态性 (HAR 假设)
    
    // 预测精度
    econometrics::MincerZarnowitzResult mz_regression; // 预测无偏性
    
    // 拟合优度
    Real r_squared;
    Real adjusted_r_squared;
    
    bool model_adequate;
};

HARDiagnosticsResult har_diagnostics(
    const std::vector<Real>& actual_rv,      // 实际 RV
    const std::vector<Real>& fitted_rv,       // HAR 拟合值
    const std::vector<Real>& residuals);      // HAR 残差

// HEAVY 模型诊断 (复用通用 volatility_diagnostics)
struct HEAVYDiagnosticsResult {
    econometrics::VolatilityDiagnosticsResult variance_equation;  // h_t 方程
    econometrics::VolatilityDiagnosticsResult measurement_equation; // RM 方程
    
    // HEAVY 特有: 两方程交叉诊断
    Real correlation_h_rm;  // h_t 与 RM_t 的相关性 (应接近 1)
    bool model_adequate;
};

HEAVYDiagnosticsResult heavy_diagnostics(
    const std::vector<Real>& rm_residuals,           // RM 方程残差
    const std::vector<Real>& rm_conditional_means,   // RM 条件均值
    const std::vector<Real>& variance_residuals,     // 方差方程残差
    const std::vector<Real>& conditional_variances); // h_t

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
```

**测试矩阵** (15 用例): HAR 残差 LB(4) + HAR MZ(4) + HEAVY 标准化残差(4) + HEAVY 交叉诊断(3)

### 5.2 `jump_test_diagnostics.hpp` - 跳跃检验多重修正 (P2)

**教材锚点**: Bonferroni 1936, Benjamini-Hochberg 1995

**接口签名**:

```cpp
#pragma once
#include <vector>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// 多重检验修正 (Bonferroni / Benjamini-Hochberg)
// 原理: 多个跳跃检验 (BNS/AJ/JO/Rank) 联合推断时需控制 FWER 或 FDR
// 排幻觉点 H23: BH 修正控制 FDR (False Discovery Rate), 非 FWER
//   Bonferroni: α' = α/m (保守, FWER ≤ α)
//   BH: 排序后 p_(i) ≤ (i/m)·α (FDR ≤ α)
struct MultipleTestCorrectionResult {
    enum class Method { Bonferroni, BenjaminiHochberg };
    
    std::vector<Real> original_p_values;
    std::vector<Real> adjusted_p_values;
    std::vector<bool> reject_null;           // 修正后是否拒绝
    Method method;
    Size n_tests;
    Size n_rejections;
};

MultipleTestCorrectionResult multiple_test_correction(
    const std::vector<Real>& p_values,
    MultipleTestCorrectionResult::Method method,
    Real alpha = 0.05);

// 跳跃检验联合诊断 (BNS + AJ + JO + Rank 四检验联合推断)
struct JumpTestDiagnosticsResult {
    std::vector<Real> test_statistics;  // [BNS, AJ, JO, Rank]
    std::vector<Real> p_values;
    MultipleTestCorrectionResult bonferroni;
    MultipleTestCorrectionResult bh;
    Size consensus_jumps;  // 多数投票的跳跃天数
    bool consistent;       // 四检验结论一致
};

JumpTestDiagnosticsResult jump_test_diagnostics(
    const std::vector<Real>& bns_stats,
    const std::vector<Real>& aj_stats,
    const std::vector<Real>& jo_stats,
    const std::vector<Real>& rank_stats,
    Real alpha = 0.05);

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
```

**测试矩阵** (10 用例): Bonferroni(3) + BH(3) + 跳跃检验联合(4)

---

## 6. 集成与验收 (Week 4)

### 6.1 集成测试 (`test_integration_phase7a.cpp`)

10 个端到端测试, 覆盖**全模块**诊断流程:

| 用例 | 模块群 | 流程 | 验证点 |
|------|--------|------|--------|
| 1 | 计量 | OLS → BP/White → HC 选择 | 异方差检测触发 HC 切换 |
| 2 | 计量 | OLS → LB/BG → HAC 选择 | 自相关检测触发 HAC 切换 |
| 3 | 高频 | HAR → 残差 LB → MZ | HAR 完整诊断流程 |
| 4 | 高频 | HEAVY → 标准化残差 → z² LB | HEAVY 完整诊断流程 |
| 5 | 风险 | VaR → Kupiec + DQ + Berkowitz + ES | 风险模型完整回测 |
| 6 | 风险 | Greeks → Analytic vs AAD vs Pathwise | Greeks 跨方法一致性 |
| 7 | 定价 | Heston → IV 拟合优度 → 价格残差 | 定价模型诊断 |
| 8 | 计量 | MLE → 信息矩阵检验 | 模型设定诊断 |
| 9 | 计量 | OLS → CUSUM → Andrews | 结构稳定性诊断 |
| 10 | 高频 | 跳跃检验 → Bonferroni/BH 修正 | 多重检验修正 |

### 6.2 验收标准

| 标准 | 要求 |
|------|------|
| 测试通过率 | 172/172 (三平台: MSVC + GCC A 站 + GCC B 站) |
| 数值基准 | 所有 18 项基准容差达标 |
| 跨平台一致性 | Philox4x64 RNG (复用 v1.5 M4 经验) |
| 无回归 | v1.5 1767 测试全部通过 (无退化) |
| 排幻觉 | 23 个排幻觉点全部验证 |
| 代码风格 | header-only, namespace 外 #include, Doxygen 注释 |

### 6.3 交付物

| 交付物 | 路径 |
|--------|------|
| 源码 | 11 个新增头文件 (5 通用 + 6 模块特定) |
| 测试 | 10 个新增测试套件 (172 用例) |
| Spec | `docs/phases/phase7/PHASE7A_FALSIFICATION_SPEC.md` (本文件) |
| 验收报告 | `docs/phases/phase7/PHASE7A_ACCEPTANCE.md` (验收后生成) |

---

## 7. 模块群覆盖矩阵

### 7.1 各模块群补齐状态

| 模块群 | 已有证伪统计量 (v1.5) | 新增证伪统计量 (v1.6) | 补齐后状态 |
|--------|---------------------|---------------------|-----------|
| **风险管理** (`risk/var/`) | Kupiec POF, Christoffersen CC/IND, Traffic Light | McNeil-Frey DQ, Berkowitz, MC收敛, ES后验 | **完整覆盖** |
| **Greeks** (`risk/greeks/`) | - | 跨方法一致性检验 | **完整覆盖** |
| **定价** (`pricing/`, `models/`) | - | IV 拟合优度, 价格残差诊断 | **基础覆盖** |
| **高频计量** (`hfecon/`) | harInsanityFilter, BNS/AJ/JO/Rank (本身是证伪统计量) | HAR 残差 LB/MZ, HEAVY 标准化残差, 跳跃检验多重修正 | **完整覆盖** |
| **经典计量** (`econometrics/`) | Wald/LR/LM, Hansen J, AR(1)/AR(2), BCa, AIC/BIC/HQ | JB/LB/BG/BP/White, 信息矩阵, MZ/DM, Cragg-Donald, CUSUM/Andrews | **完整覆盖** |

### 7.2 优先级分布

| 优先级 | 模块群 | 统计量 | 代码量估计 |
|--------|--------|--------|-----------|
| **P0** | 计量 | JB/LB/BG/BP/White | ~600 行 |
| **P0** | 计量+高频 | 标准化残差检验 (GARCH/HEAVY) | ~150 行 |
| **P0** | 风险 | MC 收敛诊断 | ~100 行 |
| **P0** | 高频 | HAR 残差 LB | ~150 行 (复用通用) |
| **P1** | 计量 | 信息矩阵/MZ/DM | ~400 行 |
| **P1** | 计量 | Cragg-Donald 弱识别 | ~250 行 |
| **P1** | 高频 | HAR MZ/HEAVY 诊断 | ~200 行 |
| **P2** | 风险 | McNeil-Frey DQ/Berkowitz/ES | ~500 行 |
| **P2** | 计量 | CUSUM/Andrews | ~550 行 |
| **P2** | Greeks | 跨方法一致性 | ~300 行 |
| **P2** | 定价 | IV 拟合优度/价格残差 | ~200 行 |
| **P2** | 高频 | 跳跃检验多重修正 | ~150 行 |
| **总计** | - | - | **~3550 行** |

---

## 8. 风险与依赖

### 8.1 依赖关系

```
Phase 7A (证伪统计量补齐)
├── 公共基础设施 (econometrics/inference/detail/) ← ADR-015 方案 B 前置
│   ├── test_result_base.hpp         ← core/ (无 Eigen3)
│   └── ols_simple.hpp               ← core/ (无 Eigen3, 参考 har_model.hpp)
│
├── 通用模块 (econometrics/inference/) ← Phase 6 (v1.5 基础设施) + detail/
│   ├── residual_diagnostics.hpp     ← detail/ (ols_simple + test_result_base)
│   ├── volatility_diagnostics.hpp   ← residual_diagnostics.hpp + Phase 7B GARCH (可选, 可用合成数据)
│   ├── specification_tests.hpp      ← detail/ (ols_simple + test_result_base)
│   ├── structural_break.hpp         ← detail/ (ols_simple + test_result_base)
│   └── conduction_metrics.hpp       ← v2.0+ 预留 (空文件)
│
├── GMM 弱识别 (econometrics/estimation/) ← ADR-015 决策点 2 (需 Eigen3)
│   └── weak_identification.hpp      ← core/linalg_dynamic.hpp (Eigen3) + detail/test_result_base
│
├── 风险模块 (risk/) ← Phase 3 (v1.0)
│   ├── risk_diagnostics.hpp         ← backtesting.hpp (已有)
│   └── greeks_consistency.hpp       ← greeks_factory.hpp (已有)
│
├── 定价模块 (pricing/) ← Phase 4 (v1.1-v1.2)
│   └── pricing_diagnostics.hpp      ← 无新依赖
│
└── 高频模块 (hfecon/) ← Phase 5 (v1.4)
    ├── hfecon_diagnostics.hpp       ← har_model.hpp, heavy_model.hpp (已有)
    │                                 + econometrics/inference/ (residual_diagnostics,
    │                                   specification_tests, volatility_diagnostics, 见 §5.1)
    │                                 [注: 首次引入 hfecon → econometrics 依赖, 见 §8.2]
    └── jump_test_diagnostics.hpp    ← bns/aj/jo/rank (已有)
```

**ADR-015 方案 B 依赖约束** (Accepted 2026-08-12, 决策点 1):
- `inference/detail/` 和 `inference/` 下通用诊断仅依赖 `core/`, 不依赖 Eigen3
- `estimation/weak_identification.hpp` 依赖 Eigen3 (Cragg-Donald 浓度矩阵需特征值分解, 决策点 2), 位于 `cpphub_econometrics` INTERFACE 库内
- `risk/` 和 `pricing/` 调用通用诊断时不引入 Eigen3 (因通用诊断不依赖 Eigen3)
- `hfecon/` 调用 `inference/` 下通用诊断 (hfecon_diagnostics), 首次引入 hfecon → econometrics/inference 依赖, 但因通用诊断无 Eigen3, 不违反 ADR-013

### 8.2 已知风险

| 风险 | 概率 | 影响 | 缓解 |
|------|------|------|------|
| Phase 7B GARCH 估计器延期 | 中 | volatility_diagnostics GARCH 测试 | 用合成数据 + HEAVY (已有) 测试 |
| R 基准环境差异 | 低 | 数值对照失败 | 硬编码 baseline (复用 v1.4 经验) |
| Andrews p 值 Hansen 表精度 | 中 | p 值容差不达标 | 模拟法替代查表, 容差 1e-4 |
| Stock-Yogo 临界值表覆盖度 | 中 | L/K 组合缺失 | 见下方 §8.3 Stock-Yogo 覆盖范围核查 |
| Greeks Pathwise/LR 随机性 | 低 | 一致性检验不稳定 | 用置信区间比较, 增大 n_paths |
| hfecon → econometrics 依赖 (首次) | 低 | hfecon_diagnostics 调用 inference/ | ADR-015 方案 B 保证通用诊断无 Eigen3, 不污染 hfecon |
| ols_simple 与 ols.hpp 重复实现 | 低 | 维护两份 OLS 代码 | detail/ols_simple 仅 ~50-80 行 (参考 har_model.hpp Gauss-Jordan 模式, ADR-015 决策点 5), 规模 N≤千级, 性能非瓶颈 |

**并行策略** (ADR-015 方案 B 更新):

- **Wave 0 (前置)**: `detail/test_result_base.hpp` + `detail/ols_simple.hpp` — 必须先于其他通用诊断完成
- **Wave 1 (P0 并行)**: `residual_diagnostics` + `volatility_diagnostics` + `risk_diagnostics` 可全部并行 (volatility_diagnostics GARCH 测试用合成数据)
- **Wave 2 (P1 并行)**: `specification_tests` + `weak_identification` (estimation/) + `hfecon_diagnostics` 可并行 (hfecon_diagnostics 依赖 Wave 1 的 residual_diagnostics/volatility_diagnostics + Wave 2 的 specification_tests, 实际启动需在 specification_tests 接口稳定后)
- **Wave 3 (P2 并行)**: `structural_break` + `greeks_consistency` + `pricing_diagnostics` + `jump_test_diagnostics` 可全部并行
- **依赖约束**: `hfecon_diagnostics` 同时依赖 Wave 1 (residual_diagnostics/volatility_diagnostics) 和 Wave 2 (specification_tests 的 MincerZarnowitzResult, 见 §5.1 L786 include + L801 使用), 因此归入 Wave 2 而非 Wave 1

### 8.3 Stock-Yogo 临界值表覆盖范围核查 (2026-08-11, 已排幻觉)

**核查依据**:
- Stock & Yogo 2005, "Testing for Weak Instruments in Linear IV Regression", in *Identification and Inference for Econometric Models* (Andrews & Stock eds.), Cambridge University Press, pp. 80-108 (论文在书中的页码范围)
- 权威覆盖范围确认来源:
  - Stata 官方论坛 statalist.org 引用 (Table 1/2 的 K 覆盖范围)
  - Stata `ivreg2` 官方输出 (显著性阈值数值)
  - Huang-Wang-Yao 2023 (arXiv:2302.14423): 明确指出 "SY obtained a table of critical values for various combinations of (p, K_n, b), where p, K_n and b are up to 3, 30 and 0.3" — 即 K(内生)≤3, L(工具)≤30, bias阈值≤30%
  - Skeels-Windmeijer 2018 (Econometrics 6(4):44, doi:10.3390/econometrics6040044)

**记号约定** (SY 原书 vs spec):
- SY `n` (或 `p`) = spec `K` = 内生变量数 (n_endogenous)
- SY `K₂` (或 `K_n`) = spec `L` = 工具变量数 (n_instruments, excluded instruments)

**原书表覆盖范围** (仅确认 TSLS 两张表, LIML 表编号未在 search 结果中明确确认):

| 表号 | 准则 | 覆盖的 K (内生) | 覆盖的 L (工具) | 显著性阈值 |
|------|------|----------------|----------------|-----------|
| Table 1 (TSLS Bias) | 最大 IV 相对偏差 | K = 1, 2, **3** | L = 1..30 (Huang et al. 2023 确认上限) | 5%, 10%, 20%, 30% |
| Table 2 (TSLS Size) | 最大 Wald 检验 size 扭曲 | K = 1, **2** (不含 K=3) | L = 3..30 (需 L>K 才过度识别) | 10%, 15%, 20%, 25% |
| LIML 表 (编号未确认) | 同 Table 1/2 但 LIML 估计 | 推断同 TSLS, 但未严格核实 | — | — |

**注**: 原 spec 草稿曾写"四张表位于 pp. 58-61", 但论文在书中为 pp. 80-108 (29页), 页码矛盾, 已删除该具体页码引用 (可能源自论坛用户误引 2003 年工作论文版本)。

**与 spec 要求 L≤8, K≤3 的对照矩阵**:

| (K, L) 组合 | Table 1 (Bias) | Table 2 (Size) | 状态 |
|------------|---------------|----------------|------|
| (K=1, L=1..8) | ✓ 全覆盖 (L≤30) | ✓ 全覆盖 (L≥2 即可) | **完全覆盖** |
| (K=2, L=2..8) | ✓ 全覆盖 (L≤30) | ✓ 全覆盖 (L≥3 即可) | **完全覆盖** |
| (K=3, L=3..8) | ✓ 全覆盖 (L≤30, Huang et al. 2023 确认) | ✗ K=3 不在原表范围 | **部分覆盖** |

**关键缺口**:
1. **K=3 时 Size 准则 (Table 2) 完全不覆盖**: 原书 Table 2 仅覆盖 K=1, 2, 未提供 K=3 的 size distortion 临界值 (Stata 论坛权威引用确认)

**缓解策略** (两层回退):

1. **优先使用 Table 1 (Bias 准则)**: 覆盖 K≤3 且 L≤30 完整, 满足 spec 主要需求。接口默认 `criterion=RelativeBias(10%)`
2. **K=3 时 Size 准则降级**: 用 Skeels-Windmeijer 2018 "On the Stock-Yogo Tables" (Econometrics 6(4):44) 的解析方法。该论文给出 SY 临界值所源自期望的封闭形式解, 以及"在多内生变量场景下可能有价值" (论文原文: "may be of value in the presence of multiple endogenous regressors") 的二阶渐近近似。容差 1e-4 (对照原表已知点验证)

   **实现约束** (K=3 Size 准则强制取值):
   - `critical_value_is_exact` 必须设为 `false` (因使用解析近似, 非原表查表)
   - `critical_value_source` 必须设为 `"Skeels-Windmeijer 2018 approx"`
   - K=1,2 Size 准则 (原表覆盖) 时 `critical_value_is_exact=true`, `critical_value_source="SY2005 Table 2"`
   - K=1,2,3 Bias 准则 (原表覆盖) 时 `critical_value_is_exact=true`, `critical_value_source="SY2005 Table 1"`

**接口设计调整** (见 §2.4):

```cpp
enum class StockYogoCriterion {
    RelativeBias,    // SY2005 TSLS Bias 表, K≤3 全覆盖 (L≤30)
    SizeDistortion   // SY2005 TSLS Size 表, K≤2 原表, K=3 用 Skeels-Windmeijer 2018 近似
};

struct WeakIdentificationResult {
    Real cragg_donald_statistic;
    Size n_instruments;
    Size n_endogenous;
    StockYogoCriterion criterion;
    Real bias_threshold;              // 0.05/0.10/0.20/0.30 (RelativeBias)
    Real size_threshold;              // 0.10/0.15/0.20/0.25 (SizeDistortion)
    Real stock_yogo_critical_value;
    bool critical_value_is_exact;     // true=原表查表, false=解析近似
    std::string critical_value_source; // "SY2005 Table 1" 或 "Skeels-Windmeijer 2018 approx"
    bool reject_weak_instruments;
    std::string interpretation;
};
```

**测试矩阵补充** (原 10 用例 → 12 用例):
- CD 统计量(4)
- Stock-Yogo TSLS Bias 表查表 (K=1,2,3 × 代表性 L 值, 6 用例, 硬编码原表值; 注: K=k 时 L≥k)
- Skeels-Windmeijer 2018 近似 (K=3, Size 准则, 2 用例, 容差 1e-4)

---

## 附录 A: 排幻觉点完整清单 (23 项)

| 编号 | 模块 | 排幻觉点 | 验证方法 |
|------|------|---------|---------|
| H1 | JB | 分母用样本标准差 s, 非 σ | 单元测试 + R 对照 |
| H2 | JB | 公式用峰度 K, K-3 才是超额峰度 | 单元测试 |
| H3 | LB | 用 N(N+2)/(N-h) 加权, 非 N | R `Box.test` 对照 |
| H4 | LB | lag 自动选择 = min(10, N/5) | 边界测试 |
| H5 | BG | 辅助回归必须包含原 X | 维度检查 |
| H6 | BP | 默认 Koenker 修正 (e², 非 e²/σ²) | R `bptest` 对照 |
| H7 | White | 高维 q = K(K+1)/2, N > q 检查 | 边界测试 |
| H8 | 波动率 | z_t² 的 LB 是关键, 非仅 z_t | ARCH 效应测试 |
| H9 | 信息矩阵 | 对 QMLE 检验均值方程, 非方差方程 | QMLE 场景测试 |
| H10 | MZ | R² 也是预测精度度量 | 完美预测测试 |
| H11 | DM | HLN 修正用 1/N, t(N-1) 非 N(0,1) | R `dm.test` 对照 |
| H12 | Cragg-Donald | CD 是 F 的矩阵推广, 非 Wald | Stata 对照 |
| H13 | Stock-Yogo | Table 2 (Size) 仅覆盖 K≤2, K=3 用 Skeels-Windmeijer 2018 近似; 接口需标注 `critical_value_is_exact` | §8.3 对照矩阵测试 + 近似容差 1e-4 |
| H14 | CUSUM | 用递归残差, 非普通残差 | 递归估计验证 |
| H15 | Andrews | p 值用 Hansen 1997 非标准分布 | 模拟法验证 |
| H16 | DQ | X_t 含 [1, Hit_{t-1}, VaR_t], 非仅 Hit | X_t 构造验证 |
| H17 | Berkowitz | 需模型 CDF, 非经验 CDF | CDF 函数测试 |
| H18 | MC收敛 | 标准误差用批次均值法, 非单次估计 | 批次划分验证 |
| H19 | ES后验 | 需在 VaR 超越条件下, 非全样本 | 条件索引验证 |
| H20 | Greeks | Pathwise/LR 用置信区间比较, 非点估计 | 随机性测试 |
| H21 | Greeks | Gamma 二阶差分 bump size 更小 (1e-4) | bump size 测试 |
| H22 | IV拟合 | 权重用 Bid-Ask 宽度, 非简单方差 | 权重验证 |
| H23 | 多重修正 | BH 控制 FDR, 非 FWER | BH vs Bonferroni 对照 |

---

## 附录 B: 教材锚点完整清单

| 检验 | 教材 | 章节 |
|------|------|------|
| Jarque-Bera | Greene 8ed | §4.8 |
| Ljung-Box | Tsay 3ed | §2 |
| Breusch-Godfrey | Greene 8ed | §13.7 |
| Breusch-Pagan | Greene 8ed | §9.5 |
| White | Greene 8ed | §9.5 |
| 波动率诊断 | Tsay 3ed | Ch.3 |
| 波动率诊断 | McNeil-Frey-Embrechts 2005 | §5.3 |
| 信息矩阵 | Davidson-MacKinnon 2004 | §11.9 |
| Mincer-Zarnowitz | Mincer-Zarnowitz 1969 | 原文 |
| Diebold-Mariano | Diebold-Mariano 1995 | 原文 |
| Diebold-Mariano (HLN) | Harvey-Leybourne-Newbold 1997 | 原文 |
| Cragg-Donald | Cragg-Donald 1993 | 原文 |
| Stock-Yogo | Stock-Yogo 2005 | 原文 (Cambridge UP, 论文 pp. 80-108) |
| Stock-Yogo (扩展) | Skeels-Windmeijer 2018, Econometrics 6(4):44 | 原文 (doi:10.3390/econometrics6040044) |
| Stock-Yogo (覆盖范围) | Huang-Wang-Yao 2023, arXiv:2302.14423 | 原文 (确认 K≤3, L≤30) |
| CUSUM | Brown-Durbin-Evans 1975 | 原文 |
| CUSUM | Greene 8ed | §9.6 |
| Andrews | Andrews 1993 | 原文 |
| Andrews (p 值) | Hansen 1997 | 原文 |
| McNeil-Frey DQ | Engle-Manganelli 2004 | 原文 |
| McNeil-Frey DQ | McNeil-Frey-Embrechts 2005 | §5.3 |
| Berkowitz | Berkowitz 2001 | 原文 |
| Greeks | Glasserman 2003 | §7 |
| Greeks | Broadie-Glasserman 1996 | 原文 |
| IV surface | Gatheral 2006 | 原文 |
| IV surface | Fengler 2009 | 原文 |
| HAR | Corsi 2009 | 原文 |
| HEAVY | Shephard-Sheppard 2010 | 原文 |
| RV 预测评估 | Patton 2011 | 原文 |
| Bonferroni | Bonferroni 1936 | 原文 |
| BH | Benjamini-Hochberg 1995 | 原文 |

---

## 附录 C: v2.0+ 预留接口

```cpp
// conduction_metrics.hpp (v2.0+ scope, 仅预留空文件)
#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/core/linalg_dynamic.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// v2.0+ 信息论事前度量 (见 INFORMATION_THEORY_METRICS_RESEARCH.md)
// 当前版本不实现, 仅预留接口空文件
// - ConductionMatrix { MatrixXD M; Real tau; Size r_eff; }
// - compute_conduction_strength(model, param_grid) → Real τ
// - fisher_mapping(M) → MatrixXD I(θ)
// - effective_rank(M, threshold) → Size

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
```

---

## 附录 D: 与双层诊断体系的定位

```
Research OS 诊断体系 (本 spec 定位: L1/L2/L3 事后)
├── L0 事前 (信息论度量, v2.0+ scope, 见 INFORMATION_THEORY_METRICS_RESEARCH.md)
│   ├── τ (传导强度): 模型选择层, 跨模型可比
│   ├── r_eff (有效秩): 独立脆弱性通道数
│   └── v_max (最脆弱方向): 参数失效定位
│
├── L1 模型内 (事后统计量, v1.6 本 spec 实现)
│   ├── 计量: JB/LB/BG/BP/White (residual_diagnostics.hpp)
│   ├── 波动率: 标准化残差/z² LB (volatility_diagnostics.hpp)
│   ├── 高频: HAR 残差 LB (hfecon_diagnostics.hpp)
│   └── 风险: MC 收敛 (risk_diagnostics.hpp)
│
├── L2 跨模型 (事后统计量, v1.5 已有 + v1.6 补齐)
│   ├── v1.5 已有: Wald/LR/LM, Hansen J, LR 检验
│   ├── v1.6 新增: MZ 回归, DM 检验 (specification_tests.hpp)
│   ├── v1.6 新增: Greeks 跨方法一致性 (greeks_consistency.hpp)
│   ├── v1.6 新增: IV 拟合优度 (pricing_diagnostics.hpp)
│   └── v1.6 新增: 跳跃检验多重修正 (jump_test_diagnostics.hpp)
│
└── L3 元诊断 (事后统计量, v1.5 已有)
    ├── BCa 置信区间 (v1.5 已实现)
    ├── 信息矩阵检验 (v1.6 新增, specification_tests.hpp)
    └── 弱识别检验 (v1.6 新增, weak_identification.hpp)
```

---

**Spec 状态**: v2.0 重新设计版 (覆盖全部已实现模块, ADR-015 已 Accepted)
**核查状态**:
- ✅ 全部 4 个模块群 (risk/pricing/hfecon/econometrics) 已覆盖
- ✅ 11 个新增头文件路径与现有代码结构对齐 (LS 验证)
- ✅ 23 个排幻觉点覆盖调研报告 §6.1-6.4 全部缺失统计量
- ✅ 与 v2.0+ 信息论度量框架边界清晰 (附录 D)
- ✅ ADR-015 (证伪统计量模块边界) 已 Accepted (2026-08-12), 方案 B 全部 5 个决策点已同步至本 spec

**依赖状态**:
1. ~~ADR-015: 证伪统计量模块边界~~ — ✅ **Accepted** (2026-08-12, 见 [ADR_INDEX.md §ADR-015](../../decisions/ADR_INDEX.md#adr-015-证伪统计量模块边界-通用-vs-模块特定)), 方案 B: 通用诊断不依赖 Eigen3
2. Stock-Yogo 临界值表的覆盖范围 — ✅ **核查完成** (见 §8.3, K=3 Size 准则用 Skeels-Windmeijer 2018 解析近似)
3. Phase 7B (金融时间序列) spec 编写 — ⏳ 待启动 (volatility_diagnostics GARCH 测试用合成数据可先行)

**下一步**: 进入实施, 按 W0-W3 里程碑推进 (Wave 0 公共基础设施 → Wave 1 P0 → Wave 2 P1 → Wave 3 P2)
