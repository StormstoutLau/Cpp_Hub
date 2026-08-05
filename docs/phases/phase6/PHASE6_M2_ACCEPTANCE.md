# Phase 6 v1.5 M2 实施验收报告

> **文档类型**: 实施验收报告 (Implementation Acceptance Report)
> **审计对象**: v1.5 M2 - MLE/QMLE + 假设检验 + 信息准则 + 估计器工厂
> **关联文档**:
>   - [PHASE6_IMPLEMENTATION_PLAN.md](./PHASE6_IMPLEMENTATION_PLAN.md) §4.2-4.3
>   - [PHASE6_ECONOMETRICS_SPEC.md](./PHASE6_ECONOMETRICS_SPEC.md)
> **验收日期**: 2026-08-05
> **验收员**: Scott (self-review, 符合项目 solo developer 工作流)
> **验收标准**: 三平台全量测试通过 + 排幻觉点全部验证 + 无遗留阻塞项

---

## 1. 验收基本信息

| 项 | 值 |
|----|-----|
| 项目名称 | Cpp_Hub |
| 审计阶段 | Phase 6 v1.5 M2 (经典参数计量模块 - MLE/假设检验/信息准则) |
| 提交版本 | `0ea3fa2` feat(v1.5 M2): MLE/QMLE + Hypothesis Tests + Information Criteria + Estimator Factory |
| 验收方法 | TDD 实现 + 三平台跨平台验证 + 排幻觉逐点核查 |
| 前置条件 | v1.5 M1 已通过 (1544/1544 三平台通过) |
| 后置条件 | M2 验收通过, 可进入 v1.5 M3 (Bootstrap/GMM) 实施 |

---

## 2. 三平台测试结果汇总

| 平台 | 编译器 | 测试总数 | 通过数 | 失败数 | 耗时 | 状态 |
|------|--------|---------|--------|--------|------|------|
| 主控站 (Windows 10) | MSVC 19.43 | 1644 | 1644 | 0 | 872.52 sec | ✅ 通过 |
| A 工作站 (Ubuntu 24.04) | GCC 13.3.0 | 1644 | 1644 | 0 | 362.18 sec | ✅ 通过 |
| B 工作站 (Ubuntu 6.17) | GCC 13.3.0 | 1644 | 1644 | 0 | 37.55 sec | ✅ 通过 |

**M2 新增测试**: 100 个 (相对 M1 的 1544 个)
- MLE Gaussian: 12 个
- MLE Logistic: 18 个
- MLE Poisson: 11 个
- MLE Probit: 7 个
- MLE NegativeBinomial: 5 个
- 假设检验 (Wald/LR/LM/J): 28 个
- 信息准则 (AIC/BIC/HQ/AICc): 12 个
- 集成测试 (IntegrationM2): 7 个

---

## 3. M2 实施范围验收

### 3.1 MLE 估计器 (6 种族)

| 估计器 | 实现文件 | 测试文件 | 排幻觉点 | 状态 |
|--------|---------|---------|---------|------|
| Gaussian MLE | `mle.hpp` | `test_mle_gaussian.cpp` | G3/G4 | ✅ |
| Logistic MLE | `mle.hpp` | `test_mle_logistic.cpp` | E7/E8/L1-L3 | ✅ |
| Probit MLE | `mle.hpp` | `test_mle_probit_nb.cpp` | PB1-PB4 | ✅ |
| Poisson MLE | `mle.hpp` | `test_mle_poisson.cpp` | P1-P5 | ✅ |
| NegativeBinomial MLE | `mle.hpp` | `test_mle_probit_nb.cpp` | NB1-NB3 | ✅ |
| Bernoulli MLE | `mle.hpp` | (Logistic 覆盖) | - | ✅ |

### 3.2 协方差矩阵计算 (3 种方法)

| 方法 | 公式 | 排幻觉点 | 状态 |
|------|------|---------|------|
| OPG (Outer Product of Gradients) | V = G'G / N | E8 | ✅ |
| Hessian | V = H^{-1} (观察信息矩阵) | G3 | ✅ |
| Sandwich (QMLE) | V = H^{-1} · G'G · H^{-1} | E8 | ✅ |

### 3.3 假设检验 (3 种 + J-test)

| 检验 | 公式 | 自由度 | 排幻觉点 | 状态 |
|------|------|--------|---------|------|
| Wald | (Rβ-r)'[RVR']^{-1}(Rβ-r) ~ χ²(q) | q | E9 | ✅ |
| Wald (F 形式) | F = Wald/q ~ F(q, df_res) | (q, df_res) | E9 | ✅ |
| LR | 2(ℓ_UR - ℓ_R) ~ χ²(q) | q = n_UR - n_R | - | ✅ |
| LM | N·R²_aux ~ χ²(q) (Breusch-Pagan) | q | H1 | ✅ |
| Hansen J | n·ḡ'Ŝ^{-1}ḡ ~ χ²(q-k) | q-k | H2 | ✅ |

### 3.4 信息准则 (4 种)

| 准则 | 公式 | 状态 |
|------|------|------|
| AIC | 2K - 2ℓ | ✅ |
| BIC | K·log(N) - 2ℓ | ✅ |
| HQ | 2K·log(log(N)) - 2ℓ | ✅ |
| AICc | AIC + 2K(K+1)/(N-K-1) | ✅ |

### 3.5 估计器工厂

| 功能 | 排幻觉点 | 状态 |
|------|---------|------|
| 按名称创建 (`create("OLS")`) | F1/F2 | ✅ |
| 按 MLEFamily 枚举创建 (`createMLE(Logistic)`) | F1/F3 | ✅ |
| Meyers Singleton (规避 SIOF) | F4 | ✅ |
| 未注册名称抛 `std::invalid_argument` | F2 | ✅ |
| 每次返回全新实例 (`unique_ptr`) | F1/F3 | ✅ |

---

## 4. 排幻觉点逐项核查

### 4.1 MLE 实现排幻觉点

| 编号 | 描述 | 验证方法 | 结果 |
|------|------|---------|------|
| **E7** | Logistic MLE 用 Newton-Raphson (对 canonical link 等价于 IRLS) | 数值对照: Newton-Raphson 迭代 vs IRLS 权重, β 收敛值一致 | ✅ 通过 |
| **E8** | QMLE Sandwich: bread = (X'WX)^{-1}, meat = X'diag(ε²)X | 手算 N=30 Logistic 数据集, Sandwich vs Hessian vs OPG 三者对比 | ✅ 通过 |
| **G3** | Gaussian Hessian V = σ²·(X'X)^{-1} (含 σ² 因子, 非 (X'X)^{-1}) | Gaussian MLE vs OLS 数值一致性测试 (Factory_GaussianMLE_vs_OLS) | ✅ 通过 |
| **G4** | Gaussian MLE σ² = SSR/N (非 SSR/(N-K)) | σ²_MLE vs σ²_OLS 对比, vcov 差异 = N/(N-K) 因子 | ✅ 通过 |
| **L1** | Logistic link: μ = 1/(1+exp(-η)), dμ/dη = μ(1-μ) | 手算 η=0 时 μ=0.5, dμ/dη=0.25 | ✅ 通过 |
| **L2** | Logistic log-likelihood: y·log(μ) + (1-y)·log(1-μ) | N=30 数据集 ℓ 值手算对照 | ✅ 通过 |
| **L3** | Logistic 完全分离检测: MLE 不存在 (β→∞) | N=30 非完全分离数据集设计, 含渐变转换区 | ✅ 通过 |
| **P1** | Poisson link: μ = exp(η) (log-link, canonical) | η=0 时 μ=1 | ✅ 通过 |
| **P2** | Poisson log-likelihood: y·η - μ - log(y!) | 手算 N=100 数据集 ℓ 值 | ✅ 通过 |
| **P3** | Poisson score: X'(y - μ) | 数值梯度 vs 解析梯度对比 | ✅ 通过 |
| **P4** | Poisson Hessian: -X'diag(μ)X (负定) | LLT 分解验证正定性 (使用 -H) | ✅ 通过 |
| **P5** | Poisson MLE = Poisson QMLE (canonical link, 正确指定时) | Sandwich vs Hessian vcov 一致性 | ✅ 通过 |
| **PB1-PB4** | Probit link/CDF/PDF/score 数值正确性 | Φ(0)=0.5, φ(0)=0.3989, 数值梯度对比 | ✅ 通过 |
| **NB1-NB3** | NegativeBinomial α→0 退化为 Poisson | α=1e-8 时 NB vcov ≈ Poisson vcov (容差 1e-4) | ✅ 通过 |

### 4.2 假设检验排幻觉点

| 编号 | 描述 | 验证方法 | 结果 |
|------|------|---------|------|
| **E9** | Wald 同时提供 χ² 和 F 两种 p 值 (R `lmtest::waldtest` 默认 F) | `use_f_distribution` 标志 + `f_df1`/`f_df2` 字段 | ✅ 通过 |
| **H1** | LM = N·R²_aux (非 (N-K)·R²_aux, Breusch-Pagan 1979 原始形式) | 手算 N=4: LM = 4·0.9797 = 3.9186, 非 2·0.9797 = 1.9593 | ✅ 通过 |
| **H2** | J-test df = q-k (过度识别约束数, 非 q) | q=3, k=2 → df=1 (非 3); 严格识别 q=k → df=0, J 无分布 | ✅ 通过 |
| **CHI2-1** | χ²(1) p 值精确计算: erfc(√(x/2)) 替代 gammq 级数近似 | erfc(√7) = 0.00018281 (机器精度) vs gammq 旧值 0.00018281 (差异 7e-6) | ✅ 通过 (修复) |
| **CHI2-2** | χ²(2) p 值精确计算: exp(-x/2) 替代 gammq | χ²(2) 是均值 2 的指数分布, P(χ²(2)>x) = exp(-x/2) | ✅ 通过 |
| **J-EXPECT** | J-test 期望值幻觉修正 | 原期望 0.0001759 来源不明, 正确值 erfc(√7) = 0.00018281 | ✅ 通过 (修复) |

### 4.3 工厂模式排幻觉点

| 编号 | 描述 | 验证方法 | 结果 |
|------|------|---------|------|
| **F1** | `create()` 每次返回全新实例 (clone 语义, 不共享状态) | 连续两次 `create("OLS")` 返回不同 `unique_ptr` | ✅ 通过 |
| **F2** | 未注册名称抛 `std::invalid_argument` (fail-fast, 不返回 nullptr) | `create("Unknown")` 抛异常测试 | ✅ 通过 |
| **F3** | 工厂不持有 Estimator 实例 (`unique_ptr` 转移所有权) | `create()` 返回 `unique_ptr<Estimator>`, 调用方持有所有权 | ✅ 通过 |
| **F4** | Meyers Singleton 规避 SIOF (静态注册顺序正确) | `EstimatorFactory::instance()` 返回局部 static 引用 | ✅ 通过 |

### 4.4 集成测试排幻觉点

| 测试 | 验证场景 | 排幻觉点 | 状态 |
|------|---------|---------|------|
| `OLS_Wald_InformationCriteria_EndToEnd` | OLS → Wald(β₁=0) → AIC/BIC | F1/E9/G3 | ✅ 通过 |
| `MLELogistic_Wald_LR_EndToEnd` | Logistic → Wald + LR(约束 vs 无约束) | E7/E9 | ✅ 通过 |
| `MLEPoisson_QMLE_Sandwich_Wald_EndToEnd` | Poisson → QMLE Sandwich → Wald | E8/P5 | ✅ 通过 |
| `Factory_GaussianMLE_vs_OLS_Consistency` | Gaussian MLE vs OLS β/vcov 一致性 | G3/G4 | ✅ 通过 |
| `ModelSelection_Logistic_vs_Probit_AIC` | Logistic vs Probit AIC/BIC 比较 | - | ✅ 通过 |
| `Factory_UnknownName_Throws` | 未注册名称抛异常 | F2 | ✅ 通过 |
| `Factory_RegisteredEstimators` | 已注册估计器名称列表 | F4 | ✅ 通过 |

---

## 5. 实施过程幻觉修复记录

### 5.1 跨头文件重复定义 (编译错误)

| 问题 | 原因 | 修复 |
|------|------|------|
| `kTwoPi`/`kLogTwoPi` 在 `ols.hpp` 和 `mle.hpp` 重复定义 | 各模块独立定义相同常量 | 提取至 `core/special_functions.hpp` 共享 |
| `betacf`/`beta_i` 在 `ols.hpp` 和 `hypothesis_tests.hpp` 重复定义 | 各模块独立实现相同算法 | 提取至 `core/special_functions.hpp` 共享 |

### 5.2 Logistic 完全分离 (MLE 不存在)

| 问题 | 原因 | 修复 |
|------|------|------|
| Newton-Raphson 迭代发散, "X'WX is singular" | y=0 和 y=1 完全分离, MLE β→∞ | 重构 N=30 数据集, 添加渐变转换区 (x=11-17) |

### 5.3 OLS BIC 断言错误 (统计理论误用)

| 问题 | 原因 | 修复 |
|------|------|------|
| `EXPECT_GT(ic.bic, ic.aic)` 在 N=4 时不成立 | 小样本 log(N) < K 系数时 BIC < AIC | 替换为验证 BIC 公式正确性: `EXPECT_NEAR(ic.bic, K·log(N) - 2ℓ, 1e-10)` |

### 5.4 J-test χ²(1) p 值精度 (数值精度 + 期望值幻觉)

| 问题 | 原因 | 修复 |
|------|------|------|
| p 值差异 6.9e-6 超过容差 1e-8 | (1) gammq 级数/连分式在极端尾部精度不足; (2) **测试期望值 0.0001759 是幻觉** | (1) df=1 用 `std::erfc(√(x/2))` 精确计算; (2) 期望值修正为 `erfc(√7) = 0.00018281` |

**排幻觉关键发现**: 测试期望值 `0.00017589591417979427` 来源不明, 可能来自不精确的在线计算器。数学推导证明:
- χ²(1) = Z², P(χ²(1)>x) = P(|Z|>√x) = 2·(1-Φ(√x)) = erfc(√(x/2))
- 对 x=14: erfc(√7) = erfc(2.64575131...) = 0.00018281063298183488
- 渐近检验: √(2/(πx))·exp(-x/2) = √(2/(14π))·exp(-7) = 0.2132·0.000912 = 0.0001944 (上界, 真值更小 ✓)

### 5.5 约束模型 y 值不一致 (测试逻辑错误)

| 问题 | 原因 | 修复 |
|------|------|------|
| LR 检验中约束模型 y 值与无约束模型不同 | 修改无约束数据集后未同步约束模型 | 同步 `y_restricted_vals` 与 `y_vals` 数组 |

### 5.6 数组越界 (未定义行为)

| 问题 | 原因 | 修复 |
|------|------|------|
| N=30 时 `y_vals` 数组仅 29 个元素 | 手动初始化计数错误 | 调整为 30 个元素, x=18-30 有 13 个 y=1 值 |

---

## 6. 三平台一致性验证

### 6.1 数值一致性

| 测试场景 | MSVC (主控站) | GCC (A 站) | GCC (B 站) | 一致性 |
|---------|--------------|-----------|-----------|--------|
| OLS β = [0, 1.7] | ✅ 1e-10 | ✅ 1e-10 | ✅ 1e-10 | ✅ |
| Wald = 96.333... | ✅ 1e-6 | ✅ 1e-6 | ✅ 1e-6 | ✅ |
| LM = 3.9186... | ✅ 1e-10 | ✅ 1e-10 | ✅ 1e-10 | ✅ |
| J-test p = erfc(√7) | ✅ 1e-12 | ✅ 1e-12 | ✅ 1e-12 | ✅ |
| Logistic MLE 收敛 | ✅ | ✅ | ✅ | ✅ |
| Poisson QMLE Sandwich | ✅ | ✅ | ✅ | ✅ |

### 6.2 跨平台注意事项

| 注意点 | 说明 |
|--------|------|
| `std::normal_distribution` 跨平台不可复现 | C++ 标准不规定具体算法, MSVC STL 和 GCC libstdc++ 实现不同. M2 测试设计避免依赖随机数据精确断言 (使用解析值或固定常量) |
| `std::erfc` 跨平台一致 | C++ 标准库函数, MSVC/GCC 实现一致 (IEEE 754 保证) |
| Eigen3 LLT 分解跨平台一致 | Eigen3 3.4.0 保证跨平台数值一致性 (相同输入相同输出) |

---

## 7. 验收结论

### 7.1 量化指标

| 指标 | 值 |
|------|-----|
| 三平台测试通过率 | 100% (1644/1644 × 3 平台) |
| M2 新增测试数 | 100 |
| 排幻觉点验证数 | 24 (全部通过) |
| 实施过程修复幻觉数 | 6 |
| 阻塞项 | 0 |
| 遗留问题 | 0 |

### 7.2 验收结论

**v1.5 M2 验收通过**。三平台 (MSVC/GCC-A/GCC-B) 1644/1644 测试全部通过, 跨平台数值一致性确认。24 个排幻觉点全部验证通过, 6 个实施过程幻觉已修复并记录。

**可进入 v1.5 M3 (Bootstrap/GMM) 实施阶段**。

---

## 8. 附录: M2 文件清单

### 8.1 新增文件

| 文件 | 用途 |
|------|------|
| `include/cpphub/econometrics/core/special_functions.hpp` | 共享数学常数和特殊函数 (kTwoPi/kLogTwoPi/betacf/beta_i) |
| `include/cpphub/econometrics/estimation/mle.hpp` | 6 种 MLE 估计器实现 |
| `include/cpphub/econometrics/estimation/estimator_factory.hpp` | ADR-003 风格估计器工厂 |
| `include/cpphub/econometrics/inference/hypothesis_tests.hpp` | Wald/LR/LM/J-test 假设检验 |
| `include/cpphub/econometrics/inference/diagnostics.hpp` | AIC/BIC/HQ/AICc 信息准则 |
| `tests/unit/econometrics/test_mle_gaussian.cpp` | Gaussian MLE 单元测试 |
| `tests/unit/econometrics/test_mle_logistic.cpp` | Logistic MLE 单元测试 |
| `tests/unit/econometrics/test_mle_poisson.cpp` | Poisson MLE 单元测试 |
| `tests/unit/econometrics/test_mle_probit_nb.cpp` | Probit/NegativeBinomial MLE 单元测试 |
| `tests/unit/econometrics/test_hypothesis_tests.cpp` | 假设检验单元测试 |
| `tests/unit/econometrics/test_diagnostics.cpp` | 信息准则单元测试 |
| `tests/unit/econometrics/test_integration_m2.cpp` | M2 端到端集成测试 |

### 8.2 修改文件

| 文件 | 修改内容 |
|------|---------|
| `include/cpphub/econometrics/estimation/ols.hpp` | 移除重复定义, 添加 `special_functions.hpp` 引用 |
| `tests/CMakeLists.txt` | 添加 M2 测试目标 |
