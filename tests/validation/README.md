# Cpp_Hub 测试基准索引

> **G4 审计项**: 本文件列出所有数值测试的基准来源、对应文件及容差,
> 满足 AUDIT_CHECKLIST.md G4 "基准索引完整" 要求。
>
> 维护规则: 新增测试必须在此登记基准来源; 无基准来源的数值测试不予合并。

---

## 1. 按基准来源分类

### 1.1 解析公式 (Black-Scholes-Merton)

| 测试文件 | 测试用例 | 基准来源 | 容差 |
|----------|----------|----------|------|
| `unit/pricing/test_bsm_analytic.cpp` | BSM 欧式期权解析解 | Black-Scholes (1973) 公式闭式解 | 1e-12 |
| `unit/risk/greeks/test_greeks_analytic_numerical.cpp` | 解析 Greeks vs 数值 Greeks | BSM Greeks 闭式公式 | 1e-10 |
| `unit/risk/greeks/test_aad_greeks.cpp` | AAD Greeks vs 解析 Greeks | BSM Greeks 闭式公式 | 1e-10 |
| `unit/risk/greeks/test_pathwise_greeks.cpp` | Pathwise Delta/Vega | BSM Greeks 闭式公式 (GBM 路径法) | 1e-2 (MC 噪声) |
| `unit/risk/greeks/test_lr_greeks.cpp` | LR Delta/Vega | BSM Greeks 闭式公式 (似然比法) | 2e-2 (MC 噪声) |
| `unit/validation/test_python_cross_validation.cpp` | Parametric VaR (Normal) | scipy.stats.norm.ppf | 1e-12 |

**公式来源**:
- Black, F., Scholes, M. (1973). "The Pricing of Options and Corporate Liabilities"
- Merton, R. C. (1973). "Theory of Rational Option Pricing"

### 1.2 论文基准表

| 测试文件 | 测试用例 | 基准来源 | 容差 |
|----------|----------|----------|------|
| `unit/pricing/test_heston_cf.cpp` | Heston CF Schoutens 表 (u=0.5, 1.0, 2.0) | Schoutens (2003) §2.5 特征函数表 | 1e-10 |
| `unit/models/test_heston_process.cpp` | Heston CF vs Schoutens 表 | Schoutens (2003) "Lévy Processes in Finance" | 1e-10 |
| `unit/pricing/test_pde_engine.cpp` | 美式看跌 Broadie-Detemple | Broadie & Detemple (1996) 美式期权基准表 | 1e-3 |
| `unit/models/vol_surface/test_ssvi.cpp` | SSVI 无套利条件 | Gatheral-Jacquier (2014) Theorem 4.2/4.4 | 布尔 (完全一致) |
| `unit/validation/test_python_cross_validation.cpp` | SSVI Power-law 公式 | Gatheral-Jacquier (2014) eq. 5.3 | 1e-12 (位精确) |

**论文引用**:
- Schoutens, W. (2003). "Lévy Processes in Finance: From Financial Models to Financial Engineering"
- Broadie, M., Detemple, J. (1996). "American Option Valuation: New Bounds, Approximations, and a Comparison of Existing Methods"
- Gatheral, J., Jacquier, A. (2014). "Arbitrage-free SVI volatility surfaces" ([arXiv:1204.0646](https://arxiv.org/abs/1204.0646))
- Heston, S. L. (1993). "A Closed-Form Solution for Options with Stochastic Volatility"

### 1.3 Python/R 开源库交叉验证

| 测试文件 | 验证模块 | 交叉验证源 | 容差 |
|----------|----------|------------|------|
| `unit/validation/test_python_cross_validation.cpp` | Historical VaR (Linear/Conservative/Empirical) | `tests/validation/python/cross_validate_var.py` (numpy) | 1e-12 (算法一致) |
| `unit/validation/test_python_cross_validation.cpp` | Parametric VaR (Normal) | scipy.stats.norm.ppf | 1e-12 |
| `unit/validation/test_python_cross_validation.cpp` | SSVI Power-law φ(θ) / 总方差 / 隐含波动率 | `tests/validation/python/cross_validate_calibration.py` (scipy.optimize) | 1e-12 (位精确) |
| `unit/validation/test_python_cross_validation.cpp` | SSVI 无套利条件 (Calendar/Butterfly) | Python `no_arbitrage_check` | 布尔 (完全一致) |
| `unit/validation/test_python_cross_validation.cpp` | Kupiec POF 统计量 | Python 标准实现 (χ²(1)) | 1e-3 |

**Python 脚本**:
- `tests/validation/python/cross_validate_var.py`: 生成 VaR 基准值 (numpy/scipy/statsmodels)
- `tests/validation/python/cross_validate_calibration.py`: 生成标定基准值 (scipy.optimize)
- `tests/validation/benchmarks_var.json`: VaR 基准值输出 (seed=42, n=5000)
- `tests/validation/benchmarks_calib.json`: 标定基准值输出 (36 点期权链)

**重新生成基准值**:
```bash
cd tests/validation
python python/cross_validate_var.py --output benchmarks_var.json
python python/cross_validate_calibration.py --output benchmarks_calib.json
```

### 1.4 内部一致性基准 (无外部源)

| 测试文件 | 验证内容 | 基准来源 | 容差 |
|----------|----------|----------|------|
| `unit/core/test_core_math.cpp` | normal_pdf/cdf/inv_cdf 往返 | Abramowitz-Stegun 7.1.26 + Acklam 逆 CDF | 1e-10 |
| `unit/core/test_rng.cpp` | Philox4x64 统计性质 | χ² 拟合优度 (均匀性) | 统计 (p>0.01) |
| `unit/core/test_linalg.cpp` | Cholesky/Gauss-Jordan 解线性系统 | 已知矩阵闭式解 | 1e-12 |
| `unit/monte_carlo/test_sobol_qmc.cpp` | Sobol 序列均匀性 | van der Corput 性质 | 1e-12 |
| `unit/monte_carlo/test_variance_reduction.cpp` | CV/MM/AV 方差缩减比 | 理论方差比 (解析) | 1e-2 (MC 噪声) |
| `unit/risk/var/test_var_engine.cpp` | VaR 引擎 (历史/参数/MC) | 内部算法一致性 | 1e-10 |
| `unit/risk/var/test_backtesting.cpp` | Kupiec/Christoffersen/Basel 回测 | 标准 χ² 临界值 | 1e-10 |
| `unit/risk/scenario/test_stress_sensitivity.cpp` | 压力测试/敏感度 | 已知场景闭式解 | 1e-10 |
| `unit/calibration/test_calibration_framework.cpp` | LM/NelderMead/DE 优化器 | Rosenbrock/Rastrigin 已知最优 | 1e-6 |
| `unit/calibration/test_svi_diagnostic.cpp` | SVI butterfly 套利 | Gatheral-Jacquier g(k)≥0 | 1e-10 |
| `unit/models/vol_surface/test_vol_surface.cpp` | 波动率曲面 (占位) | — | — |
| `unit/models/vol_surface/test_dupire_local_vol.cpp` | Dupire 局部波动率 | Dupire (1994) 公式 | 1e-10 |
| `unit/instruments/test_payoff_bridge.cpp` | Payoff 桥接模式 | 内部一致性 | 精确 |
| `unit/instruments/test_payoff_factory.cpp` | Payoff 工厂 | 内部一致性 | 精确 |
| `unit/models/test_gbm.cpp` | GBM 过程 | S_T = S_0 exp((r-σ²/2)T + σ√T Z) | 1e-12 |
| `unit/pricing/test_tree_engine.cpp` | 二叉树/三叉树 | Cox-Ross-Rubinstein (1979) | 1e-3 (收敛阶) |
| `unit/integration/test_integration_phase2.cpp` | Phase 2 端到端 | 内部 pipeline 一致性 | 1e-6 |
| `unit/integration/test_integration_phase3.cpp` | Phase 3 端到端 | 内部 pipeline 一致性 | 1e-6 |
| `unit/risk/greeks/test_greeks_factory.cpp` | Greeks Factory 分发 | Auto 策略一致性 | 1e-10 |

---

## 2. 按容差等级分类

### 2.1 位精确 (1e-12) — 解析公式直接计算
- BSM 欧式期权解析解
- Heston 特征函数 (vs Schoutens 表)
- SSVI Power-law 总方差/隐含波动率 (vs Python scipy)
- Parametric VaR Normal (vs scipy.stats.norm.ppf)
- Historical VaR 算法 (C++ vs Python 独立实现,同 PnL 输入)
- normal_cdf/inv_normal_cdf 往返
- Cholesky 分解 / Gauss-Jordan 求解

### 2.2 高精度 (1e-10) — 数值方法 vs 解析解
- AAD Greeks vs BSM 解析 Greeks (1e-10)
- 解析 Greeks vs 数值差分 Greeks (1e-10)
- Dupire 局部波动率恢复 (1e-10)
- SVI butterfly 套利 g(k)≥0 检查 (1e-10)
- Kupiec/Christoffersen 回测统计量 (1e-10)

### 2.3 MC 收敛容差 (1e-2 ~ 5e-2) — 蒙特卡洛噪声
- Pathwise Delta/Vega (1e-2, 10000 路径)
- LR Delta/Vega (2e-2 ~ 5e-2, 10000 路径)
- 方差缩减比验证 (1e-2)
- MC VaR 统计一致性 (2e-2, 不同 RNG)

### 2.4 PDE/树收敛容差 (1e-3) — 离散化误差
- 美式期权 PDE vs Broadie-Detemple (1e-3, 500×2000 网格)
- 二叉树/三叉树 vs 解析解 (1e-3, 1000 步)

---

## 3. 跨平台一致性

### 3.1 三平台验证 (MSVC Win10 + GCC-A NEX + GCC-B GTR-Pro)
- **浮点确定性**: `-ffp-contract=off` (GCC) + `/fp:precise` (MSVC) 保证 IEEE-754 严格模式
- **位精确一致**: 所有 286 测试在三平台输出完全相同的浮点位
- **性能对比**: MSVC 9.87s / GCC-A 5.89s / GCC-B 5.66s (全量回归)

### 3.2 Phase 4 LITE 新增验证
- **Python 交叉验证** (G2): 16 测试全绿,C++ 与 Python (numpy/scipy) 位精确一致
- **SSVI 跨期限** (D2): 17 测试全绿,无套利条件满足 Gatheral-Jacquier Theorem 4.2

---

## 4. 待补充基准 (Phase 4+ 整改留项)

| 模块 | 待补充基准 | 优先级 | 计划 |
|------|-----------|--------|------|
| GPU MC | CPU vs GPU 位精确/1e-12 | 中 | Phase 4 M3 |
| Heston MC | vs QuantLib/AAD 比较 | 低 | v2.0 |
| SABR 标定 | vs R `sabr` 包 | 低 | v2.0 |
| GARCH VaR | vs `arch` Python 包 | 低 | v2.0 |

---

## 5. 基准来源引用清单

### 论文
1. Black, F., Scholes, M. (1973). "The Pricing of Options and Corporate Liabilities"
2. Merton, R. C. (1973). "Theory of Rational Option Pricing"
3. Heston, S. L. (1993). "A Closed-Form Solution for Options with Stochastic Volatility"
4. Cox, J. C., Ross, S. A., Rubinstein, M. (1979). "Option Pricing: A Simplified Approach"
5. Broadie, M., Detemple, J. (1996). "American Option Valuation: New Bounds, Approximations"
6. Schoutens, W. (2003). "Lévy Processes in Finance"
7. Gatheral, J., Jacquier, A. (2014). "Arbitrage-free SVI volatility surfaces" (arXiv:1204.0646)
8. Dupire, B. (1994). "Pricing with a Smile"
9. Abramowitz, M., Stegun, I. (1964). "Handbook of Mathematical Functions" §7.1.26
10. Acklam, P. J. (2003). "An Algorithm for Computing the Inverse Normal CDF"

### 开源库
1. **NumPy** 2.x: 数组计算、分位数 (numpy.quantile)
2. **SciPy** 1.x: scipy.stats.norm (正态分布)、scipy.optimize.least_squares (LM 标定)
3. **statsmodels**: 描述性统计、分位数回归 (可选依赖)
4. **autodiff** (autodiff.github.io): 前向/反向自动微分 (C++ 第三方库,用于 Greeks 交叉验证)

### 内部基准
1. **scripts/generate_bs_benchmark.py**: BSM 基准值生成 (scipy.stats)
2. **tests/validation/python/cross_validate_var.py**: VaR 基准值生成
3. **tests/validation/python/cross_validate_calibration.py**: 标定基准值生成

---

**维护人**: Scott (鹏)
**最后更新**: 2026-07-31 (Phase 4 LITE M1-4)
**审计状态**: G4 基准索引完整 ✅
