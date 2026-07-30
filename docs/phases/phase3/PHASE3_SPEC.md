# Phase 3 执行规格书 - 风险管理与模型标定

> **版本归属**: **v1.1** (v1.1 = Phase 3 完整 scope + Phase 2 推迟内容)
> **目标**: 3 周内交付完整 Greeks 体系(AAD/Pathwise/LR/FD)、VaR/ES 引擎、生产级标定框架、波动率曲面建模
> **前置**: Phase 1-2 (v1.0) 全部通过
> **里程碑**: M1(Week 1) Greeks/AAD → M2(Week 2) VaR/ES/场景分析 → M3(Week 3) 标定/波动率曲面/验收
>
> **v1.1 Scope 声明** (2026-07-29 评审):
> - 本 Phase 包含原 Phase 3 完整内容 + 原 Phase 2 推迟的 v1.1 内容
> - **Phase 2 推迟到 v1.1 的内容** (本 Phase 一并实现):
>   - 进阶模型: SABR/Bates/VG/CEV、利率模型 (Vasicek/CIR/HW/G2++)
>   - 傅里叶引擎: COS/FFT/CONV
>   - LSMC 美式/百慕大引擎
>   - 高级方差缩减: Importance Sampling/Stratified Sampling/Conditional MC
>   - 波动率曲面: SVI/SSVI/Dupire
>   - 标定框架: LM/DE/Nelder-Mead
> - **AAD 实现**: 采用 autodiff 库 (MIT, header-only)，不自研 Tape (ADR-007 Revised)

---

## 1. 交付物清单

### 1.1 新增编译目标
| 目标 | 类型 | 关键源文件 |
|------|------|------------|
| `cpphub` (增量) | Shared Library | `src/risk/*.cpp`, `src/calibration/*.cpp`, `src/models/vol_surface/*.cpp` |
| `cpphub_risk` | Static Library | 可选单独链接的风险模块 |

### 1.2 必须通过的新增测试
| 测试套件 | 用例数 | 覆盖模块 |
|----------|--------|----------|
| `test_aad_greeks` | 20 | `aad_greeks.hpp` (BS/Heston/篮子期权) |
| `test_pathwise_greeks` | 15 | `pathwise_greeks.hpp` (欧式/障碍/亚式) |
| `test_lr_greeks` | 10 | `likelihood_ratio_greeks.hpp` |
| `test_fd_greeks` | 10 | `fd_greeks.hpp` (中心差分/复阶差分) |
| `test_greeks_consistency` | 12 | 所有方法对比一致性 (相对差 < 1e-6) |
| `test_var_engine` | 15 | `historical_var.hpp`, `parametric_var.hpp`, `mc_var.hpp` |
| `test_es_engine` | 10 | `expected_shortfall.hpp` (Cornish-Fisher/MC) |
| `test_scenario_analysis` | 10 | `stress_test.hpp`, `sensitivity.hpp` |
| `test_calibration_framework` | 20 | `lm_calibrator.hpp`, `de_calibrator.hpp`, `multi_objective.hpp` |
| `test_vol_surface` | 15 | `svi.hpp`, `ssvi.hpp`, `dupire_local_vol.hpp` |
| `test_integration_phase3` | 10 | 端到端：组合标定→全 Greeks→VaR/ES→压力测试 |

### 1.3 必须达到的数值基准
| 基准 | 容差 | 验收方式 |
|------|------|----------|
| AAD Greeks vs 解析解 (BS/Heston) | 相对误差 < 1e-10 | `test_aad_greeks.cpp` |
| Pathwise Greeks vs AAD (路径相关) | 相对误差 < 1e-6 | 同一路径集对比 |
| LR Greeks vs AAD (不连续 payoff) | 相对误差 < 1e-4 | 数字期权/障碍期权 |
| FD Greeks (中心差分 h=1e-4) | 相对误差 < 1e-6 | 光滑 payoff |
| 历史 VaR (99%, 250日) 回测覆盖率 | 1% ± 0.3% | Kupiec POF 检验 p>0.05 |
| MC VaR (1M paths) 标准误差 | < 1% VaR 值 | Bootstrap 置信区间 |
| ES 相对误差 (MC vs Cornish-Fisher) | < 2% | 正态/学生 t 分布对比 |
| Heston 联合标定 (价格+IV) | 目标函数 < 1e-6 | 市场数据重现 |
| SVI 无套利参数化 | 所有行权/期限 Butterfly/Calendar 套利检测 0 | `test_svi_no_arb.cpp` |
| Dupire 局部波动率 vs 市场价格 | 重现误差 < 1bp | `test_dupire_recovery.cpp` |

---

## 2. Week 1: Greeks 体系与自动微分 (M1)

### 2.1 自动微分核心 (`include/cpphub/risk/greeks/`)

| 文件 | 关键类/方法 | 行数 | 核心算法 |
|------|-------------|------|----------|
| `ad_tape.hpp` | `autodiff::var` 适配层 (封装 autodiff 库) | ~80 | 逆模式 AD (autodiff 库, ADR-007 Revised) |
| `ad_dual.hpp` | `autodiff::dual` 适配层 (前向模式: val, der) | ~60 | 前向模式 (autodiff 库) |
| `aad_greeks.hpp` | `AADGreeksEngine` (一次扫描全 Greeks) | ~200 | 伴随模式: 记录 → 反向传播 |

**AAD 实现规范 (采用 autodiff 库, ADR-007 Revised)**:
```cpp
#include <autodiff/forward/dual.hpp>
#include <autodiff/reverse/var.hpp>
using namespace autodiff;

// 前向模式: 单个 Greek (autodiff::dual)
// 反向模式: 全 Greeks 一次扫描 (autodiff::var)
// 不自研 Tape: std::function 类型擦除阻止内联, 性能比成熟库慢 10-50x

// 示例: BS AAD Greeks
Greeks computeBSGreeksAAD(const VanillaOption& opt) {
    var S = opt.spot(), K = opt.strike(), r = opt.rate();
    var sigma = opt.vol(), T = opt.expiry();
    var d1 = (log(S/K) + (r + sigma*sigma*0.5)*T) / (sigma*sqrt(T));
    var d2 = d1 - sigma*sqrt(T);
    var price = S * cdf(d1) - K * exp(-r*T) * cdf(d2);
    autodiff::derive(price, autodiff::wrt(S, sigma, r, T));  // 一次扫描
    return { .delta = S.derivative(), .vega = sigma.derivative(), ... };
}
```

### 2.2 Greeks 计算引擎 (`include/cpphub/risk/greeks/`)

| 文件 | 关键类 | 行数 | 适用场景 |
|------|--------|------|----------|
| `analytic_greeks.hpp` | `AnalyticGreeks` (BS/Bachelier/Black76 闭式) | ~150 | 欧式期权、光滑 payoff |
| `fd_greeks.hpp` | `FDGreeks` (中心差分/复阶差分/自动步长选择) | ~150 | 通用兜底、验证用 |
| `pathwise_greeks.hpp` | `PathwiseGreeks` (路径法: d(payoff)/dθ) | ~200 | 欧式/障碍/亚式 (连续 payoff) |
| `likelihood_ratio_greeks.hpp` | `LRGreeks` (似然比法: payoff * dlog p/dθ) | ~200 | 数字期权/障碍期权 (不连续 payoff) |
| `aad_greeks.hpp` | `AADGreeks` (伴随模式: 一次扫描全阶) | ~250 | 复杂模型/篮子期权/高阶 Greeks |
| `greeks_factory.hpp` | `GreeksFactory` (自动选择最优方法) | ~100 | 统一入口 |

**Greeks 自动选择策略**:
```cpp
enum class GreeksMethod { Auto, Analytic, Pathwise, LR, FD, AAD };

class GreeksEngine {
public:
    Greeks compute(const Option& opt, const Model& model, GreeksMethod method = Auto) {
        if (method == Auto) {
            if (opt.is_european() && model.has_analytic_greeks()) return analytic_greeks(opt, model);
            if (opt.payoff_is_smooth() && model.supports_pathwise()) return pathwise_greeks(opt, model);
            if (!opt.payoff_is_smooth() && model.supports_lr()) return lr_greeks(opt, model);
            if (model.supports_aad()) return aad_greeks(opt, model);
            return fd_greeks(opt, model);
        }
        // 显式方法分发...
    }
};
```

### 2.3 高阶 Greeks (二阶/三阶)

| Greek | 符号 | 计算方法 | 用途 |
|-------|------|----------|------|
| Gamma | Γ = ∂²V/∂S² | AAD Dual2 / FD 中心差分 | Delta 对冲再平衡 |
| Vanna | ∂²V/∂S∂σ | AAD Dual2 | 波动率微笑风险 |
| Volga | ∂²V/∂σ² | AAD Dual2 | 波动率凸性风险 |
| Charm | ∂²V/∂S∂t | FD / AAD | 时间衰减对 Delta 影响 |
| Speed | ∂³V/∂S³ | FD 三阶差分 | Gamma 变化率 |
| Zomma | ∂³V/∂S²∂σ | AAD 三阶 | Gamma 对波动率敏感度 |

---

## 3. Week 2: VaR/ES 与场景分析 (M2)

### 3.1 VaR/ES 引擎 (`include/cpphub/risk/var/`)

| 文件 | 关键类/方法 | 行数 | 核心算法 |
|------|-------------|------|----------|
| `historical_var.hpp` | `HistoricalVaR` (滚动窗口、分位数插值) | ~120 | 非参数、全估计 |
| `parametric_var.hpp` | `ParametricVaR` (正态/学生t/峰度调整 Cornish-Fisher) | ~150 | 协方差矩阵 + 分位数展开 |
| `mc_var.hpp` | `MCVaR` (完整估值/Delta-Gamma 近似/Delta 近似) | ~200 | 全估值 MC + 方差缩减 |
| `expected_shortfall.hpp` | `ExpectedShortfall` (积分定义/分位数平均/MC 平均) | ~150 | ES = E[L | L > VaR] |
| `backtesting.hpp` | `KupiecPOF`, `ChristoffenIID`, `BaselTrafficLight` | ~150 | 监管回测标准 |

**MC VaR 三种近似层级**:
```cpp
enum class VaRApproximation { Full, DeltaGamma, Delta };

class MCVaR {
    VaRApproximation approx_;
    MCConfig mc_cfg_;
public:
    Real compute_var(const Portfolio& pf, Real confidence, Size horizon_days) {
        if (approx_ == Full) {
            // 完整重估: 每条路径重新定价所有工具
            return full_revaluation_var(pf, confidence, horizon_days);
        } else if (approx_ == DeltaGamma) {
            // Delta-Gamma 近似: V ≈ V0 + Δ'S + 0.5 S' Γ S
            return delta_gamma_var(pf, confidence, horizon_days);
        } else {
            // Delta 近似: V ≈ V0 + Δ'S
            return delta_var(pf, confidence, horizon_days);
        }
    }
};
```

### 3.2 场景分析与压力测试 (`include/cpphub/risk/scenario/`)

| 文件 | 关键类 | 行数 | 场景类型 |
|------|--------|------|----------|
| `stress_test.hpp` | `StressTest` (历史/假设/相关性冲击) | ~200 | 2008/2020/自定义 |
| `sensitivity.hpp` | `SensitivityAnalysis` (单因子/多因子/朝向分析) | ~150 | Greeks 为基础的敏感度 |

**压力测试场景模板**:
```cpp
struct StressScenario {
    std::string name;
    std::map<std::string, Real> spot_shocks;      // 股价/汇率/大宗商品冲击
    std::map<std::string, Real> rate_shocks;      // 利率曲线平行/倾斜/弯曲
    std::map<std::string, Real> vol_shocks;       // 隐含波动率面冲击
    std::map<std::pair<std::string,std::string>, Real> corr_shocks; // 相关性冲击
    Real probability;                             // 发生概率 (用于加权 ES)
};

class StressTester {
public:
    PortfolioPnL run(const Portfolio& pf, const StressScenario& scenario);
    std::vector<PortfolioPnL> run_historical(const Portfolio& pf, 
                                              const std::vector<HistoricalCrisis>& crises);
    // 监管要求: Basel III FRTB, CCAR, EBA 压力测试模板
};
```

---

## 4. Week 3: 标定框架/波动率曲面/验收 (M3)

### 4.1 标定框架 (`include/cpphub/calibration/`)

| 文件 | 关键类 | 行数 | 支持模型 |
|------|--------|------|----------|
| `optimizer.hpp` | `LevenbergMarquardt` (数值雅可比+阻尼), `NelderMead`, `DifferentialEvolution` | ~300 | 全局+局部混合 |
| `calibrator.hpp` | `Calibrator` 基类 + `HestonCalibrator`, `SABRCalibrator`, `SVICalibrator` | ~250 | 多目标: 价格+IV+Greeks |
| `objective.hpp` | `ObjectiveFunction` (VEGA 加权/价格加权/相对误差) | ~100 |  |

**标定流程**:
```cpp
// Heston 标定典型流程:
// 1. DE 全局搜索 (种群 50, 200 代) 得初始猜测
// 2. LM 局部精化 (数值雅可比, 阻尼因子自适应)
// 3. 目标函数: Σ w_i (IV_model - IV_market)² + λ * (price_model - price_market)²
// 4. 约束: Feller 条件 (2κθ > σ²), 相关性 ρ ∈ [-1,1], κ,θ,σ > 0
// 5. 诊断: 雅可比条件数、参数相关性、残差分布

class HestonCalibrator : public Calibrator {
    CalibrationResult calibrate(const MarketData& market, const CalibConfig& cfg) {
        // 1. DE 全局搜索
        auto de_result = differential_evolution(objective, bounds, cfg.de_pop_size, cfg.de_generations);
        // 2. LM 局部精化
        auto lm_result = levenberg_marquardt(objective, de_result.best_params, cfg.lm_max_iter);
        // 3. 诊断
        return {lm_result.params, lm_result.objective, lm_result.jacobian, lm_result.covariance};
    }
};
```

### 4.2 波动率曲面建模 (`include/cpphub/models/vol_surface/`)

| 文件 | 关键类/函数 | 行数 | 核心特性 |
|------|-------------|------|----------|
| `vol_surface.hpp` | `VolSurface` (插值: 双线性/样条/ SVI 参数化) | ~200 | 任意行权/期限查询 |
| `svi.hpp` | `SVI` (原始/自然/跳跃翼参数化) + 无套利约束 | ~200 | 无套利保证 |
| `ssvi.hpp` | `SSVI` (可分离/无套利跨期限) | ~150 | 期限结构一致性 |
| `dupire_local_vol.hpp` | `DupireLocalVol` (有限差分/样条求导) | ~150 | 局部波动率恢复 |

**SVI 无套利约束实现**:
```cpp
// SVI 原始参数化: w(k) = a + b * (ρ*(k-m) + sqrt((k-m)² + σ²))
// 无套利充要条件 (Roger Lee / Gatheral):
// 1. a > 0, b > 0, σ > 0, |ρ| < 1
// 2. Butterfly 套利: w(k) >= 0 且 g(k) = (1 - k*w'(k)/2w)² - w'(k)²/4 * (w + 0.25) + w''(k)/2 >= 0
// 3. Calendar 套利: w_T1(k) <= w_T2(k) 对于 T1 < T2 (方差单调性)

class SVISurface {
    bool check_no_arbitrage() const {
        // 检查所有行权价和期限的 Butterfly/Calendar 条件
        // 返回第一个违反约束的 (k, T) 对
    }
    
    // 校准时使用投影法或罚函数法强制满足约束
    void calibrate_with_constraints(const MarketQuotes& quotes);
};
```

**Dupire 局部波动率恢复**:
```cpp
// Dupire 公式: σ²_loc(K,T) = (∂C/∂T + qC + (r-q)K ∂C/∂K) / (0.5 K² ∂²C/∂K²)
// 实现: 使用三次样条拟合隐含波动率面, 解析求导
// 验证: 用恢复的局部波动率跑 PDE/MC, 对比市场价格, 误差 < 1bp
```

---

## 5. 参数估计与统计推断 (v1.1+ 延伸目标)

> **状态**: 规划中，不纳入 v1.1 核心交付。v1.1 完成后评估是否下沉到 v1.2。
> **背景**: 当前 §4 的"标定"是**衍生品标定** (用市场价格反推模型参数, 工程导向)。
> 本节规划**统计估计** (用历史数据拟合模型参数, 推断导向), 两者性质不同 (见 ADR-014)。

### 5.1 标定 vs 估计的区别

| 维度 | 衍生品标定 (§4, 已规划) | 统计估计 (本节, 延伸) |
|---|---|---|
| **目标** | 模型价格匹配市场价格 | 估计统计模型参数 |
| **方法** | LM + DE (数值优化) | MLE / QMLE / GMM / Bayesian |
| **数据** | 期权报价 (IV surface) | 时间序列 (收益、因子) |
| **典型问题** | Heston 参数 → 匹配 IV 面 | GARCH 参数 → 拟合收益波动率 |
| **标准误差** | 通常不计算 | **必须计算** (Newey-West / sandwich / Hessian) |
| **模型检验** | 残差分析 | LR / Wald / LM 检验 / 信息准则 |
| **现有 C++ 库** | QuantLib 有 | **完全没有** (见调研报告 §11) |

### 5.2 规划方法清单 (按优先级)

**P0 (v1.2 候选)**:

| 方法 | 用途 | 基准来源 | 矩阵需求 |
|---|---|---|---|
| MLE / QMLE | GARCH(1,1) 参数估计 | arch + rugarch | 动态向量化, 递归滤波 |
| Hessian 标准误差 | MLE 的渐近推断 | statsmodels | 动态矩阵求逆 |
| Sandwich 估计量 | MLE 稳健标准误差 | sandwich (R) | 动态 X'X + bread + meat |

**P1 (v1.3 候选)**:

| 方法 | 用途 | 基准来源 |
|---|---|---|
| GMM (Hansen 1982) | 矩估计, IV 估计 | linearmodels |
| Block Bootstrap | 时间序列重采样 | arch.bootstrap |
| Wald / LR / LM 检验 | 模型检验 | statsmodels |

**P2 (v2.0+ 延伸)**:

| 方法 | 用途 | 基准来源 |
|---|---|---|
| MCMC (Metropolis-Hastings) | 贝叶斯估计 | PyMC / Stan |
| 状态空间 + Kalman 滤波 | 时变参数估计 | statsmodels statespace |
| SMM (Simulated Method of Moments) | 结构模型估计 | - |

### 5.3 与已有模块的共享

- **optimizer.hpp** (§4.1): LM/DE 可复用于 MLE 优化
- **econometrics/newey_west.hpp** (v1.1 计量 P0): 标准误差计算可复用
- **linalg_dynamic.hpp** (ADR-013): 动态矩阵运算, 参数估计必需

### 5.4 目录规划 (v1.2 实施时)

```
include/cpphub/econometrics/estimation/
├── mle.hpp              # 极大似然估计 (QMLE 支持)
├── gmm.hpp              # 广义矩估计
├── bootstrap.hpp        # Block / Wild bootstrap
├── standard_errors.hpp  # Hessian / Sandwich / Newey-West
└── hypothesis_tests.hpp # Wald / LR / LM 检验
```

### 5.5 对 Research OS 的意义

统计估计模块是"因子失效诊断"方向的必需基础设施:
- GARCH MLE → 估计时变波动率, 用于因子收益的风险调整
- Newey-West HAC → 因子回归的稳健推断 (已规划在计量 P0)
- GMM → IV 估计, 处理因子内生性
- Bootstrap → Romano-Wolf 多重检验 (已规划在计量 P2)

没有统计估计, Cpp_Hub 只能"定价"和"标定", 不能"推断"。完整的量化研究需要三者皆备。

详见 [ADR-014: 标定 vs 估计的分离](../decisions/ADR_INDEX.md#adr-014) 和 [ECONOMETRICS_LANDSCAPE.md §11](../research/ECONOMETRICS_LANDSCAPE.md#11-参数估计-c-生态)。

---

## 6. 验收检查表 (Phase 3)

| 类别 | 检查项 | 通过标准 |
|------|--------|----------|
| **编译** | 增量编译 < 30s | `cmake --build build --target cpphub` |
| **单测** | 新增 150+ 测试全绿 | `ctest -R "phase3"` |
| **AAD Greeks** | vs 解析解 1e-10 | `test_aad_greeks.cpp` |
| **Pathwise Greeks** | vs AAD 1e-6 | 同一路径集对比 |
| **LR Greeks** | vs AAD 1e-4 | 数字期权/障碍期权 |
| **FD Greeks** | 中心差分 1e-6 | 光滑 payoff |
| **历史 VaR** | Kupiec POF p>0.05 | 99% 250日回测 |
| **MC VaR** | 标准误差 < 1% VaR | Bootstrap CI |
| **ES** | MC vs Cornish-Fisher < 2% | 正态/t 分布对比 |
| **Heston 联合标定** | 目标函数 < 1e-6 | 市场数据重现 |
| **SVI 无套利** | Butterfly/Calendar 0 违反 | `test_svi_no_arb.cpp` |
| **Dupire 恢复** | 重现误差 < 1bp | `test_dupire_recovery.cpp` |
| **性能** | PDE 400x2000 < 100ms | `benchmarks/pde_benchmark.cpp` |
| **Python** | 新增模型/引擎绑定可用 | `cpphub.HestonModel`, `cpphub.COSEngine` |

---

**Phase 3 负责人**: _______________  
**审核人**: _______________  
**开始日期**: _______________  
**预计结束**: _______________