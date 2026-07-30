# Phase 2 执行规格书 - 进阶模型与数值方法

> **版本归属**: **v1.0 核心** (削减后，与 Phase 1 共同构成 v1.0)
> **目标**: 3 周内交付 Heston 模型、PDE/树形引擎、准蒙特卡洛完整体系
> **前置**: Phase 1 全部通过
> **里程碑**: M1(Week 1) Heston 模型层 → M2(Week 2) PDE/树形 → M3(Week 3) Sobol/QMC/方差缩减/验收
>
> **v1.0 Scope 削减声明** (2026-07-29 评审):
> - ✅ **v1.0 保留**: Heston 模型、PDE (CN+PSOR)、树形 (CRR/LR)、Sobol QMC、基础方差缩减 (Antithetic/Control Variate/Moment Matching)
> - ⏸️ **v1.1 推迟**: SABR/Bates/VG/CEV 模型、利率模型 (Vasicek/CIR/HW/G2++)、COS/FFT/CONV 引擎、LSMC、SVI/SSVI/Dupire 波动率曲面、标定框架 (LM/DE)
> - ⏸️ **v2.0 砍掉**: 无 (本 Phase 无 v2.0 内容)
>
> 削减理由: 12 周单人开发无法覆盖原 scope。v1.0 聚焦 BS/Heston + Analytic/MC/PDE/Tree + Sobol QMC，足以覆盖主流欧式/美式/路径相关期权定价需求。进阶模型 (SABR/Levy/IR)、傅里叶引擎、LSMC、标定框架推迟到 v1.1。

---

## 1. 交付物清单

### 1.1 新增编译目标
| 目标 | 类型 | 关键源文件 |
|------|------|------------|
| `cpphub` (增量) | Shared Library | 新增 `src/models/heston.cpp`, `src/pricing/{pde,tree}/*.cpp` |

### 1.2 必须通过的新增测试 (v1.0 削减后)
| 测试套件 | 用例数 | 覆盖模块 |
|----------|--------|----------|
| `test_heston_process` | 15 | `heston.hpp` (特征函数、Euler/QE/Exact 模拟) |
| `test_pde_engine` | 20 | `pde_engine.hpp` (CN/PSOR/边界/非均匀网格) |
| `test_tree_engine` | 15 | `tree_engine.hpp` (CRR/JR/Tian/LR/三叉树) |
| `test_sobol_qmc` | 12 | `sobol.hpp` (Sobol 序列、Brownian Bridge、方差缩减) |
| `test_variance_reduction` | 10 | `control_variate.hpp`, `moment_matching.hpp` |
| `test_integration_phase2` | 6 | 端到端：Heston MC → PDE 美式 → Tree 验证 |

> **v1.1 推迟的测试**: `test_sabr_process`, `test_bates_process`, `test_vg_process`, `test_ir_models`, `test_cos_engine`, `test_fft_engine`, `test_lsmc_engine`, `test_calibration`

### 1.3 必须达到的数值基准 (v1.0 削减后)
| 基准 | 容差 | 参考值来源 |
|------|------|------------|
| Heston MC (1M paths, QE) vs 解析特征函数定价 | 相对误差 < 0.5% | Heston (1993) 闭式解 |
| Heston 特征函数 vs Schoutens 表 | 1e-10 | 复数值对比 |
| 美式 Put (PDE CN+PSOR) vs Broadie-Detemple | 1e-4 | 基准文献 |
| 美式 Put (Leisen-Reimer N=200) vs 基准 | 1e-4 | 高阶树收敛 |
| 亚式几何平均 MC vs 解析解 | 1e-6 | 控制变量验证 |
| Sobol QMC 方差缩减 | > 10x vs 伪随机 | 同路径数对比 |

> **v1.1 推迟的基准**: Heston COS (N=256) 1e-8, Bates COS, VG COS, SABR Hagan IV, Heston 标定, SABR 标定

---

## 2. Week 1: Heston 模型层 (M1)

### 2.1 扩散过程模型 (`include/cpphub/models/diffusion/`)

| 文件 | 关键类/方法 | 行数 | 核心算法 |
|------|-------------|------|----------|
| `process.hpp` | `StochasticProcess` (Template Method: `dimension()`, `generate_path()`, `characteristic_function()`) | ~80 | 基类契约 |
| `gbm.hpp` | `GBM` (精确解: `S*exp((r-q-0.5σ²)t + σW)`) | ~60 | Phase 1 已有 |
| `heston.hpp` | `Heston` (v0, κ, θ, σ, ρ) | ~250 | **核心**: 特征函数 `cf(u,τ)` + QE/Exact/FullTruncation 模拟 |
| `heston_qe.hpp` | `HestonQE` (Quadratic Exponential 方案) | ~150 | Andersen (2008) QE-M/FullTruncation |

> **v1.1 推迟**: `bates.hpp`, `sabr.hpp`, `cev.hpp`, `variance_gamma.hpp`

**Heston 特征函数实现规范 (必须与 Fang-Oosterlee 一致)**:
```cpp
// pricing/analytic/heston_cf.hpp
std::complex<Real> heston_characteristic_function(
    std::complex<Real> u, Real tau,
    Real v0, Real kappa, Real theta, Real sigma, Real rho, Real r, Real q) {
    // Heston (1993) 闭式特征函数
    // 使用 Kahl-Jäckel "Rotation Count" 算法避免分支切割
    // 返回 log 特征函数 log φ(u) = C(τ,u) + D(τ,u)*v0 + i*u*(log S0 + (r-q)τ)
    // 必须通过 test_heston_cf_accuracy: 对比论文 Table 1 复数值
}
```

**Heston 路径生成三种方案 (Strategy 模式)**:
```cpp
enum class HestonScheme { Euler, QE_M, QE_Q, Exact };

class HestonPathGenerator {
    HestonScheme scheme_;
    void generate_euler(...);           // 简单但需反射/截断
    void generate_qe(Real v0, ...);     // Andersen (2008) QE 方案
    void generate_exact(Real v0, ...);  // Broadie-Kaya 精确模拟 (非中心卡方)
};
```

### 2.2 利率模型 (`include/cpphub/models/ir/`)

> ⏸️ **v1.1 推迟**: `short_rate.hpp`, `vasicek.hpp`, `cir.hpp`, `hullwhite.hpp`, `g2pp.hpp`
> v1.0 假设定常无风险利率 r，不实现利率曲线模型。

### 2.3 波动率曲面与局部波动率 (`include/cpphub/models/vol_surface/`)

> ⏸️ **v1.1 推迟**: `vol_surface.hpp`, `svi.hpp`, `sabr_calibration.hpp`, `dupire_local_vol.hpp`
> v1.0 使用单一波动率参数 σ，不实现波动率曲面。

---

## 3. Week 2: PDE/树形引擎 (M2)

### 3.1 PDE 引擎完善 (`include/cpphub/pricing/pde/`)

| 文件 | 关键类/方法 | 行数 | 关键算法 |
|------|-------------|------|----------|
| `fdm_grid.hpp` | `FDMGrid` (非均匀网格: `sinh` 变换集中 ATM, 时间网格自适应) | ~150 | `build_nonuniform_grid(S0, K, σ, T)` |
| `fdm_scheme.hpp` | `FDMScheme` 基类 + `ExplicitEuler`, `ImplicitEuler`, `CrankNicolson`, `DouglasADI`, `CraigSneydADI` | ~250 | ADI 分裂误差控制 |
| `boundary.hpp` | `BoundaryCondition` (Dirichlet/Neumann/线性互补 PSOR) | ~150 | 美式早期行使: `PSORSolver` |
| `thomas_algorithm.hpp` | `thomas_algorithm` (三对角 O(N)), `block_thomas` (块三对角) | ~100 | Phase 1 已有 |
| `pde_engine.hpp` | `PDEEngine : PricingEngine` (price, greeks, price_american) | ~200 | 统一接口 |

**Crank-Nicolson + PSOR 美式期权实现规范**:
```cpp
// 每个时间步求解线性互补问题 (LCP):
// max( A * V^{n+1} - b, V^{n+1} - payoff ) = 0
// 使用 Projected SOR (PSOR):
void psor_step(Vector& V, const Vector& payoff, const Matrix& A, const Vector& b,
               Real omega, Real tol, int max_iter) {
    // V = max(V - omega * (A*V - b), payoff) 逐分量迭代
    // 收敛判据: ||V_new - V_old||_inf < tol
    // omega 最优约 1.2-1.5 (需实验调优)
}
```

### 3.2 树形引擎 (`include/cpphub/pricing/tree/`)

| 文件 | 关键类 | 行数 | 支持类型 |
|------|--------|------|----------|
| `binomial.hpp` | `BinomialTreeEngine` (CRR/JR/Tian/Leisen-Reimer) | ~250 | 欧式/美式/百慕大/期股息 |
| `trinomial.hpp` | `TrinomialTreeEngine` (显式/隐式/显隐混合) | ~200 | 障碍/亚式/美式 |
| `tree_engine.hpp` | `TreeEngine : PricingEngine` 统一接口 | ~100 |  |

**Leisen-Reimer 高阶收敛实现**:
```cpp
// 核心: 选择概率 p 使得第三阶矩匹配
// h = (σ² + (r-q)²) * Δt;  n 取奇数
// d1 = (log(S/K) + (r-q+0.5σ²)T) / (σ√T);  d2 = d1 - σ√T
// pbar = 0.5 + 0.5 * sign(d1) * sqrt(1 - exp(-d1²/(n+1/3+0.1/(n+1))))
// p = 0.5 + 0.5 * sign(d2) * sqrt(1 - exp(-d2²/(n+1/3+0.1/(n+1))))
// 收敛阶 O(1/n²) 优于 CRR 的 O(1/n)
```

### 3.3 傅里叶引擎 (`include/cpphub/pricing/fourier/`)

> ⏸️ **v1.1 推迟**: `cos_method.hpp`, `fft_engine.hpp`, `conv_engine.hpp`, `fourier_engine.hpp`
> v1.0 不实现 COS/FFT/CONV 引擎。Heston 定价通过 MC (QE 方案) 或 PDE 完成。
> v1.1 将引入 COS 方法作为 Heston 的高精度解析引擎 (指数收敛)。

---

## 4. Week 3: Sobol QMC/方差缩减/验收 (M3)

### 4.1 LSMC (Longstaff-Schwartz) 美式/百慕大引擎

> ⏸️ **v1.1 推迟**: `lsmc_engine.hpp`, `basis_functions.hpp`, `regression.hpp`
> v1.0 美式期权通过 PDE (CN+PSOR) 或树形 (Leisen-Reimer) 定价，不实现 LSMC。
> v1.1 将引入 LSMC 处理高维美式/百慕大期权 (多因子模型)。

### 4.2 方差缩减 (v1.0 保留基础方法)

| 文件 | 关键类 | 行数 | 算法核心 | 版本 |
|------|--------|------|----------|------|
| `control_variate.hpp` | `ControlVariate` (BS 作为控制变量) | ~100 | Cov(S, Y) 优化系数 | ✅ v1.0 |
| `moment_matching.hpp` | `MomentMatchingVR` (样本矩匹配理论矩) | ~100 | 一阶/二阶矩校正 | ✅ v1.0 |
| `sobol.hpp` | `SobolSequence` (Joe-Kuo 数字网) + `BrownianBridge` | ~200 | QMC + BB 降维 | ✅ v1.0 |
| `importance_sampling.hpp` | `ImportanceSamplingVR` (Girsanov 变换漂移) | ~150 | 重要性采样 + 似然比权重 | ⏸️ v1.1 |
| `stratified_sampling.hpp` | `StratifiedSamplingVR` (分层 + 比例/最优分配) | ~120 | 分层方差最小化 | ⏸️ v1.1 |
| `conditional_mc.hpp` | `ConditionalMC` (解析积分条件期望) | ~100 | 亚式几何平均、障碍期权 BGK 修正 | ⏸️ v1.1 |

**Sobol QMC + Brownian Bridge 实现规范**:
```cpp
// Sobol 序列生成 (Joe-Kuo 数字网, 注: 数字网必须与原论文一致以保证低差异性质)
class SobolSequence {
    // 基础2激进数, 方向向量
    // generate(n, dim) -> Point in [0,1]^dim
    // 必须通过 test_sobol_property: 第一个点为 0, 数字均匀性
};

// Brownian Bridge 构造 (降低有效维度)
class BrownianBridge {
    // 将 Sobol 点映射为 Brownian 路径
    // 核心: 先生成终点, 再递归填充中间点
    // 对 Heston 等路径依赖期权, BB 可将方差缩减 10-100x
};
```

### 4.3 标定框架

> ⏸️ **v1.1 推迟**: `optimizer.hpp`, `calibrator.hpp`, `objective.hpp`
> v1.0 不实现标定框架。模型参数通过手动设置或外部输入。
> v1.1 将引入 LM + DE 标定框架，支持 Heston/SABR/SVI 标定。

---

## 5. 验收检查表 (Phase 2 - v1.0 削减后)

| 类别 | 检查项 | 通过标准 |
|------|--------|----------|
| **编译** | 增量编译 < 30s | `cmake --build build --target cpphub` |
| **单测** | 新增 78+ 测试全绿 | `ctest -R "phase2"` |
| **Heston MC** | QE 方案 vs 解析解 < 0.5% | 相同参数 1M paths |
| **Heston 特征函数** | vs Schoutens 表 1e-10 | 复数值对比 |
| **PDE 美式** | Broadie-Detemple 基准 1e-4 | `tests/validation/american_pde_broadie.cpp` |
| **树形收敛** | Leisen-Reimer O(1/n²) | log-log 斜率 -2.0 ± 0.1 |
| **Sobol QMC** | 方差缩减 > 10x vs 伪随机 | 同路径数对比 |
| **控制变量** | BS 控制变量方差缩减 > 5x | 亚式算术平均 MC |
| **性能** | PDE 400x2000 < 100ms | `benchmarks/pde_benchmark.cpp` |
| **Python** | Heston/PDE/Tree 绑定可用 | `cpphub.HestonModel`, `cpphub.PDEEngine` |

> **v1.1 推迟的验收项**: Heston COS Table 1, COS 指数收敛, LSMC 早期行使溢价, Heston/SABR 标定收敛

---

**Phase 2 负责人**: _______________
**审核人**: _______________
**开始日期**: _______________
**预计结束**: _______________