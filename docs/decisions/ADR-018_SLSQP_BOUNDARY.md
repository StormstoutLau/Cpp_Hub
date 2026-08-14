# ADR-018: SLSQP 优化器实现边界

**状态**: Accepted
**日期**: 2026-08-15
**版本归属**: v1.6 (Phase 7B 阶段 1 前置)
**关联 Phase**: 7B
**决策者**: 架构组
**调研依据**: [SLSQP_EXTENSION_RESEARCH.md](../research/SLSQP_EXTENSION_RESEARCH.md) v1.0
**设计方案**: [SLSQP_EXTENSION_DESIGN.md](../research/SLSQP_EXTENSION_DESIGN.md) v1.0 §6
**前置 ADR**: [ADR-016](ADR_INDEX.md#adr-016-金融时间序列实施边界-18-项) (G-ADR2: GARCH 参数约束方法选 SLSQP)

---

## 背景

### 约束 C1 (v1.6 M1 前置)

Phase 7B 金融时间序列模块 (v1.6 M1 GARCH) 实施前置约束 C1: 现有 optimizer 仅支持 LM/NelderMead/DE 三类无约束优化器, 无法满足 GARCH 参数约束需求 (ω>0, α≥0, β≥0, α+β<1)。

### ADR-016 G-ADR2 决策

GARCH 参数约束方法选 **SLSQP** (序列二次规划), 依据:
- `arch` 包默认用 SLSQP (scipy.optimize.minimize method='SLSQP')
- `rugarch` 默认用 solnp (增广 Lagrangian)
- log 重参数化 (ω=e^W) 配合无约束优化是备选方案, 但改变似然函数地形, 不利于数值对照

### SLSQP 算法核心 (Kraft 1988 NLPQL)

- **QP 子问题**: 每次迭代求解二次规划, 得到搜索方向 d
- **BFGS Hessian 近似**: 用 Lagrangian 梯度 (非目标梯度) 更新
- **L1 merit function**: f(x) + Σ μ_i|c_eq| + Σ μ_i max(0, -c_ineq), 用于线搜索
- **Armijo 线搜索**: 基于 merit function 回溯步长

---

## 决策

### 7 项实现边界决策

| # | 决策点 | 候选方案 | 决策 | 依据 |
|---|--------|---------|------|------|
| 1 | QP 求解策略 | (a) KKT+slack (b) active-set (c) 投影梯度 | **(b) active-set** | 数学严谨, 与 scipy 一致; 小问题 (n<10) 收敛快; KKT 系统求解等式约束 QP |
| 2 | BFGS 更新 | (a) 标准 (b) Damped (c) L-BFGS | **(b) Damped** | 曲率条件 sy >= 0.2*sBs 不满足时用 damped 版本 (Nocedal-Wright Eq 18.15), 保证 B 正定 |
| 3 | 约束表达 | (a) c_i(x) <= 0 (b) c_i(x) >= 0 | **(b) c_i(x) >= 0** | 与 scipy `f_ieqcons` 和 arch `a.dot(x)-b` 一致 (调研 I1 已确认) |
| 4 | 数值差分 | (a) 前向 (b) 中心 | **(b) 中心** | O(h²) 精度, 比 scipy 默认前向差分 O(h) 更精确; 步长 h = eps*(1+|x_j|) 自适应 |
| 5 | 边界处理 | (a) 转不等式 (b) QP 内直接处理 | **(b) QP 内直接处理** | 边界约束作为不等式加入 active-set (d >= l-x, d <= u-x), 避免 2n 个冗余约束 |
| 6 | Eigen3 依赖 | (a) 引入 (b) 不引入 | **(b) 不引入** | ADR-013 隔离原则; n<10 用 std::vector + Gauss-Jordan 足够; BFGS 矩阵用 std::vector<std::vector<Real>> |
| 7 | 等式约束支持 | (a) 不支持 (b) 支持 | **(b) 支持** | 通用性, 未来 ARIMA/MIDAS 可能需要等式约束 (如固定持久性 α+β=0.95) |

### 接口设计

```cpp
using ConstraintFn = std::function<std::vector<Real>(const std::vector<Real>&)>;
using ConstraintJacobianFn = std::function<std::vector<std::vector<Real>>(const std::vector<Real>&)>;

class SLSQP {
public:
    struct Config {
        Size max_iterations = 100;
        Real ftol = 1e-6;
        Real xtol = 1e-8;
        Real gtol = 1e-6;
        Real armijo_gamma = 0.1;
        Real armijo_beta = 0.5;
        Real epsilon = 1e-6;
        Size max_line_search = 20;
        GradientFn gradient;                  // 可选, None 时数值差分
        ConstraintJacobianFn constraint_jacobian;  // 可选, None 时数值差分
    };

    static OptimizationResult minimize(
        const ObjectiveFn& f,
        const std::vector<Real>& x0,
        const std::vector<Bounds>& bounds,
        const std::vector<ConstraintFn>& ineq_constraints = {},  // c_i(x) >= 0
        const std::vector<ConstraintFn>& eq_constraints = {},    // c_i(x) = 0
        const Config& cfg = Config{});
};
```

---

## 理由

### 为什么选 active-set 而非 KKT+slack

- **数学严谨**: active-set 是标准 QP 求解方法 (Nocedal-Wright §16.4), 正确处理不等式约束
- **小问题高效**: GARCH 应用 (n=3-5, m=4-5) active-set 通常 2-5 次迭代收敛
- **与 scipy 一致**: scipy SLSQP Fortran 实现用 active-set
- **KKT+slack 缺点**: slack 变量增加维度, 非负性仍需 active-set 处理

### 为什么选 Damped BFGS

- **曲率条件**: 标准 BFGS 要求 y^T*s > 0, 但 SLSQP 中用 Lagrangian 梯度, 可能不满足
- **Damped 保证正定**: theta = 0.8*sBs/(sBs-sy) 保证 r^T*s > 0, B 保持正定
- **标准技术**: Nocedal-Wright Eq 18.15, 工业级 SQP 实现标配

### 为什么不等式用 c_i(x) >= 0

- **scipy 一致**: `f_ieqcons(x) >= 0` (调研 I1 已确认)
- **arch 一致**: `a.dot(params) - b >= 0` (调研 §4.1 已确认)
- **直觉清晰**: "约束满足" = "约束值非负"

### 为什么不引入 Eigen3

- **ADR-013 隔离原则**: cpphub_core 不链接 Eigen3
- **规模足够**: GARCH n<10, std::vector<std::vector<Real>> 足够
- **复用基础设施**: `detail::solve_linear_system` (Gauss-Jordan) 已支持

---

## 后果

### 正面

- GARCH 参数约束需求 (C1) 满足, Phase 7B M1 可启动
- 与 scipy/arch 接口一致, 便于数值对照
- 不引入 Eigen3, 保持 core 隔离
- 支持等式约束, 未来 ARIMA/MIDAS 可复用

### 负面

- active-set QP 实现复杂 (~150 行), 但封装在 detail 命名空间
- 完整 BFGS 矩阵存储 O(n²), 大问题 (n>50) 内存开销大 (v1.7+ 可扩展 L-BFGS)
- 数值差分比解析梯度慢 (2n 次函数调用), 但 GARCH 似然函数解析梯度复杂, 数值差分可接受

### 风险缓解

| 风险 | 缓解 |
|------|------|
| QP active-set 不收敛 | 最大迭代 50 次, fallback 梯度投影 |
| BFGS 非正定 | Damped BFGS + 曲率条件检查 |
| 线搜索失败 | 最大回缩 20 次, 返回小步长 |
| KKT 系统奇异 | `solve_linear_system` 返回 false 时用单位阵 |

---

## Scope 边界

### v1.6 M1 scope (本次实现)

- SLSQP 基础功能: bounds + 不等式约束 + 等式约束
- 满足 GARCH 参数约束需求

### 不在 scope

- **L-BFGS** (大规模, v1.7+): 完整 BFGS 矩阵存储改为 limited-memory
- **trust-constr** (更稳健, v1.7+): trust-region SQP 变体
- **多起始点** (G-ADR6, Phase 7B M1 scope): 用户手动调用多次
- **warm start** (预热, v1.7+): Config 预留扩展字段

---

## 关联

- **前置**: [ADR-016](ADR_INDEX.md#adr-016-金融时间序列实施边界-18-项) G-ADR2 (GARCH 参数约束方法)
- **依赖**: [ADR-013](ADR_INDEX.md#adr-013-双层线性代数架构-固定尺寸--动态尺寸) (Eigen3 隔离原则)
- **调研**: [SLSQP_EXTENSION_RESEARCH.md](../research/SLSQP_EXTENSION_RESEARCH.md) v1.0
- **设计**: [SLSQP_EXTENSION_DESIGN.md](../research/SLSQP_EXTENSION_DESIGN.md) v1.0
- **实施**: [SLSQP_EXTENSION_IMPLEMENTATION.md](../research/SLSQP_EXTENSION_IMPLEMENTATION.md) v1.0
- **验收**: [SLSQP_EXTENSION_ACCEPTANCE_CHECKLIST.md](../research/SLSQP_EXTENSION_ACCEPTANCE_CHECKLIST.md) v1.0

---

## 参考文献

1. Kraft, D. (1988). "A Software Package for Sequential Quadratic Programming". TOIMS 733.
2. Nocedal, J. & Wright, S. (2006). Numerical Optimization, 2nd ed. Springer. Ch 18 (SQP), §16.4 (Active-Set).
3. scipy.optimize.slsqp 源码: https://github.com/scipy/scipy/blob/main/scipy/optimize/_slsqp_py.py
4. arch package 源码: https://github.com/bashtage/arch/blob/main/arch/univariate/base.py

