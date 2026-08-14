# ADR-017: 时序模块命名空间 (`cpphub::v1::timeseries`)

**状态**: Accepted
**日期**: 2026-08-15
**版本归属**: v1.6 (Phase 7B)
**关联 Phase**: 7B
**决策者**: 架构组
**调研依据**: [FINANCIAL_TIMESERIES_RESEARCH.md](../research/FINANCIAL_TIMESERIES_RESEARCH.md) v3.2 §22 (兼容性约束 C4)
**前置 ADR**: [ADR-015](ADR_INDEX.md#adr-015-证伪统计量模块边界-通用-vs-模块特定) (命名空间不对称问题)

---

## 背景

### 现有命名空间不对称

Cpp_Hub v1.0-v1.5 的命名空间布局存在不对称:

| 模块 | 命名空间 | 子命名空间 | 物理目录 |
|------|---------|-----------|----------|
| core | `cpphub::v1` | 无 | `include/cpphub/core/` |
| instruments | `cpphub::v1` | 无 | `include/cpphub/instruments/` |
| models | `cpphub::v1` | 无 | `include/cpphub/models/` |
| pricing | `cpphub::v1` | 无 | `include/cpphub/pricing/` |
| risk | `cpphub::v1` | 无 | `include/cpphub/risk/` |
| econometrics | `cpphub::v1` | `::econometrics` | `include/cpphub/econometrics/` |
| hfecon | `cpphub::v1` | `::hfecon` | `include/cpphub/hfecon/` |

**问题**: risk/pricing 落在 `cpphub::v1` (无子命名空间), econometrics/hfecon 落在子命名空间。ADR-015 决策点 4 明确指出此不对称性, 但 Phase 7A scope 不涉及命名空间调整, 留待独立 ADR 决策。

### Phase 7B 新增时序模块

v1.6+ 将新增金融时间序列模块, 包含:
- GARCH 族 (GARCH/EGARCH/GJR-GARCH/APARCH/FIGARCH/IGARCH/GARCH-M)
- 单位根检验 (ADF/DF-GLS/KPSS/PP/Ng-Perron/方差比)
- ARIMA 族 (v1.6 M3)
- MIDAS (v1.6 M4)
- VAR 族 (v1.7)
- 多元波动率 (CCC/DCC, v1.7)
- 因果检验 (Granger/Toda-Yamamoto/Transfer Entropy, v1.7+)
- 长记忆/非线性 (ARFIMA/FIGARCH/MS-AR, v1.8+)

这些模块需要明确的命名空间归属。

---

## 决策

**时序模型统一落 `cpphub::v1::timeseries` 子命名空间**, 与 `econometrics`/`hfecon` 平级。

```cpp
namespace cpphub {
inline namespace v1 {
namespace timeseries {

// GARCH 族 (v1.6 M1)
namespace garch {
    class GarchModel;
    class EGarchModel;
    class GjrGarchModel;
    // ...
}

// 单位根检验 (v1.6 M2)
namespace unit_root {
    class ADFTest;
    class DFGlsTest;
    class KPSSTest;
    class PPTest;
    // ...
}

// ARIMA 族 (v1.6 M3)
namespace arima {
    class ARModel;
    class ARIMAModel;
    // ...
}

// MIDAS (v1.6 M4)
namespace midas {
    class MIDASModel;
    // ...
}

// VAR 族 (v1.7)
namespace var {
    class VARModel;
    class VECMModel;
    // ...
}

// 多元波动率 (v1.7)
namespace multivariate_vol {
    class CCCModel;
    class DCCModel;
    // ...
}

// 因果检验 (v1.7+)
namespace causality {
    class GrangerTest;
    class TodaYamamotoTest;
    class TransferEntropy;
    // ...
}

}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
```

### 物理目录结构

```
include/cpphub/timeseries/
├── garch/                    # v1.6 M1
│   ├── garch_model.hpp
│   ├── egarch_model.hpp
│   ├── gjr_garch_model.hpp
│   ├── aparch_model.hpp
│   ├── figarch_model.hpp
│   ├── igarch_model.hpp
│   └── garch_m_model.hpp
├── unit_root/                # v1.6 M2
│   ├── adf_test.hpp
│   ├── df_gls_test.hpp
│   ├── kpss_test.hpp
│   ├── pp_test.hpp
│   ├── ng_perron_test.hpp
│   └── variance_ratio_test.hpp
├── arima/                    # v1.6 M3
│   ├── ar_model.hpp
│   └── arima_model.hpp
├── midas/                    # v1.6 M4
│   └── midas_model.hpp
├── var/                      # v1.7
│   ├── var_model.hpp
│   └── vecm_model.hpp
├── multivariate_vol/         # v1.7
│   ├── ccc_model.hpp
│   └── dcc_model.hpp
└── causality/                # v1.7+
    ├── granger_test.hpp
    ├── toda_yamamoto_test.hpp
    └── transfer_entropy.hpp
```

---

## 理由

### 1. 与 econometrics/hfecon 平级, 保持一致性

econometrics 和 hfecon 都有独立子命名空间, 时序模块作为独立的方法论体系 (GARCH/单位根/VAR/因果), 应享有同等待遇。`cpphub::v1::timeseries` 与 `cpphub::v1::econometrics`、`cpphub::v1::hfecon` 平级, 结构清晰。

### 2. 时序模块的方法论独立性

金融时间序列与横截面计量有本质差异:
- **波动率本身就是研究对象** (GARCH 族), 不是残差诊断
- **时间维度是核心** (单位根/协整/长记忆), 不是横截面维度
- **动态系统建模** (VAR/IRF/FEVD), 不是静态关系
- **信息论度量** (Transfer Entropy), 跨学科方法

独立命名空间反映这种方法论独立性。

### 3. 避免 `cpphub::v1` 顶层污染

若时序模型直接落 `cpphub::v1` (如 risk/pricing), 会导致:
- `cpphub::v1` 顶层符号爆炸 (GARCH 族 7+ 类, 单位根 6+ 类, VAR 族 5+ 类, 因果 10+ 类)
- 与 risk/pricing 的符号混淆 (如 `VARModel` vs risk 模块的 `VaR` 风险度量)
- 用户难以定位时序相关 API

### 4. 渐进式引入, 不破坏现有代码

新增 `timeseries` 子命名空间不影响现有 risk/pricing/econometrics/hfecon 代码。现有模块的时序相关功能 (如 `har_model.hpp` 在 hfecon, `volatility_diagnostics.hpp` 在 econometrics/inference) 保持原位, 不迁移。

### 5. 为 v1.7+ 多元时序预留空间

v1.7 将引入 VAR/DCC/协整等多元时序方法, v1.8+ 引入长记忆/非线性/MS-AR。`timeseries` 子命名空间为这些扩展提供清晰的归属, 避免未来再次决策命名空间。

---

## 替代方案评估

| 方案 | 优点 | 缺点 | 结论 |
|------|------|------|------|
| A. 落 `cpphub::v1` (无子命名空间, 如 risk/pricing) | 与 risk/pricing 一致 | 顶层符号爆炸, 与 VaR 混淆 | ❌ 拒绝 |
| **B. 落 `cpphub::v1::timeseries`** (本方案) | 与 econometrics/hfecon 平级, 方法论独立 | 新增子命名空间 | ✅ 采纳 |
| C. 落 `cpphub::v1::econometrics::timeseries` | 复用 econometrics 命名空间 | 时序 ≠ 计量 (波动率/因果是独立方法论), 层级过深 | ❌ 拒绝 |
| D. 统一所有模块到子命名空间 (含 risk/pricing) | 完全对称 | 需大规模迁移现有代码, 超出 v1.6 scope | ❌ 拒绝 (留待 v2.0) |

---

## 后果

### 实施约束

1. **新模块归属**: v1.6+ 所有时序相关头文件落 `include/cpphub/timeseries/`, 命名空间 `cpphub::v1::timeseries`
2. **现有模块不迁移**: hfecon 的 `har_model.hpp`、econometrics 的 `volatility_diagnostics.hpp` 保持原位, 不迁移到 timeseries
3. **CMake 目标**: 新增 `cpphub_timeseries` INTERFACE 库 (header-only), 或合并到 `cpphub_econometrics` (若依赖 Eigen3)
4. **Eigen3 隔离**: GARCH/单位根/ARIMA 的参数估计可能需 Eigen3 (矩阵运算), 按 ADR-013 原则放 `cpphub_econometrics` 链接 Eigen3; 纯序列检验 (如方差比) 可放 `cpphub_core` 不依赖 Eigen3

### 命名空间不对称的后续处理

ADR-015 决策点 4 指出的 risk/pricing 无子命名空间问题, 本 ADR **不解决** (超出 v1.6 scope)。留待 v2.0 独立 ADR 统一决策:
- 选项 1: risk/pricing 新增子命名空间 (`cpphub::v1::risk`, `cpphub::v1::pricing`)
- 选项 2: 所有模块统一到子命名空间
- 选项 3: 维持现状 (混合模式)

v1.6-v1.8 期间采用 `using` 声明过渡, 保证现有代码兼容。

### C ABI 版本化

时序模块的 C ABI 符号用 `cpphub_v1_6_*` (v1.6) / `cpphub_v1_7_*` (v1.7) 版本前缀 (约束 C5), 与 ADR-009 版本化原则一致。

---

## 关联

- **调研报告**: [FINANCIAL_TIMESERIES_RESEARCH.md](../research/FINANCIAL_TIMESERIES_RESEARCH.md) v3.2 §22 (约束 C4)
- **前置 ADR**: [ADR-013](ADR_INDEX.md#adr-013-双层线性代数架构-固定尺寸--动态尺寸) (Eigen3 隔离), [ADR-015](ADR_INDEX.md#adr-015-证伪统计量模块边界-通用-vs-模块特定) (命名空间不对称问题, 决策点 4)
- **关联 ADR**: [ADR-016](ADR-016_FINANCIAL_TIMESERIES_BOUNDARY.md) (金融时间序列实施边界, 18 项决策)
- **后续 ADR**: 待编写 ADR-018 (risk/pricing 子命名空间统一, v2.0 scope)
