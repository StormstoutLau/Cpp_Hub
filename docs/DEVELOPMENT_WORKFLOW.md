# 开发工作流：LLM 辅助开发 + 基准对齐 + 形式化验证

> **版本**: 1.0
> **日期**: 2026-07-29
> **定位**: 本文件公开记录 Cpp_Hub 的开发方法论。我们不掩饰 LLM 辅助开发，而是将其作为正向能力信号——**驾驭 LLM 生成可验证的生产级代码**本身是稀缺技能。

---

## 1. 核心立场

### 1.1 工作流定位

Cpp_Hub 采用 **AI-assisted development** 工作流，三阶段闭环：

```
LLM 生成初稿 → 严格测试 + 基准对齐 → (可选) Lean4 形式化验证
     ↓                    ↓                      ↓
  加速交付              可验证的正确性         数学证明
```

### 1.2 为什么不掩饰

| 旧叙事（掩饰） | 新叙事（坦诚） |
|---|---|
| "我独立写了这些代码" | "我用 LLM 加速生成，但每个模块通过 1e-10 基准对齐" |
| 信号：C++ 编码能力 | 信号：**AI Coding 熟练度 + 数值验证能力 + 工程判断力** |
| 风险：被发现后信任崩塌 | 风险：无（测试覆盖率是客观证据） |

### 1.3 能力信号分解

这个工作流证明的不是"我能写 C++"，而是三个稀缺组合：

1. **架构能力**：ADR + Phase Spec + 模块化设计（LLM 做不到的）
2. **验证能力**：基准对齐 + 三源交叉验证 + 1e-10 容差（LLM 做不到的）
3. **AI Coding 熟练度**：用 LLM 加速交付（传统 C++ 工程师不具备的）

---

## 2. 三阶段工作流详解

### 2.1 阶段 1: LLM 生成初稿

**输入**：
- 架构设计文档（ADR + ARCHITECTURE.md）
- Phase Spec（明确接口契约、数值基准、测试要求）
- 溯源标注（教材页码、论文公式、开源库参考实现）

**LLM 任务**：
- 根据 spec 生成 `.hpp` / `.cpp` 初稿
- 生成测试用例骨架
- 生成文档注释

**人工任务**：
- 设计模块边界、接口契约
- 选择算法（如 Heston 用 QE 而非 Euler）
- 指定基准来源（如"对比 Fang-Oosterlee 2008 Table 1"）
- Review LLM 输出的架构一致性

**产出**：可编译但未验证的初稿代码 + 测试骨架

### 2.2 阶段 2: 严格测试 + 基准对齐（核心）

**这是整个工作流的护城河**。LLM 生成的代码不信任，必须通过客观基准验证。

#### 2.2.1 三源交叉验证策略

每个数值模块必须通过三源（或以上）交叉验证：

| 验证源 | 角色 | 容差 | 工具 |
|---|---|---|---|
| **Python 开源库** | 算法基准（可读源码） | 1e-10 | statsmodels, linearmodels, arch, numpy |
| **R 包** | 算法基准（独立实现） | 1e-10 | sandwich, rugarch, wildrwolf, modelconf, qvalue |
| **Stata/EViews** | 数值验证（闭源，学术标准） | 1e-8 | Stata `newey`, `regress, vce(robust)` 等 |
| **论文基准值** | 终极真理 | 1e-8 ~ 1e-12 | Fang-Oosterlee Table 1, Broadie-Detemple, Heston 1993 |
| **Lean4 形式化** | 数学证明（核心定理） | 精确 | Lean4 定理证明 |

#### 2.2.2 基准对齐门禁（强制）

每个模块在合并到 main 分支前，必须通过 `tests/validation/` 下的基准对齐测试：

```cpp
// tests/validation/heston_cos_table1.cpp
// 基准来源: Fang-Oosterlee (2008) "A Novel Pricing Method for European Options"
// 表 1: Heston 模型 COS 方法定价基准值
// 容差: 1e-8 (论文精度)

TEST(HestonCOSTable1, MatchesPaperBenchmark) {
    // Heston 参数 (论文 Table 1 第 1 行)
    HestonParams params{v0=0.04, kappa=1.5, theta=0.04, sigma=0.3, rho=-0.9};
    double S0=100, K=100, r=0.02, q=0, T=2.0;
    
    double price = heston_cos_price(params, S0, K, r, q, T, N=256);
    
    // 论文基准值 (Table 1, 第 1 行)
    double benchmark = 4.758404404073738;  // 来自 Fang-Oosterlee 2008
    
    EXPECT_NEAR(price, benchmark, 1e-8);
}
```

#### 2.2.3 测试目录结构与基准来源标注

```
tests/
├── unit/                          # 单元测试 (功能正确性)
│   ├── test_payoff.cpp
│   ├── test_gbm.cpp
│   └── ...
├── validation/                    # 基准对齐测试 (数值正确性)
│   ├── README.md                  # 基准来源索引 (本文件)
│   ├── bs_analytic_haug.cpp       # 基准: Haug "The Complete Guide to Option Pricing Formulas" Table 1-1
│   ├── heston_cos_table1.cpp      # 基准: Fang-Oosterlee (2008) Table 1
│   ├── heston_cf_schoutens.cpp    # 基准: Schoutens "Levy Processes in Finance" 表
│   ├── american_pde_broadie.cpp   # 基准: Broadie-Detemple (1996) 表
│   ├── leisen_reimer_convergence.cpp # 基准: Leisen-Reimer (1996) 收敛阶
│   ├── sobol_variance_reduction.cpp # 基准: Joe-Kuo 数字网均匀性
│   ├── garch_rugarch.cpp          # 基准: R rugarch::ugarchfit (v1.1)
│   ├── newey_west_statsmodels.cpp # 基准: statsmodels.stats.sandwich_covariance (v1.1)
│   ├── newey_west_stata.cpp       # 基准: Stata `newey` 命令 (v1.1)
│   ├── fdr_bh_statsmodels.cpp     # 基准: statsmodels.multipletests BH (v1.1)
│   ├── fdr_by_statsmodels.cpp     # 基准: statsmodels.multipletests BY (v1.1)
│   ├── romano_wolf_rwolf.cpp      # 基准: R wildrwolf 包 (v1.1)
│   ├── mcs_arch.cpp               # 基准: Python arch.bootstrap.MCS (v1.1)
│   ├── mcs_modelconf.cpp          # 基准: R modelconf 包 (v1.1)
│   ├── spa_arch.cpp               # 基准: Python arch.bootstrap.SPA (v1.1)
│   └── qvalue_storey.cpp          # 基准: R qvalue 包 (v1.1)
└── benchmarks/                    # 性能回归测试
    ├── bs_batch_benchmark.cpp
    └── mc_throughput_benchmark.cpp
```

#### 2.2.4 基准来源索引

每个 `tests/validation/` 下的测试文件必须在头部注释标注：

```cpp
// tests/validation/heston_cos_table1.cpp
// ============================================================
// 基准来源: Fang-Oosterlee (2008) "A Novel Pricing Method for European Options"
//           Management Science, Vol. 54, No. 5
// 表/页码: Table 1, p. 11
// 基准值: 4.758404404073738 (第 1 行, S0=100, K=100, T=2y)
// 容差: 1e-8 (论文精度)
// 交叉验证: R qvalue 包, Python arch 库
// ============================================================
```

### 2.3 阶段 3: Lean4 形式化验证（可选，核心定理）

对于核心数学定理（如 BS 公式推导、Heston 特征函数、FDR 控制率），在 Lean4 中形式化证明，与 C++ 实现对照。

**当前状态**：14 个 Lean4 定理级证明（来自用户 profile）

**Lean4 与 C++ 对照策略**：
- Lean4 证明数学定理的正确性
- C++ 实现数值算法
- 测试验证 C++ 数值输出与 Lean4 推导的解析公式一致

**适用范围**（不是所有代码都需要形式化）：
- ✅ BS 公式、Greeks 解析公式
- ✅ Heston 特征函数
- ✅ FDR 控制率证明（BH/BY/Storey）
- ✅ 收敛阶证明（MC O(1/√N), LR O(1/N)）
- ❌ 工程代码（工厂模式、内存管理、Python 绑定）

---

## 3. AI Coding 熟练度作为能力信号

### 3.1 信号分解

| 能力 | 证据载体 | 市场稀缺度 |
|---|---|---|
| **LLM Prompt 工程** | 高质量初稿（无需大量返工） | 中 |
| **架构判断力** | ADR + Phase Spec（LLM 做不到） | 高 |
| **数值验证能力** | 基准对齐 + 三源交叉验证 | **极高** |
| **形式化验证** | Lean4 与 C++ 对照 | **极高** |
| **工程判断力** | 知道何时信任 LLM、何时不信任 | 高 |

### 3.2 面试叙事

> "Cpp_Hub 采用 AI-assisted development 工作流：LLM 生成初稿，我负责架构设计、模块拆分、基准对齐和形式化验证。每个模块有独立的测试套件对比论文基准值（1e-8）和开源库输出（1e-10）。这反映的不是 C++ 编码能力，而是 AI 时代的工程判断力——知道什么时候信任 LLM，什么时候不信任。"

### 3.3 防御性回答

**Q: "这是 LLM 写的吗？"**
A: "是的，初稿由 LLM 生成。但每个模块通过基准对齐验证——对比 Fang-Oosterlee 论文 1e-8、statsmodels/R/Stata 三源交叉验证 1e-10。测试覆盖率是客观证据，无论代码是谁写的。"

**Q: "为什么不直接用 QuantLib？"**
A: "三个差异化——(a) 计量模块 QuantLib 没有，(b) 确定性并行位精确复现 QuantLib 没有，(c) Lean4 形式化 QuantLib 没有。"

**Q: "LLM 生成的代码你真的理解吗？"**
A: "每个 ADR 决策都能即兴解释。比如为什么选 Philox 而不是 Mersenne Twister——因为 Philox 是计数器 RNG，支持确定性并行，相同种子不同线程数位精确相同。这是设计决策，不是 LLM 的选择。"

---

## 4. 工作流与审计清单的集成

### 4.1 基准对齐作为强制门禁

在 [AUDIT_CHECKLIST.md](../audit/AUDIT_CHECKLIST.md) 中，每个 Phase 增加一项"基准对齐验证"作为强制门禁（编号 G）：

| 编号 | 检查项 | 标准 |
|---|---|---|
| G1 | 基准来源标注 | 每个 `tests/validation/` 文件头部有 `// 基准来源:` 注释 |
| G2 | 三源交叉验证 | 每个数值模块至少对比 2 个独立来源 |
| G3 | 容差达标 | 论文基准 1e-8，开源库 1e-10，Stata 1e-8 |
| G4 | 基准索引完整 | `tests/validation/README.md` 列出所有基准来源 |

### 4.2 与溯源报告的关系

- [TRACEABILITY_REPORT.md](../TRACEABILITY_REPORT.md)：48 项技术声明的溯源（教材页码 + 开源库文件）
- 本文件：开发工作流的溯源（LLM 生成 + 基准对齐 + 形式化验证）

两者互补：溯源报告证明"技术声明有依据"，开发工作流证明"代码实现有验证"。

---

## 5. 适用范围

### 5.1 本工作流适用的模块

- ✅ 定价模块（BS, Heston, PDE, Tree, MC）
- ✅ 计量模块（Newey-West, FDR, GARCH）
- ✅ Greeks 模块（AAD, Pathwise, LR）
- ✅ 标定模块（LM, DE）

### 5.2 本工作流不适用的部分

- ❌ Lean4 形式化证明（必须人工推导，LLM 仅辅助）
- ❌ 架构设计决策（ADR 必须人工决策）
- ❌ 基准值选择（必须人工查证论文页码）
- ❌ 审计判断（必须人工 Review）

---

## 6. 工具链

### 6.1 LLM 工具

- Claude (Anthropic)：架构设计、代码生成、文档撰写
- 本地 LLM (LM Studio, 192.168.1.11/12)：隐私敏感的代码生成

### 6.2 验证工具

- GoogleTest + ctest：C++ 单元测试 + 集成测试
- statsmodels / linearmodels / arch：Python 基准
- R (sandwich, rugarch, wildrwolf, modelconf, qvalue)：R 基准
- Stata：闭源基准验证
- Lean4：形式化证明

### 6.3 工程工具

- git + GitHub Actions：版本控制 + CI
- CMake + scikit-build-core：构建系统
- nanobind：Python 绑定
- Doxygen + Sphinx：文档生成

---

## 附录: 与 Research OS 的关系

本工作流是用户 Research OS 中"建设性怀疑主义者"方法论的具体实现：

- **不信任因子** → 四路径诊断
- **不信任模型** → CDO 算子传导
- **不信任博主** → 蒸馏框架
- **不信任自己的判断** → Lean4 形式化
- **不信任 LLM** → 基准对齐 + 三源交叉验证

所有研究产出都是"如何知道"而非"知道什么"。Cpp_Hub 的开发工作流是这一方法论在工程领域的实例化。

---

## 附录: 发现日志机制 (Discoveries Log)

### 定位

Cpp_Hub 的核心价值不只是"建造完整的库", 更是"在建造过程中发现前人遗漏的裂缝"。这与用户 Research OS 中"裂缝探测器"认知架构一致: 接触模型/理论/系统时默认进入"寻找裂缝"模式, 找到后记录并建造检测工具。

### 机制

在 `docs/discoveries/` 目录维护发现日志, 每个发现是一个独立的 `.md` 文件:

```
docs/discoveries/
├── README.md                      # 发现索引
├── 001_bh_fdr_correlated_gaussian.md  # [RESOLVED] BH-FDR 在相关双侧高斯检验下失效
├── 002_heston_cf_branch_cut_edge.md   # [OPEN] Heston 特征函数分支切割极端参数失效
├── 003_psor_omega_adaptive.md         # [OPEN] PSOR 收敛性对 ω 选择的敏感性
├── 004_sobol_bb_payoff_interaction.md # [OPEN] Sobol-BB 维度分配与 payoff 交互
├── 005_aad_checkpointing_heston.md    # [OPEN] AAD checkpointing 在 Heston 路径上的最优策略
└── 006_fp_determinism_cross_compiler.md # [OPEN] 跨编译器/跨硬件浮点确定性
```

### 发现分类

| 类型 | 含义 | 处理方式 |
|---|---|---|
| **RESOLVED** | 已有反例 + 理论解释, 可能有论文 | 记录反例 + 引用 + C++ 实现时的应对策略 |
| **OPEN** | 已识别裂缝但未解决 | 记录问题描述 + 待探索方向 + 实验计划 |
| **KNOWN** | 文献有方案但开源库未实现 | 记录文献方案 + 开源库缺口 + C++ 实现计划 |

### 与工作流的集成

1. **LLM 生成初稿时**: 若遇到数值不稳定或边界 case 异常, 不是"修 bug", 而是**记录发现**
2. **基准对齐时**: 若 C++ 输出与基准不一致, 不是"调参数", 而是**分析是否为前人遗漏**
3. **形式化验证时**: 若 Lean4 证明发现 C++ 实现的边界条件未覆盖, **记录为发现**

### 价值

1. **捕捉转瞬即逝的洞察**: 实现过程中的问题往往被遗忘
2. **将工程问题转化为研究贡献**: 每个发现都是潜在论文
3. **差异化招聘信号**: 招聘方看到的不只是"代码", 而是"发现问题的能力"
4. **与 Research OS 对齐**: "裂缝探测器"认知架构的工程实例化

### 当前发现统计

- 总发现数: 6
- RESOLVED: 1 (BH-FDR 反例, GPT-5.6 构造, Dobriban 2026)
- OPEN: 5
- 潜在论文: 5 (4 篇期刊 + 1 篇技术报告)

详见 [docs/discoveries/README.md](discoveries/README.md)。
