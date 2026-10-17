# Cpp_Hub 项目文档索引

> 版本: 1.0  
> 更新: 2026-07-29  
> 状态: 文档阶段完成，可进入代码实施阶段

---

## 命名约定

| 上下文 | 名称 | 说明 |
|--------|------|------|
| 仓库目录 | `Cpp_Hub` | 保留大小写，GitHub 仓库名 |
| CMake project | `CppHub` | CMakeLists.txt 中 `project()` |
| C++ 命名空间 | `cpphub` | 全小写，无下划线 |
| Python 包名 / PyPI 包名 | `cpphub` | `import cpphub` |
| CMake target | `cpphub` / `cpphub_core` / `cpphub_python` | 小写+下划线分隔 |

## 权威文档声明

以下文档为**当前权威**，所有实施以此为准：

- `BUILD_PLAN.md` - 总体构建计划
- `docs/architecture/ARCHITECTURE.md` - 架构设计
- `docs/phases/phase{1-4}/PHASEx_SPEC.md` - 阶段执行规格书
- `docs/decisions/ADR_INDEX.md` - 架构决策记录
- `docs/audit/AUDIT_CHECKLIST.md` - 代码审计检查表

> ⚠️ `PROJECT_PLAN.md` 已标记为 **DEPRECATED**，仅作历史参考，不再维护。

---

## 版本路线图

| 版本 | 对应 Phase | 核心内容 | 预计工期 |
|------|-----------|----------|----------|
| **v1.0** | Phase 1 + Phase 2 (削减后) | BS/Heston + Analytic/MC/PDE/Tree + Sobol QMC + 基础 Greeks + Python 绑定 | 6-8 周 |
| **v1.1** | Phase 3 (+ Phase 2 推迟内容) | 进阶模型 (SABR/IR/Levy) + COS/FFT/LSMC + AAD Greeks (autodiff) + VaR/ES + 标定 + 波动率曲面 | 3-4 周 |
| **v2.0** | Phase 4 | GPU/分布式/Excel/云原生 (gRPC+Prometheus) /多平台发布 | 3-4 周 |

**开发策略**: 先构建 v1.0，后续迭代升级到 v1.1、v2.0。各版本独立可发布。

---

## 📁 文档目录结构

```
docs/
├── ARCHITECTURE.md              # 总体架构文档
├── BUILD_PLAN.md                # 构建计划总览 (含 v1.0/v1.1/v2.0 路线图)
├── PROJECT_PLAN.md              # [DEPRECATED] 早期规划草稿
├── TRACEABILITY_REPORT.md       # 源码级溯源验证报告
├── DEVELOPMENT_LOG.md           # 开发日志模板 (70天详细计划)
├── AUDIT_CHECKLIST.md           # 代码审计检查表 (4阶段全覆盖, 含 G 类基准对齐门禁)
├── DEVELOPMENT_WORKFLOW.md      # 开发工作流 (LLM 辅助开发 + 基准对齐 + 形式化验证)
├── ADR_INDEX.md                 # 架构决策记录索引 (14项 ADR, 含 ADR-013 双层 linalg, ADR-014 标定 vs 估计)
├── research/                    # [本地留存] 调研报告(含研究 IP) — .gitignore 排除, 不入公共仓库
│   ├── ECONOMETRICS_LANDSCAPE.md  # 计量经济学 C++ 开源生态调研报告 (含参数估计 §11)
│   └── MULTI_AGENT_TRIAL.md       # 多 Agent 能力验证任务 (A/B 站 qwen3-coder + GPT OSS 120B)
├── architecture/
│   └── ARCHITECTURE.md          # 详细架构设计 (分层/契约/数据流/部署)
├── phases/
│   ├── phase1/
│   │   └── PHASE1_SPEC.md       # Phase 1 执行规格书 (Week 1-3: 核心内核 MVP)
│   ├── phase2/
│   │   └── PHASE2_SPEC.md       # Phase 2 执行规格书 (Week 4-6: 进阶模型/数值方法)
│   ├── phase3/
│   │   └── PHASE3_SPEC.md       # Phase 3 执行规格书 (Week 7-9: 风险/标定)
│   └── phase4/
│       └── PHASE4_SPEC.md       # Phase 4 执行规格书 (Week 10-12: 高性能/发布)
├── audit/
│   └── AUDIT_CHECKLIST.md       # 审计清单 (架构/C++规范/数值/性能/Python/测试)
├── decisions/
│   └── ADR_INDEX.md             # 14项 ADR 详细记录
└── discoveries/                 # [本地留存] 发现日志(含潜在论文素材) — .gitignore 排除
    ├── README.md                # 发现索引 (6项: 1 RESOLVED + 5 OPEN)
    ├── 001_bh_fdr_correlated_gaussian.md  # [RESOLVED] BH-FDR 在相关双侧高斯检验下失效
    ├── 002_heston_cf_branch_cut_edge.md   # [OPEN] Heston 特征函数分支切割极端参数失效
    ├── 003_psor_omega_adaptive.md         # [OPEN] PSOR 收敛性对 ω 选择的敏感性
    ├── 004_sobol_bb_payoff_interaction.md # [OPEN] Sobol-Brownian Bridge 维度分配与 payoff 交互
    ├── 005_aad_checkpointing_heston.md    # [OPEN] AAD checkpointing 在 Heston 路径上的最优策略
    └── 006_fp_determinism_cross_compiler.md # [OPEN] 跨编译器/跨硬件浮点确定性
```

### 仓库级 feature 文档管道（Spec_Workflow 规范）

```
spec/                              # 仓库根目录 (不在 docs/ 下), 入 Git 版本控制
└── <feature>/
    ├── RESEARCH.md                # Step 2 产出: 调研文档 (A/B/C 断言分级 + §0 统计表)
    ├── DESIGN.md                  # Step 3 产出: 设计文档
    ├── IMPLEMENTATION.md          # Step 5 产出: 实施文档
    └── CHECKLIST.md               # Step 7 产出: 验收 checklist
```

- 现行实例：`spec/mcp-server/RESEARCH.md`（MCP 接口层调研，Step 1–2）
- **版本化纪律**：feature 文档入 Git（与 `docs/research/`、`docs/discoveries/` 的本地留存策略相反）；文档内**不得**引用被 `.gitignore` 排除的本地路径——须内联内容或显式标注「本地留存」

---

## 🎯 核心文档速览

| 文档 | 核心内容 | 关键指标 |
|------|----------|----------|
| **ARCHITECTURE.md** | 分层架构、模块契约、数据流、ADR 摘要 | 8层架构、Strategy/Bridge/Factory 模式 |
| **BUILD_PLAN.md** | 12周/4阶段里程碑、交付物、验收标准 | M1-M12 里程碑、性能基准 |
| **TRACEABILITY_REPORT.md** | 48项技术声明溯源：65%教材、35%开源库、0%幻觉 | 零幻觉、每项标注页码/库名 |
| **PHASE1_SPEC.md** | Week 1-3 核心内核：基础设施→PayOff/MC/解析解→Python/CI | 13个核心文件、数值基准 |
| **PHASE2_SPEC.md** | Week 4-6 Heston/SABR/IR、PDE/树形/傅里叶、LSMC/标定 | COS指数收敛、Leisen-Reimer O(1/n²) |
| **PHASE3_SPEC.md** | Week 7-9 AAD Greeks、VaR/ES、联合标定、SVI无套利 | 4种Greeks方法一致性1e-6 |
| **PHASE4_SPEC.md** | Week 10-12 GPU/分布式/Excel/云/发布 | GPU 500M paths/s、多平台Wheel |
| **PHASE5_HFE_SPEC.md** | v1.4 高频计量：TAQ reader / RV·BPV·RSV / 12 核函数 Realized Kernel / 5 多资产 Cov / HAR·HEAVY / 流动性 23 / 跳跃检验 AJ·JO·Intraday·Rank | 对标 R `highfrequency` 1.0.3，容差 1e-12，HFE 144 项 |
| **PHASE6_ECONOMETRICS_SPEC.md** | 证伪统计量 + 计量推断（Newey-West / FDR / GARCH / 多变量波动率） | 统一 `TestResult` 结构，D1–D23 排幻觉点 |
| **PHASE7{A,B,C}_*.md** | 时序与多变量：单位根六检验 / 协整 EG·Johansen / VECM / Granger·DY 溢出 / GARCH 族 / MIDAS | ADR-016~019，三平台一致 |
| **PHASE8_SPEC.md** | ARDL/PSS 边限检验 + Kalman 状态空间精确似然 | ADR-020，v1.8.0 tag |
| **spec/mcp-server/RESEARCH.md** | MCP 接口层调研（Step 1–2）：规范 2026-07-28 / 形态 B1→B2 / 诊断类工具面 | A=22 / B=5 / C=4 断言分级 |
| **AUDIT_CHECKLIST.md** | 4阶段加权审计：架构20%/C++15%/数值30%/性能15%/Python10%/测试10% | 量化评分、签名确认 |
| **ADR_INDEX.md** | 20 项决策记录（ADR-001~020）：Header-only、Bridge、Factory、RNG、SIMD、nanobind、AAD、COS、C ABI、GPU、MPI、xlOil、双层linalg、标定vs估计、证伪边界、时序命名空间、SLSQP、ARDL/Kalman | 状态跟踪、替代方案评估 |

---

## 🔗 关键追踪链路

```
教材内容 (Joshi/Duffy×2)
    ↓ 溯源验证 (TRACEABILITY_REPORT.md)
架构决策 (ADR_INDEX.md, ADR-001~020)
    ↓ 落地实施
Phase 1-8 执行规格书 (docs/phases/PHASEx_SPEC.md)
    ↓ 质量门禁
代码审计清单 (AUDIT_CHECKLIST.md)
    ↓ 过程记录
开发日志 (DEVELOPMENT_LOG.md)

新 feature 文档管道 (Spec_Workflow):
RESEARCH.md → DESIGN.md → IMPLEMENTATION.md → CHECKLIST.md   (spec/<feature>/)
```

---

## 🚀 项目状态

### 交付进度 (v1.0 → v1.8)

- ✅ Phase 1-3 (v1.0): 286/286 测试通过 (Core/PayOff/BS/Heston/MC/Sobol QMC/PDE/Tree/Greeks/VaR/Calibration/SVI)
- ✅ Phase 4 LITE (v1.1): 320/320 测试通过 (SSVI 跨期限 / Python 绑定 / GPU MC)
- ✅ Phase 5 (v1.4.0–v1.4.3): HFE 高频计量 144 项 (RV/BPV/RSV · 12 核函数 Realized Kernel · 5 多资产 Cov · HAR/HEAVY · 流动性 23 · 跳跃检验 AJ/JO/Intraday/Rank)
- ✅ Phase 6/7 (v1.5–v1.7): 证伪统计量 + 时序/多变量 (单位根 · 协整 · VECM · Granger/DY · GARCH 族 · MIDAS)
- ✅ Phase 8 (tag `v1.8.0`): ARDL/PSS 边限检验 + Kalman 状态空间精确似然
- ✅ 跨平台: MSVC (Win10) + GCC (Ubuntu 22.04/24.04) 位精确一致 (HFE 1412/1412)
- ✅ ADR: 20 项 (ADR-001~020)

### 下一步方向

| 版本 | 范围 | 状态 |
|------|------|------|
| **v1.9 (Phase 9)** | KDE/KNN + KSG/TE 信息论估计器（D-2 裁决载体） | 规划中 |
| **v2.0+** | MPI 分布式 / Excel XLL / gRPC / 云原生 / 多平台 Wheel；MCP 接口层（`spec/mcp-server/`） | 推迟 / 调研中 |

---

## 📊 验收门禁摘要

| 阶段 | 数值基准 | 性能基准 | 工程基准 |
|------|----------|----------|----------|
| **Phase 1** | BS 1e-12, MC 收敛阶 -0.5±0.05, Sobol 10x方差缩减 | 批量BSM >30M/s, MC >25M paths/s | 0警告, 覆盖率>90%, ASan/TSan 0 |
| **Phase 2** | Heston COS 1e-8, PDE 1e-4, 树形 O(1/n²) | Heston COS <5ms, PDE 400x2000 <100ms | 单测新增150+, 回归守护5% |
| **Phase 3** | AAD vs 解析 1e-10, Greeks四法一致 1e-6, VaR回测 p>0.05 | 万期权全Greeks <100ms | 标定目标函数<1e-6 |
| **Phase 4** | CPU/GPU 双精度位精确, 分布式线性加速>0.9 | GPU MC >500M paths/s, XLL UDF <5μs | 多平台Wheel, Conda/vcpkg/Conan同步 |
| **Phase 5–8** | R/论文基准对齐 (HFE 1e-12~1e-8；时序 1e-8~1e-14 三方交叉) | CI 三平台构建 + 测试全绿 | R 门禁：调研断言分级 + 排幻觉点逐点核验 |

---

## 📝 维护规范

1. **文档版本控制**: 受控 `.md` 纳入 Git，重大变更需 ADR 记录；`docs/research/`、`docs/discoveries/` 为**本地留存**（研究 IP，`.gitignore` 排除），新 feature 文档走 `spec/<feature>/` 入 Git
2. **溯源标注**: 每新增 `.hpp` 必须包含 `// SOURCE:` 头部注释
3. **日志更新**: 每日下班前更新 `DEVELOPMENT_LOG.md`，周五生成周报
4. **审计节点**: 每阶段结束按 `AUDIT_CHECKLIST.md` 逐项核对，Reviewer 签名
5. **发布冻结**: 发布前所有 ADR 必须 Accepted，无 Proposed 残留
6. **本地留存引用惯例**: 受控文档**允许**引用本地留存路径（如 `../research/xxx.md`）作为 provenance——该引用在 GitHub 上呈现为断链，属**有意取舍**（研究 IP 不入库优先于链接完整性）；引用方应在正文同时给出可独立理解的关键结论
7. **DC1 文档契约适用范围**: front-matter（`id` / `type` / `version` / `status` / `date` / `depends` / `upstream`）**仅要求** `spec/<feature>/` 下的新 feature 文档；存量 `docs/**` 文档为历史遗留，不追溯回填
8. **第三方 fixture 来源标注**: `tests/fixtures/` 引入外部出版物副本时须注明来源与可获得性（例：`timeseries/em2002_ifdp655.pdf` = FRB *International Finance Discussion Paper* #655，美联储公开出版物；`hfe/hf_source/` 按 CRAN 版权已排除入库）

---

**文档集完成度**: 核心文档已随 v1.8 更新（本索引 2026-10 刷新）  
**当前版本**: **v1.8.0**（Phase 1–8 交付，ADR-001~020）  
**下一里程碑**: v1.9 (Phase 9, KDE/KNN + KSG/TE) / MCP 接口层（`spec/mcp-server/` 调研中）