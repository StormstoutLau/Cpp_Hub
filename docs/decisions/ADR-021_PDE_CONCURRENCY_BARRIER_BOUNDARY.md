# ADR-021: PDE 引擎并发安全与障碍边界实施边界 (RISK-030/031, 6 项)

**状态**: Accepted (2026-10-24, Scott 裁决: 6 项按推荐定稿 — P2 直接删除 / P5 仅远场外推)
**日期**: 2026-10-24
**版本归属**: v1.x 维护线 (P2 风险修复, 非 feature phase)
**关联审计**: [EXTERNAL_REVIEW_2026-10-20.md](../audit/EXTERNAL_REVIEW_2026-10-20.md) RISK-030 / RISK-031 (P2)
**调研依据**: [RESEARCH_P2_RISKS.md](../../spec/p1-fixes/RESEARCH_P2_RISKS.md) §2-§3 (五路并行调研 + 主 agent 源码复核, 证据分级 A/B/C)
**前置 ADR**: [ADR-002](ADR_INDEX.md#adr-002-bridge--virtual-constructor-payoff) (PayOff Bridge + 虚拟构造) / [ADR-009](ADR_INDEX.md#adr-009-c-abi-稳定边界--版本化) (C ABI 稳定边界 — C++ 类不跨 ABI, 本变更不受 ABI 冻结约束)

---

## 背景

外部审计确认 PDE 模块两个 P2 风险, 同属 `pde_engine.hpp`, 合并为一次 ADR / 一次重编译批:

- **RISK-030 (并发安全)**: `pde_engine.hpp:259` `mutable Size last_total_iter_count_` 的写点在 `const` 方法 `price_american` 内 (:91/:125), 读点 `last_total_iterations()` (:255)。两个线程对同一 `const PDEEngine` 并发调用即 data race → UB ([intro.races]/17; Abseil "不存在安全的 data race")。
- **RISK-031 (边界正确性)**: `price_barrier` 远场边界按 `payoff.name() == "Call"` 字符串分派 (:186-204), 内置 5 种 payoff 中 Digital\* 全部落入 `bsm_put` 边界 — 调研探针实测 Up-Out 数字 put 偏差 **+963%**、Down-Out 数字 call 偏差 **−33%**; 且障碍路径完全不调用 `finalize_boundary`, `BoundaryType::Neumann` 被静默忽略 (实测 1e-10 内与 Dirichlet 相同)。对照: vanilla 路径 `set_boundary` 的 else 分支 (:300-304) 已有 `payoff(S)·df` 正确先例。

## 决策 (6 项)

| # | 决策点 | 决策 | 依据 |
|---|--------|------|------|
| P1 | 迭代次数返回方式 | 新增 `struct PricingResult { Real value; Size psor_iterations; }` 与主入口 `PricingResult price_american_result(...)`; 旧 `Real price_american(...)` 保留为 `return price_american_result(...).value` inline 包装 | RESEARCH §2.4①: 消除共享态, "本次调用的迭代数"语义成立 (同签名仅换返回类型不合法 → 改名而非重载) |
| P2 | `last_total_iterations()` 处置 | **直接删除** (含 `mutable` 成员与两处写点); 3 处 in-repo 测试调用迁移到 `price_american_result()` | 调用方仅 tests/unit/pricing/test_pde_engine.cpp:699/708/756 (grep 全库实证); ADR-009 C++ 类不跨 ABI; 保留成员 = UB 仍在, 弃用期只是延迟删除。拒绝 atomic/thread_local 方案 (§2.4②③: 只消 UB, 不消"计数归属"语义混乱) |
| P3 | payoff 类型识别 | `enum class PayOffKind { Call, Put, DigitalCall, DigitalPut, DoubleDigital, Unknown }`; `PayOff` 增非纯虚 `virtual PayOffKind kind() const { return PayOffKind::Unknown; }`; 5 个内置 vanilla 派生 override; 3 个 exotic (Asian/Barrier/Lookback) 保持默认 Unknown | QuantLib QuEP 2 先例 (enum 替代字符串分派, 边界策略与 payoff 解耦); 默认 Unknown = 自定义 payoff fail-fast, 杜绝 name() 字符串劫持 (§3.6: 自定义取名 "Call" 即静默劫持边界); 非纯虚 → 现有/用户派生类源码兼容 |
| P4 | 远场边界策略 (按 `kind()` 分派) | **`price_barrier` 远场端**: Call→`bsm_call`; Put→`bsm_put`; DigitalCall/DigitalPut/DoubleDigital→`payoff(S_far)·e^{-rτ}` (远场常值); **Unknown→`throw std::invalid_argument`** (自定义 payoff 可 override `kind()` 返回 Digital\* 类选用 discount 边界; 不允许映射到 BSM)。**`set_boundary` (vanilla/american 路径)**: Call→bsm_call 双端; Put→bsm_put 双端; 其余一切 kind (Digital\*/DoubleDigital/**Unknown**)→`payoff(S_min/S_max)·e^{-rτ}` — 识别从 `name()` 换成 `kind()`, 数值行为与现状逐位一致 (现状 else 分支), **不 throw**, 向后兼容 | §3.5③+④; `payoff(S)·df` 先例已在 set_boundary else 分支; 线性 payoff (Call/Put) 不走 discount 路径 (远场 `S` 项漏折扣, §3.6); 两条路径差异化: 障碍路径错边界代价 +963% 必须 fail-fast, vanilla 路径 else 分支本就正确无需破坏兼容 |
| P5 | 障碍路径 Neumann | **仅远场端实现线性外推** (step 后 `V[far] = 2·V[far∓1] − V[far∓2]`, 方向按 Down/Up 定位远端), 障碍端恒 Dirichlet 0; 不再静默忽略; 不 throw | §3.5: 禁止静默忽略; 与 vanilla 路径 Neumann 语义对齐; 相比 throw 保全功能 (远场外推是 Neumann 在该路径唯一合理位置; 障碍端必须保持 0, 外推会抹掉, §3.6) |
| P6 | 测试锚 | (a) 数字障碍收敛测试: 网格加倍 + 域扩大双层收敛, 对照宽域参考 (§3.3 探针单点数值不复用 — 证据局限: 单参数组/单网格族); (b) 障碍路径 Neumann≠Dirichlet 断言; (c) Unknown kind 在 `price_barrier` throw 断言 + 自定义 payoff override `kind()` 走 discount 边界断言; (d) `price_american_result().psor_iterations > 0` 且 `.value` 与旧包装位相同; (e) vanilla 路径自定义 payoff (Unknown) 行为回归断言 (不 throw, `payoff(S)·df`) | TDD 先红后绿; 探针证据 B 级仅作方向参考, 正式锚以收敛测试为准 |

## 影响面 (已核实)

- **变更文件**: `include/cpphub/instruments/payoff/payoff.hpp` (PayOffKind + kind()) / `include/cpphub/instruments/payoff/vanilla.hpp` (5 处 override) / `include/cpphub/pricing/pde/pde_engine.hpp` (P1/P2/P4/P5) / `tests/unit/pricing/test_pde_engine.cpp` (3 处迁移 + P6 新测试)
- **零影响面**: python 绑定仅暴露 price_european/price_american/greeks (`cpphub_python.cpp:1116-1141`), `price_american` 走 Real 包装签名不变; `price_barrier` 与 `last_total_iterations` 未暴露; C API 不涉及 PDE; MCP registry/golden 不涉及 PDE; `bench_pde.cpp` 不调用 `last_total_iterations` (grep 实证)
- **行为变更点** (均 fail-fast 化, 在 DEVELOPMENT_LOG 声明):
  1. 障碍路径数字期权结果大幅改变 (+963% 偏差回归正确 — 修复目的本身)
  2. 障碍路径 `BoundaryType::Neumann` 不再等同 Dirichlet
  3. 自定义 payoff (kind()==Unknown) 调 `price_barrier` 由静默错误边界改为 throw; vanilla 路径行为不变

## 后果

- **正向**: `const PDEEngine` 并发调用 UB 消除; 数字障碍从 +963% 错误回归正确; Neumann 全路径语义一致; payoff 类型识别不再可被字符串劫持
- **负向/成本**: `last_total_iterations()` 即刻删除 (若裁决选弃用期, P2 改为 `[[deprecated]]` + 成员保留一版, UB 相应延后消除); 障碍路径既有数值断言(如有)按 P6 新锚重算; `PayOff` 增加 1 个虚函数 (header-only, 无 ABI 实际影响)
- **开放问题**: 无 — 6 项均给定推荐, 待裁决生效后按 P6 以 TDD 实施

## 关联

- **调研**: [RESEARCH_P2_RISKS.md](../../spec/p1-fixes/RESEARCH_P2_RISKS.md) §2 (RISK-030) / §3 (RISK-031)
- **审计**: [EXTERNAL_REVIEW_2026-10-20.md](../audit/EXTERNAL_REVIEW_2026-10-20.md) §二 主表 RISK-030/031
- **前置 ADR**: [ADR-002](ADR_INDEX.md#adr-002-bridge--virtual-constructor-payoff) / [ADR-009](ADR_INDEX.md#adr-009-c-abi-稳定边界--版本化)
- **批次**: RESEARCH_P2_RISKS.md §6 第 2 批 (RISK-029 已单独实施为 b2339fd; RISK-035/036 为第 3 批, 数值契约修正无架构决策, 走 CHECKLIST 项不设 ADR)
