# 断言分级证据框架 (A/B/C) — 调研与审计工作流

> **版本**: v1.0 (2026-08-16)
> **来源**: Phase 7C 调研审计复盘 (126 条声明 → 8 处实质错误) + NP2 τ_T(k) MCP 学术搜索裁决
> **用途**: 所有后续调研 agent / 审计 agent 的 prompt 约束与核验流程基准
> **核心原则**: 生成端强制证据，审计端不信任引文 — 不对称配置

---

## 1. 断言三级分类定义

| 级别 | 定义 | 证据形式 | 缺证据时处理 |
|------|------|---------|-------------|
| **A 事实类** | 单点外部可验证: 版本号/公式/参数语义/章节页码/函数签名/源码行为 | **URL + 原文引文 (≤3 行, 可 grep)** | 移入"假设区", **禁止入正文** |
| **B 推断类** | 综合多源推理得出: "X 与 Y 一致/不同" / "存在/不存在 Z" / "共 N 套" | **显式推理链 + 依赖源各自附 A 级证据** | 移入假设区; 推理链缺一步即降级 |
| **C 判断类** | 决策与权衡: 对照顺序/优先级/scope 取舍 | rationale + 假设声明 (无需链接) | 照常陈述, 但标 [判断] |

**附加规则**:
- **双源规则**: spec 阻断性断言 (会影响实现正确性的) 必须 ≥2 独立来源, 否则标 `[单源-待二核]`
- **引文可核验性**: 引文须为原文精确片段 (保持拼写/连字), 不得转述; PDF 需容错 "M AIC" 式分拆排版
- **身份验证**: 任何自动下载的证据文件必须先验证首页身份 (EuropePMC PMID 错配教训: DOI 查询返回过无关 PLOS One 论文)

---

## 2. 八处错误复盘 (2026-08-16 Phase 7C 案例, 7A + 1B)

| # | 错误 | 级别 | 当时证据缺陷 | 正确证据要求 (若当时即满足则拦截) |
|---|------|------|-------------|--------------------------------|
| 1 | AR2: CSS n.cond 写 max(p,q)+1 | **A** | 凭文档印象, 未引源码行 | arima.R L158-162 链接 + 引文 `ncond = d + D*s + max(user, p + s*P)` — 摘引时即见与 q 无关 |
| 2 | V8: 称 σ_ii⁻¹ 是 "R 包误记" | **A** | 凭记号直觉, 未引 PS 原文 | PS 1998 PDF eq.(12) 引文 `θ^g_ij(n) = σ_ii⁻¹·Σ(e_i′A_lΣe_j)²/Σ e_i′A_lΣA_l′e_i` — 原文分子就是 σ_ii⁻¹ |
| 3 | CI5: 称 statsmodels 有 "三套" 临界值表 | **B** | **两处源码引用均正确, 综合推理错** (未检查 coint_johansen 是否直接调用 c_sjt/c_sja) | 推理链中 "内嵌表 ≠ c_sjt" 一环需赋值语句证据: vecm.py L724-725 `cvm[i,:]=c_sja(...)` — 赋值即证伪 "不同表" |
| 4 | midas_r 默认 Nelder-Mead | **A** | 只引了函数签名 (Ofunction="optim"), 未读函数体 | midasreg.R L858-860 引文 `## Override default method...` + `control$method <- "BFGS"` |
| 5 | zivot_andrews "0.20+" | **A** | 无 release notes 链接, 版本号凭印象 | statsmodels 0.11.0 release notes 链接 + 引文 "The Zivot-Andrews test ... has been added" |
| 6 | ZA MC c 1% = −5.83 | **A** | 引了 MC 表但未引分位数标签 | zivot_andrews.py 表源码 + 分位标签行: −5.27644@1% / −5.83192@0.1% |
| 7 | Lütkepohl Granger = §7.2.2 | **A** | 章节号凭记忆 | TOC 链接 + 引文 "3.6.1 A Wald Test for Granger-Causality ... 102" / "7.2.2 EGLS Estimation of the Cointegration Parameters ... 291" |
| 8 | fpp2 差分/drift = §8.7 | **A** | 同上 | otexts.com/fpp2 TOC + 引文 "8.1 Stationarity and differencing" |

**结论**: A 类 7/8 若生成期强制 "链接+引文" 当场拦截; 唯一 B 类 (CI5) 证据全对但推理跳步 — 必须靠 B 类专用审计 (见 §4)。

**同类案例 (非 8 处之内)**: NP2 τ_T(k) 的 λ̂−λ̃ 形式 = A 类无源断言 (全源零命中), 强制证据规则下根本无法写出, 会在生成期自曝为假设。

---

## 3. 调研 Agent Prompt 约束模板 (可直接粘贴)

```markdown
【证据纪律 — 本任务所有输出强制执行】

你产生的每一条断言必须标注级别并附规定形式的证据, 否则不得写出:

A 事实类 (单点外部可验证: 版本/公式/参数语义/章节页码/源码行为):
  → "断言 | A | URL | 引文(≤3行原文, 精确到可 grep, PDF 注意 'M AIC' 式分拆排版)"
  → 无 URL+引文 = 移入"假设区", 禁止出现在正文结论中

B 推断类 (综合多源: "X 与 Y 一致/不同"/"存在 Z"/"共 N 套"):
  → "断言 | B | 推理链(编号步骤, 每步注明依赖哪个源) | 各依赖源的 A 级证据"
  → 禁止跳步: "X 调用 Y" 类传递性结论必须有调用/赋值语句级证据
  → 存在性断言 ("无库实现 X"): 必须列出搜索过的库名清单与零命中结果

C 判断类 (决策/优先级/scope):
  → "断言 | C | rationale + 假设声明" (无需链接)

通用规则:
1. 阻断性断言 (错误会传导进 spec/实现的) 需 ≥2 独立来源, 单源标 [单源-待二核]
2. 引文禁止转述, 必须原文; 自动下载的 PDF 必须先验证首页身份 (DOI 可能错配)
3. 无法获得证据时诚实标注 ❓+查证路径, 不得以弱记忆填充
4. 输出末尾附: 断言统计表 (A/B/C 各多少条, 假设区条目列表)

自检 (提交前逐条过):
- [ ] 每条 A 类断言的引文能在所附 URL 页面内 grep 到?
- [ ] 每条 B 类断言的推理链每一步都有源, 无 "显然/众所周知"?
- [ ] 版本号/章节号/数值常量三类各有链接? (最高频错误位)
- [ ] 假设区与正文严格分离?
```

---

## 4. B 类独立重推导审计脚本逻辑

> B 类错误特征: 依赖源全部正确, 综合推理跳步 (CI5: 两处源码引用都对, "三套表"结论错)。
> 链接核验对 B 类无效 — 必须独立重推导 + 机械反证探针。

### 4.1 设计原理

```
B 类断言 = 结论 C = F(源集合 S = {s1..sn}, 算子 op)
CI5 案例解剖:
  C  = "statsmodels 有三套临界值表"
  S  = {coint_johansen 返回 cvt/cvm (对), select_coint_rank 用 c_sjt/c_sja (对)}
  op = 计数+互斥推断 ("内嵌表 ≠ c_sjt" ← 这一步无证据, 跳步)
  反证 = vecm.py L724-725 赋值语句 `cvm[i,:]=c_sja(...)` 直接证伪互斥假设
```

**核心: 按 op 类型生成机械反证探针 (falsification probe), 不依赖审计者的聪明。**

### 4.2 算子-探针对照表

| op 类型 | 典型断言形态 | 机械反证探针 | CI5 命中方式 |
|---------|-------------|-------------|-------------|
| **等价/互斥** | "X 与 Y 相同/不同/无关" | 调用图 + 赋值语句 grep: `X.*=.*Y` / `Y\(.*\)` 在 X 定义文件内 | ✅ `cvm[i,:]=c_sja(...)` |
| **存在/不存在** | "没有库实现 Z" | 否定搜索: 枚举候选库名单 → 全库 grep 关键符号 → 零命中才维持 | (NP 无 M 族即此型) |
| **计数** | "共 N 套/种/个" | 枚举再数: 生成候选全集 (非记忆), 逐个验证归属后计数 | ✅ 枚举 statsmodels 内所有 cv 表定义点 |
| **传递/依赖** | "A 基于 B" / "A 用的是 B 的结果" | import/调用链静态追踪, 缺一环即 flag | |
| **跨库一致** | "库1 与库2 公式相同" | 同输入数值 diff (构造最小样本两端跑, 比对) | |
| **因果/排序** | "因 X 所以 Y" | 反事实检查: X 不变 Y 变的对照例存在即降级 | |

### 4.3 双盲重推导流程 (脚本骨架)

> **可执行实现**: `scripts/assertion_audit.py` (本地, 不进公共仓库 — 同 `scripts/` 目录惯例) — 完整实现 5 类探针 + 双盲重推导 + STEP_GAP 检测 + 仲裁 prompt 生成; 内置 CI5/NP/计数/STEP_GAP 四个离线自检工作示例, 运行 `python scripts/assertion_audit.py demo` 验证。auditor 可插拔 (manual / OpenAI 兼容端点, 支持 LM Studio 三机)。

```python
def audit_class_b(assertion, source_evidence, auditor_agent):
    """
    assertion:       {conclusion, op_type, claimed_chain}
    source_evidence: 各依赖源的 A 级证据 (已核验为真)
    auditor_agent:   独立 agent (未见过原报告的推理文本)
    """
    # Phase 1: 机械探针 (脚本, 零 LLM, 最便宜)
    probe = PROBE_REGISTRY[assertion.op_type]        # 4.2 对照表
    falsifier = probe(assertion, source_evidence)     # 返回反证证据或 None
    if falsifier:
        return Verdict(False, kind="mechanically-falsified",
                       evidence=falsifier)             # CI5 在此终结

    # Phase 2: 双盲重推导 (剥离原推理, 只给源与命题)
    prompt = f"""
    命题: {assertion.conclusion} 是否成立?
    材料: 仅以下源证据 (无任何其他人推理):
    {source_evidence}
    要求: 给出你自己的推导链 (编号), 每步标注用到的源;
    输出 verdict: TRUE / FALSE / UNCERTAIN + 关键一步的理由。
    """
    r = auditor_agent(prompt)
    # r.chain 与 claimed_chain 逐步对比:
    #   - 结论不一致           → flag: CONFLICT (进仲裁)
    #   - 结论一致但步骤数不同   → flag: STEP_GAP (原链可能跳步, 复查差额步)
    #   - 完全一致             → PASS (标注 "双盲复核通过")

    # Phase 3: 仲裁 (仅 CONFLICT/STEP_GAP 触发)
    #   第三 agent 或人工: 只裁决分歧步, 不重跑全链
```

### 4.4 关键工程约束

1. **探针优先于 LLM**: Phase 1 脚本化 (grep/调用图/数值 diff) 零幻觉零成本, 能机械终结的不过 LLM
2. **双盲是硬条件**: Phase 2 的 auditor 不得接触原推理文本 — 否则锚定效应会把跳步复制一遍 (这正是 v1.1 审计有效的原因: "未采信报告自带引用")
3. **STEP_GAP 不是通过**: 原链 3 步、重推 5 步且结论相同时, 差额的 2 步正是 CI5 型跳步藏身处, 必须复查
4. **登记制**: B 类断言在调研产出中必须以 4.1 的结构化形式登记 {结论, op, 依赖源, 推理链}, 否则审计脚本无输入

---

## 5. 审计端不对称配置 (生成引证 ≠ 审计采信)

| 断言级 | 生成端 | 审计端 |
|--------|--------|--------|
| A 事实类 | 强制 URL+引文 | **机械核验** (可脚本化: 链接存活 + 引文在页内 grep + 首页身份) + 抽样人工 |
| B 推断类 | 强制推理链 | **§4 全流程** (探针 + 双盲重推导) — 链接核验无效 |
| C 判断类 | rationale | 只查假设是否显式声明, 不裁决对错 |
| 阻断性 | 双源 | 双源各自独立核验 + 来源间一致性 diff |

**红线**: 规则普及后最大退化风险 = 审计变成 "链接存在即通过" 的形式审计。
引文核验只证明 "生成者看过该页", 不证明 "结论从该页可推出" — B 类永远需要独立重推导。

---

## 6. 效率账 (基于 Phase 7C 实测反事实)

| 指标 | v1.0 流程 (无证据强制) | 本框架预期 |
|------|----------------------|-----------|
| 错误拦截位置 | 审计期 (8 处已入报告) | 生成期拦 ~7/8 (A 类), B 类由探针+双盲拦 |
| 审计成本 | 3 agent × 126 条全量独立重查 | 1 遍机械核验 (脚本) + 1 agent 攻 B 类重推导 + 双源抽查 |
| 审计输入 | 散文混断言 | 结构化断言表 (可直接喂 4.3 脚本) |
| 无法消除项 | — | 双方误读同一来源 (→双源规则), 推理型残余 (→STEP_GAP 复查) |

---

## 7. 调研 Agent 报告模板 Prompt 约束 (可直接粘贴)

> 目的: 让调研报告本身成为审计工具的直接输入 — `python scripts/assertion_audit.py audit --input <报告.md>`
> 自动提取附录 B 的 ```assertions 机读块 (已实现并验证)。生成端与审计端共用同一份状态词表与 ID 体系。

```markdown
【报告模板约束 — 你的 Markdown 输出必须遵循以下骨架: 逐节存在, 顺序不变, ID 稳定不复用】

# <主题> 调研报告 v<x.y> (<日期>)

## 0. 断言统计表 (必填, 审计入口)
| 级别 | 条数 | 说明 |
|------|------|------|
| A 事实类 | n | 每条附 URL+引文 |
| B 推断类 | m | 每条有 ID, 登记于附录 B |
| C 判断类 | k | 仅 rationale |
| 假设区 | h | 无证据断言, 禁止进入正文 |

正文规则 (逐条自检后才算完成):
R1 每条断言行内标注 【A】/【B#ID】/【C】 之一; 无标注句子不得含可验证事实
R2 【A】紧跟 "(源: URL; 引文: ...≤3行原文)" — 引文须能在源页 grep 到
R3 【B#ID】正文只写结论一句; 证据与推理链只出现在附录 B 对应 ID
R4 阻断性断言 (错误会传导进 spec/实现) ≥2 独立源; 单源标 [单源-待二核]
R5 B 类断言禁止 "显然/众所周知/公认"; B 类推理链禁止无源步骤
R6 假设区每条含 [H#] + 查证路径; 正文出现 FALSIFIED 断言时必须改写并记入修订记录

## 1..N 正文章节 (按调研域组织, 断言标注规则同上)

## 附录 B: 断言登记表 (机器可读 — 字段名与 assertion_audit.py 严格一致)
```assertions
[
  {
    "id": "B1",
    "conclusion": "<结论陈述一句>",
    "op": "equivalence|existence|counting|transitivity|cross_library|causal",
    "claimed_chain": [
      {"step": 1, "text": "<步骤>", "source": "<依赖源label, 无源则为null>"}
    ],
    "sources": [
      {"label": "<源名>", "path": "<本地路径|null>", "url": "<URL|null>", "quote": "<引文|null>"}
    ],
    "probe": {"type": "<与op对应>", "files": ["<文件/目录>"], "params": {"<op专属参数>"}}
  }
]
```
(op→probe.params 对应: equivalence→falsifier_pattern+direction; existence→symbols+candidates+claim;
 counting→definition_pattern+expected_count; transitivity→entry+target; cross_library→cmd_a+cmd_b+keys+tol;
 causal 无机械探针, probe 置 null)

## 附录 C: 假设区
- [H1] <断言> (查证路径: <具体到工具/URL/文献>)

审计闭环 (报告写完后必须执行并追加):
  python scripts/assertion_audit.py audit --input <本报告.md> \
      --auditor openai --base-url <LM Studio端点> --report <审计输出.md>
审计结论以 "## 审计结论 (<日期>)" 章节追加回本报告末尾, 逐断言给出:
  | ID | 最终状态 | 证据/分歧 |
状态词表 (固定, 不得自造):
  FALSIFIED 机械证伪 / SURVIVED 存活 / CONFLICT 双盲结论相反
  STEP_GAP 疑跳步 / UNCERTAIN 无法判定 / PENDING 待人工 / NO_PROBE 无机械探针
```

### 7.1 双端词汇对照 (生成端 ↔ 审计端)

| 生成端登记 | 审计端产出 | 闭环动作 |
|-----------|-----------|---------|
| 【B#ID】+ registry 行 | verdicts 按 ID 回填 | FALSIFIED → 正文改写; STEP_GAP/CONFLICT → 仲裁 |
| [单源-待二核] | 双源规则触发记录 | 补第二源或降级假设区 |
| [H#] 假设区 | 不审计 (无证据) | spec 前须转为 A/B 或清除 |

### 7.2 审计产出报告骨架 (write_report 生成, 调研 Agent 追加时保持同构)

```markdown
# B 类断言审计报告
## <ID> — <最终状态>
- 结论: <conclusion>  |  op: <op>
- **probe** [<状态>]: <detail>            (+证据代码块)
- **double_blind** [<状态>]: <detail>      (STEP_GAP 时列跳步候选)
- **arbitration** [NEEDS_ARBITRATION]: <仅分歧步的仲裁任务>
## 汇总
| ID | 最终状态 |
```
