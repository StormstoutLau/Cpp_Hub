# RESEARCH: pss_critical_values.hpp 深尾缺陷（`1−Φ` 族第 4 例）

---

id: P1-FIXES-RESEARCH-PSS-DEEPTAIL
type: research
version: 1.0
status: pending-approval
date: 2026-10-26
upstream: [RESEARCH_P2_C1128_AND_DEFERRED.md §6 遗留登记](./RESEARCH_P2_C1128_AND_DEFERRED.md) · DEVELOPMENT_LOG「P2 遗留专项批」条目
evidence-framework: docs/ASSERTION_EVIDENCE_FRAMEWORK.md（A = 本地实证/源码直读，B = 官方文档/转录，C = 社区旁证）

---

> **一句话结论**：[pss_critical_values.hpp:70](file:///f:/Cpp_Hub/include/cpphub/timeseries/ardl/pss_critical_values.hpp#L70) 的 `1.0 - normal_cdf(poly)` 是深尾 `1−P` 缺陷族第 4 例（前 3 例 = RISK-036 χ² 上尾 / SMM 联合 p / D-1 两尾正态 p）。**比前三例严重得多**：归零阈值低至现实可达的 F 统计量区间（锚键 stat≈55，全 100 键最激进 stat≈10），且 **100/100 键全部中招**。修复 = 一行上尾直算，无 golden REGEN 面，上游 statsmodels 同病（修复后本库严格优于上游）。
> **状态**：调研完成，**待用户批准后实施**（TDD 先红后绿）。

---

## 0. 结论速览

| 项 | 裁决 | 证据级 |
|---|---|---|
| 缺陷确认 | `1.0 - normal_cdf(poly)` 在 poly ≥ 8.25 时精确归零，poly ≥ 8.06 起误差 >10% | A（§2 实证） |
| 严重性 | **P1**（信息性 p 值在现实可达 stat 区间归零；三值结论判定不受污染 → 不到 P0） | A（§4） |
| 与前 3 例对比 | 前三例需 \|z\|≥8~10（VR/SMM 尺度），本例只需 stat≥52（锚键）/10（最激进键）——**F 统计量日常可达** | A（§2.3） |
| 键覆盖面 | 100/100 键在双支路（small_p + large_p）均有可表示真值被归零的点 | A（§2.4） |
| 上游语义 | statsmodels `_pss_pvalue` 源码同为 `1 - norm.cdf(...)` → 上游同病；修复为数学恒等式，非语义变更 | A（§3，probe [B] 源码实录） |
| 修复式 | `return 0.5 * std::erfc(poly / std::sqrt(2.0));`（一行） | —（§5.1） |
| golden REGEN 面 | **无**——PSS 未暴露于 C-API/Python/MCP（grep 0 命中） | A（§4.3） |
| 回归风险 | 低——`pss_pvalue_bounds` 无直接单测；间接中尾断言容差 ≫ 修复前后差（≤5e-14 rel） | A（§4.2/§2.5） |

---

## 1. 缺陷定位与机制（证据级 A：源码直读）

### 1.1 缺陷点

[pss_critical_values.hpp:70](file:///f:/Cpp_Hub/include/cpphub/timeseries/ardl/pss_critical_values.hpp#L70)，`detail::pss_response_pvalue` 末行：

```cpp
return 1.0 - normal_cdf(poly);
```

其中 `normal_cdf`（:53-55）= `0.5 * std::erfc(-x / std::sqrt(2.0))`，即标准正态 CDF `Φ(poly)`。求值链（:58-71）：

1. `x = std::log(stat)`（对数域）；
2. `stat > PSS_STAT_STAR[i]`（原始 stat 尺度比较）→ small_p（2 阶 Horner，3 系数）否则 large_p（3 阶，4 系数）；
3. `poly = Σ cᵢ·xⁱ` → `p = 1 − Φ(poly)`。

### 1.2 机制

`poly` 为大正值时（深尾小 p），`Φ(poly) = 0.5·erfc(−poly/√2)` 中 `erfc(−大值)` 舍入到 2.0，`normal_cdf` 返回精确 1.0，`1 − 1.0 = 0`。双精度灾难窗口（§2.2 实测）：

- **poly ≥ 8.06**：相对误差 > 10%（垃圾区起点，真值 p ≈ 3.8e-16）；
- **poly ≥ 8.25**：精确归零（真值 p ≈ 7.9e-17，完全可表示）。

### 1.3 求值语义为何此前漏网

Phase 8 冻结语义 probe_pss_cv.py [B] 的数值自证域是 `{stat ∈ [2, 10]}`（中尾），`1 − Φ` 与 `Φ(−poly)` 在该域相对差 ≤1e-16 —— 转录正确性自证无法暴露求值式的深尾缺陷。这与 RISK-036/D-1 的漏网模式一致：**转录验证 ≠ 求值稳健性验证**。

---

## 2. 本地实证（证据级 A：mpmath 60 dps + Python float 逐算子复刻 C++）

方法：系数取 .inc 中 double 字面量原值，`x = log(stat)`、Horner、`erfc` 全部用 IEEE double 逐算子复刻 C++ 行为（Python float = C double）；真值用 mpmath 60 dps 从同一 double 系数出发计算（锚点可复现）。复现脚本见附录 A（tmp 产物验完即删，脚本本体已删，本附录为完整自足复刻）。

### 2.1 锚键 (k=1, case=3, I0) stat 网格（small_p 支路，系数 `{-0.57887, 0.87657, 0.32929}`）

| stat | poly | true p (60 dps) | 旧实现相对误差 | 上尾直算相对误差 | 旧==0.0 |
|---:|---:|---:|---:|---:|:---:|
| 0.5 | −1.131 | 0.8709 | 4.0e-17 | 4.0e-17 | 否 |
| 2 | 0.187 | 0.4259 | 5.3e-17 | 5.3e-17 | 否 |
| 4.923 (0.1% CV) | 1.655 | 0.04898 | 4.1e-16 | 1.1e-17 | 否 |
| 10 | 3.185 | 7.2285e-4 | 2.1e-14 | 6.7e-16 | 否 |
| **100** | **10.441** | **8.027e-26** | — | 5.2e-15 | **是** |
| 1000 | 21.189 | 6.029e-100 | — | 1.2e-14 | **是** |
| 1e4 | 35.428 | 3.117e-275 | — | 1.1e-13 | **是** |
| ≥1e5 | ≥53.2 | <1e-308 | — | 真欠流归零（可接受） | 是 |

**stat=100 的 F 统计量是协整检验的日常量级**（强长期关系 + 中小样本），此区间响应面 p 值被精确归零，而上尾直算保持 5e-15 相对精度。

### 2.2 求值层直扫（poly 自变量，与键/stat 无关）

| poly | true p | 旧实现 | 相对误差 | 上尾直算相对误差 |
|---:|---:|---:|---:|---:|
| 5.0 | 2.867e-07 | 2.867e-07 | 0.00% | 2.2e-15 |
| 7.0 | 1.280e-12 | 1.280e-12 | 0.00% | 1.3e-16 |
| 7.5 | 3.191e-14 | 3.197e-14 | 0.21% | 7.4e-15 |
| 8.00 | 6.221e-16 | 6.661e-16 | 7.08% | 5.7e-15 |
| 8.25 | 7.920e-17 | **0** | 100% | 5.2e-16 |

上尾直算首个归零点 poly≈39（erfc 输入 27.6，真值 <1e-300）——纯欠流，非灾难性抵消，可接受。

### 2.3 锚键 stat 阈值（二次方程精确解）

- poly=8.06（10% 误差起点）→ **stat = 52.52**
- poly=8.25（精确归零）→ **stat = 55.44**

对比前三例：VR 两尾 p 归零需 |z|≥8.3，SMM 联合 p 需 cd≥~6，本例 **stat≥52 即归零**——F 统计量在协整检验中远比 |z|≥8 常见。

### 2.4 全 100 键扫描（缺陷签名：真值可表示 >1e-300 但旧实现返回 0.0）

- **中招键数：100/100**（small_p 支路 stat 阈值范围 [10.1, 56.8]，1e15 内全部归零）；
- **最严重单点**：键 idx=69（k=7, case=5, I1）**large_p 支路**，stat=3.162e-5 时 true p=4.141e-17 被归零（该键 large_p 系数 `{-2.08097, 2.02139, 0.5348, 0.02343}`，c₂x² 项在小 stat 侧把 poly 推到 8.33）；
- 锚键 (1,3,I0) 的 large_p 支路 poly 恒 < −0.75（导数正定），该支路无灾难窗口——**窗口是否存在取决于各键系数，不能按键抽查，必须按族修复**（即全部 100 键共用同一表达式，一行修复自动覆盖双支路）。

### 2.5 中尾等价性（回归风险评估）

锚键 stat 网格，旧 vs 新相对差：

| stat | true p | 旧−新相对差 | 旧误差 / 新误差 |
|---:|---:|---:|---|
| ≤2 | ≥0.426 | **0（逐位一致）** | 同 |
| 4.923 | 0.04898 | 4.3e-16 | 4.1e-16 / 1.1e-17（新更准 25×） |
| 10 | 7.228e-4 | 2.0e-14 | 2.1e-14 / 6.7e-16 |
| 50 | 1.514e-15 | 2.6e-02 | 2.6e-02 / 5.3e-16 |

结论：中尾 **new ≥ old**（逐位一致或更准），最大差异 5e-14 rel——远小于既有测试容差（§4.2）。

---

## 3. 上游与生态核对

### 3.1 statsmodels 上游同病（证据级 A：probe 源码实录）

[probe_pss_cv.py:33](file:///f:/Cpp_Hub/tests/fixtures/timeseries/probe_pss_cv.py#L33) 对 statsmodels 0.14.6 `model.py L2447-2455 _pss_pvalue` 的源码实录：

```python
p_val = 1 - norm.cdf(x @ p)
```

上游未用 scipy 提供的稳健上尾 `norm.sf`。即：**本缺陷系 statsmodels 原生缺陷随系数表一并转录而来**，现有「与 statsmodels ≤1e-16」自证在中尾域依然成立（数学恒等 + 中尾逐位一致），深尾域本库修复后严格优于上游。

### 3.2 修复不构成语义变更

`1 − Φ(x) ≡ Φ(−x) ≡ 0.5·erfc(x/√2)` 为数学恒等式；`.inc` 系数表、stat_star 切换、键域语义全部不动。probe [B] 冻结语义、PSS9/PSS10 排幻觉点均不受影响。

---

## 4. 影响面（证据级 A：grep + 源码直读）

### 4.1 消费方唯一且 p 值为信息性输出

[uecm_bounds.hpp:366](file:///f:/Cpp_Hub/include/cpphub/timeseries/ardl/uecm_bounds.hpp#L366)：`r.p_values = pss_pvalue_bounds(f, k, case);` 之后，三值结论判定（:368-374）用 **stat vs 临界值** 判据（`f > up[ai]` / `f < lo[ai]`），**不消费 p_values**。即缺陷污染的是 `BoundsTestResult.p_values` 报告值（深尾下报 0.0，夸大显著性量级），**结论完整性（LevelRelation/NoLevelRelation/Inconclusive）不受影响** → 定级 P1 而非 P0。

### 4.2 无既有断言冲突

- `tests/unit/timeseries/test_uecm_bounds.cpp` 全部 20 个 TEST 的 grep 中 `pss_pvalue_bounds` **0 命中**（现有测试只覆盖 CV/stat/conclusion，p 值路径无直接单测——修复同时补上该空白）；
- 间接消费者（F_K1_FiveCases_VsStatsmodelsAndR 等）用中尾 stat，`1 − Φ` 误差 ≤2e-14（§2.1 stat=10 行），远小于其对照容差（statsmodels 1e-10 / R ARDL 1e-8 口径，probe [F] 裁决）。

### 4.3 无 golden REGEN 面

`include/cpphub/c_api/**`、`python/src/**`、`src/**` 对 `bounds_f_test|uecm_bounds|pss` 大小写不敏感 grep **0 命中**——PSS 边限检验为 Phase 8 header-only 能力，未接入 C-API/Python/MCP tool 面，**无 `mcp_golden.json` REGEN 需求**。

---

## 5. 修复方案（待批后实施）

### 5.1 修复式（一行）

[pss_critical_values.hpp:70](file:///f:/Cpp_Hub/include/cpphub/timeseries/ardl/pss_critical_values.hpp#L70)：

```cpp
// 深尾上尾直算 (1−Φ 族第 4 例, RISK-036 同源): poly 大正值时 1−Φ(φ) 灾难性
// 抵消归零 (poly≥8.25 → p=0, 实测真值 7.9e-17); Φ(−x) = 0.5·erfc(x/√2) 恒等
// 式直算上尾, 中尾逐位一致或更准 (锚键 stat≤2 diff=0), 深尾保持 ≤1e-13 相对精度。
// 覆盖 small_p/large_p 双支路 (同一表达式); stat≳1e5 时真值 <1e-308 纯欠流归零, 非缺陷。
return 0.5 * std::erfc(poly / std::sqrt(2.0));
```

与 D-1（`2.0 * normal_cdf(-std::fabs(stat))`）、SMM（`-expm1`）同族风格；不改 `normal_cdf` 本身（调用点仅此一处 + 无其他消费者，避免动共享语义面）。

### 5.2 TDD 计划（先红后绿）

**测试落点**：`tests/unit/timeseries/test_uecm_bounds.cpp` 新增 1 个 suite（该文件已直接调用 `ardl::pss_*` 公共接口；**实施时注意 BOM 字节级检查**（EF BB BF，Edit 工具剥 BOM 陷阱，P2 批教训））。

**锚值**（mpmath 60 dps，从 .inc double 系数出发复算——可复现口径，附录 A）：

| 用例 | 键 | stat | 期望 p（锚） | 断言 |
|---|---|---:|---|---|
| 红（深尾） | (1,3,I0) small_p | 100 | 8.0270903435668672315e-26 | rel ≤1e-9 且 **p>0**（旧实现归零即红） |
| 红（更深尾） | (1,3,I0) small_p | 1000 | 6.0291323510070266498e-100 | 同上 |
| 红（large_p 支路） | (7,5,I1) large_p | 3.1622776601683795e-05 | 4.1407709081409948953e-17 | 同上（覆盖第二支路 + key idx 69 最严重单点） |
| 绿守卫（中尾回归） | (1,3,I0) small_p | 10 | 7.228475924500985693e-4 | rel ≤1e-12（旧实现亦过 → 纯回归守卫，非红） |

**容差论证**：double Horner poly 相对误差 ~1e-16 经 `d ln p = −poly·dpoly` 放大（poly≈8-35）→ p 相对误差 ~1e-14 量级 ≪ rtol 1e-9；旧实现 `>1e-24` 无归零守卫**在本例不适用**（最深锚 8e-26 < 1e-24，正确实现也过不了），改用 `p > 0` + rtol 双断言。

**验证矩阵**：Debug 7 目标构建 + Release 全量 ctest + `test_uecm_bounds` 定向；无 golden 面 → 无 REGEN；`.inc` 不动。

### 5.3 回填

- CHECKLIST §15（C-130 修复 / C-131 锚测试双红记录 / C-132 全量验证）；
- DEVELOPMENT_LOG 批条目（**教训**：本例推翻「深尾=理论边角」的族内惯性认知——阈值 stat≈52/10 是现实可达区间；转录验证域 ≠ 求值稳健性域）；
- 记忆更新（project_memory 深尾族洞察补第 4 例 + 阈值现实性教训）。

---

## 6. 边界与不做

- 不改 `.inc` 系数表（statsmodels 0.14.6 转录冻结，缺陷在求值式不在系数）；
- 不做 statsmodels 上游 PR（out of scope）；
- 不加 stat 上限截断/告警（真欠流区间为正确行为）；
- MC 模拟路径 / Narayan 分档表仍按 A5 不做（PSS 排幻觉点 A5 维持）。

---

## 附录 A：复现脚本（自足，≥Python 3.12 + mpmath 1.3.0）

```python
# 证据核心复刻: 逐算子 IEEE double 复刻 C++ 求值 + mpmath 60dps 真值
import math
from mpmath import mp, mpf, erfc, sqrt as msqrt, exp as mexp
mp.dps = 60
C_SMALL = (-0.57887, 0.87657, 0.32929)   # .inc 键 (1,3,I0) small_p (double 原值)
C_LARGE69 = (-2.08097, 2.02139, 0.5348, 0.02343)  # 键 (7,5,I1) large_p
SQ2 = math.sqrt(2.0)

def poly_small(x):  # Horner 升幂, 与 pss_critical_values.hpp:63-65 逐算子同构
    return C_SMALL[0] + x * (C_SMALL[1] + x * C_SMALL[2])

old_p = lambda poly: 1.0 - 0.5 * math.erfc(-poly / SQ2)   # 现实现 :70
new_p = lambda poly: 0.5 * math.erfc(poly / SQ2)          # 修复式 §5.1
true_p = lambda poly: erfc(mpf(poly) / msqrt(2)) / 2      # 60dps 真值

for stat in (10.0, 100.0, 1000.0):
    poly = poly_small(math.log(stat))
    print(stat, float(true_p(poly)), old_p(poly), new_p(poly))
# → 10:   7.228475924500986e-04  7.228475925e-04(ε2e-14)  同级
# → 100:  8.027090343566867e-26  0.0(归零)  8.027090343566867e-26
# → 1000: 6.029132351007027e-100 0.0(归零)  6.029132351007027e-100

# 归零 stat 阈值 (锚键, 二次方程精确解): poly=8.25 → stat=55.44; poly=8.06 → 52.52
c0, c1, c2 = (mpf(v) for v in C_SMALL)
x = (-c1 + msqrt(c1*c1 - 4*c2*(c0 - mpf("8.25")))) / (2*c2)
print(float(mexp(x)))  # → 55.44
```

全 100 键扫描、poly 直扫、中尾等价性表的完整脚本运行于 2026-10-26（输出已录 §2 各表），tmp 脚本按惯例验完即删；本附录含复现全部关键数字所需的系数与逻辑。

## 附录 B：证据清单

| # | 断言 | 级别 | 来源 |
|---|---|---|---|
| E1 | 缺陷点 :70 求值链语义 | A | pss_critical_values.hpp 直读 |
| E2 | 灾难窗口 poly∈[8.06,8.25]、锚键阈值 52.5/55.4 | A | §2.2/§2.3 mpmath 实测 |
| E3 | stat=100 真值 8.027e-26 被归零 | A | §2.1 实测 + 附录 A 复现 |
| E4 | 100/100 键中招、large_p 支路 key69 单点 | A | §2.4 全键扫描 |
| E5 | statsmodels `_pss_pvalue` = `1 - norm.cdf` | A | probe_pss_cv.py:33 源码实录 |
| E6 | 中尾 new≥old（≤2 逐位一致） | A | §2.5 实测 |
| E7 | p_values 信息性、结论走 stat vs CV | A | uecm_bounds.hpp:366-374 直读 |
| E8 | 无 C-API/Python/MCP 暴露 | A | grep 0 命中（2026-10-26） |
| E9 | pss_pvalue_bounds 无直接单测 | A | test_uecm_bounds.cpp grep 0 命中 |
| E10 | 上游 scipy 有 norm.sf 未被用 | B | scipy 文档语义（修复侧写） |
