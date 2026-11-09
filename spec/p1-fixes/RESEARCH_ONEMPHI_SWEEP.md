# `1−Φ` 族全库专项复核（第 5 例扫尾）

> 状态: **调研完成，未改任何源码**（等待批准后实施）
> 上游依据: RISK-036（χ² 上尾）/ SMM 联合 p / D-1（两尾正态 p）/ PSS 响应面（`1−Φ` 族第 1-4 例，均已闭环）
> 关联文档: [RESEARCH_P2_C1128_AND_DEFERRED.md](./RESEARCH_P2_C1128_AND_DEFERRED.md)、[RESEARCH_PSS_DEEPTAIL.md](./RESEARCH_PSS_DEEPTAIL.md)
> 日期: 2026-11-09

---

## §0 结论速览

全库扫描 `1−P` 型上尾写法（C++ `include/cpphub/` 581 处 `1 - …` 全量过滤 + Python 绑定层），**确认 2 处仍在族的缺陷**，其余全部良性或已修：

| # | 位置 | 形式 | 判质 | 建议 |
|---|---|---|---|---|
| **A** | `hfecon/tests/aj_jump_test.hpp:156` + `bns_jump_test.hpp:228` | `normal_cdf(x) = 0.5·(1+erf(x/√2))`，被 `2·normal_cdf(−|z|)` 用作 p 值（AJ / JO / BNS 三处） | **缺陷（新发现）** | P1，最小一行×2 |
| **B** | `pricing/analytic/cev_analytic.hpp:130` | `P1 = 1 − noncentral_chi2_cdf(c; ν₁, d)` | **缺陷（新发现，双层根因）** | P2，需成批 |
| C1 | `monte_carlo/conditional_mc.hpp:79` | `1 − exp(exponent)`（邻接族 `1−e^x`） | 量级无关紧要，形式不稳健 | P3，可选一行 |
| C2 | `cdo.hpp:320` / `basket_cds.hpp:424` | `PD = 1 − exp(−h·t)`（同邻接族） | 实践良性（PD ≥ 1e-4 时误差 ≤1e-14） | 不动 |
| — | 其余 581 处 | — | 良性 / 已修 / 非概率语义 | 不动 |

**结论**：`1−Φ` 族在库内（除上述 A、B 外）已收敛；族第 5 例 = **A（HFE 跳检验）**，B 为同族但属定价精度面。

---

## §1 扫描范围与方法（A：可复现）

### 1.1 扫描模式（ripgrep，路径 `f:/Cpp_Hub/include/cpphub`）

```
1(\.0)?\s*-\s*[A-Za-z_:.()\[\]]*(cdf|Cdf|CDF|gamma|Gamma|erf|Phi|phi|normal|chi|t_dist|lower)
2(\.0)?\s*-\s*2(\.0)?\s*\*
1(\.0)?\s*-\s*(normal_cdf|std_normal|phi|Phi|ncchi2|beta_i|F_cdf|t_cdf|...)
(1(\.0)?\s*-\s*(std::)?(erf|erfc)|-expm1|log1p\(-)
0\.5\s*\*\s*\(\s*1(\.0)?\s*[+-]\s*(std::)?erf      ← 关键: 1−Φ 的"内联"形态
-\s*2(\.0)?\s*\*\s*(std::)?(normal_cdf|erfc)
1(\.0)?\s*-\s*[A-Za-z_][A-Za-z0-9_:.]*\(           ← 581 处全量，逐条过滤
```

### 1.2 覆盖面

- `include/cpphub/` 全 100 个文件（581 处 `1 - 标识符` 命中全部过目）
- `python/` 绑定层（仅封装 + 测试，无独立 p 值计算）
- `src/`（C API 实现，仅透传）
- `include/cpphub/c_api/cpphub_tool_registry.hpp`（MCP 信封组装，无计算）

### 1.3 判质口径

- **深尾风险**：被减项是「上尾概率」而写成 `1 − 下尾 CDF`，或 `1 − 生存函数` 且结果用于小概率输出 → 灾难性抵消
- **良性**：`1 − P` 中 P 接近 1（结果近 0 但减法无损）、或 `1 − P` 是「下尾 CDF = 1 − 上尾 sf」（结果近 1，相对误差无意义）、或纯参数运算（`1 − β`、`1 − φ`、`1 − α`）

---

## §2 清单：全部 `1−P` 型判质（B，逐条）

### 2.1 良性 / 已修（42 处）

| 位置 | 形式 | 理由 |
|---|---|---|
| `core/math.hpp:22` | `normal_cdf = 0.5·erfc(−x/√2)` | ✅ 上尾直算基线（erfc） |
| `core/math.hpp:147` | `P = 1 − gcf` (x≥a+1) | P≈1，减法无损 |
| `core/math.hpp:162` | `Q = 1 − P` (x<a+1) | 该支路 Q≥0.083，相对误差 ≤3e-15（RISK-036 已论证） |
| `core/math.hpp:172` `chi2_sf` | 上尾直算 | ✅ RISK-036 |
| `econometrics/core/special_functions.hpp:151-170` | `gammp = 1−gcf` / `gammq = 1−gser` | NR 双支路，区间互补且各自安全 |
| `special_functions.hpp:192` `chi2_cdf` | `1 − chi2_sf` | 下尾，结果近 1，无风险 |
| `unit_root_common.hpp:173` | erfc 型 | ✅ |
| `estimation/mle.hpp:119` | erfc 型 | ✅ |
| `resampling/bootstrap_base.hpp:238` | erfc 型 | ✅ |
| `unit_root/variance_ratio_test.hpp:154,164` | `2·normal_cdf(−|z|)` + `−expm1(m·log1p(−q))` | ✅ D-1 / SMM 批已修 |
| `unit_root/mackinnon_cv.hpp:200-207`、`cointegration/mackinnon_coint_cv.hpp:84-92` | `p = Φ(polyval(τ))` | 直算 Φ，非 `1−Φ`；erfc 型覆盖深尾 |
| `econometrics/inference/*` 全部 p 值 | `chi2_sf` / `beta_i` / `erfc` / `exp(−c·supLR)` | ✅ 无 `1−F` 形态 |
| `risk/var/risk_diagnostics.hpp:473`、`specification_tests.hpp:262,397` | `beta_i(...)` 直算 | ✅ F/t 上尾恒等式直算 |
| `structural_break.hpp:282` | `1 − significance_level` | 置信水平，非概率 |
| 其余 ~30 处 | `1 − β`、`1 − φ`、`1 − p`（参数） | 参数运算 |

### 2.2 有意保留（AD 专用，C）

`risk/greeks/ad_tape.hpp:14`、`risk/greeks/ad_dual.hpp:17`：

```cpp
return var(0.5) * (var(1.0) + erf(x / sqrt(var(2.0))));   // 非 erfc 型
```

**不改**：autodiff reverse/forward 提供 `erf` 的导数，`erfc` 无原生导数支持；用途为 BSM Greeks，深 OTM 时 delta~1e-16（远低于任何定价容差），不产生可观察误差。若未来 AD 代码需要深尾，须自行给 `erfc` 定义导数。

### 2.3 测试内辅助实现（不在库面，记录备查）

- `tests/unit/timeseries/test_egarch_model.cpp:345`、`tests/unit/econometrics/test_integration_phase7a.cpp:237`：erf 型 normal_cdf，但只用于中尾期望值/似然，无深尾断言 → 本批不动。

---

## §3 缺陷 A：HFE 跳检验 p 值深尾归零（族第 5 例）

### 3.1 定位（A）

```cpp
// include/cpphub/hfecon/tests/aj_jump_test.hpp:155-157
inline Real normal_cdf(Real x) noexcept {
    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));   // ← 1−Φ 抵消
}
// 调用点
aj_jump_test.hpp:401   pvalue = 2.0 * detail::normal_cdf(-std::fabs(AJtest));
jo_jump_test.hpp:308   pvalue = 2.0 * detail::normal_cdf(-std::fabs(JOtest));  // 复用 AJ 的 detail
// include/cpphub/hfecon/tests/bns_jump_test.hpp:227-229 同一实现
bns_jump_test.hpp:151  p_value = 2.0 * normal_cdf(-std::fabs(z));
```

机制与族前四例同源：`erf(x)` 在 `x < −5.86` 时舍入为 `−1.0`，`0.5·(1+(−1))` 精确得 0 → p 值归零；真值仍在 double 可表示范围。

### 3.2 实证（A：mpmath 60dps + IEEE double 仿真，脚本见附录 A）

| \|z\| | 真值 `2Φ(−z)` | 现状输出 | 相对误差 |
|---|---|---|---|
| 4.6675（现有锚） | 3.048909412439e-06 | 3.048909412384e-06 | 1.80e-11 |
| 6.0 | 1.973175290075e-09 | 1.973175289827e-09 | 1.26e-10 |
| 7.0 | 2.559625087772e-12 | 2.559619183273e-12 | **2.31e-06** |
| 7.5 | 6.381783345822e-14 | 6.383782391595e-14 | **3.13e-04** |
| 8.0 | 1.244192114854e-15 | 1.221245327088e-15 | **1.84e-02** |
| 8.2 | 2.403870308547e-16 | 2.220446049250e-16 | **7.63e-02** |
| **8.25901**（二分阈值） | 1.469e-16 | **0** | 归零 |
| 9.0 | 2.257176811908e-19 | **0** | 归零 |
| 10.0 | 1.523970604832e-23 | **0** | 归零 |
| 20.0 | 5.507248237212e-89 | **0** | 归零 |

- **归零阈值 \|z\| = 8.2590072043**（二分定位）——此前已有 1.8%（\|z\|=8）/ 7.6%（\|z\|=8.2）误差，非突变
- 对照：erfc 型 `0.5·erfc(−x/√2)` 归零点在 **\|z\| = 38.4749**（真欠流，双尾 p < 4.9e-324）
- 现有基线锚 \|z\| = 0.6927 / 4.6675 均落在安全区，修复对现有测试的影响 ≤1.3e-16 绝对

### 3.3 影响面（B）

- 3 个 p 值出口：`BNSJumpTest` / `AJJumpTest` / `JOJumpTest`（`result.p_value` / `.pvalue`）
- 无 C-API / MCP / golden 面（`cpphub_tool_registry.hpp` 与 `python/tests/fixtures/mcp_golden.json` 均 0 命中）→ 纯 C++，无 REGEN
- **可实现性**：AJ 统计量 = 标准化收益，单笔 8σ 跳跃即 \|AJ\|≈8；`tests/unit/hfecon/test_aj_jump_test.cpp:181-185` 已自述"稀疏子采样使 \|z\| 极大, pvalue 可能下溢为 0"——**该注释把抵消误标为"浮点正常行为"，需一并更正**
- 决策影响有限（p 归零与 p=2e-19 都判"显著"），但属于族内已定纪律："深尾不得静默归零"

### 3.4 修复方案（A：一行 × 2，非 3）

```cpp
// aj_jump_test.hpp:155-157（JO 复用同一函数）
inline Real normal_cdf(Real x) noexcept {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}
// bns_jump_test.hpp:227-229 同改
```

- `jo_jump_test.hpp` 复用 AJ 的 `detail::normal_cdf`，无需单独改
- 中尾等价：\|z\|≤2 处与旧式逐位一致或 1 ULP 内
- `qnorm` 的 Newton 精炼（`normal_cdf(x) − p`）同步受益，alpha=0.001 处无回归
- 现有基线（z=0.6927 / 4.6675）变化 ≤1.3e-16 绝对，`TOL_STANDARD=1e-10` 与 p 容差 1e-8 均不受影响

### 3.5 TDD 计划

1. **红**：新增 `BnsJumpTest.PvalueDeepTailNoCollapse`——构造 z≈9（或直接调 `detail::normal_cdf`）断言 `p > 0` 且与 mpmath 锚 `2.257176811908e-19` 相对 1e-9 相符；旧实现 `p == 0` 必红
2. **绿**：改 `normal_cdf` 两处
3. **守卫**：`p(4.6675)` = 3.048909412439e-06 相对 1e-12 不破；AJ/JO 各加一条深尾锚
4. **更正注释**：`test_aj_jump_test.cpp:181-185` 的"下溢为 0 是浮点正常行为"改为指向本条修复

---

## §4 缺陷 B：CEV `P1 = 1 − 非中心 χ² CDF`（同族，双层根因）

### 4.1 定位（A）

```cpp
// include/cpphub/pricing/analytic/cev_analytic.hpp:123-131
// P_1 = 1 - χ'²(c; nu_1, d)  (生存概率, S_T > K)
Real chi2_1 = noncentral_chi2_cdf(c, nu_1, d);
Real P1 = 1.0 - chi2_1;
```

**双层根因**：

1. **`1 − F` 抵消**（族本体）：`F → 1` 时 P1 的相对精度崩塌
2. **`noncentral_chi2_cdf` 大 λ 分支截断偏差**（放大根因）：`core/math.hpp:226-228` 的 Poisson 窗口取 `λ/2 ± 5√(λ/2)`，只覆盖 99.99994% 质量 → **CDF 系统性偏低约 2e-7（绝对）**

即使只修 1，P1 仍被 2 的 2e-7 偏差主导（深 OTM 处形成 **3.84e-7 的 P1 地板**）。

### 4.2 实证（A：忠实镜像 C++ 算法 + mpmath 30σ 窗口真值，双路交叉验证）

参数：`S=100, T=1, r=0.05, q=0, β=0.5, σ=0.2`（ν₁=4, ν₂=2；典型标定参数）

| K | c | d(=λ) | P1 真值 | P1 现状 | 现状/真值 | 价格真值 | 价格现状 |
|---|---|---|---|---|---|---|---|
| 100 | 9.752e3 | 1.025e4 | 9.940473e-01 | 9.940475e-01 | 1.000 | 4.88097 | 4.88099（+1.9e-5） |
| 110 | 1.073e4 | 1.025e4 | 1.057153e-02 | 1.057174e-02 | 1.000 | 6.896e-03 | 6.917e-03（+2.1e-5） |
| 120 | 1.170e4 | 1.025e4 | **2.402992e-12** | **3.837701e-07** | **1.6e+05** | 6.16e-13 | 3.84e-05 |
| 150 | 1.463e4 | 1.025e4 | 1.58e-86（sf 交叉） | 3.837687e-07 | 天文 | ~0 | 3.84e-05 |
| 200 | 1.950e4 | 1.025e4 | 8.9e-323（近欠流） | 3.837687e-07 | 天文 | ~0 | 3.84e-05 |

- **大 λ 分支偏差直接测量**：K=100 时 CDF 真值 0.005952707409 vs 现状 0.005952512972 → **偏差 −1.944e-07**；K=110 → **−2.081e-07**（与 5σ 窗口漏掉的 Poisson 尾质量同量级）
- 真值双路自证：`1 − CDF(30σ 窗口)` 与 `SF(30σ 窗口)` 在 K=100/120 相对差 ≤2e-35
- **P1 地板 3.8377e-07 ⇒ 价格地板 3.84e-05**，与 K 无关（截断偏差的常数特征）

### 4.3 影响面（B）

- 唯一调用点（grep 确认 `noncentral_chi2_cdf` 仅 `cev_analytic.hpp` 使用）
- **无外部基线**：`tests/unit/models/test_cev.cpp` 全部为结构断言（parity / β=1→BSM / 单调性 / IV skew / 边界），**无绝对价格锚** → 该偏差对现有测试不可见
- 实际后果：① ATM 价格绝对偏差 ~2e-5（相对 ~4e-6）；② 深 OTM 翼 P1 触底 → **CEV 隐含波动率翼部被人为垫高**（真值应趋于 0），影响 `CEVCalibrator` 对远虚值报价的拟合形状
- 判定 **P2**（非 P1）：无决策规则依赖，绝对误差 3.8e-5；但数值上确实失效

### 4.4 修复方案（B）

**B1（族本体）**：新增 `core/math.hpp::noncentral_chi2_sf(x, k, λ)` —— 与 `noncentral_chi2_cdf` 完全镜像，把 Poisson 混合核从 `regularized_lower_gamma` 换成 `regularized_upper_gamma`（RISK-036 已保证后者上尾直算），`P1 = noncentral_chi2_sf(c, ν₁, d)`。

**B2（配套，必需）**：大 λ 分支窗口改为自适应——
- 下界：`j* − 10σ`（10σ 漏质量 ~1e-23，修掉 CDF 的 2e-7 偏差）
- 上界：向上累加至 `w_j < 1e-17·sum`（相对判据，尺度无关）或 `log w_j < −745`（欠流）
- 依据：`Q(a,x) ≤ 1`，故尾部漏项 ≤ Σw ≤ 1e-17·sum → 相对误差 ≤1e-17

**B3（验证）**：新增 mpmath 锚（K=100/110/120 的 P1 与价格，rtol 1e-9）；CEV 结构测试保持绿；因无外部基线，锚值以 mpmath 60dps 为准。

**风险**：改变 CEV 全体价格（更准）；`noncentral_chi2_cdf` 仅此一处消费，面可控；无 golden/API 影响。

---

## §5 邻接族与不动项

### 5.1 邻接族 `1 − e^x`（C：形式同源，非 `1−Φ`）

| 位置 | 形式 | 实测相对误差（x→0） | 建议 |
|---|---|---|---|
| `monte_carlo/conditional_mc.hpp:79` | `survival_prob = 1 − exp(exponent)` | x=1e-8 → 1.0e-8；x=1e-12 → 2.2e-5；x=1e-16 → 0.11 | 改 `−std::expm1(exponent)`（一行） |
| `instruments/credit/cdo.hpp:320`、`basket_cds.hpp:424` | `PD = 1 − exp(−h·t)` | x=1e-4 → 1.3e-14；x=1e-8 → 1.0e-8 | 可选同改；实践 PD ≥ 1e-4 → 不动可接受 |
| `models/ir/short_rate.hpp`、`instruments/ir/swaption.hpp` | `(1 − exp(−κτ))/κ` 等 | κτ ~1e-4 时 ~1e-14 | **不动**（物理参数下无感） |

### 5.2 明确不动

- AD 专用 erf 型 `normal_cdf`（§2.2）
- 测试内辅助 erf 型实现（§2.3）
- `path_dependent_analytic.hpp:437,514` 的 `Φ(−d1_ref) − Φ(−d1)` 型：属"两概率相减"而非 `1−Φ`，且乘在小修正项上——**不在本族范围**，如需另立课题

---

## §6 边界与不做

1. 本次复核**未修改任何源码**（仅新增本文档与 tmp 验证脚本，脚本验完即删）
2. 不覆盖 `tests/fixtures/**.py` / `verify_*.R`：属基准生成侧（用 scipy/R，本身准确）
3. 不覆盖 `1 − P` 型**下尾**补（`chi2_cdf = 1 − chi2_sf` 等）：结果近 1，深尾无风险
4. 缺陷 B 若实施，**必须**同时做 B1+B2；只做 B1 会把 P1 从 3.8e-7 降到 ~2e-7（仍错），无意义

---

## §7 实施建议（待批）

| 批 | 内容 | 预估改动 | 验证 |
|---|---|---|---|
| **批 1（建议先做）** | 缺陷 A：`aj_jump_test.hpp` + `bns_jump_test.hpp` 两处 `normal_cdf` 改 erfc 型；更正 `test_aj_jump_test.cpp:181-185` 注释；新增 3 条深尾 TDD 锚 | 2 行 + 测试 | Debug/Release 定向套件（hfecon jump 测试）+ 全量 CI |
| **批 2（可选）** | 缺陷 B：`noncentral_chi2_sf` + 大 λ 自适应窗口 + CEV P1 切换 + mpmath 锚 | ~60 行 + 测试 | CEV 套件 + Calibrator 套件 + 全量 CI |
| **批 3（可选）** | C1 `conditional_mc` 改 `expm1` | 1 行 | 障碍期权 MC 套件 |

批 1 与批 2 相互独立，可分别闭环；若只做批 1，族第 5 例即闭环。

---

## 附录 A 自足复现脚本

```python
# 需 mpmath; 用 Python312 运行
import math, mpmath as mp
mp.mp.dps = 60

def ncdf_erf(x):   # 现状 (aj/bns 局部实现)
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))
def ncdf_erfc(x):  # 稳健式
    return 0.5 * math.erfc(-x / math.sqrt(2.0))
def truth(x):
    return 0.5 * mp.erfc(-mp.mpf(x) / mp.sqrt(2))

for z in [4.6674971744735663, 7.0, 8.0, 8.259, 9.0, 10.0, 20.0]:
    t, o, n = 2*truth(-z), 2*ncdf_erf(-z), 2*ncdf_erfc(-z)
    print(f"|z|={z:8.4f} 真值={float(t):.12e} 现状={o:.12e} 修复后={n:.12e}")

# 归零阈值二分
lo, hi = 1.0, 12.0
for _ in range(80):
    mid = 0.5*(lo+hi)
    if ncdf_erf(-mid) == 0.0: hi = mid
    else: lo = mid
print("erf 型归零阈值 |z| =", hi)      # 8.2590072043
lo, hi = 1.0, 60.0
for _ in range(80):
    mid = 0.5*(lo+hi)
    if ncdf_erfc(-mid) == 0.0: hi = mid
    else: lo = mid
print("erfc 型归零阈值 |z| =", hi)     # 38.4749 (真欠流)
```

CEV 部分（`1 − F` + 大 λ 截断）的镜像脚本较长，核心见本文 §4.2 数值表：镜像 C++ `noncentral_chi2_cdf`（`gser`/`gcf` + 5σ 窗口）后与 mpmath 30σ 窗口真值比对，即可复现 `−1.944e-07`（K=100）与 P1 地板 `3.8377e-07`。

## 附录 B 证据清单

| ID | 断言 | 级别 | 证据 |
|---|---|---|---|
| E1 | `aj/bns` 局部 `normal_cdf` 为 erf 型 | A | 源码 `aj_jump_test.hpp:156`、`bns_jump_test.hpp:228` |
| E2 | erf 型 p 值在 \|z\|≥8.259 精确归零 | A | 附录 A 二分 + §3.2 表 |
| E3 | \|z\|=8 已有 1.84% 相对误差 | A | §3.2 表（double 仿真 vs mpmath 60dps） |
| E4 | 现有基线锚在安全区（≤1.3e-16 绝对变化） | A | §3.2 末行 |
| E5 | 无 MCP/golden 面 | A | `cpphub_tool_registry.hpp` 与 `mcp_golden.json` grep 0 命中 |
| E6 | `noncentral_chi2_cdf` 仅 CEV 使用 | A | grep 唯一调用点 `cev_analytic.hpp:127-128` |
| E7 | 大 λ 分支 CDF 偏差 ≈ −2e-7 | A | 镜像仿真 vs mpmath（K=100: −1.944e-07；K=110: −2.081e-07） |
| E8 | P1 地板 3.8377e-07（与 K 无关） | A | §4.2 表 K=120/150/200 同值 |
| E9 | CEV 现有测试全为结构断言，无绝对锚 | A | `tests/unit/models/test_cev.cpp:76-184` |
| E10 | `1−e^x` 相对误差曲线 | A | §5.1（double 仿真 vs mpmath） |
| E11 | AD 的 erf 型属有意选择 | C | autodiff 导数支持约束 + 深尾 delta~1e-16 无观察效应（判断） |
| E12 | 批 1 独立即可闭环族第 5 例 | C | 判质推理（§7） |

---

## §8 实施期发现与分批裁定（2026-11-09）

批 1 + 批 3 已实施闭环（提交 `b288573`）；批 2（CEV）在实施中触发更深发现，整批延后专项。

### 8.1 批 1（缺陷 A，族第 5 例）+ 批 3（邻接族 1−e^x）已闭环

- **改动**: `aj_jump_test.hpp:160` / `bns_jump_test.hpp:229-231` 的 `normal_cdf` → `0.5*erfc(-x/√2)`（bns 同时把该函数提到 `public` 供测试直调）；`conditional_mc.hpp:81` → `-std::expm1(exponent)`；更正 `test_aj_jump_test.cpp:181-185` 旧注释；新增 3 条 TDD 锚（`AJJumpTest.NormalCdfDeepTailNoCollapse` / `HFE_BNSJumpTest.PvalueDeepTailNoCollapse` / `ConditionalMC.SurvivalProbDeepSmallExponent`）。
- **TDD 红证据**（临时回退修复点实测）: `normal_cdf(-9.0)` 精确 0（真值 2.2571768119076813e-19）；`|z|=8.3` → 0（真值 1.0411139489780508e-16）；中尾 `|z|=4.6674971744735663` 相对误差 **1.8004264745741239e-11**（>1e-12 容差）；CMC `|surv−ref|/ref` = **1.0989993395467926e-6**（>1e-12）。
- **验证**: 定向套件 AJ 9/9、HFE 19/19、CMC 18/18；Release 全量 **2635/2635 passed**（265.28s）。

### 8.2 批 2（缺陷 B，CEV）延后裁定

实施批 2 时 6 个既有测试失败（`git stash` 基线确认 test_cev 21/21、test_multi_calibrator 15/15 绿 → 均为修复引入）。诊断分两层：

1. **§4.4 / E9 预测被推翻**: 现有测试虽为结构断言（`EXPECT_GT(C,0)`、IV 单调性等），但其红绿**依赖旧地板**——4 个测试（`CEVAnalyticTest.BetaZeroIsBachelierLike` / `CEVAnalyticTest.BetaBelowOneProducesSkew` / `CEVCalibratorTest.BetaBelowOneProducesSkew` / `CEVCalibratorTest.BoundaryBetaNoCrash`）在 P1 修复后语义改变：深 ITM vega→0 时反推 IV 落到搜索下界 0.001、skew 度量被饱和 ITM 腿主导、β=0.10 时 K=120 真价 ~1e-440 无法反推 IV（该测试 57.7s → 154.2s）。修复必须含测试重定基，不是容差调整。
2. **更深的既有缺陷（与 `1−Φ` 无关）**: β<0.5 时 CEV 解析价已错误——β=0.3、σ=4、K=80、r=q=0 解析价 **19.7271**（< 内在价值 20），细网格 MC（6400 步 × 20 万路径）真值 **20.7129±0.033**，差 0.986 ≈ **30 倍标准误**；同参数 β=0.5 时解析 27.0997 vs MC 27.1136（差 0.014 < MC 误差）→ **仅 β<0.5 失效**。自由度指派假说（df₁ = 2+2/(1−β) → 23.14）已排除（离 MC 更远）。IV 曲面探针：β=0.5/σ=4 时 IV 0.4235→0.3826（K=80→120，负偏斜合理）；β=0.3/σ=0.3 时 IV ~0.007 完全不可用。

**裁定**: 批 2 需先调研 β<0.5 的 Schroder 公式适用性（c/d 指派、小 β 边界、吸收壁概率质量），再连同测试重定基做独立专项批；本文档原 §7「批 2 可选」升级为「批 2 前置调研先行」。

### 8.3 对本轮扫描结论的更正

- E9 表述仍成立（无绝对锚），但「结构断言 ⇒ 与 P1 修复无关」的隐含推论错误——结构断言可依赖错误实现产生的非零地板；族修复的回归面必须按**调用链下游测试**枚举，而非按断言类型判断。
- 其余扫描结论（42 处良性/已修、AD erf 型有意保留）不受影响。
