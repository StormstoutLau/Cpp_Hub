# Phase 8 执行规格书 - ARDL 边限检验与 Kalman 状态空间模块 (v1.8 M0-M3)

> **版本**: v1.1 (2026-08-21: probe_pss_cv.py 前置完成 — 5 项修订: PSS10 新幻觉点 (statsmodels bounds_test 查表 k 偏移 +1, 本库裁决用 PSS 语义) / k=0 键不存在修正 / t 检验 p 值 NaN 化 + has_p 字段 / 转录量 1200 数字 / p 值算法冻结; v1.0 冻结)
> **版本归属**: **v1.8** (Phase 8)
> **目标**: 实现 v1.8 双免费午餐 — M1 ARDL/PSS 边限检验 + M2 Kalman/statespace 精确似然 (兑现 ADR-019 决策 1 "Kalman 留 v1.8"), 补齐 Research OS 时序推断链的水平关系检验与状态空间基础设施
> **覆盖范围**:
>   - M1: ARDL(p,q) 单方程 OLS + IC 网格阶数选择 + UECM 重参数化 + PSS bounds F/t 检验 (五 Case) + 响应面临界值
>   - M2: StateSpaceModel (DK 记号) + 常规 Kalman 滤波 + 平稳/approximate diffuse 初始化 + ARMA Harvey 形 + ARIMA 嵌入差分精确似然 MLE + innovations 常驻互检
>   - M3: 战略选项 (D-2 待裁决, 占位不阻塞)
> **前置**:
>   - Phase 7C (v1.7.1) 已闭环: 2447 测试三平台 + CI run #72/#73 全绿 (tag v1.7.1 @ 9771dd6)
>   - 可复用: SLSQP (ADR-018) + ArimaSpec/ArimaConfig (arima_model.hpp) + innovations_mle (互检对手) + cpphub_timeseries_mat 载体 + mackinnon_tables .inc 管线惯例 + residual_diagnostics (7A)
> **里程碑顺序**: M1 (纯复用, 先行热身) ∥ M2 (Eigen statespace 载体, probe 前置) → 收尾; M3 (D-2 裁决后排期)
>
> **Scope 声明** (ADR-020 冻结, 12 项决策):
> - **M1 ARDL/PSS (A1-A7)**: CV 走 statsmodels 响应面多项式路线 (PSS 原表 static_assert 锚); 模型限 ARDL(p,q) 单方程 (ardlverse 系全排除, CV 错误实锤); 阶数 IC 网格 (glob 不做); t 检验仅 Case I/III/V; 小样本靠响应面按 T 插值 (不做 Narayan 分档 / 不做 MC 模拟); UECM 含同期 Δx (causal v1.9); 仅 cov_type='nonrobust'
> - **M2 Kalman (K1-K5)**: SSModel + 常规 KF (DK Ch4) + 平稳初始化 (Lyapunov + 守卫) + approximate diffuse (κ=1e6) + ARMA Harvey 形 + 嵌入差分 ARIMA 精确似然 MLE + innovations 常驻互检; **推迟 v1.9+**: exact diffuse (de Jong) / 平滑器 RTS / 缺失数据 / 多元 / Chandrasekhar / sqrt filter / univariate filter; 优化器复用 ADR-018; σ² 非集中化 (loglik 主锚); Eigen 载体 `timeseries/statespace/` 挂 cpphub_timeseries_mat
> - **对照库版本冻结**: statsmodels 0.14.6 (ARDL + statespace 双主锚) / R ARDL 0.2.5 (JAE 验证, F 渐近交叉) / **禁用** ardlverse (CV 错误实锤)
> - **不做清单**: panel/quantile/Fourier/NARDL / glob 搜索 / Narayan 分档表 / MC 模拟路径 (asymptotic=False) / HC/HAC Wald / causal ARDL / exact diffuse / 平滑器 / 缺失数据 / 多元 statespace / SARIMA 季节
>
> **关联**:
> - [ADR-020: v1.8 ARDL/PSS + Kalman 实施边界 (12 项)](../../decisions/ADR-020_V18_ARDL_KALMAN_BOUNDARY.md) (Accepted 2026-08-21)
> - [PHASE8_RESEARCH.md](../../research/PHASE8_RESEARCH.md) v1.1 (双域调研 + 审计二轮, 3 处修正)
> - [PHASE8_SCOPE_RESEARCH.md](../../research/PHASE8_SCOPE_RESEARCH.md) v1.3 (D-1 裁决: 双免费午餐)
> - [ADR-017](../../decisions/ADR-017_TIMESERIES_NAMESPACE.md) (命名空间) / [ADR-018](../../decisions/ADR-018_SLSQP_BOUNDARY.md) (K3 复用) / [ADR-019](../../decisions/ADR-019_V17_TIMESERIES_BOUNDARY.md) (决策 1/20 本期兑现)
> - [PHASE7C_SPEC.md](../phase7/PHASE7C_SPEC.md) v1.2 (v1.7 复用基础与接口惯例来源)

---

## 0. 门禁记录 (R 清零 → spec 冻结)

| 门禁 | 状态 | 记录 |
|------|------|------|
| R1 断言统计 | ✅ (人工清点) | 域 A 9 + 域 B 10 幻觉点 + 12 项决策依据断言; 无机读 ```assertions``` 块 — **偏离 7C 全量 pilot 模式**, 理由: 断言规模 (Phase 8 ≈ 40 条 vs 7C 126 条) 且全部经主线逐条复验而非抽样 |
| R2 A 类证据可核验 | ✅ (3 处修正) | 审计二轮: ① Chandrasekhar v0.12.0 三方一手实证 (release notes + GitHub tag + 作者博客/PR #6411) ② statsmodels nsim 文档 bug 原文实锤 (源码为准) ③ Narayan k≤7 降级 B- |
| R3 B 类登记 | ✅ | Narayan k≤7 以 B- 登记并直接支撑 A5 (不采 Narayan 分档); Statalist Kripfganz 帖子作 PSS Case III k=1 锚点旁证 (F 4.94/5.73, t −2.86/−3.22) |
| R4 阻断性清零 | ✅ | 12 项决策依据无 FALSIFIED/CONFLICT; 10 项一手来源; 残余 B- 断言转 probe 前置任务 — **M1 probe_pss_cv.py 已完成 (2026-08-21, 5 项语义冻结 + 1 项新发现 PSS10 裁决落档)**; M2 probe_sarimax.py 待执行 |

**本 spec 引用的全部阻断性断言均经审计二轮复验**; 实施期发现冲突时回溯 ADR-020 修订流程, 不得静默改公式。两处 probe 前置任务为 7C "Johansen 双库 diff 前置" 同款制度。

---

## 1. 交付物清单

### 1.1 新增文件结构 (按 ADR-017 命名空间组织)

```
# M1: ARDL/PSS (Eigen3 → cpphub_timeseries_mat, 与 var/cointegration 同载体)
include/cpphub/timeseries/ardl/
├── ardl_model.hpp             # 新增 M1: ARDL(p,q) OLS + trend 四式 + hold_back
├── ardl_select.hpp            # 新增 M1: IC 网格阶数选择 (同一 hold_back 样本)
├── uecm_bounds.hpp            # 新增 M1: UECM 重参数化 + bounds F/t 检验 (五 Case, 三值结论)
├── pss_critical_values.hpp    # 新增 M1: PSS 响应面临界值 (查表 + p 值 + static_assert 锚)
└── pss_response_surface.inc   # 新增 M1: 系数表运行时 .inc (mackinnon_tables.inc 同款惯例,
                               #   生成器产出, 被 pss_critical_values.hpp 于 detail 内 include;
                               #   v1.1 勘误: 原列 tests/fixtures, 实际按 7C 运行时表惯例落 include)

# M2: Kalman/statespace (Eigen3 → cpphub_timeseries_mat, K5)
include/cpphub/timeseries/statespace/
├── ss_model.hpp               # 新增 M2: StateSpaceModel (DK 记号) + 常规 KF + 双初始化 + 守卫
└── arima_statespace.hpp       # 新增 M2: ARMA Harvey 形构造 + ARIMA 嵌入差分 + 精确似然 MLE

# verify 脚本 (沿 C6 模式: probe_*.py 源码探明 + verify_*.py/.R 基准 + gen_* 生成器; 实际入库)
tests/fixtures/timeseries/
├── probe_pss_cv.py            # M1 前置 ✅: 布局/求值/k 语义/t 存在性五项冻结 + [G] .inc 全量 diff
├── gen_phase8_pss_inc.py      # M1 Step 1 ✅: statsmodels → pss_response_surface.inc (1200 数字)
├── probe_sarimax.py           # M2 前置: loglikelihood_burn 默认公式 + 嵌入差分 T/R/Z 布局 + 变换语义探明
├── verify_ardl.py             # statsmodels ARDL.fit 系数/IC + ardl_select_order
├── verify_bounds.py           # statsmodels UECM.from_ardl + bounds_test (stat 口径; CV/p 走 R+表键直查, PSS10)
├── verify_bounds.R            # R ARDL 0.2.5 bounds_f_test/bounds_t_test 交叉 (JAE 验证包)
├── verify_statespace.py       # statsmodels SARIMAX llf/params (method='statespace')
└── gen_phase8_inc.py          # 黄金值生成 (sarimax_baseline.inc 等)
```

### 1.2 新增测试套件

| 测试套件 | 用例数 | 覆盖模块 | 里程碑 |
|----------|--------|----------|--------|
| `test_ardl_model` | 16 | 设计矩阵列序 / 系数 / trend 四式 / hold_back / IC / causal=true 拒绝 | M1 |
| `test_ardl_select` | 10 | IC 网格 / 同一样本 / argmin / 网格维度回显 | M1 |
| `test_uecm_bounds` | 18 | 重参数化 / 五 Case F 约束数 / t 仅三 Case / Case II/IV t 报错 / 三值结论 / t 方向 / CV 锚 | M1 |
| `test_statespace` | 16 | KF 递推黄金锚 / ARMA(1,1) 三矩阵 / Lyapunov + 守卫 / diffuse κ / 似然常数 / 守卫触发 | M2 |
| `test_arima_statespace` | 14 | 嵌入差分 / MLE / SARIMAX 对照 / innovations 互检 / 双口径似然 / 变换域 | M2 |
| `test_integration_phase8` | 6 | 端到端 6 场景 (§5) | 收尾 |

**新增测试总数**: ~80 (累计 ~2527)

### 1.3 必须达到的数值基准 (容差分层)

| 基准 | 容差 | 验收方式 |
|------|------|----------|
| ARDL 系数/σ̂²/IC vs statsmodels ARDL.fit | **1e-10** | `verify_ardl.py` → `test_ardl_model.cpp` |
| ardl_select_order 最优阶 + IC 值 vs statsmodels | **1e-10** | `verify_ardl.py` |
| UECM 重参数化系数 vs statsmodels UECM.from_ardl | **1e-10** | `verify_bounds.py` |
| bounds F 统计量 vs statsmodels bounds_test (stat 不受 k 偏移影响, PSS10) | **1e-10** | `verify_bounds.py` |
| bounds F 统计量 + CV + p vs R ARDL bounds_f_test (exact=FALSE, k 同口径) | 1e-8 | `verify_bounds.R` |
| bounds t 统计量 + CV vs R ARDL bounds_t_test (statsmodels 无 t, probe [D]) | 1e-8 | `verify_bounds.R` |
| CV 查表 vs statsmodels crit_vals 表键直查 (k = force 数) | **精确相等** | `probe_pss_cv.py` |
| p 值 vs statsmodels _pss_pvalue (stat 跨 stat_star 两侧) | 1e-14 | `probe_pss_cv.py` |
| 响应面系数转录 vs statsmodels pss_critical_values.py (1200 数字) | **精确相等** (全量 diff) | `probe_pss_cv.py` + `pss_response_surface.inc` |
| 响应面锚: crit_vals[(1,3,False)][1]=4.9228967 / [(1,3,True)][1]=5.7217442 | **精确相等** | static_assert (编译期) |
| PSS 原表交叉印证 (Case III k=1 5%: F 4.94/5.73, t −2.86/−3.22) | ~0.05 (不同模拟批次, 非转录断言, 文档声明) | `test_uecm_bounds.cpp` |
| ARMA(1,1) Harvey 形 T/R/Z 手推 (T=[[φ,1],[0,0]], R=[1,θ]', Z=[1,0]) | **精确相等** | `test_statespace.cpp` |
| AR(1) 平稳 P₀ = σ²/(1−φ²) 手推 | 1e-12 | `test_statespace.cpp` |
| KF 前两步 v/F/loglik 手推 (AR(1)) | 1e-12 | `test_statespace.cpp` |
| KF loglik vs statsmodels SARIMAX llf (burn 对齐后) | 1e-8 | `verify_statespace.py` |
| ARIMA statespace MLE 参数 vs statsmodels SARIMAX | 1e-6~1e-8 (优化器容差层, 7C 先例) | `verify_statespace.py` |
| innovations 常驻互检 \|KF loglik − innovations loglik\| | ≤1e-5 告警线 (T=300 实测 2.6e-8) | `test_arima_statespace.cpp` |
| Lyapunov 守卫: φ=1.01 非正定解 → throw | 行为断言 | `test_statespace.cpp` |

### 1.4 接口总则 (沿用 PHASE7C_SPEC §1.4 v1.2, 全模块强制)

1. **形态**: header-only `inline` 自由函数 + `XxxResult` 聚合体返回值 (7B/7C 惯例); 无类继承体系 — ADR-020 K1 "SSModel 基类" 落地为**聚合体载体** `StateSpaceModel` (系统矩阵 + 初始化 spec), 非继承多态
2. **Result 成员一律默认初始化** (`= 0.0` / `= false` / 空容器), 禁止裸 C 数组 (定长用 `std::array<Real, N>{}`)
3. **类型**: `Real`/`Size` 统一取自 `core/types.hpp` (Real = double, §6.6)
4. **参数化**: 封闭集合用 `enum class` (`PssCase`/`BoundsConclusion`/`SSInitMethod`); **跨库对照参数沿用字符串** (`trend` "n"/"c"/"ct"/"ctt", `ic` "aic"/"bic"/"hqic") 保持 statsmodels 直译不转译 (开放问题 #2 裁决)
5. **校验契约**: 输入边界 (最小 T/等长/无 NaN/k≤10/case-trend 相容) `throw std::invalid_argument`; "未计算"字段 = NaN + 伴随 bool 标志
6. **Eigen 边界**: M1/M2 **全部头文件** include `<Eigen/Dense>`, 挂 `cpphub_timeseries_mat` (§6.2); 不引 timeseries→econometrics 新依赖 (Wald 二次型内联实现, 理由见 §6.1)
7. **无状态**: 全部函数纯函数 (无全局可变状态), 线程策略见 §6.7
8. **随机性**: 本期两模块无随机成分 (OLS/KF 确定性; MLE 多起始用确定性初值); `seed` 字段仅 MLE 配置预留 (§6.6)

---

## 2. M1: ARDL/PSS 边限检验 (先行热身)

> **命名空间**: `cpphub::v1::timeseries::ardl` (Eigen3 → cpphub_timeseries_mat)
> **教材锚点**: Pesaran-Shin-Smith 2001 (J Appl Econ 16(3):289-326, 边限检验原始文献) / KS 2020 (OBES 82(6):1456-1481, 响应面) / DK 无关 (纯 OLS 路线)
> **对照库**: statsmodels 0.14.6 `tsa.ardl` (主锚) + R ARDL 0.2.5 (交叉, Natsiopoulos & Tzeremes 2022 JAE 验证)

### 2.1 ARDL 模型 (`ardl_model.hpp`)

**模型形式** (causal=False, A6):

```
y_t = c₀ + c₁·t (+ c₂·t²) + Σ_{i=1..p} φ_i·y_{t−i}
     + Σ_{j=1..k} Σ_{l=0..q_j} θ_{j,l}·x_{j,t−l} + ε_t
```

**接口签名**:

```cpp
namespace cpphub::v1::timeseries::ardl {

struct ArdlSpec {
    Size p = 1;                    ///< y 自身滞后阶
    std::vector<Size> q;           ///< 各 force 变量滞后阶 (含零阶; k = q.size(), 每变量独立 q_j)
    std::string trend = "c";       ///< "n"/"c"/"ct"/"ctt" (statsmodels 直译)
    bool causal = false;           ///< v1.8 恒 false (A6); 传 true → throw (v1.9 落地)
    Size hold_back = 0;            ///< 0 → 默认 max(p, max_j q_j) (statsmodels 同语义)
};

struct ArdlResult {
    std::vector<Real> params;      ///< OLS 系数 (列序 = 算法步骤 Step 2 冻结序)
    std::vector<std::vector<Real>> vcov;  ///< (X'X)⁻¹σ̂², 同列序
    Real sigma2 = 0.0;
    Real loglik = 0.0;             ///< 高斯 llf = −(T/2)(log 2πσ̂² + 1), T = n_obs_used
    Real aic = 0.0, bic = 0.0, hqic = 0.0;  ///< 镜像 statsmodels OLS 口径 (含 σ² 计数; verify_ardl.py 逐值冻结)
    std::vector<Real> fitted, residuals;
    Size n_obs_used = 0;
    Size hold_back_used = 0;
    ArdlSpec spec{};               ///< 回显 (含默认化后的 hold_back)
};

ArdlResult ardl_fit(const std::vector<Real>& y,
                    const std::vector<std::vector<Real>>& x,
                    const ArdlSpec& spec = ArdlSpec{});
///< x 为 k 个等长 force 变量; k=0 (纯 AR) 允许拟合, 但 bounds 检验拒绝
///< (probe 实证: statsmodels 响应面表无 k=0 键, k 域 = 1..10)
///< 校验: y/x 等长, 无 NaN, T > hold_back_used + 参数数, k ≤ 10 (响应面覆盖域)
}
```

**算法步骤**:

```
Step 1: hold_back_eff = (hold_back == 0) ? max(p, max_j q_j) : hold_back
Step 2: 设计矩阵列序 (冻结, params/vcov 共用):
  [trend 列 ("n"=0 列, "c"=[1], "ct"=[1,t], "ctt"=[1,t,t²])]
  + [y_{t−1..t−p}] + [x_{1,t}, x_{1,t−1}, ..., x_{1,t−q_1}]
  + ... + [x_{k,t}, ..., x_{k,t−q_k}]
  样本 t = hold_back_eff .. T−1 (causal=false 含同期 x_{j,t})
Step 3: OLS: β̂ = (X'X)⁻¹X'y (Eigen LLT); σ̂² = SSR/(T_eff − n_params);
  vcov = (X'X)⁻¹σ̂² (nonrobust, A7)
Step 4: llf/IC 镜像 statsmodels OLS 口径 (verify_ardl.py 对照冻结; IC 公式差异属阻断性断言, 冲突回溯 ADR-020)
```

**幻觉点映射**: PSS6 (k = force 变量数, 与滞后阶 q 无关) / PSS8 (样本损失 = hold_back_eff = max(p, max q_j)) / A6 (同期项含)

**测试要点**: 系数 vs statsmodels 1e-10 (trend 四式各一); hold_back 默认与显式; k=0 退化; causal=true 拒绝; NaN/不等长 throw; IC 逐值 vs statsmodels。

### 2.2 阶数选择 (`ardl_select.hpp`)

**接口签名**:

```cpp
struct ArdlSelectResult {
    std::string ic_used = "bic";
    Size p_opt = 0;
    std::vector<Size> q_opt;               ///< 最优各变量阶 (联合网格 argmin ic_used)
    Real ic_min = 0.0;
    std::vector<Real> aic_grid, bic_grid, hqic_grid;  ///< 扁平联合网格 (维度回显)
    Size n_p_grid = 0;                     ///< p ∈ 0..maxlag (n_p_grid = maxlag+1)
    std::vector<Size> q_grid_dims;         ///< 各变量维度 (maxorder+1, uniform)
    Size n_obs_used = 0;                   ///< 全网格统一 (同一样本)
    ArdlResult best_model{};
};

ArdlSelectResult ardl_select_order(const std::vector<Real>& y,
                                   const std::vector<std::vector<Real>>& x,
                                   Size maxlag, Size maxorder,
                                   const std::string& ic = "bic",
                                   const std::string& trend = "c");
///< 网格 = (maxlag+1) × (maxorder+1)^k 联合候选, glob (lag 子集搜索) 不做 (A3)
///< hold_back 全网格统一 = max(maxlag, maxorder) — 所有候选同一估计样本 (镜像 var_select offset 机制, A3)
///< maxorder 仅 uniform int (per-variable dict 不做, 手动 fit 替代)
```

**算法步骤**:

```
Step 1: hold_back_grid = max(maxlag, maxorder); 同一样本遍历全部候选
Step 2: 每候选 ardl_fit → 三 IC 落网格; argmin(ic_used) 选优
Step 3: 回显 best_model (最优模型完整 ArdlResult)
```

**幻觉点映射**: A3 (同一样本) / PSS8 (网格 hold_back 语义)

**测试要点**: 最优阶 + IC 值 vs statsmodels ardl_select_order 1e-10; 同一样本断言 (各候选 n_obs_used 相等); ic 参数三分支; k=1/k=2 网格维度。

### 2.3 UECM 重参数化与边限检验 (`uecm_bounds.hpp`)

**UECM 形式** (PSS 2001, 理论核心 [A]):

```
Δy_t = c₀ + c₁·t + π_y·y_{t−1} + Σ_{j=1..k} π_j·x_{j,t−1}
     + Σ_{i=1..p−1} ψ_{y,i}·Δy_{t−i} + Σ_{j=1..k} Σ_{l=1..q_j−1} ψ_{j,l}·Δx_{j,t−l}
     + Σ_{j=1..k} ω_j·Δx_{j,t} + ε_t
```

**接口签名**:

```cpp
enum class PssCase { I = 1, II = 2, III = 3, IV = 4, V = 5 };  // 编号跨软件一致 (PSS/statsmodels/R ARDL/pssbounds)
enum class BoundsConclusion {
    NoLevelRelation,   ///< F: stat < I(0) 下界 (t: stat > I(0) 界) — 不能拒绝 H0
    Inconclusive,      ///< 介于两界 (边限检验固有, 必须显式返回)
    LevelRelation      ///< F: stat > I(1) 上界 (t: stat < I(1) 界) — 拒绝 H0
};

struct UecmResult {
    std::vector<Real> params;          ///< UECM 系数 (列序: 截距/trend, π_y, π_1..π_k, ψ_y, ψ_j, ω)
    std::vector<std::vector<Real>> vcov;
    std::vector<Real> pi;              ///< (π_y, π_1..π_k) 水平项系数 — 检验对象
    std::vector<std::vector<Real>> pi_vcov;  ///< π 子块 vcov
    Size n_obs_used = 0;
    Real sigma2 = 0.0;
    ArdlSpec source_spec{};
};

UecmResult uecm_from_ardl(const ArdlResult& ardl);
///< ARDL → UECM 系数线性重参数化 (B 变换, 镜像 statsmodels UECM.from_ardl), 不重估

struct BoundsTestResult {
    Real stat = 0.0;                   ///< F (bounds_f_test) 或 t (bounds_t_test)
    bool is_t_test = false;
    PssCase case_used{};
    Size k = 0;                        ///< force 变量数 (PSS6)
    Size n_restrictions = 0;           ///< F 约束数: Case I/III/V = k+1, II/IV = k+2
    std::array<std::array<Real, 4>, 2> crit_vals{};  ///< [0]=I(0) 行 [1]=I(1) 行 × {10,5,1,0.1}% (statsmodels 四档, PSS9)
    std::array<Real, 2> p_values{};    ///< F: {I(0) null, I(1) null} 响应面近似 (probe [B] 冻结);
                                       ///< t: {NaN, NaN} + has_p=false (无 t 响应面, probe [D]; R ARDL 的 t p 值系 MC 自算, A5 不做)
    BoundsConclusion conclusion{};
    bool has_p = true;                 ///< false = t 检验 (无响应面 p 值)
    std::string message;
};

BoundsTestResult bounds_f_test(const UecmResult& uecm, PssCase c, Real alpha = 0.05);
BoundsTestResult bounds_t_test(const UecmResult& uecm, PssCase c, Real alpha = 0.05);
///< t 检验仅 Case I/III/V (PSS Table CII 边界, A4); Case II/IV → throw + message (镜像 pssbounds)
///< alpha 限四档 {0.10, 0.05, 0.01, 0.001} 之一, 否则 throw (无插值语义)
```

**算法步骤**:

```
Step 1: uecm_from_ardl 重参数化:
  π_y = −(1 − Σφ_i); π_j = θ_{j,0}·(level relation 代数) — B 变换镜像 statsmodels
  UECM 损失样本 = ARDL 同一样本 (重参数化不重估, 隐含 max(p, q_j) 已扣)
Step 2: F 统计量 (PSS7: F 版 = Wald 二次型 / 约束数):
  Case I:   H0: π_y = π_1..π_k = 0,               约束数 k+1
  Case II:  H0: 上 + 受限截距进水平项,             约束数 k+2
  Case III: H0: π_y = π_1..π_k = 0 (截距不限),     约束数 k+1
  Case IV:  H0: 上 + 受限趋势进水平项,             约束数 k+2
  Case V:   H0: π_y = π_1..π_k = 0 (截距/趋势均不限), 约束数 k+1
  F = (Rβ̂)' (R V R')⁻¹ (Rβ̂) / q   [Wald 二次型内联实现, Eigen LLT]
Step 3: t 统计量 = π_y / SE(π_y) (vcov 对角元开方)
Step 4: case-trend 相容校验 (镜像 R ARDL 报错):
  Case I ↔ "n"; II/III ↔ "c"; IV/V ↔ "ct"; "ctt" 不支持 bounds (throw);
  请求 case 的确定性项必须存在于已估模型, 否则 throw invalid_argument
Step 5: CV = pss_crit_values(k, case, i1); p = pss_pvalue_bounds(...); 结论:
  ⚠️ k 查表口径 = k_PSS (force 变量数, probe 裁决 PSS10): statsmodels bounds_test
  内部用 nvar = 1+k_PSS 偏移查表, 本库不镜像 — stat 对照不受影响 (Wald 不用 k 键),
  CV/p 对照走 R ARDL + statsmodels 表键直查
  F: stat > CV_I1(alpha) → LevelRelation; stat < CV_I0(alpha) → NoLevelRelation; 否则 Inconclusive
  t: 方向相反 — stat < CV_I1 (更负界, 如 Case III k=1 5%: −3.22) → LevelRelation;
     stat > CV_I0 (−2.86) → NoLevelRelation; 否则 Inconclusive
```

**幻觉点映射**: PSS3 (t 表仅 Case I/III/V) / PSS6 (k 定义) / PSS7 (F 版非 χ²) / PSS8 (样本损失) / PSS9 (statsmodels 四档 {10,5,1,0.1}%, PSS 原表 {10,5,2.5,1}% — 无 2.5% 档) / **PSS10 (statsmodels bounds_test 查表 k 偏移 +1, probe 实证, 本库用 PSS 语义)** / A4 (t 对照源 = R ARDL bounds_t_test, statsmodels 无 t 检验 — probe [D] 实证) / A7 (nonrobust)

**测试要点**: 重参数化系数 vs UECM.from_ardl 1e-10; F 统计量 vs bounds_test 1e-10 (五 case, stat 不受 k 偏移影响); t 三 case vs **R ARDL bounds_t_test** 1e-8 (statsmodels 无 t, probe [D]); Case II/IV t → throw; t 方向 (负侧拒绝, 构造已知 LevelRelation 数据); F 约束数五 case 断言; R ARDL bounds_f_test 交叉 1e-8 (CV/p 同口径); 三值结论构造数据覆盖三态; PSS 原表锚 4.94/5.73/−2.86/−3.22 交叉印证 (~0.05); **danish 案例 k 偏移回归断言** (k_PSS=3 时本库 CV = statsmodels 键 (3,3) ≠ 其 bounds_test 返回值, 防止未来误镜像)。

### 2.4 PSS 响应面临界值 (`pss_critical_values.hpp`)

**接口签名**:

```cpp
namespace cpphub::v1::timeseries::ardl::detail {
// pss_response_surface.inc: statsmodels pss_critical_values.py 全量系数转录 (probe [A] 冻结)
// 四张 constexpr 表: large_p (100 键 × 4 系数) / small_p (100 键 × 3) /
//   stat_star (100 标量) / crit_vals (100 键 × 4 档) — 共 1200 数字
// 键 (k ∈ 1..10, case ∈ 1..5, I1 ∈ {false,true}); 无 k=0 键 (probe [A] 实证)
// static_assert 双锚 (编译期):
//   crit_vals[(1, 3, false)][1] == 4.9228967
//   crit_vals[(1, 3, true)][1]  == 5.7217442
}  // namespace detail

/// 4 档临界值 {10%, 5%, 1%, 0.1%} (statsmodels 口径; PSS 原表 {10,5,2.5,1}% — PSS9 文档声明)
/// k 即 force 变量数 (表键语义直译, PSS10 裁决); k ∉ [1,10] → throw
std::array<Real, 4> pss_crit_values(Size k, Size pss_case, bool upper_i1);

/// p 值界 (双 null 近似, 镜像 statsmodels _pss_pvalue; 算法 probe [B] 冻结):
///   p = 1 − Φ(Σ_{i=0..n} cᵢ·log(stat)ⁱ), 其中 n = 3 (large) 或 2 (small);
///   系数选择: stat > stat_star ? small_p : large_p (阈值比较用原始 stat)
///   仅 F 分布覆盖 (probe [D]: statsmodels 无 t p 值路径; t 的 p 界 v1.9 评估)
std::array<Real, 2> pss_pvalue_bounds(Real stat, Size k, Size pss_case);
```

**算法步骤**:

```
Step 1 (✅ 前置已完成, probe_pss_cv.py 2026-08-21, 报告冻结):
  (a) 布局: 4 张表 × 100 键 (k=1..10 × case 1..5 × I1); large 4 系数 / small 3 / star 1 / cv 4
  (b) 路径: crit_vals = 表常量直查 (非多项式求值); p = 1−Φ(poly(log stat)), stat_star 切换
  (c) 全量 diff: probe 脚本已建框架, .inc 落地后逐键逐值 (1200 数字, 风险 #1)
Step 2: 转录落 .inc (constexpr) + static_assert 双锚 + verify 全量 diff
Step 3: 求值接口 (Φ 用 core/special_functions 正态 CDF; 与 statsmodels 数值自证
  diff ≤ 1e-16 已在 probe [B] 验证); k ∉ [1,10] → throw
```

**幻觉点映射**: PSS1 (ardlverse CV 错误实锤 — 禁作对照源, verify 脚本黑名单) / PSS2 (KS 论文出处 OBES 2020) / PSS4 (Narayan 覆盖仅 F/Case II-V/T=30-80 — A5 不采) / PSS5 (statsmodels 分位档 (90,95,99,99.9) 与 PSS (90,95,97.5,99) 错位, 第 3 列 1% 非 2.5%) / PSS10 (statsmodels bounds_test 查表 k 偏移 +1 — 表键是 PSS 语义, 调用层异类)

**测试要点**: 双锚 static_assert; 全量 diff (probe 脚本); k=1/k=10 边界; k=0/k=11 throw; 四档数值抽查 vs statsmodels crit_vals 直查; p 值 vs _pss_pvalue 抽查 (stat 跨 stat_star 两侧各一)。

---

## 3. M2: Kalman/statespace (Eigen 载体)

> **命名空间**: `cpphub::v1::timeseries::statespace` (Eigen3 → cpphub_timeseries_mat, K5)
> **教材锚点**: Durbin-Koopman 2012 *Time Series Analysis by State Space Methods* 2nd ed Ch.4 (常规 KF) / Ch.5.6-5.7 (diffuse, 本期仅 approximate) / Harvey 1989 (ARIMA 状态空间) / Francke-Koopman-De Vos 2010 JTSA (diffuse 似然不可用于估计, SS9)
> **对照库**: statsmodels 0.14.6 `tsa.statespace` (SARIMAX 主锚; ARIMA.fit method='statespace' 官方双路径之一)

### 3.1 状态空间模型与 Kalman 滤波 (`ss_model.hpp`)

**状态空间形式** (DK 记号):

```
y_t = Z_t α_t + d_t + ε_t          (ε ~ N(0, H)),  v1.8 p = 1 (单变量)
α_{t+1} = T_t α_t + c_t + R_t η_t  (η ~ N(0, Q)),  r ≤ m (k_posdef 约束)
```

**接口签名**:

```cpp
namespace cpphub::v1::timeseries::statespace {

enum class SSInitMethod { Known, Stationary, ApproximateDiffuse };

struct SSInitSpec {
    SSInitMethod method = SSInitMethod::ApproximateDiffuse;
    Eigen::VectorXd a0;          ///< Known: E[α₁] (空 = 零向量)
    Eigen::MatrixXd P0;          ///< Known: Var[α₁]
    Real kappa = 1e6;            ///< ApproximateDiffuse (statsmodels SARIMAX 非平稳默认口径, K1)
};

struct StateSpaceModel {
    Eigen::MatrixXd Z;           ///< (1×m) 观测矩阵
    Eigen::MatrixXd T;           ///< (m×m) 状态转移
    Eigen::MatrixXd R;           ///< (m×r) 状态噪声载入 (r ≤ m, SS3: R 非方阵, Q 为 r×r)
    Eigen::MatrixXd Q;           ///< (r×r) 状态噪声方差
    Eigen::MatrixXd H;           ///< (1×1) 观测噪声方差 (ARMA: 零阵, SS6)
    Eigen::VectorXd d;           ///< (1) 观测截距 (空 = 0)
    Eigen::VectorXd c;           ///< (m) 状态截距 (空 = 0)
    SSInitSpec init{};
};

struct KalmanResult {
    Real loglik = 0.0;                 ///< burn 后 Σ −½(log|F_t| + v_t'F_t⁻¹v_t + log 2π) (SS8: 单变量常数 = log 2π)
    std::vector<Real> loglik_contribs; ///< 逐 t (burn 后)
    std::vector<Real> innovations;     ///< v_t (burn 前后均记录)
    std::vector<Real> pred_var;        ///< F_t
    std::vector<Eigen::VectorXd> state_pred;      ///< a_{t+1|t}
    std::vector<Eigen::MatrixXd> state_cov_pred;  ///< P_{t+1|t}
    Size n_burn = 0;
    Size n_obs = 0;
    bool ok = false;
    std::string message;
};

/// 平稳初始化协方差: 解离散 Lyapunov 方程 T P₀ T' − P₀ = −R Q R'
/// a₀ = (I−T)⁻¹ c (c 空则 0)
/// 守卫 (SS4): 解后对称化 + LLT 正定校验; 非正定 (单位根邻域 φ≈1.01 静默负方差) → throw
Eigen::MatrixXd stationary_cov(const StateSpaceModel& m);

/// 常规 Kalman 滤波 (DK Ch4 更新形)
KalmanResult kalman_filter(const std::vector<Real>& y,
                           const StateSpaceModel& model,
                           Size loglikelihood_burn = 0);
///< init.method=Stationary 且 T 特征值 ≥1 → 内部转 ApproximateDiffuse 并 message 记录 (非静默)
}
```

**算法步骤**:

```
Step 1: 初始化 (三法):
  Known:            a₁ = a0, P₁ = P0
  Stationary:       a₁ = (I−T)⁻¹c, P₁ = Lyapunov 解 (迭代/双 Hessenberg;
                    守卫: 对称化 + LLT 正定校验, 失败 throw, SS4)
  ApproximateDiffuse: a₁ = 0, P₁ = κI (κ=1e6, statsmodels SARIMAX 非平稳默认)
Step 2: 逐 t 递推 (DK Ch4, 更新形 — 与 Joseph 形代数等价但不可混搭):
  v_t = y_t − Z a_t − d
  F_t = Z P_t Z' + H                    (单变量: 标量; F_t ≤ 0 → throw runtime_error)
  K_t = T P_t Z' F_t⁻¹
  a_{t+1} = T a_t + K_t v_t + c
  P_{t+1} = T P_t (T − K_t Z)' + R Q R'   ⚠️ DK 更新形 (非 Joseph 形)
Step 3: loglik 贡献 (t ≥ burn): −½(log F_t + v_t²/F_t + log 2π)  (SS8)
Step 4: 输出 loglik = Σ 贡献; innovations/pred_var 全程记录 (诊断用)
```

**幻觉点映射**: SS2 (T 方向与更新形) / SS3 (R m×r 非方阵, Q r×r, RQR' 秩 r 注入) / SS4 (Lyapunov 守卫) / SS6 (H=0 时 F 正定性全靠 P) / SS8 (似然常数单变量 log 2π)

**测试要点**: AR(1) 手推黄金锚 (a₁=0, P₁=σ²/(1−φ²), v₁=y₁, F₁=P₁, loglik₁=−½(log P₁ + y₁²/P₁ + log 2π)) 1e-12; Lyapunov AR(1) 精确解 1e-12; φ=1.01 守卫触发 throw; κ=1e6 diffuse 首期 F=κ 量级断言; loglik 常数项 (与 statsmodels llf 对齐时 burn 口径敏感, §3.2); F_t≤0 throw。

### 3.2 ARIMA 状态空间精确似然 MLE (`arima_statespace.hpp`)

**系统矩阵** (Harvey 形, SS2):

```
ARMA(p,q) 块 (m₁ = max(p, q+1)):
  T 第一列 = AR 系数 φ (上 Hessenberg 伴随), 次对角单位阵
  R = (1, θ₁, …, θ_{m₁−1})' (m₁×1), Q = [σ²], Z = (1, 0, …, 0), H = 0
  黄金锚 ARMA(1,1): T = [[φ, 1], [0, 0]], R = [1, θ]', Z = [1, 0]
ARIMA(p,d,q) 嵌入差分 (simple_differencing=False, K1):
  状态维 m = max(p, q+1) + d; d 个积分分量; 全部数据可用
  ⚠️ T 布局由 probe_sarimax.py 从 statsmodels sarimax.py 源码逐元素探明后冻结
    (黄金锚 ARIMA(1,1,1) 手推三矩阵落测试)
```

**接口签名**:

```cpp
struct ArimaSSConfig {
    bool simple_differencing = false;    ///< v1.8 默认 false = 嵌入差分 (K1 四件套之一); true 口径对照用例
    bool include_constant = false;       ///< d=0: 观测截距; d≥1: 差分方程常数 = levels drift (经 d_t 注入)
    bool enforce_stationarity = true;    ///< AR 平稳化变换 (Monahan 1984, statsmodels _ar_trans_params 同构)
    bool enforce_invertibility = true;   ///< MA 可逆化变换 (_ma_trans_params 同构)
    Size n_starting_points = 2;          ///< {innovations 解, CSS 解} 确定性初值 (7C 多起始惯例)
    Size seed = 42;                      ///< 预留 (本期无随机扰动)
    SLSQP::Config optimizer_config{};    ///< K3: 复用 ADR-018 (arima CSS-ML 同款)
    bool cross_check_innovations = true; ///< K1: 常驻互检
};

struct ArimaSSResult {
    std::vector<Real> phi, theta;        ///< 正号参数化 (1+θB), 与 arima/innovations 同源 (AR1)
    Real sigma2 = 0.0;                   ///< 非集中化全参数优化 (K4)
    Real loglik = 0.0;                   ///< KF loglik (burn 口径回显)
    Real aic = 0.0;
    Size loglikelihood_burn = 0;
    bool converged = false;
    Size n_iterations = 0;
    std::string message;
    Real innovations_loglik = 0.0;       ///< 互检: innovations (T−d 差分样本) loglik
    Real innovations_gap = 0.0;          ///< |KF − innovations|; > 1e-5 → message 告警 (常驻协议)
    std::vector<Real> innovations;       ///< v_t (诊断)
};

/// ARMA(p,q) Harvey 形系统矩阵 (d=0)
StateSpaceModel make_arma_statespace(const std::vector<Real>& phi,
                                     const std::vector<Real>& theta,
                                     Real sigma2);

/// ARIMA(p,d,q) 嵌入差分系统矩阵 (simple_differencing=true 时为先差分后 ARMA 形)
StateSpaceModel make_arima_statespace(const std::vector<Real>& phi,
                                      const std::vector<Real>& theta, Size d,
                                      Real sigma2, bool include_constant,
                                      bool simple_differencing = false);

/// ARIMA 精确似然 MLE (Kalman 路径)
ArimaSSResult arima_statespace_mle(const std::vector<Real>& data,
                                   const ArimaSpec& spec,          ///< 复用 arima/arima_model.hpp
                                   const ArimaSSConfig& config = ArimaSSConfig{});
///< 初始化: 平稳 (|特征根|<1) → Stationary; 否则 → ApproximateDiffuse κ=1e6 (K1)
///< burn 默认 = d (开放问题 #3 裁决); statsmodels 默认公式 probe_sarimax.py 探明, 差异文档声明
///< σ² 非集中化 (K4): 与 innovations_mle 集中化口径差异声明, 交叉验证以 loglik 为主锚
}
```

**算法步骤**:

```
Step 1 (前置, M2 第一任务): probe_sarimax.py 探明:
  (a) 嵌入差分 T/R/Z 构造 (sarimax.py _statespace_matrix 逐元素)
  (b) loglikelihood_burn 默认公式 (statsmodels None 时的解析)
  (c) transform_params 变换细节 (_ar/_ma_trans_params)
  (d) llf 口径 (含常数项与 burn 应用方式)
Step 2: 参数变换: 无约束优化变量 ⇄ (φ, θ, σ²) 经 Monahan AR 变换 + MA 可逆变换 (enforce 双 True 默认)
Step 3: 似然评估: make_arima_statespace → kalman_filter → loglik (burn=d)
Step 4: MLE: SLSQP/BFGS (K3), 初值 = {innovations_mle 解, arima CSS 解} (确定性多起始)
Step 5: 常驻互检 (K1): innovations_mle (T−d 差分样本 loglik) vs KF loglik (burn=d);
  gap > 1e-5 → message 告警 (非 throw); 7C 一次性仲裁 (T=300: 2.6e-8) 升级为常驻回归测试
Step 6: AIC = −2·loglik + 2·(p+q+d_drift+1) (口径与 statsmodels 对齐, verify 冻结)
```

**幻觉点映射**: SS1 (approximate diffuse κ=1e6 ≠ exact diffuse κ→∞ 解析极限 — 后者 v1.9) / SS2 (Harvey 形 T 第一列 AR + Z=(1,0,…,0); Hamilton 形 AR 入第一行 + Z=(0,…,0,1), 混用 = 模型级错误) / SS5 (两差分口径似然不同, 对照必须统一) / SS7 (Chandrasekhar 三限制, v1.9 不做) / SS9 (diffuse 似然不可普遍用于参数估计, FKV 2010 — 本期 approximate diffuse 仅作初始化) / SS10 (UCM 一元专属 — 多元边界) / K4 (σ² 口径差异声明)

**测试要点**: ARMA(1,1) 三矩阵手推精确断言; ARIMA(1,1,1) 嵌入布局 (probe 冻结后); KF loglik vs SARIMAX llf 1e-8 (burn 对齐); MLE 参数 vs SARIMAX 1e-6~1e-8; innovations 互检 gap ≤1e-5 (T=300 回归锚 2.6e-8); simple_differencing 双口径 loglik 差异显式 (不混比); AR 平稳化变换域测试 (优化变量 ↔ φ 映射往返); drift 用例 (d=1 + constant ↔ statsmodels trend='c'); 非平稳 → diffuse 自动切换 + message。

---

## 4. M3: 战略选项 (D-2 待裁决, 占位)

> [PHASE8_SCOPE_RESEARCH.md](../../research/PHASE8_SCOPE_RESEARCH.md) §6 D-2: (a) BVAR/TVP-VAR R 桥 [DSGE 衔接] / (b) KDE/KNN + TE [信息论铺路] — 取决于 Research OS 主线紧迫性

- **状态**: 待 Scott 裁决; 排期最晚, 不阻塞 M1/M2 实施与验收
- **占位约束**: M1/M2 收尾验收 (§8) 不含 M3; D-2 裁决后另立 spec 增补 (v1.1) 或入 v1.9
- **依赖预扫**: (a) 路线需 R 桥基础设施 (BVAR 1.0.5 / bvarsv, 7C 开放问题 e 已预扫); (b) 路线需 KDE/KNN 估计器 + 传递熵 (SEMIPARAMETRIC_LANDSCAPE.md 已预扫)

---

## 5. 端到端集成测试 (`test_integration_phase8.cpp`)

| # | 场景 | 断言要点 |
|---|------|---------|
| 1 | ARDL 完整链: select → fit → UECM → bounds F+t | 模拟水平关系数据 (β=1 协整生成) → LevelRelation 结论; k=2 双 force 变量 |
| 2 | ARIMA 精确似然三路径: CSS-ML / innovations / statespace | 同数据 loglik 落点一致 (各自口径声明下); 参数一致 1e-6 |
| 3 | 单位根预检 → 边限工作流: ADF 预检 (复用 7B) → 混合 I(0)/I(1) force 变量 → bounds 检验 | x 含 I(2) → 前置拒绝路径 (文档契约: 调用方预检责任) |
| 4 | Kalman 诊断: innovations 序列 → LB/JB (复用 7A residual_diagnostics) | statespace v_t 与 OLS 残差路径打通 |
| 5 | 双差分口径: ARIMA(0,1,1) 嵌入 vs 先差分 | 两 loglik 数值差异显式记录 (SS5), 不混比断言 |
| 6 | Lyapunov 守卫端到端: φ=1.01 模拟 → stationary_cov throw → diffuse 回退路径 | 守卫非静默 + 回退正确 |

---

## 6. 依赖关系与复用

### 6.1 模块依赖图

```
M1 ARDL/PSS
├── ardl_model ← Eigen LLT (var_model 同模式自建 OLS, 不依赖 econometrics OLSEstimator)
├── ardl_select ← ardl_fit (同一 hold_back 网格)
├── uecm_bounds ← ardl_fit 结果重参数化 + Wald 二次型内联 (Eigen, ~10 行)
│     ⚠️ 不复用 econometrics/inference/hypothesis_tests.wald_test — 其 χ²/F p 值对
│     非标准 PSS 分布无意义 (A7), 仅二次型有用; 内联实现避免 timeseries→econometrics
│     新 target 依赖 (granger_test→special_functions 为纯头文件先例, wald_test 需
│     cpphub_econometrics 链接)
└── pss_critical_values ← .inc 转录 + probe_pss_cv.py (前置)
M2 Kalman/statespace
├── ss_model ← Eigen (纯新; DK Ch4)
├── arima_statespace ← SLSQP (ADR-018, K3) + arima::ArimaSpec (头文件复用, 纯标量无 Eigen)
│     + innovations_mle (互检对手, K1) + probe_sarimax.py (前置)
└── 诊断出口 → econometrics/inference/residual_diagnostics (7A, 集成场景 4)
```

### 6.2 Eigen3 隔离 (ADR-013/ADR-017)

| 模块 | Eigen3 | 载体 |
|------|--------|------|
| M1 (ardl 全部 4 头文件) | **是** | `cpphub_timeseries_mat` (现有 INTERFACE, K5 同源) |
| M2 (statespace 全部 2 头文件) | **是** | `cpphub_timeseries_mat` |

**零新增 CMake target** (K5: 镜像 var/cointegration 隔离模式, CMakeLists L51-58 现成); ardl/statespace 头文件不被单变量纯标量头 (arima/garch/unit_root/midas) 反向依赖 (grep 断言, §8.5)。

### 6.3 C ABI (兼容性约束 4)

新模块头文件**不进 C ABI** (v1.7 现状延续); Python 绑定不扩展 (bindings.cpp 不动)。

### 6.4 构建与版本矩阵 (沿用 §8.4 冻结值)

CMake 3.25 / C++20 (`cxx_std_20`) / MSVC 19.5x + GCC 13.3 / Eigen 3.4.0 vendored (`EIGEN_MPL2_ONLY`) / googletest 1.14.0 — 全部不变; 测试注册追加 (tests/CMakeLists.txt, 唯一允许触碰的既有文件):

```cmake
# ---- Phase 8 M1: ARDL/PSS (spec PHASE8 v1.8 §1.2, Eigen3 → cpphub_timeseries_mat) ----
add_executable(test_ardl_model unit/timeseries/test_ardl_model.cpp)
target_link_libraries(test_ardl_model PRIVATE cpphub_timeseries_mat GTest::gtest_main)
cpphub_discover_tests(test_ardl_model)
# (test_ardl_select / test_uecm_bounds / test_statespace / test_arima_statespace /
#  test_integration_phase8 同模式)
```

### 6.5 兼容性: v1.7.1 零破坏条款

- **additive-only**: 仅新增头文件 (6) / tests 目录新增 / `tests/CMakeLists.txt` 追加注册 — 不修改任何 v1.7.1 既有头文件与源文件
- 既有 **2447 测试零回归**; SLSQP/arima/innovations 既有套件无退化
- 实施中若确需触碰既有文件 (含 innovations_mle.hpp 注释性互检升级) ⇒ 视同 scope 变更, 走 ADR-020 修订记录, 禁止静默

### 6.6 数值精度与可复现性政策

- `Real = double`; 编译 flag 政策不变 (MSVC `/fp:precise` + GCC `-ffp-contract=off`, 禁 `-ffast-math`) — FMA 收缩关闭 ⇒ 三平台位一致, `.inc` 黄金值可精确断言
- 本期两模块无随机成分 (OLS/KF 确定性; MLE 多起始用 innovations/CSS 确定性解, 不引入随机重启)
- 临界值系数 constexpr + static_assert 双保险 (与 mackinnon_tables.inc 管线同惯例)

### 6.7 线程安全与并行策略

全部公共接口单线程纯函数 (§1.4-7); 不引入 OpenMP/std::thread; IC 网格搜索的并行化留给调用方 (外层循环天然独立)。

### 6.8 性能预算 (单线程)

| # | 操作 | 预算 | 依据 |
|---|------|------|------|
| P1 | ardl_select maxlag=4, maxorder=3, k=2, T=250 (全网格) | < 5 s | (5×16) 候选 × 小维 LLT OLS |
| P2 | bounds_f_test T=250 | < 0.1 s | Wald 二次型 + 响应面求值 |
| P3 | KF loglik ARIMA(2,1,2) T=1000 单次 | < 50 ms | O(T·m³), m = max(2,3)+1+... ≤ 8 → ~10⁷ flop |
| P4 | arima_statespace_mle T=1000 | < 10 s | ~200 似然评估 × P3 + SLSQP 开销 |
| P5 | 全量 ctest (~2527) | < 50 min | v1.7.1 三平台基线 + 增量 |

---

## 7. 幻觉点核查清单 (实施时逐点验证, probe/verify 脚本逐条落盘)

> 编号总计: PSS 10 (含 probe 新增 PSS10) + SS 10 = **20 编号** (调研报告附录 A 19 条 + probe_pss_cv.py 实证新增 1 条)。
> 核查模式沿用 7C: probe_*.py 打印对照库源码 + 小样本手算 + .inc 黄金值断言。

### 7.1 域 A: ARDL/PSS (PSS1-PSS10)

| ID | 影响级 | 核查方式 | 容差 |
|----|--------|---------|------|
| PSS1 (ardlverse CV 错误实锤: Case 3 F 5% k=1 = 2.79/3.67 vs PSS 4.94/5.73, 方向相反) | 极高 | 对照库黑名单 (verify 脚本头部声明); 不做数值核查 (不采信即核查) | — |
| PSS2 (KS 2020 出处 = OBES 82(6), 非 Stata Journal 2019) | 低 | 文档引用核对 | — |
| PSS3 (t 界仅 Case I/III/V, PSS Table CII; 声称 II/IV 有 t 表即杜撰) | 极高 | bounds_t_test Case II/IV → throw; verify_bounds.py 三 case 逐值 | 1e-10 |
| PSS4 (Narayan 仅 F/Case II-V/T=30-80; k≤7 边界 B-) | 高 | A5 不采 Narayan — 无实现面; spec 文档声明边界 | — |
| PSS5 (statsmodels 分位档 (90,95,99,99.9): 第 3 列 1% 非 PSS 2.5%, 第 4 列 0.1% PSS 无) | 极高 | probe_pss_cv.py 逐列核对 + crit_vals 四档语义冻结 | 精确 |
| PSS6 (k = force 变量数, 不含 y 滞后; 与 q 混淆是首要错误源) | 极高 | k=1 vs q=2 构造用例 (同 k 不同 q → 同一 CV 键) | 精确 |
| PSS7 (F 版 = 二次型/约束数, 非 χ² 版) | 极高 | verify_bounds.py F 逐值对照; n_restrictions 五 case 断言 | 1e-10 |
| PSS8 (UECM 样本损失 = max(p, max q_j); 统计量与 CV 须同有效样本) | 高 | n_obs_used 断言 + select 同一样本 | 精确 |
| PSS9 (HC/HAC Wald 分布不在 PSS/KS 表覆盖内) | 高 | A7: 仅 nonrobust; API 无 cov_type 参数 (不存在即核查) | — |
| **PSS10 (statsmodels bounds_test 查表 k 偏移 +1: 表键 = PSS k 但调用用 nvar = 1+k_PSS; probe [C]/[E] danish 实证 + R ARDL 文档 null.value 实录)** | 极高 | 本库 k = force 数 (PSS 语义); danish 回归断言: 本库 CV = 键 (3,3) ≠ statsmodels bounds_test 返回 (键 (4,3)); probe_pss_cv.py [E] 段留档 | 精确 |

### 7.2 域 B: Kalman/statespace (SS1-SS10)

| ID | 影响级 | 核查方式 | 容差 |
|----|--------|---------|------|
| SS1 (approximate diffuse κ=1e6 ≠ exact diffuse κ→∞ 解析极限) | 高 | API 仅 ApproximateDiffuse (枚举无 Exact 项即核查); 文档声明 FKV 2010 | — |
| SS2 (Harvey 形 T 第一列 AR + Z=(1,0,…,0); Hamilton 形入第一行 — 混用 = 模型级错误) | 极高 | ARMA(1,1) 手推三矩阵精确断言 + probe_sarimax.py 源码逐元素 | 精确 |
| SS3 (R 非方阵 m×r, Q r×r 标量阵, RQR' 秩 r 注入) | 高 | 维度校验 throw + ARIMA R 构造断言 | 精确 |
| SS4 (Lyapunov 无守卫 → 单位根邻域 φ≈1.01 静默负方差) | 极高 | φ=1.01 触发 throw 用例 + 对称化/LLT 守卫路径 | 行为 |
| SS5 (两差分口径似然不同, simple_differencing 默认 False 嵌入式) | 极高 | 双口径测试: 两 loglik 差异显式记录, 禁止混比断言 | 1e-8 |
| SS6 (ARMA 观测 H=0, F 正定性全靠 P) | 高 | H 零阵断言 + F_t ≤ 0 throw 路径 | 行为 |
| SS7 (Chandrasekhar 三限制 + 增益 33-50% 非无条件 m³→m²; v0.12.0 新增) | 中 | v1.9 不做 — 无实现面; spec 文档留档 (审计已锚定 v0.12.0) | — |
| SS8 (似然常数: 单变量 log 2π, 多元才是 p·log 2π) | 高 | loglik 手推锚含常数项; llf vs statsmodels 1e-8 (burn 对齐) | 1e-12 |
| SS9 (diffuse likelihood 不能普遍用于参数估计, FKV 2010 JTSA) | 中 | approximate diffuse 仅作初始化 (非似然目标); 文档声明 | — |
| SS10 (UCM 一元专属; 多元 = VARMAX/DynamicFactor 边界) | 低 | scope 不做清单声明 | — |

---

## 8. 验收标准 (G 门禁)

### 8.1 编译与测试
- [ ] MSVC Release: 全量通过 (新增 ~80 + 现有 2447 ≈ 2527)
- [ ] A 站 GCC: fresh clone + rebuild + ctest (submodules: recursive)
- [ ] B 站 GCC: fresh clone + rebuild + ctest
- [ ] 三平台无数值偏差 (容差分层表 §1.3)

### 8.2 数值基准 (§1.3 全表)
- [ ] M1: statsmodels ARDL/select/UECM/bounds 1e-10 全绿; R ARDL 交叉 1e-8; 响应面转录精确相等 (全量 diff); 双锚 static_assert
- [ ] M2: ARMA(1,1) 手推精确; KF vs SARIMAX llf 1e-8; MLE 参数 1e-6~1e-8; innovations 互检 gap ≤1e-5

### 8.3 幻觉点覆盖
- [ ] §7 两域 19 编号全部核查, probe/verify 脚本逐条落盘
- [ ] 前置工程: probe_pss_cv.py (M1 第一任务) + probe_sarimax.py (M2 第一任务) 报告冻结后才冻结对应实现
- [ ] pss_response_surface.inc 全量 static_assert 可行性 (编译期内存允许时逐键, 否则双锚 + 运行时全量)

### 8.4 文档与流程对齐
- [ ] ADR-020 (12 项) 全部实施且无静默偏离 (冲突回溯 ADR 修订)
- [ ] ADR-017 命名空间 + §6.2 Eigen 隔离遵守; 兼容性约束 C1-C5 逐项满足
- [ ] DEVELOPMENT_LOG.md 逐里程碑更新; PHASE8_FINAL_ACCEPTANCE.md + checklist (实施前编制)

### 8.5 Scope 边界 (不做清单)
- [ ] 未实现: panel/quantile/Fourier/NARDL / glob 搜索 / Narayan 分档 / MC 模拟路径 / HC/HAC Wald / causal ARDL / exact diffuse / 平滑器 / 缺失数据 / 多元 statespace / Chandrasekhar / SARIMA 季节
- [ ] 单变量纯标量头 (arima/garch/unit_root/midas) 零 ardl/statespace 反向 include (grep 断言)
- [ ] M3 未实现且不阻塞验收 (D-2 待裁决声明)

---

## 9. 风险与缓解 (调研风险表全量继承 + 实施期新增)

| # | 风险 | 等级 | 缓解 |
|---|------|------|------|
| 1 | 响应面系数转录错误 (500+ 数字) | 高 | probe_pss_cv.py 全量 diff + static_assert 双锚 (Case III k=1) + PSS 原表抽样旁证 |
| 2 | Case II/IV 受限确定性项的 UECM 约束构造 | 中 | R/statsmodels 双端数值对照 (verify_bounds.py 五 case 全测); R ARDL 报错行为镜像 |
| 3 | 嵌入差分 T 矩阵布局 (Harvey 形 + d 积分块) | 高 | probe_sarimax.py 源码逐元素 + ARIMA(1,1,1) 手推黄金锚 + llf 1e-8 全局校验 |
| 4 | Lyapunov 负方差静默失效 | 高 | 守卫强制 (对称化 + LLT 正定校验) + throw 测试用例 + diffuse 回退路径 |
| 5 | burn 口径差 (本库默认 d vs statsmodels 默认公式) | 中 | K4 loglik 主锚 + probe_sarimax.py 探明 + 双口径/双 burn 测试用例 |
| 6 | MLE 多局部极值 → 对照假失败 | 中 | 确定性多起始 (innovations/CSS 解) + 容差分层 1e-6~1e-8 (7C 先例) |
| 7 | R ARDL 包未安装 (7C 环境清单无此包) | 低 | readiness 前置: `install.packages("ARDL")` 落 F:\R\win-library\4.6 镜像路径 (7C urca 先例); 未装则该交叉基准降级不阻断 |
| 8 | nsim/asymptotic 文档 bug 误导 (statsmodels 稳定版) | 低 | A5 不做 MC 路径 — 无实现面; probe 记录文档 bug 留档 |

---

## 10. 实施顺序

```
M1 ARDL/PSS (先行热身, 纯复用)
  ├── 前置: probe_pss_cv.py ✅ 已完成 (2026-08-21: 布局/求值/k 语义/t 存在性五项冻结 + PSS10 新幻觉点)
  ├── Step 1: pss_critical_values.hpp + pss_response_surface.inc ✅ 已完成 (2026-08-21:
  │     gen_phase8_pss_inc.py 机器转录 1200 数字 + 5 static_assert 锚 (四表全覆盖) +
  │     probe [G] 全量 diff 0/1200 + MSVC 编译运行 30 检查全过 — 风险 #1 钉死)
  ├── Step 2: ardl_model.hpp (设计矩阵 + OLS) + verify_ardl.py ✅ 已完成 (2026-08-22:
  │     order int/list/dict 语义 + trend 全样本 1-based + A8 计账 quirk 逐式镜像;
  │     "ctt" 扩展 (statsmodels 不支持, 手工 numpy 基准); ardl_baseline.inc 自包含)
  ├── Step 3: ardl_select.hpp (IC 网格) ✅ 已完成 (2026-08-22)
  └── Step 4: uecm_bounds.hpp (重参数化 + F/t + 三值结论) + verify_bounds.py/.R
        ✅ 已完成 (2026-08-22: B 变换 ≤4.5e-15 / F 三方 ≤8.5e-14 / t vs R ARDL ≤7.1e-15;
        t 临界值 = PSS 2001 Table CII 转录 (gen_phase8_pss_t_inc.py), 档位 {10,5,2.5,1}% 独立;
        PSS10 danish 回归断言锁死)
M2 Kalman/statespace (与 M1 后半可并行, 正交)
  ├── 前置: probe_sarimax.py (T 布局 + burn 公式 + 变换语义冻结) ✅ 已完成 (2026-08-22:
  │     [A] ARIMA(1,1,1) 黄金锚逐元素 + [B] burn=d 八组 + [C] MA 负号 + [D] DK 等价
  │     + [E] 混合初始化 + [D2] concentrate_scale 陷阱; + probe_d2_const.py 常数注入
  │     两口径 + verify_drift_init.py 混合初始化实证)
  ├── Step 1: ss_model.hpp (KF + 双初始化 + 守卫) + 手推黄金锚 ✅ 已完成 (2026-08-23)
  ├── Step 2: arima_statespace.hpp (Harvey 形 + 嵌入差分 + MLE) ✅ 已完成 (2026-08-23)
  └── Step 3: innovations 常驻互检接线 + verify_statespace.py ✅ 已完成 (2026-08-23:
        gen_phase8_ss_inc.py → ss_baseline.inc; 互检字段常驻 + 集成场景 2 K1 锚)
M3: D-2 裁决后排期 (不阻塞)
收尾: test_integration_phase8 (6 场景) ✅ 已完成 (2026-08-23, 78 新测试 MSVC 全绿)
  + 三平台验证 + PHASE8_FINAL_ACCEPTANCE.md ← 当前位置
```

顺序理由: M1 与 M2 完全正交 (OLS 路线 vs 迭代滤波路线); M1 先行 (纯复用, 风险更低) 热身, M2 引入 Eigen statespace 载体; 两 probe 前置任务分别钉死两模块的最高风险断言 (风险 #1/#3) 后才冻结对应实现。

---

## 11. 开放问题 (spec 阶段裁决记录 + 残余)

| # | 问题 | 状态 | 裁决/处置 |
|---|------|------|----------|
| 1 | bounds_test 返回结构是否含 p 值 | ✅ **已裁决 + probe 冻结**: F 含 (响应面正态近似 p = 1−Φ(poly(log stat)), stat_star 切换 large/small, probe [B] 数值自证 ≤1e-16); **t 不含** (NaN + has_p=false, 无 t 响应面, probe [D]) | §2.3/§2.4 |
| 2 | ARDL trend 参数枚举 vs 字符串 | ✅ **已裁决 (本 spec)**: 字符串 "n"/"c"/"ct"/"ctt" (statsmodels 直译, §1.4-4 惯例); PSS Case 用 `enum class PssCase`; "ctt" 不支持 bounds (throw); case-trend 相容校验镜像 R ARDL 报错 | §2.1/§2.3 |
| 3 | loglikelihood_burn 默认值 | ✅ **已裁决 (本 spec)**: 默认 = d (嵌入式差分口径), 显式参数可覆盖; statsmodels None 默认公式 probe_sarimax.py 探明 — 一致则镜像, 不一致则文档声明差异 + 双 burn 测试 | §3.2 |
| 4 | D-2 (M3 战略选项: BVAR/TVP-VAR R 桥 vs KDE/KNN+TE) | ⏳ 待 Scott 裁决 | 不阻塞 M1/M2 (§4 占位) |
| 5 | Narayan k≤7 边界 (调研 B-) | ✅ 已消解: A5 不采 Narayan 分档表 — 无实现面; 调研报告留档 | §2.4 |
| 6 | statsmodels nsim 文档 bug | ✅ 已消解: A5 不做 MC 路径; probe 留档 (风险 #8) | §9 |

---

## 12. 关联文档

- [ADR-020: v1.8 ARDL/PSS + Kalman 实施边界 (12 项)](../../decisions/ADR-020_V18_ARDL_KALMAN_BOUNDARY.md) (Accepted 2026-08-21)
- [ADR_INDEX.md](../../decisions/ADR_INDEX.md) (ADR-017/018/019/020)
- [PHASE8_RESEARCH.md](../../research/PHASE8_RESEARCH.md) v1.1 (双域调研 + 审计二轮 + 幻觉点详表 19 条)
- [PHASE8_SCOPE_RESEARCH.md](../../research/PHASE8_SCOPE_RESEARCH.md) v1.3 (D-1 裁决 + 双免费午餐支撑材料)
- [DEVELOPMENT_WORKFLOW.md](../../DEVELOPMENT_WORKFLOW.md) v1.1 §2.0/§4.2 (阶段 0 与 R/G 双门禁)
- [PHASE7C_SPEC.md](../phase7/PHASE7C_SPEC.md) v1.2 (接口总则/工程约束模板来源; v1.7 复用基础)
- [PHASE7C_M0_READINESS_CHECKLIST.md](../phase7/PHASE7C_M0_READINESS_CHECKLIST.md) (环境预检惯例; 本期 readiness 同模式编制)
