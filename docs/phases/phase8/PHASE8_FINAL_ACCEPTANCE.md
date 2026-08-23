# Phase 8 v1.8 最终验收报告 — ARDL 边限检验与 Kalman 状态空间模块

> **文档类型**: 最终验收报告 (Final Acceptance Report)
> **审计对象**: Phase 8 v1.8 完整 scope (M1 ARDL/PSS + M2 Kalman/statespace + §5 端到端集成; M3 = D-2 待裁决, 不阻塞)
> **关联文档**:
>   - [PHASE8_SPEC.md](./PHASE8_SPEC.md) (执行规格, §8 验收标准逐项对照)
>   - [ADR-020 v1.8 ARDL/PSS + Kalman 实施边界 (12 项)](../../decisions/ADR-020_V18_ARDL_KALMAN_BOUNDARY.md)
>   - [DEVELOPMENT_LOG.md](../../DEVELOPMENT_LOG.md) Phase 8 实施记录
>   - [PHASE8_RESEARCH.md](../../research/PHASE8_RESEARCH.md) v1.1 / [PHASE8_SCOPE_RESEARCH.md](../../research/PHASE8_SCOPE_RESEARCH.md) v1.3
> **验收日期**: 2026-08-23
> **验收员**: Scott (self-review, solo developer 工作流)
> **验收标准**: 三平台全量回归 + statsmodels 0.14.6 / R ARDL 0.2.5 双端基准分层容差 + 幻觉点 20 编号全覆盖 + ADR-020 12 项零偏离 + v1.7.1 零破坏

---

## 1. 验收基本信息

| 项 | 值 |
|----|-----|
| 项目名称 | Cpp_Hub |
| 审计阶段 | Phase 8 v1.8 (ARDL 边限检验 + Kalman 状态空间 — M1 边限检验四件套 + M2 KF/精确似然 MLE + 集成六场景) |
| 前置版本 | v1.7.1 (2465 主控 / 2447×2 A/B 站) |
| 后置版本 | v1.8 (**2543/2543** 主控站 MSVC / **2525/2525** ×2 A/B 站 GCC, 三平台全绿) |
| Phase 8 新增测试 | **78 个** (ardl_model 14 + ardl_select 10 + uecm_bounds 18 + statespace 16 + arima_statespace 14 + integration 6) |
| 实现 commit | 351ac98 (单 commit 交付 M1+M2+集成; ADR-020/调研文档先行 5aa7b2a) |
| 验收方法 | TDD + 双 probe 前置探明 (probe_pss_cv.py / probe_sarimax.py) + statsmodels/R ARDL 三方交叉 (1e-15~1e-8 分层) + bundle 中继三平台验证 |
| 里程碑 | 前置 probe (08-21/22) → M1 四步 (08-21/22) → M2 三步 (08-22/23) → 集成+回归 (08-23) → 三平台+CI (08-23) |

> 主控站与 GCC 18 个测试差额为平台专属用例 (Windows/MSVC-only), 7B/7C 期已知口径, 非功能差异。

---

## 2. 全量回归测试结果

### 2.1 三平台矩阵 (2026-08-23 实测)

| 平台 | 编译器 | 测试总数 | 通过数 | 失败数 | 测试耗时 | 验证方式 | 状态 |
|------|--------|---------|--------|--------|---------|---------|------|
| 主控站 (Win10) | MSVC 19.50 (VS2026) | 2543 | 2543 | 0 | ~142 sec (-j8) | 增量构建 + 全量 ctest | ✅ |
| A 工作站 (`scott-lau-NEX.local`) | GCC 13.3 | 2525 | 2525 | 0 | 54.13 sec (-j8) | bundle 增量 ff 28e3028→351ac98 + rebuild + ctest | ✅ |
| B 工作站 (`scott-lau-GTR-Pro.local`) | GCC 13.3 | 2525 | 2525 | 0 | 47.75 sec (-j8) | 同法 (mDNS 通道) | ✅ |

### 2.2 GitHub Actions CI (run #76, 351ac98)

| Job | 结果 |
|-----|------|
| Build & Test (windows-2022, MSVC) | ✅ success |
| Build & Test (ubuntu-latest, GCC) | ✅ success |
| C ABI (windows-2022) | ✅ success |
| C ABI (ubuntu-latest) | ✅ success |

> 前置 run #75 (e80617d, v1.7.1 文档) success; run #76 首次承载全量 2543 用例双平台绿 (2026-08-23 确认)。

### 2.3 测试数演进

| 轮次 | 内容 | 新增测试 | 主控站累计 | A/B 站累计 |
|------|------|---------|-----------|-----------|
| M1 | pss_cv (30 检查) + ardl_model 14 + ardl_select 10 + uecm_bounds 18 | 42 | 2507 | 2489 |
| M2 | statespace 16 + arima_statespace 14 | 30 | 2537 | 2519 |
| 收尾 | test_integration_phase8 6 | 6 | **2543** | **2525** |

### 2.4 跨平台验证执行细节

- A/B 站沿用 **bundle 中继工作流**: 主控站 `git bundle create tmp/phase8.bundle 28e3028..HEAD` (0.18MB) + scp, 站内 `git fetch <bundle> HEAD` + ff-only → 351ac98, 增量 rebuild + ctest 零外网依赖。
- v1.7.1 遗留订正 (9771dd6, ARIMA 预算余量 250×→238× 注释) 随本 bundle 一并送抵 A/B 站 (此前两站止于 28e3028)。
- Phase 8 新增 78 用例 GCC 双站全数运行, 三平台零数值偏差。

---

## 3. 实施范围汇总

### 3.1 新增头文件 (8 个, `cpphub::v1::timeseries` 命名空间, Eigen 载体 = cpphub_timeseries_mat)

| 文件 | 行数 | 用途 |
|------|------|------|
| `ardl/ardl_model.hpp` | 220 | ARDL(p,q) OLS + trend 四式 (n/c/ct/ctt) + hold_back + A8 计账镜像 |
| `ardl/ardl_select.hpp` | 163 | IC 网格阶数选择 (同一 hold_back 样本, (maxlag+1)·(maxorder+1)^k 联合网格) |
| `ardl/uecm_bounds.hpp` | 434 | UECM 重参数化 + bounds F/t 检验 (五 Case, 三值结论) |
| `ardl/pss_critical_values.hpp` | 148 | PSS 响应面临界值 (查表 + p 值正态近似 + static_assert 锚) |
| `ardl/pss_response_surface.inc` | 440 | F 表响应面系数 (机器转录 1200 数字, gen_phase8_pss_inc.py) |
| `ardl/pss_t_bounds.inc` | 146 | t 表 (PSS 2001 Table CII, gen_phase8_pss_t_inc.py; 档位 {10,5,2.5,1}%) |
| `statespace/ss_model.hpp` | 290 | StateSpaceModel (DK 记号聚合体) + KF 更新形 + 双初始化 + Lyapunov 守卫 |
| `statespace/arima_statespace.hpp` | 646 | ARMA Harvey 形 + ARIMA 嵌入差分 + 精确似然 MLE + innovations 常驻互检 |

### 3.2 入库 fixtures 脚本 (verify/probe/gen 全链落盘, 11 个 + 数据 4 个)

| 脚本 | 用途 |
|------|------|
| `probe_pss_cv.py` | M1 前置: CV 布局/求值/k 语义/t 存在性五项冻结 + [G] .inc 全量 diff |
| `gen_phase8_pss_inc.py` / `gen_phase8_pss_t_inc.py` | statsmodels/R → .inc 机器转录生成器 (1200 数字 / t 表) |
| `verify_ardl.py` | statsmodels ARDL 语义三探明 (order/trend/A8) + 手工重建对拍 → ardl_baseline.inc |
| `probe_r_ardl.R` + `gen_r_ardl_bounds.R` → `r_ardl_bounds.csv` | R ARDL 0.2.5 bounds F/t 三方交叉基准 |
| `export_pss2001_t_bounds.R` → `pss2001_t_bounds.csv` | R ARDL t 界转录 (= PSS 2001 Table CII 原文) |
| `probe_sarimax.py` | M2 前置: T/R/Z 布局 + burn=d 八组 + MA 负号 + DK 等价 + 混合初始化 + concentrate_scale 陷阱 |
| `probe_d2_const.py` / `verify_drift_init.py` | 常数注入两口径 (d=0 state_intercept vs d≥1 drift) + 混合初始化实证 |
| `gen_phase8_ss_inc.py` → `ss_baseline.inc` | SARIMAX llf/params 黄金基准 (T=120/300 双长度, 口径注释内嵌) |
| `ardl_smoke_data.csv` / `ardl_statsmodels_baselines.json` | verify_ardl.py 数据与基准落盘 |

> SPEC §1.1 计划名与实际落地的功能等价改名: verify_bounds.py/.R 并入 verify_ardl.py + probe_r_ardl.R 链路; verify_statespace.py 并入 gen_phase8_ss_inc.py + probe_sarimax.py [D]; gen_phase8_inc.py → gen_phase8_ss_inc.py。新增 pss_t_bounds.inc / probe_d2_const.py / verify_drift_init.py 三个计划外产物 (Step 4 t 表与常数注入两口径探明的直接产出)。

### 3.3 新增测试 (6 套件 78 用例 + 2 个自包含 .inc 基准)

ardl_baseline.inc (2142 行, 数据+基准全嵌) + ss_baseline.inc (180 行) — 测试自包含, 不依赖运行时外部文件 (7C 模式)。

---

## 4. 数值基准验证 (SPEC §1.3 全表实测)

| 基准 | SPEC 容差 | 实测 | 状态 |
|------|----------|------|------|
| ARDL 系数/σ̂²/IC vs statsmodels ARDL.fit | 1e-10 | 手工重建 1e-15 量级 (verify_ardl.py S1/S2 全量自证) | ✅ |
| ardl_select_order 最优阶 + IC vs statsmodels | 1e-10 | 同上 (含网格维度/argmin/同一 n_obs_used 断言) | ✅ |
| UECM 重参数化 vs statsmodels UECM.from_ardl | 1e-10 | ≤4.5e-15 | ✅ |
| bounds F 统计量 vs statsmodels bounds_test | 1e-10 | ≤8.5e-14 (三方: 本库/statsmodels/R ARDL) | ✅ |
| bounds F + CV + p vs R ARDL bounds_f_test | 1e-8 | ≤8.5e-14 | ✅ |
| bounds t 统计量 + CV vs R ARDL bounds_t_test | 1e-8 | ≤7.1e-15 | ✅ |
| CV 查表 vs statsmodels crit_vals 表键直查 | 精确相等 | probe [G] 全量 diff **0/1200** | ✅ |
| p 值 vs statsmodels _pss_pvalue (跨 stat_star) | 1e-14 | ≤1e-16 (probe [B]) | ✅ |
| 响应面系数转录 (1200 数字) | 精确相等 | 全量 diff 0 | ✅ |
| 响应面 static_assert 锚 (双锚+扩展) | 精确相等 | 5 锚 (四表全覆盖, 编译期) | ✅ |
| PSS 原表交叉印证 (Case III k=1 5%) | ~0.05 文档级 | F 4.94/5.73, t −2.86/−3.22 (Statalist Kripfganz 独立印证) | ✅ |
| ARMA(1,1) Harvey 形 T/R/Z 手推 | 精确相等 | test_statespace 黄金锚 (含 ARIMA(1,1,1) 逐元素) | ✅ |
| AR(1) 平稳 P₀ = σ²/(1−φ²) 手推 | 1e-12 | 1e-12 | ✅ |
| KF 前两步 v/F/loglik 手推 (AR(1)) | 1e-12 | 1e-12 | ✅ |
| KF loglik vs statsmodels SARIMAX llf (burn 对齐) | 1e-8 | 2.81e-9 @T=120 (probe [D]; DK 更新形等价 diff 3.7e-14) | ✅ |
| ARIMA statespace MLE 参数 vs SARIMAX | 1e-6~1e-8 | 1e-6 层 (θ 方向弱识别: statsmodels 同 DGP T=120 θ̂=0.6199 vs 本库 T=300 0.567, 向真值收敛; 三路径 1e-6 对拍锚定) | ✅ |
| innovations 常驻互检 gap | ≤1e-5 告警线 | 2.6e-8 @T=300 (集成场景 2 K1 锚) | ✅ |
| Lyapunov 守卫 φ=1.01 → throw | 行为断言 | throw + diffuse 回退逐位正确 (集成场景 6) | ✅ |

---

## 5. 幻觉点核查 (20 编号全覆盖)

### 5.1 域 A: ARDL/PSS (PSS1-PSS10)

| ID | 核查结果 | 落点 |
|----|---------|------|
| PSS1 ardlverse CV 错误黑名单 | ✅ 不采信即核查 (verify 脚本头部声明) | verify_ardl.py 头注 |
| PSS2 KS 2020 出处 = OBES 82(6) | ✅ 文档引用核对 | 调研报告 v1.1 |
| PSS3 t 界仅 Case I/III/V | ✅ Case II/IV → throw + 三 case 逐值 ≤7.1e-15 | test_uecm_bounds 18 用例 |
| PSS4 Narayan 不采 | ✅ A5 无实现面 | ADR-020 A5 |
| PSS5 分位档四档语义 (1%/0.1% ≠ PSS 2.5%) | ✅ probe 逐列核对, 四档语义冻结 | probe_pss_cv.py [B] |
| PSS6 k = force 变量数 (不含 y 滞后) | ✅ k=1 vs q=2 同 CV 键用例 | test_uecm_bounds |
| PSS7 F 版 = 二次型/约束数 (非 χ²) | ✅ F 逐值三方 ≤8.5e-14 + n_restrictions 五 case 断言 | verify_ardl.py / R 链 |
| PSS8 UECM 样本损失 = max(p, max q_j) | ✅ n_obs_used 断言 + select 同一样本 | test_ardl_select / 集成场景 1 |
| PSS9 HC/HAC 不在表覆盖内 | ✅ API 无 cov_type 参数 (不存在即核查) | uecm_bounds.hpp |
| PSS10 statsmodels bounds_test 查表 k 偏移 +1 | ✅ **danish 回归断言锁死**: 本库 CV = 键 (3,3) [3.2290,4.3223] ≠ statsmodels 返回 (键 (4,3) [2.8772,4.0105]) | probe_pss_cv.py [E] + test_uecm_bounds |

### 5.2 域 B: Kalman/statespace (SS1-SS10)

| ID | 核查结果 | 落点 |
|----|---------|------|
| SS1 approximate diffuse κ=1e6 ≠ exact | ✅ 枚举无 Exact 项; 仅作初始化 | ss_model.hpp SSInitSpec |
| SS2 Harvey 形 vs Hamilton 形 | ✅ ARMA(1,1) 手推三矩阵精确 + ARIMA(1,1,1) 黄金锚逐元素 (dT=dR=0) | test_statespace / probe_sarimax [A] |
| SS3 R 非方阵 m×r, Q r×r | ✅ 维度校验 throw + ARIMA R 构造断言 | ss_model.hpp |
| SS4 Lyapunov 无守卫 → 负方差静默失效 | ✅ 对称化 + LLT 正定校验 throw; φ=1.01 用例 + diffuse 回退 | test_statespace + 集成场景 6 |
| SS5 双差分口径似然不同 | ✅ 两 loglik 差异显式记录, 禁止混比 (gap=O(1/κ) 注释内嵌基准) | test_arima_statespace + 集成场景 5 |
| SS6 H=0, F 正定性全靠 P | ✅ H 零阵断言 + F_t ≤ 0 throw | ss_model.hpp |
| SS7 Chandrasekhar 不做 | ✅ v1.9 边界, 审计锚定 v0.12.0 | ADR-020 K2 |
| SS8 似然常数 (单变量 log 2π) | ✅ loglik 手推锚含常数项 | test_statespace |
| SS9 diffuse likelihood 不作估计目标 | ✅ approximate diffuse 仅初始化, 文档声明 | arima_statespace.hpp |
| SS10 UCM 一元专属 | ✅ scope 不做清单声明 | SPEC §8.5 |

### 5.3 实施期新发现幻觉点 (probe/verify 产出, SPEC §7 之外)

| ID | 发现 | 处置 |
|----|------|------|
| A8 | statsmodels ARDL 计账 quirk: nobs/llf/sigma2/IC 用 n_acct = T − hb_early (hb_early 在 _check_lags 固化 = 显式 hb 或 max(y_lags)=p, **exog 滞后不并入**); max(q_j) > p 时 n_acct ≠ n_act (T=250, p=1, q=2: nobs 249 vs 真实 248); bse 用真实 df | 本库逐式镜像 (llf/IC/σ² 同口径, bse 真实 df), devlog 钉死精确机制 |
| order 语义 | statsmodels order list = **显式滞后列表**套所有变量 (非逐变量最大阶!); int = 全展开; 逐变量必须 dict | 本库 `q` 向量原生 per-variable, 不镜像 list 陷阱 |
| trend 列 | 全样本 1-based 再切 hold_back (t = hb+1..T), 与 VAR 估计样本内重编号不同 — 跨模块陷阱 | verify_ardl.py S1 实证 |
| t 表档位 | t 表 {10,5,2.5,1}% ≠ F 表 {10,5,1,0.1}% (PSS9 对偶面) | pss_t_bounds.inc 独立档位 |
| MA 变换负号 | θ = −constrain_stationary_univariate(u) (statsmodels 实源, 小值区 θ≈+u 掩盖) | probe_sarimax [C] 往返 5.55e-17 |
| concentrate_scale | True 时 k_params=p+q 且 params[-1] 是 θ — 对拍必须同口径 | K4 非集中化 + probe [D2] |

---

## 6. ADR-020 对齐审计 (12 项零偏离)

| # | 决策 | 实施对照 | 状态 |
|---|------|---------|------|
| A1 | 响应面多项式路线 | pss_critical_values.hpp + 1200 数字机器转录 + p 值正态近似 (stat_star 切换) | ✅ |
| A2 | ARDL(p,q) 单方程, 排除 ardlverse 系 | ardl_model.hpp 仅单方程; 无 panel/quantile/Fourier/NARDL | ✅ |
| A3 | 网格 argmin IC, 同一 hold_back | ardl_select.hpp (hb_grid = max(maxlag, maxorder) 统一) | ✅ |
| A4 | t 仅 Case I/III/V | bounds_t_test II/IV → throw; t 表 PSS 2001 Table CII 转录 | ✅ |
| A5 | 不采 Narayan / 不做 MC | 无分档表; asymptotic 响应面天然插值 | ✅ |
| A6 | UECM 含同期 Δx_j | causal=false 恒定, true → throw (v1.9) | ✅ |
| A7 | 仅 nonrobust | API 无 cov_type | ✅ |
| K1 | DK 记号 + 常规 KF + 双初始化 + Harvey 形 + 嵌入差分 + 互检 | ss_model + arima_statespace 全落地; SARIMAX 默认口径四件套对齐 | ✅ |
| K2 | 推迟 v1.9+ (exact diffuse/平滑器/缺失/多元/Chandrasekhar/sqrt filter) | 均未实现 (scope 不做清单 §8.5) | ✅ |
| K3 | 复用 ADR-018 SLSQP + Monahan 变换 | arima_statespace_mle 多起始 {innovations 解, CSS 解} | ✅ |
| K4 | σ² 非集中化 | 全参数优化; loglik 主锚 + 口径差异声明 | ✅ |
| K5 | statespace 挂 cpphub_timeseries_mat (Eigen 隔离) | CMakeLists 同 var/cointegration 模式; 单变量纯标量头零反向 include | ✅ |

---

## 7. Scope 边界核查 (SPEC §8.5)

- ✅ 未实现清单全部遵守: panel/quantile/Fourier/NARDL / glob 搜索 / Narayan 分档 / MC 模拟 / HC/HAC Wald / causal ARDL / exact diffuse / 平滑器 / 缺失数据 / 多元 statespace / Chandrasekhar / SARIMA 季节 — 均无实现痕迹 (grep 验证 + API 面)。
- ✅ 单变量纯标量头 (arima/garch/unit_root/midas) 零 ardl/statespace 反向 include。
- ✅ M3 未实现且不阻塞 (D-2 待裁决声明, §4 占位)。
- ✅ "ctt" 扩展为 SPEC 裁决 #2 计划内 (statsmodels 不支持, 手工 numpy 基准); bounds 拒绝 ctt (case-trend 相容校验)。

## 8. 零破坏与工程收尾 (v1.7.1 兼容)

- ✅ §6.5 additive-only: Phase 7C 既有文件仅换行符差异 (gen_phase7c_m3_baseline.py / test_variance_ratio_test.cpp) 经 `git diff --ignore-cr-at-eol` 判定后 checkout 还原; 唯一实质触碰 = test_arima_model.cpp 2 行注释 (v1.7.1 遗留 250×→238× 订正, 9771dd6)。
- ✅ .gitignore 补 tmp/ scratch 目录条目。
- ✅ Eigen 隔离: 新模块挂 cpphub_timeseries_mat (INTERFACE), 零新增链接依赖。
- ✅ C ABI 不受影响 (statespace/ardl 不导出)。

## 9. 遗留问题与后续

| # | 项 | 状态 | 处置 |
|---|----|------|------|
| 1 | D-2 (M3 战略选项: BVAR/TVP-VAR R 桥 vs KDE/KNN+TE) | ✅ **已裁决 (2026-08-23, Scott): (b) KDE/KNN + KSG/TE, 载体 v1.9 (Phase 9)** | 不回填 Phase 8; 支撑材料 SCOPE_RESEARCH §5.2; BVAR 解析后验降 v1.9 候选池, TVP-VAR 挂靠平滑器前置 |
| 2 | v1.9 候选池 | 已登记 | **KDE/KNN + KSG/TE (D-2 裁决, Phase 9 M1)** + causal ARDL / glob 搜索 / HC/HAC Wald / exact diffuse / 平滑器 / 多元 statespace / BVAR 解析后验 / Narayan 重评 (一手来源出现时) |
| 3 | A8 计账 quirk 上游反馈 | 可选 | statsmodels issue 候选 (nobs 与真实估计样本不一致), 低优先 |
| 4 | tag v1.8.0 | ✅ 已完成 (2026-08-23) | 打于 351ac98 并推送 (本地+远程); CI run #76 四 job success |

## 10. 验收结论

**Phase 8 v1.8 验收通过。** 三平台 2543/2525×2 全绿; SPEC §1.3 数值基准全表达成 (含 4 项精确相等); 幻觉点 20 编号全覆盖 + 6 个实施期新发现钉死; ADR-020 12 项零偏离; v1.7.1 零破坏。M1+M2 双免费午餐兑现: statsmodels/R 双端基准成熟、UECM=OLS+Wald 与 Kalman=Eigen+innovations 复用率预判成立, 单 commit (351ac98, 11561 行净增) 周级交付。
