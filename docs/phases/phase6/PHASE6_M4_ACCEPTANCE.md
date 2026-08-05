# Phase 6 v1.5 M4 实施验收报告

> **文档类型**: 实施验收报告 (Implementation Acceptance Report)
> **审计对象**: v1.5 M4 - Bootstrap (4 种) + 端到端集成 + 跨语言验证基础设施
> **关联文档**:
>   - [PHASE6_IMPLEMENTATION_PLAN.md](./PHASE6_IMPLEMENTATION_PLAN.md) §6
>   - [PHASE6_ECONOMETRICS_SPEC.md](./PHASE6_ECONOMETRICS_SPEC.md)
> **验收日期**: 2026-08-05
> **验收员**: Scott (self-review, 符合项目 solo developer 工作流)
> **验收标准**: 主控站全量测试通过 + 排幻觉点 12 项全部验证 + R/Python 跨语言对照通过

---

## 1. 验收基本信息

| 项 | 值 |
|----|-----|
| 项目名称 | Cpp_Hub |
| 审计阶段 | Phase 6 v1.5 M4 (经典参数计量模块 - Bootstrap + 端到端) |
| 提交版本 | `e64e961` feat(v1.5 M3+M4): 偏差项补齐 - AR(1)/AR(2)检验, CUE优化, Bootstrap完整实现, R/Python跨语言验证 |
| 跨平台修复版本 | `b278151` fix(v1.5 M4): cross-platform RNG consistency - Philox4x64+box_muller 替换 std::mt19937_64+std::normal_distribution |
| 验收方法 | TDD 实现 + 三平台全量测试 + R 排幻觉点验证 + Python 跨语言对照 |
| 前置条件 | v1.5 M3 已通过 (1699/1699 主控站, 1717/1717 A/B 站) |
| 后置条件 | M4 验收通过, v1.5 经典参数计量模块 (M1+M2+M3+M4) 全部完成 |

---

## 2. 测试结果汇总

### 2.1 三平台测试结果

| 平台 | 编译器 | 测试总数 | 通过数 | 失败数 | 编译耗时 | 测试耗时 | 状态 |
|------|--------|---------|--------|--------|---------|---------|------|
| 主控站 (Windows 10) | MSVC 19.43 | 1767 | 1767 | 0 | - | - | ✅ 通过 |
| A 工作站 (Ubuntu 24.04) | GCC 13.3.0 | 1767 | 1767 | 0 | 30 sec | 41 sec | ✅ 通过 |
| B 工作站 (Ubuntu) | GCC 13.3.0 | 1767 | 1767 | 0 | 31 sec | 39 sec | ✅ 通过 |

### 2.2 v1.5 测试规模

| 项 | 值 |
|----|-----|
| v1.5 新增测试文件 | 25 个 `.cpp` |
| v1.5 新增测试用例 | 399 个 (TEST_F/TEST 宏计数) |
| spec 要求测试数 | 185 个 |
| 超出 spec | 214 个 (+115.7%) |
| 三平台编译 | ✅ 全部通过 |
| 三平台测试 | ✅ 全部通过 (1767/1767) |

**注**: A/B 站 fresh clone commit `b278151` 后 1767 个测试全部通过, 0 失败。M1-M4 累计 399 个新测试用例, 加上 Phase 3-5 遗留测试共 1767 个。

### 2.3 M4 新增测试明细

| 测试文件 | 测试数 | 覆盖范围 |
|---------|--------|---------|
| `test_bootstrap_base.cpp` | 8 | Bootstrap 基类, BCa 置信区间, 权重分布 |
| `test_paired_bootstrap.cpp` | 10 | 配对 Bootstrap, CI 覆盖概率 |
| `test_wild_bootstrap.cpp` | 10 | Wild Bootstrap, Rademacher/Mammen/Webb6 |
| `test_block_bootstrap.cpp` | 12 | Block Bootstrap, Stationary/Circular/NonOverlapping, Politis-White |
| `test_cluster_bootstrap.cpp` | 14 | Cluster Bootstrap, 面板数据 entity/time 维度 |
| `test_integration_phase6.cpp` | 17 | 端到端: OLS→HC/HAC→Wald→Bootstrap, MLE→QMLE→LR, GMM→J-test→Bootstrap |
| **M4 小计** | **71** | |

### 2.4 M3 偏差项补齐 (本次提交)

| 补齐项 | 测试文件 | 新增测试数 |
|--------|---------|-----------|
| AR(1)/AR(2) 序列相关检验 | `test_arellano_bond.cpp` (扩展) | +3 |
| CUE 完美拟合处理 | `test_gmm_cue_iterated.cpp` (修复) | 0 (修复现有) |
| **小计** | - | **3** |

---

## 3. M4 实施范围验收

### 3.1 Bootstrap 基础设施

| 组件 | 实现文件 | 排幻觉点 | 状态 |
|------|---------|---------|------|
| Bootstrap 基类 | `resampling/bootstrap_base.hpp` | E12 | ✅ |
| 配对 Bootstrap | `resampling/paired_bootstrap.hpp` | - | ✅ |
| Wild Bootstrap | `resampling/wild_bootstrap.hpp` | E12 | ✅ |
| Block Bootstrap | `resampling/block_bootstrap.hpp` | - | ✅ |
| Cluster Bootstrap | `resampling/cluster_bootstrap.hpp` | - | ✅ |

### 3.2 Bootstrap 算法验收

| 算法 | 关键参数 | 排幻觉点 | 状态 |
|------|---------|---------|------|
| 配对 Bootstrap | 有放回采样 N 对, V_boot = (1/(B-1))Σ(θ̂*_b-θ̄*)(θ̂*_b-θ̄*)' | - | ✅ |
| Wild Bootstrap (Rademacher) | y*_i = x_i'β̂ + v_i·ε̂_i, v=±1 w.p. 0.5 | E12 (CGM 2008 推荐) | ✅ |
| Wild Bootstrap (Mammen) | v=(1±√5)/2 w.p. (5±√5)/10 | E12 | ✅ |
| Wild Bootstrap (Webb6) | 6 点分布 {-√(3/2),-1,-√(1/2),√(1/2),1,√(3/2)} | E12 | ✅ |
| Block Bootstrap (Stationary) | 块长 ~ 几何分布 (Politis-Romano 1994) | - | ✅ |
| Block Bootstrap (Circular) | 环形索引, 固定块长 | - | ✅ |
| Block Bootstrap (NonOverlapping) | 固定块长不重叠 | - | ✅ |
| Politis-White 自动块长 | 基于自相关函数的自动选择 (Politis-White 2004) | - | ✅ |
| Cluster Bootstrap | 从 G 个聚类有放回采样 G 个聚类, G<20 警告 | - | ✅ |
| BCa 置信区间 | bias correction + acceleration (Davison-Hinkley §5.3) | - | ✅ |

### 3.3 M3 偏差项补齐

| 补齐项 | 实现细节 | 排幻觉点 | 状态 |
|--------|---------|---------|------|
| AR(1)/AR(2) 检验 | m_j = ΣΔε̂_it·Δε̂_{i,t-j} / sqrt(Σ(Δε̂_it)²·Σ(Δε̂_{i,t-j})²), 渐近 N(0,1), p值 via erfc | E11 | ✅ |
| CUE 数值优化 | Nelder-Mead simplex 最小化 J(β)=N·ḡ'Ŝ⁻¹ḡ, Ŝ奇异时返回 J=0 | - | ✅ |
| CUE 完美拟合 | β使ε≈0时, ḡ≈0 → J=0 (而非 max()) | - | ✅ |

### 3.4 端到端集成测试

| 场景 | 流程 | 测试数 | 状态 |
|------|------|--------|------|
| OLS → HC/HAC → Wald → Bootstrap | 因子回归 + 稳健推断 + Bootstrap CI | 4 | ✅ |
| MLE → QMLE Sandwich → LR → Bootstrap | 似然推断 + QMLE 诊断 + Bootstrap | 3 | ✅ |
| GMM → J-test → Bootstrap | 矩估计 + 过度识别检验 + Bootstrap | 3 | ✅ |
| Arellano-Bond → AR(1)/AR(2) → Sargan | 动态面板 + 序列相关 + 过度识别 | 4 | ✅ |
| Bootstrap 基础设施 | BCa/percentile CI + 权重分布 + 种子可复现 | 3 | ✅ |
| **合计** | - | **17** | ✅ |

---

## 4. 排幻觉点验证 (verify_econometrics.R)

### 4.1 验证结果汇总

| 验证项 | 模块 | 结果 | 备注 |
|--------|------|------|------|
| E1 | OLS 截距 | PASS | model.matrix 第一列全为 1 |
| E2 | HC1 自由度调整 | PASS | HC1 = N/(N-K)·HC0 (tolerance 1e-8) |
| E3 | HC2 leverage | PASS | h_i = diag(X(X'X)^{-1}X') == hatvalues() |
| E4 | Newey-West 默认 lag | PASS | lag = floor(bwNeweyWest(fm)) (NW 1994 自动带宽, 非 1987 经验法则) |
| E5 | Bartlett 权重 | PASS | w[l] = 1-l/(L+1), 非 1-l/L |
| E6 | vcovCL 小样本调整 | PASS | G/(G-1)·(N-1)/(N-K), Grunfeld G=10 (as.numeric) |
| E7 | glm IRLS vs Newton-Raphson | PASS | canonical link 等价 (check.attributes=FALSE) |
| E8 | sandwich bread/meat 分解 | PASS | bread=n·(X'WX)^{-1}, meat=crossprod(estfun)/n (OPG) |
| E9 | waldtest 默认 F 检验 | PASS | lm 对象返回 F 列 |
| E10 | gmm Ŝ tangent matrix | VERIFIED | gmm 1.6-4+ type="twoStep" (大小写敏感) |
| E11 | plm::pgmm 工具变量 | VERIFIED | plm 2.6+ multi-part formula (弃用 dynformula) |
| E12 | 完全拟合 Ŝ 奇异处理 | PASS | HC0 数值零 (abs<1e-20), R 优雅处理不报错 |

**总计**: 10 PASS + 2 VERIFIED = 12/12 全部通过

### 4.2 本次修复的排幻觉点偏差

| 排幻觉点 | 原状态 | 修复内容 |
|---------|--------|---------|
| E2 | DEVIATION | tolerance 放宽到 1e-8 (Longley 共线性浮点误差 ~1e-10) |
| E3 | DEVIATION | tolerance 放宽到 1e-8 (同上) |
| E4 | ERROR | `bwNeweyWest` API 变化; 改为验证 `NeweyWest` 默认 lag = `floor(bwNeweyWest(fm))` (NW 1994 自动带宽); Longley `prewhite=TRUE` 失败, 用 `prewhite=FALSE` |
| E5 | ERROR | `kweights` 签名变化; 改为 `kweights(x=(0:L)/(L+1), kernel="Bartlett")` |
| E6 | DEVIATION | `Grunfeld$firm` 是 data.frame, `unique()` 返回行数; 改用 `as.numeric()` 获取真实聚类数 G=10 |
| E7 | DEVIATION | `all.equal` names 不匹配; 加 `check.attributes=FALSE` |
| E8 | SKIP→PASS | car 不可用→复用 E7 的 `fm_glm`; meat 公式从 Pearson 残差改为 `crossprod(estfun)/n` (OPG); bread 含 n 因子 |
| E10 | SKIP→VERIFIED | 安装 gmm 包; `type="twostep"`→`"twoStep"` (gmm 1.6-4+ 大小写敏感) |
| E11 | ERROR→VERIFIED | `dynformula` 弃用→multi-part formula `lag(log(emp),1:2) + ... \| lag(log(emp),2:99)` |
| E12 | VERIFIED→PASS | `all(V==0)` 浮点误差失败→`all(abs(V)<1e-20)` 数值零判断 |

### 4.3 spec 同步更新

以下排幻觉点已在 [PHASE6_ECONOMETRICS_SPEC.md](./PHASE6_ECONOMETRICS_SPEC.md) §7.3 和 [PHASE6_IMPLEMENTATION_PLAN.md](./PHASE6_IMPLEMENTATION_PLAN.md) 中同步更新:

| 排幻觉点 | 更新内容 |
|---------|---------|
| E4 | NW 1994 自动带宽 (基于 s1/s0), 非 Andrews 1991; Longley prewhite=TRUE 失败 |
| E6 | cadjust=TRUE 默认, 两层调整 G/(G-1) 和 (N-1)/(N-K); Grunfeld$firm 需 as.numeric() |
| E8 | bread = (X'WX)^{-1}·n (含 n 因子); meat = crossprod(estfun)/n (OPG, 非Pearson残差外积) |
| E10 | gmm 1.6-4+ type 参数大小写敏感 ("twoStep" 非 "twostep") |
| E11 | plm 2.6+ 弃用 dynformula, 改用 multi-part formula; AR 检验为 $m1/$m2 htest 对象 |

---

## 5. 跨语言验证基础设施

### 5.1 R 基准生成脚本 (7 个)

| 脚本 | 输出 JSON | 数据集 | 状态 |
|------|----------|--------|------|
| `generate_ols_baselines.R` | `ols_hc_baseline.json` | Longley (16 obs) | ✅ EXIT 0 |
| `generate_hac_baselines.R` | `hac_baseline.json` | Longley + AR(1) | ✅ EXIT 0 |
| `generate_mle_baselines.R` | `mle_baseline.json` | Spector-Mazzeo (32 obs) + warpbreaks (48 obs) | ✅ EXIT 0 |
| `generate_gmm_baselines.R` | `gmm_baseline.json` | 合成线性 IV (Hayashi §3.5) | ✅ EXIT 0 |
| `generate_arellano_bond_baselines.R` | `arellano_bond_baseline.json` | EmplUK (140 firms, 1976-1984) | ✅ EXIT 0 |
| `generate_bootstrap_baselines.R` | `bootstrap_baseline.json` | 合成 AR(1) + 法学院数据 | ✅ EXIT 0 |
| `generate_grunfeld_baseline.R` | `grunfeld_cluster_baseline.json` | Grunfeld (10 firms × 20 years) | ✅ EXIT 0 |

**关键设计**:
- JSON 输出 `digits=17` 保证 double 全精度 (R jsonlite 默认 digits=4 严重精度丢失)
- 硬编码公开数据集 (Spector-Mazzeo Greene 8ed Table 17.1, warpbreaks) 避免包依赖
- 手动实现 HAC 协方差矩阵 (Bartlett 核, Newey-West 带宽), 不依赖 `sandwich::kweights`
- CUE 目标函数检查 ḡ≈0 返回 J=0 (完美拟合处理)

### 5.2 Python 跨语言对照脚本 (4 个)

| 脚本 | Python 库 | 对照模块 | 状态 |
|------|-----------|---------|------|
| `cross_validate_ols.py` | statsmodels `OLS` | OLS + HC0-3 | ✅ EXIT 0 |
| `cross_validate_mle.py` | statsmodels `Logit`/`Poisson` | MLE | ✅ EXIT 0 |
| `cross_validate_gmm.py` | linearmodels `IVGMM`/`IVGMMCUE` | GMM + CUE | ✅ EXIT 0 |
| `cross_validate_bootstrap.py` | arch `bootstrap` | Bootstrap (Paired/Block/BCa) | ✅ EXIT 0 |

**修复项**:
- `cross_validate_bootstrap.py`: arch 库 `random_state`→`seed` 参数; `stat_func(*args)` 签名 (新版 arch 传参方式变化)

### 5.3 verify_econometrics.R 排幻觉点验证脚本

| 项 | 值 |
|----|-----|
| 脚本 | `verify_econometrics.R` |
| 验证项 | E1-E12 (12 项排幻觉点) |
| 结果 | 10 PASS + 2 VERIFIED |
| 依赖 | sandwich, lmtest, plm, gmm (已安装) |
| 状态 | ✅ 全部通过 |

---

## 6. 已知限制与遗留项

### 6.1 已知限制 (推迟到 v1.6+)

| 限制 | 影响 | 推迟原因 |
|------|------|---------|
| 非平衡面板 | Arellano-Bond 仅支持平衡面板 | 非平衡面板 GMM 工具变量构造复杂, v1.6+ 评估 |
| 前定变量 (predetermined) | Arellano-Bond 仅支持严格外生 | 前定变量需不同工具变量集, v1.6+ |
| Wild Cluster Bootstrap | 小聚类数 (G<20) 场景 | Cameron-Gelbach-Miller 2011 推荐, v1.6+ |
| Bootstrap OpenMP 并行 | 大规模 B 时性能 | v1.5 聚焦正确性, Phase 7 性能优化 |

### 6.2 遗留项

| 项 | 优先级 | 状态 | 备注 |
|----|--------|------|------|
| ~~A/B 站跨平台验证 M4~~ | ~~P1~~ | ✅ 已完成 | 1767/1767 全部通过 (commit b278151, 2026-08-05 23:26 CST) |
| PHASE6 最终验收报告 | P1 | ⏳ 进行中 | A/B 站验证已完成, 可汇总 M1-M4 出最终报告 |

---

## 7. 交付物清单

### 7.1 代码交付 (M4)

| 类型 | 数量 | 位置 |
|------|------|------|
| Bootstrap 头文件 | 5 | `include/cpphub/econometrics/resampling/` |
| M4 测试文件 | 6 | `tests/unit/econometrics/` |
| R 基准生成脚本 | 7 | `tests/fixtures/econometrics/` |
| R 基准 JSON | 7 | `tests/fixtures/econometrics/` |
| R 排幻觉点验证脚本 | 1 | `tests/fixtures/econometrics/verify_econometrics.R` |
| Python 跨语言对照脚本 | 4 | `tests/validation/python/` |

### 7.2 v1.5 累计交付 (M1+M2+M3+M4)

| 类型 | 数量 |
|------|------|
| 计量头文件 | 22 |
| 测试文件 | 25 |
| 测试用例 | 399 (spec 要求 185, +115.7%) |
| R 脚本 | 8 (7 基准 + 1 验证) |
| Python 脚本 | 4 |

---

## 8. 验收结论

### 8.1 M4 验收标准核查

| 检查项 | 通过标准 | 结果 |
|--------|----------|------|
| 配对 Bootstrap | vs Efron-Tibshirani 1e-3 | ✅ 通过 (Python arch 对照) |
| Wild Bootstrap | vs R multiwayvcov 1e-3 | ✅ 通过 (3 种权重分布) |
| Block Bootstrap | vs R boot::tsboot 1e-3 | ✅ 通过 (3 种块类型 + Politis-White) |
| Cluster Bootstrap | vs R multiwayvcov 1e-3 | ✅ 通过 (面板数据 entity/time 维度) |
| AR(1)/AR(2) 检验 | Arellano-Bond 1991 | ✅ 通过 (M3 偏差项补齐) |
| CUE 数值优化 | Nelder-Mead + Ŝ奇异处理 | ✅ 通过 (M3 偏差项补齐) |
| 端到端集成 | 17 个集成测试 | ✅ 全部通过 |
| 排幻觉点验证 | E1-E12 全部 | ✅ 10 PASS + 2 VERIFIED |
| R 基准生成 | 7 个脚本 EXIT 0 | ✅ 全部通过 |
| Python 跨语言对照 | 4 个脚本 EXIT 0 | ✅ 全部通过 |
| 测试数 | 35+ 新增 | ✅ 71 个 (M4) + 3 个 (M3 补齐) = 74 个 |

### 8.2 验收结论

**M4 验收通过**。

v1.5 经典参数计量模块 (M1+M2+M3+M4) 全部实施完成:
- 399 个测试用例 (spec 要求 185, 超出 115.7%)
- 12 项排幻觉点全部验证 (10 PASS + 2 VERIFIED)
- 7 个 R 基准 + 4 个 Python 跨语言对照全部通过
- 三平台全量测试通过 (主控站 MSVC + A/B 站 GCC, 1767/1767)

**跨平台修复**: `b278151` 解决了 `std::normal_distribution` 在 MSVC (Box-Muller) 与 libstdc++ (Marsaglia polar) 实现差异导致的测试数据不一致问题, 用项目自研 `Philox4x64` + `box_muller` 替换, 保证跨平台 RNG 完全一致。

---

**验收状态**: ✅ M4 验收通过 (三平台)
**下一步**: PHASE6 最终验收报告 → v1.6 Research OS 对接
