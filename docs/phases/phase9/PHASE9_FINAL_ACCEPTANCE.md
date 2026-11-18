# Phase 9 v1.9 最终验收报告 — KDE/KNN + KSG/TE 信息论基础设施

> **文档类型**: 最终验收报告 (Final Acceptance Report)
> **审计对象**: Phase 9 v1.9 完整 scope (M1 1D KDE + 带宽类型 + M2 kNN 密度/KSG MI/CMI + M3 迁移熵 TE + 集成; D9 决策边界全部遵守)
> **关联文档**:
>   - [PHASE9_SPEC.md](./PHASE9_SPEC.md) v1.1 (执行规格, §7 G1-G8 门禁全勾; §2 D9-1..D9-10 决策)
>   - [PHASE9_DOMAIN_RESEARCH.md](../../research/PHASE9_DOMAIN_RESEARCH.md) v1.0 (阶段 0 调研: A 67 / B 4 / C 7 / H 5, R1-R4 清零)
>   - [DEVELOPMENT_LOG.md](../../DEVELOPMENT_LOG.md) Phase 9 实施记录 (2026-11-11 实施 + 2026-11-14 CI 修复 + 2026-11-16 bundle 验证)
>   - [PHASE8_FINAL_ACCEPTANCE.md](../phase8/PHASE8_FINAL_ACCEPTANCE.md) (前置版本验收; D-2 裁决 (b) 载体)
> **验收日期**: 2026-11-16
> **验收员**: Scott (self-review, solo developer 工作流)
> **验收标准**: 三平台全量回归 + scipy `gaussian_kde` / sklearn `_compute_mi_cc` 黄金值冻结 + Gaussian MI 20-seed 统计锚 + TE=GC/2 ±30% + D9 决策零偏离 + v1.8 零破坏

---

## 1. 验收基本信息

| 项 | 值 |
|----|-----|
| 项目名称 | Cpp_Hub |
| 审计阶段 | Phase 9 v1.9 (非参数密度 + 信息论: 1D Gaussian KDE / kNN 密度 / KSG Algorithm 1 MI+CMI / 迁移熵) |
| 前置版本 | v1.8.0 (2543 主控 / 2525×2 A/B 站, tag 打于 767fa5e) |
| 后置版本 | v1.9 (**2671/2671** 主控站 MSVC / **2609/2609** ×2 A/B 站 GCC, 三平台全绿; CI run 43709054190 十 job 全绿) |
| Phase 9 新增测试 | **41 个** (test_kde 11 + test_knn_density 6 + test_ksg 10 + test_transfer_entropy 12 + test_integration_phase9 2) |
| 实现 commit | `5476d22` (feat, 17 文件 +3074 行) + CI 修复链 `ca6da32`/`1b92321`/`04af6c8` + 验收 `8c25bd5` (GCC 补录) |
| 验收方法 | TDD + 双 probe 前置探明 (probe_sklearn_kde.py 8/8 / probe_ksg.py 9 PASS + 1 SKIP) + scipy/sklearn 黄金值冻结 (gen_phase9_baselines.py) + bundle 中继三平台验证 |
| 里程碑 | 域调研 R 清零 (11-04) → probe 双 8/9 PASS (11-08) → spec 冻结 (11-10) → M1/M2/M3 三步 (11-11) → G1-G8 门禁 (11-11) → CI 修复两轮 (11-13/14) → A/B 站 bundle (11-16) → 验收+tag (11-16) |

> 主控与 GCC 站差额 62 = 配置差额 (主控 C ABI ON; 站与 CI build-test 口径 C ABI OFF) + 平台专属用例 (7B/7C 期已知口径), 非功能差异。

---

## 2. 全量回归测试结果

### 2.1 三平台矩阵 (2026-11-11 / 2026-11-16 实测)

| 平台 | 编译器 | 测试总数 | 通过数 | 失败数 | 耗时 | 验证方式 | 状态 |
|------|--------|---------|--------|--------|------|---------|------|
| 主控站 (Win10) | MSVC 19.50 (VS2026) | 2671 | 2671 | 0 | 全量 ctest (-j) | 本地 Release 全量回归 (C ABI ON) | ✅ |
| A 工作站 (`scott-lau-NEX.local`) | GCC 13.3 (Ubuntu 24.04) | 2609 | 2609 | 0 | 构建 63s (-j16) + ctest 67.8s (-j8) | 全量 bundle (9.5MB) fresh clone `~/Cpp_Hub_p9` + 构建 + ctest | ✅ |
| B 工作站 (`scott-lau-GTR-Pro.local`) | GCC 13.3 (Ubuntu 24.04) | 2609 | 2609 | 0 | 构建 62s (-j16) + ctest 61.1s (-j8) | 同法 (mDNS 通道) | ✅ |

### 2.2 GitHub Actions CI (run 43709054190, `04af6c8`)

| Job | 结果 |
|-----|------|
| Build & Test (windows-2022, MSVC) | ✅ success |
| Build & Test (ubuntu-latest, GCC) | ✅ success (2609/2609) |
| C ABI (windows-2022) | ✅ success |
| C ABI (ubuntu-latest) | ✅ success |
| Coverage (ubuntu-latest, GCC) | ✅ success (perf 用例插桩构建 GTEST_SKIP, §2.4) |
| Debug compile guard (MSVC) | ✅ success |
| Python (windows/ubuntu × py3.11/3.12) | ✅ success ×4 |

> CI 修复链: run 43465160758 (9/10, Coverage 插桩构建下 KSG 1.3s / TE 98.9s 超预算) → `1b92321` (perf 用例 `--coverage` 构建检测 + GTEST_SKIP) → run 43479814094 (9/10, 暴露既有 `CpuFeatures.ThreadSafeConcurrentAccess` 的 `std::vector<bool>` 位压缩数据竞争) → `04af6c8` (改 `std::array<bool,8>`, 本地 30/30) → run 43709054190 (10/10)。两次失败均非 Phase 9 算法缺陷, 修复全部落地并记录 (devlog 2026-11-14)。

### 2.3 跨平台验证执行细节

- A/B 站 bundle 中继发现站侧旧基线 (`~/Cpp_Hub_7C` @ `351ac98`) **不在主控当前谱系** (同内容提交 = `767fa5e`, 主控谱系于 v1.8 验收后重写, 旧 SHA 失效) → 增量 ff 不可行, 改**全量 bundle (9.5MB) + 站侧全新克隆 `~/Cpp_Hub_p9`**; eigen 子模块簿记损坏经离线复用 7C 工作树 (同 pin `3147391d`) 绕开, 构建不受影响 (CMake 仅 include-path)。
- 站侧配置对齐 CI GCC 口径: Release + `-ffp-contract=off`; gtest 经 `FETCHCONTENT_SOURCE_DIR_GOOGLETEST=~/gtest-src` 离线覆盖, 零外网依赖。
- Phase 9 新增 41 用例 GCC 双站全数运行 (perf 用例在 Release 真实执行), 三平台零数值偏差。

---

## 3. 实施范围汇总

### 3.1 新增头文件 (5 个, `cpphub::v1::timeseries::{density,information}` 命名空间, ADR-017 载体, 纯标量 header-only)

| 文件 | 行数 | 用途 |
|------|------|------|
| `density/bandwidth.hpp` | 117 | 带宽显式类型 (BandwidthRule/Bandwidth) + 四规则解析 (kScott 1.059·A·n^-0.2 / kSilverman 0.9 / kNormalReference 冻结常数 / kUser Absolute/Factor 语义) + `A = min(sd_ddof1, IQR/1.349)` |
| `density/kde.hpp` | 147 | KernelDensity1D (Gaussian 核 pdf/log_pdf max-shift log-sum-exp) + `bandwidth_sensitivity` 网格; UCRT 标量 exp 硬下限已实测, inv_h 乘法优化 |
| `density/knn_density.hpp` | 100 | Loftsgaarden-Quesenberry `k/(n·2r_k)` + `lower_bound` 左右归并取第 k 小 (O(log n + k)) + k 边界校验 |
| `information/ksg.hpp` | 263 | KSG Algorithm 1 MI (`ψ(k)+ψ(N)−⟨ψ(nx+1)+ψ(ny+1)⟩`) + Frenzel-Pompe CMI; Chebyshev + `nextafter(ε,0)` 严格计数 (sklearn `_compute_mi_cc` 逐位对齐); 自定义 digamma (max err 1.2e-15 vs scipy); z-score + seeded jitter |
| `information/transfer_entropy.hpp` | 252 | TE = CMI(Y_t; X_past \| Y_past) + effective TE (随机置换 source) + circular block shuffle p 值 (默认块长 ceil(√T)) + 任意 k/l/delay 嵌入; 连续坐标布局 `[X_past, Y_t, Y_past]` 单次扫描四 CMI (40→20ms/评估) |

### 3.2 fixtures / probe / gen 脚本 (3 个)

| 脚本 | 用途 |
|------|------|
| `probe_sklearn_kde.py` (199 行) | M1 前置: scipy `gaussian_kde` 八项语义探明 (8/8 PASS, 实测 7.5e-16) |
| `probe_ksg.py` (238 行) | M2 前置: sklearn `_compute_mi_cc` K0-K9 十项 (9 PASS / 1 SKIP) |
| `gen_phase9_baselines.py` (80 行) | 黄金值冻结生成器: 带宽五值 + KSG MI k=3/5 + GC + TE 九锚, 复跑逐位一致 |

### 3.3 新增测试 (5 套件 41 用例)

黄金值内嵌于测试源 (自包含, 7C 模式), 不依赖运行时外部文件; 性能用例预算口径 = MSVC Release (D9-10), 覆盖率插桩构建 (`--coverage -O0`) 经 `CPPHUB_COVERAGE_BUILD` 宏 GTEST_SKIP。

---

## 4. 数值基准验证 (SPEC §1.3 全表实测)

| 基准 | SPEC 容差 | 实测 | 状态 |
|------|----------|------|------|
| KDE pdf vs scipy `gaussian_kde` (同核同 h) | 1e-14 | 黄金值逐位一致 (probe 实测 7.5e-16) | ✅ |
| 带宽规则 vs 公式 (Scott/Silverman/normal_reference) | 1e-15 | 逐位一致 (kNormalReference C = 1.059223840025653 gen 阶段冻结, Q12 关闭) | ✅ |
| kNN 密度 vs 暴力参考 (同算法) | 1e-15 | 逐位一致 + `k=n−1` 边界 | ✅ |
| KSG MI vs sklearn `_compute_mi_cc` (40 点 k=3/5) | 1e-12 | 逐位一致 (含 nextafter 严格计数与"含自身再减 1"等价性实证 K0) | ✅ |
| KSG MI vs Gaussian 解析 (n=2000, k=3, **20 seeds**) | 均值偏差 ≤0.02 nat 判红 | 达标 (审计 A3: 0.01 记录线不判红) | ✅ |
| TE vs GC/2 (Gaussian VAR φ=0.5, a=0.3) | ±30% | 5 seeds 单测 + 3 seeds 集成全达标 | ✅ |
| TE 标准化后缩放不变 (x×0.1/1/10) | 1e-12 | 达标 | ✅ |
| CMI 条件独立均值 | \|CMI\| ≤ 0.05 | 达标 | ✅ |
| 零假设校准 5a (20 条 i.i.d. p>0.05 比例 ≥60%) | spec §6.3 | 达标 (固定 seed 确定性) | ✅ |
| 零假设校准 5b (耦合 VAR n=4000, p<0.1) | spec §6.3 | 达标 (功效判据冻结 n, devlog 教训 ③) | ✅ |
| 黄金值复跑 (gen_phase9_baselines.py) | 逐位一致 | kScott=0.5689708405001331 / kMiGoldenK3=0.3176876772600816 / kGc=0.11332868530700327 等九锚逐位 | ✅ |

---

## 5. 性能预算验证 (SPEC D9-10, 单线程)

| 用例 | 预算 | MSVC Release (主控) | GCC 13.3 (A/B 站) | 状态 |
|------|------|--------------------|-------------------|------|
| KDE N=10k/M=1k | <150 ms (修订, 经用户批准 50→150) | 52-58 ms | 0.05 s | ✅ |
| KSG N=5000/k=3 | <0.5 s | 35 ms | 0.03 s | ✅ |
| TE N=2000 + ETE 100 + 块置换 100 (~200 次 CMI) | <10 s | 4.1 s | 4.62-5.11 s | ✅ |

> KDE 修订链 (2026-11-11 G 门禁, 经用户批准): 实测 52-58ms, UCRT 标量 exp ~5.5ns/次为硬下限, 暴力算法已最优 (exp2 变换实测慢 8× 禁用); 150ms 留 CI runner 余量。覆盖率 CI 插桩构建 (Debug + `--coverage -O0`) 预算口径不适用 → `CPPHUB_COVERAGE_BUILD` + GTEST_SKIP (run 43465160758 实证 KSG 1.3s/TE 98.9s), Release 三平台仍真实断言。

---

## 6. D9 决策对齐审计 (10 项零偏离, SPEC §2)

| # | 决策 | 实施对照 | 状态 |
|---|------|---------|------|
| D9-1 | 仅 1D Gaussian 核 | 无 Epanechnikov/tophat/kernel_type 痕迹 (grep 实证) | ✅ |
| D9-2 | 带宽显式类型 + 三规则 + kUser 两语义 | bandwidth.hpp 四枚举; SJ/ISJ 无实现 | ✅ |
| D9-3 | kNN 暴力 Chebyshev O(N²) | knn_density.hpp 无 k-d tree; 复杂度文档声明 | ✅ |
| D9-4 | KSG 仅 Algorithm 1 | 无 Alg2 痕迹 (grep 实证, H5 维持 P2) | ✅ |
| D9-5 | TE = KSG-CMI + 嵌入; nats 规范 + bits 视图; ETE + 块置换 | transfer_entropy.hpp 全落地 (te_bits 字段) | ✅ |
| D9-6 | KSG/TE 默认标准化 | `standardize=true` 默认; 关闭须显式且 raw 敏感性有用例记录 | ✅ |
| D9-7 | 不做多元 KDE / 多元条件 TE / 离散 TE / Rényi TE | API 面无入口 | ✅ |
| D9-8 | 不做边界修正 / 权重 KDE | reflection = 未实现 throw (显式); 无 weight 参数 | ✅ |
| D9-9 | 不新增 C ABI / Python 绑定 | src/c_api 与 python/ 零触碰 (grep 实证) | ✅ |
| D9-10 | 性能预算 + 口径 (Release; 插桩 skip) | §5 全表 + `CPPHUB_COVERAGE_BUILD` 机制 | ✅ |

**不做清单核查**: 多核/多维 KDE、带宽矩阵、SJ/ISJ、边界修正、加权 KDE、k-d tree、KSG Alg2、多元/离散/Rényi TE、排列熵、滑窗 TE、C ABI、Python — **grep + API 面双重实证零实现痕迹**。

---

## 7. 零破坏与工程收尾 (v1.8 兼容)

- ✅ 既有代码唯一实质改动 = `tests/CMakeLists.txt` (+17 注册 +8 coverage 宏, 纯增量); 既有测试修复 = `test_cpu_features.cpp` (`vector<bool>`→`array<bool>`, CI 修复批, 数据竞争消除, 非 Phase 9 模块面)。
- ✅ 零 Eigen 依赖: density/information 五头仅标准库 + core/types (grep 实证); 与 garch/unit_root 同层挂 `cpphub_timeseries` INTERFACE。
- ✅ C ABI / Python 绑定零变更 (D9-9); `python/pyproject.toml` 1.9.0 为发布轴预置, 无代码面改动。
- ✅ `VERSION` 文件 (C++ API 轴, RISK-032 两轴模型) 保持 1.0.0 冻结; 发布轴以 tag v1.9.0 承载。
- ✅ 文档链同步: PHASE9_SPEC v1.1 (§7 G1-G8 全勾 + D9-10 两轮修订链) / README Status 表 41/41 / devlog 四条目。

---

## 8. 口径对齐与实施期发现

- **口径收紧** (测试→spec §6.3): Gaussian MI 3 seeds/0.05 → **20 seeds/0.02**; TE ±35% → **±30%**; 新增 p 值校准 5a/5b 两用例 (此前仅范围断言); integration MI 同步 20 seeds。全部固定 seed, 确定性可复现。
- **实施期发现四项** (devlog 钉死): ① UCRT `std::exp2` 慢 `std::exp` 8× (性能陷阱, 禁用); ② TE 嵌入连续坐标布局单次扫描四 CMI (40→20ms); ③ TE 耦合功效须冻结 n (n=2000 seed 502 p=0.157 → n=4000 达标); ④ sklearn `kneighbors()` X=None 排除自身语义与本库等价 (K0 逐位实证)。
- **测试计数偏差记录**: spec 预估「全库 ≥2700」作废 — 实际新增 41 例, 判据以「5 套件全绿 + 全量零失败 + 计数如实记录」为准 (同 KDE 预算修订模式, 用户批准)。

---

## 9. 遗留问题与后续

| # | 项 | 状态 | 处置 |
|---|----|------|------|
| 1 | H2 JIDT TE 运行时对照 (风险 ⑭ 关闭路径) | ✅ **已关闭 (2026-11-16, tag 后 verify 批)** | probe_jidt_te.py: delay=1 结构区 12 组 max 6.53e-4 nat; JIDT 连续 TE = nats (订正 A47); JidtAssemblyAnchor 锚固化; 详见 spec §10.1 H2 + devlog |
| 2 | H4 (Kirkley reduced TE) / Q11 (Politis-White 块长) | ✅ **已关闭 (2026-11-17, P2 批)** | Q11: PW 自动块长落地 (opt-in) + arch 逐位对照 + 持久零假设误拒绝 16.7%→8.3%; H4: 离散化实证证伪 ±10% 假设 (reduced TE 正偏→强负偏, T=2e4 仍 −12.6%), 不作替代采纳; 详见 spec §10.1 |
| 2b | H1 (ks hpi/hscv) / H5 (KSG Alg2) | H1 ✅ **已关闭 (2026-11-17, P2 批, 结论: 非同一选择器)**; H5 OPEN, P2 | H1: ks 1.15.3 源码 — hpi 委托 KernSmooth::dpik (Wand-Jones plug-in), hscv 为 SCV 优化, 均非 SJ 1991 ste (见 spec §10.1); R 4.6.1 运行时复核: `ks::hpi` 与 `KernSmooth::dpik` 三组夹具**逐位相同**, 与 `bw.SJ` 显著不同; H5 维持冻结 (Alg1 已双重锚定) |
| 2c | P2-4 k-d tree 加速 (风险 ③) | ✅ **已关闭 (2026-11-17, P2 第 2 批)** | `density/kdtree.hpp` (Chebyshev, 列子集) 接入 KSG MI/CMI + TE 嵌入; 与独立暴力 oracle 逐位一致 (test_kdtree 6 例 + test_ksg MI/CMI 交叉验证), sklearn 1e-12 / JIDT 锚零漂移; 实测 KSG N=5000 37→7 ms, TE N=2000 3871→1606 ms; 详见 spec §10.2 P2-4 |
| 2d | P2-2 SJ/ISJ 带宽选择器 | ✅ **已关闭 (2026-11-17, P2 第 3 批)** | `BandwidthRule::kSj`/`kIsj` 落地 (R bw.SJ ste / KDEpy ISJ 同算法); probe_sj_isj.py 三组夹具锚, C++ 对齐 1e-10/1e-8; test_kde 11→17 例; 详见 spec §3.1/§10.2 |
| 2e | P2-1 反射边界修正 / P2-5 权重 KDE | ✅ **已关闭 (2026-11-17, P2 第 4 批)** | `KdeConfig::reflection` 启用 0 处反射 (样本严格正, ∫_0^∞=1 恒等式); `KernelDensity1D(sample, weights, cfg)` 权重仅进混合 (带宽按未加权样本, statsmodels/sklearn 同语义), `effective_sample_size()` 诊断; probe_weighted_kde.py W1-W5 + test_kde 17→28 例; 详见 spec §3.2 |
| 2f | P2-3 多变量 KDE / 条件 TE | ✅ **已关闭 (2026-11-18, P2 第 5 批)** | `KernelDensityND` (scipy 全协方差语义, Cholesky, 锚 1e-12) + `conditional_transfer_entropy` (Z_past 嵌入, 暴力 CMI 交叉验证 1e-12); probe_multivariate_kde.py W1-W3 + test_kde_nd 10 例 + 条件 TE 5 例; 详见 spec §3.4/§5.1 |
| 3 | H3 (sklearn 'scott' 裸因子) / Q12 (normal_reference C 出处) | ✅ 已关闭 | H3 不影响实现 (文档记录); Q12 gen 阶段冻结 C=1.059223840025653 |
| 4 | P2-1..P2-8 (P2-6 Alg2 / P2-7 多核 / P2-8 滑窗; **P2-1/P2-2/P2-3/P2-4/P2-5 已关闭, 见行 2c/2d/2e/2f**) | 登记 | spec §10.2 对应 D9 决策 |
| 5 | P1 工程项: RISK-033 `URL_HASH` / RISK-034 json_min 上限 / EXTERNAL_REVIEW 表回填 | ✅ **已关闭 (2026-11-16, P1 清零批)** | URL_HASH SHA256 钉死 + CI 双重守卫; json_min 入口 8 MiB + 节点预算 1e6 (既有实现复核, test_json_min 19/19 + test_c_api_tools 11/11); 审计主表 14/14 回填 (RISK-034 安全口径复核保留为独立安全评审动作) |
| 6 | tag v1.9.0 | ✅ 本验收后执行 | 打于验收 commit 并推送 origin |

## 10. 验收结论

**Phase 9 v1.9 验收通过。** 三平台 2671 (MSVC) / 2609×2 (GCC 13.3) 全绿, CI run 43709054190 十 job 全绿; SPEC §1.3 数值基准全表达成 (黄金值逐位 + 20-seed 统计锚 + p 值双向校准); D9-1..D9-10 十项决策零偏离, 不做清单 grep 双重实证零痕迹; 性能三用例双编译器达标 (预算修订链完整留痕); v1.8 零破坏 (唯一既有代码改动为测试注册与一个数据竞争修复)。信息论基础设施 (KDE/KNN/KSG/TE) 就位, 为 Factor Fingerprint 因果失效诊断与 v2.0 信息论事前度量铺路。
