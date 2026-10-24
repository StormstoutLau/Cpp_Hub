# P2 / 遗留 RISK 调研: RISK-029 / 030 / 031 / 035 / 036

> **调研日期**: 2026-10-24 · **方法**: 5 路并行检索（社区 issue / 标准草案 / 教材与论文 / 参考实现源码）+ 主 agent 源码复核
> **上游**: [EXTERNAL_REVIEW_2026-10-20.md](../../docs/audit/EXTERNAL_REVIEW_2026-10-20.md) §二 主表
> **证据分级**（沿用 [ASSERTION_EVIDENCE_FRAMEWORK.md](../../docs/ASSERTION_EVIDENCE_FRAMEWORK.md)）: A = 可点击一手来源逐字核对；B = 本地实测/推理链；C = 判断类（仅 rationale）

## 0. 结论速览

| RISK | 复核结论 | 推荐方案 | 变更面 | 建议批次 |
|---|---|---|---|---|
| 029 架构硬编码 | 成立（Release 条件化但未按架构条件化） | `option(CPPHUB_X86_V3 OFF)` + 平台/编译器守卫 + target 级选项 | 根 `CMakeLists.txt` 1 文件 | 第 1 批（独立） |
| 030 `mutable` 计数器竞态 | 成立（const 方法写共享态） | 迭代次数并入返回值（`PricingResult` 重载），旧签名保留包装 | `pde_engine.hpp` + ADR | 第 2 批（PDE 同批） |
| 031 边界字符串分派 / Neumann 失效 | 成立（数字期权静默用 put 边界；障碍路径忽略 Neumann） | `PayOff::kind()` 枚举 + 远场策略；Unknown 抛错；Neumann 显式处理 | `payoff.hpp` / `vanilla.hpp` / `pde_engine.hpp` | 第 2 批（PDE 同批） |
| 035 `chow_denning_joint` 字段语义 | 成立，但"无联合统计量"描述不精确——已有逐 horizon `Z2`，缺 `max\|Z2\|` 标量 | `VarianceRatioResult` 增 `chow_denning_statistic`，registry 分别填 statistic/p_value | 1 头文件 + 1 registry 行 + golden | 第 3 批（数值契约） |
| 036 χ² 上尾精度 | 成立；且 `regularized_upper_gamma` 已存在但为 `1-P` 且无调用者 | 原地重写为真上尾（连分式直返 Q）+ 新增 `chi2_sf` | `core/math.hpp` + `backtesting.hpp` 3 处 | 第 3 批（数值契约） |

---

## 1. RISK-029 — 架构指令集硬编码

### 1.1 事实复核（A）

`CMakeLists.txt:16-23`（Release 块 `:17-23`）在 `CMAKE_BUILD_TYPE == Release` 时全局追加 `/arch:AVX2`（MSVC）或 `-march=x86-64-v3`（GCC/Clang）。**未按 `CMAKE_SYSTEM_PROCESSOR` 条件化**：aarch64 上 GCC 直接报 unknown `-march`；产物在不支持 AVX2 的 CPU 上可 SIGILL。`python/CMakeLists.txt` 不含这些选项，wheel 不受影响；受影响面 = 根工程 Release 构建。

### 1.2 社区案例（A）

| 案例 | 结果 | 来源 |
|---|---|---|
| TensorFlow 1.6 manylinux1 轮子 | 无 AVX 的 Core i3 M330 上 `import` 即 SIGILL | [tensorflow#17411](https://github.com/tensorflow/tensorflow/issues/17411) |
| LanceDB | 默认 wheel 按 x86-64-haswell（AVX2+FMA+F16C）编译，旧 CPU import 即崩；另发 compat 包（v2 基线 + 运行时派发） | [pypi.org/project/lancedb](https://pypi.org/project/lancedb/) |
| manylinux | AlmaLinux 9 基线为 x86-64-v2，wheel tag 无法表达 v2，产物在旧 CPU SIGILL；issue 仍 open | [pypa/manylinux#1725](https://github.com/pypa/manylinux/issues/1725) |
| cryptography 44.0.0 | QEMU x86_64 上 import SIGILL | [pyca/cryptography#12069](https://github.com/pyca/cryptography/issues/12069) |
| NumPy | 构建默认 `cpu-baseline=min`、`cpu-dispatch=max` | [numpy SIMD build options](https://numpy.org/doc/stable/reference/simd/build-options.html) |

### 1.3 标准/工具链事实（A）

- x86-64 psABI 定义 v1–v4 微架构级别：<https://gitlab.com/x86-psABIs/x86-64-ABI>；GCC `-march=x86-64-v3` 即取该级：<https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html>
- MSVC `/arch` 只有 SSE2/SSE4.2/AVX/AVX2/AVX512/AVX10，**没有 v2/v3 级别**：<https://learn.microsoft.com/en-us/cpp/build/reference/arch-x64>
- GCC 函数多版本化 `target_clones` 依赖 ELF ifunc（binutils ≥2.20.1 / glibc ≥2.11.1），Windows/macOS 不可用：<https://gcc.gnu.org/onlinedocs/gcc/Function-Multiversioning.html>
- manylinux 对 RHEL9+ 的 v2 基线无法被 auditwheel 检测：<https://github.com/pypa/manylinux#caveat>；PEP 825（Provisional）方向是多基线变体而非抬高默认：<https://peps.python.org/pep-0825/>

### 1.4 方案对比（C）

| 方案 | pros | cons |
|---|---|---|
| A. 基线编译 + `target_clones` | 单二进制、按函数多版本 | 仅 ELF；MSVC 无；宏展开后样板多 |
| B. 多 TU + 手写运行时分派 | 全平台可用（xsimd 模式） | 样板多；跨 TU 边界只能标量 |
| C. Highway 式静态/动态派发 | 自动、覆盖广 | 引入依赖；静态派发遇旧 CPU 必崩 |
| **D. option 默认 OFF（推荐）** | 零依赖、零风险、消费者可 opt-in | 默认 Release 性能回退到 baseline（对量化库可接受） |

### 1.5 推荐（B）

```cmake
option(CPPHUB_X86_V3 "Enable x86-64-v3 / AVX2 codegen (requires AVX2 CPU)" OFF)
if(CPPHUB_X86_V3 AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)$")
  include(CheckCXXCompilerFlag)
  # 分别探测后 add_compile_options 或 target 级追加；失败即 FATAL_ERROR
endif()
```

要点：不覆盖 Release 默认 `-O3 /O2`；不改全局 `CMAKE_CXX_FLAGS`；wheel 路径不动；`cpu_features.hpp` 保持运行时唯一真源。

### 1.6 陷阱（A/B）

- `CMAKE_SYSTEM_PROCESSOR` 在 Windows 与交叉编译下不可靠，交叉场景必须走 toolchain 文件。
- `CheckCXXCompilerFlag` 只证"能编译"，不证"目标 CPU 支持"。
- 全局 `-march` 会让静态初始化器在 `main` 之前执行 AVX2 指令，运行时派发救不了。
- MSVC 无 `target_clones`；不同 `/arch` 的 TU 不可跨 TU 内联（LTCG 亦然）。

---

## 2. RISK-030 — `mutable` 计数器数据竞争

### 2.1 事实复核（A）

`pde_engine.hpp:259` `mutable Size last_total_iter_count_`；写点 `:91`、`:125`（均在 `const` 方法 `price_american` 内）；读点 `:255`。两个线程对同一 `const PDEEngine` 并发调用即 data race → UB。

### 2.2 标准与学术（A）

- mutable 的正文在 **[dcl.stc]/8–9**（公开草案中**不存在** `[class.mutable]` 子节）：<https://timsong-cpp.github.io/cppwp/n4618/dcl.stc>
- **[dcl.type.cv]/4**：const 对象上修改 mutable 成员 well-formed——const 是编译期访问控制，不是并发承诺：<https://timsong-cpp.github.io/cppwp/n4861/dcl.type.cv>
- **[intro.races]/17**：potentially concurrent conflicting actions + 至少一个非原子 + 无 happens-before → data race → UB；`/3` 注：relaxed 原子不构成 data race，但**不是同步操作**：<https://eel.is/c++draft/intro.races>
- **[res.on.data.races]/2–3** 与 **N2669**：标准库不得经 const 参数修改他线程可见对象——"const≈线程安全"惯例的规范出处，但**只约束标准库**：<https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2008/n2669.htm>
- 引文更正：Boehm & Adve *Foundations of the C++ Concurrency Memory Model* 发表于 **PLDI 2008**（非 POPL 2008）；Boehm *Threads Cannot Be Implemented as a Library* 为 PLDI 2005, pp.261–268。

### 2.3 社区案例（A）

| 案例 | 结果 |
|---|---|
| [libstdc++ Bug 77704](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=77704) | `std::ctype<char>::narrow() const` 读写 mutable 缓存，TSAN 报 data race；长期 NEW |
| [QuantLib LazyObject 讨论](http://quantlib.414.s1.nabble.com/Multithreading-and-LazyObject-td15310.html) | `calculate() const` 写 `calculated_`，OpenMP 下可见"已计算但未完成" |
| [JDK-8166842 讨论](https://mail.openjdk.org/archives/list/core-libs-dev@openjdk.org/message/VYHWZJWC6ZIX2KFDJIBECNKOGNPONNZH/) | `String.hashCode()` 的 hash 缓存"良性竞争"假设被否定 |
| [Abseil: Revisiting Regular Types](https://abseil.io/blog/20180531-regular-types) | mutable 成员使类型 thread-unsafe；"不存在安全的 data race" |
| [Eigen 多线程文档](https://libeigen.gitlab.io/eigen/docs-nightly/TopicMultiThreading.html) | 显式列出非线程安全入口，而非隐式承诺 const 并发安全 |

### 2.4 方案对比（B/C）

| 方案 | 结论 |
|---|---|
| ① 迭代次数并入返回值（`PricingResult{value, iterations}`） | **推荐**：消除共享态，语义"本次调用的迭代数"成立 |
| ② `std::atomic<Size>` relaxed | 只消 UB；并发调用仍互相覆盖，`last_total_iterations()` 归属不确定 |
| ③ `thread_local` | 是"每线程"而非"每实例每次调用"，线程池会串值 |
| ④ 删除访问器 | 最干净但破坏既有调用方 |
| ⑤ 仅文档声明非线程安全 | 零改动，不消 UB，转嫁风险 |

### 2.5 推荐（B）

`price_american` 内部用局部 `Size iters` 累加 → 新增返回 `PricingResult` 的重载；旧 `Real price_american(...)` 保留为 inline 包装；`last_total_iterations()` 标 deprecated 并在下一主版本删除。走 ADR 记录弃用期。**不要只换 atomic 收工**（见 2.4②）。

### 2.6 性能说明（B）

累加发生在每个 time step 一次（≈`n_time` 次），不是 PSOR 内每次迭代；若保留 atomic 也只是一次 relaxed store（x86 上等价普通 `mov`）。公开基准（如 [arXiv:2010.09852](https://arxiv.org/abs/2010.09852)）针对 CAS/FAA，**无"每步累加"专项基准**，结论为量级判断。

---

## 3. RISK-031 — 边界字符串分派与 Neumann 失效

### 3.1 事实复核（A）

- `pde_engine.hpp:186-204`：障碍侧置 0；远场侧 `payoff.name() == "Call"` 才用 `bsm_call`，否则一律 `bsm_put`。
- 内置 payoff 名称 5 种（Call/Put/DigitalCall/DigitalPut/DoubleDigital），即**数字类全部落到 put 边界**。
- `price_barrier` 循环体**不调用** `finalize_boundary`；Neumann 仅在 `:285` 与 `:308-312`（非障碍路径）生效 → `BoundaryType::Neumann` 在障碍路径被静默忽略。
- 对照：`set_boundary` 的 else 分支（:294-304）已用 `payoff(S) * df`，说明"正确做法"在库内已有先例。

### 3.2 学术与工程资料（A）

- Duffy, *Finite Difference Methods in Financial Engineering* (Wiley), Ch.14 §14.3 p.154：障碍期权是 IBVP，解由终值 + 边界数据唯一确定；put 解析解只解 put 的问题（出版社目录核对，章节页未在线获取）。
- Wilmott, Howison & Dewynne (1995), *The Mathematics of Financial Derivatives*, Ch.12 "Barrier options"。
- Fusai, Sanfelici & Tagliani (2002), *Practical Problems in the Numerical Solution of PDE's in Finance*：数字/障碍 payoff 不连续 → 振荡、离散最大原理失效。<https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/Fusai_PracticalProblems.pdf>
- QuantLib `fdblackscholesbarrierengine.cpp`：障碍边 Dirichlet(rebate) + 截断域；不支持的 payoff 直接 `QL_REQUIRE` 抛错（fail fast）。<https://github.com/lballabio/QuantLib/blob/master/ql/pricingengines/barrier/fdblackscholesbarrierengine.cpp>
- QuantLib `boundarycondition.hpp`：抽象 BC 基类；`NeumannBC`/`DirichletBC` 自 1.42 起 deprecated。<https://github.com/lballabio/QuantLib/blob/master/ql/methods/finitedifferences/boundarycondition.hpp>
- QuEP 2：enum+switch 处理 BC 被判为"严重损害扩展性"，改为多态 BC，且边界策略与 payoff 解耦。<https://www.quantlib.org/quep/quep002.html>
- Google C++ Style（RTTI 节）：运行时类型查询常意味着设计缺陷；Core Guidelines I.4：接口精确强类型。<https://google.github.io/styleguide/cppguide.html#Run-Time_Type_Information__RTTI_>

### 3.3 数值证据（B，单参数组，未独立复现）

调研 agent 用 MSVC 14.50 临时探针（400×1000，Crank-Nicolson，S0=K=100，T=0.5，r=0.05，q=0.02，σ=0.25，默认 alpha=0.2）测得：

| 场景 | 现实现 | 宽域参考 | 偏差 |
|---|---|---|---|
| Up-Out 数字 put（H=130） | 5.198 | 0.487（欧式解析 0.489） | **+963%** |
| Down-Out 数字 call（H=85） | 0.310 | 0.45–0.47 | **−33%** |
| 同参数 vanilla | 跨 alpha 差 0.002% | — | 正常 |
| Neumann vs Dirichlet（障碍路径） | 完全相同（1e-10） | — | 证实被忽略 |

**证据局限**：单一参数组、单网格族；数字期权宽域值本身有 ±2–3% 振荡散射，0.487 / 0.45–0.47 为近似参考。实施时应以收敛测试（网格加倍 + 域扩大）作为正式锚，不复用该单点数值。

### 3.4 方案对比（C）

| 方案 | pros | cons |
|---|---|---|
| ① `PayOff::boundary_value(S,t)` 虚函数 | 类型安全、可扩展 | 破坏 ABI；把 r/q/σ 塞进 payoff，职责混乱 |
| ② 独立 `BoundaryPolicy` 参数（默认 Auto） | 不改 PayOff | Auto 仍要类型识别；签名膨胀 |
| ③ `PayOff::kind()` 枚举 + switch | 无字符串、`-Wswitch` 可穷举（QuantLib 先例）；改动小 | 仍是运行时；Custom 必须 throw |
| ④ 仅线性 payoff 用解析边界，其余 `payoff(S)·e^{-rτ}` | 数字期权即正确；改动最小 | 任意 payoff 的"远场常数"假设不成立 |
| ⑤ Neumann 在障碍路径 throw | 诚实、最小 | 丢失功能；更好是仅远场端外推 |

### 3.5 推荐（B）

- **③+④**：`PayOff` 增非纯虚 `kind()`（默认 `Unknown`）；Call/Put 走 BSM；Digital*/DoubleDigital 用 `payoff(S_far)·e^{-rτ}`；`Unknown` 抛 `std::invalid_argument`（自定义 payoff 可 override `kind()` 选用该边界）。约 15 行，Python/C API 无影响。
- **Neumann**：在障碍路径要么只对远场端实现线性外推（障碍端保持 0），要么显式 throw；**禁止静默忽略**。
- 新增测试：数字障碍收敛、Neumann≠Dirichlet、Unknown 抛错。

### 3.6 陷阱（A/B）

- 默认 `alpha=0.2` 远场仅约 `1σ√T`；加大域只会掩盖错误边界，不是修复。
- 复用 `set_boundary`/`finalize_boundary` 会同时写两端，抹掉障碍端 `V=0`。
- `payoff(S)·e^{-rτ}` 只对远场为常数的 payoff 成立；线性 payoff 会漏 `S` 项。
- `name()` 可被用户 override，自定义 payoff 取名 "Call" 即静默劫持边界——这正是改用 `kind()` 的理由。

---

## 4. RISK-035 — `chow_denning_joint` 的 `statistic` 语义

### 4.1 事实复核（A，含对上游描述的更正）

- `cpphub_tool_registry.hpp:543-547` 把 `chow_denning_p_value` 同时填入 `statistic` 与 `p_value`。
- 上游称"库内无联合统计量"**不精确**：`variance_ratio_test.hpp:49` 已有逐 horizon 的 `chow_denning_stats`（有符号 `Z2`），内部 `:224-236` 已计算 `cd_abs = max|Z2|` 但**未暴露**。缺的只是"把 `max|Z2|` 存进结果结构"。

### 4.2 文献与实现（A）

- Chow & Denning, "A simple multiple variance ratio test", *Journal of Econometrics* 58(3):385–401 (1993), DOI [10.1016/0304-4076(93)90051-6](https://doi.org/10.1016/0304-4076(93)90051-6)。联合统计量 `CD = max_i |Z2(k_i)|`。
- R [vrtest::Chow.Denning](https://rdrr.io/cran/vrtest/man/Chow.Denning.html) 返回 `CD1/CD2 = max|Z2|` 与 10/5/1% 临界值，**不返回 p 值**；源码 `per=0.5*(1-(1-α)^(1/m))`、`crit=qnorm(1-per)`，代数反演即 `p = 1 - [2Φ(|CD|)-1]^m`，`m = length(kvec)`（horizon 个数）。<https://rdrr.io/cran/vrtest/src/R/Chow.Denning.R>
- 该 p 值是 Šidák / SMM(m,∞) 近似；原文有限样本用 SMM(m,T)。各 `Z2(k)` 实际正相关 → p 偏大、检验保守。
- statsmodels 的 tsa API 无 VR/Chow-Denning；Python 侧 `arch.unitroot.VarianceRatio` 只有单 horizon。

### 4.3 契约惯例（A）

- R `htest`：`statistic` 与 `p.value` 并列（[chisq.test](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/chisq.test.html)）。
- [broom::tidy.htest](https://broom.tidymodels.org/reference/tidy.htest.html)：`statistic` 列 = 检验统计量，`p.value` 独立。
- [statcheck](https://www.stats.bris.ac.uk/R/web/packages/statcheck/readme/README.html) 依赖 `statistic` 重算 p——statistic 填 p 值会使复核失效。
- 未检索到"p 值误填 statistic"的公开事故（C 级：仅约定层面佐证）。

### 4.4 方案对比（C）

| 方案 | pros | cons |
|---|---|---|
| ① 增 `chow_denning_statistic`（max\|Z2\|） | 语义正确、零数值成本、C ABI JSON shape 不变；golden 仅 1 个数变化 | 改头文件 + golden 重生成 |
| ② `statistic=null` + note | 诚实、库不改 | 丢可算信息；skeleton `<num>`→`null`，numbers 少 1 个 |
| ③ 维持现状 | 零改动 | 语义持续错，F9 常驻，无法复核 |

### 4.5 推荐（B）

方案 ①：`VarianceRatioResult` 增 `chow_denning_statistic`（单 horizon 时 = `|z2|`；多 horizon = 已算的 `cd_abs`）；registry `:544` 用新字段填 `statistic`；golden 该行对应数字改为 max|Z2|（p 不变）；skeleton 不变 → Python 测试代码不动；关闭 `spec/mcp-server/CHECKLIST.md` §1.5 F9 与 `meta.known_issues`。

### 4.6 陷阱（A/B）

- CD 非负：输出 `max|Z2|`，不要输出 argmax 的有符号 `Z2`。
- `m` 是 horizon **个数**（去重后），不是 `max(k)`。
- `1 - pow(a,m)` 在 p 极小时抵消，建议 `-expm1(m * log1p(-2*Φ(-|CD|)))`。
- p 为 SMM(m,∞)/独立近似，正相关时保守，note 应写明。

---

## 5. RISK-036 — χ² 上尾精度

### 5.1 事实复核（A，含对处置候选的更正）

- `backtesting.hpp:78/107/131` 三处 `1 - chi2_cdf(...)`。
- `core/math.hpp:144-146` **已存在** `regularized_upper_gamma`，但实现就是 `1.0 - regularized_lower_gamma(...)`，且**全仓库无调用者**（grep）。所以"新增精确上尾"实为**原地重写**，接口不变。
- `regularized_lower_gamma` 的连分式分支（`x >= a+1`）在返回前先得到 `Q` 再 `return 1.0 - Q`（math.hpp:120-140）——精确 Q 可直接从该分支返回。

### 5.2 算法/文献（A）

- Numerical Recipes 3rd ed. §6.2 `gammp`/`gammq`：`x < a+1` 用级数，否则用连分式；`gammq` 在连分式支路**直接返回 Q**，不做减法。<https://www.aquila.infn.it/pierleoni/Laboratorio_Calcolatori2/ROUTINES/Chap6.pdf>
- DiDonato & Morris, *Computation of the incomplete gamma function ratios and their inverse*, ACM TOMS 12(4):377–393, 1986, [doi:10.1145/22721.23109](https://dl.acm.org/doi/10.1145/22721.23109)。
- Cephes：`chdtr(v,x) = igamc(v/2, x/2)`（即上尾）。<https://www.netlib.org/cephes/doubldoc.html>
- Boost.Math `gamma_q` 与 `gamma_p` 是独立入口；文档明言互补量计算会"lose accuracy, perhaps totally"。<https://www.boost.org/doc/libs/1_84_0/libs/math/doc/html/math_toolkit/stat_tut/overview/complements.html>
- `std::gammaincc` **不存在**（C++17/20/26 `<cmath>` 特殊函数清单无此函数）：<https://en.cppreference.com/w/cpp/numeric/special_functions>

### 5.3 社区案例（A）

- Goldberg 1991 §2.4 Cancellation：<https://docs.oracle.com/cd/E19422-01/819-3693/ncg_goldberg.html>
- Boost "Complements are supported too"：`1-0.999999940f` 只剩 1 bit。
- R `Chisquare` 手册：非中心卡方只算下尾，"upper tail suffers from cancellation and a warning will be given"。<https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Chisquare.html>
- scipy `sf`："also defined as 1-cdf, but sf is sometimes more accurate"；1.11.0 修复 `dgamma.sf` 尾部精度：<https://scipy.github.io/devdocs/release/1.11.0-notes.html>

### 5.4 精度数据（B，本地复算：k=1 用 `erfc` 闭式；k=5 用 Γ 递推；均与 Lentz 连分式交叉到 ~1e-14）

| k | LR | 真值 Q | 现实现 `1-P` | 相对误差 |
|---|---|---|---|---|
| 1 | 64 | 1.244e-15 | 1.2e-15 | ~3% |
| 1 | 70 | 5.93e-17 | 1.11e-16 | ≈100%（剩 1 个格点） |
| 1 | 74 | 7.8117e-18 | 0 | 100% |
| 5 | 74 | 1.50403e-14 | 1.5e-14 | 3.5e-3 |
| 5 | 200 | 2.84062e-41 | 0 | 100% |
| 5 | 1000 | 6.01008e-214 | 0 | 100% |

**对上游的描述更正**：k=1 的归零点约在 `LR≈70`（`Q < ε/2 = 1.11e-16`），不是 74；74 是"取整后的近似"。

### 5.5 方案对比（C）

| 方案 | pros | cons |
|---|---|---|
| ① 重写 `regularized_upper_gamma` + 新增 `chi2_sf` | 单一来源、全 k 精度到 ~1e-306、接口不变 | 触碰 core 头文件（header-only 全量重编译一次） |
| ② 只在 `backtesting.hpp` 局部实现 | 重建面小 | 复制连分式与阈值，`midas_diagnostics.hpp:285` 同病无法复用 |
| ③ `std::gammaincc` | — | 不存在，排除 |
| ④ 维持现状 + 文档声明 | 零风险 | p 不可读；k=2 已有闭式 `exp(-LR/2)`，说明缺 sf 是真实缺口 |

"全库重编译"代价真实但有限：`core/math.hpp` 被 40 个文件 include，INTERFACE 库本就全量重编；header-only inline 无 ABI/ODR 风险，`regularized_upper_gamma` 当前无调用者 → 现有 P 数值不变。

### 5.6 推荐（B）

`math.hpp` 重写 `regularized_upper_gamma`：`x >= a+1` 直接返回连分式 Q；`x < a+1` 返回 `1-P`（该区间 Q ≥ 0.083，`1-P` 相对误差 ≤3e-15）。新增 `chi2_sf(x,k) = Q(k/2, x/2)`；`backtesting.hpp` 三处改用 `chi2_sf`。后续可把 `midas_diagnostics.hpp:285`（现为 `1 - regularized_lower_gamma`）一并迁移。同步更新 `docs/USAGE_GUIDE.md:167`（当前把 upper 描述为 `1 - P`）。

**测试锚**：k=1 可用 `std::erfc(std::sqrt(LR/2))` 作为零依赖 oracle：`LR=74 → 7.811702936622268e-18`；k=5：`74 → 1.5040254551314682e-14`、`200 → 2.8406228986415315e-41`、`1000 → 6.0100776879208042e-214`（相对容差 1e-12）。外部复核：`scipy.stats.chi2.sf` 或 `pchisq(LR,k,lower.tail=FALSE)`。

### 5.7 陷阱（A/B）

- 分支阈值必须与 NR 一致用 `x >= a+1`；挪进小 Q 区重新引入抵消。
- Q 在 LR≈1500（k=1）underflow 为 0（scipy `gammaincc(0.5,1000)=0` 同理）；要更远尾部需 log-space，本次不做。
- `x <= 0` 时 sf 返回 1，保持"负 LR → p=1"现有语义。
- P/Q 必须严格互补（P 走级数、Q 走连分式，各自相对精确）。

---

## 6. 实施顺序与待裁定项

**建议批次**：
1. **RISK-029**（独立、CMake 单文件）——先做，风险最低，立即消除 aarch64/旧 CPU 构建失败。
2. **RISK-030 + 031**（同属 `pde_engine.hpp`，一次 ADR：`PricingResult` 重载 + `PayOff::kind()` + Neumann 策略）——合并减少 ADR 与重编译次数。
3. **RISK-035 + 036**（数值契约批）——035 改结构体 + registry + golden；036 改 `math.hpp` + `backtesting.hpp` + 测试锚。

**需 Scott 裁定的点**：
- RISK-029：`CPPHUB_X86_V3` 是否默认 OFF（推荐），或保留 Release 默认为 ON 并要求显式关闭？
- RISK-030：是否接受 API 变更（新增 `PricingResult` 重载 + 弃用 `last_total_iterations()`），还是仅改 atomic（不推荐）？
- RISK-031：Neumann 在障碍路径选"仅远场外推实现"还是"显式 throw"？`PayOff::kind()` 是否允许自定义 payoff override？
- RISK-035：方案 ①（增 statistic）还是 ②（statistic=null）？
- RISK-036：是否接受触碰 `core/math.hpp`（推荐 ①），或保持现状？

**未覆盖**：RISK-029 未实测 aarch64 构建（本机无 ARM 环境，结论基于 GCC 文档）；RISK-031 数值证据未独立复现；RISK-035/036 的 golden/测试改写需实施时验证。

---

**编制**: 主 agent 汇总 5 路调研 + 源码复核 · **交叉引用**: [EXTERNAL_REVIEW_2026-10-20.md](../../docs/audit/EXTERNAL_REVIEW_2026-10-20.md) §二 · [AUDIT.md](./AUDIT.md) · [CHECKLIST.md](./CHECKLIST.md)
