# PHASE9_SPEC: v1.9 KDE/KNN + KSG/TE 信息论基础设施

> **版本**: v1.1 (2026-11-11: spec 审计修订 — A1 C++17→**C++20** / A2 normal_reference C 常数改 gen 冻结 (statsmodels 实测 ≈1.05922, 非 1.059) / A3 Gaussian MI 容差三处统一 (20 seeds, 红线 0.02 nat) / B1-B2 测试判据可执行化 / B3 n_shuffles 拆分 ETE 与块置换 / B4 MI 解析锚改纯双变量 Gaussian / C1-C7 风险表与开放问题补齐; v1.0 草拟)
> **版本归属**: **v1.9** (Phase 9)
> **目标**: 补齐被两次跳过的非参数密度与信息论基础设施 — 1D KDE + kNN 密度/搜索 + KSG 互信息 + 迁移熵 (TE), 为 Factor Fingerprint 的因果失效诊断与 v2.0 信息论事前度量铺路
> **上游裁决**: [PHASE8_SCOPE_RESEARCH.md](../../research/PHASE8_SCOPE_RESEARCH.md) §5.2.5 D-2 = (b) KDE/KNN + KSG/TE, 载体 v1.9
> **阶段 0 产出**: [PHASE9_DOMAIN_RESEARCH.md](../../research/PHASE9_DOMAIN_RESEARCH.md) v1.0 (A 67 / B 4 / C 7 / H 5; R1-R4 通过)
> **Probe 前置**: [probe_sklearn_kde.py](../../../tests/fixtures/timeseries/probe_sklearn_kde.py) 8/8 PASS; [probe_ksg.py](../../../tests/fixtures/timeseries/probe_ksg.py) 9 PASS / 1 SKIP
> **状态**: 草拟, **未冻结**; 未改动任何 C++ 代码

---

## 0. 门禁记录 (R 清零 → spec 冻结)

| 门禁 | 状态 | 记录 |
|------|------|------|
| R1 断言统计 | ✅ | A 67 / B 4 / C 7 / H 5 (合计 83); 统计表见调研报告 §0.2 |
| R2 A 类证据可核验 | ✅ | 外部 21 条 + 本地源码/运行时 46 条; 独立 agent 盲审 7 项关键事实 7/7 PASS |
| R3 B 类登记+审计 | ✅ | 附录 B 机读块; `scripts/assertion_audit.py` 执行结果: B1-B3 NO_PROBE (文档/解析类), B4 PROBE_SURVIVED |
| R4 阻断性清零 | ✅ | B1-B4 无 FALSIFIED/CONFLICT; H1-H5 保留为开放问题, 不进入阻断路径 |

**本 spec 引用的全部阻断性断言均来自调研报告 §0/附录 A**, 实施期发现冲突时回溯调研报告修订, 不得静默改公式。

---

## 1. 交付物清单

### 1.1 新增文件结构 (ADR-017: `include/cpphub/timeseries/`, namespace `cpphub::v1::timeseries`)

```
include/cpphub/timeseries/density/
├── bandwidth.hpp          # 新增: 带宽类型 (AbsoluteH / DimensionlessFactor / Rule)
├── kde.hpp                # 新增: 1D Gaussian KDE + 选择器 + 敏感性网格
└── knn_density.hpp        # 新增: 暴力 Chebyshev kNN 搜索 + Loftsgaarden-Quesenberry 密度

include/cpphub/timeseries/information/
├── ksg.hpp                # 新增: KSG Algorithm 1 MI/CMI (Chebyshev, digamma)
└── transfer_entropy.hpp   # 新增: TE = CMI(Y_t; X_past | Y_past) + effective TE + 块置换

tests/unit/timeseries/
├── test_kde.cpp
├── test_knn_density.cpp
├── test_ksg.cpp
├── test_transfer_entropy.cpp
└── test_integration_phase9.cpp

tests/fixtures/timeseries/
├── probe_sklearn_kde.py   # ✅ 已落地 (8/8 PASS)
├── probe_ksg.py           # ✅ 已落地 (9 PASS / 1 SKIP)
└── gen_phase9_baselines.py  # 新增: 生成 C++ 测试用黄金值 (.inc / JSON)
```

**依赖**: 仅 `cpphub_core` (Real/Size/异常) + 标准库; **不依赖 Eigen3** (与 garch/unit_root 同层, header-only)。不新增 CMake target, 挂入既有 `cpphub_timeseries` INTERFACE 库。

### 1.2 新增测试套件

| 套件 | 覆盖 |
|------|------|
| `test_kde.cpp` | 带宽规则/类型、pdf/log_pdf、与 scipy 黄金值、敏感性网格、边界输入 |
| `test_knn_density.cpp` | 暴力参考实现逐位、k 边界、Chebyshev 距离、密度积分近似 |
| `test_ksg.cpp` | Gaussian 解析 MI、sklearn 黄金值、缩放不变性、ties/jitter、负值截断 |
| `test_transfer_entropy.cpp` | Gaussian VAR TE=GC/2、effective TE、块置换 p 值、缩放不变性、零假设 |
| `test_integration_phase9.cpp` | KDE → KSG → TE 端到端 + 与 7C Granger 内部互检 |

### 1.3 必须达到的数值基准 (容差分层)

| 基准 | 容差 | 验收方式 |
|------|------|---------|
| KDE pdf vs scipy `gaussian_kde` (同核同 h) | **1e-14** (probe 实测 7.5e-16) | `test_kde.cpp` 黄金值 |
| 带宽规则 vs 公式 (Scott/Silverman/normal_reference) | **1e-15** | `test_kde.cpp` |
| R `bw.nrd0`/`bw.nrd` 语义交叉 | 相对差 ≤ 1e-2 (R 用 IQR/1.34, 本库用 1.349) | `gen_phase9_baselines.py` 可选 R 检查 |
| kNN 密度 vs 暴力参考 (同算法) | **1e-15** | `test_knn_density.cpp` |
| KSG MI vs sklearn (同数据, 无 jitter/standardize) | **1e-12** vs `gen_phase9_baselines.py` 冻结黄金值; sklearn 私有 API (`_compute_mi_cc`) 仅生成期使用, 不进运行时测试 (审计 C4) | `test_ksg.cpp` |
| KSG MI vs Gaussian 解析 (n=2000, k=3, **20 seeds**) | 均值偏差 > **0.02 nat** 判红; 0.01 nat 为系统性偏差记录线, 不判红 (审计 A3: probe K1 5-seed SD 0.015-0.021, 0.01 均值断言 flaky) | `test_ksg.cpp` |
| TE vs GC/2 (Gaussian VAR φ=0.5, a=0.3, n=2000) | **±30%** | `test_transfer_entropy.cpp` |
| TE 标准化后缩放不变性 (x×0.1/1/10) | **1e-12** | `test_transfer_entropy.cpp` |
| 零假设校准 (i.i.d. MI≈0; AR(1) TE 非零) | 见 §6.3 | `test_transfer_entropy.cpp` |
| 性能: KSG N=5000, k=3, 单线程 | **< 0.5 s** | 独立 perf 用例 |

### 1.4 接口总则

- 所有公共类型/函数落 `cpphub::v1::timeseries::density` / `::information`; 头文件 `#pragma once`; header-only inline。
- 输入用 `std::span<const Real>`, 不复制调用方数据; 输出 `Real`/`std::vector<Real>`/小 struct。
- 参数校验在边界 (空样本、n<2、k<1、k≥n、非有限值) 抛 `std::invalid_argument`; 内部不做重复校验。
- 异常消息统一前缀 `KDE: ` / `KNN: ` / `KSG: ` / `TE: `。
- 不引入 RNG 全局状态; 需要随机性时显式 `seed`。

---

## 2. Scope 声明 (D9 决策)

| 编号 | 决策 | 理由 |
|------|------|------|
| D9-1 | KDE **仅 1D Gaussian 核**; 其他核 (Epanechnikov/tophat 等) P2 | v1.9 目标是"基准锚 + 可视化", 1D Gaussian 足以; 多核是低收益扩展 |
| D9-2 | 带宽用显式类型: `AbsoluteH` / `DimensionlessFactor` / `Rule`; 规则仅 Scott/Silverman/normal_reference; **SJ/ISJ plug-in P2** | 调研 A1-A25 证明同名不同物; SJ 实现成本高, 先用规则 + 敏感性网格 |
| D9-3 | kNN 搜索 **暴力 Chebyshev O(N²)**; k-d tree P2 (风险 ③) | 先正确后性能; 文档显式声明复杂度 |
| D9-4 | KSG **仅 Algorithm 1** (与 sklearn 一致); Algorithm 2 P2 (H5) | probe 三变体未复现解析值, 语义未冻结 |
| D9-5 | TE = KSG-CMI + 嵌入 (k,l,delay); 单位 **nats 规范**, 提供 bits 视图; effective TE + 块置换 surrogate | 调研 §4.1/§4.5; nats 与 Barnett TE=½GC 对齐 |
| D9-6 | KSG/TE **默认标准化** (内部 z-score; 关闭须显式 `standardize=false` 且文档警告缩放敏感性) | Chebyshev 混合单位缩放敏感性 ±35% (probe K3/A41); 审计 C1 统一"强制"措辞 |
| D9-7 | 不做多变量 KDE / 多变量条件 TE / 离散 TE / Rényi TE | v1.9 聚焦连续双变量; 多变量条件化 P2 |
| D9-8 | 不做边界修正 (正支撑数据) / 权重 KDE | 边界偏倚列入开放问题; 权重 P2 |
| D9-9 | 不新增 C ABI / Python 绑定 | Phase 9 纯 C++ 内部模块; 未来按需 |
| D9-10 | 性能预算: KDE N=10k/M=1k < 50 ms; KSG N=5k < 0.5 s; TE N=2k, ETE 100 + 块置换 100 (**共 ~200 次 TE 评估**) < 10 s | 单线程, MSVC Release; 审计 C6 明确预算含 p 值置换。**[2026-11-11 G 门禁修订, 经用户批准]** KDE 预算 50→150 ms: 实测 52-58 ms (UCRT 标量 exp ~5.5 ns/次硬下限, 暴力算法已最优; exp2 变换实测慢 8× 禁用), 150 ms 留 CI runner 1.5-2× 慢速余量; KSG 实测 35 ms / TE 实测 4.1 s 均达标不动。**[2026-11-14 CI 修订]** 预算口径仅适用于优化构建: 覆盖率 CI job 为 Debug + `--coverage -O0` 插桩构建 (run 43465160758: KSG 1.3 s / TE 98.9 s 超预算) → `tests/CMakeLists.txt` 检测 `--coverage` 定义 `CPPHUB_COVERAGE_BUILD`, 三个 perf 用例在该构建下 `GTEST_SKIP()` (Release 仍断言预算) |

**不做清单**: 多核 KDE / 多维 KDE / 带宽矩阵 / SJ/ISJ / 边界修正 / 加权 KDE / k-d tree / KSG Algorithm 2 / 多变量 TE / 离散或 Rényi TE / 排列熵 / 滑窗 TE API (窗口切片由调用方实现, 审计 C5) / C ABI / Python 绑定。

---

## 3. M1: 1D KDE (`density/kde.hpp`)

### 3.1 带宽类型 (`density/bandwidth.hpp`)

```cpp
namespace cpphub::v1::timeseries::density {

enum class BandwidthRule { kUser, kScott, kSilverman, kNormalReference };

struct Bandwidth {
    BandwidthRule rule = BandwidthRule::kScott;
    Real value = 0.0;        // kUser: 绝对 h (is_factor=false) 或无量纲因子 (is_factor=true)
    bool is_factor = false;  // kUser 时有效
};

/// 规则 → 绝对 h (数据单位). 样本 sd 用 ddof=1; A = min(sd, IQR/1.349).
Real resolve_bandwidth(std::span<const Real> sample, const Bandwidth& bw);

}  // namespace density
```

**规则公式** (冻结, 对应调研 A11-A14/A19-A21):

- `kScott`: `1.059 * A * n^(-1/5)`
- `kSilverman`: `0.9 * A * n^(-1/5)`
- `kNormalReference`: `C * A * n^(-1/5)`, Gaussian 核 `C` **以 `gen_phase9_baselines.py` 从 statsmodels 实测冻结的精确值为准 (≈1.05922, Hansen 2009; 审计 A2: 并非 1.059 — 若硬编码 1.059 将与 statsmodels 产生 0.02% 系统偏差且与 kScott 重合)**; `kScott` 恒取 1.059 (statsmodels 同款硬编码)。两枚举数值仅差 0.02%, 保留独立枚举以便未来换核。
- `A = min(std(sample, ddof=1), IQR/1.349)`, `IQR = q75 - q25` (线性插值分位数)
- `kUser + is_factor=true`: `h = value * std(sample, ddof=1)` (scipy 语义; 显式)
- `kUser + is_factor=false`: `h = value` (statsmodels/R 语义)

**禁止**: 不接受裸字符串; 不提供 sklearn 的"裸 n^(-1/(d+4)) 当 h"语义 (调研 A1-A3)。

### 3.2 KDE 类

```cpp
struct KdeConfig {
    Bandwidth bandwidth{};
    bool reflection = false;  // P2 预留; v1.9 必须 false, true 抛 not-implemented
};

class KernelDensity1D {
public:
    explicit KernelDensity1D(std::span<const Real> sample, KdeConfig config = {});

    Real pdf(Real x) const;                       // 1/(n h) Σ φ((x - x_i)/h)
    Real log_pdf(Real x) const;                   // log-sum-exp, 避免下溢
    std::vector<Real> pdf(std::span<const Real> grid) const;
    Real bandwidth() const noexcept;              // 解析后的绝对 h
    BandwidthRule rule() const noexcept;
    std::size_t size() const noexcept;

private:
    std::vector<Real> sample_;
    Real h_;
    Real log_norm_;  // -log(n h sqrt(2π))
};

/// 多带宽敏感性: 对每个 h 返回 log-likelihood 与给定网格上的 pdf 摘要.
struct BandwidthSensitivityPoint {
    Real bandwidth;
    Real log_likelihood;   // Σ log pdf(x_i) (留一法可选, P2)
    Real pdf_at_mode;      // 网格最大值 (诊断用)
};
std::vector<BandwidthSensitivityPoint> bandwidth_sensitivity(
    std::span<const Real> sample,
    std::span<const Real> bandwidths,
    std::span<const Real> eval_grid);
```

**实现要点**:

- `pdf`: 直接 O(N) 求和; `log_pdf` 用 log-sum-exp (max-shift), 保证 h→0 或深尾不返回 `-inf` 以外的 NaN。
- `bandwidth_sensitivity`: 每个 h 重新构造 KDE; 返回的 `log_likelihood` 不含留一修正 (P2)。
- 输入校验: `n >= 2`, 所有值 finite, `h > 0`; 否则抛 `std::invalid_argument("KDE: ...")`。
- 不存储排序样本 (1D KDE 无需树)。

### 3.3 测试要点

1. 四个规则 vs 公式, 1e-15。
2. `kUser + factor` vs scipy `gaussian_kde(bw_method=factor)`, pdf 1e-14。
3. `kUser + absolute` vs 手写高斯和, 1e-15。
4. `log_pdf` 与 `log(pdf)` 一致 (pdf>0 时 1e-12)。
5. 敏感性网格: 3 个 h 的 `log_likelihood` 有限; 极窄/极宽 h 不 NaN。
6. 非法输入: 空样本、n=1、h=0、NaN → `invalid_argument`。
7. `reflection=true` → `not_implemented` (P2 显式)。

---

## 4. M2: kNN 密度 + KSG MI (`density/knn_density.hpp`, `information/ksg.hpp`)

### 4.1 kNN 搜索与密度

```cpp
namespace cpphub::v1::timeseries::density {

struct KnnConfig {
    std::size_t k = 5;
};

class KnnDensity1D {
public:
    explicit KnnDensity1D(std::span<const Real> sample, KnnConfig config = {});

    /// Loftsgaarden-Quesenberry: k / (n * 2 r_k(x)), r_k = 到第 k 近邻的
    /// Chebyshev 距离 (1D 即绝对距离). 边界外按 r_k 截断, 不做归一化.
    Real pdf(Real x) const;
    std::vector<Real> pdf(std::span<const Real> grid) const;
    std::size_t k() const noexcept;

private:
    std::vector<Real> sorted_;  // 升序, 用于 lower_bound 定位
    std::size_t k_;
    Real n_;
};

}  // namespace density
```

**实现要点**:

- 构造时复制并排序样本; `pdf(x)` 用 `lower_bound` 找插入位置, 在左右各取 k 个候选, 计算第 k 小绝对距离 (O(k), k 通常 ≤ 20)。
- `x` 落在样本范围外时, 半径取到端点距离; 文档声明"边界外密度不归一"。
- 与 `KnnSearch1D` 合并 (v1.9 不单独暴露搜索类; KSG 用自己的实现, 见 §4.2)。

### 4.2 KSG MI / CMI

```cpp
namespace cpphub::v1::timeseries::information {

struct KsgConfig {
    std::size_t k = 3;
    bool standardize = true;   // z-score each input (ddof=1)
    double jitter_scale = 1e-10;  // 相对 mean|x|; 0 = 关闭
    std::uint64_t seed = 0;
};

/// KSG Algorithm 1 (Eq. 8): psi(k) + psi(N) - <psi(nx+1)+psi(ny+1)>,
/// joint Chebyshev, strict < eps(i). 返回 nats.
Real ksg_mutual_info(std::span<const Real> x, std::span<const Real> y,
                     const KsgConfig& config = {});

/// Frenzel-Pompe CMI: psi(k) - <psi(nxz+1)+psi(nyz+1)-psi(nz+1)>, nats.
Real ksg_conditional_mutual_info(std::span<const Real> x,
                                 std::span<const Real> y,
                                 std::span<const Real> z,
                                 const KsgConfig& config = {});

}  // namespace information
```

**实现要点**:

- 距离: **Chebyshev** (max norm), 联合空间 `Z=(x,y)` 或 `(x,y,z)`。
- `eps(i)` = 到第 k 近邻的联合距离; 边际计数用 **严格小于** `eps(i)` (实现上用 `nextafter(eps, 0)` 或 `<` 比较)。
- 计数减 1 (排除自身); `n_x`/`n_y` 至少 0。
- 公式: `ψ(k) + ψ(N) - mean(ψ(nx+1)+ψ(ny+1))`; CMI 用 `ψ(k) - mean(ψ(nxz+1)+ψ(nyz+1)-ψ(nz+1))`。
- `standardize=true` 时内部 z-score (ddof=1); `false` 时文档警告缩放敏感性。
- `jitter_scale > 0` 时加 `jitter_scale * mean(|x|) * N(0,1)` (每维独立, seed 固定); 用于 ties。
- 返回原始值, **不截断负值**; 调用方/测试决定是否 `max(0, ·)`。
- 复杂度 O(N²) (暴力); N=5000, k=3 预算 < 0.5 s。

### 4.3 测试要点

1. 与 `gen_phase9_baselines.py` 冻结的 sklearn 黄金值逐位一致 (1e-12; 黄金值生成时关闭 jitter/standardize 以对齐 sklearn `_compute_mi_cc` 原始路径; 审计 C4)。
2. Gaussian 解析 MI: n=2000, k=3, ρ=0.5/0.8, **20 seeds**, 均值偏差 > 0.02 nat 判红 (0.01 记录线), SD 记录 (审计 A3)。
3. 标准化后 x×0.1/1/10 的 MI 不变 (1e-12); 关闭标准化时变化 > 0.002 (文档化陷阱)。
4. ties: 舍入数据 + jitter 改变估计; 不崩溃。
5. 负值: 独立小样本返回可能为负; 确认不截断 (与 sklearn 文档语义区分)。
6. 非法输入: n<k+1、k<1、非有限值、空 z。
7. kNN 密度 vs 暴力参考 (同排序/同 k) 1e-15; k=n-1 边界。

---

## 5. M3: 迁移熵 (`information/transfer_entropy.hpp`)

### 5.1 定义与接口

```cpp
namespace cpphub::v1::timeseries::information {

struct TransferEntropyConfig {
    std::size_t source_history = 1;  // k
    std::size_t target_history = 1;  // l
    std::size_t delay = 1;
    std::size_t ksg_neighbors = 3;
    bool standardize = true;
    std::size_t n_shuffles = 100;    // effective TE: 随机置换 (Marschinski-Kantz 2002)
    std::size_t n_block_boot = 100;  // p_value: 块置换 (Dimpfl-Peter 2013, 保目标自相关)
                                     // 审计 B3: 两个零分布不同, 各自独立抽样, 不得复用
    std::size_t block_length = 0;    // 0 => ceil(sqrt(T)); 块置换 surrogate
    std::uint64_t seed = 0;
};

struct TransferEntropyResult {
    Real te = 0.0;            // nats
    Real te_effective = 0.0;  // te - mean(shuffled te)
    Real te_bits = 0.0;       // te / ln 2
    Real p_value = 1.0;       // 块置换零分布
    std::size_t n_effective = 0;  // T - max(k,l)*delay
};

TransferEntropyResult transfer_entropy(std::span<const Real> source,
                                       std::span<const Real> target,
                                       const TransferEntropyConfig& config = {});

}  // namespace information
```

**定义** (Barnett 2009 Eq. 10 / Schreiber 2000):

`TE(X→Y) = CMI(Y_t ; X_{t-delay}, ..., X_{t-k*delay} | Y_{t-delay}, ..., Y_{t-l*delay})`

- 嵌入向量: 源 `k` 阶、目标 `l` 阶, 间隔 `delay`; 默认 `k=l=1, delay=1` (与 RTransferEntropy 默认一致)。
- **内部调用顺序固定为 `CMI(target_t, source_past, target_past)`** (即 `I(Y_t; X_past | Y_past)`), 避免把 source/target 写反; 单测覆盖 `TE(X→Y)` 与 `TE(Y→X)` 不等。
- 有效样本 `T - max(k,l)*delay`; 若 < `ksg_neighbors+1` 抛 `invalid_argument`。
- 标准化: 对 `source`/`target` 原始序列 z-score (ddof=1), 再构造嵌入。
- `te_effective`: 对 source 做 `n_shuffles` 次**随机置换** (Marschinski-Kantz, 破坏全部依赖含源自相关), 减去 shuffled TE 均值; **仅辅助诊断, 不作推断依据**。
- `p_value`: **块置换** (circular block shuffle, 保目标自相关) 抽样 `n_block_boot` 次, 与 `n_shuffles` 相互独立 (审计 B3); 块长默认 `ceil(sqrt(T))`; `p = (1 + #{null >= te}) / (1 + n_block_boot)`。
- `te_bits`: `te / ln 2`; 不改变 `te` 的 nats 语义。
- **不提供** "因果" 字段; 文档注明 TE 是条件 MI 的估计, 不是因果发现工具 (James-Barnett-Crutchfield 2016; Smirnov 2013)。

### 5.2 与 7C Granger 的内部互检

- 集成测试构造 Gaussian VAR(φ=0.5, a=0.3), 用 7C `granger_causality` 的 OLS/VAR 路径计算 `GC = ln(σ²_restricted/σ²_unrestricted)`。
- 断言 `TE / GC ≈ 0.5` (nats), 容差 ±30% @ n=2000 (probe K4: +25.0%)。
- 若 7C 接口只返回 F/TY/HAC-Wald 统计量, 测试侧先换算成 ln 方差比; 不在 TE 模块内做换算。

### 5.3 测试要点

1. Gaussian VAR TE=GC/2 (5 seeds, 均值 ±30%; SD 记录)。
2. `te_bits = te / ln2` (1e-15)。
3. 标准化后缩放不变性 (1e-12); 关闭标准化时变化 (文档化)。
4. effective TE: 独立 AR(1) 零假设下**记录型断言, 不判红** (审计 B1: probe K5 实测 naive/effective 同为负, 符号随种子波动, "\|eff\|<\|naive\|"不可靠); 零假设推断以 p_value 为准。
5. p 值: (a) i.i.d. 数据 20 条, 断言 p>0.05 的比例 ≥ 60%; (b) 耦合 VAR 数据 5 seeds, 断言 p < 0.1; (c) p 值分布均匀性 (KS) 留手动检查, 不进单测 (审计 B2: 100 shuffles 粒度 0.01, 均匀性检验成本 ~200×200 次 TE)。
6. 嵌入边界: k=l=1/delay=1; k=3/l=2/delay=2; T 过小抛错。
7. 零假设: i.i.d. MI≈0; AR(1) TE 非零 (不宣称正负)。

---

## 6. 集成与验证

### 6.1 端到端 (`test_integration_phase9.cpp`)

1. 生成 Gaussian VAR(φ=0.5, a=0.3, n=2000)。
2. KDE 拟合 x 的边际, 与 scipy 黄金值 1e-14。
3. KSG MI(x_t, y_t) vs 解析 `-½ln(1-ρ²)` — **ρ 取纯双变量 Gaussian fixture 的相关系数 (probe K1 同款), 不从 VAR 反解** (审计 B4: VAR 平稳同期相关 ≠ 新息相关, 有交叉滞后贡献); 偏差 ≤ 0.02 nat (20 seeds 均值)。
4. TE(x→y) vs GC/2, ±30%。
5. KDE 多带宽敏感性: 3 个 h 的 `log_likelihood` 单调性/稳定性报告。
6. 与 7C Granger 的 GC 对照 (见 §5.2)。

### 6.2 跨库验证

- Python 侧: `probe_sklearn_kde.py` (8/8) + `probe_ksg.py` (9 PASS / 1 SKIP) 已冻结语义。
- R 侧: base R `bw.nrd0/bw.nrd` 已由 probe P7 验证; `ks` 包未安装, `hpi/hscv` 留 P2 (H1)。
- JIDT 未安装; TE 的 JIDT 对照留 verify 阶段 (H2)。

### 6.3 红绿判据

- **红**: KDE 同 h pdf 差 > 1e-14; KSG vs 黄金值 > 1e-12; Gaussian MI **20 seeds 均值**偏差 > 0.02 nat; TE/GC 偏差 > 30%; 标准化后 TE 变化 > 1e-12; p 值判据 (§5.3-5a/b) 违反。
- **绿**: 上述全部通过, 且 perf 预算达标。

---

## 7. 验收标准 (G 门禁)

### 7.1 编译与测试

- [x] MSVC 零警告, **C++20** (仓库 `CMAKE_CXX_STANDARD 20`; §1.4 的 `std::span` 为 C++20 特性)。
      实测 (2026-11-11): 全量 Release 构建 EXIT=0, 全部警告均为既有项 (aad_greeks C4244 / rbf_fd
      C4996), Phase 9 文件零警告; **Phase 9 新代码定向 `/W4` 5 个 TU 全 CLEAN** (仓库全局基线
      为 /W3, spec 的 /W4 落实为新增代码定向验证, 不改全局旗标)。
      GCC 平台随 CI 双平台 job 覆盖 (本门禁批次 MSVC 实测; A/B 站 GCC 惯例随合并后 bundle 验证)。
- [x] 新增 5 个测试套件全绿: **41/41** (test_kde 11 + test_knn_density 6 + test_ksg 10 +
      test_transfer_entropy 12 + test_integration_phase9 2)。
      **全量回归**: 本地 MSVC Release **2671/2671 全绿, 零失败** (2026-11-11 17:57 终态;
      CUDA/Python OFF, C ABI ON)。**偏差记录**: spec 预估 "+60~80 例 / 全库 ≥2700" 作废 —
      实际新增 41 例, 且 ≥2700 的推算基于 CI 口径基线 2643, 未计本地 CMake 配置差异;
      判据以"5 套件全绿 + 全量零失败 + 计数如实记录"为准 (同 KDE 预算修订模式)。
- [x] 无 C ABI / Python 变更 (git status: 唯一修改的既有文件 = tests/CMakeLists.txt +17 行
      测试注册); 无 Eigen3 依赖 (density/ + information/ 头文件 grep 零命中)。

### 7.2 数值基准

- [x] §1.3 全表通过; 黄金值来源 `gen_phase9_baselines.py` (复跑输出与测试内嵌常量逐位一致,
      仅两处 ULP 级末位差来自打印端乘除顺序, 均在容差远端) + probe 输出。
- [x] 关键锚实测: sklearn/scipy KDE **1e-14** (逐位); KSG 黄金值 **1e-12** (逐位); Gaussian MI
      **20 seeds 均值** ≤ 0.02 nat (§6.1.3/§6.3 口径, 实测通过); TE=GC/2 **±30%** (5 seeds 单测 +
      3 seeds 集成, 固定 seed 确定性); 标准化缩放不变 1e-12。

### 7.3 文档与流程

- [x] `docs/DEVELOPMENT_LOG.md` 新增 Phase 9 条目。
- [x] `README.md` 模块表/特性列表更新 (KDE/KNN/KSG/TE; Status 表 + Features 各一条)。
- [x] 调研报告 §8 的 probe 结果与 spec 数值一致 (spec 冻结值均出自 probe 报告 §8 回填, 同源)。

### 7.4 Scope 边界

- [x] 不做清单 (§2) 零越界: 无 k-d tree / 多维 KDE / Alg2 / 加权 KDE / 边界修正 (reflection
      显式抛 not-implemented) / C ABI / Python 绑定; P2 项以 `[P2]` 标注 (bandwidth.hpp H1/
      kde.hpp reflection / ksg.hpp H5 / Q11 块长警告), 不进入本期测试。
- [x] TE 无 "因果" 字段 (§5.1 文档义务); effective TE 注释声明"仅诊断, 推断以 p 值为准" (审计 B1)。

### 7.5 G 门禁清单 (G1-G8, 2026-11-11 执行)

> G1-G8 由 §7.1-7.4 + §9 步 7 (性能用例) + §1.1 交付物推导 (spec 原文引用清单但未枚举, 本节
> 补齐映射 — 同 Phase 7C ACCEPTANCE_CHECKLIST 模式)。

| # | 检查项 | 来源 | 结果 | 实测 |
|---|--------|------|------|------|
| G1 | 编译零警告 (C++20) | §7.1 | [x] | 全量 Release EXIT=0, 警告全为既有项; Phase 9 定向 /W4 5 TU CLEAN |
| G2 | 5 套件全绿 + 计数 | §7.1 | [x] | 41/41; 全量本地 MSVC **2671/2671 零失败**; 预估 ≥2700 作废 (见 §7.1 偏差记录) |
| G3 | 数值基准 §1.3 全表 | §7.2 | [x] | KDE 1e-14 / KSG 1e-12 / MI 20 seeds ≤0.02 / TE ±30% / 缩放 1e-12 |
| G4 | 性能预算 (D9-10, KDE 修订 150ms) | §9 步 7 | [x] | KDE 52-58ms (<150) / KSG 35ms (<500) / TE 4.1s+校准用例 (<10s 预算含置换) |
| G5 | 交付物完整 (4 头 + 5 套件 + gen 脚本 + 2 probe) | §1.1 | [x] | gen_phase9_baselines.py 复跑黄金值逐位一致; probe 8/8 + 9/1 冻结 |
| G6 | 文档 (DEVELOPMENT_LOG + README) | §7.3 | [x] | 均已更新 (Phase 9 条目 / Status 表 + Features) |
| G7 | Scope 零越界 + P2 标注 | §7.4 | [x] | 不做清单 grep 验证; reflection 显式 not-implemented / H1 调研报告登记 / Q11 块长显式参数注记 / Alg2 未实现 (ksg.hpp 标题限定 Alg1) |
| G8 | 测试-spec 口径对齐 (审计 B1/B2/C4/C6 落实) | §5.3/§6.3 | [x] | MI 20 seeds/0.02; TE ±30%; p 值 5a 比例 ≥60% + 5b n=4000 p<0.1; sklearn 私有 API 仅生成期 |

门禁结论: **G1-G8 全绿, 允许合并** (合并前置: 全库 ctest 终态记录 + 用户批准提交)。

---

## 8. 风险与缓解 (继承调研 §5 + 新增)

| # | 风险 | 缓解 |
|---|------|------|
| ① | 带宽语义混用 | 类型系统 (§3.1) + 规则公式冻结 + 敏感性网格 |
| ② | KSG Chebyshev 缩放敏感 | 默认标准化 (D9-6) + 缩放不变性测试 |
| ③ | k-d tree 缺失导致性能 | 暴力先行 + 文档声明 O(N²); k-d tree P2 |
| ④ | TE 偏差/naive 分叉 | effective TE (辅助) + 块置换 p 值 (推断主依据) + probe K5 |
| ⑤ | TE-GC 因子/单位 | nats 规范 + `te_bits` 视图 + 集成测试 0.5 锚 |
| ⑥ | 窗口长度=平稳性声明 | `TransferEntropyConfig` 不设默认窗口; 文档要求调用方显式滑窗并记录 |
| ⑦ | KDE 带宽敏感性 | `bandwidth_sensitivity` + 集成测试报告 |
| ⑧ | statsmodels 文档/实现不一致 | 不作为依赖; 仅 probe 对照 |
| ⑨ | sklearn 字符串语义 | 不复制; 类型层拒绝裸规则字符串 |
| ⑩ | nats/bits 混用 | 结果 struct 明示单位; `te_bits` 冗余字段 |
| ⑪ | FFT 网格 vs 点评估 | 本实现直接求值, 不涉及 FFT |
| ⑫ | Alg2 语义未冻结 | 仅 Alg1; H5 留 P2 |
| ⑬ | 边界偏倚 (正支撑) | 文档警告; 边界修正 P2 (开放问题) |
| ⑭ | TE 组装层无第二实现对照 (CMI 核有 sklearn 1e-12 锚, 但嵌入/差分/有效样本裁剪仅 GC/2 ±30% 一个锚; JIDT 对照推迟 = H2) | 接受为 v1.9 残余缺口; 集成测试加 TE(X→Y)≠TE(Y→X) 方向性断言; JIDT 对照进 verify 阶段优先级最高 (审计 C2) |

---

## 9. 实施顺序

| 步 | 内容 | 依赖 | 验收 |
|----|------|------|------|
| 1 | `bandwidth.hpp` + `kde.hpp` + `test_kde.cpp` | 无 | 规则/黄金值 1e-14 |
| 2 | `knn_density.hpp` + `test_knn_density.cpp` | 无 | 暴力参考 1e-15 |
| 3 | `ksg.hpp` (MI) + `test_ksg.cpp` | 步 2 | 黄金值 1e-12; Gaussian 20 seeds 均值 ≤ 0.02 nat |
| 4 | `ksg.hpp` (CMI) + TE 骨架 | 步 3 | CMI 对称性/边界 |
| 5 | `transfer_entropy.hpp` + `test_transfer_entropy.cpp` | 步 4 | TE=GC/2 ±30%; effective/p |
| 6 | `test_integration_phase9.cpp` + `gen_phase9_baselines.py` | 步 1-5 | 端到端全绿 |
| 7 | 性能用例 + 文档 + 全量 CI | 步 6 | G 门禁 |

---

## 10. 开放问题与延后项

### 10.1 待核证据 (H 系列, 实施期/verify 期关闭)

| ID | 问题 | 处置 |
|----|------|------|
| H1 | R `ks` 的 `hpi/hscv` 与 SJ 1991 是否逐位一致 | P2; 安装 `ks` 后对照 |
| H2 | JIDT TE 运行时 `k/l/delay` 与单位行为 | verify 期最高优先 (风险 ⑭ 的关闭路径); jpype 对照 |
| H3 | sklearn `'scott'` 裸因子是否为设计意图 | 不影响本实现; 文档记录 |
| H4 | Kirkley reduced TE 是否能把 +25% 偏差压到 ±10% | P2; 实现 Eq. 12 后重跑 |
| H5 | KSG Algorithm 2 的投影半径/计数口径 | P2; 仅 Alg1 进 v1.9 |
| Q11 | 块长 ceil(√T) 在强持久序列下的 p 值校准性 (Politis-White 自动块长为文献标准) | P2; v1.9 文档警告块长为显式参数 (审计 C3) |
| Q12 | normal_reference C 的精确值 (statsmodels 实测 ≈1.05922) 与 Hansen 2009 原文常数是否逐位一致 | 实施期 gen 阶段冻结并记录来源 |

### 10.2 P2 延后项 (scope 决策的镜像, 非"问题")

> 以下各条均为 §2 已裁决 "v1.9 不做" 的登记, 供 v1.10+ 候选池检索; 语义上不是开放问题 (审计 C7)。

| ID | 延后项 | 对应决策 |
|----|--------|---------|
| P2-1 | KDE 正支撑边界修正 (反射/变换) | D9-8 / 风险 ⑬ |
| P2-2 | SJ/ISJ plug-in 带宽选择器 | D9-2 |
| P2-3 | 多变量 KDE / 多变量条件 TE | D9-7 |
| P2-4 | k-d tree / BallTree | D9-3 / 风险 ③ |
| P2-5 | 权重 KDE / neff 语义 | D9-8 |
| P2-6 | KSG Algorithm 2 | D9-4 / H5 |
| P2-7 | 多核 KDE (Epanechnikov/tophat) | D9-1 |
| P2-8 | 滑窗 TE API | 不做清单 (审计 C5) |

---

## 11. 关联文档

- [PHASE9_DOMAIN_RESEARCH.md](../../research/PHASE9_DOMAIN_RESEARCH.md) v1.0 (阶段 0 调研 + probe 结果 §8)
- [PHASE8_SCOPE_RESEARCH.md](../../research/PHASE8_SCOPE_RESEARCH.md) §5.2.5/§5.3.6 (D-2 裁决 + 7 风险)
- [PHASE8_SPEC.md](../phase8/PHASE8_SPEC.md) (章节结构与 G 门禁惯例)
- [PHASE7C_SPEC.md](../phase7/PHASE7C_SPEC.md) (Granger 内部互检接口)
- [ADR-017](../../decisions/ADR-017_TIMESERIES_NAMESPACE.md) (命名空间/物理目录)
- [DEVELOPMENT_WORKFLOW.md](../../DEVELOPMENT_WORKFLOW.md) v1.1 (阶段 0/R 门禁)
- [probe_sklearn_kde.py](../../../tests/fixtures/timeseries/probe_sklearn_kde.py) / [probe_ksg.py](../../../tests/fixtures/timeseries/probe_ksg.py)
