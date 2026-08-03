# Phase 5 执行规格书 - 高频计量经济学模块 (HFE)

> **版本归属**: v1.4 (v1.4.0 ~ v1.4.3, 分四波交付)
> **目标**: 在 Cpp_Hub 中建立工程级 C++ 高频计量经济学库, 填补 QuantLib v1.43 在 RV/RK/BPV/跳跃检验等领域的空白
> **前置**: Phase 1-4 (v1.0 + v1.4 PDE 三波) 全部通过, 当前测试基线 1268/1268
> **对标**: R `highfrequency` 0.5.3 / 1.0.0 (Boudt, Kleen, Sjørup 2022, JSS)

---

## 1. 背景与定位

### 1.1 领域空白 (来自 project_memory 实测记录)

- **QuantLib v1.43** 缺失 RV/RK/BPV/jump test 实现
- **R `highfrequency`** (80+ 函数, JSS 同行评审) 受 Rcpp 性能边界限制
- **Cpp_Hub v1.4** 定位为**首个工程级 C++ 高频计量经济学库**

### 1.2 文献基础 (严格学术依据, 排除幻觉)

本 spec 所有数学公式均来自以下已验证文献:

| 编号 | 文献 | DOI / 来源 | 用途 |
|---|---|---|---|
| [BN-S 2002] | Barndorff-Nielsen & Shephard, "Econometric Analysis of Realized Volatility" | *J. Applied Econometrics* 17, 453-475 | RV 定义与渐近理论 |
| [BN-S 2004] | Barndorff-Nielsen & Shephard, "Power and Bipower Variation with Stochastic Volatility and Jumps" | *J. Financial Econometrics* 2(1), 1-37, doi:10.1093/jjfinec/nbh001 | BPV / RV_jump 分解 |
| [BN-S 2006] | Barndorff-Nielsen & Shephard, "Econometrics of Testing for Jumps in Financial Economics Using Bipower Variation" | *J. Financial Econometrics* 4(1), 1-30 | BNS 跳跃检验统计量 |
| [ABD 2003] | Andersen, Bollerslev & Diebold, "Modeling and Forecasting Realized Volatility" | *Econometrica* 71(2), 579-625 | RV 分布与 ARFIMA 预测 |
| [H-L 2006] | Hansen & Lunde, "Realized Variance and Market Microstructure Noise" | *J. Business & Economic Statistics* 24(2), 127-161 | 微结构噪声与最优抽样 |
| [BKS 2022] | Boudt, Kleen & Sjørup, "Analyzing Intraday Financial Data in R: The highfrequency Package" | *J. Statistical Software* 104(8), 1-36, doi:10.18637/jss.v104.i08 | R 包官方 vignette, 函数定义 |
| [A-J 2009] | Aït-Sahalia & Jacod, "Testing for Jumps in a Discretely Observed Process" | *Annals of Statistics* 37(1), 184-222, doi:10.1214/07-AOS568 | AJ 跳跃检验 (对照) |

### 1.3 R `highfrequency` 包确认信息 (CRAN 实测)

- **CRAN 版本**: 0.5.3 (2018-03-03) / 1.0.0 (2022-10-26)
- **作者**: Kris Boudt, Jonathan Cornelissen, Onno Kleen, Scott Payseur, Emil Sjørup
- **License**: GPL (≥2)
- **仓库**: https://github.com/jonathancornelissen/highfrequency
- **JSS vignette**: Boudt et al. (2022), doi:10.18637/jss.v104.i08

---

## 2. 架构约束 (来自 project_memory 硬约束)

1. **模块独立**: `include/cpphub/hfecon/` 顶层模块, 与衍生品定价栈解耦, 只依赖 `core/`
2. **六层架构**: data / measures / tests / noise / models / liquidity
3. **R 对标硬约束**: 每个 HFE 函数必须与 R `highfrequency` 同名函数数值对照, 容差 1e-10 ~ 1e-12, 作为 CI gate
4. **跨平台**: MSVC + GCC (A/B 站) 三平台位精确一致
5. **性能目标**: 50+ Mtick/s (R Rcpp 的 10-50×), 复用 SIMD/OpenMP

---

## 3. 第一波 (v1.4.0) 交付项

**目标**: 建立 HFE 模块基础设施 + 核心 realized measures + BNS 跳跃检验
**预计新增测试**: ~15 个, 使总数从 1268 → ~1283

### 3.1 项 1: TAQ 数据读取器

**文件**: `include/cpphub/hfecon/data/taq_reader.hpp`

**功能**: 读取交易报价 (TAQ) 数据, 支持 CSV 和 NASDAQ ITCH 5.0 二进制格式

**ITCH 5.0 规范** (NASDAQ TotalView 官方协议):
- 二进制消息流, 每条消息 = 2 字节长度 + 消息体
- 核心消息类型 (单字符标识):
  - `S` System Event (系统事件)
  - `R` Stock Directory (股票目录)
  - `H` Stock Trading Action (交易状态)
  - `A` Add Order (添加订单, 无 MPID)
  - `F` Add Order - MPID (添加订单, 含 MPID)
  - `E` Order Executed (订单成交)
  - `C` Order Executed with Price (带价格成交)
  - `X` Order Cancel (订单取消)
  - `D` Order Delete (订单删除)
  - `U` Order Replace (订单替换)
  - `P` Trade Message (非显示交易, 含价格)
  - `Q` Cross Trade (交叉交易)
  - `B` Broken Trade (取消的交易)
- 价格字段: 4 字节有符号整数, 单位 = 1e-4 (即 price * 10000)
- 时间戳: 6 字节纳秒, 从午夜开始

**CSV 格式** (highfrequency 兼容):
```
DT,PRICE,SIZE,BID,ASK,BIDSIZE,ASKSIZE
2024-01-02 09:30:00.123456,100.05,100,100.04,100.06,500,300
```

**接口设计**:
```cpp
namespace cpphub::hfecon {

struct Trade {
    Timestamp ts;       // 纳秒时间戳 (core/types.hpp)
    Real price;
    Real size;
    std::string mpid;   // 可选, 做市商 ID
};

struct Quote {
    Timestamp ts;
    Real bid_price;
    Real ask_price;
    Real bid_size;
    Real ask_size;
};

class TaqReader {
public:
    // CSV 读取 (highfrequency 兼容格式)
    static std::vector<Trade> read_trades_csv(const std::string& path);
    static std::vector<Quote> read_quotes_csv(const std::string& path);

    // ITCH 5.0 二进制读取 (NASDAQ TotalView)
    static std::vector<Trade> read_trades_itch(const std::string& path);
    static std::vector<Quote> read_quotes_itch(const std::string& path,
                                                bool reconstruct_book = true);

    // 时间聚合 (对应 R aggregatePrice)
    static std::vector<Trade> aggregate_price(
        const std::vector<Trade>& trades,
        Size align_period,        // e.g. 5
        const std::string& align_by = "minutes",  // "seconds"/"minutes"/"ticks"
        const std::string& market_open = "09:30:00",
        const std::string& market_close = "16:00:00",
        bool fill = false);       // 是否填充缺失时段

    // 收益率生成 (对应 R makeReturns)
    static std::vector<Real> make_returns(const std::vector<Trade>& prices);
};

} // namespace cpphub::hfecon
```

**R 对照函数** (无, 基础设施): R `aggregatePrice` / `makeReturns` 作为行为对照

### 3.2 项 2: 核心 Realized Measures

**文件**: `include/cpphub/hfecon/measures/realized_measures.hpp`

**数学定义** (严格来自 [BN-S 2002, 2004], [BKS 2022] §3):

设日内对数收益率为 $r_1, r_2, \ldots, r_n$, 其中 $r_i = \log P_{t_i} - \log P_{t_{i-1}}$。

**Realized Variance (RV)** [BN-S 2002]:
$$RV = \sum_{i=1}^{n} r_i^2$$
- R 函数: `rRVar(rdata)` 或 `rCov(rdata, makeReturns=TRUE)` 的对角线
- 性质: $RV \xrightarrow{p} \int_0^T \sigma^2(u) du$ (积分方差 IV) 当 $n \to \infty$

**Realized Volatility (RVol)**:
$$RVol = \sqrt{RV}$$
- R 函数: `rRealizedVolatility(rdata)`

**Realized Quarticity (RQ)** [BN-S 2002, BN-S 2004]:
$$RQ = \frac{n}{3} \sum_{i=1}^{n} r_i^4$$
- R 函数: `rQuar(rdata)`
- 性质: $RQ \xrightarrow{p} \frac{2}{3} \int_0^T \sigma^4(u) du$ (用于跳跃检验的 IQV 估计)

**Bipower Variation (BPV)** [BN-S 2004]:
$$BPV = \frac{n}{n-1} \sum_{i=2}^{n} |r_{i-1}| |r_i|$$
- R 函数: `rBPCov(rdata, makeReturns=TRUE)`
- 性质: $BPV \xrightarrow{p} \int_0^T \sigma^2(u) du$ (IV), **对罕见跳跃稳健**

**Realized Semivariance (RSV+ / RSV-)** [Barndorff-Nielsen, Kinnebrock, Shephard 2010]:
$$RSV^+ = \sum_{i=1}^{n} r_i^2 \mathbb{1}_{\{r_i > 0\}}, \quad RSV^- = \sum_{i=1}^{n} r_i^2 \mathbb{1}_{\{r_i \leq 0\}}$$
- R 函数: `rSV(rdata)` 返回 (RSV+, RSV-)
- 性质: $RSV^+ + RSV^- = RV$, $RSV^+ \xrightarrow{p} \frac{1}{2} IV$

**接口设计**:
```cpp
namespace cpphub::hfecon {

struct RealizedMeasures {
    Real rv;       // Realized Variance
    Real rvol;     // Realized Volatility = sqrt(RV)
    Real rq;       // Realized Quarticity
    Real bpv;      // Bipower Variation
    Real rsv_pos;  // Realized Semivariance (positive)
    Real rsv_neg;  // Realized Semivariance (negative)
    Size n_obs;    // 观测数
};

class RealizedMeasuresCalculator {
public:
    // 输入: 日内对数收益率序列 (已通过 make_returns 计算)
    static RealizedMeasures compute(const std::vector<Real>& log_returns);

    // 便捷接口: 输入价格序列, 内部计算收益率
    static RealizedMeasures compute_from_prices(const std::vector<Real>& prices);

    // 多资产 RV 协方差矩阵 (对应 R rCov)
    // 输入: 每个资产一列收益率
    static Matrix realized_covariance(const std::vector<std::vector<Real>>& returns_matrix);
};

} // namespace cpphub::hfecon
```

**R 对照**:
- `rRVar(rdata)` → `RealizedMeasures::compute(returns).rv`
- `rRealizedVolatility(rdata)` → `.rvol`
- `rQuar(rdata)` → `.rq`
- `rBPCov(rdata, makeReturns=TRUE)` → `.bpv`
- `rSV(rdata)` → `.rsv_pos, .rsv_neg`
- `rCov(rdata, makeReturns=TRUE)` → `realized_covariance(...)`

### 3.3 项 3: BNS 跳跃检验

**文件**: `include/cpphub/hfecon/tests/bns_jump_test.hpp`

**数学定义** (严格来自 [BN-S 2006], [BKS 2022] §4):

**跳跃检验统计量**:
$$Z_{BNS} = \frac{RV - BPV}{\sqrt{\frac{1}{n} \left( \frac{\pi}{2} + \frac{\pi^2}{4} \right) \max\left(1, \frac{RQ}{BPV^2}\right)}}$$

等价形式 (R `BNSjumpTest` 默认, IVestimator="RV", IQVestimator="RQ"):
$$Z_{BNS} = \frac{\sqrt{n} (RV - BPV)}{\sqrt{\vartheta_{BNS}}}, \quad \vartheta_{BNS} = \left( \frac{\pi}{2} + \frac{\pi^2}{4} \right) RQ$$

其中:
- $RV$ = Realized Variance (项 2)
- $BPV$ = Bipower Variation (项 2)
- $RQ$ = Realized Quarticity (项 2)
- $n$ = 观测数

**假设检验**:
- $H_0$: 无跳跃 (价格路径连续)
- $H_1$: 存在跳跃
- $Z_{BNS} \xrightarrow{d} \mathcal{N}(0, 1)$ 在 $H_0$ 下
- 拒绝域: $Z_{BNS} > z_{1-\alpha}$ (单尾)

**跳跃贡献比**:
$$J_{ratio} = \frac{RV - BPV}{RV}$$

**R 函数签名** (来自 BKS 2022 文档):
```r
BNSjumpTest(rdata, IVestimator = "RV", IQVestimator = "RQ",
            type = "linear", linearTransformation = "Standard",
            alpha = 0.05, makeReturns = FALSE, ...)
```

**接口设计**:
```cpp
namespace cpphub::hfecon {

enum class IVEstimator { RV, BPV };           // IVestimator 参数
enum class IQVEstimator { RQ, BQV };          // IQVestimator 参数

struct BNSJumpTestResult {
    Real z_statistic;       // Z_{BNS} 检验统计量
    Real p_value;           // 单尾 p-value
    Real critical_value;    // z_{1-alpha} 临界值
    bool reject_null;       // 是否拒绝 H_0 (存在跳跃)
    Real jump_ratio;        // (RV - BPV) / RV, 跳跃贡献比
    Real rv;
    Real bpv;
    Real rq;
    Size n_obs;
};

class BNSJumpTest {
public:
    static BNSJumpTestResult test(
        const std::vector<Real>& log_returns,
        IVEstimator iv_est = IVEstimator::RV,
        IQVEstimator iqvest = IQVEstimator::RQ,
        Real alpha = 0.05);

    // 便捷接口: 输入价格序列
    static BNSJumpTestResult test_from_prices(
        const std::vector<Real>& prices,
        IVEstimator iv_est = IVEstimator::RV,
        IQVEstimator iqvest = IQVEstimator::RQ,
        Real alpha = 0.05);
};

} // namespace cpphub::hfecon
```

**R 对照**: `BNSjumpTest(rdata, IVestimator="RV", IQVestimator="RQ")` → `BNSJumpTest::test(...)`

### 3.4 项 4: 集成测试

**文件**: `tests/unit/hfecon/test_realized_measures.cpp`

**测试矩阵**:

| 测试组 | 测试数 | 验证内容 | 对照基准 |
|---|---|---|---|
| TAQ reader | 3 | CSV 读取 / 时间聚合 / make_returns | R `aggregatePrice` + `makeReturns` |
| RV/RVol/RQ | 4 | 常数序列 / 随机序列 / 与 R 对照 / 边界 | R `rRVar` / `rRealizedVolatility` / `rQuar` |
| BPV | 2 | 纯连续路径 / 含跳跃序列 | R `rBPCov` |
| RSV | 2 | 正负分离 / 与 RV 关系 | R `rSV` |
| BNS 跳跃检验 | 3 | 无跳跃 (不拒绝) / 有跳跃 (拒绝) / 与 R 对照 | R `BNSjumpTest` |
| 多资产 RV | 1 | 2 资产协方差矩阵 | R `rCov` |

**容差** (来自 project_memory 硬约束):
- 默认: `EXPECT_NEAR(cpp_val, r_val, 1e-10)`
- 严格: `EXPECT_NEAR(cpp_val, r_val, 1e-12)` (无微结构噪声的理想情形)

**R 基准生成脚本**: `tests/fixtures/hfe/generate_r_baselines.R`
```r
library(highfrequency)
data(sampleTData)
# 生成基准 JSON, 供 C++ 测试读取
source("generate_baselines.R")
```

---

## 4. 第二波 (v1.4.1) 交付项 — Realized Kernel (微结构噪声稳健 RV)

> **范围修正 (2026-08-02 严格 review)**: 原 spec 列出的 `sparseSampling`/`noiseBPM`/`noiseAC`/`optimFrequ` 在 highfrequency 1.0.3 中**不存在** (verify_v141_functions2.R/3.R 实测). 仅 `rKernelCov` 真实存在且可用. 稀疏抽样功能已由 v1.4.0 `aggregate_price(alignBy, alignPeriod)` 覆盖. 本波范围严格限定为 Realized Kernel 单资产 RV.

**目标**: 实现微结构噪声稳健的 Realized Kernel 估计量 (BNS 2008 ECTA), 对标 R `rKernelCov` 单资产模式

**预计新增测试**: 14 个 (11 核函数 + 3 R baseline 对标), 总数 1286 → 1300

### 4.1 文献基础 (零幻觉, 全部 DOI 可溯源)

| 编号 | 文献 | DOI / 来源 | 用途 |
|---|---|---|---|
| [BNS 2008] | Barndorff-Nielsen, Hansen, Lunde, Shephard, "Designing realized kernels to measure the ex post variation of equity prices in the presence of noise" | *Econometrica* 76(6), 1481-1536, doi:10.1111/j.1468-0262.2008.00837.x | Realized Kernel 定义, 11 核函数, bandwidth 选择 |
| [BNS 2011] | Barndorff-Nielsen, Hansen, Lunde, Shephard, "Multivariate realised kernels: consistent positive semi-definite estimators of the covariation of equity prices with noise and non-synchronous trading" | *Econometrica* 79(4), 1289-1314, doi:10.3982/ECTA8119 | 多资产 Kernel Cov (推迟 v1.4.2) |
| [H-L 2006] | Hansen & Lunde, "Realized Variance and Market Microstructure Noise" | *JBES* 24(2), 127-161, doi:10.1198/073500106000000072 | 噪声方差 ω² 估计 §3 |

### 4.2 项 1: 核函数库 (`measures/kernels.hpp`)

**文件**: `include/cpphub/hfecon/measures/kernels.hpp`

**R 对照**: `listAvailableKernels()` — 实测返回 11 个核 (verify_v141_functions3.R 2026-08-02)

| C++ 枚举 | R 字符串 | k(x) 定义 (highfrequency 1.0.3 源码 `KK()` 实测, realizedMeasures.cpp L16-74) | 支撑 |
|---|---|---|---|
| `Rectangular` | `"rectangular"` | `1` | 全部 (不归零) |
| `Bartlett` | `"Bartlett"` | `1 - x` | x ∈ [0,1] |
| `Second` | `"Second"` | `1 - 2x³` | x ∈ [0,1] |
| `Epanechnikov` | `"Epanechnikov"` | `1 - x²` | x ∈ [0,1] |
| `Cubic` | `"Cubic"` | `1 - 3x² + 2x³` | x ∈ [0,1] |
| `Fifth` | `"Fifth"` | `1 - 10x³ + 15x⁴ - 6x⁵` | x ∈ [0,1] |
| `Sixth` | `"Sixth"` | `1 - 15x⁴ + 24x⁵ - 10x⁶` | x ∈ [0,1] |
| `Seventh` | `"Seventh"` | `1 - 21x⁵ + 35x⁶ - 15x⁷` | x ∈ [0,1] |
| `Eighth` | `"Eighth"` | `1 - 28x⁶ + 48x⁷ - 21x⁸` | x ∈ [0,1] |
| `Parzen` | `"Parzen"` | `1 - 6x² + 6x³` (x ≤ 0.5); `2(1-x)³` (x > 0.5) | x ∈ [0,1] |
| `TukeyHanning` | `"TukeyHanning"` | `(1 + sin(π/2 - πx))/2` | x ∈ [0,1] |
| `ModifiedTukeyHanning` | `"ModifiedTukeyHanning"` | `(1 - sin(π/2 - π(1-x)²))/2` | x ∈ [0,1] |

> **公式来源 (零幻觉)**: highfrequency 1.0.3 CRAN 源码 `realizedMeasures.cpp` 函数 `KK(double x, int type)`, 通过 `download.file("https://cran.r-project.org/src/contrib/highfrequency_1.0.3.tar.gz")` 下载并解压获取 (2026-08-02 实测).
>
> **与 BNS 2008 论文的差异 (重要)**:
> 1. `Second` 核: R 实现 `1 - 2x³`, BNS 2008 Table 1 为 `1 - x²` — **R 实现与论文不符**
> 2. `Seventh`/`Eighth` 核: R 实现多项式阶数与系数均与 BNS 2008 Table 1 不同 — **R 实现与论文不符**
> 3. `Parzen` 核: R 实现是分段多项式 (Bartlett-Parzen 形式), 与 BNS 2008 §3 一致
> 4. `TukeyHanning`/`ModifiedTukeyHanning`: R 用 `sin(π/2 - πx)` 形式, 与 BNS 2008 `(1+cos(πx))/2` 等价 (sin(π/2 - θ) = cos(θ))
>
> **决策**: C++ 实现严格对标 R 源码 (而非 BNS 论文), 以保证 R baseline 数值一致. 差异记录为 "discovery" (R 实现偏离论文).

**接口签名**:

```cpp
enum class KernelType {
    Rectangular, Bartlett, Second, Epanechnikov, Cubic,
    Fifth, Sixth, Seventh, Eighth, Parzen, TukeyHanning, ModifiedTukeyHanning
};

// 核函数值 k(x), 支持 |x| > 1 时返回 0
// 异常: 未知 KernelType 抛 invalid_argument
Real kernel_value(KernelType type, Real x) noexcept(false);

// 从 R 字符串解析 KernelType (大小写敏感, 与 R listAvailableKernels() 一致)
// 异常: 未知字符串抛 invalid_argument
KernelType parse_kernel_type(const std::string& name);
```

### 4.3 项 2: 噪声方差估计 (`noise/noise_variance.hpp`)

**文件**: `include/cpphub/hfecon/noise/noise_variance.hpp`

**R 对照**: highfrequency 1.0.3 无独立导出函数, rKernelCov 内部使用 (BNS 2008 §4.4)

**算法 (BNS 2008 eq. 40, H-L 2006 §3)**:

```
ω² = (1/(2n)) * Σ_{i=1}^{n} r_i²     (n = 收益率数)
```

**接口签名**:

```cpp
struct NoiseVarianceResult {
    Real omega2;    // 噪声方差估计 ω²
    Real integrated_variance;  // IV 估计 (RV - ω²)
    Size n_obs;
};

class NoiseVarianceEstimator {
public:
    // 输入: 日内对数收益率序列
    // 异常: n < 2 抛 invalid_argument
    static NoiseVarianceResult estimate(const std::vector<Real>& log_returns);
};
```

### 4.4 项 3: 最优 Bandwidth 选择 (`noise/bandwidth.hpp`)

**文件**: `include/cpphub/hfecon/noise/bandwidth.hpp`

**R 对照**: rKernelCov 内部 `bandwidth` 选择 (BNS 2008 §4.5)

**算法 (BNS 2008 eq. 51)**:

```
H* = c × ξ^(4/5) × (ω²/IV)^(2/5) × n^(3/5)

其中:
  c = 5.74 (Bartlett 核的最优常数, BNS 2008 Table 4)
  ξ² = IV + 2ω²  (或更稳健的估计)
  ω² = 噪声方差 (项 4.3)
  IV  = RV - ω²  (积分方差估计)
  n   = 观测数

round 到最近整数, 下界 H ≥ 1
```

**接口签名**:

```cpp
// 输入: 噪声方差 ω², 积分方差 IV, 观测数 n, 核类型
// 返回: 最优 bandwidth H (整数, ≥ 1)
// 异常: omega2 ≤ 0 或 IV ≤ 0 或 n == 0 抛 invalid_argument
Size optimal_bandwidth(Real omega2, Real integrated_variance,
                       Size n_obs, KernelType type);
```

### 4.5 项 4: Realized Kernel 主估计量 (`measures/realized_kernel.hpp`)

**文件**: `include/cpphub/hfecon/measures/realized_kernel.hpp`

**R 对照**: `rKernelCov(rData, cor=FALSE, kernelType="Bartlett", kernelParam=1, kernelDOFadj=TRUE)`

**R 实测签名 (verify_v141_functions3.R)**:
```
rKernelCov(rData, cor = FALSE, alignBy = NULL, alignPeriod = NULL,
           makeReturns = FALSE, kernelType = "rectangular",
           kernelParam = 1, kernelDOFadj = TRUE, ...)
```

**算法 (highfrequency 1.0.3 源码 `kernelEstimator()` 实测, realizedMeasures.cpp L77-111)**:

```
// 输入: r (长度 n), bandwidth H (= kernelParam), 核类型 type, DOF 调整 adj
// 注意: R 实现与 BNS 2008 论文有 2 处关键差异

nab = n - 1

// Step 1: 计算自协方差 ab[h], ab2[h] (h = 0..H)
ab[h]  = Σ_{i=0}^{n-1-h} r[i] * r[i+h]    // forward lag h (γ_h)
ab2[h] = Σ_{i=h}^{n-1}   r[i] * r[i-h]    // backward lag h (γ_{-h})

// Step 2: 加权求和
ans = 0
for h = 0 to H:
    if h == 0:
        w = 1.0                              // γ_0 权重恒为 1
    else:
        w = KK((h-1)/H, type)                // 关键: (h-1)/H, 不是 h/H (半整数偏移)

    if adj == 0:
        theadj = 1.0
    else:
        theadj = n / (n - h)                 // 关键: 逐 lag 调整, 不是整体 n/(n-H)

    if h == 0:
        ans += w * theadj * ab[0]            // γ_0
    else:
        ans += w * (theadj*ab[h] + theadj*ab2[h])  // γ_h + γ_{-h}

return ans
```

**与 BNS 2008 论文的关键差异 (重要)**:
1. **权重偏移**: R 用 `w = KK((h-1)/H)`, BNS 论文用 `w = KK(h/H)`. 后果: R 实现中 h=1 时 w=KK(0)=1 (所有核 k(0)=1), 而 BNS 论文中 h=1 时 w=KK(1/H)≠1.
2. **DOF 调整**: R 用逐 lag 调整 `theadj = n/(n-h)`, BNS 论文用整体调整 `n/(n-H)`.
3. **决策**: C++ 实现严格对标 R 源码, 以保证 R baseline 数值一致 (容差 1e-12).

**接口签名**:

```cpp
struct RealizedKernelResult {
    Real rk;              // Realized Kernel 估计值 (DOF 调整后, 若启用)
    Real rv;              // γ_0 (Realized Variance, 未调整)
    Real gamma_1;         // γ_1 (一阶自协方差, 用于噪声诊断)
    Size bandwidth;       // 实际使用的 H (= kernelParam)
    Size n_obs;
    KernelType kernel;
    bool dof_adjusted;    // 是否应用 DOF 调整
};

class RealizedKernel {
public:
    // 主接口: 严格对标 R rKernelCov 单资产模式
    // 输入: 日内对数收益率序列 (R rData)
    // 参数:
    //   kernel         - 核类型 (默认 Rectangular, 与 R 默认 kernelType="rectangular" 一致)
    //   kernel_param   - bandwidth H (默认 1, 与 R 默认 kernelParam=1 一致)
    //   kernel_dof_adj - 是否应用 DOF 调整 (默认 true, 与 R 默认 kernelDOFadj=TRUE 一致)
    // 异常: n < kernel_param + 1 抛 invalid_argument (需足够观测计算 ab[H])
    //       kernel_param == 0 抛 invalid_argument
    static RealizedKernelResult estimate(
        const std::vector<Real>& log_returns,
        KernelType kernel = KernelType::Rectangular,
        Size kernel_param = 1,
        bool kernel_dof_adj = true);

    // 便捷接口: 输入价格序列 (内部 make_returns)
    static RealizedKernelResult estimate_from_prices(
        const std::vector<Real>& prices,
        KernelType kernel = KernelType::Rectangular,
        Size kernel_param = 1,
        bool kernel_dof_adj = true);
};
```

> **设计决策 (基于 R 源码实测)**:
> 1. **不自动选择 bandwidth**: R `rKernelCov` 接受用户提供的 `kernelParam`, 不内部计算最优 H. C++ 实现保持一致, `kernel_param` 为必填参数 (默认 1).
> 2. **默认 kernelType=Rectangular**: 与 R 默认一致 (`kernelType="rectangular"`). v1.4.0 spec 误写为 Bartlett, 已修正.
> 3. **噪声方差与 bandwidth 选择 (项 4.3/4.4)**: 作为 C++ 扩展工具保留 (BNS 2008 公式), 标注"非 R 对标", 供高级用户使用, 但 `RealizedKernel::estimate` 不调用它们.

### 4.6 项 5: 测试矩阵 (`tests/unit/hfecon/test_realized_kernel.cpp`)

**预计测试**: 14 个

**A. 核函数单测 (11 个, 对照 R listAvailableKernels + kfunc 源码)**:

| # | 测试名 | 描述 | 容差 |
|---|---|---|---|
| 1 | `Kernel.Rectangular` | k(0)=1, k(0.5)=1, k(1)=1, k(1.5)=0 | 1e-15 |
| 2 | `Kernel.Bartlett` | k(0)=1, k(0.5)=0.5, k(1)=0, k(1.5)=0 | 1e-15 |
| 3 | `Kernel.Second` | k(0)=1, k(0.5)=0.75, k(1)=0 | 1e-15 |
| 4 | `Kernel.Epanechnikov` | k(0)=1, k(0.5)=0.75, k(1)=0 | 1e-15 |
| 5 | `Kernel.Cubic` | k(0)=1, k(0.5)=0.5, k(1)=0 | 1e-15 |
| 6 | `Kernel.Fifth` | k(0)=1, k(0.5)=?, k(1)=0 | 1e-15 |
| 7 | `Kernel.Sixth` | k(0)=1, k(1)=0 | 1e-15 |
| 8 | `Kernel.Seventh` | k(0)=1, k(1)=0 | 1e-15 |
| 9 | `Kernel.Eighth` | k(0)=1, k(1)=0 | 1e-15 |
| 10 | `Kernel.Parzen` | k(0)=1, k(0.25)=?, k(0.5)=0.5, k(1)=0 | 1e-15 |
| 11 | `Kernel.TukeyHanning` | k(0)=1, k(0.5)=0.5, k(1)=0 | 1e-15 |

**B. R baseline 对标 (3 个, 硬编码 R 输出值)**:

| # | 测试名 | 场景 | 对照 | 容差 |
|---|---|---|---|---|
| 12 | `RBaselineExact.KernelBartlett` | R `rKernelCov(ret, kernelType="Bartlett", kernelParam=1)` 单资产 | R 输出值 | 1e-12 |
| 13 | `RBaselineExact.KernelParzen` | R `rKernelCov(ret, kernelType="Parzen", kernelParam=1)` 单资产 | R 输出值 | 1e-12 |
| 14 | `RBaselineExact.KernelNoiseRejection` | 含噪声合成数据, RK 应显著小于 RV (噪声稳健性) | 相对差 > 5% | - |

**R baseline 生成脚本**: `tests/fixtures/hfe/generate_v141_baselines.R`

---

## 5. 第三波 (v1.4.2) 交付项 — 多资产 Cov 估计 + HAR/HEAVY 预测模型

> **范围调整 (2026-08-02)**: 基于 verify_v141_functions3.R 实测, highfrequency 1.0.3 提供多资产噪声稳健协方差估计方法 (rAVGCov/rTSCov/rMRCov/rHYCov 等) 与 HAR/HEAVY 预测模型. 原计划中的 "ARFIMA 估计" 与 "预测精度评估" highfrequency 无对应函数, 推迟或剔除.
>
> **R 源码调研 (2026-08-02)**: dump 7 个函数完整源码 (`tests/fixtures/hfe/v142_source_dump.txt`, 74KB), 识别 9 项 R 源码 vs 论文差异 (幻觉源). 严格对标 R 实现, spec 显式标注差异.

**目标**: 多资产噪声稳健协方差 + HAR/HEAVY RV 预测模型

**预计新增测试**: ~31 个, 总数 1300 → ~1331

**实施波次** (基于依赖关系):
- **波次 A** (基础设施): `refreshTimeMatching` + `preAveraging` 内部工具 + rHYCov (最简单)
- **波次 B** (Cov 估计器): rAVGCov + rTSCov + rRTSCov + rMRCov
- **波次 C** (预测模型): HARmodel + HEAVYmodel

### 5.1 多资产 Cov 估计 (5 个方法, 全部 R 对标)

| 项 | 文件 | R 对照 (实测签名) | 文献 |
|---|---|---|---|
| Hayashi-Yoshida | `measures/hayashi_yoshida_cov.hpp` | `rHYCov(rData, cor=FALSE, period=1, alignBy="seconds", alignPeriod=1, makeReturns=FALSE, makePsd=TRUE)` | Hayashi & Yoshida (2005), *J. Financial Econometrics* 3(4), doi:10.1093/jjfinec/nbi013 |
| Pre-averaging Cov | `measures/preaveraged_cov.hpp` | `rAVGCov(rData, cor=FALSE, alignBy="minutes", alignPeriod=5, k=1, makeReturns=FALSE)` | Jacaud, Li, Mykland, Podolskij, Vetter (2009), *AOS* 37(1), 280-318, doi:10.1214/07-AOS568 |
| Two-scale Cov | `measures/two_scale_cov.hpp` | `rTSCov(pData, cor=FALSE, K=300, J=1, KCov=NULL, JCov=NULL, KVar=NULL, JVar=NULL, makePsd=FALSE)` | Zhang, Mykland, Aït-Sahalia (2005), *JASA* 100(472), 1394-1411, doi:10.1198/016214505000000548 |
| Robust Two-scale | `measures/robust_two_scale_cov.hpp` | `rRTSCov(pData, cor=FALSE, startIV=NULL, noisevar=NULL, K=300, J=1, ..., eta=9, makePsd=FALSE)` | Zhang (2011), *JASA* 106(495), doi:10.1198/jasa.2011.tm10384 |
| Modulated RC | `measures/modulated_realized_cov.hpp` | `rMRCov(pData, pairwise=FALSE, makePsd=FALSE, theta=0.8, crossAssetNoiseCorrection=FALSE)` | Christensen, Podolskij, Vetter (2013), *J. Econometrics* 173(1), doi:10.1016/j.jeconom.2012.08.016 |

### 5.2 HAR/HEAVY 预测模型

| 项 | 文件 | R 对照 (实测签名) | 文献 |
|---|---|---|---|
| HAR 模型 | `models/har_model.hpp` | `HARmodel(data, periods=c(1,5,22), periodsJ=c(1,5,22), periodsQ=c(1), leverage=NULL, RVest=c("rCov","rBPCov","rQuar"), type="HAR", inputType="RM", jumpTest="ABDJumptest", alpha=0.05, h=1, transform=NULL, externalRegressor=NULL, periodsExternal=c(1))` | Corsi (2009), *JFE* 4(2), 174-196, doi:10.1093/jjfinec/nbp001 |
| HEAVY 模型 | `models/heavy_model.hpp` | `HEAVYmodel(data, startingValues=NULL)` | Shephard & Sheppard (2010), *Restat* 92(2), doi:10.1162/REST_a_00017; Noureldin, Shephard, Sheppard (2012), *JAE* 27(8), doi:10.1002/jae.1260 |

### 5.3 R 源码 vs 论文差异 (排幻觉, 2026-08-02 实测)

| 编号 | 函数 | 差异 | 决策 |
|---|---|---|---|
| D1 | rTSCov | 对角线用 TSRV 公式 (`n/(n-K+1)·[Y,Y]^(K) - n/(n-J+1)·[Y,Y]^(J)`), 非对角线用 TSCov_bi 公式 (`[Y,Y]^(K) - [Y,Y]^(J)` 无系数) | C++ 分别实现两个公式 |
| D2 | rAVGCov | 单资产加 `(m+1)/m` 系数 (m = alignPeriod/k), 多资产不加 | 单/多资产分支处理 |
| D3 | rMRCov | 三分支: `crv` (单资产) / `preavbi` (双资产) / 矩阵版, 公式不兼容; `psi2kn` vs `psi2=1/12` 不同 | 按 R 三分支实现 |
| D4 | rHYCov | 用整数索引而非时间戳 (`pcovcc` 退化版), `period` 参数仅用于对齐 | C++ 实现整数索引版 |
| D5 | rHYCov | 默认 `makePsd=TRUE` (其他 Cov 函数 FALSE), 强制 PSD 投影 | 默认参数对齐 R |
| D6 | rRTSCov | `eta=9` 硬编码 `ccc=1.0415`, 其余 eta 查 30 项表 | 实现 eta=9 快速路径 + 查表 |
| D7 | HARmodel | RQ 变换用 BPQ 2016 (`sqrt(RQ) - sqrt(mean(RM3))`), 非 Corsi 2009 原始 RQ | 实现 BPQ 变换 |
| D8 | HEAVYmodel | 强制去均值 (`data - mean(data)`), 无常数项; 用 BFGS 优化 MLE | C++ 实现 BFGS + 去均值 |
| D9 | HARmodel | `har_agg` 索引边界 `[j-p, j-1]` 含 j-1 不含 j (R 1-based → C++ 0-based 需调整) | 仔细处理索引 |

### 5.4 可移植 C++ 源码 (来自 highfrequency 1.0.3 src/)

| 函数 | 位置 | 用途 |
|---|---|---|
| `refreshTimeMatching` | internals.cpp L23-72 | 4 个 Cov 函数共用的非同步时间匹配 |
| `preAveragingReturnsInternal` | internals.cpp L83-106 | rMRCov 专属 pre-averaging |
| `pcovcc` | realizedMeasures.cpp L130-157 | rHYCov 专属 (整数索引版) |
| `har_agg` | HARmodel.cpp L7-21 | HARmodel 滚动窗口聚合 |
| `calcRecVarEq` | HEAVYmodel.cpp L5-15 | HEAVYmodel 方差递归 |

### 5.5 测试矩阵 (31 个新测试)

| 模块 | 测试数 | 场景 |
|---|---|---|
| rHYCov | 5 | 单资产 (常数序列) + 双资产 (已知 ρ=0.5) + makePsd + period 参数 + 异常处理 |
| rAVGCov | 5 | 单资产 + 双资产 + k 参数 + alignBy 参数 + 异常处理 |
| rTSCov | 5 | 单资产 + 双资产 + K/J 参数 + makePsd + 异常处理 |
| rRTSCov | 4 | 单资产 + 双资产 + noisevar + eta 参数 |
| rMRCov | 4 | 单资产 + 双资产 + theta 参数 + 异常处理 |
| HARmodel | 4 | 基础 HAR + type="HARJ" + type="CHAR" + transform |
| HEAVYmodel | 4 | 基础估计 + startingValues + 残差诊断 + 异常处理 |

---

## 6. 第四波 (v1.4.3) 交付项 — 流动性度量 + 高级跳跃检验

> **范围调整 (2026-08-03 R 源码深度调研后)**: 
> - 基于 `liquidityMeasures.R`/`jumpTests.R`/`internalJumpTests.R`/`internals.cpp`/`dataHandling.R`/`spotVolAndDrift.R` 实测源码
> - `intradayJumpTest` 依赖 `spotVol` (7 种估计器) + `spotDrift`, 完整实现工作量过大. v1.4.3 仅实现 `volEstimator="RM"` 模式 (滚动窗口 rBPCov), `driftEstimator="none"`. PARM 模式推迟.
> - `rankJumpTest` 需多资产输入 (marketPrice + stockPrices list), 含 SVD + bootstrap, 完整实现.
> - 预计新增测试: ~25 个, 总数 1362 → ~1387

**目标**: 流动性度量 (3 文件) + 高级跳跃检验 (4 文件)

### 6.0 依赖分析与基础设施

| 新文件 | 依赖 (已有) | R 对照源码 |
|---|---|---|
| `liquidity/spread_cleaner.hpp` | core/ | dataHandling.R L1617/1670/3019 |
| `liquidity/liquidity_measures.hpp` | core/ | liquidityMeasures.R L231/346 |
| `liquidity/amihud.hpp` | core/ | (无 R 对照) |
| `tests/aj_jump_test.hpp` | measures/realized_measures | jumpTests.R L106 + internalJumpTests.R |
| `tests/jo_jump_test.hpp` | measures/realized_measures | jumpTests.R L446 + internals.cpp L207 |
| `tests/intraday_jump_test.hpp` | measures/realized_measures, data/ | jumpTests.R L583 + spotVolAndDrift.R L656 |
| `tests/rank_jump_test.hpp` | core/ (Jacobi SVD) | jumpTests.R L976 + internalJumpTests.R L149 |

### 6.1 流动性度量

#### 6.1.1 价差清洗工具 (`liquidity/spread_cleaner.hpp`)

**R 对照**: `rmLargeSpread(qData, maxi=50)` / `rmNegativeSpread(qData)` / `spreadPrices(data)`

**C++ 接口**:
```cpp
namespace cpphub::v1::hfecon {

// 每日计算 SPREAD = OFR - BID 中位数, 保留 SPREAD < SPREAD_MEDIAN * maxi
template<typename QuoteContainer>
QuoteContainer rm_large_spread(const QuoteContainer& qData, double maxi = 50.0);

// 保留 OFR > BID (严格大于)
template<typename QuoteContainer>
QuoteContainer rm_negative_spread(const QuoteContainer& qData);

// 长格式 (DT, SYMBOL, PRICE) → 宽格式 (DT, sym1, sym2, ...)
std::vector<std::vector<double>> spread_prices(
    const std::vector<std::chrono::system_clock::time_point>& dt,
    const std::vector<std::string>& symbols,
    const std::vector<double>& prices);

} // namespace
```

**算法 (R 源码对标)**:
- `rmLargeSpread`: 按日期分组 → 每日 SPREAD_MEDIAN = median(OFR-BID) → 保留 SPREAD < SPREAD_MEDIAN * maxi
- `rmNegativeSpread`: OFR > BID
- `spreadPrices`: split by SYMBOL → outer join on DT

#### 6.1.2 综合流动性度量 (`liquidity/liquidity_measures.hpp`)

**R 对照**: `getLiquidityMeasures(tqData, win=300)` — 23 种度量

**C++ 接口**:
```cpp
namespace cpphub::v1::hfecon {

struct LiquidityMeasures {
    std::vector<double> effectiveSpread, realizedSpread, valueTrade, signedValueTrade;
    std::vector<double> depthImbalanceDifference, depthImbalanceRatio;
    std::vector<double> proportionalEffectiveSpread, proportionalRealizedSpread;
    std::vector<double> priceImpact, proportionalPriceImpact;
    std::vector<double> halfTradedSpread, proportionalHalfTradedSpread;
    std::vector<double> squaredLogReturn, absLogReturn;
    std::vector<double> quotedSpread, proportionalQuotedSpread;
    std::vector<double> logQuotedSpread, logQuotedSize, quotedSlope, logQSlope;
    std::vector<double> midQuoteSquaredReturn, midQuoteAbsReturn, signedTradeSize;
};

LiquidityMeasures get_liquidity_measures(
    const std::vector<double>& price, const std::vector<double>& bid,
    const std::vector<double>& ofr, const std::vector<double>& size,
    const std::vector<double>& ofrsiz, const std::vector<double>& bidsiz,
    const std::optional<std::vector<int>>& direction = std::nullopt,
    int win = 300);

// Lee-Ready 交易方向推断: 返回 1 (buy) 或 -1 (sell)
std::vector<int> get_trade_direction(
    const std::vector<double>& price,
    const std::vector<double>& bid,
    const std::vector<double>& ofr);

} // namespace
```

**算法 (R 源码对标, 排幻觉)**:

`get_trade_direction` (**排幻觉 D1** — tick rule + midpoint 混合, 非纯 Lee-Ready):
1. `midpoints = (bid + ofr) / 2`
2. `rets = diff(price)`, 首元素 = 0 (R: `c(TRUE, ...)`)
3. tick rule: rets > 0 → 1, rets < 0 → -1, rets == 0 → NA + locf
4. **midpoint 覆盖**: price < mid → -1, price > mid → 1, price == mid → 保留 tick rule
5. 首观测 = 1 (buy)

`get_liquidity_measures`:
- effectiveSpread = `2 * direction * (PRICE - mid)`
- realizedSpread = `2 * direction * (PRICE - mid[t+win])` (**排幻觉 D2**: lead shift, 越界为 NaN)
- depthImbalanceRatio = `(direction * OFRSIZ / BIDSIZ) ^ direction` (**排幻觉 D3**: direction 在底数和指数)
- 其余 20 种按定义直接计算

#### 6.1.3 Amihud 流动性 (`liquidity/amihud.hpp`)

**R 对照**: 无 (highfrequency 无直接实现)

```cpp
namespace cpphub::v1::hfecon {
// ILLIQ_t = (1/N) * sum |r_t| / DVOL_t
double amihud_illiquidity(
    const std::vector<double>& dailyReturns,
    const std::vector<double>& dailyDollarVolume);
} // namespace
```

### 6.2 高级跳跃检验

#### 6.2.1 AJ 跳跃检验 (`tests/aj_jump_test.hpp`)

**R 对照**: `AJjumpTest(pData, p=4, k=2, alignBy, alignPeriod, alphaMultiplier=4, alpha=0.975)`
**文献**: Aït-Sahalia & Jacod (2009), *Annals of Statistics* 37(1), 184-222

```cpp
namespace cpphub::v1::hfecon {
struct AJJumpTestResult { double ztest, criticalLower, criticalUpper, pvalue; };
AJJumpTestResult aj_jump_test(
    const std::vector<double>& pData,
    int p = 4, int k = 2,
    const std::string& alignBy = "seconds", int alignPeriod = 1,
    double alphaMultiplier = 4.0, double alpha = 0.975);
} // namespace
```

**算法 (排幻觉 D4-D9)**:
1. **动态 alpha** (D4): `alpha = alphaMultiplier * sqrt(rCov(pData))` — 非固定, 与 RV 平方根成正比
2. `N = length(pData) - 1`, `w = 0.47`, `cvalue = alpha * (1/N)^w`
3. `h = alignPeriod * scale(alignBy)` (scale: seconds=1, minutes=60, hours=3600)
4. **整数抽样** (D5): `seq1 = seq(1, N, h)`, `seq2 = seq(1, N, h*k)`
5. `r = |makeReturns(pData)|`, pv1 = sum(r[seq1]^p), pv2 = sum(r[seq2]^p), S = pv2/pv1
6. **selection 筛选** (D6): `rse = r[|r| < cvalue]` — 只用小收益率
7. `V = calculateV(rse, p, k, N)`, `AJtest = (S - k^(p/2-1)) / sqrt(V)`

`calculateV` (D7): `Ap = (1/N)^(1-p/2)/mup * sum(rse^p)`, `A2p = (1/N)^(1-p)/mu2p * sum(rse^(2p))`, `V = calculateNpk(p,k) * A2p / (N * Ap^2)`

`calculateNpk` (D8): `npk = (1/mup^2) * (k^(p-2)*(1+k)*mu2p + k^(p-2)*(k-1)*mup^2 - 2*k^(p/2-1)*fmupk(p,k))`

`fmupk` 查表 (D9 — R 硬编码, 非论文公式):
| p\k | 2 | 3 | 4 |
|---|---|---|---|
| 2 | 4.00 | 5.00 | 6.00 |
| 3 | 24.07 | 33.63 | 43.74 |
| 4 | 204.04 | 320.26 | 455.67 |

其他 (p,k): 蒙特卡洛 `mukp(p,k,t=1e6)`, 100 rep 取均值 round 2 位

#### 6.2.2 JO 跳跃检验 (`tests/jo_jump_test.hpp`)

**R 对照**: `JOjumpTest(pData, power=4, alignBy, alignPeriod, alpha=0.975)`
**文献**: Jiang & Oomen (2008), *Mathematical Finance* 18(3), doi:10.1111/j.1467-9965.2008.00343.x

```cpp
namespace cpphub::v1::hfecon {
struct JOJumpTestResult { double ztest, criticalLower, criticalUpper, pvalue; };
JOJumpTestResult jo_jump_test(
    const std::vector<double>& pData,
    int power = 4,
    const std::string& alignBy = "seconds", int alignPeriod = 1,
    double alpha = 0.975);
} // namespace
```

**算法 (排幻觉 D10-D13)**:
1. `R = simre(pData)` 简单收益率 (D10): `R[i] = P[i]/P[i-1] - 1`, R[0]=0
2. `r = makeReturns(pData)` 对数收益率: `r[i] = log(P[i]) - log(P[i-1])`, r[0]=0
3. `N = length(pData) - 1`, `bv = RBPVar(r)` = `(pi/2)*sum(|r[0:n-1]|*|r[1:n]|)`
4. `rv = sum(r^2)`
5. **SwV = 2 * sum(R - r)** (D11): 简单与对数收益率之差
6. **mu1 = 2^3 * gamma(3.5)/gamma(0.5)** (D12): 6 阶矩 μ₆, 非 power 阶

**power=4**: `q = |rollApplyProdWrapper(r,4)|`, `mu2 = 2^(3/4)*gamma(7/4)/gamma(0.5)` (1.5 阶), `av = mu1/9 * N^3 * (mu2)^(-4) / (N-5) * sum(q^(1.5))`, `JOtest = N*bv/sqrt(av)*(1-rv/SwV)`

**power=6**: `q = |rollApplyProdWrapper(r,6)|`, `mu2 = 2^(1/2)*gamma(1)/gamma(0.5)` (1 阶), `av = mu1/9 * N^3 * (mu2)^(-6) / (N-7) * sum(q)`, `JOtest = N*bv/sqrt(av)*(1-rv/SwV)`

**rollApplyProdWrapper C++** (D13, internals.cpp L207): `m = m - 1` 后窗口 m 个元素, `out[i] = prod(x[i:i+m-1])`, 输出长度 n-m+1

#### 6.2.3 日内跳跃检验 (`tests/intraday_jump_test.hpp`)

**R 对照**: `intradayJumpTest(pData, volEstimator="RM", driftEstimator="none", alpha=0.95, ...)`
**文献**: Lee & Mykland (2008), *JFE* 6(5); Christensen, Oomen, Podolskij (2014)

```cpp
namespace cpphub::v1::hfecon {
struct IntradayJumpTestResult {
    std::vector<double> ztest, spotVol;
    double criticalValue;
    int n;
};
IntradayJumpTestResult intraday_jump_test(
    const std::vector<double>& pData,
    const std::vector<std::chrono::system_clock::time_point>& dt,
    const std::string& rmType = "rBPCov",
    int lookBackPeriod = 10,
    double alpha = 0.95,
    const std::string& alignBy = "minutes", int alignPeriod = 5,
    const std::string& marketOpen = "09:30:00",
    const std::string& marketClose = "16:00:00");
} // namespace
```

**算法 (排幻觉 D14-D16)**:
1. 聚合价格到 alignBy/alignPeriod 网格, 按日分组
2. `RETURN = log(PRICE) - log(PRICE[t-1])` (按日)
3. **spotVol RM 模式** (D14): 滚动窗口 `vol[j] = RM(returns[j-lookBack+1 : j])`, RM = rBPCov/rMinRVar/rMedRVar
4. **vol 调整** (D15): `vol = sqrt(vol^2 / (lookBackPeriod-2))`
5. drift = 0
6. `test = (return - drift) / vol`
7. **Lee-Mykland 临界值** (D16): `n = NROW(pData)` (原始观测数), `Cn = sqrt(2log(n)) - (log(pi)+log(log(n)))/(2*sqrt(2log(n)))`, `Sn = 1/sqrt(2log(n))`, `criticalValue = Cn + Sn*(-log(-log(1-alpha)))`

> v1.4.3 仅 RM 模式. PARM 模式推迟 v1.4.4.

#### 6.2.4 Rank 跳跃检验 (`tests/rank_jump_test.hpp`)

**R 对照**: `rankJumpTest(marketPrice, stockPrices, alpha=c(5,3), coarseFreq=10, localWindow=30, rank=1, BoxCox=1, quantiles, nBoot=1000, ...)`
**文献**: Bollerslev, Todorov (2011), *JFE* 9(2), doi:10.1093/jjfinec/nbr010

```cpp
namespace cpphub::v1::hfecon {
struct RankJumpTestResult {
    std::vector<double> criticalValues, testStatistic;
    std::vector<int> jumpIndices;
};
RankJumpTestResult rank_jump_test(
    const std::vector<double>& marketPrice,
    const std::vector<std::chrono::system_clock::time_point>& marketDt,
    const std::vector<std::vector<double>>& stockPrices,
    const std::vector<std::vector<std::chrono::system_clock::time_point>>& stockDts,
    std::vector<double> alpha = {5.0, 3.0},
    int coarseFreq = 10, int localWindow = 30, int rank = 1,
    std::vector<double> boxCox = {1.0},
    std::vector<double> quantiles = {0.9, 0.95, 0.99},
    int nBoot = 1000);
} // namespace
```

**算法 (排幻觉 D17-D23)**:
1. 聚合 + 对数收益率 (市场 + 个股)
2. **jumpDetection** (D17): `bpv = (pi/2)*colSums(|r[0:n-1]|*|r[1:n]|)`, `rv = colSums(r^2)`, TODadjustments (polyOrder=2), `Un = alpha*sqrt(kronecker(pmin(bpv,rv), TODfit))*(1/nRets)^0.49`, jumpIndices = which(|r| > Un)
3. **jumps 累积** (D18): `jumps = sum(stockReturns[jumpIndices+i])` for i in 0:(coarseFreq-1)
4. **SVD 全分解** (D19): `svd(jumps, nu=nrow, nv=ncol)`, U2=U[:,rank+1:], V2=V[:,rank+1:], singularValues=d[rank+1:]^2
5. **testStatistic** (D20): `sum(BoxCox__(singularValues, a))` for each BoxCox
6. **bootstrap** (D21): nBoot 次, `dxc = pmax(pmin(ret, Un), -Un)` 截断, 每次随机左右窗口 + `kappaStar = runif(1)`, `zetaStar = sqrt(kappaStar)*dxcLeft + sqrt(coarseFreq-kappaStar)*dxcRight`, `tmp = t(U2)%*%zetaStar%*%V2`, `simTestStat = sum(tmp^2)`, criticalValues = quantile(simTestStat, quantiles)
7. **BoxCox__** (D22): `lambda=0 → log(1+x)`, `lambda≠0 → ((1+x)^lambda-1)/lambda`
8. **TODadjustments** (D23): `timeOfDayScatter = 1.249531 * rowMeans(|r_i*r_{i+1}*r_{i+2}|^(2/3))`, Vandermonde OLS, 归一化均值 1

### 6.3 R 源码 vs 论文差异 (排幻觉清单, 2026-08-03 实测)

| ID | 函数 | R 源码行为 | 论文/文档 | 影响 |
|---|---|---|---|---|
| D1 | getTradeDirection | tick rule + midpoint 混合, 首观测=buy | Lee-Ready 1991 纯 tick rule | 方向推断结果不同 |
| D2 | getLiquidityMeasures | realizedSpread 用 lead shift mid[t+win] | 文档说 "t+300 秒" | 越界为 NaN |
| D3 | depthImbalanceRatio | `(direction*OFRSIZ/BIDSIZ)^direction` | 文档公式不含 direction 在底数 | direction=-1 时取倒数 |
| D4 | AJjumpTest | `alpha = alphaMultiplier*sqrt(RV)` 动态 | 论文 alpha 固定 | 阈值随波动率变化 |
| D5 | AJjumpTest | `seq(1,N,h)` 整数步长抽样 | 论文连续窗口 | 幂变差求和点不同 |
| D6 | AJjumpTest | `rse = r[|r|<cvalue]` 筛选 | 论文用全部收益率 | V 只用小收益率 |
| D7 | calculateV | `Ap = (1/N)^(1-p/2)/mup*sum(rse^p)` | 论文系数不同 | 归一化不同 |
| D8 | calculateNpk | 含 `fmupk(p,k)` 查表项 | 论文用解析公式 | (p,k) 非 (2-4,2-4) 时蒙特卡洛 |
| D9 | fmupk | 硬编码表 (p=2,3,4 × k=2,3,4) | 论文无此表 | 必须用 R 查表值 |
| D10 | JOjumpTest | R=simre 简单, r=makeReturns 对数 | 论文用同一收益率 | SwV=2*sum(R-r) |
| D11 | JOjumpTest | `SwV = 2*sum(R-r)` | 论文 SwV 定义不同 | 必须 R-r |
| D12 | JOjumpTest | `mu1 = 2^3*gamma(3.5)/gamma(0.5)` = μ₆ | 文档 "mu1" 易误解 | 6 阶矩, 非 power 阶 |
| D13 | rollApplyProdWrapper | `m = m-1` 后窗口 m 元素 | R 文档 "m 个元素乘积" | C++ 移植需注意 m-1 |
| D14 | intradayJumpTest | `vol = sqrt(vol^2/(lookBack-2))` | Lee-Mykland 原文无此调整 | RM 估计器需除以 (lookBack-2) |
| D15 | intradayJumpTest | Cn 无 sqrt(2/pi) 常数 | Lee-Mykland Eq.12 有常数 | R 去除常数使 L~N(0,1) |
| D16 | intradayJumpTest | `n = NROW(pData)` 原始观测数 | 文档说 "对齐后观测数" | 临界值用原始 n |
| D17 | jumpDetection | `Un = alpha*sqrt(kronecker(pmin(bpv,rv),TODfit))*(1/nRets)^0.49` | 论文无 TOD 调整 | 日内模式修正 |
| D18 | rankJumpTest | `jumps = sum(ret[jumpIdx+i])` i=0..coarseFreq-1 | 论文粗采样定义不同 | 累积窗口 |
| D19 | rankJumpTest | `svd(jumps, nu=nrow, nv=ncol)` 全 SVD | 标准 SVD 即可 | 需全分解取 U2/V2 |
| D20 | rankJumpTest | `testStat = sum(BoxCox__(d^2, a))` | 论文无 BoxCox | R 添加 BoxCox 变换 |
| D21 | rankJumpTest | `dxc = pmax(pmin(ret, Un), -Un)` 截断 | 论文无截断 | bootstrap 用截断收益 |
| D22 | BoxCox__ | `lambda=0 → log(1+x)` | 标准 BoxCox `log(x)` | R 用 1+x 避免 log(0) |
| D23 | timeOfDayAdjustments | `1.249531*rowMeans(|r_i*r_{i+1}*r_{i+2}|^(2/3))` | 论文无此常数 | 1.249531 来源待验证 |

### 6.4 测试矩阵

| 文件 | 测试数 | 容差 | 关键场景 |
|---|---|---|---|
| test_spread_cleaner.cpp | 4 | 1e-12 | rmLargeSpread 每日中位数; rmNegativeSpread; spreadPrices 长宽转换; 空数据 |
| test_liquidity_measures.cpp | 6 | 1e-10 | 23 种度量数值; getTradeDirection tick rule; midpoint 覆盖; 用户 DIRECTION; realizedSpread 越界; depthImbalanceRatio |
| test_amihud.cpp | 2 | 1e-12 | 基本计算; 零成交额异常 |
| test_aj_jump_test.cpp | 4 | 1e-10 | p=4,k=2 默认; p=2,k=3; fmupk 查表; calculateV 数值 |
| test_jo_jump_test.cpp | 3 | 1e-10 | power=4 默认; power=6; rollApplyProdWrapper 窗口 |
| test_intraday_jump_test.cpp | 3 | 1e-8 | RM 模式 rBPCov; Lee-Mykland 临界值; lookBackPeriod 敏感性 |
| test_rank_jump_test.cpp | 3 | 1e-8 | jumpDetection TOD; SVD 分解; bootstrap 临界值 (固定种子) |
| **合计** | **25** | | |

### 6.5 v1.4.3 任务清单

1. [ ] 实现 `liquidity/spread_cleaner.hpp` (rm_large_spread + rm_negative_spread + spread_prices)
2. [ ] 实现 `liquidity/liquidity_measures.hpp` (get_trade_direction + get_liquidity_measures 23 种度量)
3. [ ] 实现 `liquidity/amihud.hpp`
4. [ ] 实现 `tests/aj_jump_test.hpp` (含 calculateV/calculateNpk/fmupk 查表)
5. [ ] 实现 `tests/jo_jump_test.hpp` (含 roll_apply_prod_wrapper + simre + RBPVar)
6. [ ] 实现 `tests/intraday_jump_test.hpp` (RM 模式 + Lee-Mykland 临界值)
7. [ ] 实现 `tests/rank_jump_test.hpp` (含 jump_detection + SVD + bootstrap)
8. [ ] 编写 7 个测试文件 (25 测试)
9. [ ] tests/CMakeLists.txt 注册新测试目标
10. [ ] 全量回归 (目标 1387/1387)
11. [ ] A/B 站 GCC 跨平台验证
12. [ ] 更新 DEVELOPMENT_LOG + README
13. [ ] git commit + push

---

## 7. CI Gate 设计

### 7.1 R 基准生成流程

1. **环境要求**: 主控站安装 R + `highfrequency` 包 (A/B 站非必需)
2. **脚本**: `tests/fixtures/hfe/generate_r_baselines.R`
3. **输出**: `tests/fixtures/hfe/baselines.json` (JSON 格式, 供 C++ 测试读取)
4. **触发**: R 基准仅在主控站生成, A/B 站使用版本控制中的同一份 JSON

### 7.2 容差层级

| 层级 | 容差 | 适用场景 |
|---|---|---|
| 严格 | 1e-12 | 无噪声合成数据, 纯 GBM 模拟 |
| 标准 | 1e-10 | 有微结构噪声的实际数据 (默认) |
| 宽松 | 1e-8 | 迭代算法 (HAR/HEAVY 估计) |

### 7.3 回归基线

- **基线**: 1268 测试 (Phase 1-4 + v1.4 PDE 三波)
- **v1.4.0 完成后**: ~1283 测试
- **全量回归**: `ctest -C Release --parallel 8` 三平台 (主控 MSVC + A 站 GCC + B 站 GCC)

---

## 8. 风险评估

| 风险 | 概率 | 影响 | 缓解 |
|---|---|---|---|
| R `highfrequency` 安装失败 | 中 | 高 | 主控站单独验证, A/B 站使用 JSON 基准 |
| ITCH 5.0 解析错误 | 中 | 中 | v1.4.0 先实现 CSV, ITCH 推迟到 v1.4.1 |
| R 与 C++ 数值精度差异 | 低 | 高 | 容差从 1e-8 开始, 稳定后收紧到 1e-10 |
| HFE 与定价栈意外耦合 | 低 | 中 | 严格保持 `hfecon/` 独立, 只依赖 `core/` |
| 性能未达 50 Mtick/s | 中 | 中 | **推迟至 Phase 6 独立性能优化波次** (v1.4.2-v1.4.3 功能正确性已达标, SIMD/OpenMP 优化不阻塞 HFE 模块发布; 见 AUDIT_CHECKLIST F1/F2/F4/F5) |

---

## 9. 验收标准

### 9.1 v1.4.0 验收

1. `include/cpphub/hfecon/data/taq_reader.hpp` 实现完整 (CSV + ITCH 5.0)
2. `include/cpphub/hfecon/measures/realized_measures.hpp` 实现 RV/RVol/RQ/BPV/RSV
3. `include/cpphub/hfecon/tests/bns_jump_test.hpp` 实现 BNS 检验
4. `tests/unit/hfecon/test_realized_measures.cpp` 15 个测试全部通过
5. R 对照基准 JSON 生成, 所有 `EXPECT_NEAR` 在 1e-10 内通过
6. 全量回归 1283+ 测试通过, 无回归
7. 三平台 (主控 MSVC + A/B 站 GCC) 一致

### 9.2 文档要求

- 更新 `README.md` 添加 HFE 章节
- 更新 `DEVELOPMENT_LOG.md` 记录 v1.4.0 实施过程
- 更新 `project_memory.md` 记录关键决策

---

## 10. 待执行任务清单

### v1.4.0 立即可执行

1. [ ] 在主控站安装 R + `highfrequency` 包, 生成 `baselines.json`
2. [ ] 创建 `include/cpphub/hfecon/` 目录结构
3. [ ] 实现 `taq_reader.hpp` (CSV 优先, ITCH 可推迟)
4. [ ] 实现 `realized_measures.hpp` (RV/RVol/RQ/BPV/RSV)
5. [ ] 实现 `bns_jump_test.hpp`
6. [ ] 编写 `test_realized_measures.cpp` (15 测试)
7. [ ] 全量回归测试
8. [ ] 更新文档, 提交 git

---

**Spec 编写完成时间**: 2026-08-02
**学术依据**: 7 篇核心文献 (BN-S 2002/2004/2006, ABD 2003, H-L 2006, BKS 2022, A-J 2009)
**幻觉排除**: 所有数学公式均来自已验证 DOI 文献, R 函数签名来自 CRAN 官方文档与 JSS vignette
