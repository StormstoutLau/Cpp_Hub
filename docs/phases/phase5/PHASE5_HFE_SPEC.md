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

## 4. 第二波 (v1.4.1) 交付项 (待 v1.4.0 完成后细化)

**目标**: 微结构噪声分析 + realized kernel + 稀疏抽样

| 项 | 文件 | R 对照 | 文献 |
|---|---|---|---|
| 稀疏抽样 | `measures/sparse_sampling.hpp` | `sparseSampling` | [H-L 2006] |
| Realized Kernel | `measures/realized_kernel.hpp` | `rKernelCov` | Barndorff-Nielsen et al. (2008) |
| 噪声方差估计 | `noise/noise_variance.hpp` | `noiseBPM` | [H-L 2006] |
| 最优抽样频率 | `noise/optimal_frequency.hpp` | `optimFrequ` | [H-L 2006] |
| 噪声自相关 | `noise/noise_autocorr.hpp` | `noiseAC` | [H-L 2006] |

---

## 5. 第三波 (v1.4.2) 交付项 (待 v1.4.1 完成后细化)

**目标**: HAR 模型 + HEAVY 模型 + RV 预测

| 项 | 文件 | R 对照 | 文献 |
|---|---|---|---|
| HAR 模型 | `models/har_model.hpp` | `HARmodel` | Corsi (2009) |
| HEAVY 模型 | `models/heavy_model.hpp` | `HEAVYmodel` | Shephard & Sheppard (2010) |
| ARFIMA 估计 | `models/arfima.hpp` | (R `forecast`) | [ABD 2003] |
| 预测精度评估 | `models/forecast_eval.hpp` | (无直接对照) | Patton (2011) |

---

## 6. 第四波 (v1.4.3) 交付项 (待 v1.4.2 完成后细化)

**目标**: 流动性度量 + 多资产 + 高级跳跃检验

| 项 | 文件 | R 对照 | 文献 |
|---|---|---|---|
| 有效价差 | `liquidity/effective_spread.hpp` | `getLiquidityMeasures` | Hasbrouck (2009) |
| 实现价差 | `liquidity/realized_spread.hpp` | `getLiquidityMeasures` | Hasbrouck (2009) |
| Amihud 流动性 | `liquidity/amihud.hpp` | (无直接对照) | Amihud (2002) |
| AJ 跳跃检验 | `tests/aj_jump_test.hpp` | `AJjumpTest` | [A-J 2009] |
| 门限跳跃检验 | `tests/threshold_jump_test.hpp` | (无直接对照) | Corsi, Pirino, Renò (2010) |

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
| 性能未达 50 Mtick/s | 中 | 中 | v1.4.2 引入 SIMD/OpenMP 优化 |

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
