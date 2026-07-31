# P4: CDO Base Correlation 标定任务

## 任务定位
**目标**: 在现有 `include/cpphub/instruments/credit/cdo.hpp` 基础上, 实现 Base Correlation 标定与异质分券定价功能, 完成 CDO 框架收尾。
**工作站**: A 站 (`scott-lau-NEX.local`)
**模型**: `opencode/deepseek-v4-flash-free` (opencode Zen 免费模型)
**工作目录**: `/tmp/oc_codetest`

## 关键约束 (来自 project_memory)
1. A 站 opencode 路径为 `/snap/bin/opencode` (1.18.8), 非交互 SSH 必须用全路径
2. **A 站必须设置代理** `export http_proxy=http://127.0.0.1:7890` (mihomo 在 /root/clashctl/, systemd 管理)
3. 使用 `ssh -t` 模式以继承 keyring/auth 实现 PTY 交互
4. 头文件 `#include` 必须放在 namespace 外, 禁止嵌套命名空间
5. 模型选择: `opencode/deepseek-v4-flash-free` (A 站实测 100% 通过率, 16/16)
6. 使用 `--auto` 参数允许 agent 自主读写

## 现有代码状态

### `include/cpphub/instruments/credit/cdo.hpp` (已有功能)
- `CDOTrancheConfig{attachment, detachment, spread, maturity, n_premiums, ...}` — 分券配置
- `CDOTrancheResult{pv_premium_leg, pv_protection_leg, pv, par_spread, risky_pv01, expected_loss, ...}` — 定价结果
- `CDOLHPPricer(rho, pd_maturity, lgd, discount, maturity)` — LHP 大组合定价器
  - `expected_tranche_loss(t, A, D)` — 期望损失 (Gauss-Hermite 积分)
  - `expected_remaining_notional(t, A, D)`
  - `price(cfg)` — 完整分券定价
- `CDOMCPricer` — MC 定价器

### 测试文件 `tests/unit/instruments/test_cdo.cpp` (已有完整测试)

## 需要扩展的功能

### 1. Compound Correlation 反推 (单分券隐含相关系数)
**目的**: 给定一个分券的市场 par spread, 反推使模型 par spread 等于市场价的 ρ。
**算法**: Brent 或 bisection 一维求根
- 已知: market_par_spread, A, D, pd, lgd, discount, maturity
- 求: ρ ∈ [0, 0.95] 使 CDOLHPPricer(ρ,...).price(cfg).par_spread = market_par_spread
- 注意: compound correlation 可能不唯一 (smile 现象), 在 detachment > 0.10 时可能无解

**接口**:
```cpp
class CDOBaseCorrelationCalibrator {
public:
    CDOBaseCorrelationCalibrator(Real pd_maturity, Real lgd,
                                  const ZeroCurve& discount, Real maturity);

    // Compound correlation (单分券隐含 rho)
    // 返回 rho 使 CDOLHPPricer.price(cfg).par_spread = market_spread
    // 返回 nan 表示无解 (rho > 0.95 或 < 0)
    Real compound_correlation(Real A, Real D, Real market_spread,
                               Real tol = 1e-6, Size max_iter = 100) const;

    // 标定 base correlation curve (CDX/iTraxx 标准分券)
    // 输入: 标准分券 detachment points [0.03, 0.07, 0.10, 0.15, 0.30]
    //      对应 market par spreads
    // 输出: 每个 detachment 的 base correlation
    // 算法: 从 0-3% 开始, 对 detachment D_k:
    //   1. 用已校准的 base corr (D_{k-1}, rho_{k-1}) 计算 0-D_{k-1} 损失
    //   2. 找 rho_k 使 0-D_k 分券 par spread = market_spread_k
    //   3. 这里 0-D_k 分券的损失 = D_{k-1}-D_k 分券损失 + 0-D_{k-1} 损失
    std::vector<Real> calibrate_base_correlation(
        const std::vector<Real>& detachments,  // [0.03, 0.07, 0.10, 0.15, 0.30]
        const std::vector<Real>& market_spreads,
        Real attachment_first = 0.0) const;

    // 用 base correlation curve 给非标准分券定价
    // 算法: 对 detachment D 和 attachment A:
    //   1. 在 base corr curve 上插值得 rho_D 和 rho_A
    //   2. E[L_{0-D}(t)] 用 rho_D 计算
    //   3. E[L_{0-A}(t)] 用 rho_A 计算
    //   4. E[L_{A-D}(t)] = E[L_{0-D}(t)] - E[L_{0-A}(t)]
    //   5. 用此损失曲线计算 par spread
    Real price_off_market_tranche(
        Real A, Real D,
        const std::vector<Real>& base_detachments,
        const std::vector<Real>& base_correlations,
        Real spread = 0.0) const;  // spread=0 时返回 par spread

private:
    Real pd_maturity_;
    Real lgd_;
    const ZeroCurve& discount_;
    Real maturity_;
    // 线性插值 base corr curve (extrapolation flat)
    Real interpolate_base_corr(Real D,
                                const std::vector<Real>& detachments,
                                const std::vector<Real>& corrs) const;
};
```

### 2. Base Correlation 单调性校验
- base correlation 应**单调递增** ( detachment ↑ → rho ↑ )
- 违反单调性 → 标定不稳定 / 市场数据异常
- 输出: 每个分券的 base correlation + 一致性诊断

## 测试要求 (新增 ~12 测试, 加在 test_cdo.cpp 末尾)

### Compound Correlation (4 测试)
1. `CompoundCorrelationRecoversMarketSpread` — 用求得的 rho 重新 price, par_spread 应等于输入
2. `CompoundCorrelationMonotonicInSpread` — spread 越高, rho 越低 (高 spread = 高损失 = 低相关)
3. `CompoundCorrelationInRange` — ρ ∈ (0, 0.95)
4. `CompoundCorrelationEquityTrancheHighRho` — equity tranche (0-3%) 的 rho 应显著高于 senior

### Base Correlation Curve (4 测试)
5. `BaseCorrelationCurveMonotonicIncreasing` — base corr 随 detachment 单调递增
6. `BaseCorrelationRecoversMarketSpreads` — 用 base corr curve 重定价标准分券, par_spread 应等于市场价 (容差 1bp)
7. `BaseCorrelationEquityTrancheMatchesCompound` — 0-3% 分券的 base corr = compound corr
8. `BaseCorrelationInterpolationMidpoint` — 中间 detachment 的插值合理

### Off-Market Tranche Pricing (3 测试)
9. `OffMarketTrancheBetweenStandard` — 4%-6% 分券 (在 3%-7% 之间) 的 par spread 应在 3-7% 分券之间
10. `OffMarketTrancheSpreadMonotonic` — 固定 attachment, detachment ↑ → par_spread ↓
11. `OffMarketTrancheAtStandardMatches` — A=0.03, D=0.07 (标准分券) 应等于 base corr 标定输入

### 一致性诊断 (1 测试)
12. `BaseCorrelationNonMonotonicDetection` — 检测器能识别非单调的 base corr curve (模拟异常市场数据)

## 文献参考
- Li (2000) "On Default Correlation: A Copula Function Approach" J. Fixed Income
- Laurent & Gregory (2005) "Basket Default Swaps, CDOs and Factor Copulas"
- Hull & White (2004) "Valuation of a CDO and an n-th to Default CDS without Monte Carlo"
- O'Kane (2008) "Modelling Single-Name and Multi-Name Credit Derivatives" Ch.12-13
- **Base Correlation 关键文献**: 
  - McGinty & Ahluwalia (2004) "A model for base correlation calculation" JPMorgan
  - Willemann & Bicer (2005) "Base Correlation", Lehman Brothers
  - Reyfman (2004) "Valuing and Hedging Synthetic CDO Tranches Using Base Correlations", Bear Stearns

## 实施流程

1. **更新代码**: `cd /tmp/oc_codetest && git pull`
2. **opencode 启动命令** (A 站):
```bash
export http_proxy=http://127.0.0.1:7890
cd /tmp/oc_codetest
/snap/bin/opencode run -m opencode/deepseek-v4-flash-free --auto \
  "在 include/cpphub/instruments/credit/cdo.hpp 中新增 CDOBaseCorrelationCalibrator 类, 实现 P4 任务: (1) compound_correlation 单分券隐含 rho; (2) calibrate_base_correlation 标定 base correlation curve; (3) price_off_market_tranche 用 base corr 给非标准分券定价。在 tests/unit/instruments/test_cdo.cpp 末尾追加 12 个测试。所有现有测试必须仍然通过。编译: cmake --build build --config Release。运行测试: ctest --test-dir build --output-on-failure -R 'CDO|Tranche|Copula'。参考 docs/tasks/P4_CDO_BASE_CORR_TASK.md 中的接口规范和测试要求。"
```
3. **验收**: 测试通过, 跨平台编译成功 (GCC + MSVC), 提交 commit

## 验收标准
- [ ] 新增 CDOBaseCorrelationCalibrator 类 (在 cdo.hpp 末尾, namespace cpphub::v1 内)
- [ ] 新增 ≥12 个测试, 全部通过
- [ ] 全量回归测试 ≥ 986/986 通过
- [ ] 跨平台编译 (GCC + MSVC) 成功
- [ ] commit message: `feat(v1.3 P4): CDO base correlation calibration + off-market tranche pricing`
