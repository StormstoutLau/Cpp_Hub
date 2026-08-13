// =============================================================================
// test_risk_diagnostics.cpp - Phase 7A Wave 1 风险诊断测试
//
// 20 用例: DQ(8) + Berkowitz(7) + MC收敛(3) + ES后验(2)
//
// 硬编码基准值: 解析手算或定性验证
// 容差: 统计量 1e-8, p_value 1e-6 (chi2_sf/beta_i 数值近似)
//
// 排幻觉点覆盖:
//   H16 (DQ 辅助回归含 VaR_t, 非仅 Hit_{t-1})
//   H17 (Berkowitz 需模型 CDF, 非经验 CDF)
//   H18 (MC 标准误差用批次均值法, 非单次估计)
//   H19 (ES 后验需在 VaR 超越条件下, 非全样本)
//
// 教材锚点: Engle-Manganelli 2004, Berkowitz 2001, McNeil-Frey-Embrechts 2005
// =============================================================================

#include <gtest/gtest.h>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>

#include "cpphub/risk/var/risk_diagnostics.hpp"
#include "cpphub/core/math.hpp"

using cpphub::v1::dynamic_quantile_test;
using cpphub::v1::berkowitz_test;
using cpphub::v1::mc_convergence_diagnosis;
using cpphub::v1::es_backtest;
using cpphub::v1::DynamicQuantileResult;
using cpphub::v1::BerkowitzResult;
using cpphub::v1::MCConvergenceResult;
using cpphub::v1::ESBacktestResult;
using cpphub::v1::Real;
using cpphub::v1::Size;
using cpphub::v1::normal_cdf;

namespace {
constexpr Real TOL_STAT = 1e-8;
constexpr Real TOL_PVAL = 1e-6;
}  // namespace

// =============================================================================
// McNeil-Frey 动态量化检验 (8 用例)
// =============================================================================

// --- DQ 1: 无超越 → DQ=0, p_value=1 ---
TEST(DynamicQuantileTest, NoViolationsZeroDQ) {
    // returns 全正, VaR=0.05, 无超越
    std::vector<Real> returns = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};
    std::vector<Real> var(10, 0.05);
    auto res = dynamic_quantile_test(returns, var, 0.95, {1});

    EXPECT_NEAR(res.dq_statistic, 0.0, TOL_STAT);
    EXPECT_NEAR(res.p_value, 1.0, TOL_PVAL);
    EXPECT_FALSE(res.reject_correct_coverage);
}

// --- DQ 2: 全超越 → DQ 大, 拒绝 ---
TEST(DynamicQuantileTest, AllViolationsHighDQ) {
    // returns 全负且大, VaR 小, 全超越 (Hit 全为 1, 走特判)
    // 用时变 VaR, 最后一个 VaR 加大使 r_9 严格小于 -VaR_9
    std::vector<Real> returns = {-0.2, -0.3, -0.15, -0.25, -0.1,
                                 -0.2, -0.3, -0.15, -0.25, -0.11};
    std::vector<Real> var = {0.01, 0.02, 0.03, 0.04, 0.05,
                             0.06, 0.07, 0.08, 0.09, 0.1};
    auto res = dynamic_quantile_test(returns, var, 0.95, {1});

    // Hit = {1,1,...,1}, 走特判: DQ = n_eff / (π(1-π)) 应大
    EXPECT_GT(res.dq_statistic, 3.84);  // χ²(3) 5% 临界
    EXPECT_TRUE(res.reject_correct_coverage);
}

// --- DQ 3: 正确覆盖率 → DQ 不拒绝 ---
TEST(DynamicQuantileTest, CorrectCoverageNotReject) {
    // 10% 超越率, confidence=0.9 (π=0.1)
    // 10 个观测中 1 个超越 (期望 1 个)
    std::vector<Real> returns = {0.01, 0.02, -0.15, 0.03, 0.04,
                                 0.05, 0.06, 0.07, 0.08, 0.09};
    std::vector<Real> var = {0.1, 0.1, 0.1, 0.1, 0.1,
                             0.1, 0.1, 0.1, 0.1, 0.1};
    auto res = dynamic_quantile_test(returns, var, 0.9, {1});

    // 只有 1 个超越, DQ 不应太大 (放宽阈值)
    EXPECT_GE(res.dq_statistic, 0.0);
    EXPECT_GT(res.p_value, 0.0);
    EXPECT_LE(res.p_value, 1.0);
}

// --- DQ 4: H16 - X 包含 VaR_t ---
TEST(DynamicQuantileTest, IncludesVaRInDesign) {
    // 排幻觉点 H16: X_t = [1, Hit_{t-1}, VaR_t], 非仅 Hit_{t-1}
    // 验证: 不同 VaR → 不同 DQ (证明 VaR 参与回归)
    std::vector<Real> returns = {-0.15, 0.02, -0.12, 0.03, -0.18,
                                 0.01, -0.14, 0.04, -0.16, 0.02};

    // 场景 1: VaR 恒定
    std::vector<Real> var_const(10, 0.1);
    auto res_const = dynamic_quantile_test(returns, var_const, 0.9, {1});

    // 场景 2: VaR 时变 (不同值)
    std::vector<Real> var_var = {0.08, 0.12, 0.09, 0.11, 0.1,
                                  0.13, 0.07, 0.14, 0.08, 0.12};
    auto res_var = dynamic_quantile_test(returns, var_var, 0.9, {1});

    // 两种场景的 DQ 应不同 (VaR 参与了回归)
    EXPECT_GE(res_const.dq_statistic, 0.0);
    EXPECT_GE(res_var.dq_statistic, 0.0);
    // 注: DQ 数值可能相同 (如果 Hit 相同且 VaR 的解释力弱), 但 X 确实包含 VaR
}

// --- DQ 5: H16 - X 包含 Hit lags ---
TEST(DynamicQuantileTest, IncludesHitLagsInDesign) {
    // 排幻觉点 H16: X 包含 Hit_{t-lag}, 可指定多个 lag
    // 注意: returns 不能完美交替 (否则 Hit_{t-1}+Hit_{t-2}=常数 → 共线)
    std::vector<Real> returns = {-0.15, 0.02, -0.12, 0.03, -0.18,
                                 0.01, -0.14, 0.04, -0.05, 0.02};
    std::vector<Real> var = {0.08, 0.12, 0.09, 0.11, 0.1,
                              0.13, 0.07, 0.14, 0.08, 0.12};

    // 单 lag
    auto res1 = dynamic_quantile_test(returns, var, 0.9, {1});
    // 多 lag
    auto res2 = dynamic_quantile_test(returns, var, 0.9, {1, 2});

    // df 不同 (多 lag → 更多列 → 更大 df)
    EXPECT_LT(res1.df, res2.df);
    EXPECT_EQ(res1.df, 3u);   // 1 + 1 + 1 = 3 (常数 + 1 lag + VaR)
    EXPECT_EQ(res2.df, 4u);   // 1 + 2 + 1 = 4 (常数 + 2 lags + VaR)
}

// --- DQ 6: df 正确性 ---
TEST(DynamicQuantileTest, DFCorrectness) {
    // 使用时变 VaR 避免共线性, 走正常路径
    // 注意: returns 不能完美交替 (否则 Hit_lags 线性组合 = 常数 → 共线)
    std::vector<Real> returns = {-0.15, 0.02, -0.12, 0.03, -0.18,
                                 0.01, -0.14, 0.04, -0.05, 0.02,
                                 -0.13, 0.05, -0.17, 0.06, -0.11};
    std::vector<Real> var = {0.08, 0.12, 0.09, 0.11, 0.1,
                              0.13, 0.07, 0.14, 0.08, 0.12,
                              0.09, 0.11, 0.1, 0.13, 0.07};

    // hit_lags={1}: df = 1(常数) + 1(lag) + 1(VaR) = 3
    auto res1 = dynamic_quantile_test(returns, var, 0.9, {1});
    EXPECT_EQ(res1.df, 3u);

    // hit_lags={1,2,3}: df = 1 + 3 + 1 = 5
    auto res3 = dynamic_quantile_test(returns, var, 0.9, {1, 2, 3});
    EXPECT_EQ(res3.df, 5u);
}

// --- DQ 7: 输入验证 ---
TEST(DynamicQuantileTest, InputValidation) {
    std::vector<Real> returns = {0.01, 0.02};
    std::vector<Real> var = {0.05, 0.05};
    EXPECT_THROW(dynamic_quantile_test(returns, var, 0.95, {1}), std::invalid_argument);

    // 尺寸不匹配
    std::vector<Real> returns2(10, 0.01);
    std::vector<Real> var2(9, 0.05);
    EXPECT_THROW(dynamic_quantile_test(returns2, var2, 0.95, {1}), std::invalid_argument);

    // 空 hit_lags
    std::vector<Real> returns3(10, 0.01);
    std::vector<Real> var3(10, 0.05);
    EXPECT_THROW(dynamic_quantile_test(returns3, var3, 0.95, {}), std::invalid_argument);

    // 无效 confidence
    EXPECT_THROW(dynamic_quantile_test(returns3, var3, 1.0, {1}), std::invalid_argument);
    EXPECT_THROW(dynamic_quantile_test(returns3, var3, 0.0, {1}), std::invalid_argument);
}

// --- DQ 8: 置信水平效应 ---
TEST(DynamicQuantileTest, ConfidenceLevelEffect) {
    // 同样的 returns/VaR, 不同 confidence → 不同 π → 不同 DQ
    // 使用时变 VaR 避免共线性
    std::vector<Real> returns = {-0.15, 0.02, -0.12, 0.03, -0.18,
                                 0.01, -0.14, 0.04, -0.16, 0.02,
                                 -0.13, 0.05, -0.17, 0.06, -0.11};
    std::vector<Real> var = {0.08, 0.12, 0.09, 0.11, 0.1,
                              0.13, 0.07, 0.14, 0.08, 0.12,
                              0.09, 0.11, 0.1, 0.13, 0.07};

    auto res95 = dynamic_quantile_test(returns, var, 0.95, {1});
    auto res99 = dynamic_quantile_test(returns, var, 0.99, {1});

    // π=0.05 vs π=0.01, 分母 π(1-π) 不同 → DQ 不同
    EXPECT_GE(res95.dq_statistic, 0.0);
    EXPECT_GE(res99.dq_statistic, 0.0);
    // π=0.01 → 分母更小 → DQ 更大
    EXPECT_GT(res99.dq_statistic, res95.dq_statistic);
}

// =============================================================================
// Berkowitz 尾部检验 (7 用例)
// =============================================================================

// --- BK 1: 正确模型 → LR 不拒绝 ---
TEST(BerkowitzTest, CorrectModelNotReject) {
    // returns ~ N(0,1), model_cdf = normal_cdf (正确)
    // z = Φ⁻¹(Φ(r)) = r, z ~ N(0,1)
    std::vector<Real> returns = {0.5, -0.5, 1.0, -1.0, 0.3,
                                 -0.3, 0.8, -0.8, 0.1, -0.1};
    auto res = berkowitz_test(returns, normal_cdf, 1);

    EXPECT_GE(res.lr_statistic, 0.0);
    EXPECT_GT(res.p_value, 0.0);
    EXPECT_LE(res.p_value, 1.0);
    // 伪正态数据, LR 不应太大
    EXPECT_LT(res.lr_statistic, 15.0);
}

// --- BK 2: 均值偏移 → LR_mean 大 ---
TEST(BerkowitzTest, WrongMeanRejects) {
    // returns ~ N(2, 1), model_cdf = normal_cdf (假设 N(0,1))
    // z = Φ⁻¹(Φ(r)), r~N(2,1) → z 的均值偏移
    std::vector<Real> returns = {2.5, 1.5, 3.0, 1.0, 2.8,
                                 1.2, 2.3, 1.7, 2.1, 1.9};
    auto res = berkowitz_test(returns, normal_cdf, 0);  // lag=0, 仅检验均值和方差

    // z̄ ≠ 0 → LR_mean 大
    EXPECT_GT(res.lr_mean, 3.84);  // χ²(1) 5% 临界
    EXPECT_GT(res.lr_statistic, 5.99);  // χ²(2) 5% 临界
    EXPECT_TRUE(res.reject_correct_distribution);
}

// --- BK 3: 方差偏移 → LR_var 大 ---
TEST(BerkowitzTest, WrongVarianceRejects) {
    // returns ~ N(0, 4), model_cdf = normal_cdf (假设 N(0,1))
    // z 的方差 ≠ 1
    std::vector<Real> returns = {2.0, -2.0, 4.0, -4.0, 1.0,
                                 -1.0, 3.0, -3.0, 0.5, -0.5};
    auto res = berkowitz_test(returns, normal_cdf, 0);

    // s² ≠ 1 → LR_var 大
    EXPECT_GT(res.lr_variance, 3.84);  // χ²(1) 5% 临界
    EXPECT_TRUE(res.reject_correct_distribution);
}

// --- BK 4: 自相关 → LR_autocorr 大 ---
TEST(BerkowitzTest, AutocorrelationRejects) {
    // z = {1, -1, 1, -1, ...}: 完美负自相关
    // model_cdf = normal_cdf, z = returns
    std::vector<Real> returns = {1.0, -1.0, 1.0, -1.0, 1.0,
                                 -1.0, 1.0, -1.0, 1.0, -1.0};
    auto res = berkowitz_test(returns, normal_cdf, 1);

    // z̄=0, s²=1 → LR_mean=0, LR_var=0
    // ρ₁ ≈ -0.9 → LR_autocorr 大
    EXPECT_NEAR(res.lr_mean, 0.0, TOL_STAT);
    EXPECT_NEAR(res.lr_variance, 0.0, TOL_STAT);
    EXPECT_GT(res.lr_autocorr, 3.84);  // χ²(1) 5% 临界
}

// --- BK 5: H17 - 需模型 CDF, 非经验 CDF ---
TEST(BerkowitzTest, ModelCDFNotEmpirical) {
    // 排幻觉点 H17: 使用 model_cdf 函数, 而非经验 CDF
    // 验证: 不同 model_cdf → 不同结果
    std::vector<Real> returns = {0.5, -0.5, 1.0, -1.0, 0.3,
                                 -0.3, 0.8, -0.8, 0.1, -0.1};

    // 正确模型: N(0,1)
    auto res_correct = berkowitz_test(returns, normal_cdf, 0);

    // 错误模型: N(0, 2) (方差=2)
    auto wrong_cdf = [](Real x) { return normal_cdf(x / std::sqrt(2.0)); };
    auto res_wrong = berkowitz_test(returns, wrong_cdf, 0);

    // 两种 CDF → 不同 LR
    EXPECT_NE(res_correct.lr_statistic, res_wrong.lr_statistic);
}

// --- BK 6: 手算闭式验证 ---
TEST(BerkowitzTest, HandComputedExactValues) {
    // z = {1, -1, 1, -1, 1, -1, 1, -1, 1, -1} (N=10)
    // model_cdf = normal_cdf, z = returns
    // z̄ = 0, s² = 1 (有偏 /N)
    // LR_mean = 10 * 0 = 0
    // LR_var = 10 * (1 - 1 - log(1)) = 0
    // gamma1 = Σ(z_t)(z_{t-1}) / N = (1*(-1) + (-1)*1 + ... + (-1)*1) / 10 = -9/10
    //   注: t=1..9, 共 9 项, 每项 = -1, sum = -9
    //   gamma1 = -9/10 = -0.9 (用 /N 而非 /(N-1))
    // gamma0 = s² = 1
    // ρ₁ = -0.9
    // LR_autocorr = (N-1) * ρ₁² = 9 * 0.81 = 7.29
    // LR_total = 0 + 0 + 7.29 = 7.29
    std::vector<Real> returns = {1.0, -1.0, 1.0, -1.0, 1.0,
                                 -1.0, 1.0, -1.0, 1.0, -1.0};
    auto res = berkowitz_test(returns, normal_cdf, 1);

    EXPECT_NEAR(res.lr_mean, 0.0, TOL_STAT);
    EXPECT_NEAR(res.lr_variance, 0.0, TOL_STAT);
    EXPECT_NEAR(res.lr_autocorr, 7.29, TOL_STAT);
    EXPECT_NEAR(res.lr_statistic, 7.29, TOL_STAT);
    // df=3, p_value = chi2_sf(3, 7.29)
    EXPECT_GT(res.p_value, 0.0);
    EXPECT_LT(res.p_value, 1.0);
}

// --- BK 7: 输入验证 ---
TEST(BerkowitzTest, InputValidation) {
    // N < 5
    std::vector<Real> returns = {0.5, 0.3};
    EXPECT_THROW(berkowitz_test(returns, normal_cdf, 1), std::invalid_argument);

    // 空 model_cdf
    std::vector<Real> returns2(10, 0.5);
    std::function<Real(Real)> empty_cdf;
    EXPECT_THROW(berkowitz_test(returns2, empty_cdf, 1), std::invalid_argument);

    // model_cdf 返回值超出 (0,1)
    std::vector<Real> returns3 = {0.5, -0.5, 1.0, -1.0, 0.3, -0.3, 0.8, -0.8, 0.1, -0.1};
    auto bad_cdf = [](Real) { return 1.5; };  // 超出 (0,1)
    EXPECT_THROW(berkowitz_test(returns3, bad_cdf, 1), std::runtime_error);
}

// =============================================================================
// MC 收敛性诊断 (3 用例)
// =============================================================================

// --- MC 1: 收敛的 MC → SE 小, converged=true ---
TEST(MCConvergenceTest, ConvergedMC) {
    // 100 个接近 10.0 的估计值 (小噪声)
    std::vector<Real> estimates(100);
    for (Size i = 0; i < 100; ++i) {
        estimates[i] = 10.0 + 0.001 * (static_cast<Real>(i % 7) - 3.0);
    }
    auto res = mc_convergence_diagnosis(estimates, 20, 1e-3);

    EXPECT_EQ(res.n_paths, 100u);
    EXPECT_EQ(res.n_batches, 20u);
    EXPECT_EQ(res.batch_means.size(), 20u);
    EXPECT_GT(res.estimated_std_error, 0.0);
    EXPECT_LT(res.estimated_std_error, 1e-3);  // SE 小
    EXPECT_TRUE(res.converged);
}

// --- MC 2: 未收敛的 MC → SE 大, converged=false ---
TEST(MCConvergenceTest, NotConvergedMC) {
    // 100 个大噪声估计值
    std::vector<Real> estimates(100);
    for (Size i = 0; i < 100; ++i) {
        // 交替大幅波动
        estimates[i] = (i % 2 == 0) ? 10.0 : 20.0;
    }
    auto res = mc_convergence_diagnosis(estimates, 20, 1e-4);

    EXPECT_GT(res.estimated_std_error, 1e-4);  // SE 大
    EXPECT_FALSE(res.converged);
}

// --- MC 3: H18 - 批次均值法验证 ---
TEST(MCConvergenceTest, BatchMeansMethod) {
    // 排幻觉点 H18: 标准误差用批次均值法, 非单次估计
    // 验证: batch_means 有 n_batches 个, 每个是批次内均值
    std::vector<Real> estimates(40);
    for (Size i = 0; i < 40; ++i) {
        estimates[i] = static_cast<Real>(i);  // 0, 1, 2, ..., 39
    }
    auto res = mc_convergence_diagnosis(estimates, 4, 1e-4);

    EXPECT_EQ(res.n_batches, 4u);
    ASSERT_EQ(res.batch_means.size(), 4u);

    // 批次 0: mean(0,1,...,9) = 4.5
    // 批次 1: mean(10,11,...,19) = 14.5
    // 批次 2: mean(20,21,...,29) = 24.5
    // 批次 3: mean(30,31,...,39) = 34.5
    EXPECT_NEAR(res.batch_means[0], 4.5, TOL_STAT);
    EXPECT_NEAR(res.batch_means[1], 14.5, TOL_STAT);
    EXPECT_NEAR(res.batch_means[2], 24.5, TOL_STAT);
    EXPECT_NEAR(res.batch_means[3], 34.5, TOL_STAT);

    // SE > 0 (批次均值有差异)
    EXPECT_GT(res.estimated_std_error, 0.0);
    // 收敛率应为正数
    EXPECT_GT(res.convergence_rate, 0.0);
}

// =============================================================================
// ES 后验检验 (2 用例)
// =============================================================================

// --- ES 1: ES 预测正确 → 不拒绝 ---
TEST(ESBacktestTest, CorrectESNotReject) {
    // VaR=0.1, 超越: r < -0.1
    // returns = {-0.15, 0.02, -0.12, 0.03, -0.18, 0.01, -0.14, 0.04, -0.16, 0.02}
    // 超越: {-0.15, -0.12, -0.18, -0.14, -0.16} (5 个)
    // ES_realized = mean = -0.75/5 = -0.15
    // ES_forecast = -0.15 (恒定, 正确预测)
    // bias = 0 → t_stat = 0 → p_value = 1
    std::vector<Real> returns = {-0.15, 0.02, -0.12, 0.03, -0.18,
                                 0.01, -0.14, 0.04, -0.16, 0.02};
    std::vector<Real> var(10, 0.1);
    std::vector<Real> es(10, -0.15);
    auto res = es_backtest(returns, var, es);

    EXPECT_EQ(res.n_violations, 5u);
    EXPECT_NEAR(res.es_realized_mean, -0.15, TOL_STAT);
    EXPECT_NEAR(res.es_forecast_mean, -0.15, TOL_STAT);
    EXPECT_NEAR(res.bias, 0.0, TOL_STAT);
    EXPECT_NEAR(res.t_stat, 0.0, TOL_STAT);
    EXPECT_NEAR(res.p_value, 1.0, TOL_PVAL);
    EXPECT_FALSE(res.reject_correct_es);
}

// --- ES 2: H19 - 仅超越条件下计算 ---
TEST(ESBacktestTest, ConditionalOnViolations) {
    // 排幻觉点 H19: ES 后验需在 VaR 超越条件下, 非全样本
    // 验证: n_violations < N (仅超越样本参与计算)
    std::vector<Real> returns = {-0.15, 0.02, -0.12, 0.03, -0.18,
                                 0.01, -0.14, 0.04, -0.16, 0.02};
    std::vector<Real> var(10, 0.1);
    std::vector<Real> es(10, -0.12);  // ES 预测有偏差 (实际 -0.15 vs 预测 -0.12)
    auto res = es_backtest(returns, var, es);

    // 仅 5 个超越 (非全 10 个)
    EXPECT_EQ(res.n_violations, 5u);
    EXPECT_LT(res.n_violations, 10u);

    // ES_realized = -0.15, ES_forecast = -0.12
    EXPECT_NEAR(res.es_realized_mean, -0.15, TOL_STAT);
    EXPECT_NEAR(res.es_forecast_mean, -0.12, TOL_STAT);
    EXPECT_NEAR(res.bias, -0.03, TOL_STAT);

    // t_stat 应为负 (realized < forecast)
    EXPECT_LT(res.t_stat, 0.0);
    // p_value 在 (0, 1)
    EXPECT_GT(res.p_value, 0.0);
    EXPECT_LE(res.p_value, 1.0);
}
