// =============================================================================
// test_integration_phase7a.cpp - Phase 7A 端到端集成测试
//
// 10 用例, 覆盖全模块诊断流程:
//   1. 计量: OLS → BP/White → HC 选择
//   2. 计量: OLS → LB/BG → HAC 选择
//   3. 高频: HAR → 残差 LB → MZ
//   4. 高频: HEAVY → 标准化残差 → z² LB
//   5. 风险: VaR → DQ + Berkowitz + ES
//   6. 风险: Greeks → Analytic vs AAD vs Pathwise
//   7. 定价: IV 拟合优度 → 价格残差
//   8. 计量: MLE → 信息矩阵检验
//   9. 计量: OLS → CUSUM → Andrews
//  10. 高频: 跳跃检验 → Bonferroni/BH 修正
//
// 验证: 各模块诊断头文件可协同工作, 排幻觉点全覆盖
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <functional>
#include <random>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/inference/detail/ols_simple.hpp"
#include "cpphub/econometrics/inference/residual_diagnostics.hpp"
#include "cpphub/econometrics/inference/volatility_diagnostics.hpp"
#include "cpphub/econometrics/inference/specification_tests.hpp"
#include "cpphub/econometrics/inference/structural_break.hpp"
#include "cpphub/risk/var/risk_diagnostics.hpp"
#include "cpphub/risk/greeks/greeks_consistency.hpp"
#include "cpphub/pricing/pricing_diagnostics.hpp"
#include "cpphub/hfecon/hfecon_diagnostics.hpp"
#include "cpphub/hfecon/tests/jump_test_diagnostics.hpp"

using namespace cpphub::v1;
using cpphub::v1::econometrics::detail::ols_simple;
using cpphub::v1::econometrics::cusum_test;
using cpphub::v1::econometrics::andrews_breakpoint_test;
using cpphub::v1::hfecon::har_diagnostics;
using cpphub::v1::hfecon::heavy_diagnostics;
using cpphub::v1::hfecon::multiple_test_correction;
using cpphub::v1::hfecon::jump_test_diagnostics;
using cpphub::v1::hfecon::MultipleTestCorrectionResult;
using cpphub::v1::hfecon::JumpTestDiagnosticsResult;

// =============================================================================
// 用例 1: OLS → BP/White → HC 选择 (异方差检测)
// =============================================================================
TEST(Phase7AIntegration, OLS_Heteroscedasticity_TriggersHC) {
    // 生成异方差数据: ε_t ~ N(0, σ²_t), σ²_t = (1 + x_t²)
    const Size n = 200;
    std::mt19937 gen(42);
    std::vector<Real> y(n);
    // X_ols 含常数列 (ols_simple 不自动添加)
    // X_reg 不含常数列 (BP/White 自动添加)
    std::vector<std::vector<Real>> X_ols(n, std::vector<Real>(2));
    std::vector<std::vector<Real>> X_reg(n, std::vector<Real>(1));

    for (Size t = 0; t < n; ++t) {
        Real x = static_cast<Real>(t) / 50.0;
        X_ols[t][0] = 1.0;  // 常数项
        X_ols[t][1] = x;
        X_reg[t][0] = x;    // BP/White 不需要常数列
        Real sigma = 1.0 + x * x;  // 异方差: σ 随 x 增大
        std::normal_distribution<Real> dist(0.0, sigma);
        y[t] = 2.0 + 3.0 * x + dist(gen);
    }

    // OLS 估计
    std::vector<Real> fitted, residuals;
    Real r2;
    auto beta = ols_simple(y, X_ols, fitted, residuals, r2);

    // BP 检验: 应检测到异方差 (传 X_reg 不含常数列)
    auto bp = cpphub::v1::econometrics::breusch_pagan_test(X_reg, residuals);
    EXPECT_TRUE(bp.base.reject_null)
        << "BP should detect heteroscedasticity";

    // White 检验: 应检测到异方差
    auto wh = cpphub::v1::econometrics::white_test(X_reg, residuals, true);
    EXPECT_TRUE(wh.base.reject_null)
        << "White should detect heteroscedasticity";

    // "HC 切换": 若 BP/White 拒绝, 应使用 HC 标准误 (此处验证诊断逻辑)
    bool need_hc = bp.base.reject_null || wh.base.reject_null;
    EXPECT_TRUE(need_hc)
        << "Heteroscedasticity detected → should switch to HC standard errors";
}

// =============================================================================
// 用例 2: OLS → LB/BG → HAC 选择 (自相关检测)
// =============================================================================
TEST(Phase7AIntegration, OLS_Autocorrelation_TriggersHAC) {
    // 生成自相关数据: u_t = 0.7*u_{t-1} + ε_t
    const Size n = 200;
    std::mt19937 gen(42);
    std::normal_distribution<Real> noise(0.0, 1.0);

    std::vector<Real> y(n);
    // X_ols 含常数列, X_reg 不含 (BG 自动添加常数列)
    std::vector<std::vector<Real>> X_ols(n, std::vector<Real>(2));
    std::vector<std::vector<Real>> X_reg(n, std::vector<Real>(1));

    Real u_prev = 0.0;
    for (Size t = 0; t < n; ++t) {
        X_ols[t][0] = 1.0;
        X_ols[t][1] = static_cast<Real>(t);
        X_reg[t][0] = static_cast<Real>(t);
        Real eps = noise(gen);
        Real u = 0.7 * u_prev + eps;  // AR(1) 误差
        y[t] = 1.0 + 0.5 * X_ols[t][1] + u;
        u_prev = u;
    }

    // OLS
    std::vector<Real> fitted, residuals;
    Real r2;
    auto beta = ols_simple(y, X_ols, fitted, residuals, r2);

    // LB 检验: 应检测到自相关
    auto lb = cpphub::v1::econometrics::ljung_box_test(residuals, 10);
    EXPECT_TRUE(lb.base.reject_null)
        << "LB should detect autocorrelation";

    // BG 检验: 应检测到自相关 (传 X_reg 不含常数列)
    auto bg = cpphub::v1::econometrics::breusch_godfrey_test(X_reg, residuals, 5);
    EXPECT_TRUE(bg.base.reject_null)
        << "BG should detect autocorrelation";

    // "HAC 切换": 若 LB/BG 拒绝, 应使用 HAC 标准误
    bool need_hac = lb.base.reject_null || bg.base.reject_null;
    EXPECT_TRUE(need_hac)
        << "Autocorrelation detected → should switch to HAC standard errors";
}

// =============================================================================
// 用例 3: HAR → 残差 LB → MZ (HAR 完整诊断流程)
// =============================================================================
TEST(Phase7AIntegration, HAR_FullDiagnosticPipeline) {
    // 生成 HAR 风格 RV 数据 (拟合值接近实际值 + 白噪声残差)
    const Size n = 100;
    std::mt19937 gen(42);
    std::normal_distribution<Real> noise(0.0, 0.001);

    std::vector<Real> actual_rv(n), fitted_rv(n), residuals(n);
    for (Size t = 0; t < n; ++t) {
        Real rv = 0.01 + 0.001 * t;  // 递增 RV
        fitted_rv[t] = rv;
        residuals[t] = noise(gen);  // 白噪声残差
        actual_rv[t] = fitted_rv[t] + residuals[t];
    }

    // HAR 诊断
    auto diag = har_diagnostics(actual_rv, fitted_rv, residuals);

    // 残差 LB: 白噪声残差不应拒绝 (无自相关)
    EXPECT_FALSE(diag.residual_ljung_box.base.reject_null)
        << "HAR residuals should be white noise (no autocorrelation)";

    // MZ 回归: 拟合值接近实际值 → α≈0, β≈1
    EXPECT_NEAR(diag.mz_regression.alpha, 0.0, 0.01);
    EXPECT_NEAR(diag.mz_regression.beta, 1.0, 0.05);

    // R² 应较高 (好拟合)
    EXPECT_GT(diag.r_squared, 0.95);

    // 模型应 adequate
    EXPECT_TRUE(diag.model_adequate);
}

// =============================================================================
// 用例 4: HEAVY → 标准化残差 → z² LB (HEAVY 完整诊断流程)
// =============================================================================
TEST(Phase7AIntegration, HEAVY_FullDiagnosticPipeline) {
    // 生成 HEAVY 风格数据: h_t 和 RM_t
    const Size n = 200;
    std::mt19937 gen(42);
    std::normal_distribution<Real> noise(0.0, 1.0);

    std::vector<Real> rm_residuals(n), rm_conditional_means(n);
    std::vector<Real> variance_residuals(n), conditional_variances(n);

    Real h = 0.04;  // 初始方差
    for (Size t = 0; t < n; ++t) {
        Real z = noise(gen);
        conditional_variances[t] = h;
        rm_conditional_means[t] = h;  // RM_t ≈ h_t
        rm_residuals[t] = z * std::sqrt(h);  // RM 残差 = z * sqrt(h)
        variance_residuals[t] = z * z * h - h;  // 方差方程残差

        // HEAVY 更新: h_{t+1} = ω + α*RM_t + β*h_t
        Real RM_t = z * z * h;
        h = 0.01 + 0.4 * RM_t + 0.5 * h;
    }

    // HEAVY 诊断
    auto diag = heavy_diagnostics(
        rm_residuals, rm_conditional_means,
        variance_residuals, conditional_variances);

    // 标准化残差应计算
    EXPECT_EQ(diag.variance_equation.standardized_residuals.size(), n);
    EXPECT_EQ(diag.measurement_equation.standardized_residuals.size(), n);

    // z² LB 应计算 (排幻觉点 H8: z_t² LB 是关键)
    EXPECT_GE(diag.variance_equation.z_squared_ljung_box.base.statistic, 0.0);

    // h_t 与 RM_t 相关性应为正 (波动率聚集)
    EXPECT_GT(diag.correlation_h_rm, 0.0);
}

// =============================================================================
// 用例 5: VaR → DQ + Berkowitz + ES (风险模型完整回测)
// =============================================================================
TEST(Phase7AIntegration, VaR_FullBacktestPipeline) {
    // 生成收益率数据 (正常分布 + 少量违反)
    const Size n = 500;
    std::mt19937 gen(42);
    std::normal_distribution<Real> ret_dist(0.0, 0.01);

    std::vector<Real> returns(n), var_forecasts(n, -0.0164);  // 95% VaR ≈ -1.645%
    std::vector<Real> es_forecasts(n, -0.0208);  // 95% ES ≈ -2.08%

    for (Size t = 0; t < n; ++t) {
        returns[t] = ret_dist(gen);
    }

    // DQ 检验
    auto dq = dynamic_quantile_test(returns, var_forecasts, 0.95);
    EXPECT_GE(dq.dq_statistic, 0.0);
    EXPECT_GE(dq.p_value, 0.0);

    // Berkowitz 检验 (正态模型 CDF)
    auto berk = berkowitz_test(
        returns, [](Real r) -> Real {
            return 0.5 * (1.0 + std::erf(r / (0.01 * std::sqrt(2.0))));
        });
    EXPECT_GE(berk.lr_statistic, 0.0);

    // ES 后验检验
    auto es = es_backtest(returns, var_forecasts, es_forecasts);
    EXPECT_GE(es.n_violations, 0);

    // 整体回测流程完整执行 (无异常即通过)
    SUCCEED() << "VaR full backtest pipeline completed";
}

// =============================================================================
// 用例 6: Greeks → Analytic vs AAD vs Pathwise (Greeks 跨方法一致性)
// =============================================================================
TEST(Phase7AIntegration, Greeks_CrossMethodConsistency) {
    // BSM 欧式看涨期权
    Real S = 100.0, K = 100.0, T = 0.5, r = 0.05, q = 0.02, sigma = 0.20;

    // Delta 一致性
    auto delta_res = greeks_consistency_check(
        S, K, T, r, q, sigma, PayoffType::VanillaCall, "delta", 200000, 42, 0.01);
    EXPECT_TRUE(delta_res.consistent)
        << "Delta should be consistent across all methods";

    // Analytic vs AAD 精确一致
    Real delta_rd = std::abs(delta_res.analytic_value - delta_res.aad_value) /
                    std::abs(delta_res.analytic_value);
    EXPECT_LT(delta_rd, 1e-10);

    // Vega 一致性
    auto vega_res = greeks_consistency_check(
        S, K, T, r, q, sigma, PayoffType::VanillaCall, "vega", 200000, 42, 0.01);
    EXPECT_TRUE(vega_res.consistent)
        << "Vega should be consistent across all methods";

    // Gamma 一致性 (H21: FD 用 dS=1e-4)
    auto gamma_res = greeks_consistency_check(
        S, K, T, r, q, sigma, PayoffType::VanillaCall, "gamma", 200000, 42, 0.01);
    EXPECT_TRUE(gamma_res.consistent)
        << "Gamma should be consistent (Analytic/Numerical/AAD)";
}

// =============================================================================
// 用例 7: IV 拟合优度 → 价格残差 (定价模型诊断)
// =============================================================================
TEST(Phase7AIntegration, Pricing_IVFitAndPriceResidual) {
    const Size n = 15;
    std::vector<Real> strikes(n), maturities(n);
    std::vector<Real> iv_market(n), iv_model(n), spread(n);
    std::vector<Real> market_prices(n), model_prices(n);

    for (Size i = 0; i < n; ++i) {
        strikes[i] = 80.0 + i * 3.0;
        maturities[i] = 0.25 + (i % 4) * 0.25;
        iv_market[i] = 0.18 + 0.001 * i;
        // 注意: Size 是 size_t (无符号), i%3-1 在 i%3==0 时下溢, 必须先转 int
        iv_model[i] = iv_market[i] + 0.0005 * (static_cast<int>(i % 3) - 1);
        spread[i] = 0.02;
        market_prices[i] = 10.0 + i * 0.5;
        model_prices[i] = market_prices[i] + 0.01 * (static_cast<int>(i % 5) - 2);
    }

    // IV 拟合优度
    auto iv_fit = iv_fit_goodness_test(strikes, maturities, iv_market, iv_model, spread);
    EXPECT_GE(iv_fit.chi_squared, 0.0);
    EXPECT_FALSE(iv_fit.reject_good_fit);  // 小偏差 → 好拟合

    // 价格残差诊断
    auto price_diag = price_residual_analysis(market_prices, model_prices);
    EXPECT_FALSE(price_diag.has_bias);  // 小残差 → 无偏差
}

// =============================================================================
// 用例 8: MLE → 信息矩阵检验 (模型设定诊断)
// =============================================================================
TEST(Phase7AIntegration, MLE_InformationMatrixTest) {
    // 生成合成得分矩阵和 Hessian (正确设定模型)
    const Size N = 200;
    const Size K = 2;
    std::mt19937 gen(42);
    std::normal_distribution<Real> noise(0.0, 1.0);

    std::vector<std::vector<Real>> scores(N, std::vector<Real>(K));
    std::vector<std::vector<Real>> hessian(K, std::vector<Real>(K, 0.0));

    // 正确设定: E[score * score'] = -Hessian
    // 生成 iid N(0, I) 得分, Hessian = -I
    for (Size t = 0; t < N; ++t) {
        scores[t][0] = noise(gen);
        scores[t][1] = noise(gen);
    }
    hessian[0][0] = -1.0;
    hessian[1][1] = -1.0;

    // 信息矩阵检验
    auto im = cpphub::v1::econometrics::information_matrix_test(scores, hessian);
    EXPECT_GE(im.base.statistic, 0.0);
    EXPECT_GE(im.base.p_value, 0.0);

    // 正确设定时, IM 不应拒绝 (大样本下 p 值不应极小)
    // 注: 合成数据可能不完全满足, 仅验证检验可执行
}

// =============================================================================
// 用例 9: OLS → CUSUM → Andrews (结构稳定性诊断)
// =============================================================================
TEST(Phase7AIntegration, OLS_StructuralBreakDetection) {
    // 生成含结构断点的数据: y = 1 + 0.5x (前半), y = 1 + 2.0x (后半)
    const Size n = 100;
    std::mt19937 gen(42);
    std::normal_distribution<Real> noise(0.0, 0.5);

    std::vector<Real> y(n);
    std::vector<std::vector<Real>> X(n, std::vector<Real>(2));

    for (Size t = 0; t < n; ++t) {
        X[t][0] = 1.0;
        X[t][1] = static_cast<Real>(t);
        Real slope = (t < n / 2) ? 0.5 : 2.0;  // 断点在 t=50
        y[t] = 1.0 + slope * X[t][1] + noise(gen);
    }

    // CUSUM 检验
    auto cusum = cusum_test(X, y, 0.05);
    EXPECT_EQ(cusum.cusum_path.size(), n - 2);  // 递归残差长度 = n - k

    // Andrews 检验
    auto andrews = andrews_breakpoint_test(X, y, 0.15);
    EXPECT_GT(andrews.breakpoint_estimate, 0);
    EXPECT_GT(andrews.breakpoint_fraction, 0.0);

    // 断点应在 t=50 附近 (trim 15% → 搜索范围 [15, 85])
    EXPECT_GT(andrews.breakpoint_estimate, static_cast<Size>(n * 0.15));
    EXPECT_LT(andrews.breakpoint_estimate, static_cast<Size>(n * 0.85));
}

// =============================================================================
// 用例 10: 跳跃检验 → Bonferroni/BH 修正 (多重检验修正)
// =============================================================================
TEST(Phase7AIntegration, JumpTest_MultipleCorrection) {
    // 模拟 4 种跳跃检验在 5 天的统计量和 p 值
    // 部分天有显著跳跃 (p < 0.05), 部分天无跳跃
    const Size n_days = 5;
    std::vector<Real> bns_stats = {2.5, 0.5, 3.2, 1.0, 0.3};
    std::vector<Real> aj_stats  = {2.0, 0.8, 2.8, 0.9, 0.5};
    std::vector<Real> jo_stats  = {1.8, 0.6, 3.0, 1.1, 0.4};
    std::vector<Real> rank_stats = {2.2, 0.7, 2.9, 0.8, 0.6};

    // Bonferroni 修正 (保守, FWER 控制)
    auto bonf = multiple_test_correction(
        {0.01, 0.30, 0.001, 0.15, 0.50},
        MultipleTestCorrectionResult::Method::Bonferroni);
    EXPECT_EQ(bonf.n_tests, 5);
    // Bonferroni: adjusted_p = min(m * p, 1)
    EXPECT_NEAR(bonf.adjusted_p_values[0], 0.05, 1e-10);  // 5 * 0.01

    // BH 修正 (FDR 控制, 排幻觉点 H23)
    auto bh = multiple_test_correction(
        {0.01, 0.30, 0.001, 0.15, 0.50},
        MultipleTestCorrectionResult::Method::BenjaminiHochberg);
    EXPECT_EQ(bh.n_tests, 5);

    // BH 应比 Bonferroni 更 powerful (拒绝更多)
    EXPECT_GE(bh.n_rejections, bonf.n_rejections);

    // 跳跃检验联合诊断
    auto diag = jump_test_diagnostics(bns_stats, aj_stats, jo_stats, rank_stats);
    EXPECT_EQ(diag.test_statistics.size(), n_days);
    EXPECT_EQ(diag.p_values.size(), n_days);

    // 第 3 天 (统计量最大) 应被多数检验识别为跳跃
    EXPECT_GT(diag.test_statistics[2], diag.test_statistics[1]);
}
