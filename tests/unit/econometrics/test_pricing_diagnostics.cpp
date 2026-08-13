// =============================================================================
// test_pricing_diagnostics.cpp - Phase 7A Wave 3c 定价模型诊断测试
//
// 15 用例:
//   IV 拟合优度检验 (8) + 价格残差诊断 (7)
//
// 排幻觉点覆盖:
//   H22 (权重用市场 IV 的 Bid-Ask 宽度, 非简单方差)
//
// 教材锚点: Gatheral 2006 (SVI), Fengler 2009 (IV surface)
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

#include "cpphub/pricing/pricing_diagnostics.hpp"

using cpphub::v1::iv_fit_goodness_test;
using cpphub::v1::price_residual_analysis;
using cpphub::v1::IVFitGoodnessResult;
using cpphub::v1::PriceResidualDiagnostics;
using cpphub::v1::Real;
using cpphub::v1::Size;

// =============================================================================
// IV 拟合优度检验 (8 用例)
// =============================================================================

// --- IV 1: 好拟合 (模型 IV ≈ 市场 IV) → χ² 不拒绝 ---
TEST(IVFitGoodness, GoodFitNotRejected) {
    const Size n = 10;
    std::vector<Real> strikes(n), maturities(n), iv_market(n), iv_model(n), spread(n);
    for (Size i = 0; i < n; ++i) {
        strikes[i] = 90.0 + i * 2.0;
        maturities[i] = 0.25 + i * 0.05;
        iv_market[i] = 0.20 + 0.001 * i;
        iv_model[i] = iv_market[i] + 0.0001;  // 极小偏差
        spread[i] = 0.02;  // 2% bid-ask 宽度
    }
    auto res = iv_fit_goodness_test(strikes, maturities, iv_market, iv_model, spread);
    EXPECT_FALSE(res.reject_good_fit);
    EXPECT_GT(res.p_value, 0.05);
}

// --- IV 2: 差拟合 (模型 IV 偏离市场 IV) → χ² 拒绝 ---
TEST(IVFitGoodness, BadFitRejected) {
    const Size n = 10;
    std::vector<Real> strikes(n), maturities(n), iv_market(n), iv_model(n), spread(n);
    for (Size i = 0; i < n; ++i) {
        strikes[i] = 90.0 + i * 2.0;
        maturities[i] = 0.25 + i * 0.05;
        iv_market[i] = 0.20;
        iv_model[i] = 0.25;  // 大偏差 (5%)
        spread[i] = 0.02;
    }
    auto res = iv_fit_goodness_test(strikes, maturities, iv_market, iv_model, spread);
    EXPECT_TRUE(res.reject_good_fit);
    EXPECT_LT(res.p_value, 0.05);
}

// --- IV 3: 排幻觉点 H22 — 权重用 Bid-Ask 宽度 ---
// 宽度大 → 权重小 → 同样偏差 χ² 更小 (更不易拒绝)
TEST(IVFitGoodness, WeightByBidAskSpreadH22) {
    const Size n = 5;
    std::vector<Real> strikes(n), maturities(n), iv_market(n), iv_model(n);
    for (Size i = 0; i < n; ++i) {
        strikes[i] = 100.0 + i;
        maturities[i] = 0.25;
        iv_market[i] = 0.20;
        iv_model[i] = 0.21;  // 固定偏差 1%
    }

    // 小宽度 → 权重大 → χ² 大
    std::vector<Real> small_spread(n, 0.01);
    auto res_small = iv_fit_goodness_test(strikes, maturities, iv_market, iv_model, small_spread);

    // 大宽度 → 权重小 → χ² 小
    std::vector<Real> large_spread(n, 0.10);
    auto res_large = iv_fit_goodness_test(strikes, maturities, iv_market, iv_model, large_spread);

    EXPECT_GT(res_small.chi_squared, res_large.chi_squared);
    // 小宽度更可能拒绝
    EXPECT_GE(res_small.p_value, 0.0);
    EXPECT_LE(res_large.p_value, 1.0);
}

// --- IV 4: 自由度 = N - 1 ---
TEST(IVFitGoodness, DegreesOfFreedomCorrect) {
    const Size n = 15;
    std::vector<Real> strikes(n), maturities(n), iv_market(n), iv_model(n), spread(n, 0.02);
    for (Size i = 0; i < n; ++i) {
        strikes[i] = 90.0 + i;
        maturities[i] = 0.25;
        iv_market[i] = 0.20;
        iv_model[i] = 0.20;
    }
    auto res = iv_fit_goodness_test(strikes, maturities, iv_market, iv_model, spread);
    EXPECT_EQ(res.degrees_of_freedom, n - 1);
}

// --- IV 5: 残差向量正确 ---
TEST(IVFitGoodness, ResidualsCorrect) {
    const Size n = 5;
    std::vector<Real> strikes(n, 100.0), maturities(n, 0.25),
        iv_market = {0.19, 0.20, 0.21, 0.22, 0.23},
        iv_model  = {0.20, 0.20, 0.20, 0.20, 0.20},
        spread(n, 0.02);
    auto res = iv_fit_goodness_test(strikes, maturities, iv_market, iv_model, spread);

    for (Size i = 0; i < n; ++i) {
        EXPECT_NEAR(res.residuals[i], iv_model[i] - iv_market[i], 1e-10);
    }
}

// --- IV 6: RMSE 计算 ---
TEST(IVFitGoodness, RMSECalculated) {
    const Size n = 5;
    std::vector<Real> strikes(n, 100.0), maturities(n, 0.25),
        iv_market(n, 0.20), iv_model = {0.21, 0.22, 0.23, 0.24, 0.25},
        spread(n, 0.02);
    auto res = iv_fit_goodness_test(strikes, maturities, iv_market, iv_model, spread);

    // 残差 = [0.01, 0.02, 0.03, 0.04, 0.05]
    // RMSE = sqrt((0.01²+0.02²+0.03²+0.04²+0.05²)/5)
    //       = sqrt((1+4+9+16+25)*1e-4/5) = sqrt(55e-4/5) = sqrt(1.1e-3) ≈ 0.0332
    EXPECT_NEAR(res.rmse, std::sqrt(0.0011), 1e-6);
}

// --- IV 7: 最大绝对残差 ---
TEST(IVFitGoodness, MaxAbsResidualCorrect) {
    const Size n = 5;
    std::vector<Real> strikes(n, 100.0), maturities(n, 0.25),
        iv_market = {0.20, 0.20, 0.20, 0.20, 0.20},
        iv_model  = {0.21, 0.25, 0.22, 0.23, 0.20},  // 最大偏差 0.05
        spread(n, 0.02);
    auto res = iv_fit_goodness_test(strikes, maturities, iv_market, iv_model, spread);
    EXPECT_NEAR(res.max_abs_residual, 0.05, 1e-10);
}

// --- IV 8: 异常 — spread <= 0 抛异常 ---
TEST(IVFitGoodness, NonPositiveSpreadThrows) {
    const Size n = 5;
    std::vector<Real> strikes(n, 100.0), maturities(n, 0.25),
        iv_market(n, 0.20), iv_model(n, 0.20),
        spread(n, 0.0);  // 零宽度
    EXPECT_THROW(iv_fit_goodness_test(strikes, maturities, iv_market, iv_model, spread),
                 std::invalid_argument);
}

// =============================================================================
// 价格残差诊断 (7 用例)
// =============================================================================

// --- Price 1: 无偏差 (残差均值≈0) → has_bias=false ---
TEST(PriceResidual, NoBiasNotRejected) {
    const Size n = 20;
    std::vector<Real> market(n, 10.0), model(n);
    for (Size i = 0; i < n; ++i) {
        // 残差 = ±0.01 交替, 均值 = 0
        model[i] = 10.0 + (i % 2 == 0 ? 0.01 : -0.01);
    }
    auto res = price_residual_analysis(market, model);
    EXPECT_FALSE(res.has_bias);
    EXPECT_NEAR(res.mean_residual, 0.0, 1e-10);
}

// --- Price 2: 有偏差 (模型系统性高估) → has_bias=true ---
TEST(PriceResidual, BiasDetected) {
    const Size n = 50;
    std::mt19937 gen(42);
    std::normal_distribution<Real> dist(0.0, 0.01);
    std::vector<Real> market(n, 10.0), model(n);
    for (Size i = 0; i < n; ++i) {
        model[i] = 10.5 + dist(gen);  // 系统性高估 0.5 + 小噪声
    }
    auto res = price_residual_analysis(market, model);
    EXPECT_TRUE(res.has_bias);
    EXPECT_NEAR(res.mean_residual, 0.5, 0.01);
    EXPECT_LT(res.p_value_bias, 0.05);
}

// --- Price 3: 残差方向 (正=高估, 负=低估) ---
TEST(PriceResidual, ResidualDirection) {
    const Size n = 10;
    // 模型高估
    std::vector<Real> market(n, 10.0), model_high(n, 11.0);
    auto res_high = price_residual_analysis(market, model_high);
    EXPECT_GT(res_high.mean_residual, 0.0);

    // 模型低估
    std::vector<Real> model_low(n, 9.0);
    auto res_low = price_residual_analysis(market, model_low);
    EXPECT_LT(res_low.mean_residual, 0.0);
}

// --- Price 4: t 统计量计算 ---
TEST(PriceResidual, TStatisticCorrect) {
    const Size n = 100;
    std::vector<Real> market(n, 10.0), model(n);
    // 残差 = 0.1 (固定), std=0 → t 应为 0 (或处理为无偏)
    for (Size i = 0; i < n; ++i) model[i] = 10.1;
    auto res = price_residual_analysis(market, model);
    // std=0 时 t=0, p=1 (无法拒绝)
    EXPECT_NEAR(res.t_stat_bias, 0.0, 1e-10);
    EXPECT_NEAR(res.p_value_bias, 1.0, 1e-10);
    EXPECT_FALSE(res.has_bias);
}

// --- Price 5: 标准差计算 ---
TEST(PriceResidual, StdResidualCorrect) {
    const Size n = 5;
    std::vector<Real> market(n, 10.0), model = {10.1, 10.2, 10.3, 10.4, 10.5};
    auto res = price_residual_analysis(market, model);
    // 残差 = [0.1, 0.2, 0.3, 0.4, 0.5], mean=0.3
    // std = sqrt(((0.1-0.3)² + (0.2-0.3)² + (0.3-0.3)² + (0.4-0.3)² + (0.5-0.3)²) / 4)
    //     = sqrt((0.04+0.01+0+0.01+0.04)/4) = sqrt(0.1/4) = sqrt(0.025) ≈ 0.1581
    EXPECT_NEAR(res.mean_residual, 0.3, 1e-10);
    EXPECT_NEAR(res.std_residual, std::sqrt(0.025), 1e-6);
}

// --- Price 6: 残差向量正确 ---
TEST(PriceResidual, ResidualsVector) {
    const Size n = 5;
    std::vector<Real> market = {10.0, 11.0, 12.0, 13.0, 14.0};
    std::vector<Real> model  = {10.5, 11.5, 12.5, 13.5, 14.5};
    auto res = price_residual_analysis(market, model);
    for (Size i = 0; i < n; ++i) {
        EXPECT_NEAR(res.residuals[i], model[i] - market[i], 1e-10);
    }
}

// --- Price 7: 异常 — 价格太少抛异常 ---
TEST(PriceResidual, TooFewPricesThrows) {
    std::vector<Real> market(2, 10.0), model(2, 10.0);
    EXPECT_THROW(price_residual_analysis(market, model), std::invalid_argument);
}

// --- Price 8: 异常 — 大小不匹配抛异常 ---
TEST(PriceResidual, SizeMismatchThrows) {
    std::vector<Real> market(10, 10.0), model(5, 10.0);
    EXPECT_THROW(price_residual_analysis(market, model), std::invalid_argument);
}
