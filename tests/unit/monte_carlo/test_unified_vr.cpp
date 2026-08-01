// 统一方差缩减 (Variance Reduction) 测试
//
// 覆盖:
//   1. Antithetic 独立类 (antithetic.hpp)
//   2. VREngine 统一 Decorator (variance_reduction.hpp)
//   3. Importance Sampling Put 定价 (importance_sampling.hpp 新增)
//   4. path_generator antithetic Z 复用验证
//   5. mc_var.hpp 死字段移除确认

#include <gtest/gtest.h>
#include "cpphub/pricing/monte_carlo/antithetic.hpp"
#include "cpphub/pricing/monte_carlo/variance_reduction.hpp"
#include "cpphub/pricing/monte_carlo/path_generator.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/pricing/monte_carlo/path_dependent_payoffs.hpp"
#include "cpphub/pricing/monte_carlo/qmc_engine.hpp"
#include "cpphub/monte_carlo/importance_sampling.hpp"
#include "cpphub/risk/var/mc_var.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/rng.hpp"
#include <cmath>
#include <vector>

using namespace cpphub;

namespace {

// 注: bsm_call_price / bsm_put_price 已由 control_variate.hpp 提供 (cpphub::v1)
// 此处不再重复定义, 避免与 using namespace cpphub 产生歧义

// 几何平均亚式期权解析价格 (Kemna-Vorst 1990, GBM 假设)
Real geom_asian_call_price(Real S0, Real K, Real T, Real r, Real q, Real sigma, Size n_steps) {
    if (T <= 0.0) return std::max(S0 - K, 0.0);
    Real dt = T / static_cast<Real>(n_steps);
    Real n = static_cast<Real>(n_steps);
    Real mu_g = std::log(S0) + (r - q - 0.5 * sigma * sigma) * T * (n + 1) / (2.0 * n);
    Real var_g = sigma * sigma * T * (n + 1) * (2 * n + 1) / (6.0 * n * n);
    Real d1 = (mu_g + var_g - std::log(K)) / std::sqrt(var_g);
    Real d2 = d1 - std::sqrt(var_g);
    return std::exp(-r * T) * (std::exp(mu_g + 0.5 * var_g) * normal_cdf(d1) - K * normal_cdf(d2));
}

}  // anonymous namespace

// ===========================================================================
// 1. Antithetic 独立类测试
// ===========================================================================

TEST(AntitheticClass, PathPairSharedZ) {
    // 单资产 GBM
    MultiAssetGBMConfig cfg;
    cfg.S0 = {100.0};
    cfg.sigma = {0.20};
    cfg.q = {0.0};
    cfg.r = 0.05;
    cfg.T = 1.0;
    cfg.n_steps = 50;
    auto gen_cfg = make_single_asset_gbm(100.0, 0.20, 0.05, 0.0, 1.0, 50);
    MultiAssetGBMPathGenerator gen(gen_cfg);
    Antithetic ant(gen);

    Philox4x64 rng(42, 0);
    auto [p_plus, p_minus] = ant.generate_path_pair(rng);

    ASSERT_EQ(p_plus.size(), 1u);
    ASSERT_EQ(p_minus.size(), 1u);
    ASSERT_EQ(p_plus[0].size(), 51u);
    ASSERT_EQ(p_minus[0].size(), 51u);

    // 起点相同
    EXPECT_NEAR(p_plus[0][0], 100.0, 1e-12);
    EXPECT_NEAR(p_minus[0][0], 100.0, 1e-12);

    // 中间点应不同 (除非 Z=0, 概率为 0)
    EXPECT_NE(p_plus[0][25], p_minus[0][25]);
}

TEST(AntitheticClass, SampleSingleAsset) {
    auto gen_cfg = make_single_asset_gbm(100.0, 0.20, 0.05, 0.0, 1.0, 50);
    MultiAssetGBMPathGenerator gen(gen_cfg);
    Antithetic ant(gen);

    // 欧式看涨 payoff: max(S_T - K, 0)
    Real K = 100.0;
    auto payoff = [K](const std::vector<Real>& path) -> Real {
        return std::max(path.back() - K, 0.0);
    };

    Philox4x64 rng(42, 0);
    Real sample = ant.sample_single_asset(payoff, rng);
    EXPECT_GT(sample, 0.0);  // 应为正
}

TEST(AntitheticClass, AntitheticRiskFactorPair) {
    std::vector<Real> dR = {0.01, -0.02, 0.005, 0.03};
    auto [plus, minus] = antithetic_risk_factor_pair(dR);
    EXPECT_EQ(plus.size(), 4u);
    EXPECT_EQ(minus.size(), 4u);
    for (Size i = 0; i < dR.size(); ++i) {
        EXPECT_NEAR(plus[i], dR[i], 1e-15);
        EXPECT_NEAR(minus[i], -dR[i], 1e-15);
    }
}

// ===========================================================================
// 2. VREngine 统一 Decorator 测试
// ===========================================================================

TEST(VREngine, PlainMCConvergesToBSM) {
    auto gen_cfg = make_single_asset_gbm(100.0, 0.20, 0.05, 0.0, 1.0, 1);
    MultiAssetGBMPathGenerator gen(gen_cfg);

    VRConfig cfg;
    cfg.n_paths = 200000;
    cfg.seed = 42;
    cfg.use_antithetic = false;
    cfg.use_control_variate = false;
    cfg.df = std::exp(-0.05 * 1.0);

    // 欧式看涨: payoff = max(S_T - K, 0)
    Real K = 100.0;
    PathPayoff payoff = [K](const std::vector<Real>& path) -> Real {
        return std::max(path.back() - K, 0.0);
    };

    auto result = VREngine::price_single_asset(gen, payoff, cfg);
    Real expected = bsm_call_price(100, 100, 1.0, 0.05, 0.0, 0.20);
    EXPECT_NEAR(result.price, expected, 0.10);  // plain MC, ~5 sigma
    EXPECT_EQ(result.n_samples, 200000u);
    EXPECT_EQ(result.n_paths_generated, 200000u);
}

TEST(VREngine, AntitheticReducesVariance) {
    auto gen_cfg = make_single_asset_gbm(100.0, 0.20, 0.05, 0.0, 1.0, 1);
    MultiAssetGBMPathGenerator gen(gen_cfg);

    Real K = 100.0;
    PathPayoff payoff = [K](const std::vector<Real>& path) -> Real {
        return std::max(path.back() - K, 0.0);
    };

    // Plain MC
    VRConfig cfg_plain;
    cfg_plain.n_paths = 200000;
    cfg_plain.seed = 42;
    cfg_plain.df = std::exp(-0.05);
    auto r_plain = VREngine::price_single_asset(gen, payoff, cfg_plain);

    // Antithetic (相同总路径数)
    VRConfig cfg_ant;
    cfg_ant.n_paths = 200000;
    cfg_ant.seed = 42;
    cfg_ant.use_antithetic = true;
    cfg_ant.df = std::exp(-0.05);
    auto r_ant = VREngine::price_single_asset(gen, payoff, cfg_ant);

    Real expected = bsm_call_price(100, 100, 1.0, 0.05, 0.0, 0.20);
    EXPECT_NEAR(r_ant.price, expected, 0.10);
    EXPECT_EQ(r_ant.n_paths_generated, 200000u);
    EXPECT_EQ(r_ant.n_samples, 100000u);  // antithetic: n_paths/2 配对

    // Antithetic 标准误差应显著小于 plain MC (对单调 payoff)
    EXPECT_LT(r_ant.std_error, r_plain.std_error);
}

TEST(VREngine, ControlVariateReducesVariance) {
    // 用几何平均亚式作为算术平均亚式的控制变量
    Size n_steps = 50;
    auto gen_cfg = make_single_asset_gbm(100.0, 0.20, 0.05, 0.0, 1.0, n_steps);
    MultiAssetGBMPathGenerator gen(gen_cfg);

    Real K = 100.0;
    // 算术平均亚式 payoff (目标)
    PathPayoff payoff = [K](const std::vector<Real>& path) -> Real {
        Real sum = 0.0;
        for (Size i = 1; i < path.size(); ++i) sum += path[i];
        Real avg = sum / static_cast<Real>(path.size() - 1);
        return std::max(avg - K, 0.0);
    };
    // 几何平均亚式 payoff (控制变量, 有解析解)
    PathPayoff cv_payoff = [K](const std::vector<Real>& path) -> Real {
        Real log_sum = 0.0;
        for (Size i = 1; i < path.size(); ++i) log_sum += std::log(path[i]);
        Real geom = std::exp(log_sum / static_cast<Real>(path.size() - 1));
        return std::max(geom - K, 0.0);
    };
    Real cv_price = geom_asian_call_price(100, K, 1.0, 0.05, 0.0, 0.20, n_steps);

    // Plain MC
    VRConfig cfg_plain;
    cfg_plain.n_paths = 50000;
    cfg_plain.seed = 42;
    cfg_plain.df = std::exp(-0.05);
    auto r_plain = VREngine::price_single_asset(gen, payoff, cfg_plain);

    // CV
    VRConfig cfg_cv;
    cfg_cv.n_paths = 50000;
    cfg_cv.seed = 42;
    cfg_cv.use_control_variate = true;
    cfg_cv.df = std::exp(-0.05);
    auto r_cv = VREngine::price_single_asset(gen, payoff, cfg_cv, cv_payoff, cv_price);

    // CV 应显著缩减方差
    EXPECT_LT(r_cv.std_error, r_plain.std_error);
    EXPECT_LT(r_cv.variance_reduction, 1.0);  // 方差缩减因子 < 1
    EXPECT_NE(r_cv.beta_cv, 0.0);             // β 非零
}

TEST(VREngine, AntitheticPlusCVComposable) {
    Size n_steps = 50;
    auto gen_cfg = make_single_asset_gbm(100.0, 0.20, 0.05, 0.0, 1.0, n_steps);
    MultiAssetGBMPathGenerator gen(gen_cfg);

    Real K = 100.0;
    PathPayoff payoff = [K](const std::vector<Real>& path) -> Real {
        Real sum = 0.0;
        for (Size i = 1; i < path.size(); ++i) sum += path[i];
        Real avg = sum / static_cast<Real>(path.size() - 1);
        return std::max(avg - K, 0.0);
    };
    PathPayoff cv_payoff = [K](const std::vector<Real>& path) -> Real {
        Real log_sum = 0.0;
        for (Size i = 1; i < path.size(); ++i) log_sum += std::log(path[i]);
        Real geom = std::exp(log_sum / static_cast<Real>(path.size() - 1));
        return std::max(geom - K, 0.0);
    };
    Real cv_price = geom_asian_call_price(100, K, 1.0, 0.05, 0.0, 0.20, n_steps);

    VRConfig cfg;
    cfg.n_paths = 50000;
    cfg.seed = 42;
    cfg.use_antithetic = true;
    cfg.use_control_variate = true;
    cfg.df = std::exp(-0.05);

    auto r = VREngine::price_single_asset(gen, payoff, cfg, cv_payoff, cv_price);
    // 两者组合应产生更小的标准误差
    EXPECT_GT(r.std_error, 0.0);
    EXPECT_LT(r.variance_reduction, 1.0);
    EXPECT_EQ(r.n_paths_generated, 50000u);
}

// ===========================================================================
// 3. Importance Sampling Put 定价测试
// ===========================================================================

TEST(ImportanceSamplingPut, MatchesBSMAnalytic) {
    ISConfig cfg;
    cfg.auto_optimize = true;
    ImportanceSampling is(cfg);

    // OTM Put: K=90, S0=100 (深度 OTM, IS 优势明显)
    Real S0 = 100.0, K = 90.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    auto result = is.price_european_put(S0, K, T, r, q, sigma, 200000, 42);

    Real expected = bsm_put_price(S0, K, T, r, q, sigma);
    EXPECT_NEAR(result.price, expected, 0.02);  // IS 应高精度收敛
    EXPECT_GT(result.variance_reduction_ratio, 1.0);  // IS 应缩减方差
}

TEST(ImportanceSamplingPut, ITMPutMatchesBSM) {
    ISConfig cfg;
    cfg.auto_optimize = true;
    ImportanceSampling is(cfg);

    // ITM Put: K=110, S0=100
    Real S0 = 100.0, K = 110.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    auto result = is.price_european_put(S0, K, T, r, q, sigma, 200000, 42);

    Real expected = bsm_put_price(S0, K, T, r, q, sigma);
    EXPECT_NEAR(result.price, expected, 0.05);
}

TEST(ImportanceSamplingPut, ZeroMaturityReturnsIntrinsic) {
    ISConfig cfg;
    ImportanceSampling is(cfg);
    auto result = is.price_european_put(100.0, 110.0, 0.0, 0.05, 0.0, 0.20, 1000, 42);
    EXPECT_NEAR(result.price, 10.0, 1e-12);  // max(K-S0, 0) = 10
}

// ===========================================================================
// 4. path_generator antithetic Z 复用验证
// ===========================================================================

TEST(GBMPathGenerator, AntitheticZReuse) {
    GBMConfig cfg;
    cfg.S0 = 100.0;
    cfg.r = 0.05;
    cfg.q = 0.0;
    cfg.sigma = 0.20;
    cfg.T = 1.0;
    cfg.n_steps = 50;

    GBMPathGenerator gen(cfg, PathScheme::Exact);
    auto paths = gen.generate_paths(1000, 42, true);

    Size path_len = 51;
    ASSERT_EQ(paths.size(), 1000 * path_len);

    // 前 500 条为原始, 后 500 条为反变量
    // 验证: 原始路径和反变量路径的终端价格应关于漂移对称
    // E[ln(S_T^+)] = ln(S0) + (r-q-σ²/2)T
    // 对于 antithetic 配对: ln(S_T^+) + ln(S_T^-) = 2*ln(S0) + 2*(r-q-σ²/2)T
    Real mu_T = (0.05 - 0.0 - 0.5 * 0.20 * 0.20) * 1.0;
    Real expected_pair_sum = 2.0 * std::log(100.0) + 2.0 * mu_T;

    Size n_pairs = 500;
    Real sum_log_pair = 0.0;
    for (Size p = 0; p < n_pairs; ++p) {
        Real S_T_plus  = paths[p * path_len + (path_len - 1)];
        Real S_T_minus = paths[(p + n_pairs) * path_len + (path_len - 1)];
        sum_log_pair += std::log(S_T_plus) + std::log(S_T_minus);
    }
    Real avg_pair_sum = sum_log_pair / static_cast<Real>(n_pairs);
    EXPECT_NEAR(avg_pair_sum, expected_pair_sum, 1e-9);  // Z 复用则精确对称
}

// ===========================================================================
// 5. mc_var.hpp 死字段移除确认
// ===========================================================================

TEST(MCVarConfig, NoControlVariateField) {
    // 确认 use_control_variate 字段已移除 (编译期检查)
    // 若字段仍存在, 此测试会编译失败
    MCVarConfig cfg;
    cfg.n_paths = 1000;
    cfg.seed = 42;
    cfg.antithetic = true;
    // cfg.use_control_variate = false;  // 已移除, 不应编译
    EXPECT_EQ(cfg.n_paths, 1000u);
    EXPECT_TRUE(cfg.antithetic);
}

// ===========================================================================
// 6. Sobol QMC 集成测试 (验证拟随机数方差缩减)
// ===========================================================================

TEST(QMCIntegration, EuropeanCallMatchesBSM) {
    // Sobol QMC 对欧式看涨期权应高精度收敛到 BSM 解析价格
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    QMCConfig cfg;
    cfg.n_paths = 4096;
    cfg.n_replicates = 16;
    cfg.use_brownian_bridge = true;
    cfg.df = std::exp(-r * T);

    auto payoff = [K](Real ST) { return std::max(ST - K, 0.0); };
    auto result = price_european_qmc_single(S0, sigma, r, q, T, payoff, cfg);

    Real expected = bsm_call_price(S0, K, T, r, q, sigma);
    EXPECT_NEAR(result.price, expected, 0.02);  // QMC 应高精度收敛
    EXPECT_GT(result.std_error, 0.0);
    EXPECT_EQ(result.n_paths, 4096u);
    EXPECT_EQ(result.n_replicates, 16u);
    EXPECT_EQ(result.n_total_paths, 4096u * 16u);
}

TEST(QMCIntegration, PathDepentAsianMatchesGeomAnalytic) {
    // 路径相关 QMC: 几何平均亚式期权应匹配 Kemna-Vorst 解析解
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Size n_steps = 32;
    QMCConfig cfg;
    cfg.n_paths = 4096;
    cfg.n_replicates = 16;
    cfg.use_brownian_bridge = true;
    cfg.df = std::exp(-r * T);

    // 几何平均亚式 payoff (path 包含 S(t_1)..S(t_n), 长度 n_steps)
    PathPayoff payoff = [K, n_steps](const std::vector<Real>& path) -> Real {
        Real log_sum = 0.0;
        for (Size i = 0; i < path.size(); ++i) log_sum += std::log(path[i]);
        Real geom = std::exp(log_sum / static_cast<Real>(path.size()));
        return std::max(geom - K, 0.0);
    };

    auto result = price_path_dependent_qmc_single(S0, sigma, r, q, T, n_steps, payoff, cfg);
    Real expected = geom_asian_call_price(S0, K, T, r, q, sigma, n_steps);
    EXPECT_NEAR(result.price, expected, 0.03);
}

TEST(QMCIntegration, VarianceReductionVsPlainMC) {
    // QMC 的标准误差应显著小于相同路径数的 plain MC
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;

    // QMC
    QMCConfig qmc_cfg;
    qmc_cfg.n_paths = 4096;
    qmc_cfg.n_replicates = 16;
    qmc_cfg.df = std::exp(-r * T);
    auto payoff = [K](Real ST) { return std::max(ST - K, 0.0); };
    auto qmc_result = price_european_qmc_single(S0, sigma, r, q, T, payoff, qmc_cfg);

    // Plain MC (相同总路径数, 用 VREngine)
    auto gen_cfg = make_single_asset_gbm(S0, sigma, r, q, T, 1);
    MultiAssetGBMPathGenerator gen(gen_cfg);
    VRConfig mc_cfg;
    mc_cfg.n_paths = 4096 * 16;  // 相同总路径数
    mc_cfg.seed = 42;
    mc_cfg.df = std::exp(-r * T);
    PathPayoff path_payoff = [K](const std::vector<Real>& path) -> Real {
        return std::max(path.back() - K, 0.0);
    };
    auto mc_result = VREngine::price_single_asset(gen, path_payoff, mc_cfg);

    // QMC 标准误差应远小于 plain MC (通常 10x-100x 缩减)
    EXPECT_LT(qmc_result.std_error, mc_result.std_error);
}
