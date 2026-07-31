// SOURCE: Glasserman (2003) Ch.4 — 高级方差缩减方法测试
// 测试模块: Importance Sampling / Conditional MC / Stratified Sampling
//
// 容差策略: MC 测试用 4σ 或 5σ 容差 (99.99% CI); 方差缩减比 > 1.0 即可
#include <gtest/gtest.h>
#include "cpphub/monte_carlo/importance_sampling.hpp"
#include "cpphub/monte_carlo/conditional_mc.hpp"
#include "cpphub/monte_carlo/stratified_sampling.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"
#include "cpphub/pricing/monte_carlo/sobol.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/rng.hpp"
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace cpphub;

namespace {
// 辅助: 从 Philox4x64 生成一个 N(0,1) 样本
Real generate_normal(Philox4x64& rng) {
    uint64_t r1 = rng();
    uint64_t r2 = rng();
    Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
    Real u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
    return box_muller(u1, u2).first;
}

// 标准 MC 欧式 Call 定价 (返回价格和标准误差)
struct MCPrice { Real price; Real std_error; Real variance; };
MCPrice standard_mc_european_call(Real S0, Real K, Real T, Real r, Real q,
                                    Real sigma, Size n_paths, uint64_t seed) {
    Real mu_T = (r - q - 0.5 * sigma * sigma) * T;
    Real sqrtT = std::sqrt(T);
    Real df = std::exp(-r * T);
    Real sum = 0.0, sum_sq = 0.0;
    for (Size i = 0; i < n_paths; ++i) {
        Philox4x64 rng(seed, static_cast<uint64_t>(i) + 1);
        Real Z = generate_normal(rng);
        Real S_T = S0 * std::exp(mu_T + sigma * sqrtT * Z);
        Real payoff = std::max(S_T - K, 0.0);
        sum += payoff;
        sum_sq += payoff * payoff;
    }
    Real inv_n = 1.0 / static_cast<Real>(n_paths);
    Real mean = sum * inv_n;
    Real var = sum_sq * inv_n - mean * mean;
    if (var < 0.0) var = 0.0;
    return { df * mean, df * std::sqrt(var * inv_n), var };
}

// 标准 MC 障碍期权 (Up-and-Out Call) — 路径模拟, 用 n_steps 个时间步
MCPrice standard_mc_barrier_uoc(Real S0, Real K, Real B, Real T, Real r, Real q,
                                  Real sigma, Size n_steps, Size n_paths, uint64_t seed) {
    Real mu = r - q - 0.5 * sigma * sigma;
    Real dt = T / static_cast<Real>(n_steps);
    Real df = std::exp(-r * T);
    Real sum = 0.0, sum_sq = 0.0;
    for (Size p = 0; p < n_paths; ++p) {
        Philox4x64 rng(seed, static_cast<uint64_t>(p) + 1);
        Real S = S0;
        bool knocked_out = false;
        for (Size i = 0; i < n_steps; ++i) {
            Real Z = generate_normal(rng);
            S = S * std::exp(mu * dt + sigma * std::sqrt(dt) * Z);
            if (S >= B) { knocked_out = true; break; }
        }
        Real payoff = knocked_out ? 0.0 : std::max(S - K, 0.0);
        sum += payoff;
        sum_sq += payoff * payoff;
    }
    Real inv_n = 1.0 / static_cast<Real>(n_paths);
    Real mean = sum * inv_n;
    Real var = sum_sq * inv_n - mean * mean;
    if (var < 0.0) var = 0.0;
    return { df * mean, df * std::sqrt(var * inv_n), var };
}

// 标准 MC 亚式算术平均 Call — 路径模拟
MCPrice standard_mc_asian_arith_call(Real S0, Real K, Real T, Real r, Real q,
                                       Real sigma, Size n_steps, Size n_paths, uint64_t seed) {
    Real mu = r - q - 0.5 * sigma * sigma;
    Real dt = T / static_cast<Real>(n_steps);
    Real df = std::exp(-r * T);
    Real inv_n_steps = 1.0 / static_cast<Real>(n_steps);
    Real sum = 0.0, sum_sq = 0.0;
    for (Size p = 0; p < n_paths; ++p) {
        Philox4x64 rng(seed, static_cast<uint64_t>(p) + 1);
        Real S = S0;
        Real sum_S = 0.0;
        for (Size i = 0; i < n_steps; ++i) {
            Real Z = generate_normal(rng);
            S = S * std::exp(mu * dt + sigma * std::sqrt(dt) * Z);
            sum_S += S;
        }
        Real avg = sum_S * inv_n_steps;
        Real payoff = std::max(avg - K, 0.0);
        sum += payoff;
        sum_sq += payoff * payoff;
    }
    Real inv_n = 1.0 / static_cast<Real>(n_paths);
    Real mean = sum * inv_n;
    Real var = sum_sq * inv_n - mean * mean;
    if (var < 0.0) var = 0.0;
    return { df * mean, df * std::sqrt(var * inv_n), var };
}
} // anonymous namespace

// ============================================================
// ============ Importance Sampling 测试 (6 用例) ============
// ============================================================

// 1. 最优 theta 自动计算: theta* = (ln(K/S0) - (r-q-sigma^2/2)T) / (sigma*sqrt(T)) = -d2
TEST(ImportanceSampling, OptimalThetaCall) {
    Real S0 = 100, K = 120, T = 1.0, r = 0.05, q = 0.02, sigma = 0.3;
    Real theta = ImportanceSampling::optimal_theta_call(S0, K, T, r, q, sigma);
    // 验证: theta* = -d2 (BSM d2)
    Real d1 = (std::log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    Real d2 = d1 - sigma * std::sqrt(T);
    Real expected_theta = -d2;
    EXPECT_NEAR(theta, expected_theta, 1e-12);
    // OTM Call (K > S0) 应有 theta > 0 (推高 S_T)
    EXPECT_GT(theta, 0.0);
}

// 2. 最优 theta Put: OTM Put (K < S0) 应有 theta < 0 (推低 S_T)
TEST(ImportanceSampling, OptimalThetaPut) {
    Real S0 = 100, K = 80, T = 1.0, r = 0.05, q = 0.0, sigma = 0.25;
    Real theta = ImportanceSampling::optimal_theta_put(S0, K, T, r, q, sigma);
    EXPECT_NEAR(theta, ImportanceSampling::optimal_theta_call(S0, K, T, r, q, sigma), 1e-15);
    // OTM Put (K < S0): theta < 0
    EXPECT_LT(theta, 0.0);
}

// 3. 似然比正确性: L = exp(theta*Z + 0.5*theta^2)
TEST(ImportanceSampling, LikelihoodRatioCorrectness) {
    ISConfig cfg;
    cfg.theta = 0.5;
    cfg.auto_optimize = false;
    ImportanceSampling is(cfg);
    Real T = 1.0;
    Real Z = 1.23;
    Real L = 0.0;
    Real W_T = is.sample_shifted_WT(T, Z, L);
    // W_T = sqrt(T) * (Z + theta)
    EXPECT_NEAR(W_T, std::sqrt(T) * (Z + 0.5), 1e-12);
    // L = exp(theta*Z + 0.5*theta^2)
    Real expected_L = std::exp(0.5 * Z + 0.5 * 0.25);
    EXPECT_NEAR(L, expected_L, 1e-12);
    // 验证: dP/dQ = 1/L
    // (Z + theta) 是 Q 下的样本, 由密度比推导得 dP/dQ = exp(-theta*Z - 0.5*theta^2)
    Real expected_weight = std::exp(-0.5 * Z - 0.5 * 0.25);
    EXPECT_NEAR(1.0 / L, expected_weight, 1e-12);
}

// 4. IS 价格匹配 BSM (auto_optimize)
TEST(ImportanceSampling, PriceMatchesBSM) {
    Real S0 = 100, K = 100, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Size n_paths = 20000;
    uint64_t seed = 42;
    ISConfig cfg;
    cfg.auto_optimize = true;
    ImportanceSampling is(cfg);
    auto result = is.price_european_call(S0, K, T, r, q, sigma, n_paths, seed);
    Real bsm = bsm_call_price(S0, K, T, r, q, sigma);
    // 5σ 容差
    Real tol = std::max(0.05, 5.0 * result.std_error);
    EXPECT_NEAR(result.price, bsm, tol);
}

// 5. IS 方差缩减 (vs 标准 MC): OTM Call (K > S0)
TEST(ImportanceSampling, VarianceReductionOTM) {
    Real S0 = 100, K = 130, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Size n_paths = 10000;
    uint64_t seed = 42;
    ISConfig cfg;
    cfg.auto_optimize = true;
    ImportanceSampling is(cfg);
    auto result = is.price_european_call(S0, K, T, r, q, sigma, n_paths, seed);
    // 对 OTM Call, theta* > 0 显著, IS 应实现方差缩减
    EXPECT_GT(result.variance_reduction_ratio, 1.0);
    // 验证价格匹配 BSM
    Real bsm = bsm_call_price(S0, K, T, r, q, sigma);
    Real tol = std::max(0.05, 5.0 * result.std_error);
    EXPECT_NEAR(result.price, bsm, tol);
}

// 6. 稀有事件 (deep OTM) 方差缩减显著
TEST(ImportanceSampling, DeepOTMVarianceReduction) {
    Real S0 = 100, K = 200, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Size n_paths = 20000;
    uint64_t seed = 7;
    ISConfig cfg;
    cfg.auto_optimize = true;
    ImportanceSampling is(cfg);
    auto result = is.price_european_call(S0, K, T, r, q, sigma, n_paths, seed);
    // Deep OTM (K=2*S0): IS 应有显著方差缩减 (>> 1)
    EXPECT_GT(result.variance_reduction_ratio, 2.0);
    // 验证 IS 价格与 BSM 一致 (5σ 容差)
    Real bsm = bsm_call_price(S0, K, T, r, q, sigma);
    Real tol = std::max(0.02, 5.0 * result.std_error);
    EXPECT_NEAR(result.price, bsm, tol);
    // 标准 MC 在 deep OTM 下应有大标准误差
    auto mc = standard_mc_european_call(S0, K, T, r, q, sigma, n_paths, seed);
    // IS 标准误差应远小于标准 MC
    EXPECT_LT(result.std_error, mc.std_error);
}

// ============================================================
// ============ Conditional MC 测试 (5 用例) ============
// ============================================================

// 1. 障碍 CMC 匹配标准 MC (CMC 使用 Bachelier-Lévy 线性障碍公式, 精确连续监控)
TEST(ConditionalMC, BarrierUpOutCallMatchesMC) {
    Real S0 = 100, K = 100, B = 130, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Size n_paths_cmc = 50000;
    Size n_paths_mc = 30000;
    Size n_steps_mc = 1000;  // 高分辨率以减小离散监控偏差
    uint64_t seed = 42;
    Real df = std::exp(-r * T);
    Real sqrtT = std::sqrt(T);
    // CMC 估计 (精确连续监控)
    Real sum = 0.0, sum_sq = 0.0;
    for (Size i = 0; i < n_paths_cmc; ++i) {
        Philox4x64 rng(seed, static_cast<uint64_t>(i) + 1);
        Real Z = generate_normal(rng);
        Real W_T = sqrtT * Z;
        Real surv;
        Real payoff = ConditionalMC::barrier_up_out_call_cmc(S0, K, B, T, r, q, sigma, W_T, surv);
        sum += payoff;
        sum_sq += payoff * payoff;
    }
    Real inv_n = 1.0 / static_cast<Real>(n_paths_cmc);
    Real cmc_price = df * (sum * inv_n);
    Real cmc_var = sum_sq * inv_n - (sum * inv_n) * (sum * inv_n);
    Real cmc_se = df * std::sqrt(cmc_var * inv_n);
    // 标准 MC (高分辨率, 仍有少量离散监控偏差)
    auto mc = standard_mc_barrier_uoc(S0, K, B, T, r, q, sigma, n_steps_mc, n_paths_mc, seed + 100);
    // 容差: 5σ (CMC) + 5σ (MC) + 0.2 (残余离散监控偏差)
    // CMC 给出精确连续监控价格; 标准 MC 离散监控会有系统性高估 (Broadie-Glasserman-Kou 1997)
    Real tol = std::max(0.3, 5.0 * cmc_se + 5.0 * mc.std_error + 0.2);
    EXPECT_NEAR(cmc_price, mc.price, tol);
    // CMC 应有显著方差缩减 (CMC 单步 vs MC 多步路径)
    EXPECT_LT(cmc_se, mc.std_error);
}

// 2. 亚式 CMC 方差缩减 vs 标准 MC
TEST(ConditionalMC, AsianArithmeticVarianceReduction) {
    Real S0 = 100, K = 100, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Size n_steps = 12;
    Size n_paths = 20000;
    uint64_t seed = 42;
    Real df = std::exp(-r * T);
    Real sqrtT = std::sqrt(T);
    // CMC 估计
    Real sum = 0.0, sum_sq = 0.0;
    for (Size i = 0; i < n_paths; ++i) {
        Philox4x64 rng(seed, static_cast<uint64_t>(i) + 1);
        Real Z = generate_normal(rng);
        Real W_T = sqrtT * Z;
        Real cond_mean;
        Real payoff = ConditionalMC::asian_arithmetic_call_cmc(S0, K, T, r, q, sigma,
                                                                  n_steps, W_T, cond_mean);
        sum += payoff;
        sum_sq += payoff * payoff;
    }
    Real inv_n = 1.0 / static_cast<Real>(n_paths);
    Real cmc_price = df * (sum * inv_n);
    Real cmc_var = sum_sq * inv_n - (sum * inv_n) * (sum * inv_n);
    Real cmc_se = df * std::sqrt(cmc_var * inv_n);
    // 标准 MC
    auto mc = standard_mc_asian_arith_call(S0, K, T, r, q, sigma, n_steps, n_paths, seed + 100);
    // CMC 估计与标准 MC 应在统计上一致 (5σ + 5σ 合成容差)
    Real tol = std::max(0.1, 5.0 * cmc_se + 5.0 * mc.std_error);
    EXPECT_NEAR(cmc_price, mc.price, tol);
    // CMC 应有方差缩减
    EXPECT_LT(cmc_se, mc.std_error);
}

// 3. survival_prob 范围 [0, 1]
TEST(ConditionalMC, SurvivalProbRange) {
    Real S0 = 100, K = 100, B = 130, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Real mu = r - q - 0.5 * sigma * sigma;
    Real b_0 = std::log(B / S0) / sigma;             // t=0 W-空间障碍
    Real b_T = (std::log(B / S0) - mu * T) / sigma;  // t=T W-空间障碍
    // 情形 1: W_T << b_T (远离障碍) → survival_prob 接近 1
    Real surv1;
    ConditionalMC::barrier_up_out_call_cmc(S0, K, B, T, r, q, sigma, b_T - 5.0, surv1);
    EXPECT_GE(surv1, 0.0);
    EXPECT_LE(surv1, 1.0);
    EXPECT_NEAR(surv1, 1.0, 1e-3);
    // 情形 2: W_T 接近 b_T (从下方) → survival_prob 在 (0,1) 之间
    Real surv2;
    ConditionalMC::barrier_up_out_call_cmc(S0, K, B, T, r, q, sigma, b_T - 0.05, surv2);
    EXPECT_GT(surv2, 0.0);
    EXPECT_LT(surv2, 1.0);
    // 情形 3: W_T >= b_T → survival_prob = 0 (终端已越过 t=T 障碍)
    Real surv3;
    Real payoff3 = ConditionalMC::barrier_up_out_call_cmc(S0, K, B, T, r, q, sigma, b_T + 0.5, surv3);
    EXPECT_NEAR(surv3, 0.0, 1e-15);
    EXPECT_NEAR(payoff3, 0.0, 1e-15);
}

// 4. 条件期望正确性: E[A_T | W_T] 与直接计算一致
TEST(ConditionalMC, ConditionalMeanCorrectness) {
    Real S0 = 100, K = 100, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Size n_steps = 12;
    Real W_T = 0.3;
    Real mu = r - q - 0.5 * sigma * sigma;
    Real dt = T / static_cast<Real>(n_steps);
    Real inv_n = 1.0 / static_cast<Real>(n_steps);
    // 直接计算 E[A_T | W_T=w]
    Real expected_mean = 0.0;
    for (Size i = 0; i < n_steps; ++i) {
        Real t_i = static_cast<Real>(i + 1) * dt;
        Real mean_W = (t_i / T) * W_T;
        Real var_W = t_i * (T - t_i) / T;
        expected_mean += S0 * std::exp(mu * t_i + sigma * mean_W + 0.5 * sigma * sigma * var_W);
    }
    expected_mean *= inv_n;
    Real cond_mean;
    ConditionalMC::asian_arithmetic_call_cmc(S0, K, T, r, q, sigma, n_steps, W_T, cond_mean);
    EXPECT_NEAR(cond_mean, expected_mean, 1e-10 * expected_mean);
}

// 5. 障碍 CMC 条件期望非负且与 S_T 关系合理
TEST(ConditionalMC, BarrierPayoffConsistency) {
    Real S0 = 100, K = 100, B = 150, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Real mu_T = (r - q - 0.5 * sigma * sigma) * T;
    Real b = std::log(B / S0) / sigma;
    // W_T = b - 1 (远离障碍, survival_prob 高)
    Real W_T = b - 1.0;
    Real S_T = S0 * std::exp(mu_T + sigma * W_T);
    Real surv;
    Real payoff = ConditionalMC::barrier_up_out_call_cmc(S0, K, B, T, r, q, sigma, W_T, surv);
    // payoff = max(S_T - K, 0) * survival_prob
    Real expected_payoff = std::max(S_T - K, 0.0) * surv;
    EXPECT_NEAR(payoff, expected_payoff, 1e-12);
    EXPECT_GE(payoff, 0.0);
    EXPECT_LE(surv, 1.0);
    EXPECT_GE(surv, 0.0);
}

// ============================================================
// ============ Stratified Sampling 测试 (6 用例) ============
// ============================================================

// 1. 分层 U 分布正确: 每层 [i/K, (i+1)/K) 内有 n/K 个样本
TEST(StratifiedSampling, StratifiedUniformDistribution) {
    StratifiedConfig cfg;
    cfg.n_strata = 10;
    cfg.neyman_allocation = false;
    StratifiedSampling ss(cfg);
    Size n_paths = 1000;
    Philox4x64 rng(42);
    auto uniforms = ss.generate_stratified_uniforms(n_paths, rng);
    EXPECT_EQ(uniforms.size(), n_paths);
    Size K = cfg.n_strata;
    // 每层应有 n_paths / K = 100 个样本
    std::vector<Size> count(K, 0);
    for (Real u : uniforms) {
        EXPECT_GE(u, 0.0);
        EXPECT_LT(u, 1.0);
        Size k = static_cast<Size>(u * K);
        if (k >= K) k = K - 1;  // 边界保护
        count[k]++;
    }
    for (Size k = 0; k < K; ++k) {
        EXPECT_EQ(count[k], n_paths / K);
    }
    // 每层均值应接近 (k + 0.5) / K
    std::vector<Real> sum(K, 0.0);
    for (Real u : uniforms) {
        Size k = static_cast<Size>(u * K);
        if (k >= K) k = K - 1;
        sum[k] += u;
    }
    for (Size k = 0; k < K; ++k) {
        Real layer_mean = sum[k] / static_cast<Real>(count[k]);
        Real expected = (static_cast<Real>(k) + 0.5) / static_cast<Real>(K);
        EXPECT_NEAR(layer_mean, expected, 0.02);
    }
}

// 2. Neyman 分配合理: 高方差层 (深度 ITM/OTM) 样本数应 >= 低方差层
TEST(StratifiedSampling, NeymanAllocationReasonable) {
    StratifiedConfig cfg;
    cfg.n_strata = 10;
    cfg.neyman_allocation = true;
    cfg.pilot_samples = 200;
    StratifiedSampling ss(cfg);
    Real S0 = 100, K = 100, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Size n_paths = 1000;
    auto n_k = ss.neyman_allocation(S0, K, T, r, q, sigma, n_paths, 42);
    // 总样本数应为 n_paths
    Size total = std::accumulate(n_k.begin(), n_k.end(), Size(0));
    EXPECT_EQ(total, n_paths);
    // 各层样本数 > 0
    for (Size k = 0; k < n_k.size(); ++k) {
        EXPECT_GT(n_k[k], 0u);
    }
    // 对 BSM Call, 中间层 (S_T ≈ K 附近) 方差最大, 应分配更多样本
    // 第 K/2 层 (中间) 应比第 0 层 (深度 ITM, payoff 几乎确定) 多
    Size mid = n_k.size() / 2;
    EXPECT_GE(n_k[mid], n_k[0]);
}

// 3. 分层 MC 匹配 BSM
TEST(StratifiedSampling, StratifiedMCMatchesBSM) {
    Real S0 = 100, K = 100, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Size n_paths = 10000;
    uint64_t seed = 42;
    StratifiedConfig cfg;
    cfg.n_strata = 20;
    cfg.neyman_allocation = false;  // 等比例分配
    StratifiedSampling ss(cfg);
    auto result = ss.price_european_call_stratified(S0, K, T, r, q, sigma, n_paths, seed);
    Real bsm = bsm_call_price(S0, K, T, r, q, sigma);
    // 5σ 容差
    Real tol = std::max(0.05, 5.0 * result.std_error);
    EXPECT_NEAR(result.price, bsm, tol);
}

// 4. Neyman 分配下方差缩减 (vs 标准 MC)
TEST(StratifiedSampling, VarianceReductionNeyman) {
    Real S0 = 100, K = 100, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Size n_paths = 10000;
    uint64_t seed = 42;
    StratifiedConfig cfg;
    cfg.n_strata = 20;
    cfg.neyman_allocation = true;
    cfg.pilot_samples = 100;
    StratifiedSampling ss(cfg);
    auto result = ss.price_european_call_stratified(S0, K, T, r, q, sigma, n_paths, seed);
    EXPECT_GT(result.variance_reduction_ratio, 1.0);
    // 标准 MC 应有更大标准误差
    auto mc = standard_mc_european_call(S0, K, T, r, q, sigma, n_paths, seed);
    EXPECT_LT(result.std_error, mc.std_error);
}

// 5. 等比例分配下方差缩减 (vs 标准 MC)
TEST(StratifiedSampling, VarianceReductionProportional) {
    Real S0 = 100, K = 110, T = 1.0, r = 0.05, q = 0.0, sigma = 0.25;
    Size n_paths = 10000;
    uint64_t seed = 99;
    StratifiedConfig cfg;
    cfg.n_strata = 16;
    cfg.neyman_allocation = false;
    StratifiedSampling ss(cfg);
    auto result = ss.price_european_call_stratified(S0, K, T, r, q, sigma, n_paths, seed);
    // 等比例分配应也有方差缩减 (分层永远不增方差)
    EXPECT_GT(result.variance_reduction_ratio, 1.0);
    Real bsm = bsm_call_price(S0, K, T, r, q, sigma);
    Real tol = std::max(0.05, 5.0 * result.std_error);
    EXPECT_NEAR(result.price, bsm, tol);
}

// 6. 与 Sobol QMC 对比: 分层 MC 应与 Sobol QMC 精度相当 (都在 BSM 附近)
TEST(StratifiedSampling, CompareWithSobolQMC) {
    Real S0 = 100, K = 100, T = 1.0, r = 0.05, q = 0.0, sigma = 0.2;
    Real df = std::exp(-r * T);
    Real bsm = bsm_call_price(S0, K, T, r, q, sigma);
    Size n_paths = 8192;
    // 分层 MC (等比例分配)
    StratifiedConfig cfg;
    cfg.n_strata = 32;
    cfg.neyman_allocation = false;
    StratifiedSampling ss(cfg);
    auto strat = ss.price_european_call_stratified(S0, K, T, r, q, sigma, n_paths, 42);
    // Sobol QMC (1 维)
    SobolSequence seq(1);
    Real sum_sobol = 0.0;
    Real mu_T = (r - q - 0.5 * sigma * sigma) * T;
    Real sqrtT = std::sqrt(T);
    for (Size i = 1; i <= n_paths; ++i) {
        auto u = seq(static_cast<uint64_t>(i));
        Real Z = inv_normal_cdf(u[0]);
        Real S_T = S0 * std::exp(mu_T + sigma * sqrtT * Z);
        sum_sobol += std::max(S_T - K, 0.0);
    }
    Real sobol_price = df * sum_sobol / static_cast<Real>(n_paths);
    // 两者都应在 BSM 附近
    Real strat_tol = std::max(0.1, 5.0 * strat.std_error);
    EXPECT_NEAR(strat.price, bsm, strat_tol);
    EXPECT_NEAR(sobol_price, bsm, std::max(0.1, bsm * 0.01));
    // 两者之间应一致 (互在 5σ 内)
    Real diff_tol = std::max(0.1, 5.0 * strat.std_error + bsm * 0.01);
    EXPECT_NEAR(strat.price, sobol_price, diff_tol);
}
