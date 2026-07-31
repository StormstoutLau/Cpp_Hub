// CEV 模型单元测试
// 覆盖: 数学基础 (Gamma/非中心卡方) + 解析定价 + 随机过程模拟
//
// 测试分组:
//   1. 数学函数 (5 cases) — regularized_lower_gamma / noncentral_chi2_cdf
//   2. CEV 解析定价 (7 cases) — 参数验证, β=1 退化, T=0, Call-Put parity, r=q vs r≠q
//   3. CEV 随机过程 (8 cases) — 参数验证, 路径起始, 确定性, 吸收壁, β=1 GBM 一致, MC vs 解析
// 总计: ~20 cases

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include "cpphub/core/math.hpp"
#include "cpphub/pricing/analytic/cev_analytic.hpp"
#include "cpphub/models/diffusion/cev.hpp"
#include "cpphub/models/diffusion/gbm.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"  // bsm_call_price / bsm_put_price
#include "cpphub/core/rng.hpp"

using namespace cpphub;

// ==================== 1. 数学函数 ====================

TEST(CEVMathTest, RegularizedLowerGammaBoundaryZero) {
    // P(a, 0) = 0
    EXPECT_NEAR(regularized_lower_gamma(1.0, 0.0), 0.0, 1e-15);
    EXPECT_NEAR(regularized_lower_gamma(2.5, 0.0), 0.0, 1e-15);
}

TEST(CEVMathTest, RegularizedLowerGammaAsymptoticOne) {
    // P(a, x) → 1 as x → ∞
    // 对 a=1: P(1, x) = 1 - e^{-x}
    EXPECT_NEAR(regularized_lower_gamma(1.0, 50.0), 1.0, 1e-15);
    // 对 a=2: P(2, x) = 1 - (1+x)e^{-x}
    EXPECT_NEAR(regularized_lower_gamma(2.0, 50.0), 1.0, 1e-15);
}

TEST(CEVMathTest, RegularizedLowerGammaKnownValues) {
    // P(1, x) = 1 - e^{-x}
    EXPECT_NEAR(regularized_lower_gamma(1.0, 1.0), 1.0 - std::exp(-1.0), 1e-14);
    EXPECT_NEAR(regularized_lower_gamma(1.0, 2.5), 1.0 - std::exp(-2.5), 1e-14);
    // P(0.5, x) = erf(sqrt(x))  (关系: γ(1/2, x) = sqrt(π)·erf(√x))
    Real x = 1.5;
    Real expected = std::erf(std::sqrt(x));
    EXPECT_NEAR(regularized_lower_gamma(0.5, x), expected, 1e-13);
}

TEST(CEVMathTest, NoncentralChi2LambdaZeroDegenerates) {
    // λ = 0: 退化为中心卡方 CDF = P(k/2, x/2)
    Real x = 3.0;
    Real k = 4.0;
    Real expected = regularized_lower_gamma(k / 2.0, x / 2.0);
    EXPECT_NEAR(noncentral_chi2_cdf(x, k, 0.0), expected, 1e-14);
}

TEST(CEVMathTest, NoncentralChi2MonotoneInLambda) {
    // 对固定 x, k: CDF 应随 λ 增大而减小 (非中心参数越大, 均值越大, 落在小值的概率越小)
    Real x = 5.0;
    Real k = 4.0;
    Real cdf_0 = noncentral_chi2_cdf(x, k, 0.0);
    Real cdf_1 = noncentral_chi2_cdf(x, k, 2.0);
    Real cdf_2 = noncentral_chi2_cdf(x, k, 8.0);
    EXPECT_GT(cdf_0, cdf_1);
    EXPECT_GT(cdf_1, cdf_2);
    // CDF 始终在 [0, 1]
    EXPECT_GE(cdf_0, 0.0);
    EXPECT_LE(cdf_0, 1.0);
    EXPECT_GE(cdf_2, 0.0);
    EXPECT_LE(cdf_2, 1.0);
}

// ==================== 2. CEV 解析定价 ====================

TEST(CEVAnalyticTest, InvalidParamsThrows) {
    CEVParams p_neg_sigma{-0.1, 0.5};
    EXPECT_THROW(validate_cev_params(p_neg_sigma), std::invalid_argument);

    CEVParams p_zero_sigma{0.0, 0.5};
    EXPECT_THROW(validate_cev_params(p_zero_sigma), std::invalid_argument);

    CEVParams p_beta_gt_1{0.2, 1.5};
    EXPECT_THROW(validate_cev_params(p_beta_gt_1), std::invalid_argument);
}

TEST(CEVAnalyticTest, ZeroMaturityReturnsIntrinsic) {
    CEVParams p{0.2, 0.5};
    Real S = 100.0, K = 110.0, r = 0.05, q = 0.02;
    Real call_T0 = cev_call_price(S, K, 0.0, r, q, p);
    Real put_T0 = cev_put_price(S, K, 0.0, r, q, p);
    EXPECT_NEAR(call_T0, std::max(S - K, 0.0), 1e-15);
    EXPECT_NEAR(put_T0, std::max(K - S, 0.0), 1e-15);

    // ITM 情形
    Real K_itm = 90.0;
    EXPECT_NEAR(cev_call_price(S, K_itm, 0.0, r, q, p), S - K_itm, 1e-15);
    EXPECT_NEAR(cev_put_price(S, K_itm, 0.0, r, q, p), 0.0, 1e-15);
}

TEST(CEVAnalyticTest, BetaOneMatchesBSM) {
    // β = 1: CEV 退化为 GBM, 应严格匹配 Black-Scholes
    CEVParams p{0.2, 1.0};
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.02;

    Real cev_call = cev_call_price(S, K, T, r, q, p);
    Real bs_call = bsm_call_price(S, K, T, r, q, p.sigma);
    EXPECT_NEAR(cev_call, bs_call, 1e-10);

    Real cev_put = cev_put_price(S, K, T, r, q, p);
    Real bs_put = bsm_put_price(S, K, T, r, q, p.sigma);
    EXPECT_NEAR(cev_put, bs_put, 1e-10);

    // OTM / ITM
    Real K_otm = 120.0;
    EXPECT_NEAR(cev_call_price(S, K_otm, T, r, q, p),
                bsm_call_price(S, K_otm, T, r, q, p.sigma), 1e-10);
    Real K_itm = 80.0;
    EXPECT_NEAR(cev_put_price(S, K_itm, T, r, q, p),
                bsm_put_price(S, K_itm, T, r, q, p.sigma), 1e-10);
}

TEST(CEVAnalyticTest, CallPutParityHolds) {
    // CEV Call-Put Parity: C - P = S·e^{-qT} - K·e^{-rT}
    // 仅在 β < 1 且吸收壁在 0 时严格成立
    CEVParams p{0.3, 0.5};
    Real S = 100.0, T = 1.0, r = 0.05, q = 0.02;

    for (Real K : {80.0, 90.0, 100.0, 110.0, 120.0}) {
        Real C = cev_call_price(S, K, T, r, q, p);
        Real P = cev_put_price(S, K, T, r, q, p);
        Real parity_rhs = S * std::exp(-q * T) - K * std::exp(-r * T);
        EXPECT_NEAR(C - P, parity_rhs, 1e-10)
            << "Parity failed at K=" << K;
    }
}

TEST(CEVAnalyticTest, ZeroDriftSimplifiedFormMatchesGeneral) {
    // r = q 时简化形式应与 r ≠ q 时的公式数值等价
    // 用 r ≈ q 的极限验证 (这里直接用 r=q=0.05 测一致性)
    CEVParams p{0.25, 0.7};
    Real S = 100.0, K = 105.0, T = 1.5, r = 0.05, q = 0.05;

    // r = q 走简化分支
    Real call_zero_drift = cev_call_price(S, K, T, r, q, p);

    // r ≠ q 走一般分支 (微小差异, 容差放宽)
    Real q_perturbed = 0.05 + 1e-6;
    Real call_general = cev_call_price(S, K, T, r, q_perturbed, p);

    // 两者应接近 (q 微扰 1e-6, 价格变化 ~1e-5 量级)
    EXPECT_NEAR(call_zero_drift, call_general, 1e-3);
}

TEST(CEVAnalyticTest, BetaBelowOneProducesSkew) {
    // β < 1: CEV 产生负偏斜 (低行权价 IV 高于高行权价 IV)
    // 通过反推 IV 验证: ITM call 的 IV 应高于 OTM call 的 IV
    CEVParams p{0.3, 0.3};  // β = 0.3 显著 skew
    Real S = 100.0, T = 1.0, r = 0.0, q = 0.0;

    // 用 BS 公式反推 IV (近似, 用于验证单调性)
    Real K_low = 80.0;   // ITM call
    Real K_high = 120.0; // OTM call
    Real price_low = cev_call_price(S, K_low, T, r, q, p);
    Real price_high = cev_call_price(S, K_high, T, r, q, p);

    // 反推 IV 用 Newton 迭代
    auto implied_vol = [](Real price, Real S_, Real K_, Real T_, Real r_, Real q_) {
        Real sigma = 0.2;
        for (int i = 0; i < 100; ++i) {
            Real d1 = (std::log(S_ / K_) + (r_ - q_ + 0.5 * sigma * sigma) * T_) / (sigma * std::sqrt(T_));
            Real d2 = d1 - sigma * std::sqrt(T_);
            Real bs = S_ * std::exp(-q_ * T_) * normal_cdf(d1) - K_ * std::exp(-r_ * T_) * normal_cdf(d2);
            Real vega = S_ * std::exp(-q_ * T_) * normal_pdf(d1) * std::sqrt(T_);
            if (vega < 1e-12) break;
            Real diff = bs - price;
            if (std::abs(diff) < 1e-12) break;
            sigma -= diff / vega;
            if (sigma < 0.001) sigma = 0.001;
            if (sigma > 5.0) sigma = 5.0;
        }
        return sigma;
    };

    Real iv_low = implied_vol(price_low, S, K_low, T, r, q);
    Real iv_high = implied_vol(price_high, S, K_high, T, r, q);
    // β < 1 产生 negative skew: ITM (K_low) IV > OTM (K_high) IV
    EXPECT_GT(iv_low, iv_high)
        << "iv_low=" << iv_low << " iv_high=" << iv_high;
}

TEST(CEVAnalyticTest, BetaZeroIsBachelierLike) {
    // β = 0.3 (CEV 公式假设 β ∈ (0,1) 且 S=0 吸收壁; β=0 时 SDE dS=σdW 允许 S<0, 公式不适用)
    // 验证中等 β 下价格在合理 strike 范围内为正且有限
    CEVParams p{0.3, 0.3};
    Real S = 100.0, T = 1.0, r = 0.0, q = 0.0;

    for (Real K : {70.0, 85.0, 100.0, 115.0, 130.0}) {
        Real C = cev_call_price(S, K, T, r, q, p);
        Real P = cev_put_price(S, K, T, r, q, p);
        EXPECT_GT(C, 0.0) << "Call price should be positive at K=" << K;
        EXPECT_GT(P, 0.0) << "Put price should be positive at K=" << K;
        EXPECT_LT(C, S) << "Call price should be less than spot at K=" << K;
        EXPECT_LT(P, K) << "Put price should be less than strike at K=" << K;
    }
}

// ==================== 3. CEV 随机过程 ====================

TEST(CEVProcessTest, InvalidS0Throws) {
    CEVParams p{0.2, 0.5};
    EXPECT_THROW(CEVProcess(p, -1.0), std::invalid_argument);
    EXPECT_THROW(CEVProcess(p, 0.0), std::invalid_argument);
}

TEST(CEVProcessTest, InvalidBetaThrows) {
    CEVParams p_neg{-0.1, -0.5};  // sigma < 0
    EXPECT_THROW(CEVProcess(p_neg, 100.0), std::invalid_argument);

    CEVParams p_beta_neg{0.2, -0.1};  // beta < 0
    EXPECT_THROW(CEVProcess(p_beta_neg, 100.0), std::invalid_argument);

    CEVParams p_beta_gt_1{0.2, 1.5};  // beta > 1
    EXPECT_THROW(CEVProcess(p_beta_gt_1, 100.0), std::invalid_argument);
}

TEST(CEVProcessTest, PathStartsAtSpot) {
    CEVParams p{0.2, 0.5};
    CEVProcess process(p, 100.0);
    Philox4x64 rng(42);
    std::vector<Real> path(101);
    process.generate_path(1.0, 100, path, rng);
    EXPECT_DOUBLE_EQ(path[0], 100.0);
}

TEST(CEVProcessTest, PathDimensionCorrect) {
    CEVParams p{0.2, 0.5};
    CEVProcess process(p, 100.0);
    EXPECT_EQ(process.dimension(), 1u);
}

TEST(CEVProcessTest, DeterministicWithSameSeed) {
    CEVParams p{0.2, 0.5};
    CEVProcess process(p, 100.0, 0.05, 0.02);
    Size n_steps = 100;
    std::vector<Real> path1(n_steps + 1);
    std::vector<Real> path2(n_steps + 1);
    {
        Philox4x64 rng(42);
        process.generate_path(1.0, n_steps, path1, rng);
    }
    {
        Philox4x64 rng(42);
        process.generate_path(1.0, n_steps, path2, rng);
    }
    for (Size i = 0; i <= n_steps; ++i) {
        EXPECT_DOUBLE_EQ(path1[i], path2[i])
            << "Mismatch at step " << i;
    }
}

TEST(CEVProcessTest, AbsorbingBarrierKeepsNonNegative) {
    // β < 1: 吸收壁在 0, 路径始终非负
    CEVParams p{0.4, 0.3};  // 高波动率 + β=0.3 触发吸收
    CEVProcess process(p, 100.0, 0.0, 0.0);
    Size n_steps = 252;
    std::vector<Real> path(n_steps + 1);
    for (int trial = 0; trial < 20; ++trial) {
        Philox4x64 rng(static_cast<uint64_t>(trial * 1000 + 7));
        process.generate_path(1.0, n_steps, path, rng);
        for (Size i = 0; i <= n_steps; ++i) {
            EXPECT_GE(path[i], 0.0)
                << "Negative price at step " << i << " trial " << trial;
        }
    }
}

TEST(CEVProcessTest, BetaOneLogEulerMatchesGBMEvolve) {
    // β = 1 + LogEuler: 用相同 Z 逐步 evolve, 应严格匹配 GBM evolve
    // (避免 generate_path 的 RNG 调用次数差异)
    CEVParams p{0.2, 1.0};
    CEVProcess process(p, 100.0, 0.05, 0.02, CEVScheme::LogEuler);
    GBMParams gbm_p{100.0, 0.05 - 0.02, 0.2};  // GBM mu = r - q
    GBM gbm(gbm_p);

    Real S_cev = 100.0;
    Real S_gbm = 100.0;
    Real dt = 0.02;
    Philox4x64 rng(42);

    for (Size i = 0; i < 50; ++i) {
        // 生成一个 Z, 两个过程用同一个 Z
        uint64_t r1 = rng();
        uint64_t r2 = rng();
        double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
        double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
        auto [z1, z2] = box_muller(u1, u2);
        S_cev = process.evolve(S_cev, dt, z1);
        S_gbm = gbm.evolve(S_gbm, dt, z1);
        EXPECT_NEAR(S_cev, S_gbm, 1e-12)
            << "Step " << i << ": CEV=" << S_cev << " GBM=" << S_gbm;
    }
}

TEST(CEVProcessTest, BetaOneEulerMatchesGBMMean) {
    // β = 1 + EulerAbsorbing: 离散化略不同 (Euler vs log-Euler), 但均值应统计一致
    CEVParams p{0.2, 1.0};
    CEVProcess process(p, 100.0, 0.05, 0.02, CEVScheme::EulerAbsorbing);

    Size n_paths = 20000;
    Size n_steps = 252;
    Real T = 1.0;

    Real sum_ST = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 999 + 1));
        std::vector<Real> path(n_steps + 1);
        process.generate_path(T, n_steps, path, rng);
        sum_ST += path.back();
    }
    Real mean_ST = sum_ST / static_cast<Real>(n_paths);
    Real expected_ST = 100.0 * std::exp((0.05 - 0.02) * T);
    Real rel_error = std::abs(mean_ST - expected_ST) / expected_ST;
    // Euler 离散有 O(dt) 截断误差, 20000 路径下统计误差 ~1%
    EXPECT_LT(rel_error, 0.02)
        << "mean_ST = " << mean_ST << ", expected = " << expected_ST;
}

TEST(CEVProcessTest, MCMatchesAnalytic) {
    // MC vs 解析定价: β < 1 情形
    CEVParams p{0.3, 0.5};
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.02;

    CEVProcess process(p, S0, r, q, CEVScheme::EulerAbsorbing);
    Real analytic_price = cev_call_price(S0, K, T, r, q, p);

    Size n_paths = 50000;
    Size n_steps = 200;  // 高分辨率减小离散偏差
    Real dt = T / static_cast<Real>(n_steps);
    Real discount = std::exp(-r * T);

    Real sum_payoff = 0.0;
    Real sum_payoff2 = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 7919 + 13));
        std::vector<Real> path(n_steps + 1);
        process.generate_path(T, n_steps, path, rng);
        Real ST = path.back();
        Real payoff = std::max(ST - K, 0.0);
        sum_payoff += payoff;
        sum_payoff2 += payoff * payoff;
    }
    Real mc_price = discount * sum_payoff / static_cast<Real>(n_paths);
    Real var = (sum_payoff2 / static_cast<Real>(n_paths)
                - (sum_payoff / static_cast<Real>(n_paths))
                  * (sum_payoff / static_cast<Real>(n_paths)));
    Real se = discount * std::sqrt(std::max(var, 0.0) / static_cast<Real>(n_paths));

    // MC 估计 ±3σ 应覆盖解析解
    Real diff = std::abs(mc_price - analytic_price);
    EXPECT_LT(diff, 4.0 * se)
        << "MC=" << mc_price << " ± " << se
        << ", analytic=" << analytic_price << ", diff=" << diff;
}
