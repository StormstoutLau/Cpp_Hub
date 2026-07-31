// Bates 模型单元测试
// 覆盖: 特征函数 (CF) + 随机过程模拟 + COS 半解析定价
//
// 测试分组:
//   1. Bates CF (10 cases) — 单位模, λ=0 退化, 跳跃 CF 性质, 衰减, 平滑性
//   2. Bates 过程 (10 cases) — 参数验证, 路径性质, λ=0 匹配 Heston, Poisson/LogNormal 采样
//   3. Bates 定价 (5 cases) — COS 定价, λ=0 匹配 Heston, Call-Put parity, 跳跃效应, MC vs COS
// 总计: 25 cases

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include "cpphub/pricing/analytic/bates_cf.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include "cpphub/models/diffusion/bates.hpp"
#include "cpphub/models/diffusion/heston.hpp"
#include "cpphub/pricing/fourier/cos_method.hpp"
#include "cpphub/core/rng.hpp"

using namespace cpphub;

// ============ 1. Bates 特征函数 ============

TEST(BatesCFTest, AtZeroReturnsOne) {
    BatesCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.2, 0.03, 0.0};
    Complex phi = bates_characteristic_function(Complex(0, 0), 1.0, 100.0, p);
    EXPECT_NEAR(std::real(phi), 1.0, 1e-12);
    EXPECT_NEAR(std::imag(phi), 0.0, 1e-12);
}

TEST(BatesCFTest, UnitModulusForRealU) {
    BatesCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.2, 0.03, 0.0};
    for (Real u_val = 0.1; u_val <= 5.0; u_val += 0.1) {
        Complex u(u_val, 0);
        Complex phi = bates_characteristic_function(u, 1.0, 100.0, p);
        EXPECT_LE(std::abs(phi), 1.0 + 1e-10)
            << "Failed at u=" << u_val << ", |phi|=" << std::abs(phi);
    }
}

TEST(BatesCFTest, LambdaZeroMatchesHeston) {
    // λ=0: Bates CF 应严格等于 Heston CF (相同 Heston 参数, 无 drift 调整)
    BatesCFParams bp{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, -0.1, 0.2, 0.03, 0.0};
    HestonCFParams hp{0.04, 1.5, 0.04, 0.3, -0.5, 0.03, 0.0};
    Real S0 = 100.0, tau = 1.0;
    for (Real u_val : {0.1, 0.5, 1.0, 2.0, 5.0}) {
        Complex u(u_val, 0);
        Complex bates = bates_characteristic_function(u, tau, S0, bp);
        Complex heston = heston_characteristic_function(u, tau, S0, hp);
        EXPECT_NEAR(std::real(bates), std::real(heston), 1e-12)
            << "Real mismatch at u=" << u_val;
        EXPECT_NEAR(std::imag(bates), std::imag(heston), 1e-12)
            << "Imag mismatch at u=" << u_val;
    }
}

TEST(BatesCFTest, MertonJumpCFAtZeroIsOne) {
    // φ_J(0, τ) = exp(λτ(e^0 - 1)) = exp(0) = 1
    Complex phi_J = merton_jump_cf(Complex(0, 0), 1.0, 0.5, -0.1, 0.2);
    EXPECT_NEAR(std::real(phi_J), 1.0, 1e-14);
    EXPECT_NEAR(std::imag(phi_J), 0.0, 1e-14);
}

TEST(BatesCFTest, MertonJumpCFReducesModulus) {
    // λ>0 时 |φ_J(u)| < 1 for u != 0 (跳跃增加方差, CF 衰减更快)
    Real tau = 1.0, lambda = 0.5, mu_J = -0.1, sigma_J = 0.2;
    for (Real u_val : {0.5, 1.0, 2.0, 5.0}) {
        Complex u(u_val, 0);
        Complex phi_J = merton_jump_cf(u, tau, lambda, mu_J, sigma_J);
        EXPECT_LT(std::abs(phi_J), 1.0 - 1e-10)
            << "Failed at u=" << u_val << ", |phi_J|=" << std::abs(phi_J);
    }
}

TEST(BatesCFTest, JumpCompensationFormula) {
    // m = E[J-1] = exp(μ_J + σ_J²/2) - 1
    EXPECT_NEAR(bates_jump_compensation(0.0, 0.0), 0.0, 1e-15);
    EXPECT_NEAR(bates_jump_compensation(0.0, 0.2), std::exp(0.02) - 1.0, 1e-15);
    EXPECT_NEAR(bates_jump_compensation(-0.1, 0.3), std::exp(-0.1 + 0.045) - 1.0, 1e-15);
}

TEST(BatesCFTest, HigherLambdaDecaysFaster) {
    // 更大的 λ → CF 衰减更快 (更重的尾部)
    Real tau = 1.0, S0 = 100.0;
    BatesCFParams p_low{0.04, 1.5, 0.04, 0.3, -0.5, 0.1, -0.1, 0.2, 0.03, 0.0};
    BatesCFParams p_high{0.04, 1.5, 0.04, 0.3, -0.5, 1.0, -0.1, 0.2, 0.03, 0.0};
    Complex u(2.0, 0);
    Real abs_low = std::abs(bates_characteristic_function(u, tau, S0, p_low));
    Real abs_high = std::abs(bates_characteristic_function(u, tau, S0, p_high));
    EXPECT_LT(abs_high, abs_low)
        << "Higher lambda should decay faster: low=" << abs_low << " high=" << abs_high;
}

TEST(BatesCFTest, CFDecaysAtInfinity) {
    BatesCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.2, 0.03, 0.0};
    Complex u(50.0, 0);
    Complex phi = bates_characteristic_function(u, 1.0, 100.0, p);
    EXPECT_LT(std::abs(phi), 0.01);
}

TEST(BatesCFTest, NoBranchCutDiscontinuity) {
    BatesCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.2, 0.03, 0.0};
    Real tau = 1.0;
    const Real eps = 1e-8;
    for (Real u_val = 0.1; u_val <= 10.0; u_val += 0.1) {
        Complex phi_fwd = bates_characteristic_function(Complex(u_val + eps, 0), tau, 100.0, p);
        Complex phi_bwd = bates_characteristic_function(Complex(u_val - eps, 0), tau, 100.0, p);
        Complex phi_mid = bates_characteristic_function(Complex(u_val, 0), tau, 100.0, p);
        Real smoothness = std::abs(phi_fwd + phi_bwd - Real(2) * phi_mid);
        EXPECT_LT(smoothness, 1e-6)
            << "Discontinuity at u=" << u_val << ", smoothness=" << smoothness;
    }
}

TEST(BatesCFTest, SymmetricJumpsZeroMuJ) {
    // μ_J = 0: 跳跃对称, 不引入偏度 (仅增加峰度)
    // 检查: 对实数 u, φ_J 的虚部来自 iμ_J u 项, μ_J=0 时该项消失
    // 但 σ_J²u²/2 仍为实数, 故 φ_J 在 μ_J=0 时对实数 u 为实数且 < 1
    Real tau = 1.0, lambda = 0.5, sigma_J = 0.2;
    for (Real u_val : {0.5, 1.0, 2.0}) {
        Complex u(u_val, 0);
        Complex phi_J = merton_jump_cf(u, tau, lambda, 0.0, sigma_J);
        // μ_J=0: exponent = -σ_J²u²/2 (纯实数), exp(λτ(e^{-σ_J²u²/2} - 1)) 为实数
        EXPECT_NEAR(std::imag(phi_J), 0.0, 1e-14)
            << "Imaginary part should be ~0 for symmetric jumps at u=" << u_val;
        EXPECT_LT(std::real(phi_J), 1.0)
            << "Real part < 1 at u=" << u_val;
    }
}

// ============ 2. Bates 随机过程 ============

TEST(BatesProcessTest, InvalidParamsThrows) {
    BatesParams p_neg_sigma{100, 0.04, 1.5, 0.04, -0.3, -0.5, 0.5, -0.1, 0.2, 0.03, 0.0};
    EXPECT_THROW(BatesProcess{p_neg_sigma}, std::invalid_argument);

    BatesParams p_neg_lambda{100, 0.04, 1.5, 0.04, 0.3, -0.5, -0.1, -0.1, 0.2, 0.03, 0.0};
    EXPECT_THROW(BatesProcess{p_neg_lambda}, std::invalid_argument);

    BatesParams p_zero_sigmaJ{100, 0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.0, 0.03, 0.0};
    EXPECT_THROW(BatesProcess{p_zero_sigmaJ}, std::invalid_argument);

    BatesParams p_bad_rho{100, 0.04, 1.5, 0.04, 0.3, 1.5, 0.5, -0.1, 0.2, 0.03, 0.0};
    EXPECT_THROW(BatesProcess{p_bad_rho}, std::invalid_argument);
}

TEST(BatesProcessTest, PathStartsAtSpot) {
    BatesParams p{100, 0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.2, 0.03, 0.0};
    BatesProcess process(p);
    Philox4x64 rng(42);
    std::vector<Real> path(101);
    process.generate_path(1.0, 100, path, rng);
    EXPECT_DOUBLE_EQ(path[0], 100.0);
}

TEST(BatesProcessTest, PathDimensionCorrect) {
    BatesParams p{100, 0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.2, 0.03, 0.0};
    BatesProcess process(p);
    EXPECT_EQ(process.dimension(), 2u);
}

TEST(BatesProcessTest, DeterministicWithSameSeed) {
    BatesParams p{100, 0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.2, 0.03, 0.0};
    BatesProcess process(p);
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

TEST(BatesProcessTest, LambdaZeroMatchesHestonPath) {
    // λ=0: Bates 路径应与 Heston 路径完全一致 (相同 seed, 相同 RNG 调用序列)
    BatesParams bp{100, 0.04, 1.5, 0.04, 0.3, -0.5, 0.0, -0.1, 0.2, 0.03, 0.0};
    HestonParams hp{100, 0.04, 1.5, 0.04, 0.3, -0.5, 0.03, 0.0};

    BatesProcess bates(bp, HestonScheme::FullTruncation);
    Heston heston(hp, HestonScheme::FullTruncation);

    Size n_steps = 100;
    std::vector<Real> path_bates(n_steps + 1);
    std::vector<Real> path_heston(n_steps + 1);
    {
        Philox4x64 rng(42);
        bates.generate_path(1.0, n_steps, path_bates, rng);
    }
    {
        Philox4x64 rng(42);
        heston.generate_path(1.0, n_steps, path_heston, rng);
    }
    for (Size i = 0; i <= n_steps; ++i) {
        EXPECT_NEAR(path_bates[i], path_heston[i], 1e-12)
            << "Path mismatch at step " << i
            << ": Bates=" << path_bates[i] << " Heston=" << path_heston[i];
    }
}

TEST(BatesProcessTest, PositiveSpotAlways) {
    // S 通过 exp 更新, 始终为正
    BatesParams p{100, 0.04, 1.5, 0.04, 0.3, -0.5, 2.0, -0.2, 0.3, 0.03, 0.0};
    BatesProcess process(p);
    Size n_steps = 252;
    std::vector<Real> path(n_steps + 1);
    for (int trial = 0; trial < 20; ++trial) {
        Philox4x64 rng(static_cast<uint64_t>(trial * 1000 + 7));
        process.generate_path(1.0, n_steps, path, rng);
        for (Size i = 0; i <= n_steps; ++i) {
            EXPECT_GT(path[i], 0.0)
                << "Non-positive spot at step " << i << " trial " << trial;
        }
    }
}

TEST(BatesProcessTest, JumpIncreasesPathVariance) {
    // 有跳跃的路径终端值方差应 > 无跳跃 (λ=0) 的方差
    BatesParams p_no_jump{100, 0.04, 1.5, 0.04, 0.3, -0.5, 0.0, -0.1, 0.2, 0.03, 0.0};
    BatesParams p_with_jump{100, 0.04, 1.5, 0.04, 0.3, -0.5, 3.0, -0.1, 0.2, 0.03, 0.0};

    BatesProcess proc_no(p_no_jump);
    BatesProcess proc_yes(p_with_jump);

    Size n_paths = 5000;
    Size n_steps = 100;
    Real T = 1.0;

    auto compute_var = [&](BatesProcess& proc) -> Real {
        Real sum = 0.0, sum2 = 0.0;
        for (Size j = 0; j < n_paths; ++j) {
            Philox4x64 rng(static_cast<uint64_t>(j * 7919 + 13));
            std::vector<Real> path(n_steps + 1);
            proc.generate_path(T, n_steps, path, rng);
            Real ST = path.back();
            sum += ST;
            sum2 += ST * ST;
        }
        Real mean = sum / n_paths;
        return sum2 / n_paths - mean * mean;
    };

    Real var_no = compute_var(proc_no);
    Real var_yes = compute_var(proc_yes);
    EXPECT_GT(var_yes, var_no * 1.1)
        << "Jumps should increase variance: no_jump=" << var_no
        << " with_jump=" << var_yes;
}

TEST(BatesProcessTest, PoissonSamplerMean) {
    // Poisson(λdt) 的样本均值应 ≈ λdt
    Real lambda_dt = 0.3;  // 期望 0.3 个跳跃/步
    Real sum = 0.0;
    int N = 100000;
    Philox4x64 rng(42);
    for (int i = 0; i < N; ++i) {
        sum += BatesProcess::sample_poisson(lambda_dt, rng);
    }
    Real mean = sum / N;
    // 大数定律: SE ≈ sqrt(λdt/N) ≈ 0.0017
    EXPECT_NEAR(mean, lambda_dt, 0.01)
        << "Poisson mean: expected=" << lambda_dt << " got=" << mean;
}

TEST(BatesProcessTest, LognormalJumpMean) {
    // E[J] = exp(μ_J + σ_J²/2)
    BatesParams p{100, 0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.2, 0.03, 0.0};
    BatesProcess process(p);
    Real expected_mean = std::exp(-0.1 + 0.5 * 0.04);
    Real sum = 0.0;
    int N = 100000;
    Philox4x64 rng(42);
    for (int i = 0; i < N; ++i) {
        sum += process.sample_lognormal_jump(rng);
    }
    Real mean = sum / N;
    // SE ≈ std::exp(μ_J + σ_J²/2) * std::sqrt((exp(σ_J²)-1)/N)  ≈ 0.0027
    EXPECT_NEAR(mean, expected_mean, 0.01)
        << "Lognormal jump mean: expected=" << expected_mean << " got=" << mean;
}

// ============ 3. Bates COS 定价 ============

TEST(BatesPricingTest, LambdaZeroMatchesHestonPrice) {
    // λ=0: Bates COS 定价应与 Heston COS 定价一致
    Real S0 = 100.0, r = 0.03, q = 0.0, T = 1.0, K = 100.0;

    BatesCFParams bp{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, -0.1, 0.2, r, q};
    HestonCFParams hp{0.04, 1.5, 0.04, 0.3, -0.5, r, q};

    auto bates_phi = make_bates_cf(S0, r, q, bp, T);
    auto heston_phi = make_heston_cf(S0, r, q, hp, T);

    COSEngine bates_engine(bates_phi, S0, r, q, T);
    COSEngine heston_engine(heston_phi, S0, r, q, T);

    Real bates_call = bates_engine.price_call(K);
    Real heston_call = heston_engine.price_call(K);

    EXPECT_NEAR(bates_call, heston_call, 1e-8)
        << "Bates=" << bates_call << " Heston=" << heston_call;
}

TEST(BatesPricingTest, CallPricePositive) {
    Real S0 = 100.0, r = 0.03, q = 0.0, T = 1.0;
    BatesCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.2, r, q};
    auto phi = make_bates_cf(S0, r, q, p, T);
    COSEngine engine(phi, S0, r, q, T);
    for (Real K : {80.0, 90.0, 100.0, 110.0, 120.0}) {
        Real C = engine.price_call(K);
        EXPECT_GT(C, 0.0) << "Call price should be positive at K=" << K;
        EXPECT_LT(C, S0) << "Call price should be less than spot at K=" << K;
    }
}

TEST(BatesPricingTest, CallPutParityHolds) {
    // Bates Call-Put Parity: C - P = S e^{-qT} - K e^{-rT}
    // (跳跃补偿已含在 drift 中, parity 严格成立)
    Real S0 = 100.0, r = 0.03, q = 0.0, T = 1.0;
    BatesCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.2, r, q};
    auto phi = make_bates_cf(S0, r, q, p, T);
    COSEngine engine(phi, S0, r, q, T);
    for (Real K : {80.0, 90.0, 100.0, 110.0, 120.0}) {
        Real C = engine.price_call(K);
        Real P = engine.price_put(K);
        Real parity_rhs = S0 * std::exp(-q * T) - K * std::exp(-r * T);
        EXPECT_NEAR(C - P, parity_rhs, 1e-6)
            << "Parity failed at K=" << K << ": C-P=" << (C - P)
            << " rhs=" << parity_rhs;
    }
}

TEST(BatesPricingTest, JumpIncreasesOTMCallPrice) {
    // 跳跃增加尾部, OTM call 应更贵 (λ>0 vs λ=0)
    Real S0 = 100.0, r = 0.03, q = 0.0, T = 1.0, K = 120.0;  // OTM call
    BatesCFParams p_no{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, -0.1, 0.2, r, q};
    BatesCFParams p_yes{0.04, 1.5, 0.04, 0.3, -0.5, 2.0, -0.1, 0.2, r, q};

    auto phi_no = make_bates_cf(S0, r, q, p_no, T);
    auto phi_yes = make_bates_cf(S0, r, q, p_yes, T);
    COSEngine engine_no(phi_no, S0, r, q, T);
    COSEngine engine_yes(phi_yes, S0, r, q, T);

    Real price_no = engine_no.price_call(K);
    Real price_yes = engine_yes.price_call(K);
    EXPECT_GT(price_yes, price_no)
        << "Jumps should increase OTM call: no_jump=" << price_no
        << " with_jump=" << price_yes;
}

TEST(BatesPricingTest, MCMatchesCOS) {
    // MC vs COS 半解析定价
    Real S0 = 100.0, r = 0.03, q = 0.0, T = 1.0, K = 100.0;
    BatesParams p{S0, 0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.2, r, q};
    BatesProcess process(p);

    // COS 半解析
    BatesCFParams cf_p{0.04, 1.5, 0.04, 0.3, -0.5, 0.5, -0.1, 0.2, r, q};
    auto phi = make_bates_cf(S0, r, q, cf_p, T);
    COSEngine engine(phi, S0, r, q, T);
    Real cos_price = engine.price_call(K);

    // MC
    Size n_paths = 30000;
    Size n_steps = 200;
    Real discount = std::exp(-r * T);

    Real sum_payoff = 0.0, sum_payoff2 = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 7919 + 13));
        std::vector<Real> path(n_steps + 1);
        process.generate_path(T, n_steps, path, rng);
        Real ST = path.back();
        Real payoff = std::max(ST - K, 0.0);
        sum_payoff += payoff;
        sum_payoff2 += payoff * payoff;
    }
    Real mc_price = discount * sum_payoff / n_paths;
    Real var = sum_payoff2 / n_paths - (sum_payoff / n_paths) * (sum_payoff / n_paths);
    Real se = discount * std::sqrt(std::max(var, 0.0) / n_paths);

    // MC ±4σ 应覆盖 COS 解 (Full Truncation 离散偏差 + MC 误差)
    Real diff = std::abs(mc_price - cos_price);
    EXPECT_LT(diff, 4.0 * se + 0.15)  // +0.15 容差给离散偏差
        << "MC=" << mc_price << " ± " << se
        << ", COS=" << cos_price << ", diff=" << diff;
}
