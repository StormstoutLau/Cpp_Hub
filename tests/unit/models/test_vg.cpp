// Variance Gamma (VG) 模型单元测试
// 覆盖: 特征函数 + 累积量 + 随机过程模拟 + COS 半解析定价
//
// 测试分组:
//   1. VG CF (8 cases) — 单位模, 衰减, omega, Feller, CF 一致性, 累积量
//   2. VG 过程 (7 cases) — 参数验证, 路径性质, Gamma 采样, MC 矩匹配
//   3. VG 定价 (5 cases) — COS 定价, Call-Put parity, ν→0 退化 BS, MC vs COS
// 总计: 20 cases

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include "cpphub/pricing/analytic/vg_analytic.hpp"
#include "cpphub/pricing/fourier/characteristic_functions.hpp"
#include "cpphub/pricing/fourier/cos_method.hpp"
#include "cpphub/models/diffusion/variance_gamma.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"  // bsm_call_price
#include "cpphub/calibration/calibrator.hpp"        // bsm_implied_vol
#include "cpphub/core/rng.hpp"

using namespace cpphub;

// ============ 1. VG 特征函数 ============

TEST(VGCFTest, AtZeroReturnsOne) {
    VGParams p{0.2, 0.5, -0.1};
    Complex phi = vg_characteristic_function(Complex(0, 0), 1.0, 100.0, 0.03, 0.0, p);
    EXPECT_NEAR(std::real(phi), 1.0, 1e-12);
    EXPECT_NEAR(std::imag(phi), 0.0, 1e-12);
}

TEST(VGCFTest, UnitModulusForRealU) {
    VGParams p{0.2, 0.5, -0.1};
    for (Real u_val = 0.1; u_val <= 5.0; u_val += 0.1) {
        Complex u(u_val, 0);
        Complex phi = vg_characteristic_function(u, 1.0, 100.0, 0.03, 0.0, p);
        EXPECT_LE(std::abs(phi), 1.0 + 1e-10)
            << "Failed at u=" << u_val << ", |phi|=" << std::abs(phi);
    }
}

TEST(VGCFTest, CFDecaysAtInfinity) {
    VGParams p{0.2, 0.5, -0.1};
    Complex u(50.0, 0);
    Complex phi = vg_characteristic_function(u, 1.0, 100.0, 0.03, 0.0, p);
    EXPECT_LT(std::abs(phi), 0.01);
}

TEST(VGCFTest, OmegaFormula) {
    // ω = (1/ν) ln(1 - θν - σ²ν/2)
    Real sigma = 0.2, nu = 0.5, theta = -0.1;
    Real expected = std::log(1.0 - theta * nu - 0.5 * sigma * sigma * nu) / nu;
    EXPECT_NEAR(vg_omega(sigma, nu, theta), expected, 1e-15);
}

TEST(VGCFTest, FellerConditionThrows) {
    // 1 - θν - σ²ν/2 > 0 必须满足
    VGParams p_bad{0.5, 2.0, 0.5};  // 1 - 0.5*2 - 0.5*0.25*2 = 1 - 1 - 0.25 = -0.25 < 0
    EXPECT_THROW(validate_vg_params(p_bad), std::invalid_argument);
    EXPECT_THROW(vg_omega(0.5, 2.0, 0.5), std::invalid_argument);

    VGParams p_zero_sigma{0.0, 0.5, -0.1};
    EXPECT_THROW(validate_vg_params(p_zero_sigma), std::invalid_argument);

    VGParams p_zero_nu{0.2, 0.0, -0.1};
    EXPECT_THROW(validate_vg_params(p_zero_nu), std::invalid_argument);
}

TEST(VGCFTest, CFMatchesExistingMakeVGCf) {
    // vg_characteristic_function 应与 characteristic_functions.hpp 的 make_vg_cf 一致
    Real S0 = 100.0, r = 0.03, q = 0.0, T = 1.0;
    VGParams p{0.2, 0.5, -0.1};
    auto phi_factory = make_vg_cf(S0, r, q, p.sigma, p.nu, p.theta, T);
    for (Real u_val : {0.1, 0.5, 1.0, 2.0, 5.0}) {
        Complex u(u_val, 0);
        Complex phi_direct = vg_characteristic_function(u, T, S0, r, q, p);
        Complex phi_from_factory = phi_factory(u);
        EXPECT_NEAR(std::real(phi_direct), std::real(phi_from_factory), 1e-12)
            << "Real mismatch at u=" << u_val;
        EXPECT_NEAR(std::imag(phi_direct), std::imag(phi_from_factory), 1e-12)
            << "Imag mismatch at u=" << u_val;
    }
}

TEST(VGCFTest, NoBranchCutDiscontinuity) {
    VGParams p{0.2, 0.5, -0.1};
    Real tau = 1.0;
    const Real eps = 1e-8;
    for (Real u_val = 0.1; u_val <= 10.0; u_val += 0.1) {
        Complex phi_fwd = vg_characteristic_function(Complex(u_val + eps, 0), tau, 100.0, 0.03, 0.0, p);
        Complex phi_bwd = vg_characteristic_function(Complex(u_val - eps, 0), tau, 100.0, 0.03, 0.0, p);
        Complex phi_mid = vg_characteristic_function(Complex(u_val, 0), tau, 100.0, 0.03, 0.0, p);
        Real smoothness = std::abs(phi_fwd + phi_bwd - Real(2) * phi_mid);
        EXPECT_LT(smoothness, 1e-6)
            << "Discontinuity at u=" << u_val << ", smoothness=" << smoothness;
    }
}

TEST(VGCFTest, CumulantsCorrect) {
    // 验证 VG 累积量公式
    VGParams p{0.2, 0.5, -0.1};
    Real T = 1.0;

    // E[X_T] = θT
    EXPECT_NEAR(vg_cumulant_mean(T, p), p.theta * T, 1e-15);

    // Var[X_T] = (σ² + θ²ν)T
    Real expected_var = (p.sigma * p.sigma + p.theta * p.theta * p.nu) * T;
    EXPECT_NEAR(vg_cumulant_variance(T, p), expected_var, 1e-15);

    // θ = 0 时偏度 = 0 (对称 VG)
    VGParams p_sym{0.2, 0.5, 0.0};
    EXPECT_NEAR(vg_cumulant_skewness(T, p_sym), 0.0, 1e-15);

    // 超额峰度 > 0 (VG 有厚尾)
    EXPECT_GT(vg_cumulant_kurtosis_excess(T, p), 0.0);
    // θ = 0 时超额峰度 = 3σ⁴νT / (σ²T)² = 3ν/T > 0
    Real expected_kurt_excess_sym = 3.0 * p_sym.nu / T;
    EXPECT_NEAR(vg_cumulant_kurtosis_excess(T, p_sym), expected_kurt_excess_sym, 1e-15);
}

// ============ 2. VG 随机过程 ============

TEST(VGProcessTest, InvalidParamsThrows) {
    VGParams p_bad{0.5, 2.0, 0.5};  // Feller 违反
    EXPECT_THROW((VarianceGammaProcess{p_bad, 100.0}), std::invalid_argument);

    VGParams p_neg_S0{0.2, 0.5, -0.1};
    EXPECT_THROW((VarianceGammaProcess{p_neg_S0, -100.0}), std::invalid_argument);
    EXPECT_THROW((VarianceGammaProcess{p_neg_S0, 0.0}), std::invalid_argument);
}

TEST(VGProcessTest, PathStartsAtSpot) {
    VGParams p{0.2, 0.5, -0.1};
    VarianceGammaProcess process(p, 100.0, 0.03, 0.0);
    Philox4x64 rng(42);
    std::vector<Real> path(101);
    process.generate_path(1.0, 100, path, rng);
    EXPECT_DOUBLE_EQ(path[0], 100.0);
}

TEST(VGProcessTest, PathDimensionCorrect) {
    VGParams p{0.2, 0.5, -0.1};
    VarianceGammaProcess process(p, 100.0);
    EXPECT_EQ(process.dimension(), 1u);
}

TEST(VGProcessTest, DeterministicWithSameSeed) {
    VGParams p{0.2, 0.5, -0.1};
    VarianceGammaProcess process(p, 100.0, 0.03, 0.0);
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

TEST(VGProcessTest, PositiveSpotAlways) {
    VGParams p{0.3, 0.5, -0.2};
    VarianceGammaProcess process(p, 100.0, 0.03, 0.0);
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

TEST(VGProcessTest, GammaIncrementMean) {
    // Gamma(dt/ν, ν) 的样本均值应 ≈ dt (因为 shape*scale = (dt/ν)*ν = dt)
    VGParams p{0.2, 0.5, -0.1};
    VarianceGammaProcess process(p, 100.0);
    Real dt = 0.01;
    Real expected_mean = dt;  // shape*scale = dt
    Real sum = 0.0;
    int N = 100000;
    Philox4x64 rng(42);
    for (int i = 0; i < N; ++i) {
        sum += process.sample_gamma_increment(dt, rng);
    }
    Real mean = sum / N;
    // Gamma(dt/ν, ν) 的方差 = dt*ν, SE = sqrt(dt*ν/N) ≈ 0.0022
    EXPECT_NEAR(mean, expected_mean, 0.01)
        << "Gamma increment mean: expected=" << expected_mean << " got=" << mean;
}

TEST(VGProcessTest, MCMomentsMatchCumulants) {
    // MC 采样终端 log-return, 验证均值和方差匹配累积量
    VGParams p{0.2, 0.5, -0.1};
    Real S0 = 100.0, r = 0.03, q = 0.0, T = 1.0;
    VarianceGammaProcess process(p, S0, r, q);

    // 理论 log-return 矩
    Real omega = vg_omega(p.sigma, p.nu, p.theta);
    Real expected_mean_lr = std::log(S0) + (r - q + omega) * T + vg_cumulant_mean(T, p);
    Real expected_var_lr = vg_cumulant_variance(T, p);

    Size n_paths = 50000;
    Size n_steps = 1;  // 直接采样终端值 (VG 有精确增量)
    Real sum_lr = 0.0, sum_lr2 = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 7919 + 13));
        std::vector<Real> path(n_steps + 1);
        process.generate_path(T, n_steps, path, rng);
        Real lr = std::log(path.back());
        sum_lr += lr;
        sum_lr2 += lr * lr;
    }
    Real mean_lr = sum_lr / n_paths;
    Real var_lr = sum_lr2 / n_paths - mean_lr * mean_lr;

    // MC ±5% 容差 (log-return 厚尾, 需要较大容差)
    EXPECT_NEAR(mean_lr, expected_mean_lr, 0.05)
        << "Log-return mean: expected=" << expected_mean_lr << " got=" << mean_lr;
    EXPECT_NEAR(var_lr, expected_var_lr, 0.1 * expected_var_lr)
        << "Log-return var: expected=" << expected_var_lr << " got=" << var_lr;
}

// ============ 3. VG COS 定价 ============

TEST(VGPricingTest, CallPricePositive) {
    Real S0 = 100.0, r = 0.03, q = 0.0, T = 1.0;
    VGParams p{0.2, 0.5, -0.1};
    auto phi = make_vg_cf_direct(S0, r, q, p, T);
    COSEngine engine(phi, S0, r, q, T);
    for (Real K : {80.0, 90.0, 100.0, 110.0, 120.0}) {
        Real C = engine.price_call(K);
        EXPECT_GT(C, 0.0) << "Call price should be positive at K=" << K;
        EXPECT_LT(C, S0) << "Call price should be less than spot at K=" << K;
    }
}

TEST(VGPricingTest, CallPutParityHolds) {
    // VG Call-Put Parity: C - P = S e^{-qT} - K e^{-rT}
    // (ω 鞅修正保证 E[S_T] = S_0 e^{(r-q)T})
    Real S0 = 100.0, r = 0.03, q = 0.0, T = 1.0;
    VGParams p{0.2, 0.5, -0.1};
    auto phi = make_vg_cf_direct(S0, r, q, p, T);
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

TEST(VGPricingTest, NuSmallApproachesBSM) {
    // ν → 0 且 θ = 0 时, VG 退化为 BS (vol = σ)
    // 此时 ω → -σ²/2, drift = r-q-σ²/2, X_T = σW(T)
    Real S0 = 100.0, r = 0.03, q = 0.0, T = 1.0, K = 100.0;
    Real sigma = 0.2;
    VGParams p{sigma, 0.001, 0.0};  // ν 很小, θ=0

    auto phi = make_vg_cf_direct(S0, r, q, p, T);
    COSEngine engine(phi, S0, r, q, T);
    Real vg_call = engine.price_call(K);
    Real bs_call = bsm_call_price(S0, K, T, r, q, sigma);

    // ν=0.001 时 VG 应接近 BS (容差 0.05)
    EXPECT_NEAR(vg_call, bs_call, 0.05)
        << "VG (nu small) should approach BSM: VG=" << vg_call << " BSM=" << bs_call;
}

TEST(VGPricingTest, ThetaZeroSymmetricSmile) {
    // θ = 0: 对称 VG, IV 微笑对称中心在 K* = F·exp(ωT) (非 F, 因 omega 鞅修正)
    // log(S_T/F) = ωT + X_T, X_T 对称于 0 ⇒ 分布对称于 ωT ⇒ 微笑对称于 log(K/F)=ωT
    Real S0 = 100.0, r = 0.0, q = 0.0, T = 1.0;
    Real F = S0 * std::exp((r - q) * T);  // = 100
    VGParams p{0.2, 0.5, 0.0};  // θ=0 对称
    Real omega = vg_omega(p.sigma, p.nu, p.theta);
    Real K_star = F * std::exp(omega * T);  // 微笑对称中心

    auto phi = make_vg_cf_direct(S0, r, q, p, T);
    COSEngine engine(phi, S0, r, q, T);

    // log-moneyness 对称于 K*: K1 * K2 = K_star^2
    Real K1 = 90.0;
    Real K2 = K_star * K_star / K1;
    Real C1 = engine.price_call(K1);
    Real C2 = engine.price_call(K2);
    Real iv1 = bsm_implied_vol(C1, S0, K1, T, r, q, true);
    Real iv2 = bsm_implied_vol(C2, S0, K2, T, r, q, true);
    // IV 微笑对称 (容差 0.8% vol, 含 COS 数值误差)
    EXPECT_NEAR(iv1, iv2, 0.008)
        << "Symmetric VG IV smile: IV(K1=" << K1 << ")=" << iv1
        << " IV(K2=" << K2 << ")=" << iv2;
}

TEST(VGPricingTest, MCMatchesCOS) {
    // MC vs COS 半解析定价
    Real S0 = 100.0, r = 0.03, q = 0.0, T = 1.0, K = 100.0;
    VGParams p{0.2, 0.5, -0.1};
    VarianceGammaProcess process(p, S0, r, q);

    // COS 半解析
    auto phi = make_vg_cf_direct(S0, r, q, p, T);
    COSEngine engine(phi, S0, r, q, T);
    Real cos_price = engine.price_call(K);

    // MC (VG 有精确增量, 用 1 步直接采样终端值)
    Size n_paths = 50000;
    Size n_steps = 1;  // VG 是 Levy 过程, 1 步精确
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

    // MC ±4σ 应覆盖 COS 解 (VG 1 步精确, 无离散偏差)
    Real diff = std::abs(mc_price - cos_price);
    EXPECT_LT(diff, 4.0 * se)
        << "MC=" << mc_price << " ± " << se
        << ", COS=" << cos_price << ", diff=" << diff;
}
