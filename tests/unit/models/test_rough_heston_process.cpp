#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <complex>
#include "cpphub/models/diffusion/rough_heston.hpp"
#include "cpphub/models/diffusion/heston.hpp"
#include "cpphub/pricing/fourier/cos_method.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"  // bsm_call_price / bsm_put_price

using namespace cpphub;

// H = 0.5 (α = 1) 退化基准参数 (Feller 条件满足, 几乎不发生截断)
static RoughHestonParams make_degenerate_params() {
    return RoughHestonParams{0.5, 2.0, 0.04, 0.2, -0.5, 0.04, 100.0, 0.0, 0.0};
}

// 二分法 BSM 隐含波动率 (avoid calibrator.hpp's conflicting HestonParams)
static Real bsm_implied_vol_local(Real price, Real S, Real K, Real T,
                                  Real r, Real q, bool is_call) {
    auto bs = [&](Real sigma) {
        return is_call ? bsm_call_price(S, K, T, r, q, sigma)
                       : bsm_put_price(S, K, T, r, q, sigma);
    };
    Real lo = 1e-4, hi = 5.0;
    for (int i = 0; i < 80; ++i) {
        const Real mid = 0.5 * (lo + hi);
        if (bs(mid) > price) hi = mid; else lo = mid;
    }
    return 0.5 * (lo + hi);
}

// 标准 Heston Euler (Full Truncation) MC 定价, 与 rough MC 使用相同 seed 序列
static Real heston_euler_mc_price(const HestonParams& hp, Real T, Real K, bool is_call,
                                  Size n_steps, Size n_paths, uint64_t seed, Real& se_out) {
    Heston h(hp, HestonScheme::FullTruncation);
    Real sum = 0.0, sum2 = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(seed + j);
        std::vector<Real> path(n_steps + 1);
        h.generate_path(T, n_steps, path, rng);
        const Real ST = path.back();
        const Real payoff = is_call ? std::max(ST - K, 0.0) : std::max(K - ST, 0.0);
        sum += payoff;
        sum2 += payoff * payoff;
    }
    const Real mean = sum / static_cast<Real>(n_paths);
    Real var = sum2 / static_cast<Real>(n_paths) - mean * mean;
    if (var < 0.0) var = 0.0;
    se_out = std::exp(-hp.r * T) * std::sqrt(var / static_cast<Real>(n_paths));
    return std::exp(-hp.r * T) * mean;
}

// 1. 核矩阵: 对角元 K[j][j] = Δt^α/α, 上三角 (i>j) 为 0
TEST(RoughHestonKernel, DiagonalAndUpperZero) {
    const Real T = 1.0;
    const Size N = 50;
    const Real alpha = 0.8;
    const Real dt = T / static_cast<Real>(N);
    auto K = rough_heston_kernel(T, N, alpha);
    ASSERT_EQ(K.size(), N);
    for (Size j = 0; j < N; ++j) {
        EXPECT_NEAR(K[j][j], std::pow(dt, alpha) / alpha, 1e-12)
            << "diagonal mismatch at j = " << j;
        for (Size i = j + 1; i < N; ++i) {
            EXPECT_DOUBLE_EQ(K[j][i], 0.0)
                << "upper triangle not zero at (j,i) = (" << j << "," << i << ")";
        }
    }
}

// 2. 核积分总和: Σᵢ K[N-1][i] = T^α/α
TEST(RoughHestonKernel, IntegralSumEqualsTPowAlphaOverAlpha) {
    const Real T = 1.0;
    const Size N = 64;
    const Real alpha = 0.75;
    auto K = rough_heston_kernel(T, N, alpha);
    Real sum = 0.0;
    for (Size i = 0; i < N; ++i) sum += K[N - 1][i];
    EXPECT_NEAR(sum, std::pow(T, alpha) / alpha, 1e-12)
        << "kernel integral sum = " << sum;
}

// 3. α = 1 退化: 核为常数 Δt
TEST(RoughHestonKernel, AlphaOneConstantDt) {
    const Real T = 1.0;
    const Size N = 32;
    const Real alpha = 1.0;
    const Real dt = T / static_cast<Real>(N);
    auto K = rough_heston_kernel(T, N, alpha);
    for (Size j = 0; j < N; ++j) {
        for (Size i = 0; i <= j; ++i) {
            EXPECT_NEAR(K[j][i], dt, 1e-12)
                << "K[" << j << "][" << i << "] = " << K[j][i];
        }
    }
}

// 4. Full truncation 保证方差非负 (1000 条路径)
TEST(RoughHestonProcess, VarianceStaysNonNegative) {
    RoughHestonParams p{0.3, 1.0, 0.04, 0.5, -0.5, 0.04, 100.0, 0.0, 0.0};
    RoughHestonProcess proc(p, 1.0, 100);
    for (Size trial = 0; trial < 1000; ++trial) {
        Philox4x64 rng(static_cast<uint64_t>(1000 + trial * 7));
        auto paths = proc.generate_path(rng, 1.0);
        for (Real v : paths[1]) {
            EXPECT_GE(v, 0.0) << "negative variance at trial " << trial;
            EXPECT_TRUE(std::isfinite(v));
        }
    }
}

// 5. H = 0.5 (α=1) 退化为标准 Heston Euler, 价格一致 (容差 2×std_error)
TEST(RoughHestonProcess, AlphaOneMatchesHestonEuler) {
    auto p = make_degenerate_params();
    const Size n_steps = 52, n_paths = 10000;
    const Real T = 1.0, K = 100.0;
    const uint64_t seed = 12345;

    auto rough = rough_heston_price_european(p, T, K, true, n_steps, n_paths, seed);

    HestonParams hp{p.S0, p.v0, p.kappa, p.theta, p.sigma, p.rho, p.r, p.q};
    Real h_se = 0.0;
    const Real heston = heston_euler_mc_price(hp, T, K, true, n_steps, n_paths, seed, h_se);

    const Real se_combined =
        std::sqrt(rough.std_error * rough.std_error + h_se * h_se);
    EXPECT_NEAR(rough.price, heston, 2.0 * se_combined)
        << "rough MC = " << rough.price << ", heston Euler MC = " << heston
        << ", combined se = " << se_combined;
}

// 6. H = 0.5 (α=1) MC 价格与 Heston COS (CF) 价格一致 (容差 5%)
TEST(RoughHestonProcess, AlphaOneMCvsHestonCOS) {
    auto p = make_degenerate_params();
    const Size n_steps = 100, n_paths = 50000;
    const Real T = 1.0, K = 100.0;

    auto rough = rough_heston_price_european(p, T, K, true, n_steps, n_paths, 777);

    HestonCFParams hp{p.v0, p.kappa, p.theta, p.sigma, p.rho, p.r, p.q};
    const Real cos_price = cos_call_heston(p.S0, K, T, p.r, p.q, hp, 512, 12.0);

    EXPECT_NEAR(rough.price, cos_price, 0.05 * cos_price)
        << "rough MC = " << rough.price << ", Heston COS = " << cos_price;
}

// 7. 杠杆效应: ρ < 0 时 OTM put 的 IV > OTM call 的 IV (负 skew)
TEST(RoughHestonProcess, LeverageSkew) {
    RoughHestonParams p{0.2, 1.0, 0.04, 0.4, -0.7, 0.04, 100.0, 0.0, 0.0};
    const Real T = 0.5;
    const Size n_steps = 100, n_paths = 60000;
    const Real K_put = 90.0, K_call = 110.0;
    const uint64_t seed = 4242;

    auto put_res = rough_heston_price_european(p, T, K_put, false, n_steps, n_paths, seed);
    auto call_res = rough_heston_price_european(p, T, K_call, true, n_steps, n_paths, seed);

    const Real iv_put = bsm_implied_vol_local(put_res.price, p.S0, K_put, T, p.r, p.q, false);
    const Real iv_call = bsm_implied_vol_local(call_res.price, p.S0, K_call, T, p.r, p.q, true);

    EXPECT_GT(iv_put, iv_call)
        << "iv_put = " << iv_put << ", iv_call = " << iv_call
        << ", put price = " << put_res.price << " (se " << put_res.std_error << ")"
        << ", call price = " << call_res.price << " (se " << call_res.std_error << ")";
}

// 8. 路径确定性: 相同 seed 生成完全相同路径 (可复现性)
TEST(RoughHestonProcess, DeterministicWithSameSeed) {
    auto p = make_degenerate_params();
    RoughHestonProcess proc(p, 1.0, 50);
    std::vector<std::vector<Real>> path1, path2;
    {
        Philox4x64 rng(42);
        path1 = proc.generate_path(rng, 1.0);
    }
    {
        Philox4x64 rng(42);
        path2 = proc.generate_path(rng, 1.0);
    }
    for (Size i = 0; i <= 50; ++i) {
        EXPECT_DOUBLE_EQ(path1[0][i], path2[0][i])
            << "S mismatch at step " << i;
        EXPECT_DOUBLE_EQ(path1[1][i], path2[1][i])
            << "v mismatch at step " << i;
    }
}

// 9. 价格鞅性: E[S_T] = S0·exp((r-q)·T) (容差 1%)
TEST(RoughHestonProcess, Martingale) {
    RoughHestonParams p{0.3, 1.5, 0.04, 0.3, -0.5, 0.04, 100.0, 0.05, 0.02};
    const Size n_steps = 52, n_paths = 50000;
    const Real T = 1.0;
    RoughHestonProcess proc(p, T, n_steps);

    Real sum_ST = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 12345 + 7));
        auto paths = proc.generate_path(rng, 1.0);
        sum_ST += paths[0].back();
    }
    const Real mean_ST = sum_ST / static_cast<Real>(n_paths);
    const Real expected = p.S0 * std::exp((p.r - p.q) * T);
    const Real rel = std::abs(mean_ST - expected) / expected;
    EXPECT_LT(rel, 0.01)
        << "mean_ST = " << mean_ST << ", expected = " << expected
        << ", rel_error = " << rel;
}

// 10. 数值稳定性: 极端粗糙 (H=0.05) + 长期限 (T=2) 不产生 NaN/Inf
TEST(RoughHestonProcess, NumericalStabilityExtremeRoughness) {
    RoughHestonParams p{0.05, 2.0, 0.04, 0.4, -0.6, 0.04, 100.0, 0.0, 0.0};
    const Real T = 2.0;
    const Size n_steps = 200;
    RoughHestonProcess proc(p, T, n_steps);
    for (Size trial = 0; trial < 200; ++trial) {
        Philox4x64 rng(static_cast<uint64_t>(trial * 31 + 1));
        auto paths = proc.generate_path(rng, 1.0);
        for (Real s : paths[0]) {
            EXPECT_TRUE(std::isfinite(s) && s > 0.0)
                << "non-finite S at trial " << trial;
        }
        for (Real v : paths[1]) {
            EXPECT_TRUE(std::isfinite(v) && v >= 0.0)
                << "non-finite v at trial " << trial;
        }
    }
}

// 附加: 接口契约 (dimension / spot / alpha / 路径起点)
TEST(RoughHestonProcess, DimensionAndSpotContract) {
    auto p = make_degenerate_params();
    RoughHestonProcess proc(p, 1.0, 52);
    EXPECT_EQ(proc.dimension(), 2u);
    EXPECT_DOUBLE_EQ(proc.spot(), 100.0);
    EXPECT_DOUBLE_EQ(proc.alpha(), 1.0);
    EXPECT_DOUBLE_EQ(proc.T(), 1.0);
    EXPECT_EQ(proc.n_steps(), 52u);

    Philox4x64 rng(1);
    std::vector<Real> path(2 * 53);
    proc.generate_path(1.0, 52, path, rng);
    EXPECT_DOUBLE_EQ(path[0], 100.0);          // S0
    EXPECT_DOUBLE_EQ(path[53], 0.04);          // v0 (双维展平第二段起点)
}

// 附加: 不同 seed 生成不同路径 (随机性)
TEST(RoughHestonProcess, DifferentSeedDifferentPath) {
    auto p = make_degenerate_params();
    RoughHestonProcess proc(p, 1.0, 25);
    std::vector<Real> s1, s2;
    {
        Philox4x64 rng(1);
        auto paths = proc.generate_path(rng, 1.0);
        s1 = paths[0];
    }
    {
        Philox4x64 rng(2);
        auto paths = proc.generate_path(rng, 1.0);
        s2 = paths[0];
    }
    bool any_diff = false;
    for (Size i = 0; i <= 25; ++i) {
        if (std::abs(s1[i] - s2[i]) > 1e-12) { any_diff = true; break; }
    }
    EXPECT_TRUE(any_diff);
}
