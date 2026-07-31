#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/models/diffusion/rough_bergomi.hpp"
#include "cpphub/models/diffusion/rbergomi_hybrid_scheme.hpp"

using namespace cpphub::v1;

namespace {

// 用 Philox4x64 + box_muller 生成 n 个独立 N(0,1) (与 rough_bergomi.hpp 一致的取法)
std::vector<Real> make_normal(Philox4x64& rng, Size n) {
    std::vector<Real> z(n, 0.0);
    for (Size i = 0; i < n; i += 2) {
        uint64_t r1 = rng();
        uint64_t r2 = rng();
        double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
        double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
        auto [z1, z2] = box_muller(u1, u2);
        z[i] = z1;
        if (i + 1 < n) z[i + 1] = z2;
    }
    return z;
}

// 最大相对误差 (忽略 |B| < 1e-6 * max|B| 的微小元素, 避免无意义的大相对误差)
Real max_rel_error(const std::vector<std::vector<Real>>& A,
                   const std::vector<std::vector<Real>>& B) {
    const Size n = B.size();
    Real m = 0.0;
    for (const auto& row : B)
        for (Real v : row) m = std::max(m, std::abs(v));
    const Real floor_v = 1e-6 * m;
    Real mr = 0.0;
    for (Size i = 0; i < n; ++i) {
        for (Size j = 0; j < n; ++j) {
            if (std::abs(B[i][j]) > floor_v) {
                mr = std::max(mr, std::abs(A[i][j] - B[i][j]) / std::abs(B[i][j]));
            }
        }
    }
    return mr;
}

Real max_abs_error(const std::vector<std::vector<Real>>& A,
                   const std::vector<std::vector<Real>>& B) {
    const Size n = B.size();
    Real m = 0.0;
    for (Size i = 0; i < n; ++i)
        for (Size j = 0; j < n; ++j) m = std::max(m, std::abs(A[i][j] - B[i][j]));
    return m;
}

// 与 rough_bergomi.hpp 的 generate_path_with_sampler 完全一致的 log-Euler 定价,
// 但 W̃^H 路径由任意采样器生成 (wh: Z -> W̃^H)。
Real price_with_fbm(const RoughBergomiParams& p, Real K, Real T, Size n_steps,
                    Size n_paths, uint64_t seed,
                    const std::function<std::vector<Real>(const std::vector<Real>&)>& wh) {
    const Real dt = T / static_cast<Real>(n_steps);
    const Real sqrt_dt = std::sqrt(dt);
    const Real eta = p.eta, rho = p.rho, xi0 = p.xi0, H = p.H;
    const Real sqrt_1m_rho2 = std::sqrt(1.0 - rho * rho);

    Real sum_payoff = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(seed + j);
        std::vector<Real> Z(n_steps);
        std::vector<Real> dW_perp(n_steps);
        for (Size i = 0; i < n_steps; ++i) {
            uint64_t r1 = rng();
            uint64_t r2 = rng();
            double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            auto [z1, z2] = box_muller(u1, u2);
            Z[i] = z1;
            dW_perp[i] = z2 * sqrt_dt;
        }

        std::vector<Real> WH = wh(Z);

        Real log_S = std::log(p.S0);
        for (Size i = 0; i < n_steps; ++i) {
            Real v_step;
            if (i == 0) {
                v_step = xi0;
            } else {
                Real t_prev = static_cast<Real>(i) * dt;
                Real t_prev_2H = std::pow(t_prev, 2.0 * H);
                Real log_v_prev = std::log(xi0) + eta * WH[i - 1] - 0.5 * eta * eta * t_prev_2H;
                v_step = std::exp(log_v_prev);
            }
            Real dW_i = Z[i] * sqrt_dt;
            Real dZ = rho * dW_i + sqrt_1m_rho2 * dW_perp[i];
            log_S += -0.5 * v_step * dt + std::sqrt(v_step) * dZ;
        }
        Real ST = std::exp(log_S);
        sum_payoff += std::max(ST - K, 0.0);
    }
    return std::exp(-p.r * T) * sum_payoff / static_cast<Real>(n_paths);
}

}  // namespace

// ============ HybridSchemeTest ============

TEST(HybridSchemeTest, ParameterValidation) {
    // H 越界 / T 非正 / n_steps 为 0 / b 越界 都应抛异常
    EXPECT_THROW(RLFbmHybridSampler(1.0, 64, 0.0), std::invalid_argument);
    EXPECT_THROW(RLFbmHybridSampler(1.0, 64, 0.5), std::invalid_argument);
    EXPECT_THROW(RLFbmHybridSampler(1.0, 64, 1.0), std::invalid_argument);
    EXPECT_THROW(RLFbmHybridSampler(1.0, 64, -0.1), std::invalid_argument);
    EXPECT_THROW(RLFbmHybridSampler(0.0, 64, 0.1), std::invalid_argument);
    EXPECT_THROW(RLFbmHybridSampler(-1.0, 64, 0.1), std::invalid_argument);
    EXPECT_THROW(RLFbmHybridSampler(1.0, 0, 0.1), std::invalid_argument);
    EXPECT_THROW(RLFbmHybridSampler(1.0, 64, 0.1, HybridSchemeConfig{0, false}),
                 std::invalid_argument);
    EXPECT_THROW(RLFbmHybridSampler(1.0, 64, 0.1, HybridSchemeConfig{65, false}),
                 std::invalid_argument);
    // 合法构造
    EXPECT_NO_THROW(RLFbmHybridSampler(1.0, 64, 0.1));
    EXPECT_NO_THROW(RLFbmHybridSampler(1.0, 64, 0.1, HybridSchemeConfig{64, false}));
}

TEST(HybridSchemeTest, CovarianceDiagonal) {
    // Var(W̃^H_{t_i}) = (2H+1)/(2H) * t_i^{2H} (真实解析值)
    // H=0.49 时核平坦, 平顶远端近似误差最小, 容差 10%
    const Real T = 1.0, H = 0.49;
    const Size n = 64;
    RLFbmHybridSampler sampler(T, n, H);  // b=1 默认
    auto C = sampler.implied_covariance();

    const Real dt = T / static_cast<Real>(n);
    for (Size i = 0; i < n; ++i) {
        Real t = static_cast<Real>(i + 1) * dt;
        Real var_true = (2.0 * H + 1.0) / (2.0 * H) * std::pow(t, 2.0 * H);
        EXPECT_LT(std::abs(C[i][i] - var_true) / var_true, 0.10)
            << "diagonal element " << i << " deviates from analytic variance";
    }
}

TEST(HybridSchemeTest, CovarianceOffDiagonal) {
    // 非对角线协方差与 Cholesky 精确值 rbergomi_fbm_covariance 对比, 容差 10%
    const Real T = 1.0, H = 0.49;
    const Size n = 128;
    RLFbmHybridSampler sampler(T, n, H);
    auto implied = sampler.implied_covariance();
    auto exact = rbergomi_fbm_covariance(T, n, H);
    EXPECT_LT(max_rel_error(implied, exact), 0.10);
}

TEST(HybridSchemeTest, PathMeanZero) {
    // MC 验证 E[W̃^H_t] = 0 (10000 路径, 均值 |mean| < 0.05)
    const Real T = 1.0, H = 0.1;
    const Size n = 64, n_paths = 10000;
    RLFbmHybridSampler sampler(T, n, H);

    Real sum = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(42 + j);
        auto Z = make_normal(rng, n);
        auto W = sampler.sample(Z);
        sum += W.back();
    }
    Real mean = sum / static_cast<Real>(n_paths);
    EXPECT_LT(std::abs(mean), 0.05);
}

TEST(HybridSchemeTest, PathVariance) {
    // MC 验证 Var(W̃^H_t) ≈ (2H+1)/(2H) * t^{2H} (10000 路径, 相对误差 < 10%)
    const Real T = 1.0, H = 0.49;
    const Size n = 64, n_paths = 10000;
    RLFbmHybridSampler sampler(T, n, H);

    std::vector<Real> Wt(n_paths);
    Real sum = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(7 + j);
        auto Z = make_normal(rng, n);
        auto W = sampler.sample(Z);
        Wt[j] = W.back();
        sum += Wt[j];
    }
    Real mean = sum / static_cast<Real>(n_paths);
    Real var = 0.0;
    for (Real v : Wt) var += (v - mean) * (v - mean);
    var /= static_cast<Real>(n_paths - 1);

    Real var_true = (2.0 * H + 1.0) / (2.0 * H) * std::pow(T, 2.0 * H);
    EXPECT_LT(std::abs(var - var_true) / var_true, 0.10);
}

TEST(HybridSchemeTest, VsCholeskyStats) {
    // 同 Z 下 Hybrid vs Cholesky 路径统计一致性: 均值差 < 0.1, 方差差 < 10%
    const Real T = 1.0, H = 0.49;
    const Size n = 64, n_paths = 10000;
    RLFbmHybridSampler hybrid(T, n, H);
    RLFbmSampler chol(T, n, H);

    Real sum_h = 0.0, sum_c = 0.0;
    std::vector<Real> Wh(n_paths), Wc(n_paths);
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(100 + j);
        auto Z = make_normal(rng, n);
        Wh[j] = hybrid.sample(Z).back();
        Wc[j] = chol.sample(Z).back();
        sum_h += Wh[j];
        sum_c += Wc[j];
    }
    Real mean_h = sum_h / static_cast<Real>(n_paths);
    Real mean_c = sum_c / static_cast<Real>(n_paths);
    Real var_h = 0.0, var_c = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        var_h += (Wh[j] - mean_h) * (Wh[j] - mean_h);
        var_c += (Wc[j] - mean_c) * (Wc[j] - mean_c);
    }
    var_h /= static_cast<Real>(n_paths - 1);
    var_c /= static_cast<Real>(n_paths - 1);

    EXPECT_LT(std::abs(mean_h - mean_c), 0.1);
    EXPECT_LT(std::abs(var_h - var_c) / var_c, 0.10);
}

TEST(HybridSchemeTest, LargeNPerformance) {
    // N=256 时 Hybrid vs Cholesky 采样时间比 < 0.5 (Hybrid 更快)
    const Real T = 1.0, H = 0.1;
    const Size N = 256, P = 1000;
    RLFbmHybridSampler hybrid(T, N, H);
    RLFbmSampler chol(T, N, H);
    Philox4x64 rng(3);
    auto Z = make_normal(rng, N);

    using clock = std::chrono::steady_clock;
    // 预热
    for (Size p = 0; p < 20; ++p) { volatile auto w1 = hybrid.sample(Z); volatile auto w2 = chol.sample(Z); (void)w1; (void)w2; }

    auto time_sample = [&](auto& sampler) {
        double best = std::numeric_limits<double>::max();
        for (int trial = 0; trial < 3; ++trial) {
            auto t0 = clock::now();
            for (Size p = 0; p < P; ++p) { auto W = sampler.sample(Z); (void)W; }
            auto t1 = clock::now();
            best = std::min(best, std::chrono::duration<double>(t1 - t0).count());
        }
        return best;
    };

    double t_hyb = time_sample(hybrid);
    double t_chol = time_sample(chol);
    double ratio = t_hyb / t_chol;
    EXPECT_LT(ratio, 0.5);
    std::cout << "[ perf ] N=256, P=1000: Hybrid=" << (t_hyb * 1e3) << " ms, "
              << "Cholesky=" << (t_chol * 1e3) << " ms, ratio=" << ratio << std::endl;
}

TEST(HybridSchemeTest, B1Default) {
    // b=1 (默认) 时, 协方差与 Cholesky 差异 < 10%
    const Real T = 1.0, H = 0.49;
    const Size n = 128;
    RLFbmHybridSampler sampler(T, n, H, HybridSchemeConfig{1, false});
    EXPECT_EQ(sampler.b(), 1);
    auto implied = sampler.implied_covariance();
    auto exact = rbergomi_fbm_covariance(T, n, H);
    EXPECT_LT(max_rel_error(implied, exact), 0.10);
}

TEST(HybridSchemeTest, B3Accuracy) {
    // b=3 时, 协方差与 Cholesky 差异 < 5% (b 越大越精确)
    const Real T = 1.0, H = 0.49;
    const Size n = 128;
    RLFbmHybridSampler sampler(T, n, H, HybridSchemeConfig{3, false});
    auto implied = sampler.implied_covariance();
    auto exact = rbergomi_fbm_covariance(T, n, H);
    EXPECT_LT(max_rel_error(implied, exact), 0.05);
}

TEST(HybridSchemeTest, BEqualsN) {
    // b = n_steps 时 Hybrid 退化为全近端, 与 Cholesky 等价 (差异 < 1e-10)
    const Real T = 1.0, H = 0.1;
    const Size n = 64;
    RLFbmHybridSampler sampler(T, n, H, HybridSchemeConfig{n, false});
    auto implied = sampler.implied_covariance();
    auto exact = rbergomi_fbm_covariance(T, n, H);
    EXPECT_LT(max_abs_error(implied, exact), 1e-10);
}

TEST(HybridSchemeTest, ExtremeH005) {
    // H=0.05 (极端粗糙) 数值稳定, 无 NaN/Inf
    const Real T = 1.0, H = 0.05;
    const Size n = 256;
    RLFbmHybridSampler sampler(T, n, H);
    Philox4x64 rng(11);
    auto Z = make_normal(rng, n);
    auto W = sampler.sample(Z);
    for (Real v : W) {
        EXPECT_TRUE(std::isfinite(v));
    }
    auto C = sampler.implied_covariance();
    for (Size i = 0; i < n; ++i) {
        EXPECT_TRUE(std::isfinite(C[i][i]));
        EXPECT_GT(C[i][i], 0.0);
    }
}

TEST(HybridSchemeTest, H049) {
    // H=0.49 (接近标准 Brownian) 时, Hybrid 与 Cholesky 差异 < 5%
    const Real T = 1.0, H = 0.49;
    const Size n = 128;
    RLFbmHybridSampler sampler(T, n, H);
    auto implied = sampler.implied_covariance();
    auto exact = rbergomi_fbm_covariance(T, n, H);
    EXPECT_LT(max_rel_error(implied, exact), 0.05);
}

TEST(HybridSchemeTest, Determinism) {
    // 同 seed 同 Z 序列 → 同路径 (位精确)
    const Real T = 1.0, H = 0.1;
    const Size n = 64;
    RLFbmHybridSampler sampler(T, n, H);
    Philox4x64 rng(2024);
    auto Z = make_normal(rng, n);
    auto W1 = sampler.sample(Z);
    auto W2 = sampler.sample(Z);
    ASSERT_EQ(W1.size(), W2.size());
    for (Size i = 0; i < n; ++i) EXPECT_EQ(W1[i], W2[i]);
}

TEST(HybridSchemeTest, PathLength) {
    // 输出路径长度 = n_steps
    const Real T = 1.0, H = 0.2;
    const Size n = 37;
    RLFbmHybridSampler sampler(T, n, H);
    Philox4x64 rng(5);
    auto Z = make_normal(rng, n);
    auto W = sampler.sample(Z);
    EXPECT_EQ(W.size(), n);
    EXPECT_EQ(sampler.n_steps(), n);
}

TEST(HybridSchemeTest, IntegrationWithBergomi) {
    // 用 RLFbmHybridSampler 替换 RLFbmSampler 组装 rBergomi 路径,
    // MC 定价与 Cholesky 版一致 (容差 5%)。同 seed → 同 Z, 差异仅来自采样器离散化。
    const Real T = 1.0, H = 0.1;
    const Size n_steps = 64, n_paths = 3000;
    const uint64_t seed = 12345;

    RoughBergomiParams p;
    p.H = 0.1; p.eta = 0.3; p.rho = -0.7; p.xi0 = 0.04; p.S0 = 100.0; p.r = 0.03; p.q = 0.0;

    RLFbmHybridSampler hybrid(T, n_steps, H);
    RLFbmSampler chol(T, n_steps, H);

    Real price_h = price_with_fbm(p, 100.0, T, n_steps, n_paths, seed,
                                  [&](const std::vector<Real>& Z) { return hybrid.sample(Z); });
    Real price_c = price_with_fbm(p, 100.0, T, n_steps, n_paths, seed,
                                  [&](const std::vector<Real>& Z) { return chol.sample(Z); });
    std::cout << "[ integ ] ATM call: Hybrid=" << price_h << ", Cholesky=" << price_c
              << ", rel diff=" << std::abs(price_h - price_c) / price_c << std::endl;
    EXPECT_LT(std::abs(price_h - price_c) / price_c, 0.05);
}

TEST(HybridSchemeTest, ImpliedCovarianceMatchesMC) {
    // sample() 的统计行为与 implied_covariance() 一致 (MC 交叉验证)
    const Real T = 1.0, H = 0.49;
    const Size n = 16, n_paths = 20000;
    RLFbmHybridSampler sampler(T, n, H);
    auto implied = sampler.implied_covariance();

    std::vector<std::vector<Real>> acc(n, std::vector<Real>(n, 0.0));
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(99 + j);
        auto Z = make_normal(rng, n);
        auto W = sampler.sample(Z);
        for (Size i = 0; i < n; ++i)
            for (Size k = 0; k <= i; ++k) {
                acc[i][k] += W[i] * W[k];
            }
    }
    Real maxrel = 0.0;
    for (Size i = 0; i < n; ++i)
        for (Size k = 0; k <= i; ++k) {
            Real mc = acc[i][k] / static_cast<Real>(n_paths);
            maxrel = std::max(maxrel, std::abs(mc - implied[i][k]) / std::abs(implied[i][k]));
        }
    EXPECT_LT(maxrel, 0.15);  // MC 噪声容差
}

TEST(HybridSchemeTest, KernelStructure) {
    // 核系数结构与公式一致: 近端 = sqrt(2H+1)*sqrt(dt)*(m*dt)^alpha,
    // 远端平顶权重 = sqrt(2H+1)*sqrt(dt)*(max(i-b,1)*dt)^alpha
    const Real T = 1.0, H = 0.2;
    const Size n = 32;
    const Size b = 3;
    RLFbmHybridSampler sampler(T, n, H, HybridSchemeConfig{b, false});
    const Real dt = T / static_cast<Real>(n);
    const Real alpha = H - 0.5;
    const Real sc = std::sqrt(2.0 * H + 1.0) * std::sqrt(dt);

    for (Size m = 1; m <= b; ++m) {
        Real expect = sc * std::pow(static_cast<Real>(m) * dt, alpha);
        EXPECT_NEAR(sampler.near_kernel(m), expect, 1e-15);
    }
    for (Size i = b; i < n; ++i) {
        // 实现用 max(i-b, 1) 规避 0^alpha = inf
        Size j = std::max(Size{1}, i - b);
        Real expect = sc * std::pow(static_cast<Real>(j) * dt, alpha);
        EXPECT_NEAR(sampler.far_weight(i), expect, 1e-15);
    }
    EXPECT_EQ(sampler.far_weight(b - 1), 0.0);  // 远端为空
    EXPECT_EQ(sampler.T(), T);
    EXPECT_EQ(sampler.n_steps(), n);
    EXPECT_EQ(sampler.H(), H);
    EXPECT_EQ(sampler.b(), b);
    EXPECT_FALSE(sampler.config().use_fft);
}

TEST(HybridSchemeTest, SymmetricPositiveDefinite) {
    // 隐含协方差对称且对角线为正
    const Real T = 1.0, H = 0.3;
    const Size n = 32;
    RLFbmHybridSampler sampler(T, n, H);
    auto C = sampler.implied_covariance();
    for (Size i = 0; i < n; ++i) {
        EXPECT_GT(C[i][i], 0.0);
        for (Size j = 0; j < n; ++j) EXPECT_DOUBLE_EQ(C[i][j], C[j][i]);
    }
}
