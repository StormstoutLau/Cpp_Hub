#include <gtest/gtest.h>
// 路径相关解析解测试 (Kemna-Vorst / Reiner-Rubinstein / Goldman-Sosin-Gatto)
// 验证策略:
//   1. 解析公式内部一致性 (put-call parity, barrier parity)
//   2. 与 MC 模拟数值对照 (容差取决于 MC 路径数)
//   3. 边界情况 (H→S, m=S, T→0)

#include "cpphub/pricing/analytic/path_dependent_analytic.hpp"
#include "cpphub/pricing/monte_carlo/path_generator.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/pricing/monte_carlo/path_dependent_payoffs.hpp"
#include "cpphub/pricing/monte_carlo/variance_reduction.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/core/math.hpp"
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace cpphub::v1;

namespace {

// BSM 闭式解 (用于 barrier parity 验证)
Real bsm_call(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(S - K, 0.0);
    Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    Real d2 = d1 - sigma * std::sqrt(T);
    return S * std::exp(-q * T) * normal_cdf(d1) - K * std::exp(-r * T) * normal_cdf(d2);
}

Real bsm_put(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(K - S, 0.0);
    Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    Real d2 = d1 - sigma * std::sqrt(T);
    return K * std::exp(-r * T) * normal_cdf(-d2) - S * std::exp(-q * T) * normal_cdf(-d1);
}

// MC 障碍期权定价 (用于数值验证)
Real mc_barrier(Real S, Real K, Real H, Real T, Real r, Real q, Real sigma,
                BarrierDirection dir, BarrierKnock knock, OptionType opt,
                Size n_paths = 500000, uint64_t seed = 42) {
    auto cfg = make_single_asset_gbm(S, sigma, r, q, T, 200);
    MultiAssetGBMPathGenerator gen(cfg);
    Philox4x64 rng(seed);
    Real df = std::exp(-r * T);
    Real sum = 0.0, sum2 = 0.0;
    for (Size i = 0; i < n_paths; ++i) {
        auto paths = gen.generate_path(rng);
        const auto& path = paths[0];
        bool knocked = false;
        for (Size j = 0; j < path.size(); ++j) {
            if (dir == BarrierDirection::Up && path[j] >= H) { knocked = true; break; }
            if (dir == BarrierDirection::Down && path[j] <= H) { knocked = true; break; }
        }
        bool active = (knock == BarrierKnock::Out) ? !knocked : knocked;
        if (!active) { sum += 0.0; continue; }
        Real ST = path.back();
        Real payoff = (opt == OptionType::Call) ? std::max(ST - K, 0.0) : std::max(K - ST, 0.0);
        sum += payoff;
        sum2 += payoff * payoff;
    }
    Real price = df * sum / n_paths;
    return price;
}

// MC Lookback 浮动 call 定价 (payoff = S_T - m_T)
// 使用 Brownian Bridge 极值采样: 在 [t_i, t_{i+1}] 间已知端点时,
// P(min > m | X_i, X_{i+1}) = 1 - exp(-2*(X_i-m)*(X_{i+1}-m)/(σ²Δt))
// 反变换: y = (-z + sqrt(z² + 2σ²Δt*(-log(u)))) / 2, m = X_i - y, z = X_{i+1}-X_i
Real mc_lookback_call_float(Real S, Real T, Real r, Real q, Real sigma,
                             Size n_paths = 500000, uint64_t seed = 42) {
    auto cfg = make_single_asset_gbm(S, sigma, r, q, T, 50);  // 50 步 + BB 精确化
    MultiAssetGBMPathGenerator gen(cfg);
    Philox4x64 rng(seed);
    Real df = std::exp(-r * T);
    Real dt = T / 50.0;
    Real sigma_sq_dt = sigma * sigma * dt;
    Real sum = 0.0;
    for (Size i = 0; i < n_paths; ++i) {
        auto paths = gen.generate_path(rng);
        const auto& path = paths[0];
        Real m = path[0];
        for (Size j = 0; j + 1 < path.size(); ++j) {
            Real Xi = std::log(path[j]);
            Real Xj = std::log(path[j + 1]);
            Real z = Xj - Xi;
            // u ~ U(0,1) via Philox
            Real u = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
            if (u < 1e-15) u = 1e-15;
            Real disc = z * z + 2.0 * sigma_sq_dt * (-std::log(u));
            if (disc < 0.0) disc = 0.0;
            Real y = 0.5 * (-z + std::sqrt(disc));  // y ≥ 0 (取正根)
            Real X_min_bb = Xi - y;
            Real S_min_bb = std::exp(X_min_bb);
            if (S_min_bb < m) m = S_min_bb;
            if (path[j + 1] < m) m = path[j + 1];
        }
        sum += (path.back() - m);
    }
    return df * sum / n_paths;
}

// MC Lookback 浮动 put 定价 (payoff = M_T - S_T) — Brownian Bridge 改进
// 最大值采样: M = -min(-X), 对称变换 y_max = (z + sqrt(z² + 2σ²Δt*(-log(u))))/2, M = X_i + y_max
Real mc_lookback_put_float(Real S, Real T, Real r, Real q, Real sigma,
                            Size n_paths = 500000, uint64_t seed = 42) {
    auto cfg = make_single_asset_gbm(S, sigma, r, q, T, 50);
    MultiAssetGBMPathGenerator gen(cfg);
    Philox4x64 rng(seed);
    Real df = std::exp(-r * T);
    Real dt = T / 50.0;
    Real sigma_sq_dt = sigma * sigma * dt;
    Real sum = 0.0;
    for (Size i = 0; i < n_paths; ++i) {
        auto paths = gen.generate_path(rng);
        const auto& path = paths[0];
        Real M = path[0];
        for (Size j = 0; j + 1 < path.size(); ++j) {
            Real Xi = std::log(path[j]);
            Real Xj = std::log(path[j + 1]);
            Real z = Xj - Xi;
            Real u = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
            if (u < 1e-15) u = 1e-15;
            Real disc = z * z + 2.0 * sigma_sq_dt * (-std::log(u));
            if (disc < 0.0) disc = 0.0;
            Real y_max = 0.5 * (z + std::sqrt(disc));  // y_max ≥ 0
            Real X_max_bb = Xi + y_max;
            Real S_max_bb = std::exp(X_max_bb);
            if (S_max_bb > M) M = S_max_bb;
            if (path[j + 1] > M) M = path[j + 1];
        }
        sum += (M - path.back());
    }
    return df * sum / n_paths;
}

// MC 几何亚式 call 定价
Real mc_geom_asian_call(Real S0, Real K, Real T, Real r, Real q, Real sigma, Size n_steps,
                         Size n_paths = 500000, uint64_t seed = 42) {
    auto cfg = make_single_asset_gbm(S0, sigma, r, q, T, n_steps);
    MultiAssetGBMPathGenerator gen(cfg);
    Philox4x64 rng(seed);
    Real df = std::exp(-r * T);
    Real sum = 0.0;
    for (Size i = 0; i < n_paths; ++i) {
        auto paths = gen.generate_path(rng);
        const auto& path = paths[0];
        Real sum_log = 0.0;
        for (Size j = 1; j <= n_steps; ++j) sum_log += std::log(path[j]);
        Real geom_avg = std::exp(sum_log / n_steps);
        sum += std::max(geom_avg - K, 0.0);
    }
    return df * sum / n_paths;
}

}  // anonymous namespace

// ================================================================
// ============ 1. Geometric Asian (Kemna-Vorst) 测试 ============
// ================================================================

// 1.1 离散 vs 连续: n_steps → ∞ 时离散 → 连续
TEST(GeomAsian, DiscreteConvergesToContinuous) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real cont = geom_asian_call_continuous(S0, K, T, r, q, sigma);
    Real disc_50  = geom_asian_call_price(S0, K, T, r, q, sigma, 50);
    Real disc_500 = geom_asian_call_price(S0, K, T, r, q, sigma, 500);
    Real disc_2000 = geom_asian_call_price(S0, K, T, r, q, sigma, 2000);
    // n=2000 应非常接近连续
    EXPECT_NEAR(disc_2000, cont, 0.05);
    // n 越大越接近连续
    EXPECT_LT(std::abs(disc_500 - cont), std::abs(disc_50 - cont));
}

// 1.2 put-call parity: C - P = e^{-rT}(E[G] - K)
TEST(GeomAsian, PutCallParity) {
    Real S0 = 100.0, K = 105.0, T = 1.0, r = 0.05, q = 0.02, sigma = 0.30;
    Size n = 50;
    Real c = geom_asian_call_price(S0, K, T, r, q, sigma, n);
    Real p = geom_asian_put_price(S0, K, T, r, q, sigma, n);
    // E[G] = exp(mu_g + var_g/2)
    Real b = r - q;
    Real nn = static_cast<Real>(n);
    Real mu_g = std::log(S0) + (b - 0.5 * sigma * sigma) * T * (nn + 1.0) / (2.0 * nn);
    Real var_g = sigma * sigma * T * (nn + 1.0) * (2.0 * nn + 1.0) / (6.0 * nn * nn);
    Real EG = std::exp(mu_g + 0.5 * var_g);
    Real parity = std::exp(-r * T) * (EG - K);
    EXPECT_NEAR(c - p, parity, 1e-10);
}

// 1.3 vs MC (离散几何亚式 call)
TEST(GeomAsian, CallVsMC) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Size n = 50;
    Real analytic = geom_asian_call_price(S0, K, T, r, q, sigma, n);
    Real mc = mc_geom_asian_call(S0, K, T, r, q, sigma, n, 500000);
    // MC 标准误差约 0.02-0.05, 容差 0.15 (3σ)
    EXPECT_NEAR(analytic, mc, 0.15);
}

// 1.4 vs MC (ITM case)
TEST(GeomAsian, CallITMvsMC) {
    Real S0 = 110.0, K = 100.0, T = 0.5, r = 0.03, q = 0.01, sigma = 0.25;
    Size n = 252;
    Real analytic = geom_asian_call_price(S0, K, T, r, q, sigma, n);
    Real mc = mc_geom_asian_call(S0, K, T, r, q, sigma, n, 500000);
    EXPECT_NEAR(analytic, mc, 0.15);
}

// 1.5 连续 put vs MC (大 n 近似)
TEST(GeomAsian, ContinuousPutVsMC) {
    Real S0 = 100.0, K = 95.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real cont = geom_asian_put_continuous(S0, K, T, r, q, sigma);
    // 用大 n 离散 MC 近似连续
    auto cfg = make_single_asset_gbm(S0, sigma, r, q, T, 1000);
    MultiAssetGBMPathGenerator gen(cfg);
    Philox4x64 rng(42);
    Real df = std::exp(-r * T);
    Real sum = 0.0;
    Size n_paths = 500000;
    for (Size i = 0; i < n_paths; ++i) {
        auto paths = gen.generate_path(rng);
        const auto& path = paths[0];
        Real sum_log = 0.0;
        for (Size j = 1; j <= 1000; ++j) sum_log += std::log(path[j]);
        Real geom = std::exp(sum_log / 1000.0);
        sum += std::max(K - geom, 0.0);
    }
    Real mc = df * sum / n_paths;
    EXPECT_NEAR(cont, mc, 0.10);
}

// ================================================================
// ============ 2. Barrier Options (Reiner-Rubinstein) 测试 ============
// ================================================================

// 2.1 Barrier parity: In + Out = Vanilla (所有 8 种类型)
TEST(Barrier, ParityAllTypes) {
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real H_down = 85.0;  // down barrier
    Real H_up = 115.0;   // up barrier

    Real vc = bsm_call(S, K, T, r, q, sigma);
    Real vp = bsm_put(S, K, T, r, q, sigma);

    // Down barriers (H < S)
    Real doc = barrier_option_price(S, K, H_down, T, r, q, sigma, BarrierType::DownOutCall);
    Real dic = barrier_option_price(S, K, H_down, T, r, q, sigma, BarrierType::DownInCall);
    EXPECT_NEAR(doc + dic, vc, 1e-10);

    Real dop = barrier_option_price(S, K, H_down, T, r, q, sigma, BarrierType::DownOutPut);
    Real dip = barrier_option_price(S, K, H_down, T, r, q, sigma, BarrierType::DownInPut);
    EXPECT_NEAR(dop + dip, vp, 1e-10);

    // Up barriers (H > S)
    Real uoc = barrier_option_price(S, K, H_up, T, r, q, sigma, BarrierType::UpOutCall);
    Real uic = barrier_option_price(S, K, H_up, T, r, q, sigma, BarrierType::UpInCall);
    EXPECT_NEAR(uoc + uic, vc, 1e-10);

    Real uop = barrier_option_price(S, K, H_up, T, r, q, sigma, BarrierType::UpOutPut);
    Real uip = barrier_option_price(S, K, H_up, T, r, q, sigma, BarrierType::UpInPut);
    EXPECT_NEAR(uop + uip, vp, 1e-10);
}

// 2.2 Out 期权价格 ≤ Vanilla
TEST(Barrier, OutLeqVanilla) {
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real vc = bsm_call(S, K, T, r, q, sigma);
    Real vp = bsm_put(S, K, T, r, q, sigma);

    Real doc = barrier_option_price(S, K, 85.0, T, r, q, sigma, BarrierType::DownOutCall);
    Real dop = barrier_option_price(S, K, 85.0, T, r, q, sigma, BarrierType::DownOutPut);
    Real uoc = barrier_option_price(S, K, 115.0, T, r, q, sigma, BarrierType::UpOutCall);
    Real uop = barrier_option_price(S, K, 115.0, T, r, q, sigma, BarrierType::UpOutPut);

    EXPECT_LE(doc, vc + 1e-10);
    EXPECT_LE(dop, vp + 1e-10);
    EXPECT_LE(uoc, vc + 1e-10);
    EXPECT_LE(uop, vp + 1e-10);
}

// 2.3 In 期权价格 ≥ 0
TEST(Barrier, InNonNegative) {
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real dic = barrier_option_price(S, K, 85.0, T, r, q, sigma, BarrierType::DownInCall);
    Real dip = barrier_option_price(S, K, 85.0, T, r, q, sigma, BarrierType::DownInPut);
    Real uic = barrier_option_price(S, K, 115.0, T, r, q, sigma, BarrierType::UpInCall);
    Real uip = barrier_option_price(S, K, 115.0, T, r, q, sigma, BarrierType::UpInPut);
    EXPECT_GE(dic, 0.0);
    EXPECT_GE(dip, 0.0);
    EXPECT_GE(uic, 0.0);
    EXPECT_GE(uip, 0.0);
}

// 2.4 vs MC: Down-and-Out Call (H < K)
TEST(Barrier, DownOutCallVsMC) {
    Real S = 100.0, K = 100.0, H = 85.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real analytic = barrier_option_price(S, K, H, T, r, q, sigma, BarrierType::DownOutCall);
    Real mc = mc_barrier(S, K, H, T, r, q, sigma,
                          BarrierDirection::Down, BarrierKnock::Out, OptionType::Call);
    // MC 离散监控误差较大, 容差 0.5
    EXPECT_NEAR(analytic, mc, 0.50);
}

// 2.5 vs MC: Up-and-Out Call (K < H)
TEST(Barrier, UpOutCallVsMC) {
    Real S = 100.0, K = 100.0, H = 115.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real analytic = barrier_option_price(S, K, H, T, r, q, sigma, BarrierType::UpOutCall);
    Real mc = mc_barrier(S, K, H, T, r, q, sigma,
                          BarrierDirection::Up, BarrierKnock::Out, OptionType::Call);
    EXPECT_NEAR(analytic, mc, 0.50);
}

// 2.6 vs MC: Down-and-Out Put (H < K)
TEST(Barrier, DownOutPutVsMC) {
    Real S = 100.0, K = 110.0, H = 85.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.25;
    Real analytic = barrier_option_price(S, K, H, T, r, q, sigma, BarrierType::DownOutPut);
    Real mc = mc_barrier(S, K, H, T, r, q, sigma,
                          BarrierDirection::Down, BarrierKnock::Out, OptionType::Put);
    EXPECT_NEAR(analytic, mc, 0.50);
}

// 2.7 vs MC: Up-and-Out Put
TEST(Barrier, UpOutPutVsMC) {
    Real S = 100.0, K = 95.0, H = 115.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real analytic = barrier_option_price(S, K, H, T, r, q, sigma, BarrierType::UpOutPut);
    Real mc = mc_barrier(S, K, H, T, r, q, sigma,
                          BarrierDirection::Up, BarrierKnock::Out, OptionType::Put);
    EXPECT_NEAR(analytic, mc, 0.50);
}

// 2.8 vs MC: Down-and-In Call
TEST(Barrier, DownInCallVsMC) {
    Real S = 100.0, K = 100.0, H = 85.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real analytic = barrier_option_price(S, K, H, T, r, q, sigma, BarrierType::DownInCall);
    Real mc = mc_barrier(S, K, H, T, r, q, sigma,
                          BarrierDirection::Down, BarrierKnock::In, OptionType::Call);
    EXPECT_NEAR(analytic, mc, 0.50);
}

// 2.9 H → S: 障碍很远时 Out 期权 ≈ Vanilla (障碍几乎不触发)
TEST(Barrier, BarrierFarFromS) {
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real vc = bsm_call(S, K, T, r, q, sigma);
    Real vp = bsm_put(S, K, T, r, q, sigma);
    // H 非常低 → DOC/DOP ≈ Vanilla (down barrier 几乎不触发)
    Real doc_far = barrier_option_price(S, K, 10.0, T, r, q, sigma, BarrierType::DownOutCall);
    EXPECT_NEAR(doc_far, vc, 1e-6);
    // H 非常高 → UOC/UOP ≈ Vanilla (up barrier 几乎不触发)
    Real uoc_far = barrier_option_price(S, K, 1000.0, T, r, q, sigma, BarrierType::UpOutCall);
    EXPECT_NEAR(uoc_far, vc, 1e-6);
    Real uop_far = barrier_option_price(S, K, 1000.0, T, r, q, sigma, BarrierType::UpOutPut);
    EXPECT_NEAR(uop_far, vp, 1e-6);
}

// 2.10 H = K 边界情况 (两种公式分支在 H=K 处连续)
TEST(Barrier, HEqualsK) {
    // Down barrier: H < S 必须满足, 取 H = K = 95 < S = 100
    Real S = 100.0, K = 95.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real H = K;  // H = K = 95 < S = 100, 合法的 down barrier
    Real doc = barrier_option_price(S, K, H, T, r, q, sigma, BarrierType::DownOutCall);
    EXPECT_GT(doc, 0.0);
    EXPECT_LT(doc, bsm_call(S, K, T, r, q, sigma));
    Real dic = barrier_option_price(S, K, H, T, r, q, sigma, BarrierType::DownInCall);
    EXPECT_NEAR(doc + dic, bsm_call(S, K, T, r, q, sigma), 1e-10);

    // Up barrier: H > S 必须满足, 取 H = K = 105 > S = 100
    Real S2 = 100.0, K2 = 105.0;
    Real H2 = K2;  // H = K = 105 > S = 100, 合法的 up barrier
    Real uoc = barrier_option_price(S2, K2, H2, T, r, q, sigma, BarrierType::UpOutCall);
    // K = H: max < H ≤ K ⟹ S_T < K, call payoff = 0
    EXPECT_EQ(uoc, 0.0);
    Real uop = barrier_option_price(S2, K2, H2, T, r, q, sigma, BarrierType::UpOutPut);
    EXPECT_GT(uop, 0.0);
}

// 2.11 带股息率 q > 0
TEST(Barrier, WithDividend) {
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.03, sigma = 0.20;
    Real H = 85.0;
    Real vc = bsm_call(S, K, T, r, q, sigma);
    Real vp = bsm_put(S, K, T, r, q, sigma);
    Real doc = barrier_option_price(S, K, H, T, r, q, sigma, BarrierType::DownOutCall);
    Real dic = barrier_option_price(S, K, H, T, r, q, sigma, BarrierType::DownInCall);
    Real dop = barrier_option_price(S, K, H, T, r, q, sigma, BarrierType::DownOutPut);
    Real dip = barrier_option_price(S, K, H, T, r, q, sigma, BarrierType::DownInPut);
    EXPECT_NEAR(doc + dic, vc, 1e-10);
    EXPECT_NEAR(dop + dip, vp, 1e-10);
}

// ================================================================
// ============ 3. Lookback Options (Goldman-Sosin-Gatto) 测试 ============
// ================================================================

// 3.1 浮动 call vs MC (m = S, 刚启动)
TEST(Lookback, FloatingCallVsMC) {
    Real S = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real m = S;  // 刚启动, m = S
    Real analytic = lookback_call_floating(S, m, T, r, q, sigma);
    Real mc = mc_lookback_call_float(S, T, r, q, sigma, 500000);
    // Brownian Bridge 改进后 MC 标准误差约 0.3-0.5, 容差 1.5 容忍 BB+MC 噪声
    EXPECT_NEAR(analytic, mc, 1.5);
}

// 3.2 浮动 put vs MC (M = S, 刚启动)
TEST(Lookback, FloatingPutVsMC) {
    Real S = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real M = S;
    Real analytic = lookback_put_floating(S, M, T, r, q, sigma);
    Real mc = mc_lookback_put_float(S, T, r, q, sigma, 500000);
    EXPECT_NEAR(analytic, mc, 1.0);
}

// 3.3 浮动 call + 浮动 put = forward + lookback premium
// Call_float + Put_float = (S_T - m) + (M - S_T) = M - m (极差)
// E[M - m] 的折现 = Call + Put (无折现对应)
// 但这不是标准 parity, 跳过

// 3.4 浮动 call 价格 > vanilla call (lookback 有额外价值)
TEST(Lookback, FloatingCallGreaterThanVanilla) {
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real m = S;
    Real lb = lookback_call_floating(S, m, T, r, q, sigma);
    Real vc = bsm_call(S, K, T, r, q, sigma);
    EXPECT_GT(lb, vc);
}

// 3.5 浮动 put 价格 > vanilla put
TEST(Lookback, FloatingPutGreaterThanVanilla) {
    Real S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real M = S;
    Real lb = lookback_put_floating(S, M, T, r, q, sigma);
    Real vp = bsm_put(S, K, T, r, q, sigma);
    EXPECT_GT(lb, vp);
}

// 3.6 浮动 call: m < S (已运行最小值低于现货)
TEST(Lookback, FloatingCallWithRunningMin) {
    Real S = 100.0, T = 0.5, r = 0.05, q = 0.0, sigma = 0.20;
    Real m = 90.0;  // 已有最小值
    Real analytic = lookback_call_floating(S, m, T, r, q, sigma);
    // m < S 时, call = S_T - m_T, m_T ≤ m
    // 至少 payoff ≥ S_T - m (因为 m_T ≤ m)
    // 所以 call ≥ e^{-rT}(E[S_T] - m) = S*e^{-qT} - m*e^{-rT}
    Real lower_bound = S * std::exp(-q * T) - m * std::exp(-r * T);
    EXPECT_GT(analytic, lower_bound);
}

// 3.7 浮动 put: M > S (已运行最大值高于现货)
TEST(Lookback, FloatingPutWithRunningMax) {
    Real S = 100.0, T = 0.5, r = 0.05, q = 0.0, sigma = 0.20;
    Real M = 110.0;
    Real analytic = lookback_put_floating(S, M, T, r, q, sigma);
    // M_T ≥ M, 所以 payoff = M_T - S_T ≥ M - S_T
    // put ≥ e^{-rT}(M - E[S_T]) = M*e^{-rT} - S*e^{-qT}
    Real lower_bound = M * std::exp(-r * T) - S * std::exp(-q * T);
    EXPECT_GT(analytic, lower_bound);
}

// 3.8 b = 0 (r = q) 极限情况
TEST(Lookback, FloatingCallBZero) {
    Real S = 100.0, T = 1.0, r = 0.03, q = 0.03, sigma = 0.20;  // b = 0
    Real m = S;
    // 不应抛异常, 应给出有限值
    Real analytic = lookback_call_floating(S, m, T, r, q, sigma);
    EXPECT_TRUE(std::isfinite(analytic));
    EXPECT_GT(analytic, 0.0);
    Real mc = mc_lookback_call_float(S, T, r, q, sigma, 500000);
    EXPECT_NEAR(analytic, mc, 1.5);
}

// 3.9 固定 call: M > K (已 ITM)
TEST(Lookback, FixedCallITM) {
    Real S = 100.0, K = 95.0, M = 110.0, T = 0.5, r = 0.05, q = 0.0, sigma = 0.20;
    Real analytic = lookback_call_fixed(S, K, M, T, r, q, sigma);
    // payoff = M_T - K, M_T ≥ M > K
    // 至少 ≥ (M - K) * e^{-rT}
    Real lower = (M - K) * std::exp(-r * T);
    EXPECT_GT(analytic, lower);
}

// 3.10 固定 put: m < K (已 ITM)
TEST(Lookback, FixedPutITM) {
    Real S = 100.0, K = 105.0, m = 90.0, T = 0.5, r = 0.05, q = 0.0, sigma = 0.20;
    Real analytic = lookback_put_fixed(S, K, m, T, r, q, sigma);
    // payoff = K - m_T, m_T ≤ m < K
    // 至少 ≥ (K - m) * e^{-rT}
    Real lower = (K - m) * std::exp(-r * T);
    EXPECT_GT(analytic, lower);
}

// 3.11 固定 call: M ≤ K (刚启动, M = S ≤ K) vs MC
TEST(Lookback, FixedCallATMvsMC) {
    Real S = 100.0, K = 100.0, M = S, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real analytic = lookback_call_fixed(S, K, M, T, r, q, sigma);
    // MC: payoff = max(M_T - K, 0)
    auto cfg = make_single_asset_gbm(S, sigma, r, q, T, 200);
    MultiAssetGBMPathGenerator gen(cfg);
    Philox4x64 rng(42);
    Real df = std::exp(-r * T);
    Real sum = 0.0;
    Size n_paths = 500000;
    for (Size i = 0; i < n_paths; ++i) {
        auto paths = gen.generate_path(rng);
        const auto& path = paths[0];
        Real MT = *std::max_element(path.begin(), path.end());
        sum += std::max(MT - K, 0.0);
    }
    Real mc = df * sum / n_paths;
    EXPECT_NEAR(analytic, mc, 1.5);
}
