// v1.2: HJM (Heath-Jarrow-Morton) + LMM (LIBOR Market Model / BGM) 测试
// - HJM: 初始曲线复原, drift 解析, Markovian (HW 等价) 一致性, 模拟路径合理性
// - LMM: 初始债券复原, caplet Black 76 严格性, cap-floor parity, Rebonato swaption, 模拟路径
#include <gtest/gtest.h>
#include "cpphub/models/ir/hjm.hpp"
#include "cpphub/models/ir/lmm.hpp"
#include "cpphub/models/ir/short_rate.hpp"
#include "cpphub/core/rng.hpp"
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace cpphub;

namespace {
constexpr Real kTol = 1e-10;
constexpr Real kMCTol = 5e-3;  // MC 3σ 容差

// 平坦收益率曲线工厂 (用于 HJM/HW 对比)
std::pair<std::vector<Real>, std::vector<Real>> flat_curve(Real y, Real Tmax, Real dt) {
    std::vector<Real> mats, bonds;
    for (Real T = dt; T <= Tmax + 1e-10; T += dt) {
        mats.push_back(T);
        bonds.push_back(std::exp(-y * T));
    }
    return {mats, bonds};
}
}  // namespace

// ============================================================
// 1. HJM 测试
// ============================================================

TEST(HJMTest, InitialForwardCurveMatchesInput) {
    // 初始远期曲线在 tenor 节点处应等于输入
    auto cfg = make_flat_hjm(0.04, 10.0, 0.5,
                              [](Real /*tau*/) { return 0.01; });
    HJM hjm(cfg);
    for (Size i = 0; i < cfg.tenors.size(); ++i) {
        EXPECT_NEAR(hjm.forward_rate(cfg.tenors[i]), 0.04, kTol);
    }
}

TEST(HJMTest, ZeroCouponBondFlatCurve) {
    // 平坦远期曲线 f(0,T)=y 下, P(0,T) = exp(-y*T)
    auto cfg = make_flat_hjm(0.05, 10.0, 0.25,
                              [](Real /*tau*/) { return 0.01; });
    HJM hjm(cfg);
    for (Real T = 0.25; T <= 10.0; T += 0.5) {
        Real P = hjm.zero_coupon_bond(T);
        Real expected = std::exp(-0.05 * T);
        EXPECT_NEAR(P, expected, 1e-6) << "T=" << T;
    }
}

TEST(HJMTest, DriftConstantVol) {
    // σ(τ)=σ (常数) 时, drift μ(τ) = σ * ∫_0^τ σ du = σ² τ
    Real sigma = 0.01;
    auto cfg = make_flat_hjm(0.04, 5.0, 0.5,
                              [sigma](Real /*tau*/) { return sigma; });
    HJM hjm(cfg);
    for (Real tau = 0.5; tau <= 5.0; tau += 0.5) {
        Real mu = hjm.drift(tau);
        Real expected = sigma * sigma * tau;
        EXPECT_NEAR(mu, expected, 1e-8) << "tau=" << tau;
    }
}

TEST(HJMTest, DriftExponentialVolHW) {
    // σ(τ) = σ_hw exp(-κτ) 时, drift μ(τ) = σ_hw exp(-κτ) * (σ_hw/κ)(1 - exp(-κτ))
    //                                     = (σ_hw²/κ) exp(-κτ) (1 - exp(-κτ))
    Real sigma_hw = 0.01, kappa = 0.1;
    auto cfg = make_hw_equivalent_hjm(0.04, kappa, sigma_hw, 10.0, 0.25);
    HJM hjm(cfg);
    for (Real tau = 0.5; tau <= 10.0; tau += 0.5) {
        Real mu = hjm.drift(tau);
        Real expected = (sigma_hw * sigma_hw / kappa) * std::exp(-kappa * tau)
                      * (1.0 - std::exp(-kappa * tau));
        EXPECT_NEAR(mu, expected, 1e-8) << "tau=" << tau;
    }
}

TEST(HJMTest, RejectsInvalidConfig) {
    HJMConfig cfg;
    cfg.tenors = {0.0, 1.0};
    cfg.initial_forwards = {0.04};  // size mismatch
    cfg.vol_functions = {[](Real) { return 0.01; }};
    EXPECT_THROW(HJM hjm(cfg), std::invalid_argument);

    cfg.tenors = {1.0, 2.0};  // tenors[0] != 0
    cfg.initial_forwards = {0.04, 0.04};
    EXPECT_THROW(HJM hjm(cfg), std::invalid_argument);

    cfg.tenors = {0.0, 1.0};
    cfg.vol_functions = {};  // no factor
    EXPECT_THROW(HJM hjm(cfg), std::invalid_argument);
}

TEST(HJMTest, SimulatePathPreservesForwardShape) {
    // 零波动率下, 模拟路径上 forward 应保持初始曲线 (无扩散, drift 也不变累积)
    // 但 drift 仍存在; 此处仅验证零波动率时 forward 增量 = drift * dt
    auto cfg = make_flat_hjm(0.04, 5.0, 0.5,
                              [](Real) { return 0.0; });  // 零波动率
    HJM hjm(cfg);
    Philox4x64 rng(42, 0);
    hjm.simulate_path(2.0, 100, rng);
    // 零波动率 + 平坦曲线: forward 演化为 f(0,T) + drift*dt 累积, drift=0 时严格不变
    for (Size i = 0; i < cfg.tenors.size(); ++i) {
        if (cfg.tenors[i] < 2.0) continue;  // 仅检查未到期 tenor
        Real f_end = hjm.forward_rate_at_step(100, cfg.tenors[i]);
        EXPECT_NEAR(f_end, 0.04, 1e-10) << "i=" << i << " T=" << cfg.tenors[i];
    }
}

TEST(HJMTest, SimulatedShortRateMatchesInitial) {
    // 零波动率 + 平坦曲线下, r(t) = f(t, t) ≈ 0.04 (因 drift 严格为 0)
    auto cfg = make_flat_hjm(0.04, 5.0, 0.25,
                              [](Real) { return 0.0; });
    HJM hjm(cfg);
    Philox4x64 rng(42, 0);
    hjm.simulate_path(2.0, 200, rng);
    Real r = hjm.short_rate_at(200);
    EXPECT_NEAR(r, 0.04, 1e-6);
}

TEST(HJMTest, HWEquivalentDriftMatches) {
    // HJM (σ=σ_hw*exp(-κτ)) 的 drift 与 HW bond option 隐含的 drift 一致
    // HW bond option 已在 short_rate.hpp 中验证, 此处仅验证 HJM drift 公式正确
    Real sigma_hw = 0.01, kappa = 0.1;
    auto cfg = make_hw_equivalent_hjm(0.04, kappa, sigma_hw, 10.0, 0.25);
    HJM hjm(cfg);
    // μ(τ) = (σ²/κ) exp(-κτ) (1 - exp(-κτ))
    Real tau = 5.0;
    Real mu = hjm.drift(tau);
    Real expected = (sigma_hw * sigma_hw / kappa) * std::exp(-kappa * tau)
                  * (1.0 - std::exp(-kappa * tau));
    // 梯形积分 64 步, 误差 ~1e-9
    EXPECT_NEAR(mu, expected, 1e-8);
}

TEST(HJMTest, MultiFactorRejectsBadCorrelation) {
    HJMConfig cfg;
    cfg.tenors = {0.0, 1.0, 2.0};
    cfg.initial_forwards = {0.04, 0.04, 0.04};
    cfg.vol_functions = {[](Real) { return 0.01; }, [](Real) { return 0.005; }};
    // 未提供 correlation → 应抛
    EXPECT_THROW(HJM hjm(cfg), std::invalid_argument);

    // 非单位对角线
    cfg.correlation = {{0.5, 0.3}, {0.3, 0.5}};
    EXPECT_THROW(HJM hjm(cfg), std::invalid_argument);
}

TEST(HJMTest, MultiFactorSimulateRuns) {
    // 双因子 HJM 模拟应能运行 (相关性矩阵正定)
    HJMConfig cfg;
    cfg.tenors = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
    cfg.initial_forwards = {0.04, 0.04, 0.042, 0.043, 0.044, 0.045, 0.046};
    cfg.vol_functions = {
        [](Real tau) { return 0.01 * std::exp(-0.1 * tau); },
        [](Real tau) { return 0.005 * std::exp(-0.2 * tau); }
    };
    cfg.correlation = {{1.0, 0.5}, {0.5, 1.0}};
    HJM hjm(cfg);
    Philox4x64 rng(123, 0);
    hjm.simulate_path(1.0, 50, rng);
    // 简单检查: 模拟后 forward 在合理范围内 (非负, 非 NaN)
    auto f = hjm.current_forward_state();
    for (Real x : f) {
        EXPECT_TRUE(std::isfinite(x));
        EXPECT_GT(x, -0.5) << "forward rate should not explode";
        EXPECT_LT(x, 0.5);
    }
}

// ============================================================
// 2. LMM 测试
// ============================================================

TEST(LMMTest, InitialBondsFromFlatLibors) {
    // 平坦 LIBOR L 下, P(0, T_i) = 1/(1+L*τ)^i
    Real L = 0.05, tau = 0.5, sigma = 0.2;
    auto cfg = make_flat_lmm(L, 5.0, tau, sigma);
    LMM lmm(cfg);
    Size N = cfg.n_libors();
    for (Size i = 0; i <= N; ++i) {
        Real T_i = cfg.tenors[i];
        Real P = lmm.zero_coupon_bond(T_i);
        Real expected = std::pow(1.0 / (1.0 + L * tau), static_cast<Real>(i));
        EXPECT_NEAR(P, expected, 1e-10) << "i=" << i;
    }
}

TEST(LMMTest, InitialLiborFromBonds) {
    // L(0; T_a, T_b) = (P(0,T_a)/P(0,T_b) - 1) / (T_b - T_a)
    Real L = 0.05, tau = 0.5, sigma = 0.2;
    auto cfg = make_flat_lmm(L, 5.0, tau, sigma);
    LMM lmm(cfg);
    for (Size i = 0; i < cfg.n_libors(); ++i) {
        Real L_i = lmm.initial_libor(static_cast<Real>(cfg.tenors[i]),
                                       static_cast<Real>(cfg.tenors[i + 1]));
        EXPECT_NEAR(L_i, L, 1e-10) << "i=" << i;
    }
}

TEST(LMMTest, CapletBlack76Formula) {
    // LMM caplet 应等于 Black 76 公式 (在 L_i 自身 forward measure 下严格成立)
    Real L = 0.05, tau = 0.5, sigma = 0.2;
    auto cfg = make_flat_lmm(L, 5.0, tau, sigma);
    LMM lmm(cfg);

    Real K = 0.045;
    Size i = 2;  // L_2, reset at T_2=1.0
    Real price = lmm.caplet_price(i, K, 100.0);

    // 手算 Black 76
    Real F = L;
    Real T_i = cfg.tenors[i];
    Real P = lmm.zero_coupon_bond(cfg.tenors[i + 1]);
    Real d1 = (std::log(F / K) + 0.5 * sigma * sigma * T_i) / (sigma * std::sqrt(T_i));
    Real d2 = d1 - sigma * std::sqrt(T_i);
    Real expected = tau * 100.0 * P * (F * normal_cdf(d1) - K * normal_cdf(d2));
    EXPECT_NEAR(price, expected, 1e-10);
}

TEST(LMMTest, CapletFloorletParityAtParRate) {
    // caplet_i + K*τ*P = floorlet_i + F*τ*P (forward synthetic)
    // 等价: caplet - floorlet = (F - K) * τ * P * N
    Real L = 0.05, tau = 0.5, sigma = 0.2;
    auto cfg = make_flat_lmm(L, 5.0, tau, sigma);
    LMM lmm(cfg);
    Real K = 0.045;
    Size i = 1;
    Real caplet = lmm.caplet_price(i, K, 1.0);
    Real floorlet = lmm.floorlet_price(i, K, 1.0);
    Real F = L;
    Real P = lmm.zero_coupon_bond(cfg.tenors[i + 1]);
    Real expected_diff = (F - K) * tau * P;
    EXPECT_NEAR(caplet - floorlet, expected_diff, 1e-10);
}

TEST(LMMTest, CapFloorParity) {
    // Cap - Floor = Payer Forward Swap
    Real L = 0.05, tau = 0.5, sigma = 0.2;
    auto cfg = make_flat_lmm(L, 5.0, tau, sigma);
    LMM lmm(cfg);
    Real K = 0.05;
    Real cap_pv = lmm.cap_price(K, 100.0, 0);
    Real floor_pv = lmm.floor_price(K, 100.0, 0);
    Real payer_swap = lmm.cap_floor_parity_payer_swap(K, 100.0, 0);
    EXPECT_NEAR(cap_pv - floor_pv, payer_swap, 1e-10);
}

TEST(LMMTest, CapIncreasesWithVol) {
    Real L = 0.05, tau = 0.5;
    auto cfg_low = make_flat_lmm(L, 5.0, tau, 0.1);
    auto cfg_high = make_flat_lmm(L, 5.0, tau, 0.4);
    LMM lmm_low(cfg_low), lmm_high(cfg_high);
    Real K = 0.055;
    EXPECT_LT(lmm_low.cap_price(K, 100.0), lmm_high.cap_price(K, 100.0));
}

TEST(LMMTest, ZeroVolCapletIntrinsic) {
    // σ=0 时, caplet = 内在价值 (折现)
    Real L = 0.05, tau = 0.5;
    auto cfg = make_flat_lmm(L, 5.0, tau, 0.0);
    LMM lmm(cfg);
    Real K = 0.045;
    Size i = 2;
    Real price = lmm.caplet_price(i, K, 100.0);
    Real F = L;
    Real P = lmm.zero_coupon_bond(cfg.tenors[i + 1]);
    Real intrinsic = tau * 100.0 * P * std::max(F - K, 0.0);
    EXPECT_NEAR(price, intrinsic, 1e-10);
}

TEST(LMMTest, RejectsInvalidConfig) {
    LMMConfig cfg;
    cfg.tenors = {0.0, 1.0};
    cfg.initial_libors = {0.05, 0.05};  // size mismatch (N=2 but tenors N+1=2)
    cfg.volatilities = {0.2, 0.2};
    cfg.correlation = {{1.0, 0.5}, {0.5, 1.0}};
    EXPECT_THROW(LMM lmm(cfg), std::invalid_argument);

    cfg.tenors = {0.0, 1.0, 2.0};
    cfg.initial_libors = {0.05, 0.05};
    cfg.volatilities = {0.2, 0.2};
    cfg.correlation = {{1.0, 0.5}, {0.5, 1.0}};
    cfg.measure = "bad";
    EXPECT_THROW(LMM lmm(cfg), std::invalid_argument);

    cfg.measure = "spot";
    cfg.correlation = {{1.0, 1.2}, {1.2, 1.0}};  // |ρ|>1
    EXPECT_THROW(LMM lmm(cfg), std::invalid_argument);
}

TEST(LMMTest, SwaptionRebonatoPayerReceiverSymmetry) {
    // ATM swaption (K = S(0)) 时, payer + receiver ≈ 2 * ATM price
    // Payer - Receiver = S - K (forward swap value)
    Real L = 0.05, tau = 0.5, sigma = 0.2;
    auto cfg = make_flat_lmm(L, 5.0, tau, sigma);
    LMM lmm(cfg);

    // T_ex 必须 > 0, 否则 sigma_swap = sqrt(var * 0) = 0, swaption 退化为内在价值
    Size a = 2, b = cfg.n_libors() - 1;  // T_ex = T_2 = 1.0, swap [T_2, T_N]
    Real T_ex = cfg.tenors[a];
    EXPECT_GT(T_ex, 0.0);
    // par swap rate S = (P_a - P_b)/A
    Real P_a = lmm.zero_coupon_bond(cfg.tenors[a]);
    Real P_b = lmm.zero_coupon_bond(cfg.tenors[b + 1]);
    Real A = 0.0;
    for (Size i = a; i <= b; ++i) A += tau * lmm.zero_coupon_bond(cfg.tenors[i + 1]);
    Real S = (P_a - P_b) / A;

    Real payer = lmm.swaption_price_rebonato(T_ex, S, true, a, b, 1.0);
    Real receiver = lmm.swaption_price_rebonato(T_ex, S, false, a, b, 1.0);
    // ATM 时 payer ≈ receiver (高斯偏差微小)
    Real avg = 0.5 * (payer + receiver);
    EXPECT_NEAR(payer, receiver, 0.05 * avg + 1e-10);

    // payer - receiver = S - K = 0 (ATM)
    EXPECT_NEAR(payer - receiver, 0.0, 1e-10);

    // OTM payer (K > S) 应小于 ATM
    Real payer_otm = lmm.swaption_price_rebonato(T_ex, S + 0.01, true, a, b, 1.0);
    EXPECT_LT(payer_otm, payer);
}

TEST(LMMTest, SimulatePathKeepsLiborsPositive) {
    // Euler on log L 保证 L > 0
    Real L = 0.05, tau = 0.5, sigma = 0.3;
    auto cfg = make_flat_lmm(L, 5.0, tau, sigma);
    LMM lmm(cfg);
    Philox4x64 rng(42, 0);
    lmm.simulate_path(2.0, 200, rng);
    // 检查每一步每个 LIBOR 为正
    for (Size step = 0; step <= 200; ++step) {
        for (Size i = 0; i < cfg.n_libors(); ++i) {
            Real L_i = lmm.libor_at_step(step, i);
            EXPECT_GT(L_i, 0.0) << "step=" << step << " i=" << i;
            EXPECT_TRUE(std::isfinite(L_i));
        }
    }
}

TEST(LMMTest, SimulateZeroVolPreservesLibors) {
    // σ=0 时, L_i(t) = L_i(0) (drift 也为 0, 因为 drift 含 σ_j)
    Real L = 0.05, tau = 0.5;
    auto cfg = make_flat_lmm(L, 5.0, tau, 0.0);
    LMM lmm(cfg);
    Philox4x64 rng(42, 0);
    lmm.simulate_path(2.0, 200, rng);
    for (Size step = 0; step <= 200; ++step) {
        for (Size i = 0; i < cfg.n_libors(); ++i) {
            EXPECT_NEAR(lmm.libor_at_step(step, i), L, 1e-12)
                << "step=" << step << " i=" << i;
        }
    }
}

TEST(LMMTest, SimulatedBondMatchesInitialAtZeroVol) {
    // σ=0 时, P(t, T) 沿模拟路径应保持等于 P(0, T)/P(0, t) (无随机演化)
    // 简化验证: t=0 时 P(0, T) 严格复原
    Real L = 0.05, tau = 0.5;
    auto cfg = make_flat_lmm(L, 5.0, tau, 0.0);
    LMM lmm(cfg);
    // 模拟前 P(0, T) = 初始债券
    for (Size i = 0; i <= cfg.n_libors(); ++i) {
        Real P = lmm.zero_coupon_bond(cfg.tenors[i]);
        Real expected = std::pow(1.0 / (1.0 + L * tau), static_cast<Real>(i));
        EXPECT_NEAR(P, expected, 1e-10);
    }
}

// ============================================================
// 3. HJM-LMM 跨模型一致性
// ============================================================

TEST(HJMLMMConsistencyTest, FlatCurveZeroCouponBondsMatch) {
    // 平坦曲线下, HJM (f=0.05) 和 LMM (L=0.05) 给出的 P(0, T) 应一致
    // (连续复利 vs 简单复利, 差异 ~ y²τ/2)
    Real y = 0.05, tau = 0.5, sigma = 0.01;
    auto hjm_cfg = make_flat_hjm(y, 5.0, 0.25, [sigma](Real) { return sigma; });
    auto lmm_cfg = make_flat_lmm(y, 5.0, tau, sigma);
    HJM hjm(hjm_cfg);
    LMM lmm(lmm_cfg);

    // 在各自 tenor 节点上比较
    // HJM P(0, T) = exp(-y*T)
    // LMM P(0, T_i) = 1/(1+y*τ)^i  (简单复利)
    // 差异 = convexity, 应在 y²τ/2 量级
    for (Real T = 0.5; T <= 5.0; T += 0.5) {
        Real P_hjm = hjm.zero_coupon_bond(T);
        Real P_lmm = lmm.zero_coupon_bond(T);
        // 相对差异应小于 1%
        Real rel_diff = std::abs(P_hjm - P_lmm) / P_hjm;
        EXPECT_LT(rel_diff, 0.01) << "T=" << T;
    }
}

TEST(HJMLMMConsistencyTest, HWEquivalentHJMMatchesHWBonds) {
    // HJM (σ_hw exp(-κτ)) 给出的 P(0, T) 应等于 HW 给出的 P(0, T)
    // (因 HW 是 HJM 的 Markovian 特例, 两者校准到同一初始曲线)
    Real r0 = 0.04, kappa = 0.1, sigma_hw = 0.01;
    auto [mats, bonds] = flat_curve(r0, 10.0, 0.25);
    HullWhiteParams hw_p{r0, kappa, sigma_hw};
    HullWhite hw(hw_p, mats, bonds);

    auto hjm_cfg = make_hw_equivalent_hjm(r0, kappa, sigma_hw, 10.0, 0.25);
    HJM hjm(hjm_cfg);

    for (Real T = 0.5; T <= 10.0; T += 0.5) {
        Real P_hw = hw.zero_coupon_bond(T);
        Real P_hjm = hjm.zero_coupon_bond(T);
        // 两者都是 exp(-r0*T) (平坦曲线)
        EXPECT_NEAR(P_hjm, P_hw, 1e-6) << "T=" << T;
        EXPECT_NEAR(P_hjm, std::exp(-r0 * T), 1e-6);
    }
}
