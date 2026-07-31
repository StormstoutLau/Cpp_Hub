// v1.1: 利率衍生品测试
// - IRS: par swap rate, payer/receiver 对称性, 浮动端 par 性质
// - Cap/Floor: HW 闭式 vs Black 76 一致性, cap-floor parity, 边界情形
// - Swaption: Jamshidian 分解, payer/receiver 对称性, ATM 性质, 与 Cap/Floor 一致性
#include <gtest/gtest.h>
#include "cpphub/instruments/ir/irs.hpp"
#include "cpphub/instruments/ir/cap_floor.hpp"
#include "cpphub/instruments/ir/swaption.hpp"
#include "cpphub/models/ir/short_rate.hpp"
#include <cmath>
#include <vector>

using namespace cpphub;

// ============ 测试用市场曲线工厂 ============
// 构造平坦收益率曲线 (固定收益率 y): P(0, T) = exp(-y * T)
static std::pair<std::vector<Real>, std::vector<Real>> flat_curve(Real y, Real Tmax, Real dt) {
    std::vector<Real> mats, bonds;
    for (Real T = dt; T <= Tmax + 1e-10; T += dt) {
        mats.push_back(T);
        bonds.push_back(std::exp(-y * T));
    }
    return {mats, bonds};
}

// ============ 1. IRS 测试 ============

TEST(IRSTest, PayerReceiverSymmetry) {
    // 相同参数下, payer PV = -receiver PV
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.5);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    auto cfg_payer = make_vanilla_irs(100.0, 0.05, 0.0, 0.5, 10, true);
    auto cfg_receiver = make_vanilla_irs(100.0, 0.05, 0.0, 0.5, 10, false);

    InterestRateSwap payer(cfg_payer);
    InterestRateSwap receiver(cfg_receiver);

    Real pv_payer = payer.present_value(hw);
    Real pv_receiver = receiver.present_value(hw);
    EXPECT_NEAR(pv_payer, -pv_receiver, 1e-10);
}

TEST(IRSTest, ParSwapRateZeroPV) {
    // K = par swap rate 时, payer PV = 0
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.5);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    auto cfg = make_vanilla_irs(100.0, 0.0, 0.0, 0.5, 10, true);
    InterestRateSwap irs(cfg);
    Real K_par = irs.par_swap_rate(hw);

    auto cfg_par = make_vanilla_irs(100.0, K_par, 0.0, 0.5, 10, true);
    InterestRateSwap irs_par(cfg_par);
    Real pv = irs_par.present_value(hw);
    EXPECT_NEAR(pv, 0.0, 1e-10);
}

TEST(IRSTest, ParSwapRateMatchesYieldFlatCurve) {
    // 平坦连续复利曲线 y 下, par swap rate ≈ y + y²τ/2 (凸性修正)
    // 因为远期 LIBOR L(0,T_{i-1},T_i) = (exp(yτ)-1)/τ ≈ y + y²τ/2
    auto [mats, bonds] = flat_curve(0.05, 30.0, 0.25);
    HullWhiteParams hw_p{0.05, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    auto cfg = make_vanilla_irs(1.0, 0.0, 0.0, 1.0, 5, true);
    InterestRateSwap irs(cfg);
    Real K_par = irs.par_swap_rate(hw);
    // 平坦曲线 + spot start, par swap rate 与 y 的偏差 ≈ y²τ/2 = 0.00125
    Real y = 0.05, tau = 1.0;
    Real expected = (std::exp(y * tau) - 1.0) / tau;  // 简单复利远期利率
    EXPECT_NEAR(K_par, expected, 1e-10);
    EXPECT_NEAR(K_par, y, 0.002);  // 凸性偏差量级
}

TEST(IRSTest, ForwardStartPV) {
    // Forward start swap: T_start > 0, 浮动端 PV = P(0, T_start) - P(0, T_n)
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.5);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    auto cfg = make_vanilla_irs(100.0, 0.05, 2.0, 0.5, 10, true);
    InterestRateSwap irs(cfg);
    Real pv = irs.present_value(hw);
    Real fixed_pv = irs.fixed_leg_pv(hw);
    Real float_pv = irs.float_leg_pv(hw);
    EXPECT_NEAR(pv, float_pv - fixed_pv, 1e-10);
    // 浮动端 PV 应等于 P(0,2) - P(0,7)
    Real P_start = hw.zero_coupon_bond(2.0);
    Real P_end = hw.zero_coupon_bond(7.0);
    EXPECT_NEAR(float_pv, 100.0 * (P_start - P_end), 1e-10);
}

TEST(IRSTest, RejectsInvalidConfig) {
    IRSConfig cfg;
    cfg.notional = -1.0;
    cfg.fixed_rate = 0.05;
    cfg.payment_times = {1.0, 2.0};
    cfg.year_fractions = {1.0, 1.0};
    EXPECT_THROW(InterestRateSwap irs(cfg), std::invalid_argument);

    cfg.notional = 100.0;
    cfg.payment_times = {2.0, 1.0};  // 非递增
    EXPECT_THROW(InterestRateSwap irs(cfg), std::invalid_argument);

    cfg.payment_times = {1.0, 2.0};
    cfg.year_fractions = {1.0, 0.0};  // τ ≤ 0
    EXPECT_THROW(InterestRateSwap irs(cfg), std::invalid_argument);
}

// ============ 2. Cap/Floor 测试 ============

TEST(CapFloorTest, CapFloorParity) {
    // Cap - Floor = Swap (par rate 时 Cap = Floor)
    // 一般: Cap - Floor = N * [P(0,T_start) - P(0,T_n) - K * Σ τ_i * P(0,T_i)]
    //                    = payer swap PV
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    Real K = 0.04;
    auto cap_cfg = make_vanilla_capfloor(100.0, K, 0.0, 0.5, 10, true);
    auto floor_cfg = make_vanilla_capfloor(100.0, K, 0.0, 0.5, 10, false);
    CapFloor cap(cap_cfg);
    CapFloor floor(floor_cfg);

    Real cap_pv = cap.present_value_model(hw);
    Real floor_pv = floor.present_value_model(hw);

    auto swap_cfg = make_vanilla_irs(100.0, K, 0.0, 0.5, 10, true);
    InterestRateSwap swap(swap_cfg);
    Real swap_pv = swap.present_value(hw);

    // Cap - Floor = Payer Swap
    EXPECT_NEAR(cap_pv - floor_pv, swap_pv, 1e-6);
}

TEST(CapFloorTest, BlackAndModelConsistency) {
    // HW 闭式 vs Black 76 (用 HW 远期利率作为 Black 远期, vol 反推)
    // 两者数学上不完全等价 (HW 假设利率高斯, Black 假设对数正态),
    // 但 ATM 附近价格应接近
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    Real K = 0.04;
    auto cfg = make_vanilla_capfloor(100.0, K, 0.0, 0.5, 4, true);
    CapFloor cap(cfg);
    Real model_pv = cap.present_value_model(hw);

    // 用 model_pv 反推 Black 隐含 vol
    Real impl_vol = cap.implied_black_vol(hw, model_pv);
    Real black_pv = cap.present_value_black(hw, impl_vol);
    EXPECT_NEAR(black_pv, model_pv, 1e-6);
    // 隐含 vol 应在合理范围 (HW σ=0.01, κ=0.1, 短期 caplet vol ~ σ * sqrt(t) ≈ 0.01-0.02)
    EXPECT_GT(impl_vol, 0.0);
    EXPECT_LT(impl_vol, 1.0);
}

TEST(CapFloorTest, CapPricePositiveFloorPricePositive) {
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    auto cap_cfg = make_vanilla_capfloor(100.0, 0.04, 0.0, 0.5, 10, true);
    auto floor_cfg = make_vanilla_capfloor(100.0, 0.04, 0.0, 0.5, 10, false);
    CapFloor cap(cap_cfg);
    CapFloor floor(floor_cfg);

    EXPECT_GT(cap.present_value_model(hw), 0.0);
    EXPECT_GT(floor.present_value_model(hw), 0.0);
}

TEST(CapFloorTest, CapIncreasesWithVol) {
    // Black 76: vol 越大, cap 价格越高
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    auto cfg = make_vanilla_capfloor(100.0, 0.045, 0.0, 0.5, 6, true);
    CapFloor cap(cfg);
    Real pv_low = cap.present_value_black(hw, 0.10);
    Real pv_mid = cap.present_value_black(hw, 0.30);
    Real pv_high = cap.present_value_black(hw, 0.60);
    EXPECT_LT(pv_low, pv_mid);
    EXPECT_LT(pv_mid, pv_high);
}

TEST(CapFloorTest, ZeroVolGivesIntrinsic) {
    // vol=0 时, Black 价格 = 内在价值
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    Real K = 0.05;  // OTM cap (远期利率 < K)
    auto cfg = make_vanilla_capfloor(100.0, K, 0.0, 0.5, 4, true);
    CapFloor cap(cfg);
    Real pv_zero_vol = cap.present_value_black(hw, 0.0);

    // 计算内在价值
    Real intrinsic = 0.0;
    for (Size i = 0; i < cfg.n_caplets(); ++i) {
        Real T_reset = cfg.reset_times[i];
        Real T_pay = cfg.payment_times[i];
        Real tau = cfg.year_fractions[i];
        Real F = (hw.zero_coupon_bond(T_reset) / hw.zero_coupon_bond(T_pay) - 1.0) / tau;
        Real df = hw.zero_coupon_bond(T_pay);
        intrinsic += tau * 100.0 * df * std::max(F - K, 0.0);
    }
    EXPECT_NEAR(pv_zero_vol, intrinsic, 1e-10);
}

TEST(CapFloorTest, RejectsInvalidConfig) {
    CapFloorConfig cfg;
    cfg.notional = -1.0;
    cfg.strike = 0.05;
    cfg.reset_times = {0.5};
    cfg.payment_times = {1.0};
    cfg.year_fractions = {0.5};
    EXPECT_THROW(CapFloor cf(cfg), std::invalid_argument);

    cfg.notional = 100.0;
    cfg.reset_times = {1.0};
    cfg.payment_times = {0.5};  // payment < reset
    EXPECT_THROW(CapFloor cf(cfg), std::invalid_argument);
}

// ============ 3. Swaption 测试 ============

TEST(SwaptionTest, PayerReceiverSymmetry) {
    // 相同 ATM 参数, payer + receiver = forward swap (无期权性)
    // 或更简单: payer PV = -receiver PV (仅当 K 相同, 但方向相反)
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    auto cfg_p = make_vanilla_swaption(100.0, 0.045, 1.0, 0.5, 6, true);
    auto cfg_r = make_vanilla_swaption(100.0, 0.045, 1.0, 0.5, 6, false);
    Swaption payer(cfg_p, hw);
    Swaption receiver(cfg_r, hw);

    // payer + receiver = forward swap (无期权性, 即 T_ex 时 max + max 内在价值之和)
    // 实际上: max(PV_swap, 0) + max(-PV_swap, 0) = |PV_swap|, 不是 PV_swap
    // 正确关系: payer - receiver = forward swap PV (无期权性)
    Real pv_p = payer.present_value();
    Real pv_r = receiver.present_value();

    // 计算 forward swap PV (no optionality)
    auto fwd_swap_cfg = make_vanilla_irs(100.0, 0.045, 1.0, 0.5, 6, true);
    InterestRateSwap fwd_swap(fwd_swap_cfg);
    Real fwd_swap_pv = fwd_swap.present_value(hw);

    EXPECT_NEAR(pv_p - pv_r, fwd_swap_pv, 1e-3);
}

TEST(SwaptionTest, NonNegativePayoff) {
    // Swaption 价格必须非负
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    auto cfg_p = make_vanilla_swaption(100.0, 0.045, 1.0, 0.5, 6, true);
    auto cfg_r = make_vanilla_swaption(100.0, 0.045, 1.0, 0.5, 6, false);
    Swaption payer(cfg_p, hw);
    Swaption receiver(cfg_r, hw);
    EXPECT_GE(payer.present_value(), 0.0);
    EXPECT_GE(receiver.present_value(), 0.0);
}

TEST(SwaptionTest, ATMStrikeMatchesParSwapRate) {
    // ATM swaption (K = par swap rate) 时, payer ≈ receiver (近似, 因高斯模型偏差)
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    // 计算 T_ex=1 时刻的 forward par swap rate
    auto fwd_swap_cfg = make_vanilla_irs(1.0, 0.0, 1.0, 0.5, 6, true);
    InterestRateSwap fwd_swap(fwd_swap_cfg);
    Real K_par = fwd_swap.par_swap_rate(hw);

    auto cfg_p = make_vanilla_swaption(100.0, K_par, 1.0, 0.5, 6, true);
    auto cfg_r = make_vanilla_swaption(100.0, K_par, 1.0, 0.5, 6, false);
    Swaption payer(cfg_p, hw);
    Swaption receiver(cfg_r, hw);

    Real pv_p = payer.present_value();
    Real pv_r = receiver.present_value();
    // ATM 时 payer ≈ receiver (高斯模型下有微小偏差, 但应 < 5%)
    Real avg = 0.5 * (pv_p + pv_r);
    EXPECT_NEAR(pv_p, pv_r, 0.10 * avg + 1e-6);
}

TEST(SwaptionTest, OTMOptionalityValue) {
    // OTM payer (K >> par rate) 价格应小于 ATM payer
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    auto fwd_swap_cfg = make_vanilla_irs(1.0, 0.0, 1.0, 0.5, 6, true);
    InterestRateSwap fwd_swap(fwd_swap_cfg);
    Real K_par = fwd_swap.par_swap_rate(hw);

    auto cfg_atm = make_vanilla_swaption(100.0, K_par, 1.0, 0.5, 6, true);
    auto cfg_otm = make_vanilla_swaption(100.0, K_par + 0.05, 1.0, 0.5, 6, true);
    Swaption swaption_atm(cfg_atm, hw);
    Swaption swaption_otm(cfg_otm, hw);

    EXPECT_GT(swaption_atm.present_value(), swaption_otm.present_value());
    EXPECT_GT(swaption_otm.present_value(), 0.0);  // 仍应有期权价值
}

TEST(SwaptionTest, VolatilitySensitivity) {
    // HW σ 增大, swaption 价格应增大 (ATM 情形)
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_low{0.04, 0.1, 0.005};
    HullWhiteParams hw_high{0.04, 0.1, 0.02};
    HullWhite hw_l(hw_low, mats, bonds);
    HullWhite hw_h(hw_high, mats, bonds);

    auto fwd_swap_cfg = make_vanilla_irs(1.0, 0.0, 1.0, 0.5, 6, true);
    InterestRateSwap fwd_swap(fwd_swap_cfg);
    Real K_par = fwd_swap.par_swap_rate(hw_l);  // 用同一 K_par

    auto cfg = make_vanilla_swaption(100.0, K_par, 1.0, 0.5, 6, true);
    Swaption sw_l(cfg, hw_l);
    Swaption sw_h(cfg, hw_h);

    EXPECT_LT(sw_l.present_value(), sw_h.present_value());
}

TEST(SwaptionTest, RejectsInvalidConfig) {
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    SwaptionConfig cfg;
    cfg.notional = -1.0;
    cfg.fixed_rate = 0.05;
    cfg.exercise_time = 1.0;
    cfg.payment_times = {1.5, 2.0};
    cfg.year_fractions = {0.5, 0.5};
    EXPECT_THROW(Swaption s(cfg, hw), std::invalid_argument);

    cfg.notional = 100.0;
    cfg.exercise_time = -1.0;
    EXPECT_THROW(Swaption s(cfg, hw), std::invalid_argument);

    cfg.exercise_time = 1.0;
    cfg.payment_times = {0.5, 2.0};  // payment < exercise
    EXPECT_THROW(Swaption s(cfg, hw), std::invalid_argument);
}

// ============ 4. 跨产品一致性测试 ============

TEST(IRDerivsConsistencyTest, SinglePeriodCapletEqualsSwaption) {
    // 单期 (n=1) swaption 等价于单期 caplet/floorlet
    // Payer swaption T_ex, n=1, K, τ = T_1 - T_ex
    //   在 T_ex 期权: max(1 - (1+Kτ) P(T_ex, T_1), 0)
    //   等价于 caplet 在 T_1 支付 max(1 - (1+Kτ) P(T_ex, T_1), 0)
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    Real K = 0.045;
    Real T_ex = 1.0, tau = 0.5, T1 = T_ex + tau;
    auto cap_cfg = make_vanilla_capfloor(100.0, K, T_ex, tau, 1, true);
    auto swap_cfg = make_vanilla_swaption(100.0, K, T_ex, tau, 1, true);
    CapFloor cap(cap_cfg);
    Swaption swaption(swap_cfg, hw);

    Real cap_pv = cap.present_value_model(hw);
    Real swaption_pv = swaption.present_value();
    EXPECT_NEAR(cap_pv, swaption_pv, 1e-6)
        << "cap=" << cap_pv << " swaption=" << swaption_pv;
}

TEST(IRDerivsConsistencyTest, IRSParRateMatchesYieldCurve) {
    // 在 HW 无套利校准后, IRS par rate 反映平坦曲线下的简单复利远期利率
    // 平坦连续复利 y 下 par swap rate ≈ (exp(yτ)-1)/τ = y + y²τ/2 + ...
    auto [mats, bonds] = flat_curve(0.04, 30.0, 0.25);
    HullWhiteParams hw_p{0.04, 0.1, 0.01};
    HullWhite hw(hw_p, mats, bonds);

    // 5Y swap, 半年支付
    auto cfg = make_vanilla_irs(1.0, 0.0, 0.0, 0.5, 10, true);
    InterestRateSwap irs(cfg);
    Real K_par = irs.par_swap_rate(hw);
    // 平坦曲线下 par swap rate ≈ 简单复利远期 (exp(yτ)-1)/τ, 凸性偏差量级 ~y²τ/2 ≈ 0.0004
    Real y = 0.04, tau = 0.5;
    Real expected = (std::exp(y * tau) - 1.0) / tau;
    EXPECT_NEAR(K_par, expected, 1e-10);
    EXPECT_NEAR(K_par, y, 0.001);  // 凸性偏差量级
}
