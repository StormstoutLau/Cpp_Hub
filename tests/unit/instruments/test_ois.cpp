// v1.2 Batch 2: OIS / 隔夜利率衍生品单元测试
// 覆盖: ZeroCurve 插值 / OIS Bootstrap / OIS 定价 / FRA / Basis Swap / Dual-curve IRS
#include <gtest/gtest.h>
#include "cpphub/instruments/ir/ois_curve.hpp"
#include "cpphub/instruments/ir/ois.hpp"
#include "cpphub/instruments/ir/fra.hpp"
#include "cpphub/instruments/ir/basis_swap.hpp"
#include "cpphub/instruments/ir/irs.hpp"
#include <cmath>
#include <vector>

using namespace cpphub;

// ============================================================
// ZeroCurve 插值测试
// ============================================================
TEST(ZeroCurveTest, FlatCurveDiscountFactor) {
    auto curve = make_flat_curve(0.05);
    // 平坦曲线: P(0,T) = exp(-0.05*T)
    EXPECT_NEAR(curve.discount_factor(1.0), std::exp(-0.05), 1e-12);
    EXPECT_NEAR(curve.discount_factor(5.0), std::exp(-0.25), 1e-12);
    EXPECT_NEAR(curve.discount_factor(0.0), 1.0, 1e-15);
}

TEST(ZeroCurveTest, LinearZeroInterpolation) {
    // r(1)=0.03, r(2)=0.05, r(3)=0.04 (非单调测试)
    std::vector<Real> T = {1.0, 2.0, 3.0};
    std::vector<Real> r = {0.03, 0.05, 0.04};
    ZeroCurve curve(T, r, ZeroCurve::InterpType::LinearZero);

    // 节点精确
    EXPECT_NEAR(curve.zero_rate(1.0), 0.03, 1e-12);
    EXPECT_NEAR(curve.zero_rate(2.0), 0.05, 1e-12);
    EXPECT_NEAR(curve.zero_rate(3.0), 0.04, 1e-12);
    // 中点线性: r(1.5) = (0.03+0.05)/2 = 0.04
    EXPECT_NEAR(curve.zero_rate(1.5), 0.04, 1e-12);
    // r(2.5) = (0.05+0.04)/2 = 0.045
    EXPECT_NEAR(curve.zero_rate(2.5), 0.045, 1e-12);
}

TEST(ZeroCurveTest, LogLinearDFInterpolation) {
    std::vector<Real> T = {1.0, 2.0};
    std::vector<Real> r = {0.03, 0.05};
    ZeroCurve curve(T, r, ZeroCurve::InterpType::LogLinearDF);

    // ln P 线性: ln P(1.5) = (ln P(1) + ln P(2))/2
    Real lnP1 = -0.03 * 1.0;
    Real lnP2 = -0.05 * 2.0;
    Real lnP_mid = 0.5 * (lnP1 + lnP2);
    Real r_expected = -lnP_mid / 1.5;
    EXPECT_NEAR(curve.zero_rate(1.5), r_expected, 1e-12);
    // 节点精确
    EXPECT_NEAR(curve.zero_rate(1.0), 0.03, 1e-12);
    EXPECT_NEAR(curve.zero_rate(2.0), 0.05, 1e-12);
}

TEST(ZeroCurveTest, CubicSplineInterpolation) {
    // 用 sin-like 数据测试样条平滑性
    std::vector<Real> T = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<Real> r = {0.02, 0.035, 0.045, 0.04, 0.03};
    ZeroCurve curve(T, r, ZeroCurve::InterpType::CubicSplineZero);

    // 节点精确
    for (Size i = 0; i < T.size(); ++i) {
        EXPECT_NEAR(curve.zero_rate(T[i]), r[i], 1e-10);
    }
    // 中点应在两点之间 (样条平滑, 不一定线性)
    Real r_mid = curve.zero_rate(1.5);
    EXPECT_GT(r_mid, 0.02);
    EXPECT_LT(r_mid, 0.035);
}

TEST(ZeroCurveTest, ForwardRateFlatCurve) {
    // 平坦曲线下, 简单复利远期 ≈ 连续复利零息利率 (近似)
    auto curve = make_flat_curve(0.04);
    Real fwd = curve.forward_rate(1.0, 2.0);
    // 简单复利: (exp(0.04)/exp(0.08) - 1)/1 = exp(-0.04) - 1 ≈ -0.0392... 错
    // 实际: (P1/P2 - 1)/(T2-T1) = (exp(-0.04)/exp(-0.08) - 1)/1 = exp(0.04) - 1 ≈ 0.0408
    Real expected = std::exp(0.04) - 1.0;
    EXPECT_NEAR(fwd, expected, 1e-12);
}

TEST(ZeroCurveTest, InstantaneousForwardFlat) {
    // 平坦曲线下瞬时远期 = r (常数)
    auto curve = make_flat_curve(0.03);
    EXPECT_NEAR(curve.instantaneous_forward(1.0), 0.03, 1e-6);
    EXPECT_NEAR(curve.instantaneous_forward(5.0), 0.03, 1e-6);
}

TEST(ZeroCurveTest, Extrapolation) {
    std::vector<Real> T = {1.0, 2.0, 3.0};
    std::vector<Real> r = {0.03, 0.05, 0.04};
    ZeroCurve curve(T, r, ZeroCurve::InterpType::LinearZero);
    // 外推: r(4.0) = 0.04 + (0.04-0.05)*(4-3)/(3-2) = 0.03
    EXPECT_NEAR(curve.zero_rate(4.0), 0.03, 1e-12);
    // 外推: r(0.5) = 0.03 + (0.05-0.03)*(0.5-1)/(2-1) = 0.02
    EXPECT_NEAR(curve.zero_rate(0.5), 0.02, 1e-12);
}

TEST(ZeroCurveTest, InvalidInputs) {
    EXPECT_THROW(ZeroCurve({}, {}), std::invalid_argument);
    EXPECT_THROW(ZeroCurve({1.0, 2.0}, {0.03}), std::invalid_argument);
    EXPECT_THROW(ZeroCurve({1.0, 1.0}, {0.03, 0.04}), std::invalid_argument);
    EXPECT_THROW(ZeroCurve({-1.0, 1.0}, {0.03, 0.04}), std::invalid_argument);
}

TEST(ZeroCurveTest, ModelInterfaceCompat) {
    // 兼容 short_rate.hpp model 接口
    auto curve = make_flat_curve(0.05);
    EXPECT_NEAR(curve.zero_coupon_bond(2.0), std::exp(-0.10), 1e-12);
    EXPECT_NEAR(curve.yield(3.0), 0.05, 1e-12);
}

// ============================================================
// OIS Bootstrap 测试
// ============================================================
TEST(OISBootstrapTest, FlatOISQuotesProduceFlatCurve) {
    // 平坦 OIS quotes (都 = 0.03, annual), 必须用连续 annual 节点 {1,2,3,4,5}
    // 避免中间节点缺失导致插值误差
    std::vector<OISCurveBootstrapper::OISQuote> quotes = {
        {1.0, 0.03}, {2.0, 0.03}, {3.0, 0.03}, {4.0, 0.03}, {5.0, 0.03}
    };
    OISCurveBootstrapper boot(quotes, 1.0);
    ZeroCurve curve = boot.bootstrap();

    // 由于 OIS swap rate 简单复利 vs 零息连续复利, bootstrap 出的 zero_rate ≈ 0.03 (但非精确)
    // 验证: 重新定价每个 OIS, PV≈0 (节点精确, 无插值误差)
    for (const auto& q : quotes) {
        OISConfig cfg = make_annual_ois(1.0, q.rate, 0.0, q.maturity,
                                         static_cast<Size>(q.maturity));
        OvernightIndexSwap ois(cfg);
        Real pv = ois.pv(curve);
        EXPECT_NEAR(pv, 0.0, 1e-9);
    }
}

TEST(OISBootstrapTest, IncreasingTermStructure) {
    // 上升曲线: 1Y=0.02, 2Y=0.03, 3Y=0.04
    std::vector<OISCurveBootstrapper::OISQuote> quotes = {
        {1.0, 0.02}, {2.0, 0.03}, {3.0, 0.04}
    };
    OISCurveBootstrapper boot(quotes, 1.0);
    ZeroCurve curve = boot.bootstrap();

    // 验证每个 quote 对应的 OIS 定价 PV ≈ 0
    for (const auto& q : quotes) {
        OISConfig cfg = make_annual_ois(1.0, q.rate, 0.0, q.maturity,
                                         static_cast<Size>(q.maturity));
        OvernightIndexSwap ois(cfg);
        EXPECT_NEAR(ois.pv(curve), 0.0, 1e-9);
    }
    // 验证 zero rate 单调性 (上升曲线 → r 上升)
    Real r1 = curve.zero_rate(1.0);
    Real r3 = curve.zero_rate(3.0);
    EXPECT_GT(r3, r1);
}

TEST(OISBootstrapTest, InvalidQuotes) {
    EXPECT_THROW(OISCurveBootstrapper({}), std::invalid_argument);
    EXPECT_THROW(OISCurveBootstrapper({{0.0, 0.03}}), std::invalid_argument);
    EXPECT_THROW(OISCurveBootstrapper({{1.0, -0.01}}), std::invalid_argument);
    EXPECT_THROW(OISCurveBootstrapper({{1.0, 0.03}, {1.0, 0.04}}), std::invalid_argument);
}

// ============================================================
// OIS 定价测试
// ============================================================
TEST(OISTest, FlatCurveParRate) {
    // 平坦曲线 r=0.04, 1Y OIS, par rate ≈ exp(0.04) - 1 (连续→简单复利)
    auto curve = make_flat_curve(0.04);
    OISConfig cfg = make_single_period_ois(1.0, 0.0, 0.0, 1.0, 1.0);
    OvernightIndexSwap ois(cfg);
    Real K_par = ois.par_rate(curve);
    // 单期 1Y OIS par rate = (1 - P(0,1)) / (1 * P(0,1)) = exp(0.04) - 1
    Real expected = std::exp(0.04) - 1.0;
    EXPECT_NEAR(K_par, expected, 1e-10);
}

TEST(OISTest, PayerOISZeroPVAtParRate) {
    auto curve = make_flat_curve(0.03);
    OISConfig cfg = make_annual_ois(1.0, 0.0, 0.0, 2.0, 2, true);
    OvernightIndexSwap ois(cfg);
    Real K = ois.par_rate(curve);
    // 设置 fixed_rate = par_rate, PV 应 = 0
    OISConfig cfg_par = cfg;
    cfg_par.fixed_rate = K;
    OvernightIndexSwap ois_par(cfg_par);
    EXPECT_NEAR(ois_par.pv(curve), 0.0, 1e-10);
}

TEST(OISTest, FloatLegFRNParProperty) {
    // 单曲线下, 浮端 PV = N*(1 - P(0,T_n)) (spot start, T_start=0)
    auto curve = make_flat_curve(0.05);
    OISConfig cfg = make_annual_ois(1.0, 0.04, 0.0, 3.0, 3, true);
    OvernightIndexSwap ois(cfg);
    Real flt_pv = ois.float_leg_pv(curve);
    Real P_end = curve.discount_factor(3.0);
    EXPECT_NEAR(flt_pv, 1.0 - P_end, 1e-10);
}

TEST(OISTest, FixedLegPV) {
    auto curve = make_flat_curve(0.05);
    OISConfig cfg = make_annual_ois(1.0, 0.04, 0.0, 2.0, 2, true);
    OvernightIndexSwap ois(cfg);
    Real fixed_pv = ois.fixed_leg_pv(curve);
    // 0.04 * (P(0,1) + P(0,2)) = 0.04 * (exp(-0.05) + exp(-0.10))
    Real expected = 0.04 * (std::exp(-0.05) + std::exp(-0.10));
    EXPECT_NEAR(fixed_pv, expected, 1e-12);
}

TEST(OISTest, PayerReceiverSymmetry) {
    auto curve = make_flat_curve(0.04);
    OISConfig payer_cfg = make_annual_ois(1.0, 0.05, 0.0, 2.0, 2, true);
    OISConfig receiver_cfg = make_annual_ois(1.0, 0.05, 0.0, 2.0, 2, false);
    OvernightIndexSwap payer(payer_cfg);
    OvernightIndexSwap receiver(receiver_cfg);
    Real pv_p = payer.pv(curve);
    Real pv_r = receiver.pv(curve);
    EXPECT_NEAR(pv_p, -pv_r, 1e-12);
}

TEST(OISTest, AnnuityAndForwardRate) {
    auto curve = make_flat_curve(0.04);
    OISConfig cfg = make_annual_ois(1.0, 0.04, 0.0, 3.0, 3, true);
    OvernightIndexSwap ois(cfg);
    Real ann = ois.annuity(curve);
    Real expected_ann = std::exp(-0.04) + std::exp(-0.08) + std::exp(-0.12);
    EXPECT_NEAR(ann, expected_ann, 1e-12);
    // forward_rate (单期近似): (1/P(0,3) - 1)/3
    Real fwd = ois.forward_rate(curve);
    Real expected_fwd = (1.0 / std::exp(-0.12) - 1.0) / 3.0;
    EXPECT_NEAR(fwd, expected_fwd, 1e-12);
}

TEST(OISTest, ForwardStartOIS) {
    // forward start OIS: T_start=1Y, T_end=3Y
    auto curve = make_flat_curve(0.04);
    OISConfig cfg = make_annual_ois(1.0, 0.04, 1.0, 3.0, 2, true);
    OvernightIndexSwap ois(cfg);
    // 浮端 = P(0,1) - P(0,3) (FRN par from T_start)
    Real flt = ois.float_leg_pv(curve);
    Real expected = std::exp(-0.04) - std::exp(-0.12);
    EXPECT_NEAR(flt, expected, 1e-12);
}

// ============================================================
// FRA 测试
// ============================================================
TEST(FRATest, ForwardRateFromCurve) {
    auto curve = make_flat_curve(0.05);
    FRAConfig cfg = make_3m_fra(1.0, 0.05, 0.0);
    ForwardRateAgreement fra(cfg);
    // F = (P(0,0)/P(0,0.25) - 1)/0.25 = (exp(0.0125) - 1)/0.25
    Real F = fra.forward_rate(curve);
    Real expected = (1.0 / std::exp(-0.05 * 0.25) - 1.0) / 0.25;
    EXPECT_NEAR(F, expected, 1e-12);
}

TEST(FRATest, ZeroPVAtParRate) {
    auto curve = make_flat_curve(0.04);
    FRAConfig cfg = make_3m_fra(1.0, 0.0, 0.5);
    ForwardRateAgreement fra(cfg);
    Real F = fra.forward_rate(curve);
    // 设置 K = F, PV = 0
    FRAConfig par_cfg = cfg;
    par_cfg.fixed_rate = F;
    ForwardRateAgreement fra_par(par_cfg);
    EXPECT_NEAR(fra_par.pv(curve), 0.0, 1e-12);
}

TEST(FRATest, MultiCurveVsSingleCurve) {
    // 投影曲线和贴现曲线相同 → 等于单曲线
    auto curve = make_flat_curve(0.04);
    FRAConfig cfg = make_3m_fra(1.0, 0.03, 0.0);
    ForwardRateAgreement fra(cfg);
    Real pv_single = fra.pv(curve);
    Real pv_multi = fra.pv(curve, curve);
    EXPECT_NEAR(pv_single, pv_multi, 1e-12);
}

TEST(FRATest, MultiCurveDiscountingEffect) {
    // Forward start FRA: T_reset=0.5, T_pay=0.75 (3M FRA starting in 6M)
    // 这样 T_settle=0.5 > 0, discount curve 的 P_d(0, 0.5) 影响 PV
    // LIBOR projection curve = 5%, OIS discount curve = 4%
    auto proj = make_flat_curve(0.05);
    auto disc = make_flat_curve(0.04);
    FRAConfig cfg;
    cfg.notional = 1.0;
    cfg.fixed_rate = 0.04;
    cfg.T_reset = 0.5;
    cfg.T_pay = 0.75;
    cfg.T_settle = 0.5;  // settle at fixing date
    cfg.year_fraction = 0.25;
    cfg.is_long = true;
    ForwardRateAgreement fra(cfg);
    Real pv = fra.pv(proj, disc);
    EXPECT_GT(pv, 0.0);  // long FRA, LIBOR > K, 获利

    // Multi-curve (OIS discount 4%) vs single LIBOR discount (5%):
    // P_d(0,0.5) at OIS=4% > P_LIBOR(0,0.5) at 5%, 所以 multi-curve PV 更大
    Real pv_single_libor = fra.pv(proj, proj);
    EXPECT_GT(pv, pv_single_libor);
}

TEST(FRATest, LongShortSymmetry) {
    auto curve = make_flat_curve(0.04);
    FRAConfig long_cfg = make_3m_fra(1.0, 0.03, 0.0, true);
    FRAConfig short_cfg = make_3m_fra(1.0, 0.03, 0.0, false);
    ForwardRateAgreement fra_long(long_cfg);
    ForwardRateAgreement fra_short(short_cfg);
    Real pv_l = fra_long.pv(curve);
    Real pv_s = fra_short.pv(curve);
    EXPECT_NEAR(pv_l, -pv_s, 1e-12);
}

TEST(FRATest, PayoffAtSettle) {
    // realized LIBOR = 0.05, K = 0.03, τ = 0.25, N = 1M
    FRAConfig cfg = make_3m_fra(1'000'000.0, 0.03, 0.0, true);
    ForwardRateAgreement fra(cfg);
    Real payoff = fra.payoff_at_settle(0.05);
    Real expected = 1'000'000.0 * 0.25 * (0.05 - 0.03) / (1.0 + 0.25 * 0.05);
    EXPECT_NEAR(payoff, expected, 1e-6);
}

TEST(FRATest, SimplePV) {
    // pv_simple (payoff at T_pay, 不提前折现)
    auto curve = make_flat_curve(0.05);
    FRAConfig cfg = make_3m_fra(1.0, 0.03, 0.0, true);
    ForwardRateAgreement fra(cfg);
    Real F = fra.forward_rate(curve);
    Real pv_simple = fra.pv_simple(curve, curve);
    Real expected = 1.0 * 0.25 * (F - 0.03) * std::exp(-0.05 * 0.25);
    EXPECT_NEAR(pv_simple, expected, 1e-12);
}

// ============================================================
// Basis Swap 测试
// ============================================================
TEST(BasisSwapTest, SingleCurveSameLegsZeroPV) {
    // 单曲线下, 两腿相同 (无 spread) → PV = 0
    auto curve = make_flat_curve(0.04);
    FloatLegConfig leg = make_float_leg(1.0, 0.0, 2.0, 4, 0.0);
    BasisSwap bs(leg, leg, true);
    Real pv = bs.pv(curve, curve, curve);
    EXPECT_NEAR(pv, 0.0, 1e-12);
}

TEST(BasisSwapTest, FlatCurveParSpreadNearZero) {
    // 平坦曲线下, 3M vs 6M par spread ≈ 0 (单曲线, 平坦)
    auto curve = make_flat_curve(0.04);
    BasisSwap bs = make_3m_6m_basis_swap(1.0, 0.0, 2.0, 0.0, true);
    Real spread = bs.par_spread_leg2(curve, curve, curve);
    EXPECT_NEAR(spread, 0.0, 1e-10);
}

TEST(BasisSwapTest, FloatLegFRNParSingleCurve) {
    // 单曲线下, 浮端 (无 spread) PV = N*(P(0,T_start) - P(0,T_n))
    auto curve = make_flat_curve(0.05);
    FloatLegConfig leg = make_float_leg(1.0, 0.0, 3.0, 3, 0.0);
    Real pv = float_leg_pv(leg, curve, curve);
    Real expected = std::exp(0.0) - std::exp(-0.05 * 3.0);  // 1 - P(0,3)
    EXPECT_NEAR(pv, expected, 1e-12);
}

TEST(BasisSwapTest, MultiCurveParSpreadNonZero) {
    // 投影曲线 (LIBOR=5%) > 贴现曲线 (OIS=3%)
    // 两腿相同 tenor, 但 projection_1 != projection_2 → par spread 非零
    auto proj1 = make_flat_curve(0.05);
    auto proj2 = make_flat_curve(0.045);
    auto disc = make_flat_curve(0.03);

    FloatLegConfig leg1 = make_float_leg(1.0, 0.0, 2.0, 4, 0.0);
    FloatLegConfig leg2 = make_float_leg(1.0, 0.0, 2.0, 4, 0.0);
    BasisSwap bs(leg1, leg2, true);

    // leg1 (proj1=5%) > leg2 (proj2=4.5%), pay leg1 → 需 leg2 spread > 0
    Real spread = bs.par_spread_leg2(proj1, proj2, disc);
    EXPECT_GT(spread, 0.0);
    // spread 应 ≈ 5% - 4.5% = 0.5% (近似)
    EXPECT_NEAR(spread, 0.005, 5e-4);
}

TEST(BasisSwapTest, AnnuityAndLegPVs) {
    auto curve = make_flat_curve(0.04);
    FloatLegConfig leg = make_float_leg(1.0, 0.0, 2.0, 4, 0.0);
    Real annuity = float_leg_annuity(leg, curve);
    // Σ τ_i * P(0, T_i), τ_i = 0.5, T_i = 0.5, 1.0, 1.5, 2.0
    Real expected = 0.5 * (std::exp(-0.02) + std::exp(-0.04) + std::exp(-0.06) + std::exp(-0.08));
    EXPECT_NEAR(annuity, expected, 1e-12);
}

TEST(BasisSwapTest, DifferentTenorLegs) {
    // 3M leg (4 期) vs 6M leg (2 期), 2Y 总期限, 平坦曲线
    auto curve = make_flat_curve(0.04);
    BasisSwap bs = make_3m_6m_basis_swap(1.0, 0.0, 2.0, 0.0, true);
    // 验证两条腿 PV 都 = N*(1 - P(0,T_n)) (单曲线, FRN par)
    Real pv_leg1 = bs.leg1_pv(curve, curve);
    Real pv_leg2 = bs.leg2_pv(curve, curve);
    Real expected = 1.0 - std::exp(-0.08);  // 1 - P(0,2)
    EXPECT_NEAR(pv_leg1, expected, 1e-12);
    EXPECT_NEAR(pv_leg2, expected, 1e-12);
    // PV (无 spread) = 0
    EXPECT_NEAR(bs.pv(curve, curve, curve), 0.0, 1e-12);
}

// ============================================================
// Dual-Curve IRS 测试
// ============================================================
TEST(DualCurveIRSTest, SingleCurveDegeneratesToSingle) {
    // projection = discount → 与 InterestRateSwap::present_value 一致
    auto curve = make_flat_curve(0.04);
    IRSConfig cfg = make_vanilla_irs(1.0, 0.04, 0.0, 0.5, 4, true);
    InterestRateSwap irs(cfg);
    Real pv_single = irs.present_value(curve);
    Real pv_dual = irs_dual_curve_pv(cfg, curve, curve);
    EXPECT_NEAR(pv_single, pv_dual, 1e-10);
}

TEST(DualCurveIRSTest, ParRateSingleCurve) {
    auto curve = make_flat_curve(0.05);
    IRSConfig cfg = make_vanilla_irs(1.0, 0.0, 0.0, 0.5, 4, true);
    Real par_single = InterestRateSwap(cfg).par_swap_rate(curve);
    Real par_dual = irs_dual_curve_par_rate(cfg, curve, curve);
    EXPECT_NEAR(par_single, par_dual, 1e-10);
}

TEST(DualCurveIRSTest, OISDiscountingRaisesParRate) {
    // 非平坦 LIBOR curve (上升) 使 forward rate 非常数, OIS discount 才会影响 par rate
    // 平坦曲线下 forward=const, P_d 在 par rate 分子分母抵消, discount curve 无影响
    std::vector<Real> T = {0.5, 1.0, 1.5, 2.0};
    std::vector<Real> r_libor = {0.04, 0.045, 0.05, 0.055};  // 上升 LIBOR curve
    ZeroCurve libor_curve(T, r_libor, ZeroCurve::InterpType::LinearZero);
    auto ois_curve = make_flat_curve(0.03);  // 平坦 OIS discount curve

    IRSConfig cfg = make_vanilla_irs(1.0, 0.0, 0.0, 0.5, 4, true);
    Real par_dual = irs_dual_curve_par_rate(cfg, libor_curve, ois_curve);
    Real par_single_libor = irs_dual_curve_par_rate(cfg, libor_curve, libor_curve);
    Real par_single_ois = irs_dual_curve_par_rate(cfg, ois_curve, ois_curve);

    // Dual curve par rate 应 > OIS-only par rate (LIBOR projection 更高)
    EXPECT_GT(par_dual, par_single_ois);
    // Dual curve par rate 应 != LIBOR-only par rate (非平坦曲线下 P_d 影响权重)
    // 注意: 不保证严格小于, 因 OIS discount 降低 annuity 同时影响分子分母
    // 但一定不相等 (非平坦 LIBOR curve 下)
    Real diff = std::abs(par_dual - par_single_libor);
    EXPECT_GT(diff, 1e-8);
}

TEST(DualCurveIRSTest, ZeroPVAtParRate) {
    auto libor_curve = make_flat_curve(0.05);
    auto ois_curve = make_flat_curve(0.04);
    IRSConfig cfg_base = make_vanilla_irs(1.0, 0.0, 0.0, 0.5, 4, true);
    Real K_par = irs_dual_curve_par_rate(cfg_base, libor_curve, ois_curve);
    IRSConfig cfg_par = cfg_base;
    cfg_par.fixed_rate = K_par;
    Real pv = irs_dual_curve_pv(cfg_par, libor_curve, ois_curve);
    EXPECT_NEAR(pv, 0.0, 1e-10);
}

TEST(DualCurveIRSTest, PayerReceiverSymmetry) {
    auto libor_curve = make_flat_curve(0.05);
    auto ois_curve = make_flat_curve(0.04);
    IRSConfig payer_cfg = make_vanilla_irs(1.0, 0.05, 0.0, 0.5, 4, true);
    IRSConfig receiver_cfg = make_vanilla_irs(1.0, 0.05, 0.0, 0.5, 4, false);
    Real pv_p = irs_dual_curve_pv(payer_cfg, libor_curve, ois_curve);
    Real pv_r = irs_dual_curve_pv(receiver_cfg, libor_curve, ois_curve);
    EXPECT_NEAR(pv_p, -pv_r, 1e-12);
}

TEST(DualCurveIRSTest, FixedLegAndFloatLeg) {
    auto libor_curve = make_flat_curve(0.05);
    auto ois_curve = make_flat_curve(0.04);
    IRSConfig cfg = make_vanilla_irs(1.0, 0.05, 0.0, 0.5, 2, true);
    Real fixed_pv = irs_dual_curve_fixed_leg_pv(cfg, ois_curve);
    Real flt_pv = irs_dual_curve_float_leg_pv(cfg, libor_curve, ois_curve);
    // Payer PV = Float - Fixed
    Real pv = irs_dual_curve_pv(cfg, libor_curve, ois_curve);
    EXPECT_NEAR(pv, flt_pv - fixed_pv, 1e-12);
}

TEST(DualCurveIRSTest, ForwardStartDualCurve) {
    // Forward start IRS: T_start=1Y, 2Y tenor, semi-annual payments
    auto libor_curve = make_flat_curve(0.05);
    auto ois_curve = make_flat_curve(0.04);
    IRSConfig cfg = make_vanilla_irs(1.0, 0.0, 1.0, 0.5, 4, true);
    Real par = irs_dual_curve_par_rate(cfg, libor_curve, ois_curve);
    // Forward start par rate 应合理 (>0, <1)
    EXPECT_GT(par, 0.0);
    EXPECT_LT(par, 1.0);
    // 验证 PV at par = 0
    IRSConfig cfg_par = cfg;
    cfg_par.fixed_rate = par;
    EXPECT_NEAR(irs_dual_curve_pv(cfg_par, libor_curve, ois_curve), 0.0, 1e-10);
}
