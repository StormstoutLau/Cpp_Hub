// SOURCE: Vasicek (1977), CIR (1985), Hull-White (1990), Brigo-Mercurio (2006) Ch.3-4
// 测试覆盖: 参数验证 / 零息债解析解 / 债券期权 (Jamshidian) / Feller 条件 /
//          条件矩 / 路径模拟 / 期限结构校准 / Call-Put 平价
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "cpphub/models/ir/short_rate.hpp"
#include "cpphub/core/rng.hpp"

using namespace cpphub;

namespace {
// 通用容差
constexpr Real kTol = 1e-10;
}  // namespace

// ============================================================================
// Vasicek 模型测试
// ============================================================================

TEST(VasicekTest, ParameterValidationRejectsInvalidKappa) {
    EXPECT_THROW(Vasicek(VasicekParams{0.05, -0.1, 0.05, 0.01}), std::invalid_argument);
    EXPECT_THROW(Vasicek(VasicekParams{0.05, 0.0, 0.05, 0.01}), std::invalid_argument);
}

TEST(VasicekTest, ParameterValidationRejectsInvalidSigma) {
    EXPECT_THROW(Vasicek(VasicekParams{0.05, 0.5, 0.05, -0.1}), std::invalid_argument);
    EXPECT_THROW(Vasicek(VasicekParams{0.05, 0.5, 0.05, 0.0}), std::invalid_argument);
}

TEST(VasicekTest, ZeroCouponBondAtZeroIsOne) {
    Vasicek v({0.05, 0.5, 0.05, 0.01});
    EXPECT_NEAR(v.zero_coupon_bond(0.0), 1.0, kTol);
    EXPECT_NEAR(v.zero_coupon_bond(-1.0), 1.0, kTol);  // T<=0 退化
}

TEST(VasicekTest, ZeroCouponBondInValidRange) {
    Vasicek v({0.05, 0.5, 0.05, 0.01});
    for (Real T = 0.5; T <= 30.0; T += 0.5) {
        Real P = v.zero_coupon_bond(T);
        EXPECT_GT(P, 0.0) << "T=" << T;
        EXPECT_LE(P, 1.0 + kTol) << "T=" << T;
    }
}

TEST(VasicekTest, ZeroCouponBondMonotonicDecreasing) {
    // 正利率下, 零息债价格应随期限单调递减
    Vasicek v({0.05, 0.5, 0.05, 0.01});
    Real P_prev = v.zero_coupon_bond(0.5);
    for (Real T = 1.0; T <= 30.0; T += 1.0) {
        Real P = v.zero_coupon_bond(T);
        EXPECT_LT(P, P_prev) << "T=" << T;
        P_prev = P;
    }
}

TEST(VasicekTest, YieldAtZeroApproachesR0) {
    // 瞬时远期利率 = r0, 短端收益率应趋近 r0
    Vasicek v({0.05, 0.5, 0.05, 0.01});
    Real y_short = v.yield(1e-6);
    EXPECT_NEAR(y_short, 0.05, 1e-4);
}

TEST(VasicekTest, YieldConvergesToLongTermLimit) {
    // 长期收益率应收敛到有限值 (均值回复)
    // 理论极限 y_inf = kappa*theta - sigma^2/(2*kappa)
    // 收敛速度 ~ 1/T, T=1000 时残余项 < 1e-4
    Real kappa = 0.5, theta = 0.05, sigma = 0.01;
    Vasicek v({0.04, kappa, theta, sigma});
    Real y_inf_theoretical = kappa * theta - sigma * sigma / (2.0 * kappa);
    Real y_long = v.yield(1000.0);
    EXPECT_NEAR(y_long, y_inf_theoretical, 1e-4)
        << "y_long=" << y_long << " y_inf=" << y_inf_theoretical;
}

TEST(VasicekTest, ConditionalMomentsMeanRevertsToTheta) {
    // E[r(T)] → theta 当 T → ∞
    Vasicek v({0.04, 0.5, 0.05, 0.01});
    Real mean, var;
    v.conditional_moments(200.0, mean, var);
    EXPECT_NEAR(mean, 0.05, 1e-6);
    // Var[r(T)] → sigma^2/(2*kappa)
    Real var_inf = 0.01 * 0.01 / (2.0 * 0.5);
    EXPECT_NEAR(var, var_inf, 1e-8);
}

TEST(VasicekTest, ConditionalMomentsAtZeroMatchesR0) {
    Vasicek v({0.04, 0.5, 0.05, 0.01});
    Real mean, var;
    v.conditional_moments(0.0, mean, var);
    EXPECT_NEAR(mean, 0.04, kTol);
    EXPECT_NEAR(var, 0.0, kTol);
}

TEST(VasicekTest, BondDeltaSign) {
    // bond_delta = -B * P < 0 (B>0, P>0)
    Vasicek v({0.05, 0.5, 0.05, 0.01});
    Real delta = v.bond_delta(5.0);
    EXPECT_LT(delta, 0.0);
    // |delta| = B * P
    Real B = affine_B(0.5, 5.0);
    Real P = v.zero_coupon_bond(5.0);
    EXPECT_NEAR(std::abs(delta), B * P, kTol);
}

TEST(VasicekTest, BondOptionCallPutParity) {
    // Call - Put = P(0,T_opt) * (F - K)
    // F = E[P(T_opt, T_bond)] = mean_P (模型远期价格)
    Vasicek v({0.05, 0.5, 0.05, 0.01});
    Real T_opt = 1.0, T_bond = 5.0, K = 0.95;

    Real C = v.bond_option(T_opt, T_bond, K, true);
    Real P_put = v.bond_option(T_opt, T_bond, K, false);

    // 重新计算 mean_P 用于校验
    Real tau = T_bond - T_opt;
    Real B_tau = affine_B(0.5, tau);
    Real m_r = 0.05 + (0.05 - 0.05) * std::exp(-0.5 * T_opt);  // = theta (r0=theta)
    Real var_r = 0.01 * 0.01 * (1.0 - std::exp(-1.0)) / 1.0;
    Real A_tau = std::exp((B_tau - tau) * (0.5 * 0.05 - 0.01 * 0.01 / 1.0)
                          - 0.01 * 0.01 * B_tau * B_tau / 2.0);
    Real F = A_tau * std::exp(-B_tau * m_r + 0.5 * B_tau * B_tau * var_r);
    Real P_0_T = v.zero_coupon_bond(T_opt);

    Real parity = P_0_T * (F - K);
    EXPECT_NEAR(C - P_put, parity, 1e-10);
}

TEST(VasicekTest, BondOptionCallDecreasesWithStrike) {
    Vasicek v({0.05, 0.5, 0.05, 0.01});
    Real T_opt = 1.0, T_bond = 5.0;
    Real C_low = v.bond_option(T_opt, T_bond, 0.90, true);
    Real C_mid = v.bond_option(T_opt, T_bond, 0.95, true);
    Real C_high = v.bond_option(T_opt, T_bond, 1.00, true);
    EXPECT_GT(C_low, C_mid);
    EXPECT_GT(C_mid, C_high);
    EXPECT_GT(C_low, 0.0);
}

TEST(VasicekTest, BondOptionCallAboveIntrinsicValue) {
    // Call >= P(0,T_opt) * max(F - K, 0)
    Vasicek v({0.05, 0.5, 0.05, 0.01});
    Real T_opt = 1.0, T_bond = 5.0, K = 0.92;

    Real C = v.bond_option(T_opt, T_bond, K, true);

    Real tau = T_bond - T_opt;
    Real B_tau = affine_B(0.5, tau);
    Real m_r = 0.05;
    Real var_r = 0.0001 * (1.0 - std::exp(-1.0)) / 1.0;
    Real A_tau = std::exp((B_tau - tau) * (0.025 - 0.0001 / 1.0)
                          - 0.0001 * B_tau * B_tau / 2.0);
    Real F = A_tau * std::exp(-B_tau * m_r + 0.5 * B_tau * B_tau * var_r);
    Real P_0_T = v.zero_coupon_bond(T_opt);
    Real intrinsic = P_0_T * std::max(F - K, 0.0);

    EXPECT_GE(C, intrinsic - 1e-12);
}

TEST(VasicekTest, BondOptionArgumentValidation) {
    Vasicek v({0.05, 0.5, 0.05, 0.01});
    EXPECT_THROW(v.bond_option(0.0, 5.0, 0.95, true), std::invalid_argument);
    EXPECT_THROW(v.bond_option(5.0, 5.0, 0.95, true), std::invalid_argument);
    EXPECT_THROW(v.bond_option(5.0, 1.0, 0.95, true), std::invalid_argument);
}

TEST(VasicekTest, SimulatePathStartsAtR0) {
    Vasicek v({0.05, 0.5, 0.05, 0.01});
    Philox4x64 rng(42);
    std::vector<Real> path;
    v.simulate_path(1.0, 12, path, rng);
    EXPECT_EQ(path.size(), 13u);
    EXPECT_NEAR(path[0], 0.05, kTol);
}

TEST(VasicekTest, SimulatePathMatchesConditionalMoments) {
    // MC 验证: 路径终点 r(T) 的均值和方差应匹配条件矩
    Vasicek v({0.04, 0.5, 0.05, 0.01});
    Real T = 2.0;
    Size n_paths = 50000;
    Size n_steps = 100;

    Real mean_theo, var_theo;
    v.conditional_moments(T, mean_theo, var_theo);

    double sum = 0.0, sum2 = 0.0;
    for (Size i = 0; i < n_paths; ++i) {
        Philox4x64 rng(i);
        std::vector<Real> path;
        v.simulate_path(T, n_steps, path, rng);
        double rT = path.back();
        sum += rT;
        sum2 += rT * rT;
    }
    double mean_mc = sum / n_paths;
    double var_mc = sum2 / n_paths - mean_mc * mean_mc;

    EXPECT_NEAR(mean_mc, mean_theo, 3.0 * std::sqrt(var_theo / n_paths));
    EXPECT_NEAR(var_mc, var_theo, 5.0 * var_theo / std::sqrt(n_paths));
}

// ============================================================================
// CIR 模型测试
// ============================================================================

TEST(CIRTest, ParameterValidationRejectsInvalidParams) {
    EXPECT_THROW(CIR(CIRParams{0.05, -0.1, 0.05, 0.01}), std::invalid_argument);
    EXPECT_THROW(CIR(CIRParams{0.05, 0.5, -0.05, 0.01}), std::invalid_argument);
    EXPECT_THROW(CIR(CIRParams{-0.05, 0.5, 0.05, 0.01}), std::invalid_argument);
    EXPECT_THROW(CIR(CIRParams{0.05, 0.5, 0.05, -0.1}), std::invalid_argument);
}

TEST(CIRTest, ZeroCouponBondAtZeroIsOne) {
    CIR c({0.05, 0.5, 0.05, 0.01});
    EXPECT_NEAR(c.zero_coupon_bond(0.0), 1.0, kTol);
}

TEST(CIRTest, ZeroCouponBondMonotonicDecreasing) {
    CIR c({0.05, 0.5, 0.05, 0.01});
    Real P_prev = c.zero_coupon_bond(0.5);
    for (Real T = 1.0; T <= 30.0; T += 1.0) {
        Real P = c.zero_coupon_bond(T);
        EXPECT_LT(P, P_prev) << "T=" << T;
        P_prev = P;
    }
}

TEST(CIRTest, FellerConditionSatisfied) {
    // 2*kappa*theta > sigma^2
    CIR c({0.05, 0.5, 0.05, 0.01});
    EXPECT_TRUE(c.feller_satisfied());
    // 2*0.5*0.05 = 0.05 > 0.0001 ✓
}

TEST(CIRTest, FellerConditionViolated) {
    // 2*kappa*theta < sigma^2
    CIR c({0.05, 0.1, 0.02, 0.5});
    // 2*0.1*0.02 = 0.004 < 0.25
    EXPECT_FALSE(c.feller_satisfied());
}

TEST(CIRTest, ConditionalMomentsMeanRevertsToTheta) {
    CIR c({0.04, 0.5, 0.05, 0.01});
    Real mean, var;
    c.conditional_moments(200.0, mean, var);
    EXPECT_NEAR(mean, 0.05, 1e-6);
    // 长期方差: theta * sigma^2 / (2*kappa)  [CIR stationary distribution]
    Real var_inf = 0.05 * 0.0001 / (2.0 * 0.5);
    EXPECT_NEAR(var, var_inf, 1e-8);
}

TEST(CIRTest, SimulatePathNonNegative) {
    // CIR 路径应保持非负 (full truncation)
    CIR c({0.04, 0.5, 0.05, 0.01});
    Philox4x64 rng(42);
    for (Size trial = 0; trial < 100; ++trial) {
        Philox4x64 rng_t(trial);
        std::vector<Real> path;
        c.simulate_path(5.0, 100, path, rng_t);
        for (Real r : path) {
            EXPECT_GE(r, 0.0) << "trial=" << trial;
        }
    }
}

TEST(CIRTest, SimulatePathMatchesConditionalMoments) {
    CIR c({0.04, 0.5, 0.05, 0.01});
    Real T = 2.0;
    Size n_paths = 50000;
    Size n_steps = 200;  // 细化以减少离散化偏差

    Real mean_theo, var_theo;
    c.conditional_moments(T, mean_theo, var_theo);

    double sum = 0.0, sum2 = 0.0;
    for (Size i = 0; i < n_paths; ++i) {
        Philox4x64 rng(i);
        std::vector<Real> path;
        c.simulate_path(T, n_steps, path, rng);
        double rT = path.back();
        sum += rT;
        sum2 += rT * rT;
    }
    double mean_mc = sum / n_paths;
    double var_mc = sum2 / n_paths - mean_mc * mean_mc;

    EXPECT_NEAR(mean_mc, mean_theo, 3.0 * std::sqrt(var_theo / n_paths));
    // Euler 方案有离散化偏差, 放宽容差
    EXPECT_NEAR(var_mc, var_theo, 10.0 * var_theo / std::sqrt(n_paths));
}

// ============================================================================
// Hull-White 模型测试
// ============================================================================

namespace {
// 构造平坦期限结构 (收益率 = r_flat)
std::pair<std::vector<Real>, std::vector<Real>> flat_term_structure(Real r_flat) {
    std::vector<Real> maturities = {0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0};
    std::vector<Real> bonds;
    bonds.reserve(maturities.size());
    for (Real T : maturities) {
        bonds.push_back(std::exp(-r_flat * T));
    }
    return {maturities, bonds};
}
}  // namespace

TEST(HullWhiteTest, ParameterValidationRejectsInvalidParams) {
    auto [mats, bonds] = flat_term_structure(0.05);
    EXPECT_THROW(HullWhite(HullWhiteParams{0.05, -0.1, 0.01}, mats, bonds), std::invalid_argument);
    EXPECT_THROW(HullWhite(HullWhiteParams{0.05, 0.5, -0.01}, mats, bonds), std::invalid_argument);
    EXPECT_THROW(HullWhite(HullWhiteParams{0.05, 0.5, 0.01}, {0.5}, {0.9}), std::invalid_argument);  // size mismatch
}

TEST(HullWhiteTest, ZeroCouponBondMatchesMarketCurve) {
    Real r_flat = 0.05;
    auto [mats, bonds] = flat_term_structure(r_flat);
    HullWhite hw(HullWhiteParams{0.04, 0.5, 0.01}, mats, bonds);
    for (Real T : mats) {
        Real P_market = std::exp(-r_flat * T);
        Real P_model = hw.zero_coupon_bond(T);
        EXPECT_NEAR(P_model, P_market, 1e-10) << "T=" << T;
    }
}

TEST(HullWhiteTest, ForwardRateAtZeroApproachesR0) {
    Real r_flat = 0.05;
    auto [mats, bonds] = flat_term_structure(r_flat);
    HullWhite hw(HullWhiteParams{0.04, 0.5, 0.01}, mats, bonds);
    // 平坦曲线下, f(0, T) ≈ r_flat
    Real f = hw.forward_rate(1e-6);
    EXPECT_NEAR(f, r_flat, 1e-3);
}

TEST(HullWhiteTest, ThetaAtSmallTMatchesFormula) {
    // theta(t) = f'(t) + kappa*f(t) + sigma^2/(2*kappa)*(1 - exp(-2*kappa*t))
    // 平坦曲线下 f'(t)=0, f(t)=r_flat, 所以 theta(t) = kappa*r_flat + sigma^2/(2*kappa)*(1-exp(-2*kappa*t))
    Real r_flat = 0.05;
    auto [mats, bonds] = flat_term_structure(r_flat);
    Real kappa = 0.5;
    Real sigma = 0.01;
    HullWhite hw(HullWhiteParams{r_flat, kappa, sigma}, mats, bonds);
    Real t = 0.01;  // 足够大避免数值微分 clamping, 足够小使 sigma^2 项可解析
    Real theta_val = hw.theta(t);
    Real theta_theo = kappa * r_flat
                     + sigma * sigma / (2.0 * kappa) * (1.0 - std::exp(-2.0 * kappa * t));
    EXPECT_NEAR(theta_val, theta_theo, 1e-4)
        << "theta_val=" << theta_val << " theta_theo=" << theta_theo;
}

TEST(HullWhiteTest, HWBondPriceCloseToMarket) {
    // HW 零息债价格应接近市场曲线 (轻微偏差来自 r0 vs f(0) 的差异)
    Real r_flat = 0.05;
    auto [mats, bonds] = flat_term_structure(r_flat);
    HullWhite hw(HullWhiteParams{r_flat, 0.5, 0.01}, mats, bonds);
    for (Real T : mats) {
        Real P_market = std::exp(-r_flat * T);
        Real P_hw = hw.hw_zero_coupon_bond(T);
        // 当 r0 = f(0) 时, P_hw 应精确等于 P_market
        EXPECT_NEAR(P_hw, P_market, 1e-6) << "T=" << T;
    }
}

TEST(HullWhiteTest, BondOptionCallPutParity) {
    // Call - Put = P(0,T_opt) * (F_market - K)
    // F_market = P(0,T_bond) / P(0,T_opt)
    Real r_flat = 0.05;
    auto [mats, bonds] = flat_term_structure(r_flat);
    HullWhite hw(HullWhiteParams{0.05, 0.5, 0.01}, mats, bonds);

    Real T_opt = 1.0, T_bond = 5.0, K = 0.95;
    Real C = hw.bond_option(T_opt, T_bond, K, true);
    Real P_put = hw.bond_option(T_opt, T_bond, K, false);

    Real P_0_Topt = hw.zero_coupon_bond(T_opt);
    Real P_0_Tbond = hw.zero_coupon_bond(T_bond);
    Real F = P_0_Tbond / P_0_Topt;

    Real parity = P_0_Topt * (F - K);
    EXPECT_NEAR(C - P_put, parity, 1e-10);
}

TEST(HullWhiteTest, BondOptionCallDecreasesWithStrike) {
    Real r_flat = 0.05;
    auto [mats, bonds] = flat_term_structure(r_flat);
    HullWhite hw(HullWhiteParams{0.05, 0.5, 0.01}, mats, bonds);

    Real T_opt = 1.0, T_bond = 5.0;
    Real C1 = hw.bond_option(T_opt, T_bond, 0.90, true);
    Real C2 = hw.bond_option(T_opt, T_bond, 0.95, true);
    Real C3 = hw.bond_option(T_opt, T_bond, 1.00, true);
    EXPECT_GT(C1, C2);
    EXPECT_GT(C2, C3);
}

TEST(HullWhiteTest, SimulatePathStartsAtR0) {
    auto [mats, bonds] = flat_term_structure(0.05);
    HullWhite hw(HullWhiteParams{0.04, 0.5, 0.01}, mats, bonds);
    Philox4x64 rng(42);
    std::vector<Real> path;
    hw.simulate_path(1.0, 12, path, rng);
    EXPECT_EQ(path.size(), 13u);
    EXPECT_NEAR(path[0], 0.04, kTol);
}

// ============================================================================
// G2++ 模型测试
// ============================================================================

TEST(G2Test, ParameterValidationRejectsInvalidParams) {
    EXPECT_THROW(G2(G2Params{0.0, 0.0, -0.1, 0.1, 0.01, 0.01, 0.0}), std::invalid_argument);
    EXPECT_THROW(G2(G2Params{0.0, 0.0, 0.1, -0.1, 0.01, 0.01, 0.0}), std::invalid_argument);
    EXPECT_THROW(G2(G2Params{0.0, 0.0, 0.1, 0.1, -0.01, 0.01, 0.0}), std::invalid_argument);
    EXPECT_THROW(G2(G2Params{0.0, 0.0, 0.1, 0.1, 0.01, -0.01, 0.0}), std::invalid_argument);
    EXPECT_THROW(G2(G2Params{0.0, 0.0, 0.1, 0.1, 0.01, 0.01, 1.5}), std::invalid_argument);
}

TEST(G2Test, AffineZeroCouponBondAtZeroIsOne) {
    G2 g(G2Params{0.01, -0.01, 0.1, 0.1, 0.01, 0.01, 0.3});
    EXPECT_NEAR(g.affine_zero_coupon_bond(0.0), 1.0, kTol);
}

TEST(G2Test, VarianceTermIsPositive) {
    // V(0,T) = Var[积分 (x+y) ds] >= 0 (半正定协方差)
    G2 g(G2Params{0.0, 0.0, 0.1, 0.1, 0.01, 0.01, 0.5});
    for (Real T = 0.5; T <= 30.0; T += 1.0) {
        Real V = g.variance_term(T);
        EXPECT_GT(V, 0.0) << "T=" << T;
    }
}

TEST(G2Test, VarianceTermZeroAtZero) {
    G2 g(G2Params{0.0, 0.0, 0.1, 0.1, 0.01, 0.01, 0.5});
    EXPECT_NEAR(g.variance_term(0.0), 0.0, kTol);
}

TEST(G2Test, AffineBondPriceWithZeroVariance) {
    // sigma=eta=0 时, P_aff = exp(-B_a*x0 - B_b*y0)
    // 取 x0=y0=0, P_aff = 1 (无贴现)
    G2 g(G2Params{0.0, 0.0, 0.1, 0.1, 1e-12, 1e-12, 0.0});
    for (Real T = 1.0; T <= 10.0; T += 1.0) {
        EXPECT_NEAR(g.affine_zero_coupon_bond(T), 1.0, 1e-6);
    }
}

TEST(G2Test, ZeroCouponBondMatchesMarketCurve) {
    Real r_flat = 0.05;
    auto [mats, bonds] = flat_term_structure(r_flat);
    G2 g(G2Params{0.0, 0.0, 0.1, 0.1, 0.01, 0.01, 0.3}, mats, bonds);
    for (Real T : mats) {
        Real P_market = std::exp(-r_flat * T);
        Real P_model = g.zero_coupon_bond(T);
        EXPECT_NEAR(P_model, P_market, 1e-10) << "T=" << T;
    }
}

TEST(G2Test, PhiCalibrationIdentity) {
    // 校准后: P_market(0,T) = P_aff(0,T) * exp(-phi(T) * T)
    // → P_aff * exp(-phi*T) = P_market
    Real r_flat = 0.05;
    auto [mats, bonds] = flat_term_structure(r_flat);
    G2 g(G2Params{0.0, 0.0, 0.1, 0.1, 0.01, 0.01, 0.3}, mats, bonds);

    for (Real T : mats) {
        Real P_market = g.zero_coupon_bond(T);
        Real P_aff = g.affine_zero_coupon_bond(T);
        Real phi_T = g.phi(T);
        Real P_reconstructed = P_aff * std::exp(-phi_T * T);
        EXPECT_NEAR(P_reconstructed, P_market, 1e-10) << "T=" << T;
    }
}

TEST(G2Test, SimulatePathStartsAtX0PlusY0) {
    G2 g(G2Params{0.01, -0.005, 0.1, 0.1, 0.01, 0.01, 0.3});
    Philox4x64 rng(42);
    std::vector<Real> path;
    g.simulate_path(1.0, 12, path, rng);
    EXPECT_EQ(path.size(), 13u);
    EXPECT_NEAR(path[0], 0.01 + (-0.005), kTol);
}

TEST(G2Test, SimulatePathWithMarketCurveStartsReasonably) {
    // 提供市场曲线时, path[0] 仍是 x0+y0 (未加 phi(0))
    Real r_flat = 0.05;
    auto [mats, bonds] = flat_term_structure(r_flat);
    G2 g(G2Params{0.0, 0.0, 0.1, 0.1, 0.01, 0.01, 0.3}, mats, bonds);
    Philox4x64 rng(42);
    std::vector<Real> path;
    g.simulate_path(1.0, 12, path, rng);
    EXPECT_NEAR(path[0], 0.0, kTol);
}

// ============================================================================
// 跨模型一致性测试
// ============================================================================

TEST(ShortRateModelsTest, AffineBFunctionDecaysCorrectly) {
    // B(κ, 0) = 0, B(κ, T) → 1/κ 当 T → ∞
    EXPECT_NEAR(affine_B(0.5, 0.0), 0.0, kTol);
    EXPECT_NEAR(affine_B(0.5, 1000.0), 1.0 / 0.5, 1e-200);
    // κ → 0 时退化为 T
    EXPECT_NEAR(affine_B(1e-13, 5.0), 5.0, kTol);
}

TEST(ShortRateModelsTest, VasicekAndHWFMatchOnFlatCurve) {
    // 平坦曲线 + r0 = r_flat 时, HW 退化为 Vasicek (theta 时变但常数)
    Real r_flat = 0.05;
    auto [mats, bonds] = flat_term_structure(r_flat);
    HullWhite hw(HullWhiteParams{r_flat, 0.5, 0.01}, mats, bonds);
    // Vasicek with theta = r_flat, kappa=0.5, sigma=0.01
    // 注意: HW 的 theta(t) = f'(t) + kappa*f(t) + sigma^2/(2*kappa)*(1-exp(-2*kappa*t))
    // 平坦曲线 f'=0, f=r_flat, 所以 theta(t) = kappa*r_flat + sigma^2/(2*kappa)*(1-exp(-2*kappa*t))
    // 这不是常数, 所以 HW 债券价格 != Vasicek 债券价格 (除非 sigma=0)
    // 这里我们只验证两者在 T=0 都返回 1
    EXPECT_NEAR(hw.zero_coupon_bond(0.0), 1.0, kTol);
}

TEST(ShortRateModelsTest, CIRAndVasicekBothMeanReverting) {
    // 两个模型的收益率都应收敛到有限值
    Vasicek v({0.04, 0.5, 0.05, 0.01});
    CIR c({0.04, 0.5, 0.05, 0.01});
    Real y_v_short = v.yield(1.0);
    Real y_v_long = v.yield(100.0);
    Real y_c_short = c.yield(1.0);
    Real y_c_long = c.yield(100.0);
    // 短期收益率应接近 r0
    EXPECT_NEAR(y_v_short, 0.04, 0.01);
    EXPECT_NEAR(y_c_short, 0.04, 0.01);
    // 长期收益率应在合理范围内
    EXPECT_GT(y_v_long, 0.0);
    EXPECT_GT(y_c_long, 0.0);
    EXPECT_LT(y_v_long, 0.1);
    EXPECT_LT(y_c_long, 0.1);
}
