// SOURCE: Discovery 002 — Heston 特征函数分支切割在极端参数下失效
// SOURCE: Albrecher et al. (2007) "The little Heston trap" (主分支选择)
// SOURCE: Paper A (Peng Liu) "Conduction Intensity and Structural Fragility of Heston Volatility Surfaces"
//         — Feller 边界 δ_F = 2κθ - σ_v² 已证明为 soft crossover (Proposition 3)
//
// 测试目标 (Discovery 002 修正后):
//   - ρ → ±1: 接近 "no moment explosion" 边界 2κ = ρσ_v
//   - σ_v → 0: Heston 退化为 Black-Scholes
//   - τ → 0: 短到期 Riccati 解的渐近行为
//   - δ_F < 0: Feller 违反区域 (Paper A 明确排除, 但实际标定可能进入)
//
// 关键: 当前实现使用 log-of-ratio form (Albrecher 2007 Little Trap), 避免分支切割不连续
//       Paper A Proposition 3 已证明 Feller 边界是 soft crossover, 不构成裂缝
//       本测试聚焦其他极端参数是否产生数值不稳定
//
// 容差: 1e-10 (常规), 1e-8 (极端参数)
#include <gtest/gtest.h>
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include <cmath>
#include <complex>
#include <vector>

using namespace cpphub;

namespace {
constexpr Real kTolRegular = 1e-10;
constexpr Real kTolExtreme = 1e-8;
constexpr Real kTolDegenerate = 1e-6;  // σ_v → 0 退化场景

// Schoutens 标准参数集 (基准对照)
HestonCFParams schoutens_params() {
    return HestonCFParams{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, 0.0};
}

// σ_v → 0 时 Heston 的极限: v(t) 确定性均值回复
//   v(t) = θ + (v0 - θ) e^{-κt}
//   ∫₀^τ v(s) ds = θτ + (v0-θ)(1-e^{-κτ})/κ  (确定性积分方差)
//   φ(u) = exp(i*u*(ln S0 + (r-q)τ) - 0.5*(u² + iu) * ∫v ds)
Complex bs_characteristic_function(Complex u, Real tau, Real S0, const HestonCFParams& p) {
    Real integral_v = p.theta * tau
                       + (p.v0 - p.theta) * (1.0 - std::exp(-p.kappa * tau)) / p.kappa;
    Complex i(0, 1);
    Complex exponent = i * u * (std::log(S0) + (p.r - p.q) * tau)
                       - 0.5 * (u * u + i * u) * integral_v;
    return std::exp(exponent);
}

// 检查 cf 在 u 扫描下的连续性 (双侧差分)
Real cf_continuity_residual(Complex u, Real tau, Real S0, const HestonCFParams& p, Real eps = 1e-6) {
    Complex phi_fwd = heston_characteristic_function(u + eps, tau, S0, p);
    Complex phi_bwd = heston_characteristic_function(u - eps, tau, S0, p);
    Complex phi_mid = heston_characteristic_function(u, tau, S0, p);
    return std::abs(phi_fwd + phi_bwd - Real(2) * phi_mid);
}
}  // namespace

// ============================================================================
// 1. ρ → +1 极限 (接近 "no moment explosion" 边界 2κ = ρσ_v)
// ============================================================================
// Paper A Assumption 2: 2κ > ρσ_v 限制了 ρ 的上界
// 当 ρ → +1, 边界条件 2κ > σ_v 接近违反 (需 κ > σ_v/2)
// 测试: cf 在 ρ → +1 时是否仍连续可计算

TEST(HestonCFEdgeCases, RhoNearPlusOneCFComputable) {
    // ρ = 0.999: 接近 +1 但仍满足 2κ = 3.0 > ρσ_v = 0.999*0.3 = 0.2997
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, 0.999, 0.0, 0.0};
    Real tau = 1.0;
    Real S0 = 100.0;
    for (Real u_val = 0.1; u_val <= 10.0; u_val += 0.5) {
        Complex u(u_val, 0);
        Complex phi = heston_characteristic_function(u, tau, S0, p);
        EXPECT_TRUE(std::isfinite(std::real(phi))) << "u=" << u_val;
        EXPECT_TRUE(std::isfinite(std::imag(phi))) << "u=" << u_val;
        EXPECT_LE(std::abs(phi), 1.0 + kTolExtreme) << "u=" << u_val;
    }
}

TEST(HestonCFEdgeCases, RhoNearPlusOneNoBranchCutDiscontinuity) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, 0.999, 0.0, 0.0};
    Real tau = 1.0;
    Real S0 = 100.0;
    for (Real u_val = 0.1; u_val <= 20.0; u_val += 0.1) {
        Complex u(u_val, 0);
        Real residual = cf_continuity_residual(u, tau, S0, p);
        EXPECT_LT(residual, kTolExtreme) << "u=" << u_val;
    }
}

TEST(HestonCFEdgeCases, RhoNearPlusOneContinuousInRho) {
    // ρ 从 0.99 → 0.999 → 0.9999, cf 应连续变化
    Real tau = 1.0;
    Real S0 = 100.0;
    Complex u(1.0, 0);
    std::vector<Real> rhos = {0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999};
    std::vector<Complex> phis;
    for (Real rho : rhos) {
        HestonCFParams p{0.04, 1.5, 0.04, 0.3, rho, 0.0, 0.0};
        phis.push_back(heston_characteristic_function(u, tau, S0, p));
    }
    // 相邻 ρ 的 cf 变化应单调 (无突变)
    for (Size i = 1; i < phis.size(); ++i) {
        Real delta = std::abs(phis[i] - phis[i - 1]);
        EXPECT_LT(delta, 0.1) << "rho transition " << rhos[i - 1] << " -> " << rhos[i];
    }
}

// ============================================================================
// 2. ρ → -1 极限
// ============================================================================

TEST(HestonCFEdgeCases, RhoNearMinusOneCFComputable) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, -0.999, 0.0, 0.0};
    Real tau = 1.0;
    Real S0 = 100.0;
    for (Real u_val = 0.1; u_val <= 10.0; u_val += 0.5) {
        Complex u(u_val, 0);
        Complex phi = heston_characteristic_function(u, tau, S0, p);
        EXPECT_TRUE(std::isfinite(std::real(phi))) << "u=" << u_val;
        EXPECT_TRUE(std::isfinite(std::imag(phi))) << "u=" << u_val;
        EXPECT_LE(std::abs(phi), 1.0 + kTolExtreme) << "u=" << u_val;
    }
}

TEST(HestonCFEdgeCases, RhoNearMinusOneNoBranchCutDiscontinuity) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, -0.999, 0.0, 0.0};
    Real tau = 1.0;
    Real S0 = 100.0;
    for (Real u_val = 0.1; u_val <= 20.0; u_val += 0.1) {
        Complex u(u_val, 0);
        Real residual = cf_continuity_residual(u, tau, S0, p);
        EXPECT_LT(residual, kTolExtreme) << "u=" << u_val;
    }
}

// ============================================================================
// 3. σ_v → 0 退化 (Heston → Black-Scholes)
// ============================================================================
// 当 σ_v → 0, v(t) = v0 常数, Heston 退化为 BS with σ² = v0
// 测试: cf 是否平滑收敛到 BS 极限

TEST(HestonCFEdgeCases, SigmaVNearZeroConvergesToBS) {
    // σ_v → 0 时 Heston 应单调趋近 BS 极限 (确定性方差)
    // 不验证严格相等 (σ_v² 在分母导致数值精度损失), 而验证收敛趋势
    Real tau = 1.0;
    Real S0 = 100.0;
    Complex u(2.0, 0);  // 固定 u, 扫描 σ_v
    HestonCFParams p_base{0.04, 1.5, 0.04, 0.3, -0.5, 0.03, 0.01};
    Complex phi_bs = bs_characteristic_function(u, tau, S0, p_base);

    std::vector<Real> sigma_vs = {0.3, 0.1, 0.03, 0.01, 0.003};
    std::vector<Real> errors;
    for (Real sigma_v : sigma_vs) {
        HestonCFParams p = p_base;
        p.sigma = sigma_v;
        Complex phi_heston = heston_characteristic_function(u, tau, S0, p);
        Real err = std::abs(phi_heston - phi_bs);
        errors.push_back(err);
        EXPECT_TRUE(std::isfinite(err)) << "sigma_v=" << sigma_v;
    }
    // 误差应随 σ_v 减小而单调减小 (趋近 BS 极限)
    for (Size i = 1; i < errors.size(); ++i) {
        EXPECT_LT(errors[i], errors[i - 1])
            << "sigma_v transition " << sigma_vs[i - 1] << " -> " << sigma_vs[i];
    }
    // σ_v = 0.003 时误差应已较小 (< 0.05)
    EXPECT_LT(errors.back(), 0.05);
}

TEST(HestonCFEdgeCases, SigmaVNearZeroCFContinuous) {
    // σ_v 从 0.3 → 0.0001, cf 应连续变化 (无突变)
    Real tau = 1.0;
    Real S0 = 100.0;
    Complex u(1.0, 0);
    std::vector<Real> sigma_vs = {0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001};
    std::vector<Complex> phis;
    for (Real sigma_v : sigma_vs) {
        HestonCFParams p{0.04, 1.5, 0.04, sigma_v, -0.5, 0.0, 0.0};
        phis.push_back(heston_characteristic_function(u, tau, S0, p));
        EXPECT_TRUE(std::isfinite(std::real(phis.back()))) << "sigma_v=" << sigma_v;
        EXPECT_TRUE(std::isfinite(std::imag(phis.back()))) << "sigma_v=" << sigma_v;
    }
    // 相邻 σ_v 的 cf 变化应有限 (趋近 BS 极限, 无突变)
    for (Size i = 1; i < phis.size(); ++i) {
        Real delta = std::abs(phis[i] - phis[i - 1]);
        EXPECT_LT(delta, 0.5) << "sigma_v transition " << sigma_vs[i - 1] << " -> " << sigma_vs[i];
    }
    // σ_v = 1e-4 应接近 BS 极限 (σ_v² = 1e-8 不会下溢)
    HestonCFParams p_bs{0.04, 1.5, 0.04, 1e-4, -0.5, 0.0, 0.0};
    Complex phi_final = heston_characteristic_function(u, tau, S0, p_bs);
    Complex phi_bs = bs_characteristic_function(u, tau, S0, p_bs);
    EXPECT_NEAR(std::real(phi_final), std::real(phi_bs), kTolDegenerate);
    EXPECT_NEAR(std::imag(phi_final), std::imag(phi_bs), kTolDegenerate);
}

// ============================================================================
// 4. τ → 0 短到期
// ============================================================================
// τ → 0 时, Riccati 解的渐近行为
// 极限: cf(u, 0) = 1 (确定性极限)

TEST(HestonCFEdgeCases, ShortMaturityCFConvergesToSpot) {
    // τ → 0 时 S_T → S0, cf(u, 0) = E[e^{iu ln S_T}] → e^{iu ln S0} = S0^{iu}
    HestonCFParams p = schoutens_params();
    Real S0 = 100.0;
    Complex u(1.0, 0);
    Complex phi_limit = std::exp(Complex(0, 1) * u * std::log(S0));  // S0^{iu}
    for (Real tau : {1.0, 0.1, 0.01, 0.001, 1e-6}) {
        Complex phi = heston_characteristic_function(u, tau, S0, p);
        Real dist_to_limit = std::abs(phi - phi_limit);
        EXPECT_TRUE(std::isfinite(dist_to_limit)) << "tau=" << tau;
        // τ 减小时 cf 应趋近极限 (单调性不保证, 但距离应减小)
        EXPECT_LT(dist_to_limit, 1.0) << "tau=" << tau;
    }
    // 极短到期应非常接近 S0^{iu}
    Complex phi_tiny = heston_characteristic_function(u, 1e-8, S0, p);
    EXPECT_NEAR(std::real(phi_tiny), std::real(phi_limit), kTolExtreme);
    EXPECT_NEAR(std::imag(phi_tiny), std::imag(phi_limit), kTolExtreme);
}

TEST(HestonCFEdgeCases, ShortMaturityNoBranchCutDiscontinuity) {
    HestonCFParams p = schoutens_params();
    Real S0 = 100.0;
    Real tau = 1e-4;  // 短到期
    for (Real u_val = 0.1; u_val <= 20.0; u_val += 0.1) {
        Complex u(u_val, 0);
        Real residual = cf_continuity_residual(u, tau, S0, p);
        EXPECT_LT(residual, kTolExtreme) << "u=" << u_val << " tau=" << tau;
    }
}

// ============================================================================
// 5. Feller 违反区域 (δ_F < 0)
// ============================================================================
// Paper A 明确排除此区域 (Assumption 2: δ_F > 0)
// 但实际标定可能进入此区域, 测试 cf 是否仍可计算 (Albrecher 2007 算法仍适用)
// δ_F = 2κθ - σ_v²; 违反时 δ_F < 0

TEST(HestonCFEdgeCases, FellerViolationCFComputable) {
    // 2*0.5*0.04 = 0.04 < 0.3² = 0.09 → δ_F = -0.05 (Feller 违反)
    HestonCFParams p{0.04, 0.5, 0.04, 0.3, -0.7, 0.0, 0.0};
    Real tau = 1.0;
    Real S0 = 100.0;
    for (Real u_val = 0.1; u_val <= 10.0; u_val += 0.5) {
        Complex u(u_val, 0);
        Complex phi = heston_characteristic_function(u, tau, S0, p);
        EXPECT_TRUE(std::isfinite(std::real(phi))) << "u=" << u_val;
        EXPECT_TRUE(std::isfinite(std::imag(phi))) << "u=" << u_val;
        // |cf| ≤ 1 仍应满足 (cf 仍是特征函数)
        EXPECT_LE(std::abs(phi), 1.0 + kTolExtreme) << "u=" << u_val;
    }
}

TEST(HestonCFEdgeCases, FellerViolationNoBranchCutDiscontinuity) {
    HestonCFParams p{0.04, 0.5, 0.04, 0.3, -0.7, 0.0, 0.0};
    Real tau = 1.0;
    Real S0 = 100.0;
    for (Real u_val = 0.1; u_val <= 20.0; u_val += 0.1) {
        Complex u(u_val, 0);
        Real residual = cf_continuity_residual(u, tau, S0, p);
        EXPECT_LT(residual, kTolExtreme) << "u=" << u_val;
    }
}

TEST(HestonCFEdgeCases, FellerBoundarySoftCrossover) {
    // Paper A Proposition 3: Feller 边界 δ_F = 0 是 soft crossover
    // 测试: δ_F 从 +0.01 → 0 → -0.01, cf 应连续变化 (无突变)
    Real tau = 1.0;
    Real S0 = 100.0;
    Complex u(1.0, 0);
    Real kappa = 0.5;
    Real theta = 0.04;
    // δ_F = 2*kappa*theta - sigma_v² = 0.04 - sigma_v²
    // δ_F = +0.01 → sigma_v = sqrt(0.03) ≈ 0.1732
    // δ_F = 0     → sigma_v = sqrt(0.04) = 0.2
    // δ_F = -0.01 → sigma_v = sqrt(0.05) ≈ 0.2236
    std::vector<Real> sigma_vs = {0.15, 0.17, 0.19, 0.20, 0.21, 0.23};
    std::vector<Complex> phis;
    for (Real sigma_v : sigma_vs) {
        HestonCFParams p{0.04, kappa, theta, sigma_v, -0.5, 0.0, 0.0};
        phis.push_back(heston_characteristic_function(u, tau, S0, p));
    }
    // 相邻 σ_v 的 cf 变化应平滑 (跨过 Feller 边界无突变)
    for (Size i = 1; i < phis.size(); ++i) {
        Real delta = std::abs(phis[i] - phis[i - 1]);
        EXPECT_LT(delta, 0.1) << "sigma_v transition " << sigma_vs[i - 1] << " -> " << sigma_vs[i];
    }
}

// ============================================================================
// 6. 跨极端参数的连续性测试
// ============================================================================

TEST(HestonCFEdgeCases, CFContinuousInRhoScan) {
    // ρ 从 -0.999 → +0.999 扫描, cf 应连续变化
    HestonCFParams base = schoutens_params();
    Real tau = 1.0;
    Real S0 = 100.0;
    Complex u(1.0, 0);
    Complex phi_prev;
    bool first = true;
    for (Real rho = -0.999; rho <= 0.999; rho += 0.001) {
        HestonCFParams p = base;
        p.rho = rho;
        Complex phi = heston_characteristic_function(u, tau, S0, p);
        if (!first) {
            Real delta = std::abs(phi - phi_prev);
            // 步长 0.001 下, 相邻 cf 变化应 < 0.01 (平滑)
            EXPECT_LT(delta, 0.01) << "rho=" << rho;
        }
        phi_prev = phi;
        first = false;
    }
}

TEST(HestonCFEdgeCases, CFContinuousInSigmaVScan) {
    // σ_v 从 0.001 → 1.0 扫描, cf 应连续变化
    HestonCFParams base = schoutens_params();
    Real tau = 1.0;
    Real S0 = 100.0;
    Complex u(1.0, 0);
    Complex phi_prev;
    bool first = true;
    for (Real sigma_v = 0.001; sigma_v <= 1.0; sigma_v += 0.001) {
        HestonCFParams p = base;
        p.sigma = sigma_v;
        Complex phi = heston_characteristic_function(u, tau, S0, p);
        if (!first) {
            Real delta = std::abs(phi - phi_prev);
            EXPECT_LT(delta, 0.01) << "sigma_v=" << sigma_v;
        }
        phi_prev = phi;
        first = false;
    }
}

TEST(HestonCFEdgeCases, CFContinuousInTauScan) {
    // τ 从 1e-6 → 10.0 扫描, cf 应连续变化
    HestonCFParams p = schoutens_params();
    Real S0 = 100.0;
    Complex u(1.0, 0);
    Complex phi_prev;
    bool first = true;
    // 对数等间隔扫描
    for (Real tau = 1e-6; tau <= 10.0; tau *= 1.1) {
        Complex phi = heston_characteristic_function(u, tau, S0, p);
        if (!first) {
            Real delta = std::abs(phi - phi_prev);
            // 对数步长 1.1 下, 相邻 cf 变化应 < 0.1
            EXPECT_LT(delta, 0.1) << "tau=" << tau;
        }
        phi_prev = phi;
        first = false;
    }
}

// ============================================================================
// 7. 概率积分校验 (cf 基本性质)
// ============================================================================
// 特征函数基本性质: φ(0) = 1, |φ(u)| ≤ 1

TEST(HestonCFEdgeCases, CFAtZeroIsOneAllExtremeParams) {
    // 所有极端参数下 φ(0) = 1
    Real tau = 1.0;
    Real S0 = 100.0;
    std::vector<HestonCFParams> extreme_params = {
        {0.04, 1.5, 0.04, 0.3, 0.999, 0.0, 0.0},    // ρ → +1
        {0.04, 1.5, 0.04, 0.3, -0.999, 0.0, 0.0},   // ρ → -1
        {0.04, 1.5, 0.04, 1e-6, -0.5, 0.0, 0.0},    // σ_v → 0
        {0.04, 0.5, 0.04, 0.3, -0.7, 0.0, 0.0},     // Feller 违反
        {0.04, 1.5, 0.04, 0.3, -0.5, 0.0, 0.0}      // Schoutens 基准
    };
    for (const auto& p : extreme_params) {
        Complex phi = heston_characteristic_function(Complex(0, 0), tau, S0, p);
        EXPECT_NEAR(std::real(phi), 1.0, kTolRegular);
        EXPECT_NEAR(std::imag(phi), 0.0, kTolRegular);
    }
}

TEST(HestonCFEdgeCases, CFModulusLeqOneAllExtremeParams) {
    // 所有极端参数下 |φ(u)| ≤ 1 (概率测度下的特征函数性质)
    Real tau = 1.0;
    Real S0 = 100.0;
    std::vector<HestonCFParams> extreme_params = {
        {0.04, 1.5, 0.04, 0.3, 0.999, 0.0, 0.0},
        {0.04, 1.5, 0.04, 0.3, -0.999, 0.0, 0.0},
        {0.04, 1.5, 0.04, 1e-6, -0.5, 0.0, 0.0},
        {0.04, 0.5, 0.04, 0.3, -0.7, 0.0, 0.0}
    };
    for (const auto& p : extreme_params) {
        for (Real u_val = 0.01; u_val <= 50.0; u_val += 0.5) {
            Complex u(u_val, 0);
            Complex phi = heston_characteristic_function(u, tau, S0, p);
            EXPECT_LE(std::abs(phi), 1.0 + kTolExtreme)
                << "u=" << u_val << " rho=" << p.rho << " sigma_v=" << p.sigma;
        }
    }
}

// ============================================================================
// 8. 复平面虚部位移测试 (Carr-Madan P1 积分需要 u - i)
// ============================================================================
// Carr-Madan P1 积分使用 u_shifted = u - i (虚部位移)
// 分支切割问题常在此处出现 (原版 Heston trap)

TEST(HestonCFEdgeCases, ComplexUScanNoBranchCut) {
    // 复平面 u = u_real + i*u_imag 扫描, 验证无分支切割
    HestonCFParams p = schoutens_params();
    Real tau = 1.0;
    Real S0 = 100.0;
    // Carr-Madan P1 用 u - i, P2 用 u (实轴)
    for (Real u_real = 0.1; u_real <= 10.0; u_real += 0.5) {
        for (Real u_imag = -1.5; u_imag <= 1.5; u_imag += 0.5) {
            Complex u(u_real, u_imag);
            // 双侧差分连续性 (实部和虚部方向)
            Real eps = 1e-6;
            Complex phi_real_fwd = heston_characteristic_function(u + eps, tau, S0, p);
            Complex phi_real_bwd = heston_characteristic_function(u - eps, tau, S0, p);
            Complex phi_imag_fwd = heston_characteristic_function(u + Complex(0, eps), tau, S0, p);
            Complex phi_imag_bwd = heston_characteristic_function(u - Complex(0, eps), tau, S0, p);
            Complex phi_mid = heston_characteristic_function(u, tau, S0, p);

            Real resid_real = std::abs(phi_real_fwd + phi_real_bwd - Real(2) * phi_mid);
            Real resid_imag = std::abs(phi_imag_fwd + phi_imag_bwd - Real(2) * phi_mid);
            // 复平面扫描容差放宽到 1e-7 (数值微分误差 + 浮点精度)
            Real kTolComplex = 1e-7;

            EXPECT_LT(resid_real, kTolComplex)
                << "u_real=" << u_real << " u_imag=" << u_imag << " (real direction)";
            EXPECT_LT(resid_imag, kTolComplex)
                << "u_real=" << u_real << " u_imag=" << u_imag << " (imag direction)";
        }
    }
}

TEST(HestonCFEdgeCases, CarrMadanP1ImagShiftStable) {
    // Carr-Madan P1 积分核心: u_shifted = u - i (固定虚部 -1)
    // 在极端参数下验证此位移仍稳定
    std::vector<HestonCFParams> extreme_params = {
        {0.04, 1.5, 0.04, 0.3, 0.999, 0.0, 0.0},
        {0.04, 1.5, 0.04, 0.3, -0.999, 0.0, 0.0},
        {0.04, 1.5, 0.04, 1e-4, -0.5, 0.0, 0.0},
        {0.04, 0.5, 0.04, 0.3, -0.7, 0.0, 0.0}  // Feller 违反
    };
    Real tau = 1.0;
    Real S0 = 100.0;
    for (const auto& p : extreme_params) {
        for (Real u_real = 0.01; u_real <= 20.0; u_real += 0.1) {
            Complex u_shifted(u_real, -1.0);  // Carr-Madan P1
            Complex phi = heston_characteristic_function(u_shifted, tau, S0, p);
            EXPECT_TRUE(std::isfinite(std::real(phi)))
                << "u_real=" << u_real << " rho=" << p.rho;
            EXPECT_TRUE(std::isfinite(std::imag(phi)))
                << "u_real=" << u_real << " rho=" << p.rho;
        }
    }
}

// ============================================================================
// 9. 与 Schoutens 基准值交叉验证 (确保极端参数测试不破坏常规正确性)
// ============================================================================

TEST(HestonCFEdgeCases, SchoutensBenchmarkPreserved) {
    // 确保新增极端参数测试后, 常规 Schoutens 基准值仍正确
    HestonCFParams p = schoutens_params();
    Real tau = 1.0;
    Real S0 = 100.0;

    // Schoutens 表 u=0.5
    Complex phi_05 = heston_characteristic_function(Complex(0.5, 0), tau, S0, p);
    EXPECT_NEAR(std::real(phi_05), -0.6573742159702644, kTolRegular);
    EXPECT_NEAR(std::imag(phi_05), 0.7466039816575294, kTolRegular);

    // Schoutens 表 u=1.0
    Complex phi_10 = heston_characteristic_function(Complex(1.0, 0), tau, S0, p);
    EXPECT_NEAR(std::real(phi_10), -0.1231706804396436, kTolRegular);
    EXPECT_NEAR(std::imag(phi_10), -0.9715285746424274, kTolRegular);

    // Schoutens 表 u=2.0
    Complex phi_20 = heston_characteristic_function(Complex(2.0, 0), tau, S0, p);
    EXPECT_NEAR(std::real(phi_20), -0.8932404438080451, kTolRegular);
    EXPECT_NEAR(std::imag(phi_20), 0.2240614476780658, kTolRegular);
}

TEST(HestonCFEdgeCases, HestonVarianceCFEdgeCases) {
    // heston_variance_cf (方差特征函数) 在极端参数下也应稳定
    std::vector<HestonCFParams> extreme_params = {
        {0.04, 1.5, 0.04, 0.3, 0.999, 0.0, 0.0},
        {0.04, 1.5, 0.04, 0.3, -0.999, 0.0, 0.0},
        {0.04, 0.5, 0.04, 0.3, -0.7, 0.0, 0.0}  // Feller 违反
    };
    Real tau = 1.0;
    for (const auto& p : extreme_params) {
        for (Real u_val = 0.1; u_val <= 10.0; u_val += 0.5) {
            Complex u(u_val, 0);
            Complex phi = heston_variance_cf(u, tau, p);
            EXPECT_TRUE(std::isfinite(std::real(phi))) << "u=" << u_val;
            EXPECT_TRUE(std::isfinite(std::imag(phi))) << "u=" << u_val;
            // 方差 cf 的 |φ| 应 ≤ 1 (在方差测度下)
            EXPECT_LE(std::abs(phi), 1.0 + kTolExtreme) << "u=" << u_val;
        }
    }
}

// ============================================================================
// 10. 组合极端参数 (多重 stress)
// ============================================================================

TEST(HestonCFEdgeCases, CombinedExtremeParamsStable) {
    // 多重极端参数组合: ρ→+1 + 短到期 + 高波动率
    std::vector<HestonCFParams> combined_extreme = {
        {0.04, 1.5, 0.04, 0.3, 0.999, 0.0, 0.0},       // ρ→+1
        {0.04, 1.5, 0.04, 0.3, -0.999, 0.0, 0.0},      // ρ→-1
        {0.04, 0.5, 0.04, 0.5, 0.9, 0.0, 0.0},         // Feller 违反 + 高 ρ
        {0.01, 0.3, 0.01, 0.2, 0.95, 0.0, 0.0},        // 低方差 + 高 ρ + Feller 边界
        {0.09, 1.0, 0.09, 0.4, -0.95, 0.05, 0.02}      // 高方差 + 高 |ρ| + 含 r, q
    };
    Real S0 = 100.0;
    for (const auto& p : combined_extreme) {
        for (Real tau : {0.001, 0.01, 0.1, 1.0, 5.0}) {
            for (Real u_val = 0.01; u_val <= 10.0; u_val += 0.5) {
                Complex u(u_val, 0);
                Complex phi = heston_characteristic_function(u, tau, S0, p);
                EXPECT_TRUE(std::isfinite(std::real(phi)))
                    << "u=" << u_val << " tau=" << tau << " rho=" << p.rho;
                EXPECT_TRUE(std::isfinite(std::imag(phi)))
                    << "u=" << u_val << " tau=" << tau << " rho=" << p.rho;
                EXPECT_LE(std::abs(phi), 1.0 + kTolExtreme)
                    << "u=" << u_val << " tau=" << tau << " rho=" << p.rho;
            }
        }
    }
}
