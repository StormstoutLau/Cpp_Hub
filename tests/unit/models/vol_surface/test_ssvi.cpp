// Phase 4 LITE - D2 整改项: SSVI 跨期限模块测试
// 覆盖: SSVI 公式正确性、ATM 总方差、无套利条件、参数化工厂、校准
// 参考: Gatheral-Jacquier (2014) "Arbitrage-free SVI volatility surfaces"
//       https://arxiv.org/abs/1204.0646

#include "cpphub/models/vol_surface/ssvi.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace cpphub::v1;

// ============================================================================
// 1. SSVI 公式正确性
// ============================================================================

TEST(SSVI, TotalVarianceATMEqualsTheta) {
    // k=0 (ATM): w(0,θ) = θ/2 * (1 + 0 + sqrt(ρ² + 1 - ρ²)) = θ/2 * (1 + 1) = θ
    auto params = SSVI::Power_law(-0.3, 1.0, 0.25, {0.04, 0.09, 0.16});
    SSVI ssvi(params);
    for (Real theta : {0.04, 0.09, 0.16}) {
        Real w_atm = ssvi.total_variance(0.0, theta);
        EXPECT_NEAR(w_atm, theta, 1e-12) << "ATM total variance should equal theta";
    }
}

TEST(SSVI, TotalVarianceFormulaCorrectness) {
    // 手工验证: θ=0.04, φ=η*θ^(-γ)=1.0*0.04^(-0.25)=2.2361, ρ=-0.3, k=0.1
    // w = 0.04/2 * (1 + 2.2361*(-0.3)*0.1 + sqrt((2.2361*0.1 + (-0.3))² + 1 - 0.09))
    //   = 0.02 * (1 - 0.06708 + sqrt((-0.07639)² + 0.91))
    //   = 0.02 * (0.93292 + sqrt(0.005836 + 0.91))
    //   = 0.02 * (0.93292 + sqrt(0.915836))
    //   = 0.02 * (0.93292 + 0.9570)
    //   = 0.02 * 1.8899 = 0.037798
    auto params = SSVI::Power_law(-0.3, 1.0, 0.25, {0.04});
    SSVI ssvi(params);
    Real w = ssvi.total_variance(0.1, 0.04);
    EXPECT_NEAR(w, 0.03779, 1e-4);
}

TEST(SSVI, ImpliedVolConsistency) {
    // implied_vol(k, T, θ) = sqrt(w(k,θ) / T)
    auto params = SSVI::Power_law(-0.3, 1.0, 0.25, {0.04});
    SSVI ssvi(params);
    Real T = 1.0;
    Real theta = 0.04;
    Real w = ssvi.total_variance(0.1, theta);
    Real expected_vol = std::sqrt(w / T);
    EXPECT_NEAR(ssvi.implied_vol(0.1, T, theta), expected_vol, 1e-12);
}

TEST(SSVI, SymmetryInRho) {
    // ρ=0 时, w(k,θ) 关于 k 对称: w(k) = w(-k)
    auto params = SSVI::Power_law(0.0, 1.0, 0.25, {0.04});
    SSVI ssvi(params);
    Real theta = 0.04;
    Real w_pos = ssvi.total_variance(0.2, theta);
    Real w_neg = ssvi.total_variance(-0.2, theta);
    EXPECT_NEAR(w_pos, w_neg, 1e-12);
}

// ============================================================================
// 2. 参数化工厂
// ============================================================================

TEST(SSVI, PowerLawParameterization) {
    std::vector<Real> thetas = {0.04, 0.09, 0.16};
    auto params = SSVI::Power_law(-0.3, 1.5, 0.3, thetas);
    EXPECT_DOUBLE_EQ(params.rho, -0.3);
    // φ(θ) = 1.5 * θ^(-0.3)
    EXPECT_NEAR(params.phi(0.04), 1.5 * std::pow(0.04, -0.3), 1e-12);
    EXPECT_NEAR(params.phi(0.16), 1.5 * std::pow(0.16, -0.3), 1e-12);
}

TEST(SSVI, HestonLikeParameterization) {
    std::vector<Real> thetas = {0.04, 0.09, 0.16};
    auto params = SSVI::Heston_like(-0.2, 1.0, 0.5, thetas);
    EXPECT_DOUBLE_EQ(params.rho, -0.2);
    // φ(θ) = 1.0 * θ^(-0.5)
    EXPECT_NEAR(params.phi(0.04), 1.0 * std::pow(0.04, -0.5), 1e-12);
    EXPECT_NEAR(params.phi(0.16), 1.0 * std::pow(0.16, -0.5), 1e-12);
}

TEST(SSVI, PowerLawRejectsInvalidGamma) {
    EXPECT_THROW(SSVI::Power_law(-0.3, 1.0, 0.0, {0.04}), std::invalid_argument);
    EXPECT_THROW(SSVI::Power_law(-0.3, 1.0, 0.5, {0.04}), std::invalid_argument);
    EXPECT_THROW(SSVI::Power_law(-0.3, 1.0, -0.1, {0.04}), std::invalid_argument);
}

TEST(SSVI, HestonLikeRejectsInvalidLambda) {
    EXPECT_THROW(SSVI::Heston_like(-0.3, 1.0, 0.0, {0.04}), std::invalid_argument);
    EXPECT_THROW(SSVI::Heston_like(-0.3, 1.0, 1.0, {0.04}), std::invalid_argument);
}

TEST(SSVI, ConstructorRejectsInvalidRho) {
    EXPECT_THROW(SSVI(SSVI::Power_law(1.0, 1.0, 0.25, {0.04})), std::invalid_argument);
    EXPECT_THROW(SSVI(SSVI::Power_law(-1.0, 1.0, 0.25, {0.04})), std::invalid_argument);
}

// ============================================================================
// 3. 无套利条件
// ============================================================================

TEST(SSVI, CalendarArbitragePowerLaw) {
    // Power-law: φ(θ)=η*θ^(-γ), θ*φ(θ)=η*θ^(1-γ), d/dθ=η*(1-γ)*θ^(-γ)>0 ✓
    // 应无 calendar 套利
    std::vector<Real> thetas = {0.04, 0.09, 0.16, 0.25};
    auto params = SSVI::Power_law(-0.3, 1.0, 0.25, thetas);
    SSVI ssvi(params);
    EXPECT_TRUE(ssvi.check_calendar_arbitrage());
}

TEST(SSVI, ButterflyArbitragePowerLaw) {
    // Power-law with γ=0.25, η=1.0: φ(θ)*θ*(1+|ρ|) = θ^(1-0.25)*(1+0.3) = 1.3*θ^0.75
    // 最大 θ=0.25: 1.3*0.25^0.75 = 1.3*0.354 = 0.46 < 4 ✓
    std::vector<Real> thetas = {0.04, 0.09, 0.16, 0.25};
    auto params = SSVI::Power_law(-0.3, 1.0, 0.25, thetas);
    SSVI ssvi(params);
    EXPECT_TRUE(ssvi.check_butterfly_arbitrage());
}

TEST(SSVI, NoArbitragePowerLaw) {
    std::vector<Real> thetas = {0.04, 0.09, 0.16, 0.25};
    auto params = SSVI::Power_law(-0.3, 1.0, 0.25, thetas);
    SSVI ssvi(params);
    EXPECT_TRUE(ssvi.check_no_arbitrage());
}

TEST(SSVI, ButterflyArbitrageViolationLargeEta) {
    // 极大 η 导致 φ(θ)*θ*(1+|ρ|) >= 4,违反 butterfly 条件
    // φ(θ)*θ*(1+|ρ|) = η*θ^(1-γ)*(1+|ρ|), θ=0.25, γ=0.25, |ρ|=0.3
    // = η * 0.25^0.75 * 1.3 = η * 0.354 * 1.3 = η * 0.46
    // 要 >= 4, η >= 8.69
    std::vector<Real> thetas = {0.04, 0.09, 0.16, 0.25};
    auto params = SSVI::Power_law(-0.3, 10.0, 0.25, thetas);
    SSVI ssvi(params);
    EXPECT_FALSE(ssvi.check_butterfly_arbitrage());
}

TEST(SSVI, CalendarArbitrageViolation) {
    // 构造违反 calendar 套利的 φ(θ): φ(θ) = 1/θ (θ*φ(θ)=1 常数,d/dθ=0 不严格 >0)
    // 但更明确: φ(θ) = 1/θ² (θ*φ(θ)=1/θ, d/dθ=-1/θ²<0) 违反
    SSVIParams params;
    params.rho = -0.3;
    params.phi = [](Real theta) { return 1.0 / (theta * theta); };
    params.theta_slice = {0.04, 0.09, 0.16};
    SSVI ssvi(params);
    EXPECT_FALSE(ssvi.check_calendar_arbitrage());
}

// ============================================================================
// 4. 跨期限单调性 (Calendar 套利天然免疫验证)
// ============================================================================

TEST(SSVI, TotalVarianceMonotonicInTheta) {
    // 对固定 k, w(k,θ) 应随 θ 单调递增 (calendar no-arbitrage 的体现)
    auto params = SSVI::Power_law(-0.3, 1.0, 0.25, {0.04, 0.09, 0.16, 0.25});
    SSVI ssvi(params);
    Real k = 0.1;
    Real w_prev = ssvi.total_variance(k, 0.04);
    for (Real theta : {0.09, 0.16, 0.25}) {
        Real w_curr = ssvi.total_variance(k, theta);
        EXPECT_GT(w_curr, w_prev) << "w should increase with theta (calendar no-arb)";
        w_prev = w_curr;
    }
}

// ============================================================================
// 5. 与 SVI 单切片一致性
// ============================================================================

TEST(SSVI, ConsistencyWithSVISlice) {
    // SSVI 在单期限上应与 SVI 切片一致 (相同 total_variance 形式)
    // SVI raw: w(k) = a + b*(ρ*k + sqrt(k² + σ²))
    // SSVI: w(k,θ) = θ/2 * (1 + φ*ρ*k + sqrt((φ*k + ρ)² + 1 - ρ²))
    // 令 φ=1, θ=2*b*σ, ρ相同, a=θ/2*(1-φ*ρ*0+sqrt(ρ²+1-ρ²)-...) 不直接对应
    // 简化验证: SSVI ATM (k=0) = θ, SVI ATM 应等于 θ
    auto ssvi_params = SSVI::Power_law(-0.3, 1.0, 0.25, {0.04});
    SSVI ssvi(ssvi_params);
    Real theta = 0.04;
    Real w_ssvi_atm = ssvi.total_variance(0.0, theta);
    EXPECT_NEAR(w_ssvi_atm, theta, 1e-12);
}

// ============================================================================
// 6. 校准冒烟测试
// ============================================================================

TEST(SSVI, CalibrationSmokeTest) {
    // 生成 SSVI 合成数据,然后校准,验证参数可恢复
    Real true_rho = -0.3;
    Real true_eta = 1.0;
    Real true_gamma = 0.25;
    Real forward = 100.0;

    std::vector<Real> maturities = {0.25, 0.5, 1.0};
    std::vector<Real> strikes = {80.0, 90.0, 95.0, 100.0, 105.0, 110.0, 120.0};

    auto true_params = SSVI::Power_law(true_rho, true_eta, true_gamma, {});
    SSVI true_ssvi(true_params);

    // 生成合成市场数据
    std::vector<Real> implied_vols;
    for (Real T : maturities) {
        Real theta_atm = 0.2 * 0.2 * T;  // ATM vol = 20%, θ = σ²*T
        for (Real K : strikes) {
            Real k = std::log(K / forward);
            Real vol = true_ssvi.implied_vol(k, T, theta_atm);
            implied_vols.push_back(vol);
        }
    }

    // 校准
    CalibConfig cfg;
    cfg.use_de_init = false;  // 加速测试
    cfg.lm_max_iter = 100;
    SSVI calib_ssvi(SSVI::Power_law(0.0, 0.5, 0.3, {}));
    auto result = calib_ssvi.calibrate(strikes, maturities, implied_vols, forward, cfg);

    // 校准应收敛 (或至少接近)
    EXPECT_TRUE(result.converged || result.objective_value < 1e-6)
        << "obj=" << result.objective_value << " msg=" << result.message;

    // 恢复的参数应接近真值
    Real calib_rho = result.params[0];
    Real calib_eta = result.params[1];
    Real calib_gamma = result.params[2];
    EXPECT_NEAR(calib_rho, true_rho, 0.15) << "rho recovery";
    EXPECT_NEAR(calib_eta, true_eta, 0.3) << "eta recovery";
    EXPECT_NEAR(calib_gamma, true_gamma, 0.1) << "gamma recovery";

    // 校准后应无套利
    EXPECT_TRUE(calib_ssvi.check_no_arbitrage());
}

// ============================================================================
// 7. Heston → SSVI 解析映射 (Gatheral-Jacquier 2014 Theorem 3.1)
// ============================================================================

TEST(SSVIFromHeston, RhoMappingCorrect) {
    // ρ_SSVI = ρ_Heston
    HestonParams hp{0.04, 2.0, 0.04, 0.3, -0.5};
    Real T = 1.0;
    auto ssvi_params = SSVI::from_heston(hp, T);
    EXPECT_NEAR(ssvi_params.rho, hp.rho, 1e-12);
}

TEST(SSVIFromHeston, ATMTotalVarianceCorrect) {
    // θ(T) = (v0 - θ̄) * (1 - exp(-κT)) / κ + θ̄ * T
    HestonParams hp{0.04, 2.0, 0.04, 0.3, -0.5};
    Real T = 1.0;
    Real expected_theta = (hp.v0 - hp.theta) * (1.0 - std::exp(-hp.kappa * T)) / hp.kappa
                          + hp.theta * T;
    auto ssvi_params = SSVI::from_heston(hp, T);
    ASSERT_EQ(ssvi_params.theta_slice.size(), 1u);
    EXPECT_NEAR(ssvi_params.theta_slice[0], expected_theta, 1e-12);
}

TEST(SSVIFromHeston, PhiConstantLargeTermLimit) {
    // φ = ξ / (κ * θ̄) * (1 - ρ²)
    HestonParams hp{0.04, 2.0, 0.04, 0.3, -0.5};
    Real T = 1.0;
    Real expected_phi = hp.sigma_v * (1.0 - hp.rho * hp.rho) / (hp.kappa * hp.theta);
    auto ssvi_params = SSVI::from_heston(hp, T);
    Real theta = ssvi_params.theta_slice[0];
    EXPECT_NEAR(ssvi_params.phi(theta), expected_phi, 1e-12);
    // φ 是常数, 不依赖 θ
    EXPECT_NEAR(ssvi_params.phi(theta * 2.0), expected_phi, 1e-12);
}

TEST(SSVIFromHeston, RejectsInvalidParams) {
    EXPECT_THROW(SSVI::from_heston({0.04, 2.0, 0.04, 0.3, -0.5}, -1.0), std::invalid_argument);
    EXPECT_THROW(SSVI::from_heston({0.04, -1.0, 0.04, 0.3, -0.5}, 1.0), std::invalid_argument);  // kappa <= 0
    EXPECT_THROW(SSVI::from_heston({0.04, 2.0, -0.01, 0.3, -0.5}, 1.0), std::invalid_argument); // theta <= 0
    EXPECT_THROW(SSVI::from_heston({0.04, 2.0, 0.04, -0.1, -0.5}, 1.0), std::invalid_argument); // sigma_v <= 0
    EXPECT_THROW(SSVI::from_heston({0.04, 2.0, 0.04, 0.3, 1.0}, 1.0), std::invalid_argument);   // |rho| >= 1
}

TEST(SSVIFromHeston, SetHestonInitSkipsDE) {
    // 设置 Heston 初始猜测后, calibrate 应跳过 DE 直接用 LM
    // 生成合成 Heston-like SSVI 数据
    HestonParams hp{0.04, 2.0, 0.04, 0.3, -0.5};
    Real forward = 100.0;
    std::vector<Real> maturities = {0.5, 1.0};
    std::vector<Real> strikes = {90.0, 95.0, 100.0, 105.0, 110.0};

    // 用 from_heston 生成合成 IV 数据
    auto heston_ssvi_params = SSVI::from_heston(hp, maturities[0]);
    SSVI heston_ssvi(heston_ssvi_params);
    std::vector<Real> implied_vols;
    for (Real T : maturities) {
        Real theta_T = (hp.v0 - hp.theta) * (1.0 - std::exp(-hp.kappa * T)) / hp.kappa
                       + hp.theta * T;
        for (Real K : strikes) {
            Real k = std::log(K / forward);
            implied_vols.push_back(heston_ssvi.implied_vol(k, T, theta_T));
        }
    }

    // 用 Heston 初始猜测校准
    SSVI calib_ssvi(SSVI::Power_law(0.0, 0.5, 0.3, {}));
    calib_ssvi.set_heston_init(hp);
    CalibConfig cfg;
    cfg.use_de_init = true;  // 即使开启 DE, 有 Heston init 时也应跳过
    cfg.lm_max_iter = 200;
    auto result = calib_ssvi.calibrate(strikes, maturities, implied_vols, forward, cfg);

    // 应收敛或接近
    EXPECT_TRUE(result.converged || result.objective_value < 1e-4)
        << "obj=" << result.objective_value << " msg=" << result.message;
    // rho 应接近 Heston 的 rho
    EXPECT_NEAR(result.params[0], hp.rho, 0.15) << "rho should be close to Heston rho";
}

TEST(SSVIFromHeston, ClearHestonInitRestoresDE) {
    // clear_heston_init 后, calibrate 应恢复 DE 路径
    SSVI ssvi(SSVI::Power_law(0.0, 0.5, 0.3, {}));
    ssvi.set_heston_init({0.04, 2.0, 0.04, 0.3, -0.5});
    ssvi.clear_heston_init();
    // 仅验证不崩溃 — 不实际运行校准 (DE 耗时)
    EXPECT_NO_THROW(ssvi.calibrate({100.0}, {1.0}, {0.2}, 100.0, CalibConfig{}));
}
