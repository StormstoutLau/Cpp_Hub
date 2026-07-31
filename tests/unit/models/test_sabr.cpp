// SABR 模型单元测试
// 覆盖: Hagan IV 闭式公式 + SABR 过程模拟
//
// 测试分组:
//   1. 参数验证 (8 cases)
//   2. Hagan IV 数值正确性 (5 cases)
//   3. Hagan IV 极限与退化 (5 cases)
//   4. Hagan IV 微笑性质 (5 cases)
//   5. ATM IV 反推 alpha (3 cases)
//   6. SABR 过程模拟 (7 cases)
// 总计: ~33 cases

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include "cpphub/pricing/analytic/sabr_hagan.hpp"
#include "cpphub/models/diffusion/sabr.hpp"
#include "cpphub/core/rng.hpp"

using namespace cpphub;

// ==================== 1. 参数验证 ====================

TEST(SABRParamsTest, NegativeAlphaThrows) {
    SABRParams p{-0.1, 0.5, -0.3, 0.4};
    EXPECT_THROW(validate_sabr_params(p), std::invalid_argument);
}

TEST(SABRParamsTest, ZeroAlphaThrows) {
    SABRParams p{0.0, 0.5, -0.3, 0.4};
    EXPECT_THROW(validate_sabr_params(p), std::invalid_argument);
}

TEST(SABRParamsTest, BetaOutOfRangeThrows) {
    SABRParams p1{0.2, -0.1, -0.3, 0.4};
    SABRParams p2{0.2, 1.1, -0.3, 0.4};
    EXPECT_THROW(validate_sabr_params(p1), std::invalid_argument);
    EXPECT_THROW(validate_sabr_params(p2), std::invalid_argument);
}

TEST(SABRParamsTest, RhoBoundaryThrows) {
    SABRParams p1{0.2, 0.5, -1.0, 0.4};
    SABRParams p2{0.2, 0.5, 1.0, 0.4};
    EXPECT_THROW(validate_sabr_params(p1), std::invalid_argument);
    EXPECT_THROW(validate_sabr_params(p2), std::invalid_argument);
}

TEST(SABRParamsTest, NegativeNuThrows) {
    SABRParams p{0.2, 0.5, -0.3, -0.1};
    EXPECT_THROW(validate_sabr_params(p), std::invalid_argument);
}

TEST(SABRParamsTest, ValidParamsNoThrow) {
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    EXPECT_NO_THROW(validate_sabr_params(p));
}

TEST(SABRParamsTest, BetaBoundaryValuesValid) {
    SABRParams p0{0.2, 0.0, -0.3, 0.4};   // β = 0 (normal SABR)
    SABRParams p1{0.2, 1.0, -0.3, 0.4};   // β = 1 (lognormal SABR)
    EXPECT_NO_THROW(validate_sabr_params(p0));
    EXPECT_NO_THROW(validate_sabr_params(p1));
}

TEST(SABRParamsTest, ZeroNuValid) {
    SABRParams p{0.2, 0.5, -0.3, 0.0};    // ν = 0 (CEV 局部波动率极限)
    EXPECT_NO_THROW(validate_sabr_params(p));
}

// ==================== 2. Hagan IV 数值正确性 ====================

TEST(SABRHaganTest, NonNegativeStrikeThrows) {
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    EXPECT_THROW(sabr_implied_vol_hagan(-0.01, 0.05, 1.0, p), std::invalid_argument);
    EXPECT_THROW(sabr_implied_vol_hagan(0.0, 0.05, 1.0, p), std::invalid_argument);
}

TEST(SABRHaganTest, NonPositiveForwardThrows) {
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    EXPECT_THROW(sabr_implied_vol_hagan(0.05, 0.0, 1.0, p), std::invalid_argument);
    EXPECT_THROW(sabr_implied_vol_hagan(0.05, -0.01, 1.0, p), std::invalid_argument);
}

TEST(SABRHaganTest, NegativeTimeThrows) {
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    EXPECT_THROW(sabr_implied_vol_hagan(0.05, 0.05, -0.1, p), std::invalid_argument);
}

TEST(SABRHaganTest, ATMImpliedVolPositiveAndFinite) {
    // 典型利率衍生品参数: α=0.2, β=0.5, ρ=-0.3, ν=0.4, F=0.05, T=1
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    Real F = 0.05, T = 1.0;
    Real iv = sabr_implied_vol_hagan(F, F, T, p);
    EXPECT_GT(iv, 0.0);
    EXPECT_LT(iv, 1.0);  // 波动率 < 100%
    EXPECT_TRUE(std::isfinite(iv));
}

TEST(SABRHaganTest, OTMImpliedVolPositiveAndFinite) {
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    Real F = 0.05, K = 0.06, T = 2.0;  // OTM call
    Real iv = sabr_implied_vol_hagan(K, F, T, p);
    EXPECT_GT(iv, 0.0);
    EXPECT_LT(iv, 2.0);
    EXPECT_TRUE(std::isfinite(iv));
}

TEST(SABRHaganTest, ITMImpliedVolPositiveAndFinite) {
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    Real F = 0.05, K = 0.04, T = 2.0;  // ITM call
    Real iv = sabr_implied_vol_hagan(K, F, T, p);
    EXPECT_GT(iv, 0.0);
    EXPECT_LT(iv, 2.0);
    EXPECT_TRUE(std::isfinite(iv));
}

TEST(SABRHaganTest, ATMViaGenericFormulaMatchesATMSpecialized) {
    // K=F 时, 完整公式 (通过分支判断) 应与 ATM 闭式公式一致
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    Real F = 0.05, T = 1.0;
    Real iv_general = sabr_implied_vol_hagan(F, F, T, p);
    Real iv_atm = sabr_implied_vol_atm(F, T, p);
    EXPECT_NEAR(iv_general, iv_atm, 1e-12);
}

// ==================== 3. Hagan IV 极限与退化 ====================

TEST(SABRHaganTest, ZeroTimeToMaturityReturnsInstantaneousVol) {
    // T → 0: σ_B → α / F^(1-β)
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    Real F = 0.05;
    Real iv = sabr_implied_vol_hagan(F, F, 0.0, p);
    Real expected = p.alpha / std::pow(F, 1.0 - p.beta);
    EXPECT_NEAR(iv, expected, 1e-12);
}

TEST(SABRHaganTest, ZeroNuReturnsCEVVolatility) {
    // ν = 0: σ 不随机, SABR 退化为 CEV, IV = α / F^(1-β) (对所有 K)
    // 实际上 ν=0 时 z=0, x(z)=0, 需用 0/0 极限 → 1, IV = α / D · (1 + c·T)
    // 但当 β=1 时 D=1, c 仅含 ν² 项 (为 0), 故 IV = α (常数)
    SABRParams p{0.2, 1.0, -0.3, 0.0};  // β=1, ν=0
    Real F = 100.0, T = 1.0;
    Real iv_atm = sabr_implied_vol_atm(F, T, p);
    // β=1: σ_ATM = α/1 · (1 + 0 + 0 + 0) = α
    EXPECT_NEAR(iv_atm, 0.2, 1e-12);
}

TEST(SABRHaganTest, BetaOneGivesLognormalBehavior) {
    // β=1: (FK)^((1-β)/2) = 1, IV 不依赖 F·K 乘积的幂
    SABRParams p{0.2, 1.0, -0.0, 0.3};  // ρ=0 对称微笑
    Real F = 100.0, T = 1.0;
    Real K_otm = 110.0;
    Real K_itm = 100.0 * 100.0 / K_otm;  // F²/K_otm, 保证 ln(F/K) 对称
    Real iv_otm = sabr_implied_vol_hagan(K_otm, F, T, p);
    Real iv_itm = sabr_implied_vol_hagan(K_itm, F, T, p);
    // β=1, ρ=0 时 IV 对 K 对称 (关于 F)
    EXPECT_NEAR(iv_otm, iv_itm, 1e-9);
}

TEST(SABRHaganTest, BetaZeroNormalSABR) {
    // β=0: normal SABR, (FK)^((1-β)/2) = √(FK)
    SABRParams p{0.005, 0.0, -0.3, 0.4};  // 小 alpha, 因为 normal vol 量级不同
    Real F = 0.03, T = 1.0;
    Real iv_atm = sabr_implied_vol_atm(F, T, p);
    // β=0: σ_ATM = α/F · (1 + (α²/F²)/24 · T + 0 + ((2-3ρ²)/24)ν²T)
    Real expected = (p.alpha / F) * (1.0
                   + (p.alpha * p.alpha / (F * F) / 24.0) * T
                   + ((2.0 - 3.0 * p.rho * p.rho) / 24.0) * p.nu * p.nu * T);
    EXPECT_NEAR(iv_atm, expected, 1e-12);
}

TEST(SABRHaganTest, NearATMFormulaContinuous) {
    // K 接近 F 时, 完整公式与 ATM 公式连续过渡
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    Real F = 0.05, T = 1.0;
    Real iv_atm = sabr_implied_vol_atm(F, T, p);

    Real eps = 1e-6;
    Real iv_near1 = sabr_implied_vol_hagan(F * (1 - eps), F, T, p);
    Real iv_near2 = sabr_implied_vol_hagan(F * (1 + eps), F, T, p);
    EXPECT_NEAR(iv_near1, iv_atm, 1e-6);
    EXPECT_NEAR(iv_near2, iv_atm, 1e-6);
}

// ==================== 4. Hagan IV 微笑性质 ====================

TEST(SABRHaganTest, SmileIsUShaped) {
    // 微笑: OTM 和 ITM 的 IV 应高于 ATM 的 IV (典型微笑形状)
    // 注: 仅当 β=1 (lognormal SABR, (FK)^0=1) 且 ρ=0 时严格对称微笑
    // β<1 时 (FK)^((1-β)/2) 项随 K 变化, 会引入单调倾斜掩盖微笑
    SABRParams p{0.2, 1.0, 0.0, 0.4};  // β=1, ρ=0 严格对称微笑
    Real F = 100.0, T = 1.0;
    Real K_otm = 110.0;
    Real K_itm = F * F / K_otm;  // 对数对称: ln(F/K_otm) = -ln(F/K_itm)
    Real iv_atm = sabr_implied_vol_atm(F, T, p);
    Real iv_otm = sabr_implied_vol_hagan(K_otm, F, T, p);
    Real iv_itm = sabr_implied_vol_hagan(K_itm, F, T, p);
    EXPECT_GT(iv_otm, iv_atm);
    EXPECT_GT(iv_itm, iv_atm);
    // β=1, ρ=0 且 K 对数对称时 IV 严格对称
    EXPECT_NEAR(iv_otm, iv_itm, 1e-12);
}

TEST(SABRHaganTest, NegativeRhoTiltsSmileDownwardInK) {
    // ρ < 0: K 越大 IV 越小 (skew, 股票市场典型形态)
    // 用 β=1 避免 (FK)^((1-β)/2) 项干扰, 单独检验 ρ 的倾斜效果
    SABRParams p{0.2, 1.0, -0.5, 0.4};
    Real F = 100.0, T = 1.0;
    Real iv_low_K = sabr_implied_vol_hagan(90.0, F, T, p);
    Real iv_high_K = sabr_implied_vol_hagan(110.0, F, T, p);
    EXPECT_GT(iv_low_K, iv_high_K);
}

TEST(SABRHaganTest, PositiveRhoTiltsSmileUpwardInK) {
    // ρ > 0: K 越大 IV 越大 (反向 skew, 商品市场常见)
    // 用 β=1 避免 (FK)^((1-β)/2) 项干扰
    SABRParams p{0.2, 1.0, 0.5, 0.4};
    Real F = 100.0, T = 1.0;
    Real iv_low_K = sabr_implied_vol_hagan(90.0, F, T, p);
    Real iv_high_K = sabr_implied_vol_hagan(110.0, F, T, p);
    EXPECT_LT(iv_low_K, iv_high_K);
}

TEST(SABRHaganTest, HigherNuGivesMoreCurvedSmile) {
    // ν 越大, 微笑越陡 (OTM IV 相对 ATM 上升更多)
    Real F = 0.05, T = 1.0;
    Real K_otm = 0.06;
    SABRParams p_low{0.2, 0.5, 0.0, 0.2};
    SABRParams p_high{0.2, 0.5, 0.0, 0.8};
    Real iv_atm_low = sabr_implied_vol_atm(F, T, p_low);
    Real iv_atm_high = sabr_implied_vol_atm(F, T, p_high);
    Real iv_otm_low = sabr_implied_vol_hagan(K_otm, F, T, p_low);
    Real iv_otm_high = sabr_implied_vol_hagan(K_otm, F, T, p_high);
    // OTM IV - ATM IV 的差应随 ν 增大而增大
    Real skew_low = iv_otm_low - iv_atm_low;
    Real skew_high = iv_otm_high - iv_atm_high;
    EXPECT_GT(skew_high, skew_low);
}

TEST(SABRHaganTest, AlphaScalesATMVolProportionally) {
    // σ_ATM ≈ α / F^(1-β) (主项), α 增大 → σ_ATM 增大
    Real F = 0.05, T = 1.0;
    SABRParams p1{0.1, 0.5, -0.3, 0.4};
    SABRParams p2{0.3, 0.5, -0.3, 0.4};
    Real iv1 = sabr_implied_vol_atm(F, T, p1);
    Real iv2 = sabr_implied_vol_atm(F, T, p2);
    EXPECT_GT(iv2, iv1);
    // 主项比例: iv2/iv1 ≈ α2/α1 = 3 (因高阶项略有偏差)
    EXPECT_NEAR(iv2 / iv1, 3.0, 0.05);
}

// ==================== 5. ATM IV 反推 alpha ====================

TEST(SABRSolveAlphaTest, RoundTripRecoversAlpha) {
    // 给定 α 计算 σ_ATM, 再从 σ_ATM 反推 α, 应一致
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    Real F = 0.05, T = 1.0;
    Real sigma_atm = sabr_implied_vol_atm(F, T, p);
    Real alpha_recovered = sabr_solve_alpha_from_atm(sigma_atm, F, T, p.beta, p.rho, p.nu);
    EXPECT_NEAR(alpha_recovered, p.alpha, 1e-10);
}

TEST(SABRSolveAlphaTest, RoundTripAtBetaOne) {
    // β=1 时反推
    SABRParams p{0.25, 1.0, -0.2, 0.3};
    Real F = 100.0, T = 2.0;
    Real sigma_atm = sabr_implied_vol_atm(F, T, p);
    Real alpha_recovered = sabr_solve_alpha_from_atm(sigma_atm, F, T, p.beta, p.rho, p.nu);
    EXPECT_NEAR(alpha_recovered, p.alpha, 1e-10);
}

TEST(SABRSolveAlphaTest, RoundTripAtBetaZero) {
    // β=0 时反推
    SABRParams p{0.005, 0.0, -0.3, 0.5};
    Real F = 0.03, T = 1.0;
    Real sigma_atm = sabr_implied_vol_atm(F, T, p);
    Real alpha_recovered = sabr_solve_alpha_from_atm(sigma_atm, F, T, p.beta, p.rho, p.nu);
    EXPECT_NEAR(alpha_recovered, p.alpha, 1e-10);
}

// ==================== 6. SABR 过程模拟 ====================

TEST(SABRProcessTest, ConstructorValidation) {
    SABRParams valid{0.2, 0.5, -0.3, 0.4};
    EXPECT_NO_THROW(SABRProcess(valid, 0.05));

    SABRParams bad_alpha{-0.1, 0.5, -0.3, 0.4};
    EXPECT_THROW(SABRProcess(bad_alpha, 0.05), std::invalid_argument);

    SABRParams bad_beta{0.2, 1.5, -0.3, 0.4};
    EXPECT_THROW(SABRProcess(bad_beta, 0.05), std::invalid_argument);

    EXPECT_THROW(SABRProcess(valid, -0.01), std::invalid_argument);  // F0 <= 0
    EXPECT_THROW(SABRProcess(valid, 0.0), std::invalid_argument);
}

TEST(SABRProcessTest, DimensionIsTwo) {
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    SABRProcess sabr(p, 0.05);
    EXPECT_EQ(sabr.dimension(), 2u);
}

TEST(SABRProcessTest, SpotReturnsForward0) {
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    SABRProcess sabr(p, 0.05);
    EXPECT_DOUBLE_EQ(sabr.spot(), 0.05);
    EXPECT_DOUBLE_EQ(sabr.forward0(), 0.05);
}

TEST(SABRProcessTest, PathLengthAndInitialValue) {
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    SABRProcess sabr(p, 0.05);
    Size n_steps = 100;
    std::vector<Real> path(n_steps + 1);
    Philox4x64 rng(42);
    sabr.generate_path(1.0, n_steps, path, rng);
    EXPECT_DOUBLE_EQ(path[0], 0.05);
    EXPECT_EQ(path.size(), n_steps + 1);
}

TEST(SABRProcessTest, PathNonNegative) {
    // Euler 吸收壁方案: F(t) >= 0
    SABRParams p{0.3, 0.5, -0.5, 0.5};
    SABRProcess sabr(p, 0.05);
    Size n_steps = 500;
    std::vector<Real> path(n_steps + 1);
    Philox4x64 rng(42);
    sabr.generate_path(2.0, n_steps, path, rng);
    for (Real x : path) {
        EXPECT_GE(x, 0.0);
    }
}

TEST(SABRProcessTest, BetaOneMeanMatchesLognormal) {
    // β=1, ν=0: SABR 退化为 GBM (σ = α 常数)
    // dF = α F dW_1
    // E[F(T)] = F_0 (drift-free, 因 SABR 无 r)
    // 但 MC 误差较大, 容差放宽
    SABRParams p{0.2, 1.0, 0.0, 0.0};  // β=1, ν=0 → GBM with σ=α
    Real F0 = 100.0;
    SABRProcess sabr(p, F0);
    Size n_paths = 20000;
    Size n_steps = 50;
    Real T = 1.0;
    Real alpha = p.alpha;

    Real sum_F = 0.0;
    Real sum_F2 = 0.0;
    Philox4x64 rng(12345);
    std::vector<Real> path(n_steps + 1);

    for (Size i = 0; i < n_paths; ++i) {
        sabr.generate_path(T, n_steps, path, rng);
        Real F_T = path.back();
        sum_F += F_T;
        sum_F2 += F_T * F_T;
    }

    Real mean_F = sum_F / n_paths;
    // GBM (无 drift): E[F(T)] = F_0
    EXPECT_NEAR(mean_F, F0, 0.02 * F0);  // 2% 容差 (MC 误差)

    // Var[F(T)] = F_0² (exp(σ²T) - 1), σ = α
    Real var_F = sum_F2 / n_paths - mean_F * mean_F;
    Real expected_var = F0 * F0 * (std::exp(alpha * alpha * T) - 1.0);
    EXPECT_NEAR(var_F, expected_var, 0.10 * expected_var);  // 10% 容差
}

TEST(SABRProcessTest, NegativeRhoGivesNegativeCorrelation) {
    // ρ < 0: F 上升时 σ 倾向下降, 反之亦然
    // 检验: 路径上 ΔF 和 Δσ 的样本相关性符号
    SABRParams p{0.2, 0.5, -0.8, 0.4};
    Real F0 = 100.0;
    SABRProcess sabr(p, F0);
    Size n_steps = 200;
    Real T = 1.0;
    Real dt = T / n_steps;

    std::vector<Real> dF_arr, dsigma_arr;
    dF_arr.reserve(n_steps * 100);
    dsigma_arr.reserve(n_steps * 100);

    Philox4x64 rng(42);
    std::vector<Real> path(n_steps + 1);

    Size n_paths = 100;
    for (Size i = 0; i < n_paths; ++i) {
        sabr.generate_path(T, n_steps, path, rng);
        // 需要单独跟踪 sigma, generate_path 不输出 sigma 路径
        // 重新模拟一次并记录 sigma
        Real F = F0;
        Real sigma = p.alpha;
        for (Size k = 0; k < n_steps; ++k) {
            Real Z1, Z2;
            // 复用 rng 生成相同路径 (但 rng 已被 generate_path 消耗, 这里用独立 rng)
            uint64_t r1 = rng();
            uint64_t r2 = rng();
            double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            auto [w1, w2] = box_muller(u1, u2);
            Z1 = w1;
            Z2 = p.rho * w1 + std::sqrt(1.0 - p.rho * p.rho) * w2;
            Real F_old = F;
            Real sigma_old = sigma;
            F = sabr.evolve(F, sigma, dt, Z1, Z2);
            dF_arr.push_back(F - F_old);
            // σ 用 log-Euler: Δln σ ≈ Δσ/σ
            dsigma_arr.push_back(std::log(sigma / sigma_old));
        }
    }

    // 样本相关性
    Size n = dF_arr.size();
    Real mean_dF = std::accumulate(dF_arr.begin(), dF_arr.end(), 0.0) / n;
    Real mean_ds = std::accumulate(dsigma_arr.begin(), dsigma_arr.end(), 0.0) / n;
    Real cov = 0.0, var_dF = 0.0, var_ds = 0.0;
    for (Size i = 0; i < n; ++i) {
        Real df = dF_arr[i] - mean_dF;
        Real ds = dsigma_arr[i] - mean_ds;
        cov += df * ds;
        var_dF += df * df;
        var_ds += ds * ds;
    }
    Real corr = cov / std::sqrt(var_dF * var_ds);
    // ρ = -0.8, 样本相关性应为负
    EXPECT_LT(corr, -0.3);  // 放宽容差, 因 Δln σ 与 Δσ 不完全相同
}

TEST(SABRProcessTest, ZeroRhoGivesUncorrelatedFAndSigma) {
    // ρ = 0: F 和 σ 独立, 样本相关性接近 0
    SABRParams p{0.2, 0.5, 0.0, 0.4};
    Real F0 = 100.0;
    SABRProcess sabr(p, F0);
    Size n_steps = 200;
    Real T = 1.0;
    Real dt = T / n_steps;

    std::vector<Real> dF_arr, dsigma_arr;
    Philox4x64 rng(42);

    Size n_paths = 200;
    for (Size i = 0; i < n_paths; ++i) {
        Real F = F0;
        Real sigma = p.alpha;
        for (Size k = 0; k < n_steps; ++k) {
            uint64_t r1 = rng();
            uint64_t r2 = rng();
            double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            auto [w1, w2] = box_muller(u1, u2);
            Real Z1 = w1;
            Real Z2 = w2;  // ρ=0, Z2 = W2
            Real F_old = F;
            Real sigma_old = sigma;
            F = sabr.evolve(F, sigma, dt, Z1, Z2);
            dF_arr.push_back(F - F_old);
            dsigma_arr.push_back(std::log(sigma / sigma_old));
        }
    }

    Size n = dF_arr.size();
    Real mean_dF = std::accumulate(dF_arr.begin(), dF_arr.end(), 0.0) / n;
    Real mean_ds = std::accumulate(dsigma_arr.begin(), dsigma_arr.end(), 0.0) / n;
    Real cov = 0.0, var_dF = 0.0, var_ds = 0.0;
    for (Size i = 0; i < n; ++i) {
        Real df = dF_arr[i] - mean_dF;
        Real ds = dsigma_arr[i] - mean_ds;
        cov += df * ds;
        var_dF += df * df;
        var_ds += ds * ds;
    }
    Real corr = cov / std::sqrt(var_dF * var_ds);
    EXPECT_NEAR(corr, 0.0, 0.10);  // 10% 容差
}

TEST(SABRProcessTest, DeterministicSeedReproducible) {
    SABRParams p{0.2, 0.5, -0.3, 0.4};
    SABRProcess sabr(p, 0.05);
    Size n_steps = 50;
    std::vector<Real> path1(n_steps + 1);
    std::vector<Real> path2(n_steps + 1);

    Philox4x64 rng1(42);
    Philox4x64 rng2(42);
    sabr.generate_path(1.0, n_steps, path1, rng1);
    sabr.generate_path(1.0, n_steps, path2, rng2);

    for (Size i = 0; i <= n_steps; ++i) {
        EXPECT_DOUBLE_EQ(path1[i], path2[i]);
    }
}

TEST(SABRProcessTest, ZeroNuKeepsSigmaConstant) {
    // ν = 0: σ 不随机, 保持 σ(0) = α
    // 但 generate_path 不输出 σ 路径, 用 evolve 直接验证
    SABRParams p{0.2, 0.5, 0.0, 0.0};  // ν = 0
    SABRProcess sabr(p, 100.0);
    Real F = 100.0;
    Real sigma = p.alpha;
    Real dt = 0.01;
    Philox4x64 rng(42);
    for (int i = 0; i < 100; ++i) {
        Real Z1, Z2;
        uint64_t r1 = rng();
        uint64_t r2 = rng();
        double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
        double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
        auto [w1, w2] = box_muller(u1, u2);
        Z1 = w1;
        Z2 = w2;
        F = sabr.evolve(F, sigma, dt, Z1, Z2);
        EXPECT_NEAR(sigma, p.alpha, 1e-12);  // σ 保持不变
    }
}
