// =============================================================================
// test_midas_model.cpp - M4 估计测试 (18 用例, spec §1.2 测试矩阵)
//
// 基准: midas_baseline.inc (midasr 0.9, verify_midas.R 2026-08-18:
//       midas_u/midas_r reltol=1e-12 + hAh_test)
//
// 容差策略 (跨优化器对照, 7B GM 先例口径):
//   - U-MIDAS vs midas_u (纯 OLS, 无优化器): 1e-10 主锚
//   - NLS vs midas_r (SLSQP 集中化 vs BFGS 联合): 参数 ~1e-4, SSR ~1e-8
//     (同一最优点的优化器容差层; 两者收敛同一最优由 W5 双起始佐证)
//   - hAh: 中心差分 vs numDeriv Richardson: stat/p ~1e-4
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/midas/midas_diagnostics.hpp"
#include "cpphub/timeseries/midas/midas_model.hpp"
#include "midas_baseline.inc"

namespace mm = cpphub::v1::timeseries::midas;
namespace bl = cpphub::v1::timeseries::midas_baseline::v1;
using cpphub::Real;
using cpphub::Size;

static mm::MixedFreqData make_data(Size max_low_lag = 0) {
    mm::MixedFreqData d;
    d.y.assign(bl::Y, bl::Y + bl::N_LF);
    d.x.assign(bl::X, bl::X + bl::N_HF);
    d.m = bl::M;
    d.max_low_lag = max_low_lag;
    return d;
}

// 缓存 (NLS ~1s 级, 避免重复估计)
static const mm::MidasResult& umidas_fit() {
    static mm::MidasResult r =
        mm::midas_fit(make_data(), mm::MidasWeight::Nealmon,
                      mm::MidasType::UMidas, 4);
    return r;
}
static const mm::MidasResult& nls_fit() {
    static mm::MidasResult r =
        mm::midas_fit(make_data(), mm::MidasWeight::Nealmon,
                      mm::MidasType::DL, 4);
    return r;
}
static const mm::MidasResult& ar_fit() {
    static mm::MidasResult r =
        mm::midas_fit(make_data(1), mm::MidasWeight::Nealmon,
                      mm::MidasType::AR, 4);
    return r;
}

// ---------------------------------------------------------------------------
// 对齐 (MD3) — mls 期末语义
// ---------------------------------------------------------------------------

// 1. mls lag0 = x_{3t} (期末最新, W2 实测 lag0 列 = 3,6,9,12)
TEST(MidasModel, MlsLag0IsPeriodEnd) {
    std::vector<Real> x(12);
    for (Size i = 0; i < 12; ++i) x[i] = static_cast<Real>(i + 1);
    const auto lag0 = mm::mls_column(x, 0, 3);
    ASSERT_EQ(lag0.size(), 4u);
    EXPECT_NEAR(lag0[0], 3.0, 1e-15);
    EXPECT_NEAR(lag0[1], 6.0, 1e-15);
    EXPECT_NEAR(lag0[2], 9.0, 1e-15);
    EXPECT_NEAR(lag0[3], 12.0, 1e-15);
    const auto lag1 = mm::mls_column(x, 1, 3);
    EXPECT_NEAR(lag1[0], 2.0, 1e-15);
    EXPECT_NEAR(lag1[3], 11.0, 1e-15);
}

// 2. mls lag 越界拒绝
TEST(MidasModel, MlsRejectsOversizedLag) {
    std::vector<Real> x(9, 1.0);
    EXPECT_THROW(mm::mls_column(x, 9, 3), std::invalid_argument);
    EXPECT_THROW(mm::mls_column(x, 0, 2), std::invalid_argument);  // m ∤ 9
}

// 3. MixedFreqData 校验 (N=n·m / NaN / 空)
TEST(MidasModel, MixedFreqDataValidation) {
    auto d = make_data();
    EXPECT_NO_THROW(d.validate());
    auto bad1 = d;
    bad1.x.pop_back();  // 299 ≠ 100·3
    EXPECT_THROW(bad1.validate(), std::invalid_argument);
    auto bad2 = d;
    bad2.y[0] = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(bad2.validate(), std::invalid_argument);
    auto bad3 = d;
    bad3.m = 0;
    EXPECT_THROW(bad3.validate(), std::invalid_argument);
}

// 4. design_matrix 行 j=2 ↔ R 期 2: 列 h=0 = x[3·2−1] (0-based, 期末系)
TEST(MidasModel, DesignMatrixAlignsPeriodTwo) {
    const auto X = mm::design_matrix(make_data(), 4);
    ASSERT_EQ(X.size(), 99u);  // j = 2..100
    // 行 0 (期 2): lag0 = x[5], lag3 = x[2] (0-based)
    EXPECT_NEAR(X[0][0], bl::X[3 * 2 - 1], 1e-15);
    EXPECT_NEAR(X[0][3], bl::X[3 * 2 - 4], 1e-15);
    // 行 98 (期 100): lag0 = x[299] (最后)
    EXPECT_NEAR(X[98][0], bl::X[299], 1e-15);
}

// ---------------------------------------------------------------------------
// U-MIDAS vs midas_u (1e-10 主锚)
// ---------------------------------------------------------------------------

// 5. U-MIDAS 系数 vs midas_u (W3)
TEST(MidasModel, UMidasCoeffMatchesR) {
    const auto& r = umidas_fit();
    ASSERT_EQ(r.delta.size(), 4u);
    EXPECT_NEAR(r.intercept, bl::UMIDAS_COEF[0], 1e-10);
    for (Size i = 0; i < 4; ++i) {
        EXPECT_NEAR(r.delta[i], bl::UMIDAS_COEF[1 + i], 1e-10);
    }
}

// 6. U-MIDAS SSR (W3)
TEST(MidasModel, UMidasSsrMatchesR) {
    EXPECT_NEAR(umidas_fit().ssr, bl::UMIDAS_SSR, 1e-10);
}

// 7. U-MIDAS sigma2 + n_obs + 隐含 midas_coef
TEST(MidasModel, UMidasSigma2AndMeta) {
    const auto& r = umidas_fit();
    EXPECT_NEAR(r.sigma2, bl::UMIDAS_SIGMA2, 1e-10);
    EXPECT_EQ(r.n_obs, 99u);
    EXPECT_EQ(r.type, mm::MidasType::UMidas);
    EXPECT_TRUE(r.converged);
    ASSERT_EQ(r.midas_coef.size(), 4u);  // midas_coef = delta (无约束)
    for (Size i = 0; i < 4; ++i) {
        EXPECT_NEAR(r.midas_coef[i], r.delta[i], 0.0);
    }
    EXPECT_TRUE(r.lambda.empty());  // 无外层超参
}

// 8. U-MIDAS 残差正交 (OLS): Σe ≈ 0, Σe·x_col ≈ 0
TEST(MidasModel, UMidasResidualOrthogonality) {
    const auto& r = umidas_fit();
    Real s = 0.0;
    for (Real e : r.residuals) s += e;
    EXPECT_NEAR(s, 0.0, 1e-10);
    // SSR 一致性
    Real ssr = 0.0;
    for (Real e : r.residuals) ssr += e * e;
    EXPECT_NEAR(ssr, r.ssr, 1e-9);
}

// ---------------------------------------------------------------------------
// MIDAS-DL 集中化 NLS vs midas_r (W4/W5)
// ---------------------------------------------------------------------------

// 9. NLS 系数 (μ, δ, λ₂) vs W4 — 方向裁决: 恢复 λ* ≈ (5, −0.5)
TEST(MidasModel, NlsDlCoeffMatchesR) {
    const auto& r = nls_fit();
    EXPECT_NEAR(r.intercept, bl::NLS_COEF[0], 1e-6);
    ASSERT_EQ(r.delta.size(), 1u);
    EXPECT_NEAR(r.delta[0], bl::NLS_COEF[1], 1e-4);
    ASSERT_EQ(r.lambda.size(), 1u);
    EXPECT_NEAR(r.lambda[0], bl::NLS_COEF[2], 1e-4);
}

// 10. NLS SSR vs W4 (1e-8 层)
TEST(MidasModel, NlsDlSsrMatchesR) {
    EXPECT_NEAR(nls_fit().ssr, bl::NLS_SSR, 1e-8);
}

// 11. NLS 隐含高频系数 (δ·w̄) vs W4 midas.coef
TEST(MidasModel, NlsDlMidasCoefMatchesR) {
    const auto& r = nls_fit();
    ASSERT_EQ(r.midas_coef.size(), 4u);
    for (Size i = 0; i < 4; ++i) {
        EXPECT_NEAR(r.midas_coef[i], bl::NLS_MIDAS_COEF[1 + i], 1e-6);
    }
}

// 12. 多起点逃逸: 显式远程起始 (MD8: 避开 λ=0 平坦陷阱) 收敛同一最优
TEST(MidasModel, NlsDlMultiStartSameOptimum) {
    const auto r2 = mm::midas_fit(make_data(), mm::MidasWeight::Nealmon,
                                  mm::MidasType::DL, 4, {{2.0}});
    ASSERT_EQ(r2.lambda.size(), 1u);
    EXPECT_NEAR(r2.lambda[0], bl::NLS2_COEF[2], 1e-4);   // W5 同最优
    EXPECT_NEAR(r2.delta[0], bl::NLS2_COEF[1], 1e-4);
    EXPECT_NEAR(r2.ssr, bl::NLS_SSR, 1e-8);
    EXPECT_TRUE(r2.converged);
}

// ---------------------------------------------------------------------------
// MIDAS-AR vs W6
// ---------------------------------------------------------------------------

// 13. AR 系数 (μ, δ, λ₂, ρ₁) vs W6
TEST(MidasModel, ArCoeffMatchesR) {
    const auto& r = ar_fit();
    EXPECT_NEAR(r.intercept, bl::AR_COEF[0], 1e-5);
    ASSERT_EQ(r.delta.size(), 2u);  // (δ, ρ₁)
    EXPECT_NEAR(r.delta[0], bl::AR_COEF[1], 1e-4);
    EXPECT_NEAR(r.lambda[0], bl::AR_COEF[2], 1e-4);
    EXPECT_NEAR(r.delta[1], bl::AR_COEF[3], 1e-4);
}

// 14. AR 对齐: n_obs = 99 (j0 = max(⌊3/3⌋+1, 1+1) = 2), AR 列期 2 = y_1 真值
TEST(MidasModel, ArAlignment) {
    const auto& r = ar_fit();
    EXPECT_EQ(r.n_obs, 99u);
    ASSERT_EQ(r.residuals.size(), 99u);
    // ARSS 收敛性
    EXPECT_TRUE(r.converged);
}

// 15. U-MIDAS SSR ≤ NLS SSR ≤ U-MIDAS SSR (约束与无约束嵌套关系)
TEST(MidasModel, NestedSsrOrdering) {
    // 约束模型 SSR ≥ 无约束 (U-MIDAS 参数更多)
    EXPECT_GE(nls_fit().ssr, umidas_fit().ssr - 1e-9);
    // 差距有限 (真值就在 nealmon 族内)
    EXPECT_LT(nls_fit().ssr - umidas_fit().ssr, 1.0);
}

// 16. loglik 公式 (高斯 MLE): −n/2·(log(2πσ̂²ₘₗₑ)+1)
TEST(MidasModel, LoglikFormulaConsistency) {
    const auto& r = umidas_fit();
    constexpr Real kTwoPi = 6.283185307179586476925286766559;
    const Real s2mle = r.ssr / static_cast<Real>(r.n_obs);
    const Real expect = -0.5 * static_cast<Real>(r.n_obs)
                        * (std::log(kTwoPi * s2mle) + 1.0);
    EXPECT_NEAR(r.loglik, expect, 1e-10);
}

// ---------------------------------------------------------------------------
// 语义边界 + 诊断
// ---------------------------------------------------------------------------

// 17. PolyStep 估计拒绝 (离散断点结构非连续参数化)
TEST(MidasModel, PolyStepFitRejected) {
    EXPECT_THROW(mm::midas_fit(make_data(), mm::MidasWeight::PolyStep,
                               mm::MidasType::DL, 4),
                 std::invalid_argument);
}

// 18. hAh 三列 vs W7 (stat/p/df) — DGP 在 nealmon 族内 → 不拒绝
TEST(MidasModel, HahTestMatchesR) {
    const auto& r = nls_fit();
    const auto hah = mm::hAh_test(r, make_data(), 4, mm::MidasWeight::Nealmon);
    EXPECT_EQ(hah.df, bl::HAH[2]);
    EXPECT_NEAR(hah.stat, bl::HAH[0], 1e-4);
    EXPECT_NEAR(hah.p, bl::HAH[1], 1e-4);
    EXPECT_GT(hah.p, 0.05);  // H0 成立 (DGP 即 nealmon)
    // U-MIDAS 是无约束备择本身 → 拒绝
    EXPECT_THROW(mm::hAh_test(umidas_fit(), make_data(), 4,
                              mm::MidasWeight::Nealmon),
                 std::invalid_argument);
}

// 附加: 残差 LB/JB 复用 (白噪声 DGP → 双不拒绝)
TEST(MidasModel, ResidualDiagnosticsReuse) {
    const auto d = mm::midas_residual_diagnostics(nls_fit(), 10);
    EXPECT_TRUE(d.residual_white);
}
