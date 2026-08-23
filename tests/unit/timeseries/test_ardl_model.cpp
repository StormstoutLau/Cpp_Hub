// =============================================================================
// test_ardl_model.cpp - ARDL 设计矩阵/OLS/A8 计账镜像 (Phase 8 M1 Step 2, 14 用例)
//
// 基准: ardl_baseline.inc (statsmodels 0.14.6 ARDL.fit 主基准 1e-10 +
// ctt 手工 numpy OLS 基准; verify_ardl.py [S1]/[S2] 逐值自证)
// 幻觉点: A6 (causal=false 同期项) / A8 (n_acct = T − hb_early 计账:
//         llf/sigma2/IC ≠ 真实估计样本口径, 显式 hb 或 max q > p 时暴露) /
//         PSS8 (hold_back 默认 = max(p, max q_j)) / trend 列全样本 1-based
//         (≠ VAR 估计样本内重编号, 系数对齐即隐式验证)
// =============================================================================

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/ardl/ardl_model.hpp"
#include "ardl_baseline.inc"

namespace ardl = cpphub::v1::timeseries::ardl;
namespace ab = cpphub::v1::ardl_baseline;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

template <Size N>
std::vector<Real> vec(const double (&a)[N]) {
    return std::vector<Real>(std::begin(a), std::end(a));
}

std::vector<std::vector<Real>> exog_k1() {
    return {vec(ab::DATA_K1_X1)};
}

std::vector<std::vector<Real>> exog_k2() {
    return {vec(ab::DATA_K2_X1), vec(ab::DATA_K2_X2)};
}

/// 通用基准对照: params + bse(√diag vcov) + SCALARS 九元组 + 回显
/// hb_early: A8 计账早期 hold_back = 显式 hb 或 max(y_lags)=p (仅 y 滞后)
void check_against(const ardl::ArdlResult& r, const double params[],
                   Size np, const double bse[], const double scalars[],
                   Real tol, Size hb_early) {
    ASSERT_EQ(r.params.size(), np);
    for (Size i = 0; i < np; ++i) {
        EXPECT_NEAR(r.params[i], params[i], tol) << "param " << i;
        EXPECT_NEAR(std::sqrt(r.vcov[i][i]), bse[i], tol) << "bse " << i;
    }
    // SCALARS = [llf, aic, bic, hqic, sigma2_df, sigma2_acct, hold_back,
    //            nobs(=n_acct), K]
    EXPECT_NEAR(r.loglik, scalars[0], tol);
    EXPECT_NEAR(r.aic, scalars[1], tol);
    EXPECT_NEAR(r.bic, scalars[2], tol);
    EXPECT_NEAR(r.hqic, scalars[3], tol);
    EXPECT_NEAR(r.sigma2, scalars[5], tol);      // A8 口径 (= res.sigma2)
    EXPECT_EQ(r.hold_back_used, static_cast<Size>(scalars[6]));
    // A8: n_acct = T − hb_early (T = 250 夹具)
    EXPECT_EQ(static_cast<Size>(scalars[7]), 250u - hb_early);
    EXPECT_EQ(r.params.size(), static_cast<Size>(scalars[8]));
    // sigma2_df 自洽: ssr 由 residuals 重建
    Real ssr = 0.0;
    for (Real e : r.residuals) ssr += e * e;
    const Size n_act = r.n_obs_used;
    EXPECT_NEAR(ssr / static_cast<Real>(n_act - r.params.size()),
                scalars[4], tol);
    // llf 镜像自证 (A8 分母)
    const Real n_a = scalars[7];
    const Real llf = -n_a / 2.0 * (std::log(2.0 * std::acos(-1.0)) +
                                   std::log(ssr / n_a) + 1.0);
    EXPECT_NEAR(llf, r.loglik, tol);
}

}  // namespace

// 1. trend="c" p=2 q=[1]: 系数/bse/IC vs statsmodels 1e-10
TEST(ArdlModel, CoefficientsC_P2Q1_VsStatsmodels) {
    ardl::ArdlSpec s;
    s.p = 2;
    s.q = {1};
    auto r = ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s);
    check_against(r, ab::PARAMS_C_P2Q1, 5, ab::BSE_C_P2Q1, ab::SCALARS_C_P2Q1,
                  1e-10, /*hb_early=*/2);
}

// 2. trend="n": 无确定性列
TEST(ArdlModel, CoefficientsN_P2Q1_VsStatsmodels) {
    ardl::ArdlSpec s;
    s.p = 2;
    s.q = {1};
    s.trend = "n";
    auto r = ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s);
    check_against(r, ab::PARAMS_N_P2Q1, 4, ab::BSE_N_P2Q1, ab::SCALARS_N_P2Q1,
                  1e-10, /*hb_early=*/2);
}

// 3. trend="ct" p=1 q=[2]: A8 暴露用例 — n_act=248 但 n_acct=249
TEST(ArdlModel, CoefficientsCt_P1Q2_VsStatsmodels_A8Exposed) {
    ardl::ArdlSpec s;
    s.p = 1;
    s.q = {2};
    s.trend = "ct";
    auto r = ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s);
    check_against(r, ab::PARAMS_CT_P1Q2, 6, ab::BSE_CT_P1Q2, ab::SCALARS_CT_P1Q2,
                  1e-10, /*hb_early=*/1);
    // A8 显式: 真实估计样本 248, statsmodels 计账 249 = T − p
    EXPECT_EQ(r.n_obs_used, 248u);
    EXPECT_EQ(r.hold_back_used, 2u);
    // llf 若误用 n_act (248) 将偏离基准 ~1.4 → A8 镜像已被 check_against 锁死
}

// 4. trend="ctt": statsmodels ARDL 无此档 — 手工 numpy 基准
//    (语义正确性由构造式 [1, t, t²] 直接保证, verify_ardl.py)
TEST(ArdlModel, CoefficientsCtt_P1Q1_VsManualBaseline) {
    ardl::ArdlSpec s;
    s.p = 1;
    s.q = {1};
    s.trend = "ctt";
    auto r = ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s);
    check_against(r, ab::PARAMS_CTT_P1Q1, 6, ab::BSE_CTT_P1Q1,
                  ab::SCALARS_CTT_P1Q1, 1e-10, /*hb_early=*/1);
}

// 5. k=2 双变量逐阶 q=[1,2]
TEST(ArdlModel, CoefficientsC_K2_VsStatsmodels) {
    ardl::ArdlSpec s;
    s.p = 2;
    s.q = {1, 2};
    auto r = ardl::ardl_fit(vec(ab::DATA_K2_Y), exog_k2(), s);
    check_against(r, ab::PARAMS_C_K2, 8, ab::BSE_C_K2, ab::SCALARS_C_K2, 1e-10,
                  /*hb_early=*/2);
    EXPECT_EQ(r.n_obs_used, 248u);
}

// 6. 边界阶: p=0 (纯 DL) 与 q=0 (仅同期)
TEST(ArdlModel, BoundaryOrders_P0Q1_And_P2Q0) {
    ardl::ArdlSpec a;
    a.p = 0;
    a.q = {1};
    auto ra = ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), a);
    check_against(ra, ab::PARAMS_C_P0Q1, 3, ab::BSE_C_P0Q1, ab::SCALARS_C_P0Q1,
                  1e-10, /*hb_early=*/0);  // p=0 → hb_early=0, n_acct=T=250
    EXPECT_EQ(ra.n_obs_used, 249u);   // hb = max(0, 1) = 1
    EXPECT_EQ(ra.hold_back_used, 1u);

    ardl::ArdlSpec b;
    b.p = 2;
    b.q = {0};
    auto rb = ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), b);
    check_against(rb, ab::PARAMS_C_P2Q0, 4, ab::BSE_C_P2Q0, ab::SCALARS_C_P2Q0,
                  1e-10, /*hb_early=*/2);
}

// 7. 显式 hold_back=5: n_acct = T − hb (hb_early = 显式值)
TEST(ArdlModel, HoldBackExplicit_HB5) {
    ardl::ArdlSpec s;
    s.p = 2;
    s.q = {1};
    s.hold_back = 5;
    auto r = ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s);
    check_against(r, ab::PARAMS_C_HB5, 5, ab::BSE_C_HB5, ab::SCALARS_C_HB5,
                  1e-10, /*hb_early=*/5);
    EXPECT_EQ(r.n_obs_used, 245u);
    EXPECT_EQ(r.hold_back_used, 5u);
    EXPECT_EQ(r.spec.hold_back, 5u);
}

// 8. k=0 纯 AR 退化 (bounds 侧将拒绝, 拟合允许)
TEST(ArdlModel, PureAR_K0) {
    ardl::ArdlSpec s;
    s.p = 2;
    s.q = {};
    auto r = ardl::ardl_fit(vec(ab::DATA_K1_Y), {}, s);
    check_against(r, ab::PARAMS_AR2, 3, ab::BSE_AR2, ab::SCALARS_AR2, 1e-10,
                  /*hb_early=*/2);
}

// 9. fitted/residuals: 长度 + y 重建 + OLS 正交性
TEST(ArdlModel, FittedAndResiduals) {
    ardl::ArdlSpec s;
    s.p = 2;
    s.q = {1};
    auto r = ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s);
    const auto& y = ab::DATA_K1_Y;
    ASSERT_EQ(r.fitted.size(), r.n_obs_used);
    ASSERT_EQ(r.residuals.size(), r.n_obs_used);
    for (Size i = 0; i < r.n_obs_used; ++i) {
        EXPECT_NEAR(r.fitted[i] + r.residuals[i], y[r.hold_back_used + i],
                    1e-12);
    }
    // 残差与设计矩阵各列正交 (OLS 一阶条件): 用 X 列重建核验 fitted 亦可,
    // 此处核验 E[resid] ≈ 0 与 resid·y_lag 正交 (数值抽样)
    Real dot = 0.0;
    for (Size i = 0; i < r.n_obs_used; ++i) {
        dot += r.residuals[i] * y[r.hold_back_used + i - 1];  // vs y.L1
    }
    EXPECT_NEAR(dot, 0.0, 1e-8);
}

// 10. 默认 hold_back = max(p, max_j q_j) 回显 (PSS8)
TEST(ArdlModel, HoldBackDefaultEcho) {
    ardl::ArdlSpec s;
    s.p = 2;
    s.q = {1, 3};
    auto r = ardl::ardl_fit(vec(ab::DATA_K2_Y), exog_k2(), s);
    EXPECT_EQ(r.hold_back_used, 3u);
    EXPECT_EQ(r.spec.hold_back, 3u);
    EXPECT_EQ(r.n_obs_used, 247u);
}

// 11. A8 计账结构断言: n_acct 与 n_obs_used 的偏离量 = hb_early − hb_eff
TEST(ArdlModel, A8AccountingDivergence) {
    // p=1, q=2 (默认 hb): hb_eff=2, hb_early=1 → n_acct = n_act + 1
    ardl::ArdlSpec s;
    s.p = 1;
    s.q = {2};
    auto r = ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s);
    const Real n_acct = 250.0 - 1.0;
    Real ssr = 0.0;
    for (Real e : r.residuals) ssr += e * e;
    const Real llf_a8 = -n_acct / 2.0 * (std::log(2.0 * std::acos(-1.0)) +
                                         std::log(ssr / n_acct) + 1.0);
    EXPECT_NEAR(r.loglik, llf_a8, 1e-10);
    // 对照: 若用 n_act (248) 计账, llf 差异 ~1.0 (阻断性; 阈值 0.5 防脆弱)
    const Real n_act = static_cast<Real>(r.n_obs_used);
    const Real llf_naive =
        -n_act / 2.0 * (std::log(2.0 * std::acos(-1.0)) +
                        std::log(ssr / n_act) + 1.0);
    EXPECT_GT(std::abs(llf_a8 - llf_naive), 0.5);
}

// 12. 校验: causal / trend 非法 / 不等长 / NaN / q 不匹配 / k>10 / 样本过小
TEST(ArdlModel, ValidationThrows) {
    const auto y = vec(ab::DATA_K1_Y);
    const auto x = exog_k1();
    ardl::ArdlSpec s;
    s.p = 1;
    s.q = {1};

    ardl::ArdlSpec bad = s;
    bad.causal = true;
    EXPECT_THROW(ardl::ardl_fit(y, x, bad), std::invalid_argument);

    bad = s;
    bad.trend = "cttt";
    EXPECT_THROW(ardl::ardl_fit(y, x, bad), std::invalid_argument);

    bad = s;
    bad.q = {1, 1};  // q.size() != x.size()
    EXPECT_THROW(ardl::ardl_fit(y, x, bad), std::invalid_argument);

    auto x2 = exog_k1();
    x2[0].pop_back();  // 不等长
    EXPECT_THROW(ardl::ardl_fit(y, x2, s), std::invalid_argument);

    auto xn = exog_k1();
    xn[0][10] = std::nan("");
    EXPECT_THROW(ardl::ardl_fit(y, xn, s), std::invalid_argument);

    bad = s;
    bad.q.assign(11, 1);  // k > 10 (响应面覆盖域)
    std::vector<std::vector<Real>> x11(11, x[0]);
    EXPECT_THROW(ardl::ardl_fit(y, x11, bad), std::invalid_argument);

    bad = s;
    bad.p = 249;  // 样本过小 (hb=249, npar>剩余)
    EXPECT_THROW(ardl::ardl_fit(y, x, bad), std::invalid_argument);
}

// 13. vcov 对称性与 bse/vcov 一致性 (A7 nonrobust)
TEST(ArdlModel, VcovSymmetry) {
    ardl::ArdlSpec s;
    s.p = 2;
    s.q = {1};
    auto r = ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s);
    for (Size i = 0; i < r.params.size(); ++i) {
        for (Size j = i + 1; j < r.params.size(); ++j) {
            EXPECT_NEAR(r.vcov[i][j], r.vcov[j][i], 1e-14);
        }
    }
}

// 14. IC 三式惩罚项自洽 (σ² 计入 → 恒 K+1)
TEST(ArdlModel, ICFormulaSelfConsistency) {
    ardl::ArdlSpec s;
    s.p = 2;
    s.q = {1};
    auto r = ardl::ardl_fit(vec(ab::DATA_K1_Y), exog_k1(), s);
    const Real m2llf = -2.0 * r.loglik;
    const Real Kp1 = static_cast<Real>(r.params.size()) + 1.0;
    EXPECT_NEAR(r.aic - m2llf, 2.0 * Kp1, 1e-9);
    EXPECT_NEAR(r.bic - m2llf, std::log(248.0) * Kp1, 1e-9);
    EXPECT_NEAR(r.hqic - m2llf,
                2.0 * std::log(std::log(248.0)) * Kp1, 1e-9);
}
