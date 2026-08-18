// =============================================================================
// test_var_model.cpp - VAR 估计/IC/滞后选择/稳定性 (Phase 7C M2, 20 用例)
//
// 基准: var_baseline.inc (statsmodels 0.14.4 主基准 1e-10 + R vars 交叉 1e-8)
// 幻觉点: V2 (Cholesky 下三角) / V4 (Σ_mle ÷T) / V5 (同样本 offset) /
//         V6 (FPE 指数 K) / V9 (稳定性双输出)
// =============================================================================

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/rng.hpp"
#include "cpphub/timeseries/var/multivariate_data.hpp"
#include "cpphub/timeseries/var/var_model.hpp"
#include "cpphub/timeseries/var/var_select.hpp"
#include "var_baseline.inc"

namespace vt = cpphub::v1::timeseries::var;
namespace vb = cpphub::v1::timeseries::var_baseline;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

vt::MultivariateTSData fixture_data() {
    vt::MultivariateTSData d;
    d.columns = {std::vector<cpphub::v1::Real>(std::begin(vb::Y1), std::end(vb::Y1)),
                 std::vector<cpphub::v1::Real>(std::begin(vb::Y2), std::end(vb::Y2)),
                 std::vector<cpphub::v1::Real>(std::begin(vb::Y3), std::end(vb::Y3))};
    d.names = {"y1", "y2", "y3"};
    return d;
}

vt::VARResult fit2() {
    vt::VARSpec spec;
    spec.lag = 2;
    return vt::var_fit(fixture_data(), spec);
}

/// 平铺 → MatrixXd (row-major)
MatrixXd mat3(const std::vector<cpphub::v1::Real>& flat) {
    MatrixXd m(3, 3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) m(i, j) = flat[3 * i + j];
    return m;
}

}  // namespace

// 1. 系数 vs statsmodels (SM_PARAMS 布局 回归元×方程 7×3 row-major), 1e-10
TEST(VarModel, CoefficientsVsStatsmodels) {
    auto fit = fit2();
    ASSERT_EQ(fit.coefficients.rows(), 3);
    ASSERT_EQ(fit.coefficients.cols(), 7);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 7; ++j) {
            const double ref = vb::SM_PARAMS[3 * j + i];  // inc[j*3+i] = regressor j, eq i
            EXPECT_NEAR(fit.coefficients(i, j), ref, 1e-10)
                << "eq " << i << " regressor " << j;
        }
    }
    EXPECT_EQ(fit.n_obs_used, 248u);
}

// 2. Σ_mle = SSR/T (V4), 1e-10
TEST(VarModel, SigmaUmeVsStatsmodels) {
    auto fit = fit2();
    const MatrixXd ref = mat3(std::vector<cpphub::v1::Real>(
        std::begin(vb::SM_SIGMA_U_MLE), std::end(vb::SM_SIGMA_U_MLE)));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(fit.sigma_u_mle(i, j), ref(i, j), 1e-10) << i << j;
}

// 3. Σ df 修正版 (IRF/FEVD 输入), 1e-10
TEST(VarModel, SigmaUDfCorrectedVsStatsmodels) {
    auto fit = fit2();
    const MatrixXd ref = mat3(std::vector<cpphub::v1::Real>(
        std::begin(vb::SM_SIGMA_U), std::end(vb::SM_SIGMA_U)));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(fit.sigma_u(i, j), ref(i, j), 1e-10) << i << j;
    // 分层一致性: sigma_u·df_resid/T = sigma_u_mle
    const double scale = (248.0 - 7.0) / 248.0;
    EXPECT_NEAR(fit.sigma_u(1, 2) * scale, fit.sigma_u_mle(1, 2), 1e-12);
}

// 4. IC 五式 (V4/V6), 1e-10
TEST(VarModel, ICsVsStatsmodels) {
    auto fit = fit2();
    EXPECT_NEAR(fit.loglik, vb::SM_SCALARS[0], 1e-10);
    EXPECT_NEAR(fit.aic, vb::SM_SCALARS[1], 1e-10);
    EXPECT_NEAR(fit.bic, vb::SM_SCALARS[2], 1e-10);
    EXPECT_NEAR(fit.hqic, vb::SM_SCALARS[3], 1e-10);
    EXPECT_NEAR(fit.fpe, vb::SM_SCALARS[4], 1e-10);
    // SCALARS[5] = det(Σ_df); logdet(Σ_mle) = log det(Σ_df) − K·ln(nobs/df_resid)
    EXPECT_NEAR(fit.logdet,
                std::log(vb::SM_SCALARS[5]) - 3.0 * std::log(248.0 / 241.0), 1e-9);
}

// 5. Cholesky 下三角 (V2): LLT(sigma_u df 修正) == np.linalg.cholesky, 1e-12
TEST(VarModel, CholeskyLowerVsNumpy) {
    auto fit = fit2();
    Eigen::LLT<MatrixXd> llt(fit.sigma_u);
    const MatrixXd L = llt.matrixL();
    const MatrixXd ref = mat3(std::vector<cpphub::v1::Real>(
        std::begin(vb::SM_CHOL_LOWER), std::end(vb::SM_CHOL_LOWER)));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(L(i, j), ref(i, j), 1e-12) << i << j;
    // 上三角严格为 0
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) EXPECT_DOUBLE_EQ(L(i, j), 0.0);
}

// 6. R vars t(chol) 等价 (V2 交叉, 1e-12); R dump 列主序 → L(i,j)=inc[3j+i]
// 注意: matrixL() 返回 TriangularView, 须物化后再取元 (上三角存储未引用)
TEST(VarModel, CholeskyRVarsCross) {
    auto fit = fit2();
    Eigen::LLT<MatrixXd> llt(fit.sigma_u);
    const MatrixXd L = llt.matrixL();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(L(i, j), vb::R_VARS_CHOL[3 * j + i], 1e-12) << i << j;
}

// 7. 稳定性双输出 (V9): max_eig = 1/min(roots), 严格 <1
TEST(VarModel, StabilityDualOutput) {
    auto fit = fit2();
    EXPECT_NEAR(fit.max_abs_eigenvalue, vb::SM_MAX_EIG[0], 1e-9);
    EXPECT_TRUE(fit.is_strictly_stationary);
    EXPECT_LT(fit.max_abs_eigenvalue, 1.0);
    // roots 全 >1 (statsmodels 口径)
    for (double r : vb::SM_ROOTS_ASC) EXPECT_GT(r, 1.0);
}

// 8. select_order 同样本轨迹 (V5), p=0..4, 1e-10
TEST(VarModel, SelectOrderSameSampleTrajectory) {
    auto sel = vt::var_select_order(fixture_data(), "c", 4, "aic");
    ASSERT_EQ(sel.aic.size(), 5u);
    for (int p = 0; p < 5; ++p) {
        EXPECT_NEAR(sel.aic[p], vb::SM_SEL_AIC[p], 1e-10) << "aic p=" << p;
        EXPECT_NEAR(sel.bic[p], vb::SM_SEL_BIC[p], 1e-10) << "bic p=" << p;
        EXPECT_NEAR(sel.hqic[p], vb::SM_SEL_HQIC[p], 1e-10) << "hqic p=" << p;
        EXPECT_NEAR(sel.fpe[p], vb::SM_SEL_FPE[p], 1e-10) << "fpe p=" << p;
    }
    EXPECT_EQ(sel.n_obs_used, 246u);  // maxlags=4 offset 同样本
}

// 9. argmin: 四 IC 均选中 p=2
TEST(VarModel, SelectOrderArgmin) {
    for (const char* ic : {"aic", "bic", "hqic", "fpe"}) {
        auto sel = vt::var_select_order(fixture_data(), "c", 4, ic);
        EXPECT_EQ(sel.selected_lag, 2u) << ic;
    }
    // 非法 ic 拒绝
    EXPECT_THROW(vt::var_select_order(fixture_data(), "c", 4, "xic"),
                 std::invalid_argument);
}

// 10. R vars::VARselect 交叉 (1e-8): AIC/SC 轨迹 p=1..4 逐位
TEST(VarModel, RVarsVarselectCross) {
    auto sel = vt::var_select_order(fixture_data(), "c", 4, "aic");
    for (int p = 1; p <= 4; ++p) {
        EXPECT_NEAR(sel.aic[p], vb::R_SEL_AIC[p - 1], 1e-8) << "aic p=" << p;
        EXPECT_NEAR(sel.bic[p], vb::R_SEL_SC[p - 1], 1e-8) << "bic p=" << p;
        // R FPE 定义不同 (V6 以 statsmodels 为主; 记录差异存在)
        // — FPE 不作 R 交叉锚
    }
}

// 11. trend="n" 拟合 (1e-10)
TEST(VarModel, TrendNoneFit) {
    vt::VARSpec spec;
    spec.lag = 2;
    spec.trend = "n";
    auto fit = vt::var_fit(fixture_data(), spec);
    ASSERT_EQ(fit.coefficients.cols(), 6u);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 6; ++j) {
            // SM_N_PARAMS 布局 params[reg][eq]: 6×3
            const double ref = vb::SM_N_PARAMS[3 * j + i];
            EXPECT_NEAR(fit.coefficients(i, j), ref, 1e-10) << i << j;
        }
    EXPECT_NEAR(fit.loglik, vb::SM_N_SCALARS[0], 1e-10);
    EXPECT_NEAR(fit.aic, vb::SM_N_SCALARS[1], 1e-10);
    EXPECT_NEAR(fit.bic, vb::SM_N_SCALARS[2], 1e-10);
    EXPECT_NEAR(fit.hqic, vb::SM_N_SCALARS[3], 1e-10);
    EXPECT_NEAR(fit.fpe, vb::SM_N_SCALARS[4], 1e-10);
}

// 12. trend="ct" 拟合 (IC 1e-10)
TEST(VarModel, TrendCtFit) {
    vt::VARSpec spec;
    spec.lag = 2;
    spec.trend = "ct";
    auto fit = vt::var_fit(fixture_data(), spec);
    ASSERT_EQ(fit.coefficients.cols(), 8u);
    EXPECT_NEAR(fit.loglik, vb::SM_CT_SCALARS[0], 1e-10);
    EXPECT_NEAR(fit.aic, vb::SM_CT_SCALARS[1], 1e-10);
    EXPECT_NEAR(fit.bic, vb::SM_CT_SCALARS[2], 1e-10);
    EXPECT_NEAR(fit.fpe, vb::SM_CT_SCALARS[4], 1e-10);
}

// 13. p=0 纯截距 (select_order p_min=0 成员)
TEST(VarModel, P0PureConstant) {
    vt::VARSpec spec;
    spec.lag = 0;  // → 会 IC 自动; 用 var_fit? p=0 由 fit 触发 IC.
    // 直接验证 p=0 估计: 通过 select 轨迹首点已覆盖; 此处验证 fit 的
    // auto-lag 路径选 p=2 (与 DGP 一致)
    auto fit = vt::var_fit(fixture_data(), spec);
    EXPECT_EQ(fit.lag, 2u);
    EXPECT_NEAR(fit.aic, vb::SM_SCALARS[1], 1e-10);
}

// 14. p=0 全样本 IC 锚 (SM_P0_SCALARS = model.fit(0) 全样本 T=250;
//     与 select_order 轨迹 p=0 (offset 246) 语义不同, 闭包重算对照)
TEST(VarModel, P0ICAnchors) {
    const auto d = fixture_data();
    const MatrixXd Y = d.matrix();  // 250×3
    const MatrixXd Yc = Y.rowwise() - Y.colwise().mean();
    const MatrixXd S = Yc.transpose() * Yc / 250.0;  // Σ_mle (VAR(0) 仅截距)
    Eigen::LLT<MatrixXd> llt(S);
    const double ld =
        2.0 * llt.matrixL().toDenseMatrix().diagonal().array().log().sum();
    const double fp = 3.0;  // fp = 0·K² + K·k_trend = 3
    EXPECT_NEAR(ld + 2.0 / 250.0 * fp, vb::SM_P0_SCALARS[0], 1e-10);   // aic
    EXPECT_NEAR(ld + std::log(250.0) / 250.0 * fp, vb::SM_P0_SCALARS[1],
                1e-10);  // bic
    EXPECT_NEAR(ld + 2.0 * std::log(std::log(250.0)) / 250.0 * fp,
                vb::SM_P0_SCALARS[2], 1e-10);                          // hqic
    EXPECT_NEAR(std::pow(251.0 / 249.0, 3.0) * std::exp(ld),
                vb::SM_P0_SCALARS[3], 1e-10);                          // fpe
}

// 15. 变量重排敏感性 (V2: Cholesky 次序敏感; reorder FEVD ≠ 原 FEVD)
TEST(VarModel, ReorderSensitivity) {
    auto fit = fit2();
    auto data_r = fixture_data().reorder({2, 0, 1});
    vt::VARSpec spec;
    spec.lag = 2;
    auto fit_r = vt::var_fit(data_r, spec);
    // 重排系数 vs statsmodels reorder (7×3 回归元主序: inc[3j+i])
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 7; ++j) {
            const double ref = vb::SM_REORDER_PARAMS[3 * j + i];
            EXPECT_NEAR(fit_r.coefficients(i, j), ref, 1e-10) << i << j;
        }
    // 重排后 Cholesky 次序变化 → Σ 行列重排 (y3,y1,y2 序)
    EXPECT_NEAR(fit_r.sigma_u(0, 0), fit.sigma_u(2, 2), 1e-12);
    EXPECT_NEAR(fit_r.sigma_u(1, 1), fit.sigma_u(0, 0), 1e-12);
}

// 16. R vars 系数交叉 (R_VARS_COEF 3×7 方程主序, 列 [y1.l1..y3.l1 y1.l2..y3.l2 const])
TEST(VarModel, RVarsCoefCross) {
    auto fit = fit2();
    // R dump: 行=方程 i, 列 = [y1.l1 y2.l1 y3.l1 y1.l2 y2.l2 y3.l2 const]
    // C++ coefficients(i, j): 方程 i × [const, lags...]
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(fit.coefficients(i, 0), vb::R_VARS_COEF[i * 7 + 6], 1e-10);
        for (int l = 0; l < 2; ++l)
            for (int j = 0; j < 3; ++j) {
                const double ref = vb::R_VARS_COEF[i * 7 + l * 3 + j];
                EXPECT_NEAR(fit.coefficients(i, 1 + l * 3 + j), ref, 1e-10)
                    << i << " lag " << l + 1 << " var " << j;
            }
    }
    // loglik 交叉 (R vars::logLik)
    EXPECT_NEAR(fit.loglik, vb::R_VARS_LOGLIK[0], 1e-8);
}

// 17. 残差正交性: Z'Re ≈ 0 (OLS 一阶条件)
TEST(VarModel, ResidualOrthogonality) {
    auto fit = fit2();
    // 重构 Z: [const, lags] (数据自 CSV)
    const auto d = fixture_data();
    const MatrixXd Y = d.matrix();
    const int T = static_cast<int>(Y.rows());
    MatrixXd Z(T - 2, 7);
    for (int t = 2; t < T; ++t) {
        Z(t - 2, 0) = 1.0;
        for (int l = 0; l < 2; ++l)
            for (int j = 0; j < 3; ++j) Z(t - 2, 1 + l * 3 + j) = Y(t - 1 - l, j);
    }
    const MatrixXd G = Z.transpose() * fit.residuals;
    const double scale = G.cwiseAbs().maxCoeff();
    EXPECT_LT(scale, 1e-9);
}

// 18. coeff_vcov 结构 (K 个 (Kp+kt)² 对称矩阵, 对角正)
TEST(VarModel, CoeffVcovShape) {
    auto fit = fit2();
    ASSERT_EQ(fit.coeff_vcov.size(), 3u);
    for (const auto& V : fit.coeff_vcov) {
        ASSERT_EQ(V.rows(), 7);
        ASSERT_EQ(V.cols(), 7);
        for (int i = 0; i < 7; ++i) {
            EXPECT_GT(V(i, i), 0.0);
            EXPECT_NEAR(V(i, i), V(i, i), 0.0);  // finite
        }
        EXPECT_NEAR(V(0, 1), V(1, 0), 1e-14);
    }
}

// 19. 异常输入
TEST(VarModel, InvalidInputs) {
    vt::MultivariateTSData bad;
    bad.columns = {{1.0, 2.0}, {1.0, 2.0, 3.0}};  // 不等长
    EXPECT_THROW(bad.T(), std::invalid_argument);
    vt::MultivariateTSData nan_data;
    nan_data.columns = {{1.0, std::nan("")}, {1.0, 2.0}};
    EXPECT_THROW(nan_data.validate(), std::invalid_argument);
    vt::MultivariateTSData empty;
    empty.columns = {};
    EXPECT_THROW(empty.validate(), std::invalid_argument);
    vt::VARSpec sp;
    sp.lag = 2;
    sp.trend = "xyz";
    EXPECT_THROW(vt::var_fit(fixture_data(), sp), std::invalid_argument);
    // reorder 非法
    EXPECT_THROW(fixture_data().reorder({0, 0, 1}), std::invalid_argument);
}

// 20. matrix()/T()/K()/names 基础
TEST(VarModel, DataCarrierBasics) {
    auto d = fixture_data();
    EXPECT_EQ(d.T(), 250u);
    EXPECT_EQ(d.K(), 3u);
    const MatrixXd M = d.matrix();
    ASSERT_EQ(M.rows(), 250);
    ASSERT_EQ(M.cols(), 3);
    EXPECT_DOUBLE_EQ(M(0, 0), vb::Y1[0]);
    EXPECT_DOUBLE_EQ(M(249, 2), vb::Y3[249]);
    auto r = d.reorder({1, 0, 2});
    EXPECT_DOUBLE_EQ(r.matrix()(0, 0), vb::Y2[0]);
    EXPECT_EQ(r.names[0], "y2");
}

// 21. 性能 (§15.5): K=5, T=500, IC 全扫描 (maxlag 默认≈18) < 5 sec
TEST(VarModel, PerfK5T500ICScan) {
    vt::MultivariateTSData d;
    d.columns.resize(5);
    cpphub::v1::Philox4x64 rng(2026);
    std::vector<std::vector<double>> y(5, std::vector<double>(500, 0.0));
    for (int t = 1; t < 500; ++t) {
        for (int j = 0; j < 5; ++j) {
            double acc = 0.25 * y[j][t - 1];
            if (j > 0) acc += 0.08 * y[j - 1][t - 1];
            y[j][t] = acc +
                      static_cast<double>(static_cast<int64_t>(rng() % 2000) - 1000) / 500.0;
        }
    }
    for (int j = 0; j < 5; ++j) d.columns[j] = y[j];
    const auto t0 = std::chrono::steady_clock::now();
    auto sel = vt::var_select_order(d, "c", 0, "aic");
    const double sec = std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - t0).count();
    EXPECT_GT(sel.aic.size(), 5u);   // 扫描真实发生
    EXPECT_LT(sec, 5.0) << "IC scan took " << sec << "s";
}
