// =============================================================================
// test_var_irf_fevd.cpp - IRF/FEVD 双轨/GFEVD 框架 (Phase 7C M2, 18 用例)
//
// 基准: var_baseline.inc (statsmodels 0.14.4 orth_ma_rep/fevd 1e-10~1e-12 +
//       R Spillover 0.1.1 g.fevd 1e-8 + vars fevd 交叉)
// 幻觉点: V3 (Θ 行=响应列=冲击) / V7 (行归一) / V8 (DY σ_jj⁻¹ vs PS σ_ii⁻¹)
//       / V12 (不稳定拦截) / V13 (bootstrap 带不进容差)
// =============================================================================

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/timeseries/var/fevd.hpp"
#include "cpphub/timeseries/var/irf.hpp"
#include "cpphub/timeseries/var/multivariate_data.hpp"
#include "cpphub/timeseries/var/var_model.hpp"
#include "var_baseline.inc"

namespace vt = cpphub::v1::timeseries::var;
namespace vb = cpphub::v1::timeseries::var_baseline;
using Eigen::MatrixXd;

namespace {

vt::MultivariateTSData fixture_data() {
    vt::MultivariateTSData d;
    d.columns = {std::vector<cpphub::v1::Real>(std::begin(vb::Y1), std::end(vb::Y1)),
                 std::vector<cpphub::v1::Real>(std::begin(vb::Y2), std::end(vb::Y2)),
                 std::vector<cpphub::v1::Real>(std::begin(vb::Y3), std::end(vb::Y3))};
    return d;
}

vt::VARResult fit2() {
    vt::VARSpec spec;
    spec.lag = 2;
    return vt::var_fit(fixture_data(), spec);
}

MatrixXd mat3(const double* flat) {
    MatrixXd m(3, 3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) m(i, j) = flat[3 * i + j];
    return m;
}

}  // namespace

// 1. Φ_0 = I (精确)
TEST(VarIrfFevd, PhiH0Identity) {
    auto fit = fit2();
    auto irf = vt::var_irf(fit, 3);
    ASSERT_EQ(irf.phi.size(), 3u);
    EXPECT_DOUBLE_EQ(irf.phi[0](0, 0), 1.0);
    EXPECT_DOUBLE_EQ(irf.phi[0](1, 1), 1.0);
    EXPECT_DOUBLE_EQ(irf.phi[0](0, 1), 0.0);
}

// 2. Φ_1 vs statsmodels ma_rep (未正交化), 1e-12
TEST(VarIrfFevd, PhiH1VsStatsmodels) {
    auto fit = fit2();
    auto irf = vt::var_irf(fit, 3);
    const MatrixXd ref = mat3(vb::SM_PHI_H1.data());
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(irf.phi[1](i, j), ref(i, j), 1e-12) << i << j;
    // Φ_1 = A_1 (伴随递推退化)
    EXPECT_NEAR(irf.phi[1](0, 0), fit.coefficients(0, 1), 1e-14);
}

// 2b. Φ 递推恒等式: Φ_2 = Φ_1·A_1 + A_2 (Lütkepohl Ch.2, p=2 情形)
TEST(VarIrfFevd, PhiRecursionIdentity) {
    auto fit = fit2();
    auto irf = vt::var_irf(fit, 3);
    const auto A = vt::detail::var_companion_blocks(fit);
    const MatrixXd expect = irf.phi[1] * A[0] + A[1];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(irf.phi[2](i, j), expect(i, j), 1e-13) << i << j;
}

// 2c. horizon 边界: IRF/FEVD horizon=0 均拒绝
TEST(VarIrfFevd, HorizonZeroThrows) {
    auto fit = fit2();
    EXPECT_THROW(vt::var_irf(fit, 0), std::invalid_argument);
    EXPECT_THROW(vt::var_fevd(fit, 0), std::invalid_argument);
}

// 3. Ψ_0 = P 下三角 (V3: 正交化零期 = Cholesky 因子)
TEST(VarIrfFevd, OrthIrfH0IsCholesky) {
    auto fit = fit2();
    auto irf = vt::var_irf(fit, 2);
    const MatrixXd ref = mat3(vb::SM_IRF_H0.data());
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(irf.theta[0](i, j), ref(i, j), 1e-12) << i << j;
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) EXPECT_DOUBLE_EQ(irf.theta[0](i, j), 0.0);
}

// 4-6. Ψ_h vs statsmodels orth_ma_rep (V3), h=1/2/10, 1e-12
TEST(VarIrfFevd, OrthIrfVsStatsmodels) {
    auto fit = fit2();
    auto irf = vt::var_irf(fit, 11);
    ASSERT_EQ(irf.theta.size(), 11u);
    const struct { int h; const std::array<double, 9>& ref; } cases[] = {
        {1, vb::SM_IRF_H1}, {2, vb::SM_IRF_H2}, {10, vb::SM_IRF_H10}};
    for (const auto& c : cases) {
        const MatrixXd ref = mat3(c.ref.data());
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                EXPECT_NEAR(irf.theta[c.h](i, j), ref(i, j), 1e-12)
                    << "h=" << c.h << " " << i << j;
    }
}

// 7. V3 方向: Θ[i,j] 行=响应 列=冲击 (非对称断言 + 与 P 列一致性)
TEST(VarIrfFevd, IRFDirection) {
    auto fit = fit2();
    auto irf = vt::var_irf(fit, 2);
    // Ψ_1(1,0) = Σ_k Φ_1(1,k)·P(k,0) (P 下三角 ⇒ k≤0 全列参与)
    double expect = 0.0;
    for (int k = 0; k < 3; ++k)
        expect += fit.coefficients(1, 1 + k) * irf.P(k, 0);
    EXPECT_NEAR(irf.theta[1](1, 0), expect, 1e-13);
    // Θ 非对称 (若转置则对角块错)
    EXPECT_NE(irf.theta[1](0, 2), irf.theta[1](2, 0));
}

// 8. P 注入 (决策 11): 下三角注入生效; 非下三角拒绝
TEST(VarIrfFevd, PInjection) {
    auto fit = fit2();
    MatrixXd Pinj = MatrixXd::Zero(3, 3);
    Pinj(0, 0) = 1.0;
    Pinj(1, 1) = 1.0;
    Pinj(2, 2) = 1.0;
    Pinj(1, 0) = 0.5;
    vt::VARSpec spec;
    spec.lag = 2;
    spec.identification_P = Pinj;
    // 注入通过 detail::orthogonal_P (由 var_irf 内部走默认 — 注入路径经
    // fevd/irf 公共接口暂为 Cholesky; 此处直接测 detail 层)
    const MatrixXd P = vt::detail::orthogonal_P(fit, Pinj);
    EXPECT_DOUBLE_EQ(P(1, 0), 0.5);
    EXPECT_DOUBLE_EQ(P(2, 2), 1.0);
    MatrixXd bad = MatrixXd::Identity(3, 3);
    bad(0, 2) = 0.3;  // 上三角非零
    EXPECT_THROW(vt::detail::orthogonal_P(fit, bad), std::invalid_argument);
    // 空 (0×0) → Cholesky 默认
    const MatrixXd Pd = vt::detail::orthogonal_P(fit, MatrixXd{});
    EXPECT_NEAR(Pd(0, 0), vb::SM_CHOL_LOWER[0], 1e-12);
}

// 9. FEVD Cholesky H=10 vs statsmodels (1e-10)
TEST(VarIrfFevd, FEVDCholeskyH10) {
    auto fit = fit2();
    auto f = vt::var_fevd(fit, 10, vt::FevdFramework::Cholesky);
    const MatrixXd ref = mat3(vb::SM_FEVD_H10.data());
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(f.fevd(i, j), ref(i, j), 1e-10) << i << j;
    EXPECT_EQ(f.framework, vt::FevdFramework::Cholesky);
    EXPECT_EQ(f.horizon, 10u);
}

// 10. FEVD Cholesky H=1 vs statsmodels (1e-10)
TEST(VarIrfFevd, FEVDCholeskyH1) {
    auto fit = fit2();
    auto f = vt::var_fevd(fit, 1, vt::FevdFramework::Cholesky);
    const MatrixXd ref = mat3(vb::SM_FEVD_H1.data());
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(f.fevd(i, j), ref(i, j), 1e-10) << i << j;
}

// 11. V7: 行和精确 = 1 (Cholesky + DY), 1e-12
TEST(VarIrfFevd, FEVDRowSumOne) {
    auto fit = fit2();
    for (auto fw : {vt::FevdFramework::Cholesky, vt::FevdFramework::GeneralizedDY}) {
        auto f = vt::var_fevd(fit, 10, fw);
        for (int i = 0; i < 3; ++i)
            EXPECT_NEAR(f.fevd.row(i).sum(), 1.0, 1e-12) << static_cast<int>(fw);
    }
}

// 12. R vars fevd 交叉 (1e-10); R dump 列主序 → fevd(i,j)=inc[3j+i]
TEST(VarIrfFevd, FEVDRVarsCross) {
    auto fit = fit2();
    auto f = vt::var_fevd(fit, 10, vt::FevdFramework::Cholesky);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(f.fevd(i, j), vb::R_VARS_FEVD_H10[3 * j + i], 1e-10)
                << i << j;
}

// 13. GFEVD DY vs Spillover g.fevd (主基准 1e-8)
TEST(VarIrfFevd, GFEVDDYVsSpillover) {
    auto fit = fit2();
    auto f = vt::var_fevd(fit, 10, vt::FevdFramework::GeneralizedDY);
    const MatrixXd ref = mat3(vb::SP_GFEVD_H10.data());
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(f.fevd(i, j), ref(i, j), 1e-8) << i << j;
    EXPECT_EQ(f.framework, vt::FevdFramework::GeneralizedDY);
}

// 14. 未归一 GFEVD 结构: raw 行归一 == DY 归一版 (R Spillover raw/normalized
//     行比值恒等: c = [1.1122, 1.2462, 1.1468] > 1, 分母为加权总方差)
TEST(VarIrfFevd, GFEVDRawStructure) {
    auto fit = fit2();
    auto f = vt::var_fevd(fit, 10, vt::FevdFramework::GeneralizedDY);
    const MatrixXd raw = mat3(vb::SP_GFEVD_H10_RAW.data());
    for (int i = 0; i < 3; ++i) {
        const double row_sum = raw.row(i).sum();
        EXPECT_GT(row_sum, 1.0);   // R raw 行和 > 1 (归一分母 < 行和)
        EXPECT_LT(row_sum, 2.0);
        // normalized = raw / row_sum
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(f.fevd(i, j), raw(i, j) / row_sum, 1e-9) << i << j;
    }
}

// 15. V8: PS 框架 ≠ DY 框架 (σ 不等时); PS 未归一 (行和 ≠ 1)
TEST(VarIrfFevd, PSvsDYFrameworkDifference) {
    auto fit = fit2();
    auto fdy = vt::var_fevd(fit, 10, vt::FevdFramework::GeneralizedDY);
    auto fps = vt::var_fevd(fit, 10, vt::FevdFramework::GeneralizedPS);
    // PS 行和 ≠ 1 (未归一)
    for (int i = 0; i < 3; ++i)
        EXPECT_NE(fps.fevd.row(i).sum(), 1.0);
    // PS 归一化后与 DY 不同 (Σ 对角不等: 0.868/1.420/0.848)
    double max_diff = 0.0;
    for (int i = 0; i < 3; ++i) {
        const double rs = fps.fevd.row(i).sum();
        for (int j = 0; j < 3; ++j)
            max_diff = std::max(max_diff,
                                std::abs(fps.fevd(i, j) / rs - fdy.fevd(i, j)));
    }
    EXPECT_GT(max_diff, 1e-3);  // 数值显著不同
    EXPECT_EQ(fps.framework, vt::FevdFramework::GeneralizedPS);
}

// 16. 数学恒等: raw GFEVD 行归一 == DY 归一版 (逐行, 非仅首行; Cholesky
//     FEVD 与 GFEVD 首行无此恒等 — 分母定义不同, 不可混用)
TEST(VarIrfFevd, RawGFEVDRowNormalizationIdentity) {
    auto fit = fit2();
    auto fdy = vt::var_fevd(fit, 10, vt::FevdFramework::GeneralizedDY);
    const MatrixXd raw = mat3(vb::SP_GFEVD_H10_RAW.data());
    for (int i = 0; i < 3; ++i) {
        const double row_sum = raw.row(i).sum();
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(fdy.fevd(i, j), raw(i, j) / row_sum, 1e-9) << i << j;
    }
}

// 17. V12: 不稳定 VAR 拦截 FEVD (爆炸性 ρ=1.05 ⇒ 特征值模 >1 确定性;
//     随机游走有限样本 OLS 常得 ρ̂<1, 不能稳定触发拦截)
TEST(VarIrfFevd, UnstableVarThrows) {
    // 构造爆炸 VAR(1): y_t = 1.05·y_{t-1} + e
    vt::MultivariateTSData d;
    d.columns.resize(2);
    cpphub::v1::Philox4x64 rng(7);
    double a = 0.0, b = 0.0;
    for (int t = 0; t < 200; ++t) {
        a = 1.05 * a + static_cast<double>(static_cast<int64_t>(rng() % 2000) - 1000) / 1000.0;
        b = 1.05 * b + static_cast<double>(static_cast<int64_t>(rng() % 2000) - 1000) / 1000.0;
        d.columns[0].push_back(a);
        d.columns[1].push_back(b);
    }
    vt::VARSpec spec;
    spec.lag = 1;
    auto fit = vt::var_fit(d, spec);
    ASSERT_GT(fit.max_abs_eigenvalue, 1.0);
    ASSERT_FALSE(fit.is_strictly_stationary);
    EXPECT_THROW(vt::var_fevd(fit, 10, vt::FevdFramework::GeneralizedDY),
                 std::invalid_argument);
    EXPECT_THROW(vt::var_fevd(fit, 10, vt::FevdFramework::Cholesky),
                 std::invalid_argument);
}

// 18. V13: bootstrap 带存在且仅描述不确定性 (带不进容差断言 — 结构检查)
TEST(VarIrfFevd, BootstrapBands) {
    auto fit = fit2();
    auto data = fixture_data();
    auto irf = vt::var_irf_bootstrap(data, fit, 3, 200, 42);
    EXPECT_TRUE(irf.has_bands);
    ASSERT_EQ(irf.irf_lower.size(), 3u);
    ASSERT_EQ(irf.irf_upper.size(), 3u);
    // 结构: lower <= upper 且均为有限值
    for (int h = 0; h < 3; ++h)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                EXPECT_LE(irf.irf_lower[h](i, j), irf.irf_upper[h](i, j) + 1e-12);
                EXPECT_TRUE(std::isfinite(irf.irf_lower[h](i, j)));
                EXPECT_TRUE(std::isfinite(irf.irf_upper[h](i, j)));
            }
    // h=0 对角: P 对角为正, 带应包含正区间 (弱结构约束)
    for (int i = 0; i < 3; ++i) EXPECT_GT(irf.irf_upper[0](i, i), 0.0);
    // 无数据重载: NaN 占位
    auto plain = vt::var_irf(fit, 3, /*bootstrap=*/true, 100, 42);
    EXPECT_FALSE(plain.has_bands);
    EXPECT_TRUE(std::isnan(plain.irf_lower[0](0, 0)));
}
