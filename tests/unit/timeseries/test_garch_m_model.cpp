// =============================================================================
// test_garch_m_model.cpp - GARCH(1,1)-M 测试 (16 用例, spec §1.2/§9.5)
//
// 基准: tests/unit/timeseries/gm_baseline.inc (verify_gm.py 自动生成:
//       arch 8.0.0 ARCHInMean, rescale=False, Normal, 三 form)
//
// 容差 (7B test_garch_model.cpp 头注先例, 可达精度):
//   - params/llf: 1e-4 (实测 1e-5~1e-6; scipy SLSQP 解析梯度 vs 自研数值
//     梯度落点差; spec 名义 1e-8~1e-10 仅同优化器同梯度模式可达)
//   - robust SE: 5e-2 (实测 ≤2e-2; 数值 Hessian/OPG 步长舍入噪声)
//   - 确定性函数 (gm_g/filter/等变映射): 1e-12~1e-15
//
// 幻觉点覆盖 (spec §9.5):
//   GM1 form 映射 (5/6/7/11) / GM2 arch 8.0 才有 (基准即证) / GM3 log 有基准
//   (3) / GM4 sandwich SE (4/13) / GM5 耦合递归 (8/9/10)
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/garch/garch_m_model.hpp"
#include "gm_baseline.inc"

namespace ts = cpphub::v1::timeseries;
using ts::garch::estimate_garch_m;
using ts::garch::GarchConfig;
using ts::garch::GarchDist;
using ts::garch::GarchMForm;
using cpphub::Real;
using cpphub::Size;

static std::vector<Real> y_gm() {
    return {ts::garch::gm_baseline::Y_GM,
            ts::garch::gm_baseline::Y_GM + ts::garch::gm_baseline::T};
}

// 三 form 参数/llf 对照 (容差 1e-4, 名义可达精度)
static void check_form(const ts::garch::GarchMResult& r,
                       const ts::garch::gm_baseline::GmCase& c) {
    EXPECT_NEAR(r.params.mu, c.mu, 1e-4);
    EXPECT_NEAR(r.params.lambda, c.lambda, 1e-4);
    EXPECT_NEAR(r.params.omega, c.omega, 1e-4);
    EXPECT_NEAR(r.params.alpha, c.alpha, 1e-4);
    EXPECT_NEAR(r.params.beta, c.beta, 1e-4);
    EXPECT_NEAR(r.log_likelihood, c.llf, 1e-5);
    EXPECT_TRUE(r.converged);
}

// ---------------------------------------------------------------------------
// 1-3: arch 基准三 form (params/llf, GM2/GM3 — arch 8.0 ARCHInMean 存在性即证)
// ---------------------------------------------------------------------------
TEST(GarchMTest, ArchBaselineVol) {
    check_form(estimate_garch_m(y_gm(), GarchMForm::Volatility),
               ts::garch::gm_baseline::ARCH_VOL);
}

TEST(GarchMTest, ArchBaselineVar) {
    check_form(estimate_garch_m(y_gm(), GarchMForm::Variance),
               ts::garch::gm_baseline::ARCH_VAR);
}

TEST(GarchMTest, ArchBaselineLog) {
    check_form(estimate_garch_m(y_gm(), GarchMForm::LogVariance),
               ts::garch::gm_baseline::ARCH_LOG);
}

// ---------------------------------------------------------------------------
// 4: GM4 — robust SE (BW QMLE) vs arch std_err (5e-2, 数值三明治噪声层)
// ---------------------------------------------------------------------------
TEST(GarchMTest, ArchBaselineRobustSe) {
    struct Case {
        GarchMForm form;
        const ts::garch::gm_baseline::GmCase* ref;
    };
    for (const auto& c : std::vector<Case>{
             {GarchMForm::Volatility, &ts::garch::gm_baseline::ARCH_VOL},
             {GarchMForm::Variance, &ts::garch::gm_baseline::ARCH_VAR},
             {GarchMForm::LogVariance, &ts::garch::gm_baseline::ARCH_LOG}}) {
        const auto r = estimate_garch_m(y_gm(), c.form);
        ASSERT_EQ(r.std_errors.size(), 5u);
        EXPECT_GT(r.std_errors[0], 0.0);
        EXPECT_GT(r.std_errors[1], 0.0);
        EXPECT_NEAR(r.std_errors[0], c.ref->se_mu, 5e-2);
        EXPECT_NEAR(r.std_errors[1], c.ref->se_lambda, 5e-2);
    }
}

// ---------------------------------------------------------------------------
// 5: GM1 — g(h) 形式函数精确值
// ---------------------------------------------------------------------------
TEST(GarchMTest, GmFormFunctions) {
    EXPECT_DOUBLE_EQ(ts::garch::detail::gm_g(0.25, GarchMForm::Variance), 0.25);
    EXPECT_DOUBLE_EQ(ts::garch::detail::gm_g(0.25, GarchMForm::Volatility), 0.5);
    EXPECT_DOUBLE_EQ(ts::garch::detail::gm_g(0.25, GarchMForm::LogVariance),
                     std::log(0.25));
}

// ---------------------------------------------------------------------------
// 6: GM1 — λ 等变映射往返恒等 + Jacobian 值 (vol:1 / var:s / log:1/s)
// ---------------------------------------------------------------------------
TEST(GarchMTest, LambdaScaleMappingRoundtrip) {
    const Real s = 3.7;
    const Real lam = -0.42;
    for (const auto form : {GarchMForm::Variance, GarchMForm::Volatility,
                            GarchMForm::LogVariance}) {
        const Real ls = ts::garch::detail::gm_lambda_to_scaled(lam, s, form);
        EXPECT_NEAR(ts::garch::detail::gm_lambda_from_scaled(ls, s, form), lam,
                    1e-15);
    }
    EXPECT_DOUBLE_EQ(ts::garch::detail::gm_lambda_jacobian(s, GarchMForm::Variance), s);
    EXPECT_DOUBLE_EQ(ts::garch::detail::gm_lambda_jacobian(s, GarchMForm::Volatility), 1.0);
    EXPECT_DOUBLE_EQ(ts::garch::detail::gm_lambda_jacobian(s, GarchMForm::LogVariance),
                     1.0 / s);
}

// ---------------------------------------------------------------------------
// 7: log 形 μ 映射交叉项往返 (μ' = s·μ − 2sλ·lns)
// ---------------------------------------------------------------------------
TEST(GarchMTest, MuLogFormCrossTerm) {
    const Real s = 5.0, mu = 0.11, lam = -0.3;
    const Real mu_s = ts::garch::detail::gm_mu_to_scaled(mu, lam, s,
                                                         GarchMForm::LogVariance);
    const Real lam_s = ts::garch::detail::gm_lambda_to_scaled(
        lam, s, GarchMForm::LogVariance);
    EXPECT_NEAR(ts::garch::detail::gm_mu_from_scaled(mu_s, lam_s, s,
                                                     GarchMForm::LogVariance),
                mu, 1e-15);
    // vol/var 无交叉项: 纯线性
    EXPECT_DOUBLE_EQ(
        ts::garch::detail::gm_mu_to_scaled(mu, lam, s, GarchMForm::Volatility),
        mu * s);
}

// ---------------------------------------------------------------------------
// 8: GM5 退化性 — λ=0 时耦合递归 ≡ 7B filter_garch11 (h 逐位, eps = y − μ)
// ---------------------------------------------------------------------------
TEST(GarchMTest, DegenerateLambdaZeroMatchesGarch11) {
    const auto y = y_gm();
    const Size T = y.size();
    const Real mu = 0.05, om = 0.05, al = 0.10, be = 0.85;
    std::vector<Real> eps0(T);
    for (Size t = 0; t < T; ++t) eps0[t] = y[t] - mu;
    const Real s0 = ts::garch::backcast_variance(eps0, 0.94, 0);
    const auto p = ts::garch::detail::filter_garch_m(
        GarchMForm::Volatility, mu, 0.0, om, al, be, y, s0);
    const ts::garch::GarchParams gp{om, al, be, 0.0};
    const auto h11 = ts::garch::filter_garch11(gp, eps0, s0);
    ASSERT_EQ(p.h.size(), T);
    for (Size t = 0; t < T; ++t) {
        EXPECT_NEAR(p.h[t], h11[t], 1e-15);
        EXPECT_NEAR(p.eps[t], eps0[t], 1e-15);
    }
}

// ---------------------------------------------------------------------------
// 9: GM5 (issue #269) — λ 扰动重写整条 h 路径: h[0] 不变, t ≥ 1 全变
// ---------------------------------------------------------------------------
TEST(GarchMTest, CouplingLambdaPerturbsPath) {
    const auto y = y_gm();
    const Size T = y.size();
    std::vector<Real> eps0(T);
    for (Size t = 0; t < T; ++t) eps0[t] = y[t] - 0.05;
    const Real s0 = ts::garch::backcast_variance(eps0, 0.94, 0);
    const auto p0 = ts::garch::detail::filter_garch_m(
        GarchMForm::Volatility, 0.05, 0.5, 0.05, 0.10, 0.85, y, s0);
    const auto p1 = ts::garch::detail::filter_garch_m(
        GarchMForm::Volatility, 0.05, 0.51, 0.05, 0.10, 0.85, y, s0);
    EXPECT_DOUBLE_EQ(p0.h[0], p1.h[0]);  // h₁ = f(backcast) 与 λ 无关
    Real maxdh = 0.0;
    for (Size t = 1; t < T; ++t) {
        maxdh = std::max(maxdh, std::fabs(p0.h[t] - p1.h[t]));
    }
    EXPECT_GT(maxdh, 1e-6);  // 耦合: λ+0.01 传播到全路径
}

// ---------------------------------------------------------------------------
// 10: GM5 — 耦合 vs 无耦合反例: 用 λ=0 残差驱动的 GARCH 路径 ≠ 耦合路径
// ---------------------------------------------------------------------------
TEST(GarchMTest, CoupledVsUncoupledDivergence) {
    const auto y = y_gm();
    const Size T = y.size();
    const Real mu = 0.05, lam = 0.5, om = 0.05, al = 0.10, be = 0.85;
    std::vector<Real> eps0(T);
    for (Size t = 0; t < T; ++t) eps0[t] = y[t] - mu;
    const Real s0 = ts::garch::backcast_variance(eps0, 0.94, 0);
    // 两步解耦 (错): 固定残差 eps0 驱动方差
    const auto h_unc = ts::garch::filter_garch11(ts::garch::GarchParams{om, al, be, 0.0},
                                                 eps0, s0);
    // 耦合 (对): ε 含 λ·√h 反馈
    const auto p = ts::garch::detail::filter_garch_m(
        GarchMForm::Volatility, mu, lam, om, al, be, y, s0);
    Real maxd = 0.0;
    for (Size t = 0; t < T; ++t) {
        maxd = std::max(maxd, std::fabs(p.h[t] - h_unc[t]));
    }
    EXPECT_GT(maxd, 1e-3);  // 显著分离 (GM5: 不可两步解耦)
}

// ---------------------------------------------------------------------------
// 11: 等变性实测 — 数据×2: vol λ 不变 / var λ 减半 / μ×2 / ω×4
// ---------------------------------------------------------------------------
TEST(GarchMTest, ScaleInvarianceEmpirical) {
    const auto y = y_gm();
    std::vector<Real> y2(y.size());
    for (Size t = 0; t < y.size(); ++t) y2[t] = 2.0 * y[t];
    const auto r1v = estimate_garch_m(y, GarchMForm::Volatility);
    const auto r2v = estimate_garch_m(y2, GarchMForm::Volatility);
    EXPECT_NEAR(r2v.params.lambda, r1v.params.lambda, 1e-6);
    EXPECT_NEAR(r2v.params.mu, 2.0 * r1v.params.mu, 1e-6);
    EXPECT_NEAR(r2v.params.omega, 4.0 * r1v.params.omega, 1e-5);
    const auto r1a = estimate_garch_m(y, GarchMForm::Variance);
    const auto r2a = estimate_garch_m(y2, GarchMForm::Variance);
    // 容差 1e-5: λ' 边界随 s 变化 (lam_hi = 10/s) 致落点差 ~2e-6 (实测)
    EXPECT_NEAR(r2a.params.lambda, 0.5 * r1a.params.lambda, 1e-5);
}

// ---------------------------------------------------------------------------
// 12: 默认 form = Volatility (arch 默认 'vol' 对齐) + 收敛
// ---------------------------------------------------------------------------
TEST(GarchMTest, DefaultFormIsVolatility) {
    const auto rd = estimate_garch_m(y_gm());
    const auto rv = estimate_garch_m(y_gm(), GarchMForm::Volatility);
    EXPECT_DOUBLE_EQ(rd.params.lambda, rv.params.lambda);
    EXPECT_DOUBLE_EQ(rd.params.omega, rv.params.omega);
    EXPECT_TRUE(rd.converged);
    EXPECT_GT(rd.n_iterations, 0u);
}

// ---------------------------------------------------------------------------
// 13: GM4 — vcov 形状/对称性 (log 形 J 交叉项不破坏对称)
// ---------------------------------------------------------------------------
TEST(GarchMTest, SandwichShapeSymmetry) {
    for (const auto form : {GarchMForm::Variance, GarchMForm::Volatility,
                            GarchMForm::LogVariance}) {
        const auto r = estimate_garch_m(y_gm(), form);
        ASSERT_EQ(r.vcov.size(), 5u);
        for (Size a = 0; a < 5; ++a) {
            ASSERT_EQ(r.vcov[a].size(), 5u);
            EXPECT_GT(r.vcov[a][a], 0.0);
            for (Size b = 0; b < 5; ++b) {
                EXPECT_NEAR(r.vcov[a][b], r.vcov[b][a], 1e-10);
            }
        }
        for (Size a = 0; a < 5; ++a) {
            EXPECT_NEAR(r.std_errors[a], std::sqrt(r.vcov[a][a]), 1e-15);
        }
    }
}

// ---------------------------------------------------------------------------
// 14: 恒等式 — z_t = ε_t/√h_t + AIC/BIC 公式 (k=5)
// ---------------------------------------------------------------------------
TEST(GarchMTest, StdResidualsAndIcIdentity) {
    const auto r = estimate_garch_m(y_gm(), GarchMForm::Variance);
    const Size T = ts::garch::gm_baseline::T;
    ASSERT_EQ(r.std_residuals.size(), T);
    ASSERT_EQ(r.conditional_variances.size(), T);
    Real zmean = 0.0;
    for (Size t = 0; t < T; ++t) {
        EXPECT_NEAR(r.std_residuals[t],
                    r.residuals[t] / std::sqrt(r.conditional_variances[t]), 1e-12);
        zmean += r.std_residuals[t];
    }
    EXPECT_LT(std::fabs(zmean / static_cast<Real>(T)), 0.1);
    EXPECT_NEAR(r.aic, -2.0 * r.log_likelihood + 10.0, 1e-9);
    EXPECT_NEAR(r.bic,
                -2.0 * r.log_likelihood
                    + 5.0 * std::log(static_cast<Real>(T)),
                1e-9);
}

// ---------------------------------------------------------------------------
// 15: t 分布 pass-through (ν 联合估计收敛, vcov 6×6)
// ---------------------------------------------------------------------------
TEST(GarchMTest, StudentTPassThrough) {
    GarchConfig cfg;
    cfg.dist = GarchDist::StudentT;
    const auto r = estimate_garch_m(y_gm(), GarchMForm::Volatility, cfg);
    EXPECT_TRUE(r.converged);
    EXPECT_GT(r.params.omega, 0.0);
    // spec 冻结 GarchMResult 无 ν 字段 (7B GarchResult 有, GM 不回显);
    // ν 进入似然与 sandwich 维度: vcov/std_errors 6 维
    ASSERT_EQ(r.vcov.size(), 6u);
    ASSERT_EQ(r.std_errors.size(), 6u);
}

// ---------------------------------------------------------------------------
// 16: 异常输入 (T<10 / NaN / 常数列)
// ---------------------------------------------------------------------------
TEST(GarchMTest, InvalidInputThrows) {
    EXPECT_THROW(estimate_garch_m(std::vector<Real>(9, 0.01)),
                 std::invalid_argument);
    auto y = y_gm();
    y[7] = std::nan("");
    EXPECT_THROW(estimate_garch_m(y), std::invalid_argument);
    EXPECT_THROW(estimate_garch_m(std::vector<Real>(200, 3.14)),
                 std::invalid_argument);
}
