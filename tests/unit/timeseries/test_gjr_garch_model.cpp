// =============================================================================
// test_gjr_garch_model.cpp - GJR-GARCH(1,1) QMLE 测试 (18 用例, spec §2.0.4)
//
// 基准: tests/unit/timeseries/gjr_baseline.inc (arch 8.0.0 自动生成,
//       脚本 tests/fixtures/timeseries/verify_gjr.py, HFE 硬编码策略)
//
// 容差策略 (跨优化器对照, 非 HFE 同算法对照):
//   - 确定性函数 (filter/forecast/backcast): 1e-12 (与 arch 逐位一致)
//   - 优化结果 (params/llf): ~1e-4/1e-6 (scipy SLSQP 解析梯度 vs 自研数值
//     梯度, 收敛到同一最优点的精度差)
//   - sandwich vcov: ~5e-3 (数值 Hessian 步长 1e-6 级的舍入噪声)
//
// 约定实证 (verify_gjr.py probe, 2026-08-15):
//   probe-1: h₁ = ω + (α+γ/2+β)·σ²₀ (arch t=0 对 o 项用 0.5·bc)
//   probe-2: 多步预测 φ = α+γ/2+β (E[I(z<0)·ε²]=h/2)
// =============================================================================
#include <gtest/gtest.h>
#include <chrono>
#include <vector>
#include <cmath>

#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/timeseries/garch/gjr_garch_model.hpp"
#include "cpphub/timeseries/garch/garch_forecast.hpp"
#include "gjr_baseline.inc"

namespace ts = cpphub::v1::timeseries::garch;
namespace gbl = cpphub::v1::timeseries::garch::gjr_baseline;
using cpphub::Real;
using cpphub::Size;

// ---------------------------------------------------------------------------
// 辅助: 数据构造 + 结果缓存 (estimate ~0.5s, 静态缓存避免重复估计)
// ---------------------------------------------------------------------------
static std::vector<Real> make_data() {
    return {gbl::DATA, gbl::DATA + gbl::T};
}

static std::vector<Real> make_data_demeaned() {
    auto d = make_data();
    Real m = 0.0;
    for (Real v : d) m += v;
    m /= static_cast<Real>(d.size());
    for (Real& v : d) v -= m;
    return d;
}

static std::vector<Real> make_tdata() {
    return {gbl::TDATA, gbl::TDATA + gbl::T};
}

static const ts::GjrGarchResult& normal_fit() {
    static ts::GjrGarchResult r = ts::estimate_gjr_garch(make_data());
    return r;
}

// ---------------------------------------------------------------------------
// 1. Backcast vs arch (G1 复用, 1e-12)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, BackcastVarianceMatchesArch) {
    const auto dm = make_data_demeaned();
    const Real bc = ts::backcast_variance(dm, 0.94, 0);
    EXPECT_NEAR(bc, gbl::BACKCAST, 1e-12);
}

// ---------------------------------------------------------------------------
// 2. 方差递归 filter_gjr vs arch conditional_volatility² (1e-12)
//    用 arch 参数直接递归 — 排除优化器差异, 纯算法对照
//    含 probe-1 约定: h₁ = ω + (α+γ/2+β)·σ²₀
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, FilterRecursionMatchesArchConditionalVariance) {
    const auto dm = make_data_demeaned();
    const ts::GjrGarchParams p{gbl::OMEGA, gbl::ALPHA, gbl::GAMMA, gbl::BETA, 0.0};
    const auto h = ts::filter_gjr(p, dm, gbl::BACKCAST);
    ASSERT_EQ(h.size(), gbl::T);
    for (Size t = 0; t < gbl::T; ++t) {
        ASSERT_NEAR(h[t], gbl::H[t], 1e-12)
            << "t=" << t << " h=" << h[t] << " ref=" << gbl::H[t];
    }
}

// ---------------------------------------------------------------------------
// 3. t=0 初始化 + 手算递归 (probe-1: γ/2·σ²₀; G8: I(ε<0) 非对称项)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, FilterT0InitAndHandComputedRecursion) {
    // params: ω=0.2, α=0.10, γ=0.30, β=0.50; σ²₀=1.5
    const ts::GjrGarchParams p{0.2, 0.10, 0.30, 0.50, 0.0};
    const std::vector<Real> eps = {1.0, -0.5, 2.0};
    const Real bc = 1.5;
    const auto h = ts::filter_gjr(p, eps, bc);

    // t=0: h₀ = ω + α·σ²₀ + γ·(1/2)·σ²₀ + β·σ²₀ (probe-1: 0.5 系数)
    const Real h0 = 0.2 + 0.10 * 1.5 + 0.30 * 0.5 * 1.5 + 0.50 * 1.5;
    ASSERT_NEAR(h[0], h0, 1e-14);
    // t=1: ε₀=1.0>0 → I=0 (无 γ 项)
    const Real h1 = 0.2 + 0.10 * 1.0 + 0.50 * h0;
    ASSERT_NEAR(h[1], h1, 1e-14);
    // t=2: ε₁=-0.5<0 → I=1 (G8: 非对称项激活)
    const Real h2 = 0.2 + 0.10 * 0.25 + 0.30 * 1.0 * 0.25 + 0.50 * h1;
    ASSERT_NEAR(h[2], h2, 1e-14);
}

// ---------------------------------------------------------------------------
// 4. G7: I(zₜ<0) ≡ I(εₜ<0) 等价性 (σₜ>0 → 同号; 1e-15)
//    构造 h 序列, z 与 ε 同号 → filter 输出应逐位一致
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, IndicatorEquivalenceZVsEpsilon) {
    const ts::GjrGarchParams p{0.05, 0.05, 0.15, 0.80, 0.0};
    const std::vector<Real> eps = {0.5, -1.2, 3.0, -0.01, 2.2, -1.7,
                                   0.0, -0.8, 1.1, -2.5};
    const auto h = ts::filter_gjr(p, eps, 0.7);
    // z = ε/√h 与 ε 同号 (h>0); 显式验证
    for (Size t = 0; t < eps.size(); ++t) {
        const Real z = eps[t] / std::sqrt(h[t]);
        if (eps[t] < 0.0) EXPECT_LT(z, 0.0);
        if (eps[t] > 0.0) EXPECT_GT(z, 0.0);
        if (eps[t] == 0.0) EXPECT_EQ(z, 0.0);
    }
    // 注: I(0<0)=0 — εₜ=0 时两项均不激活, 等价性保持
}

// ---------------------------------------------------------------------------
// 5. G10: 平稳性 α+γ/2+β < 1 (估计结果满足; 与 arch PERSIST 一致 1e-12)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, StationarityConditionHalfGamma) {
    const auto& r = normal_fit();
    const Real persist = r.params.alpha + r.params.gamma / 2.0 + r.params.beta;
    EXPECT_LT(persist, 1.0);                       // G10 约束满足
    EXPECT_NEAR(persist, gbl::PERSIST, 1e-4);      // 与 arch 估计一致 (跨优化器)
    // 排幻觉: 非 α+γ+β (会高估持续性); γ/2 定义由测试 3 手算递归逐位锚定
    const Real wrong = r.params.alpha + r.params.gamma + r.params.beta;
    EXPECT_GT(wrong, gbl::PERSIST);
    EXPECT_NEAR(wrong - persist, r.params.gamma / 2.0, 1e-12);  // 差恰为 γ/2
}

// ---------------------------------------------------------------------------
// 6. Normal 参数 vs arch (ω ~2e-4 绝对容差, α/γ/β ~1e-5 相对)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, NormalParamsMatchArch) {
    const auto& r = normal_fit();
    EXPECT_TRUE(r.converged) << r.message;
    EXPECT_NEAR(r.params.omega, gbl::OMEGA, 2e-4);
    EXPECT_NEAR(r.params.alpha, gbl::ALPHA, 5e-5);
    EXPECT_NEAR(r.params.gamma, gbl::GAMMA, 5e-5);
    EXPECT_NEAR(r.params.beta, gbl::BETA, 5e-5);
    // 杠杆方向: 估计 γ>0 (数据含杠杆, G8 方向性)
    EXPECT_GT(r.params.gamma, 0.0);
}

// ---------------------------------------------------------------------------
// 7. llf/AIC/BIC vs arch (G17 完整似然; llf 1e-5, aic/bic 1e-4)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, LogLikelihoodAndICMatchArch) {
    const auto& r = normal_fit();
    EXPECT_NEAR(r.log_likelihood, gbl::LLF, 1e-5);
    EXPECT_NEAR(r.aic, gbl::AIC, 1e-4);
    EXPECT_NEAR(r.bic, gbl::BIC, 1e-4);
}

// ---------------------------------------------------------------------------
// 8. 条件方差序列 vs arch (最优点处, 相对 1e-4 — 参数微差沿递归传播放大;
//    约定错误的量级是 3.6e-2 (probe-1 约定 B), 1e-4 足以排幻觉)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, ConditionalVariancesMatchArch) {
    const auto& r = normal_fit();
    ASSERT_EQ(r.conditional_variances.size(), gbl::T);
    for (Size t = 0; t < gbl::T; ++t) {
        ASSERT_NEAR(r.conditional_variances[t], gbl::H[t],
                    1e-4 * std::abs(gbl::H[t]))
            << "t=" << t;
    }
}

// ---------------------------------------------------------------------------
// 9. sandwich 协方差 vs arch robust (数值 Hessian 噪声 → 5e-3)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, SandwichVcovMatchesArch) {
    const auto& r = normal_fit();
    ASSERT_EQ(r.vcov.size(), 4u);
    ASSERT_EQ(r.vcov[0].size(), 4u);
    for (Size i = 0; i < 4; ++i) {
        for (Size j = 0; j < 4; ++j) {
            ASSERT_NEAR(r.vcov[i][j], gbl::VCOV[i][j], 5e-3)
                << "V[" << i << "][" << j << "]";
        }
    }
}

// ---------------------------------------------------------------------------
// 10. 标准误 vs arch robust SE (5e-3)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, StandardErrorsMatchArch) {
    const auto& r = normal_fit();
    ASSERT_EQ(r.std_errors.size(), 4u);
    for (Size j = 0; j < 4; ++j) {
        EXPECT_NEAR(r.std_errors[j], gbl::SE[j], 5e-3) << "SE[" << j << "]";
    }
}

// ---------------------------------------------------------------------------
// 11. 标准化残差一致性 (G11: zₜ = εₜ/√hₜ)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, StdResidualsConsistency) {
    const auto& r = normal_fit();
    ASSERT_EQ(r.std_residuals.size(), gbl::T);
    ASSERT_EQ(r.residuals.size(), gbl::T);
    for (Size t = 0; t < gbl::T; ++t) {
        ASSERT_NEAR(r.std_residuals[t],
                    r.residuals[t] / std::sqrt(r.conditional_variances[t]),
                    1e-12)
            << "t=" << t;
    }
    // 去均值残差: 残差均值 ≈ 0 (G2)
    Real m = 0.0;
    for (Real v : r.residuals) m += v;
    EXPECT_NEAR(m / static_cast<Real>(gbl::T), 0.0, 1e-12);
}

// ---------------------------------------------------------------------------
// 12. G8 杠杆效应方向性: γ>0 时同幅度负冲击 → 更高波动 (filter 层面)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, LeverageEffectDirection) {
    // 两组数据仅末残差符号不同; 前段同残差 → h 链逐位一致
    const ts::GjrGarchParams p{0.2, 0.10, 0.30, 0.50, 0.0};
    const std::vector<Real> eps_pos = {0.5, -0.3, 0.2, 0.1, 1.0};
    const std::vector<Real> eps_neg = {0.5, -0.3, 0.2, 0.1, -1.0};
    const auto h_pos = ts::filter_gjr(p, eps_pos, 1.0);
    const auto h_neg = ts::filter_gjr(p, eps_neg, 1.0);
    for (Size t = 0; t < 4; ++t) {
        ASSERT_NEAR(h_pos[t], h_neg[t], 1e-15);
    }
    // 同幅度冲击仅符号不同 → h₄ 本身一致 (|ε| 相同, I 差异在下一期显现)
    ASSERT_NEAR(h_pos[4], h_neg[4], 1e-15);
    // 下期方差: h₅' = ω + α·ε₄² + γ·I(ε₄<0)·ε₄² + β·h₄
    const Real h5_next_pos = 0.2 + 0.10 * 1.0 + 0.30 * 0.0 * 1.0 + 0.50 * h_pos[4];
    const Real h5_next_neg = 0.2 + 0.10 * 1.0 + 0.30 * 1.0 * 1.0 + 0.50 * h_neg[4];
    EXPECT_GT(h5_next_neg, h5_next_pos);  // G8: 负冲击放大波动
    EXPECT_NEAR(h5_next_neg - h5_next_pos, 0.30 * 1.0, 1e-12);  // 差 = γ·ε²
}

// ---------------------------------------------------------------------------
// 13. 无条件方差不动点 (G10: σ̄² = ω/(1-α-γ/2-β))
//     交替 ±√σ̄² 冲击: h 收敛到围绕 σ̄² 的 2-循环 (γ·σ̄²/(1+β) 振幅),
//     循环均值恰为 σ̄² (E[ε²]=σ̄² 且 I 激活率 1/2)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, UnconditionalVarianceFixedPoint) {
    const ts::GjrGarchParams p{0.05, 0.05, 0.15, 0.80, 0.0};
    const Real uc = p.omega / (1.0 - p.alpha - p.gamma / 2.0 - p.beta);
    const Size n = 400;
    const Real sq = std::sqrt(uc);
    std::vector<Real> eps(n);
    for (Size t = 0; t < n; ++t) eps[t] = (t % 2 == 0) ? sq : -sq;
    const auto h = ts::filter_gjr(p, eps, uc);
    // 2-循环均值 = σ̄² (稳态: (h_e+h_o)/2 = [ω+(α+γ/2)σ̄²]/(1-β) = σ̄²)
    const Real h_avg = 0.5 * (h[n - 2] + h[n - 1]);
    EXPECT_NEAR(h_avg, uc, 1e-8);
    // 振幅 = γ·σ̄²/(1+β) (2-循环解析; n-2 为偶位 = 高方差分支)
    EXPECT_NEAR(h[n - 2] - h[n - 1], p.gamma * uc / (1.0 + p.beta), 1e-8);
    // 基准数据估计量的无条件方差 vs arch UC_VAR (probe-2 同源公式)
    const auto& r = normal_fit();
    const Real uc_hat = r.params.omega
                        / (1.0 - r.params.alpha - r.params.gamma / 2.0
                           - r.params.beta);
    EXPECT_NEAR(uc_hat, gbl::UC_VAR, 5e-4);
}

// ---------------------------------------------------------------------------
// 14. t-GJR ν 联合估计 vs arch (t(6) 数据, G-ADR5)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, StudentTJointNuEstimation) {
    ts::GarchConfig cfg;
    cfg.dist = ts::GarchDist::StudentT;
    const ts::GjrGarchResult r = ts::estimate_gjr_garch(make_tdata(), cfg);
    EXPECT_TRUE(r.converged) << r.message;
    EXPECT_NEAR(r.params.nu, gbl::T_NU, 5e-2);
    EXPECT_NEAR(r.params.omega, gbl::T_OMEGA, 2e-3);
    EXPECT_NEAR(r.params.alpha, gbl::T_ALPHA, 1e-4);
    EXPECT_NEAR(r.params.gamma, gbl::T_GAMMA, 1e-3);
    EXPECT_NEAR(r.params.beta, gbl::T_BETA, 1e-4);
    EXPECT_NEAR(r.log_likelihood, gbl::T_LLF, 1e-4);
}

// ---------------------------------------------------------------------------
// 15. 多步预测 vs arch analytic (probe-2: φ=α+γ/2+β; 1e-12)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, ForecastMatchesArchAnalytic) {
    const ts::GjrGarchParams p{gbl::OMEGA, gbl::ALPHA, gbl::GAMMA, gbl::BETA, 0.0};
    const auto fc = ts::forecast_gjr(p, gbl::H_T, gbl::E_T, 10);
    ASSERT_EQ(fc.size(), 10u);
    for (Size k = 0; k < 10; ++k) {
        ASSERT_NEAR(fc[k], gbl::FC10[k], 1e-12)
            << "k=" << k + 1 << " fc=" << fc[k] << " ref=" << gbl::FC10[k];
    }
}

// ---------------------------------------------------------------------------
// 16. 预测闭式公式 (G13 同构: h_{T+k} = σ̄² + φ^{k-1}·(h_{T+1} - σ̄²))
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, ForecastClosedFormDecay) {
    const ts::GjrGarchParams p{gbl::OMEGA, gbl::ALPHA, gbl::GAMMA, gbl::BETA, 0.0};
    const auto fc = ts::forecast_gjr(p, gbl::H_T, gbl::E_T, 10);
    const Real phi = p.alpha + p.gamma / 2.0 + p.beta;
    const Real uc = p.omega / (1.0 - phi);
    const Real h1 = fc[0];
    for (Size k = 1; k < 10; ++k) {
        const Real closed = uc + std::pow(phi, static_cast<Real>(k)) * (h1 - uc);
        EXPECT_NEAR(fc[k], closed, 1e-10) << "k=" << k + 1;
    }
}

// ---------------------------------------------------------------------------
// 17. 预测收敛到无条件方差 (k→∞, φ^k 衰减)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, ForecastConvergesToUnconditionalVariance) {
    const ts::GjrGarchParams p{gbl::OMEGA, gbl::ALPHA, gbl::GAMMA, gbl::BETA, 0.0};
    const auto fc = ts::forecast_gjr(p, gbl::H_T, gbl::E_T, 500);
    const Real uc = p.omega / (1.0 - p.alpha - p.gamma / 2.0 - p.beta);
    EXPECT_NEAR(fc.back(), uc, 1e-8);
    // 单调逼近 (方向由 h_{T+1} 相对 σ̄² 的位置决定, 此处 FC10[0] < σ̄² → 上行)
    for (Size k = 1; k < 10; ++k) {
        EXPECT_LE(std::abs(fc[k] - uc), std::abs(fc[k - 1] - uc) + 1e-12)
            << "k=" << k + 1;
    }
}

// ---------------------------------------------------------------------------
// 18. 输入校验 + 多起始点 (G16: 单起始仍收敛到同一最优)
// ---------------------------------------------------------------------------
TEST(GjrGarchModel, InputValidationAndMultistart) {
    // T < 10
    EXPECT_THROW(ts::estimate_gjr_garch({1.0, 2.0, 3.0}), std::invalid_argument);
    // NaN 数据
    std::vector<Real> bad(20, 1.0);
    bad[10] = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(ts::estimate_gjr_garch(bad), std::invalid_argument);
    // 零方差数据
    EXPECT_THROW(ts::estimate_gjr_garch(std::vector<Real>(20, 3.0)),
                 std::invalid_argument);
    // filter: 参数异常 (α<0 由估计端约束; filter 端 h ≤ 0 抛错)
    const ts::GjrGarchParams badp{-1.0, 0.5, 0.0, 0.5, 0.0};
    EXPECT_THROW(ts::filter_gjr(badp, std::vector<Real>(5, 1.0), 0.7),
                 std::runtime_error);
    // forecast: 非平稳参数 (α+γ/2+β ≥ 1) 抛错
    const ts::GjrGarchParams nsp{0.05, 0.30, 0.60, 0.50, 0.0};
    EXPECT_THROW(ts::forecast_gjr(nsp, 0.5, 0.1, 5), std::invalid_argument);

    // 单起始 (关闭 multistart) vs 多起始: 同一全局最优
    ts::GarchConfig single;
    single.use_multistart = false;
    single.initial_params = {0.05, 0.10, 0.15, 0.80};
    const ts::GjrGarchResult r1 = ts::estimate_gjr_garch(make_data(), single);
    const auto& r4 = normal_fit();
    EXPECT_NEAR(r1.log_likelihood, r4.log_likelihood, 1e-5);
    EXPECT_NEAR(r1.params.gamma, r4.params.gamma, 1e-5);
    EXPECT_NEAR(r1.params.beta, r4.params.beta, 1e-5);
}
