// =============================================================================
// test_egarch_model.cpp - EGARCH(1,1) QMLE 测试 (18 用例, spec §2.0.3 测试矩阵)
//
// 基准: tests/unit/timeseries/egarch_baseline.inc (arch 8.0.0 自动生成,
//       脚本 tests/fixtures/timeseries/verify_egarch.py, HFE 硬编码策略)
//
// 容差策略 (跨优化器对照, 非 HFE 同算法对照):
//   - 确定性函数 (filter/forecast/backcast): 1e-12 (与 arch 逐位一致)
//   - 优化结果 (params/llf): ~1e-4/1e-6 (scipy SLSQP 解析梯度 vs 自研数值
//     梯度, 收敛到同一最优点的精度差)
//   - sandwich vcov: ~3e-2/5e-2 (数值 Hessian 歮差方案噪声)
//   - simulation 多步预测中位数: ~10% (Monte Carlo 噪声 + ln h 递归中位数
//     vs 均值的偏度修正, γ 小时可忽略)
//
// 参数映射 (G23): spec {omega, alpha(非对称 z 项), beta, gamma(对称项)};
//   arch 顺序 {omega, alpha[1](对称), gamma[1](非对称), beta[1]} — 基准
//   已在 Python 端映射为 spec 顺序
// =============================================================================
#include <gtest/gtest.h>
#include <chrono>
#include <vector>
#include <cmath>

#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/timeseries/garch/egarch_model.hpp"
#include "cpphub/timeseries/garch/garch_forecast.hpp"
#include "egarch_baseline.inc"

namespace ts = cpphub::v1::timeseries::garch;
namespace ebl = cpphub::v1::timeseries::garch::egarch_baseline;
using cpphub::Real;
using cpphub::Size;

// ---------------------------------------------------------------------------
// 辅助: 数据构造 + 结果缓存 (estimate ~1s, 静态缓存避免重复估计)
// ---------------------------------------------------------------------------
static std::vector<Real> make_data() {
    return {ebl::DATA, ebl::DATA + ebl::T};
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
    return {ebl::TDATA, ebl::TDATA + ebl::T};
}

static const ts::EGarchResult& normal_fit() {
    static ts::EGarchResult r = ts::estimate_egarch(make_data());
    return r;
}

// ---------------------------------------------------------------------------
// 1. E|z| = √(2/π) 常数 (G5, 1e-15; arch 源码 np.sqrt(2/np.pi))
// ---------------------------------------------------------------------------
TEST(EgarchModel, EAbsZConstantMatchesSqrt2OverPi) {
    const Real ref = std::sqrt(2.0 / 3.14159265358979323846);
    EXPECT_NEAR(ts::EGARCH_E_ABS_Z, ref, 1e-15);
    // 排幻觉: 非 2/π ≈ 0.6366
    EXPECT_GT(ts::EGARCH_E_ABS_Z, 0.79);
    EXPECT_LT(ts::EGARCH_E_ABS_Z, 0.80);
}

// ---------------------------------------------------------------------------
// 2. log backcast vs arch EGARCH.backcast (G-q: log 尺度, 负值)
// ---------------------------------------------------------------------------
TEST(EgarchModel, LogBackcastMatchesArch) {
    const auto dm = make_data_demeaned();
    const Real bc_raw = ts::backcast_variance(dm, 0.94, 0);
    // arch EGARCH.backcast() = log(super().backcast()) → log 尺度 (可为负)
    EXPECT_NEAR(std::log(bc_raw), ebl::BACKCAST, 1e-12);
    // 印证: 该数据下 backcast 为负 (log 尺度), raw ≈ exp(负)
    EXPECT_LT(ebl::BACKCAST, 0.0);
    EXPECT_NEAR(bc_raw, std::exp(ebl::BACKCAST), 1e-12);
}

// ---------------------------------------------------------------------------
// 3. 方差递归 filter_egarch vs arch conditional_volatility² (1e-12)
//    用 arch 参数直接递归 — 排除优化器差异, 纯算法对照 (含 G-q t=0 初始化)
// ---------------------------------------------------------------------------
TEST(EgarchModel, FilterRecursionMatchesArchConditionalVariance) {
    const auto dm = make_data_demeaned();
    const ts::EGarchParams p{ebl::OMEGA, ebl::ALPHA, ebl::BETA, ebl::GAMMA, 0.0};
    const auto h = ts::filter_egarch(p, dm, ebl::BACKCAST);
    ASSERT_EQ(h.size(), ebl::T);
    for (Size t = 0; t < ebl::T; ++t) {
        ASSERT_NEAR(h[t], ebl::H[t], 1e-12)
            << "t=" << t << " h=" << h[t] << " ref=" << ebl::H[t];
    }
}

// ---------------------------------------------------------------------------
// 4. t=0 初始化 + 手算递归 (G-q + G6: 非对称项用标准化 zₜ₋₁)
// ---------------------------------------------------------------------------
TEST(EgarchModel, FilterT0InitAndHandComputedRecursion) {
    // params: ω=0.1, α=-0.3, β=0.8, γ=0.2; bc_log=-1.5
    const ts::EGarchParams p{0.1, -0.3, 0.8, 0.2, 0.0};
    const std::vector<Real> eps = {1.0, 0.5, -0.2};
    const Real bc_log = -1.5;
    const auto h = ts::filter_egarch(p, eps, bc_log);

    // t=0: ln h₀ = ω + β·bc_log (G-q)
    const Real ln_h0 = 0.1 + 0.8 * (-1.5);
    ASSERT_NEAR(h[0], std::exp(ln_h0), 1e-12);
    // t=1: z₀ = ε₀/√h₀ (G6 标准化; h₀=exp(-1.1)≈0.333 → z₀≈3 ≠ ε₀=1)
    const Real z0 = 1.0 / std::sqrt(std::exp(ln_h0));
    const Real ln_h1 = 0.1 + 0.8 * ln_h0 + (-0.3) * z0
                       + 0.2 * (std::abs(z0) - ts::EGARCH_E_ABS_Z);
    ASSERT_NEAR(h[1], std::exp(ln_h1), 1e-12);
    // t=2: 链式递归
    const Real z1 = 0.5 / std::sqrt(std::exp(ln_h1));
    const Real ln_h2 = 0.1 + 0.8 * ln_h1 + (-0.3) * z1
                       + 0.2 * (std::abs(z1) - ts::EGARCH_E_ABS_Z);
    ASSERT_NEAR(h[2], std::exp(ln_h2), 1e-12);
    // G6 排幻觉: 若误用未标准化 ε₀ (=1.0), h₁ 偏离 ~7.5% (α+γ 均作用于差)
    const Real ln_h1_wrong = 0.1 + 0.8 * ln_h0 + (-0.3) * 1.0
                             + 0.2 * (1.0 - ts::EGARCH_E_ABS_Z);
    EXPECT_GT(std::abs(h[1] - std::exp(ln_h1_wrong)),
              0.05 * std::abs(h[1]));
}

// ---------------------------------------------------------------------------
// 5. G23 参数映射: spec 顺序 filter ≡ arch 顺序手写递归 (1e-14)
// ---------------------------------------------------------------------------
TEST(EgarchModel, G23ParameterMappingConsistency) {
    // arch 顺序递归: th = {ω, alpha_sym(乘 |z|-E|z|), gamma_asym(乘 z), β}
    auto ref_arch = [](const std::vector<Real>& th,
                       const std::vector<Real>& eps, Real bc_log) {
        const Size T = eps.size();
        std::vector<Real> h(T);
        Real ln_prev = th[0] + th[3] * bc_log;
        for (Size t = 0; t < T; ++t) {
            if (t > 0) {
                const Real z = eps[t - 1] / std::sqrt(h[t - 1]);
                ln_prev = th[0] + th[3] * ln_prev + th[2] * z
                          + th[1] * (std::abs(z) - ts::EGARCH_E_ABS_Z);
            }
            h[t] = std::exp(ln_prev);
        }
        return h;
    };
    const ts::EGarchParams p{-0.02, -0.08, 0.91, 0.17, 0.0};
    const std::vector<Real> eps = {0.5, -1.2, 0.8, -0.3, 1.5, -0.9, 0.2};
    const auto h_spec = ts::filter_egarch(p, eps, -0.4);
    // 映射: arch alpha_sym = spec gamma, arch gamma_asym = spec alpha
    const auto h_arch =
        ref_arch({p.omega, p.gamma, p.alpha, p.beta}, eps, -0.4);
    ASSERT_EQ(h_spec.size(), h_arch.size());
    for (Size t = 0; t < h_spec.size(); ++t) {
        ASSERT_NEAR(h_spec[t], h_arch[t], 1e-14) << "t=" << t;
    }
}

// ---------------------------------------------------------------------------
// 6. 参数 vs arch + |β|<1 约束 + 确定性
// ---------------------------------------------------------------------------
TEST(EgarchModel, ParamsMatchArchAndConstraintsHold) {
    const auto& r = normal_fit();
    EXPECT_TRUE(r.converged) << r.message;
    EXPECT_NEAR(r.params.omega, ebl::OMEGA, 1e-3);
    EXPECT_NEAR(r.params.alpha, ebl::ALPHA, 1e-4);
    EXPECT_NEAR(r.params.beta, ebl::BETA, 2e-4);
    EXPECT_NEAR(r.params.gamma, ebl::GAMMA, 1e-4);
    // 平稳性 |β| < 1 (spec §2.0.3 Step 3.2, bounds β ∈ [0,1))
    EXPECT_GE(r.params.beta, 0.0);
    EXPECT_LT(r.params.beta, 1.0);
    // 确定性: 重复运行逐位一致 (Philox 种子固定)
    const auto r2 = ts::estimate_egarch(make_data());
    EXPECT_EQ(r2.params.omega, r.params.omega);
    EXPECT_EQ(r2.params.alpha, r.params.alpha);
    EXPECT_EQ(r2.params.beta, r.params.beta);
    EXPECT_EQ(r2.params.gamma, r.params.gamma);
    EXPECT_EQ(r2.log_likelihood, r.log_likelihood);
}

// ---------------------------------------------------------------------------
// 7. 对数似然 vs arch (G3: 含常数项)
// ---------------------------------------------------------------------------
TEST(EgarchModel, LogLikelihoodMatchesArch) {
    const auto& r = normal_fit();
    EXPECT_NEAR(r.log_likelihood, ebl::LLF, 1e-6);
}

// ---------------------------------------------------------------------------
// 8. AIC/BIC vs arch (G17: 完整似然, k=4)
// ---------------------------------------------------------------------------
TEST(EgarchModel, AicBicMatchArch) {
    const auto& r = normal_fit();
    EXPECT_NEAR(r.aic, ebl::AIC, 1e-5);
    EXPECT_NEAR(r.bic, ebl::BIC, 1e-5);
    // 公式恒等式: AIC = -2ℓ + 2k (k=4), BIC = -2ℓ + k·log(T)
    EXPECT_NEAR(r.aic, -2.0 * r.log_likelihood + 8.0, 1e-9);
    EXPECT_NEAR(r.bic,
                -2.0 * r.log_likelihood
                    + 4.0 * std::log(static_cast<Real>(ebl::T)),
                1e-9);
}

// ---------------------------------------------------------------------------
// 9. 条件方差: 与公开 filter 自洽 (无尺度变换 → 逐位) + vs arch (1e-3)
// ---------------------------------------------------------------------------
TEST(EgarchModel, ConditionalVariancesSelfConsistentAndMatchArch) {
    const auto& r = normal_fit();
    ASSERT_EQ(r.conditional_variances.size(), ebl::T);
    // 自洽: 估计器内部与公开 filter 同一路径 (无尺度变换, 逐位一致)
    const auto dm = make_data_demeaned();
    const Real bc_log = std::log(ts::backcast_variance(dm, 0.94, 0));
    const auto h_ref = ts::filter_egarch(r.params, dm, bc_log);
    for (Size t = 0; t < ebl::T; ++t) {
        ASSERT_NEAR(r.conditional_variances[t], h_ref[t], 1e-14) << "t=" << t;
    }
    // vs arch (参数差 ~1e-4 传播到 h)
    for (Size t = 0; t < ebl::T; ++t) {
        ASSERT_NEAR(r.conditional_variances[t], ebl::H[t], 1e-3)
            << "t=" << t;
    }
}

// ---------------------------------------------------------------------------
// 10. sandwich vcov vs arch robust cov (G9, spec 顺序 4×4)
// ---------------------------------------------------------------------------
TEST(EgarchModel, SandwichVcovMatchesArchRobust) {
    const auto& r = normal_fit();
    ASSERT_EQ(r.vcov.size(), 4u);
    for (Size i = 0; i < 4; ++i) {
        for (Size j = 0; j < 4; ++j) {
            const Real ref = ebl::VCOV[i][j];
            if (i == j) {
                ASSERT_NEAR(r.vcov[i][j], ref, 3e-2 * std::abs(ref))
                    << "vcov[" << i << "][" << j << "]";
            } else {
                const Real sc = std::sqrt(std::abs(ebl::VCOV[i][i]
                                                   * ebl::VCOV[j][j]));
                ASSERT_NEAR(r.vcov[i][j], ref, 5e-2 * sc)
                    << "vcov[" << i << "][" << j << "]";
            }
        }
    }
    for (Size i = 0; i < 4; ++i)
        for (Size j = i + 1; j < 4; ++j)
            EXPECT_NEAR(r.vcov[i][j], r.vcov[j][i], 1e-12);
}

// ---------------------------------------------------------------------------
// 11. 标准误 = sqrt(diag(V)) 恒等式 + vs arch
// ---------------------------------------------------------------------------
TEST(EgarchModel, StdErrorsIdentityAndMatchArch) {
    const auto& r = normal_fit();
    ASSERT_EQ(r.std_errors.size(), 4u);
    for (Size i = 0; i < 4; ++i) {
        EXPECT_NEAR(r.std_errors[i], std::sqrt(r.vcov[i][i]), 1e-12);
        EXPECT_NEAR(r.std_errors[i], ebl::SE[i], 2e-2 * ebl::SE[i]);
    }
}

// ---------------------------------------------------------------------------
// 12. 标准化残差恒等式 zₜ = εₜ/√hₜ (G11) + E[z²]≈1
// ---------------------------------------------------------------------------
TEST(EgarchModel, StdResidualsIdentity) {
    const auto& r = normal_fit();
    ASSERT_EQ(r.std_residuals.size(), ebl::T);
    for (Size t = 0; t < ebl::T; ++t) {
        const Real z = r.residuals[t] / std::sqrt(r.conditional_variances[t]);
        ASSERT_NEAR(r.std_residuals[t], z, 1e-10) << "t=" << t;
    }
    Real m2 = 0.0;
    for (Real z : r.std_residuals) m2 += z * z;
    m2 /= static_cast<Real>(r.std_residuals.size());
    EXPECT_NEAR(m2, 1.0, 0.1);  // QMLE 一阶条件
}

// ---------------------------------------------------------------------------
// 13. 杠杆效应方向性 (α<0: 负冲击放大波动) + 精确公式
// ---------------------------------------------------------------------------
TEST(EgarchModel, LeverageEffectDirectionAndFormula) {
    // α = -0.1 (Nelson θ < 0 → 杠杆), γ=0.2; 同幅度冲击: z=-2 vs z=+2
    const ts::EGarchParams p{0.0, -0.1, 0.9, 0.2, 0.0};
    const auto fc_neg = ts::forecast_egarch(p, 0.0, -2.0, 1);
    const auto fc_pos = ts::forecast_egarch(p, 0.0, 2.0, 1);
    ASSERT_EQ(fc_neg.size(), 1u);
    // ln h_{T+1}(z=-2) - ln h_{T+1}(z=+2) = α·(-2) - α·(+2) = -4α > 0
    EXPECT_GT(fc_neg[0], fc_pos[0]);
    EXPECT_NEAR(std::log(fc_neg[0] / fc_pos[0]), -4.0 * p.alpha, 1e-12);
    // size 效应 (γ>0): |z|=2 vs z=0 同等方向下 |z| 大 → 波动大
    const auto fc_zero = ts::forecast_egarch(p, 0.0, 0.0, 1);
    EXPECT_GT(fc_pos[0], fc_zero[0]);
    // 杠杆非对称量级: 负冲击 ln h = α·(-2) 与 size 项 γ·(2-E|z|)
    const Real ln_neg = 0.0 + 0.9 * 0.0 + (-0.1) * (-2.0)
                        + 0.2 * (2.0 - ts::EGARCH_E_ABS_Z);
    EXPECT_NEAR(std::log(fc_neg[0]), ln_neg, 1e-12);
}

// ---------------------------------------------------------------------------
// 14. 1-step 预测 vs arch analytic (纯算术, 1e-10)
// ---------------------------------------------------------------------------
TEST(EgarchModel, ForecastOneStepMatchesArchAnalytic) {
    const ts::EGarchParams p{ebl::OMEGA, ebl::ALPHA, ebl::BETA, ebl::GAMMA, 0.0};
    const auto fc = ts::forecast_egarch(p, ebl::LN_H_T, ebl::Z_T, 1);
    ASSERT_EQ(fc.size(), 1u);
    EXPECT_NEAR(fc[0], ebl::FC1, 1e-10);
}

// ---------------------------------------------------------------------------
// 15. 多步预测收敛到 exp(ω/(1-β)) (1e-6) + 单调性
// ---------------------------------------------------------------------------
TEST(EgarchModel, ForecastConvergesToExpOmegaOverOneMinusBeta) {
    const ts::EGarchParams p{-0.02, -0.05, 0.94, 0.15, 0.0};
    const Real limit = std::exp(p.omega / (1.0 - p.beta));
    const auto fc = ts::forecast_egarch(p, 1.2, -1.5, 3000);
    EXPECT_NEAR(fc.back(), limit, 1e-6 * limit);
    // ln 空间单调逼近 (|ln h_{T+k} - limit_ln| 递减)
    const Real lim_ln = p.omega / (1.0 - p.beta);
    Real prev_gap = std::numeric_limits<Real>::infinity();
    for (Size k : {0u, 10u, 100u, 1000u}) {
        const Real gap = std::abs(std::log(fc[k]) - lim_ln);
        EXPECT_LT(gap, prev_gap);
        prev_gap = gap;
    }
}

// ---------------------------------------------------------------------------
// 16. 多步预测: 矩母函数闭式 E[h_{T+k}] vs arch simulation 路径均值 (3%)
//
// 数学 (probe_egarch_sim.py 验证): ln h_{T+k} = c_k + Σ_{j=1}^{k-1} β^{k-1-j}ξⱼ,
// ξⱼ = α·zⱼ + γ·(|zⱼ|-c) iid, z ~ N(0,1); c_k = 确定性递归 (本实现)。
// E[h_{T+k}] = e^{c_k}·Πⱼ M(tⱼ), tⱼ = β^{k-1-j}, M 为 ξ 的矩母函数:
//   E[e^{tξ}] = e^{-tγc}·[e^{t²(α-γ)²/2}·Φ(-t(α-γ)) + e^{t²(α+γ)²/2}·Φ(t(α+γ))]
// arch simulation 输出 fc.variance = 路径 h 均值 (E[h]); 手算 E[h_{T+2}]
// = 0.6501 vs arch 0.6477 (1000 路径 MC 噪声 0.4%)。同时验证 Jensen 方向:
// E[h] ≥ exp(c_k) = 本实现返回值 (中位数式预测)
// ---------------------------------------------------------------------------
TEST(EgarchModel, ForecastMultiStepMomentClosedFormVsArchSimulation) {
    const ts::EGarchParams p{ebl::OMEGA, ebl::ALPHA, ebl::BETA, ebl::GAMMA, 0.0};
    const auto fc = ts::forecast_egarch(p, ebl::LN_H_T, ebl::Z_T, 10);
    ASSERT_EQ(fc.size(), 10u);

    const auto phi_cdf = [](Real x) {
        return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
    };
    // M(t) = E[e^{tξ}], ξ = αz + γ(|z|-c)
    const auto mgf = [&](Real t) {
        const Real um = t * (p.alpha - p.gamma);
        const Real up = t * (p.alpha + p.gamma);
        const Real m = std::exp(0.5 * um * um) * phi_cdf(-um)
                       + std::exp(0.5 * up * up) * phi_cdf(up);
        return std::exp(-t * p.gamma * ts::EGARCH_E_ABS_Z) * m;
    };

    for (Size k = 0; k < 10; ++k) {
        // E[h_{T+k}] 闭式: e^{c_k}·Π_{j=1}^{k-1} M(β^{k-1-j})
        Real prod = 1.0;
        for (Size j = 1; j <= k; ++j) {
            prod *= mgf(std::pow(p.beta, static_cast<Real>(k - j)));
        }
        const Real eh = fc[k] * prod;  // fc[k] = exp(c_k)
        EXPECT_NEAR(eh, ebl::FC_SIM_VAR[k], 0.03 * ebl::FC_SIM_VAR[k])
            << "k=" << k + 1 << " closed=" << eh
            << " arch=" << ebl::FC_SIM_VAR[k];
        // Jensen 方向: E[h] ≥ exp(E ln h) (方差项使矩母函数乘积 ≥ 1)
        EXPECT_GE(eh, fc[k] * (1.0 - 1e-12));
    }
    // 渐近极限 (验证用, 不参与断言数值): E[h_∞] = exp(ω/(1-β) + Var(ξ)/(2(1-β²))),
    // Var(ξ) = α² + γ²(1-2/π) > 确定性极限 exp(ω/(1-β)); k=10 时 E[h] 尚未
    // 收敛 (0.835 < 0.971, 偏差以 β^k 衰减), 故仅验证序列单调上行逼近
    for (Size k = 1; k < 10; ++k) {
        EXPECT_GT(ebl::FC_SIM_VAR[k], ebl::FC_SIM_VAR[k - 1] - 0.02);
    }
}

// ---------------------------------------------------------------------------
// 17. t-EGARCH ν 联合估计 vs arch dist='t' (G-ADR5)
// ---------------------------------------------------------------------------
TEST(EgarchModel, TEgarchNuJointEstimation) {
    ts::GarchConfig cfg;
    cfg.dist = ts::GarchDist::StudentT;
    const auto r = ts::estimate_egarch(make_tdata(), cfg);
    EXPECT_TRUE(r.converged) << r.message;
    EXPECT_NEAR(r.params.nu, ebl::T_NU, 5e-2);
    EXPECT_NEAR(r.params.omega, ebl::T_OMEGA, 2e-3);
    EXPECT_NEAR(r.params.alpha, ebl::T_ALPHA, 2e-3);
    EXPECT_NEAR(r.params.beta, ebl::T_BETA, 5e-3);
    EXPECT_NEAR(r.params.gamma, ebl::T_GAMMA, 2e-3);
    EXPECT_NEAR(r.log_likelihood, ebl::T_LLF, 1e-6);
    EXPECT_GE(r.params.nu, 2.05);
    EXPECT_LE(r.params.nu, 500.0);
}

// ---------------------------------------------------------------------------
// 18. 输入校验 + 多起始不劣于单起始 + forecast 边界
// ---------------------------------------------------------------------------
TEST(EgarchModel, InputValidationAndMultistart) {
    // T < 10
    EXPECT_THROW(ts::estimate_egarch({1.0, 2.0, 3.0}),
                 std::invalid_argument);
    // NaN 数据
    auto d = make_data();
    d[5] = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(ts::estimate_egarch(d), std::invalid_argument);
    // 零方差 (全同值)
    EXPECT_THROW(ts::estimate_egarch(std::vector<Real>(100, 3.14)),
                 std::invalid_argument);
    // forecast 边界
    const ts::EGarchParams p{-0.02, -0.05, 0.94, 0.15, 0.0};
    EXPECT_TRUE(ts::forecast_egarch(p, 0.1, 0.5, 0).empty());
    EXPECT_THROW(ts::forecast_egarch(p, std::numeric_limits<Real>::quiet_NaN(),
                                     0.5, 5),
                 std::invalid_argument);
    EXPECT_THROW(ts::forecast_egarch(ts::EGarchParams{-0.02, -0.05, 1.2, 0.15, 0.0},
                                     0.1, 0.5, 5),
                 std::invalid_argument);
    // 多起始不劣于单起始 (G16)
    ts::GarchConfig single;
    single.use_multistart = false;
    const auto r1 = ts::estimate_egarch(make_data(), single);
    const auto& r4 = normal_fit();
    EXPECT_TRUE(r1.converged);
    EXPECT_GE(r4.log_likelihood, r1.log_likelihood - 1e-8);
}
