// =============================================================================
// test_garch_diagnostics.cpp - GARCH 标准化残差诊断测试 (15 用例, spec §2.0.6)
//
// 基准策略: 统计量层与 Phase 7A jarque_bera_test/ljung_box_test 逐位一致
// (复用即正确性), 行为层用确定性模拟数据 (Philox 固定 seed) 检验检验力方向。
//
// 幻觉点防护:
//   G11: 标准化责任在调用方 — 原始 ε (含 ARCH) 必须被诊断出 ARCH 效应,
//        标准化 z = ε/√h 后必须通过 (用例 8 直接验证)
//   G12: 自适应 lag = floor(log(T)), 非固定 10
//   G-ADR4: Bootstrap p 值确定性可复现, 退回模式与渐近逐位一致
// =============================================================================
#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>

#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/econometrics/inference/residual_diagnostics.hpp"
#include "cpphub/timeseries/garch/garch_diagnostics.hpp"

namespace ts = cpphub::v1::timeseries::garch;
namespace ec = cpphub::v1::econometrics;
using cpphub::Real;
using cpphub::Size;

// ---------------------------------------------------------------------------
// 辅助: 确定性数据生成 (Philox 固定 seed)
// ---------------------------------------------------------------------------
static std::vector<Real> make_normal(Size n, uint64_t seed) {
    cpphub::Philox4x64 rng(seed, 0xD1A9ULL);
    const Real inv_2p53 = 1.0 / 9007199254740992.0;
    std::vector<Real> z(n);
    for (Size i = 0; i < n; i += 2) {
        const Real u1 = (static_cast<Real>(rng() >> 11) + 0.5) * inv_2p53;
        const Real u2 = (static_cast<Real>(rng() >> 11) + 0.5) * inv_2p53;
        const auto [z1, z2] = cpphub::box_muller(u1, u2);
        z[i] = z1;
        if (i + 1 < n) z[i + 1] = z2;
    }
    return z;
}

// AR(1): z_t = φ z_{t-1} + w_t (强自相关)
static std::vector<Real> make_ar1(Size n, Real phi, uint64_t seed) {
    const auto w = make_normal(n, seed);
    std::vector<Real> z(n);
    z[0] = w[0] / std::sqrt(1.0 - phi * phi);
    for (Size t = 1; t < n; ++t) z[t] = phi * z[t - 1] + w[t];
    return z;
}

// 异方差 z: z_t = w_t·√h_t, h_t = 0.5 + 0.9·z²_{t-1} → z² 强自相关 (伪 ARCH)
static std::vector<Real> make_arch_z(Size n, uint64_t seed) {
    const auto w = make_normal(n, seed);
    std::vector<Real> z(n);
    z[0] = w[0];
    for (Size t = 1; t < n; ++t) {
        const Real h = 0.5 + 0.9 * z[t - 1] * z[t - 1];
        z[t] = w[t] * std::sqrt(h);
    }
    return z;
}

// Student-t(3): t = Z₁·√3/√(Z₂²+Z₃²+Z₄²) (厚尾, JB 应拒绝)
static std::vector<Real> make_t3(Size n, uint64_t seed) {
    cpphub::Philox4x64 rng(seed, 0x7C0ULL);
    const Real inv_2p53 = 1.0 / 9007199254740992.0;
    std::vector<Real> z(n);
    for (Size i = 0; i < n; ++i) {
        Real nz[4];
        for (int k = 0; k < 4; k += 2) {
            const Real u1 = (static_cast<Real>(rng() >> 11) + 0.5) * inv_2p53;
            const Real u2 = (static_cast<Real>(rng() >> 11) + 0.5) * inv_2p53;
            const auto [z1, z2] = cpphub::box_muller(u1, u2);
            nz[k] = z1;
            nz[k + 1] = z2;
        }
        z[i] = nz[0] * std::sqrt(3.0) /
               std::sqrt(nz[1] * nz[1] + nz[2] * nz[2] + nz[3] * nz[3]);
    }
    return z;
}

// ---------------------------------------------------------------------------
// 1. iid 正态 z: 三项全通过 + summary 非空
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, NormalZAllPass) {
    const auto z = make_normal(500, 101);
    const auto r = ts::diagnose_garch_residuals(z, 500);
    EXPECT_TRUE(r.passes_normality);
    EXPECT_TRUE(r.passes_no_autocorr);
    EXPECT_TRUE(r.passes_no_arch_effect);
    EXPECT_FALSE(r.summary.empty());
    EXPECT_NE(r.summary.find("LB(z^2)"), std::string::npos);
}

// ---------------------------------------------------------------------------
// 2. G12 自适应 lag = floor(log(T)), z 与 z² 同 lag
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, AutoLagFloorLogT) {
    EXPECT_EQ(ts::diagnostics_auto_lag(1000), 6);  // log(1000)=6.907
    EXPECT_EQ(ts::diagnostics_auto_lag(500), 6);   // log(500)=6.215
    EXPECT_EQ(ts::diagnostics_auto_lag(100), 4);   // log(100)=4.605
    EXPECT_EQ(ts::diagnostics_auto_lag(10), 2);    // log(10)=2.303
    EXPECT_EQ(ts::diagnostics_auto_lag(3), 1);     // 下限保护

    const auto z = make_normal(500, 102);
    const auto r = ts::diagnose_garch_residuals(z, 0);
    EXPECT_EQ(r.lb_z.lag, 6);
    EXPECT_EQ(r.lb_z_squared.lag, 6);  // G12: 同 lag 保持子检验可比
}

// ---------------------------------------------------------------------------
// 3. 显式 lb_lag 覆盖
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, ExplicitLagOverride) {
    const auto z = make_normal(300, 103);
    const auto r = ts::diagnose_garch_residuals(z, 0, 7);
    EXPECT_EQ(r.lb_z.lag, 7);
    EXPECT_EQ(r.lb_z_squared.lag, 7);
}

// ---------------------------------------------------------------------------
// 4. JB 统计量与 Phase 7A 逐位一致 (复用即正确性)
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, JBMatchesPhase7A) {
    const auto z = make_normal(400, 104);
    const auto r = ts::diagnose_garch_residuals(z, 0);
    const auto jb_ref = ec::jarque_bera_test(z);
    EXPECT_DOUBLE_EQ(r.jb_test.base.statistic, jb_ref.base.statistic);
    EXPECT_DOUBLE_EQ(r.jb_test.skewness, jb_ref.skewness);
    EXPECT_DOUBLE_EQ(r.jb_test.kurtosis, jb_ref.kurtosis);
    EXPECT_DOUBLE_EQ(r.jb_test.base.p_value, jb_ref.base.p_value);
    EXPECT_EQ(r.jb_test.base.method_name, "Jarque-Bera");  // 无 bootstrap 后缀
}

// ---------------------------------------------------------------------------
// 5. LB 统计量与 Phase 7A 逐位一致
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, LBMatchesPhase7A) {
    const auto z = make_normal(400, 105);
    const auto r = ts::diagnose_garch_residuals(z, 0);
    const auto lb_ref = ec::ljung_box_test(z, r.lb_z.lag);
    EXPECT_DOUBLE_EQ(r.lb_z.base.statistic, lb_ref.base.statistic);
    EXPECT_DOUBLE_EQ(r.lb_z_squared.base.statistic,
                     ec::ljung_box_test([&] {
                         std::vector<Real> z2(z.size());
                         for (Size i = 0; i < z.size(); ++i) z2[i] = z[i] * z[i];
                         return z2;
                     }(),
                                        r.lb_z_squared.lag)
                         .base.statistic);
}

// ---------------------------------------------------------------------------
// 6. 自相关检验力: AR(1) φ=0.5 → LB(z) 拒绝
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, AutocorrelationDetection) {
    const auto z = make_ar1(500, 0.5, 106);
    const auto r = ts::diagnose_garch_residuals(z, 0);
    EXPECT_FALSE(r.passes_no_autocorr);
    EXPECT_TRUE(r.lb_z.base.reject_null);
    EXPECT_LT(r.lb_z.base.p_value, 0.05);
}

// ---------------------------------------------------------------------------
// 7. ARCH 效应检验力: 异方差 z → LB(z²) 拒绝 (Li-Mak)
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, ArchEffectDetection) {
    const auto z = make_arch_z(500, 107);
    const auto r = ts::diagnose_garch_residuals(z, 0);
    EXPECT_FALSE(r.passes_no_arch_effect);
    EXPECT_TRUE(r.lb_z_squared.base.reject_null);
    EXPECT_LT(r.lb_z_squared.base.p_value, 0.05);
}

// ---------------------------------------------------------------------------
// 8. G11 标准化责任: 原始 ε (ARCH) 必须检出; z = ε/√h 后必须通过
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, G11StandardizationResponsibility) {
    // 构造 GARCH(1,1) 数据 ε_t = w_t·√h_t, h_t = 0.1 + 0.15ε² + 0.80h_{t-1}
    // (α=0.15: ε² 的 ρ₁≈0.30, LB(6) 检验力充足; α=0.05 时 ρ₁≈0.07 过弱)
    const auto w = make_normal(600, 108);
    const Size n = w.size();
    std::vector<Real> eps(n), h(n);
    h[0] = 0.1 / (1.0 - 0.15 - 0.80);  // 无条件方差
    eps[0] = w[0] * std::sqrt(h[0]);
    for (Size t = 1; t < n; ++t) {
        h[t] = 0.1 + 0.15 * eps[t - 1] * eps[t - 1] + 0.80 * h[t - 1];
        eps[t] = w[t] * std::sqrt(h[t]);
    }

    // 原始 ε: 未标准化 → ARCH 效应必须被检出
    const auto r_raw = ts::diagnose_garch_residuals(eps, 0);
    EXPECT_FALSE(r_raw.passes_no_arch_effect);

    // 标准化 z = ε/√h: 条件异方差已消除 → ARCH 效应必须通过
    std::vector<Real> zstd(n);
    for (Size t = 0; t < n; ++t) zstd[t] = eps[t] / std::sqrt(h[t]);
    const auto r_std = ts::diagnose_garch_residuals(zstd, 0);
    EXPECT_TRUE(r_std.passes_no_arch_effect);
    EXPECT_TRUE(r_std.passes_no_autocorr);  // w iid → z iid
}

// ---------------------------------------------------------------------------
// 9. 非正态检验力: t(3) 厚尾 → JB 拒绝 (渐近 + Bootstrap 双模式)
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, NonNormalityDetection) {
    const auto z = make_t3(500, 109);
    // 渐近
    const auto r_asy = ts::diagnose_garch_residuals(z, 0);
    EXPECT_FALSE(r_asy.passes_normality);
    EXPECT_GT(r_asy.jb_test.kurtosis, 6.0);  // t(3) 峰度 = 9 (非超额 9? 超额 6)
    // Bootstrap
    const auto r_boot = ts::diagnose_garch_residuals(z, 500);
    EXPECT_FALSE(r_boot.passes_normality);
    EXPECT_LT(r_boot.jb_test.base.p_value, 0.05);
}

// ---------------------------------------------------------------------------
// 10. G-ADR4: 正态 z 的 Bootstrap p ≈ 渐近 p (MC 误差内)
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, BootstrapPvalueCloseToAsymptotic) {
    const auto z = make_normal(500, 110);
    const Real p_asy = ec::jarque_bera_test(z).base.p_value;
    const auto r = ts::diagnose_garch_residuals(z, 1000);
    EXPECT_NEAR(r.jb_test.base.p_value, p_asy, 0.05);
    EXPECT_TRUE(r.passes_normality);
}

// ---------------------------------------------------------------------------
// 11. Bootstrap 确定性: 同输入两次 → p 逐位一致
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, BootstrapDeterministic) {
    const auto z = make_normal(300, 111);
    const Real p1 = ts::jb_bootstrap_pvalue(z, 400);
    const Real p2 = ts::jb_bootstrap_pvalue(z, 400);
    EXPECT_DOUBLE_EQ(p1, p2);
    const auto r1 = ts::diagnose_garch_residuals(z, 400);
    const auto r2 = ts::diagnose_garch_residuals(z, 400);
    EXPECT_DOUBLE_EQ(r1.jb_test.base.p_value, r2.jb_test.base.p_value);
}

// ---------------------------------------------------------------------------
// 12. n_bootstrap = 0: 退回渐近 (p 逐位一致, method 无后缀)
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, BootstrapZeroRepsFallsBackToAsymptotic) {
    const auto z = make_normal(300, 112);
    const auto r = ts::diagnose_garch_residuals(z, 0);
    EXPECT_EQ(r.jb_bootstrap_reps, 0u);
    EXPECT_DOUBLE_EQ(r.jb_test.base.p_value, ec::jarque_bera_test(z).base.p_value);
    EXPECT_EQ(r.jb_test.base.method_name, "Jarque-Bera");
}

// ---------------------------------------------------------------------------
// 13. Bootstrap 元数据: reps 记录 + method 后缀
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, BootstrapRepsRecorded) {
    const auto z = make_normal(200, 113);
    const auto r = ts::diagnose_garch_residuals(z, 250);
    EXPECT_EQ(r.jb_bootstrap_reps, 250u);
    EXPECT_EQ(r.jb_test.base.method_name, "Jarque-Bera (bootstrap)");
    // 统计量仍是观测样本的 JB (bootstrap 只替换 p 值)
    EXPECT_DOUBLE_EQ(r.jb_test.base.statistic, ec::jarque_bera_test(z).base.statistic);
    // p ∈ [1/(B+1), 1] (Davison-Hinkley 加一修正下界)
    EXPECT_GE(r.jb_test.base.p_value, 1.0 / 251.0);
    EXPECT_LE(r.jb_test.base.p_value, 1.0);
}

// ---------------------------------------------------------------------------
// 14. 异常输入: 长度不足 / NaN / 非法 lag
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, InvalidInput) {
    const std::vector<Real> short_z(9, 0.1);
    EXPECT_THROW(ts::diagnose_garch_residuals(short_z), std::invalid_argument);

    auto z = make_normal(200, 114);
    z[100] = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(ts::diagnose_garch_residuals(z), std::invalid_argument);

    const auto z2 = make_normal(200, 115);
    EXPECT_THROW(ts::diagnose_garch_residuals(z2, 0, 200), std::invalid_argument);
    // 注: lb_lag=0 合法 (spec: 0 = 自动 floor(log(T))), 不应抛出
    EXPECT_NO_THROW(ts::diagnose_garch_residuals(z2, 0, 0));
    EXPECT_THROW(ts::jb_bootstrap_pvalue(z2, 0), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 15. 退化输入: 零方差 → Phase 7A jarque_bera_test 抛 runtime_error
// ---------------------------------------------------------------------------
TEST(GarchDiagnostics, DegenerateZeroVarianceThrows) {
    const std::vector<Real> flat(100, 0.5);
    EXPECT_THROW(ts::diagnose_garch_residuals(flat), std::runtime_error);
}
