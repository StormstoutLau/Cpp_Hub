// =============================================================================
// garch_diagnostics.hpp - GARCH 标准化残差诊断 (spec §2.0.6)
//
// Phase 7B v1.6 M1 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 教材锚点: Tsay 3ed Ch 5 / Li-Mak 1994 / Phase 7A residual_diagnostics.hpp
// 复用: jarque_bera_test() / ljung_box_test() (Phase 7A, ADR-015 方案 B,
//       不重复实现)
//
// 幻觉点防护 (spec §2.0.6):
//   G11:   输入是标准化残差 zₜ = εₜ/√hₜ (非原始 εₜ, 非 εₜ/σₜ 标准差符号)
//   G12:   z²ₜ LB 滞后自适应 floor(log(T)) (非固定 10; Tsay Ch5 诊断惯例)
//   G-ADR4: JB 检验小样本用参数化 Bootstrap p 值 (渐近 χ²(2) 在 N<200
//          过度拒绝; 抽样 iid N(0,1) — H0 下 z 的总体即标准正态)
//
// Bootstrap 实现说明: spec 名义复用 v1.5 block_bootstrap.hpp, 但该组件
// 面向估计器对象 + Eigen3 (linalg_dynamic), 而此处只需 iid 正态参数化
// bootstrap — 按 ADR-015 方案 B 精神直接用 core Philox, 避免 timeseries
// 模块引入 Eigen3 依赖。固定 seed 保证可复现 (同输入 → 同 p 值, 逐位)。
//
// p_boot = (1 + #{JB* ≥ JB_obs}) / (B + 1)  (Davison-Hinkley 加一修正,
//          避免 p = 0; B = n_bootstrap, 0 = 退回渐近 χ²)
// =============================================================================
#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/econometrics/inference/residual_diagnostics.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace garch {

// GARCH 标准化残差诊断结果 (聚合多个子检验, 复合诊断不组合 base — ADR-015 决策点 3)
struct GarchDiagnosticsResult {
    // L1 模型内诊断
    econometrics::JarqueBeraResult jb_test;       ///< 正态性 (Bootstrap p 值, G-ADR4)
    econometrics::LjungBoxResult lb_z;            ///< zₜ 自相关 (G11)
    econometrics::LjungBoxResult lb_z_squared;    ///< z²ₜ ARCH 效应 (G12/Li-Mak)

    // 诊断结论
    bool passes_normality;                        ///< JB 检验是否通过 (p > 0.05)
    bool passes_no_autocorr;                      ///< LB zₜ 检验是否通过
    bool passes_no_arch_effect;                   ///< LB z²ₜ 检验是否通过
    Size jb_bootstrap_reps;                       ///< JB Bootstrap 次数 (0 = 渐近)
    std::string summary;                          ///< 诊断摘要
};

/// @brief G12 自适应 LB 滞后数: floor(log(T)) (自然对数, 下限 1, 上限 (T-1)/2)
inline Size diagnostics_auto_lag(Size n) {
    const Size m = static_cast<Size>(std::floor(std::log(static_cast<Real>(n))));
    return std::max<Size>(1, std::min<Size>(m, (n - 1) / 2));
}

/// @brief JB 参数化 Bootstrap p 值 (G-ADR4)
///
/// H0: z ~ iid N(0,1)。重复 B 次: 抽同长度 iid 正态样本 → 计算其 JB 统计量,
/// p = (1 + #{JB* ≥ JB_obs}) / (B + 1)。固定 seed, 同输入逐位可复现。
///
/// @param z       标准化残差 (长度 N ≥ 10)
/// @param B       Bootstrap 次数 (≥ 1)
/// @param seed    Philox 种子 (默认固定, 保证可复现)
inline Real jb_bootstrap_pvalue(const std::vector<Real>& z, Size B,
                                uint64_t seed = 0x5EED5EED5EED5EEDULL) {
    const Size n = z.size();
    if (B == 0) {
        throw std::invalid_argument("jb_bootstrap_pvalue: B must be >= 1");
    }

    // 观测 JB 统计量 (复用 Phase 7A 实现)
    const Real jb_obs = econometrics::jarque_bera_test(z).base.statistic;

    Philox4x64 rng(seed, 0x4A42ULL);  // stream 标识 JB-bootstrap
    const Real inv_2p53 = 1.0 / 9007199254740992.0;
    Size count = 0;
    std::vector<Real> zs(n);
    for (Size b = 0; b < B; ++b) {
        for (Size i = 0; i < n; i += 2) {
            // +0.5 nudge → u ∈ (0,1), 避免 box_muller 的 log(0)
            const Real u1 = (static_cast<Real>(rng() >> 11) + 0.5) * inv_2p53;
            const Real u2 = (static_cast<Real>(rng() >> 11) + 0.5) * inv_2p53;
            const auto [z1, z2] = box_muller(u1, u2);
            zs[i] = z1;
            if (i + 1 < n) zs[i + 1] = z2;
        }
        if (econometrics::jarque_bera_test(zs).base.statistic >= jb_obs) {
            ++count;
        }
    }
    return (static_cast<Real>(count) + 1.0) / (static_cast<Real>(B) + 1.0);
}

/// @brief 执行 GARCH 标准化残差诊断 (spec §2.0.6)
///
/// 三类检验: JB(zₜ) 正态性 (Bootstrap p 值)、LB(zₜ) 自相关、
/// LB(z²ₜ) ARCH 效应 (Li-Mak; GARCH 未充分捕捉条件异方差时拒绝)。
///
/// @param std_residuals 标准化残差 zₜ = εₜ/√hₜ (G11: 需调用方完成标准化)
/// @param n_bootstrap   JB Bootstrap 次数 (G-ADR4; 0 = 渐近 χ²(2))
/// @param lb_lag        LB 滞后数 (0 = 自适应 floor(log(T)), G12)
inline GarchDiagnosticsResult diagnose_garch_residuals(
    const std::vector<Real>& std_residuals,
    Size n_bootstrap = 1000,
    Size lb_lag = 0) {

    const Size n = std_residuals.size();
    if (n < 10) {
        throw std::invalid_argument(
            "diagnose_garch_residuals: need at least 10 observations");
    }
    for (Size t = 0; t < n; ++t) {
        if (!std::isfinite(std_residuals[t])) {
            throw std::invalid_argument(
                "diagnose_garch_residuals: non-finite value in std_residuals");
        }
    }

    // G12: 自适应滞后 (zₜ 与 z²ₜ 使用同一 lag, 保持子检验可比)
    const Size m = (lb_lag == 0) ? diagnostics_auto_lag(n) : lb_lag;
    if (m == 0 || m >= n) {
        throw std::invalid_argument("diagnose_garch_residuals: invalid lb_lag");
    }

    GarchDiagnosticsResult r;

    // zₜ 自相关
    r.lb_z = econometrics::ljung_box_test(std_residuals, m);

    // z²ₜ 自相关 (ARCH 效应, Li-Mak — G12 关键: z² 而非 z)
    std::vector<Real> z2(n);
    for (Size t = 0; t < n; ++t) z2[t] = std_residuals[t] * std_residuals[t];
    r.lb_z_squared = econometrics::ljung_box_test(z2, m);

    // JB 正态性 + Bootstrap p 值 (G-ADR4)
    r.jb_test = econometrics::jarque_bera_test(std_residuals);
    r.jb_bootstrap_reps = n_bootstrap;
    if (n_bootstrap > 0) {
        r.jb_test.base.p_value = jb_bootstrap_pvalue(std_residuals, n_bootstrap);
        r.jb_test.base.method_name = "Jarque-Bera (bootstrap)";
        r.jb_test.base.reject_null = (r.jb_test.base.p_value < 0.05);
    }

    // 诊断结论
    r.passes_normality = !r.jb_test.base.reject_null;
    r.passes_no_autocorr = !r.lb_z.base.reject_null;
    r.passes_no_arch_effect = !r.lb_z_squared.base.reject_null;

    // 摘要
    r.summary = std::string("GARCH diagnostics: JB ") +
                (r.passes_normality ? "pass" : "FAIL") +
                " (p=" + std::to_string(r.jb_test.base.p_value) + ")" +
                ", LB(z) " + (r.passes_no_autocorr ? "pass" : "FAIL") +
                " (lag=" + std::to_string(m) + ")" +
                ", LB(z^2) " + (r.passes_no_arch_effect ? "pass" : "FAIL");
    return r;
}

}  // namespace garch
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
