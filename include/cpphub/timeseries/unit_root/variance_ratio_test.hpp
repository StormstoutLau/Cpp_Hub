// =============================================================================
// variance_ratio_test.hpp - 方差比检验 (spec §3.5)
//
// Phase 7B v1.6 M2 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 教材锚点: Lo-MacKinlay 1988 / Chow-Denning 1993 / Campbell-Lo-MacKinlay 1997
// 对照库: Python arch 8.0.0 arch/unitroot/unitroot.py:1704-1765 VarianceRatio
//
// 幻觉点防护 (spec §6.2 U14-U17 + 本实现新实测):
//   U14-vr: VR(k) = sigma2_k / sigma2_1, sigma2_k 已含 1/k 因子, 不再除 k
//   U15: Z2 三重修正 — 前置因子 4 (非 2); delta_j 为 4 阶矩 (非 2 阶自协方差);
//        stat = sqrt(T)·(VR-1)/sqrt(theta), theta 是 O(1) 不含 1/T
//   U16: Chow-Denning 用 Z2 (异方差稳健) + SMM(m,inf) 联合分布 p=1-[2Phi(|CD|)-1]^m
//   U17: debiased 是 CLM 1997 重叠块偏差修正 (非 Chen-Deo 2006):
//        sigma2_1 *= nq/(nq-1); m = k·(nq-k+1)·(1-k/nq); sigma2_k *= nq·k/m
//        (等价 sigma2_k *= T^2/[(T-k+1)(T-k)])
//   U-log (spec Step 1.1 残留幻觉, 本实现以 arch 源码为准):
//        arch 用 delta_y = diff(y) 简单差分, 非对数收益率!
//        mu = (y[-1]-y[0])/(nobs-1) (telescoping 等价 sum(dy)/nq)
//   方差公式 (overlap=True, arch 默认):
//        同方差: var_h = 2(2k-1)(k-1)/(3k)            (CLM 2.4.39, 无 1/T)
//        稳健:   theta = sum_j 4(1-j/k)^2 · delta_j    (CLM 2.4.43)
//        delta_j = nq·sum_t z2_t·z2_{t-j} / (sum_t z2_t)^2,  z2_t = (dy_t-mu)^2
//   p 值 = 2 - 2·Phi(|stat|) (标准正态, arch 无 MacKinnon 表)
// =============================================================================
#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace unit_root {

// ---------------------------------------------------------------------------
// 结果 (spec §3.5 接口签名)
// ---------------------------------------------------------------------------
struct VarianceRatioResult {
    Real vr_statistic = 0.0;  ///< VR(k) = sigma2_k / sigma2_1 (含 debias 因子)
    Real z1_statistic = 0.0;  ///< Lo-MacKinlay Z1 (同方差, U14)
    Real z1_p_value = 0.0;
    Real z2_statistic = 0.0;  ///< Lo-MacKinlay Z2 (异方差稳健, U15)
    Real z2_p_value = 0.0;
    std::vector<Real> chow_denning_stats;  ///< 各 horizon 的 Z2 (U16)
    Real chow_denning_p_value = 0.0;       ///< SMM(m,inf) 联合 p 值 (U16-smm)
    Real clm_debiased_statistic = 0.0;     ///< CLM debiased Z2 (U17, 恒 debias)
    Real clm_debiased_p_value = 0.0;
    std::vector<Real> vr_list;   ///< 各 horizon 的 VR(k) (multi 时填充)
    std::vector<Size> k_values;  ///< horizon 列表 (multi 时填充)
    Size k = 0;                  ///< 主 horizon (single; multi 时为首元素)
    bool reject_null = false;    ///< H0: 随机游走 (5%, 基于 Z2/chow_denning)
    std::string summary;
};

// ---------------------------------------------------------------------------
// 内部核心: 单 horizon VR + Z1 + Z2 (arch _compute_statistic 逐行对照)
// ---------------------------------------------------------------------------
namespace detail {

struct VrCore {
    Real vr = 0.0;
    Real z1 = 0.0;
    Real z2 = 0.0;
    Size nq = 0;  ///< 收益样本量 = nobs - 1
};

inline VrCore vr_core(const std::vector<Real>& data, Size k, bool debiased,
                      const std::string& trend) {
    const Size nobs = data.size();
    // arch _check_specification: required = lags + 1 + int(debiased) (overlap)
    const Size required = k + 1 + (debiased ? 1 : 0);
    if (nobs < required) {
        throw std::invalid_argument(
            "variance_ratio_test: sample too small for horizon " +
            std::to_string(k));
    }
    const Real nq_r = static_cast<Real>(nobs - 1);

    // Step 1: 漂移 (U17-demean; trend="n" 零漂移)
    const Real mu = (trend == "n")
                        ? 0.0
                        : (data[nobs - 1] - data[0]) / static_cast<Real>(nobs - 1);

    // Step 2: 一阶差分 + 1 期方差 (简单差分, 非 log! U-log)
    std::vector<Real> dy(nobs - 1);
    for (Size t = 0; t < nobs - 1; ++t) dy[t] = data[t + 1] - data[t];
    Real sigma2_1 = 0.0;
    for (Real d : dy) sigma2_1 += (d - mu) * (d - mu);
    sigma2_1 /= nq_r;

    // k 期重叠块方差: dyq = y[t+k]-y[t], 分母 nq·k (U14-vr, 已含 1/k)
    Real sigma2_q = 0.0;
    for (Size t = 0; t + k < nobs; ++t) {
        const Real d = (data[t + k] - data[t]) - static_cast<Real>(k) * mu;
        sigma2_q += d * d;
    }
    sigma2_q /= nq_r * static_cast<Real>(k);

    // Step 5: CLM debiased (U17, 重叠块偏差修正)
    if (debiased) {
        sigma2_1 *= nq_r / (nq_r - 1.0);
        const Real m = static_cast<Real>(k) * (nq_r - static_cast<Real>(k) + 1.0) *
                       (1.0 - static_cast<Real>(k) / nq_r);
        sigma2_q *= nq_r * static_cast<Real>(k) / m;
    }

    if (!(sigma2_1 > 0.0)) {
        throw std::runtime_error(
            "variance_ratio_test: zero 1-period variance (constant series?)");
    }
    VrCore c;
    c.vr = sigma2_q / sigma2_1;
    c.nq = nobs - 1;

    // Z1 同方差 (CLM 2.4.39, 无 1/T 因子 — sqrt(nq) 在分子)
    const Real var_h =
        2.0 * (2.0 * static_cast<Real>(k) - 1.0) * (static_cast<Real>(k) - 1.0) /
        (3.0 * static_cast<Real>(k));
    c.z1 = std::sqrt(nq_r) * (c.vr - 1.0) / std::sqrt(var_h);

    // Z2 异方差稳健 (U15: 4 阶矩 delta_j + 前置 4 + sqrt(nq))
    Real scale = 0.0;
    std::vector<Real> z2v(nobs - 1);
    for (Size t = 0; t < nobs - 1; ++t) {
        z2v[t] = (dy[t] - mu) * (dy[t] - mu);
        scale += z2v[t];
    }
    scale *= scale;
    Real theta = 0.0;
    for (Size j = 1; j < k; ++j) {
        Real s = 0.0;
        for (Size t = j; t < nobs - 1; ++t) s += z2v[t] * z2v[t - j];
        const Real delta = nq_r * s / scale;
        const Real w = 1.0 - static_cast<Real>(j) / static_cast<Real>(k);
        theta += 4.0 * w * w * delta;
    }
    if (!(theta > 0.0)) {
        throw std::runtime_error(
            "variance_ratio_test: non-positive robust variance (theta <= 0)");
    }
    c.z2 = std::sqrt(nq_r) * (c.vr - 1.0) / std::sqrt(theta);
    return c;
}

// 两尾 p 值: 2 - 2·Phi(|stat|) (arch unitroot.py:1765)
inline Real vr_two_sided_p(Real stat) noexcept {
    return 2.0 - 2.0 * normal_cdf(std::fabs(stat));
}

}  // namespace detail

// ---------------------------------------------------------------------------
// 单 horizon 方差比检验 (U-ADR10, U14-U17)
//
// @param data 水平序列 (arch 语义: 内部一阶差分; 传价格即简单收益)
// @param k 检验 horizon (>= 2)
// @param use_debiased 是否用 CLM 小样本 debiased 修正 (默认 true, arch 默认)
// @param trend "c" (去漂移, arch 默认) / "n" (零漂移)
// ---------------------------------------------------------------------------
inline VarianceRatioResult variance_ratio_test(
    const std::vector<Real>& data, Size k = 2, bool use_debiased = true,
    const std::string& trend = "c") {
    if (k < 2) {
        throw std::invalid_argument("variance_ratio_test: k must be >= 2");
    }
    if (trend != "c" && trend != "n") {
        throw std::invalid_argument("variance_ratio_test: trend must be c or n");
    }

    VarianceRatioResult r;
    r.k = k;
    const detail::VrCore c = detail::vr_core(data, k, use_debiased, trend);
    r.vr_statistic = c.vr;
    r.z1_statistic = c.z1;
    r.z1_p_value = detail::vr_two_sided_p(c.z1);
    r.z2_statistic = c.z2;
    r.z2_p_value = detail::vr_two_sided_p(c.z2);

    // CLM debiased 参考统计量 (恒用 debias, 即 arch 默认 robust+debiased 组合)
    const detail::VrCore cd = use_debiased ? c : detail::vr_core(data, k, true, trend);
    r.clm_debiased_statistic = cd.z2;
    r.clm_debiased_p_value = detail::vr_two_sided_p(cd.z2);

    // 单 horizon Chow-Denning 退化为 m=1: p = 2 - 2·Phi(|Z2|) = z2_p
    r.chow_denning_stats = {c.z2};
    r.chow_denning_p_value = r.z2_p_value;
    r.reject_null = r.z2_p_value < 0.05;

    r.summary = "Variance-Ratio Test: H0 random walk (k=" + std::to_string(k) +
                ", trend=" + trend + (use_debiased ? ", debiased" : "") +
                "); VR=" + std::to_string(r.vr_statistic) +
                ", Z2=" + std::to_string(r.z2_statistic) +
                ", p=" + std::to_string(r.z2_p_value);
    return r;
}

// ---------------------------------------------------------------------------
// 多 horizon 方差比 + Chow-Denning 联合检验 (U16)
//
// @param k_list horizon 列表 (如 {2, 5, 10, 20}, 全部 >= 2, 非空)
// 注: CD 用 Z2 (异方差稳健); 联合 p 值 SMM(m,inf): 1 - [2·Phi(|CD|)-1]^m
// ---------------------------------------------------------------------------
inline VarianceRatioResult variance_ratio_test_multi(
    const std::vector<Real>& data, const std::vector<Size>& k_list = {2, 5, 10, 20},
    bool use_debiased = true, const std::string& trend = "c") {
    if (k_list.empty()) {
        throw std::invalid_argument("variance_ratio_test_multi: empty k_list");
    }
    if (trend != "c" && trend != "n") {
        throw std::invalid_argument("variance_ratio_test_multi: trend must be c or n");
    }

    VarianceRatioResult r;
    r.k = k_list[0];
    r.k_values = k_list;
    r.chow_denning_stats.reserve(k_list.size());
    r.vr_list.reserve(k_list.size());

    Real cd_abs = 0.0;
    for (Size k : k_list) {
        const detail::VrCore c = detail::vr_core(data, k, use_debiased, trend);
        r.chow_denning_stats.push_back(c.z2);
        r.vr_list.push_back(c.vr);
        if (k == k_list[0]) {  // 首 horizon 填充标量字段
            r.vr_statistic = c.vr;
            r.z1_statistic = c.z1;
            r.z1_p_value = detail::vr_two_sided_p(c.z1);
            r.z2_statistic = c.z2;
            r.z2_p_value = detail::vr_two_sided_p(c.z2);
        }
        if (std::fabs(c.z2) > cd_abs) cd_abs = std::fabs(c.z2);
    }

    // CLM debiased 参考统计量 (首 horizon, 恒 debias)
    const detail::VrCore cd0 = detail::vr_core(data, k_list[0], true, trend);
    r.clm_debiased_statistic = cd0.z2;
    r.clm_debiased_p_value = detail::vr_two_sided_p(cd0.z2);

    // SMM(m, inf) 联合 p 值 (U16-smm)
    const Real one_sided = 2.0 * detail::normal_cdf(cd_abs) - 1.0;
    Real powm = 1.0;
    for (Size i = 0; i < k_list.size(); ++i) powm *= one_sided;
    r.chow_denning_p_value = 1.0 - powm;
    r.reject_null = r.chow_denning_p_value < 0.05;

    r.summary = "Chow-Denning joint VR test: H0 random walk (m=" +
                std::to_string(k_list.size()) + " horizons, Z2-based); CD=" +
                std::to_string(cd_abs) + ", p=" +
                std::to_string(r.chow_denning_p_value);
    return r;
}

}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
