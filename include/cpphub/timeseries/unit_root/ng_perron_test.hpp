// =============================================================================
// ng_perron_test.hpp - Ng-Perron 2001 M 检验族 + MAIC 滞后选择 (spec §2.1)
//
// Phase 7C v1.7 M0 (PHASE7C_SPEC.md v1.2)
//
// 教材锚点: Ng-Perron 2001 Econometrica 69(6):1519-1554
//   (DOI 10.1111/1468-0262.00256; ⚠️ 非 J.Econometrics, NP1)
//   原文公式 2026-08-17 PDF 复核 (AU wp 副本, Section 2/3.1 逐字):
//     MZα = (T⁻¹y_T² − s²_AR)·(2·T⁻²Σỹ²_{t−1})⁻¹
//     MSB = (T⁻²Σỹ²_{t−1}/s²_AR)^{1/2};  MZt = MZα×MSB (恒等式)
//     MPT p=0: [c̄²T⁻²Σỹ² − c̄·T⁻¹ỹ_T²]/s²_AR   (c̄=−7 ⇒ 末项 +7)
//     MPT p=1: [c̄²T⁻²Σỹ² + (1−c̄)T⁻¹ỹ_T²]/s²_AR (c̄=−13.5 ⇒ 末项 +14.5, NP4)
//     s²_AR = σ̂²(k)/(1−β̂(1))²                (AR 谱, 对 Δỹ 拟合, NP3)
//
// 对照生态 (NP6): 无成熟开源库输出 M 族 (statsmodels/arch/urca 均无);
//   Stata dfgls 可编程访问未验证 (readiness C6) ⇒ 本实现基准 =
//   原文公式钉死 (上) + 恒等式自检 (1e-12) + 文献 Table 1 临界值精确相等 +
//   模拟方向断言; Stata 可用后补逐 k MAIC 1e-10 对照 (verify_np_stata.py 占位)
//
// 口径决策 (实施记录, Stata 验证前冻结):
//   - 固定样本: 所有 k 的辅助回归统一 t = k_max+1..T (1-based), 观测 n = T−k_max;
//     σ̂²(k) = SSR(k)/(n−1) (spec §2.1 Step2.2, 4 源裁决 NP2)
//   - Σỹ² 与 τ_T 同口径 (spec Step2.3): Σ_{t=k_max+1}^{T} ỹ²_{t−1}
//   - 统计量 T 因子取 T_eff = n (固定样本一致性; ỹ_T = 全样本最后观测)
//   - MBIC/seq-t: 冻结签名无输出字段 (v1.2 §1.4), 推迟; MAIC 为主选择准则
//
// 幻觉点防护 (spec §9.5):
//   NP1 (Econometrica) / NP2 (τ_T 4 源裁决, λ̂−λ̃ 形零命中) /
//   NP3 (AR 谱对 Δỹ 差分拟合) / NP4 (MPT 分情形 +7/+14.5) / NP5 (方向+表) /
//   NP6 (Stata 非基准)
//
// B1 复用处置 (readiness 预批): ERS GLS 去势重实现 (df_gls_test.hpp 内联不可
//   include; 提取公共函数 = 触碰 v1.6 文件, 违 §8.5 additive-only)
// =============================================================================
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/np_tables.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace unit_root {

/// Ng-Perron 检验结果 (spec §2.1 v1.2 冻结签名, 默认初始化)
struct NgPerronResult {
    // 四统计量 (NP 2001 §3.1, GLS 去势后)
    Real mz_alpha = 0.0;  ///< MZα = (T⁻¹ỹ_T² − s²)/(2·T⁻²Σỹ²)
    Real mz_t = 0.0;      ///< MZt = MZα × MSB (恒等式, 自检 1e-12)
    Real msb = 0.0;       ///< MSB = √(T⁻²Σỹ²/s²)
    Real mpt = 0.0;       ///< MPT 分情形 (NP4): c̄=−7 末项 +7 / c̄=−13.5 末项 +14.5
    // 临界值与结论 — 顺序 {MZα, MZt, MSB, MPT}; std::array (§1.4-2)
    std::array<Real, 4> cv_1pct{};
    std::array<Real, 4> cv_5pct{};
    std::array<Real, 4> cv_10pct{};
    std::array<bool, 4> reject_5pct{};  ///< NP5: 统一 stat < cv 拒绝
    // 滞后选择
    Size selected_lag = 0;         ///< MAIC 最优 k*
    std::vector<Real> maic;        ///< 逐 k MAIC(k) 轨迹 (Stata dfgls 对照点)
    std::vector<Real> sigma2_k;    ///< 逐 k σ̂²(k) (对照点, 1e-10 层)
    std::string trend_spec;        ///< "c" / "ct"
    std::string summary;
};

namespace detail {

/// ERS GLS 去势 (B1 处置: 重实现, 与 df_gls_test.hpp 内联变换逐步一致;
/// 公开以便测试与去势等价性对照)
/// @param data 原始序列 (0-based, T 个)
/// @param trend_spec "c" (c̄=−7, z=[1]) / "ct" (c̄=−13.5, z=[1,t])
/// @return 去势序列 ỹ = y − z·ψ̂ (ψ̂ 最小化准差分平方和)
inline std::vector<Real> gls_detrend(const std::vector<Real>& data,
                                     const std::string& trend_spec) {
    const Size T = data.size();
    const Size tc = (trend_spec == "c") ? 1u : 2u;
    const Real c_bar = (trend_spec == "c") ? -7.0 : -13.5;
    const Real rho = 1.0 + c_bar / static_cast<Real>(T);

    // z 趋势列 (c → [1]; ct → [1, t])
    std::vector<std::vector<Real>> z(T, std::vector<Real>(tc, 1.0));
    if (tc == 2) {
        for (Size t = 0; t < T; ++t) z[t][1] = static_cast<Real>(t);
    }
    // GLS 准差分 (第一项不变换, U9)
    std::vector<std::vector<Real>> dz(T, std::vector<Real>(tc, 0.0));
    dz[0] = z[0];
    for (Size t = 1; t < T; ++t) {
        for (Size j = 0; j < tc; ++j) dz[t][j] = z[t][j] - rho * z[t - 1][j];
    }
    std::vector<Real> dy(T);
    dy[0] = data[0];
    for (Size t = 1; t < T; ++t) dy[t] = data[t] - rho * data[t - 1];

    const auto fit = ols_fit(dy, dz);
    std::vector<Real> yd(T);
    for (Size t = 0; t < T; ++t) {
        Real pred = 0.0;
        for (Size j = 0; j < tc; ++j) pred += fit.beta[j] * z[t][j];
        yd[t] = data[t] - pred;
    }
    return yd;
}

}  // namespace detail

/// @brief Ng-Perron M 检验族 (MZα/MZt/MSB/MPT) + MAIC 滞后选择
///
/// H0: 单位根 (α=1); NP5 拒绝方向: 四统计量统一 statistic < cv
///
/// @param data 原始序列 (水平值, 无 NaN)
/// @param trend_spec "c" (常数, c̄=−7) / "ct" (趋势, c̄=−13.5)
/// @param max_lag 0 => Schwert 上限; >0 => 用户固定 k_max (固定样本口径)
/// @return NgPerronResult (含逐 k MAIC/σ̂² 轨迹)
/// @throws std::invalid_argument 样本过小/含 NaN/零方差/trend 非法/lag 过大
inline NgPerronResult ng_perron_test(const std::vector<Real>& data,
                                     const std::string& trend_spec = "ct",
                                     Size max_lag = 0) {
    const Size T = data.size();
    if (T < 10) {
        throw std::invalid_argument("ng_perron_test: sample too small (T < 10)");
    }
    if (trend_spec != "c" && trend_spec != "ct") {
        throw std::invalid_argument("ng_perron_test: trend_spec must be c/ct");
    }
    Real vmin = data[0], vmax = data[0];
    for (Real v : data) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("ng_perron_test: NaN/Inf in data");
        }
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }
    if (vmax == vmin) {
        throw std::invalid_argument("ng_perron_test: zero variance (constant)");
    }

    // Step 1: GLS 去势 (B1 重实现 ERS 变换)
    const std::vector<Real> yd = detail::gls_detrend(data, trend_spec);
    const Real c_bar = (trend_spec == "c") ? -7.0 : -13.5;

    // Step 2: 固定样本辅助回归 (逐 k = 0..k_max)
    //   坐标勘误 (2026-08-18 实施期发现, 越界实锤): 原冻结口径
    //   "t = k_max+1..T (1-based), n = T−k_max" 在 Δỹ_t 索引上偏 1 —
    //   Δỹ_T 需 dyd[T−1] 越界 (dyd 合法 0..T−2), 堆越界读 ⇒ 跨进程概率性
    //   mz_alpha 漂移 (主控站 3/8 复现, B 站 GCC 偶未触发);
    //   且最深滞后 Δỹ_{t−k_max} 在 t=k_max+1 处需 Δỹ₁ (依赖不存在的 y₀)。
    //   修正口径 (与 statsmodels adfuller 族一致): 回归观测取 0-based
    //   Δỹ 下标 idx ∈ [k_max, T−2], n = T−1−k_max; 最深滞后 dyd[idx−k] ≥ 0 恒合法
    const Size k_max = (max_lag > 0) ? max_lag : schwert_lag(T);
    if (T < 2 || k_max + 3 > T) {
        throw std::invalid_argument(
            "ng_perron_test: lag too large for sample (n < 3)");
    }
    const Size n = T - 1 - k_max;  // 固定样本观测数 (修正口径)

    // Σỹ² (τ_T 与统计量同口径): 回归 ỹ_{idx} 列 (= 1-based ỹ_{t−1}) 平方和,
    //   0-based idx ∈ [k_max, T−2]
    Real sum_y2 = 0.0;
    for (Size s = k_max; s + 2 <= T; ++s) {  // s = k_max .. T−2
        sum_y2 += yd[s] * yd[s];
    }
    const Real n_r = static_cast<Real>(n);
    const Real t_eff = n_r;  // 口径决策: T 因子 = 固定样本观测数

    // 预备 Δỹ (dyd[j] = ỹ_{j+1} − ỹ_j, j = 0..T−2)
    std::vector<Real> dyd(T - 1);
    for (Size t = 1; t < T; ++t) dyd[t - 1] = yd[t] - yd[t - 1];

    std::vector<Real> maic(k_max + 1);
    std::vector<Real> sigma2(k_max + 1);
    std::vector<std::vector<Real>> betas(k_max + 1);  // 逐 k 系数 (s² 用 k*)

    for (Size k = 0; k <= k_max; ++k) {
        // 回归: Δỹ[idx] = β₀·ỹ[idx] + Σ_{j=1..k} β_j·Δỹ[idx−j]
        //   idx (0-based Δỹ 下标) ∈ [k_max, T−2]; ỹ 列与 Δỹ 同时刻基 (adfuller 族)
        std::vector<Real> lhs(n);
        std::vector<std::vector<Real>> X(n, std::vector<Real>(1 + k));
        for (Size i = 0; i < n; ++i) {
            const Size idx = k_max + i;   // 0-based Δỹ 下标
            lhs[i] = dyd[idx];            // Δỹ_{idx+1} (1-based t = idx+1)
            X[i][0] = yd[idx];            // ỹ_idx = 1-based ỹ_{t−1} ✓
            for (Size j = 1; j <= k; ++j) {
                X[i][j] = dyd[idx - j];   // Δỹ_{idx+1−j} (1-based), idx−j ≥ 0 ✓
            }
        }
        Real ssr = 0.0;
        std::vector<Real> beta(1 + k, 0.0);
        try {
            const auto fit = detail::ols_fit(lhs, X);
            ssr = fit.ssr;
            beta = fit.beta;
        } catch (const std::invalid_argument&) {
            // 奇异 (共线): 大 σ̂² ⇒ MAIC 排除该 k
            sigma2[k] = std::numeric_limits<Real>::infinity();
            maic[k] = std::numeric_limits<Real>::infinity();
            betas[k] = beta;
            continue;
        }
        sigma2[k] = ssr / (n_r - 1.0);  // spec Step2.2: SSR/(T−k_max−1)
        const Real tau_T =
            beta[0] * beta[0] * sum_y2 / sigma2[k];  // NP2 4 源裁决形
        maic[k] = std::log(sigma2[k]) + 2.0 * (tau_T + static_cast<Real>(k)) / n_r;
        betas[k] = beta;
    }

    // Step 2.5: k* = argmin MAIC (首个)
    Size k_star = 0;
    Real best = std::numeric_limits<Real>::infinity();
    for (Size k = 0; k <= k_max; ++k) {
        if (maic[k] < best) {
            best = maic[k];
            k_star = k;
        }
    }

    // Step 3: 长期方差 s² (NP3: AR 谱对 Δỹ, k* 系数)
    Real sum_beta = 0.0;
    for (Size j = 1; j <= k_star; ++j) sum_beta += betas[k_star][j];
    const Real s2 = sigma2[k_star] / ((1.0 - sum_beta) * (1.0 - sum_beta));
    if (!(s2 > 0.0) || !std::isfinite(s2)) {
        throw std::runtime_error("ng_perron_test: invalid spectral density");
    }

    // Step 4: 四统计量 (口径: T_eff = n, ỹ_T = 最后观测, Σỹ² 固定样本)
    const Real yT = yd[T - 1];
    NgPerronResult res;
    res.mz_alpha = (yT * yT / t_eff - s2) / (2.0 * sum_y2 / (t_eff * t_eff));
    res.msb = std::sqrt(sum_y2 / (t_eff * t_eff) / s2);
    res.mz_t = res.mz_alpha * res.msb;  // 恒等式 (测试断言 1e-12)
    // NP4 分情形末项: p=0 → −c̄ = +7; p=1 → (1−c̄) = +14.5
    const Real tail =
        (trend_spec == "c") ? -c_bar : (1.0 - c_bar);
    res.mpt =
        (c_bar * c_bar * sum_y2 / (t_eff * t_eff) + tail * yT * yT / t_eff) / s2;

    // Step 5: 临界值 (np_tables 双源转录) + NP5 统一拒绝方向
    const std::array<Real, 4> stats{res.mz_alpha, res.mz_t, res.msb, res.mpt};
    for (Size s = 0; s < 4; ++s) {
        const auto st = static_cast<NPStat>(s);
        res.cv_1pct[s] = np_critical_value(st, trend_spec, 0.01);
        res.cv_5pct[s] = np_critical_value(st, trend_spec, 0.05);
        res.cv_10pct[s] = np_critical_value(st, trend_spec, 0.10);
        res.reject_5pct[s] = stats[s] < res.cv_5pct[s];
    }

    res.selected_lag = k_star;
    res.maic = std::move(maic);
    res.sigma2_k = std::move(sigma2);
    res.trend_spec = trend_spec;
    res.summary = "Ng-Perron M tests (trend=" + trend_spec + ", k*=" +
                  std::to_string(k_star) + " MAIC): H0: unit root; "
                  "MZa=" + std::to_string(res.mz_alpha) +
                  " MZt=" + std::to_string(res.mz_t) +
                  " MSB=" + std::to_string(res.msb) +
                  " MPT=" + std::to_string(res.mpt);
    return res;
}

}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
