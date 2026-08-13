// =============================================================================
// risk_diagnostics.hpp - VaR/ES 后验诊断 (P0-P2)
//
// Phase 7A Wave 1: 模块特定证伪统计量
//
// 包含 4 个检验:
//   1. McNeil-Frey 动态量化检验 (Engle-Manganelli 2004 DQ)
//   2. Berkowitz 尾部检验 (Berkowitz 2001)
//   3. MC 收敛性诊断 (batch means)
//   4. ES 后验检验 (Expected Shortfall backtesting)
//
// ADR-015 决策点 4: 命名空间维持 cpphub::v1 (与 backtesting.hpp 一致)
// 依赖: core/math.hpp (normal_cdf, inv_normal_cdf), special_functions.hpp (chi2_sf, beta_i)
//       detail/ols_simple.hpp (DQ 辅助回归)
//
// 教材锚点: Engle-Manganelli 2004, Berkowitz 2001, McNeil-Frey-Embrechts 2005
// 排幻觉点: H16(DQ含VaR)/H17(Berkowitz需模型CDF)/H18(MC批次均值)/H19(ES条件超越)
// =============================================================================
#pragma once

#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/econometrics/core/special_functions.hpp"
#include "cpphub/econometrics/inference/detail/ols_simple.hpp"

namespace cpphub {
inline namespace v1 {

// =============================================================================
// McNeil-Frey 动态量化检验 (Engle-Manganelli 2004)
// =============================================================================

// DQ = Hit'X(X'X)^{-1}X'Hit / (π(1-π)) ~ χ²(q)
// 排幻觉点 H16: X_t 包含 [1, Hit_{t-1}, ..., VaR_t], 非仅 Hit_{t-1}
struct DynamicQuantileResult {
    Real dq_statistic;
    Size df;
    Real p_value;
    bool reject_correct_coverage;
};

inline DynamicQuantileResult dynamic_quantile_test(
    const std::vector<Real>& returns,
    const std::vector<Real>& var_forecasts,
    Real confidence_level,
    const std::vector<Size>& hit_lags = {1, 2, 3, 4}) {

    const Size n = returns.size();
    if (n < 5) {
        throw std::invalid_argument("dynamic_quantile_test: need at least 5 observations");
    }
    if (var_forecasts.size() != n) {
        throw std::invalid_argument("dynamic_quantile_test: returns and var_forecasts size mismatch");
    }
    if (hit_lags.empty()) {
        throw std::invalid_argument("dynamic_quantile_test: hit_lags is empty");
    }

    const Real pi = 1.0 - confidence_level;  // 超越概率 (如 95% VaR → π=0.05)
    if (pi <= 0.0 || pi >= 1.0) {
        throw std::invalid_argument("dynamic_quantile_test: invalid confidence_level");
    }

    // 构造 Hit_t (超越指示变量)
    // Hit_t = 1 if r_t < -VaR_t else 0
    std::vector<Real> hit(n);
    for (Size t = 0; t < n; ++t) {
        hit[t] = (returns[t] < -var_forecasts[t]) ? 1.0 : 0.0;
    }

    // 最大滞后阶
    Size max_lag = 0;
    for (Size lag : hit_lags) {
        if (lag > max_lag) max_lag = lag;
    }
    if (max_lag >= n) {
        throw std::invalid_argument("dynamic_quantile_test: max_lag >= N");
    }

    // 构造辅助回归设计矩阵 X (排幻觉点 H16)
    // X_t = [1, Hit_{t-1}, Hit_{t-2}, ..., VaR_t]
    // 有效行数 = n - max_lag
    const Size n_eff = n - max_lag;
    const Size k_lags = hit_lags.size();
    const Size p_cols = 1 + k_lags + 1;  // 常数 + hit lags + VaR_t

    if (n_eff <= p_cols) {
        throw std::invalid_argument("dynamic_quantile_test: insufficient observations after lag");
    }

    // 特判: Hit 全为 0 或全为 1 → 设计矩阵共线 (Hit_lag 列 = 常数列)
    // 此时 OLS 无解, 但 DQ 有明确含义:
    //   - Hit 全为 0: 无超越, fitted=0, DQ=0, p_value=1 (不拒绝)
    //   - Hit 全为 1: 全超越, fitted≈1, DQ=n_eff/(π(1-π)) (应拒绝)
    Size hit_sum = 0;
    for (Size t = 0; t < n; ++t) hit_sum += static_cast<Size>(hit[t]);
    if (hit_sum == 0 || hit_sum == n) {
        DynamicQuantileResult result;
        result.df = p_cols;
        if (hit_sum == 0) {
            result.dq_statistic = 0.0;
            result.p_value = 1.0;
            result.reject_correct_coverage = false;
        } else {
            // Hit 全为 1: fitted = 1 (常数拟合), DQ = n_eff / (π(1-π))
            result.dq_statistic = static_cast<Real>(n_eff) / (pi * (1.0 - pi));
            result.p_value = cpphub::v1::econometrics::detail::chi2_sf(
                static_cast<Real>(p_cols), result.dq_statistic);
            result.reject_correct_coverage = (result.p_value < 0.05);
        }
        return result;
    }

    std::vector<std::vector<Real>> X(n_eff, std::vector<Real>(p_cols));
    std::vector<Real> y(n_eff);

    for (Size i = 0; i < n_eff; ++i) {
        const Size t = i + max_lag;  // 原始时间索引
        Size col = 0;
        X[i][col++] = 1.0;  // 常数项
        for (Size lag : hit_lags) {
            X[i][col++] = hit[t - lag];  // Hit_{t-lag}
        }
        X[i][col++] = var_forecasts[t];  // VaR_t (排幻觉点 H16)
        y[i] = hit[t];
    }

    // OLS: Hit on X → fitted values
    std::vector<Real> fitted, resid;
    Real r2;
    try {
        econometrics::detail::ols_simple(y, X, fitted, resid, r2);
    } catch (const std::runtime_error&) {
        // 设计矩阵奇异 (如 VaR 恒定导致共线), 降级为仅用常数+Hit_lags
        // 移除 VaR 列, 重新构造 X
        const Size p_cols_no_var = 1 + k_lags;
        std::vector<std::vector<Real>> X2(n_eff, std::vector<Real>(p_cols_no_var));
        for (Size i = 0; i < n_eff; ++i) {
            const Size t = i + max_lag;
            Size col = 0;
            X2[i][col++] = 1.0;
            for (Size lag : hit_lags) {
                X2[i][col++] = hit[t - lag];
            }
        }
        econometrics::detail::ols_simple(y, X2, fitted, resid, r2);
        // 降级后 df 减少 1 (移除 VaR 列)
        const Real dq = [&] {
            Real s = 0.0;
            for (Size i = 0; i < n_eff; ++i) s += fitted[i] * fitted[i];
            return s / (pi * (1.0 - pi));
        }();
        DynamicQuantileResult result;
        result.dq_statistic = dq;
        result.df = p_cols_no_var;
        result.p_value = cpphub::v1::econometrics::detail::chi2_sf(
            static_cast<Real>(p_cols_no_var), dq);
        result.reject_correct_coverage = (result.p_value < 0.05);
        return result;
    }

    // DQ = Σ fitted_t² / (π(1-π)) = Hit'X(X'X)^{-1}X'Hit / (π(1-π))
    Real ss_fitted = 0.0;
    for (Size i = 0; i < n_eff; ++i) {
        ss_fitted += fitted[i] * fitted[i];
    }
    const Real dq = ss_fitted / (pi * (1.0 - pi));
    const Size df = p_cols;
    const Real p_value = cpphub::v1::econometrics::detail::chi2_sf(
        static_cast<Real>(df), dq);

    DynamicQuantileResult result;
    result.dq_statistic = dq;
    result.df = df;
    result.p_value = p_value;
    result.reject_correct_coverage = (p_value < 0.05);
    return result;
}

// =============================================================================
// Berkowitz 尾部检验 (Berkowitz 2001)
// =============================================================================

// z_t = Φ⁻¹(F_model(r_t)), 检验 z ~ N(0,1) 且无自相关
// LR = LR_mean + LR_var + LR_autocorr ~ χ²(3)
// 排幻觉点 H17: 需模型 CDF, 非经验 CDF
struct BerkowitzResult {
    Real lr_statistic;    // ~ χ²(3)
    Real lr_mean;         // 均值检验分量 ~ χ²(1)
    Real lr_variance;     // 方差检验分量 ~ χ²(1)
    Real lr_autocorr;     // 自相关检验分量 ~ χ²(1)
    Real p_value;
    bool reject_correct_distribution;
};

inline BerkowitzResult berkowitz_test(
    const std::vector<Real>& returns,
    std::function<Real(Real)> model_cdf,
    Size lag = 1) {

    const Size n = returns.size();
    if (n < 5) {
        throw std::invalid_argument("berkowitz_test: need at least 5 observations");
    }
    if (!model_cdf) {
        throw std::invalid_argument("berkowitz_test: model_cdf is empty");
    }

    // z_t = Φ⁻¹(F_model(r_t)) (排幻觉点 H17: 模型 CDF, 非经验 CDF)
    std::vector<Real> z(n);
    for (Size t = 0; t < n; ++t) {
        const Real u = model_cdf(returns[t]);
        if (u <= 0.0 || u >= 1.0) {
            throw std::runtime_error("berkowitz_test: model_cdf returned value outside (0,1)");
        }
        z[t] = inv_normal_cdf(u);
    }

    // 均值和方差 (有偏 /N)
    Real z_mean = 0.0;
    for (Size t = 0; t < n; ++t) z_mean += z[t];
    z_mean /= static_cast<Real>(n);

    Real z_var = 0.0;
    for (Size t = 0; t < n; ++t) {
        const Real d = z[t] - z_mean;
        z_var += d * d;
    }
    z_var /= static_cast<Real>(n);

    if (z_var <= 0.0) {
        throw std::runtime_error("berkowitz_test: zero variance in z");
    }

    // lag-1 自相关 (如果 lag >= 1)
    Real rho1 = 0.0;
    if (lag >= 1 && n > 1) {
        Real gamma0 = z_var * static_cast<Real>(n);  // Σ(z-z̄)²
        Real gamma1 = 0.0;
        for (Size t = 1; t < n; ++t) {
            gamma1 += (z[t] - z_mean) * (z[t - 1] - z_mean);
        }
        rho1 = gamma1 / gamma0;
    }

    // LR 分量 (Berkowitz 2001)
    // LR_mean = N·z̄² (假设 σ²=1, 检验 μ=0)
    const Real lr_mean = static_cast<Real>(n) * z_mean * z_mean;

    // LR_var = N·(s² - 1 - log(s²)) (检验 σ²=1)
    const Real lr_var = static_cast<Real>(n) * (z_var - 1.0 - std::log(z_var));

    // LR_autocorr = (N-1)·ρ₁² (检验无自相关)
    Real lr_autocorr = 0.0;
    if (lag >= 1) {
        lr_autocorr = static_cast<Real>(n - 1) * rho1 * rho1;
    }

    // LR_total ~ χ²(3) (如果 lag >= 1), 否则 ~ χ²(2)
    const Real lr_total = lr_mean + lr_var + lr_autocorr;
    const Size df = (lag >= 1) ? 3 : 2;
    const Real p_value = cpphub::v1::econometrics::detail::chi2_sf(
        static_cast<Real>(df), lr_total);

    BerkowitzResult result;
    result.lr_statistic = lr_total;
    result.lr_mean = lr_mean;
    result.lr_variance = lr_var;
    result.lr_autocorr = lr_autocorr;
    result.p_value = p_value;
    result.reject_correct_distribution = (p_value < 0.05);
    return result;
}

// =============================================================================
// MC 收敛性诊断 (batch means)
// =============================================================================

// 原理: MC 标准误差应随 √N 衰减
// 排幻觉点 H18: 标准误差用批次均值法 (batch means), 非单次估计
struct MCConvergenceResult {
    Real estimated_std_error;    // 估计的 MC 标准误差
    Size n_paths;
    Size n_batches;              // 批次数
    std::vector<Real> batch_means;
    Real convergence_rate;       // 应接近 1/√2 ≈ 0.707 (全样本 SE / 半样本 SE)
    bool converged;              // SE < tolerance
};

inline MCConvergenceResult mc_convergence_diagnosis(
    const std::vector<Real>& mc_estimates,
    Size n_batches = 20,
    Real tolerance = 1e-4) {

    const Size n = mc_estimates.size();
    if (n < 10) {
        throw std::invalid_argument("mc_convergence_diagnosis: need at least 10 estimates");
    }
    if (n_batches < 2 || n_batches > n) {
        throw std::invalid_argument("mc_convergence_diagnosis: invalid n_batches");
    }

    // 批次大小 (每批 M 个样本)
    const Size batch_size = n / n_batches;
    if (batch_size < 1) {
        throw std::invalid_argument("mc_convergence_diagnosis: batch_size < 1");
    }

    // 计算批次均值 (排幻觉点 H18)
    std::vector<Real> batch_means(n_batches);
    for (Size b = 0; b < n_batches; ++b) {
        Real sum = 0.0;
        const Size start = b * batch_size;
        const Size end = (b == n_batches - 1) ? n : start + batch_size;
        const Size count = end - start;
        for (Size i = start; i < end; ++i) {
            sum += mc_estimates[i];
        }
        batch_means[b] = sum / static_cast<Real>(count);
    }

    // 总均值
    Real grand_mean = 0.0;
    for (Size b = 0; b < n_batches; ++b) grand_mean += batch_means[b];
    grand_mean /= static_cast<Real>(n_batches);

    // 批次均值的样本方差 (无偏, /(B-1))
    Real batch_var = 0.0;
    for (Size b = 0; b < n_batches; ++b) {
        const Real d = batch_means[b] - grand_mean;
        batch_var += d * d;
    }
    batch_var /= static_cast<Real>(n_batches - 1);

    // SE(x̄) = sqrt(s²_batch / B) (排幻觉点 H18: 批次均值法)
    const Real se_full = std::sqrt(batch_var / static_cast<Real>(n_batches));

    // 收敛率: 全样本 SE / 半样本 SE, 应接近 1/√2 ≈ 0.707
    // 半样本: 用前 n/2 个估计值, n_batches/2 批
    Real convergence_rate = std::sqrt(0.5);  // 默认理论值 1/√2
    const Size n_half = n / 2;
    const Size n_batches_half = n_batches / 2;
    if (n_half >= 10 && n_batches_half >= 2) {
        const Size batch_size_half = n_half / n_batches_half;
        if (batch_size_half >= 1) {
            std::vector<Real> batch_means_half(n_batches_half);
            for (Size b = 0; b < n_batches_half; ++b) {
                Real sum = 0.0;
                const Size start = b * batch_size_half;
                const Size end = (b == n_batches_half - 1) ? n_half : start + batch_size_half;
                const Size count = end - start;
                for (Size i = start; i < end; ++i) {
                    sum += mc_estimates[i];
                }
                batch_means_half[b] = sum / static_cast<Real>(count);
            }
            Real grand_mean_half = 0.0;
            for (Size b = 0; b < n_batches_half; ++b) grand_mean_half += batch_means_half[b];
            grand_mean_half /= static_cast<Real>(n_batches_half);
            Real batch_var_half = 0.0;
            for (Size b = 0; b < n_batches_half; ++b) {
                const Real d = batch_means_half[b] - grand_mean_half;
                batch_var_half += d * d;
            }
            batch_var_half /= static_cast<Real>(n_batches_half - 1);
            const Real se_half = std::sqrt(batch_var_half / static_cast<Real>(n_batches_half));
            if (se_half > 0.0) {
                convergence_rate = se_full / se_half;
            }
        }
    }

    MCConvergenceResult result;
    result.estimated_std_error = se_full;
    result.n_paths = n;
    result.n_batches = n_batches;
    result.batch_means = batch_means;
    result.convergence_rate = convergence_rate;
    result.converged = (se_full < tolerance);
    return result;
}

// =============================================================================
// ES 后验检验 (Expected Shortfall backtesting)
// =============================================================================

// 原理: 超越条件均值应等于 ES 预测
// 排幻觉点 H19: ES 后验需在 VaR 超越条件下, 非全样本
//   ES_realized = mean(r_t | r_t < -VaR_t), 应接近 ES_forecast
struct ESBacktestResult {
    Real es_forecast_mean;
    Real es_realized_mean;      // 超越条件下的实际均值
    Real bias;                  // es_realized - es_forecast
    Size n_violations;
    Real t_stat;                // (es_realized - es_forecast) / SE
    Real p_value;               // 双侧 t 检验 p 值
    bool reject_correct_es;
};

inline ESBacktestResult es_backtest(
    const std::vector<Real>& returns,
    const std::vector<Real>& var_forecasts,
    const std::vector<Real>& es_forecasts) {

    const Size n = returns.size();
    if (n < 5) {
        throw std::invalid_argument("es_backtest: need at least 5 observations");
    }
    if (var_forecasts.size() != n || es_forecasts.size() != n) {
        throw std::invalid_argument("es_backtest: size mismatch");
    }

    // 找到 VaR 超越 (排幻觉点 H19: 仅超越条件下)
    std::vector<Real> realized_exceedances;
    std::vector<Real> forecast_es_at_violations;
    for (Size t = 0; t < n; ++t) {
        if (returns[t] < -var_forecasts[t]) {
            realized_exceedances.push_back(returns[t]);
            forecast_es_at_violations.push_back(es_forecasts[t]);
        }
    }

    const Size n_viol = realized_exceedances.size();
    if (n_viol < 2) {
        // 不足 2 个超越, 无法计算 t 统计量
        ESBacktestResult result;
        result.es_forecast_mean = (n_viol == 1) ? forecast_es_at_violations[0] : 0.0;
        result.es_realized_mean = (n_viol == 1) ? realized_exceedances[0] : 0.0;
        result.bias = result.es_realized_mean - result.es_forecast_mean;
        result.n_violations = n_viol;
        result.t_stat = 0.0;
        result.p_value = 1.0;  // 无法拒绝
        result.reject_correct_es = false;
        return result;
    }

    // ES_realized = mean(r_t | r_t < -VaR_t) (排幻觉点 H19)
    Real es_realized = 0.0;
    for (Size i = 0; i < n_viol; ++i) es_realized += realized_exceedances[i];
    es_realized /= static_cast<Real>(n_viol);

    // ES_forecast_mean = mean(ES_forecast_t | t is violation)
    Real es_forecast = 0.0;
    for (Size i = 0; i < n_viol; ++i) es_forecast += forecast_es_at_violations[i];
    es_forecast /= static_cast<Real>(n_viol);

    const Real bias = es_realized - es_forecast;

    // SE = std(realized) / sqrt(n_viol) (样本标准差, 无偏 /N-1)
    Real var_realized = 0.0;
    for (Size i = 0; i < n_viol; ++i) {
        const Real d = realized_exceedances[i] - es_realized;
        var_realized += d * d;
    }
    var_realized /= static_cast<Real>(n_viol - 1);  // 无偏
    const Real se = std::sqrt(var_realized / static_cast<Real>(n_viol));

    // t 统计量 (双侧检验)
    Real t_stat = 0.0;
    Real p_value = 1.0;
    if (se > 0.0) {
        t_stat = bias / se;
        // 双侧 p_value: p = beta_i(df/2, 0.5, df/(df+t²))
        // 其中 df = n_viol - 1
        const Real df = static_cast<Real>(n_viol - 1);
        const Real x = df / (df + t_stat * t_stat);
        p_value = cpphub::v1::econometrics::detail::beta_i(df / 2.0, 0.5, x);
    }

    ESBacktestResult result;
    result.es_forecast_mean = es_forecast;
    result.es_realized_mean = es_realized;
    result.bias = bias;
    result.n_violations = n_viol;
    result.t_stat = t_stat;
    result.p_value = p_value;
    result.reject_correct_es = (p_value < 0.05);
    return result;
}

}  // namespace v1
}  // namespace cpphub
