// =============================================================================
// pricing_diagnostics.hpp - 定价模型诊断 (P2)
//
// Phase 7A Wave 3c: IV 微笑拟合优度 + 价格残差诊断
//
// 教材锚点:
//   - Gatheral 2006 (SVI, IV surface)
//   - Fengler 2009 (IV surface 动态)
//
// ADR-015 方案 B: 仅依赖 core/, 不依赖 Eigen3
//   χ² p 值用 econometrics/core/special_functions.hpp::chi2_sf
//   t 检验 p 值用 std::erfc (正态近似, 大样本)
//
// 排幻觉点:
//   H22: 权重用市场 IV 的 Bid-Ask 宽度 (1/spread²), 非简单方差
//        Bid-Ask 宽度反映市场对该行权价/到期日的定价不确定性
// =============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/special_functions.hpp"

namespace cpphub {
inline namespace v1 {

// =============================================================================
// IVFitGoodnessResult - IV 微笑拟合优度检验结果
// =============================================================================
struct IVFitGoodnessResult {
    Real chi_squared;
    Size degrees_of_freedom;
    Real p_value;
    Real rmse;                   // 均方根误差
    Real max_abs_residual;       // 最大绝对残差
    std::vector<Real> residuals;
    bool reject_good_fit;        // p < alpha → 拒绝好拟合
};

// =============================================================================
// iv_fit_goodness_test - IV 微笑拟合优度检验
//
// 原理: 模型 IV vs 市场 IV 的加权残差应无系统性偏差
//   χ² = Σ w_i * (IV_model_i - IV_market_i)²
//   其中 w_i = 1 / (spread_i / 2)²  (排幻觉点 H22)
//
// 排幻觉点 H22: 权重用市场 IV 的 Bid-Ask 宽度, 非简单方差
//   Bid-Ask 宽度大 → 权重小 (该点定价不确定性高, 容忍更大偏差)
//   Bid-Ask 宽度小 → 权重大 (该点定价精确, 偏差应小)
//
// H0: 模型 IV 与市场 IV 无显著差异 (好拟合)
// H1: 模型 IV 与市场 IV 有显著差异 (差拟合)
// χ² ~ χ²(N-1) under H0 (N-1 自由度, 因为 IV 曲面通常有 N-1 个独立参数)
//
// @param strikes 行权价
// @param maturities 到期日
// @param iv_market 市场 IV
// @param iv_model 模型 IV
// @param iv_bid_ask_spread IV 的 Bid-Ask 宽度 (权重 = 1/(spread/2)²)
// =============================================================================
inline IVFitGoodnessResult iv_fit_goodness_test(
    const std::vector<Real>& strikes,
    const std::vector<Real>& maturities,
    const std::vector<Real>& iv_market,
    const std::vector<Real>& iv_model,
    const std::vector<Real>& iv_bid_ask_spread,
    Real alpha = 0.05) {

    const Size n = iv_market.size();
    if (n < 3) {
        throw std::invalid_argument("iv_fit_goodness_test: need at least 3 points");
    }
    if (iv_model.size() != n || iv_bid_ask_spread.size() != n) {
        throw std::invalid_argument("iv_fit_goodness_test: size mismatch");
    }
    if (strikes.size() != n || maturities.size() != n) {
        throw std::invalid_argument("iv_fit_goodness_test: strikes/maturities size mismatch");
    }

    IVFitGoodnessResult result;
    result.residuals.resize(n);
    result.chi_squared = 0.0;
    result.max_abs_residual = 0.0;
    Real sum_sq_residual = 0.0;

    for (Size i = 0; i < n; ++i) {
        // 检查 spread 有效
        if (iv_bid_ask_spread[i] <= 0.0) {
            throw std::invalid_argument(
                "iv_fit_goodness_test: bid-ask spread must be positive");
        }

        result.residuals[i] = iv_model[i] - iv_market[i];
        const Real abs_res = std::fabs(result.residuals[i]);
        if (abs_res > result.max_abs_residual) {
            result.max_abs_residual = abs_res;
        }

        // 权重: 1 / (spread/2)² (排幻觉点 H22)
        // spread/2 是半宽度, 代表市场定价的标准差近似
        const Real half_spread = iv_bid_ask_spread[i] / 2.0;
        const Real weight = 1.0 / (half_spread * half_spread);

        result.chi_squared += weight * result.residuals[i] * result.residuals[i];
        sum_sq_residual += result.residuals[i] * result.residuals[i];
    }

    // 自由度: N - 1 (IV 曲面拟合通常有 N-1 个独立参数)
    result.degrees_of_freedom = n - 1;

    // p 值: P(χ² > observed) = chi2_sf(df, chi_squared)
    result.p_value = econometrics::detail::chi2_sf(
        static_cast<Real>(result.degrees_of_freedom), result.chi_squared);

    // RMSE
    result.rmse = std::sqrt(sum_sq_residual / static_cast<Real>(n));

    // 拒绝好拟合: p < alpha
    result.reject_good_fit = (result.p_value < alpha);

    return result;
}

// =============================================================================
// PriceResidualDiagnostics - 价格残差诊断结果
// =============================================================================
struct PriceResidualDiagnostics {
    Real mean_residual;
    Real std_residual;
    Real t_stat_bias;           // H0: mean = 0
    Real p_value_bias;          // 双尾 p 值
    bool has_bias;              // p < alpha → 有偏差
    std::vector<Real> residuals;
};

// =============================================================================
// price_residual_analysis - 模型价格 vs 市场价格残差诊断
//
// 原理: 检验残差无偏 (mean=0), 无系统性偏差
//   t = mean / (std / sqrt(N)) ~ t(N-1) (大样本用正态近似)
//
// @param market_prices 市场价格
// @param model_prices 模型价格
// =============================================================================
inline PriceResidualDiagnostics price_residual_analysis(
    const std::vector<Real>& market_prices,
    const std::vector<Real>& model_prices,
    Real alpha = 0.05) {

    const Size n = market_prices.size();
    if (n < 3) {
        throw std::invalid_argument("price_residual_analysis: need at least 3 prices");
    }
    if (model_prices.size() != n) {
        throw std::invalid_argument("price_residual_analysis: size mismatch");
    }

    PriceResidualDiagnostics result;
    result.residuals.resize(n);

    // 残差 = 模型价格 - 市场价格 (正值=模型高估, 负值=模型低估)
    Real sum = 0.0;
    for (Size i = 0; i < n; ++i) {
        result.residuals[i] = model_prices[i] - market_prices[i];
        sum += result.residuals[i];
    }

    // 均值
    result.mean_residual = sum / static_cast<Real>(n);

    // 标准差 (样本标准差, Bessel 校正)
    Real sum_sq = 0.0;
    for (Size i = 0; i < n; ++i) {
        const Real d = result.residuals[i] - result.mean_residual;
        sum_sq += d * d;
    }
    result.std_residual = std::sqrt(sum_sq / static_cast<Real>(n - 1));

    // t 统计量: t = mean / (std / sqrt(N))
    if (result.std_residual < 1e-300) {
        result.t_stat_bias = 0.0;
        result.p_value_bias = 1.0;
    } else {
        result.t_stat_bias = result.mean_residual /
            (result.std_residual / std::sqrt(static_cast<Real>(n)));
        // 大样本正态近似: p = 2 * Phi(-|t|)
        result.p_value_bias = 2.0 * std::erfc(std::fabs(result.t_stat_bias) /
                                               std::sqrt(2.0));
        if (result.p_value_bias > 1.0) result.p_value_bias = 1.0;
    }

    result.has_bias = (result.p_value_bias < alpha);

    return result;
}

}  // namespace v1
}  // namespace cpphub
