// =============================================================================
// greeks_consistency.hpp - Greeks 跨方法一致性检验 (P2)
//
// Phase 7A Wave 3d: Analytic vs Numerical vs AAD vs Pathwise vs LR
//
// 教材锚点:
//   - Glasserman 2003 §7 (Monte Carlo Greeks)
//   - Broadie-Glasserman 1996
//
// 依赖: 已有 risk/greeks/greeks_factory.hpp (Analytic/Numerical/AAD/Pathwise/LR)
//
// 排幻觉点:
//   H20: Pathwise/LR 有随机性, 用置信区间比较而非点估计
//        |mean_pathwise - analytic| < z_{0.975} * SE_pathwise
//        仅当 |diff| > z*SE 且 rel_diff > tolerance 时才判不一致
//   H21: Gamma 的 Numerical 用二阶差分, bump size 需更小 (1e-4 而非 1e-2)
//        二阶差分截断误差 O(dS²), dS=0.01 → 误差 ~1e-4, dS=1e-4 → 误差 ~1e-8
//
// 方法适用性矩阵:
//   | Method    | delta | gamma | vega | theta | rho |
//   |-----------|-------|-------|------|-------|-----|
//   | Analytic  |   Y   |   Y   |  Y   |   Y   |  Y  |
//   | Numerical |   Y   |   Y*  |  Y   |   Y   |  Y  |
//   | AAD       |   Y   |   Y   |  Y   |   Y   |  Y  |
//   | Pathwise  |   Y   |   N   |  Y   |   N   |  N  |
//   | LR        |   Y   |   N   |  Y   |   N   |  N  |
//   * Gamma 用 dS=1e-4 (H21), 其余用默认 dS=0.01
//
// 注: Pathwise/LR 不支持 gamma/theta/rho (gamma 需二阶 score O(1/T) 方差,
//     theta/rho 需路径对 T/r 的导数, Pathwise/LR 仅实现 S/σ)
// =============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/risk/greeks/greeks_factory.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"
#include "cpphub/risk/greeks/greeks_numerical.hpp"

namespace cpphub {
inline namespace v1 {

// =============================================================================
// GreeksConsistencyResult - Greeks 跨方法一致性检验结果
// =============================================================================
struct GreeksConsistencyResult {
    Real analytic_value = 0.0;       // Analytic (闭式解, 精确)
    Real numerical_value = 0.0;      // Numerical (FD, 中心差分)
    Real aad_value = 0.0;            // AAD (自动微分, 解析导数)
    Real pathwise_mean = 0.0;        // Pathwise (MC, 批次均值)
    Real pathwise_std_error = 0.0;   // Pathwise 标准误差 (批次均值法)
    Real lr_mean = 0.0;              // LR (MC, 批次均值)
    Real lr_std_error = 0.0;         // LR 标准误差
    Real max_discrepancy = 0.0;      // 最大相对偏差
    bool consistent = true;          // 所有方法在容差内一致
    std::vector<std::string> warnings;  // 不一致的方法对或不适用的方法
};

// =============================================================================
// greeks_consistency_check - Greeks 跨方法一致性检验
//
// 原理: 同一 Greeks 用不同方法计算应数值一致 (容差内)
//   - Analytic vs AAD: 应精确一致 (AAD 是解析导数, 容差 1e-10)
//   - Analytic vs Numerical: FD 精度受 bump size 影响 (容差 relative_tolerance)
//   - Analytic vs Pathwise: H20 用置信区间比较 (|diff| < z*SE 或 rel_diff < tol)
//   - Analytic vs LR: H20 用置信区间比较
//
// 批次均值法 (batch means):
//   将 n_paths 分成 n_batches=10 批, 每批 paths_per_batch = n_paths/10 条路径
//   每批独立计算 Greeks 估计, mean = 批均值, SE = 批标准差 / sqrt(n_batches)
//   这是 MC 标准误差的稳健估计 (排幻觉点 H18 同理)
//
// @param S 标的现价
// @param K 行权价
// @param T 到期时间 (年)
// @param r 无风险利率
// @param q 股息率
// @param sigma 波动率
// @param payoff 期权类型 (VanillaCall/VanillaPut/DigitalCall/DigitalPut)
// @param greek_name Greeks 名称 ("delta"/"gamma"/"vega"/"theta"/"rho")
// @param n_paths MC 路径总数 (分 10 批)
// @param seed 随机种子
// @param relative_tolerance 相对容差 (默认 1%)
// =============================================================================
inline GreeksConsistencyResult greeks_consistency_check(
    Real S, Real K, Real T, Real r, Real q, Real sigma,
    PayoffType payoff,
    const std::string& greek_name,
    Size n_paths = 100000,
    uint64_t seed = 42,
    Real relative_tolerance = 0.01) {

    // 参数校验
    if (greek_name != "delta" && greek_name != "gamma" &&
        greek_name != "vega" && greek_name != "theta" && greek_name != "rho") {
        throw std::invalid_argument(
            "greeks_consistency_check: unknown greek '" + greek_name +
            "' (expected: delta/gamma/vega/theta/rho)");
    }
    if (S <= 0.0 || K <= 0.0 || T <= 0.0 || sigma <= 0.0) {
        throw std::invalid_argument(
            "greeks_consistency_check: S, K, T, sigma must be positive");
    }
    if (n_paths < 100) {
        throw std::invalid_argument(
            "greeks_consistency_check: n_paths must be >= 100");
    }

    GreeksConsistencyResult result;

    // 辅助: 从 UnifiedGreeks 提取指定 Greeks 值
    auto extract = [](const UnifiedGreeks& g, const std::string& name) -> Real {
        if (name == "delta") return g.delta;
        if (name == "gamma") return g.gamma;
        if (name == "vega")  return g.vega;
        if (name == "theta") return g.theta;
        return g.rho;  // "rho"
    };

    const bool is_vanilla =
        (payoff == PayoffType::VanillaCall || payoff == PayoffType::VanillaPut);
    const bool is_call =
        (payoff == PayoffType::VanillaCall || payoff == PayoffType::DigitalCall);

    // --- 1. Analytic (闭式解, 精确) ---
    // 注: digital 时 factory 回退到 LR, analytic_value 仍是有效估计
    auto g_analytic = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, payoff, GreeksMethod::Analytic);
    result.analytic_value = extract(g_analytic, greek_name);

    // --- 2. AAD (自动微分, 解析导数) ---
    auto g_aad = GreeksFactory::compute_bsm(
        S, K, T, r, q, sigma, payoff, GreeksMethod::AAD);
    result.aad_value = extract(g_aad, greek_name);

    // --- 3. Numerical (FD, 中心差分) ---
    // H21: Gamma 用更小的 bump size (dS=1e-4, 二阶差分截断误差 O(dS²))
    if (greek_name == "gamma" && is_vanilla) {
        // 直接调用 NumericalGreeksEngine, 用 dS=1e-4 (排幻觉点 H21)
        NumericalGreeksEngine::PriceFn price_fn =
            [](Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call) -> Real {
                return AnalyticGreeksEngine::bsm_european(
                    S, K, T, r, q, sigma, is_call).price;
            };
        auto g_fd = NumericalGreeksEngine::bsm_european(
            S, K, T, r, q, sigma, is_call, price_fn,
            1e-4,       // dS — H21: 二阶差分需更小 bump
            0.0001,     // dSigma
            1.0 / 365.0, // dT
            0.0001);    // dR
        result.numerical_value = g_fd.gamma;
    } else {
        auto g_fd = GreeksFactory::compute_bsm(
            S, K, T, r, q, sigma, payoff, GreeksMethod::FD);
        result.numerical_value = extract(g_fd, greek_name);
    }

    // --- 4. Pathwise (MC, 批次均值法估计 SE) ---
    // Pathwise 仅支持 delta/vega, 且仅 vanilla (光滑 payoff)
    bool pathwise_applicable =
        is_vanilla && (greek_name == "delta" || greek_name == "vega");

    if (pathwise_applicable) {
        const Size n_batches = 10;
        const Size paths_per_batch =
            std::max(n_paths / n_batches, static_cast<Size>(1000));
        std::vector<Real> batch_estimates(n_batches);
        for (Size b = 0; b < n_batches; ++b) {
            auto g = GreeksFactory::compute_bsm(
                S, K, T, r, q, sigma, payoff,
                GreeksMethod::Pathwise,
                paths_per_batch, seed + b);
            batch_estimates[b] = extract(g, greek_name);
        }
        // 批次均值
        Real sum = 0.0;
        for (Real v : batch_estimates) sum += v;
        result.pathwise_mean = sum / static_cast<Real>(n_batches);
        // SE = sample_std / sqrt(n_batches)
        Real sum_sq = 0.0;
        for (Real v : batch_estimates) {
            Real d = v - result.pathwise_mean;
            sum_sq += d * d;
        }
        Real sample_var = sum_sq / static_cast<Real>(n_batches - 1);
        result.pathwise_std_error =
            std::sqrt(sample_var / static_cast<Real>(n_batches));
    } else {
        result.pathwise_mean = 0.0;
        result.pathwise_std_error = 0.0;
        if (!is_vanilla) {
            result.warnings.push_back(
                "Pathwise not applicable: discontinuous payoff (digital)");
        } else {
            result.warnings.push_back(
                "Pathwise does not support " + greek_name +
                " (only delta/vega for smooth payoff)");
        }
    }

    // --- 5. LR (MC, 批次均值法估计 SE) ---
    // LR 支持 delta/vega (vanilla 和 digital 均可)
    bool lr_applicable = (greek_name == "delta" || greek_name == "vega");

    if (lr_applicable) {
        const Size n_batches = 10;
        const Size paths_per_batch =
            std::max(n_paths / n_batches, static_cast<Size>(1000));
        std::vector<Real> batch_estimates(n_batches);
        for (Size b = 0; b < n_batches; ++b) {
            auto g = GreeksFactory::compute_bsm(
                S, K, T, r, q, sigma, payoff,
                GreeksMethod::LR,
                paths_per_batch, seed + 1000 + b);
            batch_estimates[b] = extract(g, greek_name);
        }
        Real sum = 0.0;
        for (Real v : batch_estimates) sum += v;
        result.lr_mean = sum / static_cast<Real>(n_batches);
        Real sum_sq = 0.0;
        for (Real v : batch_estimates) {
            Real d = v - result.lr_mean;
            sum_sq += d * d;
        }
        Real sample_var = sum_sq / static_cast<Real>(n_batches - 1);
        result.lr_std_error =
            std::sqrt(sample_var / static_cast<Real>(n_batches));
    } else {
        result.lr_mean = 0.0;
        result.lr_std_error = 0.0;
        result.warnings.push_back(
            "LR does not support " + greek_name +
            " (only delta/vega; gamma needs 2nd-order score)");
    }

    // --- 一致性检验 ---
    result.max_discrepancy = 0.0;
    result.consistent = true;

    // z_{0.975} = 1.959963984540054 (标准正态 97.5% 分位数)
    constexpr Real z_975 = 1.959963984540054;

    auto relative_diff = [](Real a, Real b) -> Real {
        Real ref = std::max(std::abs(a), std::abs(b));
        if (ref < 1e-15) return 0.0;
        return std::abs(a - b) / ref;
    };

    auto update_max = [&](Real rd) {
        if (rd > result.max_discrepancy) result.max_discrepancy = rd;
    };

    // (a) Analytic vs Numerical (确定性方法, 相对容差)
    {
        Real rd = relative_diff(result.analytic_value, result.numerical_value);
        update_max(rd);
        if (rd > relative_tolerance) {
            result.consistent = false;
            result.warnings.push_back(
                "Analytic vs Numerical: rel_diff=" + std::to_string(rd) +
                " > tol=" + std::to_string(relative_tolerance));
        }
    }

    // (b) Analytic vs AAD (应精确一致, AAD 是解析导数)
    //     容差用 max(relative_tolerance * 0.01, 1e-10) — 远比整体容差严格
    {
        Real rd = relative_diff(result.analytic_value, result.aad_value);
        update_max(rd);
        Real aad_tol = std::max(relative_tolerance * 0.01,
                                static_cast<Real>(1e-10));
        if (rd > aad_tol) {
            result.consistent = false;
            result.warnings.push_back(
                "Analytic vs AAD: rel_diff=" + std::to_string(rd) +
                " > tol=" + std::to_string(aad_tol) +
                " (should be machine precision)");
        }
    }

    // (c) Analytic vs Pathwise (随机方法, H20: 置信区间比较)
    if (pathwise_applicable) {
        Real diff = std::abs(result.analytic_value - result.pathwise_mean);
        Real ci = z_975 * result.pathwise_std_error;
        Real rd = relative_diff(result.analytic_value, result.pathwise_mean);
        update_max(rd);
        // H20: 仅当 |diff| > z*SE 且 rel_diff > tolerance 时判不一致
        //      |diff| <= z*SE → 差异在采样噪声范围内, 一致
        if (diff > ci && rd > relative_tolerance) {
            result.consistent = false;
            result.warnings.push_back(
                "Analytic vs Pathwise: |diff|=" + std::to_string(diff) +
                " > z*SE=" + std::to_string(ci) +
                ", rel_diff=" + std::to_string(rd));
        }
    }

    // (d) Analytic vs LR (随机方法, H20: 置信区间比较)
    if (lr_applicable) {
        Real diff = std::abs(result.analytic_value - result.lr_mean);
        Real ci = z_975 * result.lr_std_error;
        Real rd = relative_diff(result.analytic_value, result.lr_mean);
        update_max(rd);
        if (diff > ci && rd > relative_tolerance) {
            result.consistent = false;
            result.warnings.push_back(
                "Analytic vs LR: |diff|=" + std::to_string(diff) +
                " > z*SE=" + std::to_string(ci) +
                ", rel_diff=" + std::to_string(rd));
        }
    }

    return result;
}

}  // namespace v1
}  // namespace cpphub
