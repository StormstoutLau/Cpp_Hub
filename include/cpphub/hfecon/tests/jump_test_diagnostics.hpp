// =============================================================================
// jump_test_diagnostics.hpp - 跳跃检验多重修正与联合诊断 (P2)
//
// Phase 7A Wave 3: 高频跳跃检验的多重检验修正
//
// 教材锚点:
//   - Bonferroni 1936 (FWER 控制)
//   - Benjamini-Hochberg 1995 (FDR 控制)
//
// 原理:
//   多个跳跃检验 (BNS/AJ/JO/Rank) 联合推断时, 单个检验 α=0.05 会导致
//   总体第一类错误膨胀。需用多重检验修正控制 FWER 或 FDR。
//
// 排幻觉点:
//   H23: BH 修正控制 FDR (False Discovery Rate), 非 FWER
//     - Bonferroni: α' = α/m (保守, FWER ≤ α)
//     - BH: 排序后 p_(i) ≤ (i/m)·α (FDR ≤ α)
//
// ADR-015 方案 B: 仅依赖 core/, 不依赖 Eigen3
// =============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// MultipleTestCorrectionResult - 多重检验修正结果
// =============================================================================
struct MultipleTestCorrectionResult {
    enum class Method { Bonferroni, BenjaminiHochberg };

    std::vector<Real> original_p_values;
    std::vector<Real> adjusted_p_values;  // 修正后的 p 值
    std::vector<bool> reject_null;        // 修正后是否拒绝 H0
    Method method;
    Size n_tests;
    Size n_rejections;
};

// =============================================================================
// multiple_test_correction - 多重检验修正
//
// Bonferroni 修正:
//   adjusted_p = min(m * p_i, 1.0)
//   reject if adjusted_p < α
//   特点: 保守, 控制 FWER (Family-Wise Error Rate) ≤ α
//
// Benjamini-Hochberg 修正 (排幻觉点 H23):
//   1. 将 p 值升序排列: p_(1) ≤ p_(2) ≤ ... ≤ p_(m)
//   2. 找最大 k 使得 p_(k) ≤ (k/m) * α
//   3. 拒绝 H0_(1), ..., H0_(k)
//   4. adjusted_p_(i) = min_{j>=i} (m * p_(j) / j), 取累积最小值后截断到 [0,1]
//   特点: 控制 FDR (False Discovery Rate) ≤ α, 比 Bonferroni 更强力
//
// @param p_values 原始 p 值向量
// @param method 修正方法 (Bonferroni 或 BenjaminiHochberg)
// @param alpha 显著性水平 (默认 0.05)
// @return MultipleTestCorrectionResult
// =============================================================================
inline MultipleTestCorrectionResult multiple_test_correction(
    const std::vector<Real>& p_values,
    MultipleTestCorrectionResult::Method method,
    Real alpha = 0.05) {

    if (p_values.empty()) {
        throw std::invalid_argument("multiple_test_correction: p_values is empty");
    }
    if (alpha <= 0.0 || alpha >= 1.0) {
        throw std::invalid_argument("multiple_test_correction: alpha must be in (0, 1)");
    }

    const Size m = p_values.size();
    MultipleTestCorrectionResult result;
    result.method = method;
    result.n_tests = m;
    result.original_p_values = p_values;
    result.adjusted_p_values.resize(m);
    result.reject_null.resize(m);

    if (method == MultipleTestCorrectionResult::Method::Bonferroni) {
        // Bonferroni: adjusted_p = min(m * p_i, 1.0)
        for (Size i = 0; i < m; ++i) {
            result.adjusted_p_values[i] = std::min(static_cast<Real>(m) * p_values[i], 1.0);
            result.reject_null[i] = (result.adjusted_p_values[i] < alpha);
        }
    } else {
        // Benjamini-Hochberg (排幻觉点 H23: 控制 FDR, 非 FWER)
        // 1. 创建排序索引
        std::vector<Size> order(m);
        for (Size i = 0; i < m; ++i) order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](Size a, Size b) { return p_values[a] < p_values[b]; });

        // 2. 计算排序后的 adjusted p: p_adj_(k) = min_{j>=k} (m * p_(j) / j)
        //    从大到小 (k=m-1 到 0) 取累积最小值
        std::vector<Real> sorted_adjusted(m);
        Real cum_min = 1.0;  // 截断到 [0, 1]
        for (int k = static_cast<int>(m) - 1; k >= 0; --k) {
            const Size kk = static_cast<Size>(k);
            const Real rank = static_cast<Real>(kk + 1);  // 1-indexed rank
            const Real raw_adj = static_cast<Real>(m) * p_values[order[kk]] / rank;
            cum_min = std::min(cum_min, raw_adj);
            sorted_adjusted[kk] = std::min(cum_min, 1.0);  // 截断到 [0, 1]
        }

        // 3. 将排序后的 adjusted p 映射回原始顺序
        for (Size i = 0; i < m; ++i) {
            result.adjusted_p_values[order[i]] = sorted_adjusted[i];
        }

        // 4. BH 拒绝规则: 找最大 k 使 p_(k) <= (k/m) * alpha, 拒绝 1..k
        //    等价于 adjusted_p < alpha
        for (Size i = 0; i < m; ++i) {
            result.reject_null[i] = (result.adjusted_p_values[i] < alpha);
        }
    }

    // 统计拒绝数
    result.n_rejections = 0;
    for (Size i = 0; i < m; ++i) {
        if (result.reject_null[i]) ++result.n_rejections;
    }

    return result;
}

// =============================================================================
// JumpTestDiagnosticsResult - 跳跃检验联合诊断结果
// =============================================================================
struct JumpTestDiagnosticsResult {
    std::vector<Real> test_statistics;  // [BNS, AJ, JO, Rank] 每天的统计量
    std::vector<Real> p_values;         // 对应 p 值
    MultipleTestCorrectionResult bonferroni;
    MultipleTestCorrectionResult bh;
    Size consensus_jumps;  // 多数投票的跳跃天数 (>=3/4 检验拒绝)
    bool consistent;       // 四检验结论一致 (全拒绝或全不拒绝)
};

// =============================================================================
// jump_test_diagnostics - 跳跃检验联合诊断
//
// 对每天的跳跃检验结果进行多重修正和联合推断:
//   1. 收集四类跳跃检验 (BNS/AJ/JO/Rank) 的统计量和 p 值
//   2. 分别用 Bonferroni 和 BH 修正
//   3. 共识跳跃: >=3/4 检验拒绝 → 确认为跳跃
//   4. 一致性: 四检验全拒绝或全不拒绝 → consistent=true
//
// @param bns_stats BNS 检验统计量向量 (每天一个)
// @param aj_stats AJ 检验统计量向量
// @param jo_stats JO 检验统计量向量
// @param rank_stats Rank 检验统计量向量
// @param alpha 显著性水平
// @return JumpTestDiagnosticsResult
// =============================================================================
inline JumpTestDiagnosticsResult jump_test_diagnostics(
    const std::vector<Real>& bns_stats,
    const std::vector<Real>& aj_stats,
    const std::vector<Real>& jo_stats,
    const std::vector<Real>& rank_stats,
    Real alpha = 0.05) {

    const Size n_days = bns_stats.size();
    if (aj_stats.size() != n_days || jo_stats.size() != n_days ||
        rank_stats.size() != n_days) {
        throw std::invalid_argument("jump_test_diagnostics: size mismatch");
    }
    if (n_days == 0) {
        throw std::invalid_argument("jump_test_diagnostics: empty input");
    }

    JumpTestDiagnosticsResult result;

    // 收集统计量和 p 值 (双尾 Z 检验: p = 2 * Phi(-|Z|))
    result.test_statistics.resize(n_days);
    result.p_values.resize(n_days);
    for (Size d = 0; d < n_days; ++d) {
        // 取四检验统计量的平均值作为代表
        result.test_statistics[d] =
            (bns_stats[d] + aj_stats[d] + jo_stats[d] + rank_stats[d]) / 4.0;

        // 每天有 4 个检验, 收集 4 个 p 值
        std::vector<Real> daily_p_values(4);
        daily_p_values[0] = 2.0 * std::erfc(std::fabs(bns_stats[d]) / std::sqrt(2.0));
        daily_p_values[1] = 2.0 * std::erfc(std::fabs(aj_stats[d]) / std::sqrt(2.0));
        daily_p_values[2] = 2.0 * std::erfc(std::fabs(jo_stats[d]) / std::sqrt(2.0));
        daily_p_values[3] = 2.0 * std::erfc(std::fabs(rank_stats[d]) / std::sqrt(2.0));

        // 取最小 p 值作为当天的代表 (最显著的检验)
        result.p_values[d] = *std::min_element(daily_p_values.begin(),
                                                daily_p_values.end());
    }

    // 多重检验修正 (对 n_days * 4 个检验)
    // 展开 4 * n_days 个 p 值
    std::vector<Real> all_p_values;
    all_p_values.reserve(4 * n_days);
    for (Size d = 0; d < n_days; ++d) {
        all_p_values.push_back(2.0 * std::erfc(std::fabs(bns_stats[d]) / std::sqrt(2.0)));
        all_p_values.push_back(2.0 * std::erfc(std::fabs(aj_stats[d]) / std::sqrt(2.0)));
        all_p_values.push_back(2.0 * std::erfc(std::fabs(jo_stats[d]) / std::sqrt(2.0)));
        all_p_values.push_back(2.0 * std::erfc(std::fabs(rank_stats[d]) / std::sqrt(2.0)));
    }

    result.bonferroni = multiple_test_correction(
        all_p_values, MultipleTestCorrectionResult::Method::Bonferroni, alpha);
    result.bh = multiple_test_correction(
        all_p_values, MultipleTestCorrectionResult::Method::BenjaminiHochberg, alpha);

    // 共识跳跃: 每天检查 >=3/4 检验拒绝 (Bonferroni 修正后)
    result.consensus_jumps = 0;
    bool all_consistent = true;
    for (Size d = 0; d < n_days; ++d) {
        Size daily_rejections = 0;
        for (Size t = 0; t < 4; ++t) {
            if (result.bonferroni.reject_null[4 * d + t]) ++daily_rejections;
        }
        // 共识: >=3/4 检验拒绝
        if (daily_rejections >= 3) ++result.consensus_jumps;
        // 一致性: 全拒绝 (4/4) 或全不拒绝 (0/4)
        if (daily_rejections != 0 && daily_rejections != 4) {
            all_consistent = false;
        }
    }
    result.consistent = all_consistent;

    return result;
}

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
