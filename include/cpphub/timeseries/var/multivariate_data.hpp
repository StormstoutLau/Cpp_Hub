// =============================================================================
// multivariate_data.hpp - 多变量时序数据载体 MultivariateTSData (spec §4.1, C4)
//
// Phase 7C v1.7 M2 (PHASE7C_SPEC.md v1.2 §4.1)
//
// 契约 (spec §4.1 + §1.4-5):
//   - K 列等长 T, 不等长/NaN/K<1/T<2 ⇒ throw std::invalid_argument
//   - matrix(): T×K 行=时间 (v1.2)
//   - reorder(order): 变量重排副本 — Cholesky 排序敏感性检验用 (V2)
//   - names 可空 (DY 输出用)
//
// Eigen 边界 (ADR-013/§8.2): 本头文件属 cpphub_timeseries_mat 载体,
//   公共接口使用 Eigen::MatrixXd; 单变量 timeseries 头不得 include 本文件
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace var {

/// 多变量时序数据 (C4): K 列等长
struct MultivariateTSData {
    std::vector<std::vector<Real>> columns;  ///< K 列, 每列 T
    std::vector<std::string> names;          ///< 变量名 (可空, DY 输出/重排用)

    /// 有效样本数 (等长校验, 不等长 ⇒ invalid_argument)
    Size T() const {
        validate();
        return columns.empty() ? 0 : columns[0].size();
    }

    /// 变量数 (≥1)
    Size K() const {
        validate();
        return columns.size();
    }

    /// T×K 矩阵, 行=时间 行=变量序 (v1.2)
    Eigen::MatrixXd matrix() const {
        validate();
        const Size t = T(), k = K();
        Eigen::MatrixXd m(t, k);
        for (Size j = 0; j < k; ++j) {
            for (Size i = 0; i < t; ++i) m(i, j) = columns[j][i];
        }
        return m;
    }

    /// 变量重排副本: order[j] = 原变量索引 → 新位置 j (排序敏感性检验, V2)
    MultivariateTSData reorder(const std::vector<Size>& order) const {
        validate();
        if (order.size() != columns.size()) {
            throw std::invalid_argument("reorder: order size mismatch with K");
        }
        std::vector<bool> seen(columns.size(), false);
        for (Size idx : order) {
            if (idx >= columns.size() || seen[idx]) {
                throw std::invalid_argument("reorder: order must be a permutation of 0..K-1");
            }
            seen[idx] = true;
        }
        MultivariateTSData out;
        out.columns.resize(columns.size());
        for (Size j = 0; j < order.size(); ++j) out.columns[j] = columns[order[j]];
        if (!names.empty() && names.size() == columns.size()) {
            out.names.resize(order.size());
            for (Size j = 0; j < order.size(); ++j) out.names[j] = names[order[j]];
        }
        return out;
    }

    /// 校验: K≥1, T≥2, 等长, 无 NaN
    void validate() const {
        if (columns.empty()) throw std::invalid_argument("data: no columns (K=0)");
        const Size t = columns[0].size();
        if (t < 2) throw std::invalid_argument("data: T must be >= 2");
        for (const auto& col : columns) {
            if (col.size() != t) {
                throw std::invalid_argument("data: unequal column lengths");
            }
            for (Real v : col) {
                if (std::isnan(v)) throw std::invalid_argument("data: NaN in column");
            }
        }
        if (!names.empty() && names.size() != columns.size()) {
            throw std::invalid_argument("data: names size mismatch with K");
        }
    }
};

}  // namespace var
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
