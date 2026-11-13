// =============================================================================
// knn_density.hpp - 1D 暴力 kNN 密度 (Phase 9 M2, PHASE9_SPEC v1.1 §4.1)
//
// Loftsgaarden-Quesenberry:
//   pdf(x) = k / (n · 2 · r_k(x)),  r_k = 到第 k 近邻的绝对距离
//   1D 下 Chebyshev 距离即 |x_j − x_i|.
//
// 约定:
//   - 近邻计数包含查询点自身 (dist=0 占第 1 位, 与 sklearn kneighbors
//     排除自身的语义不同; 黄金值按含自身口径冻结)
//   - 构造时复制并升序排序; pdf 用 lower_bound 定位后左右归并取第 k 小,
//     O(log n + k) 查询
//   - x 落在样本范围外时半径取到端点距离; 边界外密度不归一 (D9-8)
// =============================================================================

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace density {

/// kNN 密度配置
struct KnnConfig {
    std::size_t k = 5;
};

/// 一维 kNN 密度估计 (header-only)
class KnnDensity1D {
public:
    /// @throws std::invalid_argument 空样本 / 非有限值 / k < 1 / k >= n
    explicit KnnDensity1D(std::span<const Real> sample, KnnConfig config = {}) {
        if (sample.empty()) {
            throw std::invalid_argument("KNN: sample must not be empty");
        }
        for (Real v : sample) {
            if (!std::isfinite(v)) {
                throw std::invalid_argument("KNN: sample contains non-finite value");
            }
        }
        if (config.k < 1) {
            throw std::invalid_argument("KNN: k must be >= 1");
        }
        if (config.k >= sample.size()) {
            throw std::invalid_argument("KNN: k must be < sample size");
        }
        sorted_.assign(sample.begin(), sample.end());
        std::sort(sorted_.begin(), sorted_.end());
        k_ = config.k;
        n_ = static_cast<Real>(sample.size());
    }

    /// Loftsgaarden-Quesenberry 密度 k/(n·2r_k); 样本范围外不归一 (文档化)
    Real pdf(Real x) const {
        const Real r = kth_radius(x);
        return static_cast<Real>(k_) / (n_ * 2.0 * r);
    }

    /// 网格批量求值 (等价于逐点 pdf)
    std::vector<Real> pdf(std::span<const Real> grid) const {
        std::vector<Real> out;
        out.reserve(grid.size());
        for (Real x : grid) {
            out.push_back(pdf(x));
        }
        return out;
    }

    std::size_t k() const noexcept { return k_; }

private:
    /// 第 k 小绝对距离 (含自身): lower_bound 定位后左右归并 k 步
    Real kth_radius(Real x) const {
        const Size n = sorted_.size();
        const Size idx = static_cast<Size>(
            std::lower_bound(sorted_.begin(), sorted_.end(), x) - sorted_.begin());
        Size left = idx;   // 下一个左候选: left − 1
        Size right = idx;  // 下一个右候选: right
        constexpr Real inf = std::numeric_limits<Real>::infinity();
        Real r = 0.0;
        for (Size step = 0; step < k_; ++step) {
            const Real dl = (left > 0) ? (x - sorted_[left - 1]) : inf;
            const Real dr = (right < n) ? (sorted_[right] - x) : inf;
            if (dl <= dr) {
                r = dl;
                --left;
            } else {
                r = dr;
                ++right;
            }
        }
        return r;
    }

    std::vector<Real> sorted_;  ///< 升序样本
    std::size_t k_ = 1;
    Real n_ = 0.0;
};

}  // namespace density
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
