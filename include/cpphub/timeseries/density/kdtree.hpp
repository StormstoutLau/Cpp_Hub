// =============================================================================
// kdtree.hpp - 固定维 Chebyshev (L∞) k-d tree (Phase 9 P2-4, D9-3 延后项落地)
//
// 用途: 加速 KSG MI/CMI 与 TE 的联合空间 k 近邻与边际空间半径计数.
//       所有查询均为精确结果 (整数计数 / 顺序统计量), 与暴力参考逐位一致 —
//       树只改变访问路径, 不改变任何距离或计数的数值语义.
//
// 语义 (与 ksg.hpp / transfer_entropy.hpp 的暴力口径严格对齐):
//   - 点在行主序缓冲区 points[n × stride] 中; 仅 dims 指定的列参与距离,
//     距离 = Chebyshev (各维绝对差的最大值, 从 0.0 起逐维 max)
//   - kth_distance_including_self(q, k): 含查询点自身 (距离 0) 的第 k 小
//     Chebyshev 距离; k = 1 恒为 0 (自身). KSG 用 k+1 复现"排除自身"语义
//     (多重集去掉自身的一个 0 后, 第 k 小 = 含自身的第 k+1 小)
//   - count_within(q, r): 满足 cheb(q, j) <= r 的点数 (含自身), 精确计数;
//     r = 0 时统计与 q 完全重合的点 (含自身), 与 KSG 的 nextafter(ε,0) 口径
//     一致 (重复点场景不引入近似)
//   - 构建: 每层取展宽最大维的中位数切分, O(n log n · d) / O(n · d) 空间
//   - 查询: 典型 O(log n + k) 与 O(log n + 命中数); 高维/全重合数据退化 O(n)
//
// 生命周期: 非拥有视图 (仅存点索引与包围盒), points 缓冲区须在树对象存活
//           期内保持有效 (KSG/TE 在单次估计调用内部构造并立即使用).
// =============================================================================

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace density {

/// 固定维 Chebyshev k-d tree (header-only, 非拥有 points 视图)
class KdTree {
public:
    KdTree() = default;

    /// @throws std::invalid_argument n == 0 / dims 为空 / 列号越界 /
    ///         points 长度不足 / 选中列含非有限值
    KdTree(std::span<const Real> points, Size n, Size stride,
           std::span<const Size> dims) {
        build(points, n, stride, dims);
    }

    void build(std::span<const Real> points, Size n, Size stride,
               std::span<const Size> dims) {
        if (n == 0) {
            throw std::invalid_argument("KDTREE: n must be >= 1");
        }
        if (dims.empty()) {
            throw std::invalid_argument("KDTREE: dims must not be empty");
        }
        if (stride == 0 || n > points.size() / stride) {
            throw std::invalid_argument("KDTREE: points buffer is too small");
        }
        for (Size c : dims) {
            if (c >= stride) {
                throw std::invalid_argument("KDTREE: dims index out of range");
            }
        }
        for (Size i = 0; i < n; ++i) {
            for (Size c : dims) {
                if (!std::isfinite(points[i * stride + c])) {
                    throw std::invalid_argument(
                        "KDTREE: points contain non-finite value");
                }
            }
        }

        points_ = points;
        n_ = n;
        stride_ = stride;
        dims_.assign(dims.begin(), dims.end());
        idx_.resize(n);
        std::iota(idx_.begin(), idx_.end(), Size{0});
        nodes_.clear();
        nodes_.reserve(n);
        lo_.assign(n * dims_.size(), 0.0);
        hi_.assign(n * dims_.size(), 0.0);
        root_ = build_node(0, n);
    }

    /// @throws std::logic_error 未构建 (默认构造)
    Size size() const { return n_; }

    /// @throws std::logic_error 未构建 (默认构造)
    Size dims() const { return dims_.size(); }

    /// 含自身的第 k 小 Chebyshev 距离 (k ∈ [1, n])
    /// @throws std::invalid_argument 查询点越界 / k 越界
    Real kth_distance_including_self(Size query, Size k) const {
        require_ready();
        if (query >= n_) {
            throw std::invalid_argument("KDTREE: query index out of range");
        }
        if (k < 1 || k > n_) {
            throw std::invalid_argument("KDTREE: k must be in [1, n]");
        }
        std::vector<Real> buf(k, std::numeric_limits<Real>::infinity());
        knn_search(root_, query, buf);
        return buf[k - 1];
    }

    /// 半径内点数 (含自身): cheb(query, j) <= radius
    /// @throws std::invalid_argument 查询点越界 / radius 为负或 NaN
    Size count_within(Size query, Real radius) const {
        require_ready();
        if (query >= n_) {
            throw std::invalid_argument("KDTREE: query index out of range");
        }
        if (!(radius >= 0.0)) {
            throw std::invalid_argument("KDTREE: radius must be >= 0");
        }
        return count_node(root_, query, radius);
    }

private:
    static constexpr Size kNil = static_cast<Size>(-1);

    struct Node {
        Size point = 0;   ///< 本节点存储的点索引
        Size axis = 0;    ///< 切分轴 (dims_ 下标)
        Size left = kNil;  ///< 左子树 (< 切分值)
        Size right = kNil; ///< 右子树 (> 切分值)
        Size count = 0;   ///< 子树点数 (含本节点)
        Size box = 0;     ///< 包围盒偏移 (每节点 dims_.size())
    };

    void require_ready() const {
        if (n_ == 0) {
            throw std::logic_error("KDTREE: tree is not built");
        }
    }

    Size build_node(Size lo, Size hi) {
        const Size id = nodes_.size();
        nodes_.push_back(Node{});
        const Size d = dims_.size();
        const Size box = id * d;

        for (Size a = 0; a < d; ++a) {
            Real l = std::numeric_limits<Real>::infinity();
            Real h = -std::numeric_limits<Real>::infinity();
            for (Size t = lo; t < hi; ++t) {
                const Real v = points_[idx_[t] * stride_ + dims_[a]];
                l = std::min(l, v);
                h = std::max(h, v);
            }
            lo_[box + a] = l;
            hi_[box + a] = h;
        }

        Size axis = 0;
        Real best = -1.0;
        for (Size a = 0; a < d; ++a) {
            const Real spread = hi_[box + a] - lo_[box + a];
            if (spread > best) {
                best = spread;
                axis = a;
            }
        }

        const Size mid = lo + (hi - lo) / 2;
        const Size dim = dims_[axis];
        std::nth_element(idx_.begin() + static_cast<std::ptrdiff_t>(lo),
                         idx_.begin() + static_cast<std::ptrdiff_t>(mid),
                         idx_.begin() + static_cast<std::ptrdiff_t>(hi),
                         [this, dim](Size a, Size b) {
                             return points_[a * stride_ + dim] <
                                    points_[b * stride_ + dim];
                         });

        nodes_[id].point = idx_[mid];
        nodes_[id].axis = axis;
        nodes_[id].count = hi - lo;
        nodes_[id].box = box;
        if (mid > lo) {
            nodes_[id].left = build_node(lo, mid);
        }
        if (mid + 1 < hi) {
            nodes_[id].right = build_node(mid + 1, hi);
        }
        return id;
    }

    /// 查询点到点 b 的 Chebyshev 距离 (逐维 max, 与暴力口径逐位一致)
    Real cheb_distance(Size query, Size b) const {
        Real d = 0.0;
        for (Size a = 0; a < dims_.size(); ++a) {
            d = std::max(d, std::fabs(points_[query * stride_ + dims_[a]] -
                                      points_[b * stride_ + dims_[a]]));
        }
        return d;
    }

    /// 查询点到包围盒的最小 Chebyshev 距离 (0 = 盒内)
    Real box_min_distance(const Node& node, Size query) const {
        Real d = 0.0;
        for (Size a = 0; a < dims_.size(); ++a) {
            const Real v = points_[query * stride_ + dims_[a]];
            const Real l = lo_[node.box + a];
            const Real h = hi_[node.box + a];
            if (v < l) {
                d = std::max(d, l - v);
            } else if (v > h) {
                d = std::max(d, v - h);
            }
        }
        return d;
    }

    /// 查询点到包围盒的最大 Chebyshev 距离 (盒内所有点的上界)
    Real box_max_distance(const Node& node, Size query) const {
        Real d = 0.0;
        for (Size a = 0; a < dims_.size(); ++a) {
            const Real v = points_[query * stride_ + dims_[a]];
            d = std::max(d, std::max(std::fabs(v - lo_[node.box + a]),
                                     std::fabs(v - hi_[node.box + a])));
        }
        return d;
    }

    /// 定长有序缓冲插入 (与暴力实现的严格 < 语义一致)
    static void insert_sorted(std::vector<Real>& buf, Real d) {
        if (!(d < buf.back())) {
            return;
        }
        Size pos = buf.size() - 1;
        while (pos > 0 && buf[pos - 1] > d) {
            buf[pos] = buf[pos - 1];
            --pos;
        }
        buf[pos] = d;
    }

    void knn_search(Size node, Size query, std::vector<Real>& buf) const {
        const Node& nd = nodes_[node];
        if (buf.back() < box_min_distance(nd, query)) {
            return;
        }
        const Real d = cheb_distance(query, nd.point);
        insert_sorted(buf, d);

        Size first = nd.left;
        Size second = nd.right;
        if (nd.left != kNil && nd.right != kNil &&
            points_[query * stride_ + dims_[nd.axis]] >
                points_[nd.point * stride_ + dims_[nd.axis]]) {
            first = nd.right;
            second = nd.left;
        }
        if (first != kNil) {
            knn_search(first, query, buf);
        }
        if (second != kNil) {
            knn_search(second, query, buf);
        }
    }

    Size count_node(Size node, Size query, Real radius) const {
        const Node& nd = nodes_[node];
        if (box_max_distance(nd, query) <= radius) {
            return nd.count;
        }
        if (box_min_distance(nd, query) > radius) {
            return 0;
        }
        Size c = (cheb_distance(query, nd.point) <= radius) ? 1 : 0;
        if (nd.left != kNil) {
            c += count_node(nd.left, query, radius);
        }
        if (nd.right != kNil) {
            c += count_node(nd.right, query, radius);
        }
        return c;
    }

    std::span<const Real> points_;  ///< 非拥有, 行主序 n × stride
    Size n_ = 0;
    Size stride_ = 0;
    std::vector<Size> dims_;        ///< 参与距离的列
    std::vector<Size> idx_;         ///< 节点切分顺序 (build 期重排)
    std::vector<Node> nodes_;
    std::vector<Real> lo_;          ///< 包围盒下界 (每节点 dims_.size())
    std::vector<Real> hi_;          ///< 包围盒上界
    Size root_ = kNil;
};

}  // namespace density
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
