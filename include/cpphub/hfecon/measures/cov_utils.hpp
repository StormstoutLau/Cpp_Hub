// =============================================================================
// cov_utils.hpp
// Phase 5 v1.4.2 - 多资产 Cov 估计共享工具
//
// 功能:
//   - make_psd: 对称矩阵 PSD 投影 (R highfrequency makePsd 对标)
//   - refresh_time_matching: 非同步时间匹配 (R refreshTimeMatching 对标)
//
// SOURCE: PHASE5_HFE_SPEC §5.4
//   R highfrequency 1.0.3 src/internals.cpp (refreshTimeMatching L23-72)
//   R highfrequency 1.0.3 R/makePsd.R (eigen 投影)
//
// R 对照:
//   makePsd(sigma)              -> make_psd(mat, n)
//   refreshTimeMatching(x, idx) -> refresh_time_matching(x, idx)
// =============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// make_psd: 对称矩阵 PSD 投影
//
// R highfrequency makePsd 实现:
//   eig <- eigen(sigma, symmetric = TRUE)
//   eig$values <- pmax(eig$values, 0)
//   sigma <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
//
// C++ 实现: Jacobi 特征值分解 (对称矩阵, 自包含, 无外部线性代数依赖)
// 输入: mat (n x n row-major 对称矩阵), n
// 输出: PSD 投影后的矩阵 (n x n row-major)
// =============================================================================
inline std::vector<Real> make_psd(const std::vector<Real>& mat, Size n) {
    if (n == 0) return {};
    if (n == 1) {
        return {std::max(mat[0], 0.0)};
    }

    // Jacobi 特征值分解: A = Q * diag(lambda) * Q^T
    std::vector<Real> A = mat;  // 工作副本, 会被旋转为对角阵
    std::vector<Real> Q(n * n, 0.0);
    for (Size i = 0; i < n; ++i) Q[i * n + i] = 1.0;

    const int max_iter = 100;
    const Real tol = 1e-15;

    for (int iter = 0; iter < max_iter; ++iter) {
        // 寻找最大非对角元素
        Size p = 0, q = 1;
        Real max_val = 0.0;
        for (Size i = 0; i < n; ++i) {
            for (Size j = i + 1; j < n; ++j) {
                Real val = std::fabs(A[i * n + j]);
                if (val > max_val) {
                    max_val = val;
                    p = i;
                    q = j;
                }
            }
        }
        if (max_val < tol) break;

        const Real app = A[p * n + p];
        const Real aqq = A[q * n + q];
        const Real apq = A[p * n + q];
        if (std::fabs(apq) < 1e-300) break;

        // Jacobi 旋转参数
        const Real tau = (aqq - app) / (2.0 * apq);
        const Real t = std::copysign(
            1.0 / (std::fabs(tau) + std::sqrt(tau * tau + 1.0)), tau);
        const Real c = 1.0 / std::sqrt(t * t + 1.0);
        const Real s = t * c;

        // 更新对角元素
        A[p * n + p] = app - t * apq;
        A[q * n + q] = aqq + t * apq;
        A[p * n + q] = 0.0;
        A[q * n + p] = 0.0;

        // 更新非对角元素 (i != p, q)
        for (Size i = 0; i < n; ++i) {
            if (i != p && i != q) {
                const Real aip = A[i * n + p];
                const Real aiq = A[i * n + q];
                A[i * n + p] = c * aip - s * aiq;
                A[p * n + i] = A[i * n + p];
                A[i * n + q] = s * aip + c * aiq;
                A[q * n + i] = A[i * n + q];
            }
        }

        // 更新特征向量矩阵 Q
        for (Size i = 0; i < n; ++i) {
            const Real qip = Q[i * n + p];
            const Real qiq = Q[i * n + q];
            Q[i * n + p] = c * qip - s * qiq;
            Q[i * n + q] = s * qip + c * qiq;
        }
    }

    // 重建: A_new = Q * diag(max(lambda, 0)) * Q^T
    std::vector<Real> result(n * n, 0.0);
    for (Size k = 0; k < n; ++k) {
        const Real lambda = std::max(A[k * n + k], 0.0);
        if (lambda == 0.0) continue;
        for (Size i = 0; i < n; ++i) {
            const Real qik = Q[i * n + k];
            if (qik == 0.0) continue;
            for (Size j = 0; j < n; ++j) {
                result[i * n + j] += qik * lambda * Q[j * n + k];
            }
        }
    }

    // 对称化 (消除数值不对称)
    for (Size i = 0; i < n; ++i) {
        for (Size j = i + 1; j < n; ++j) {
            const Real avg = (result[i * n + j] + result[j * n + i]) / 2.0;
            result[i * n + j] = avg;
            result[j * n + i] = avg;
        }
    }

    return result;
}

// =============================================================================
// refresh_time_matching: 非同步时间匹配 (R refreshTimeMatching 对标)
//
// R highfrequency 1.0.3 src/internals.cpp L23-72:
//   输入: x (N x D 矩阵, 可含 NA), idx (时间索引)
//   输出: 对齐后的矩阵 (使用 refresh time scheme)
//
// 算法:
//   1. 第一个 refresh time = 所有资产第一个观测时间的最大值
//   2. 后续 refresh time = 每个资产下一个观测时间的最小值
//   3. 在每个 refresh time, 使用各资产最近的观测值
//
// 输入: x (D 个资产, 每个长度 N, NaN 表示缺失), idx (N 个时间戳)
// 输出: pair (对齐矩阵 D x M, 对齐时间戳 M)
//
// 注意: 此函数为 rAVGCov/rTSCov/rMRCov 等多资产 Cov 估计器的预处理步骤.
//       rHYCov 不使用此函数 (它用 pcovcc 整数索引版).
// =============================================================================
inline std::pair<std::vector<std::vector<Real>>, std::vector<Real>>
refresh_time_matching(const std::vector<std::vector<Real>>& x,
                      const std::vector<Real>& idx) {
    const Size N = x.size();        // 资产数 (注意: 输入按资产分行)
    if (N == 0) return {{}, {}};
    const Size M = x[0].size();     // 观测数
    for (Size k = 0; k < N; ++k) {
        if (x[k].size() != M) {
            throw std::invalid_argument(
                "refresh_time_matching: all assets must have equal length");
        }
    }

    // 转置: 按时间行存储 (time x asset), 便于按行处理
    // x_t[a] = x[a][t]
    // NaN 表示该时刻该资产无观测

    // refresh time 方案:
    // 1. 找到每个资产的第一个有效观测时间
    // 2. 第一个 refresh time = max(各资产第一个观测时间)
    // 3. 后续 refresh time = min(各资产在当前 refresh time 之后的下一个观测时间)
    // 4. 在每个 refresh time, 取各资产最近的 (<= refresh time) 观测值

    // 简化实现: 假设 idx 已排序, x 中的 NaN 表示缺失
    // 找到每个资产的第一个有效索引
    std::vector<Size> first_valid(N, M);
    for (Size a = 0; a < N; ++a) {
        for (Size t = 0; t < M; ++t) {
            if (!std::isnan(x[a][t])) {
                first_valid[a] = t;
                break;
            }
        }
    }

    // 检查是否所有资产都有至少一个有效观测
    for (Size a = 0; a < N; ++a) {
        if (first_valid[a] >= M) {
            return {{}, {}};  // 某资产全 NaN, 返回空
        }
    }

    // refresh time 序列
    std::vector<Real> refresh_idx;
    std::vector<std::vector<Real>> refresh_data(N);

    // 第一个 refresh time = max(各资产第一个观测的时间戳)
    Real current_refresh = idx[first_valid[0]];
    for (Size a = 1; a < N; ++a) {
        current_refresh = std::max(current_refresh, idx[first_valid[a]]);
    }

    // 每个资产的当前位置指针
    std::vector<Size> pos(N, 0);
    for (Size a = 0; a < N; ++a) {
        pos[a] = first_valid[a];
    }

    while (true) {
        // 在当前 refresh time, 取各资产最近的观测值 (<= current_refresh)
        std::vector<Real> row(N);
        bool all_valid = true;
        for (Size a = 0; a < N; ++a) {
            // 前进到不超过 current_refresh 的最后一个观测
            while (pos[a] + 1 < M && !std::isnan(x[a][pos[a] + 1]) &&
                   idx[pos[a] + 1] <= current_refresh) {
                ++pos[a];
            }
            if (std::isnan(x[a][pos[a]])) {
                all_valid = false;
                break;
            }
            row[a] = x[a][pos[a]];
        }

        if (!all_valid) break;

        refresh_idx.push_back(current_refresh);
        for (Size a = 0; a < N; ++a) {
            refresh_data[a].push_back(row[a]);
        }

        // 下一个 refresh time = min(各资产在当前 refresh time 之后的下一个观测时间)
        Real next_refresh = std::numeric_limits<Real>::infinity();
        for (Size a = 0; a < N; ++a) {
            Size next_pos = pos[a] + 1;
            while (next_pos < M && std::isnan(x[a][next_pos])) ++next_pos;
            if (next_pos < M) {
                next_refresh = std::min(next_refresh, idx[next_pos]);
            }
        }
        if (!std::isfinite(next_refresh)) break;  // 所有资产已耗尽
        current_refresh = next_refresh;

        // 前进各资产指针到 >= current_refresh
        for (Size a = 0; a < N; ++a) {
            while (pos[a] + 1 < M && !std::isnan(x[a][pos[a] + 1]) &&
                   idx[pos[a] + 1] < current_refresh) {
                ++pos[a];
            }
        }
    }

    return {refresh_data, refresh_idx};
}

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
