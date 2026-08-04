// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.4 Day 7-8 任务 1.10 - 聚类稳健协方差 (排幻觉点 E6)
// 教材锚点:
//   - Liang-Zeger (1986) Biometrika (原始聚类 SE)
//   - Cameron-Gelbach-Miller (2011) 双向聚类 (REStat)
//   - Cameron-Miller (2015) "A Practitioner's Guide to Cluster-Robust Inference"
//   - R sandwich::vcovCL 源码 (L80-100, type="HC1" 默认含小样本调整)
//
// 排幻觉点 (R sandwich 对照):
//   E6: 小样本调整 G/(G-1)·(N-1)/(N-K), 非 1/(G-1) 或 G/(G-1)
//       R vcovCL L80-100 实测: type="HC1" 默认含此调整, 分母 (G-1)·(N-K)
//       文档锚点: ?sandwich::vcovCL "type: ... HC1 ... applies the cluster
//       correction (G/(G-1)) and the degrees-of-freedom correction ((N-1)/(N-K))"
//   双向聚类: V_twoway = V(g1) + V(g2) - V(g1∩g2), 每个分量独立做小样本调整
//       (Cameron-Gelbach-Miller 2011 §2.3, R vcovCL cluster=list(g1,g2))
//
// 公式 (Liang-Zeger 1986, 严格遵循):
//   V_cluster = (X'X)^{-1} · meat · (X'X)^{-1} · [G/(G-1)·(N-1)/(N-K)]
//   meat = Σ_g X_g' · ε_g · ε_g' · X_g = Σ_g (X_g' ε_g)(X_g' ε_g)'
//   其中 X_g, ε_g 为第 g 个聚类的子矩阵/子向量, G = 聚类数
//
// 双向聚类 (Cameron-Gelbach-Miller 2011):
//   V_twoway = V_cluster(g1) + V_cluster(g2) - V_cluster(g1∩g2)
//   g1∩g2 = 按 (g1, g2) 组合的合并聚类, 每个分量独立做小样本调整
//
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训)
#pragma once

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/covariance_type.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

using linalg::dynamic::MatrixXD;
using linalg::dynamic::VectorXD;

namespace detail {

/// @brief 计算单向聚类 meat 矩阵 (无小样本调整)
/// @param X 设计矩阵 (N × k)
/// @param residuals 残差向量 (N)
/// @param cluster_id 0-based 聚类 ID (长度 N)
/// @param G 聚类数 (唯一聚类数)
/// @return meat 矩阵 (k × k) = Σ_g (X_g' ε_g)(X_g' ε_g)'
inline MatrixXD compute_cluster_meat(const MatrixXD& X, const VectorXD& residuals,
                                       const std::vector<Index>& cluster_id, Index G) {
    const Size N = X.rows();
    const Size k = X.cols();

    // 按聚类分组观测索引
    std::unordered_map<Index, std::vector<Size>> groups;
    for (Size i = 0; i < N; ++i) {
        groups[cluster_id[i]].push_back(i);
    }

    // meat = Σ_g (X_g' ε_g)(X_g' ε_g)'  (秩 1 累加)
    Eigen::MatrixXd meat = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(k),
                                                   static_cast<Eigen::Index>(k));
    for (const auto& kv : groups) {
        const std::vector<Size>& idx = kv.second;
        const Size ng = idx.size();
        // X_g (ng × k), ε_g (ng)
        // v_g = X_g' · ε_g (k 维列向量)
        Eigen::VectorXd v_g = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(k));
        for (Size i = 0; i < ng; ++i) {
            const Size row = idx[i];
            for (Size j = 0; j < k; ++j) {
                v_g(static_cast<Eigen::Index>(j)) += X(row, j) * residuals(row);
            }
        }
        // meat += v_g · v_g'  (k×k outer product)
        meat += v_g * v_g.transpose();
    }
    (void)G;  // G 不参与 meat 计算, 仅用于小样本调整
    return MatrixXD(meat);
}

/// @brief 单向聚类协方差核心计算 (含小样本调整 E6)
/// @throws std::invalid_argument 若 G < 2 或 N <= K 或维度不匹配
inline MatrixXD compute_cluster_vcov_oneway(const MatrixXD& X, const VectorXD& residuals,
                                              const MatrixXD& XtX_inv,
                                              const std::vector<Index>& cluster_id) {
    const Size N = X.rows();
    const Size k = X.cols();

    // 维度校验
    if (residuals.size() != N) {
        throw std::invalid_argument(
            "compute_cluster_vcov: residuals size (" + std::to_string(residuals.size()) +
            ") != X.rows() (" + std::to_string(N) + ")");
    }
    if (XtX_inv.rows() != k || XtX_inv.cols() != k) {
        throw std::invalid_argument("compute_cluster_vcov: XtX_inv must be k×k");
    }
    if (cluster_id.size() != N) {
        throw std::invalid_argument(
            "compute_cluster_vcov: cluster_id size (" + std::to_string(cluster_id.size()) +
            ") != X.rows() (" + std::to_string(N) + ")");
    }
    if (N <= k) {
        throw std::invalid_argument(
            "compute_cluster_vcov: requires N > K (degrees of freedom)");
    }

    // 计算唯一聚类数 G
    std::unordered_map<Index, Index> unique_clusters;
    for (Index id : cluster_id) {
        unique_clusters[id] = 1;
    }
    const Index G = static_cast<Index>(unique_clusters.size());

    // 排幻觉点 E6: 分母 (G-1), G=1 时不可计算
    if (G < 2) {
        throw std::invalid_argument(
            "compute_cluster_vcov: requires G >= 2 clusters (got G=1, "
            "small-sample adjustment G/(G-1) undefined)");
    }

    // meat = Σ_g (X_g' ε_g)(X_g' ε_g)'
    const MatrixXD meat = compute_cluster_meat(X, residuals, cluster_id, G);

    // V_unscaled = (X'X)^{-1} · meat · (X'X)^{-1}
    const Eigen::MatrixXd V_unscaled =
        XtX_inv.eigen() * meat.eigen() * XtX_inv.eigen();

    // 排幻觉点 E6: 小样本调整 G/(G-1)·(N-1)/(N-K)
    // R sandwich::vcovCL type="HC1" 默认含此调整
    const Real scale = static_cast<Real>(G) / static_cast<Real>(G - 1) *
                       static_cast<Real>(N - 1) / static_cast<Real>(N - k);

    Eigen::MatrixXd V = V_unscaled * scale;
    return MatrixXD(V);
}

}  // namespace detail

/// @brief 合并双向聚类 ID (g1, g2) → 唯一组合的 0-based 紧凑 ID
/// @param cluster_id1 第一维聚类 ID (长度 N, 0-based)
/// @param cluster_id2 第二维聚类 ID (长度 N, 0-based)
/// @return 合并后的 0-based 紧凑 ID (长度 N), 每个唯一 (g1,g2) 组合映射到一个整数
/// @throws std::invalid_argument 若两向量长度不一致
inline std::vector<Index> merge_twoway_cluster_ids(const std::vector<Index>& cluster_id1,
                                                     const std::vector<Index>& cluster_id2) {
    if (cluster_id1.size() != cluster_id2.size()) {
        throw std::invalid_argument(
            "merge_twoway_cluster_ids: size mismatch (" +
            std::to_string(cluster_id1.size()) + " vs " +
            std::to_string(cluster_id2.size()) + ")");
    }
    std::vector<Index> merged;
    merged.reserve(cluster_id1.size());
    std::map<std::pair<Index, Index>, Index> mapping;
    Index next_id = 0;
    for (Size i = 0; i < cluster_id1.size(); ++i) {
        const auto key = std::make_pair(cluster_id1[i], cluster_id2[i]);
        auto it = mapping.find(key);
        if (it == mapping.end()) {
            mapping[key] = next_id;
            merged.push_back(next_id);
            ++next_id;
        } else {
            merged.push_back(it->second);
        }
    }
    return merged;
}

/// @brief 计算聚类稳健协方差矩阵 (Liang-Zeger 1986, 含双向聚类)
/// @param X 设计矩阵 (N × k)
/// @param residuals 残差向量 (N)
/// @param XtX_inv (X'X)^{-1} (k × k)
/// @param cluster_id 聚类 ID (长度 N, 0-based, 可非连续)
/// @param twoway 是否双向聚类 (Cameron-Gelbach-Miller 2011)
/// @param cluster_id2 第二维聚类 ID (仅 twoway=true 时使用, 长度 N)
/// @return 聚类稳健协方差矩阵 (k × k)
/// @throws std::invalid_argument 若 G < 2 / N <= K / 维度不匹配 / twoway=true 但 cluster_id2 长度不匹配
///
/// 公式 (单向, Liang-Zeger 1986):
///   V_cluster = (X'X)^{-1} · [Σ_g X_g' ε_g ε_g' X_g] · (X'X)^{-1} · [G/(G-1)·(N-1)/(N-K)]
/// 公式 (双向, Cameron-Gelbach-Miller 2011):
///   V_twoway = V(g1) + V(g2) - V(g1∩g2)
///   每个分量独立做小样本调整 (G1, G2, G12 各自的 G/(G-1)·(N-1)/(N-K))
inline MatrixXD compute_cluster_vcov(const MatrixXD& X, const VectorXD& residuals,
                                       const MatrixXD& XtX_inv,
                                       const std::vector<Index>& cluster_id,
                                       bool twoway = false,
                                       const std::vector<Index>& cluster_id2 = {}) {
    if (!twoway) {
        // 单向聚类
        return detail::compute_cluster_vcov_oneway(X, residuals, XtX_inv, cluster_id);
    }

    // 双向聚类 (Cameron-Gelbach-Miller 2011)
    const Size N = X.rows();
    if (cluster_id2.size() != N) {
        throw std::invalid_argument(
            "compute_cluster_vcov: cluster_id2 size (" +
            std::to_string(cluster_id2.size()) +
            ") != X.rows() (" + std::to_string(N) + ") in twoway mode");
    }

    // V(g1)
    const MatrixXD V_g1 = detail::compute_cluster_vcov_oneway(X, residuals, XtX_inv, cluster_id);
    // V(g2)
    const MatrixXD V_g2 = detail::compute_cluster_vcov_oneway(X, residuals, XtX_inv, cluster_id2);
    // V(g1∩g2) - 合并聚类
    const std::vector<Index> g12 = merge_twoway_cluster_ids(cluster_id, cluster_id2);
    const MatrixXD V_g12 = detail::compute_cluster_vcov_oneway(X, residuals, XtX_inv, g12);

    // V_twoway = V(g1) + V(g2) - V(g1∩g2)
    const Eigen::MatrixXd V_twoway =
        V_g1.eigen() + V_g2.eigen() - V_g12.eigen();
    return MatrixXD(V_twoway);
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
