// SOURCE: PHASE6_IMPLEMENTATION_PLAN §6.2 任务 4.5 - Cluster Bootstrap (Liu 1988 / Cameron-Gelbach-Miller 2008)
// 教材锚点:
//   - Liu (1988) "Bootstrap procedures under some non-i.i.d. models" Ann. Stat.
//   - Cameron-Gelbach-Miller (2008) REStat "Bootstrap-based improvements for inference
//     with clustered errors" (CGM 2008, cluster bootstrap vs wild cluster bootstrap)
//   - Cameron-Miller (2015) "A Practitioner's Guide to Cluster-Robust Inference"
//   - Davison-Hinkley (1997) §3.9 (cluster bootstrap)
//   - R multiwayvcov::cluster.boot(boot_type='xy') (pairs cluster bootstrap)
//
// 算法 (Pairs Cluster Bootstrap, CGM 2008):
//   1. 原样本: (y_i, x_i), i = 1..N, 分属 G 个聚类
//      聚类 g 包含观测集 I_g, |I_g| = n_g, Σ n_g = N
//   2. 原估计: θ̂ = estimator.estimate(data)
//   3. 对 b = 1..B:
//      a. 有放回采样 G 个聚类索引: g*_1, ..., g*_G ~ Uniform{1, ..., G}
//      b. 构造重采样数据: 按采样顺序拼接各聚类的全部观测
//         (X*_b, y*_b) = concat( (X_{g*_1}, y_{g*_1}), ..., (X_{g*_G}, y_{g*_G}) )
//      c. 重采样观测数 N*_b = Σ_{j=1}^{G} n_{g*_j} (不等大小聚类时 N*_b ≠ N)
//      d. 估计 θ̂*_b = estimator.estimate(data*_b)
//   4. Bootstrap 协方差: V_boot = (1/(B-1)) Σ (θ̂*_b - θ̄*)(θ̂*_b - θ̄*)'
//   5. 置信区间: percentile CI (默认)
//
// 排幻觉点:
//   H-023: Cluster Bootstrap 采样聚类 (G 个有放回), 不是观测 (N 个独立采样)
//          (独立采样观测会破坏聚类内相关性, 导致 V_boot 严重偏小)
//   H-024: 同一聚类的全部观测必须整体保留 (cluster g 被选中 → 所有 |I_g| 个观测一起进入重采样)
//          (拆分聚类内的观测会破坏聚类结构)
//   H-025: G < 20 时 Cluster Bootstrap 有限样本性质差 (Cameron-Miller 2015 推荐 G >= 20)
//          (G < 10 时推荐 Wild Cluster Bootstrap, v1.6+ 实现)
//          (当前实现: G < 20 时设置 warning flag, 不阻止运行)
//   H-026: 不等大小聚类时 N*_b ≠ N (重采样观测数为 Σ n_{g*_j}, 非 N)
//          (强制 N*_b = N 会错误截断或重复某些聚类的观测)
//
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训)
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/data_types.hpp"
#include "cpphub/econometrics/core/estimation_result.hpp"
#include "cpphub/econometrics/core/estimator_base.hpp"
#include "cpphub/econometrics/resampling/bootstrap_base.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

using linalg::dynamic::MatrixXD;
using linalg::dynamic::VectorXD;

// =============================================================================
// ClusterBootstrap - 聚类 (Pairs) Bootstrap 引擎 (CGM 2008 boot_type='xy')
//
// 适用场景:
//   - 面板数据 (PanelData, 聚类=entity)
//   - 横截面数据带聚类结构 (CrossSectionData + 显式 cluster_id)
//   - 聚类数 G >= 20 (Cameron-Miller 2015 推荐)
//
// 不适用场景:
//   - 时间序列 (用 BlockBootstrap)
//   - G < 10 (推荐 Wild Cluster Bootstrap, v1.6+)
//   - 无聚类结构 (用 PairedBootstrap)
//
// 用法 1 (PanelData, 聚类=entity):
//   ClusterBootstrap engine;
//   OLSEstimator ols;
//   BootstrapResult r = engine.resample(ols, panel_data, 999, 42);
//
// 用法 2 (CrossSectionData + 显式 cluster_id):
//   ClusterBootstrap engine;
//   engine.setClusterId(cluster_ids);  // 0-based 或任意整数 ID
//   OLSEstimator ols;
//   BootstrapResult r = engine.resample(ols, cs_data, 999, 42);
//
// 用法 3 (PanelData, 聚类=time):
//   ClusterBootstrap engine;
//   engine.setClusterDimension("time");
//   BootstrapResult r = engine.resample(ols, panel_data, 999, 42);
// =============================================================================
class ClusterBootstrap : public BootstrapEngine {
public:
    /// @brief 默认构造 (聚类维度="entity", alpha=0.05, seed=42)
    ClusterBootstrap() = default;

    /// @brief 设置显式聚类 ID (用于 CrossSectionData)
    /// @param cluster_id 聚类标识向量 (长度 N, 任意整数 ID, 内部重映射为 0-based 紧凑索引)
    void setClusterId(const std::vector<Index>& cluster_id) {
        explicit_cluster_id_ = cluster_id;
        has_explicit_cluster_id_ = true;
    }

    /// @brief 设置聚类维度 (仅对 PanelData 有效)
    /// @param dim "entity" (默认) 或 "time"
    /// @throws std::invalid_argument 若 dim 非 "entity"/"time"
    void setClusterDimension(const std::string& dim) {
        if (dim != "entity" && dim != "time") {
            throw std::invalid_argument(
                "ClusterBootstrap::setClusterDimension: dim must be 'entity' or 'time', got '"
                + dim + "'");
        }
        cluster_dim_ = dim;
    }

    /// @brief 获取聚类数 G (resample 后填充, resample 前为 0)
    Size nClusters() const { return n_clusters_; }

    /// @brief 是否触发了 G < 20 警告 (resample 后填充)
    bool smallClusterWarning() const { return small_cluster_warning_; }

    /// @brief 执行 Cluster Bootstrap 重采样
    /// @param estimator 估计器 (调用其 estimate() 方法)
    /// @param data 原始数据 (PanelData 或 CrossSectionData)
    /// @param n_replicates 重抽样次数 B (默认 999)
    /// @param seed 随机种子 (默认 42)
    /// @return Bootstrap 结果
    BootstrapResult resample(Estimator& estimator, const EconData& data,
                              Size n_replicates = 999,
                              std::uint64_t seed = 42) override {
        // -----------------------------------------------------------------
        // 1. 提取 (X, y) 和聚类 ID
        // -----------------------------------------------------------------
        MatrixXD X_orig;
        VectorXD y_orig;
        std::vector<Index> cluster_id;
        std::vector<std::string> x_names;
        std::string y_name;

        if (std::holds_alternative<PanelData>(data)) {
            const auto& panel = std::get<PanelData>(data);
            X_orig = panel.X;
            y_orig = panel.y;
            x_names = panel.x_names;
            y_name = panel.y_name;

            if (has_explicit_cluster_id_) {
                cluster_id = explicit_cluster_id_;
            } else {
                // 从 PanelData 提取聚类 ID (entity 或 time)
                const std::vector<Index>& src = (cluster_dim_ == "entity")
                    ? panel.entity_id : panel.time_id;
                // 重映射为 0-based 紧凑索引 (排幻觉点 P2: 不假设原始 ID 从 0 开始)
                cluster_id = remap_to_compact(src);
            }
        } else if (std::holds_alternative<CrossSectionData>(data)) {
            const auto& cs = std::get<CrossSectionData>(data);
            X_orig = cs.X;
            y_orig = cs.y;
            x_names = cs.x_names;
            y_name = cs.y_name;

            if (has_explicit_cluster_id_) {
                cluster_id = explicit_cluster_id_;
            } else {
                throw std::invalid_argument(
                    "ClusterBootstrap::resample: CrossSectionData requires explicit "
                    "cluster_id (call setClusterId() first)");
            }
        } else {
            throw std::invalid_argument(
                "ClusterBootstrap::resample: only PanelData/CrossSectionData supported");
        }

        const Size N = X_orig.rows();
        const Size k = X_orig.cols();

        if (N < 2) {
            throw std::invalid_argument(
                "ClusterBootstrap::resample: requires N >= 2 observations");
        }
        if (k == 0) {
            throw std::invalid_argument(
                "ClusterBootstrap::resample: design matrix has 0 columns");
        }
        if (cluster_id.size() != N) {
            throw std::invalid_argument(
                "ClusterBootstrap::resample: cluster_id size ("
                + std::to_string(cluster_id.size())
                + ") != N (" + std::to_string(N) + ")");
        }
        if (n_replicates < 2) {
            throw std::invalid_argument(
                "ClusterBootstrap::resample: requires n_replicates >= 2");
        }

        // -----------------------------------------------------------------
        // 2. 按聚类分组观测索引
        // -----------------------------------------------------------------
        // groups[g] = {观测索引列表}, g = 0, 1, ..., G-1
        std::unordered_map<Index, std::vector<Size>> group_map;
        std::vector<Index> cluster_order;  // 按首次出现顺序保留聚类 ID
        for (Size i = 0; i < N; ++i) {
            const Index g = cluster_id[i];
            if (group_map.find(g) == group_map.end()) {
                cluster_order.push_back(g);
            }
            group_map[g].push_back(i);
        }
        const Size G = cluster_order.size();
        n_clusters_ = G;

        if (G < 2) {
            throw std::invalid_argument(
                "ClusterBootstrap::resample: requires G >= 2 clusters (got G=1, "
                "cannot resample with only one cluster)");
        }

        // H-025: G < 20 警告 (Cameron-Miller 2015)
        small_cluster_warning_ = (G < 20);

        // 构建聚类大小向量和累计偏移 (用于快速拼接)
        std::vector<Size> cluster_sizes(G);
        for (Size g = 0; g < G; ++g) {
            cluster_sizes[g] = group_map[cluster_order[g]].size();
        }

        // -----------------------------------------------------------------
        // 3. Bootstrap 重采样主循环 (H-023, H-024, H-026)
        // -----------------------------------------------------------------
        std::vector<VectorXD> samples;
        samples.reserve(n_replicates);
        Size n_failed = 0;

        for (Size b = 0; b < n_replicates; ++b) {
            Philox4x64 rng(seed, b);

            // H-023: 采样 G 个聚类索引 (有放回, 聚类级采样)
            std::vector<Size> sampled_clusters(G);
            for (Size g = 0; g < G; ++g) {
                sampled_clusters[g] = randint(rng, G);
            }

            // H-026: 计算重采样总观测数 N*_b = Σ n_{g*_j}
            Size N_b = 0;
            for (Size g = 0; g < G; ++g) {
                N_b += cluster_sizes[sampled_clusters[g]];
            }

            // H-024: 拼接各聚类的全部观测 (整体保留)
            MatrixXD X_b(N_b, k);
            VectorXD y_b(N_b);
            Size pos = 0;
            for (Size g = 0; g < G; ++g) {
                const Size src_cluster = sampled_clusters[g];
                const std::vector<Size>& obs = group_map[cluster_order[src_cluster]];
                for (Size idx : obs) {
                    y_b(pos) = y_orig(idx);
                    for (Size j = 0; j < k; ++j) {
                        X_b(pos, j) = X_orig(idx, j);
                    }
                    ++pos;
                }
            }

            // 估计 θ̂*_b
            try {
                EconData data_b = make_cross_section(X_b, y_b, x_names, y_name);
                EstimationResult r_b = estimator.estimate(data_b);

                if (r_b.coefficients.size() != k ||
                    !r_b.coefficients.eigen().allFinite()) {
                    ++n_failed;
                    continue;
                }

                samples.push_back(r_b.coefficients);
            } catch (const std::exception&) {
                ++n_failed;
            }
        }

        // -----------------------------------------------------------------
        // 4. 计算 Bootstrap 统计量
        // -----------------------------------------------------------------
        BootstrapResult result;
        result.n_replicates = n_replicates;
        result.n_failed = n_failed;
        result.bootstrap_samples = samples;

        if (samples.empty()) {
            result.coef_mean = VectorXD(k);
            result.coef_std = VectorXD(k);
            result.coef_vcov = MatrixXD(k, k);
            for (Size i = 0; i < k; ++i) {
                result.coef_mean(i) = std::numeric_limits<Real>::quiet_NaN();
                result.coef_std(i) = std::numeric_limits<Real>::quiet_NaN();
            }
            result.lower_ci = std::numeric_limits<Real>::quiet_NaN();
            result.upper_ci = std::numeric_limits<Real>::quiet_NaN();
            return result;
        }

        compute_stats(samples, result.coef_mean, result.coef_vcov);

        result.coef_std = VectorXD(k);
        for (Size i = 0; i < k; ++i) {
            result.coef_std(i) = std::sqrt(result.coef_vcov(i, i));
        }

        if (k >= 1) {
            std::vector<Real> ci_samples;
            ci_samples.reserve(samples.size());
            for (const auto& s : samples) {
                ci_samples.push_back(s(0));
            }
            auto [lo, hi] = percentileCI(ci_samples, alpha_);
            result.lower_ci = lo;
            result.upper_ci = hi;
        }

        return result;
    }

    std::string name() const override { return "ClusterBootstrap"; }

private:
    std::vector<Index> explicit_cluster_id_;
    bool has_explicit_cluster_id_ = false;
    std::string cluster_dim_ = "entity";
    Size n_clusters_ = 0;
    bool small_cluster_warning_ = false;

    /// @brief 将任意整数 ID 重映射为 0-based 紧凑索引
    /// (排幻觉点 P2: 不假设原始 ID 从 0 开始或连续)
    static std::vector<Index> remap_to_compact(const std::vector<Index>& src) {
        std::vector<Index> result;
        result.reserve(src.size());
        std::unordered_map<Index, Index> mapping;
        Index next_id = 0;
        for (Index v : src) {
            auto it = mapping.find(v);
            if (it == mapping.end()) {
                mapping[v] = next_id;
                result.push_back(next_id);
                ++next_id;
            } else {
                result.push_back(it->second);
            }
        }
        return result;
    }
};

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
