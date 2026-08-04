// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.4 Day 7-8 任务 1.11 - 面板数据工具
// 用途: 面板数据长表/宽表转换, 平衡检测, 聚类标识提取
// 依赖: core/data_types.hpp (PanelData), linalg_dynamic (MatrixXD/VectorXD)
// 约定: 头文件 #include 必须位于 namespace 外
#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/data_types.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

using linalg::dynamic::MatrixXD;
using linalg::dynamic::VectorXD;

/// @brief 面板平衡检测结果
struct PanelBalance {
    bool is_balanced;       ///< 是否平衡面板
    Index n_entities;       ///< 唯一 entity 数
    Index n_periods;        ///< 平衡时为 T (每个 entity 的时期数), 非平衡时为 0
};

/// @brief 聚类 ID 提取结果
struct ClusterIdResult {
    std::vector<Index> ids;        ///< 0-based 紧凑聚类 ID (长度 N)
    Index n_clusters;              ///< 唯一聚类数 G
    std::vector<Index> original;   ///< 原始 ID (去重, 按首次出现顺序)
};

/// @brief 检测面板是否平衡
/// @param panel 面板数据 (含 entity_id, time_id)
/// @return PanelBalance
/// @throws std::invalid_argument 若 panel 为空 (N=0)
///
/// 平衡判定: 每个 entity 的 time_id 集合完全相同, 且 size = T
/// (排幻觉点 P1: 非平衡面板也可能 N == E*T 但 time_id 集合不同)
inline PanelBalance detect_balance(const PanelData& panel) {
    const Size N = panel.entity_id.size();
    if (N == 0) {
        throw std::invalid_argument("detect_balance: empty panel (N=0)");
    }
    if (panel.time_id.size() != N) {
        throw std::invalid_argument("detect_balance: entity_id/time_id size mismatch");
    }

    // 收集每个 entity 的 time_id 集合 (按首次出现顺序保留 entity)
    std::unordered_map<Index, std::unordered_set<Index>> entity_times;
    std::vector<Index> entity_order;  // 按首次出现顺序
    for (Size i = 0; i < N; ++i) {
        Index eid = panel.entity_id[i];
        Index tid = panel.time_id[i];
        if (entity_times.find(eid) == entity_times.end()) {
            entity_order.push_back(eid);
        }
        entity_times[eid].insert(tid);
    }

    const Index E = static_cast<Index>(entity_order.size());
    if (E == 0) {
        return {true, 0, 0};
    }

    // 单 entity: 平衡 (只要 time_id 无重复)
    if (E == 1) {
        const Index T = static_cast<Index>(entity_times[entity_order[0]].size());
        return {true, 1, T};
    }

    // 检查所有 entity 的 time_id 集合是否相同
    const std::unordered_set<Index>& ref_times = entity_times[entity_order[0]];
    for (Size i = 1; i < entity_order.size(); ++i) {
        const std::unordered_set<Index>& times = entity_times[entity_order[i]];
        if (times.size() != ref_times.size()) {
            return {false, E, 0};
        }
        // 检查集合相等
        for (Index tid : ref_times) {
            if (times.find(tid) == times.end()) {
                return {false, E, 0};
            }
        }
    }

    const Index T = static_cast<Index>(ref_times.size());
    return {true, E, T};
}

/// @brief 提取聚类标识 (entity_id 或 time_id 转为 0-based 紧凑索引)
/// @param panel 面板数据
/// @param type "entity" 或 "time"
/// @return ClusterIdResult
/// @throws std::invalid_argument 若 type 无效或 panel 为空
///
/// 排幻觉点 P2: 不能假设原始 ID 从 0 开始或连续, 必须重新映射
inline ClusterIdResult extract_cluster_ids(const PanelData& panel,
                                            const std::string& type) {
    const Size N = panel.entity_id.size();
    if (N == 0) {
        return {{}, 0, {}};
    }

    const std::vector<Index>& src = (type == "entity") ? panel.entity_id : panel.time_id;
    if (type != "entity" && type != "time") {
        throw std::invalid_argument(
            "extract_cluster_ids: type must be 'entity' or 'time', got '" + type + "'");
    }
    if (src.size() != N) {
        throw std::invalid_argument("extract_cluster_ids: source size mismatch");
    }

    std::vector<Index> ids;
    ids.reserve(N);
    std::unordered_map<Index, Index> mapping;
    std::vector<Index> original;
    Index next_id = 0;
    for (Size i = 0; i < N; ++i) {
        Index orig = src[i];
        auto it = mapping.find(orig);
        if (it == mapping.end()) {
            mapping[orig] = next_id;
            ids.push_back(next_id);
            original.push_back(orig);
            ++next_id;
        } else {
            ids.push_back(it->second);
        }
    }
    return {ids, next_id, original};
}

/// @brief 长表 → 宽表 (y 变量)
/// @param panel 面板数据
/// @return MatrixXD (n_entities × n_periods), 行=entity, 列=time
/// @throws std::invalid_argument 若 panel 为空或非平衡 (无法构造完整宽表)
///
/// 排幻觉点 P4: entity/time 顺序按首次出现, 不排序
/// 排幻觉点 P5: 非平衡面板的缺失值填 NaN, 不用 0/-1 哨兵值
inline MatrixXD reshape_y_to_wide(const PanelData& panel) {
    const Size N = panel.entity_id.size();
    if (N == 0) {
        throw std::invalid_argument("reshape_y_to_wide: empty panel");
    }
    const PanelBalance bal = detect_balance(panel);
    if (!bal.is_balanced) {
        throw std::invalid_argument(
            "reshape_y_to_wide: non-balanced panel not supported (use NaN-aware version)");
    }
    const Size E = static_cast<Size>(bal.n_entities);
    const Size T = static_cast<Size>(bal.n_periods);

    // 收集 entity/time 的首次出现顺序
    std::unordered_map<Index, Size> entity_row;
    std::unordered_map<Index, Size> time_col;
    std::vector<Index> entity_order;
    std::vector<Index> time_order;
    for (Size i = 0; i < N; ++i) {
        Index eid = panel.entity_id[i];
        Index tid = panel.time_id[i];
        if (entity_row.find(eid) == entity_row.end()) {
            entity_row[eid] = entity_order.size();
            entity_order.push_back(eid);
        }
        if (time_col.find(tid) == time_col.end()) {
            time_col[tid] = time_order.size();
            time_order.push_back(tid);
        }
    }

    MatrixXD wide(E, T);
    for (Size i = 0; i < N; ++i) {
        const Size r = entity_row[panel.entity_id[i]];
        const Size c = time_col[panel.time_id[i]];
        wide(r, c) = panel.y(i);
    }
    return wide;
}

/// @brief 长表 → 宽表 (X 变量, 指定列)
/// @param panel 面板数据
/// @param col_index X 的列索引 (0-based)
/// @return MatrixXD (n_entities × n_periods)
/// @throws std::invalid_argument 若 panel 为空或非平衡或 col_index 越界
inline MatrixXD reshape_x_to_wide(const PanelData& panel, Index col_index) {
    const Size N = panel.entity_id.size();
    if (N == 0) {
        throw std::invalid_argument("reshape_x_to_wide: empty panel");
    }
    if (col_index < 0 || static_cast<Size>(col_index) >= panel.X.cols()) {
        throw std::invalid_argument(
            "reshape_x_to_wide: col_index out of range");
    }
    const PanelBalance bal = detect_balance(panel);
    if (!bal.is_balanced) {
        throw std::invalid_argument(
            "reshape_x_to_wide: non-balanced panel not supported");
    }
    const Size E = static_cast<Size>(bal.n_entities);
    const Size T = static_cast<Size>(bal.n_periods);

    std::unordered_map<Index, Size> entity_row;
    std::unordered_map<Index, Size> time_col;
    std::vector<Index> entity_order;
    std::vector<Index> time_order;
    for (Size i = 0; i < N; ++i) {
        Index eid = panel.entity_id[i];
        Index tid = panel.time_id[i];
        if (entity_row.find(eid) == entity_row.end()) {
            entity_row[eid] = entity_order.size();
            entity_order.push_back(eid);
        }
        if (time_col.find(tid) == time_col.end()) {
            time_col[tid] = time_order.size();
            time_order.push_back(tid);
        }
    }

    MatrixXD wide(E, T);
    for (Size i = 0; i < N; ++i) {
        const Size r = entity_row[panel.entity_id[i]];
        const Size c = time_col[panel.time_id[i]];
        wide(r, c) = panel.X(i, col_index);
    }
    return wide;
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
