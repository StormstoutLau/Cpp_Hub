// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.1 任务 1.3 - 数据容器类型 (ADR-013 双层 linalg)
// 用途: 计量模型的横截面/面板/时间序列数据结构 + 便利构造函数
// 依赖: linalg::dynamic::MatrixXD/VectorXD (ADR-013) + types.hpp
// 约定: 头文件 #include 必须位于 namespace 外
#pragma once

#include <string>
#include <variant>
#include <vector>

#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// 别名引用 linalg::dynamic 层类型 (ADR-013)
// 在 econometrics 命名空间内通过 linalg::dynamic::MatrixXD 使用
using linalg::dynamic::MatrixXD;
using linalg::dynamic::VectorXD;

// =============================================================================
// CrossSectionData - 横截面数据
//   X: n×k 解释变量矩阵, y: n 维响应向量
// =============================================================================
struct CrossSectionData {
    MatrixXD X;                    ///< 设计矩阵 (n × k)
    VectorXD y;                    ///< 响应向量 (n)
    std::vector<std::string> x_names;  ///< 解释变量名称 (k)
    std::string y_name;                ///< 响应变量名称
};

// =============================================================================
// PanelData - 面板 (纵向) 数据
// =============================================================================
struct PanelData {
    MatrixXD X;                    ///< 设计矩阵 (N × k)
    VectorXD y;                    ///< 响应向量 (N)
    std::vector<Index> entity_id;  ///< 实体标识 (每个观测所属截面单元)
    std::vector<Index> time_id;    ///< 时间标识 (每个观测所属时期)
    std::vector<std::string> x_names;  ///< 解释变量名称 (k)
    std::string y_name;                ///< 响应变量名称
    bool balanced = false;          ///< 是否平衡面板 (各实体时期数相同)
};

// =============================================================================
// TimeSeriesData - 时间序列数据
// =============================================================================
struct TimeSeriesData {
    VectorXD y;                    ///< 响应向量 (n)
    MatrixXD X;                    ///< 设计矩阵 (n × k)
    std::vector<Real> timestamps;  ///< 时间戳 (n)
    std::vector<std::string> x_names;  ///< 解释变量名称 (k)
    std::string y_name;                ///< 响应变量名称
};

// =============================================================================
// EconData - 统一数据载体 (variant 分发)
// =============================================================================
using EconData = std::variant<CrossSectionData, PanelData, TimeSeriesData>;

/// @brief 构造 CrossSectionData
/// @param X 设计矩阵 (n × k)
/// @param y 响应向量 (n)
/// @param x_names 解释变量名称 (k)
/// @param y_name 响应变量名称
/// @return 填充完成的 CrossSectionData
inline CrossSectionData make_cross_section(const MatrixXD& X, const VectorXD& y,
                                           std::vector<std::string> x_names,
                                           std::string y_name) {
    CrossSectionData d;
    d.X = X;
    d.y = y;
    d.x_names = std::move(x_names);
    d.y_name = std::move(y_name);
    return d;
}

/// @brief 构造 PanelData
/// @param X 设计矩阵 (N × k)
/// @param y 响应向量 (N)
/// @param entity_id 实体标识 (N)
/// @param time_id 时间标识 (N)
/// @param x_names 解释变量名称 (k)
/// @param y_name 响应变量名称
/// @param balanced 是否平衡面板
/// @return 填充完成的 PanelData
inline PanelData make_panel(const MatrixXD& X, const VectorXD& y,
                            std::vector<Index> entity_id, std::vector<Index> time_id,
                            std::vector<std::string> x_names, std::string y_name,
                            bool balanced) {
    PanelData d;
    d.X = X;
    d.y = y;
    d.entity_id = std::move(entity_id);
    d.time_id = std::move(time_id);
    d.x_names = std::move(x_names);
    d.y_name = std::move(y_name);
    d.balanced = balanced;
    return d;
}

/// @brief 构造 TimeSeriesData
/// @param y 响应向量 (n)
/// @param X 设计矩阵 (n × k)
/// @param timestamps 时间戳 (n)
/// @param x_names 解释变量名称 (k)
/// @param y_name 响应变量名称
/// @return 填充完成的 TimeSeriesData
inline TimeSeriesData make_time_series(const VectorXD& y, const MatrixXD& X,
                                       std::vector<Real> timestamps,
                                       std::vector<std::string> x_names,
                                       std::string y_name) {
    TimeSeriesData d;
    d.y = y;
    d.X = X;
    d.timestamps = std::move(timestamps);
    d.x_names = std::move(x_names);
    d.y_name = std::move(y_name);
    return d;
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub