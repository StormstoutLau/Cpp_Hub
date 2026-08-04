// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.1 任务 1.2 - CovarianceType 枚举
// 用途: 计量估计结果与 Estimator 协方差类型标记
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训, 避免嵌套 namespace 错误)
#pragma once

#include <string>

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// =============================================================================
// CovarianceType - 协方差矩阵估计类型
// =============================================================================
enum class CovarianceType {
    Classical,                    ///< 经典 OLS (X'X)^{-1} s²
    HC0,                          ///< White (1980) 异方差稳健
    HC1,                          ///< HC0 * n/(n-k)
    HC2,                          ///< MacKinnon & White (1985)
    HC3,                          ///< jackknife 近似
    HC4,                          ///< Cribari-Neto (2004)
    HC5,                          ///< Cribari-Neto (2007)
    HAC_Bartlett,                 ///< Newey-West (Bartlett 核)
    HAC_QuadraticSpectral,        ///< Andrews (1991) QS 核
    HAC_Parzen,                   ///< Parzen 核
    HAC_TukeyHanning,             ///< Tukey-Hanning 核
    Cluster_OneWay,               ///< 单维聚类稳健
    Cluster_TwoWay,               ///< 二维聚类稳健
    OPG,                          ///< Outer Product of Gradients
    Hessian,                      ///< Hessian 逆
    Sandwich,                     ///< 三明治估计
    Bootstrap,                    ///< Bootstrap 协方差
    Custom                        ///< 用户自定义
};

/// @brief 将 CovarianceType 枚举转换为字符串名
/// @param type 协方差类型枚举值
/// @return 对应的字符串名 (如 "HC1")
inline std::string to_string(CovarianceType type) {
    switch (type) {
        case CovarianceType::Classical:            return "Classical";
        case CovarianceType::HC0:                  return "HC0";
        case CovarianceType::HC1:                  return "HC1";
        case CovarianceType::HC2:                  return "HC2";
        case CovarianceType::HC3:                  return "HC3";
        case CovarianceType::HC4:                  return "HC4";
        case CovarianceType::HC5:                  return "HC5";
        case CovarianceType::HAC_Bartlett:         return "HAC_Bartlett";
        case CovarianceType::HAC_QuadraticSpectral: return "HAC_QuadraticSpectral";
        case CovarianceType::HAC_Parzen:           return "HAC_Parzen";
        case CovarianceType::HAC_TukeyHanning:     return "HAC_TukeyHanning";
        case CovarianceType::Cluster_OneWay:       return "Cluster_OneWay";
        case CovarianceType::Cluster_TwoWay:       return "Cluster_TwoWay";
        case CovarianceType::OPG:                  return "OPG";
        case CovarianceType::Hessian:              return "Hessian";
        case CovarianceType::Sandwich:             return "Sandwich";
        case CovarianceType::Bootstrap:            return "Bootstrap";
        case CovarianceType::Custom:               return "Custom";
    }
    return "Unknown";
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub