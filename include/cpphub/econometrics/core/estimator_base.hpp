// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.1 任务 1.5 + ADR-002 (Bridge + Virtual Constructor)
// 用途: 计量估计器抽象基类 + 估计器类别枚举
// 依赖: covariance_type / data_types / estimation_result + <memory> + <string>
// 约定: 头文件 #include 必须位于 namespace 外
#pragma once

#include <memory>
#include <string>

#include "cpphub/econometrics/core/covariance_type.hpp"
#include "cpphub/econometrics/core/data_types.hpp"
#include "cpphub/econometrics/core/estimation_result.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// =============================================================================
// EstimatorClass - 估计器类别
// =============================================================================
enum class EstimatorClass {
    Parametric,        ///< 参数化估计器 (OLS/MLE/GLS)
    Semiparametric,    ///< 半参数估计器 (GMM/两步法)
    Nonparametric,     ///< 非参数估计器 (核回归)
    MachineLearning    ///< 机器学习估计器 (LASSO/随机森林)
};

// =============================================================================
// Estimator - 计量估计器抽象基类 (ADR-002: 虚构造 + 接口桥接)
//   派生类必须实现 estimate / name / clone
// =============================================================================
class Estimator {
public:
    /// 虚析构 (多态删除安全)
    virtual ~Estimator() = default;

    /// @brief 对给定数据执行估计
    /// @param data 统一数据载体 (横截面/面板/时间序列)
    /// @return 估计结果
    virtual EstimationResult estimate(const EconData& data) = 0;

    /// @brief 估计器名称
    /// @return 名称字符串
    virtual std::string name() const = 0;

    /// @brief 估计器类别
    /// @return 默认 Parametric
    virtual EstimatorClass estimatorClass() const { return EstimatorClass::Parametric; }

    /// @brief 是否参数化估计器
    /// @return 默认 true
    virtual bool isParametric() const { return true; }

    /// @brief 是否半参数估计器
    /// @return 默认 false
    virtual bool isSemiparametric() const { return false; }

    /// @brief 是否非参数估计器
    /// @return 默认 false
    virtual bool isNonparametric() const { return false; }

    /// @brief 设置协方差类型
    /// @param type 协方差类型
    virtual void setCovarianceType(CovarianceType type) { cov_type_ = type; }

    /// @brief 获取当前协方差类型
    /// @return 协方差类型
    CovarianceType covarianceType() const { return cov_type_; }

    /// @brief 虚构造函数 (ADR-002): 克隆当前估计器配置
    /// @return 与当前对象等价的独立副本
    virtual std::unique_ptr<Estimator> clone() const = 0;

protected:
    CovarianceType cov_type_ = CovarianceType::Classical;  ///< 当前协方差类型
};

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub