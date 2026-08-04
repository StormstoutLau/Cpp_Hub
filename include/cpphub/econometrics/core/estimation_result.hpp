// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.1 任务 1.4 - 估计结果类型
// 用途: 估计器输出/推断结果/Bootstrap 结果的数据载体
// 依赖: linalg::dynamic::MatrixXD/VectorXD (ADR-013) + covariance_type.hpp
// 约定: 头文件 #include 必须位于 namespace 外
#pragma once

#include <vector>

#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/covariance_type.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

using linalg::dynamic::MatrixXD;
using linalg::dynamic::VectorXD;

// =============================================================================
// EstimationResult - 基础估计结果
// =============================================================================
struct EstimationResult {
    VectorXD coefficients;        ///< 系数估计 (k)
    VectorXD std_errors;          ///< 标准误 (k)
    VectorXD t_statistics;        ///< t 统计量 (k)
    VectorXD p_values;            ///< p 值 (k)
    MatrixXD vcov;                ///< 协方差矩阵 (k × k)
    Real log_likelihood = 0.0;    ///< 对数似然值
    Real r_squared = 0.0;         ///< R²
    Real adj_r_squared = 0.0;     ///< 调整后 R²
    Size n_obs = 0;               ///< 观测数 n
    Size n_params = 0;            ///< 参数个数 k
    Size df_residual = 0;         ///< 残差自由度 n - k
    CovarianceType cov_type = CovarianceType::Classical;  ///< 协方差类型
};

// =============================================================================
// InferenceResult - 推断扩展结果 (继承 EstimationResult)
//   增加 Wald / LR / LM 三大检验统计量
// =============================================================================
struct InferenceResult : EstimationResult {
    Real wald_statistic = 0.0;    ///< Wald 检验统计量
    Real lr_statistic = 0.0;      ///< 似然比检验统计量
    Real lm_statistic = 0.0;      ///< 拉格朗日乘数 (score) 检验统计量
    Real wald_pvalue = 0.0;       ///< Wald 检验 p 值
    Real lr_pvalue = 0.0;         ///< LR 检验 p 值
    Real lm_pvalue = 0.0;         ///< LM 检验 p 值
    Size df_test = 0;             ///< 检验自由度
};

// =============================================================================
// BootstrapResult - Bootstrap 重抽样结果
// =============================================================================
struct BootstrapResult {
    VectorXD coef_mean;                          ///< Bootstrap 系数均值 (k)
    VectorXD coef_std;                           ///< Bootstrap 系数标准差 (k)
    MatrixXD coef_vcov;                          ///< Bootstrap 系数协方差 (k × k)
    std::vector<VectorXD> bootstrap_samples;     ///< 每次重抽样的系数样本
    Real lower_ci = 0.0;                         ///< 置信区间下界
    Real upper_ci = 0.0;                         ///< 置信区间上界
    Size n_replicates = 0;                       ///< 重抽样次数
    Size n_failed = 0;                           ///< 失败的抽样次数
};

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub