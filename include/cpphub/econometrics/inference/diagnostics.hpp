// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.2 任务 2.4 - 信息准则 (AIC/BIC/HQ)
// 教材锚点: Greene 8ed §5.3 (信息准则), Burnham-Anderson 2002 (模型选择)
// 公式:
//   AIC = 2K - 2ℓ                    (Akaike 1973)
//   BIC = K·log(N) - 2ℓ              (Schwarz 1978)
//   HQ  = 2K·log(log(N)) - 2ℓ        (Hannan-Quinn 1979)
// 注: AIC 小样本修正 AICc = AIC + 2K(K+1)/(N-K-1) (Sugiura 1978)
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训)
#pragma once

#include <cmath>
#include <stdexcept>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// =============================================================================
// InformationCriteria - 模型选择信息准则
// =============================================================================
struct InformationCriteria {
    Real aic;    ///< Akaike 信息准则 (越小越好)
    Real bic;    ///< Bayesian 信息准则 (越小越好)
    Real hq;     ///< Hannan-Quinn 信息准则 (越小越好)
    Real aicc;   ///< AIC 小样本修正 (Sugiura 1978)
};

/// @brief 计算信息准则 (AIC/BIC/HQ/AICc)
/// @param log_likelihood 对数似然值 ℓ
/// @param n_params 参数个数 K (含截距)
/// @param n_obs 观测数 N
/// @return InformationCriteria 结构体
/// @throws std::invalid_argument 若 n_obs == 0 或 n_params == 0
inline InformationCriteria compute_information_criteria(Real log_likelihood,
                                                         Size n_params,
                                                         Size n_obs) {
    if (n_obs == 0) {
        throw std::invalid_argument("compute_information_criteria: n_obs must be > 0");
    }
    if (n_params == 0) {
        throw std::invalid_argument("compute_information_criteria: n_params must be > 0");
    }

    InformationCriteria ic;
    const Real two_ell = 2.0 * log_likelihood;
    ic.aic = 2.0 * static_cast<Real>(n_params) - two_ell;
    ic.bic = static_cast<Real>(n_params) * std::log(static_cast<Real>(n_obs)) - two_ell;
    ic.hq = 2.0 * static_cast<Real>(n_params) * std::log(std::log(static_cast<Real>(n_obs))) - two_ell;

    // AICc 小样本修正 (Sugiura 1978): AICc = AIC + 2K(K+1)/(N-K-1)
    // 当 N - K - 1 <= 0 时未定义, 返回 AIC (并附 NaN 标记通过 aicc 为 AIC 值)
    if (n_obs > n_params + 1) {
        const Real correction = 2.0 * static_cast<Real>(n_params) *
                                static_cast<Real>(n_params + 1) /
                                static_cast<Real>(n_obs - n_params - 1);
        ic.aicc = ic.aic + correction;
    } else {
        // 小样本未定义, 返回 AIC (调用方应检查 n_obs > n_params + 1)
        ic.aicc = ic.aic;
    }

    return ic;
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
