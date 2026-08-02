// SOURCE: PHASE5_HFE_SPEC §4.3
//   [BNS 2008] Barndorff-Nielsen, Hansen, Lunde, Shephard,
//              Econometrica 76(6), 1481-1536, doi:10.1111/j.1468-0262.2008.00837.x
//   [H-L 2006] Hansen & Lunde, JBES 24(2), 127-161, doi:10.1198/073500106000000072
// R 对照: highfrequency 1.0.3 无独立导出函数, rKernelCov 内部使用 (BNS 2008 §4.4)
#pragma once

#include <vector>
#include <stdexcept>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// NoiseVarianceEstimator: 微结构噪声方差 ω² 估计
// =============================================================================
// 算法 (BNS 2008 eq. 40, H-L 2006 §3):
//   ω² = (1/(2n)) * Σ_{i=1}^{n} r_i²
//
// 其中 n = 收益率数, r_i = 日内对数收益率
//
// 性质: 在噪声独立同分布假设下, ω² 是噪声方差的一致估计量
// 用途: rKernelCov 内部 bandwidth 选择 (BNS 2008 §4.5)
//
// 注意: 此模块为 BNS 2008 公式实现, 非 R 对标 (R 无独立导出函数)
//       RealizedKernel::estimate 不调用此模块
// =============================================================================

struct NoiseVarianceResult {
    Real omega2;               // 噪声方差估计 ω²
    Real integrated_variance;  // IV 估计 = RV - ω²
    Size n_obs;                // 观测数
};

class NoiseVarianceEstimator {
public:
    // 输入: 日内对数收益率序列
    // 异常: n < 2 抛 invalid_argument
    static NoiseVarianceResult estimate(const std::vector<Real>& log_returns) {
        const Size n = log_returns.size();
        if (n < 2) {
            throw std::invalid_argument(
                "NoiseVarianceEstimator::estimate requires n >= 2 returns");
        }

        Real sum_r2 = 0.0;
        for (Size i = 0; i < n; ++i) {
            sum_r2 += log_returns[i] * log_returns[i];
        }

        NoiseVarianceResult result;
        result.omega2 = sum_r2 / (2.0 * static_cast<Real>(n));
        result.integrated_variance = sum_r2 - result.omega2;  // IV = RV - ω²
        result.n_obs = n;
        return result;
    }
};

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
