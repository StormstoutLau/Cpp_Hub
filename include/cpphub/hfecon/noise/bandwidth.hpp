// SOURCE: PHASE5_HFE_SPEC §4.4
//   [BNS 2008] Barndorff-Nielsen, Hansen, Lunde, Shephard,
//              Econometrica 76(6), 1481-1536, doi:10.1111/j.1468-0262.2008.00837.x
// R 对照: rKernelCov 内部 bandwidth 选择 (BNS 2008 §4.5)
#pragma once

#include <cmath>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/measures/kernels.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// optimal_bandwidth: BNS 2008 最优 bandwidth 选择
// =============================================================================
// 算法 (BNS 2008 eq. 51):
//   H* = c × ξ^(4/5) × (ω²/IV)^(2/5) × n^(3/5)
//
// 其中:
//   c  = 5.74 (Bartlett 核的最优常数, BNS 2008 Table 4)
//        其他核的 c 值不同, 但 R rKernelCov 默认使用 Bartlett
//   ξ² = IV + 2ω²  (噪声修正的 IV)
//   ω² = 噪声方差 (NoiseVarianceEstimator)
//   IV  = RV - ω²  (积分方差估计)
//   n   = 观测数
//
// round 到最近整数, 下界 H ≥ 1
//
// 注意: 此模块为 BNS 2008 公式实现, 非 R 对标
//       RealizedKernel::estimate 不调用此模块 (R rKernelCov 接受用户提供 kernelParam)
// =============================================================================

inline Size optimal_bandwidth(Real omega2, Real integrated_variance,
                               Size n_obs, KernelType type) {
    if (omega2 <= 0.0) {
        throw std::invalid_argument(
            "optimal_bandwidth: omega2 must be positive");
    }
    if (integrated_variance <= 0.0) {
        throw std::invalid_argument(
            "optimal_bandwidth: integrated_variance must be positive");
    }
    if (n_obs == 0) {
        throw std::invalid_argument(
            "optimal_bandwidth: n_obs must be positive");
    }

    // 最优常数 c (BNS 2008 Table 4)
    // Bartlett: c = 5.74 (R 默认)
    // 其他核的理论 c 值不同, 但 BNS 2008 仅给出 Bartlett 的精确值
    // 对于非 Bartlett 核, 使用 Bartlett 的 c 作为近似 (BNS 2008 §4.5 建议)
    Real c = 5.74;  // Bartlett 默认

    // ξ² = IV + 2ω²
    Real xi2 = integrated_variance + 2.0 * omega2;
    Real xi = std::sqrt(xi2);

    // H* = c × ξ^(4/5) × (ω²/IV)^(2/5) × n^(3/5)
    Real n_real = static_cast<Real>(n_obs);
    Real H_star = c *
                  std::pow(xi, 4.0 / 5.0) *
                  std::pow(omega2 / integrated_variance, 2.0 / 5.0) *
                  std::pow(n_real, 3.0 / 5.0);

    // round 到最近整数, 下界 H ≥ 1
    Size H = static_cast<Size>(std::round(H_star));
    if (H < 1) H = 1;
    return H;
}

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
