// SOURCE: PHASE5_HFE_SPEC §6.1.3
//   [Amihud 2002] Yakov Amihud, "Illiquidity and stock returns: cross-section
//                  and time-series effects", J. Financial Markets 52(5), 2002
// R对照: highfrequency 无直接实现; 公式为 Amihud (2002) 标准定义
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// amihud_illiquidity — Amihud (2002) 非流动性度量
// =============================================================================
//
// 公式:
//   ILLIQ = (1/N) * sum_{t=1}^{N} ( |r_t| / DVOL_t )
//
// 其中:
//   r_t     = 日对数收益率 (输入序列)
//   DVOL_t  = 日美元成交额 (dollar volume, 通常 = price * volume)
//   N       = 观测数 (输入序列长度)
//
// 性质:
//   - ILLIQ 越大, 流动性越差 (单位成交额引起的价格变动越大)
//   - 单位: 1/货币 (如 1/USD), 通常乘以 1e6 以 (basis points per million USD) 表示
//   - 对零成交额敏感 (除零), C++ 实现显式抛异常
//
// 异常:
//   - 空输入序列: invalid_argument
//   - returns.size() != dollarVolume.size(): invalid_argument
//   - 任一 dollarVolume[t] <= 0: invalid_argument (避免除零和负成交额)
//
// 容差: 1e-12 (纯算术, 无浮点累积)
// =============================================================================
inline double amihud_illiquidity(
    const std::vector<double>& dailyReturns,
    const std::vector<double>& dailyDollarVolume) {
    const Size n = dailyReturns.size();
    if (n == 0) {
        throw std::invalid_argument(
            "amihud_illiquidity: empty input sequences");
    }
    if (n != dailyDollarVolume.size()) {
        throw std::invalid_argument(
            "amihud_illiquidity: dailyReturns and dailyDollarVolume "
            "must have the same length");
    }

    double sum = 0.0;
    for (Size t = 0; t < n; ++t) {
        if (dailyDollarVolume[t] <= 0.0) {
            throw std::invalid_argument(
                "amihud_illiquidity: dailyDollarVolume must be > 0 "
                "(zero or negative dollar volume at index " +
                std::to_string(t) + ")");
        }
        sum += std::fabs(dailyReturns[t]) / dailyDollarVolume[t];
    }

    return sum / static_cast<double>(n);
}

} // namespace hfecon
} // namespace v1
} // namespace cpphub
