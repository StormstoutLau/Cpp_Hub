// SOURCE: PHASE5_HFE_SPEC §6.1.2
//   [Hasbrouck 2009] Hasbrouck, "Trading Costs and Returns for U.S. Equities"
//   [Lee-Ready 1991] Lee, Ready, "Inferring trade direction from intraday data"
// R对照: liquidityMeasures.R L231 (getLiquidityMeasures) + L346 (getTradeDirection)
// R 源码实测: 2026-08-03 (hf_src/highfrequency/R/liquidityMeasures.R)
#pragma once

#include <vector>
#include <optional>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// LiquidityMeasures — 23 种流动性度量 (R getLiquidityMeasures L268-312)
// =============================================================================
struct LiquidityMeasures {
    std::vector<double> effectiveSpread;                 // 1.  2*D*(P-mid)
    std::vector<double> realizedSpread;                  // 2.  2*D*(P-mid[t+win])
    std::vector<double> valueTrade;                      // 3.  SIZE*P
    std::vector<double> signedValueTrade;                // 4.  D*SIZE*P
    std::vector<double> depthImbalanceDifference;        // 5.  D*(OFRSIZ-BIDSIZ)/(OFRSIZ+BIDSIZ)
    std::vector<double> depthImbalanceRatio;             // 6.  (D*OFRSIZ/BIDSIZ)^D  [D3: D=-1 负倒数]
    std::vector<double> proportionalEffectiveSpread;     // 7.  effectiveSpread/mid
    std::vector<double> proportionalRealizedSpread;      // 8.  realizedSpread/mid
    std::vector<double> priceImpact;                     // 9.  (effectiveSpread-realizedSpread)/2
    std::vector<double> proportionalPriceImpact;         // 10. priceImpact/mid
    std::vector<double> halfTradedSpread;                // 11. D*(P-mid)
    std::vector<double> proportionalHalfTradedSpread;    // 12. halfTradedSpread/mid
    std::vector<double> squaredLogReturn;                // 13. (log(P)-log(P[t-1]))^2  [首元素 NaN]
    std::vector<double> absLogReturn;                    // 14. |log(P)-log(P[t-1])|     [首元素 NaN]
    std::vector<double> quotedSpread;                    // 15. OFR-BID
    std::vector<double> proportionalQuotedSpread;        // 16. quotedSpread/mid
    std::vector<double> logQuotedSpread;                 // 17. log(OFR/BID)
    std::vector<double> logQuotedSize;                   // 18. log(OFRSIZ)+log(BIDSIZ)
    std::vector<double> quotedSlope;                     // 19. quotedSpread/logQuotedSize
    std::vector<double> logQSlope;                       // 20. logQuotedSpread/logQuotedSize
    std::vector<double> midQuoteSquaredReturn;           // 21. (log(mid)-log(mid[t-1]))^2  [首元素 NaN]
    std::vector<double> midQuoteAbsReturn;               // 22. |log(mid)-log(mid[t-1])|     [首元素 NaN]
    std::vector<double> signedTradeSize;                 // 23. D*SIZE
};

// =============================================================================
// get_trade_direction — R getTradeDirection (L346-380)
// =============================================================================
//
// 算法 (R 源码实测, 排幻觉 D1):
//   1. midpoints = (bid + ofr) / 2
//   2. rets = diff(price)  (长度 n-1)
//   3. tick rule: c(TRUE, fifelse(rets>0,TRUE,fifelse(rets<0,FALSE,NA))) * 2 - 1
//      - 首元素 = TRUE → 1 (buy), **非 0**
//      - rets>0 → 1, rets<0 → -1, rets==0 → NA
//   4. LOCF: NA → 上一个非 NA 值 (首元素恒为 1, 故无 "无前值" 情况)
//   5. midpoint 覆盖 (在 LOCF 之后):
//      - price < mid → -1 (sell)
//      - price > mid → 1 (buy)
//      - price == mid → 保留 tick rule 结果
//
// 返回: 1 (buy) 或 -1 (sell), 长度 n
// 异常: 空输入或 n < 2 抛 invalid_argument
// =============================================================================
inline std::vector<int> get_trade_direction(
    const std::vector<double>& price,
    const std::vector<double>& bid,
    const std::vector<double>& ofr) {
    const Size n = price.size();
    if (n == 0) {
        throw std::invalid_argument("get_trade_direction: empty input");
    }
    if (n != bid.size() || n != ofr.size()) {
        throw std::invalid_argument(
            "get_trade_direction: price, bid, ofr must have equal length");
    }

    // midpoints
    std::vector<double> midpoints(n);
    for (Size i = 0; i < n; ++i) {
        midpoints[i] = 0.5 * (bid[i] + ofr[i]);
    }

    // tick rule + LOCF
    // R: c(TRUE, fifelse(rets>0,TRUE,fifelse(rets<0,FALSE,NA))) * 2 - 1, then nafill locf
    std::vector<int> tick_dir(n);
    tick_dir[0] = 1;  // 首元素 = TRUE → 1 (buy), 非 0
    int last_valid = 1;
    for (Size i = 1; i < n; ++i) {
        double ret = price[i] - price[i - 1];  // diff(price)
        int dir;
        if (ret > 0.0) {
            dir = 1;   // uptick → buy
        } else if (ret < 0.0) {
            dir = -1;  // downtick → sell
        } else {
            dir = last_valid;  // zero-tick → LOCF
        }
        tick_dir[i] = dir;
        last_valid = dir;
    }

    // midpoint 覆盖 (在 LOCF 之后)
    // R: fifelse(price<midpoints, -1, fifelse(price>midpoints, 1, buys))
    std::vector<int> direction(n);
    for (Size i = 0; i < n; ++i) {
        if (price[i] < midpoints[i]) {
            direction[i] = -1;  // price < mid → sell (覆盖)
        } else if (price[i] > midpoints[i]) {
            direction[i] = 1;   // price > mid → buy (覆盖)
        } else {
            direction[i] = tick_dir[i];  // price == mid → 保留 tick rule
        }
    }
    return direction;
}

// =============================================================================
// get_liquidity_measures — R getLiquidityMeasures (L231-319)
// =============================================================================
//
// 输入:
//   price, bid, ofr, size, ofrsiz, bidsiz — 等长 (n) 的 TAQ 数据
//   direction — 可选, 用户指定方向 (仅 ±1); 若 nullopt 则调用 get_trade_direction
//   win — realizedSpread 的 lead shift 窗口 (默认 300)
//
// 输出: LiquidityMeasures 结构体, 23 种度量
//
// 异常:
//   - 输入长度不一致: invalid_argument
//   - 用户 direction 含非 ±1 值: invalid_argument
//   - 空输入: invalid_argument
//
// 排幻觉 (R 源码实测 2026-08-03):
//   D1: direction 推断为 tick rule + midpoint 覆盖混合, 首元素 = 1 (buy)
//   D2: realizedSpread 用 lead shift mid[t+win], 越界为 NaN
//   D3: depthImbalanceRatio = (D*OFRSIZ/BIDSIZ)^D, D=-1 时结果为负倒数
//
// NaN 传播规则:
//   - realizedSpread: t+win >= n 时 NaN
//   - proportionalRealizedSpread, priceImpact, proportionalPriceImpact:
//     依赖 realizedSpread, NaN 传播
//   - squaredLogReturn, absLogReturn: 首元素 NaN (无前值)
//   - midQuoteSquaredReturn, midQuoteAbsReturn: 首元素 NaN
// =============================================================================
inline LiquidityMeasures get_liquidity_measures(
    const std::vector<double>& price, const std::vector<double>& bid,
    const std::vector<double>& ofr, const std::vector<double>& size,
    const std::vector<double>& ofrsiz, const std::vector<double>& bidsiz,
    const std::optional<std::vector<int>>& direction = std::nullopt,
    int win = 300) {
    const Size n = price.size();
    if (n == 0) {
        throw std::invalid_argument("get_liquidity_measures: empty input");
    }
    if (n != bid.size() || n != ofr.size() || n != size.size() ||
        n != ofrsiz.size() || n != bidsiz.size()) {
        throw std::invalid_argument(
            "get_liquidity_measures: all input vectors must have equal length");
    }
    if (win < 0) {
        throw std::invalid_argument(
            "get_liquidity_measures: win must be >= 0");
    }

    // 方向推断或用户输入
    std::vector<int> dir;
    if (direction.has_value()) {
        dir = direction.value();
        if (dir.size() != n) {
            throw std::invalid_argument(
                "get_liquidity_measures: direction size mismatch");
        }
        for (Size i = 0; i < n; ++i) {
            if (dir[i] != 1 && dir[i] != -1) {
                throw std::invalid_argument(
                    "get_liquidity_measures: DIRECTION must be -1 or 1 "
                    "(invalid value at index " + std::to_string(i) + ")");
            }
        }
    } else {
        dir = get_trade_direction(price, bid, ofr);
    }

    // midpoints
    std::vector<double> midpoints(n);
    for (Size i = 0; i < n; ++i) {
        midpoints[i] = 0.5 * (bid[i] + ofr[i]);
    }

    // NaN 常量
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

    LiquidityMeasures m;
    m.effectiveSpread.resize(n);
    m.realizedSpread.resize(n);
    m.valueTrade.resize(n);
    m.signedValueTrade.resize(n);
    m.depthImbalanceDifference.resize(n);
    m.depthImbalanceRatio.resize(n);
    m.proportionalEffectiveSpread.resize(n);
    m.proportionalRealizedSpread.resize(n);
    m.priceImpact.resize(n);
    m.proportionalPriceImpact.resize(n);
    m.halfTradedSpread.resize(n);
    m.proportionalHalfTradedSpread.resize(n);
    m.squaredLogReturn.resize(n);
    m.absLogReturn.resize(n);
    m.quotedSpread.resize(n);
    m.proportionalQuotedSpread.resize(n);
    m.logQuotedSpread.resize(n);
    m.logQuotedSize.resize(n);
    m.quotedSlope.resize(n);
    m.logQSlope.resize(n);
    m.midQuoteSquaredReturn.resize(n);
    m.midQuoteAbsReturn.resize(n);
    m.signedTradeSize.resize(n);

    for (Size t = 0; t < n; ++t) {
        const double D = static_cast<double>(dir[t]);
        const double P = price[t];
        const double mid = midpoints[t];
        const double S = size[t];
        const double osiz = ofrsiz[t];
        const double bsiz = bidsiz[t];

        // 1. effectiveSpread = 2*D*(P-mid)
        m.effectiveSpread[t] = 2.0 * D * (P - mid);

        // 2. realizedSpread = 2*D*(P - mid[t+win])  [D2: lead shift, 越界 NaN]
        if (t + static_cast<Size>(win) < n) {
            m.realizedSpread[t] = 2.0 * D * (P - midpoints[t + win]);
        } else {
            m.realizedSpread[t] = NaN;
        }

        // 3. valueTrade = SIZE*P
        m.valueTrade[t] = S * P;

        // 4. signedValueTrade = D*valueTrade
        m.signedValueTrade[t] = D * m.valueTrade[t];

        // 5. depthImbalanceDifference = D*(OFRSIZ-BIDSIZ)/(OFRSIZ+BIDSIZ)
        m.depthImbalanceDifference[t] = D * (osiz - bsiz) / (osiz + bsiz);

        // 6. depthImbalanceRatio = (D*OFRSIZ/BIDSIZ)^D  [D3: D=-1 负倒数]
        // R: (direction * OFRSIZ / BIDSIZ) ^ direction
        // C++ std::pow 支持负底数 + 整数指数 (C++11 起)
        m.depthImbalanceRatio[t] = std::pow(D * osiz / bsiz, D);

        // 7. proportionalEffectiveSpread = effectiveSpread/mid
        m.proportionalEffectiveSpread[t] = m.effectiveSpread[t] / mid;

        // 8. proportionalRealizedSpread = realizedSpread/mid  [NaN 传播]
        m.proportionalRealizedSpread[t] = m.realizedSpread[t] / mid;

        // 9. priceImpact = (effectiveSpread - realizedSpread)/2  [NaN 传播]
        m.priceImpact[t] = (m.effectiveSpread[t] - m.realizedSpread[t]) / 2.0;

        // 10. proportionalPriceImpact = priceImpact/mid  [NaN 传播]
        m.proportionalPriceImpact[t] = m.priceImpact[t] / mid;

        // 11. halfTradedSpread = D*(P-mid)
        m.halfTradedSpread[t] = D * (P - mid);

        // 12. proportionalHalfTradedSpread = halfTradedSpread/mid
        m.proportionalHalfTradedSpread[t] = m.halfTradedSpread[t] / mid;

        // 13. squaredLogReturn = (log(P)-log(P[t-1]))^2  [首元素 NaN]
        if (t >= 1) {
            double lr = std::log(P) - std::log(price[t - 1]);
            m.squaredLogReturn[t] = lr * lr;
        } else {
            m.squaredLogReturn[t] = NaN;
        }

        // 14. absLogReturn = |log(P)-log(P[t-1])|  [首元素 NaN]
        if (t >= 1) {
            m.absLogReturn[t] = std::fabs(std::log(P) - std::log(price[t - 1]));
        } else {
            m.absLogReturn[t] = NaN;
        }

        // 15. quotedSpread = OFR-BID
        m.quotedSpread[t] = ofr[t] - bid[t];

        // 16. proportionalQuotedSpread = quotedSpread/mid
        m.proportionalQuotedSpread[t] = m.quotedSpread[t] / mid;

        // 17. logQuotedSpread = log(OFR/BID)
        m.logQuotedSpread[t] = std::log(ofr[t] / bid[t]);

        // 18. logQuotedSize = log(OFRSIZ)+log(BIDSIZ)
        m.logQuotedSize[t] = std::log(osiz) + std::log(bsiz);

        // 19. quotedSlope = quotedSpread/logQuotedSize
        m.quotedSlope[t] = m.quotedSpread[t] / m.logQuotedSize[t];

        // 20. logQSlope = logQuotedSpread/logQuotedSize
        m.logQSlope[t] = m.logQuotedSpread[t] / m.logQuotedSize[t];

        // 21. midQuoteSquaredReturn = (log(mid)-log(mid[t-1]))^2  [首元素 NaN]
        if (t >= 1) {
            double mlr = std::log(mid) - std::log(midpoints[t - 1]);
            m.midQuoteSquaredReturn[t] = mlr * mlr;
        } else {
            m.midQuoteSquaredReturn[t] = NaN;
        }

        // 22. midQuoteAbsReturn = |log(mid)-log(mid[t-1])|  [首元素 NaN]
        if (t >= 1) {
            m.midQuoteAbsReturn[t] =
                std::fabs(std::log(mid) - std::log(midpoints[t - 1]));
        } else {
            m.midQuoteAbsReturn[t] = NaN;
        }

        // 23. signedTradeSize = D*SIZE
        m.signedTradeSize[t] = D * S;
    }

    return m;
}

} // namespace hfecon
} // namespace v1
} // namespace cpphub
