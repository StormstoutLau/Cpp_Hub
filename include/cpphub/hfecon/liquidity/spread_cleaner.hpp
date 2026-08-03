// SOURCE: PHASE5_HFE_SPEC §6.1.1
//   [BKS 2022] Boudt, Kleen, Sjoerup, JSS 104(8), 1-36, doi:10.18637/jss.v104.i08
// R对照: rmLargeSpread / rmNegativeSpread / spreadPrices
//   rmLargeSpread:  dataHandling.R L1617 — 每日 SPREAD 中位数筛选
//   rmNegativeSpread: dataHandling.R L1670 — 保留 OFR > BID
//   spreadPrices:   dataHandling.R L3019 — 长格式 (DT,SYMBOL,PRICE) → 宽格式
// R 源码实测: 2026-08-03
#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <cmath>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// 数据结构: 报价观测 (用于 spread cleaner)
// =============================================================================
struct QuoteObservation {
    int64_t ts;      // 纳秒时间戳
    Real bid;        // BID 价格
    Real ask;        // OFR 价格
};

// 宽格式行: DT + 各 SYMBOL 价格
struct SpreadPriceRow {
    int64_t dt;
    std::vector<Real> prices;  // 按 symbols 顺序排列, 缺失为 NaN
};

// =============================================================================
// rm_negative_spread — 保留 OFR > BID (严格大于)
// R: qData[OFR > BID]  (dataHandling.R L1686)
// =============================================================================
inline std::vector<QuoteObservation> rm_negative_spread(
    const std::vector<QuoteObservation>& qData) {
    std::vector<QuoteObservation> result;
    result.reserve(qData.size());
    for (const auto& q : qData) {
        if (q.ask > q.bid) {  // R: OFR > BID (严格)
            result.push_back(q);
        }
    }
    return result;
}

// =============================================================================
// rm_large_spread — 每日 SPREAD 中位数筛选
// R (dataHandling.R L1645-1647):
//   qData[, DATE := as.Date(DT, tz = tz)]
//   qData[, SPREAD := OFR - BID]
//   qData[, SPREAD_MEDIAN := median(SPREAD), by = "DATE"]
//   qData[SPREAD < (SPREAD_MEDIAN * maxi)]
//
// 日期分组: 按 ts / 86400e9 取整作为 DATE
// 中位数: std::nth_element (偶数取中间两数均值, R median 默认 type=7)
// 筛选: SPREAD < SPREAD_MEDIAN * maxi (严格小于)
// =============================================================================
inline std::vector<QuoteObservation> rm_large_spread(
    const std::vector<QuoteObservation>& qData, double maxi = 50.0) {
    if (qData.empty()) return {};
    if (maxi <= 0.0) {
        throw std::invalid_argument("rm_large_spread: maxi must be > 0");
    }

    // 按日期分组 (epoch 天数)
    const int64_t one_day_ns = 86400LL * 1000000000LL;
    std::map<int64_t, std::vector<Size>> day_groups;  // day_idx -> [original indices]
    for (Size i = 0; i < qData.size(); ++i) {
        int64_t day_idx = qData[i].ts / one_day_ns;
        if (qData[i].ts < 0 && qData[i].ts % one_day_ns != 0) --day_idx;
        day_groups[day_idx].push_back(i);
    }

    // 每日计算中位数并筛选
    std::vector<QuoteObservation> result;
    result.reserve(qData.size());
    for (const auto& [day_idx, indices] : day_groups) {
        // 计算每日 SPREAD
        std::vector<Real> spreads;
        spreads.reserve(indices.size());
        for (Size idx : indices) {
            spreads.push_back(qData[idx].ask - qData[idx].bid);
        }
        // 中位数 (R type=7: n 个元素, 取 (n-1)/2 和 n/2 的均值, 即标准中位数)
        std::vector<Real> sorted_spreads = spreads;
        std::sort(sorted_spreads.begin(), sorted_spreads.end());
        Real median;
        Size n = sorted_spreads.size();
        if (n % 2 == 1) {
            median = sorted_spreads[n / 2];
        } else {
            median = 0.5 * (sorted_spreads[n / 2 - 1] + sorted_spreads[n / 2]);
        }
        // 筛选: SPREAD < SPREAD_MEDIAN * maxi
        Real threshold = median * maxi;
        for (Size idx : indices) {
            Real spread = qData[idx].ask - qData[idx].bid;
            if (spread < threshold) {  // R: SPREAD < (SPREAD_MEDIAN * maxi)
                result.push_back(qData[idx]);
            }
        }
    }
    return result;
}

// =============================================================================
// spread_prices — 长格式 (DT, SYMBOL, PRICE) → 宽格式
// R (dataHandling.R L3040-3047):
//   splitted <- split(data[,list(DT, PRICE, SYMBOL)], by = 'SYMBOL')
//   collected <- Reduce(merge, lapply(splitted, ...))  # outer join on DT
//
// 实现:
//   1. 收集所有唯一 DT (排序)
//   2. 收集所有唯一 SYMBOL (按首次出现顺序)
//   3. 构建行: 每个 DT 一行, prices[i] = 对应 SYMBOL 的价格 (缺失为 NaN)
// =============================================================================
inline std::vector<SpreadPriceRow> spread_prices(
    const std::vector<int64_t>& dt,
    const std::vector<std::string>& symbols,
    const std::vector<Real>& prices) {
    if (dt.size() != symbols.size() || dt.size() != prices.size()) {
        throw std::invalid_argument(
            "spread_prices: dt, symbols, prices must have same length");
    }
    if (dt.empty()) return {};

    // 收集唯一 SYMBOL (按首次出现顺序)
    std::vector<std::string> uniq_symbols;
    std::map<std::string, Size> sym_col;  // symbol -> column index
    for (const auto& s : symbols) {
        if (sym_col.find(s) == sym_col.end()) {
            sym_col[s] = uniq_symbols.size();
            uniq_symbols.push_back(s);
        }
    }
    Size n_cols = uniq_symbols.size();

    // 收集唯一 DT (排序)
    std::set<int64_t> uniq_dt_set(dt.begin(), dt.end());
    std::vector<int64_t> uniq_dt(uniq_dt_set.begin(), uniq_dt_set.end());
    std::sort(uniq_dt.begin(), uniq_dt.end());

    // 构建 DT -> row_idx 映射
    std::map<int64_t, Size> dt_row;
    for (Size i = 0; i < uniq_dt.size(); ++i) {
        dt_row[uniq_dt[i]] = i;
    }

    // 初始化结果矩阵 (NaN)
    std::vector<std::vector<Real>> matrix(uniq_dt.size(),
                                           std::vector<Real>(n_cols,
                                               std::numeric_limits<Real>::quiet_NaN()));

    // 填充: (DT, SYMBOL, PRICE) -> matrix[row][col]
    for (Size i = 0; i < dt.size(); ++i) {
        Size row = dt_row[dt[i]];
        Size col = sym_col[symbols[i]];
        matrix[row][col] = prices[i];
    }

    // 转换为 SpreadPriceRow
    std::vector<SpreadPriceRow> result;
    result.reserve(uniq_dt.size());
    for (Size i = 0; i < uniq_dt.size(); ++i) {
        result.push_back({uniq_dt[i], std::move(matrix[i])});
    }
    return result;
}

} // namespace hfecon
} // namespace v1
} // namespace cpphub
