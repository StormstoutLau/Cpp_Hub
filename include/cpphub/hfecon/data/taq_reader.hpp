// SOURCE: PHASE5_HFE_SPEC §3.1
//   [BKS 2022] Boudt, Kleen, Sjoerup, JSS 104(8), 1-36, doi:10.18637/jss.v104.i08
// R对照: aggregatePrice / makeReturns
// ITCH 5.0: NASDAQ TotalView 官方协议 (v1.4.0 推迟, 仅实现 CSV)
#pragma once

// MSVC CRT 安全警告抑制 (sscanf 用于数值解析, 无缓冲区溢出风险)
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable: 4996)  // 'sscanf': This function or variable may be unsafe
#endif

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <ctime>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// TAQ 数据结构 (spec §3.1)
// =============================================================================
// Timestamp: 纳秒时间戳 (int64, core/types.hpp 未定义, 此处用 int64_t)
// Trade/Quote: 对应 R highfrequency sampleTData / sampleQData 结构
// =============================================================================

using Timestamp = int64_t;  // 纳秒, 1970-01-01 UTC

struct Trade {
    Timestamp ts;       // 纳秒时间戳
    Real price;
    Real size;
    std::string mpid;   // 可选, 做市商 ID (v1.4.0 暂不解析)
};

struct Quote {
    Timestamp ts;
    Real bid_price;
    Real ask_price;
    Real bid_size;
    Real ask_size;
};

// =============================================================================
// TaqReader: CSV 读取 + 时间聚合 (spec §3.1)
// =============================================================================
// CSV 格式 (highfrequency 兼容):
//   DT,PRICE,SIZE,BID,ASK,BIDSIZE,ASKSIZE
//   2024-01-02 09:30:00.123456,100.05,100,100.04,100.06,500,300
//
// aggregatePrice 语义 (R highfrequency):
//   - 按 alignBy ("seconds"/"minutes"/"ticks") 对齐到 alignPeriod 桶
//   - 每桶取最后成交价 (last-tick aggregation)
//   - 时间戳为桶的结束时刻
//   - market_open/market_close 限定交易时段 (格式 "HH:MM:SS")
// =============================================================================

class TaqReader {
public:
    // CSV 读取 (highfrequency 兼容格式, 至少含 DT,PRICE,SIZE 列)
    // 抛 invalid_argument: 文件无法打开或格式错误
    static std::vector<Trade> read_trades_csv(const std::string& path) {
        std::ifstream ifs(path);
        if (!ifs) {
            throw std::invalid_argument("read_trades_csv: cannot open " + path);
        }
        std::string line;
        std::vector<Trade> trades;
        // 解析表头, 定位列索引
        std::getline(ifs, line);
        auto hdr = parse_csv_line(line);
        Size idx_dt = SIZE_MAX, idx_price = SIZE_MAX, idx_size = SIZE_MAX;
        for (Size i = 0; i < hdr.size(); ++i) {
            std::string h = to_lower(trim(hdr[i]));
            if (h == "dt")               idx_dt = i;
            else if (h == "price")       idx_price = i;
            else if (h == "size")        idx_size = i;
        }
        if (idx_dt == SIZE_MAX || idx_price == SIZE_MAX) {
            throw std::invalid_argument(
                "read_trades_csv: missing DT or PRICE column in " + path);
        }

        while (std::getline(ifs, line)) {
            if (line.empty()) continue;
            auto fields = parse_csv_line(line);
            if (fields.size() <= std::max(idx_dt, idx_price)) continue;
            Trade t;
            t.ts = parse_timestamp(trim(fields[idx_dt]));
            t.price = std::stod(trim(fields[idx_price]));
            t.size = (idx_size != SIZE_MAX && idx_size < fields.size())
                         ? std::stod(trim(fields[idx_size])) : 0.0;
            trades.push_back(std::move(t));
        }
        return trades;
    }

    // CSV 读取报价 (至少含 DT,BID,ASK 列; BIDSIZE/ASKSIZE 可选)
    static std::vector<Quote> read_quotes_csv(const std::string& path) {
        std::ifstream ifs(path);
        if (!ifs) {
            throw std::invalid_argument("read_quotes_csv: cannot open " + path);
        }
        std::string line;
        std::vector<Quote> quotes;
        std::getline(ifs, line);
        auto hdr = parse_csv_line(line);
        Size idx_dt = SIZE_MAX, idx_bid = SIZE_MAX, idx_ask = SIZE_MAX;
        Size idx_bsize = SIZE_MAX, idx_asize = SIZE_MAX;
        for (Size i = 0; i < hdr.size(); ++i) {
            std::string h = to_lower(trim(hdr[i]));
            if (h == "dt")           idx_dt = i;
            else if (h == "bid")     idx_bid = i;
            else if (h == "ask")     idx_ask = i;
            else if (h == "bidsize") idx_bsize = i;
            else if (h == "asksize") idx_asize = i;
        }
        if (idx_dt == SIZE_MAX || idx_bid == SIZE_MAX || idx_ask == SIZE_MAX) {
            throw std::invalid_argument(
                "read_quotes_csv: missing DT/BID/ASK column in " + path);
        }
        while (std::getline(ifs, line)) {
            if (line.empty()) continue;
            auto fields = parse_csv_line(line);
            Quote q;
            q.ts = parse_timestamp(trim(fields[idx_dt]));
            q.bid_price = std::stod(trim(fields[idx_bid]));
            q.ask_price = std::stod(trim(fields[idx_ask]));
            q.bid_size = (idx_bsize != SIZE_MAX && idx_bsize < fields.size())
                             ? std::stod(trim(fields[idx_bsize])) : 0.0;
            q.ask_size = (idx_asize != SIZE_MAX && idx_asize < fields.size())
                             ? std::stod(trim(fields[idx_asize])) : 0.0;
            quotes.push_back(std::move(q));
        }
        return quotes;
    }

    // 时间聚合 (对应 R aggregatePrice, last-tick aggregation)
    // align_by: "seconds" / "minutes" / "ticks"
    // 返回: 聚合后的 Trade 序列 (price = 每桶最后成交价, size = 0)
    // 异常: 空输入或未知 align_by 抛 invalid_argument
    static std::vector<Trade> aggregate_price(
        const std::vector<Trade>& trades,
        Size align_period,
        const std::string& align_by = "minutes",
        const std::string& market_open = "09:30:00",
        const std::string& market_close = "16:00:00",
        bool fill = false) {
        (void)fill;  // v1.4.0 暂不实现 fill
        if (trades.empty()) {
            throw std::invalid_argument("aggregate_price: empty trades");
        }
        if (align_period == 0) {
            throw std::invalid_argument("aggregate_price: align_period must be > 0");
        }

        const Timestamp open_ns = parse_hhmmss(market_open);
        const Timestamp close_ns = parse_hhmmss(market_close);
        const Timestamp bucket_ns = bucket_to_ns(align_by, align_period);

        std::vector<Trade> result;
        if (align_by == "ticks") {
            // 每 align_period 个 tick 取最后一个
            for (Size i = 0; i < trades.size(); i += align_period) {
                Size last = std::min(i + align_period - 1, trades.size() - 1);
                result.push_back(trades[last]);
            }
            return result;
        }

        // 时间桶聚合: 按 (ts - open) / bucket 分桶, 每桶取最后成交价
        // 算法: 单遍扫描, 维护 (bucket_idx, last_price) 向量, 桶切换时新增条目
        std::vector<std::pair<Timestamp, Real>> buckets;  // (bucket_idx, last_price)
        bool in_market = false;
        for (const auto& tr : trades) {
            Timestamp t = ts_to_intraday_ns(tr.ts);
            if (t < open_ns || t > close_ns) continue;
            in_market = true;
            Timestamp bucket_idx = (t - open_ns) / bucket_ns;
            if (buckets.empty() || bucket_idx != buckets.back().first) {
                buckets.emplace_back(bucket_idx, tr.price);
            } else {
                buckets.back().second = tr.price;
            }
        }
        if (!in_market) {
            throw std::invalid_argument(
                "aggregate_price: no trades within market hours");
        }
        for (const auto& b : buckets) {
            Trade t;
            t.ts = open_ns + (b.first + 1) * bucket_ns;  // 桶结束时刻
            t.price = b.second;
            t.size = 0;
            result.push_back(t);
        }
        return result;
    }

    // 收益率生成 (对应 R makeReturns, 数值向量行为)
    // 复用 RealizedMeasuresCalculator::make_returns 避免重复
    static std::vector<Real> make_returns(const std::vector<Trade>& trades) {
        std::vector<Real> prices;
        prices.reserve(trades.size());
        for (const auto& t : trades) prices.push_back(t.price);
        // 委托给 realized_measures.hpp 的实现 (等长, 首0)
        // 为避免循环依赖, 此处内联实现
        const Size n = prices.size();
        if (n == 0) return {};
        std::vector<Real> ret(n, 0.0);
        for (Size i = 1; i < n; ++i) {
            ret[i] = std::log(prices[i] / prices[i - 1]);
        }
        return ret;
    }

    // ITCH 5.0 二进制读取 (v1.4.0 推迟, 抛 not_implemented)
    static std::vector<Trade> read_trades_itch(const std::string& /*path*/) {
        throw std::runtime_error(
            "read_trades_itch: ITCH 5.0 parser deferred to v1.4.1");
    }
    static std::vector<Quote> read_quotes_itch(const std::string& /*path*/,
                                                bool /*reconstruct_book*/ = true) {
        throw std::runtime_error(
            "read_quotes_itch: ITCH 5.0 parser deferred to v1.4.1");
    }

private:
    // CSV 行解析 (不处理引号内逗号, TAQ 数据通常无引号)
    static std::vector<std::string> parse_csv_line(const std::string& line) {
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, ',')) {
            fields.push_back(field);
        }
        return fields;
    }

    static std::string trim(const std::string& s) {
        Size a = 0, b = s.size();
        while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
        while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
        return s.substr(a, b - a);
    }

    static std::string to_lower(const std::string& s) {
        std::string r = s;
        std::transform(r.begin(), r.end(), r.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        return r;
    }

    // 解析 "YYYY-MM-DD HH:MM:SS.ffffff" 为纳秒时间戳 (UTC epoch)
    // 仅取日内时间部分用于聚合 (返回值含完整 epoch ns)
    static Timestamp parse_timestamp(const std::string& s) {
        // 简化: 期望格式 "YYYY-MM-DD HH:MM:SS" 或带 ".ffffff"
        // 使用 std::tm + timegm
        std::tm tm{};
        int year, mon, day, hh, mm, ss;
        char dot;
        double frac = 0.0;
        int n = std::sscanf(s.c_str(), "%d-%d-%d %d:%d:%d",
                            &year, &mon, &day, &hh, &mm, &ss);
        if (n < 6) {
            throw std::invalid_argument(
                "parse_timestamp: bad format '" + s + "'");
        }
        // 尝试解析小数秒
        std::sscanf(s.c_str(), "%*d-%*d-%*d %*d:%*d:%d%c%lf",
                    &ss, &dot, &frac);
        tm.tm_year = year - 1900;
        tm.tm_mon = mon - 1;
        tm.tm_mday = day;
        tm.tm_hour = hh;
        tm.tm_min = mm;
        tm.tm_sec = ss;
        // 跨平台 UTC epoch: MSVC _mkgmtime / POSIX timegm
        std::time_t t;
#if defined(_MSC_VER)
        t = _mkgmtime(&tm);
        if (t == static_cast<std::time_t>(-1)) {
            t = std::mktime(&tm);  // 回退本地时区
        }
#else
        // POSIX timegm (GCC/Clang); 若不可用回退 mktime
        extern time_t timegm(struct tm*);
        t = timegm(&tm);
        if (t == static_cast<std::time_t>(-1)) {
            t = std::mktime(&tm);
        }
#endif
        Timestamp ns = static_cast<Timestamp>(t) * 1000000000LL;
        ns += static_cast<Timestamp>(frac * 1e9);
        return ns;
    }

    // 提取日内纳秒 (从 epoch ns 中提取 HH:MM:SS.ns)
    static Timestamp ts_to_intraday_ns(Timestamp epoch_ns) {
        // 取模一天的纳秒数
        const Timestamp one_day = 86400LL * 1000000000LL;
        Timestamp t = epoch_ns % one_day;
        if (t < 0) t += one_day;
        return t;
    }

    // "HH:MM:SS" -> 日内纳秒
    static Timestamp parse_hhmmss(const std::string& s) {
        int hh, mm, ss;
        if (std::sscanf(s.c_str(), "%d:%d:%d", &hh, &mm, &ss) != 3) {
            throw std::invalid_argument(
                "parse_hhmmss: bad format '" + s + "'");
        }
        return (static_cast<Timestamp>(hh) * 3600 +
                static_cast<Timestamp>(mm) * 60 +
                static_cast<Timestamp>(ss)) * 1000000000LL;
    }

    // align_by + align_period -> 桶大小 (纳秒)
    static Timestamp bucket_to_ns(const std::string& align_by, Size align_period) {
        Timestamp unit;
        if (align_by == "seconds" || align_by == "secs") {
            unit = 1000000000LL;
        } else if (align_by == "minutes" || align_by == "mins") {
            unit = 60LL * 1000000000LL;
        } else if (align_by == "hours") {
            unit = 3600LL * 1000000000LL;
        } else if (align_by == "ticks") {
            return 1;  // ticks 模式不使用 ns 桶
        } else {
            throw std::invalid_argument(
                "bucket_to_ns: unknown align_by '" + align_by + "'");
        }
        return unit * static_cast<Timestamp>(align_period);
    }
};

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub

#if defined(_MSC_VER)
#pragma warning(pop)  // 恢复 C4996 警告
#endif
