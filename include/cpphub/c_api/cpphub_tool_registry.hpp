// 工具注册表：7 个诊断类 MCP 工具的 name / description / inputSchema / handler。
//
// 依据 spec/mcp-server/IMPLEMENTATION.md §3.3、DESIGN.md §4.1（F7 订正后入参）、
// §6.1（逐方法信封）、§6.2（ToolResult 容器 + Conclusion）。
//
// 设计要点（ADR-015 / DESIGN §8 不变式）：
//   - 单一登记点：新增工具只改本文件 → 所有前端（B1 Python 壳 / B2 原生 exe）自动获得
//   - handler 为纯函数：仅依赖本次请求参数，无服务端会话态（不变式 1）
//   - 只读：不写文件、不发网络（不变式 4）
//   - description 内嵌文献陷阱注释（GR6 方向 / CI2 p 与 cv 分列 / U12 KPSS 方向等）
//   - handler 抛异常不逃逸：在 handler 内逐测试捕获，转为该行的 notes + statistic=null
#pragma once

#include "cpphub/c_api/detail/json_min.hpp"

#include "cpphub/econometrics/inference/detail/test_result_base.hpp"
#include "cpphub/timeseries/arima/granger_test.hpp"
#include "cpphub/timeseries/cointegration/engle_granger.hpp"
#include "cpphub/timeseries/cointegration/johansen_test.hpp"
#include "cpphub/timeseries/cointegration/vecm_model.hpp"
#include "cpphub/timeseries/garch/egarch_model.hpp"
#include "cpphub/timeseries/garch/garch_diagnostics.hpp"
#include "cpphub/timeseries/garch/garch_m_model.hpp"
#include "cpphub/timeseries/garch/garch_model.hpp"
#include "cpphub/timeseries/garch/gjr_garch_model.hpp"
#include "cpphub/timeseries/unit_root/adf_test.hpp"
#include "cpphub/timeseries/unit_root/df_gls_test.hpp"
#include "cpphub/timeseries/unit_root/kpss_test.hpp"
#include "cpphub/timeseries/unit_root/ng_perron_test.hpp"
#include "cpphub/timeseries/unit_root/pp_test.hpp"
#include "cpphub/timeseries/unit_root/variance_ratio_test.hpp"
#include "cpphub/timeseries/unit_root/zivot_andrews_test.hpp"
#include "cpphub/timeseries/var/dy_spillover.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cpphub {
inline namespace v1 {
namespace c_api {

namespace json = json_min;

// 序列长度上限（DESIGN §7「序列超长（>10k）拒绝」）
inline constexpr std::size_t kMaxSeriesLength = 10000;
// 多变量列数上限（K）
inline constexpr std::size_t kMaxColumns = 64;
// args_json 字节上限（RISK-034）：属工具面输入策略（与序列/列数上限同层），
// 由 c_api.cpp tool_call 入口在 strlen 后前置检查，防超大输入在解析前占用内存。
inline constexpr std::size_t kMaxArgsJsonBytes = 8u * 1024u * 1024u;  // 8 MiB

// ---------------------------------------------------------------------------
// handler 结果
// ---------------------------------------------------------------------------
// 错误分类（映射到 C ABI 错误码，见 DESIGN §7 / CHECKLIST §5）
//   InvalidArg -> CPPHUB_ERR_INVALID_ARG（参数缺失/类型/范围/长度/非有限）
//   Runtime    -> CPPHUB_ERR_RUNTIME（内核抛 std::runtime_error，如退化方差）
enum class ToolErrorKind { InvalidArg, Runtime };

struct ToolOutcome {
    bool ok = false;
    json::Value result;  // ok=true：ToolResult（DESIGN §6.2）
    std::string error;   // ok=false：错误信息（经 cpphub_v1_get_last_error 下发）
    ToolErrorKind kind = ToolErrorKind::InvalidArg;
};

using ToolHandler = ToolOutcome (*)(const json::Value& args);

struct ToolSpec {
    const char* name;
    const char* description;   // 含文献陷阱注释
    const char* input_schema;  // JSON 文本（DESIGN §4.1 订正后入参）
    ToolHandler handler;
};

namespace detail {

using Real = cpphub::Real;
using Size = cpphub::Size;

inline ToolOutcome fail(std::string msg) {
    ToolOutcome o;
    o.ok = false;
    o.error = std::move(msg);
    o.kind = ToolErrorKind::InvalidArg;
    return o;
}

// 内核运行期失败（std::runtime_error 等）→ CPPHUB_ERR_RUNTIME
inline ToolOutcome fail_runtime(std::string msg) {
    ToolOutcome o;
    o.ok = false;
    o.error = std::move(msg);
    o.kind = ToolErrorKind::Runtime;
    return o;
}

inline ToolOutcome succeed(json::Value v) {
    ToolOutcome o;
    o.ok = true;
    o.result = std::move(v);
    return o;
}

// ---- 信封与容器装配（DESIGN §6.1 / §6.2）----

inline json::Value opt_number_or_null(double x) {
    return std::isfinite(x) ? json::Value::number(x) : json::Value::null();
}

inline json::Value make_envelope(const char* tool, const std::string& method,
                                 double statistic, double p_value,
                                 double cv_5pct, bool reject_null,
                                 std::vector<std::string> notes) {
    json::Value env = json::Value::object();
    env.set("tool", json::Value::string(tool));
    env.set("method", json::Value::string(method));
    env.set("statistic", opt_number_or_null(statistic));
    env.set("p_value", opt_number_or_null(p_value));
    env.set("cv_5pct", opt_number_or_null(cv_5pct));
    env.set("reject_null", json::Value::boolean(reject_null));
    json::Value n = json::Value::array();
    for (std::string& s : notes) n.push(json::Value::string(std::move(s)));
    env.set("notes", std::move(n));
    return env;
}

// 主检验行查找：按 method 标签定位；缺失则退化为首行（并补一条 note）
inline json::Value make_tool_result(const char* tool, std::vector<json::Value> rows,
                                    const std::string& primary_label) {
    const json::Value* primary = nullptr;
    for (const json::Value& r : rows) {
        const json::Value* m = r.find("method");
        if (m != nullptr && m->is_string() && m->as_string() == primary_label) {
            primary = &r;
            break;
        }
    }
    if (primary == nullptr && !rows.empty()) primary = &rows.front();

    bool reject_null = false;
    std::string method_name = primary_label;
    if (primary != nullptr) {
        const json::Value* rn = primary->find("reject_null");
        const json::Value* mn = primary->find("method");
        if (rn != nullptr && rn->is_bool()) reject_null = rn->as_bool();
        if (mn != nullptr && mn->is_string()) method_name = mn->as_string();
    }

    json::Value conclusion = json::Value::object();
    conclusion.set("primary_method", json::Value::string(method_name));
    conclusion.set("reject_null", json::Value::boolean(reject_null));
    conclusion.set("n_total", json::Value::number(static_cast<double>(rows.size())));

    json::Value results = json::Value::array();
    for (json::Value& r : rows) results.push(std::move(r));

    json::Value out = json::Value::object();
    out.set("tool", json::Value::string(tool));
    out.set("results", std::move(results));
    out.set("conclusion", std::move(conclusion));
    return out;
}

// ---- 参数读取（缺省 vs 类型错误严格区分）----

// 数值（整数/浮点同源）。缺省返回 false 且不置 err；类型错误返回 false 并置 err。
inline bool opt_number(const json::Value& args, const char* key, double& out,
                       std::string& err) {
    const json::Value* v = args.find(key);
    if (v == nullptr || v->is_null()) return false;
    if (!v->is_number()) {
        err = std::string("parameter '") + key + "' must be a number";
        return false;
    }
    out = v->as_number();
    return true;
}

inline bool opt_int(const json::Value& args, const char* key, long long& out,
                    std::string& err) {
    double d = 0.0;
    const json::Value* v = args.find(key);
    if (v == nullptr || v->is_null()) return false;
    if (!v->is_number()) {
        err = std::string("parameter '") + key + "' must be an integer";
        return false;
    }
    d = v->as_number();
    if (d != std::floor(d) || !std::isfinite(d)) {
        err = std::string("parameter '") + key + "' must be an integer";
        return false;
    }
    out = static_cast<long long>(d);
    return true;
}

inline bool opt_bool(const json::Value& args, const char* key, bool& out,
                     std::string& err) {
    const json::Value* v = args.find(key);
    if (v == nullptr || v->is_null()) return false;
    if (!v->is_bool()) {
        err = std::string("parameter '") + key + "' must be a boolean";
        return false;
    }
    out = v->as_bool();
    return true;
}

inline bool opt_string(const json::Value& args, const char* key, std::string& out,
                       std::string& err) {
    const json::Value* v = args.find(key);
    if (v == nullptr || v->is_null()) return false;
    if (!v->is_string()) {
        err = std::string("parameter '") + key + "' must be a string";
        return false;
    }
    out = v->as_string();
    return true;
}

// 有限性 + 长度校验（CHECKLIST §2.5：非有限值拒绝、长度超限拒绝）
inline bool check_finite(const std::vector<Real>& v, const char* what,
                         std::string& err) {
    if (v.size() > kMaxSeriesLength) {
        err = std::string(what) + " exceeds max length (" +
              std::to_string(kMaxSeriesLength) + ")";
        return false;
    }
    for (Real x : v) {
        if (!std::isfinite(x)) {
            err = std::string(what) + " contains non-finite value (NaN/Inf)";
            return false;
        }
    }
    return true;
}

inline bool req_number_array(const json::Value& args, const char* key,
                             std::vector<Real>& out, std::string& err) {
    const json::Value* v = args.find(key);
    if (v == nullptr || v->is_null()) {
        err = std::string("required parameter '") + key + "' is missing";
        return false;
    }
    if (!v->is_array()) {
        err = std::string("parameter '") + key + "' must be an array of numbers";
        return false;
    }
    out.clear();
    out.reserve(v->as_array().size());
    for (const json::Value& e : v->as_array()) {
        if (!e.is_number()) {
            err = std::string("parameter '") + key + "' must be an array of numbers";
            return false;
        }
        out.push_back(e.as_number());
    }
    return check_finite(out, key, err);
}

// number[][]（列优先：外层为变量列，与 MultivariateTSData.columns 一致）
inline bool req_matrix(const json::Value& args, const char* key,
                       std::vector<std::vector<Real>>& out, std::string& err) {
    const json::Value* v = args.find(key);
    if (v == nullptr || v->is_null()) {
        err = std::string("required parameter '") + key + "' is missing";
        return false;
    }
    if (!v->is_array()) {
        err = std::string("parameter '") + key + "' must be an array of arrays";
        return false;
    }
    out.clear();
    if (v->as_array().size() > kMaxColumns) {
        err = std::string(key) + " has too many columns";
        return false;
    }
    for (const json::Value& col : v->as_array()) {
        if (!col.is_array()) {
            err = std::string("parameter '") + key + "' must be an array of arrays";
            return false;
        }
        std::vector<Real> c;
        c.reserve(col.as_array().size());
        for (const json::Value& e : col.as_array()) {
            if (!e.is_number()) {
                err = std::string("parameter '") + key +
                      "' must be an array of arrays of numbers";
                return false;
            }
            c.push_back(e.as_number());
        }
        if (!check_finite(c, key, err)) return false;
        out.push_back(std::move(c));
    }
    return true;
}

inline bool opt_int_array(const json::Value& args, const char* key,
                          std::vector<Size>& out, std::string& err) {
    const json::Value* v = args.find(key);
    if (v == nullptr || v->is_null()) return false;
    if (!v->is_array()) {
        err = std::string("parameter '") + key + "' must be an array of integers";
        return false;
    }
    out.clear();
    for (const json::Value& e : v->as_array()) {
        if (!e.is_number()) {
            err = std::string("parameter '") + key + "' must be an array of integers";
            return false;
        }
        const double d = e.as_number();
        if (d != std::floor(d) || !std::isfinite(d) || d < 0.0) {
            err = std::string("parameter '") + key +
                  "' must be an array of non-negative integers";
            return false;
        }
        out.push_back(static_cast<Size>(d));
    }
    return true;
}

inline bool opt_string_array(const json::Value& args, const char* key,
                             std::vector<std::string>& out, std::string& err) {
    const json::Value* v = args.find(key);
    if (v == nullptr || v->is_null()) return false;
    if (!v->is_array()) {
        err = std::string("parameter '") + key + "' must be an array of strings";
        return false;
    }
    out.clear();
    for (const json::Value& e : v->as_array()) {
        if (!e.is_string()) {
            err = std::string("parameter '") + key + "' must be an array of strings";
            return false;
        }
        out.push_back(e.as_string());
    }
    return true;
}

inline bool contains(const std::vector<std::string>& v, const char* s) {
    for (const std::string& e : v) {
        if (e == s) return true;
    }
    return false;
}

// ===========================================================================
// Handler 1：cpphub_unit_root_diagnose（ADF/PP/KPSS/DF-GLS/NP/ZA 一揽子）
// ===========================================================================
inline ToolOutcome handle_unit_root(const json::Value& args) {
    static const char* kTool = "cpphub_unit_root_diagnose";
    std::vector<Real> series;
    std::string err;
    if (!req_number_array(args, "series", series, err)) return fail(err);

    std::string trend = "auto";
    opt_string(args, "trend", trend, err);
    if (!err.empty()) return fail(err);
    if (trend != "auto" && trend != "nc" && trend != "n" && trend != "c" &&
        trend != "ct") {
        return fail("parameter 'trend' must be auto/nc/n/c/ct");
    }
    long long max_lag = 0;
    if (opt_int(args, "max_lag", max_lag, err) && max_lag < 0) {
        return fail("parameter 'max_lag' must be >= 0");
    }
    if (!err.empty()) return fail(err);

    std::string za_model = "C";
    opt_string(args, "za_model", za_model, err);
    if (!err.empty()) return fail(err);

    std::vector<std::string> tests = {"adf", "pp", "kpss", "df_gls",
                                     "ng_perron", "za"};
    opt_string_array(args, "tests", tests, err);
    if (!err.empty()) return fail(err);
    for (const std::string& t : tests) {
        if (t != "adf" && t != "pp" && t != "kpss" && t != "df_gls" &&
            t != "ng_perron" && t != "za") {
            return fail("parameter 'tests' has unknown test: " + t);
        }
    }

    using namespace cpphub::v1::timeseries::unit_root;
    const Size ml = static_cast<Size>(max_lag);
    // KPSS/DF-GLS/NP 仅接受 c/ct；trend=="auto" 时沿用各检验自身默认
    const bool c_ct_only = (trend == "c" || trend == "ct");
    const std::string t_kpss = c_ct_only ? trend : "c";
    const std::string t_gls = c_ct_only ? trend : "ct";
    const std::string t_np = c_ct_only ? trend : "ct";

    std::vector<json::Value> rows;

    // 单检验执行：库内抛异常（样本过小 / 含 NaN / 参数非法）转为该行
    // statistic=null + notes 失败原因，保持"一揽子"形态不中断整次调用
    const double nan = std::numeric_limits<double>::quiet_NaN();
    auto run = [&](const char* label, auto&& fn) {
        try {
            rows.push_back(fn());
        } catch (const std::exception& e) {
            rows.push_back(make_envelope(kTool, label, nan, nan, nan, false,
                                         {std::string("failed: ") + e.what()}));
        }
    };

    if (contains(tests, "adf")) {
        run("adf", [&] {
            const ADFResult r = adf_test(series, trend, ml, false);
            return make_envelope(
                kTool, "adf", r.statistic, r.p_value, r.critical_value_5pct,
                r.reject_null,
                {"H0=单位根(非平稳); H1=平稳 (U4: tau 非标准分布)",
                 "n_lags=" + std::to_string(r.n_lags) + "; trend=" + r.trend_spec});
        });
    }
    if (contains(tests, "pp")) {
        run("pp", [&] {
            const PPResult r = pp_test(series, trend, 0);
            return make_envelope(
                kTool, "pp", r.statistic, r.p_value, r.critical_value_5pct,
                r.reject_null,
                {"H0=单位根(非平稳); H1=平稳",
                 "bandwidth=" + std::to_string(r.bandwidth) +
                     "; trend=" + r.trend_spec});
        });
    }
    if (contains(tests, "kpss")) {
        run("kpss", [&] {
            const KPSSResult r = kpss_test(series, t_kpss, 0);
            return make_envelope(
                kTool, "kpss", r.statistic, r.p_value, r.critical_value_5pct,
                r.reject_null,
                {"U12 陷阱: KPSS 的 H0 是【平稳性】, 与 ADF/PP 的 H0(单位根)【相反】",
                 "reject_null=true 表示拒绝平稳性(即支持单位根)",
                 "bandwidth=" + std::to_string(r.bandwidth) +
                     "; trend=" + r.trend_spec});
        });
    }
    if (contains(tests, "df_gls")) {
        run("df_gls", [&] {
            const DFGlsResult r = df_gls_test(series, t_gls, ml);
            return make_envelope(
                kTool, "df_gls", r.statistic, r.p_value, r.critical_value_5pct,
                r.reject_null,
                {"H0=单位根(非平稳); H1=平稳 (U21: 独立 response surface)",
                 "c_bar=" + std::to_string(r.c_bar) + "; n_lags=" +
                     std::to_string(r.n_lags) + "; trend=" + r.trend_spec});
        });
    }
    if (contains(tests, "ng_perron")) {
        run("ng_perron", [&] {
            const NgPerronResult r = ng_perron_test(series, t_np, ml);
            return make_envelope(
                kTool, "ng_perron", r.mz_t, std::numeric_limits<double>::quiet_NaN(),
                r.cv_5pct[1], r.reject_5pct[1],
                {"Ng-Perron 四统计量同 H0(单位根); 本行报 MZt(索引 1)",
                 "cv 数组顺序 = {MZalpha, MZt, MSB, MPT}; MZalpha=" +
                     std::to_string(r.mz_alpha) + "; MSB=" + std::to_string(r.msb) +
                     "; MPT=" + std::to_string(r.mpt),
                 "selected_lag(MAIC)=" + std::to_string(r.selected_lag) +
                     "; trend=" + r.trend_spec});
        });
    }
    if (contains(tests, "za")) {
        run("za", [&] {
            ZAModel model = ZAModel::C;
            if (za_model == "A" || za_model == "a") model = ZAModel::A;
            else if (za_model == "B" || za_model == "b") model = ZAModel::B;
            else if (za_model != "C" && za_model != "c") {
                throw std::invalid_argument("za_model must be A/B/C");
            }
            const ZAResult r = zivot_andrews_test(series, model, ml, false, 0.15);
            return make_envelope(
                kTool, "za", r.statistic, r.p_value_baum, r.critical_5pct,
                r.reject_null,
                {"H0=带单断点单位根 (ZA1992); 临界值取论文表, p 取 MC 插值(不同源)",
                 "break_index(0-based)=" + std::to_string(r.break_index) +
                     "; model=" + za_model + "; trim=" + std::to_string(r.trim) +
                     "; n_lags=" + std::to_string(r.n_lags)});
        });
    }

    if (rows.empty()) return fail("parameter 'tests' selected no test");
    return succeed(make_tool_result(kTool, std::move(rows), "adf"));
}

// ===========================================================================
// Handler 2：cpphub_variance_ratio_test
// ===========================================================================
inline ToolOutcome handle_variance_ratio(const json::Value& args) {
    static const char* kTool = "cpphub_variance_ratio_test";
    std::vector<Real> series;
    std::string err;
    if (!req_number_array(args, "series", series, err)) return fail(err);

    std::vector<Size> k_list = {2, 5, 10, 20};
    opt_int_array(args, "k", k_list, err);
    if (!err.empty()) return fail(err);
    if (k_list.empty()) return fail("parameter 'k' must be a non-empty array");
    for (Size k : k_list) {
        if (k < 2) return fail("parameter 'k' must have all horizons >= 2");
    }
    bool use_debiased = true;
    opt_bool(args, "use_debiased", use_debiased, err);
    if (!err.empty()) return fail(err);
    std::string trend = "c";
    opt_string(args, "trend", trend, err);
    if (!err.empty()) return fail(err);
    if (trend != "c" && trend != "n") return fail("parameter 'trend' must be c/n");

    try {
        using namespace cpphub::v1::timeseries::unit_root;
        const VarianceRatioResult r =
            variance_ratio_test_multi(series, k_list, use_debiased, trend);

        std::string horizon_note = "VR(k) 逐 horizon: ";
        for (std::size_t i = 0; i < r.vr_list.size(); ++i) {
            if (i != 0) horizon_note += ", ";
            horizon_note += "k=" + std::to_string(r.k_values[i]) + " -> " +
                            std::to_string(r.vr_list[i]);
        }

        std::vector<json::Value> rows;
        rows.push_back(make_envelope(
            kTool, "lo_mackinlay_z1", r.z1_statistic, r.z1_p_value,
            std::numeric_limits<double>::quiet_NaN(), r.z1_p_value < 0.05,
            {"H0=随机游走 (Lo-MacKinlay 1988 同方差 Z1, U14)", horizon_note}));
        rows.push_back(make_envelope(
            kTool, "lo_mackinlay_z2", r.z2_statistic, r.z2_p_value,
            std::numeric_limits<double>::quiet_NaN(), r.reject_null,
            {"H0=随机游走 (异方差稳健 Z2, U15; 主检验)",
             "trend=" + trend + "; debiased=" + (use_debiased ? "true" : "false"),
             horizon_note}));
        rows.push_back(make_envelope(
            kTool, "chow_denning_joint", r.chow_denning_statistic,
            r.chow_denning_p_value, std::numeric_limits<double>::quiet_NaN(),
            r.chow_denning_p_value < 0.05,
            {"多 horizon 联合 SMM(m,inf) (U16: 多重检验校正; horizon 间正相关 → p 偏保守)",
             horizon_note}));
        rows.push_back(make_envelope(
            kTool, "clm_debiased", r.clm_debiased_statistic, r.clm_debiased_p_value,
            std::numeric_limits<double>::quiet_NaN(), r.clm_debiased_p_value < 0.05,
            {"CLM 小样本 debiased Z2 (U17, 恒 debias)", horizon_note}));

        return succeed(make_tool_result(kTool, std::move(rows), "lo_mackinlay_z2"));
    } catch (const std::exception& e) {
        return fail_runtime(std::string("variance_ratio_test failed: ") + e.what());
    }
}

// ===========================================================================
// Handler 3：cpphub_cointegration_scan（EG + Johansen）
// ===========================================================================
inline ToolOutcome handle_cointegration(const json::Value& args) {
    static const char* kTool = "cpphub_cointegration_scan";
    std::vector<std::vector<Real>> cols;
    std::string err;
    if (!req_matrix(args, "series", cols, err)) return fail(err);
    if (cols.size() < 2) return fail("parameter 'series' needs at least 2 columns");

    std::string eg_trend = "c";
    opt_string(args, "eg_trend", eg_trend, err);
    if (!err.empty()) return fail(err);
    long long det_order = 0;
    if (opt_int(args, "det_order", det_order, err) &&
        (det_order < -1 || det_order > 1)) {
        return fail("parameter 'det_order' must be -1/0/1");
    }
    if (!err.empty()) return fail(err);
    long long k_ar_diff = 0;
    if (opt_int(args, "k_ar_diff", k_ar_diff, err) && k_ar_diff < 0) {
        return fail("parameter 'k_ar_diff' must be >= 0");
    }
    if (!err.empty()) return fail(err);

    using namespace cpphub::v1::timeseries;
    const Size T0 = cols[0].size();
    for (const std::vector<Real>& c : cols) {
        if (c.size() != T0) return fail("parameter 'series' columns must be equal length");
    }

    std::vector<json::Value> rows;

    if (cols.size() == 2) {
        try {
            const cointegration::EGResult r =
                cointegration::engle_granger(cols[0], cols[1], eg_trend);
            rows.push_back(make_envelope(
                kTool, "engle_granger", r.statistic, r.p_value, r.cv_5pct,
                r.reject_null,
                {"CI3 方向依赖: 参数顺序 (y0,y1) 中 y0 为被解释变量; 交换列会改变结果",
                 "CI2: p 值(MacKinnon 1994) 与临界值(MacKinnon 2010) 不同源",
                 "trend=" + eg_trend + "; n_obs=" + std::to_string(r.n_obs)}));
        } catch (const std::exception& e) {
            rows.push_back(make_envelope(
                kTool, "engle_granger",
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {std::string("failed: ") + e.what()}));
        }
    } else {
        rows.push_back(make_envelope(
            kTool, "engle_granger", std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(), false,
            {"K != 2: Engle-Granger 两步法仅适用于双变量, 本行未计算"}));
    }

    try {
        var::MultivariateTSData data;
        data.columns = cols;
        const cointegration::JohansenResult jr = cointegration::coint_johansen(
            data, static_cast<int>(det_order), static_cast<Size>(k_ar_diff));
        const Size K = static_cast<Size>(jr.lr1.size());
        for (Size r = 0; r < K; ++r) {
            const double cv = jr.cvt(static_cast<Eigen::Index>(r), 1);
            rows.push_back(make_envelope(
                kTool, "johansen_trace_r" + std::to_string(r), jr.lr1(r), 
                std::numeric_limits<double>::quiet_NaN(), cv, jr.lr1(r) > cv,
                {"H0: rank <= r; H1: rank > r (迹统计量; 拒绝方向 stat > cv)",
                 "cv_source=" + jr.cv_source + "; cv 列序 {90,95,99} -> 取 95",
                 "det_order=" + std::to_string(det_order) +
                     "; k_ar_diff=" + std::to_string(k_ar_diff) +
                     "; n_obs=" + std::to_string(jr.n_obs)}));
            const double cvm95 = jr.cvm(static_cast<Eigen::Index>(r), 1);
            rows.push_back(make_envelope(
                kTool, "johansen_maxeig_r" + std::to_string(r), jr.lr2(r),
                std::numeric_limits<double>::quiet_NaN(), cvm95,
                jr.lr2(r) > cvm95,
                {"H0: rank <= r; H1: rank > r (最大特征值; 拒绝方向 stat > cv)",
                 "CI5: cvt/cvm 双表源 (OL1992 / MHM96)", "cv_source=" + jr.cv_source}));
        }
    } catch (const std::exception& e) {
        rows.push_back(make_envelope(
            kTool, "johansen_trace_r0",
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(), false,
            {std::string("failed: ") + e.what()}));
    }

    const std::string primary = (cols.size() == 2) ? "engle_granger"
                                                  : "johansen_trace_r0";
    return succeed(make_tool_result(kTool, std::move(rows), primary));
}

// ===========================================================================
// Handler 4：cpphub_vecm_fit
// ===========================================================================
inline ToolOutcome handle_vecm(const json::Value& args) {
    static const char* kTool = "cpphub_vecm_fit";
    std::vector<std::vector<Real>> cols;
    std::string err;
    if (!req_matrix(args, "series", cols, err)) return fail(err);
    if (cols.size() < 2) return fail("parameter 'series' needs at least 2 columns");

    long long rank = 0;
    if (args.find("rank") == nullptr || args.find("rank")->is_null()) {
        return fail("required parameter 'rank' is missing");
    }
    if (!opt_int(args, "rank", rank, err)) return fail(err);
    if (rank < 1 || static_cast<std::size_t>(rank) > cols.size()) {
        return fail("parameter 'rank' must be in [1, K]");
    }
    if (!err.empty()) return fail(err);

    long long k_ar_diff = 0;
    if (opt_int(args, "k_ar_diff", k_ar_diff, err) && k_ar_diff < 0) {
        return fail("parameter 'k_ar_diff' must be >= 0");
    }
    if (!err.empty()) return fail(err);

    std::string det = "n";
    opt_string(args, "det", det, err);
    if (!err.empty()) return fail(err);
    if (det != "n" && det != "co" && det != "ci" && det != "lo" && det != "li") {
        return fail("parameter 'det' must be n/co/ci/lo/li");
    }
    bool urca_norm = false;
    opt_bool(args, "urca_normalization", urca_norm, err);
    if (!err.empty()) return fail(err);

    for (const std::vector<Real>& c : cols) {
        if (c.size() != cols[0].size()) {
            return fail("parameter 'series' columns must be equal length");
        }
    }

    std::vector<json::Value> rows;
    try {
        using namespace cpphub::v1::timeseries;
        var::MultivariateTSData data;
        data.columns = cols;
        const cointegration::VECMResult r = cointegration::vecm_fit(
            data, static_cast<Size>(rank), static_cast<Size>(k_ar_diff), det,
            urca_norm);

        if (r.has_ect_t) {
            const Size K = static_cast<Size>(r.ect_t_stat.size());
            for (Size i = 0; i < K; ++i) {
                const double t = r.ect_t_stat[i];
                const double cv = (i < r.ect_cv_5pct.size())
                                      ? r.ect_cv_5pct[i]
                                      : std::numeric_limits<double>::quiet_NaN();
                rows.push_back(make_envelope(
                    kTool, "ect_t_" + std::to_string(i),
                    t, std::numeric_limits<double>::quiet_NaN(), cv,
                    std::isfinite(cv) ? (std::fabs(t) > cv) : false,
                    {"ECT t 统计量 (EM2002 5% 临界值, 双尾 |t| > cv)",
                     "仅 rank=1 时可用; det=" + det + "; rank=" + std::to_string(rank)}));
            }
        }
        rows.push_back(make_envelope(
            kTool, "log_likelihood", r.loglik,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(), false,
            {"非检验行: VECM ML 对数似然 (供模型比较)",
             "rank=" + std::to_string(r.rank) + "; det=" + r.det +
                 "; k_ar_diff=" + std::to_string(r.k_ar_diff) +
                 "; n_obs=" + std::to_string(r.n_obs) +
                 "; PI=alpha*beta' (CI9; alpha/beta 见库内 VECMResult)"}));

        const std::string primary = r.has_ect_t ? "ect_t_0" : "log_likelihood";
        return succeed(make_tool_result(kTool, std::move(rows), primary));
    } catch (const std::exception& e) {
        return fail_runtime(std::string("vecm_fit failed: ") + e.what());
    }
}

// ===========================================================================
// Handler 5：cpphub_garch_diagnose
// ===========================================================================
inline ToolOutcome handle_garch(const json::Value& args) {
    static const char* kTool = "cpphub_garch_diagnose";
    std::vector<Real> series;
    std::string err;
    if (!req_number_array(args, "series", series, err)) return fail(err);

    std::string model = "";
    opt_string(args, "model", model, err);
    if (!err.empty()) return fail(err);
    if (model.empty()) return fail("required parameter 'model' is missing");
    if (model != "GARCH" && model != "EGARCH" && model != "GJR" &&
        model != "GARCH-M") {
        return fail("parameter 'model' must be GARCH/EGARCH/GJR/GARCH-M");
    }
    std::string dist = "normal";
    opt_string(args, "dist", dist, err);
    if (!err.empty()) return fail(err);
    std::string form = "vol";
    opt_string(args, "form", form, err);
    if (!err.empty()) return fail(err);

    using namespace cpphub::v1::timeseries::garch;
    GarchDist d = GarchDist::Normal;
    if (dist == "t" || dist == "studentt") d = GarchDist::StudentT;
    else if (dist == "ged") d = GarchDist::GED;
    else if (dist != "normal") return fail("parameter 'dist' must be normal/t/ged");

    GarchConfig cfg;
    cfg.dist = d;
    cfg.compute_sandwich = true;
    cfg.compute_diagnostics = true;

    std::vector<Real> std_res;
    std::vector<json::Value> rows;
    std::string primary = "omega";  // 主检验行标签（GARCH-M 时为 lambda）
    try {
        if (model == "GARCH") {
            const GarchResult r = estimate_garch11(series, cfg);
            std_res = r.std_residuals;
            const double persistence = r.params.alpha + r.params.beta;
            rows.push_back(make_envelope(
                kTool, "omega", r.params.omega, 
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"参数估计 (QMLE); std_error=" + std::to_string(
                     r.std_errors.size() > 0 ? r.std_errors[0] : 0.0) +
                     " (sandwich)"}));
            rows.push_back(make_envelope(
                kTool, "alpha", r.params.alpha,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"ARCH 系数; std_error=" + std::to_string(
                     r.std_errors.size() > 1 ? r.std_errors[1] : 0.0)}));
            rows.push_back(make_envelope(
                kTool, "beta", r.params.beta,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"GARCH 系数; std_error=" + std::to_string(
                     r.std_errors.size() > 2 ? r.std_errors[2] : 0.0),
                 "稳定性: alpha+beta=" + std::to_string(persistence) +
                     (persistence < 1.0 ? " (<1, 协方差平稳)" : " (>=1, 非平稳!)")}));
            rows.push_back(make_envelope(
                kTool, "log_likelihood", r.log_likelihood,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"非检验行; AIC=" + std::to_string(r.aic) +
                     "; BIC=" + std::to_string(r.bic) +
                     "; converged=" + (r.converged ? "true" : "false") +
                     "; dist=" + dist}));
        } else if (model == "EGARCH") {
            const EGarchResult r = estimate_egarch(series, cfg);
            std_res = r.std_residuals;
            rows.push_back(make_envelope(
                kTool, "omega", r.params.omega,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"EGARCH 参数 (log 方差方程); std_error=" + std::to_string(
                     r.std_errors.size() > 0 ? r.std_errors[0] : 0.0)}));
            rows.push_back(make_envelope(
                kTool, "alpha", r.params.alpha,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"非对称系数 (乘 z_{t-1}); std_error=" + std::to_string(
                     r.std_errors.size() > 1 ? r.std_errors[1] : 0.0)}));
            rows.push_back(make_envelope(
                kTool, "gamma", r.params.gamma,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"对称 size 系数 (乘 |z|-E|z|); std_error=" + std::to_string(
                     r.std_errors.size() > 2 ? r.std_errors[2] : 0.0)}));
            rows.push_back(make_envelope(
                kTool, "beta", r.params.beta,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"持续性 beta (log 方差); |beta|<1 平稳",
                 "稳定性: |beta|=" + std::to_string(std::fabs(r.params.beta))}));
            rows.push_back(make_envelope(
                kTool, "log_likelihood", r.log_likelihood,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"非检验行; AIC=" + std::to_string(r.aic) +
                     "; BIC=" + std::to_string(r.bic) + "; dist=" + dist}));
        } else if (model == "GJR") {
            const GjrGarchResult r = estimate_gjr_garch(series, cfg);
            std_res = r.std_residuals;
            const double persistence =
                r.params.alpha + 0.5 * r.params.gamma + r.params.beta;
            rows.push_back(make_envelope(
                kTool, "omega", r.params.omega,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"GJR-GARCH 参数; std_error=" + std::to_string(
                     r.std_errors.size() > 0 ? r.std_errors[0] : 0.0)}));
            rows.push_back(make_envelope(
                kTool, "alpha", r.params.alpha,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"ARCH 系数; std_error=" + std::to_string(
                     r.std_errors.size() > 1 ? r.std_errors[1] : 0.0)}));
            rows.push_back(make_envelope(
                kTool, "gamma", r.params.gamma,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"非对称(杠杆)系数: gamma>0 负冲击放大波动 (G8); std_error=" +
                     std::to_string(r.std_errors.size() > 2 ? r.std_errors[2] : 0.0)}));
            rows.push_back(make_envelope(
                kTool, "beta", r.params.beta,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"GARCH 系数; std_error=" + std::to_string(
                     r.std_errors.size() > 3 ? r.std_errors[3] : 0.0),
                 "稳定性: alpha+gamma/2+beta=" + std::to_string(persistence)}));
            rows.push_back(make_envelope(
                kTool, "log_likelihood", r.log_likelihood,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"非检验行; AIC=" + std::to_string(r.aic) +
                     "; BIC=" + std::to_string(r.bic) + "; dist=" + dist}));
        } else {  // GARCH-M
            GarchMForm f = GarchMForm::Volatility;
            if (form == "var") f = GarchMForm::Variance;
            else if (form == "vol") f = GarchMForm::Volatility;
            else if (form == "log") f = GarchMForm::LogVariance;
            else return fail("parameter 'form' must be var/vol/log");
            const GarchMResult r = estimate_garch_m(series, f, cfg);
            std_res = r.std_residuals;
            primary = "lambda";
            rows.push_back(make_envelope(
                kTool, "lambda", r.params.lambda,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"风险溢价系数 (GARCH-M 关键参数, GM1: g(h) 由 form 决定)",
                 "form=" + form + "; std_error=" + std::to_string(
                     r.std_errors.size() > 1 ? r.std_errors[1] : 0.0)}));
            rows.push_back(make_envelope(
                kTool, "mu", r.params.mu,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"均值方程截距"}));
            rows.push_back(make_envelope(
                kTool, "log_likelihood", r.log_likelihood,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false,
                {"非检验行; AIC=" + std::to_string(r.aic) +
                     "; BIC=" + std::to_string(r.bic) +
                     "; alpha+beta=" + std::to_string(r.params.alpha + r.params.beta)}));
        }

        // L2 标准化残差诊断（LB(z) / LB(z^2) / JB）
        if (!std_res.empty()) {
            const GarchDiagnosticsResult g =
                diagnose_garch_residuals(std_res, 1000, 0);
            rows.push_back(make_envelope(
                kTool, "ljung_box_z", g.lb_z.base.statistic, g.lb_z.base.p_value,
                std::numeric_limits<double>::quiet_NaN(),
                g.lb_z.base.reject_null,
                {"H0=标准化残差无自相关 (LB(z), G11); lag=" +
                     std::to_string(g.lb_z.lag),
                 std::string("passes_no_autocorr=") +
                     (g.passes_no_autocorr ? "true" : "false")}));
            rows.push_back(make_envelope(
                kTool, "ljung_box_z_squared", g.lb_z_squared.base.statistic,
                g.lb_z_squared.base.p_value,
                std::numeric_limits<double>::quiet_NaN(),
                g.lb_z_squared.base.reject_null,
                {"H0=无 ARCH 效应 (LB(z^2), G12/Li-Mak); lag=" +
                     std::to_string(g.lb_z_squared.lag),
                 std::string("passes_no_arch_effect=") +
                     (g.passes_no_arch_effect ? "true" : "false")}));
            rows.push_back(make_envelope(
                kTool, "jarque_bera", g.jb_test.base.statistic,
                g.jb_test.base.p_value, std::numeric_limits<double>::quiet_NaN(),
                g.jb_test.base.reject_null,
                {"H0=标准化残差正态 (JB; G-ADR4: p 经 Bootstrap, reps=" +
                     std::to_string(g.jb_bootstrap_reps) + ")",
                 "skewness=" + std::to_string(g.jb_test.skewness) +
                     "; kurtosis(非超额)=" + std::to_string(g.jb_test.kurtosis),
                 std::string("passes_normality=") +
                     (g.passes_normality ? "true" : "false")}));
        }
    } catch (const std::exception& e) {
        return fail_runtime(std::string("garch_diagnose failed: ") + e.what());
    }

    if (rows.empty()) return fail("garch_diagnose produced no results");
    return succeed(make_tool_result(kTool, std::move(rows), primary));
}

// ===========================================================================
// Handler 6：cpphub_granger_causality
// ===========================================================================
inline ToolOutcome handle_granger(const json::Value& args) {
    static const char* kTool = "cpphub_granger_causality";
    std::vector<Real> cause, effect;
    std::string err;
    if (!req_number_array(args, "cause", cause, err)) return fail(err);
    if (!req_number_array(args, "effect", effect, err)) return fail(err);
    if (cause.size() != effect.size()) {
        return fail("parameters 'cause' and 'effect' must have equal length");
    }

    long long lag = 1;
    if (opt_int(args, "lag", lag, err) && lag < 1) {
        return fail("parameter 'lag' must be >= 1");
    }
    if (!err.empty()) return fail(err);
    long long d_max = 0;
    if (opt_int(args, "d_max", d_max, err) && d_max < 0) {
        return fail("parameter 'd_max' must be >= 0");
    }
    if (!err.empty()) return fail(err);
    bool with_hac = true;
    opt_bool(args, "with_hac", with_hac, err);
    if (!err.empty()) return fail(err);
    long long bw = 0;
    if (opt_int(args, "hac_bandwidth", bw, err) && bw < 0) {
        return fail("parameter 'hac_bandwidth' must be >= 0");
    }
    if (!err.empty()) return fail(err);

    try {
        using namespace cpphub::v1::timeseries::granger;
        const GrangerResult r = granger_test(
            cause, effect, static_cast<Size>(lag), static_cast<Size>(d_max),
            with_hac, static_cast<Size>(bw));

        std::vector<json::Value> rows;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        rows.push_back(make_envelope(
            kTool, "ssr_f", r.f_stat, r.f_p, nan, r.f_p < 0.05,
            {"H0=cause 不 Granger-cause effect (statsmodels ssr_ftest, 主检验)",
             "GR6 方向陷阱: 参数顺序显式 (cause, effect); 与 statsmodels 列序相反",
             "df1=" + std::to_string(r.df1) + "; df2=" + std::to_string(r.df2)}));
        rows.push_back(make_envelope(
            kTool, "ssr_chi2", r.chi2_stat, r.chi2_p, nan, r.chi2_p < 0.05,
            {"H0 同上 (ssr_chi2test, 渐近 chi2)"}));
        rows.push_back(make_envelope(
            kTool, "lr", r.lr_stat, r.lr_p, nan, r.lr_p < 0.05,
            {"H0 同上 (likelihood-ratio)"}));
        if (r.has_ty) {
            rows.push_back(make_envelope(
                kTool, "ty_wald", r.ty_wald_stat, r.ty_wald_p, nan,
                r.ty_wald_p < 0.05,
                {"Toda-Yamamoto 增广 Wald (d_max=" + std::to_string(d_max) +
                     "); GR4: 增广项不进入约束"}));
        }
        if (r.has_hac) {
            rows.push_back(make_envelope(
                kTool, "hac_wald", r.hac_wald_stat, r.hac_wald_p, nan,
                r.hac_wald_p < 0.05,
                {"HAC 稳健 Wald (Newey-West); GR5: 其余统计量均非异方差稳健",
                 "hac_bandwidth=" + std::to_string(bw) + " (0=NW 经验法则)"}));
        }
        return succeed(make_tool_result(kTool, std::move(rows), "ssr_f"));
    } catch (const std::exception& e) {
        return fail_runtime(std::string("granger_test failed: ") + e.what());
    }
}

// ===========================================================================
// Handler 7：cpphub_dy_spillover
// ===========================================================================
inline ToolOutcome handle_dy(const json::Value& args) {
    static const char* kTool = "cpphub_dy_spillover";
    std::vector<std::vector<Real>> cols;
    std::string err;
    if (!req_matrix(args, "series", cols, err)) return fail(err);
    const std::size_t K = cols.size();
    if (K < 2) return fail("parameter 'series' needs at least 2 columns");
    for (const std::vector<Real>& c : cols) {
        if (c.size() != cols[0].size()) {
            return fail("parameter 'series' columns must be equal length");
        }
    }

    long long horizon = 10;
    if (opt_int(args, "horizon", horizon, err) && horizon < 1) {
        return fail("parameter 'horizon' must be >= 1");
    }
    if (!err.empty()) return fail(err);
    long long lag = 0;
    if (opt_int(args, "lag", lag, err) && lag < 0) {
        return fail("parameter 'lag' must be >= 0");
    }
    if (!err.empty()) return fail(err);
    std::string trend = "c";
    opt_string(args, "trend", trend, err);
    if (!err.empty()) return fail(err);
    if (trend != "c" && trend != "n") return fail("parameter 'trend' must be c/n");
    long long rolling = 0;
    if (opt_int(args, "rolling", rolling, err) && rolling < 0) {
        return fail("parameter 'rolling' must be >= 0");
    }
    if (!err.empty()) return fail(err);
    if (rolling > 0 && static_cast<std::size_t>(rolling) <= 2 * K) {
        return fail("parameter 'rolling' (window) must be > 2*K (K=" +
                    std::to_string(K) + ")");
    }

    try {
        using namespace cpphub::v1::timeseries::var;
        MultivariateTSData data;
        data.columns = cols;
        const DYResult r =
            (rolling > 0)
                ? dy_spillover(data, static_cast<Size>(rolling),
                               static_cast<Size>(horizon), trend,
                               static_cast<Size>(lag))
                : dy_spillover_static(data, static_cast<Size>(horizon), trend,
                                      static_cast<Size>(lag));

        std::vector<json::Value> rows;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        std::string layout = "变量序 = series 列序 (TO/FROM/NET 按变量)";
        rows.push_back(make_envelope(
            kTool, "tci", r.tci, nan, nan, false,
            {"总溢出指数 TCI (Diebold-Yilmaz, 全部 x100/K 口径; 非检验, 无 H0)",
             layout, "K=" + std::to_string(data.K()) +
                 "; T=" + std::to_string(data.T()) +
                 "; horizon=" + std::to_string(horizon) +
                 "; lag=" + std::to_string(r.lag) +
                 (rolling > 0 ? "; window=" + std::to_string(r.window) : "; 全样本静态")}));
        for (std::size_t v = 0; v < r.to_spillover.size(); ++v) {
            const std::string idx = std::to_string(v);
            rows.push_back(make_envelope(kTool, "to_v" + idx, r.to_spillover[v],
                                         nan, nan, false,
                                         {"变量 " + idx + " 的 TO 溢出 (对外)"}));
            rows.push_back(make_envelope(kTool, "from_v" + idx,
                                         r.from_spillover[v], nan, nan, false,
                                         {"变量 " + idx + " 的 FROM 溢出 (接收)"}));
            rows.push_back(make_envelope(kTool, "net_v" + idx,
                                         r.net_spillover[v], nan, nan, false,
                                         {"NET = TO - FROM (>0 净输出)"}));
        }
        return succeed(make_tool_result(kTool, std::move(rows), "tci"));
    } catch (const std::exception& e) {
        return fail_runtime(std::string("dy_spillover failed: ") + e.what());
    }
}

}  // namespace detail

// ---------------------------------------------------------------------------
// 注册表（单一登记点；惰性构建 + 只读缓存）
// ---------------------------------------------------------------------------
inline const std::vector<ToolSpec>& tool_specs() {
    static const std::vector<ToolSpec> kSpecs = {
        {"cpphub_unit_root_diagnose",
         "单位根六检验一揽子诊断 (ADF/PP/KPSS/DF-GLS/Ng-Perron/Zivot-Andrews) + "
         "逐方法信封与综合结论。陷阱注释: KPSS 的 H0 是平稳性, 与 ADF/PP 相反 (U12), "
         "勿把 reject_null 跨方法直接相加; ZA 为带单断点 H0。",
         R"json({"type":"object","properties":{"series":{"type":"array","items":{"type":"number"},"description":"价格/收益水平序列"},"trend":{"type":"string","enum":["auto","nc","c","ct"],"default":"auto","description":"趋势设定; KPSS/DF-GLS/Ng-Perron 仅接受 c/ct"},"tests":{"type":"array","items":{"type":"string","enum":["adf","pp","kpss","df_gls","ng_perron","za"]},"description":"要运行的检验子集, 缺省为全部"},"za_model":{"type":"string","enum":["A","B","C"],"default":"C"},"max_lag":{"type":"integer","minimum":0,"default":0,"description":"0=自动 (Schwert/AIC)"}},"required":["series"],"additionalProperties":false})json",
         detail::handle_unit_root},
        {"cpphub_variance_ratio_test",
         "方差比检验 (Lo-MacKinlay Z1 同方差 / Z2 异方差稳健 / Chow-Denning 多 horizon "
         "联合 / CLM debiased)。H0=随机游走。",
         R"json({"type":"object","properties":{"series":{"type":"array","items":{"type":"number"}},"k":{"type":"array","items":{"type":"integer","minimum":2},"default":[2,5,10,20],"description":"horizon 列表, 逐项 >= 2"},"use_debiased":{"type":"boolean","default":true},"trend":{"type":"string","enum":["c","n"],"default":"c"}},"required":["series"],"additionalProperties":false})json",
         detail::handle_variance_ratio},
        {"cpphub_cointegration_scan",
         "协整扫描: Engle-Granger 两步法 (K=2) + Johansen 迹/最大特征值 + rank 判定。"
         "陷阱注释: EG 方向依赖 (CI3, y0=被解释变量, 交换列会变结果); EG 的 p 值与临界值不同源 (CI2)。",
         R"json({"type":"object","properties":{"series":{"type":"array","items":{"type":"array","items":{"type":"number"}},"description":"K 列等长序列 (K>=2), 外层为变量列"},"eg_trend":{"type":"string","enum":["n","c","ct","ctt"],"default":"c"},"det_order":{"type":"integer","enum":[-1,0,1],"default":0,"description":"Johansen 确定项: -1 无 / 0 常数 / 1 线性趋势"},"k_ar_diff":{"type":"integer","minimum":0,"default":0}},"required":["series"],"additionalProperties":false})json",
         detail::handle_cointegration},
        {"cpphub_vecm_fit",
         "VECM (Johansen ML) 估计: alpha/beta + ECT t 统计量 (EM2002 5% 临界值) + "
         "PI=alpha*beta'。陷阱注释: ECT t 仅 rank=1 可用 (CI4 五情形 det)。",
         R"json({"type":"object","properties":{"series":{"type":"array","items":{"type":"array","items":{"type":"number"}},"description":"K 列等长水平序列 (K>=2)"},"rank":{"type":"integer","minimum":1,"description":"协整秩 r, 须在 [1,K]"},"k_ar_diff":{"type":"integer","minimum":0,"default":0},"det":{"type":"string","enum":["n","co","ci","lo","li"],"default":"n"},"urca_normalization":{"type":"boolean","default":false}},"required":["series","rank"],"additionalProperties":false})json",
         detail::handle_vecm},
        {"cpphub_garch_diagnose",
         "GARCH 族诊断: 参数估计 + sandwich SE + LB(z)/LB(z^2)/JB 标准化残差诊断 + "
         "稳定性。陷阱注释: JB 的 p 值经 Bootstrap (G-ADR4)。",
         R"json({"type":"object","properties":{"series":{"type":"array","items":{"type":"number"},"description":"收益序列"},"model":{"type":"string","enum":["GARCH","EGARCH","GJR","GARCH-M"]},"dist":{"type":"string","enum":["normal","t","ged"],"default":"normal"},"form":{"type":"string","enum":["var","vol","log"],"default":"vol","description":"仅 GARCH-M 有效"}},"required":["series","model"],"additionalProperties":false})json",
         detail::handle_garch},
        {"cpphub_granger_causality",
         "Granger 因果: F/chi2/LR + Toda-Yamamoto 增广 Wald + HAC 稳健 Wald。"
         "陷阱注释: 参数顺序显式 (cause, effect), 与 statsmodels 列序相反 (GR6); "
         "除 HAC 版外统计量均非异方差稳健 (GR5)。",
         R"json({"type":"object","properties":{"cause":{"type":"array","items":{"type":"number"},"description":"原因序列 (进入 effect 方程的滞后项)"},"effect":{"type":"array","items":{"type":"number"},"description":"效应序列 (被解释变量)"},"lag":{"type":"integer","minimum":1,"default":1},"d_max":{"type":"integer","minimum":0,"default":0,"description":"TY 增广阶; 0=不做 TY"},"with_hac":{"type":"boolean","default":true},"hac_bandwidth":{"type":"integer","minimum":0,"default":0}},"required":["cause","effect"],"additionalProperties":false})json",
         detail::handle_granger},
        {"cpphub_dy_spillover",
         "Diebold-Yilmaz 溢出指数: TCI + 逐变量 TO/FROM/NET (全部 x100/K 口径), "
         "可选滚动窗口。陷阱注释: 滚动窗口 rolling 须 > 2*K。",
         R"json({"type":"object","properties":{"series":{"type":"array","items":{"type":"array","items":{"type":"number"}},"description":"K 列等长序列 (K>=2)"},"horizon":{"type":"integer","minimum":1,"default":10},"lag":{"type":"integer","minimum":0,"default":0,"description":"0=每窗口 IC 自动选阶"},"trend":{"type":"string","enum":["c","n"],"default":"c"},"rolling":{"type":"integer","minimum":0,"default":0,"description":"滚动窗口大小; 0=全样本静态; 须 > 2*K"}},"required":["series"],"additionalProperties":false})json",
         detail::handle_dy},
    };
    return kSpecs;
}

// 按名查找（不存在返回 nullptr）
inline const ToolSpec* find_tool(const std::string& name) {
    for (const ToolSpec& s : tool_specs()) {
        if (name == s.name) return &s;
    }
    return nullptr;
}

// 工具面规模（不变式 3：<= 10）
inline std::size_t tool_count() { return tool_specs().size(); }

}  // namespace c_api
}  // namespace v1
}  // namespace cpphub
