// S2 单元测试：cpphub_tool_registry.hpp（7 工具注册表）
// 依据 CHECKLIST.md §2.2：
//   1. 7 工具注册完整（名称集合 = DESIGN §4.1 七个）
//   2. 名称唯一 + 字符集合规（[A-Za-z0-9_.-]{1,128}）+ 工具面 <= 10（不变式 3）
//   3. schema 哈希快照锁定（防"工具投毒式"静默变更）
//   4. inputSchema 为合法 JSON object 且含 properties/required
//   5. 逐工具功能冒烟：handler 返回 ToolResult 容器（DESIGN §6.2）
#include "cpphub/c_api/cpphub_tool_registry.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <regex>
#include <set>
#include <string>
#include <vector>

namespace c_api = cpphub::v1::c_api;
namespace json = cpphub::v1::json_min;

namespace {

// ---- schema 哈希（FNV-1a 64 位，跨平台稳定）----
std::uint64_t fnv1a(const std::string& s,
                    std::uint64_t h = 1469598103934665603ULL) {
    for (const unsigned char c : s) {
        h ^= c;
        h *= 1099511628211ULL;
    }
    return h;
}

std::string spec_hash(const c_api::ToolSpec& s) {
    std::uint64_t h = fnv1a(std::string(s.name));
    h = fnv1a(std::string(s.description), h);
    h = fnv1a(std::string(s.input_schema), h);
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%016llx",
                  static_cast<unsigned long long>(h));
    return std::string(buf);
}

// ---- 确定性序列（minstd_rand 为全标准规定，跨平台逐位一致）----
std::vector<double> make_series(std::size_t n, double phi, unsigned seed) {
    std::minstd_rand rng(seed);
    std::uniform_real_distribution<double> u(-0.5, 0.5);
    std::vector<double> x(n);
    double v = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        v = phi * v + u(rng);
        x[i] = v;
    }
    return x;
}

json::Value num_array(const std::vector<double>& v) {
    json::Value a = json::Value::array();
    for (double x : v) a.push(json::Value::number(x));
    return a;
}

json::Value matrix_value(const std::vector<std::vector<double>>& cols) {
    json::Value a = json::Value::array();
    for (const std::vector<double>& c : cols) a.push(num_array(c));
    return a;
}

// 调用工具并返回 ToolResult；失败时经断言暴露错误
json::Value call_tool(const char* name, const json::Value& args) {
    const c_api::ToolSpec* spec = c_api::find_tool(name);
    EXPECT_NE(spec, nullptr) << name;
    if (spec == nullptr) return json::Value();
    const c_api::ToolOutcome out = spec->handler(args);
    EXPECT_TRUE(out.ok) << name << " error: " << out.error;
    return out.result;
}

// ToolResult 容器形状校验（DESIGN §6.2）+ 每行信封 7 字段（§6.1）
void check_container_shape(const json::Value& r, const char* tool) {
    ASSERT_TRUE(r.is_object()) << tool;
    const json::Value* t = r.find("tool");
    ASSERT_NE(t, nullptr) << tool;
    ASSERT_TRUE(t->is_string());
    EXPECT_EQ(t->as_string(), tool);

    const json::Value* results = r.find("results");
    ASSERT_NE(results, nullptr) << tool;
    ASSERT_TRUE(results->is_array()) << tool;
    EXPECT_FALSE(results->as_array().empty()) << tool;

    for (const json::Value& row : results->as_array()) {
        ASSERT_TRUE(row.is_object()) << tool;
        for (const char* f : {"tool", "method", "statistic", "p_value",
                              "cv_5pct", "reject_null", "notes"}) {
            EXPECT_NE(row.find(f), nullptr) << tool << " row missing field " << f;
        }
        const json::Value* m = row.find("method");
        ASSERT_TRUE(m->is_string());
        // statistic 可为 null（该检验失败或未计算），但不允许缺字段
    }

    const json::Value* c = r.find("conclusion");
    ASSERT_NE(c, nullptr) << tool;
    ASSERT_TRUE(c->is_object()) << tool;
    EXPECT_NE(c->find("primary_method"), nullptr) << tool;
    EXPECT_NE(c->find("reject_null"), nullptr) << tool;
    EXPECT_NE(c->find("n_total"), nullptr) << tool;
    EXPECT_TRUE(c->find("reject_null")->is_bool()) << tool;

    // 主检验行必须存在于 results 中
    const std::string primary = c->find("primary_method")->as_string();
    bool found = false;
    for (const json::Value& row : results->as_array()) {
        if (row.find("method")->as_string() == primary) found = true;
    }
    EXPECT_TRUE(found) << tool << " primary_method not in results: " << primary;

    // 规范化稳定：dump -> parse -> dump 逐字节相同
    const std::string text = json::dump(r);
    json::Value back;
    std::string err;
    ASSERT_TRUE(json::parse(text, back, err)) << tool << " : " << err;
    EXPECT_EQ(json::dump(back), text) << tool;
}

}  // namespace

// ---------------------------------------------------------------------------
// 1. 注册完整性
// ---------------------------------------------------------------------------
TEST(ToolRegistry, SevenToolsRegisteredInDesignOrder) {
    const std::vector<std::string> expected = {
        "cpphub_unit_root_diagnose", "cpphub_variance_ratio_test",
        "cpphub_cointegration_scan", "cpphub_vecm_fit",
        "cpphub_garch_diagnose",     "cpphub_granger_causality",
        "cpphub_dy_spillover"};
    std::vector<std::string> actual;
    for (const c_api::ToolSpec& s : c_api::tool_specs()) {
        actual.push_back(s.name);
    }
    EXPECT_EQ(actual, expected);
    EXPECT_EQ(c_api::tool_count(), 7u);
}

TEST(ToolRegistry, NamesUniqueCharsetCompliantAndToolFaceCapped) {
    const std::regex re("^[A-Za-z0-9_.-]{1,128}$");
    std::set<std::string> seen;
    for (const c_api::ToolSpec& s : c_api::tool_specs()) {
        EXPECT_TRUE(std::regex_match(s.name, re)) << s.name;
        EXPECT_TRUE(seen.insert(s.name).second) << "duplicate name: " << s.name;
        EXPECT_NE(s.handler, nullptr) << s.name;
        ASSERT_NE(s.description, nullptr) << s.name;
        EXPECT_GT(std::strlen(s.description), 0u) << s.name;
    }
    EXPECT_LE(c_api::tool_count(), 10u);  // 不变式 3
    EXPECT_EQ(c_api::find_tool("no_such_tool"), nullptr);
}

// ---------------------------------------------------------------------------
// 2. schema 快照锁定（防投毒式变更）
// ---------------------------------------------------------------------------
TEST(ToolRegistry, SchemaHashSnapshotLocked) {
    // 快照 = FNV-1a64(name || description || inputSchema)
    // 变更工具名/描述/schema 必须显式更新此表（即"投毒式静默变更"不可通过）
    const std::vector<std::pair<std::string, std::string>> expected = {
        {"cpphub_unit_root_diagnose", "e14c201b5ae55741"},
        {"cpphub_variance_ratio_test", "038c3fcd32181479"},
        {"cpphub_cointegration_scan", "0ab5c133754fe8b3"},
        {"cpphub_vecm_fit", "5dd27ed7aa72a0fd"},
        {"cpphub_garch_diagnose", "d457fe4cef29fe9c"},
        {"cpphub_granger_causality", "f956c7c3a80911ac"},
        {"cpphub_dy_spillover", "9c2f7816784078c5"},
    };
    ASSERT_EQ(c_api::tool_specs().size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const c_api::ToolSpec& s = c_api::tool_specs()[i];
        EXPECT_EQ(s.name, expected[i].first);
        EXPECT_EQ(spec_hash(s), expected[i].second) << "tool=" << s.name;
    }
}

TEST(ToolRegistry, InputSchemasAreValidJsonObjects) {
    for (const c_api::ToolSpec& s : c_api::tool_specs()) {
        std::string err;
        json::Value v;
        ASSERT_TRUE(json::parse(s.input_schema, v, err)) << s.name << ": " << err;
        ASSERT_TRUE(v.is_object()) << s.name;
        const json::Value* t = v.find("type");
        ASSERT_NE(t, nullptr) << s.name;
        EXPECT_EQ(t->as_string(), "object") << s.name;
        EXPECT_NE(v.find("properties"), nullptr) << s.name;
        EXPECT_TRUE(v.find("properties")->is_object()) << s.name;
        const json::Value* req = v.find("required");
        ASSERT_NE(req, nullptr) << s.name;
        EXPECT_TRUE(req->is_array()) << s.name;
        // 陷阱注释须随 description 下发（DESIGN §4.1 末段）
        EXPECT_GT(strlen(s.description), 30u) << s.name;
    }
}

// ---------------------------------------------------------------------------
// 3. 逐工具功能冒烟（容器形状 + 规范化稳定）
// ---------------------------------------------------------------------------
TEST(ToolRegistry, UnitRootDiagnoseSmoke) {
    json::Value args = json::Value::object();
    args.set("series", num_array(make_series(200, 0.98, 7)));
    args.set("trend", json::Value::string("ct"));
    const json::Value r = call_tool("cpphub_unit_root_diagnose", args);
    check_container_shape(r, "cpphub_unit_root_diagnose");
    // 默认 6 个检验各一行
    EXPECT_EQ(r.find("results")->as_array().size(), 6u);
    // 主检验 = adf（DESIGN §4.1）
    EXPECT_EQ(r.find("conclusion")->find("primary_method")->as_string(), "adf");
    EXPECT_EQ(r.find("conclusion")->find("n_total")->as_number(), 6.0);
    // KPSS 行必须携带 U12 方向陷阱注释（消费方免疫）
    for (const json::Value& row : r.find("results")->as_array()) {
        if (row.find("method")->as_string() == "kpss") {
            bool has_trap = false;
            for (const json::Value& n : row.find("notes")->as_array()) {
                if (n.as_string().find("U12") != std::string::npos) has_trap = true;
            }
            EXPECT_TRUE(has_trap) << "KPSS row must carry U12 trap note";
        }
    }
}

TEST(ToolRegistry, UnitRootShortSeriesDegradesPerRow) {
    // T=8：ADF/PP/KPSS/DF-GLS 可跑，Ng-Perron/ZA 需 T>=10 → 该行 statistic=null 且有 failed 注
    json::Value args = json::Value::object();
    args.set("series", num_array(make_series(8, 0.9, 11)));
    const json::Value r = call_tool("cpphub_unit_root_diagnose", args);
    ASSERT_TRUE(r.is_object());
    EXPECT_EQ(r.find("results")->as_array().size(), 6u);
    std::size_t null_stat = 0;
    for (const json::Value& row : r.find("results")->as_array()) {
        if (row.find("statistic")->is_null()) ++null_stat;
    }
    EXPECT_GE(null_stat, 2u) << "Ng-Perron / ZA should degrade (T<10)";
}

TEST(ToolRegistry, VarianceRatioSmoke) {
    json::Value args = json::Value::object();
    args.set("series", num_array(make_series(240, 0.0, 13)));
    json::Value k = json::Value::array();
    k.push(json::Value::number(2));
    k.push(json::Value::number(5));
    args.set("k", std::move(k));
    const json::Value r = call_tool("cpphub_variance_ratio_test", args);
    check_container_shape(r, "cpphub_variance_ratio_test");
    EXPECT_EQ(r.find("conclusion")->find("primary_method")->as_string(),
              "lo_mackinlay_z2");
    EXPECT_EQ(r.find("results")->as_array().size(), 4u);
}

TEST(ToolRegistry, CointegrationScanSmoke) {
    const std::vector<double> y1 = make_series(200, 0.97, 17);
    std::vector<double> y0(200);
    std::minstd_rand rng(19);
    std::uniform_real_distribution<double> u(-0.2, 0.2);
    for (std::size_t i = 0; i < 200; ++i) y0[i] = 0.9 * y1[i] + u(rng);

    json::Value args = json::Value::object();
    args.set("series", matrix_value({y0, y1}));
    args.set("det_order", json::Value::number(0));
    args.set("k_ar_diff", json::Value::number(1));
    const json::Value r = call_tool("cpphub_cointegration_scan", args);
    check_container_shape(r, "cpphub_cointegration_scan");
    // K=2：EG 1 行 + Johansen 迹/最大特征值 各 2 行 = 5 行
    EXPECT_EQ(r.find("results")->as_array().size(), 5u);
    EXPECT_EQ(r.find("conclusion")->find("primary_method")->as_string(),
              "engle_granger");
}

TEST(ToolRegistry, VECMFitSmoke) {
    const std::vector<double> y1 = make_series(200, 0.97, 23);
    std::vector<double> y0(200);
    std::minstd_rand rng(29);
    std::uniform_real_distribution<double> u(-0.2, 0.2);
    for (std::size_t i = 0; i < 200; ++i) y0[i] = 0.9 * y1[i] + u(rng);

    json::Value args = json::Value::object();
    args.set("series", matrix_value({y0, y1}));
    args.set("rank", json::Value::number(1));
    args.set("k_ar_diff", json::Value::number(1));
    args.set("det", json::Value::string("co"));
    const json::Value r = call_tool("cpphub_vecm_fit", args);
    check_container_shape(r, "cpphub_vecm_fit");
    ASSERT_GT(r.find("results")->as_array().size(), 1u);
    EXPECT_EQ(r.find("conclusion")->find("primary_method")->as_string(), "ect_t_0");
}

TEST(ToolRegistry, GarchDiagnoseSmoke) {
    json::Value args = json::Value::object();
    args.set("series", num_array(make_series(300, 0.0, 31)));
    args.set("model", json::Value::string("GARCH"));
    const json::Value r = call_tool("cpphub_garch_diagnose", args);
    check_container_shape(r, "cpphub_garch_diagnose");
    bool has_lb = false, has_jb = false;
    for (const json::Value& row : r.find("results")->as_array()) {
        const std::string m = row.find("method")->as_string();
        if (m == "ljung_box_z_squared") has_lb = true;
        if (m == "jarque_bera") has_jb = true;
    }
    EXPECT_TRUE(has_lb);
    EXPECT_TRUE(has_jb);
}

TEST(ToolRegistry, GrangerCausalitySmoke) {
    const std::vector<double> effect = make_series(200, 0.9, 37);
    std::vector<double> cause(effect.size());
    for (std::size_t i = 0; i < effect.size(); ++i) {
        cause[i] = (i == 0) ? 0.0 : effect[i - 1] * 0.7;
    }
    json::Value args = json::Value::object();
    args.set("cause", num_array(cause));
    args.set("effect", num_array(effect));
    args.set("lag", json::Value::number(1));
    args.set("d_max", json::Value::number(1));
    const json::Value r = call_tool("cpphub_granger_causality", args);
    check_container_shape(r, "cpphub_granger_causality");
    EXPECT_EQ(r.find("conclusion")->find("primary_method")->as_string(), "ssr_f");
    // TY 与 HAC 行须存在
    std::set<std::string> methods;
    for (const json::Value& row : r.find("results")->as_array()) {
        methods.insert(row.find("method")->as_string());
    }
    EXPECT_TRUE(methods.count("ty_wald") == 1u);
    EXPECT_TRUE(methods.count("hac_wald") == 1u);
}

TEST(ToolRegistry, DYSpilloverSmoke) {
    const std::vector<double> a = make_series(200, 0.9, 41);
    const std::vector<double> b = make_series(200, 0.9, 43);
    json::Value args = json::Value::object();
    args.set("series", matrix_value({a, b}));
    args.set("horizon", json::Value::number(5));
    const json::Value r = call_tool("cpphub_dy_spillover", args);
    check_container_shape(r, "cpphub_dy_spillover");
    EXPECT_EQ(r.find("conclusion")->find("primary_method")->as_string(), "tci");
    // tci 1 行 + 每变量 TO/FROM/NET 各 1 行 = 1 + 3*2 = 7
    EXPECT_EQ(r.find("results")->as_array().size(), 7u);
}

// ---------------------------------------------------------------------------
// 4. 输入校验（CHECKLIST §2.5：超长/非有限/类型错误）
// ---------------------------------------------------------------------------
TEST(ToolRegistry, RejectsInvalidInput) {
    const c_api::ToolSpec* s = c_api::find_tool("cpphub_unit_root_diagnose");
    ASSERT_NE(s, nullptr);

    // 缺必填 series
    {
        const c_api::ToolOutcome o = s->handler(json::Value::object());
        EXPECT_FALSE(o.ok);
        EXPECT_FALSE(o.error.empty());
    }
    // 类型错误：series 不是数组
    {
        json::Value a = json::Value::object();
        a.set("series", json::Value::string("not an array"));
        const c_api::ToolOutcome o = s->handler(a);
        EXPECT_FALSE(o.ok);
    }
    // 非有限值
    {
        json::Value a = json::Value::object();
        json::Value arr = json::Value::array();
        arr.push(json::Value::number(1.0));
        arr.push(json::Value::number(std::numeric_limits<double>::quiet_NaN()));
        a.set("series", std::move(arr));
        const c_api::ToolOutcome o = s->handler(a);
        EXPECT_FALSE(o.ok);
        EXPECT_NE(o.error.find("non-finite"), std::string::npos);
    }
    // 超长（> kMaxSeriesLength）
    {
        json::Value a = json::Value::object();
        json::Value arr = json::Value::array();
        for (std::size_t i = 0; i < c_api::kMaxSeriesLength + 1; ++i) {
            arr.push(json::Value::number(1.0));
        }
        a.set("series", std::move(arr));
        const c_api::ToolOutcome o = s->handler(a);
        EXPECT_FALSE(o.ok);
        EXPECT_NE(o.error.find("max length"), std::string::npos);
    }
    // 非法 trend
    {
        json::Value a = json::Value::object();
        a.set("series", num_array(make_series(50, 0.9, 3)));
        a.set("trend", json::Value::string("nonsense"));
        const c_api::ToolOutcome o = s->handler(a);
        EXPECT_FALSE(o.ok);
    }
}

TEST(ToolRegistry, ToolCallIsStatelessAndRepeatable) {
    // 不变式 1：同一请求连续两次调用结果逐字节相同
    json::Value args = json::Value::object();
    args.set("series", num_array(make_series(120, 0.95, 5)));
    const c_api::ToolSpec* s = c_api::find_tool("cpphub_unit_root_diagnose");
    ASSERT_NE(s, nullptr);
    const c_api::ToolOutcome a1 = s->handler(args);
    const c_api::ToolOutcome a2 = s->handler(args);
    ASSERT_TRUE(a1.ok);
    ASSERT_TRUE(a2.ok);
    EXPECT_EQ(json::dump(a1.result), json::dump(a2.result));
}
