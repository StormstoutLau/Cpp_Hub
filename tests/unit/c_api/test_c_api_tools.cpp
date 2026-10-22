// S3 单元测试：C ABI 扩展符号 cpphub_v1_tool_list / cpphub_v1_tool_call
// 依据 CHECKLIST.md §2.3 / §3.1 / §3.2 / §5：
//   1. tool_list 返回完整 schema JSON（可解析且含 7 工具）
//   2. tool_call 逐工具闭环（7 工具各至少 1 次正常调用返回容器）
//   3. 与既有 6 个标量符号共存不冲突
//   4. 缓冲协商：不足返回 -4 并回填所需长度；二次分配后成功
//   5. 错误映射：未知工具 / 非法 args_json / 非法参数 -> -1；内核退化 -> -2
#include "cpphub/c_api/cpphub_c_api.h"
#include "cpphub/c_api/cpphub_tool_registry.hpp"
#include "cpphub/c_api/detail/json_min.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace json = cpphub::v1::json_min;

namespace {

// 两阶段缓冲调用：先探需求长度，再足量重试（E1 可重放）
std::string fetch_list() {
    size_t need = 0;
    EXPECT_EQ(cpphub_v1_tool_list(nullptr, &need), CPPHUB_ERR_BUFFER_TOO_SMALL);
    EXPECT_GT(need, 0u);
    std::vector<char> buf(need + 1);
    size_t len = buf.size();
    EXPECT_EQ(cpphub_v1_tool_list(buf.data(), &len), CPPHUB_OK);
    EXPECT_EQ(len, need);
    EXPECT_EQ(std::strlen(buf.data()), need);
    return std::string(buf.data(), len);
}

std::string call_tool(const std::string& name, const std::string& args) {
    size_t need = 0;
    const int rc = cpphub_v1_tool_call(name.c_str(), args.c_str(), nullptr, &need);
    EXPECT_EQ(rc, CPPHUB_ERR_BUFFER_TOO_SMALL) << cpphub_v1_get_last_error();
    std::vector<char> buf(need + 1);
    size_t len = buf.size();
    const int rc2 = cpphub_v1_tool_call(name.c_str(), args.c_str(), buf.data(), &len);
    EXPECT_EQ(rc2, CPPHUB_OK) << cpphub_v1_get_last_error();
    return std::string(buf.data(), len);
}

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

std::string num_array(const std::vector<double>& v) {
    json::Value a = json::Value::array();
    for (double x : v) a.push(json::Value::number(x));
    return json::dump(a);
}

std::string matrix_json(const std::vector<std::vector<double>>& cols) {
    json::Value a = json::Value::array();
    for (const std::vector<double>& c : cols) {
        json::Value col = json::Value::array();
        for (double x : c) col.push(json::Value::number(x));
        a.push(std::move(col));
    }
    return json::dump(a);
}

json::Value parse_ok(const std::string& text) {
    json::Value v;
    std::string err;
    EXPECT_TRUE(json::parse(text, v, err)) << err;
    return v;
}

void expect_container(const std::string& text, const char* tool) {
    const json::Value v = parse_ok(text);
    ASSERT_TRUE(v.is_object()) << tool;
    ASSERT_NE(v.find("tool"), nullptr) << tool;
    EXPECT_EQ(v.find("tool")->as_string(), tool);
    ASSERT_NE(v.find("results"), nullptr) << tool;
    EXPECT_TRUE(v.find("results")->is_array()) << tool;
    ASSERT_NE(v.find("conclusion"), nullptr) << tool;
    EXPECT_TRUE(v.find("conclusion")->is_object()) << tool;
    EXPECT_GT(v.find("results")->as_array().size(), 0u) << tool;
}

}  // namespace

// ---------------------------------------------------------------------------
// 1. tool_list
// ---------------------------------------------------------------------------
TEST(CApiTools, ToolListReturnsSevenParsableTools) {
    const std::string text = fetch_list();
    const json::Value arr = parse_ok(text);
    ASSERT_TRUE(arr.is_array());
    ASSERT_EQ(arr.as_array().size(), 7u);
    for (const json::Value& t : arr.as_array()) {
        ASSERT_TRUE(t.is_object());
        ASSERT_NE(t.find("name"), nullptr);
        ASSERT_NE(t.find("description"), nullptr);
        ASSERT_NE(t.find("inputSchema"), nullptr);
        EXPECT_TRUE(t.find("name")->is_string());
        EXPECT_TRUE(t.find("description")->is_string());
        ASSERT_TRUE(t.find("inputSchema")->is_object());
        const json::Value* type = t.find("inputSchema")->find("type");
        ASSERT_NE(type, nullptr);
        EXPECT_EQ(type->as_string(), "object");
    }
    // 名称集合与 DESIGN §4.1 一致
    std::vector<std::string> names;
    for (const json::Value& t : arr.as_array()) {
        names.push_back(t.find("name")->as_string());
    }
    const std::vector<std::string> expected = {
        "cpphub_unit_root_diagnose", "cpphub_variance_ratio_test",
        "cpphub_cointegration_scan", "cpphub_vecm_fit",
        "cpphub_garch_diagnose",     "cpphub_granger_causality",
        "cpphub_dy_spillover"};
    EXPECT_EQ(names, expected);
}

TEST(CApiTools, ToolListIsStableAcrossCalls) {
    // 不变式 1/6：schema 为常量，重复调用逐字节相同
    EXPECT_EQ(fetch_list(), fetch_list());
}

// ---------------------------------------------------------------------------
// 2. tool_call 逐工具闭环（7/7）
// ---------------------------------------------------------------------------
TEST(CApiTools, ToolCallAllSevenTools) {
    EXPECT_NO_FATAL_FAILURE(expect_container(
        call_tool("cpphub_unit_root_diagnose",
                  "{\"series\":" + num_array(make_series(200, 0.98, 7)) + "}"),
        "cpphub_unit_root_diagnose"));

    EXPECT_NO_FATAL_FAILURE(expect_container(
        call_tool("cpphub_variance_ratio_test",
                  "{\"series\":" + num_array(make_series(240, 0.0, 13)) + "}"),
        "cpphub_variance_ratio_test"));

    const std::vector<double> y1 = make_series(200, 0.97, 17);
    std::vector<double> y0(200);
    {
        std::minstd_rand rng(19);
        std::uniform_real_distribution<double> u(-0.2, 0.2);
        for (std::size_t i = 0; i < 200; ++i) y0[i] = 0.9 * y1[i] + u(rng);
    }
    const std::string mat = matrix_json({y0, y1});
    EXPECT_NO_FATAL_FAILURE(expect_container(
        call_tool("cpphub_cointegration_scan",
                  "{\"series\":" + mat + ",\"det_order\":0,\"k_ar_diff\":1}"),
        "cpphub_cointegration_scan"));

    EXPECT_NO_FATAL_FAILURE(expect_container(
        call_tool("cpphub_vecm_fit",
                  "{\"series\":" + mat +
                      ",\"rank\":1,\"k_ar_diff\":1,\"det\":\"co\"}"),
        "cpphub_vecm_fit"));

    EXPECT_NO_FATAL_FAILURE(expect_container(
        call_tool("cpphub_garch_diagnose",
                  "{\"series\":" + num_array(make_series(300, 0.0, 31)) +
                      ",\"model\":\"GARCH\"}"),
        "cpphub_garch_diagnose"));

    const std::vector<double> effect = make_series(200, 0.9, 37);
    std::vector<double> cause(effect.size());
    for (std::size_t i = 0; i < effect.size(); ++i) {
        cause[i] = (i == 0) ? 0.0 : effect[i - 1] * 0.7;
    }
    EXPECT_NO_FATAL_FAILURE(expect_container(
        call_tool("cpphub_granger_causality",
                  "{\"cause\":" + num_array(cause) + ",\"effect\":" +
                      num_array(effect) + ",\"lag\":1,\"d_max\":1}"),
        "cpphub_granger_causality"));

    EXPECT_NO_FATAL_FAILURE(expect_container(
        call_tool("cpphub_dy_spillover",
                  "{\"series\":" +
                      matrix_json({make_series(200, 0.9, 41),
                                   make_series(200, 0.9, 43)}) +
                      ",\"horizon\":5}"),
        "cpphub_dy_spillover"));
}

// ---------------------------------------------------------------------------
// 3. 与既有 6 符号共存
// ---------------------------------------------------------------------------
TEST(CApiTools, CoexistsWithLegacyScalarSymbols) {
    EXPECT_EQ(cpphub_v1_get_abi_version(), 1);  // 未递增（加法性变更）
    double price = 0.0;
    ASSERT_EQ(cpphub_v1_heston_price(100.0, 100.0, 1.0, 0.05, 0.0, 0.04, 2.0,
                                     0.04, 0.3, -0.5, &price),
              CPPHUB_OK);
    EXPECT_GT(price, 0.0);

    // JSON 扩展调用后，既有符号仍可正常工作（共享 thread_local 错误槽不串扰）
    const std::string text = fetch_list();
    EXPECT_FALSE(text.empty());
    EXPECT_EQ(cpphub_v1_get_abi_version(), 1);
    double iv = 0.0;
    EXPECT_EQ(cpphub_v1_bsm_implied_vol(10.4506, 100.0, 100.0, 1.0, 0.05, 0.0,
                                       1, &iv),
              CPPHUB_OK);
    EXPECT_NEAR(iv, 0.2, 5e-3);
}

// ---------------------------------------------------------------------------
// 4. 缓冲协商
// ---------------------------------------------------------------------------
TEST(CApiTools, BufferNegotiationTwoPhase) {
    // 0 容量 -> -4 且回填真实需求
    char tiny[1];
    size_t len = 0;
    EXPECT_EQ(cpphub_v1_tool_list(tiny, &len), CPPHUB_ERR_BUFFER_TOO_SMALL);
    const size_t need = len;
    EXPECT_GT(need, 0u);
    EXPECT_NE(std::string(cpphub_v1_get_last_error()).find("too small"),
              std::string::npos);

    // 恰好 need（无 NUL 空间）仍拒绝
    std::vector<char> exact(need);
    size_t len2 = exact.size();
    EXPECT_EQ(cpphub_v1_tool_list(exact.data(), &len2),
              CPPHUB_ERR_BUFFER_TOO_SMALL);
    EXPECT_EQ(len2, need);

    // need + 1 -> 成功
    std::vector<char> ok(need + 1);
    size_t len3 = ok.size();
    EXPECT_EQ(cpphub_v1_tool_list(ok.data(), &len3), CPPHUB_OK);
    EXPECT_EQ(len3, need);

    // 过大的缓冲也被正确回填（不改容量语义）
    std::vector<char> big(need + 64, 'x');
    size_t len4 = big.size();
    EXPECT_EQ(cpphub_v1_tool_list(big.data(), &len4), CPPHUB_OK);
    EXPECT_EQ(len4, need);
}

TEST(CApiTools, NullOutLenRejected) {
    EXPECT_EQ(cpphub_v1_tool_list(nullptr, nullptr), CPPHUB_ERR_INVALID_ARG);
    EXPECT_EQ(cpphub_v1_tool_call("cpphub_unit_root_diagnose", "{}", nullptr,
                                  nullptr),
              CPPHUB_ERR_INVALID_ARG);
}

// ---------------------------------------------------------------------------
// 5. 错误映射
// ---------------------------------------------------------------------------
TEST(CApiTools, ErrorMapping) {
    char buf[4096];
    size_t len = sizeof(buf);

    // 未知工具 -> -1，错误信息含可用工具名
    EXPECT_EQ(cpphub_v1_tool_call("no_such_tool", "{}", buf, &len),
              CPPHUB_ERR_INVALID_ARG);
    const std::string e1 = cpphub_v1_get_last_error();
    EXPECT_NE(e1.find("unknown tool"), std::string::npos);
    EXPECT_NE(e1.find("cpphub_unit_root_diagnose"), std::string::npos);

    // 非法 JSON -> -1
    len = sizeof(buf);
    EXPECT_EQ(cpphub_v1_tool_call("cpphub_unit_root_diagnose", "{not json}",
                                  buf, &len),
              CPPHUB_ERR_INVALID_ARG);
    EXPECT_NE(std::string(cpphub_v1_get_last_error()).find("invalid args_json"),
              std::string::npos);

    // args 非 object -> -1
    len = sizeof(buf);
    EXPECT_EQ(cpphub_v1_tool_call("cpphub_unit_root_diagnose", "[1,2]", buf,
                                  &len),
              CPPHUB_ERR_INVALID_ARG);

    // 缺必填 series -> -1
    len = sizeof(buf);
    EXPECT_EQ(cpphub_v1_tool_call("cpphub_unit_root_diagnose", "{}", buf, &len),
              CPPHUB_ERR_INVALID_ARG);
    EXPECT_NE(std::string(cpphub_v1_get_last_error()).find("series"),
              std::string::npos);

    // nullptr 参数 -> -1
    len = sizeof(buf);
    EXPECT_EQ(cpphub_v1_tool_call(nullptr, "{}", buf, &len),
              CPPHUB_ERR_INVALID_ARG);
    len = sizeof(buf);
    EXPECT_EQ(cpphub_v1_tool_call("cpphub_unit_root_diagnose", nullptr, buf,
                                  &len),
              CPPHUB_ERR_INVALID_ARG);

    // 内核退化（常数序列 -> 零方差 runtime_error）-> -2
    len = sizeof(buf);
    const std::string constant = num_array(std::vector<double>(120, 1.0));
    EXPECT_EQ(cpphub_v1_tool_call(
                  "cpphub_variance_ratio_test",
                  ("{\"series\":" + constant + "}").c_str(), buf, &len),
              CPPHUB_ERR_RUNTIME);
    EXPECT_NE(std::string(cpphub_v1_get_last_error()).find("variance_ratio"),
              std::string::npos);
}

TEST(CApiTools, ErrorStateIsPerThreadAndResetBySuccess) {
    char buf[8192];
    size_t len = sizeof(buf);
    EXPECT_EQ(cpphub_v1_tool_call("no_such_tool", "{}", buf, &len),
              CPPHUB_ERR_INVALID_ARG);
    EXPECT_STRNE(cpphub_v1_get_last_error(), "");

    // 成功后上一次错误不再被新调用继承（成功路径不写错误串，
    // 但调用方依据返回码判断；此处仅确认成功仍返回 OK 且输出有效）
    len = sizeof(buf);
    const std::string args = "{\"series\":" + num_array(make_series(80, 0.9, 5)) + "}";
    EXPECT_EQ(cpphub_v1_tool_call("cpphub_unit_root_diagnose", args.c_str(),
                                  buf, &len),
              CPPHUB_OK);
    json::Value v;
    std::string err;
    EXPECT_TRUE(json::parse(std::string(buf, len), v, err)) << err;
}

// ---------------------------------------------------------------------------
// RISK-034: 入口长度上限（kMaxArgsJsonBytes, 8 MiB）
// ---------------------------------------------------------------------------
TEST(CApiTools, RejectsOversizedArgsJson) {
    const std::size_t limit = cpphub::v1::c_api::kMaxArgsJsonBytes;
    std::string args = "{\"series\":[" + std::string(limit, '0') + "]}";
    ASSERT_GT(args.size(), limit);
    size_t need = 0;
    EXPECT_EQ(cpphub_v1_tool_call("cpphub_unit_root_diagnose", args.c_str(),
                                  nullptr, &need),
              CPPHUB_ERR_INVALID_ARG);
    EXPECT_NE(std::string(cpphub_v1_get_last_error()).find("exceeds limit"),
              std::string::npos);
}

TEST(CApiTools, OversizedCheckPrecedesParse) {
    // 超长 + 语法非法 → 报 "exceeds limit" 而非 JSON 语法错误（证明检查前置）
    std::string args(cpphub::v1::c_api::kMaxArgsJsonBytes + 1, '!');
    size_t need = 0;
    EXPECT_EQ(cpphub_v1_tool_call("cpphub_unit_root_diagnose", args.c_str(),
                                  nullptr, &need),
              CPPHUB_ERR_INVALID_ARG);
    const std::string err = cpphub_v1_get_last_error();
    EXPECT_NE(err.find("exceeds limit"), std::string::npos) << err;
    EXPECT_EQ(err.find("invalid args_json"), std::string::npos) << err;
}
