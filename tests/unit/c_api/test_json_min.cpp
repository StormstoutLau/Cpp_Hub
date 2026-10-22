// S1 单元测试：json_min.hpp（最小 JSON 读写）
// 依据 IMPLEMENTATION.md §3.2 / CHECKLIST.md §2.1：
//   1. double 往返位精度（IEEE-754 位级相同）
//   2. 嵌套对象/数组（结构与值一致）
//   3. 字符串转义 / UTF-8（\u 转义、代理对、多字节直通）
//   4. null / 布尔 / 空数组边界
//   5. 非法输入拒绝（严格 JSON 语法）
#include "cpphub/c_api/detail/json_min.hpp"

#include <gtest/gtest.h>

#include <bit>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace jm = cpphub::v1::json_min;

namespace {

std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

// 往返：dump 后逐字节回读，返回重解析得到的 double
double round_trip_number(double x, const char** text_out = nullptr) {
    static std::string text;
    text = jm::dump(jm::Value::number(x));
    if (text_out) *text_out = text.c_str();
    std::string err;
    jm::Value v;
    EXPECT_TRUE(jm::parse(text, v, err)) << "text=" << text << " err=" << err;
    EXPECT_TRUE(v.is_number()) << "text=" << text;
    return v.is_number() ? v.as_number() : 0.0;
}

// 往返一个字符串值
std::string round_trip_string(const std::string& s) {
    const std::string text = jm::dump(jm::Value::string(s));
    std::string err;
    jm::Value v;
    EXPECT_TRUE(jm::parse(text, v, err)) << "text=" << text << " err=" << err;
    EXPECT_TRUE(v.is_string()) << "text=" << text;
    return v.is_string() ? v.as_string() : std::string();
}

}  // namespace

// ---------------------------------------------------------------------------
// 1. double 往返位精度
// ---------------------------------------------------------------------------
TEST(JsonMinNumber, DoubleRoundTripIsBitExact) {
    const std::vector<double> cases = {
        0.0,
        -0.0,
        1.0,
        -1.0,
        0.1,                       // 不可精确表示的十进制
        2.0 / 3.0,
        1.0 / 3.0,
        1e-300,
        1e300,
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::min(),        // 最小正规数
        std::numeric_limits<double>::denorm_min(), // 最小非正规数
        3.14159265358979323846,
        -123456.789012345,
    };
    for (double x : cases) {
        EXPECT_EQ(bits(round_trip_number(x)), bits(x)) << "value=" << x;
    }
}

TEST(JsonMinNumber, ShortestRoundTripForm) {
    // 最短往返表示：整数型 double 不带小数点/指数
    const char* text = nullptr;
    round_trip_number(1.0, &text);
    EXPECT_EQ(std::string(text), "1");
    round_trip_number(-2.0, &text);
    EXPECT_EQ(std::string(text), "-2");
    // 0.1 的最短往返表示为 "0.1"（而非 0.1000000000000000055…）
    round_trip_number(0.1, &text);
    EXPECT_EQ(std::string(text), "0.1");
}

// 外部 oracle 校验：与 Python repr()（最短往返表示）逐字符比对。
// 独立于本实现——不依赖自证往返，直接锚定第三方格式化结果。
TEST(JsonMinNumber, MatchesPythonShortestReprOracle) {
    struct Case {
        double value;
        const char* python_repr;  // 由 python -c "print(repr(v))" 取得
    };
    const Case cases[] = {
        {0.1, "0.1"},
        {1.0 / 3.0, "0.3333333333333333"},
        {2.0 / 3.0, "0.6666666666666666"},
        {1e-300, "1e-300"},
        {1e300, "1e+300"},
        {3.14159265358979323846, "3.141592653589793"},
        {-123456.789012345, "-123456.789012345"},
        {5e-324, "5e-324"},
    };
    for (const Case& c : cases) {
        EXPECT_EQ(jm::dump(jm::Value::number(c.value)), std::string(c.python_repr))
            << "value=" << c.value;
        // oracle 输出亦须被本解析器逐位还原
        std::string err;
        jm::Value v;
        ASSERT_TRUE(jm::parse(c.python_repr, v, err)) << err;
        EXPECT_EQ(bits(v.as_number()), bits(c.value)) << "oracle=" << c.python_repr;
    }
    // 整值：JSON 数值无 int/float 之分，最短形式省略 ".0"（Python json.loads 亦可解析）
    EXPECT_EQ(jm::dump(jm::Value::number(1.0)), "1");
    EXPECT_EQ(jm::dump(jm::Value::number(-2.0)), "-2");
    EXPECT_EQ(jm::dump(jm::Value::number(0.0)), "0");
    EXPECT_EQ(jm::dump(jm::Value::number(-0.0)), "-0");
}

TEST(JsonMinNumber, NonFiniteSerializesAsNull) {
    // JSON 无法表示 NaN/Inf → 落为 null（信封 p_value 可缺失的既定约定）
    EXPECT_EQ(jm::dump(jm::Value::number(std::numeric_limits<double>::quiet_NaN())), "null");
    EXPECT_EQ(jm::dump(jm::Value::number(std::numeric_limits<double>::infinity())), "null");
    EXPECT_EQ(jm::dump(jm::Value::number(-std::numeric_limits<double>::infinity())), "null");
}

// ---------------------------------------------------------------------------
// 2. 嵌套对象 / 数组
// ---------------------------------------------------------------------------
TEST(JsonMinStructure, NestedObjectAndArray) {
    jm::Value root = jm::Value::object({
        {"tool", jm::Value::string("cpphub_unit_root_diagnose")},
        {"reject_null", jm::Value::boolean(false)},
        {"p_value", jm::Value::number(0.03125)},
        {"series", jm::Value::array({jm::Value::number(1.0), jm::Value::number(2.5)})},
        {"notes", jm::Value::array()},
        {"cv", jm::Value::object({{"adf", jm::Value::number(-2.86)}})},
        {"meta", jm::Value::null()},
    });

    const std::string text = jm::dump(root);
    std::string err;
    jm::Value back;
    ASSERT_TRUE(jm::parse(text, back, err)) << "text=" << text << " err=" << err;
    ASSERT_TRUE(back.is_object());

    const jm::Value* tool = back.find("tool");
    ASSERT_NE(tool, nullptr);
    ASSERT_TRUE(tool->is_string());
    EXPECT_EQ(tool->as_string(), "cpphub_unit_root_diagnose");

    const jm::Value* rn = back.find("reject_null");
    ASSERT_NE(rn, nullptr);
    ASSERT_TRUE(rn->is_bool());
    EXPECT_FALSE(rn->as_bool());

    const jm::Value* pv = back.find("p_value");
    ASSERT_NE(pv, nullptr);
    ASSERT_TRUE(pv->is_number());
    EXPECT_EQ(bits(pv->as_number()), bits(0.03125));

    const jm::Value* series = back.find("series");
    ASSERT_NE(series, nullptr);
    ASSERT_TRUE(series->is_array());
    ASSERT_EQ(series->as_array().size(), 2u);
    EXPECT_EQ(bits(series->as_array()[0].as_number()), bits(1.0));
    EXPECT_EQ(bits(series->as_array()[1].as_number()), bits(2.5));

    const jm::Value* notes = back.find("notes");
    ASSERT_NE(notes, nullptr);
    ASSERT_TRUE(notes->is_array());
    EXPECT_TRUE(notes->as_array().empty());

    const jm::Value* cv = back.find("cv");
    ASSERT_NE(cv, nullptr);
    ASSERT_TRUE(cv->is_object());
    const jm::Value* adf = cv->find("adf");
    ASSERT_NE(adf, nullptr);
    EXPECT_EQ(bits(adf->as_number()), bits(-2.86));

    const jm::Value* meta = back.find("meta");
    ASSERT_NE(meta, nullptr);
    EXPECT_TRUE(meta->is_null());

    // 未知键返回 nullptr
    EXPECT_EQ(back.find("no_such_key"), nullptr);

    // 规范化稳定性：dump → parse → dump 逐字节相同
    EXPECT_EQ(jm::dump(back), text);
}

TEST(JsonMinStructure, DeepNestingWithinLimit) {
    std::string text = "[[[1]]]";
    std::string err;
    jm::Value v;
    ASSERT_TRUE(jm::parse(text, v, err)) << err;
    ASSERT_TRUE(v.is_array());
    EXPECT_EQ(jm::dump(v), "[[[1]]]");
}

// ---------------------------------------------------------------------------
// 3. 字符串转义 / UTF-8
// ---------------------------------------------------------------------------
TEST(JsonMinString, EscapesRoundTrip) {
    const std::string raw = "quote:\" back:\\ slash:/ \b\f\n\r\t ctrl:\x01 end";
    EXPECT_EQ(round_trip_string(raw), raw);
}

TEST(JsonMinString, Utf8Passthrough) {
    const std::string raw = "中文 € → ✓";
    EXPECT_EQ(round_trip_string(raw), raw);
    // 非 ASCII 字节不经 \u 转义（保持 UTF-8 直通）
    const std::string dumped = jm::dump(jm::Value::string(raw));
    EXPECT_EQ(dumped.find("\\u"), std::string::npos);
}

TEST(JsonMinString, UnicodeEscapeDecoded) {
    std::string err;
    jm::Value v;
    // \u0041 → 'A'，\n → 换行
    ASSERT_TRUE(jm::parse("\"a\\u0041\\n\"", v, err)) << err;
    ASSERT_TRUE(v.is_string());
    EXPECT_EQ(v.as_string(), std::string("aA\n"));
    // \u4e2d → U+4E2D（中），UTF-8 为 E4 B8 AD
    ASSERT_TRUE(jm::parse("\"\\u4e2d\"", v, err)) << err;
    ASSERT_TRUE(v.is_string());
    EXPECT_EQ(v.as_string(), std::string("\xE4\xB8\xAD"));
    // 代理对 U+1F600 → F0 9F 98 80
    ASSERT_TRUE(jm::parse("\"\\uD83D\\uDE00\"", v, err)) << err;
    ASSERT_TRUE(v.is_string());
    EXPECT_EQ(v.as_string(), std::string("\xF0\x9F\x98\x80"));
}

TEST(JsonMinString, RejectsBadEscapes) {
    std::string err;
    jm::Value v;
    EXPECT_FALSE(jm::parse("\"\\x\"", v, err));       // 非法转义字符
    EXPECT_FALSE(jm::parse("\"\\u12\"", v, err));     // 十六进制位不足
    EXPECT_FALSE(jm::parse("\"\\uD83D\"", v, err));   // 孤立高代理
    EXPECT_FALSE(jm::parse("\"\\uDE00\"", v, err));   // 孤立低代理
    EXPECT_FALSE(jm::parse("\"unterminated", v, err));
    EXPECT_FALSE(jm::parse("\"raw\ncontrol\"", v, err));  // 未转义控制字符
}

// ---------------------------------------------------------------------------
// 4. null / 布尔 / 空数组边界
// ---------------------------------------------------------------------------
TEST(JsonMinEdge, NullBoolEmptyArray) {
    std::string err;
    jm::Value v;

    ASSERT_TRUE(jm::parse("null", v, err)) << err;
    EXPECT_TRUE(v.is_null());
    EXPECT_EQ(jm::dump(v), "null");

    ASSERT_TRUE(jm::parse("true", v, err)) << err;
    ASSERT_TRUE(v.is_bool());
    EXPECT_TRUE(v.as_bool());
    EXPECT_EQ(jm::dump(v), "true");

    ASSERT_TRUE(jm::parse("false", v, err)) << err;
    ASSERT_TRUE(v.is_bool());
    EXPECT_FALSE(v.as_bool());
    EXPECT_EQ(jm::dump(v), "false");

    ASSERT_TRUE(jm::parse("[]", v, err)) << err;
    ASSERT_TRUE(v.is_array());
    EXPECT_TRUE(v.as_array().empty());
    EXPECT_EQ(jm::dump(v), "[]");

    ASSERT_TRUE(jm::parse("{}", v, err)) << err;
    ASSERT_TRUE(v.is_object());
    EXPECT_TRUE(v.as_object().empty());
    EXPECT_EQ(jm::dump(v), "{}");

    // 空数组嵌在对象中
    ASSERT_TRUE(jm::parse("{\"a\":[]}", v, err)) << err;
    EXPECT_EQ(jm::dump(v), "{\"a\":[]}");

    // 单元素数组 / 空白容忍
    ASSERT_TRUE(jm::parse("  [ 1 ,  2 ]  ", v, err)) << err;
    EXPECT_EQ(jm::dump(v), "[1,2]");
}

// ---------------------------------------------------------------------------
// 5. 非法输入拒绝（严格 JSON 语法）
// ---------------------------------------------------------------------------
TEST(JsonMinParse, RejectsMalformed) {
    const std::vector<std::string> bad = {
        "",          // 空输入
        "   ",       // 仅空白
        "{",         // 对象未闭合
        "[",         // 数组未闭合
        "{\"a\":}",  // 缺值
        "[1,]",      // 尾随逗号
        "[,1]",      // 前导逗号
        "{\"a\":1,}", // 对象尾随逗号
        "tru",       // 字面量不完整
        "nul",
        "NULL",      // 关键字大小写敏感
        "01",        // 前导零
        "1.",        // 小数点后无数字
        ".5",        // 缺整数部分
        "+1",        // 不允许前导正号
        "1e",        // 指数无数字
        "1e+",       // 指数符号后无数字
        "\"a\"\"b\"", // 值后跟多余内容
        "{} []",     // 顶层多余 token
        "[1 2]",     // 元素间缺逗号
        "{a:1}",     // 键未加引号
        "{'a':1}",   // 单引号键
        "NaN",       // 非 JSON 字面量
        "Infinity",
    };
    for (const std::string& text : bad) {
        std::string err;
        jm::Value v;
        EXPECT_FALSE(jm::parse(text, v, err)) << "should reject: " << text;
        EXPECT_FALSE(err.empty()) << "error message empty for: " << text;
    }
}

TEST(JsonMinParse, RejectsDepthOverflow) {
    // 深度上限保护：超限输入被拒（不递归爆栈）
    std::string deep;
    for (int i = 0; i < 200; ++i) deep.push_back('[');
    deep += "1";
    for (int i = 0; i < 200; ++i) deep.push_back(']');
    std::string err;
    jm::Value v;
    EXPECT_FALSE(jm::parse(deep, v, err));
    EXPECT_FALSE(err.empty());
}

TEST(JsonMinParse, ErrorMessageOnFailureIsSet) {
    std::string err = "sentinel";
    jm::Value v;
    EXPECT_FALSE(jm::parse("{\"a\" 1}", v, err));
    EXPECT_NE(err, "sentinel");
    EXPECT_FALSE(err.empty());
    // 成功时清空错误信息
    ASSERT_TRUE(jm::parse("{\"a\":1}", v, err));
    EXPECT_TRUE(err.empty());
}

// ---------------------------------------------------------------------------
// RISK-034: 节点预算（kMaxParseNodes）
// 注意: 文件顶部 (:18) 已有 `namespace jm = cpphub::v1::json_min;`,
// 同一作用域重复声明 namespace alias = 编译错误 (Step 6 勘误 R-1), 此处不得再写。
// ---------------------------------------------------------------------------

TEST(JsonMinBudget, AcceptsAtBoundary) {
    // N 元素数组 = 1(数组) + N(元素) 节点 → kMaxParseNodes-1 个元素恰为 1e6 节点
    std::string text = "[";
    for (std::size_t i = 0; i < jm::kMaxParseNodes - 1; ++i) {
        if (i) text += ',';
        text += '0';
    }
    text += ']';
    std::string err;
    jm::Value v;
    ASSERT_TRUE(jm::parse(text, v, err)) << err;
    ASSERT_TRUE(v.is_array());
    EXPECT_EQ(v.as_array().size(), jm::kMaxParseNodes - 1);
}

TEST(JsonMinBudget, RejectsWideArrayOverBudget) {
    // kMaxParseNodes+2 个元素 → 1'000'003 节点, 超限拒绝且错误含 "node"
    std::string text = "[";
    for (std::size_t i = 0; i < jm::kMaxParseNodes + 2; ++i) {
        if (i) text += ',';
        text += '0';
    }
    text += ']';
    std::string err;
    jm::Value v;
    EXPECT_FALSE(jm::parse(text, v, err));
    EXPECT_NE(std::string(err).find("node"), std::string::npos) << err;
}

TEST(JsonMinBudget, DeepNestingStillLimitedByDepth) {
    // 回归锁: 节点预算不遮蔽深度限制（小节点数 + 深嵌套仍报 nesting depth）
    std::string deep;
    for (int i = 0; i < 200; ++i) deep.push_back('[');
    deep += "1";
    for (int i = 0; i < 200; ++i) deep.push_back(']');
    std::string err;
    jm::Value v;
    EXPECT_FALSE(jm::parse(deep, v, err));
    EXPECT_NE(std::string(err).find("nesting depth"), std::string::npos) << err;
}
