// 最小 JSON 读写（无第三方依赖）
//
// 依据 spec/mcp-server/IMPLEMENTATION.md §3.2（U-6 定案）：
//   - 范围限定：仅支持本 feature 所需子集（对象 / 数组 / 字符串 / 数字 / 布尔 / null）
//   - 不做通用 JSON 库；header-only；解析单遍无回溯
//   - double 采用最短往返表示（std::to_chars / std::from_chars），保证 IEEE-754 位精度往返
//   - 非有限值（NaN/±Inf）JSON 无法表示 → 序列化为 null
//
// 使用方：cpphub_v1_tool_list / cpphub_v1_tool_call（src/c_api/c_api.cpp）。
#ifndef CPPHUB_C_API_DETAIL_JSON_MIN_HPP
#define CPPHUB_C_API_DETAIL_JSON_MIN_HPP

#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace cpphub {
inline namespace v1 {
namespace json_min {

// 解析深度上限：递归下降解析器对深层嵌套输入的保护（避免栈溢出）。
inline constexpr int kMaxParseDepth = 64;

// JSON 值：可判别联合（duck-typed 树）。
class Value {
public:
    enum class Type { Null, Bool, Number, String, Array, Object };

    using Array = std::vector<Value>;
    using Member = std::pair<std::string, Value>;
    using Object = std::vector<Member>;

    Value() = default;  // 默认 null

    static Value null() { return Value(); }
    static Value boolean(bool b) {
        Value v;
        v.type_ = Type::Bool;
        v.bool_ = b;
        return v;
    }
    static Value number(double x) {
        Value v;
        v.type_ = Type::Number;
        v.num_ = x;
        return v;
    }
    static Value string(std::string s) {
        Value v;
        v.type_ = Type::String;
        v.str_ = std::move(s);
        return v;
    }
    static Value array(Array a = {}) {
        Value v;
        v.type_ = Type::Array;
        v.arr_ = std::move(a);
        return v;
    }
    static Value object(Object o = {}) {
        Value v;
        v.type_ = Type::Object;
        v.obj_ = std::move(o);
        return v;
    }

    Type type() const noexcept { return type_; }
    bool is_null() const noexcept { return type_ == Type::Null; }
    bool is_bool() const noexcept { return type_ == Type::Bool; }
    bool is_number() const noexcept { return type_ == Type::Number; }
    bool is_string() const noexcept { return type_ == Type::String; }
    bool is_array() const noexcept { return type_ == Type::Array; }
    bool is_object() const noexcept { return type_ == Type::Object; }

    bool as_bool() const noexcept { return bool_; }
    double as_number() const noexcept { return num_; }
    const std::string& as_string() const noexcept { return str_; }
    const Array& as_array() const noexcept { return arr_; }
    const Object& as_object() const noexcept { return obj_; }
    Array& as_array() noexcept { return arr_; }
    Object& as_object() noexcept { return obj_; }

    // 对象成员查找（不存在返回 nullptr）。
    const Value* find(std::string_view key) const noexcept {
        for (const Member& m : obj_) {
            if (m.first == key) return &m.second;
        }
        return nullptr;
    }

    // 追加对象成员（保持插入顺序，不覆盖同名键）。
    void set(std::string key, Value v) { obj_.emplace_back(std::move(key), std::move(v)); }
    // 追加数组元素。
    void push(Value v) { arr_.push_back(std::move(v)); }

private:
    Type type_ = Type::Null;
    bool bool_ = false;
    double num_ = 0.0;
    std::string str_;
    Array arr_;
    Object obj_;
};

namespace detail {

inline bool is_digit(char c) noexcept { return c >= '0' && c <= '9'; }

// 码点 → UTF-8 追加。
inline void append_utf8(std::string& out, std::uint32_t cp) {
    if (cp <= 0x7Fu) {
        out.push_back(static_cast<char>(cp));
    } else if (cp <= 0x7FFu) {
        out.push_back(static_cast<char>(0xC0u | (cp >> 6)));
        out.push_back(static_cast<char>(0x80u | (cp & 0x3Fu)));
    } else if (cp <= 0xFFFFu) {
        out.push_back(static_cast<char>(0xE0u | (cp >> 12)));
        out.push_back(static_cast<char>(0x80u | ((cp >> 6) & 0x3Fu)));
        out.push_back(static_cast<char>(0x80u | (cp & 0x3Fu)));
    } else {
        out.push_back(static_cast<char>(0xF0u | (cp >> 18)));
        out.push_back(static_cast<char>(0x80u | ((cp >> 12) & 0x3Fu)));
        out.push_back(static_cast<char>(0x80u | ((cp >> 6) & 0x3Fu)));
        out.push_back(static_cast<char>(0x80u | (cp & 0x3Fu)));
    }
}

// 读取 4 位十六进制（\uXXXX）；不足或非法返回 false。
inline bool hex4(std::string_view s, std::size_t& i, std::uint32_t& out) noexcept {
    if (i + 4 > s.size()) return false;
    std::uint32_t v = 0;
    for (int k = 0; k < 4; ++k) {
        const char c = s[i + static_cast<std::size_t>(k)];
        std::uint32_t d = 0;
        if (c >= '0' && c <= '9') {
            d = static_cast<std::uint32_t>(c - '0');
        } else if (c >= 'a' && c <= 'f') {
            d = static_cast<std::uint32_t>(c - 'a' + 10);
        } else if (c >= 'A' && c <= 'F') {
            d = static_cast<std::uint32_t>(c - 'A' + 10);
        } else {
            return false;
        }
        v = (v << 4) | d;
    }
    i += 4;
    out = v;
    return true;
}

// 单遍递归下降解析器。
class Parser {
public:
    Parser(std::string_view text, std::string& err) : s_(text), err_(err) {}

    bool parse_document(Value& out) {
        skip_ws();
        Value v;
        if (!parse_value(v, 0)) return false;
        skip_ws();
        if (pos_ != s_.size()) return fail("trailing characters after JSON value");
        out = std::move(v);
        return true;
    }

private:
    std::string_view s_;
    std::size_t pos_ = 0;
    std::string& err_;

    bool fail(const char* msg) {
        err_ = msg;
        return false;
    }
    char peek() const noexcept { return pos_ < s_.size() ? s_[pos_] : '\0'; }
    void skip_ws() noexcept {
        while (pos_ < s_.size()) {
            const char c = s_[pos_];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
                ++pos_;
            } else {
                break;
            }
        }
    }

    bool parse_value(Value& out, int depth) {
        switch (peek()) {
            case '{':
                return parse_object(out, depth);
            case '[':
                return parse_array(out, depth);
            case '"': {
                std::string str;
                if (!parse_string(str)) return false;
                out = Value::string(std::move(str));
                return true;
            }
            case 't':
                if (s_.substr(pos_, 4) == "true") {
                    pos_ += 4;
                    out = Value::boolean(true);
                    return true;
                }
                return fail("invalid literal (expected true)");
            case 'f':
                if (s_.substr(pos_, 5) == "false") {
                    pos_ += 5;
                    out = Value::boolean(false);
                    return true;
                }
                return fail("invalid literal (expected false)");
            case 'n':
                if (s_.substr(pos_, 4) == "null") {
                    pos_ += 4;
                    out = Value::null();
                    return true;
                }
                return fail("invalid literal (expected null)");
            case '\0':
                return fail("unexpected end of input");
            default:
                return parse_number(out);
        }
    }

    bool parse_object(Value& out, int depth) {
        if (depth >= kMaxParseDepth) return fail("maximum nesting depth exceeded");
        ++pos_;  // '{'
        Value::Object obj;
        skip_ws();
        if (peek() == '}') {
            ++pos_;
            out = Value::object(std::move(obj));
            return true;
        }
        for (;;) {
            skip_ws();
            if (peek() != '"') return fail("expected string key in object");
            std::string key;
            if (!parse_string(key)) return false;
            skip_ws();
            if (peek() != ':') return fail("expected ':' after object key");
            ++pos_;
            skip_ws();
            Value val;
            if (!parse_value(val, depth + 1)) return false;
            obj.emplace_back(std::move(key), std::move(val));
            skip_ws();
            const char c = peek();
            if (c == ',') {
                ++pos_;
                continue;
            }
            if (c == '}') {
                ++pos_;
                out = Value::object(std::move(obj));
                return true;
            }
            return fail("expected ',' or '}' in object");
        }
    }

    bool parse_array(Value& out, int depth) {
        if (depth >= kMaxParseDepth) return fail("maximum nesting depth exceeded");
        ++pos_;  // '['
        Value::Array arr;
        skip_ws();
        if (peek() == ']') {
            ++pos_;
            out = Value::array(std::move(arr));
            return true;
        }
        for (;;) {
            skip_ws();
            Value val;
            if (!parse_value(val, depth + 1)) return false;
            arr.push_back(std::move(val));
            skip_ws();
            const char c = peek();
            if (c == ',') {
                ++pos_;
                continue;
            }
            if (c == ']') {
                ++pos_;
                out = Value::array(std::move(arr));
                return true;
            }
            return fail("expected ',' or ']' in array");
        }
    }

    bool parse_string(std::string& out) {
        if (peek() != '"') return fail("expected string");
        ++pos_;
        out.clear();
        for (;;) {
            if (pos_ >= s_.size()) return fail("unterminated string");
            const unsigned char c = static_cast<unsigned char>(s_[pos_]);
            if (c == '"') {
                ++pos_;
                return true;
            }
            if (c == '\\') {
                ++pos_;
                if (pos_ >= s_.size()) return fail("unterminated escape sequence");
                switch (s_[pos_]) {
                    case '"': out.push_back('"'); ++pos_; break;
                    case '\\': out.push_back('\\'); ++pos_; break;
                    case '/': out.push_back('/'); ++pos_; break;
                    case 'b': out.push_back('\b'); ++pos_; break;
                    case 'f': out.push_back('\f'); ++pos_; break;
                    case 'n': out.push_back('\n'); ++pos_; break;
                    case 'r': out.push_back('\r'); ++pos_; break;
                    case 't': out.push_back('\t'); ++pos_; break;
                    case 'u': {
                        ++pos_;
                        std::uint32_t cp = 0;
                        if (!hex4(s_, pos_, cp)) return fail("invalid \\u escape");
                        if (cp >= 0xD800u && cp <= 0xDBFFu) {
                            // 高代理：必须跟随 \uDC00-\uDFFF 低代理
                            if (pos_ + 1 >= s_.size() || s_[pos_] != '\\' || s_[pos_ + 1] != 'u') {
                                return fail("lone high surrogate");
                            }
                            pos_ += 2;
                            std::uint32_t lo = 0;
                            if (!hex4(s_, pos_, lo)) return fail("invalid \\u escape");
                            if (lo < 0xDC00u || lo > 0xDFFFu) return fail("invalid low surrogate");
                            cp = 0x10000u + ((cp - 0xD800u) << 10) + (lo - 0xDC00u);
                        } else if (cp >= 0xDC00u && cp <= 0xDFFFu) {
                            return fail("lone low surrogate");
                        }
                        append_utf8(out, cp);
                        break;
                    }
                    default:
                        return fail("invalid escape character");
                }
                continue;
            }
            if (c < 0x20u) return fail("unescaped control character in string");
            out.push_back(s_[pos_]);
            ++pos_;
        }
    }

    bool parse_number(Value& out) {
        const std::size_t start = pos_;
        if (peek() == '-') ++pos_;
        if (pos_ >= s_.size()) return fail("unexpected end of number");
        const char c = s_[pos_];
        if (c == '0') {
            ++pos_;
        } else if (c >= '1' && c <= '9') {
            while (pos_ < s_.size() && is_digit(s_[pos_])) ++pos_;
        } else {
            return fail("invalid number");
        }
        if (pos_ < s_.size() && s_[pos_] == '.') {
            ++pos_;
            if (pos_ >= s_.size() || !is_digit(s_[pos_])) {
                return fail("invalid number: missing fraction digits");
            }
            while (pos_ < s_.size() && is_digit(s_[pos_])) ++pos_;
        }
        if (pos_ < s_.size() && (s_[pos_] == 'e' || s_[pos_] == 'E')) {
            ++pos_;
            if (pos_ < s_.size() && (s_[pos_] == '+' || s_[pos_] == '-')) ++pos_;
            if (pos_ >= s_.size() || !is_digit(s_[pos_])) {
                return fail("invalid number: missing exponent digits");
            }
            while (pos_ < s_.size() && is_digit(s_[pos_])) ++pos_;
        }
        // 语法已校验；由此处做一次正确舍入的十进制 → double 转换
        const char* const first = s_.data() + start;
        const char* const last = s_.data() + pos_;
        double value = 0.0;
        const std::from_chars_result res = std::from_chars(first, last, value);
        if (res.ec != std::errc() || res.ptr != last) {
            return fail("invalid number: out of range or unparsable");
        }
        out = Value::number(value);
        return true;
    }
};

inline void dump_string(std::string& out, const std::string& s) {
    out.push_back('"');
    for (const unsigned char c : s) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20u) {
                    static const char* const kHex = "0123456789abcdef";
                    out += "\\u00";
                    out.push_back(kHex[(c >> 4) & 0x0Fu]);
                    out.push_back(kHex[c & 0x0Fu]);
                } else {
                    // 含非 ASCII：UTF-8 字节直通（不做 \u 转义）
                    out.push_back(static_cast<char>(c));
                }
        }
    }
    out.push_back('"');
}

inline void dump_number(std::string& out, double x) {
    if (!std::isfinite(x)) {
        // NaN / ±Inf 无 JSON 表示 → null
        out += "null";
        return;
    }
    char buf[32];
    const std::to_chars_result res = std::to_chars(buf, buf + sizeof(buf), x);
    out.append(buf, res.ptr);
}

inline void dump_value(std::string& out, const Value& v) {
    switch (v.type()) {
        case Value::Type::Null:
            out += "null";
            break;
        case Value::Type::Bool:
            out += (v.as_bool() ? "true" : "false");
            break;
        case Value::Type::Number:
            dump_number(out, v.as_number());
            break;
        case Value::Type::String:
            dump_string(out, v.as_string());
            break;
        case Value::Type::Array: {
            const Value::Array& a = v.as_array();
            out.push_back('[');
            for (std::size_t i = 0; i < a.size(); ++i) {
                if (i != 0) out.push_back(',');
                dump_value(out, a[i]);
            }
            out.push_back(']');
            break;
        }
        case Value::Type::Object: {
            const Value::Object& o = v.as_object();
            out.push_back('{');
            for (std::size_t i = 0; i < o.size(); ++i) {
                if (i != 0) out.push_back(',');
                dump_string(out, o[i].first);
                out.push_back(':');
                dump_value(out, o[i].second);
            }
            out.push_back('}');
            break;
        }
    }
}

}  // namespace detail

// 解析 JSON 文本。成功返回 true 并写入 `out`；
// 失败返回 false，`error` 置为诊断信息（成功时清空）。
inline bool parse(std::string_view text, Value& out, std::string& error) {
    error.clear();
    detail::Parser parser(text, error);
    if (!parser.parse_document(out)) {
        out = Value();
        return false;
    }
    return true;
}

// 序列化 JSON 值（紧凑形式，无多余空白）。
inline std::string dump(const Value& v) {
    std::string out;
    out.reserve(64);
    detail::dump_value(out, v);
    return out;
}

}  // namespace json_min
}  // namespace v1
}  // namespace cpphub

#endif  // CPPHUB_C_API_DETAIL_JSON_MIN_HPP
