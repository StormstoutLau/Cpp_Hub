// =============================================================================
// test_result_base.hpp - 通用检验结果基结构
//
// Phase 7A Wave 0: detail/ 公共基础设施 (ADR-015 方案 B)
//
// 设计依据: R htest S3 类 (statistic/p.value/method 三核心字段)
// 复用方式: 组合 (非继承), 避免 vtable 开销, 保持 POD 友好
// 例外: 复合诊断 (VolatilityDiagnosticsResult 等) 不组合 base (ADR-015 决策点 3)
//
// 教材锚点: ADR-015 §4.1 (R htest 字段核实)
// =============================================================================
#pragma once

#include <string>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {
namespace detail {

// 通用检验结果基结构 (组合方式复用)
// 对应 R htest S3 类的 statistic/p.value/method 三个核心字段
// reject_null 为 Cpp_Hub 扩展字段, 便于业务层快速判断
struct TestResultBase {
    Real statistic;           // 检验统计量 (JB/LB/CD/DM/...)
    Real p_value;             // p 值 (无 p 值时为 NaN, 如 Cragg-Donald 非标准分布)
    std::string method_name;  // 方法名 ("Jarque-Bera"/"Ljung-Box"/...)
    bool reject_null;         // 在 significance_level 下是否拒绝 H0
};

}  // namespace detail
}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
