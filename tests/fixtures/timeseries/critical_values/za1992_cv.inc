// =============================================================================
// za1992_cv.inc - Zivot-Andrews 1992 论文临界值黄金值 (spec §1.2.21 / readiness D2)
//
// Phase 7C v1.7 M0 (PHASE7C_SPEC.md v1.2; 消费方: zivot_andrews_test.hpp 与
// test_zivot_andrews — 测试将本黄金副本与头文件运行时表逐值断言相等, 双录互校)
//
// 溯源 (2026-08-17 三源一致, 9/9 值零差异):
//   [S1] urca 1.3-4 (已装 F:/R/win-library/4.6) ur.za 源码 deparse 逐字:
//        L39  intercept: cval <- c(-5.34, -4.8, -4.58)
//        L53  trend:     cval <- c(-4.93, -4.42, -4.11)
//        L69  both:      cval <- c(-5.57, -5.08, -4.82)
//        (各为 {1%, 5%, 10%}; 模型参数名 intercept/trend/both, lag 为数值/NULL 无 trim)
//   [S2] Baum zandrews.ado v1.0.5 (fmwww.bc.edu/repec/bocode/z/zandrews.ado,
//        statsmodels Baum 近似的正主) — 注释标明论文表号:
//        "Table 2, Model A" -5.34/-4.80/-4.58; "Table 3, Model B" -4.93/-4.42/-4.11;
//        "Table 4, Model C" -5.57/-5.08/-4.82
//   [S3] 讲义转录 (randall-romero.github.io/econometria, "Fuente: Zivot1992"):
//        A {-5.34,-4.80,-4.58} / B {-4.93,-4.42,-4.11} / C {-5.57,-5.08,-4.82}
//   论文锚: Zivot-Andrews 1992 JBES 10(3):251-270, Tables 2/3/4 (c 1% = −5.34,
//   Phase 7C 调研 ZA4 已裁决: 论文 −5.34 vs MC −5.27644)
//
// 排幻觉点 (spec §9.5 ZA4):
//   ⚠️ 本表 = ZA1992 论文表 (urca/Baum 采用); 与 statsmodels/arch 共用的 MC 表
//      (c 1% = −5.27644, 100,000 次模拟) 是两套表, 数值不得混用;
//      MC 表的 −5.83 为 0.1% 分位值, 非 1% (v1.1 审计修正)。MC 表另录 za_mc_cv.inc
//
// 用法: 在测试 TU 全局作用域 #include 本文件 (无外部依赖, 内含 <array>;
//       黄金值用裸 double, 与 Real=double §8.6 一致, 保持零依赖独立性)
// =============================================================================
#pragma once

#include <array>

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace unit_root {
namespace detail {

/// ZA1992 论文临界值黄金副本 — 行序 {Model A, Model B, Model C},
/// 列序 {1%, 5%, 10%}; Model A=崩溃均值(intercept/DU), B=断裂趋势(trend/DT),
/// C=两者(both) — 与 ZAModel enum (zivot_andrews_test.hpp) 一一对应
constexpr std::array<std::array<double, 3>, 3> ZA1992_CV_GOLDEN{{
    {{-5.34, -4.80, -4.58}},  // Model A (ZA1992 Table 2; urca "intercept")
    {{-4.93, -4.42, -4.11}},  // Model B (ZA1992 Table 3; urca "trend")
    {{-5.57, -5.08, -4.82}}   // Model C (ZA1992 Table 4; urca "both")
}};

// ---------------------------------------------------------------------------
// 编译期守护: 三源转录身份锚 (9 值逐一) + 结构单调性 + MC 混表防呆
// ---------------------------------------------------------------------------
static_assert(ZA1992_CV_GOLDEN[0][0] == -5.34, "ZA1992 A 1% (S1 L39/S2 Table2/S3)");
static_assert(ZA1992_CV_GOLDEN[0][1] == -4.80, "ZA1992 A 5%");
static_assert(ZA1992_CV_GOLDEN[0][2] == -4.58, "ZA1992 A 10%");
static_assert(ZA1992_CV_GOLDEN[1][0] == -4.93, "ZA1992 B 1% (S1 L53/S2 Table3/S3)");
static_assert(ZA1992_CV_GOLDEN[1][1] == -4.42, "ZA1992 B 5%");
static_assert(ZA1992_CV_GOLDEN[1][2] == -4.11, "ZA1992 B 10%");
static_assert(ZA1992_CV_GOLDEN[2][0] == -5.57, "ZA1992 C 1% (S1 L69/S2 Table4/S3)");
static_assert(ZA1992_CV_GOLDEN[2][1] == -5.08, "ZA1992 C 5%");
static_assert(ZA1992_CV_GOLDEN[2][2] == -4.82, "ZA1992 C 10%");
// 单调性: 1% < 5% < 10% (越显著越负)
static_assert(ZA1992_CV_GOLDEN[0][0] < ZA1992_CV_GOLDEN[0][1] &&
              ZA1992_CV_GOLDEN[0][1] < ZA1992_CV_GOLDEN[0][2], "ZA1992 A monotone");
static_assert(ZA1992_CV_GOLDEN[1][0] < ZA1992_CV_GOLDEN[1][1] &&
              ZA1992_CV_GOLDEN[1][1] < ZA1992_CV_GOLDEN[1][2], "ZA1992 B monotone");
static_assert(ZA1992_CV_GOLDEN[2][0] < ZA1992_CV_GOLDEN[2][1] &&
              ZA1992_CV_GOLDEN[2][1] < ZA1992_CV_GOLDEN[2][2], "ZA1992 C monotone");
// 防混表: 论文表 A 1% = −5.34, 不得混入 MC 表值 −5.27644 (ZA4 两套表)
static_assert(ZA1992_CV_GOLDEN[0][0] != -5.27644,
              "ZA4 trap: paper table A 1% = -5.34, MC -5.27644 belongs za_mc_cv.inc");
static_assert(ZA1992_CV_GOLDEN[0][0] != -5.83,
              "ZA4 trap: -5.83 is MC 0.1% percentile, never a 1% value");

}  // namespace detail
}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
