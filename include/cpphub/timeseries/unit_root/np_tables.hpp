// =============================================================================
// np_tables.hpp - Ng-Perron 2001 Table 1 渐近临界值 (spec §2.1 / readiness D1)
//
// Phase 7C v1.7 M0 (PHASE7C_SPEC.md v1.2 §1.1.3)
//
// 教材锚点: Ng-Perron 2001 Econometrica 69(6):1519-1554, Table I
//
// 溯源 (2026-08-17 双源取证, 24/24 值零差异):
//   [S1] AU OA 副本 (BC wp369 工作稿, PDF p.30 Table 1):
//        https://w.american.edu/cas/economics/gaussres/authors/perron_pierre/ng_perron00.pdf
//        逐字: "Case: p=0,c=-7.0 / .01 -13.8 -2.58 .174 1.78 / .05 -8.1 -1.98 .233 3.17 /
//               .10 -5.7 -1.62 .275 4.45 / Case: p=1,c=-13.5 /
//               .01 -23.8 -3.42 .143 4.03 / .05 -17.3 -2.91 .168 5.48 / .10 -14.2 -2.62 .185 6.67"
//   [S2] Econometrica 发表版 (作者存档, p.1524 TABLE I):
//        https://www.columbia.edu/~sn2294/pub/ecta01.pdf — 同上逐字一致
//   [S3] MetricGate (第三源, 仅 5% 八值) — 一致
//   脚注: p=0 的 MZα/MZt/ADFGLS 三列取自 Fuller (1976); 其余条目
//         = 20,000 次模拟 × 5,000 步 Wiener 逼近
//
// 排幻觉点 (spec §9.5):
//   NP5: 拒绝方向 — 四统计量统一 "statistic < cv 拒绝 H0(单位根)":
//        MZα/MZt 越负越拒绝, MSB/MPT 越小越拒绝, 代数上同为 stat < cv
//   ⚠️ 转录陷阱 (双源取证拦截, 2026-08-17): 趋势情形 MSB 1% = 0.143 (非 0.121),
//        MPT 1% = 4.03 (非 4.47) — 0.121/4.47 变体系无源讹传, 两版原文均无此值
//
// 关联: ng_perron_test.hpp (消费方); static_assert 锚 = spec §2.1 冻结 5% 值
// =============================================================================
#pragma once

#include <array>
#include <stdexcept>
#include <string>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace unit_root {

/// NP 四统计量 (Ng-Perron 2001 §3.1)
enum class NPStat {
    MZa,  ///< MZα = (T⁻¹ỹ_T² − s²_AR) / (2·T⁻²Σỹ²_{t−1})
    MZt,  ///< MZt = MZα × MSB (恒等式, 实现侧自检)
    MSB,  ///< MSB = √(T⁻²Σỹ²_{t−1} / s²_AR)
    MPT   ///< MPT 分情形 (NP4): p=0 末项 −c̄; p=1 末项 (1−c̄)
};

namespace detail {

/// Table 1 布局: [stat][level], level 顺序 {1%, 5%, 10%}
/// 常数情形 (p=0, c̄=−7.0) — 溯源 [S1]/[S2] "Case: p=0,c=-7.0" 三行逐字
constexpr std::array<std::array<Real, 3>, 4> NP_CV_CONSTANT{{
    {{-13.8, -8.1, -5.7}},   // MZα  (1%/5%/10%; 5% 取自 Fuller 1976)
    {{-2.58, -1.98, -1.62}}, // MZt  (同上)
    {{0.174, 0.233, 0.275}}, // MSB  (20,000 模拟)
    {{1.78, 3.17, 4.45}}     // MPT
}};

/// 趋势情形 (p=1, c̄=−13.5) — 溯源 [S1]/[S2] "Case: p=1,c=-13.5" 三行逐字
/// ⚠️ MSB 1%=0.143 / MPT 1%=4.03 为原文值 (0.121/4.47 讹传变体, 见文件头)
constexpr std::array<std::array<Real, 3>, 4> NP_CV_TREND{{
    {{-23.8, -17.3, -14.2}}, // MZα
    {{-3.42, -2.91, -2.62}}, // MZt
    {{0.143, 0.168, 0.185}}, // MSB
    {{4.03, 5.48, 6.67}}     // MPT
}};

// ---------------------------------------------------------------------------
// 编译期守护: 5% 八值 = spec §2.1 冻结锚 (PHASE7C_SPEC v1.2, R 门禁审计链)
// ---------------------------------------------------------------------------
static_assert(NP_CV_CONSTANT[0][1] == -8.10, "NP5 anchor: MZa const 5% = -8.10");
static_assert(NP_CV_CONSTANT[1][1] == -1.98, "NP5 anchor: MZt const 5% = -1.98");
static_assert(NP_CV_CONSTANT[2][1] == 0.233, "NP5 anchor: MSB const 5% = 0.233");
static_assert(NP_CV_CONSTANT[3][1] == 3.17, "NP5 anchor: MPT const 5% = 3.17");
static_assert(NP_CV_TREND[0][1] == -17.30, "NP5 anchor: MZa trend 5% = -17.30");
static_assert(NP_CV_TREND[1][1] == -2.91, "NP5 anchor: MZt trend 5% = -2.91");
static_assert(NP_CV_TREND[2][1] == 0.168, "NP5 anchor: MSB trend 5% = 0.168");
static_assert(NP_CV_TREND[3][1] == 5.48, "NP5 anchor: MPT trend 5% = 5.48");
// 防呆: 拦截常见讹传变体混入 (双源原文均无 0.121/4.47)
static_assert(NP_CV_TREND[2][0] == 0.143, "NP5 trap: MSB trend 1% = 0.143 (非 0.121)");
static_assert(NP_CV_TREND[3][0] == 4.03, "NP5 trap: MPT trend 1% = 4.03 (非 4.47)");
// 单调性 (临界值表内在结构: MZα/MZt 随水平升而增, MSB/MPT 亦增)
static_assert(NP_CV_CONSTANT[0][0] < NP_CV_CONSTANT[0][1] &&
              NP_CV_CONSTANT[0][1] < NP_CV_CONSTANT[0][2], "MZa const monotone");
static_assert(NP_CV_TREND[3][0] < NP_CV_TREND[3][1] &&
              NP_CV_TREND[3][1] < NP_CV_TREND[3][2], "MPT trend monotone");

}  // namespace detail

/// NP 渐近临界值查表 (Ng-Perron 2001 Table I, 双源转录)
/// @param stat 四统计量之一
/// @param trend_spec "c" (常数, p=0, c̄=−7.0) / "ct" (趋势, p=1, c̄=−13.5)
/// @param p 显著性水平 (0.01 / 0.05 / 0.10; 表仅此三档)
/// @return 临界值 (渐近; 精确到原文 3 位有效数字)
/// @throws std::invalid_argument 若 trend_spec/p 不在支持档位
/// NP5 拒绝方向: 四统计量统一 statistic < cv 即拒绝 H0 (单位根)
inline Real np_critical_value(NPStat stat, const std::string& trend_spec, Real p) {
    const auto& table = (trend_spec == "c") ? detail::NP_CV_CONSTANT
                        : (trend_spec == "ct") ? detail::NP_CV_TREND
                                               : throw std::invalid_argument(
                                                     "np_critical_value: trend_spec must be c/ct");
    Size lvl = 0;
    if (p == 0.10) {
        lvl = 2;
    } else if (p == 0.05) {
        lvl = 1;
    } else if (p == 0.01) {
        lvl = 0;
    } else {
        throw std::invalid_argument(
            "np_critical_value: p must be 0.01/0.05/0.10 (Table 1 仅三档)");
    }
    return table[static_cast<Size>(stat)][lvl];
}

}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
