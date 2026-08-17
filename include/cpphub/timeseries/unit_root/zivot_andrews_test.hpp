// =============================================================================
// zivot_andrews_test.hpp - Zivot-Andrews 结构断点单位根检验 (spec §2.2)
//
// Phase 7C v1.7 M0 (PHASE7C_SPEC.md v1.2)
//
// 教材锚点: Zivot-Andrews 1992 JBES 10(3):251-270 (Tables 2/3/4)
//   ⚠️ statsmodels 文档误写刊名 "J. Business & Economic Studies" (实为 Statistics),
//      勿照抄 (调研 ZA 勘误 + 2026-08-17 复核仍误)
// 对照库: statsmodels 0.14.6 zivot_andrews (Baum 近似, 主对照 1e-10)
//       + R urca 1.3-4 ur.za (固定 lag, 1e-8)
//
// 幻觉点防护 (spec §9.5):
//   ZA1: 三库 lag 策略差异 → 双模式 API (固定 lag 主 / Baum 预选对照)
//   ZA2: trim 参数化 (默认 0.15, 上限 1/3; urca 无 trim → 复现须小 trim 放开网格)
//   ZA3: DU_t = 1(t > Tb); DT_t = (t − Tb)·1(t > Tb) 断点后重新计时 (非 DU·t 全局)
//   ZA4: 双临界值表 — 主 = ZA1992 论文表 (za1992_cv.inc 同源), MC 表 (za_mc_cv.inc
//        同源, statsmodels/arch 共用) 仅作 p 值插值与 statsmodels 对照;
//        ⚠️ MC c 1% = −5.27644 (非论文 −5.34); −5.83 是 MC 0.1% 分位值
//   ZA5: 统计量 = min_{Tb} t(α̂) (最负), 非 max
//   spec v1.2 勘误: Model A 恒含趋势项 (三源实证: statsmodels L2700 exog[:,2]=trend
//        于 regression="c" 亦在 / urca datmat trend 列 / Baum .ado baseline 含 trend)
//
// statsmodels 对齐语义 (本地 stattools.py L2590-2738 一手, 2026-08-17):
//   - Baum 模式 lag = adfuller(x, regression="ct", autolag="AIC") 语义
//     → 复用 select_lag_by_ic(data, "ct", 0, "aic") (7B 已实证同口径)
//   - 网格: trimcnt = int(T·trim), 候选 b ∈ [trimcnt+1, T−trimcnt]
//   - DU(0-based t) = 1(t ≥ b) ⟺ 1-based 1(t > b) — 与 urca z 约定一致
//   - DT 断点后自 1 重新计时 (urca 形); statsmodels "ct" 形自 2 起 — 差为 DU 的
//     常数倍, DU 在 Model C span 内 ⇒ C 模型 t-stat 不受影响 (代数等价);
//     Model B ("t") statsmodels 另有 cutoff−1 边界 quirk, 与 urca 形不同,
//     verify_za.py 量化并以 urca 为主对照
//   - 报告约定: statsmodels bpidx = argmin − 1 ⇒ 本实现 break_index = b* − 1 (0-based)
//   - p 值: np.interp(stat, cv 升序, pct)/100, 端点 clamp (L2529 复刻)
//
// 复用 (readiness B2, additive-only §8.5): detail::ols_fit (bse/SSR),
//   schwert_lag, select_lag_by_ic (unit_root_common.hpp, 不修改)
// =============================================================================
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace unit_root {

/// ZA 三模型 (ZA1992 §2; urca: intercept/trend/both; statsmodels: c/t/ct)
enum class ZAModel {
    A,  ///< 崩溃均值: c + trend + DU + y_{t−1} + lags (urca "intercept", sm "c")
    B,  ///< 断裂趋势: c + trend + DT + y_{t−1} + lags (urca "trend", sm "t")
    C   ///< 两者:     c + trend + DU + DT + y_{t−1} + lags (urca "both", sm "ct")
};

/// ZA 检验结果 (spec §2.2 签名 + §1.4 默认初始化)
struct ZAResult {
    Real statistic = 0.0;         ///< min_{Tb} t(α̂) (ZA5: 最负)
    Size break_index = 0;         ///< 0-based y 索引 (statsmodels bpidx 约定: b* − 1)
    Real critical_1pct = 0.0;     ///< ZA1992 论文表 (主, za1992_cv.inc 同源)
    Real critical_5pct = 0.0;
    Real critical_10pct = 0.0;
    Real p_value_baum = 0.0;      ///< MC 表插值 p (statsmodels 对照; 论文表无 p)
    ZAModel model = ZAModel::C;
    Size n_lags = 0;
    Real trim = 0.15;
    std::vector<Real> t_stats_path;  ///< 逐断点 t(α̂) 轨迹 (诊断)
    bool reject_null = false;     ///< H0: 带单断点单位根; stat < cv_5pct 拒绝
    std::string summary;
};

namespace detail {

/// ZA1992 论文临界值 runtime 副本 — 与 tests/.../za1992_cv.inc 双录互校
/// (溯源: urca ur.za 源码 + Baum zandrews.ado Table 2/3/4 + 讲义, 三源 9/9)
constexpr std::array<std::array<Real, 3>, 3> ZA1992_CV{{
    {{-5.34, -4.80, -4.58}},  // Model A (Table 2)
    {{-4.93, -4.42, -4.11}},  // Model B (Table 3)
    {{-5.57, -5.08, -4.82}}   // Model C (Table 4)
}};
static_assert(ZA1992_CV[0][0] == -5.34 && ZA1992_CV[0][1] == -4.80 &&
                  ZA1992_CV[0][2] == -4.58, "ZA1992 A row");
static_assert(ZA1992_CV[1][0] == -4.93 && ZA1992_CV[2][0] == -5.57, "ZA1992 B/C 1%");
static_assert(ZA1992_CV[2][2] == -4.82, "ZA1992 C 10%");

/// MC 百分位点 (pct: 0.001~99.9; cv: 统计量) — 与 za_mc_cv.inc 双录互校
/// (程序化生成: statsmodels 0.14.6 与 arch 8.0.0 本地双库 144/144 零差异,
///  gen_za_mc_tables.py; MC = 100,000 模拟 × 2,000 点)
struct ZaMcPoint {
    Real pct;
    Real cv;
};

constexpr std::array<ZaMcPoint, 48> ZA_MC_C{{
    {0.001, -6.78442},
    {0.100, -5.83192},
    {0.200, -5.68139},
    {0.300, -5.58461},
    {0.400, -5.51308},
    {0.500, -5.45043},
    {0.600, -5.39924},
    {0.700, -5.36023},
    {0.800, -5.33219},
    {0.900, -5.30294},
    {1.000, -5.27644},
    {2.500, -5.0334},
    {5.000, -4.81067},
    {7.500, -4.67636},
    {10.000, -4.56618},
    {12.500, -4.4813},
    {15.000, -4.40507},
    {17.500, -4.33947},
    {20.000, -4.28155},
    {22.500, -4.22683},
    {25.000, -4.1783},
    {27.500, -4.13101},
    {30.000, -4.08586},
    {32.500, -4.04455},
    {35.000, -4.0038},
    {37.500, -3.96144},
    {40.000, -3.92078},
    {42.500, -3.88178},
    {45.000, -3.84503},
    {47.500, -3.80549},
    {50.000, -3.77031},
    {52.500, -3.73209},
    {55.000, -3.696},
    {57.500, -3.65985},
    {60.000, -3.62126},
    {65.000, -3.5458},
    {70.000, -3.46848},
    {75.000, -3.38533},
    {80.000, -3.29112},
    {85.000, -3.17832},
    {90.000, -3.04165},
    {92.500, -2.95146},
    {95.000, -2.83179},
    {96.000, -2.76465},
    {97.000, -2.68624},
    {98.000, -2.57884},
    {99.000, -2.40044},
    {99.900, -1.88932},
}};

constexpr std::array<ZaMcPoint, 48> ZA_MC_T{{
    {0.001, -83.9094},
    {0.100, -13.8837},
    {0.200, -9.13205},
    {0.300, -6.32564},
    {0.400, -5.60803},
    {0.500, -5.38794},
    {0.600, -5.26585},
    {0.700, -5.18734},
    {0.800, -5.12756},
    {0.900, -5.07984},
    {1.000, -5.03421},
    {2.500, -4.65634},
    {5.000, -4.4058},
    {7.500, -4.25214},
    {10.000, -4.13678},
    {12.500, -4.03765},
    {15.000, -3.95185},
    {17.500, -3.87945},
    {20.000, -3.81295},
    {22.500, -3.75273},
    {25.000, -3.69836},
    {27.500, -3.64785},
    {30.000, -3.59819},
    {32.500, -3.55146},
    {35.000, -3.50522},
    {37.500, -3.45987},
    {40.000, -3.41672},
    {42.500, -3.37465},
    {45.000, -3.33394},
    {47.500, -3.29393},
    {50.000, -3.25316},
    {52.500, -3.21244},
    {55.000, -3.17124},
    {57.500, -3.13211},
    {60.000, -3.09204},
    {65.000, -3.01135},
    {70.000, -2.92897},
    {75.000, -2.83614},
    {80.000, -2.73893},
    {85.000, -2.6284},
    {90.000, -2.49611},
    {92.500, -2.41337},
    {95.000, -2.3082},
    {96.000, -2.25797},
    {97.000, -2.19648},
    {98.000, -2.1132},
    {99.000, -1.99138},
    {99.900, -1.67466},
}};

constexpr std::array<ZaMcPoint, 48> ZA_MC_CT{{
    {0.001, -38.178},
    {0.100, -6.43107},
    {0.200, -6.07279},
    {0.300, -5.95496},
    {0.400, -5.86254},
    {0.500, -5.77081},
    {0.600, -5.72541},
    {0.700, -5.68406},
    {0.800, -5.65163},
    {0.900, -5.60419},
    {1.000, -5.57556},
    {2.500, -5.29704},
    {5.000, -5.07332},
    {7.500, -4.93003},
    {10.000, -4.82668},
    {12.500, -4.73711},
    {15.000, -4.6602},
    {17.500, -4.5897},
    {20.000, -4.52855},
    {22.500, -4.471},
    {25.000, -4.42011},
    {27.500, -4.37387},
    {30.000, -4.32705},
    {32.500, -4.28126},
    {35.000, -4.23793},
    {37.500, -4.19822},
    {40.000, -4.158},
    {42.500, -4.11946},
    {45.000, -4.08064},
    {47.500, -4.04286},
    {50.000, -4.00489},
    {52.500, -3.96837},
    {55.000, -3.932},
    {57.500, -3.89496},
    {60.000, -3.85577},
    {65.000, -3.77795},
    {70.000, -3.69794},
    {75.000, -3.61852},
    {80.000, -3.52485},
    {85.000, -3.41665},
    {90.000, -3.28527},
    {92.500, -3.19724},
    {95.000, -3.08769},
    {96.000, -3.03088},
    {97.000, -2.96091},
    {98.000, -2.85581},
    {99.000, -2.71015},
    {99.900, -2.28767},
}};

// runtime 侧锚 (全量互校在测试层: GOLDEN == runtime 逐值)
static_assert(ZA_MC_C.size() == 48 && ZA_MC_T.size() == 48 && ZA_MC_CT.size() == 48,
              "MC tables 48 points");
static_assert(ZA_MC_C[10].pct == 1.0 && ZA_MC_C[10].cv == -5.27644, "MC c 1%");
static_assert(ZA_MC_C[12].cv == -4.81067 && ZA_MC_C[14].cv == -4.56618, "MC c 5/10%");
static_assert(ZA_MC_T[10].cv == -5.03421 && ZA_MC_CT[10].cv == -5.57556, "MC t/ct 1%");
static_assert(ZA_MC_C[1].pct == 0.1 && ZA_MC_C[1].cv == -5.83192,
              "ZA4 trap: -5.83 is 0.1-percentile, never 1%");

/// MC p 值插值 (statsmodels L2529 复刻: np.interp(stat, cv 升序, pct)/100, 端点 clamp)
inline Real za_mc_pvalue(Real stat, ZAModel m) {
    const std::array<ZaMcPoint, 48>& tbl =
        (m == ZAModel::A) ? ZA_MC_C : (m == ZAModel::B) ? ZA_MC_T : ZA_MC_CT;
    if (stat <= tbl.front().cv) {
        return tbl.front().pct / 100.0;
    }
    if (stat >= tbl.back().cv) {
        return tbl.back().pct / 100.0;
    }
    for (Size i = 1; i < tbl.size(); ++i) {
        if (stat < tbl[i].cv) {
            const ZaMcPoint& a = tbl[i - 1];
            const ZaMcPoint& b = tbl[i];
            const Real w = (stat - a.cv) / (b.cv - a.cv);
            return (a.pct + w * (b.pct - a.pct)) / 100.0;
        }
    }
    return tbl.back().pct / 100.0;  // 不可达
}

}  // namespace detail

/// Zivot-Andrews 结构断点单位根检验 (H0: 带单断点单位根)
///
/// @param data 原始序列 (水平值, 无 NaN)
/// @param model A/B/C (崩溃均值/断裂趋势/两者; 三模型均含趋势项 — spec v1.2 勘误)
/// @param fixed_lag 固定滞后数: 0 => Schwert 自动; ≥1 => 用户固定 (urca 对齐)
/// @param baum_preselect true => Baum 式一次性预选 (statsmodels 对齐:
///        adfuller(regression="ct", autolag="AIC") 语义, 全样本选 lag 后冻结)
/// @param trim 网格修剪比例 (0 < trim ≤ 1/3; 默认 0.15; urca 复现须放开)
/// @return ZAResult (论文表 cv 主 + MC 插值 p 对照)
/// @throws std::invalid_argument 样本过小/含 NaN/零方差/trim 越界/网格空/自由度不足
inline ZAResult zivot_andrews_test(const std::vector<Real>& data,
                                   ZAModel model = ZAModel::C,
                                   Size fixed_lag = 0,
                                   bool baum_preselect = false,
                                   Real trim = 0.15) {
    const Size T = data.size();
    if (T < 10) {
        throw std::invalid_argument("zivot_andrews_test: sample too small (T < 10)");
    }
    if (!(trim > 0.0 && trim <= 1.0 / 3.0)) {
        throw std::invalid_argument("zivot_andrews_test: trim must be in (0, 1/3]");
    }
    Real ymin = data[0], ymax = data[0];
    for (Real v : data) {
        if (std::isnan(v)) {
            throw std::invalid_argument("zivot_andrews_test: NaN in data");
        }
        if (v < ymin) ymin = v;
        if (v > ymax) ymax = v;
    }
    if (ymax == ymin) {
        throw std::invalid_argument("zivot_andrews_test: zero variance (constant)");
    }

    // ---- Step 1: 滞后选择 (ZA1 双模式) ----
    Size k = 0;
    if (baum_preselect) {
        // statsmodels: adfuller(x, maxlag=None, regression="ct", autolag="AIC")
        k = select_lag_by_ic(data, "ct", 0, "aic");
    } else {
        k = (fixed_lag > 0) ? fixed_lag : schwert_lag(T);  // 主模式: 固定 lag
    }

    // ---- Step 2: 断点网格 (statsmodels 边界: trimcnt = int(T·trim), b ∈ [trimcnt+1, T−trimcnt]) ----
    // 秩保护 (等价 statsmodels L2718 rank-check): DU 需行集内变异 ⇒ b ≥ k+2 (行 t 最小 = k+1);
    // DT 需至少一个非零行 ⇒ b ≤ T−1 (行 t 最大 = T−1)。极小 trim 时边界内收, doc 注明。
    const Size trimcnt = static_cast<Size>(T * trim);
    const Size b_lo = std::max<Size>(trimcnt + 1, k + 2);
    const Size b_hi = std::min<Size>(T - trimcnt, T - 1);
    if (b_hi <= b_lo) {
        throw std::invalid_argument("zivot_andrews_test: empty break grid (trim too large)");
    }

    // ---- Step 3: 回归数据 (行 r ↔ 0-based t = r + k + 1; Δy_t = y[t] − y[t−1]) ----
    const Size R = T - k - 1;                  // 有效行数
    // 列: 模型公共 [1, DU?, trend, DT?, ylag, Δy 滞后 ×k]
    const bool has_du = (model != ZAModel::B);
    const bool has_dt = (model != ZAModel::A);
    const Size ncol = 1 + (has_du ? 1 : 0) + 1 + (has_dt ? 1 : 0) + 1 + k;
    const Size gamma_col = 1 + (has_du ? 1 : 0) + 1 + (has_dt ? 1 : 0);  // ylag 位置
    if (R <= ncol) {
        throw std::invalid_argument(
            "zivot_andrews_test: insufficient degrees of freedom (reduce lags/trim)");
    }

    // 预备: Δy 与其滞后 (行对齐)
    std::vector<Real> dy(R);        // Δy_t, 行 r ↔ t = r + k + 1
    std::vector<Real> dy_lags(R * k, 0.0);  // 行 r 的第 j 滞后 (j=1..k): Δy_{t−j}
    for (Size r = 0; r < R; ++r) {
        const Size t = r + k + 1;   // 0-based y 索引
        dy[r] = data[t] - data[t - 1];
        for (Size j = 1; j <= k; ++j) {
            dy_lags[r * k + (j - 1)] = data[t - j] - data[t - j - 1];
        }
    }

    // ---- Step 4: 断点搜索: 统计量 = min t(α̂) (ZA5), 首个 argmin ----
    ZAResult res;
    res.model = model;
    res.n_lags = k;
    res.trim = trim;
    Real best_stat = std::numeric_limits<Real>::infinity();
    Size best_b = b_lo;
    res.t_stats_path.reserve(b_hi - b_lo + 1);

    std::vector<std::vector<Real>> X(R, std::vector<Real>(ncol, 0.0));
    for (Size b = b_lo; b <= b_hi; ++b) {
        for (Size r = 0; r < R; ++r) {
            const Size t = r + k + 1;             // 0-based
            Size c = 0;
            X[r][c++] = 1.0;                      // 常数
            if (has_du) {
                // ZA3: DU = 1(t ≥ b) 0-based ⟺ 1-based 1(t > Tb), Tb = b
                X[r][c++] = (t >= b) ? 1.0 : 0.0;
            }
            // trend/DT 按 T 缩放 (列乘常数不改 t-stat; 对齐 statsmodels 数值稳定化,
            // 正规方程条件数 ↓, urca QR 精度可达)
            X[r][c++] = static_cast<Real>(t) / static_cast<Real>(T);
            if (has_dt) {
                // ZA3: DT = (t − b + 1)·1(t ≥ b) — 断点后自 1 重新计时 (urca 形)
                X[r][c++] =
                    (t >= b) ? static_cast<Real>(t - b + 1) / static_cast<Real>(T) : 0.0;
            }
            X[r][c++] = data[t - 1];              // y_{t−1} (γ 所在列)
            for (Size j = 0; j < k; ++j) {
                X[r][c++] = dy_lags[r * k + j];
            }
        }
        // 奇异/近奇异候选点 (如极小 b 下 DU/DT 与 trend 近共线): 跳过, 不进 min
        // (语义同 urca lm 产生 NA 不进 which.min / statsmodels quick_ols 不检查)
        Real tstat = std::numeric_limits<Real>::infinity();
        try {
            const detail::OlsFull ols = detail::ols_fit(dy, X);
            tstat = ols.beta[gamma_col] / ols.bse[gamma_col];
        } catch (const std::invalid_argument&) {
            ;  // +inf 保留
        }
        res.t_stats_path.push_back(tstat);
        if (tstat < best_stat) {                  // 严格 < ⇒ 首个 argmin (np.argmin/which.min)
            best_stat = tstat;
            best_b = b;
        }
    }

    // ---- Step 5: 结果 ----
    const auto& cvrow = detail::ZA1992_CV[static_cast<Size>(model)];
    res.statistic = best_stat;
    res.break_index = best_b - 1;  // statsmodels bpidx 约定 (0-based; urca z = best_b)
    res.critical_1pct = cvrow[0];
    res.critical_5pct = cvrow[1];
    res.critical_10pct = cvrow[2];
    res.p_value_baum = detail::za_mc_pvalue(best_stat, model);
    res.reject_null = best_stat < res.critical_5pct;
    res.summary = "Zivot-Andrews (model " + std::to_string(static_cast<int>(model) + 1) +
                  ", lag=" + std::to_string(k) + ", trim=" + std::to_string(trim) +
                  "): H0 = unit root with single structural break; "
                  "min t(alpha) = " + std::to_string(best_stat) +
                  " at break idx " + std::to_string(res.break_index);
    return res;
}

}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
