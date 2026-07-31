#pragma once
// SOURCE: Basel Committee on Banking Supervision (2014) "The standardised approach for
//         measuring counterparty credit risk exposures" (BCBS d279)
// SOURCE: BCBS (2019) SA-CCR consolidating document
// SOURCE: Zhu, Pykhtin (2007) "A Guide to Modeling Counterparty Credit Risk"
// SOURCE: Canabarro, Duffie (2003) "Measuring and Marking Counterparty Risk"
// 模块: PFE (Potential Future Exposure) + SA-CCR (Standardised Approach for Counterparty Credit Risk)
//
// ==================== PFE 数学框架 ====================
//
// 设 V(t) 为衍生品组合在 t 时刻的市值 (从银行视角), 沿 MC 路径获得样本
// V_samples[path][time]. 暴露:
//   E(t) = max(V(t), 0)    (正值暴露, bank 收钱)
//   NE(t) = max(-V(t), 0)  (负值暴露, bank 付钱)
//
// 1. Expected Exposure (EE):
//      EE(t) = E[ max(V(t), 0) ]
//
// 2. Expected Positive Exposure (EPE):
//      EPE = (1/T) ∫_0^T EE(t) dt  (时间平均)
//
// 3. Effective Expected Exposure (EEE): EE 的 running maximum
//      EEE(t) = max_{s≤t} EE(s)  (体现暴露单调累积的监管保守假设)
//
// 4. Effective Expected Positive Exposure (EEPE): Basel 用作 SA-CCR 输入
//      EEPE = (1/T) ∫_0^T EEE(t) dt
//
// 5. Potential Future Exposure (PFE): 暴露分布的 α 分位数
//      PFE_α(t) = inf{x : P(E(t) ≤ x) ≥ α}
//    Basel α = 99% (1 年期, 99% PFE)
//
// ==================== SA-CCR 数学框架 (BCBS d279) ====================
//
// Replacement Cost (RC) — 当前替换成本:
//   RC = max(V - C, 0)             (V = 组合市值, C = 抵押品净额)
//   若 V - C ≤ 0, RC = 0 (无当前暴露)
//
// Potential Future Exposure (PFE_addon) — addon 反映未来潜在暴露:
//   PFE_addon = mul_addon * aggregate_addon
//   aggregate_addon = Σ_h (SF_h × EF_h × DF_h × MF_h)  (各 hedging set 汇总)
//   mul_addon = 1.4 (supervisory multiplier, 反映 maturity 与 correlation)
//
// 简化版 (单 hedging set, 5 个资产类别):
//   PFE_addon = Σ_i (notional_i × RF_i × MF_i)   (RF = supervisory factor)
//
// Supervisory Factors (RF) per asset class:
//   Interest Rate:       0.5%   (0.005)
//   FX:                  4.0%   (0.04)
//   Credit (QQ):        3.6%~6.0% (取决 rating)
//   Equity:             12%~32%  (取决 index/single)
//   Commodity:          18%~40%
//
// Maturity Factor (MF) — 期限因子 (剩余期限 < 1Y 用 under-1Y 公式):
//   MF_i = min{1, (M_i - t_now)/1Y}   若剩余期限 ≥ 1Y
//   MF_i = sqrt((M_i - t_now)/1Y)     若剩余期限 < 1Y
//   (Basel 2019 修订: < 1Y 用 √(T/1Y), 体现短期高频展期)
//
// Hedging Set Level Aggregation:
//   SF_h = supervisory factor × (1 + 0.4 × net_to_gross_ratio_h)
//   其中 net_to_gross_ratio = Σ|net notional| / Σ|gross notional|
//
// SA-CCR Exposure at Default (EAD):
//   EAD = 1.4 × (RC + PFE_addon)
//   1.4 为 supervisory factor (Basel 默认, 监管保守系数)
//
// ==================== 与 PFE (MC) 的关系 ====================
//
// SA-CCR 是监管公式 (保守, 解析, 资本计提), 与 MC PFE 在概念上对应:
//   - SA-CCR PFE_addon ≈ 模型 PFE_99% × mul_factor
//   - Basel 允许 IMM (Internal Model Method) 银行用 MC PFE 替代 SA-CCR addon
//   - 但 SA-CCR 为标准方法, 不依赖 MC, 用于非 IMM 银行或 fallback

#include "cpphub/core/types.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include "cpphub/instruments/credit/credit_curve.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/core/rng.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <functional>

namespace cpphub {
inline namespace v1 {

// ============ Exposure Stats (沿 MC 路径的暴露统计) ============
// 从 V_samples[path][time] 计算各时刻的暴露分布统计量
struct ExposureStats {
    std::vector<Real> times;        // 时间网格, n+1 点
    std::vector<Real> ee;           // Expected Exposure: E[max(V(t),0)]
    std::vector<Real> ene;          // Expected Negative Exposure: E[max(-V(t),0)]
    std::vector<Real> eee;          // Effective EE: running max of EE
    std::vector<Real> pfe;          // PFE (α 分位数)
    Real epe = 0.0;                 // Expected Positive Exposure (时间平均 EE)
    Real eepe = 0.0;                // Effective EPE (时间平均 EEE)
    Real max_pfe = 0.0;             // 最大 PFE (沿时间)
    Real max_ee = 0.0;              // 最大 EE
};

// ============ PFE 配置 ============
struct PFEConfig {
    Real confidence_level = 0.99;   // Basel 默认 99%
    bool use_running_max_ee = true; // true: EEE = running max(EE); false: EEE = EE
};

// ============ 从 MC 样本计算暴露统计 ============
// V_samples[path_idx][time_idx] = 衍品组合在第 path 条路径 time_idx 处的市值 V(t)
// (折现或未折现均可, 内部按 discount_curve 折现)
// times: 时间网格, 长度 = V_samples 列数
inline ExposureStats compute_exposure_stats(
        const std::vector<std::vector<Real>>& V_samples,
        const std::vector<Real>& times,
        const ZeroCurve& discount_curve,
        const PFEConfig& cfg = {}) {
    if (V_samples.empty()) {
        throw std::invalid_argument("compute_exposure_stats: V_samples empty");
    }
    if (cfg.confidence_level <= 0.0 || cfg.confidence_level >= 1.0) {
        throw std::invalid_argument("compute_exposure_stats: confidence_level must be in (0, 1)");
    }
    const Size n_paths = V_samples.size();
    const Size n_times = times.size();
    for (const auto& v : V_samples) {
        if (v.size() != n_times) {
            throw std::invalid_argument("compute_exposure_stats: column size mismatch");
        }
    }

    ExposureStats stats;
    stats.times = times;
    stats.ee.assign(n_times, 0.0);
    stats.ene.assign(n_times, 0.0);
    stats.eee.assign(n_times, 0.0);
    stats.pfe.assign(n_times, 0.0);

    // 排序索引, 用于分位数计算
    Real alpha = cfg.confidence_level;
    Real running_max_ee = 0.0;
    Real sum_ee = 0.0;
    Real sum_eee = 0.0;
    stats.max_pfe = 0.0;
    stats.max_ee = 0.0;

    for (Size j = 0; j < n_times; ++j) {
        Real sum_pos = 0.0;
        Real sum_neg = 0.0;
        std::vector<Real> exposures_pos;
        exposures_pos.reserve(n_paths);
        Real P_d = discount_curve.discount_factor(times[j]);

        for (Size p = 0; p < n_paths; ++p) {
            Real v = V_samples[p][j] * P_d;  // 折现暴露
            Real e = std::max(v, 0.0);
            Real ne = std::max(-v, 0.0);
            sum_pos += e;
            sum_neg += ne;
            exposures_pos.push_back(e);
        }
        Real ee = sum_pos / static_cast<Real>(n_paths);
        Real ene = sum_neg / static_cast<Real>(n_paths);
        stats.ee[j] = ee;
        stats.ene[j] = ene;

        // PFE (α 分位数)
        std::sort(exposures_pos.begin(), exposures_pos.end());
        // 经验 α 分位数: index = α * (n - 1)
        Real idx_f = alpha * static_cast<Real>(n_paths - 1);
        Size idx_lo = static_cast<Size>(std::floor(idx_f));
        Size idx_hi = static_cast<Size>(std::ceil(idx_f));
        if (idx_hi >= n_paths) idx_hi = n_paths - 1;
        Real w = idx_f - static_cast<Real>(idx_lo);
        Real pfe = (1.0 - w) * exposures_pos[idx_lo] + w * exposures_pos[idx_hi];
        stats.pfe[j] = pfe;

        // EEE (running max)
        if (cfg.use_running_max_ee) {
            running_max_ee = std::max(running_max_ee, ee);
        } else {
            running_max_ee = ee;
        }
        stats.eee[j] = running_max_ee;

        // 累积时间平均 (trapezoidal 需要区间, 这里用简单平均的等价形式)
        if (j > 0) {
            Real dt = times[j] - times[j - 1];
            if (dt > 0.0) {
                sum_ee += 0.5 * (stats.ee[j - 1] + stats.ee[j]) * dt;
                sum_eee += 0.5 * (stats.eee[j - 1] + stats.eee[j]) * dt;
            }
        }
        stats.max_ee = std::max(stats.max_ee, ee);
        stats.max_pfe = std::max(stats.max_pfe, pfe);
    }

    // EPE/EEPE: 时间平均 (除以总时长)
    Real T = times.back() - times.front();
    if (T > 0.0) {
        stats.epe = sum_ee / T;
        stats.eepe = sum_eee / T;
    }
    return stats;
}

// ============ 端到端 MC PFE ============
// value_fn(t, S_t) 返回衍生品组合在 t 时刻的市值
inline ExposureStats compute_pfe_mc(
        const MultiAssetGBMPathGenerator& gen,
        std::function<Real(Real, const std::vector<Real>&)> value_fn,
        const std::vector<Real>& exposure_times,
        const ZeroCurve& discount_curve,
        Size n_paths,
        uint64_t seed,
        const PFEConfig& cfg = {}) {
    if (n_paths == 0) {
        throw std::invalid_argument("compute_pfe_mc: n_paths must be positive");
    }
    if (exposure_times.empty()) {
        throw std::invalid_argument("compute_pfe_mc: exposure_times empty");
    }
    const Real dt_gen = gen.dt();
    const Real T_gen = gen.config().T;
    const Size n_steps = gen.config().n_steps;

    std::vector<Size> step_indices;
    step_indices.reserve(exposure_times.size());
    for (Real t : exposure_times) {
        if (t < 0.0 || t > T_gen + 1e-9) {
            throw std::invalid_argument("compute_pfe_mc: exposure_time out of range");
        }
        Size idx = static_cast<Size>(std::round(t / dt_gen));
        if (idx > n_steps) idx = n_steps;
        step_indices.push_back(idx);
    }

    std::vector<std::vector<Real>> V_samples(n_paths,
                                                std::vector<Real>(exposure_times.size()));
    for (Size p = 0; p < n_paths; ++p) {
        Philox4x64 rng(seed, static_cast<uint64_t>(p));
        auto paths = gen.generate_path(rng, 1.0);
        for (Size j = 0; j < exposure_times.size(); ++j) {
            Size step = step_indices[j];
            std::vector<Real> S_at_t(paths.size());
            for (Size a = 0; a < paths.size(); ++a) {
                S_at_t[a] = paths[a][step];
            }
            V_samples[p][j] = value_fn(exposure_times[j], S_at_t);
        }
    }

    return compute_exposure_stats(V_samples, exposure_times, discount_curve, cfg);
}

// ============ SA-CCR (Basel III 标准化方法) ============

// 资产类别 (5 类 Basel 分类)
enum class AssetClass {
    InterestRate,   // 利率
    FX,             // 外汇
    Credit,         // 信用 (单名/指数)
    Equity,         // 股票
    Commodity       // 商品
};

// Supervisory Factor (RF) 默认值 (BCBS d279 表 1)
// 注意: 信用类需按 rating 调整, 股票需按 index/single 调整
inline Real supervisory_factor(AssetClass ac) {
    switch (ac) {
        case AssetClass::InterestRate: return 0.005;   // 0.5%
        case AssetClass::FX:           return 0.04;    // 4.0%
        case AssetClass::Credit:       return 0.05;    // 5.0% (IG 默认, 评级调整见下)
        case AssetClass::Equity:       return 0.32;    // 32% (single name 默认)
        case AssetClass::Commodity:    return 0.18;    // 18% (默认)
    }
    return 0.05;
}

// 信用类 supervisory factor 按 rating (BCBS d279 表 1)
// rating: AAA/A: 0.0038; BBB: 0.01; BB: 0.06; B: 0.10; CCC: 0.15; 非违约 unspecified IG: 0.05
inline Real credit_supervisory_factor_by_rating(const std::string& rating) {
    if (rating == "AAA" || rating == "AA" || rating == "A") return 0.0038;
    if (rating == "BBB") return 0.01;
    if (rating == "BB")  return 0.06;
    if (rating == "B")   return 0.10;
    if (rating == "CCC") return 0.15;
    return 0.05;  // default IG
}

// 股票类 supervisory factor: index 12%, single name 32%
inline Real equity_supervisory_factor(bool is_index) {
    return is_index ? 0.12 : 0.32;
}

// ============ SA-CCR 交易 (单笔) ============
struct SACCRTrade {
    Real notional = 0.0;            // 名义本金
    Real remaining_maturity = 1.0;  // 剩余期限 (年)
    AssetClass asset_class = AssetClass::InterestRate;
    Real price = 1.0;               // 当前价格 (用于 notional = contract_notional * price, FX/股票类)
    bool is_index_equity = false;   // 股票类: 是否为 index
    std::string credit_rating;      // 信用类: rating (AAA/AA/A/BBB/BB/B/CCC)
    Real long_short_sign = 1.0;     // +1: 多头, -1: 空头 (用于 netting)
};

// ============ SA-CCR 投资组合 (单 hedging set) ============
struct SACCRPortfolio {
    std::vector<SACCRTrade> trades;
    Real net_collateral = 0.0;       // 净抵押品 (C, 银行已收抵押 - 已付)
    Real portfolio_value = 0.0;      // 组合市值 V (银行视角, 正=收钱)
    Real correlation = 0.50;         // supervisory correlation (默认 50%)
    // mul_addon: supervisory multiplier, 默认 1.4 (BCBS d279 §143)
    Real supervisory_multiplier = 1.4;

    // Net-to-Gross Ratio (NGR) per hedging set
    // NGR = Σ|net notional| / Σ|gross notional|
    Real compute_ngr() const {
        Real gross = 0.0, net = 0.0;
        for (const auto& t : trades) {
            gross += std::abs(t.notional * t.price);
        }
        // net: 按 asset_class 分组求和后取绝对值 (简化: 不分组)
        std::vector<Real> by_class_long(5, 0.0), by_class_short(5, 0.0);
        for (const auto& t : trades) {
            int idx = static_cast<int>(t.asset_class);
            Real n = t.notional * t.price * t.long_short_sign;
            if (n >= 0.0) by_class_long[idx] += n;
            else by_class_short[idx] += -n;
        }
        for (int i = 0; i < 5; ++i) {
            Real l = by_class_long[i];
            Real s = by_class_short[i];
            net += std::abs(l - s);
        }
        if (gross <= 0.0) return 0.0;
        return net / gross;
    }
};

// ============ Maturity Factor (MF) ============
inline Real maturity_factor(Real remaining_maturity) {
    if (remaining_maturity <= 0.0) return 0.0;
    if (remaining_maturity >= 1.0) return 1.0;
    return std::sqrt(remaining_maturity);  // BCBS 2019: < 1Y 用 √(T/1Y)
}

// ============ SA-CCR 计算 ============
struct SACCRResult {
    Real rc = 0.0;             // Replacement Cost
    Real pfe_addon = 0.0;      // PFE Addon
    Real ead = 0.0;            // Exposure at Default
    Real mul_addon = 0.0;      // Supervisory Multiplier (1.4)
    Real aggregate_addon = 0.0; // 加总 addon (mul 之前)
    Real ngr = 0.0;            // Net-to-Gross Ratio
};

// 简化 SA-CCR: 单 hedging set, 不区分 hedging set level aggregation
// 适用于教学/快速估算; 完整 SA-CCR 见 BCBS d279 §157-187
inline SACCRResult compute_sa_ccr(const SACCRPortfolio& portfolio) {
    SACCRResult r;
    r.mul_addon = portfolio.supervisory_multiplier;
    r.ngr = portfolio.compute_ngr();

    // Replacement Cost
    Real V = portfolio.portfolio_value;
    Real C = portfolio.net_collateral;
    r.rc = std::max(V - C, 0.0);

    // Aggregate Addon: Σ_i (notional_i × RF_i × MF_i × EF_i)
    // EF (Effective Notional) = d * SF * MF (简化版, 不区分 hedging set)
    // 这里按 asset_class 分组, 每组内做 netting (long-short 抵消)
    std::vector<Real> addon_by_class(5, 0.0);
    std::vector<Real> gross_long(5, 0.0), gross_short(5, 0.0);
    for (const auto& t : portfolio.trades) {
        int idx = static_cast<int>(t.asset_class);
        Real effective_notional = t.notional * t.price * maturity_factor(t.remaining_maturity);
        Real signed_n = effective_notional * t.long_short_sign;
        if (signed_n >= 0.0) gross_long[idx] += signed_n;
        else gross_short[idx] += -signed_n;
    }

    for (int i = 0; i < 5; ++i) {
        if (gross_long[i] == 0.0 && gross_short[i] == 0.0) continue;
        AssetClass ac = static_cast<AssetClass>(i);
        Real RF = supervisory_factor(ac);
        if (ac == AssetClass::Credit) {
            // 信用类: 取第一个信用交易的 rating (简化)
            for (const auto& t : portfolio.trades) {
                if (t.asset_class == AssetClass::Credit && !t.credit_rating.empty()) {
                    RF = credit_supervisory_factor_by_rating(t.credit_rating);
                    break;
                }
            }
        }
        if (ac == AssetClass::Equity) {
            bool any_index = false;
            for (const auto& t : portfolio.trades) {
                if (t.asset_class == AssetClass::Equity && t.is_index_equity) {
                    any_index = true;
                    break;
                }
            }
            RF = equity_supervisory_factor(any_index);
        }

        // Hedging set 内 netting: addon = RF × max(net, gross × 0.4) (supervisory floor)
        // BCBS §158: addon = SF × |net + gross × 0.4| (simplified)
        // 完整公式: addon = SF × sqrt((sum long² + sum short² + 2·ρ·sum_long·sum_short))
        // 这里用 supervisory floor 简化版: addon = RF × max(|L-S|, 0.4 × (L+S))
        Real L = gross_long[i];
        Real S = gross_short[i];
        Real net = std::abs(L - S);
        Real gross = L + S;
        Real effective_notional = std::max(net, 0.4 * gross);
        addon_by_class[i] = RF * effective_notional;
    }

    // Aggregate across asset classes (假设不相关, 直接求和)
    r.aggregate_addon = std::accumulate(addon_by_class.begin(), addon_by_class.end(), 0.0);

    // Apply NGR adjustment (supervisory factor with net/gross)
    // BCBS §149: PFE_addon = mul_addon × aggregate_addon × (1 + 0.4 × NGR) / (1 + 0.4)
    // (NGR 已隐含在 hedging set addon 中, 这里使用 1.4× 直接计算)
    r.pfe_addon = r.mul_addon * r.aggregate_addon;

    // BCBS §161: EAD = α × (RC + PFE_addon), α = 1.4 (supervisory factor)
    r.ead = 1.4 * (r.rc + r.pfe_addon);

    return r;
}

// ============ 便捷工厂: IRS Trade (用于 SA-CCR) ============
inline SACCRTrade make_sa_ccr_irs_trade(Real notional, Real remaining_maturity,
                                          Real long_short_sign = 1.0) {
    SACCRTrade t;
    t.notional = notional;
    t.remaining_maturity = remaining_maturity;
    t.asset_class = AssetClass::InterestRate;
    t.price = 1.0;
    t.long_short_sign = long_short_sign;
    return t;
}

inline SACCRTrade make_sa_ccr_fx_trade(Real notional, Real remaining_maturity,
                                          Real long_short_sign = 1.0) {
    SACCRTrade t;
    t.notional = notional;
    t.remaining_maturity = remaining_maturity;
    t.asset_class = AssetClass::FX;
    t.price = 1.0;
    t.long_short_sign = long_short_sign;
    return t;
}

inline SACCRTrade make_sa_ccr_credit_trade(Real notional, Real remaining_maturity,
                                              const std::string& rating,
                                              Real long_short_sign = 1.0) {
    SACCRTrade t;
    t.notional = notional;
    t.remaining_maturity = remaining_maturity;
    t.asset_class = AssetClass::Credit;
    t.price = 1.0;
    t.credit_rating = rating;
    t.long_short_sign = long_short_sign;
    return t;
}

inline SACCRTrade make_sa_ccr_equity_trade(Real notional, Real remaining_maturity,
                                              bool is_index, Real long_short_sign = 1.0) {
    SACCRTrade t;
    t.notional = notional;
    t.remaining_maturity = remaining_maturity;
    t.asset_class = AssetClass::Equity;
    t.price = 1.0;
    t.is_index_equity = is_index;
    t.long_short_sign = long_short_sign;
    return t;
}

inline SACCRTrade make_sa_ccr_commodity_trade(Real notional, Real remaining_maturity,
                                                  Real long_short_sign = 1.0) {
    SACCRTrade t;
    t.notional = notional;
    t.remaining_maturity = remaining_maturity;
    t.asset_class = AssetClass::Commodity;
    t.price = 1.0;
    t.long_short_sign = long_short_sign;
    return t;
}

}  // namespace v1
}  // namespace cpphub
