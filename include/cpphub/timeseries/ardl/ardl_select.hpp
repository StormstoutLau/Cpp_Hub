// =============================================================================
// ardl_select.hpp - ARDL IC 网格阶数选择, 同一 hold_back 样本 (spec §2.2)
//
// Phase 8 v1.8 M1 Step 3 (PHASE8_SPEC.md v1.1 §2.2)
//
// 一手源码语义落档 (statsmodels 0.14.6 tsa/ardl/model.py ardl_select_order
// L1452-1570 实录, verify_ardl.py sel_k1/sel_k2 基准逐值自证):
//   - 网格 = (maxlag+1) × (maxorder+1)^k 联合候选 (p ∈ 0..maxlag,
//     q_j ∈ 0..maxorder 每变量独立, uniform int; per-variable dict 不做,
//     手动 ardl_fit 替代)
//   - ⚠️ statsmodels 搜索域更宽: 其 iter_orders = range(block_cols + 1)
//     含 io_j = 0 候选 = 变量 j 完全排除 (连 L0 都不含, q_j=None) →
//     实际域 (maxlag+1)·(maxorder+2)^k ⊃ 本库域. A3 裁决: 变量排除不做
//     (需要排除时手动 ardl_fit); 若 statsmodels 最优落在排除区, 本库
//     argmin 结果可能与其偏离 — 已知且接受
//   - 同一样本 (A3 核心): 全网格统一 hold_back_grid = max(maxlag,
//     maxorder) — 各候选仅列数不同, 估计样本恒 T − hb_grid (可比性,
//     镜像 var_select offset 机制)
//   - 网格 IC 计账 (compute_ics L1454-1470 实录): lstsq OLS,
//     n = n_act = T − hb_grid, σ² = SSR/n, llf = −n/2·(ln2π + lnσ² + 1),
//     IC 复用 ARDLResults.aic/bic/hqic → 惩罚恒 (K+1) (σ² 计入).
//     本库实现: ardl_fit 显式 hold_back = hb_grid → A8 计账 hb_early =
//     显式值 → n_acct = T − hb_grid = n_act, IC 字段逐值 = 网格口径
//     (ardl_model.hpp A8 机制在此恰与 compute_ics 重合)
//   - argmin: 严格 <, 先到先得; 扫描序 = product 序 (p 最外层, 最后
//     变量阶最快) — 平局语义与 statsmodels dict 插入序一致
//   - refit: sel.model 无显式 hb → fit() 走 A8 计账 n_acct = T − p_best,
//     hb_eff = max(p_best, max q_best) — 与网格口径 (n = T − hb_grid) 的
//     差异 = "最优模型自有样本 vs 统一网格样本" (sel_k1: 249 vs 247,
//     ic 值不同; sel_k2: p_best=2=hb_grid → 两口径重合). 本库 best_model
//     = ardl_fit(默认 hb) 逐值镜像 ic_value_refit 基准
//   - hb_grid=0 edge (maxlag=maxorder=0): 显式 0 在 ardl_fit 走默认分支
//     max(p, maxq) = 0, hb_early = p = 0 → n_acct = T = n_act, 网格口径
//     仍一致 (该 edge 下网格唯一候选 p=0, q_j=0)
// =============================================================================

#pragma once

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/timeseries/ardl/ardl_model.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace ardl {

/// IC 网格阶数选择结果 (spec §2.2 逐字段)
struct ArdlSelectResult {
    std::string ic_used = "bic";
    Size p_opt = 0;                        ///< 最优 y 滞后阶
    std::vector<Size> q_opt;               ///< 最优各变量阶 (联合网格 argmin)
    Real ic_min = 0.0;                     ///< 网格口径最优 IC (统一 hb, = statsmodels sel.<ic> 首值)
    std::vector<Real> aic_grid, bic_grid, hqic_grid;  ///< 扁平联合网格 (p 外层, 末变量最快)
    Size n_p_grid = 0;                     ///< maxlag + 1
    std::vector<Size> q_grid_dims;         ///< 各变量维度 (maxorder+1, uniform)
    Size n_obs_used = 0;                   ///< 全网格统一 = T − hb_grid
    ArdlResult best_model{};               ///< refit (默认 hb → A8 计账口径)
};

/// ARDL IC 网格阶数选择 (同一 hold_back 样本, A3)
///
/// 网格 = (maxlag+1) × (maxorder+1)^k; glob (lag 子集搜索) 不做 (A3);
/// 变量完全排除候选不做 (statsmodels 独有, 见头注释 — 需要时手动 ardl_fit)
inline ArdlSelectResult ardl_select_order(
    const std::vector<Real>& y,
    const std::vector<std::vector<Real>>& x,
    Size maxlag, Size maxorder, const std::string& ic = "bic",
    const std::string& trend = "c") {
    const bool use_aic = (ic == "aic");
    const bool use_bic = (ic == "bic");
    const bool use_hqic = (ic == "hqic");
    if (!(use_aic || use_bic || use_hqic)) {
        throw std::invalid_argument(
            "ardl_select_order: ic must be one of aic/bic/hqic");
    }
    detail::trend_order(trend);  // trend 非法 → throw

    const Size T = y.size();
    const Size k = x.size();
    if (T == 0) {
        throw std::invalid_argument("ardl_select_order: empty endog");
    }
    for (Size j = 0; j < k; ++j) {
        if (x[j].size() != T) {
            throw std::invalid_argument(
                "ardl_select_order: x[j] length != y length");
        }
    }

    // 同一样本 (A3): 全网格统一 hb (含 ardl_fit 内 NaN/k>10/样本过小校验,
    // 任一候选不合法 → 整个选择 throw)
    const Size hb_grid = std::max(maxlag, maxorder);
    const Size nq = maxorder + 1;

    ArdlSelectResult res;
    res.ic_used = ic;
    res.n_p_grid = maxlag + 1;
    res.q_grid_dims.assign(k, nq);
    res.n_obs_used = T - hb_grid;

    Size total = maxlag + 1;
    for (Size j = 0; j < k; ++j) total *= nq;
    res.aic_grid.reserve(total);
    res.bic_grid.reserve(total);
    res.hqic_grid.reserve(total);

    // 扁平索引 → (p, q_0..q_{k-1}): 末变量最低位 (变化最快), p 最高位
    std::vector<Size> qs(k);
    Real best = std::numeric_limits<Real>::infinity();
    Size best_idx = 0;
    for (Size idx = 0; idx < total; ++idx) {
        Size rem = idx;
        for (Size j = k; j-- > 0;) {
            qs[j] = rem % nq;
            rem /= nq;
        }
        const Size p = rem;  // rem ∈ 0..maxlag

        ArdlSpec s;
        s.p = p;
        s.q = qs;
        s.trend = trend;
        s.hold_back = hb_grid;  // A8: 显式 hb → n_acct = n_act (网格口径)
        const ArdlResult r = ardl_fit(y, x, s);
        res.aic_grid.push_back(r.aic);
        res.bic_grid.push_back(r.bic);
        res.hqic_grid.push_back(r.hqic);

        const Real val = use_aic ? r.aic : (use_bic ? r.bic : r.hqic);
        if (val < best) {  // 严格 <, 先到先得 (product 序, 平局同 statsmodels)
            best = val;
            best_idx = idx;
        }
    }

    Size rem = best_idx;
    res.q_opt.assign(k, 0);
    for (Size j = k; j-- > 0;) {
        res.q_opt[j] = rem % nq;
        rem /= nq;
    }
    res.p_opt = rem;
    res.ic_min = best;

    // refit: 默认 hold_back → A8 计账 (n_acct = T − p_best) = statsmodels
    // ic_value_refit 口径 (hb_eff = max(p_best, max q_best))
    ArdlSpec s;
    s.p = res.p_opt;
    s.q = res.q_opt;
    s.trend = trend;
    res.best_model = ardl_fit(y, x, s);
    return res;
}

}  // namespace ardl
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
