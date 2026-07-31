#pragma once
// SOURCE: Hull (2018) Ch.27 (Multi-Asset Options)
// SOURCE: Margrabe (1978) "The Value of an Option to Exchange One Asset for Another"
// SOURCE: Stulz (1982) "Options on the Minimum or the Maximum of Two Risky Assets"
// 模块: 多资产 (篮子/彩虹/价差) payoff
//
// 约定:
//   paths[asset_idx][step_idx] - n_assets × (n_steps+1) 矩阵
//   terminal price S_i(T) = paths[i].back()
//   仅 terminal payoff (path-independent multi-asset); 路径相关篮子由组合方式构造

#include "cpphub/core/types.hpp"
#include "cpphub/pricing/monte_carlo/path_dependent_payoffs.hpp"  // OptionType
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <numeric>

namespace cpphub {
inline namespace v1 {

// ============ 多资产 payoff 函数签名 ============
// 接收 [n_assets][n_steps+1] 路径矩阵, 返回 payoff (未折现)
using MultiAssetPayoff = std::function<Real(const std::vector<std::vector<Real>>&)>;

// ============ 1. Basket Option ============
// 篮子期权: max(Σ w_i S_i(T) - K, 0) (call) / max(K - Σ w_i S_i(T), 0) (put)
inline MultiAssetPayoff make_basket_payoff(const std::vector<Real>& weights,
                                             Real K, OptionType opt) {
    if (weights.empty()) throw std::invalid_argument("Basket: weights empty");
    if (K < 0.0) throw std::invalid_argument("Basket: K must be non-negative");
    Real wsum = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (std::abs(wsum - 1.0) > 1e-6)
        throw std::invalid_argument("Basket: weights must sum to 1");
    return [weights, K, opt](const std::vector<std::vector<Real>>& paths) -> Real {
        if (paths.size() != weights.size())
            throw std::invalid_argument("Basket: paths/weights size mismatch");
        Real basket = 0.0;
        for (Size i = 0; i < weights.size(); ++i) {
            if (paths[i].empty()) throw std::invalid_argument("Basket: path empty");
            basket += weights[i] * paths[i].back();
        }
        return (opt == OptionType::Call)
                ? std::max(basket - K, 0.0)
                : std::max(K - basket, 0.0);
    };
}

// ============ 2. Rainbow Options (Stulz 1982) ============
enum class RainbowType { BestOf, WorstOf, Nth };  // max / min / n-th

// Rainbow call on the max/min of n assets
// BestOf call: max(max_i S_i(T) - K, 0)
// WorstOf call: max(min_i S_i(T) - K, 0)
inline MultiAssetPayoff make_rainbow_payoff(Real K, OptionType opt,
                                              RainbowType rtype, Size nth = 0) {
    if (K < 0.0) throw std::invalid_argument("Rainbow: K must be non-negative");
    return [K, opt, rtype, nth](const std::vector<std::vector<Real>>& paths) -> Real {
        if (paths.empty()) throw std::invalid_argument("Rainbow: paths empty");
        std::vector<Real> terminals;
        terminals.reserve(paths.size());
        for (const auto& p : paths) {
            if (p.empty()) throw std::invalid_argument("Rainbow: path empty");
            terminals.push_back(p.back());
        }
        Real extremum;
        if (rtype == RainbowType::BestOf) {
            extremum = *std::max_element(terminals.begin(), terminals.end());
        } else if (rtype == RainbowType::WorstOf) {
            extremum = *std::min_element(terminals.begin(), terminals.end());
        } else {
            // Nth: sort and pick (0-indexed)
            std::sort(terminals.begin(), terminals.end());
            Size idx = std::min(nth, paths.size() - 1);
            extremum = terminals[idx];
        }
        return (opt == OptionType::Call)
                ? std::max(extremum - K, 0.0)
                : std::max(K - extremum, 0.0);
    };
}

// ============ 3. Spread Option (Margrabe 1978 generalized) ============
// Spread call: max(S_1(T) - S_2(T) - K, 0)
// K=0 退化为 Margrabe exchange option (S_1 换 S_2)
inline MultiAssetPayoff make_spread_payoff(Real K, OptionType opt,
                                             Size idx1 = 0, Size idx2 = 1) {
    if (K < 0.0) throw std::invalid_argument("Spread: K must be non-negative");
    return [K, opt, idx1, idx2](const std::vector<std::vector<Real>>& paths) -> Real {
        if (paths.size() <= std::max(idx1, idx2))
            throw std::invalid_argument("Spread: index out of range");
        Real S1 = paths[idx1].back();
        Real S2 = paths[idx2].back();
        Real spread = S1 - S2;
        return (opt == OptionType::Call)
                ? std::max(spread - K, 0.0)
                : std::max(K - spread, 0.0);
    };
}

// ============ 4. Multi-asset digital (all-itm / any-itm) ============
enum class DigitalCombo { AllITM, AnyITM };

// 多资产 digital: AllITM = 所有 asset 都 ITM 才付 1; AnyITM = 任一 ITM 付 1
inline MultiAssetPayoff make_multi_digital_payoff(const std::vector<Real>& strikes,
                                                    const std::vector<OptionType>& opts,
                                                    DigitalCombo combo,
                                                    Real payoff = 1.0) {
    if (strikes.size() != opts.size())
        throw std::invalid_argument("MultiDigital: size mismatch");
    return [strikes, opts, combo, payoff](const std::vector<std::vector<Real>>& paths) -> Real {
        if (paths.size() != strikes.size())
            throw std::invalid_argument("MultiDigital: paths/strikes size mismatch");
        bool all_itm = true, any_itm = false;
        for (Size i = 0; i < paths.size(); ++i) {
            Real S_T = paths[i].back();
            bool itm = (opts[i] == OptionType::Call) ? (S_T > strikes[i]) : (S_T < strikes[i]);
            if (itm) any_itm = true; else all_itm = false;
        }
        bool pay = (combo == DigitalCombo::AllITM) ? all_itm : any_itm;
        return pay ? payoff : 0.0;
    };
}

}  // namespace v1
}  // namespace cpphub
