#pragma once
// SOURCE: Hull (2018) "Options, Futures, and Other Derivatives" Ch.29
// SOURCE: Brigo & Mercurio (2006) "Interest Rate Models - Theory and Practice" Ch.1, Ch.3
// 模块: 利率上限/下限 (Interest Rate Cap / Floor)
//
// Cap = 一系列 caplet 的组合; 每个 caplet 在 T_i 支付 τ_i * max(L(T_{i-1}, T_i) - K, 0) * N
// Floor = 一系列 floorlet 的组合; 每个 floorlet 在 T_i 支付 τ_i * max(K - L(T_{i-1}, T_i), 0) * N
//
// LIBOR 远期利率 (简单复利):
//   L(T_{i-1}, T_i) = (1/τ_i) * (1/P(T_{i-1}, T_i) - 1)
//
// Caplet/Floorlet 等价于零息债期权:
//   caplet(T_i) = max(1 - (1+Kτ) * P(T_{i-1}, T_i), 0)         [payoff at T_i]
//   floorlet(T_i) = max((1+Kτ) * P(T_{i-1}, T_i) - 1, 0)       [payoff at T_i]
//
//   等价零息债期权 (strike K_eff = (1+Kτ)^-1):
//   caplet = (1+Kτ) * Put_on_P(T_opt=T_{i-1}, T_bond=T_i, K=(1+Kτ)^-1)
//   floorlet = (1+Kτ) * Call_on_P(T_opt=T_{i-1}, T_bond=T_i, K=(1+Kτ)^-1)
//
// 两种定价方法:
//   1. 模型闭式 (HW/Vasicek): 用 bond_option(T_opt, T_bond, K_eff, is_call)
//   2. Black 76 市场报价 (假设远期利率对数正态):
//      caplet = τ * N * P(0, T_i) * (F * N(d1) - K * N(d2))
//      d1 = (log(F/K) + 0.5 σ² T_{i-1}) / (σ sqrt(T_{i-1}))
//      d2 = d1 - σ sqrt(T_{i-1})
//      F = (P(0,T_{i-1})/P(0,T_i) - 1) / τ   (远期 LIBOR)

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/models/ir/short_rate.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ Cap/Floor 配置 ============
struct CapFloorConfig {
    Real notional = 1.0;
    Real strike;                       // K (年化上限/下限利率)
    std::vector<Real> reset_times;       // T_0, T_1, ..., T_{n-1} (重置日, 严格递增)
    std::vector<Real> payment_times;     // T_1, T_2, ..., T_n (支付日, 严格递增)
    std::vector<Real> year_fractions;    // τ_i = T_i - T_{i-1}
    bool is_cap = true;                // true: Cap, false: Floor

    Size n_caplets() const { return reset_times.size(); }

    void validate() const {
        if (notional <= 0.0) throw std::invalid_argument("CapFloor: notional must be positive");
        if (strike < 0.0) throw std::invalid_argument("CapFloor: strike must be non-negative");
        if (reset_times.size() != payment_times.size()) {
            throw std::invalid_argument("CapFloor: reset/payment size mismatch");
        }
        if (reset_times.size() != year_fractions.size()) {
            throw std::invalid_argument("CapFloor: reset/year_fractions size mismatch");
        }
        if (reset_times.empty()) {
            throw std::invalid_argument("CapFloor: at least 1 caplet required");
        }
        for (Size i = 0; i < reset_times.size(); ++i) {
            if (reset_times[i] < 0.0) throw std::invalid_argument("CapFloor: reset_times must be non-negative");
            if (payment_times[i] <= reset_times[i]) {
                throw std::invalid_argument("CapFloor: payment_times must be > reset_times");
            }
            if (year_fractions[i] <= 0.0) {
                throw std::invalid_argument("CapFloor: year_fractions must be positive");
            }
            if (i > 0) {
                if (reset_times[i] <= reset_times[i - 1]) {
                    throw std::invalid_argument("CapFloor: reset_times must be strictly increasing");
                }
                if (payment_times[i] <= payment_times[i - 1]) {
                    throw std::invalid_argument("CapFloor: payment_times must be strictly increasing");
                }
            }
        }
    }
};

// ============ Cap/Floor 定价器 ============
class CapFloor {
public:
    explicit CapFloor(CapFloorConfig cfg) : cfg_(std::move(cfg)) {
        cfg_.validate();
    }

    // 模型闭式 (HW/Vasicek): 每个 caplet/floorlet 用 bond_option
    // caplet = (1+Kτ) * bond_option(T_{i-1}, T_i, (1+Kτ)^-1, is_call=false)
    // floorlet = (1+Kτ) * bond_option(T_{i-1}, T_i, (1+Kτ)^-1, is_call=true)
    // bond_option 接口在 Vasicek 和 HullWhite 中都有
    template <typename Model>
    Real present_value_model(const Model& model) const {
        Real pv = 0.0;
        for (Size i = 0; i < cfg_.n_caplets(); ++i) {
            Real T_reset = cfg_.reset_times[i];
            Real T_pay = cfg_.payment_times[i];
            Real tau = cfg_.year_fractions[i];
            Real K_eff = 1.0 / (1.0 + cfg_.strike * tau);  // (1+Kτ)^-1
            Real mult = 1.0 + cfg_.strike * tau;            // (1+Kτ)
            // Cap = put on bond; Floor = call on bond
            bool is_call = !cfg_.is_cap;
            Real caplet_pv;
            if (T_reset <= 0.0) {
                // 边界: T_reset=0 时远期利率已确定, caplet = 内在价值
                // Payoff at T_pay = max(±(1/P(0,T_pay) - (1+Kτ)), 0)
                // PV = P(0,T_pay) * payoff = max(±(1 - (1+Kτ)*P(0,T_pay)), 0)  (P>0 提出)
                Real P_pay = model.zero_coupon_bond(T_pay);
                caplet_pv = cfg_.is_cap
                    ? std::max(1.0 - mult * P_pay, 0.0)
                    : std::max(mult * P_pay - 1.0, 0.0);
            } else {
                caplet_pv = mult * model.bond_option(T_reset, T_pay, K_eff, is_call);
            }
            pv += cfg_.notional * caplet_pv;
        }
        return pv;
    }

    // Black 76 形式 (市场报价, 假设远期利率对数正态)
    // caplet = τ * N * P(0, T_i) * (F * N(d1) - K * N(d2))
    // floorlet = τ * N * P(0, T_i) * (K * N(-d2) - F * N(-d1))
    // vol: caplet/floorlet 波动率 (年化)
    template <typename Model>
    Real present_value_black(const Model& model, Real vol) const {
        if (vol < 0.0) throw std::invalid_argument("CapFloor::black: vol must be non-negative");
        Real pv = 0.0;
        for (Size i = 0; i < cfg_.n_caplets(); ++i) {
            Real T_reset = cfg_.reset_times[i];
            Real T_pay = cfg_.payment_times[i];
            Real tau = cfg_.year_fractions[i];
            Real P_reset = model.zero_coupon_bond(T_reset);
            Real P_pay = model.zero_coupon_bond(T_pay);
            // 远期 LIBOR (简单复利)
            Real F = (P_reset / P_pay - 1.0) / tau;
            Real df = P_pay;  // 折现到 t=0 (payoff 在 T_pay 支付)
            Real caplet;
            if (vol < 1e-12 || T_reset < 1e-12) {
                // 退化: 直接内在价值
                Real intrinsic = cfg_.is_cap ? std::max(F - cfg_.strike, 0.0)
                                              : std::max(cfg_.strike - F, 0.0);
                caplet = tau * cfg_.notional * df * intrinsic;
            } else {
                Real sqrt_T = std::sqrt(T_reset);
                Real d1 = (std::log(F / cfg_.strike) + 0.5 * vol * vol * T_reset)
                          / (vol * sqrt_T);
                Real d2 = d1 - vol * sqrt_T;
                if (cfg_.is_cap) {
                    caplet = tau * cfg_.notional * df
                             * (F * normal_cdf(d1) - cfg_.strike * normal_cdf(d2));
                } else {
                    caplet = tau * cfg_.notional * df
                             * (cfg_.strike * normal_cdf(-d2) - F * normal_cdf(-d1));
                }
            }
            pv += caplet;
        }
        return pv;
    }

    // 计算隐含 Black 波动率 (给定模型闭式价格, 反推 vol)
    // 用二分法
    template <typename Model>
    Real implied_black_vol(const Model& model, Real target_price,
                            Real vol_low = 1e-4, Real vol_high = 5.0,
                            Real tol = 1e-8, Size max_iter = 100) const {
        if (target_price <= 0.0) return 0.0;
        Real p_low = present_value_black(model, vol_low);
        Real p_high = present_value_black(model, vol_high);
        if (target_price <= p_low) return vol_low;
        if (target_price >= p_high) return vol_high;
        for (Size iter = 0; iter < max_iter; ++iter) {
            Real vol_mid = 0.5 * (vol_low + vol_high);
            Real p_mid = present_value_black(model, vol_mid);
            if (std::abs(p_mid - target_price) < tol) return vol_mid;
            if (p_mid < target_price) {
                vol_low = vol_mid;
                p_low = p_mid;
            } else {
                vol_high = vol_mid;
                p_high = p_mid;
            }
        }
        return 0.5 * (vol_low + vol_high);
    }

    const CapFloorConfig& config() const { return cfg_; }
    Real notional() const { return cfg_.notional; }
    Real strike() const { return cfg_.strike; }
    bool is_cap() const { return cfg_.is_cap; }

private:
    CapFloorConfig cfg_;
};

// ============ 便捷工厂: 等间隔 Cap/Floor ============
inline CapFloorConfig make_vanilla_capfloor(Real notional, Real strike, Real start_time,
                                              Real tau, Size n_caplets, bool is_cap = true) {
    if (tau <= 0.0) throw std::invalid_argument("make_vanilla_capfloor: tau must be positive");
    if (n_caplets == 0) throw std::invalid_argument("make_vanilla_capfloor: n_caplets must be > 0");
    CapFloorConfig cfg;
    cfg.notional = notional;
    cfg.strike = strike;
    cfg.is_cap = is_cap;
    cfg.reset_times.resize(n_caplets);
    cfg.payment_times.resize(n_caplets);
    cfg.year_fractions.resize(n_caplets, tau);
    for (Size i = 0; i < n_caplets; ++i) {
        cfg.reset_times[i] = start_time + static_cast<Real>(i) * tau;
        cfg.payment_times[i] = start_time + static_cast<Real>(i + 1) * tau;
    }
    return cfg;
}

}  // namespace v1
}  // namespace cpphub
