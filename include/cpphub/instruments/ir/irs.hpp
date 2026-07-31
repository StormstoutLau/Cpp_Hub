#pragma once
// SOURCE: Hull (2018) "Options, Futures, and Other Derivatives" Ch.7
// SOURCE: Brigo & Mercurio (2006) "Interest Rate Models - Theory and Practice" Ch.1
// 模块: 利率互换 (Interest Rate Swap, IRS)
//
// Vanilla Fixed-Floating Swap:
//   - 固定端: 在支付日 T_1, ..., T_n 支付 K * τ_i * N
//   - 浮动端: 基于参考利率 (LIBOR/term rate), 在同样支付日支付 τ_i * L(T_{i-1}, T_i) * N
//   - 在 t=0 进入 swap (spot start) 或 T_start > 0 进入 (forward start)
//
// 定价 (无套利, 假设浮动端基于零息债):
//   Float leg PV (t=0) = N * [P(0, T_start) - P(0, T_n)]
//     (floating rate note 在 reset 时刻等于 par, 即 PV=1 at reset)
//   Fixed leg PV (t=0) = N * K * Σ τ_i * P(0, T_i)
//
// Payer swap (pay fixed, receive float) PV = Float - Fixed
// Receiver swap (receive fixed, pay float) PV = Fixed - Float
//
// 互换利率 (par swap rate): 让 PV_payer = 0 的 K 值
//   K_par = [P(0, T_start) - P(0, T_n)] / Σ τ_i * P(0, T_i)

#include "cpphub/core/types.hpp"
#include "cpphub/models/ir/short_rate.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"  // ZeroCurve (dual-curve 支持)
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ IRS 配置 ============
struct IRSConfig {
    Real notional = 1.0;          // 名义本金 N
    Real fixed_rate;              // 固定利率 K (年化)
    Real start_time = 0.0;        // T_start (swap 起始/reset 日, 默认 spot start)
    std::vector<Real> payment_times;   // T_1, ..., T_n (固定+浮动端支付日, 严格递增)
    std::vector<Real> year_fractions;  // τ_1, ..., τ_n (年化分数, 对应 T_i - T_{i-1})
    bool is_payer = true;         // true: payer (付固定收浮动); false: receiver

    void validate() const {
        if (notional <= 0.0) throw std::invalid_argument("IRS: notional must be positive");
        if (fixed_rate < 0.0) throw std::invalid_argument("IRS: fixed_rate must be non-negative");
        if (payment_times.size() != year_fractions.size()) {
            throw std::invalid_argument("IRS: payment_times and year_fractions size mismatch");
        }
        if (payment_times.size() < 1) {
            throw std::invalid_argument("IRS: at least 1 payment required");
        }
        if (start_time < 0.0) throw std::invalid_argument("IRS: start_time must be non-negative");
        // 支付时间必须严格递增且 ≥ start_time
        for (Size i = 0; i < payment_times.size(); ++i) {
            if (payment_times[i] <= start_time) {
                throw std::invalid_argument("IRS: payment_times must be > start_time");
            }
            if (i > 0 && payment_times[i] <= payment_times[i - 1]) {
                throw std::invalid_argument("IRS: payment_times must be strictly increasing");
            }
            if (year_fractions[i] <= 0.0) {
                throw std::invalid_argument("IRS: year_fractions must be positive");
            }
        }
    }
};

// ============ IRS 定价器 ============
// 接受任意提供 zero_coupon_bond(T) 接口的利率模型 (Vasicek/CIR/HW/G2++)
class InterestRateSwap {
public:
    explicit InterestRateSwap(IRSConfig cfg) : cfg_(std::move(cfg)) {
        cfg_.validate();
    }

    // 固定端现值: N * K * Σ τ_i * P(0, T_i)
    template <typename Model>
    Real fixed_leg_pv(const Model& model) const {
        Real pv = 0.0;
        for (Size i = 0; i < cfg_.payment_times.size(); ++i) {
            Real P = model.zero_coupon_bond(cfg_.payment_times[i]);
            pv += cfg_.year_fractions[i] * P;
        }
        return cfg_.notional * cfg_.fixed_rate * pv;
    }

    // 浮动端现值: N * [P(0, T_start) - P(0, T_n)]
    // 假设浮动端 reset 在 T_start, 最后支付在 T_n
    // (floating rate note 在 reset 日等于 par 的无套利结果)
    template <typename Model>
    Real float_leg_pv(const Model& model) const {
        Real P_start = model.zero_coupon_bond(cfg_.start_time);
        Real P_end = model.zero_coupon_bond(cfg_.payment_times.back());
        return cfg_.notional * (P_start - P_end);
    }

    // 互换现值
    template <typename Model>
    Real present_value(const Model& model) const {
        Real fixed_pv = fixed_leg_pv(model);
        Real float_pv = float_leg_pv(model);
        return cfg_.is_payer ? (float_pv - fixed_pv) : (fixed_pv - float_pv);
    }

    // 互换利率 (par swap rate): 让 payer PV = 0 的 K
    // K_par = [P(0, T_start) - P(0, T_n)] / Σ τ_i * P(0, T_i)
    template <typename Model>
    Real par_swap_rate(const Model& model) const {
        Real annuity = 0.0;  // Σ τ_i * P(0, T_i)
        for (Size i = 0; i < cfg_.payment_times.size(); ++i) {
            annuity += cfg_.year_fractions[i] * model.zero_coupon_bond(cfg_.payment_times[i]);
        }
        if (annuity <= 0.0) {
            throw std::runtime_error("IRS: annuity must be positive");
        }
        Real P_start = model.zero_coupon_bond(cfg_.start_time);
        Real P_end = model.zero_coupon_bond(cfg_.payment_times.back());
        return (P_start - P_end) / annuity;
    }

    // 年金因子 (annuity factor): Σ τ_i * P(0, T_i)
    template <typename Model>
    Real annuity(const Model& model) const {
        Real a = 0.0;
        for (Size i = 0; i < cfg_.payment_times.size(); ++i) {
            a += cfg_.year_fractions[i] * model.zero_coupon_bond(cfg_.payment_times[i]);
        }
        return a;
    }

    const IRSConfig& config() const { return cfg_; }
    Real notional() const { return cfg_.notional; }
    Real fixed_rate() const { return cfg_.fixed_rate; }
    bool is_payer() const { return cfg_.is_payer; }

private:
    IRSConfig cfg_;
};

// ============ 便捷工厂: 等间隔 IRS ============
// 构造 T_start 开始, 每 τ 支付一次, 共 n_payments 期的 vanilla IRS
inline IRSConfig make_vanilla_irs(Real notional, Real fixed_rate, Real start_time,
                                    Real tau, Size n_payments, bool is_payer = true) {
    if (tau <= 0.0) throw std::invalid_argument("make_vanilla_irs: tau must be positive");
    if (n_payments == 0) throw std::invalid_argument("make_vanilla_irs: n_payments must be > 0");
    IRSConfig cfg;
    cfg.notional = notional;
    cfg.fixed_rate = fixed_rate;
    cfg.start_time = start_time;
    cfg.is_payer = is_payer;
    cfg.payment_times.resize(n_payments);
    cfg.year_fractions.resize(n_payments, tau);
    for (Size i = 0; i < n_payments; ++i) {
        cfg.payment_times[i] = start_time + static_cast<Real>(i + 1) * tau;
    }
    return cfg;
}

// ============ Multi-Curve IRS (OIS Discounting) ============
// Dual-curve: LIBOR projection curve (P_f) + OIS discount curve (P_d)
// Post-crisis standard: LIBOR 用于 forward rate projection, OIS 用于折现.
//
// Fixed leg PV = N * K * Σ τ_i * P_d(0, T_i)
// Float leg PV = N * Σ τ_i * F_i * P_d(0, T_i)
//   where F_i = (P_f(0, T_{i-1})/P_f(0, T_i) - 1) / τ_i  (LIBOR forward)
//
// 单曲线退化 (projection = discount): Float PV = N * [P(0, T_start) - P(0, T_n)] (FRN par)

// 固定端 PV (OIS discounting)
inline Real irs_dual_curve_fixed_leg_pv(const IRSConfig& cfg, const ZeroCurve& discount_curve) {
    Real pv = 0.0;
    for (Size i = 0; i < cfg.payment_times.size(); ++i) {
        pv += cfg.year_fractions[i] * discount_curve.discount_factor(cfg.payment_times[i]);
    }
    return cfg.notional * cfg.fixed_rate * pv;
}

// 浮动端 PV (multi-curve: LIBOR projection + OIS discount)
inline Real irs_dual_curve_float_leg_pv(const IRSConfig& cfg,
                                         const ZeroCurve& projection_curve,
                                         const ZeroCurve& discount_curve) {
    Real pv = 0.0;
    for (Size i = 0; i < cfg.payment_times.size(); ++i) {
        Real T_prev = (i == 0) ? cfg.start_time : cfg.payment_times[i - 1];
        Real T_curr = cfg.payment_times[i];
        Real tau = cfg.year_fractions[i];
        Real P_f_prev = projection_curve.discount_factor(T_prev);
        Real P_f_curr = projection_curve.discount_factor(T_curr);
        Real F = (P_f_prev / P_f_curr - 1.0) / tau;  // LIBOR forward rate
        Real P_d = discount_curve.discount_factor(T_curr);
        pv += tau * F * P_d;
    }
    return cfg.notional * pv;
}

// IRS PV (multi-curve): payer = Float - Fixed, receiver = Fixed - Float
inline Real irs_dual_curve_pv(const IRSConfig& cfg,
                               const ZeroCurve& projection_curve,
                               const ZeroCurve& discount_curve) {
    Real fixed = irs_dual_curve_fixed_leg_pv(cfg, discount_curve);
    Real flt = irs_dual_curve_float_leg_pv(cfg, projection_curve, discount_curve);
    return cfg.is_payer ? (flt - fixed) : (fixed - flt);
}

// Annuity (OIS discounting): Σ τ_i * P_d(0, T_i)
inline Real irs_dual_curve_annuity(const IRSConfig& cfg, const ZeroCurve& discount_curve) {
    Real a = 0.0;
    for (Size i = 0; i < cfg.payment_times.size(); ++i) {
        a += cfg.year_fractions[i] * discount_curve.discount_factor(cfg.payment_times[i]);
    }
    return a;
}

// Par swap rate (dual-curve): K_par = Float_PV_per_N / annuity
inline Real irs_dual_curve_par_rate(const IRSConfig& cfg,
                                     const ZeroCurve& projection_curve,
                                     const ZeroCurve& discount_curve) {
    // Float PV per unit notional
    Real float_pv_per_N = 0.0;
    for (Size i = 0; i < cfg.payment_times.size(); ++i) {
        Real T_prev = (i == 0) ? cfg.start_time : cfg.payment_times[i - 1];
        Real T_curr = cfg.payment_times[i];
        Real tau = cfg.year_fractions[i];
        Real P_f_prev = projection_curve.discount_factor(T_prev);
        Real P_f_curr = projection_curve.discount_factor(T_curr);
        Real F = (P_f_prev / P_f_curr - 1.0) / tau;
        Real P_d = discount_curve.discount_factor(T_curr);
        float_pv_per_N += tau * F * P_d;
    }
    Real annuity = irs_dual_curve_annuity(cfg, discount_curve);
    if (annuity <= 0.0) {
        throw std::runtime_error("irs_dual_curve_par_rate: annuity must be positive");
    }
    return float_pv_per_N / annuity;
}

}  // namespace v1
}  // namespace cpphub
