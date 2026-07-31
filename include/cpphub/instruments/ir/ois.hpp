#pragma once
// SOURCE: Hull (2018) Ch.7, Henrard (2014) "Multi-Curve Framework"
// SOURCE: Clarke (2010) "OIS Discounting" (Barclays Capital)
// 模块: 隔夜指数互换 (Overnight Index Swap, OIS)
//
// OIS 结构:
//   - 固定端: 在 T_1,...,T_n 支付 K * τ_i * N (K = OIS swap rate, 简单复利)
//   - 浮动端: 基于隔夜指数 (SOFR/ESTR/SONIA) 复合平均, 在同样日期支付 τ_i * L_OIS(T_{i-1},T_i) * N
//   - 单曲线框架 (post-crisis 标准): OIS curve 既作 projection (浮端远期) 又作 discount
//
// 定价 (单曲线, OIS curve = discount = projection):
//   Fixed leg PV = N * K * Σ τ_i * P(0, T_i)
//   Float leg PV = N * [P(0, T_start) - P(0, T_n)]  (FRN par 性质)
//   Payer OIS (pay fixed, receive float) PV = Float - Fixed
//
// Par OIS rate: 让 PV_payer = 0
//   K_par = [P(0, T_start) - P(0, T_n)] / Σ τ_i * P(0, T_i)
//
// 与 IRS 的关系: 数学结构相同, 但 OIS 基于 OIS curve, IRS (传统) 基于 LIBOR curve.
// 多曲线框架下 IRS 浮端用 LIBOR projection, 贴现用 OIS curve, 两者分离.

#include "cpphub/core/types.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ OIS 配置 ============
struct OISConfig {
    Real notional = 1.0;          // 名义本金 N
    Real fixed_rate;              // 固定利率 K (年化, 简单复利)
    Real start_time = 0.0;        // T_start (swap 起始, 默认 spot start)
    std::vector<Real> payment_times;   // T_1, ..., T_n (固定+浮动端支付日, 严格递增)
    std::vector<Real> year_fractions;  // τ_1, ..., τ_n (年化分数)
    bool is_payer = true;         // true: payer (付固定收浮动); false: receiver

    void validate() const {
        if (notional <= 0.0) throw std::invalid_argument("OIS: notional must be positive");
        if (fixed_rate < 0.0) throw std::invalid_argument("OIS: fixed_rate must be non-negative");
        if (payment_times.size() != year_fractions.size()) {
            throw std::invalid_argument("OIS: payment_times and year_fractions size mismatch");
        }
        if (payment_times.size() < 1) {
            throw std::invalid_argument("OIS: at least 1 payment required");
        }
        if (start_time < 0.0) throw std::invalid_argument("OIS: start_time must be non-negative");
        for (Size i = 0; i < payment_times.size(); ++i) {
            if (payment_times[i] <= start_time) {
                throw std::invalid_argument("OIS: payment_times must be > start_time");
            }
            if (i > 0 && payment_times[i] <= payment_times[i - 1]) {
                throw std::invalid_argument("OIS: payment_times must be strictly increasing");
            }
            if (year_fractions[i] <= 0.0) {
                throw std::invalid_argument("OIS: year_fractions must be positive");
            }
        }
    }
};

// ============ OIS 定价器 (单曲线框架) ============
class OvernightIndexSwap {
public:
    explicit OvernightIndexSwap(OISConfig cfg) : cfg_(std::move(cfg)) {
        cfg_.validate();
    }

    // 固定端现值: N * K * Σ τ_i * P(0, T_i)
    Real fixed_leg_pv(const ZeroCurve& curve) const {
        Real pv = 0.0;
        for (Size i = 0; i < cfg_.payment_times.size(); ++i) {
            Real P = curve.discount_factor(cfg_.payment_times[i]);
            pv += cfg_.year_fractions[i] * P;
        }
        return cfg_.notional * cfg_.fixed_rate * pv;
    }

    // 浮动端现值 (单曲线, FRN par 性质): N * [P(0, T_start) - P(0, T_n)]
    // 若 start_time == 0: P(0, 0) = 1, 浮端 = N * [1 - P(0, T_n)]
    Real float_leg_pv(const ZeroCurve& curve) const {
        Real P_start = curve.discount_factor(cfg_.start_time);
        Real P_end = curve.discount_factor(cfg_.payment_times.back());
        return cfg_.notional * (P_start - P_end);
    }

    // OIS PV: payer = Float - Fixed, receiver = Fixed - Float
    Real pv(const ZeroCurve& curve) const {
        Real fixed = fixed_leg_pv(curve);
        Real flt = float_leg_pv(curve);
        return cfg_.is_payer ? (flt - fixed) : (fixed - flt);
    }

    // Par OIS rate: 让 PV_payer = 0
    //   K_par = [P(0, T_start) - P(0, T_n)] / Σ τ_i * P(0, T_i)
    Real par_rate(const ZeroCurve& curve) const {
        Real P_start = curve.discount_factor(cfg_.start_time);
        Real P_end = curve.discount_factor(cfg_.payment_times.back());
        Real annuity = 0.0;
        for (Size i = 0; i < cfg_.payment_times.size(); ++i) {
            annuity += cfg_.year_fractions[i] * curve.discount_factor(cfg_.payment_times[i]);
        }
        if (annuity <= 0.0) {
            throw std::runtime_error("OIS::par_rate: annuity must be positive");
        }
        return (P_start - P_end) / annuity;
    }

    // Annuity (固定端单位利率的 PV): Σ τ_i * P(0, T_i)
    Real annuity(const ZeroCurve& curve) const {
        Real a = 0.0;
        for (Size i = 0; i < cfg_.payment_times.size(); ++i) {
            a += cfg_.year_fractions[i] * curve.discount_factor(cfg_.payment_times[i]);
        }
        return a;
    }

    // Forward OIS rate (从 T_start 到 T_n 的远期 OIS rate, 单期近似)
    // F = (P(0,T_start)/P(0,T_n) - 1) / Σ τ_i  (仅对单期 OIS 精确, 多期为平均)
    Real forward_rate(const ZeroCurve& curve) const {
        Real P_start = curve.discount_factor(cfg_.start_time);
        Real P_end = curve.discount_factor(cfg_.payment_times.back());
        Real total_tau = 0.0;
        for (Size i = 0; i < cfg_.year_fractions.size(); ++i) {
            total_tau += cfg_.year_fractions[i];
        }
        if (total_tau <= 0.0) return 0.0;
        return (P_start / P_end - 1.0) / total_tau;
    }

    const OISConfig& config() const noexcept { return cfg_; }
    Real notional() const noexcept { return cfg_.notional; }
    Real fixed_rate() const noexcept { return cfg_.fixed_rate; }
    bool is_payer() const noexcept { return cfg_.is_payer; }

private:
    OISConfig cfg_;
};

// ============ 工厂函数: 构造单期 OIS (短期 OIS, ≤1Y) ============
inline OISConfig make_single_period_ois(Real notional, Real fixed_rate, Real T_start, Real T_end,
                                         Real year_fraction, bool is_payer = true) {
    OISConfig cfg;
    cfg.notional = notional;
    cfg.fixed_rate = fixed_rate;
    cfg.start_time = T_start;
    cfg.payment_times = {T_end};
    cfg.year_fractions = {year_fraction > 0 ? year_fraction : (T_end - T_start)};
    cfg.is_payer = is_payer;
    return cfg;
}

// ============ 工厂函数: 构造等距支付 OIS (annual coupon) ============
inline OISConfig make_annual_ois(Real notional, Real fixed_rate, Real T_start, Real T_end,
                                  Size n_payments, bool is_payer = true) {
    OISConfig cfg;
    cfg.notional = notional;
    cfg.fixed_rate = fixed_rate;
    cfg.start_time = T_start;
    cfg.payment_times.resize(n_payments);
    cfg.year_fractions.resize(n_payments);
    Real dt = (T_end - T_start) / static_cast<Real>(n_payments);
    for (Size i = 0; i < n_payments; ++i) {
        cfg.payment_times[i] = T_start + (i + 1) * dt;
        cfg.year_fractions[i] = dt;
    }
    cfg.is_payer = is_payer;
    return cfg;
}

}  // namespace v1
}  // namespace cpphub
