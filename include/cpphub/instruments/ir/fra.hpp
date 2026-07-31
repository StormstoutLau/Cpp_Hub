#pragma once
// SOURCE: Hull (2018) "Options, Futures, and Other Derivatives" Ch.4 (FRA)
// SOURCE: Henrard (2014) "Interest Rate Modelling in the Multi-Curve Framework" Ch.3
// 模块: 远期利率协议 (Forward Rate Agreement, FRA) - Multi-Curve OIS Discounting
//
// FRA 结构:
//   - 在 T_reset (reset/fixing date) 重置参考利率 L(T_reset, T_pay)
//   - 在 T_settle (settlement date, 通常 = T_reset) 时刻提前结算
//   - Payoff at T_settle: N * τ * (L(T_reset, T_pay) - K) / (1 + τ * L)   (market standard)
//     其中 τ = T_pay - T_reset (年化)
//   - 持有 FRA (long, pay fixed K receive floating L) 在 L > K 时获利
//
// Multi-Curve 定价 (post-crisis standard):
//   Forward rate F (LIBOR projection):
//     F = (P_f(0, T_reset) / P_f(0, T_pay) - 1) / τ   (从 projection curve 推导)
//   Payoff PV (折现到 t=0):
//     PV = N * τ * (F - K) / (1 + τ * F) * P_d(0, T_settle)
//     其中 P_d 是 discount curve (OIS), T_settle 通常 = T_reset
//
// 单曲线退化 (传统 pre-crisis, projection = discount):
//   F = (P(0,T_reset)/P(0,T_pay) - 1) / τ
//   PV = N * τ * (F - K) / (1 + τ * F) * P(0, T_reset)
//
// FRA rate (par forward rate, 让 PV=0): K_par = F

#include "cpphub/core/types.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include <cmath>
#include <stdexcept>

namespace cpphub {
inline namespace v1 {

// ============ FRA 配置 ============
struct FRAConfig {
    Real notional = 1.0;        // 名义本金 N
    Real fixed_rate;            // 固定利率 K (FRA rate, 年化简单复利)
    Real T_reset;               // 重置日 (fixing date)
    Real T_pay;                 // 支付日 (maturity of underlying LIBOR)
    Real T_settle = -1.0;       // 结算日 (默认 = T_reset, 提前结算)
    Real year_fraction = -1.0;  // τ = T_pay - T_reset (默认自动计算)
    bool is_long = true;        // true: long (付 K 收 L); false: short (收 K 付 L)

    void validate() const {
        if (notional <= 0.0) throw std::invalid_argument("FRA: notional must be positive");
        if (fixed_rate < 0.0) throw std::invalid_argument("FRA: fixed_rate must be non-negative");
        if (T_reset < 0.0) throw std::invalid_argument("FRA: T_reset must be non-negative");
        if (T_pay <= T_reset) throw std::invalid_argument("FRA: T_pay must be > T_reset");
        Real Ts = (T_settle < 0.0) ? T_reset : T_settle;
        if (Ts < 0.0) throw std::invalid_argument("FRA: T_settle must be non-negative");
        if (Ts > T_pay) throw std::invalid_argument("FRA: T_settle must be <= T_pay");
    }

    Real settle_time() const { return (T_settle < 0.0) ? T_reset : T_settle; }
    Real tau() const { return (year_fraction > 0.0) ? year_fraction : (T_pay - T_reset); }
};

// ============ FRA 定价器 (Multi-Curve) ============
class ForwardRateAgreement {
public:
    explicit ForwardRateAgreement(FRAConfig cfg) : cfg_(std::move(cfg)) {
        cfg_.validate();
    }

    // Forward LIBOR rate from projection curve:
    //   F = (P_f(0, T_reset) / P_f(0, T_pay) - 1) / τ
    Real forward_rate(const ZeroCurve& projection_curve) const {
        Real Pr = projection_curve.discount_factor(cfg_.T_reset);
        Real Pp = projection_curve.discount_factor(cfg_.T_pay);
        return (Pr / Pp - 1.0) / cfg_.tau();
    }

    // FRA PV (multi-curve OIS discounting, market standard settlement)
    //   PV = ±N * τ * (F - K) / (1 + τ * F) * P_d(0, T_settle)
    //   +: long (pay K, receive L); -: short
    Real pv(const ZeroCurve& projection_curve, const ZeroCurve& discount_curve) const {
        Real F = forward_rate(projection_curve);
        Real Pd = discount_curve.discount_factor(cfg_.settle_time());
        Real payoff = cfg_.notional * cfg_.tau() * (F - cfg_.fixed_rate) / (1.0 + cfg_.tau() * F);
        return cfg_.is_long ? payoff * Pd : -payoff * Pd;
    }

    // 单曲线便捷接口 (projection = discount)
    Real pv(const ZeroCurve& curve) const {
        return pv(curve, curve);
    }

    // Par FRA rate (让 PV = 0): K_par = F
    Real par_rate(const ZeroCurve& projection_curve) const {
        return forward_rate(projection_curve);
    }

    // 持有方 P&L at T_settle (未折现): N * τ * (L - K) / (1 + τ * L)
    Real payoff_at_settle(Real realized_libor) const {
        Real p = cfg_.notional * cfg_.tau() * (realized_libor - cfg_.fixed_rate)
                 / (1.0 + cfg_.tau() * realized_libor);
        return cfg_.is_long ? p : -p;
    }

    // 简单 PV (不提前折现, payoff at T_pay)
    //   PV = ±N * τ * (F - K) * P_d(0, T_pay)
    // 用于多期 IRS 中 caplet/floorlet 的等价计算
    Real pv_simple(const ZeroCurve& projection_curve, const ZeroCurve& discount_curve) const {
        Real F = forward_rate(projection_curve);
        Real Pd = discount_curve.discount_factor(cfg_.T_pay);
        Real payoff = cfg_.notional * cfg_.tau() * (F - cfg_.fixed_rate);
        return cfg_.is_long ? payoff * Pd : -payoff * Pd;
    }

    const FRAConfig& config() const noexcept { return cfg_; }
    Real notional() const noexcept { return cfg_.notional; }
    Real fixed_rate() const noexcept { return cfg_.fixed_rate; }

private:
    FRAConfig cfg_;
};

// ============ 工厂函数: 标准 3M FRA (T_reset, T_pay = T_reset + 0.25) ============
inline FRAConfig make_3m_fra(Real notional, Real fixed_rate, Real T_reset,
                              bool is_long = true) {
    FRAConfig cfg;
    cfg.notional = notional;
    cfg.fixed_rate = fixed_rate;
    cfg.T_reset = T_reset;
    cfg.T_pay = T_reset + 0.25;
    cfg.T_settle = T_reset;  // standard settlement at fixing date
    cfg.year_fraction = 0.25;
    cfg.is_long = is_long;
    return cfg;
}

// ============ 工厂函数: 6M FRA ============
inline FRAConfig make_6m_fra(Real notional, Real fixed_rate, Real T_reset,
                              bool is_long = true) {
    FRAConfig cfg;
    cfg.notional = notional;
    cfg.fixed_rate = fixed_rate;
    cfg.T_reset = T_reset;
    cfg.T_pay = T_reset + 0.5;
    cfg.T_settle = T_reset;
    cfg.year_fraction = 0.5;
    cfg.is_long = is_long;
    return cfg;
}

}  // namespace v1
}  // namespace cpphub
