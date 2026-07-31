#pragma once
// SOURCE: Henrard (2014) "Interest Rate Modelling in the Multi-Curve Framework" Ch.5-6
// SOURCE: Mercurio (2010) "Interest Rates and FVA Calculations" (multi-curve bootstrapping)
// 模块: 基差互换 (Basis Swap) - Multi-Curve 框架
//
// Basis Swap 结构:
//   两条浮动腿交换, 不同 tenor (e.g. 3M LIBOR vs 6M LIBOR) 或 不同基准 (LIBOR vs OIS)
//   Leg 1: N * τ_i * L^{(1)}(T_{i-1}, T_i) at payment dates T_i  (e.g. 3M LIBOR)
//   Leg 2: N * τ_j * (L^{(2)}(T_{j-1}, T_j) + spread) at T_j       (e.g. 6M LIBOR + spread)
//   spread = tenor basis spread (通常 6M > 3M, 信用/流动性风险补偿)
//
// Multi-Curve 定价 (post-crisis standard):
//   Forward rate on leg k from projection curve k:
//     F_k(T_{i-1}, T_i) = (P_f^{(k)}(0, T_{i-1}) / P_f^{(k)}(0, T_i) - 1) / τ_i
//   Leg k PV (浮端 + spread):
//     PV_k = N * Σ_i τ_i * (F_k + spread_k) * P_d(0, T_i)
//          = N * Σ_i (P_f^{(k)}(0,T_{i-1})/P_f^{(k)}(0,T_i) - 1) * P_d(0, T_i)
//            + N * spread_k * Σ_i τ_i * P_d(0, T_i)
//   其中 P_d = discount curve (OIS), P_f^{(k)} = projection curve k (LIBOR tenor k)
//
// 单曲线退化 (projection = discount):
//   Leg PV (no spread) = N * [P(0, T_start) - P(0, T_n)]  (FRN par 性质)
//
// Par basis spread: 让 PV_leg1 = PV_leg2 解出 leg2 上的 spread
//   若 leg1 上 spread=0, leg2 上 spread=s, 则 s_par = (PV_leg1 - PV_leg2_no_spread) / annuity_leg2

#include "cpphub/core/types.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ 浮动腿配置 ============
struct FloatLegConfig {
    Real notional = 1.0;
    Real start_time = 0.0;
    std::vector<Real> payment_times;   // T_1, ..., T_n (严格递增)
    std::vector<Real> year_fractions;  // τ_1, ..., τ_n
    Real spread = 0.0;                 // 浮端利差 (加到 forward rate)

    void validate() const {
        if (notional <= 0.0) throw std::invalid_argument("FloatLeg: notional must be positive");
        if (payment_times.size() != year_fractions.size()) {
            throw std::invalid_argument("FloatLeg: payment_times/year_fractions size mismatch");
        }
        if (payment_times.empty()) throw std::invalid_argument("FloatLeg: at least 1 payment required");
        if (start_time < 0.0) throw std::invalid_argument("FloatLeg: start_time must be non-negative");
        for (Size i = 0; i < payment_times.size(); ++i) {
            if (payment_times[i] <= start_time) {
                throw std::invalid_argument("FloatLeg: payment_times must be > start_time");
            }
            if (i > 0 && payment_times[i] <= payment_times[i - 1]) {
                throw std::invalid_argument("FloatLeg: payment_times must be strictly increasing");
            }
            if (year_fractions[i] <= 0.0) {
                throw std::invalid_argument("FloatLeg: year_fractions must be positive");
            }
        }
    }

    // reset times: T_0=start_time, T_i=payment_times[i-1] for i>=1
    // 即第 i 期 [T_{i-1}, T_i], reset 在 T_{i-1}
    Real reset_time(Size i) const {
        if (i == 0) return start_time;
        return payment_times[i - 1];
    }
};

// ============ 浮动腿 PV (Multi-Curve) ============
// PV = N * Σ_i τ_i * (F_i + spread) * P_d(0, T_i)
//   其中 F_i = (P_f(0, T_{i-1})/P_f(0, T_i) - 1) / τ_i
inline Real float_leg_pv(const FloatLegConfig& leg,
                         const ZeroCurve& projection_curve,
                         const ZeroCurve& discount_curve) {
    Real pv = 0.0;
    for (Size i = 0; i < leg.payment_times.size(); ++i) {
        Real T_prev = leg.reset_time(i);
        Real T_curr = leg.payment_times[i];
        Real tau = leg.year_fractions[i];
        Real P_f_prev = projection_curve.discount_factor(T_prev);
        Real P_f_curr = projection_curve.discount_factor(T_curr);
        Real F = (P_f_prev / P_f_curr - 1.0) / tau;  // forward LIBOR
        Real P_d = discount_curve.discount_factor(T_curr);
        pv += tau * (F + leg.spread) * P_d;
    }
    return leg.notional * pv;
}

// 浮动腿 annuity (单位 spread 的 PV): Σ τ_i * P_d(0, T_i)
inline Real float_leg_annuity(const FloatLegConfig& leg, const ZeroCurve& discount_curve) {
    Real a = 0.0;
    for (Size i = 0; i < leg.payment_times.size(); ++i) {
        a += leg.year_fractions[i] * discount_curve.discount_factor(leg.payment_times[i]);
    }
    return leg.notional * a;
}

// ============ Basis Swap (两条浮动腿交换) ============
class BasisSwap {
public:
    BasisSwap(FloatLegConfig leg1, FloatLegConfig leg2, bool pay_leg1 = true)
        : leg1_(std::move(leg1)), leg2_(std::move(leg2)), pay_leg1_(pay_leg1) {
        leg1_.validate();
        leg2_.validate();
        if (std::abs(leg1_.notional - leg2_.notional) > 1e-10) {
            throw std::invalid_argument("BasisSwap: leg1/leg2 notional must match");
        }
    }

    // Leg 1 PV (基于 projection_curve_1)
    Real leg1_pv(const ZeroCurve& proj1, const ZeroCurve& discount) const {
        return float_leg_pv(leg1_, proj1, discount);
    }

    // Leg 2 PV (基于 projection_curve_2)
    Real leg2_pv(const ZeroCurve& proj2, const ZeroCurve& discount) const {
        return float_leg_pv(leg2_, proj2, discount);
    }

    // Basis Swap PV: pay_leg1 ? PV_leg1 - PV_leg2 : PV_leg2 - PV_leg1
    // 注意: leg2 的 spread 已包含在 leg2_pv 中
    Real pv(const ZeroCurve& proj1, const ZeroCurve& proj2,
            const ZeroCurve& discount) const {
        Real pv1 = leg1_pv(proj1, discount);
        Real pv2 = leg2_pv(proj2, discount);
        return pay_leg1_ ? (pv2 - pv1) : (pv1 - pv2);
    }

    // Par basis spread (让 PV = 0 时 leg2 上的 spread)
    // 解 spread: PV_leg2(s) = PV_leg1
    // PV_leg2(s) = PV_leg2(0) + s * annuity_leg2
    // s_par = (PV_leg1 - PV_leg2(0)) / annuity_leg2
    Real par_spread_leg2(const ZeroCurve& proj1, const ZeroCurve& proj2,
                         const ZeroCurve& discount) const {
        // 临时将 leg2.spread 设为 0
        FloatLegConfig leg2_no_spread = leg2_;
        leg2_no_spread.spread = 0.0;
        Real pv1 = float_leg_pv(leg1_, proj1, discount);
        Real pv2_0 = float_leg_pv(leg2_no_spread, proj2, discount);
        Real annuity2 = float_leg_annuity(leg2_, discount);
        if (annuity2 <= 0.0) {
            throw std::runtime_error("BasisSwap::par_spread: annuity must be positive");
        }
        // PV(pay leg1) = PV_leg2(s) - PV_leg1 = 0
        // PV_leg2(0) + s * annuity2 = PV_leg1
        return (pv1 - pv2_0) / annuity2;
    }

    // 等价的 par spread on leg1 (让 PV=0)
    Real par_spread_leg1(const ZeroCurve& proj1, const ZeroCurve& proj2,
                         const ZeroCurve& discount) const {
        FloatLegConfig leg1_no_spread = leg1_;
        leg1_no_spread.spread = 0.0;
        Real pv1_0 = float_leg_pv(leg1_no_spread, proj1, discount);
        Real pv2 = float_leg_pv(leg2_, proj2, discount);
        Real annuity1 = float_leg_annuity(leg1_, discount);
        if (annuity1 <= 0.0) {
            throw std::runtime_error("BasisSwap::par_spread: annuity must be positive");
        }
        return (pv2 - pv1_0) / annuity1;
    }

    const FloatLegConfig& leg1() const noexcept { return leg1_; }
    const FloatLegConfig& leg2() const noexcept { return leg2_; }
    bool pay_leg1() const noexcept { return pay_leg1_; }

private:
    FloatLegConfig leg1_;
    FloatLegConfig leg2_;
    bool pay_leg1_;
};

// ============ 工厂函数: 等额支付浮动腿 ============
inline FloatLegConfig make_float_leg(Real notional, Real start_time, Real end_time,
                                     Size n_payments, Real spread = 0.0) {
    FloatLegConfig leg;
    leg.notional = notional;
    leg.start_time = start_time;
    leg.payment_times.resize(n_payments);
    leg.year_fractions.resize(n_payments);
    Real dt = (end_time - start_time) / static_cast<Real>(n_payments);
    for (Size i = 0; i < n_payments; ++i) {
        leg.payment_times[i] = start_time + (i + 1) * dt;
        leg.year_fractions[i] = dt;
    }
    leg.spread = spread;
    return leg;
}

// ============ 工厂函数: 3M vs 6M LIBOR basis swap ============
// leg1: 3M LIBOR (4 期/年), leg2: 6M LIBOR (2 期/年)
inline BasisSwap make_3m_6m_basis_swap(Real notional, Real T_start, Real T_end,
                                       Real spread_on_6m = 0.0, bool pay_3m = true) {
    Real total_tau = T_end - T_start;
    Size n_3m = static_cast<Size>(std::round(total_tau / 0.25));
    Size n_6m = static_cast<Size>(std::round(total_tau / 0.5));
    FloatLegConfig leg_3m = make_float_leg(notional, T_start, T_end, n_3m, 0.0);
    FloatLegConfig leg_6m = make_float_leg(notional, T_start, T_end, n_6m, spread_on_6m);
    return BasisSwap(leg_3m, leg_6m, pay_3m);
}

}  // namespace v1
}  // namespace cpphub
