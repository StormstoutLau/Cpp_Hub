#pragma once
// SOURCE: Jamshidian (1989) "An Exact Bond Option Formula"
// SOURCE: Brigo & Mercurio (2006) "Interest Rate Models - Theory and Practice" Ch.3
// 模块: 互换期权 (Swaption) — Hull-White 模型下的 Jamshidian 分解
//
// Swaption = 进入利率互换的期权
//   - Payer swaption: 有权在 T_ex 进入 payer swap (付固定, 收浮动)
//   - Receiver swaption: 有权在 T_ex 进入 receiver swap (收固定, 付浮动)
//
// 假设: swap 在 T_ex 起始 (spot start after exercise), 浮动端在 T_ex reset
//
// Jamshidian (1989) 分解 (单因子仿射模型, 如 HW/Vasicek):
//   1. T_ex 时刻 payer swap 价值 (per unit notional):
//      PV(T_ex, r) = 1 - P(T_ex, T_n, r) - K * Σ τ_i * P(T_ex, T_i, r)
//                   = c_0 + Σ_i c_i * P(T_ex, T_i, r)
//      c_0 = +1 (T_ex, 浮动端 par, P=1 常数)
//      c_i = -K * τ_i (固定端支付, i=1..n)
//      c_n += -1 (浮动端末端, 即若 T_n 同时是固定+浮动支付日, c_n = -K*τ_n - 1)
//   2. 二分搜索 r* 使 PV(r*) = 0:
//      P(T_ex, T_n, r*) + K * Σ τ_i * P(T_ex, T_i, r*) = 1
//   3. 对每个 i ≥ 1 (跳过 T_ex, 因 P=1 是常数):
//      P_i* = P(T_ex, T_i, r*) (Jamshidian strike)
//      分解为零息债期权:
//        if c_i > 0: call on P(T_ex, T_i), strike P_i*
//        if c_i < 0: put on P(T_ex, T_i), strike P_i*
//   4. Swaption_payer = N * Σ_i |c_i| * bond_option(T_ex, T_i, P_i*, is_call=(c_i>0))
//
// HW 模型解析公式 (Brigo-Mercurio eq 3.30):
//   P(t, T) = A(t, T) * exp(-B(t, T) * r(t))
//   B(t, T) = (1 - exp(-κ(T-t))) / κ
//   A(t, T) = P_M(0, T) / P_M(0, t) * exp(-B(t,T) * f(0, t) - 0.5 * σ² * B(t,T)² * V²(t))
//   V²(t) = (1 - exp(-2κ t)) / (2κ)
//
// 注: 仅支持 Hull-White 模型 (无套利校准, 市场标准).
//     Vasicek 的 Jamshidian 分解类似但 A(t,T) 公式不同, 需单独实现.

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/models/ir/short_rate.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ Swaption 配置 ============
struct SwaptionConfig {
    Real notional = 1.0;
    Real fixed_rate;                // K (年化, swap 固定利率)
    Real exercise_time;             // T_ex (期权到期, 同时是 swap 起始)
    std::vector<Real> payment_times;    // T_1, ..., T_n (swap 支付日, 全部 > T_ex)
    std::vector<Real> year_fractions;   // τ_1, ..., τ_n
    bool is_payer = true;           // true: payer swaption; false: receiver

    Size n_payments() const { return payment_times.size(); }

    void validate() const {
        if (notional <= 0.0) throw std::invalid_argument("Swaption: notional must be positive");
        if (fixed_rate < 0.0) throw std::invalid_argument("Swaption: fixed_rate must be non-negative");
        if (exercise_time <= 0.0) throw std::invalid_argument("Swaption: exercise_time must be positive");
        if (payment_times.size() != year_fractions.size()) {
            throw std::invalid_argument("Swaption: payment_times/year_fractions size mismatch");
        }
        if (payment_times.empty()) {
            throw std::invalid_argument("Swaption: at least 1 payment required");
        }
        for (Size i = 0; i < payment_times.size(); ++i) {
            if (payment_times[i] <= exercise_time) {
                throw std::invalid_argument("Swaption: payment_times must be > exercise_time");
            }
            if (year_fractions[i] <= 0.0) {
                throw std::invalid_argument("Swaption: year_fractions must be positive");
            }
            if (i > 0 && payment_times[i] <= payment_times[i - 1]) {
                throw std::invalid_argument("Swaption: payment_times must be strictly increasing");
            }
        }
    }
};

// ============ Swaption 定价器 (Hull-White Jamshidian 分解) ============
class Swaption {
public:
    Swaption(SwaptionConfig cfg, const HullWhite& model)
        : cfg_(std::move(cfg)), model_(model) {
        cfg_.validate();
    }

    // Jamshidian 分解定价 (HW 闭式)
    // 数学步骤:
    //   1. 构造 swap 现金流系数 c_i (per unit notional)
    //   2. 二分搜索 r* 使 Σ c_i * P(T_ex, T_i, r*) = -c_0 (即 PV=0)
    //   3. 对每个 i, 计算 P_i* = P(T_ex, T_i, r*)
    //   4. 分解为零息债期权: |c_i| * bond_option(T_ex, T_i, P_i*, is_call)
    Real present_value() const {
        const Size n = cfg_.n_payments();
        // 现金流系数 c_i (i=0..n), 其中 i=0 对应 T_ex (P=1 常数, 不参与期权分解)
        // Payer swap: c_0 = +1 (float par at T_ex), c_i = -K*τ_i (i=1..n-1), c_n = -K*τ_n - 1
        // Receiver swap: c_0 = -1, c_i = +K*τ_i (i=1..n-1), c_n = +K*τ_n + 1
        std::vector<Real> c(n + 1);
        Real sign = cfg_.is_payer ? -1.0 : +1.0;  // 固定端方向: payer 付 (负), receiver 收 (正)
        c[0] = cfg_.is_payer ? +1.0 : -1.0;       // 浮动端 par
        for (Size i = 1; i <= n; ++i) {
            c[i] = sign * cfg_.fixed_rate * cfg_.year_fractions[i - 1];
            if (i == n) {
                // 浮动末端: payer 付 1 (即 -1), receiver 收 1 (+1)
                c[i] += cfg_.is_payer ? -1.0 : +1.0;
            }
        }

        // Step 2: 二分搜索 r* 使 PV(r*) = c_0 + Σ_{i=1..n} c_i * P(T_ex, T_i, r*) = 0
        // 即 Σ_{i=1..n} c_i * P(T_ex, T_i, r*) = -c_0
        Real target = -c[0];  // = +1 for payer (浮端 par 抵消固定端+末端), -1 for receiver
        auto pv_at_r = [&](Real r) -> Real {
            Real sum = 0.0;
            for (Size i = 1; i <= n; ++i) {
                sum += c[i] * hw_bond_price(cfg_.exercise_time, cfg_.payment_times[i - 1], r);
            }
            return sum;
        };

        // 二分搜索 r* (典型 r* 在 [-0.2, 0.5] 范围内)
        Real r_lo = -0.5, r_hi = 1.0;
        Real pv_lo = pv_at_r(r_lo) - target;
        Real pv_hi = pv_at_r(r_hi) - target;
        if (pv_lo * pv_hi > 0.0) {
            // 无解或退化, 回退到 par swap rate 作为 r*
            // 回退方案: 直接用 PV 内在价值
            Real par_K = compute_par_swap_rate_at_exercise();
            Real intrinsic = cfg_.is_payer
                ? std::max(par_K - cfg_.fixed_rate, 0.0)
                : std::max(cfg_.fixed_rate - par_K, 0.0);
            // 用 annuity 近似
            Real annuity = 0.0;
            for (Size i = 0; i < n; ++i) {
                annuity += cfg_.year_fractions[i] * model_.zero_coupon_bond(cfg_.payment_times[i]);
            }
            return cfg_.notional * annuity * intrinsic;
        }
        // 二分法
        for (Size iter = 0; iter < 100; ++iter) {
            Real r_mid = 0.5 * (r_lo + r_hi);
            Real pv_mid = pv_at_r(r_mid) - target;
            if (std::abs(pv_mid) < 1e-12) {
                r_lo = r_hi = r_mid;
                break;
            }
            if (pv_lo * pv_mid < 0.0) {
                r_hi = r_mid;
                pv_hi = pv_mid;
            } else {
                r_lo = r_mid;
                pv_lo = pv_mid;
            }
        }
        Real r_star = 0.5 * (r_lo + r_hi);

        // Step 3: 对每个 i ≥ 1, 计算 P_i* 并分解为零息债期权
        Real pv = 0.0;
        for (Size i = 1; i <= n; ++i) {
            Real P_star = hw_bond_price(cfg_.exercise_time, cfg_.payment_times[i - 1], r_star);
            if (P_star <= 0.0) continue;  // 防御性
            // 现金流方向决定看涨/看跌:
            //   c_i > 0: 持有 P 的多头, payoff = c_i * max(P - P*, 0) → call on bond
            //   c_i < 0: 持有 P 的空头, payoff = |c_i| * max(P* - P, 0) → put on bond
            bool is_call = (c[i] > 0.0);
            Real weight = std::abs(c[i]);
            // HW.bond_option 返回 P(0, T_ex) * E[max(... , 0)]
            Real opt = model_.bond_option(cfg_.exercise_time, cfg_.payment_times[i - 1],
                                          P_star, is_call);
            pv += weight * opt;
        }
        return cfg_.notional * pv;
    }

    // 计算 T_ex 时刻的 par swap rate (用于估算内在价值)
    // K_par(T_ex) = [1 - P(T_ex, T_n)] / Σ τ_i * P(T_ex, T_i)
    // 注: P(T_ex, T_i) 由 HW 模型 E[P(T_ex, T_i)] 给出 (期望价格), 非某个 r* 下价格
    Real compute_par_swap_rate_at_exercise() const {
        const Size n = cfg_.n_payments();
        Real numerator = 1.0 - expected_bond_price(cfg_.exercise_time, cfg_.payment_times[n - 1]);
        Real annuity = 0.0;
        for (Size i = 0; i < n; ++i) {
            annuity += cfg_.year_fractions[i] * expected_bond_price(cfg_.exercise_time, cfg_.payment_times[i]);
        }
        if (annuity <= 0.0) return 0.0;
        return numerator / annuity;
    }

    const SwaptionConfig& config() const { return cfg_; }
    const HullWhite& model() const { return model_; }

private:
    SwaptionConfig cfg_;
    HullWhite model_;  // 拷贝模型 (避免引用生命周期问题)

    // HW 模型下, 给定 r(t)=r 时的 P(t, T) 解析表达 (Brigo-Mercurio eq 3.30)
    // P(t, T, r) = A(t, T) * exp(-B(t, T) * r)
    // B(t, T) = (1 - exp(-κ(T-t))) / κ
    // A(t, T) = P_M(0, T) / P_M(0, t) * exp(-B * f(0, t) - 0.5 * σ² * B² * V²(t))
    // V²(t) = (1 - exp(-2κ t)) / (2κ)
    Real hw_bond_price(Real t, Real T, Real r) const {
        if (T <= t) return 1.0;
        Real kappa = model_.params().kappa;
        Real sigma = model_.params().sigma;
        Real tau = T - t;
        Real B = affine_B(kappa, tau);
        // P_M(0, T), P_M(0, t) 来自市场曲线 (model 的 zero_coupon_bond)
        Real P_M_T = model_.zero_coupon_bond(T);
        Real P_M_t = model_.zero_coupon_bond(t);
        // f(0, t) = -∂ ln P(0, t) / ∂t (瞬时远期利率)
        Real f_0_t = model_.forward_rate(t);
        // V²(t) = σ² * (1 - exp(-2κ t)) / (2κ)  [Brigo-Mercurio (3.26), HW 条件方差]
        Real var_r = sigma * sigma * (1.0 - std::exp(-2.0 * kappa * t)) / (2.0 * kappa);
        // A(t, T) = (P_M(0,T) / P_M(0,t)) * exp(-B * f(0,t) - 0.5 * B² * var_r)
        Real A = (P_M_T / P_M_t)
                 * std::exp(-B * f_0_t - 0.5 * B * B * var_r);
        return A * std::exp(-B * r);
    }

    // HW 模型下 E[P(t, T)] (期望价格, 用于 par swap rate 估算)
    // E[P(t, T)] = A(t, T) * E[exp(-B * r(t))] = A(t, T) * exp(-B * E[r(t)] + 0.5 * B² * Var[r(t)])
    // E[r(t)] = r0 (假设 r0 = f(0,0) 校准后) 或 forward_rate(0)
    // Var[r(t)] = σ² * (1 - exp(-2κ t)) / (2κ)
    Real expected_bond_price(Real t, Real T) const {
        if (T <= t) return 1.0;
        Real kappa = model_.params().kappa;
        Real sigma = model_.params().sigma;
        Real tau = T - t;
        Real B = affine_B(kappa, tau);
        Real P_M_T = model_.zero_coupon_bond(T);
        Real P_M_t = model_.zero_coupon_bond(t);
        Real f_0_t = model_.forward_rate(t);
        Real var_r = sigma * sigma * (1.0 - std::exp(-2.0 * kappa * t)) / (2.0 * kappa);
        Real A = (P_M_T / P_M_t)
                 * std::exp(-B * f_0_t - 0.5 * B * B * var_r);
        // E[r(t)] 在 HW 下: E[r(t)] = r0*exp(-κt) + ∫₀^t exp(-κ(t-s)) θ(s) ds
        // 对无套利 HW, E[r(t)] = f(0, t) (远期利率)
        Real mean_r = model_.forward_rate(t);
        Real exp_factor = std::exp(-B * mean_r + 0.5 * B * B * var_r);
        return A * exp_factor;
    }
};

// ============ 便捷工厂: 等间隔 Swaption ============
// T_ex 行权, swap 在 T_ex 起始, 每 τ 支付一次, 共 n_payments 期
inline SwaptionConfig make_vanilla_swaption(Real notional, Real fixed_rate, Real exercise_time,
                                              Real tau, Size n_payments, bool is_payer = true) {
    if (tau <= 0.0) throw std::invalid_argument("make_vanilla_swaption: tau must be positive");
    if (n_payments == 0) throw std::invalid_argument("make_vanilla_swaption: n_payments must be > 0");
    SwaptionConfig cfg;
    cfg.notional = notional;
    cfg.fixed_rate = fixed_rate;
    cfg.exercise_time = exercise_time;
    cfg.is_payer = is_payer;
    cfg.payment_times.resize(n_payments);
    cfg.year_fractions.resize(n_payments, tau);
    for (Size i = 0; i < n_payments; ++i) {
        cfg.payment_times[i] = exercise_time + static_cast<Real>(i + 1) * tau;
    }
    return cfg;
}

}  // namespace v1
}  // namespace cpphub
