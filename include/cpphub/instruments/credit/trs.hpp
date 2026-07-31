#pragma once
// SOURCE: O'Kane (2008) "Modelling Single-Name and Multi-Name Credit Derivatives" Ch.10 (TRS)
// SOURCE: Schönbucher (2003) "Credit Derivatives Pricing Models" Ch.2, Ch.9
// SOURCE: Hull & White (2003) "The Valuation of Credit-Sensitive Notes"
// SOURCE: Das (2005) "Credit Derivatives: CDOs, Structured Products, and Synthetic Strategies" Ch.6
// 模块: 总回报互换 (Total Return Swap, TRS)
//
// ==================== TRS 结构 ====================
//
// TRS 是一种信用衍生品: 一方 (总回报支付方, payer) 将参考资产的全部回报 (价格变动 + 票息/股息)
// 转给另一方 (总回报接收方, receiver), 换取融资利率 (floating + spread) 的支付.
//
// 现金流 (从接收方视角, 正=收到):
//   Asset leg (总回报端, 收入):
//     - 价格回报: [S(t_i) - S(t_{i-1})] / S(t_{i-1}) * N  (可负, 价格下跌时支付方"收"负回报)
//     - 票息/股息: q * τ_i * N  (q 为 asset yield, 年化)
//   Funding leg (融资端, 支出):
//     - (L_i + spread) * τ_i * N  (L_i 为远期利率, spread 为融资利差)
//
// 违约处理 (参考实体违约发生在 (t_{i-1}, t_i] 时):
//   - 接收方支付 LGD = (1 - R) * N (信用损失)
//   - 合约在违约时终止, 假设违约在区间中点 t_mid
//   - 应计项 (asset return, funding) 按半期计算 (mid-point convention, 与 CDS 一致)
//
// ==================== 定价公式 (从接收方视角) ====================
//
// 确定性远期定价 (deterministic forward pricing):
//   远期价格 F(t) = S0 * exp((r - q) * t)  (风险中性测度下)
//   区间总回报: R_i = F(t_i)/F(t_{i-1}) - 1 + q * τ_i
//                    = exp((r-q)*τ_i) - 1 + q*τ_i
//
// Asset leg PV (receiver):
//   PV_asset = Σ_i {
//       Q(0, t_i) * R_i * N * P(0, t_i)                 // 生存: 完整回报
//     - PD(t_{i-1}, t_i) * (1-R) * N * P(0, t_mid_i)    // 违约: 信用损失
//   }
//
// Funding leg PV (receiver, 支出为正):
//   PV_funding = Σ_i {
//       Q(0, t_i) * (L_i + spread) * τ_i * N * P(0, t_i)             // 生存
//     + PD(t_{i-1}, t_i) * (L_i + spread) * τ_i/2 * N * P(0, t_mid_i) // 违约: 半期应计
//   }
//
// TRS PV (receiver) = PV_asset - PV_funding
//
// ==================== 约定 ====================
//
// - 支付频率: 季度 (τ ≈ 0.25) 或月度, 与 IRS 类似
// - 违约时点假设: 区间中点 (mid-point convention)
// - 折现: OIS discount curve
// - 信用曲线: CreditCurve (hazard rate + recovery)
// - 远期利率: L(T1, T2) = (P(T1)/P(T2) - 1) / (T2 - T1)  (简单复利)

#include "cpphub/core/types.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include "cpphub/instruments/credit/credit_curve.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>

namespace cpphub {
inline namespace v1 {

// ============ TRS 配置 ============
struct TRSConfig {
    Real notional = 1.0;                // 名义本金
    Real S0 = 100.0;                    // 初始资产价格
    Real asset_yield = 0.0;             // 资产票息/股息收益率 (年化, 如 0.03 = 3%)
    Real funding_spread = 0.005;        // 融资利差 (年化, 如 0.005 = 50bp)
    Real start_time = 0.0;              // 起始时间 (通常 0 = spot)
    Real maturity = 5.0;                // 到期时间 (年)
    Size n_payments = 20;               // 支付次数 (如 5Y 季度 = 20)
    bool is_receiver = true;            // true = 接收方视角, false = 支付方视角
    // 支付时间表 (若为空则按等间隔生成)
    std::vector<Real> payment_times;    // 若非空则用此, 否则按等间隔生成

    // 等间隔支付时间表生成
    void generate_schedule() {
        if (!payment_times.empty()) return;
        if (n_payments == 0) {
            throw std::invalid_argument("TRSConfig: n_payments must be positive");
        }
        if (maturity <= start_time) {
            throw std::invalid_argument("TRSConfig: maturity must be > start_time");
        }
        Real dt = (maturity - start_time) / static_cast<Real>(n_payments);
        payment_times.reserve(n_payments);
        for (Size i = 1; i <= n_payments; ++i) {
            payment_times.push_back(start_time + static_cast<Real>(i) * dt);
        }
    }

    // 验证配置
    // 注: funding_spread 允许为负 — 当 receiver 承担参考资产信用风险时,
    // par_spread 可能为负 (receiver 应获得信用补偿, 即支付负 spread).
    // 这与 "receiver takes credit loss" 约定一致 (见 CreditRiskReducesReceiverPV 测试).
    void validate() const {
        if (notional <= 0.0) throw std::invalid_argument("TRSConfig: notional must be positive");
        if (S0 <= 0.0) throw std::invalid_argument("TRSConfig: S0 must be positive");
        if (asset_yield < 0.0) throw std::invalid_argument("TRSConfig: asset_yield must be non-negative");
        if (maturity <= start_time) throw std::invalid_argument("TRSConfig: maturity must be > start_time");
        if (n_payments == 0) throw std::invalid_argument("TRSConfig: n_payments must be positive");
        if (!payment_times.empty()) {
            if (payment_times.size() != n_payments) {
                throw std::invalid_argument("TRSConfig: payment_times size != n_payments");
            }
            for (Size i = 0; i < payment_times.size(); ++i) {
                if (payment_times[i] <= 0.0) {
                    throw std::invalid_argument("TRSConfig: payment_times must be positive");
                }
                if (i > 0 && payment_times[i] <= payment_times[i - 1]) {
                    throw std::invalid_argument("TRSConfig: payment_times must be strictly increasing");
                }
            }
        }
    }
};

// ============ TRS 结果 ============
struct TRSResult {
    Real pv_asset_leg = 0.0;            // 总回报端 PV (receiver 视角, 正=收入)
    Real pv_funding_leg = 0.0;          // 融资端 PV (receiver 视角, 正=支出)
    Real pv = 0.0;                       // TRS PV = PV_asset - PV_funding (receiver 视角)
    Real pv_credit_loss = 0.0;           // 信用损失项 PV (receiver 支付的 LGD)
    Real pv_asset_return = 0.0;          // 资产回报项 PV (不含信用损失)
    // 各期现金流明细 (用于诊断与敏感性分析)
    std::vector<Real> asset_leg_cashflows;   // 各期 asset leg 现金流 (未折现)
    std::vector<Real> funding_leg_cashflows;  // 各期 funding leg 现金流 (未折现)
    std::vector<Real> payment_times;         // 支付时间表
};

// ============ TRS 定价 (确定性远期) ============
// 使用远期价格 F(t) = S0 * exp((r-q)*t) 进行确定性定价
// 信用风险通过 CreditCurve (hazard rate + recovery) 反映
inline TRSResult price_trs(
        const TRSConfig& cfg,
        const ZeroCurve& discount_curve,
        const CreditCurve& credit_curve) {
    TRSConfig cfg_validated = cfg;
    cfg_validated.generate_schedule();
    cfg_validated.validate();

    TRSResult result;
    result.payment_times = cfg_validated.payment_times;

    const Size n = cfg_validated.n_payments;
    const Real N = cfg_validated.notional;
    const Real S0 = cfg_validated.S0;
    const Real q = cfg_validated.asset_yield;
    const Real spread = cfg_validated.funding_spread;
    const Real R = credit_curve.recovery_rate();

    result.asset_leg_cashflows.assign(n, 0.0);
    result.funding_leg_cashflows.assign(n, 0.0);

    // 起始时间 t_0 (通常为 0)
    Real t_prev = cfg_validated.start_time;
    Real S_prev = S0;  // F(t_0) = S0 (远期价格在 t=0 等于现货)

    // 如果 start_time > 0, 需要用远期价格
    if (t_prev > 0.0) {
        // F(t_0) = S0 * exp((r - q) * t_0), 但 r 需从曲线获取
        // 使用 zero_rate(t_0) 作为 r 的代理
        Real r0 = discount_curve.zero_rate(t_prev);
        S_prev = S0 * std::exp((r0 - q) * t_prev);
    }

    Real pv_asset = 0.0;
    Real pv_funding = 0.0;
    Real pv_credit_loss = 0.0;
    Real pv_asset_return = 0.0;

    for (Size i = 0; i < n; ++i) {
        Real t_i = cfg_validated.payment_times[i];
        Real tau_i = t_i - t_prev;  // 年份分数
        if (tau_i <= 0.0) {
            throw std::invalid_argument("TRS: non-positive year fraction at period " + std::to_string(i));
        }

        // 远期价格 F(t_i) = S0 * exp((r - q) * t_i)
        // 但 r 沿曲线变化, 用 t_i 处的零息率
        Real r_i = discount_curve.zero_rate(t_i);
        Real F_i = S0 * std::exp((r_i - q) * t_i);

        // 区间总回报 (价格回报 + 票息)
        Real price_return = (F_i / S_prev) - 1.0;
        Real income_return = q * tau_i;
        Real total_return = price_return + income_return;

        // 远期利率 L(t_prev, t_i) = (P(t_prev)/P(t_i) - 1) / tau_i
        Real L_i = discount_curve.forward_rate(t_prev, t_i);

        // 折现因子
        Real P_i = discount_curve.discount_factor(t_i);
        Real t_mid = 0.5 * (t_prev + t_i);
        Real P_mid = discount_curve.discount_factor(t_mid);

        // 生存概率与违约概率
        Real Q_i = credit_curve.survival_prob(t_i);          // Q(0, t_i)
        Real Q_prev = credit_curve.survival_prob(t_prev);    // Q(0, t_{i-1})
        Real PD_interval = Q_prev - Q_i;                      // PD(t_{i-1}, t_i)
        if (PD_interval < 0.0) PD_interval = 0.0;             // 数值保护

        // ===== Asset leg (receiver: 正=收入) =====
        // 生存: 完整回报
        Real cf_asset_survival = total_return * N;
        // 违约: 信用损失 (receiver 支付 LGD)
        Real cf_credit_loss = -(1.0 - R) * N;  // 负号: receiver 支出

        Real pv_asset_survival = Q_i * cf_asset_survival * P_i;
        Real pv_default_loss = PD_interval * cf_credit_loss * P_mid;

        result.asset_leg_cashflows[i] = cf_asset_survival + PD_interval * cf_credit_loss;
        pv_asset_return += pv_asset_survival;
        pv_credit_loss += pv_default_loss;
        pv_asset += pv_asset_survival + pv_default_loss;

        // ===== Funding leg (receiver: 正=支出) =====
        // 生存: 完整融资支付
        Real cf_funding_survival = (L_i + spread) * tau_i * N;
        // 违约: 半期应计融资
        Real cf_funding_default = (L_i + spread) * (tau_i * 0.5) * N;

        Real pv_funding_survival = Q_i * cf_funding_survival * P_i;
        Real pv_funding_default = PD_interval * cf_funding_default * P_mid;

        result.funding_leg_cashflows[i] = cf_funding_survival + PD_interval * cf_funding_default;
        pv_funding += pv_funding_survival + pv_funding_default;

        // 更新前一期变量
        t_prev = t_i;
        S_prev = F_i;
    }

    result.pv_asset_leg = pv_asset;
    result.pv_funding_leg = pv_funding;
    result.pv_credit_loss = pv_credit_loss;
    result.pv_asset_return = pv_asset_return;
    result.pv = pv_asset - pv_funding;

    // 支付方视角: 取反
    if (!cfg_validated.is_receiver) {
        result.pv_asset_leg = -result.pv_asset_leg;
        result.pv_funding_leg = -result.pv_funding_leg;
        result.pv = -result.pv;
        result.pv_credit_loss = -result.pv_credit_loss;
        result.pv_asset_return = -result.pv_asset_return;
        for (auto& cf : result.asset_leg_cashflows) cf = -cf;
        for (auto& cf : result.funding_leg_cashflows) cf = -cf;
    }

    return result;
}

// ============ TRS 等价融资利差 (Par Spread) ============
// 求解使 TRS PV = 0 的 funding_spread
// PV(spread) = [PV_asset - PV_funding_excluding_spread] - spread * funding_pv01 = 0
//   funding_pv01 = Σ_i [Q(0,t_i)*τ_i*N*P(0,t_i) + PD(t_{i-1},t_i)*τ_i/2*N*P(0,t_mid)]
//   par_spread = [PV_asset - PV_funding_excluding_spread] / funding_pv01
inline Real trs_par_spread(
        const TRSConfig& cfg,
        const ZeroCurve& discount_curve,
        const CreditCurve& credit_curve) {
    TRSConfig cfg_validated = cfg;
    cfg_validated.generate_schedule();
    cfg_validated.validate();

    const Size n = cfg_validated.n_payments;
    const Real N = cfg_validated.notional;
    const Real S0 = cfg_validated.S0;
    const Real q = cfg_validated.asset_yield;
    const Real R = credit_curve.recovery_rate();

    Real t_prev = cfg_validated.start_time;
    Real S_prev = S0;
    if (t_prev > 0.0) {
        Real r0 = discount_curve.zero_rate(t_prev);
        S_prev = S0 * std::exp((r0 - q) * t_prev);
    }

    Real pv_asset = 0.0;
    Real pv_funding_base = 0.0;  // 不含 spread 的 funding PV
    Real funding_pv01 = 0.0;     // unit spread 的 funding PV (risky annuity)

    for (Size i = 0; i < n; ++i) {
        Real t_i = cfg_validated.payment_times[i];
        Real tau_i = t_i - t_prev;
        if (tau_i <= 0.0) {
            throw std::invalid_argument("TRS par_spread: non-positive year fraction");
        }

        Real r_i = discount_curve.zero_rate(t_i);
        Real F_i = S0 * std::exp((r_i - q) * t_i);
        Real price_return = (F_i / S_prev) - 1.0;
        Real total_return = price_return + q * tau_i;

        Real L_i = discount_curve.forward_rate(t_prev, t_i);
        Real P_i = discount_curve.discount_factor(t_i);
        Real t_mid = 0.5 * (t_prev + t_i);
        Real P_mid = discount_curve.discount_factor(t_mid);

        Real Q_i = credit_curve.survival_prob(t_i);
        Real Q_prev = credit_curve.survival_prob(t_prev);
        Real PD_interval = std::max(Q_prev - Q_i, 0.0);

        // Asset leg (同 price_trs)
        pv_asset += Q_i * total_return * N * P_i
                  + PD_interval * (-(1.0 - R)) * N * P_mid;

        // Funding leg (不含 spread)
        pv_funding_base += Q_i * L_i * tau_i * N * P_i
                        + PD_interval * L_i * (tau_i * 0.5) * N * P_mid;

        // funding_pv01 (unit spread 的 risky annuity)
        funding_pv01 += Q_i * tau_i * N * P_i
                      + PD_interval * (tau_i * 0.5) * N * P_mid;

        t_prev = t_i;
        S_prev = F_i;
    }

    if (funding_pv01 <= 0.0) {
        throw std::runtime_error("trs_par_spread: funding_pv01 <= 0");
    }

    // par_spread = (PV_asset - PV_funding_base) / funding_pv01
    // 注意符号: funding 是 receiver 的支出, spread 增加支出 → PV 降低
    Real par_spread = (pv_asset - pv_funding_base) / funding_pv01;
    return par_spread;
}

}  // namespace v1
}  // namespace cpphub
