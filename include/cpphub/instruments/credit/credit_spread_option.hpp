#pragma once
// SOURCE: Schönbucher (2003) "Credit Derivatives Pricing Models" Ch.10-11 (Spread Options)
// SOURCE: O'Kane (2008) "Modelling Single-Name and Multi-Name Credit Derivatives" Ch.8 (Options)
// SOURCE: Hull & White (2003) "The Valuation of Credit-Sensitive Notes"
// SOURCE: Black (1976) "The Pricing of Commodity Contracts" J. Financial Economics 3, 167-179
// 模块: 信用利差期权 (Credit Spread Option, CSO)
//
// ==================== CSO 结构 ====================
//
// 信用利差期权赋予持有人在期权到期 T 时, 以预先约定利差 K 进入一笔 CDS 的权利 (而非义务).
//
// 两类期权:
//   - Payer (Call): 买入保护权 — 有权以 K 支付保费换取保护
//     若市场利差 S_T > K, 期权持有人获利 (以低于市场的 K 买入保护)
//     Payoff = max(S_T - K, 0) * N * τ   (τ 为基础 CDS 的 risk period)
//   - Receiver (Put): 卖出保护权 — 有权以 K 收取保费提供保护
//     若市场利差 S_T < K, 期权持有人获利 (以高于市场的 K 卖出保护)
//     Payoff = max(K - S_T, 0) * N * τ
//
// ==================== 定价公式 (Black 1976 模型) ====================
//
// 假设远期利差 F_T 在 T-forward 测度下服从对数正态分布:
//   ln(F_T) ~ N(ln(F_0) - 0.5σ²T, σ²T)
//
// 非敲出 (non-knock-out) 期权:
//   PV_call = exp(-rT) * [F_0 * N(d1) - K * N(d2)] * N * τ
//   PV_put  = exp(-rT) * [K * N(-d2) - F_0 * N(-d1)] * N * τ
//   d1 = [ln(F_0/K) + 0.5σ²T] / (σ√T)
//   d2 = d1 - σ√T
//
// 敲出 (knock-out) 期权 — 参考实体在 T 前违约则期权作废:
//   PV_call = Q(0,T) * exp(-rT) * [F_0 * N(d1) - K * N(d2)] * N * τ
//   PV_put  = Q(0,T) * exp(-rT) * [K * N(-d2) - F_0 * N(-d1)] * N * τ
//   其中 Q(0,T) = 生存概率 (从 CreditCurve 获取)
//
// ==================== Greeks ====================
//
// Delta (∂PV/∂F_0):
//   call: Q * exp(-rT) * N(d1) * N_notional * τ
//   put:  -Q * exp(-rT) * N(-d1) * N_notional * τ
//
// Gamma (∂²PV/∂F_0²):
//   call/put: Q * exp(-rT) * n(d1) / (F_0 * σ * √T) * N_notional * τ
//
// Vega (∂PV/∂σ):
//   call/put: Q * exp(-rT) * F_0 * n(d1) * √T * N_notional * τ
//
// 其中 n(x) 为标准正态密度函数.
//
// ==================== 约定 ====================
//
// - 远期利差 F_0: 由市场报价或从信用曲线推得 (此处作为输入参数)
// - 波动率 σ: 年化, 对数正态波动率 (可从市场报价或历史数据估计)
// - τ (cds_tenor): 基础 CDS 期限, 用于放大 payoff (近似 risky annuity)
// - 折现: OIS discount curve
// - 信用曲线: CreditCurve (提供生存概率 Q(0,T))

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include "cpphub/instruments/credit/credit_curve.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>

namespace cpphub {
inline namespace v1 {

// ============ 信用利差期权配置 ============
struct CSOConfig {
    Real notional = 1.0;                // 名义本金
    Real strike_spread = 0.01;          // K: 执行利差 (年化, 如 0.01 = 100bp)
    Real option_maturity = 1.0;         // T: 期权到期时间 (年)
    Real cds_tenor = 5.0;               // τ: 基础 CDS 期限 (年, 近似 risky annuity)
    Real forward_spread = 0.01;         // F_0: 远期利差 (年化)
    Real spread_volatility = 0.3;       // σ: 利差波动率 (年化, 对数正态)
    bool is_payer = true;               // true = payer (call), false = receiver (put)
    bool is_knock_out = true;           // true = 标的违约后期权作废

    // 验证配置
    void validate() const {
        if (notional <= 0.0) {
            throw std::invalid_argument("CSOConfig: notional must be positive");
        }
        if (strike_spread < 0.0) {
            throw std::invalid_argument("CSOConfig: strike_spread must be non-negative");
        }
        if (option_maturity <= 0.0) {
            throw std::invalid_argument("CSOConfig: option_maturity must be positive");
        }
        if (cds_tenor <= 0.0) {
            throw std::invalid_argument("CSOConfig: cds_tenor must be positive");
        }
        if (forward_spread < 0.0) {
            throw std::invalid_argument("CSOConfig: forward_spread must be non-negative");
        }
        if (spread_volatility < 0.0) {
            throw std::invalid_argument("CSOConfig: spread_volatility must be non-negative");
        }
    }
};

// ============ 信用利差期权结果 ============
struct CSOResult {
    Real pv = 0.0;                  // 期权现值
    Real d1 = 0.0;                  // Black d1
    Real d2 = 0.0;                  // Black d2
    Real delta = 0.0;               // ∂PV/∂F_0
    Real gamma = 0.0;               // ∂²PV/∂F_0²
    Real vega = 0.0;                // ∂PV/∂σ (per 1.0 vol change, 即 100%)
    Real theta = 0.0;               // ∂PV/∂T (per year)
    Real survival_prob = 1.0;       // Q(0, T) (敲出期权使用)
    Real discount_factor = 1.0;     // P(0, T)
    Real intrinsic_value = 0.0;     // 内在价值 (未折现 payoff 的现值)
};

// ============ 信用利差期权定价 (Black 1976) ============
// 使用 Black 模型对远期利差定价, 支持敲出 (knock-out) 与非敲出两种形式.
inline CSOResult price_credit_spread_option(
        const CSOConfig& cfg,
        const ZeroCurve& discount_curve,
        const CreditCurve& credit_curve) {
    CSOConfig cfg_valid = cfg;
    cfg_valid.validate();

    CSOResult result;

    const Real N = cfg_valid.notional;
    const Real K = cfg_valid.strike_spread;
    const Real T = cfg_valid.option_maturity;
    const Real tau = cfg_valid.cds_tenor;
    const Real F = cfg_valid.forward_spread;
    const Real sigma = cfg_valid.spread_volatility;

    // 折现因子与生存概率
    result.discount_factor = discount_curve.discount_factor(T);
    result.survival_prob = cfg_valid.is_knock_out ? credit_curve.survival_prob(T) : 1.0;

    const Real df = result.discount_factor;
    const Real Q = result.survival_prob;

    // 处理退化情况: F=0 或 K=0 或 sigma=0
    if (F <= 0.0 && K <= 0.0) {
        result.pv = 0.0;
        return result;
    }

    // 内在价值 (未折现)
    if (cfg_valid.is_payer) {
        result.intrinsic_value = std::max(F - K, 0.0) * N * tau;
    } else {
        result.intrinsic_value = std::max(K - F, 0.0) * N * tau;
    }

    // sigma=0: 期权等于内在价值的折现
    if (sigma < 1e-12) {
        result.pv = Q * df * result.intrinsic_value;
        return result;
    }

    // 计算 d1, d2
    const Real sqrtT = std::sqrt(T);
    const Real sig_sqrtT = sigma * sqrtT;

    if (F <= 0.0) {
        // F=0: call 无价值, put = K * exp(-rT) * Q * N * τ
        if (cfg_valid.is_payer) {
            result.pv = 0.0;
            result.d1 = -std::numeric_limits<Real>::infinity();
            result.d2 = -std::numeric_limits<Real>::infinity();
        } else {
            result.d1 = -std::numeric_limits<Real>::infinity();
            result.d2 = -std::numeric_limits<Real>::infinity();
            result.pv = Q * df * K * N * tau;
        }
        return result;
    }

    if (K <= 0.0) {
        // K=0: call = F * exp(-rT) * Q * N * τ, put 无价值
        if (cfg_valid.is_payer) {
            result.d1 = std::numeric_limits<Real>::infinity();
            result.d2 = std::numeric_limits<Real>::infinity();
            result.pv = Q * df * F * N * tau;
        } else {
            result.pv = 0.0;
            result.d1 = std::numeric_limits<Real>::infinity();
            result.d2 = std::numeric_limits<Real>::infinity();
        }
        return result;
    }

    result.d1 = (std::log(F / K) + 0.5 * sigma * sigma * T) / sig_sqrtT;
    result.d2 = result.d1 - sig_sqrtT;

    const Real N_d1 = normal_cdf(result.d1);
    const Real N_d2 = normal_cdf(result.d2);
    const Real N_minus_d1 = normal_cdf(-result.d1);
    const Real N_minus_d2 = normal_cdf(-result.d2);
    const Real n_d1 = normal_pdf(result.d1);

    // payoff 放大因子: N * τ (notional * cds_tenor)
    const Real scale = N * tau;

    if (cfg_valid.is_payer) {
        // Call: PV = Q * df * [F * N(d1) - K * N(d2)] * N * τ
        result.pv = Q * df * (F * N_d1 - K * N_d2) * scale;
        // Greeks
        result.delta = Q * df * N_d1 * scale;
        result.gamma = Q * df * n_d1 / (F * sig_sqrtT) * scale;
        result.vega = Q * df * F * n_d1 * sqrtT * scale;
        // Theta: ∂PV/∂T (注意: 此处 theta 为 ∂PV/∂T, 正值表示随时间增加 PV 增加)
        // 对 Black 模型: ∂PV/∂T = Q*df * [F*n(d1)*σ/(2√T) - r*(F*N(d1) - K*N(d2)) + h*(F*N(d1) - K*N(d2))] * scale
        // (h 为 hazard rate, 来自 Q 的导数)
        Real h_T = credit_curve.hazard_rate(T);
        Real r_T = discount_curve.zero_rate(T);
        Real black_core = F * N_d1 - K * N_d2;
        Real theta_black = F * n_d1 * sigma / (2.0 * sqrtT);
        Real theta_disc = -r_T * black_core;
        Real theta_surv = cfg_valid.is_knock_out ? -h_T * black_core : 0.0;
        // 注意: df = exp(-rT), d(df)/dT = -r * df
        // Q = exp(-H(T)), dQ/dT = -h(T) * Q
        result.theta = Q * df * (theta_black + theta_disc) * scale
                     + (cfg_valid.is_knock_out ? -h_T * result.pv : 0.0);
    } else {
        // Put: PV = Q * df * [K * N(-d2) - F * N(-d1)] * N * τ
        result.pv = Q * df * (K * N_minus_d2 - F * N_minus_d1) * scale;
        // Greeks
        result.delta = -Q * df * N_minus_d1 * scale;
        result.gamma = Q * df * n_d1 / (F * sig_sqrtT) * scale;  // 同 call
        result.vega = Q * df * F * n_d1 * sqrtT * scale;  // 同 call
        // Theta
        Real h_T = credit_curve.hazard_rate(T);
        Real r_T = discount_curve.zero_rate(T);
        Real black_core_put = K * N_minus_d2 - F * N_minus_d1;
        Real theta_black_put = -F * n_d1 * sigma / (2.0 * sqrtT);  // put 的 theta_black 为负
        Real theta_disc_put = -r_T * black_core_put;
        result.theta = Q * df * (theta_black_put + theta_disc_put) * scale
                     + (cfg_valid.is_knock_out ? -h_T * result.pv : 0.0);
    }

    return result;
}

// ============ 隐含波动率求解 (Bisection) ============
// 给定市场价格, 反解 Black 模型的隐含波动率
// 适用于非退化情形 (F>0, K>0, T>0)
inline Real credit_spread_implied_vol(
        const CSOConfig& cfg,
        const ZeroCurve& discount_curve,
        const CreditCurve& credit_curve,
        Real market_price,
        Real vol_tol = 1e-8,
        Size max_iter = 100) {
    CSOConfig cfg_valid = cfg;
    cfg_valid.validate();

    if (cfg_valid.forward_spread <= 0.0 || cfg_valid.strike_spread <= 0.0) {
        throw std::invalid_argument("implied_vol: require F>0 and K>0");
    }
    if (market_price < 0.0) {
        throw std::invalid_argument("implied_vol: market_price must be non-negative");
    }

    // 检查市场价是否在无套利区间内
    CSOConfig cfg_lower = cfg_valid;
    cfg_lower.spread_volatility = 0.0;
    Real price_lower = price_credit_spread_option(cfg_lower, discount_curve, credit_curve).pv;

    Real price_upper = 0.0;
    if (cfg_valid.is_payer) {
        // call 上界: Q * df * F * N * τ
        Real Q = cfg_valid.is_knock_out ? credit_curve.survival_prob(cfg_valid.option_maturity) : 1.0;
        Real df = discount_curve.discount_factor(cfg_valid.option_maturity);
        price_upper = Q * df * cfg_valid.forward_spread * cfg_valid.notional * cfg_valid.cds_tenor;
    } else {
        // put 上界: Q * df * K * N * τ
        Real Q = cfg_valid.is_knock_out ? credit_curve.survival_prob(cfg_valid.option_maturity) : 1.0;
        Real df = discount_curve.discount_factor(cfg_valid.option_maturity);
        price_upper = Q * df * cfg_valid.strike_spread * cfg_valid.notional * cfg_valid.cds_tenor;
    }

    if (market_price < price_lower - 1e-12 || market_price > price_upper + 1e-12) {
        throw std::runtime_error("implied_vol: market_price outside no-arbitrage bounds");
    }

    // 二分法
    Real vol_lo = 1e-8;
    Real vol_hi = 5.0;  // 500% vol 上限

    CSOConfig cfg_tmp = cfg_valid;
    cfg_tmp.spread_volatility = vol_lo;
    Real pv_lo = price_credit_spread_option(cfg_tmp, discount_curve, credit_curve).pv;
    if (std::abs(pv_lo - market_price) < vol_tol) return vol_lo;

    cfg_tmp.spread_volatility = vol_hi;
    Real pv_hi = price_credit_spread_option(cfg_tmp, discount_curve, credit_curve).pv;
    if (std::abs(pv_hi - market_price) < vol_tol) return vol_hi;

    for (Size iter = 0; iter < max_iter; ++iter) {
        Real vol_mid = 0.5 * (vol_lo + vol_hi);
        cfg_tmp.spread_volatility = vol_mid;
        Real pv_mid = price_credit_spread_option(cfg_tmp, discount_curve, credit_curve).pv;

        if (std::abs(pv_mid - market_price) < vol_tol) return vol_mid;

        // PV 是 σ 的单调递增函数 (对 call 和 put 都成立)
        if (pv_mid < market_price) {
            vol_lo = vol_mid;
        } else {
            vol_hi = vol_mid;
        }
        if (vol_hi - vol_lo < vol_tol) break;
    }

    return 0.5 * (vol_lo + vol_hi);
}

// ============ 信用利差期权平价关系 (Put-Call Parity) ============
// 对于相同 K, T, τ, F, σ 的 payer 和 receiver:
//   PV_call - PV_put = Q * df * (F - K) * N * τ
// 此函数验证平价关系是否成立 (返回差值, 应接近 0)
inline Real credit_spread_parity_check(
        const CSOConfig& cfg,
        const ZeroCurve& discount_curve,
        const CreditCurve& credit_curve) {
    CSOConfig cfg_payer = cfg;
    cfg_payer.is_payer = true;
    CSOConfig cfg_receiver = cfg;
    cfg_receiver.is_payer = false;

    Real pv_call = price_credit_spread_option(cfg_payer, discount_curve, credit_curve).pv;
    Real pv_put = price_credit_spread_option(cfg_receiver, discount_curve, credit_curve).pv;

    Real Q = cfg.is_knock_out ? credit_curve.survival_prob(cfg.option_maturity) : 1.0;
    Real df = discount_curve.discount_factor(cfg.option_maturity);
    Real parity = Q * df * (cfg.forward_spread - cfg.strike_spread) * cfg.notional * cfg.cds_tenor;

    return (pv_call - pv_put) - parity;
}

}  // namespace v1
}  // namespace cpphub
