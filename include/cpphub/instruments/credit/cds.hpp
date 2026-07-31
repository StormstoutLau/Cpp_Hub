#pragma once
// SOURCE: O'Kane (2008) "Modelling Single-Name and Multi-Name Credit Derivatives" Ch.5-7
// SOURCE: Schönbucher (2003) "Credit Derivatives Pricing Models" Ch.6
// SOURCE: Hull & White (2000) "Valuing Credit Default Swaps I: No Counterparty Default Risk"
// SOURCE: ISDA (2003) "Credit Derivatives Definitions"
// 模块: 信用违约互换 (Credit Default Swap, CDS)
//
// ==================== CDS 结构 ====================
//
// CDS 是信用衍生品中最基础的产品: 保护买方 (protection buyer) 定期支付保费 (premium),
// 换取参考实体 (reference entity) 违约时获得赔付 (protection).
//
// 现金流 (从保护买方视角):
//   - Premium leg (付费端): 在支付日 t_1, ..., t_n 支付 s * τ_i * N, s 为 CDS spread
//     若违约发生在 (t_{i-1}, t_i], 还需支付应计保费 s * (τ_default - τ_{i-1}) * N
//   - Protection leg (保护端): 违约时获得 (1 - R) * N, R 为回收率
//
// ==================== 定价公式 (从保护买方视角) ====================
//
// Premium leg PV (含应计):
//   PV_premium = s * N * Σ_i τ_i * [Q(0,t_i) * P(0,t_i) + 0.5 * PD(t_{i-1}, t_i) * P(0, t_mid_i)]
//   其中 Q 为生存概率, PD 为区间违约概率, P 为折现因子
//   第一项: 存活到 t_i 支付完整保费
//   第二项: 在 (t_{i-1}, t_i] 违约, 假设违约在区间中点, 支付半期应计保费
//
// Protection leg PV:
//   PV_protection = (1 - R) * N * Σ_i PD(t_{i-1}, t_i) * P(0, t_mid_i)
//   假设违约在区间中点, 折现到 t_mid_i
//
// CDS PV (buyer) = PV_protection - PV_premium
//   buyer 买入保护, 支付保费, 获得赔付 → PV > 0 表示保护买方获益
//
// Par spread (让 PV = 0 的 s):
//   s_par = PV_protection / Risky_pv01
//   Risky_pv01 = Σ_i τ_i * [Q(0,t_i) * P(0,t_i) + 0.5 * PD(t_{i-1}, t_i) * P(0, t_mid_i)]
//   (即 unit spread 的 premium leg PV, 也称 risky annuity)
//
// ==================== CDS Curve Bootstrap ====================
//
// 给定一组 CDS quotes (maturity, spread), 逐点 bootstrap hazard rate term structure:
//   第 j 个 quote: 已知前 j-1 个 hazard rate 段, 求第 j 段 h_j 使 CDS PV = 0
//   PV_protection(h_j) = s_j * Risky_pv01(h_j)
//   由于 h_j 影响第 j 段的 Q 和 PD, 可用二分法求解 (h_j ≥ 0)
//
// ==================== 约定 ====================
//
// - 支付频率: 季度 (ACT/365, τ ≈ 0.25) 或 半年度 (τ ≈ 0.5), ISDA 标准
// - 违约时点假设: 区间中点 (mid-point convention)
// - 折现: OIS discount curve (post-crisis 标准)
// - 应计保费: 全额应计 (full first coupon), 假设违约在区间中点支付半期

#include "cpphub/core/types.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include "cpphub/instruments/credit/credit_curve.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ CDS 配置 ============
struct CDSConfig {
    Real notional = 1.0;                // 名义本金
    Real spread = 0.01;                  // CDS spread (年化, 如 0.01 = 100bp)
    Real start_time = 0.0;               // 起始时间 (0 = spot)
    Real maturity = 5.0;                 // 到期时间 (年)
    Size n_premiums = 20;                // 保费支付次数 (如 5Y 季度 = 20)
    bool is_buyer = true;                // true = 保护买方, false = 保护卖方
    // 保费支付时间表 (自动生成等间隔), 也可外部指定
    std::vector<Real> payment_times;     // 若非空则用此, 否则按等间隔生成
    std::vector<Real> year_fractions;    // 若非空则用此, 否则按 1/n_premiums_per_year

    // 等间隔支付时间表生成 (若 payment_times 为空)
    void generate_schedule() {
        if (!payment_times.empty()) return;
        if (n_premiums == 0) {
            throw std::invalid_argument("CDSConfig: n_premiums must be positive");
        }
        Real tau = (maturity - start_time) / static_cast<Real>(n_premiums);
        payment_times.resize(n_premiums);
        year_fractions.assign(n_premiums, tau);
        for (Size i = 0; i < n_premiums; ++i) {
            payment_times[i] = start_time + static_cast<Real>(i + 1) * tau;
        }
    }
};

// ============ CDS 定价器 ============
class CreditDefaultSwap {
public:
    explicit CreditDefaultSwap(CDSConfig cfg) : cfg_(std::move(cfg)) {
        cfg_.generate_schedule();
        validate();
    }

    // Premium leg PV (含应计, 从保护买方视角为负, 因为买方支付)
    // PV_premium = s * N * Σ τ_i * [Q(t_i)*P(t_i) + 0.5*PD(t_{i-1},t_i)*P(t_mid_i)]
    Real premium_leg_pv(const CreditCurve& credit, const ZeroCurve& discount) const {
        Real pv = 0.0;
        Real s = cfg_.spread;
        Real N = cfg_.notional;
        for (Size i = 0; i < cfg_.payment_times.size(); ++i) {
            Real t_i = cfg_.payment_times[i];
            Real t_prev = (i == 0) ? cfg_.start_time : cfg_.payment_times[i - 1];
            Real tau = cfg_.year_fractions[i];
            Real t_mid = 0.5 * (t_prev + t_i);
            Real Q_i = credit.survival_prob(t_i);
            Real PD_interval = credit.default_prob(t_prev, t_i);
            Real P_i = discount.discount_factor(t_i);
            Real P_mid = discount.discount_factor(t_mid);
            // 完整支付 (存活到 t_i) + 半期应计 (区间内违约)
            pv += s * N * tau * (Q_i * P_i + 0.5 * PD_interval * P_mid);
        }
        return pv;
    }

    // Protection leg PV (从保护买方视角为正, 因为买方获得赔付)
    // PV_protection = (1-R) * N * Σ PD(t_{i-1},t_i) * P(t_mid_i)
    Real protection_leg_pv(const CreditCurve& credit, const ZeroCurve& discount) const {
        Real pv = 0.0;
        Real LGD = credit.lgd();
        Real N = cfg_.notional;
        for (Size i = 0; i < cfg_.payment_times.size(); ++i) {
            Real t_i = cfg_.payment_times[i];
            Real t_prev = (i == 0) ? cfg_.start_time : cfg_.payment_times[i - 1];
            Real t_mid = 0.5 * (t_prev + t_i);
            Real PD_interval = credit.default_prob(t_prev, t_i);
            Real P_mid = discount.discount_factor(t_mid);
            pv += LGD * N * PD_interval * P_mid;
        }
        return pv;
    }

    // Risky PV01 (unit spread 的 premium leg PV, 也称 risky annuity)
    // = Σ τ_i * [Q(t_i)*P(t_i) + 0.5*PD(t_{i-1},t_i)*P(t_mid_i)]
    Real risky_pv01(const CreditCurve& credit, const ZeroCurve& discount) const {
        Real pv = 0.0;
        for (Size i = 0; i < cfg_.payment_times.size(); ++i) {
            Real t_i = cfg_.payment_times[i];
            Real t_prev = (i == 0) ? cfg_.start_time : cfg_.payment_times[i - 1];
            Real tau = cfg_.year_fractions[i];
            Real t_mid = 0.5 * (t_prev + t_i);
            Real Q_i = credit.survival_prob(t_i);
            Real PD_interval = credit.default_prob(t_prev, t_i);
            Real P_i = discount.discount_factor(t_i);
            Real P_mid = discount.discount_factor(t_mid);
            pv += tau * (Q_i * P_i + 0.5 * PD_interval * P_mid);
        }
        return pv;
    }

    // Par spread (让 PV = 0 的 spread): s_par = PV_protection / Risky_pv01
    Real par_spread(const CreditCurve& credit, const ZeroCurve& discount) const {
        Real rpv01 = risky_pv01(credit, discount);
        if (rpv01 <= 0.0) {
            throw std::runtime_error("CDS::par_spread: risky_pv01 must be positive");
        }
        return protection_leg_pv(credit, discount) / rpv01;
    }

    // CDS PV (从保护买方视角): PV = PV_protection - PV_premium
    // buyer: PV > 0 表示保护买方获益 (买入保护划算)
    // seller: PV = -buyer_PV
    Real pv(const CreditCurve& credit, const ZeroCurve& discount) const {
        Real net = protection_leg_pv(credit, discount) - premium_leg_pv(credit, discount);
        return cfg_.is_buyer ? net : -net;
    }

    // 访问器
    const CDSConfig& config() const noexcept { return cfg_; }
    Real notional() const noexcept { return cfg_.notional; }
    Real spread() const noexcept { return cfg_.spread; }
    Real maturity() const noexcept { return cfg_.maturity; }
    bool is_buyer() const noexcept { return cfg_.is_buyer; }

private:
    CDSConfig cfg_;

    void validate() const {
        if (cfg_.notional <= 0.0) {
            throw std::invalid_argument("CDS: notional must be positive");
        }
        if (cfg_.spread < 0.0) {
            throw std::invalid_argument("CDS: spread must be non-negative");
        }
        if (cfg_.maturity <= cfg_.start_time) {
            throw std::invalid_argument("CDS: maturity must be > start_time");
        }
        if (cfg_.payment_times.empty()) {
            throw std::invalid_argument("CDS: payment_times empty");
        }
        if (cfg_.year_fractions.size() != cfg_.payment_times.size()) {
            throw std::invalid_argument("CDS: year_fractions size mismatch");
        }
        for (Real t : cfg_.payment_times) {
            if (t <= 0.0) {
                throw std::invalid_argument("CDS: payment times must be positive");
            }
        }
    }
};

// ============ CDS Curve Bootstrapper ============
// 给定一组 CDS quotes (maturity, spread), 逐点 bootstrap hazard rate term structure.
// 假设: 每个 quote 对应一个新的 hazard rate 段 (分段常数), 支付频率按等间隔.
//
// 输入: CDS quotes + 回收率 + 折现曲线 + 每年支付频率
// 输出: CreditCurve (PDCurve + recovery), hazard rate 分段常数
class CDSCurveBootstrapper {
public:
    struct CDSQuote {
        Real maturity;   // CDS 期限 (年)
        Real spread;     // CDS spread (年化, 如 0.01 = 100bp)
    };

    CDSCurveBootstrapper(std::vector<CDSQuote> quotes, Real recovery_rate,
                           const ZeroCurve& discount_curve, Size premiums_per_year = 4)
        : quotes_(std::move(quotes)),
          recovery_(recovery_rate),
          discount_(discount_curve),
          freq_(premiums_per_year) {
        validate();
    }

    // Bootstrap 信用曲线
    // 逐点求解: 第 j 个 quote, 已知前 j 段 hazard rate, 二分法求 h_j 使 spot CDS (0→mat) PV = 0
    // 关键: 每次求解的是 spot CDS (从 0 到 mat), 而非 forward CDS (prev_mat → mat)
    //       这样 bootstrap 出的曲线对每个 quote 期限的 spot CDS 都有 PV ≈ 0
    CreditCurve bootstrap() const {
        const Real h_tol = 1e-12;   // hazard rate 收敛容差
        const Real pv_tol = 1e-10;  // PV 收敛容差 (单位 notional)
        const Size max_iter = 100;

        std::vector<Real> times;        // hazard rate 时间节点 (T_1, ..., T_n)
        std::vector<Real> hazard_rates; // 分段常数 hazard rate (h_0, ..., h_{n-1})

        for (Size j = 0; j < quotes_.size(); ++j) {
            Real mat = quotes_[j].maturity;
            Real s = quotes_[j].spread;

            // 二分法求第 j 段 hazard rate h_j ∈ [0, h_max]
            // h_j 影响 (prev_mat, mat] 区间的 Q 和 PD, 从而影响 spot CDS (0→mat) PV
            Real h_lo = 0.0;
            Real h_hi = 1.0;  // 100% hazard rate 上限 (足够大)
            // 扩展上限直到 spot CDS PV 在 h_hi 下 < 0 (protection > premium)
            for (Size iter = 0; iter < 50; ++iter) {
                Real pv_hi = spot_cds_pv_with_hazard_(h_hi, times, hazard_rates, mat, s);
                if (pv_hi < 0.0) break;
                h_hi *= 2.0;
                if (h_hi > 1e6) break;  // 防止无限扩展
            }

            Real h_mid = 0.0;
            for (Size iter = 0; iter < max_iter; ++iter) {
                h_mid = 0.5 * (h_lo + h_hi);
                Real pv_mid = spot_cds_pv_with_hazard_(h_mid, times, hazard_rates, mat, s);
                if (std::abs(pv_mid) < pv_tol) break;
                // PV > 0 (protection > premium) → h 太大 (违约多), 降低 h
                // PV < 0 (protection < premium) → h 太小, 提高 h
                if (pv_mid > 0.0) {
                    h_hi = h_mid;
                } else {
                    h_lo = h_mid;
                }
                if (h_hi - h_lo < h_tol) break;
            }

            // 记录第 j 段 hazard rate
            times.push_back(mat);
            hazard_rates.push_back(h_mid);
        }

        // 构造 CreditCurve: PDCurve 分段约定 hazard_rates_[i] 对应 [times_[i], times_[i+1})
        // 目标分段:
        //   [0.01, T_1): h_0
        //   [T_1, T_2): h_1
        //   ...
        //   [T_{n-1}, T_n): h_{n-1}
        //   [T_n, ∞): h_{n-1}  (最后一段延伸)
        std::vector<Real> pd_times;
        std::vector<Real> pd_hrs;
        pd_times.push_back(0.01);
        pd_hrs.push_back(hazard_rates[0]);  // 段 [0.01, T_1)
        const Size n_segments = times.size();
        for (Size i = 0; i < n_segments; ++i) {
            pd_times.push_back(times[i]);  // T_{i+1}
            // 段 [T_{i+1}, T_{i+2}) 用 h_{i+1}; 最后一段 [T_n, ∞) 用 h_{n-1}
            if (i + 1 < n_segments) {
                pd_hrs.push_back(hazard_rates[i + 1]);
            } else {
                pd_hrs.push_back(hazard_rates[i]);  // 无穷段重复最后 hazard
            }
        }
        PDCurve pd(std::move(pd_times), std::move(pd_hrs));
        return CreditCurve(std::move(pd), recovery_);
    }

    const std::vector<CDSQuote>& quotes() const noexcept { return quotes_; }
    Real recovery_rate() const noexcept { return recovery_; }
    Size premiums_per_year() const noexcept { return freq_; }

private:
    std::vector<CDSQuote> quotes_;
    Real recovery_;
    const ZeroCurve& discount_;
    Size freq_;

    void validate() const {
        if (quotes_.empty()) {
            throw std::invalid_argument("CDSCurveBootstrapper: at least 1 quote required");
        }
        for (Size i = 0; i < quotes_.size(); ++i) {
            if (quotes_[i].maturity <= 0.0) {
                throw std::invalid_argument("CDSCurveBootstrapper: maturity must be positive");
            }
            if (quotes_[i].spread < 0.0) {
                throw std::invalid_argument("CDSCurveBootstrapper: spread must be non-negative");
            }
            if (i > 0 && quotes_[i].maturity <= quotes_[i - 1].maturity) {
                throw std::invalid_argument("CDSCurveBootstrapper: maturities must be strictly increasing");
            }
        }
        if (recovery_ < 0.0 || recovery_ >= 1.0) {
            throw std::invalid_argument("CDSCurveBootstrapper: recovery must be in [0, 1)");
        }
        if (freq_ == 0) {
            throw std::invalid_argument("CDSCurveBootstrapper: premiums_per_year must be positive");
        }
    }

    // 计算给定第 j 段 hazard rate h_j 下的 spot CDS (0→mat) PV (buyer, notional=1)
    // 已知前 j 段 hazard rate (prev_times=[T_1,...,T_j], prev_hrs=[h_0,...,h_{j-1}])
    // 第 j 段: (T_j, mat], hazard = h_j (新段)
    //
    // PDCurve 分段约定: hazard_rates_[i] 对应段 [times_[i], times_[i+1})
    // 目标分段:
    //   [0.01, T_1): h_0
    //   [T_1, T_2): h_1
    //   ...
    //   [T_j, mat): h_j   (新段)
    //   [mat, ∞): h_j     (无穷段, 不影响有限区间计算)
    //
    // CDS 配置: spot start (0 → mat), n_payments = round(mat * freq_), spread = s
    Real spot_cds_pv_with_hazard_(Real h_j, const std::vector<Real>& prev_times,
                                       const std::vector<Real>& prev_hrs,
                                       Real mat, Real s) const {
        const Size j = prev_times.size();  // 已知前 j 段
        // 构造临时 PDCurve
        std::vector<Real> pd_times;
        std::vector<Real> pd_hrs;
        pd_times.push_back(0.01);
        pd_hrs.push_back(prev_hrs.empty() ? h_j : prev_hrs[0]);  // 段 [0.01, T_1)
        for (Size i = 0; i < j; ++i) {
            pd_times.push_back(prev_times[i]);  // T_{i+1}
            // 段 [T_{i+1}, T_{i+2}) 用 h_{i+1}; 若 i+1 == j (即 T_{i+1}=T_j), 段 [T_j, mat) 用 h_j
            if (i + 1 < j) {
                pd_hrs.push_back(prev_hrs[i + 1]);
            } else {
                pd_hrs.push_back(h_j);
            }
        }
        // 添加 mat 节点 (无穷段 hazard = h_j)
        pd_times.push_back(mat);
        pd_hrs.push_back(h_j);

        PDCurve pd(std::move(pd_times), std::move(pd_hrs));
        CreditCurve credit(std::move(pd), recovery_);

        // 构造 spot CDS (0 → mat), 支付频率 freq_, spread = s, notional = 1
        Size n_payments = std::max<Size>(1, static_cast<Size>(std::round(mat * freq_)));
        CDSConfig cfg;
        cfg.notional = 1.0;
        cfg.spread = s;
        cfg.start_time = 0.0;
        cfg.maturity = mat;
        cfg.n_premiums = n_payments;
        cfg.is_buyer = true;
        cfg.payment_times.clear();
        cfg.year_fractions.clear();
        cfg.generate_schedule();

        CreditDefaultSwap cds(cfg);
        return cds.pv(credit, discount_);
    }
};

// ============ 便捷工厂函数 ============
inline CDSConfig make_cds(Real notional, Real spread, Real maturity,
                            Size premiums_per_year = 4, bool is_buyer = true) {
    CDSConfig cfg;
    cfg.notional = notional;
    cfg.spread = spread;
    cfg.start_time = 0.0;
    cfg.maturity = maturity;
    cfg.n_premiums = static_cast<Size>(std::round(maturity * premiums_per_year));
    cfg.is_buyer = is_buyer;
    cfg.payment_times.clear();
    cfg.year_fractions.clear();
    cfg.generate_schedule();
    return cfg;
}

}  // namespace v1
}  // namespace cpphub
