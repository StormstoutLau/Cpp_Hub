#pragma once
// SOURCE: Schönbucher (2003) "Credit Derivatives Pricing Models" Ch.2
// SOURCE: O'Kane (2008) "Modelling Single-Name and Multi-Name Credit Derivatives"
// SOURCE: Brigo & Mercurio (2006) "Interest Rate Models" Ch.21 (Credit)
// 模块: 信用曲线 (CreditCurve) + 违约概率曲线 (PDCurve)
//
// ==================== 数学框架 ====================
//
// 违约时间 τ 的建模: 生存概率 Q(0,T) = P(τ > T) = exp(-H(T))
//   H(T) = ∫_0^T h(s) ds  (累积 hazard function)
//   h(s) 为瞬时 hazard rate (违约强度), 分段常数假设
//
// 区间违约概率: PD(t1, t2) = Q(0,t1) - Q(0,t2) = P(t1 < τ ≤ t2)
// 累积违约概率: PD(0, T) = 1 - Q(0, T)
//
// 从 CDS spread 单点近似 (flat curve, 已知 recovery R):
//   s ≈ h * (1 - R)  →  h = s / (1 - R)
//   (Hazard rate 与 CDS spread 的线性近似, 忽略折现与应计效应)
//
// CreditCurve = PDCurve + recovery rate + 名字 (reference entity)
//   - 损失给定违约 LGD = 1 - R
//   - 预期损失 EL(T) = LGD * PD(0, T)

#include "cpphub/core/types.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>

namespace cpphub {
inline namespace v1 {

// ============ 违约概率曲线 (PDCurve) ============
// 分段常数 hazard rate, 提供生存/违约概率 + 区间违约概率
// (从 xva.hpp 迁移至此, 作为信用模块的基础设施)
class PDCurve {
public:
    PDCurve() = default;

    PDCurve(std::vector<Real> times, std::vector<Real> hazard_rates)
        : times_(std::move(times)), hazard_rates_(std::move(hazard_rates)) {
        validate();
    }

    // 生存概率 Q(0, T) = exp(-H(T)), H(T) = ∫_0^T h(s) ds
    Real survival_prob(Real T) const {
        if (T <= 0.0) return 1.0;
        return std::exp(-cumulative_hazard_(T));
    }

    // 累积违约概率 PD(0, T) = 1 - Q(0, T)
    Real default_prob(Real T) const {
        return 1.0 - survival_prob(T);
    }

    // 区间违约概率 PD(t1, t2) = Q(0,t1) - Q(0,t2), 假设 t1 < t2
    Real default_prob(Real t1, Real t2) const {
        if (t2 <= t1) return 0.0;
        return survival_prob(t1) - survival_prob(t2);
    }

    // 瞬时 hazard rate h(T) (分段常数)
    Real hazard_rate(Real T) const {
        if (T <= 0.0) return hazard_rates_.empty() ? 0.0 : hazard_rates_[0];
        Size idx = find_segment_(T);
        return hazard_rates_[idx];
    }

    // 累积 hazard H(T) = ∫_0^T h(s) ds
    Real cumulative_hazard(Real T) const {
        return cumulative_hazard_(T);
    }

    // 从 CDS spread 单点近似构造 (flat hazard rate): h = s / (1 - R)
    static PDCurve from_cds_spread(Real cds_spread, Real recovery_rate,
                                     Real max_maturity = 30.0, Size n_points = 31) {
        if (recovery_rate >= 1.0) {
            throw std::invalid_argument("PDCurve::from_cds_spread: recovery must be < 1");
        }
        Real h = cds_spread / (1.0 - recovery_rate);
        return flat(h, max_maturity, n_points);
    }

    // 平坦 hazard rate 工厂
    static PDCurve flat(Real hazard_rate, Real max_maturity = 30.0, Size n_points = 31) {
        std::vector<Real> times, hrs;
        times.reserve(n_points);
        hrs.reserve(n_points);
        for (Size i = 0; i < n_points; ++i) {
            Real T = (i == 0) ? 0.01 : static_cast<Real>(i);
            if (T > max_maturity) break;
            times.push_back(T);
            hrs.push_back(hazard_rate);
        }
        return PDCurve(std::move(times), std::move(hrs));
    }

    const std::vector<Real>& times() const noexcept { return times_; }
    const std::vector<Real>& hazard_rates() const noexcept { return hazard_rates_; }
    bool empty() const noexcept { return times_.empty(); }
    Size size() const noexcept { return times_.size(); }

private:
    std::vector<Real> times_;          // 严格递增, times_[0] > 0
    std::vector<Real> hazard_rates_;   // hazard_rates_[i] 在 [times_[i], times_[i+1}) 内常数

    void validate() const {
        if (times_.size() != hazard_rates_.size()) {
            throw std::invalid_argument("PDCurve: times/hazard_rates size mismatch");
        }
        if (times_.empty()) {
            throw std::invalid_argument("PDCurve: at least 1 point required");
        }
        for (Size i = 0; i < times_.size(); ++i) {
            if (times_[i] <= 0.0) {
                throw std::invalid_argument("PDCurve: times must be positive");
            }
            if (hazard_rates_[i] < 0.0) {
                throw std::invalid_argument("PDCurve: hazard rates must be non-negative");
            }
            if (i > 0 && times_[i] <= times_[i - 1]) {
                throw std::invalid_argument("PDCurve: times must be strictly increasing");
            }
        }
    }

    // 累积 hazard H(T) = ∫_0^T h(s) ds, 分段常数 → 分段线性
    Real cumulative_hazard_(Real T) const {
        if (T <= 0.0) return 0.0;
        if (times_.size() == 1) {
            return hazard_rates_[0] * T;  // flat
        }
        if (T <= times_[0]) return hazard_rates_[0] * T;
        Real H = hazard_rates_[0] * times_[0];
        for (Size i = 0; i + 1 < times_.size(); ++i) {
            Real t_lo = times_[i];
            Real t_hi = times_[i + 1];
            Real h = hazard_rates_[i];
            if (T <= t_hi) {
                H += h * (T - t_lo);
                return H;
            }
            H += h * (t_hi - t_lo);
        }
        H += hazard_rates_.back() * (T - times_.back());
        return H;
    }

    Size find_segment_(Real T) const {
        if (T <= times_[0]) return 0;
        if (T >= times_.back()) return hazard_rates_.size() - 1;
        Size lo = 0, hi = times_.size() - 1;
        while (hi - lo > 1) {
            Size mid = (lo + hi) / 2;
            if (times_[mid] <= T) lo = mid;
            else hi = mid;
        }
        return lo;
    }
};

// ============ 信用曲线 (CreditCurve) ============
// 组合 PDCurve + recovery rate, 提供损失给定违约 (LGD) 与预期损失 (EL)
class CreditCurve {
public:
    CreditCurve() = default;

    CreditCurve(PDCurve pd_curve, Real recovery_rate)
        : pd_(std::move(pd_curve)), recovery_(recovery_rate) {
        if (recovery_ < 0.0 || recovery_ >= 1.0) {
            throw std::invalid_argument("CreditCurve: recovery must be in [0, 1)");
        }
    }

    // 从 flat CDS spread 构造 (h = s/(1-R) 单点近似)
    static CreditCurve from_flat_cds(Real cds_spread, Real recovery_rate,
                                       Real max_maturity = 30.0, Size n_points = 31) {
        return CreditCurve(
            PDCurve::from_cds_spread(cds_spread, recovery_rate, max_maturity, n_points),
            recovery_rate);
    }

    // 从 flat hazard rate 构造
    static CreditCurve flat(Real hazard_rate, Real recovery_rate,
                              Real max_maturity = 30.0, Size n_points = 31) {
        return CreditCurve(
            PDCurve::flat(hazard_rate, max_maturity, n_points),
            recovery_rate);
    }

    // 生存/违约概率 (委托给 PDCurve)
    Real survival_prob(Real T) const { return pd_.survival_prob(T); }
    Real default_prob(Real T) const { return pd_.default_prob(T); }
    Real default_prob(Real t1, Real t2) const { return pd_.default_prob(t1, t2); }
    Real hazard_rate(Real T) const { return pd_.hazard_rate(T); }
    Real cumulative_hazard(Real T) const { return pd_.cumulative_hazard(T); }

    // 回收率与损失给定违约
    Real recovery_rate() const noexcept { return recovery_; }
    Real lgd() const noexcept { return 1.0 - recovery_; }  // Loss Given Default

    // 预期损失 EL(0, T) = LGD * PD(0, T)
    Real expected_loss(Real T) const {
        return lgd() * default_prob(T);
    }

    // 区间预期损失 EL(t1, t2) = LGD * PD(t1, t2)
    Real expected_loss(Real t1, Real t2) const {
        return lgd() * default_prob(t1, t2);
    }

    const PDCurve& pd_curve() const noexcept { return pd_; }

private:
    PDCurve pd_;
    Real recovery_ = 0.40;  // 默认 40% 回收率
};

}  // namespace v1
}  // namespace cpphub
