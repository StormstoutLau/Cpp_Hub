#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

namespace detail {

/// χ²(k) 的 CDF: P(χ²_k ≤ x) = P(k/2, x/2) (正则化下不完全 Gamma, NR §6.2)
/// RISK-028: 原 Wilson-Hilferty 三次根近似在 k=1 临界带偏差 ~3e-3 (可翻转结论),
/// 改用 core/math.hpp 的精确实现 (精度 ~1e-14)。
/// 注意: 上尾按 1 - P 计算, 故 p < ~1e-16 时归零 (不影响拒绝/不拒绝判定)。
inline Real chi2_cdf(Real x, Size k) {
    if (x <= 0) return 0;  // 数值上略负的 LR 亦归零, 不抛异常
    return regularized_lower_gamma(static_cast<Real>(k) / Real(2), x / Real(2));
}

/// x·log(y), 约定 0·log(0) = 0 (RISK-024: 计数为 0 的似然项贡献 0 而非 NaN)
inline Real xlogy(Real x, Real y) {
    return (x == Real(0)) ? Real(0) : x * std::log(y);
}

/// x·log(1-y), 约定 0·log(0) = 0; y ∈ [0, 1]
inline Real xlog1m(Real x, Real y) {
    return (x == Real(0)) ? Real(0) : x * std::log1p(-y);
}

/// 二项对数似然 logL = (N - x)·log(1-p) + x·log(p)
/// RISK-024: 逐项对数化替代连乘 (原形式在 N ≳ 3850 时下溢为 0)
inline Real log_lik_binomial(Size n_obs, Size n_success, Real p_success) {
    return xlog1m(static_cast<Real>(n_obs - n_success), p_success)
         + xlogy(static_cast<Real>(n_success), p_success);
}

/// 校验 confidence ∈ (0, 1) 并返回期望违约率 p = 1 - confidence
/// (confidence ≥ 1 会让 log(p) 取 -inf, 产生静默 inf/NaN)
inline Real expected_rate_from_confidence(Real confidence, const char* who) {
    if (!(confidence > Real(0) && confidence < Real(1))) {
        throw std::invalid_argument(std::string(who) + ": confidence 必须 ∈ (0, 1)");
    }
    return Real(1) - confidence;
}

} // namespace detail

struct BacktestResult {
    Size n_violations;
    Size n_observations;
    Real violation_rate;
    Real expected_rate;
    Real p_value;
    bool reject_null;
};

class KupiecPOF {
public:
    /// Kupiec POF 检验 (无条件覆盖率), LR ~ χ²(1)
    /// @throws std::invalid_argument n_observations == 0 / n_violations > n_observations / confidence ∉ (0,1)
    static BacktestResult test(Size n_violations, Size n_observations, Real confidence) {
        if (n_observations == 0 || n_violations > n_observations) {
            throw std::invalid_argument(
                "KupiecPOF::test: 需要 0 <= n_violations <= n_observations 且 n_observations > 0");
        }
        Real p = detail::expected_rate_from_confidence(confidence, "KupiecPOF::test");
        Real violation_rate = static_cast<Real>(n_violations) / static_cast<Real>(n_observations);

        // 对数空间 LR = 2·(logL_alt - logL_null) (RISK-024)
        // 0·log0 = 0 约定天然覆盖 v=0 与 v=N 两个边界, 原特判分支由此删除
        Real logL0 = detail::log_lik_binomial(n_observations, n_violations, p);
        Real logL1 = detail::log_lik_binomial(n_observations, n_violations, violation_rate);
        Real LR = Real(2) * (logL1 - logL0);
        Real p_value = Real(1) - detail::chi2_cdf(LR, Size(1));

        return {n_violations, n_observations, violation_rate, p, p_value, p_value < Real(0.05)};
    }

    static BacktestResult test(const std::vector<Real>& var_series,
                                const std::vector<Real>& realized_losses,
                                Real confidence) {
        Size n_violations = 0;
        Size n = std::min(var_series.size(), realized_losses.size());
        for (Size i = 0; i < n; ++i) {
            if (realized_losses[i] > var_series[i]) {
                ++n_violations;
            }
        }
        return test(n_violations, n, confidence);
    }
};

class ChristoffersenIID {
public:
    /// 独立性检验 (Christoffersen 1998 §3), LR_ind ~ χ²(1)
    /// @throws std::invalid_argument 观测数 < 2 (独立性基于 n-1 个转移)
    static BacktestResult test(const std::vector<Real>& var_series,
                                const std::vector<Real>& realized_losses) {
        TransitionCounts t = count_transitions(var_series, realized_losses);
        Real pi = static_cast<Real>(t.v_lag) / static_cast<Real>(t.n - 1);
        Real LR_ind = independence_lr(t);

        Real p_value = Real(1) - detail::chi2_cdf(LR_ind, Size(1));
        return {t.v, t.n, static_cast<Real>(t.v) / static_cast<Real>(t.n),
                pi, p_value, p_value < Real(0.05)};
    }

    /// 联合条件覆盖检验, LR_cc = LR_pof + LR_ind ~ χ²(2) (Christoffersen 1998 式 (10))
    ///
    /// RISK-023: 原实现 LR_cc_total = LR_cc 只含覆盖率分量, 独立性证据被丢弃
    /// (且按 χ²(2) 取值 ⇒ 检验过于保守)。此处按定义把两个对数似然比相加。
    /// @throws std::invalid_argument 观测数 < 2 / confidence ∉ (0,1)
    static BacktestResult joint_test(const std::vector<Real>& var_series,
                                      const std::vector<Real>& realized_losses,
                                      Real confidence) {
        Real p = detail::expected_rate_from_confidence(confidence, "ChristoffersenIID::joint_test");
        TransitionCounts t = count_transitions(var_series, realized_losses);

        // LR_pof: H0 = 期望违约率 p (全样本 n 个观测)
        Real logL0 = detail::log_lik_binomial(t.n, t.v, p);
        Real logL1 = detail::log_lik_binomial(t.n, t.v, static_cast<Real>(t.v) / static_cast<Real>(t.n));
        Real LR_pof = Real(2) * (logL1 - logL0);

        // LR_ind: 与 test() 共用同一实现 (单一来源)
        Real LR_cc = LR_pof + independence_lr(t);

        Real p_value = Real(1) - detail::chi2_cdf(LR_cc, Size(2));
        return {t.v, t.n, static_cast<Real>(t.v) / static_cast<Real>(t.n),
                p, p_value, p_value < Real(0.05)};
    }

private:
    /// 违约指示与 Christoffersen 转移计数 (1998 §3)
    struct TransitionCounts {
        Size n = 0;       ///< 观测数 (含首期)
        Size v = 0;       ///< 全样本违约数 (与 Kupiec/联合检验同口径)
        Size v_lag = 0;   ///< 独立性样本内的违约数 = n01 + n11
        Size n00 = 0, n01 = 0, n10 = 0, n11 = 0;

        Size n0() const { return n00 + n01; }
        Size n1() const { return n10 + n11; }
    };

    static TransitionCounts count_transitions(const std::vector<Real>& var_series,
                                              const std::vector<Real>& realized_losses) {
        TransitionCounts t;
        t.n = std::min(var_series.size(), realized_losses.size());
        if (t.n < 2) {
            throw std::invalid_argument(
                "ChristoffersenIID: 至少需要 2 个观测 (独立性检验基于 n-1 个转移)");
        }
        bool prev = realized_losses[0] > var_series[0];
        if (prev) ++t.v;
        for (Size i = 1; i < t.n; ++i) {
            bool cur = realized_losses[i] > var_series[i];
            if (cur) ++t.v;
            if (!prev && !cur)      ++t.n00;
            else if (!prev &&  cur) ++t.n01;
            else if ( prev && !cur) ++t.n10;
            else                    ++t.n11;
            prev = cur;
        }
        t.v_lag = t.n01 + t.n11;
        return t;
    }

    /// 独立性对数似然比 (单一实现, test 与 joint_test 共用)
    ///
    /// RISK-024: 原连乘形式在 n ≳ 3750 时下溢 (L_IID → 0 / 次正规), 旧守卫
    /// (L_IID>0 && L_nonIID>0) 会命中并把 LR_ind 置 0 → p=1 (静默不拒绝聚簇)。
    /// 对数空间 + 0·log0 = 0 约定后, "完美聚集" 正确给出 logL_non = 0 → LR_ind > 0;
    /// 守卫不再需要, 已删除。
    static Real independence_lr(const TransitionCounts& t) {
        Real pi  = static_cast<Real>(t.v_lag) / static_cast<Real>(t.n - 1);
        Real pi0 = (t.n0() > 0) ? static_cast<Real>(t.n01) / static_cast<Real>(t.n0()) : Real(0);
        Real pi1 = (t.n1() > 0) ? static_cast<Real>(t.n11) / static_cast<Real>(t.n1()) : Real(0);

        Real logL_iid = detail::log_lik_binomial(t.n - 1, t.v_lag, pi);
        Real logL_non = detail::xlog1m(static_cast<Real>(t.n00), pi0)
                      + detail::xlogy(static_cast<Real>(t.n01), pi0)
                      + detail::xlog1m(static_cast<Real>(t.n10), pi1)
                      + detail::xlogy(static_cast<Real>(t.n11), pi1);

        return Real(2) * (logL_non - logL_iid);
    }
};

enum class BaselZone { Green, Yellow, Red };

struct BaselTrafficLightResult {
    BaselZone zone;
    Size n_violations;
    Size n_observations;
    Real capital_multiplier;
    std::string description;
};

class BaselTrafficLight {
public:
    static BaselTrafficLightResult assess(Size n_violations, Size n_observations = 250) {
        BaselZone zone;
        Real multiplier;
        std::string desc;

        if (n_violations <= 4) {
            zone = BaselZone::Green;
            multiplier = Real(3.0);
            desc = "Green zone: model OK";
        } else if (n_violations <= 9) {
            zone = BaselZone::Yellow;
            static const Real multipliers[] = {Real(3.4), Real(3.5), Real(3.65), Real(3.75), Real(3.85)};
            multiplier = multipliers[n_violations - 5];
            desc = "Yellow zone: model may be inaccurate";
        } else {
            zone = BaselZone::Red;
            multiplier = Real(4.0);
            desc = "Red zone: model is inaccurate";
        }

        return {zone, n_violations, n_observations, multiplier, desc};
    }
};

} // namespace v1
} // namespace cpphub
