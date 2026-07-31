#pragma once
// SOURCE: Li (2000) "On Default Correlation: A Copula Function Approach" J. Fixed Income 9(4)
// SOURCE: Hull & White (2004) "Valuation of a CDO and an n-th to Default CDS without Monte Carlo Simulation"
// SOURCE: Laurent & Gregory (2005) "Basket Default Swaps, CDOs and Factor Copulas"
// SOURCE: O'Kane (2008) "Modelling Single-Name and Multi-Name Credit Derivatives" Ch.11 (Basket CDS)
// SOURCE: Schönbucher (2003) "Credit Derivatives Pricing Models" Ch.10 (Multi-name)
// 模块: 篮子 CDS (Basket CDS / Nth-to-Default CDS)
//
// ==================== Basket CDS 结构 ====================
//
// Basket CDS 是基于多个参考实体 (reference entities) 的信用衍生品:
//   - Nth-to-Default CDS: 当第 N 个名字违约时, 保护卖方支付 (1-R) * N_tranche
//     之后合约终止 (或继续到下一违约, 取决于合约)
//   - M-th to N-th Default CDS: 覆盖第 M 到第 N 个违约的损失
//   - 全部违约 (All-to-Default): 所有名字都违约才赔付
//
// 本模块实现:
//   - NthToDefaultCDS: 第 N 个违约触发赔付, 之后合约终止
//   - MTMToDefaultCDS: M-th to N-th 违约篮子 (覆盖第 M 到第 N 个违约)
//
// ==================== 定价公式 (从保护买方视角) ====================
//
// 设 τ_(1) ≤ τ_(2) ≤ ... ≤ τ_(n) 为顺序统计量 (按违约时间排序)
// Nth-to-Default: 在 τ_(N) 时刻支付 (1-R) * N, 之后合约终止
//
// Premium leg (买方支付保费, 直到第 N 个违约或到期):
//   PV_premium = s * N * Σ_i τ_i * P(τ_(N) > t_i) * P(0, t_i)
//              + s * N * 0.5 * Σ_i τ_i * P(t_{i-1} < τ_(N) ≤ t_i) * P(0, t_mid_i)
//   第一项: 第 N 个违约在 t_i 之后 (合约存活), 支付完整保费
//   第二项: 第 N 个违约在 (t_{i-1}, t_i] 内, 支付半期应计
//
// Protection leg (买方获得赔付):
//   PV_protection = (1-R) * N * E[P(0, τ_(N))] * 1{τ_(N) ≤ T}
//                  = (1-R) * N * Σ_i P(τ_(N) ∈ (t_{i-1}, t_i]) * P(0, t_mid_i)
//
// Par spread:
//   s_par = PV_protection / Risky_pv01_basket
//   Risky_pv01_basket = Σ_i τ_i * [P(τ_(N) > t_i) * P(0,t_i)
//                                  + 0.5 * P(t_{i-1} < τ_(N) ≤ t_i) * P(0, t_mid_i)]
//
// ==================== 损失分布 (单因子 Gaussian Copula, 半解析) ====================
//
// 给定 M=m, 各名字条件独立:
//   PD_i(t | m) = Φ((Φ^{-1}(PD_i(t)) - √ρ * m) / √(1-ρ))
//   P(τ_i ≤ t | m) = PD_i(t | m), 独立
//
// 第 k 个顺序统计量条件分布 (异质组合):
//   P(τ_(k) ≤ t | m) = P(at least k of n default by t | m)
//   = Σ_{j=k}^{n} Σ_{S: |S|=j} Π_{i∈S} PD_i(t|m) * Π_{i∉S} (1 - PD_i(t|m))
//   (用 Poisson-binomial 分布)
//
// 同质组合 (homogeneous, 所有 PD_i 相同):
//   P(K ≥ k | m) = Σ_{j=k}^{n} C(n,j) * p(t|m)^j * (1-p(t|m))^{n-j}
//   其中 p(t|m) = PD(t|m), K|n,m ~ Binomial(n, p(t|m))
//
// 无条件: P(τ_(k) ≤ t) = ∫ P(τ_(k) ≤ t | m) φ(m) dm (Gauss-Hermite 积分)
//
// ==================== MC 模拟 (通用) ====================
//
// 适用于任意 Copula (Gaussian / t) 和任意异质组合:
//   1. 采样违约时间 τ_1, ..., τ_n
//   2. 排序得到 τ_(1) ≤ ... ≤ τ_(n)
//   3. 统计 τ_(N) 分布, 计算 premium / protection leg
//
// ==================== 约定 ====================
//
// - 输入: n 个 CreditCurve (允许异质), Copula (相关结构)
// - 第 N 个违约: N ∈ [1, n], N=1 为 first-to-default
// - 支付频率: 季度 (与 CDS 一致)
// - 违约时点: 区间中点 (mid-point convention)
// - 折现: OIS discount curve
// - 赔付: (1-R) * N_notional (假设所有名字 LGD 相同, 否则用各名字 LGD)

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include "cpphub/instruments/credit/credit_curve.hpp"
#include "cpphub/instruments/credit/copula.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <numeric>

namespace cpphub {
inline namespace v1 {

// ============ Basket CDS 配置 ============
struct BasketCDSConfig {
    Real notional = 1.0;            // 名义本金 (每个名字的本金)
    Real spread = 0.01;             // CDS spread (年化)
    Real maturity = 5.0;            // 到期时间 (年)
    Size n_premiums = 20;           // 保费支付次数
    Size nth_default = 1;           // 第 N 个违约 (1 = first-to-default)
    bool is_buyer = true;           // true = 保护买方视角
    std::vector<Real> payment_times;   // 自定义支付时间表
    std::vector<Real> year_fractions;  // 自定义年份分数

    void generate_schedule() {
        if (!payment_times.empty()) return;
        if (n_premiums == 0) {
            throw std::invalid_argument("BasketCDS: n_premiums must be positive");
        }
        Real tau = maturity / static_cast<Real>(n_premiums);
        payment_times.resize(n_premiums);
        year_fractions.assign(n_premiums, tau);
        for (Size i = 0; i < n_premiums; ++i) {
            payment_times[i] = static_cast<Real>(i + 1) * tau;
        }
    }

    void validate() const {
        if (notional <= 0.0) {
            throw std::invalid_argument("BasketCDS: notional must be positive");
        }
        if (spread < 0.0) {
            throw std::invalid_argument("BasketCDS: spread must be non-negative");
        }
        if (maturity <= 0.0) {
            throw std::invalid_argument("BasketCDS: maturity must be positive");
        }
        if (n_premiums == 0) {
            throw std::invalid_argument("BasketCDS: n_premiums must be positive");
        }
        if (nth_default == 0) {
            throw std::invalid_argument("BasketCDS: nth_default must be >= 1");
        }
    }
};

// ============ Basket CDS 结果 ============
struct BasketCDSResult {
    Real pv_premium_leg = 0.0;       // 保费腿 PV (买方视角, 正=支出)
    Real pv_protection_leg = 0.0;    // 保护腿 PV (买方视角, 正=收入)
    Real pv = 0.0;                   // Basket PV = PV_protection - PV_premium (买方视角)
    Real par_spread = 0.0;           // 让 PV = 0 的 spread
    Real risky_pv01 = 0.0;           // Risky annuity
    Real prob_nth_default_by_maturity = 0.0;  // P(τ_(N) ≤ T)
    std::vector<Real> survival_prob_at_t;  // P(τ_(N) > t_i) 各时点
    std::vector<Real> default_prob_at_t;   // P(τ_(N) ∈ (t_{i-1}, t_i]) 各区间
};

// ============ MC Basket CDS 定价器 ============
// 使用 Copula 模型采样违约时间, 统计第 N 个违约时间分布
// 适用于任意 Copula (Gaussian / t) 和任意异质组合
class BasketCDSMCPricer {
public:
    BasketCDSMCPricer(std::vector<CreditCurve> credit_curves,
                        const ZeroCurve& discount)
        : credit_curves_(std::move(credit_curves)),
          discount_(discount) {
        if (credit_curves_.empty()) {
            throw std::invalid_argument("BasketCDSMC: credit_curves empty");
        }
    }

    // MC 模拟 + basket CDS 定价
    // copula: 已配置好的 Copula
    // cfg: basket CDS 配置
    // n_paths: MC 路径数
    // seed: RNG 种子
    template <typename Copula>
    BasketCDSResult price(const Copula& copula,
                             const BasketCDSConfig& cfg,
                             Size n_paths, uint64_t seed = 42) const {
        BasketCDSConfig cfg_v = cfg;
        cfg_v.generate_schedule();
        cfg_v.validate();

        if (copula.n_names() != credit_curves_.size()) {
            throw std::invalid_argument("BasketCDSMC: copula.n_names != credit_curves size");
        }
        if (cfg_v.nth_default > credit_curves_.size()) {
            throw std::invalid_argument(
                "BasketCDSMC: nth_default > n_names (不可能违约)");
        }

        const Size n = cfg_v.n_premiums;
        const Size n_names = credit_curves_.size();
        const Real N_notional = cfg_v.notional;
        const Real s = cfg_v.spread;
        const Size Nth = cfg_v.nth_default;

        // 累积各时点统计量
        std::vector<Real> sum_survival(n, 0.0);  // 1{τ_(N) > t_i}
        std::vector<Real> sum_default_in_interval(n, 0.0);  // 1{t_{i-1} < τ_(N) ≤ t_i}
        Real sum_premium = 0.0;
        Real sum_protection = 0.0;
        Real sum_rpv01 = 0.0;
        Real sum_default_by_T = 0.0;  // 1{τ_(N) ≤ T}

        // 计算 LGD (假设所有名字 LGD 相同, 取第一个; 否则用加权平均)
        Real lgd = 0.0;
        for (Size i = 0; i < n_names; ++i) lgd += credit_curves_[i].lgd();
        lgd /= static_cast<Real>(n_names);

        Philox4x64 rng(seed);

        for (Size path = 0; path < n_paths; ++path) {
            // 采样违约时间
            auto tau = copula.sample_default_times(rng, credit_curves_);

            // 排序得到顺序统计量
            std::vector<Real> sorted_tau(tau.begin(), tau.end());
            std::sort(sorted_tau.begin(), sorted_tau.end());

            // 第 N 个违约时间
            Real tau_N = sorted_tau[Nth - 1];  // 0-indexed → Nth-1

            // 检查是否在到期前违约
            bool default_before_T = (tau_N <= cfg_v.maturity) && !std::isinf(tau_N);

            // 遍历各支付时点
            for (Size i = 0; i < n; ++i) {
                Real t_i = cfg_v.payment_times[i];
                Real t_prev = (i == 0) ? 0.0 : cfg_v.payment_times[i - 1];
                Real tau_i = cfg_v.year_fractions[i];
                Real t_mid = 0.5 * (t_prev + t_i);
                Real P_i = discount_.discount_factor(t_i);
                Real P_mid = discount_.discount_factor(t_mid);

                // 1{τ_(N) > t_i}: 合约在 t_i 存活
                bool alive_at_t_i = (tau_N > t_i);
                // 1{t_{i-1} < τ_(N) ≤ t_i}: 第 N 个违约在此区间
                bool default_in_interval = (tau_N > t_prev) && (tau_N <= t_i);

                if (alive_at_t_i) sum_survival[i] += 1.0;
                if (default_in_interval) sum_default_in_interval[i] += 1.0;

                // Premium leg: 完整支付 + 半期应计
                Real cf_premium = 0.0;
                if (alive_at_t_i) {
                    cf_premium += s * N_notional * tau_i * P_i;
                } else if (default_in_interval) {
                    cf_premium += s * N_notional * 0.5 * tau_i * P_mid;
                }
                sum_premium += cf_premium;
                sum_rpv01 += (alive_at_t_i)
                    ? tau_i * N_notional * P_i
                    : (default_in_interval ? 0.5 * tau_i * N_notional * P_mid : 0.0);

                // Protection leg: 区间内违约赔付 (折现到 t_mid)
                if (default_in_interval) {
                    sum_protection += lgd * N_notional * P_mid;
                }
            }

            if (default_before_T) sum_default_by_T += 1.0;
        }

        // 取平均
        Real inv_n = 1.0 / static_cast<Real>(n_paths);
        BasketCDSResult result;
        result.survival_prob_at_t.resize(n);
        result.default_prob_at_t.resize(n);
        for (Size i = 0; i < n; ++i) {
            result.survival_prob_at_t[i] = sum_survival[i] * inv_n;
            result.default_prob_at_t[i] = sum_default_in_interval[i] * inv_n;
        }
        result.pv_premium_leg = sum_premium * inv_n;
        result.pv_protection_leg = sum_protection * inv_n;
        result.risky_pv01 = sum_rpv01 * inv_n;
        result.prob_nth_default_by_maturity = sum_default_by_T * inv_n;
        result.pv = result.pv_protection_leg - result.pv_premium_leg;
        result.par_spread = (result.risky_pv01 > 0.0)
                              ? result.pv_protection_leg / result.risky_pv01 : 0.0;

        if (!cfg_v.is_buyer) {
            result.pv = -result.pv;
            result.pv_premium_leg = -result.pv_premium_leg;
            result.pv_protection_leg = -result.pv_protection_leg;
        }

        return result;
    }

    Size n_names() const noexcept { return credit_curves_.size(); }

private:
    std::vector<CreditCurve> credit_curves_;
    const ZeroCurve& discount_;
};

// ============ 同质组合 + 单因子 Gaussian Copula 半解析定价器 ============
// 适用于所有名字同质 (相同 PD, LGD) + 单因子相关结构
// 用 Gauss-Hermite 积分计算 P(τ_(N) ≤ t)
class NthToDefaultHomogeneousPricer {
public:
    NthToDefaultHomogeneousPricer(Real rho, Real pd_at_maturity, Real lgd,
                                     Size n_names, const ZeroCurve& discount,
                                     Real maturity)
        : rho_(rho), pd_(pd_at_maturity), lgd_(lgd),
          n_names_(n_names), discount_(discount), maturity_(maturity) {
        if (rho_ < 0.0 || rho_ >= 1.0) {
            throw std::invalid_argument("NthToDefault: rho must be in [0, 1)");
        }
        if (pd_ < 0.0 || pd_ > 1.0) {
            throw std::invalid_argument("NthToDefault: pd must be in [0, 1]");
        }
        if (lgd_ < 0.0 || lgd_ > 1.0) {
            throw std::invalid_argument("NthToDefault: lgd must be in [0, 1]");
        }
        if (n_names_ == 0) {
            throw std::invalid_argument("NthToDefault: n_names must be positive");
        }
        sqrt_rho_ = std::sqrt(rho_);
        sqrt_1mrho_ = std::sqrt(1.0 - rho_);
    }

    // 计算 P(τ_(N) ≤ t), 即第 N 个违约在 t 前发生的概率
    // = E_M[P(K ≥ N | M)], K|M ~ Binomial(n_names, PD(t|M))
    // = ∫ Σ_{k=N}^{n} C(n,k) * p(t|m)^k * (1-p(t|m))^{n-k} φ(m) dm
    Real prob_nth_default_by(Real t, Size N) const {
        if (N > n_names_) return 0.0;  // 不可能
        if (t <= 0.0) return 0.0;
        Real pd_t = pd_at_time_(t);
        if (pd_t <= 0.0) return 0.0;

        auto integrand = [&](Real m) -> Real {
            Real p_cond = conditional_pd_(pd_t, m);
            // P(K ≥ N | m) = Σ_{k=N}^{n} C(n,k) p^k (1-p)^{n-k}
            return binomial_tail_(n_names_, N, p_cond);
        };

        return gauss_hermite_integrate_(integrand);
    }

    // 计算 P(τ_(N) > t) = 1 - P(τ_(N) ≤ t)
    Real survival_prob(Real t, Size N) const {
        return 1.0 - prob_nth_default_by(t, N);
    }

    // 完整定价
    BasketCDSResult price(const BasketCDSConfig& cfg) const {
        BasketCDSConfig cfg_v = cfg;
        cfg_v.generate_schedule();
        cfg_v.validate();

        const Size N = cfg_v.nth_default;
        if (N > n_names_) {
            throw std::invalid_argument(
                "NthToDefault: nth_default > n_names");
        }

        const Size n = cfg_v.n_premiums;
        const Real N_notional = cfg_v.notional;
        const Real s = cfg_v.spread;

        BasketCDSResult result;
        result.survival_prob_at_t.resize(n);
        result.default_prob_at_t.resize(n);

        Real pv_premium = 0.0;
        Real pv_protection = 0.0;
        Real rpv01 = 0.0;
        Real prev_cum_default = 0.0;  // P(τ_(N) ≤ t_{i-1})

        for (Size i = 0; i < n; ++i) {
            Real t_i = cfg_v.payment_times[i];
            Real t_prev = (i == 0) ? 0.0 : cfg_v.payment_times[i - 1];
            Real tau_i = cfg_v.year_fractions[i];
            Real t_mid = 0.5 * (t_prev + t_i);
            Real P_i = discount_.discount_factor(t_i);
            Real P_mid = discount_.discount_factor(t_mid);

            Real cum_default_i = prob_nth_default_by(t_i, N);
            Real survival_i = 1.0 - cum_default_i;
            Real default_in_interval = std::max(0.0, cum_default_i - prev_cum_default);

            result.survival_prob_at_t[i] = survival_i;
            result.default_prob_at_t[i] = default_in_interval;

            // Premium leg: 完整支付 + 半期应计
            pv_premium += s * N_notional * tau_i * survival_i * P_i;
            pv_premium += s * N_notional * 0.5 * tau_i * default_in_interval * P_mid;
            rpv01 += tau_i * survival_i * P_i + 0.5 * tau_i * default_in_interval * P_mid;

            // Protection leg
            pv_protection += lgd_ * N_notional * default_in_interval * P_mid;

            prev_cum_default = cum_default_i;
        }

        result.pv_premium_leg = pv_premium * N_notional / N_notional;  // 保持单位
        result.pv_protection_leg = pv_protection;
        result.risky_pv01 = rpv01 * N_notional;  // unit: N_notional * time
        result.prob_nth_default_by_maturity = prev_cum_default;
        result.pv = pv_protection - pv_premium;
        result.par_spread = (result.risky_pv01 > 0.0)
                              ? result.pv_protection_leg / result.risky_pv01 : 0.0;

        if (!cfg_v.is_buyer) {
            result.pv = -result.pv;
            result.pv_premium_leg = -result.pv_premium_leg;
            result.pv_protection_leg = -result.pv_protection_leg;
        }

        return result;
    }

    Real rho() const noexcept { return rho_; }
    Size n_names() const noexcept { return n_names_; }
    Real pd_at_maturity() const noexcept { return pd_; }
    Real lgd() const noexcept { return lgd_; }

private:
    Real rho_;
    Real pd_;
    Real lgd_;
    Size n_names_;
    const ZeroCurve& discount_;
    Real maturity_;
    Real sqrt_rho_;
    Real sqrt_1mrho_;

    // PD(t) 假设 flat hazard
    Real pd_at_time_(Real t) const {
        if (t <= 0.0) return 0.0;
        if (t >= maturity_) return pd_;
        if (pd_ >= 1.0) return 1.0;
        Real h = -std::log(1.0 - pd_) / maturity_;
        return 1.0 - std::exp(-h * t);
    }

    Real conditional_pd_(Real pd_t, Real m) const {
        if (pd_t <= 0.0) return 0.0;
        if (pd_t >= 1.0) return 1.0;
        Real z = inv_normal_cdf(pd_t);
        Real x = (z - sqrt_rho_ * m) / sqrt_1mrho_;
        return normal_cdf(x);
    }

    // 二项分布尾概率 P(K ≥ k | n, p) = Σ_{j=k}^{n} C(n,j) p^j (1-p)^{n-j}
    static Real binomial_tail_(Size n, Size k, Real p) {
        if (k > n) return 0.0;
        if (p <= 0.0) return 0.0;
        if (p >= 1.0) return 1.0;
        // 用对数计算避免溢出
        Real q = 1.0 - p;
        Real log_p = std::log(p);
        Real log_q = std::log(q);
        Real sum = 0.0;
        for (Size j = k; j <= n; ++j) {
            Real log_term = log_binomial_(n, j) + static_cast<Real>(j) * log_p
                            + static_cast<Real>(n - j) * log_q;
            sum += std::exp(log_term);
        }
        return std::min(sum, 1.0);
    }

    // log(C(n, k))
    static Real log_binomial_(Size n, Size k) {
        if (k > n) return -std::numeric_limits<Real>::infinity();
        if (k == 0 || k == n) return 0.0;
        // C(n,k) = n! / (k! (n-k)!)
        // log C(n,k) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
        return std::lgamma(static_cast<Real>(n + 1))
               - std::lgamma(static_cast<Real>(k + 1))
               - std::lgamma(static_cast<Real>(n - k + 1));
    }

    // Gauss-Hermite 积分: ∫_{-∞}^{∞} f(m) φ(m) dm = (1/√π) Σ_{正+负} w_k f(√2 * x_k)
    // 使用 n=20 GH (e^{-x²} 权重) 的 10 个正节点, 镜像得 20 节点
    // 节点/权重来源: numpy.polynomial.hermite.hermgauss(20)
    // 验证: Σ_{正} w_k = √π/2 ≈ 0.8862, 完整 Σ w_k = √π
    template <typename F>
    Real gauss_hermite_integrate_(F&& f) const {
        // n=20 Gauss-Hermite 正节点 (physicist's, e^{-x²} 权重)
        static const Real nodes[] = {
            0.2453407083009, 0.7374737285454, 1.2340762153953, 1.7385377121166,
            2.2549740020893, 2.7888060584281, 3.3478545673832, 3.9447640401156,
            4.6036824495507, 5.3874808900112
        };
        static const Real weights[] = {
            0.4622436696006, 0.2866755053628, 0.1090172060200, 0.0248105208875,
            0.0032437733422, 0.0002283386360, 0.0000078025565, 0.0000001086069,
            0.0000000004399, 0.000000000000223
        };
        const Size n_nodes = 10;

        Real sum = 0.0;
        for (Size i = 0; i < n_nodes; ++i) {
            Real m_pos = SQRT_2 * nodes[i];
            Real m_neg = -m_pos;
            sum += weights[i] * (f(m_pos) + f(m_neg));
        }
        static const Real sqrt_pi = std::sqrt(PI);
        return sum / sqrt_pi;
    }
};

// ============ 便捷工厂 ============

// 构造 N 个同质信用曲线 (相同 hazard rate, recovery)
inline std::vector<CreditCurve> make_homogeneous_credit_curves(
        Size n_names, Real hazard_rate, Real recovery_rate,
        Real max_maturity = 30.0, Size n_points = 31) {
    std::vector<CreditCurve> curves;
    curves.reserve(n_names);
    for (Size i = 0; i < n_names; ++i) {
        curves.push_back(CreditCurve::flat(hazard_rate, recovery_rate,
                                              max_maturity, n_points));
    }
    return curves;
}

// 构造 N 个异质信用曲线 (从 hazard rate 向量)
inline std::vector<CreditCurve> make_heterogeneous_credit_curves(
        const std::vector<Real>& hazard_rates, Real recovery_rate,
        Real max_maturity = 30.0, Size n_points = 31) {
    std::vector<CreditCurve> curves;
    curves.reserve(hazard_rates.size());
    for (Real h : hazard_rates) {
        curves.push_back(CreditCurve::flat(h, recovery_rate,
                                              max_maturity, n_points));
    }
    return curves;
}

}  // namespace v1
}  // namespace cpphub
