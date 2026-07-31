#pragma once
// SOURCE: Li (2000) "On Default Correlation: A Copula Function Approach" J. Fixed Income 9(4), 43-54
// SOURCE: Laurent & Gregory (2005) "Basket Default Swaps, CDOs and Factor Copulas"
// SOURCE: Andersen, Sidenius & Basu (2003) "All Your Hedges in One Basket" Risk 16(11), 67-72
// SOURCE: Hull & White (2004) "Valuation of a CDO and an n-th to Default CDS without Monte Carlo Simulation"
// SOURCE: Vasicek (1987) "Probability of Loss on Loan Portfolio" KMV technical report
// SOURCE: O'Kane (2008) "Modelling Single-Name and Multi-Name Credit Derivatives" Ch.12-13 (CDO)
// 模块: 合成 CDO (Synthetic CDO) 定价 — Tranche 损失分布 + 保费/保护腿
//
// ==================== CDO 结构 ====================
//
// 合成 CDO: 将一组 CDS (reference portfolio) 的信用风险重新分配到不同 tranches.
// 每个 tranche 由 attachment point (A) 和 detachment point (D) 定义:
//   - Tranche 损失 L(t) = max(0, min(L_port(t), D) - A)
//     其中 L_port(t) 为组合累计损失
//   - Tranche 名义本金 N(t) = (D - A) - L(t) = max(0, D - max(A, L_port(t)))
//
// 标准化 CDX/iTraxx tranches (示例):
//   CDX IG: 0-3%, 3-7%, 7-10%, 10-15%, 15-30%, 30-100%
//   iTraxx Europe: 0-3%, 3-6%, 6-9%, 9-12%, 12-22%, 22-100%
//
// ==================== 定价公式 (从保护卖方视角) ====================
//
// Premium leg (卖方收取保费):
//   PV_premium = s * N_tranche * Σ_i τ_i * E[(D-A-L(t_i))^+] * P(0, t_i)
//              + s * N_tranche * 0.5 * Σ_i τ_i * E[L_tranche(t_{i-1}) - L_tranche(t_i)] * P(0, t_mid_i)
//   第一项: 生存到 t_i 支付完整保费 (基于剩余名义)
//   第二项: 区间内违约的应计 (基于损失变化)
//
// Protection leg (卖方支付赔付):
//   PV_protection = E[Σ_i (L_tranche(t_i) - L_tranche(t_{i-1})) * P(0, t_mid_i)]
//                 = Σ_i E[L_tranche(t_i) - L_tranche(t_{i-1})] * P(0, t_mid_i)
//   即期望损失增量的折现值, mid-point convention
//
// Par spread (让 PV = 0 的 s):
//   s_par = PV_protection / Risky_pv01_tranche
//   Risky_pv01_tranche = Σ_i τ_i * [E[(D-A-L(t_i))^+] * P(0,t_i)
//                                  + 0.5 * E[ΔL_tranche(t_i)] * P(0, t_mid_i)]
//
// ==================== 损失分布计算 (单因子 Gaussian Copula) ====================
//
// 组合损失 L(t) = Σ_i LGD_i * 1{τ_i ≤ t}
// 在单因子模型下, 给定 M=m, 各名字条件独立:
//   PD_i(t | m) = Φ((Φ^{-1}(PD_i(t)) - √ρ * m) / √(1-ρ))
//   L(t | m) = Σ_i LGD_i * Bernoulli(PD_i(t | m))
//
// 等权重组合 (homogeneous): LGD_i = LGD, PD_i(t) = PD(t)
//   条件违约数 K | M=m ~ Binomial(n, PD(t|m))
//   条件损失 L(t|m) = K * LGD
//   无条件 E[L_tranche(t)] = ∫ L_tranche(t|m) φ(m) dm (用 Gaussian quadrature 数值积分)
//
// 异质组合 (heterogeneous): 用 Andersen-Sidenius-Basu (ASB) 递归
//   或 large-homogeneous-portfolio (LHA) Vasicek 渐近近似
//
// ==================== Large Homogeneous Portfolio (LHP/Vasicek) 近似 ====================
//
// 当 n → ∞ 且组合同质, 由大数定律:
//   L(t|m)/N → LGD * PD(t|m)  (a.s.)
//   即条件损失率 = LGD * Φ((Φ^{-1}(PD(t)) - √ρ * m) / √(1-ρ))
// Tranche 条件损失率:
//   l_tranche(t|m) = max(0, min(L(t|m)/N, D) - A)
// 无条件 E[L_tranche(t)/N] = ∫ l_tranche(t|m) φ(m) dm
//
// 本模块实现:
//   1. LHP 近似 (半解析, 快速, 适合同质大组合)
//   2. MC 模拟 (精确, 处理任意异质组合)
//
// ==================== 约定 ====================
//
// - Tranche: attachment/detachment 以 fraction 表示 (如 0.03, 0.07)
// - 组合损失: 以 fraction of total notional 表示
// - 支付频率: 季度 (与 CDS 一致)
// - 违约时点: 区间中点 (mid-point convention)
// - 折现: OIS discount curve
// - 信用曲线: 各名字的 CreditCurve (允许异质)

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
#include <functional>
#include <limits>

namespace cpphub {
inline namespace v1 {

// ============ CDO Tranche 配置 ============
struct CDOTrancheConfig {
    Real attachment = 0.0;          // A: attachment point (fraction, e.g. 0.03)
    Real detachment = 0.03;         // D: detachment point (fraction, e.g. 0.07)
    Real spread = 0.05;             // Tranche spread (年化, e.g. 0.05 = 500bp for equity)
    Real maturity = 5.0;            // 到期时间 (年)
    Size n_premiums = 20;           // 保费支付次数 (5Y 季度 = 20)
    bool is_protection_seller = true;  // true = 卖方视角 (收保费, 付赔付)
    std::vector<Real> payment_times;   // 自定义支付时间表 (空则等间隔)
    std::vector<Real> year_fractions;  // 自定义年份分数 (空则等间隔)

    void generate_schedule() {
        if (!payment_times.empty()) return;
        if (n_premiums == 0) {
            throw std::invalid_argument("CDOTranche: n_premiums must be positive");
        }
        Real tau = maturity / static_cast<Real>(n_premiums);
        payment_times.resize(n_premiums);
        year_fractions.assign(n_premiums, tau);
        for (Size i = 0; i < n_premiums; ++i) {
            payment_times[i] = static_cast<Real>(i + 1) * tau;
        }
    }

    void validate() const {
        if (attachment < 0.0 || attachment >= 1.0) {
            throw std::invalid_argument("CDOTranche: attachment must be in [0, 1)");
        }
        if (detachment <= attachment || detachment > 1.0) {
            throw std::invalid_argument("CDOTranche: detachment must be in (attachment, 1]");
        }
        if (spread < 0.0) {
            throw std::invalid_argument("CDOTranche: spread must be non-negative");
        }
        if (maturity <= 0.0) {
            throw std::invalid_argument("CDOTranche: maturity must be positive");
        }
        if (n_premiums == 0) {
            throw std::invalid_argument("CDOTranche: n_premiums must be positive");
        }
        if (!payment_times.empty() && payment_times.size() != n_premiums) {
            throw std::invalid_argument("CDOTranche: payment_times size mismatch");
        }
        if (!year_fractions.empty() && year_fractions.size() != n_premiums) {
            throw std::invalid_argument("CDOTranche: year_fractions size mismatch");
        }
    }

    Real tranche_width() const noexcept { return detachment - attachment; }
};

// ============ CDO Tranche 结果 ============
struct CDOTrancheResult {
    Real pv_premium_leg = 0.0;       // 保费腿 PV (卖方视角, 正=收入)
    Real pv_protection_leg = 0.0;    // 保护腿 PV (卖方视角, 正=支出)
    Real pv = 0.0;                   // Tranche PV = PV_premium - PV_protection (卖方视角)
    Real par_spread = 0.0;           // 让 PV = 0 的 spread
    Real risky_pv01 = 0.0;           // Tranche risky annuity (unit spread PV)
    Real expected_loss = 0.0;        // E[L_tranche(T)] (Tranche 总期望损失 fraction)
    std::vector<Real> expected_tranche_loss;  // 各支付时点 E[L_tranche(t_i)]
    std::vector<Real> expected_remaining_notional;  // E[(D-A-L(t_i))^+]
};

// ============ LHP (Large Homogeneous Portfolio) 定价器 ============
// 基于 Vasicek (1987) 单因子模型大组合渐近近似
// 假设: 所有名字同质 (相同 PD, LGD), 数量足够大
// 优点: 半解析, 计算快速, 适合标准 CDX/iTraxx 定价
// 局限: 不能处理异质组合 (不同信用曲线)
class CDOLHPPricer {
public:
    CDOLHPPricer(Real rho, Real pd_maturity, Real lgd,
                   const ZeroCurve& discount, Real maturity)
        : rho_(rho), pd_(pd_maturity), lgd_(lgd),
          discount_(discount), maturity_(maturity) {
        if (rho_ < 0.0 || rho_ >= 1.0) {
            throw std::invalid_argument("CDOLHP: rho must be in [0, 1)");
        }
        if (pd_ < 0.0 || pd_ > 1.0) {
            throw std::invalid_argument("CDOLHP: pd must be in [0, 1]");
        }
        if (lgd_ < 0.0 || lgd_ > 1.0) {
            throw std::invalid_argument("CDOLHP: lgd must be in [0, 1]");
        }
        if (maturity_ <= 0.0) {
            throw std::invalid_argument("CDOLHP: maturity must be positive");
        }
        sqrt_rho_ = std::sqrt(rho_);
        sqrt_1mrho_ = std::sqrt(1.0 - rho_);
    }

    // 计算 t 时刻组合损失率分布的期望损失率 E[L(t)/N]
    // E[L(t)/N] = LGD * PD(t)
    // 注: 这是无条件期望, 与 rho 无关 (线性可加性)
    Real expected_portfolio_loss(Real t) const {
        // 若 pd_ 对应 maturity_, 则按比例缩放到 t
        // 简化: 假设 flat hazard, PD(t) = 1 - exp(-h*t), h = -ln(1-pd_)/maturity_
        Real pd_t = pd_at_time_(t);
        return lgd_ * pd_t;
    }

    // 计算 t 时刻 tranche 期望损失率 E[L_tranche(t) / (D-A)]
    // LHP: L(t|m)/N → LGD * PD(t|m) (a.s.)
    //   PD(t|m) = Φ((Φ^{-1}(PD(t)) - √ρ * m) / √(1-ρ))
    //   L_tranche(t|m) = max(0, min(L(t|m)/N, D) - A)
    //   E[L_tranche(t)/N] = ∫ max(0, min(lgd*PD(t|m), D) - A) φ(m) dm
    // 用 Gaussian quadrature 数值积分 (Hermite 节点)
    Real expected_tranche_loss(Real t, Real A, Real D) const {
        if (t <= 0.0) return 0.0;
        Real pd_t = pd_at_time_(t);
        if (pd_t <= 0.0) return 0.0;

        // 被积函数 f(m) = max(0, min(lgd * PD(t|m), D) - A)
        auto integrand = [&](Real m) -> Real {
            Real pd_cond = conditional_pd_lhp_(pd_t, m);
            Real loss_rate = lgd_ * pd_cond;  // L(t|m)/N
            Real tranche_loss = std::max(0.0, std::min(loss_rate, D) - A);
            return tranche_loss;
        };

        // Gauss-Hermite 积分: ∫ f(m) φ(m) dm = (1/√π) Σ w_k f(√2 * x_k)
        // 使用 32 点 Gauss-Hermite 节点
        return gauss_hermite_integrate_(integrand);
    }

    // 计算 t 时刻剩余名义期望 E[(D-A-L(t))^+] / (D-A)
    // = (D-A - E[L_tranche(t)]) / (D-A) = 1 - E[L_tranche(t)]/(D-A)
    Real expected_remaining_notional(Real t, Real A, Real D) const {
        Real width = D - A;
        if (width <= 0.0) return 0.0;
        Real el = expected_tranche_loss(t, A, D);
        return std::max(0.0, width - el) / width;
    }

    // ============ 完整 Tranche 定价 ============
    CDOTrancheResult price(const CDOTrancheConfig& cfg) const {
        CDOTrancheConfig cfg_v = cfg;
        cfg_v.generate_schedule();
        cfg_v.validate();

        CDOTrancheResult result;
        const Real A = cfg_v.attachment;
        const Real D = cfg_v.detachment;
        const Real width = cfg_v.tranche_width();
        const Size n = cfg_v.n_premiums;
        const Real s = cfg_v.spread;

        result.expected_tranche_loss.assign(n, 0.0);
        result.expected_remaining_notional.assign(n, 0.0);

        Real pv_premium = 0.0;
        Real pv_protection = 0.0;
        Real rpv01 = 0.0;
        Real prev_tranche_loss = 0.0;  // E[L_tranche(0)] = 0

        for (Size i = 0; i < n; ++i) {
            Real t_i = cfg_v.payment_times[i];
            Real tau = cfg_v.year_fractions[i];
            Real t_prev = (i == 0) ? 0.0 : cfg_v.payment_times[i - 1];
            Real t_mid = 0.5 * (t_prev + t_i);

            // E[L_tranche(t_i)] (fraction of total portfolio notional)
            Real el_i = expected_tranche_loss(t_i, A, D);
            result.expected_tranche_loss[i] = el_i;

            // E[remaining notional at t_i] = (width - el_i) / width (fraction of tranche notional)
            Real rn_i = std::max(0.0, width - el_i) / width;
            result.expected_remaining_notional[i] = rn_i;

            // 折现因子
            Real P_i = discount_.discount_factor(t_i);
            Real P_mid = discount_.discount_factor(t_mid);

            // Premium leg: s * width * [tau * rn_i * P_i + 0.5 * tau * (el_i - prev) * P_mid]
            // 注: width * rn_i = (width - el_i) (剩余名义 fraction)
            //     width * (el_i - prev) = ΔL_tranche(t_i) (损失变化 fraction)
            Real remaining_notional = std::max(0.0, width - el_i);
            Real delta_loss = std::max(0.0, el_i - prev_tranche_loss);
            Real cf_premium = s * tau * (remaining_notional + 0.5 * delta_loss);
            pv_premium += cf_premium * P_i;  // 用 t_i 折现 (近似)
            rpv01 += tau * (remaining_notional + 0.5 * delta_loss) * P_i;

            // Protection leg: ΔL_tranche(t_i) 折现到 t_mid
            pv_protection += delta_loss * P_mid;

            prev_tranche_loss = el_i;
        }

        result.pv_premium_leg = pv_premium;
        result.pv_protection_leg = pv_protection;
        result.risky_pv01 = rpv01;
        result.expected_loss = prev_tranche_loss;  // E[L_tranche(T)]
        result.pv = pv_premium - pv_protection;
        result.par_spread = (rpv01 > 0.0) ? pv_protection / rpv01 : 0.0;

        // 买方视角: 取反
        if (!cfg_v.is_protection_seller) {
            result.pv = -result.pv;
            result.pv_premium_leg = -result.pv_premium_leg;
            result.pv_protection_leg = -result.pv_protection_leg;
        }

        return result;
    }

    Real rho() const noexcept { return rho_; }
    Real pd_maturity() const noexcept { return pd_; }
    Real lgd() const noexcept { return lgd_; }
    Real maturity() const noexcept { return maturity_; }

private:
    Real rho_;
    Real pd_;          // PD(0, maturity)
    Real lgd_;
    const ZeroCurve& discount_;
    Real maturity_;
    Real sqrt_rho_;
    Real sqrt_1mrho_;

    // PD(t) 假设 flat hazard: PD(t) = 1 - exp(-h*t), h = -ln(1-pd_)/maturity_
    Real pd_at_time_(Real t) const {
        if (t <= 0.0) return 0.0;
        if (t >= maturity_) return pd_;
        // h = -ln(1-pd_)/maturity_
        if (pd_ >= 1.0) return 1.0;
        Real h = -std::log(1.0 - pd_) / maturity_;
        return 1.0 - std::exp(-h * t);
    }

    // 条件违约概率 PD(t|m) = Φ((Φ^{-1}(PD(t)) - √ρ * m) / √(1-ρ))
    Real conditional_pd_lhp_(Real pd_t, Real m) const {
        if (pd_t <= 0.0) return 0.0;
        if (pd_t >= 1.0) return 1.0;
        Real z = inv_normal_cdf(pd_t);
        Real x = (z - sqrt_rho_ * m) / sqrt_1mrho_;
        return normal_cdf(x);
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
            // 对称节点 ±√2 * x_k
            Real m_pos = SQRT_2 * nodes[i];
            Real m_neg = -m_pos;
            sum += weights[i] * (f(m_pos) + f(m_neg));
        }
        // 归一化: (1/√π)
        static const Real sqrt_pi = std::sqrt(PI);
        return sum / sqrt_pi;
    }
};

// ============ MC CDO 定价器 (精确, 处理异质组合) ============
// 使用 Copula 模型采样违约时间, 统计 tranche 损失分布
// 适用于任意异质组合, 但需要 MC 模拟
class CDOMCPricer {
public:
    CDOMCPricer(std::vector<CreditCurve> credit_curves,
                  const ZeroCurve& discount,
                  Real total_notional)
        : credit_curves_(std::move(credit_curves)),
          discount_(discount),
          total_notional_(total_notional) {
        if (credit_curves_.empty()) {
            throw std::invalid_argument("CDOMC: credit_curves empty");
        }
        if (total_notional_ <= 0.0) {
            throw std::invalid_argument("CDOMC: total_notional must be positive");
        }
        // 各名字等权重名义
        per_name_notional_ = total_notional_ / static_cast<Real>(credit_curves_.size());
    }

    // MC 模拟 + tranche 定价
    // copula: 已配置好的 Copula (Gaussian / OneFactor / t)
    // cfg: tranche 配置
    // n_paths: MC 路径数
    // seed: RNG 种子
    template <typename Copula>
    CDOTrancheResult price(const Copula& copula,
                              const CDOTrancheConfig& cfg,
                              Size n_paths, uint64_t seed = 42) const {
        CDOTrancheConfig cfg_v = cfg;
        cfg_v.generate_schedule();
        cfg_v.validate();

        if (copula.n_names() != credit_curves_.size()) {
            throw std::invalid_argument("CDOMC: copula.n_names != credit_curves size");
        }

        const Real A = cfg_v.attachment;
        const Real D = cfg_v.detachment;
        const Real width = cfg_v.tranche_width();
        const Size n = cfg_v.n_premiums;
        const Real s = cfg_v.spread;
        const Size n_names = credit_curves_.size();

        // 累积各时点的 E[L_tranche(t_i)] 和 E[remaining notional]
        std::vector<Real> sum_tranche_loss(n, 0.0);
        std::vector<Real> sum_remaining(n, 0.0);
        Real sum_premium = 0.0;
        Real sum_protection = 0.0;
        Real sum_rpv01 = 0.0;

        Philox4x64 rng(seed);

        for (Size path = 0; path < n_paths; ++path) {
            // 采样违约时间
            auto tau = copula.sample_default_times(rng, credit_curves_);

            // 计算各支付时点的累计损失率 (fraction of total notional)
            std::vector<Real> cum_loss_at_t(n, 0.0);
            for (Size i = 0; i < n; ++i) {
                Real t_i = cfg_v.payment_times[i];
                Real port_loss = 0.0;  // fraction
                for (Size k = 0; k < n_names; ++k) {
                    if (tau[k] <= t_i) {
                        // 名字 k 在 t_i 前违约, 损失 = LGD_k * notional_k / total_notional
                        Real loss_k = credit_curves_[k].lgd() * per_name_notional_ / total_notional_;
                        port_loss += loss_k;
                    }
                }
                cum_loss_at_t[i] = port_loss;
            }

            // 计算各时点 tranche 损失
            Real prev_tranche_loss = 0.0;
            for (Size i = 0; i < n; ++i) {
                Real port_loss = cum_loss_at_t[i];
                Real tranche_loss = std::max(0.0, std::min(port_loss, D) - A);
                Real remaining = std::max(0.0, width - tranche_loss);
                Real delta_loss = std::max(0.0, tranche_loss - prev_tranche_loss);

                Real t_i = cfg_v.payment_times[i];
                Real tau_i = cfg_v.year_fractions[i];
                Real t_prev = (i == 0) ? 0.0 : cfg_v.payment_times[i - 1];
                Real t_mid = 0.5 * (t_prev + t_i);
                Real P_i = discount_.discount_factor(t_i);
                Real P_mid = discount_.discount_factor(t_mid);

                sum_tranche_loss[i] += tranche_loss;
                sum_remaining[i] += remaining;

                // Premium leg (per unit notional * width)
                Real cf_premium = s * tau_i * (remaining + 0.5 * delta_loss);
                sum_premium += cf_premium * P_i;
                sum_rpv01 += tau_i * (remaining + 0.5 * delta_loss) * P_i;

                // Protection leg
                sum_protection += delta_loss * P_mid;

                prev_tranche_loss = tranche_loss;
            }
        }

        // 取平均
        Real inv_n = 1.0 / static_cast<Real>(n_paths);
        CDOTrancheResult result;
        result.expected_tranche_loss.resize(n);
        result.expected_remaining_notional.resize(n);
        for (Size i = 0; i < n; ++i) {
            result.expected_tranche_loss[i] = sum_tranche_loss[i] * inv_n;
            result.expected_remaining_notional[i] = sum_remaining[i] * inv_n / width;
        }
        result.pv_premium_leg = sum_premium * inv_n;
        result.pv_protection_leg = sum_protection * inv_n;
        result.risky_pv01 = sum_rpv01 * inv_n;
        result.expected_loss = result.expected_tranche_loss.back();
        result.pv = result.pv_premium_leg - result.pv_protection_leg;
        result.par_spread = (result.risky_pv01 > 0.0)
                              ? result.pv_protection_leg / result.risky_pv01 : 0.0;

        if (!cfg_v.is_protection_seller) {
            result.pv = -result.pv;
            result.pv_premium_leg = -result.pv_premium_leg;
            result.pv_protection_leg = -result.pv_protection_leg;
        }

        return result;
    }

    Size n_names() const noexcept { return credit_curves_.size(); }
    Real total_notional() const noexcept { return total_notional_; }

private:
    std::vector<CreditCurve> credit_curves_;
    const ZeroCurve& discount_;
    Real total_notional_;
    Real per_name_notional_;
};

// ============ 便捷工厂 ============

// 构造标准 CDX IG 5Y tranches (0-3%, 3-7%, 7-10%, 10-15%, 15-30%, 30-100%)
inline std::vector<CDOTrancheConfig> make_cdx_ig_tranches(Real maturity = 5.0,
                                                             Size premiums_per_year = 4) {
    std::vector<std::pair<Real, Real>> ads = {
        {0.00, 0.03}, {0.03, 0.07}, {0.07, 0.10},
        {0.10, 0.15}, {0.15, 0.30}, {0.30, 1.00}
    };
    std::vector<CDOTrancheConfig> tranches;
    tranches.reserve(ads.size());
    Size n_prem = static_cast<Size>(std::round(maturity * premiums_per_year));
    for (auto [A, D] : ads) {
        CDOTrancheConfig cfg;
        cfg.attachment = A;
        cfg.detachment = D;
        cfg.maturity = maturity;
        cfg.n_premiums = n_prem;
        tranches.push_back(cfg);
    }
    return tranches;
}

// ============ Base Correlation 标定器 ============
// SOURCE: McGinty & Ahluwalia (2004) "A model for base correlation calculation" JPMorgan
// SOURCE: O'Kane (2008) "Modelling Single-Name and Multi-Name Credit Derivatives" Ch.12-13
// 用于: (1) compound correlation 反推 (单分券隐含 rho)
//       (2) base correlation curve 标定 (CDX/iTraxx 标准分券)
//       (3) 非标准分券定价 (base corr 插值)
//
// Base correlation 方法 (McGinty-Ahluwalia):
//   每个标准分券 (A_k, D_k) 对应一个 0-D_k 的"基准分券".
//   用递归方式标定: 对 D_k, 找 rho_k 使
//     E[L_{0-D_k}(t; rho_k)] - E[L_{0-D_{k-1}}(t; rho_{k-1})]
//   构成的分券 par spread = market spread_k.
//   得到 base corr curve {(D_k, rho_k)}, 对任意 (A, D) 分券:
//     rho_A = interp(A), rho_D = interp(D)
//     E[L_{A-D}(t)] = E[L_{0-D}(t; rho_D)] - E[L_{0-A}(t; rho_A)]
class CDOBaseCorrelationCalibrator {
public:
    CDOBaseCorrelationCalibrator(Real pd_maturity, Real lgd,
                                 const ZeroCurve& discount, Real maturity)
        : pd_maturity_(pd_maturity), lgd_(lgd),
          discount_(discount), maturity_(maturity) {
        if (pd_maturity_ < 0.0 || pd_maturity_ > 1.0) {
            throw std::invalid_argument("CDOBaseCorr: pd must be in [0, 1]");
        }
        if (lgd_ < 0.0 || lgd_ > 1.0) {
            throw std::invalid_argument("CDOBaseCorr: lgd must be in [0, 1]");
        }
        if (maturity_ <= 0.0) {
            throw std::invalid_argument("CDOBaseCorr: maturity must be positive");
        }
    }

    // Compound correlation (单分券隐含 rho)
    // 返回 rho 使 CDOLHPPricer(ρ,...).price(cfg).par_spread = market_spread
    // 返回 nan 表示无解 (rho > 0.95 或 < 0)
    Real compound_correlation(Real A, Real D, Real market_spread,
                              Real tol = 1e-6, Size max_iter = 100) const {
        if (!(A >= 0.0 && D > A && D <= 1.0)) {
            throw std::invalid_argument("CDOBaseCorr: need 0 <= A < D <= 1");
        }
        auto par_spread_fn = [&](Real rho) -> Real {
            return lhp_par_spread_(rho, A, D);
        };
        return root_find_spread_(par_spread_fn, market_spread, tol, max_iter);
    }

    // 标定 base correlation curve (CDX/iTraxx 标准分券)
    // 输入: 标准分券 detachment points [0.03, 0.07, 0.10, 0.15, 0.30]
    //      对应 market par spreads (每个 (attachment, detachment) 分券)
    // 输出: 每个 detachment 的 base correlation
    std::vector<Real> calibrate_base_correlation(
        const std::vector<Real>& detachments,
        const std::vector<Real>& market_spreads,
        Real attachment_first = 0.0) const {
        if (detachments.size() != market_spreads.size()) {
            throw std::invalid_argument("CDOBaseCorr: detachments/market_spreads size mismatch");
        }
        if (detachments.empty()) return {};

        std::vector<Real> base_dets;
        std::vector<Real> base_corrs;
        base_dets.reserve(detachments.size());
        base_corrs.reserve(detachments.size());

        Real A = attachment_first;
        for (Size k = 0; k < detachments.size(); ++k) {
            Real D = detachments[k];
            if (D <= A) {
                throw std::invalid_argument("CDOBaseCorr: detachments must be strictly increasing");
            }
            base_dets.push_back(D);

            // 对 0-D_k 基准分券找 rho_k, 使 (A, D_k) 分券 par spread = market_spread_k
            // 其中 0-A 损失用已标定的 (D_{k-1}, rho_{k-1}) 插值得到 rho_A
            auto par_spread_fn = [&](Real rho) -> Real {
                base_corrs.push_back(rho);
                Real ps = base_corr_par_spread_(A, D, base_dets, base_corrs);
                base_corrs.pop_back();
                return ps;
            };
            Real rho_k = root_find_spread_(par_spread_fn, market_spreads[k], 1e-6, 100);
            base_corrs.push_back(rho_k);
            A = D;
        }
        return base_corrs;
    }

    // 用 base correlation curve 给非标准分券定价
    // 算法: 对 D 和 A 在 base corr curve 上插值得 rho_D, rho_A
    //   E[L_{A-D}(t)] = E[L_{0-D}(t; rho_D)] - E[L_{0-A}(t; rho_A)]
    // 返回: spread=0 时 par spread; 否则该 spread 下的 PV (卖方视角)
    Real price_off_market_tranche(
        Real A, Real D,
        const std::vector<Real>& base_detachments,
        const std::vector<Real>& base_correlations,
        Real spread = 0.0) const {
        if (!(A >= 0.0 && D > A && D <= 1.0)) {
            throw std::invalid_argument("CDOBaseCorr: need 0 <= A < D <= 1");
        }
        Real rho_D = interpolate_base_corr(D, base_detachments, base_correlations);
        Real rho_A = interpolate_base_corr(A, base_detachments, base_correlations);

        CDOLHPPricer pricer_D(rho_D, pd_maturity_, lgd_, discount_, maturity_);
        CDOLHPPricer pricer_A(rho_A, pd_maturity_, lgd_, discount_, maturity_);

        Real pv_protection = 0.0;
        Real rpv01 = 0.0;
        Real prev_el = 0.0;
        for (Size i = 0; i < n_premiums_(); ++i) {
            Real t_i = payment_time_(i);
            Real t_prev = (i == 0) ? 0.0 : payment_time_(i - 1);
            Real t_mid = 0.5 * (t_prev + t_i);
            Real tau = tau_();

            // E[L_{0-D}(t_i)] - E[L_{0-A}(t_i)]
            Real el_D = pricer_D.expected_tranche_loss(t_i, 0.0, D);
            Real el_A = pricer_A.expected_tranche_loss(t_i, 0.0, A);
            Real el = std::max(0.0, el_D - el_A);
            Real delta_loss = std::max(0.0, el - prev_el);
            Real remaining = std::max(0.0, (D - A) - el);

            Real P_i = discount_.discount_factor(t_i);
            Real P_mid = discount_.discount_factor(t_mid);

            rpv01 += tau * (remaining + 0.5 * delta_loss) * P_i;
            pv_protection += delta_loss * P_mid;
            prev_el = el;
        }

        if (spread == 0.0) {
            return (rpv01 > 0.0) ? pv_protection / rpv01 : 0.0;
        }
        return spread * rpv01 - pv_protection;
    }

    // 一致性诊断: base correlation 是否单调递增 (detachment ↑ → rho ↑)
    static bool is_base_correlation_monotonic(const std::vector<Real>& corrs) noexcept {
        for (Size i = 1; i < corrs.size(); ++i) {
            if (corrs[i] < corrs[i - 1]) return false;
        }
        return true;
    }

    Real pd_maturity() const noexcept { return pd_maturity_; }
    Real lgd() const noexcept { return lgd_; }
    Real maturity() const noexcept { return maturity_; }

private:
    Real pd_maturity_;
    Real lgd_;
    const ZeroCurve& discount_;
    Real maturity_;

    // 季度支付频率
    Size n_premiums_() const noexcept {
        Size n = static_cast<Size>(std::round(maturity_ * 4.0));
        return (n > 0) ? n : 1;
    }

    Real tau_() const noexcept { return maturity_ / static_cast<Real>(n_premiums_()); }

    Real payment_time_(Size i) const noexcept {
        return static_cast<Real>(i + 1) * tau_();
    }

    // LHP par spread (单一 rho)
    Real lhp_par_spread_(Real rho, Real A, Real D) const {
        CDOLHPPricer pricer(rho, pd_maturity_, lgd_, discount_, maturity_);
        CDOTrancheConfig cfg;
        cfg.attachment = A;
        cfg.detachment = D;
        cfg.spread = 0.0;
        cfg.maturity = maturity_;
        cfg.n_premiums = n_premiums_();
        return pricer.price(cfg).par_spread;
    }

    // base corr par spread: E[L_{A-D}(t)] = E[L_{0-D}(t; rho_D)] - E[L_{0-A}(t; rho_A)]
    Real base_corr_par_spread_(Real A, Real D,
                               const std::vector<Real>& dets,
                               const std::vector<Real>& corrs) const {
        return price_off_market_tranche(A, D, dets, corrs, 0.0);
    }

    // 一维求根: 找 rho ∈ [0, 0.95] 使 par_spread_fn(rho) = market_spread
    // par_spread_fn 对 rho 单调 (方向不限), 用 bisection
    // 无解 (端点同号) 返回 nan
    Real root_find_spread_(const std::function<Real(Real)>& par_spread_fn,
                           Real market_spread, Real tol, Size max_iter) const {
        const Real rho_lo = 0.0;
        const Real rho_hi = 0.95;
        Real f_lo = par_spread_fn(rho_lo) - market_spread;
        Real f_hi = par_spread_fn(rho_hi) - market_spread;

        if (std::abs(f_lo) <= tol) return rho_lo;
        if (std::abs(f_hi) <= tol) return rho_hi;
        if (!std::isfinite(f_lo) || !std::isfinite(f_hi)) {
            return std::numeric_limits<Real>::quiet_NaN();
        }
        // 端点同号 → 区间内无根 (单调函数)
        if ((f_lo > 0.0) == (f_hi > 0.0)) {
            return std::numeric_limits<Real>::quiet_NaN();
        }

        Real lo = rho_lo;
        Real hi = rho_hi;
        for (Size it = 0; it < max_iter; ++it) {
            Real mid = 0.5 * (lo + hi);
            Real f_mid = par_spread_fn(mid) - market_spread;
            if (std::abs(f_mid) <= tol || (hi - lo) < 1e-12) return mid;
            if ((f_mid > 0.0) == (f_lo > 0.0)) {
                lo = mid;
                f_lo = f_mid;
            } else {
                hi = mid;
                f_hi = f_mid;
            }
        }
        return 0.5 * (lo + hi);
    }

    // 线性插值 base corr curve (extrapolation flat)
    Real interpolate_base_corr(Real D,
                               const std::vector<Real>& detachments,
                               const std::vector<Real>& corrs) const {
        if (detachments.empty() || detachments.size() != corrs.size()) {
            return std::numeric_limits<Real>::quiet_NaN();
        }
        if (D <= detachments.front()) return corrs.front();
        if (D >= detachments.back()) return corrs.back();
        for (Size i = 0; i + 1 < detachments.size(); ++i) {
            if (D >= detachments[i] && D <= detachments[i + 1]) {
                Real width = detachments[i + 1] - detachments[i];
                if (width <= 0.0) return corrs[i];
                Real w = (D - detachments[i]) / width;
                return corrs[i] + w * (corrs[i + 1] - corrs[i]);
            }
        }
        return corrs.back();
    }
};

}  // namespace v1
}  // namespace cpphub
