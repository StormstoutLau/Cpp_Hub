#pragma once
// SOURCE: Vasicek (1977) "An equilibrium characterization of the term structure"
// SOURCE: Cox-Ingersoll-Ross (1985) "A Theory of the Term Structure of Interest Rates"
// SOURCE: Hull & White (1990) "Pricing interest-rate derivative securities"
// SOURCE: Brigo & Mercurio (2006) "Interest Rate Models - Theory and Practice" Ch.3-4
// 模块: 仿射短期利率模型族 (Vasicek / CIR / Hull-White / G2++)
//
// 仿射期限结构 (ATSM) 统一形式:
//   dr(t) = (θ(t) - κ*r(t)) dt + σ(t) * r(t)^β dW(t)
//   β=0 → Vasicek/Hull-White (高斯, ORNSTEIN-UHLENBECK)
//   β=1/2 → CIR (平方根, Feller 条件 2κθ > σ²)
//
// 零息债价格 P(t,T) = A(t,T) * exp(-B(t,T) * r(t))
//   B(κ,T-t) 仿射求解, A(κ,T-t) 依赖波动率结构
//
// G2++ (Brigo-Mercurio 2-factor Gaussian):
//   r(t) = x(t) + y(t) + φ(t)
//   dx = -a*x dt + σ dW1
//   dy = -b*y dt + η dW2,  d<W1,W2> = ρ dt
//   P(0,T) = exp(-φ_T) * A_2F(a,b,ρ,σ,η,T) * exp(-B_a(T)*x0 - B_b(T)*y0)
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/models/diffusion/process.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ 单因子短期利率模型参数 ============

struct VasicekParams {
    Real r0;       // 当前短期利率
    Real kappa;    // 均值回复速度 (>0)
    Real theta;    // 长期均值
    Real sigma;    // 波动率 (>0)
};

struct CIRParams {
    Real r0;
    Real kappa;
    Real theta;
    Real sigma;
    // Feller 条件: 2*kappa*theta > sigma^2 (保证 r > 0)
};

struct HullWhiteParams {
    Real r0;
    Real kappa;    // a(t) 通常常数, 但可时变
    Real sigma;
    // θ(t) 由初始期限结构匹配决定 (calibration to P(0,T))
};

struct G2Params {
    Real x0 = 0.0;
    Real y0 = 0.0;
    Real a = 0.1;     // x 因子回复速度
    Real b = 0.1;     // y 因子回复速度
    Real sigma = 0.01;  // x 因子波动率
    Real eta = 0.01;    // y 因子波动率
    Real rho = 0.0;     // 两因子相关性
    // φ(t) 由初始期限结构匹配
};

// ============ B 函数 (零息债价格指数项) ============

// Vasicek/CIR/HW: B(κ, τ) = (1 - exp(-κ*τ)) / κ
// 当 κ → 0 时退化为 τ
inline Real affine_B(Real kappa, Real tau) {
    if (kappa < 1e-12) return tau;
    return (1.0 - std::exp(-kappa * tau)) / kappa;
}

// ============ Vasicek 模型 (Ornstein-Uhlenbeck) ============
// dr = κ(θ - r) dt + σ dW
// r(t) ~ N(r0*exp(-κt) + θ*(1-exp(-κt)), σ²*(1-exp(-2κt))/(2κ))
// P(0,T) = A(τ)*exp(-B(τ)*r0)
// A(τ) = exp( (B(τ) - τ)*(θ - σ²/(2κ²))*κ - σ²*B(τ)²/(4κ) )

class Vasicek {
public:
    explicit Vasicek(VasicekParams p) : params_(p) {
        validate();
    }

    // 零息债价格 P(0,T)
    Real zero_coupon_bond(Real T) const {
        if (T <= 0.0) return 1.0;
        Real tau = T;
        Real B = affine_B(params_.kappa, tau);
        Real A = compute_A(tau, B);
        return A * std::exp(-B * params_.r0);
    }

    // 零息债收益率 (连续复利)
    Real yield(Real T) const {
        if (T <= 0.0) return params_.r0;
        return -std::log(zero_coupon_bond(T)) / T;
    }

    // 零息债价格关于 r 的偏导 (用于构建利率二叉树)
    Real bond_delta(Real T) const {
        return -affine_B(params_.kappa, T) * zero_coupon_bond(T);
    }

    // 债券期权 (Jamshidian 1989): 对 Vasicek 有解析解
    // Call on P(T,T+S): max(P(T,T+S) - K, 0)
    // P(t,T) = A(t,T) exp(-B(t,T) r(t))
    // P(T,T+S) 在 T 时刻是 r(T) 的指数函数 → Jamshidian 分解
    Real bond_option(Real T_opt, Real T_bond, Real K, bool is_call) const {
        if (T_opt <= 0.0 || T_bond <= T_opt) {
            throw std::invalid_argument("Vasicek::bond_option: require T_opt>0 and T_bond>T_opt");
        }
        Real tau = T_bond - T_opt;  // 债券在 T_opt 时点的剩余期限
        Real B_tau = affine_B(params_.kappa, tau);
        Real A_tau = compute_A(tau, B_tau);

        // T_opt 时刻 r(T_opt) 的分布
        Real B_T = affine_B(params_.kappa, T_opt);
        Real A_T = compute_A(T_opt, B_T);
        // r(T) 的均值和方差
        Real m_r = params_.theta + (params_.r0 - params_.theta) * std::exp(-params_.kappa * T_opt);
        Real var_r = params_.sigma * params_.sigma * (1.0 - std::exp(-2.0 * params_.kappa * T_opt))
                     / (2.0 * params_.kappa);
        Real sigma_r = std::sqrt(var_r);

        // Jamshidian: 找 r* 使得 A_tau * exp(-B_tau * r*) = K
        // → r* = -log(K / A_tau) / B_tau
        Real r_star = -std::log(K / A_tau) / B_tau;

        // P(T_opt, T_bond) 的均值 (在风险中性测度下)
        // E[P(T,T+S)] = A_tau * E[exp(-B_tau * r(T))]
        //             = A_tau * exp(-B_tau * m_r + 0.5 * B_tau^2 * var_r)
        Real mean_P = A_tau * std::exp(-B_tau * m_r + 0.5 * B_tau * B_tau * var_r);

        // 零息债期权 (Jamshidian 分解为多个零息债期权之和, 这里单一零息债)
        // Call = P(0,T_opt) * E_Q[max(P(T_opt,T_bond) - K, 0) | F_0]
        // 由于 P(T,T+S) = A_tau * exp(-B_tau * r(T)), 且 r(T) 服从高斯分布,
        // 可视为对高斯变量 r(T) 的指数函数期权, 用类似 Black-Scholes 公式
        Real P_0_T = zero_coupon_bond(T_opt);  // 折现因子
        // 等价: 把 P(T_opt, T_bond) 视为标的资产, 其对数正态分布的均值和方差
        // ln P(T,T+S) = ln A_tau - B_tau * r(T)
        // E[ln P] = ln A_tau - B_tau * m_r
        // Var[ln P] = B_tau^2 * var_r
        Real mu_lnP = std::log(A_tau) - B_tau * m_r;
        Real var_lnP = B_tau * B_tau * var_r;
        Real sigma_P = std::sqrt(var_lnP);
        // 等价 Black 公式 (P 为对数正态):
        // Forward F = E[P(T,T+S)] = mean_P
        // d1 = (ln(F/K) + 0.5*sigma_P^2) / sigma_P
        // d2 = d1 - sigma_P
        Real F = mean_P;
        if (F <= 0.0 || K <= 0.0 || sigma_P < 1e-15) {
            // 退化情形: 直接用内在价值
            Real payoff = is_call ? std::max(F - K, 0.0) : std::max(K - F, 0.0);
            return P_0_T * payoff;
        }
        Real d1 = (std::log(F / K) + 0.5 * var_lnP) / sigma_P;
        Real d2 = d1 - sigma_P;
        if (is_call) {
            return P_0_T * (F * normal_cdf(d1) - K * normal_cdf(d2));
        } else {
            return P_0_T * (K * normal_cdf(-d2) - F * normal_cdf(-d1));
        }
    }

    // 模拟路径 (Euler 方案)
    void simulate_path(Real T, Size n_steps, std::vector<Real>& path, Philox4x64& rng) const {
        path.resize(n_steps + 1);
        path[0] = params_.r0;
        Real dt = T / static_cast<Real>(n_steps);
        Real sqrt_dt = std::sqrt(dt);
        for (Size i = 1; i <= n_steps; ++i) {
            // Box-Muller
            Real u1 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
            if (u1 < 1e-300) u1 = 1e-300;
            Real u2 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
            Real z = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586476925286766559 * u2);
            Real r_prev = path[i - 1];
            path[i] = r_prev + params_.kappa * (params_.theta - r_prev) * dt
                     + params_.sigma * sqrt_dt * z;
        }
    }

    // 短期利率的条件期望与方差 (解析解, 用于验证模拟)
    void conditional_moments(Real T, Real& mean, Real& variance) const {
        Real e_kt = std::exp(-params_.kappa * T);
        mean = params_.theta + (params_.r0 - params_.theta) * e_kt;
        variance = params_.sigma * params_.sigma * (1.0 - e_kt * e_kt) / (2.0 * params_.kappa);
    }

    const VasicekParams& params() const { return params_; }

private:
    VasicekParams params_;

    void validate() const {
        if (params_.kappa <= 0.0) throw std::invalid_argument("Vasicek: kappa must be positive");
        if (params_.sigma <= 0.0) throw std::invalid_argument("Vasicek: sigma must be positive");
    }

    Real compute_A(Real tau, Real B) const {
        Real kappa = params_.kappa;
        Real theta = params_.theta;
        Real sigma = params_.sigma;
        // A(τ) = exp( (B - τ)*(κθ - σ²/(2κ)) - σ²*B²/(4κ) )
        // 推导: A = exp[ (B-τ)*(κθ - σ²/2κ²)*κ - σ²*B²/(4κ) ]
        //       = exp[ (B-τ)*(κθ - σ²/(2κ)) - σ²*B²/(4κ) ]
        Real term1 = (B - tau) * (kappa * theta - sigma * sigma / (2.0 * kappa));
        Real term2 = sigma * sigma * B * B / (4.0 * kappa);
        return std::exp(term1 - term2);
    }
};

// ============ CIR 模型 (平方根过程) ============
// dr = κ(θ - r) dt + σ*sqrt(r) dW
// Feller 条件: 2κθ > σ² (保证 r > 0 a.s.)
// P(0,T) = A(τ)*exp(-B(τ)*r0)
// B(τ) = 2*(exp(γτ) - 1) / ((κ+γ)*(exp(γτ) - 1) + 2γ)
// A(τ) = [ 2γ*exp((κ+γ)τ/2) / ((κ+γ)*(exp(γτ)-1) + 2γ) ]^(2κθ/σ²)
// γ = sqrt(κ² + 2σ²)

class CIR {
public:
    explicit CIR(CIRParams p) : params_(p) {
        validate();
    }

    Real zero_coupon_bond(Real T) const {
        if (T <= 0.0) return 1.0;
        Real gamma = std::sqrt(params_.kappa * params_.kappa
                               + 2.0 * params_.sigma * params_.sigma);
        Real tau = T;
        Real exp_gt = std::exp(gamma * tau);
        Real denom = (params_.kappa + gamma) * (exp_gt - 1.0) + 2.0 * gamma;
        Real B = 2.0 * (exp_gt - 1.0) / denom;
        Real A_base = 2.0 * gamma * std::exp((params_.kappa + gamma) * tau / 2.0) / denom;
        Real exponent = 2.0 * params_.kappa * params_.theta
                       / (params_.sigma * params_.sigma);
        Real A = std::pow(A_base, exponent);
        return A * std::exp(-B * params_.r0);
    }

    Real yield(Real T) const {
        if (T <= 0.0) return params_.r0;
        return -std::log(zero_coupon_bond(T)) / T;
    }

    bool feller_satisfied() const {
        return 2.0 * params_.kappa * params_.theta > params_.sigma * params_.sigma;
    }

    // 模拟路径 (Euler 全截断方案, Glasserman Ch.3)
    void simulate_path(Real T, Size n_steps, std::vector<Real>& path, Philox4x64& rng) const {
        path.resize(n_steps + 1);
        path[0] = params_.r0;
        Real dt = T / static_cast<Real>(n_steps);
        Real sqrt_dt = std::sqrt(dt);
        for (Size i = 1; i <= n_steps; ++i) {
            Real u1 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
            if (u1 < 1e-300) u1 = 1e-300;
            Real u2 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
            Real z = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586476925286766559 * u2);
            Real r_prev = path[i - 1];
            // Full truncation: 用 max(r, 0) 计算 sqrt
            Real r_pos = std::max(r_prev, 0.0);
            path[i] = r_prev + params_.kappa * (params_.theta - r_prev) * dt
                     + params_.sigma * std::sqrt(r_pos) * sqrt_dt * z;
            // 强制非负
            if (path[i] < 0.0) path[i] = 0.0;
        }
    }

    // 条件期望 (无解析分布的高阶矩, 但一阶矩有闭式)
    void conditional_moments(Real T, Real& mean, Real& variance) const {
        Real e_kt = std::exp(-params_.kappa * T);
        mean = params_.theta + (params_.r0 - params_.theta) * e_kt;
        // 方差: σ² r0 (e^{-κT} - e^{-2κT})/κ + θσ² (1-e^{-κT})²/(2κ)
        Real sigma2 = params_.sigma * params_.sigma;
        Real kappa = params_.kappa;
        variance = sigma2 * params_.r0 * (e_kt - e_kt * e_kt) / kappa
                 + params_.theta * sigma2 * (1.0 - e_kt) * (1.0 - e_kt) / (2.0 * kappa);
    }

    const CIRParams& params() const { return params_; }

private:
    CIRParams params_;

    void validate() const {
        if (params_.kappa <= 0.0) throw std::invalid_argument("CIR: kappa must be positive");
        if (params_.theta <= 0.0) throw std::invalid_argument("CIR: theta must be positive");
        if (params_.sigma <= 0.0) throw std::invalid_argument("CIR: sigma must be positive");
        if (params_.r0 < 0.0) throw std::invalid_argument("CIR: r0 must be non-negative");
    }
};

// ============ Hull-White 模型 (时变 θ(t), 无套利校准) ============
// dr = (θ(t) - κ*r) dt + σ dW
// θ(t) = ∂f(0,t)/∂t + κ*f(0,t) + σ²/(2κ)*(1 - exp(-2κt))
// 其中 f(0,t) = -∂ln P(0,t)/∂t 为初始远期利率曲线
//
// 给定初始零息债曲线 P(0,T), HW 可精确匹配当前期限结构

class HullWhite {
public:
    HullWhite(HullWhiteParams p, std::vector<Real> maturities, std::vector<Real> bond_prices)
        : params_(p), maturities_(std::move(maturities)), bonds_(std::move(bond_prices)) {
        validate();
    }

    // 零息债价格 P(0,T) — 从初始曲线插值
    Real zero_coupon_bond(Real T) const {
        if (T <= 0.0) return 1.0;
        return interpolate_bond(T);
    }

    // 远期利率 f(0,T) = -∂ln P(0,T)/∂T
    Real forward_rate(Real T) const {
        if (T <= 0.0) return params_.r0;
        // 数值微分 (中心差分)
        Real h = 1e-6;
        Real P_plus = interpolate_bond(T + h);
        Real P_minus = interpolate_bond(std::max(T - h, 1e-10));
        return -(std::log(P_plus) - std::log(P_minus)) / (2.0 * h);
    }

    // θ(t) — Hull-White 时变漂移
    Real theta(Real t) const {
        if (t <= 0.0) return params_.kappa * params_.r0;
        Real f_t = forward_rate(t);
        Real f_prime;  // ∂f/∂t
        Real h = 1e-5;
        Real f_plus = forward_rate(t + h);
        Real f_minus = forward_rate(std::max(t - h, 1e-10));
        f_prime = (f_plus - f_minus) / (2.0 * h);
        Real kappa = params_.kappa;
        Real sigma = params_.sigma;
        return f_prime + kappa * f_t + sigma * sigma / (2.0 * kappa) * (1.0 - std::exp(-2.0 * kappa * t));
    }

    // Hull-White 零息债 P(0,T) (应等于输入曲线, 用于验证校准)
    // P_HW(0,T) = A(0,T) * exp(-B(0,T) * r0)
    // B(0,T) = (1 - exp(-κT))/κ
    // A(0,T) = P_M(0,T) * exp(B(0,T) * f(0,0))  [Brigo-Mercurio eq 3.33, 在 t=0 方差项消失]
    // → P(0,T) = P_M(0,T) * exp(B(0,T) * (f(0,0) - r0))
    // 当 r0 = f(0,0) (校准后) 时, P(0,T) = P_M(0,T) 恒等
    Real hw_zero_coupon_bond(Real T) const {
        if (T <= 0.0) return 1.0;
        Real B = affine_B(params_.kappa, T);
        Real P_market = interpolate_bond(T);
        // t=0 时方差项 (1 - exp(-2κ*0)) = 0, 仅剩余短期利率匹配项
        Real P_hw = P_market * std::exp(-B * (params_.r0 - forward_rate(0.0)));
        return P_hw;
    }

    // 债券期权 (Black 形式, 类似 Vasicek 但 θ 时变)
    // Call on P(T_opt, T_bond): 用 HW 闭式解
    Real bond_option(Real T_opt, Real T_bond, Real K, bool is_call) const {
        if (T_opt <= 0.0 || T_bond <= T_opt) {
            throw std::invalid_argument("HullWhite::bond_option: require T_opt>0 and T_bond>T_opt");
        }
        Real tau = T_bond - T_opt;
        Real B_tau = affine_B(params_.kappa, tau);
        Real B_T = affine_B(params_.kappa, T_opt);

        // P(T_opt, T_bond) = A(T_opt, T_bond) * exp(-B_tau * r(T_opt))
        // r(T_opt) 服从高斯分布 (HW 仍是高斯过程)
        // 远期债券价格 F = P(0,T_bond) / P(0,T_opt)
        Real P_0_Tbond = interpolate_bond(T_bond);
        Real P_0_Topt = interpolate_bond(T_opt);
        Real F = P_0_Tbond / P_0_Topt;

        // P(T_opt, T_bond) 的对数方差
        // σ_P² = σ² * (1-exp(-2κT_opt))/(2κ) * B_tau²
        Real var_r = params_.sigma * params_.sigma
                    * (1.0 - std::exp(-2.0 * params_.kappa * T_opt))
                    / (2.0 * params_.kappa);
        Real var_lnP = B_tau * B_tau * var_r;
        Real sigma_P = std::sqrt(var_lnP);

        if (sigma_P < 1e-15 || F <= 0.0 || K <= 0.0) {
            Real payoff = is_call ? std::max(F - K, 0.0) : std::max(K - F, 0.0);
            return P_0_Topt * payoff;
        }

        Real d1 = (std::log(F / K) + 0.5 * var_lnP) / sigma_P;
        Real d2 = d1 - sigma_P;
        if (is_call) {
            return P_0_Topt * (F * normal_cdf(d1) - K * normal_cdf(d2));
        } else {
            return P_0_Topt * (K * normal_cdf(-d2) - F * normal_cdf(-d1));
        }
    }

    // Jamshidian (1989) 分解: 求 r* 使 Σ_i c_i * A(T_opt, t_i) * exp(-B(T_opt, t_i) * r*) = K
    // P(T_opt, t_i) = A(T_opt, t_i) * exp(-B(T_opt, t_i) * r(T_opt))
    // 用带括号保护的 Newton 迭代求解 (通常在 20 步内二次收敛)
    Real jamshidian_r_star(Real T_opt,
                           const std::vector<Real>& payment_times,
                           const std::vector<Real>& cashflows,
                           Real K) const {
        if (payment_times.size() != cashflows.size()) {
            throw std::invalid_argument("HullWhite::jamshidian_r_star: size mismatch");
        }
        if (payment_times.empty()) {
            throw std::invalid_argument("HullWhite::jamshidian_r_star: empty schedule");
        }
        if (T_opt <= 0.0) {
            throw std::invalid_argument("HullWhite::jamshidian_r_star: require T_opt > 0");
        }
        if (K <= 0.0) {
            throw std::invalid_argument("HullWhite::jamshidian_r_star: require K > 0");
        }
        std::vector<Real> times, cfs;
        filter_payments(T_opt, payment_times, cashflows, times, cfs);
        if (times.empty()) {
            throw std::invalid_argument("HullWhite::jamshidian_r_star: no payments after T_opt");
        }
        std::vector<Real> A_terms, B_terms;
        build_affine_terms(T_opt, times, A_terms, B_terms);
        return solve_r_star(K, A_terms, B_terms, cfs);
    }

    // Coupon-bearing bond option via Jamshidian 分解
    // cashflows[i] 在 payment_times[i] 时支付; call = Σ c_i * ZBP(T_opt, t_i, K_i)
    Real coupon_bond_option(Real T_opt,
                            const std::vector<Real>& payment_times,
                            const std::vector<Real>& cashflows,
                            Real K, bool is_call) const {
        if (payment_times.size() != cashflows.size()) {
            throw std::invalid_argument("HullWhite::coupon_bond_option: size mismatch");
        }
        if (payment_times.empty()) {
            throw std::invalid_argument("HullWhite::coupon_bond_option: empty schedule");
        }
        if (T_opt <= 0.0) {
            throw std::invalid_argument("HullWhite::coupon_bond_option: require T_opt > 0");
        }
        if (K <= 0.0) {
            throw std::invalid_argument("HullWhite::coupon_bond_option: require K > 0");
        }
        std::vector<Real> times, cfs;
        filter_payments(T_opt, payment_times, cashflows, times, cfs);
        if (times.empty()) {
            // T_opt 之后无剩余现金流: 期权无内在价值
            return is_call ? 0.0 : K * interpolate_bond(T_opt);
        }
        Real r_star = jamshidian_r_star(T_opt, times, cfs, K);
        std::vector<Real> A_terms, B_terms;
        build_affine_terms(T_opt, times, A_terms, B_terms);
        Real sum = 0.0;
        for (Size i = 0; i < times.size(); ++i) {
            Real K_i = A_terms[i] * std::exp(-B_terms[i] * r_star);
            sum += cfs[i] * bond_option(T_opt, times[i], K_i, is_call);
        }
        return sum;
    }

    // Swaption (payer/receiver) via Jamshidian 分解
    // 把 payer swap 视为组合: +1 at T_start, -K*τ 在中间支付日, -(1+K*τ) at T_end
    // V_payer(r) 关于 r 单调递增, V_payer(r*)=0 ⟹
    //   payer   = -Σ_i c_i * put(ZCB_i, K_i),   receiver = -Σ_i c_i * call(ZCB_i, K_i)
    Real swaption(Real T_opt, Real T_start, Real T_end,
                  Size n_periods, Real K, bool is_payer) const {
        if (T_opt <= 0.0 || T_opt > T_start) {
            throw std::invalid_argument("HullWhite::swaption: require 0 < T_opt <= T_start");
        }
        if (T_start >= T_end) {
            throw std::invalid_argument("HullWhite::swaption: require T_start < T_end");
        }
        if (n_periods == 0) {
            throw std::invalid_argument("HullWhite::swaption: n_periods must be positive");
        }
        if (K < 0.0) {
            throw std::invalid_argument("HullWhite::swaption: require K >= 0");
        }
        Real tau = (T_end - T_start) / static_cast<Real>(n_periods);
        std::vector<Real> times;
        std::vector<Real> coeffs;
        times.reserve(n_periods + 1);
        coeffs.reserve(n_periods + 1);
        times.push_back(T_start);
        coeffs.push_back(1.0);
        for (Size i = 1; i <= n_periods; ++i) {
            times.push_back(T_start + static_cast<Real>(i) * tau);
            coeffs.push_back((i == n_periods) ? -(1.0 + K * tau) : -K * tau);
        }
        std::vector<Real> A_terms, B_terms;
        build_affine_terms(T_opt, times, A_terms, B_terms);
        Real r_star = solve_r_star(0.0, A_terms, B_terms, coeffs);
        Real sum = 0.0;
        for (Size i = 0; i < times.size(); ++i) {
            if (times[i] <= T_opt + 1e-12) continue;  // 到期时确定, 期权价值为 0
            Real K_i = A_terms[i] * std::exp(-B_terms[i] * r_star);
            sum += -coeffs[i] * bond_option(T_opt, times[i], K_i, !is_payer);
        }
        return sum;
    }

    // HW 下 caplet = put on ZCB: caplet = (1 + τ*K) * put(T_start, T_end, 1/(1 + τ*K))
    // 推导: τ*max(L-K,0) = (1/P)*max(1 - (1+τK)*P, 0), 在 Q^{t_start} 测度下取期望
    Real caplet(Real t_start, Real t_end, Real K_cap) const {
        if (t_end <= t_start) return 0.0;
        if (K_cap < 0.0) {
            throw std::invalid_argument("HullWhite::caplet: require K_cap >= 0");
        }
        Real tau = t_end - t_start;
        if (t_start <= 0.0) {
            // 起始 LIBOR 已知, 直接取内在价值
            Real P_end = interpolate_bond(t_end);
            Real L0 = (1.0 / tau) * (1.0 / P_end - 1.0);
            return P_end * tau * std::max(L0 - K_cap, 0.0);
        }
        Real factor = 1.0 + tau * K_cap;
        Real put = bond_option(t_start, t_end, 1.0 / factor, false);
        return factor * put;
    }

    // Floorlet = call on ZCB (对称)
    Real floorlet(Real t_start, Real t_end, Real K_floor) const {
        if (t_end <= t_start) return 0.0;
        if (K_floor < 0.0) {
            throw std::invalid_argument("HullWhite::floorlet: require K_floor >= 0");
        }
        Real tau = t_end - t_start;
        if (t_start <= 0.0) {
            Real P_end = interpolate_bond(t_end);
            Real L0 = (1.0 / tau) * (1.0 / P_end - 1.0);
            return P_end * tau * std::max(K_floor - L0, 0.0);
        }
        Real factor = 1.0 + tau * K_floor;
        Real call = bond_option(t_start, t_end, 1.0 / factor, true);
        return factor * call;
    }

    // Cap = Σ caplets, 每个 caplet 覆盖 [reset_times[i], reset_times[i+1])
    Real cap(const std::vector<Real>& reset_times, Real K_cap) const {
        if (reset_times.size() < 2) return 0.0;
        Real sum = 0.0;
        for (Size i = 0; i + 1 < reset_times.size(); ++i) {
            sum += caplet(reset_times[i], reset_times[i + 1], K_cap);
        }
        return sum;
    }

    // Floor = Σ floorlets
    Real floor(const std::vector<Real>& reset_times, Real K_floor) const {
        if (reset_times.size() < 2) return 0.0;
        Real sum = 0.0;
        for (Size i = 0; i + 1 < reset_times.size(); ++i) {
            sum += floorlet(reset_times[i], reset_times[i + 1], K_floor);
        }
        return sum;
    }

    // 模拟路径 (Euler, 使用 θ(t))
    void simulate_path(Real T, Size n_steps, std::vector<Real>& path, Philox4x64& rng) const {
        path.resize(n_steps + 1);
        path[0] = params_.r0;
        Real dt = T / static_cast<Real>(n_steps);
        Real sqrt_dt = std::sqrt(dt);
        for (Size i = 1; i <= n_steps; ++i) {
            Real t = static_cast<Real>(i - 1) * dt;
            Real u1 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
            if (u1 < 1e-300) u1 = 1e-300;
            Real u2 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
            Real z = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586476925286766559 * u2);
            Real r_prev = path[i - 1];
            Real theta_t = theta(t);
            path[i] = r_prev + (theta_t - params_.kappa * r_prev) * dt
                     + params_.sigma * sqrt_dt * z;
        }
    }

    const HullWhiteParams& params() const { return params_; }
    const std::vector<Real>& maturities() const { return maturities_; }
    const std::vector<Real>& bond_prices() const { return bonds_; }

private:
    HullWhiteParams params_;
    std::vector<Real> maturities_;  // 初始期限结构: 到期时间
    std::vector<Real> bonds_;        // 对应的零息债价格 P(0,T)

    void validate() const {
        if (params_.kappa <= 0.0) throw std::invalid_argument("HullWhite: kappa must be positive");
        if (params_.sigma <= 0.0) throw std::invalid_argument("HullWhite: sigma must be positive");
        if (maturities_.size() != bonds_.size()) {
            throw std::invalid_argument("HullWhite: maturities and bonds size mismatch");
        }
        // interpolate_bond 需要至少 2 个点用于线性插值/外推
        if (maturities_.size() < 2) {
            throw std::invalid_argument("HullWhite: term structure requires at least 2 points");
        }
        for (Size i = 0; i < maturities_.size(); ++i) {
            if (bonds_[i] <= 0.0 || bonds_[i] > 1.0 + 1e-10) {
                throw std::invalid_argument("HullWhite: bond price must be in (0, 1]");
            }
        }
        // 期限必须严格递增 (插值要求)
        for (Size i = 1; i < maturities_.size(); ++i) {
            if (maturities_[i] <= maturities_[i - 1]) {
                throw std::invalid_argument("HullWhite: maturities must be strictly increasing");
            }
        }
    }

    // 线性插值零息债价格 (简单方案; 实际可使用对数插值或三次样条)
    Real interpolate_bond(Real T) const {
        if (T <= maturities_.front()) {
            // 短端线性外推 (基于瞬时远期利率)
            Real T0 = maturities_[0];
            Real T1 = maturities_[1];
            Real P0 = bonds_[0];
            Real P1 = bonds_[1];
            Real y0 = -std::log(P0) / T0;
            Real y1 = -std::log(P1) / T1;
            Real y = y0 + (y1 - y0) * (T - T0) / (T1 - T0);
            return std::exp(-y * T);
        }
        if (T >= maturities_.back()) {
            // 长端对数线性外推
            Real Tn = maturities_.back();
            Real Pn = bonds_.back();
            Real Tn_1 = maturities_[maturities_.size() - 2];
            Real Pn_1 = bonds_[bonds_.size() - 2];
            Real y_n = -std::log(Pn) / Tn;
            Real y_n_1 = -std::log(Pn_1) / Tn_1;
            Real y = y_n + (y_n - y_n_1) * (T - Tn) / (Tn - Tn_1);
            return std::exp(-y * T);
        }
        // 中间: 对数线性插值 (保证 P 单调递减)
        auto it = std::upper_bound(maturities_.begin(), maturities_.end(), T);
        Size idx = static_cast<Size>(it - maturities_.begin());
        Real T1 = maturities_[idx - 1];
        Real T2 = maturities_[idx];
        Real P1 = bonds_[idx - 1];
        Real P2 = bonds_[idx];
        Real w = (T - T1) / (T2 - T1);
        Real lnP = std::log(P1) + w * (std::log(P2) - std::log(P1));
        return std::exp(lnP);
    }

    // P(t,T) = A(t,T) * exp(-B(t,T) * r(t)), 校准到初始期限结构
    // A(t,T) = P(0,T)/P(0,t) * exp( B(t,T) f(0,t) - σ²/(4κ)(1-e^{-2κt}) B(t,T)² )
    // (Brigo-Mercurio 2006 eq 3.33; 与 bond_option 的 Black 形式一致)
    Real hw_A(Real t, Real T) const {
        if (T <= t) return 1.0;
        Real B = affine_B(params_.kappa, T - t);
        Real P0T = interpolate_bond(T);
        Real P0t = interpolate_bond(t);
        Real f0t = forward_rate(t);
        Real kappa = params_.kappa;
        Real sigma = params_.sigma;
        Real correction = sigma * sigma / (4.0 * kappa)
                        * (1.0 - std::exp(-2.0 * kappa * t)) * B * B;
        return (P0T / P0t) * std::exp(B * f0t - correction);
    }

    // 仅保留 T_opt 之后的支付
    void filter_payments(Real T_opt,
                         const std::vector<Real>& payment_times,
                         const std::vector<Real>& cashflows,
                         std::vector<Real>& times,
                         std::vector<Real>& cfs) const {
        times.clear();
        cfs.clear();
        for (Size i = 0; i < payment_times.size(); ++i) {
            if (payment_times[i] > T_opt + 1e-12) {
                times.push_back(payment_times[i]);
                cfs.push_back(cashflows[i]);
            }
        }
    }

    void build_affine_terms(Real T_opt, const std::vector<Real>& times,
                            std::vector<Real>& A_terms,
                            std::vector<Real>& B_terms) const {
        A_terms.resize(times.size());
        B_terms.resize(times.size());
        for (Size i = 0; i < times.size(); ++i) {
            A_terms[i] = hw_A(T_opt, times[i]);
            B_terms[i] = affine_B(params_.kappa, times[i] - T_opt);
        }
    }

    // F(r) = Σ c_i * A_i * exp(-B_i * r) - target 及其导数
    void value_and_deriv(Real r,
                         const std::vector<Real>& A_terms,
                         const std::vector<Real>& B_terms,
                         const std::vector<Real>& cashflows,
                         Real target,
                         Real& val, Real& deriv) const {
        val = -target;
        deriv = 0.0;
        for (Size i = 0; i < cashflows.size(); ++i) {
            Real e = A_terms[i] * std::exp(-B_terms[i] * r);
            val += cashflows[i] * e;
            deriv -= cashflows[i] * B_terms[i] * e;
        }
    }

    // 带括号保护的 Newton 迭代 (coupon bond: F 单调递减; swaption: F 单调递增)
    Real solve_r_star(Real target,
                      const std::vector<Real>& A_terms,
                      const std::vector<Real>& B_terms,
                      const std::vector<Real>& cashflows) const {
        if (cashflows.empty()) {
            throw std::invalid_argument("HullWhite::solve_r_star: empty schedule");
        }
        auto F = [&](Real r) {
            Real v, d;
            value_and_deriv(r, A_terms, B_terms, cashflows, target, v, d);
            return v;
        };
        auto dF = [&](Real r) {
            Real v, d;
            value_and_deriv(r, A_terms, B_terms, cashflows, target, v, d);
            return d;
        };
        // 初始括号, 双向倍增外扩直至符号相反
        Real r_lo = -0.25;
        Real r_hi = 0.25;
        Real f_lo = F(r_lo);
        Real f_hi = F(r_hi);
        for (Size i = 0; i < 200 && f_lo * f_hi > 0.0; ++i) {
            Real width = r_hi - r_lo;
            r_lo -= width;
            r_hi += width;
            f_lo = F(r_lo);
            f_hi = F(r_hi);
        }
        if (f_lo * f_hi > 0.0) {
            throw std::runtime_error("HullWhite::solve_r_star: failed to bracket root");
        }
        // Newton 迭代 (20 步内二次收敛), 越界时退化为二分
        Real r = 0.5 * (r_lo + r_hi);
        for (Size iter = 0; iter < 60; ++iter) {
            Real f = F(r);
            if (std::abs(f) < 1e-13) break;
            Real df = dF(r);
            Real r_new = (std::abs(df) < 1e-18) ? 0.5 * (r_lo + r_hi) : r - f / df;
            if (!(r_new > r_lo && r_new < r_hi)) r_new = 0.5 * (r_lo + r_hi);
            Real f_new = F(r_new);
            if (f_lo * f_new <= 0.0) {
                r_hi = r_new;
                f_hi = f_new;
            } else {
                r_lo = r_new;
                f_lo = f_new;
            }
            r = r_new;
            if (std::abs(r_hi - r_lo) < 1e-14 * std::max(1.0, std::abs(r))) break;
        }
        // 二分抛光保证精度
        for (Size i = 0; i < 60; ++i) {
            if (r_hi - r_lo < 1e-15 * std::max(1.0, std::abs(r_lo))) break;
            Real mid = 0.5 * (r_lo + r_hi);
            Real f_mid = F(mid);
            if (f_lo * f_mid <= 0.0) r_hi = mid;
            else r_lo = mid;
        }
        return 0.5 * (r_lo + r_hi);
    }
};

// ============ G2++ (两因子高斯模型, Brigo-Mercurio) ============
// r(t) = x(t) + y(t) + φ(t)
// dx = -a*x dt + σ dW1,  x(0) = x0
// dy = -b*y dt + η dW2,  y(0) = y0
// d<W1,W2> = ρ dt
// φ(t) 由初始期限结构校准 (类似 HW 的 θ(t))
//
// 零息债 (忽略 φ 校准, 仅计算仿射部分):
// P_aff(0,T) = exp(-B_a(T)*x0 - B_b(T)*y0 + V(0,T)/2)   [Brigo-Mercurio eq 4.15-4.16]
// V(0,T) = σ²/a²*(T + (1-exp(-2aT))/(2a) - 2*(1-exp(-aT))/a)
//        + η²/b²*(T + (1-exp(-2bT))/(2b) - 2*(1-exp(-bT))/b)
//        + 2ρση/(a*b)*(T - (1-exp(-aT))/a - (1-exp(-bT))/b + (1-exp(-(a+b)T))/(a+b))
//
// 校准: P_market(0,T) = P_aff(0,T) * exp(-φ_T * T)
//   → φ_T = -log(P_market / P_aff) / T

class G2 {
public:
    G2(G2Params p, std::vector<Real> maturities, std::vector<Real> bond_prices)
        : params_(p), maturities_(std::move(maturities)), bonds_(std::move(bond_prices)) {
        validate();
    }
    G2(G2Params p) : params_(p) {
        // 无市场曲线时, 用平坦 φ(t)=r0 假设 (仅用于测试仿射部分)
        validate();
    }

    Real affine_B_a(Real T) const { return affine_B(params_.a, T); }
    Real affine_B_b(Real T) const { return affine_B(params_.b, T); }

    // 方差项 V(0,T) = Var[∫₀ᵀ (x(s)+y(s)) ds]  [Brigo-Mercurio (2006) eq 4.16]
    Real variance_term(Real T) const {
        Real a = params_.a, b = params_.b;
        Real sigma = params_.sigma, eta = params_.eta, rho = params_.rho;

        Real Vx = sigma * sigma / (a * a)
                  * (T + (1.0 - std::exp(-2.0 * a * T)) / (2.0 * a)
                     - 2.0 * (1.0 - std::exp(-a * T)) / a);
        Real Vy = eta * eta / (b * b)
                  * (T + (1.0 - std::exp(-2.0 * b * T)) / (2.0 * b)
                     - 2.0 * (1.0 - std::exp(-b * T)) / b);
        Real Vxy = 2.0 * rho * sigma * eta / (a * b)
                   * (T - (1.0 - std::exp(-a * T)) / a
                      - (1.0 - std::exp(-b * T)) / b
                      + (1.0 - std::exp(-(a + b) * T)) / (a + b));
        return Vx + Vy + Vxy;
    }

    // 仿射零息债 (未校准 φ)
    Real affine_zero_coupon_bond(Real T) const {
        if (T <= 0.0) return 1.0;
        Real B_a = affine_B_a(T);
        Real B_b = affine_B_b(T);
        Real V = variance_term(T);
        return std::exp(-B_a * params_.x0 - B_b * params_.y0 + 0.5 * V);
    }

    // 校准 φ(t) 使得 P_model(0,T) = P_market(0,T)
    // P_model = exp(-φ_T * T) * P_aff
    Real phi(Real T) const {
        if (maturities_.empty()) return 0.0;
        if (T <= 0.0) return 0.0;
        Real P_market = interpolate_bond(T);
        Real P_aff = affine_zero_coupon_bond(T);
        if (P_aff <= 0.0) return 0.0;
        return -std::log(P_market / P_aff) / T;
    }

    // 校准后零息债 (应等于市场曲线)
    Real zero_coupon_bond(Real T) const {
        if (T <= 0.0) return 1.0;
        if (maturities_.empty()) return affine_zero_coupon_bond(T);
        return interpolate_bond(T);
    }

    Real yield(Real T) const {
        if (T <= 0.0) return params_.x0 + params_.y0;
        return -std::log(zero_coupon_bond(T)) / T;
    }

    // 模拟两因子路径 (Euler)
    void simulate_path(Real T, Size n_steps, std::vector<Real>& path, Philox4x64& rng) const {
        path.resize(n_steps + 1);
        path[0] = params_.x0 + params_.y0;
        Real dt = T / static_cast<Real>(n_steps);
        Real sqrt_dt = std::sqrt(dt);
        Real x = params_.x0, y = params_.y0;
        for (Size i = 1; i <= n_steps; ++i) {
            // 生成相关正态
            Real u1 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
            if (u1 < 1e-300) u1 = 1e-300;
            Real u2 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
            Real z1 = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586476925286766559 * u2);
            u1 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
            if (u1 < 1e-300) u1 = 1e-300;
            u2 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
            Real z2_indep = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586476925286766559 * u2);
            Real z2 = params_.rho * z1 + std::sqrt(1.0 - params_.rho * params_.rho) * z2_indep;
            x = x - params_.a * x * dt + params_.sigma * sqrt_dt * z1;
            y = y - params_.b * y * dt + params_.eta * sqrt_dt * z2;
            // 加上 φ(t) (校准后)
            Real t_i = static_cast<Real>(i) * dt;
            Real phi_t = phi(t_i);
            path[i] = x + y + phi_t;
        }
    }

    const G2Params& params() const { return params_; }

private:
    G2Params params_;
    std::vector<Real> maturities_;
    std::vector<Real> bonds_;

    void validate() const {
        if (params_.a <= 0.0) throw std::invalid_argument("G2: a must be positive");
        if (params_.b <= 0.0) throw std::invalid_argument("G2: b must be positive");
        if (params_.sigma <= 0.0) throw std::invalid_argument("G2: sigma must be positive");
        if (params_.eta <= 0.0) throw std::invalid_argument("G2: eta must be positive");
        if (params_.rho < -1.0 || params_.rho > 1.0) {
            throw std::invalid_argument("G2: rho must be in [-1, 1]");
        }
        if (!maturities_.empty()) {
            if (maturities_.size() != bonds_.size()) {
                throw std::invalid_argument("G2: maturities and bonds size mismatch");
            }
            if (maturities_.size() < 2) {
                throw std::invalid_argument("G2: term structure requires at least 2 points");
            }
            for (Size i = 1; i < maturities_.size(); ++i) {
                if (maturities_[i] <= maturities_[i - 1]) {
                    throw std::invalid_argument("G2: maturities must be strictly increasing");
                }
            }
        }
    }

    Real interpolate_bond(Real T) const {
        if (maturities_.empty()) return affine_zero_coupon_bond(T);
        if (T <= maturities_.front()) {
            Real T0 = maturities_[0], T1 = maturities_[1];
            Real P0 = bonds_[0], P1 = bonds_[1];
            Real y0 = -std::log(P0) / T0;
            Real y1 = -std::log(P1) / T1;
            Real y = y0 + (y1 - y0) * (T - T0) / (T1 - T0);
            return std::exp(-y * T);
        }
        if (T >= maturities_.back()) {
            Real Tn = maturities_.back();
            Real Pn = bonds_.back();
            Real Tn_1 = maturities_[maturities_.size() - 2];
            Real Pn_1 = bonds_[bonds_.size() - 2];
            Real y_n = -std::log(Pn) / Tn;
            Real y_n_1 = -std::log(Pn_1) / Tn_1;
            Real y = y_n + (y_n - y_n_1) * (T - Tn) / (Tn - Tn_1);
            return std::exp(-y * T);
        }
        auto it = std::upper_bound(maturities_.begin(), maturities_.end(), T);
        Size idx = static_cast<Size>(it - maturities_.begin());
        Real T1 = maturities_[idx - 1];
        Real T2 = maturities_[idx];
        Real P1 = bonds_[idx - 1];
        Real P2 = bonds_[idx];
        Real w = (T - T1) / (T2 - T1);
        Real lnP = std::log(P1) + w * (std::log(P2) - std::log(P1));
        return std::exp(lnP);
    }
};

}  // namespace v1
}  // namespace cpphub
