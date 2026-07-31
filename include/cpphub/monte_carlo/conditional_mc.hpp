#pragma once
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" Ch.4.3
// SOURCE: Curran (1994) "Valuing Asian and Portfolio Options by Conditioning on the Geometric Mean Price"
// SOURCE: Hull (2018) Ch.26 (Path-Dependent Derivatives, Barrier options)
// 模块: 条件蒙特卡洛 (Conditional Monte Carlo, CMC) — 通过条件期望缩减方差
//
// 数学:
//   CMC 原理: 设 X 为待估随机变量, Y 为条件变量. 则
//     E[X] = E[E[X | Y]],  Var(X) = Var(E[X|Y]) + E[Var(X|Y)]
//   用 E[X|Y] 替代 X 作为估计量, 消除 E[Var(X|Y)] 项, 方差严格不增.
//
//   1. 障碍期权 CMC (条件化 W_T, 用反射原理计算跨越概率):
//      - Up-and-Out Call, 障碍 B (S 空间), 终端 Brownian W_T = w
//      - 反射原理: P(max_{0<=t<=T} W_t >= b | W_T=w) = exp(-2*b*(b-w)/T),  w < b
//      - 生存概率: P(survive | W_T=w) = 1 - exp(-2*b*(b-w)/T),  w < b; 0, w >= b
//      - 其中 b 为 W 空间障碍 (近似常数: b = ln(B/S0)/sigma)
//      - 条件期望 payoff = max(S_T - K, 0) * survival_prob
//        其中 S_T = S0 * exp((r-q-0.5*sigma^2)*T + sigma*w)
//
//   2. 亚式期权 CMC (条件化 W_T, 用对数正态条件期望, Curran 1994 思路):
//      - A_T = (1/n) * sum_{i=1}^n S(t_i),  S(t_i) = S0*exp(mu*t_i + sigma*W(t_i))
//      - Brownian bridge: W(t_i) | W_T=w ~ N((t_i/T)*w, t_i*(T-t_i)/T)
//      - 条件均值 m = E[A_T | W_T=w] = (1/n) * sum E[S(t_i)|W_T=w]
//      - 条件方差 v = Var(A_T|W_T=w) 用 Brownian bridge 协方差计算
//      - 用对数正态近似 A_T | W_T, 计算 Call 价格 (moment matching)

#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/math.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cpphub {
inline namespace v1 {

// ============ Conditional MC 工具集 ============
class ConditionalMC {
public:
    // 障碍期权 CMC (条件化 W_T, 用 Bachelier-Lévy 公式计算线性障碍跨越概率)
    // Up-and-Out Call: 障碍 B > S0, B > K (典型情形)
    // 输入: W_T (条件化的终端 Brownian), 输出 survival_prob (条件生存概率)
    // 返回: 条件期望 payoff = max(S_T - K, 0) * survival_prob (未折现)
    //
    // 数学: S_t = S0*exp(mu*t + sigma*W_t), 障碍 B 在 W 空间为线性时变:
    //       b(t) = (ln(B/S0) - mu*t)/sigma = b_0 - (mu/sigma)*t
    //       其中 b_0 = ln(B/S0)/sigma, b_T = b_0 - mu*T/sigma
    // Bachelier-Lévy 公式 (Brownian bridge 跨越线性边界):
    //   P(max_{t<=T} W_t >= b(t) | W_T=w) = exp(-2*b_0*(b_T - w)/T),  w < b_T
    //   survival_prob = 1 - exp(-2*b_0*(b_T - w)/T),  w < b_T; 0, w >= b_T
    // 注意: 此公式对常数 S-空间障碍是精确的 (连续监控); 离散监控会有微小偏差.
    static Real barrier_up_out_call_cmc(Real S0, Real K, Real B, Real T,
                                         Real r, Real q, Real sigma,
                                         Real W_T,
                                         Real& survival_prob) {
        if (T <= 0.0) throw std::invalid_argument("barrier_cmc: T must be positive");
        if (sigma <= 0.0) throw std::invalid_argument("barrier_cmc: sigma must be positive");
        if (S0 <= 0.0 || K <= 0.0 || B <= 0.0)
            throw std::invalid_argument("barrier_cmc: S0, K, B must be positive");

        Real mu = r - q - 0.5 * sigma * sigma;
        Real mu_T = mu * T;
        Real S_T = S0 * std::exp(mu_T + sigma * W_T);

        // W-空间线性障碍: b(t) = b_0 - (mu/sigma)*t
        Real b_0 = std::log(B / S0) / sigma;             // t=0 障碍
        Real b_T = (std::log(B / S0) - mu_T) / sigma;    // t=T 障碍

        // Bachelier-Lévy: 跨越概率 = exp(-2*b_0*(b_T - w)/T),  w < b_T
        if (W_T >= b_T) {
            survival_prob = 0.0;
            return 0.0;  // 终端已经越过 t=T 障碍, 必然跨越
        }
        Real exponent = -2.0 * b_0 * (b_T - W_T) / T;
        if (exponent < -700.0) {
            survival_prob = 1.0;  // 跨越概率下溢为 0, 生存概率 = 1
        } else {
            survival_prob = 1.0 - std::exp(exponent);
        }

        Real payoff = std::max(S_T - K, 0.0);
        return payoff * survival_prob;
    }

    // 亚式期权 CMC (条件化 W_T, 用对数正态条件期望近似)
    // A_T = (1/n) * sum_{i=1}^n S(t_i), t_i = i*T/n
    // 输出 conditional_mean = E[A_T | W_T] (exact)
    // 返回: E[max(A_T - K, 0) | W_T] (lognormal moment matching 近似, 未折现)
    static Real asian_arithmetic_call_cmc(Real S0, Real K, Real T,
                                           Real r, Real q, Real sigma,
                                           Size n_steps, Real W_T,
                                           Real& conditional_mean) {
        if (T <= 0.0) throw std::invalid_argument("asian_cmc: T must be positive");
        if (sigma <= 0.0) throw std::invalid_argument("asian_cmc: sigma must be positive");
        if (S0 <= 0.0 || K <= 0.0) throw std::invalid_argument("asian_cmc: S0, K must be positive");
        if (n_steps == 0) throw std::invalid_argument("asian_cmc: n_steps must be positive");

        Real mu = r - q - 0.5 * sigma * sigma;
        Real dt = T / static_cast<Real>(n_steps);
        Real inv_n = 1.0 / static_cast<Real>(n_steps);

        // 计算每个 S(t_i) 的条件期望 E[S(t_i) | W_T=w]
        // W(t_i) | W_T=w ~ N((t_i/T)*w, t_i*(T-t_i)/T)
        // E[S(t_i)|W_T=w] = S0 * exp(mu*t_i + sigma*(t_i/T)*w + 0.5*sigma^2*t_i*(1-t_i/T))
        std::vector<Real> cond_E_S(n_steps);
        conditional_mean = 0.0;
        for (Size i = 0; i < n_steps; ++i) {
            Real t_i = static_cast<Real>(i + 1) * dt;
            Real mean_W = (t_i / T) * W_T;
            Real var_W = t_i * (T - t_i) / T;
            cond_E_S[i] = S0 * std::exp(mu * t_i + sigma * mean_W + 0.5 * sigma * sigma * var_W);
            conditional_mean += cond_E_S[i];
        }
        conditional_mean *= inv_n;

        // 计算条件方差 Var(A_T | W_T=w)
        // Cov(W(t_i), W(t_j) | W_T) = min(t_i, t_j) - t_i*t_j/T  (Brownian bridge)
        // Cov(S_i, S_j | W_T) = E[S_i|W_T]*E[S_j|W_T] * (exp(sigma^2 * Cov(W_i,W_j|W_T)) - 1)
        Real cond_var = 0.0;
        for (Size i = 0; i < n_steps; ++i) {
            Real t_i = static_cast<Real>(i + 1) * dt;
            for (Size j = 0; j < n_steps; ++j) {
                Real t_j = static_cast<Real>(j + 1) * dt;
                Real cov_W = std::min(t_i, t_j) - t_i * t_j / T;
                cond_var += cond_E_S[i] * cond_E_S[j]
                            * (std::exp(sigma * sigma * cov_W) - 1.0);
            }
        }
        cond_var *= inv_n * inv_n;
        if (cond_var < 0.0) cond_var = 0.0;

        // 用对数正态近似 A_T | W_T: moment matching
        // 若 A ~ Lognormal 使得 E[A]=m, Var(A)=v, 则
        //   sigma_L^2 = ln(1 + v/m^2),  mu_L = ln(m) - 0.5*sigma_L^2
        //   E[max(A-K,0)] = m*N(d1) - K*N(d2)
        //   d1 = (ln(m/K) + 0.5*sigma_L^2) / sigma_L,  d2 = d1 - sigma_L
        if (conditional_mean <= 0.0 || cond_var <= 0.0) {
            return std::max(conditional_mean - K, 0.0);
        }
        Real sigma_L_sq = std::log(1.0 + cond_var / (conditional_mean * conditional_mean));
        if (sigma_L_sq <= 0.0) {
            return std::max(conditional_mean - K, 0.0);
        }
        Real sigma_L = std::sqrt(sigma_L_sq);
        Real d1 = (std::log(conditional_mean / K) + 0.5 * sigma_L_sq) / sigma_L;
        Real d2 = d1 - sigma_L;
        return conditional_mean * normal_cdf(d1) - K * normal_cdf(d2);
    }
};

}  // namespace v1
}  // namespace cpphub
