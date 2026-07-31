#pragma once
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" Ch.4.5-4.6
// SOURCE: Glasserman, Heidelberger, Shahabuddin (1999) "Asymptotically optimal importance sampling
//         and stratification for pricing path-dependent options"
// 模块: 重要性抽样 (Importance Sampling, IS) — 通过 Girsanov 测度变换缩减方差
//
// 数学:
//   设定 (Z-shift 约定, 与 optimal_theta 公式一致):
//     - 原测度 P: Z ~ N(0,1), W_T = sqrt(T) * Z, S_T = S0 * exp((r-q-sigma^2/2)T + sigma*W_T)
//     - 新测度 Q: Z ~ N(theta, 1) (漂移调整), 通过 Girsanov 变换实现
//     - Radon-Nikodym 导数 (Q 关于 P, 在 Q-样本 z=Z+theta 处):
//         L = dQ/dP = exp(theta*Z + 0.5*theta^2)
//       (Z 为 N(0,1) 抽样, Q-样本 z = Z + theta)
//     - IS 估计量: X_IS = X * L^{-1} = X * exp(-theta*Z - 0.5*theta^2)
//       使得 E_P[X] = E_Q[X * L^{-1}] (无偏)
//     - 最优 theta (BSM 欧式 Call/Put, 使 S_T 集中在 K 附近):
//         theta* = (ln(K/S0) - (r-q-sigma^2/2)*T) / (sigma*sqrt(T))
//       即把 Z 的均值偏移 -d2 (BSM d2), 让 S_T 在 K 附近高概率取值
//
//   稀有事件优势: 对深度 OTM 期权, 标准 MC 大部分样本 payoff=0, 方差大;
//   IS 将样本推向 K 附近, 显著增加有效样本数, 方差缩减可达数十倍

#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cpphub {
inline namespace v1 {

// ============ IS 配置 ============
struct ISConfig {
    Real theta = 0.0;             // Girsanov 漂移调整 (Z-shift 约定)
    bool auto_optimize = true;    // 自动计算最优 theta (仅 BSM 欧式 Call/Put)
};

// ============ Importance Sampling 引擎 ============
class ImportanceSampling {
public:
    explicit ImportanceSampling(ISConfig cfg)
        : cfg_(cfg), theta_(cfg.theta) {
        if (std::isnan(theta_)) {
            // NaN 表示自动优化 (但需要参数, 在 price_european_call 中设置)
            theta_ = 0.0;
            auto_optimize_pending_ = true;
        } else {
            auto_optimize_pending_ = false;
        }
    }

    // 采样一条路径的漂移调整 W_T (返回 W_T under Q, 以及似然比 L = dQ/dP)
    // 输入 Z ~ N(0,1) (P 下的标准正态), 输出 W_T = sqrt(T)*(Z + theta) (Q 下的终端 Brownian)
    // likelihood_ratio = L = exp(theta*Z + 0.5*theta^2)
    // IS 估计: X_IS = X * L^{-1} = X * exp(-theta*Z - 0.5*theta^2)
    Real sample_shifted_WT(Real T, Real Z, Real& likelihood_ratio) const {
        Real sqrtT = std::sqrt(T);
        Real W_T = sqrtT * (Z + theta_);
        // L = dQ/dP 在 Q-样本 z = Z+theta 处: exp(theta*Z + 0.5*theta^2)
        // (Z 是 P-N(0,1), Q-样本为 Z+theta; 由密度比推导得 dP/dQ = exp(-theta*Z - 0.5*theta^2))
        likelihood_ratio = std::exp(theta_ * Z + 0.5 * theta_ * theta_);
        return W_T;
    }

    // 最优 theta (BSM 欧式 Call): theta* = (ln(K/S0) - (r-q-sigma^2/2)T) / (sigma*sqrt(T)) = -d2
    static Real optimal_theta_call(Real S0, Real K, Real T, Real r, Real q, Real sigma) {
        if (T <= 0.0 || sigma <= 0.0 || S0 <= 0.0 || K <= 0.0)
            throw std::invalid_argument("optimal_theta_call: invalid parameters");
        Real mu_T = (r - q - 0.5 * sigma * sigma) * T;
        return (std::log(K / S0) - mu_T) / (sigma * std::sqrt(T));
    }

    // 最优 theta (BSM 欧式 Put): 同 Call 公式 (使 S_T 集中在 K 附近)
    // 对 OTM Put (K < S0), theta* < 0, 将 S_T 向下推至 K
    static Real optimal_theta_put(Real S0, Real K, Real T, Real r, Real q, Real sigma) {
        return optimal_theta_call(S0, K, T, r, q, sigma);
    }

    // IS 结果
    struct ISResult {
        Real price = 0.0;                  // 折现后期权价格
        Real std_error = 0.0;              // 标准误差
        Real variance_reduction_ratio = 1.0;  // 方差缩减比 (vs 标准 MC)
    };

    // 批量 IS 估计: BSM 欧式 Call (使用 Girsanov 测度变换)
    ISResult price_european_call(Real S0, Real K, Real T, Real r, Real q, Real sigma,
                                  Size n_paths, uint64_t seed) const {
        if (n_paths == 0) throw std::invalid_argument("IS: n_paths must be positive");
        if (T <= 0.0) {
            ISResult res;
            res.price = std::max(S0 - K, 0.0);
            return res;
        }

        // 解析最优 theta (若 auto_optimize)
        Real theta_use = theta_;
        if (cfg_.auto_optimize || auto_optimize_pending_) {
            theta_use = optimal_theta_call(S0, K, T, r, q, sigma);
        }

        Real mu_T = (r - q - 0.5 * sigma * sigma) * T;
        Real sqrtT = std::sqrt(T);
        Real df = std::exp(-r * T);

        // IS 估计: E_P[payoff] = E_Q[payoff(S_T^Q) * exp(-theta*Z - 0.5*theta^2)]
        // 其中 S_T^Q = S0 * exp(mu_T + sigma*sqrt(T)*(Z + theta))
        Real sum_x = 0.0, sum_xx = 0.0;
        // 标准 MC 估计 (用于方差比较)
        Real sum_x_mc = 0.0, sum_xx_mc = 0.0;

        for (Size i = 0; i < n_paths; ++i) {
            Philox4x64 rng(seed, static_cast<uint64_t>(i));
            uint64_t r1 = rng();
            uint64_t r2 = rng();
            Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            Real u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            Real Z = box_muller(u1, u2).first;

            // IS 路径: S_T^Q = S0 * exp(mu_T + sigma*sqrt(T)*(Z + theta))
            Real W_T_Q = sqrtT * (Z + theta_use);
            Real S_T_Q = S0 * std::exp(mu_T + sigma * W_T_Q);
            Real payoff_is = std::max(S_T_Q - K, 0.0);
            // IS 权重 (dP/dQ): exp(-theta*Z - 0.5*theta^2)
            Real weight = std::exp(-theta_use * Z - 0.5 * theta_use * theta_use);
            Real X_is = payoff_is * weight;

            // 标准 MC 路径 (P 下): S_T = S0 * exp(mu_T + sigma*sqrt(T)*Z)
            Real S_T_P = S0 * std::exp(mu_T + sigma * sqrtT * Z);
            Real X_mc = std::max(S_T_P - K, 0.0);

            sum_x += X_is;
            sum_xx += X_is * X_is;
            sum_x_mc += X_mc;
            sum_xx_mc += X_mc * X_mc;
        }

        Real inv_n = 1.0 / static_cast<Real>(n_paths);
        Real mean_is = sum_x * inv_n;
        Real var_is = sum_xx * inv_n - mean_is * mean_is;
        if (var_is < 0.0) var_is = 0.0;

        Real mean_mc = sum_x_mc * inv_n;
        Real var_mc = sum_xx_mc * inv_n - mean_mc * mean_mc;
        if (var_mc < 0.0) var_mc = 0.0;

        ISResult res;
        res.price = df * mean_is;
        res.std_error = df * std::sqrt(var_is * inv_n);
        if (var_is > 1e-30) {
            res.variance_reduction_ratio = var_mc / var_is;
        } else {
            res.variance_reduction_ratio = std::numeric_limits<Real>::max();
        }
        return res;
    }

    ISConfig config() const { return cfg_; }

    Real theta() const { return theta_; }

private:
    ISConfig cfg_;
    Real theta_;
    bool auto_optimize_pending_;
};

}  // namespace v1
}  // namespace cpphub
