#pragma once
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" Ch.4.4
// SOURCE: Neyman (1934) "On the two different aspects of the representative method"
// 模块: 分层抽样 (Stratified Sampling) — 通过将样本空间分层缩减方差
//
// 数学:
//   设将 [0,1] 分为 K 层, 第 k 层区间 I_k = [(k-1)/K, k/K), 概率 p_k = 1/K.
//   从第 k 层抽 n_k 个样本, 样本均值 x_bar_k.
//   分层估计量: X_strat = sum_k p_k * x_bar_k  (无偏)
//   分层方差: Var(X_strat) = sum_k p_k^2 * sigma_k^2 / n_k
//
//   分配策略:
//     - 等比例分配: n_k = n * p_k = n / K  (每层样本数相同)
//     - Neyman 最优分配: n_k ∝ p_k * sigma_k, 即方差大的层多采样
//       n_k = n * (p_k * sigma_k) / sum_j (p_j * sigma_j)
//       最优分配下 Var(X_strat) = (sum_k p_k * sigma_k)^2 / n  (Cauchy-Schwarz)
//
//   与朴素 MC 对比: 朴素 MC 方差 = Var(X) / n, 分层方差 ≤ Var(X) / n (always).
//   分层样本生成: 第 k 层第 j 个样本 U_kj = (k + V_kj) / K, V_kj ~ U(0,1)

#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"
#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ 分层抽样配置 ============
struct StratifiedConfig {
    Size n_strata = 10;                // 分层数 K
    bool neyman_allocation = true;     // Neyman 最优分配 (否则等比例)
    Size pilot_samples = 100;          // Neyman 预采样数 (估计各层方差)
};

// ============ Stratified Sampling 引擎 ============
class StratifiedSampling {
public:
    explicit StratifiedSampling(StratifiedConfig cfg)
        : cfg_(cfg) {
        if (cfg_.n_strata == 0)
            throw std::invalid_argument("StratifiedSampling: n_strata must be positive");
        if (cfg_.pilot_samples == 0)
            throw std::invalid_argument("StratifiedSampling: pilot_samples must be positive");
    }

    // 生成分层采样的 U(0,1) 样本 (返回 n_paths 个样本)
    // 等比例分配: 每层 n_paths / K 个样本 (余数分配到前几层)
    // 每个样本 u_ij = (i + V_ij) / K, i=0..K-1, j=1..n_i
    std::vector<Real> generate_stratified_uniforms(Size n_paths, Philox4x64& rng) const {
        if (n_paths == 0) return {};
        Size K = cfg_.n_strata;
        std::vector<Real> uniforms;
        uniforms.reserve(n_paths);

        // 等比例分配: 每层 base = n_paths / K, 余数 r 个层各 +1
        Size base = n_paths / K;
        Size remainder = n_paths % K;
        for (Size i = 0; i < K; ++i) {
            Size n_i = base + (i < remainder ? 1 : 0);
            for (Size j = 0; j < n_i; ++j) {
                // V_ij ~ U(0,1) from rng
                uint64_t r = rng();
                Real V = (r >> 11) * (1.0 / 9007199254740992.0);
                // 防止 V = 0 导致 inv_normal_cdf(-inf); clip 到 (0,1)
                if (V < 1e-15) V = 1e-15;
                if (V > 1.0 - 1e-15) V = 1.0 - 1e-15;
                Real u = (static_cast<Real>(i) + V) / static_cast<Real>(K);
                uniforms.push_back(u);
            }
        }
        return uniforms;
    }

    // 分层结果
    struct StratResult {
        Real price = 0.0;                  // 折现后期权价格
        Real std_error = 0.0;              // 标准误差
        Real variance_reduction_ratio = 1.0;  // 方差缩减比 (vs 标准 MC)
    };

    // 分层 MC 估计: BSM 欧式 Call
    // W_T = sqrt(T) * Phi^{-1}(u), 然后标准 BSM payoff
    StratResult price_european_call_stratified(Real S0, Real K, Real T, Real r, Real q,
                                                 Real sigma, Size n_paths, uint64_t seed) const {
        if (n_paths == 0) throw std::invalid_argument("Strat: n_paths must be positive");
        if (T <= 0.0) {
            StratResult res;
            res.price = std::max(S0 - K, 0.0);
            return res;
        }

        Real mu_T = (r - q - 0.5 * sigma * sigma) * T;
        Real sqrtT = std::sqrt(T);
        Real df = std::exp(-r * T);
        Size K_strata = cfg_.n_strata;

        std::vector<Size> n_k;
        if (cfg_.neyman_allocation) {
            n_k = neyman_allocation(S0, K, T, r, q, sigma, n_paths, seed);
        } else {
            n_k.assign(K_strata, n_paths / K_strata);
            Size rem = n_paths % K_strata;
            for (Size i = 0; i < rem; ++i) n_k[i] += 1;
        }

        // 在每层内生成样本, 计算层内均值和方差
        // 同时计算标准 MC 方差 (合并所有样本的样本方差)
        Real price_sum = 0.0;          // Σ p_k * x̄_k = Σ (n_k/n) * x̄_k
        Real var_sum = 0.0;            // Σ p_k^2 * s_k^2 / n_k (stratified 方差)
        Real all_sum = 0.0, all_sum_sq = 0.0;  // 标准 MC 方差用
        Size actual_total = 0;

        // 用独立的 rng 流以避免与 pilot 重叠
        Philox4x64 main_rng(seed, 999983);

        for (Size k = 0; k < K_strata; ++k) {
            if (n_k[k] == 0) continue;
            Real sum_k = 0.0, sum_sq_k = 0.0;
            for (Size j = 0; j < n_k[k]; ++j) {
                uint64_t r1 = main_rng();
                Real V = (r1 >> 11) * (1.0 / 9007199254740992.0);
                if (V < 1e-15) V = 1e-15;
                if (V > 1.0 - 1e-15) V = 1.0 - 1e-15;
                Real u = (static_cast<Real>(k) + V) / static_cast<Real>(K_strata);
                Real Z = inv_normal_cdf(u);
                Real W_T = sqrtT * Z;
                Real S_T = S0 * std::exp(mu_T + sigma * W_T);
                Real payoff = std::max(S_T - K, 0.0);
                sum_k += payoff;
                sum_sq_k += payoff * payoff;
                all_sum += payoff;
                all_sum_sq += payoff * payoff;
                ++actual_total;
            }
            Real mean_k = sum_k / static_cast<Real>(n_k[k]);
            Real var_k = (n_k[k] >= 2)
                          ? (sum_sq_k - sum_k * sum_k / static_cast<Real>(n_k[k]))
                              / static_cast<Real>(n_k[k] - 1)
                          : 0.0;
            if (var_k < 0.0) var_k = 0.0;
            Real p_k = static_cast<Real>(n_k[k]) / static_cast<Real>(n_paths);
            price_sum += p_k * mean_k;
            var_sum += p_k * p_k * var_k / static_cast<Real>(n_k[k]);
        }

        // 标准 MC 方差估计 (基于全部样本)
        Real inv_total = 1.0 / static_cast<Real>(actual_total);
        Real mc_mean = all_sum * inv_total;
        Real mc_var = (actual_total >= 2)
                       ? (all_sum_sq - all_sum * all_sum * inv_total)
                           / static_cast<Real>(actual_total - 1)
                       : 0.0;
        if (mc_var < 0.0) mc_var = 0.0;
        // 标准 MC 估计量方差 = mc_var / n
        Real mc_est_var = mc_var * inv_total;

        StratResult res;
        res.price = df * price_sum;
        res.std_error = df * std::sqrt(var_sum);
        if (var_sum > 1e-30) {
            res.variance_reduction_ratio = mc_est_var / var_sum;
        } else {
            res.variance_reduction_ratio = std::numeric_limits<Real>::max();
        }
        return res;
    }

    // Neyman 分配 (预采样估计各层方差, 返回各层样本数)
    // 总样本数 = n_paths, 各层分配 n_k ∝ p_k * sigma_k (p_k = 1/K 等概率分层)
    std::vector<Size> neyman_allocation(Real S0, Real K, Real T, Real r, Real q,
                                         Real sigma, Size n_paths, uint64_t seed) const {
        Size K_strata = cfg_.n_strata;
        Size pilot = cfg_.pilot_samples;
        std::vector<Size> n_k(K_strata, 0);
        if (n_paths == 0) return n_k;

        Real mu_T = (r - q - 0.5 * sigma * sigma) * T;
        Real sqrtT = std::sqrt(T);

        // 用独立 stream 做 pilot, 不影响主采样
        Philox4x64 pilot_rng(seed, 1234567);

        std::vector<Real> sigma_k(K_strata, 0.0);
        for (Size k = 0; k < K_strata; ++k) {
            Real sum = 0.0, sum_sq = 0.0;
            for (Size j = 0; j < pilot; ++j) {
                uint64_t r1 = pilot_rng();
                Real V = (r1 >> 11) * (1.0 / 9007199254740992.0);
                if (V < 1e-15) V = 1e-15;
                if (V > 1.0 - 1e-15) V = 1.0 - 1e-15;
                Real u = (static_cast<Real>(k) + V) / static_cast<Real>(K_strata);
                Real Z = inv_normal_cdf(u);
                Real W_T = sqrtT * Z;
                Real S_T = S0 * std::exp(mu_T + sigma * W_T);
                Real payoff = std::max(S_T - K, 0.0);
                sum += payoff;
                sum_sq += payoff * payoff;
            }
            Real mean = sum / static_cast<Real>(pilot);
            Real var = (pilot >= 2)
                        ? (sum_sq - sum * sum / static_cast<Real>(pilot))
                            / static_cast<Real>(pilot - 1)
                        : 0.0;
            if (var < 0.0) var = 0.0;
            sigma_k[k] = std::sqrt(var);
        }

        // Neyman: n_k = n * (p_k * sigma_k) / sum(p_j * sigma_j), p_k = 1/K
        Real total_weight = 0.0;
        for (Size k = 0; k < K_strata; ++k) total_weight += sigma_k[k];
        if (total_weight <= 0.0) {
            // 全部 sigma_k = 0 (退化), 等比例分配
            Size base = n_paths / K_strata;
            Size rem = n_paths % K_strata;
            for (Size k = 0; k < K_strata; ++k) {
                n_k[k] = base + (k < rem ? 1 : 0);
            }
            return n_k;
        }

        // 分配: 比例后取整, 余数按四舍五入分配给最大权重层
        std::vector<Real> ideal(K_strata);
        std::vector<Real> frac(K_strata);
        Size allocated = 0;
        for (Size k = 0; k < K_strata; ++k) {
            ideal[k] = static_cast<Real>(n_paths) * sigma_k[k] / total_weight;
            n_k[k] = static_cast<Size>(std::floor(ideal[k]));
            frac[k] = ideal[k] - static_cast<Real>(n_k[k]);
            allocated += n_k[k];
        }
        // 剩余样本按 frac 降序分配
        Size remaining = n_paths - allocated;
        std::vector<Size> order(K_strata);
        for (Size k = 0; k < K_strata; ++k) order[k] = k;
        std::sort(order.begin(), order.end(),
                   [&](Size a, Size b) { return frac[a] > frac[b]; });
        for (Size i = 0; i < remaining; ++i) {
            n_k[order[i % K_strata]] += 1;
        }
        // 保证至少 1 个样本 (若某层 sigma_k = 0 且 n_paths < K_strata)
        for (Size k = 0; k < K_strata; ++k) {
            if (n_k[k] == 0 && n_paths >= K_strata) {
                // 从最多的层挪 1 个过来
                Size max_idx = 0;
                for (Size kk = 1; kk < K_strata; ++kk) {
                    if (n_k[kk] > n_k[max_idx]) max_idx = kk;
                }
                if (n_k[max_idx] > 1) {
                    n_k[max_idx] -= 1;
                    n_k[k] = 1;
                }
            }
        }
        return n_k;
    }

    StratifiedConfig config() const { return cfg_; }

private:
    StratifiedConfig cfg_;
};

}  // namespace v1
}  // namespace cpphub
