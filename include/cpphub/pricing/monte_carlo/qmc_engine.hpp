#pragma once
// SOURCE: Caflisch, Morokoff, Owen (1997) "Valuation of mortgage-backed securities
//         using Brownian bridges to reduce effective dimension" J. Comput. Finance 1, 27-46.
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" Ch.5
// SOURCE: L'Ecuyer, Lemieux (2002) "A survey of randomized quasi-Monte Carlo methods"
//         in "Uniform Distribution Theory" Ch.12
// 模块: Quasi-Monte Carlo (QMC) 引擎 — Sobol + Brownian Bridge + 误差估计
//
// ==================== QMC 数学框架 ====================
//
// 标准 MC: 估计 E[f(X)] ≈ (1/N) Σ f(X_i), X_i i.i.d.
//   收敛速度: O(N^{-1/2}) (中心极限定理)
//
// Quasi-MC: 用低偏差序列 (Sobol) 替代伪随机数
//   估计 E[f(X)] ≈ (1/N) Σ f(X_i), X_i 为 Sobol 点 (确定性)
//   收敛速度: O(N^{-1} (log N)^d) (Koksma-Hlawka 不等式)
//   对低有效维度问题, QMC 比 MC 快 ~10x-1000x
//
// Randomized QMC (RQMC): 用 Owen 随机化 Sobol 估计误差
//   生成 R 个独立 scrambled Sobol 序列, 每个长度 N
//   估计 θ̂ = (1/R) Σ_r θ̂_r, 其中 θ̂_r = (1/N) Σ_i f(X_{r,i})
//   标准误差: SE = std(θ̂_r) / sqrt(R)
//   保持 QMC 收敛速度, 同时提供无偏误差估计
//
// Brownian Bridge 路径构造 (与 Sobol 配合):
//   标准 MC 路径生成: W(t_i) = W(t_{i-1}) + sqrt(dt) * Z_i (维度同等重要)
//   BB 路径生成: 先 W(T), 再递归中点 (最重要的维度控制 W(T))
//   有效维度从 n 降到 ~log2(n), 与 Sobol 配合显著提升收敛速度
//
// ==================== 与现有 MC 引擎集成 ====================
//
// QMC 引擎复用:
//   - MultiAssetGBMPathGenerator (路径生成, GBM Exact scheme)
//   - PathPayoff / MultiAssetPayoff (payoff 函数)
//   - MCResult / MCConfig (结果与配置)
//
// 新增:
//   - QMCConfig (n_paths, n_replicates for RQMC, use_bb, scramble_seed)
//   - price_european_qmc (欧式期权, 维度 = n_assets)
//   - price_path_dependent_qmc (路径相关, 维度 = n_steps * n_assets)

#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/pricing/monte_carlo/sobol.hpp"
#include "cpphub/pricing/monte_carlo/brownian_bridge.hpp"
#include "cpphub/pricing/monte_carlo/mc_engine.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/pricing/monte_carlo/path_dependent_payoffs.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_payoffs.hpp"
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ QMC 配置 ============
struct QMCConfig {
    Size n_paths = 4096;            // 每个 Sobol 序列的点数 (建议 2^k)
    Size n_replicates = 16;         // RQMC 独立重复次数 (用于误差估计, ≥1)
    uint64_t base_seed = 42;        // 基础 seed (每个 replicate 用不同 scramble_seed)
    bool use_brownian_bridge = true; // 是否使用 Brownian Bridge 路径构造
    Real df = 1.0;                  // 折现因子 exp(-rT)
    bool skip_first = true;         // 跳过 n=0 (Sobol 第一个点是原点)
};

// ============ QMC 结果 (复用 MCResult + 额外信息) ============
struct QMCResult {
    Real price = 0.0;          // QMC 估计 (所有 replicate 的均值)
    Real std_error = 0.0;      // 标准误差 (RQMC, replicate 间样本标准差 / sqrt(R))
    Real ci_lower = 0.0;       // 95% 置信区间下界
    Real ci_upper = 0.0;       // 95% 置信区间上界
    Size n_paths = 0;          // 每个 replicate 的路径数
    Size n_replicates = 0;     // RQMC 重复次数
    Size n_total_paths = 0;    // 总路径数 = n_paths * n_replicates
    Real variance_reduction = 0.0;  // vs 标准 MC 的方差缩减 (可选, 0 表示未计算)
    std::vector<Real> replicate_estimates;  // 各 replicate 的独立估计 (for 诊断)
};

// ============ 单资产欧式期权 QMC ============
// 维度 = 1 (只需 W(T))
// uniforms[0] → W(T) = sqrt(T) * inv_normal_cdf(uniforms[0])
// S(T) = S0 * exp((r-q-0.5σ²)T + σW(T))
inline QMCResult price_european_qmc_single(
        Real S0, Real sigma, Real r, Real q, Real T,
        const std::function<Real(Real ST)>& payoff,  // 接收 S(T), 返回未折现 payoff
        const QMCConfig& cfg = {}) {
    if (S0 <= 0.0) throw std::invalid_argument("QMC: S0 must be positive");
    if (sigma < 0.0) throw std::invalid_argument("QMC: sigma must be non-negative");
    if (T <= 0.0) throw std::invalid_argument("QMC: T must be positive");
    if (cfg.n_paths == 0) throw std::invalid_argument("QMC: n_paths must be positive");
    if (cfg.n_replicates == 0) throw std::invalid_argument("QMC: n_replicates must be positive");

    Real drift_T = (r - q - 0.5 * sigma * sigma) * T;
    Real sqrt_T = std::sqrt(T);
    // 注意: W(T) = sqrt(T)*Z 已含 sqrt(T) 因子, 故 vol 系数为 σ (而非 σ*sqrt(T))
    // S(T) = S0 * exp(drift_T + σ * W(T))

    QMCResult result;
    result.n_paths = cfg.n_paths;
    result.n_replicates = cfg.n_replicates;
    result.n_total_paths = cfg.n_paths * cfg.n_replicates;
    result.replicate_estimates.reserve(cfg.n_replicates);

    uint64_t start_n = cfg.skip_first ? 1 : 0;

    for (Size rep = 0; rep < cfg.n_replicates; ++rep) {
        uint64_t scramble_seed = (cfg.base_seed == 0 && rep == 0)
            ? 0  // 第一个 replicate 可以用未 scrambled 版本
            : cfg.base_seed + static_cast<uint64_t>(rep) * 0x9E3779B97F4A7C15ULL + 1;

        SobolSequence sobol(1, scramble_seed);
        Real sum = 0.0;
        for (uint64_t i = 0; i < cfg.n_paths; ++i) {
            auto u = sobol(start_n + i);
            Real u0 = u[0];
            if (u0 <= 0.0) u0 = 1e-15;
            if (u0 >= 1.0) u0 = 1.0 - 1e-15;
            Real WT = sqrt_T * inv_normal_cdf(u0);
            Real ST = S0 * std::exp(drift_T + sigma * WT);
            sum += payoff(ST);
        }
        Real est = sum / static_cast<Real>(cfg.n_paths);
        result.replicate_estimates.push_back(est);
    }

    // RQMC 误差估计
    Real mean_est = std::accumulate(result.replicate_estimates.begin(),
                                      result.replicate_estimates.end(), 0.0)
                    / static_cast<Real>(cfg.n_replicates);
    result.price = mean_est * cfg.df;

    if (cfg.n_replicates >= 2) {
        Real var = 0.0;
        for (Real e : result.replicate_estimates) {
            var += (e - mean_est) * (e - mean_est);
        }
        var /= static_cast<Real>(cfg.n_replicates - 1);
        result.std_error = std::sqrt(var / static_cast<Real>(cfg.n_replicates)) * cfg.df;
    }
    result.ci_lower = result.price - 1.96 * result.std_error;
    result.ci_upper = result.price + 1.96 * result.std_error;
    return result;
}

// ============ 单资产路径相关 QMC (Brownian Bridge) ============
// 维度 = n_steps
// Sobol uniforms[n_steps] → Brownian Bridge → W(t_1)..W(t_n) → S(t_1)..S(t_n)
inline QMCResult price_path_dependent_qmc_single(
        Real S0, Real sigma, Real r, Real q, Real T, Size n_steps,
        const PathPayoff& payoff,  // 接收 path[0..n_steps-1] = S(t_1)..S(t_n), 返回未折现 payoff
        const QMCConfig& cfg = {}) {
    if (S0 <= 0.0) throw std::invalid_argument("QMC: S0 must be positive");
    if (sigma < 0.0) throw std::invalid_argument("QMC: sigma must be non-negative");
    if (T <= 0.0) throw std::invalid_argument("QMC: T must be positive");
    if (n_steps == 0) throw std::invalid_argument("QMC: n_steps must be positive");
    if (!payoff) throw std::invalid_argument("QMC: payoff function is null");

    Real dt = T / static_cast<Real>(n_steps);
    Real drift_dt = (r - q - 0.5 * sigma * sigma) * dt;
    // dW = W(t+dt) - W(t) ~ N(0, dt) 已含 sqrt(dt) 因子, 故 vol 系数为 σ
    // S(t+dt) = S(t) * exp(drift_dt + σ * dW)

    BrownianBridge bb(n_steps, T);

    QMCResult result;
    result.n_paths = cfg.n_paths;
    result.n_replicates = cfg.n_replicates;
    result.n_total_paths = cfg.n_paths * cfg.n_replicates;
    result.replicate_estimates.reserve(cfg.n_replicates);

    uint64_t start_n = cfg.skip_first ? 1 : 0;

    for (Size rep = 0; rep < cfg.n_replicates; ++rep) {
        uint64_t scramble_seed = (cfg.base_seed == 0 && rep == 0)
            ? 0
            : cfg.base_seed + static_cast<uint64_t>(rep) * 0x9E3779B97F4A7C15ULL + 1;

        SobolSequence sobol(n_steps, scramble_seed);
        Real sum = 0.0;
        for (uint64_t i = 0; i < cfg.n_paths; ++i) {
            auto uniforms = sobol(start_n + i);
            std::vector<Real> W_path;
            if (cfg.use_brownian_bridge) {
                W_path = bb.generate_path(uniforms);  // W(t_1)..W(t_n), 长度 n_steps
            } else {
                // 标准 incremental 构造
                W_path.resize(n_steps);
                Real W = 0.0;
                for (Size s = 0; s < n_steps; ++s) {
                    Real u = uniforms[s];
                    if (u <= 0.0) u = 1e-15;
                    if (u >= 1.0) u = 1.0 - 1e-15;
                    // 增量构造: dW = sqrt(dt) * Z, 累积得 W(t)
                    W += std::sqrt(dt) * inv_normal_cdf(u);
                    W_path[s] = W;
                }
            }
            // 从 W(t_1)..W(t_n) 构造 S(t_1)..S(t_n) (GBM Exact scheme)
            // S(t_{s+1}) = S(t_s) * exp(drift_dt + σ * (W(t_{s+1}) - W(t_s)))
            std::vector<Real> S_path(n_steps);
            Real S_prev = S0;
            for (Size s = 0; s < n_steps; ++s) {
                Real dW = (s == 0) ? W_path[0] : (W_path[s] - W_path[s - 1]);
                S_prev = S_prev * std::exp(drift_dt + sigma * dW);
                S_path[s] = S_prev;
            }
            sum += payoff(S_path);
        }
        Real est = sum / static_cast<Real>(cfg.n_paths);
        result.replicate_estimates.push_back(est);
    }

    Real mean_est = std::accumulate(result.replicate_estimates.begin(),
                                      result.replicate_estimates.end(), 0.0)
                    / static_cast<Real>(cfg.n_replicates);
    result.price = mean_est * cfg.df;

    if (cfg.n_replicates >= 2) {
        Real var = 0.0;
        for (Real e : result.replicate_estimates) {
            var += (e - mean_est) * (e - mean_est);
        }
        var /= static_cast<Real>(cfg.n_replicates - 1);
        result.std_error = std::sqrt(var / static_cast<Real>(cfg.n_replicates)) * cfg.df;
    }
    result.ci_lower = result.price - 1.96 * result.std_error;
    result.ci_upper = result.price + 1.96 * result.std_error;
    return result;
}

// ============ 多资产欧式期权 QMC ============
// 维度 = n_assets
// Sobol uniforms[n_assets] → W_a(T) = sqrt(T) * inv_normal_cdf(uniforms[a])
// Cholesky 相关化: W_correlated = L * W_independent
// S_a(T) = S0_a * exp((r-q_a-0.5σ_a²)T + σ_a * W_correlated_a)
inline QMCResult price_european_qmc_multi(
        const std::vector<Real>& S0,
        const std::vector<Real>& sigma,
        const std::vector<Real>& q,
        Real r, Real T,
        const std::vector<std::vector<Real>>& correlation,
        const MultiAssetPayoff& payoff,
        const QMCConfig& cfg = {}) {
    const Size n_assets = S0.size();
    if (n_assets == 0) throw std::invalid_argument("QMC: S0 empty");
    if (sigma.size() != n_assets) throw std::invalid_argument("QMC: sigma size mismatch");
    if (q.size() != n_assets && !q.empty()) throw std::invalid_argument("QMC: q size mismatch");
    if (correlation.size() != n_assets) throw std::invalid_argument("QMC: correlation size mismatch");
    if (!payoff) throw std::invalid_argument("QMC: payoff is null");

    // Cholesky 分解
    std::vector<std::vector<Real>> L(n_assets, std::vector<Real>(n_assets, 0.0));
    for (Size i = 0; i < n_assets; ++i) {
        for (Size j = 0; j <= i; ++j) {
            Real s = correlation[i][j];
            for (Size k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
            if (i == j) {
                if (s < -1e-12) throw std::invalid_argument("QMC: correlation not PSD");
                L[i][j] = std::sqrt(std::max(s, 0.0));
            } else {
                L[i][j] = (L[j][j] > 1e-15) ? s / L[j][j] : 0.0;
            }
        }
    }

    std::vector<Real> drift_T(n_assets);
    for (Size a = 0; a < n_assets; ++a) {
        Real qa = q.empty() ? 0.0 : q[a];
        drift_T[a] = (r - qa - 0.5 * sigma[a] * sigma[a]) * T;
    }
    Real sqrt_T = std::sqrt(T);

    QMCResult result;
    result.n_paths = cfg.n_paths;
    result.n_replicates = cfg.n_replicates;
    result.n_total_paths = cfg.n_paths * cfg.n_replicates;
    result.replicate_estimates.reserve(cfg.n_replicates);

    uint64_t start_n = cfg.skip_first ? 1 : 0;

    for (Size rep = 0; rep < cfg.n_replicates; ++rep) {
        uint64_t scramble_seed = (cfg.base_seed == 0 && rep == 0)
            ? 0
            : cfg.base_seed + static_cast<uint64_t>(rep) * 0x9E3779B97F4A7C15ULL + 1;

        SobolSequence sobol(n_assets, scramble_seed);
        Real sum = 0.0;
        for (uint64_t i = 0; i < cfg.n_paths; ++i) {
            auto uniforms = sobol(start_n + i);
            // 1. 生成 n_assets 个独立 W(T)
            std::vector<Real> W_indep(n_assets);
            for (Size a = 0; a < n_assets; ++a) {
                Real u = uniforms[a];
                if (u <= 0.0) u = 1e-15;
                if (u >= 1.0) u = 1.0 - 1e-15;
                W_indep[a] = sqrt_T * inv_normal_cdf(u);
            }
            // 2. Cholesky 相关化
            std::vector<Real> W_corr(n_assets, 0.0);
            for (Size a = 0; a < n_assets; ++a) {
                for (Size b = 0; b <= a; ++b) {
                    W_corr[a] += L[a][b] * W_indep[b];
                }
            }
            // 3. 计算每个资产的 S(T) — W_corr[a] = W_a(T) 已含 sqrt(T) 因子
            std::vector<std::vector<Real>> paths(n_assets, std::vector<Real>(1));
            for (Size a = 0; a < n_assets; ++a) {
                paths[a][0] = S0[a] * std::exp(drift_T[a] + sigma[a] * W_corr[a]);
            }
            sum += payoff(paths);
        }
        Real est = sum / static_cast<Real>(cfg.n_paths);
        result.replicate_estimates.push_back(est);
    }

    Real mean_est = std::accumulate(result.replicate_estimates.begin(),
                                      result.replicate_estimates.end(), 0.0)
                    / static_cast<Real>(cfg.n_replicates);
    result.price = mean_est * cfg.df;

    if (cfg.n_replicates >= 2) {
        Real var = 0.0;
        for (Real e : result.replicate_estimates) {
            var += (e - mean_est) * (e - mean_est);
        }
        var /= static_cast<Real>(cfg.n_replicates - 1);
        result.std_error = std::sqrt(var / static_cast<Real>(cfg.n_replicates)) * cfg.df;
    }
    result.ci_lower = result.price - 1.96 * result.std_error;
    result.ci_upper = result.price + 1.96 * result.std_error;
    return result;
}

// ============ 多资产路径相关 QMC (Brownian Bridge) ============
// 维度 = n_steps * n_assets
// 使用 MultiAssetBrownianBridge 生成相关 Brownian 路径
inline QMCResult price_path_dependent_qmc_multi(
        const std::vector<Real>& S0,
        const std::vector<Real>& sigma,
        const std::vector<Real>& q,
        Real r, Real T, Size n_steps,
        const std::vector<std::vector<Real>>& correlation,
        const MultiAssetPayoff& payoff,
        const QMCConfig& cfg = {}) {
    const Size n_assets = S0.size();
    if (n_assets == 0) throw std::invalid_argument("QMC: S0 empty");
    if (sigma.size() != n_assets) throw std::invalid_argument("QMC: sigma size mismatch");
    if (q.size() != n_assets && !q.empty()) throw std::invalid_argument("QMC: q size mismatch");
    if (correlation.size() != n_assets) throw std::invalid_argument("QMC: correlation size mismatch");
    if (n_steps == 0) throw std::invalid_argument("QMC: n_steps must be positive");
    if (!payoff) throw std::invalid_argument("QMC: payoff is null");

    MultiAssetBrownianBridge mabb(n_steps, T, n_assets, correlation);

    Real dt = T / static_cast<Real>(n_steps);
    std::vector<Real> drift_dt(n_assets);
    for (Size a = 0; a < n_assets; ++a) {
        Real qa = q.empty() ? 0.0 : q[a];
        drift_dt[a] = (r - qa - 0.5 * sigma[a] * sigma[a]) * dt;
    }

    QMCResult result;
    result.n_paths = cfg.n_paths;
    result.n_replicates = cfg.n_replicates;
    result.n_total_paths = cfg.n_paths * cfg.n_replicates;
    result.replicate_estimates.reserve(cfg.n_replicates);

    uint64_t start_n = cfg.skip_first ? 1 : 0;
    Size sobol_dim = n_steps * n_assets;

    for (Size rep = 0; rep < cfg.n_replicates; ++rep) {
        uint64_t scramble_seed = (cfg.base_seed == 0 && rep == 0)
            ? 0
            : cfg.base_seed + static_cast<uint64_t>(rep) * 0x9E3779B97F4A7C15ULL + 1;

        SobolSequence sobol(sobol_dim, scramble_seed);
        Real sum = 0.0;
        for (uint64_t i = 0; i < cfg.n_paths; ++i) {
            auto uniforms = sobol(start_n + i);
            // 1. Brownian Bridge 生成相关 Brownian 路径 W_a(t_1)..W_a(t_n)
            auto W_paths = mabb.generate_paths(uniforms);  // [n_assets][n_steps]
            // 2. GBM Exact scheme: dW = W(t+dt) - W(t) 已含 sqrt(dt) 因子, vol 系数为 σ_a
            std::vector<std::vector<Real>> S_paths(n_assets, std::vector<Real>(n_steps));
            for (Size a = 0; a < n_assets; ++a) {
                Real S_prev = S0[a];
                for (Size s = 0; s < n_steps; ++s) {
                    Real dW = (s == 0) ? W_paths[a][0] : (W_paths[a][s] - W_paths[a][s - 1]);
                    S_prev = S_prev * std::exp(drift_dt[a] + sigma[a] * dW);
                    S_paths[a][s] = S_prev;
                }
            }
            sum += payoff(S_paths);
        }
        Real est = sum / static_cast<Real>(cfg.n_paths);
        result.replicate_estimates.push_back(est);
    }

    Real mean_est = std::accumulate(result.replicate_estimates.begin(),
                                      result.replicate_estimates.end(), 0.0)
                    / static_cast<Real>(cfg.n_replicates);
    result.price = mean_est * cfg.df;

    if (cfg.n_replicates >= 2) {
        Real var = 0.0;
        for (Real e : result.replicate_estimates) {
            var += (e - mean_est) * (e - mean_est);
        }
        var /= static_cast<Real>(cfg.n_replicates - 1);
        result.std_error = std::sqrt(var / static_cast<Real>(cfg.n_replicates)) * cfg.df;
    }
    result.ci_lower = result.price - 1.96 * result.std_error;
    result.ci_upper = result.price + 1.96 * result.std_error;
    return result;
}

}  // namespace v1
}  // namespace cpphub
