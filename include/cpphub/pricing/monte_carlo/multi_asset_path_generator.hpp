#pragma once
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" Ch.3
// SOURCE: Hull (2018) "Options, Futures, and Other Derivatives" Ch.27 (Multi-Asset)
// 模块: 多资产 GBM 路径生成器
//
// 多资产 GBM (Black-Scholes 多维):
//   dS_i(t)/S_i(t) = (r - q_i) dt + σ_i dW_i(t)
//   dW_i dW_j = ρ_ij dt
//
// 通过 Cholesky 分解相关矩阵 R = LL^T, 将独立 Brownian Z 转换为相关 Brownian:
//   dW = L * dZ  (Z 为独立标准正态)
//
// Exact scheme (multiplicative):
//   S_i(t+dt) = S_i(t) * exp((r - q_i - 0.5 σ_i²) dt + σ_i sqrt(dt) (L Z)_i)

#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/pricing/monte_carlo/path_generator.hpp"  // next_normal
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ 多资产 GBM 配置 ============
struct MultiAssetGBMConfig {
    std::vector<Real> S0;            // 初始价格 [n_assets]
    std::vector<Real> sigma;         // 波动率 [n_assets]
    std::vector<Real> q;             // 股息率 [n_assets] (可空, 默认 0)
    Real r = 0.05;                   // 无风险利率 (共享)
    Real T = 1.0;
    Size n_steps = 50;
    std::vector<std::vector<Real>> correlation;  // n_assets × n_assets 相关矩阵

    void validate() const {
        const Size n = S0.size();
        if (n == 0) throw std::invalid_argument("MultiAssetGBM: S0 empty");
        if (sigma.size() != n) throw std::invalid_argument("MultiAssetGBM: sigma size mismatch");
        if (!q.empty() && q.size() != n) throw std::invalid_argument("MultiAssetGBM: q size mismatch");
        if (T <= 0.0) throw std::invalid_argument("MultiAssetGBM: T must be positive");
        if (n_steps == 0) throw std::invalid_argument("MultiAssetGBM: n_steps must be positive");
        if (correlation.size() != n) throw std::invalid_argument("MultiAssetGBM: correlation rows mismatch");
        for (Size i = 0; i < n; ++i) {
            if (S0[i] <= 0.0) throw std::invalid_argument("MultiAssetGBM: S0 must be positive");
            if (sigma[i] < 0.0) throw std::invalid_argument("MultiAssetGBM: sigma must be non-negative");
            if (correlation[i].size() != n) throw std::invalid_argument("MultiAssetGBM: correlation cols mismatch");
        }
        // 对角元 = 1, 对称, 半正定
        for (Size i = 0; i < n; ++i) {
            if (std::abs(correlation[i][i] - 1.0) > 1e-10)
                throw std::invalid_argument("MultiAssetGBM: diagonal must be 1");
            for (Size j = 0; j < i; ++j) {
                if (std::abs(correlation[i][j] - correlation[j][i]) > 1e-10)
                    throw std::invalid_argument("MultiAssetGBM: correlation must be symmetric");
            }
        }
    }

    Size n_assets() const { return S0.size(); }
    Real dividend(Size i) const { return q.empty() ? 0.0 : q[i]; }
};

// ============ 多资产 GBM 路径生成器 ============
// 路径存储: paths[asset_idx][step_idx]
// step 0 = S0, step n_steps = S_T
class MultiAssetGBMPathGenerator {
public:
    explicit MultiAssetGBMPathGenerator(MultiAssetGBMConfig cfg)
        : cfg_(std::move(cfg)) {
        cfg_.validate();
        dt_ = cfg_.T / static_cast<Real>(cfg_.n_steps);
        sqrt_dt_ = std::sqrt(dt_);
        L_ = cholesky_semi_definite(cfg_.correlation);
        // 预计算 drift 和 vol
        drift_.resize(cfg_.n_assets());
        vol_.resize(cfg_.n_assets());
        for (Size i = 0; i < cfg_.n_assets(); ++i) {
            drift_[i] = (cfg_.r - cfg_.dividend(i) - 0.5 * cfg_.sigma[i] * cfg_.sigma[i]) * dt_;
            vol_[i] = cfg_.sigma[i] * sqrt_dt_;
        }
    }

    // 生成一条路径, 返回 [n_assets][n_steps+1] 矩阵
    // sign = +1 普通路径, sign = -1 反变量路径 (Z -> -Z)
    std::vector<std::vector<Real>> generate_path(Philox4x64& rng, Real sign = 1.0) const {
        const Size n = cfg_.n_assets();
        std::vector<std::vector<Real>> paths(n, std::vector<Real>(cfg_.n_steps + 1));
        for (Size i = 0; i < n; ++i) paths[i][0] = cfg_.S0[i];

        for (Size step = 1; step <= cfg_.n_steps; ++step) {
            // 1. 生成 n 个独立 Z
            std::vector<Real> Z(n);
            for (Size i = 0; i < n; ++i) Z[i] = sign * next_normal(rng);

            // 2. 相关化: W = L * Z
            std::vector<Real> W(n, 0.0);
            for (Size i = 0; i < n; ++i) {
                for (Size j = 0; j <= i; ++j) {
                    W[i] += L_[i][j] * Z[j];
                }
            }

            // 3. 更新每个 asset (Exact scheme)
            for (Size i = 0; i < n; ++i) {
                paths[i][step] = paths[i][step - 1] * std::exp(drift_[i] + vol_[i] * W[i]);
            }
        }
        return paths;
    }

    // 从预生成的 Z 矩阵生成路径 (用于真正的 antithetic: 同一 Z 序列, sign=±1)
    // Z_matrix[step-1][asset_idx] = 标准正态 (step=1..n_steps)
    std::vector<std::vector<Real>> generate_path_from_Z(
            const std::vector<std::vector<Real>>& Z_matrix, Real sign = 1.0) const {
        const Size n = cfg_.n_assets();
        std::vector<std::vector<Real>> paths(n, std::vector<Real>(cfg_.n_steps + 1));
        for (Size i = 0; i < n; ++i) paths[i][0] = cfg_.S0[i];

        for (Size step = 1; step <= cfg_.n_steps; ++step) {
            std::vector<Real> Z(n);
            for (Size i = 0; i < n; ++i) Z[i] = sign * Z_matrix[step - 1][i];

            std::vector<Real> W(n, 0.0);
            for (Size i = 0; i < n; ++i) {
                for (Size j = 0; j <= i; ++j) {
                    W[i] += L_[i][j] * Z[j];
                }
            }
            for (Size i = 0; i < n; ++i) {
                paths[i][step] = paths[i][step - 1] * std::exp(drift_[i] + vol_[i] * W[i]);
            }
        }
        return paths;
    }

    // 生成 Z 矩阵 (n_steps × n_assets 独立标准正态)
    std::vector<std::vector<Real>> generate_Z_matrix(Philox4x64& rng) const {
        std::vector<std::vector<Real>> Z(cfg_.n_steps, std::vector<Real>(cfg_.n_assets()));
        for (Size s = 0; s < cfg_.n_steps; ++s) {
            for (Size i = 0; i < cfg_.n_assets(); ++i) {
                Z[s][i] = next_normal(rng);
            }
        }
        return Z;
    }

    // 生成一条单资产路径的辅助函数 (n=1 时仍走多资产逻辑, 测试一致性)
    std::vector<Real> generate_single_path(Philox4x64& rng, Real sign = 1.0) const {
        auto paths = generate_path(rng, sign);
        return paths[0];
    }

    const MultiAssetGBMConfig& config() const { return cfg_; }
    Real dt() const { return dt_; }
    const std::vector<std::vector<Real>>& cholesky_L() const { return L_; }

private:
    MultiAssetGBMConfig cfg_;
    Real dt_;
    Real sqrt_dt_;
    std::vector<Real> drift_;
    std::vector<Real> vol_;
    std::vector<std::vector<Real>> L_;  // Cholesky 下三角

    // 半正定 Cholesky (允许秩亏, 如完美相关 ρ=±1)
    static std::vector<std::vector<Real>> cholesky_semi_definite(
            const std::vector<std::vector<Real>>& A) {
        const Size n = A.size();
        std::vector<std::vector<Real>> L(n, std::vector<Real>(n, 0.0));
        for (Size i = 0; i < n; ++i) {
            for (Size j = 0; j <= i; ++j) {
                Real s = A[i][j];
                for (Size k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
                if (i == j) {
                    if (s < -1e-12)
                        throw std::invalid_argument("MultiAssetGBM: correlation not positive semidefinite");
                    L[i][j] = std::sqrt(std::max(s, 0.0));
                } else {
                    L[i][j] = (L[j][j] > 1e-15) ? s / L[j][j] : 0.0;
                }
            }
        }
        return L;
    }
};

// ============ 便捷工厂 ============
inline MultiAssetGBMConfig make_multi_asset_gbm(
        const std::vector<Real>& S0,
        const std::vector<Real>& sigma,
        Real r, Real T, Size n_steps,
        const std::vector<std::vector<Real>>& correlation,
        const std::vector<Real>& q = {}) {
    MultiAssetGBMConfig cfg;
    cfg.S0 = S0;
    cfg.sigma = sigma;
    cfg.r = r;
    cfg.T = T;
    cfg.n_steps = n_steps;
    cfg.correlation = correlation;
    cfg.q = q;
    return cfg;
}

// 单资产 GBM 配置 (n=1, correlation = {{1.0}})
inline MultiAssetGBMConfig make_single_asset_gbm(
        Real S0, Real sigma, Real r, Real q, Real T, Size n_steps) {
    MultiAssetGBMConfig cfg;
    cfg.S0 = {S0};
    cfg.sigma = {sigma};
    cfg.q = {q};
    cfg.r = r;
    cfg.T = T;
    cfg.n_steps = n_steps;
    cfg.correlation = {{1.0}};
    return cfg;
}

}  // namespace v1
}  // namespace cpphub
