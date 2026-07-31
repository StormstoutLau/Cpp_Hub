#pragma once
// SOURCE: Longstaff & Schwartz (2001) "Valuing American Options by Simulation"
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" Ch.3-4
// 模块: GBM 路径生成器，支持精确解 (exact) 和 Euler 方案，供 LSMC 使用
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/rng.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>

namespace cpphub {
inline namespace v1 {

enum class PathScheme {
    Exact,   // GBM 闭式解: S_{t+dt} = S_t * exp((r-q-0.5*sigma^2)*dt + sigma*sqrt(dt)*Z)
    Euler,   // Euler-Maruyama: S_{t+dt} = S_t * (1 + (r-q)*dt + sigma*sqrt(dt)*Z)
    Milstein // Milstein: Euler + 0.5*sigma^2*Z*(Z-1)*dt
};

struct GBMConfig {
    Real S0 = 100.0;
    Real r = 0.05;
    Real q = 0.0;
    Real sigma = 0.20;
    Real T = 1.0;
    Size n_steps = 50;  // 行使点数 (不含 t=0)
};

// 从 Philox4x64 生成标准正态随机数 (Box-Muller)
inline Real next_normal(Philox4x64& rng) {
    // 使用两个 uniform 生成一个 normal (Box-Muller)
    constexpr Real two_pi = 6.283185307179586476925286766559;
    Real u1, u2;
    // 均匀分布 (0,1]，避免 log(0)
    u1 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
    if (u1 < 1e-300) u1 = 1e-300;
    u2 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
    Real r = std::sqrt(-2.0 * std::log(u1));
    return r * std::cos(two_pi * u2);
}

// GBM 路径生成器
// 路径存储: paths[path_idx * (n_steps+1) + step_idx] = S_{step_idx}
// step 0 = S0, step n_steps = S_T
class GBMPathGenerator {
public:
    explicit GBMPathGenerator(GBMConfig cfg, PathScheme scheme = PathScheme::Exact)
        : cfg_(cfg), scheme_(scheme) {
        if (cfg_.T <= 0.0) throw std::invalid_argument("GBMPathGenerator: T must be positive");
        if (cfg_.n_steps == 0) throw std::invalid_argument("GBMPathGenerator: n_steps must be positive");
        if (cfg_.sigma < 0.0) throw std::invalid_argument("GBMPathGenerator: sigma must be non-negative");
        if (cfg_.S0 <= 0.0) throw std::invalid_argument("GBMPathGenerator: S0 must be positive");
        dt_ = cfg_.T / static_cast<Real>(cfg_.n_steps);
        sqrt_dt_ = std::sqrt(dt_);
    }

    // 生成单条路径，返回 (n_steps+1) 个价格点
    std::vector<Real> generate_path(Philox4x64& rng) const {
        std::vector<Real> path(cfg_.n_steps + 1);
        path[0] = cfg_.S0;
        Real drift, vol;
        switch (scheme_) {
            case PathScheme::Exact:
                drift = (cfg_.r - cfg_.q - 0.5 * cfg_.sigma * cfg_.sigma) * dt_;
                vol = cfg_.sigma * sqrt_dt_;
                for (Size i = 1; i <= cfg_.n_steps; ++i) {
                    Real Z = next_normal(rng);
                    path[i] = path[i - 1] * std::exp(drift + vol * Z);
                }
                break;
            case PathScheme::Euler:
                drift = (cfg_.r - cfg_.q) * dt_;
                vol = cfg_.sigma * sqrt_dt_;
                for (Size i = 1; i <= cfg_.n_steps; ++i) {
                    Real Z = next_normal(rng);
                    path[i] = path[i - 1] * (1.0 + drift + vol * Z);
                    // 防止 Euler 出现负价格
                    if (path[i] <= 0.0) path[i] = 1e-10;
                }
                break;
            case PathScheme::Milstein:
                drift = (cfg_.r - cfg_.q) * dt_;
                vol = cfg_.sigma * sqrt_dt_;
                for (Size i = 1; i <= cfg_.n_steps; ++i) {
                    Real Z = next_normal(rng);
                    path[i] = path[i - 1] * (1.0 + drift + vol * Z
                              + 0.5 * cfg_.sigma * cfg_.sigma * dt_ * (Z * Z - 1.0));
                    if (path[i] <= 0.0) path[i] = 1e-10;
                }
                break;
        }
        return path;
    }

    // 批量生成 n_paths 条路径，返回 (n_paths x (n_steps+1)) 平铺数组
    // 路径 i 的第 j 个点: paths[i * (n_steps+1) + j]
    std::vector<Real> generate_paths(Size n_paths, uint64_t seed, bool antithetic = false) const {
        Size path_len = cfg_.n_steps + 1;
        Size total_paths = antithetic ? n_paths : n_paths;
        // antithetic: 生成 n_paths/2 原始 + n_paths/2 反变量
        Size n_generate = antithetic ? (n_paths + 1) / 2 : n_paths;
        std::vector<Real> paths(total_paths * path_len);
        // 使用确定性分块: 每条路径独立 stream = path_idx
        for (Size p = 0; p < n_generate; ++p) {
            Philox4x64 rng(seed, p);  // key=seed, stream=path_idx
            std::vector<Real> path = generate_path(rng);
            for (Size j = 0; j < path_len; ++j) {
                paths[p * path_len + j] = path[j];
            }
            if (antithetic && p + n_generate < n_paths) {
                // 反变量路径: 翻转所有 Z 的符号
                // 重新生成，但用 -Z (等效于用相同的 uniforms 但取反)
                Philox4x64 rng2(seed, p);
                std::vector<Real> anti_path(cfg_.n_steps + 1);
                anti_path[0] = cfg_.S0;
                Real drift, vol;
                switch (scheme_) {
                    case PathScheme::Exact:
                        drift = (cfg_.r - cfg_.q - 0.5 * cfg_.sigma * cfg_.sigma) * dt_;
                        vol = cfg_.sigma * sqrt_dt_;
                        for (Size i = 1; i <= cfg_.n_steps; ++i) {
                            Real Z = next_normal(rng2);
                            anti_path[i] = anti_path[i - 1] * std::exp(drift - vol * Z);
                        }
                        break;
                    case PathScheme::Euler:
                        drift = (cfg_.r - cfg_.q) * dt_;
                        vol = cfg_.sigma * sqrt_dt_;
                        for (Size i = 1; i <= cfg_.n_steps; ++i) {
                            Real Z = next_normal(rng2);
                            anti_path[i] = anti_path[i - 1] * (1.0 + drift - vol * Z);
                            if (anti_path[i] <= 0.0) anti_path[i] = 1e-10;
                        }
                        break;
                    case PathScheme::Milstein:
                        drift = (cfg_.r - cfg_.q) * dt_;
                        vol = cfg_.sigma * sqrt_dt_;
                        for (Size i = 1; i <= cfg_.n_steps; ++i) {
                            Real Z = next_normal(rng2);
                            anti_path[i] = anti_path[i - 1] * (1.0 + drift - vol * Z
                                      + 0.5 * cfg_.sigma * cfg_.sigma * dt_ * (Z * Z - 1.0));
                            if (anti_path[i] <= 0.0) anti_path[i] = 1e-10;
                        }
                        break;
                }
                for (Size j = 0; j < path_len; ++j) {
                    paths[(p + n_generate) * path_len + j] = anti_path[j];
                }
            }
        }
        return paths;
    }

    const GBMConfig& config() const { return cfg_; }
    Real dt() const { return dt_; }
    Size path_length() const { return cfg_.n_steps + 1; }
    PathScheme scheme() const { return scheme_; }

private:
    GBMConfig cfg_;
    PathScheme scheme_;
    Real dt_ = 0.0;
    Real sqrt_dt_ = 0.0;
};

}  // namespace v1
}  // namespace cpphub
