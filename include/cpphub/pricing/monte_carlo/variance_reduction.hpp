#pragma once
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" Ch.4
// 模块: 统一方差缩减 Decorator
//
// 组合三种核心方差缩减技术:
//   1. Antithetic (对偶变量): Z 与 -Z 配对, 取均值作为独立样本
//   2. Control Variate (控制变量): X_cv = X - β(Y - μ_Y), β* = Cov(X,Y)/Var(Y)
//   3. (可扩展) Importance Sampling: 通过 Girsanov 测度变换
//
// 设计原则:
//   - 单一入口: price_with_vr() 根据 VRConfig 自动分派
//   - 可组合: Antithetic + CV 可同时启用 (Glasserman 4.3)
//   - 诊断: VRResult 返回 beta, variance_reduction_ratio, effective_n_paths
//   - 不重新生成 Z: Antithetic 共享 Z 矩阵 (与 mc_engine.hpp 一致)
//
// 与现有代码的关系:
//   - 本头文件是 mc_engine.hpp 中 antithetic/CV 逻辑的统一封装
//   - mc_engine.hpp 的 price_path_dependent / price_multi_asset 保留向后兼容
//   - 新代码推荐使用 price_with_vr

#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/pricing/monte_carlo/path_dependent_payoffs.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_payoffs.hpp"
#include "cpphub/pricing/monte_carlo/antithetic.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ 统一 VR 配置 ============
struct VRConfig {
    Size n_paths = 10000;             // 目标独立样本数
                                      // antithetic 下内部生成 n_paths/2 对
    uint64_t seed = 42;
    bool use_antithetic = false;      // 对偶变量
    bool use_control_variate = false; // 控制变量 (需提供 cv_payoff + cv_analytic_price)
    Real df = 1.0;                    // 折现因子 exp(-rT)
};

// ============ 统一 VR 结果 ============
struct VRResult {
    Real price = 0.0;                 // 折现后期权价格
    Real std_error = 0.0;             // 标准误差
    Real ci_lower = 0.0;              // 95% 置信区间下界
    Real ci_upper = 0.0;              // 95% 置信区间上界
    Size n_samples = 0;               // 独立样本数 (antithetic 配对均值数)
    Size n_paths_generated = 0;       // 实际路径生成数 (antithetic = 2*n_samples)
    Real beta_cv = 0.0;               // 控制变量最优 β
    Real variance_reduction = 1.0;    // 方差缩减因子 (vs 无缩减, <1 表示缩减)
};

// ============ 统一 VR 引擎 ============
class VREngine {
public:
    // 统计工具 (与 MCEngine 一致)
    static Real mean(const std::vector<Real>& xs) {
        if (xs.empty()) return 0.0;
        return std::accumulate(xs.begin(), xs.end(), 0.0) / static_cast<Real>(xs.size());
    }

    static Real variance(const std::vector<Real>& xs, Real mu = NaN_) {
        if (xs.size() < 2) return 0.0;
        if (std::isnan(mu)) mu = mean(xs);
        Real s = 0.0;
        for (Real x : xs) s += (x - mu) * (x - mu);
        return s / static_cast<Real>(xs.size() - 1);
    }

    static Real covariance(const std::vector<Real>& xs, const std::vector<Real>& ys,
                            Real mu_x = NaN_, Real mu_y = NaN_) {
        if (xs.size() != ys.size() || xs.size() < 2) return 0.0;
        if (std::isnan(mu_x)) mu_x = mean(xs);
        if (std::isnan(mu_y)) mu_y = mean(ys);
        Real s = 0.0;
        for (Size i = 0; i < xs.size(); ++i) s += (xs[i] - mu_x) * (ys[i] - mu_y);
        return s / static_cast<Real>(xs.size() - 1);
    }

    // ============ 单资产路径相关定价 (统一 VR 入口) ============
    // gen: 单资产 GBM 路径生成器
    // payoff: 路径相关 payoff (接收单条路径 vector)
    // cv_payoff: 控制变量 payoff (可选, 需与 use_control_variate=true 配合)
    // cv_analytic_price: 控制变量的解析价格 (折现前)
    static VRResult price_single_asset(const MultiAssetGBMPathGenerator& gen,
                                        const PathPayoff& payoff,
                                        const VRConfig& cfg,
                                        const PathPayoff& cv_payoff = nullptr,
                                        Real cv_analytic_price = 0.0) {
        if (cfg.n_paths == 0) throw std::invalid_argument("VR: n_paths must be positive");
        if (gen.config().n_assets() != 1)
            throw std::invalid_argument("price_single_asset: generator must be single-asset");

        const bool use_cv = cfg.use_control_variate && cv_payoff;
        const bool use_ant = cfg.use_antithetic;
        // antithetic: 生成 n_paths/2 对, 每对取均值为 1 个独立样本
        const Size n_pairs = use_ant ? std::max<Size>(cfg.n_paths / 2, 1) : cfg.n_paths;

        std::vector<Real> Xs, Ys;
        Xs.reserve(n_pairs);
        if (use_cv) Ys.reserve(n_pairs);

        Antithetic ant(gen);
        for (Size p = 0; p < n_pairs; ++p) {
            Philox4x64 rng(cfg.seed, static_cast<uint64_t>(p));
            if (use_ant) {
                // 共享 Z, sign=±1 生成配对路径, 取均值
                auto [p_plus, p_minus] = ant.generate_path_pair(rng);
                const auto& path1 = p_plus[0];
                const auto& path2 = p_minus[0];
                Real X1 = payoff(path1);
                Real X2 = payoff(path2);
                Xs.push_back(0.5 * (X1 + X2));
                if (use_cv) {
                    Real Y1 = cv_payoff(path1);
                    Real Y2 = cv_payoff(path2);
                    Ys.push_back(0.5 * (Y1 + Y2));
                }
            } else {
                auto path_mat = gen.generate_path(rng, 1.0);
                const auto& path = path_mat[0];
                Xs.push_back(payoff(path));
                if (use_cv) Ys.push_back(cv_payoff(path));
            }
        }

        // 应用控制变量
        Real beta = 0.0;
        Real var_reduction = 1.0;
        if (use_cv) {
            std::vector<Real> Xs_orig = Xs;
            Real mu_X = mean(Xs);
            Real mu_Y = mean(Ys);
            Real var_Y = variance(Ys, mu_Y);
            Real cov_XY = covariance(Xs, Ys, mu_X, mu_Y);
            if (var_Y > 1e-15) {
                beta = cov_XY / var_Y;
                for (Size i = 0; i < Xs.size(); ++i) {
                    Xs[i] -= beta * (Ys[i] - cv_analytic_price);
                }
                Real var_X_orig = variance(Xs_orig, mu_X);
                Real var_X_cv = variance(Xs, mean(Xs));
                if (var_X_orig > 1e-15) {
                    var_reduction = var_X_cv / var_X_orig;
                }
            }
        }

        // 折现
        for (auto& x : Xs) x *= cfg.df;

        VRResult r;
        r.n_samples = Xs.size();
        r.n_paths_generated = use_ant ? 2 * r.n_samples : r.n_samples;
        r.price = mean(Xs);
        Real var = variance(Xs, r.price);
        r.std_error = std::sqrt(var / static_cast<Real>(r.n_samples));
        r.ci_lower = r.price - 1.96 * r.std_error;
        r.ci_upper = r.price + 1.96 * r.std_error;
        r.beta_cv = beta;
        r.variance_reduction = var_reduction;
        return r;
    }

    // ============ 多资产定价 (统一 VR 入口) ============
    static VRResult price_multi_asset(const MultiAssetGBMPathGenerator& gen,
                                       const MultiAssetPayoff& payoff,
                                       const VRConfig& cfg,
                                       const MultiAssetPayoff& cv_payoff = nullptr,
                                       Real cv_analytic_price = 0.0) {
        if (cfg.n_paths == 0) throw std::invalid_argument("VR: n_paths must be positive");
        if (gen.config().n_assets() < 1)
            throw std::invalid_argument("price_multi_asset: generator must have >=1 asset");

        const bool use_cv = cfg.use_control_variate && cv_payoff;
        const bool use_ant = cfg.use_antithetic;
        const Size n_pairs = use_ant ? std::max<Size>(cfg.n_paths / 2, 1) : cfg.n_paths;

        std::vector<Real> Xs, Ys;
        Xs.reserve(n_pairs);
        if (use_cv) Ys.reserve(n_pairs);

        Antithetic ant(gen);
        for (Size p = 0; p < n_pairs; ++p) {
            Philox4x64 rng(cfg.seed, static_cast<uint64_t>(p));
            if (use_ant) {
                auto [paths1, paths2] = ant.generate_path_pair(rng);
                Real X1 = payoff(paths1);
                Real X2 = payoff(paths2);
                Xs.push_back(0.5 * (X1 + X2));
                if (use_cv) {
                    Real Y1 = cv_payoff(paths1);
                    Real Y2 = cv_payoff(paths2);
                    Ys.push_back(0.5 * (Y1 + Y2));
                }
            } else {
                auto paths = gen.generate_path(rng, 1.0);
                Xs.push_back(payoff(paths));
                if (use_cv) Ys.push_back(cv_payoff(paths));
            }
        }

        // 应用控制变量
        Real beta = 0.0;
        Real var_reduction = 1.0;
        if (use_cv) {
            std::vector<Real> Xs_orig = Xs;
            Real mu_X = mean(Xs);
            Real mu_Y = mean(Ys);
            Real var_Y = variance(Ys, mu_Y);
            Real cov_XY = covariance(Xs, Ys, mu_X, mu_Y);
            if (var_Y > 1e-15) {
                beta = cov_XY / var_Y;
                for (Size i = 0; i < Xs.size(); ++i) {
                    Xs[i] -= beta * (Ys[i] - cv_analytic_price);
                }
                Real var_X_orig = variance(Xs_orig, mu_X);
                Real var_X_cv = variance(Xs, mean(Xs));
                if (var_X_orig > 1e-15) {
                    var_reduction = var_X_cv / var_X_orig;
                }
            }
        }

        for (auto& x : Xs) x *= cfg.df;

        VRResult r;
        r.n_samples = Xs.size();
        r.n_paths_generated = use_ant ? 2 * r.n_samples : r.n_samples;
        r.price = mean(Xs);
        Real var = variance(Xs, r.price);
        r.std_error = std::sqrt(var / static_cast<Real>(r.n_samples));
        r.ci_lower = r.price - 1.96 * r.std_error;
        r.ci_upper = r.price + 1.96 * r.std_error;
        r.beta_cv = beta;
        r.variance_reduction = var_reduction;
        return r;
    }

private:
    static constexpr Real NaN_ = std::numeric_limits<Real>::quiet_NaN();
};

}  // namespace v1
}  // namespace cpphub
