#pragma once
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" Ch.4 (Variance Reduction)
// SOURCE: Boyle, Broadie, Glasserman (1997) "Monte Carlo methods for security pricing"
// 模块: 通用 MC 引擎
//
// 功能:
//   1. 欧式期权定价 (单/多资产, 路径相关/无关)
//   2. 反变量方差缩减 (antithetic variates): Z 和 -Z 配对
//   3. 控制变量方差缩减 (control variates): 已知 closed form 的相关 payoff
//   4. 标准误差 + 95% 置信区间
//
// 数学:
//   Antithetic: X_av = (X(Z) + X(-Z))/2, Var(X_av) = (Var(X) + Cov(X(Z),X(-Z)))/2
//   Control Variate: X_cv = X - β(Y - μ_Y), β* = Cov(X,Y)/Var(Y), Var(X_cv) = Var(X) - β² Var(Y)
//   两者组合: 在 antithetic 后应用 CV (Glasserman 4.3)

#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/pricing/monte_carlo/path_dependent_payoffs.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_payoffs.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>

namespace cpphub {
inline namespace v1 {

// ============ MC 配置 ============
struct MCConfig {
    Size n_paths = 10000;           // 路径数 (antithetic 下实际生成 n_paths/2 对)
    uint64_t seed = 42;
    bool use_antithetic = false;    // 反变量
    bool use_control_variate = false; // 控制变量
    Real df = 1.0;                  // 折现因子 exp(-rT)
};

// ============ MC 结果 ============
struct MCResult {
    Real price = 0.0;         // 折现后期权价格
    Real std_error = 0.0;     // 标准误差 (price 的)
    Real ci_lower = 0.0;      // 95% 置信区间下界
    Real ci_upper = 0.0;      // 95% 置信区间上界
    Size n_paths = 0;         // 实际路径数 (antithetic 下为 n_paths)
    Real beta_cv = 0.0;       // 控制变量最优 β (若启用)
    Real variance_reduction = 1.0;  // 方差缩减因子 (vs 无缩减)
};

// ============ MC 引擎 (统计工具) ============
class MCEngine {
public:
    // 计算均值和标准误差
    static Real mean(const std::vector<Real>& xs) {
        if (xs.empty()) return 0.0;
        Real s = std::accumulate(xs.begin(), xs.end(), 0.0);
        return s / static_cast<Real>(xs.size());
    }

    // 样本方差 (无偏估计, 除以 n-1)
    static Real variance(const std::vector<Real>& xs, Real mu = std::numeric_limits<Real>::quiet_NaN()) {
        if (xs.size() < 2) return 0.0;
        if (std::isnan(mu)) mu = mean(xs);
        Real s = 0.0;
        for (Real x : xs) s += (x - mu) * (x - mu);
        return s / static_cast<Real>(xs.size() - 1);
    }

    // 协方差 (无偏)
    static Real covariance(const std::vector<Real>& xs, const std::vector<Real>& ys,
                            Real mu_x = std::numeric_limits<Real>::quiet_NaN(),
                            Real mu_y = std::numeric_limits<Real>::quiet_NaN()) {
        if (xs.size() != ys.size() || xs.size() < 2) return 0.0;
        if (std::isnan(mu_x)) mu_x = mean(xs);
        if (std::isnan(mu_y)) mu_y = mean(ys);
        Real s = 0.0;
        for (Size i = 0; i < xs.size(); ++i) s += (xs[i] - mu_x) * (ys[i] - mu_y);
        return s / static_cast<Real>(xs.size() - 1);
    }

    // 从样本构造结果 (含置信区间)
    static MCResult make_result(const std::vector<Real>& discounted_payoffs,
                                  Real df, Real variance_reduction = 1.0,
                                  Real beta_cv = 0.0) {
        MCResult r;
        r.n_paths = discounted_payoffs.size();
        r.price = mean(discounted_payoffs);
        Real var = variance(discounted_payoffs, r.price);
        r.std_error = std::sqrt(var / static_cast<Real>(r.n_paths));
        // 95% CI: ±1.96 * se
        r.ci_lower = r.price - 1.96 * r.std_error;
        r.ci_upper = r.price + 1.96 * r.std_error;
        r.beta_cv = beta_cv;
        r.variance_reduction = variance_reduction;
        return r;
    }
};

// ============ 单资产 MC 定价 (使用 MultiAssetGBMPathGenerator 单资产) ============
// 接收单资产 path generator, 路径相关 payoff
inline MCResult price_path_dependent(const MultiAssetGBMPathGenerator& gen,
                                       const PathPayoff& payoff,
                                       const MCConfig& cfg,
                                       // 控制变量: cv_payoff 路径相关, cv_analytic_price 已知
                                       const PathPayoff& cv_payoff = nullptr,
                                       Real cv_analytic_price = 0.0) {
    if (cfg.n_paths == 0) throw std::invalid_argument("MC: n_paths must be positive");
    if (gen.config().n_assets() != 1)
        throw std::invalid_argument("price_path_dependent: generator must be single-asset");

    bool use_cv = cfg.use_control_variate && cv_payoff;
    // Antithetic 下生成 n_paths/2 对, 每对取均值作为独立样本
    Size n_pairs = cfg.use_antithetic ? (cfg.n_paths / 2) : cfg.n_paths;
    if (n_pairs == 0) n_pairs = 1;
    Size n_independent = n_pairs;  // 独立样本数

    std::vector<Real> Xs, Ys;
    Xs.reserve(n_independent);
    if (use_cv) Ys.reserve(n_independent);

    for (Size p = 0; p < n_pairs; ++p) {
        Philox4x64 rng(cfg.seed, static_cast<uint64_t>(p));
        if (cfg.use_antithetic) {
            // 真正的 antithetic: 同一 Z 序列, sign=±1 生成配对路径
            auto Z = gen.generate_Z_matrix(rng);
            auto path_mat1 = gen.generate_path_from_Z(Z, 1.0);
            auto path_mat2 = gen.generate_path_from_Z(Z, -1.0);
            const auto& path1 = path_mat1[0];
            const auto& path2 = path_mat2[0];
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
        Real mu_X = MCEngine::mean(Xs);
        Real mu_Y = MCEngine::mean(Ys);
        Real var_Y = MCEngine::variance(Ys, mu_Y);
        Real cov_XY = MCEngine::covariance(Xs, Ys, mu_X, mu_Y);
        if (var_Y > 1e-15) {
            beta = cov_XY / var_Y;
            for (Size i = 0; i < Xs.size(); ++i) {
                Xs[i] -= beta * (Ys[i] - cv_analytic_price);
            }
            Real var_X_orig = MCEngine::variance(Xs_orig, mu_X);
            Real var_X_cv = MCEngine::variance(Xs, MCEngine::mean(Xs));
            if (var_X_orig > 1e-15) {
                var_reduction = var_X_cv / var_X_orig;
            }
        }
    }

    // 折现
    for (auto& x : Xs) x *= cfg.df;
    return MCEngine::make_result(Xs, cfg.df, var_reduction, beta);
}

// ============ 多资产 MC 定价 ============
// 接收多资产 path generator, 多资产 payoff
inline MCResult price_multi_asset(const MultiAssetGBMPathGenerator& gen,
                                    const MultiAssetPayoff& payoff,
                                    const MCConfig& cfg,
                                    const MultiAssetPayoff& cv_payoff = nullptr,
                                    Real cv_analytic_price = 0.0) {
    if (cfg.n_paths == 0) throw std::invalid_argument("MC: n_paths must be positive");
    if (gen.config().n_assets() < 1)
        throw std::invalid_argument("price_multi_asset: generator must have >=1 asset");

    bool use_cv = cfg.use_control_variate && cv_payoff;
    Size n_pairs = cfg.use_antithetic ? (cfg.n_paths / 2) : cfg.n_paths;
    if (n_pairs == 0) n_pairs = 1;

    std::vector<Real> Xs, Ys;
    Xs.reserve(n_pairs);
    if (use_cv) Ys.reserve(n_pairs);

    for (Size p = 0; p < n_pairs; ++p) {
        Philox4x64 rng(cfg.seed, static_cast<uint64_t>(p));
        if (cfg.use_antithetic) {
            auto Z = gen.generate_Z_matrix(rng);
            auto paths1 = gen.generate_path_from_Z(Z, 1.0);
            auto paths2 = gen.generate_path_from_Z(Z, -1.0);
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
        Real mu_X = MCEngine::mean(Xs);
        Real mu_Y = MCEngine::mean(Ys);
        Real var_Y = MCEngine::variance(Ys, mu_Y);
        Real cov_XY = MCEngine::covariance(Xs, Ys, mu_X, mu_Y);
        if (var_Y > 1e-15) {
            beta = cov_XY / var_Y;
            for (Size i = 0; i < Xs.size(); ++i) {
                Xs[i] -= beta * (Ys[i] - cv_analytic_price);
            }
            Real var_X_orig = MCEngine::variance(Xs_orig, mu_X);
            Real var_X_cv = MCEngine::variance(Xs, MCEngine::mean(Xs));
            if (var_X_orig > 1e-15) {
                var_reduction = var_X_cv / var_X_orig;
            }
        }
    }

    for (auto& x : Xs) x *= cfg.df;
    return MCEngine::make_result(Xs, cfg.df, var_reduction, beta);
}

}  // namespace v1
}  // namespace cpphub
