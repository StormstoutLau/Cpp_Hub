#pragma once
// SOURCE: Glasserman (2003) "Monte Carlo Methods in Financial Engineering" Ch.4.2
// SOURCE: Hammersley, Handscomb (1964) "Monte Carlo Methods" — antithetic variates
// 模块: 统一的 Antithetic 对偶变量方差缩减
//
// 数学:
//   对称分布 (如 N(0,1)): 若 Z ~ N(0,1), 则 -Z ~ N(0,1) 同分布
//   Antithetic 估计量: X_av = (X(Z) + X(-Z)) / 2
//   Var(X_av) = (Var(X) + Cov(X(Z), X(-Z))) / 2
//   当 X(Z) 与 X(-Z) 强负相关 (单调 payoff) 时, Cov < 0, 方差显著缩减
//
// 设计:
//   本类封装"生成 Z 一次, 用 ±sign 生成配对路径"的模式,
//   消除 mc_engine / path_generator / mc_var 中各自手写 antithetic 的重复.
//   核心原则: 同一组随机数 Z, 通过 sign 翻转产生配对, 不重新生成 Z.

#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include <vector>
#include <span>
#include <utility>
#include <stdexcept>

namespace cpphub {
inline namespace v1 {

// ============ Antithetic 对偶变量采样器 ============
// 包装 MultiAssetGBMPathGenerator, 提供"Z 共享 + sign 翻转"的配对路径生成.
//
// 使用方式:
//   Antithetic ant(gen);
//   Philox4x64 rng(seed, stream);
//   auto [path_plus, path_minus] = ant.generate_path_pair(rng);
//   Real X1 = payoff(path_plus);
//   Real X2 = payoff(path_minus);
//   Real sample = 0.5 * (X1 + X2);  // 一个独立样本 (antithetic 配对均值)
class Antithetic {
public:
    using PathMatrix = std::vector<std::vector<Real>>;  // [asset][step]
    using ZMatrix    = std::vector<std::vector<Real>>;  // [step][asset]

    explicit Antithetic(const MultiAssetGBMPathGenerator& gen)
        : gen_(gen) {}

    // 生成一对 antithetic 路径 (共享 Z)
    // 返回 (paths(sign=+1), paths(sign=-1)), 每个是 [asset][step] 矩阵
    std::pair<PathMatrix, PathMatrix> generate_path_pair(Philox4x64& rng) const {
        ZMatrix Z = gen_.generate_Z_matrix(rng);
        PathMatrix p_plus  = gen_.generate_path_from_Z(Z,  1.0);
        PathMatrix p_minus = gen_.generate_path_from_Z(Z, -1.0);
        return {std::move(p_plus), std::move(p_minus)};
    }

    // 生成 Z 矩阵 (供外部复用, 例如多 payoff 评估)
    ZMatrix generate_Z(Philox4x64& rng) const {
        return gen_.generate_Z_matrix(rng);
    }

    // 从给定 Z 生成 sign=±1 的路径
    PathMatrix path_from_Z(const ZMatrix& Z, Real sign) const {
        return gen_.generate_path_from_Z(Z, sign);
    }

    // 单次 antithetic 采样: 对一个 payoff 返回配对均值
    // X_av = 0.5 * (payoff(path(+Z)) + payoff(path(-Z)))
    template<typename PayoffFn>
    Real sample(const PayoffFn& payoff, Philox4x64& rng) const {
        auto [p_plus, p_minus] = generate_path_pair(rng);
        Real X1 = payoff(p_plus);
        Real X2 = payoff(p_minus);
        return 0.5 * (X1 + X2);
    }

    // 单资产便捷版: payoff 接收单条路径 vector
    // 要求 gen 为单资产配置
    template<typename SinglePathPayoff>
    Real sample_single_asset(const SinglePathPayoff& payoff, Philox4x64& rng) const {
        auto [p_plus, p_minus] = generate_path_pair(rng);
        // 单资产: paths[0] 是唯一资产的路径
        Real X1 = payoff(p_plus[0]);
        Real X2 = payoff(p_minus[0]);
        return 0.5 * (X1 + X2);
    }

    // 同时对主 payoff 和控制变量 payoff 采样
    // 返回 (X_av, Y_av) = (0.5*(X1+X2), 0.5*(Y1+Y2))
    template<typename PayoffFn, typename CvPayoffFn>
    std::pair<Real, Real> sample_with_cv(const PayoffFn& payoff,
                                          const CvPayoffFn& cv_payoff,
                                          Philox4x64& rng) const {
        auto [p_plus, p_minus] = generate_path_pair(rng);
        Real X1 = payoff(p_plus), X2 = payoff(p_minus);
        Real Y1 = cv_payoff(p_plus), Y2 = cv_payoff(p_minus);
        return {0.5 * (X1 + X2), 0.5 * (Y1 + Y2)};
    }

    const MultiAssetGBMPathGenerator& generator() const { return gen_; }

private:
    const MultiAssetGBMPathGenerator& gen_;
};

// ============ Antithetic 用于 VaR 场景 (风险因子向量) ============
// mc_var.hpp 的 antithetic 是对 dR (风险因子扰动) 取反, 而非路径 Z.
// 本函数提供统一的"对一组正态扰动取反"逻辑.
//
// 给定 dR ~ N(0, Σ) (通过 Cholesky 生成), antithetic 配对为 (dR, -dR).
// 对线性/二次 payoff, 这等价于路径 antithetic.
inline std::pair<std::vector<Real>, std::vector<Real>>
antithetic_risk_factor_pair(const std::vector<Real>& dR) {
    std::vector<Real> neg_dR(dR.size());
    for (Size i = 0; i < dR.size(); ++i) neg_dR[i] = -dR[i];
    return {dR, neg_dR};
}

}  // namespace v1
}  // namespace cpphub
