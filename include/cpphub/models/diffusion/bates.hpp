#pragma once
// SOURCE: Bates (1996) "Jumps and Stochastic Volatility: Exchange Rate Processes
//         Implicit in Deutsche Mark Options"
// SOURCE: Merton (1976) "Option Pricing When Underlying Stock Returns Are Discontinuous"
//
// 模块: Bates 随机过程 (Heston 随机波动率 + Merton 对数正态跳跃)
//
// ==================== 过程离散化 ====================
//
// 风险中性 SDE:
//   dS/S = (r - q - λm) dt + sqrt(v) dW^S + dJ
//   dv   = κ(θ - v) dt + σ sqrt(v) dW^v
//   J = Σ_{k=1}^{N_t} (J_k - 1),  N_t ~ Poisson(λt),  J_k ~ LogNormal(μ_J, σ_J²)
//   m = E[J_k - 1] = exp(μ_J + σ_J²/2) - 1
//
// 离散化 (Full Truncation Euler + 跳跃):
//   每步 dt:
//     1. 方差更新 (Full Truncation, 同 Heston):
//        v_pos = max(v, 0)
//        v_{t+dt} = v + κ(θ - v)dt + σ sqrt(v_pos) Z_v sqrt(dt)
//     2. 跳跃采样:
//        N ~ Poisson(λ dt),  跳跃对数收益和 = Σ log(J_k) ~ N(N*μ_J, N*σ_J²)
//     3. 价格更新 (log 形式, 含跳跃补偿):
//        drift = (r - q - λm - 0.5 v_pos) dt
//        S_{t+dt} = S * exp(drift + sqrt(v_pos) Z_S sqrt(dt) + jump_log_sum)
//
// 注: 跳跃补偿 λm 在 drift 中扣除, 保证 E[S_T] = S_0 exp((r-q)T) (风险中性)

#include <cmath>
#include <complex>
#include <span>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/models/diffusion/process.hpp"
#include "cpphub/models/diffusion/heston.hpp"
#include "cpphub/pricing/analytic/bates_cf.hpp"

namespace cpphub {
inline namespace v1 {

// ============ Bates 过程参数 ============
struct BatesParams {
    Real S0;       // 初始价格 > 0
    // Heston 部分
    Real v0;       // 初始方差 > 0
    Real kappa;    // 均值回归速度 > 0
    Real theta;    // 长期方差 > 0
    Real sigma;    // 波动率波动率 > 0
    Real rho;      // 相关性 [-1, 1]
    // Merton 跳跃部分
    Real lambda;   // 跳跃强度 >= 0
    Real mu_J;     // log(J) 的均值
    Real sigma_J;  // log(J) 的标准差 > 0
    // 市场
    Real r;        // 无风险利率
    Real q;        // 股息率
};

// ============ Bates 过程类 ============
class BatesProcess : public StochasticProcess {
public:
    explicit BatesProcess(BatesParams p,
                           HestonScheme scheme = HestonScheme::FullTruncation)
        : params_(p), scheme_(scheme) {
        validate_bates_params(p);
    }

    Size dimension() const override { return 2; }
    Real spot() const override { return params_.S0; }

    Complex characteristic_function(Complex u, Real tau) const override;
    void generate_path(Real T, Size n_steps,
                       std::span<Real> path, Philox4x64& rng) const override;

    const BatesParams& params() const { return params_; }
    HestonScheme scheme() const { return scheme_; }

    // 单步演化 (供测试逐步对比)
    // 输入: 当前 S, v, dt, 两个相关正态 Z1/Z2, 跳跃对数收益和 jump_log_sum
    // 输出: 新的 S (v 通过引用更新)
    Real evolve(Real S, Real& v, Real dt,
                Real z1, Real z2, Real jump_log_sum) const noexcept;

    // Poisson 采样 (Knuth 算法, 适合小 λdt)
    static int sample_poisson(Real lambda_dt, Philox4x64& rng);

    // 对数正态跳跃采样 J = exp(μ_J + σ_J * Z)
    Real sample_lognormal_jump(Philox4x64& rng) const;

    // 参数验证
    static void validate_bates_params(const BatesParams& p);

private:
    BatesParams params_;
    HestonScheme scheme_;

    void generate_correlated_normals(Philox4x64& rng, Real& z1, Real& z2) const;
};

inline void BatesProcess::validate_bates_params(const BatesParams& p) {
    if (p.S0 <= 0.0) throw std::invalid_argument("Bates: S0 must be positive");
    if (p.v0 <= 0.0) throw std::invalid_argument("Bates: v0 must be positive");
    if (p.kappa <= 0.0) throw std::invalid_argument("Bates: kappa must be positive");
    if (p.theta <= 0.0) throw std::invalid_argument("Bates: theta must be positive");
    if (p.sigma <= 0.0) throw std::invalid_argument("Bates: sigma must be positive");
    if (p.rho < -1.0 || p.rho > 1.0) throw std::invalid_argument("Bates: rho must be in [-1, 1]");
    if (p.lambda < 0.0) throw std::invalid_argument("Bates: lambda must be non-negative");
    if (p.sigma_J <= 0.0) throw std::invalid_argument("Bates: sigma_J must be positive");
}

inline void BatesProcess::generate_correlated_normals(Philox4x64& rng,
                                                        Real& z1, Real& z2) const {
    uint64_t r1 = rng();
    uint64_t r2 = rng();
    double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
    double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
    auto [w1, w2] = box_muller(u1, u2);
    Real rho = params_.rho;
    z1 = w1;
    z2 = rho * w1 + std::sqrt(1.0 - rho * rho) * w2;
}

inline int BatesProcess::sample_poisson(Real lambda_dt, Philox4x64& rng) {
    if (lambda_dt <= 0.0) return 0;
    // Knuth 算法 (适合 λdt < 30, 典型期权定价 dt 很小)
    Real L = std::exp(-lambda_dt);
    Real p = 1.0;
    int k = 0;
    do {
        k++;
        uint64_t r = rng();
        double u = (r >> 11) * (1.0 / 9007199254740992.0);
        p *= u;
    } while (p > L);
    return k - 1;
}

inline Real BatesProcess::sample_lognormal_jump(Philox4x64& rng) const {
    uint64_t r1 = rng();
    uint64_t r2 = rng();
    double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
    double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
    auto [z, _] = box_muller(u1, u2);
    (void)_;
    return std::exp(params_.mu_J + params_.sigma_J * z);
}

inline Real BatesProcess::evolve(Real S, Real& v, Real dt,
                                   Real z1, Real z2, Real jump_log_sum) const noexcept {
    Real v_old = v;
    Real v_pos = (v_old > 0.0) ? v_old : 0.0;
    Real sqrt_v = std::sqrt(v_pos);

    Real kappa = params_.kappa;
    Real theta = params_.theta;
    Real sigma = params_.sigma;
    Real r = params_.r;
    Real q = params_.q;
    Real lambda = params_.lambda;
    Real m = bates_jump_compensation(params_.mu_J, params_.sigma_J);

    // 方差更新 (Full Truncation)
    v = v_old + kappa * (theta - v_old) * dt + sigma * sqrt_v * z2 * std::sqrt(dt);
    if (v < 0.0) v = 0.0;

    // 价格更新 (含跳跃补偿 λm)
    Real drift_adj = r - q - lambda * m;
    Real drift = (drift_adj - 0.5 * v_pos) * dt;
    Real diffusion = sqrt_v * z1 * std::sqrt(dt);
    return S * std::exp(drift + diffusion + jump_log_sum);
}

inline void BatesProcess::generate_path(Real T, Size n_steps,
                                          std::span<Real> path, Philox4x64& rng) const {
    Real dt = T / static_cast<Real>(n_steps);
    path[0] = params_.S0;

    Real S = params_.S0;
    Real v = params_.v0;
    Real lambda_dt = params_.lambda * dt;

    for (Size i = 0; i < n_steps; ++i) {
        Real z1, z2;
        generate_correlated_normals(rng, z1, z2);

        // 跳跃采样
        Real jump_log_sum = 0.0;
        if (params_.lambda > 0.0) {
            int n_jumps = sample_poisson(lambda_dt, rng);
            for (int j = 0; j < n_jumps; ++j) {
                Real J = sample_lognormal_jump(rng);
                jump_log_sum += std::log(J);
            }
        }

        S = evolve(S, v, dt, z1, z2, jump_log_sum);
        path[i + 1] = S;
    }
}

inline Complex BatesProcess::characteristic_function(Complex u, Real tau) const {
    BatesCFParams p{params_.v0, params_.kappa, params_.theta,
                    params_.sigma, params_.rho,
                    params_.lambda, params_.mu_J, params_.sigma_J,
                    params_.r, params_.q};
    return bates_characteristic_function(u, tau, params_.S0, p);
}

}  // namespace v1
}  // namespace cpphub
