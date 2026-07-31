#pragma once
// SOURCE: Heath, Jarrow, Morton (1992) "Bond Pricing and the Term Structure of Interest Rates"
// SOURCE: Brigo & Mercurio (2006) "Interest Rate Models - Theory and Practice" Ch.5
// 模块: HJM (Heath-Jarrow-Morton) 远期利率框架
//
// HJM 直接对瞬时远期利率 f(t,T) 建模 (T≥t):
//   df(t,T) = μ(t,T) dt + Σ_d σ_d(t,T) dW_d(t)
//
// 无套利条件 (HJM drift theorem, risk-neutral measure Q):
//   μ(t,T) = Σ_d σ_d(t,T) ∫_t^T σ_d(t,s) ds    (单因子: μ = σ(t,T) ∫_t^T σ(t,s) ds)
//
// 零息债价格:
//   P(t,T) = exp(-∫_t^T f(t,s) ds)
//
// 短期利率 r(t) = f(t,t), P(t,T) 也满足 dP/P = r dt + (bond vol) dW
//
// 实现: 离散 tenor grid {T_0=0, T_1, ..., T_N}
//   演化 forward vector f(t, T_i) for T_i > t (frozen tenor 方法, Brace-Gatarek-Musiela)
//   Euler scheme: f(t+dt, T_i) = f(t, T_i) + μ(t,T_i) dt + σ(t,T_i) dW(t)
//
// Markovian HJM: 当波动率结构为 σ(t,T) = σ exp(-κ(T-t)) 时, HJM 退化为 Hull-White
//   即整个 forward curve 可由 (r(t), 时间 t) 有限维 Markov 状态表示
//   验证: f(t,T) = f(0,T) + ∫_0^t α(s,T) ds + σ exp(-κ(T-t)) ∫_0^t exp(-κ(t-s)) dW(s)
//        其中 α(s,T) = σ²/κ [exp(-κ(T-s)) - exp(-2κ(T-s))] (无套利 drift 积分)

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/rng.hpp"
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>

namespace cpphub {
inline namespace v1 {

// ============ HJM 配置 ============
struct HJMConfig {
    // 初始远期曲线: f(0, T_i) 在 tenor grid 上取值 (T_0=0, T_N=T_max)
    std::vector<Real> tenors;            // T_0=0, T_1, ..., T_N (严格递增)
    std::vector<Real> initial_forwards;  // f(0, T_0), f(0, T_1), ..., f(0, T_N) (长度 N+1)

    // 波动率结构: σ_d(t, T) → Real
    // 单因子: 一个 VolFn; 多因子: 多个 VolFn (各因子独立, 通过 dW 相关性耦合)
    // 简化: 时间齐次 (σ 仅依赖 τ=T-t), VolFn 接受 (τ) 返回 σ
    using VolFn = std::function<Real(Real)>;  // σ(τ=T-t)
    std::vector<VolFn> vol_functions;    // 因子数 D = vol_functions.size()

    // 因子间相关性 (D×D 矩阵, Cholesky 分解用于生成相关 Brownian)
    // 单因子时忽略; 多因子需提供正定相关矩阵
    std::vector<std::vector<Real>> correlation;  // 空或 D×D 单位矩阵

    Real dt = 0.01;  // Euler 步长 (与 tenor grid 解耦, 用于演化 forward curve)
    Size n_factors() const { return vol_functions.size(); }
    Size n_tenors() const { return tenors.size(); }

    void validate() const {
        if (tenors.empty()) throw std::invalid_argument("HJM: tenors must not be empty");
        if (tenors.size() != initial_forwards.size()) {
            throw std::invalid_argument("HJM: tenors/initial_forwards size mismatch");
        }
        if (tenors[0] != 0.0) {
            throw std::invalid_argument("HJM: tenors[0] must be 0 (spot)");
        }
        for (Size i = 1; i < tenors.size(); ++i) {
            if (tenors[i] <= tenors[i - 1]) {
                throw std::invalid_argument("HJM: tenors must be strictly increasing");
            }
        }
        if (vol_functions.empty()) {
            throw std::invalid_argument("HJM: at least 1 factor required");
        }
        for (auto& f : vol_functions) {
            if (!f) throw std::invalid_argument("HJM: null vol function");
        }
        Size D = n_factors();
        if (D > 1) {
            if (correlation.size() != D || correlation[0].size() != D) {
                throw std::invalid_argument("HJM: correlation must be D×D for multi-factor");
            }
            // 简单检查对角线
            for (Size i = 0; i < D; ++i) {
                if (std::abs(correlation[i][i] - 1.0) > 1e-10) {
                    throw std::invalid_argument("HJM: correlation diagonal must be 1");
                }
            }
        }
        if (dt <= 0.0) throw std::invalid_argument("HJM: dt must be positive");
    }
};

// ============ HJM 模型 ============
// 通过离散 tenor grid 演化 forward curve, 实现 zero_coupon_bond(t,T) 查询
class HJM {
public:
    explicit HJM(HJMConfig cfg) : cfg_(std::move(cfg)) {
        cfg_.validate();
        // 预计算 Cholesky 分解 (多因子)
        Size D = cfg_.n_factors();
        if (D > 1) {
            chol_L_ = cholesky(cfg_.correlation);
        } else {
            chol_L_.assign(D, std::vector<Real>(D, 0.0));
            chol_L_[0][0] = 1.0;
        }
    }

    // 初始远期利率 f(0, T)
    Real forward_rate(Real T) const {
        return interpolate_forward(0.0, T, cfg_.initial_forwards);
    }

    // 零息债价格 P(0, T) = exp(-∫_0^T f(0,s) ds)
    Real zero_coupon_bond(Real T) const {
        if (T <= 0.0) return 1.0;
        // 梯形积分 f(0,·) over [0, T]
        Real integral = integrate_forward(0.0, T, cfg_.initial_forwards);
        return std::exp(-integral);
    }

    // 无套利 drift μ(0, T) = Σ_d σ_d(T) ∫_0^T σ_d(s) ds  (t=0 时刻)
    // 注: 时间齐次 σ(τ) 下, μ(t,T) = Σ_d σ_d(T-t) ∫_t^T σ_d(s-t) ds
    //                          = Σ_d σ_d(τ) ∫_0^τ σ_d(u) du   (τ=T-t)
    Real drift(Real tau) const {
        Real mu = 0.0;
        for (auto& sigma : cfg_.vol_functions) {
            Real s_tau = sigma(tau);
            // ∫_0^τ σ(u) du (梯形积分, 自适应步长)
            Real integ = integrate_vol(sigma, 0.0, tau);
            mu += s_tau * integ;
        }
        return mu;
    }

    // 短期利率 r(t) = f(t, t)
    // 模拟后, 当前路径在 t=r_t 的 forward curve 上 f(t, t) 即 r(t)
    Real short_rate() const {
        if (!has_simulated_) {
            // 未模拟时返回 r(0) = f(0, 0)
            return cfg_.initial_forwards[0];
        }
        // 模拟后: forward_state_ 存储当前 t 下的 f(t, T_i) for T_i ≥ t
        // 短期利率 = f(t, t) = 在 t 处插值 (通常 t 落在 tenor grid 上)
        return forward_state_[current_step_];
    }

    // 模拟 forward curve 演化: t: 0 → T_horizon
    // 路径存储: paths_[step][tenor_idx] = f(t_step, T_tenor_idx)
    // 仅保留 T_tenor_idx ≥ t_step 的有效 forward (frozen tenor)
    void simulate_path(Real T_horizon, Size n_steps, Philox4x64& rng) {
        if (T_horizon <= 0.0) throw std::invalid_argument("HJM::simulate: T_horizon must be positive");
        if (n_steps == 0) throw std::invalid_argument("HJM::simulate: n_steps must be positive");
        Real dt = T_horizon / static_cast<Real>(n_steps);
        Size D = cfg_.n_factors();
        Size N = cfg_.n_tenors();  // tenor 数 (含 T_0=0)

        // 初始化 forward state: f(0, T_i) for all i
        std::vector<Real> f_state = cfg_.initial_forwards;
        paths_.assign(n_steps + 1, std::vector<Real>(N));
        paths_[0] = f_state;

        Size current_tenor_idx = 0;  // 下一个未到期的 tenor 索引 (T_i ≥ t)

        for (Size step = 1; step <= n_steps; ++step) {
            Real t = static_cast<Real>(step) * dt;
            // 更新 current_tenor_idx: 跳过已到期 tenor (T_i < t)
            while (current_tenor_idx < N && cfg_.tenors[current_tenor_idx] < t - 1e-12) {
                current_tenor_idx++;
            }

            // 生成相关 Brownian 增量 dW_d ~ N(0, dt), 相关性 via Cholesky
            std::vector<Real> dW(D, 0.0);
            std::vector<Real> z(D);
            for (Size d = 0; d < D; ++d) {
                z[d] = next_normal(rng) * std::sqrt(dt);
            }
            for (Size d = 0; d < D; ++d) {
                for (Size k = 0; k <= d; ++k) {
                    dW[d] += chol_L_[d][k] * z[k];
                }
            }

            // 演化所有未到期 forward f(t, T_i) for T_i ≥ t (i.e., i ≥ current_tenor_idx)
            // 实际上 HJM 演化 forward f(t-, T_i) → f(t+, T_i), 所有 T_i > t- 都更新
            // 简化: 对所有 i ≥ current_tenor_idx, τ_i = T_i - (t-dt)
            Real t_prev = static_cast<Real>(step - 1) * dt;
            for (Size i = current_tenor_idx; i < N; ++i) {
                Real T_i = cfg_.tenors[i];
                Real tau = T_i - t_prev;  // 在 t_prev 时刻的 τ = T_i - t_prev
                if (tau <= 0.0) continue;
                // drift: μ(t_prev, T_i) = Σ_d σ_d(τ) ∫_0^τ σ_d(u) du
                Real mu = 0.0;
                for (auto& sigma : cfg_.vol_functions) {
                    Real s = sigma(tau);
                    Real integ = integrate_vol(sigma, 0.0, tau);
                    mu += s * integ;
                }
                // diffusion: Σ_d σ_d(τ) dW_d
                Real diff = 0.0;
                for (Size d = 0; d < D; ++d) {
                    diff += cfg_.vol_functions[d](tau) * dW[d];
                }
                f_state[i] = f_state[i] + mu * dt + diff;
            }

            paths_[step] = f_state;
        }

        // 保存最终状态 (在 t=T_horizon 时刻的 forward curve)
        forward_state_ = f_state;
        current_step_ = n_steps;
        sim_dt_ = dt;
        sim_T_ = T_horizon;
        has_simulated_ = true;
    }

    // 模拟后查询: f(t_step, T) — 通过路径 paths_[step] 在 tenor grid 上插值
    Real forward_rate_at_step(Size step, Real T) const {
        if (!has_simulated_) throw std::runtime_error("HJM: must simulate first");
        if (step >= paths_.size()) throw std::out_of_range("HJM: step out of range");
        Real t = static_cast<Real>(step) * sim_dt_;
        return interpolate_forward(t, T, paths_[step]);
    }

    // 模拟后查询: P(t_step, T) = exp(-∫_{t_step}^T f(t_step, s) ds)
    Real zero_coupon_bond_at_step(Size step, Real T) const {
        if (!has_simulated_) throw std::runtime_error("HJM: must simulate first");
        if (step >= paths_.size()) throw std::out_of_range("HJM: step out of range");
        Real t = static_cast<Real>(step) * sim_dt_;
        if (T <= t) return 1.0;
        Real integral = integrate_forward(t, T, paths_[step]);
        return std::exp(-integral);
    }

    // 模拟后的当前 forward curve (最后一步)
    const std::vector<Real>& current_forward_state() const {
        if (!has_simulated_) throw std::runtime_error("HJM: must simulate first");
        return forward_state_;
    }

    Real simulated_dt() const { return sim_dt_; }
    Real simulated_T() const { return sim_T_; }
    Size n_steps() const { return paths_.size() - 1; }
    const HJMConfig& config() const { return cfg_; }
    bool has_simulated() const { return has_simulated_; }

    // 模拟后当前时刻 t 的短期利率 r(t) = f(t, t)
    // 注: 若 t 不在 tenor grid 上, 用线性插值
    Real short_rate_at(Size step) const {
        if (!has_simulated_) throw std::runtime_error("HJM: must simulate first");
        Real t = static_cast<Real>(step) * sim_dt_;
        return interpolate_forward(t, t, paths_[step]);
    }

private:
    HJMConfig cfg_;
    std::vector<std::vector<Real>> chol_L_;  // Cholesky 下三角
    std::vector<std::vector<Real>> paths_;   // paths_[step][i] = f(t_step, T_i)
    std::vector<Real> forward_state_;        // 最后一步的 forward state
    Real sim_dt_ = 0.0;
    Real sim_T_ = 0.0;
    Size current_step_ = 0;
    bool has_simulated_ = false;

    // Box-Muller 生成 N(0,1)
    static Real next_normal(Philox4x64& rng) {
        Real u1 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
        if (u1 < 1e-300) u1 = 1e-300;
        Real u2 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
        return std::sqrt(-2.0 * std::log(u1))
             * std::cos(6.283185307179586476925286766559 * u2);
    }

    // 在 tenor grid 上插值 f(t, T)
    // state[i] 存储 f(t, T_i) (tenors_[i] 是 T_i)
    // t 是当前模拟时间; 当 T < t 时无意义, 返回 r(t) 近似
    Real interpolate_forward(Real t, Real T, const std::vector<Real>& state) const {
        // 找到 T 在 tenors_ 中的位置
        const auto& tenors = cfg_.tenors;
        if (T <= t) {
            // 短端: 返回 f(t, t) 即 r(t) (在 t 处插值)
            // 找 t 在 tenors 中的位置
            if (t <= tenors.front()) return state[0];
            if (t >= tenors.back()) return state.back();
            // 线性插值
            auto it = std::upper_bound(tenors.begin(), tenors.end(), t);
            Size idx = static_cast<Size>(it - tenors.begin());
            Real T1 = tenors[idx - 1], T2 = tenors[idx];
            Real f1 = state[idx - 1], f2 = state[idx];
            return f1 + (f2 - f1) * (t - T1) / (T2 - T1);
        }
        // T > t: 在 tenor grid 上插值
        if (T <= tenors.front()) return state[0];
        if (T >= tenors.back()) return state.back();
        auto it = std::upper_bound(tenors.begin(), tenors.end(), T);
        Size idx = static_cast<Size>(it - tenors.begin());
        Real T1 = tenors[idx - 1], T2 = tenors[idx];
        Real f1 = state[idx - 1], f2 = state[idx];
        return f1 + (f2 - f1) * (T - T1) / (T2 - T1);
    }

    // ∫_t^T f(t, s) ds (state 存 f(t, T_i), 在 [t, T] 上梯形积分)
    Real integrate_forward(Real t, Real T, const std::vector<Real>& state) const {
        if (T <= t) return 0.0;
        const auto& tenors = cfg_.tenors;
        // 找到所有落在 [t, T] 内的 tenor 节点
        // 端点 t 和 T 处用插值, 内部节点用 state 直接取值
        // 简化: 用密集子区间梯形积分 (使用 tenor grid 节点 + 端点)
        std::vector<std::pair<Real, Real>> pts;  // (s, f(t, s))
        pts.emplace_back(t, interpolate_forward(t, t, state));
        for (Size i = 0; i < tenors.size(); ++i) {
            if (tenors[i] > t && tenors[i] < T) {
                pts.emplace_back(tenors[i], state[i]);
            }
        }
        pts.emplace_back(T, interpolate_forward(t, T, state));
        // 梯形积分
        Real integral = 0.0;
        for (Size i = 1; i < pts.size(); ++i) {
            integral += 0.5 * (pts[i].second + pts[i - 1].second) * (pts[i].first - pts[i - 1].first);
        }
        return integral;
    }

    // ∫_0^tau σ(u) du (自适应梯形, 因 σ 可能任意函数)
    static Real integrate_vol(const HJMConfig::VolFn& sigma, Real a, Real b) {
        if (b <= a) return 0.0;
        // 简单梯形: 用 64 个子区间 (对常见 σ(τ)=σ·exp(-κτ) 足够精确)
        Size n = 64;
        Real h = (b - a) / static_cast<Real>(n);
        Real sum = 0.5 * (sigma(a) + sigma(b));
        for (Size i = 1; i < n; ++i) {
            sum += sigma(a + static_cast<Real>(i) * h);
        }
        return sum * h;
    }

    // Cholesky 分解 (对称半正定矩阵, 允许秩亏)
    static std::vector<std::vector<Real>> cholesky(const std::vector<std::vector<Real>>& A) {
        Size n = A.size();
        std::vector<std::vector<Real>> L(n, std::vector<Real>(n, 0.0));
        for (Size i = 0; i < n; ++i) {
            for (Size j = 0; j <= i; ++j) {
                Real s = A[i][j];
                for (Size k = 0; k < j; ++k) {
                    s -= L[i][k] * L[j][k];
                }
                if (i == j) {
                    if (s < -1e-12) {
                        throw std::invalid_argument("HJM: correlation not positive semidefinite");
                    }
                    L[i][j] = std::sqrt(std::max(s, 0.0));
                } else {
                    L[i][j] = (L[j][j] > 1e-15) ? s / L[j][j] : 0.0;
                }
            }
        }
        return L;
    }
};

// ============ 便捷工厂: 平坦远期曲线 ============
inline HJMConfig make_flat_hjm(Real r_flat, Real T_max, Real dt_tenor,
                                  HJMConfig::VolFn vol_fn,
                                  Size n_factors = 1) {
    HJMConfig cfg;
    for (Real T = 0.0; T <= T_max + 1e-10; T += dt_tenor) {
        cfg.tenors.push_back(T);
        cfg.initial_forwards.push_back(r_flat);
    }
    cfg.vol_functions.assign(n_factors, vol_fn);
    return cfg;
}

// ============ 便捷工厂: Hull-White 等价 HJM ============
// σ(τ) = σ_hw * exp(-κ * τ), 此时 HJM 退化为 HW
inline HJMConfig make_hw_equivalent_hjm(Real r0, Real kappa, Real sigma_hw,
                                           Real T_max, Real dt_tenor) {
    return make_flat_hjm(r0, T_max, dt_tenor,
                          [sigma_hw, kappa](Real tau) {
                              return sigma_hw * std::exp(-kappa * tau);
                          }, 1);
}

}  // namespace v1
}  // namespace cpphub
