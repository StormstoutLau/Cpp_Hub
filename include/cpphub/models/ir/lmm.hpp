#pragma once
// SOURCE: Brace, Gatarek, Musiela (1997) "The Market Model of Interest Rate Dynamics"
// SOURCE: Jamshidian (1997) "LIBOR and Swap Market Models and Measures"
// SOURCE: Brigo & Mercurio (2006) "Interest Rate Models - Theory and Practice" Ch.6
// 模块: LIBOR 市场模型 (LMM / BGM)
//
// LMM 直接对离散远期 LIBOR L(t; T_i, T_{i+1}) 建模 (而非瞬时远期 f(t,T)):
//   dL_i(t) = μ_i(t) L_i(t) dt + σ_i(t) L_i(t) dW_i(t)
//   其中 L_i(t) = L(t; T_i, T_{i+1}), τ_i = T_{i+1} - T_i
//
// 测度选择:
//   (1) Spot LIBOR measure (numeraire = rolling bond): drift 紧凑, 推荐使用
//       μ_i(t) = Σ_{j: T_j ≤ T_i} ρ_ij σ_j(t) L_j(t) τ_j / (1 + L_j(t) τ_j)
//       (求和下标 j 从当前 reset index 到 i, 含 i)
//   (2) Terminal measure (numeraire = P(t, T_N)): μ_i(t) = -Σ_{j: T_j > T_i} ρ_ij σ_j(t) L_j(t) τ_j / (1 + L_j(t) τ_j)
//
// 每个 L_i 在其自身 forward measure Q^{T_{i+1}} 下是 driftless martingale (对数正态 → Black 76)
//   → caplet_i 定价 Black 76 严格成立: caplet = τ_i N P(0,T_{i+1}) (F_i N(d1) - K N(d2))
//
// Swaption: Rebonato 近似公式 (swap rate 近似对数正态)
//   σ_swap² T_ex ≈ Σ_i Σ_j w_i w_j ρ_ij σ_i σ_j  (vol 假设常数)
//   swaption ≈ Black 76 on swap rate with σ_swap
//
// 实现: Euler on log L (保证 L > 0)
//   d log L_i = (μ_i - 0.5 σ_i²) dt + σ_i dW_i (Itô)

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

// ============ LMM 配置 ============
struct LMMConfig {
    // Tenor 结构: T_0, T_1, ..., T_N (严格递增, T_0 通常为 spot 0)
    // 共 N 个 LIBOR: L_0, L_1, ..., L_{N-1}, 其中 L_i 对应 [T_i, T_{i+1}]
    std::vector<Real> tenors;            // T_0, ..., T_N (长度 N+1)
    std::vector<Real> initial_libors;    // L_0(0), ..., L_{N-1}(0) (长度 N)

    // 波动率: σ_i(t) (常数或时间依赖); 简化为 σ_i (常数, 或即时 t→σ_i)
    // 这里用即时 t 的函数: σ_i(t) → Real, 由用户决定参数化 (piecewise constant 等)
    // 简化接口: 假设 σ_i 为常数 (用 σ_i 即可); 若需时变可扩展
    std::vector<Real> volatilities;      // σ_0, ..., σ_{N-1}

    // 相关性矩阵 ρ_ij (N×N 对称正定, 对角为 1)
    std::vector<std::vector<Real>> correlation;

    // 选用的测度: "spot" (rolling bond) 或 "terminal" (P(t, T_N))
    std::string measure = "spot";

    Size n_libors() const { return initial_libors.size(); }
    Size n_tenors() const { return tenors.size(); }

    void validate() const {
        if (tenors.size() != initial_libors.size() + 1) {
            throw std::invalid_argument("LMM: tenors size must be N+1 = n_libors + 1");
        }
        if (initial_libors.size() != volatilities.size()) {
            throw std::invalid_argument("LMM: initial_libors/volatilities size mismatch");
        }
        if (tenors.empty()) throw std::invalid_argument("LMM: tenors empty");
        for (Size i = 1; i < tenors.size(); ++i) {
            if (tenors[i] <= tenors[i - 1]) {
                throw std::invalid_argument("LMM: tenors must be strictly increasing");
            }
        }
        Size N = n_libors();
        if (N == 0) throw std::invalid_argument("LMM: at least 1 LIBOR required");
        if (correlation.size() != N) {
            throw std::invalid_argument("LMM: correlation must be N×N");
        }
        for (Size i = 0; i < N; ++i) {
            if (correlation[i].size() != N) {
                throw std::invalid_argument("LMM: correlation must be N×N");
            }
            if (std::abs(correlation[i][i] - 1.0) > 1e-10) {
                throw std::invalid_argument("LMM: correlation diagonal must be 1");
            }
            for (Size j = 0; j < N; ++j) {
                if (std::abs(correlation[i][j] - correlation[j][i]) > 1e-10) {
                    throw std::invalid_argument("LMM: correlation must be symmetric");
                }
                if (std::abs(correlation[i][j]) > 1.0 + 1e-10) {
                    throw std::invalid_argument("LMM: |correlation| must be ≤ 1");
                }
            }
        }
        for (Real v : volatilities) {
            if (v < 0.0) throw std::invalid_argument("LMM: vol must be non-negative");
        }
        for (Real L : initial_libors) {
            if (L <= 0.0) throw std::invalid_argument("LMM: initial LIBOR must be positive");
        }
        if (measure != "spot" && measure != "terminal") {
            throw std::invalid_argument("LMM: measure must be 'spot' or 'terminal'");
        }
    }

    // 便利: 第 i 个 LIBOR 的年化分数 τ_i = T_{i+1} - T_i
    Real tau(Size i) const {
        if (i + 1 >= tenors.size()) {
            throw std::out_of_range("LMM::tau: index out of range");
        }
        return tenors[i + 1] - tenors[i];
    }
};

// ============ LMM 模型 ============
class LMM {
public:
    explicit LMM(LMMConfig cfg) : cfg_(std::move(cfg)) {
        cfg_.validate();
        chol_L_ = cholesky(cfg_.correlation);
        // 预计算 t=0 时的零息债 P(0, T_i) = Π_{k<i} 1/(1 + L_k(0) τ_k)
        Size N = cfg_.n_libors();
        Size M = cfg_.n_tenors();  // N+1
        init_bonds_.assign(M, 1.0);
        init_bonds_[0] = 1.0;  // P(0, T_0) = 1
        for (Size i = 0; i < N; ++i) {
            Real L0 = cfg_.initial_libors[i];
            Real tau_i = cfg_.tau(i);
            init_bonds_[i + 1] = init_bonds_[i] / (1.0 + L0 * tau_i);
        }
    }

    // t=0 时的零息债价格 P(0, T_i)
    Real zero_coupon_bond(Real T) const {
        // 找到 T 对应的 tenor 索引; 超出网格则线性外推
        const auto& tenors = cfg_.tenors;
        if (T <= tenors.front()) {
            // 短端外推 (用第 1 个 LIBOR)
            return std::exp(-cfg_.initial_libors[0] * T);
        }
        if (T >= tenors.back()) {
            // 长端外推 (用最后一个 LIBOR)
            Real tau_last = cfg_.tau(cfg_.n_libors() - 1);
            Real r_last = cfg_.initial_libors.back();
            Real T_back = tenors.back();
            Real P_back = init_bonds_.back();
            return P_back * std::exp(-r_last * (T - T_back));
        }
        // 在网格内: 线性插值 bond 价格
        auto it = std::upper_bound(tenors.begin(), tenors.end(), T);
        Size idx = static_cast<Size>(it - tenors.begin());
        Real T1 = tenors[idx - 1], T2 = tenors[idx];
        Real P1 = init_bonds_[idx - 1], P2 = init_bonds_[idx];
        // log-linear 插值 (保证 P 单调)
        Real y1 = -std::log(P1) / T1;
        Real y2 = -std::log(P2) / T2;
        Real y = y1 + (y2 - y1) * (T - T1) / (T2 - T1);
        return std::exp(-y * T);
    }

    // t=0 时的 LIBOR L(0; T_i, T_{i+1})
    Real initial_libor(Size i) const {
        if (i >= cfg_.n_libors()) throw std::out_of_range("LMM::initial_libor");
        return cfg_.initial_libors[i];
    }

    // 初始 LIBOR L(0; T_a, T_b) for any T_a < T_b on the grid
    Real initial_libor(Real T_a, Real T_b) const {
        if (T_b <= T_a) throw std::invalid_argument("LMM::libor: T_b must > T_a");
        Real P_a = zero_coupon_bond(T_a);
        Real P_b = zero_coupon_bond(T_b);
        Real tau = T_b - T_a;
        return (P_a / P_b - 1.0) / tau;
    }

    // ============ Caplet 定价 (Black 76, 严格解析) ============
    // 第 i 个 caplet: payoff at T_{i+1} = τ_i N max(L_i(T_i) - K, 0)
    // L_i 在 Q^{T_{i+1}} 下对数正态且 driftless, 故 caplet = Black 76
    //   caplet_i = τ_i N P(0, T_{i+1}) (F_i N(d1) - K N(d2))
    //   F_i = L_i(0), d1 = (log(F/K) + 0.5 σ_i² T_i) / (σ_i √T_i), d2 = d1 - σ_i √T_i
    Real caplet_price(Size i, Real K, Real N = 1.0) const {
        if (i >= cfg_.n_libors()) throw std::out_of_range("LMM::caplet_price");
        Real F = cfg_.initial_libors[i];
        Real sigma = cfg_.volatilities[i];
        Real T_i = cfg_.tenors[i];
        Real tau_i = cfg_.tau(i);
        Real P = init_bonds_[i + 1];  // P(0, T_{i+1})
        if (sigma < 1e-15 || T_i < 1e-15) {
            Real intrinsic = std::max(F - K, 0.0);
            return tau_i * N * P * intrinsic;
        }
        Real sqrt_T = std::sqrt(T_i);
        Real d1 = (std::log(F / K) + 0.5 * sigma * sigma * T_i) / (sigma * sqrt_T);
        Real d2 = d1 - sigma * sqrt_T;
        return tau_i * N * P * (F * normal_cdf(d1) - K * normal_cdf(d2));
    }

    // Floorlet (Black 76): payoff = τ_i N max(K - L_i, 0)
    Real floorlet_price(Size i, Real K, Real N = 1.0) const {
        if (i >= cfg_.n_libors()) throw std::out_of_range("LMM::floorlet_price");
        Real F = cfg_.initial_libors[i];
        Real sigma = cfg_.volatilities[i];
        Real T_i = cfg_.tenors[i];
        Real tau_i = cfg_.tau(i);
        Real P = init_bonds_[i + 1];
        if (sigma < 1e-15 || T_i < 1e-15) {
            Real intrinsic = std::max(K - F, 0.0);
            return tau_i * N * P * intrinsic;
        }
        Real sqrt_T = std::sqrt(T_i);
        Real d1 = (std::log(F / K) + 0.5 * sigma * sigma * T_i) / (sigma * sqrt_T);
        Real d2 = d1 - sigma * sqrt_T;
        return tau_i * N * P * (K * normal_cdf(-d2) - F * normal_cdf(-d1));
    }

    // ============ Cap/Floor 价格 (Caplet/Floorlet 之和) ============
    // cap covering L_a, ..., L_b (reset at T_a, ..., T_b; pay at T_{a+1}, ..., T_{b+1})
    Real cap_price(Real K, Real N = 1.0, Size start_idx = 0) const {
        Real pv = 0.0;
        for (Size i = start_idx; i < cfg_.n_libors(); ++i) {
            pv += caplet_price(i, K, N);
        }
        return pv;
    }

    Real floor_price(Real K, Real N = 1.0, Size start_idx = 0) const {
        Real pv = 0.0;
        for (Size i = start_idx; i < cfg_.n_libors(); ++i) {
            pv += floorlet_price(i, K, N);
        }
        return pv;
    }

    // ============ Cap-Floor Parity (无套利验证) ============
    // Cap - Floor = Payer Forward Swap (no optionality)
    // Cap(K) - Floor(K) = N * [P(0,T_start) - P(0,T_end) - K * Σ τ_i P(0,T_{i+1})]
    Real cap_floor_parity_payer_swap(Real K, Real N = 1.0, Size start_idx = 0) const {
        Size end_idx = cfg_.n_libors() - 1;  // 最后一个 LIBOR 索引
        Real P_start = init_bonds_[start_idx];
        Real P_end = init_bonds_[end_idx + 1];
        Real annuity = 0.0;
        for (Size i = start_idx; i <= end_idx; ++i) {
            annuity += cfg_.tau(i) * init_bonds_[i + 1];
        }
        return N * (P_start - P_end - K * annuity);
    }

    // ============ Swaption 近似定价 (Rebonato 公式) ============
    // Payer swaption on swap [T_a, T_b], expiry T_a, strike K
    // Rebonato: swap rate 近似对数正态, σ_swap² T_a ≈ Σ_{i,j=a..b-1} w_i w_j ρ_ij σ_i σ_j
    //   w_i = τ_i P(0, T_{i+1}) / A, A = Σ τ_i P(0, T_{i+1})  (annuity weights)
    //   swap rate S(0) = (P(0,T_a) - P(0,T_b)) / A
    //   swaption = A * (S N(d1) - K N(d2)), d1 = (log(S/K) + 0.5 σ_swap² T_a)/(σ_swap √T_a)
    Real swaption_price_rebonato(Real T_ex, Real K, bool is_payer,
                                    Size start_idx, Size end_idx, Real N = 1.0) const {
        // T_ex 应等于 cfg_.tenors[start_idx]
        if (start_idx >= end_idx) throw std::invalid_argument("LMM::swaption: start_idx < end_idx required");
        if (end_idx >= cfg_.n_libors()) throw std::out_of_range("LMM::swaption: end_idx");

        // Annuity A and weights w_i
        Real A = 0.0;
        for (Size i = start_idx; i <= end_idx; ++i) {
            A += cfg_.tau(i) * init_bonds_[i + 1];
        }
        if (A <= 0.0) throw std::runtime_error("LMM::swaption: annuity non-positive");

        // Swap rate S(0) = (P(0,T_a) - P(0,T_{b+1})) / A
        Real P_a = init_bonds_[start_idx];
        Real P_b = init_bonds_[end_idx + 1];
        Real S = (P_a - P_b) / A;

        // Rebonato effective variance: σ_swap² T_ex = Σ_{i,j} w_i w_j ρ_ij σ_i σ_j T_ex
        // (假设 σ_i 为常数, 故 T_ex 因子提出)
        Real var_swap = 0.0;
        for (Size i = start_idx; i <= end_idx; ++i) {
            for (Size j = start_idx; j <= end_idx; ++j) {
                Real w_i = cfg_.tau(i) * init_bonds_[i + 1] / A;
                Real w_j = cfg_.tau(j) * init_bonds_[j + 1] / A;
                var_swap += w_i * w_j * cfg_.correlation[i][j]
                          * cfg_.volatilities[i] * cfg_.volatilities[j];
            }
        }
        // var_swap 已是 σ_swap² (单位时间), 乘 T_ex 得到方差
        Real total_var = var_swap * T_ex;
        Real sigma_swap = std::sqrt(total_var);

        if (sigma_swap < 1e-15) {
            Real intrinsic = is_payer ? std::max(S - K, 0.0) : std::max(K - S, 0.0);
            return N * A * intrinsic;
        }
        Real d1 = (std::log(S / K) + 0.5 * total_var) / sigma_swap;
        Real d2 = d1 - sigma_swap;
        Real price;
        if (is_payer) {
            price = N * A * (S * normal_cdf(d1) - K * normal_cdf(d2));
        } else {
            price = N * A * (K * normal_cdf(-d2) - S * normal_cdf(-d1));
        }
        return price;
    }

    // ============ 模拟演化 (Euler on log L) ============
    // 演化所有 LIBOR 到 T_horizon; n_steps 步长; 路径存储在 paths_[step][i]
    // 注: 模拟只覆盖 t ≤ T_horizon 内的 LIBOR; 已到期的 LIBOR 不再演化
    void simulate_path(Real T_horizon, Size n_steps, Philox4x64& rng) {
        if (T_horizon <= 0.0) throw std::invalid_argument("LMM::simulate: T_horizon must be positive");
        if (n_steps == 0) throw std::invalid_argument("LMM::simulate: n_steps must be positive");
        Size N = cfg_.n_libors();
        Real dt = T_horizon / static_cast<Real>(n_steps);

        // 当前 LIBOR state
        std::vector<Real> L = cfg_.initial_libors;
        paths_.assign(n_steps + 1, std::vector<Real>(N));
        paths_[0] = L;

        // Spot measure: 当前 active 的 LIBOR 索引 = 第一个未到期 (T_i ≥ t)
        Size active_idx = 0;

        for (Size step = 1; step <= n_steps; ++step) {
            Real t = static_cast<Real>(step) * dt;
            Real t_prev = static_cast<Real>(step - 1) * dt;

            // 更新 active_idx: 跳过已到期 LIBOR (T_i < t_prev)
            while (active_idx < N && cfg_.tenors[active_idx] < t_prev - 1e-12) {
                active_idx++;
            }

            // 生成相关 Brownian dW_i ~ N(0, dt)
            std::vector<Real> z(N, 0.0);
            std::vector<Real> dW(N, 0.0);
            for (Size i = 0; i < N; ++i) {
                z[i] = next_normal(rng) * std::sqrt(dt);
            }
            for (Size i = 0; i < N; ++i) {
                for (Size k = 0; k <= i; ++k) {
                    dW[i] += chol_L_[i][k] * z[k];
                }
            }

            // 演化所有 active LIBOR (i ≥ active_idx)
            for (Size i = active_idx; i < N; ++i) {
                Real L_i = L[i];
                Real sigma_i = cfg_.volatilities[i];
                if (sigma_i < 1e-15) continue;  // 零波动率不演化
                // drift μ_i (spot measure):
                //   μ_i = Σ_{j: active_idx ≤ j ≤ i} ρ_ij σ_j L_j τ_j / (1 + L_j τ_j)
                // (求和从 active_idx 到 i, 含 i; 上方 LIBOR 的 drift 为 0)
                Real mu_i = 0.0;
                if (cfg_.measure == "spot") {
                    for (Size j = active_idx; j <= i; ++j) {
                        Real L_j = L[j];
                        Real sigma_j = cfg_.volatilities[j];
                        Real tau_j = cfg_.tau(j);
                        mu_i += cfg_.correlation[i][j] * sigma_j * L_j * tau_j
                              / (1.0 + L_j * tau_j);
                    }
                } else {
                    // Terminal measure: μ_i = -Σ_{j > i} ρ_ij σ_j L_j τ_j / (1 + L_j τ_j)
                    for (Size j = i + 1; j < N; ++j) {
                        Real L_j = L[j];
                        Real sigma_j = cfg_.volatilities[j];
                        Real tau_j = cfg_.tau(j);
                        mu_i -= cfg_.correlation[i][j] * sigma_j * L_j * tau_j
                              / (1.0 + L_j * tau_j);
                    }
                }
                // Euler on log L: d log L_i = (μ_i - 0.5 σ_i²) dt + σ_i dW_i
                Real d_log_L = (mu_i - 0.5 * sigma_i * sigma_i) * dt + sigma_i * dW[i];
                L[i] = L_i * std::exp(d_log_L);
            }

            paths_[step] = L;
        }

        L_state_ = L;
        sim_dt_ = dt;
        sim_T_ = T_horizon;
        has_simulated_ = true;
    }

    // 模拟后查询: L_i(t_step)
    Real libor_at_step(Size step, Size i) const {
        if (!has_simulated_) throw std::runtime_error("LMM: must simulate first");
        if (step >= paths_.size()) throw std::out_of_range("LMM::libor_at_step: step");
        if (i >= cfg_.n_libors()) throw std::out_of_range("LMM::libor_at_step: i");
        return paths_[step][i];
    }

    // 模拟后查询: P(t_step, T) — 由模拟后的 LIBOR 重建
    // P(t, T_{k+1}) = Π_{i=k_start..k} 1/(1 + L_i(t) τ_i), 其中 k_start = 第一个 T_i > t
    Real zero_coupon_bond_at_step(Size step, Real T) const {
        if (!has_simulated_) throw std::runtime_error("LMM: must simulate first");
        if (step >= paths_.size()) throw std::out_of_range("LMM::bond_at_step: step");
        Real t = static_cast<Real>(step) * sim_dt_;
        if (T <= t) return 1.0;
        const auto& tenors = cfg_.tenors;
        const auto& L = paths_[step];
        // 找到 t 落在哪个 LIBOR 区间: t ∈ [T_k, T_{k+1}) → 从 LIBOR k 开始累积
        Size k_start = 0;
        while (k_start < cfg_.n_libors() && tenors[k_start + 1] <= t + 1e-12) {
            k_start++;
        }
        // 累积 P(t, T_k_start+1) = 1/(1 + L_{k_start}(t) * (T_{k_start+1} - t))
        // 简化: 假设 t 落在 tenor 节点上 (dt 与 tenor grid 对齐), 则 t = T_{k_start}
        Real P = 1.0;
        // 如果 t 在网格节点 T_{k_start} 上, 直接用 LIBOR 累积
        if (std::abs(t - tenors[k_start]) < 1e-9) {
            for (Size i = k_start; i < cfg_.n_libors() && tenors[i + 1] <= T + 1e-12; ++i) {
                P /= (1.0 + L[i] * cfg_.tau(i));
            }
            // 若 T 不在 tenor 节点上, 线性外推最后一段
            if (T > tenors.back()) {
                Real r_last = L.back();
                P *= std::exp(-r_last * (T - tenors.back()));
            }
            return P;
        }
        // t 不在节点: 简化处理, 用初始 LIBOR 估算
        // (此分支主要用于非对齐模拟, 测试中应保证对齐)
        Real partial_tau = tenors[k_start + 1] - t;
        if (k_start < cfg_.n_libors()) {
            P /= (1.0 + L[k_start] * partial_tau);
        }
        for (Size i = k_start + 1; i < cfg_.n_libors() && tenors[i + 1] <= T + 1e-12; ++i) {
            P /= (1.0 + L[i] * cfg_.tau(i));
        }
        return P;
    }

    const LMMConfig& config() const { return cfg_; }
    const std::vector<Real>& initial_bonds() const { return init_bonds_; }
    Real simulated_dt() const { return sim_dt_; }
    Real simulated_T() const { return sim_T_; }
    bool has_simulated() const { return has_simulated_; }
    const std::vector<Real>& current_libors() const { return L_state_; }

private:
    LMMConfig cfg_;
    std::vector<std::vector<Real>> chol_L_;  // Cholesky 下三角
    std::vector<Real> init_bonds_;           // P(0, T_0..T_N)
    std::vector<std::vector<Real>> paths_;   // paths_[step][i] = L_i(t_step)
    std::vector<Real> L_state_;              // 最后一步的 LIBOR state
    Real sim_dt_ = 0.0;
    Real sim_T_ = 0.0;
    bool has_simulated_ = false;

    static Real next_normal(Philox4x64& rng) {
        Real u1 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
        if (u1 < 1e-300) u1 = 1e-300;
        Real u2 = static_cast<Real>(rng()) / static_cast<Real>(UINT64_MAX);
        return std::sqrt(-2.0 * std::log(u1))
             * std::cos(6.283185307179586476925286766559 * u2);
    }

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
                    // 允许半正定 (s=0): 单因子全 1 矩阵秩亏时 L 对角元为 0
                    if (s < -1e-12) {
                        throw std::invalid_argument("LMM: correlation not positive semidefinite");
                    }
                    L[i][j] = std::sqrt(std::max(s, 0.0));
                } else {
                    // L[j][j]=0 时 (秩亏), L[i][j]=0 避免除零
                    L[i][j] = (L[j][j] > 1e-15) ? s / L[j][j] : 0.0;
                }
            }
        }
        return L;
    }
};

// ============ 便捷工厂: 平坦 LIBOR 曲线 + 单因子 (ρ=1) ============
inline LMMConfig make_flat_lmm(Real L_flat, Real T_max, Real dt_tenor,
                                  Real sigma_uniform, Size n_factors = 1) {
    LMMConfig cfg;
    for (Real T = 0.0; T <= T_max + 1e-10; T += dt_tenor) {
        cfg.tenors.push_back(T);
    }
    Size N = cfg.tenors.size() - 1;
    cfg.initial_libors.assign(N, L_flat);
    cfg.volatilities.assign(N, sigma_uniform);
    // 单因子: ρ_ij = 1
    cfg.correlation.assign(N, std::vector<Real>(N, 1.0));
    // 多因子: 简化用指数衰减相关性 ρ_ij = exp(-β |T_i - T_j|)
    if (n_factors > 1) {
        Real beta = 0.1;
        for (Size i = 0; i < N; ++i) {
            for (Size j = 0; j < N; ++j) {
                cfg.correlation[i][j] = std::exp(-beta * std::abs(cfg.tenors[i] - cfg.tenors[j]));
            }
        }
    }
    return cfg;
}

}  // namespace v1
}  // namespace cpphub
