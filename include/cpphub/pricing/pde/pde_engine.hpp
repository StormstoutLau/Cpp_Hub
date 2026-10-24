#pragma once
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <string>
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/pricing/pde/fdm_grid.hpp"
#include "cpphub/pricing/pde/fdm_scheme.hpp"
#include "cpphub/pricing/pde/thomas_solver.hpp"
#include "cpphub/instruments/payoff/payoff.hpp"

namespace cpphub {
inline namespace v1 {

// 边界条件类型
enum class BoundaryType {
    Dirichlet,  // 强加边界值 (V(S_min), V(S_max) 由解析公式或 payoff 给定)
    Neumann     // Γ=0 线性外推: V[0] = 2*V[1] - V[2], V[n-1] = 2*V[n-2] - V[n-3]
};

struct PDEEngineConfig {
    Size n_spatial = 400;
    Size n_time = 1000;
    Real alpha = 0.2;
    FDMSchemeType scheme = FDMSchemeType::CrankNicolson;
    Real s_multiplier = 5.0;
    // Rannacher smoothing 参数 (仅当 scheme = RannacherSmoothing 时生效)
    Size rannacher_warmup = 4;
    // 边界条件类型 (默认 Dirichlet 保持向后兼容)
    BoundaryType boundary = BoundaryType::Dirichlet;
    // RISK-006: PSOR 松弛因子配置
    // psor_omega > 0.0: 用户指定固定值 (默认 1.5 保持向后兼容)
    // psor_omega == 0.0: 自适应估计 (Gershgorin 上界 + Young 公式)
    Real psor_omega = 1.5;
    Size psor_max_iter = 2000;
    Real psor_tol = 1e-8;
};

class PDEEngine {
public:
    explicit PDEEngine(PDEEngineConfig config = PDEEngineConfig{})
        : config_(config) {}

    Real price_european(const PayOff& payoff, Real S0, Real K, Real T,
                         Real r, Real q, Real sigma) const {
        PDEParams params{r, q, sigma, T, K, S0};
        FDMGrid grid(config_.n_spatial, S0, K, sigma, T, config_.alpha);
        TimeGrid time(config_.n_time, T);
        Size n = grid.size();
        Size n_steps = config_.n_time;
        Real dt = time.dt();

        std::vector<Real> V(n), V_new(n);
        for (Size i = 0; i < n; ++i) {
            V[i] = payoff(grid.s(i));
        }

        auto scheme = create_scheme();

        for (Size step = 0; step < n_steps; ++step) {
            Real tau = static_cast<Real>(step + 1) * dt;
            V_new = V;
            set_boundary(V_new, tau, grid, payoff, r, q, K, sigma, false);
            scheme->step(V, V_new, dt, grid, params);
            finalize_boundary(V_new);
            V = V_new;
        }

        return interpolate_at_S0(V, grid, S0);
    }

    // RISK-030 (ADR-021 P1): 主入口返回 PricingResult — "本次调用的 PSOR 迭代数"
    // 随结果一起返回, 取代 mutable 共享态 (const 方法内写 mutable 成员, 同一
    // const PDEEngine 并发调用即 data race UB)。
    struct PricingResult {
        Real value;
        Size psor_iterations;
    };

    PricingResult price_american_result(const PayOff& payoff, Real S0, Real K, Real T,
                                          Real r, Real q, Real sigma) const {
        PDEParams params{r, q, sigma, T, K, S0};
        FDMGrid grid(config_.n_spatial, S0, K, sigma, T, config_.alpha);
        TimeGrid time(config_.n_time, T);
        Size n = grid.size();
        Size n_steps = config_.n_time;
        Real dt = time.dt();

        std::vector<Real> V(n), V_new(n);
        for (Size i = 0; i < n; ++i) {
            V[i] = payoff(grid.s(i));
        }

        // RISK-006: 用户指定 omega 或自适应估计
        Real tol = config_.psor_tol;
        Size max_iter = config_.psor_max_iter;
        Size total_iter = 0;

        for (Size step = 0; step < n_steps; ++step) {
            Real tau = static_cast<Real>(step + 1) * dt;

            std::vector<Real> payoff_vals(n);
            for (Size i = 0; i < n; ++i) {
                payoff_vals[i] = payoff(grid.s(i));
            }

            std::vector<Real> a, b, c, d;
            build_matrix(grid, params, dt, a, b, c);

            // RISK-006: 每步根据当前矩阵估计最优 omega (自适应模式)
            Real omega = (config_.psor_omega > 0.0)
                           ? config_.psor_omega
                           : estimate_optimal_omega(a, b, c);

            d.assign(n, 0.0);
            if (config_.scheme == FDMSchemeType::ImplicitEuler) {
                for (Size i = 1; i < n - 1; ++i) d[i] = V[i];
            } else {
                std::vector<Real> a_op, b_op, c_op;
                compute_operator_coeffs(grid, params, a_op, b_op, c_op);
                Real dt_half = dt * 0.5;
                for (Size i = 1; i < n - 1; ++i) {
                    d[i] = V[i] + dt_half * (a_op[i] * V[i - 1] +
                                              b_op[i] * V[i] +
                                              c_op[i] * V[i + 1]);
                }
            }

            V_new = V;
            set_boundary(V_new, tau, grid, payoff, r, q, K, sigma, true);
            total_iter += psor_solve(V_new, payoff_vals, a, b, c, d, omega, tol, max_iter);
            finalize_boundary(V_new);
            V = V_new;
        }

        return {interpolate_at_S0(V, grid, S0), total_iter};
    }

    // 旧签名包装 (ADR-021 P1): 只取 value, 数值路径与主入口完全一致
    Real price_american(const PayOff& payoff, Real S0, Real K, Real T,
                         Real r, Real q, Real sigma) const {
        return price_american_result(payoff, S0, K, T, r, q, sigma).value;
    }

    struct Greeks {
        Real delta;
        Real gamma;
        Real theta;
    };

    // 障碍侧 (H < S0 为 down barrier, H > S0 为 up barrier)
    enum class BarrierSide { Down, Up };

    // 障碍期权 FDM 定价 (连续监控, Out 期权)
    // 网格在障碍处截断, V[barrier] = 0 (Dirichlet)
    // 另一端 (远场) 边界按 payoff.kind() 分派 (ADR-021 P4):
    //   Call/Put → BSM 解析式 (线性 payoff, 含 S 项的折扣语义);
    //   DigitalCall/DigitalPut/DoubleDigital → payoff(S_far)·e^{-rτ} (远场常值);
    //   Unknown → throw std::invalid_argument (fail-fast, 不静默映射 BSM)。
    // 远场域宽倍数取 config_.s_multiplier (原硬编码 5.0, 默认行为不变;
    // P6(a) 测试锚依赖该配置实现默认域/宽域对比)。
    // BoundaryType::Neumann 仅作用于远场端 (ADR-021 P5): step 后线性外推,
    // 障碍端恒 Dirichlet 0 (外推会抹掉 V=0, 故不复用 finalize_boundary)。
    Real price_barrier(const PayOff& payoff, Real S0, Real K, Real H, Real T,
                        Real r, Real q, Real sigma, BarrierSide side) const {
        // RISK-031 (ADR-021 P4): Unknown kind 无法确定安全远场边界
        // (BSM 公式对非线性 payoff 不适用, name() 字符串分派可被冒名劫持)
        const PayOffKind kind = payoff.kind();
        if (kind == PayOffKind::Unknown) {
            throw std::invalid_argument(
                "price_barrier: payoff kind() is Unknown — override "
                "PayOff::kind() in your PayOff-derived class to select a "
                "far-field boundary (Call/Put map to BSM, Digital*/DoubleDigital "
                "map to discounted payoff)");
        }

        PDEParams params{r, q, sigma, T, K, S0};
        Real c = config_.alpha;
        Real m = config_.s_multiplier;
        Real s_far = S0 * std::exp(c * sigma * std::sqrt(T) * m);

        // 构建截断网格: 障碍端为 H, 另一端为正常远场边界
        Real s_min, s_max;
        if (side == BarrierSide::Down) {
            // Down barrier: H < S0, 网格 [H, s_far]
            if (H >= S0) throw std::invalid_argument("price_barrier: Down requires H < S0");
            s_min = H;
            s_max = s_far;
        } else {
            // Up barrier: H > S0, 网格 [s_far_inv, H]
            if (H <= S0) throw std::invalid_argument("price_barrier: Up requires H > S0");
            s_min = S0 * std::exp(-c * sigma * std::sqrt(T) * m);
            s_max = H;
        }

        FDMGrid grid(custom_boundary, config_.n_spatial, s_min, s_max, S0, config_.alpha);
        TimeGrid time(config_.n_time, T);
        Size n = grid.size();
        Size n_steps = config_.n_time;
        Real dt = time.dt();

        // 远场端 Dirichlet 边界值 (kind() 分派; kind 已守卫, default 只剩 Digital*)
        auto far_boundary = [&](Real S_far, Real tau) -> Real {
            switch (kind) {
                case PayOffKind::Call:
                    return bsm_call(S_far, K, tau, r, q, sigma);
                case PayOffKind::Put:
                    return bsm_put(S_far, K, tau, r, q, sigma);
                default:  // DigitalCall/DigitalPut/DoubleDigital: 远场常值 discount
                    return payoff(S_far) * std::exp(-r * tau);
            }
        };

        std::vector<Real> V(n), V_new(n);
        for (Size i = 0; i < n; ++i) {
            V[i] = payoff(grid.s(i));
            // 障碍侧 payoff 为 0 (Out 期权触及障碍即失效)
            if (side == BarrierSide::Down && grid.s(i) <= H) V[i] = 0.0;
            if (side == BarrierSide::Up && grid.s(i) >= H) V[i] = 0.0;
        }

        auto scheme = create_scheme();
        const bool neumann_far = (config_.boundary == BoundaryType::Neumann);

        for (Size step = 0; step < n_steps; ++step) {
            Real tau = static_cast<Real>(step + 1) * dt;
            V_new = V;

            // 障碍侧: V = 0 (Dirichlet, 恒定)
            // 远场侧: Dirichlet 时强加 kind() 分派边界; Neumann 时保持旧内部值
            // (step 后由远场外推更新 — 与 set_boundary 的 Neumann 语义对齐)
            if (side == BarrierSide::Down) {
                V_new[0] = 0.0;
                if (!neumann_far) V_new[n - 1] = far_boundary(grid.s_max(), tau);
            } else {
                V_new[n - 1] = 0.0;
                if (!neumann_far) V_new[0] = far_boundary(grid.s_min(), tau);
            }

            scheme->step(V, V_new, dt, grid, params);

            // Neumann 远场端后处理: Γ=0 线性外推 (方向按 Down/Up 定位远端)
            if (neumann_far) {
                if (side == BarrierSide::Down) {
                    V_new[n - 1] = 2.0 * V_new[n - 2] - V_new[n - 3];
                } else {
                    V_new[0] = 2.0 * V_new[1] - V_new[2];
                }
            }
            V = V_new;
        }

        return interpolate_at_S0(V, grid, S0);
    }

    Greeks greeks(const PayOff& payoff, Real S0, Real K, Real T,
                   Real r, Real q, Real sigma, bool american = false) const {
        Real h = 0.01;
        Real V0 = price_either(payoff, S0, K, T, r, q, sigma, american);
        Real Vp = price_either(payoff, S0 + h, K, T, r, q, sigma, american);
        Real Vm = price_either(payoff, S0 - h, K, T, r, q, sigma, american);

        Real delta = (Vp - Vm) / (2.0 * h);
        Real gamma = (Vp - 2.0 * V0 + Vm) / (h * h);

        Real dt_small = 1.0 / 365.0;
        Real V_forward = price_either(payoff, S0, K, T + dt_small, r, q, sigma, american);
        Real theta = -(V_forward - V0) / dt_small;

        return {delta, gamma, theta};
    }

    std::string name() const { return "PDEEngine"; }
    const PDEEngineConfig& config() const { return config_; }

    // RISK-006: Gershgorin 上界估计 Jacobi 谱半径, Young 公式计算最优 ω
    // 纯函数, 暴露为 public static 便于测试
    // ρ(B) ≤ max_i (|a_i/b_i| + |c_i/b_i|)
    // ω* = 2/(1+√(1-ρ²)), 裁剪到 [1.0, 1.95]
    static Real estimate_optimal_omega(const std::vector<Real>& a,
                                          const std::vector<Real>& b,
                                          const std::vector<Real>& c) {
        Size n = b.size();
        Real rho_upper = 0.0;
        for (Size i = 1; i < n - 1; ++i) {
            if (std::abs(b[i]) < 1e-15) continue;
            Real ratio = (std::abs(a[i]) + std::abs(c[i])) / std::abs(b[i]);
            rho_upper = std::max(rho_upper, ratio);
        }
        // 裁剪到 [0, 0.999] 避免 sqrt 负数
        rho_upper = std::min(rho_upper, 0.999);
        Real omega_opt = 2.0 / (1.0 + std::sqrt(1.0 - rho_upper * rho_upper));
        // 裁剪到安全范围 [1.0, 1.95]
        return std::max(1.0, std::min(omega_opt, 1.95));
    }

    // RISK-030 (ADR-021 P2): last_total_iterations() 与 mutable 成员已删除 —
    // "迭代数"语义改由 price_american_result() 返回值承载

private:
    PDEEngineConfig config_;

    Real price_either(const PayOff& payoff, Real S0, Real K, Real T,
                       Real r, Real q, Real sigma, bool american) const {
        return american ? price_american(payoff, S0, K, T, r, q, sigma)
                        : price_european(payoff, S0, K, T, r, q, sigma);
    }

    static Real bsm_call(Real S, Real K, Real tau, Real r, Real q, Real sigma) {
        if (tau <= 0.0) return std::max(S - K, 0.0);
        Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * tau) / (sigma * std::sqrt(tau));
        Real d2 = d1 - sigma * std::sqrt(tau);
        return S * std::exp(-q * tau) * normal_cdf(d1) - K * std::exp(-r * tau) * normal_cdf(d2);
    }

    static Real bsm_put(Real S, Real K, Real tau, Real r, Real q, Real sigma) {
        if (tau <= 0.0) return std::max(K - S, 0.0);
        Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * tau) / (sigma * std::sqrt(tau));
        Real d2 = d1 - sigma * std::sqrt(tau);
        return K * std::exp(-r * tau) * normal_cdf(-d2) - S * std::exp(-q * tau) * normal_cdf(-d1);
    }

    void set_boundary(std::vector<Real>& V, Real tau, const FDMGrid& grid,
                       const PayOff& payoff, Real r, Real q, Real K, Real sigma,
                       bool) const {
        Size n = V.size();
        if (config_.boundary == BoundaryType::Neumann) {
            // Neumann: Γ=0 线性外推 (初始边界用旧内部值, step 后由 finalize_boundary 更新)
            V[0] = 2.0 * V[1] - V[2];
            V[n - 1] = 2.0 * V[n - 2] - V[n - 3];
            return;
        }
        // Dirichlet: 解析边界 (RISK-031 ADR-021 P4: 识别从 name() 换成 kind(),
        // else 分支 payoff(S)·e^{-rτ} 对 Digital*/DoubleDigital/Unknown 依然适用,
        // 数值行为与旧实现逐位一致, 不 throw)
        Real S_min = grid.s_min();
        Real S_max = grid.s_max();
        const PayOffKind kind = payoff.kind();
        if (kind == PayOffKind::Call) {
            V[0] = bsm_call(S_min, K, tau, r, q, sigma);
            V[n - 1] = bsm_call(S_max, K, tau, r, q, sigma);
        } else if (kind == PayOffKind::Put) {
            V[0] = bsm_put(S_min, K, tau, r, q, sigma);
            V[n - 1] = bsm_put(S_max, K, tau, r, q, sigma);
        } else {
            Real df = std::exp(-r * tau);
            V[0] = payoff(S_min) * df;
            V[n - 1] = payoff(S_max) * df;
        }
    }

    // Neumann 边界后处理: step 后用新内部值外推边界
    void finalize_boundary(std::vector<Real>& V) const {
        if (config_.boundary != BoundaryType::Neumann) return;
        Size n = V.size();
        if (n < 4) return;
        V[0] = 2.0 * V[1] - V[2];
        V[n - 1] = 2.0 * V[n - 2] - V[n - 3];
    }

    // RISK-006: 返回实际迭代次数 (用于性能对比)
    Size psor_solve(std::vector<Real>& V, const std::vector<Real>& payoff,
                     const std::vector<Real>& a, const std::vector<Real>& b,
                     const std::vector<Real>& c, const std::vector<Real>& d,
                     Real omega, Real tol, Size max_iter) const {
        Size n = V.size();
        for (Size iter = 0; iter < max_iter; ++iter) {
            Real max_diff = 0.0;
            for (Size i = 1; i < n - 1; ++i) {
                Real old_val = V[i];
                Real y = (d[i] - a[i] * V[i - 1] - c[i] * V[i + 1]) / b[i];
                V[i] = std::max(payoff[i], old_val + omega * (y - old_val));
                max_diff = std::max(max_diff, std::abs(V[i] - old_val));
            }
            if (max_diff < tol) return iter + 1;
        }
        return max_iter;
    }

    void build_matrix(const FDMGrid& grid, const PDEParams& params, Real dt,
                       std::vector<Real>& a, std::vector<Real>& b,
                       std::vector<Real>& c) const {
        Size n = grid.size();
        std::vector<Real> a_op, b_op, c_op;
        compute_operator_coeffs(grid, params, a_op, b_op, c_op);

        a.assign(n, 0.0); b.assign(n, 0.0); c.assign(n, 0.0);

        Real theta = (config_.scheme == FDMSchemeType::ImplicitEuler) ? 1.0 : 0.5;
        Real dt_theta = dt * theta;

        for (Size i = 1; i < n - 1; ++i) {
            a[i] = -dt_theta * a_op[i];
            b[i] = 1.0 - dt_theta * b_op[i];
            c[i] = -dt_theta * c_op[i];
        }

        b[0] = 1.0; c[0] = 0.0;
        a[n - 1] = 0.0; b[n - 1] = 1.0;
    }

    std::unique_ptr<FDMScheme> create_scheme() const {
        switch (config_.scheme) {
            case FDMSchemeType::ExplicitEuler:
                return std::make_unique<ExplicitEuler>();
            case FDMSchemeType::ImplicitEuler:
                return std::make_unique<ImplicitEuler>();
            case FDMSchemeType::RannacherSmoothing:
                return std::make_unique<RannacherSmoothing>(config_.rannacher_warmup, 0.5);
            case FDMSchemeType::CrankNicolson:
            default:
                return std::make_unique<CrankNicolson>(0.5);
        }
    }

    Real interpolate_at_S0(const std::vector<Real>& V, const FDMGrid& grid, Real S0) const {
        Size n = grid.size();
        for (Size i = 0; i < n - 1; ++i) {
            if (grid.s(i) <= S0 && S0 <= grid.s(i + 1)) {
                Real w = (S0 - grid.s(i)) / (grid.s(i + 1) - grid.s(i));
                return (1.0 - w) * V[i] + w * V[i + 1];
            }
        }
        return V[n - 1];
    }
};

}  // namespace v1
}  // namespace cpphub
