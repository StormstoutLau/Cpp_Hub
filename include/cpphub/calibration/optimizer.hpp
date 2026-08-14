#pragma once
// SOURCE: PHASE3_SPEC §4.1 - Optimization algorithms for calibration
// Implemented on main station (MSVC) - 2026-07-30
// LM: numerical Jacobian + damped normal equations + Gauss-Jordan solve
// NelderMead: standard simplex (reflect/expand/contract/shrink)
// DE: rand/1/bin mutation + binomial crossover
// NOTE: uses inline LCG instead of <random> to avoid MSVC ICE (C1001) in <random>
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {

using ObjectiveFn = std::function<Real(const std::vector<Real>&)>;
using GradientFn = std::function<std::vector<Real>(const std::vector<Real>&)>;
using ResidualFn = std::function<std::vector<Real>(const std::vector<Real>&)>;

// 约束函数: 返回约束值向量
// 不等式约束: c_i(x) >= 0 (与 scipy f_ieqcons 一致)
// 等式约束: c_i(x) = 0 (与 scipy f_eqcons 一致)
using ConstraintFn = std::function<std::vector<Real>(const std::vector<Real>&)>;

// 约束 Jacobian: 返回 (m, n) 矩阵, J[i][j] = dc_i/dx_j
// 可选: 用户未提供时用数值差分估计
using ConstraintJacobianFn = std::function<std::vector<std::vector<Real>>(const std::vector<Real>&)>;

struct OptimizationResult {
    std::vector<Real> x;
    Real fx;
    Size n_iterations;
    Size n_function_evaluations;
    bool converged;
    std::string message;
};

struct Bounds {
    Real lower;
    Real upper;
};

// Calibration common types (moved here from calibrator.hpp to break circular deps)
struct CalibrationResult {
    std::vector<Real> params;
    Real objective_value;
    std::vector<std::vector<Real>> jacobian;
    std::vector<std::vector<Real>> covariance;
    Size n_iterations;
    bool converged;
    std::string message;
    std::vector<Real> residuals;
};

struct CalibConfig {
    Size de_pop_size = 50;
    Size de_generations = 200;
    Size lm_max_iter = 200;
    Real ftol = 1e-8;
    Real xtol = 1e-8;
    uint64_t seed = 42;
    bool use_de_init = true;
    bool compute_diagnostics = true;
    // --- v1.1 校准稳定性增强 (Task 5) ---
    // Tikhonov 正则化: 目标函数加 0.5 * lambda_reg * ||x - params_prior||^2
    // 当 params_prior 为空或 lambda_reg <= 0 时不启用正则化
    Real lambda_reg = 0.0;
    std::vector<Real> params_prior;
    // 早停: 当残差 RMSE = sqrt(2*fx / m) < early_stop_rmse 时停止迭代
    // 0 表示不启用早停。实务中可设为 bid-ask spread 的一半
    Real early_stop_rmse = 0.0;
};

namespace detail {

// Gauss-Jordan elimination with partial pivoting for small dense systems.
// Solves A x = b in place. A is n×n row-major, b is n. Returns true if solvable.
inline bool solve_linear_system(std::vector<std::vector<Real>>& A,
                                std::vector<Real>& b,
                                Size n) {
    for (Size k = 0; k < n; ++k) {
        // Partial pivot
        Size piv = k;
        Real maxv = std::abs(A[k][k]);
        for (Size i = k + 1; i < n; ++i) {
            Real v = std::abs(A[i][k]);
            if (v > maxv) { maxv = v; piv = i; }
        }
        if (maxv < 1e-15) return false;  // singular
        if (piv != k) {
            std::swap(A[piv], A[k]);
            std::swap(b[piv], b[k]);
        }
        // Eliminate
        Real akk = A[k][k];
        for (Size i = k + 1; i < n; ++i) {
            Real f = A[i][k] / akk;
            if (f == 0.0) continue;
            for (Size j = k; j < n; ++j) A[i][j] -= f * A[k][j];
            b[i] -= f * b[k];
        }
    }
    // Back substitution
    for (Size i = n; i-- > 0;) {
        Real s = b[i];
        for (Size j = i + 1; j < n; ++j) s -= A[i][j] * b[j];
        b[i] = s / A[i][i];
    }
    return true;
}

// Numerical Jacobian via central differences for residual vector r(x).
// J[i][j] = dr_i / dx_j, shape (m, n)
// Central difference: O(h^2) accuracy, less cancellation error than forward diff.
// Note: ConstraintFn is an alias for ResidualFn, so this also serves constraints.
inline std::vector<std::vector<Real>> numerical_jacobian(
        const ResidualFn& r,
        const std::vector<Real>& x,
        Real eps = 1e-6) {
    Size n = x.size();
    auto r0 = r(x);
    Size m = r0.size();
    std::vector<std::vector<Real>> J(m, std::vector<Real>(n, 0.0));
    for (Size j = 0; j < n; ++j) {
        std::vector<Real> xp = x, xm = x;
        Real h = eps * (1.0 + std::abs(x[j]));
        xp[j] += h;
        xm[j] -= h;
        auto rp = r(xp);
        auto rm = r(xm);
        Real inv_2h = 0.5 / h;
        for (Size i = 0; i < m; ++i) J[i][j] = (rp[i] - rm[i]) * inv_2h;
    }
    return J;
}

// Numerical gradient via central differences for scalar objective f(x).
// g[j] = (f(x + h*e_j) - f(x - h*e_j)) / (2*h)
// O(h^2) accuracy, more precise than scipy default forward difference.
inline std::vector<Real> numerical_gradient(
        const ObjectiveFn& f,
        const std::vector<Real>& x,
        Real eps = 1e-6) {
    Size n = x.size();
    std::vector<Real> g(n, 0.0);
    for (Size j = 0; j < n; ++j) {
        std::vector<Real> xp = x, xm = x;
        Real h = eps * (1.0 + std::abs(x[j]));
        xp[j] += h;
        xm[j] -= h;
        g[j] = (f(xp) - f(xm)) / (2.0 * h);
    }
    return g;
}

// Solve equality-constrained QP via KKT system:
//   min 0.5*d^T B d + g^T d  s.t.  A*d = b
// KKT:
//   [B  A^T] [d]   [-g]
//   [A  0  ] [l] = [ b]
// Returns true if KKT matrix is non-singular.
inline bool solve_eqp_kkt(
        const std::vector<std::vector<Real>>& B,
        const std::vector<Real>& g,
        const std::vector<std::vector<Real>>& A,  // (m_eq, n)
        const std::vector<Real>& b,                // (m_eq,)
        std::vector<Real>& d,
        std::vector<Real>& lambda) {
    Size n = g.size();
    Size m = b.size();

    if (m == 0) {
        // Unconstrained: B*d = -g
        std::vector<std::vector<Real>> Bcopy = B;
        std::vector<Real> rhs(n);
        for (Size i = 0; i < n; ++i) rhs[i] = -g[i];
        if (!solve_linear_system(Bcopy, rhs, n)) return false;
        d = rhs;
        lambda.clear();
        return true;
    }

    // Build (n+m) x (n+m) KKT system
    Size dim = n + m;
    std::vector<std::vector<Real>> K(dim, std::vector<Real>(dim, 0.0));
    std::vector<Real> rhs(dim, 0.0);

    // Top-left n x n: B
    for (Size i = 0; i < n; ++i)
        for (Size j = 0; j < n; ++j)
            K[i][j] = B[i][j];
    // Top-right n x m: A^T
    for (Size i = 0; i < n; ++i)
        for (Size j = 0; j < m; ++j)
            K[i][n + j] = A[j][i];
    // Bottom-left m x n: A
    for (Size i = 0; i < m; ++i)
        for (Size j = 0; j < n; ++j)
            K[n + i][j] = A[i][j];
    // Bottom-right m x m: 0 (already zero)
    // RHS: [-g; b]
    for (Size i = 0; i < n; ++i) rhs[i] = -g[i];
    for (Size i = 0; i < m; ++i) rhs[n + i] = b[i];

    if (!solve_linear_system(K, rhs, dim)) return false;

    d.assign(rhs.begin(), rhs.begin() + n);
    lambda.assign(rhs.begin() + n, rhs.end());
    return true;
}

// Active-set QP solver for small problems (n<10, m<20).
//   min 0.5*d^T B d + g^T d
//   s.t.  A_ineq[i]*d >= b_ineq[i]   (m_ineq inequalities)
//         A_eq[i]*d   =  b_eq[i]     (m_eq equalities)
// Returns true if a solution was found; fills d and multipliers.
// Note: caller must ensure b_ineq encodes the linearized constraint violation
//       such that d=0 is feasible when b_ineq[i] <= 0 (i.e. constraint satisfied at x).
inline bool solve_qp_active_set(
        const std::vector<std::vector<Real>>& B,
        const std::vector<Real>& g,
        const std::vector<std::vector<Real>>& A_ineq,
        const std::vector<Real>& b_ineq,
        const std::vector<std::vector<Real>>& A_eq,
        const std::vector<Real>& b_eq,
        std::vector<Real>& d,
        std::vector<Real>& lambda_ineq,
        std::vector<Real>& lambda_eq,
        Size max_iter = 50) {
    Size n = g.size();
    Size m_ineq = b_ineq.size();
    Size m_eq = b_eq.size();

    d.assign(n, 0.0);
    lambda_ineq.assign(m_ineq, 0.0);
    lambda_eq.assign(m_eq, 0.0);

    if (m_ineq == 0 && m_eq == 0) {
        // Unconstrained QP: B*d = -g
        std::vector<std::vector<Real>> Bcopy = B;
        std::vector<Real> rhs(n);
        for (Size i = 0; i < n; ++i) rhs[i] = -g[i];
        if (!solve_linear_system(Bcopy, rhs, n)) return false;
        d = rhs;
        return true;
    }

    // Initial working set: active constraints at d=0 (b_ineq[i] >= -tol means active/violated)
    std::vector<bool> active(m_ineq, false);
    for (Size i = 0; i < m_ineq; ++i) {
        if (b_ineq[i] >= -1e-12) active[i] = true;
    }

    for (Size iter = 0; iter < max_iter; ++iter) {
        // Assemble active constraint matrix (equalities + active inequalities)
        std::vector<Size> active_idx;
        for (Size i = 0; i < m_ineq; ++i) {
            if (active[i]) active_idx.push_back(i);
        }
        Size m_active = m_eq + active_idx.size();
        std::vector<std::vector<Real>> A_all(m_active, std::vector<Real>(n, 0.0));
        std::vector<Real> b_all(m_active, 0.0);
        for (Size i = 0; i < m_eq; ++i) {
            A_all[i] = A_eq[i];
            b_all[i] = b_eq[i];
        }
        for (Size k = 0; k < active_idx.size(); ++k) {
            A_all[m_eq + k] = A_ineq[active_idx[k]];
            b_all[m_eq + k] = b_ineq[active_idx[k]];
        }

        std::vector<Real> d_new;
        std::vector<Real> lam_all;
        if (!solve_eqp_kkt(B, g, A_all, b_all, d_new, lam_all)) {
            return false;
        }

        // Compute step length alpha along (d_new - d)
        Real alpha = 1.0;
        Size blocking = m_ineq;  // invalid index
        bool has_blocking = false;
        for (Size i = 0; i < m_ineq; ++i) {
            if (active[i]) continue;
            // Constraint: A_ineq[i]*d >= b_ineq[i]
            // Along p = d_new - d: A_ineq[i]*(d + alpha*p) >= b_ineq[i]
            Real ad_dir = 0.0;   // A_ineq[i]*p
            Real ad_curr = 0.0;  // A_ineq[i]*d
            for (Size j = 0; j < n; ++j) {
                ad_dir += A_ineq[i][j] * (d_new[j] - d[j]);
                ad_curr += A_ineq[i][j] * d[j];
            }
            if (ad_dir < -1e-15) {
                // Direction violates this constraint; compute max step
                Real alpha_i = (b_ineq[i] - ad_curr) / ad_dir;
                if (alpha_i < alpha) {
                    alpha = alpha_i;
                    blocking = i;
                    has_blocking = true;
                }
            }
        }
        if (alpha < 0.0) alpha = 0.0;

        // Update d
        for (Size j = 0; j < n; ++j) d[j] += alpha * (d_new[j] - d[j]);

        if (has_blocking) {
            active[blocking] = true;
            continue;
        }

        // No blocking constraint: check multipliers of active inequalities.
        // KKT convention used here: L = f + λ^T*(A*d - b), so for inequality
        // A_ineq*d >= b_ineq the optimal multiplier must satisfy λ <= 0.
        // If any active inequality has λ > 0, remove it (constraint should be inactive).
        Real most_positive = 1e-10;
        Size remove_idx = m_ineq;
        for (Size k = 0; k < active_idx.size(); ++k) {
            Real lam = lam_all[m_eq + k];
            if (lam > most_positive) {
                most_positive = lam;
                remove_idx = active_idx[k];
            }
        }
        if (remove_idx == m_ineq) {
            // All inequality multipliers <= 0: optimal
            for (Size k = 0; k < active_idx.size(); ++k) {
                lambda_ineq[active_idx[k]] = lam_all[m_eq + k];
            }
            for (Size i = 0; i < m_eq; ++i) lambda_eq[i] = lam_all[i];
            return true;
        }
        active[remove_idx] = false;
    }
    return false;  // did not converge
}

// Damped BFGS update (Nocedal-Wright Eq 18.15).
// Ensures B stays positive definite even when curvature condition s^T y > 0 fails.
//   s = x_{k+1} - x_k,  y = gradL_{k+1} - gradL_k
//   sy = s^T y,  sBs = s^T B s
//   if sy >= 0.2 * sBs: theta = 1
//   else:               theta = 0.8 * sBs / (sBs - sy)
//   r = theta*y + (1-theta)*B*s
//   B <- B + r*r^T/(r^T s) - (B*s)*(B*s)^T/(s^T B s)
inline void bfgs_damped_update(
        std::vector<std::vector<Real>>& B,
        const std::vector<Real>& s,
        const std::vector<Real>& y) {
    Size n = s.size();
    Real sy = 0.0;
    for (Size i = 0; i < n; ++i) sy += s[i] * y[i];

    std::vector<Real> Bs(n, 0.0);
    for (Size i = 0; i < n; ++i)
        for (Size j = 0; j < n; ++j)
            Bs[i] += B[i][j] * s[j];

    Real sBs = 0.0;
    for (Size i = 0; i < n; ++i) sBs += s[i] * Bs[i];

    if (sBs < 1e-15) return;  // degenerate, skip update

    Real theta;
    if (sy >= 0.2 * sBs) {
        theta = 1.0;
    } else {
        Real denom = sBs - sy;
        if (std::abs(denom) < 1e-15) return;
        theta = 0.8 * sBs / denom;
    }

    std::vector<Real> r(n);
    for (Size i = 0; i < n; ++i) r[i] = theta * y[i] + (1.0 - theta) * Bs[i];

    Real rs = 0.0;
    for (Size i = 0; i < n; ++i) rs += r[i] * s[i];
    if (std::abs(rs) < 1e-15) return;

    for (Size i = 0; i < n; ++i) {
        for (Size j = 0; j < n; ++j) {
            B[i][j] += r[i] * r[j] / rs - Bs[i] * Bs[j] / sBs;
        }
    }
}

// L1 merit function: f(x) + sum mu_i*|c_eq_i| + sum mu_i*max(0, -c_ineq_i)
inline Real l1_merit(
        const ObjectiveFn& f,
        const std::vector<Real>& x,
        const std::vector<ConstraintFn>& ineq_constraints,
        const std::vector<ConstraintFn>& eq_constraints,
        const std::vector<Real>& mu_ineq,
        const std::vector<Real>& mu_eq) {
    Real merit = f(x);
    for (Size i = 0; i < ineq_constraints.size(); ++i) {
        auto c = ineq_constraints[i](x);
        for (Size k = 0; k < c.size(); ++k) {
            if (c[k] < 0.0) merit += mu_ineq[i] * (-c[k]);
        }
    }
    for (Size i = 0; i < eq_constraints.size(); ++i) {
        auto c = eq_constraints[i](x);
        for (Size k = 0; k < c.size(); ++k) {
            merit += mu_eq[i] * std::abs(c[k]);
        }
    }
    return merit;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Levenberg-Marquardt: minimize 0.5 * sum r_i(x)^2
// Iteration: (J^T J + lambda I) dx = -J^T r
// ---------------------------------------------------------------------------
class LevenbergMarquardt {
public:
    struct Config {
        Size max_iterations = 200;
        Real ftol = 1e-8;
        Real xtol = 1e-8;
        Real gtol = 1e-8;
        Real lambda_init = 1e-3;
        Real lambda_up = 10.0;
        Real lambda_down = 0.1;
        // --- v1.1 Tikhonov 正则化 (Task 5) ---
        // 目标函数: 0.5 * sum(r_i^2) + 0.5 * lambda_reg * ||x - params_prior||^2
        // 实现方式: 扩展残差向量 r_ext = [r, sqrt(lambda_reg)*(x - prior)]
        // 当 params_prior.size() != x0.size() 或 lambda_reg <= 0 时不启用
        Real lambda_reg = 0.0;
        std::vector<Real> params_prior;
        // 早停: RMSE = sqrt(2*fx_orig / m_orig) < early_stop_rmse 时停止
        Real early_stop_rmse = 0.0;
    };

    // GCC: 在 class 内使用 Config{} 作为默认参数会触发 "default member initializer
    // required before the end of its enclosing class" 错误,因为嵌套 struct Config
    // 在外层 class LevenbergMarquardt 定义结束前不算 complete。
    // 解决方案:类内声明不带默认参数,类外定义带默认参数。
    static OptimizationResult minimize(
            const ResidualFn& residuals,
            const std::vector<Real>& x0,
            const Config& cfg);
};

inline OptimizationResult LevenbergMarquardt::minimize(
        const ResidualFn& residuals,
        const std::vector<Real>& x0,
        const Config& cfg = LevenbergMarquardt::Config{}) {
        // Default argument defined out-of-class to satisfy GCC completeness rules
        // (nested Config is complete only after the enclosing class definition ends)
        OptimizationResult result{};  // RISK-014: 显式零初始化,避免 converged 字段为垃圾值
        std::vector<Real> x = x0;
        Size n = x.size();
        Real lambda = cfg.lambda_init;
        Size n_evals = 0;

        // v1.1: Tikhonov 正则化 — 通过扩展残差向量实现
        // r_ext = [r_orig(x), sqrt(lambda_reg) * (x - prior)]
        // 这样 J_ext = [J_orig; sqrt(lambda_reg)*I], J^T J 自动加 lambda_reg*I,
        // J^T r 自动加 lambda_reg*(x - prior), 严格等价于 Tikhonov 正则化 LM
        bool use_reg = (cfg.lambda_reg > 0.0) &&
                       (cfg.params_prior.size() == n);
        Real sqrt_lambda_reg = use_reg ? std::sqrt(cfg.lambda_reg) : 0.0;

        // 扩展残差函数: 在原始残差后追加 n 个正则化项
        ResidualFn r_ext = [&residuals, use_reg, sqrt_lambda_reg, &cfg, n](
                const std::vector<Real>& xx) -> std::vector<Real> {
            auto r = residuals(xx);
            if (!use_reg) return r;
            Size m = r.size();
            std::vector<Real> extended(m + n);
            for (Size i = 0; i < m; ++i) extended[i] = r[i];
            for (Size i = 0; i < n; ++i) {
                extended[m + i] = sqrt_lambda_reg * (xx[i] - cfg.params_prior[i]);
            }
            return extended;
        };

        // 原始 cost (用于早停 RMSE 判断, 不含正则化项)
        auto compute_cost_orig = [&](const std::vector<Real>& xx) -> Real {
            auto r = residuals(xx);
            ++n_evals;
            Real s = 0.0;
            for (Real v : r) s += v * v;
            return 0.5 * s;
        };
        // 扩展 cost (实际优化目标, 含正则化项)
        auto compute_cost_ext = [&](const std::vector<Real>& xx) -> Real {
            auto r = r_ext(xx);
            ++n_evals;
            Real s = 0.0;
            for (Real v : r) s += v * v;
            return 0.5 * s;
        };

        Real fx = compute_cost_ext(x);
        Real fx_orig = use_reg ? compute_cost_orig(x) : fx;
        result.n_iterations = 0;

        // 早停检查: RMSE = sqrt(2*fx_orig / m_orig) < early_stop_rmse
        auto check_early_stop = [&](Real fx_o) -> bool {
            if (cfg.early_stop_rmse <= 0.0) return false;
            auto r = residuals(x);  // 不计 n_evals (复用)
            Size m = r.size();
            if (m == 0) return false;
            Real rmse = std::sqrt(2.0 * fx_o / static_cast<Real>(m));
            return rmse < cfg.early_stop_rmse;
        };

        for (Size iter = 0; iter < cfg.max_iterations; ++iter) {
            result.n_iterations = iter + 1;
            auto r = r_ext(x); ++n_evals;
            auto J = detail::numerical_jacobian(r_ext, x);
            Size m = r.size();

            // Build J^T J (n×n) and J^T r (n)
            std::vector<std::vector<Real>> JtJ(n, std::vector<Real>(n, 0.0));
            std::vector<Real> Jtr(n, 0.0);
            for (Size i = 0; i < m; ++i) {
                for (Size j = 0; j < n; ++j) {
                    Jtr[j] += J[i][j] * r[i];
                    for (Size k = 0; k < n; ++k) JtJ[j][k] += J[i][j] * J[i][k];
                }
            }

            // Gradient norm convergence
            Real gnorm = 0.0;
            for (Real g : Jtr) gnorm += g * g;
            gnorm = std::sqrt(gnorm);
            if (gnorm < cfg.gtol) {
                result.converged = true;
                result.message = "gtol satisfied";
                // v1.1: 早停优先报告 (若 RMSE 已低于阈值, 报告 early_stop 而非 gtol)
                if (cfg.early_stop_rmse > 0.0) {
                    Real fx_o = use_reg ? compute_cost_orig(x) : fx;
                    if (check_early_stop(fx_o)) {
                        result.message = "early_stop_rmse satisfied";
                    }
                }
                break;
            }

            // Damped solve: (J^T J + lambda * I) dx = -J^T r  (Levenberg scaling)
            bool step_accepted = false;
            for (Size inner = 0; inner < 30; ++inner) {
                auto A = JtJ;
                auto b = Jtr;
                for (Size i = 0; i < n; ++i) {
                    A[i][i] += lambda;
                    b[i] = -b[i];
                }
                std::vector<Real> dx;
                if (detail::solve_linear_system(A, b, n)) {
                    dx = b;
                } else {
                    lambda *= cfg.lambda_up;
                    continue;
                }
                std::vector<Real> x_new(n);
                for (Size i = 0; i < n; ++i) x_new[i] = x[i] + dx[i];
                Real fx_new = compute_cost_ext(x_new);
                if (fx_new < fx) {
                    Real dxnorm = 0.0;
                    for (Real d : dx) dxnorm += d * d;
                    dxnorm = std::sqrt(dxnorm);
                    Real rel_impr = (fx - fx_new) / std::max(fx, 1e-30);
                    x = x_new;
                    fx = fx_new;
                    if (use_reg) fx_orig = compute_cost_orig(x);
                    lambda = std::max(lambda * cfg.lambda_down, 1e-12);
                    step_accepted = true;
                    if (dxnorm < cfg.xtol) {
                        result.converged = true;
                        result.message = "xtol satisfied";
                    }
                    if (rel_impr < cfg.ftol) {
                        result.converged = true;
                        result.message = "ftol satisfied";
                    }
                    // v1.1: 早停检查 (基于原始 RMSE)
                    if (cfg.early_stop_rmse > 0.0 && check_early_stop(fx_orig)) {
                        result.converged = true;
                        result.message = "early_stop_rmse satisfied";
                    }
                    break;
                } else {
                    lambda *= cfg.lambda_up;
                    if (lambda > 1e12) break;
                }
            }
            if (!step_accepted) {
                result.converged = false;
                result.message = "no progress (lambda exploded)";
                break;
            }
            if (result.converged) break;
        }

        result.x = x;
        result.fx = fx;
        result.n_function_evaluations = n_evals;
        if (result.message.empty()) {
            result.converged = false;
            result.message = "max iterations reached";
        }
        return result;
}

// ---------------------------------------------------------------------------
// Nelder-Mead simplex (for non-smooth or noisy objectives)
// ---------------------------------------------------------------------------
class NelderMead {
public:
    struct Config {
        Size max_iterations = 500;
        Real ftol = 1e-10;
        Real xtol = 1e-10;
        Real alpha = 1.0;   // reflection
        Real gamma = 2.0;   // expansion
        Real rho = 0.5;     // contraction
        Real sigma = 0.5;   // shrink
    };

    static OptimizationResult minimize(
            ObjectiveFn f,
            const std::vector<Real>& x0,
            const Config& cfg);
};

inline OptimizationResult NelderMead::minimize(
        ObjectiveFn f,
        const std::vector<Real>& x0,
        const Config& cfg = NelderMead::Config{}) {
        OptimizationResult result{};  // RISK-014: 显式零初始化
        Size n = x0.size();
        if (n == 0) {
            result.converged = true;
            result.message = "empty x0";
            return result;
        }

        // Build initial simplex: x0 plus n perturbed vertices
        std::vector<std::vector<Real>> simplex(n + 1, x0);
        std::vector<Real> fval(n + 1);
        fval[0] = f(simplex[0]);
        for (Size i = 0; i < n; ++i) {
            Real step = (x0[i] != 0.0) ? 0.05 * x0[i] : 0.00025;
            simplex[i + 1][i] += step;
            fval[i + 1] = f(simplex[i + 1]);
        }
        Size n_evals = n + 1;
        result.n_iterations = 0;

        for (Size iter = 0; iter < cfg.max_iterations; ++iter) {
            result.n_iterations = iter + 1;
            // Sort by fval ascending
            std::vector<Size> idx(n + 1);
            for (Size i = 0; i <= n; ++i) idx[i] = i;
            std::sort(idx.begin(), idx.end(), [&](Size a, Size b) { return fval[a] < fval[b]; });

            // Convergence: range of fval
            Real frange = fval[idx[n]] - fval[idx[0]];
            if (frange < cfg.ftol) {
                result.converged = true;
                result.message = "ftol satisfied";
                break;
            }
            // Convergence: max vertex distance from best
            Real xdist = 0.0;
            for (Size i = 1; i <= n; ++i) {
                for (Size j = 0; j < n; ++j) {
                    Real d = std::abs(simplex[idx[i]][j] - simplex[idx[0]][j]);
                    if (d > xdist) xdist = d;
                }
            }
            if (xdist < cfg.xtol) {
                result.converged = true;
                result.message = "xtol satisfied";
                break;
            }

            // Centroid of all but worst
            std::vector<Real> centroid(n, 0.0);
            for (Size i = 0; i < n; ++i) {
                for (Size j = 0; j < n; ++j) centroid[j] += simplex[idx[i]][j];
            }
            for (Real& c : centroid) c /= static_cast<Real>(n);

            // Reflection
            std::vector<Real> xr(n);
            for (Size j = 0; j < n; ++j)
                xr[j] = centroid[j] + cfg.alpha * (centroid[j] - simplex[idx[n]][j]);
            Real fr = f(xr); ++n_evals;

            if (fr >= fval[idx[0]] && fr < fval[idx[n - 1]]) {
                simplex[idx[n]] = xr;
                fval[idx[n]] = fr;
            } else if (fr < fval[idx[0]]) {
                // Expansion
                std::vector<Real> xe(n);
                for (Size j = 0; j < n; ++j)
                    xe[j] = centroid[j] + cfg.gamma * (xr[j] - centroid[j]);
                Real fe = f(xe); ++n_evals;
                if (fe < fr) {
                    simplex[idx[n]] = xe;
                    fval[idx[n]] = fe;
                } else {
                    simplex[idx[n]] = xr;
                    fval[idx[n]] = fr;
                }
            } else {
                // Contraction
                std::vector<Real> xc(n);
                if (fr < fval[idx[n]]) {
                    // outside contraction
                    for (Size j = 0; j < n; ++j)
                        xc[j] = centroid[j] + cfg.rho * (xr[j] - centroid[j]);
                    Real fc = f(xc); ++n_evals;
                    if (fc <= fr) {
                        simplex[idx[n]] = xc;
                        fval[idx[n]] = fc;
                    } else {
                        // Shrink
                        for (Size i = 1; i <= n; ++i) {
                            for (Size j = 0; j < n; ++j)
                                simplex[idx[i]][j] = simplex[idx[0]][j] +
                                    cfg.sigma * (simplex[idx[i]][j] - simplex[idx[0]][j]);
                            fval[idx[i]] = f(simplex[idx[i]]); ++n_evals;
                        }
                    }
                } else {
                    // inside contraction
                    for (Size j = 0; j < n; ++j)
                        xc[j] = centroid[j] - cfg.rho * (centroid[j] - simplex[idx[n]][j]);
                    Real fc = f(xc); ++n_evals;
                    if (fc < fval[idx[n]]) {
                        simplex[idx[n]] = xc;
                        fval[idx[n]] = fc;
                    } else {
                        // Shrink
                        for (Size i = 1; i <= n; ++i) {
                            for (Size j = 0; j < n; ++j)
                                simplex[idx[i]][j] = simplex[idx[0]][j] +
                                    cfg.sigma * (simplex[idx[i]][j] - simplex[idx[0]][j]);
                            fval[idx[i]] = f(simplex[idx[i]]); ++n_evals;
                        }
                    }
                }
            }
        }

        // Return best vertex
        Size best = 0;
        for (Size i = 1; i <= n; ++i) if (fval[i] < fval[best]) best = i;
        result.x = simplex[best];
        result.fx = fval[best];
        result.n_function_evaluations = n_evals;
        if (result.message.empty()) {
            result.converged = false;
            result.message = "max iterations reached";
        }
        return result;
}

// ---------------------------------------------------------------------------
// Differential Evolution (rand/1/bin) - global optimizer
// Uses inline LCG (xorshift64*) instead of <random> to avoid MSVC ICE (C1001)
// ---------------------------------------------------------------------------
namespace detail {

// xorshift64* PRNG - fast, deterministic, no <random> dependency
inline uint64_t xorshift64(uint64_t& state) {
    uint64_t x = state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    state = x;
    return x * 0x2545F4914F6CDD1DULL;
}

// Uniform real in [0, 1) using xorshift64*
inline Real uniform01(uint64_t& state) {
    // Use top 53 bits for full double precision
    return static_cast<Real>(xorshift64(state) >> 11) * (1.0 / 9007199254740992.0);
}

}  // namespace detail

class DifferentialEvolution {
public:
    struct Config {
        Size population_size = 50;
        Size max_generations = 200;
        Real F = 0.8;       // mutation scale
        Real CR = 0.9;       // crossover probability
        Real tol = 1e-8;
        uint64_t seed = 42;
        // --- v1.1 Tikhonov 正则化 (Task 5) ---
        // 目标函数: f(x) + 0.5 * lambda_reg * ||x - params_prior||^2
        // 当 params_prior.size() != bounds.size() 或 lambda_reg <= 0 时不启用
        Real lambda_reg = 0.0;
        std::vector<Real> params_prior;
        // 早停: 当原始目标 f(x_best) < early_stop_rmse^2 * m / 2 时停止
        // (近似: DE 用 objective 而非 residual, 早停基于 best_fit)
        Real early_stop_rmse = 0.0;
    };

    static OptimizationResult minimize(
            ObjectiveFn f,
            const std::vector<Bounds>& bounds,
            const Config& cfg);
};

inline OptimizationResult DifferentialEvolution::minimize(
        ObjectiveFn f,
        const std::vector<Bounds>& bounds,
        const Config& cfg = DifferentialEvolution::Config{}) {
        OptimizationResult result;
        Size n = bounds.size();
        if (n == 0) {
            result.converged = true;
            result.message = "empty bounds";
            return result;
        }

        // v1.1: Tikhonov 正则化 — 包装 objective 函数
        bool use_reg = (cfg.lambda_reg > 0.0) &&
                       (cfg.params_prior.size() == n);
        ObjectiveFn f_wrapped = [&f, use_reg, &cfg, n](const std::vector<Real>& x) -> Real {
            Real base = f(x);
            if (!use_reg) return base;
            Real reg = 0.0;
            for (Size i = 0; i < n; ++i) {
                Real d = x[i] - cfg.params_prior[i];
                reg += d * d;
            }
            return base + 0.5 * cfg.lambda_reg * reg;
        };

        uint64_t rng_state = cfg.seed;
        if (rng_state == 0) rng_state = 1;  // xorshift needs nonzero state

        Size pop = cfg.population_size;
        std::vector<std::vector<Real>> population(pop, std::vector<Real>(n));
        std::vector<Real> fitness(pop);

        for (Size i = 0; i < pop; ++i) {
            for (Size j = 0; j < n; ++j) {
                population[i][j] = bounds[j].lower +
                    detail::uniform01(rng_state) * (bounds[j].upper - bounds[j].lower);
            }
            fitness[i] = f_wrapped(population[i]);
        }
        Size n_evals = pop;

        Real best_fit = fitness[0];
        Size best_idx = 0;
        for (Size i = 1; i < pop; ++i) {
            if (fitness[i] < best_fit) { best_fit = fitness[i]; best_idx = i; }
        }

        result.n_iterations = 0;
        for (Size gen = 0; gen < cfg.max_generations; ++gen) {
            result.n_iterations = gen + 1;
            bool improved = false;
            for (Size i = 0; i < pop; ++i) {
                // Pick r1 != r2 != r3 != i
                Size r1, r2, r3;
                do { r1 = static_cast<Size>(detail::uniform01(rng_state) * pop); } while (r1 == i);
                do { r2 = static_cast<Size>(detail::uniform01(rng_state) * pop); } while (r2 == i || r2 == r1);
                do { r3 = static_cast<Size>(detail::uniform01(rng_state) * pop); } while (r3 == i || r3 == r1 || r3 == r2);

                // Mutant vector v = x_r1 + F * (x_r2 - x_r3)
                std::vector<Real> mutant(n);
                for (Size j = 0; j < n; ++j) {
                    mutant[j] = population[r1][j] + cfg.F * (population[r2][j] - population[r3][j]);
                    // Clamp to bounds
                    mutant[j] = std::max(bounds[j].lower, std::min(bounds[j].upper, mutant[j]));
                }

                // Binomial crossover
                std::vector<Real> trial = population[i];
                Size j_rand = static_cast<Size>(detail::uniform01(rng_state) * n);
                for (Size j = 0; j < n; ++j) {
                    if (j == j_rand || detail::uniform01(rng_state) < cfg.CR) trial[j] = mutant[j];
                }

                Real f_trial = f_wrapped(trial); ++n_evals;
                if (f_trial <= fitness[i]) {
                    population[i] = trial;
                    fitness[i] = f_trial;
                    if (f_trial < best_fit) {
                        best_fit = f_trial;
                        best_idx = i;
                        improved = true;
                    }
                }
            }
            // v1.1: 早停检查 (基于原始 objective, 不含正则化项)
            if (cfg.early_stop_rmse > 0.0) {
                Real base_fit = f(population[best_idx]); ++n_evals;
                if (base_fit < cfg.early_stop_rmse * cfg.early_stop_rmse) {
                    result.converged = true;
                    result.message = "early_stop_rmse satisfied";
                    break;
                }
            }
            // Stagnation check
            if (!improved && gen > 20) {
                // Compute population spread
                Real spread = 0.0;
                for (Size i = 0; i < pop; ++i) {
                    for (Size j = 0; j < n; ++j) {
                        Real d = population[i][j] - population[best_idx][j];
                        spread += d * d;
                    }
                }
                spread = std::sqrt(spread / pop);
                if (spread < cfg.tol) {
                    result.converged = true;
                    result.message = "population converged";
                    break;
                }
            }
        }

        result.x = population[best_idx];
        // 返回原始 objective 值 (不含正则化), 便于跨校准器比较
        result.fx = use_reg ? f(result.x) : best_fit;
        result.n_function_evaluations = n_evals;
        if (result.message.empty()) {
            result.converged = false;
            result.message = "max generations reached";
        }
        return result;
}

// ---------------------------------------------------------------------------
// SLSQP: Sequential Least Squares Programming (Kraft 1988 NLPQL)
// Supports bounds + inequality constraints (c_i(x) >= 0) + equality constraints (c_i(x) = 0)
// Algorithm: QP subproblem (active-set) + Damped BFGS + L1 merit + Armijo line search
// ADR-018: implementation boundary decisions (7 items)
// ---------------------------------------------------------------------------
class SLSQP {
public:
    struct Config {
        Size max_iterations = 100;
        Real ftol = 1e-6;
        Real xtol = 1e-8;
        Real gtol = 1e-6;
        Real armijo_gamma = 0.1;
        Real armijo_beta = 0.5;
        Real epsilon = 1e-6;
        Size max_line_search = 20;
        GradientFn gradient;                  // optional, None => numerical diff
        ConstraintJacobianFn constraint_jacobian;  // optional, None => numerical diff
    };

    static OptimizationResult minimize(
            const ObjectiveFn& f,
            const std::vector<Real>& x0,
            const std::vector<Bounds>& bounds,
            const std::vector<ConstraintFn>& ineq_constraints,
            const std::vector<ConstraintFn>& eq_constraints,
            const Config& cfg);
};

inline OptimizationResult SLSQP::minimize(
        const ObjectiveFn& f,
        const std::vector<Real>& x0,
        const std::vector<Bounds>& bounds,
        const std::vector<ConstraintFn>& ineq_constraints = {},
        const std::vector<ConstraintFn>& eq_constraints = {},
        const SLSQP::Config& cfg = SLSQP::Config{}) {
    OptimizationResult result{};  // zero-init (RISK-014)
    Size n = x0.size();
    if (n == 0 || bounds.size() != n) {
        result.converged = false;
        result.message = "invalid input dimensions";
        return result;
    }

    // Project x0 onto bounds
    std::vector<Real> x = x0;
    for (Size i = 0; i < n; ++i) {
        x[i] = std::max(bounds[i].lower, std::min(bounds[i].upper, x[i]));
    }

    // BFGS Hessian approximation B = I
    std::vector<std::vector<Real>> B(n, std::vector<Real>(n, 0.0));
    for (Size i = 0; i < n; ++i) B[i][i] = 1.0;

    // Penalty parameters (L1 merit)
    std::vector<Real> mu_ineq(ineq_constraints.size(), 10.0);
    std::vector<Real> mu_eq(eq_constraints.size(), 10.0);

    Size n_evals = 0;
    Real fx = f(x); ++n_evals;
    result.n_iterations = 0;

    // Track previous Lagrangian gradient for BFGS update
    std::vector<Real> gradL_prev;
    bool has_gradL_prev = false;

    for (Size iter = 0; iter < cfg.max_iterations; ++iter) {
        result.n_iterations = iter + 1;

        // Evaluate objective gradient
        std::vector<Real> grad = cfg.gradient ? cfg.gradient(x)
                                              : detail::numerical_gradient(f, x, cfg.epsilon);
        if (!cfg.gradient) n_evals += 2 * n;

        // Evaluate inequality constraints and Jacobians
        // Flatten: each ConstraintFn returns a vector; we concatenate all.
        std::vector<std::vector<Real>> A_ineq;  // (m_ineq_total, n)
        std::vector<Real> c_ineq_flat;          // (m_ineq_total,)
        for (Size i = 0; i < ineq_constraints.size(); ++i) {
            auto cvals = ineq_constraints[i](x);
            auto J = detail::numerical_jacobian(ineq_constraints[i], x, cfg.epsilon);
            n_evals += 2 * n * cvals.size();
            for (Size k = 0; k < cvals.size(); ++k) {
                A_ineq.push_back(J[k]);
                c_ineq_flat.push_back(cvals[k]);
            }
        }

        // Evaluate equality constraints and Jacobians
        std::vector<std::vector<Real>> A_eq;
        std::vector<Real> c_eq_flat;
        for (Size i = 0; i < eq_constraints.size(); ++i) {
            auto cvals = eq_constraints[i](x);
            auto J = detail::numerical_jacobian(eq_constraints[i], x, cfg.epsilon);
            n_evals += 2 * n * cvals.size();
            for (Size k = 0; k < cvals.size(); ++k) {
                A_eq.push_back(J[k]);
                c_eq_flat.push_back(cvals[k]);
            }
        }

        // Add bound constraints as inequalities: d >= l - x,  -d >= -(u - x)
        // i.e.  I*d >= (l - x)  and  -I*d >= -(u - x)
        Size m_ineq_user = A_ineq.size();
        Size m_ineq_bounds = 2 * n;
        Size m_ineq_total = m_ineq_user + m_ineq_bounds;

        std::vector<std::vector<Real>> A_ineq_full(m_ineq_total, std::vector<Real>(n, 0.0));
        std::vector<Real> b_ineq_full(m_ineq_total, 0.0);
        // User inequalities: A_ineq[i]*d >= -c_ineq[i]  (linearization of c(x+d) >= 0)
        for (Size i = 0; i < m_ineq_user; ++i) {
            A_ineq_full[i] = A_ineq[i];
            b_ineq_full[i] = -c_ineq_flat[i];
        }
        // Bound lower: d_j >= l_j - x_j  =>  e_j * d >= (l - x)
        for (Size j = 0; j < n; ++j) {
            A_ineq_full[m_ineq_user + j][j] = 1.0;
            b_ineq_full[m_ineq_user + j] = bounds[j].lower - x[j];
        }
        // Bound upper: -d_j >= -(u_j - x_j)  =>  -e_j * d >= x_j - u_j
        for (Size j = 0; j < n; ++j) {
            A_ineq_full[m_ineq_user + n + j][j] = -1.0;
            b_ineq_full[m_ineq_user + n + j] = x[j] - bounds[j].upper;
        }

        // Equality: A_eq[i]*d = -c_eq[i]
        std::vector<Real> b_eq_full(c_eq_flat.size());
        for (Size i = 0; i < c_eq_flat.size(); ++i) b_eq_full[i] = -c_eq_flat[i];

        // Solve QP subproblem
        std::vector<Real> d;
        std::vector<Real> lam_ineq;
        std::vector<Real> lam_eq;
        bool qp_ok = detail::solve_qp_active_set(
            B, grad, A_ineq_full, b_ineq_full, A_eq, b_eq_full,
            d, lam_ineq, lam_eq);

        if (!qp_ok) {
            // QP failed: use steepest descent direction projected to bounds
            d.assign(n, 0.0);
            for (Size j = 0; j < n; ++j) d[j] = -grad[j];
        }

        // Convergence: step too small
        Real dxnorm = 0.0;
        for (Real v : d) dxnorm += v * v;
        dxnorm = std::sqrt(dxnorm);
        if (dxnorm < cfg.xtol) {
            result.converged = true;
            result.message = "xtol satisfied";
            break;
        }

        // Convergence: gradient (projected) small
        Real gnorm = 0.0;
        for (Real g : grad) gnorm += g * g;
        gnorm = std::sqrt(gnorm);
        // KKT-style: if step tiny and gradient small, converged
        if (gnorm < cfg.gtol && dxnorm < 10.0 * cfg.xtol) {
            result.converged = true;
            result.message = "gtol satisfied";
            break;
        }

        // Update penalty parameters: mu_i = max(mu_i, 1.5 * |lambda_i|) (Boggs-Tolle heuristic)
        for (Size i = 0; i < mu_ineq.size() && i < m_ineq_user; ++i) {
            // lam_ineq layout: [m_ineq_user | n lower-bounds | n upper-bounds]
            // But user inequality lambda for constraint function i corresponds to
            // flattened index i (since we flattened user constraints first).
            // Note: each ConstraintFn may return multiple values; mapping is complex.
            // Simplified: use average of relevant lambdas.
            Real lam_abs = std::abs(lam_ineq[i]);
            mu_ineq[i] = std::max(mu_ineq[i], 1.5 * lam_abs);
        }
        for (Size i = 0; i < mu_eq.size() && i < lam_eq.size(); ++i) {
            mu_eq[i] = std::max(mu_eq[i], 1.5 * std::abs(lam_eq[i]));
        }

        // Armijo line search on L1 merit
        Real merit_x = detail::l1_merit(f, x, ineq_constraints, eq_constraints, mu_ineq, mu_eq);
        // Directional derivative of merit at x along d
        Real D = 0.0;
        for (Size j = 0; j < n; ++j) D += grad[j] * d[j];
        // Subtract current constraint violation contribution (merit decreases when violation decreases)
        for (Size i = 0; i < ineq_constraints.size(); ++i) {
            auto cvals = ineq_constraints[i](x);
            for (Size k = 0; k < cvals.size(); ++k) {
                if (cvals[k] < 0.0) D -= mu_ineq[i] * (-cvals[k]);
            }
        }
        for (Size i = 0; i < eq_constraints.size(); ++i) {
            auto cvals = eq_constraints[i](x);
            for (Size k = 0; k < cvals.size(); ++k) {
                D -= mu_eq[i] * std::abs(cvals[k]);
            }
        }

        Real alpha = 1.0;
        bool line_search_ok = false;
        for (Size ls = 0; ls < cfg.max_line_search; ++ls) {
            std::vector<Real> x_new(n);
            for (Size j = 0; j < n; ++j) x_new[j] = x[j] + alpha * d[j];
            // Project to bounds
            for (Size j = 0; j < n; ++j) {
                x_new[j] = std::max(bounds[j].lower, std::min(bounds[j].upper, x_new[j]));
            }
            Real merit_new = detail::l1_merit(f, x_new, ineq_constraints, eq_constraints, mu_ineq, mu_eq);
            if (merit_new <= merit_x + cfg.armijo_gamma * alpha * D) {
                line_search_ok = true;
                // Accept step
                Real fx_new = f(x_new); ++n_evals;
                // Compute s and y for BFGS update
                std::vector<Real> s(n);
                for (Size j = 0; j < n; ++j) s[j] = x_new[j] - x[j];

                // Lagrangian gradient at x_new.
                // Our KKT convention: L = f + sum λ_eq * c_eq + sum λ_ineq * c_ineq
                // (because KKT system uses [B A^T; A 0][d;λ]=[-g;b], i.e. λ has opposite
                //  sign to the standard convention L = f - λ^T c).
                // So gradL = grad_f + sum λ_eq * grad_c_eq + sum λ_ineq * grad_c_ineq.
                std::vector<Real> gradL_new(n, 0.0);
                auto grad_new = cfg.gradient ? cfg.gradient(x_new)
                                             : detail::numerical_gradient(f, x_new, cfg.epsilon);
                if (!cfg.gradient) n_evals += 2 * n;
                for (Size j = 0; j < n; ++j) gradL_new[j] = grad_new[j];

                // Add equality constraint contributions
                Size eq_offset = 0;
                for (Size i = 0; i < eq_constraints.size(); ++i) {
                    auto J = detail::numerical_jacobian(eq_constraints[i], x_new, cfg.epsilon);
                    n_evals += 2 * n * J.size();
                    for (Size k = 0; k < J.size(); ++k) {
                        for (Size j = 0; j < n; ++j) {
                            gradL_new[j] += lam_eq[eq_offset + k] * J[k][j];
                        }
                    }
                    eq_offset += J.size();
                }
                // Add inequality constraint contributions (only active ones: λ < 0 in our convention)
                Size ineq_offset = 0;
                for (Size i = 0; i < ineq_constraints.size(); ++i) {
                    auto J = detail::numerical_jacobian(ineq_constraints[i], x_new, cfg.epsilon);
                    n_evals += 2 * n * J.size();
                    for (Size k = 0; k < J.size(); ++k) {
                        if (lam_ineq[ineq_offset + k] < -1e-10) {
                            for (Size j = 0; j < n; ++j) {
                                gradL_new[j] += lam_ineq[ineq_offset + k] * J[k][j];
                            }
                        }
                    }
                    ineq_offset += J.size();
                }

                // BFGS update
                if (has_gradL_prev) {
                    std::vector<Real> y(n);
                    for (Size j = 0; j < n; ++j) y[j] = gradL_new[j] - gradL_prev[j];
                    detail::bfgs_damped_update(B, s, y);
                }
                gradL_prev = gradL_new;
                has_gradL_prev = true;

                // Check convergence using merit function improvement (not raw fx),
                // because when starting from an infeasible point the first step may
                // reduce constraint violation without reducing f, and using fx would
                // falsely trigger ftol convergence.
                Real rel_merit_impr = (merit_x - merit_new) / std::max(std::abs(merit_x), 1e-30);
                x = x_new;
                fx = fx_new;

                if (rel_merit_impr < cfg.ftol) {
                    result.converged = true;
                    result.message = "ftol satisfied";
                }
                break;
            }
            alpha *= cfg.armijo_beta;
            if (alpha < 1e-15) break;
        }

        if (!line_search_ok) {
            // Line search failed: try tiny step or declare convergence
            result.converged = false;
            result.message = "line search failed";
            break;
        }

        if (result.converged) break;
    }

    result.x = x;
    result.fx = fx;
    result.n_function_evaluations = n_evals;
    if (result.message.empty()) {
        result.converged = false;
        result.message = "max iterations reached";
    }
    return result;
}

}  // inline namespace v1
}  // namespace cpphub
