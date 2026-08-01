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

}  // inline namespace v1
}  // namespace cpphub
