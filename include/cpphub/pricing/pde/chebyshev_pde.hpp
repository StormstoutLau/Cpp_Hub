// Chebyshev collocation spectral method for 1D BSM PDE
//
// Discretizes the spatial domain using Chebyshev-Gauss-Lobatto points and
// approximates derivatives via the Chebyshev differentiation matrix. For
// smooth solutions this yields spectral (exponential) convergence, far
// superior to finite-difference methods which converge only polynomially.
//
// BSM PDE in log-price y = ln(S/K):
//   dV/dtau = 0.5*sigma^2 * d2V/dy2 + (r-q-0.5*sigma^2) * dV/dy - r*V
// with V(y, 0) = payoff(K*exp(y)) and Dirichlet far-field BCs.
//
// Time discretization: theta-scheme (theta=0.5 => Crank-Nicolson, second-order).
// The LHS matrix M2 = I - theta*dt*A is LU-factored once and reused per step.
//
// Reference: Trefethen (2000) "Spectral Methods in MATLAB", Chapter 6.
//            Boyd (2001) "Chebyshev and Fourier Spectral Methods".
#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {

struct ChebyshevPDEConfig {
    Size n_points = 64;       // Number of Chebyshev intervals; total nodes = n_points + 1
    Size n_time = 200;        // Number of time steps
    Real x_range = 5.0;       // Half-width of log-S domain (y in [y0-L, y0+L])
    bool is_call = true;      // true = call, false = put
    Real theta = 0.5;         // 0.5 = Crank-Nicolson, 1.0 = implicit Euler
    // Rannacher smoothing: first n_rannacher_warmup steps use implicit Euler
    // (theta=1.0) to dampen Gibbs oscillations from the non-smooth payoff
    // kink at S=K, then switch to theta-scheme (typically CN) for accuracy.
    // Set to 0 to disable. Default 4 is standard for BSM.
    Size n_rannacher_warmup = 4;
};

class ChebyshevPDEEngine {
public:
    explicit ChebyshevPDEEngine(ChebyshevPDEConfig config)
        : config_(config) {
        if (config_.n_points < 4) {
            throw std::invalid_argument("ChebyshevPDEEngine: n_points must be >= 4");
        }
        if (config_.x_range <= 0.0) {
            throw std::invalid_argument("ChebyshevPDEEngine: x_range must be positive");
        }
    }

    // Price a European vanilla option under BSM via Chebyshev collocation.
    // The spot S0 maps to the center of the domain (x = 0 in [-1, 1]).
    Real price(Real S0, Real K, Real T, Real r, Real q, Real sigma) const {
        if (T <= 0.0) {
            Real S = S0;
            return config_.is_call ? std::max(S - K, 0.0) : std::max(K - S, 0.0);
        }
        if (sigma <= 0.0) {
            // Degenerate: deterministic forward
            Real fwd = S0 * std::exp((r - q) * T);
            Real disc = std::exp(-r * T);
            return config_.is_call
                ? std::max(fwd - K, 0.0) * disc
                : std::max(K - fwd, 0.0) * disc;
        }

        const Size N = config_.n_points;
        const Size n = N + 1;  // total nodes

        const Real y0 = std::log(S0 / K);
        const Real L  = config_.x_range;

        // --- Chebyshev-Gauss-Lobatto points on [-1, 1] ---
        // x[k] = cos(pi * k / N), k = 0, ..., N
        // Note: x[0] = +1 (=> y_max => S_max), x[N] = -1 (=> y_min => S_min)
        std::vector<Real> x(n);
        const Real pi_val = pi();
        for (Size k = 0; k < n; ++k) {
            x[k] = std::cos(pi_val * static_cast<Real>(k) / static_cast<Real>(N));
        }

        // Map to y-domain: y[i] = y0 + L * x[i]
        std::vector<Real> y(n);
        for (Size i = 0; i < n; ++i) {
            y[i] = y0 + L * x[i];
        }

        // --- Chebyshev differentiation matrix D on [-1, 1] ---
        // Trefethen (2000) Program 17.
        auto D = build_diff_matrix(N, x);

        // D2 = D * D  (second-derivative matrix on [-1, 1])
        auto D2 = mat_mul(D, D);

        // --- BSM spatial operator in y-domain ---
        // dy/dx = L  =>  d/dy = (1/L) d/dx,  d2/dy2 = (1/L^2) d2/dx2
        // A = (sigma^2 / (2 L^2)) D2 + ((r - q - 0.5 sigma^2) / L) D - r I
        const Real sigma2 = sigma * sigma;
        const Real drift  = r - q - 0.5 * sigma2;
        const Real coef_D2 = sigma2 / (2.0 * L * L);
        const Real coef_D  = drift / L;

        std::vector<std::vector<Real>> A(n, std::vector<Real>(n, 0.0));
        for (Size i = 0; i < n; ++i) {
            for (Size j = 0; j < n; ++j) {
                A[i][j] = coef_D2 * D2[i][j] + coef_D * D[i][j];
            }
            A[i][i] -= r;
        }

        // --- Time-stepping matrices (theta-scheme) ---
        // (I - theta*dt*A) V^{n+1} = (I + (1-theta)*dt*A) V^n
        const Real dt = T / static_cast<Real>(config_.n_time);
        const Real th = config_.theta;

        // Rannacher smoothing: first n_warmup steps use implicit Euler (theta=1.0)
        // to dampen Gibbs oscillations from the payoff kink at S=K.
        // Implicit Euler amplification factor 1/(1-dt*lambda) strongly damps
        // high-frequency modes (large |lambda|), while CN factor
        // (1+0.5*dt*lambda)/(1-0.5*dt*lambda) only flips their sign without
        // significant decay — causing persistent oscillations.
        Size n_warmup = std::min(config_.n_rannacher_warmup, config_.n_time);

        // Build M2_main = I - theta*dt*A, M1_main = I + (1-theta)*dt*A
        auto build_matrices = [&](Real theta_val) {
            std::vector<std::vector<Real>> M1m(n, std::vector<Real>(n, 0.0));
            std::vector<std::vector<Real>> M2m(n, std::vector<Real>(n, 0.0));
            for (Size i = 0; i < n; ++i) {
                for (Size j = 0; j < n; ++j) {
                    M1m[i][j] = (1.0 - theta_val) * dt * A[i][j];
                    M2m[i][j] = -theta_val * dt * A[i][j];
                }
                M1m[i][i] += 1.0;
                M2m[i][i] += 1.0;
            }
            // Replace boundary rows of M2 with unit rows (Dirichlet BCs)
            for (Size j = 0; j < n; ++j) {
                M2m[0][j]     = (j == 0)     ? 1.0 : 0.0;
                M2m[n - 1][j] = (j == n - 1) ? 1.0 : 0.0;
            }
            return std::make_pair(std::move(M1m), std::move(M2m));
        };

        // Warmup matrices (theta=1.0 => M1=I, M2=I-dt*A)
        auto [M1_warm, M2_warm] = (n_warmup > 0)
            ? build_matrices(1.0)
            : std::make_pair(std::vector<std::vector<Real>>{},
                             std::vector<std::vector<Real>>{});
        auto lu_warm = (n_warmup > 0) ? lu_decompose(M2_warm) : LUFactor{};

        // Main matrices (theta=config_.theta, typically 0.5 = CN)
        auto [M1_main, M2_main] = build_matrices(th);
        auto lu_main = lu_decompose(std::move(M2_main));

        // --- Initial condition: payoff at tau = 0 ---
        std::vector<Real> V(n);
        for (Size i = 0; i < n; ++i) {
            Real S = K * std::exp(y[i]);
            V[i] = config_.is_call ? std::max(S - K, 0.0)
                                   : std::max(K - S, 0.0);
        }

        // Helper: set boundary entries of b given tau
        auto set_bc = [&](std::vector<Real>& bvec, Real tau) {
            Real disc_r = std::exp(-r * tau);
            if (config_.is_call) {
                Real S_max = K * std::exp(y[0]);
                bvec[0]     = S_max - K * disc_r;
                bvec[n - 1] = 0.0;
            } else {
                bvec[0]     = 0.0;
                bvec[n - 1] = K * disc_r;
            }
        };

        // Helper: mat-vec b = M * V (full n x n)
        auto mat_vec = [&](const std::vector<std::vector<Real>>& M,
                           const std::vector<Real>& v,
                           std::vector<Real>& bvec) {
            for (Size i = 0; i < n; ++i) {
                Real s = 0.0;
                const auto& row = M[i];
                for (Size j = 0; j < n; ++j) {
                    s += row[j] * v[j];
                }
                bvec[i] = s;
            }
        };

        // --- Time stepping ---
        std::vector<Real> b(n, 0.0);
        for (Size step = 0; step < config_.n_time; ++step) {
            Real tau = static_cast<Real>(step + 1) * dt;

            if (step < n_warmup) {
                // Warmup: implicit Euler. M1 = I, so b = V.
                b = V;
            } else {
                // Main: theta-scheme. b = M1_main * V.
                mat_vec(M1_main, V, b);
            }

            // Apply Dirichlet boundary conditions
            set_bc(b, tau);

            // Solve
            const LUFactor& lu = (step < n_warmup) ? lu_warm : lu_main;
            V = lu_solve(lu, b);
        }

        // --- Interpolate at S0 (x = 0) via barycentric formula ---
        return barycentric_interp(x, V, 0.0);
    }

    const ChebyshevPDEConfig& config() const { return config_; }

private:
    ChebyshevPDEConfig config_;

    static Real pi() {
        // std::acos(-1.0) is portable across MSVC / GCC / Clang
        return std::acos(-1.0);
    }

    // Build Chebyshev differentiation matrix D on [-1, 1].
    // Trefethen (2000) "Spectral Methods in MATLAB", Program 17.
    // x[k] = cos(pi * k / N) for k = 0..N.
    // c_k = 2*(-1)^k for k in {0, N}, 1*(-1)^k otherwise.
    // D[i][j] = c_i / c_j / (x_i - x_j)  for i != j
    // D[i][i] = -sum_{j != i} D[i][j]
    static std::vector<std::vector<Real>> build_diff_matrix(
            Size N, const std::vector<Real>& x) {
        Size n = N + 1;
        std::vector<Real> c(n);
        for (Size i = 0; i < n; ++i) {
            Real base = (i == 0 || i == n - 1) ? 2.0 : 1.0;
            c[i] = base * ((i % 2 == 0) ? 1.0 : -1.0);
        }

        std::vector<std::vector<Real>> D(n, std::vector<Real>(n, 0.0));
        for (Size i = 0; i < n; ++i) {
            for (Size j = 0; j < n; ++j) {
                if (i != j) {
                    D[i][j] = c[i] / c[j] / (x[i] - x[j]);
                }
            }
        }
        // Diagonal: negative sum of off-diagonal entries (row-sum = 0
        // for differentiation of a constant).
        for (Size i = 0; i < n; ++i) {
            Real s = 0.0;
            for (Size j = 0; j < n; ++j) {
                if (j != i) s += D[i][j];
            }
            D[i][i] = -s;
        }
        return D;
    }

    // Dense matrix-matrix product C = A * B (both n x n).
    static std::vector<std::vector<Real>> mat_mul(
            const std::vector<std::vector<Real>>& A,
            const std::vector<std::vector<Real>>& B) {
        Size n = A.size();
        std::vector<std::vector<Real>> C(n, std::vector<Real>(n, 0.0));
        for (Size i = 0; i < n; ++i) {
            for (Size k = 0; k < n; ++k) {
                Real aik = A[i][k];
                if (aik == 0.0) continue;
                for (Size j = 0; j < n; ++j) {
                    C[i][j] += aik * B[k][j];
                }
            }
        }
        return C;
    }

    // LU factorization with partial pivoting (in-place on a copy of A).
    // Stores both L (unit lower, multipliers in strictly-lower part of U)
    // and U (upper triangular) in a single matrix for cache efficiency.
    struct LUFactor {
        std::vector<std::vector<Real>> U;  // combined L\U storage
        std::vector<Size> piv;             // row permutation
    };

    static LUFactor lu_decompose(std::vector<std::vector<Real>> A) {
        Size n = A.size();
        LUFactor lu;
        lu.U = std::move(A);
        lu.piv.resize(n);
        for (Size i = 0; i < n; ++i) lu.piv[i] = i;

        for (Size k = 0; k < n; ++k) {
            // Partial pivoting: find row p with max |U[p][k]| in column k
            Size p = k;
            Real max_val = std::abs(lu.U[k][k]);
            for (Size i = k + 1; i < n; ++i) {
                Real v = std::abs(lu.U[i][k]);
                if (v > max_val) { max_val = v; p = i; }
            }
            if (p != k) {
                std::swap(lu.U[p], lu.U[k]);
                std::swap(lu.piv[p], lu.piv[k]);
            }

            Real pivot = lu.U[k][k];
            if (std::abs(pivot) < 1e-300) {
                throw std::runtime_error("ChebyshevPDE: LU singular matrix");
            }

            // Eliminate below pivot; store multipliers in lower part of U.
            for (Size i = k + 1; i < n; ++i) {
                Real m = lu.U[i][k] / pivot;
                lu.U[i][k] = m;  // multiplier (L[i][k])
                for (Size j = k + 1; j < n; ++j) {
                    lu.U[i][j] -= m * lu.U[k][j];
                }
            }
        }
        return lu;
    }

    // Solve A x = b given LU factorization. O(n^2).
    static std::vector<Real> lu_solve(const LUFactor& lu,
                                       const std::vector<Real>& b) {
        Size n = b.size();
        // Apply permutation: y = Pb
        std::vector<Real> y(n);
        for (Size i = 0; i < n; ++i) y[i] = b[lu.piv[i]];

        // Forward substitution: L y' = y  (L unit-diagonal, lower part of U)
        for (Size i = 1; i < n; ++i) {
            Real s = y[i];
            for (Size j = 0; j < i; ++j) {
                s -= lu.U[i][j] * y[j];
            }
            y[i] = s;
        }
        // Backward substitution: U x = y'
        std::vector<Real> x(n);
        for (Size i = n; i-- > 0; ) {
            Real s = y[i];
            for (Size j = i + 1; j < n; ++j) {
                s -= lu.U[i][j] * x[j];
            }
            x[i] = s / lu.U[i][i];
        }
        return x;
    }

    // Barycentric interpolation at point x0 for data (x_i, y_i) on
    // Chebyshev-Gauss-Lobatto nodes. Exact for polynomials of degree <= N.
    static Real barycentric_interp(const std::vector<Real>& x,
                                    const std::vector<Real>& y,
                                    Real x0) {
        Size n = x.size();

        // If x0 coincides with a node, return the node value directly.
        for (Size i = 0; i < n; ++i) {
            if (std::abs(x0 - x[i]) < 1e-14) return y[i];
        }

        // Barycentric weights for CGL points:
        //   w_k = (-1)^k           for 1 <= k <= N-1
        //   w_k = (-1)^k / 2        for k = 0 or k = N
        std::vector<Real> w(n);
        for (Size i = 0; i < n; ++i) {
            w[i] = (i % 2 == 0) ? 1.0 : -1.0;
            if (i == 0 || i == n - 1) w[i] *= 0.5;
        }

        Real num = 0.0, den = 0.0;
        for (Size i = 0; i < n; ++i) {
            Real diff = x0 - x[i];
            Real t = w[i] / diff;
            num += t * y[i];
            den += t;
        }
        return num / den;
    }
};

}  // namespace v1
}  // namespace cpphub
