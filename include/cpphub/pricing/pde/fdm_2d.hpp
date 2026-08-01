// ADI (Alternating Direction Implicit) schemes for 2D PDEs
// Reference: in 't Hout & Foulon (2010) "ADI finite difference discretization
//            for the Heston PDE", Craig & Sneyd (1988), Hundsdorfer & Verwer (2003)
//
// Heston PDE (log-S transform, x = ln(S/K)):
//   dV/dt = 0.5 v d2V/dx2 + (r-q-0.5v) dV/dx
//         + 0.5 xi^2 v d2V/dv2 + kappa(theta-v) dV/dv
//         + rho xi v d2V/(dx dv) - r V
//
// Operator splitting: L = L_x + L_v + L_xv + L_0
//   L_x   = 0.5 v d2/dx2 + (r-q-0.5v) d/dx           (x-direction convection-diffusion)
//   L_v   = 0.5 xi^2 v d2/dv2 + kappa(theta-v) d/dv  (v-direction convection-diffusion)
//   L_xv  = rho xi v d2/(dx dv)                       (mixed derivative, treated explicitly)
//   L_0   = -r                                         (discount)
//
// ADI schemes (in 't Hout & Foulon 2010, with A0=L_x, A1=L_v, A2=L_xv, A=L_x+L_v+L_xv+L_0):
//   Y0 = V_n + dt*A*V_n                                (explicit predictor)
//   (I - theta*dt*A0) Y1 = Y0 - theta*dt*A0*V_n        (implicit x)
//   (I - theta*dt*A1) Y2 = Y1 - theta*dt*A1*V_n        (implicit v)
//   CS:  Yhat = Y2 + theta*dt*A2*V_n + (theta-0.5)*dt*A*Y0
//   HV:  Yhat = Y2                    + (theta-0.5)*dt*A*Y0
//   ytilde = Yhat - Y2 + Y1
//   (I - theta*dt*A0) ytilde1 = ytilde - theta*dt*A0*Yhat
//   (I - theta*dt*A1) ytilde2 = ytilde1 - theta*dt*A1*Yhat
//   CS:  V_{n+1} = ytilde2 + theta*dt*A2*Yhat + (theta-0.5)*dt*A*Yhat
//   HV:  V_{n+1} = ytilde2                    + (theta-0.5)*dt*A*Yhat
//
// Note: CS explicitly adds the mixed-derivative term A2=L_xv in both correction stages;
//       HV folds A2 into the full operator A (handled via L_full application).
// Both schemes are second-order accurate in time when theta in [0.5, 1].
#pragma once
#include <vector>
#include <cmath>
#include "cpphub/core/types.hpp"
#include "cpphub/pricing/pde/thomas_solver.hpp"

namespace cpphub {
inline namespace v1 {

// Heston model parameters for PDE
struct HestonPDEParams {
    Real kappa;    // mean reversion speed
    Real theta;    // long-term variance
    Real xi;       // vol of vol
    Real rho;      // correlation
    Real r;        // risk-free rate
    Real q;        // dividend yield
    Real T;        // maturity
    Real K;        // strike
    Real S0;       // spot
    Real v0;       // initial variance
};

// 2D grid (x = log-S, v = variance) supporting uniform and sinh-transformed grids.
// The sinh transformation concentrates grid points near:
//   x = x_center (ATM, where payoff has kink)
//   v = v_min    (where Feller boundary layer causes rapid solution variation)
// When alpha=0, the sinh grid degenerates to uniform (L'Hôpital).
class FDMGrid2D {
public:
    // Uniform grid constructor (backward compatible)
    FDMGrid2D(Size n_x, Size n_v, Real x_min, Real x_max,
              Real v_min, Real v_max)
        : n_x_(n_x), n_v_(n_v),
          x_min_(x_min), x_max_(x_max), v_min_(v_min), v_max_(v_max),
          dx_((x_max - x_min) / static_cast<Real>(n_x - 1)),
          dv_((v_max - v_min) / static_cast<Real>(n_v - 1)),
          is_uniform_(true) {
        x_points_.resize(n_x);
        v_points_.resize(n_v);
        for (Size i = 0; i < n_x; ++i) x_points_[i] = x_min + static_cast<Real>(i) * dx_;
        for (Size j = 0; j < n_v; ++j) v_points_[j] = v_min + static_cast<Real>(j) * dv_;
    }

    // Non-uniform grid constructor (sinh transformation)
    // alpha_x > 0: concentrates x points near x_center = (x_min + x_max)/2
    // alpha_v > 0: concentrates v points near v_min
    // alpha = 0: degenerates to uniform (via L'Hôpital)
    FDMGrid2D(Size n_x, Size n_v, Real x_min, Real x_max,
              Real v_min, Real v_max, Real alpha_x, Real alpha_v)
        : n_x_(n_x), n_v_(n_v),
          x_min_(x_min), x_max_(x_max), v_min_(v_min), v_max_(v_max),
          is_uniform_(false) {
        x_points_.resize(n_x);
        v_points_.resize(n_v);

        // x: symmetric sinh around center
        Real x_center = (x_min + x_max) / 2.0;
        Real x_range = (x_max - x_min) / 2.0;
        if (std::abs(alpha_x) < 1e-14) {
            for (Size i = 0; i < n_x; ++i)
                x_points_[i] = x_min + static_cast<Real>(i) * (x_max - x_min) / static_cast<Real>(n_x - 1);
        } else {
            Real inv_sinh_ax = 1.0 / std::sinh(alpha_x);
            for (Size i = 0; i < n_x; ++i) {
                Real xi = -1.0 + 2.0 * static_cast<Real>(i) / static_cast<Real>(n_x - 1);
                x_points_[i] = x_center + x_range * std::sinh(alpha_x * xi) * inv_sinh_ax;
            }
        }

        // v: asymmetric sinh concentrating near v_min (j=0)
        Real v_range = v_max - v_min;
        if (std::abs(alpha_v) < 1e-14) {
            for (Size j = 0; j < n_v; ++j)
                v_points_[j] = v_min + static_cast<Real>(j) * v_range / static_cast<Real>(n_v - 1);
        } else {
            Real inv_sinh_av = 1.0 / std::sinh(alpha_v);
            for (Size j = 0; j < n_v; ++j) {
                Real eta = static_cast<Real>(j) / static_cast<Real>(n_v - 1);
                v_points_[j] = v_min + v_range * std::sinh(alpha_v * eta) * inv_sinh_av;
            }
        }

        // Average spacing (for diagnostics like CFL check)
        dx_ = (x_max - x_min) / static_cast<Real>(n_x - 1);
        dv_ = (v_max - v_min) / static_cast<Real>(n_v - 1);
    }

    Size n_x() const noexcept { return n_x_; }
    Size n_v() const noexcept { return n_v_; }
    Size size() const noexcept { return n_x_ * n_v_; }

    Real x(Size i) const noexcept { return x_points_[i]; }
    Real v(Size j) const noexcept { return v_points_[j]; }

    Real x_min() const noexcept { return x_min_; }
    Real x_max() const noexcept { return x_max_; }
    Real v_min() const noexcept { return v_min_; }
    Real v_max() const noexcept { return v_max_; }

    // Average spacing (exact for uniform, approximate for non-uniform)
    // Used for CFL diagnostics only. Use dx_plus/dx_minus for FD coefficients.
    Real dx() const noexcept { return dx_; }
    Real dv() const noexcept { return dv_; }

    // Non-uniform spacing: forward and backward differences at point i
    Real dx_plus(Size i) const noexcept { return x_points_[i + 1] - x_points_[i]; }
    Real dx_minus(Size i) const noexcept { return x_points_[i] - x_points_[i - 1]; }
    Real dv_plus(Size j) const noexcept { return v_points_[j + 1] - v_points_[j]; }
    Real dv_minus(Size j) const noexcept { return v_points_[j] - v_points_[j - 1]; }

    bool is_uniform() const noexcept { return is_uniform_; }

    // Index: V[i + j * n_x] (v varies slowly, x varies fast)
    Size idx(Size i, Size j) const noexcept { return i + j * n_x_; }

private:
    Size n_x_, n_v_;
    Real x_min_, x_max_, v_min_, v_max_;
    Real dx_, dv_;  // average spacing
    bool is_uniform_;
    std::vector<Real> x_points_, v_points_;
};

// 2D Heston operator coefficients for ADI splitting.
// Each sub-operator is discretized via central differences on the grid
// (supports both uniform and non-uniform spacing).
// Storage: all arrays indexed [j][i] (v index outer, x index inner).
struct HestonOperatorCoeffs {
    std::vector<std::vector<Real>> a_x, b_x, c_x;  // L_x tridiagonal (per v_j row)
    std::vector<std::vector<Real>> a_v, b_v, c_v;  // L_v tridiagonal (per x_i column)
    std::vector<std::vector<Real>> m_xv;           // L_xv coefficient rho*xi*v_j
    Real L0;                                       // L_0 = -r (constant discount)
};

inline HestonOperatorCoeffs compute_heston_coeffs(const FDMGrid2D& grid,
                                                    const HestonPDEParams& p) {
    Size nx = grid.n_x(), nv = grid.n_v();
    Real rho_xi = p.rho * p.xi;

    HestonOperatorCoeffs c;
    c.a_x.assign(nv, std::vector<Real>(nx, 0.0));
    c.b_x.assign(nv, std::vector<Real>(nx, 0.0));
    c.c_x.assign(nv, std::vector<Real>(nx, 0.0));
    c.a_v.assign(nv, std::vector<Real>(nx, 0.0));
    c.b_v.assign(nv, std::vector<Real>(nx, 0.0));
    c.c_v.assign(nv, std::vector<Real>(nx, 0.0));
    c.m_xv.assign(nv, std::vector<Real>(nx, 0.0));
    c.L0 = -p.r;

    for (Size j = 0; j < nv; ++j) {
        Real v_j = grid.v(j);
        Real drift_x = p.r - p.q - 0.5 * v_j;
        Real drift_v = p.kappa * (p.theta - v_j);
        Real xi2_v = p.xi * p.xi * v_j;
        for (Size i = 0; i < nx; ++i) {
            if (i > 0 && i < nx - 1) {
                // Non-uniform central differences for L_x = 0.5*v*d²/dx² + drift_x*d/dx
                // With hp = dx_plus(i), hm = dx_minus(i):
                //   a = [v - drift_x*hp] / (hm*(hp+hm))
                //   b = [-v + drift_x*(hp-hm)] / (hp*hm)
                //   c = [v + drift_x*hm] / (hp*(hp+hm))
                // Reduces to uniform formula when hp=hm=h:
                //   a = v/(2h²) - drift_x/(2h), b = -v/h², c = v/(2h²) + drift_x/(2h)
                Real hp_x = grid.dx_plus(i);
                Real hm_x = grid.dx_minus(i);
                Real sum_x = hp_x + hm_x;
                Real prod_x = hp_x * hm_x;
                c.a_x[j][i] = (v_j - drift_x * hp_x) / (hm_x * sum_x);
                c.b_x[j][i] = (-v_j + drift_x * (hp_x - hm_x)) / prod_x;
                c.c_x[j][i] = (v_j + drift_x * hm_x) / (hp_x * sum_x);
            }
            if (j > 0 && j < nv - 1) {
                // Non-uniform central differences for L_v = 0.5*xi²*v*d²/dv² + drift_v*d/dv
                Real hp_v = grid.dv_plus(j);
                Real hm_v = grid.dv_minus(j);
                Real sum_v = hp_v + hm_v;
                Real prod_v = hp_v * hm_v;
                c.a_v[j][i] = (xi2_v - drift_v * hp_v) / (hm_v * sum_v);
                c.b_v[j][i] = (-xi2_v + drift_v * (hp_v - hm_v)) / prod_v;
                c.c_v[j][i] = (xi2_v + drift_v * hm_v) / (hp_v * sum_v);
            }
            c.m_xv[j][i] = rho_xi * v_j;
        }
    }
    return c;
}

// Apply mixed derivative operator L_xv (central 5-point stencil):
// Uniform: (L_xv V)[i,j] = rho*xi*v_j/(4*dx*dv) * (V[i+1,j+1] - V[i-1,j+1] - V[i+1,j-1] + V[i-1,j-1])
// Non-uniform: denominator = (dx_plus+dx_minus)*(dv_plus+dv_minus)
// (reduces to 4*dx*dv when uniform)
inline void apply_mixed_derivative(const std::vector<Real>& V,
                                    const FDMGrid2D& grid,
                                    const HestonOperatorCoeffs& c,
                                    std::vector<Real>& out) {
    Size nx = grid.n_x(), nv = grid.n_v();
    out.assign(V.size(), 0.0);
    for (Size j = 1; j < nv - 1; ++j) {
        Real dvp = grid.dv_plus(j);
        Real dvm = grid.dv_minus(j);
        Real denom_v = dvp + dvm;
        for (Size i = 1; i < nx - 1; ++i) {
            Real dxp = grid.dx_plus(i);
            Real dxm = grid.dx_minus(i);
            Real denom_x = dxp + dxm;
            Real coef = c.m_xv[j][i] / (denom_x * denom_v);
            out[grid.idx(i, j)] = coef * (V[grid.idx(i + 1, j + 1)]
                                        - V[grid.idx(i - 1, j + 1)]
                                        - V[grid.idx(i + 1, j - 1)]
                                        + V[grid.idx(i - 1, j - 1)]);
        }
    }
}

// Apply full operator L = L_x + L_v + L_xv + L_0 to V
// Boundary points (i=0, i=nx-1, j=0, j=nv-1) return 0 so that the predictor
// Y0 = V + dt*L*V keeps V unchanged at boundaries. The actual boundary evolution
// is handled externally by apply_boundary_conditions() in PDEEngine2D.
inline void apply_heston_operator(const std::vector<Real>& V,
                                   const FDMGrid2D& grid,
                                   const HestonOperatorCoeffs& c,
                                   std::vector<Real>& out) {
    Size nx = grid.n_x(), nv = grid.n_v();
    out.assign(V.size(), 0.0);
    for (Size j = 1; j < nv - 1; ++j) {
        for (Size i = 1; i < nx - 1; ++i) {
            Size k = grid.idx(i, j);
            Real val = c.L0 * V[k];
            val += c.a_x[j][i] * V[grid.idx(i - 1, j)]
                 + c.b_x[j][i] * V[k]
                 + c.c_x[j][i] * V[grid.idx(i + 1, j)];
            val += c.a_v[j][i] * V[grid.idx(i, j - 1)]
                 + c.b_v[j][i] * V[k]
                 + c.c_v[j][i] * V[grid.idx(i, j + 1)];
            out[k] = val;
        }
    }
    std::vector<Real> mixed;
    apply_mixed_derivative(V, grid, c, mixed);
    for (Size k = 0; k < out.size(); ++k) out[k] += mixed[k];
}

// Apply only L_x to V (interior points only; boundary remains 0)
// Loop range j=1..nv-2, i=1..nx-2 ensures v-boundary rows and x-boundary cols
// are NOT computed, keeping LxVn=0 at all boundaries. This is critical for ADI
// stability: if LxVn≠0 at v-boundary, Y1=Y0-theta*dt*LxVn corrupts the boundary,
// which propagates into interior via subsequent solve_v_direction.
inline void apply_x_operator(const std::vector<Real>& V,
                              const FDMGrid2D& grid,
                              const HestonOperatorCoeffs& c,
                              std::vector<Real>& out) {
    Size nx = grid.n_x(), nv = grid.n_v();
    out.assign(V.size(), 0.0);
    for (Size j = 1; j < nv - 1; ++j) {
        for (Size i = 1; i < nx - 1; ++i) {
            Size k = grid.idx(i, j);
            out[k] = c.a_x[j][i] * V[grid.idx(i - 1, j)]
                   + c.b_x[j][i] * V[k]
                   + c.c_x[j][i] * V[grid.idx(i + 1, j)];
        }
    }
}

// Apply only L_v to V (interior points only; boundary remains 0)
// Loop range j=1..nv-2, i=1..nx-2 ensures x-boundary cols and v-boundary rows
// are NOT computed, keeping LvVn=0 at all boundaries.
inline void apply_v_operator(const std::vector<Real>& V,
                              const FDMGrid2D& grid,
                              const HestonOperatorCoeffs& c,
                              std::vector<Real>& out) {
    Size nx = grid.n_x(), nv = grid.n_v();
    out.assign(V.size(), 0.0);
    for (Size j = 1; j < nv - 1; ++j) {
        for (Size i = 1; i < nx - 1; ++i) {
            Size k = grid.idx(i, j);
            out[k] = c.a_v[j][i] * V[grid.idx(i, j - 1)]
                   + c.b_v[j][i] * V[k]
                   + c.c_v[j][i] * V[grid.idx(i, j + 1)];
        }
    }
}

// Solve (I - theta_dt*L_x) Y = RHS along x-direction (Thomas algorithm per v_j row).
// Boundary rows (i=0, i=nx-1) keep identity (V unchanged).
inline void solve_x_direction(std::vector<Real>& V,
                               const FDMGrid2D& grid,
                               const HestonOperatorCoeffs& c,
                               Real theta_dt) {
    Size nx = grid.n_x(), nv = grid.n_v();
    std::vector<Real> a(nx, 0.0), b(nx, 0.0), cc(nx, 0.0), d(nx, 0.0);
    for (Size j = 0; j < nv; ++j) {
        a[0] = 0.0; b[0] = 1.0; cc[0] = 0.0; d[0] = V[grid.idx(0, j)];
        a[nx - 1] = 0.0; b[nx - 1] = 1.0; cc[nx - 1] = 0.0; d[nx - 1] = V[grid.idx(nx - 1, j)];
        for (Size i = 1; i < nx - 1; ++i) {
            a[i] = -theta_dt * c.a_x[j][i];
            b[i] = 1.0 - theta_dt * c.b_x[j][i];
            cc[i] = -theta_dt * c.c_x[j][i];
            d[i] = V[grid.idx(i, j)];
        }
        auto x = thomas_solve(a, b, cc, d);
        for (Size i = 0; i < nx; ++i) V[grid.idx(i, j)] = x[i];
    }
}

// Solve (I - theta_dt*L_v) Y = RHS along v-direction (Thomas algorithm per x_i column).
inline void solve_v_direction(std::vector<Real>& V,
                               const FDMGrid2D& grid,
                               const HestonOperatorCoeffs& c,
                               Real theta_dt) {
    Size nx = grid.n_x(), nv = grid.n_v();
    std::vector<Real> a(nv, 0.0), b(nv, 0.0), cc(nv, 0.0), d(nv, 0.0);
    for (Size i = 0; i < nx; ++i) {
        a[0] = 0.0; b[0] = 1.0; cc[0] = 0.0; d[0] = V[grid.idx(i, 0)];
        a[nv - 1] = 0.0; b[nv - 1] = 1.0; cc[nv - 1] = 0.0; d[nv - 1] = V[grid.idx(i, nv - 1)];
        for (Size j = 1; j < nv - 1; ++j) {
            a[j] = -theta_dt * c.a_v[j][i];
            b[j] = 1.0 - theta_dt * c.b_v[j][i];
            cc[j] = -theta_dt * c.c_v[j][i];
            d[j] = V[grid.idx(i, j)];
        }
        auto x = thomas_solve(a, b, cc, d);
        for (Size j = 0; j < nv; ++j) V[grid.idx(i, j)] = x[j];
    }
}

// Polymorphic base for ADI schemes
class ADISchemeBase {
public:
    virtual ~ADISchemeBase() = default;
    virtual void step(const std::vector<Real>& V_old, std::vector<Real>& V_new,
                      Real dt, const FDMGrid2D& grid,
                      const HestonOperatorCoeffs& c) const = 0;
    // Reset internal state (e.g., Rannacher step counter) for a new price() call.
    // Default: no-op (stateless schemes). RannacherSmoothing2D overrides this.
    virtual void reset() const {}
};

// Craig-Sneyd scheme (second-order, theta = 0.5 default)
// Explicitly treats mixed derivative L_xv in both correction stages.
class CraigSneydScheme : public ADISchemeBase {
public:
    explicit CraigSneydScheme(Real theta = 0.5) : theta_(theta) {}

    void step(const std::vector<Real>& V_old, std::vector<Real>& V_new,
              Real dt, const FDMGrid2D& grid,
              const HestonOperatorCoeffs& c) const override {
        Real theta_dt = theta_ * dt;
        Real corr_coeff = (theta_ - 0.5) * dt;
        Size total = V_old.size();

        // L_x(V_n), L_v(V_n), L_xv(V_n) — needed for implicit corrections
        std::vector<Real> LxVn, LvVn, LxvVn;
        apply_x_operator(V_old, grid, c, LxVn);
        apply_v_operator(V_old, grid, c, LvVn);
        apply_mixed_derivative(V_old, grid, c, LxvVn);

        // Y0 = V_n + dt * L_full * V_n
        // Use apply_heston_operator (not component sum) to keep boundary handling
        // consistent: boundaries return 0 so Y0 = V_old at boundaries, matching
        // the external boundary-condition evolution in PDEEngine2D.
        // Component sum (Lx+Lv+Lxv+L0) would evolve boundaries inconsistently
        // because apply_x_operator / apply_v_operator populate boundary rows,
        // leading to numerical blow-up in subsequent L_full(Y0) applications.
        std::vector<Real> LVn;
        apply_heston_operator(V_old, grid, c, LVn);
        std::vector<Real> Y0(total);
        for (Size k = 0; k < total; ++k) Y0[k] = V_old[k] + dt * LVn[k];

        // (I - theta_dt L_x) Y1 = Y0 - theta_dt L_x V_n
        std::vector<Real> Y1 = Y0;
        for (Size k = 0; k < total; ++k) Y1[k] -= theta_dt * LxVn[k];
        solve_x_direction(Y1, grid, c, theta_dt);

        // (I - theta_dt L_v) Y2 = Y1 - theta_dt L_v V_n
        std::vector<Real> Y2 = Y1;
        for (Size k = 0; k < total; ++k) Y2[k] -= theta_dt * LvVn[k];
        solve_v_direction(Y2, grid, c, theta_dt);

        // Yhat = Y2 + theta_dt L_xv V_n + (theta - 0.5) dt L_full Y0
        std::vector<Real> LY0;
        apply_heston_operator(Y0, grid, c, LY0);
        std::vector<Real> Yhat = Y2;
        for (Size k = 0; k < total; ++k) {
            Yhat[k] += theta_dt * LxvVn[k] + corr_coeff * LY0[k];
        }

        // ytilde = Yhat - Y2 + Y1
        std::vector<Real> ytilde(total);
        for (Size k = 0; k < total; ++k) ytilde[k] = Yhat[k] - Y2[k] + Y1[k];

        // L_x(Yhat), L_v(Yhat), L_xv(Yhat)
        std::vector<Real> LxYhat, LvYhat, LxvYhat;
        apply_x_operator(Yhat, grid, c, LxYhat);
        apply_v_operator(Yhat, grid, c, LvYhat);
        apply_mixed_derivative(Yhat, grid, c, LxvYhat);

        // (I - theta_dt L_x) ytilde1 = ytilde - theta_dt L_x Yhat
        std::vector<Real> ytilde1 = ytilde;
        for (Size k = 0; k < total; ++k) ytilde1[k] -= theta_dt * LxYhat[k];
        solve_x_direction(ytilde1, grid, c, theta_dt);

        // (I - theta_dt L_v) ytilde2 = ytilde1 - theta_dt L_v Yhat
        std::vector<Real> ytilde2 = ytilde1;
        for (Size k = 0; k < total; ++k) ytilde2[k] -= theta_dt * LvYhat[k];
        solve_v_direction(ytilde2, grid, c, theta_dt);

        // V_{n+1} = ytilde2 + theta_dt L_xv Yhat + (theta - 0.5) dt L_full Yhat
        std::vector<Real> LYhat;
        apply_heston_operator(Yhat, grid, c, LYhat);
        V_new.assign(total, 0.0);
        for (Size k = 0; k < total; ++k) {
            V_new[k] = ytilde2[k] + theta_dt * LxvYhat[k] + corr_coeff * LYhat[k];
        }
    }

private:
    Real theta_;
};

// Hundsdorfer-Verwer scheme (second-order, theta = 0.5 default; optimal theta = 0.5 + 1/sqrt(12))
// Folds mixed derivative into full operator L (no explicit A2 term in correction).
class HundsdorferVerwerScheme : public ADISchemeBase {
public:
    explicit HundsdorferVerwerScheme(Real theta = 0.5) : theta_(theta) {}

    void step(const std::vector<Real>& V_old, std::vector<Real>& V_new,
              Real dt, const FDMGrid2D& grid,
              const HestonOperatorCoeffs& c) const override {
        Real theta_dt = theta_ * dt;
        Real corr_coeff = (theta_ - 0.5) * dt;
        Size total = V_old.size();

        std::vector<Real> LxVn, LvVn;
        apply_x_operator(V_old, grid, c, LxVn);
        apply_v_operator(V_old, grid, c, LvVn);

        // Y0 = V_n + dt L_full V_n
        std::vector<Real> LVn;
        apply_heston_operator(V_old, grid, c, LVn);
        std::vector<Real> Y0(total);
        for (Size k = 0; k < total; ++k) Y0[k] = V_old[k] + dt * LVn[k];

        // (I - theta_dt L_x) Y1 = Y0 - theta_dt L_x V_n
        std::vector<Real> Y1 = Y0;
        for (Size k = 0; k < total; ++k) Y1[k] -= theta_dt * LxVn[k];
        solve_x_direction(Y1, grid, c, theta_dt);

        // (I - theta_dt L_v) Y2 = Y1 - theta_dt L_v V_n
        std::vector<Real> Y2 = Y1;
        for (Size k = 0; k < total; ++k) Y2[k] -= theta_dt * LvVn[k];
        solve_v_direction(Y2, grid, c, theta_dt);

        // Yhat = Y2 + (theta - 0.5) dt L_full Y0
        std::vector<Real> LY0;
        apply_heston_operator(Y0, grid, c, LY0);
        std::vector<Real> Yhat = Y2;
        for (Size k = 0; k < total; ++k) Yhat[k] += corr_coeff * LY0[k];

        // ytilde = Yhat - Y2 + Y1
        std::vector<Real> ytilde(total);
        for (Size k = 0; k < total; ++k) ytilde[k] = Yhat[k] - Y2[k] + Y1[k];

        std::vector<Real> LxYhat, LvYhat;
        apply_x_operator(Yhat, grid, c, LxYhat);
        apply_v_operator(Yhat, grid, c, LvYhat);

        // (I - theta_dt L_x) ytilde1 = ytilde - theta_dt L_x Yhat
        std::vector<Real> ytilde1 = ytilde;
        for (Size k = 0; k < total; ++k) ytilde1[k] -= theta_dt * LxYhat[k];
        solve_x_direction(ytilde1, grid, c, theta_dt);

        // (I - theta_dt L_v) ytilde2 = ytilde1 - theta_dt L_v Yhat
        std::vector<Real> ytilde2 = ytilde1;
        for (Size k = 0; k < total; ++k) ytilde2[k] -= theta_dt * LvYhat[k];
        solve_v_direction(ytilde2, grid, c, theta_dt);

        // V_{n+1} = ytilde2 + (theta - 0.5) dt L_full Yhat
        std::vector<Real> LYhat;
        apply_heston_operator(Yhat, grid, c, LYhat);
        V_new.assign(total, 0.0);
        for (Size k = 0; k < total; ++k) {
            V_new[k] = ytilde2[k] + corr_coeff * LYhat[k];
        }
    }

private:
    Real theta_;
};

// Modified Craig-Sneyd scheme (second-order, theta = 0.5 default; optimal theta = 0.5 + 1/sqrt(12))
//
// Identical to CS EXCEPT in the final correction step:
//   CS:  V_{n+1} = ytilde2 + theta*dt*L_xv*Yhat + (theta-0.5)*dt*L_full*Yhat
//   MCS: V_{n+1} = ytilde2 + theta*dt*L_xv*Yhat + (theta-0.5)*dt*L_full*ytilde2
// At theta=0.5, (theta-0.5)=0 => CS == MCS (identical).
// At theta!=0.5, MCS replaces L_full(Yhat) with L_full(ytilde2).
// Reference: in 't Hout & Foulon (2010) Table 2, in 't Hout & Welfert (2007).
// Unconditionally stable and second-order for theta in [0.5, 1].
class ModifiedCraigSneydScheme : public ADISchemeBase {
public:
    explicit ModifiedCraigSneydScheme(Real theta = 0.5) : theta_(theta) {}

    void step(const std::vector<Real>& V_old, std::vector<Real>& V_new,
              Real dt, const FDMGrid2D& grid,
              const HestonOperatorCoeffs& c) const override {
        Real theta_dt = theta_ * dt;
        Real corr_coeff = (theta_ - 0.5) * dt;
        Size total = V_old.size();

        // L_x(V_n), L_v(V_n), L_xv(V_n)
        std::vector<Real> LxVn, LvVn, LxvVn;
        apply_x_operator(V_old, grid, c, LxVn);
        apply_v_operator(V_old, grid, c, LvVn);
        apply_mixed_derivative(V_old, grid, c, LxvVn);

        // Y0 = V_n + dt * L_full * V_n
        std::vector<Real> LVn;
        apply_heston_operator(V_old, grid, c, LVn);
        std::vector<Real> Y0(total);
        for (Size k = 0; k < total; ++k) Y0[k] = V_old[k] + dt * LVn[k];

        // (I - theta_dt L_x) Y1 = Y0 - theta_dt L_x V_n
        std::vector<Real> Y1 = Y0;
        for (Size k = 0; k < total; ++k) Y1[k] -= theta_dt * LxVn[k];
        solve_x_direction(Y1, grid, c, theta_dt);

        // (I - theta_dt L_v) Y2 = Y1 - theta_dt L_v V_n
        std::vector<Real> Y2 = Y1;
        for (Size k = 0; k < total; ++k) Y2[k] -= theta_dt * LvVn[k];
        solve_v_direction(Y2, grid, c, theta_dt);

        // Yhat = Y2 + theta_dt L_xv V_n + (theta - 0.5) dt L_full Y0
        std::vector<Real> LY0;
        apply_heston_operator(Y0, grid, c, LY0);
        std::vector<Real> Yhat = Y2;
        for (Size k = 0; k < total; ++k) {
            Yhat[k] += theta_dt * LxvVn[k] + corr_coeff * LY0[k];
        }

        // ytilde = Yhat - Y2 + Y1
        std::vector<Real> ytilde(total);
        for (Size k = 0; k < total; ++k) ytilde[k] = Yhat[k] - Y2[k] + Y1[k];

        // L_x(Yhat), L_v(Yhat), L_xv(Yhat)
        std::vector<Real> LxYhat, LvYhat, LxvYhat;
        apply_x_operator(Yhat, grid, c, LxYhat);
        apply_v_operator(Yhat, grid, c, LvYhat);
        apply_mixed_derivative(Yhat, grid, c, LxvYhat);

        // (I - theta_dt L_x) ytilde1 = ytilde - theta_dt L_x Yhat
        std::vector<Real> ytilde1 = ytilde;
        for (Size k = 0; k < total; ++k) ytilde1[k] -= theta_dt * LxYhat[k];
        solve_x_direction(ytilde1, grid, c, theta_dt);

        // (I - theta_dt L_v) ytilde2 = ytilde1 - theta_dt L_v Yhat
        std::vector<Real> ytilde2 = ytilde1;
        for (Size k = 0; k < total; ++k) ytilde2[k] -= theta_dt * LvYhat[k];
        solve_v_direction(ytilde2, grid, c, theta_dt);

        // MCS final correction: replace L_full(Yhat) with L_full(ytilde2)
        // CS:  V_{n+1} = ytilde2 + theta_dt*L_xv*Yhat + corr_coeff*L_full*Yhat
        // MCS: V_{n+1} = ytilde2 + theta_dt*L_xv*Yhat + corr_coeff*L_full*ytilde2
        V_new.assign(total, 0.0);
        if (corr_coeff == 0.0) {
            // theta == 0.5: corr term vanishes, identical to CS (skip L_full(ytilde2))
            for (Size k = 0; k < total; ++k) {
                V_new[k] = ytilde2[k] + theta_dt * LxvYhat[k];
            }
        } else {
            std::vector<Real> Lytilde2;
            apply_heston_operator(ytilde2, grid, c, Lytilde2);
            for (Size k = 0; k < total; ++k) {
                V_new[k] = ytilde2[k] + theta_dt * LxvYhat[k] + corr_coeff * Lytilde2[k];
            }
        }
    }

private:
    Real theta_;
};

// Rannacher smoothing wrapper for 2D ADI schemes.
//
// First n_warmup time steps use the warmup scheme (typically ADI with theta=1.0,
// which is L-stable and strongly damps high-frequency oscillations from non-smooth
// payoff kinks at the strike). After warmup, switches to the main scheme
// (typically theta=0.5 for second-order accuracy).
//
// This mirrors the 1D RannacherSmoothing in fdm_scheme.hpp: the warmup phase
// damps the spurious oscillations that Crank-Nicolson (theta=0.5) produces
// near non-smooth regions (kinks, discontinuities) at coarse time grids.
//
// Reference: Rannacher (1984), Giles & Carter (1988); ADI adaptation follows
// the theta-switching approach (theta=1 warmup → theta=0.5 main).
//
// Note: The wrapper is stateful (tracks step count via mutable member) because
// the ADISchemeBase::step interface is const. The step counter resets on each
// new price() call via reset().
class RannacherSmoothing2D : public ADISchemeBase {
public:
    RannacherSmoothing2D(std::unique_ptr<ADISchemeBase> warmup_scheme,
                         std::unique_ptr<ADISchemeBase> main_scheme,
                         Size n_warmup)
        : warmup_scheme_(std::move(warmup_scheme)),
          main_scheme_(std::move(main_scheme)),
          n_warmup_(n_warmup),
          step_count_(0) {}

    void step(const std::vector<Real>& V_old, std::vector<Real>& V_new,
              Real dt, const FDMGrid2D& grid,
              const HestonOperatorCoeffs& c) const override {
        if (step_count_ < n_warmup_ && warmup_scheme_) {
            warmup_scheme_->step(V_old, V_new, dt, grid, c);
        } else {
            main_scheme_->step(V_old, V_new, dt, grid, c);
        }
        ++step_count_;
    }

    // Reset step counter for a new price() call.
    void reset() const override { step_count_ = 0; }

    Size warmup_steps() const noexcept { return n_warmup_; }
    Size current_step() const noexcept { return step_count_; }

private:
    std::unique_ptr<ADISchemeBase> warmup_scheme_;
    std::unique_ptr<ADISchemeBase> main_scheme_;
    Size n_warmup_;
    mutable Size step_count_;
};

}  // namespace v1
}  // namespace cpphub
