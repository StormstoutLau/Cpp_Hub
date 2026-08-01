#pragma once
#include <vector>
#include <functional>
#include <cmath>
#include "cpphub/core/types.hpp"
#include "cpphub/pricing/pde/fdm_grid.hpp"
#include "cpphub/pricing/pde/thomas_solver.hpp"

namespace cpphub {
inline namespace v1 {

enum class FDMSchemeType {
    ExplicitEuler,
    ImplicitEuler,
    CrankNicolson,
    RannacherSmoothing
};

struct PDEParams {
    Real r;
    Real q;
    Real sigma;
    Real T;
    Real K;
    Real S0;
};

class FDMScheme {
public:
    virtual ~FDMScheme() = default;
    virtual void step(const std::vector<Real>& V_old,
                       std::vector<Real>& V_new,
                       Real dt,
                       const FDMGrid& grid,
                       const PDEParams& params) const = 0;
};

inline void compute_operator_coeffs(const FDMGrid& grid, const PDEParams& params,
                                     std::vector<Real>& a, std::vector<Real>& b,
                                     std::vector<Real>& c) {
    Size n = grid.size();
    a.resize(n, 0.0); b.resize(n, 0.0); c.resize(n, 0.0);
    Real sigma2 = params.sigma * params.sigma;
    Real r_minus_q = params.r - params.q;

    for (Size i = 1; i < n - 1; ++i) {
        Real ds_p = grid.s(i + 1) - grid.s(i);
        Real ds_m = grid.s(i) - grid.s(i - 1);
        Real ds_av = grid.ds_avg(i);
        Real Si = grid.s(i);

        Real d1_left = -1.0 / (2.0 * ds_av);
        Real d1_right = 1.0 / (2.0 * ds_av);

        Real d2_left = 1.0 / (ds_av * ds_m);
        Real d2_mid = -2.0 / (ds_p * ds_m);
        Real d2_right = 1.0 / (ds_av * ds_p);

        Real half_sigma2_S2 = 0.5 * sigma2 * Si * Si;
        Real rq_S = r_minus_q * Si;

        a[i] = half_sigma2_S2 * d2_left + rq_S * d1_left;
        b[i] = half_sigma2_S2 * d2_mid - params.r;
        c[i] = half_sigma2_S2 * d2_right + rq_S * d1_right;
    }
}

class ExplicitEuler : public FDMScheme {
public:
    void step(const std::vector<Real>& V_old, std::vector<Real>& V_new, Real dt,
              const FDMGrid& grid, const PDEParams& params) const override {
        Size n = grid.size();
        std::vector<Real> a, b_op, c;
        compute_operator_coeffs(grid, params, a, b_op, c);

        for (Size i = 1; i < n - 1; ++i) {
            V_new[i] = V_old[i] + dt * (a[i] * V_old[i - 1] +
                                         b_op[i] * V_old[i] +
                                         c[i] * V_old[i + 1]);
        }
    }
};

class ImplicitEuler : public FDMScheme {
public:
    void step(const std::vector<Real>& V_old, std::vector<Real>& V_new, Real dt,
              const FDMGrid& grid, const PDEParams& params) const override {
        Size n = grid.size();
        std::vector<Real> a_op, b_op, c_op;
        compute_operator_coeffs(grid, params, a_op, b_op, c_op);

        std::vector<Real> a(n, 0.0), b(n, 0.0), c(n, 0.0), d(n, 0.0);

        b[0] = 1.0; c[0] = 0.0; d[0] = V_new[0];
        a[n - 1] = 0.0; b[n - 1] = 1.0; d[n - 1] = V_new[n - 1];

        for (Size i = 1; i < n - 1; ++i) {
            a[i] = -dt * a_op[i];
            b[i] = 1.0 - dt * b_op[i];
            c[i] = -dt * c_op[i];
            d[i] = V_old[i];
        }

        V_new = thomas_solve(a, b, c, d);
    }
};

class CrankNicolson : public FDMScheme {
public:
    explicit CrankNicolson(Real theta = 0.5) : theta_(theta) {}
    void step(const std::vector<Real>& V_old, std::vector<Real>& V_new, Real dt,
              const FDMGrid& grid, const PDEParams& params) const override {
        Size n = grid.size();
        std::vector<Real> a_op, b_op, c_op;
        compute_operator_coeffs(grid, params, a_op, b_op, c_op);

        std::vector<Real> a(n, 0.0), b(n, 0.0), c(n, 0.0), d(n, 0.0);

        b[0] = 1.0; c[0] = 0.0; d[0] = V_new[0];
        a[n - 1] = 0.0; b[n - 1] = 1.0; d[n - 1] = V_new[n - 1];

        Real dt_theta = dt * theta_;
        Real dt_1mtheta = dt * (1.0 - theta_);

        for (Size i = 1; i < n - 1; ++i) {
            a[i] = -dt_theta * a_op[i];
            b[i] = 1.0 - dt_theta * b_op[i];
            c[i] = -dt_theta * c_op[i];
            d[i] = V_old[i] + dt_1mtheta * (a_op[i] * V_old[i - 1] +
                                             b_op[i] * V_old[i] +
                                             c_op[i] * V_old[i + 1]);
        }

        V_new = thomas_solve(a, b, c, d);
    }
private:
    Real theta_;
};

// Rannacher smoothing: 前 n_warmup 步用 ImplicitEuler (L 稳定, damping 高频振荡),
// 之后切换到 CrankNicolson. 用于处理非光滑 payoff (如行权价处 kink, 数字期权)
// 导致的 CrankNicolson 振荡. 参考 Rannacher (1984), Giles & Carter (1988).
// 简化版: warmup 阶段用 ImplicitEuler(dt) 全步, 保持总时间 = n_steps * dt = T.
// 严格版会用 ImplicitEuler(dt/2) 半步 (需额外步数补偿).
class RannacherSmoothing : public FDMScheme {
public:
    explicit RannacherSmoothing(Size n_warmup = 4, Real theta = 0.5)
        : n_warmup_(n_warmup), theta_(theta), step_count_(0) {}

    void reset() const { step_count_ = 0; }

    void step(const std::vector<Real>& V_old, std::vector<Real>& V_new, Real dt,
              const FDMGrid& grid, const PDEParams& params) const override {
        if (step_count_ < n_warmup_) {
            // warmup: ImplicitEuler 全步 (L 稳定, damping 非光滑初值的高频分量)
            ImplicitEuler impl;
            impl.step(V_old, V_new, dt, grid, params);
        } else {
            // 主循环: CrankNicolson 二阶精度
            CrankNicolson cn(theta_);
            cn.step(V_old, V_new, dt, grid, params);
        }
        ++step_count_;
    }

    Size warmup_steps() const noexcept { return n_warmup_; }
    Size current_step() const noexcept { return step_count_; }

private:
    Size n_warmup_;
    Real theta_;
    mutable Size step_count_;
};

}  // namespace v1
}  // namespace cpphub
