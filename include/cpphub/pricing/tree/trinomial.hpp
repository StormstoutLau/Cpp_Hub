#pragma once
#include <vector>
#include <memory>
#include <cmath>
#include <string>
#include <stdexcept>
#include <algorithm>
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/instruments/payoff/payoff.hpp"

namespace cpphub {
inline namespace v1 {

enum class TrinomialType {
    Explicit,
    Implicit,
    Hybrid
};

struct TrinomialParams {
    Real S0;
    Real K;
    Real T;
    Real r;
    Real q;
    Real sigma;
    Size n_steps;
};

class TrinomialTreeEngine {
public:
    TrinomialTreeEngine(TrinomialParams params, TrinomialType type = TrinomialType::Explicit)
        : params_(params), type_(type) {
        compute_parameters();
    }

    Real price_european(const PayOff& payoff) const {
        Size n = params_.n_steps;
        if (n == 0) return 0.0;
        std::vector<Real> V(2 * n + 1, 0.0);
        for (Index j = -static_cast<Index>(n); j <= static_cast<Index>(n); ++j) {
            Real S = params_.S0 * std::pow(u_, j);
            V[j + n] = payoff(S);
        }
        for (Index i = n - 1; i >= 0; --i) {
            std::vector<Real> Vn(2 * n + 1, 0.0);
            Size cntr = static_cast<Size>(i + 1);
            for (Index j = -i; j <= i; ++j) {
                Real cont = disc_ * (pu_ * V[j + cntr + 1] + pm_ * V[j + cntr] + pd_ * V[j + cntr - 1]);
                Vn[j + i] = cont;
            }
            V.swap(Vn);
        }
        return V[0];
    }

    Real price_american(const PayOff& payoff) const {
        Size n = params_.n_steps;
        if (n == 0) return 0.0;
        std::vector<Real> V(2 * n + 1, 0.0);
        for (Index j = -static_cast<Index>(n); j <= static_cast<Index>(n); ++j) {
            Real S = params_.S0 * std::pow(u_, j);
            V[j + n] = payoff(S);
        }
        for (Index i = n - 1; i >= 0; --i) {
            std::vector<Real> Vn(2 * n + 1, 0.0);
            Size cntr = static_cast<Size>(i + 1);
            for (Index j = -i; j <= i; ++j) {
                Real S = params_.S0 * std::pow(u_, j);
                Real cont = disc_ * (pu_ * V[j + cntr + 1] + pm_ * V[j + cntr] + pd_ * V[j + cntr - 1]);
                Vn[j + i] = std::max(cont, payoff(S));
            }
            V.swap(Vn);
        }
        return V[0];
    }

    Real price_bermudan(const PayOff& payoff, const std::vector<Size>& exercise_steps) const {
        Size n = params_.n_steps;
        if (n == 0) return 0.0;
        std::vector<bool> is_exercise(n, false);
        for (auto s : exercise_steps) {
            if (s < n) is_exercise[s] = true;
        }
        std::vector<Real> V(2 * n + 1, 0.0);
        for (Index j = -static_cast<Index>(n); j <= static_cast<Index>(n); ++j) {
            Real S = params_.S0 * std::pow(u_, j);
            V[j + n] = payoff(S);
        }
        for (Index i = n - 1; i >= 0; --i) {
            std::vector<Real> Vn(2 * n + 1, 0.0);
            Size cntr = static_cast<Size>(i + 1);
            for (Index j = -i; j <= i; ++j) {
                Real S = params_.S0 * std::pow(u_, j);
                Real cont = disc_ * (pu_ * V[j + cntr + 1] + pm_ * V[j + cntr] + pd_ * V[j + cntr - 1]);
                if (is_exercise[i]) {
                    Vn[j + i] = std::max(cont, payoff(S));
                } else {
                    Vn[j + i] = cont;
                }
            }
            V.swap(Vn);
        }
        return V[0];
    }

    struct Greeks {
        Real delta;
        Real gamma;
        Real theta;
    };

    Greeks greeks(const PayOff& payoff, bool american = false) const {
        Real dS = params_.S0 * 0.01;
        Real V0 = american ? price_american(payoff) : price_european(payoff);

        TrinomialParams params_up = params_;
        params_up.S0 = params_.S0 + dS;
        TrinomialTreeEngine engine_up(params_up, type_);
        Real V_up = american ? engine_up.price_american(payoff) : engine_up.price_european(payoff);

        TrinomialParams params_down = params_;
        params_down.S0 = params_.S0 - dS;
        TrinomialTreeEngine engine_down(params_down, type_);
        Real V_down = american ? engine_down.price_american(payoff) : engine_down.price_european(payoff);

        Greeks g;
        g.delta = (V_up - V_down) / (2.0 * dS);
        g.gamma = (V_up - 2.0 * V0 + V_down) / (dS * dS);

        Real dt = params_.T / static_cast<Real>(params_.n_steps);
        if (params_.n_steps > 1) {
            TrinomialParams params_theta = params_;
            params_theta.T = params_.T - dt;
            params_theta.n_steps = params_.n_steps - 1;
            TrinomialTreeEngine engine_theta(params_theta, type_);
            Real V_theta = american ? engine_theta.price_american(payoff) : engine_theta.price_european(payoff);
            g.theta = (V_theta - V0) / dt;
        } else {
            g.theta = 0.0;
        }
        return g;
    }

    std::string name() const { return "TrinomialTreeEngine"; }
    TrinomialType type() const { return type_; }

private:
    TrinomialParams params_;
    TrinomialType type_;
    Real u_ = 1.0;
    Real pu_ = 0.0;
    Real pm_ = 0.0;
    Real pd_ = 0.0;
    Real disc_ = 1.0;
    Real dt_ = 0.0;

    void compute_parameters() {
        switch (type_) {
            case TrinomialType::Explicit: compute_explicit(); break;
            case TrinomialType::Implicit: compute_implicit(); break;
            case TrinomialType::Hybrid: compute_hybrid(); break;
        }
    }

    void compute_explicit() {
        dt_ = params_.T / static_cast<Real>(params_.n_steps);
        Real dx = params_.sigma * std::sqrt(dt_);
        u_ = std::exp(dx);
        Real dx2 = dx * dx;
        Real mu = params_.r - params_.q - 0.5 * params_.sigma * params_.sigma;
        pu_ = dt_ * (0.5 * params_.sigma * params_.sigma / dx2 + mu / (2.0 * dx));
        pd_ = dt_ * (0.5 * params_.sigma * params_.sigma / dx2 - mu / (2.0 * dx));
        pm_ = 1.0 - dt_ * params_.sigma * params_.sigma / dx2;
        disc_ = std::exp(-params_.r * dt_);
    }

    void compute_implicit() {
        dt_ = params_.T / static_cast<Real>(params_.n_steps);
        Real dx = params_.sigma * std::sqrt(3.0 * dt_);
        u_ = std::exp(dx);
        Real dx2 = dx * dx;
        pu_ = dt_ * (0.5 * params_.sigma * params_.sigma / dx2 + (params_.r - params_.q - 0.5 * params_.sigma * params_.sigma) / (2.0 * dx));
        pd_ = dt_ * (0.5 * params_.sigma * params_.sigma / dx2 - (params_.r - params_.q - 0.5 * params_.sigma * params_.sigma) / (2.0 * dx));
        pm_ = 1.0 - pu_ - pd_;
        disc_ = std::exp(-params_.r * dt_);
    }

    void compute_hybrid() {
        compute_explicit();
    }
};

}  // namespace v1
}  // namespace cpphub
