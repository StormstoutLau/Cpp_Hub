#pragma once
#include <vector>
#include <memory>
#include <cmath>
#include <string>
#include <functional>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/instruments/payoff/payoff.hpp"

namespace cpphub {
inline namespace v1 {

enum class BinomialType {
    CRR,
    JarrowRudd,
    Tian,
    LeisenReimer
};

struct BinomialParams {
    Real S0;
    Real K;
    Real T;
    Real r;
    Real q;
    Real sigma;
    Size n_steps;
};

class BinomialTreeEngine {
public:
    BinomialTreeEngine(BinomialParams params, BinomialType type = BinomialType::CRR)
        : params_(params), type_(type) {
        compute_tree_parameters();
    }

    Real price_european(const PayOff& payoff) const {
        Size n = params_.n_steps;
        std::vector<Real> V(n + 1);
        for (Size j = 0; j <= n; ++j) {
            Real S = params_.S0 * std::pow(u_, n - j) * std::pow(d_, j);
            V[j] = payoff(S);
        }
        for (Size i = n; i > 0; --i) {
            for (Size j = 0; j < i; ++j) {
                V[j] = disc_ * (p_ * V[j] + (1.0 - p_) * V[j + 1]);
            }
        }
        return V[0];
    }

    Real price_american(const PayOff& payoff) const {
        Size n = params_.n_steps;
        std::vector<Real> V(n + 1);
        for (Size j = 0; j <= n; ++j) {
            Real S = params_.S0 * std::pow(u_, n - j) * std::pow(d_, j);
            V[j] = payoff(S);
        }
        for (Size i = n; i > 0; --i) {
            for (Size j = 0; j < i; ++j) {
                Real S = params_.S0 * std::pow(u_, i - 1 - j) * std::pow(d_, j);
                Real continuation = disc_ * (p_ * V[j] + (1.0 - p_) * V[j + 1]);
                V[j] = std::max(continuation, payoff(S));
            }
        }
        return V[0];
    }

    Real price_bermudan(const PayOff& payoff, const std::vector<Size>& exercise_steps) const {
        Size n = params_.n_steps;
        std::vector<bool> is_exercise(n, false);
        for (auto s : exercise_steps) {
            if (s < n) is_exercise[s] = true;
        }
        std::vector<Real> V(n + 1);
        for (Size j = 0; j <= n; ++j) {
            Real S = params_.S0 * std::pow(u_, n - j) * std::pow(d_, j);
            V[j] = payoff(S);
        }
        for (Size i = n; i > 0; --i) {
            for (Size j = 0; j < i; ++j) {
                Real S = params_.S0 * std::pow(u_, i - 1 - j) * std::pow(d_, j);
                Real continuation = disc_ * (p_ * V[j] + (1.0 - p_) * V[j + 1]);
                if (is_exercise[i - 1]) {
                    V[j] = std::max(continuation, payoff(S));
                } else {
                    V[j] = continuation;
                }
            }
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

        BinomialParams params_up = params_;
        params_up.S0 = params_.S0 + dS;
        BinomialTreeEngine engine_up(params_up, type_);
        Real V_up = american ? engine_up.price_american(payoff) : engine_up.price_european(payoff);

        BinomialParams params_down = params_;
        params_down.S0 = params_.S0 - dS;
        BinomialTreeEngine engine_down(params_down, type_);
        Real V_down = american ? engine_down.price_american(payoff) : engine_down.price_european(payoff);

        Greeks g;
        g.delta = (V_up - V_down) / (2.0 * dS);
        g.gamma = (V_up - 2.0 * V0 + V_down) / (dS * dS);

        Real dt = params_.T / static_cast<Real>(params_.n_steps);
        if (params_.n_steps > 1) {
            BinomialParams params_theta = params_;
            params_theta.T = params_.T - dt;
            params_theta.n_steps = params_.n_steps - 1;
            BinomialTreeEngine engine_theta(params_theta, type_);
            Real V_theta = american ? engine_theta.price_american(payoff) : engine_theta.price_european(payoff);
            g.theta = (V_theta - V0) / dt;
        } else {
            g.theta = 0.0;
        }

        return g;
    }

    std::string name() const { return "BinomialTreeEngine"; }
    BinomialType type() const { return type_; }
    const BinomialParams& params() const { return params_; }

private:
    BinomialParams params_;
    BinomialType type_;
    Real u_ = 1.0;
    Real d_ = 1.0;
    Real p_ = 0.5;
    Real disc_ = 1.0;

    void compute_tree_parameters() {
        Real dt = params_.T / static_cast<Real>(params_.n_steps);
        if (dt <= 0.0 || params_.n_steps == 0) {
            throw std::invalid_argument("BinomialTreeEngine: n_steps must be positive");
        }
        disc_ = std::exp(-params_.r * dt);
        switch (type_) {
            case BinomialType::CRR: compute_crr(); break;
            case BinomialType::JarrowRudd: compute_jarrow_rudd(); break;
            case BinomialType::Tian: compute_tian(); break;
            case BinomialType::LeisenReimer: compute_leisen_reimer(); break;
        }
    }

    void compute_crr() {
        Real dt = params_.T / static_cast<Real>(params_.n_steps);
        u_ = std::exp(params_.sigma * std::sqrt(dt));
        d_ = 1.0 / u_;
        Real a = std::exp((params_.r - params_.q) * dt);
        p_ = (a - d_) / (u_ - d_);
    }

    void compute_jarrow_rudd() {
        Real dt = params_.T / static_cast<Real>(params_.n_steps);
        Real drift = (params_.r - params_.q - 0.5 * params_.sigma * params_.sigma) * dt;
        Real vol = params_.sigma * std::sqrt(dt);
        u_ = std::exp(drift + vol);
        d_ = std::exp(drift - vol);
        p_ = 0.5;
    }

    void compute_tian() {
        Real dt = params_.T / static_cast<Real>(params_.n_steps);
        Real R = std::exp((params_.r - params_.q) * dt);
        Real V = std::exp(params_.sigma * params_.sigma * dt);
        Real sqrt_term = std::sqrt(V * V + 2.0 * V - 3.0);
        u_ = 0.5 * R * V * (V + 1.0 + sqrt_term);
        d_ = 0.5 * R * V * (V + 1.0 - sqrt_term);
        p_ = (R - d_) / (u_ - d_);
    }

    void compute_leisen_reimer() {
        Size n = params_.n_steps;
        if (n % 2 == 0) n += 1;
        Real dt = params_.T / static_cast<Real>(n);

        Real d1 = (std::log(params_.S0 / params_.K) +
                   (params_.r - params_.q + 0.5 * params_.sigma * params_.sigma) * params_.T) /
                  (params_.sigma * std::sqrt(params_.T));
        Real d2 = d1 - params_.sigma * std::sqrt(params_.T);

        Real pbar = peizer_pratt_inversion(d1, n);
        p_ = peizer_pratt_inversion(d2, n);

        Real a = std::exp((params_.r - params_.q) * dt);
        u_ = a * pbar / p_;
        d_ = a * (1.0 - pbar) / (1.0 - p_);
        disc_ = std::exp(-params_.r * dt);
    }

    static Real peizer_pratt_inversion(Real z, Size n) {
        if (z == 0.0) return 0.5;
        Real n_real = static_cast<Real>(n);
        Real a = n_real + 1.0 / 3.0 + 0.1 / (n_real + 1.0);
        Real sign = (z > 0.0) ? 1.0 : -1.0;
        Real num = z * z * (n_real + 1.0 / 6.0);
        Real den = a * a;
        Real exp_arg = -num / den;
        if (exp_arg < -700.0) exp_arg = -700.0;
        return 0.5 + 0.5 * sign * std::sqrt(-std::expm1(exp_arg));
    }
};

}  // namespace v1
}  // namespace cpphub
