#pragma once
#include <span>
#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/core/constants.hpp"

namespace cpphub {
inline namespace v1 {

struct GBMParams {
    Real S0;
    Real mu;
    Real sigma;
};

class GBM {
public:
    explicit GBM(GBMParams p) : params_(p) {}

    void generate_path(Real T, Size n_steps,
                       std::span<Real> path,
                       Philox4x64& rng) const {
        Real dt = T / static_cast<Real>(n_steps);
        Real drift = (params_.mu - Real(0.5) * params_.sigma * params_.sigma) * dt;
        Real vol = params_.sigma * std::sqrt(dt);
        path[0] = params_.S0;
        for (Size i = 1; i <= n_steps; ++i) {
            uint64_t r1 = rng();
            uint64_t r2 = rng();
            double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            auto [z1, z2] = box_muller(u1, u2);
            Real z = (i < n_steps) ? z1 : z2;
            path[i] = path[i - 1] * std::exp(drift + vol * z);
        }
    }

    Real evolve(Real S, Real dt, Real Z) const noexcept {
        return S * std::exp((params_.mu - Real(0.5) * params_.sigma * params_.sigma) * dt
                            + params_.sigma * std::sqrt(dt) * Z);
    }

    Complex log_characteristic_function(Complex u, Real t) const {
        Complex iu = Complex(0, 1) * u;
        Real log_S0 = std::log(params_.S0);
        Real mean = log_S0 + (params_.mu - Real(0.5) * params_.sigma * params_.sigma) * t;
        Real var = params_.sigma * params_.sigma * t;
        return std::exp(iu * mean - Real(0.5) * u * u * var);
    }

    const GBMParams& params() const { return params_; }
    Real spot() const { return params_.S0; }

private:
    GBMParams params_;
};

}  // namespace v1
}  // namespace cpphub
