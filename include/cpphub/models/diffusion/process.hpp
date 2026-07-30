#pragma once
#include <complex>
#include <vector>
#include <span>
#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"

namespace cpphub {
inline namespace v1 {

class StochasticProcess {
public:
    virtual ~StochasticProcess() = default;
    virtual Size dimension() const = 0;
    virtual void generate_path(Real T, Size n_steps,
                                std::span<Real> path,
                                Philox4x64& rng) const = 0;
    virtual Complex characteristic_function(Complex u, Real tau) const {
        (void)u; (void)tau;
        return Complex{0.0, 0.0};
    }
    virtual Real spot() const = 0;
};

}  // namespace v1
}  // namespace cpphub
