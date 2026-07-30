#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include <vector>
#include <cmath>
#include <functional>

namespace cpphub {
inline namespace v1 {

class BrownianBridge {
public:
    explicit BrownianBridge(Size n_steps, Real T = 1.0)
        : n_steps_(n_steps), T_(T), dt_(T / static_cast<Real>(n_steps))
    {
        build_schedule();
    }

    std::vector<Real> generate(const std::vector<Real>& uniforms) const
    {
        std::vector<Real> W(n_steps_ + 1, 0.0);
        for (Size s = 0; s < schedule_.size(); ++s) {
            Real Z = inv_normal_cdf(uniforms[s]);
            const auto& step = schedule_[s];
            Real mean = (static_cast<Real>(step.right - step.index) * W[step.left]
                       + static_cast<Real>(step.index - step.left) * W[step.right])
                      / static_cast<Real>(step.right - step.left);
            W[step.index] = mean + step.sigma * Z;
        }
        std::vector<Real> increments(n_steps_);
        for (Size i = 0; i < n_steps_; ++i) {
            increments[i] = W[i + 1] - W[i];
        }
        return increments;
    }

    std::vector<Real> generate_path(const std::vector<Real>& uniforms) const
    {
        std::vector<Real> W(n_steps_ + 1, 0.0);
        for (Size s = 0; s < schedule_.size(); ++s) {
            Real Z = inv_normal_cdf(uniforms[s]);
            const auto& step = schedule_[s];
            Real mean = (static_cast<Real>(step.right - step.index) * W[step.left]
                       + static_cast<Real>(step.index - step.left) * W[step.right])
                      / static_cast<Real>(step.right - step.left);
            W[step.index] = mean + step.sigma * Z;
        }
        std::vector<Real> path(n_steps_);
        for (Size i = 0; i < n_steps_; ++i) {
            path[i] = W[i + 1];
        }
        return path;
    }

    Size n_steps() const noexcept { return n_steps_; }

private:
    Size n_steps_;
    Real T_;
    Real dt_;

    struct BridgeStep {
        Size index;
        Size left;
        Size right;
        Real sigma;
    };

    std::vector<BridgeStep> schedule_;

    void build_schedule()
    {
        schedule_.clear();
        // First step: generate W(T) with direct Normal
        schedule_.push_back({n_steps_, 0, n_steps_, std::sqrt(T_)});

        // Recursively fill midpoints
        std::function<void(Size, Size)> bridge = [&](Size left, Size right) {
            if (right - left <= 1) return;
            Size mid = (left + right) / 2;
            Real sigma = std::sqrt(dt_ * static_cast<Real>(mid - left)
                                 * static_cast<Real>(right - mid)
                                 / static_cast<Real>(right - left));
            schedule_.push_back({mid, left, right, sigma});
            bridge(left, mid);
            bridge(mid, right);
        };
        bridge(0, n_steps_);
    }
};

}  // namespace v1
}  // namespace cpphub
