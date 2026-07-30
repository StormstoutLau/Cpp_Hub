#pragma once
#include <vector>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {

inline std::vector<Real> thomas_solve(const std::vector<Real>& a,
                                       const std::vector<Real>& b,
                                       const std::vector<Real>& c,
                                       const std::vector<Real>& d) {
    Size n = b.size();
    std::vector<Real> cp(n, 0.0), dp(n, 0.0), x(n, 0.0);
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for (Size i = 1; i < n; ++i) {
        Real denom = b[i] - a[i] * cp[i - 1];
        if (i < n - 1) cp[i] = c[i] / denom;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }
    x[n - 1] = dp[n - 1];
    for (Size i = n - 1; i-- > 0; ) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }
    return x;
}

}  // namespace v1
}  // namespace cpphub
