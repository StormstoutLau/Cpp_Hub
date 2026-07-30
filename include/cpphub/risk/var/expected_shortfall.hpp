#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace cpphub {
inline namespace v1 {

class ExpectedShortfall {
public:
    Real from_losses(const std::vector<Real>& losses, Real confidence) const {
        if (losses.empty()) return 0;
        std::vector<Real> sorted = losses;
        std::sort(sorted.begin(), sorted.end());
        Size n = sorted.size();
        Size tail_start = static_cast<Size>(confidence * n);
        if (tail_start >= n) tail_start = n - 1;
        Real sum = 0;
        for (Size i = tail_start; i < n; ++i) {
            sum += sorted[i];
        }
        Size tail_count = n - tail_start;
        if (tail_count == 0) return sorted[n - 1];
        return sum / static_cast<Real>(tail_count);
    }

    Real normal_es(Real mean, Real sigma, Real confidence) const {
        Real z = inv_normal_cdf(1.0 - confidence);
        Real phi_z = normal_pdf(z);
        return -mean + sigma * phi_z / (1.0 - confidence);
    }

    Real student_t_es(Real mean, Real sigma, Real dof, Real confidence) const {
        if (dof <= 1) dof = 2.1;
        Real z = inv_normal_cdf(1.0 - confidence);
        Real x = z;
        for (int iter = 0; iter < 100; ++iter) {
            bool neg = (x < 0);
            Real x_beta = dof / (dof + x * x);
            Real ib = betainc(x_beta, dof / 2.0, 0.5);
            Real cdf_val = neg ? 0.5 * ib : 1.0 - 0.5 * ib;
            Real f = cdf_val - (1.0 - confidence);
            if (std::abs(f) < 1e-15) break;
            Real pdf_val = std::exp(std::lgamma((dof + 1.0) / 2.0) - std::lgamma(dof / 2.0))
                           * std::pow(1.0 + x * x / dof, -(dof + 1.0) / 2.0)
                           / std::sqrt(dof * PI);
            if (std::abs(pdf_val) < 1e-30) break;
            x -= f / pdf_val;
        }
        Real t_q = x;
        Real pdf_q = std::exp(std::lgamma((dof + 1.0) / 2.0) - std::lgamma(dof / 2.0))
                     * std::pow(1.0 + t_q * t_q / dof, -(dof + 1.0) / 2.0)
                     / std::sqrt(dof * PI);
        Real es_factor = pdf_q / (1.0 - confidence) * (dof + t_q * t_q) / (dof - 1.0);
        Real scale = std::sqrt(sigma * sigma * (dof - 2.0) / dof);
        return -mean + scale * es_factor;
    }

    Real cornish_fisher_es(Real mean, Real sigma, Real skew, Real kurt,
                             Real confidence) const {
        Real z = inv_normal_cdf(1.0 - confidence);
        Real z2 = z * z;
        Real z3 = z2 * z;
        Real z_cf = z + (z2 - 1.0) * skew / 6.0
                     + (z3 - 3.0 * z) * kurt / 24.0
                     - (2.0 * z3 - 5.0 * z) * skew * skew / 36.0;
        Real phi_z = normal_pdf(z_cf);
        return -mean + sigma * phi_z / (1.0 - confidence);
    }

    Real from_mc_paths(const std::vector<Real>& pnl_paths, Real confidence) const {
        return from_losses(pnl_paths, confidence);
    }

    Real tail_average(const std::vector<Real>& losses, Real confidence,
                       Size n_tail_points = 100) const {
        if (losses.empty()) return 0;
        std::vector<Real> sorted = losses;
        std::sort(sorted.begin(), sorted.end());
        Size n = sorted.size();
        Real sum = 0;
        for (Size k = 1; k <= n_tail_points; ++k) {
            Real p = confidence + (1.0 - confidence) * static_cast<Real>(k) / static_cast<Real>(n_tail_points + 1);
            Real index = p * (n - 1);
            if (index >= n - 1) index = static_cast<Real>(n - 1);
            Size lo = static_cast<Size>(std::floor(index));
            Size hi = static_cast<Size>(std::ceil(index));
            if (lo >= n) lo = n - 1;
            if (hi >= n) hi = n - 1;
            Real val;
            if (lo == hi) {
                val = sorted[lo];
            } else {
                Real frac = index - lo;
                val = sorted[lo] + frac * (sorted[hi] - sorted[lo]);
            }
            sum += val;
        }
        return sum / static_cast<Real>(n_tail_points);
    }

private:
    static Real betainc(Real x, Real a, Real b) {
        if (x <= 0) return 0;
        if (x >= 1) return 1;
        if (x > (a + 1.0) / (a + b + 2.0)) {
            return 1.0 - betainc(1.0 - x, b, a);
        }
        Real lbeta = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
        Real front = std::exp(a * std::log(x) + b * std::log(1.0 - x) - lbeta) / a;
        Real f = 1.0;
        Real C = 1.0;
        Real D = 0.0;
        Real tiny = 1e-30;
        for (int m = 1; m <= 200; ++m) {
            Real a_m;
            if (m % 2 == 0) {
                int k = m / 2;
                a_m = static_cast<Real>(k) * (b - static_cast<Real>(k)) * x
                      / ((a + static_cast<Real>(m) - 1.0) * (a + static_cast<Real>(m)));
            } else {
                int k = (m - 1) / 2;
                a_m = -(a + static_cast<Real>(k)) * (a + b + static_cast<Real>(k)) * x
                      / ((a + static_cast<Real>(m) - 1.0) * (a + static_cast<Real>(m)));
            }
            D = 1.0 + a_m * D;
            if (D == 0) D = tiny;
            C = 1.0 + a_m / C;
            if (C == 0) C = tiny;
            D = 1.0 / D;
            Real delta = C * D;
            f *= delta;
            if (std::abs(delta - 1.0) < 1e-14) break;
        }
        return front / f;
    }
};

}  // namespace v1
}  // namespace cpphub
