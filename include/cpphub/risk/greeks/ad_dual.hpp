#pragma once
#include <autodiff/forward/dual.hpp>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"

namespace cpphub {
inline namespace v1 {

using autodiff::dual;
using autodiff::wrt;
using autodiff::derivative;
using autodiff::val;

using dual2nd = autodiff::dual2nd;

inline dual normal_cdf_dual(dual x) {
    return dual(0.5) * (dual(1.0) + erf(x / sqrt(dual(2.0))));
}

inline dual normal_pdf_dual(dual x) {
    return exp(-x * x / dual(2.0)) / sqrt(dual(2.0 * PI));
}

inline dual bsm_d1_dual(dual S, dual K, dual T, dual r, dual q, dual sigma) {
    return (log(S / K) + (r - q + sigma * sigma / dual(2.0)) * T) / (sigma * sqrt(T));
}

inline dual bsm_d2_dual(dual d1, dual sigma, dual T) {
    return d1 - sigma * sqrt(T);
}

}  // namespace v1
}  // namespace cpphub
