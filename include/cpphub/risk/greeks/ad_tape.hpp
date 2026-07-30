#pragma once
#include <autodiff/reverse/var.hpp>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"

namespace cpphub::v1 {

using autodiff::var;
using autodiff::wrt;
using autodiff::derivatives;
using autodiff::val;

inline var normal_cdf_var(var x) {
    return var(0.5) * (var(1.0) + erf(x / sqrt(var(2.0))));
}

inline var normal_pdf_var(var x) {
    return exp(-x * x / var(2.0)) / sqrt(var(2.0 * PI));
}

inline var bsm_d1_var(var S, var K, var T, var r, var q, var sigma) {
    return (log(S / K) + (r - q + sigma * sigma / var(2.0)) * T) / (sigma * sqrt(T));
}

inline var bsm_d2_var(var d1, var sigma, var T) {
    return d1 - sigma * sqrt(T);
}

}  // namespace cpphub::v1
