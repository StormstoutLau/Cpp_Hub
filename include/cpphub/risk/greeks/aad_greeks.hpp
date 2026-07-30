#pragma once
#include "cpphub/risk/greeks/ad_tape.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include <random>
#include <cmath>

namespace cpphub::v1 {

struct AADGreeks {
    Real price;
    Real delta;
    Real vega;
    Real rho;
    Real theta;
    Real gamma;
};

class AADGreeksEngine {
public:
    static AADGreeks bsm_european(Real S, Real K, Real T, Real r, Real q,
                                   Real sigma, bool is_call);
    static AADGreeks heston_mc(Real S, Real K, Real T, Real r, Real q,
                                Real v0, Real kappa, Real theta, Real sigma_v, Real rho,
                                bool is_call, Size n_paths = 10000, uint64_t seed = 42);
};

inline AADGreeks AADGreeksEngine::bsm_european(Real S, Real K, Real T, Real r,
                                                Real q, Real sigma, bool is_call) {
    var vS = S, vK = K, vT = T, vr = r, vq = q, vsigma = sigma;
    var d1 = bsm_d1_var(vS, vK, vT, vr, vq, vsigma);
    var d2 = bsm_d2_var(d1, vsigma, vT);
    var Nd1 = normal_cdf_var(d1);
    var Nd2 = normal_cdf_var(d2);
    var price;
    if (is_call) {
        price = vS * exp(-vq * vT) * Nd1 - vK * exp(-vr * vT) * Nd2;
    } else {
        price = vK * exp(-vr * vT) * (var(1.0) - Nd2) - vS * exp(-vq * vT) * (var(1.0) - Nd1);
    }
    auto [dprice_dS, dprice_dK, dprice_dT, dprice_dr, dprice_dq, dprice_dsigma] =
        derivatives(price, wrt(vS, vK, vT, vr, vq, vsigma));

    Real d1_val = val(d1);
    Real gamma = normal_pdf(d1_val) * std::exp(-q * T) / (S * sigma * std::sqrt(T));

    return {
        .price = val(price),
        .delta = dprice_dS,
        .vega = dprice_dsigma,
        .rho = dprice_dr,
        .theta = -dprice_dT,
        .gamma = gamma
    };
}

inline AADGreeks AADGreeksEngine::heston_mc(Real S, Real K, Real T, Real r, Real q,
                                              Real v0, Real kappa, Real theta, Real sigma_v, Real rho,
                                              bool is_call, Size n_paths, uint64_t seed) {
    std::mt19937 gen(seed);
    std::normal_distribution<Real> norm(0.0, 1.0);

    Size n_steps = 1;
    Real dt = T;

    // 逐路径 AAD: 每条路径独立计算图 + Real 累加梯度。
    // 数学等价于整体 AAD (Leibniz: E[d/dθ Payoff] = d/dθ E[Payoff],路径独立时成立),
    // 但栈深度从 O(n_paths) 降到 O(path_length),避免 MSVC 1MB 栈溢出 (RISK-011)。
    // 注意: vv0/vS 必须在循环内声明,确保每路径计算图出作用域后被 shared_ptr 释放,
    // 不会跨路径累积节点。
    Real sum_price = 0.0;
    Real sum_delta = 0.0;
    Real sum_vega  = 0.0;  // d Price / d v0

    for (Size path = 0; path < n_paths; ++path) {
        var vS = S;
        var vv0 = v0;

        var logS = log(vS);
        var vsigma = sqrt(vv0);
        for (Size step = 0; step < n_steps; ++step) {
            Real Z = norm(gen);
            logS = logS + (r - q - vv0 / var(2.0)) * dt + vsigma * std::sqrt(dt) * Z;
        }
        var ST = exp(logS);
        var payoff;
        if (is_call) {
            payoff = max(ST - K, var(0.0));
        } else {
            payoff = max(K - ST, var(0.0));
        }
        var discounted = payoff * exp(-var(r) * var(T));

        auto [dS, dv0] = derivatives(discounted, wrt(vS, vv0));

        sum_price += val(discounted);
        sum_delta += dS;
        sum_vega  += dv0;
    }

    Real inv_n = 1.0 / static_cast<Real>(n_paths);
    return {
        .price = sum_price * inv_n,
        .delta = sum_delta * inv_n,
        .vega  = sum_vega  * inv_n,
        .rho = 0.0,
        .theta = 0.0,
        .gamma = 0.0
    };
}

}  // namespace cpphub::v1
