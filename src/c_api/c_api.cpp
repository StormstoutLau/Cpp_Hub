// CppHub C ABI implementation.
// ADR-009: All C++ exceptions are caught at the C boundary and converted to
// error codes. The last error message is stored in thread_local storage.
#include "cpphub/c_api/cpphub_c_api.h"

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"  // bsm_call_price / bsm_put_price
#include "cpphub/calibration/calibrator.hpp"        // bsm_implied_vol, detail::heston_call_price_cf, HestonParams

#include <cmath>
#include <stdexcept>
#include <string>
#include <cstring>

namespace {

// Thread-local error message buffer.
thread_local std::string g_last_error;

void set_error(const char* msg) {
    g_last_error = msg;
}

void set_error(const std::string& msg) {
    g_last_error = msg;
}

int translate_exception() {
    try {
        throw;
    } catch (const std::invalid_argument& e) {
        set_error(e.what());
        return CPPHUB_ERR_INVALID_ARG;
    } catch (const std::runtime_error& e) {
        set_error(e.what());
        return CPPHUB_ERR_RUNTIME;
    } catch (const std::exception& e) {
        set_error(e.what());
        return CPPHUB_ERR_UNKNOWN;
    } catch (...) {
        set_error("unknown C++ exception");
        return CPPHUB_ERR_UNKNOWN;
    }
}

}  // namespace

extern "C" {

int cpphub_v1_get_abi_version(void) {
    return CPPHUB_ABI_VERSION;
}

const char* cpphub_v1_get_last_error(void) {
    if (g_last_error.empty()) return "";
    return g_last_error.c_str();
}

int cpphub_v1_bsm_price_batch(const double* spots,
                              const double* strikes,
                              const double* rates,
                              const double* vols,
                              const double* expiries,
                              double* prices,
                              size_t n,
                              char opt_type) {
    try {
        if (spots == nullptr || strikes == nullptr || rates == nullptr ||
            vols == nullptr || expiries == nullptr || prices == nullptr) {
            set_error("null pointer argument");
            return CPPHUB_ERR_INVALID_ARG;
        }
        if (n == 0) {
            set_error("n must be positive");
            return CPPHUB_ERR_INVALID_ARG;
        }
        char t = opt_type;
        if (t == 'c') t = 'C';
        if (t == 'p') t = 'P';
        if (t != 'C' && t != 'P') {
            set_error("opt_type must be 'C' or 'P'");
            return CPPHUB_ERR_INVALID_ARG;
        }
        using cpphub::v1::Real;
        using cpphub::v1::bsm_call_price;
        using cpphub::v1::bsm_put_price;
        for (size_t i = 0; i < n; ++i) {
            Real S = spots[i];
            Real K = strikes[i];
            Real r = rates[i];
            Real sigma = vols[i];
            Real T = expiries[i];
            Real q = 0.0;  // batch API does not take dividend; caller discounts manually
            prices[i] = (t == 'C') ? bsm_call_price(S, K, T, r, q, sigma)
                                   : bsm_put_price(S, K, T, r, q, sigma);
        }
        return CPPHUB_OK;
    } catch (...) {
        return translate_exception();
    }
}

int cpphub_v1_heston_price(double S, double K, double T, double r, double q,
                           double v0, double kappa, double theta,
                           double sigma_v, double rho, double* price) {
    try {
        if (price == nullptr) {
            set_error("null output pointer");
            return CPPHUB_ERR_INVALID_ARG;
        }
        if (T <= 0.0) {
            set_error("T must be positive");
            return CPPHUB_ERR_INVALID_ARG;
        }
        if (S <= 0.0 || K <= 0.0) {
            set_error("S and K must be positive");
            return CPPHUB_ERR_INVALID_ARG;
        }
        using cpphub::v1::Real;
        using cpphub::v1::HestonParams;
        using cpphub::v1::detail::heston_call_price_cf;
        HestonParams hp{v0, kappa, theta, sigma_v, rho};
        *price = heston_call_price_cf(S, K, T, r, q, hp);
        return CPPHUB_OK;
    } catch (...) {
        return translate_exception();
    }
}

int cpphub_v1_bsm_implied_vol(double C_market, double S, double K, double T,
                              double r, double q, int is_call, double* iv) {
    try {
        if (iv == nullptr) {
            set_error("null output pointer");
            return CPPHUB_ERR_INVALID_ARG;
        }
        if (T <= 0.0 || C_market <= 0.0 || S <= 0.0 || K <= 0.0) {
            set_error("T, C_market, S, K must be positive");
            return CPPHUB_ERR_INVALID_ARG;
        }
        using cpphub::v1::bsm_implied_vol;
        *iv = bsm_implied_vol(C_market, S, K, T, r, q, is_call != 0);
        return CPPHUB_OK;
    } catch (...) {
        return translate_exception();
    }
}

int cpphub_v1_mc_price(double S, double K, double T, double r, double q,
                       double sigma, int is_call, size_t n_paths,
                       unsigned long long seed, cpphub_mc_result_t* result) {
    try {
        if (result == nullptr) {
            set_error("null output pointer");
            return CPPHUB_ERR_INVALID_ARG;
        }
        if (T <= 0.0 || S <= 0.0 || K <= 0.0 || sigma <= 0.0) {
            set_error("T, S, K, sigma must be positive");
            return CPPHUB_ERR_INVALID_ARG;
        }
        if (n_paths == 0) {
            set_error("n_paths must be positive");
            return CPPHUB_ERR_INVALID_ARG;
        }
        using cpphub::v1::Real;
        using cpphub::v1::Philox4x64;
        using cpphub::v1::box_muller;
        using cpphub::v1::Size;

        Philox4x64 rng(static_cast<uint64_t>(seed));
        Real drift = (r - q - 0.5 * sigma * sigma) * T;
        Real vol = sigma * std::sqrt(T);
        Real disc = std::exp(-r * T);

        // Welford-style accumulation for mean and variance of payoff.
        double sum = 0.0;
        double sum_sq = 0.0;
        Size n = static_cast<Size>(n_paths);
        Size i = 0;
        while (i < n) {
            // Generate up to 2 normals per RNG draw (Box-Muller pair).
            uint64_t r1 = rng();
            uint64_t r2 = rng();
            double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            // Guard u1 == 0 to avoid log(0)
            if (u1 < 1e-300) u1 = 1e-300;
            auto [z1, z2] = box_muller(u1, u2);

            Real S_T1 = S * std::exp(drift + vol * z1);
            Real p1 = is_call ? std::max(S_T1 - K, 0.0) : std::max(K - S_T1, 0.0);
            double pv1 = disc * p1;
            sum += pv1;
            sum_sq += pv1 * pv1;
            ++i;

            if (i < n) {
                Real S_T2 = S * std::exp(drift + vol * z2);
                Real p2 = is_call ? std::max(S_T2 - K, 0.0) : std::max(K - S_T2, 0.0);
                double pv2 = disc * p2;
                sum += pv2;
                sum_sq += pv2 * pv2;
                ++i;
            }
        }

        double mean = sum / static_cast<double>(n);
        double var = (n > 1)
            ? (sum_sq - static_cast<double>(n) * mean * mean) / static_cast<double>(n - 1)
            : 0.0;
        if (var < 0.0) var = 0.0;
        double std_err = std::sqrt(var / static_cast<double>(n));

        result->price = mean;
        result->std_err = std_err;
        return CPPHUB_OK;
    } catch (...) {
        return translate_exception();
    }
}

}  // extern "C"
