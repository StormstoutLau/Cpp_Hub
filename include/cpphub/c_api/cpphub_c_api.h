// CppHub C ABI - Stable C interface for cross-language bindings.
// ADR-009: versioned symbols `cpphub_v1_*`, POD types only, error codes.
// Implemented in src/c_api/c_api.cpp.
//
// Usage from C:
//   #include "cpphub/c_api/cpphub_c_api.h"
//   double price = 0;
//   int rc = cpphub_v1_heston_price(100, 100, 1, 0.05, 0,
//                                    0.04, 2.0, 0.04, 0.3, -0.5, &price);
//   if (rc != 0) fprintf(stderr, "error: %s\n", cpphub_v1_get_last_error());
//
// Usage from Python (ctypes):
//   lib = ctypes.CDLL("libcpphub_c_api.so")
//   lib.cpphub_v1_get_abi_version.restype = ctypes.c_int
//   assert lib.cpphub_v1_get_abi_version() == 1
#ifndef CPPHUB_C_API_H
#define CPPHUB_C_API_H

#include <stddef.h>  // size_t

#ifdef __cplusplus
extern "C" {
#endif

#define CPPHUB_ABI_VERSION 1

// Error codes returned by all cpphub_v1_* functions.
// 0 = success; negative = failure (use cpphub_v1_get_last_error for details).
#define CPPHUB_OK                  0
#define CPPHUB_ERR_INVALID_ARG    -1
#define CPPHUB_ERR_RUNTIME        -2
#define CPPHUB_ERR_UNKNOWN        -3

// Greeks struct (POD, C-compatible).
typedef struct {
    double delta;
    double gamma;
    double vega;
    double theta;
    double rho;
} cpphub_greeks_t;

// Monte Carlo result (POD, C-compatible).
typedef struct {
    double price;
    double std_err;
} cpphub_mc_result_t;

// ---------------------------------------------------------------------------
// Batch BSM option pricing.
// All input arrays must have at least `n` elements; `prices` must have `n`
// slots for output. `opt_type` is 'C' for call, 'P' for put (case-insensitive).
// Returns CPPHUB_OK on success, error code otherwise.
// ---------------------------------------------------------------------------
int cpphub_v1_bsm_price_batch(const double* spots,
                              const double* strikes,
                              const double* rates,
                              const double* vols,
                              const double* expiries,
                              double* prices,
                              size_t n,
                              char opt_type);

// ---------------------------------------------------------------------------
// Single Heston call price via Carr-Madan Fourier inversion.
// params: v0, kappa, theta, sigma_v, rho (Heston 1993 parameters).
// `price` receives the call price on success.
// ---------------------------------------------------------------------------
int cpphub_v1_heston_price(double S, double K, double T, double r, double q,
                           double v0, double kappa, double theta,
                           double sigma_v, double rho, double* price);

// ---------------------------------------------------------------------------
// BSM implied volatility inversion (Newton-Raphson with vega).
// `is_call` = 1 for call, 0 for put.
// `iv` receives the implied vol on success.
// ---------------------------------------------------------------------------
int cpphub_v1_bsm_implied_vol(double C_market, double S, double K, double T,
                              double r, double q, int is_call, double* iv);

// ---------------------------------------------------------------------------
// Monte Carlo option price under geometric Brownian motion (BSM dynamics).
// Uses Philox4x64-10 RNG for reproducibility. `seed` initializes the RNG key.
// `result` receives {price, std_err} on success.
// ---------------------------------------------------------------------------
int cpphub_v1_mc_price(double S, double K, double T, double r, double q,
                       double sigma, int is_call, size_t n_paths,
                       unsigned long long seed, cpphub_mc_result_t* result);

// ---------------------------------------------------------------------------
// Retrieve the last error message for the calling thread.
// Returns a pointer to a thread-local NUL-terminated string. The pointer is
// valid until the next cpphub_v1_* call from the same thread. Returns "" if
// no error has occurred.
// ---------------------------------------------------------------------------
const char* cpphub_v1_get_last_error(void);

// Returns the ABI version (currently 1).
int cpphub_v1_get_abi_version(void);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // CPPHUB_C_API_H
