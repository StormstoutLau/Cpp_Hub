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
// Output buffer too small (JSON extension only, spec/mcp-server DESIGN §4.2).
// The required length is written back through `out_len`; caller reallocates
// and retries. Added after the original 4 codes — existing symbols never
// return it, so CPPHUB_ABI_VERSION stays 1 (additive change only).
#define CPPHUB_ERR_BUFFER_TOO_SMALL -4

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

// ---------------------------------------------------------------------------
// MCP tool registry extension (spec/mcp-server).
//
// These two symbols expose the diagnostic tool surface as JSON, so any
// language binding can act as an MCP server without re-implementing the
// protocol or the tool schemas.
//
// Buffer contract (NEW for these symbols; the original scalar symbols above
// have no buffer-query precedent):
//   * On entry, `*out_len` is the caller-allocated capacity of `out_json`.
//   * On return, `*out_len` is the number of bytes REQUIRED (excluding the
//     terminating NUL).
//   * If capacity is insufficient, returns CPPHUB_ERR_BUFFER_TOO_SMALL and
//     still writes the required length, so the caller can reallocate and
//     retry. No partial JSON is written in that case.
//   * On success, `out_json` holds a NUL-terminated JSON document and
//     `*out_len` is its length (excluding the NUL).
// ---------------------------------------------------------------------------

// Enumerate registered tools: a JSON array of
//   {"name": string, "description": string, "inputSchema": object}
// Returns CPPHUB_OK on success.
int cpphub_v1_tool_list(char* out_json, size_t* out_len);

// Invoke a tool by name. `args_json` is a JSON object of arguments.
// On success `out_json` receives a JSON object (the ToolResult container:
// tool / results[] / conclusion, see DESIGN §6.2).
// Unknown tool name -> CPPHUB_ERR_INVALID_ARG (message lists available tools).
// Malformed args_json or invalid argument values -> CPPHUB_ERR_INVALID_ARG.
// Kernel failure -> CPPHUB_ERR_RUNTIME.
int cpphub_v1_tool_call(const char* tool_name, const char* args_json,
                        char* out_json, size_t* out_len);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // CPPHUB_C_API_H
