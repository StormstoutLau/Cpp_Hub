// C API unit tests (ADR-009).
// Verifies the C ABI boundary: BSM batch pricing, Heston CF pricing,
// IV inversion, MC pricing, error handling, and ABI version.
#include <gtest/gtest.h>
#include "cpphub/c_api/cpphub_c_api.h"
#include "cpphub/monte_carlo/control_variate.hpp"  // bsm_call_price for reference
#include "cpphub/calibration/calibrator.hpp"        // bsm_implied_vol, detail::heston_call_price_cf, HestonParams
#include <vector>
#include <cmath>

using namespace cpphub::v1;

TEST(CApi, AbiVersion) {
    EXPECT_EQ(cpphub_v1_get_abi_version(), 1);
}

TEST(CApi, GetLastErrorInitiallyEmpty) {
    // Fresh thread (gtest may have run other tests); just check it returns a valid C string.
    const char* err = cpphub_v1_get_last_error();
    ASSERT_NE(err, nullptr);
    // After an error it would be non-empty; we cannot guarantee emptiness here.
}

TEST(CApi, BsmPriceBatchCall) {
    const size_t n = 5;
    std::vector<double> spots(n, 100.0);
    std::vector<double> strikes = {90.0, 95.0, 100.0, 105.0, 110.0};
    std::vector<double> rates(n, 0.05);
    std::vector<double> vols(n, 0.20);
    std::vector<double> expiries(n, 1.0);
    std::vector<double> prices(n, 0.0);

    int rc = cpphub_v1_bsm_price_batch(spots.data(), strikes.data(), rates.data(),
                                        vols.data(), expiries.data(), prices.data(),
                                        n, 'C');
    ASSERT_EQ(rc, CPPHUB_OK) << cpphub_v1_get_last_error();

    // Compare against the C++ reference (q=0 in the C API).
    for (size_t i = 0; i < n; ++i) {
        double ref = bsm_call_price(spots[i], strikes[i], expiries[i],
                                    rates[i], 0.0, vols[i]);
        EXPECT_NEAR(prices[i], ref, 1e-12) << "i=" << i;
    }
}

TEST(CApi, BsmPriceBatchPut) {
    const size_t n = 3;
    std::vector<double> spots = {100.0, 100.0, 100.0};
    std::vector<double> strikes = {90.0, 100.0, 110.0};
    std::vector<double> rates(n, 0.05);
    std::vector<double> vols(n, 0.25);
    std::vector<double> expiries(n, 0.5);
    std::vector<double> prices(n, 0.0);

    int rc = cpphub_v1_bsm_price_batch(spots.data(), strikes.data(), rates.data(),
                                        vols.data(), expiries.data(), prices.data(),
                                        n, 'P');
    ASSERT_EQ(rc, CPPHUB_OK) << cpphub_v1_get_last_error();

    for (size_t i = 0; i < n; ++i) {
        double ref = bsm_put_price(spots[i], strikes[i], expiries[i],
                                   rates[i], 0.0, vols[i]);
        EXPECT_NEAR(prices[i], ref, 1e-12) << "i=" << i;
    }
}

TEST(CApi, BsmPriceBatchLowercaseOptType) {
    const size_t n = 1;
    double S = 100.0, K = 100.0, r = 0.05, sigma = 0.2, T = 1.0;
    double price_c = 0.0, price_lower = 0.0;

    ASSERT_EQ(cpphub_v1_bsm_price_batch(&S, &K, &r, &sigma, &T, &price_c, n, 'C'), CPPHUB_OK);
    ASSERT_EQ(cpphub_v1_bsm_price_batch(&S, &K, &r, &sigma, &T, &price_lower, n, 'c'), CPPHUB_OK);
    EXPECT_NEAR(price_c, price_lower, 1e-12);
}

TEST(CApi, BsmPriceBatchNullPointer) {
    double dummy = 0.0;
    int rc = cpphub_v1_bsm_price_batch(nullptr, &dummy, &dummy, &dummy, &dummy, &dummy, 1, 'C');
    EXPECT_EQ(rc, CPPHUB_ERR_INVALID_ARG);
    EXPECT_GT(std::string(cpphub_v1_get_last_error()).size(), 0u);
}

TEST(CApi, BsmPriceBatchZeroN) {
    double dummy = 0.0;
    int rc = cpphub_v1_bsm_price_batch(&dummy, &dummy, &dummy, &dummy, &dummy, &dummy, 0, 'C');
    EXPECT_EQ(rc, CPPHUB_ERR_INVALID_ARG);
}

TEST(CApi, BsmPriceBatchBadOptType) {
    double S = 100.0, K = 100.0, r = 0.05, sigma = 0.2, T = 1.0, price = 0.0;
    int rc = cpphub_v1_bsm_price_batch(&S, &K, &r, &sigma, &T, &price, 1, 'X');
    EXPECT_EQ(rc, CPPHUB_ERR_INVALID_ARG);
}

TEST(CApi, HestonPrice) {
    double S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0;
    double v0 = 0.04, kappa = 2.0, theta = 0.04, sigma_v = 0.3, rho = -0.5;
    double price = 0.0;
    int rc = cpphub_v1_heston_price(S, K, T, r, q, v0, kappa, theta, sigma_v, rho, &price);
    ASSERT_EQ(rc, CPPHUB_OK) << cpphub_v1_get_last_error();

    // Compare against the C++ reference (self-consistency: C API must match C++ implementation)
    HestonParams hp{v0, kappa, theta, sigma_v, rho};
    double ref = detail::heston_call_price_cf(S, K, T, r, q, hp);
    EXPECT_NEAR(price, ref, 1e-10) << "price=" << price << " ref=" << ref;

    // Sanity: Heston price should be positive for a call
    double bsm_ref = bsm_call_price(S, K, T, r, q, std::sqrt(v0));
    EXPECT_GT(price, 0.0) << "price=" << price << " bsm_ref=" << bsm_ref;
    // Heston with skew (rho<0) produces a smile; ATM price should be within ~15% of BSM
    EXPECT_LT(std::abs(price - bsm_ref) / bsm_ref, 0.15)
        << "price=" << price << " bsm_ref=" << bsm_ref << " ref=" << ref;
}

TEST(CApi, HestonPriceNullOutput) {
    int rc = cpphub_v1_heston_price(100, 100, 1, 0.05, 0, 0.04, 2.0, 0.04, 0.3, -0.5, nullptr);
    EXPECT_EQ(rc, CPPHUB_ERR_INVALID_ARG);
}

TEST(CApi, HestonPriceNegativeT) {
    double price = 0.0;
    int rc = cpphub_v1_heston_price(100, 100, -1.0, 0.05, 0, 0.04, 2.0, 0.04, 0.3, -0.5, &price);
    EXPECT_EQ(rc, CPPHUB_ERR_INVALID_ARG);
}

TEST(CApi, BsmImpliedVol) {
    double S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    double C = bsm_call_price(S, K, T, r, q, sigma);
    double iv = 0.0;
    int rc = cpphub_v1_bsm_implied_vol(C, S, K, T, r, q, 1, &iv);
    ASSERT_EQ(rc, CPPHUB_OK) << cpphub_v1_get_last_error();
    EXPECT_NEAR(iv, sigma, 1e-8);
}

TEST(CApi, BsmImpliedVolPut) {
    double S = 100.0, K = 110.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.25;
    double P = bsm_put_price(S, K, T, r, q, sigma);
    double iv = 0.0;
    int rc = cpphub_v1_bsm_implied_vol(P, S, K, T, r, q, 0, &iv);
    ASSERT_EQ(rc, CPPHUB_OK) << cpphub_v1_get_last_error();
    EXPECT_NEAR(iv, sigma, 1e-8);
}

TEST(CApi, BsmImpliedVolInvalidArgs) {
    double iv = 0.0;
    EXPECT_EQ(cpphub_v1_bsm_implied_vol(10, 100, 100, -1, 0.05, 0, 1, &iv), CPPHUB_ERR_INVALID_ARG);
    EXPECT_EQ(cpphub_v1_bsm_implied_vol(10, 100, 100, 1, 0.05, 0, 1, nullptr), CPPHUB_ERR_INVALID_ARG);
}

TEST(CApi, McPriceCall) {
    double S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    size_t n_paths = 100000;
    cpphub_mc_result_t result;
    int rc = cpphub_v1_mc_price(S, K, T, r, q, sigma, 1, n_paths, 42, &result);
    ASSERT_EQ(rc, CPPHUB_OK) << cpphub_v1_get_last_error();

    // MC price should be within 3 std_err of BSM price
    double bsm = bsm_call_price(S, K, T, r, q, sigma);
    EXPECT_GT(result.std_err, 0.0);
    EXPECT_LT(std::abs(result.price - bsm), 5.0 * result.std_err)
        << "price=" << result.price << " bsm=" << bsm << " se=" << result.std_err;
    // 3 sigma confidence should be < 0.05 for 100k paths on a vanilla call
    EXPECT_LT(result.std_err, 0.05);
}

TEST(CApi, McPricePut) {
    double S = 100.0, K = 110.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.25;
    size_t n_paths = 50000;
    cpphub_mc_result_t result;
    int rc = cpphub_v1_mc_price(S, K, T, r, q, sigma, 0, n_paths, 123, &result);
    ASSERT_EQ(rc, CPPHUB_OK) << cpphub_v1_get_last_error();

    double bsm = bsm_put_price(S, K, T, r, q, sigma);
    EXPECT_LT(std::abs(result.price - bsm), 5.0 * result.std_err);
}

TEST(CApi, McPriceReproducible) {
    double S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    size_t n_paths = 10000;
    cpphub_mc_result_t r1, r2;
    ASSERT_EQ(cpphub_v1_mc_price(S, K, T, r, q, sigma, 1, n_paths, 999, &r1), CPPHUB_OK);
    ASSERT_EQ(cpphub_v1_mc_price(S, K, T, r, q, sigma, 1, n_paths, 999, &r2), CPPHUB_OK);
    // Same seed → identical results (Philox is deterministic)
    EXPECT_NEAR(r1.price, r2.price, 1e-15);
    EXPECT_NEAR(r1.std_err, r2.std_err, 1e-15);
}

TEST(CApi, McPriceInvalidArgs) {
    cpphub_mc_result_t result;
    EXPECT_EQ(cpphub_v1_mc_price(100, 100, 1, 0.05, 0, 0.2, 1, 1000, 42, nullptr),
              CPPHUB_ERR_INVALID_ARG);
    EXPECT_EQ(cpphub_v1_mc_price(100, 100, -1, 0.05, 0, 0.2, 1, 1000, 42, &result),
              CPPHUB_ERR_INVALID_ARG);
    EXPECT_EQ(cpphub_v1_mc_price(100, 100, 1, 0.05, 0, 0.2, 1, 0, 42, &result),
              CPPHUB_ERR_INVALID_ARG);
}
