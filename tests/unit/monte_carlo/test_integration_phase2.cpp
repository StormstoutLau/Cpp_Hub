#include <gtest/gtest.h>
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include "cpphub/pricing/pde/pde_engine.hpp"
#include "cpphub/pricing/tree/binomial.hpp"
#include "cpphub/pricing/tree/tree_engine.hpp"
#include "cpphub/models/diffusion/heston.hpp"
#include "cpphub/models/diffusion/heston_qe.hpp"
#include "cpphub/instruments/payoff/vanilla.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"
#include <cmath>
#include <complex>
#include <vector>

using namespace cpphub;

namespace {

Real bsm_call(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(S - K, 0.0);
    Real d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    Real d2 = d1 - sigma * std::sqrt(T);
    return S * std::exp(-q * T) * normal_cdf(d1) - K * std::exp(-r * T) * normal_cdf(d2);
}

Real generate_normal(Philox4x64& rng) {
    uint64_t r1 = rng();
    uint64_t r2 = rng();
    Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
    Real u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
    return box_muller(u1, u2).first;
}

Real heston_cf_price(Real S0, Real K, Real T, Real r, Real q, const HestonCFParams& p) {
    Complex i(0, 1);
    Real logK = std::log(K);

    auto phi = [&](Complex u) -> Complex {
        if (std::abs(u) < Real(1e-15)) return Complex(1, 0);
        return heston_characteristic_function(u, T, S0, p);
    };

    Real f_neg_i = S0 * std::exp((r - q) * T);
    Size N = 200000;
    Real du = 100.0 / static_cast<Real>(N);

    // P2 = 1/2 + 1/π * ∫₀^∞ Re[φ(u) * exp(-iu*ln(K)) / (iu)] du
    Real int_P2 = 0.0;
    for (Size j = 1; j < N; ++j) {
        Real u = static_cast<Real>(j) * du;
        Complex cu(u, 0);
        Complex phi_u = phi(cu);
        Complex num = phi_u * std::exp(-i * u * logK);
        Complex term = num / (Complex(0, 1) * u);
        Real val = std::real(term);
        if (std::isfinite(val)) int_P2 += val * du;
    }
    Real P2 = 0.5 + int_P2 / PI;

    // P1 = 1/2 + 1/π * ∫₀^∞ Re[φ(u-i) * exp(-iu*ln(K)) / (iu * φ(-i))] du
    Real int_P1 = 0.0;
    for (Size j = 1; j < N; ++j) {
        Real u = static_cast<Real>(j) * du;
        Complex cu(u, -1);
        Complex phi_ui = phi(cu);
        Complex num = phi_ui * std::exp(-i * u * logK);
        Complex term = num / (Complex(0, 1) * u * f_neg_i);
        Real val = std::real(term);
        if (std::isfinite(val)) int_P1 += val * du;
    }
    Real P1 = 0.5 + int_P1 / PI;

    return S0 * std::exp(-q * T) * P1 - K * std::exp(-r * T) * P2;
}

Real geom_asian_price(Real S0, Real K, Real T, Real r, Real q, Real sigma, Size n_steps) {
    if (T <= 0.0) return std::max(S0 - K, 0.0);
    Real m_geom = std::log(S0) + (r - q - 0.5 * sigma * sigma) * T
        * static_cast<Real>(n_steps + 1) / (2.0 * static_cast<Real>(n_steps));
    Real v_geom = sigma * sigma * T
        * static_cast<Real>(n_steps + 1) * static_cast<Real>(2 * n_steps + 1)
        / (6.0 * static_cast<Real>(n_steps * n_steps));
    Real d1 = (m_geom + v_geom - std::log(K)) / std::sqrt(v_geom);
    Real d2 = d1 - std::sqrt(v_geom);
    return std::exp(-r * T) * (std::exp(m_geom + 0.5 * v_geom) * normal_cdf(d1) - K * normal_cdf(d2));
}

// Halton sequence for QMC
class HaltonSeq {
public:
    explicit HaltonSeq(Size dim) : dim_(dim), count_(0) {
        primes_[0] = 2; primes_[1] = 3; primes_[2] = 5; primes_[3] = 7;
        primes_[4] = 11; primes_[5] = 13; primes_[6] = 17; primes_[7] = 19;
        primes_[8] = 23; primes_[9] = 29; primes_[10] = 31; primes_[11] = 37;
    }

    std::vector<Real> next() {
        std::vector<Real> result(dim_);
        uint64_t n = count_++;
        for (Size d = 0; d < dim_; ++d) {
            uint64_t base = primes_[d];
            uint64_t i = n + 1;
            Real f = 1.0 / static_cast<Real>(base);
            Real val = 0.0;
            uint64_t ival = i;
            while (ival > 0) {
                val += f * static_cast<Real>(ival % base);
                ival /= base;
                f /= static_cast<Real>(base);
            }
            result[d] = val;
        }
        return result;
    }

    Size dimension() const noexcept { return dim_; }

private:
    Size dim_;
    uint64_t count_;
    uint64_t primes_[12];
};

} // anonymous namespace

TEST(Integration, HestonMCMatchesCharacteristicFunction) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0;

    // Verify CF price against BSM for near-deterministic volatility
    // Use small sigma_vol so Heston approaches BSM; handle numerical limit
    HestonCFParams p_near_bsm{0.04, 0.5, 0.04, 0.01, 0.0, r, q};
    Real cf_near_bsm = heston_cf_price(S0, K, T, r, q, p_near_bsm);
    Real bsm_ref = bsm_call(S0, K, T, r, q, std::sqrt(0.04));
    EXPECT_TRUE(std::isfinite(cf_near_bsm));
    EXPECT_NEAR(cf_near_bsm, bsm_ref, bsm_ref * 0.05);

    HestonCFParams p{0.04, 2.0, 0.04, 0.5, -0.7, r, q};

    HestonParams hp{S0, p.v0, p.kappa, p.theta, p.sigma, p.rho, p.r, p.q};
    Heston h_heston(hp);

    Complex i(0,1);
    for (Real u_test = 0.5; u_test <= 10.0; u_test += 0.5) {
        Complex cu(u_test, 0);
        Complex cf1 = heston_characteristic_function(cu, T, S0, p);
        Complex cf2 = h_heston.characteristic_function(cu, T);
        EXPECT_NEAR(std::real(cf1), std::real(cf2), 1e-8);
        EXPECT_NEAR(std::imag(cf1), std::imag(cf2), 1e-8);
    }

    Real cf_price = heston_cf_price(S0, K, T, r, q, p);
    EXPECT_TRUE(std::isfinite(cf_price) && cf_price > 0.0);

    HestonQE heston_qe(hp);

    Size n_paths = 500000;
    Size n_steps = 1;
    Real sum_qe = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(static_cast<uint64_t>(j * 12345 + 1));
        std::vector<Real> path(n_steps + 1);
        heston_qe.generate_path(T, n_steps, path, rng);
        sum_qe += std::max(path.back() - K, 0.0);
    }
    Real qe_price = std::exp(-r * T) * sum_qe / static_cast<Real>(n_paths);

    {
        Heston h_euler(hp, HestonScheme::FullTruncation);
        n_steps = 52;
        Real sum_euler = 0.0;
        for (Size j = 0; j < n_paths; ++j) {
            Philox4x64 rng(static_cast<uint64_t>(j * 12345 + 1));
            std::vector<Real> path(n_steps + 1);
            h_euler.generate_path(T, n_steps, path, rng);
            sum_euler += std::max(path.back() - K, 0.0);
        }
        Real euler_price = std::exp(-r * T) * sum_euler / static_cast<Real>(n_paths);
        EXPECT_NEAR(euler_price, cf_price, cf_price * 0.03)
            << "euler_price=" << euler_price << " cf_price=" << cf_price;
    }

    EXPECT_NEAR(qe_price, cf_price, cf_price * 0.03)
        << "qe_price=" << qe_price << " cf_price=" << cf_price;
}

TEST(Integration, HestonQEVsEuler) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0;
    HestonParams hp{S0, 0.04, 1.5, 0.04, 0.3, -0.5, r, q};
    Heston heston_euler(hp, HestonScheme::FullTruncation);
    HestonQE heston_qe(hp);

    Size n_paths = 10000;
    Size n_steps = 52;

    auto compute_se = [&](auto& model) -> Real {
        Real sum_p = 0, sum_p2 = 0;
        for (Size j = 0; j < n_paths; ++j) {
            Philox4x64 rng(static_cast<uint64_t>(j * 9999 + 1));
            std::vector<Real> path(n_steps + 1);
            model.generate_path(T, n_steps, path, rng);
            Real ST = path.back();
            Real payoff = std::max(ST - K, 0.0);
            sum_p += payoff;
            sum_p2 += payoff * payoff;
        }
        Real mean = sum_p / static_cast<Real>(n_paths);
        Real var = sum_p2 / static_cast<Real>(n_paths) - mean * mean;
        return std::sqrt(std::max(var, 0.0) / static_cast<Real>(n_paths));
    };

    Real se_euler = compute_se(heston_euler);
    Real se_qe = compute_se(heston_qe);
    EXPECT_LT(se_qe, se_euler * 1.5);
}

TEST(Integration, PDEAmericanPut) {
    PDEEngineConfig cfg;
    cfg.n_spatial = 500;
    cfg.n_time = 2000;
    cfg.scheme = FDMSchemeType::CrankNicolson;
    cfg.alpha = 0.5;
    PDEEngine engine(cfg);

    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    PutPayOff payoff(K);
    Real price = engine.price_american(payoff, S0, K, T, r, q, sigma);
    EXPECT_NEAR(price, 6.0909, 1e-3);
}

TEST(Integration, PDEEuropeanMatchesBSM) {
    PDEEngineConfig cfg;
    cfg.n_spatial = 400;
    cfg.n_time = 1000;
    cfg.scheme = FDMSchemeType::CrankNicolson;
    PDEEngine engine(cfg);

    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0.02, sigma = 0.2;
    CallPayOff payoff(K);
    Real pde_price = engine.price_european(payoff, S0, K, T, r, q, sigma);
    Real bsm = bsm_call(S0, K, T, r, q, sigma);
    EXPECT_NEAR(pde_price, bsm, 1e-4);
}

TEST(Integration, TreeAmericanPut) {
    BinomialParams bp{100, 100, 0.25, 0.10, 0.0, 0.40, 201};
    BinomialTreeEngine bte(bp, BinomialType::LeisenReimer);
    PutPayOff payoff(100);
    Real tree_price = bte.price_american(payoff);

    PDEEngineConfig cfg;
    cfg.n_spatial = 500;
    cfg.n_time = 2000;
    cfg.scheme = FDMSchemeType::CrankNicolson;
    cfg.alpha = 0.5;
    PDEEngine engine(cfg);
    Real pde_price = engine.price_american(payoff, 100, 100, 0.25, 0.10, 0, 0.40);
    EXPECT_NEAR(tree_price, pde_price, 0.01);
}

TEST(Integration, QMCAsianOption) {
    Real S0 = 100, K = 100, T = 1, r = 0.05, q = 0, sigma = 0.2;
    Size n_steps = 12;
    Size n_paths = 1024;
    Real dt = T / static_cast<Real>(n_steps);

    // Reference: large-sample Halton + CV
    HaltonSeq halton_ref(n_steps);
    ControlVariate cv_ref(geom_asian_price(S0, K, T, r, q, sigma, n_steps));
    Size n_ref = 102400;
    for (Size j = 0; j < n_ref; ++j) {
        std::vector<Real> uniforms = halton_ref.next();
        Real S = S0;
        Real sum_arith = 0.0;
        Real prod_geom = 1.0;
        for (Size i = 0; i < n_steps; ++i) {
            Real u = uniforms[i];
            if (u < 1e-15) u = 1e-15;
            if (u > 1.0 - 1e-15) u = 1.0 - 1e-15;
            Real Z = inv_normal_cdf(u);
            S = S * std::exp((r - q - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * Z);
            sum_arith += S;
            prod_geom *= S;
        }
        Real avg_arith = sum_arith / static_cast<Real>(n_steps);
        Real avg_geom = std::pow(prod_geom, 1.0 / static_cast<Real>(n_steps));
        cv_ref.add_sample(std::max(avg_arith - K, 0.0), std::max(avg_geom - K, 0.0));
    }
    Real ref_price = std::exp(-r * T) * cv_ref.estimate();

    // Test: 1024 QMC + CV with Brownian bridge
    HaltonSeq halton(n_steps);
    ControlVariate cv(geom_asian_price(S0, K, T, r, q, sigma, n_steps));
    for (Size j = 0; j < n_paths; ++j) {
        std::vector<Real> uniforms = halton.next();
        Real S = S0;
        Real sum_arith = 0.0;
        Real prod_geom = 1.0;
        for (Size i = 0; i < n_steps; ++i) {
            Real u = uniforms[i];
            if (u < 1e-15) u = 1e-15;
            if (u > 1.0 - 1e-15) u = 1.0 - 1e-15;
            Real Z = inv_normal_cdf(u);
            S = S * std::exp((r - q - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * Z);
            sum_arith += S;
            prod_geom *= S;
        }
        Real avg_arith = sum_arith / static_cast<Real>(n_steps);
        Real avg_geom = std::pow(prod_geom, 1.0 / static_cast<Real>(n_steps));
        cv.add_sample(std::max(avg_arith - K, 0.0), std::max(avg_geom - K, 0.0));
    }

    Real cv_price = std::exp(-r * T) * cv.estimate();
    Real ratio = cv.variance_reduction_ratio();
    EXPECT_GT(ratio, 1.0);
    EXPECT_NEAR(cv_price, ref_price, 0.02);
}
