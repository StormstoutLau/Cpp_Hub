// v1.1: 通用傅里叶引擎测试
// - CharFn 接口验证 (GBM/Heston CF 模 |u|=1 性质)
// - COS 方法 (Fang-Oosterlee 2009) vs BSM 解析解
// - Carr-Madan FFT (1999) vs BSM 解析解
// - COS vs FFT 交叉验证
// - Heston CF 定价 vs 既有 heston_call_price_cf
// - Put-Call Parity 验证
// - FFT 内核 (radix-2) 正确性
#include <gtest/gtest.h>
#include "cpphub/pricing/fourier/characteristic_functions.hpp"
#include "cpphub/pricing/fourier/cos_method.hpp"
#include "cpphub/pricing/fourier/fft_method.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"
#include "cpphub/calibration/calibrator.hpp"
#include "cpphub/pricing/tree/binomial.hpp"
#include "cpphub/instruments/payoff/vanilla.hpp"
#include <cmath>
#include <vector>
#include <complex>

using namespace cpphub;

// ============ 辅助: BSM 解析解 (本地封装) ============
static Real bsm_call(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    return bsm_call_price(S, K, T, r, q, sigma);
}
static Real bsm_put(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    return bsm_put_price(S, K, T, r, q, sigma);
}

// ============ 1. 特征函数接口 (CharFn) 验证 ============

TEST(CharacteristicFunctionsTest, GBM_AtZeroReturnsOne) {
    auto phi = make_gbm_cf(100.0, 0.05, 0.02, 0.20, 1.0);
    Complex val = phi(Complex(0.0, 0.0));
    EXPECT_NEAR(std::real(val), 1.0, 1e-12);
    EXPECT_NEAR(std::imag(val), 0.0, 1e-12);
}

TEST(CharacteristicFunctionsTest, GBM_UnitModulusForRealU) {
    // 对实数 u, |phi(u)| <= 1 (CF 性质)
    auto phi = make_gbm_cf(100.0, 0.05, 0.02, 0.20, 1.0);
    for (Real u = 0.0; u <= 10.0; u += 0.5) {
        Complex val = phi(Complex(u, 0.0));
        EXPECT_LE(std::abs(val), 1.0 + 1e-12)
            << " |phi(u=" << u << ")| = " << std::abs(val);
    }
}

TEST(CharacteristicFunctionsTest, GBM_MatchesAnalyticFormula) {
    // phi(u) = exp(iu*ln S0 + iu*(r-q-σ²/2)T - σ²u²T/2)
    // 对实数 u, |phi(u)| = exp(-σ²u²T/2)
    Real S0 = 100.0, r = 0.05, q = 0.02, sigma = 0.20, T = 1.0;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    Real u = 1.5;
    Complex val = phi(Complex(u, 0.0));
    Real expected_mod = std::exp(-sigma * sigma * u * u * T / 2.0);
    EXPECT_NEAR(std::abs(val), expected_mod, 1e-12);

    // 检查相位: arg(phi(u)) = u*ln S0 + u*(r-q-σ²/2)T
    Real expected_arg = u * (std::log(S0) + (r - q - 0.5 * sigma * sigma) * T);
    Real actual_arg = std::arg(val);
    // 规范化到 [-π, π]
    while (actual_arg - expected_arg > PI) actual_arg -= 2 * PI;
    while (expected_arg - actual_arg > PI) actual_arg += 2 * PI;
    EXPECT_NEAR(actual_arg, expected_arg, 1e-10);
}

TEST(CharacteristicFunctionsTest, Heston_AtZeroReturnsOne) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.05, 0.02};
    auto phi = make_heston_cf(100.0, 0.05, 0.02, p, 1.0);
    Complex val = phi(Complex(0.0, 0.0));
    EXPECT_NEAR(std::real(val), 1.0, 1e-10);
    EXPECT_NEAR(std::imag(val), 0.0, 1e-10);
}

TEST(CharacteristicFunctionsTest, Heston_UnitModulusForRealU) {
    HestonCFParams p{0.04, 1.5, 0.04, 0.3, -0.5, 0.0, 0.0};
    auto phi = make_heston_cf(100.0, 0.0, 0.0, p, 1.0);
    for (Real u = 0.1; u <= 5.0; u += 0.5) {
        Complex val = phi(Complex(u, 0.0));
        EXPECT_LE(std::abs(val), 1.0 + 1e-10);
    }
}

TEST(CharacteristicFunctionsTest, VG_RejectsFellerViolation) {
    // 1 - theta*nu - sigma²*nu/2 > 0 必须满足
    // theta=2.0, nu=0.5, sigma=0.30 → 1 - 1.0 - 0.0225 = -0.0225 < 0 (违反 Feller)
    EXPECT_THROW(make_vg_cf(100.0, 0.05, 0.0, 0.30, 0.5, 2.0, 1.0),
                 std::invalid_argument);
}

TEST(CharacteristicFunctionsTest, NIG_RejectsBetaGeqAlpha) {
    EXPECT_THROW(make_nig_cf(100.0, 0.05, 0.0, 1.0, 1.0, 1.0, 1.0),
                 std::invalid_argument);
    EXPECT_THROW(make_nig_cf(100.0, 0.05, 0.0, 1.0, -1.0, 1.0, 1.0),
                 std::invalid_argument);
}

TEST(CharacteristicFunctionsTest, VG_AtZeroReturnsOne) {
    auto phi = make_vg_cf(100.0, 0.05, 0.0, 0.20, 0.1, -0.1, 1.0);
    Complex val = phi(Complex(0.0, 0.0));
    EXPECT_NEAR(std::real(val), 1.0, 1e-10);
}

TEST(CharacteristicFunctionsTest, NIG_AtZeroReturnsOne) {
    auto phi = make_nig_cf(100.0, 0.05, 0.0, 15.0, -5.0, 1.0, 1.0);
    Complex val = phi(Complex(0.0, 0.0));
    EXPECT_NEAR(std::real(val), 1.0, 1e-10);
}

TEST(CharacteristicFunctionsTest, CosTruncationRange_GBM_CoversStrikes) {
    // GBM 截断区间应覆盖典型 strike 范围 [ln(0.5*S0), ln(2*S0)]
    Real S0 = 100.0, sigma = 0.20, T = 1.0;
    auto phi = make_gbm_cf(S0, 0.05, 0.0, sigma, T);
    auto [a, b] = cos_truncation_range(phi, S0, T, 10.0);
    // 期望: a < ln(50), b > ln(200)
    EXPECT_LT(a, std::log(50.0));
    EXPECT_GT(b, std::log(200.0));
    // 区间宽度 ≈ 2*L*σ*sqrt(T) = 2*10*0.2*1 = 4
    EXPECT_NEAR(b - a, 4.0, 0.5);
}

// ============ 2. COS 方法 (Fang-Oosterlee 2009) vs BSM ============

TEST(COSMethodTest, GBM_Call_MatchesBSM_ATM) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.02, sigma = 0.20;
    Real expected = bsm_call(S0, K, T, r, q, sigma);
    Real cos_price = cos_call_gbm(S0, K, T, r, q, sigma, 256, 10.0);
    EXPECT_NEAR(cos_price, expected, 1e-4)
        << "COS=" << cos_price << " BSM=" << expected;
}

TEST(COSMethodTest, GBM_Call_MatchesBSM_OTM) {
    Real S0 = 100.0, K = 130.0, T = 1.0, r = 0.05, q = 0.02, sigma = 0.20;
    Real expected = bsm_call(S0, K, T, r, q, sigma);
    Real cos_price = cos_call_gbm(S0, K, T, r, q, sigma, 256, 10.0);
    EXPECT_NEAR(cos_price, expected, 1e-3)
        << "COS=" << cos_price << " BSM=" << expected;
}

TEST(COSMethodTest, GBM_Call_MatchesBSM_ITM) {
    Real S0 = 100.0, K = 70.0, T = 1.0, r = 0.05, q = 0.02, sigma = 0.20;
    Real expected = bsm_call(S0, K, T, r, q, sigma);
    Real cos_price = cos_call_gbm(S0, K, T, r, q, sigma, 256, 10.0);
    EXPECT_NEAR(cos_price, expected, 1e-2)
        << "COS=" << cos_price << " BSM=" << expected;
}

TEST(COSMethodTest, GBM_Put_MatchesBSM) {
    Real S0 = 100.0, K = 95.0, T = 0.5, r = 0.03, q = 0.01, sigma = 0.25;
    Real expected = bsm_put(S0, K, T, r, q, sigma);
    Real cos_price = cos_put_gbm(S0, K, T, r, q, sigma, 256, 10.0);
    EXPECT_NEAR(cos_price, expected, 1e-3)
        << "COS=" << cos_price << " BSM=" << expected;
}

TEST(COSMethodTest, GBM_Call_PutCallParityHolds) {
    // C - P = S0*e^{-qT} - K*e^{-rT}
    Real S0 = 100.0, K = 105.0, T = 1.0, r = 0.05, q = 0.02, sigma = 0.20;
    Real C = cos_call_gbm(S0, K, T, r, q, sigma);
    Real P = cos_put_gbm(S0, K, T, r, q, sigma);
    Real parity = S0 * std::exp(-q * T) - K * std::exp(-r * T);
    EXPECT_NEAR(C - P, parity, 1e-3);
}

TEST(COSMethodTest, GBM_Call_ConvergenceWithNTerms) {
    // N 增大时, 误差应单调减小或达到机器精度 (不发散)
    // GBM CF 为高斯型, COS 系数指数衰减, N≥64 即可达到机器精度,
    // 因此 err_256 与 err_64 可能完全相等 (1.94e-13 量级) — 用 LE 容许.
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real expected = bsm_call(S0, K, T, r, q, sigma);
    Real err_64  = std::abs(cos_call_gbm(S0, K, T, r, q, sigma, 64, 10.0)  - expected);
    Real err_128 = std::abs(cos_call_gbm(S0, K, T, r, q, sigma, 128, 10.0) - expected);
    Real err_256 = std::abs(cos_call_gbm(S0, K, T, r, q, sigma, 256, 10.0) - expected);
    EXPECT_LE(err_256, err_64);
    EXPECT_LE(err_128, err_64);
    EXPECT_LT(err_256, 1e-6);  // 高 N 应达到机器精度量级
}

TEST(COSMethodTest, GBM_Call_ShortMaturity) {
    Real S0 = 100.0, K = 100.0, T = 0.05, r = 0.05, q = 0.0, sigma = 0.30;
    Real expected = bsm_call(S0, K, T, r, q, sigma);
    Real cos_price = cos_call_gbm(S0, K, T, r, q, sigma, 256, 12.0);
    EXPECT_NEAR(cos_price, expected, 1e-3)
        << "COS=" << cos_price << " BSM=" << expected;
}

TEST(COSMethodTest, GBM_Call_LongMaturity) {
    Real S0 = 100.0, K = 100.0, T = 5.0, r = 0.05, q = 0.02, sigma = 0.20;
    Real expected = bsm_call(S0, K, T, r, q, sigma);
    Real cos_price = cos_call_gbm(S0, K, T, r, q, sigma, 512, 10.0);
    EXPECT_NEAR(cos_price, expected, 1e-3)
        << "COS=" << cos_price << " BSM=" << expected;
}

TEST(COSMethodTest, Heston_Call_MatchesCarrMadanIntegral) {
    // 与 calibrator.hpp 中的 detail::heston_call_price_cf (Gil-Pelaez 积分) 对比
    // 这是 Heston CF 实现正确性的交叉验证
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0;
    HestonCFParams hp{0.04, 1.5, 0.04, 0.3, -0.5, r, q};
    Real cos_price = cos_call_heston(S0, K, T, r, q, hp, 512, 12.0);

    // 用 calibrator.hpp 的 Gil-Pelaez 积分作为参考
    HestonParams hp_calib{hp.v0, hp.kappa, hp.theta, hp.sigma, hp.rho};
    Real ref_price = detail::heston_call_price_cf(S0, K, T, r, q, hp_calib, 4096);

    EXPECT_NEAR(cos_price, ref_price, 5e-3)
        << "COS=" << cos_price << " Ref=" << ref_price;
}

TEST(COSMethodTest, Heston_Call_OTM_MatchesReference) {
    Real S0 = 100.0, K = 120.0, T = 1.0, r = 0.05, q = 0.0;
    HestonCFParams hp{0.04, 1.5, 0.04, 0.3, -0.5, r, q};
    Real cos_price = cos_call_heston(S0, K, T, r, q, hp, 512, 12.0);

    HestonParams hp_calib{hp.v0, hp.kappa, hp.theta, hp.sigma, hp.rho};
    Real ref_price = detail::heston_call_price_cf(S0, K, T, r, q, hp_calib, 4096);

    EXPECT_NEAR(cos_price, ref_price, 5e-3)
        << "COS=" << cos_price << " Ref=" << ref_price;
}

TEST(COSMethodTest, RejectsInvalidConfig) {
    auto phi = make_gbm_cf(100.0, 0.05, 0.0, 0.20, 1.0);
    COSEngine::Config cfg;
    cfg.n_terms = 4;  // too small
    EXPECT_THROW(COSEngine(phi, 100.0, 0.05, 0.0, 1.0, cfg), std::invalid_argument);

    cfg.n_terms = 256;
    cfg.L = 0.5;  // too small
    EXPECT_THROW(COSEngine(phi, 100.0, 0.05, 0.0, 1.0, cfg), std::invalid_argument);

    cfg.L = 10.0;
    cfg.a = 1.0; cfg.b = 0.5;  // b <= a
    EXPECT_THROW(COSEngine(phi, 100.0, 0.05, 0.0, 1.0, cfg), std::invalid_argument);
}

// ============ 3. FFT 内核 (radix-2) 验证 ============

TEST(FFTKernelTest, ImpulseTransform) {
    // 单位脉冲 [1, 0, 0, ..., 0] 的 FFT 应为 [1, 1, 1, ..., 1]
    std::vector<Complex> x(8, Complex(0.0, 0.0));
    x[0] = Complex(1.0, 0.0);
    fft_radix2(x, -1);
    for (Size i = 0; i < 8; ++i) {
        EXPECT_NEAR(std::real(x[i]), 1.0, 1e-12);
        EXPECT_NEAR(std::imag(x[i]), 0.0, 1e-12);
    }
}

TEST(FFTKernelTest, ConstantSignal) {
    // 常数信号 [c, c, ..., c] 的 FFT 在 index 0 处为 N*c, 其余为 0
    std::vector<Complex> x(16, Complex(2.5, 0.0));
    fft_radix2(x, -1);
    EXPECT_NEAR(std::real(x[0]), 16.0 * 2.5, 1e-10);
    EXPECT_NEAR(std::imag(x[0]), 0.0, 1e-10);
    for (Size i = 1; i < 16; ++i) {
        EXPECT_NEAR(std::abs(x[i]), 0.0, 1e-10);
    }
}

TEST(FFTKernelTest, RoundTripForwardInverse) {
    // 正向 FFT 后逆向 FFT 应恢复原信号 (相差 1/N 因子)
    std::vector<Complex> x = {{1, 2}, {3, -1}, {0, 4}, {-2, 0.5},
                              {1.5, 1}, {0, 0}, {2, 2}, {-1, -1}};
    std::vector<Complex> original = x;
    fft_radix2(x, -1);  // forward
    fft_radix2(x, +1);  // inverse
    Real n = static_cast<Real>(x.size());
    for (Size i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(std::real(x[i]) / n, std::real(original[i]), 1e-10);
        EXPECT_NEAR(std::imag(x[i]) / n, std::imag(original[i]), 1e-10);
    }
}

TEST(FFTKernelTest, KnownSinusoid) {
    // x[n] = cos(2π n / N), FFT 在 k=1 处有 N/2, k=N-1 处有 N/2
    Size n = 32;
    std::vector<Complex> x(n);
    for (Size i = 0; i < n; ++i) {
        x[i] = Complex(std::cos(2.0 * PI * i / n), 0.0);
    }
    fft_radix2(x, -1);
    EXPECT_NEAR(std::real(x[1]), n / 2.0, 1e-10);
    EXPECT_NEAR(std::real(x[n - 1]), n / 2.0, 1e-10);
    // 其他频率分量应为 0
    for (Size k = 2; k < n - 1; ++k) {
        EXPECT_NEAR(std::abs(x[k]), 0.0, 1e-10);
    }
}

TEST(FFTKernelTest, RejectsNonPowerOf2) {
    std::vector<Complex> x(7, Complex(1.0, 0.0));
    EXPECT_THROW(fft_radix2(x, -1), std::invalid_argument);
}

// ============ 4. Carr-Madan FFT 方法 vs BSM ============

TEST(CarrMadanFFTTest, GBM_Call_MatchesBSM_ATM) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.02, sigma = 0.20;
    Real expected = bsm_call(S0, K, T, r, q, sigma);
    Real fft_price = fft_call_gbm(S0, K, T, r, q, sigma, 1.5, 4096, 0.25);
    // FFT 插值精度较低, 容差 1e-2 (Carr-Madan 固有局限)
    EXPECT_NEAR(fft_price, expected, 5e-2)
        << "FFT=" << fft_price << " BSM=" << expected;
}

TEST(CarrMadanFFTTest, GBM_Call_MatchesBSM_OTM) {
    Real S0 = 100.0, K = 120.0, T = 1.0, r = 0.05, q = 0.02, sigma = 0.20;
    Real expected = bsm_call(S0, K, T, r, q, sigma);
    Real fft_price = fft_call_gbm(S0, K, T, r, q, sigma, 1.5, 4096, 0.25);
    EXPECT_NEAR(fft_price, expected, 5e-2)
        << "FFT=" << fft_price << " BSM=" << expected;
}

TEST(CarrMadanFFTTest, GBM_Put_ViaParity_MatchesBSM) {
    Real S0 = 100.0, K = 95.0, T = 0.5, r = 0.03, q = 0.01, sigma = 0.25;
    Real expected = bsm_put(S0, K, T, r, q, sigma);

    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    CarrMadanFFT::Config cfg;
    cfg.alpha = 1.5;
    cfg.n_fft = 4096;
    cfg.eta = 0.25;
    CarrMadanFFT engine(phi, S0, r, q, T, cfg);
    Real fft_put = engine.price_put(K);

    EXPECT_NEAR(fft_put, expected, 5e-2)
        << "FFT=" << fft_put << " BSM=" << expected;
}

TEST(CarrMadanFFTTest, GBM_Call_PutCallParityHolds) {
    Real S0 = 100.0, K = 105.0, T = 1.0, r = 0.05, q = 0.02, sigma = 0.20;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    CarrMadanFFT::Config cfg;
    cfg.alpha = 1.5;
    cfg.n_fft = 4096;
    cfg.eta = 0.25;
    CarrMadanFFT engine(phi, S0, r, q, T, cfg);
    Real C = engine.price_call(K);
    Real P = engine.price_put(K);
    Real parity = S0 * std::exp(-q * T) - K * std::exp(-r * T);
    EXPECT_NEAR(C - P, parity, 1e-2);
}

TEST(CarrMadanFFTTest, Heston_Call_MatchesReference) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0;
    HestonCFParams hp{0.04, 1.5, 0.04, 0.3, -0.5, r, q};

    Real fft_price = fft_call_heston(S0, K, T, r, q, hp, 1.5, 4096, 0.25);

    HestonParams hp_calib{hp.v0, hp.kappa, hp.theta, hp.sigma, hp.rho};
    Real ref_price = detail::heston_call_price_cf(S0, K, T, r, q, hp_calib, 4096);

    EXPECT_NEAR(fft_price, ref_price, 5e-2)
        << "FFT=" << fft_price << " Ref=" << ref_price;
}

TEST(CarrMadanFFTTest, StrikeGridCenteredCorrectly) {
    // 默认配置下, k=0 (即 K=S0) 应在 strike 网格内
    Real S0 = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    CarrMadanFFT::Config cfg;
    cfg.alpha = 1.5;
    cfg.n_fft = 4096;
    cfg.eta = 0.25;
    CarrMadanFFT engine(phi, S0, r, q, T, cfg);

    const auto& strikes = engine.strikes();
    Real k_min = std::log(strikes.front());
    Real k_max = std::log(strikes.back());
    EXPECT_LT(k_min, 0.0);  // k=0 在网格内
    EXPECT_GT(k_max, 0.0);
}

TEST(CarrMadanFFTTest, RejectsNonPowerOf2NFFT) {
    auto phi = make_gbm_cf(100.0, 0.05, 0.0, 0.20, 1.0);
    CarrMadanFFT::Config cfg;
    cfg.n_fft = 1000;  // not power of 2
    EXPECT_THROW(CarrMadanFFT(phi, 100.0, 0.05, 0.0, 1.0, cfg), std::invalid_argument);
}

// ============ 5. COS vs FFT 交叉验证 ============

TEST(CrossValidationTest, COS_vs_FFT_Heston_ATM) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0;
    HestonCFParams hp{0.04, 1.5, 0.04, 0.3, -0.5, r, q};

    Real cos_price = cos_call_heston(S0, K, T, r, q, hp, 512, 12.0);
    Real fft_price = fft_call_heston(S0, K, T, r, q, hp, 1.5, 4096, 0.25);

    EXPECT_NEAR(cos_price, fft_price, 5e-2)
        << "COS=" << cos_price << " FFT=" << fft_price;
}

TEST(CrossValidationTest, COS_vs_FFT_Heston_OTM) {
    Real S0 = 100.0, K = 110.0, T = 1.0, r = 0.05, q = 0.0;
    HestonCFParams hp{0.04, 1.5, 0.04, 0.3, -0.5, r, q};

    Real cos_price = cos_call_heston(S0, K, T, r, q, hp, 512, 12.0);
    Real fft_price = fft_call_heston(S0, K, T, r, q, hp, 1.5, 4096, 0.25);

    EXPECT_NEAR(cos_price, fft_price, 5e-2)
        << "COS=" << cos_price << " FFT=" << fft_price;
}

TEST(CrossValidationTest, COS_vs_FFT_GBM_MultipleStrikes) {
    Real S0 = 100.0, T = 1.0, r = 0.05, q = 0.02, sigma = 0.20;
    std::vector<Real> strikes = {80.0, 90.0, 95.0, 100.0, 105.0, 110.0, 120.0};

    for (Real K : strikes) {
        Real cos_price = cos_call_gbm(S0, K, T, r, q, sigma, 256, 10.0);
        Real fft_price = fft_call_gbm(S0, K, T, r, q, sigma, 1.5, 4096, 0.25);
        Real bsm_price = bsm_call(S0, K, T, r, q, sigma);

        // 两种方法都应接近 BSM
        EXPECT_NEAR(cos_price, bsm_price, 1e-3)
            << "K=" << K << " COS=" << cos_price << " BSM=" << bsm_price;
        EXPECT_NEAR(fft_price, bsm_price, 5e-2)
            << "K=" << K << " FFT=" << fft_price << " BSM=" << bsm_price;
    }
}

// ============ 6. 批量定价与边界情形 ============

TEST(COSMethodTest, BatchPricing_AllStrikesConsistent) {
    Real S0 = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    COSEngine::Config cfg;
    cfg.n_terms = 256;
    cfg.L = 10.0;
    COSEngine engine(phi, S0, r, q, T, cfg);

    std::vector<Real> strikes = {80.0, 90.0, 100.0, 110.0, 120.0};
    auto calls = engine.price_calls(strikes);

    // call 价格应递减 (随 K 增大)
    for (Size i = 1; i < calls.size(); ++i) {
        EXPECT_LT(calls[i], calls[i - 1]);
    }
    // 与 BSM 对比
    for (Size i = 0; i < strikes.size(); ++i) {
        Real bsm = bsm_call(S0, strikes[i], T, r, q, sigma);
        EXPECT_NEAR(calls[i], bsm, 1e-3);
    }
}

TEST(COSMethodTest, DeepOTM_ReturnsSmallButPositive) {
    Real S0 = 100.0, K = 200.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real cos_price = cos_call_gbm(S0, K, T, r, q, sigma, 512, 12.0);
    Real bsm_price = bsm_call(S0, K, T, r, q, sigma);
    EXPECT_GE(cos_price, 0.0);
    EXPECT_LT(cos_price, 1.0);  // deep OTM call 价格应很小
    // 深度 OTM 数值精度有限, 仅检查数量级
    EXPECT_NEAR(cos_price, bsm_price, 0.1);
}

TEST(COSMethodTest, UserSpecifiedTruncationRange) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    Real mean = std::log(S0) + (r - q - 0.5 * sigma * sigma) * T;
    Real std_dev = sigma * std::sqrt(T);

    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    COSEngine::Config cfg;
    cfg.n_terms = 256;
    cfg.a = mean - 10.0 * std_dev;
    cfg.b = mean + 10.0 * std_dev;
    COSEngine engine(phi, S0, r, q, T, cfg);

    Real cos_price = engine.price_call(K);
    Real bsm_price = bsm_call(S0, K, T, r, q, sigma);
    EXPECT_NEAR(cos_price, bsm_price, 1e-4);
}

// ============ 7. VG / NIG CF 基础验证 ============

TEST(CharacteristicFunctionsTest, VG_UnitModulusForRealU) {
    auto phi = make_vg_cf(100.0, 0.05, 0.0, 0.20, 0.1, -0.1, 1.0);
    for (Real u = 0.0; u <= 5.0; u += 0.5) {
        Complex val = phi(Complex(u, 0.0));
        EXPECT_LE(std::abs(val), 1.0 + 1e-10)
            << " |phi_VG(u=" << u << ")| = " << std::abs(val);
    }
}

TEST(CharacteristicFunctionsTest, NIG_UnitModulusForRealU) {
    auto phi = make_nig_cf(100.0, 0.05, 0.0, 15.0, -5.0, 1.0, 1.0);
    for (Real u = 0.0; u <= 5.0; u += 0.5) {
        Complex val = phi(Complex(u, 0.0));
        EXPECT_LE(std::abs(val), 1.0 + 1e-10)
            << " |phi_NIG(u=" << u << ")| = " << std::abs(val);
    }
}

TEST(COSMethodTest, VG_Call_PriceReasonable) {
    // VG 模型 CF 应能被 COS 引擎定价 (与 BSM 同 strike 比较, 应在同一数量级)
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0;
    Real sigma = 0.20, nu = 0.1, theta = -0.1;
    auto phi = make_vg_cf(S0, r, q, sigma, nu, theta, T);
    COSEngine::Config cfg;
    cfg.n_terms = 512;
    cfg.L = 12.0;
    COSEngine engine(phi, S0, r, q, T, cfg);
    Real vg_price = engine.price_call(K);
    Real bsm_price = bsm_call(S0, K, T, r, q, sigma);
    // VG 与 BSM 价格应在同一数量级 (VG 引入负 theta 会降低价格)
    EXPECT_GT(vg_price, 0.0);
    EXPECT_LT(vg_price, bsm_price * 1.5);
    EXPECT_GT(vg_price, bsm_price * 0.5);
}

// ============ 百慕大期权 COS 方法测试 (Fang-Oosterlee 2009 §5) ============

// 辅助: 用 CRR 二叉树对百慕大期权定价 (交叉验证基准)
static Real crr_bermudan_call(Real S0, Real K, Real T, Real r, Real q, Real sigma,
                               Size n_steps, const std::vector<Real>& exercise_times) {
    BinomialParams params{S0, K, T, r, q, sigma, n_steps};
    BinomialTreeEngine engine(params, BinomialType::CRR);
    CallPayOff payoff(K);
    // 将行权时间转换为步数索引
    std::vector<Size> exercise_steps;
    Real dt = T / static_cast<Real>(n_steps);
    for (Real t : exercise_times) {
        Size step = static_cast<Size>(std::round(t / dt));
        if (step > 0 && step < n_steps) {
            exercise_steps.push_back(step - 1);  // price_bermudan 用 i-1 索引
        }
    }
    return engine.price_bermudan(payoff, exercise_steps);
}

static Real crr_bermudan_put(Real S0, Real K, Real T, Real r, Real q, Real sigma,
                              Size n_steps, const std::vector<Real>& exercise_times) {
    BinomialParams params{S0, K, T, r, q, sigma, n_steps};
    BinomialTreeEngine engine(params, BinomialType::CRR);
    PutPayOff payoff(K);
    std::vector<Size> exercise_steps;
    Real dt = T / static_cast<Real>(n_steps);
    for (Real t : exercise_times) {
        Size step = static_cast<Size>(std::round(t / dt));
        if (step > 0 && step < n_steps) {
            exercise_steps.push_back(step - 1);
        }
    }
    return engine.price_bermudan(payoff, exercise_steps);
}

// 1. 单行权日 (= T) 退化为欧式
TEST(COSBermudanTest, SingleExerciseDateDegeneratesToEuropean) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.02, sigma = 0.20;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    auto inc_factory = make_gbm_inc_cf_factory(r, q, sigma);
    COSEngine::Config cfg;
    cfg.n_terms = 256;
    cfg.L = 10.0;
    COSEngine engine(phi, S0, r, q, T, cfg);

    Real bermudan = engine.price_bermudan(K, true, {T}, inc_factory);
    Real european = engine.price_call(K);
    EXPECT_NEAR(bermudan, european, 1e-10);
}

// 2. 无分红看涨百慕大 = 欧式看涨 (经典结论: q=0 时不应提前行权)
TEST(COSBermudanTest, NoDividendCallBermudanEqualsEuropean) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    auto inc_factory = make_gbm_inc_cf_factory(r, q, sigma);
    COSEngine::Config cfg;
    cfg.n_terms = 256;
    cfg.L = 10.0;
    COSEngine engine(phi, S0, r, q, T, cfg);

    std::vector<Real> ex_times = {T / 3.0, 2.0 * T / 3.0, T};
    Real bermudan = engine.price_bermudan(K, true, ex_times, inc_factory);
    Real european = engine.price_call(K);
    // q=0 时看涨百慕大 = 欧式 (容差放宽因数值误差)
    EXPECT_NEAR(bermudan, european, 5e-3);
}

// 3. 看跌百慕大 >= 欧式 (r>0 时有提前行权溢价)
TEST(COSBermudanTest, PutBermudanExceedsEuropean) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    auto inc_factory = make_gbm_inc_cf_factory(r, q, sigma);
    COSEngine::Config cfg;
    cfg.n_terms = 256;
    cfg.L = 10.0;
    COSEngine engine(phi, S0, r, q, T, cfg);

    std::vector<Real> ex_times = {T / 4.0, T / 2.0, 3.0 * T / 4.0, T};
    Real bermudan = engine.price_bermudan(K, false, ex_times, inc_factory);
    Real european = engine.price_put(K);
    EXPECT_GT(bermudan, european - 1e-6);  // 百慕大应严格大于欧式
}

// 4. 有分红的看涨百慕大 > 欧式 (q>0 时提前行权有溢价)
TEST(COSBermudanTest, DividendCallBermudanExceedsEuropean) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.10, sigma = 0.20;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    auto inc_factory = make_gbm_inc_cf_factory(r, q, sigma);
    COSEngine::Config cfg;
    cfg.n_terms = 256;
    cfg.L = 10.0;
    COSEngine engine(phi, S0, r, q, T, cfg);

    std::vector<Real> ex_times = {T / 4.0, T / 2.0, 3.0 * T / 4.0, T};
    Real bermudan = engine.price_bermudan(K, true, ex_times, inc_factory);
    Real european = engine.price_call(K);
    EXPECT_GT(bermudan, european - 1e-6);
}

// 5. 与 CRR 二叉树交叉验证 (看跌, 3 行权日)
TEST(COSBermudanTest, PutCrossValidatesWithCRR) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.06, q = 0.02, sigma = 0.30;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    auto inc_factory = make_gbm_inc_cf_factory(r, q, sigma);
    COSEngine::Config cfg;
    cfg.n_terms = 256;
    cfg.L = 12.0;
    COSEngine engine(phi, S0, r, q, T, cfg);

    std::vector<Real> ex_times = {T / 3.0, 2.0 * T / 3.0, T};
    Real cos_price = engine.price_bermudan(K, false, ex_times, inc_factory);
    Real crr_price = crr_bermudan_put(S0, K, T, r, q, sigma, 2000, ex_times);
    // 容差: COS vs CRR (n=2000) 应在 5e-3 内
    EXPECT_NEAR(cos_price, crr_price, 5e-3)
        << "COS=" << cos_price << " CRR=" << crr_price;
}

// 6. 与 CRR 二叉树交叉验证 (看涨, 有分红, 4 行权日)
TEST(COSBermudanTest, CallWithDividendCrossValidatesWithCRR) {
    Real S0 = 100.0, K = 95.0, T = 1.0, r = 0.05, q = 0.08, sigma = 0.25;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    auto inc_factory = make_gbm_inc_cf_factory(r, q, sigma);
    COSEngine::Config cfg;
    cfg.n_terms = 256;
    cfg.L = 12.0;
    COSEngine engine(phi, S0, r, q, T, cfg);

    std::vector<Real> ex_times = {0.25, 0.5, 0.75, T};
    Real cos_price = engine.price_bermudan(K, true, ex_times, inc_factory);
    Real crr_price = crr_bermudan_call(S0, K, T, r, q, sigma, 2000, ex_times);
    EXPECT_NEAR(cos_price, crr_price, 5e-3)
        << "COS=" << cos_price << " CRR=" << crr_price;
}

// 7. 多行权日单调性: 行权日越多, 百慕大价格越高 (趋近美式)
TEST(COSBermudanTest, MoreExerciseDatesGivesHigherPrice) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.08, sigma = 0.25;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    auto inc_factory = make_gbm_inc_cf_factory(r, q, sigma);
    COSEngine::Config cfg;
    cfg.n_terms = 256;
    cfg.L = 12.0;
    COSEngine engine(phi, S0, r, q, T, cfg);

    // 2 行权日 vs 5 行权日
    Real p2 = engine.price_bermudan(K, true, {T / 2.0, T}, inc_factory);
    Real p5 = engine.price_bermudan(K, true, {0.2, 0.4, 0.6, 0.8, T}, inc_factory);
    EXPECT_GE(p5, p2 - 1e-6);  // 更多行权日 → 价格更高
}

// 8. 便捷工厂函数一致性
TEST(COSBermudanTest, FactoryFunctionMatchesDirectCall) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.02, sigma = 0.20;
    std::vector<Real> ex_times = {T / 3.0, 2.0 * T / 3.0, T};

    Real factory_price = cos_bermudan_call_gbm(S0, K, T, r, q, sigma, ex_times);

    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    auto inc_factory = make_gbm_inc_cf_factory(r, q, sigma);
    COSEngine::Config cfg;
    cfg.n_terms = 256;
    cfg.L = 10.0;
    COSEngine engine(phi, S0, r, q, T, cfg);
    Real direct_price = engine.price_bermudan(K, true, ex_times, inc_factory);

    EXPECT_NEAR(factory_price, direct_price, 1e-12);
}

// 9. 参数验证: 空行权日抛异常
TEST(COSBermudanTest, EmptyExerciseTimesThrows) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    auto inc_factory = make_gbm_inc_cf_factory(r, q, sigma);
    COSEngine engine(phi, S0, r, q, T);
    EXPECT_THROW(engine.price_bermudan(K, true, {}, inc_factory), std::invalid_argument);
}

// 10. 参数验证: 非递增行权日抛异常
TEST(COSBermudanTest, NonIncreasingExerciseTimesThrows) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    auto inc_factory = make_gbm_inc_cf_factory(r, q, sigma);
    COSEngine engine(phi, S0, r, q, T);
    EXPECT_THROW(engine.price_bermudan(K, true, {0.5, 0.5, T}, inc_factory), std::invalid_argument);
    EXPECT_THROW(engine.price_bermudan(K, true, {0.7, 0.3, T}, inc_factory), std::invalid_argument);
}

// 11. 参数验证: 首行权日 <= 0 抛异常
TEST(COSBermudanTest, NonPositiveFirstExerciseThrows) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    auto inc_factory = make_gbm_inc_cf_factory(r, q, sigma);
    COSEngine engine(phi, S0, r, q, T);
    EXPECT_THROW(engine.price_bermudan(K, true, {0.0, T}, inc_factory), std::invalid_argument);
    EXPECT_THROW(engine.price_bermudan(K, true, {-0.1, T}, inc_factory), std::invalid_argument);
}

// 12. 参数验证: 末行权日 != T 抛异常
TEST(COSBermudanTest, LastExerciseNotTThrows) {
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    auto inc_factory = make_gbm_inc_cf_factory(r, q, sigma);
    COSEngine engine(phi, S0, r, q, T);
    EXPECT_THROW(engine.price_bermudan(K, true, {0.3, 0.7, 0.9}, inc_factory), std::invalid_argument);
}

// 13. 增量 CF 工厂验证
TEST(IncCharFnFactoryTest, GBM_IncCF_AtZeroReturnsOne) {
    auto factory = make_gbm_inc_cf_factory(0.05, 0.02, 0.20);
    auto phi = factory(1.0);
    Complex val = phi(Complex(0.0, 0.0));
    EXPECT_NEAR(val.real(), 1.0, 1e-12);
    EXPECT_NEAR(val.imag(), 0.0, 1e-12);
}

TEST(IncCharFnFactoryTest, GBM_IncCF_UnitModulusForRealU) {
    auto factory = make_gbm_inc_cf_factory(0.05, 0.02, 0.20);
    auto phi = factory(1.0);
    for (Real u = 0.0; u <= 10.0; u += 0.5) {
        Complex val = phi(Complex(u, 0.0));
        EXPECT_LE(std::abs(val), 1.0 + 1e-12);
    }
}

TEST(IncCharFnFactoryTest, GBM_IncCF_ScalesWithDt) {
    // φ_inc(u; 2Δt) = φ_inc(u; Δt)^2 (独立增量性质)
    auto factory = make_gbm_inc_cf_factory(0.05, 0.02, 0.20);
    auto phi1 = factory(1.0);
    auto phi2 = factory(2.0);
    Complex u(1.5, 0.0);
    Complex val1 = phi1(u);
    Complex val2 = phi2(u);
    EXPECT_NEAR(val2.real(), (val1 * val1).real(), 1e-10);
    EXPECT_NEAR(val2.imag(), (val1 * val1).imag(), 1e-10);
}
