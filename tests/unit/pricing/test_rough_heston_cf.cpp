// Rough Heston 特征函数单元测试 (rough_heston_cf.hpp)
// 覆盖: 归一性 / 模长有界 / 共轭对称 / H→0.5 与 H=0.5 退化 (vs heston_characteristic_function)
//       / COS Call-Put 平价 / COS vs Heston 一致性 / ATM skew 粗糙度 /
//       rho 单调性 / 参数校验 / 大 |u| 稳定性 / 网格收敛 / 鞅性质
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <complex>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/pricing/analytic/rough_heston_cf.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include "cpphub/pricing/fourier/characteristic_functions.hpp"
#include "cpphub/pricing/fourier/cos_method.hpp"
#include "cpphub/calibration/calibrator.hpp"  // bsm_implied_vol

using namespace cpphub;

namespace {

RoughHestonCFParams base_params() {
    return {0.1, 1.5, 0.04, 0.3, -0.5, 0.04, 100.0, 0.03, 0.0, 1.0, 500};
}

HestonCFParams heston_params(const RoughHestonCFParams& p) {
    return {p.v0, p.kappa, p.theta, p.sigma, p.rho, p.r, p.q};
}

Complex heston_cf(Complex u, const RoughHestonCFParams& p) {
    return heston_characteristic_function(u, p.T, p.S0, heston_params(p));
}

// COS 价格 → BSM 隐含波动率
// n_terms=256 为 COS 引擎默认 (Fang-Oosterlee 2009, 256 项 ~0.1ms/solve);
// 大 u 处分数阶 Riccati 数值稳定性受 |2c·h|·(Δt^α/Γ(α+2)) 限制, 256 项时
// u_max ≈ 155, 位于稳定区 (n_terms=512 时 u 可达 ~306, 超出稳定阈值 ~205).
Real call_iv(const RoughHestonCFParams& p, Real K, Size n_terms = 256, Real L = 12.0) {
    auto phi = make_rough_heston_cf(p);
    COSEngine::Config cfg;
    cfg.n_terms = n_terms;
    cfg.L = L;
    COSEngine engine(phi, p.S0, p.r, p.q, p.T, cfg);
    Real price = engine.price_call(K);
    return bsm_implied_vol(price, p.S0, K, p.T, p.r, p.q, true);
}

Real skew_measure(const RoughHestonCFParams& p, Real K_low, Real K_high) {
    return call_iv(p, K_low) - call_iv(p, K_high);
}

}  // namespace

// 1. u=0 归一: φ(0) = 1
TEST(RoughHestonCFTest, PhiAtZeroIsOne) {
    auto p = base_params();
    Complex phi = rough_heston_characteristic_function(Complex(0, 0), p);
    EXPECT_NEAR(std::real(phi), 1.0, 1e-10);
    EXPECT_NEAR(std::imag(phi), 0.0, 1e-10);
}

// 2. 模长 ≤ 1 对多个 u
TEST(RoughHestonCFTest, ModulusBoundedByOne) {
    auto p = base_params();
    for (Real u : {0.5, 1.0, 2.0, 5.0, 10.0}) {
        Complex phi = rough_heston_characteristic_function(Complex(u, 0), p);
        EXPECT_LE(std::abs(phi), 1.0 + 1e-4) << "u=" << u;
        EXPECT_TRUE(std::isfinite(std::real(phi))) << "u=" << u;
    }
}

// 3. 共轭对称: φ(-u) = conj(φ(u))
TEST(RoughHestonCFTest, ConjugateSymmetry) {
    auto p = base_params();
    for (Real u : {0.5, 1.5, 3.0}) {
        Complex pp = rough_heston_characteristic_function(Complex(u, 0), p);
        Complex pn = rough_heston_characteristic_function(Complex(-u, 0), p);
        EXPECT_NEAR(std::real(pn), std::real(pp), 1e-5) << "u=" << u;
        EXPECT_NEAR(std::imag(pn), -std::imag(pp), 1e-5) << "u=" << u;
    }
}

// 4. H→0.5 退化 (H=0.499, α≈1): 与标准 Heston 一致 (容差 1e-3, 因 α≠1 有小偏差)
TEST(RoughHestonCFTest, NearHalfDegeneration) {
    auto p = base_params();
    p.H = 0.499;
    p.n_steps = 1000;
    for (Real u : {0.5, 1.0, 2.0}) {
        Complex rough = rough_heston_characteristic_function(Complex(u, 0), p);
        Complex heston = heston_cf(Complex(u, 0), p);
        EXPECT_NEAR(std::real(rough), std::real(heston), 1e-3) << "u=" << u;
        EXPECT_NEAR(std::imag(rough), std::imag(heston), 1e-3) << "u=" << u;
    }
}

// 5. H=0.5 精确退化 (α=1): 退化为标准 Heston, 容差 1e-4
TEST(RoughHestonCFTest, ExactHalfDegeneration) {
    auto p = base_params();
    p.H = 0.5;
    p.n_steps = 1000;
    for (Real u : {0.5, 1.0, 2.0, 3.0}) {
        Complex rough = rough_heston_characteristic_function(Complex(u, 0), p);
        Complex heston = heston_cf(Complex(u, 0), p);
        EXPECT_NEAR(std::real(rough), std::real(heston), 1e-4) << "u=" << u;
        EXPECT_NEAR(std::imag(rough), std::imag(heston), 1e-4) << "u=" << u;
    }
}

// 6. COS Call-Put 平价: C - P = S₀e^{-qT} - Ke^{-rT} (容差 1e-4)
TEST(RoughHestonCFTest, CallPutParity) {
    auto p = base_params();  // H=0.1
    auto phi = make_rough_heston_cf(p);
    COSEngine::Config cfg;
    cfg.n_terms = 256;
    cfg.L = 12.0;
    COSEngine engine(phi, p.S0, p.r, p.q, p.T, cfg);
    for (Real K : {90.0, 100.0, 110.0}) {
        Real C = engine.price_call(K);
        Real P = engine.price_put(K);
        Real parity = p.S0 * std::exp(-p.q * p.T) - K * std::exp(-p.r * p.T);
        EXPECT_NEAR(C - P, parity, 1e-4) << "K=" << K;
    }
}

// 7. COS vs Heston 一致性: H=0.5 时 Rough Heston COS 价格 == Heston COS 价格 (容差 1e-4)
TEST(RoughHestonCFTest, COSMatchHeston) {
    auto p = base_params();
    p.H = 0.5;
    p.n_steps = 1000;
    auto heston_phi = make_heston_cf(p.S0, p.r, p.q, heston_params(p), p.T);
    auto rng = cos_truncation_range(heston_phi, p.S0, p.T, 12.0);
    COSEngine::Config cfg;
    cfg.n_terms = 256;
    cfg.L = 12.0;
    cfg.a = rng.first;  // 两引擎使用相同截断区间, 仅 CF 不同
    cfg.b = rng.second;
    COSEngine hest_engine(heston_phi, p.S0, p.r, p.q, p.T, cfg);
    COSEngine rough_engine(make_rough_heston_cf(p), p.S0, p.r, p.q, p.T, cfg);
    for (Real K : {95.0, 100.0, 105.0}) {
        EXPECT_NEAR(rough_engine.price_call(K), hest_engine.price_call(K), 1e-4)
            << "K=" << K;
    }
}

// 8. ATM smile 倾斜: H<0.5 时 IV skew 比 Heston (H=0.5) 更陡 (rough volatility 特征)
TEST(RoughHestonCFTest, RoughSkewSteeperThanHeston) {
    auto rough = base_params();
    rough.T = 0.5;
    auto heston = base_params();
    heston.H = 0.5;
    heston.n_steps = 1000;
    heston.T = 0.5;

    Real skew_rough = skew_measure(rough, 90.0, 110.0);
    Real skew_heston = skew_measure(heston, 90.0, 110.0);
    // rho<0 → 下行微笑 (低行权价 IV 更高), skew = IV(90)-IV(110) > 0
    EXPECT_GT(skew_rough, 0.0) << "rough skew=" << skew_rough;
    EXPECT_GT(skew_heston, 0.0) << "heston skew=" << skew_heston;
    EXPECT_GT(skew_rough, skew_heston)
        << "rough skew=" << skew_rough << " heston skew=" << skew_heston;
}

// 9. rho 单调性: rho→-1 时 IV skew 加深
TEST(RoughHestonCFTest, RhoDeepensSkew) {
    auto mild = base_params();
    mild.rho = -0.3;
    auto strong = base_params();
    strong.rho = -0.9;

    Real skew_mild = skew_measure(mild, 90.0, 110.0);
    Real skew_strong = skew_measure(strong, 90.0, 110.0);
    EXPECT_GT(skew_strong, 0.0) << "strong skew=" << skew_strong;
    EXPECT_GT(skew_strong, skew_mild)
        << "mild skew=" << skew_mild << " strong skew=" << skew_strong;
}

// 10. 参数校验: H>0.5 / H≤0 / sigma≤0 / rho 越界 / S0,T,v0 非法 抛异常
TEST(RoughHestonCFTest, ParameterValidation) {
    auto good = base_params();
    EXPECT_NO_THROW(validate_rough_heston_cf_params(good));
    // H=0.5 合法 (α=1 退化验证需要)
    auto half = good;
    half.H = 0.5;
    EXPECT_NO_THROW(validate_rough_heston_cf_params(half));

    auto bad = [](Real H, Real sigma, Real rho, Real S0, Real T, Real v0) {
        return RoughHestonCFParams{H, 1.5, 0.04, sigma, rho, v0, S0, 0.03, 0.0, T, 200};
    };
    EXPECT_THROW(validate_rough_heston_cf_params(bad(0.6, 0.3, -0.5, 100, 1, 0.04)),
                 std::invalid_argument);   // H > 0.5
    EXPECT_THROW(validate_rough_heston_cf_params(bad(0.0, 0.3, -0.5, 100, 1, 0.04)),
                 std::invalid_argument);   // H = 0
    EXPECT_THROW(validate_rough_heston_cf_params(bad(-0.1, 0.3, -0.5, 100, 1, 0.04)),
                 std::invalid_argument);   // H < 0
    EXPECT_THROW(validate_rough_heston_cf_params(bad(0.1, 0.0, -0.5, 100, 1, 0.04)),
                 std::invalid_argument);   // sigma = 0
    EXPECT_THROW(validate_rough_heston_cf_params(bad(0.1, -1.0, -0.5, 100, 1, 0.04)),
                 std::invalid_argument);   // sigma < 0
    EXPECT_THROW(validate_rough_heston_cf_params(bad(0.1, 0.3, 1.1, 100, 1, 0.04)),
                 std::invalid_argument);   // rho > 1
    EXPECT_THROW(validate_rough_heston_cf_params(bad(0.1, 0.3, -1.5, 100, 1, 0.04)),
                 std::invalid_argument);   // rho < -1
    EXPECT_THROW(validate_rough_heston_cf_params(bad(0.1, 0.3, -0.5, 0.0, 1, 0.04)),
                 std::invalid_argument);   // S0 = 0
    EXPECT_THROW(validate_rough_heston_cf_params(bad(0.1, 0.3, -0.5, 100, 0.0, 0.04)),
                 std::invalid_argument);   // T = 0
    EXPECT_THROW(validate_rough_heston_cf_params(bad(0.1, 0.3, -0.5, 100, 1, -0.1)),
                 std::invalid_argument);   // v0 < 0

    // 非法参数在 CF / 工厂入口同样抛异常
    auto invalid = bad(0.6, 0.3, -0.5, 100, 1, 0.04);
    EXPECT_THROW(make_rough_heston_cf(invalid), std::invalid_argument);
    EXPECT_THROW(rough_heston_characteristic_function(Complex(1, 0), invalid),
                 std::invalid_argument);
}

// 11. 数值稳定性: 大 |u| (|u|=100) 不产生 NaN/Inf
// 注: 仅测试实值 u. u=100i (即 E[S_T^{-100}]) 对应 Heston 矩爆炸
// (E[S_T^p] 对 p<-1 发散), 特征函数本应发散, 不属于稳定性范畴.
TEST(RoughHestonCFTest, LargeUFinite) {
    auto p = base_params();
    for (Complex u : {Complex(100, 0), Complex(-100, 0)}) {
        Complex phi = rough_heston_characteristic_function(u, p);
        EXPECT_TRUE(std::isfinite(std::real(phi))) << "u=" << u;
        EXPECT_TRUE(std::isfinite(std::imag(phi))) << "u=" << u;
        EXPECT_GE(std::abs(phi), 0.0);
        EXPECT_LT(std::abs(phi), 1.0);
    }
}

// 12. 网格收敛: n_steps 100→500, CF 值收敛 (差值 < 1e-4)
TEST(RoughHestonCFTest, GridConvergence) {
    auto p = base_params();
    Complex u(1.0, 0);
    p.n_steps = 100;
    Complex phi100 = rough_heston_characteristic_function(u, p);
    p.n_steps = 500;
    Complex phi500 = rough_heston_characteristic_function(u, p);
    EXPECT_NEAR(std::abs(phi100 - phi500), 0.0, 1e-4)
        << "phi100=" << phi100 << " phi500=" << phi500;
}

// 13. 鞅性质: φ(-i) = E[S_T] = S₀e^{(r-q)T}
TEST(RoughHestonCFTest, MartingaleSpot) {
    auto p = base_params();
    Complex phi = rough_heston_characteristic_function(Complex(0, -1), p);
    Real expected = p.S0 * std::exp((p.r - p.q) * p.T);
    EXPECT_NEAR(std::real(phi), expected, 1e-6);
    EXPECT_NEAR(std::imag(phi), 0.0, 1e-6);
}
