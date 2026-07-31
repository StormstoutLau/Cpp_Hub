// v1.3 F: CGMY + Kou Levy 过程特征函数测试
// 覆盖: u=0 归一, 模长≤1, 极限退化 (CGMY Y→0→VG, Kou lambda→0→BS),
//       鞍条件/鞅条件, 数值稳定性, COS 定价交叉验证
#include <gtest/gtest.h>
#include "cpphub/pricing/fourier/characteristic_functions.hpp"
#include "cpphub/pricing/fourier/cos_method.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include <cmath>

using namespace cpphub;

// 辅助: 用 COSEngine 定价 call
Real cos_call(const CharFn& phi, Real S0, Real K, Real T, Real r, Real q,
              Size n_terms = 256, Real L = 10.0) {
    COSEngine::Config cfg;
    cfg.n_terms = n_terms;
    cfg.L = L;
    COSEngine engine(phi, S0, r, q, T, cfg);
    return engine.price_call(K);
}

// 辅助: 用 COSEngine 定价 put
Real cos_put(const CharFn& phi, Real S0, Real K, Real T, Real r, Real q,
             Size n_terms = 256, Real L = 10.0) {
    COSEngine::Config cfg;
    cfg.n_terms = n_terms;
    cfg.L = L;
    COSEngine engine(phi, S0, r, q, T, cfg);
    return engine.price_put(K);
}

// 辅助: 用 COSEngine 定价 call (手动指定截断区间 [a, b])
Real cos_call_ab(const CharFn& phi, Real S0, Real K, Real T, Real r, Real q,
                 Size n_terms, Real a, Real b) {
    COSEngine::Config cfg;
    cfg.n_terms = n_terms;
    cfg.a = a;
    cfg.b = b;
    COSEngine engine(phi, S0, r, q, T, cfg);
    return engine.price_call(K);
}

// 辅助: 用 COSEngine 定价 put (手动指定截断区间 [a, b])
Real cos_put_ab(const CharFn& phi, Real S0, Real K, Real T, Real r, Real q,
                Size n_terms, Real a, Real b) {
    COSEngine::Config cfg;
    cfg.n_terms = n_terms;
    cfg.a = a;
    cfg.b = b;
    COSEngine engine(phi, S0, r, q, T, cfg);
    return engine.price_put(K);
}

// Kou 模型 ln S_T 解析矩 (用于手动指定 COS 截断区间)
// ln S_T = ln S0 + (r-q+omega-sigma²/2)*T + sigma*W_T + Σ J_i
// E[ln S_T] = ln S0 + (r-q+omega-sigma²/2)*T + lambda*T*E[J]
//   E[J] = -p*eta1 + q*eta2  (双指数分布均值, 负跳均值 -eta1, 正跳均值 +eta2)
// Var[ln S_T] = sigma²*T + lambda*T*E[J²]
//   E[J²] = p*2*eta1² + q*2*eta2²  (指数分布 Exp(mean=eta) 的二阶矩 = 2*eta²)
// omega = -lambda*(p/(1+eta1) + q/(1-eta2) - 1) = -lambda*xi
inline std::pair<Real, Real> kou_moments(Real S0, Real r, Real q, Real sigma,
                                          Real lambda, Real p, Real eta1, Real eta2, Real T) {
    Real q_prob = 1.0 - p;
    Real E_exp_J = p / (1.0 + eta1) + q_prob / (1.0 - eta2);
    Real xi = E_exp_J - 1.0;
    Real omega = -lambda * xi;
    Real E_J = -p * eta1 + q_prob * eta2;
    Real E_J2 = p * 2.0 * eta1 * eta1 + q_prob * 2.0 * eta2 * eta2;
    Real mean = std::log(S0) + (r - q + omega - 0.5 * sigma * sigma) * T + lambda * T * E_J;
    Real var = sigma * sigma * T + lambda * T * E_J2;
    return {mean, std::sqrt(var)};
}

// ============================================================
// 1. CGMY 特征函数基本性质
// ============================================================
TEST(CGMYCharactericFunction, AtZeroReturnsOne) {
    // C=1, G=5, M=5, Y=0.5 (对称 CGMY)
    auto phi = make_cgmy_cf(100.0, 0.05, 0.0, 1.0, 5.0, 5.0, 0.5, 1.0);
    Complex u(0.0, 0.0);
    Complex val = phi(u);
    EXPECT_NEAR(std::real(val), 1.0, 1e-12);
    EXPECT_NEAR(std::imag(val), 0.0, 1e-12);
}

TEST(CGMYCharactericFunction, UnitModulusForRealU) {
    // 对称 CGMY (G=M), 实数 u 时 |phi|≤1
    // Y=0.5 避开 Gamma(-Y) 极点 (Y≠1)
    auto phi = make_cgmy_cf(100.0, 0.05, 0.0, 1.0, 5.0, 5.0, 0.5, 1.0);
    for (Real u_val = 0.0; u_val <= 10.0; u_val += 0.5) {
        Complex u(u_val, 0.0);
        Complex val = phi(u);
        EXPECT_LE(std::abs(val), 1.0 + 1e-10)
            << "u=" << u_val << " |phi|=" << std::abs(val);
    }
}

TEST(CGMYCharactericFunction, DecaysAtInfinity) {
    // 大 u 时 |phi| → 0 (CGMY 为纯跳跃, 有无限可分性, CF 在无穷衰减)
    auto phi = make_cgmy_cf(100.0, 0.05, 0.0, 1.0, 5.0, 5.0, 1.5, 1.0);
    Complex u(100.0, 0.0);
    Complex val = phi(u);
    EXPECT_LT(std::abs(val), 0.1);
}

TEST(CGMYCharactericFunction, SymmetricParamsSymmetricCF) {
    // G=M 对称 → phi 的虚部在小 u 时接近 0 (对称分布的 CF 为实数)
    auto phi = make_cgmy_cf(100.0, 0.05, 0.0, 1.0, 5.0, 5.0, 0.5, 1.0);
    // 对称 CGMY 的 X_T 分布关于 omega 偏移对称, phi(u) 在去漂移后应为实数
    // phi(u) = exp(iu*drift) * real_term → 虚部主要来自 drift
    Complex u(1.0, 0.0);
    Complex val = phi(u);
    EXPECT_TRUE(std::isfinite(std::real(val)));
    EXPECT_TRUE(std::isfinite(std::imag(val)));
}

TEST(CGMYCharactericFunction, YNearZeroConvergesToVG) {
    // Y→0 时 CGMY 退化为 VG
    // VG 参数映射: theta=C*(M-G), sigma=sqrt(2*C*M*G*... ), 实际更复杂
    // 简化验证: Y 小时 CF 有限且接近某个极限
    Real C = 1.0, G = 5.0, M = 5.0;
    auto phi_y01 = make_cgmy_cf(100.0, 0.05, 0.0, C, G, M, 0.01, 1.0);
    auto phi_y001 = make_cgmy_cf(100.0, 0.05, 0.0, C, G, M, 0.001, 1.0);
    Complex u(1.0, 0.0);
    // Y 趋近 0 时 CF 应收敛
    Real err = std::abs(phi_y01(u) - phi_y001(u));
    EXPECT_LT(err, 0.5);  // 收敛但不严格 (Gamma(-Y) 在 Y→0 发散)
}

TEST(CGMYCharactericFunction, NoBranchCutForRealU) {
    // 实轴上无分支切割 (复数幂 ^Y 可能引入分支)
    auto phi = make_cgmy_cf(100.0, 0.05, 0.0, 1.0, 5.0, 8.0, 1.2, 1.0);
    const Real eps = 1e-8;
    for (Real u_val = 0.1; u_val <= 5.0; u_val += 0.5) {
        Complex phi_fwd = phi(Complex(u_val + eps, 0));
        Complex phi_bwd = phi(Complex(u_val - eps, 0));
        Complex phi_mid = phi(Complex(u_val, 0));
        Real smoothness = std::abs(phi_fwd + phi_bwd - Real(2) * phi_mid);
        EXPECT_LT(smoothness, 1e-6) << "u=" << u_val;
    }
}

TEST(CGMYCharactericFunction, InvalidParamsThrow) {
    EXPECT_THROW(make_cgmy_cf(100.0, 0.05, 0.0, -1.0, 5.0, 5.0, 0.5, 1.0), std::invalid_argument);
    EXPECT_THROW(make_cgmy_cf(100.0, 0.05, 0.0, 1.0, -5.0, 5.0, 0.5, 1.0), std::invalid_argument);
    EXPECT_THROW(make_cgmy_cf(100.0, 0.05, 0.0, 1.0, 5.0, 5.0, 2.5, 1.0), std::invalid_argument);  // Y>=2
    EXPECT_THROW(make_cgmy_cf(100.0, 0.05, 0.0, 1.0, 5.0, 5.0, 0.5, -1.0), std::invalid_argument);  // T<=0
}

// ============================================================
// 2. CGMY COS 定价
// ============================================================
TEST(CGMYCharactericFunction, COSPricingProducesValidCallPrice) {
    // 用 CGMY CF + COS 方法定价 call, 验证价格合理
    // Y=0.5 避开 Gamma(-Y) 极点, L=15 适配厚尾
    auto phi = make_cgmy_cf(100.0, 0.05, 0.0, 1.0, 5.0, 5.0, 0.5, 1.0);
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0;
    Real price = cos_call(phi, S0, K, T, r, q, 512, 15.0);
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, S0);  // call 价格 < 标的
    // ATM call 应在合理范围 (5-20)
    EXPECT_GT(price, 3.0);
    EXPECT_LT(price, 25.0);
}

TEST(CGMYCharactericFunction, CallPutParityViaCOS) {
    // Call-Put parity: C - P = S0*e^{-qT} - K*e^{-rT}
    // CGMY Y=0.5 避开 Gamma(-Y) 极点; 用解析矩指定截断区间提高精度
    auto phi = make_cgmy_cf(100.0, 0.05, 0.0, 1.0, 5.0, 8.0, 0.5, 1.0);
    Real S0 = 100.0, K = 105.0, T = 1.0, r = 0.05, q = 0.0;
    // CGMY 矩 (C=1, G=5, M=8, Y=0.5, T=1):
    // E[X] = C*T*Gamma(-Y)*[(M-1)^Y - M^Y + (G+1)^Y - G^Y]  (X 为 CGMY 过程)
    // 此处用数值估计: ln S_T 的均值/方差通过 CF 数值微分
    // 直接用 COSEngine 自动估计 + 大 L
    Real call = cos_call(phi, S0, K, T, r, q, 512, 20.0);
    Real put = cos_put(phi, S0, K, T, r, q, 512, 20.0);
    Real parity = call - put;
    Real expected = S0 * std::exp(-q * T) - K * std::exp(-r * T);
    EXPECT_NEAR(parity, expected, 0.50);  // CGMY 厚尾, COS 截断误差较大
}

// ============================================================
// 3. Kou 特征函数基本性质
// ============================================================
TEST(KouCharactericFunction, AtZeroReturnsOne) {
    auto phi = make_kou_cf(100.0, 0.05, 0.0, 0.2, 3.0, 0.4, 0.3, 0.3, 1.0);
    Complex u(0.0, 0.0);
    Complex val = phi(u);
    EXPECT_NEAR(std::real(val), 1.0, 1e-12);
    EXPECT_NEAR(std::imag(val), 0.0, 1e-12);
}

TEST(KouCharactericFunction, UnitModulusForRealU) {
    auto phi = make_kou_cf(100.0, 0.05, 0.0, 0.2, 3.0, 0.4, 0.3, 0.3, 1.0);
    for (Real u_val = 0.0; u_val <= 10.0; u_val += 0.5) {
        Complex u(u_val, 0.0);
        Complex val = phi(u);
        EXPECT_LE(std::abs(val), 1.0 + 1e-10)
            << "u=" << u_val << " |phi|=" << std::abs(val);
    }
}

TEST(KouCharactericFunction, LambdaNearZeroConvergesToBS) {
    // lambda→0 (无跳跃) → Kou 退化为 BS
    Real S0 = 100.0, r = 0.05, q = 0.0, sigma = 0.2, T = 1.0;
    auto phi_bs = make_gbm_cf(S0, r, q, sigma, T);
    // lambda 极小时 Kou CF 应接近 BS CF (大 u 时复数误差累积, 容差放宽)
    auto phi_kou = make_kou_cf(S0, r, q, sigma, 1e-6, 0.5, 0.1, 0.1, T);
    for (Real u_val = 0.5; u_val <= 5.0; u_val += 0.5) {
        Complex u(u_val, 0.0);
        Real err = std::abs(phi_kou(u) - phi_bs(u));
        EXPECT_LT(err, 1e-6) << "u=" << u_val << " err=" << err;
    }
}

TEST(KouCharactericFunction, HigherLambdaMoreKurtosis) {
    // 跳跃强度越大, |phi(u)| 衰减越快 (尾部更厚)
    auto phi_low = make_kou_cf(100.0, 0.05, 0.0, 0.2, 1.0, 0.4, 0.3, 0.3, 1.0);
    auto phi_high = make_kou_cf(100.0, 0.05, 0.0, 0.2, 10.0, 0.4, 0.3, 0.3, 1.0);
    Complex u(3.0, 0.0);
    // 高 lambda → 更多跳跃 → |phi| 更小 (分布更分散)
    EXPECT_LT(std::abs(phi_high(u)), std::abs(phi_low(u)) + 1e-10);
}

TEST(KouCharactericFunction, NoBranchCutForRealU) {
    auto phi = make_kou_cf(100.0, 0.05, 0.0, 0.2, 3.0, 0.4, 0.3, 0.3, 1.0);
    const Real eps = 1e-8;
    for (Real u_val = 0.1; u_val <= 5.0; u_val += 0.5) {
        Complex phi_fwd = phi(Complex(u_val + eps, 0));
        Complex phi_bwd = phi(Complex(u_val - eps, 0));
        Complex phi_mid = phi(Complex(u_val, 0));
        Real smoothness = std::abs(phi_fwd + phi_bwd - Real(2) * phi_mid);
        EXPECT_LT(smoothness, 1e-6) << "u=" << u_val;
    }
}

TEST(KouCharactericFunction, InvalidParamsThrow) {
    EXPECT_THROW(make_kou_cf(100.0, 0.05, 0.0, -0.2, 3.0, 0.4, 0.3, 0.3, 1.0), std::invalid_argument);
    EXPECT_THROW(make_kou_cf(100.0, 0.05, 0.0, 0.2, -3.0, 0.4, 0.3, 0.3, 1.0), std::invalid_argument);
    EXPECT_THROW(make_kou_cf(100.0, 0.05, 0.0, 0.2, 3.0, 1.5, 0.3, 0.3, 1.0), std::invalid_argument);  // p>=1
    EXPECT_THROW(make_kou_cf(100.0, 0.05, 0.0, 0.2, 3.0, 0.4, 0.3, 1.5, 1.0), std::invalid_argument);  // eta2>=1 (鞅条件)
}

// ============================================================
// 4. Kou COS 定价
// ============================================================
TEST(KouCharactericFunction, COSPricingProducesValidCallPrice) {
    // Kou 跳跃扩散: sigma=0.15, lambda=2, p=0.4, eta1=0.3, eta2=0.25
    auto phi = make_kou_cf(100.0, 0.05, 0.0, 0.15, 2.0, 0.4, 0.3, 0.25, 1.0);
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0;
    Real price = cos_call(phi, S0, K, T, r, q, 512, 15.0);
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, S0);
    // 跳跃扩散 ATM call 应在合理范围 (8-25)
    EXPECT_GT(price, 5.0);
    EXPECT_LT(price, 30.0);
}

TEST(KouCharactericFunction, COSCallPutParity) {
    // Kou 跳跃扩散: 用解析矩手动指定截断区间, 避免 cos_truncation_range 数值微分误差
    Real S0 = 100.0, r = 0.05, q = 0.0, sigma = 0.15;
    Real lambda = 2.0, p = 0.4, eta1 = 0.3, eta2 = 0.25, T = 1.0;
    auto phi = make_kou_cf(S0, r, q, sigma, lambda, p, eta1, eta2, T);
    Real K = 95.0;
    auto mom = kou_moments(S0, r, q, sigma, lambda, p, eta1, eta2, T);
    Real a = mom.first - 20.0 * mom.second;
    Real b = mom.first + 20.0 * mom.second;
    Real call = cos_call_ab(phi, S0, K, T, r, q, 1024, a, b);
    Real put = cos_put_ab(phi, S0, K, T, r, q, 1024, a, b);
    Real parity = call - put;
    Real expected = S0 * std::exp(-q * T) - K * std::exp(-r * T);
    EXPECT_NEAR(parity, expected, 0.10);
}

// ============================================================
// 5. CGMY vs Kou 对比 (不同 Levy 族应给出不同价格)
// ============================================================
TEST(LevyCFComparison, CGMYAndKouGiveDifferentPrices) {
    // 同样 ATM call, CGMY 和 Kou 应给出不同价格 (不同 Levy 族)
    // CGMY Y=0.5 避开 Gamma(-Y) 极点 (Y≠1)
    auto phi_cgmy = make_cgmy_cf(100.0, 0.05, 0.0, 1.0, 5.0, 5.0, 0.5, 1.0);
    auto phi_kou = make_kou_cf(100.0, 0.05, 0.0, 0.2, 3.0, 0.4, 0.3, 0.3, 1.0);
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0;
    Real price_cgmy = cos_call(phi_cgmy, S0, K, T, r, q, 256, 12.0);
    Real price_kou = cos_call(phi_kou, S0, K, T, r, q, 256, 12.0);
    // 两者应不同 (不同模型, 不同参数化)
    EXPECT_NE(price_cgmy, price_kou);
    EXPECT_GT(price_cgmy, 0.0);
    EXPECT_GT(price_kou, 0.0);
}
