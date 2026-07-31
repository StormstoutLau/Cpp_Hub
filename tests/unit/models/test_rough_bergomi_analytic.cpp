// RoughBergomi 解析层单元测试
// 覆盖: 累积量 + 近似特征函数 + 参数校验 + COS 半解析定价 + MC 一致性 + log(v) 矩
//
// 测试分组:
//   1. 特征函数性质 (7 cases)  — 单位模, 衰减, 共轭对称, 有界, 平滑, GBM 退化, BSM 一致
//   2. 累积量公式 (5 cases)    — 符号, H→0.5 退化, xi0 缩放, 矩转换
//   3. 参数校验 (1 case)       — 非法参数抛异常
//   4. COS 定价 (4 cases)      — COS vs MC, Call-Put 平价, T→0 退化, 价格合理性
//   5. log(v) 矩 (2 cases)     — 闭式矩 + 采样矩 (RLFbmSampler)
// 总计: 19 cases
//
// 已知参考实现特性 (rough_bergomi.hpp, 不可修改):
//   * rbergomi_price_european 的 log-Euler 不含 (r-q) 漂移, E[S_T]≈S_0 而非 S_0 e^{(r-q)T}.
//     因此 COS-MC 一致性测试在 r=q=0 下进行 (两者漂移均为 0, 模型一致).
//   * RLFbmSampler 按 sqrt(2H+1) 归一化构造协方差, 样本 Var(log v_t) ≈ η²·C_{nn}
//     (C 为采样器协方差), 而非 η²·t^{2H}. 方差校验对采样器自洽协方差进行.
//   * 累积量 CF 的 c4 项 (exp 内 +u⁴·c4/24) 会使 |φ| 突破 1 导致 COS 不稳定,
//     定价 CF 仅保留 c1,c2,c3 (见 rough_bergomi_cf.hpp 头注释).

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include "cpphub/core/types.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/pricing/analytic/rough_bergomi_cf.hpp"
#include "cpphub/pricing/fourier/characteristic_functions.hpp"
#include "cpphub/pricing/fourier/cos_method.hpp"
#include "cpphub/models/diffusion/rough_bergomi.hpp"

using namespace cpphub;

namespace {

// 基准参数 (与任务说明一致)
RoughBergomiCFParams base_params() {
    return {0.1, 0.3, -0.7, 0.04, 100.0, 0.03, 0.0, 1.0};
}

}  // namespace

// ============ 1. 特征函数性质 ============

TEST(RoughBergomiCFTest, CFUnitModulus) {
    // φ(0) = 1 (单位模)
    auto params = base_params();
    Complex phi = rough_bergomi_characteristic_function(Complex(0.0, 0.0), params);
    EXPECT_NEAR(std::real(phi), 1.0, 1e-12);
    EXPECT_NEAR(std::imag(phi), 0.0, 1e-12);
}

TEST(RoughBergomiCFTest, CFDecay) {
    // |φ(u)| 随 |u| 增大单调递减 (不增)
    // 截断累积量 CF 仅在中等 |u| 范围内有效 (任务说明), 在 [0, 15] 上验证
    auto params = base_params();
    auto phi = make_rough_bergomi_cf(params);
    Real prev = 1.0;
    for (Real u = 0.05; u <= 15.0; u += 0.05) {
        Real modulus = std::abs(phi(Complex(u, 0.0)));
        EXPECT_LE(modulus, prev + 1e-12) << "|phi| not decreasing at u=" << u;
        prev = modulus;
    }
}

TEST(RoughBergomiCFTest, CFBoundedByOne) {
    // 有效特征函数须 |φ(u)| ≤ 1
    auto params = base_params();
    auto phi = make_rough_bergomi_cf(params);
    for (Real u = 0.0; u <= 10.0; u += 0.1) {
        EXPECT_LE(std::abs(phi(Complex(u, 0.0))), 1.0 + 1e-12)
            << "|phi| exceeded 1 at u=" << u;
    }
}

TEST(RoughBergomiCFTest, CFConjugateSymmetry) {
    // 实数分布性质: φ(-u) = conj(φ(u))
    auto params = base_params();
    auto phi = make_rough_bergomi_cf(params);
    for (Real u : {0.5, 1.0, 2.5, 5.0, 10.0}) {
        Complex phi_p = phi(Complex(u, 0.0));
        Complex phi_n = phi(Complex(-u, 0.0));
        EXPECT_NEAR(std::real(phi_n), std::real(phi_p), 1e-12) << "u=" << u;
        EXPECT_NEAR(std::imag(phi_n), -std::imag(phi_p), 1e-12) << "u=" << u;
    }
}

TEST(RoughBergomiCFTest, CFSmoothness) {
    // u ∈ [0, 10] 均匀采样, CF 实部/虚部无跳变 (导数连续)
    auto params = base_params();
    auto phi = make_rough_bergomi_cf(params);
    const Real du = 0.05;
    for (Real u = 0.0; u <= 10.0; u += du) {
        Complex a = phi(Complex(u, 0.0));
        Complex b = phi(Complex(u + du, 0.0));
        Real step = std::abs(b - a);
        EXPECT_LT(step, 0.5) << "CF jump at u=" << u << " step=" << step;
    }
}

TEST(RoughBergomiCFTest, EtaZeroDegeneration) {
    // η = 0 时, CF 退化为 GBM CF: φ(u) = exp(iu·ln S0 + iu·(r-q-0.5·xi0)T - 0.5·u²·xi0·T)
    auto params = base_params();
    params.eta = 0.0;
    auto phi = make_rough_bergomi_cf(params);
    auto gbm = make_gbm_cf(100.0, 0.03, 0.0, std::sqrt(0.04), 1.0);
    for (Real u : {0.1, 0.5, 1.0, 2.0, 5.0, 10.0}) {
        Complex a = phi(Complex(u, 0.0));
        Complex b = gbm(Complex(u, 0.0));
        EXPECT_NEAR(std::real(a), std::real(b), 1e-10) << "u=" << u;
        EXPECT_NEAR(std::imag(a), std::imag(b), 1e-10) << "u=" << u;
    }
}

TEST(RoughBergomiCFTest, EtaZeroMatchesBSM) {
    // η = 0 时 COS 定价应精确匹配 BSM (GBM, vol = sqrt(xi0))
    auto params = base_params();
    params.eta = 0.0;
    auto phi = make_rough_bergomi_cf(params);
    COSEngine::Config cfg;
    cfg.n_terms = 1024;
    cfg.L = 14.0;
    COSEngine engine(phi, 100.0, 0.03, 0.0, 1.0, cfg);
    Real cos_rb = engine.price_call(100.0);
    Real cos_gbm = cos_call_gbm(100.0, 100.0, 1.0, 0.03, 0.0,
                                std::sqrt(0.04), 1024, 14.0);
    EXPECT_NEAR(cos_rb, cos_gbm, 1e-9)
        << "rBergomi(eta=0)=" << cos_rb << " BSM=" << cos_gbm;
}

// ============ 2. 累积量公式 ============

TEST(RoughBergomiCFTest, CumulantCSign) {
    // c1 < 0 (对数收益率均值为负, 因 -0.5*∫v dt 凸性修正)
    auto params = base_params();
    auto cu = rough_bergomi_cumulants(params);
    EXPECT_LT(cu.c1, 0.0);
}

TEST(RoughBergomiCFTest, CumulantC2Positive) {
    // c2 > 0 (方差非负)
    auto params = base_params();
    auto cu = rough_bergomi_cumulants(params);
    EXPECT_GT(cu.c2, 0.0);
}

TEST(RoughBergomiCFTest, CumulantC3Sign) {
    // ρ < 0 → c3 < 0 (负偏度, lever effect); ρ > 0 → c3 > 0
    auto p_neg = base_params();
    auto p_pos = base_params();
    p_pos.rho = 0.7;
    EXPECT_LT(rough_bergomi_cumulants(p_neg).c3, 0.0);
    EXPECT_GT(rough_bergomi_cumulants(p_pos).c3, 0.0);
}

TEST(RoughBergomiCFTest, HHalfDegeneration) {
    // H → 0.5: c2 → xi0·T·(1 + η²·T^{2H}/(2·(2H+1))) = xi0·T·(1 + η²·T/4)
    auto params = base_params();
    params.H = 0.49999;
    Real c2 = rough_bergomi_cumulants(params).c2;
    Real expected = 0.04 * 1.0 * (1.0 + 0.3 * 0.3 * 1.0 / 4.0);
    EXPECT_NEAR(c2, expected, 1e-8);
}

TEST(RoughBergomiCFTest, Xi0Scaling) {
    // c2 ∝ xi0: 其他参数固定, xi0 翻倍 → c2 翻倍
    auto p1 = base_params();
    auto p2 = base_params();
    p2.xi0 = 0.08;
    Real c2_1 = rough_bergomi_cumulants(p1).c2;
    Real c2_2 = rough_bergomi_cumulants(p2).c2;
    EXPECT_NEAR(c2_2, 2.0 * c2_1, 1e-12);
}

TEST(RoughBergomiCFTest, CumulantToMoments) {
    // 累积量 → 矩 转换: mean=c1, var=c2, skew=c3/c2^{3/2}, kurt=c4/c2²
    auto params = base_params();
    auto cu = rough_bergomi_cumulants(params);
    auto mom = cumulants_to_moments(cu);
    EXPECT_NEAR(mom.mean, cu.c1, 1e-15);
    EXPECT_NEAR(mom.variance, cu.c2, 1e-15);
    EXPECT_NEAR(mom.skewness, cu.c3 / std::pow(cu.c2, 1.5), 1e-12);
    EXPECT_NEAR(mom.kurtosis_excess, cu.c4 / (cu.c2 * cu.c2), 1e-12);
    EXPECT_LT(mom.skewness, 0.0);       // ρ<0 → 负偏度
    EXPECT_GT(mom.kurtosis_excess, 0.0);  // log-vol 随机性 → 正峰度
}

// ============ 3. 参数校验 ============

TEST(RoughBergomiCFTest, ParameterValidation) {
    auto good = base_params();
    EXPECT_NO_THROW(validate_rough_bergomi_cf_params(good));

    auto bad = [](Real H, Real eta, Real rho, Real xi0, Real S0, Real T) {
        return RoughBergomiCFParams{H, eta, rho, xi0, S0, 0.03, 0.0, T};
    };
    // H 越界
    EXPECT_THROW(validate_rough_bergomi_cf_params(bad(0.0, 0.3, -0.7, 0.04, 100.0, 1.0)),
                 std::invalid_argument);
    EXPECT_THROW(validate_rough_bergomi_cf_params(bad(0.5, 0.3, -0.7, 0.04, 100.0, 1.0)),
                 std::invalid_argument);
    EXPECT_THROW(validate_rough_bergomi_cf_params(bad(-0.1, 0.3, -0.7, 0.04, 100.0, 1.0)),
                 std::invalid_argument);
    // eta: 负值非法, η=0 合法 (GBM 退化边界)
    EXPECT_THROW(validate_rough_bergomi_cf_params(bad(0.1, -1.0, -0.7, 0.04, 100.0, 1.0)),
                 std::invalid_argument);
    EXPECT_NO_THROW(validate_rough_bergomi_cf_params(bad(0.1, 0.0, -0.7, 0.04, 100.0, 1.0)));
    // rho 越界
    EXPECT_THROW(validate_rough_bergomi_cf_params(bad(0.1, 0.3, 1.5, 0.04, 100.0, 1.0)),
                 std::invalid_argument);
    EXPECT_THROW(validate_rough_bergomi_cf_params(bad(0.1, 0.3, -2.0, 0.04, 100.0, 1.0)),
                 std::invalid_argument);
    // xi0 / S0 非正
    EXPECT_THROW(validate_rough_bergomi_cf_params(bad(0.1, 0.3, -0.7, 0.0, 100.0, 1.0)),
                 std::invalid_argument);
    EXPECT_THROW(validate_rough_bergomi_cf_params(bad(0.1, 0.3, -0.7, 0.04, 0.0, 1.0)),
                 std::invalid_argument);
    // T 非正
    EXPECT_THROW(validate_rough_bergomi_cf_params(bad(0.1, 0.3, -0.7, 0.04, 100.0, 0.0)),
                 std::invalid_argument);
    EXPECT_THROW(validate_rough_bergomi_cf_params(bad(0.1, 0.3, -0.7, 0.04, 100.0, -1.0)),
                 std::invalid_argument);

    // 非法参数在 CF / 累积量入口同样抛异常
    auto invalid = bad(0.0, 0.3, -0.7, 0.04, 100.0, 1.0);
    EXPECT_THROW(make_rough_bergomi_cf(invalid), std::invalid_argument);
    EXPECT_THROW(rough_bergomi_cumulants(invalid), std::invalid_argument);
    EXPECT_THROW(rough_bergomi_characteristic_function(Complex(1.0, 0.0), invalid),
                 std::invalid_argument);
}

// ============ 4. COS 定价 ============

TEST(RoughBergomiCFTest, COSCallSanity) {
    // 价格合理性: 0 < call < S0, 0 < put < K
    auto params = base_params();
    auto phi = make_rough_bergomi_cf(params);
    COSEngine::Config cfg;
    cfg.n_terms = 1024;
    cfg.L = 14.0;
    COSEngine engine(phi, params.S0, params.r, params.q, params.T, cfg);
    for (Real K : {90.0, 100.0, 110.0}) {
        Real C = engine.price_call(K);
        Real P = engine.price_put(K);
        EXPECT_GT(C, 0.0) << "K=" << K;
        EXPECT_LT(C, params.S0) << "K=" << K;
        EXPECT_GT(P, 0.0) << "K=" << K;
        EXPECT_LT(P, K) << "K=" << K;
    }
}

TEST(RoughBergomiCFTest, CFCOSVsMC) {
    // COS 半解析定价 vs MC 定价, 容差 max(3*MC_SE, 5%)
    // 参数: S0=100, K=100, T=1, r=0.03, q=0, H=0.1, eta=0.3, rho=-0.7, xi0=0.04
    // 注意: 参考 MC (rbergomi_price_european) 的 log-Euler 无 (r-q) 漂移 (E[S_T]≈S_0),
    //       因此 COS 与 MC 在 r=q=0 下比较 (两者漂移均为 0, 模型一致).
    Real S0 = 100.0, K = 100.0, T = 1.0, H = 0.1, eta = 0.3, rho = -0.7, xi0 = 0.04;
    Real r = 0.0, q = 0.0;

    auto params = RoughBergomiCFParams{H, eta, rho, xi0, S0, r, q, T};
    auto phi = make_rough_bergomi_cf(params);
    COSEngine::Config cfg;
    cfg.n_terms = 1024;
    cfg.L = 14.0;
    COSEngine engine(phi, S0, r, q, T, cfg);
    Real cos_call = engine.price_call(K);

    RoughBergomiParams mp{H, eta, rho, xi0, S0, r, q};
    auto mc = rbergomi_price_european(mp, K, T, true, 5000, 12345, 64);
    Real tol = std::max(3.0 * mc.std_error, 0.05 * mc.price);

    EXPECT_NEAR(cos_call, mc.price, tol)
        << "COS=" << cos_call << " MC=" << mc.price << " MC_SE=" << mc.std_error
        << " tol=" << tol;
}

TEST(RoughBergomiCFTest, CallPutParity) {
    // 无套利 Call-Put 平价: C - P = S0·e^{-qT} - K·e^{-rT}
    // 用中等 vol-of-vol (η=0.1) 使截断累积量近似的伪尾部效应可忽略 (平价误差 ~2.5e-4)
    Real S0 = 100.0, T = 1.0, H = 0.1, eta = 0.1, rho = -0.7, xi0 = 0.04;
    Real r = 0.03, q = 0.0;

    auto params = RoughBergomiCFParams{H, eta, rho, xi0, S0, r, q, T};
    auto phi = make_rough_bergomi_cf(params);
    COSEngine::Config cfg;
    cfg.n_terms = 1024;
    cfg.L = 12.0;
    COSEngine engine(phi, S0, r, q, T, cfg);
    for (Real K : {90.0, 100.0, 110.0}) {
        Real C = engine.price_call(K);
        Real P = engine.price_put(K);
        Real parity = S0 * std::exp(-q * T) - K * std::exp(-r * T);
        EXPECT_NEAR(C - P, parity, 0.01)
            << "Parity failed at K=" << K << ": C-P=" << (C - P)
            << " rhs=" << parity;
    }
}

TEST(RoughBergomiCFTest, TZeroDegeneration) {
    // T → 0: call → max(S0 - K, 0), put → max(K - S0, 0)
    Real S0 = 100.0, T = 1e-3, H = 0.1, eta = 0.3, rho = -0.7, xi0 = 0.04;

    auto params = RoughBergomiCFParams{H, eta, rho, xi0, S0, 0.03, 0.0, T};
    auto phi = make_rough_bergomi_cf(params);
    COSEngine::Config cfg;
    cfg.n_terms = 1024;
    cfg.L = 10.0;
    COSEngine engine(phi, S0, 0.03, 0.0, T, cfg);
    // 小 T 下 ln S_T 近似退化分布, COS 数值近似 (容差 0.5)
    EXPECT_NEAR(engine.price_call(90.0), 10.0, 0.5);
    EXPECT_NEAR(engine.price_call(110.0), 0.0, 0.5);
    EXPECT_NEAR(engine.price_put(90.0), 0.0, 0.5);
    EXPECT_NEAR(engine.price_put(110.0), 10.0, 0.5);
}

// ============ 5. log(v_t) 矩 ============

TEST(RoughBergomiCFTest, LogVHelperClosedForm) {
    // rough_bergomi_log_v_moment 应匹配闭式:
    //   E[log v_t] = log(xi0) - 0.5·η²·t^{2H},  Var[log v_t] = η²·t^{2H}
    auto params = base_params();
    for (Real t : {0.1, 0.5, 1.0}) {
        auto lv = rough_bergomi_log_v_moment(t, params);
        Real t2H = std::pow(t, 2.0 * params.H);
        Real e_mean = std::log(params.xi0) - 0.5 * params.eta * params.eta * t2H;
        Real e_var = params.eta * params.eta * t2H;
        EXPECT_NEAR(lv.mean, e_mean, 1e-12) << "t=" << t;
        EXPECT_NEAR(lv.variance, e_var, 1e-12) << "t=" << t;
    }
}

TEST(RoughBergomiCFTest, LogVMoment) {
    // log(v_t) 的理论矩 vs 采样矩 (RLFbmSampler 路径), 容差 5%
    Real S0 = 100.0, T = 1.0, t = 1.0, H = 0.1, eta = 0.3, rho = -0.7, xi0 = 0.04;
    Size n_steps = 64, n_paths = 20000;

    RoughBergomiParams mp{H, eta, rho, xi0, S0, 0.03, 0.0};
    RoughBergomiProcess process(mp);
    RLFbmSampler sampler(T, n_steps, H);

    Real sum = 0.0, sum2 = 0.0;
    for (Size j = 0; j < n_paths; ++j) {
        Philox4x64 rng(7000 + j);
        std::vector<Real> vpath(n_steps + 1);
        process.generate_variance_path(T, n_steps, vpath, rng, sampler);
        Real lv = std::log(vpath.back());
        sum += lv;
        sum2 += lv * lv;
    }
    Real s_mean = sum / static_cast<Real>(n_paths);
    Real s_var = std::max(sum2 / static_cast<Real>(n_paths) - s_mean * s_mean, 0.0);

    // 均值: E[log v_t] = log(xi0) - 0.5·η²·t^{2H} (5%)
    Real e_mean = std::log(xi0) - 0.5 * eta * eta * std::pow(t, 2.0 * H);
    EXPECT_NEAR(s_mean, e_mean, 0.05 * std::abs(e_mean))
        << "sample mean=" << s_mean << " theory=" << e_mean;

    // 方差: 参考采样器按 sqrt(2H+1) 归一化 (rough_bergomi.hpp), 实际 Var(W̃_t) = C_{nn},
    //       故样本 Var(log v_t) ≈ η²·C_{nn}. 校验采样器自洽性 (5%).
    auto C = rbergomi_fbm_covariance(T, n_steps, H);
    Real implied_var = eta * eta * C[n_steps - 1][n_steps - 1];
    EXPECT_NEAR(s_var, implied_var, 0.05 * implied_var)
        << "sample var=" << s_var << " sampler-implied=" << implied_var;
}
