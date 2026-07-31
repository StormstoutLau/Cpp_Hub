// v1.2 Batch 8: Copula 模型 (Gaussian / One-Factor / t) 单元测试
// 覆盖: GaussianCopula / OneFactorGaussianCopula / TCopula
//       + cholesky_dynamic + make_equicorrelation
// 测试维度: 配置验证/Cholesky正确性/均匀变量转换/违约时间采样/
//          条件违约概率/尾部相关性/MC统计性质
#include <gtest/gtest.h>
#include "cpphub/instruments/credit/copula.hpp"
#include "cpphub/instruments/credit/credit_curve.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/core/math.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace cpphub;

namespace {
// 构造 n×n 单位矩阵 (独立 Copula 等价)
std::vector<std::vector<Real>> identity_matrix(Size n) {
    std::vector<std::vector<Real>> m(n, std::vector<Real>(n, 0.0));
    for (Size i = 0; i < n; ++i) m[i][i] = 1.0;
    return m;
}

// 构造等相关矩阵: Σ_ii = 1, Σ_ij = ρ (i≠j)
std::vector<std::vector<Real>> equicorr(Size n, Real rho) {
    std::vector<std::vector<Real>> m(n, std::vector<Real>(n, rho));
    for (Size i = 0; i < n; ++i) m[i][i] = 1.0;
    return m;
}
}  // namespace

// ============================================================
// 1. cholesky_dynamic 基础测试
// ============================================================
TEST(CholeskyDynamicTest, IdentityMatrix) {
    auto A = identity_matrix(3);
    auto L = cholesky_dynamic(A);
    // 单位矩阵的 Cholesky = 单位矩阵
    for (Size i = 0; i < 3; ++i) {
        for (Size j = 0; j < 3; ++j) {
            Real expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(L[i][j], expected, 1e-12);
        }
    }
}

TEST(CholeskyDynamicTest, ReconstructsOriginalMatrix) {
    // 3×3 正定矩阵
    std::vector<std::vector<Real>> A = {
        {1.0, 0.5, 0.3},
        {0.5, 1.0, 0.4},
        {0.3, 0.4, 1.0}
    };
    auto L = cholesky_dynamic(A);
    // 验证 L * L^T = A
    for (Size i = 0; i < 3; ++i) {
        for (Size j = 0; j < 3; ++j) {
            Real sum = 0.0;
            for (Size k = 0; k < 3; ++k) sum += L[i][k] * L[j][k];
            EXPECT_NEAR(sum, A[i][j], 1e-10);
        }
    }
}

TEST(CholeskyDynamicTest, ThrowsOnNonPositiveSemidefinite) {
    std::vector<std::vector<Real>> A = {
        {1.0, 2.0},
        {2.0, 1.0}  // 特征值 -1, 3 (非半正定)
    };
    EXPECT_THROW(cholesky_dynamic(A), std::invalid_argument);
}

TEST(CholeskyDynamicTest, HandlesSingularMatrix) {
    // ρ=1 的 2×2 等相关矩阵 (半正定, 特征值 0, 2)
    auto A = equicorr(2, 1.0);
    EXPECT_NO_THROW(cholesky_dynamic(A));
    auto L = cholesky_dynamic(A);
    // L 应为 [[1, 0], [1, 0]] (退化情形)
    EXPECT_NEAR(L[0][0], 1.0, 1e-10);
    EXPECT_NEAR(L[1][0], 1.0, 1e-10);
    EXPECT_NEAR(L[1][1], 0.0, 1e-10);
}

// ============================================================
// 2. GaussianCopula 配置与验证
// ============================================================
TEST(GaussianCopulaTest, ConstructorValidates) {
    // 空矩阵
    std::vector<std::vector<Real>> empty_mat;
    EXPECT_THROW(GaussianCopula{empty_mat}, std::invalid_argument);

    // 非方阵
    std::vector<std::vector<Real>> non_square = {{1.0, 0.5}, {0.5, 1.0, 0.3}};
    EXPECT_THROW(GaussianCopula{non_square}, std::invalid_argument);

    // 对角线不为 1
    std::vector<std::vector<Real>> bad_diag = {{0.5, 0.3}, {0.3, 1.0}};
    EXPECT_THROW(GaussianCopula{bad_diag}, std::invalid_argument);

    // 不对称
    std::vector<std::vector<Real>> asym = {{1.0, 0.5}, {0.3, 1.0}};
    EXPECT_THROW(GaussianCopula{asym}, std::invalid_argument);

    // |相关系数| > 1
    std::vector<std::vector<Real>> big_corr = {{1.0, 1.5}, {1.5, 1.0}};
    EXPECT_THROW(GaussianCopula{big_corr}, std::invalid_argument);
}

TEST(GaussianCopulaTest, NNamesAndAccessors) {
    GaussianCopula gc(equicorr(4, 0.3));
    EXPECT_EQ(gc.n_names(), 4u);
    EXPECT_EQ(gc.correlation().size(), 4u);
    EXPECT_EQ(gc.cholesky_L().size(), 4u);
    EXPECT_NEAR(gc.correlation()[0][1], 0.3, 1e-15);
    EXPECT_NEAR(gc.correlation()[1][0], 0.3, 1e-15);
    EXPECT_NEAR(gc.correlation()[0][0], 1.0, 1e-15);
}

// ============================================================
// 3. GaussianCopula 均匀变量转换
// ============================================================
TEST(GaussianCopulaTest, TransformUniformsOutputInRange) {
    GaussianCopula gc(equicorr(3, 0.5));
    std::vector<Real> u_in = {0.3, 0.6, 0.8};
    auto u_out = gc.transform_uniforms(u_in);
    ASSERT_EQ(u_out.size(), 3u);
    for (Real u : u_out) {
        EXPECT_GT(u, 0.0);
        EXPECT_LT(u, 1.0);
    }
}

TEST(GaussianCopulaTest, TransformUniformsSizeMismatch) {
    GaussianCopula gc(equicorr(3, 0.5));
    std::vector<Real> u_wrong_size = {0.3, 0.6};
    EXPECT_THROW(gc.transform_uniforms(u_wrong_size), std::invalid_argument);
}

TEST(GaussianCopulaTest, IdentityCorrelationPreservesUniforms) {
    // ρ=0 (单位矩阵): 输出的相关均匀变量应等于输入 (容差范围内)
    GaussianCopula gc(identity_matrix(3));
    std::vector<Real> u_in = {0.3, 0.6, 0.8};
    auto u_out = gc.transform_uniforms(u_in);
    // ρ=0 时, L = I, X = Z = Φ^{-1}(u), U = Φ(X) = u
    for (Size i = 0; i < 3; ++i) {
        EXPECT_NEAR(u_out[i], u_in[i], 1e-10);
    }
}

TEST(GaussianCopulaTest, PerfectCorrelationProducesIdenticalUniforms) {
    // ρ=1: 所有名字的均匀变量应相同
    GaussianCopula gc(equicorr(3, 1.0));
    std::vector<Real> u_in = {0.3, 0.6, 0.8};
    auto u_out = gc.transform_uniforms(u_in);
    // ρ=1 时, 所有 X_i 相同 → 所有 U_i 相同
    // 注意: 由于 Cholesky 在半正定情形下的实现, 输出取决于 L 的第一列
    Real u0 = u_out[0];
    for (Size i = 1; i < 3; ++i) {
        EXPECT_NEAR(u_out[i], u0, 1e-10);
    }
}

TEST(GaussianCopulaTest, MonotonicTransform) {
    // 固定其他输入, 单个输入增大时, 输出也应单调递增
    GaussianCopula gc(equicorr(2, 0.5));
    Real u2_fixed = 0.5;
    Real prev = -1.0;
    for (Size k = 0; k <= 10; ++k) {
        Real u1 = 0.05 + 0.09 * k;  // 0.05 → 0.95
        auto u_out = gc.transform_uniforms({u1, u2_fixed});
        EXPECT_GT(u_out[0], prev);
        prev = u_out[0];
    }
}

// ============================================================
// 4. GaussianCopula MC 统计性质
// ============================================================
TEST(GaussianCopulaTest, SampleUniformsMeanApproxHalf) {
    // 大量采样后, 每个名字的均匀变量均值应接近 0.5
    GaussianCopula gc(equicorr(3, 0.3));
    Philox4x64 rng(42);
    const Size n_paths = 50000;
    std::vector<Real> sum(3, 0.0);
    for (Size p = 0; p < n_paths; ++p) {
        auto u = gc.sample_uniforms(rng);
        for (Size i = 0; i < 3; ++i) sum[i] += u[i];
    }
    for (Size i = 0; i < 3; ++i) {
        Real mean = sum[i] / static_cast<Real>(n_paths);
        EXPECT_NEAR(mean, 0.5, 0.02);  // MC 误差 ~ 1/sqrt(N)
    }
}

TEST(GaussianCopulaTest, SampleUniformsCorrelationMatchesInput) {
    // 采样相关均匀变量, 转换回正态后, 经验相关系数应接近输入 ρ
    Real rho_input = 0.4;
    GaussianCopula gc(equicorr(2, rho_input));
    Philox4x64 rng(42);
    const Size n_paths = 50000;
    Real sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0, sum_yy = 0.0;
    for (Size p = 0; p < n_paths; ++p) {
        auto u = gc.sample_uniforms(rng);
        // 转换回正态
        Real x = inv_normal_cdf(u[0]);
        Real y = inv_normal_cdf(u[1]);
        sum_x += x; sum_y += y;
        sum_xy += x * y;
        sum_xx += x * x; sum_yy += y * y;
    }
    Real n = static_cast<Real>(n_paths);
    Real mean_x = sum_x / n, mean_y = sum_y / n;
    Real cov = sum_xy / n - mean_x * mean_y;
    Real var_x = sum_xx / n - mean_x * mean_x;
    Real var_y = sum_yy / n - mean_y * mean_y;
    Real corr = cov / std::sqrt(var_x * var_y);
    EXPECT_NEAR(corr, rho_input, 0.02);
}

// ============================================================
// 5. GaussianCopula 违约时间采样
// ============================================================
TEST(GaussianCopulaTest, UniformsToDefaultTimesMonotonic) {
    // U = PD(0, τ) = 1 - exp(-H(τ)), U 越大 → H(τ) 越大 → τ 越大 (越晚违约)
    // 因为 H(τ) 单调递增, 所以 U 与 τ 同向单调
    auto cc = CreditCurve::flat(0.05, 0.4);  // 5% hazard
    GaussianCopula gc(identity_matrix(1));
    std::vector<CreditCurve> curves = {cc};

    Real prev_tau = 0.0;
    for (Size k = 0; k <= 10; ++k) {
        Real u = 0.05 + 0.09 * k;  // 0.05 → 0.95
        std::vector<Real> taus = gc.uniforms_to_default_times({u}, curves);
        // τ 应单调递增
        EXPECT_GE(taus[0], prev_tau);
        prev_tau = taus[0];
    }
}

TEST(GaussianCopulaTest, UniformsToDefaultTimesMatchesFlatHazard) {
    // flat hazard h: τ = -ln(1-U) / h
    Real h = 0.05;
    auto cc = CreditCurve::flat(h, 0.4);
    GaussianCopula gc(identity_matrix(1));
    std::vector<CreditCurve> curves = {cc};

    Real u = 0.3;
    Real expected_tau = -std::log(1.0 - u) / h;
    auto taus = gc.uniforms_to_default_times({u}, curves);
    EXPECT_NEAR(taus[0], expected_tau, 1e-6);
}

TEST(GaussianCopulaTest, UniformsToDefaultTimesSizeMismatch) {
    GaussianCopula gc(identity_matrix(2));
    auto cc = CreditCurve::flat(0.05, 0.4);
    std::vector<CreditCurve> curves = {cc};  // 只有 1 条
    EXPECT_THROW(gc.uniforms_to_default_times({0.3, 0.5}, curves),
                  std::invalid_argument);
}

TEST(GaussianCopulaTest, SampleDefaultTimesProducesFiniteTimes) {
    // 中等 hazard 下, 大多数采样违约时间应有限
    auto cc = CreditCurve::flat(0.1, 0.4);
    GaussianCopula gc(equicorr(3, 0.3));
    std::vector<CreditCurve> curves = {cc, cc, cc};
    Philox4x64 rng(42);
    const Size n_paths = 1000;
    Size finite_count = 0;
    for (Size p = 0; p < n_paths; ++p) {
        auto taus = gc.sample_default_times(rng, curves);
        for (Real t : taus) {
            if (std::isfinite(t)) ++finite_count;
        }
    }
    // 至少 99% 应有限
    EXPECT_GT(finite_count, static_cast<Size>(0.99 * n_paths * 3));
}

// ============================================================
// 6. OneFactorGaussianCopula
// ============================================================
TEST(OneFactorGaussianCopulaTest, ConstructorValidates) {
    EXPECT_THROW(OneFactorGaussianCopula(-0.1, 5), std::invalid_argument);
    EXPECT_THROW(OneFactorGaussianCopula(1.0, 5), std::invalid_argument);
    EXPECT_THROW(OneFactorGaussianCopula(1.5, 5), std::invalid_argument);
    EXPECT_THROW(OneFactorGaussianCopula(0.3, 0), std::invalid_argument);
}

TEST(OneFactorGaussianCopulaTest, Accessors) {
    OneFactorGaussianCopula ofgc(0.3, 5);
    EXPECT_NEAR(ofgc.rho(), 0.3, 1e-15);
    EXPECT_EQ(ofgc.n_names(), 5u);
}

TEST(OneFactorGaussianCopulaTest, CorrelationMatrixMatchesParameter) {
    OneFactorGaussianCopula ofgc(0.4, 4);
    auto corr = ofgc.correlation_matrix();
    ASSERT_EQ(corr.size(), 4u);
    for (Size i = 0; i < 4; ++i) {
        EXPECT_NEAR(corr[i][i], 1.0, 1e-15);
        for (Size j = 0; j < 4; ++j) {
            if (i != j) EXPECT_NEAR(corr[i][j], 0.4, 1e-15);
        }
    }
}

TEST(OneFactorGaussianCopulaTest, TransformUniformsSizeMismatch) {
    OneFactorGaussianCopula ofgc(0.3, 3);
    EXPECT_THROW(ofgc.transform_uniforms(0.5, {0.2, 0.4}),
                  std::invalid_argument);
}

TEST(OneFactorGaussianCopulaTest, ZeroRhoProducesIndependentUniforms) {
    // ρ=0: X_i = Z_i (完全独立), U_i = Φ(Z_i) 应与输入相同
    OneFactorGaussianCopula ofgc(0.0, 3);
    Real u_sys = 0.5;
    std::vector<Real> u_id = {0.3, 0.6, 0.8};
    auto u_out = ofgc.transform_uniforms(u_sys, u_id);
    // ρ=0 时, U_i 与 u_sys 无关, 应等于 u_id[i]
    for (Size i = 0; i < 3; ++i) {
        EXPECT_NEAR(u_out[i], u_id[i], 1e-10);
    }
}

TEST(OneFactorGaussianCopulaTest, ConditionalPDMatchesFormula) {
    // 验证 conditional_pd 与手算公式一致
    // PD_i(t|m) = Φ((Φ^{-1}(PD) - √ρ*m) / √(1-ρ))
    OneFactorGaussianCopula ofgc(0.3, 5);
    Real pd = 0.1;
    Real m = 0.5;
    Real expected = normal_cdf((inv_normal_cdf(pd) - std::sqrt(0.3) * m)
                                 / std::sqrt(1.0 - 0.3));
    Real actual = ofgc.conditional_pd(pd, m);
    EXPECT_NEAR(actual, expected, 1e-12);
}

TEST(OneFactorGaussianCopulaTest, ConditionalPDMonotonicInM) {
    // 给定 PD, m 越大 (系统因子越好), 条件违约概率应越小
    OneFactorGaussianCopula ofgc(0.3, 5);
    Real pd = 0.1;
    Real prev = 1.0;
    for (Size k = 0; k <= 10; ++k) {
        Real m = -2.0 + 0.4 * k;  // -2.0 → 2.0
        Real p = ofgc.conditional_pd(pd, m);
        EXPECT_LE(p, prev);
        prev = p;
    }
}

TEST(OneFactorGaussianCopulaTest, ConditionalPDIncreasesWithPD) {
    // 给定 m, PD 越大, 条件 PD 越大
    OneFactorGaussianCopula ofgc(0.3, 5);
    Real m = 0.0;
    Real prev = -1.0;
    for (Size k = 0; k <= 10; ++k) {
        Real pd = 0.01 + 0.09 * k;
        Real p = ofgc.conditional_pd(pd, m);
        EXPECT_GT(p, prev);
        prev = p;
    }
}

TEST(OneFactorGaussianCopulaTest, SampleUniformsMeanApproxHalf) {
    OneFactorGaussianCopula ofgc(0.3, 3);
    Philox4x64 rng(42);
    const Size n_paths = 50000;
    std::vector<Real> sum(3, 0.0);
    for (Size p = 0; p < n_paths; ++p) {
        auto u = ofgc.sample_uniforms(rng);
        for (Size i = 0; i < 3; ++i) sum[i] += u[i];
    }
    for (Size i = 0; i < 3; ++i) {
        Real mean = sum[i] / static_cast<Real>(n_paths);
        EXPECT_NEAR(mean, 0.5, 0.02);
    }
}

// ============================================================
// 7. TCopula
// ============================================================
TEST(TCopulaTest, ConstructorValidates) {
    // 自由度 ≤ 2 应拒绝
    EXPECT_THROW(TCopula(equicorr(3, 0.5), 1.0), std::invalid_argument);
    EXPECT_THROW(TCopula(equicorr(3, 0.5), 2.0), std::invalid_argument);
    // 自由度 > 2 应通过
    EXPECT_NO_THROW(TCopula(equicorr(3, 0.5), 3.0));
    EXPECT_NO_THROW(TCopula(equicorr(3, 0.5), 10.0));
}

TEST(TCopulaTest, Accessors) {
    TCopula tc(equicorr(4, 0.4), 5.0);
    EXPECT_EQ(tc.n_names(), 4u);
    EXPECT_NEAR(tc.degrees_of_freedom(), 5.0, 1e-15);
    EXPECT_EQ(tc.correlation().size(), 4u);
}

TEST(TCopulaTest, TransformUniformsOutputInRange) {
    TCopula tc(equicorr(3, 0.5), 5.0);
    std::vector<Real> u_in = {0.3, 0.6, 0.8};
    Real w = 5.0;  // χ²(5) 的均值
    auto u_out = tc.transform_uniforms(u_in, w);
    ASSERT_EQ(u_out.size(), 3u);
    for (Real u : u_out) {
        EXPECT_GT(u, 0.0);
        EXPECT_LT(u, 1.0);
    }
}

TEST(TCopulaTest, TransformUniformsSizeMismatch) {
    TCopula tc(equicorr(3, 0.5), 5.0);
    EXPECT_THROW(tc.transform_uniforms({0.3, 0.5}, 5.0),
                  std::invalid_argument);
}

TEST(TCopulaTest, TransformUniformsRejectsNonPositiveW) {
    TCopula tc(equicorr(3, 0.5), 5.0);
    std::vector<Real> u_in = {0.3, 0.6, 0.8};
    EXPECT_THROW(tc.transform_uniforms(u_in, 0.0), std::invalid_argument);
    EXPECT_THROW(tc.transform_uniforms(u_in, -1.0), std::invalid_argument);
}

TEST(TCopulaTest, SampleUniformsMeanApproxHalf) {
    TCopula tc(equicorr(3, 0.3), 5.0);
    Philox4x64 rng(42);
    const Size n_paths = 50000;
    std::vector<Real> sum(3, 0.0);
    for (Size p = 0; p < n_paths; ++p) {
        auto u = tc.sample_uniforms(rng);
        for (Size i = 0; i < 3; ++i) sum[i] += u[i];
    }
    for (Size i = 0; i < 3; ++i) {
        Real mean = sum[i] / static_cast<Real>(n_paths);
        EXPECT_NEAR(mean, 0.5, 0.03);
    }
}

TEST(TCopulaTest, LargeNuApproachesGaussian) {
    // ν → ∞ 时, t-Copula 退化为 Gaussian Copula
    // 用很大的 ν 检验: 输出应接近 Gaussian Copula 的输出
    Real rho = 0.4;
    GaussianCopula gc(equicorr(2, rho));
    TCopula tc(equicorr(2, rho), 1000.0);  // ν 很大
    std::vector<Real> u_in = {0.3, 0.7};
    Real w = 1000.0;  // ≈ ν (均值)
    auto u_gauss = gc.transform_uniforms(u_in);
    auto u_t = tc.transform_uniforms(u_in, w);
    // 应非常接近 (差异来自 1/√(W/ν) ≈ 1)
    EXPECT_NEAR(u_t[0], u_gauss[0], 0.05);
    EXPECT_NEAR(u_t[1], u_gauss[1], 0.05);
}

TEST(TCopulaTest, UniformsToDefaultTimesMonotonic) {
    // U 越大 → 违约时间越晚 (同 GaussianCopula, 见上文数学说明)
    Real h = 0.05;
    auto cc = CreditCurve::flat(h, 0.4);
    TCopula tc(identity_matrix(1), 5.0);
    std::vector<CreditCurve> curves = {cc};
    Real prev = 0.0;
    for (Size k = 0; k <= 10; ++k) {
        Real u = 0.05 + 0.09 * k;
        auto taus = tc.uniforms_to_default_times({u}, curves);
        EXPECT_GE(taus[0], prev);
        prev = taus[0];
    }
}

// ============================================================
// 8. 便捷工厂函数
// ============================================================
TEST(CopulaFactoryTest, MakeEquicorrelation) {
    auto corr = make_equicorrelation(4, 0.3);
    ASSERT_EQ(corr.size(), 4u);
    for (Size i = 0; i < 4; ++i) {
        EXPECT_NEAR(corr[i][i], 1.0, 1e-15);
        for (Size j = 0; j < 4; ++j) {
            if (i != j) EXPECT_NEAR(corr[i][j], 0.3, 1e-15);
        }
    }
}

TEST(CopulaFactoryTest, MakeGaussianCopula) {
    auto corr = make_equicorrelation(3, 0.5);
    auto gc = make_gaussian_copula(corr);
    EXPECT_EQ(gc.n_names(), 3u);
}

TEST(CopulaFactoryTest, MakeTCopula) {
    auto corr = make_equicorrelation(3, 0.5);
    auto tc = make_t_copula(corr, 5.0);
    EXPECT_EQ(tc.n_names(), 3u);
    EXPECT_NEAR(tc.degrees_of_freedom(), 5.0, 1e-15);
}

TEST(CopulaFactoryTest, MakeEquicorrelationInvalidRho) {
    // ρ ≥ 1 非法
    EXPECT_THROW(make_equicorrelation(3, 1.0), std::invalid_argument);
    EXPECT_THROW(make_equicorrelation(3, 1.5), std::invalid_argument);
}
