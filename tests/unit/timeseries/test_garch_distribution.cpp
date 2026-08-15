// =============================================================================
// test_garch_distribution.cpp - GARCH 分布族似然函数测试 (12 用例)
//
// Phase 7B v1.6 M1 (PHASE7B_FINANCIAL_TS_SPEC.md §2.0.1 测试矩阵)
//
// 幻觉点覆盖:
//   G3: Normal 似然含 -0.5·log(2π) 常数项 (硬编码手算值)
//   G14: t 分布似然公式 (Bollerslev 1987)
//   G-GED1: GED 指数 ν + 缩放常数 c 含 2^(-1/ν) (arch 源码实测公式)
//   ν→∞ 退化 / GED ν=2 退化 / 批量求和一致性 / 异常输入
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <stdexcept>

#include "cpphub/timeseries/garch/garch_distribution.hpp"

using namespace cpphub::v1::timeseries::garch;
using cpphub::v1::Real;
using cpphub::v1::Size;

// 独立重算 (转录核查: 测试内独立写公式, 与实现解耦)
namespace ref {

constexpr double kPi = 3.14159265358979323846;

double normal_ll(double eps, double h) {
    return -0.5 * (std::log(2.0 * kPi) + std::log(h) + eps * eps / h);
}

double t_ll(double eps, double h, double nu) {
    const double hp = 0.5 * (nu + 1.0);
    const double u = eps * eps / ((nu - 2.0) * h);
    return std::lgamma(hp) - std::lgamma(0.5 * nu) - 0.5 * std::log(nu - 2.0)
           - 0.5 * std::log(kPi) - 0.5 * std::log(h) - hp * std::log(1.0 + u);
}

double ged_c(double nu) {
    return std::pow(2.0, -1.0 / nu)
           * std::sqrt(std::tgamma(1.0 / nu) / std::tgamma(3.0 / nu));
}

double ged_ll(double eps, double h, double nu) {
    const double c = ged_c(nu);
    const double a = std::abs(eps) / (std::sqrt(h) * c);
    return std::log(nu) - std::log(c) - std::lgamma(1.0 / nu)
           - (1.0 + 1.0 / nu) * std::log(2.0) - 0.5 * std::log(h)
           - 0.5 * std::pow(a, nu);
}

}  // namespace ref

// ---------------------------------------------------------------------------
// 1. Normal 单观测似然 vs 手算 (ε=0.1, h=0.04) — G3 常数项
// ---------------------------------------------------------------------------
TEST(GarchDistribution, NormalSingleTermHandComputed) {
    // 手算: -0.5·[log(2π) + log(0.04) + 0.01/0.04]
    //     = -0.5·[1.8378770664093453 - 3.2188758248682006 + 0.25]
    //     = 0.56549937922942765
    const Real ll = log_likelihood_term(0.1, 0.04, GarchDist::Normal);
    EXPECT_NEAR(ll, 0.56549937922942765, 1e-15);
    EXPECT_NEAR(ll, ref::normal_ll(0.1, 0.04), 1e-15);
}

// ---------------------------------------------------------------------------
// 2. Normal ε=0, h=1 → -0.5·log(2π) (纯常数项验证)
// ---------------------------------------------------------------------------
TEST(GarchDistribution, NormalZeroResidualUnitVariance) {
    const Real ll = log_likelihood_term(0.0, 1.0, GarchDist::Normal);
    EXPECT_NEAR(ll, -0.5 * std::log(2.0 * ref::kPi), 1e-15);
    EXPECT_NEAR(ll, -0.91893853320467274, 1e-15);
}

// ---------------------------------------------------------------------------
// 3. Student-t ν=5 单观测 vs 手算 (G14 公式)
// ---------------------------------------------------------------------------
TEST(GarchDistribution, StudentTSingleTermHandComputed) {
    // 手算 (ν=5, ε=0.1, h=0.04):
    //   lgamma(3) - lgamma(2.5) - 0.5·log(3) - 0.5·log(π) - 0.5·log(0.04)
    //   - 3·log(1 + 0.01/0.12) = 0.6561030122417623
    const Real ll = log_likelihood_term(0.1, 0.04, GarchDist::StudentT, 5.0);
    EXPECT_NEAR(ll, 0.6561030122417623, 1e-12);
    EXPECT_NEAR(ll, ref::t_ll(0.1, 0.04, 5.0), 1e-15);
}

// ---------------------------------------------------------------------------
// 4. t 分布 ν→∞ 退化到 Normal (渐近行为, 收敛率 O(1/ν))
// ---------------------------------------------------------------------------
TEST(GarchDistribution, StudentTLargeNuDegeneratesToNormal) {
    // 排幻觉修正: t(ν)→Normal 收敛率为 O(1/ν), ν=1000 时残差 ~1e-3 达不到
    // 1e-6; 用 ν=1e8 验证 spec 的 1e-6 容差 (渐近极限的数学事实)
    const Real ll_t = log_likelihood_term(0.3, 0.02, GarchDist::StudentT, 1e8);
    const Real ll_n = log_likelihood_term(0.3, 0.02, GarchDist::Normal);
    EXPECT_NEAR(ll_t, ll_n, 1e-6);
}

// ---------------------------------------------------------------------------
// 5. GED ν=2 (正态特殊情况) 退化到 Normal
// ---------------------------------------------------------------------------
TEST(GarchDistribution, GEDNu2EqualsNormal) {
    for (double eps : {-1.2, -0.3, 0.0, 0.7, 2.5}) {
        for (double h : {0.01, 0.5, 3.0}) {
            EXPECT_NEAR(log_likelihood_term(eps, h, GarchDist::GED, 2.0),
                        log_likelihood_term(eps, h, GarchDist::Normal), 1e-10)
                << "eps=" << eps << " h=" << h;
        }
    }
}

// ---------------------------------------------------------------------------
// 6. GED ε=0 时似然 = 纯常数项 (G-GED1 常数项形态手算)
// ---------------------------------------------------------------------------
TEST(GarchDistribution, GEDZeroResidualConstantTerm) {
    // ν=2.5, ε=0, h=1: ℓ = log(ν) - log(c) - log(Γ(1/ν)) - (1+1/ν)·log(2)
    const Real nu = 2.5;
    const Real c = ref::ged_c(nu);
    const Real expected = std::log(nu) - std::log(c) - std::lgamma(1.0 / nu)
                          - (1.0 + 1.0 / nu) * std::log(2.0);
    EXPECT_NEAR(log_likelihood_term(0.0, 1.0, GarchDist::GED, nu), expected, 1e-14);
}

// ---------------------------------------------------------------------------
// 6b. GED/t ν 低于 arch 下界 2.05 → 抛异常 (边界行为)
// ---------------------------------------------------------------------------
TEST(GarchDistribution, GEDNuBelowArchBoundThrows) {
    EXPECT_THROW(log_likelihood_term(0.0, 1.0, GarchDist::GED, 1.0),
                 std::invalid_argument);
    EXPECT_THROW(log_likelihood_term(0.0, 1.0, GarchDist::StudentT, 2.0),
                 std::invalid_argument);
    EXPECT_NO_THROW(log_likelihood_term(0.0, 1.0, GarchDist::GED, 2.051));
}

// ---------------------------------------------------------------------------
// 7. GED ν=2.5 单观测 vs 独立公式重算 (G-GED1 公式形态)
// ---------------------------------------------------------------------------
TEST(GarchDistribution, GEDMatchesReferenceFormula) {
    for (double eps : {-0.8, 0.0, 1.5}) {
        for (double h : {0.1, 2.0}) {
            EXPECT_NEAR(log_likelihood_term(eps, h, GarchDist::GED, 2.5),
                        ref::ged_ll(eps, h, 2.5), 1e-14)
                << "eps=" << eps << " h=" << h;
        }
    }
    // 手算锚点 (ν=2.5, ε=1.5, h=1):
    const Real c = ref::ged_c(2.5);
    const Real expected = std::log(2.5) - std::log(c) - std::lgamma(0.4)
                          - 1.4 * std::log(2.0) - 0.5 * std::pow(1.5 / c, 2.5);
    EXPECT_NEAR(log_likelihood_term(1.5, 1.0, GarchDist::GED, 2.5), expected, 1e-14);
}

// ---------------------------------------------------------------------------
// 8. 批量似然 = 逐项求和 (1e-15)
// ---------------------------------------------------------------------------
TEST(GarchDistribution, BatchEqualsSumOfTerms) {
    const std::vector<Real> eps = {0.01, -0.02, 0.03, -0.01, 0.005, -0.04};
    const std::vector<Real> h = {0.0001, 0.00011, 0.00012, 0.00013, 0.00014, 0.00015};
    Real manual = 0.0;
    for (Size i = 0; i < eps.size(); ++i) {
        manual += log_likelihood_term(eps[i], h[i], GarchDist::StudentT, 7.0);
    }
    EXPECT_NEAR(log_likelihood(eps, h, GarchDist::StudentT, 7.0), manual, 1e-15);
}

// ---------------------------------------------------------------------------
// 9. 长度不一致 / variance<=0 → 抛异常
// ---------------------------------------------------------------------------
TEST(GarchDistribution, InvalidInputsThrow) {
    const std::vector<Real> eps = {0.01, 0.02};
    const std::vector<Real> h_bad = {0.01};
    EXPECT_THROW(log_likelihood(eps, h_bad, GarchDist::Normal),
                 std::invalid_argument);
    EXPECT_THROW(log_likelihood_term(0.1, 0.0, GarchDist::Normal),
                 std::invalid_argument);
    EXPECT_THROW(log_likelihood_term(0.1, -0.01, GarchDist::Normal),
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 10. t 分布 ν 增大 → 似然单调趋近 Normal (ε²/h 较大处差异明显)
// ---------------------------------------------------------------------------
TEST(GarchDistribution, StudentTConvergesMonotonically) {
    const Real ll_n = ref::normal_ll(2.0, 1.0);
    Real prev_diff = std::abs(ref::t_ll(2.0, 1.0, 5.0) - ll_n);
    for (double nu : {10.0, 50.0, 200.0, 1000.0}) {
        const Real diff = std::abs(ref::t_ll(2.0, 1.0, nu) - ll_n);
        EXPECT_LT(diff, prev_diff + 1e-12);
        prev_diff = diff;
    }
    // O(1/ν) 收敛率: ν=1e8 时残差 < 1e-6 (排幻觉: ν=1000 仅 ~1e-3)
    EXPECT_NEAR(ref::t_ll(2.0, 1.0, 1e8), ll_n, 1e-6);
}

// ---------------------------------------------------------------------------
// 11. 梯度: dℓ/dν (t 解析) vs 数值中心差分
// ---------------------------------------------------------------------------
TEST(GarchDistribution, StudentTGradientMatchesNumerical) {
    const std::vector<Real> eps = {0.5, -1.0, 0.2, -0.3, 1.5, -0.8, 0.0, 0.9};
    std::vector<Real> h(eps.size(), 1.0);
    const Real nu = 8.0;
    const auto grad = log_likelihood_gradient(eps, h, GarchDist::StudentT, nu);

    // 前 3 分量 (ω,α,β) 恒为 0 (静态序列无递归路径, 见头文件说明)
    EXPECT_EQ(grad.size(), 4u);
    EXPECT_EQ(grad[0], 0.0);
    EXPECT_EQ(grad[1], 0.0);
    EXPECT_EQ(grad[2], 0.0);

    const Real delta = 1e-6;
    const Real num = (log_likelihood(eps, h, GarchDist::StudentT, nu + delta)
                      - log_likelihood(eps, h, GarchDist::StudentT, nu - delta))
                     / (2.0 * delta);
    EXPECT_NEAR(grad[3], num, 1e-6);
}

// ---------------------------------------------------------------------------
// 12. 梯度: Normal → 全 0; GED → 数值差分
// ---------------------------------------------------------------------------
TEST(GarchDistribution, GradientNormalAndGED) {
    const std::vector<Real> eps = {0.4, -0.9, 1.1, -0.2};
    std::vector<Real> h(eps.size(), 0.5);

    const auto g_n = log_likelihood_gradient(eps, h, GarchDist::Normal);
    EXPECT_EQ(g_n, (std::vector<Real>{0.0, 0.0, 0.0, 0.0}));

    const Real nu = 1.8;
    const auto g_g = log_likelihood_gradient(eps, h, GarchDist::GED, nu);
    const Real delta = 1e-6;
    const Real num = (log_likelihood(eps, h, GarchDist::GED, nu + delta)
                      - log_likelihood(eps, h, GarchDist::GED, nu - delta))
                     / (2.0 * delta);
    EXPECT_NEAR(g_g[3], num, 1e-5);
}
