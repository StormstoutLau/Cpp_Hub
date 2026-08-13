// Phase 7A Wave 0: detail/ 公共基础设施测试
//
// 测试目标:
//   1. TestResultBase (test_result_base.hpp): 聚合初始化 / 字段访问 / R htest 对应
//   2. ols_simple (ols_simple.hpp): 与 har_model.hpp::ols_estimate 对齐 + 有意差异验证
//
// 教材锚点: ADR-015 §2.0, spec §2.0.1/§2.0.2
// 排幻觉点: ols_simple 不自动加常数列 (调用方决定), 不算 adj_r_squared/llh
//
// 容差: 1e-12 (与 har_model.hpp ols_estimate 对齐, Gauss-Jordan 双精度)

#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>

#include "cpphub/econometrics/inference/detail/test_result_base.hpp"
#include "cpphub/econometrics/inference/detail/ols_simple.hpp"
#include "cpphub/hfecon/models/har_model.hpp"

using cpphub::v1::econometrics::detail::TestResultBase;
using cpphub::v1::econometrics::detail::ols_simple;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
constexpr Real TOL = 1e-12;
}  // namespace

// =============================================================================
// TestResultBase 测试
// =============================================================================

TEST(TestResultBase, AggregateInitialization) {
    TestResultBase r{2.5, 0.03, "Jarque-Bera", true};
    EXPECT_DOUBLE_EQ(r.statistic, 2.5);
    EXPECT_DOUBLE_EQ(r.p_value, 0.03);
    EXPECT_EQ(r.method_name, "Jarque-Bera");
    EXPECT_TRUE(r.reject_null);
}

TEST(TestResultBase, DefaultConstruction) {
    TestResultBase r{};
    EXPECT_DOUBLE_EQ(r.statistic, 0.0);
    EXPECT_DOUBLE_EQ(r.p_value, 0.0);
    EXPECT_TRUE(r.method_name.empty());
    EXPECT_FALSE(r.reject_null);
}

TEST(TestResultBase, NanPValueForNonStandardDistribution) {
    // Cragg-Donald 等非标准分布无 p_value, 用 NaN 表示
    TestResultBase r{15.3, std::numeric_limits<Real>::quiet_NaN(),
                     "Cragg-Donald", false};
    EXPECT_TRUE(std::isnan(r.p_value));
    EXPECT_FALSE(r.reject_null);
}

TEST(TestResultBase, CompositionPattern) {
    // 验证组合模式: JarqueBeraResult 内嵌 TestResultBase
    struct JarqueBeraResult {
        TestResultBase base;
        Real skewness;
        Real kurtosis;
    };
    JarqueBeraResult jb{{5.2, 0.074, "Jarque-Bera", false}, 0.3, 3.5};
    EXPECT_DOUBLE_EQ(jb.base.statistic, 5.2);
    EXPECT_DOUBLE_EQ(jb.skewness, 0.3);
    EXPECT_DOUBLE_EQ(jb.kurtosis, 3.5);
}

// =============================================================================
// ols_simple 测试
// =============================================================================

TEST(OlsSimple, PerfectFitLinear) {
    // y = -1 + 3x, 完美拟合
    // 调用方需自行添加常数列 (与 ols_estimate 的有意差异)
    std::vector<Real> y = {2.0, 5.0, 8.0, 11.0};  // -1+3*1=2, -1+3*2=5, ...
    std::vector<std::vector<Real>> X = {
        {1.0, 1.0},  // [常数, x]
        {1.0, 2.0},
        {1.0, 3.0},
        {1.0, 4.0}
    };
    std::vector<Real> fitted, resid;
    Real r2;
    std::vector<Real> beta = ols_simple(y, X, fitted, resid, r2);

    ASSERT_EQ(beta.size(), 2u);
    EXPECT_NEAR(beta[0], -1.0, TOL);  // 截距
    EXPECT_NEAR(beta[1], 3.0, TOL);   // 斜率
    EXPECT_NEAR(r2, 1.0, TOL);

    for (Size i = 0; i < y.size(); ++i) {
        EXPECT_NEAR(fitted[i], y[i], TOL);
        EXPECT_NEAR(resid[i], 0.0, TOL);
    }
}

TEST(OlsSimple, NoIntercept) {
    // y = 3x (无常数项), 验证 ols_simple 不自动添加常数列
    std::vector<Real> y = {3.0, 6.0, 9.0, 12.0};
    std::vector<std::vector<Real>> X = {
        {1.0},
        {2.0},
        {3.0},
        {4.0}
    };
    std::vector<Real> fitted, resid;
    Real r2;
    std::vector<Real> beta = ols_simple(y, X, fitted, resid, r2);

    ASSERT_EQ(beta.size(), 1u);
    EXPECT_NEAR(beta[0], 3.0, TOL);
    EXPECT_NEAR(r2, 1.0, TOL);
}

TEST(OlsSimple, AlignWithOlsEstimate) {
    // 与 har_model.hpp::ols_estimate 对齐验证
    // ols_estimate 自动加常数列, ols_simple 不加
    // 若给 ols_simple 传入 [1, X], 两者结果应一致 (beta/fitted/resid/r²)
    std::vector<Real> y = {1.0, 3.0, 2.0, 5.0, 4.0};
    std::vector<std::vector<Real>> X_raw = {
        {0.5},
        {1.5},
        {1.0},
        {2.5},
        {2.0}
    };

    // ols_estimate (自动加常数列)
    std::vector<Real> fitted_ref, resid_ref;
    Real r2_ref, adj_r2_ref, llh_ref;
    std::vector<Real> beta_ref = cpphub::v1::hfecon::ols_estimate(
        y, X_raw, fitted_ref, resid_ref, r2_ref, adj_r2_ref, llh_ref);

    // ols_simple (手动加常数列)
    std::vector<std::vector<Real>> X_with_const(X_raw.size(),
                                                  std::vector<Real>(2));
    for (Size i = 0; i < X_raw.size(); ++i) {
        X_with_const[i][0] = 1.0;
        X_with_const[i][1] = X_raw[i][0];
    }
    std::vector<Real> fitted, resid;
    Real r2;
    std::vector<Real> beta = ols_simple(y, X_with_const, fitted, resid, r2);

    // beta 应一致 (ols_estimate 返回 [截距, 斜率], ols_simple 返回 [常数列系数, X列系数])
    ASSERT_EQ(beta.size(), beta_ref.size());
    for (Size i = 0; i < beta.size(); ++i) {
        EXPECT_NEAR(beta[i], beta_ref[i], TOL);
    }
    EXPECT_NEAR(r2, r2_ref, TOL);
    for (Size i = 0; i < y.size(); ++i) {
        EXPECT_NEAR(fitted[i], fitted_ref[i], TOL);
        EXPECT_NEAR(resid[i], resid_ref[i], TOL);
    }
}

TEST(OlsSimple, NoisyDataRSquaredLessThanOne) {
    // 非完美拟合, R² < 1
    // y = [1, 2, 3, 4] + noise, X = [1, 2, 3, 4]
    // 用 y = [1.1, 1.9, 3.2, 3.8] (有噪声)
    std::vector<Real> y = {1.1, 1.9, 3.2, 3.8};
    std::vector<std::vector<Real>> X = {
        {1.0, 1.0},
        {1.0, 2.0},
        {1.0, 3.0},
        {1.0, 4.0}
    };
    std::vector<Real> fitted, resid;
    Real r2;
    std::vector<Real> beta = ols_simple(y, X, fitted, resid, r2);

    EXPECT_GT(r2, 0.0);
    EXPECT_LT(r2, 1.0);
    // 残差和应接近 0 (OLS 正规方程保证)
    Real sum_resid = 0.0;
    for (Real r : resid) sum_resid += r;
    EXPECT_NEAR(sum_resid, 0.0, 1e-10);
}

TEST(OlsSimple, MultivariateRegression) {
    // 多元回归: y = 1 + 2*x1 + 3*x2
    std::vector<Real> y = {1.0 + 2.0*1 + 3.0*2,
                           1.0 + 2.0*2 + 3.0*3,
                           1.0 + 2.0*3 + 3.0*1,
                           1.0 + 2.0*4 + 3.0*0};
    // y = [9, 14, 10, 9]
    std::vector<std::vector<Real>> X = {
        {1.0, 1.0, 2.0},
        {1.0, 2.0, 3.0},
        {1.0, 3.0, 1.0},
        {1.0, 4.0, 0.0}
    };
    std::vector<Real> fitted, resid;
    Real r2;
    std::vector<Real> beta = ols_simple(y, X, fitted, resid, r2);

    ASSERT_EQ(beta.size(), 3u);
    EXPECT_NEAR(beta[0], 1.0, TOL);
    EXPECT_NEAR(beta[1], 2.0, TOL);
    EXPECT_NEAR(beta[2], 3.0, TOL);
    EXPECT_NEAR(r2, 1.0, TOL);
}

TEST(OlsSimple, EmptyInputThrows) {
    std::vector<Real> y;
    std::vector<std::vector<Real>> X;
    std::vector<Real> fitted, resid;
    Real r2;
    EXPECT_THROW(ols_simple(y, X, fitted, resid, r2), std::invalid_argument);
}

TEST(OlsSimple, InsufficientObservationsThrows) {
    // N < K (3 个观测, 4 个参数)
    std::vector<Real> y = {1.0, 2.0, 3.0};
    std::vector<std::vector<Real>> X = {
        {1.0, 0.0, 0.0, 0.0},
        {1.0, 1.0, 0.0, 0.0},
        {1.0, 0.0, 1.0, 0.0}
    };
    std::vector<Real> fitted, resid;
    Real r2;
    EXPECT_THROW(ols_simple(y, X, fitted, resid, r2), std::invalid_argument);
}

TEST(OlsSimple, SingularMatrixThrows) {
    // 完全共线性: 两列相同
    std::vector<Real> y = {1.0, 2.0, 3.0, 4.0};
    std::vector<std::vector<Real>> X = {
        {1.0, 1.0},
        {1.0, 1.0},
        {1.0, 1.0},
        {1.0, 1.0}
    };
    std::vector<Real> fitted, resid;
    Real r2;
    EXPECT_THROW(ols_simple(y, X, fitted, resid, r2), std::runtime_error);
}

TEST(OlsSimple, NoAdjRSquaredOrLlh) {
    // 验证 ols_simple 不计算 adj_r_squared / llh (与 ols_estimate 的有意差异)
    // 接口签名不含这两个输出参数
    std::vector<Real> y = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<std::vector<Real>> X = {
        {1.0, 1.0},
        {1.0, 2.0},
        {1.0, 3.0},
        {1.0, 4.0},
        {1.0, 5.0}
    };
    std::vector<Real> fitted, resid;
    Real r2;
    std::vector<Real> beta = ols_simple(y, X, fitted, resid, r2);

    // 接口仅返回 beta + fitted + resid + r_squared
    // 若需要 adj_r_squared / llh, 调用方应使用 estimation/ols.hpp (Eigen3)
    EXPECT_NEAR(beta[0], 0.0, TOL);
    EXPECT_NEAR(beta[1], 1.0, TOL);
    EXPECT_NEAR(r2, 1.0, TOL);
}
