// =============================================================================
// test_har_model.cpp
// Phase 5 v1.4.2 Wave C - HAR 模型测试
//
// 对标: R highfrequency 1.0.3 HARmodel
// 容差: 1e-8 (迭代算法, 宽松容差)
//
// SOURCE: PHASE5_HFE_SPEC §5.2, §5.3 D7/D9, §5.5
//   R highfrequency 1.0.3 src/HARmodel.cpp har_agg (L7-21)
//   R highfrequency 1.0.3 R/HARmodel.R HARmodel (L189-514)
//
// 关键幻觉排除 (spec §5.3 D7/D9):
//   D7: RQ 变换用 BPQ 2016 (sqrt(RQ) - sqrt(mean(RM3)))
//   D9: har_agg 索引边界 [j-p, j-1] 含 j-1 不含 j
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/models/har_model.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace cpphub::v1::hfecon;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
constexpr Real TOL_STRICT = 1e-12;
constexpr Real TOL_LOOSE  = 1e-8;
}  // namespace

// =============================================================================
// 辅助函数: 生成 RV 序列 (固定种子)
// =============================================================================
namespace {
std::vector<Real> make_rv_series(Size n, Real base, Real vol, uint64_t seed) {
    std::vector<Real> rv(n);
    for (Size i = 0; i < n; ++i) {
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g = (u - 0.5) * 2.0;  // [-1, 1]
        rv[i] = std::max(base * (1.0 + vol * g), 1e-8);
    }
    return rv;
}
}  // namespace

// =============================================================================
// TEST 1: har_agg 基础功能 — 滚动窗口平均
//   RM = [1, 2, 3, 4, 5], periods = [2]
//   har_agg[0] = NaN (p=2, row=0 < p-1=1)
//   har_agg[1] = (1+2)/2 = 1.5
//   har_agg[2] = (2+3)/2 = 2.5
//   har_agg[3] = (3+4)/2 = 3.5
//   har_agg[4] = (4+5)/2 = 4.5
// =============================================================================
TEST(HarAggTest, BasicRollingMean) {
    std::vector<Real> RM = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<Size> periods = {2};
    auto result = har_agg(RM, periods);

    ASSERT_EQ(result.size(), 5u);
    ASSERT_EQ(result[0].size(), 1u);
    EXPECT_TRUE(std::isnan(result[0][0]));  // row=0 < p-1=1
    EXPECT_NEAR(result[1][0], 1.5, TOL_STRICT);
    EXPECT_NEAR(result[2][0], 2.5, TOL_STRICT);
    EXPECT_NEAR(result[3][0], 3.5, TOL_STRICT);
    EXPECT_NEAR(result[4][0], 4.5, TOL_STRICT);
}

// =============================================================================
// TEST 2: har_agg 多周期 — periods = [1, 2, 3]
// =============================================================================
TEST(HarAggTest, MultiPeriod) {
    std::vector<Real> RM = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<Size> periods = {1, 2, 3};
    auto result = har_agg(RM, periods);

    ASSERT_EQ(result.size(), 6u);
    ASSERT_EQ(result[0].size(), 3u);

    // period=1: 直接复制
    EXPECT_NEAR(result[0][0], 1.0, TOL_STRICT);
    EXPECT_NEAR(result[5][0], 6.0, TOL_STRICT);

    // period=2: (RM[i-1]+RM[i])/2
    EXPECT_TRUE(std::isnan(result[0][1]));
    EXPECT_NEAR(result[1][1], 1.5, TOL_STRICT);
    EXPECT_NEAR(result[5][1], 5.5, TOL_STRICT);

    // period=3: (RM[i-2]+RM[i-1]+RM[i])/3
    EXPECT_TRUE(std::isnan(result[0][2]));
    EXPECT_TRUE(std::isnan(result[1][2]));
    EXPECT_NEAR(result[2][2], 2.0, TOL_STRICT);
    EXPECT_NEAR(result[5][2], 5.0, TOL_STRICT);
}

// =============================================================================
// TEST 3: har_agg_single — 单周期聚合
// =============================================================================
TEST(HarAggTest, SinglePeriod) {
    std::vector<Real> RM = {1.0, 2.0, 3.0, 4.0, 5.0};
    auto result = har_agg_single(RM, 3);

    ASSERT_EQ(result.size(), 5u);
    EXPECT_TRUE(std::isnan(result[0]));
    EXPECT_TRUE(std::isnan(result[1]));
    EXPECT_NEAR(result[2], 2.0, TOL_STRICT);
    EXPECT_NEAR(result[3], 3.0, TOL_STRICT);
    EXPECT_NEAR(result[4], 4.0, TOL_STRICT);
}

// =============================================================================
// TEST 4: har_insanity_filter — BPQ insanity filter
// =============================================================================
TEST(HarAggTest, InsanityFilter) {
    std::vector<Real> fitted = {1.0, -5.0, 3.0, 100.0, 2.0};
    Real lower = 0.0, upper = 10.0, replacement = 2.0;
    auto result = har_insanity_filter(fitted, lower, upper, replacement);

    EXPECT_NEAR(result[0], 1.0, TOL_STRICT);     // 在范围内
    EXPECT_NEAR(result[1], 2.0, TOL_STRICT);     // < lower → replacement
    EXPECT_NEAR(result[2], 3.0, TOL_STRICT);     // 在范围内
    EXPECT_NEAR(result[3], 2.0, TOL_STRICT);     // > upper → replacement
    EXPECT_NEAR(result[4], 2.0, TOL_STRICT);     // 在范围内
}

// =============================================================================
// TEST 5: 基础 HAR-RV 模型 — type="HAR"
//   构造 100 天 RV 序列, periods=[1,5,22], h=1
//   验证系数估计合理 (非 NaN, 有限)
// =============================================================================
TEST(HarModelTest, BasicHAR) {
    auto RM1 = make_rv_series(100, 0.01, 0.3, 42);
    auto result = HarModel::estimate_har(RM1, {1, 5, 22}, 1, HarTransform::None);

    // 验证模型结构
    EXPECT_EQ(result.type, HarType::HAR);
    EXPECT_EQ(result.transform, HarTransform::None);
    EXPECT_EQ(result.h, 1u);
    EXPECT_EQ(result.maxp, 22u);

    // 系数: [beta0, RV1, RV5, RV22]
    ASSERT_EQ(result.coefficients.size(), 4u);
    EXPECT_EQ(result.coef_names.size(), 4u);
    EXPECT_EQ(result.coef_names[0], "beta0");
    EXPECT_EQ(result.coef_names[1], "RV1");
    EXPECT_EQ(result.coef_names[2], "RV5");
    EXPECT_EQ(result.coef_names[3], "RV22");

    // 所有系数应有限
    for (Real c : result.coefficients) {
        EXPECT_TRUE(std::isfinite(c));
    }

    // R^2 应在 [0, 1]
    EXPECT_GE(result.r_squared, 0.0);
    EXPECT_LE(result.r_squared, 1.0);

    // 拟合值和残差应有限
    EXPECT_EQ(result.fitted_values.size(), result.n_obs);
    EXPECT_EQ(result.residuals.size(), result.n_obs);
    for (Size i = 0; i < result.n_obs; ++i) {
        EXPECT_TRUE(std::isfinite(result.fitted_values[i]));
        EXPECT_TRUE(std::isfinite(result.residuals[i]));
    }
}

// =============================================================================
// TEST 6: HAR-RV-J 模型 — type="HARJ"
//   需要提供 RM1 (RV) 和 RM2 (BPV)
//   注意: RM2 必须有独立随机性, 否则 J = RM1 - RM2 与 RM1 完美共线性
//         导致 OLS 设计矩阵奇异 (实测: RM2 = RM1*0.95 会触发 singular matrix)
// =============================================================================
TEST(HarModelTest, HARJ) {
    auto RM1 = make_rv_series(100, 0.01, 0.3, 42);
    // RM2 (BPV) 用独立随机序列, 略小于 RM1 (因跳跃)
    auto RM2 = make_rv_series(100, 0.0095, 0.3, 99);  // 独立种子

    auto result = HarModel::estimate_harj(RM1, RM2, {1, 5, 22}, {1, 5, 22});

    EXPECT_EQ(result.type, HarType::HARJ);
    // 系数: [beta0, RV1, RV5, RV22, J1, J5, J22]
    ASSERT_EQ(result.coefficients.size(), 7u);
    EXPECT_EQ(result.coef_names[0], "beta0");
    EXPECT_EQ(result.coef_names[1], "RV1");
    EXPECT_EQ(result.coef_names[4], "J1");

    for (Real c : result.coefficients) {
        EXPECT_TRUE(std::isfinite(c));
    }
}

// =============================================================================
// TEST 7: HAR-Q 模型 — type="HARQ" (含 BPQ 四阶矩变换)
// =============================================================================
TEST(HarModelTest, HARQ) {
    auto RM1 = make_rv_series(100, 0.01, 0.3, 42);
    auto RM3 = make_rv_series(100, 0.0001, 0.5, 99);  // RQ

    auto result = HarModel::estimate_harq(RM1, RM3, {1, 5, 22}, {1});

    EXPECT_EQ(result.type, HarType::HARQ);
    // 系数: [beta0, RV1, RV5, RV22, RQ1]
    ASSERT_EQ(result.coefficients.size(), 5u);
    EXPECT_EQ(result.coef_names[0], "beta0");
    EXPECT_EQ(result.coef_names[1], "RV1");
    EXPECT_EQ(result.coef_names[4], "RQ1");

    for (Real c : result.coefficients) {
        EXPECT_TRUE(std::isfinite(c));
    }
}

// =============================================================================
// TEST 8: CHAR 模型 — type="CHAR" (连续成分)
// =============================================================================
TEST(HarModelTest, CHAR) {
    auto RM1 = make_rv_series(100, 0.01, 0.3, 42);
    std::vector<Real> RM2(100);
    for (Size i = 0; i < 100; ++i) {
        RM2[i] = RM1[i] * 0.95;
    }

    auto result = HarModel::estimate_char(RM1, RM2, {1, 5, 22});

    EXPECT_EQ(result.type, HarType::CHAR);
    // 系数: [beta0, RV1, RV5, RV22] (用 BPV 替代 RV)
    ASSERT_EQ(result.coefficients.size(), 4u);

    for (Real c : result.coefficients) {
        EXPECT_TRUE(std::isfinite(c));
    }
}

// =============================================================================
// TEST 9: transform="log" — 对数变换
// =============================================================================
TEST(HarModelTest, LogTransform) {
    auto RM1 = make_rv_series(100, 0.01, 0.3, 42);
    auto result = HarModel::estimate_har(RM1, {1, 5, 22}, 1, HarTransform::Log);

    EXPECT_EQ(result.transform, HarTransform::Log);

    // 变换后 y 和 x 应为对数
    ASSERT_FALSE(result.y.empty());
    for (Size i = 0; i < result.y.size(); ++i) {
        EXPECT_TRUE(std::isfinite(result.y[i]));
    }

    // 预测应通过 exp backtransform
    Real pred = HarModel::predict_last(result);
    EXPECT_GT(pred, 0.0);  // exp 总是正
    EXPECT_TRUE(std::isfinite(pred));
}

// =============================================================================
// TEST 10: transform="sqrt" — 平方根变换
// =============================================================================
TEST(HarModelTest, SqrtTransform) {
    auto RM1 = make_rv_series(100, 0.01, 0.3, 42);
    auto result = HarModel::estimate_har(RM1, {1, 5, 22}, 1, HarTransform::Sqrt);

    EXPECT_EQ(result.transform, HarTransform::Sqrt);

    // 预测应通过 x^2 backtransform
    Real pred = HarModel::predict_last(result);
    EXPECT_TRUE(std::isfinite(pred));
}

// =============================================================================
// TEST 11: h 参数 — 多步聚合
// =============================================================================
TEST(HarModelTest, MultiStepH) {
    auto RM1 = make_rv_series(100, 0.01, 0.3, 42);

    auto r1 = HarModel::estimate_har(RM1, {1, 5, 22}, 1);
    auto r5 = HarModel::estimate_har(RM1, {1, 5, 22}, 5);

    EXPECT_EQ(r1.h, 1u);
    EXPECT_EQ(r5.h, 5u);

    // h=5 时 maxp 应 >= max(periods, h) = max(22, 5) = 22
    EXPECT_GE(r5.maxp, 22u);

    // h=5 的估计样本应小于 h=1
    EXPECT_LE(r5.n_obs, r1.n_obs);
}

// =============================================================================
// TEST 12: 预测 — predict_last
// =============================================================================
TEST(HarModelTest, PredictLast) {
    auto RM1 = make_rv_series(100, 0.01, 0.3, 42);
    auto result = HarModel::estimate_har(RM1, {1, 5, 22}, 1);

    Real pred = HarModel::predict_last(result);
    EXPECT_TRUE(std::isfinite(pred));

    // 验证预测值 = beta0 + sum(beta_i * x_last_i)
    Real expected = result.coefficients[0];
    ASSERT_FALSE(result.x.empty());
    for (Size j = 0; j < result.x.back().size(); ++j) {
        expected += result.coefficients[j + 1] * result.x.back()[j];
    }
    EXPECT_NEAR(pred, expected, TOL_STRICT);
}

// =============================================================================
// TEST 13: 异常处理 — 数据不足 / 空输入
// =============================================================================
TEST(HarModelTest, ExceptionHandling) {
    // 空输入
    std::vector<Real> empty;
    EXPECT_THROW(HarModel::estimate_har(empty), std::invalid_argument);

    // 数据不足 (n < maxp + h + 1)
    std::vector<Real> short_data(20, 0.01);
    EXPECT_THROW(HarModel::estimate_har(short_data, {1, 5, 22}, 1),
                 std::invalid_argument);

    // HARJ 缺少 RM2
    std::vector<Real> RM1(100, 0.01);
    EXPECT_THROW(HarModel::estimate(RM1, {}, {}, {1,5,22}, {1,5,22}, {},
                                     HarType::HARJ, 1, HarTransform::None, 0.05),
                 std::invalid_argument);

    // HARQ 缺少 RM3
    EXPECT_THROW(HarModel::estimate(RM1, {}, {}, {1,5,22}, {}, {1},
                                     HarType::HARQ, 1, HarTransform::None, 0.05),
                 std::invalid_argument);

    // CHAR 缺少 RM2
    EXPECT_THROW(HarModel::estimate(RM1, {}, {}, {1,5,22}, {}, {},
                                     HarType::CHAR, 1, HarTransform::None, 0.05),
                 std::invalid_argument);
}

// =============================================================================
// TEST 14: OLS 估计验证 — 已知线性关系
//   构造 y = 1 + 2*x, 验证 OLS 能精确恢复系数
// =============================================================================
TEST(HarModelTest, OlsRecovery) {
    // 构造 50 个观测, y = 1 + 2*x
    std::vector<Real> y(50);
    std::vector<std::vector<Real>> X(50, std::vector<Real>(1));
    for (Size i = 0; i < 50; ++i) {
        X[i][0] = static_cast<Real>(i + 1);
        y[i] = 1.0 + 2.0 * X[i][0];
    }

    std::vector<Real> fitted, residuals;
    Real r_sq, adj_r_sq, llh;
    auto beta = ols_estimate(y, X, fitted, residuals, r_sq, adj_r_sq, llh);

    ASSERT_EQ(beta.size(), 2u);
    EXPECT_NEAR(beta[0], 1.0, TOL_STRICT);
    EXPECT_NEAR(beta[1], 2.0, TOL_STRICT);
    EXPECT_NEAR(r_sq, 1.0, TOL_LOOSE);  // 完美线性拟合
}
