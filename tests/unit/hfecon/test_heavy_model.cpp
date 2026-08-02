// =============================================================================
// test_heavy_model.cpp
// Phase 5 v1.4.2 Wave C - HEAVY 模型测试
//
// 对标: R highfrequency 1.0.3 HEAVYmodel
// 容差: 1e-8 (数值优化, 宽松容差)
//
// SOURCE: PHASE5_HFE_SPEC §5.2, §5.3 D8, §5.5
//   R highfrequency 1.0.3 src/HEAVYmodel.cpp calcRecVarEq (L5-15)
//   R highfrequency 1.0.3 R/HEAVYmodel.R HEAVYmodel (L72-137)
//   R highfrequency 1.0.3 R/internalHEAVY.R heavyLLH (L32-39)
//
// 关键幻觉排除 (spec §5.3 D8):
//   强制去均值 (data - mean(data))
//   g[0] = mean(rm), 不是 mean(ret^2)
//   方差方程 condVar = calcRecVarEq(par, rm)
//   RM 方程 condVar = calcRecVarEq(par, rm) (相同递归, 不同似然)
// =============================================================================
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/models/heavy_model.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace cpphub::v1::hfecon;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {
constexpr Real TOL_STRICT = 1e-12;
constexpr Real TOL_LOOSE  = 1e-6;
constexpr Real TOL_VERY_LOOSE = 1e-3;
}  // namespace

// =============================================================================
// 辅助函数: 生成 ret/rm 序列 (固定种子 LCG, 跨平台一致)
// =============================================================================
namespace {
void make_ret_rm_series(Size n, Real base_ret, Real base_rm,
                        Real vol_ret, Real vol_rm,
                        uint64_t seed,
                        std::vector<Real>& ret,
                        std::vector<Real>& rm) {
    ret.resize(n);
    rm.resize(n);
    for (Size i = 0; i < n; ++i) {
        // PCG-like LCG (Numerical Recipes 3rd ed.)
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u1 = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g1 = (u1 - 0.5) * 2.0;  // [-1, 1]
        seed = (6364136223846793005ULL * seed + 1442695040888963407ULL);
        const Real u2 = static_cast<Real>((seed >> 11) & 0xFFFFFF) / 16777216.0;
        const Real g2 = (u2 - 0.5) * 2.0;

        // ret ~ N(0, base_ret^2) 近似 (Box-Muller 简化)
        // 用 g1 * base_ret 作为收益率
        ret[i] = g1 * base_ret;
        // rm > 0, 接近 base_rm
        rm[i] = std::max(base_rm * (1.0 + vol_rm * g2), 1e-8);
    }
}
}  // namespace

// =============================================================================
// TEST 1: calc_rec_var_eq 基础递归 — 已知参数手算
//   par = [0.1, 0.4, 0.5], rm = [1, 2, 3, 4]
//   g[0] = mean(rm) = 2.5
//   g[1] = 0.1 + 0.4*rm[0] + 0.5*g[0] = 0.1 + 0.4*1 + 0.5*2.5 = 1.75
//   g[2] = 0.1 + 0.4*rm[1] + 0.5*g[1] = 0.1 + 0.4*2 + 0.5*1.75 = 1.775
//   g[3] = 0.1 + 0.4*rm[2] + 0.5*g[2] = 0.1 + 0.4*3 + 0.5*1.775 = 2.1875
// =============================================================================
TEST(HeavyRecVarEqTest, BasicRecursion) {
    std::vector<Real> par = {0.1, 0.4, 0.5};
    std::vector<Real> rm = {1.0, 2.0, 3.0, 4.0};

    auto g = calc_rec_var_eq(par, rm);

    ASSERT_EQ(g.size(), 4u);
    EXPECT_NEAR(g[0], 2.5, TOL_STRICT);     // mean(rm) = (1+2+3+4)/4
    EXPECT_NEAR(g[1], 1.75, TOL_STRICT);    // 0.1 + 0.4*1 + 0.5*2.5
    EXPECT_NEAR(g[2], 1.775, TOL_STRICT);   // 0.1 + 0.4*2 + 0.5*1.75
    EXPECT_NEAR(g[3], 2.1875, TOL_STRICT);  // 0.1 + 0.4*3 + 0.5*1.775
}

// =============================================================================
// TEST 2: calc_rec_var_eq g[0] = mean(rm) — 关键幻觉排除
//   R 源码 (HEAVYmodel.cpp L7): g[0] = mean(rm), 不是 mean(ret^2)
// =============================================================================
TEST(HeavyRecVarEqTest, G0IsMeanRM) {
    std::vector<Real> par = {0.2, 0.3, 0.4};
    std::vector<Real> rm = {2.0, 4.0, 6.0, 8.0};

    auto g = calc_rec_var_eq(par, rm);

    // mean(rm) = 5.0, 不是 mean(ret^2) 或其他
    EXPECT_NEAR(g[0], 5.0, TOL_STRICT);
}

// =============================================================================
// TEST 3: calc_rec_var_eq 稳定性 — g[i] <= 0 时截断为 1e-10
//   构造极端参数使 g[i] < 0, 验证截断
// =============================================================================
TEST(HeavyRecVarEqTest, NegativeVarianceClipping) {
    // omega = -10 (违反约束但函数不强制检查), 使 g[i] 可能为负
    std::vector<Real> par = {-10.0, 0.5, 0.5};
    std::vector<Real> rm = {0.1, 0.1, 0.1};

    auto g = calc_rec_var_eq(par, rm);

    ASSERT_EQ(g.size(), 3u);
    EXPECT_NEAR(g[0], 0.1, TOL_STRICT);  // mean(rm)
    // g[1] = -10 + 0.5*0.1 + 0.5*0.1 = -9.9 → 截断为 1e-10
    EXPECT_GT(g[1], 0.0);
    EXPECT_NEAR(g[1], 1e-10, TOL_LOOSE);
}

// =============================================================================
// TEST 4: calc_rec_var_eq 空输入 / 参数不足
// =============================================================================
TEST(HeavyRecVarEqTest, ExceptionHandling) {
    std::vector<Real> empty;
    std::vector<Real> par = {0.1, 0.2, 0.3};

    // 空输入返回空向量
    auto g_empty = calc_rec_var_eq(par, empty);
    EXPECT_TRUE(g_empty.empty());

    // 参数不足抛异常
    std::vector<Real> rm = {1.0, 2.0, 3.0};
    std::vector<Real> par_short = {0.1, 0.2};
    EXPECT_THROW(calc_rec_var_eq(par_short, rm), std::invalid_argument);
}

// =============================================================================
// TEST 5: heavy_llh 方差方程似然 — RMEq=FALSE
//   验证似然值公式: -1/2 * log(2pi) - 1/2 * (log(g) + ret^2/g)
//
//   R 源码 (internalHEAVY.R L36): condVar <- calcRecVarEq(par, ret^2)
//   注意: 方差方程 MLE 用 ret^2 递归, 不是 rm (关键幻觉排除)
// =============================================================================
TEST(HeavyLlhTest, VarEquationLikelihood) {
    std::vector<Real> par = {0.1, 0.4, 0.5};
    std::vector<Real> ret = {0.01, -0.02, 0.015, 0.005};
    std::vector<Real> rm = {1.0, 2.0, 3.0, 4.0};

    auto llh = heavy_llh(par, ret, rm, false);
    ASSERT_EQ(llh.size(), 4u);

    // 手算第一个观测的似然
    // R 源码: condVar = calcRecVarEq(par, ret^2), g[0] = mean(ret^2)
    // ret^2 = [0.0001, 0.0004, 0.000225, 0.000025]
    // mean(ret^2) = 0.00075 / 4 = 0.0001875
    // llh[0] = -0.5*log(2pi) - 0.5*(log(g[0]) + ret[0]^2/g[0])
    const Real log_2pi = std::log(2.0 * cpphub::v1::PI);
    const Real g0 = 0.0001875;  // mean(ret^2)
    const Real expected_llh0 = -0.5 * log_2pi - 0.5 * (std::log(g0) + 0.01*0.01 / g0);
    EXPECT_NEAR(llh[0], expected_llh0, TOL_STRICT);

    // 所有似然值应有限
    for (Real v : llh) {
        EXPECT_TRUE(std::isfinite(v));
    }
}

// =============================================================================
// TEST 6: heavy_llh RM 方程似然 — RMEq=TRUE
//   验证似然值公式: -1/2 * log(2pi) - 1/2 * (log(g) + rm/g)
// =============================================================================
TEST(HeavyLlhTest, RMEquationLikelihood) {
    std::vector<Real> par = {0.1, 0.4, 0.5};
    std::vector<Real> ret = {0.01, -0.02, 0.015, 0.005};
    std::vector<Real> rm = {1.0, 2.0, 3.0, 4.0};

    auto llh = heavy_llh(par, ret, rm, true);
    ASSERT_EQ(llh.size(), 4u);

    // 手算第一个观测的似然
    // g[0] = mean(rm) = 2.5
    // llh[0] = -0.5*log(2pi) - 0.5*(log(2.5) + rm[0]/2.5)
    //       = -0.5*log(2pi) - 0.5*(log(2.5) + 1.0/2.5)
    const Real log_2pi = std::log(2.0 * cpphub::v1::PI);
    const Real expected_llh0 = -0.5 * log_2pi - 0.5 * (std::log(2.5) + 1.0 / 2.5);
    EXPECT_NEAR(llh[0], expected_llh0, TOL_STRICT);
}

// =============================================================================
// TEST 7: HEAVY 模型估计 — 基本估计收敛
//   构造 200 天 ret/rm 数据, 验证估计能收敛
// =============================================================================
TEST(HeavyModelTest, BasicEstimation) {
    std::vector<Real> ret, rm;
    make_ret_rm_series(200, 0.01, 0.001, 0.5, 0.3, 42, ret, rm);

    auto result = HeavyModel::estimate(ret, rm);

    // 验证模型结构
    EXPECT_EQ(result.n_obs, 200u);
    EXPECT_TRUE(result.converged);

    // 系数: [omega, alpha, beta, omegaR, alphaR, betaR]
    ASSERT_EQ(result.coefficients.size(), 6u);
    ASSERT_EQ(result.var_coefficients.size(), 3u);
    ASSERT_EQ(result.rm_coefficients.size(), 3u);
    EXPECT_EQ(result.coef_names.size(), 6u);
    EXPECT_EQ(result.coef_names[0], "omega");
    EXPECT_EQ(result.coef_names[1], "alpha");
    EXPECT_EQ(result.coef_names[2], "beta");
    EXPECT_EQ(result.coef_names[3], "omegaR");
    EXPECT_EQ(result.coef_names[4], "alphaR");
    EXPECT_EQ(result.coef_names[5], "betaR");

    // 所有系数应有限
    for (Real c : result.coefficients) {
        EXPECT_TRUE(std::isfinite(c)) << "Coefficient is not finite";
    }

    // 约束检查: omega > 0, alpha > 0, beta > 0
    // NelderMead penalty 不保证严格约束, 允许小范围违反 (TOL_VERY_LOOSE)
    EXPECT_GE(result.var_coefficients[0], -TOL_VERY_LOOSE);  // omega
    EXPECT_GE(result.var_coefficients[1], -TOL_VERY_LOOSE);  // alpha
    EXPECT_GE(result.var_coefficients[2], -TOL_VERY_LOOSE);  // beta
    EXPECT_GE(result.rm_coefficients[0], -TOL_VERY_LOOSE);   // omegaR
    EXPECT_GE(result.rm_coefficients[1], -TOL_VERY_LOOSE);   // alphaR
    EXPECT_GE(result.rm_coefficients[2], -TOL_VERY_LOOSE);   // betaR

    // 对数似然应为有限实数
    EXPECT_TRUE(std::isfinite(result.llh_var));
    EXPECT_TRUE(std::isfinite(result.llh_rm));
}

// =============================================================================
// TEST 8: HEAVY 模型 — 条件方差序列
//   验证 var_cond_variances 和 rm_cond_variances 长度正确
// =============================================================================
TEST(HeavyModelTest, ConditionalVariances) {
    std::vector<Real> ret, rm;
    make_ret_rm_series(100, 0.01, 0.001, 0.5, 0.3, 77, ret, rm);

    auto result = HeavyModel::estimate(ret, rm);

    ASSERT_EQ(result.var_cond_variances.size(), 100u);
    ASSERT_EQ(result.rm_cond_variances.size(), 100u);
    ASSERT_EQ(result.residuals.size(), 100u);

    // 条件方差应为正
    for (Size i = 0; i < 100; ++i) {
        EXPECT_GT(result.var_cond_variances[i], 0.0)
            << "var_cond_variances[" << i << "] <= 0";
        EXPECT_GT(result.rm_cond_variances[i], 0.0)
            << "rm_cond_variances[" << i << "] <= 0";
        EXPECT_TRUE(std::isfinite(result.residuals[i]));
    }

    // 残差 = ret / sqrt(h_t) (去均值后的 ret)
    for (Size i = 0; i < 100; ++i) {
        Real expected_res = result.ret[i] / std::sqrt(result.var_cond_variances[i]);
        EXPECT_NEAR(result.residuals[i], expected_res, TOL_LOOSE);
    }
}

// =============================================================================
// TEST 9: HEAVY 模型 — 去均值验证
//   R 源码 (HEAVYmodel.R L76): ret <- ret - mean(ret)
//   传入非零均值的 ret, 验证内部去均值
// =============================================================================
TEST(HeavyModelTest, DemeanRet) {
    std::vector<Real> ret, rm;
    make_ret_rm_series(100, 0.01, 0.001, 0.5, 0.3, 42, ret, rm);

    // 故意加一个非零均值
    Real bias = 0.001;
    std::vector<Real> ret_biased(100);
    for (Size i = 0; i < 100; ++i) ret_biased[i] = ret[i] + bias;

    auto result = HeavyModel::estimate(ret_biased, rm);

    // 内部去均值后, result.ret 的均值应接近 0
    Real mean_ret = 0.0;
    for (Size i = 0; i < 100; ++i) mean_ret += result.ret[i];
    mean_ret /= 100.0;
    EXPECT_NEAR(mean_ret, 0.0, TOL_LOOSE);
}

// =============================================================================
// TEST 10: HEAVY 模型 — 自定义起始值
//   传入 starting_values, 验证使用而非自动计算
// =============================================================================
TEST(HeavyModelTest, CustomStartingValues) {
    std::vector<Real> ret, rm;
    make_ret_rm_series(150, 0.01, 0.001, 0.5, 0.3, 42, ret, rm);

    // 自定义起始值
    std::vector<Real> start = {0.0001, 0.3, 0.5, 0.0001, 0.6, 0.3};
    auto result = HeavyModel::estimate(ret, rm, start);

    EXPECT_EQ(result.n_obs, 150u);
    ASSERT_EQ(result.coefficients.size(), 6u);
    for (Real c : result.coefficients) {
        EXPECT_TRUE(std::isfinite(c));
    }
}

// =============================================================================
// TEST 11: HEAVY 模型 — 预测功能 (1 步)
//   predict_one_step 应返回 (h_{T+1}, mu_{T+1})
// =============================================================================
TEST(HeavyModelTest, PredictOneStep) {
    std::vector<Real> ret, rm;
    make_ret_rm_series(100, 0.01, 0.001, 0.5, 0.3, 42, ret, rm);

    auto result = HeavyModel::estimate(ret, rm);

    auto pred = HeavyModel::predict_one_step(result);

    // 预测值应有限
    EXPECT_TRUE(std::isfinite(pred.first));   // h_{T+1}
    EXPECT_TRUE(std::isfinite(pred.second));  // mu_{T+1}

    // 手算预测值
    const Size n = ret.size();
    const Real omega = result.var_coefficients[0];
    const Real alpha = result.var_coefficients[1];
    const Real beta = result.var_coefficients[2];
    const Real omega_r = result.rm_coefficients[0];
    const Real alpha_r = result.rm_coefficients[1];
    const Real beta_r = result.rm_coefficients[2];

    const Real last_rm = result.rm[n - 1];
    const Real last_h = result.var_cond_variances[n - 1];
    const Real last_mu = result.rm_cond_variances[n - 1];

    Real expected_h = omega + alpha * last_rm + beta * last_h;
    Real expected_mu = omega_r + alpha_r * last_rm + beta_r * last_mu;

    EXPECT_NEAR(pred.first, expected_h, TOL_LOOSE);
    EXPECT_NEAR(pred.second, expected_mu, TOL_LOOSE);
}

// =============================================================================
// TEST 12: HEAVY 模型 — 多步预测
//   predict_multi_step 应返回 steps_ahead 个预测
// =============================================================================
TEST(HeavyModelTest, PredictMultiStep) {
    std::vector<Real> ret, rm;
    make_ret_rm_series(100, 0.01, 0.001, 0.5, 0.3, 42, ret, rm);

    auto result = HeavyModel::estimate(ret, rm);

    Size steps = 5;
    auto preds = HeavyModel::predict_multi_step(result, steps);

    ASSERT_EQ(preds.size(), steps);
    for (Size k = 0; k < steps; ++k) {
        EXPECT_TRUE(std::isfinite(preds[k].first));
        EXPECT_TRUE(std::isfinite(preds[k].second));
    }

    // 1 步预测应与 predict_one_step 一致
    auto pred1 = HeavyModel::predict_one_step(result);
    EXPECT_NEAR(preds[0].first, pred1.first, TOL_LOOSE);
    EXPECT_NEAR(preds[0].second, pred1.second, TOL_LOOSE);
}

// =============================================================================
// TEST 13: HEAVY 模型 — 异常处理
//   数据不足、长度不匹配、空输入
// =============================================================================
TEST(HeavyModelTest, ExceptionHandling) {
    // 空输入
    std::vector<Real> empty;
    EXPECT_THROW(HeavyModel::estimate(empty, empty), std::invalid_argument);

    // 数据不足 (n < 5)
    std::vector<Real> ret_short = {0.01, 0.02, 0.03};
    std::vector<Real> rm_short = {0.001, 0.002, 0.003};
    EXPECT_THROW(HeavyModel::estimate(ret_short, rm_short),
                 std::invalid_argument);

    // 长度不匹配
    std::vector<Real> ret100(100, 0.01);
    std::vector<Real> rm90(90, 0.001);
    EXPECT_THROW(HeavyModel::estimate(ret100, rm90), std::invalid_argument);
}

// =============================================================================
// TEST 14: HEAVY 模型 — 起始值边界
//   starting_values 长度不足 6 应自动用默认起始值
// =============================================================================
TEST(HeavyModelTest, ShortStartingValues) {
    std::vector<Real> ret, rm;
    make_ret_rm_series(100, 0.01, 0.001, 0.5, 0.3, 42, ret, rm);

    // 长度不足 6, 应自动用默认起始值
    std::vector<Real> short_start = {0.0001, 0.3};
    EXPECT_NO_THROW({
        auto result = HeavyModel::estimate(ret, rm, short_start);
        EXPECT_EQ(result.n_obs, 100u);
        ASSERT_EQ(result.coefficients.size(), 6u);
    });
}
