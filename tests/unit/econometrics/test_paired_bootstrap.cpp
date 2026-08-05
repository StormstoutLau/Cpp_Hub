// SOURCE: PHASE6_IMPLEMENTATION_PLAN §6.1 任务 4.2 - PairedBootstrap 测试
// 验证:
//   - 配对 Bootstrap 基本功能 (Efron-Tibshirani 1993 法学院 15 obs)
//   - 固定种子可复现性 (Philox 分块, ADR-004)
//   - 大样本下 Bootstrap 均值/SE 接近原估计 (大数定律)
//   - CI 覆盖率 (蒙特卡洛, 95% CI 应覆盖真值 ~93-97%)
//   - 失败 replicate 计数 (H-009)
//   - 排幻觉点 H-006 (配对采样), H-008 (B-1 分母), H-012 (返回 samples)

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/econometrics/estimation/ols.hpp"
#include "cpphub/econometrics/resampling/paired_bootstrap.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;

// =============================================================================
// 辅助: 构造横截面数据
// =============================================================================
namespace {

// 构造简单线性回归数据 (y = β0 + β1·x + ε)
EconData make_linear_data(Size N, Real beta0, Real beta1, Real sigma,
                          std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<Real> noise(0.0, sigma);
    std::uniform_real_distribution<Real> xdist(0.0, 10.0);

    MatrixXD X(N, 2);
    VectorXD y(N);
    for (Size i = 0; i < N; ++i) {
        const Real xi = xdist(rng);
        const Real yi = beta0 + beta1 * xi + noise(rng);
        X(i, 0) = 1.0;
        X(i, 1) = xi;
        y(i) = yi;
    }
    return make_cross_section(X, y, {"intercept", "x"}, "y");
}

// Efron-Tibshirani (1993) 法学院数据 (15 obs, 法学院 LSAT vs GPA)
// SOURCE: Efron-Tibshirani (1993) Table 3.1, 法学院 15 学校
// LSAT:  576, 635, 558, 578, 666, 580, 555, 661, 651, 605, 653, 575, 545, 572, 594
// GPA:   3.39, 3.30, 2.81, 3.03, 3.44, 3.07, 3.00, 3.43, 3.36, 3.13, 3.12, 2.74, 2.76, 2.88, 2.96
EconData make_law_school_data() {
    static const std::vector<Real> lsat = {
        576, 635, 558, 578, 666, 580, 555, 661, 651, 605,
        653, 575, 545, 572, 594
    };
    static const std::vector<Real> gpa = {
        3.39, 3.30, 2.81, 3.03, 3.44, 3.07, 3.00, 3.43, 3.36, 3.13,
        3.12, 2.74, 2.76, 2.88, 2.96
    };
    const Size N = lsat.size();
    MatrixXD X(N, 2);
    VectorXD y(N);
    for (Size i = 0; i < N; ++i) {
        X(i, 0) = 1.0;
        X(i, 1) = lsat[i];
        y(i) = gpa[i];
    }
    return make_cross_section(X, y, {"intercept", "LSAT"}, "GPA");
}

}  // namespace

// =============================================================================
// 测试 1: 基本 resample - 不抛异常, 返回正确结构
// =============================================================================
TEST(PairedBootstrap, BasicResample_ReturnsValidResult) {
    auto data = make_linear_data(50, 1.0, 2.0, 0.5, 42);
    OLSEstimator ols;
    PairedBootstrap engine;

    BootstrapResult r = engine.resample(ols, data, 199, 42);

    EXPECT_EQ(r.n_replicates, 199u);
    EXPECT_EQ(r.n_failed, 0u);
    EXPECT_EQ(r.bootstrap_samples.size(), 199u);
    EXPECT_EQ(r.coef_mean.size(), 2u);
    EXPECT_EQ(r.coef_std.size(), 2u);
    EXPECT_EQ(r.coef_vcov.rows(), 2u);
    EXPECT_EQ(r.coef_vcov.cols(), 2u);
}

// =============================================================================
// 测试 2: 维度匹配 - 系数维度与估计器一致
// =============================================================================
TEST(PairedBootstrap, ResultDimensions_MatchEstimator) {
    MatrixXD X(30, 3);  // K=3 (截距 + 2 个解释变量)
    VectorXD y(30);
    std::mt19937_64 rng(123);
    std::normal_distribution<Real> n(0.0, 1.0);
    for (Size i = 0; i < 30; ++i) {
        X(i, 0) = 1.0;
        X(i, 1) = n(rng);
        X(i, 2) = n(rng);
        y(i) = 1.0 + 0.5 * X(i, 1) - 0.3 * X(i, 2) + n(rng) * 0.5;
    }
    EconData data = make_cross_section(X, y, {"b0", "x1", "x2"}, "y");

    OLSEstimator ols;
    PairedBootstrap engine;
    BootstrapResult r = engine.resample(ols, data, 99, 7);

    EXPECT_EQ(r.coef_mean.size(), 3u);
    EXPECT_EQ(r.coef_std.size(), 3u);
    EXPECT_EQ(r.coef_vcov.rows(), 3u);
    EXPECT_EQ(r.coef_vcov.cols(), 3u);
    for (const auto& s : r.bootstrap_samples) {
        EXPECT_EQ(s.size(), 3u);
    }
}

// =============================================================================
// 测试 3: 固定种子可复现 - 同 (seed, B) 产出位精确相同结果
//   排幻觉点 H-013: Philox4x64(seed, b) 分块, 确定性可复现
// =============================================================================
TEST(PairedBootstrap, Reproducible_SameSeed) {
    auto data = make_linear_data(40, 0.5, 1.5, 0.3, 100);
    OLSEstimator ols;
    PairedBootstrap engine;

    BootstrapResult r1 = engine.resample(ols, data, 199, 42);
    BootstrapResult r2 = engine.resample(ols, data, 199, 42);

    EXPECT_EQ(r1.bootstrap_samples.size(), r2.bootstrap_samples.size());
    for (Size b = 0; b < r1.bootstrap_samples.size(); ++b) {
        for (Size j = 0; j < r1.bootstrap_samples[b].size(); ++j) {
            EXPECT_NEAR(r1.bootstrap_samples[b](j),
                        r2.bootstrap_samples[b](j), 1e-15);
        }
    }
    EXPECT_NEAR(r1.lower_ci, r2.lower_ci, 1e-15);
    EXPECT_NEAR(r1.upper_ci, r2.upper_ci, 1e-15);
}

// =============================================================================
// 测试 4: 不同种子产生不同结果 (统计独立)
// =============================================================================
TEST(PairedBootstrap, DifferentSeeds_DifferentResults) {
    auto data = make_linear_data(100, 1.0, 2.0, 0.5, 200);
    OLSEstimator ols;
    PairedBootstrap engine;

    BootstrapResult r1 = engine.resample(ols, data, 499, 42);
    BootstrapResult r2 = engine.resample(ols, data, 499, 999);

    // 均值相近 (大样本), 但具体值不同
    EXPECT_NEAR(r1.coef_mean(1), r2.coef_mean(1), 0.2);  // 同分布
    // 95% CI 不同
    EXPECT_NE(r1.lower_ci, r2.lower_ci);
}

// =============================================================================
// 测试 5: 全部估计成功 - n_failed=0 (OLS 不会失败)
// =============================================================================
TEST(PairedBootstrap, NoFailedReplicates_ValidEstimator) {
    auto data = make_linear_data(50, 1.0, 2.0, 0.5, 300);
    OLSEstimator ols;
    PairedBootstrap engine;

    BootstrapResult r = engine.resample(ols, data, 299, 42);
    EXPECT_EQ(r.n_failed, 0u);
    EXPECT_EQ(r.bootstrap_samples.size(), 299u);
}

// =============================================================================
// 测试 6: 大样本 Bootstrap 均值接近原估计
//   大数定律: θ̄* → θ̂ 当 B → ∞
// =============================================================================
TEST(PairedBootstrap, BootstrapMean_CloseToOriginal) {
    auto data = make_linear_data(200, 1.0, 2.0, 0.5, 400);
    OLSEstimator ols;
    EstimationResult orig = ols.estimate(data);

    PairedBootstrap engine;
    BootstrapResult r = engine.resample(ols, data, 999, 42);

    // 大样本 B=999, Bootstrap 均值应接近原估计 (容差 ~2·SE)
    EXPECT_NEAR(r.coef_mean(0), orig.coefficients(0), 0.1);
    EXPECT_NEAR(r.coef_mean(1), orig.coefficients(1), 0.05);
}

// =============================================================================
// 测试 7: 大样本 Bootstrap SE 接近 OLS Classical SE
//   同方差下 Bootstrap SE ≈ OLS SE (n^{-1/2} 比例)
// =============================================================================
TEST(PairedBootstrap, BootstrapSE_CloseToOLSSe) {
    auto data = make_linear_data(300, 1.0, 2.0, 0.5, 500);
    OLSEstimator ols(CovarianceType::Classical);
    EstimationResult orig = ols.estimate(data);

    PairedBootstrap engine;
    BootstrapResult r = engine.resample(ols, data, 999, 42);

    // Bootstrap SE 应在 OLS SE 的 ±20% 范围内 (Bootstrap 自身有 ~1/√B 波动)
    const Real se0_ratio = r.coef_std(0) / orig.std_errors(0);
    const Real se1_ratio = r.coef_std(1) / orig.std_errors(1);
    EXPECT_GT(se0_ratio, 0.80);
    EXPECT_LT(se0_ratio, 1.20);
    EXPECT_GT(se1_ratio, 0.80);
    EXPECT_LT(se1_ratio, 1.20);
}

// =============================================================================
// 测试 8: 法学院数据 (Efron-Tibshirani 1993, 15 obs)
//   - 原样本估计稳定 (β1 > 0, GPA 与 LSAT 正相关)
//   - Bootstrap 估计有限, n_failed=0
//   - Bootstrap SE > 0 (有效变异)
// =============================================================================
TEST(PairedBootstrap, LawSchoolData_EfronTibshirani) {
    auto data = make_law_school_data();
    OLSEstimator ols;
    EstimationResult orig = ols.estimate(data);

    // LSAT 系数 (β1) 应为正 (GPA 随 LSAT 增加而增加)
    EXPECT_GT(orig.coefficients(1), 0.0);

    PairedBootstrap engine;
    BootstrapResult r = engine.resample(ols, data, 999, 42);

    EXPECT_EQ(r.n_failed, 0u);
    EXPECT_EQ(r.bootstrap_samples.size(), 999u);

    // Bootstrap 均值接近原估计
    EXPECT_NEAR(r.coef_mean(0), orig.coefficients(0), 0.1);
    EXPECT_NEAR(r.coef_mean(1), orig.coefficients(1), 0.005);

    // Bootstrap SE > 0
    EXPECT_GT(r.coef_std(0), 0.0);
    EXPECT_GT(r.coef_std(1), 0.0);
    EXPECT_TRUE(std::isfinite(r.lower_ci));
    EXPECT_TRUE(std::isfinite(r.upper_ci));
    EXPECT_LT(r.lower_ci, r.upper_ci);
}

// =============================================================================
// 测试 9: CI 覆盖率 - 蒙特卡洛验证
//   生成 M 个数据集, 每个 Bootstrap 估计 95% CI
//   真值 β1=2.0 应在 ~95% 的 CI 内 (覆盖率 0.93-0.97, 容许抽样波动)
// =============================================================================
namespace {
// 暴露 percentileCI 用于测试 (基类中为 protected)
class TestablePairedBootstrap : public PairedBootstrap {
public:
    using PairedBootstrap::percentileCI;
};
}  // namespace

TEST(PairedBootstrap, CICoverage_LargeSample) {
    const Real beta1_true = 2.0;
    const Size M = 50;  // 蒙特卡洛次数 (受测试时间限制)
    const Size B = 199;  // Bootstrap 次数
    Size covered = 0;

    TestablePairedBootstrap engine;
    for (Size m = 0; m < M; ++m) {
        auto data = make_linear_data(100, 1.0, beta1_true, 0.5, 1000 + m);
        OLSEstimator ols;

        BootstrapResult r = engine.resample(ols, data, B, 42 + m);

        // 提取 β1 (第 1 个系数) 的 CI
        std::vector<Real> beta1_samples;
        for (const auto& s : r.bootstrap_samples) {
            beta1_samples.push_back(s(1));
        }
        if (beta1_samples.empty()) continue;

        auto ci = engine.percentileCI(beta1_samples, 0.05);
        if (ci.first <= beta1_true && beta1_true <= ci.second) {
            ++covered;
        }
    }

    const Real coverage = static_cast<Real>(covered) / static_cast<Real>(M);
    // 95% CI 覆盖率应在 [0.80, 1.00] 范围内 (M=50 抽样波动较大)
    // 注: 严格 [0.93, 0.97] 需要 M >= 1000, 测试时间不允许
    EXPECT_GT(coverage, 0.80);
    EXPECT_LT(coverage, 1.01);
}

// =============================================================================
// 测试 10: 失败 replicate 计数
//   构造一个总抛异常的估计器, 验证 n_failed = B, samples 为空
//   排幻觉点 H-009: 失败必须计入, 不默默忽略
// =============================================================================
namespace {

// 总是抛异常的估计器 (用于测试失败处理)
class AlwaysFailEstimator : public Estimator {
public:
    EstimationResult estimate(const EconData& /*data*/) override {
        throw std::runtime_error("AlwaysFailEstimator: forced failure");
    }
    std::string name() const override { return "AlwaysFail"; }
    EstimatorClass estimatorClass() const override { return EstimatorClass::Parametric; }
    std::unique_ptr<Estimator> clone() const override {
        return std::make_unique<AlwaysFailEstimator>(*this);
    }
};

}  // namespace

TEST(PairedBootstrap, FailedReplicates_Counted) {
    auto data = make_linear_data(20, 1.0, 2.0, 0.5, 42);
    AlwaysFailEstimator fail_est;
    PairedBootstrap engine;

    BootstrapResult r = engine.resample(fail_est, data, 99, 42);

    // 所有 replicate 都失败, n_failed = 99, samples 为空
    EXPECT_EQ(r.n_failed, 99u);
    EXPECT_EQ(r.bootstrap_samples.size(), 0u);
    EXPECT_TRUE(std::isnan(r.lower_ci));
    EXPECT_TRUE(std::isnan(r.upper_ci));
}
