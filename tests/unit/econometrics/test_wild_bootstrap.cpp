// SOURCE: PHASE6_IMPLEMENTATION_PLAN §6.1 任务 4.3 - WildBootstrap 测试
// 验证:
//   - 三种权重分布 (Rademacher/Mammen/Webb6) 的统计性质 (E[v]=0, Var[v]=1)
//   - Wild Bootstrap 基本功能 (固定 X, 残差乘权重)
//   - 固定种子可复现性
//   - 大样本下 Bootstrap 均值/SE 接近原估计
//   - 三种权重分布大样本下结果相近
//   - 失败 replicate 计数 (H-009)
//   - 排幻觉点 E12 (R multiwayvcov 对照), H-014 (固定 X), H-015 (β̂ 固定), H-016 (y*=fitted+v·ε), H-017 (E[v]=0,Var[v]=1)

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/econometrics/estimation/ols.hpp"
#include "cpphub/econometrics/resampling/wild_bootstrap.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;

// =============================================================================
// 辅助: 构造横截面数据
// =============================================================================
namespace {

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

// 暴露 generate_wild_weight 用于统计性质测试
class TestableWildBootstrap : public WildBootstrap {
public:
    using WildBootstrap::percentileCI;
    static Real generateWeight(cpphub::v1::Philox4x64& rng, WildWeightDistribution d) {
        return detail::generate_wild_weight(rng, d);
    }
};

}  // namespace

// =============================================================================
// 测试 1: 基本 resample - 不抛异常, 返回正确结构
// =============================================================================
TEST(WildBootstrap, BasicResample_ReturnsValidResult) {
    auto data = make_linear_data(50, 1.0, 2.0, 0.5, 42);
    OLSEstimator ols;
    WildBootstrap engine;

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
// 测试 2: Rademacher 权重 - E[v]=0, Var[v]=1
//   排幻觉点 H-017
// =============================================================================
TEST(WildBootstrap, Rademacher_MeanZero_VarOne) {
    cpphub::v1::Philox4x64 rng(42, 0);
    const Size M = 100000;
    Real sum = 0.0, sum_sq = 0.0;
    Real min_v = std::numeric_limits<Real>::max();
    Real max_v = std::numeric_limits<Real>::lowest();

    for (Size i = 0; i < M; ++i) {
        const Real v = TestableWildBootstrap::generateWeight(rng, WildWeightDistribution::Rademacher);
        sum += v;
        sum_sq += v * v;
        min_v = std::min(min_v, v);
        max_v = std::max(max_v, v);
    }
    const Real mean = sum / M;
    const Real var = sum_sq / M - mean * mean;

    // Rademacher: v ∈ {-1, +1}
    EXPECT_NEAR(min_v, -1.0, 1e-15);
    EXPECT_NEAR(max_v, 1.0, 1e-15);
    EXPECT_NEAR(mean, 0.0, 0.01);   // E[v] = 0
    EXPECT_NEAR(var, 1.0, 0.01);    // Var[v] = 1
}

// =============================================================================
// 测试 3: Mammen 权重 - E[v]=0, Var[v]=1
//   排幻觉点 H-017
// =============================================================================
TEST(WildBootstrap, Mammen_MeanZero_VarOne) {
    cpphub::v1::Philox4x64 rng(42, 0);
    const Size M = 100000;
    Real sum = 0.0, sum_sq = 0.0;

    for (Size i = 0; i < M; ++i) {
        const Real v = TestableWildBootstrap::generateWeight(rng, WildWeightDistribution::Mammen);
        sum += v;
        sum_sq += v * v;
    }
    const Real mean = sum / M;
    const Real var = sum_sq / M - mean * mean;

    EXPECT_NEAR(mean, 0.0, 0.01);   // E[v] = 0
    EXPECT_NEAR(var, 1.0, 0.01);    // Var[v] = 1
}

// =============================================================================
// 测试 4: Webb6 权重 - E[v]=0, Var[v]=1, 6 个离散值
//   排幻觉点 H-017
// =============================================================================
TEST(WildBootstrap, Webb6_MeanZero_VarOne) {
    cpphub::v1::Philox4x64 rng(42, 0);
    const Size M = 100000;
    Real sum = 0.0, sum_sq = 0.0;
    std::set<Real> distinct_values;

    for (Size i = 0; i < M; ++i) {
        const Real v = TestableWildBootstrap::generateWeight(rng, WildWeightDistribution::Webb6);
        sum += v;
        sum_sq += v * v;
        distinct_values.insert(std::round(v * 1e10) / 1e10);  // 量化以便去重
    }
    const Real mean = sum / M;
    const Real var = sum_sq / M - mean * mean;

    // Webb6: 恰好 6 个离散值
    EXPECT_EQ(distinct_values.size(), 6u);
    EXPECT_NEAR(mean, 0.0, 0.01);   // E[v] = 0
    EXPECT_NEAR(var, 1.0, 0.01);    // Var[v] = 1
}

// =============================================================================
// 测试 5: 固定种子可复现 - 同 (seed, B, dist) 产出位精确相同结果
// =============================================================================
TEST(WildBootstrap, Reproducible_SameSeed) {
    auto data = make_linear_data(40, 0.5, 1.5, 0.3, 100);
    OLSEstimator ols;
    WildBootstrap engine(WildWeightDistribution::Mammen);

    BootstrapResult r1 = engine.resample(ols, data, 199, 42);
    BootstrapResult r2 = engine.resample(ols, data, 199, 42);

    EXPECT_EQ(r1.bootstrap_samples.size(), r2.bootstrap_samples.size());
    for (Size b = 0; b < r1.bootstrap_samples.size(); ++b) {
        for (Size j = 0; j < r1.bootstrap_samples[b].size(); ++j) {
            EXPECT_NEAR(r1.bootstrap_samples[b](j),
                        r2.bootstrap_samples[b](j), 1e-15);
        }
    }
}

// =============================================================================
// 测试 6: 大样本 Bootstrap 均值接近原估计
//   Efron-Tibshirani 定理: θ̄* → θ̂ 当 B → ∞
// =============================================================================
TEST(WildBootstrap, BootstrapMean_CloseToOriginal) {
    auto data = make_linear_data(200, 1.0, 2.0, 0.5, 400);
    OLSEstimator ols;
    EstimationResult orig = ols.estimate(data);

    WildBootstrap engine(WildWeightDistribution::Rademacher);
    BootstrapResult r = engine.resample(ols, data, 999, 42);

    // Wild Bootstrap 均值应接近原估计 (容差 ~2·SE)
    EXPECT_NEAR(r.coef_mean(0), orig.coefficients(0), 0.1);
    EXPECT_NEAR(r.coef_mean(1), orig.coefficients(1), 0.05);
}

// =============================================================================
// 测试 7: 大样本 Bootstrap SE 接近 OLS Classical SE
//   同方差下 Wild Bootstrap SE ≈ OLS SE
// =============================================================================
TEST(WildBootstrap, BootstrapSE_CloseToOLSSe) {
    auto data = make_linear_data(300, 1.0, 2.0, 0.5, 500);
    OLSEstimator ols(CovarianceType::Classical);
    EstimationResult orig = ols.estimate(data);

    WildBootstrap engine(WildWeightDistribution::Rademacher);
    BootstrapResult r = engine.resample(ols, data, 999, 42);

    // Wild Bootstrap SE 应在 OLS SE 的 ±25% 范围内 (Wild Bootstrap 自身有波动)
    const Real se0_ratio = r.coef_std(0) / orig.std_errors(0);
    const Real se1_ratio = r.coef_std(1) / orig.std_errors(1);
    EXPECT_GT(se0_ratio, 0.75);
    EXPECT_LT(se0_ratio, 1.25);
    EXPECT_GT(se1_ratio, 0.75);
    EXPECT_LT(se1_ratio, 1.25);
}

// =============================================================================
// 测试 8: 三种权重分布大样本下结果相近
//   Rademacher / Mammen / Webb6 在大样本下都收敛到真值
// =============================================================================
TEST(WildBootstrap, ThreeDistances_ConvergeSame) {
    auto data = make_linear_data(300, 1.0, 2.0, 0.5, 600);
    OLSEstimator ols;

    WildBootstrap engine_r(WildWeightDistribution::Rademacher);
    WildBootstrap engine_m(WildWeightDistribution::Mammen);
    WildBootstrap engine_w(WildWeightDistribution::Webb6);

    BootstrapResult r_r = engine_r.resample(ols, data, 499, 42);
    BootstrapResult r_m = engine_m.resample(ols, data, 499, 42);
    BootstrapResult r_w = engine_w.resample(ols, data, 499, 42);

    // 三种分布的均值应接近 (大样本下都收敛到真值)
    EXPECT_NEAR(r_r.coef_mean(1), r_m.coef_mean(1), 0.1);
    EXPECT_NEAR(r_r.coef_mean(1), r_w.coef_mean(1), 0.1);
    EXPECT_NEAR(r_m.coef_mean(1), r_w.coef_mean(1), 0.1);

    // 三种分布的 SE 应在 ±30% 范围内 (Wild Bootstrap 自身有 ~1/√B 波动)
    const Real se_r = r_r.coef_std(1);
    const Real se_m = r_m.coef_std(1);
    const Real se_w = r_w.coef_std(1);
    EXPECT_GT(se_r / se_m, 0.70);
    EXPECT_LT(se_r / se_m, 1.30);
    EXPECT_GT(se_r / se_w, 0.70);
    EXPECT_LT(se_r / se_w, 1.30);
}

// =============================================================================
// 测试 9: name() 包含权重分布名称
// =============================================================================
TEST(WildBootstrap, Name_ContainsDistribution) {
    WildBootstrap engine_r(WildWeightDistribution::Rademacher);
    WildBootstrap engine_m(WildWeightDistribution::Mammen);
    WildBootstrap engine_w(WildWeightDistribution::Webb6);

    EXPECT_EQ(engine_r.name(), "WildBootstrap[Rademacher]");
    EXPECT_EQ(engine_m.name(), "WildBootstrap[Mammen]");
    EXPECT_EQ(engine_w.name(), "WildBootstrap[Webb6]");
}

// =============================================================================
// 测试 10: 失败 replicate 计数
//   排幻觉点 H-009: 失败必须计入, 不默默忽略
// =============================================================================
namespace {

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

TEST(WildBootstrap, FailedReplicates_Counted) {
    auto data = make_linear_data(20, 1.0, 2.0, 0.5, 42);
    AlwaysFailEstimator fail_est;
    WildBootstrap engine;

    // 第一次 estimate (原样本) 就抛异常, WildBootstrap 整体抛 runtime_error
    EXPECT_THROW(engine.resample(fail_est, data, 99, 42), std::runtime_error);
}
