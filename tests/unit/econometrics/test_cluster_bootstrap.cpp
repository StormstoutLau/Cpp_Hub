// SOURCE: PHASE6_IMPLEMENTATION_PLAN §6.2 任务 4.5 - ClusterBootstrap 测试
// 验证:
//   - Cluster Bootstrap 基本功能 (PanelData, 聚类=entity)
//   - 显式 cluster_id (CrossSectionData + setClusterId)
//   - 固定种子可复现性
//   - 大样本下 Bootstrap 均值接近原估计
//   - G < 20 警告标志
//   - 不等大小聚类处理 (H-026)
//   - 排幻觉点 H-023 (聚类级采样), H-024 (整体保留), H-025 (G<20 警告), H-026 (N*≠N)

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
#include "cpphub/econometrics/resampling/cluster_bootstrap.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;
using cpphub::v1::Index;

// =============================================================================
// 辅助: 构造平衡面板数据 (G 个 entity, T 期, y = β0 + β1·x + ε, 聚类内相关)
// =============================================================================
namespace {

// 构造平衡面板: G entities × T periods = N obs
// 聚类内误差相关 (entity 固定效应 α_g ~ N(0, σ_α²))
EconData make_panel_data(Size G, Size T, Real beta0, Real beta1,
                          Real sigma_eps, Real sigma_alpha, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<Real> eps(0.0, sigma_eps);
    std::normal_distribution<Real> alpha(0.0, sigma_alpha);
    std::uniform_real_distribution<Real> xdist(0.0, 10.0);

    const Size N = G * T;
    MatrixXD X(N, 2);
    VectorXD y(N);
    std::vector<Index> entity_id(N);
    std::vector<Index> time_id(N);

    Size row = 0;
    for (Size g = 0; g < G; ++g) {
        const Real ag = alpha(rng);  // entity 固定效应
        for (Size t = 0; t < T; ++t) {
            const Real xi = xdist(rng);
            X(row, 0) = 1.0;
            X(row, 1) = xi;
            y(row) = beta0 + beta1 * xi + ag + eps(rng);
            entity_id[row] = static_cast<Index>(g);
            time_id[row] = static_cast<Index>(t);
            ++row;
        }
    }

    return make_panel(X, y, entity_id, time_id, {"intercept", "x"}, "y", true);
}

// 构造横截面数据 + 显式 cluster_id (G 个聚类, 不等大小)
EconData make_cs_with_clusters(Size N, Real beta0, Real beta1, Real sigma,
                                 const std::vector<Index>& cluster_id,
                                 std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<Real> noise(0.0, sigma);
    std::uniform_real_distribution<Real> xdist(0.0, 10.0);

    MatrixXD X(N, 2);
    VectorXD y(N);
    for (Size i = 0; i < N; ++i) {
        const Real xi = xdist(rng);
        X(i, 0) = 1.0;
        X(i, 1) = xi;
        y(i) = beta0 + beta1 * xi + noise(rng);
    }
    return make_cross_section(X, y, {"intercept", "x"}, "y");
}

}  // namespace

// =============================================================================
// 测试 1: 基本 resample - PanelData, 聚类=entity, 不抛异常
// =============================================================================
TEST(ClusterBootstrap, BasicResample_PanelData_ReturnsValidResult) {
    auto data = make_panel_data(30, 5, 1.0, 2.0, 0.5, 0.3, 42);  // G=30, T=5, N=150
    OLSEstimator ols;
    ClusterBootstrap engine;

    BootstrapResult r = engine.resample(ols, data, 199, 42);

    EXPECT_EQ(r.n_replicates, 199u);
    EXPECT_EQ(r.bootstrap_samples.size(), 199u);
    EXPECT_EQ(r.coef_mean.size(), 2u);
    EXPECT_EQ(r.coef_vcov.rows(), 2u);
    EXPECT_EQ(r.coef_vcov.cols(), 2u);
    EXPECT_EQ(engine.nClusters(), 30u);
}

// =============================================================================
// 测试 2: CrossSectionData + 显式 cluster_id
// =============================================================================
TEST(ClusterBootstrap, CrossSection_WithExplicitClusterId) {
    const Size N = 100;
    // 20 个聚类, 每个聚类 5 个观测
    std::vector<Index> cluster_id(N);
    for (Size i = 0; i < N; ++i) {
        cluster_id[i] = static_cast<Index>(i / 5);  // 聚类 0..19
    }

    auto data = make_cs_with_clusters(N, 1.0, 2.0, 0.5, cluster_id, 42);
    OLSEstimator ols;
    ClusterBootstrap engine;
    engine.setClusterId(cluster_id);

    BootstrapResult r = engine.resample(ols, data, 99, 42);

    EXPECT_EQ(r.n_replicates, 99u);
    EXPECT_EQ(r.bootstrap_samples.size(), 99u);
    EXPECT_EQ(engine.nClusters(), 20u);
}

// =============================================================================
// 测试 3: CrossSectionData 无 cluster_id 抛异常
// =============================================================================
TEST(ClusterBootstrap, CrossSection_NoClusterId_Throws) {
    auto data = make_cs_with_clusters(50, 1.0, 2.0, 0.5, {}, 42);
    OLSEstimator ols;
    ClusterBootstrap engine;

    EXPECT_THROW(engine.resample(ols, data, 99, 42), std::invalid_argument);
}

// =============================================================================
// 测试 4: 固定种子可复现性 (H-013, ADR-004)
// =============================================================================
TEST(ClusterBootstrap, FixedSeed_Reproducible) {
    auto data = make_panel_data(25, 4, 1.0, 2.0, 0.5, 0.3, 42);
    OLSEstimator ols;
    ClusterBootstrap engine;

    BootstrapResult r1 = engine.resample(ols, data, 99, 42);
    BootstrapResult r2 = engine.resample(ols, data, 99, 42);

    ASSERT_EQ(r1.bootstrap_samples.size(), r2.bootstrap_samples.size());
    for (Size b = 0; b < r1.bootstrap_samples.size(); ++b) {
        for (Size i = 0; i < r1.bootstrap_samples[b].size(); ++i) {
            EXPECT_NEAR(r1.bootstrap_samples[b](i),
                        r2.bootstrap_samples[b](i), 1e-12);
        }
    }
}

// =============================================================================
// 测试 5: 大样本下 Bootstrap 均值接近原估计 (H-023, H-024)
// =============================================================================
TEST(ClusterBootstrap, LargeSample_MeanConvergesToOriginal) {
    // G=50, T=10, N=500, 聚类数足够大
    auto data = make_panel_data(50, 10, 1.0, 2.0, 0.3, 0.2, 42);
    // OLS 仅接受 CrossSectionData, 从 PanelData 提取
    const auto& panel = std::get<PanelData>(data);
    EconData cs_data = make_cross_section(panel.X, panel.y, panel.x_names, panel.y_name);
    OLSEstimator ols;
    EstimationResult orig = ols.estimate(cs_data);

    ClusterBootstrap engine;
    BootstrapResult r = engine.resample(ols, data, 999, 42);

    // Bootstrap 均值应接近原估计
    EXPECT_NEAR(r.coef_mean(0), orig.coefficients(0), 0.2);
    EXPECT_NEAR(r.coef_mean(1), orig.coefficients(1), 0.1);
}

// =============================================================================
// 测试 6: G < 20 警告标志 (H-025)
// =============================================================================
TEST(ClusterBootstrap, SmallClusterCount_SetsWarning) {
    // G=10 < 20, 应触发警告
    auto data = make_panel_data(10, 5, 1.0, 2.0, 0.5, 0.3, 42);
    OLSEstimator ols;
    ClusterBootstrap engine;

    EXPECT_FALSE(engine.smallClusterWarning());  // resample 前为 false

    BootstrapResult r = engine.resample(ols, data, 99, 42);

    EXPECT_TRUE(engine.smallClusterWarning());   // resample 后为 true
    EXPECT_EQ(engine.nClusters(), 10u);
}

// =============================================================================
// 测试 7: G >= 20 无警告
// =============================================================================
TEST(ClusterBootstrap, LargeClusterCount_NoWarning) {
    auto data = make_panel_data(30, 5, 1.0, 2.0, 0.5, 0.3, 42);
    OLSEstimator ols;
    ClusterBootstrap engine;

    BootstrapResult r = engine.resample(ols, data, 99, 42);

    EXPECT_FALSE(engine.smallClusterWarning());
    EXPECT_EQ(engine.nClusters(), 30u);
}

// =============================================================================
// 测试 8: 不等大小聚类 (H-026: N* ≠ N)
// =============================================================================
TEST(ClusterBootstrap, UnequalClusterSizes_HandlesCorrectly) {
    // 构造不等大小聚类: 聚类 0 有 3 个观测, 聚类 1 有 7 个观测, ... 共 20 个聚类
    const Size N = 100;
    std::vector<Index> cluster_id(N);
    std::vector<Size> cluster_sizes(20);
    for (Size i = 0; i < 20; ++i) cluster_sizes[i] = 3 + (i % 5);  // 大小 3,4,5,6,7 循环
    Size pos = 0;
    for (Size g = 0; g < 20; ++g) {
        for (Size j = 0; j < cluster_sizes[g]; ++j) {
            cluster_id[pos++] = static_cast<Index>(g);
        }
    }

    auto data = make_cs_with_clusters(N, 1.0, 2.0, 0.5, cluster_id, 42);
    OLSEstimator ols;
    ClusterBootstrap engine;
    engine.setClusterId(cluster_id);

    BootstrapResult r = engine.resample(ols, data, 99, 42);

    // H-026: 不等大小聚类时应正常运行, 失败数应合理
    EXPECT_EQ(r.n_replicates, 99u);
    // 大部分 replicate 应成功 (可能因不等大小导致个别失败)
    EXPECT_GT(r.bootstrap_samples.size(), 80u);
    EXPECT_EQ(engine.nClusters(), 20u);
}

// =============================================================================
// 测试 9: G=1 抛异常 (无法重采样)
// =============================================================================
TEST(ClusterBootstrap, SingleCluster_Throws) {
    const Size N = 20;
    std::vector<Index> cluster_id(N, 0);  // 全部属于聚类 0
    auto data = make_cs_with_clusters(N, 1.0, 2.0, 0.5, cluster_id, 42);
    OLSEstimator ols;
    ClusterBootstrap engine;
    engine.setClusterId(cluster_id);

    EXPECT_THROW(engine.resample(ols, data, 99, 42), std::invalid_argument);
}

// =============================================================================
// 测试 10: cluster_id 长度不匹配抛异常
// =============================================================================
TEST(ClusterBootstrap, ClusterIdSizeMismatch_Throws) {
    auto data = make_cs_with_clusters(50, 1.0, 2.0, 0.5, {}, 42);
    OLSEstimator ols;
    ClusterBootstrap engine;
    engine.setClusterId({0, 1, 2});  // 长度 3 != N=50

    EXPECT_THROW(engine.resample(ols, data, 99, 42), std::invalid_argument);
}

// =============================================================================
// 测试 11: 聚类维度="time" (PanelData)
// =============================================================================
TEST(ClusterBootstrap, ClusterDimension_Time) {
    // G=10 entities × T=25 periods, 按 time 聚类 → G=25 个聚类
    auto data = make_panel_data(10, 25, 1.0, 2.0, 0.5, 0.3, 42);
    OLSEstimator ols;
    ClusterBootstrap engine;
    engine.setClusterDimension("time");

    BootstrapResult r = engine.resample(ols, data, 99, 42);

    EXPECT_EQ(engine.nClusters(), 25u);
    EXPECT_FALSE(engine.smallClusterWarning());  // 25 >= 20
    EXPECT_EQ(r.bootstrap_samples.size(), 99u);
}

// =============================================================================
// 测试 12: 无效聚类维度抛异常
// =============================================================================
TEST(ClusterBootstrap, InvalidClusterDimension_Throws) {
    ClusterBootstrap engine;
    EXPECT_THROW(engine.setClusterDimension("invalid"), std::invalid_argument);
}

// =============================================================================
// 测试 13: 非法 cluster_id 值 (非 0-based) 也能正确处理 (排幻觉点 P2)
// =============================================================================
TEST(ClusterBootstrap, NonZeroBasedClusterId_HandledCorrectly) {
    // 原始 ID 从 100 开始 (非 0-based), 应被重映射为 0-based 紧凑索引
    const Size N = 100;
    std::vector<Index> cluster_id(N);
    for (Size i = 0; i < N; ++i) {
        cluster_id[i] = static_cast<Index>(100 + i / 5);  // 聚类 100..119
    }

    auto data = make_cs_with_clusters(N, 1.0, 2.0, 0.5, cluster_id, 42);
    OLSEstimator ols;
    ClusterBootstrap engine;
    engine.setClusterId(cluster_id);

    BootstrapResult r = engine.resample(ols, data, 99, 42);

    // 应正确识别 20 个聚类 (100..119 重映射为 0..19)
    EXPECT_EQ(engine.nClusters(), 20u);
    EXPECT_EQ(r.bootstrap_samples.size(), 99u);
}

// =============================================================================
// 测试 14: name() 返回正确
// =============================================================================
TEST(ClusterBootstrap, Name_Correct) {
    ClusterBootstrap engine;
    EXPECT_EQ(engine.name(), "ClusterBootstrap");
}
