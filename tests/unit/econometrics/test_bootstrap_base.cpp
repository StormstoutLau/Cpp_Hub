// SOURCE: PHASE6_IMPLEMENTATION_PLAN §6.1 任务 4.1 - Bootstrap 基类编译验证
// 验证: BootstrapEngine 抽象基类可编译 + 工具方法正确性
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

#include "cpphub/econometrics/resampling/bootstrap_base.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;

// =============================================================================
// 辅助: 具体引擎 (用于测试基类工具方法)
// =============================================================================
class TestEngine : public BootstrapEngine {
public:
    BootstrapResult resample(Estimator& /*estimator*/, const EconData& /*data*/,
                              Size /*n_replicates*/, std::uint64_t /*seed*/) override {
        BootstrapResult r;
        r.n_replicates = 0;
        return r;
    }
    std::string name() const override { return "TestEngine"; }

    // 暴露基类工具方法用于测试
    using BootstrapEngine::percentileCI;
    using BootstrapEngine::bcaCI;
    using BootstrapEngine::randint;
    using BootstrapEngine::compute_stats;
    using BootstrapEngine::extract_cross_section;
    using BootstrapEngine::make_cross_section_data;
};

// =============================================================================
// 测试 1: percentileCI - 已知样本的分位数
// =============================================================================
TEST(BootstrapBase, PercentileCI_KnownValues) {
    TestEngine engine;
    // 样本: 1, 2, 3, ..., 100
    std::vector<Real> samples;
    for (Size i = 1; i <= 100; ++i) samples.push_back(static_cast<Real>(i));

    auto [lo, hi] = engine.percentileCI(samples, 0.05);
    // 95% CI: 2.5% 和 97.5% 分位数
    // R type 7: h = 99 * 0.025 = 2.475 → lo = 3.475
    //           h = 99 * 0.975 = 96.525 → hi = 97.525
    EXPECT_NEAR(lo, 3.475, 1e-6);
    EXPECT_NEAR(hi, 97.525, 1e-6);
}

// =============================================================================
// 测试 2: percentileCI - 空样本返回 NaN
// =============================================================================
TEST(BootstrapBase, PercentileCI_Empty) {
    TestEngine engine;
    auto [lo, hi] = engine.percentileCI({}, 0.05);
    EXPECT_TRUE(std::isnan(lo));
    EXPECT_TRUE(std::isnan(hi));
}

// =============================================================================
// 测试 3: randint - 均匀分布范围 [0, n)
// =============================================================================
TEST(BootstrapBase, RandInt_Range) {
    TestEngine engine;
    cpphub::v1::Philox4x64 rng(42);
    const Size n = 100;
    const Size trials = 10000;
    for (Size i = 0; i < trials; ++i) {
        Size idx = engine.randint(rng, n);
        EXPECT_LT(idx, n);
    }
}

// =============================================================================
// 测试 4: randint - 固定种子可复现
// =============================================================================
TEST(BootstrapBase, RandInt_Reproducible) {
    TestEngine engine;
    cpphub::v1::Philox4x64 rng1(42);
    cpphub::v1::Philox4x64 rng2(42);
    const Size n = 1000;
    for (Size i = 0; i < 100; ++i) {
        EXPECT_EQ(engine.randint(rng1, n), engine.randint(rng2, n));
    }
}

// =============================================================================
// 测试 5: compute_stats - 均值和协方差
// =============================================================================
TEST(BootstrapBase, ComputeStats_MeanAndCov) {
    TestEngine engine;
    // 3 个样本, 2 维
    std::vector<cpphub::v1::linalg::dynamic::VectorXD> samples;
    cpphub::v1::linalg::dynamic::VectorXD s1(2), s2(2), s3(2);
    s1(0) = 1.0; s1(1) = 2.0;
    s2(0) = 3.0; s2(1) = 4.0;
    s3(0) = 5.0; s3(1) = 6.0;
    samples.push_back(s1);
    samples.push_back(s2);
    samples.push_back(s3);

    cpphub::v1::linalg::dynamic::VectorXD mean;
    cpphub::v1::linalg::dynamic::MatrixXD vcov;
    engine.compute_stats(samples, mean, vcov);

    // 均值 = [3.0, 4.0]
    EXPECT_NEAR(mean(0), 3.0, 1e-10);
    EXPECT_NEAR(mean(1), 4.0, 1e-10);

    // 协方差: (1/(B-1)) Σ (s_b - mean)(s_b - mean)'
    // B-1 = 2, 差值: [-2,-2], [0,0], [2,2]
    // V = (1/2) * ([4,4;4,4] + [0,0;0,0] + [4,4;4,4]) = [4,4;4,4]
    EXPECT_NEAR(vcov(0, 0), 4.0, 1e-10);
    EXPECT_NEAR(vcov(0, 1), 4.0, 1e-10);
    EXPECT_NEAR(vcov(1, 0), 4.0, 1e-10);
    EXPECT_NEAR(vcov(1, 1), 4.0, 1e-10);
}

// =============================================================================
// 测试 6: BCa CI - 对称分布近似 percentile CI
// =============================================================================
TEST(BootstrapBase, BCaCI_Symmetric) {
    TestEngine engine;
    // 对称正态样本: z0 ≈ 0, a ≈ 0 → BCa ≈ percentile
    std::vector<Real> samples;
    for (Size i = 1; i <= 1000; ++i) {
        samples.push_back(static_cast<Real>(i));
    }
    Real theta_hat = 500.5;  // 中位数

    // Jackknife: 对称样本, a ≈ 0
    std::vector<Real> jk;
    for (Size i = 0; i < 100; ++i) {
        jk.push_back(500.5 + (i % 2 == 0 ? 0.1 : -0.1));
    }

    auto [p_lo, p_hi] = engine.percentileCI(samples, 0.05);
    auto [b_lo, b_hi] = engine.bcaCI(samples, 0.05, jk, theta_hat);

    // BCa 应接近 percentile (对称分布)
    EXPECT_NEAR(b_lo, p_lo, 5.0);
    EXPECT_NEAR(b_hi, p_hi, 5.0);
}

// =============================================================================
// 测试 7: extract_cross_section - 正确提取
// =============================================================================
TEST(BootstrapBase, ExtractCrossSection) {
    using cpphub::v1::linalg::dynamic::MatrixXD;
    using cpphub::v1::linalg::dynamic::VectorXD;
    TestEngine engine;
    MatrixXD X(3, 2);
    VectorXD y(3);
    X(0, 0) = 1.0; X(0, 1) = 2.0;
    X(1, 0) = 3.0; X(1, 1) = 4.0;
    X(2, 0) = 5.0; X(2, 1) = 6.0;
    y(0) = 10.0; y(1) = 20.0; y(2) = 30.0;
    EconData data = make_cross_section(X, y, {"x1", "x2"}, "y");

    MatrixXD X_out;
    VectorXD y_out;
    engine.extract_cross_section(data, X_out, y_out);

    EXPECT_EQ(X_out.rows(), 3);
    EXPECT_EQ(X_out.cols(), 2);
    EXPECT_EQ(y_out.size(), 3);
    EXPECT_NEAR(X_out(0, 0), 1.0, 1e-10);
    EXPECT_NEAR(y_out(2), 30.0, 1e-10);
}

// =============================================================================
// 测试 8: extract_cross_section - 非 CrossSectionData 抛异常
// =============================================================================
TEST(BootstrapBase, ExtractCrossSection_WrongType) {
    TestEngine engine;
    PanelData panel;
    panel.X = cpphub::v1::linalg::dynamic::MatrixXD(2, 1);
    panel.y = cpphub::v1::linalg::dynamic::VectorXD(2);
    panel.entity_id = {0, 0};
    panel.time_id = {0, 1};
    panel.balanced = true;
    EconData data = panel;

    MatrixXD X_out;
    VectorXD y_out;
    EXPECT_THROW(engine.extract_cross_section(data, X_out, y_out),
                 std::invalid_argument);
}
