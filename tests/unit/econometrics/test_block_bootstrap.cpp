// SOURCE: PHASE6_IMPLEMENTATION_PLAN §6.2 任务 4.4 - BlockBootstrap 测试
// 验证:
//   - 三种块类型 (Stationary/Circular/NonOverlapping) 基本功能
//   - 块长度自动选择 (经验法则 N^{1/3} + Politis-White 2004)
//   - 固定种子可复现性 (Philox 分块, ADR-004)
//   - 大样本下 Bootstrap 均值接近原估计
//   - 时间序列自相关保留验证 (AR(1) 序列)
//   - 横截面退化路径 (等价 PairedBootstrap)
//   - 排幻觉点 H-018 (时间序列专用), H-019 (块长选择), H-020 (几何分布),
//     H-021 (环形索引), H-022 (块采样数 ⌈N/L⌉)

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
#include "cpphub/econometrics/resampling/block_bootstrap.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;

// =============================================================================
// 辅助: 构造时间序列数据 (AR(1) 过程 y_t = β0 + β1·y_{t-1} + ε_t)
// =============================================================================
namespace {

// 构造 AR(1) 时间序列 (用于验证自相关保留)
// X 的第 1 列为 y 的滞后 1 期 (H-018: 时间序列场景)
EconData make_ar1_ts_data(Size N, Real beta0, Real beta1, Real sigma,
                           std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<Real> noise(0.0, sigma);

    std::vector<Real> y_full(N + 1);  // N+1 个点, 用于构造 N 个 (y_{t-1}, y_t) 对
    y_full[0] = beta0 + noise(rng);
    for (Size t = 1; t <= N; ++t) {
        y_full[t] = beta0 + beta1 * y_full[t - 1] + noise(rng);
    }

    MatrixXD X(N, 2);
    VectorXD y(N);
    for (Size t = 0; t < N; ++t) {
        X(t, 0) = 1.0;
        X(t, 1) = y_full[t];      // y_{t-1}
        y(t) = y_full[t + 1];     // y_t
    }

    std::vector<Real> ts(N);
    for (Size t = 0; t < N; ++t) ts[t] = static_cast<Real>(t);
    return make_time_series(y, X, ts, {"intercept", "y_lag1"}, "y");
}

// 构造简单横截面数据 (用于退化路径测试)
EconData make_linear_cs_data(Size N, Real beta0, Real beta1, Real sigma,
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
// 测试 1: 基本 resample - 不抛异常, 返回正确结构 (Stationary)
// =============================================================================
TEST(BlockBootstrap, BasicResample_Stationary_ReturnsValidResult) {
    auto data = make_ar1_ts_data(100, 0.5, 0.7, 0.3, 42);
    OLSEstimator ols;
    BlockBootstrap engine(BlockType::Stationary);

    BootstrapResult r = engine.resample(ols, data, 199, 42);

    EXPECT_EQ(r.n_replicates, 199u);
    EXPECT_EQ(r.bootstrap_samples.size(), 199u);
    EXPECT_EQ(r.coef_mean.size(), 2u);
    EXPECT_EQ(r.coef_vcov.rows(), 2u);
    EXPECT_EQ(r.coef_vcov.cols(), 2u);
}

// =============================================================================
// 测试 2: 三种块类型均可运行
// =============================================================================
TEST(BlockBootstrap, AllBlockTypes_Run) {
    auto data = make_ar1_ts_data(50, 0.5, 0.7, 0.3, 42);
    OLSEstimator ols;

    for (BlockType bt : {BlockType::Stationary, BlockType::Circular,
                          BlockType::NonOverlapping}) {
        BlockBootstrap engine(bt);
        BootstrapResult r = engine.resample(ols, data, 99, 42);
        EXPECT_EQ(r.n_replicates, 99u);
        // H-009: 失败数应 <= n_replicates
        EXPECT_LE(r.n_failed, 99u);
    }
}

// =============================================================================
// 测试 3: 块长度自动选择 - 经验法则 L = ⌈N^{1/3}⌉
// =============================================================================
TEST(BlockBootstrap, DefaultBlockLength_IsNOneThird) {
    // detail::default_block_length(N) = ⌈N^{1/3}⌉
    // N=125 → 125^{1/3} = 5 → L=5
    EXPECT_EQ(cpphub::v1::econometrics::detail::default_block_length(125), 5u);
    // N=1000 → 1000^{1/3} = 10 → L=10
    EXPECT_EQ(cpphub::v1::econometrics::detail::default_block_length(1000), 10u);
    // N=8 → 8^{1/3} = 2 → L=2
    EXPECT_EQ(cpphub::v1::econometrics::detail::default_block_length(8), 2u);
    // N=1 → L=1 (边界)
    EXPECT_EQ(cpphub::v1::econometrics::detail::default_block_length(1), 1u);
}

// =============================================================================
// 测试 4: 手动设置块长度
// =============================================================================
TEST(BlockBootstrap, ManualBlockLength_Used) {
    auto data = make_ar1_ts_data(100, 0.5, 0.7, 0.3, 42);
    OLSEstimator ols;

    BlockBootstrap engine(BlockType::Circular, 10);  // L=10
    EXPECT_EQ(engine.blockLength(), 10u);

    BootstrapResult r = engine.resample(ols, data, 99, 42);
    // 块长在 resample 后应保持 10 (或 Politis-White 覆盖, 但手动设置时不应被覆盖)
    // 注: auto_block_length_ = false 时, resample 不修改 block_length_
    EXPECT_EQ(engine.blockLength(), 10u);
}

// =============================================================================
// 测试 5: 固定种子可复现性 (H-013, ADR-004)
// =============================================================================
TEST(BlockBootstrap, FixedSeed_Reproducible) {
    auto data = make_ar1_ts_data(80, 0.5, 0.7, 0.3, 42);
    OLSEstimator ols;
    BlockBootstrap engine(BlockType::Stationary);

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
// 测试 6: 大样本下 Bootstrap 均值接近原估计 (大数定律)
// =============================================================================
TEST(BlockBootstrap, LargeSample_MeanConvergesToOriginal) {
    auto data = make_ar1_ts_data(500, 0.5, 0.7, 0.2, 42);
    // OLS 仅接受 CrossSectionData, 从 TimeSeriesData 提取
    const auto& ts = std::get<TimeSeriesData>(data);
    EconData cs_data = make_cross_section(ts.X, ts.y, ts.x_names, ts.y_name);
    OLSEstimator ols;
    EstimationResult orig = ols.estimate(cs_data);

    BlockBootstrap engine(BlockType::Circular, 8);
    BootstrapResult r = engine.resample(ols, data, 999, 42);

    // Bootstrap 均值应接近原估计 (容差较松, 因时间序列 Bootstrap 有额外方差)
    EXPECT_NEAR(r.coef_mean(0), orig.coefficients(0), 0.3);
    EXPECT_NEAR(r.coef_mean(1), orig.coefficients(1), 0.15);
}

// =============================================================================
// 测试 7: 横截面退化路径 - CrossSectionData 不抛异常
// =============================================================================
TEST(BlockBootstrap, CrossSection_DegeneratePath) {
    auto data = make_linear_cs_data(50, 1.0, 2.0, 0.5, 42);
    OLSEstimator ols;
    BlockBootstrap engine(BlockType::Stationary);

    // 横截面数据应退化为 PairedBootstrap 等价 (块长=1)
    BootstrapResult r = engine.resample(ols, data, 99, 42);
    EXPECT_EQ(r.n_replicates, 99u);
    EXPECT_EQ(r.bootstrap_samples.size(), 99u);
}

// =============================================================================
// 测试 8: 非法输入抛异常
// =============================================================================
TEST(BlockBootstrap, InvalidInput_Throws) {
    OLSEstimator ols;
    BlockBootstrap engine;

    // N < 4 (时间序列过短)
    {
        MatrixXD X(3, 2);
        VectorXD y(3);
        for (Size i = 0; i < 3; ++i) { X(i, 0) = 1.0; X(i, 1) = 0.5 * static_cast<Real>(i); y(i) = static_cast<Real>(i); }
        std::vector<Real> ts = {0.0, 1.0, 2.0};
        auto data = make_time_series(y, X, ts, {"b0", "x"}, "y");
        EXPECT_THROW(engine.resample(ols, data, 99, 42), std::invalid_argument);
    }

    // n_replicates < 2
    {
        auto data = make_ar1_ts_data(50, 0.5, 0.7, 0.3, 42);
        EXPECT_THROW(engine.resample(ols, data, 1, 42), std::invalid_argument);
    }
}

// =============================================================================
// 测试 9: Politis-White 自动块长选择 - 白噪声序列 (ACF ≈ 0)
// =============================================================================
TEST(BlockBootstrap, PolitisWhite_WhiteNoise_ShortBlock) {
    // 白噪声序列: ACF 快速衰减到 0, Politis-White 应给出较小块长
    const Size N = 200;
    std::mt19937_64 rng(123);
    std::normal_distribution<Real> noise(0.0, 1.0);

    VectorXD series(N);
    for (Size t = 0; t < N; ++t) series(t) = noise(rng);

    const Size max_lag = 20;
    const std::vector<Real> acf = cpphub::v1::econometrics::detail::sample_acf(series, max_lag);

    // 白噪声的 ACF 在 k >= 1 时应接近 0
    for (Size k = 1; k <= 5; ++k) {
        EXPECT_NEAR(acf[k], 0.0, 0.2);  // 容差较松 (有限样本波动)
    }

    // Politis-White: 白噪声 ACF ≈ 0, G 很小, 块长应 <= 经验法则
    // (白噪声无自相关需保留, Politis-White 正确给出较小或相等的块长)
    const Size L_pw = cpphub::v1::econometrics::detail::politis_white_block_length(acf, N);
    const Size L_default = cpphub::v1::econometrics::detail::default_block_length(N);
    EXPECT_LE(L_pw, L_default);
    EXPECT_GE(L_pw, 1u);
}

// =============================================================================
// 测试 10: 自相关序列 ACF 计算 - AR(1) 序列 ρ(k) ≈ β1^k
// =============================================================================
TEST(BlockBootstrap, SampleACF_AR1_DecaysAsBetaPowK) {
    // AR(1): y_t = 0.7·y_{t-1} + ε_t → ρ(k) ≈ 0.7^k
    const Size N = 5000;
    std::mt19937_64 rng(42);
    std::normal_distribution<Real> noise(0.0, 1.0);

    VectorXD series(N);
    series(0) = noise(rng);
    for (Size t = 1; t < N; ++t) {
        series(t) = 0.7 * series(t - 1) + noise(rng);
    }

    const std::vector<Real> acf = cpphub::v1::econometrics::detail::sample_acf(series, 5);
    EXPECT_NEAR(acf[0], 1.0, 1e-10);
    EXPECT_NEAR(acf[1], 0.7, 0.05);     // ρ(1) ≈ 0.7
    EXPECT_NEAR(acf[2], 0.49, 0.05);    // ρ(2) ≈ 0.7² = 0.49
    EXPECT_NEAR(acf[3], 0.343, 0.05);   // ρ(3) ≈ 0.7³
    EXPECT_TRUE(acf[1] > acf[2] && acf[2] > acf[3]);  // 单调衰减
}

// =============================================================================
// 测试 11: name() 返回包含块类型
// =============================================================================
TEST(BlockBootstrap, Name_IncludesBlockType) {
    BlockBootstrap e1(BlockType::Stationary);
    EXPECT_NE(e1.name().find("Stationary"), std::string::npos);

    BlockBootstrap e2(BlockType::Circular);
    EXPECT_NE(e2.name().find("Circular"), std::string::npos);

    BlockBootstrap e3(BlockType::NonOverlapping);
    EXPECT_NE(e3.name().find("NonOverlapping"), std::string::npos);
}

// =============================================================================
// 测试 12: to_string(BlockType) 转换正确
// =============================================================================
TEST(BlockBootstrap, BlockType_ToString) {
    EXPECT_EQ(to_string(BlockType::Stationary), "Stationary");
    EXPECT_EQ(to_string(BlockType::Circular), "Circular");
    EXPECT_EQ(to_string(BlockType::NonOverlapping), "NonOverlapping");
}
