// =============================================================================
// test_integration_phase9.cpp - Phase 9 端到端: KDE -> KSG -> TE + 7C Granger 互检
// (TDD 前置)
//
// 状态: Phase 9 M1-M3 已实现, 已接入 tests/CMakeLists.txt.
//
// 验收:
//   1. KDE(Scott) 密度积分 ≈ 1
//   2. KSG MI vs Gaussian 解析 (rho=0.8) 偏差 ≤ 0.05 nat
//   3. TE(x->y) vs GC/2, 3 seeds 均值 ±35%
//   4. TE/GC_hat ≈ 0.5, 容差 ±0.15 (0.35..0.65)
//      GC_hat = ln(1 + F/df2) 由 7C Granger 的 SSR-F 换算
// =============================================================================

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <span>
#include <vector>

#include "cpphub/timeseries/arima/granger_test.hpp"
#include "cpphub/timeseries/density/kde.hpp"
#include "cpphub/timeseries/information/ksg.hpp"
#include "cpphub/timeseries/information/transfer_entropy.hpp"

namespace de = cpphub::v1::timeseries::density;
namespace inf = cpphub::v1::timeseries::information;
namespace gr = cpphub::v1::timeseries::granger;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

struct VarSample {
    std::vector<Real> x;
    std::vector<Real> y;
};

VarSample generate_var(Size n, Real phi, Real a, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<Real> nd(0.0, 1.0);
    VarSample s;
    s.x.resize(n);
    s.y.resize(n);
    s.x[0] = nd(rng) / std::sqrt(1.0 - phi * phi);
    s.y[0] = nd(rng) / std::sqrt(1.0 - phi * phi);
    for (Size t = 1; t < n; ++t) {
        s.x[t] = phi * s.x[t - 1] + nd(rng);
        s.y[t] = phi * s.y[t - 1] + a * s.x[t - 1] + nd(rng);
    }
    return s;
}

struct Pair {
    std::vector<Real> x;
    std::vector<Real> y;
};

Pair gaussian_pair(Size n, Real rho, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<Real> nd(0.0, 1.0);
    Pair p;
    p.x.resize(n);
    p.y.resize(n);
    const Real s = std::sqrt(1.0 - rho * rho);
    for (Size i = 0; i < n; ++i) {
        const Real z1 = nd(rng);
        const Real z2 = nd(rng);
        p.x[i] = z1;
        p.y[i] = rho * z1 + s * z2;
    }
    return p;
}

de::KdeConfig scott_config() {
    de::KdeConfig cfg;
    cfg.bandwidth.rule = de::BandwidthRule::kScott;
    return cfg;
}

Real integrate_pdf(const de::KernelDensity1D& kde, Real lo, Real hi, Size steps) {
    const Real dx = (hi - lo) / static_cast<Real>(steps);
    Real sum = 0.0;
    for (Size i = 0; i <= steps; ++i) {
        const Real x = lo + dx * static_cast<Real>(i);
        const Real w = (i == 0 || i == steps) ? 0.5 : 1.0;
        sum += w * kde.pdf(x);
    }
    return sum * dx;
}

}  // namespace

TEST(Phase9Integration, KdeKsgTePipeline) {
    constexpr Size n = 2000;
    const VarSample v = generate_var(n, 0.5, 0.3, 900);

    // 1. KDE(Scott) 密度积分 ≈ 1
    const de::KernelDensity1D kde(std::span<const Real>(v.x), scott_config());
    EXPECT_GT(kde.bandwidth(), 0.0);
    EXPECT_NEAR(integrate_pdf(kde, -6.0, 6.0, 2400), 1.0, 1e-3);

    // 2. KSG MI vs Gaussian 解析 (纯双变量 Gaussian fixture, 不从 VAR 反解 rho;
    //    spec §6.1.3: 20 seeds 均值偏差 ≤ 0.02 nat, 审计 B4)
    Real mi_sum = 0.0;
    constexpr Size kMiSeeds = 20;
    for (Size s = 0; s < kMiSeeds; ++s) {
        const Pair g = gaussian_pair(n, 0.8, 901 + s);
        inf::KsgConfig ksg_cfg;
        ksg_cfg.k = 3;
        ksg_cfg.standardize = true;
        ksg_cfg.jitter_scale = 0.0;
        mi_sum += inf::ksg_mutual_info(
            std::span<const Real>(g.x), std::span<const Real>(g.y), ksg_cfg);
    }
    const Real mi = mi_sum / static_cast<Real>(kMiSeeds);
    EXPECT_NEAR(mi, 0.5108256237659907, 0.02) << "mean MI=" << mi;

    // 3. TE(x->y) vs GC/2 (3 seeds 均值)
    constexpr Real gc = 0.11332868530700327;
    constexpr Real te_true = gc / 2.0;
    inf::TransferEntropyConfig te_cfg;
    te_cfg.source_history = 1;
    te_cfg.target_history = 1;
    te_cfg.delay = 1;
    te_cfg.ksg_neighbors = 3;
    te_cfg.standardize = true;
    te_cfg.n_shuffles = 10;
    te_cfg.n_block_boot = 10;
    te_cfg.seed = 902;

    Real te_sum = 0.0;
    for (Size s = 0; s < 3; ++s) {
        const VarSample vs = generate_var(n, 0.5, 0.3, 910 + s);
        te_sum += inf::transfer_entropy(
                      std::span<const Real>(vs.x), std::span<const Real>(vs.y), te_cfg)
                      .te;
    }
    const Real te = te_sum / 3.0;
    // spec §6.1.4/§6.3: ±30% (probe K4 bias +25%, 固定 seed 确定性)
    EXPECT_NEAR(te, te_true, 0.30 * te_true) << "TE=" << te << " GC/2=" << te_true;

    // 4. 7C Granger 内部互检: GC_hat = ln(1 + F/df2)
    const auto gres = gr::granger_test(v.x, v.y, 1, 0, false);
    ASSERT_TRUE(std::isfinite(gres.f_stat));
    ASSERT_GT(gres.df2, 0u);
    const Real gc_hat = std::log(1.0 + gres.f_stat / static_cast<Real>(gres.df2));
    const Real ratio = te / gc_hat;
    EXPECT_NEAR(ratio, 0.5, 0.15) << "TE/GC_hat=" << ratio << " (expect 0.35..0.65)";
}

TEST(Phase9Integration, KdeBandwidthSensitivity) {
    const VarSample v = generate_var(800, 0.5, 0.3, 920);
    de::KdeConfig scott = scott_config();
    de::KdeConfig silver;
    silver.bandwidth.rule = de::BandwidthRule::kSilverman;

    const de::KernelDensity1D kde_scott(std::span<const Real>(v.x), scott);
    const de::KernelDensity1D kde_silver(std::span<const Real>(v.x), silver);
    EXPECT_GT(kde_scott.bandwidth(), kde_silver.bandwidth());
    EXPECT_NEAR(integrate_pdf(kde_scott, -6.0, 6.0, 1200), 1.0, 5e-3);
    EXPECT_NEAR(integrate_pdf(kde_silver, -6.0, 6.0, 1200), 1.0, 5e-3);

    const std::vector<Real> bandwidths = {kde_scott.bandwidth(), kde_silver.bandwidth()};
    const std::vector<Real> grid = {-2.0, -1.0, 0.0, 1.0, 2.0};
    const auto pts = de::bandwidth_sensitivity(
        std::span<const Real>(v.x), std::span<const Real>(bandwidths), std::span<const Real>(grid));
    ASSERT_EQ(pts.size(), 2u);
    for (const auto& p : pts) {
        EXPECT_TRUE(std::isfinite(p.log_likelihood));
        EXPECT_TRUE(std::isfinite(p.pdf_at_mode));
    }
}
