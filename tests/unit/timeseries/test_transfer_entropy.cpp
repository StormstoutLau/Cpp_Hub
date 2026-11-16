// =============================================================================
// test_transfer_entropy.cpp - Phase 9 M3: KSG-CMI 迁移熵 (TDD 前置)
//
// 状态: M3 已实现 (include/cpphub/timeseries/information/transfer_entropy.hpp),
//       已接入 tests/CMakeLists.txt.
//
// 锚:
//   - Gaussian VAR(phi=0.5, a=0.3): GC = ln(1+a^2/(1-phi^2)) = 0.113329,
//     TE = GC/2 = 0.056664 (nats); KSG 估计容差 ±35% (probe K4: +25.0%)
//   - te_bits = te / ln2
//   - 默认 standardize=true 时 x*10 不变
//   - 方向性: TE(x->y) > TE(y->x) (a>0)
// =============================================================================

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

#include "cpphub/core/rng.hpp"
#include "cpphub/timeseries/information/transfer_entropy.hpp"

namespace inf = cpphub::v1::timeseries::information;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

struct VarSample {
    std::vector<Real> x;
    std::vector<Real> y;
};

// x_t = phi*x_{t-1} + e_x; y_t = phi*y_{t-1} + a*x_{t-1} + e_y
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

inf::TransferEntropyConfig te_config(Size n_shuffles, Size n_block_boot) {
    inf::TransferEntropyConfig cfg;
    cfg.source_history = 1;
    cfg.target_history = 1;
    cfg.delay = 1;
    cfg.ksg_neighbors = 3;
    cfg.standardize = true;
    cfg.n_shuffles = n_shuffles;
    cfg.n_block_boot = n_block_boot;
    cfg.block_length = 0;
    cfg.seed = 42;
    return cfg;
}

constexpr Real kGc = 0.11332868530700327;  // ln(1 + 0.3^2/(1-0.5^2))
constexpr Real kTe = 0.056664342653501635; // GC/2

// H2/风险⑭ 关闭锚: JIDT 1.6.1 TransferEntropyCalculatorKraskov 运行时对照黄金值
// (probe_jidt_te.py 实测; NOISE_LEVEL_TO_ADD=0, DYN_CORR_EXCL=0; 数据 = Philox
// VAR 与本测试 jidt_anchor_var() 位精确同源)。延迟=1 结构区 12 组 max 偏差
// 6.53e-4 nat (组装层=嵌入/差分/有效样本裁剪 与 JIDT 同构), 锚容差放宽 3×。
constexpr Real kJidtBaseSeed3000 = 0.072790;  // (k,l,d,K)=(1,1,1,3)
constexpr Real kJidtDeepSeed3000 = 0.044291;  // (k,l,d,K)=(3,2,1,3)
constexpr Real kJidtAnchorTol = 0.002;

// Philox4x64 + normal_simd: 跨平台位精确的 VAR(1) (与驱动/探针同源, 避开
// std::normal_distribution 的标准库实现差异)
VarSample jidt_anchor_var(Size n, Real phi, Real a, std::uint64_t seed) {
    cpphub::v1::Philox4x64 rng(seed);
    VarSample s;
    s.x.resize(n);
    s.y.resize(n);
    // f64x4 AVX2 分支无 operator[] — 跨平台取元素统一走 store (test_rng.cpp 同惯例)
    alignas(32) double z0v[4];
    cpphub::v1::store(z0v, cpphub::v1::normal_simd(rng));
    s.x[0] = z0v[0] / std::sqrt(1.0 - phi * phi);
    s.y[0] = z0v[1] / std::sqrt(1.0 - phi * phi);
    for (Size t = 1; t < n; ++t) {
        alignas(32) double zv[4];
        cpphub::v1::store(zv, cpphub::v1::normal_simd(rng));
        s.x[t] = phi * s.x[t - 1] + zv[0];
        s.y[t] = phi * s.y[t - 1] + a * s.x[t - 1] + zv[1];
    }
    return s;
}

// JIDT 对照口径的纯点估计配置 (无置换, 本库默认标准化对齐探针 z-score)
inf::TransferEntropyConfig jidt_anchor_config(Size k, Size l, Size delay, Size knn) {
    inf::TransferEntropyConfig cfg;
    cfg.source_history = k;
    cfg.target_history = l;
    cfg.delay = delay;
    cfg.ksg_neighbors = knn;
    cfg.standardize = true;
    cfg.n_shuffles = 0;
    cfg.n_block_boot = 0;
    return cfg;
}

}  // namespace

TEST(TransferEntropy, GaussianVarMatchesGcHalf) {
    // spec §5.3.1/§6.3: 5 seeds 均值 ±30% (probe K4 bias +25%, 固定 seed 确定性)
    constexpr Size kSeeds = 5;
    Real sum = 0.0;
    for (Size s = 0; s < kSeeds; ++s) {
        const VarSample v = generate_var(2000, 0.5, 0.3, 200 + s);
        sum += inf::transfer_entropy(std::span<const Real>(v.x), std::span<const Real>(v.y),
                                     te_config(10, 10))
                   .te;
    }
    const Real mean = sum / static_cast<Real>(kSeeds);
    EXPECT_NEAR(mean, kTe, 0.30 * kTe) << "mean=" << mean << " GC/2=" << kTe;
}

TEST(TransferEntropy, BitsConversion) {
    const VarSample v = generate_var(500, 0.5, 0.3, 201);
    const auto r = inf::transfer_entropy(std::span<const Real>(v.x), std::span<const Real>(v.y),
                                         te_config(5, 5));
    EXPECT_NEAR(r.te_bits, r.te / std::log(2.0), 1e-15);
}

TEST(TransferEntropy, Directionality) {
    const VarSample v = generate_var(1000, 0.5, 0.3, 202);
    const auto fwd = inf::transfer_entropy(std::span<const Real>(v.x), std::span<const Real>(v.y),
                                           te_config(5, 5));
    const auto rev = inf::transfer_entropy(std::span<const Real>(v.y), std::span<const Real>(v.x),
                                           te_config(5, 5));
    EXPECT_GT(fwd.te, rev.te + 0.005) << "fwd=" << fwd.te << " rev=" << rev.te;
}

TEST(TransferEntropy, StandardizeRemovesScaleSensitivity) {
    VarSample v = generate_var(600, 0.5, 0.3, 203);
    std::vector<Real> x10 = v.x;
    for (Real& value : x10) {
        value *= 10.0;
    }
    const auto base = inf::transfer_entropy(std::span<const Real>(v.x), std::span<const Real>(v.y),
                                            te_config(5, 5));
    const auto scaled = inf::transfer_entropy(std::span<const Real>(x10), std::span<const Real>(v.y),
                                              te_config(5, 5));
    EXPECT_NEAR(base.te, scaled.te, 1e-10);
}

TEST(TransferEntropy, ResultFields) {
    constexpr Size n = 500;
    const VarSample v = generate_var(n, 0.5, 0.3, 204);
    const auto r = inf::transfer_entropy(std::span<const Real>(v.x), std::span<const Real>(v.y),
                                         te_config(5, 5));
    EXPECT_EQ(r.n_effective, n - 1);
    EXPECT_TRUE(std::isfinite(r.te));
    EXPECT_TRUE(std::isfinite(r.te_effective));
    EXPECT_TRUE(std::isfinite(r.te_bits));
    EXPECT_GE(r.p_value, 0.0);
    EXPECT_LE(r.p_value, 1.0);
}

TEST(TransferEntropy, IidPValueInRange) {
    std::mt19937_64 rng(205);
    std::normal_distribution<Real> nd(0.0, 1.0);
    std::vector<Real> x(300), y(300);
    for (Size i = 0; i < x.size(); ++i) {
        x[i] = nd(rng);
        y[i] = nd(rng);
    }
    const auto r = inf::transfer_entropy(std::span<const Real>(x), std::span<const Real>(y),
                                         te_config(10, 10));
    EXPECT_GT(r.p_value, 0.0);
    EXPECT_LE(r.p_value, 1.0);
}

// spec §5.3.5a: 20 条 i.i.d. 序列, p>0.05 的比例 ≥ 60% (块置换 p 值零假设校准)
TEST(TransferEntropy, IidPValueCalibration) {
    constexpr Size kSeries = 20;
    constexpr Size n = 300;
    Size not_rejected = 0;
    for (Size s = 0; s < kSeries; ++s) {
        std::mt19937_64 rng(400 + s);
        std::normal_distribution<Real> nd(0.0, 1.0);
        std::vector<Real> x(n), y(n);
        for (Size i = 0; i < n; ++i) {
            x[i] = nd(rng);
            y[i] = nd(rng);
        }
        inf::TransferEntropyConfig cfg = te_config(10, 50);
        cfg.seed = 400 + s;
        const auto r = inf::transfer_entropy(std::span<const Real>(x),
                                             std::span<const Real>(y), cfg);
        if (r.p_value > 0.05) {
            ++not_rejected;
        }
    }
    EXPECT_GE(not_rejected * 100 / kSeries, 60u)
        << "p>0.05 proportion = " << not_rejected << "/" << kSeries;
}

// spec §5.3.5b: 耦合 VAR 5 seeds, p < 0.1 (真因果应被检出)
// n=4000 (spec 未冻结 n): n=2000 实测 seed 502 p=0.157 — TE 估计方差
// (SD~0.02 nat) 与块置换零分布重叠; 加倍 n 收窄方差提高功效。
// n_shuffles=0 跳过 effective TE 评估 (本用例不断言 te_effective);
// B=50 保 p 分辨率 1/51。运行成本 ~30s (5 seeds × 51 次 CMI @ n=4000)
TEST(TransferEntropy, CoupledVarPValueSignificant) {
    for (Size s = 0; s < 5; ++s) {
        const VarSample v = generate_var(4000, 0.5, 0.3, 500 + s);
        inf::TransferEntropyConfig cfg = te_config(0, 50);
        cfg.seed = 500 + s;
        const auto r = inf::transfer_entropy(std::span<const Real>(v.x),
                                             std::span<const Real>(v.y), cfg);
        EXPECT_LT(r.p_value, 0.1) << "seed " << 500 + s << " p=" << r.p_value;
    }
}

TEST(TransferEntropy, MultivariateEmbedding) {
    // spec §5.3.6: k=3/l=2/delay=2 通用嵌入路径 (非 k=l=1 特例)
    constexpr Size n = 400;
    const VarSample v = generate_var(n, 0.5, 0.3, 207);
    inf::TransferEntropyConfig cfg = te_config(3, 3);
    cfg.source_history = 3;
    cfg.target_history = 2;
    cfg.delay = 2;
    const auto r = inf::transfer_entropy(std::span<const Real>(v.x),
                                         std::span<const Real>(v.y), cfg);
    EXPECT_EQ(r.n_effective, n - 3 * 2);
    EXPECT_TRUE(std::isfinite(r.te));
    EXPECT_TRUE(std::isfinite(r.te_effective));
    EXPECT_NEAR(r.te_bits, r.te / std::log(2.0), 1e-15);
    EXPECT_GE(r.p_value, 0.0);
    EXPECT_LE(r.p_value, 1.0);
}

TEST(TransferEntropy, EmbeddingBoundaryThrows) {
    const VarSample v = generate_var(8, 0.5, 0.3, 206);
    inf::TransferEntropyConfig cfg = te_config(2, 2);
    cfg.source_history = 4;
    cfg.target_history = 4;
    cfg.delay = 2;
    EXPECT_THROW(inf::transfer_entropy(std::span<const Real>(v.x), std::span<const Real>(v.y), cfg),
                 std::invalid_argument);
}

TEST(TransferEntropy, Validation) {
    std::vector<Real> empty;
    std::vector<Real> short_series = {0.0, 1.0, 2.0};
    std::vector<Real> mismatch = {0.0, 1.0, 2.0, 3.0};

    EXPECT_THROW(inf::transfer_entropy(std::span<const Real>(empty), std::span<const Real>(empty),
                                       te_config(2, 2)),
                 std::invalid_argument);
    EXPECT_THROW(inf::transfer_entropy(std::span<const Real>(short_series), std::span<const Real>(short_series),
                                       te_config(2, 2)),
                 std::invalid_argument);
    EXPECT_THROW(inf::transfer_entropy(std::span<const Real>(short_series), std::span<const Real>(mismatch),
                                       te_config(2, 2)),
                 std::invalid_argument);
}

// spec D9-10 性能预算: TE N=2000, ETE 100 + 块置换 100 (共 ~200 次 CMI 评估)
// < 10 s (MSVC Release, 单线程; 审计 C6: 预算含 p 值置换)
// 覆盖率 CI 为 Debug -O0 + --coverage 插桩构建, 预算口径不适用 → skip
TEST(TransferEntropy, PerformanceBudgetN2000) {
#if defined(CPPHUB_COVERAGE_BUILD)
    GTEST_SKIP() << "perf budget (PHASE9_SPEC D9-10) asserted on optimized "
                    "builds; coverage job is Debug -O0 + --coverage";
#endif
    constexpr Size n = 2000;
    const VarSample v = generate_var(n, 0.5, 0.3, 3000);
    const auto cfg = te_config(100, 100);

    const auto t0 = std::chrono::steady_clock::now();
    const auto r = inf::transfer_entropy(std::span<const Real>(v.x),
                                         std::span<const Real>(v.y), cfg);
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - t0)
                        .count();
    EXPECT_TRUE(std::isfinite(r.te));
    EXPECT_GE(r.p_value, 0.0);
    EXPECT_LE(r.p_value, 1.0);
    EXPECT_LT(static_cast<double>(ms), 10000.0)
        << "TE N=2000 with ETE+block boot took " << ms << " ms (budget 10 s)";
}

// H2/风险⑭ 关闭锚 (PHASE9_SPEC §8/§10.1): TE 组装层第二实现对照
// JIDT 1.6.1 TransferEntropyCalculatorKraskov 运行时对照 (probe_jidt_te.py):
// 嵌入映射 (JIDT kHist=本库 target_history, lHist=source_history, delay=1 同构),
// 单位 (JIDT 连续 TE = nats, 订正调研 A47 "bits 默认"), 有效样本裁剪一致。
// 数据跨平台位精确 (Philox), 黄金值为 JIDT 实测值; 容差 0.002 覆盖两实现的
// KSG 内部差异 (实测 max 6.53e-4) 与 GCC/MSVC 浮点次序差。
TEST(TransferEntropy, JidtAssemblyAnchor) {
    const VarSample v = jidt_anchor_var(2000, 0.5, 0.3, 3000);
    const auto base = inf::transfer_entropy(std::span<const Real>(v.x),
                                            std::span<const Real>(v.y),
                                            jidt_anchor_config(1, 1, 1, 3));
    EXPECT_NEAR(base.te, kJidtBaseSeed3000, kJidtAnchorTol)
        << "base (k,l,d,K)=(1,1,1,3): te=" << base.te << " vs JIDT " << kJidtBaseSeed3000;
    const auto deep = inf::transfer_entropy(std::span<const Real>(v.x),
                                            std::span<const Real>(v.y),
                                            jidt_anchor_config(3, 2, 1, 3));
    EXPECT_NEAR(deep.te, kJidtDeepSeed3000, kJidtAnchorTol)
        << "deep (k,l,d,K)=(3,2,1,3): te=" << deep.te << " vs JIDT " << kJidtDeepSeed3000;
    // 组装层语义不变式: 双边历史加深后估计仍有限、bits 换算仍成立
    EXPECT_TRUE(std::isfinite(deep.te));
    EXPECT_NEAR(deep.te_bits, deep.te / std::log(2.0), 1e-15);
}
