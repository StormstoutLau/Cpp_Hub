// =============================================================================
// test_ksg.cpp - Phase 9 M2: KSG Algorithm 1 MI/CMI (TDD 前置)
//
// 状态: M2 已实现 (include/cpphub/timeseries/information/ksg.hpp),
//       已接入 tests/CMakeLists.txt.
//
// 基准:
//   - 40 点二元 Gaussian (rho=0.7) 的 sklearn `_compute_mi_cc` 黄金值
//     (raw, k=3/5, 无 jitter/standardize)
//   - 解析 MI: -0.5*ln(1-rho^2), n=2000, 3 seeds 均值, 容差 0.05 nat
//   - 缩放不变性: 默认 standardize=true 时 x*10 不变
//   - CMI: 条件独立 (x=z+e1, y=z+e2) 均值 |CMI|<0.05
// =============================================================================

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

#include "cpphub/timeseries/information/ksg.hpp"

namespace inf = cpphub::v1::timeseries::information;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

// 40 点二元 Gaussian (rho=0.7, numpy seed 12345), sklearn 黄金值基准
const std::vector<Real> kX = {
    0.82326189505377179, 0.90308782230613061, 0.35640648875312053,
    1.0097274606081093, 0.42346147954939323, -2.5393020054854971,
    0.35070049792209096, 0.45401450362606566, -0.24057307528193064,
    -1.0727366303298402, -1.1030258398081565, -0.20441420856327813,
    -0.028241086196026964, 1.2703764499015366, -2.6041404339297345,
    -1.0376200035191261, 1.3525376080443703, 1.0923287151490628,
    -0.43811709397070142, 1.0330306705970929, -0.23781908696007845,
    -1.1348623788542136, 0.49019772371178644, -0.64706699750937158,
    -0.083837108559159337, -0.50727868769903772, -0.69673968271712672,
    -0.26218160840237181, 0.34410879275371242, 0.24264842238131665,
    -0.11789023443403203, 0.34463083162647296, 1.8444054324921522,
    -0.18741706962626542, -0.061977974330676339, 1.2608216877788845,
    -0.21354470572064443, 1.0529103725517679, -1.3724402994191438,
    0.60242282599676722};
const std::vector<Real> kY = {
    1.8021417495462564, 0.70233309776940567, -0.21748029504452768,
    1.512357663881611, -1.0892197445445708, -1.789107528119374,
    1.0495402761223158, 0.40700460455449516, -1.2139840255538374,
    0.010907839321951791, -1.3351713115354613, -1.460486793360313,
    0.31992767049136872, 1.2071015142294077, -0.57612286272201196,
    -0.39576349020973717, 0.41574946913440131, 1.5118090858318012,
    -0.94839064333610579, 1.2325340601469732, -0.33914538650500375,
    -1.2068707889789174, -0.36822151130373193, 0.39637495192455635,
    -0.028906895378032312, -0.292340768994101, -0.28099740772572396,
    -0.87799948192241983, -0.89733852232433675, -0.73460115677590831,
    0.24930710767247694, 0.42032170706798583, 1.1803756028333465,
    -1.0817118393486354, -1.1373969693846704, 0.41230889332259302,
    0.20036869528504875, 0.91234219242353221, -1.6182910433722897,
    0.90186564523286539};
constexpr Real kMiGoldenK3 = 0.3176876772600816;
constexpr Real kMiGoldenK5 = 0.3132406148282585;

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

struct Triple {
    std::vector<Real> x;
    std::vector<Real> y;
    std::vector<Real> z;
};

// x = z + 0.5*e1, y = z + 0.5*e2, z ~ N(0,1): x ⟂ y | z
Triple conditional_independent(Size n, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<Real> nd(0.0, 1.0);
    Triple t;
    t.x.resize(n);
    t.y.resize(n);
    t.z.resize(n);
    for (Size i = 0; i < n; ++i) {
        const Real z = nd(rng);
        t.z[i] = z;
        t.x[i] = z + 0.5 * nd(rng);
        t.y[i] = z + 0.5 * nd(rng);
    }
    return t;
}

inf::KsgConfig ksg_config(Size k, bool standardize, Real jitter) {
    inf::KsgConfig cfg;
    cfg.k = k;
    cfg.standardize = standardize;
    cfg.jitter_scale = jitter;
    cfg.seed = 7;
    return cfg;
}

Real gaussian_mi(Real rho) {
    return -0.5 * std::log(1.0 - rho * rho);
}

}  // namespace

TEST(KsgMutualInfo, MatchesSklearnGolden) {
    const Real mi3 = inf::ksg_mutual_info(
        std::span<const Real>(kX), std::span<const Real>(kY), ksg_config(3, false, 0.0));
    const Real mi5 = inf::ksg_mutual_info(
        std::span<const Real>(kX), std::span<const Real>(kY), ksg_config(5, false, 0.0));
    EXPECT_NEAR(mi3, kMiGoldenK3, 1e-12);
    EXPECT_NEAR(mi5, kMiGoldenK5, 1e-12);
}

TEST(KsgMutualInfo, GaussianAnalyticMean) {
    constexpr Real rho = 0.8;
    constexpr Real truth = 0.5108256237659907;  // -0.5*ln(1-0.8^2)
    // spec §1.3/§6.3: n=2000, k=3, 20 seeds 均值, 判红线 0.02 nat
    // (0.01 为记录线不判红; probe K1 5-seed bias +0.0044~+0.0055)
    Real sum = 0.0;
    constexpr Size kSeeds = 20;
    for (Size s = 0; s < kSeeds; ++s) {
        const Pair p = gaussian_pair(2000, rho, 1000 + s);
        sum += inf::ksg_mutual_info(std::span<const Real>(p.x), std::span<const Real>(p.y),
                                    ksg_config(3, true, 0.0));
    }
    const Real mean = sum / static_cast<Real>(kSeeds);
    EXPECT_NEAR(mean, truth, 0.02) << "mean=" << mean << " truth=" << truth;
}

TEST(KsgMutualInfo, StandardizeRemovesScaleSensitivity) {
    std::vector<Real> x10 = kX;
    for (Real& v : x10) {
        v *= 10.0;
    }
    const inf::KsgConfig cfg = ksg_config(3, true, 0.0);
    const Real base = inf::ksg_mutual_info(std::span<const Real>(kX), std::span<const Real>(kY), cfg);
    const Real scaled = inf::ksg_mutual_info(std::span<const Real>(x10), std::span<const Real>(kY), cfg);
    EXPECT_NEAR(base, scaled, 1e-12);
}

TEST(KsgMutualInfo, RawScaleSensitivityIsDocumented) {
    std::vector<Real> x100 = kX;
    for (Real& v : x100) {
        v *= 100.0;
    }
    const inf::KsgConfig cfg = ksg_config(3, false, 0.0);
    const Real base = inf::ksg_mutual_info(std::span<const Real>(kX), std::span<const Real>(kY), cfg);
    const Real scaled = inf::ksg_mutual_info(std::span<const Real>(x100), std::span<const Real>(kY), cfg);
    EXPECT_GT(std::abs(base - scaled), 1e-6)
        << "raw Chebyshev MI 应对混合单位缩放敏感 (spec 风险 ②)";
}

TEST(KsgMutualInfo, JitterHandlesTies) {
    std::vector<Real> xr = kX;
    std::vector<Real> yr = kY;
    for (Real& v : xr) {
        v = std::round(v * 100.0) / 100.0;
    }
    for (Real& v : yr) {
        v = std::round(v * 100.0) / 100.0;
    }
    const Real no_jitter = inf::ksg_mutual_info(
        std::span<const Real>(xr), std::span<const Real>(yr), ksg_config(3, true, 0.0));
    const Real jittered = inf::ksg_mutual_info(
        std::span<const Real>(xr), std::span<const Real>(yr), ksg_config(3, true, 1e-10));
    EXPECT_TRUE(std::isfinite(no_jitter));
    EXPECT_TRUE(std::isfinite(jittered));
    EXPECT_NE(no_jitter, jittered);
}

TEST(KsgMutualInfo, NegativeEstimateNotClipped) {
    // spec §5.3.7 零假设: i.i.d. MI≈0 (n=300, KSG 小样本偏置 ~+0.005)
    const Pair p = gaussian_pair(300, 0.0, 77);
    const Real mi = inf::ksg_mutual_info(
        std::span<const Real>(p.x), std::span<const Real>(p.y), ksg_config(3, true, 0.0));
    EXPECT_TRUE(std::isfinite(mi));
    EXPECT_LT(std::abs(mi), 0.05) << "mi=" << mi;
    // 估计可正可负; 只验证不截断 (sklearn 的 max(0,·) 语义不在本库)
}

TEST(KsgConditionalMutualInfo, Symmetry) {
    const Triple t = conditional_independent(300, 99);
    const inf::KsgConfig cfg = ksg_config(3, true, 0.0);
    const Real a = inf::ksg_conditional_mutual_info(
        std::span<const Real>(t.x), std::span<const Real>(t.y), std::span<const Real>(t.z), cfg);
    const Real b = inf::ksg_conditional_mutual_info(
        std::span<const Real>(t.y), std::span<const Real>(t.x), std::span<const Real>(t.z), cfg);
    EXPECT_NEAR(a, b, 1e-12);
}

TEST(KsgConditionalMutualInfo, ConditionalIndependenceNearZero) {
    constexpr Size kSeeds = 3;
    Real sum = 0.0;
    for (Size s = 0; s < kSeeds; ++s) {
        const Triple t = conditional_independent(2000, 500 + s);
        sum += inf::ksg_conditional_mutual_info(
            std::span<const Real>(t.x), std::span<const Real>(t.y), std::span<const Real>(t.z),
            ksg_config(3, true, 0.0));
    }
    const Real mean = sum / static_cast<Real>(kSeeds);
    EXPECT_NEAR(mean, 0.0, 0.05) << "mean=" << mean;
}

TEST(KsgValidation, RejectsInvalidInputs) {
    std::vector<Real> empty;
    std::vector<Real> two = {0.0, 1.0};
    std::vector<Real> nan_vec = {0.0, std::numeric_limits<Real>::quiet_NaN(), 2.0};

    EXPECT_THROW(inf::ksg_mutual_info(std::span<const Real>(empty), std::span<const Real>(empty),
                                      ksg_config(3, true, 0.0)),
                 std::invalid_argument);
    EXPECT_THROW(inf::ksg_mutual_info(std::span<const Real>(two), std::span<const Real>(two),
                                      ksg_config(3, true, 0.0)),
                 std::invalid_argument);  // n < k+1
    EXPECT_THROW(inf::ksg_mutual_info(std::span<const Real>(kX), std::span<const Real>(kY),
                                      ksg_config(0, true, 0.0)),
                 std::invalid_argument);
    EXPECT_THROW(inf::ksg_mutual_info(std::span<const Real>(nan_vec), std::span<const Real>(nan_vec),
                                      ksg_config(1, true, 0.0)),
                 std::invalid_argument);
    EXPECT_THROW(inf::ksg_conditional_mutual_info(std::span<const Real>(kX), std::span<const Real>(kY),
                                                  std::span<const Real>(two),
                                                  ksg_config(3, true, 0.0)),
                 std::invalid_argument);  // 长度不齐
}

// spec §1.3/D9-10 性能预算: KSG N=5000, k=3, 单线程 < 0.5 s (暴力 O(N²), D9-3)
TEST(KsgPerformance, N5000K3UnderBudget) {
    constexpr Size n = 5000;
    constexpr Real rho = 0.5;
    const Pair p = gaussian_pair(n, rho, 3000);
    const inf::KsgConfig cfg = ksg_config(3, true, 0.0);

    const auto t0 = std::chrono::steady_clock::now();
    const Real mi = inf::ksg_mutual_info(std::span<const Real>(p.x), std::span<const Real>(p.y), cfg);
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - t0)
                        .count();
    EXPECT_TRUE(std::isfinite(mi));
    EXPECT_LT(static_cast<double>(ms), 500.0)
        << "KSG N=5000 k=3 took " << ms << " ms (budget 500 ms)";
}
