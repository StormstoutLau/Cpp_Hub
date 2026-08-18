// =============================================================================
// test_midas_weights.cpp - M4 权重族测试 (16 用例, spec §1.2 测试矩阵)
//
// 基准: tests/unit/timeseries/midas_baseline.inc (midasr 0.9 源函数直调,
//       verify_midas.R W1, 2026-08-18) — 逐点 1e-12 (同公式复刻, 无优化器)
//
// 覆盖幻觉点 (spec §6.2/MD1/MD2/MD7):
//   MD1: nealmon i 从 1 起 / nbeta xi 从 0 起 (两套起点)
//   MD2: Σw = δ 独立尺度
//   MD7: log-sum-exp 防溢出, 非溢出区间与裸公式差 <1e-14
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/midas/midas_weights.hpp"
#include "midas_baseline.inc"

namespace mw = cpphub::v1::timeseries::midas;
namespace bl = cpphub::v1::timeseries::midas_baseline::v1;
using cpphub::Real;
using cpphub::Size;

// 1. nealmon 2 参数逐点 vs midasr (1e-12, MD1 i 从 1 起)
TEST(MidasWeights, Nealmon2ParamMatchesR) {
    const auto w = mw::nealmon_weights({1.0, -0.5}, 4);
    ASSERT_EQ(w.size(), 4u);
    for (Size i = 0; i < 4; ++i) EXPECT_NEAR(w[i], bl::NEALMON1[i], 1e-12);
}

// 2. nealmon 3 参数逐点 (二次形状, 驼峰)
TEST(MidasWeights, Nealmon3ParamMatchesR) {
    const auto w = mw::nealmon_weights({0.8, 0.1, -0.02}, 6);
    ASSERT_EQ(w.size(), 6u);
    for (Size i = 0; i < 6; ++i) EXPECT_NEAR(w[i], bl::NEALMON2[i], 1e-12);
}

// 3. MD2: Σw = δ (首参数独立尺度, 非归一化 1)
TEST(MidasWeights, NealmonSumEqualsDelta) {
    const Real deltas[] = {1.0, 2.5, -0.7, 10.0};
    for (Real d : deltas) {
        const auto w = mw::nealmon_weights({d, -0.3}, 5);
        Real s = 0.0;
        for (Real v : w) s += v;
        EXPECT_NEAR(s, d, 1e-12);
    }
}

// 4. MD7: log-sum-exp 与裸公式差 < 1e-14 (非溢出区间)
TEST(MidasWeights, NealmonLogSumExpMatchesNaive) {
    const std::vector<Real> lam{2.0, -0.4, 0.02};
    const Size d = 8;
    const auto w = mw::nealmon_weights(lam, d);
    // 裸公式: δ·exp(η)/Σexp(η)
    std::vector<Real> eta(d);
    for (Size i = 1; i <= d; ++i) {
        eta[i - 1] = lam[1] * static_cast<Real>(i)
                     + lam[2] * static_cast<Real>(i) * static_cast<Real>(i);
    }
    Real se = 0.0;
    for (Real e : eta) se += std::exp(e);
    for (Size i = 0; i < d; ++i) {
        EXPECT_NEAR(w[i], lam[0] * std::exp(eta[i]) / se, 1e-14);
    }
}

// 5. MD7: 溢出防护 (裸公式 exp(λ₂·d) 溢出区间, lse 仍有限)
TEST(MidasWeights, NealmonOverflowProtected) {
    // λ₂ = −0.2, d = 200: η_min = −40 (裸 exp 下溢 0 → 0/0 NaN 风险)
    const auto w = mw::nealmon_weights({1.0, -0.2}, 200);
    Real s = 0.0;
    for (Real v : w) {
        EXPECT_TRUE(std::isfinite(v));
        s += v;
    }
    EXPECT_NEAR(s, 1.0, 1e-10);
}

// 6. nbetaMT θ₀=0 逐点 (xi 从 0 起, 端点 eps 防护)
TEST(MidasWeights, NbetaThetaZeroMatchesR) {
    const auto w = mw::nbeta_weights({1.0, 2.0, 3.0, 0.0}, 5);
    ASSERT_EQ(w.size(), 5u);
    for (Size i = 0; i < 5; ++i) {
        EXPECT_NEAR(w[i], bl::NBETA1[i], 1e-12);
    }
}

// 7. nbetaMT θ₀=0.1 逐点 (均匀混合)
TEST(MidasWeights, NbetaThetaMixMatchesR) {
    const auto w = mw::nbeta_weights({0.5, 1.5, 1.5, 0.1}, 7);
    ASSERT_EQ(w.size(), 7u);
    for (Size i = 0; i < 7; ++i) {
        EXPECT_NEAR(w[i], bl::NBETA2[i], 1e-12);
    }
}

// 8. nbeta Σw = δ
TEST(MidasWeights, NbetaSumEqualsDelta) {
    const auto w = mw::nbeta_weights({3.0, 1.5, 2.5, 0.1}, 9);
    Real s = 0.0;
    for (Real v : w) s += v;
    EXPECT_NEAR(s, 3.0, 1e-12);
}

// 9. nbeta d=1 拒绝 (midasr NaN 同)
TEST(MidasWeights, NbetaRejectsD1) {
    EXPECT_THROW(mw::nbeta_weights({1.0, 1.5, 1.5, 0.0}, 1),
                 std::invalid_argument);
}

// 10. nbeta 参数个数校验 (4 参数)
TEST(MidasWeights, NbetaRejectsWrongParamCount) {
    EXPECT_THROW(mw::nbeta_weights({1.0, 1.5, 1.5}, 5),
                 std::invalid_argument);
    EXPECT_THROW(mw::nbeta_weights({1.0, 1.5, 1.5, 0.0, 0.2}, 5),
                 std::invalid_argument);
}

// 11. almonp 逐点 (不归一化, 可负)
TEST(MidasWeights, AlmonpMatchesR) {
    const auto w = mw::almonp_weights({0.1, 0.05, -0.01}, 5);
    ASSERT_EQ(w.size(), 5u);
    for (Size i = 0; i < 5; ++i) {
        EXPECT_NEAR(w[i], bl::ALMONP[i], 1e-12);
    }
}

// 12. almonp 手算: w_i = c + a·i (raw poly)
TEST(MidasWeights, AlmonpHandComputed) {
    const auto w = mw::almonp_weights({2.0, -0.5}, 4);
    // w_i = 2 − 0.5·i, i = 1..4 → 1.5, 1.0, 0.5, 0.0
    EXPECT_NEAR(w[0], 1.5, 1e-15);
    EXPECT_NEAR(w[1], 1.0, 1e-15);
    EXPECT_NEAR(w[2], 0.5, 1e-15);
    EXPECT_NEAR(w[3], 0.0, 1e-15);
    // Σw ≠ 1 (不归一化)
    Real s = 0.0;
    for (Real v : w) s += v;
    EXPECT_NE(s, 1.0);
}

// 13. polystep 逐点 (断点 {2,5}, 段值 3 段)
TEST(MidasWeights, PolystepMatchesR) {
    const auto w = mw::polystep_weights({0.5, 0.2, 0.1}, {2, 5}, 8);
    ASSERT_EQ(w.size(), 8u);
    for (Size i = 0; i < 8; ++i) {
        EXPECT_NEAR(w[i], bl::POLYSTEP[i], 1e-12);
    }
}

// 14. polystep 断点越界/非严格递增拒绝
TEST(MidasWeights, PolystepRejectsBadSteps) {
    EXPECT_THROW(mw::polystep_weights({0.5, 0.2}, {0, 5}, 8),
                 std::invalid_argument);  // 断点 ≤ 0
    EXPECT_THROW(mw::polystep_weights({0.5, 0.2}, {5, 8}, 8),
                 std::invalid_argument);  // 断点 ≥ d
    EXPECT_THROW(mw::polystep_weights({0.5, 0.2, 0.1}, {5, 2}, 8),
                 std::invalid_argument);  // 非递增
    EXPECT_THROW(mw::polystep_weights({0.5}, {2}, 8),
                 std::invalid_argument);  // 段数不匹配
}

// 15. harstep 逐点 (d=20, HAR(3) 结构)
TEST(MidasWeights, HarstepMatchesR) {
    const auto w = mw::harstep_weights({0.6, 0.3, 0.1}, 20);
    ASSERT_EQ(w.size(), 20u);
    EXPECT_NEAR(w[0], bl::HARSTEP[0], 1e-12);  // w1 = p1+p2/5+p3/20
    EXPECT_NEAR(w[1], bl::HARSTEP[1], 1e-12);  // w2 = p2/5+p3/20
    EXPECT_NEAR(w[5], bl::HARSTEP[2], 1e-12);  // w6 = p3/20
    EXPECT_NEAR(w[19], bl::HARSTEP[3], 1e-12); // w20 = p3/20
}

// 16. harstep d≠20 拒绝 + nealmon 空/单参数
TEST(MidasWeights, HarstepAndNealmonEdgeCases) {
    EXPECT_THROW(mw::harstep_weights({0.6, 0.3, 0.1}, 19),
                 std::invalid_argument);
    EXPECT_THROW(mw::harstep_weights({0.6, 0.3}, 20), std::invalid_argument);
    // nealmon 仅 δ: 均匀 δ/d
    const auto w = mw::nealmon_weights({2.0}, 4);
    for (Real v : w) EXPECT_NEAR(v, 0.5, 1e-15);
    EXPECT_THROW(mw::nealmon_weights({}, 4), std::invalid_argument);
    EXPECT_THROW(mw::nealmon_weights({1.0, -0.5}, 0), std::invalid_argument);
}
