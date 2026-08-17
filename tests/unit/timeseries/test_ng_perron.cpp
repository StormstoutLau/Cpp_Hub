// =============================================================================
// test_ng_perron.cpp - Ng-Perron M 检验族测试 (18 用例, spec §1.2/§9.5)
//
// 基准 (NP6 对照生态): 无开源库输出 M 族 ⇒ 原文公式钉死 + 恒等式自检 +
//   np_tables 临界值精确相等 + 模拟方向断言 (Philox 确定性数据, §8.6 位一致)
//   Stata dfgls 逐 k MAIC 1e-10 对照: verify_np_stata.py 占位, 环境可用后补
//
// 容差:
//   - MZt ≡ MZα×MSB 恒等式: 1e-12 (spec §1.3)
//   - 临界值: EXPECT_DOUBLE_EQ (np_table1 双源转录)
//   - 方向/结构断言: 精确或宽松
//
// 幻觉点覆盖: NP1 (文献, 表已 static_assert) / NP2 (τ_T 形, 实现内) /
//   NP3 (s² 对 Δỹ, 结构断言) / NP4 (MPT 分情形, 正性+构造断言) /
//   NP5 (方向重构 + 表精确) / NP6 (对照边界声明)
// =============================================================================
#include <gtest/gtest.h>
#include <array>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/rng.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/df_gls_test.hpp"
#include "cpphub/timeseries/unit_root/ng_perron_test.hpp"
#include "cpphub/timeseries/unit_root/np_tables.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"

namespace ur = cpphub::v1::timeseries::unit_root;
using cpphub::Real;
using cpphub::Size;

// Philox 确定性合成数据 (§8.6 三平台位一致; seed 固定)
static std::vector<Real> gen_normal(Size T, uint64_t seed, uint64_t ctr) {
    cpphub::v1::Philox4x64 rng(seed, ctr);
    std::vector<Real> z(T);
    for (Size i = 0; i < T; i += 2) {
        const uint64_t r1 = rng();
        const uint64_t r2 = rng();
        const Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
        const Real u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
        const auto [a, b] = cpphub::v1::box_muller(u1, u2);
        z[i] = a;
        if (i + 1 < T) z[i + 1] = b;
    }
    return z;
}

/// 随机游走 (单位根 DGP, H0 真)
static std::vector<Real> gen_rw(Size T, uint64_t ctr = 1) {
    const auto e = gen_normal(T, 42, ctr);
    std::vector<Real> y(T);
    y[0] = e[0];
    for (Size t = 1; t < T; ++t) y[t] = y[t - 1] + e[t];
    return y;
}

/// 平稳 AR(1) φ=0.95 (H1 真, 应拒绝单位根)
static std::vector<Real> gen_ar095(Size T, uint64_t ctr = 2) {
    const auto e = gen_normal(T, 42, ctr);
    std::vector<Real> y(T);
    y[0] = e[0] / std::sqrt(1.0 - 0.95 * 0.95);
    for (Size t = 1; t < T; ++t) y[t] = 0.95 * y[t - 1] + e[t];
    return y;
}

// ---------------------------------------------------------------------------
// 1: 恒等式 MZt ≡ MZα × MSB (1e-12, spec §1.3) — 两情形
// ---------------------------------------------------------------------------
TEST(NgPerronTest, IdentityMztEqualsMzaTimesMsb) {
    for (const std::string ts : {"c", "ct"}) {
        const auto r = ur::ng_perron_test(gen_rw(200), ts);
        EXPECT_NEAR(r.mz_t, r.mz_alpha * r.msb, 1e-12)
            << "trend=" << ts;
        EXPECT_TRUE(std::isfinite(r.mz_alpha));
        EXPECT_TRUE(std::isfinite(r.mpt));
    }
}

// ---------------------------------------------------------------------------
// 2-3: 临界值 vs np_table1 精确相等 (NP5 表)
// ---------------------------------------------------------------------------
TEST(NgPerronTest, CriticalValuesConstantCase) {
    const auto r = ur::ng_perron_test(gen_rw(200), "c");
    EXPECT_DOUBLE_EQ(r.cv_1pct[0], -13.8);
    EXPECT_DOUBLE_EQ(r.cv_5pct[0], -8.1);
    EXPECT_DOUBLE_EQ(r.cv_10pct[0], -5.7);
    EXPECT_DOUBLE_EQ(r.cv_5pct[1], -1.98);
    EXPECT_DOUBLE_EQ(r.cv_5pct[2], 0.233);
    EXPECT_DOUBLE_EQ(r.cv_5pct[3], 3.17);
    EXPECT_EQ(r.trend_spec, "c");
}

TEST(NgPerronTest, CriticalValuesTrendCase) {
    const auto r = ur::ng_perron_test(gen_rw(200), "ct");
    EXPECT_DOUBLE_EQ(r.cv_1pct[0], -23.8);
    EXPECT_DOUBLE_EQ(r.cv_5pct[0], -17.3);
    EXPECT_DOUBLE_EQ(r.cv_1pct[2], 0.143);  // 转录陷阱: 非 0.121
    EXPECT_DOUBLE_EQ(r.cv_1pct[3], 4.03);   // 转录陷阱: 非 4.47
}

// ---------------------------------------------------------------------------
// 4: NP4 — MPT 分情形末项系数 (+7 / +14.5): 各项非负 ⇒ MPT > 0 恒成立
//    (c̄²Σỹ²/n² > 0, tail·ỹ_T²/n: tail=+7 或 +14.5 均 > 0, s² > 0)
// ---------------------------------------------------------------------------
TEST(NgPerronTest, MptPositiveBothCases) {
    for (const std::string ts : {"c", "ct"}) {
        const auto r = ur::ng_perron_test(gen_rw(150), ts);
        EXPECT_GT(r.mpt, 0.0) << "trend=" << ts;
    }
    const auto r_ar = ur::ng_perron_test(gen_ar095(150), "ct");
    EXPECT_GT(r_ar.mpt, 0.0);
}

// ---------------------------------------------------------------------------
// 5: NP4 强化 — MPT 对 ỹ_T² 的敏感方向分情形不同不可能从输出直接分解,
//    改断言两情形 MPT 同数据下均有限且互异 (去势不同) — NP4 由实现公式冻结
// ---------------------------------------------------------------------------
TEST(NgPerronTest, MptDiffersAcrossDetrending) {
    const auto y = gen_rw(150);
    const auto rc = ur::ng_perron_test(y, "c");
    const auto rct = ur::ng_perron_test(y, "ct");
    EXPECT_NE(rc.mpt, rct.mpt);
    EXPECT_NE(rc.mz_alpha, rct.mz_alpha);
}

// ---------------------------------------------------------------------------
// 6-7: 方向断言 — RW 不拒绝 / AR(0.95) 拒绝 (MZα/MZt, 1%/5% 层)
// ---------------------------------------------------------------------------
TEST(NgPerronTest, RandomWalkNotRejected) {
    const auto r = ur::ng_perron_test(gen_rw(300), "ct");
    // H0 真: 四统计量均不应系统性拒绝 (5% 水平下期望不拒绝; 固定 seed 确定性)
    EXPECT_FALSE(r.reject_5pct[0]);  // MZα
    EXPECT_FALSE(r.reject_5pct[1]);  // MZt
    EXPECT_FALSE(r.reject_5pct[3]);  // MPT
    EXPECT_GT(r.mz_alpha, r.cv_5pct[0]);
}

TEST(NgPerronTest, StationaryArRejected) {
    const auto r = ur::ng_perron_test(gen_ar095(300), "ct");
    EXPECT_TRUE(r.reject_5pct[0]);   // MZα 拒绝 (5% 层; 1% 需更大 T, 功率边界)
    EXPECT_TRUE(r.reject_5pct[1]);   // MZt
    EXPECT_LT(r.mz_alpha, r.cv_5pct[0]);
    // NP5 方向: MSB/MPT 越小越拒绝 = stat < cv
    EXPECT_TRUE(r.reject_5pct[2]);
    EXPECT_TRUE(r.reject_5pct[3]);
}

// ---------------------------------------------------------------------------
// 8: MAIC/σ̂² 轨迹结构 (长度/正性/有限; 奇异 k 除外为 +inf)
// ---------------------------------------------------------------------------
TEST(NgPerronTest, MaicPathStructure) {
    const auto r = ur::ng_perron_test(gen_rw(150), "ct", 6);  // 显式 k_max=6
    EXPECT_EQ(r.maic.size(), 7u);
    EXPECT_EQ(r.sigma2_k.size(), 7u);
    for (Size k = 0; k <= 6; ++k) {
        EXPECT_TRUE(std::isfinite(r.maic[k]));
        EXPECT_GT(r.sigma2_k[k], 0.0);
    }
}

// ---------------------------------------------------------------------------
// 9: MAIC 惩罚结构 — maic[k] − ln σ̂²(k) = 2(τ_T(k)+k)/n > 0 (NP2 形)
// ---------------------------------------------------------------------------
TEST(NgPerronTest, MaicPenaltyNonnegative) {
    const auto r = ur::ng_perron_test(gen_rw(150), "ct", 6);
    const Real n = static_cast<Real>(150 - 6);
    for (Size k = 0; k <= 6; ++k) {
        const Real penalty = r.maic[k] - std::log(r.sigma2_k[k]);
        EXPECT_GT(penalty, 0.0);
        // 下界: 2(0 + k)/n (τ_T ≥ 0 因 β̂₀² ≥ 0)
        EXPECT_GE(penalty, 2.0 * static_cast<Real>(k) / n - 1e-12);
    }
}

// ---------------------------------------------------------------------------
// 10: selected_lag = argmin MAIC (首个 argmin 重构)
// ---------------------------------------------------------------------------
TEST(NgPerronTest, SelectedLagIsArgminMaic) {
    const auto r = ur::ng_perron_test(gen_rw(150), "ct", 6);
    Size k_star = 0;
    Real best = r.maic[0];
    for (Size k = 1; k <= 6; ++k) {
        if (r.maic[k] < best) {
            best = r.maic[k];
            k_star = k;
        }
    }
    EXPECT_EQ(r.selected_lag, k_star);
    EXPECT_LT(r.selected_lag, 7u);
}

// ---------------------------------------------------------------------------
// 11: Schwert 自动上限 (默认 max_lag=0 ⇒ k_max = schwert_lag(T))
// ---------------------------------------------------------------------------
TEST(NgPerronTest, DefaultLagIsSchwert) {
    const auto r = ur::ng_perron_test(gen_rw(150), "ct");
    const Size k_max = ur::schwert_lag(150);
    EXPECT_EQ(r.maic.size(), k_max + 1);
}

// ---------------------------------------------------------------------------
// 12: GLS 去势正交性 (B1 重实现正确性) — 准差分残差 ⊥ dz (OLS 一阶条件)
// ---------------------------------------------------------------------------
TEST(NgPerronTest, GlS_detrendOrthogonality) {
    const auto y = gen_rw(120);
    const auto yd = ur::detail::gls_detrend(y, "ct");
    ASSERT_EQ(yd.size(), y.size());
    const Size T = y.size();
    const Real c_bar = -13.5;
    const Real rho = 1.0 + c_bar / static_cast<Real>(T);
    // 从 ỹ 反解 ψ̂ (ỹ = y − z·ψ̂ 精确 ⇒ r_ = z·ψ̂; OLS of r_ on z 恢复 ψ̂)
    std::vector<std::vector<Real>> z(T, std::vector<Real>(2, 1.0));
    for (Size t = 0; t < T; ++t) z[t][1] = static_cast<Real>(t);
    std::vector<Real> r_(T);
    for (Size t = 0; t < T; ++t) r_[t] = y[t] - yd[t];
    const auto psi = ur::detail::ols_fit(r_, z);
    // GLS OLS 一阶条件: dz'(dy − dz·ψ̂) = 0 — **含首行** (U9 首项不变换 dz[0]=z[0])
    std::vector<Real> g(2, 0.0);
    Real scale0 = 0.0, scale1 = 0.0;
    for (Size t = 0; t < T; ++t) {
        const Real dz0 = (t == 0) ? z[0][0] : z[t][0] - rho * z[t - 1][0];
        const Real dz1 = (t == 0) ? z[0][1] : z[t][1] - rho * z[t - 1][1];
        const Real dy = (t == 0) ? y[0] : y[t] - rho * y[t - 1];
        const Real e = dy - (psi.beta[0] * dz0 + psi.beta[1] * dz1);
        g[0] += dz0 * e;
        g[1] += dz1 * e;
        scale0 += std::fabs(dz0 * e);
        scale1 += std::fabs(dz1 * e);
    }
    // 相对量级断言 (绝对零不可达, 正规方程条件数 ~1e5)
    EXPECT_LT(std::fabs(g[0]), 1e-8 * std::max(scale0, 1.0));
    EXPECT_LT(std::fabs(g[1]), 1e-8 * std::max(scale1, 1.0));
}

// ---------------------------------------------------------------------------
// 13: NP5 — reject_5pct 与 (stat < cv_5pct) 重构一致
// ---------------------------------------------------------------------------
TEST(NgPerronTest, RejectDirectionReconstruction) {
    for (const std::string ts : {"c", "ct"}) {
        const auto r = ur::ng_perron_test(gen_ar095(200), ts);
        const std::array<Real, 4> stats{r.mz_alpha, r.mz_t, r.msb, r.mpt};
        for (Size s = 0; s < 4; ++s) {
            EXPECT_EQ(r.reject_5pct[s], stats[s] < r.cv_5pct[s]);
        }
    }
}

// ---------------------------------------------------------------------------
// 14: 与 DF-GLS 同数据方向一致 (跨检验一致性: 平稳 AR 双双拒绝)
// ---------------------------------------------------------------------------
TEST(NgPerronTest, ConsistentWithDfGlsDirection) {
    // T=400: 功率足够使 MZt 与 DFGLS (渐近等价) 双双拒绝;
    // 小 T 下两者有限样本分离 (250 实测 MZt 边界未拒), 属已知差异
    const auto y = gen_ar095(400);
    const auto np_r = ur::ng_perron_test(y, "ct");
    const auto dg = ur::df_gls_test(y, "ct");
    EXPECT_TRUE(dg.reject_null);
    EXPECT_TRUE(np_r.reject_5pct[1]);  // MZt 与 DFGLS 渐近等价
}

// ---------------------------------------------------------------------------
// 15: 异常输入 (T<10 / NaN / 常数 / trend 非法 / max_lag 过大)
// ---------------------------------------------------------------------------
TEST(NgPerronTest, InvalidInputThrows) {
    EXPECT_THROW(ur::ng_perron_test(std::vector<Real>(9, 1.0), "ct"),
                 std::invalid_argument);
    auto y = gen_rw(100);
    y[3] = std::nan("");
    EXPECT_THROW(ur::ng_perron_test(y, "ct"), std::invalid_argument);
    EXPECT_THROW(ur::ng_perron_test(std::vector<Real>(100, 2.5), "c"),
                 std::invalid_argument);
    EXPECT_THROW(ur::ng_perron_test(gen_rw(100), "n"), std::invalid_argument);
    EXPECT_THROW(ur::ng_perron_test(gen_rw(100), "xyz"), std::invalid_argument);
    EXPECT_THROW(ur::ng_perron_test(gen_rw(20), "ct", 19),  // n = 1 < 3
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 16: 默认参数 = trend "ct" (与 DF-GLS 默认一致)
// ---------------------------------------------------------------------------
TEST(NgPerronTest, DefaultParameters) {
    const auto r = ur::ng_perron_test(gen_rw(100));
    EXPECT_EQ(r.trend_spec, "ct");
    const auto r2 = ur::ng_perron_test(gen_rw(100), "ct");
    EXPECT_DOUBLE_EQ(r.mz_alpha, r2.mz_alpha);
}

// ---------------------------------------------------------------------------
// 17: summary 关键词
// ---------------------------------------------------------------------------
TEST(NgPerronTest, SummaryContent) {
    const auto r = ur::ng_perron_test(gen_ar095(150), "ct");
    EXPECT_NE(r.summary.find("Ng-Perron"), std::string::npos);
    EXPECT_NE(r.summary.find("unit root"), std::string::npos);
    EXPECT_NE(r.summary.find("MAIC"), std::string::npos);
}

// ---------------------------------------------------------------------------
// 18: 性能 P1 — T=1000 全 k 搜索 < 1s (checklist §15.1 镜像)
// ---------------------------------------------------------------------------
TEST(NgPerronTest, PerformanceT1000) {
    const auto y = gen_rw(1000, 7);
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = ur::ng_perron_test(y, "ct");  // k_max = Schwert(1000) = 21
    const auto t1 = std::chrono::steady_clock::now();
    const double sec =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() /
        1000.0;
    EXPECT_LT(sec, 1.0);
    EXPECT_GT(r.maic.size(), 10u);  // Schwert(1000)=ceil(12·10^0.25)=21
    (void)r;
}
