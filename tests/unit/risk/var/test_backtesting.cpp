#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/risk/var/backtesting.hpp"

namespace cpphub {
inline namespace v1 {

TEST(KupiecPOF, CorrectCoverage) {
    auto result = KupiecPOF::test(2, 250, 0.99);
    EXPECT_EQ(result.n_violations, 2u);
    EXPECT_EQ(result.n_observations, 250u);
    EXPECT_NEAR(result.expected_rate, 0.01, 1e-10);
    EXPECT_GT(result.p_value, 0.05);
    EXPECT_FALSE(result.reject_null);
}

TEST(KupiecPOF, UnderestimatedRisk) {
    auto result = KupiecPOF::test(10, 250, 0.99);
    EXPECT_EQ(result.n_violations, 10u);
    EXPECT_LT(result.p_value, 0.01);
    EXPECT_TRUE(result.reject_null);
}

TEST(KupiecPOF, OverestimatedRisk) {
    auto result = KupiecPOF::test(0, 250, 0.99);
    EXPECT_EQ(result.n_violations, 0u);
    EXPECT_EQ(result.violation_rate, 0.0);
    EXPECT_GT(result.p_value, 0.01);
    EXPECT_LT(result.p_value, 0.05);
}

TEST(KupiecPOF, FromSeries) {
    std::vector<Real> var_series(250, 100.0);
    std::vector<Real> losses(250, 50.0);
    losses[0] = 110.0;
    losses[1] = 120.0;

    auto result = KupiecPOF::test(var_series, losses, 0.99);
    EXPECT_EQ(result.n_violations, 2u);

    auto direct = KupiecPOF::test(2, 250, 0.99);
    EXPECT_NEAR(result.p_value, direct.p_value, 1e-10);
    EXPECT_EQ(result.reject_null, direct.reject_null);
}

TEST(ChristoffersenIID, IndependentBreaches) {
    std::vector<Real> var_series(500, 100.0);
    std::vector<Real> losses(500, 50.0);
    for (Size i = 0; i < 500; i += 50) {
        losses[i] = 110.0;
    }

    auto result = ChristoffersenIID::test(var_series, losses);
    EXPECT_GT(result.p_value, 0.05);
}

TEST(ChristoffersenIID, ClusteredBreaches) {
    std::vector<Real> var_series(500, 100.0);
    std::vector<Real> losses(500, 50.0);
    for (Size i = 200; i < 210; ++i) {
        losses[i] = 110.0;
    }

    auto result = ChristoffersenIID::test(var_series, losses);
    EXPECT_LT(result.p_value, 0.05);
}

TEST(ChristoffersenIID, JointTest) {
    std::vector<Real> var_series(500, 100.0);
    std::vector<Real> losses(500, 50.0);
    for (Size i = 200; i < 210; ++i) {
        losses[i] = 110.0;
    }

    auto ind = ChristoffersenIID::test(var_series, losses);
    auto joint = ChristoffersenIID::joint_test(var_series, losses, 0.99);

    // 收紧 (RISK-023): 联合检验 = 覆盖率 + 独立性, 聚簇序列两者都必须拒绝;
    // LR_cc = LR_pof + LR_ind >= LR_ind ⇒ 联合 p 不应大于独立性 p
    EXPECT_TRUE(ind.reject_null);
    EXPECT_TRUE(joint.reject_null);
    EXPECT_LE(joint.p_value, ind.p_value + 1e-12);
}

TEST(BaselTrafficLight, GreenZone) {
    for (Size v = 0; v <= 4; ++v) {
        auto result = BaselTrafficLight::assess(v, 250);
        EXPECT_EQ(result.zone, BaselZone::Green);
        EXPECT_EQ(result.capital_multiplier, 3.0);
    }
}

TEST(BaselTrafficLight, YellowZone) {
    Real expected_mult[] = {3.4, 3.5, 3.65, 3.75, 3.85};
    for (Size v = 5; v <= 9; ++v) {
        auto result = BaselTrafficLight::assess(v, 250);
        EXPECT_EQ(result.zone, BaselZone::Yellow);
        EXPECT_NEAR(result.capital_multiplier, expected_mult[v - 5], 1e-10);
    }
}

TEST(BaselTrafficLight, RedZone) {
    for (Size v = 10; v <= 15; ++v) {
        auto result = BaselTrafficLight::assess(v, 250);
        EXPECT_EQ(result.zone, BaselZone::Red);
        EXPECT_EQ(result.capital_multiplier, 4.0);
    }
}

// ============================================================================
// RISK-023 / RISK-024 / RISK-028 回归
//   023: 联合检验 LR_cc = LR_pof + LR_ind ~ χ²(2) (独立性分量此前被丢弃)
//   024: 似然连乘下溢 → 对数似然 (0·log0 = 0) + 守卫重写
//   028: χ² CDF 精确实现 (原 Wilson-Hilferty 近似临界带偏差 ~3e-3, 可翻转结论)
//
// 参考值策略 (双重来源, 避免用实现验证实现):
//   - 对数似然: 测试内按教科书公式独立重建 (对数空间)
//   - χ² 上尾闭式: k=1 → Q(x) = erfc(√(x/2)); k=2 → Q(x) = exp(-x/2)
// ============================================================================

namespace {

/// 影子参考: 0·log(y) = 0 约定
Real ref_xlog1m(Real x, Real y) { return x == 0.0 ? 0.0 : x * std::log1p(-y); }
Real ref_xlogy(Real x, Real y)  { return x == 0.0 ? 0.0 : x * std::log(y); }

/// 转移计数 (含首期) 参考
struct RefCounts {
    Size n = 0, v = 0, n00 = 0, n01 = 0, n10 = 0, n11 = 0;
};

RefCounts ref_counts(const std::vector<Real>& var_series,
                     const std::vector<Real>& losses) {
    RefCounts c;
    c.n = std::min(var_series.size(), losses.size());
    for (Size i = 0; i < c.n; ++i) {
        if (losses[i] > var_series[i]) ++c.v;
        if (i == 0) continue;
        bool prev = losses[i-1] > var_series[i-1];
        bool cur  = losses[i] > var_series[i];
        if (!prev && !cur) ++c.n00;
        if (!prev &&  cur) ++c.n01;
        if ( prev && !cur) ++c.n10;
        if ( prev &&  cur) ++c.n11;
    }
    return c;
}

/// Kupiec LR (对数空间, 独立实现)
Real ref_kupiec_lr(Size v, Size N, Real confidence) {
    Real p = 1.0 - confidence;
    Real vr = static_cast<Real>(v) / static_cast<Real>(N);
    Real logL0 = ref_xlog1m(static_cast<Real>(N - v), p) + ref_xlogy(static_cast<Real>(v), p);
    Real logL1 = ref_xlog1m(static_cast<Real>(N - v), vr) + ref_xlogy(static_cast<Real>(v), vr);
    return 2.0 * (logL1 - logL0);
}

/// Christoffersen 独立性 LR (对数空间, 独立实现)
Real ref_ind_lr(const RefCounts& c) {
    Real m = static_cast<Real>(c.n - 1);
    Real pi  = static_cast<Real>(c.n01 + c.n11) / m;
    Size n0 = c.n00 + c.n01, n1 = c.n10 + c.n11;
    Real pi0 = (n0 > 0) ? static_cast<Real>(c.n01) / static_cast<Real>(n0) : 0.0;
    Real pi1 = (n1 > 0) ? static_cast<Real>(c.n11) / static_cast<Real>(n1) : 0.0;
    Real log_iid = ref_xlog1m(m - static_cast<Real>(c.n01 + c.n11), pi)
                 + ref_xlogy(static_cast<Real>(c.n01 + c.n11), pi);
    Real log_non = ref_xlog1m(static_cast<Real>(c.n00), pi0)
                 + ref_xlogy(static_cast<Real>(c.n01), pi0)
                 + ref_xlog1m(static_cast<Real>(c.n10), pi1)
                 + ref_xlogy(static_cast<Real>(c.n11), pi1);
    return 2.0 * (log_non - log_iid);
}

/// 构造 n 期序列: 违约位置由 breach_idx 指定 (违约损失 110 > VaR 100)
void build_series(Size n, const std::vector<Size>& breach_idx,
                  std::vector<Real>& var_series, std::vector<Real>& losses) {
    var_series.assign(n, 100.0);
    losses.assign(n, 50.0);
    for (Size i : breach_idx) losses[i] = 110.0;
}

} // namespace

// ---------------------------------------------------------------------------
// RISK-028: χ² CDF 精确性
// ---------------------------------------------------------------------------

TEST(Chi2Cdf, MatchesExactGammaAndClosedForm) {
    // k=1 有闭式: CDF(x) = erf(√(x/2))
    EXPECT_NEAR(detail::chi2_cdf(2.0, 1), std::erf(std::sqrt(1.0)), 1e-14);
    EXPECT_NEAR(detail::chi2_cdf(0.5, 1), std::erf(std::sqrt(0.25)), 1e-14);
    // 标准临界值: CDF(x_crit) = 1 - α (旧 Wilson-Hilferty 偏差 ~3e-3 → 必红)
    EXPECT_NEAR(detail::chi2_cdf(3.841459, 1), 0.95, 1e-6);
    EXPECT_NEAR(detail::chi2_cdf(2.705543, 1), 0.90, 1e-6);
    EXPECT_NEAR(detail::chi2_cdf(6.634897, 1), 0.99, 1e-6);
    EXPECT_NEAR(detail::chi2_cdf(5.991465, 2), 0.95, 1e-6);
    EXPECT_NEAR(detail::chi2_cdf(9.487729, 4), 0.95, 1e-6);
    EXPECT_EQ(detail::chi2_cdf(0.0, 1), 0.0);
    EXPECT_EQ(detail::chi2_cdf(-1.0, 2), 0.0);  // 数值负 LR → 0 (不抛异常)
}

TEST(KupiecPOF, PValueMatchesExactChi2Tail) {
    // RISK-028: 6/250 @99% → LR ≈ 3.556; 旧近似 p≈0.0569 vs 精确 p≈0.0593
    const Size v = 6, N = 250;
    auto result = KupiecPOF::test(v, N, 0.99);
    double lr = ref_kupiec_lr(v, N, 0.99);
    double expected = std::erfc(std::sqrt(lr / 2.0));  // χ²(1) 上尾闭式
    EXPECT_NEAR(result.p_value, expected, 1e-10) << "LR=" << lr;
}

// ---------------------------------------------------------------------------
// RISK-024: 下溢与静默失效
// ---------------------------------------------------------------------------

TEST(KupiecPOF, LargeSampleNoUnderflow) {
    // L0 = 0.99^2500·0.01^2500 ≈ 1e-5011 = 0 (精确), L1 = 0.5^5000 ≈ 1e-1505 = 0
    // 旧实现: log(0/0) = NaN → p=NaN → reject_null 恒 false (50% 违约率都不拒绝)
    auto result = KupiecPOF::test(2500, 5000, 0.99);
    EXPECT_FALSE(std::isnan(result.p_value)) << "旧实现返回 p=NaN (静默不拒绝)";
    EXPECT_TRUE(std::isfinite(result.p_value));
    EXPECT_TRUE(result.reject_null);
    EXPECT_LT(result.p_value, 1e-10);
}

TEST(ChristoffersenIID, LargeSampleNoUnderflow) {
    // n=4000, 200 连簇: pi≈0.05 → L_IID ≈ 1e-343 → 下溢为 0 (次正规/精确 0)
    // 旧实现: 守卫 (L_IID>0 && L_nonIID>0) 命中 → LR_ind 置 0 → p=1 (静默不拒绝聚簇)
    std::vector<Size> idx;
    for (Size i = 1000; i < 1200; ++i) idx.push_back(i);
    std::vector<Real> var_series, losses;
    build_series(4000, idx, var_series, losses);

    auto result = ChristoffersenIID::test(var_series, losses);
    EXPECT_TRUE(std::isfinite(result.p_value));
    EXPECT_TRUE(result.reject_null) << "200 连簇在 4000 观测下必须拒绝独立性";
    EXPECT_LT(result.p_value, 1e-10);
}

TEST(ChristoffersenIID, PValueMatchesExactChi2Tail) {
    // 温和聚簇: LR_ind ≈ 3.71 → p ≈ 0.054 (旧 Wilson-Hilferty 给 0.0513, 偏差 3e-3)
    // 尾部对比需挑 p 未被浮点下限 (≈1e-16) 截断的区间, 否则新旧都会被归零而失去区分力
    std::vector<Size> idx{100, 101, 200, 300, 400, 450};  // 一对相邻 + 4 孤立
    std::vector<Real> var_series, losses;
    build_series(500, idx, var_series, losses);

    auto result = ChristoffersenIID::test(var_series, losses);
    double lr = ref_ind_lr(ref_counts(var_series, losses));
    double expected = std::erfc(std::sqrt(lr / 2.0));
    EXPECT_GT(expected, 1e-3) << "用例自检: 期望 p 必须远离浮点下限 (LR_ind=" << lr << ")";
    EXPECT_NEAR(result.p_value, expected, 1e-10) << "LR_ind=" << lr;
}

// ---------------------------------------------------------------------------
// RISK-023: 联合检验分量
// ---------------------------------------------------------------------------

TEST(ChristoffersenIID, JointTestSumsLRComponents) {
    // LR_cc = LR_pof + LR_ind ~ χ²(2); χ²(2) 上尾闭式 Q(x) = exp(-x/2)
    // 旧实现: LR_cc_total = LR_cc (仅覆盖率) → 独立性证据完全丢失
    std::vector<Size> idx{100, 101, 200, 300, 400, 450};  // 一对相邻 + 4 孤立
    std::vector<Real> var_series, losses;
    build_series(500, idx, var_series, losses);

    auto joint = ChristoffersenIID::joint_test(var_series, losses, 0.99);
    RefCounts c = ref_counts(var_series, losses);
    double lr_cc = ref_kupiec_lr(c.v, c.n, 0.99) + ref_ind_lr(c);
    EXPECT_NEAR(joint.p_value, std::exp(-lr_cc / 2.0), 1e-12) << "LR_cc=" << lr_cc;
    EXPECT_FALSE(joint.reject_null) << "温和聚簇 (LR_cc≈3.9) 不应拒绝";
}

TEST(ChristoffersenIID, JointTestClusteredAtExpectedRate) {
    // 评审失败用例: 违约率恰等于期望 5%, 但集中成 5 个 5 连簇
    // LR_pof ≈ 0 (覆盖率无证据), 全部证据在 LR_ind
    // 旧实现: LR_cc = LR_pof ≈ 0 → p≈1 → 静默不拒绝 (正确检验应拒绝)
    std::vector<Size> idx;
    for (Size b = 0; b < 5; ++b) {
        for (Size k = 0; k < 5; ++k) idx.push_back(50 + 100 * b + k);
    }
    std::vector<Real> var_series, losses;
    build_series(500, idx, var_series, losses);

    auto joint = ChristoffersenIID::joint_test(var_series, losses, 0.95);
    EXPECT_EQ(joint.n_violations, 25u) << "25/500 = 5% = 期望违约率";
    EXPECT_TRUE(joint.reject_null) << "等率聚集序列必须拒绝 (证据全在 LR_ind)";
    EXPECT_LT(joint.p_value, 1e-6);

    auto ind = ChristoffersenIID::test(var_series, losses);
    EXPECT_TRUE(ind.reject_null);
}

TEST(ChristoffersenIID, JointTestIndependentCorrectCoverage) {
    // 对照: 5 个孤立违约 (5/500 = 期望 1%@99% ≈ 0 偏差), 独立性成立 → 不拒绝
    std::vector<Size> idx{100, 200, 300, 400, 450};
    std::vector<Real> var_series, losses;
    build_series(500, idx, var_series, losses);

    auto joint = ChristoffersenIID::joint_test(var_series, losses, 0.99);
    RefCounts c = ref_counts(var_series, losses);
    double lr_cc = ref_kupiec_lr(c.v, c.n, 0.99) + ref_ind_lr(c);
    EXPECT_NEAR(joint.p_value, std::exp(-lr_cc / 2.0), 1e-12) << "LR_cc=" << lr_cc;
    EXPECT_FALSE(joint.reject_null);
}

TEST(ChristoffersenIID, ViolationCountsConsistentAcrossTests) {
    // 首期即违约: 旧 test() 返回 n01+n11 (漏计首期) 与 Kupiec/joint 口径不一致
    std::vector<Size> idx{0, 1, 50};
    std::vector<Real> var_series, losses;
    build_series(100, idx, var_series, losses);

    auto k = KupiecPOF::test(var_series, losses, 0.99);
    auto ind = ChristoffersenIID::test(var_series, losses);
    auto joint = ChristoffersenIID::joint_test(var_series, losses, 0.99);
    EXPECT_EQ(k.n_violations, 3u);
    EXPECT_EQ(ind.n_violations, 3u) << "旧实现返回 2 (漏计首期违约)";
    EXPECT_EQ(joint.n_violations, 3u);
    EXPECT_NEAR(ind.violation_rate, k.violation_rate, 1e-15);
}

// ---------------------------------------------------------------------------
// 输入边界 (原实现静默产出 NaN)
// ---------------------------------------------------------------------------

TEST(ChristoffersenIID, RejectsTooFewObservations) {
    std::vector<Real> var1{100.0}, loss1{50.0};
    EXPECT_THROW(ChristoffersenIID::test(var1, loss1), std::invalid_argument);
    EXPECT_THROW(ChristoffersenIID::joint_test(var1, loss1, 0.99), std::invalid_argument);
    std::vector<Real> empty;
    EXPECT_THROW(ChristoffersenIID::test(empty, empty), std::invalid_argument);
}

TEST(KupiecPOF, RejectsInvalidInputs) {
    EXPECT_THROW(KupiecPOF::test(3, 2, 0.99), std::invalid_argument);   // v > N
    EXPECT_THROW(KupiecPOF::test(0, 0, 0.99), std::invalid_argument);   // N = 0
    EXPECT_THROW(KupiecPOF::test(1, 100, 1.0), std::invalid_argument);  // confidence ∉ (0,1)
}

} // namespace v1
} // namespace cpphub
