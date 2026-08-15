// =============================================================================
// test_variance_ratio_test.cpp - 方差比检验测试 (20 用例, spec §3.5 测试矩阵)
//
// 基准: tests/unit/timeseries/unit_root_baseline.inc (arch 8.0.0 自动生成,
//       脚本 tests/fixtures/timeseries/verify_unit_root.py, 硬编码策略)
//
// VR_CASES 布局 (12 行, 顺序 series×k×(rob,deb)):
//   idx 0/1/2:   rw/k=2 × (rob=T,deb=T)/(rob=T,deb=F)/(rob=F,deb=T)
//   idx 3/4/5:   rw/k=5 × 同上
//   idx 6/7/8:   ar/k=2 × 同上;  idx 9/10/11: ar/k=5 × 同上
// 映射: 我们的接口一次算 Z1+Z2 — rob=T 行对照 z2, rob=F 行对照 z1,
//       vr 只依赖 deb (rob=T/rob=F 同 deb 行的 vr 相同)
//
// 容差 (spec §7.1):
//   - vr/statistic vs arch: 1e-10; p 值: 1e-12
//   - Chow-Denning 联合检验: 1e-8 (R vrtest 未装, 用 arch Z2 分量 + SMM 公式重构)
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/unit_root/variance_ratio_test.hpp"
#include "unit_root_baseline.inc"

namespace ur = cpphub::v1::timeseries::unit_root;
using cpphub::Real;
using cpphub::Size;

// VR 基准在价格序列 P_RW/P_AR (TP=501) 上生成 — 非 Y_RW/Y_AR (T=250)!
static std::vector<Real> p_rw() {
    return {ur::baseline::P_RW, ur::baseline::P_RW + ur::baseline::TP};
}
static std::vector<Real> p_ar() {
    return {ur::baseline::P_AR, ur::baseline::P_AR + ur::baseline::TP};
}
static std::vector<Real> tiny_p() {
    return {ur::baseline::TINY_P, ur::baseline::TINY_P + ur::baseline::TINY_P_N};
}

// rob=True 行 → z2 分量; vr 一致
static void check_robust(const ur::VarianceRatioResult& r,
                         const ur::baseline::VrCase& c) {
    EXPECT_NEAR(r.vr_statistic, c.vr, 1e-10);
    EXPECT_NEAR(r.z2_statistic, c.stat, 1e-10);
    EXPECT_NEAR(r.z2_p_value, c.p, 1e-12);
}
// rob=False 行 → z1 分量; vr 一致
static void check_homosk(const ur::VarianceRatioResult& r,
                         const ur::baseline::VrCase& c) {
    EXPECT_NEAR(r.vr_statistic, c.vr, 1e-10);
    EXPECT_NEAR(r.z1_statistic, c.stat, 1e-10);
    EXPECT_NEAR(r.z1_p_value, c.p, 1e-12);
}

// ---------------------------------------------------------------------------
// 1-6: rw 序列 (H0 真, 不应拒绝), k=2/5 × deb/rob 组合
// ---------------------------------------------------------------------------
TEST(VarianceRatioTest, RwK2DebiasedRobust) {
    check_robust(ur::variance_ratio_test(p_rw(), 2, true), ur::baseline::VR_CASES[0]);
}

TEST(VarianceRatioTest, RwK2UndebiasedRobust) {
    check_robust(ur::variance_ratio_test(p_rw(), 2, false), ur::baseline::VR_CASES[1]);
}

TEST(VarianceRatioTest, RwK2DebiasedHomoskedastic) {
    check_homosk(ur::variance_ratio_test(p_rw(), 2, true), ur::baseline::VR_CASES[2]);
}

TEST(VarianceRatioTest, RwK5DebiasedRobust) {
    check_robust(ur::variance_ratio_test(p_rw(), 5, true), ur::baseline::VR_CASES[3]);
}

TEST(VarianceRatioTest, RwK5UndebiasedRobust) {
    check_robust(ur::variance_ratio_test(p_rw(), 5, false), ur::baseline::VR_CASES[4]);
}

TEST(VarianceRatioTest, RwK5DebiasedHomoskedastic) {
    check_homosk(ur::variance_ratio_test(p_rw(), 5, true), ur::baseline::VR_CASES[5]);
}

// ---------------------------------------------------------------------------
// 7-12: ar 序列 (H0 假) — VR 显著 < 1 (负自相关)
// ---------------------------------------------------------------------------
TEST(VarianceRatioTest, ArK2DebiasedRobust) {
    const auto r = ur::variance_ratio_test(p_ar(), 2, true);
    check_robust(r, ur::baseline::VR_CASES[6]);
    // p=0.0899 > 0.05: Z2 稳健统计量单 horizon 不拒绝 (小样本稳健代价)
    EXPECT_FALSE(r.reject_null);
}

TEST(VarianceRatioTest, ArK2UndebiasedRobust) {
    check_robust(ur::variance_ratio_test(p_ar(), 2, false), ur::baseline::VR_CASES[7]);
}

TEST(VarianceRatioTest, ArK2DebiasedHomoskedastic) {
    check_homosk(ur::variance_ratio_test(p_ar(), 2, true), ur::baseline::VR_CASES[8]);
}

TEST(VarianceRatioTest, ArK5DebiasedRobust) {
    const auto r = ur::variance_ratio_test(p_ar(), 5, true);
    check_robust(r, ur::baseline::VR_CASES[9]);
    EXPECT_TRUE(r.reject_null);  // p=0.0159 < 0.05
}

TEST(VarianceRatioTest, ArK5UndebiasedRobust) {
    check_robust(ur::variance_ratio_test(p_ar(), 5, false), ur::baseline::VR_CASES[10]);
}

TEST(VarianceRatioTest, ArK5DebiasedHomoskedastic) {
    check_homosk(ur::variance_ratio_test(p_ar(), 5, true), ur::baseline::VR_CASES[11]);
}

// ---------------------------------------------------------------------------
// 13-14: trend="n" 零漂移 (mu=0, U17-demean)
// ---------------------------------------------------------------------------
TEST(VarianceRatioTest, TrendN_K2) {
    check_robust(ur::variance_ratio_test(p_rw(), 2, true, "n"),
                 ur::baseline::VR_N_CASES[0]);
}

TEST(VarianceRatioTest, TrendN_K5) {
    check_robust(ur::variance_ratio_test(p_rw(), 5, true, "n"),
                 ur::baseline::VR_N_CASES[1]);
}

// ---------------------------------------------------------------------------
// 15-16: tiny 手算分量白盒 (P=11 点带噪声序列, k=2)
//   deb=T/rob=T: vr=VRD, z2=STATD (与 arch 交叉验证值逐位一致)
//   deb=F/homosk: vr=VR, z1=STATH
//   theta 一致性: sqrt(nq)·(VRD-1)/sqrt(THETA) == ARCHSTAT
// ---------------------------------------------------------------------------
TEST(VarianceRatioTest, TinyDebiasedRobust) {
    const auto r = ur::variance_ratio_test(tiny_p(), 2, true);
    EXPECT_NEAR(r.vr_statistic, ur::baseline::VR_TINY_VRD, 1e-15);
    EXPECT_NEAR(r.vr_statistic, ur::baseline::VR_TINY_ARCHVR, 1e-15);
    EXPECT_NEAR(r.z2_statistic, ur::baseline::VR_TINY_STATD, 1e-12);
    EXPECT_NEAR(r.z2_statistic, ur::baseline::VR_TINY_ARCHSTAT, 1e-12);
    // 内部一致性: 基准 theta 与 arch 统计量满足同一公式
    const Real nq = static_cast<Real>(ur::baseline::TINY_P_N - 1);
    const Real recomputed =
        std::sqrt(nq) * (ur::baseline::VR_TINY_VRD - 1.0) /
        std::sqrt(ur::baseline::VR_TINY_THETA);
    EXPECT_NEAR(recomputed, ur::baseline::VR_TINY_ARCHSTAT, 1e-12);
}

TEST(VarianceRatioTest, TinyUndebiasedHomoskedastic) {
    const auto r = ur::variance_ratio_test(tiny_p(), 2, false);
    EXPECT_NEAR(r.vr_statistic, ur::baseline::VR_TINY_VR, 1e-15);
    EXPECT_NEAR(r.z1_statistic, ur::baseline::VR_TINY_STATH, 1e-12);
    EXPECT_NEAR(r.z1_statistic, ur::baseline::VR_TINY_ARCHSTAT_H, 1e-12);
}

// ---------------------------------------------------------------------------
// 17: CLM debiased 参考字段 — 与 use_debiased 无关, 恒为 arch 默认组合
//     (robust+debiased) 的 Z2
// ---------------------------------------------------------------------------
TEST(VarianceRatioTest, ClmFieldAlwaysDebiased) {
    const auto r_undeb = ur::variance_ratio_test(p_rw(), 2, false);
    EXPECT_NEAR(r_undeb.clm_debiased_statistic, ur::baseline::VR_CASES[0].stat,
                1e-10);
    EXPECT_NEAR(r_undeb.clm_debiased_p_value, ur::baseline::VR_CASES[0].p, 1e-12);
    const auto r_deb = ur::variance_ratio_test(p_rw(), 2, true);
    EXPECT_NEAR(r_deb.clm_debiased_statistic, r_deb.z2_statistic, 0.0);
    EXPECT_DOUBLE_EQ(r_deb.clm_debiased_p_value, r_deb.z2_p_value);
}

// ---------------------------------------------------------------------------
// 18: 单 horizon Chow-Denning 退化 (m=1): p = 2-2·Phi(|Z2|) = z2_p
// ---------------------------------------------------------------------------
TEST(VarianceRatioTest, ChowDenningSingleKDegenerate) {
    const auto r = ur::variance_ratio_test(p_ar(), 5, true);
    ASSERT_EQ(r.chow_denning_stats.size(), 1u);
    EXPECT_NEAR(r.chow_denning_stats[0], ur::baseline::VR_CASES[9].stat, 1e-10);
    EXPECT_NEAR(r.chow_denning_p_value, r.z2_p_value, 1e-15);
}

// ---------------------------------------------------------------------------
// 19-20: 多 horizon Chow-Denning 联合检验 (U16: Z2 + SMM(m,inf))
// ---------------------------------------------------------------------------
TEST(VarianceRatioTest, MultiChowDenningArRejects) {
    const std::vector<Size> ks = {2, 5};
    const auto r = ur::variance_ratio_test_multi(p_ar(), ks, true);
    ASSERT_EQ(r.chow_denning_stats.size(), 2u);
    EXPECT_NEAR(r.chow_denning_stats[0], ur::baseline::VR_CASES[6].stat, 1e-10);
    EXPECT_NEAR(r.chow_denning_stats[1], ur::baseline::VR_CASES[9].stat, 1e-10);
    // SMM(2,inf) 重构: p = 1 - [2·Phi(|CD|)-1]^2 (1e-8, spec §7.1)
    const Real cd = std::max(std::fabs(ur::baseline::VR_CASES[6].stat),
                             std::fabs(ur::baseline::VR_CASES[9].stat));
    const Real phi = 0.5 * std::erfc(-cd / std::sqrt(2.0));
    const Real p_ref = 1.0 - (2.0 * phi - 1.0) * (2.0 * phi - 1.0);
    EXPECT_NEAR(r.chow_denning_p_value, p_ref, 1e-8);
    EXPECT_TRUE(r.reject_null);
    // 标量字段 = 首 horizon (k=2); vr_list/k_values 填充
    EXPECT_EQ(r.k, 2u);
    EXPECT_EQ(r.k_values.size(), 2u);
    ASSERT_EQ(r.vr_list.size(), 2u);
    EXPECT_NEAR(r.vr_list[0], ur::baseline::VR_CASES[6].vr, 1e-10);
    EXPECT_NEAR(r.vr_list[1], ur::baseline::VR_CASES[9].vr, 1e-10);
    EXPECT_NEAR(r.z2_statistic, ur::baseline::VR_CASES[6].stat, 1e-10);
}

TEST(VarianceRatioTest, MultiChowDenningRwNotRejected) {
    const auto r = ur::variance_ratio_test_multi(p_rw(), {2, 5}, true);
    const Real cd = std::max(std::fabs(ur::baseline::VR_CASES[0].stat),
                             std::fabs(ur::baseline::VR_CASES[3].stat));
    const Real phi = 0.5 * std::erfc(-cd / std::sqrt(2.0));
    const Real p_ref = 1.0 - (2.0 * phi - 1.0) * (2.0 * phi - 1.0);
    EXPECT_NEAR(r.chow_denning_p_value, p_ref, 1e-8);
    EXPECT_FALSE(r.reject_null);
}

// ---------------------------------------------------------------------------
// 21: 平移不变性 — 差分消除常数水平
// ---------------------------------------------------------------------------
TEST(VarianceRatioTest, InvarianceToLevelShift) {
    std::vector<Real> shifted = p_rw();
    for (Real& v : shifted) v += 1000.0;
    const auto a = ur::variance_ratio_test(p_rw(), 2, true);
    const auto b = ur::variance_ratio_test(shifted, 2, true);
    EXPECT_NEAR(b.vr_statistic, a.vr_statistic, 1e-10);
    EXPECT_NEAR(b.z2_statistic, a.z2_statistic, 1e-8);
    EXPECT_NEAR(b.z1_statistic, a.z1_statistic, 1e-8);
}

// ---------------------------------------------------------------------------
// 22: 异常输入
// ---------------------------------------------------------------------------
TEST(VarianceRatioTest, InvalidInputs) {
    EXPECT_THROW(ur::variance_ratio_test(p_rw(), 1), std::invalid_argument);   // k<2
    EXPECT_THROW(ur::variance_ratio_test(p_rw(), 2, true, "ct"),
                 std::invalid_argument);                                        // trend
    EXPECT_THROW(ur::variance_ratio_test(tiny_p(), 20), std::invalid_argument); // 样本不足
    EXPECT_THROW(ur::variance_ratio_test({}, 2), std::invalid_argument);        // 空序列
    EXPECT_THROW(ur::variance_ratio_test(std::vector<Real>(20, 5.0), 2),
                 std::runtime_error);  // 常数序列 → 零方差
    EXPECT_THROW(ur::variance_ratio_test_multi(p_rw(), {}),
                 std::invalid_argument);  // 空 k_list
}
