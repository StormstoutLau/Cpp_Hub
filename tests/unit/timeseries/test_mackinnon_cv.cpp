// =============================================================================
// test_mackinnon_cv.cpp - MacKinnon 2010 临界值/p 值 + KPSS 插值测试
// (12 用例, spec §3.0.2 测试矩阵)
//
// Phase 7B v1.6 M2 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 基准: tests/unit/timeseries/unit_root_baseline.inc MKCV/MKP/KPSSCRIT
//       (arch 8.0.0 mackinnoncrit/mackinnonp/kpss_crit 自动生成)
//
// 幻觉点覆盖:
//   U3: 4 系数 3 次多项式 CV = c3/T³ + c2/T² + c1/T + c0 (T=0 → c0)
//   U13: ADF/PP p 值分 smallp(二次)/largep(三次) 段, 系数升幂
//   U7: DF-GLS 独立临界值表 (非 ERS 1996 原表)
//   U6: KPSS 线性插值 + np.interp 越界 clamp
// =============================================================================
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/timeseries/unit_root/mackinnon_cv.hpp"
#include "unit_root_baseline.inc"

namespace ur = cpphub::v1::timeseries::unit_root;
using cpphub::Real;
using cpphub::Size;

// ---------------------------------------------------------------------------
// 1-3. ADF 临界值: nc/c/ct × T={233, 100, 0(∞)} × 3 水平 (1e-12)
// ---------------------------------------------------------------------------
TEST(MackinnonCv, AdfCriticalValueN) {
    for (int i = 0; i < 3; ++i) {
        const Real p = i == 0 ? 0.01 : (i == 1 ? 0.05 : 0.10);
        EXPECT_NEAR(ur::mackinnon_critical_value("adf", "nc", 233, 1, p),
                    ur::baseline::MKCV_ADF_N_233[i], 1e-12);
        EXPECT_NEAR(ur::mackinnon_critical_value("adf", "n", 100, 1, p),
                    ur::baseline::MKCV_ADF_N_100[i], 1e-12);
        EXPECT_NEAR(ur::mackinnon_critical_value("adf", "nc", 0, 1, p),
                    ur::baseline::MKCV_ADF_N_INF[i], 1e-12);
    }
}

TEST(MackinnonCv, AdfCriticalValueC) {
    for (int i = 0; i < 3; ++i) {
        const Real p = i == 0 ? 0.01 : (i == 1 ? 0.05 : 0.10);
        EXPECT_NEAR(ur::mackinnon_critical_value("adf", "c", 233, 1, p),
                    ur::baseline::MKCV_ADF_C_233[i], 1e-12);
        EXPECT_NEAR(ur::mackinnon_critical_value("adf", "c", 100, 1, p),
                    ur::baseline::MKCV_ADF_C_100[i], 1e-12);
        EXPECT_NEAR(ur::mackinnon_critical_value("adf", "c", 0, 1, p),
                    ur::baseline::MKCV_ADF_C_INF[i], 1e-12);
    }
}

TEST(MackinnonCv, AdfCriticalValueCt) {
    for (int i = 0; i < 3; ++i) {
        const Real p = i == 0 ? 0.01 : (i == 1 ? 0.05 : 0.10);
        EXPECT_NEAR(ur::mackinnon_critical_value("adf", "ct", 233, 1, p),
                    ur::baseline::MKCV_ADF_CT_233[i], 1e-12);
        EXPECT_NEAR(ur::mackinnon_critical_value("adf", "ct", 100, 1, p),
                    ur::baseline::MKCV_ADF_CT_100[i], 1e-12);
        EXPECT_NEAR(ur::mackinnon_critical_value("adf", "ct", 0, 1, p),
                    ur::baseline::MKCV_ADF_CT_INF[i], 1e-12);
    }
}

// ---------------------------------------------------------------------------
// 4. PP 共享 ADF 表 (同为 DF 分布, U13)
// ---------------------------------------------------------------------------
TEST(MackinnonCv, PpSharesAdfTable) {
    for (const std::string tr : {"nc", "c", "ct"}) {
        for (Size T : {233u, 100u, 0u}) {
            for (Real p : {0.01, 0.05, 0.10}) {
                EXPECT_DOUBLE_EQ(
                    ur::mackinnon_critical_value("pp", tr, T, 1, p),
                    ur::mackinnon_critical_value("adf", tr, T, 1, p));
                EXPECT_DOUBLE_EQ(
                    ur::mackinnon_p_value(-2.87, "pp", tr, T, 1),
                    ur::mackinnon_p_value(-2.87, "adf", tr, T, 1));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 5. DF-GLS 临界值 (arch 独立模拟表, U7): c/ct × 3 T × 3 水平 (1e-12)
// ---------------------------------------------------------------------------
TEST(MackinnonCv, DfglsCriticalValues) {
    for (int i = 0; i < 3; ++i) {
        const Real p = i == 0 ? 0.01 : (i == 1 ? 0.05 : 0.10);
        EXPECT_NEAR(ur::mackinnon_critical_value("df_gls", "c", 233, 1, p),
                    ur::baseline::MKCV_DFGLS_C_233[i], 1e-12);
        EXPECT_NEAR(ur::mackinnon_critical_value("df_gls", "c", 100, 1, p),
                    ur::baseline::MKCV_DFGLS_C_100[i], 1e-12);
        EXPECT_NEAR(ur::mackinnon_critical_value("df_gls", "c", 0, 1, p),
                    ur::baseline::MKCV_DFGLS_C_INF[i], 1e-12);
        EXPECT_NEAR(ur::mackinnon_critical_value("df_gls", "ct", 233, 1, p),
                    ur::baseline::MKCV_DFGLS_CT_233[i], 1e-12);
        EXPECT_NEAR(ur::mackinnon_critical_value("df_gls", "ct", 100, 1, p),
                    ur::baseline::MKCV_DFGLS_CT_100[i], 1e-12);
        EXPECT_NEAR(ur::mackinnon_critical_value("df_gls", "ct", 0, 1, p),
                    ur::baseline::MKCV_DFGLS_CT_INF[i], 1e-12);
    }
}

// ---------------------------------------------------------------------------
// 6. 非法输入抛异常 (p/trend/test_type/n_params)
// ---------------------------------------------------------------------------
TEST(MackinnonCv, ThrowsOnInvalidInputs) {
    EXPECT_THROW(ur::mackinnon_critical_value("adf", "c", 100, 1, 0.025),
                 std::invalid_argument);
    EXPECT_THROW(ur::mackinnon_critical_value("adf", "ctt", 100, 1, 0.05),
                 std::invalid_argument);
    EXPECT_THROW(ur::mackinnon_critical_value("kpss", "c", 100, 1, 0.05),
                 std::invalid_argument);
    EXPECT_THROW(ur::mackinnon_critical_value("df_gls", "nc", 100, 1, 0.05),
                 std::invalid_argument);
    EXPECT_THROW(ur::mackinnon_critical_value("adf", "c", 100, 2, 0.05),
                 std::invalid_argument);
    EXPECT_THROW(ur::mackinnon_p_value(-2.0, "adf", "bad", 100, 1),
                 std::invalid_argument);
    EXPECT_THROW(ur::mackinnon_p_value(-2.0, "df_gls", "nc", 100, 1),
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 7. ADF p 值 smallp 段: t={-3.5, -2.5} (< tau_star=-1.61) (1e-12)
// ---------------------------------------------------------------------------
TEST(MackinnonCv, AdfPValueSmallPRegion) {
    EXPECT_NEAR(ur::mackinnon_p_value(-3.5, "adf", "c", 233, 1),
                ur::baseline::MKP_ADF_C_T233_M3p5, 1e-12);
    EXPECT_NEAR(ur::mackinnon_p_value(-2.5, "adf", "c", 233, 1),
                ur::baseline::MKP_ADF_C_T233_M2p5, 1e-12);
}

// ---------------------------------------------------------------------------
// 8. ADF p 值 largep 段: t=-1 (> tau_star) — 升幂方向性验证 (降幂给 1-p)
// ---------------------------------------------------------------------------
TEST(MackinnonCv, AdfPValueLargePRegion) {
    EXPECT_NEAR(ur::mackinnon_p_value(-1.0, "adf", "c", 233, 1),
                ur::baseline::MKP_ADF_C_T233_M1, 1e-12);
}

// ---------------------------------------------------------------------------
// 9. DF-GLS p 值 (arch 独立系数, U7): t={-3.5, -2.0} (1e-12)
// ---------------------------------------------------------------------------
TEST(MackinnonCv, DfglsPValue) {
    EXPECT_NEAR(ur::mackinnon_p_value(-3.5, "df_gls", "c", 233, 1),
                ur::baseline::MKP_DFGLS_C_M3p5, 1e-12);
    EXPECT_NEAR(ur::mackinnon_p_value(-2.0, "df_gls", "c", 233, 1),
                ur::baseline::MKP_DFGLS_C_M2, 1e-12);
}

// ---------------------------------------------------------------------------
// 10. p 值越界 clamp: stat < tau_min → 0; stat > tau_max → 1
// ---------------------------------------------------------------------------
TEST(MackinnonCv, PValueClampsOutsideSupport) {
    EXPECT_DOUBLE_EQ(ur::mackinnon_p_value(-25.0, "adf", "c", 233, 1), 0.0);
    EXPECT_DOUBLE_EQ(ur::mackinnon_p_value(5.0, "adf", "c", 233, 1), 1.0);
}

// ---------------------------------------------------------------------------
// 11. KPSS p 值线性插值 + 越界 clamp (U6)
// ---------------------------------------------------------------------------
TEST(MackinnonCv, KpssPValueInterpolation) {
    EXPECT_NEAR(ur::kpss_p_value(0.5, "c"),
                ur::baseline::KPSSCRIT_C_P05, 1e-14);
    EXPECT_NEAR(ur::kpss_p_value(0.3, "ct"),
                ur::baseline::KPSSCRIT_CT_P03, 1e-14);
    // 超上界 → clamp 到 y_min/100 = 0.01/100 = 0.0001
    EXPECT_NEAR(ur::kpss_p_value(10.0, "c"),
                ur::baseline::KPSSCRIT_C_P10, 1e-14);
    EXPECT_THROW(ur::kpss_p_value(0.5, "nc"), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 12. KPSS 临界值 {1%, 5%, 10%}: 表内精确匹配点 (1e-14)
// ---------------------------------------------------------------------------
TEST(MackinnonCv, KpssCriticalValues) {
    const auto cv_c = ur::kpss_critical_values("c");
    ASSERT_EQ(cv_c.size(), 3u);
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(cv_c[i], ur::baseline::KPSSCRIT_C_CV[i], 1e-14);
    }
    const auto cv_ct = ur::kpss_critical_values("ct");
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(cv_ct[i], ur::baseline::KPSSCRIT_CT_CV[i], 1e-14);
    }
    EXPECT_THROW(ur::kpss_critical_values("n"), std::invalid_argument);
}
