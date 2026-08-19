// =============================================================================
// test_vecm_model.cpp - VECM 估计 (Phase 7C M3, 16 用例)
//
// 基准: coint_baseline.inc (statsmodels 0.14.4 VECM ML, 1e-10)
// 幻觉点: CI4 (5 情形) / CI8 (β 双归一 + 投影空间对照) / CI9 (Π=αβ') /
//         CI10 (ECT t: EM2002 表, 非标准 t)
// =============================================================================

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/timeseries/cointegration/ericsson_mackinnon_cv.hpp"
#include "cpphub/timeseries/cointegration/vecm_model.hpp"
#include "cpphub/timeseries/var/multivariate_data.hpp"
#include "coint_baseline.inc"

namespace cb = cpphub::v1::timeseries::coint_baseline;
namespace co = cpphub::v1::timeseries::cointegration;
namespace vt = cpphub::v1::timeseries::var;
using cpphub::v1::Real;
using cpphub::v1::Size;
using Eigen::MatrixXd;

namespace {

vt::MultivariateTSData fixture3() {
    vt::MultivariateTSData d;
    d.columns = {std::vector<Real>(std::begin(cb::Y1), std::end(cb::Y1)),
                 std::vector<Real>(std::begin(cb::Y2), std::end(cb::Y2)),
                 std::vector<Real>(std::begin(cb::Y3), std::end(cb::Y3))};
    d.names = {"y1", "y2", "y3"};
    return d;
}

// 基准数组选取: VECM_<DET>_R<r>_K<k>_<WHAT>
const double* va(const char* det, int r, int k, const char* what) {
    const std::string key =
        std::string(det) + "_r" + std::to_string(r) + "_k" +
        std::to_string(k) + "_" + what;
    if (key == "n_r1_k1_alpha") return cb::VECM_N_R1_K1_ALPHA;
    if (key == "n_r1_k1_stderr_alpha") return cb::VECM_N_R1_K1_STDERR_ALPHA;
    if (key == "n_r1_k1_beta") return cb::VECM_N_R1_K1_BETA;
    if (key == "n_r1_k1_gamma") return cb::VECM_N_R1_K1_GAMMA;
    if (key == "n_r1_k1_sigma_u") return cb::VECM_N_R1_K1_SIGMA_U;
    if (key == "n_r1_k1_llf") return cb::VECM_N_R1_K1_LLF;
    if (key == "n_r1_k1_resid_head") return cb::VECM_N_R1_K1_RESID_HEAD;
    if (key == "n_r1_k1_nobs") return cb::VECM_N_R1_K1_NOBS;
    if (key == "co_r1_k1_alpha") return cb::VECM_CO_R1_K1_ALPHA;
    if (key == "co_r1_k1_beta") return cb::VECM_CO_R1_K1_BETA;
    if (key == "co_r1_k1_gamma") return cb::VECM_CO_R1_K1_GAMMA;
    if (key == "co_r1_k1_det_coef") return cb::VECM_CO_R1_K1_DET_COEF;
    if (key == "co_r1_k1_sigma_u") return cb::VECM_CO_R1_K1_SIGMA_U;
    if (key == "co_r1_k1_llf") return cb::VECM_CO_R1_K1_LLF;
    if (key == "co_r1_k1_stderr_alpha") return cb::VECM_CO_R1_K1_STDERR_ALPHA;
    if (key == "ci_r1_k1_alpha") return cb::VECM_CI_R1_K1_ALPHA;
    if (key == "ci_r1_k1_beta") return cb::VECM_CI_R1_K1_BETA;
    if (key == "ci_r1_k1_gamma") return cb::VECM_CI_R1_K1_GAMMA;
    if (key == "ci_r1_k1_det_coef_coint") return cb::VECM_CI_R1_K1_DET_COEF_COINT;
    if (key == "ci_r1_k1_llf") return cb::VECM_CI_R1_K1_LLF;
    if (key == "lo_r1_k1_alpha") return cb::VECM_LO_R1_K1_ALPHA;
    if (key == "lo_r1_k1_beta") return cb::VECM_LO_R1_K1_BETA;
    if (key == "lo_r1_k1_gamma") return cb::VECM_LO_R1_K1_GAMMA;
    if (key == "lo_r1_k1_det_coef") return cb::VECM_LO_R1_K1_DET_COEF;
    if (key == "lo_r1_k1_llf") return cb::VECM_LO_R1_K1_LLF;
    if (key == "li_r1_k1_alpha") return cb::VECM_LI_R1_K1_ALPHA;
    if (key == "li_r1_k1_beta") return cb::VECM_LI_R1_K1_BETA;
    if (key == "li_r1_k1_gamma") return cb::VECM_LI_R1_K1_GAMMA;
    if (key == "li_r1_k1_det_coef_coint") return cb::VECM_LI_R1_K1_DET_COEF_COINT;
    if (key == "li_r1_k1_llf") return cb::VECM_LI_R1_K1_LLF;
    if (key == "n_r1_k2_alpha") return cb::VECM_N_R1_K2_ALPHA;
    if (key == "n_r1_k2_beta") return cb::VECM_N_R1_K2_BETA;
    if (key == "n_r1_k2_gamma") return cb::VECM_N_R1_K2_GAMMA;
    if (key == "n_r1_k2_llf") return cb::VECM_N_R1_K2_LLF;
    if (key == "n_r2_k1_alpha") return cb::VECM_N_R2_K1_ALPHA;
    if (key == "n_r2_k1_beta") return cb::VECM_N_R2_K1_BETA;
    if (key == "n_r2_k1_gamma") return cb::VECM_N_R2_K1_GAMMA;
    if (key == "n_r2_k1_llf") return cb::VECM_N_R2_K1_LLF;
    if (key == "n_r2_k1_stderr_alpha") return cb::VECM_N_R2_K1_STDERR_ALPHA;
    if (key == "ci_r1_k1_sigma_u") return cb::VECM_CI_R1_K1_SIGMA_U;
    ADD_FAILURE() << "unknown vecm baseline: " << key;
    return nullptr;
}

}  // namespace

// 1. alpha vs statsmodels — 5 情形 × r=1 × k=1 (CI4), 1e-10
TEST(VecmModel, AlphaAllDetCases) {
    const char* dets[] = {"n", "co", "ci", "lo", "li"};
    for (const char* det : dets) {
        const auto res = co::vecm_fit(fixture3(), 1, 1, det);
        const double* ref = va(det, 1, 1, "alpha");
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(res.alpha(i, 0), ref[i], 1e-10) << "det=" << det;
        }
        EXPECT_EQ(res.det, det);
        EXPECT_EQ(res.rank, 1u);
        EXPECT_EQ(res.n_obs, 248u);
    }
}

// 2. beta vs statsmodels (r=1: 前 1 行 = I → β = [1, b₂, b₃]'; 归一确定性)
TEST(VecmModel, BetaAllDetCases) {
    const char* dets[] = {"n", "co", "ci", "lo", "li"};
    for (const char* det : dets) {
        const auto res = co::vecm_fit(fixture3(), 1, 1, det);
        const double* ref = va(det, 1, 1, "beta");
        EXPECT_DOUBLE_EQ(res.beta(0, 0), 1.0);  // Phillips 归一 (决策 21)
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(res.beta(i, 0), ref[i], 1e-10) << "det=" << det;
        }
    }
}

// 3. gamma vs statsmodels (K×K·k 行主序), 1e-10
TEST(VecmModel, GammaAllDetCases) {
    const char* dets[] = {"n", "co", "ci", "lo", "li"};
    for (const char* det : dets) {
        const auto res = co::vecm_fit(fixture3(), 1, 1, det);
        const double* ref = va(det, 1, 1, "gamma");
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                EXPECT_NEAR(res.gamma_flat(i, j), ref[3 * i + j], 1e-10)
                    << "det=" << det;
            }
        }
        // 逐滞后块 = flat 块
        ASSERT_EQ(res.gamma.size(), 1u);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                EXPECT_EQ(res.gamma[0](i, j), res.gamma_flat(i, j));
            }
        }
    }
}

// 4. det_coef / det_coef_coint vs statsmodels (co/lo 外部, ci/li 内部)
TEST(VecmModel, DeterministicCoefficients) {
    {
        const auto res = co::vecm_fit(fixture3(), 1, 1, "co");
        const double* ref = va("co", 1, 1, "det_coef");
        ASSERT_EQ(res.det_coef.rows(), 3);
        ASSERT_EQ(res.det_coef.cols(), 1);
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(res.det_coef(i, 0), ref[i], 1e-10);
        }
        EXPECT_EQ(res.det_coef_coint.size(), 0);
    }
    {
        const auto res = co::vecm_fit(fixture3(), 1, 1, "ci");
        const double* ref = va("ci", 1, 1, "det_coef_coint");
        ASSERT_EQ(res.det_coef_coint.rows(), 1);
        ASSERT_EQ(res.det_coef_coint.cols(), 1);
        EXPECT_NEAR(res.det_coef_coint(0, 0), ref[0], 1e-10);
        EXPECT_EQ(res.det_coef.size(), 0);
    }
    {
        const auto res = co::vecm_fit(fixture3(), 1, 1, "lo");
        const double* ref = va("lo", 1, 1, "det_coef");
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(res.det_coef(i, 0), ref[i], 1e-10);
        }
    }
    {
        const auto res = co::vecm_fit(fixture3(), 1, 1, "li");
        const double* ref = va("li", 1, 1, "det_coef_coint");
        EXPECT_NEAR(res.det_coef_coint(0, 0), ref[0], 1e-10);
    }
}

// 5. sigma_u vs statsmodels (ML: ÷T), 1e-10
TEST(VecmModel, SigmaUVsStatsmodels) {
    for (const char* det : {"n", "co", "ci"}) {
        const auto res = co::vecm_fit(fixture3(), 1, 1, det);
        const double* ref = va(det, 1, 1, "sigma_u");
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                EXPECT_NEAR(res.sigma_u(i, j), ref[3 * i + j], 1e-10)
                    << "det=" << det;
            }
        }
    }
}

// 6. llf vs statsmodels (Lütkepohl 7.2.20), 1e-10
TEST(VecmModel, LogLikelihoodAllCases) {
    const char* dets[] = {"n", "co", "ci", "lo", "li"};
    for (const char* det : dets) {
        const auto res = co::vecm_fit(fixture3(), 1, 1, det);
        EXPECT_NEAR(res.loglik, va(det, 1, 1, "llf")[0], 1e-10)
            << "det=" << det;
    }
}

// 7. 残差 vs statsmodels (前 4 个时点)
//    基准 resid_head 为 statsmodels resid (T×K) 前 4 行的行主序:
//    ref[3j+i] = 时间 j 变量 i; 我们的 resid 为 K×T_eff
TEST(VecmModel, ResidualsVsStatsmodels) {
    const auto res = co::vecm_fit(fixture3(), 1, 1, "n");
    const double* ref = va("n", 1, 1, "resid_head");
    ASSERT_EQ(res.resid.cols(), 248);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_NEAR(res.resid(i, j), ref[3 * j + i], 1e-10);
        }
    }
}

// 8. k_ar_diff=2 vs statsmodels
TEST(VecmModel, KarDiff2) {
    const auto res = co::vecm_fit(fixture3(), 1, 2, "n");
    const double* a = va("n", 1, 2, "alpha");
    const double* b = va("n", 1, 2, "beta");
    const double* g = va("n", 1, 2, "gamma");
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(res.alpha(i, 0), a[i], 1e-10);
        EXPECT_NEAR(res.beta(i, 0), b[i], 1e-10);
        for (int j = 0; j < 6; ++j) {
            EXPECT_NEAR(res.gamma_flat(i, j), g[6 * i + j], 1e-10);
        }
    }
    EXPECT_NEAR(res.loglik, va("n", 1, 2, "llf")[0], 1e-10);
    ASSERT_EQ(res.gamma.size(), 2u);
    EXPECT_EQ(res.n_obs, 247u);
}

// 9. rank=2 vs statsmodels (β 前 2 行 = I₂)
TEST(VecmModel, Rank2Case) {
    const auto res = co::vecm_fit(fixture3(), 2, 1, "n");
    ASSERT_EQ(res.beta.rows(), 3);
    ASSERT_EQ(res.beta.cols(), 2);
    EXPECT_DOUBLE_EQ(res.beta(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(res.beta(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(res.beta(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(res.beta(1, 1), 1.0);
    const double* a = va("n", 2, 1, "alpha");
    const double* b = va("n", 2, 1, "beta");
    const double* g = va("n", 2, 1, "gamma");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            EXPECT_NEAR(res.alpha(i, j), a[2 * i + j], 1e-10);
        }
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(res.gamma_flat(i, j), g[3 * i + j], 1e-10);
        }
    }
    for (int i = 0; i < 3; ++i) {
        // β 第 2 列: K×r 行主序 → b[2i+1] (b[2i] 为第 1 列)
        EXPECT_NEAR(res.beta(i, 1), b[2 * i + 1], 1e-10);
    }
    EXPECT_NEAR(res.loglik, va("n", 2, 1, "llf")[0], 1e-10);
}

// 10. CI8: β 投影空间对照 (P = β(β'β)⁻¹β', r=2)
TEST(VecmModel, BetaProjectionSpace) {
    const auto res = co::vecm_fit(fixture3(), 2, 1, "n");
    const double* bref = va("n", 2, 1, "beta");
    MatrixXd b_ref(3, 2);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) b_ref(i, j) = bref[2 * i + j];
    }
    auto proj = [](const MatrixXd& b) {
        return b * (b.transpose() * b).inverse() * b.transpose();
    };
    const MatrixXd p1 = proj(res.beta);
    const MatrixXd p2 = proj(b_ref);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(p1(i, j), p2(i, j), 1e-10);
        }
    }
}

// 11. urca 归一开关 (决策 21): β 首变量 = 1; Π = αβ' 不变 (CI9)
TEST(VecmModel, UrcaNormalizationSwitch) {
    const auto dflt = co::vecm_fit(fixture3(), 2, 1, "n");
    const auto ur = co::vecm_fit(fixture3(), 2, 1, "n", true);
    EXPECT_TRUE(ur.urca_normalization);
    // 首变量归一: β[0, :] = 1
    EXPECT_DOUBLE_EQ(ur.beta(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(ur.beta(0, 1), 1.0);
    // Π 不变 (α·β' 逐元素)
    const MatrixXd pi1 = dflt.alpha * dflt.beta.transpose();
    const MatrixXd pi2 = ur.alpha * ur.beta.transpose();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(pi1(i, j), pi2(i, j), 1e-10);  // CI9
        }
    }
    // Γ/Σ 与归一无关
    EXPECT_TRUE(ur.gamma_flat.isApprox(dflt.gamma_flat, 1e-12));
    EXPECT_TRUE(ur.sigma_u.isApprox(dflt.sigma_u, 1e-12));
}

// 12. CI9: Π = αβ' 的秩 = r (SVD 奇异值)
TEST(VecmModel, PiRankEqualsCointRank) {
    const auto r1 = co::vecm_fit(fixture3(), 1, 1, "n");
    Eigen::JacobiSVD<MatrixXd> svd1(r1.alpha * r1.beta.transpose());
    EXPECT_GT(svd1.singularValues()(0), 1e-8);
    EXPECT_LT(svd1.singularValues()(1), 1e-8);  // rank 1
    const auto r2 = co::vecm_fit(fixture3(), 2, 1, "n");
    Eigen::JacobiSVD<MatrixXd> svd2(r2.alpha * r2.beta.transpose());
    EXPECT_GT(svd2.singularValues()(1), 1e-8);
    EXPECT_LT(svd2.singularValues()(2), 1e-8);  // rank 2
}

// 13. ECT t vs statsmodels (t = α_j / stderr_alpha_j; r=1), 1e-10
TEST(VecmModel, EctTStatVsStatsmodels) {
    const auto res = co::vecm_fit(fixture3(), 1, 1, "co");
    const double* a = va("co", 1, 1, "alpha");
    const double* se = va("co", 1, 1, "stderr_alpha");
    ASSERT_TRUE(res.has_ect_t);
    for (int i = 0; i < 3; ++i) {
        const double t_ref = a[i] / se[i];
        EXPECT_NEAR(res.ect_t_stat[i], t_ref, 1e-10) << "eq " << i;
    }
}

// 14. rank≥2: ECT t = NaN + has_ect_t = false (§1.4-5)
TEST(VecmModel, EctTNaNWhenRank2) {
    const auto res = co::vecm_fit(fixture3(), 2, 1, "n");
    EXPECT_FALSE(res.has_ect_t);
    for (Real v : res.ect_t_stat) EXPECT_TRUE(std::isnan(v));
    for (Real v : res.ect_cv_5pct) EXPECT_TRUE(std::isnan(v));
}

// 15. ECT 临界值 = EM2002 查表 (CI10; det→case 映射 n→n, co/ci→c, lo/li→ct)
TEST(VecmModel, EctCriticalValuesEM2002) {
    struct Case { const char* det; const char* em; };
    for (const Case& c : {Case{"n", "n"}, Case{"co", "c"}, Case{"ci", "c"},
                          Case{"lo", "ct"}, Case{"li", "ct"}}) {
        const auto res = co::vecm_fit(fixture3(), 1, 1, c.det);
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(res.ect_cv_5pct[i],
                        co::em2002_ect_critical_value(3, c.em, 5.0, 248), 0.0)
                << "det=" << c.det;
        }
    }
    // EM2002 表锚 (转录保真): c n=1 1% 渐近 −3.4307
    EXPECT_DOUBLE_EQ(co::em2002_ect_critical_value_asymptotic(1, "c", 1.0),
                     -3.4307);
    // 有限样本公式: ctt n=3 1% T=51 → −5.0860 (论文表 7 复现, 三重验证)
    EXPECT_NEAR(co::em2002_ect_critical_value(3, "ctt", 1.0, 51), -5.0860,
                5e-4);
    EXPECT_THROW(co::em2002_ect_critical_value(13, "c", 5.0, 100),
                 std::invalid_argument);
    EXPECT_THROW(co::em2002_ect_critical_value(2, "c", 2.5, 100),
                 std::invalid_argument);
}

// 16. 输入校验
TEST(VecmModel, InputValidation) {
    EXPECT_THROW(co::vecm_fit(fixture3(), 0, 1, "n"), std::invalid_argument);
    EXPECT_THROW(co::vecm_fit(fixture3(), 4, 1, "n"), std::invalid_argument);
    EXPECT_THROW(co::vecm_fit(fixture3(), 1, 1, "colo"),
                 std::invalid_argument);
    EXPECT_THROW(co::vecm_fit(fixture3(), 1, 300, "n"),
                 std::invalid_argument);
}
