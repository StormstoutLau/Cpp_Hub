// =============================================================================
// test_johansen_test.cpp - Johansen 协整秩检验 (Phase 7C M3, 18 用例)
//
// 基准: coint_baseline.inc (statsmodels 0.14.4 1e-10 主基准 +
//       urca 1.3-4 det0↔none 映射交叉 1e-8, JOHANSEN_DUAL_LIB_DIFF.md §2)
// 幻觉点: CI4 (3 情形) / CI5 (双表源) / CI6 (迹公式+λ降序+有效T) /
//         CI7 (spec 恒等, 报告 §3) / CI8 (evec 列符号任意)
// =============================================================================

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/timeseries/cointegration/johansen_test.hpp"
#include "cpphub/timeseries/cointegration/osterwald_lenum_cv.hpp"
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

vt::MultivariateTSData fixture2() {  // y1, y3 (无协整对)
    vt::MultivariateTSData d;
    d.columns = {std::vector<Real>(std::begin(cb::Y1), std::end(cb::Y1)),
                 std::vector<Real>(std::begin(cb::Y3), std::end(cb::Y3))};
    d.names = {"y1", "y3"};
    return d;
}

const double* jo_arr(int det, int k, const char* what) {
    // JO_DET<d>_K<k>_<WHAT>
    if (det == -1 && k == 1) {
        if (std::string(what) == "EIG") return cb::JO_DETM1_K1_EIG;
        if (std::string(what) == "LR1") return cb::JO_DETM1_K1_LR1;
        if (std::string(what) == "LR2") return cb::JO_DETM1_K1_LR2;
        if (std::string(what) == "CVT") return cb::JO_DETM1_K1_CVT;
        if (std::string(what) == "CVM") return cb::JO_DETM1_K1_CVM;
        if (std::string(what) == "EVEC") return cb::JO_DETM1_K1_EVEC;
    }
    if (det == -1 && k == 2) {
        if (std::string(what) == "EIG") return cb::JO_DETM1_K2_EIG;
        if (std::string(what) == "LR1") return cb::JO_DETM1_K2_LR1;
        if (std::string(what) == "LR2") return cb::JO_DETM1_K2_LR2;
        if (std::string(what) == "CVT") return cb::JO_DETM1_K2_CVT;
        if (std::string(what) == "CVM") return cb::JO_DETM1_K2_CVM;
        if (std::string(what) == "EVEC") return cb::JO_DETM1_K2_EVEC;
    }
    if (det == 0 && k == 1) {
        if (std::string(what) == "EIG") return cb::JO_DET0_K1_EIG;
        if (std::string(what) == "LR1") return cb::JO_DET0_K1_LR1;
        if (std::string(what) == "LR2") return cb::JO_DET0_K1_LR2;
        if (std::string(what) == "CVT") return cb::JO_DET0_K1_CVT;
        if (std::string(what) == "CVM") return cb::JO_DET0_K1_CVM;
        if (std::string(what) == "EVEC") return cb::JO_DET0_K1_EVEC;
    }
    if (det == 0 && k == 2) {
        if (std::string(what) == "EIG") return cb::JO_DET0_K2_EIG;
        if (std::string(what) == "LR1") return cb::JO_DET0_K2_LR1;
        if (std::string(what) == "LR2") return cb::JO_DET0_K2_LR2;
        if (std::string(what) == "CVT") return cb::JO_DET0_K2_CVT;
        if (std::string(what) == "CVM") return cb::JO_DET0_K2_CVM;
        if (std::string(what) == "EVEC") return cb::JO_DET0_K2_EVEC;
    }
    if (det == 1 && k == 1) {
        if (std::string(what) == "EIG") return cb::JO_DET1_K1_EIG;
        if (std::string(what) == "LR1") return cb::JO_DET1_K1_LR1;
        if (std::string(what) == "LR2") return cb::JO_DET1_K1_LR2;
        if (std::string(what) == "CVT") return cb::JO_DET1_K1_CVT;
        if (std::string(what) == "CVM") return cb::JO_DET1_K1_CVM;
        if (std::string(what) == "EVEC") return cb::JO_DET1_K1_EVEC;
    }
    if (det == 1 && k == 2) {
        if (std::string(what) == "EIG") return cb::JO_DET1_K2_EIG;
        if (std::string(what) == "LR1") return cb::JO_DET1_K2_LR1;
        if (std::string(what) == "LR2") return cb::JO_DET1_K2_LR2;
        if (std::string(what) == "CVT") return cb::JO_DET1_K2_CVT;
        if (std::string(what) == "CVM") return cb::JO_DET1_K2_CVM;
        if (std::string(what) == "EVEC") return cb::JO_DET1_K2_EVEC;
    }
    ADD_FAILURE() << "unknown jo baseline";
    return nullptr;
}

}  // namespace

// 1. 特征值 λ 降序 vs statsmodels (CI6), 全 det × k
TEST(JohansenTest, EigenvaluesVsStatsmodels) {
    for (int det : {-1, 0, 1}) {
        for (int k : {1, 2}) {
            const auto r = co::coint_johansen(fixture3(), det, k);
            const double* ref = jo_arr(det, k, "EIG");
            for (int i = 0; i < 3; ++i) {
                EXPECT_NEAR(r.eig(i), ref[i], 1e-10)
                    << "det=" << det << " k=" << k << " i=" << i;
            }
            // λ 降序 (CI6)
            EXPECT_GE(r.eig(0), r.eig(1));
            EXPECT_GE(r.eig(1), r.eig(2));
        }
    }
}

// 2. 迹统计量 vs statsmodels (CI6: −T·Σ_{i>r} ln(1−λᵢ))
TEST(JohansenTest, TraceStatVsStatsmodels) {
    for (int det : {-1, 0, 1}) {
        for (int k : {1, 2}) {
            const auto r = co::coint_johansen(fixture3(), det, k);
            const double* ref = jo_arr(det, k, "LR1");
            for (int i = 0; i < 3; ++i) {
                EXPECT_NEAR(r.lr1(i), ref[i], 1e-10)
                    << "det=" << det << " k=" << k;
            }
        }
    }
}

// 3. 最大特征值统计量 vs statsmodels
TEST(JohansenTest, MaxEigStatVsStatsmodels) {
    for (int det : {-1, 0, 1}) {
        for (int k : {1, 2}) {
            const auto r = co::coint_johansen(fixture3(), det, k);
            const double* ref = jo_arr(det, k, "LR2");
            for (int i = 0; i < 3; ++i) {
                EXPECT_NEAR(r.lr2(i), ref[i], 1e-10)
                    << "det=" << det << " k=" << k;
            }
        }
    }
}

// 4. cvt/cvm vs statsmodels (MHM96 表, 90/95/99)
TEST(JohansenTest, CriticalValuesVsStatsmodels) {
    for (int det : {-1, 0, 1}) {
        for (int k : {1, 2}) {
            const auto r = co::coint_johansen(fixture3(), det, k);
            const double* cvt = jo_arr(det, k, "CVT");
            const double* cvm = jo_arr(det, k, "CVM");
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    EXPECT_NEAR(r.cvt(i, j), cvt[3 * i + j], 1e-10);
                    EXPECT_NEAR(r.cvm(i, j), cvm[3 * i + j], 1e-10);
                }
            }
            EXPECT_EQ(r.cv_source, "MHM96");
        }
    }
}

// 5. evec vs statsmodels (CI8: 列符号任意 — 逐列符号对齐后比较)
TEST(JohansenTest, EigenvectorsColumnAligned) {
    for (int det : {-1, 0, 1}) {
        for (int k : {1, 2}) {
            const auto r = co::coint_johansen(fixture3(), det, k);
            const double* ref = jo_arr(det, k, "EVEC");
            for (int col = 0; col < 3; ++col) {
                // 符号对齐: 用最大绝对值元素
                int big = 0;
                for (int i = 1; i < 3; ++i) {
                    if (std::fabs(r.evec(i, col)) >
                        std::fabs(r.evec(big, col))) {
                        big = i;
                    }
                }
                const double s =
                    (r.evec(big, col) * ref[3 * big + col] < 0) ? -1.0 : 1.0;
                for (int i = 0; i < 3; ++i) {
                    EXPECT_NEAR(s * r.evec(i, col), ref[3 * i + col], 1e-10)
                        << "det=" << det << " k=" << k << " col=" << col;
                }
            }
        }
    }
}

// 6. 迹/最大特征值恒等式: lr1[r] − lr1[r+1] == lr2[r] (1e-12 自洽)
TEST(JohansenTest, TraceMaxEigIdentity) {
    const auto r = co::coint_johansen(fixture3(), 0, 1);
    EXPECT_NEAR(r.lr1(0) - r.lr1(1), r.lr2(0), 1e-12);
    EXPECT_NEAR(r.lr1(1) - r.lr1(2), r.lr2(1), 1e-12);
    // lr2[r] = −T·ln(1−λ_{r+1}) 手算
    const Real T = static_cast<Real>(r.n_obs);
    EXPECT_NEAR(r.lr2(0), -T * std::log(1.0 - r.eig(0)), 1e-10);
}

// 7. 有效样本 T−1−k (CI6)
TEST(JohansenTest, EffectiveSampleSize) {
    EXPECT_EQ(co::coint_johansen(fixture3(), 0, 1).n_obs, 248u);
    EXPECT_EQ(co::coint_johansen(fixture3(), 0, 2).n_obs, 247u);
    EXPECT_EQ(co::coint_johansen(fixture3(), -1, 0).n_obs, 249u);
}

// 8. rank 选择 (trace) vs statsmodels select_coint_rank 全组合
TEST(JohansenTest, RankSelectionTrace) {
    struct Case { int det; int k; };
    const double sigs[3] = {0.1, 0.05, 0.01};
    for (const Case& c : {Case{-1, 1}, Case{-1, 2}, Case{0, 1}, Case{0, 2},
                          Case{1, 1}, Case{1, 2}}) {
        for (double s : sigs) {
            const Size rk = co::select_coint_rank(fixture3(), c.det, c.k,
                                                  "trace", s);
            // 基准命名 JO_DET<d>_K<k>_RANK_TRACE_<sig>
            double ref;
            if (c.det == -1 && c.k == 1)
                ref = (s == 0.1) ? cb::JO_DETM1_K1_RANK_TRACE_0_1[0]
                      : (s == 0.05) ? cb::JO_DETM1_K1_RANK_TRACE_0_05[0]
                                    : cb::JO_DETM1_K1_RANK_TRACE_0_01[0];
            else if (c.det == -1 && c.k == 2)
                ref = (s == 0.1) ? cb::JO_DETM1_K2_RANK_TRACE_0_1[0]
                      : (s == 0.05) ? cb::JO_DETM1_K2_RANK_TRACE_0_05[0]
                                    : cb::JO_DETM1_K2_RANK_TRACE_0_01[0];
            else if (c.det == 0 && c.k == 1)
                ref = (s == 0.1) ? cb::JO_DET0_K1_RANK_TRACE_0_1[0]
                      : (s == 0.05) ? cb::JO_DET0_K1_RANK_TRACE_0_05[0]
                                    : cb::JO_DET0_K1_RANK_TRACE_0_01[0];
            else if (c.det == 0 && c.k == 2)
                ref = (s == 0.1) ? cb::JO_DET0_K2_RANK_TRACE_0_1[0]
                      : (s == 0.05) ? cb::JO_DET0_K2_RANK_TRACE_0_05[0]
                                    : cb::JO_DET0_K2_RANK_TRACE_0_01[0];
            else if (c.det == 1 && c.k == 1)
                ref = (s == 0.1) ? cb::JO_DET1_K1_RANK_TRACE_0_1[0]
                      : (s == 0.05) ? cb::JO_DET1_K1_RANK_TRACE_0_05[0]
                                    : cb::JO_DET1_K1_RANK_TRACE_0_01[0];
            else
                ref = (s == 0.1) ? cb::JO_DET1_K2_RANK_TRACE_0_1[0]
                      : (s == 0.05) ? cb::JO_DET1_K2_RANK_TRACE_0_05[0]
                                    : cb::JO_DET1_K2_RANK_TRACE_0_01[0];
            EXPECT_EQ(rk, static_cast<Size>(ref))
                << "det=" << c.det << " k=" << c.k << " sig=" << s;
        }
    }
}

// 9. rank 选择 (maxeig) vs statsmodels
TEST(JohansenTest, RankSelectionMaxEig) {
    EXPECT_EQ(co::select_coint_rank(fixture3(), 0, 1, "maxeig", 0.05),
              static_cast<Size>(cb::JO_DET0_K1_RANK_MAXEIG_0_05[0]));
    EXPECT_EQ(co::select_coint_rank(fixture3(), 1, 2, "maxeig", 0.1),
              static_cast<Size>(cb::JO_DET1_K2_RANK_MAXEIG_0_1[0]));
    EXPECT_EQ(co::select_coint_rank(fixture3(), -1, 1, "maxeig", 0.01),
              static_cast<Size>(cb::JO_DETM1_K1_RANK_MAXEIG_0_01[0]));
}

// 10. rank 逐级检验语义 (B4: lr[r] < cv[r, signif] 即接受)
TEST(JohansenTest, RankSelectionSemantics) {
    const auto r = co::coint_johansen(fixture3(), 0, 1);
    // r=0: 96.69 > cv 29.80 → 拒绝; r=1: 14.51 < 15.49 → 接受 ⇒ rank 1
    EXPECT_GT(r.lr1(0), r.cvt(0, 1));
    EXPECT_LT(r.lr1(1), r.cvt(1, 1));
    EXPECT_EQ(co::select_coint_rank(fixture3(), 0, 1, "trace", 0.05), 1u);
    // 非协整对 (y1, y3): K=2
    const auto r2 = co::coint_johansen(fixture2(), 0, 1);
    Size expect = r2.lr1(0) < r2.cvt(0, 1) ? 0 : 1;
    EXPECT_EQ(co::select_coint_rank(fixture2(), 0, 1, "trace", 0.05), expect);
}

// 11. OL1992 查表 (urca 转录; CI5 双表)
TEST(JohansenTest, OL1992Lookups) {
    const auto t3 = co::ol1992_trace_cv("none", 3);
    EXPECT_DOUBLE_EQ(t3[0], 28.71);
    EXPECT_DOUBLE_EQ(t3[1], 31.52);
    EXPECT_DOUBLE_EQ(t3[2], 37.22);
    const auto m2 = co::ol1992_maxeig_cv("trend", 2);
    EXPECT_DOUBLE_EQ(m2[1], 18.96);
    EXPECT_THROW(co::ol1992_trace_cv("none", 12), std::invalid_argument);
    EXPECT_THROW(co::ol1992_trace_cv("zz", 1), std::invalid_argument);
}

// 12. MHM96 查表 (statsmodels 转录)
TEST(JohansenTest, MHM96Lookups) {
    const auto t = co::mhm96_trace_cv(0, 3);
    EXPECT_DOUBLE_EQ(t[1], 29.7961);
    const auto m = co::mhm96_maxeig_cv(1, 3);
    EXPECT_DOUBLE_EQ(m[0], 21.8731);
    EXPECT_THROW(co::mhm96_trace_cv(0, 13), std::invalid_argument);
    EXPECT_THROW(co::mhm96_trace_cv(2, 3), std::invalid_argument);
}

// 13. CI5: 双表源差异断言 (OL1992 vs MHM96 名义对应情形数值不同 — 已冻结的
//     双库差异; 不得混用)
TEST(JohansenTest, DualTableSourceDifference) {
    const auto ol = co::ol1992_trace_cv("none", 3);   // urca ecdet=none
    const auto mh = co::mhm96_trace_cv(0, 3);         // SM det_order=0
    // 同一映射情形 (det0 ↔ none) 的两套表: 95% 值 31.52 vs 29.7961 (MC 裁决:
    // OL1992 匹配统计量分布, MHM96 不匹配 — 见 JOHANSEN_DUAL_LIB_DIFF.md §5)
    EXPECT_NE(ol[1], mh[1]);
    EXPECT_NEAR(ol[1], 31.52, 0.0);
    EXPECT_NEAR(mh[1], 29.7961, 0.0);
}

// 14. urca 交叉: det_order=0 ↔ ecdet="none" (k_ar_diff = K−1), 1e-8
TEST(JohansenTest, UrcaCrossCheckDet0) {
    const auto k1 = co::coint_johansen(fixture3(), 0, 1);  // ↔ urca K=2
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(k1.eig(i), cb::UR_K2_NONE_EIG[i], 1e-8);
        EXPECT_NEAR(k1.lr1(i), cb::UR_K2_NONE_LR1[i], 1e-8);
        EXPECT_NEAR(k1.lr2(i), cb::UR_K2_NONE_LR2[i], 1e-8);
    }
    EXPECT_EQ(k1.n_obs, static_cast<Size>(cb::UR_K2_NONE_NOBS[0]));
    const auto k2 = co::coint_johansen(fixture3(), 0, 2);  // ↔ urca K=3
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(k2.eig(i), cb::UR_K3_NONE_EIG[i], 1e-8);
        EXPECT_NEAR(k2.lr1(i), cb::UR_K3_NONE_LR1[i], 1e-8);
        EXPECT_NEAR(k2.lr2(i), cb::UR_K3_NONE_LR2[i], 1e-8);
    }
}

// 15. 变量重排不变性: RRR 谱对变量置换不变 (eig/lr1/lr2 逐值相同; evec 行同步置换)
TEST(JohansenTest, ReorderInvariance) {
    const auto base = co::coint_johansen(fixture3(), 0, 1);
    const auto perm = co::coint_johansen(fixture3().reorder({2, 0, 1}), 0, 1);
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(base.eig(i), perm.eig(i), 1e-12);
        EXPECT_NEAR(base.lr1(i), perm.lr1(i), 1e-10);
        EXPECT_NEAR(base.lr2(i), perm.lr2(i), 1e-10);
    }
    // evec 行按置换同步: reorder({2,0,1}) ⇒ perm 行 j = base 行 perm_map[j]
    const int perm_map[3] = {2, 0, 1};
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            // 符号可能逐列不同 — 用绝对值比较 (列符号任意, CI8)
            EXPECT_NEAR(std::fabs(perm.evec(row, col)),
                        std::fabs(base.evec(perm_map[row], col)), 1e-10);
        }
    }
}

// 16. 输入校验
TEST(JohansenTest, InputValidation) {
    EXPECT_THROW(co::coint_johansen(fixture3(), 2, 1), std::invalid_argument);
    EXPECT_THROW(co::coint_johansen(fixture3(), -2, 1), std::invalid_argument);
    EXPECT_THROW(co::coint_johansen(fixture3(), 0, 400), std::invalid_argument);
    vt::MultivariateTSData one;
    one.columns = {std::vector<Real>(std::begin(cb::Y1),
                                     std::begin(cb::Y1) + 50)};
    EXPECT_THROW(co::coint_johansen(one, 0, 1), std::invalid_argument);  // K<2
    EXPECT_THROW(co::select_coint_rank(fixture3(), 0, 1, "zz", 0.05),
                 std::invalid_argument);
    EXPECT_THROW(co::select_coint_rank(fixture3(), 0, 1, "trace", 0.025),
                 std::invalid_argument);
}

// 17. 结果字段回显 (det_order/k_ar_diff/n_obs/维度)
TEST(JohansenTest, ResultFieldsEcho) {
    const auto r = co::coint_johansen(fixture3(), 1, 2);
    EXPECT_EQ(r.det_order, 1);
    EXPECT_EQ(r.k_ar_diff, 2u);
    EXPECT_EQ(r.n_obs, 247u);
    EXPECT_EQ(r.eig.size(), 3);
    EXPECT_EQ(r.cvt.rows(), 3);
    EXPECT_EQ(r.cvt.cols(), 3);
}

// 18. 性能预算 (P6: Johansen N=3, T=500 < 2s 单线程)
TEST(JohansenTest, PerformanceBudget) {
    vt::MultivariateTSData big;
    big.columns = {std::vector<Real>(), std::vector<Real>(),
                   std::vector<Real>()};
    for (int rep = 0; rep < 2; ++rep) {  // T=500
        big.columns[0].insert(big.columns[0].end(), std::begin(cb::Y1),
                              std::end(cb::Y1));
        big.columns[1].insert(big.columns[1].end(), std::begin(cb::Y2),
                              std::end(cb::Y2));
        big.columns[2].insert(big.columns[2].end(), std::begin(cb::Y3),
                              std::end(cb::Y3));
    }
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = co::coint_johansen(big, 0, 2);
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - t0)
                        .count();
    EXPECT_EQ(r.n_obs, 497u);
    EXPECT_LT(ms, 2000);  // P6 预算 (单线程)
}
