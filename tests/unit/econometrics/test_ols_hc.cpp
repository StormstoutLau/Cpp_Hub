// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.2 Day 3-4 任务 1.6-1.7 (排幻觉点 E1/E2/E3)
// OLS 估计器 + HC0-HC5 异方差稳健协方差 (手算解析值对照, 容差 1e-12)
//
// 手算数据集 (N=4, K=2, 含截距; 排幻觉点 E1: 截距由用户显式构造, 不自动添加):
//   X = [1,1; 1,2; 1,3; 1,4],  y = [2; 3; 5; 7]
//   X'X = [4,10; 10,30],  det = 20
//   (X'X)^{-1} = (1/20)[30,-10; -10,4] = [3/2,-1/2; -1/2,1/5]
//   X'y = [17; 51]   (Σy=17, Σxy=1·2+2·3+3·5+4·7=2+6+15+28=51)
//   β = (X'X)^{-1} X'y = [0; 17/10] = [0; 1.7]
//   fitted = X·β = [1.7; 3.4; 5.1; 6.8]
//   residuals u = y - fitted = [0.3; -0.4; -0.1; 0.2]
//   u² = [0.09; 0.16; 0.01; 0.04],  SSR = 0.30
//   ȳ = 4.25,  SST = 14.75
//   R² = 1 - 0.30/14.75 = 289/295
//   adj R² = 1 - (1-R²)·(N-1)/(N-K) = 286/295
//   F = (R²/K)/((1-R²)/(N-K)) = 289/6
//   leverage h = [0.7; 0.3; 0.3; 0.7],  Σh = K = 2
//   Classical: σ²=SSR/(N-K)=0.15, vcov=0.15·(X'X)^{-1}=[9/40,-3/40; -3/40,3/100]
//   HC0 = [7/50, -41/1000; -41/1000, 67/5000]
//   HC1 = N/(N-K)·HC0 = 2·HC0 = [7/25, -41/500; -41/500, 67/2500]  (排幻觉点 E2: 分母 N-K)
//   HC2 = [41/105, -17/140; -17/140, 29/700]                       (排幻觉点 E3: leverage h_i)
//   HC3 = [526/441, -563/1470; -563/1470, 327/2450]
//   HC4/HC5 含无理幂 (0.7^1.4 等), 用独立参考实现 (显式循环) 交叉验证

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/econometrics/estimation/ols.hpp"
#include "cpphub/econometrics/inference/hc_standard_errors.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;

// MSVC 默认不定义 M_PI, 用字面量避免 _USE_MATH_DEFINES 依赖
constexpr Real kTwoPi = 6.2831853071795864769252867665590057683943387987502;

// =============================================================================
// 测试夹具: 构造 N=4, K=2 含截距数据集 + 预计算 (X'X)^{-1}, β, residuals
// =============================================================================
namespace {

// 构造设计矩阵 X (含截距列, 排幻觉点 E1: 显式构造, 不自动添加)
MatrixXD make_X() {
    MatrixXD X(4, 2);
    X(0, 0) = 1.0; X(0, 1) = 1.0;
    X(1, 0) = 1.0; X(1, 1) = 2.0;
    X(2, 0) = 1.0; X(2, 1) = 3.0;
    X(3, 0) = 1.0; X(3, 1) = 4.0;
    return X;
}

VectorXD make_y() {
    VectorXD y(4);
    y(0) = 2.0;
    y(1) = 3.0;
    y(2) = 5.0;
    y(3) = 7.0;
    return y;
}

EconData make_data() {
    return make_cross_section(make_X(), make_y(), {"intercept", "x"}, "y");
}

// (X'X)^{-1} (手算: [3/2,-1/2; -1/2,1/5]) — 通过 inverse_symmetric 计算, 复用受信任路径
MatrixXD make_XtX_inv() {
    MatrixXD X = make_X();
    MatrixXD XtX = MatrixXD(X.eigen().transpose() * X.eigen());
    return inverse_symmetric(XtX);
}

// β (手算: [0; 1.7]) — 通过正规方程计算
VectorXD make_beta() {
    MatrixXD X = make_X();
    VectorXD y = make_y();
    MatrixXD A = make_XtX_inv();
    Eigen::VectorXd Xty = X.eigen().transpose() * y.eigen();
    Eigen::VectorXd beta = A.eigen() * Xty;
    return VectorXD(beta);
}

// residuals (手算: [0.3; -0.4; -0.1; 0.2])
VectorXD make_residuals() {
    MatrixXD X = make_X();
    VectorXD y = make_y();
    VectorXD beta = make_beta();
    Eigen::VectorXd fitted = X.eigen() * beta.eigen();
    Eigen::VectorXd res = y.eigen() - fitted;
    return VectorXD(res);
}

// 独立参考实现: HC4/HC5 协方差 (显式 4 重循环, 与生产代码 Eigen 矩阵乘路径完全独立)
// 公式严格遵循任务描述:
//   HC4: δ_i = min(N·h_i/K, 4),            w_i = u_i²/(1-h_i)^δ_i
//   HC5: δ_i = min(N·h_i/K, max(4, 0.7·N·h_max/K)),  w_i = u_i²/(1-h_i)^δ_i,  h_max=max(h_i)
MatrixXD reference_hc45(const MatrixXD& X, const VectorXD& u, const MatrixXD& A, CovarianceType type) {
    const Size n = X.rows();
    const Size k = X.cols();
    // leverage h_i = Σ_{a,b} X(i,a)·A(a,b)·X(i,b)  (独立于生产的 P=X·A·X' 路径)
    std::vector<Real> h(n, 0.0);
    for (Size i = 0; i < n; ++i) {
        Real s = 0.0;
        for (Size a = 0; a < k; ++a)
            for (Size b = 0; b < k; ++b)
                s += X(i, a) * A(a, b) * X(i, b);
        h[i] = s;
    }
    Real hmax = *std::max_element(h.begin(), h.end());
    // 权重 w_i
    std::vector<Real> w(n, 0.0);
    for (Size i = 0; i < n; ++i) {
        Real u2 = u(i) * u(i);
        Real one_minus_h = 1.0 - h[i];
        Real delta = 0.0;
        if (type == CovarianceType::HC4) {
            delta = std::min(static_cast<Real>(n) * h[i] / static_cast<Real>(k), static_cast<Real>(4));
        } else {  // HC5
            delta = std::min(static_cast<Real>(n) * h[i] / static_cast<Real>(k),
                             std::max(static_cast<Real>(4),
                                      static_cast<Real>(0.7) * static_cast<Real>(n) * hmax / static_cast<Real>(k)));
        }
        w[i] = u2 / std::pow(one_minus_h, delta);
    }
    // meat = X' diag(w) X (显式求和)
    MatrixXD meat(k, k);
    for (Size a = 0; a < k; ++a)
        for (Size b = 0; b < k; ++b) {
            Real s = 0.0;
            for (Size i = 0; i < n; ++i) s += w[i] * X(i, a) * X(i, b);
            meat(a, b) = s;
        }
    // V = A · meat · A (显式 4 重循环, 独立于 Eigen)
    MatrixXD V(k, k);
    for (Size a = 0; a < k; ++a)
        for (Size b = 0; b < k; ++b) {
            Real s = 0.0;
            for (Size p = 0; p < k; ++p)
                for (Size q = 0; q < k; ++q)
                    s += A(a, p) * meat(p, q) * A(q, b);
            V(a, b) = s;
        }
    return V;
}

}  // namespace

// =============================================================================
// 1. OLS 系数 β (手算: [0; 1.7])
// =============================================================================
TEST(OlsHcTest, BetaCoefficients_MatchHandCalc) {
    OLSEstimator ols;
    EstimationResult r = ols.estimate(make_data());
    ASSERT_EQ(r.coefficients.size(), 2u);
    EXPECT_NEAR(r.coefficients(0), 0.0, 1e-12);
    EXPECT_NEAR(r.coefficients(1), 1.7, 1e-12);
    EXPECT_EQ(r.n_obs, 4u);
    EXPECT_EQ(r.n_params, 2u);
    EXPECT_EQ(r.df_residual, 2u);
}

// =============================================================================
// 2. 拟合值 (手算: [1.7; 3.4; 5.1; 6.8])
// =============================================================================
TEST(OlsHcTest, FittedValues_MatchHandCalc) {
    MatrixXD X = make_X();
    VectorXD beta = make_beta();
    VectorXD fitted = OLSEstimator{}.computeFittedValues(X, beta);
    ASSERT_EQ(fitted.size(), 4u);
    EXPECT_NEAR(fitted(0), 1.7, 1e-12);
    EXPECT_NEAR(fitted(1), 3.4, 1e-12);
    EXPECT_NEAR(fitted(2), 5.1, 1e-12);
    EXPECT_NEAR(fitted(3), 6.8, 1e-12);
}

// =============================================================================
// 3. 残差 (手算: [0.3; -0.4; -0.1; 0.2])
// =============================================================================
TEST(OlsHcTest, Residuals_MatchHandCalc) {
    MatrixXD X = make_X();
    VectorXD y = make_y();
    VectorXD beta = make_beta();
    VectorXD res = OLSEstimator{}.computeResiduals(X, y, beta);
    ASSERT_EQ(res.size(), 4u);
    EXPECT_NEAR(res(0), 0.3, 1e-12);
    EXPECT_NEAR(res(1), -0.4, 1e-12);
    EXPECT_NEAR(res(2), -0.1, 1e-12);
    EXPECT_NEAR(res(3), 0.2, 1e-12);
}

// =============================================================================
// 4. R² (手算: 1 - 0.30/14.75 = 289/295)
// =============================================================================
TEST(OlsHcTest, RSquared_MatchHandCalc) {
    OLSEstimator ols;
    EstimationResult r = ols.estimate(make_data());
    EXPECT_NEAR(r.r_squared, 289.0 / 295.0, 1e-12);
}

// =============================================================================
// 5. 调整后 R² (手算: 1 - (1-R²)·3/2 = 286/295)
// =============================================================================
TEST(OlsHcTest, AdjustedRSquared_MatchHandCalc) {
    OLSEstimator ols;
    EstimationResult r = ols.estimate(make_data());
    EXPECT_NEAR(r.adj_r_squared, 286.0 / 295.0, 1e-12);
    // 直接方法对照
    EXPECT_NEAR(ols.computeAdjustedRSquared(289.0 / 295.0, 4, 2), 286.0 / 295.0, 1e-12);
}

// =============================================================================
// 6. F 统计量 (手算: R²/(1-R²)·(N-K)/K = 289/6)
// =============================================================================
TEST(OlsHcTest, FStatistic_MatchHandCalc) {
    OLSEstimator ols;
    EXPECT_NEAR(ols.computeFStatistic(289.0 / 295.0, 4, 2), 289.0 / 6.0, 1e-12);
}

// =============================================================================
// 7. 投影矩阵对角 = leverage h_i (手算: [0.7; 0.3; 0.3; 0.7])
// =============================================================================
TEST(OlsHcTest, ProjectionMatrix_DiagonalIsLeverage) {
    MatrixXD X = make_X();
    MatrixXD A = make_XtX_inv();
    MatrixXD P = OLSEstimator{}.computeProjectionMatrix(X, A);
    ASSERT_EQ(P.rows(), 4u);
    ASSERT_EQ(P.cols(), 4u);
    EXPECT_NEAR(P(0, 0), 0.7, 1e-12);
    EXPECT_NEAR(P(1, 1), 0.3, 1e-12);
    EXPECT_NEAR(P(2, 2), 0.3, 1e-12);
    EXPECT_NEAR(P(3, 3), 0.7, 1e-12);
}

// =============================================================================
// 8. 投影矩阵对称 + 幂等 (P²=P) + 迹=K
// =============================================================================
TEST(OlsHcTest, ProjectionMatrix_SymmetricIdempotent_TraceK) {
    MatrixXD X = make_X();
    MatrixXD A = make_XtX_inv();
    MatrixXD P = OLSEstimator{}.computeProjectionMatrix(X, A);
    // 对称
    for (Size i = 0; i < 4; ++i)
        for (Size j = 0; j < 4; ++j)
            EXPECT_NEAR(P(i, j), P(j, i), 1e-12);
    // 幂等: P*P = P
    MatrixXD PP = MatrixXD(P.eigen() * P.eigen());
    for (Size i = 0; i < 4; ++i)
        for (Size j = 0; j < 4; ++j)
            EXPECT_NEAR(PP(i, j), P(i, j), 1e-12);
    // 迹 = K = 2
    Real tr = P(0, 0) + P(1, 1) + P(2, 2) + P(3, 3);
    EXPECT_NEAR(tr, 2.0, 1e-12);
}

// =============================================================================
// 9. Classical 协方差 + 标准误 + t 统计量 + p 值
//   σ² = SSR/(N-K) = 0.15,  vcov = 0.15·(X'X)^{-1} = [9/40,-3/40; -3/40,3/100]
//   se = [√(9/40); √(3/100)],  t = [0; 17/√3]  (β_0=0 → t_0=0)
//   p (df=2 闭式): p = 1 - |t|/√(2+t²);  t_0=0 → p_0=1;  t_1=17/√3 → p_1=1-17/√295
// =============================================================================
TEST(OlsHcTest, ClassicalVcov_StdErrors_T_PValues) {
    OLSEstimator ols(CovarianceType::Classical);
    EstimationResult r = ols.estimate(make_data());
    EXPECT_EQ(r.cov_type, CovarianceType::Classical);
    // vcov = [9/40, -3/40; -3/40, 3/100] = [0.225, -0.075; -0.075, 0.03]
    EXPECT_NEAR(r.vcov(0, 0), 9.0 / 40.0, 1e-12);
    EXPECT_NEAR(r.vcov(0, 1), -3.0 / 40.0, 1e-12);
    EXPECT_NEAR(r.vcov(1, 1), 3.0 / 100.0, 1e-12);
    // 标准误
    EXPECT_NEAR(r.std_errors(0), std::sqrt(9.0 / 40.0), 1e-12);
    EXPECT_NEAR(r.std_errors(1), std::sqrt(3.0 / 100.0), 1e-12);
    // t 统计量 (手算: β_0=0 → t_0=0; β_1=1.7, se_1=√(3/100)=√3/10 → t_1=17/√3)
    EXPECT_NEAR(r.t_statistics(0), 0.0, 1e-12);
    EXPECT_NEAR(r.t_statistics(1), 17.0 / std::sqrt(3.0), 1e-10);
    // p 值 (df=2 解析闭式: p = 1 - |t|/√(2+t²), 独立于 incomplete beta 实现)
    // t_0=0 → p_0=1.0;  t_1=17/√3 → t_1²=289/3, p_1=1-(17/√3)/√(2+289/3)=1-17/√295
    Real t0 = 0.0;
    Real t1 = 17.0 / std::sqrt(3.0);
    Real p0 = 1.0 - t0 / std::sqrt(2.0 + t0 * t0);
    Real p1 = 1.0 - t1 / std::sqrt(2.0 + t1 * t1);
    EXPECT_NEAR(r.p_values(0), p0, 1e-10);
    EXPECT_NEAR(r.p_values(1), p1, 1e-10);
    EXPECT_TRUE(r.p_values(0) >= 0.0 && r.p_values(0) <= 1.0);
    EXPECT_TRUE(r.p_values(1) >= 0.0 && r.p_values(1) <= 1.0);
    // 高斯对数似然: σ²_MLE = SSR/N = 0.075, LL = -N/2·(log(2π)+1+log(σ²))
    Real ll = -0.5 * 4.0 * (std::log(kTwoPi) + 1.0 + std::log(0.075));
    EXPECT_NEAR(r.log_likelihood, ll, 1e-10);
}

// =============================================================================
// 10. HC0 (White 1980): V = (X'X)^{-1} (X' diag(u²) X) (X'X)^{-1}
//     手算 = [7/50, -41/1000; -41/1000, 67/5000]
// =============================================================================
TEST(OlsHcTest, HC0_Vcov_MatchHandCalc) {
    MatrixXD X = make_X();
    MatrixXD A = make_XtX_inv();
    VectorXD u = make_residuals();
    MatrixXD V = compute_hc_vcov(X, u, A, CovarianceType::HC0);
    EXPECT_NEAR(V(0, 0), 7.0 / 50.0, 1e-12);
    EXPECT_NEAR(V(0, 1), -41.0 / 1000.0, 1e-12);
    EXPECT_NEAR(V(1, 0), -41.0 / 1000.0, 1e-12);
    EXPECT_NEAR(V(1, 1), 67.0 / 5000.0, 1e-12);
    // estimate() 用 HC0 应一致
    OLSEstimator ols(CovarianceType::HC0);
    EstimationResult r = ols.estimate(make_data());
    EXPECT_EQ(r.cov_type, CovarianceType::HC0);
    EXPECT_NEAR(r.vcov(0, 0), 7.0 / 50.0, 1e-12);
    EXPECT_NEAR(r.vcov(1, 1), 67.0 / 5000.0, 1e-12);
}

// =============================================================================
// 11. HC1 (MacKinnon-White 1985): V = N/(N-K) · HC0
//     排幻觉点 E2: 分母是 N-K (非 n-k 笔误), N=4,K=2 → 2·HC0
//     手算 = [7/25, -41/500; -41/500, 67/2500]
// =============================================================================
TEST(OlsHcTest, HC1_Vcov_MatchHandCalc) {
    MatrixXD X = make_X();
    MatrixXD A = make_XtX_inv();
    VectorXD u = make_residuals();
    MatrixXD V = compute_hc_vcov(X, u, A, CovarianceType::HC1);
    EXPECT_NEAR(V(0, 0), 7.0 / 25.0, 1e-12);
    EXPECT_NEAR(V(0, 1), -41.0 / 500.0, 1e-12);
    EXPECT_NEAR(V(1, 0), -41.0 / 500.0, 1e-12);
    EXPECT_NEAR(V(1, 1), 67.0 / 2500.0, 1e-12);
    // HC1 = N/(N-K) · HC0 关系
    MatrixXD V0 = compute_hc_vcov(X, u, A, CovarianceType::HC0);
    EXPECT_NEAR(V(0, 0), (4.0 / 2.0) * V0(0, 0), 1e-12);
    EXPECT_NEAR(V(1, 1), (4.0 / 2.0) * V0(1, 1), 1e-12);
}

// =============================================================================
// 12. HC2 (MacKinnon-White 1985): V = (X'X)^{-1} (X' diag(u²/(1-h)) X) (X'X)^{-1}
//     排幻觉点 E3: leverage h_i = x_i'(X'X)^{-1}x_i
//     手算 = [41/105, -17/140; -17/140, 29/700]
// =============================================================================
TEST(OlsHcTest, HC2_Vcov_MatchHandCalc) {
    MatrixXD X = make_X();
    MatrixXD A = make_XtX_inv();
    VectorXD u = make_residuals();
    MatrixXD V = compute_hc_vcov(X, u, A, CovarianceType::HC2);
    EXPECT_NEAR(V(0, 0), 41.0 / 105.0, 1e-12);
    EXPECT_NEAR(V(0, 1), -17.0 / 140.0, 1e-12);
    EXPECT_NEAR(V(1, 0), -17.0 / 140.0, 1e-12);
    EXPECT_NEAR(V(1, 1), 29.0 / 700.0, 1e-12);
}

// =============================================================================
// 13. HC3 (Jackknife 近似): V = (X'X)^{-1} (X' diag(u²/(1-h)²) X) (X'X)^{-1}
//     手算 = [526/441, -563/1470; -563/1470, 327/2450]
// =============================================================================
TEST(OlsHcTest, HC3_Vcov_MatchHandCalc) {
    MatrixXD X = make_X();
    MatrixXD A = make_XtX_inv();
    VectorXD u = make_residuals();
    MatrixXD V = compute_hc_vcov(X, u, A, CovarianceType::HC3);
    EXPECT_NEAR(V(0, 0), 526.0 / 441.0, 1e-12);
    EXPECT_NEAR(V(0, 1), -563.0 / 1470.0, 1e-12);
    EXPECT_NEAR(V(1, 0), -563.0 / 1470.0, 1e-12);
    EXPECT_NEAR(V(1, 1), 327.0 / 2450.0, 1e-12);
}

// =============================================================================
// 14. HC4 (Cribari-Neto 2004): δ_i = min(N·h_i/K, 4), w_i = u_i²/(1-h_i)^δ_i
//     无理幂 → 独立参考实现 (显式循环) 交叉验证
// =============================================================================
TEST(OlsHcTest, HC4_Vcov_MatchesIndependentReference) {
    MatrixXD X = make_X();
    MatrixXD A = make_XtX_inv();
    VectorXD u = make_residuals();
    MatrixXD Vref = reference_hc45(X, u, A, CovarianceType::HC4);
    MatrixXD V = compute_hc_vcov(X, u, A, CovarianceType::HC4);
    for (Size i = 0; i < 2; ++i)
        for (Size j = 0; j < 2; ++j)
            EXPECT_NEAR(V(i, j), Vref(i, j), 1e-12);
    // 对称性
    EXPECT_NEAR(V(0, 1), V(1, 0), 1e-12);
    // δ_i 手算: N·h/K = 2·h = [1.4, 0.6, 0.6, 1.4], 均 < 4 → δ = [1.4,0.6,0.6,1.4]
    // HC4 对角 > 0 (正半定)
    EXPECT_GT(V(0, 0), 0.0);
    EXPECT_GT(V(1, 1), 0.0);
}

// =============================================================================
// 15. HC5 (Cribari-Neto-Souza 2007): δ_i = min(N·h_i/K, max(4, 0.7·N·h_max/K))
//     公式同 HC4 但 δ_i 不同 (任务描述; 注: spec §7.3 的 (δ/2) 为 R sandwich HC4m 行为, 非 HC5)
// =============================================================================
TEST(OlsHcTest, HC5_Vcov_MatchesIndependentReference) {
    MatrixXD X = make_X();
    MatrixXD A = make_XtX_inv();
    VectorXD u = make_residuals();
    MatrixXD Vref = reference_hc45(X, u, A, CovarianceType::HC5);
    MatrixXD V = compute_hc_vcov(X, u, A, CovarianceType::HC5);
    for (Size i = 0; i < 2; ++i)
        for (Size j = 0; j < 2; ++j)
            EXPECT_NEAR(V(i, j), Vref(i, j), 1e-12);
    EXPECT_NEAR(V(0, 1), V(1, 0), 1e-12);
    // h_max = 0.7, 0.7·N·h_max/K = 0.7·4·0.7/2 = 0.98 < 4 → max(4, 0.98) = 4
    // δ_i = min(2·h, 4) = min([1.4,0.6,0.6,1.4], 4) = [1.4,0.6,0.6,1.4] (与 HC4 相同)
    // 故 HC5 应等于 HC4 (本数据集 h_max 较小, 0.7·N·h_max/K < 4)
    MatrixXD V4 = compute_hc_vcov(X, u, A, CovarianceType::HC4);
    EXPECT_NEAR(V(0, 0), V4(0, 0), 1e-12);
    EXPECT_NEAR(V(1, 1), V4(1, 1), 1e-12);
}

// =============================================================================
// 16. 非法协方差类型 (Classical/HAC/Cluster 等) 传给 compute_hc_vcov 抛 invalid_argument
// =============================================================================
TEST(OlsHcTest, HC_InvalidType_Throws) {
    MatrixXD X = make_X();
    MatrixXD A = make_XtX_inv();
    VectorXD u = make_residuals();
    EXPECT_THROW(compute_hc_vcov(X, u, A, CovarianceType::Classical), std::invalid_argument);
    EXPECT_THROW(compute_hc_vcov(X, u, A, CovarianceType::HAC_Bartlett), std::invalid_argument);
    EXPECT_THROW(compute_hc_vcov(X, u, A, CovarianceType::Cluster_OneWay), std::invalid_argument);
}

// =============================================================================
// 17. N < K 报错 (自由度不足)
// =============================================================================
TEST(OlsHcTest, FewerObsThanParams_Throws) {
    MatrixXD X(2, 3);  // n=2 < k=3
    X(0, 0) = 1.0; X(0, 1) = 1.0; X(0, 2) = 2.0;
    X(1, 0) = 1.0; X(1, 1) = 2.0; X(1, 2) = 3.0;
    VectorXD y(2);
    y(0) = 1.0; y(1) = 2.0;
    EconData data = make_cross_section(X, y, {"a", "b", "c"}, "y");
    OLSEstimator ols;
    EXPECT_THROW(ols.estimate(data), std::invalid_argument);
}

// =============================================================================
// 18. X 奇异 (列共线) 报错 — (X'X) 非 SPD, inverse_symmetric(LLT) 抛 runtime_error
// =============================================================================
TEST(OlsHcTest, SingularX_Throws) {
    MatrixXD X(3, 2);  // 第二列 = 第一列 (秩 1, X'X 奇异)
    X(0, 0) = 1.0; X(0, 1) = 1.0;
    X(1, 0) = 1.0; X(1, 1) = 1.0;
    X(2, 0) = 1.0; X(2, 1) = 1.0;
    VectorXD y(3);
    y(0) = 1.0; y(1) = 2.0; y(2) = 3.0;
    EconData data = make_cross_section(X, y, {"a", "b"}, "y");
    OLSEstimator ols;
    EXPECT_THROW(ols.estimate(data), std::runtime_error);
}

// =============================================================================
// 19. PanelData / TimeSeriesData 输入报错 (OLS 仅支持 CrossSectionData)
// =============================================================================
TEST(OlsHcTest, NonCrossSectionData_Throws) {
    MatrixXD X = make_X();
    VectorXD y = make_y();
    OLSEstimator ols;
    // PanelData
    PanelData pd;
    pd.X = X; pd.y = y;
    pd.entity_id = {0, 0, 1, 1};
    pd.time_id = {0, 1, 0, 1};
    pd.x_names = {"intercept", "x"};
    pd.y_name = "y";
    pd.balanced = true;
    EXPECT_THROW(ols.estimate(EconData{pd}), std::invalid_argument);
    // TimeSeriesData
    TimeSeriesData ts;
    ts.y = y; ts.X = X;
    ts.timestamps = {0.0, 1.0, 2.0, 3.0};
    ts.x_names = {"intercept", "x"};
    ts.y_name = "y";
    EXPECT_THROW(ols.estimate(EconData{ts}), std::invalid_argument);
}

// =============================================================================
// 20. clone() 多态 — 配置保留, 通过基类指针估计
// =============================================================================
TEST(OlsHcTest, Clone_Polymorphic) {
    OLSEstimator ols(CovarianceType::HC3);
    std::unique_ptr<Estimator> clone = ols.clone();
    EXPECT_EQ(clone->name(), "OLS");
    EXPECT_EQ(clone->estimatorClass(), EstimatorClass::Parametric);
    EXPECT_EQ(clone->covarianceType(), CovarianceType::HC3);
    // 通过基类指针执行估计 (多态)
    EstimationResult r = clone->estimate(make_data());
    EXPECT_EQ(r.cov_type, CovarianceType::HC3);
    EXPECT_NEAR(r.coefficients(0), 0.0, 1e-12);
    EXPECT_NEAR(r.coefficients(1), 1.7, 1e-12);
    // HC3 协方差一致
    EXPECT_NEAR(r.vcov(0, 0), 526.0 / 441.0, 1e-12);
    EXPECT_NEAR(r.vcov(1, 1), 327.0 / 2450.0, 1e-12);
    // clone 独立于原对象
    clone->setCovarianceType(CovarianceType::HC0);
    EXPECT_EQ(ols.covarianceType(), CovarianceType::HC3);
    EXPECT_EQ(clone->covarianceType(), CovarianceType::HC0);
}
