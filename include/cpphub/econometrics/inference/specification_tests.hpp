// =============================================================================
// specification_tests.hpp - 模型设定与预测检验 (P1)
//
// Phase 7A Wave 2: 通用证伪统计量
//
// 包含 3 个检验:
//   1. 信息矩阵等式检验 (White 1982)
//   2. Mincer-Zarnowitz 回归 (Mincer-Zarnowitz 1969)
//   3. Diebold-Mariano 检验 (Diebold-Mariano 1995, HLN 1997 修正)
//
// ADR-015 方案 B: 仅依赖 core/, 不依赖 linalg_dynamic.hpp (Eigen3)
// MZ 回归用 detail/ols_simple.hpp, 信息矩阵用 Gauss-Jordan 自实现
//
// 教材锚点: White 1982, Mincer-Zarnowitz 1969, Diebold-Mariano 1995, HLN 1997
// 排幻觉点: H9(IM对QMLE均值方程)/H10(MZ R²也是预测精度)/H11(HLN修正1/N, t分布)
// =============================================================================
#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/special_functions.hpp"
#include "cpphub/econometrics/inference/detail/test_result_base.hpp"
#include "cpphub/econometrics/inference/detail/ols_simple.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// =============================================================================
// 信息矩阵等式检验 (White 1982)
// =============================================================================

// IM = N · vech(D)' [Â]^{-1} vech(D) ~ χ²(q)
// 其中 D = I_outer - I_inner, q = K(K+1)/2
//   I_outer = (1/N) Σ s_t s_t'  (外积形式信息矩阵)
//   I_inner = -H/N               (Hessian 形式信息矩阵)
// 排幻觉点 H9: 对 QMLE 检验均值方程正确性, 不是方差方程
struct InformationMatrixResult {
    detail::TestResultBase base;  // statistic=IM, p_value, method="Information Matrix", reject_null
    Size df;  // q = K(K+1)/2
};

inline InformationMatrixResult information_matrix_test(
    const std::vector<std::vector<Real>>& scores,   // N×K 得分矩阵
    const std::vector<std::vector<Real>>& hessian) { // K×K Hessian 矩阵

    const Size N = scores.size();
    if (N < 5) {
        throw std::invalid_argument("information_matrix_test: need at least 5 observations");
    }
    const Size K = scores[0].size();
    if (K == 0) {
        throw std::invalid_argument("information_matrix_test: empty scores");
    }
    if (hessian.size() != K || hessian[0].size() != K) {
        throw std::invalid_argument("information_matrix_test: hessian dimension mismatch");
    }

    // 验证 scores 每行长度一致
    for (Size t = 0; t < N; ++t) {
        if (scores[t].size() != K) {
            throw std::invalid_argument("information_matrix_test: scores row size mismatch");
        }
    }

    // I_outer = (1/N) Σ s_t s_t'
    std::vector<std::vector<Real>> I_outer(K, std::vector<Real>(K, 0.0));
    for (Size t = 0; t < N; ++t) {
        for (Size i = 0; i < K; ++i) {
            for (Size j = 0; j < K; ++j) {
                I_outer[i][j] += scores[t][i] * scores[t][j];
            }
        }
    }
    for (Size i = 0; i < K; ++i) {
        for (Size j = 0; j < K; ++j) {
            I_outer[i][j] /= static_cast<Real>(N);
        }
    }

    // I_inner = -H / N
    std::vector<std::vector<Real>> I_inner(K, std::vector<Real>(K, 0.0));
    for (Size i = 0; i < K; ++i) {
        for (Size j = 0; j < K; ++j) {
            I_inner[i][j] = -hessian[i][j] / static_cast<Real>(N);
        }
    }

    // D = I_outer - I_inner
    std::vector<std::vector<Real>> D(K, std::vector<Real>(K, 0.0));
    for (Size i = 0; i < K; ++i) {
        for (Size j = 0; j < K; ++j) {
            D[i][j] = I_outer[i][j] - I_inner[i][j];
        }
    }

    // vech(D) - 半向量化 (下三角, 含对角线)
    const Size q = K * (K + 1) / 2;
    std::vector<Real> vech_D(q);
    Size idx = 0;
    for (Size i = 0; i < K; ++i) {
        for (Size j = 0; j <= i; ++j) {
            vech_D[idx++] = D[i][j];
        }
    }

    // Â = (1/N) Σ [vech(s_t s_t') - vech(I_outer)] [vech(s_t s_t') - vech(I_outer)]'
    // (White 1982 协方差矩阵估计, 中心化用 I_outer 因 E[vech(s_t s_t')]=vech(I_outer))
    std::vector<Real> vech_I_outer(q);
    idx = 0;
    for (Size i = 0; i < K; ++i) {
        for (Size j = 0; j <= i; ++j) {
            vech_I_outer[idx++] = I_outer[i][j];
        }
    }
    std::vector<std::vector<Real>> A_hat(q, std::vector<Real>(q, 0.0));
    std::vector<Real> vech_st(q);
    for (Size t = 0; t < N; ++t) {
        idx = 0;
        for (Size i = 0; i < K; ++i) {
            for (Size j = 0; j <= i; ++j) {
                vech_st[idx++] = scores[t][i] * scores[t][j];
            }
        }
        for (Size a = 0; a < q; ++a) {
            const Real diff_a = vech_st[a] - vech_I_outer[a];
            for (Size b = a; b < q; ++b) {
                A_hat[a][b] += diff_a * (vech_st[b] - vech_I_outer[b]);
            }
        }
    }
    for (Size a = 0; a < q; ++a) {
        for (Size b = a; b < q; ++b) {
            A_hat[a][b] /= static_cast<Real>(N);
            A_hat[b][a] = A_hat[a][b];  // 对称化
        }
    }

    // 求解 Â x = vech_D (Gauss-Jordan, partial pivoting)
    std::vector<std::vector<Real>> A_solve = A_hat;
    std::vector<Real> b_solve = vech_D;
    for (Size col = 0; col < q; ++col) {
        Size piv = col;
        Real maxv = std::fabs(A_solve[col][col]);
        for (Size r = col + 1; r < q; ++r) {
            Real v = std::fabs(A_solve[r][col]);
            if (v > maxv) { maxv = v; piv = r; }
        }
        if (maxv < 1e-15) {
            // Â 奇异, 降级: 用对角元素
            throw std::runtime_error("information_matrix_test: singular covariance matrix");
        }
        if (piv != col) {
            std::swap(A_solve[piv], A_solve[col]);
            std::swap(b_solve[piv], b_solve[col]);
        }
        Real akk = A_solve[col][col];
        for (Size r = col + 1; r < q; ++r) {
            Real f = A_solve[r][col] / akk;
            if (f == 0.0) continue;
            for (Size c = col; c < q; ++c) A_solve[r][c] -= f * A_solve[col][c];
            b_solve[r] -= f * b_solve[col];
        }
    }
    std::vector<Real> x(q);
    for (Size i = q; i-- > 0;) {
        Real s = b_solve[i];
        for (Size j = i + 1; j < q; ++j) s -= A_solve[i][j] * x[j];
        x[i] = s / A_solve[i][i];
    }

    // IM = N · vech_D' x = N · vech_D' Â^{-1} vech_D
    Real im = 0.0;
    for (Size a = 0; a < q; ++a) {
        im += vech_D[a] * x[a];
    }
    im *= static_cast<Real>(N);

    const Real p_value = detail::chi2_sf(static_cast<Real>(q), im);

    InformationMatrixResult result;
    result.base.statistic = im;
    result.base.p_value = p_value;
    result.base.method_name = "Information Matrix";
    result.base.reject_null = (p_value < 0.05);
    result.df = q;
    return result;
}

// =============================================================================
// Mincer-Zarnowitz 回归 (Mincer-Zarnowitz 1969)
// =============================================================================

// y_t = α + β·ŷ_t + ε_t, H0: α=0 且 β=1
// joint F = ((RSS_r - RSS_u) / 2) / (RSS_u / (N-2)) ~ F(2, N-2)
// 排幻觉点 H10: R² 也是预测精度度量 (越大预测越准)
struct MincerZarnowitzResult {
    detail::TestResultBase base;  // statistic=joint F, p_value, method="Mincer-Zarnowitz", reject_null
    Real alpha;
    Real beta;
    Real alpha_t_stat;
    Real beta_t_stat;
    Real r_squared;
};

inline MincerZarnowitzResult mincer_zarnowitz_regression(
    const std::vector<Real>& actual,
    const std::vector<Real>& forecast) {

    const Size N = actual.size();
    if (N < 5) {
        throw std::invalid_argument("mincer_zarnowitz_regression: need at least 5 observations");
    }
    if (forecast.size() != N) {
        throw std::invalid_argument("mincer_zarnowitz_regression: size mismatch");
    }

    // OLS: actual = α + β·forecast + ε
    // 设计矩阵 X = [1, forecast]
    std::vector<std::vector<Real>> X(N, std::vector<Real>(2));
    for (Size t = 0; t < N; ++t) {
        X[t][0] = 1.0;
        X[t][1] = forecast[t];
    }
    std::vector<Real> fitted, resid;
    Real r2;
    std::vector<Real> beta_hat = detail::ols_simple(actual, X, fitted, resid, r2);

    const Real alpha = beta_hat[0];
    const Real beta = beta_hat[1];

    // RSS_unrestricted = Σ resid²
    Real rss_u = 0.0;
    for (Size t = 0; t < N; ++t) rss_u += resid[t] * resid[t];

    // RSS_restricted = Σ (actual - forecast)² (强制 α=0, β=1)
    Real rss_r = 0.0;
    for (Size t = 0; t < N; ++t) {
        const Real d = actual[t] - forecast[t];
        rss_r += d * d;
    }

    // joint F = ((RSS_r - RSS_u) / 2) / (RSS_u / (N-2)) ~ F(2, N-2)
    Real f_stat = 0.0;
    if (rss_u > 1e-300 && N > 2) {
        f_stat = ((rss_r - rss_u) / 2.0) / (rss_u / static_cast<Real>(N - 2));
    }

    // p_value = F_sf(2, N-2, f_stat) = 1 - F_cdf(2, N-2, f_stat)
    const Real df1 = 2.0;
    const Real df2 = static_cast<Real>(N - 2);
    // F cdf: P(F <= f) = I_{df1*f/(df1*f+df2)}(df1/2, df2/2)
    // F_sf = 1 - F_cdf = I_{df2/(df1*f+df2)}(df2/2, df1/2) = beta_i(df2/2, df1/2, df2/(df1*f+df2))
    Real p_value = 1.0;
    if (f_stat > 0.0) {
        const Real x = df2 / (df1 * f_stat + df2);
        p_value = detail::beta_i(df2 / 2.0, df1 / 2.0, x);
    }

    // t 统计量 (α=0, β=1 的单参数检验)
    // SE(α) = sqrt(σ² · (Σf²)/(N·Σ(f-f̄)²))
    // SE(β) = sqrt(σ² / Σ(f-f̄)²)
    // σ² = RSS_u / (N-2)
    const Real sigma2 = rss_u / static_cast<Real>(N - 2);
    Real f_mean = 0.0;
    for (Size t = 0; t < N; ++t) f_mean += forecast[t];
    f_mean /= static_cast<Real>(N);
    Real ss_f = 0.0, sf2 = 0.0;
    for (Size t = 0; t < N; ++t) {
        const Real d = forecast[t] - f_mean;
        ss_f += d * d;
        sf2 += forecast[t] * forecast[t];
    }
    Real se_alpha = 0.0, se_beta = 0.0;
    if (ss_f > 1e-300) {
        se_beta = std::sqrt(sigma2 / ss_f);
        // SE(α) = sqrt(σ² · (1/N + f̄²/SS_f))
        se_alpha = std::sqrt(sigma2 * (1.0 / static_cast<Real>(N) + f_mean * f_mean / ss_f));
    }
    const Real alpha_t = (se_alpha > 0.0) ? alpha / se_alpha : 0.0;
    const Real beta_t = (se_beta > 0.0) ? (beta - 1.0) / se_beta : 0.0;  // H0: β=1

    MincerZarnowitzResult result;
    result.base.statistic = f_stat;
    result.base.p_value = p_value;
    result.base.method_name = "Mincer-Zarnowitz";
    result.base.reject_null = (p_value < 0.05);
    result.alpha = alpha;
    result.beta = beta;
    result.alpha_t_stat = alpha_t;
    result.beta_t_stat = beta_t;
    result.r_squared = r2;
    return result;
}

// =============================================================================
// Diebold-Mariano 检验 (Diebold-Mariano 1995, HLN 1997 修正)
// =============================================================================

// DM = d̄ / √(V̂/N) ~ t(N-1)
// d_t = L(e1_t) - L(e2_t)
// V̂ = γ̂_0 + 2 Σ_{h=1}^{h-1} γ̂_h  (HLN 修正)
// γ̂_h = (1/N) Σ (d_t - d̄)(d_{t-h} - d̄)  (排幻觉点 H11: 用 1/N 而非 1/(N-h))
// 排幻觉点 H11: HLN 修正用 1/N 计算 γ̂_h, DM ~ t(N-1) 非 N(0,1)
struct DieboldMarianoResult {
    detail::TestResultBase base;  // statistic=DM, p_value, method="Diebold-Mariano", reject_null
    Real mean_loss_diff;
};

inline DieboldMarianoResult diebold_mariano_test(
    const std::vector<Real>& actual,
    const std::vector<Real>& forecast1,
    const std::vector<Real>& forecast2,
    const std::string& loss_function = "mse",  // "mse"/"mae"/"hmahe"
    Size h = 1) {

    const Size N = actual.size();
    if (N < 5) {
        throw std::invalid_argument("diebold_mariano_test: need at least 5 observations");
    }
    if (forecast1.size() != N || forecast2.size() != N) {
        throw std::invalid_argument("diebold_mariano_test: size mismatch");
    }
    if (h < 1 || h >= N) {
        throw std::invalid_argument("diebold_mariano_test: invalid h");
    }

    // 计算损失差 d_t = L(e1_t) - L(e2_t)
    std::vector<Real> d(N);
    for (Size t = 0; t < N; ++t) {
        const Real e1 = actual[t] - forecast1[t];
        const Real e2 = actual[t] - forecast2[t];
        Real L1, L2;
        if (loss_function == "mse") {
            L1 = e1 * e1;
            L2 = e2 * e2;
        } else if (loss_function == "mae") {
            L1 = std::fabs(e1);
            L2 = std::fabs(e2);
        } else if (loss_function == "hmahe") {
            // HMAE (Harvey-Leybourne-Newbold): |e|/(actual²+forecast²) 的变体
            // 标准 HMAE: 2|e|/(|actual|+|forecast|)
            L1 = 2.0 * std::fabs(e1) / (std::fabs(actual[t]) + std::fabs(forecast1[t]) + 1e-300);
            L2 = 2.0 * std::fabs(e2) / (std::fabs(actual[t]) + std::fabs(forecast2[t]) + 1e-300);
        } else {
            throw std::invalid_argument("diebold_mariano_test: unknown loss_function");
        }
        d[t] = L1 - L2;
    }

    // d̄
    Real d_mean = 0.0;
    for (Size t = 0; t < N; ++t) d_mean += d[t];
    d_mean /= static_cast<Real>(N);

    // γ̂_h (排幻觉点 H11: 用 1/N 而非 1/(N-h))
    // V̂ = γ̂_0 + 2 Σ_{h=1}^{h-1} γ̂_h
    Real gamma0 = 0.0;
    for (Size t = 0; t < N; ++t) {
        const Real diff = d[t] - d_mean;
        gamma0 += diff * diff;
    }
    gamma0 /= static_cast<Real>(N);

    Real v_hat = gamma0;
    for (Size lag = 1; lag < h; ++lag) {
        Real gamma_h = 0.0;
        for (Size t = lag; t < N; ++t) {
            gamma_h += (d[t] - d_mean) * (d[t - lag] - d_mean);
        }
        gamma_h /= static_cast<Real>(N);  // 排幻觉点 H11: 1/N
        v_hat += 2.0 * gamma_h;
    }

    if (v_hat <= 0.0) {
        // V̂ = 0: 两个预测完全相同 (d_t 全等于 d̄) → DM=0, p_value=1
        DieboldMarianoResult result;
        result.base.statistic = 0.0;
        result.base.p_value = 1.0;
        result.base.method_name = "Diebold-Mariano";
        result.base.reject_null = false;
        result.mean_loss_diff = d_mean;
        return result;
    }

    // DM = d̄ / √(V̂/N) ~ t(N-1) (排幻觉点 H11: t 分布非正态)
    const Real dm = d_mean / std::sqrt(v_hat / static_cast<Real>(N));

    // p_value (双侧, t 分布): p = beta_i(df/2, 0.5, df/(df+t²))
    const Real df = static_cast<Real>(N - 1);
    const Real x = df / (df + dm * dm);
    const Real p_value = detail::beta_i(df / 2.0, 0.5, x);

    DieboldMarianoResult result;
    result.base.statistic = dm;
    result.base.p_value = p_value;
    result.base.method_name = "Diebold-Mariano";
    result.base.reject_null = (p_value < 0.05);
    result.mean_loss_diff = d_mean;
    return result;
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
