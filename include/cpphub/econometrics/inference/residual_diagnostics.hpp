// =============================================================================
// residual_diagnostics.hpp - 残差诊断 (P0)
//
// Phase 7A Wave 1: 通用证伪统计量 (ADR-015 方案 B)
//
// 包含 5 个检验:
//   1. Jarque-Bera 正态性检验 (Jarque-Bera 1987)
//   2. Ljung-Box 残差自相关检验 (Ljung-Box 1978)
//   3. Breusch-Godfrey LM 自相关检验 (Breusch 1978, Godfrey 1978)
//   4. Breusch-Pagan 异方差检验 (Breusch-Pagan 1979, Koenker 1981 修正)
//   5. White 异方差检验 (White 1980)
//
// ADR-015 方案 B: 仅依赖 core/, 不依赖 linalg_dynamic.hpp (Eigen3)
// 回归检验 (BG/BP/White) 辅助回归用 detail/ols_simple.hpp
//
// 教材锚点: Greene 8ed §4.8/§9.5/§13.7, Tsay 3ed §2
// 排幻觉点: H1(σ有偏)/H2(峰度非超额)/H3(LB加权)/H4(lag自动)/H5(BG含X)/H6(BP用Koenker)/H7(White高维检查)
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
// Jarque-Bera 正态性检验
// =============================================================================

// JB = N * [S²/6 + (K-3)²/24] ~ χ²(2)
// 排幻觉点 H1: 分母用样本标准差 s (有偏, 除以 N), 非 σ (除以 N-1)
// 排幻觉点 H2: 公式用峰度 K, K-3 才是超额峰度
struct JarqueBeraResult {
    detail::TestResultBase base;  // statistic=JB, p_value, method="Jarque-Bera", reject_null
    Real skewness;
    Real kurtosis;        // 非超额峰度
};

inline JarqueBeraResult jarque_bera_test(const std::vector<Real>& residuals) {
    const Size n = residuals.size();
    if (n < 3) {
        throw std::invalid_argument("jarque_bera_test: need at least 3 observations");
    }

    // 均值
    Real mean = 0.0;
    for (Size i = 0; i < n; ++i) mean += residuals[i];
    mean /= static_cast<Real>(n);

    // 有偏方差 (除以 N, 排幻觉点 H1)
    Real var = 0.0;
    for (Size i = 0; i < n; ++i) {
        const Real d = residuals[i] - mean;
        var += d * d;
    }
    var /= static_cast<Real>(n);

    if (var <= 0.0) {
        throw std::runtime_error("jarque_bera_test: zero variance");
    }
    const Real s = std::sqrt(var);

    // 偏度 S = (1/N) Σ ((x_i - x̄)/s)³
    Real skew = 0.0;
    for (Size i = 0; i < n; ++i) {
        const Real z = (residuals[i] - mean) / s;
        skew += z * z * z;
    }
    skew /= static_cast<Real>(n);

    // 峰度 K = (1/N) Σ ((x_i - x̄)/s)⁴ (非超额, 排幻觉点 H2)
    Real kurt = 0.0;
    for (Size i = 0; i < n; ++i) {
        const Real z = (residuals[i] - mean) / s;
        kurt += z * z * z * z;
    }
    kurt /= static_cast<Real>(n);

    // JB = N * [S²/6 + (K-3)²/24]
    const Real jb = static_cast<Real>(n) *
                    (skew * skew / 6.0 + (kurt - 3.0) * (kurt - 3.0) / 24.0);
    const Real p_value = detail::chi2_sf(2.0, jb);

    JarqueBeraResult result;
    result.base.statistic = jb;
    result.base.p_value = p_value;
    result.base.method_name = "Jarque-Bera";
    result.base.reject_null = (p_value < 0.05);
    result.skewness = skew;
    result.kurtosis = kurt;
    return result;
}

// =============================================================================
// Ljung-Box 残差自相关检验
// =============================================================================

// LB = N(N+2) Σ ρ_h²/(N-h) ~ χ²(m)
// 排幻觉点 H3: Ljung-Box 用 N(N+2)/(N-h) 加权, 非 Box-Pierce 的 N
// 排幻觉点 H4: lag=0 时自动选择 m = min(10, N/5) (Hyndman 推荐)
struct LjungBoxResult {
    detail::TestResultBase base;  // statistic=LB, p_value, method="Ljung-Box", reject_null
    Size lag;
    std::vector<Real> autocorrelations;
};

inline LjungBoxResult ljung_box_test(const std::vector<Real>& residuals, Size lag = 0) {
    const Size n = residuals.size();
    if (n < 5) {
        throw std::invalid_argument("ljung_box_test: need at least 5 observations");
    }

    // 排幻觉点 H4: lag 自动选择
    Size m = lag;
    if (m == 0) {
        m = std::min(static_cast<Size>(10), n / 5);
        if (m == 0) m = 1;
    }
    if (m >= n) {
        throw std::invalid_argument("ljung_box_test: lag must be less than N");
    }

    // 均值
    Real mean = 0.0;
    for (Size i = 0; i < n; ++i) mean += residuals[i];
    mean /= static_cast<Real>(n);

    // 自协方差 γ_h = (1/N) Σ (x_t - x̄)(x_{t-h} - x̄)
    // 自相关 ρ_h = γ_h / γ_0
    std::vector<Real> acf(m);
    Real gamma0 = 0.0;
    for (Size t = 0; t < n; ++t) {
        const Real d = residuals[t] - mean;
        gamma0 += d * d;
    }
    gamma0 /= static_cast<Real>(n);

    if (gamma0 <= 0.0) {
        throw std::runtime_error("ljung_box_test: zero variance");
    }

    for (Size h = 1; h <= m; ++h) {
        Real gamma_h = 0.0;
        for (Size t = h; t < n; ++t) {
            gamma_h += (residuals[t] - mean) * (residuals[t - h] - mean);
        }
        gamma_h /= static_cast<Real>(n);
        acf[h - 1] = gamma_h / gamma0;
    }

    // LB = N(N+2) Σ ρ_h²/(N-h) (排幻觉点 H3)
    Real lb = 0.0;
    for (Size h = 1; h <= m; ++h) {
        lb += acf[h - 1] * acf[h - 1] / static_cast<Real>(n - h);
    }
    lb *= static_cast<Real>(n) * static_cast<Real>(n + 2);

    const Real p_value = detail::chi2_sf(static_cast<Real>(m), lb);

    LjungBoxResult result;
    result.base.statistic = lb;
    result.base.p_value = p_value;
    result.base.method_name = "Ljung-Box";
    result.base.reject_null = (p_value < 0.05);
    result.lag = m;
    result.autocorrelations = acf;
    return result;
}

// =============================================================================
// Breusch-Godfrey LM 自相关检验
// =============================================================================

// 辅助回归 e_t = X_t γ + Σ δ_h e_{t-h} + u, LM = N'·R²_aux ~ χ²(p)
// 排幻觉点 H5: 辅助回归必须包含原 X
struct BreuschGodfreyResult {
    detail::TestResultBase base;  // statistic=LM, p_value, method="Breusch-Godfrey", reject_null
    Size lag;
};

// ADR-015 方案 B: X 用 std::vector<std::vector<Real>>, 不依赖 Eigen3
// 辅助回归内部用 detail/ols_simple()
inline BreuschGodfreyResult breusch_godfrey_test(
    const std::vector<std::vector<Real>>& X,  // N×K 解释变量矩阵 (不含常数列)
    const std::vector<Real>& residuals,
    Size lag) {

    const Size n = residuals.size();
    if (n < 5) {
        throw std::invalid_argument("breusch_godfrey_test: need at least 5 observations");
    }
    if (X.size() != n) {
        throw std::invalid_argument("breusch_godfrey_test: X and residuals size mismatch");
    }
    if (lag == 0 || lag >= n) {
        throw std::invalid_argument("breusch_godfrey_test: invalid lag");
    }

    const Size k = X.empty() ? 0 : X[0].size();
    const Size n_eff = n - lag;  // 有效观测数 (丢弃前 lag 个)

    // 构造辅助回归设计矩阵 [1, X, e_{t-1}, ..., e_{t-lag}]
    // 行数 = n_eff, 列数 = 1 + k + lag
    const Size p_cols = 1 + k + lag;
    if (n_eff <= p_cols) {
        throw std::invalid_argument("breusch_godfrey_test: insufficient observations after lag");
    }

    std::vector<std::vector<Real>> design(n_eff, std::vector<Real>(p_cols));
    std::vector<Real> y_eff(n_eff);

    for (Size i = 0; i < n_eff; ++i) {
        const Size t = i + lag;  // 原始时间索引
        design[i][0] = 1.0;  // 常数项
        for (Size j = 0; j < k; ++j) {
            design[i][1 + j] = X[t][j];  // 排幻觉点 H5: 必须包含原 X
        }
        for (Size h = 1; h <= lag; ++h) {
            design[i][1 + k + (h - 1)] = residuals[t - h];  // 滞后残差
        }
        y_eff[i] = residuals[t];
    }

    std::vector<Real> fitted, resid;
    Real r2;
    detail::ols_simple(y_eff, design, fitted, resid, r2);

    // LM = N' * R²_aux ~ χ²(lag)
    const Real lm = static_cast<Real>(n_eff) * r2;
    const Real p_value = detail::chi2_sf(static_cast<Real>(lag), lm);

    BreuschGodfreyResult result;
    result.base.statistic = lm;
    result.base.p_value = p_value;
    result.base.method_name = "Breusch-Godfrey";
    result.base.reject_null = (p_value < 0.05);
    result.lag = lag;
    return result;
}

// =============================================================================
// Breusch-Pagan 异方差检验 (Koenker 1981 修正)
// =============================================================================

// 辅助回归 e²_t = 常数 + X_t γ + u, LM = N·R²_aux ~ χ²(K)
// 排幻觉点 H6: 默认 Koenker 修正 (用 e², 非 e²/σ²), R bptest 默认也是 studentized
struct BreuschPaganResult {
    detail::TestResultBase base;  // statistic=LM, p_value, method="Breusch-Pagan", reject_null
};

inline BreuschPaganResult breusch_pagan_test(
    const std::vector<std::vector<Real>>& X,  // ADR-015 方案 B: 不依赖 Eigen3
    const std::vector<Real>& residuals) {

    const Size n = residuals.size();
    if (n < 5) {
        throw std::invalid_argument("breusch_pagan_test: need at least 5 observations");
    }
    if (X.size() != n) {
        throw std::invalid_argument("breusch_pagan_test: X and residuals size mismatch");
    }

    const Size k = X.empty() ? 0 : X[0].size();
    if (k == 0) {
        throw std::invalid_argument("breusch_pagan_test: empty design matrix");
    }

    // 构造辅助回归: e²_t = 常数 + X_t γ + u (排幻觉点 H6: 用 e²)
    std::vector<std::vector<Real>> design(n, std::vector<Real>(1 + k));
    std::vector<Real> e2(n);
    for (Size i = 0; i < n; ++i) {
        design[i][0] = 1.0;  // 常数项
        for (Size j = 0; j < k; ++j) {
            design[i][1 + j] = X[i][j];
        }
        e2[i] = residuals[i] * residuals[i];
    }

    std::vector<Real> fitted, resid;
    Real r2;
    detail::ols_simple(e2, design, fitted, resid, r2);

    // LM = N * R²_aux ~ χ²(K) (K = X 的列数, 不含常数)
    const Real lm = static_cast<Real>(n) * r2;
    const Real p_value = detail::chi2_sf(static_cast<Real>(k), lm);

    BreuschPaganResult result;
    result.base.statistic = lm;
    result.base.p_value = p_value;
    result.base.method_name = "Breusch-Pagan";
    result.base.reject_null = (p_value < 0.05);
    return result;
}

// =============================================================================
// White 异方差检验
// =============================================================================

// 辅助回归 e² = Zγ + u, Z = [常数, X, X²交叉项, X²]
// 排幻觉点 H7: 高维 q = K(K+1)/2, N > q 强制检查
struct WhiteResult {
    detail::TestResultBase base;  // statistic=LM, p_value, method="White", reject_null
};

inline WhiteResult white_test(
    const std::vector<std::vector<Real>>& X,  // ADR-015 方案 B: 不依赖 Eigen3
    const std::vector<Real>& residuals,
    bool include_cross_terms = true) {

    const Size n = residuals.size();
    if (n < 5) {
        throw std::invalid_argument("white_test: need at least 5 observations");
    }
    if (X.size() != n) {
        throw std::invalid_argument("white_test: X and residuals size mismatch");
    }

    const Size k = X.empty() ? 0 : X[0].size();
    if (k == 0) {
        throw std::invalid_argument("white_test: empty design matrix");
    }

    // 构造 Z = [常数, X, 交叉项, X²]
    // 列数: 1 + k + (include_cross_terms ? k*(k-1)/2 : 0) + k
    Size n_cross = 0;
    if (include_cross_terms && k >= 2) {
        n_cross = k * (k - 1) / 2;
    }
    const Size q = 1 + k + n_cross + k;  // 总列数 (含常数)

    // 排幻觉点 H7: N > q 强制检查
    if (n <= q) {
        throw std::invalid_argument(
            "white_test: insufficient observations (N must exceed q)");
    }

    std::vector<std::vector<Real>> design(n, std::vector<Real>(q));
    std::vector<Real> e2(n);
    for (Size i = 0; i < n; ++i) {
        Size col = 0;
        design[i][col++] = 1.0;  // 常数项
        for (Size j = 0; j < k; ++j) {
            design[i][col++] = X[i][j];  // 原始 X
        }
        if (include_cross_terms && k >= 2) {
            for (Size j = 0; j < k; ++j) {
                for (Size l = j + 1; l < k; ++l) {
                    design[i][col++] = X[i][j] * X[i][l];  // 交叉项
                }
            }
        }
        for (Size j = 0; j < k; ++j) {
            design[i][col++] = X[i][j] * X[i][j];  // 平方项
        }
        e2[i] = residuals[i] * residuals[i];
    }

    std::vector<Real> fitted, resid;
    Real r2;
    detail::ols_simple(e2, design, fitted, resid, r2);

    // LM = N * R²_aux ~ χ²(q-1) (减1因常数列不计)
    const Real lm = static_cast<Real>(n) * r2;
    const Real df = static_cast<Real>(q - 1);
    const Real p_value = detail::chi2_sf(df, lm);

    WhiteResult result;
    result.base.statistic = lm;
    result.base.p_value = p_value;
    result.base.method_name = "White";
    result.base.reject_null = (p_value < 0.05);
    return result;
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
