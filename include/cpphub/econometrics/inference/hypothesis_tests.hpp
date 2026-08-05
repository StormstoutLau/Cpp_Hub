// SOURCE: PHASE6_IMPLEMENTATION_PLAN §4.2 任务 2.3 - 假设检验 (Wald/LR/LM)
// 教材锚点:
//   - Greene 8ed Ch.5 (Wald/LR/LM 显式公式)
//   - Wooldridge CS 2ed Ch.12 (渐近等价性证明)
//   - Engle 1984 (Wald/LR/LM 对偶关系)
//
// 公式 (Greene 8ed §5.3):
//   Wald = (Rβ̂ - r)' [R V R']^{-1} (Rβ̂ - r)  ~ χ²(q)
//   LR   = 2 (ℓ_UR - ℓ_R)                       ~ χ²(q)
//   LM   = N · R²_aux                           ~ χ²(q)  (Breusch-Pagan form)
//          等价于 score' · I^{-1} · score
//
// 排幻觉点 E9: R `lmtest::waldtest` 默认 F 检验 (小样本), C++ 同时提供 χ² 和 F
//   F = Wald / q ~ F(q, df_residual), 小样本更稳健
//
// 排幻觉点 H1: LM = N·R²_aux (非 (N-K)·R²_aux), Breusch-Pagan 1979 原始形式
// 排幻觉点 H2: J-test df = q - k (过度识别约束数), 非 q (矩条件数)
//   q = 矩条件数, k = 参数数; 严格识别 q=k 时 df=0, J 恒为 0
//
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训)
#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/estimation_result.hpp"
#include "cpphub/econometrics/core/special_functions.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

using linalg::dynamic::MatrixXD;
using linalg::dynamic::VectorXD;

// =============================================================================
// detail: 不完全 Gamma 函数 (Numerical Recipes gammp/gammq)
//   用于 χ² 分布 p 值计算: P(χ² > x) = gammq(df/2, x/2)
//
// 注: gammln 改用 std::lgamma (C++标准, 精度更高), 替代 Lanczos 近似
//   原 Lanczos 近似对小参数 (如 0.5) 有 1e-4 级误差, std::lgamma 达机器精度
// =============================================================================
namespace detail {

constexpr Real kGammaEuler = 0.57721566490153286060651209008240243104215933593992;

/// @brief 对数 Gamma 函数 (委托 std::lgamma, C++标准保证精度)
inline Real gammln(Real xx) {
    return std::lgamma(xx);
}

// 不完全 Gamma 函数 P(a,x) = γ(a,x)/Γ(a) (级数展开, 适用于 x < a+1)
inline Real gser(Real a, Real x) {
    const Real EPS = 3e-16;
    const int ITMAX = 300;
    if (x <= 0.0) return 0.0;
    Real gln = gammln(a);
    Real ap = a;
    Real sum = 1.0 / a;
    Real del = sum;
    for (int n = 0; n < ITMAX; ++n) {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (std::abs(del) < std::abs(sum) * EPS) break;
    }
    return sum * std::exp(-x + a * std::log(x) - gln);
}

// 不完全 Gamma 函数 Q(a,x) = 1 - P(a,x) (连分式展开, 适用于 x >= a+1)
inline Real gcf(Real a, Real x) {
    const Real EPS = 3e-16;
    const Real FPMIN = 1e-300;
    const int ITMAX = 300;
    Real gln = gammln(a);
    Real b = x + 1.0 - a;
    Real c = 1.0 / FPMIN;
    Real d = 1.0 / b;
    Real h = d;
    for (int i = 1; i <= ITMAX; ++i) {
        const Real an = -static_cast<Real>(i) * (static_cast<Real>(i) - a);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        const Real del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < EPS) break;
    }
    return std::exp(-x + a * std::log(x) - gln) * h;
}

/// @brief 不完全 Gamma 函数 P(a,x) = γ(a,x)/Γ(a)
inline Real gammp(Real a, Real x) {
    if (x < 0.0 || a <= 0.0) {
        throw std::domain_error("gammp: invalid arguments (x<0 or a<=0)");
    }
    if (x < a + 1.0) {
        return gser(a, x);
    }
    return 1.0 - gcf(a, x);
}

/// @brief 不完全 Gamma 函数 Q(a,x) = 1 - P(a,x) (上尾)
inline Real gammq(Real a, Real x) {
    if (x < 0.0 || a <= 0.0) {
        throw std::domain_error("gammq: invalid arguments (x<0 or a<=0)");
    }
    if (x < a + 1.0) {
        return 1.0 - gser(a, x);
    }
    return gcf(a, x);
}

/// @brief χ² 分布上尾概率 P(χ²_df > x)
///   通用形式: gammq(df/2, x/2)
///   df=1 精确形式: P(χ²(1) > x) = erfc(√(x/2))  (因 χ²(1) = Z², P(|Z|>√x) = erfc(√x/√2))
///   df=2 精确形式: P(χ²(2) > x) = exp(-x/2)     (因 χ²(2) 是均值 2 的指数分布)
///
/// 排幻觉点: df=1/2 使用 std::erfc/std::exp (机器精度), 避免 gammq 级数/连分式
///   在极端尾部 (如 x=14, p≈1.8e-4) 的精度损失. SciPy 同样对 df∈{1,2} 特殊处理.
inline Real chi2_sf(Real df, Real x) {
    if (df <= 0.0) return std::numeric_limits<Real>::quiet_NaN();
    if (x <= 0.0) return 1.0;
    // df=1 精确: χ²(1) = Z², P(χ²(1)>x) = erfc(√(x/2))
    if (df == 1.0) {
        return std::erfc(std::sqrt(0.5 * x));
    }
    // df=2 精确: χ²(2) ~ Exp(1/2), P(χ²(2)>x) = exp(-x/2)
    if (df == 2.0) {
        return std::exp(-0.5 * x);
    }
    return gammq(0.5 * df, 0.5 * x);
}

/// @brief χ² 分布下尾概率 P(χ²_df ≤ x) = 1 - chi2_sf(df, x)
inline Real chi2_cdf(Real df, Real x) {
    if (df <= 0.0) return std::numeric_limits<Real>::quiet_NaN();
    if (x <= 0.0) return 0.0;
    return 1.0 - chi2_sf(df, x);
}

/// @brief χ² 分布分位函数 (临界值, Wilson-Hilferty 近似)
///   仅用于临界值显示, 不用于精确推断
inline Real chi2_ppf_approx(Real df, Real p) {
    // Wilson-Hilferty 1931: χ²_p ≈ df·(1 - 2/(9df) + z_p·√(2/(9df)))³
    // z_p 为标准正态分位数
    // 简化: 用查表 + 线性插值 (p ∈ {0.95, 0.99})
    (void)p;
    // 对 p=0.95 和 p=0.99 用 Wilson-Hilferty
    // z_0.95 = 1.6449, z_0.99 = 2.3263
    const Real z = (p > 0.975) ? 2.3263 : 1.6449;
    const Real t = 1.0 - 2.0 / (9.0 * df) + z * std::sqrt(2.0 / (9.0 * df));
    return df * t * t * t;
}

// betacf / beta_i 已移至 core/special_functions.hpp (共享定义点, 消除跨头文件重复)

/// @brief F 分布上尾概率 P(F_{df1,df2} > f)
///   F CDF: P(F ≤ f) = I_{x'}(df1/2, df2/2), x' = df1·f / (df1·f + df2)
///   F SF  = 1 - CDF = I_{1-x'}(df2/2, df1/2) = I_x(df2/2, df1/2)
///   其中 x = df2 / (df2 + df1·f) = 1 - x'
///
/// 排幻觉点: 参数顺序为 (df2/2, df1/2), 不是 (df1/2, df2/2)
///   I_x(a,b) ≠ I_x(b,a) — 正则化不完全贝塔函数参数不对称
inline Real f_sf(Real df1, Real df2, Real f) {
    if (df1 <= 0.0 || df2 <= 0.0) return std::numeric_limits<Real>::quiet_NaN();
    if (f <= 0.0) return 1.0;
    const Real x = df2 / (df2 + df1 * f);
    return beta_i(0.5 * df2, 0.5 * df1, x);  // 排幻觉: (df2/2, df1/2) 非对称
}

}  // namespace detail

// =============================================================================
// HypothesisTestResult - 假设检验结果
// =============================================================================
struct HypothesisTestResult {
    std::string test_name;             ///< "Wald" / "LR" / "LM" / "J"
    Real statistic = 0.0;              ///< 检验统计量
    Real p_value = 0.0;                ///< p 值 (χ² 或 F 上尾)
    Size df = 0;                       ///< 自由度 (约束个数 q 或 q-k)
    Real critical_value_95 = 0.0;      ///< 5% 临界值 (近似)
    Real critical_value_99 = 0.0;      ///< 1% 临界值 (近似)
    bool reject_null_95 = false;       ///< 5% 显著性下是否拒绝 H0
    bool reject_null_99 = false;       ///< 1% 显著性下是否拒绝 H0
    bool use_f_distribution = false;   ///< 是否使用 F 分布 (小样本修正)
    Real f_df1 = 0.0;                  ///< F 分布分子自由度
    Real f_df2 = 0.0;                  ///< F 分布分母自由度
};

// =============================================================================
// Wald 检验 (线性约束 Rβ = r)
//   H0: Rβ = r,  H1: Rβ ≠ r
//   Wald = (Rβ̂ - r)' [R V R']^{-1} (Rβ̂ - r)  ~ χ²(q)
//   F 形式: F = Wald / q ~ F(q, df_residual)  (小样本修正)
//
// 排幻觉点 E9: 同时返回 χ² 和 F 两种 p 值
//   R `lmtest::waldtest` 默认 F (小样本), C++ 主调方按需选择
// =============================================================================
inline HypothesisTestResult wald_test(const VectorXD& beta,
                                       const MatrixXD& vcov,
                                       const MatrixXD& R,
                                       const VectorXD& r,
                                       Real df_residual = -1.0) {
    const Size k = beta.size();
    const Size q = R.rows();  // 约束个数

    if (R.cols() != k) {
        throw std::invalid_argument("wald_test: R.cols() must equal beta.size()");
    }
    if (r.size() != q) {
        throw std::invalid_argument("wald_test: r.size() must equal R.rows()");
    }
    if (vcov.rows() != k || vcov.cols() != k) {
        throw std::invalid_argument("wald_test: vcov must be k×k");
    }
    if (q == 0) {
        throw std::invalid_argument("wald_test: q (number of restrictions) must be > 0");
    }

    // Rβ̂ - r
    const Eigen::VectorXd Rb_minus_r = R.eigen() * beta.eigen() - r.eigen();

    // R V R'  (q×q)
    const Eigen::MatrixXd RVR = R.eigen() * vcov.eigen() * R.eigen().transpose();

    // (R V R')^{-1} via LLT (假设正定)
    Eigen::LLT<Eigen::MatrixXd> llt(RVR);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("wald_test: RVR' is not positive definite (LLT failed)");
    }

    // Wald = (Rβ̂-r)' (RVR')^{-1} (Rβ̂-r)
    const Eigen::VectorXd solved = llt.solve(Rb_minus_r);
    const Real wald = static_cast<Real>(Rb_minus_r.dot(solved));

    HypothesisTestResult result;
    result.test_name = "Wald";
    result.statistic = wald;
    result.df = q;

    // χ² p 值
    const Real chi2_pvalue = detail::chi2_sf(static_cast<Real>(q), wald);

    // F p 值 (小样本修正): F = Wald/q ~ F(q, df_residual)
    Real f_pvalue = std::numeric_limits<Real>::quiet_NaN();
    if (df_residual > 0.0) {
        const Real f_stat = wald / static_cast<Real>(q);
        f_pvalue = detail::f_sf(static_cast<Real>(q), df_residual, f_stat);
        result.use_f_distribution = true;
        result.f_df1 = static_cast<Real>(q);
        result.f_df2 = df_residual;
        result.p_value = f_pvalue;  // 默认用 F (与 R 一致, 更保守)
    } else {
        result.p_value = chi2_pvalue;  // 无 df_residual 时用 χ²
    }

    // 临界值 (近似)
    result.critical_value_95 = detail::chi2_ppf_approx(static_cast<Real>(q), 0.95);
    result.critical_value_99 = detail::chi2_ppf_approx(static_cast<Real>(q), 0.99);
    result.reject_null_95 = (result.p_value < 0.05);
    result.reject_null_99 = (result.p_value < 0.01);

    return result;
}

// =============================================================================
// LR 检验 (似然比)
//   H0: 约束成立
//   LR = 2 (ℓ_UR - ℓ_R)  ~ χ²(q)
//   q = n_params_UR - n_params_R (约束个数)
// =============================================================================
inline HypothesisTestResult lr_test(const EstimationResult& unrestricted,
                                     const EstimationResult& restricted) {
    const Real ll_ur = unrestricted.log_likelihood;
    const Real ll_r = restricted.log_likelihood;

    if (ll_ur < ll_r - 1e-10) {
        // 数值误差容许微小负值, 但显著负值是错误
        throw std::invalid_argument(
            "lr_test: unrestricted log_likelihood < restricted (should be >=)");
    }

    const Real lr = 2.0 * (ll_ur - ll_r);
    const Size q = unrestricted.n_params - restricted.n_params;

    if (restricted.n_params >= unrestricted.n_params) {
        throw std::invalid_argument(
            "lr_test: restricted.n_params must be < unrestricted.n_params");
    }

    HypothesisTestResult result;
    result.test_name = "LR";
    result.statistic = lr;
    result.df = q;
    result.p_value = detail::chi2_sf(static_cast<Real>(q), lr);
    result.critical_value_95 = detail::chi2_ppf_approx(static_cast<Real>(q), 0.95);
    result.critical_value_99 = detail::chi2_ppf_approx(static_cast<Real>(q), 0.99);
    result.reject_null_95 = (result.p_value < 0.05);
    result.reject_null_99 = (result.p_value < 0.01);

    return result;
}

// =============================================================================
// LM 检验 (score test, Breusch-Pagan 1979 形式)
//   H0: 约束成立 (使用有约束估计的残差)
//   LM = N · R²_aux  ~ χ²(q)
//   其中 R²_aux 是有约束残差对无约束 X 回归的 R²
//
// 排幻觉点 H1: LM = N·R²_aux (非 (N-K)·R²_aux)
//   Breusch-Pagan 1979 原始形式, 大样本下 N 与 N-K 渐近等价
//
// 矩阵形式: LM = ε_R' X (X'X)^{-1} X' ε_R / σ²_R
//   其中 σ²_R = ε_R'ε_R / N (MLE), 化简后 LM = N · R²_aux
// =============================================================================
inline HypothesisTestResult lm_test(const MatrixXD& X,
                                     const VectorXD& residuals_restricted,
                                     const MatrixXD& XtX_inv) {
    const Size N = X.rows();
    const Size K = X.cols();

    if (residuals_restricted.size() != N) {
        throw std::invalid_argument("lm_test: residuals_restricted.size() != X.rows()");
    }
    if (XtX_inv.rows() != K || XtX_inv.cols() != K) {
        throw std::invalid_argument("lm_test: XtX_inv must be K×K");
    }

    // σ²_R = ε_R'ε_R / N (MLE)
    const Real ssr = residuals_restricted.eigen().squaredNorm();
    if (ssr <= 0.0) {
        throw std::runtime_error("lm_test: zero residual sum of squares (perfect fit)");
    }
    const Real sigma2_r = ssr / static_cast<Real>(N);

    // score = X' ε_R  (K×1)
    const Eigen::VectorXd score = X.eigen().transpose() * residuals_restricted.eigen();

    // LM = score' (X'X)^{-1} score / σ²_R = ε_R'X(X'X)^{-1}X'ε_R / σ²_R
    const Eigen::VectorXd tmp = XtX_inv.eigen() * score;
    const Real lm = static_cast<Real>(score.dot(tmp)) / sigma2_r;

    // q = K (无约束模型参数数, 约束模型参数为 0 的极端情况)
    // 注: 实际 q 应由调用方根据约束个数确定, 这里默认 q = K
    const Size q = K;

    HypothesisTestResult result;
    result.test_name = "LM";
    result.statistic = lm;
    result.df = q;
    result.p_value = detail::chi2_sf(static_cast<Real>(q), lm);
    result.critical_value_95 = detail::chi2_ppf_approx(static_cast<Real>(q), 0.95);
    result.critical_value_99 = detail::chi2_ppf_approx(static_cast<Real>(q), 0.99);
    result.reject_null_95 = (result.p_value < 0.05);
    result.reject_null_99 = (result.p_value < 0.01);

    return result;
}

// =============================================================================
// 过度识别检验 (Hansen J-test, Hansen 1982)
//   J = n · ḡ' Ŝ^{-1} ḡ  ~ χ²(q - k)
//   q = 矩条件数, k = 参数数
//   严格识别 (q=k): df=0, J 恒为 0 (无过度识别约束)
//
// 排幻觉点 H2: J-test df = q - k (过度识别约束数), 非 q
//
// 注: moments 假设为 ḡ (均值矩条件), weighting_matrix = n·Ŝ^{-1} (已含 n)
//   则 J = moments' · weighting_matrix · moments
// =============================================================================
inline HypothesisTestResult overidentification_test(const VectorXD& moments,
                                                     const MatrixXD& weighting_matrix,
                                                     Size n_params) {
    const Size q = moments.size();  // 矩条件数
    const Size k = n_params;        // 参数数

    if (weighting_matrix.rows() != q || weighting_matrix.cols() != q) {
        throw std::invalid_argument(
            "overidentification_test: weighting_matrix dimension mismatch");
    }
    if (q < k) {
        throw std::invalid_argument(
            "overidentification_test: q (moments) must be >= k (params)");
    }

    // J = moments' W moments
    const Eigen::VectorXd tmp = weighting_matrix.eigen() * moments.eigen();
    const Real J = static_cast<Real>(moments.eigen().dot(tmp));

    const Size df = q - k;  // 过度识别约束数

    HypothesisTestResult result;
    result.test_name = "J";
    result.statistic = J;
    result.df = df;

    if (df == 0) {
        // 严格识别: J 无分布, p 值无意义
        result.p_value = std::numeric_limits<Real>::quiet_NaN();
        result.critical_value_95 = 0.0;
        result.critical_value_99 = 0.0;
        result.reject_null_95 = false;
        result.reject_null_99 = false;
    } else {
        result.p_value = detail::chi2_sf(static_cast<Real>(df), J);
        result.critical_value_95 = detail::chi2_ppf_approx(static_cast<Real>(df), 0.95);
        result.critical_value_99 = detail::chi2_ppf_approx(static_cast<Real>(df), 0.99);
        result.reject_null_95 = (result.p_value < 0.05);
        result.reject_null_99 = (result.p_value < 0.01);
    }

    return result;
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
