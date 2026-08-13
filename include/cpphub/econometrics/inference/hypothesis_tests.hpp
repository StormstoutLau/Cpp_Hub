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
// detail: 假设检验辅助函数
//   chi2_sf/gammp/gammq/gser/gcf/gammln/chi2_cdf/chi2_ppf_approx/f_sf/beta_i
//   已移至 core/special_functions.hpp (Phase 7A: 供 ADR-015 方案 B 头文件复用, 不引入 Eigen3)
// =============================================================================
namespace detail {

// betacf / beta_i / gammln / gser / gcf / gammp / gammq / chi2_sf / chi2_cdf
// / chi2_ppf_approx / f_sf 均定义在 core/special_functions.hpp

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
