// =============================================================================
// uecm_bounds.hpp - UECM 重参数化 + PSS 边限检验 F/t (spec §2.3)
//
// Phase 8 v1.8 M1 Step 4 (PHASE8_SPEC.md v1.1 §2.3)
//
// 理论锚点: Pesaran-Shin-Smith 2001 (JAE 16(3):289-326) bounds approach
//
// 一手源码语义落档 (statsmodels 0.14.6 tsa/ardl/model.py UECM._construct_
// regressors L1848-1925 / bounds_test L2267-2444 实录; R ARDL 0.2.5
// parse_case/bounds_f_test 交叉, probe_r_ardl.R + verify_ardl.py [S4]
// 逐值自证: B 变换 ≤4.5e-15, F 三方 ≤8.5e-14, t 五组 ≤7.1e-15):
//
//   - UECM 列序 (statsmodels blocks 序冻结): [det(kt), π_y, π_1..π_k,
//     ψ_1..ψ_{p−1}, γ_{1,0}..γ_{1,q_1−1}, ..., γ_{k,·}] — Δx 块把同期
//     差分 (γ_{j,0}=ω_j) 与滞后差分连续排列 (与 spec UECM 形式一致)
//   - B 变换 (不重估, OLS 下与 statsmodels UECM.from_ardl().fit() 等价):
//       π_y = Σφ_i − 1 (仿射, b0 = −1);  π_j = Σ_{l=0..q_j} θ_{j,l}
//       ψ_i = −Σ_{r=i+1..p} φ_r (i=1..p−1)
//       γ_{j,0} = θ_{j,0};  γ_{j,l} = −Σ_{r=l+1..q_j} θ_{j,r} (l ≥ 1)
//     vcov: V_u = B V_a Bᵀ (仿射项不影响)
//   - UECM 样本 = ARDL 同一样本 (重参数化不改估计样本); 校验镜像
//     statsmodels UECM 构造: hold_back < max(p, max q_j, 1) → throw
//     (差分至少损失 1 观测, UECM _maxlag 下限 1, L1909)
//   - F 统计量 (PSS7: Wald 二次型/约束数, statsmodels bounds_test
//     L2414 直译): F = (Rβ)'(RVR')⁻¹(Rβ)/q, 约束行 (UECM 列序上):
//       I:   [π_y, π_j]                    (k+1)
//       II:  [const, π_y, π_j]             (k+2)
//       III: [π_y, π_j]                    (k+1)
//       IV:  [trend, π_y, π_j]             (k+2)
//       V:   [π_y, π_j]                    (k+1)
//   - t 统计量: t = π_y/SE(π_y) (vcov 对角元开方); 仅 Case I/III/V
//     (PSS Table CII 边界, A4; statsmodels 无 t — probe [D]; R ARDL
//     bounds_t_test 交叉 1e-8, 实测 ≤7.1e-15)
//   - t 临界值表 (PSS_T_CRIT): R ARDL 0.2.5 crit_val_bounds_pss2001$t
//     转录 = PSS 2001 Table CII 原文 (两位小数); ⚠️ 档位 {10,5,2.5,1}%
//     ≠ F 表 {10,5,1,0.1}% (PSS9 对偶面: 两表第 3/4 档语义不同) —
//     bounds_t_test 的 alpha 域因此与 F 不同
//   - case-trend 相容 (镜像 R ARDL parse_case throw): I↔"n", II/III↔"c",
//     IV/V↔"ct"; "ctt" 不支持 bounds; 不匹配 → invalid_argument
//   - ⚠️ PSS10: statsmodels bounds_test 查表键 k = nvar = 1+k_PSS (偏移
//     +1, danish 实证: k_PSS=3 → 键 (4,3) CV [2.8772,4.0105] vs 本库键
//     (3,3) [3.2290,4.3223]); 本库 k = force 变量数 (PSS/R ARDL 语义),
//     CV/p 走 pss_critical_values 表键直查 — 回归断言锁死防误镜像
//   - 结论三值 (与 p 值对偶): F: stat > CV_I1(α) → LevelRelation;
//     stat < CV_I0(α) → NoLevelRelation; 否则 Inconclusive (边限检验
//     固有, 必须显式). t: 方向相反 (负侧拒绝)
//   - alpha 限四档 {0.10, 0.05, 0.01, 0.001} (PSS9: statsmodels 分位档
//     (90,95,99,99.9), PSS 原表 (90,95,97.5,99) — 第 3 档是 1% 非 2.5%)
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/ardl/ardl_model.hpp"
#include "cpphub/timeseries/ardl/pss_critical_values.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace ardl {

// PSS 2001 Table CII t 界 (R ARDL 0.2.5 转录, gen_phase8_pss_t_inc.py 生成)
namespace detail {
#include "cpphub/timeseries/ardl/pss_t_bounds.inc"
}  // namespace detail

// 编译期锚 (PSS 2001 原文 Table CII(iii), R ARDL 转录逐位一致):
// case III k=1 5%: I(0) = −2.86, I(1) = −3.22
static_assert(detail::PSS_T_CRIT[1][0][0][1] == -2.86,
              "PSS t bounds anchor (case III, k=1, I0, 5%) = -2.86");
static_assert(detail::PSS_T_CRIT[1][0][1][1] == -3.22,
              "PSS t bounds anchor (case III, k=1, I1, 5%) = -3.22");
// 辅锚: case I k=1 5% I0 = −1.95 (原文 Table CII(i))
static_assert(detail::PSS_T_CRIT[0][0][0][1] == -1.95,
              "PSS t bounds anchor (case I, k=1, I0, 5%) = -1.95");

/// PSS Case 编号 (跨 PSS/statsmodels/R ARDL/pssbounds 一致)
enum class PssCase { I = 1, II = 2, III = 3, IV = 4, V = 5 };

/// 边限检验三值结论
enum class BoundsConclusion {
    NoLevelRelation,   ///< F: stat < I(0) 下界 (t: stat > I(0) 界) — 不能拒绝 H0
    Inconclusive,      ///< 介于两界 (边限检验固有, 必须显式返回)
    LevelRelation      ///< F: stat > I(1) 上界 (t: stat < I(1) 界) — 拒绝 H0
};

/// UECM 重参数化结果 (spec §2.3 逐字段)
struct UecmResult {
    std::vector<Real> params;          ///< [det, π_y, π_1..π_k, ψ_y, ψ_j(=γ)] (列序见头注释)
    std::vector<std::vector<Real>> vcov;  ///< B V_a Bᵀ, 同列序
    std::vector<Real> pi;              ///< (π_y, π_1..π_k) 水平项系数 — 检验对象
    std::vector<std::vector<Real>> pi_vcov;  ///< π 子块 vcov
    Size n_obs_used = 0;               ///< = ARDL 同一样本 (T − hold_back_used)
    Real sigma2 = 0.0;                 ///< = ardl.sigma2 (SSR 与 n_acct 均不变)
    ArdlSpec source_spec{};            ///< 源 ARDL 设定回显
};

namespace detail {

/// case → 期望 trend ("ctt" 无对应 case, 由调用方拒绝)
inline const char* case_trend(PssCase c) {
    switch (c) {
        case PssCase::I: return "n";
        case PssCase::II: return "c";
        case PssCase::III: return "c";
        case PssCase::IV: return "ct";
        case PssCase::V: return "ct";
    }
    throw std::invalid_argument("bounds: case must be I..V");
}

/// alpha → F 表四档索引 {0.10, 0.05, 0.01, 0.001} → {0,1,2,3} (PSS9)
inline Size alpha_index(Real alpha) {
    if (std::abs(alpha - 0.10) < 1e-12) return 0;
    if (std::abs(alpha - 0.05) < 1e-12) return 1;
    if (std::abs(alpha - 0.01) < 1e-12) return 2;
    if (std::abs(alpha - 0.001) < 1e-12) return 3;
    throw std::invalid_argument(
        "bounds test: alpha must be one of {0.10, 0.05, 0.01, 0.001} "
        "(statsmodels percentile grid; no interpolation semantics)");
}

/// alpha → t 表四档索引 {0.10, 0.05, 0.025, 0.01} → {0,1,2,3}
/// (PSS 2001 Table CII 档位 — ≠ F 表档, PSS9 对偶面)
inline Size t_alpha_index(Real alpha) {
    if (std::abs(alpha - 0.10) < 1e-12) return 0;
    if (std::abs(alpha - 0.05) < 1e-12) return 1;
    if (std::abs(alpha - 0.025) < 1e-12) return 2;
    if (std::abs(alpha - 0.01) < 1e-12) return 3;
    throw std::invalid_argument(
        "bounds_t_test: alpha must be one of {0.10, 0.05, 0.025, 0.01} "
        "(PSS 2001 Table CII grid — note: NOT the F-test grid, whose 3rd/4th "
        "levels are 1%/0.1% while the t-table's are 2.5%/1%)");
}

/// PSS 2001 Table CII t 临界值: {I(0) 四档, I(1) 四档} (均为负, I1 更负)
/// 仅 case I/III/V (调用方保证); k ∈ [1,10]
inline std::array<std::array<Real, 4>, 2> t_crit_rows(Size k, PssCase c) {
    Size ci;
    switch (c) {
        case PssCase::I:   ci = 0; break;
        case PssCase::III: ci = 1; break;
        case PssCase::V:   ci = 2; break;
        default:
            throw std::invalid_argument(
                "t_crit_rows: t bounds table covers Case I/III/V only");
    }
    if (k < 1 || k > 10) {
        throw std::invalid_argument(
            "t_crit_rows: k must be in [1, 10] (PSS Table CII coverage)");
    }
    std::array<std::array<Real, 4>, 2> out{};
    for (Size i1 = 0; i1 < 2; ++i1)
        for (Size a = 0; a < 4; ++a)
            out[i1][a] = PSS_T_CRIT[ci][k - 1][i1][a];
    return out;
}

}  // namespace detail

/// ARDL → UECM 系数 B 变换 (不重估; 镜像 statsmodels UECM.from_ardl 数学语义)
///
/// \throw std::invalid_argument hold_back < max(p, max q_j, 1) (差分样本
///        不足, 镜像 statsmodels UECM 构造校验) 或 ARDL(0,·) k=0 退化
inline UecmResult uecm_from_ardl(const ArdlResult& ardl) {
    const ArdlSpec& sp = ardl.spec;
    const Size p = sp.p;
    const Size k = static_cast<Size>(sp.q.size());
    if (p == 0 && k == 0) {
        throw std::invalid_argument(
            "uecm_from_ardl: degenerate ARDL(0,0) has no UECM");
    }
    Size maxq = 0;
    for (Size qj : sp.q) maxq = std::max(maxq, qj);
    if (ardl.hold_back_used < std::max(std::max(p, maxq), Size{1})) {
        throw std::invalid_argument(
            "uecm_from_ardl: ARDL hold_back < max(p, max q_j, 1) — UECM "
            "requires one extra observation for differencing");
    }
    const Size kt = detail::trend_order(sp.trend);

    const Size n_a = ardl.params.size();
    Size n_u = kt + 1 + k + (p >= 1 ? p - 1 : 0);
    for (Size qj : sp.q) n_u += qj;
    if (n_u != n_a) {
        throw std::invalid_argument(
            "uecm_from_ardl: internal layout mismatch (n_u != n_a)");
    }

    // ---- B 变换 (verify_ardl.py manual_uecm_from_ardl 逐行镜像) ----
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n_u, n_a);
    Eigen::VectorXd b0 = Eigen::VectorXd::Zero(n_u);
    std::vector<Size> off(k);
    Size cur = kt + p;
    for (Size j = 0; j < k; ++j) {
        off[j] = cur;
        cur += sp.q[j] + 1;
    }
    for (Size i = 0; i < kt; ++i) B(i, i) = 1.0;        // det 恒等
    for (Size i = 0; i < p; ++i) B(kt, kt + i) = 1.0;   // π_y = Σφ − 1
    b0(kt) = -1.0;
    for (Size j = 0; j < k; ++j)                         // π_j = Σθ_{j,·}
        for (Size l = 0; l <= sp.q[j]; ++l) B(kt + 1 + j, off[j] + l) = 1.0;
    for (Size i = 1; i < p; ++i)                         // ψ_i = −Σ_{r>i} φ_r
        for (Size r = i + 1; r <= p; ++r)
            B(kt + 1 + k + (i - 1), kt + r - 1) = -1.0;
    const Size gam0 = kt + 1 + k + (p >= 1 ? p - 1 : 0);  // γ 块起始行
    Size gbase = gam0;
    for (Size j = 0; j < k; ++j) {                       // γ_{j,l}
        for (Size l = 0; l < sp.q[j]; ++l) {
            const Size row = gbase + l;
            if (l == 0) {
                B(row, off[j]) = 1.0;                    // γ_{j,0} = θ_{j,0}
            } else {
                for (Size r = l + 1; r <= sp.q[j]; ++r)
                    B(row, off[j] + r) = -1.0;  // −Σ_{r>l} θ_{j,r}
            }
        }
        gbase += sp.q[j];
    }

    Eigen::Map<const Eigen::VectorXd> pa(ardl.params.data(),
                                         static_cast<Eigen::Index>(n_a));
    Eigen::MatrixXd Va = Eigen::MatrixXd::Zero(n_a, n_a);
    for (Size i = 0; i < n_a; ++i)
        for (Size j = 0; j < n_a; ++j) Va(i, j) = ardl.vcov[i][j];
    const Eigen::VectorXd pu = B * pa + b0;
    const Eigen::MatrixXd Vu = B * Va * B.transpose();

    UecmResult u;
    u.params.resize(n_u);
    u.vcov.assign(n_u, std::vector<Real>(n_u, 0.0));
    for (Size i = 0; i < n_u; ++i) {
        u.params[i] = pu(i);
        for (Size j = 0; j < n_u; ++j) u.vcov[i][j] = Vu(i, j);
    }
    u.pi.assign(k + 1, 0.0);
    u.pi_vcov.assign(k + 1, std::vector<Real>(k + 1, 0.0));
    for (Size i = 0; i <= k; ++i) {
        u.pi[i] = u.params[kt + i];  // π 块 = params[kt .. kt+k]
        for (Size j = 0; j <= k; ++j)
            u.pi_vcov[i][j] = u.vcov[kt + i][kt + j];
    }
    u.n_obs_used = ardl.n_obs_used;
    u.sigma2 = ardl.sigma2;
    u.source_spec = ardl.spec;
    return u;
}

/// 边限检验结果 (spec §2.3 逐字段)
struct BoundsTestResult {
    Real stat = 0.0;                   ///< F (bounds_f_test) 或 t (bounds_t_test)
    bool is_t_test = false;
    PssCase case_used{};
    Size k = 0;                        ///< force 变量数 (PSS6/PSS10)
    Size n_restrictions = 0;           ///< F 约束数: I/III/V = k+1, II/IV = k+2
    std::array<std::array<Real, 4>, 2> crit_vals{};  ///< [0]=I(0) 行 [1]=I(1) 行 × {10,5,1,0.1}% (PSS9)
    std::array<Real, 2> p_values{};    ///< F: {I(0) null, I(1) null} 响应面近似;
                                       ///< t: {NaN, NaN} + has_p=false (无 t 响应面, A5)
    BoundsConclusion conclusion{};
    bool has_p = true;                 ///< false = t 检验 (无响应面 p 值)
    std::string message;
};

namespace detail {

/// case-trend 相容校验 + k 域校验 (bounds_f/t_test 公共前置)
inline void bounds_check_compat(const UecmResult& u, PssCase c) {
    const std::string& tr = u.source_spec.trend;
    if (tr == "ctt") {
        throw std::invalid_argument(
            "bounds test: trend='ctt' has no PSS case mapping "
            "(PSS tables cover n/c/ct only)");
    }
    const char* want = case_trend(c);
    if (tr != want) {
        throw std::invalid_argument(
            "bounds test: case " + std::to_string(static_cast<int>(c)) +
            " requires trend='" + want + "' but the underlying ARDL has "
            "trend='" + tr + "' (mirror of R ARDL parse_case)");
    }
    const Size k = static_cast<Size>(u.source_spec.q.size());
    if (k < 1 || k > 10) {
        throw std::invalid_argument(
            "bounds test: k (force variables) must be in [1, 10] "
            "(PSS response surface coverage)");
    }
}

/// 约束行索引 (UECM params 列序; statsmodels bounds_test rest 直译):
/// I: [0..k]; II: [0..k+1]; III: [1..k+1]; IV: [1..k+2]; V: [2..k+2]
inline std::vector<Size> bounds_restriction_rows(PssCase c, Size k) {
    const Size nvar = 1 + k;  // π_y + π_j (statsmodels nvar 语义)
    Size lo, hi;
    switch (c) {
        case PssCase::I:   lo = 0; hi = nvar;     break;
        case PssCase::II:  lo = 0; hi = nvar + 1; break;
        case PssCase::III: lo = 1; hi = nvar + 1; break;
        case PssCase::IV:  lo = 1; hi = nvar + 2; break;
        case PssCase::V:   lo = 2; hi = nvar + 2; break;
        default: throw std::invalid_argument("bounds: case must be I..V");
    }
    std::vector<Size> rows;
    for (Size i = lo; i < hi; ++i) rows.push_back(i);
    return rows;
}

}  // namespace detail

/// PSS 边限 F 检验 (Wald 二次型/约束数, PSS7)
///
/// F > CV_I1(α) → LevelRelation; F < CV_I0(α) → NoLevelRelation;
/// 否则 Inconclusive. CV/p 查表键 k = force 变量数 (PSS10: 不镜像
/// statsmodels bounds_test 的 +1 偏移)
///
/// \throw std::invalid_argument case-trend 不相容 / k ∉ [1,10] /
///        alpha ∉ 四档 / 约束子矩阵奇异 (vcov 退化)
inline BoundsTestResult bounds_f_test(const UecmResult& uecm, PssCase c,
                                      Real alpha = 0.05) {
    detail::bounds_check_compat(uecm, c);
    const Size k = static_cast<Size>(uecm.source_spec.q.size());
    const Size ai = detail::alpha_index(alpha);

    const std::vector<Size> rows = detail::bounds_restriction_rows(c, k);
    const Size q = static_cast<Size>(rows.size());
    const Size n = uecm.params.size();
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(q, n);
    for (Size i = 0; i < q; ++i) R(i, rows[i]) = 1.0;
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(n);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, n);
    for (Size i = 0; i < n; ++i) {
        beta(i) = uecm.params[i];
        for (Size j = 0; j < n; ++j) V(i, j) = uecm.vcov[i][j];
    }
    const Eigen::VectorXd rb = R * beta;
    const Eigen::MatrixXd rvr = R * V * R.transpose();
    const Eigen::LLT<Eigen::MatrixXd> llt(rvr);
    if (llt.info() != Eigen::Success) {
        throw std::invalid_argument(
            "bounds_f_test: restriction covariance not positive definite");
    }
    const Real wald = rb.dot(llt.solve(rb));
    const Real f = wald / static_cast<Real>(q);

    BoundsTestResult r;
    r.stat = f;
    r.is_t_test = false;
    r.case_used = c;
    r.k = k;
    r.n_restrictions = q;
    const auto lo = pss_crit_values(k, static_cast<Size>(c), false);
    const auto up = pss_crit_values(k, static_cast<Size>(c), true);
    for (Size i = 0; i < 4; ++i) {
        r.crit_vals[0][i] = lo[i];
        r.crit_vals[1][i] = up[i];
    }
    r.p_values = pss_pvalue_bounds(f, k, static_cast<Size>(c));
    r.has_p = true;
    if (f > up[ai]) {
        r.conclusion = BoundsConclusion::LevelRelation;
    } else if (f < lo[ai]) {
        r.conclusion = BoundsConclusion::NoLevelRelation;
    } else {
        r.conclusion = BoundsConclusion::Inconclusive;
    }
    r.message = "Bounds F-test (Wald) for no cointegration";
    return r;
}

/// PSS 边限 t 检验 (仅 Case I/III/V, PSS Table CII 边界, A4)
///
/// t = π_y/SE(π_y); 方向与 F 相反 (负侧拒绝):
/// t < CV_I1(α) (更负界) → LevelRelation; t > CV_I0(α) → NoLevelRelation;
/// 否则 Inconclusive. 无 p 值 (statsmodels 无 t 响应面路径, probe [D];
/// R ARDL 的 t p 值系 MC 自算, A5 不做)
///
/// \throw std::invalid_argument Case II/IV / case-trend 不相容 /
///        k ∉ [1,10] / alpha ∉ 四档
inline BoundsTestResult bounds_t_test(const UecmResult& uecm, PssCase c,
                                      Real alpha = 0.05) {
    if (c == PssCase::II || c == PssCase::IV) {
        throw std::invalid_argument(
            "bounds_t_test: t-test not available for Case II/IV "
            "(PSS Table CII covers Case I/III/V only)");
    }
    detail::bounds_check_compat(uecm, c);
    const Size k = static_cast<Size>(uecm.source_spec.q.size());
    const Size ai = detail::t_alpha_index(alpha);  // t 表档 {10,5,2.5,1}%

    const Real se = std::sqrt(uecm.pi_vcov[0][0]);
    if (!(se > 0.0)) {
        throw std::invalid_argument("bounds_t_test: SE(pi_y) is zero/NaN");
    }
    const Real t = uecm.pi[0] / se;

    BoundsTestResult r;
    r.stat = t;
    r.is_t_test = true;
    r.case_used = c;
    r.k = k;
    r.n_restrictions = 1;
    // t CV = PSS 2001 Table CII (R ARDL 转录), ≠ F 表 (档位/数值均不同)
    const auto cv = detail::t_crit_rows(k, c);
    for (Size i = 0; i < 4; ++i) {
        r.crit_vals[0][i] = cv[0][i];  // I(0) 下界 (较不负)
        r.crit_vals[1][i] = cv[1][i];  // I(1) 上界 (更负)
    }
    r.p_values = {std::numeric_limits<Real>::quiet_NaN(),
                  std::numeric_limits<Real>::quiet_NaN()};
    r.has_p = false;
    if (t < cv[1][ai]) {       // 负侧: 比更负的 I(1) 界还负 → 拒绝
        r.conclusion = BoundsConclusion::LevelRelation;
    } else if (t > cv[0][ai]) {
        r.conclusion = BoundsConclusion::NoLevelRelation;
    } else {
        r.conclusion = BoundsConclusion::Inconclusive;
    }
    r.message = "Bounds t-test for no cointegration";
    return r;
}

}  // namespace ardl
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
