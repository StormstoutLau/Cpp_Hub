// =============================================================================
// unit_root_common.hpp - 单位根检验共享工具 (spec §3.0.1)
//
// Phase 7B v1.6 M2 (PHASE7B_FINANCIAL_TS_SPEC.md)
//
// 教材锚点: Schwert 1989 / Ng-Perron 1995 (lag 选择) / Newey-West 1987 (LRV)
// 对照库: Python arch 8.0.0 arch/unitroot/unitroot.py (U-ADR9)
//
// 幻觉点防护 (spec §6.2):
//   U1: Schwert 规则 ceil(12·(T/100)^0.25) 向上取整 (arch unitroot.py:416,
//       1128 实测), ADF/DFGLS 有上限 min(lag, max((T-1)//2-1-trend_cols, 0)),
//       PP 无上限 (直接赋值), KPSS legacy 无上限
//   U2: 自动方程形式基于趋势显著性
//   U5/U11: 长期方差 Newey-West Bartlett 核 (cov_nw 复刻, arch cov.py)
//   AIC 惩罚项只数 lag 列数 (arch _select_best_ic: crit = -2·llf + 2·arange,
//       不含 y_{t-1}/常数列; llf 用 σ²=SSR/n 的 MLE 高斯似然)
// =============================================================================
#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace unit_root {

// ---------------------------------------------------------------------------
// 内部: 完整 OLS (系数 + 标准误 + 残差) — 正规方程 + Gauss-Jordan 求逆
// 与 econometrics/detail/ols_simple 的差异: 单位根检验需要 bse/SSR/XtX⁻¹
// ---------------------------------------------------------------------------
namespace detail {

struct OlsFull {
    std::vector<Real> beta;       ///< 系数 (K)
    std::vector<Real> bse;        ///< 标准误 (K): √(s²·(X'X)⁻¹_jj), s²=SSR/(n-K)
    std::vector<Real> resid;      ///< 残差 (n)
    Real ssr = 0.0;               ///< 残差平方和
    Size nobs = 0;                ///< 样本量
    Size n_params = 0;            ///< 参数数 K
    Real sigma2 = 0.0;            ///< s² = SSR/(n-K) (df-corrected)
};

// Gauss-Jordan 矩阵求逆 (partial pivoting), A 为 n×n, 原地求逆
inline void invert_matrix(std::vector<std::vector<Real>>& A) {
    const Size n = A.size();
    std::vector<std::vector<Real>> I(n, std::vector<Real>(n, 0.0));
    for (Size i = 0; i < n; ++i) I[i][i] = 1.0;
    for (Size col = 0; col < n; ++col) {
        Size piv = col;
        Real maxv = std::fabs(A[col][col]);
        for (Size r = col + 1; r < n; ++r) {
            if (std::fabs(A[r][col]) > maxv) {
                maxv = std::fabs(A[r][col]);
                piv = r;
            }
        }
        if (maxv < 1e-300) {
            throw std::runtime_error("invert_matrix: singular matrix");
        }
        if (piv != col) {
            std::swap(A[piv], A[col]);
            std::swap(I[piv], I[col]);
        }
        const Real d = A[col][col];
        for (Size c = 0; c < n; ++c) {
            A[col][c] /= d;
            I[col][c] /= d;
        }
        for (Size r = 0; r < n; ++r) {
            if (r == col) continue;
            const Real f = A[r][col];
            if (f == 0.0) continue;
            for (Size c = 0; c < n; ++c) {
                A[r][c] -= f * A[col][c];
                I[r][c] -= f * I[col][c];
            }
        }
    }
    A = std::move(I);
}

// OLS: y (n) = X (n×K) beta + u; 返回系数/标准误/残差/SSR
inline OlsFull ols_fit(const std::vector<Real>& y,
                       const std::vector<std::vector<Real>>& X) {
    const Size n = y.size();
    if (n == 0 || X.empty()) {
        throw std::invalid_argument("ols_fit: empty inputs");
    }
    const Size k = X[0].size();
    if (n <= k) {
        throw std::invalid_argument("ols_fit: insufficient observations");
    }

    // 正规方程
    std::vector<std::vector<Real>> XtX(k, std::vector<Real>(k, 0.0));
    std::vector<Real> Xty(k, 0.0);
    for (Size i = 0; i < k; ++i) {
        for (Size j = i; j < k; ++j) {
            Real s = 0.0;
            for (Size t = 0; t < n; ++t) s += X[t][i] * X[t][j];
            XtX[i][j] = s;
            XtX[j][i] = s;
        }
        Real s = 0.0;
        for (Size t = 0; t < n; ++t) s += X[t][i] * y[t];
        Xty[i] = s;
    }

    // Gauss-Jordan 消元解 beta (增广矩阵)
    std::vector<std::vector<Real>> A = XtX;
    std::vector<Real> b = Xty;
    for (Size col = 0; col < k; ++col) {
        Size piv = col;
        Real maxv = std::fabs(A[col][col]);
        for (Size r = col + 1; r < k; ++r) {
            if (std::fabs(A[r][col]) > maxv) {
                maxv = std::fabs(A[r][col]);
                piv = r;
            }
        }
        if (maxv < 1e-300) {
            throw std::runtime_error("ols_fit: singular design matrix");
        }
        if (piv != col) {
            std::swap(A[piv], A[col]);
            std::swap(b[piv], b[col]);
        }
        const Real akk = A[col][col];
        for (Size r = col + 1; r < k; ++r) {
            const Real f = A[r][col] / akk;
            if (f == 0.0) continue;
            for (Size c = col; c < k; ++c) A[r][c] -= f * A[col][c];
            b[r] -= f * b[col];
        }
    }
    OlsFull res;
    res.beta.resize(k);
    for (Size i = k; i-- > 0;) {
        Real s = b[i];
        for (Size j = i + 1; j < k; ++j) s -= A[i][j] * res.beta[j];
        res.beta[i] = s / A[i][i];
    }

    // 残差与 SSR
    res.resid.resize(n);
    res.ssr = 0.0;
    for (Size t = 0; t < n; ++t) {
        Real pred = 0.0;
        for (Size j = 0; j < k; ++j) pred += res.beta[j] * X[t][j];
        res.resid[t] = y[t] - pred;
        res.ssr += res.resid[t] * res.resid[t];
    }

    // 标准误: s²·(X'X)⁻¹ 对角
    invert_matrix(XtX);
    res.nobs = n;
    res.n_params = k;
    res.sigma2 = res.ssr / static_cast<Real>(n - k);
    res.bse.resize(k);
    for (Size j = 0; j < k; ++j) {
        res.bse[j] = std::sqrt(res.sigma2 * XtX[j][j]);
    }
    return res;
}

// 正态 CDF: Φ(x) = 0.5·erfc(-x/√2)
inline Real normal_cdf(Real x) noexcept {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// trend 字符串 → trend 列数 ("n"/"nc"=0, "c"=1, "ct"=2)
inline Size trend_cols(const std::string& trend) {
    if (trend == "n" || trend == "nc") return 0;
    if (trend == "c") return 1;
    if (trend == "ct") return 2;
    throw std::invalid_argument("trend must be one of n/nc/c/ct");
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Schwert lag 选择规则 (U1, U-ADR1)
// ceil(12·(T/100)^0.25), 向上取整 (非 floor; arch unitroot.py:416 实测)
// ---------------------------------------------------------------------------
inline Size schwert_lag(Size T) {
    if (T == 0) {
        throw std::invalid_argument("schwert_lag: empty sample");
    }
    return static_cast<Size>(
        std::ceil(12.0 * std::pow(static_cast<Real>(T) / 100.0, 0.25)));
}

// ---------------------------------------------------------------------------
// ADF 类回归的 lag 自动选择 (AIC/BIC) — 复刻 arch _df_select_lags +
// _autolag_ols + _select_best_ic (unitroot.py:147-204, 299-371, 372-445)
//
// 关键约定:
//   - 样本固定截断: nobs = (T-1) - max_lag (所有候选 lag 用同一样本)
//   - exog = [trend 列..., y_{t-1}, dy_{t-1}, ..., dy_{t-max_lag}]
//     (arch lagmat trim="both" original="in", rhs[:,0] 替换为 y_{t-1},
//      add_trend prepend=True → trend 列在前)
//   - 候选 p = 0..max_lag: 用前 (trend_cols+1+p) 列回归, σ² = SSR/nobs (MLE)
//   - llf = -nobs/2·(log 2π + log σ² + 1)
//   - AIC: crit = -2·llf + 2·p (惩罚只数 lag 列, 不数 y_{t-1}/trend)
//   - BIC: crit = -2·llf + log(nobs)·p
//   - 选 crit 最小; tie 时取小 p (顺序遍历严格小于)
//
// @param y 原始序列 (水平值)
// @param trend "n"/"c"/"ct"
// @param max_lag 最大候选 lag (0 => Schwert 自动 + arch 上限保护)
// @param criterion "aic" 或 "bic"
// @return 最优 lag 数 p
inline Size select_lag_by_ic(const std::vector<Real>& y,
                             const std::string& trend = "c",
                             Size max_lag = 0,
                             const std::string& criterion = "aic") {
    const Size T = y.size();
    if (T < 3) {
        throw std::invalid_argument("select_lag_by_ic: sample too small");
    }
    bool use_bic = false;
    if (criterion == "aic") {
        use_bic = false;
    } else if (criterion == "bic") {
        use_bic = true;
    } else {
        throw std::invalid_argument("criterion must be aic or bic");
    }

    // arch 上限保护 (unitroot.py:412-417)
    if (max_lag == 0) {
        Size max_max = (T - 1) / 2 > 1 ? (T - 1) / 2 - 1 : 0;
        const Size tc = detail::trend_cols(trend);
        if (tc > max_max) max_max = 0; else max_max -= tc;
        max_lag = schwert_lag(T);
        if (max_lag > max_max) max_lag = max_max;
    }
    if (T - 1 <= max_lag) {
        throw std::invalid_argument("select_lag_by_ic: max_lag too large");
    }

    // 固定截断样本: dy[i] 对应 lhs, rhs = [y_{t-1}, dy_{t-1}, ..., dy_{t-max_lag}]
    const Size nobs = (T - 1) - max_lag;
    std::vector<std::vector<Real>> X(nobs);
    for (Size i = 0; i < nobs; ++i) {
        // 时间索引: t = max_lag + i + 1 (dy_t 的 t), y_{t-1} = y[max_lag+i]
        std::vector<Real> row;
        row.reserve(max_lag + 3);
        // trend 列在前 (arch add_trend prepend=True)
        const Size tc = detail::trend_cols(trend);
        if (tc >= 1) row.push_back(1.0);
        if (tc >= 2) row.push_back(static_cast<Real>(i));  // t 从 0 起
        row.push_back(y[max_lag + i]);                     // y_{t-1}
        for (Size l = 1; l <= max_lag; ++l) {
            row.push_back(y[max_lag + i + 1 - l] - y[max_lag + i - l]);  // dy_{t-l}
        }
        X[i] = std::move(row);
    }
    std::vector<Real> lhs(nobs);
    for (Size i = 0; i < nobs; ++i) {
        lhs[i] = y[max_lag + i + 1] - y[max_lag + i];  // dy_t
    }

    const Size tc = detail::trend_cols(trend);
    const Size base_cols = tc + 1;  // trend + y_{t-1}
    const Real n = static_cast<Real>(nobs);
    const Real log_2pi = std::log(2.0 * 3.14159265358979323846);

    Real best_crit = std::numeric_limits<Real>::infinity();
    Size best_lag = 0;
    for (Size p = 0; p <= max_lag; ++p) {
        // 用前 base_cols + p 列
        std::vector<std::vector<Real>> Xp(nobs);
        for (Size i = 0; i < nobs; ++i) {
            Xp[i].assign(X[i].begin(), X[i].begin() + static_cast<std::ptrdiff_t>(base_cols + p));
        }
        // 正规方程求解 (只用 XtX, 避免 ols_fit 的额外输出)
        const Size k = base_cols + p;
        std::vector<std::vector<Real>> XtX(k, std::vector<Real>(k, 0.0));
        std::vector<Real> Xty(k, 0.0);
        for (Size a = 0; a < k; ++a) {
            for (Size b = a; b < k; ++b) {
                Real s = 0.0;
                for (Size t = 0; t < nobs; ++t) s += Xp[t][a] * Xp[t][b];
                XtX[a][b] = s;
                XtX[b][a] = s;
            }
            Real s = 0.0;
            for (Size t = 0; t < nobs; ++t) s += Xp[t][a] * lhs[t];
            Xty[a] = s;
        }
        // Gauss-Jordan
        for (Size col = 0; col < k; ++col) {
            Size piv = col;
            Real maxv = std::fabs(XtX[col][col]);
            for (Size r = col + 1; r < k; ++r) {
                if (std::fabs(XtX[r][col]) > maxv) {
                    maxv = std::fabs(XtX[r][col]);
                    piv = r;
                }
            }
            if (piv != col) {
                std::swap(XtX[piv], XtX[col]);
                std::swap(Xty[piv], Xty[col]);
            }
            const Real akk = XtX[col][col];
            for (Size r = col + 1; r < k; ++r) {
                const Real f = XtX[r][col] / akk;
                if (f == 0.0) continue;
                for (Size c = col; c < k; ++c) XtX[r][c] -= f * XtX[col][c];
                Xty[r] -= f * Xty[col];
            }
        }
        std::vector<Real> beta(k);
        for (Size i = k; i-- > 0;) {
            Real s = Xty[i];
            for (Size j = i + 1; j < k; ++j) s -= XtX[i][j] * beta[j];
            beta[i] = s / XtX[i][i];
        }
        // SSR (arch: sigma2 = (ypy - b'Xpx b)/n, 代数等价)
        Real ssr = 0.0;
        for (Size t = 0; t < nobs; ++t) {
            Real e = lhs[t];
            for (Size a = 0; a < k; ++a) e -= beta[a] * Xp[t][a];
            ssr += e * e;
        }
        const Real sigma2 = ssr / n;
        const Real llf = -n / 2.0 * (log_2pi + std::log(sigma2) + 1.0);
        const Real crit = use_bic
                              ? -2.0 * llf + std::log(n) * static_cast<Real>(p)
                              : -2.0 * llf + 2.0 * static_cast<Real>(p);
        if (crit < best_crit) {
            best_crit = crit;
            best_lag = p;
        }
    }
    return best_lag;
}

// ---------------------------------------------------------------------------
// 自动方程形式选择 (U2, U-ADR2)
// 基于趋势显著性: ct 回归趋势项 |t|>1.96 → "ct"; c 回归常数 |t|>1.96 → "c";
// 否则 "nc"。 (arch 无对应功能 — arch trend 由用户显式指定; 此处为 spec U2
// 的工程化实现, 测试做方向性验证而非数值对照)
// ---------------------------------------------------------------------------
inline std::string select_trend_spec(const std::vector<Real>& data) {
    const Size T = data.size();
    if (T < 5) {
        throw std::invalid_argument("select_trend_spec: sample too small");
    }
    // dy_t = a + b·t + g·y_{t-1} (ct 形式, 无 lag)
    const Size n = T - 1;
    std::vector<Real> lhs(n);
    std::vector<std::vector<Real>> Xc(n, std::vector<Real>(3));
    for (Size i = 0; i < n; ++i) {
        const Size t = i + 1;  // dy_t 时刻
        lhs[i] = data[t] - data[t - 1];
        Xc[i][0] = 1.0;                     // a
        Xc[i][1] = static_cast<Real>(i);    // b (0-based)
        Xc[i][2] = data[t - 1];             // g
    }
    // ct 回归: y_{t-1} 与 t 强共线时 (近完美趋势) 设计矩阵奇异,
    // 视为趋势显著直接返回 "ct" (趋势方向明确的病态输入)
    try {
        const auto rc = detail::ols_fit(lhs, Xc);
        if (rc.bse[1] > 0.0 &&
            std::fabs(rc.beta[1] / rc.bse[1]) > 1.959963984540054) {
            return "ct";
        }
    } catch (const std::runtime_error&) {
        return "ct";  // 奇异: 完美趋势序列
    }
    // dy_t = a + g·y_{t-1}
    std::vector<std::vector<Real>> Xm(n, std::vector<Real>(2));
    for (Size i = 0; i < n; ++i) {
        Xm[i][0] = 1.0;
        Xm[i][1] = data[i];
    }
    const auto rm = detail::ols_fit(lhs, Xm);
    if (rm.bse[0] > 0.0 && std::fabs(rm.beta[0] / rm.bse[0]) > 1.959963984540054) {
        return "c";
    }
    return "nc";
}

// ---------------------------------------------------------------------------
// 长期方差估计 (Newey-West Bartlett 核, U5/U11/U-ADR7/U-ADR8)
// 复刻 arch utility/cov.py::cov_nw(x, lags, demean=False):
//   γ_j = (1/n)·Σₜ uₜ·uₜ₋ⱼ  (分母恒为 n, 非 n-j)
//   λ = γ₀ + 2·Σ_{j=1}^{L} (1 - j/(L+1))·γ_j
// bandwidth 语义: 显式滞后阶数 (0 = 仅 γ₀, 即 σ²·(n-1)/n)。
// 调用方 (PP/KPSS) 总是先算好带宽再传入 — arch cov_nw 同样要求显式 lags
// 手算基准 (tests/unit/timeseries/unit_root_baseline.inc):
//   COVNW_U5_L0: u=[1..5], L=0 → 11.0
//   COVNW_U5_L2: u=[1..5], L=2 → 25.133333333333336
// ---------------------------------------------------------------------------
inline Real long_run_variance(const std::vector<Real>& residuals,
                              Size bandwidth = 0,
                              const std::string& kernel = "Bartlett") {
    if (kernel != "Bartlett") {
        throw std::invalid_argument(
            "long_run_variance: only Bartlett kernel is supported "
            "(arch cov_nw, U6/U11-kernel)");
    }
    const Size n = residuals.size();
    if (n == 0) {
        throw std::invalid_argument("long_run_variance: empty residuals");
    }
    const Size L = bandwidth;  // 显式带宽 (0 = 仅 γ₀)
    if (L >= n) {
        throw std::invalid_argument("long_run_variance: bandwidth >= n");
    }

    const Real nn = static_cast<Real>(n);
    Real gamma0 = 0.0;
    for (Real u : residuals) gamma0 += u * u;
    gamma0 /= nn;

    Real lam = gamma0;
    for (Size j = 1; j <= L; ++j) {
        Real gj = 0.0;
        for (Size t = j; t < n; ++t) gj += residuals[t] * residuals[t - j];
        gj /= nn;
        const Real w = 1.0 - static_cast<Real>(j) / static_cast<Real>(L + 1);
        lam += 2.0 * w * gj;
    }
    return lam;
}

}  // namespace unit_root
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
