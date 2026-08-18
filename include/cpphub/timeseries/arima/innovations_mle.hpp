// =============================================================================
// innovations_mle.hpp - ARMA 精确 MLE (innovations 算法, B&D 2016 §5.2)
//
// Phase 7C v1.7 M1 (PHASE7C_SPEC.md v1.2 §3.1 Step5 / AR8)
//
// 教材锚点: Brockwell-Davis 2016 ITSF 2nd ed §5.2 (innovations 算法)
// 对照库: statsmodels 0.14.6 tsa/arima/estimators/innovations.py +
//   tsa/_innovations.pyx (一手源码 2026-08-18 GitHub main 实录)
//
// 一手源码语义落档 (与 statsmodels 逐位对齐的依据):
//   innovations_algo(acov, nobs) [_innovations.pyx L70-86 实录]:
//     max_lag = max{h : acov[h] ≠ 0}
//     v[0] = acov[0]
//     for i = 1..n−1:
//       for k = max(i−max_lag, 0)..i−1:
//         sub = Σ_{j=max(i−max_lag,0)}^{k−1} θ[k,k−j]·θ[i,i−j]·v[j]
//         θ[i,i−k] = (acov[i−k] − sub) / v[k]
//       v[i] = acov[0] − Σ_{j=max(i−max_lag,0)}^{i−1} θ[i,i−j]²·v[j]
//     返回 θ[:n, 1:] (对角偏移坐标 → 行 i 列 j = 内部 θ[i,j+1])
//   innovations_filter(z, θ) [同文件实录]:
//     u[0] = z[0]; u[i] = z[i] − Σ_{j<min(i,k)} θ[i][j]·u[i−1−j]
//   arma_acovf (arima_process.py L128-201): B&D 2009 eq 3.3.8/3.3.9
//     m = max(p,q)+1; ψ = arma2ma 截断 m 项;
//     A[k, :k+1] = ar[:k+1][::-1], A[k, 1:m−k] += ar[k+1:m]
//     b[k] = σ²·Σ_j ma[k+j]·ψ[j] (j = 0..q−k)
//     γ[:m] = A⁻¹b; γ_h = Σφᵢγ_{h−i} (h ≥ m, B&D 3.3.9)
//   黄金锚 (test_stattools.py L1315-1321, B&D MA(1) θ=−0.9):
//     acov=[1.81,−0.9] → θ 行 [0]/[−0.4972]/[−0.6606]/[−0.7404],
//     v = [1.81, 1.3625, 1.2155, 1.1436] — 本实现 4 位逐位复现 (测试断言)
//
// 幻觉点防护 (spec §9.1):
//   AR8: Innovations 仅限无缺失 + 无季节 + q ≥ 0; d > 0 时差分后估计
//        (statsmodels innovations_mle L155-161 同语义); exog/trend 不支持
//   AR1: MA 正号参数化 (1 + θ₁B + …) — statsmodels ma_params 同号直传
//
// σ² 集中化 (与 statsmodels 差异声明): statsmodels 优化全参数 (φ,θ,σ²)
//   (其 TODO L137-139 明示未集中化); 本实现解析集中化 — θ/u/v 与 σ² 无关
//   (γ/σ² 归一化后 σ² 只进入 log 因子), σ̂² = (Σu²/v)/T 解析最优,
//   同一似然面的最优落点等价, 落点差在优化器容差层 (实测落档)
// =============================================================================

#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/calibration/optimizer.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/arima/hannan_rissanen.hpp"
#include "cpphub/timeseries/unit_root/unit_root_common.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace arima {

/// innovations MLE 结果 (spec §3.1 Step5; ArimaParams 子集 + 回显)
struct InnovationsResult {
    std::vector<Real> phi;    ///< AR (1 − φ₁B − …)
    std::vector<Real> theta;  ///< MA (1 + θ₁B + …) 正号 (AR1)
    Real sigma2 = 0.0;
    Real loglik = 0.0;        ///< 精确高斯 loglik (demean 后 T 个观测)
    std::vector<Real> innovations;  ///< u_t (一步预测误差)
    bool demean_used = true;
    bool converged = false;
    Size n_iterations = 0;
    std::string message;
};

namespace detail {

/// 参数拼接 (φ, θ) → 优化向量 (前置声明供起始点构造)
inline std::vector<Real> concat_params(const std::vector<Real>& a,
                                       const std::vector<Real>& b) {
    std::vector<Real> out(a);
    out.insert(out.end(), b.begin(), b.end());
    return out;
}

/// ar 多项式第 i 个 φ 系数 (ar[i] = −φ_i; 前置声明供 arma_acovf 用)
inline Real phi_coef(const std::vector<Real>& ar, Size i) { return -ar[i]; }

/// AR → MA(∞) 截断系数 (arma2ma 语义): ψ_0 = 1;
/// ψ_j = Σ_{i≤min(j,p)} φ_i·ψ_{j−i} + θ_j (j ≤ q, 否则 θ_j = 0)
/// ar/ma 均为"含零滞后"多项式: ar = (1, −φ₁,…), ma = (1, θ₁,…)
inline std::vector<Real> arma_psi(const std::vector<Real>& ar,
                                  const std::vector<Real>& ma, Size lags) {
    const Size p = ar.size() - 1;
    const Size q = ma.size() - 1;
    std::vector<Real> psi(lags + 1, 0.0);
    psi[0] = 1.0;
    for (Size j = 1; j <= lags; ++j) {
        Real s = 0.0;
        for (Size i = 1; i <= std::min<Size>(j, p); ++i) {
            s -= ar[i] * psi[j - i];  // ar[i] = −φ_i
        }
        if (j <= q) s += ma[j];
        psi[j] = s;
    }
    return psi;
}

/// ARMA 理论自协方差 (statsmodels arma_acovf 逐字语义, B&D 3.3.8/3.3.9)
/// ar = (1, −φ…), ma = (1, θ…); 返回 nobs 项 γ_0..γ_{nobs−1}
/// @throws std::runtime_error 非平稳 (A 奇异) 或 σ² ≤ 0
inline std::vector<Real> arma_acovf(const std::vector<Real>& ar,
                                    const std::vector<Real>& ma,
                                    Real sigma2, Size nobs) {
    const Size p = ar.size() - 1;
    const Size q = ma.size() - 1;
    const Size m = std::max(p, q) + 1;
    if (sigma2 <= 0.0) throw std::runtime_error("arma_acovf: sigma2 <= 0");

    if (p == 0 && q == 0) {
        std::vector<Real> g(nobs, 0.0);
        g[0] = sigma2;
        return g;
    }
    // ψ 权重截断 m 项 (arma2ma lags=m → ψ_0..ψ_m)
    const std::vector<Real> psi = arma_psi(ar, ma, m);

    // A·γ = b (B&D eq. 3.3.8)
    std::vector<std::vector<Real>> A(m, std::vector<Real>(m, 0.0));
    std::vector<Real> b(m, 0.0);
    std::vector<Real> tmp_ar(m, 0.0);
    for (Size i = 0; i <= p && i < m; ++i) tmp_ar[i] = ar[i];
    for (Size k = 0; k < m; ++k) {
        // A[k, :k+1] = tmp_ar[:k+1][::-1]
        for (Size c = 0; c <= k; ++c) {
            A[k][c] = tmp_ar[k - c];
        }
        // A[k, 1:m−k] += tmp_ar[k+1:m]
        for (Size c = k + 1; c < m; ++c) {
            A[k][c] += tmp_ar[c];
        }
        // b[k] = σ²·Σ_{j=0}^{q−k} ma[k+j]·ψ[j]
        Real s = 0.0;
        for (Size j = 0; k + j <= q; ++j) {
            s += ma[k + j] * psi[j];
        }
        b[k] = sigma2 * s;
    }
    // Gauss-Jordan 解 A·x = b (复用 unit_root::detail::invert_matrix;
    // 奇异 throw = 非平稳防护, 由上层 catch 返回惩罚)
    auto Ainv = A;
    try {
        unit_root::detail::invert_matrix(Ainv);
    } catch (const std::runtime_error&) {
        throw std::runtime_error("arma_acovf: non-stationary (singular A)");
    }
    const Size n_out = std::max(nobs, m);
    std::vector<Real> g(n_out, 0.0);
    for (Size r = 0; r < m; ++r) {
        Real s = 0.0;
        for (Size c = 0; c < m; ++c) s += Ainv[r][c] * b[c];
        g[r] = s;
    }
    // B&D eq. 3.3.9: γ_h = Σ_{i=1}^{p} φ_i·γ_{h−i} (h ≥ m)
    for (Size h = m; h < n_out; ++h) {
        Real s = 0.0;
        for (Size i = 1; i <= p; ++i) s += phi_coef(ar, i) * g[h - i];
        g[h] = s;
    }
    g.resize(nobs);
    return g;
}

/// innovations 算法 (statsmodels _innovations.pyx 逐字; γ 可含零尾)
/// 返回 {theta_out (T×max_lag, 行 i 列 j = 内部 θ[i,j+1]), v (T)}
inline void innovations_algo(const std::vector<Real>& acov, Size nobs,
                             std::vector<std::vector<Real>>& theta_out,
                             std::vector<Real>& v) {
    const Size L = acov.size();
    Size max_lag = 0;
    for (Size h = 0; h < L; ++h) {
        if (acov[h] != 0.0) max_lag = h;
    }
    // 内部 θ (nobs+1) × (max_lag+1); 仅用 [i, i−k] 条目 (k ≤ max_lag)
    std::vector<std::vector<Real>> th(nobs, std::vector<Real>(max_lag + 2, 0.0));
    v.assign(nobs, 0.0);
    v[0] = acov[0];
    for (Size i = 1; i < nobs; ++i) {
        const Size k_lo = (i > max_lag) ? (i - max_lag) : 0;
        for (Size k = k_lo; k < i; ++k) {
            Real sub = 0.0;
            for (Size j = k_lo; j < k; ++j) {
                sub += th[k][k - j] * th[i][i - j] * v[j];
            }
            th[i][i - k] = (acov[i - k] - sub) / v[k];
        }
        Real vi = acov[0];
        for (Size j = k_lo; j < i; ++j) {
            vi -= th[i][i - j] * th[i][i - j] * v[j];
        }
        v[i] = vi;
        if (!(vi > 0.0)) {
            throw std::runtime_error(
                "innovations_algo: non-positive prediction variance "
                "(non-stationary/invertibility violation)");
        }
    }
    // 输出裁剪: theta_out[i][j] = th[i][j+1], j = 0..max_lag−1
    theta_out.assign(nobs, std::vector<Real>(max_lag, 0.0));
    for (Size i = 0; i < nobs; ++i) {
        for (Size j = 0; j + 1 <= max_lag; ++j) {
            theta_out[i][j] = th[i][j + 1];
        }
    }
}

/// innovations 滤波 (statsmodels innovations_filter 逐字)
/// u[0] = z[0]; u[i] = z[i] − Σ_{j<min(i,k)} θ[i][j]·u[i−1−j]
inline std::vector<Real> innovations_filter(
        const std::vector<Real>& z,
        const std::vector<std::vector<Real>>& theta) {
    const Size T = z.size();
    const Size k = theta.empty() ? 0 : theta[0].size();
    std::vector<Real> u(T, 0.0);
    u[0] = z[0];
    for (Size i = 1; i < T; ++i) {
        Real hat = 0.0;
        for (Size j = 0; j < std::min<Size>(i, k); ++j) {
            hat += theta[i][j] * u[i - j - 1];
        }
        u[i] = z[i] - hat;
    }
    return u;
}

/// ARMA 精确 loglik 逐观测值 (归一化 γ 口径: σ² 由外部单独进入)
/// γ_norm = arma_acovf(φ,θ,σ²=1); θ/v 由 γ_norm 决定 (σ² 无关);
/// ℓ = −0.5·Σ_t [log(2π·σ²·v_t) + u_t²/(σ²·v_t)]
/// 返回 {S = Σ u²/v (σ² 集中化用), Σ log v}
struct LoglikePieces { Real s_uv = 0.0; Real sum_logv = 0.0; };
inline LoglikePieces arma_innovations_pieces(
        const std::vector<Real>& z, const std::vector<Real>& phi,
        const std::vector<Real>& theta) {
    const Size T = z.size();
    std::vector<Real> ar(1, 1.0);
    for (Real p : phi) ar.push_back(-p);
    std::vector<Real> ma(1, 1.0);
    for (Real t : theta) ma.push_back(t);
    const std::vector<Real> gamma = arma_acovf(ar, ma, 1.0, T);
    std::vector<std::vector<Real>> th;
    std::vector<Real> v;
    innovations_algo(gamma, T, th, v);
    const std::vector<Real> u = innovations_filter(z, th);
    LoglikePieces pc;
    for (Size t = 0; t < T; ++t) {
        pc.s_uv += u[t] * u[t] / v[t];
        pc.sum_logv += std::log(v[t]);
    }
    if (!std::isfinite(pc.s_uv) || !std::isfinite(pc.sum_logv)) {
        throw std::runtime_error("arma loglik: non-finite pieces");
    }
    return pc;
}

/// 集中化负 loglik (σ² 解析最优): 参数 x = {φ₁..φ_p, θ₁..θ_q}
inline Real concentrated_nll(const std::vector<Real>& x,
                             const std::vector<Real>& z, Size p, Size q) {
    try {
        std::vector<Real> phi(x.begin(), x.begin() + p);
        std::vector<Real> theta(x.begin() + p, x.end());
        const auto pc = arma_innovations_pieces(z, phi, theta);
        const Real T = static_cast<Real>(z.size());
        const Real s2 = pc.s_uv / T;
        if (!(s2 > 0.0)) return 1e100;
        return 0.5 * T * std::log(s2) + 0.5 * pc.sum_logv;
    } catch (...) {
        return 1e100;  // 非平稳/发散参数惩罚
    }
}

}  // namespace detail

/// @brief ARMA(p,q) 精确 MLE — innovations 算法 (B&D 2016 §5.2)
///
/// 数据预处理 (statsmodels innovations_mle L155-163 语义):
///   d > 0 → 差分 d 次; demean=true → 减样本均值 (AR8: d>0 时 statsmodels
///   对差分序列 demean; 均值不可识别时设 false)
/// 起始: HR (hannan_rissanen) + 零向量 + 随机扰动 (Philox, seed 固定);
///   非平稳/不可逆起始按 statsmodels L198-204 置零处理
/// 优化: SLSQP 无约束 (bounds ±0.999 防护) on 集中化 nll; σ̂² 解析
///
/// @param z 序列 (差分由本函数按 d 执行)
/// @param p, q AR/MA 阶
/// @param d 差分阶 (≥ 0)
/// @param demean 差分后是否减均值 (默认 true = statsmodels)
/// @param seed 随机扰动起始 (§1.4-8)
/// @throws std::invalid_argument T 过小 / NaN / AR8 违例 (NaN 输入)
inline InnovationsResult innovations_mle(const std::vector<Real>& z_input,
                                         Size p, Size q, Size d = 0,
                                         bool demean = true, Size seed = 42) {
    // 输入校验 (AR8: 无缺失)
    for (Real v : z_input) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument(
                "innovations_mle: NaN/Inf not supported (AR8)");
        }
    }
    if (p + q == 0) {
        throw std::invalid_argument("innovations_mle: p=q=0 trivial");
    }
    // 差分 d
    std::vector<Real> z = z_input;
    for (Size i = 0; i < d; ++i) {
        std::vector<Real> dz(z.size() - 1);
        for (Size t = 1; t < z.size(); ++t) dz[t - 1] = z[t] - z[t - 1];
        z = dz;
    }
    Real mean = 0.0;
    if (demean) {
        for (Real v : z) mean += v;
        mean /= static_cast<Real>(z.size());
        for (Real& v : z) v -= mean;
    }
    const Size T = z.size();
    if (T < p + q + 5) {
        throw std::invalid_argument("innovations_mle: sample too small");
    }

    const Size k = p + q;
    // 起始点集: {HR, 零, 零+扰动×2} (AR7 精神; statsmodels 仅 HR)
    std::vector<std::vector<Real>> starts;
    starts.push_back(std::vector<Real>(k, 0.0));
    try {
        const auto hr = hannan_rissanen(z, p, q);
        starts.push_back(detail::concat_params(hr.phi, hr.theta));
    } catch (...) {
        // HR 失败 (退化) → 仅零起始
    }
    {
        Philox4x64 rng(seed, 11);
        for (Size rep = 0; rep < 2; ++rep) {
            std::vector<Real> x0(k, 0.0);
            for (Size i = 0; i < k; ++i) {
                const uint64_t r1 = rng();
                const uint64_t r2 = rng();
                const Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
                const Real u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
                const auto [a, b] = box_muller(u1, u2);
                x0[i] = 0.2 * (rep == 0 ? a : b);
            }
            starts.push_back(x0);
        }
    }

    SLSQP::Config cfg;
    cfg.max_iterations = 300;
    cfg.ftol = 1e-11;
    cfg.xtol = 1e-10;

    std::vector<Bounds> bounds(k, Bounds{-0.999, 0.999});
    Real best_f = std::numeric_limits<Real>::infinity();
    std::vector<Real> best_x;
    bool conv = false;
    Size niter = 0;
    std::string msg;
    for (const auto& x0 : starts) {
        OptimizationResult r = SLSQP::minimize(
            [&](const std::vector<Real>& x) {
                return detail::concentrated_nll(x, z, p, q);
            },
            x0, bounds, {}, {}, cfg);
        if (r.fx < best_f) {
            best_f = r.fx;
            best_x = r.x;
            conv = r.converged;
            niter = r.n_iterations;
            msg = r.message;
        }
    }

    InnovationsResult res;
    res.phi.assign(best_x.begin(), best_x.begin() + p);
    res.theta.assign(best_x.begin() + p, best_x.end());
    res.demean_used = demean;
    res.converged = conv;
    res.n_iterations = niter;
    res.message = msg;
    // σ̂² 解析 + 完整 loglik (含 demean 后 T 观测)
    const auto pc = detail::arma_innovations_pieces(z, res.phi, res.theta);
    const Real Tn = static_cast<Real>(T);
    res.sigma2 = pc.s_uv / Tn;
    constexpr Real kTwoPi = 6.283185307179586476925286766559;
    res.loglik = -0.5 * (Tn * std::log(kTwoPi * res.sigma2)
                         + Tn + pc.sum_logv);
    // innovations (一步预测误差) 输出
    {
        std::vector<Real> ar(1, 1.0);
        for (Real p : res.phi) ar.push_back(-p);
        std::vector<Real> ma(1, 1.0);
        for (Real t : res.theta) ma.push_back(t);
        const auto gamma = detail::arma_acovf(ar, ma, 1.0, T);
        std::vector<std::vector<Real>> th;
        std::vector<Real> v;
        detail::innovations_algo(gamma, T, th, v);
        res.innovations = detail::innovations_filter(z, th);
    }
    return res;
}

}  // namespace arima
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
