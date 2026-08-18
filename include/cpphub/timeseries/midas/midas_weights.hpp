// =============================================================================
// midas_weights.hpp - MIDAS 权重族五件 (spec §6.2 / 决策 22)
//
// Phase 7C v1.7 M4 (PHASE7C_SPEC.md v1.2 §6.2)
//
// 对照库: R midasr 0.9 (唯一主数值基准, B5) — 源码一手实录
//   2026-08-18 print(midasr::nealmon/nbetaMT/almonp/polystep/harstep):
//
//   nealmon(p, d, m):
//     i <- 1:d; plc <- poly(i, degree=len(p)-1, raw=TRUE) %*% p[-1]
//     p[1] * exp(plc)/sum(exp(plc))
//     → w_i = δ·exp(Σ_k λ_{k+1}·i^k)/Σexp(·), i = 1..d (MD1 i 从 1 起);
//       返回 Σw = δ = p[1] (MD2 δ 独立线性参数, 非归一化 1)
//     ⚠️ 裸 exp 在 plc 大值溢出 — midasr 无防护; 本实现 log-sum-exp
//       (决策 25/MD7), 非溢出区间与裸公式差 < 1e-14 (测试断言)
//
//   nbetaMT(p, d, m)  [4 参数 (δ, κ₁, κ₂, θ₀) — spec §6.2 仅记 κ₁κ₂,
//     实施期勘误: 第 4 参数 = 零假设混合权重 θ₀, Ghysels 2016 SEM']:
//     xi <- (1:d−1)/(d−1); xi[1]+=eps; xi[d]−=eps   (xi 从 0 起, 端点防护)
//     nb <- xi^(κ₁−1)·(1−xi)^(κ₂−1)
//     sum(nb)<eps: |θ₀|<eps → 零向量; 否则 δ·rep(1/d)
//     w <- nb/sum(nb) + θ₀; δ·w/sum(w)                (Σw = δ)
//     ⚠️ d = 1 → (0)/(0) = NaN (midasr 同), 本实现拒绝 d < 2
//
//   almonp(p, d, m):
//     w_i = p₁ + Σ_k p_{k+1}·i^k (raw poly, i = 1..d) — 不归一化, 可负
//
//   polystep(p, d, m, a): 阶梯 — 断点 a ⊂ (1,d), 长度 len(p)−1;
//     w = rep(p, times = diff(c(0, a, d)))
//     ⚠️ spec §6.2 签名 (a,b,d) 二步特例 → 通用 vector 断点 (实施勘误,
//     与 midasr 语义一致); a 越界 (≤0 或 ≥d) 拒绝
//
//   harstep(p, d, m): HAR(3)-RV — d = 20 硬编码 (midasr 同):
//     w₁ = p₁+p₂/5+p₃/20; w₂..₅ = p₂/5+p₃/20; w₆..₂₀ = p₃/20
//
// 参数化总则 (MD2): 所有权重函数含 δ 总尺度于首参数; 回归中 δ 消去到
//   内层线性参数 (集中化 NLS 的内层 OLS 吸收), midas_weights 本身按
//   midasr 逐点复刻 (含 δ) 以支撑 1e-12 对照
// =============================================================================

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace midas {

/// 权重族选择 (spec §6.2) — 语义见各 weights 函数
enum class MidasWeight { Nealmon, NBeta, AlmonP, PolyStep, HarStep };

/// log-sum-exp (决策 25 / MD7): logΣexp(η) = max(η) + log Σ exp(η − max η)
/// 防溢出: 裸 exp 在 λ₂·d² ≳ 709 溢出; 非溢出区间与裸公式差 < 1e-14
inline Real log_sum_exp(const std::vector<Real>& eta) {
    if (eta.empty()) throw std::invalid_argument("log_sum_exp: empty");
    Real mx = eta[0];
    for (Real e : eta) mx = std::max(mx, e);
    if (!std::isfinite(mx)) return mx;  // −inf 输入 → −inf
    Real s = 0.0;
    for (Real e : eta) s += std::exp(e - mx);
    return mx + std::log(s);
}

/// @brief nealmon 权重 (指数 Almon; MD1: i = 1..d 从 1 起)
/// @param lambda {δ, λ₂, λ₃, …} — λ₁ = δ 总尺度 (Σw = δ, MD2)
/// @param d 高频滞后数 (≥ 1)
/// @return w_i = δ·exp(Σ_k λ_{k+1}·i^k)/Σexp(·), i = 1..d (log-sum-exp 实现)
inline std::vector<Real> nealmon_weights(const std::vector<Real>& lambda,
                                         Size d) {
    if (lambda.size() < 2) {
        // 仅 δ: midasr poly degree=0 → plc = 0 → w = δ/d
        if (lambda.size() == 1) {
            return std::vector<Real>(d, lambda[0] / static_cast<Real>(d));
        }
        throw std::invalid_argument("nealmon: lambda = {δ, ...} ≥ 1");
    }
    if (d == 0) throw std::invalid_argument("nealmon: d ≥ 1");
    std::vector<Real> eta(d);
    for (Size i = 1; i <= d; ++i) {
        Real e = 0.0;
        Real ipow = 1.0;
        for (Size k = 1; k < lambda.size(); ++k) {
            ipow *= static_cast<Real>(i);  // i^k
            e += lambda[k] * ipow;
        }
        eta[i - 1] = e;
    }
    const Real lse = log_sum_exp(eta);
    std::vector<Real> w(d);
    for (Size i = 0; i < d; ++i) {
        w[i] = lambda[0] * std::exp(eta[i] - lse);
    }
    return w;
}

/// @brief nbetaMT 权重 (β 密度 + 零假设混合; 实施勘误: 4 参数)
/// @param lambda {δ, κ₁, κ₂, θ₀} — θ₀ = 均匀混合权重 (midasr p[4])
/// @param d ≥ 2 (d=1 → midasr NaN, 拒绝)
/// @return Σw = δ; xi = (i−1)/(d−1) 从 0 起 (与 nealmon 不同!)
inline std::vector<Real> nbeta_weights(const std::vector<Real>& lambda,
                                        Size d) {
    if (lambda.size() != 4) {
        throw std::invalid_argument("nbeta: lambda = {δ,κ₁,κ₂,θ₀} (4)");
    }
    if (d < 2) {
        throw std::invalid_argument("nbeta: d ≥ 2 (d=1 → NaN in midasr)");
    }
    constexpr Real kEps = std::numeric_limits<double>::epsilon();
    const Real d_r = static_cast<Real>(d);
    std::vector<Real> xi(d);
    for (Size i = 1; i <= d; ++i) {
        xi[i - 1] = static_cast<Real>(i - 1) / static_cast<Real>(d - 1);
    }
    xi[0] += kEps;
    xi[d - 1] -= kEps;
    std::vector<Real> nb(d);
    Real sum_nb = 0.0;
    for (Size i = 0; i < d; ++i) {
        // log 域防溢出: xi^(κ₁−1)·(1−xi)^(κ₂−1)
        const Real lg = (lambda[1] - 1.0) * std::log(xi[i])
                        + (lambda[2] - 1.0) * std::log(1.0 - xi[i]);
        nb[i] = std::exp(lg);
        sum_nb += nb[i];
    }
    if (sum_nb < kEps) {
        if (std::fabs(lambda[3]) < kEps) {
            return std::vector<Real>(d, 0.0);
        }
        return std::vector<Real>(d, lambda[0] / d_r);
    }
    std::vector<Real> w(d);
    Real sum_w = 0.0;
    for (Size i = 0; i < d; ++i) {
        w[i] = nb[i] / sum_nb + lambda[3];
        sum_w += w[i];
    }
    for (Size i = 0; i < d; ++i) {
        w[i] = lambda[0] * w[i] / sum_w;
    }
    return w;
}

/// @brief almonp 权重 (raw 多项式; 不归一化可负)
/// @param lambda {c, a₁, a₂, …} — w_i = c + Σ_k a_k·i^k, i = 1..d
inline std::vector<Real> almonp_weights(const std::vector<Real>& lambda,
                                        Size d) {
    if (lambda.empty()) throw std::invalid_argument("almonp: empty");
    if (d == 0) throw std::invalid_argument("almonp: d ≥ 1");
    std::vector<Real> w(d);
    for (Size i = 1; i <= d; ++i) {
        Real v = lambda[0];
        Real ipow = 1.0;
        for (Size k = 1; k < lambda.size(); ++k) {
            ipow *= static_cast<Real>(i);
            v += lambda[k] * ipow;
        }
        w[i - 1] = v;
    }
    return w;
}

/// @brief polystep 阶梯权重 (通用断点, 实施勘误: midasr 语义非二步特例)
/// @param values 每段常值 {v₁, …, v_{s+1}} (s = steps.size())
/// @param steps 断点 a ⊂ (0, d) 整数, 严格递增; 段 r 覆盖 i ∈ (a_{r−1}, a_r]
/// @param d 总滞后数; w = rep(values, times = diff(c(0, a, d)))
inline std::vector<Real> polystep_weights(const std::vector<Real>& values,
                                          const std::vector<Size>& steps,
                                          Size d) {
    if (steps.size() + 1 != values.size()) {
        throw std::invalid_argument(
            "polystep: len(values) = len(steps) + 1");
    }
    Size prev = 0;
    std::vector<Real> w;
    w.reserve(d);
    for (Size s = 0; s < steps.size(); ++s) {
        if (steps[s] <= prev || steps[s] >= d) {
            throw std::invalid_argument(
                "polystep: steps strictly inside (0, d)");
        }
        w.insert(w.end(), steps[s] - prev, values[s]);
        prev = steps[s];
    }
    w.insert(w.end(), d - prev, values.back());
    return w;
}

/// @brief harstep 权重 (HAR(3)-RV; midasr 硬编码 d = 20)
/// @param lambda {p₁, p₂, p₃} (日/周/月系数)
/// @param d 必须为 20 (midasr 同, 其他 d 拒绝)
inline std::vector<Real> harstep_weights(const std::vector<Real>& lambda,
                                         Size d) {
    if (d != 20) {
        throw std::invalid_argument("harstep: d = 20 (HAR(3)-RV, midasr)");
    }
    if (lambda.size() != 3) {
        throw std::invalid_argument("harstep: lambda = {p1,p2,p3}");
    }
    std::vector<Real> w(20, 0.0);
    w[0] = lambda[0] + lambda[1] / 5.0 + lambda[2] / 20.0;
    for (Size i = 1; i <= 4; ++i) {
        w[i] = lambda[1] / 5.0 + lambda[2] / 20.0;
    }
    for (Size i = 5; i < 20; ++i) {
        w[i] = lambda[2] / 20.0;
    }
    return w;
}

}  // namespace midas
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
