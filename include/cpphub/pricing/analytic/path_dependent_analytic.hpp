#pragma once
// SOURCE: Kemna & Vorst (1990) "A pricing method for options based on average asset values"
//         J. Banking & Finance 14, 113-129.
// SOURCE: Reiner & Rubinstein (1991) "Unscrambling the Binary Code" Risk 4, 75-83.
// SOURCE: Goldman, Sosin, Gatto (1979) "Path dependent options: buy at the low, sell at the high"
//         J. Finance 34, 1111-1127.
// SOURCE: Haug (2007) "The Complete Guide to Option Pricing Formulas" Ch.4 (Barrier/Lookback)
// SOURCE: Hull (2018) "Options, Futures, and Other Derivatives" Ch.26 (Path-Dependent)
// SOURCE: Conze & Viswanathan (1991) "Path dependent options: The case of lookback options"
//         J. Finance 46, 1893-1907.
//
// 模块: 路径相关衍生品解析定价 (GBM 假设)
//
// 数学概要:
//   1. Geometric Asian (Kemna-Vorst):
//      GBM 下几何平均 G = (∏ S(tᵢ))^{1/n} 服从对数正态分布.
//      连续: ln G ~ N(ln S₀ + (b-σ²/2)T/2, σ²T/3)
//      离散: ln G ~ N(ln S₀ + (b-σ²/2)T(n+1)/(2n), σ²T(n+1)(2n+1)/(6n²))
//      其中 b = r - q. BSM-like 公式用调整后的漂移和波动率.
//
//   2. Barrier (Reiner-Rubinstein):
//      反射原理: 若 S(t) 触及 H, 则 S'(t) = H²/S(t) 也是合法 GBM (测度变换).
//      反射因子: (H/S)^{2μ}, 其中 μ = (b - σ²/2)/σ², b = r - q.
//      反射现货: S' = H²/S.
//      障碍期权 = BSM - 反射项, 障碍 Parity: In = Vanilla - Out.
//
//   3. Lookback (Goldman-Sosin-Gatto):
//      浮动执行价 Call (payoff = S_T - m_T): 利用 min 的反射原理.
//      浮动执行价 Put  (payoff = M_T - S_T): 利用 max 的反射原理.
//      公式含 σ²/(2b) 项, 反映极值过程的分布.

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ 内部辅助函数 ============
namespace detail {

// BSM d1: (ln(S/K) + (b + σ²/2)T) / (σ√T), b = r - q (cost of carry)
inline Real bsm_d1(Real S, Real K, Real T, Real b, Real sigma) {
    return (std::log(S / K) + (b + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
}

inline Real bsm_d2(Real S, Real K, Real T, Real b, Real sigma) {
    return bsm_d1(S, K, T, b, sigma) - sigma * std::sqrt(T);
}

// BSM call (discounted): S*e^{-qT}*N(d1) - K*e^{-rT}*N(d2)
inline Real bsm_call(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(S - K, 0.0);
    if (sigma <= 0.0) return S * std::exp(-q * T) - K * std::exp(-r * T);
    Real b = r - q;
    Real d1 = bsm_d1(S, K, T, b, sigma);
    Real d2 = d1 - sigma * std::sqrt(T);
    return S * std::exp(-q * T) * normal_cdf(d1) - K * std::exp(-r * T) * normal_cdf(d2);
}

// BSM put (discounted): K*e^{-rT}*N(-d2) - S*e^{-qT}*N(-d1)
inline Real bsm_put(Real S, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(K - S, 0.0);
    if (sigma <= 0.0) return K * std::exp(-r * T) - S * std::exp(-q * T);
    Real b = r - q;
    Real d1 = bsm_d1(S, K, T, b, sigma);
    Real d2 = d1 - sigma * std::sqrt(T);
    return K * std::exp(-r * T) * normal_cdf(-d2) - S * std::exp(-q * T) * normal_cdf(-d1);
}

// 反射因子 μ = (b - σ²/2)/σ²
inline Real reflection_mu(Real b, Real sigma) {
    return (b - 0.5 * sigma * sigma) / (sigma * sigma);
}

// 反射现货 S' = H²/S
inline Real reflected_spot(Real S, Real H) {
    return H * H / S;
}

// ============ Down-barrier building blocks (H < S) ============
// P_down(X) = P(S_T > X, min > H) = N(d2(S,X)) - (H/S)^{2μ} * N(d2(S',X))
inline Real P_down(Real S, Real X, Real H, Real T, Real b, Real sigma) {
    Real mu = reflection_mu(b, sigma);
    Real Sp = reflected_spot(S, H);
    Real factor = std::pow(H / S, 2.0 * mu);
    return normal_cdf(bsm_d2(S, X, T, b, sigma)) - factor * normal_cdf(bsm_d2(Sp, X, T, b, sigma));
}

// A_down(X) = E[S_T * I(S_T > X, min > H)] = S*e^{bT} * [N(d1(S,X)) - (H/S)^{2μ+2} * N(d1(S',X))]
inline Real A_down(Real S, Real X, Real H, Real T, Real b, Real sigma) {
    Real mu = reflection_mu(b, sigma);
    Real Sp = reflected_spot(S, H);
    Real factor2 = std::pow(H / S, 2.0 * mu + 2.0);
    return S * std::exp(b * T) * (normal_cdf(bsm_d1(S, X, T, b, sigma))
                                  - factor2 * normal_cdf(bsm_d1(Sp, X, T, b, sigma)));
}

// ============ Up-barrier building blocks (H > S) ============
// P_up(X) = P(S_T < X, max < H) = N(-d2(S,X)) - (H/S)^{2μ} * N(-d2(S',X))   for X ≤ H
inline Real P_up(Real S, Real X, Real H, Real T, Real b, Real sigma) {
    Real mu = reflection_mu(b, sigma);
    Real Sp = reflected_spot(S, H);
    Real factor = std::pow(H / S, 2.0 * mu);
    return normal_cdf(-bsm_d2(S, X, T, b, sigma)) - factor * normal_cdf(-bsm_d2(Sp, X, T, b, sigma));
}

// A_up(X) = E[S_T * I(S_T < X, max < H)] = S*e^{bT} * [N(-d1(S,X)) - (H/S)^{2μ+2} * N(-d1(S',X))]
inline Real A_up(Real S, Real X, Real H, Real T, Real b, Real sigma) {
    Real mu = reflection_mu(b, sigma);
    Real Sp = reflected_spot(S, H);
    Real factor2 = std::pow(H / S, 2.0 * mu + 2.0);
    return S * std::exp(b * T) * (normal_cdf(-bsm_d1(S, X, T, b, sigma))
                                  - factor2 * normal_cdf(-bsm_d1(Sp, X, T, b, sigma)));
}

}  // namespace detail

// ================================================================
// ============ 1. Geometric Asian Options (Kemna-Vorst) ============
// ================================================================

// 离散几何平均亚式期权 (n_steps 个等距观察点 t_i = i*dt, dt = T/n, i=1..n)
// ln G = (1/n) Σ ln S(t_i), G 服从对数正态分布
// μ_G = ln S₀ + (b-σ²/2) * T(n+1)/(2n)
// σ_G² = σ² * T(n+1)(2n+1)/(6n²)
//
// Call = e^{-rT} * [E[G]*N(d1) - K*N(d2)]
// Put  = e^{-rT} * [K*N(-d2) - E[G]*N(-d1)]
// d1 = (μ_G - ln K + σ_G²) / σ_G,  d2 = d1 - σ_G
// E[G] = exp(μ_G + σ_G²/2)

inline Real geom_asian_call_price(Real S0, Real K, Real T, Real r, Real q, Real sigma,
                                   Size n_steps) {
    if (T <= 0.0) return std::max(S0 - K, 0.0);
    if (n_steps == 0) throw std::invalid_argument("geom_asian: n_steps must be positive");
    if (sigma <= 0.0) {
        Real b = r - q;
        Real n = static_cast<Real>(n_steps);
        Real mu_g = std::log(S0) + (b - 0.5 * sigma * sigma) * T * (n + 1.0) / (2.0 * n);
        Real G = std::exp(mu_g);
        return std::exp(-r * T) * std::max(G - K, 0.0);
    }
    Real b = r - q;
    Real n = static_cast<Real>(n_steps);
    Real mu_g = std::log(S0) + (b - 0.5 * sigma * sigma) * T * (n + 1.0) / (2.0 * n);
    Real var_g = sigma * sigma * T * (n + 1.0) * (2.0 * n + 1.0) / (6.0 * n * n);
    Real sg = std::sqrt(var_g);
    Real EG = std::exp(mu_g + 0.5 * var_g);
    Real d1 = (mu_g - std::log(K) + var_g) / sg;
    Real d2 = d1 - sg;
    return std::exp(-r * T) * (EG * normal_cdf(d1) - K * normal_cdf(d2));
}

inline Real geom_asian_put_price(Real S0, Real K, Real T, Real r, Real q, Real sigma,
                                  Size n_steps) {
    if (T <= 0.0) return std::max(K - S0, 0.0);
    if (n_steps == 0) throw std::invalid_argument("geom_asian: n_steps must be positive");
    if (sigma <= 0.0) {
        Real b = r - q;
        Real n = static_cast<Real>(n_steps);
        Real mu_g = std::log(S0) + (b - 0.5 * sigma * sigma) * T * (n + 1.0) / (2.0 * n);
        Real G = std::exp(mu_g);
        return std::exp(-r * T) * std::max(K - G, 0.0);
    }
    Real b = r - q;
    Real n = static_cast<Real>(n_steps);
    Real mu_g = std::log(S0) + (b - 0.5 * sigma * sigma) * T * (n + 1.0) / (2.0 * n);
    Real var_g = sigma * sigma * T * (n + 1.0) * (2.0 * n + 1.0) / (6.0 * n * n);
    Real sg = std::sqrt(var_g);
    Real EG = std::exp(mu_g + 0.5 * var_g);
    Real d1 = (mu_g - std::log(K) + var_g) / sg;
    Real d2 = d1 - sg;
    return std::exp(-r * T) * (K * normal_cdf(-d2) - EG * normal_cdf(-d1));
}

// 连续几何平均亚式期权 (积分形式, G = exp((1/T) ∫₀ᵀ ln S(t) dt))
// μ_G = ln S₀ + (b-σ²/2) * T/2
// σ_G² = σ² * T/3
inline Real geom_asian_call_continuous(Real S0, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(S0 - K, 0.0);
    if (sigma <= 0.0) {
        Real b = r - q;
        Real mu_g = std::log(S0) + (b - 0.5 * sigma * sigma) * T * 0.5;
        return std::exp(-r * T) * std::max(std::exp(mu_g) - K, 0.0);
    }
    Real b = r - q;
    Real mu_g = std::log(S0) + (b - 0.5 * sigma * sigma) * T * 0.5;
    Real var_g = sigma * sigma * T / 3.0;
    Real sg = std::sqrt(var_g);
    Real EG = std::exp(mu_g + 0.5 * var_g);
    Real d1 = (mu_g - std::log(K) + var_g) / sg;
    Real d2 = d1 - sg;
    return std::exp(-r * T) * (EG * normal_cdf(d1) - K * normal_cdf(d2));
}

inline Real geom_asian_put_continuous(Real S0, Real K, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(K - S0, 0.0);
    if (sigma <= 0.0) {
        Real b = r - q;
        Real mu_g = std::log(S0) + (b - 0.5 * sigma * sigma) * T * 0.5;
        return std::exp(-r * T) * std::max(K - std::exp(mu_g), 0.0);
    }
    Real b = r - q;
    Real mu_g = std::log(S0) + (b - 0.5 * sigma * sigma) * T * 0.5;
    Real var_g = sigma * sigma * T / 3.0;
    Real sg = std::sqrt(var_g);
    Real EG = std::exp(mu_g + 0.5 * var_g);
    Real d1 = (mu_g - std::log(K) + var_g) / sg;
    Real d2 = d1 - sg;
    return std::exp(-r * T) * (K * normal_cdf(-d2) - EG * normal_cdf(-d1));
}

// ================================================================
// ============ 2. Barrier Options (Reiner-Rubinstein) ============
// ================================================================

enum class BarrierType {
    DownInCall, DownOutCall, DownInPut, DownOutPut,
    UpInCall,   UpOutCall,   UpInPut,   UpOutPut
};

// 障碍期权统一定价接口
// S: 现货, K: 行权价, H: 障碍水平, T: 到期, r: 无风险利率, q: 股息率, sigma: 波动率
// 要求: H < S (down barrier) 或 H > S (up barrier)
//
// 公式推导基于反射原理 building blocks:
//   Down (H < S): P_down(X) = P(S_T > X, min > H),  A_down(X) = E[S_T * I(S_T > X, min > H)]
//   Up   (H > S): P_up(X)   = P(S_T < X, max < H),  A_up(X)   = E[S_T * I(S_T < X, max < H)]
//
// Out 期权 (生存至到期):
//   DOC (H≤K) = e^{-rT}[A_down(K) - K*P_down(K)]
//   DOC (H>K) = e^{-rT}[A_down(H) - K*P_down(H)]   (S_T > H > K, 无期权性)
//   DOP (H<K) = e^{-rT}[K*(P_down(H)-P_down(K)) - (A_down(H)-A_down(K))]
//   DOP (H≥K) = 0                                    (min > H ≥ K ⟹ S_T > K)
//   UOC (K<H) = e^{-rT}[(A_up(H)-A_up(K)) - K*(P_up(H)-P_up(K))]
//   UOC (K≥H) = 0                                    (max < H ≤ K ⟹ S_T < K)
//   UOP (K<H) = e^{-rT}[K*P_up(K) - A_up(K)]
//   UOP (K≥H) = e^{-rT}[K*P_up(H) - A_up(H)]        (S_T < H ≤ K, 无期权性)
// In = Vanilla - Out (parity)

inline Real barrier_option_price(Real S, Real K, Real H, Real T, Real r, Real q,
                                  Real sigma, BarrierType type) {
    if (T <= 0.0) {
        // 到期: 检查是否触及障碍
        // 简化: 假设未触及 (连续监控下概率为 0)
        switch (type) {
            case BarrierType::DownOutCall: return std::max(S - K, 0.0);
            case BarrierType::DownOutPut:  return std::max(K - S, 0.0);
            case BarrierType::UpOutCall:   return std::max(S - K, 0.0);
            case BarrierType::UpOutPut:    return std::max(K - S, 0.0);
            case BarrierType::DownInCall:  return 0.0;
            case BarrierType::DownInPut:   return 0.0;
            case BarrierType::UpInCall:    return 0.0;
            case BarrierType::UpInPut:     return 0.0;
        }
    }
    if (sigma <= 0.0) {
        // 零波动率: 路径确定, 是否触及障碍取决于漂移方向
        // 简化处理: 用 BSM 公式
        Real vc = detail::bsm_call(S, K, T, r, q, sigma);
        Real vp = detail::bsm_put(S, K, T, r, q, sigma);
        bool is_call = (type == BarrierType::DownInCall || type == BarrierType::DownOutCall ||
                        type == BarrierType::UpInCall || type == BarrierType::UpOutCall);
        return is_call ? vc : vp;
    }
    if (S <= 0.0 || K < 0.0 || H <= 0.0) throw std::invalid_argument("barrier: S, H must be positive, K non-negative");

    Real b = r - q;
    Real df = std::exp(-r * T);
    Real vanilla = 0.0;

    // Vanilla for parity (compute once, needed for In options)
    bool need_vanilla = (type == BarrierType::DownInCall || type == BarrierType::DownInPut ||
                         type == BarrierType::UpInCall || type == BarrierType::UpInPut);
    if (need_vanilla) {
        bool is_call = (type == BarrierType::DownInCall || type == BarrierType::UpInCall);
        vanilla = is_call ? detail::bsm_call(S, K, T, r, q, sigma)
                          : detail::bsm_put(S, K, T, r, q, sigma);
    }

    switch (type) {
        // ---- Down-and-Out Call (H < S) ----
        case BarrierType::DownOutCall: {
            if (H >= S) throw std::invalid_argument("DownOutCall: H must be < S");
            Real out;
            if (H <= K) {
                out = df * (detail::A_down(S, K, H, T, b, sigma) - K * detail::P_down(S, K, H, T, b, sigma));
            } else {
                out = df * (detail::A_down(S, H, H, T, b, sigma) - K * detail::P_down(S, H, H, T, b, sigma));
            }
            return out;
        }
        // ---- Down-and-In Call ----
        case BarrierType::DownInCall: {
            if (H >= S) throw std::invalid_argument("DownInCall: H must be < S");
            Real out;
            if (H <= K) {
                out = df * (detail::A_down(S, K, H, T, b, sigma) - K * detail::P_down(S, K, H, T, b, sigma));
            } else {
                out = df * (detail::A_down(S, H, H, T, b, sigma) - K * detail::P_down(S, H, H, T, b, sigma));
            }
            return vanilla - out;
        }
        // ---- Down-and-Out Put (H < S) ----
        case BarrierType::DownOutPut: {
            if (H >= S) throw std::invalid_argument("DownOutPut: H must be < S");
            Real out;
            if (H < K) {
                // (K-S_T)^+ * I(min > H) = (K-S_T) * I(H < S_T < K, min > H)
                out = df * (K * (detail::P_down(S, H, H, T, b, sigma) - detail::P_down(S, K, H, T, b, sigma))
                          - (detail::A_down(S, H, H, T, b, sigma) - detail::A_down(S, K, H, T, b, sigma)));
            } else {
                out = 0.0;  // H >= K: min > H >= K ⟹ S_T > K, put payoff = 0
            }
            return out;
        }
        // ---- Down-and-In Put ----
        case BarrierType::DownInPut: {
            if (H >= S) throw std::invalid_argument("DownInPut: H must be < S");
            Real out;
            if (H < K) {
                out = df * (K * (detail::P_down(S, H, H, T, b, sigma) - detail::P_down(S, K, H, T, b, sigma))
                          - (detail::A_down(S, H, H, T, b, sigma) - detail::A_down(S, K, H, T, b, sigma)));
            } else {
                out = 0.0;
            }
            return vanilla - out;
        }
        // ---- Up-and-Out Call (H > S) ----
        case BarrierType::UpOutCall: {
            if (H <= S) throw std::invalid_argument("UpOutCall: H must be > S");
            Real out;
            if (K < H) {
                // (S_T-K)^+ * I(max < H) = (S_T-K) * I(K < S_T < H, max < H)
                out = df * ((detail::A_up(S, H, H, T, b, sigma) - detail::A_up(S, K, H, T, b, sigma))
                          - K * (detail::P_up(S, H, H, T, b, sigma) - detail::P_up(S, K, H, T, b, sigma)));
            } else {
                out = 0.0;  // K >= H: max < H <= K ⟹ S_T < K, call payoff = 0
            }
            return out;
        }
        // ---- Up-and-In Call ----
        case BarrierType::UpInCall: {
            if (H <= S) throw std::invalid_argument("UpInCall: H must be > S");
            Real out;
            if (K < H) {
                out = df * ((detail::A_up(S, H, H, T, b, sigma) - detail::A_up(S, K, H, T, b, sigma))
                          - K * (detail::P_up(S, H, H, T, b, sigma) - detail::P_up(S, K, H, T, b, sigma)));
            } else {
                out = 0.0;
            }
            return vanilla - out;
        }
        // ---- Up-and-Out Put (H > S) ----
        case BarrierType::UpOutPut: {
            if (H <= S) throw std::invalid_argument("UpOutPut: H must be > S");
            Real out;
            if (K < H) {
                // (K-S_T)^+ * I(max < H) = (K-S_T) * I(S_T < K, max < H)
                out = df * (K * detail::P_up(S, K, H, T, b, sigma) - detail::A_up(S, K, H, T, b, sigma));
            } else {
                // K >= H: max < H <= K, S_T < H <= K, put always in the money
                out = df * (K * detail::P_up(S, H, H, T, b, sigma) - detail::A_up(S, H, H, T, b, sigma));
            }
            return out;
        }
        // ---- Up-and-In Put ----
        case BarrierType::UpInPut: {
            if (H <= S) throw std::invalid_argument("UpInPut: H must be > S");
            Real out;
            if (K < H) {
                out = df * (K * detail::P_up(S, K, H, T, b, sigma) - detail::A_up(S, K, H, T, b, sigma));
            } else {
                out = df * (K * detail::P_up(S, H, H, T, b, sigma) - detail::A_up(S, H, H, T, b, sigma));
            }
            return vanilla - out;
        }
    }
    throw std::invalid_argument("barrier: unknown type");
}

// ================================================================
// ============ 3. Lookback Options (Goldman-Sosin-Gatto) ============
// ================================================================
// 参考: Goldman, Sosin, Gatto (1979) J. Finance 34, 1111-1127.
//       Conze & Viswanathan (1991) J. Finance 46, 1893-1907.
//       Haug (2007) "The Complete Guide to Option Pricing Formulas" §4.5
//
// 浮动执行价 Call (payoff = S_T - m_T, m_T = min path)
// 浮动执行价 Put  (payoff = M_T - S_T, M_T = max path)
// 固定执行价 Call (payoff = (M_T - K)^+)
// 固定执行价 Put  (payoff = (K - m_T)^+)
//
// 统一记号 (b = r - q, μ = (b-σ²/2)/σ²):
//   d1     = (ln(S/X) + (b+σ²/2)T) / (σ√T)     [X = m 或 M 或 K]
//   d2     = d1 - σ√T
//   d1_ref = (ln(S/X) - (b-σ²/2)T) / (σ√T) = d1 - 2b√T/σ  [reflected drift]
//
// 浮动 Call (m ≤ S, payoff = S_T - m_T):
//   LC = S*e^{-qT}*N(d1) - m*e^{-rT}*N(d2)
//      + S*e^{-qT}*(σ²/(2b)) * [(m/S)^{2μ+1} * N(-d1_ref) - N(-d1)]
//   注: 2μ+1 = 2b/σ²
//
// 浮动 Put (M ≥ S, payoff = M_T - S_T):
//   LP = M*e^{-rT}*N(-d2) - S*e^{-qT}*N(-d1)
//      + S*e^{-qT}*(σ²/(2b)) * [N(d1) - (M/S)^{2μ+1} * N(d1_ref)]
//
// b=0 (r=q) 极限 (L'Hôpital, σ²/(2b) * bracket → finite):
//   Call: third → S*e^{-qT} * [ln(m/S)*N(γ) + σ√T*φ(γ)],  γ = (ln(m/S) - σ²T/2)/(σ√T)
//   Put:  third → S*e^{-qT} * [σ√T*φ(d1_0) - ln(M/S)*N(d1_0)],  d1_0 = (ln(S/M) + σ²T/2)/(σ√T)

// 浮动执行价 Lookback Call (payoff = S_T - m_T)
// m: 当前运行最小值 (m ≤ S, 若期权刚启动则 m = S)
inline Real lookback_call_floating(Real S, Real m, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(S - m, 0.0);
    if (S <= 0.0 || m <= 0.0 || sigma <= 0.0) throw std::invalid_argument("lookback_call: invalid params");
    if (m > S) throw std::invalid_argument("lookback_call: m must be <= S");
    Real b = r - q;
    Real sqT = sigma * std::sqrt(T);
    Real d1 = detail::bsm_d1(S, m, T, b, sigma);
    Real d2 = d1 - sqT;
    Real d1_ref = (std::log(S / m) - (b - 0.5 * sigma * sigma) * T) / sqT;
    Real base = S * std::exp(-q * T) * normal_cdf(d1) - m * std::exp(-r * T) * normal_cdf(d2);

    Real third;
    if (std::abs(b) < 1e-10) {
        // b=0 极限: γ = (ln(m/S) - σ²T/2)/(σ√T) = -d1(at b=0)
        Real gamma = (std::log(m / S) - sigma * sigma * T * 0.5) / sqT;
        third = S * std::exp(-q * T) * (std::log(m / S) * normal_cdf(gamma) + sigma * std::sqrt(T) * normal_pdf(gamma));
    } else {
        Real mu = detail::reflection_mu(b, sigma);
        Real ratio = std::pow(m / S, 2.0 * mu + 1.0);  // (m/S)^{2b/σ²}
        Real bracket = ratio * normal_cdf(-d1_ref) - normal_cdf(-d1);
        third = S * std::exp(-q * T) * (sigma * sigma / (2.0 * b)) * bracket;
    }
    return base + third;
}

// 浮动执行价 Lookback Put (payoff = M_T - S_T)
// M: 当前运行最大值 (M ≥ S, 若期权刚启动则 M = S)
inline Real lookback_put_floating(Real S, Real M, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(M - S, 0.0);
    if (S <= 0.0 || M <= 0.0 || sigma <= 0.0) throw std::invalid_argument("lookback_put: invalid params");
    if (M < S) throw std::invalid_argument("lookback_put: M must be >= S");
    Real b = r - q;
    Real sqT = sigma * std::sqrt(T);
    Real d1 = detail::bsm_d1(S, M, T, b, sigma);
    Real d2 = d1 - sqT;
    Real d1_ref = (std::log(S / M) - (b - 0.5 * sigma * sigma) * T) / sqT;
    Real base = M * std::exp(-r * T) * normal_cdf(-d2) - S * std::exp(-q * T) * normal_cdf(-d1);

    Real third;
    if (std::abs(b) < 1e-10) {
        // b=0 极限: d1_0 = (ln(S/M) + σ²T/2)/(σ√T)
        Real d1_0 = (std::log(S / M) + sigma * sigma * T * 0.5) / sqT;
        third = S * std::exp(-q * T) * (sigma * std::sqrt(T) * normal_pdf(d1_0) - std::log(M / S) * normal_cdf(d1_0));
    } else {
        Real mu = detail::reflection_mu(b, sigma);
        Real ratio = std::pow(M / S, 2.0 * mu + 1.0);  // (M/S)^{2b/σ²}
        Real bracket = normal_cdf(d1) - ratio * normal_cdf(d1_ref);
        third = S * std::exp(-q * T) * (sigma * sigma / (2.0 * b)) * bracket;
    }
    return base + third;
}

// 固定执行价 Lookback Call (payoff = (M_T - K)^+)
// M: 当前运行最大值 (M ≥ S)
//
// 若 M ≥ K: payoff = M_T - K (无条件 ITM, 因 M_T ≥ M ≥ K)
//   C_fix = lookback_put_floating(S, M, ...) + S*e^{-qT} - K*e^{-rT}
//   (由 E[M_T - K] = E[M_T - S_T] + E[S_T - K] = LP + (S*e^{bT} - K)*e^{-rT})
//
// 若 M < K (假设 M = S, 刚启动且 OTM): Conze-Viswanathan 公式
//   d1 = bsm_d1(S, K), d1_ref = d1 - 2b√T/σ
//   C_fix = S*e^{-qT}*N(d1) - K*e^{-rT}*N(d2)
//         + S*e^{-qT}*(σ²/(2b)) * [(S/K)^{-2μ} * N(-d1_ref) - N(-d1)]
inline Real lookback_call_fixed(Real S, Real K, Real M, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(M - K, 0.0);
    if (S <= 0.0 || K <= 0.0 || sigma <= 0.0) throw std::invalid_argument("lookback_call_fixed: invalid params");
    if (M < S) throw std::invalid_argument("lookback_call_fixed: M must be >= S");
    Real b = r - q;
    Real sqT = sigma * std::sqrt(T);

    if (M >= K) {
        // M_T ≥ M ≥ K ⟹ payoff = M_T - K (无条件)
        return lookback_put_floating(S, M, T, r, q, sigma)
             + S * std::exp(-q * T) - K * std::exp(-r * T);
    } else {
        // M < K (假设 M = S): Conze-Viswanathan
        Real d1 = detail::bsm_d1(S, K, T, b, sigma);
        Real d2 = d1 - sqT;
        Real d1_ref = (std::log(S / K) - (b - 0.5 * sigma * sigma) * T) / sqT;
        Real base = S * std::exp(-q * T) * normal_cdf(d1) - K * std::exp(-r * T) * normal_cdf(d2);

        Real third;
        if (std::abs(b) < 1e-10) {
            // b=0 极限: 与 floating call 相同形式, γ = (ln(S/K) - σ²T/2)/(σ√T)... 
            // 实际上 fixed call (M<K) 的 b=0 极限 = floating call (m=S) 的 b=0 极限 (用 K 代 m)
            Real gamma = (std::log(S / K) - sigma * sigma * T * 0.5) / sqT;  // = -d1 at b=0
            // 注意 ln(S/K) > 0 当 S > K; 但 M < K 意味 M=S < K, 所以 S < K, ln(S/K) < 0
            third = S * std::exp(-q * T) * (std::log(S / K) * normal_cdf(gamma) + sigma * std::sqrt(T) * normal_pdf(gamma));
            // 修正: floating call 极限用 ln(m/S), 这里 m→K, 所以用 ln(K/S) = -ln(S/K)
            // third = S*e^{-qT} * [ln(K/S)*N(γ') + σ√T*φ(γ')], γ' = (ln(K/S)-σ²T/2)/(σ√T)
            // 但 γ' = -γ - σ√T... 让我直接用正确形式:
            Real gamma2 = (std::log(K / S) - sigma * sigma * T * 0.5) / sqT;
            third = S * std::exp(-q * T) * (std::log(K / S) * normal_cdf(gamma2) + sigma * std::sqrt(T) * normal_pdf(gamma2));
        } else {
            Real mu = detail::reflection_mu(b, sigma);
            Real ratio = std::pow(S / K, -2.0 * mu);  // (S/K)^{-2μ} = (S/K)^{(σ²-2b)/σ²}
            Real bracket = ratio * normal_cdf(-d1_ref) - normal_cdf(-d1);
            third = S * std::exp(-q * T) * (sigma * sigma / (2.0 * b)) * bracket;
        }
        (void)M;  // M < K 时假设 M = S
        return base + third;
    }
}

// 固定执行价 Lookback Put (payoff = (K - m_T)^+)
// m: 当前运行最小值 (m ≤ S)
//
// 若 m ≤ K: payoff = K - m_T (无条件 ITM, 因 m_T ≤ m ≤ K)
//   P_fix = lookback_call_floating(S, m, ...) + K*e^{-rT} - S*e^{-qT}
//   (由 E[K - m_T] = E[K - S_T] + E[S_T - m_T] = (K - S*e^{bT})*e^{-rT} + LC)
//
// 若 m > K (假设 m = S, 刚启动且 OTM): Conze-Viswanathan 公式
//   d1 = bsm_d1(S, K), d1_ref = d1 - 2b√T/σ
//   P_fix = K*e^{-rT}*N(-d2) - S*e^{-qT}*N(-d1)
//         + S*e^{-qT}*(σ²/(2b)) * [N(d1) - (S/K)^{-2μ} * N(d1_ref)]
inline Real lookback_put_fixed(Real S, Real K, Real m, Real T, Real r, Real q, Real sigma) {
    if (T <= 0.0) return std::max(K - m, 0.0);
    if (S <= 0.0 || K <= 0.0 || sigma <= 0.0) throw std::invalid_argument("lookback_put_fixed: invalid params");
    if (m > S) throw std::invalid_argument("lookback_put_fixed: m must be <= S");
    Real b = r - q;
    Real sqT = sigma * std::sqrt(T);

    if (m <= K) {
        // m_T ≤ m ≤ K ⟹ payoff = K - m_T (无条件)
        return K * std::exp(-r * T) - S * std::exp(-q * T)
             + lookback_call_floating(S, m, T, r, q, sigma);
    } else {
        // m > K (假设 m = S): Conze-Viswanathan
        Real d1 = detail::bsm_d1(S, K, T, b, sigma);
        Real d2 = d1 - sqT;
        Real d1_ref = (std::log(S / K) - (b - 0.5 * sigma * sigma) * T) / sqT;
        Real base = K * std::exp(-r * T) * normal_cdf(-d2) - S * std::exp(-q * T) * normal_cdf(-d1);

        Real third;
        if (std::abs(b) < 1e-10) {
            // b=0 极限: 与 floating put 相同形式, 用 K 代 M
            // floating put 极限: σ√T*φ(d1_0) - ln(M/S)*N(d1_0), 这里 M→K
            Real d1_0 = (std::log(S / K) + sigma * sigma * T * 0.5) / sqT;
            third = S * std::exp(-q * T) * (sigma * std::sqrt(T) * normal_pdf(d1_0) - std::log(K / S) * normal_cdf(d1_0));
        } else {
            Real mu = detail::reflection_mu(b, sigma);
            Real ratio = std::pow(S / K, -2.0 * mu);  // (S/K)^{-2μ}
            Real bracket = normal_cdf(d1) - ratio * normal_cdf(d1_ref);
            third = S * std::exp(-q * T) * (sigma * sigma / (2.0 * b)) * bracket;
        }
        (void)m;  // m > K 时假设 m = S
        return base + third;
    }
}

}  // namespace v1
}  // namespace cpphub
