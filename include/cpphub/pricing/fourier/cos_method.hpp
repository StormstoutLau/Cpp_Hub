#pragma once
// SOURCE: Fang & Oosterlee (2009) "A novel pricing method for options based on
//         Fourier-cosine series expansion" (J. Comp. Finance, 12(2), 209-227).
//
// COS 方法核心公式 (Fang-Oosterlee 2009, eq. 19):
//   v(x_0, t_0) ≈ e^{-rT} * (b-a)/2 * Σ'_{k=0}^{N-1} Re[φ(kπ/(b-a)) * e^{i kπ (x_0-a)/(b-a)}] * V_k
//
// 其中:
//   - φ(u) = E[exp(iu * ln S_T)] 为 ln S_T 的特征函数 (T 已绑定)
//   - x_0 = ln S_0 (初始状态)
//   - [a, b] 为 ln S_T 的截断区间 (典型 L=10 标准差)
//   - V_k 为 payoff 函数 v(y) 的傅里叶余弦系数
//   - Σ' 表示首项 (k=0) 权重减半
//
// Payoff 系数 (Fang-Oosterlee 2009, eq. 22-23):
//   Call: V_k = (2/(b-a)) * [χ_k(c_1, b) - K * ψ_k(c_1, b)],  c_1 = ln K
//   Put:  V_k = (2/(b-a)) * [K * ψ_k(a, c_1) - χ_k(a, c_1)],  c_1 = ln K
//
//   χ_k(s, e) = (1/(1+α²)) * [e^e * (cos(α(e-a)) + α sin(α(e-a))) - e^s * (cos(α(s-a)) + α sin(α(s-a)))]
//   ψ_k(s, e) = (sin(α(e-a)) - sin(α(s-a))) / α    (k ≥ 1);  ψ_0(s, e) = e - s
//   其中 α = kπ/(b-a)
//
// 数值稳定性:
//   - 当 e = b 时, α(b-a) = kπ, sin(kπ)=0, cos(kπ)=(-1)^k → 简化 χ_k 公式
//   - 当 s = a 时, sin(0)=0, cos(0)=1 → 进一步简化
//   - k=0 单独处理避免除零

#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/pricing/fourier/characteristic_functions.hpp"
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <utility>

namespace cpphub {
inline namespace v1 {

// ============ COS 方法引擎 ============
class COSEngine {
public:
    struct Config {
        Size n_terms = 256;       // 余弦级数项数 N (典型 128-512)
        Real L = 10.0;            // 截断区间倍数 (L 标准差, 典型 10-15)
        // 可选: 用户指定 [a, b] (覆盖自动估计)
        Real a = std::numeric_limits<Real>::quiet_NaN();
        Real b = std::numeric_limits<Real>::quiet_NaN();
    };

    // 构造: 绑定特征函数 (T 已封装在 phi 内)
    // 注意: phi 必须是 ln S_T 的 CF, 即 phi(u) = E[exp(iu * ln S_T)]
    // GCC 13 兼容: 嵌套 Config 在类体结束前不完整, 默认参数移到类外定义
    // (与 optimizer.hpp 中 LevenbergMarquardt::minimize 的处理方式一致)
    explicit COSEngine(CharFn phi, Real S0, Real r, Real q, Real T, Config cfg);

    // 欧式看涨期权价格
    Real price_call(Real K) const {
        if (K <= 0.0) throw std::invalid_argument("COSEngine: K must be positive");
        return price_vanilla(K, true);
    }

    // 欧式看跌期权价格
    Real price_put(Real K) const {
        if (K <= 0.0) throw std::invalid_argument("COSEngine: K must be positive");
        return price_vanilla(K, false);
    }

    // 批量定价: 同一到期 T, 多个行权价 K
    std::vector<Real> price_calls(const std::vector<Real>& strikes) const {
        std::vector<Real> out;
        out.reserve(strikes.size());
        for (Real K : strikes) out.push_back(price_call(K));
        return out;
    }
    std::vector<Real> price_puts(const std::vector<Real>& strikes) const {
        std::vector<Real> out;
        out.reserve(strikes.size());
        for (Real K : strikes) out.push_back(price_put(K));
        return out;
    }

    // 访问器
    Real a() const { return a_; }
    Real b() const { return b_; }
    Real S0() const { return S0_; }
    Real T() const { return T_; }
    Size n_terms() const { return cfg_.n_terms; }

private:
    CharFn phi_;
    Real S0_, r_, q_, T_;
    Config cfg_;
    Real a_, b_, x0_;

    // 计算 χ_k(s, e): 见 Fang-Oosterlee (2009) eq. 29
    // χ_k(s, e) = (1/(1+α²)) * [e^e * (cos(α(e-a)) + α sin(α(e-a)))
    //                            - e^s * (cos(α(s-a)) + α sin(α(s-a)))]
    // 其中 α = kπ/(b-a)
    static inline Real chi_k(Real s, Real e, Real a, Real alpha, Real b_minus_a) {
        (void)b_minus_a;
        Real alpha_e_a = alpha * (e - a);
        Real alpha_s_a = alpha * (s - a);
        Real denom = 1.0 + alpha * alpha;
        Real term_e = std::exp(e) * (std::cos(alpha_e_a) + alpha * std::sin(alpha_e_a));
        Real term_s = std::exp(s) * (std::cos(alpha_s_a) + alpha * std::sin(alpha_s_a));
        return (term_e - term_s) / denom;
    }

    // 计算 ψ_k(s, e): 见 Fang-Oosterlee (2009) eq. 30
    // k=0:  ψ_0(s, e) = e - s
    // k≥1:  ψ_k(s, e) = (sin(α(e-a)) - sin(α(s-a))) / α
    static inline Real psi_k(Real s, Real e, Real a, Real alpha, Size k) {
        if (k == 0) return e - s;
        Real alpha_e_a = alpha * (e - a);
        Real alpha_s_a = alpha * (s - a);
        return (std::sin(alpha_e_a) - std::sin(alpha_s_a)) / alpha;
    }

    // 核心: COS 定价
    // 数学推导:
    //   f(y|x0) ≈ (2/h) Σ' A_k cos(kπ(y-a)/h),  A_k = ∫_a^b f(y|x0) cos(kπ(y-a)/h) dy
    //   A_k ≈ Re[φ_XT(kπ/h) * e^{-i kπ a/h}]    (φ_XT = E[exp(iu X_T)], X_T=ln S_T)
    //   v(x0) = e^{-rT} ∫ v(y) f(y|x0) dy
    //         ≈ e^{-rT} Σ' A_k V_k                (V_k = (2/h) ∫ v(y) cos(kπ(y-a)/h) dy)
    //   注: V_k 已含 (2/h) 归一化, 因此最终公式无 (b-a)/2 因子
    Real price_vanilla(Real K, bool is_call) const {
        const Real h = b_ - a_;                  // 截断区间长度
        const Real c1 = std::log(K);             // 行权价对数
        const Real pi_over_h = PI / h;
        const Complex i_unit(0.0, 1.0);

        // 处理 K 超出截断区间的退化情形
        Real s_chi, e_chi, s_psi, e_psi;
        if (is_call) {
            // Call payoff: (e^y - K)^+ 对 y ∈ [max(c1, a), b] 非零
            s_chi = std::max(c1, a_);
            e_chi = b_;
            s_psi = std::max(c1, a_);
            e_psi = b_;
            if (c1 >= b_) return 0.0;  // K 太大, call 价值 ≈ 0
        } else {
            // Put payoff: (K - e^y)^+ 对 y ∈ [a, min(c1, b)] 非零
            s_chi = a_;
            e_chi = std::min(c1, b_);
            s_psi = a_;
            e_psi = std::min(c1, b_);
            if (c1 <= a_) return 0.0;  // K 太小, put 价值 ≈ 0
        }

        Real sum = 0.0;
        for (Size k = 0; k < cfg_.n_terms; ++k) {
            Real alpha_k = static_cast<Real>(k) * pi_over_h;
            // V_k = (2/h) * (chi_k - K * psi_k)  for call
            // V_k = (2/h) * (K * psi_k - chi_k)  for put
            Real chi = chi_k(s_chi, e_chi, a_, alpha_k, h);
            Real psi = psi_k(s_psi, e_psi, a_, alpha_k, k);
            Real Vk;
            if (is_call) {
                Vk = (2.0 / h) * (chi - K * psi);
            } else {
                Vk = (2.0 / h) * (K * psi - chi);
            }

            // φ_XT(kπ/h) * exp(-i kπ a/h)
            // 注意: phi_ 是 X_T=ln S_T 的 CF (不是 increment).
            // 用 phase = exp(-i kπ a/h) 等价于转换为 increment CF φ_inc = φ_XT * e^{-iu x0}
            // 再用论文的 phase = exp(i kπ (x0-a)/h). 两者数学等价.
            Real u_k = alpha_k;  // kπ/h
            Complex phi_val = phi_(Complex(u_k, 0.0));
            Complex phase = std::exp(-i_unit * alpha_k * a_);
            Complex integrand = phi_val * phase;

            Real F_k = integrand.real();  // Re[φ * phase] = A_k (余弦系数近似)

            // Σ' : 首项权重 1/2
            Real weight = (k == 0) ? 0.5 : 1.0;
            sum += weight * F_k * Vk;
        }

        // 注: V_k 已含 (2/h) 归一化, 最终公式无需 (b-a)/2 因子
        Real price = std::exp(-r_ * T_) * sum;
        return std::max(price, 0.0);  // 期权价格非负
    }
};

// 构造 (类外定义, 默认参数在 Config 完整后解析, 兼容 GCC 13)
inline COSEngine::COSEngine(CharFn phi, Real S0, Real r, Real q, Real T,
                            Config cfg = COSEngine::Config{})
    : phi_(std::move(phi)), S0_(S0), r_(r), q_(q), T_(T), cfg_(cfg) {
    if (S0 <= 0.0) throw std::invalid_argument("COSEngine: S0 must be positive");
    if (T <= 0.0) throw std::invalid_argument("COSEngine: T must be positive");
    if (cfg_.n_terms < 8) throw std::invalid_argument("COSEngine: n_terms too small (>=8 required)");
    if (cfg_.L < 1.0) throw std::invalid_argument("COSEngine: L too small (>=1.0 required)");
    // 自动估计截断区间 [a, b] (若用户未指定)
    if (std::isnan(cfg_.a) || std::isnan(cfg_.b)) {
        auto rng = cos_truncation_range(phi_, S0_, T_, cfg_.L);
        a_ = rng.first;
        b_ = rng.second;
    } else {
        a_ = cfg_.a;
        b_ = cfg_.b;
        if (b_ <= a_) throw std::invalid_argument("COSEngine: b must be > a");
    }
    x0_ = std::log(S0_);
}

// ============ 便捷工厂函数 ============
// 使用 BSM CF 的 COS 定价 (用于验证 COS 实现的正确性)
inline Real cos_call_gbm(Real S0, Real K, Real T, Real r, Real q, Real sigma,
                          Size n_terms = 256, Real L = 10.0) {
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    COSEngine::Config cfg;
    cfg.n_terms = n_terms;
    cfg.L = L;
    COSEngine engine(phi, S0, r, q, T, cfg);
    return engine.price_call(K);
}

inline Real cos_put_gbm(Real S0, Real K, Real T, Real r, Real q, Real sigma,
                         Size n_terms = 256, Real L = 10.0) {
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    COSEngine::Config cfg;
    cfg.n_terms = n_terms;
    cfg.L = L;
    COSEngine engine(phi, S0, r, q, T, cfg);
    return engine.price_put(K);
}

// 使用 Heston CF 的 COS 定价
inline Real cos_call_heston(Real S0, Real K, Real T, Real r, Real q,
                             const HestonCFParams& hp,
                             Size n_terms = 512, Real L = 12.0) {
    auto phi = make_heston_cf(S0, r, q, hp, T);
    COSEngine::Config cfg;
    cfg.n_terms = n_terms;
    cfg.L = L;
    COSEngine engine(phi, S0, r, q, T, cfg);
    return engine.price_call(K);
}

inline Real cos_put_heston(Real S0, Real K, Real T, Real r, Real q,
                            const HestonCFParams& hp,
                            Size n_terms = 512, Real L = 12.0) {
    auto phi = make_heston_cf(S0, r, q, hp, T);
    COSEngine::Config cfg;
    cfg.n_terms = n_terms;
    cfg.L = L;
    COSEngine engine(phi, S0, r, q, T, cfg);
    return engine.price_put(K);
}

}  // namespace v1
}  // namespace cpphub
