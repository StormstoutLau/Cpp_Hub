#pragma once
// SOURCE: PHASE3_SPEC §2.2 - GreeksFactory: unified entry point with method dispatch
// Auto-selection strategy (PHASE3_SPEC §2.2):
//   - Vanilla European + BSM has analytic → Analytic (closed-form, exact, fastest)
//   - Discontinuous payoff (digital/barrier) → LR (Pathwise returns 0, AAD biased)
//   - Smooth payoff + no analytic → Pathwise (lower variance than LR/AAD)
//   - Path-dependent / basket → AAD (general, but slower)
//   - Fallback → FD (centered difference, universal but slow & noisy)
//
// Design note: Cpp_Hub does not yet have Option/Model abstract base classes,
// so Factory takes BSM European parameters + PayoffType directly. When Option/
// Model abstractions are added (v1.2+), Factory can be extended to dispatch on
// model.has_analytic_greeks() etc. without breaking this API.
#include "cpphub/core/types.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"
#include "cpphub/risk/greeks/greeks_numerical.hpp"
#include "cpphub/risk/greeks/aad_greeks.hpp"
#include "cpphub/risk/greeks/pathwise_greeks.hpp"
#include "cpphub/risk/greeks/likelihood_ratio_greeks.hpp"
#include <functional>

namespace cpphub {
inline namespace v1 {

enum class GreeksMethod {
    Auto,       // 自动选择最优方法
    Analytic,   // BSM 闭式解 (最快、最准,仅适用欧式 vanilla)
    Pathwise,   // 路径法 (光滑 payoff,方差最小)
    LR,         // 似然比法 (不连续 payoff,如数字/障碍)
    FD,         // 有限差分 (通用兜底,中心差分)
    AAD         // 自动微分 (复杂模型/篮子,精确但慢)
};

enum class PayoffType {
    VanillaCall,    // max(S_T - K, 0)
    VanillaPut,     // max(K - S_T, 0)
    DigitalCall,    // 1{S_T > K}
    DigitalPut      // 1{S_T < K}
};

struct UnifiedGreeks {
    Real price = 0.0;
    Real delta = 0.0;   // d Price / d S
    Real gamma = 0.0;   // d² Price / d S²
    Real vega  = 0.0;   // d Price / d σ
    Real theta = 0.0;   // -d Price / d T
    Real rho   = 0.0;   // d Price / d r
    GreeksMethod method_used = GreeksMethod::Auto;  // 实际使用的方法
    std::string note;   // 可选诊断信息
};

class GreeksFactory {
public:
    // BSM European option Greeks — unified entry point
    // 默认 method=Auto: vanilla → Analytic, digital → LR
    static UnifiedGreeks compute_bsm(
        Real S, Real K, Real T, Real r, Real q, Real sigma,
        PayoffType payoff,
        GreeksMethod method = GreeksMethod::Auto,
        Size n_paths = 100000, uint64_t seed = 42) {

        UnifiedGreeks out;

        switch (method) {
            case GreeksMethod::Analytic:
                return compute_analytic(S, K, T, r, q, sigma, payoff);
            case GreeksMethod::Pathwise:
                return compute_pathwise(S, K, T, r, q, sigma, payoff, n_paths, seed);
            case GreeksMethod::LR:
                return compute_lr(S, K, T, r, q, sigma, payoff, n_paths, seed);
            case GreeksMethod::FD:
                return compute_fd(S, K, T, r, q, sigma, payoff);
            case GreeksMethod::AAD:
                return compute_aad(S, K, T, r, q, sigma, payoff);
            case GreeksMethod::Auto:
            default:
                return compute_auto(S, K, T, r, q, sigma, payoff, n_paths, seed);
        }
    }

private:
    static bool is_vanilla(PayoffType p) {
        return p == PayoffType::VanillaCall || p == PayoffType::VanillaPut;
    }

    static bool is_call(PayoffType p) {
        return p == PayoffType::VanillaCall || p == PayoffType::DigitalCall;
    }

    // Auto dispatch:
    //   vanilla + BSM → Analytic (closed-form available)
    //   digital       → LR (discontinuous payoff, Pathwise fails)
    static UnifiedGreeks compute_auto(
        Real S, Real K, Real T, Real r, Real q, Real sigma,
        PayoffType payoff, Size n_paths, uint64_t seed) {

        if (is_vanilla(payoff)) {
            // BSM 有闭式解,直接用 Analytic
            return compute_analytic(S, K, T, r, q, sigma, payoff);
        }
        // 数字/障碍期权 payoff 不连续,Pathwise 失效 (delta=0 a.e.),
        // AAD 在 max() kink 处有偏,LR 是首选。
        return compute_lr(S, K, T, r, q, sigma, payoff, n_paths, seed);
    }

    static UnifiedGreeks compute_analytic(
        Real S, Real K, Real T, Real r, Real q, Real sigma, PayoffType payoff) {
        UnifiedGreeks out;
        out.method_used = GreeksMethod::Analytic;
        if (!is_vanilla(payoff)) {
            out.note = "Analytic not implemented for digital; use LR/AAD";
            // Fallback to LR for digital
            return compute_lr(S, K, T, r, q, sigma, payoff, 200000, 42);
        }
        auto g = AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, is_call(payoff));
        out.price = g.price; out.delta = g.delta; out.gamma = g.gamma;
        out.vega = g.vega;   out.theta = g.theta; out.rho = g.rho;
        return out;
    }

    static UnifiedGreeks compute_pathwise(
        Real S, Real K, Real T, Real r, Real q, Real sigma,
        PayoffType payoff, Size n_paths, uint64_t seed) {
        UnifiedGreeks out;
        out.method_used = GreeksMethod::Pathwise;
        if (!is_vanilla(payoff)) {
            // Fallback to LR for discontinuous payoff
            out = compute_lr(S, K, T, r, q, sigma, payoff, n_paths, seed);
            out.note = "Pathwise not applicable to discontinuous payoff; fell back to LR";
            return out;
        }
        auto g = PathwiseGreeksEngine::bsm_european(
            S, K, T, r, q, sigma, is_call(payoff), n_paths, seed);
        out.price = g.price; out.delta = g.delta; out.vega = g.vega;
        out.gamma = 0.0;  // pathwise 二阶不直接支持
        return out;
    }

    static UnifiedGreeks compute_lr(
        Real S, Real K, Real T, Real r, Real q, Real sigma,
        PayoffType payoff, Size n_paths, uint64_t seed) {
        UnifiedGreeks out;
        out.method_used = GreeksMethod::LR;
        if (is_vanilla(payoff)) {
            auto g = LRGreeksEngine::bsm_european(
                S, K, T, r, q, sigma, is_call(payoff), n_paths, seed);
            out.price = g.price; out.delta = g.delta; out.vega = g.vega;
        } else {
            auto g = LRGreeksEngine::digital_european(
                S, K, T, r, q, sigma, is_call(payoff), n_paths, seed);
            out.price = g.price; out.delta = g.delta; out.vega = g.vega;
        }
        out.gamma = 0.0;  // LR gamma 需要 2nd-order score, O(1/T) 方差,不实现
        return out;
    }

    static UnifiedGreeks compute_aad(
        Real S, Real K, Real T, Real r, Real q, Real sigma, PayoffType payoff) {
        UnifiedGreeks out;
        out.method_used = GreeksMethod::AAD;
        if (!is_vanilla(payoff)) {
            out.note = "AAD not implemented for digital; use LR";
            return compute_lr(S, K, T, r, q, sigma, payoff, 200000, 42);
        }
        auto g = AADGreeksEngine::bsm_european(S, K, T, r, q, sigma, is_call(payoff));
        out.price = g.price; out.delta = g.delta; out.gamma = g.gamma;
        out.vega = g.vega;   out.theta = g.theta; out.rho = g.rho;
        return out;
    }

    static UnifiedGreeks compute_fd(
        Real S, Real K, Real T, Real r, Real q, Real sigma, PayoffType payoff) {
        UnifiedGreeks out;
        out.method_used = GreeksMethod::FD;
        // Price function for FD: vanilla uses BSM closed-form; digital uses closed-form too
        NumericalGreeksEngine::PriceFn price_fn;
        if (is_vanilla(payoff)) {
            price_fn = [](Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call) -> Real {
                return AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, is_call).price;
            };
        } else {
            // Digital: cash-or-nothing, 1 unit cash if ITM
            price_fn = [](Real S, Real K, Real T, Real r, Real q, Real sigma, bool is_call) -> Real {
                Real d2 = (std::log(S / K) + (r - q - 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
                Real disc = std::exp(-r * T);
                return is_call ? disc * normal_cdf(d2) : disc * normal_cdf(-d2);
            };
        }
        auto g = NumericalGreeksEngine::bsm_european(
            S, K, T, r, q, sigma, is_call(payoff), price_fn);
        out.delta = g.delta; out.gamma = g.gamma; out.vega = g.vega;
        out.theta = g.theta; out.rho = g.rho;
        out.price = price_fn(S, K, T, r, q, sigma, is_call(payoff));
        return out;
    }
};

}  // inline namespace v1
}  // namespace cpphub
