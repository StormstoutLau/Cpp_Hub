#pragma once
// SOURCE: Hull (2018) Ch.26 (Path-Dependent Derivatives)
// SOURCE: Glasserman (2003) Ch.3 (Path-dependent options)
// SOURCE: Goldman, Sosin, Gatto (1979) "Path dependent options: buy at the low, sell at the high"
// SOURCE: Reiner & Rubinstein (1991) "Unscrambling the Binary Code"
// 模块: 路径相关 payoff (Asian / Lookback / Barrier)
//
// 约定:
//   path 为 (n_steps+1) 长度的向量, path[0]=S0, path[n_steps]=S_T
//   Asian 观察点: 通常为 path[1..n_steps] (不含 t=0)
//   Lookback 观察点: 通常为 path[0..n_steps] (含 t=0)
//   Barrier 监控: 通常为 path[0..n_steps] (含 t=0, 连续监控的离散近似)

#include "cpphub/core/types.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <functional>

namespace cpphub {
inline namespace v1 {

// ============ Payoff 函数签名 ============
// 路径相关 payoff: 接收完整路径, 返回 payoff (未折现)
using PathPayoff = std::function<Real(const std::vector<Real>&)>;

enum class OptionType { Call, Put };

// ============ 1. Asian Options ============
enum class AsianAverageType { Arithmetic, Geometric };

// Asian option payoff
// 观察点: path[1..n] (跳过 t=0)
inline PathPayoff make_asian_payoff(Real K, OptionType opt,
                                      AsianAverageType avg = AsianAverageType::Arithmetic) {
    if (K < 0.0) throw std::invalid_argument("Asian: K must be non-negative");
    return [K, opt, avg](const std::vector<Real>& path) -> Real {
        if (path.size() < 2) throw std::invalid_argument("Asian: path too short");
        Real sum = 0.0, sum_log = 0.0;
        Size n = path.size() - 1;  // 观察点数 (跳过 t=0)
        for (Size i = 1; i <= n; ++i) {
            if (path[i] <= 0.0) throw std::invalid_argument("Asian: path must be positive");
            sum += path[i];
            sum_log += std::log(path[i]);
        }
        Real avg_price = (avg == AsianAverageType::Arithmetic)
                          ? sum / static_cast<Real>(n)
                          : std::exp(sum_log / static_cast<Real>(n));
        if (opt == OptionType::Call) {
            return std::max(avg_price - K, 0.0);
        } else {
            return std::max(K - avg_price, 0.0);
        }
    };
}

// ============ 2. Lookback Options ============
enum class LookbackType { FixedStrike, FloatingStrike };

// Lookback option payoff
// 观察点: path[0..n] (含 t=0)
inline PathPayoff make_lookback_payoff(Real K, OptionType opt,
                                         LookbackType type = LookbackType::FixedStrike) {
    if (K < 0.0) throw std::invalid_argument("Lookback: K must be non-negative");
    return [K, opt, type](const std::vector<Real>& path) -> Real {
        if (path.empty()) throw std::invalid_argument("Lookback: path empty");
        Real M = path[0], m = path[0];
        for (Size i = 1; i < path.size(); ++i) {
            M = std::max(M, path[i]);
            m = std::min(m, path[i]);
        }
        Real S_T = path.back();
        if (type == LookbackType::FixedStrike) {
            // Fixed strike: max(M - K, 0) (call), max(K - m, 0) (put)
            if (opt == OptionType::Call) {
                return std::max(M - K, 0.0);
            } else {
                return std::max(K - m, 0.0);
            }
        } else {
            // Floating strike: max(S_T - m, 0) (call), max(M - S_T, 0) (put)
            if (opt == OptionType::Call) {
                return std::max(S_T - m, 0.0);
            } else {
                return std::max(M - S_T, 0.0);
            }
        }
    };
}

// ============ 3. Barrier Options ============
enum class BarrierDirection { Up, Down };
enum class BarrierKnock { In, Out };

struct BarrierSpec {
    Real barrier;
    BarrierDirection dir;
    BarrierKnock knock;
    // inner payoff at T (vanilla call/put on S_T): max(S_T - K, 0) or max(K - S_T, 0)
    Real K;
    OptionType inner_opt = OptionType::Call;
};

// Barrier option payoff
// 监控点: path[0..n] (含 t=0, 连续监控的离散近似)
inline PathPayoff make_barrier_payoff(BarrierSpec spec) {
    if (spec.barrier <= 0.0) throw std::invalid_argument("Barrier: barrier must be positive");
    if (spec.K < 0.0) throw std::invalid_argument("Barrier: K must be non-negative");
    return [spec](const std::vector<Real>& path) -> Real {
        if (path.empty()) throw std::invalid_argument("Barrier: path empty");
        // 检查是否触及障碍
        bool knocked = false;
        for (Size i = 0; i < path.size(); ++i) {
            if (spec.dir == BarrierDirection::Up && path[i] >= spec.barrier) {
                knocked = true; break;
            }
            if (spec.dir == BarrierDirection::Down && path[i] <= spec.barrier) {
                knocked = true; break;
            }
        }
        // Up-and-Out: 触及归零; Up-and-In: 触及才有效
        bool active = (spec.knock == BarrierKnock::Out) ? !knocked : knocked;
        if (!active) return 0.0;
        Real S_T = path.back();
        return (spec.inner_opt == OptionType::Call)
                ? std::max(S_T - spec.K, 0.0)
                : std::max(spec.K - S_T, 0.0);
    };
}

// ============ 4. Vanilla payoff (用于对比测试) ============
inline PathPayoff make_vanilla_payoff(Real K, OptionType opt) {
    if (K < 0.0) throw std::invalid_argument("Vanilla: K must be non-negative");
    return [K, opt](const std::vector<Real>& path) -> Real {
        Real S_T = path.back();
        return (opt == OptionType::Call)
                ? std::max(S_T - K, 0.0)
                : std::max(K - S_T, 0.0);
    };
}

// ============ 5. Digital (binary) payoff ============
inline PathPayoff make_digital_payoff(Real K, OptionType opt, Real payoff = 1.0) {
    if (K < 0.0) throw std::invalid_argument("Digital: K must be non-negative");
    return [K, opt, payoff](const std::vector<Real>& path) -> Real {
        Real S_T = path.back();
        if (opt == OptionType::Call) {
            return (S_T > K) ? payoff : 0.0;
        } else {
            return (S_T < K) ? payoff : 0.0;
        }
    };
}

}  // namespace v1
}  // namespace cpphub
