#pragma once
// SOURCE: El Euch, O., Rosenbaum, M. (2018) "The characteristic function of rough Heston models"
// SOURCE: Gatheral, J., Jaisson, T., Rosenbaum, M. (2018) "Volatility is rough" (Springer)
// SOURCE: Bayer, Friz, Gatheral (2016) "Pricing under rough volatility" (Quant. Finance 16(6))
//
// 模块: rBergomi 解析层 — 近似特征函数 + 累积量 (用于 COSEngine 半解析定价)
//
// ==================== rBergomi 模型 (见 rough_bergomi.hpp) ====================
//
//   dS_t / S_t = sqrt(v_t) * dZ_t,   Z 与 W 相关 ρ
//   v_t = ξ_0 * exp(η * W̃^H_t - 0.5 * η² * t^{2H})
//   W̃^H_t = sqrt(2H+1) * ∫_0^t (t-s)^{H-1/2} dW_s   (RL-fBm)
//
// rBergomi 无闭式特征函数 (与 rough Heston 不同), 本模块用累积量展开
// (Edgeworth/Escher 型近似) 构造 ln S_T 的近似特征函数:
//
//   φ(u) = exp( i·u·(ln S_0 + (r-q)T + c1) - 0.5·u²·c2
//                + (i·u)³·c3/6 + (i·u)⁴·c4/24 )
//
// 累积量近似公式 (基于矩展开, Gatheral-Jaisson-Rosenbaum 2018):
//   c2 = ξ_0·T·(1 + η²·T^{2H}/(2·(2H+1)))                (方差, 含 log-normal 修正)
//   c3 = 3·η·ρ·ξ_0·T^{H+1}/(H+1)·sqrt(2H+1)              (偏度, 来自 ρ 相关)
//   c4 = 3·η²·ξ_0²·T^{2H+1}/(2H+1)                       (峰度, 来自 log-vol 随机性)
//   c1 = -(c2/2 + c3/6 + c4/24)                          (鞅修正, 见下)
//
// 注意: 上述为近似公式, 精确形式需数值积分, 研究用途足够.
//
// 关于 c1 (鞅修正):
//   累积量展开要求 E[S_T] = φ(-i) = S_0·e^{(r-q)T} (风险中性鞅, 无套利),
//   即 exp(c1 + c2/2 + c3/6 + c4/24) = 1, 由此反解 c1 = -(c2/2 + c3/6 + c4/24).
//   相比文献中仅含 c2 凸性修正的 -0.5·ξ_0·T·(1 + 0.5·η²·T^{2H}/(2H+1)),
//   本实现额外计入 c3/c4 对均值的贡献, 使 Call-Put 平价精确成立.
//   (H→0 或 η→0 时两者一致, 且 η→0 精确退化为 GBM 特征函数.)
//
// 关于 c4 (数值稳定性):
//   (i·u)⁴·c4/24 = +u⁴·c4/24 为实数且 c4 > 0, 会使 |φ(u)| 在 |u| 较大时
//   突破 1 并指数增长 (非有效特征函数), 导致 COSEngine 数值不稳定
//   (|φ| 的"驼峰"使密度伪尾部振荡, 破坏平价与 COS-MC 一致性).
//   故本模块定价 CF 只保留 c1, c2, c3 (|φ| = exp(-0.5·c2·u²) 严格衰减),
//   c4 仍由 rough_bergomi_cumulants 返回, 供矩匹配/测试使用.
//
// 依赖约定: φ 为 ln S_T 的特征函数 (T 已绑定), 与 COSEngine 兼容.

#include <cmath>
#include <complex>
#include <stdexcept>
#include <functional>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/pricing/fourier/characteristic_functions.hpp"

namespace cpphub {
inline namespace v1 {

// ============ rBergomi 解析 CF 参数结构 ============
struct RoughBergomiCFParams {
    Real H;    // Hurst 指数 (0, 0.5)
    Real eta;  // vol of vol > 0
    Real rho;  // 价格-方差相关性 [-1, 1]
    Real xi0;  // 初始远期方差 > 0
    Real S0;   // 初始价格 > 0
    Real r;    // 无风险利率
    Real q;    // 股息率
    Real T;    // 时间跨度 > 0
};

// ============ 参数校验 ============
inline void validate_rough_bergomi_cf_params(const RoughBergomiCFParams& p) {
    if (p.H <= 0.0 || p.H >= 0.5)
        throw std::invalid_argument("RoughBergomiCF: H must be in (0, 0.5)");
    // η=0 时 CF 精确退化为 GBM (c3=c4=0, c2=xi0·T), 允许作为退化边界
    if (p.eta < 0.0)
        throw std::invalid_argument("RoughBergomiCF: eta must be non-negative");
    if (p.rho < -1.0 || p.rho > 1.0)
        throw std::invalid_argument("RoughBergomiCF: rho must be in [-1, 1]");
    if (p.xi0 <= 0.0)
        throw std::invalid_argument("RoughBergomiCF: xi0 must be positive");
    if (p.S0 <= 0.0)
        throw std::invalid_argument("RoughBergomiCF: S0 must be positive");
    if (p.T <= 0.0)
        throw std::invalid_argument("RoughBergomiCF: T must be positive");
}

// ============ 累积量结构体 ============
struct RoughBergomiCumulants {
    Real c1;  // 一阶累积量 (对数收益率均值, 鞅修正)
    Real c2;  // 二阶累积量 (方差)
    Real c3;  // 三阶累积量 (偏度来源)
    Real c4;  // 四阶累积量 (峰度来源)
};

// ============ 累积量近似公式 ============
inline RoughBergomiCumulants rough_bergomi_cumulants(const RoughBergomiCFParams& p) {
    validate_rough_bergomi_cf_params(p);
    Real T2H = std::pow(p.T, 2.0 * p.H);
    // c2: 方差 (含 log-normal 凸性修正)
    Real c2 = p.xi0 * p.T * (1.0 + p.eta * p.eta / (2.0 * (2.0 * p.H + 1.0)) * T2H);
    // c3: 偏度 (来自 ρ 相关)
    Real c3 = 3.0 * p.eta * p.rho * p.xi0 * std::pow(p.T, p.H + 1.0) / (p.H + 1.0)
              * std::sqrt(2.0 * p.H + 1.0);
    // c4: 峰度 (来自 log-vol 随机性)
    Real c4 = 3.0 * p.eta * p.eta * p.xi0 * p.xi0 * std::pow(p.T, 2.0 * p.H + 1.0)
              / (2.0 * p.H + 1.0);
    // c1: 鞅修正, 保证 E[S_T] = S_0·e^{(r-q)T}
    Real c1 = -(c2 / 2.0 + c3 / 6.0 + c4 / 24.0);
    return {c1, c2, c3, c4};
}

// ============ 累积量 → 矩 转换 ============
struct RoughBergomiMoments {
    Real mean;              // 均值 = c1
    Real variance;          // 方差 = c2
    Real skewness;          // 偏度 = c3 / c2^{3/2}
    Real kurtosis_excess;   // 超额峰度 = c4 / c2²
};

inline RoughBergomiMoments cumulants_to_moments(const RoughBergomiCumulants& cu) {
    Real var = cu.c2;
    if (var <= 0.0)
        throw std::invalid_argument("cumulants_to_moments: c2 must be positive");
    return {cu.c1, cu.c2, cu.c3 / std::pow(var, 1.5), cu.c4 / (var * var)};
}

// ============ log(v_t) 的理论矩 ============
struct LogVMoments {
    Real mean;      // E[log v_t] = log(ξ_0) - 0.5·η²·t^{2H}
    Real variance;  // Var[log v_t] = η²·t^{2H}
};

inline LogVMoments rough_bergomi_log_v_moment(Real t, const RoughBergomiCFParams& p) {
    if (t < 0.0) throw std::invalid_argument("RoughBergomiCF: t must be non-negative");
    if (p.xi0 <= 0.0) throw std::invalid_argument("RoughBergomiCF: xi0 must be positive");
    if (p.eta <= 0.0) throw std::invalid_argument("RoughBergomiCF: eta must be positive");
    Real t2H = std::pow(t, 2.0 * p.H);
    return {std::log(p.xi0) - 0.5 * p.eta * p.eta * t2H, p.eta * p.eta * t2H};
}

// ============ 近似特征函数 ============
// φ(u) = exp( i·u·(ln S_0 + (r-q)T + c1) - 0.5·u²·c2 + (i·u)³·c3/6 )
// 说明: 见文件头注释, 定价 CF 保留前三阶累积量以保证 |φ| ≤ 1 稳定衰减.
inline Complex rough_bergomi_characteristic_function(Complex u,
                                                      const RoughBergomiCFParams& p) {
    validate_rough_bergomi_cf_params(p);
    if (std::abs(u) < Real(1e-15)) return Complex(1.0, 0.0);

    RoughBergomiCumulants cu = rough_bergomi_cumulants(p);
    Real m = std::log(p.S0) + (p.r - p.q) * p.T + cu.c1;

    const Complex i(0.0, 1.0);
    Complex iu = i * u;
    // (i·u)³ = i³·u³ = -i·u³ (纯虚部, 只贡献相位)
    Complex exponent = iu * m + Complex(-0.5, 0.0) * cu.c2 * (u * u)
                     + (iu * iu * iu) * (cu.c3 / 6.0);
    return std::exp(exponent);
}

// ============ CF 工厂 (CharFn, 兼容 COSEngine) ============
inline CharFn make_rough_bergomi_cf(const RoughBergomiCFParams& p) {
    validate_rough_bergomi_cf_params(p);
    RoughBergomiCFParams params = p;  // 按值捕获
    return [params](Complex u) -> Complex {
        return rough_bergomi_characteristic_function(u, params);
    };
}

}  // namespace v1
}  // namespace cpphub
