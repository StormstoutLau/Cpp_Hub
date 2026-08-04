#pragma once
// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.3 任务 1.8 - HAC 内核
// 排幻觉点: E5 (Bartlett w[l] = 1 - l/(L+1), 非 1 - l/L)
// 审计修复: Parzen K(u) = 1 - 6u^2 + 6|u|^3 (含 |u|, 对称函数)
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include "cpphub/core/constants.hpp"
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

// HAC 内核类型枚举
enum class HacKernel {
    Bartlett,           // Newey-West 1987, K(u) = 1 - |u| for |u| <= 1
    QuadraticSpectral,  // Andrews 1991, K(u) = 25/(12*pi^2*u^2) * [sin(6*pi*u/5)/(6*pi*u/5) - cos(6*pi*u/5)]
    Parzen,             // Gallant 1987, K(u) = 1 - 6u^2 + 6|u|^3 for 0 <= |u| <= 0.5; 2(1-|u|)^3 for 0.5 < |u| <= 1
    TukeyHanning        // K(u) = (1 + cos(pi*u)) / 2 for |u| <= 1
};

// 内核权重函数
// @param kernel 内核类型
// @param u 参数, u = lag / (L+1), 其中 L 是最大滞后
// @return 内核权重 K(u), 满足 K(0)=1, K(u)>=0 for |u|<=1, K(u)=0 for |u|>1
// @throws std::invalid_argument 如果 kernel 无效
Real kernel_weight(HacKernel kernel, Real u) {
    const Real au = std::fabs(u);
    switch (kernel) {
        case HacKernel::Bartlett:
            return (au <= 1.0) ? (1.0 - au) : 0.0;
        case HacKernel::QuadraticSpectral: {
            if (u == 0.0) return 1.0;  // 极限
            const Real z = 6.0 * PI * u / 5.0;
            return 25.0 / (12.0 * PI * PI * u * u) *
                   (std::sin(z) / z - std::cos(z));
        }
        case HacKernel::Parzen:
            if (au > 1.0) return 0.0;
            if (au <= 0.5) return 1.0 - 6.0 * au * au + 6.0 * au * au * au;
            return 2.0 * (1.0 - au) * (1.0 - au) * (1.0 - au);
        case HacKernel::TukeyHanning:
            return (au <= 1.0) ? ((1.0 + std::cos(PI * u)) / 2.0) : 0.0;
    }
    throw std::invalid_argument("kernel_weight: unknown HacKernel");
}

// 批量计算滞后权重向量
// @param kernel 内核类型
// @param max_lag 最大滞后 L
// @return 权重向量 w[0..L], w[0]=1.0, w[l] = kernel_weight(kernel, l/(L+1))
// 排幻觉点 E5: Bartlett w[l] = 1 - l/(L+1), 非 1 - l/L (R sandwich::kweights 实测)
std::vector<Real> kernel_weights(HacKernel kernel, Size max_lag) {
    std::vector<Real> w(max_lag + 1);
    for (Size l = 0; l <= max_lag; ++l) {
        w[l] = kernel_weight(kernel, static_cast<Real>(l) / (max_lag + 1));
    }
    return w;
}

// 最优滞后选择 (Andrews 1991 自动带宽)
// @param n_obs 样本数 T
// @param kernel 内核类型
// @param ar1_coef AR(1) 系数 alpha(1) (可选, 默认 0)
// @return 最优最大滞后 L = floor(1.1447 * (alpha(1) * T)^(1/3)) (Andrews 1991 for Bartlett)
//        对于 QS 内核, 返回最优带宽 b* = 1.3221 * (alpha(2) * T)^(1/5) (Andrews 1991)
// 排幻觉点 E4: 默认自动带宽基于 AR(1) 拟合的 Andrews (1991) 公式, 非 NW 1987 经验法则 floor(4*(T/100)^(2/9))
Size select_max_lag(Size n_obs, HacKernel kernel, Real ar1_coef = 0.0) {
    const Real T = static_cast<Real>(n_obs);
    if (ar1_coef == 0.0) {
        // 退化 (Alpha=0), fallback 到 NW 经验法则
        return static_cast<Size>(std::floor(4.0 * std::pow(T / 100.0, 2.0 / 9.0)));
    }
    const Real rho = ar1_coef;
    const Real rho2 = rho * rho;
    if (kernel == HacKernel::QuadraticSpectral) {
        // alpha(2) = 4*rho^2 / (1-rho^2)^4 (Andrews 1991)
        const Real alpha2 = 4.0 * rho2 / std::pow(1.0 - rho2, 4);
        return static_cast<Size>(std::floor(1.3221 * std::pow(alpha2 * T, 1.0 / 5.0)));
    }
    // Bartlett (及其余截断内核): alpha(1) = 4*rho^2 / (1-rho^2)^2
    const Real alpha1 = 4.0 * rho2 / ((1.0 - rho2) * (1.0 - rho2));
    return static_cast<Size>(std::floor(1.1447 * std::pow(alpha1 * T, 1.0 / 3.0)));
}

// 内核名称字符串
std::string to_string(HacKernel kernel) {
    switch (kernel) {
        case HacKernel::Bartlett: return "Bartlett";
        case HacKernel::QuadraticSpectral: return "QuadraticSpectral";
        case HacKernel::Parzen: return "Parzen";
        case HacKernel::TukeyHanning: return "TukeyHanning";
    }
    throw std::invalid_argument("to_string: unknown HacKernel");
}

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub