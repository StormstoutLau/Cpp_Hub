// SOURCE: PHASE5_HFE_SPEC §4.5
//   [BNS 2008] Barndorff-Nielsen, Hansen, Lunde, Shephard,
//              Econometrica 76(6), 1481-1536, doi:10.1111/j.1468-0262.2008.00837.x
// R 对照: highfrequency 1.0.3 rKernelCov() → kernelEstimator()
//         realizedMeasures.cpp L77-111 (CRAN 源码实测 2026-08-02)
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/measures/kernels.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// RealizedKernel: 微结构噪声稳健 RV 估计 (BNS 2008 ECTA)
// =============================================================================
//
// 算法 (R highfrequency 1.0.3 kernelEstimator() 源码实测):
//
// 输入: r[0..n-1] (对数收益率), H (bandwidth), kernel type, DOF adj
// nab = n - 1
//
// Step 1: 计算自协方差
//   ab[h]  = Σ_{i=0}^{n-1-h} r[i] * r[i+h]    (forward lag h = γ_h)
//   ab2[h] = Σ_{i=h}^{n-1}   r[i] * r[i-h]    (backward lag h = γ_{-h})
//   对于单资产: ab[h] = ab2[h] = γ_h
//
// Step 2: 加权求和
//   ans = 0
//   for h = 0 to H:
//     w = (h == 0) ? 1.0 : KK((h-1)/H, type)     // 关键: (h-1)/H, 不是 h/H
//     theadj = adj ? n/(n-h) : 1.0                // 关键: 逐 lag 调整, 不是 n/(n-H)
//     if h == 0:
//       ans += w * theadj * γ_0
//     else:
//       ans += w * theadj * (γ_h + γ_{-h})        // 单资产: 2 * γ_h
//
// 与 BNS 2008 论文的关键差异 (重要):
//   1. 权重偏移: R 用 (h-1)/H, 论文用 h/H
//      后果: h=1 时 R w=KK(0)=1 (所有核 k(0)=1), 论文 w=KK(1/H)≠1
//   2. DOF 调整: R 用逐 lag n/(n-h), 论文用整体 n/(n-H)
//   决策: C++ 严格对标 R 源码, 保证 R baseline 数值一致 (容差 1e-12)
//
// R 对照:
//   rKernelCov(rData, kernelType, kernelParam, kernelDOFadj)
//     → RealizedKernel::estimate(returns, kernel, kernel_param, kernel_dof_adj)
//
// 容差: 1e-12 (无噪声合成数据, R 对标硬约束)
// =============================================================================

struct RealizedKernelResult {
    Real rk;              // Realized Kernel 估计值 (DOF 调整后, 若启用)
    Real rv;              // γ_0 (Realized Variance, 未调整)
    Real gamma_1;         // γ_1 (一阶自协方差, 用于噪声诊断)
    Size bandwidth;       // 实际使用的 H (= kernel_param)
    Size n_obs;           // 观测数 (收益率序列长度)
    KernelType kernel;    // 核类型
    bool dof_adjusted;    // 是否应用 DOF 调整
};

class RealizedKernel {
public:
    // 主接口: 严格对标 R rKernelCov 单资产模式
    // 输入: 日内对数收益率序列 (R rData)
    // 参数:
    //   kernel         - 核类型 (默认 Rectangular, 与 R kernelType="rectangular" 一致)
    //   kernel_param   - bandwidth H (默认 1, 与 R kernelParam=1 一致)
    //   kernel_dof_adj - 是否应用 DOF 调整 (默认 true, 与 R kernelDOFadj=TRUE 一致)
    // 异常: n < kernel_param + 1 抛 invalid_argument (需足够观测计算 γ_H)
    //       kernel_param == 0 抛 invalid_argument
    static RealizedKernelResult estimate(
        const std::vector<Real>& log_returns,
        KernelType kernel = KernelType::Rectangular,
        Size kernel_param = 1,
        bool kernel_dof_adj = true) {

        const Size n = log_returns.size();
        const Size H = kernel_param;

        if (H == 0) {
            throw std::invalid_argument(
                "RealizedKernel::estimate: kernel_param (bandwidth H) must be >= 1");
        }
        if (n < H + 1) {
            throw std::invalid_argument(
                "RealizedKernel::estimate: requires n >= H+1 returns "
                "(need at least H+1 observations to compute gamma_H)");
        }

        const Size nab = n - 1;  // R kernelEstimator: nab = na - 1

        // Step 1: 计算自协方差 ab[h] (forward) 和 ab2[h] (backward)
        // R 源码:
        //   for (j = 0; j <= lags; j++) {
        //     for(i = 0; i <= nab-j; i++)  ab[j]  += a[i] * b[i+j];
        //     for(i = j; i <= nab; i++)    ab2[j] += a[i] * b[i-j];
        //   }
        // 单资产: a = b = r, 所以 ab[h] = ab2[h] = γ_h
        std::vector<Real> ab(H + 1, 0.0);
        std::vector<Real> ab2(H + 1, 0.0);

        for (Size j = 0; j <= H; ++j) {
            // ab[j] = Σ_{i=0}^{n-1-j} r[i] * r[i+j]  (i from 0 to nab-j)
            for (Size i = 0; i + j <= nab; ++i) {
                ab[j] += log_returns[i] * log_returns[i + j];
            }
            // ab2[j] = Σ_{i=j}^{n-1} r[i] * r[i-j]  (i from j to nab)
            for (Size i = j; i <= nab; ++i) {
                ab2[j] += log_returns[i] * log_returns[i - j];
            }
        }

        // Step 2: 加权求和
        // R 源码:
        //   for (i = 0; i <= lags; i++) {
        //     w = (i == 0) ? 1 : KK((i-1)/lags, type)
        //     theadj = adj ? (nab+1)/((nab+1)-i) : 1
        //     if (i == 0) ans += w * theadj * ab[i]
        //     else        ans += w * (theadj*ab[i] + theadj*ab2[i])
        //   }
        Real ans = 0.0;
        for (Size i = 0; i <= H; ++i) {
            Real w;
            if (i == 0) {
                w = 1.0;
            } else {
                // 关键: (h-1)/H, 不是 h/H (半整数偏移)
                w = kernel_value(kernel, static_cast<Real>(i - 1) / static_cast<Real>(H));
            }

            Real theadj;
            if (!kernel_dof_adj) {
                theadj = 1.0;
            } else {
                // 关键: 逐 lag 调整 n/(n-h), 不是整体 n/(n-H)
                // R 源码: (nab+1) / ((nab+1) - i) = n / (n - i)
                theadj = static_cast<Real>(nab + 1) /
                         static_cast<Real>(nab + 1 - i);
            }

            if (i == 0) {
                ans += w * theadj * ab[0];
            } else {
                ans += w * (theadj * ab[i] + theadj * ab2[i]);
            }
        }

        RealizedKernelResult result;
        result.rk            = ans;
        result.rv            = ab[0];          // γ_0 = Realized Variance
        result.gamma_1       = ab[1];          // γ_1 (一阶自协方差)
        result.bandwidth     = H;
        result.n_obs         = n;
        result.kernel        = kernel;
        result.dof_adjusted  = kernel_dof_adj;
        return result;
    }

    // 便捷接口: 输入价格序列 (内部 make_returns)
    // make_returns 语义: 返回等长序列, 首元素 = 0 (对齐 R makeReturns)
    static RealizedKernelResult estimate_from_prices(
        const std::vector<Real>& prices,
        KernelType kernel = KernelType::Rectangular,
        Size kernel_param = 1,
        bool kernel_dof_adj = true) {
        return estimate(make_returns(prices), kernel, kernel_param, kernel_dof_adj);
    }

private:
    // make_returns: 价格 -> 对数收益率
    // R makeReturns 数值向量行为: 返回等长序列, ret[0] = 0, ret[i] = log(p[i]/p[i-1])
    static std::vector<Real> make_returns(const std::vector<Real>& prices) {
        const Size n = prices.size();
        if (n == 0) {
            throw std::invalid_argument("make_returns: empty prices");
        }
        std::vector<Real> ret(n, 0.0);
        for (Size i = 1; i < n; ++i) {
            ret[i] = std::log(prices[i] / prices[i - 1]);
        }
        return ret;
    }
};

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
