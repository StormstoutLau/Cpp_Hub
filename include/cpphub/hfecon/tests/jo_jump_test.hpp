// =============================================================================
// jo_jump_test.hpp
// Phase 5 v1.4.3 - HFE JO Jump Test (Jiang & Oomen 2008)
//
// 对标 R highfrequency 1.0.3:
//   JOjumpTest(pData, power=4, alignBy, alignPeriod, alpha=0.975)
//   — jumpTests.R L446 + internalJumpTests.R + internals.cpp L207
//
// 排幻觉点 (R 源码实测 2026-08-03):
//   D10: R = simre(pData) 简单收益率, r = makeReturns(pData) 对数收益率
//        SwV = 2*sum(R - r)  (论文用同一收益率, R 源码用两种)
//   D11: SwV = 2*sum(R-r), 不是 sum(2*(R-r)) — 数值相同但语义明确
//   D12: mu1 = 2^(6/2)*gamma(3.5)/gamma(0.5) = μ₆ = 15.0 (6 阶矩, 非 power 阶)
//   D13: rollApplyProdWrapper: m = m-1 后窗口 m 个元素, 输出长度 n-m+1
//        internals.cpp L207-216 实测: out[i] = prod(x[i:i+m]), 共 m+1 个元素
//
// 容差: 1e-10 (R 对标)
//
// SOURCE:
//   [JO 2008] Jiang & Oomen, Mathematical Finance 18(3), doi:10.1111/j.1467-9965.2008.00343.x
// R 源码: tests/fixtures/hfe/hf_src/highfrequency/R/jumpTests.R L446
//         tests/fixtures/hfe/hf_src/highfrequency/R/internalJumpTests.R (simre, RBPVar)
//         tests/fixtures/hfe/hf_src/highfrequency/src/internals.cpp L207 (rollApplyProdWrapper)
// =============================================================================
#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "cpphub/core/types.hpp"
// 复用 aj_jump_test.hpp 中的 detail::normal_cdf 和 detail::qnorm
#include "cpphub/hfecon/tests/aj_jump_test.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// detail namespace — JO jump test 特有辅助函数
// (normal_cdf / qnorm 复用自 aj_jump_test.hpp 的 detail 命名空间)
// =============================================================================
namespace detail {

// -----------------------------------------------------------------------------
// simre — 简单收益率 (D10)
// R: internalJumpTests.R L79-101
//   R[0] = 0
//   R[i] = P[i]/P[i-1] - 1
// -----------------------------------------------------------------------------
inline std::vector<Real> simre(const std::vector<Real>& pData) {
    Size n = pData.size();
    if (n == 0) return {};
    std::vector<Real> R(n);
    R[0] = 0.0;
    for (Size i = 1; i < n; ++i) {
        R[i] = pData[i] / pData[i - 1] - 1.0;
    }
    return R;
}

// -----------------------------------------------------------------------------
// make_returns — 对数收益率 (D10)
// R: makeReturns(pData)
//   r[0] = 0
//   r[i] = log(P[i]) - log(P[i-1])
// -----------------------------------------------------------------------------
inline std::vector<Real> make_returns(const std::vector<Real>& pData) {
    Size n = pData.size();
    if (n == 0) return {};
    std::vector<Real> r(n);
    r[0] = 0.0;
    for (Size i = 1; i < n; ++i) {
        r[i] = std::log(pData[i]) - std::log(pData[i - 1]);
    }
    return r;
}

// -----------------------------------------------------------------------------
// roll_apply_prod_wrapper — 滚动乘积 (D13)
// R: internals.cpp L207-216
//   m = m - 1 (关键偏移!)
//   out 长度 = n - m = n - (orig_m - 1) = n - orig_m + 1
//   out[i] = prod(x[i : i+m]), 共 m+1 = orig_m 个元素
//
// 手算验证: r = [1,2,3,4,5], orig_m = 4
//   m = 3, out_len = 2
//   out[0] = prod(r[0:3]) = 1*2*3*4 = 24
//   out[1] = prod(r[1:4]) = 2*3*4*5 = 120
// -----------------------------------------------------------------------------
inline std::vector<Real> roll_apply_prod_wrapper(const std::vector<Real>& x, int orig_m) {
    if (orig_m < 1) {
        throw std::invalid_argument("roll_apply_prod_wrapper: m must be >= 1");
    }
    int m = orig_m - 1;  // D13: m = m - 1
    int n = static_cast<int>(x.size());
    if (n <= m) {
        throw std::invalid_argument("roll_apply_prod_wrapper: n must be > m-1");
    }
    int out_len = n - m;
    std::vector<Real> out(static_cast<Size>(out_len));
    for (int i = 0; i < out_len; ++i) {
        Real prod = 1.0;
        for (int j = 0; j <= m; ++j) {  // j: 0..m, 共 m+1 = orig_m 个元素
            prod *= x[static_cast<Size>(i + j)];
        }
        out[static_cast<Size>(i)] = prod;
    }
    return out;
}

// -----------------------------------------------------------------------------
// rbp_var — 双幂变差 (RBPVar)
// R: internalRealizedMeasures.R L256-262
//   bv = (pi/2) * sum(|r[0:n-2]| * |r[1:n-1]|)
//   注意: r[0]=0, 所以第一项 |r[0]|*|r[1]| = 0
// -----------------------------------------------------------------------------
inline Real rbp_var(const std::vector<Real>& r) {
    Size n = r.size();
    if (n < 2) return 0.0;
    Real sum = 0.0;
    for (Size i = 0; i + 1 < n; ++i) {
        sum += std::fabs(r[i]) * std::fabs(r[i + 1]);
    }
    const Real pi = std::acos(-1.0);
    return (pi / 2.0) * sum;
}

// -----------------------------------------------------------------------------
// mu1_jo — 6 阶矩 μ₆ (D12)
// R: jumpTests.R L495
//   mu1 = 2^(6/2) * gamma(1/2*(6+1)) / gamma(1/2)
//       = 2^3 * gamma(3.5) / gamma(0.5)
//       = 8 * (15/8 * sqrt(pi)) / sqrt(pi)
//       = 15.0
//
// 注意: R 源码注释 "mu1" 易误解为 power 阶矩, 实际是固定的 6 阶矩 μ₆
// -----------------------------------------------------------------------------
inline Real mu1_jo() {
    return std::pow(2.0, 6.0 / 2.0) * std::tgamma((6.0 + 1.0) / 2.0) / std::tgamma(0.5);
}

} // namespace detail

// =============================================================================
// JO Jump Test Result
// =============================================================================
struct JOJumpTestResult {
    Real ztest;           // JO 检验统计量
    Real criticalLower;   // qnorm(1 - alpha)
    Real criticalUpper;   // qnorm(alpha)
    Real pvalue;          // 2 * Phi(-|z|)  (双尾)
};

// =============================================================================
// JO Jump Test (Jiang & Oomen 2008)
// =============================================================================
//
// 检验统计量 (R highfrequency 1.0.3 JOjumpTest 源码实测):
//   JO = N * bv / sqrt(av) * (1 - rv/SwV)
//
// 其中:
//   R  = simre(pData)          简单收益率 (D10)
//   r  = makeReturns(pData)    对数收益率 (D10)
//   N  = length(pData) - 1
//   bv = RBPVar(r)             双幂变差
//   rv = sum(r^2)              实现方差
//   SwV = 2 * sum(R - r)       Swap variance (D11)
//   mu1 = μ₆ = 15.0            6 阶矩 (D12, 非 power 阶)
//
// power=4:
//   q   = |rollApplyProdWrapper(r, 4)|     滚动乘积 (D13)
//   mu2 = 2^(3/4) * gamma(7/4) / gamma(1/2)
//   av  = mu1/9 * N^3 * mu2^(-4) / (N-5) * sum(q^(3/2))
//
// power=6:
//   q   = |rollApplyProdWrapper(r, 6)|     滚动乘积 (D13)
//   mu2 = 2^(1/2) * gamma(1) / gamma(1/2) = sqrt(2/pi)
//   av  = mu1/9 * N^3 * mu2^(-6) / (N-7) * sum(q)
//
// 假设检验:
//   H0: 无跳跃 → SwV ≈ rv (简单与对数收益率之差 ≈ 0)
//   H1: 存在跳跃 → SwV > rv (跳跃放大简单-对数收益差)
//   p_value = 2 * Phi(-|JO|)
//
// 注意:
//   - R 源码只支持 power=4 和 power=6, 其他值无定义
//   - rollApplyProdWrapper 的 m=m-1 偏移是 R C++ 源码的 bug-for-bug 行为 (D13)
//   - r[0]=0 导致 q[0]=0, 但不影响 sum (0^x = 0)
// =============================================================================
inline JOJumpTestResult jo_jump_test(
    const std::vector<Real>& pData,
    int power = 4,
    const std::string& alignBy = "seconds", int alignPeriod = 1,
    double alpha = 0.975) {

    // --- 输入校验 ---
    if (pData.size() < 2) {
        throw std::invalid_argument("jo_jump_test: pData must have at least 2 elements");
    }
    if (power != 4 && power != 6) {
        throw std::invalid_argument("jo_jump_test: power must be 4 or 6");
    }
    if (alignPeriod < 1) {
        throw std::invalid_argument("jo_jump_test: alignPeriod must be >= 1");
    }
    (void)alignBy;  // alignBy 仅用于接口兼容, 子采样基于 alignPeriod

    // --- 子采样 (当 alignPeriod > 1, 近似 R fastTickAggregation 对等间隔数据) ---
    std::vector<Real> data = pData;
    if (alignPeriod > 1) {
        std::vector<Real> aggregated;
        for (Size i = 0; i < data.size(); i += static_cast<Size>(alignPeriod)) {
            aggregated.push_back(data[i]);
        }
        data = std::move(aggregated);
    }

    // --- 核心计算 (D10-D13) ---
    auto R = detail::simre(data);         // 简单收益率 (D10)
    auto r = detail::make_returns(data);  // 对数收益率 (D10)
    Size n = data.size();
    Size N = n - 1;                       // R 源码: N = length(pData) - 1

    // N 必须足够大 (分母 N - power - 1 > 0)
    if (N <= static_cast<Size>(power + 1)) {
        throw std::invalid_argument("jo_jump_test: N must be > power + 1");
    }

    Real bv = detail::rbp_var(r);         // 双幂变差 (RBPVar)

    // rv = sum(r^2) — R 源码: rv <- rRVar(r); rv <- rv[[length(rv)]]
    Real rv = 0.0;
    for (Size i = 0; i < n; ++i) {
        rv += r[i] * r[i];
    }

    // SwV = 2 * sum(R - r) (D11)
    Real SwV = 0.0;
    for (Size i = 0; i < n; ++i) {
        SwV += R[i] - r[i];
    }
    SwV *= 2.0;

    if (SwV == 0.0) {
        throw std::invalid_argument("jo_jump_test: SwV is zero (degenerate input)");
    }

    // mu1 = μ₆ = 15.0 (D12)
    Real mu1 = detail::mu1_jo();

    Real N_real = static_cast<Real>(N);
    Real JOtest;

    if (power == 4) {
        // q = |roll_apply_prod_wrapper(r, 4)| (D13)
        auto q_vec = detail::roll_apply_prod_wrapper(r, 4);
        Real sum_q = 0.0;
        for (Real q : q_vec) {
            sum_q += std::pow(std::fabs(q), 6.0 / 4.0);  // q^(6/4) = q^1.5
        }

        // mu2 = 2^((6/4)/2) * gamma(1/2*(6/4+1)) / gamma(1/2)
        //      = 2^(3/4) * gamma(7/4) / gamma(1/2)
        Real mu2 = std::pow(2.0, (6.0 / 4.0) / 2.0) *
                   std::tgamma(0.5 * (6.0 / 4.0 + 1.0)) / std::tgamma(0.5);

        // av = mu1/9 * N^3 * (mu2)^(-4) / (N-5) * sum(q^(6/4))
        Real av = mu1 / 9.0 * std::pow(N_real, 3.0) * std::pow(mu2, -4.0) /
                  (N_real - 5.0) * sum_q;

        if (av <= 0.0) {
            throw std::invalid_argument("jo_jump_test: av is non-positive (power=4)");
        }

        // JOtest = N * bv / sqrt(av) * (1 - rv/SwV)
        JOtest = N_real * bv / std::sqrt(av) * (1.0 - rv / SwV);
    } else {
        // power == 6
        auto q_vec = detail::roll_apply_prod_wrapper(r, 6);
        Real sum_q = 0.0;
        for (Real q : q_vec) {
            sum_q += std::fabs(q);  // sum(q) (q^(6/6) = q^1)
        }

        // mu2 = 2^((6/6)/2) * gamma(1/2*(6/6+1)) / gamma(1/2)
        //      = 2^(1/2) * gamma(1) / gamma(1/2) = sqrt(2/pi)
        Real mu2 = std::pow(2.0, (6.0 / 6.0) / 2.0) *
                   std::tgamma(0.5 * (6.0 / 6.0 + 1.0)) / std::tgamma(0.5);

        // av = mu1/9 * N^3 * (mu2)^(-6) / (N-7) * sum(q)
        Real av = mu1 / 9.0 * std::pow(N_real, 3.0) * std::pow(mu2, -6.0) /
                  (N_real - 7.0) * sum_q;

        if (av <= 0.0) {
            throw std::invalid_argument("jo_jump_test: av is non-positive (power=6)");
        }

        JOtest = N_real * bv / std::sqrt(av) * (1.0 - rv / SwV);
    }

    // critical.value = qnorm(c(1-alpha, alpha))
    Real criticalLower = detail::qnorm(1.0 - alpha);
    Real criticalUpper = detail::qnorm(alpha);

    // pvalue = 2 * pnorm(-|JOtest|)
    Real pvalue = 2.0 * detail::normal_cdf(-std::fabs(JOtest));

    return {JOtest, criticalLower, criticalUpper, pvalue};
}

} // namespace hfecon
} // inline namespace v1
} // namespace cpphub
