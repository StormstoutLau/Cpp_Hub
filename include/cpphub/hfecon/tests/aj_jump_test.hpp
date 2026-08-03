// =============================================================================
// aj_jump_test.hpp
// Phase 5 v1.4.3 - HFE AJ Jump Test (Aït-Sahalia & Jacod 2009)
//
// 对标 R highfrequency 1.0.3:
//   AJjumpTest(pData, p=4, k=2, alignBy, alignPeriod, alphaMultiplier=4, alpha=0.975)
//   — jumpTests.R L106 + internalJumpTests.R
//
// 排幻觉点 (R 源码实测 2026-08-03):
//   D4: alpha 动态 = alphaMultiplier * sqrt(RV), 覆盖入参 alpha
//       (R bug: critical value 用动态 alpha, qnorm(c(1-alpha_dyn, alpha_dyn)))
//   D5: seq(1,N,h) 整数步长抽样, h = alignPeriod * scale(alignBy)
//   D6: rse = abs(makeReturns(pData[selection])), 过滤价格后重算收益 (非 r[selection])
//       R 回收: pData[N+1] 复用 selection[1] (R 1-based)
//   D7: Ap = (1/N)^(1-p/2)/mup * sum(rse^p)
//   D8: calculateNpk 含 fmupk(p,k) 查表项
//   D9: fmupk 硬编码表 (p=2,3,4 × k=2,3,4), 其他 (p,k) 走 MC (v1.4.3 未实现, 抛异常)
//
// 容差: 1e-10 (R 对标)
//
// SOURCE:
//   [AS-J 2009] Aït-Sahalia & Jacod, Annals of Statistics 37(1), 184-222
// R 源码: tests/fixtures/hfe/hf_src/highfrequency/R/jumpTests.R L106
//         tests/fixtures/hfe/hf_src/highfrequency/R/internalJumpTests.R
// =============================================================================
#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// detail namespace — 内部辅助函数 (对标 R internalJumpTests.R)
// =============================================================================
namespace detail {

// -----------------------------------------------------------------------------
// fmupk — R 硬编码查表 (D9)
// R: internalJumpTests.R L18-48
// 表: (p=2,3,4) × (k=2,3,4), 其他 (p,k) 走 MC (v1.4.3 未实现)
// -----------------------------------------------------------------------------
inline Real fmupk(Real p, Real k) {
    const int pi = static_cast<int>(p);
    const int ki = static_cast<int>(k);

    if (pi == 2) {
        switch (ki) {
            case 2: return 4.00;
            case 3: return 5.00;
            case 4: return 6.00;
        }
    }
    if (pi == 3) {
        switch (ki) {
            case 2: return 24.07;
            case 3: return 33.63;
            case 4: return 43.74;
        }
    }
    if (pi == 4) {
        switch (ki) {
            case 2: return 204.04;
            case 3: return 320.26;
            case 4: return 455.67;
        }
    }
    // R fallback: MC mukp(p,k,t=1e6) × 100 rep, round 2 位
    // v1.4.3 未实现 MC, 抛异常
    throw std::invalid_argument(
        "fmupk: (p,k) not in table {2,3,4}x{2,3,4}; "
        "MC fallback not implemented in v1.4.3");
}

// -----------------------------------------------------------------------------
// mup — 标准正态 p 阶绝对矩
// R: mup = 2^(p/2) * gamma((p+1)/2) / gamma(0.5)
//    internalJumpTests.R L53 (calculateNpk 内联)
//    internalJumpTests.R L61 (calculateV 内联)
// -----------------------------------------------------------------------------
inline Real mup(Real p) {
    return std::pow(2.0, p / 2.0) * std::tgamma((p + 1.0) / 2.0) / std::tgamma(0.5);
}

// -----------------------------------------------------------------------------
// mu2p — 标准正态 2p 阶绝对矩
// R: mu2p = 2^p * gamma((2p+1)/2) / gamma(0.5)
//    internalJumpTests.R L54 (calculateNpk: 2^((2*p)/2) = 2^p)
//    internalJumpTests.R L62 (calculateV: 2^p)
// -----------------------------------------------------------------------------
inline Real mu2p(Real p) {
    return std::pow(2.0, p) * std::tgamma((2.0 * p + 1.0) / 2.0) / std::tgamma(0.5);
}

// -----------------------------------------------------------------------------
// calculate_npk — N(p,k) 系数 (D8)
// R: internalJumpTests.R L52-56
//   npk = (1/mup^2) * (k^(p-2)*(1+k)*mu2p + k^(p-2)*(k-1)*mup^2
//                       - 2*k^(p/2-1)*fmupk(p,k))
// -----------------------------------------------------------------------------
inline Real calculate_npk(Real p, Real k) {
    const Real mup_val = mup(p);
    const Real mu2p_val = mu2p(p);
    const Real fmupk_val = fmupk(p, k);

    const Real kp2 = std::pow(k, p - 2.0);           // k^(p-2)
    const Real khalf = std::pow(k, p / 2.0 - 1.0);   // k^(p/2-1)

    return (1.0 / (mup_val * mup_val)) *
           (kp2 * (1.0 + k) * mu2p_val +
            kp2 * (k - 1.0) * mup_val * mup_val -
            2.0 * khalf * fmupk_val);
}

// -----------------------------------------------------------------------------
// calculate_v — V 统计量 (D7)
// R: internalJumpTests.R L60-66
//   Ap  = (1/N)^(1-p/2) / mup * sum(rse^p)
//   A2p = (1/N)^(1-p)   / mu2p * sum(rse^(2p))
//   V   = calculateNpk(p,k) * A2p / (N * Ap^2)
// -----------------------------------------------------------------------------
inline Real calculate_v(const std::vector<Real>& rse, Real p, Real k, Size N) {
    const Real mup_val = mup(p);
    const Real mu2p_val = mu2p(p);
    const Real N_real = static_cast<Real>(N);

    Real sum_rse_p = 0.0;
    Real sum_rse_2p = 0.0;
    for (const Real r : rse) {
        sum_rse_p += std::pow(r, p);
        sum_rse_2p += std::pow(r, 2.0 * p);
    }

    const Real Ap = std::pow(1.0 / N_real, 1.0 - p / 2.0) / mup_val * sum_rse_p;
    const Real A2p = std::pow(1.0 / N_real, 1.0 - p) / mu2p_val * sum_rse_2p;

    if (Ap == 0.0) {
        return 0.0;  // 避免除零
    }

    const Real npk = calculate_npk(p, k);
    return npk * A2p / (N_real * Ap * Ap);
}

// -----------------------------------------------------------------------------
// normal_cdf — 标准正态 CDF Phi(x)
// Phi(x) = 0.5 * (1 + erf(x / sqrt(2)))
// -----------------------------------------------------------------------------
inline Real normal_cdf(Real x) noexcept {
    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

// -----------------------------------------------------------------------------
// qnorm — 逆标准正态 CDF (Acklam's algorithm + Newton 精炼, 精度 ~1e-12)
// R qnorm(0.975) = 1.9599639845400536
// 对 alpha_dyn 超出 (0,1) 返回 NaN/Inf (匹配 R qnorm 行为)
// -----------------------------------------------------------------------------
inline Real qnorm(Real p) noexcept {
    if (p <= 0.0 || p >= 1.0) {
        if (p == 0.0) return -std::numeric_limits<Real>::infinity();
        if (p == 1.0) return  std::numeric_limits<Real>::infinity();
        return std::numeric_limits<Real>::quiet_NaN();
    }

    // Acklam's coefficients
    static const Real a[] = {
        -3.969683028665376e+01,  2.209460984245205e+02,
        -2.759285104469687e+02,  1.383577518672690e+02,
        -3.066479806614716e+01,  2.506628277459239e+00
    };
    static const Real b[] = {
        -5.447609879822406e+01,  1.615858368580409e+02,
        -1.556989798598866e+02,  6.680131188771972e+01,
        -1.328068155288572e+01
    };
    static const Real c[] = {
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00,  2.938163982698783e+00
    };
    static const Real d[] = {
         7.784695709041462e-03,  3.224671290700398e-01,
         2.445134137142996e+00,  3.754408661907416e+00
    };
    static const Real plow  = 0.02425;
    static const Real phigh = 1.0 - plow;

    Real x;
    if (p < plow) {
        Real q = std::sqrt(-2.0 * std::log(p));
        x = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
            ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    } else if (p <= phigh) {
        Real q = p - 0.5;
        Real r = q * q;
        x = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q /
            (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);
    } else {
        Real q = std::sqrt(-2.0 * std::log(1.0 - p));
        x = -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
             ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    }

    // Newton 精炼一步
    Real e = normal_cdf(x) - p;
    Real u = e * std::sqrt(2.0 * 3.14159265358979323846) * std::exp(x * x / 2.0);
    x = x - u / (1.0 + x * u / 2.0);
    return x;
}

} // namespace detail

// =============================================================================
// AJ Jump Test Result
// =============================================================================
struct AJJumpTestResult {
    Real ztest;           // AJ 检验统计量
    Real criticalLower;   // qnorm(1 - alpha_dyn)  [D4: 动态 alpha]
    Real criticalUpper;   // qnorm(alpha_dyn)
    Real pvalue;          // 2 * Phi(-|z|)  (双尾)
};

// =============================================================================
// AJ Jump Test (Aït-Sahalia & Jacod 2009)
// =============================================================================
//
// 检验统计量 (R highfrequency 1.0.3 AJjumpTest 源码实测):
//   AJ = (S - k^(p/2-1)) / sqrt(V)
//
// 其中:
//   S = PV(p, kh) / PV(p, h)  = sum(|r_slow|^p) / sum(|r_fast|^p)
//   V = calculateNpk(p,k) * A2p / (N * Ap^2)
//
// 假设检验:
//   H0: 无跳跃 → S → k^(p/2-1)
//   H1: 存在跳跃 → S → 1
//   p_value = 2 * Phi(-|AJ|)
//
// 注意 (D4): alpha 在函数内部被覆盖为 alphaMultiplier * sqrt(RV),
//            critical value 使用此动态 alpha (R bug, bug-for-bug 复现)
// =============================================================================
inline AJJumpTestResult aj_jump_test(
    const std::vector<Real>& pData,
    int p = 4, int k = 2,
    const std::string& alignBy = "seconds", int alignPeriod = 1,
    double alphaMultiplier = 4.0, double alpha = 0.975) {

    // --- 输入校验 ---
    const Size n_prices = pData.size();
    if (n_prices < 2) {
        throw std::invalid_argument(
            "aj_jump_test: pData must have at least 2 prices");
    }
    if (alignPeriod <= 0) {
        throw std::invalid_argument(
            "aj_jump_test: alignPeriod must be positive");
    }

    // scale(alignBy): seconds=1, minutes=60, hours=3600
    // R scale() 仅支持 "seconds"/"minutes"/"hours", 其他返回 NULL (导致后续报错)
    int scale_factor;
    if (alignBy == "seconds") {
        scale_factor = 1;
    } else if (alignBy == "minutes") {
        scale_factor = 60;
    } else if (alignBy == "hours") {
        scale_factor = 3600;
    } else {
        throw std::invalid_argument(
            "aj_jump_test: alignBy must be 'seconds', 'minutes', or 'hours'");
    }

    const Real p_real = static_cast<Real>(p);
    const Real k_real = static_cast<Real>(k);

    // --- makeReturns: 对数收益率 ---
    // R: makeReturns(pData) = diff(log(pData)), 长度 N
    const Size N = n_prices - 1;
    std::vector<Real> r(N);
    for (Size i = 0; i < N; ++i) {
        r[i] = std::fabs(std::log(pData[i + 1]) - std::log(pData[i]));
    }

    // --- RV (realized variance) ---
    Real rv = 0.0;
    for (Size i = 0; i < N; ++i) {
        rv += r[i] * r[i];
    }

    // --- D4: 动态 alpha = alphaMultiplier * sqrt(RV) ---
    // R: alpha <- alphaMultiplier * sqrt(rCov(pData, makeReturns=TRUE))
    //    rCov(makeReturns=TRUE) = sum(r^2) = RV
    const Real alpha_dyn = alphaMultiplier * std::sqrt(rv);

    // --- cvalue = alpha * (1/N)^w, w = 0.47 ---
    const Real w = 0.47;
    const Real cvalue = alpha_dyn * std::pow(1.0 / static_cast<Real>(N), w);

    // --- D5: 整数步长抽样 ---
    // R: seq1 = seq(1, N, h), seq2 = seq(1, N, h*k)
    //    h = alignPeriod * scale(alignBy)
    // C++ (0-based): seq1 = {0, h, 2h, ...} while idx < N
    const int h = alignPeriod * scale_factor;
    const int hk = h * k;

    std::vector<Size> seq1;
    for (int idx = 0; idx < static_cast<int>(N); idx += h) {
        seq1.push_back(static_cast<Size>(idx));
    }

    std::vector<Size> seq2;
    for (int idx = 0; idx < static_cast<int>(N); idx += hk) {
        seq2.push_back(static_cast<Size>(idx));
    }

    // --- 子采样价格 ---
    std::vector<Real> pData1, pData2;
    pData1.reserve(seq1.size());
    for (Size idx : seq1) {
        pData1.push_back(pData[idx]);
    }
    pData2.reserve(seq2.size());
    for (Size idx : seq2) {
        pData2.push_back(pData[idx]);
    }

    // --- 子采样收益率 + 幂变差 ---
    // R: r1 = abs(makeReturns(pData1)), pv1 = sum(r1^p)
    Real pv1 = 0.0, pv2 = 0.0;
    if (pData1.size() >= 2) {
        for (Size i = 0; i + 1 < pData1.size(); ++i) {
            Real r1 = std::fabs(std::log(pData1[i + 1]) - std::log(pData1[i]));
            pv1 += std::pow(r1, p_real);
        }
    }
    if (pData2.size() >= 2) {
        for (Size i = 0; i + 1 < pData2.size(); ++i) {
            Real r2 = std::fabs(std::log(pData2[i + 1]) - std::log(pData2[i]));
            pv2 += std::pow(r2, p_real);
        }
    }

    // S = pv2 / pv1
    Real S = 0.0;
    if (pv1 > 0.0) {
        S = pv2 / pv1;
    }

    // --- D6: rse = abs(makeReturns(pData[selection])) ---
    // R: selection = abs(r) < cvalue (r 已是 abs, 故 selection = r < cvalue)
    //    rse = abs(makeReturns(pData[selection]))
    //    pData[selection]: R 回收 selection (长度 N) 到 pData (长度 N+1)
    //    pData[N+1] (R 1-based) = pData[N] (C++ 0-based) 复用 selection[1] = selection[0]
    std::vector<Real> filtered_prices;
    for (Size i = 0; i < N; ++i) {
        if (r[i] < cvalue) {
            filtered_prices.push_back(pData[i]);
        }
    }
    // R 回收: pData[N] 复用 selection[0]
    if (N > 0 && r[0] < cvalue) {
        filtered_prices.push_back(pData[N]);
    }

    // rse = abs(makeReturns(filtered_prices))
    std::vector<Real> rse;
    if (filtered_prices.size() >= 2) {
        rse.reserve(filtered_prices.size() - 1);
        for (Size i = 0; i + 1 < filtered_prices.size(); ++i) {
            rse.push_back(std::fabs(
                std::log(filtered_prices[i + 1]) - std::log(filtered_prices[i])));
        }
    }

    // --- V = calculateV(rse, p, k, N) ---
    Real V = 0.0;
    if (!rse.empty()) {
        V = detail::calculate_v(rse, p_real, k_real, N);
    }

    // --- AJ test statistic ---
    // AJ = (S - k^(p/2-1)) / sqrt(V)
    Real AJtest = 0.0;
    const Real k_term = std::pow(k_real, p_real / 2.0 - 1.0);
    if (V > 0.0) {
        AJtest = (S - k_term) / std::sqrt(V);
    }

    // --- critical value (D4: 使用动态 alpha, R bug-for-bug 复现) ---
    // R: qnorm(c(1-alpha, alpha)) where alpha = alpha_dyn
    const Real crit_lower = detail::qnorm(1.0 - alpha_dyn);
    const Real crit_upper = detail::qnorm(alpha_dyn);

    // --- p-value: 2 * Phi(-|AJ|) ---
    const Real pvalue = 2.0 * detail::normal_cdf(-std::fabs(AJtest));

    AJJumpTestResult result;
    result.ztest = AJtest;
    result.criticalLower = crit_lower;
    result.criticalUpper = crit_upper;
    result.pvalue = pvalue;
    return result;
}

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
