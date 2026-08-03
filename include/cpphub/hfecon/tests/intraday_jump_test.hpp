// =============================================================================
// intraday_jump_test.hpp
// Phase 5 v1.4.3 - HFE Intraday Jump Test (Lee & Mykland 2008)
//
// 对标 R highfrequency 1.0.3:
//   intradayJumpTest(pData, volEstimator="RM", driftEstimator="none", alpha=0.95, ...)
//   — jumpTests.R L583 (RM 模式简化版, PARM 推迟 v1.4.4)
//
// 排幻觉点 (R 源码实测 2026-08-03):
//   D14: vol 调整 = sqrt(vol^2 / (lookBackPeriod-2))
//        R: spot = sqrt(sigma2hat) = sqrt(RBPVar), 然后 vol$spot = sqrt(spot^2/(K-2)) = sqrt(RBPVar/(K-2))
//        Lee-Mykland 原文无此 (K-2) 除法, R 特有调整
//   D15: Cn 无 sqrt(2/pi) 常数
//        R: Cn = sqrt(2*log(n)) - (log(pi)+log(log(n)))/(2*sqrt(2*log(n)))
//        Lee-Mykland Eq.12 有 sqrt(2/pi) 常数, R 去除使 L~N(0,1)
//   D16: n = NROW(pData) 原始观测数 (非对齐后观测数)
//
// 简化假设 (v1.4.3):
//   - 输入为单日等间隔价格 vector (无需 aggregatePrice/按日分组)
//   - 仅支持 rBPCov RM 估计器 (rMinRVar/rMedRVar 推迟 v1.4.4)
//   - drift = 0 (driftEstimator="none")
//
// 容差: 1e-8 (R 对标, 滚动窗口浮点累积)
//
// SOURCE:
//   [LM 2008] Lee & Mykland, JFE 6(5)
//   [COP 2014] Christensen, Oomen, Podolskij, JFE 144, 576-599
// R 源码: tests/fixtures/hfe/hf_src/highfrequency/R/jumpTests.R L583
//         tests/fixtures/hfe/hf_src/highfrequency/R/internalSpotVolAndDrift.R L940
//         tests/fixtures/hfe/hf_src/highfrequency/R/internalRealizedMeasures.R L256
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
// detail namespace — 内部辅助函数
// =============================================================================
namespace detail {

// -----------------------------------------------------------------------------
// rolling_rbp_var — 滚动窗口 RBPVar (Realized Bipower Variation)
// R: internalSpotVolAndDrift.R L967-982 (realizedMeasureSpotVol, RM="rBPCov")
//   RBPVar(r) = (pi/2) * sum(|r[0:n-2]| * |r[1:n-1]|)  (n-1 项)
//   滚动: 对每个大小为 lookBack 的窗口计算 RBPVar
//   窗口数 = r.size() - lookBack + 1
//
// R 源码 (internalRealizedMeasures.R L256-262):
//   RBPVar <- function(rData) {
//     returns <- as.vector(as.numeric(rData))
//     n <- length(returns)
//     rbpvar <- (pi/2) * sum(abs(returns[1:(n-1)]) * abs(returns[2:n]))
//     return(rbpvar)
//   }
// -----------------------------------------------------------------------------
inline std::vector<Real> rolling_rbp_var(const std::vector<Real>& r, int lookBack) {
    if (lookBack <= 1) {
        throw std::invalid_argument(
            "rolling_rbp_var: lookBack must be at least 2");
    }
    const Size n = r.size();
    if (n < static_cast<Size>(lookBack)) {
        throw std::invalid_argument(
            "rolling_rbp_var: r size must be >= lookBack");
    }

    const Size n_windows = n - static_cast<Size>(lookBack) + 1;
    std::vector<Real> vol(n_windows);
    const Real pi_half = std::acos(-1.0) / 2.0;

    for (Size w = 0; w < n_windows; ++w) {
        Real sum = 0.0;
        for (int i = 0; i < lookBack - 1; ++i) {
            sum += std::fabs(r[w + i]) * std::fabs(r[w + i + 1]);
        }
        vol[w] = pi_half * sum;
    }
    return vol;
}

// -----------------------------------------------------------------------------
// lee_mykland_critical_value — Lee-Mykland 临界值 (D15, D16)
// R: jumpTests.R L707-722
//   n = NROW(pData)  (D16: 原始观测数)
//   Cn = sqrt(2*log(n)) - (log(pi)+log(log(n)))/(2*sqrt(2*log(n)))  (D15: 无 sqrt(2/pi))
//   Sn = 1/sqrt(2*log(n))
//   betastar = -log(-log(1-alpha))
//   criticalValue = Cn + Sn*betastar
// -----------------------------------------------------------------------------
inline Real lee_mykland_critical_value(Size n, Real alpha) {
    if (n <= 1) {
        throw std::invalid_argument(
            "lee_mykland_critical_value: n must be > 1");
    }
    if (alpha <= 0.0 || alpha >= 1.0) {
        throw std::invalid_argument(
            "lee_mykland_critical_value: alpha must be in (0, 1)");
    }

    const Real n_real = static_cast<Real>(n);
    const Real log_n = std::log(n_real);
    const Real two_log_n = 2.0 * log_n;
    const Real sqrt_two_log_n = std::sqrt(two_log_n);

    const Real Cn = sqrt_two_log_n
        - (std::log(std::acos(-1.0)) + std::log(log_n)) / (2.0 * sqrt_two_log_n);
    const Real Sn = 1.0 / sqrt_two_log_n;
    const Real betastar = -std::log(-std::log(1.0 - alpha));

    return Cn + Sn * betastar;
}

} // namespace detail

// =============================================================================
// Intraday Jump Test Result
// =============================================================================
struct IntradayJumpTestResult {
    std::vector<Real> ztest;        // Lee-Mykland 检验统计量
    std::vector<Real> spotVol;      // D14 调整后的 spot volatility
    Real criticalValue;             // Lee-Mykland 临界值
    Size n;                         // 原始观测数 (D16)
};

// =============================================================================
// Intraday Jump Test (Lee & Mykland 2008)
// =============================================================================
//
// 检验统计量 (R highfrequency 1.0.3 intradayJumpTest 源码实测):
//   L[i] = returns[i] / spotVol[i]
//
// 其中:
//   returns[i]   = log(P[i+1] / P[i])
//   spotVol[i]   = sqrt(RBPVar_window / (lookBackPeriod - 2))   [D14]
//   RBPVar       = (pi/2) * sum(|r[j]| * |r[j+1]|)  over window
//
// 临界值 (D15, D16):
//   Cn = sqrt(2*log(n)) - (log(pi)+log(log(n)))/(2*sqrt(2*log(n)))
//   Sn = 1/sqrt(2*log(n))
//   betastar = -log(-log(1-alpha))
//   criticalValue = Cn + Sn*betastar
//   n = NROW(pData)  原始观测数
//
// 假设检验:
//   H0: 无跳跃 → L ~ N(0, 1)
//   H1: 存在跳跃 → |L| > criticalValue
//
// 简化假设 (v1.4.3):
//   - 单日等间隔价格 vector, 无需 aggregatePrice
//   - 仅 rBPCov RM 估计器
//   - drift = 0
// =============================================================================
inline IntradayJumpTestResult intraday_jump_test(
    const std::vector<Real>& pData,
    int lookBackPeriod = 10,
    double alpha = 0.95,
    const std::string& RM = "rBPCov") {

    // --- 输入校验 ---
    const Size n = pData.size();
    if (n == 0) {
        throw std::invalid_argument(
            "intraday_jump_test: pData is empty");
    }
    if (lookBackPeriod <= 2) {
        throw std::invalid_argument(
            "intraday_jump_test: lookBackPeriod must be > 2 (division by lookBackPeriod-2)");
    }
    if (RM != "rBPCov") {
        throw std::invalid_argument(
            "intraday_jump_test: only 'rBPCov' RM estimator is supported in v1.4.3");
    }
    if (n <= static_cast<Size>(lookBackPeriod)) {
        throw std::invalid_argument(
            "intraday_jump_test: pData size must be > lookBackPeriod");
    }

    // --- 计算对数收益率 ---
    // R: returns = log(PRICE) - shift(log(PRICE), type="lag"), 去除第一个 NA
    // C++: returns[i] = log(P[i+1]/P[i]), i = 0..n-2, 长度 = n-1
    std::vector<Real> returns(n - 1);
    for (Size i = 0; i + 1 < n; ++i) {
        returns[i] = std::log(pData[i + 1] / pData[i]);
    }

    // --- 滚动 RBPVar (RM="rBPCov") ---
    // R: realizedMeasureSpotVol, sigma2hat = RBPVar(window)
    auto rbp_vars = detail::rolling_rbp_var(returns, lookBackPeriod);
    // rbp_vars.size() = (n-1) - lookBack + 1 = n - lookBack

    // --- D14: spot volatility 调整 ---
    // R: vol$spot <- sqrt((vol$spot^2)/(lookBackPeriod-2))
    //    vol$spot = sqrt(sigma2hat) = sqrt(RBPVar)
    //    vol$spot^2 = RBPVar
    //    vol$spot (new) = sqrt(RBPVar / (lookBackPeriod-2))
    std::vector<Real> spotVol(rbp_vars.size());
    for (Size i = 0; i < rbp_vars.size(); ++i) {
        spotVol[i] = std::sqrt(rbp_vars[i] / (lookBackPeriod - 2));
    }

    // --- Lee-Mykland 检验统计量 ---
    // R: tests <- (returns - drift) / vol$spot, drift = 0
    // ztest[i] = returns[i + lookBack - 1] / spotVol[i]
    // (窗口末尾对应的收益率除以该窗口的 spot volatility)
    std::vector<Real> ztest(spotVol.size());
    for (Size i = 0; i < spotVol.size(); ++i) {
        if (spotVol[i] == 0.0) {
            ztest[i] = std::numeric_limits<Real>::quiet_NaN();
        } else {
            ztest[i] = returns[i + static_cast<Size>(lookBackPeriod) - 1] / spotVol[i];
        }
    }

    // --- 临界值 (D15, D16) ---
    Real cv = detail::lee_mykland_critical_value(n, alpha);

    return IntradayJumpTestResult{ztest, spotVol, cv, n};
}

} // namespace hfecon
} // namespace v1
} // namespace cpphub
