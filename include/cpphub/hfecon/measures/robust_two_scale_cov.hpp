// =============================================================================
// robust_two_scale_cov.hpp
// Phase 5 v1.4.2 Wave B - Robust Two-Scale Covariance (Zhang 2011)
//
// R 对照: rRTSCov(pData, cor=FALSE, startIV=NULL, noisevar=NULL,
//                  K=300, J=1, KCov=NULL, JCov=NULL, KVar=NULL, JVar=NULL,
//                  eta=9, makePsd=FALSE)
//
// 文献: Zhang (2011), JASA 106(495), doi:10.1198/jasa.2011.tm10384
//
// SOURCE: PHASE5_HFE_SPEC §5.1, §5.3 D6
//   R highfrequency 1.0.3 R/internalRealizedMeasures.R
//     RTSRV (L490-538): 单资产 robust TSRV (迭代截断)
//     RTSCov_bi (L360-487): 双资产 robust TSCov (迭代截断)
//     cfactor_RTSCV (L28-34): eta -> ccc 查表 (30 项)
//
// 关键幻觉排除 (spec §5.3 D6):
//   eta=9 时 R 硬编码 ccc=1.0415, 不走查表 (查表值为 1.04146535666802)
//   其余 eta 走 30 项查表
//   zeta = 1/pchisq(eta, 3) (单资产), ccc = cfactor(eta) (双资产)
//   zeta 与 ccc 是不同的校正因子!
//
// 算法 (R 源码实测):
//   单资产 RTSRV (迭代 20 次):
//     logprices = log(prices)
//     nbarK = (n-K+1)/K, nbarJ = (n-J+1)/J
//     adj = 1/(1 - nbarK/nbarJ)
//     zeta = 1/pchisq(eta, 3)
//     对每个 k=1..K: logreturns_K = diff(logprices[seq(k, n, K)])
//     对每个 j=1..J: logreturns_J = diff(logprices[seq(j, n, J)])
//     vdelta_K[h] = (seconds[h+K] - seconds[h]) / secday  (等间隔时 = K/(n-1))
//     vdelta_J[h] = (seconds[h+J] - seconds[h]) / secday  (等间隔时 = J/(n-1))
//     noisevar = max(0, 1/(2*nbarJ) * (sum(lr_J^2)/J - TSRV))
//     初始: RTSRV = rMedRVar(diff(logprices[seq(1,n,K)]))
//     迭代 20 次:
//       I_K = (lr_K^2 <= eta * (RTSRV * vdelta_K + 2*noisevar))
//       I_J = (lr_J^2 <= eta * (RTSRV * vdelta_J + 2*noisevar))
//       若 sum(I_K)==0, I_K = rep(1, length)
//       若 sum(I_J)==0, I_J = rep(1, length)
//       RTSRV = adj * (zeta/K * sum(lr_K^2 * I_K)/mean(I_K)
//                      - (nbarK/nbarJ) * zeta/J * sum(lr_J^2 * I_J)/mean(I_J))
//
//   双资产 RTSCov_bi:
//     先用 RTSRV 估计各资产 IV (作为 startIV)
//     refreshTime 对齐 -> newprice1, newprice2
//     ccc = (eta==9) ? 1.0415 : cfactor(eta)
//     对每个 k=1..K: lr_K1, lr_K2, vdelta_K
//     对每个 j=1..J: lr_J1, lr_J2, vdelta_J
//     I_K1 = (lr_K1^2 <= eta*(RTSRV1*vdelta_K + 2*noisevar1))
//     I_K2 = (lr_K2^2 <= eta*(RTSRV2*vdelta_K + 2*noisevar2))
//     I_J1, I_J2 同理
//     adj_cov = n / ((K-J) * nbarK)
//     RTSCV = adj_cov * (ccc/K * sum(lr_K1 * I_K1 * lr_K2 * I_K2) / mean(I_K1 * I_K2)
//                        - (nbarK/nbarJ) * ccc/J * sum(lr_J1 * lr_J2 * I_J1 * I_J2) / mean(I_J1 * I_J2))
//
// 简化假设 (与 rTSCov 一致):
//   - 等间隔观测: seconds[i] = i, secday = n-1
//   - vdelta_K[h] = K / (n-1), vdelta_J[h] = J / (n-1)
//   - 对同步数据, refreshTime 不改变序列
// =============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/hfecon/measures/cov_utils.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

struct RobustTwoScaleCovResult {
    std::vector<Real> cov;       // d x d 协方差矩阵 (row-major)
    Size n_assets;               // 资产数 d
    Size n_obs;                  // 观测数 N
    Size K;                      // K 参数
    Size J;                      // J 参数
    Real eta;                    // eta 参数
};

class RobustTwoScaleCov {
public:
    // =====================================================================
    // 主接口: 多资产 rRTSCov
    // 输入: prices_per_asset (d 个资产的价格序列, 等长)
    // 参数:
    //   K       - 快尺度参数 (默认 300, R 默认)
    //   J       - 慢尺度参数 (默认 1, R 默认)
    //   eta     - 截断阈值 (默认 9, R 默认)
    //   make_psd_flag - 是否 PSD 投影 (默认 false)
    //   noisevar - 各资产噪声方差 (空向量表示自动估计)
    //   start_iv - 各资产初始 IV 估计 (空向量表示自动估计)
    // 异常: n < 10*K 抛 invalid_argument (R 源码硬约束)
    //       K <= J 抛 invalid_argument
    // =====================================================================
    static RobustTwoScaleCovResult estimate(
        const std::vector<std::vector<Real>>& prices_per_asset,
        Size K = 300,
        Size J = 1,
        Real eta = 9.0,
        bool make_psd_flag = false,
        std::vector<Real> noisevar = {},
        std::vector<Real> start_iv = {}) {

        if (prices_per_asset.empty()) {
            throw std::invalid_argument(
                "RobustTwoScaleCov: at least one asset required");
        }
        const Size d = prices_per_asset.size();
        const Size N = prices_per_asset[0].size();
        for (Size a = 1; a < d; ++a) {
            if (prices_per_asset[a].size() != N) {
                throw std::invalid_argument(
                    "RobustTwoScaleCov: all assets must have equal length");
            }
        }
        if (N < 10 * K) {
            throw std::invalid_argument(
                "RobustTwoScaleCov: need at least 10*K observations");
        }
        if (K <= J) {
            throw std::invalid_argument(
                "RobustTwoScaleCov: K must be greater than J");
        }

        // 计算各资产的 RTSRV (单资产 robust TSRV)
        std::vector<Real> iv(d);
        std::vector<Real> nv(d);
        for (Size a = 0; a < d; ++a) {
            Real auto_nv = 0.0;
            if (a < start_iv.size()) {
                iv[a] = start_iv[a];
                // 若 startIV 已指定, 仍需估计 noisevar
                rtsrv(prices_per_asset[a], K, J, eta, auto_nv);
            } else {
                iv[a] = rtsrv(prices_per_asset[a], K, J, eta, auto_nv);
            }
            nv[a] = (a < noisevar.size()) ? noisevar[a] : auto_nv;
        }

        std::vector<Real> cov(d * d, 0.0);

        // 对角线: RTSRV
        for (Size i = 0; i < d; ++i) {
            cov[i * d + i] = iv[i];
        }

        // 非对角线: RTSCov_bi
        for (Size i = 1; i < d; ++i) {
            for (Size j = 0; j < i; ++j) {
                const Real c = rtscov_bi(prices_per_asset[i],
                                         prices_per_asset[j],
                                         iv[i], iv[j],
                                         nv[i], nv[j],
                                         K, J, eta);
                cov[i * d + j] = c;
                cov[j * d + i] = c;
            }
        }

        if (make_psd_flag) {
            cov = make_psd(cov, d);
        }

        return RobustTwoScaleCovResult{cov, d, N, K, J, eta};
    }

    // =====================================================================
    // 单资产便捷接口
    // =====================================================================
    static Real estimate_univariate(const std::vector<Real>& prices,
                                     Size K = 300, Size J = 1,
                                     Real eta = 9.0) {
        Real nv = 0.0;
        return rtsrv(prices, K, J, eta, nv);
    }

private:
    // =====================================================================
    // chi-squared CDF for df=3: P(X <= x) for X ~ chi^2_3
    //   P(X <= x) = erf(sqrt(x/2)) - sqrt(2x/pi) * exp(-x/2)
    //   1/pchisq(eta, 3) = 1 / P(X <= eta)
    // =====================================================================
    static Real inv_pchisq3(Real x) {
        if (x <= 0.0) return std::numeric_limits<Real>::infinity();
        const Real s = std::sqrt(x / 2.0);
        const Real cdf = std::erf(s) - std::sqrt(2.0 * x / M_PI) * std::exp(-x / 2.0);
        return 1.0 / std::max(cdf, 1e-300);
    }

    // =====================================================================
    // cfactor_RTSCV: eta -> ccc 查表 (R internalRealizedMeasures.R L28-34)
    // 注意: R 源码 eta=9 时硬编码 ccc=1.0415, 不走查表
    // =====================================================================
    static Real cfactor_rtscv(Real eta) {
        // R 源码: if (eta == 9) ccc <- 1.0415 else ccc <- cfactor_RTSCV(eta)
        // 查表值是整数索引 1..30, eta=9 对应 table[9] = 1.04146535666802
        // 但 R 代码用 eta==9 ? 1.0415 : table[eta], 即硬编码值略大于查表值
        if (std::fabs(eta - 9.0) < 1e-10) {
            return 1.0415;
        }
        // 查表 (R 源码 internalRealizedMeasures.R L30-34)
        static const Real table[30] = {
            2.24524337411497, 1.51316853553965, 1.281198454641, 1.17397943904068,
            1.11618284593878, 1.08253306733799, 1.06211117548291, 1.04943239440651,
            1.04146535666802, 1.03642979300709, 1.03324015884075, 1.03121967913847,
            1.02994125279178, 1.02913375893126, 1.02862474610108, 1.02830455056679,
            1.02810353991678, 1.02797759342303, 1.02789882104643, 1.02784963478995,
            1.02781896890344, 1.02779987627651, 1.0277880042143, 1.0277806305394,
            1.02777605564142, 1.02777321996772, 1.02777146389586, 1.027770377296,
            1.02776970545685, 1.02776929035619
        };
        const int idx = static_cast<int>(std::round(eta));
        if (idx < 1 || idx > 30) {
            throw std::invalid_argument(
                "RobustTwoScaleCov: eta must be in [1, 30]");
        }
        return table[idx - 1];
    }

    // =====================================================================
    // rMedRVar: 中位数已实现方差 (R rMedRVar)
    //   rMedRVar = N * (N/(N-1)) * (pi/(6-4*sqrt(3)+pi)) * sum(|r_i| * |r_{i-1}|)
    //   这是 RTSRV 的初始 IV 估计 (R 源码 L517-518)
    //
    // R rMedRVar 实测公式 (realizedMeasures.R):
    //   z <- abs(as.matrix(rData))
    //   z <- rollApplyMedianWrapper(z)  // 中位数过滤
    //   N <- nrow(z) + 2
    //   const <- N * (N/(N-2)) * (pi/(6 - 4*sqrt(3) + pi))
    //   return(const * colSums(z^2))
    //
    // 简化: 对单资产, 直接用 sum(|r_i| * |r_{i-1}|) 近似
    // =====================================================================
    // 3 元素中位数 (避免 std::median C++20 跨编译器兼容性问题)
    static inline Real median3(Real a, Real b, Real c) {
        if (a > b) std::swap(a, b);
        if (b > c) std::swap(b, c);
        if (a > b) std::swap(a, b);
        return b;
    }

    static Real rmedrvar(const std::vector<Real>& returns) {
        const Size n = returns.size();
        if (n < 3) return 0.0;
        // R 源码: rollApplyMedianWrapper 取 3 元素中位数, 输出长度 n-2
        // 然后 N = n-2+2 = n, const = N*(N/(N-2))*(pi/(6-4*sqrt(3)+pi))
        // 这里简化: 直接用中位数过滤后的 sum^2
        std::vector<Real> med(n - 2);
        for (Size i = 1; i < n - 1; ++i) {
            med[i - 1] = median3(std::fabs(returns[i - 1]),
                                  std::fabs(returns[i]),
                                  std::fabs(returns[i + 1]));
        }
        Real sum_sq = 0.0;
        for (Size i = 0; i < med.size(); ++i) {
            sum_sq += med[i] * med[i];
        }
        const Real N = static_cast<Real>(n);
        const Real c = N * (N / (N - 2.0)) * (M_PI / (6.0 - 4.0 * std::sqrt(3.0) + M_PI));
        return c * sum_sq;
    }

    // =====================================================================
    // TSRV: 非稳健 TSRV (用于 noisevar 估计, R internalRealizedMeasures.R L603-621)
    //   adj = 1/(1 - nbarK/nbarJ)
    //   TSRV = adj * ((1/K)*sum(lr_K^2) - (nbarK/nbarJ)*(1/J)*sum(lr_J^2))
    // =====================================================================
    static Real tsrv(const std::vector<Real>& logprices, Size K, Size J) {
        const Size n = logprices.size();
        const Real nbarK = static_cast<Real>(n - K + 1) / static_cast<Real>(K);
        const Real nbarJ = static_cast<Real>(n - J + 1) / static_cast<Real>(J);
        const Real adj = 1.0 / (1.0 - nbarK / nbarJ);

        Real lr_K = 0.0;
        for (Size k = 0; k < K; ++k) {
            for (Size i = k; i + K < n; i += K) {
                const Real r = logprices[i + K] - logprices[i];
                lr_K += r * r;
            }
        }
        Real lr_J = 0.0;
        for (Size j = 0; j < J; ++j) {
            for (Size i = j; i + J < n; i += J) {
                const Real r = logprices[i + J] - logprices[i];
                lr_J += r * r;
            }
        }
        return adj * ((1.0 / K) * lr_K - (nbarK / nbarJ) * (1.0 / J) * lr_J);
    }

    // =====================================================================
    // RTSRV: 单资产 robust TSRV (R internalRealizedMeasures.R L490-538)
    //   输出: RTSRV 估计值, 同时通过引用返回 noisevar
    // =====================================================================
    static Real rtsrv(const std::vector<Real>& prices,
                      Size K, Size J, Real eta, Real& noisevar) {
        std::vector<Real> logprices(prices.size());
        for (Size i = 0; i < prices.size(); ++i) {
            logprices[i] = std::log(prices[i]);
        }
        const Size n = logprices.size();
        const Real nbarK = static_cast<Real>(n - K + 1) / static_cast<Real>(K);
        const Real nbarJ = static_cast<Real>(n - J + 1) / static_cast<Real>(J);
        const Real adj = 1.0 / (1.0 - nbarK / nbarJ);
        const Real zeta = inv_pchisq3(eta);

        // 等间隔假设: seconds[i] = i, secday = n-1
        const Real secday = static_cast<Real>(n - 1);

        // K-spaced log returns + vdelta
        std::vector<Real> lr_K, vdelta_K;
        for (Size k = 0; k < K; ++k) {
            for (Size i = k; i + K < n; i += K) {
                lr_K.push_back(logprices[i + K] - logprices[i]);
                // vdelta_K = (seconds[i+K] - seconds[i]) / secday = K / (n-1)
                vdelta_K.push_back(static_cast<Real>(K) / secday);
            }
        }

        // J-spaced log returns + vdelta
        std::vector<Real> lr_J, vdelta_J;
        for (Size j = 0; j < J; ++j) {
            for (Size i = j; i + J < n; i += J) {
                lr_J.push_back(logprices[i + J] - logprices[i]);
                vdelta_J.push_back(static_cast<Real>(J) / secday);
            }
        }

        // noisevar 估计 (R 源码 L510-512):
        //   noisevar = max(0, 1/(2*nbarJ) * (sum(lr_J^2)/J - TSRV))
        const Real tsrv_est = tsrv(logprices, K, J);
        Real sum_lr_J_sq = 0.0;
        for (Real v : lr_J) sum_lr_J_sq += v * v;
        noisevar = std::max(0.0, 1.0 / (2.0 * nbarJ) *
                            (sum_lr_J_sq / static_cast<Real>(J) - tsrv_est));

        // 初始 IV: rMedRVar(diff(logprices[seq(1, n, K)]))
        //   R 源码: sel <- seq(1, n, K); RTSRV <- rMedRVar(as.matrix(diff(logprices[sel])))
        std::vector<Real> init_returns;
        for (Size i = 0; i + K < n; i += K) {
            init_returns.push_back(logprices[i + K] - logprices[i]);
        }
        Real RTSRV = rmedrvar(init_returns);

        // 迭代 20 次
        const int max_iter = 20;
        for (int iter = 0; iter < max_iter; ++iter) {
            // I_K = (lr_K^2 <= eta * (RTSRV * vdelta_K + 2*noisevar))
            std::vector<Real> I_K(lr_K.size()), I_J(lr_J.size());
            Real sum_I_K = 0.0, sum_I_J = 0.0;
            for (Size i = 0; i < lr_K.size(); ++i) {
                const Real thresh = eta * (RTSRV * vdelta_K[i] + 2.0 * noisevar);
                I_K[i] = (lr_K[i] * lr_K[i] <= thresh) ? 1.0 : 0.0;
                sum_I_K += I_K[i];
            }
            for (Size i = 0; i < lr_J.size(); ++i) {
                const Real thresh = eta * (RTSRV * vdelta_J[i] + 2.0 * noisevar);
                I_J[i] = (lr_J[i] * lr_J[i] <= thresh) ? 1.0 : 0.0;
                sum_I_J += I_J[i];
            }

            // 若 sum == 0, I = rep(1, length)
            if (sum_I_K == 0.0) {
                std::fill(I_K.begin(), I_K.end(), 1.0);
                sum_I_K = static_cast<Real>(I_K.size());
            }
            if (sum_I_J == 0.0) {
                std::fill(I_J.begin(), I_J.end(), 1.0);
                sum_I_J = static_cast<Real>(I_J.size());
            }

            // sum(lr_K^2 * I_K) / mean(I_K)
            Real sum_lrK_sq_IK = 0.0;
            for (Size i = 0; i < lr_K.size(); ++i) {
                sum_lrK_sq_IK += lr_K[i] * lr_K[i] * I_K[i];
            }
            const Real mean_I_K = sum_I_K / static_cast<Real>(I_K.size());

            Real sum_lrJ_sq_IJ = 0.0;
            for (Size i = 0; i < lr_J.size(); ++i) {
                sum_lrJ_sq_IJ += lr_J[i] * lr_J[i] * I_J[i];
            }
            const Real mean_I_J = sum_I_J / static_cast<Real>(I_J.size());

            // RTSRV = adj * (zeta/K * sum(lr_K^2 * I_K)/mean(I_K)
            //                - (nbarK/nbarJ) * zeta/J * sum(lr_J^2 * I_J)/mean(I_J))
            RTSRV = adj * (zeta / static_cast<Real>(K) * sum_lrK_sq_IK / mean_I_K
                           - (nbarK / nbarJ) * zeta / static_cast<Real>(J) *
                             sum_lrJ_sq_IJ / mean_I_J);
        }

        return RTSRV;
    }

    // =====================================================================
    // RTSCov_bi: 双资产 robust TSCov (R internalRealizedMeasures.R L360-487)
    //   输入: prices1, prices2, startIV1, startIV2, noisevar1, noisevar2, K, J, eta
    //   输出: RTSCV 估计值
    // =====================================================================
    static Real rtscov_bi(const std::vector<Real>& prices1,
                          const std::vector<Real>& prices2,
                          Real startIV1, Real startIV2,
                          Real noisevar1, Real noisevar2,
                          Size K, Size J, Real eta) {
        // 假设已同步 (refreshTime 对同步数据不改变序列)
        std::vector<Real> logp1(prices1.size()), logp2(prices2.size());
        for (Size i = 0; i < prices1.size(); ++i) {
            logp1[i] = std::log(prices1[i]);
            logp2[i] = std::log(prices2[i]);
        }
        const Size n = logp1.size();
        const Real nbarK = static_cast<Real>(n - K + 1) / static_cast<Real>(K);
        const Real nbarJ = static_cast<Real>(n - J + 1) / static_cast<Real>(J);
        const Real adj = static_cast<Real>(n) /
            (static_cast<Real>(K - J) * nbarK);
        const Real ccc = cfactor_rtscv(eta);
        const Real secday = static_cast<Real>(n - 1);

        // K-spaced log returns + vdelta
        std::vector<Real> lr_K1, lr_K2, vdelta_K;
        for (Size k = 0; k < K; ++k) {
            for (Size i = k; i + K < n; i += K) {
                lr_K1.push_back(logp1[i + K] - logp1[i]);
                lr_K2.push_back(logp2[i + K] - logp2[i]);
                vdelta_K.push_back(static_cast<Real>(K) / secday);
            }
        }

        // J-spaced log returns + vdelta
        std::vector<Real> lr_J1, lr_J2, vdelta_J;
        for (Size j = 0; j < J; ++j) {
            for (Size i = j; i + J < n; i += J) {
                lr_J1.push_back(logp1[i + J] - logp1[i]);
                lr_J2.push_back(logp2[i + J] - logp2[i]);
                vdelta_J.push_back(static_cast<Real>(J) / secday);
            }
        }

        // I_K1, I_K2, I_J1, I_J2
        std::vector<Real> I_K1(lr_K1.size()), I_K2(lr_K2.size());
        std::vector<Real> I_J1(lr_J1.size()), I_J2(lr_J2.size());
        for (Size i = 0; i < lr_K1.size(); ++i) {
            const Real t1 = eta * (startIV1 * vdelta_K[i] + 2.0 * noisevar1);
            const Real t2 = eta * (startIV2 * vdelta_K[i] + 2.0 * noisevar2);
            I_K1[i] = (lr_K1[i] * lr_K1[i] <= t1) ? 1.0 : 0.0;
            I_K2[i] = (lr_K2[i] * lr_K2[i] <= t2) ? 1.0 : 0.0;
        }
        for (Size i = 0; i < lr_J1.size(); ++i) {
            const Real t1 = eta * (startIV1 * vdelta_J[i] + 2.0 * noisevar1);
            const Real t2 = eta * (startIV2 * vdelta_J[i] + 2.0 * noisevar2);
            I_J1[i] = (lr_J1[i] * lr_J1[i] <= t1) ? 1.0 : 0.0;
            I_J2[i] = (lr_J2[i] * lr_J2[i] <= t2) ? 1.0 : 0.0;
        }

        // sum(lr_K1 * I_K1 * lr_K2 * I_K2) / mean(I_K1 * I_K2)
        Real sum_K = 0.0, sum_IKK = 0.0;
        for (Size i = 0; i < lr_K1.size(); ++i) {
            sum_K += lr_K1[i] * I_K1[i] * lr_K2[i] * I_K2[i];
            sum_IKK += I_K1[i] * I_K2[i];
        }
        const Real mean_IKK = (sum_IKK > 0.0) ?
            (sum_IKK / static_cast<Real>(I_K1.size())) : 1.0;

        Real sum_J = 0.0, sum_IJJ = 0.0;
        for (Size i = 0; i < lr_J1.size(); ++i) {
            sum_J += lr_J1[i] * I_J1[i] * lr_J2[i] * I_J2[i];
            sum_IJJ += I_J1[i] * I_J2[i];
        }
        const Real mean_IJJ = (sum_IJJ > 0.0) ?
            (sum_IJJ / static_cast<Real>(I_J1.size())) : 1.0;

        // RTSCV = adj * (ccc/K * sum_K / mean_IKK
        //                 - (nbarK/nbarJ) * ccc/J * sum_J / mean_IJJ)
        return adj * (ccc / static_cast<Real>(K) * sum_K / mean_IKK
                      - (nbarK / nbarJ) * ccc / static_cast<Real>(J) *
                        sum_J / mean_IJJ);
    }
};

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
