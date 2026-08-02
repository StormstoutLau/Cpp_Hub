// =============================================================================
// preaveraged_cov.hpp
// Phase 5 v1.4.2 Wave B - Subsampled (Averaged) Realized Covariance
//
// R 对照: rAVGCov(rData, cor=FALSE, alignBy="minutes", alignPeriod=5, k=1,
//                 makeReturns=FALSE, ...)
//
// 文献: Liu, Patton, Sheppard (2015), J. Econometrics 187, 293-311,
//       doi:10.1016/j.jeconom.2015.02.005
//       Zhang, Mykland, Aït-Sahalia (2005), JASA 100(472), 1394-1411,
//       doi:10.1198/016214505000000548
//
// SOURCE: PHASE5_HFE_SPEC §5.1, §5.3 D2
//   R highfrequency 1.0.3 R/realizedMeasures.R rAVGCov (L790-951)
//
// 关键幻觉排除 (spec §5.3 D2):
//   单资产加 (m+1)/m 系数 (m = alignPeriod/k = scalingFraction)
//   多资产不加该校正系数
//   scalingFraction = alignPeriod/k 必须为正整数
//
// 算法 (R 源码等间隔数据等价实现):
//   scalingFraction s = alignPeriod / k
//   单资产:
//     主项 (j=0): 按 s 分组求和 -> 平方和
//     偏移项 (j=1..s-1): 去掉前 j 和后 s-j 个观测,
//                        按 s 分组求和 -> 平方和,
//                        乘以校正因子 (nrow/frac + 1) / (nrow/frac)
//                        其中 nrow = N - s, frac = s
//     RV = (主项 + sum(偏移项)) / s
//   多资产:
//     对角线: 单资产 rAVGCov
//     非对角线: subsampled cross-product, 不加校正
//
// 注意: R 实现使用 data.table 时间戳聚合, 本实现针对等间隔数据简化,
//       使用整数索引分组 (数学等价于 R 的 ceiling 时间聚合对等间隔数据).
// =============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/measures/cov_utils.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

struct PreaveragedCovResult {
    std::vector<Real> cov;       // d x d 协方差矩阵 (row-major)
    Size n_assets;               // 资产数 d
    Size n_obs;                  // 观测数 N
    Size scaling_fraction;       // s = alignPeriod / k
};

class PreaveragedCov {
public:
    // =====================================================================
    // 主接口: 多资产 rAVGCov
    // 输入: prices_per_asset (d 个资产的价格序列, 等长)
    // 参数:
    //   align_period - 慢尺度 (默认 5, R 默认 alignPeriod=5)
    //   k            - 快尺度 (默认 1, R 默认 k=1)
    //   cor          - 是否返回相关系数矩阵 (默认 false)
    //   make_psd_flag- 是否 PSD 投影 (默认 false, R rAVGCov 无 makePsd 参数,
    //                  但保持接口一致性)
    // 异常: align_period % k != 0 抛 invalid_argument
    //       align_period / k <= 0 抛 invalid_argument
    //       价格序列长度不一致或为空抛 invalid_argument
    // =====================================================================
    static PreaveragedCovResult estimate(
        const std::vector<std::vector<Real>>& prices_per_asset,
        Size align_period = 5,
        Size k = 1,
        bool cor = false,
        bool make_psd_flag = false) {

        if (prices_per_asset.empty()) {
            throw std::invalid_argument(
                "PreaveragedCov: at least one asset required");
        }
        if (k == 0) {
            throw std::invalid_argument(
                "PreaveragedCov: k must be positive");
        }
        if (align_period % k != 0) {
            throw std::invalid_argument(
                "PreaveragedCov: alignPeriod/k must be an integer");
        }
        const Size s = align_period / k;
        if (s == 0) {
            throw std::invalid_argument(
                "PreaveragedCov: alignPeriod must be >= k");
        }

        const Size d = prices_per_asset.size();
        const Size N = prices_per_asset[0].size();
        for (Size a = 1; a < d; ++a) {
            if (prices_per_asset[a].size() != N) {
                throw std::invalid_argument(
                    "PreaveragedCov: all assets must have equal length");
            }
        }
        if (N < s + 1) {
            throw std::invalid_argument(
                "PreaveragedCov: need at least s+1 observations");
        }

        std::vector<Real> cov(d * d, 0.0);

        // 对角线: 单资产 rAVGCov
        for (Size i = 0; i < d; ++i) {
            cov[i * d + i] = ravg_univariate(prices_per_asset[i], s);
        }

        // 非对角线: subsampled cross-product
        if (d > 1) {
            for (Size i = 1; i < d; ++i) {
                for (Size j = 0; j < i; ++j) {
                    const Real c = ravg_bivariate(
                        prices_per_asset[i], prices_per_asset[j], s);
                    cov[i * d + j] = c;
                    cov[j * d + i] = c;
                }
            }
        }

        // PSD 投影
        if (!cor && make_psd_flag && d > 1) {
            cov = make_psd(cov, d);
        }

        // 相关矩阵
        if (cor) {
            std::vector<Real> diag_sqrt(d);
            for (Size i = 0; i < d; ++i) {
                diag_sqrt[i] = std::sqrt(cov[i * d + i]);
            }
            for (Size i = 0; i < d; ++i) {
                for (Size j = 0; j < d; ++j) {
                    if (diag_sqrt[i] > 1e-300 && diag_sqrt[j] > 1e-300) {
                        cov[i * d + j] /= (diag_sqrt[i] * diag_sqrt[j]);
                    } else {
                        cov[i * d + j] = 0.0;
                    }
                }
            }
            for (Size i = 0; i < d; ++i) cov[i * d + i] = 1.0;
            if (make_psd_flag && d > 1) {
                cov = make_psd(cov, d);
                for (Size i = 0; i < d; ++i) cov[i * d + i] = 1.0;
            }
        }

        return PreaveragedCovResult{cov, d, N, s};
    }

    // =====================================================================
    // 单资产便捷接口
    // =====================================================================
    static Real estimate_univariate(const std::vector<Real>& prices,
                                     Size align_period = 5,
                                     Size k = 1) {
        if (k == 0 || align_period % k != 0) {
            throw std::invalid_argument(
                "PreaveragedCov: alignPeriod/k must be a positive integer");
        }
        return ravg_univariate(prices, align_period / k);
    }

private:
    // =====================================================================
    // make_returns: log price diff
    // =====================================================================
    static std::vector<Real> make_returns(const std::vector<Real>& prices) {
        if (prices.size() < 2) return {};
        std::vector<Real> ret(prices.size() - 1);
        for (Size i = 0; i + 1 < prices.size(); ++i) {
            ret[i] = std::log(prices[i + 1]) - std::log(prices[i]);
        }
        return ret;
    }

    // =====================================================================
    // bucket_sum_sq: 按 s 分组求和, 返回每组的平方和
    //   输入: returns[0..n-1], 分组大小 s
    //   输出: sum_{buckets} (sum(bucket))^2
    //   最后一个不满 s 的桶也会被计算 (R ceiling 行为)
    // =====================================================================
    static Real bucket_sum_sq(const std::vector<Real>& returns, Size s) {
        const Size n = returns.size();
        if (n == 0 || s == 0) return 0.0;
        Real total = 0.0;
        Size i = 0;
        while (i < n) {
            const Size end = std::min(i + s, n);
            Real bucket = 0.0;
            for (Size j = i; j < end; ++j) bucket += returns[j];
            total += bucket * bucket;
            i = end;
        }
        return total;
    }

    // =====================================================================
    // bucket_sum_cross: 按 s 分组求和, 返回两组对应桶的乘积和
    //   sum_{buckets} (sum(bucket1)) * (sum(bucket2))
    // =====================================================================
    static Real bucket_sum_cross(const std::vector<Real>& r1,
                                  const std::vector<Real>& r2, Size s) {
        const Size n = std::min(r1.size(), r2.size());
        if (n == 0 || s == 0) return 0.0;
        Real total = 0.0;
        Size i = 0;
        while (i < n) {
            const Size end = std::min(i + s, n);
            Real b1 = 0.0, b2 = 0.0;
            for (Size j = i; j < end; ++j) {
                b1 += r1[j];
                b2 += r2[j];
            }
            total += b1 * b2;
            i = end;
        }
        return total;
    }

    // =====================================================================
    // ravg_univariate: 单资产 subsampled RV (含 R 的 (m+1)/m 校正)
    //
    // R 源码逻辑 (L884-893):
    //   rvavg = bucket_sum_sq(rData, s)  // 主项 j=0, 全部 N 个观测
    //   for j = 1..s-1:
    //     rdatasub = rData[j .. N-s+j-1]  // 去掉前 j 和后 s-j 个, 剩 N-s 个
    //     nrow = N - s, frac = s, nrow/frac = (N-s)/s
    //     校正 = (nrow/frac + 1) / nrow/frac = (N/s) / (N/s - 1)
    //     rvavg += bucket_sum_sq(rdatasub, s) * 校正
    //   return rvavg / s
    //
    // 边界: 当 N == s 时, 偏移项无观测 (N-s=0), 跳过校正
    // =====================================================================
    static Real ravg_univariate(const std::vector<Real>& prices, Size s) {
        const std::vector<Real> ret = make_returns(prices);
        const Size N = ret.size();
        if (N == 0) return 0.0;
        if (s == 1) {
            // s=1: 退化为标准 RV
            Real rv = 0.0;
            for (Real v : ret) rv += v * v;
            return rv;
        }

        // 主项: j=0, 全部 N 个观测按 s 分组
        Real rvavg = bucket_sum_sq(ret, s);

        // 偏移项: j=1..s-1
        if (N > s) {
            // nrow = N - s, frac = s
            // 校正 = (nrow/frac + 1) / nrow/frac = (N/s) / (N/s - 1)
            const Real nrow_over_frac = static_cast<Real>(N - s) /
                                         static_cast<Real>(s);
            const Real correction = (nrow_over_frac > 0.0)
                ? (nrow_over_frac + 1.0) / nrow_over_frac
                : 1.0;

            for (Size j = 1; j < s; ++j) {
                // rdatasub = ret[j .. N-s+j-1], 长度 N-s
                std::vector<Real> sub(N - s);
                for (Size i = 0; i < N - s; ++i) {
                    sub[i] = ret[j + i];
                }
                rvavg += bucket_sum_sq(sub, s) * correction;
            }
        }

        return rvavg / static_cast<Real>(s);
    }

    // =====================================================================
    // ravg_bivariate: 双资产 subsampled cross-product (无 (m+1)/m 校正)
    //
    // R 源码逻辑 (L929-942):
    //   covavg = bucket_sum_cross(r1, r2, s)  // 主项 j=0
    //   for kk = 1..s-1:
    //     returns_sub = returns[kk .. N-s+kk-1]  // 双资产同步去除
    //     covavg += bucket_sum_cross(sub1, sub2, s)  // 注意: 无校正因子
    //   rdatamatrix[i,j] = covavg / s
    // =====================================================================
    static Real ravg_bivariate(const std::vector<Real>& prices1,
                                const std::vector<Real>& prices2, Size s) {
        const std::vector<Real> r1 = make_returns(prices1);
        const std::vector<Real> r2 = make_returns(prices2);
        const Size N = std::min(r1.size(), r2.size());
        if (N == 0) return 0.0;
        if (s == 1) {
            // s=1: 退化为标准 cross-product
            Real c = 0.0;
            for (Size i = 0; i < N; ++i) c += r1[i] * r2[i];
            return c;
        }

        // 主项: j=0
        Real covavg = bucket_sum_cross(r1, r2, s);

        // 偏移项: j=1..s-1 (无校正)
        if (N > s) {
            for (Size j = 1; j < s; ++j) {
                std::vector<Real> sub1(N - s), sub2(N - s);
                for (Size i = 0; i < N - s; ++i) {
                    sub1[i] = r1[j + i];
                    sub2[i] = r2[j + i];
                }
                covavg += bucket_sum_cross(sub1, sub2, s);
            }
        }

        return covavg / static_cast<Real>(s);
    }
};

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
