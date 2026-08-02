// =============================================================================
// two_scale_cov.hpp
// Phase 5 v1.4.2 Wave B - Two-Scale 协方差估计 (Zhang, Mykland, Aït-Sahalia 2005)
//
// R 对照: rTSCov(pData, cor=FALSE, K=300, J=1, KCov=NULL, JCov=NULL,
//                 KVar=NULL, JVar=NULL, makePsd=FALSE)
//
// 文献: Zhang, Mykland, Aït-Sahalia (2005), JASA 100(472), 1394-1411,
//       doi:10.1198/016214505000000548
//
// SOURCE: PHASE5_HFE_SPEC §5.1, §5.3 D1
//   R highfrequency 1.0.3 rTSCov (v142_source_dump.txt L387-464)
//   R highfrequency 1.0.3 R/internalRealizedMeasures.R TSRV (L603-621)
//   R highfrequency 1.0.3 R/internalRealizedMeasures.R TSCov_bi (L567-600)
//
// 关键幻觉排除 (spec §5.3 D1):
//   对角线用 TSRV 公式 (adj = 1/(1 - nbarK/nbarJ))
//   非对角线用 TSCov_bi 公式 (adj = n/((K-J)*nbarK))
//   两个 adj 系数不同!
//
// 算法:
//   单资产 (d=1): TSRV(prices, K, J)
//     nbarK = (n-K+1)/K, nbarJ = (n-J+1)/J
//     adj = 1/(1 - nbarK/nbarJ)
//     lr_K = sum over k=1..K of sum(diff(logprices[seq(k,n,K)])^2)
//     lr_J = sum over j=1..J of sum(diff(logprices[seq(j,n,J)])^2)
//     TSRV = adj * ((1/K)*lr_K - (nbarK/nbarJ)*(1/J)*lr_J)
//
//   多资产 (d>1):
//     对角线: TSRV(prices_i, K, J)
//     非对角线: TSCov_bi(prices_i, prices_j, K, J)
//       adj = n/((K-J)*nbarK)
//       lr_K = sum over k=1..K of sum(diff(logp1[seq])*diff(logp2[seq]))
//       lr_J = sum over j=1..J of sum(diff(logp1[seq])*diff(logp2[seq]))
//       TSCov = adj * ((1/K)*lr_K - (nbarK/nbarJ)*(1/J)*lr_J)
// =============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/measures/cov_utils.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

struct TwoScaleCovResult {
    std::vector<Real> cov;
    Size n_assets;
    Size n_obs;
};

class TwoScaleCov {
public:
    // R rTSCov 等价
    // 输入: price_matrix[k] = 资产 k 的价格序列 (等长 n)
    // K: 慢尺度 (subsample 间隔), J: 快尺度 (通常 J=1)
    // make_psd: PSD 投影 (R makePsd, rTSCov 默认 FALSE)
    // cor: 相关系数矩阵
    // 异常: K<=J 或 n < 10*K 时抛 invalid_argument
    static TwoScaleCovResult estimate(
        const std::vector<std::vector<Real>>& price_matrix,
        Size K = 300,
        Size J = 1,
        bool make_psd_flag = false,
        bool cor = false) {

        const Size d = price_matrix.size();
        if (d == 0) {
            throw std::invalid_argument("TwoScaleCov: empty price_matrix");
        }
        const Size n = price_matrix[0].size();
        for (Size k = 0; k < d; ++k) {
            if (price_matrix[k].size() != n) {
                throw std::invalid_argument(
                    "TwoScaleCov: all assets must have equal length");
            }
        }
        if (n < 10 * K) {
            throw std::invalid_argument(
                "TwoScaleCov: need at least 10*K observations");
        }
        if (K <= J) {
            throw std::invalid_argument(
                "TwoScaleCov: K must be greater than J");
        }

        std::vector<Real> cov(d * d, 0.0);

        // 对角线: TSRV
        for (Size i = 0; i < d; ++i) {
            cov[i * d + i] = tsrv(price_matrix[i], K, J);
        }

        // 非对角线: TSCov_bi
        if (d > 1) {
            for (Size i = 1; i < d; ++i) {
                for (Size j = 0; j < i; ++j) {
                    const Real val = tscov_bi(
                        price_matrix[i], price_matrix[j], K, J);
                    cov[i * d + j] = val;
                    cov[j * d + i] = val;
                }
            }
        }

        // makePsd
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

        TwoScaleCovResult result;
        result.cov = cov;
        result.n_assets = d;
        result.n_obs = n;
        return result;
    }

private:
    // TSRV: Two-Scale Realized Volatility (单资产方差)
    // R highfrequency 1.0.3 R/internalRealizedMeasures.R TSRV (L603-621)
    static Real tsrv(const std::vector<Real>& prices, Size K, Size J) {
        const Size n = prices.size();
        std::vector<Real> logp(n);
        for (Size i = 0; i < n; ++i) logp[i] = std::log(prices[i]);

        const Real nbarK = static_cast<Real>(n - K + 1) / static_cast<Real>(K);
        const Real nbarJ = static_cast<Real>(n - J + 1) / static_cast<Real>(J);
        const Real adj = 1.0 / (1.0 - nbarK / nbarJ);

        // K-subsampling: sum of squared returns
        Real lr_K = 0.0;
        for (Size k = 1; k <= K; ++k) {
            for (Size i = k - 1; i + K < n; i += K) {
                const Real ret = logp[i + K] - logp[i];
                lr_K += ret * ret;
            }
        }

        // J-subsampling: sum of squared returns
        Real lr_J = 0.0;
        for (Size j = 1; j <= J; ++j) {
            for (Size i = j - 1; i + J < n; i += J) {
                const Real ret = logp[i + J] - logp[i];
                lr_J += ret * ret;
            }
        }

        return adj * ((1.0 / K) * lr_K - (nbarK / nbarJ) * (1.0 / J) * lr_J);
    }

    // TSCov_bi: Two-Scale Covariance (双资产交叉协方差)
    // R highfrequency 1.0.3 R/internalRealizedMeasures.R TSCov_bi (L567-600)
    // 注意: adj 与 TSRV 不同! (spec §5.3 D1)
    static Real tscov_bi(const std::vector<Real>& prices1,
                         const std::vector<Real>& prices2,
                         Size K, Size J) {
        const Size n = prices1.size();
        std::vector<Real> logp1(n), logp2(n);
        for (Size i = 0; i < n; ++i) {
            logp1[i] = std::log(prices1[i]);
            logp2[i] = std::log(prices2[i]);
        }

        const Real nbarK = static_cast<Real>(n - K + 1) / static_cast<Real>(K);
        const Real nbarJ = static_cast<Real>(n - J + 1) / static_cast<Real>(J);
        // TSCov_bi 的 adj 与 TSRV 不同!
        const Real adj = static_cast<Real>(n) /
                         (static_cast<Real>(K - J) * nbarK);

        // K-subsampling: cross-product
        Real lr_K = 0.0;
        for (Size k = 1; k <= K; ++k) {
            for (Size i = k - 1; i + K < n; i += K) {
                lr_K += (logp1[i + K] - logp1[i]) * (logp2[i + K] - logp2[i]);
            }
        }

        // J-subsampling: cross-product
        Real lr_J = 0.0;
        for (Size j = 1; j <= J; ++j) {
            for (Size i = j - 1; i + J < n; i += J) {
                lr_J += (logp1[i + J] - logp1[i]) * (logp2[i + J] - logp2[i]);
            }
        }

        return adj * ((1.0 / K) * lr_K - (nbarK / nbarJ) * (1.0 / J) * lr_J);
    }
};

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
