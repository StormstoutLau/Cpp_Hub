// =============================================================================
// modulated_realized_cov.hpp
// Phase 5 v1.4.2 Wave B - Modulated Realized Covariance (Christensen, Podolskij, Vetter 2013)
//
// R 对照: rMRCov(pData, pairwise=FALSE, makePsd=FALSE, theta=0.8,
//                 crossAssetNoiseCorrection=FALSE)
//
// 文献: Christensen, Podolskij, Vetter (2013), J. Econometrics 173(1),
//       doi:10.1016/j.jeconom.2012.08.016
//
// SOURCE: PHASE5_HFE_SPEC §5.1, §5.3 D3
//   R highfrequency 1.0.3 R/internalPreaveringEstimators.R
//     crv (L2-16): 单资产 pre-averaged RV + 噪声校正
//     hatreturn (L19-34): pre-averaged returns
//     gfunction (L39-44): g(x) = min(x, 1-x)
//     preavbi (L47-72): 双资产 pre-averaged Cov
//   R highfrequency 1.0.3 R/realizedMeasures.R rMRCov (L642-721)
//   R highfrequency 1.0.3 src/internals.cpp preAveragingReturnsInternal (L83-106)
//
// 关键幻觉排除 (spec §5.3 D3):
//   三分支公式不兼容:
//   - 单资产 (crv): 1/(sqrt(N)*theta*psi2kn) * sum(r1^2) - 噪声校正
//   - 双资产 (preavbi): N/(N-kn+2) * 1/(psi2*kn) * sum(r1*r2) [- 噪声校正]
//   - 矩阵版: N/(N-kn+2) * 1/(psi2*kn) * S - 噪声校正对角线
//   注意: psi2kn (有限样本) 与 psi2=1/12 (渐近) 在不同分支使用
//
// 算法 (R 源码实测):
//   kn = floor(theta * sqrt(N))
//   g(x) = min(x, 1-x)
//   psi1kn = kn * sum_{h=1..kn} (g(h/kn) - g((h-1)/kn))^2
//   psi2kn = (1/kn) * sum_{h=1..kn} g(h/kn)^2
//   psi2 = 1/12  (渐近值, 用于双资产/矩阵版)
//
//   Pre-averaged returns (hatreturn):
//     weights[h] = g(h/kn) = min(h/kn, 1-h/kn), h=1..kn-1
//     hatre[i] = sum_{h=1..kn-1} weights[h] * r[i+h-1], i=0..N-kn
//
//   单资产 (crv):
//     crv = 1/(sqrt(N)*theta*psi2kn) * sum(hatre^2)
//           - psi1kn * (1/N) / (2*theta^2*psi2kn) * sum(r^2)
//
//   双资产 (preavbi, crossAssetNoiseCorrection=FALSE):
//     mrc = N/(N-kn+2) * 1/(psi2*kn) * sum(r1*r2)
//   双资产 (preavbi, crossAssetNoiseCorrection=TRUE):
//     mrc = N/(N-kn+2) * 1/(psi2*kn) * sum(r1*r2)
//           - psi1kn * (1/N) / (2*theta^2*psi2kn) * 1/(2*N) * sum(ret1*ret2)
//
//   矩阵版 (pairwise=FALSE):
//     refresh_time_matching -> preavreturn (N-kn+1 x d)
//     S = t(preavreturn) %*% preavreturn  (d x d)
//     对角线: - psi1kn*(1/N)/(2*theta^2*psi2kn) * 1/(2*N) * sum(r_i^2)
//     非对角线: 不校正 (crossAssetNoiseCorrection=FALSE)
//              或 - psi1kn*(1/N)/(2*theta^2*psi2kn) * 1/(2*N) * sum(r_i*r_j)
//     mrc = N/(N-kn+2) * 1/(psi2*kn) * S + 噪声校正项
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

struct ModulatedRealizedResult {
    std::vector<Real> cov;       // d x d 协方差矩阵 (row-major)
    Size n_assets;               // 资产数 d
    Size n_obs;                  // 观测数 N (refresh time 对齐后)
    Size kn;                     // pre-averaging 窗口
    Real theta;                  // theta 参数
};

class ModulatedRealizedCov {
public:
    // =====================================================================
    // 主接口: 多资产 rMRCov
    // 输入: prices_per_asset (d 个资产的价格序列, 等长)
    // 参数:
    //   theta                        - pre-averaging 参数 (默认 0.8, R 默认)
    //   pairwise                     - 是否使用 pairwise 模式 (默认 false)
    //   make_psd_flag                - 是否 PSD 投影 (默认 false, R 默认)
    //   cross_asset_noise_correction - 是否对非对角线做噪声校正 (默认 false)
    // 异常: 资产数为 0 或价格序列长度不一致抛 invalid_argument
    //       kn < 2 抛 invalid_argument (需足够观测计算 pre-averaged returns)
    // =====================================================================
    static ModulatedRealizedResult estimate(
        const std::vector<std::vector<Real>>& prices_per_asset,
        Real theta = 0.8,
        bool pairwise = false,
        bool make_psd_flag = false,
        bool cross_asset_noise_correction = false) {

        if (prices_per_asset.empty()) {
            throw std::invalid_argument(
                "ModulatedRealizedCov: at least one asset required");
        }
        const Size d = prices_per_asset.size();
        const Size N = prices_per_asset[0].size();
        for (Size a = 1; a < d; ++a) {
            if (prices_per_asset[a].size() != N) {
                throw std::invalid_argument(
                    "ModulatedRealizedCov: all assets must have equal length");
            }
        }

        // 单资产分支: crv
        if (d == 1) {
            const Real mrc = crv(prices_per_asset[0], theta);
            return ModulatedRealizedResult{
                {mrc}, 1, N, compute_kn(N, theta), theta
            };
        }

        // 多资产分支
        std::vector<Real> cov(d * d, 0.0);

        if (pairwise) {
            // Pairwise 模式: 对角线用 crv, 非对角线用 preavbi
            for (Size i = 0; i < d; ++i) {
                cov[i * d + i] = crv(prices_per_asset[i], theta);
            }
            for (Size i = 1; i < d; ++i) {
                for (Size j = 0; j < i; ++j) {
                    const Real c = preavbi(prices_per_asset[i],
                                           prices_per_asset[j],
                                           theta,
                                           cross_asset_noise_correction);
                    cov[i * d + j] = c;
                    cov[j * d + i] = c;
                }
            }
        } else {
            // 矩阵模式: refresh time 对齐 + pre-averaged returns + S
            // 注意: 对同步数据, refresh time = 原始时间, 无变化
            //       对非同步数据, 这里使用简化版 (假设已同步)
            const Size kn = compute_kn(N, theta);
            if (kn < 2) {
                throw std::invalid_argument(
                    "ModulatedRealizedCov: kn < 2, need more observations");
            }

            // 计算各资产的 pre-averaged returns
            // hatreturn_i: 长度 N-kn+1
            std::vector<std::vector<Real>> preav_returns(d);
            for (Size a = 0; a < d; ++a) {
                preav_returns[a] = hatreturn(prices_per_asset[a], kn);
            }
            const Size M = preav_returns[0].size();  // N - kn + 1

            // S = t(preavreturn) %*% preavreturn (d x d)
            for (Size i = 0; i < d; ++i) {
                for (Size j = i; j < d; ++j) {
                    Real s = 0.0;
                    for (Size t = 0; t < M; ++t) {
                        s += preav_returns[i][t] * preav_returns[j][t];
                    }
                    cov[i * d + j] = s;
                    cov[j * d + i] = s;
                }
            }

            // 噪声校正项
            const Real psi1kn = compute_psi1kn(kn);
            const Real psi2kn = compute_psi2kn(kn);
            const Real psi2 = 1.0 / 12.0;
            const Real noise_factor =
                psi1kn * (1.0 / static_cast<Real>(N)) /
                (2.0 * theta * theta * psi2kn) / (2.0 * static_cast<Real>(N));

            // 主系数: N/(N-kn+2) * 1/(psi2*kn)
            const Real main_factor = static_cast<Real>(N) /
                (static_cast<Real>(N) - static_cast<Real>(kn) + 2.0) /
                (psi2 * static_cast<Real>(kn));

            // 计算各资产的 RV (makeReturns 后的 sum(r^2))
            std::vector<Real> rv(d, 0.0);
            std::vector<std::vector<Real>> rets(d);
            for (Size a = 0; a < d; ++a) {
                rets[a] = make_returns(prices_per_asset[a]);
                for (Size t = 0; t < rets[a].size(); ++t) {
                    rv[a] += rets[a][t] * rets[a][t];
                }
            }

            for (Size i = 0; i < d; ++i) {
                for (Size j = i; j < d; ++j) {
                    Real mrc = main_factor * cov[i * d + j];
                    if (i == j) {
                        // 对角线总是做噪声校正
                        mrc -= noise_factor * rv[i];
                    } else if (cross_asset_noise_correction) {
                        // 非对角线条件做噪声校正
                        Real cross_rv = 0.0;
                        const Size n_ret = std::min(rets[i].size(), rets[j].size());
                        for (Size t = 0; t < n_ret; ++t) {
                            cross_rv += rets[i][t] * rets[j][t];
                        }
                        mrc -= noise_factor * cross_rv;
                    }
                    cov[i * d + j] = mrc;
                    cov[j * d + i] = mrc;
                }
            }
        }

        if (make_psd_flag) {
            cov = make_psd(cov, d);
        }

        return ModulatedRealizedResult{
            cov, d, N, compute_kn(N, theta), theta
        };
    }

    // =====================================================================
    // 单资产便捷接口
    // =====================================================================
    static Real estimate_univariate(const std::vector<Real>& prices,
                                     Real theta = 0.8) {
        return crv(prices, theta);
    }

private:
    // =====================================================================
    // gfunction: g(x) = min(x, 1-x), R internalPreaveringEstimators.R L39-44
    // =====================================================================
    static inline Real gfunction(Real x) {
        return std::min(x, 1.0 - x);
    }

    // =====================================================================
    // kn = floor(theta * sqrt(N))
    // =====================================================================
    static inline Size compute_kn(Size N, Real theta) {
        return static_cast<Size>(std::floor(theta * std::sqrt(static_cast<Real>(N))));
    }

    // =====================================================================
    // psi1kn = kn * sum_{h=1..kn} (g(h/kn) - g((h-1)/kn))^2
    // =====================================================================
    static Real compute_psi1kn(Size kn) {
        Real sum = 0.0;
        for (Size h = 1; h <= kn; ++h) {
            const Real g1 = gfunction(static_cast<Real>(h) / static_cast<Real>(kn));
            const Real g0 = gfunction(static_cast<Real>(h - 1) / static_cast<Real>(kn));
            const Real diff = g1 - g0;
            sum += diff * diff;
        }
        return static_cast<Real>(kn) * sum;
    }

    // =====================================================================
    // psi2kn = (1/kn) * sum_{h=1..kn} g(h/kn)^2
    // =====================================================================
    static Real compute_psi2kn(Size kn) {
        Real sum = 0.0;
        for (Size h = 1; h <= kn; ++h) {
            const Real g = gfunction(static_cast<Real>(h) / static_cast<Real>(kn));
            sum += g * g;
        }
        return sum / static_cast<Real>(kn);
    }

    // =====================================================================
    // make_returns: log price diff (R makeReturns)
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
    // hatreturn: pre-averaged returns (R hatreturn + preAveragingReturnsInternal)
    //
    // R preAveragingReturnsInternal (internals.cpp L83-106):
    //   weights = linspace(1, kn-1, kn-1) / kn  =  [1/kn, 2/kn, ..., (kn-1)/kn]
    //   foo = 1 - weights
    //   weights[weights > 1-weights] = foo[weights > 1-weights]  // = min(w, 1-w)
    //   out[i] = sum(ret[i..i+kn-2] * weights), i=0..N-kn
    //
    // 即: hatre[i] = sum_{h=1..kn-1} g(h/kn) * r[i+h-1], i=0..N-kn
    // =====================================================================
    static std::vector<Real> hatreturn(const std::vector<Real>& prices, Size kn) {
        const std::vector<Real> ret = make_returns(prices);
        const Size N = ret.size();  // 注意: ret 长度 = prices.size() - 1

        if (kn < 2) {
            // kn == 1 时, R 直接返回 rData
            return ret;
        }

        // 权重: w[h] = min(h/kn, 1-h/kn), h=1..kn-1
        std::vector<Real> weights(kn - 1);
        for (Size h = 1; h <= kn - 1; ++h) {
            const Real x = static_cast<Real>(h) / static_cast<Real>(kn);
            weights[h - 1] = gfunction(x);
        }

        // out[i] = sum_{h=0..kn-2} weights[h] * ret[i+h], i=0..N-kn
        // R: N = ret.n_rows + 1, out 长度 = N - kn + 1 = ret.size() - kn + 1
        const Size out_len = (N + 1 >= kn) ? (N + 1 - kn) : 0;
        std::vector<Real> out(out_len, 0.0);
        for (Size i = 0; i < out_len; ++i) {
            Real s = 0.0;
            for (Size h = 0; h < kn - 1; ++h) {
                s += weights[h] * ret[i + h];
            }
            out[i] = s;
        }
        return out;
    }

    // =====================================================================
    // crv: 单资产 pre-averaged RV + 噪声校正 (R internalPreaveringEstimators.R L2-16)
    //
    // crv = 1/(sqrt(N)*theta*psi2kn) * sum(r1^2)
    //       - psi1kn * (1/N) / (2*theta^2*psi2kn) * sum(rData^2)
    //
    // 其中 N = nrow(pData) = prices.size(), rData = makeReturns(prices)
    // =====================================================================
    static Real crv(const std::vector<Real>& prices, Real theta) {
        const Size N = prices.size();
        if (N < 4) {
            throw std::invalid_argument(
                "ModulatedRealizedCov::crv: need at least 4 prices");
        }
        const Size kn = compute_kn(N, theta);
        if (kn < 2) {
            throw std::invalid_argument(
                "ModulatedRealizedCov::crv: kn < 2, need more observations");
        }

        const Real psi1kn = compute_psi1kn(kn);
        const Real psi2kn = compute_psi2kn(kn);

        // r1 = pre-averaged returns
        const std::vector<Real> r1 = hatreturn(prices, kn);
        Real sum_r1_sq = 0.0;
        for (Real v : r1) sum_r1_sq += v * v;

        // rData = makeReturns(pData)
        const std::vector<Real> rData = make_returns(prices);
        Real sum_r_sq = 0.0;
        for (Real v : rData) sum_r_sq += v * v;

        const Real sqrt_N = std::sqrt(static_cast<Real>(N));
        const Real main_term = sum_r1_sq / (sqrt_N * theta * psi2kn);
        const Real noise_term = psi1kn * (1.0 / static_cast<Real>(N)) /
                                 (2.0 * theta * theta * psi2kn) * sum_r_sq;

        return main_term - noise_term;
    }

    // =====================================================================
    // preavbi: 双资产 pre-averaged Cov (R internalPreaveringEstimators.R L47-72)
    //
    // refreshTime 对齐后:
    //   kn = floor(theta * sqrt(N))
    //   psi2 = 1/12
    //   r1 = hatreturn(newprice1, kn), r2 = hatreturn(newprice2, kn)
    //
    //   crossAssetNoiseCorrection = FALSE:
    //     mrc = N/(N-kn+2) * 1/(psi2*kn) * sum(r1*r2)
    //   TRUE:
    //     mrc = N/(N-kn+2) * 1/(psi2*kn) * sum(r1*r2)
    //           - psi1kn*(1/N)/(2*theta^2*psi2kn) * 1/(2*N) * sum(ret1*ret2)
    //
    // 注意: 对同步数据, refreshTime 不改变序列, 直接使用原 prices
    // =====================================================================
    static Real preavbi(const std::vector<Real>& prices1,
                        const std::vector<Real>& prices2,
                        Real theta,
                        bool cross_asset_noise_correction) {
        const Size N = prices1.size();
        if (prices2.size() != N) {
            throw std::invalid_argument(
                "ModulatedRealizedCov::preavbi: prices must have equal length");
        }
        const Size kn = compute_kn(N, theta);
        if (kn < 2) {
            throw std::invalid_argument(
                "ModulatedRealizedCov::preavbi: kn < 2");
        }

        const Real psi2 = 1.0 / 12.0;
        const std::vector<Real> r1 = hatreturn(prices1, kn);
        const std::vector<Real> r2 = hatreturn(prices2, kn);

        Real sum_r1r2 = 0.0;
        const Size M = std::min(r1.size(), r2.size());
        for (Size t = 0; t < M; ++t) {
            sum_r1r2 += r1[t] * r2[t];
        }

        const Real main_factor = static_cast<Real>(N) /
            (static_cast<Real>(N) - static_cast<Real>(kn) + 2.0) /
            (psi2 * static_cast<Real>(kn));
        Real mrc = main_factor * sum_r1r2;

        if (cross_asset_noise_correction) {
            const Real psi1kn = compute_psi1kn(kn);
            const Real psi2kn = compute_psi2kn(kn);
            const std::vector<Real> ret1 = make_returns(prices1);
            const std::vector<Real> ret2 = make_returns(prices2);
            Real sum_ret = 0.0;
            const Size n_ret = std::min(ret1.size(), ret2.size());
            for (Size t = 0; t < n_ret; ++t) {
                sum_ret += ret1[t] * ret2[t];
            }
            const Real noise_factor = psi1kn * (1.0 / static_cast<Real>(N)) /
                (2.0 * theta * theta * psi2kn) / (2.0 * static_cast<Real>(N));
            mrc -= noise_factor * sum_ret;
        }

        return mrc;
    }
};

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
