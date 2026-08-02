// =============================================================================
// hayashi_yoshida_cov.hpp
// Phase 5 v1.4.2 - Hayashi-Yoshida 非同步协方差估计
//
// R 对照: rHYCov(rData, cor=FALSE, period=1, alignBy="seconds",
//                 alignPeriod=1, makeReturns=FALSE, makePsd=TRUE)
//
// 文献: Hayashi & Yoshida (2005), J. Financial Econometrics 3(4),
//       doi:10.1093/jjfinec/nbi013
//
// SOURCE: PHASE5_HFE_SPEC §5.1, §5.3 D4/D5
//   R highfrequency 1.0.3 rHYCov (v142_source_dump.txt L895-1035)
//   R highfrequency 1.0.3 src/realizedMeasures.cpp pcovcc (L130-157)
//
// 关键幻觉排除 (spec §5.3):
//   D4: rHYCov 用整数索引而非时间戳 (pcovcc 退化版), period 仅用于聚合
//   D5: 默认 makePsd=TRUE (其他 Cov 函数 FALSE), 强制 PSD 投影
//
// 算法:
//   - 单资产 (d=1): 返回 rCov(aggrdata) = RV = sum(r^2)
//   - 多资产 (d>1):
//     对角线: diagonal[i] = rCov(aggrdata[,i]) = sum(r_i^2)
//     非对角线: cov[i,j] = sum(pcovcc(r_i, r_j, ...))
//     pcovcc 使用整数时间索引 (1..n) 进行周期聚合和时间匹配
//   - 默认 makePsd=TRUE: 对称矩阵 PSD 投影
//   - cor=TRUE: 转换为相关矩阵 (在 makePsd 之后, 对相关矩阵再做 makePsd)
//
// pcovcc 算法 (R highfrequency 1.0.3 src/realizedMeasures.cpp):
//   输入: a (资产i收益率), b (资产j收益率), at/bt (1-based 时间索引)
//         na, nap(=na/period), nb, period
//   1. 周期聚合: ap[i/period] += a[i], atp[i/period] = at[i]
//   2. 时间匹配: 对每个聚合区间 i:
//      tmpRet = sum(b[j] for j from prevj while bt[j] <= atp[i])
//      ans[i] = ap[i] * tmpRet
//   3. HY 协方差 = sum(ans)
//
// 默认参数 (period=1): 退化为标准同步协方差 sum(r_i[t]*r_j[t])
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

struct HayashiYoshidaResult {
    std::vector<Real> cov;   // d x d row-major 协方差矩阵
    Size n_assets;
    Size n_obs;
};

class HayashiYoshidaCov {
public:
    // R rHYCov 等价 (整数索引版, 对标 pcovcc)
    // 输入: returns_matrix[k] = 资产 k 的收益率序列 (等长 n)
    // period: 聚合周期 (R period * alignPeriod)
    // make_psd_flag: 是否投影到 PSD (R makePsd, rHYCov 默认 TRUE)
    // cor: 是否返回相关系数矩阵
    // 异常: 空输入或资产长度不一致时抛 invalid_argument
    static HayashiYoshidaResult estimate(
        const std::vector<std::vector<Real>>& returns_matrix,
        Size period = 1,
        bool make_psd_flag = true,
        bool cor = false) {

        const Size d = returns_matrix.size();
        if (d == 0) {
            throw std::invalid_argument(
                "HayashiYoshidaCov::estimate: empty returns_matrix");
        }
        const Size n = returns_matrix[0].size();
        for (Size k = 0; k < d; ++k) {
            if (returns_matrix[k].size() != n) {
                throw std::invalid_argument(
                    "HayashiYoshidaCov::estimate: all assets must have equal length");
            }
        }
        if (n == 0) {
            throw std::invalid_argument(
                "HayashiYoshidaCov::estimate: empty returns");
        }
        if (period == 0) {
            throw std::invalid_argument(
                "HayashiYoshidaCov::estimate: period must be >= 1");
        }

        std::vector<Real> cov(d * d, 0.0);

        // 对角线: RV of each asset (R rCov)
        for (Size i = 0; i < d; ++i) {
            Real rv = 0.0;
            for (Size t = 0; t < n; ++t) {
                rv += returns_matrix[i][t] * returns_matrix[i][t];
            }
            cov[i * d + i] = rv;
        }

        // 非对角线: pcovcc (Hayashi-Yoshida)
        if (d > 1) {
            const int na  = static_cast<int>(n);
            const int nb  = static_cast<int>(n);
            const int per = static_cast<int>(period);
            const int nap = na / per;

            for (Size i = 1; i < d; ++i) {
                for (Size j = 0; j < i; ++j) {
                    const Real hy = pcovcc_sum(
                        returns_matrix[i], returns_matrix[j],
                        na, nap, nb, per);
                    cov[i * d + j] = hy;
                    cov[j * d + i] = hy;
                }
            }
        }

        // makePsd (仅在 cor=FALSE 时应用于协方差矩阵)
        if (!cor && make_psd_flag && d > 1) {
            cov = make_psd(cov, d);
        }

        // 相关矩阵转换
        if (cor) {
            // 保存对角线平方根 (避免原地修改导致除数错误)
            std::vector<Real> diag_sqrt(d);
            for (Size i = 0; i < d; ++i) {
                diag_sqrt[i] = std::sqrt(cov[i * d + i]);
            }
            // rcor = D^{-1/2} * cov * D^{-1/2}
            for (Size i = 0; i < d; ++i) {
                for (Size j = 0; j < d; ++j) {
                    if (diag_sqrt[i] > 1e-300 && diag_sqrt[j] > 1e-300) {
                        cov[i * d + j] /= (diag_sqrt[i] * diag_sqrt[j]);
                    } else {
                        cov[i * d + j] = 0.0;
                    }
                }
            }
            // 对角线强制为 1
            for (Size i = 0; i < d; ++i) {
                cov[i * d + i] = 1.0;
            }
            // 对相关矩阵做 makePsd
            if (make_psd_flag && d > 1) {
                cov = make_psd(cov, d);
                for (Size i = 0; i < d; ++i) {
                    cov[i * d + i] = 1.0;
                }
            }
        }

        HayashiYoshidaResult result;
        result.cov = cov;
        result.n_assets = d;
        result.n_obs = n;
        return result;
    }

private:
    // pcovcc: Hayashi-Yoshida 交叉协方差 (周期聚合 + 时间匹配)
    // 对标 R highfrequency 1.0.3 src/realizedMeasures.cpp pcovcc()
    //
    // 输入: a (资产i收益率), b (资产j收益率)
    //       na (a 长度), nap (聚合区间数 = na/period), nb (b 长度), period (聚合周期)
    // 时间索引: 1-based (匹配 R c(1:n))
    // 返回: sum(ans) = HY 协方差
    static Real pcovcc_sum(
        const std::vector<Real>& a,
        const std::vector<Real>& b,
        int na, int nap, int nb, int period) {

        // 周期聚合缓冲区
        std::vector<Real> ap(nap, 0.0);
        std::vector<Real> atp(nap, 0.0);

        // 步骤 1: 周期聚合
        // ap[i/period] += a[i], atp[i/period] = at[i] = i+1 (1-based)
        for (int i = 0; i < na; ++i) {
            ap[i / period]  += a[i];
            atp[i / period]  = static_cast<Real>(i + 1);
        }

        // 步骤 2: 时间匹配 + 乘积求和
        Real total = 0.0;
        int prevj = 0;
        for (int i = 0; i < nap; ++i) {
            Real tmp_ret = 0.0;
            for (int j = prevj; j < nb; ++j) {
                tmp_ret += b[j];
                const Real bt_j = static_cast<Real>(j + 1);  // 1-based
                if (bt_j > atp[i]) {
                    prevj = j;
                    break;
                } else if (bt_j == atp[i]) {
                    prevj = j + 1;
                    break;
                }
            }
            total += ap[i] * tmp_ret;
        }
        return total;
    }
};

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
