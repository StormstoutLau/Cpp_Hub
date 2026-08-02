// SOURCE: PHASE5_HFE_SPEC §3.2
//   [BN-S 2002] Barndorff-Nielsen & Shephard, J. Applied Econometrics 17, 453-475
//   [BN-S 2004] Barndorff-Nielsen & Shephard, J. Financial Econometrics 2(1), 1-37, doi:10.1093/jjfinec/nbh001
//   [BKS 2022]  Boudt, Kleen, Sjoerup, JSS 104(8), 1-36, doi:10.18637/jss.v104.i08
// R对照: rRVar / rRealizedVolatility / rQuar / rBPCov / rSVar / rCov
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// RealizedMeasures: RV / RVol / RQ / BPV / RSV+/- (BN-S 2002, 2004; BKS 2022 §3)
// =============================================================================
//
// 设日内对数收益率 r_1, ..., r_n (n = 观测数 - 1, 或直接输入收益率序列长度).
//
// Realized Variance:      RV  = sum_i r_i^2                           [BN-S 2002]
// Realized Volatility:    RVol = sqrt(RV)
// Realized Quarticity:    RQ  = ((n+1)/3) * sum_i r_i^4               [R highfrequency 1.0.3 rQuar]
//    注意: R highfrequency rQuar 源码实测 (verify_rq.R 2026-08-02):
//          N <- nrow(q) + 1; rQuar <- N/3 * colSums(q^4)
//          即 R 的 N = n+1 (而非 BN-S 2004 原始 (n/3) 系数)
//          case2 验证: (5+1)/3 * 1.15e-6 = 2.3e-6, 与 R 1.0.3 一致
// Bipower Variation:      BPV = (pi/2) * sum_{i=2}^n |r_{i-1}|*|r_i|  [R highfrequency rBPCov 1.0.3]
//    注意: BN-S 2004 原始定义含 n/(n-1) 系数, 但 R highfrequency rBPCov 实现省略 (渐近等价)
//    实测 (2026-08-02): case2 sum|r_{i-1}*r_i| = 0.0013, (pi/2)*0.0013 = 0.0020420352,
//                       与 R 1.0.3 baseline 一致 (若加 n/(n-1)=(5/4) 系数则得 0.002553, 不匹配)
// Realized Semivariance+: RSV+ = sum_i r_i^2 * 1{r_i > 0}
// Realized Semivariance-: RSV- = sum_i r_i^2 * 1{r_i <= 0}
//
// 性质: RV -> IV, BPV -> IV (跳跃稳健), RSV+ + RSV- = RV, RQ -> (2/3)*int sigma^4
//
// R 对照:
//   rRVar(rdata)                 -> rv
//   rRealizedVolatility(rdata)   -> rvol
//   rQuar(rdata)                 -> rq
//   rBPCov(rdata, makeReturns=F) -> bpv
//   rSVar(rdata)                 -> (rSVarupside=rsv_pos, rSVardownside=rsv_neg)
//   rCov(R, makeReturns=F)       -> realized_covariance(returns_matrix)
//
// 容差: 1e-12 (无噪声合成数据), 1e-10 (默认, R 对标硬约束, project_memory)
// =============================================================================

struct RealizedMeasures {
    Real rv;       // Realized Variance
    Real rvol;     // Realized Volatility = sqrt(RV)
    Real rq;       // Realized Quarticity
    Real bpv;      // Bipower Variation
    Real rsv_pos;  // Realized Semivariance (positive, r_i > 0)
    Real rsv_neg;  // Realized Semivariance (negative, r_i <= 0)
    Size n_obs;    // 观测数 (= 输入收益率序列长度)
};

class RealizedMeasuresCalculator {
public:
    // 输入: 日内对数收益率序列 (已通过 make_returns 计算)
    // 异常: n < 2 时 BPV 系数 n/(n-1) 未定义, 抛 invalid_argument
    static RealizedMeasures compute(const std::vector<Real>& log_returns) {
        const Size n = log_returns.size();
        if (n < 2) {
            throw std::invalid_argument(
                "RealizedMeasuresCalculator::compute requires n >= 2 returns "
                "(BPV coefficient n/(n-1) undefined for n < 2)");
        }

        Real sum_r2  = 0.0;
        Real sum_r4  = 0.0;
        Real sum_abs_prod = 0.0;  // BPV sum |r_{i-1}|*|r_i|, i=2..n
        Real rsv_pos = 0.0;
        Real rsv_neg = 0.0;

        for (Size i = 0; i < n; ++i) {
            const Real r = log_returns[i];
            const Real r2 = r * r;
            sum_r2 += r2;
            sum_r4 += r2 * r2;
            if (r > 0.0) {
                rsv_pos += r2;
            } else {
                rsv_neg += r2;  // r <= 0 (含 0, 与 R rSVar 一致: 1{r_i <= 0})
            }
            if (i >= 1) {
                sum_abs_prod += std::fabs(log_returns[i - 1]) * std::fabs(r);
            }
        }

        RealizedMeasures m;
        m.rv      = sum_r2;
        m.rvol    = std::sqrt(sum_r2);
        // RQ = ((n+1)/3) * sum(r^4)  [R highfrequency 1.0.3 rQuar 源码实测]
        // verify_rq.R (2026-08-02): N <- nrow(q) + 1; rQuar <- N/3 * colSums(q^4)
        // case2 验证: (5+1)/3 * 1.15e-6 = 2.3e-6, 与 R 1.0.3 一致
        m.rq      = (static_cast<Real>(n + 1) / 3.0) * sum_r4;
        // BPV = (pi/2) * sum |r_{i-1}|*|r_i|  [R highfrequency rBPCov 1.0.3 实测]
        // 注意: R 实现省略了 BN-S 2004 原始的 n/(n-1) 因子 (渐近等价)
        // 实测验证: case2 bpv = (pi/2)*0.0013 = 0.0020420352, 与 R 1.0.3 一致
        constexpr Real PI_HALF = 3.14159265358979323846 / 2.0;
        m.bpv     = PI_HALF * sum_abs_prod;
        m.rsv_pos = rsv_pos;
        m.rsv_neg = rsv_neg;
        m.n_obs   = n;
        return m;
    }

    // 便捷接口: 输入价格序列, 内部计算收益率
    // make_returns 语义: 返回等长序列, 首元素 = 0 (对齐 R makeReturns 数值向量行为)
    static RealizedMeasures compute_from_prices(const std::vector<Real>& prices) {
        return compute(make_returns(prices));
    }

    // 多资产 RV 协方差矩阵 (对应 R rCov)
    // 输入: returns_matrix[k] = 第 k 个资产的收益率序列 (各资产等长 n)
    // 输出: d x d 矩阵, 以 row-major std::vector<Real> 存储, 元素 (i,j) = sum_t r_i[t]*r_j[t]
    // 异常: 空输入或资产长度不一致时抛 invalid_argument
    static std::vector<Real> realized_covariance(
        const std::vector<std::vector<Real>>& returns_matrix) {
        const Size d = returns_matrix.size();
        if (d == 0) {
            throw std::invalid_argument(
                "realized_covariance: empty returns_matrix");
        }
        const Size n = returns_matrix[0].size();
        for (Size k = 0; k < d; ++k) {
            if (returns_matrix[k].size() != n) {
                throw std::invalid_argument(
                    "realized_covariance: all assets must have equal length");
            }
        }

        std::vector<Real> cov(d * d, 0.0);
        for (Size t = 0; t < n; ++t) {
            for (Size i = 0; i < d; ++i) {
                for (Size j = i; j < d; ++j) {
                    cov[i * d + j] += returns_matrix[i][t] * returns_matrix[j][t];
                }
            }
        }
        // 对称填充下三角
        for (Size i = 0; i < d; ++i) {
            for (Size j = 0; j < i; ++j) {
                cov[i * d + j] = cov[j * d + i];
            }
        }
        return cov;
    }

    // make_returns: 价格 -> 对数收益率
    // R makeReturns 数值向量行为: 返回等长序列, ret[0] = 0, ret[i] = log(p[i]/p[i-1])
    // 异常: prices 长度 < 1 时抛 invalid_argument
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
