// SOURCE: PHASE5_HFE_SPEC §3.3
//   [BN-S 2006] Barndorff-Nielsen & Shephard, J. Financial Econometrics 4(1), 1-30
//   [BKS 2022]  Boudt, Kleen, Sjoerup, JSS 104(8), 1-36, doi:10.18637/jss.v104.i08
// R对照: BNSjumpTest(rdata, IVestimator="BV", IQestimator="TP", alpha=0.975)
//
// 关键 (2026-08-02 R 1.0.3 源码实测, verify_bns_internals.R):
//   - theta = tt("BV") = pi^2/4 + pi - 3  (非 pi^2/4)
//   - theta - 2 = pi^2/4 + pi - 5 ≈ 0.608994
//   - rTPQuar 实际公式: N*(N/(N-2))*C*sum(|r_{i-2}*r_{i-1}*r_i|^(4/3))
//     其中 C = (gamma(0.5)/(2^(2/3)*gamma(7/6)))^3
//     (非文档中的 "3σ 截断 + 4 次幂", 实际是 Truncated Power Variance)
//   - p-value = 2*Phi(-|z|)  (双尾, 非单尾)
//   - 拒绝域: |z| > z_{1-alpha}  (双尾)
//   - critical.value = c(-z_{1-alpha}, +z_{1-alpha})
//
// 实测验证 (2026-08-02):
//   CASE2 (n=5):   R z=-0.18669, C++ verify=-0.18669 ✓
//   CASE3 (n=100): R z=0.69272,  C++ verify=0.69272 ✓
//   CASE5 baseline: z=0.69271791587237797, p=0.48848659354781576
//   CASE6 baseline: z=4.6674971744735663, p=3.0489094124391184e-06
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/hfecon/measures/realized_measures.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// BNS Jump Test (Barndorff-Nielsen & Shephard 2006)
// =============================================================================
//
// 检验统计量 (R highfrequency 1.0.3 BNSjumpTest 源码实测):
//   Z_BNS = sqrt(N) * (RV - IV) / sqrt((theta - 2) * IQ)
//
// 其中 (IVestimator="BV", IQestimator="TP"):
//   IV  = hatIV("BV") = rBPCov = (pi/2) * sum |r_{i-1}*r_i|     [跳跃稳健]
//   IQ  = hatIQ("TP") = rTPQuar                                 [跳跃稳健]
//   theta = tt("BV") = pi^2/4 + pi - 3                          [R 1.0.3 tt 函数]
//   theta - 2 = pi^2/4 + pi - 5 ≈ 0.608994
//
// rTPQuar 实际公式 (R highfrequency 1.0.3 rTPQuar 源码, 非文档):
//   q[i] = |r_{i-2}| * |r_{i-1}| * |r_i|,  i = 2..n-1 (0-indexed, 长度 n-2)
//   N = n  (N = nrow(q) + 2 = (n-2) + 2 = n)
//   C = (gamma(0.5) / (2^(2/3) * gamma(7/6)))^3
//   rTPQuar = N * (N/(N-2)) * C * sum(q^(4/3))
//
// 注意: rTPQuar 函数名暗示 "Truncated Power Quarticity", 但实际是
//       "Truncated Power Variance" (4/3 次幂, 非 4 次幂; 滚动乘积, 非截断)
//
// 假设检验 (R 1.0.3 BNSjumpTest, type="linear", logTransform=FALSE, max=FALSE):
//   H0: 无跳跃 (价格路径连续)
//   H1: 存在跳跃
//   Z_BNS ->d N(0,1) under H0
//   p_value = 2 * Phi(-|Z|)                      [双尾]
//   critical.value = c(-z_{1-alpha}, +z_{1-alpha})
//   拒绝域: |Z| > z_{1-alpha}                    [双尾]
// =============================================================================

enum class IVEstimator { RV, BPV };           // IVestimator 参数 (R: "RV"/"BV")
enum class IQVEstimator { RQ, TPQ };          // IQestimator 参数 (R: "RQ"/"TP")

struct BNSJumpTestResult {
    Real z_statistic;       // Z_{BNS} 检验统计量
    Real p_value;           // 双尾 p-value = 2 * Phi(-|Z|)
    Real critical_value;    // z_{1-alpha} 临界值 (正侧, 负侧为 -critical_value)
    bool reject_null;       // 是否拒绝 H_0 (|Z| > critical_value)
    Real jump_ratio;        // (RV - BPV) / RV, 跳跃贡献比
    Real rv;
    Real bpv;
    Real rq;
    Real tpq;               // rTPQuar (Truncated Power Variance, R 1.0.3 实际公式)
    Size n_obs;
};

class BNSJumpTest {
public:
    // BNS 跳跃检验
    // 输入: 日内对数收益率序列
    // 参数:
    //   iv_est  - IV 估计量 (BPV 对应 R "BV", 跳跃稳健; RV 对应 R "RV", 含跳跃)
    //   iqvest  - IQ 估计量 (TPQ 对应 R "TP", 跳跃稳健; RQ 对应 R "RQ", 含跳跃)
    //   alpha   - 分位点 (R 1.0.3 默认 0.975, 即 97.5% -> z=1.96)
    // 异常: n < 4 时抛 invalid_argument (需足够观测估计 BPV/rTPQuar)
    //       IVestimator::RV 时抛 invalid_argument (R 1.0.3 tt 函数未定义 "RV" 的 theta)
    static BNSJumpTestResult test(
        const std::vector<Real>& log_returns,
        IVEstimator iv_est = IVEstimator::BPV,
        IQVEstimator iqvest = IQVEstimator::TPQ,
        Real alpha = 0.975) {

        const Size n = log_returns.size();
        if (n < 4) {
            throw std::invalid_argument(
                "BNSJumpTest::test requires n >= 4 returns");
        }
        if (iv_est == IVEstimator::RV) {
            // R 1.0.3 tt("RV") 未定义 (switch 无 "RV" case), BNSjumpTest 会出错
            throw std::invalid_argument(
                "BNSJumpTest::test: IVestimator::RV not supported "
                "(R 1.0.3 tt() does not define theta for 'RV')");
        }

        // 计算 realized measures (复用 RealizedMeasuresCalculator)
        const RealizedMeasures m = RealizedMeasuresCalculator::compute(log_returns);
        const Real rv  = m.rv;
        const Real bpv = m.bpv;
        const Real rq  = m.rq;

        // IV 估计量: BPV (R "BV"), 跳跃稳健
        const Real iv = bpv;

        // theta = tt("BV") = pi^2/4 + pi - 3  [R 1.0.3 tt 函数源码实测]
        // theta - 2 = pi^2/4 + pi - 5
        constexpr Real PI       = 3.14159265358979323846;
        constexpr Real THETA_BV = (PI * PI) / 4.0 + PI - 3.0;
        constexpr Real THETA_BV_MINUS_2 = THETA_BV - 2.0;  // ≈ 0.608994

        // IQ 估计量 + vartheta_BNS
        // R 公式: vartheta = (theta - 2) * IQ
        Real tpq = 0.0;
        Real iq;
        if (iqvest == IQVEstimator::TPQ) {
            // rTPQuar (R highfrequency 1.0.3 实际公式, 非 "Truncated Quarticity")
            // q[i] = |r_{i-2}|*|r_{i-1}|*|r_i|, i=2..n-1 (0-indexed)
            // rTPQuar = N * (N/(N-2)) * C * sum(q^(4/3))
            //   N = n, C = (gamma(0.5)/(2^(2/3)*gamma(7/6)))^3
            tpq = compute_rtpquar(log_returns, n);
            iq = tpq;
        } else {
            // RQ 模式: IQ = rQuar
            iq = rq;
        }
        const Real vartheta = THETA_BV_MINUS_2 * iq;

        // Z 统计量: Z = sqrt(N) * (RV - IV) / sqrt(vartheta)
        const Real numerator = rv - iv;
        Real z = 0.0;
        if (vartheta > 0.0) {
            z = std::sqrt(static_cast<Real>(n)) * numerator / std::sqrt(vartheta);
        }

        // 临界值: z_{1-alpha} (R qnorm(c(1-alpha, alpha)))
        const Real z_crit = inverse_normal_cdf(alpha);

        // p-value: 双尾 = 2 * Phi(-|Z|)  [R 1.0.3 BNSjumpTest 源码]
        const Real p_value = 2.0 * normal_cdf(-std::fabs(z));

        // 拒绝域: |Z| > z_{1-alpha}  (双尾)
        const bool reject = (std::fabs(z) > z_crit);

        // 跳跃贡献比
        Real jump_ratio = 0.0;
        if (rv > 0.0) {
            jump_ratio = (rv - bpv) / rv;
        }

        BNSJumpTestResult result;
        result.z_statistic    = z;
        result.p_value        = p_value;
        result.critical_value = z_crit;
        result.reject_null    = reject;
        result.jump_ratio     = jump_ratio;
        result.rv             = rv;
        result.bpv            = bpv;
        result.rq             = rq;
        result.tpq            = tpq;
        result.n_obs          = n;
        return result;
    }

    // 便捷接口: 输入价格序列
    static BNSJumpTestResult test_from_prices(
        const std::vector<Real>& prices,
        IVEstimator iv_est = IVEstimator::BPV,
        IQVEstimator iqvest = IQVEstimator::TPQ,
        Real alpha = 0.975) {
        return test(RealizedMeasuresCalculator::make_returns(prices),
                    iv_est, iqvest, alpha);
    }

private:
    // rTPQuar (R highfrequency 1.0.3 实际公式, verify_bns_internals.R 2026-08-02)
    // 注意: 函数名 "rTPQuar" 暗示 Quarticity, 实际是 Truncated Power Variance
    //
    // 算法:
    //   q = |r|  (绝对收益率)
    //   q_rolled[i] = q[i-2] * q[i-1] * q[i],  i = 2..n-1 (0-indexed, 长度 n-2)
    //   N = n  (N = nrow(q_rolled) + 2 = (n-2) + 2)
    //   C = (gamma(0.5) / (2^(2/3) * gamma(7/6)))^3
    //   rTPQuar = N * (N/(N-2)) * C * sum(q_rolled^(4/3))
    //
    // 异常: n < 3 时滚动乘积为空, 抛 invalid_argument (但 test() 已要求 n>=4)
    static Real compute_rtpquar(const std::vector<Real>& log_returns,
                                Size n) {
        if (n < 4) return 0.0;  // test() 已检查, 此处防御

        // 滚动 3 元素乘积 |r_{i-2}|*|r_{i-1}|*|r_i|, i=2..n-1
        Real sum_q_4_3 = 0.0;
        for (Size i = 2; i < n; ++i) {
            const Real q = std::fabs(log_returns[i - 2]) *
                           std::fabs(log_returns[i - 1]) *
                           std::fabs(log_returns[i]);
            sum_q_4_3 += std::pow(q, 4.0 / 3.0);
        }

        // 常数 C = (gamma(0.5)/(2^(2/3)*gamma(7/6)))^3
        // gamma(0.5) = sqrt(pi), gamma(7/6) 用 std::tgamma 计算
        const Real gamma_half = std::tgamma(0.5);      // = sqrt(pi)
        const Real gamma_7_6  = std::tgamma(7.0 / 6.0);
        const Real two_2_3    = std::pow(2.0, 2.0 / 3.0);
        const Real C = std::pow(gamma_half / (two_2_3 * gamma_7_6), 3.0);

        // N = n, N/(N-2) = n/(n-2)
        const Real N = static_cast<Real>(n);
        const Real N_over_N_minus_2 = N / static_cast<Real>(n - 2);

        return N * N_over_N_minus_2 * C * sum_q_4_3;
    }

    // 标准正态 CDF Phi(x) via std::erf
    // Phi(x) = 0.5 * (1 + erf(x / sqrt(2)))
    static Real normal_cdf(Real x) noexcept {
        return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
    }

    // 逆标准正态 CDF (Acklam's algorithm + Newton 精炼, 精度 ~1e-12)
    // 对于 alpha=0.975, 返回 1.959964; alpha=0.95 返回 1.644854
    // R qnorm(0.975) = 1.9599639845400536
    static Real inverse_normal_cdf(Real p) noexcept {
        // Acklam's coefficients
        static const Real a[] = {
            -3.969683028665376e+01,
             2.209460984245205e+02,
            -2.759285104469687e+02,
             1.383577518672690e+02,
            -3.066479806614716e+01,
             2.506628277459239e+00
        };
        static const Real b[] = {
            -5.447609879822406e+01,
             1.615858368580409e+02,
            -1.556989798598866e+02,
             6.680131188771972e+01,
            -1.328068155288572e+01
        };
        static const Real c[] = {
            -7.784894002430293e-03,
            -3.223964580411365e-01,
            -2.400758277161838e+00,
            -2.549732539343734e+00,
             4.374664141464968e+00,
             2.938163982698783e+00
        };
        static const Real d[] = {
             7.784695709041462e-03,
             3.224671290700398e-01,
             2.445134137142996e+00,
             3.754408661907416e+00
        };
        static const Real plow  = 0.02425;
        static const Real phigh = 1.0 - plow;

        Real x;
        if (p < plow) {
            // 左尾
            Real q = std::sqrt(-2.0 * std::log(p));
            x = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
                ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
        } else if (p <= phigh) {
            // 主体
            Real q = p - 0.5;
            Real r = q * q;
            x = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q /
                (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);
        } else {
            // 右尾
            Real q = std::sqrt(-2.0 * std::log(1.0 - p));
            x = -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
                 ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
        }

        // Newton 精炼一步 (提升到 ~1e-12)
        Real e = normal_cdf(x) - p;
        Real u = e * std::sqrt(2.0 * 3.14159265358979323846) * std::exp(x * x / 2.0);
        x = x - u / (1.0 + x * u / 2.0);
        return x;
    }
};

}  // namespace hfecon
}  // namespace v1
}  // namespace cpphub
