// =============================================================================
// rank_jump_test.hpp
// Phase 5 v1.4.3 - HFE Rank Jump Test (Bollerslev & Todorov 2011)
//
// 对标 R highfrequency 1.0.3:
//   rankJumpTest(marketPrice, stockPrices, alpha=c(5,3), coarseFreq=10,
//               localWindow=30, rank=1, BoxCox=1, quantiles, nBoot=1000, ...)
//   — jumpTests.R L976 + internalJumpTests.R L149 (jumpDetection) + L115 (BoxCox__)
//
// 排幻觉点 (R 源码实测 2026-08-03):
//   D17: jumpDetection — Un = alpha*sqrt(kronecker(pmin(bpv,rv), TODfit))*(1/nRets)^0.49
//        论文无 TOD 调整, R 特有日内模式修正
//   D18: rankJumpTest — jumps = sum(ret[jumpIdx+i]) for i=0..coarseFreq-1 (累积窗口)
//   D19: rankJumpTest — svd(jumps, nu=nrow, nv=ncol) 全 SVD (非 reduced)
//   D20: rankJumpTest — testStat = sum(BoxCox__(d^2, a)) (论文无 BoxCox)
//   D21: rankJumpTest — dxc = pmax(pmin(ret, Un), -Un) 截断, bootstrap 用截断收益
//   D22: BoxCox__ — lambda=0 → log(1+x) (标准 BoxCox 是 log(x), R 用 1+x 避免 log(0))
//   D23: timeOfDayAdjustments — 1.249531*rowMeans(|r_i*r_{i+1}*r_{i+2}|^(2/3))
//        1.249531 = (2^(2/3)*gamma(7/6)/gamma(1/2))^2 (tripower quarticity 归一化)
//
// 简化假设 (v1.4.3):
//   - 输入为已聚合的对数收益率序列 (跳过 aggregatePrice + makeReturns)
//   - 多日数据通过 nRets x nDays 扁平 vector 传入 (按列存储, 日为列)
//   - bootstrap 用 std::mt19937_64 固定种子, 不与 R runif 数值对标 (仅验证可复现性)
//
// 容差: 1e-8 (SVD Jacobi 收敛 + bootstrap 浮点累积)
//
// SOURCE:
//   [BT 2011] Bollerslev & Todorov, JFE 9(2), doi:10.1093/jjfinec/nbr010
// R 源码: tests/fixtures/hfe/hf_src/highfrequency/R/jumpTests.R L976
//         tests/fixtures/hfe/hf_src/highfrequency/R/internalJumpTests.R L115 (BoxCox__)
//         tests/fixtures/hfe/hf_src/highfrequency/R/internalJumpTests.R L125 (TOD)
//         tests/fixtures/hfe/hf_src/highfrequency/R/internalJumpTests.R L149 (jumpDetection)
// =============================================================================
#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace hfecon {

// =============================================================================
// detail namespace — 内部辅助函数
// =============================================================================
namespace detail {

// -----------------------------------------------------------------------------
// BoxCox 变换 (D22)
// R: internalJumpTests.R L115-121
//   BoxCox__(x, lambda):
//     lambda == 0 → log(1+x)        (标准 BoxCox 是 log(x), R 用 1+x 避免 log(0))
//     lambda != 0 → ((1+x)^lambda - 1)/lambda
// -----------------------------------------------------------------------------
inline Real box_cox(Real x, Real lambda) {
    if (lambda == 0.0) {
        return std::log(1.0 + x);
    }
    return (std::pow(1.0 + x, lambda) - 1.0) / lambda;
}

// -----------------------------------------------------------------------------
// TOD 调整结果 (D23)
// -----------------------------------------------------------------------------
struct TODAdjustment {
    std::vector<Real> timeOfDayScatter;  // nRets, 归一化 mean=1
    std::vector<Real> timeOfDayFit;      // nRets, 归一化 mean=1
    std::vector<Real> timeOfDayBeta;     // polyOrder+1
};

// -----------------------------------------------------------------------------
// time_of_day_adjustments — 日内模式调整 (D23)
// R: internalJumpTests.R L125-145
//   timePolyMatrix: Vandermonde (1:nRets)^{0:polyOrder}, nRets x (polyOrder+1)
//   timeOfDayScatter = 1.249531 * rowMeans(|r[,1:(m-2)]| * |r[,2:(m-1)]| * |r[,3:m]|^(2/3))
//   timeOfDayBeta = (X^T X)^{-1} X^T scatter (OLS)
//   timeOfDayFit = X %*% beta, 归一化 mean=1
//   timeOfDayScatter 归一化 mean=1 (R L140)
//
// 输入: returns_flat (nRets*nDays, 按列存储, 日为列), nRets, nDays, polyOrder
// 输出: TODAdjustment
//
// 边界: nDays < 3 时无法计算三元组, 返回全 1 的 TODfit (退化模式)
// -----------------------------------------------------------------------------
inline TODAdjustment time_of_day_adjustments(
        const std::vector<Real>& returns_flat,
        Size nRets, Size nDays, int polyOrder = 2) {

    TODAdjustment adj;
    const Size ncol = static_cast<Size>(polyOrder + 1);

    if (nDays < 3) {
        // R: returns[,1:(m-2)] 当 m<3 时为空矩阵, rowMeans=NaN
        // 退化: TODfit = 全 1, scatter = 全 1
        adj.timeOfDayScatter.assign(nRets, 1.0);
        adj.timeOfDayFit.assign(nRets, 1.0);
        adj.timeOfDayBeta.assign(ncol, 0.0);
        return adj;
    }

    // timePolyMatrix: nRets x ncol, 每行 [1, t, t^2, ...] (t = 1..nRets, R 1-based)
    std::vector<std::vector<Real>> X(nRets, std::vector<Real>(ncol));
    for (Size t = 0; t < nRets; ++t) {
        Real tt = static_cast<Real>(t + 1);
        Real pow_val = 1.0;
        for (Size p = 0; p < ncol; ++p) {
            X[t][p] = pow_val;
            pow_val *= tt;
        }
    }

    // timeOfDayScatter = 1.249531 * rowMeans(|r_d1 * r_d2 * r_d3|^(2/3))
    // n_triples = nDays - 2 (三元组数: (d, d+1, d+2) for d=0..nDays-3)
    const Size n_triples = nDays - 2;
    const Real const_1249531 = 1.249531;
    adj.timeOfDayScatter.resize(nRets);

    for (Size i = 0; i < nRets; ++i) {
        Real sum = 0.0;
        for (Size d = 0; d < n_triples; ++d) {
            Real r1 = returns_flat[d * nRets + i];
            Real r2 = returns_flat[(d + 1) * nRets + i];
            Real r3 = returns_flat[(d + 2) * nRets + i];
            Real product = std::fabs(r1) * std::fabs(r2) * std::fabs(r3);
            sum += std::pow(product, 2.0 / 3.0);
        }
        adj.timeOfDayScatter[i] = const_1249531 * sum / static_cast<Real>(n_triples);
    }

    // OLS: beta = (X^T X)^{-1} X^T scatter
    std::vector<std::vector<Real>> XtX(ncol, std::vector<Real>(ncol, 0.0));
    std::vector<Real> Xty(ncol, 0.0);
    for (Size i = 0; i < ncol; ++i) {
        for (Size j = 0; j < ncol; ++j) {
            Real s = 0.0;
            for (Size t = 0; t < nRets; ++t) {
                s += X[t][i] * X[t][j];
            }
            XtX[i][j] = s;
        }
        Real s = 0.0;
        for (Size t = 0; t < nRets; ++t) {
            s += X[t][i] * adj.timeOfDayScatter[t];
        }
        Xty[i] = s;
    }

    // Gauss-Jordan 求解 (partial pivoting)
    std::vector<std::vector<Real>> A = XtX;
    std::vector<Real> b = Xty;
    for (Size col = 0; col < ncol; ++col) {
        Size piv = col;
        Real maxv = std::fabs(A[col][col]);
        for (Size r = col + 1; r < ncol; ++r) {
            if (std::fabs(A[r][col]) > maxv) {
                maxv = std::fabs(A[r][col]);
                piv = r;
            }
        }
        if (maxv < 1e-15) {
            throw std::runtime_error(
                "time_of_day_adjustments: singular X^T X matrix");
        }
        if (piv != col) {
            std::swap(A[piv], A[col]);
            std::swap(b[piv], b[col]);
        }
        Real akk = A[col][col];
        for (Size r = col + 1; r < ncol; ++r) {
            Real f = A[r][col] / akk;
            if (f == 0.0) continue;
            for (Size c = col; c < ncol; ++c) A[r][c] -= f * A[col][c];
            b[r] -= f * b[col];
        }
    }
    std::vector<Real> beta(ncol);
    for (Size i = ncol; i-- > 0;) {
        Real s = b[i];
        for (Size j = i + 1; j < ncol; ++j) s -= A[i][j] * beta[j];
        beta[i] = s / A[i][i];
    }
    adj.timeOfDayBeta = beta;

    // timeOfDayFit = X %*% beta
    adj.timeOfDayFit.resize(nRets);
    for (Size t = 0; t < nRets; ++t) {
        Real fit = 0.0;
        for (Size p = 0; p < ncol; ++p) {
            fit += X[t][p] * beta[p];
        }
        adj.timeOfDayFit[t] = fit;
    }

    // 归一化 mean=1 (R L138)
    Real mean_fit = 0.0;
    for (Size t = 0; t < nRets; ++t) mean_fit += adj.timeOfDayFit[t];
    mean_fit /= static_cast<Real>(nRets);
    if (std::fabs(mean_fit) > 1e-15) {
        for (Size t = 0; t < nRets; ++t) adj.timeOfDayFit[t] /= mean_fit;
    }

    // 归一化 scatter mean=1 (R L140)
    Real mean_scatter = 0.0;
    for (Size t = 0; t < nRets; ++t) mean_scatter += adj.timeOfDayScatter[t];
    mean_scatter /= static_cast<Real>(nRets);
    if (std::fabs(mean_scatter) > 1e-15) {
        for (Size t = 0; t < nRets; ++t) adj.timeOfDayScatter[t] /= mean_scatter;
    }

    return adj;
}

// -----------------------------------------------------------------------------
// 跳跃检测结果 (D17)
// -----------------------------------------------------------------------------
struct JumpDetection {
    std::vector<int> jumpIndices;   // 0-based 全局索引
    std::vector<Real> Un;           // nRets*nDays, 阈值
    std::vector<Real> todFit;       // nRets, 日内模式 (归一化)
};

// -----------------------------------------------------------------------------
// jump_detection — 跳跃检测 (D17)
// R: internalJumpTests.R L149-161
//   returns matrix: nRets x nDays (按列, 日为列)
//   bpv = (pi/2) * colSums(|r[0:n-1]| * |r[1:n]|)  (每列一个 BPV)
//   rv = colSums(r^2)  (每列一个 RV)
//   TODadjustments(returns, n=nRets, m=nDays, polyOrder=2)
//   Un = alpha * sqrt(kronecker(pmin(bpv,rv), TODfit)) * (1/nRets)^0.49
//   jumpIndices = which(|r| > Un)
//
// 输入: returns_flat (nRets*nDays, 按列存储), alpha, nRets, nDays
// 输出: JumpDetection
// -----------------------------------------------------------------------------
inline JumpDetection jump_detection(
        const std::vector<Real>& returns_flat,
        Real alpha, Size nRets, Size nDays) {

    JumpDetection jd;
    const Real pi_half = std::acos(-1.0) / 2.0;
    const Size total = nRets * nDays;

    // bpv, rv: 每日一个
    std::vector<Real> bpv(nDays), rv(nDays);
    for (Size d = 0; d < nDays; ++d) {
        Real bpv_sum = 0.0, rv_sum = 0.0;
        for (Size i = 0; i + 1 < nRets; ++i) {
            Real r1 = returns_flat[d * nRets + i];
            Real r2 = returns_flat[d * nRets + i + 1];
            bpv_sum += std::fabs(r1) * std::fabs(r2);
            rv_sum += r1 * r1;
        }
        // rv 包含最后一个观测
        rv_sum += returns_flat[d * nRets + nRets - 1]
                * returns_flat[d * nRets + nRets - 1];
        bpv[d] = pi_half * bpv_sum;
        rv[d] = rv_sum;
    }

    // TOD 调整
    auto adj = time_of_day_adjustments(returns_flat, nRets, nDays, 2);
    jd.todFit = adj.timeOfDayFit;

    // Un = alpha * sqrt(kronecker(pmin(bpv,rv), TODfit)) * (1/nRets)^0.49
    // kronecker(pmin(bpv,rv), TODfit): 长度 = nDays * nRets
    // 顺序: (day1_tod1,...,day1_todN, day2_tod1,...,dayD_todN) — 与 returns_flat 相同
    jd.Un.resize(total);
    const Real scale = std::pow(1.0 / static_cast<Real>(nRets), 0.49);
    for (Size d = 0; d < nDays; ++d) {
        Real pmin = std::min(bpv[d], rv[d]);
        for (Size t = 0; t < nRets; ++t) {
            jd.Un[d * nRets + t] = alpha * std::sqrt(pmin * adj.timeOfDayFit[t]) * scale;
        }
    }

    // jumpIndices = which(|r| > Un)
    for (Size i = 0; i < total; ++i) {
        if (std::fabs(returns_flat[i]) > jd.Un[i]) {
            jd.jumpIndices.push_back(static_cast<int>(i));
        }
    }

    return jd;
}

// -----------------------------------------------------------------------------
// SVD 全分解结果 (D19)
// -----------------------------------------------------------------------------
struct SVDResult {
    std::vector<Real> U;    // m x m, row-major
    std::vector<Real> d;    // min(m, n) 个奇异值 (降序)
    std::vector<Real> V;    // n x n, row-major
    Size m = 0, n = 0;
};

// -----------------------------------------------------------------------------
// svd_full — 全 SVD 分解 (D19)
// R: svd(A, nu=nrow, nv=ncol) → U (m x m), d (min(m,n)), V (n x n)
// 实现: one-sided Jacobi SVD + Gram-Schmidt 补全
//
// 算法:
//   1. 如果 m < n, 对 A^T 做 SVD, 交换 U, V
//   2. one-sided Jacobi: 迭代正交化列, A = U_red * diag(d) * V^T
//   3. 归一化 U_red 列, d[i] = ||U_red[:,i]||
//   4. 按 d 降序排列 U, V, d
//   5. Gram-Schmidt 补全 U_red 到 m x m (如果 m > n)
// -----------------------------------------------------------------------------
inline SVDResult svd_full(const std::vector<Real>& A, Size m, Size n) {
    SVDResult result;
    result.m = m;
    result.n = n;

    if (m == 0 || n == 0) return result;

    // 转置处理 (one-sided Jacobi 要求 m >= n)
    bool transposed = (m < n);
    Size M, N;
    std::vector<Real> At;
    if (!transposed) {
        M = m; N = n;
        At = A;
    } else {
        M = n; N = m;
        At.resize(M * N);
        for (Size i = 0; i < m; ++i)
            for (Size j = 0; j < n; ++j)
                At[j * m + i] = A[i * n + j];
    }

    // one-sided Jacobi SVD: At = Up * diag(d) * Vp^T
    // Up: M x N (正交列), Vp: N x N (正交)
    std::vector<Real> Up = At;  // M x N, row-major (Up[i*N + j])
    std::vector<Real> Vp(N * N, 0.0);
    for (Size i = 0; i < N; ++i) Vp[i * N + i] = 1.0;

    const int max_iter = 100;
    const Real tol = 1e-14;

    for (int iter = 0; iter < max_iter; ++iter) {
        bool converged = true;
        for (Size i = 0; i < N; ++i) {
            for (Size j = i + 1; j < N; ++j) {
                // 计算列 i, j 的内积
                Real alpha = 0.0;
                for (Size k = 0; k < M; ++k) {
                    alpha += Up[k * N + i] * Up[k * N + j];
                }
                // 相对容差: |alpha| / sqrt(||ci||^2 * ||cj||^2)
                Real norm_i = 0.0, norm_j = 0.0;
                for (Size k = 0; k < M; ++k) {
                    norm_i += Up[k * N + i] * Up[k * N + i];
                    norm_j += Up[k * N + j] * Up[k * N + j];
                }
                Real threshold = tol * std::sqrt(norm_i * norm_j);
                if (threshold < 1e-300) threshold = 1e-300;

                if (std::fabs(alpha) > threshold) {
                    converged = false;
                    // Jacobi 旋转角
                    Real a = norm_i, c = norm_j;
                    Real theta;
                    if (std::fabs(a - c) < 1e-300) {
                        theta = std::acos(-1.0) / 4.0;
                    } else {
                        theta = 0.5 * std::atan2(2.0 * alpha, a - c);
                    }
                    Real cos_t = std::cos(theta), sin_t = std::sin(theta);
                    // 旋转 Up 列 i, j
                    for (Size k = 0; k < M; ++k) {
                        Real u_i = Up[k * N + i];
                        Real u_j = Up[k * N + j];
                        Up[k * N + i] = cos_t * u_i + sin_t * u_j;
                        Up[k * N + j] = -sin_t * u_i + cos_t * u_j;
                    }
                    // 旋转 Vp 列 i, j
                    for (Size k = 0; k < N; ++k) {
                        Real v_i = Vp[k * N + i];
                        Real v_j = Vp[k * N + j];
                        Vp[k * N + i] = cos_t * v_i + sin_t * v_j;
                        Vp[k * N + j] = -sin_t * v_i + cos_t * v_j;
                    }
                }
            }
        }
        if (converged) break;
    }

    // d[i] = ||Up[:,i]||, 归一化 Up[:,i]
    std::vector<Real> d(N);
    for (Size i = 0; i < N; ++i) {
        Real norm = 0.0;
        for (Size k = 0; k < M; ++k) {
            norm += Up[k * N + i] * Up[k * N + i];
        }
        norm = std::sqrt(norm);
        d[i] = norm;
        if (norm > 1e-300) {
            for (Size k = 0; k < M; ++k) {
                Up[k * N + i] /= norm;
            }
        }
    }

    // 按 d 降序排列
    std::vector<Size> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](Size a, Size b) { return d[a] > d[b]; });

    std::vector<Real> Up_sorted(M * N), Vp_sorted(N * N), d_sorted(N);
    for (Size i = 0; i < N; ++i) {
        d_sorted[i] = d[idx[i]];
        for (Size k = 0; k < M; ++k) {
            Up_sorted[k * N + i] = Up[k * N + idx[i]];
        }
        for (Size k = 0; k < N; ++k) {
            Vp_sorted[k * N + i] = Vp[k * N + idx[i]];
        }
    }

    // 补全 Up_sorted 到 M x M (Gram-Schmidt, 如果 M > N)
    std::vector<Real> Up_full(M * M, 0.0);
    for (Size i = 0; i < N; ++i) {
        for (Size k = 0; k < M; ++k) {
            Up_full[k * M + i] = Up_sorted[k * N + i];
        }
    }
    Size n_complete = N;
    for (Size k = 0; k < M && n_complete < M; ++k) {
        std::vector<Real> e(M, 0.0);
        e[k] = 1.0;
        // Modified Gram-Schmidt: 减去在已有列上的投影
        for (Size i = 0; i < n_complete; ++i) {
            Real proj = 0.0;
            for (Size j = 0; j < M; ++j) {
                proj += e[j] * Up_full[j * M + i];
            }
            for (Size j = 0; j < M; ++j) {
                e[j] -= proj * Up_full[j * M + i];
            }
        }
        Real norm = 0.0;
        for (Size j = 0; j < M; ++j) norm += e[j] * e[j];
        norm = std::sqrt(norm);
        if (norm > 1e-10) {
            for (Size j = 0; j < M; ++j) Up_full[j * M + n_complete] = e[j] / norm;
            ++n_complete;
        }
    }

    // 设置返回值 (处理转置)
    if (!transposed) {
        // A = At = Up_full * diag(d_ext) * Vp_sorted^T
        // U = Up_full (m x m), V = Vp_sorted (n x n), d = d_sorted (n 个)
        result.U = Up_full;
        result.V = Vp_sorted;
        result.d = d_sorted;
    } else {
        // A = At^T = (Up_full * diag(d_ext) * Vp_sorted^T)^T
        //         = Vp_sorted * diag(d_ext) * Up_full^T
        // U_A = Vp_sorted (m x m = m x m, 因为 Vp_sorted 是 N x N = m x m)
        // V_A = Up_full (n x n = n x n, 因为 Up_full 是 M x M = n x n)
        // d_A = d_sorted (m 个, min(m,n) = m)
        result.U = Vp_sorted;
        result.V = Up_full;
        result.d = d_sorted;
    }

    return result;
}

} // namespace detail

// =============================================================================
// Rank Jump Test Result
// =============================================================================
struct RankJumpTestResult {
    std::vector<Real> criticalValues;   // 对应 quantiles 的临界值
    std::vector<Real> testStatistic;    // 对应 boxCox 参数的检验统计量
    std::vector<int> jumpIndices;       // 0-based 全局跳跃索引
    Size nAssets = 0;
    Size nJumps = 0;
};

// =============================================================================
// rank_jump_test — Rank 跳跃检验 (Bollerslev & Todorov 2011)
//
// 算法 (R 源码 jumpTests.R L976-1130 实测):
//   1. 数据准备: 输入为已聚合收益率 (marketReturns + stockReturns)
//   2. jumpDetection (D17): 市场跳跃检测, Un 阈值
//   3. 个股 jumpDetection: 每个股票的 Un (用于 bootstrap 截断)
//   4. jumps 累积 (D18): sum(stockReturns[jumpIdx+i]) for i=0..coarseFreq-1
//   5. SVD 全分解 (D19): U2=U[:,rank:], V2=V[:,rank:], singularValues=d[rank:]^2
//   6. testStatistic (D20): sum(BoxCox__(singularValues, a)) per boxCox
//   7. bootstrap (D21): nBoot 次, dxc 截断 + zetaStar + tmp = t(U2)%*%zetaStar%*%V2
//   8. criticalValues = quantile(simTestStat, quantiles) (type 7)
//
// 输入:
//   marketReturns: 长度 nRets*nDays, 按列存储 (日为列)
//   stockReturns: nAssets 个资产, 每个长度 nRets*nDays
//   nRets: 每日观测数, nDays: 天数
//   alphaMarket/alphaStock: jumpDetection 阈值 (标准差倍数)
//   coarseFreq: 跳跃累积窗口长度
//   localWindow: bootstrap 局部窗口
//   rank: 截断秩
//   boxCox: BoxCox 变换参数列表
//   quantiles: 临界值分位数列表
//   nBoot: bootstrap 次数
//   dontTestAtBoundaries: 是否避免跨日采样
//   seed: 随机种子 (可复现)
// =============================================================================
inline RankJumpTestResult rank_jump_test(
        const std::vector<Real>& marketReturns,
        const std::vector<std::vector<Real>>& stockReturns,
        Size nRets, Size nDays,
        Real alphaMarket = 5.0,
        Real alphaStock = 3.0,
        int coarseFreq = 10,
        int localWindow = 30,
        int rank = 1,
        const std::vector<Real>& boxCox = {1.0},
        const std::vector<Real>& quantiles = {0.9, 0.95, 0.99},
        int nBoot = 1000,
        bool dontTestAtBoundaries = true,
        std::uint64_t seed = 12345ULL) {

    RankJumpTestResult result;

    // --- 输入校验 ---
    const Size total = nRets * nDays;
    if (marketReturns.empty() || stockReturns.empty()) {
        throw std::invalid_argument(
            "rank_jump_test: marketReturns and stockReturns must be non-empty");
    }
    if (marketReturns.size() != total) {
        throw std::invalid_argument(
            "rank_jump_test: marketReturns size != nRets*nDays");
    }
    const Size nAssets = stockReturns.size();
    for (Size a = 0; a < nAssets; ++a) {
        if (stockReturns[a].size() != total) {
            throw std::invalid_argument(
                "rank_jump_test: stockReturns[a] size != nRets*nDays");
        }
    }
    if (nRets == 0 || nDays == 0) {
        throw std::invalid_argument(
            "rank_jump_test: nRets and nDays must be > 0");
    }
    if (coarseFreq < 1) {
        throw std::invalid_argument(
            "rank_jump_test: coarseFreq must be >= 1");
    }

    result.nAssets = nAssets;

    // --- 步骤 1: 市场跳跃检测 (D17) ---
    auto marketJD = detail::jump_detection(marketReturns, alphaMarket, nRets, nDays);
    result.jumpIndices = marketJD.jumpIndices;

    if (marketJD.jumpIndices.empty()) {
        // R L1051-1054: 无跳跃, 返回空
        return result;
    }

    // 过滤: 跳跃索引 + coarseFreq - 1 必须不越界
    std::vector<int> validJumps;
    validJumps.reserve(marketJD.jumpIndices.size());
    for (int idx : marketJD.jumpIndices) {
        if (static_cast<Size>(idx) + static_cast<Size>(coarseFreq) - 1 < total) {
            validJumps.push_back(idx);
        }
    }
    if (validJumps.empty()) {
        result.jumpIndices.clear();
        return result;
    }
    result.jumpIndices = validJumps;
    result.nJumps = validJumps.size();

    // --- 步骤 2: 个股 jumpDetection (取 Un) ---
    // R: stockJumpDetections[,j] <- jumpDetection(stockReturns[,j], alpha[2], ...)[["Un"]]
    std::vector<std::vector<Real>> stockUn(nAssets, std::vector<Real>(total));
    for (Size a = 0; a < nAssets; ++a) {
        auto jd = detail::jump_detection(stockReturns[a], alphaStock, nRets, nDays);
        stockUn[a] = jd.Un;
    }

    // --- 步骤 3: jumps 累积 (D18) ---
    // R: jumps = sum(stockReturns[jumpIndices+i,]) for i=0..coarseFreq-1
    //   jumps 初始 = stockReturns[jumpIndices,] (nJumps x nAssets)
    //   累加 stockReturns[jumpIndices+i,] for i=1..coarseFreq-1
    //   最后 t(jumps) → nAssets x nJumps
    const Size nJumps = validJumps.size();
    std::vector<Real> jumps(nAssets * nJumps, 0.0);  // row-major: jumps[a*nJumps + j]
    for (Size a = 0; a < nAssets; ++a) {
        for (Size j = 0; j < nJumps; ++j) {
            Real sum = 0.0;
            for (int i = 0; i < coarseFreq; ++i) {
                Size idx = static_cast<Size>(validJumps[j]) + static_cast<Size>(i);
                sum += stockReturns[a][idx];
            }
            jumps[a * nJumps + j] = sum;
        }
    }

    // --- 步骤 4: SVD 全分解 (D19) ---
    // R: decomp <- svd(jumps, nu=nrow, nv=ncol)
    //   jumps: nAssets x nJumps, m=nAssets, n=nJumps
    if (nAssets == 1 && nJumps == 1) {
        throw std::runtime_error(
            "rank_jump_test: SVD cannot be calculated for 1x1 matrix");
    }

    auto svd = detail::svd_full(jumps, nAssets, nJumps);

    // U2 = U[:, rank:nAssets] (0-based, R (rank+1):n 1-based)
    // V2 = V[:, rank:nJumps]
    // singularValues = d[rank:min(nAssets,nJumps)]^2
    const Size svd_k = svd.d.size();  // min(nAssets, nJumps)
    if (static_cast<Size>(rank) >= svd_k) {
        throw std::invalid_argument(
            "rank_jump_test: rank must be < min(nAssets, nJumps)");
    }

    const Size u2_cols = nAssets - static_cast<Size>(rank);
    const Size v2_cols = nJumps - static_cast<Size>(rank);
    const Size sv_count = svd_k - static_cast<Size>(rank);

    // U2: nAssets x u2_cols (从 U 的第 rank 列开始)
    std::vector<Real> U2(nAssets * u2_cols);
    for (Size i = 0; i < nAssets; ++i) {
        for (Size j = 0; j < u2_cols; ++j) {
            U2[i * u2_cols + j] = svd.U[i * nAssets + (static_cast<Size>(rank) + j)];
        }
    }
    // V2: nJumps x v2_cols
    std::vector<Real> V2(nJumps * v2_cols);
    for (Size i = 0; i < nJumps; ++i) {
        for (Size j = 0; j < v2_cols; ++j) {
            V2[i * v2_cols + j] = svd.V[i * nJumps + (static_cast<Size>(rank) + j)];
        }
    }

    // singularValues = d[rank:svd_k]^2
    std::vector<Real> singularValues(sv_count);
    for (Size i = 0; i < sv_count; ++i) {
        singularValues[i] = svd.d[static_cast<Size>(rank) + i]
                          * svd.d[static_cast<Size>(rank) + i];
    }

    // --- 步骤 5: testStatistic (D20) ---
    // R: testStatistic[i] = sum(BoxCox__(singularValues, a))
    result.testStatistic.resize(boxCox.size());
    for (Size i = 0; i < boxCox.size(); ++i) {
        Real sum = 0.0;
        for (Size k = 0; k < sv_count; ++k) {
            sum += detail::box_cox(singularValues[k], boxCox[i]);
        }
        result.testStatistic[i] = sum;
    }

    // --- 步骤 6: bootstrap (D21) ---
    // R: dxc = pmax(pmin(stockReturns, stockUn), -stockUn)
    std::vector<std::vector<Real>> dxc(nAssets, std::vector<Real>(total));
    for (Size a = 0; a < nAssets; ++a) {
        for (Size i = 0; i < total; ++i) {
            Real val = stockReturns[a][i];
            Real un = stockUn[a][i];
            dxc[a][i] = std::max(std::min(val, un), -un);
        }
    }

    // R bootstrap: p = ncol(jumps) = nJumps
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    std::vector<Real> simTestStats(static_cast<Size>(nBoot));

    for (int b = 0; b < nBoot; ++b) {
        // zetaStar: nAssets x nJumps (row-major)
        std::vector<Real> zetaStar(nAssets * nJumps, 0.0);

        for (Size i = 0; i < nJumps; ++i) {
            int jmp0 = validJumps[static_cast<size_t>(i)];  // 0-based global index

            Size pos0;  // 0-based intraday position
            Size leftKN, rightKN;

            if (dontTestAtBoundaries) {
                // R: pos = ((jmp-1) %% nRets) + 1  (1-based)
                // C++ 0-based: pos0 = jmp0 % nRets
                pos0 = static_cast<Size>(jmp0) % nRets;
                // R: leftKN = min(localWindow, pos-1), pos-1 = pos0 (0-based)
                leftKN = std::min(static_cast<Size>(localWindow), pos0);
                // R: rightKN = min(localWindow, nRets - pos), nRets-pos = nRets-pos0-1
                rightKN = std::min(static_cast<Size>(localWindow),
                                   nRets - 1 - pos0);
            } else {
                leftKN = static_cast<Size>(localWindow);
                rightKN = static_cast<Size>(localWindow);
            }

            // jumpLeft, jumpRight (R: ceiling(runif(1) * KN))
            // R 1-based, jumpLeft in 1..leftKN (if leftKN > 0), else 0
            int jumpLeft = 0;
            int jumpRight = 0;
            if (leftKN > 0) {
                Real u = unif(rng);
                jumpLeft = 1 + static_cast<int>(u * static_cast<Real>(leftKN));
                if (jumpLeft > static_cast<int>(leftKN))
                    jumpLeft = static_cast<int>(leftKN);
                if (jumpLeft < 1) jumpLeft = 1;
            }
            if (rightKN > 0) {
                Real u = unif(rng);
                jumpRight = 1 + static_cast<int>(u * static_cast<Real>(rightKN));
                if (jumpRight > static_cast<int>(rightKN))
                    jumpRight = static_cast<int>(rightKN);
                if (jumpRight < 1) jumpRight = 1;
            }

            // dxcLeft[a] = dxc[a][jmp0 - jumpLeft], dxcRight[a] = dxc[a][jmp0 + jumpRight]
            Size idxLeft = static_cast<Size>(jmp0 - jumpLeft);
            Size idxRight = static_cast<Size>(jmp0 + jumpRight);

            // kappaStar = runif(1)
            Real kappaStar = unif(rng);

            // zetaStar[:,i] = sqrt(kappa)*dxcLeft + sqrt(C-kappa)*dxcRight
            Real sqrt_kappa = std::sqrt(kappaStar);
            Real sqrt_C_minus_kappa = std::sqrt(
                static_cast<Real>(coarseFreq) - kappaStar);

            for (Size a = 0; a < nAssets; ++a) {
                Real dxcL = dxc[a][idxLeft];
                Real dxcR = dxc[a][idxRight];
                zetaStar[a * nJumps + i] = sqrt_kappa * dxcL + sqrt_C_minus_kappa * dxcR;
            }
        }

        // tmp = t(U2) %*% zetaStar %*% V2
        // t(U2): u2_cols x nAssets
        // zetaStar: nAssets x nJumps
        // V2: nJumps x v2_cols
        // tmp: u2_cols x v2_cols
        // simTestStat = sum(tmp^2)

        // 中间: M = t(U2) %*% zetaStar  (u2_cols x nJumps)
        std::vector<Real> M(u2_cols * nJumps, 0.0);
        for (Size r = 0; r < u2_cols; ++r) {
            for (Size c = 0; c < nJumps; ++c) {
                Real s = 0.0;
                for (Size k = 0; k < nAssets; ++k) {
                    // t(U2)[r][k] = U2[k * u2_cols + r]
                    s += U2[k * u2_cols + r] * zetaStar[k * nJumps + c];
                }
                M[r * nJumps + c] = s;
            }
        }

        // tmp = M %*% V2  (u2_cols x v2_cols)
        Real simStat = 0.0;
        for (Size r = 0; r < u2_cols; ++r) {
            for (Size c = 0; c < v2_cols; ++c) {
                Real s = 0.0;
                for (Size k = 0; k < nJumps; ++k) {
                    s += M[r * nJumps + k] * V2[k * v2_cols + c];
                }
                simStat += s * s;
            }
        }
        simTestStats[static_cast<Size>(b)] = simStat;
    }

    // --- 步骤 7: criticalValues = quantile(simTestStat, quantiles) ---
    // R quantile type 7 (默认): h = (n-1)*p, Q = x[floor(h)] + (h-floor(h))*(x[floor(h)+1]-x[floor(h)])
    std::sort(simTestStats.begin(), simTestStats.end());
    const Real n_real = static_cast<Real>(nBoot);

    result.criticalValues.resize(quantiles.size());
    for (Size q = 0; q < quantiles.size(); ++q) {
        Real p = quantiles[q];
        if (nBoot == 1) {
            result.criticalValues[q] = simTestStats[0];
            continue;
        }
        Real h = (n_real - 1.0) * p;
        Size lo = static_cast<Size>(std::floor(h));
        Size hi = (lo + 1 < static_cast<Size>(nBoot)) ? lo + 1 : lo;
        Real frac = h - static_cast<Real>(lo);
        result.criticalValues[q] = simTestStats[lo]
                                 + frac * (simTestStats[hi] - simTestStats[lo]);
    }

    return result;
}

} // namespace hfecon
} // namespace v1
} // namespace cpphub
