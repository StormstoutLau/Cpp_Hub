// =============================================================================
// mixed_freq_data.hpp - 混频对齐 mls + 设计矩阵 (spec §6.1 / MD3 / C4)
//
// Phase 7C v1.7 M4 (PHASE7C_SPEC.md v1.2 §6.1)
//
// 对照库: R midasr 0.9 mls — 源码一手实录 (2026-08-18 print(mls)):
//   mls(x, k, m):
//     n <- n.x %/% m; 要求 m | n.x ("Incomplete high frequency data")
//     kmax <- max(k); mk <- min(kmax+1, 0)→0 (k 正则)
//     idx <- m·(((kmax−1)%/%m + 1):(n − 0))        # R: (k−1)%/%m with k=kmax+1
//     行 r (1-based) ↔ j = ⌊kmax/m⌋ + r; 列 h 的元素 = x[m·j − h]
//     输出 (n − ⌊kmax/m⌋) 行 × length(k) 列 (选取 res[, k+1])
//   → lag 0 列 = x[m·j] = 低频期 j 内最后一个高频观测 (期末对齐, MD3);
//     "期初可得信息" 预测须从 h ≥ m 起窗 (spec §6.1)
//
//   实测算例 (R 实测 2026-08-18): mls(x=1:12, 0:2, 3) →
//     行 1..4, 列 lag0 = x[3,6,9,12] (期末 ✓)
//     kmax=3 (lags 0:3) → 首期丢弃 (行 j 从 2 起, 滞后窗完整原则)
//
// ⚠️ spec §6.3 签名 Eigen::MatrixXd design_matrix — 与 §1.4-6
//   "M0/M1/M4 头文件零 Eigen include" 矛盾 → 按总则执行:
//   std::vector<std::vector<Real>> (spec 勘误, 总则优先)
// =============================================================================

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace midas {

/// 混频数据载体 (spec §6.1 / C4)
struct MixedFreqData {
    std::vector<Real> y;    ///< 低频响应 (n 期)
    std::vector<Real> x;    ///< 高频解释变量 (N = n·m)
    Size m = 1;             ///< 频率比 (m ∤ len(x) ⇒ invalid_argument)
    Size max_low_lag = 0;   ///< 低频 AR 最大滞后 (MIDAS-AR 用)

    /// 构造校验: 等长约束 N = n·m、无 NaN、m ≥ 1
    void validate() const {
        if (m < 1) throw std::invalid_argument("MixedFreqData: m >= 1");
        if (x.empty() || y.empty()) {
            throw std::invalid_argument("MixedFreqData: empty series");
        }
        if (x.size() % m != 0) {
            throw std::invalid_argument(
                "MixedFreqData: m must divide len(x) (midasr mls 同)");
        }
        if (x.size() / m != y.size()) {
            throw std::invalid_argument(
                "MixedFreqData: len(x) = len(y)*m required");
        }
        for (Real v : y) {
            if (!std::isfinite(v)) {
                throw std::invalid_argument("MixedFreqData: NaN in y");
            }
        }
        for (Real v : x) {
            if (!std::isfinite(v)) {
                throw std::invalid_argument("MixedFreqData: NaN in x");
            }
        }
    }
};

/// @brief mls 单列 (midasr 语义): 列 = x[m·t − h], t = 1..n (期末对齐系)
/// @param x 高频序列 (m | len(x))
/// @param h 滞后 (0 = 期末 x_tm; MD3)
/// @param m 频率比
/// @return 长度 n 的列 (元素 t−1 = x[m·t − h − 1] 0-based)
/// @throws std::invalid_argument m ∤ len(x) 或越界
inline std::vector<Real> mls_column(const std::vector<Real>& x, Size h,
                                    Size m) {
    if (m == 0 || x.size() % m != 0) {
        throw std::invalid_argument("mls_column: m must divide len(x)");
    }
    const Size n = x.size() / m;
    std::vector<Real> col(n);
    // 1-based: 元素 = x[m·t − h]; 0-based 索引 = m·t − h − 1 ≥ 0 需校验
    for (Size t = 1; t <= n; ++t) {
        const Size idx = m * t - h;  // 1-based
        if (idx < 1) {
            throw std::invalid_argument(
                "mls_column: lag exceeds available pre-sample (h too big)");
        }
        col[t - 1] = x[idx - 1];
    }
    return col;
}

/// @brief DL/U-MIDAS 设计矩阵 (滞后窗列 + 有效行对齐)
/// @param data 已 validate 的混频数据
/// @param k_high 高频滞后数 (列 h = 0..k_high−1, lag0 = 期末, MD3)
/// @param h_start 列起始滞后 (默认 0; "期初信息" 预测用 h_start ≥ m)
/// @return n_eff × k_high 矩阵; 行 r ↔ 低频期 t = j0 + r (1-based, 0-based y 索引 j0+r−1)
///   其中 j0 = ⌊(h_start+k_high−1)/m⌋ + 1 (midasr idx 起始; kmax 滞后窗完整)
inline std::vector<std::vector<Real>> design_matrix(const MixedFreqData& data,
                                                    Size k_high,
                                                    Size h_start = 0) {
    data.validate();
    if (k_high == 0) throw std::invalid_argument("design_matrix: k_high ≥ 1");
    const Size kmax = h_start + k_high - 1;
    const Size n = data.y.size();
    const Size j0 = kmax / data.m + 1;
    if (j0 > n) {
        throw std::invalid_argument(
            "design_matrix: burn-in exceeds sample (k_high too large)");
    }
    const Size n_eff = n - j0 + 1;
    std::vector<std::vector<Real>> X(n_eff,
                                     std::vector<Real>(k_high, 0.0));
    for (Size r = 0; r < n_eff; ++r) {
        const Size j = j0 + r;  // 1-based 低频期
        for (Size c = 0; c < k_high; ++c) {
            const Size h = h_start + c;
            const Size idx = data.m * j - h;  // 1-based
            X[r][c] = data.x[idx - 1];
        }
    }
    return X;
}

/// @brief 设计矩阵对齐的低频 y 切片 (行 r ↔ y[j0+r], 0-based 索引 j0+r−1)
inline std::vector<Real> aligned_y(const MixedFreqData& data, Size k_high,
                                   Size h_start = 0, Size low_lags = 0) {
    data.validate();
    const Size kmax_all = h_start + k_high - 1 + low_lags * data.m;
    const Size n = data.y.size();
    const Size j0 = kmax_all / data.m + 1;
    if (j0 > n) {
        throw std::invalid_argument("aligned_y: burn-in exceeds sample");
    }
    const Size n_eff = n - j0 + 1;
    std::vector<Real> out(n_eff);
    for (Size r = 0; r < n_eff; ++r) {
        out[r] = data.y[j0 + r - 1];
    }
    return out;
}

/// @brief 低频自身滞后列 (MIDAS-AR 用; midasr mls(y, k, 1) 语义)
/// @param y 低频序列; @param lag 滞后阶 (≥ 1)
/// @return 长度 n−lag 的列 (元素 r = y[r] 即 0-based 前移 lag)
inline std::vector<Real> ar_lag_column(const std::vector<Real>& y, Size lag) {
    if (lag == 0 || lag >= y.size()) {
        throw std::invalid_argument("ar_lag_column: 1 ≤ lag < n");
    }
    return std::vector<Real>(y.begin(), y.end() - static_cast<std::ptrdiff_t>(lag));
}

}  // namespace midas
}  // namespace timeseries
}  // namespace v1
}  // namespace cpphub
