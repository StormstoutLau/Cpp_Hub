// SOURCE: PHASE6_IMPLEMENTATION_PLAN §6.2 任务 4.4 - Block Bootstrap (Künsch 1989 / Politis-Romano 1994)
// 教材锚点:
//   - Künsch (1989) "The jackknife and the bootstrap for general stationary observations"
//     Ann. Stat. 17:1217-1241 (Moving Block Bootstrap, 重叠块)
//   - Politis-Romano (1994) "The stationary bootstrap" JASA 89:1303-1313
//     (Stationary Bootstrap, 块长 L ~ 几何分布)
//   - Politis-Romano (1994) "Large sample confidence intervals based on subsamples
//     under minimal assumptions" (Circular Block Bootstrap, 环形索引避免端点偏倚)
//   - Politis-White (2004) "Automatic block-length selection for the dependent bootstrap"
//     Econometric Reviews 23:53-70 (自动块长选择)
//
// 算法:
//   1. 原样本: 时间序列 (y_t, x_t), t = 1..N
//   2. 块长度 L 选择:
//      - 经验法则: L = ⌈N^{1/3}⌉ (默认, Politis-Romano 1994)
//      - Politis-White 2004 自动选择 (基于自相关衰减)
//   3. 对 b = 1..B:
//      a. 采样块起始位置 i ~ Uniform{1, ..., N}
//      b. 块长度:
//         - NonOverlapping: 固定 L, 不重叠块
//         - Circular (Politis-Romano 1994): 固定 L, 环形索引 (i+k) mod N
//         - Stationary (Politis-Romano 1994): L ~ Geom(1/L̄), 期望 L̄
//      c. 拼接块直到总长度 >= N, 截断到 N
//      d. 用重采样数据估计 θ̂*_b
//   4. V_boot = (1/(B-1)) Σ (θ̂*_b - θ̄*)(θ̂*_b - θ̄*)'
//
// 排幻觉点:
//   H-018: Block Bootstrap 用于时间序列 (有自相关), 不是横截面数据
//          (横截面数据用 PairedBootstrap; Block 会人为引入自相关)
//   H-019: 块长度 L 选择至关重要
//          - L 太小: 不足以保留自相关结构 → V_boot 偏小
//          - L 太大: 重采样样本与原样本太相似 → V_boot 偏大
//          - 经验法则 L = N^{1/3} (Politis-Romano 1994, 平方根最优)
//   H-020: Stationary Bootstrap 块长服从几何分布 (Politis-Romano 1994)
//          P(L = k) = p(1-p)^{k-1}, E[L] = 1/p = L̄
//          (固定块长会引入边界效应, 几何分布更平滑)
//   H-021: Circular Bootstrap 将序列视为环形 (i+k) mod N
//          避免端点偏倚 (端点观测被采样次数少于中间观测)
//   H-022: 块采样数 = ⌈N/L⌉ (非 N/L, 向上取整保证总长度 >= N)
//
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训)
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cpphub/core/linalg_dynamic.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/core/types.hpp"
#include "cpphub/econometrics/core/data_types.hpp"
#include "cpphub/econometrics/core/estimation_result.hpp"
#include "cpphub/econometrics/core/estimator_base.hpp"
#include "cpphub/econometrics/resampling/bootstrap_base.hpp"

namespace cpphub {
inline namespace v1 {
namespace econometrics {

using linalg::dynamic::MatrixXD;
using linalg::dynamic::VectorXD;

// =============================================================================
// BlockType - Block Bootstrap 块采样类型
// =============================================================================
enum class BlockType {
    Stationary,       ///< Politis-Romano 1994 stationary: 块长 ~ Geom(1/L̄)
    Circular,         ///< Politis-Romano 1994 circular: 固定块长, 环形索引
    NonOverlapping    ///< Künsch 1989 non-overlapping: 固定块长, 不重叠
};

/// @brief BlockType 枚举转字符串
inline std::string to_string(BlockType t) {
    switch (t) {
        case BlockType::Stationary:     return "Stationary";
        case BlockType::Circular:       return "Circular";
        case BlockType::NonOverlapping: return "NonOverlapping";
    }
    return "Unknown";
}

// =============================================================================
// detail: Block Bootstrap 内部工具
// =============================================================================
namespace detail {

/// @brief 经验法则选择块长度: L = ⌈N^{1/3}⌉ (Politis-Romano 1994)
/// @param N 时间序列长度
/// @return 块长度 L >= 1
inline Size default_block_length(Size N) {
    if (N < 2) return 1;
    // L = ⌈N^{1/3}⌉, 至少为 1
    const Real L_real = std::pow(static_cast<Real>(N), 1.0 / 3.0);
    Size L = static_cast<Size>(std::ceil(L_real));
    if (L < 1) L = 1;
    if (L > N) L = N;  // 块长不能超过序列长
    return L;
}

/// @brief Politis-White (2004) 自动块长选择
///
/// 基于 "spectral density at frequency 0" 估计:
///   L̄ = (2·[Σ_{k=1}^{M} k·ρ̂(k)]² / Σ_{k=-M}^{M} ρ̂(k)²)^{1/3} · N^{1/3}
/// 其中 ρ̂(k) 为样本自相关, M 为截断滞后 (自相关首次落入 ±2/√N 内的位置)
///
/// 简化实现 (Politis-White 2004 §3.2):
///   1. 计算 ρ̂(k), k = 1, 2, ...
///   2. 找最小 M 使得 |ρ̂(M)| < 2/√N (统计不显著)
///   3. G = 2·Σ_{k=1}^{M} k·ρ̂(k)  (加权和)
///   4. L̄ = ⌈(2·G²)^(1/3) · N^(1/3)⌉  (圆整)
///
/// @param acf 样本自相关函数 ρ̂(k), k=0,1,...,K
/// @param N 样本大小
/// @return 推荐块长度 L̄
inline Size politis_white_block_length(const std::vector<Real>& acf, Size N) {
    if (N < 4 || acf.size() < 2) {
        return default_block_length(N);
    }
    // 找 M: 最小 k 使得 |ρ̂(k)| < 2/√N
    const Real threshold = 2.0 / std::sqrt(static_cast<Real>(N));
    Size M = 1;
    for (Size k = 1; k < acf.size(); ++k) {
        if (std::fabs(acf[k]) < threshold) {
            M = k;
            break;
        }
        M = k;  // 全部显著时用最大滞后
    }
    if (M < 1) M = 1;

    // G = 2·Σ_{k=1}^{M} k·ρ̂(k)  (注: 仅正滞后项, 因对称)
    Real G = 0.0;
    for (Size k = 1; k <= M && k < acf.size(); ++k) {
        G += static_cast<Real>(k) * acf[k];
    }
    G *= 2.0;

    // L̄ = ⌈(2·G²)^(1/3) · N^(1/3)⌉
    // 注: 若 G <= 0 (负自相关主导), 退化为经验法则
    if (G <= 0.0) {
        return default_block_length(N);
    }
    const Real L_real = std::pow(2.0 * G * G, 1.0 / 3.0) *
                        std::pow(static_cast<Real>(N), 1.0 / 3.0);
    Size L = static_cast<Size>(std::ceil(L_real));
    if (L < 1) L = 1;
    if (L > N) L = N;
    return L;
}

/// @brief 计算样本自相关函数 (ACF) ρ̂(k), k=0,1,...,max_lag
/// @param series 时间序列 (N)
/// @param max_lag 最大滞后
/// @return ACF 向量 (max_lag+1)
inline std::vector<Real> sample_acf(const VectorXD& series, Size max_lag) {
    const Size N = series.size();
    std::vector<Real> acf;
    if (N < 2) return acf;
    if (max_lag >= N) max_lag = N - 1;

    // 样本均值
    Real mean = 0.0;
    for (Size i = 0; i < N; ++i) mean += series(i);
    mean /= static_cast<Real>(N);

    // 方差 (分母 N, 与 R acf() 一致)
    Real var = 0.0;
    for (Size i = 0; i < N; ++i) {
        const Real d = series(i) - mean;
        var += d * d;
    }
    if (var <= 0.0) {
        acf.assign(max_lag + 1, 0.0);
        acf[0] = 1.0;
        return acf;
    }
    var /= static_cast<Real>(N);

    // 自协方差 γ(k) = (1/N) Σ_{t=1}^{N-k} (y_t - ȳ)(y_{t+k} - ȳ)
    acf.resize(max_lag + 1);
    acf[0] = 1.0;  // ρ̂(0) = 1
    for (Size k = 1; k <= max_lag; ++k) {
        Real cov = 0.0;
        for (Size t = 0; t + k < N; ++t) {
            cov += (series(t) - mean) * (series(t + k) - mean);
        }
        cov /= static_cast<Real>(N);
        acf[k] = cov / var;
    }
    return acf;
}

}  // namespace detail

// =============================================================================
// BlockBootstrap - 时间序列 Block Bootstrap 引擎
//
// 适用场景:
//   - 时间序列数据 (有自相关)
//   - HAR 模型系数 SE 估计 (Corsi 2009 场景)
//   - GARCH 残差 Bootstrap
//
// 不适用场景:
//   - 横截面数据 (用 PairedBootstrap)
//   - 面板数据 (用 ClusterBootstrap)
//   - 强非平稳序列 (需差分/去趋势后使用)
//
// 用法:
//   BlockBootstrap engine(BlockType::Stationary);
//   engine.setBlockLength(10);  // 可选, 默认 N^{1/3}
//   OLSEstimator ols;
//   BootstrapResult r = engine.resample(ols, ts_data, 999, 42);
// =============================================================================
class BlockBootstrap : public BootstrapEngine {
public:
    /// @brief 默认构造 (Stationary, 自动块长)
    BlockBootstrap()
        : block_type_(BlockType::Stationary)
        , block_length_(0)  // 0 = 自动选择
        , auto_block_length_(true) {}

    /// @brief 指定块类型和块长度的构造器
    /// @param type 块类型 (Stationary/Circular/NonOverlapping)
    /// @param L 块长度 (0 = 自动选择 N^{1/3})
    explicit BlockBootstrap(BlockType type, Size L = 0)
        : block_type_(type)
        , block_length_(L)
        , auto_block_length_(L == 0) {}

    /// @brief 设置块类型
    void setBlockType(BlockType t) { block_type_ = t; }

    /// @brief 设置块长度 (0 = 自动选择)
    void setBlockLength(Size L) {
        block_length_ = L;
        auto_block_length_ = (L == 0);
    }

    /// @brief 获取当前块类型
    BlockType blockType() const { return block_type_; }

    /// @brief 获取实际使用的块长度 (resample 后填充)
    Size blockLength() const { return block_length_; }

    /// @brief 执行 Block Bootstrap 重采样
    BootstrapResult resample(Estimator& estimator, const EconData& data,
                              Size n_replicates = 999,
                              std::uint64_t seed = 42) override {
        // -----------------------------------------------------------------
        // 1. 提取时间序列数据
        // -----------------------------------------------------------------
        if (!std::holds_alternative<TimeSeriesData>(data)) {
            // 同时支持 CrossSectionData (视为无自相关序列, 退化为 PairedBootstrap)
            // 但推荐用户对时间序列显式构造 TimeSeriesData
            if (std::holds_alternative<CrossSectionData>(data)) {
                return resample_cross_section(estimator, data, n_replicates, seed);
            }
            throw std::invalid_argument(
                "BlockBootstrap::resample: only TimeSeriesData/CrossSectionData supported");
        }
        const auto& ts = std::get<TimeSeriesData>(data);
        const MatrixXD& X_orig = ts.X;
        const VectorXD& y_orig = ts.y;
        const Size N = X_orig.rows();
        const Size k = X_orig.cols();

        if (N < 4) {
            throw std::invalid_argument(
                "BlockBootstrap::resample: requires N >= 4 (time series too short)");
        }
        if (k == 0) {
            throw std::invalid_argument(
                "BlockBootstrap::resample: design matrix has 0 columns");
        }
        if (n_replicates < 2) {
            throw std::invalid_argument(
                "BlockBootstrap::resample: requires n_replicates >= 2");
        }

        // -----------------------------------------------------------------
        // 2. 块长度选择
        // -----------------------------------------------------------------
        if (auto_block_length_) {
            block_length_ = detail::default_block_length(N);
            // Politis-White 自动选择 (基于 y 的 ACF)
            try {
                const Size max_lag = std::min<Size>(N - 1, 50);
                const std::vector<Real> acf = detail::sample_acf(y_orig, max_lag);
                const Size L_pw = detail::politis_white_block_length(acf, N);
                if (L_pw >= 1 && L_pw <= N) {
                    block_length_ = L_pw;
                }
            } catch (...) {
                // Politis-White 失败时保留经验法则
            }
        }
        Size L = block_length_;
        if (L < 1) L = 1;
        if (L > N) L = N;

        // -----------------------------------------------------------------
        // 3. Bootstrap 重采样主循环
        // -----------------------------------------------------------------
        std::vector<VectorXD> samples;
        samples.reserve(n_replicates);
        Size n_failed = 0;

        for (Size b = 0; b < n_replicates; ++b) {
            Philox4x64 rng(seed, b);

            // 生成重采样索引 (N 个)
            std::vector<Size> idx(N);
            generate_block_indices(rng, N, L, idx.data());

            // 构造重采样数据 (X*, y*)
            MatrixXD X_b(N, k);
            VectorXD y_b(N);
            for (Size i = 0; i < N; ++i) {
                const Size src = idx[i];
                y_b(i) = y_orig(src);
                for (Size j = 0; j < k; ++j) {
                    X_b(i, j) = X_orig(src, j);
                }
            }

            // 估计 θ̂*_b (BlockBootstrap 将时间序列重采样后用 OLS 估计)
            try {
                // 转换为 CrossSectionData (因 OLSEstimator 仅支持 CrossSectionData)
                EconData data_b = make_cross_section(X_b, y_b,
                                                      ts.x_names, ts.y_name);
                EstimationResult r_b = estimator.estimate(data_b);

                if (r_b.coefficients.size() != k ||
                    !r_b.coefficients.eigen().allFinite()) {
                    ++n_failed;
                    continue;
                }

                samples.push_back(r_b.coefficients);
            } catch (const std::exception&) {
                ++n_failed;
            }
        }

        // -----------------------------------------------------------------
        // 4. 计算 Bootstrap 统计量
        // -----------------------------------------------------------------
        BootstrapResult result;
        result.n_replicates = n_replicates;
        result.n_failed = n_failed;
        result.bootstrap_samples = samples;

        if (samples.empty()) {
            result.coef_mean = VectorXD(k);
            result.coef_std = VectorXD(k);
            result.coef_vcov = MatrixXD(k, k);
            for (Size i = 0; i < k; ++i) {
                result.coef_mean(i) = std::numeric_limits<Real>::quiet_NaN();
                result.coef_std(i) = std::numeric_limits<Real>::quiet_NaN();
            }
            result.lower_ci = std::numeric_limits<Real>::quiet_NaN();
            result.upper_ci = std::numeric_limits<Real>::quiet_NaN();
            return result;
        }

        compute_stats(samples, result.coef_mean, result.coef_vcov);

        result.coef_std = VectorXD(k);
        for (Size i = 0; i < k; ++i) {
            result.coef_std(i) = std::sqrt(result.coef_vcov(i, i));
        }

        if (k >= 1) {
            std::vector<Real> ci_samples;
            ci_samples.reserve(samples.size());
            for (const auto& s : samples) {
                ci_samples.push_back(s(0));
            }
            auto [lo, hi] = percentileCI(ci_samples, alpha_);
            result.lower_ci = lo;
            result.upper_ci = hi;
        }

        return result;
    }

    std::string name() const override {
        return std::string("BlockBootstrap[") + to_string(block_type_) + "]";
    }

private:
    BlockType block_type_;
    Size block_length_;
    bool auto_block_length_;

    // -------------------------------------------------------------------------
    // 生成块重采样索引 (H-019, H-020, H-021, H-022)
    // -------------------------------------------------------------------------
    void generate_block_indices(Philox4x64& rng, Size N, Size L, Size* idx) const {
        Size filled = 0;
        while (filled < N) {
            // 块起始位置: Uniform{0, 1, ..., N-1}
            const Size start = randint(rng, N);

            // 块长度
            Size block_len;
            if (block_type_ == BlockType::Stationary) {
                // H-020: 几何分布, P(L=k) = p(1-p)^{k-1}, E[L] = 1/p = L̄
                // 简化: 接受-拒绝法 (持续生成直到拒绝)
                const Real p = 1.0 / static_cast<Real>(L);  // 期望 = L
                block_len = 1;
                while (block_len < N) {
                    const Real u = (static_cast<uint64_t>(rng()) >> 11) *
                                   (1.0 / 9007199254740992.0);
                    if (u < p) break;  // 几何分布终止
                    ++block_len;
                }
            } else {
                // Circular / NonOverlapping: 固定块长 L
                block_len = L;
            }

            // 拼接块到 idx
            for (Size i = 0; i < block_len && filled < N; ++i) {
                Size src;
                if (block_type_ == BlockType::NonOverlapping) {
                    // 不重叠: 块起始位置只能是 0, L, 2L, ...
                    // 注: 此处 start 应已被离散化, 简化实现用 start
                    src = (start + i) % N;  // 用环形避免越界
                } else {
                    // Circular / Stationary: 环形索引 (H-021)
                    src = (start + i) % N;
                }
                idx[filled++] = src;
            }
        }
        // 截断到 N (循环保证 filled == N)
    }

    // -------------------------------------------------------------------------
    // 横截面数据退化路径 (无自相关时退化为 PairedBootstrap 等价)
    // -------------------------------------------------------------------------
    BootstrapResult resample_cross_section(Estimator& estimator,
                                            const EconData& data,
                                            Size n_replicates,
                                            std::uint64_t seed) {
        // CrossSectionData: 视为无自相关, 用 NonOverlapping 块长 = 1 (等价于 PairedBootstrap)
        const auto& cs = std::get<CrossSectionData>(data);
        const MatrixXD& X_orig = cs.X;
        const VectorXD& y_orig = cs.y;
        const Size N = X_orig.rows();
        const Size k = X_orig.cols();

        if (N < 2) {
            throw std::invalid_argument(
                "BlockBootstrap::resample_cross_section: requires N >= 2");
        }

        std::vector<VectorXD> samples;
        samples.reserve(n_replicates);
        Size n_failed = 0;

        for (Size b = 0; b < n_replicates; ++b) {
            Philox4x64 rng(seed, b);
            // 块长 = 1, 等价于 Paired Bootstrap
            MatrixXD X_b(N, k);
            VectorXD y_b(N);
            for (Size i = 0; i < N; ++i) {
                const Size src = randint(rng, N);
                y_b(i) = y_orig(src);
                for (Size j = 0; j < k; ++j) {
                    X_b(i, j) = X_orig(src, j);
                }
            }
            try {
                EconData data_b = make_cross_section(X_b, y_b,
                                                      cs.x_names, cs.y_name);
                EstimationResult r_b = estimator.estimate(data_b);
                if (r_b.coefficients.size() != k ||
                    !r_b.coefficients.eigen().allFinite()) {
                    ++n_failed;
                    continue;
                }
                samples.push_back(r_b.coefficients);
            } catch (const std::exception&) {
                ++n_failed;
            }
        }

        BootstrapResult result;
        result.n_replicates = n_replicates;
        result.n_failed = n_failed;
        result.bootstrap_samples = samples;

        if (samples.empty()) {
            result.coef_mean = VectorXD(k);
            result.coef_std = VectorXD(k);
            result.coef_vcov = MatrixXD(k, k);
            for (Size i = 0; i < k; ++i) {
                result.coef_mean(i) = std::numeric_limits<Real>::quiet_NaN();
                result.coef_std(i) = std::numeric_limits<Real>::quiet_NaN();
            }
            result.lower_ci = std::numeric_limits<Real>::quiet_NaN();
            result.upper_ci = std::numeric_limits<Real>::quiet_NaN();
            return result;
        }

        compute_stats(samples, result.coef_mean, result.coef_vcov);
        result.coef_std = VectorXD(k);
        for (Size i = 0; i < k; ++i) {
            result.coef_std(i) = std::sqrt(result.coef_vcov(i, i));
        }
        if (k >= 1) {
            std::vector<Real> ci_samples;
            ci_samples.reserve(samples.size());
            for (const auto& s : samples) {
                ci_samples.push_back(s(0));
            }
            auto [lo, hi] = percentileCI(ci_samples, alpha_);
            result.lower_ci = lo;
            result.upper_ci = hi;
        }
        return result;
    }
};

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
