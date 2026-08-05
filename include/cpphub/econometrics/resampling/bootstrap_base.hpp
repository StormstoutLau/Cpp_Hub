// SOURCE: PHASE6_IMPLEMENTATION_PLAN §6.1 任务 4.1 - Bootstrap 引擎基类
// 教材锚点:
//   - Davison-Hinkley (1997) "Bootstrap Methods and Their Application"
//   - Cameron-Trivedi Ch.11 (配对/非参数/残差/Wild 完整方法)
//   - Efron-Tibshirani (1993) "An Introduction to the Bootstrap"
//
// 排幻觉点:
//   E12: Wild Bootstrap 权重分布严格按 R multiwayvcov::cluster.boot 文档
//        (默认 Rademacher, Cameron-Gelbach-Miller 2008 推荐)
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

namespace cpphub {
inline namespace v1 {
namespace econometrics {

using linalg::dynamic::MatrixXD;
using linalg::dynamic::VectorXD;

// =============================================================================
// BootstrapEngine - Bootstrap 引擎抽象基类
//
// 设计:
//   - resample() 为纯虚, 由 PairedBootstrap/WildBootstrap/BlockBootstrap/
//     ClusterBootstrap 实现
//   - percentileCI / bcaCI 为通用工具方法, 供派生类使用
//   - RNG: 复用 core/rng.hpp Philox4x64 (ADR-004, 确定性分块并行)
//
// 注: Estimator::estimate() 为非 const, 故 resample() 接受 Estimator& (非 const)
// =============================================================================
class BootstrapEngine {
public:
    virtual ~BootstrapEngine() = default;

    /// @brief 核心: 执行 Bootstrap 重采样
    /// @param estimator 估计器 (调用其 estimate() 方法)
    /// @param data 原始数据
    /// @param n_replicates 重抽样次数 B (默认 999)
    /// @param seed 随机种子 (默认 42)
    /// @return Bootstrap 结果 (系数均值/标准差/协方差/CI)
    virtual BootstrapResult resample(Estimator& estimator, const EconData& data,
                                      Size n_replicates = 999,
                                      std::uint64_t seed = 42) = 0;

    /// @brief 引擎名称
    virtual std::string name() const = 0;

    /// @brief 设置置信水平
    /// @param alpha 显著性水平 (默认 0.05 → 95% CI)
    void setConfidenceLevel(Real alpha) { alpha_ = alpha; }

    /// @brief 设置随机种子
    void setRNG(std::uint64_t seed) { seed_ = seed; }

protected:
    Real alpha_ = 0.05;            ///< 显著性水平
    std::uint64_t seed_ = 42;      ///< 随机种子

    // -------------------------------------------------------------------------
    // 百分位置信区间 (Efron's percentile method)
    //   排序样本, 取 alpha/2 和 1-alpha/2 分位数
    //   使用 R type 7 分位数 (线性插值, R 默认方法)
    // -------------------------------------------------------------------------
    std::pair<Real, Real> percentileCI(const std::vector<Real>& samples,
                                        Real alpha) const {
        if (samples.empty()) {
            return {std::numeric_limits<Real>::quiet_NaN(),
                    std::numeric_limits<Real>::quiet_NaN()};
        }
        std::vector<Real> sorted = samples;
        std::sort(sorted.begin(), sorted.end());

        const Size B = sorted.size();
        // R type 7: h = (B-1) * p, 线性插值
        auto quantile = [&](Real p) -> Real {
            if (p <= 0.0) return sorted.front();
            if (p >= 1.0) return sorted.back();
            Real h = static_cast<Real>(B - 1) * p;
            Size lo = static_cast<Size>(std::floor(h));
            Size hi = std::min(lo + 1, B - 1);
            Real frac = h - static_cast<Real>(lo);
            return sorted[lo] + frac * (sorted[hi] - sorted[lo]);
        };

        return {quantile(alpha / 2.0), quantile(1.0 - alpha / 2.0)};
    }

    // -------------------------------------------------------------------------
    // BCa 置信区间 (bias-corrected and accelerated, Davison-Hinkley §5.3)
    //   z0 = Φ⁻¹(#{θ̂*_b < θ̂} / B)         (bias correction)
    //   a  = Σ(θ̄*_j - θ̂*_j)³ / [6·(Σ(θ̄*_j - θ̂*_j)²)^{3/2}]  (acceleration)
    //   α1 = Φ(z0 + (z0 + z_{α/2})/(1 - a·(z0 + z_{α/2})))
    //   α2 = Φ(z0 + (z0 + z_{1-α/2})/(1 - a·(z0 + z_{1-α/2})))
    // -------------------------------------------------------------------------
    std::pair<Real, Real> bcaCI(const std::vector<Real>& samples, Real alpha,
                                 const std::vector<Real>& jackknife_samples,
                                 Real theta_hat) const {
        if (samples.size() < 2 || jackknife_samples.empty()) {
            return percentileCI(samples, alpha);
        }

        // Bias correction z0
        Size count_below = 0;
        for (Real s : samples) {
            if (s < theta_hat) ++count_below;
        }
        Real prop = static_cast<Real>(count_below) / static_cast<Real>(samples.size());
        // 边界: prop=0 或 1 时 z0 → ±∞, 使用截断
        prop = std::max(prop, 1e-10);
        prop = std::min(prop, 1.0 - 1e-10);
        Real z0 = inverse_normal_cdf(prop);

        // Acceleration a (jackknife)
        const Size n_jk = jackknife_samples.size();
        Real jk_mean = 0.0;
        for (Real s : jackknife_samples) jk_mean += s;
        jk_mean /= static_cast<Real>(n_jk);

        Real num = 0.0, den = 0.0;
        for (Real s : jackknife_samples) {
            Real diff = jk_mean - s;
            num += diff * diff * diff;
            den += diff * diff;
        }
        Real a = (den > 0.0)
            ? num / (6.0 * std::pow(den, 1.5))
            : 0.0;

        // BCa 分位数
        Real z_lo = inverse_normal_cdf(alpha / 2.0);
        Real z_hi = inverse_normal_cdf(1.0 - alpha / 2.0);

        auto bca_quantile = [&](Real z_alpha) -> Real {
            Real denom = 1.0 - a * (z0 + z_alpha);
            if (std::fabs(denom) < 1e-15) denom = 1e-15;
            Real z_adj = z0 + (z0 + z_alpha) / denom;
            return normal_cdf(z_adj);
        };

        Real alpha1 = bca_quantile(z_lo);
        Real alpha2 = bca_quantile(z_hi);

        // 用 BCa 分位数重新计算百分位
        return percentileCI_at(samples, alpha1, alpha2);
    }

    // -------------------------------------------------------------------------
    // 工具: 生成 [0, n) 均匀随机整数
    //   使用 53 位精度, 避免 % 的模偏差
    // -------------------------------------------------------------------------
    static Size randint(Philox4x64& rng, Size n) {
        if (n <= 1) return 0;
        uint64_t r = rng();
        double u = (r >> 11) * (1.0 / 9007199254740992.0);  // [0, 1) 53 bits
        return static_cast<Size>(u * static_cast<double>(n));
    }

    // -------------------------------------------------------------------------
    // 工具: 计算 Bootstrap 样本的均值和协方差
    //   V_boot = (1/(B-1)) Σ (θ̂*_b - θ̄*)(θ̂*_b - θ̄*)'
    // -------------------------------------------------------------------------
    static void compute_stats(const std::vector<VectorXD>& samples,
                               VectorXD& mean, MatrixXD& vcov) {
        const Size B = samples.size();
        if (B == 0) return;
        const Size k = samples[0].size();
        mean = VectorXD(k);
        vcov = MatrixXD(k, k);

        // 均值
        for (Size i = 0; i < k; ++i) mean(i) = 0.0;
        for (const auto& s : samples) {
            for (Size i = 0; i < k; ++i) mean(i) += s(i);
        }
        for (Size i = 0; i < k; ++i) mean(i) /= static_cast<Real>(B);

        // 协方差
        for (Size i = 0; i < k; ++i)
            for (Size j = 0; j < k; ++j) vcov(i, j) = 0.0;
        for (const auto& s : samples) {
            for (Size i = 0; i < k; ++i) {
                for (Size j = 0; j < k; ++j) {
                    vcov(i, j) += (s(i) - mean(i)) * (s(j) - mean(j));
                }
            }
        }
        if (B > 1) {
            const Real scale = 1.0 / static_cast<Real>(B - 1);
            for (Size i = 0; i < k; ++i)
                for (Size j = 0; j < k; ++j) vcov(i, j) *= scale;
        }
    }

    // -------------------------------------------------------------------------
    // 工具: 从 EconData 提取 CrossSectionData 的 (X, y)
    //   若 data 非 CrossSectionData, 抛 invalid_argument
    // -------------------------------------------------------------------------
    static void extract_cross_section(const EconData& data,
                                       MatrixXD& X, VectorXD& y) {
        if (!std::holds_alternative<CrossSectionData>(data)) {
            throw std::invalid_argument(
                "Bootstrap: only CrossSectionData supported by this engine");
        }
        const auto& cs = std::get<CrossSectionData>(data);
        X = cs.X;
        y = cs.y;
    }

    // -------------------------------------------------------------------------
    // 工具: 从 (X, y) 构造 CrossSectionData → EconData
    // -------------------------------------------------------------------------
    static EconData make_cross_section_data(const MatrixXD& X, const VectorXD& y,
                                             const std::vector<std::string>& x_names = {},
                                             const std::string& y_name = "") {
        return make_cross_section(X, y, x_names, y_name);
    }

private:
    // 标准正态 CDF: Φ(x) = 0.5 * erfc(-x/√2)
    static Real normal_cdf(Real x) {
        return 0.5 * std::erfc(-x / std::sqrt(2.0));
    }

    // 标准正态逆 CDF: Φ⁻¹(p), 二分法 (精度 ~1e-12)
    static Real inverse_normal_cdf(Real p) {
        if (p <= 0.0) return -std::numeric_limits<Real>::infinity();
        if (p >= 1.0) return std::numeric_limits<Real>::infinity();
        Real lo = -10.0, hi = 10.0;
        for (int i = 0; i < 100; ++i) {
            Real mid = 0.5 * (lo + hi);
            if (normal_cdf(mid) < p) lo = mid;
            else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    // 指定分位数位置的百分位 CI (BCa 内部使用)
    std::pair<Real, Real> percentileCI_at(const std::vector<Real>& samples,
                                           Real p_lo, Real p_hi) const {
        if (samples.empty()) {
            return {std::numeric_limits<Real>::quiet_NaN(),
                    std::numeric_limits<Real>::quiet_NaN()};
        }
        std::vector<Real> sorted = samples;
        std::sort(sorted.begin(), sorted.end());

        const Size B = sorted.size();
        auto quantile = [&](Real p) -> Real {
            if (p <= 0.0) return sorted.front();
            if (p >= 1.0) return sorted.back();
            Real h = static_cast<Real>(B - 1) * p;
            Size lo = static_cast<Size>(std::floor(h));
            Size hi = std::min(lo + 1, B - 1);
            Real frac = h - static_cast<Real>(lo);
            return sorted[lo] + frac * (sorted[hi] - sorted[lo]);
        };

        return {quantile(p_lo), quantile(p_hi)};
    }
};

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
