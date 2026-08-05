// SOURCE: PHASE6_IMPLEMENTATION_PLAN §6.1 任务 4.2 - 配对 (Paired/Nonparametric) Bootstrap
// 教材锚点:
//   - Efron-Tibshirani (1993) "An Introduction to the Bootstrap" Ch.7-8
//   - Davison-Hinkley (1997) §3.2 (配对 Bootstrap, 非参数重采样)
//   - Cameron-Trivedi (2005) §11.2 (配对 Bootstrap vs 残差 Bootstrap)
//
// 算法:
//   1. 原样本: (y_i, x_i), i = 1..N
//   2. 原估计: θ̂ = estimator.estimate(data)
//   3. 对 b = 1..B:
//      a. 有放回采样 N 个索引 i*_1, ..., i*_N ~ Uniform{1, ..., N}
//      b. 构造重采样数据: y*_b = (y_{i*_1}, ..., y_{i*_N}), X*_b 同理
//      c. 估计 θ̂*_b = estimator.estimate(data*_b)
//   4. Bootstrap 协方差: V_boot = (1/(B-1)) Σ (θ̂*_b - θ̄*)(θ̂*_b - θ̄*)'
//   5. 置信区间: percentile CI (默认) / BCa (可选, 需 jackknife)
//
// 排幻觉点:
//   H-006: 配对 Bootstrap 是对 (x_i, y_i) 配对采样, 不是分别采样 x 和 y
//          (分别采样会破坏 x-y 之间的相关性, 导致 V_boot 严重偏小)
//   H-008: V_boot 的分母是 (B-1) 而非 B (样本方差的无偏估计)
//          (1/B 会导致 V_boot 偏小约 1/B 比例, B=999 时偏误 ~0.1%)
//   H-009: 失败的 replicate (estimate 抛异常或返回 NaN) 必须计入 n_failed,
//          不参与均值/协方差计算, 但 n_replicates 仍为请求的 B
//   H-012: bootstrap_samples 必须返回每次重抽样的系数向量,
//          不能只返回汇总统计量 (用于后续 BCa/plot/诊断)
//   H-013: RNG 使用 Philox4x64(seed, b) 分块, 保证每个 replicate 独立可复现
//          (ADR-004, 即使并行也是确定性结果)
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
// PairedBootstrap - 配对 (非参数) Bootstrap
//
// 适用场景:
//   - 横截面数据 (CrossSectionData)
//   - 估计量为光滑函数 (μ, β, σ² 等可微统计量)
//   - 大样本 (N >= 20, 经验法则; 小样本考虑 Wild Bootstrap)
//
// 不适用场景:
//   - 时间序列 (有自相关, 用 BlockBootstrap)
//   - 面板数据 (有聚类结构, 用 ClusterBootstrap)
//   - 边界估计 (分位数/极值, 用 m-out-of-n Bootstrap)
//
// 用法:
//   PairedBootstrap engine;
//   engine.setConfidenceLevel(0.05);
//   OLSEstimator ols(CovarianceType::HC1);
//   BootstrapResult r = engine.resample(ols, data, 999, 42);
//
// 注: 为保证可复现性, 同 (seed, B, data, estimator) 组合应产出位精确相同结果
// =============================================================================
class PairedBootstrap : public BootstrapEngine {
public:
    /// @brief 默认构造 (alpha=0.05, seed=42)
    PairedBootstrap() = default;

    /// @brief 执行配对 Bootstrap 重采样
    /// @param estimator 估计器 (调用其 estimate() 方法)
    /// @param data 原始数据 (仅 CrossSectionData)
    /// @param n_replicates 重抽样次数 B (默认 999)
    /// @param seed 随机种子 (默认 42)
    /// @return Bootstrap 结果 (系数均值/标准差/协方差/CI/samples)
    BootstrapResult resample(Estimator& estimator, const EconData& data,
                              Size n_replicates = 999,
                              std::uint64_t seed = 42) override {
        // -----------------------------------------------------------------
        // 1. 提取原数据 (X, y)
        // -----------------------------------------------------------------
        MatrixXD X_orig;
        VectorXD y_orig;
        extract_cross_section(data, X_orig, y_orig);

        const Size N = X_orig.rows();
        const Size k = X_orig.cols();
        if (N < 2) {
            throw std::invalid_argument(
                "PairedBootstrap::resample: requires N >= 2 observations");
        }
        if (k == 0) {
            throw std::invalid_argument(
                "PairedBootstrap::resample: design matrix has 0 columns");
        }
        if (n_replicates < 2) {
            throw std::invalid_argument(
                "PairedBootstrap::resample: requires n_replicates >= 2");
        }

        // -----------------------------------------------------------------
        // 2. 原样本估计 θ̂ (用于 BCa bias correction, 不参与 V_boot 计算)
        // -----------------------------------------------------------------
        EstimationResult orig_result;
        try {
            orig_result = estimator.estimate(data);
        } catch (const std::exception&) {
            // 原样本估计失败: 无法做 bias correction, 但仍可计算 V_boot
            orig_result.coefficients = VectorXD(k);
            for (Size i = 0; i < k; ++i) orig_result.coefficients(i) = 0.0;
        }

        // 提取 x_names / y_name 用于构造重采样数据 (保持元信息)
        std::vector<std::string> x_names;
        std::string y_name;
        const auto* cs_ptr = std::get_if<CrossSectionData>(&data);
        if (cs_ptr) {
            x_names = cs_ptr->x_names;
            y_name = cs_ptr->y_name;
        }

        // -----------------------------------------------------------------
        // 3. Bootstrap 重采样主循环
        // -----------------------------------------------------------------
        std::vector<VectorXD> samples;
        samples.reserve(n_replicates);
        Size n_failed = 0;

        for (Size b = 0; b < n_replicates; ++b) {
            // H-013: 每个 replicate 独立种子 (Philox 分块, 确定性可复现)
            Philox4x64 rng(seed, b);

            // 采样 N 个索引 (有放回, 配对采样)
            std::vector<Size> idx(N);
            for (Size i = 0; i < N; ++i) {
                idx[i] = randint(rng, N);
            }

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

            // 估计 θ̂*_b (捕获异常, 计入 n_failed)
            try {
                EconData data_b = make_cross_section(X_b, y_b, x_names, y_name);
                EstimationResult r_b = estimator.estimate(data_b);

                // 验证结果有效 (系数非 NaN/Inf, 维度匹配)
                if (r_b.coefficients.size() != k ||
                    !r_b.coefficients.eigen().allFinite()) {
                    ++n_failed;
                    continue;
                }

                samples.push_back(r_b.coefficients);
            } catch (const std::exception&) {
                // H-009: 失败计入 n_failed, 不参与 V_boot
                ++n_failed;
            }
        }

        // -----------------------------------------------------------------
        // 4. 计算 Bootstrap 统计量 (均值/标准差/协方差)
        // -----------------------------------------------------------------
        BootstrapResult result;
        result.n_replicates = n_replicates;
        result.n_failed = n_failed;
        result.bootstrap_samples = samples;

        if (samples.empty()) {
            // 全部失败: 返回 NaN
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

        // H-008: V_boot = (1/(B-1)) Σ (θ̂*_b - θ̄*)(θ̂*_b - θ̄*)'
        compute_stats(samples, result.coef_mean, result.coef_vcov);

        // 标准差 = sqrt(diag(V_boot))
        result.coef_std = VectorXD(k);
        for (Size i = 0; i < k; ++i) {
            result.coef_std(i) = std::sqrt(result.coef_vcov(i, i));
        }

        // -----------------------------------------------------------------
        // 5. 置信区间 (percentile CI, 对每个系数单独计算)
        //    注: BootstrapResult.lower_ci/upper_ci 为第 0 个系数的 CI
        //        (向后兼容, 多系数 CI 通过 bootstrap_samples 自行计算)
        // -----------------------------------------------------------------
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

    std::string name() const override { return "PairedBootstrap"; }
};

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
