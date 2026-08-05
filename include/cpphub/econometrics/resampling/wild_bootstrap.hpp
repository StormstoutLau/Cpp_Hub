// SOURCE: PHASE6_IMPLEMENTATION_PLAN §6.1 任务 4.3 - Wild Bootstrap (Wu 1986 / Liu 1988)
// 教材锚点:
//   - Cameron-Trivedi (2005) §11.3 (Wild Bootstrap 异方差稳健)
//   - Davidson-Flachaire (2008) "Wild bootstrap with extra wild clamping" (Rademacher 推荐)
//   - Cameron-Gelbach-Miller (2008) REStat (CGM 2008, Rademacher 默认, 6 种权重对比)
//   - Mammen (1993) Ann. Stat. (两点分布, E[v]=0, Var[v]=1, 三阶矩非零)
//   - Webb (2018) Econometrics Journal (6 点分布, 小样本稳健)
//   - R multiwayvcov::cluster.boot (默认 Rademacher, CGM 2008 推荐)
//
// 算法 (Wild Bootstrap, Wu 1986):
//   1. 原样本估计: β̂ = estimator.estimate(data)
//   2. 残差: ε̂_i = y_i - x_i'β̂
//   3. 对 b = 1..B:
//      a. 对每个观测 i, 生成权重 v_i ~ F (Rademacher / Mammen / Webb6)
//      b. 构造 y*_i = x_i'β̂ + v_i · ε̂_i (固定 X, 残差乘权重)
//      c. 估计 θ̂*_b = estimator.estimate((X, y*))
//   4. V_boot = (1/(B-1)) Σ (θ̂*_b - θ̄*)(θ̂*_b - θ̄*)'
//
// 排幻觉点:
//   E12: Wild Bootstrap 权重分布严格按 R multiwayvcov::cluster.boot 文档
//        (默认 Rademacher, Cameron-Gelbach-Miller 2008 推荐)
//   H-014: Wild Bootstrap 固定 X (不重采样 X), 仅对残差加权
//          (重采样 X 会破坏条件期望结构, Wild Bootstrap 的前提是 E[ε|X]=0)
//   H-015: 残差 ε̂_i = y_i - x_i'β̂ (用原样本估计的 β̂, 不是每次重估)
//          (β̂ 固定保证 v_i 仅作用于残差, 不污染 X-β 路径)
//   H-016: y*_i = x_i'β̂ + v_i · ε̂_i (不是 v_i · y_i 或 y_i + v_i · ε̂_i)
//          (v_i · y_i 会破坏 X-y 关系; y_i + v_i · ε̂_i 会让均值漂移)
//   H-017: 三种权重分布必须满足 E[v]=0, Var[v]=1
//          (保证 y* 的条件均值 = x'β̂, 条件方差 ≈ Var[ε], 与原模型一致)
//
// 权重分布 (R multiwayvcov 1.2.3 对照):
//   Rademacher (默认): v ∈ {-1, +1} 各 w.p. 0.5
//     E[v] = 0, Var[v] = 1, 三阶矩 = 0 (对称)
//     优点: 简单, 数值稳定, 对异方差最稳健 (Davidson-Flachaire 2008)
//     缺点: 仅两点, 小样本下离散性高
//
//   Mammen (1993): v ∈ {(1-√5)/2, (1+√5)/2} w.p. {(5+√5)/10, (5-√5)/10}
//     即 v ≈ -0.6180 w.p. 0.7236, v ≈ 1.6180 w.p. 0.2764
//     E[v] = 0, Var[v] = 1, 三阶矩 = 2 (非零, 模拟偏度)
//     优点: 理论最优 (Mammen 1993 证明在 Edgeworth 展开下二阶精确)
//     缺点: 非对称, 小样本下可能引入偏误
//
//   Webb6 (2018): v ∈ {-√(3/2), -1, -√(1/2), +√(1/2), +1, +√(3/2)} 各 w.p. 1/6
//     E[v] = 0, Var[v] = (3/2 + 1 + 1/2 + 1/2 + 1 + 3/2)/6 = 6/6 = 1
//     优点: 6 点比 Rademacher 更平滑, 极少聚类 (G < 10) 时最稳健
//     缺点: 计算略贵 (6 路分支)
//     注: R multiwayvcov 1.2.3 未内置 Webb6, 需通过 wild_type=function() 自定义
//
// 约定: 头文件 #include 必须位于 namespace 外 (project_memory 教训)
#pragma once

#include <algorithm>
#include <array>
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
// WildWeightDistribution - Wild Bootstrap 权重分布枚举
// =============================================================================
enum class WildWeightDistribution {
    Rademacher,   ///< v ∈ {-1, +1} 各 w.p. 0.5 (默认, CGM 2008 推荐)
    Mammen,       ///< v ∈ {(1∓√5)/2} w.p. {(5±√5)/10} (Mammen 1993)
    Webb6         ///< v ∈ {±√(3/2), ±1, ±√(1/2)} 各 w.p. 1/6 (Webb 2018)
};

/// @brief WildWeightDistribution 枚举转字符串
inline std::string to_string(WildWeightDistribution d) {
    switch (d) {
        case WildWeightDistribution::Rademacher: return "Rademacher";
        case WildWeightDistribution::Mammen:     return "Mammen";
        case WildWeightDistribution::Webb6:       return "Webb6";
    }
    return "Unknown";
}

// =============================================================================
// detail: Wild Bootstrap 权重生成器
//
// 数值常量 (高精度预计算, 避免运行时 sqrt):
//   kSqrt5    = √5    ≈ 2.2360679774997896964091736687312762354406183596115
//   kMammenLo = (1-√5)/2 ≈ -0.61803398874989484820458683436563811772030917980576
//   kMammenHi = (1+√5)/2 ≈  1.6180339887498948482045868343656381177203091798058
//   kMammenPLo = (5+√5)/10 ≈ 0.72360679774997896964091736687312762354406183596115
//   kMammenPHi = (5-√5)/10 ≈ 0.27639320225002103035908263312687237645593816403885
//   kSqrt3_2   = √(3/2)   ≈ 1.2247448713915890490986420373529456959839616437987
//   kSqrt1_2   = √(1/2)   ≈ 0.70710678118654752440084436210484903928483593768847
// =============================================================================
namespace detail {

inline constexpr Real kSqrt5     = 2.2360679774997896964091736687312762354406183596115;
inline constexpr Real kMammenLo  = (1.0 - kSqrt5) / 2.0;   // (1-√5)/2 ≈ -0.618
inline constexpr Real kMammenHi  = (1.0 + kSqrt5) / 2.0;   // (1+√5)/2 ≈ +1.618
inline constexpr Real kMammenPLo = (5.0 + kSqrt5) / 10.0;  // (5+√5)/10 ≈ 0.7236
inline constexpr Real kMammenPHi = (5.0 - kSqrt5) / 10.0;  // (5-√5)/10 ≈ 0.2764
inline constexpr Real kSqrt3_2   = 1.2247448713915890490986420373529456959839616437987;
inline constexpr Real kSqrt1_2   = 0.70710678118654752440084436210484903928483593768847;

/// @brief 生成一个 Wild Bootstrap 权重 v
/// @param rng Philox4x64 RNG
/// @param dist 权重分布
/// @return 权重 v, 满足 E[v]=0, Var[v]=1
inline Real generate_wild_weight(Philox4x64& rng, WildWeightDistribution dist) {
    // 53 位精度均匀 [0, 1)
    auto uniform = [&rng]() -> double {
        uint64_t r = rng();
        return (r >> 11) * (1.0 / 9007199254740992.0);
    };

    switch (dist) {
        case WildWeightDistribution::Rademacher: {
            // v = +1 w.p. 0.5, v = -1 w.p. 0.5
            return (uniform() < 0.5) ? 1.0 : -1.0;
        }
        case WildWeightDistribution::Mammen: {
            // v = (1-√5)/2 w.p. (5+√5)/10 ≈ 0.7236
            // v = (1+√5)/2 w.p. (5-√5)/10 ≈ 0.2764
            return (uniform() < kMammenPLo) ? kMammenLo : kMammenHi;
        }
        case WildWeightDistribution::Webb6: {
            // 6 点等概率分布, 每 1/6
            // {-√(3/2), -1, -√(1/2), +√(1/2), +1, +√(3/2)}
            const double u = uniform() * 6.0;
            const int idx = static_cast<int>(u);  // [0, 5]
            switch (idx) {
                case 0: return -kSqrt3_2;
                case 1: return -1.0;
                case 2: return -kSqrt1_2;
                case 3: return  kSqrt1_2;
                case 4: return  1.0;
                default: return kSqrt3_2;  // idx = 5 (u 接近 6.0 时)
            }
        }
    }
    return 0.0;  // 不可达
}

}  // namespace detail

// =============================================================================
// WildBootstrap - Wild Bootstrap 引擎 (Wu 1986 / Liu 1988)
//
// 适用场景:
//   - 异方差稳健推断 (无需指定异方差形式)
//   - 小样本 N < 50 (比配对 Bootstrap 更稳健)
//   - 非线性模型 (Logit/Probit 等, 比 PairedBootstrap 数值更稳定)
//   - 聚类数 G < 10 时改用 Wild Cluster Bootstrap (v1.6+)
//
// 不适用场景:
//   - 时间序列 (有自相关, 用 BlockBootstrap)
//   - 内生性 (IV/GMM 场景, 需特殊 Wild Bootstrap 变体)
//
// 用法:
//   WildBootstrap engine(WildWeightDistribution::Rademacher);
//   engine.setConfidenceLevel(0.05);
//   OLSEstimator ols;
//   BootstrapResult r = engine.resample(ols, data, 999, 42);
// =============================================================================
class WildBootstrap : public BootstrapEngine {
public:
    /// @brief 默认构造 (Rademacher, CGM 2008 推荐)
    WildBootstrap() : weight_dist_(WildWeightDistribution::Rademacher) {}

    /// @brief 指定权重分布的构造器
    /// @param dist 权重分布 (Rademacher/Mammen/Webb6)
    explicit WildBootstrap(WildWeightDistribution dist) : weight_dist_(dist) {}

    /// @brief 设置权重分布
    void setWeightDistribution(WildWeightDistribution dist) { weight_dist_ = dist; }

    /// @brief 获取当前权重分布
    WildWeightDistribution weightDistribution() const { return weight_dist_; }

    /// @brief 执行 Wild Bootstrap 重采样
    /// @param estimator 估计器 (调用其 estimate() 方法)
    /// @param data 原始数据 (仅 CrossSectionData)
    /// @param n_replicates 重抽样次数 B (默认 999)
    /// @param seed 随机种子 (默认 42)
    /// @return Bootstrap 结果
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
                "WildBootstrap::resample: requires N >= 2 observations");
        }
        if (k == 0) {
            throw std::invalid_argument(
                "WildBootstrap::resample: design matrix has 0 columns");
        }
        if (n_replicates < 2) {
            throw std::invalid_argument(
                "WildBootstrap::resample: requires n_replicates >= 2");
        }

        // -----------------------------------------------------------------
        // 2. 原样本估计 β̂ (H-015: 用原样本估计的 β̂, 不重估)
        // -----------------------------------------------------------------
        EstimationResult orig_result = estimator.estimate(data);
        const VectorXD& beta_hat = orig_result.coefficients;
        if (beta_hat.size() != k || !beta_hat.eigen().allFinite()) {
            throw std::runtime_error(
                "WildBootstrap::resample: original estimate failed or NaN");
        }

        // 计算残差 ε̂_i = y_i - x_i'β̂ (H-015)
        const VectorXD resid = VectorXD(y_orig.eigen() - X_orig.eigen() * beta_hat.eigen());

        // 计算拟合值 x_i'β̂ (H-016: y*_i = fitted + v_i · resid)
        const VectorXD fitted = VectorXD(X_orig.eigen() * beta_hat.eigen());

        // 提取 x_names / y_name
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

            // 生成权重 v_i, 构造 y*_i = fitted_i + v_i · ε̂_i (H-014, H-016)
            VectorXD y_b(N);
            for (Size i = 0; i < N; ++i) {
                const Real v = detail::generate_wild_weight(rng, weight_dist_);
                y_b(i) = fitted(i) + v * resid(i);
            }

            // 估计 θ̂*_b (X 固定, 仅 y 变)
            try {
                // H-014: X 固定, 仅 y 重生成
                EconData data_b = make_cross_section(X_orig, y_b, x_names, y_name);
                EstimationResult r_b = estimator.estimate(data_b);

                // 验证结果有效
                if (r_b.coefficients.size() != k ||
                    !r_b.coefficients.eigen().allFinite()) {
                    ++n_failed;
                    continue;
                }

                samples.push_back(r_b.coefficients);
            } catch (const std::exception&) {
                // H-009: 失败计入 n_failed
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

        // 第 0 个系数的 CI (向后兼容)
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
        return std::string("WildBootstrap[") + to_string(weight_dist_) + "]";
    }

private:
    WildWeightDistribution weight_dist_;
};

}  // namespace econometrics
}  // namespace v1
}  // namespace cpphub
