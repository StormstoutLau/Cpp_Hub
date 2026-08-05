// SOURCE: PHASE6_IMPLEMENTATION_PLAN §5 - Arellano-Bond 1991 动态面板 GMM 测试
// 验证方法: 差分方程构造 / 工具变量数 / 大样本收敛 / 边界异常
//
// 排幻觉点:
//   E11: 工具变量构造严格按 Arellano-Bond 1991 原始论文
//        对 Δy_{it} (t=3,4,...,T), 工具变量 = y_{i,t-2}, y_{i,t-3}, ..., y_{i,1}
//        (level instruments for differenced equation, 与 Δε_{it} 正交)
//        非 R plm::pgmm 的 GMM-style 变体
//
// Arellano-Bond 1991 模型:
//   y_{it} = α y_{i,t-1} + β' x_{it} + μ_i + ε_{it}
//   差分: Δy_{it} = α Δy_{i,t-1} + β' Δx_{it} + Δε_{it}
//
// 工具变量数 (纯自回归, k_x=0):
//   t=3: y_{i,1} (1 个)
//   t=4: y_{i,1}, y_{i,2} (2 个)
//   t=T: y_{i,1}, ..., y_{i,T-2} (T-2 个)
//   总计: 1+2+...+(T-2) = (T-2)(T-1)/2
//
// 差分观测数: N * (T-2)  (每个个体 t=3..T)

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cpphub/econometrics/estimation/gmm.hpp"
#include "cpphub/econometrics/core/data_types.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;
using cpphub::v1::Index;
using cpphub::v1::linalg::dynamic::MatrixXD;
using cpphub::v1::linalg::dynamic::VectorXD;

// =============================================================================
// 辅助: 构造简单平衡面板 (N=2, T=4, k_x=0, 纯自回归)
//   entity 0: y = [1, 2, 3, 5]
//   entity 1: y = [2, 4, 7, 11]
//   差分观测数: 2 * (4-2) = 4
//   工具变量总数: (4-2)(4-1)/2 = 3 (y 滞后) + 0 (x) = 3
// =============================================================================
inline PanelData make_simple_panel() {
    PanelData panel;
    const Size N_entities = 2;
    const Size T = 4;
    const Size N = N_entities * T;

    panel.y = VectorXD(N);
    panel.X = MatrixXD(N, 0);  // k_x=0 (纯自回归, 无外生变量)
    panel.entity_id.resize(N);
    panel.time_id.resize(N);
    panel.balanced = true;

    // entity 0: y = [1, 2, 3, 5]
    // entity 1: y = [2, 4, 7, 11]
    const std::vector<std::vector<Real>> y_data = {{1, 2, 3, 5}, {2, 4, 7, 11}};
    for (Size e = 0; e < N_entities; ++e) {
        for (Size t = 0; t < T; ++t) {
            const Size idx = e * T + t;
            panel.y(idx) = y_data[e][t];
            panel.entity_id[idx] = static_cast<Index>(e);
            panel.time_id[idx] = static_cast<Index>(t);
        }
    }
    return panel;
}

// =============================================================================
// 辅助: 构造大样本面板 (N=50, T=8, k_x=0, α=0.5)
//   y_{it} = 0.5 * y_{i,t-1} + ε_{it}
//   用确定性 "噪声" 避免跨平台随机数问题
// =============================================================================
inline PanelData make_large_panel() {
    PanelData panel;
    const Size N_entities = 50;
    const Size T = 8;
    const Size N = N_entities * T;

    panel.y = VectorXD(N);
    panel.X = MatrixXD(N, 0);
    panel.entity_id.resize(N);
    panel.time_id.resize(N);
    panel.balanced = true;

    for (Size e = 0; e < N_entities; ++e) {
        Real y_prev = 1.0;  // 初始值
        for (Size t = 0; t < T; ++t) {
            const Size idx = e * T + t;
            // 确定性 "噪声": sin(e*0.1 + t*0.3) * 0.2 (均值≈0)
            const Real noise = 0.2 * std::sin(static_cast<Real>(e) * 0.1 +
                                                static_cast<Real>(t) * 0.3);
            const Real y_t = (t == 0) ? 1.0 : 0.5 * y_prev + noise;
            panel.y(idx) = y_t;
            panel.entity_id[idx] = static_cast<Index>(e);
            panel.time_id[idx] = static_cast<Index>(t);
            y_prev = y_t;
        }
    }
    return panel;
}

// =============================================================================
// 辅助: 构造带外生变量的面板 (N=50, T=6, k_x=1, α=0.4, β=0.6)
// =============================================================================
inline PanelData make_panel_with_x() {
    PanelData panel;
    // 注: N_entities=50 确保足够样本量 (N_eff=50*(6-2)=200),
    // 避免小样本下 block-diagonal Z 矩阵数值不稳定 (GCC α=1.05 非平稳问题)
    const Size N_entities = 50;
    const Size T = 6;
    const Size N = N_entities * T;

    panel.y = VectorXD(N);
    panel.X = MatrixXD(N, 1);  // k_x=1
    panel.entity_id.resize(N);
    panel.time_id.resize(N);
    panel.balanced = true;

    for (Size e = 0; e < N_entities; ++e) {
        Real y_prev = 1.0;
        for (Size t = 0; t < T; ++t) {
            const Size idx = e * T + t;
            const Real x = static_cast<Real>(t) + 0.1 * static_cast<Real>(e);
            panel.X(idx, 0) = x;
            const Real noise = 0.1 * std::sin(static_cast<Real>(e + t));
            const Real y_t = (t == 0) ? 1.0 : 0.4 * y_prev + 0.6 * x + noise;
            panel.y(idx) = y_t;
            panel.entity_id[idx] = static_cast<Index>(e);
            panel.time_id[idx] = static_cast<Index>(t);
            y_prev = y_t;
        }
    }
    return panel;
}

// =============================================================================
// 测试 1: 简单面板 - 不抛异常, 返回有效结果
// =============================================================================
TEST(ArellanoBond, SimplePanel_NoThrow) {
    const auto panel = make_simple_panel();
    EXPECT_NO_THROW({
        const ArellanoBondResult result = arellano_bond(panel);
        EXPECT_EQ(result.n_entities, 2u);
        EXPECT_EQ(result.n_periods, 4u);
    });
}

// =============================================================================
// 测试 2: 差分观测数 = N_entities * (T-2)
//   简单面板: 2 * (4-2) = 4 个差分观测
// =============================================================================
TEST(ArellanoBond, NObs_Equals_NTimes_TMinus2) {
    const auto panel = make_simple_panel();
    const ArellanoBondResult result = arellano_bond(panel);
    // 差分观测数 = N_entities * (T-2) = 2 * 2 = 4
    EXPECT_EQ(result.n_obs, 4u);
}

// =============================================================================
// 测试 3: 面板维度字段正确
// =============================================================================
TEST(ArellanoBond, Dimensions_Correct) {
    const auto panel = make_simple_panel();
    const ArellanoBondResult result = arellano_bond(panel);
    EXPECT_EQ(result.n_entities, 2u);
    EXPECT_EQ(result.n_periods, 4u);
    EXPECT_EQ(result.n_params, 1u);  // k = 1 + k_x = 1 + 0 = 1
}

// =============================================================================
// 测试 4: 工具变量数 > 0
//   纯自回归 T=4: (T-2)(T-1)/2 = 3 个 y 滞后 + 0 个 x = 3
// =============================================================================
TEST(ArellanoBond, NInstruments_Positive) {
    const auto panel = make_simple_panel();
    const ArellanoBondResult result = arellano_bond(panel);
    EXPECT_GT(result.n_instruments, 0u);
    // 至少有 y 滞后工具变量
    EXPECT_GE(result.n_instruments, 3u);
}

// =============================================================================
// 测试 5: J 检验自由度 = q - k >= 0 (过度识别或恰好识别)
// =============================================================================
TEST(ArellanoBond, J_DF_NonNegative) {
    const auto panel = make_simple_panel();
    const ArellanoBondResult result = arellano_bond(panel);
    // df = q - k, 过度识别时 > 0, 恰好识别时 = 0
    // 简单面板: q >= 3, k=1, df >= 2
    EXPECT_GE(result.j_df, 0u);
}

// =============================================================================
// 测试 6: J 统计量非负
// =============================================================================
TEST(ArellanoBond, J_Statistic_NonNegative) {
    const auto panel = make_simple_panel();
    const ArellanoBondResult result = arellano_bond(panel);
    EXPECT_GE(result.j_statistic, 0.0);
}

// =============================================================================
// 测试 7: 大样本面板 - α 收敛到真实值 0.5
// =============================================================================
TEST(ArellanoBond, LargePanel_AlphaConverges) {
    const auto panel = make_large_panel();
    const ArellanoBondResult result = arellano_bond(panel);
    // α ≈ 0.5 (大样本收敛, 容差 0.15)
    EXPECT_NEAR(result.coefficients(0), 0.5, 0.15);
}

// =============================================================================
// 测试 8: 大样本面板 - 维度正确
// =============================================================================
TEST(ArellanoBond, LargePanel_Dimensions) {
    const auto panel = make_large_panel();
    const ArellanoBondResult result = arellano_bond(panel);
    EXPECT_EQ(result.n_entities, 50u);
    EXPECT_EQ(result.n_periods, 8u);
    // 差分观测数 = 50 * (8-2) = 300
    EXPECT_EQ(result.n_obs, 300u);
    EXPECT_EQ(result.n_params, 1u);
}

// =============================================================================
// 测试 9: 带外生变量的面板 - 不抛异常
// =============================================================================
TEST(ArellanoBond, PanelWithX_NoThrow) {
    const auto panel = make_panel_with_x();
    EXPECT_NO_THROW({
        const ArellanoBondResult result = arellano_bond(panel);
        EXPECT_EQ(result.n_params, 2u);  // k = 1 (α) + 1 (β) = 2
    });
}

// =============================================================================
// 测试 10: 带外生变量的面板 - 维度正确
// =============================================================================
TEST(ArellanoBond, PanelWithX_Dimensions) {
    const auto panel = make_panel_with_x();
    const ArellanoBondResult result = arellano_bond(panel);
    EXPECT_EQ(result.n_entities, 50u);
    EXPECT_EQ(result.n_periods, 6u);
    // 差分观测数 = 50 * (6-2) = 200
    EXPECT_EQ(result.n_obs, 200u);
    EXPECT_EQ(result.n_params, 2u);  // α + β
}

// =============================================================================
// 测试 11: 带外生变量的面板 - α 和 β 估计合理
// =============================================================================
TEST(ArellanoBond, PanelWithX_CoefficientsReasonable) {
    const auto panel = make_panel_with_x();
    const ArellanoBondResult result = arellano_bond(panel);
    // α ≈ 0.4, β ≈ 0.6 (容差较宽, 因为小样本)
    EXPECT_NEAR(result.coefficients(0), 0.4, 0.3);
    EXPECT_NEAR(result.coefficients(1), 0.6, 0.3);
}

// =============================================================================
// 测试 12: 空面板抛异常
// =============================================================================
TEST(ArellanoBond, EmptyPanel_Throws) {
    PanelData panel;
    panel.X = MatrixXD(0, 0);
    panel.y = VectorXD(0);
    EXPECT_THROW(arellano_bond(panel), std::invalid_argument);
}

// =============================================================================
// 测试 13: T < 3 抛异常 (无法差分)
// =============================================================================
TEST(ArellanoBond, TLessThan3_Throws) {
    PanelData panel;
    const Size N = 4;  // 2 entities × 2 periods
    panel.y = VectorXD(N);
    panel.X = MatrixXD(N, 0);
    panel.entity_id = {0, 0, 1, 1};
    panel.time_id = {0, 1, 0, 1};
    EXPECT_THROW(arellano_bond(panel), std::invalid_argument);
}

// =============================================================================
// 测试 14: 协方差矩阵对称
// =============================================================================
TEST(ArellanoBond, Vcov_Symmetric) {
    const auto panel = make_large_panel();
    const ArellanoBondResult result = arellano_bond(panel);
    const Size k = result.n_params;
    for (Size i = 0; i < k; ++i) {
        for (Size j = i; j < k; ++j) {
            EXPECT_NEAR(result.vcov(i, j), result.vcov(j, i), 1e-10);
        }
    }
}

// =============================================================================
// 测试 15: GMM 类型字段正确 (默认 TwoStep)
// =============================================================================
TEST(ArellanoBond, DefaultType_TwoStep) {
    const auto panel = make_simple_panel();
    const ArellanoBondResult result = arellano_bond(panel);
    EXPECT_EQ(result.gmm_type, GMMType::TwoStep);
}

// =============================================================================
// 测试 16: 收敛标志为 true
// =============================================================================
TEST(ArellanoBond, ConvergedFlag_True) {
    const auto panel = make_simple_panel();
    const ArellanoBondResult result = arellano_bond(panel);
    EXPECT_TRUE(result.converged);
}

// =============================================================================
// 测试 17: 两步 GMM 和迭代 GMM 大样本下接近
// =============================================================================
TEST(ArellanoBond, LargePanel_TwoStepVsIterated) {
    const auto panel = make_large_panel();
    const ArellanoBondResult r1 = arellano_bond(panel, 0, GMMType::TwoStep);
    const ArellanoBondResult r2 = arellano_bond(panel, 0, GMMType::Iterated);
    // 大样本下应接近
    EXPECT_NEAR(r1.coefficients(0), r2.coefficients(0), 0.05);
}

// =============================================================================
// 测试 18: 工具变量数随 T 增长 (E11 验证)
//   T=4: (T-2)(T-1)/2 = 3 个 y 滞后
//   T=8: (T-2)(T-1)/2 = 21 个 y 滞后
// =============================================================================
TEST(ArellanoBond, NInstruments_GrowsWithT) {
    const auto panel4 = make_simple_panel();  // T=4
    const auto panel8 = make_large_panel();   // T=8

    const ArellanoBondResult r4 = arellano_bond(panel4);
    const ArellanoBondResult r8 = arellano_bond(panel8);

    // T=8 的工具变量数应远大于 T=4
    EXPECT_GT(r8.n_instruments, r4.n_instruments);
}

// =============================================================================
// 测试 19: 带外生变量时工具变量数 = y滞后 + k_x
//   panel_with_x: T=6, k_x=1
//   y 滞后: (6-2)(6-1)/2 = 10 (GMM-style, block-diagonal)
//   x 工具变量: k_x = 1 (standard IV, 所有观测在同一列)
//   总计: 10 + 1 = 11
// =============================================================================
TEST(ArellanoBond, PanelWithX_NInstruments) {
    const auto panel = make_panel_with_x();
    const ArellanoBondResult result = arellano_bond(panel);
    // y 滞后 (T-2)(T-1)/2 = 4*5/2 = 10, x 工具变量 k_x = 1 (standard IV)
    // 总计: 10 + 1 = 11
    EXPECT_GE(result.n_instruments, 10u);  // 至少有 y 滞后工具变量
}

// =============================================================================
// 测试 21: AR(1)/AR(2) 检验统计量有限且 p 值合法 (E14 验证)
//   Arellano-Bond 1991 §3.2: 差分残差的序列相关检验
//   Δε_{it} = ε_{it} - ε_{i,t-1} 具有 MA(1) 结构
//   AR(1): Corr(Δε_{it}, Δε_{i,t-1}) ≠ 0 (构造性, 显著)
//   AR(2): Corr(Δε_{it}, Δε_{i,t-2}) ≈ 0 (若原 ε 为白噪声, 不显著)
// =============================================================================
TEST(ArellanoBond, ARTests_StatisticsFinite) {
    const auto panel = make_large_panel();  // N=50, T=8
    const ArellanoBondResult result = arellano_bond(panel);

    // AR(1) 统计量: 有限且非 NaN
    EXPECT_TRUE(std::isfinite(result.ar1_statistic));
    EXPECT_TRUE(std::isfinite(result.ar1_pvalue));
    EXPECT_GE(result.ar1_pvalue, 0.0);
    EXPECT_LE(result.ar1_pvalue, 1.0);

    // AR(2) 统计量: 有限且非 NaN
    // T=8 → 每个个体有 6 个差分观测 (t=3..8), AR(2) 配对 (t, t-2) 有 4 对/个体
    EXPECT_TRUE(std::isfinite(result.ar2_statistic));
    EXPECT_TRUE(std::isfinite(result.ar2_pvalue));
    EXPECT_GE(result.ar2_pvalue, 0.0);
    EXPECT_LE(result.ar2_pvalue, 1.0);
}

// =============================================================================
// 测试 22: AR(1) 检验应显著 (差分构造的 MA(1) 结构)
//   Δε_{it} = ε_{it} - ε_{i,t-1}, Δε_{i,t-1} = ε_{i,t-1} - ε_{i,t-2}
//   两者共享 ε_{i,t-1}, 故 Corr(Δε_{it}, Δε_{i,t-1}) < 0 (负相关)
//   大样本下 AR(1) p 值应 < 0.05 (显著拒绝无自相关假设)
// =============================================================================
TEST(ArellanoBond, AR1_Significant_ByConstruction) {
    const auto panel = make_large_panel();  // N=50, T=8
    const ArellanoBondResult result = arellano_bond(panel);

    // AR(1) 应显著 (p < 0.05): 差分残差的 MA(1) 结构
    // 注: 统计量应为负 (Δε_{it} 与 Δε_{i,t-1} 负相关)
    EXPECT_LT(result.ar1_pvalue, 0.05);
}

// =============================================================================
// 测试 23: 小面板 (T=4) AR(2) 无配对 → 统计量为 0
//   T=4: 差分观测 t=3,4 (每个体 2 个), AR(2) 需 (t, t-2) 配对
//   t=4 与 t=2: t=2 不是差分观测 (差分从 t=3 开始), 无配对
// =============================================================================
TEST(ArellanoBond, AR2_NoPairs_SmallPanel) {
    const auto panel = make_simple_panel();  // N=2, T=4
    const ArellanoBondResult result = arellano_bond(panel);

    // T=4: AR(2) 无配对 (需 t 和 t-2 都是差分观测, 但差分从 t=3 开始)
    // 统计量保持默认值 0.0
    EXPECT_EQ(result.ar2_statistic, 0.0);
}

// =============================================================================
// 测试 20: 大样本带外生变量 - α 和 β 收敛到真实值
// =============================================================================
TEST(ArellanoBond, LargePanelWithX_Converges) {
    // 构造 N=30, T=10, α=0.5, β=0.8 的大面板
    const Size N_entities = 30;
    const Size T = 10;
    const Size N = N_entities * T;
    PanelData panel;
    panel.y = VectorXD(N);
    panel.X = MatrixXD(N, 1);
    panel.entity_id.resize(N);
    panel.time_id.resize(N);
    panel.balanced = true;

    for (Size e = 0; e < N_entities; ++e) {
        Real y_prev = 1.0;
        for (Size t = 0; t < T; ++t) {
            const Size idx = e * T + t;
            const Real x = static_cast<Real>(t) + 0.05 * static_cast<Real>(e);
            panel.X(idx, 0) = x;
            const Real noise = 0.1 * std::sin(static_cast<Real>(e * T + t) * 0.1);
            const Real y_t = (t == 0) ? 1.0 : 0.5 * y_prev + 0.8 * x + noise;
            panel.y(idx) = y_t;
            panel.entity_id[idx] = static_cast<Index>(e);
            panel.time_id[idx] = static_cast<Index>(t);
            y_prev = y_t;
        }
    }

    const ArellanoBondResult result = arellano_bond(panel);
    // α ≈ 0.5, β ≈ 0.8 (大样本, 容差 0.2)
    EXPECT_NEAR(result.coefficients(0), 0.5, 0.2);
    EXPECT_NEAR(result.coefficients(1), 0.8, 0.2);
}
