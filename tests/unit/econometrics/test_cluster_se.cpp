// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.4 Day 7-8 任务 1.10 - 聚类稳健协方差 (排幻觉点 E6)
// 教材锚点:
//   - Liang-Zeger (1986) Biometrika (原始聚类 SE)
//   - Cameron-Gelbach-Miller (2011) 双向聚类
//   - Cameron-Miller (2015) "A Practitioner's Guide to Cluster-Robust Inference"
//
// 排幻觉点 (R sandwich::vcovCL 对照):
//   E6: 小样本调整 G/(G-1)·(N-1)/(N-K), 非 1/(G-1) 或 G/(G-1)
//       R vcovCL L80-100 实测: type="HC1" 默认含此调整, 分母 (G-1)·(N-K)
//   双向聚类: V_twoway = V(g1) + V(g2) - V(g1∩g2), 每个分量独立做小样本调整
//
// 手算数据集 1 (N=4, K=1, 无截距, G=2):
//   X=[1;2;3;4], y=[1;2;4;5], cluster_id=[0,0,1,1]
//   X'X=30, (X'X)^{-1}=1/30, β=37/30
//   u=[-7/30; -7/15; 3/10; 1/15]
//   聚类 0: X_0=[1;2], u_0=[-7/30;-7/15], v_0=X_0'u_0=-7/6, meat_0=49/36
//   聚类 1: X_1=[3;4], u_1=[3/10;1/15],  v_1=X_1'u_1=7/6,  meat_1=49/36
//   meat = 49/36 + 49/36 = 49/18
//   V_unscaled = (1/30)·(49/18)·(1/30) = 49/16200
//   小样本调整 = G/(G-1)·(N-1)/(N-K) = 2/1·3/3 = 2
//   V_cluster = 2·49/16200 = 49/8100
//
// 手算数据集 2 (N=4, K=2, 含截距, G=2):
//   X=[[1,1];[1,2];[1,3];[1,4]], y=[2;3;5;7], cluster_id=[0,0,1,1]
//   (X'X)^{-1}=[[3/2,-1/2];[-1/2,1/5]], β=[0;1.7], u=[0.3;-0.4;-0.1;0.2]
//   聚类 0: v_0=[-0.1;-0.5], meat_0=[[0.01,0.05];[0.05,0.25]]
//   聚类 1: v_1=[0.1;0.5],  meat_1=[[0.01,0.05];[0.05,0.25]]
//   meat = [[0.02,0.10];[0.10,0.50]]
//   V_unscaled = A·meat·A = [[1/50,-1/100];[-1/100,1/200]]
//   小样本调整 = 2/1·3/2 = 3
//   V_cluster = [[3/50,-3/100];[-3/100,3/200]]
//
// 手算数据集 3 (双向聚类, N=4, K=1, G1=2, G2=2):
//   X=[1;2;3;4], y=[1;2;4;5]
//   g1=[0,0,1,1] (entity), g2=[0,1,0,1] (time)
//   V(g1) = 49/8100 (同数据集 1)
//   V(g2): 聚类 0=[obs0,obs2], 聚类 1=[obs1,obs3]
//     v_0 = 1·(-7/30)+3·(3/10) = 2/3,  meat_0 = 4/9
//     v_1 = 2·(-7/15)+4·(1/15) = -2/3, meat_1 = 4/9
//     meat = 8/9, V(g2)_unscaled = 8/8100, scale=2, V(g2) = 16/8100
//   V(g1∩g2): 每观测一聚类, G12=4
//     meat = Σ x_i²u_i² = 49/900 + 196/225 + 81/100 + 16/225 = 271/150
//     V(g1∩g2)_unscaled = 271/135000
//     scale = 4/3·3/3 = 4/3
//     V(g1∩g2) = 4/3·271/135000 = 271/101250
//   V_twoway = 49/8100 + 16/8100 - 271/101250
//            = 1225/202500 + 400/202500 - 542/202500
//            = 1083/202500 = 361/67500

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cpphub/econometrics/inference/cluster_vcov.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;
using cpphub::v1::Index;

// =============================================================================
// 辅助: 数据集 1 (K=1, 单向聚类)
// =============================================================================
struct K1ClusterDataset {
    MatrixXD X;
    VectorXD u;
    MatrixXD XtX_inv;
    std::vector<Index> cluster_id;
    K1ClusterDataset() : X(4, 1), u(4), XtX_inv(1, 1), cluster_id{0, 0, 1, 1} {
        X(0, 0) = 1.0; X(1, 0) = 2.0; X(2, 0) = 3.0; X(3, 0) = 4.0;
        u(0) = -7.0 / 30.0; u(1) = -7.0 / 15.0;
        u(2) = 3.0 / 10.0;  u(3) = 1.0 / 15.0;
        XtX_inv(0, 0) = 1.0 / 30.0;
    }
};

// =============================================================================
// 辅助: 数据集 2 (K=2 含截距, 单向聚类)
// =============================================================================
struct K2ClusterDataset {
    MatrixXD X;
    VectorXD u;
    MatrixXD XtX_inv;
    std::vector<Index> cluster_id;
    K2ClusterDataset() : X(4, 2), u(4), XtX_inv(2, 2), cluster_id{0, 0, 1, 1} {
        X(0, 0) = 1.0; X(0, 1) = 1.0;
        X(1, 0) = 1.0; X(1, 1) = 2.0;
        X(2, 0) = 1.0; X(2, 1) = 3.0;
        X(3, 0) = 1.0; X(3, 1) = 4.0;
        u(0) = 0.3;  u(1) = -0.4; u(2) = -0.1; u(3) = 0.2;
        XtX_inv(0, 0) = 3.0 / 2.0;  XtX_inv(0, 1) = -1.0 / 2.0;
        XtX_inv(1, 0) = -1.0 / 2.0; XtX_inv(1, 1) = 1.0 / 5.0;
    }
};

// =============================================================================
// 辅助: 数据集 3 (K=1, 双向聚类)
// =============================================================================
struct K1TwoWayDataset {
    MatrixXD X;
    VectorXD u;
    MatrixXD XtX_inv;
    std::vector<Index> g1;  // entity
    std::vector<Index> g2;  // time
    K1TwoWayDataset() : X(4, 1), u(4), XtX_inv(1, 1),
                        g1{0, 0, 1, 1}, g2{0, 1, 0, 1} {
        X(0, 0) = 1.0; X(1, 0) = 2.0; X(2, 0) = 3.0; X(3, 0) = 4.0;
        u(0) = -7.0 / 30.0; u(1) = -7.0 / 15.0;
        u(2) = 3.0 / 10.0;  u(3) = 1.0 / 15.0;
        XtX_inv(0, 0) = 1.0 / 30.0;
    }
};

// =============================================================================
// 1. 单向聚类 K=1 精确值 (V = 49/8100)
// =============================================================================
TEST(ClusterSETest, OneWay_K1_ExactValue) {
    K1ClusterDataset d;
    MatrixXD V = compute_cluster_vcov(d.X, d.u, d.XtX_inv, d.cluster_id);
    EXPECT_NEAR(V(0, 0), 49.0 / 8100.0, 1e-12);
}

// =============================================================================
// 2. 单向聚类 K=2 含截距 精确值 (V = [[3/50, -3/100];[-3/100, 3/200]])
// =============================================================================
TEST(ClusterSETest, OneWay_K2_WithIntercept_ExactValue) {
    K2ClusterDataset d;
    MatrixXD V = compute_cluster_vcov(d.X, d.u, d.XtX_inv, d.cluster_id);
    EXPECT_NEAR(V(0, 0), 3.0 / 50.0, 1e-12);   // 0.06
    EXPECT_NEAR(V(0, 1), -3.0 / 100.0, 1e-12); // -0.03
    EXPECT_NEAR(V(1, 0), -3.0 / 100.0, 1e-12); // -0.03
    EXPECT_NEAR(V(1, 1), 3.0 / 200.0, 1e-12);  // 0.015
}

// =============================================================================
// 3. 单向聚类 V 对称 (K=2)
// =============================================================================
TEST(ClusterSETest, OneWay_K2_Symmetric) {
    K2ClusterDataset d;
    MatrixXD V = compute_cluster_vcov(d.X, d.u, d.XtX_inv, d.cluster_id);
    EXPECT_NEAR(V(0, 1), V(1, 0), 1e-12);
}

// =============================================================================
// 4. 排幻觉点 E6: 小样本调整 G/(G-1)·(N-1)/(N-K)
//    验证: V_cluster = scale · V_unscaled, scale = 2·1 = 2 (K=1 数据集)
//    若错误使用 1/(G-1) 调整, V 会小 4 倍; 若不调整, V 会小 2 倍
// =============================================================================
TEST(ClusterSETest, SmallSampleAdjustment_E6_GOverGMinus1) {
    K1ClusterDataset d;
    MatrixXD V = compute_cluster_vcov(d.X, d.u, d.XtX_inv, d.cluster_id);
    // V_unscaled = 49/16200, scale = 2, V = 49/8100
    Real V_unscaled = 49.0 / 16200.0;
    Real expected_scale = 2.0;  // G/(G-1)·(N-1)/(N-K) = 2/1·3/3 = 2
    EXPECT_NEAR(V(0, 0), expected_scale * V_unscaled, 1e-12);
    // 排幻觉: 不应是 1/(G-1) = 1.0 (那样 V = V_unscaled, 错)
    EXPECT_NE(V(0, 0), V_unscaled);
    // 排幻觉: 不应是 G/(G-1) = 2.0 不含 (N-1)/(N-K) (K=1 时恰好相同, 用 K=2 验证)
}

// =============================================================================
// 5. 排幻觉点 E6: K=2 时小样本调整 = G/(G-1)·(N-1)/(N-K) = 2/1·3/2 = 3
//    若错误使用 G/(G-1) = 2, V 会小 1.5 倍
// =============================================================================
TEST(ClusterSETest, SmallSampleAdjustment_E6_K2_WithNMinusK) {
    K2ClusterDataset d;
    MatrixXD V = compute_cluster_vcov(d.X, d.u, d.XtX_inv, d.cluster_id);
    // V_unscaled[0][0] = 1/50 = 0.02
    // 正确 scale = 3, V = 3/50 = 0.06
    // 错误 scale = 2 (漏 (N-1)/(N-K)), V = 2/50 = 0.04
    Real V_unscaled_00 = 1.0 / 50.0;
    Real correct_scale = 3.0;  // 2/1 · 3/2
    Real wrong_scale = 2.0;    // 仅 G/(G-1), 漏 (N-1)/(N-K)
    EXPECT_NEAR(V(0, 0), correct_scale * V_unscaled_00, 1e-12);
    EXPECT_NE(V(0, 0), wrong_scale * V_unscaled_00);
}

// =============================================================================
// 6. 双向聚类 K=1 精确值 (V_twoway = 361/67500)
//    V_twoway = V(g1) + V(g2) - V(g1∩g2) = 49/8100 + 16/8100 - 271/101250
// =============================================================================
TEST(ClusterSETest, TwoWay_K1_ExactValue) {
    K1TwoWayDataset d;
    MatrixXD V = compute_cluster_vcov(d.X, d.u, d.XtX_inv, d.g1, true, d.g2);
    EXPECT_NEAR(V(0, 0), 361.0 / 67500.0, 1e-12);
}

// =============================================================================
// 7. 双向聚类 = V(g1) + V(g2) - V(g1∩g2) 数值验证
// =============================================================================
TEST(ClusterSETest, TwoWay_EqualsSumMinusIntersect) {
    K1TwoWayDataset d;
    MatrixXD V_twoway = compute_cluster_vcov(d.X, d.u, d.XtX_inv, d.g1, true, d.g2);
    MatrixXD V_g1 = compute_cluster_vcov(d.X, d.u, d.XtX_inv, d.g1);

    // 构造 g1∩g2 合并聚类 ID
    std::vector<Index> g12 = merge_twoway_cluster_ids(d.g1, d.g2);
    MatrixXD V_g12 = compute_cluster_vcov(d.X, d.u, d.XtX_inv, g12);

    MatrixXD V_g2 = compute_cluster_vcov(d.X, d.u, d.XtX_inv, d.g2);

    // V_twoway = V(g1) + V(g2) - V(g1∩g2)
    Real expected = V_g1(0, 0) + V_g2(0, 0) - V_g12(0, 0);
    EXPECT_NEAR(V_twoway(0, 0), expected, 1e-12);
    EXPECT_NEAR(V_twoway(0, 0), 361.0 / 67500.0, 1e-12);
}

// =============================================================================
// 8. 边界: G=1 抛 std::invalid_argument (分母 G-1=0)
// =============================================================================
TEST(ClusterSETest, SingleCluster_Throws) {
    K1ClusterDataset d;
    std::vector<Index> single_cluster{0, 0, 0, 0};  // G=1
    EXPECT_THROW(compute_cluster_vcov(d.X, d.u, d.XtX_inv, single_cluster),
                 std::invalid_argument);
}

// =============================================================================
// 9. 边界: 聚类内观测数不等 (G=2, 但聚类 0 有 3 个, 聚类 1 有 1 个)
// =============================================================================
TEST(ClusterSETest, UnbalancedClusters_OK) {
    K1ClusterDataset d;
    std::vector<Index> unbalanced{0, 0, 0, 1};  // 聚类 0: 3 obs, 聚类 1: 1 obs
    // 应正常运行, 不抛异常
    MatrixXD V = compute_cluster_vcov(d.X, d.u, d.XtX_inv, unbalanced);
    EXPECT_GT(V(0, 0), 0.0);
    // V 应为正数 (方差非负)
}

// =============================================================================
// 10. 维度校验: cluster_id 长度与 N 不匹配抛异常
// =============================================================================
TEST(ClusterSETest, DimensionMismatch_Throws) {
    K1ClusterDataset d;
    std::vector<Index> short_id{0, 0, 1};  // 长度 3, 但 N=4
    EXPECT_THROW(compute_cluster_vcov(d.X, d.u, d.XtX_inv, short_id),
                 std::invalid_argument);
}

// =============================================================================
// 11. 双向聚类 g1∩g2 合并 ID 函数 (独立测试)
// =============================================================================
TEST(ClusterSETest, MergeTwoWayClusterIds_Correct) {
    std::vector<Index> g1{0, 0, 1, 1};
    std::vector<Index> g2{0, 1, 0, 1};
    std::vector<Index> merged = merge_twoway_cluster_ids(g1, g2);
    // 应得到 4 个唯一组合: (0,0)→0, (0,1)→1, (1,0)→2, (1,1)→3
    ASSERT_EQ(merged.size(), 4u);
    EXPECT_EQ(merged[0], 0);  // (0,0)
    EXPECT_EQ(merged[1], 1);  // (0,1)
    EXPECT_EQ(merged[2], 2);  // (1,0)
    EXPECT_EQ(merged[3], 3);  // (1,1)
}

// =============================================================================
// 12. 输出维度正确 (k × k)
// =============================================================================
TEST(ClusterSETest, ResultDimension_Correct) {
    K1ClusterDataset d1;
    MatrixXD V1 = compute_cluster_vcov(d1.X, d1.u, d1.XtX_inv, d1.cluster_id);
    EXPECT_EQ(V1.rows(), 1u);
    EXPECT_EQ(V1.cols(), 1u);

    K2ClusterDataset d2;
    MatrixXD V2 = compute_cluster_vcov(d2.X, d2.u, d2.XtX_inv, d2.cluster_id);
    EXPECT_EQ(V2.rows(), 2u);
    EXPECT_EQ(V2.cols(), 2u);
}

// =============================================================================
// 13. twoway=false 时 cluster_id2 被忽略 (单向聚类)
// =============================================================================
TEST(ClusterSETest, TwoWayFalse_IgnoresClusterId2) {
    K1ClusterDataset d;
    std::vector<Index> dummy_g2{5, 5, 5, 5};  // 不应被使用
    MatrixXD V1 = compute_cluster_vcov(d.X, d.u, d.XtX_inv, d.cluster_id, false, dummy_g2);
    MatrixXD V2 = compute_cluster_vcov(d.X, d.u, d.XtX_inv, d.cluster_id);
    EXPECT_NEAR(V1(0, 0), V2(0, 0), 1e-12);
}

// =============================================================================
// 14. 双向聚类 K=2 对称
// =============================================================================
TEST(ClusterSETest, TwoWay_K2_Symmetric) {
    // 复用 K2 数据集, 但需要扩展 g1/g2 到 N=4
    K2ClusterDataset d;
    std::vector<Index> g1{0, 0, 1, 1};
    std::vector<Index> g2{0, 1, 0, 1};
    MatrixXD V = compute_cluster_vcov(d.X, d.u, d.XtX_inv, g1, true, g2);
    EXPECT_NEAR(V(0, 1), V(1, 0), 1e-12);
}

// =============================================================================
// 15. 双向聚类 g1=g2 (完全重合) 时 V_twoway = V(g1) + V(g1) - V(g1) = V(g1)
// =============================================================================
TEST(ClusterSETest, TwoWay_G1EqualsG2_EqualsOneWay) {
    K1ClusterDataset d;
    std::vector<Index> g1{0, 0, 1, 1};
    std::vector<Index> g2{0, 0, 1, 1};  // 与 g1 相同
    MatrixXD V_twoway = compute_cluster_vcov(d.X, d.u, d.XtX_inv, g1, true, g2);
    MatrixXD V_oneway = compute_cluster_vcov(d.X, d.u, d.XtX_inv, g1);
    // V_twoway = V(g1) + V(g1) - V(g1) = V(g1)
    EXPECT_NEAR(V_twoway(0, 0), V_oneway(0, 0), 1e-12);
}
