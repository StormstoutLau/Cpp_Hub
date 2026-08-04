// SOURCE: PHASE6_IMPLEMENTATION_PLAN §3.4 任务 1.11 - 面板数据工具测试
// 覆盖: detect_balance / extract_cluster_ids / reshape_y_to_wide / reshape_x_to_wide
// TDD: 手算期望值, 不依赖外部基准
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cpphub/econometrics/data/panel_data.hpp"

using namespace cpphub::v1::econometrics;
using cpphub::v1::Real;
using cpphub::v1::Size;
using cpphub::v1::Index;

// =============================================================================
// 辅助: 构造平衡面板 (N=6, E=3, T=2)
// entity=[0,0,1,1,2,2], time=[0,1,0,1,0,1]
// y=[10,11,20,21,30,31], X[:,0]=[1,2,3,4,5,6]
// =============================================================================
PanelData make_balanced_panel() {
    MatrixXD X(6, 1);
    X(0, 0) = 1.0; X(1, 0) = 2.0; X(2, 0) = 3.0;
    X(3, 0) = 4.0; X(4, 0) = 5.0; X(5, 0) = 6.0;
    VectorXD y(6);
    y(0) = 10.0; y(1) = 11.0; y(2) = 20.0;
    y(3) = 21.0; y(4) = 30.0; y(5) = 31.0;
    return make_panel(X, y, {0, 0, 1, 1, 2, 2}, {0, 1, 0, 1, 0, 1},
                       {"x"}, "y", true);
}

// =============================================================================
// 辅助: 构造非平衡面板 (entity 0 有 3 期, entity 1 有 2 期)
// entity=[0,0,0,1,1], time=[0,1,2,0,1]
// =============================================================================
PanelData make_unbalanced_panel() {
    MatrixXD X(5, 1);
    X(0, 0) = 1.0; X(1, 0) = 2.0; X(2, 0) = 3.0;
    X(3, 0) = 4.0; X(4, 0) = 5.0;
    VectorXD y(5);
    y(0) = 10.0; y(1) = 11.0; y(2) = 12.0; y(3) = 20.0; y(4) = 21.0;
    return make_panel(X, y, {0, 0, 0, 1, 1}, {0, 1, 2, 0, 1},
                       {"x"}, "y", false);
}

// =============================================================================
// 1. detect_balance 平衡面板 → is_balanced=true, E=3, T=2
// =============================================================================
TEST(PanelDataTest, DetectBalance_Balanced) {
    PanelData p = make_balanced_panel();
    PanelBalance b = detect_balance(p);
    EXPECT_TRUE(b.is_balanced);
    EXPECT_EQ(b.n_entities, 3);
    EXPECT_EQ(b.n_periods, 2);
}

// =============================================================================
// 2. detect_balance 非平衡面板 → is_balanced=false, n_periods=0
// =============================================================================
TEST(PanelDataTest, DetectBalance_Unbalanced) {
    PanelData p = make_unbalanced_panel();
    PanelBalance b = detect_balance(p);
    EXPECT_FALSE(b.is_balanced);
    EXPECT_EQ(b.n_entities, 2);
    EXPECT_EQ(b.n_periods, 0);
}

// =============================================================================
// 3. detect_balance 空面板 → 抛异常
// =============================================================================
TEST(PanelDataTest, DetectBalance_Empty_Throws) {
    PanelData p;
    EXPECT_THROW(detect_balance(p), std::invalid_argument);
}

// =============================================================================
// 4. detect_balance 单 entity 多时期 → is_balanced=true, E=1
// =============================================================================
TEST(PanelDataTest, DetectBalance_SingleEntity) {
    MatrixXD X(3, 1);
    X(0, 0) = 1.0; X(1, 0) = 2.0; X(2, 0) = 3.0;
    VectorXD y(3);
    y(0) = 10.0; y(1) = 11.0; y(2) = 12.0;
    PanelData p = make_panel(X, y, {0, 0, 0}, {0, 1, 2}, {"x"}, "y", true);
    PanelBalance b = detect_balance(p);
    EXPECT_TRUE(b.is_balanced);
    EXPECT_EQ(b.n_entities, 1);
    EXPECT_EQ(b.n_periods, 3);
}

// =============================================================================
// 5. reshape_y_to_wide 平衡面板 → 宽表 (3×2), 数值正确
// =============================================================================
TEST(PanelDataTest, ReshapeY_Balanced) {
    PanelData p = make_balanced_panel();
    MatrixXD wide = reshape_y_to_wide(p);
    EXPECT_EQ(wide.rows(), 3u);
    EXPECT_EQ(wide.cols(), 2u);
    // entity 0: y=[10, 11] (time 0, 1)
    EXPECT_NEAR(wide(0, 0), 10.0, 1e-12);
    EXPECT_NEAR(wide(0, 1), 11.0, 1e-12);
    // entity 1: y=[20, 21]
    EXPECT_NEAR(wide(1, 0), 20.0, 1e-12);
    EXPECT_NEAR(wide(1, 1), 21.0, 1e-12);
    // entity 2: y=[30, 31]
    EXPECT_NEAR(wide(2, 0), 30.0, 1e-12);
    EXPECT_NEAR(wide(2, 1), 31.0, 1e-12);
}

// =============================================================================
// 6. reshape_y_to_wide 非平衡面板 → 抛异常 (当前实现不支持)
// =============================================================================
TEST(PanelDataTest, ReshapeY_Unbalanced_Throws) {
    PanelData p = make_unbalanced_panel();
    EXPECT_THROW(reshape_y_to_wide(p), std::invalid_argument);
}

// =============================================================================
// 7. reshape_y_to_wide entity 顺序按首次出现 (entity=[2,2,1,1,0,0])
//    行序应为 [2, 1, 0]
// =============================================================================
TEST(PanelDataTest, ReshapeY_EntityOrder_FirstOccurrence) {
    MatrixXD X(6, 1);
    for (int i = 0; i < 6; ++i) X(i, 0) = static_cast<Real>(i + 1);
    VectorXD y(6);
    y(0) = 10.0; y(1) = 11.0; y(2) = 20.0;
    y(3) = 21.0; y(4) = 30.0; y(5) = 31.0;
    // entity 顺序: 2, 1, 0 (不是 0, 1, 2)
    PanelData p = make_panel(X, y, {2, 2, 1, 1, 0, 0}, {0, 1, 0, 1, 0, 1},
                               {"x"}, "y", true);
    MatrixXD wide = reshape_y_to_wide(p);
    EXPECT_EQ(wide.rows(), 3u);
    EXPECT_EQ(wide.cols(), 2u);
    // 第一行应为 entity 2 的数据 (y=[10, 11])
    EXPECT_NEAR(wide(0, 0), 10.0, 1e-12);
    EXPECT_NEAR(wide(0, 1), 11.0, 1e-12);
    // 第二行应为 entity 1 的数据 (y=[20, 21])
    EXPECT_NEAR(wide(1, 0), 20.0, 1e-12);
    EXPECT_NEAR(wide(1, 1), 21.0, 1e-12);
    // 第三行应为 entity 0 的数据 (y=[30, 31])
    EXPECT_NEAR(wide(2, 0), 30.0, 1e-12);
    EXPECT_NEAR(wide(2, 1), 31.0, 1e-12);
}

// =============================================================================
// 8. extract_cluster_ids entity 类型, 原始 ID [3,3,5,5,7,7] → 紧凑 [0,0,1,1,2,2]
// =============================================================================
TEST(PanelDataTest, ExtractClusterIds_Entity) {
    MatrixXD X(6, 1);
    VectorXD y(6);
    PanelData p = make_panel(X, y, {3, 3, 5, 5, 7, 7}, {0, 1, 0, 1, 0, 1},
                               {"x"}, "y", true);
    ClusterIdResult r = extract_cluster_ids(p, "entity");
    EXPECT_EQ(r.n_clusters, 3);
    ASSERT_EQ(r.ids.size(), 6u);
    EXPECT_EQ(r.ids[0], 0); EXPECT_EQ(r.ids[1], 0);
    EXPECT_EQ(r.ids[2], 1); EXPECT_EQ(r.ids[3], 1);
    EXPECT_EQ(r.ids[4], 2); EXPECT_EQ(r.ids[5], 2);
    // original 保留首次出现顺序
    ASSERT_EQ(r.original.size(), 3u);
    EXPECT_EQ(r.original[0], 3);
    EXPECT_EQ(r.original[1], 5);
    EXPECT_EQ(r.original[2], 7);
}

// =============================================================================
// 9. extract_cluster_ids time 类型, 原始 ID [10,20,10,20] → 紧凑 [0,1,0,1]
// =============================================================================
TEST(PanelDataTest, ExtractClusterIds_Time) {
    MatrixXD X(4, 1);
    VectorXD y(4);
    PanelData p = make_panel(X, y, {0, 0, 1, 1}, {10, 20, 10, 20},
                               {"x"}, "y", true);
    ClusterIdResult r = extract_cluster_ids(p, "time");
    EXPECT_EQ(r.n_clusters, 2);
    ASSERT_EQ(r.ids.size(), 4u);
    EXPECT_EQ(r.ids[0], 0); EXPECT_EQ(r.ids[1], 1);
    EXPECT_EQ(r.ids[2], 0); EXPECT_EQ(r.ids[3], 1);
    EXPECT_EQ(r.original[0], 10);
    EXPECT_EQ(r.original[1], 20);
}

// =============================================================================
// 10. extract_cluster_ids 空输入 → n_clusters=0
// =============================================================================
TEST(PanelDataTest, ExtractClusterIds_Empty) {
    PanelData p;
    ClusterIdResult r = extract_cluster_ids(p, "entity");
    EXPECT_EQ(r.n_clusters, 0);
    EXPECT_TRUE(r.ids.empty());
}

// =============================================================================
// 11. extract_cluster_ids 单聚类 (所有观测同一 entity) → n_clusters=1
// =============================================================================
TEST(PanelDataTest, ExtractClusterIds_SingleCluster) {
    MatrixXD X(3, 1);
    VectorXD y(3);
    PanelData p = make_panel(X, y, {7, 7, 7}, {0, 1, 2}, {"x"}, "y", true);
    ClusterIdResult r = extract_cluster_ids(p, "entity");
    EXPECT_EQ(r.n_clusters, 1);
    ASSERT_EQ(r.ids.size(), 3u);
    EXPECT_EQ(r.ids[0], 0);
    EXPECT_EQ(r.ids[1], 0);
    EXPECT_EQ(r.ids[2], 0);
}

// =============================================================================
// 12. reshape_x_to_wide 平衡面板, X 第 0 列 → 宽表 (3×2)
// =============================================================================
TEST(PanelDataTest, ReshapeX_Balanced_Col0) {
    PanelData p = make_balanced_panel();
    MatrixXD wide = reshape_x_to_wide(p, 0);
    EXPECT_EQ(wide.rows(), 3u);
    EXPECT_EQ(wide.cols(), 2u);
    // entity 0: X[:,0]=[1, 2] (time 0, 1)
    EXPECT_NEAR(wide(0, 0), 1.0, 1e-12);
    EXPECT_NEAR(wide(0, 1), 2.0, 1e-12);
    // entity 2: X[:,0]=[5, 6]
    EXPECT_NEAR(wide(2, 0), 5.0, 1e-12);
    EXPECT_NEAR(wide(2, 1), 6.0, 1e-12);
}

// =============================================================================
// 13. extract_cluster_ids 无效 type → 抛异常
// =============================================================================
TEST(PanelDataTest, ExtractClusterIds_InvalidType_Throws) {
    PanelData p = make_balanced_panel();
    EXPECT_THROW(extract_cluster_ids(p, "invalid"), std::invalid_argument);
}

// =============================================================================
// 14. reshape_x_to_wide col_index 越界 → 抛异常
// =============================================================================
TEST(PanelDataTest, ReshapeX_ColIndexOutOfRange_Throws) {
    PanelData p = make_balanced_panel();
    EXPECT_THROW(reshape_x_to_wide(p, 1), std::invalid_argument);  // X 只有 1 列
    EXPECT_THROW(reshape_x_to_wide(p, -1), std::invalid_argument);
}

// =============================================================================
// 15. detect_balance entity 重复 time (同 entity 同时段 2 个观测) → 非平衡
// 排幻觉点 P1: 不是简单 N == E*T, 必须 time_id 集合相同
// =============================================================================
TEST(PanelDataTest, DetectBalance_DuplicateTime_Unbalanced) {
    // entity=[0,0,0,1,1,1], time=[0,1,1,0,1,2]
    // entity 0: time {0, 1} (但 time=1 重复, set={0,1}, size=2)
    // entity 1: time {0, 1, 2}, size=3
    // 集合大小不同 → 非平衡
    MatrixXD X(6, 1);
    VectorXD y(6);
    PanelData p = make_panel(X, y, {0, 0, 0, 1, 1, 1}, {0, 1, 1, 0, 1, 2},
                               {"x"}, "y", false);
    PanelBalance b = detect_balance(p);
    EXPECT_FALSE(b.is_balanced);
    EXPECT_EQ(b.n_entities, 2);
    EXPECT_EQ(b.n_periods, 0);
}
