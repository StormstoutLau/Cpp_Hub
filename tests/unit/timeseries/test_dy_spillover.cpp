// =============================================================================
// test_dy_spillover.cpp - DY 溢出指数 TCI/TO/FROM/NET + 滚动 (Phase 7C M2, 12 用例)
//
// 主基准: R Spillover 0.1.1 G.spillover/g.fevd/net/roll.spillover (1e-8)
// 幻觉点: V10 (H 可配敏感性) / §13-a (window 必填 > 2·K, 无默认)
// =============================================================================

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/timeseries/var/dy_spillover.hpp"
#include "cpphub/timeseries/var/fevd.hpp"
#include "cpphub/timeseries/var/multivariate_data.hpp"
#include "cpphub/timeseries/var/var_model.hpp"
#include "var_baseline.inc"

namespace vt = cpphub::v1::timeseries::var;
namespace vb = cpphub::v1::timeseries::var_baseline;

namespace {

vt::MultivariateTSData fixture_data() {
    vt::MultivariateTSData d;
    d.columns = {std::vector<cpphub::v1::Real>(std::begin(vb::Y1), std::end(vb::Y1)),
                 std::vector<cpphub::v1::Real>(std::begin(vb::Y2), std::end(vb::Y2)),
                 std::vector<cpphub::v1::Real>(std::begin(vb::Y3), std::end(vb::Y3))};
    d.names = {"y1", "y2", "y3"};
    return d;
}

}  // namespace

// 1. θ̃·100/K 表 vs G.spillover standardized (1e-8); R dump 列主序 → [3j+i]
TEST(DySpillover, TableVsRSpillover) {
    auto res = vt::dy_spillover_static(fixture_data(), 10, "c", 2);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(res.fevd(i, j) * 100.0 / 3.0, vb::SP_TABLE[3 * j + i], 1e-8)
                << i << j;
}

// 2. TCI vs G.spillover 表尾 (1e-8)
TEST(DySpillover, TCIvsSpillover) {
    auto res = vt::dy_spillover_static(fixture_data(), 10, "c", 2);
    EXPECT_NEAR(res.tci, vb::SP_TCI[0], 1e-8);
    EXPECT_GT(res.tci, 0.0);
    EXPECT_LT(res.tci, 100.0);
}

// 3. TO_j vs G.spillover (1e-8)
TEST(DySpillover, TOvsSpillover) {
    auto res = vt::dy_spillover_static(fixture_data(), 10, "c", 2);
    ASSERT_EQ(res.to_spillover.size(), 3u);
    for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(res.to_spillover[j], vb::SP_TO[j], 1e-8) << j;
}

// 4. FROM_j vs G.spillover (1e-8)
TEST(DySpillover, FROMvsSpillover) {
    auto res = vt::dy_spillover_static(fixture_data(), 10, "c", 2);
    for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(res.from_spillover[j], vb::SP_FROM[j], 1e-8) << j;
}

// 5. NET_j vs Spillover::net (1e-8)
TEST(DySpillover, NETvsSpillover) {
    auto res = vt::dy_spillover_static(fixture_data(), 10, "c", 2);
    for (int j = 0; j < 3; ++j) {
        EXPECT_NEAR(res.net_spillover[j], vb::SP_NET[j], 1e-8) << j;
        // 恒等式: NET = TO − FROM
        EXPECT_NEAR(res.net_spillover[j],
                    res.to_spillover[j] - res.from_spillover[j], 1e-12);
    }
    // 净方向: y1 发送方 (TRUE), y2 接收方
    EXPECT_GT(res.net_spillover[0], 0.0);
    EXPECT_LT(res.net_spillover[1], 0.0);
}

// 6. V10: H=10 vs H=50 敏感性 (均对照 R; 差异存在但小)
TEST(DySpillover, HSensitivity) {
    auto r10 = vt::dy_spillover_static(fixture_data(), 10, "c", 2);
    auto r50 = vt::dy_spillover_static(fixture_data(), 50, "c", 2);
    // H=50 主锚
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(r50.fevd(i, j), vb::SP_GFEVD_H50[3 * i + j], 1e-8) << i << j;
    EXPECT_NEAR(r50.tci, vb::SP_TCI_H50[0], 1e-8);
    // 敏感性: TCI 差异小 (收敛中) 但非零
    const double diff = std::abs(r50.tci - r10.tci);
    EXPECT_LT(diff, 0.1);
    EXPECT_GT(diff, 0.0);
}

// 7. 滚动 TCI 路径 vs R roll.spillover (w=150, 101 窗口)
//    R 侧 roll.spillover 未传 p ⇒ vars::VAR 默认 p=1 (trend const),
//    align="right" 窗口集与 C++ start=0..T−w 相同 → C++ 用 lag=1 对照
TEST(DySpillover, RollingTCIPathVsR) {
    auto res = vt::dy_spillover(fixture_data(), 150, 10, "c", 1);
    ASSERT_EQ(res.tci_path.size(), 101u);
    ASSERT_EQ(vb::SP_ROLL_TCI.size(), 101u);
    double max_err = 0.0;
    for (int w = 0; w < 101; ++w) {
        const double err = std::abs(res.tci_path[w] - vb::SP_ROLL_TCI[w]);
        max_err = std::max(max_err, err);
    }
    // R rollapply 每窗口 zoo 子样本 → VAR(p=1) → G.spillover; C++ 同路径
    EXPECT_LT(max_err, 1e-6) << "max_err=" << max_err;
}

// 8. §13-a: window 必填校验 (≤2K 拒绝; > T 拒绝; 可估计边界)
TEST(DySpillover, WindowValidation) {
    EXPECT_THROW(vt::dy_spillover(fixture_data(), 6, 10, "c", 2),
                 std::invalid_argument);   // 2K = 6 不允许
    EXPECT_THROW(vt::dy_spillover(fixture_data(), 251, 10, "c", 2),
                 std::invalid_argument);   // > T
    // 可估计性: Σ_mle 满秩需 w ≥ 12; 再加首窗平稳性 (V12) 为数据依赖,
    // 此夹具取 w=20 稳定通过 (w=12 首窗伴随特征值 ≥1 触发 V12 拦截)
    EXPECT_NO_THROW(vt::dy_spillover(fixture_data(), 20, 10, "c", 2));
    EXPECT_THROW(vt::dy_spillover(fixture_data(), 11, 10, "c", 2),
                 std::invalid_argument);
}

// 9. 滚动路径结构 (有限值 + 数量 = T − window + 1)
TEST(DySpillover, RollingPathStructure) {
    auto res = vt::dy_spillover(fixture_data(), 100, 10, "c", 2);
    EXPECT_EQ(res.tci_path.size(), 151u);
    EXPECT_EQ(res.net_path.size(), 151u);
    for (const auto& tci : res.tci_path) {
        EXPECT_TRUE(std::isfinite(tci));
        EXPECT_GT(tci, -1.0);
        EXPECT_LT(tci, 101.0);
    }
    for (const auto& net : res.net_path) {
        ASSERT_EQ(net.size(), 3u);
        for (double v : net) EXPECT_TRUE(std::isfinite(v));
    }
    EXPECT_EQ(res.window, 100u);
    EXPECT_EQ(res.horizon, 10u);
    // 末窗口 = 标量输出
    EXPECT_NEAR(res.tci, res.tci_path.back(), 0.0);
}

// 10. 末窗口一致性: 滚动末窗口 == 手动子样本 static (T=150 子样本)
TEST(DySpillover, RollingLastWindowConsistency) {
    auto roll = vt::dy_spillover(fixture_data(), 150, 10, "c", 2);
    auto d = fixture_data();
    vt::MultivariateTSData sub;
    sub.columns.resize(3);
    for (int j = 0; j < 3; ++j) {
        sub.columns[j].assign(d.columns[j].begin() + 100, d.columns[j].end());
    }
    auto stat = vt::dy_spillover_static(sub, 10, "c", 2);
    EXPECT_NEAR(roll.tci, stat.tci, 1e-10);
    for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(roll.net_spillover[j], stat.net_spillover[j], 1e-10);
}

// 11. lag=0 每窗口 IC 自动选阶 (滚动路径仍产出)
TEST(DySpillover, RollingAutoLag) {
    auto res = vt::dy_spillover(fixture_data(), 150, 10, "c", 0);
    EXPECT_EQ(res.tci_path.size(), 101u);
    EXPECT_GT(res.lag, 0u);
    // 与固定 lag=2 路径相近但可能不同 (IC 每窗口可能选不同阶)
    auto fixed = vt::dy_spillover(fixture_data(), 150, 10, "c", 2);
    EXPECT_NEAR(res.tci, fixed.tci, 2.0);  // 方向性 (同数据同 H)
}

// 12. 频率默认表知识锚 (§13-a): {日 200/周 200/月 60} × H {10/10/12}
//     实现为 API 层无默认 (window 必填), 表值记录于注释 — 断言文档性行为:
//     月频 T=60 月窗 + H=12 可配执行
TEST(DySpillover, MonthlyConventionUsable) {
    auto d = fixture_data();
    vt::MultivariateTSData sub;
    sub.columns.resize(3);
    for (int j = 0; j < 3; ++j) {
        sub.columns[j].assign(d.columns[j].begin(), d.columns[j].begin() + 60);
    }
    // 月频惯例: window=60 → sub.T()=60 全样本; H=12
    auto res = vt::dy_spillover_static(sub, 12, "c", 1);
    EXPECT_TRUE(std::isfinite(res.tci));
    EXPECT_EQ(res.horizon, 12u);
}
