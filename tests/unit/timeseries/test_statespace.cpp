// =============================================================================
// test_statespace.cpp - 状态空间模型 + 常规 Kalman 滤波 (Phase 8 M2 Step 2,
// 16 用例)
//
// 基准: 全部手推黄金锚 (spec §1.3: AR(1) 前两步 v/F/loglik 1e-12 /
//   AR(1) 平稳 P₀ = σ²/(1−φ²) 1e-12 / Lyapunov 守卫行为断言);
//   statsmodels SARIMAX llf 1e-8 对照属 §3.2 (arima_statespace + verify_
//   statespace.py), 本套件不依赖外部基准
//
// 手推锚推导 (全部闭式, 非循环):
//   AR(1) φ=0.5 σ²=0.9 (stationary): P₁ = 1.2, a₁ = 0
//     t=0: v=1.0, F=1.2, K=0.5, a₂=0.5, P₂=φ²(P₁−P₁²/F₁)+σ²=0.9 (H=0 稳态)
//     contrib₁ = −½(log 1.2 + 1/1.2 + log 2π) = −1.4267659782683167
//     contrib₂ = −½(log 0.9 + 1/0.9 + log 2π) = −1.4218138309313151
//   ARMA(1,1) φ=0.5 θ=0.3 σ²=1 (stationary, T=[[φ,1],[0,0]], R=[1,θ]'):
//     Lyapunov 闭式 P₁ = [[σ²(1+2φθ+θ²)/(1−φ²), σ²θ],[σ²θ, σ²θ²]]
//                     = [[139/75, 3/10],[3/10, 9/100]]
//     t=0: K₁=[92/139,0], a₂=[92/139,0], P₂=[[3619/3475, 3/10],[3/10, 9/100]]
//     t=1: v₂ = −323/278, F₂ = 3619/3475
//   Known+截距 m=1 (T=0.5, Q=0.9, H=0.2, c=0.1, d=−0.05, a₀=0.3, P₀=0.7):
//     t=0: v=0.75, F=0.9, K=7/18, a'=13/24, P'=169/180
//     t=1: v=−11/120, F=41/36
//
// 幻觉点: SS2 (DK 更新形增益含 T 前乘) / SS3 (R m×r 维度校验) / SS4
//   (Lyapunov 守卫: φ=1.01 非正定 throw + φ=1.0 奇异 throw) / SS6 (F≤0
//   throw) / SS8 (单变量常数 log 2π) / SS1 (κ=1e6 首期 F=κ) /
//   burn = 贡献求和起点 (滤波全程运行)
// =============================================================================

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/statespace/ss_model.hpp"

namespace ss = cpphub::v1::timeseries::statespace;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

const Real kLog2Pi = std::log(2.0 * 3.14159265358979323846);

/// AR(1) 状态空间: y_t = φ y_{t−1} + ε_t (m=1, r=1)
ss::StateSpaceModel ar1_model(Real phi, Real sigma2) {
    ss::StateSpaceModel m;
    m.Z = Eigen::MatrixXd(1, 1);
    m.Z << 1.0;
    m.T = Eigen::MatrixXd(1, 1);
    m.T << phi;
    m.R = Eigen::MatrixXd(1, 1);
    m.R << 1.0;
    m.Q = Eigen::MatrixXd(1, 1);
    m.Q << sigma2;
    m.H = Eigen::MatrixXd(1, 1);
    m.H << 0.0;
    return m;
}

/// ARMA(1,1) Harvey 形: T=[[φ,1],[0,0]], R=[1,θ]', Z=[1,0], H=0 (SS2)
ss::StateSpaceModel arma11_model(Real phi, Real theta, Real sigma2) {
    ss::StateSpaceModel m;
    m.Z = Eigen::MatrixXd(1, 2);
    m.Z << 1.0, 0.0;
    m.T = Eigen::MatrixXd(2, 2);
    m.T << phi, 1.0, 0.0, 0.0;
    m.R = Eigen::MatrixXd(2, 1);
    m.R << 1.0, theta;
    m.Q = Eigen::MatrixXd(1, 1);
    m.Q << sigma2;
    m.H = Eigen::MatrixXd(1, 1);
    m.H << 0.0;
    return m;
}

/// Known 初始化 + 双截距模型 (m=1): T=0.5, Q=0.9, H=0.2, c=0.1, d=−0.05,
/// a₀=0.3, P₀=0.7
ss::StateSpaceModel known_intercept_model() {
    ss::StateSpaceModel m = ar1_model(0.5, 0.9);
    m.H << 0.2;
    m.d = Eigen::VectorXd(1);
    m.d << -0.05;
    m.c = Eigen::VectorXd(1);
    m.c << 0.1;
    m.init.method = ss::SSInitMethod::Known;
    m.init.a0 = Eigen::VectorXd(1);
    m.init.a0 << 0.3;
    m.init.P0 = Eigen::MatrixXd(1, 1);
    m.init.P0 << 0.7;
    return m;
}

}  // namespace

// --- AR(1) 黄金锚 (spec §3.1 测试要点, 1e-12) ---

TEST(Statespace, Ar1GoldenAnchorTwoSteps) {
    ss::StateSpaceModel m = ar1_model(0.5, 0.9);
    m.init.method = ss::SSInitMethod::Stationary;
    const std::vector<Real> y = {1.0, -0.5};

    const ss::KalmanResult r = ss::kalman_filter(y, m);
    ASSERT_TRUE(r.ok);
    ASSERT_EQ(r.innovations.size(), 2u);
    ASSERT_EQ(r.pred_var.size(), 2u);
    ASSERT_EQ(r.state_pred.size(), 2u);
    ASSERT_EQ(r.state_cov_pred.size(), 2u);

    // t=0: a₁=0 (平稳均值, c 空), P₁=σ²/(1−φ²)=1.2 → v₁=y₁, F₁=P₁
    EXPECT_NEAR(r.state_pred[0](0), 0.0, 1e-12);
    EXPECT_NEAR(r.state_cov_pred[0](0, 0), 1.2, 1e-12);
    EXPECT_NEAR(r.innovations[0], 1.0, 1e-12);
    EXPECT_NEAR(r.pred_var[0], 1.2, 1e-12);

    // t=1: a₂ = K₁v₁ = 0.5·1 = 0.5; P₂ = 0.9 (H=0 时 P 即条件预测方差)
    EXPECT_NEAR(r.state_pred[1](0), 0.5, 1e-12);
    EXPECT_NEAR(r.state_cov_pred[1](0, 0), 0.9, 1e-12);
    EXPECT_NEAR(r.innovations[1], -1.0, 1e-12);
    EXPECT_NEAR(r.pred_var[1], 0.9, 1e-12);

    // 全程 loglik (burn=0) 手推字面量 (非公式重建)
    EXPECT_NEAR(r.loglik, -2.8485798091996318, 1e-12);
}

TEST(Statespace, Ar1ContribsExact) {
    ss::StateSpaceModel m = ar1_model(0.5, 0.9);
    m.init.method = ss::SSInitMethod::Stationary;
    const ss::KalmanResult r = ss::kalman_filter({1.0, -0.5}, m);

    ASSERT_EQ(r.loglik_contribs.size(), 2u);
    EXPECT_NEAR(r.loglik_contribs[0], -1.4267659782683167, 1e-12);
    EXPECT_NEAR(r.loglik_contribs[1], -1.4218138309313151, 1e-12);
    // 公式形交叉: −½(log F_t + v_t²/F_t + log 2π)
    EXPECT_NEAR(r.loglik_contribs[0],
                -0.5 * (std::log(1.2) + 1.0 / 1.2 + kLog2Pi), 1e-15);
    EXPECT_NEAR(r.loglik_contribs[1],
                -0.5 * (std::log(0.9) + 1.0 / 0.9 + kLog2Pi), 1e-15);
}

TEST(Statespace, Ar1BurnSemantics) {
    ss::StateSpaceModel m = ar1_model(0.5, 0.9);
    m.init.method = ss::SSInitMethod::Stationary;
    const ss::KalmanResult r = ss::kalman_filter({1.0, -0.5}, m, 1);

    EXPECT_EQ(r.n_burn, 1u);
    EXPECT_EQ(r.n_obs, 2u);
    EXPECT_EQ(r.loglik_contribs.size(), 1u);
    // burn 前后均记录 (滤波全程运行, 非丢观测 — probe [B 补充])
    EXPECT_EQ(r.innovations.size(), 2u);
    EXPECT_EQ(r.pred_var.size(), 2u);
    // loglik = burn 后贡献和 = 仅第二期
    EXPECT_NEAR(r.loglik, -1.4218138309313151, 1e-12);
    EXPECT_NEAR(r.loglik, r.loglik_contribs[0], 1e-15);
}

TEST(Statespace, Ar1StationaryCovExact) {
    const ss::StateSpaceModel m = ar1_model(0.5, 0.9);
    const Eigen::MatrixXd P = ss::stationary_cov(m);
    ASSERT_EQ(P.rows(), 1);
    ASSERT_EQ(P.cols(), 1);
    EXPECT_NEAR(P(0, 0), 0.9 / (1.0 - 0.5 * 0.5), 1e-12);
}

TEST(Statespace, Arma11StationaryCovExact) {
    const ss::StateSpaceModel m = arma11_model(0.5, 0.3, 1.0);
    const Eigen::MatrixXd P = ss::stationary_cov(m);
    ASSERT_EQ(P.rows(), 2);
    ASSERT_EQ(P.cols(), 2);
    // 手推闭式: [[σ²(1+2φθ+θ²)/(1−φ²), σ²θ],[σ²θ, σ²θ²]]
    EXPECT_NEAR(P(0, 0), (1.0 + 2.0 * 0.5 * 0.3 + 0.09) / (1.0 - 0.25),
                1e-12);
    EXPECT_NEAR(P(0, 1), 0.3, 1e-12);
    EXPECT_NEAR(P(1, 0), 0.3, 1e-12);
    EXPECT_NEAR(P(1, 1), 0.09, 1e-12);
}

// --- Lyapunov 守卫 (SS4) ---

TEST(Statespace, LyapunovGuardExplosiveThrows) {
    // φ=1.01: 方程有解但非正定 (静默负方差风险 — 守卫必须 throw)
    const ss::StateSpaceModel m = ar1_model(1.01, 0.9);
    EXPECT_THROW(ss::stationary_cov(m), std::runtime_error);
}

TEST(Statespace, LyapunovGuardUnitRootThrows) {
    // φ=1.0: I − T⊗T 奇异 (单位根)
    const ss::StateSpaceModel m = ar1_model(1.0, 0.9);
    EXPECT_THROW(ss::stationary_cov(m), std::runtime_error);
}

// --- F 正定性守卫 (SS6) 与 diffuse (SS1) ---

TEST(Statespace, NonPositiveFThrows) {
    ss::StateSpaceModel m = ar1_model(0.5, 0.9);
    m.init.method = ss::SSInitMethod::Known;
    m.init.P0 = Eigen::MatrixXd::Zero(1, 1);  // P₁=0 且 H=0 → F₁=0
    m.init.a0 = Eigen::VectorXd::Zero(1);
    EXPECT_THROW(ss::kalman_filter({1.0, 0.5}, m), std::runtime_error);
}

TEST(Statespace, DiffuseFirstFIsKappa) {
    ss::StateSpaceModel m = arma11_model(0.5, 0.3, 1.0);
    // 默认 ApproximateDiffuse, κ=1e6 (statsmodels 口径, K1)
    const ss::KalmanResult r = ss::kalman_filter({1.0, -0.5}, m);
    ASSERT_TRUE(r.ok);
    EXPECT_TRUE(r.message.empty());
    // a₁=0, P₁=κI, Z=[1,0] → F₁=κ (精确); v₁=y₁
    EXPECT_DOUBLE_EQ(r.pred_var[0], 1e6);
    EXPECT_DOUBLE_EQ(r.innovations[0], 1.0);
    EXPECT_DOUBLE_EQ(r.state_cov_pred[0](0, 0), 1e6);
    EXPECT_DOUBLE_EQ(r.state_cov_pred[0](1, 1), 1e6);
    EXPECT_DOUBLE_EQ(r.state_cov_pred[0](0, 1), 0.0);
}

TEST(Statespace, StationaryFallsBackToDiffuse) {
    ss::StateSpaceModel m = ar1_model(1.01, 0.9);
    m.init.method = ss::SSInitMethod::Stationary;
    const ss::KalmanResult r = ss::kalman_filter({1.0, 0.5}, m);
    // 非静默回退 (spec §3.1: message 记录)
    EXPECT_TRUE(r.ok);
    EXPECT_FALSE(r.message.empty());
    EXPECT_NE(r.message.find("ApproximateDiffuse"), std::string::npos);
    EXPECT_DOUBLE_EQ(r.pred_var[0], 1e6);  // 回退后 F₁ = κ
    EXPECT_DOUBLE_EQ(r.innovations[0], 1.0);
}

// --- Known 初始化 + 截距 (d/c) 递推 ---

TEST(Statespace, KnownInitAndIntercepts) {
    const ss::StateSpaceModel m = known_intercept_model();
    const std::vector<Real> y = {1.0, 0.4};
    const ss::KalmanResult r = ss::kalman_filter(y, m);
    ASSERT_TRUE(r.ok);

    // t=0: v = y − Za₀ − d = 1.0 − 0.3 + 0.05 = 0.75; F = P₀ + H = 0.9
    EXPECT_NEAR(r.innovations[0], 0.75, 1e-12);
    EXPECT_NEAR(r.pred_var[0], 0.9, 1e-12);
    // 更新 (K = TP₀Z/F = 7/18): a' = Ta₀ + Kv + c = 13/24;
    // P' = TP₀(T−KZ)' + RQR' = 169/180
    EXPECT_NEAR(r.state_pred[1](0), 13.0 / 24.0, 1e-12);
    EXPECT_NEAR(r.state_cov_pred[1](0, 0), 169.0 / 180.0, 1e-12);
    // t=1: v = 0.4 − 13/24 + 0.05 = −11/120; F = 169/180 + 1/5 = 41/36
    EXPECT_NEAR(r.innovations[1], -11.0 / 120.0, 1e-12);
    EXPECT_NEAR(r.pred_var[1], 41.0 / 36.0, 1e-12);

    // loglik 闭式 (burn=0, SS8 常数含 log 2π)
    const Real c0 =
        -0.5 * (std::log(0.9) + 0.75 * 0.75 / 0.9 + kLog2Pi);
    const Real c1 = -0.5 * (std::log(41.0 / 36.0) +
                            (11.0 / 120.0) * (11.0 / 120.0) /
                                (41.0 / 36.0) +
                            kLog2Pi);
    ASSERT_EQ(r.loglik_contribs.size(), 2u);
    EXPECT_NEAR(r.loglik_contribs[0], c0, 1e-14);
    EXPECT_NEAR(r.loglik_contribs[1], c1, 1e-14);
    EXPECT_NEAR(r.loglik, c0 + c1, 1e-12);
}

TEST(Statespace, ContribContainsLog2Pi) {
    // SS8: 单变量常数 = log 2π — contrib 与去常数部分之差精确等于 −½log 2π
    ss::StateSpaceModel m = ar1_model(0.5, 0.9);
    m.init.method = ss::SSInitMethod::Stationary;
    const ss::KalmanResult r = ss::kalman_filter({1.0, -0.5}, m);
    ASSERT_GE(r.loglik_contribs.size(), 1u);
    const Real no_const = -0.5 * (std::log(1.2) + 1.0 / 1.2);
    EXPECT_NEAR(r.loglik_contribs[0] - no_const, -0.5 * kLog2Pi, 1e-15);
}

TEST(Statespace, LoglikEqualsSumOfContribs) {
    // burn 语义: 贡献求和起点; 同一滤波轨道 (innovations 逐值一致)
    const ss::StateSpaceModel m = arma11_model(0.5, 0.3, 1.0);
    const std::vector<Real> y = {1.0, -0.5, 0.8, 0.2};

    const ss::KalmanResult r0 = ss::kalman_filter(y, m, 0);
    Real s = 0.0;
    for (Real c : r0.loglik_contribs) s += c;
    EXPECT_NEAR(r0.loglik, s, 1e-12);
    EXPECT_EQ(r0.loglik_contribs.size(), 4u);

    const ss::KalmanResult r2 = ss::kalman_filter(y, m, 2);
    EXPECT_EQ(r2.loglik_contribs.size(), 2u);
    EXPECT_EQ(r2.innovations.size(), 4u);
    EXPECT_EQ(r2.pred_var.size(), 4u);
    EXPECT_EQ(r2.n_burn, 2u);
    EXPECT_EQ(r2.n_obs, 4u);
    for (Size t = 0; t < 4u; ++t) {
        EXPECT_DOUBLE_EQ(r0.innovations[t], r2.innovations[t]);
        EXPECT_DOUBLE_EQ(r0.pred_var[t], r2.pred_var[t]);
    }
}

TEST(Statespace, StatePredAlignment) {
    // state_pred[0] = 初始先验; innovations[t] = y[t] − Z·state_pred[t] − d
    const ss::StateSpaceModel m = known_intercept_model();
    const std::vector<Real> y = {1.0, 0.4};
    const ss::KalmanResult r = ss::kalman_filter(y, m);
    ASSERT_EQ(r.state_pred.size(), 2u);
    EXPECT_NEAR(r.state_pred[0](0), 0.3, 1e-12);  // = a0
    for (Size t = 0; t < 2u; ++t) {
        const Real v_rec =
            y[t] - (m.Z * r.state_pred[t])(0, 0) - m.d(0);
        EXPECT_NEAR(r.innovations[t], v_rec, 1e-12) << "t=" << t;
    }
}

TEST(Statespace, Arma11TwoStepTrace) {
    // m=2 黄金锚: Harvey 形 ARMA(1,1) 平稳初始化两步递推 (全部闭式分数)
    ss::StateSpaceModel m = arma11_model(0.5, 0.3, 1.0);
    m.init.method = ss::SSInitMethod::Stationary;
    const std::vector<Real> y = {1.0, -0.5};

    const ss::KalmanResult r = ss::kalman_filter(y, m);
    ASSERT_TRUE(r.ok);
    EXPECT_TRUE(r.message.empty());  // ρ(T)=0.5 < 1, 无回退

    // t=0: v₁ = y₁ (a₁=0); F₁ = P₁[0,0] = 139/75
    EXPECT_NEAR(r.innovations[0], 1.0, 1e-12);
    EXPECT_NEAR(r.pred_var[0], 139.0 / 75.0, 1e-12);
    // t=1: a₂ = K₁v₁ = [92/139, 0];
    //      P₂ = [[3619/3475, 3/10],[3/10, 9/100]]
    ASSERT_EQ(r.state_pred[1].size(), 2);
    EXPECT_NEAR(r.state_pred[1](0), 92.0 / 139.0, 1e-12);
    EXPECT_NEAR(r.state_pred[1](1), 0.0, 1e-12);
    EXPECT_NEAR(r.state_cov_pred[1](0, 0), 3619.0 / 3475.0, 1e-12);
    EXPECT_NEAR(r.state_cov_pred[1](0, 1), 0.3, 1e-12);
    EXPECT_NEAR(r.state_cov_pred[1](1, 0), 0.3, 1e-12);
    EXPECT_NEAR(r.state_cov_pred[1](1, 1), 0.09, 1e-12);
    // v₂ = −0.5 − 92/139 = −323/278; F₂ = 3619/3475
    EXPECT_NEAR(r.innovations[1], -323.0 / 278.0, 1e-12);
    EXPECT_NEAR(r.pred_var[1], 3619.0 / 3475.0, 1e-12);

    // loglik 闭式
    const Real c0 =
        -0.5 * (std::log(139.0 / 75.0) + 1.0 / (139.0 / 75.0) + kLog2Pi);
    const Real c1 = -0.5 * (std::log(3619.0 / 3475.0) +
                            (323.0 / 278.0) * (323.0 / 278.0) /
                                (3619.0 / 3475.0) +
                            kLog2Pi);
    EXPECT_NEAR(r.loglik, c0 + c1, 1e-12);
}

// --- 输入校验 (SS3 维度 + Known init 契约) ---

TEST(Statespace, ValidationThrows) {
    const ss::StateSpaceModel base = arma11_model(0.5, 0.3, 1.0);
    const std::vector<Real> y = {1.0, -0.5};

    ss::StateSpaceModel bad_z = base;
    bad_z.Z = Eigen::MatrixXd::Zero(2, 2);  // Z 非 (1×m)
    EXPECT_THROW(ss::kalman_filter(y, bad_z), std::invalid_argument);

    ss::StateSpaceModel bad_r = base;
    bad_r.R = Eigen::MatrixXd::Zero(2, 3);  // r=3 > m=2 (SS3)
    EXPECT_THROW(ss::kalman_filter(y, bad_r), std::invalid_argument);

    ss::StateSpaceModel bad_r2 = base;
    bad_r2.R = Eigen::MatrixXd::Zero(1, 1);  // R.rows ≠ m
    EXPECT_THROW(ss::kalman_filter(y, bad_r2), std::invalid_argument);

    ss::StateSpaceModel bad_q = base;
    bad_q.Q = Eigen::MatrixXd::Zero(2, 2);  // Q 非 (r×r)
    EXPECT_THROW(ss::kalman_filter(y, bad_q), std::invalid_argument);

    ss::StateSpaceModel bad_h = base;
    bad_h.H = Eigen::MatrixXd::Zero(2, 2);  // H 非 (1×1)
    EXPECT_THROW(ss::kalman_filter(y, bad_h), std::invalid_argument);

    ss::StateSpaceModel bad_t = base;
    bad_t.T = Eigen::MatrixXd::Zero(2, 3);  // T 非方阵
    EXPECT_THROW(ss::kalman_filter(y, bad_t), std::invalid_argument);

    ss::StateSpaceModel bad_d = base;
    bad_d.d = Eigen::VectorXd::Zero(2);  // d 非 size 1
    EXPECT_THROW(ss::kalman_filter(y, bad_d), std::invalid_argument);

    ss::StateSpaceModel bad_c = base;
    bad_c.c = Eigen::VectorXd::Zero(3);  // c 非 size m
    EXPECT_THROW(ss::kalman_filter(y, bad_c), std::invalid_argument);

    EXPECT_THROW(ss::kalman_filter(std::vector<Real>{}, base),
                 std::invalid_argument);  // 空 y

    const std::vector<Real> y_nan = {1.0, std::nan("")};
    EXPECT_THROW(ss::kalman_filter(y_nan, base),
                 std::invalid_argument);  // y 含 NaN

    // Known init 契约: P0 必填 m×m; a0 空或 m
    ss::StateSpaceModel k0 = base;
    k0.init.method = ss::SSInitMethod::Known;
    EXPECT_THROW(ss::kalman_filter(y, k0), std::invalid_argument);  // P0 空
    k0.init.P0 = Eigen::MatrixXd::Zero(1, 1);  // P0 尺寸错
    EXPECT_THROW(ss::kalman_filter(y, k0), std::invalid_argument);
    k0.init.P0 = Eigen::MatrixXd::Identity(2, 2);
    k0.init.a0 = Eigen::VectorXd::Zero(3);  // a0 尺寸错
    EXPECT_THROW(ss::kalman_filter(y, k0), std::invalid_argument);

    // 系统矩阵含 NaN
    ss::StateSpaceModel bad_nan = base;
    bad_nan.T(0, 0) = std::nan("");
    EXPECT_THROW(ss::kalman_filter(y, bad_nan), std::invalid_argument);
}
