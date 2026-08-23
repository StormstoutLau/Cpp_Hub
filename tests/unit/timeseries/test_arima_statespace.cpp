// =============================================================================
// test_arima_statespace.cpp - ARIMA Harvey 形 + 嵌入差分 + 精确似然 MLE
// (Phase 8 M2 Step 3, 14 用例)
//
// 基准: ss_baseline.inc (gen_phase8_ss_inc.py 自动生成, statsmodels 0.14.6
//   对照 2026-08-23): 黄金锚 llf 五组 (固定参数, concentrate_scale=False,
//   burn=d 口径) + MLE 四组 (参数/llf/aic) + Monahan 变换锚 + 互检 gap 锚
//
// 语义冻结来源 (tests/fixtures/timeseries/):
//   probe_sarimax.py — Harvey 形布局/burn 公式/初始化/参数变换
//   probe_d2_const.py — d=0+trend='c' 常数直入 state_intercept; d=2 布局
//   verify_drift_init.py — drift 混合初始化 (差分块 diffuse ⊕ ARMA 块
//     stationary + a₀ = (I−T_sub)⁻¹c_sub)
//
// 幻觉点: SS2 (Harvey 形唯一 — AR 入 T 第一列 + Z 前置) / SS5 (双差分口径
//   似然不同, 禁混比; burn=0 含 diffuse 伪贡献) / SS9 (κ=1e6 仅初始化) /
//   K1 (innovations 常驻互检 ≤1e-5) / K4 (σ² 非集中化) / Monahan MA 负号
// =============================================================================

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/statespace/arima_statespace.hpp"
#include "ss_baseline.inc"

namespace ss = cpphub::v1::timeseries::statespace;
namespace arima_ns = cpphub::v1::timeseries::arima;
namespace ssb = cpphub::v1::timeseries::ss_baseline;
using cpphub::v1::Real;
using cpphub::v1::Size;

namespace {

template <Size N>
std::vector<Real> vec(const double (&a)[N]) {
    return std::vector<Real>(std::begin(a), std::end(a));
}

std::vector<Real> data120() { return vec(ssb::DATA); }
std::vector<Real> data300() { return vec(ssb::DATA300); }

std::vector<Real> diff1(const std::vector<Real>& y) {
    std::vector<Real> dy(y.size() - 1);
    for (Size t = 1; t < y.size(); ++t) dy[t - 1] = y[t] - y[t - 1];
    return dy;
}

}  // namespace

// ---- 1. ARMA(1,1) Harvey 形三矩阵手推黄金锚 (SS2) + p=q=0 边界 ----
TEST(ArimaStatespace, Arma11HarveyMatrices) {
    const ss::StateSpaceModel m = ss::make_arma_statespace({0.55}, {0.35}, 0.9);
    ASSERT_EQ(m.T.rows(), 2);
    ASSERT_EQ(m.T.cols(), 2);
    // T = [[φ, 1], [0, 0]] (AR 入第一列 + 次对角单位, 上 Hessenberg)
    EXPECT_NEAR(m.T(0, 0), 0.55, 1e-14);
    EXPECT_NEAR(m.T(0, 1), 1.0, 1e-14);
    EXPECT_NEAR(m.T(1, 0), 0.0, 1e-14);
    EXPECT_NEAR(m.T(1, 1), 0.0, 1e-14);
    // R = (1, θ)', Z = (1, 0), Q = [σ²], H = 0 (SS6)
    EXPECT_NEAR(m.R(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(m.R(1, 0), 0.35, 1e-14);
    EXPECT_NEAR(m.Z(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(m.Z(0, 1), 0.0, 1e-14);
    EXPECT_NEAR(m.Q(0, 0), 0.9, 1e-14);
    EXPECT_NEAR(m.H(0, 0), 0.0, 1e-14);
    EXPECT_EQ(m.init.method, ss::SSInitMethod::Stationary);

    // v1.8 边界: 纯白噪声/纯随机游走不支持 (throw 非静默)
    EXPECT_THROW(ss::make_arma_statespace({}, {}, 0.9),
                 std::invalid_argument);
    EXPECT_THROW(ss::make_arma_statespace({}, {}, -0.1),
                 std::invalid_argument);  // σ² ≤ 0
}

// ---- 2. ARIMA(1,1,1) 嵌入差分黄金锚 (probe [A2] 逐元素冻结) ----
TEST(ArimaStatespace, Arima111EmbeddedMatrices) {
    const ss::StateSpaceModel m =
        ss::make_arima_statespace({0.55}, {0.35}, 1, 0.9, false, false);
    ASSERT_EQ(m.T.rows(), 3);
    // T = [[1, 1, 0], [0, φ, 1], [0, 0, 0]]
    //   (差分块 triu + 链接列 T[:d,d]=1 + ARMA companion)
    EXPECT_NEAR(m.T(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(m.T(0, 1), 1.0, 1e-14);
    EXPECT_NEAR(m.T(0, 2), 0.0, 1e-14);
    EXPECT_NEAR(m.T(1, 0), 0.0, 1e-14);
    EXPECT_NEAR(m.T(1, 1), 0.55, 1e-14);
    EXPECT_NEAR(m.T(1, 2), 1.0, 1e-14);
    EXPECT_NEAR(m.T(2, 0), 0.0, 1e-14);
    EXPECT_NEAR(m.T(2, 1), 0.0, 1e-14);
    EXPECT_NEAR(m.T(2, 2), 0.0, 1e-14);
    // R = (0, 1, θ)', Z = (1, 1, 0)
    EXPECT_NEAR(m.R(0, 0), 0.0, 1e-14);
    EXPECT_NEAR(m.R(1, 0), 1.0, 1e-14);
    EXPECT_NEAR(m.R(2, 0), 0.35, 1e-14);
    EXPECT_NEAR(m.Z(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(m.Z(0, 1), 1.0, 1e-14);
    EXPECT_NEAR(m.Z(0, 2), 0.0, 1e-14);
    // 平稳 AR → 混合初始化 (Known 注入: 差分块 diffuse κ ⊕ ARMA 块平稳)
    EXPECT_EQ(m.init.method, ss::SSInitMethod::Known);
    ASSERT_EQ(m.init.P0.rows(), 3);
    EXPECT_NEAR(m.init.P0(0, 0), 1e6, 1.0);  // κ (SS1)
}

// ---- 3. d=2 嵌入差分布局 (probe_d2_const [P2] 三合同款冻结) ----
TEST(ArimaStatespace, D2EmbeddedMatrices) {
    // ARIMA(0,2,1): m1 = max(0, 2) = 2, m = 4
    // T = triu(2×2 全 1) ⊕ 链接列 T[:2,2]=1 ⊕ ARMA companion(2)
    //   [1, 1, 1, 0]
    //   [0, 1, 1, 0]
    //   [0, 0, 0, 1]  ← companion 超对角 (无 AR → T[2,2]=0)
    //   [0, 0, 0, 0]
    const ss::StateSpaceModel m01 =
        ss::make_arima_statespace({}, {0.4}, 2, 1.0, false, false);
    ASSERT_EQ(m01.T.rows(), 4);
    ASSERT_EQ(m01.T.cols(), 4);
    EXPECT_NEAR(m01.T(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(m01.T(0, 1), 1.0, 1e-14);
    EXPECT_NEAR(m01.T(1, 0), 0.0, 1e-14);
    EXPECT_NEAR(m01.T(1, 1), 1.0, 1e-14);
    EXPECT_NEAR(m01.T(0, 2), 1.0, 1e-14);  // 链接列
    EXPECT_NEAR(m01.T(1, 2), 1.0, 1e-14);
    EXPECT_NEAR(m01.T(2, 2), 0.0, 1e-14);  // 无 AR
    EXPECT_NEAR(m01.T(2, 3), 1.0, 1e-14);  // companion 超对角
    EXPECT_NEAR(m01.T(3, 3), 0.0, 1e-14);
    // R = (0, 0, 1, θ)', Z = (1, 1, 1, 0)
    EXPECT_NEAR(m01.R(0, 0), 0.0, 1e-14);
    EXPECT_NEAR(m01.R(1, 0), 0.0, 1e-14);
    EXPECT_NEAR(m01.R(2, 0), 1.0, 1e-14);
    EXPECT_NEAR(m01.R(3, 0), 0.4, 1e-14);
    EXPECT_NEAR(m01.Z(0, 2), 1.0, 1e-14);
    EXPECT_NEAR(m01.Z(0, 3), 0.0, 1e-14);
    // p=q=0 边界 (纯随机游走 throw)
    EXPECT_THROW(ss::make_arima_statespace({}, {}, 2, 1.0, false, false),
                 std::invalid_argument);

    // ARIMA(1,2,1): m1 = 2, m = 4; AR 首列 + companion 超对角
    const ss::StateSpaceModel m11 =
        ss::make_arima_statespace({0.5}, {0.3}, 2, 1.0, false, false);
    ASSERT_EQ(m11.T.rows(), 4);
    EXPECT_NEAR(m11.T(2, 2), 0.5, 1e-14);  // AR 首列 (d=2 偏移)
    EXPECT_NEAR(m11.T(2, 3), 1.0, 1e-14);  // companion 超对角
    EXPECT_NEAR(m11.T(3, 2), 0.0, 1e-14);
    EXPECT_NEAR(m11.T(3, 3), 0.0, 1e-14);
    EXPECT_NEAR(m11.R(2, 0), 1.0, 1e-14);
    EXPECT_NEAR(m11.R(3, 0), 0.3, 1e-14);
    EXPECT_NEAR(m11.Z(0, 3), 0.0, 1e-14);
}

// ---- 4. KF llf vs statsmodels 黄金锚: ARMA(1,1) burn=0 (1e-8) ----
TEST(ArimaStatespace, KfLlfArma11Golden) {
    const ss::StateSpaceModel m = ss::make_arma_statespace({0.55}, {0.35}, 0.9);
    const ss::KalmanResult r = ss::kalman_filter(data120(), m, 0);
    ASSERT_TRUE(r.ok);
    // statsmodels SARIMAX(y, order=(1,0,1)).loglike([.55,.35,.9]), burn=0
    EXPECT_NEAR(r.loglik, ssb::GOLD_LL_ARMA11, 1e-8);
}

// ---- 5. KF llf vs statsmodels 黄金锚: ARIMA(1,1,1) embedded burn=d ----
TEST(ArimaStatespace, KfLlfArima111Golden) {
    const ss::StateSpaceModel m =
        ss::make_arima_statespace({0.55}, {0.35}, 1, 0.9, false, false);
    const ss::KalmanResult r = ss::kalman_filter(data120(), m, 1);
    ASSERT_TRUE(r.ok);
    // statsmodels SARIMAX(y, order=(1,1,1)).loglike([...]), burn=1
    // (开放问题 #3 裁决: 默认 burn = d = statsmodels None 公式)
    //
    // 阈值 1e-7 (非 1e-8) 仲裁依据 (tmp/arbiter_llf.py 2026-08-23):
    //   本库 KF 与独立 manual KF (probe [D] 同构) 逐位一致 (diff=0),
    //   2.81e-8 偏差归属 statsmodels — 其 Cython KF 稳态捷径在 P 收敛后
    //   冻结增益 (simple 口径 t≥60 后 |ΔF_t|<1e-14 全 58 期), 截断误差
    //   累积 ~2.8e-8; κ=1e6 diffuse 放大 (probe [D] 同结论, 阈值同 1e-7)
    EXPECT_NEAR(r.loglik, ssb::GOLD_LL_ARIMA111, 1e-7);
}

// ---- 6. KF llf vs statsmodels 黄金锚: ARIMA(0,1,1)+drift 混合初始化 ----
TEST(ArimaStatespace, KfLlfDriftGolden) {
    // statsmodels params = [drift=0.2, θ=0.5, σ²=0.9] (trend='c' 在前)
    const ss::StateSpaceModel m =
        ss::make_arima_statespace({}, {0.5}, 1, 0.9, true, false, 0.2);
    const ss::KalmanResult r = ss::kalman_filter(data120(), m, 1);
    ASSERT_TRUE(r.ok);
    EXPECT_NEAR(r.loglik, ssb::GOLD_LL_ARIMA011_DRIFT, 1e-8);

    // 常数注入位置守卫: c[d]=drift (ARMA 块首行), 非观测截距 d_t
    ASSERT_EQ(m.c.size(), 3);
    EXPECT_NEAR(m.c(0), 0.0, 1e-14);
    EXPECT_NEAR(m.c(1), 0.2, 1e-14);
    EXPECT_NEAR(m.c(2), 0.0, 1e-14);
    ASSERT_TRUE(m.d.size() == 0 || std::abs(m.d(0)) < 1e-14);
    // 混合初始化 a₀^ARMA = (I−T_sub)⁻¹c_sub (verify_drift_init 冻结)
    // ARMA(0,1) 子块: T_sub=[0], c_sub=[0.2] → a₀[1] = 0.2
    ASSERT_EQ(m.init.a0.size(), 3);
    EXPECT_NEAR(m.init.a0(1), 0.2, 1e-12);
}

// ---- 7. SS5 双差分口径: embedded vs simple 时间平移等价 + burn=0 伪贡献 ----
TEST(ArimaStatespace, SimpleDifferencingDualCaliber) {
    const std::vector<Real> y = data120();
    const std::vector<Real> dy = diff1(y);

    // embedded (levels, burn=1) vs simple (dy 样本, burn=0)
    const ss::StateSpaceModel me =
        ss::make_arima_statespace({0.55}, {0.35}, 1, 0.9, false, false);
    const ss::KalmanResult re = ss::kalman_filter(y, me, 1);
    const ss::StateSpaceModel ms =
        ss::make_arima_statespace({0.55}, {0.35}, 0, 0.9, false, true);
    const ss::KalmanResult rs = ss::kalman_filter(dy, ms, 0);
    ASSERT_TRUE(re.ok);
    ASSERT_TRUE(rs.ok);

    // 双口径锚: gap = O(1/κ) ≈ 5.3e-7 (κ=1e6 有限 diffuse 尾差) —
    // SS5: 两口径似然不同 (非零 gap), 数值差异显式记录, 禁混比断言
    // 阈值 1e-7 同 KfLlfArima111Golden 仲裁 (statsmodels 稳态捷径
    // 2.8e-8; 本库 = 精确递推, 与 manual KF 逐位一致)
    EXPECT_NEAR(re.loglik, ssb::GOLD_LL_ARIMA111, 1e-7);
    EXPECT_NEAR(rs.loglik, ssb::GOLD_LL_ARIMA111_SIMPLE, 1e-7);
    const Real gap_dual = std::abs(re.loglik - rs.loglik);
    EXPECT_GT(gap_dual, 0.0);       // 非恒等 (口径不同)
    EXPECT_LT(gap_dual, 1e-5);      // 但时间平移等价 (O(1/κ))

    // burn=0 embedded: 含 diffuse 伪贡献 (t=0 的 −½log κ ≈ −7.8) —
    // 与差分样本口径不可混比 (SS5)
    const ss::KalmanResult rb = ss::kalman_filter(y, me, 0);
    ASSERT_TRUE(rb.ok);
    EXPECT_NEAR(rb.loglik, ssb::GOLD_LL_ARIMA111_BURN0, 1e-6);
    EXPECT_NEAR(rb.loglik - re.loglik, -7.826695, 1e-3);  // 伪贡献显式
}

// ---- 8. Monahan 变换锚: constrain 逐位 + 往返 (probe + gen [T]) ----
TEST(ArimaStatespace, MonahanTransformAnchors) {
    for (Size i = 0; i < 4; ++i) {
        const Size n = ssb::MONAHAN_LEN[i];
        std::vector<Real> u(ssb::MONAHAN_IN[i],
                            ssb::MONAHAN_IN[i] + n);
        const std::vector<Real> c = ss::detail::constrain_stationary(u);
        ASSERT_EQ(c.size(), n);
        for (Size j = 0; j < n; ++j)
            EXPECT_NEAR(c[j], ssb::MONAHAN_OUT[i][j], 1e-13)
                << "case " << i << " coef " << j;
        // 往返: unconstrain(constrain(u)) = u (整组联合, 非逐系数)
        const std::vector<Real> back = ss::detail::unconstrain_stationary(c);
        for (Size j = 0; j < n; ++j)
            EXPECT_NEAR(back[j], u[j], 1e-12)
                << "roundtrip case " << i << " coef " << j;
    }
}

// ---- 9. Monahan MA 负号 + unpack/params_to_u 域往返 (probe [C]) ----
TEST(ArimaStatespace, MonahanMaNegativeSign) {
    // MA: θ = −constrain(u) (负号实锤; 小值区 θ ≈ +u — u=−0.2 ⇒ θ≈−0.2)
    // constrain(−0.2) = +0.2/√1.04 (输出带负号: −y) ⇒ θ = −0.2/√1.04
    const auto up = ss::detail::unpack_params({-0.2, 1.0}, 0, 1, true, true,
                                              false);
    ASSERT_EQ(up.theta.size(), 1);
    EXPECT_NEAR(up.theta[0], -0.2 / std::sqrt(1.04), 1e-13);  // −constrain(−0.2)

    // AR: φ = constrain(u) 直传 (statsmodels AR 系数同口径)
    const auto up2 = ss::detail::unpack_params({0.3, 1.0}, 1, 0, true, true,
                                               false);
    EXPECT_NEAR(up2.phi[0], -0.3 / std::sqrt(1.09), 1e-13);

    // σ² = u_σ² (K4 非集中化)
    EXPECT_NEAR(up.sigma2, 1.0, 1e-14);

    // 域往返: params_to_u(phi, θ) 恢复无约束 u (enforce 双开)
    std::vector<Real> u_back;
    ss::detail::params_to_u(up2.phi, up.theta, true, true, u_back);
    ASSERT_EQ(u_back.size(), 2);
    EXPECT_NEAR(u_back[0], 0.3, 1e-12);
    EXPECT_NEAR(u_back[1], -0.2, 1e-12);
}

// ---- 10. MLE vs statsmodels: ARMA(1,1) + ARMA(2,2) 高阶 Monahan ----
TEST(ArimaStatespace, MleArmaVsSarimax) {
    const std::vector<Real> y = data120();

    arima_ns::ArimaSpec s11;  // p=1, d=0, q=1
    s11.p = 1; s11.q = 1;
    const ss::ArimaSSResult r11 = ss::arima_statespace_mle(y, s11);
    EXPECT_TRUE(r11.converged) << r11.message;
    // llf/aic 严格 (同一似然面最优); 参数宽松 (SLSQP vs L-BFGS-B 落点差)
    EXPECT_NEAR(r11.loglik, ssb::MLE_ARMA11_LLF, 1e-5);
    EXPECT_NEAR(r11.aic, ssb::MLE_ARMA11_AIC, 1e-4);
    ASSERT_EQ(r11.phi.size(), 1);
    ASSERT_EQ(r11.theta.size(), 1);
    EXPECT_NEAR(r11.phi[0], ssb::MLE_ARMA11_PARAMS[0], 5e-3);
    EXPECT_NEAR(r11.theta[0], ssb::MLE_ARMA11_PARAMS[1], 5e-3);
    EXPECT_NEAR(r11.sigma2, ssb::MLE_ARMA11_PARAMS[2], 5e-2);
    EXPECT_EQ(r11.loglikelihood_burn, 0);
    // K1 常驻互检: d=0 全样本恒等 (锚 GAP_T300_D0 ~1e-9 量级)
    EXPECT_LT(r11.innovations_gap, 1e-5);
    EXPECT_EQ(r11.innovations.size(), y.size());

    // ARMA(2,2): 高阶 Monahan (5 参数)
    arima_ns::ArimaSpec s22;
    s22.p = 2; s22.q = 2;
    const ss::ArimaSSResult r22 = ss::arima_statespace_mle(y, s22);
    EXPECT_TRUE(r22.converged) << r22.message;
    EXPECT_NEAR(r22.loglik, ssb::MLE_ARMA22_LLF, 1e-4);
    EXPECT_NEAR(r22.aic, ssb::MLE_ARMA22_AIC, 1e-3);
    ASSERT_EQ(r22.phi.size(), 2);
    ASSERT_EQ(r22.theta.size(), 2);
    for (Size j = 0; j < 2; ++j) {
        EXPECT_NEAR(r22.phi[j], ssb::MLE_ARMA22_PARAMS[j], 5e-2);
        EXPECT_NEAR(r22.theta[j], ssb::MLE_ARMA22_PARAMS[2 + j], 5e-2);
    }
}

// ---- 11. MLE vs statsmodels: ARIMA(1,1,1) 嵌入差分 ----
TEST(ArimaStatespace, MleArima111VsSarimax) {
    arima_ns::ArimaSpec s;
    s.p = 1; s.d = 1; s.q = 1;
    const ss::ArimaSSResult r = ss::arima_statespace_mle(data120(), s);
    EXPECT_TRUE(r.converged) << r.message;
    EXPECT_NEAR(r.loglik, ssb::MLE_ARIMA111_LLF, 1e-5);
    EXPECT_NEAR(r.aic, ssb::MLE_ARIMA111_AIC, 1e-4);
    ASSERT_EQ(r.phi.size(), 1);
    ASSERT_EQ(r.theta.size(), 1);
    EXPECT_NEAR(r.phi[0], ssb::MLE_ARIMA111_PARAMS[0], 5e-3);
    EXPECT_NEAR(r.theta[0], ssb::MLE_ARIMA111_PARAMS[1], 5e-3);
    EXPECT_NEAR(r.sigma2, ssb::MLE_ARIMA111_PARAMS[2], 5e-2);
    EXPECT_EQ(r.loglikelihood_burn, 1);  // burn = d (开放问题 #3 裁决)
    // K1: embedded burn=d vs innovations 时间平移等价, gap O(1/κ)
    EXPECT_LT(r.innovations_gap, 1e-5);
}

// ---- 12. MLE vs statsmodels: ARIMA(0,1,1)+drift (常数参数化) ----
TEST(ArimaStatespace, MleDriftVsSarimax) {
    arima_ns::ArimaSpec s;
    s.d = 1; s.q = 1; s.include_drift = true;
    const ss::ArimaSSResult r = ss::arima_statespace_mle(data120(), s);
    EXPECT_TRUE(r.converged) << r.message;
    // statsmodels params = [drift, θ, σ²] (trend='c' 在前)
    EXPECT_NEAR(r.loglik, ssb::MLE_ARIMA011D_LLF, 1e-5);
    EXPECT_NEAR(r.aic, ssb::MLE_ARIMA011D_AIC, 1e-4);
    EXPECT_NEAR(r.constant, ssb::MLE_ARIMA011D_PARAMS[0], 2e-3);
    ASSERT_EQ(r.theta.size(), 1);
    EXPECT_NEAR(r.theta[0], ssb::MLE_ARIMA011D_PARAMS[1], 5e-3);
    EXPECT_NEAR(r.sigma2, ssb::MLE_ARIMA011D_PARAMS[2], 5e-2);
    EXPECT_EQ(r.loglikelihood_burn, 1);
    // K1: 中心化口径 (Δy − μ, μ = c/(1−Σφ)) 同参数等价
    EXPECT_LT(r.innovations_gap, 1e-5);
}

// ---- 13. K1 互检 gap: T=300 d=0 恒等 / d=1 时间平移 (回归锚) ----
TEST(ArimaStatespace, CrossCheckInnovationsGap) {
    const std::vector<Real> y300 = data300();
    const std::vector<Real> dy300 = diff1(y300);

    // d=0: KF(steady) vs innovations (零均值) — 数学恒等
    const ss::StateSpaceModel m0 = ss::make_arma_statespace({0.6}, {0.4}, 1.0);
    const ss::KalmanResult k0 = ss::kalman_filter(y300, m0, 0);
    ASSERT_TRUE(k0.ok);
    const Real ll_inv_d0 =
        ss::detail::innovations_loglik_at(y300, {0.6}, {0.4}, 1.0);
    const Real gap_d0 = std::abs(k0.loglik - ll_inv_d0);
    // 锚 GAP_T300_D0 = 8.05e-10 (statsmodels vs 手工 DL); C++ 双路径同量级
    EXPECT_LT(gap_d0, 1e-8);
    EXPECT_NEAR(gap_d0, ssb::GAP_T300_D0, 1e-7);

    // d=1: embedded(burn=1) vs innovations(Δy) — 时间平移等价, O(1/κ)
    const ss::StateSpaceModel m1 =
        ss::make_arima_statespace({0.6}, {0.4}, 1, 1.0, false, false);
    const ss::KalmanResult k1 = ss::kalman_filter(y300, m1, 1);
    ASSERT_TRUE(k1.ok);
    const Real ll_inv_d1 =
        ss::detail::innovations_loglik_at(dy300, {0.6}, {0.4}, 1.0);
    const Real gap_d1 = std::abs(k1.loglik - ll_inv_d1);
    // 锚 GAP_T300_D1 = 7.87e-7 (κ=1e6 有限尾差)
    EXPECT_LT(gap_d1, 1e-5);                       // K1 告警线
    EXPECT_NEAR(gap_d1, ssb::GAP_T300_D1, 1e-6);   // 锚量级
}

// ---- 14. 非平稳 → 全 diffuse 自动切换 (非静默) + 混合初始化结构 ----
TEST(ArimaStatespace, NonstationaryDiffuseFallback) {
    // φ=1.05 爆炸 AR: 混合初始化跳过 → 全 ApproximateDiffuse (κ=1e6)
    const ss::StateSpaceModel m_exp =
        ss::make_arima_statespace({1.05}, {0.3}, 1, 0.9, false, false);
    EXPECT_EQ(m_exp.init.method, ss::SSInitMethod::ApproximateDiffuse);

    // 全 diffuse KF: burn = d + m₁ = 2 跑通 (丢首 2 期贡献)
    const ss::KalmanResult r_exp =
        ss::kalman_filter(data120(), m_exp, 2);
    EXPECT_TRUE(r_exp.ok);
    EXPECT_EQ(r_exp.loglik_contribs.size(),
              data120().size() - 2);

    // φ=0.5 平稳: 混合初始化结构断言 (Known 注入)
    //   P₀ = diag(κ, P_ARMA), a₀ = (0, a₀^ARMA)
    const ss::StateSpaceModel m_st =
        ss::make_arima_statespace({0.5}, {0.3}, 1, 0.9, false, false);
    ASSERT_EQ(m_st.init.method, ss::SSInitMethod::Known);
    ASSERT_EQ(m_st.init.P0.rows(), 3);
    EXPECT_NEAR(m_st.init.P0(0, 0), 1e6, 1.0);      // 差分块 diffuse
    EXPECT_NEAR(m_st.init.P0(0, 1), 0.0, 1e-14);    // 块间零交叉
    EXPECT_GT(m_st.init.P0(1, 1), 0.0);             // ARMA 块 Lyapunov 正定
    ASSERT_EQ(m_st.init.a0.size(), 3);
    EXPECT_NEAR(m_st.init.a0(0), 0.0, 1e-14);       // 差分块零初始化

    // arima_statespace_mle 非平稳输入: 随机游走数据 d=1 差分后平稳,
    // 互检正常路径 (message 无 WARNING) — 全 diffuse MLE 分支经
    // enforce_stationarity=false 直测:
    arima_ns::ArimaSpec s;
    s.p = 1; s.d = 1; s.q = 0;
    ss::ArimaSSConfig cfg;
    cfg.enforce_stationarity = false;  // 允许非平稳域
    const ss::ArimaSSResult r = ss::arima_statespace_mle(data120(), s, cfg);
    // 收敛或平稳落点均可; 关键: 无 throw + message 协议有效
    if (r.innovations_gap > 1e-5) {
        EXPECT_NE(r.message.find("cross-check"), std::string::npos);
    }
}
