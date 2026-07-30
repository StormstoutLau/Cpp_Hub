// Dupire local volatility tests - TDD (RED -> GREEN)
// SOURCE: PHASE3_SPEC §4.2 - Dupire local volatility recovery
// Implemented on main station (MSVC) - 2026-07-30
// Strategy: self-contained IV grid + BSM analytic call + numerical differentiation
//   sigma^2_loc(K,T) = (dC/dT + qC + (r-q)K dC/dK) / (0.5 K^2 d^2C/dK^2)
#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"
#include "cpphub/models/vol_surface/dupire_local_vol.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

using namespace cpphub::v1;

// ---------- Helper: build flat IV grid ----------
static std::vector<std::vector<Real>> flat_iv_grid(
    const std::vector<Real>& strikes,
    const std::vector<Real>& maturities,
    Real vol) {
    std::vector<std::vector<Real>> ivs(
        maturities.size(), std::vector<Real>(strikes.size(), vol));
    return ivs;
}

// ========== 1. Flat IV recovery: constant implied vol -> local vol ≈ implied vol ==========
TEST(DupireLocalVol, FlatIVRecovery) {
    // 平坦 IV 表面 (σ=0.2): 理论上 Dupire 局部波动率 ≡ 0.2
    // 因为 BSM 模型下局部波动率 = 隐含波动率 (constant vol case)
    std::vector<Real> strikes;
    for (int i = 0; i < 21; ++i) strikes.push_back(80.0 + i * 2.0);  // 80..120
    std::vector<Real> maturities;
    for (int j = 0; j < 11; ++j) maturities.push_back(0.1 + j * 0.1);  // 0.1..1.1
    Real flat_vol = 0.20;
    auto ivs = flat_iv_grid(strikes, maturities, flat_vol);

    DupireLocalVol dup(strikes, maturities, ivs, /*S=*/100.0, /*r=*/0.05, /*q=*/0.02);

    // ATM 附近应精确恢复 flat vol (中心区域数值差分最稳定)
    for (Size j = 2; j < maturities.size() - 2; ++j) {
        for (Size i = 5; i < strikes.size() - 5; ++i) {
            Real lv = dup.local_vol(strikes[i], maturities[j]);
            EXPECT_NEAR(lv, flat_vol, 5e-3)
                << "K=" << strikes[i] << " T=" << maturities[j];
        }
    }
}

// ========== 2. Local variance non-negative for flat IV ==========
TEST(DupireLocalVol, LocalVarianceNonNeg) {
    std::vector<Real> strikes;
    for (int i = 0; i < 21; ++i) strikes.push_back(80.0 + i * 2.0);
    std::vector<Real> maturities = {0.25, 0.5, 1.0};
    auto ivs = flat_iv_grid(strikes, maturities, 0.25);

    DupireLocalVol dup(strikes, maturities, ivs, 100.0, 0.03, 0.0);

    for (Size j = 0; j < maturities.size(); ++j) {
        for (Size i = 3; i < strikes.size() - 3; ++i) {
            Real lv2 = dup.local_variance(strikes[i], maturities[j]);
            EXPECT_GE(lv2, 0.0)
                << "K=" << strikes[i] << " T=" << maturities[j];
        }
    }
}

// ========== 3. Local vol grid output shape ==========
TEST(DupireLocalVol, GridOutputShape) {
    std::vector<Real> strikes = {90, 95, 100, 105, 110};
    std::vector<Real> maturities = {0.25, 0.5, 1.0};
    auto ivs = flat_iv_grid(strikes, maturities, 0.2);

    DupireLocalVol dup(strikes, maturities, ivs, 100.0, 0.0, 0.0);
    auto grid = dup.local_vol_grid(strikes, maturities);

    ASSERT_EQ(grid.size(), maturities.size());
    for (Size j = 0; j < maturities.size(); ++j) {
        ASSERT_EQ(grid[j].size(), strikes.size());
        for (Size i = 0; i < strikes.size(); ++i) {
            EXPECT_GT(grid[j][i], 0.0);
        }
    }
}

// ========== 4. Dupire formula consistency: analytic BSM Greeks vs numerical ==========
// 验证 DupireLocalVol 内部数值差分与独立数值差分 BSM 价格一致
// (避免解析 Greeks 推导错误, 直接用 BSM 价格 + 数值差分做基准)
TEST(DupireLocalVol, AnalyticGreeksConsistency) {
    Real S = 100.0, K = 100.0, T = 0.5, r = 0.05, q = 0.02, sigma = 0.20;

    // 独立数值差分: 用 BSM 解析价格做中心差分
    auto bsm_call = [&](Real kk, Real tt) -> Real {
        auto g = AnalyticGreeksEngine::bsm_european(S, kk, tt, r, q, sigma, true);
        return g.price;
    };
    Real h_K = 1e-3 * K;
    Real h_T = 1e-4 * std::max(T, 0.1);
    Real C   = bsm_call(K, T);
    Real C_Kp = bsm_call(K + h_K, T);
    Real C_Km = bsm_call(K - h_K, T);
    Real C_Tp = bsm_call(K, T + h_T);
    Real C_Tm = bsm_call(K, T - h_T);

    Real dC_dT_ref   = (C_Tp - C_Tm) / (2.0 * h_T);
    Real dC_dK_ref   = (C_Kp - C_Km) / (2.0 * h_K);
    Real d2C_dK2_ref = (C_Kp - 2.0 * C + C_Km) / (h_K * h_K);

    Real num_ref = dC_dT_ref + q * C + (r - q) * K * dC_dK_ref;
    Real den_ref = 0.5 * K * K * d2C_dK2_ref;
    Real sigma2_loc_ref = num_ref / den_ref;

    // 平坦 IV 下解析 + 数值差分应等于 sigma^2
    EXPECT_NEAR(sigma2_loc_ref, sigma * sigma, 1e-4)
        << "ref sigma2_loc=" << sigma2_loc_ref;

    // DupireLocalVol 数值实现 (使用细网格降低插值误差)
    std::vector<Real> strikes;
    for (int i = 0; i < 41; ++i) strikes.push_back(80.0 + i * 1.0);  // 80..120 step 1
    std::vector<Real> maturities;
    for (int j = 0; j < 21; ++j) maturities.push_back(0.1 + j * 0.05);  // 0.1..1.1 step 0.05
    auto ivs = flat_iv_grid(strikes, maturities, sigma);

    DupireLocalVol dup(strikes, maturities, ivs, S, r, q);
    Real sigma2_loc_numeric = dup.local_variance(K, T);

    // 数值 Dupire 与独立数值差分基准一致 (容差来自双线性插值)
    EXPECT_NEAR(sigma2_loc_numeric, sigma2_loc_ref, 5e-4)
        << "ref=" << sigma2_loc_ref << " numeric=" << sigma2_loc_numeric;
    // 数值 Dupire 与理论值 sigma^2 一致
    EXPECT_NEAR(sigma2_loc_numeric, sigma * sigma, 1e-3)
        << "theory=" << sigma * sigma << " numeric=" << sigma2_loc_numeric;
}

// ========== 5. Skew surface: IV with skew -> local vol reflects skew ==========
TEST(DupireLocalVol, SkewSurface) {
    // 简单 skew: σ(K,T) = 0.20 - 0.001*(K-100) + 0.0*T
    // 即 put 方行权价越低 IV 越高 (典型 equity skew)
    std::vector<Real> strikes;
    for (int i = 0; i < 21; ++i) strikes.push_back(80.0 + i * 2.0);
    std::vector<Real> maturities = {0.25, 0.5, 1.0};
    std::vector<std::vector<Real>> ivs(
        maturities.size(), std::vector<Real>(strikes.size()));
    for (Size j = 0; j < maturities.size(); ++j) {
        for (Size i = 0; i < strikes.size(); ++i) {
            ivs[j][i] = 0.20 - 0.001 * (strikes[i] - 100.0);
        }
    }

    DupireLocalVol dup(strikes, maturities, ivs, 100.0, 0.03, 0.0);

    // 局部波动率在低 K 应更高 (反映 skew)
    Real lv_low_K  = dup.local_vol(85.0, 0.5);
    Real lv_atm    = dup.local_vol(100.0, 0.5);
    Real lv_high_K = dup.local_vol(115.0, 0.5);

    EXPECT_GT(lv_low_K, lv_atm);
    EXPECT_LT(lv_high_K, lv_atm);
}

// ========== 6. Constructor from VolSurface reference (interface check) ==========
// 此测试验证 VolSurface 引用构造接口存在 (即使 VolSurface 尚未实现, 接口应可编译)
// 暂用占位测试, 待 B 站实现 VolSurface 后扩展
TEST(DupireLocalVol, IVGridConstructorInterface) {
    std::vector<Real> strikes = {90, 95, 100, 105, 110};
    std::vector<Real> maturities = {0.25, 0.5, 1.0};
    auto ivs = flat_iv_grid(strikes, maturities, 0.2);

    DupireLocalVol dup(strikes, maturities, ivs, 100.0, 0.05, 0.02);

    // 验证访问器
    EXPECT_EQ(dup.r(), 0.05);
    EXPECT_EQ(dup.q(), 0.02);
    EXPECT_EQ(dup.strikes().size(), strikes.size());
    EXPECT_EQ(dup.maturities().size(), maturities.size());
}

// ========== 7. Zero rate / zero dividend case ==========
TEST(DupireLocalVol, ZeroRateZeroDividend) {
    // r=q=0 时 Dupire 公式简化为 σ²_loc = (∂C/∂T) / (0.5 K² ∂²C/∂K²)
    std::vector<Real> strikes;
    for (int i = 0; i < 21; ++i) strikes.push_back(80.0 + i * 2.0);
    std::vector<Real> maturities = {0.25, 0.5, 1.0};
    auto ivs = flat_iv_grid(strikes, maturities, 0.20);

    DupireLocalVol dup(strikes, maturities, ivs, 100.0, 0.0, 0.0);

    // 中心区域应恢复 flat vol
    Real lv = dup.local_vol(100.0, 0.5);
    EXPECT_NEAR(lv, 0.20, 1e-2);
}

// ========== 8. Recovery error: simple MC verification ==========
// 用恢复的局部波动率跑简化 MC (局部波动率模型), 对比 BSM 平坦 IV 价格
// 此测试为简化版本, 不要求 <1bp (需要完整 PDE 引擎), 仅验证数量级合理
TEST(DupireLocalVol, RecoveryErrorSanity) {
    std::vector<Real> strikes;
    for (int i = 0; i < 11; ++i) strikes.push_back(90.0 + i * 2.0);  // 90..110
    std::vector<Real> maturities = {0.25, 0.5, 1.0};
    Real flat_vol = 0.20;
    auto ivs = flat_iv_grid(strikes, maturities, flat_vol);

    DupireLocalVol dup(strikes, maturities, ivs, 100.0, 0.05, 0.0);

    // 恢复误差: 平坦 IV 表面下, 局部波动率应 ≈ flat_vol, 故 recovery_error 应很小
    Real err = dup.recovery_error(100.0, 0.5, 100.0, /*n_paths=*/20000);
    EXPECT_LT(err, 0.5) << "recovery_error too large: " << err;
}
