// v1.2 Batch 5: PFE (Potential Future Exposure) + SA-CCR 单元测试
// 覆盖: ExposureStats (EE/EPE/EEE/EEPE/PFE) + compute_exposure_stats + compute_pfe_mc +
//       SA-CCR (RC/PFE_addon/EAD/NGR/MF/supervisory factor)
#include <gtest/gtest.h>
#include "cpphub/risk/pfe_sa_ccr.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"  // bsm_call_price
#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>

using namespace cpphub;

namespace {
ZeroCurve flat_zero_curve(Real rate, Real max_T = 30.0) {
    std::vector<Real> Ts, rs;
    for (Size i = 1; i <= 30; ++i) {
        Real T = static_cast<Real>(i);
        if (T > max_T) break;
        Ts.push_back(T);
        rs.push_back(rate);
    }
    return ZeroCurve(Ts, rs, ZeroCurve::InterpType::LinearZero);
}
}  // namespace

// ============================================================
// compute_exposure_stats: 基础统计量测试
// ============================================================
TEST(ExposureStatsTest, ConstantPositiveExposure) {
    // 所有路径 V(t) = 1.0 (常数正值), 折现到 t=0
    // EE = 1.0, ENE = 0, PFE = 1.0
    auto discount = flat_zero_curve(0.0);
    std::vector<Real> times = {0.0, 1.0, 2.0, 3.0};
    std::vector<std::vector<Real>> V_samples(100, std::vector<Real>(4, 1.0));
    PFEConfig cfg;
    cfg.confidence_level = 0.99;
    auto stats = compute_exposure_stats(V_samples, times, discount, cfg);
    EXPECT_NEAR(stats.ee[0], 1.0, 1e-12);
    EXPECT_NEAR(stats.ene[0], 0.0, 1e-12);
    EXPECT_NEAR(stats.pfe[0], 1.0, 1e-12);
    EXPECT_NEAR(stats.pfe[3], 1.0, 1e-12);
    EXPECT_NEAR(stats.max_ee, 1.0, 1e-12);
    EXPECT_NEAR(stats.max_pfe, 1.0, 1e-12);
    EXPECT_GT(stats.epe, 0.0);
}

TEST(ExposureStatsTest, ConstantNegativeExposure) {
    // 所有路径 V(t) = -1.0 (负值), EE=0, ENE=1
    auto discount = flat_zero_curve(0.0);
    std::vector<Real> times = {0.0, 1.0, 2.0};
    std::vector<std::vector<Real>> V_samples(100, std::vector<Real>(3, -1.0));
    auto stats = compute_exposure_stats(V_samples, times, discount);
    EXPECT_NEAR(stats.ee[0], 0.0, 1e-12);
    EXPECT_NEAR(stats.ene[0], 1.0, 1e-12);
    EXPECT_NEAR(stats.pfe[0], 0.0, 1e-12);
    EXPECT_NEAR(stats.epe, 0.0, 1e-12);
}

TEST(ExposureStatsTest, SymmetricExposure) {
    // 50% 路径 V=1, 50% 路径 V=-1 → EE=0.5, ENE=0.5
    auto discount = flat_zero_curve(0.0);
    std::vector<Real> times = {0.0};
    std::vector<std::vector<Real>> V_samples(100, std::vector<Real>(1));
    for (Size i = 0; i < 100; ++i) {
        V_samples[i][0] = (i < 50) ? 1.0 : -1.0;
    }
    auto stats = compute_exposure_stats(V_samples, times, discount);
    EXPECT_NEAR(stats.ee[0], 0.5, 1e-12);
    EXPECT_NEAR(stats.ene[0], 0.5, 1e-12);
    // PFE 99%: 排序后最大 1% = 全部 1.0
    EXPECT_NEAR(stats.pfe[0], 1.0, 1e-12);
}

TEST(ExposureStatsTest, PFESinglePathConstant) {
    // 单路径 V=2.0, EE=2.0, PFE=2.0
    auto discount = flat_zero_curve(0.0);
    std::vector<Real> times = {0.0, 1.0};
    std::vector<std::vector<Real>> V_samples(1, std::vector<Real>(2, 2.0));
    auto stats = compute_exposure_stats(V_samples, times, discount);
    EXPECT_NEAR(stats.ee[0], 2.0, 1e-12);
    EXPECT_NEAR(stats.pfe[0], 2.0, 1e-12);
}

TEST(ExposureStatsTest, EmptySamplesThrows) {
    auto discount = flat_zero_curve(0.0);
    std::vector<Real> times = {0.0, 1.0};
    std::vector<std::vector<Real>> V_samples;
    EXPECT_THROW(compute_exposure_stats(V_samples, times, discount), std::invalid_argument);
}

TEST(ExposureStatsTest, InvalidConfidenceThrows) {
    auto discount = flat_zero_curve(0.0);
    std::vector<Real> times = {0.0};
    std::vector<std::vector<Real>> V_samples(1, std::vector<Real>(1, 1.0));
    PFEConfig cfg;
    cfg.confidence_level = 0.0;
    EXPECT_THROW(compute_exposure_stats(V_samples, times, discount, cfg), std::invalid_argument);
    cfg.confidence_level = 1.5;
    EXPECT_THROW(compute_exposure_stats(V_samples, times, discount, cfg), std::invalid_argument);
}

TEST(ExposureStatsTest, ColumnSizeMismatchThrows) {
    auto discount = flat_zero_curve(0.0);
    std::vector<Real> times = {0.0, 1.0};
    std::vector<std::vector<Real>> V_samples;
    V_samples.push_back({1.0, 2.0});
    V_samples.push_back({1.0});  // 长度不一致
    EXPECT_THROW(compute_exposure_stats(V_samples, times, discount), std::invalid_argument);
}

// ============================================================
// EEE (Effective EE) running max 测试
// ============================================================
TEST(ExposureStatsTest, EEERunningMax) {
    // EE(t1)=1, EE(t2)=2, EE(t3)=0.5 → EEE(t3) 应 = 2 (running max)
    auto discount = flat_zero_curve(0.0);
    std::vector<Real> times = {0.0, 1.0, 2.0, 3.0};
    Size n_paths = 1000;
    std::vector<std::vector<Real>> V_samples(n_paths, std::vector<Real>(4));
    // t=0: V=0; t=1: V=1 (正); t=2: V=2 (正); t=3: V=0.5 (正, 但比 t=2 小)
    for (Size p = 0; p < n_paths; ++p) {
        V_samples[p][0] = 0.0;
        V_samples[p][1] = 1.0;
        V_samples[p][2] = 2.0;
        V_samples[p][3] = 0.5;
    }
    PFEConfig cfg;
    cfg.use_running_max_ee = true;
    auto stats = compute_exposure_stats(V_samples, times, discount, cfg);
    EXPECT_NEAR(stats.ee[0], 0.0, 1e-12);
    EXPECT_NEAR(stats.ee[1], 1.0, 1e-12);
    EXPECT_NEAR(stats.ee[2], 2.0, 1e-12);
    EXPECT_NEAR(stats.ee[3], 0.5, 1e-12);
    EXPECT_NEAR(stats.eee[0], 0.0, 1e-12);
    EXPECT_NEAR(stats.eee[1], 1.0, 1e-12);
    EXPECT_NEAR(stats.eee[2], 2.0, 1e-12);
    EXPECT_NEAR(stats.eee[3], 2.0, 1e-12);  // running max 不下降
}

TEST(ExposureStatsTest, EEEDisabledWhenFlagOff) {
    auto discount = flat_zero_curve(0.0);
    std::vector<Real> times = {0.0, 1.0, 2.0};
    Size n_paths = 100;
    std::vector<std::vector<Real>> V_samples(n_paths, std::vector<Real>(3));
    for (Size p = 0; p < n_paths; ++p) {
        V_samples[p][0] = 0.0;
        V_samples[p][1] = 5.0;
        V_samples[p][2] = 1.0;
    }
    PFEConfig cfg;
    cfg.use_running_max_ee = false;
    auto stats = compute_exposure_stats(V_samples, times, discount, cfg);
    EXPECT_NEAR(stats.eee[2], 1.0, 1e-12);  // 不取 running max, EEE=EE
    EXPECT_NE(stats.eee[2], stats.eee[1]);
}

// ============================================================
// Discount Effect 测试
// ============================================================
TEST(ExposureStatsTest, DiscountReducesExposure) {
    // V(t)=1, r=5% → EE(t) = 1 * exp(-r*t), 随 t 衰减
    auto discount = flat_zero_curve(0.05);
    std::vector<Real> times = {0.0, 1.0, 5.0, 10.0};
    Size n_paths = 100;
    std::vector<std::vector<Real>> V_samples(n_paths, std::vector<Real>(4, 1.0));
    auto stats = compute_exposure_stats(V_samples, times, discount);
    EXPECT_NEAR(stats.ee[0], 1.0, 1e-12);                          // exp(0) = 1
    EXPECT_NEAR(stats.ee[1], std::exp(-0.05 * 1.0), 1e-12);
    EXPECT_NEAR(stats.ee[2], std::exp(-0.05 * 5.0), 1e-12);
    EXPECT_NEAR(stats.ee[3], std::exp(-0.05 * 10.0), 1e-12);
    // EE 随 t 递减 (折现效应)
    EXPECT_GT(stats.ee[0], stats.ee[1]);
    EXPECT_GT(stats.ee[1], stats.ee[2]);
}

// ============================================================
// PFE 分位数精度测试 (均匀分布)
// ============================================================
TEST(ExposureStatsTest, PFEQuantileUniform) {
    // V(t) ~ U[0, 1] (10000 条路径), PFE_99% 应 ≈ 0.99
    auto discount = flat_zero_curve(0.0);
    std::vector<Real> times = {0.0};
    Size n_paths = 10000;
    std::vector<std::vector<Real>> V_samples(n_paths, std::vector<Real>(1));
    for (Size p = 0; p < n_paths; ++p) {
        V_samples[p][0] = static_cast<Real>(p) / static_cast<Real>(n_paths - 1);
    }
    PFEConfig cfg;
    cfg.confidence_level = 0.99;
    auto stats = compute_exposure_stats(V_samples, times, discount, cfg);
    EXPECT_NEAR(stats.pfe[0], 0.99, 1e-3);
}

TEST(ExposureStatsTest, PFE50Quantile) {
    // α=0.5 应 ≈ 中位数
    auto discount = flat_zero_curve(0.0);
    std::vector<Real> times = {0.0};
    Size n_paths = 1001;  // 奇数, 中位数清晰
    std::vector<std::vector<Real>> V_samples(n_paths, std::vector<Real>(1));
    for (Size p = 0; p < n_paths; ++p) {
        V_samples[p][0] = static_cast<Real>(p);  // V = 0, 1, 2, ..., 1000
    }
    PFEConfig cfg;
    cfg.confidence_level = 0.5;
    auto stats = compute_exposure_stats(V_samples, times, discount, cfg);
    EXPECT_NEAR(stats.pfe[0], 500.0, 1.0);
}

// ============================================================
// compute_pfe_mc 端到端 (vanilla call 的暴露)
// ============================================================
TEST(PFEMCTest, VanillaCallExposure) {
    // BSM call 的暴露: V(t) = call(t, S(t))
    // 在 ATM 附近, EE(t) 应先上升后下降 (典型 PFE profile)
    auto discount = flat_zero_curve(0.05);
    Real S0 = 100.0, sigma = 0.20, T = 1.0, K = 100.0, r = 0.05;
    auto cfg_gbm = make_single_asset_gbm(S0, sigma, r, 0.0, T, 50);
    MultiAssetGBMPathGenerator gen(cfg_gbm);
    std::vector<Real> exposure_times = {0.0, 0.25, 0.5, 0.75, 1.0};

    auto value_fn = [K, r, sigma, T](Real t, const std::vector<Real>& S) -> Real {
        if (S.empty()) return 0.0;
        Real tau = T - t;
        if (tau <= 0.0) return std::max(S[0] - K, 0.0);
        return bsm_call_price(S[0], K, r, 0.0, sigma, tau);
    };

    PFEConfig cfg;
    cfg.confidence_level = 0.95;
    auto stats = compute_pfe_mc(gen, value_fn, exposure_times, discount, 5000, 42, cfg);
    EXPECT_EQ(stats.times.size(), 5u);
    EXPECT_GE(stats.ee[0], 0.0);
    // t=0 的 call value 应 ≈ BSM call price
    Real bs_price = bsm_call_price(S0, K, r, 0.0, sigma, T);
    EXPECT_NEAR(stats.ee[0], bs_price, 0.5);  // MC 误差
    // PFE 应随时间增加 (未来不确定性)
    EXPECT_GE(stats.pfe[3], stats.pfe[1]);
}

TEST(PFEMCTest, ZeroVolGivesConstantExposure) {
    // sigma=0 → 路径确定 S(t)=S0*exp(rt), V(t) 确定
    auto discount = flat_zero_curve(0.05);
    Real S0 = 100.0, sigma = 0.0, T = 1.0, K = 100.0, r = 0.05;
    auto cfg_gbm = make_single_asset_gbm(S0, sigma, r, 0.0, T, 12);
    MultiAssetGBMPathGenerator gen(cfg_gbm);
    std::vector<Real> exposure_times = {0.0, 0.5, 1.0};

    auto value_fn = [K, r, sigma, T](Real t, const std::vector<Real>& S) -> Real {
        Real tau = T - t;
        if (tau <= 0.0) return std::max(S[0] - K, 0.0);
        return bsm_call_price(S[0], K, r, 0.0, sigma, tau);
    };

    auto stats = compute_pfe_mc(gen, value_fn, exposure_times, discount, 100, 42);
    // 零波动率: 所有路径相同, PFE = EE
    EXPECT_NEAR(stats.pfe[0], stats.ee[0], 1e-9);
    EXPECT_NEAR(stats.pfe[1], stats.ee[1], 1e-9);
}

// ============================================================
// SA-CCR: Supervisory Factor 测试
// ============================================================
TEST(SACCRSupervisoryFactorTest, AssetClassFactors) {
    EXPECT_NEAR(supervisory_factor(AssetClass::InterestRate), 0.005, 1e-15);
    EXPECT_NEAR(supervisory_factor(AssetClass::FX), 0.04, 1e-15);
    EXPECT_NEAR(supervisory_factor(AssetClass::Credit), 0.05, 1e-15);
    EXPECT_NEAR(supervisory_factor(AssetClass::Equity), 0.32, 1e-15);
    EXPECT_NEAR(supervisory_factor(AssetClass::Commodity), 0.18, 1e-15);
}

TEST(SACCRSupervisoryFactorTest, CreditByRating) {
    EXPECT_NEAR(credit_supervisory_factor_by_rating("AAA"), 0.0038, 1e-15);
    EXPECT_NEAR(credit_supervisory_factor_by_rating("AA"), 0.0038, 1e-15);
    EXPECT_NEAR(credit_supervisory_factor_by_rating("A"), 0.0038, 1e-15);
    EXPECT_NEAR(credit_supervisory_factor_by_rating("BBB"), 0.01, 1e-15);
    EXPECT_NEAR(credit_supervisory_factor_by_rating("BB"), 0.06, 1e-15);
    EXPECT_NEAR(credit_supervisory_factor_by_rating("B"), 0.10, 1e-15);
    EXPECT_NEAR(credit_supervisory_factor_by_rating("CCC"), 0.15, 1e-15);
    EXPECT_NEAR(credit_supervisory_factor_by_rating(""), 0.05, 1e-15);
}

TEST(SACCRSupervisoryFactorTest, EquityIndexVsSingleName) {
    EXPECT_NEAR(equity_supervisory_factor(true), 0.12, 1e-15);
    EXPECT_NEAR(equity_supervisory_factor(false), 0.32, 1e-15);
}

// ============================================================
// SA-CCR: Maturity Factor (MF) 测试
// ============================================================
TEST(SACCRMaturityFactorTest, ZeroMaturity) {
    EXPECT_NEAR(maturity_factor(0.0), 0.0, 1e-15);
    EXPECT_NEAR(maturity_factor(-1.0), 0.0, 1e-15);
}

TEST(SACCRMaturityFactorTest, BelowOneYearSqrt) {
    EXPECT_NEAR(maturity_factor(0.25), std::sqrt(0.25), 1e-15);  // = 0.5
    EXPECT_NEAR(maturity_factor(0.5), std::sqrt(0.5), 1e-15);
    EXPECT_NEAR(maturity_factor(0.75), std::sqrt(0.75), 1e-15);
}

TEST(SACCRMaturityFactorTest, AtOrAboveOneYear) {
    EXPECT_NEAR(maturity_factor(1.0), 1.0, 1e-15);
    EXPECT_NEAR(maturity_factor(5.0), 1.0, 1e-15);
    EXPECT_NEAR(maturity_factor(30.0), 1.0, 1e-15);
}

// ============================================================
// SA-CCR: 单笔交易测试
// ============================================================
TEST(SACCRTest, SingleIRSTrade) {
    // 单笔 IRS, notional=10M, 5Y → addon = 10M × 0.005 × 1 = 50000
    // PFE_addon = 1.4 × 50000 = 70000
    // RC = 0 (假设 V=0, C=0)
    // EAD = 1.4 × (0 + 70000) = 98000
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_irs_trade(10'000'000.0, 5.0));
    portfolio.portfolio_value = 0.0;
    portfolio.net_collateral = 0.0;

    auto r = compute_sa_ccr(portfolio);
    EXPECT_NEAR(r.rc, 0.0, 1e-9);
    EXPECT_NEAR(r.aggregate_addon, 10'000'000.0 * 0.005, 1e-6);
    EXPECT_NEAR(r.pfe_addon, 1.4 * 10'000'000.0 * 0.005, 1e-6);
    EXPECT_NEAR(r.ead, 1.4 * (0.0 + 1.4 * 10'000'000.0 * 0.005), 1e-6);
    EXPECT_NEAR(r.ngr, 1.0, 1e-9);  // 单笔, NGR=1
}

TEST(SACCRTest, SingleFXTrade) {
    // 单笔 FX, notional=1M, 1Y → addon = 1M × 0.04 × 1 = 40000
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_fx_trade(1'000'000.0, 1.0));
    portfolio.portfolio_value = 0.0;
    auto r = compute_sa_ccr(portfolio);
    EXPECT_NEAR(r.aggregate_addon, 1'000'000.0 * 0.04, 1e-6);
    EXPECT_NEAR(r.pfe_addon, 1.4 * 1'000'000.0 * 0.04, 1e-6);
}

TEST(SACCRTest, SingleCreditTradeAAA) {
    // AAA 信用, notional=5M, 3Y → addon = 5M × 0.0038 × 1 = 19000
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_credit_trade(5'000'000.0, 3.0, "AAA"));
    portfolio.portfolio_value = 0.0;
    auto r = compute_sa_ccr(portfolio);
    EXPECT_NEAR(r.aggregate_addon, 5'000'000.0 * 0.0038, 1e-6);
}

TEST(SACCRTest, SingleEquityTradeIndex) {
    // Index 股票, notional=2M, 2Y → addon = 2M × 0.12 × 1 = 240000
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_equity_trade(2'000'000.0, 2.0, true));
    portfolio.portfolio_value = 0.0;
    auto r = compute_sa_ccr(portfolio);
    EXPECT_NEAR(r.aggregate_addon, 2'000'000.0 * 0.12, 1e-6);
}

TEST(SACCRTest, SingleCommodityTrade) {
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_commodity_trade(1'000'000.0, 1.0));
    portfolio.portfolio_value = 0.0;
    auto r = compute_sa_ccr(portfolio);
    EXPECT_NEAR(r.aggregate_addon, 1'000'000.0 * 0.18, 1e-6);
}

// ============================================================
// SA-CCR: Replacement Cost (RC) 测试
// ============================================================
TEST(SACCRTest, PositivePortfolioValueGivesRC) {
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 5.0));
    portfolio.portfolio_value = 100'000.0;   // 组合有正市值
    portfolio.net_collateral = 0.0;
    auto r = compute_sa_ccr(portfolio);
    EXPECT_NEAR(r.rc, 100'000.0, 1e-9);
}

TEST(SACCRTest, NegativePortfolioValueZeroRC) {
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 5.0));
    portfolio.portfolio_value = -200'000.0;   // 负市值
    portfolio.net_collateral = 0.0;
    auto r = compute_sa_ccr(portfolio);
    EXPECT_NEAR(r.rc, 0.0, 1e-9);  // RC = max(V-C, 0) = 0
}

TEST(SACCRTest, CollateralReducesRC) {
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 5.0));
    portfolio.portfolio_value = 100'000.0;
    portfolio.net_collateral = 60'000.0;   // 银行已收抵押
    auto r = compute_sa_ccr(portfolio);
    EXPECT_NEAR(r.rc, 40'000.0, 1e-9);
}

TEST(SACCRTest, ExcessCollateralNoNegativeRC) {
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 5.0));
    portfolio.portfolio_value = 100'000.0;
    portfolio.net_collateral = 200'000.0;   // 抵押 > V
    auto r = compute_sa_ccr(portfolio);
    EXPECT_NEAR(r.rc, 0.0, 1e-9);  // 不能为负
}

// ============================================================
// SA-CCR: Maturity Factor 应用测试
// ============================================================
TEST(SACCRTest, ShortMaturityReducesAddon) {
    // 同样 notional, 0.25Y vs 1Y → 0.25Y addon = 0.5 × 1Y addon
    SACCRPortfolio p_short, p_long;
    p_short.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 0.25));
    p_long.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 1.0));
    auto r_short = compute_sa_ccr(p_short);
    auto r_long = compute_sa_ccr(p_long);
    EXPECT_NEAR(r_short.aggregate_addon, 0.5 * r_long.aggregate_addon, 1e-6);
}

TEST(SACCRTest, AboveOneYearMaturityCappedAtOne) {
    // 5Y 和 10Y MF 都 = 1.0
    SACCRPortfolio p5, p10;
    p5.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 5.0));
    p10.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 10.0));
    auto r5 = compute_sa_ccr(p5);
    auto r10 = compute_sa_ccr(p10);
    EXPECT_NEAR(r5.aggregate_addon, r10.aggregate_addon, 1e-9);
}

// ============================================================
// SA-CCR: Netting 测试 (long-short 抵消)
// ============================================================
TEST(SACCRTest, OppositeDirectionReducesAddon) {
    // 两笔 IRS, 一多一空, 完全抵消 → net=0, 但 supervisory floor = 0.4 × gross
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 5.0, 1.0));
    portfolio.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 5.0, -1.0));
    portfolio.portfolio_value = 0.0;
    auto r = compute_sa_ccr(portfolio);
    // L=1M, S=1M, net=0, gross=2M, addon = 0.005 × max(0, 0.4×2M) = 0.005 × 800000 = 4000
    Real expected_addon = 0.005 * 0.4 * 2'000'000.0;
    EXPECT_NEAR(r.aggregate_addon, expected_addon, 1e-6);
}

TEST(SACCRTest, SameDirectionAddsNotional) {
    // 两笔 IRS 同方向 → net = 2M, gross = 2M, addon = 0.005 × max(2M, 0.4×2M) = 0.005 × 2M = 10000
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 5.0, 1.0));
    portfolio.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 5.0, 1.0));
    portfolio.portfolio_value = 0.0;
    auto r = compute_sa_ccr(portfolio);
    Real expected_addon = 0.005 * 2'000'000.0;
    EXPECT_NEAR(r.aggregate_addon, expected_addon, 1e-6);
}

TEST(SACCRTest, NGROneForSingleDirection) {
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 5.0, 1.0));
    EXPECT_NEAR(portfolio.compute_ngr(), 1.0, 1e-9);
}

TEST(SACCRTest, NGRZeroForPerfectOffset) {
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 5.0, 1.0));
    portfolio.trades.push_back(make_sa_ccr_irs_trade(1'000'000.0, 5.0, -1.0));
    // 同 asset class 完全抵消 → net=0, gross=2M → NGR=0
    EXPECT_NEAR(portfolio.compute_ngr(), 0.0, 1e-9);
}

// ============================================================
// SA-CCR: 多资产类别组合测试
// ============================================================
TEST(SACCRTest, MultiAssetClassPortfolio) {
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_irs_trade(10'000'000.0, 5.0));       // IR
    portfolio.trades.push_back(make_sa_ccr_fx_trade(2'000'000.0, 1.0));         // FX
    portfolio.trades.push_back(make_sa_ccr_equity_trade(1'000'000.0, 2.0, false)); // Equity
    portfolio.portfolio_value = 0.0;
    auto r = compute_sa_ccr(portfolio);
    // 各类 addon 单独计算后求和 (无 netting 跨类)
    Real ir_addon = 0.005 * 10'000'000.0;
    Real fx_addon = 0.04 * 2'000'000.0;
    Real eq_addon = 0.32 * 1'000'000.0;
    Real expected = ir_addon + fx_addon + eq_addon;
    EXPECT_NEAR(r.aggregate_addon, expected, 1e-6);
    EXPECT_NEAR(r.pfe_addon, 1.4 * expected, 1e-6);
    EXPECT_NEAR(r.ead, 1.4 * 1.4 * expected, 1e-6);  // RC=0
}

TEST(SACCRTest, EADFormula) {
    // EAD = 1.4 × (RC + PFE_addon)
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_irs_trade(5'000'000.0, 3.0));
    portfolio.portfolio_value = 50'000.0;
    portfolio.net_collateral = 10'000.0;
    auto r = compute_sa_ccr(portfolio);
    Real expected_rc = 40'000.0;  // max(50000 - 10000, 0)
    Real expected_ead = 1.4 * (expected_rc + r.pfe_addon);
    EXPECT_NEAR(r.ead, expected_ead, 1e-6);
}

TEST(SACCRTest, EmptyPortfolio) {
    SACCRPortfolio portfolio;
    auto r = compute_sa_ccr(portfolio);
    EXPECT_NEAR(r.rc, 0.0, 1e-15);
    EXPECT_NEAR(r.aggregate_addon, 0.0, 1e-15);
    EXPECT_NEAR(r.pfe_addon, 0.0, 1e-15);
    EXPECT_NEAR(r.ead, 0.0, 1e-15);
    EXPECT_NEAR(r.ngr, 0.0, 1e-15);
}

// ============================================================
// SA-CCR: Rating 影响 addon 测试
// ============================================================
TEST(SACCRTest, LowerRatingIncreasesAddon) {
    // 同样 notional, AAA vs CCC → CCC addon 远大于 AAA
    SACCRPortfolio p_aaa, p_ccc;
    p_aaa.trades.push_back(make_sa_ccr_credit_trade(1'000'000.0, 5.0, "AAA"));
    p_ccc.trades.push_back(make_sa_ccr_credit_trade(1'000'000.0, 5.0, "CCC"));
    auto r_aaa = compute_sa_ccr(p_aaa);
    auto r_ccc = compute_sa_ccr(p_ccc);
    EXPECT_GT(r_ccc.aggregate_addon, r_aaa.aggregate_addon);
    // 比例应 = 0.15 / 0.0038
    Real ratio = r_ccc.aggregate_addon / r_aaa.aggregate_addon;
    EXPECT_NEAR(ratio, 0.15 / 0.0038, 0.01);
}

// ============================================================
// SA-CCR vs MC PFE 概念对照 (一致性验证)
// ============================================================
TEST(SACCROverlapTest, SingleIRSExposureScale) {
    // 单笔 IRS notional=10M, 5Y, SA-CCR EAD 与 MC PFE 99% 数量级一致
    // SA-CCR 通常更保守 (监管公式)
    SACCRPortfolio portfolio;
    portfolio.trades.push_back(make_sa_ccr_irs_trade(10'000'000.0, 5.0));
    portfolio.portfolio_value = 0.0;
    auto r = compute_sa_ccr(portfolio);
    Real sa_ccr_ead = r.ead;

    // IRS 的 MC PFE 99% 1Y 一般在 notional 的 1%-5% 量级
    // 这里仅检查 SA-CCR EAD 在合理范围 (notional 的 1%-10%)
    Real pct = sa_ccr_ead / 10'000'000.0;
    EXPECT_GT(pct, 0.005);  // > 0.5%
    EXPECT_LT(pct, 0.10);   // < 10%
}
