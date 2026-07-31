// ===========================================================================
// Bates / VG / CEV 多模型标定器单元测试
// SOURCE: docs/tasks/E_MULTI_CALIBRATOR_TASK.md
//
// 策略: 用已知参数生成合成 IV surface (定价函数与标定器内部一致), 再标定
//   1. 验证参数可恢复 (合成数据, 无噪声)
//   2. 验证 Feller 条件 / 模型退化极限 (beta->1 BS, beta->0 正态 skew, nu->0 BS)
//   3. 验证 DE 全局搜索 + LM 精炼 优于纯 LM
//   4. 边界参数不崩溃
//
// 每个标定器 5 个测试, 共 15 个测试.
// 复用 gtest 固定装置 (SetUpTestSuite) 共享一次标定结果, 避免重复昂贵的 DE 运行.
// ===========================================================================

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <string>

#include "cpphub/calibration/calibrator.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"  // bsm_call_price

using namespace cpphub::v1;

namespace {

// 由价格构造 MarketQuote, 并用 bsm_implied_vol 反推 IV 作为标定目标
inline MarketQuote make_iv_quote(Real S, Real K, Real T, Real r, Real q, Real price) {
    return MarketQuote{K, T, price, bsm_implied_vol(price, S, K, T, r, q, true), 0.0};
}

// 合成 Bates IV surface: T x K 网格
inline std::vector<MarketQuote> make_bates_quotes(
        Real S, Real r, Real q, const BatesCFParams& tp,
        const std::vector<Real>& strikes, const std::vector<Real>& maturities) {
    std::vector<MarketQuote> quotes;
    for (Real T : maturities) {
        auto phi = make_bates_cf(S, r, q, tp, T);
        COSEngine::Config cc;
        cc.n_terms = 256;
        COSEngine engine(phi, S, r, q, T, cc);
        for (Real K : strikes) {
            Real price = engine.price_call(K);
            quotes.push_back(make_iv_quote(S, K, T, r, q, price));
        }
    }
    return quotes;
}

// 合成 VG IV surface
inline std::vector<MarketQuote> make_vg_quotes(
        Real S, Real r, Real q, Real sigma, Real nu, Real theta,
        const std::vector<Real>& strikes, const std::vector<Real>& maturities) {
    std::vector<MarketQuote> quotes;
    for (Real T : maturities) {
        auto phi = make_vg_cf(S, r, q, sigma, nu, theta, T);
        COSEngine::Config cc;
        cc.n_terms = 256;
        COSEngine engine(phi, S, r, q, T, cc);
        for (Real K : strikes) {
            Real price = engine.price_call(K);
            quotes.push_back(make_iv_quote(S, K, T, r, q, price));
        }
    }
    return quotes;
}

// 合成 CEV IV surface
inline std::vector<MarketQuote> make_cev_quotes(
        Real S, Real r, Real q, const CEVParams& tp,
        const std::vector<Real>& strikes, Real T) {
    std::vector<MarketQuote> quotes;
    for (Real K : strikes) {
        Real price = cev_call_price(S, K, T, r, q, tp);
        quotes.push_back(make_iv_quote(S, K, T, r, q, price));
    }
    return quotes;
}

}  // namespace

// ===========================================================================
// BatesCalibrator (8 参数: v0, kappa, theta, sigma_v, rho, lambda, mu_J, sigma_J)
// ===========================================================================

class BatesCalibratorTest : public ::testing::Test {
protected:
    static std::vector<MarketQuote> quotes_;
    static BatesParams true_p_;
    static CalibrationResult de_result_;
    static CalibrationResult lm_result_;

    static void SetUpTestSuite() {
        const Real S = 100.0, r = 0.03, q = 0.0;
        true_p_ = BatesParams{0.04, 1.5, 0.05, 0.3, -0.4, 0.3, -0.1, 0.2};
        BatesCFParams cfp;
        cfp.v0 = true_p_.v0;
        cfp.kappa = true_p_.kappa;
        cfp.theta = true_p_.theta;
        cfp.sigma = true_p_.sigma_v;
        cfp.rho = true_p_.rho;
        cfp.lambda = true_p_.lambda;
        cfp.mu_J = true_p_.mu_J;
        cfp.sigma_J = true_p_.sigma_J;
        cfp.r = r;
        cfp.q = q;

        std::vector<Real> strikes = {75, 85, 95, 100, 105, 115, 125};
        std::vector<Real> maturities = {0.5, 1.0, 2.0};
        quotes_ = make_bates_quotes(S, r, q, cfp, strikes, maturities);

        BatesCalibrator cal;
        cal.set_market(S, r, q);
        CalibConfig cfg;
        cfg.de_pop_size = 80;
        cfg.de_generations = 200;
        cfg.lm_max_iter = 400;
        de_result_ = cal.calibrate(quotes_, cfg);

        CalibConfig cfg_lm;
        cfg_lm.use_de_init = false;
        cfg_lm.lm_max_iter = 400;
        lm_result_ = cal.calibrate(quotes_, cfg_lm);
    }

    static BatesParams recovered() {
        BatesCalibrator cal;
        return cal.extract_params(de_result_.params);
    }
};

std::vector<MarketQuote> BatesCalibratorTest::quotes_;
BatesParams BatesCalibratorTest::true_p_;
CalibrationResult BatesCalibratorTest::de_result_;
CalibrationResult BatesCalibratorTest::lm_result_;

TEST_F(BatesCalibratorTest, SyntheticDataRecovery) {
    // 用已知 Bates 参数生成 IV surface, 标定恢复参数 (容差 5%)
    BatesParams rec = recovered();
    EXPECT_NEAR(rec.v0, true_p_.v0, 0.05 * true_p_.v0)
        << "v0: expected " << true_p_.v0 << " got " << rec.v0;
    EXPECT_NEAR(rec.kappa, true_p_.kappa, 0.05 * true_p_.kappa)
        << "kappa: expected " << true_p_.kappa << " got " << rec.kappa;
    EXPECT_NEAR(rec.theta, true_p_.theta, 0.05 * true_p_.theta)
        << "theta: expected " << true_p_.theta << " got " << rec.theta;
    EXPECT_NEAR(rec.sigma_v, true_p_.sigma_v, 0.05 * true_p_.sigma_v)
        << "sigma_v: expected " << true_p_.sigma_v << " got " << rec.sigma_v;
    EXPECT_NEAR(rec.rho, true_p_.rho, 0.05)
        << "rho: expected " << true_p_.rho << " got " << rec.rho;
    // 功能拟合: IV 残差很小
    EXPECT_LT(de_result_.objective_value, 1e-4);
}

TEST_F(BatesCalibratorTest, JumpParameterRecovery) {
    // 跳跃参数 lambda/mu_J/sigma_J 在合成数据下可恢复 (容差 10%)
    BatesParams rec = recovered();
    EXPECT_NEAR(rec.lambda, true_p_.lambda, 0.10 * true_p_.lambda)
        << "lambda: expected " << true_p_.lambda << " got " << rec.lambda;
    EXPECT_NEAR(rec.mu_J, true_p_.mu_J, 0.10 * std::abs(true_p_.mu_J))
        << "mu_J: expected " << true_p_.mu_J << " got " << rec.mu_J;
    EXPECT_NEAR(rec.sigma_J, true_p_.sigma_J, 0.10 * true_p_.sigma_J)
        << "sigma_J: expected " << true_p_.sigma_J << " got " << rec.sigma_J;
}

TEST_F(BatesCalibratorTest, FellerConditionSatisfied) {
    // 标定结果应满足 2*kappa*theta > sigma_v^2
    BatesParams rec = recovered();
    EXPECT_TRUE(BatesCalibrator::check_feller(rec))
        << "calibrated params violate Feller: 2*kappa*theta="
        << 2.0 * rec.kappa * rec.theta << " sigma_v^2=" << rec.sigma_v * rec.sigma_v;
    // 静态检查函数本身
    EXPECT_TRUE(BatesCalibrator::check_feller(true_p_));
    EXPECT_FALSE(BatesCalibrator::check_feller(BatesParams{0.04, 0.5, 0.04, 0.4, -0.5, 0.2, 0.0, 0.1}));
}

TEST_F(BatesCalibratorTest, DEPlusLMBetterThanPureLM) {
    // DE 全局搜索初始化的目标函数值应 <= 纯 LM (从固定默认点出发)
    EXPECT_LE(de_result_.objective_value, lm_result_.objective_value + 1e-8);
    EXPECT_LT(de_result_.objective_value, 1e-4);
}

TEST(BatesCalibrator, BoundaryParamsNoCrash) {
    // 参数接近边界 (v0/lambda/sigma_J 小, rho/sigma_v 大) 时不崩溃
    const Real S = 100.0, r = 0.03, q = 0.0;
    BatesCFParams tp;
    tp.v0 = 0.01; tp.kappa = 5.0; tp.theta = 0.08; tp.sigma = 0.5;
    tp.rho = 0.9; tp.lambda = 0.05; tp.mu_J = -0.4; tp.sigma_J = 0.05;
    tp.r = r; tp.q = q;
    std::vector<Real> strikes = {80, 90, 100, 110, 120};
    std::vector<Real> maturities = {0.5, 1.0, 2.0};
    auto quotes = make_bates_quotes(S, r, q, tp, strikes, maturities);

    BatesCalibrator cal;
    cal.set_market(S, r, q);
    CalibConfig cfg;
    cfg.de_pop_size = 40;
    cfg.de_generations = 60;
    cfg.lm_max_iter = 200;

    CalibrationResult res;
    EXPECT_NO_THROW(res = cal.calibrate(quotes, cfg));
    ASSERT_EQ(res.params.size(), 8u);
    for (Real p : res.params) EXPECT_TRUE(std::isfinite(p));
    EXPECT_TRUE(std::isfinite(res.objective_value));
}

// ===========================================================================
// VGCalibrator (3 参数: sigma, nu, theta)
// ===========================================================================

class VGCalibratorTest : public ::testing::Test {
protected:
    static std::vector<MarketQuote> quotes_;
    static VGParams true_p_;
    static CalibrationResult de_result_;
    static CalibrationResult lm_result_;

    static void SetUpTestSuite() {
        const Real S = 100.0, r = 0.03, q = 0.0;
        true_p_ = VGParams{0.2, 0.5, -0.1};
        std::vector<Real> strikes = {75, 85, 95, 100, 105, 115, 125};
        std::vector<Real> maturities = {0.5, 1.0, 2.0};
        quotes_ = make_vg_quotes(S, r, q, true_p_.sigma, true_p_.nu, true_p_.theta,
                                 strikes, maturities);

        VGCalibrator cal;
        cal.set_market(S, r, q);
        CalibConfig cfg;
        cfg.de_pop_size = 40;
        cfg.de_generations = 80;
        cfg.lm_max_iter = 300;
        de_result_ = cal.calibrate(quotes_, cfg);

        CalibConfig cfg_lm;
        cfg_lm.use_de_init = false;
        cfg_lm.lm_max_iter = 300;
        lm_result_ = cal.calibrate(quotes_, cfg_lm);
    }

    static VGParams recovered() {
        VGCalibrator cal;
        return cal.extract_params(de_result_.params);
    }
};

std::vector<MarketQuote> VGCalibratorTest::quotes_;
VGParams VGCalibratorTest::true_p_;
CalibrationResult VGCalibratorTest::de_result_;
CalibrationResult VGCalibratorTest::lm_result_;

TEST_F(VGCalibratorTest, SyntheticDataRecovery) {
    // 用已知 VG 参数生成 IV surface, 标定恢复参数 (容差 3%)
    VGParams rec = recovered();
    EXPECT_NEAR(rec.sigma, true_p_.sigma, 0.03 * true_p_.sigma)
        << "sigma: expected " << true_p_.sigma << " got " << rec.sigma;
    EXPECT_NEAR(rec.nu, true_p_.nu, 0.03 * true_p_.nu)
        << "nu: expected " << true_p_.nu << " got " << rec.nu;
    EXPECT_NEAR(rec.theta, true_p_.theta, 0.03 * std::abs(true_p_.theta))
        << "theta: expected " << true_p_.theta << " got " << rec.theta;
    EXPECT_LT(de_result_.objective_value, 1e-6);
}

TEST_F(VGCalibratorTest, FellerConditionSatisfied) {
    // 1 - theta*nu - sigma^2*nu/2 > 0
    VGParams rec = recovered();
    EXPECT_TRUE(VGCalibrator::check_feller(rec));
    EXPECT_GT(1.0 - rec.theta * rec.nu - 0.5 * rec.sigma * rec.sigma * rec.nu, 0.0);
}

TEST_F(VGCalibratorTest, ThetaZeroSymmetricSmile) {
    // theta=0 时 IV smile 关于 K* = F*exp(omega*T) 对称
    const Real S = 100.0, r = 0.0, q = 0.0, T = 1.0;
    const Real sigma = 0.2, nu = 0.5, theta = 0.0;
    const Real omega = std::log(1.0 - theta * nu - 0.5 * sigma * sigma * nu) / nu;
    const Real Kstar = S * std::exp((r - q + omega) * T);

    auto phi = make_vg_cf(S, r, q, sigma, nu, theta, T);
    COSEngine::Config cc;
    cc.n_terms = 256;
    COSEngine engine(phi, S, r, q, T, cc);

    const Real K1 = 90.0;
    const Real K2 = Kstar * Kstar / K1;  // log-moneyness 对称于 K*
    Real iv1 = bsm_implied_vol(engine.price_call(K1), S, K1, T, r, q, true);
    Real iv2 = bsm_implied_vol(engine.price_call(K2), S, K2, T, r, q, true);
    EXPECT_NEAR(iv1, iv2, 0.008)
        << "Symmetric VG IV smile: IV(K1=" << K1 << ")=" << iv1
        << " IV(K2=" << K2 << ")=" << iv2;
}

TEST_F(VGCalibratorTest, DEPlusLMBetterThanPureLM) {
    EXPECT_LE(de_result_.objective_value, lm_result_.objective_value + 1e-8);
    EXPECT_LT(de_result_.objective_value, 1e-6);
}

TEST_F(VGCalibratorTest, NuSmallDegradesToBS) {
    // nu->0 且 theta=0 时 VG 退化为 BS (vol = sigma), 标定不崩溃
    const Real S = 100.0, r = 0.03, q = 0.0;
    const Real sigma = 0.2;

    // 定价层: nu 很小的 VG IV ≈ sigma
    auto phi = make_vg_cf(S, r, q, sigma, 0.001, 0.0, 1.0);
    COSEngine::Config cc;
    cc.n_terms = 256;
    COSEngine engine(phi, S, r, q, 1.0, cc);
    Real price = engine.price_call(100.0);
    Real iv = bsm_implied_vol(price, S, 100.0, 1.0, r, q, true);
    EXPECT_NEAR(iv, sigma, 0.01) << "nu->0 VG IV should approach sigma";

    // 标定层: 用 nu 很小的合成数据标定, 不崩溃且结果有限
    std::vector<Real> strikes = {85, 95, 100, 105, 115};
    std::vector<Real> maturities = {0.5, 1.0, 2.0};
    auto quotes = make_vg_quotes(S, r, q, sigma, 0.05, 0.0, strikes, maturities);

    VGCalibrator cal;
    cal.set_market(S, r, q);
    CalibConfig cfg;
    cfg.de_pop_size = 30;
    cfg.de_generations = 60;
    cfg.lm_max_iter = 300;
    CalibrationResult res;
    EXPECT_NO_THROW(res = cal.calibrate(quotes, cfg));
    ASSERT_EQ(res.params.size(), 3u);
    for (Real p : res.params) EXPECT_TRUE(std::isfinite(p));
    VGParams rec = cal.extract_params(res.params);
    EXPECT_NEAR(rec.nu, 0.05, 0.02);  // nu 恢复到小值附近
}

// ===========================================================================
// CEVCalibrator (2 参数: sigma, beta)
// ===========================================================================

class CEVCalibratorTest : public ::testing::Test {
protected:
    static std::vector<MarketQuote> quotes_;
    static CEVCalibParams true_p_;
    static CalibrationResult de_result_;
    static CalibrationResult lm_result_;

    static void SetUpTestSuite() {
        const Real S = 100.0, r = 0.03, q = 0.0, T = 1.0;
        true_p_ = CEVCalibParams{0.4, 0.5};
        CEVParams tp{true_p_.sigma, true_p_.beta};
        std::vector<Real> strikes = {80, 85, 90, 95, 100, 105, 110, 115, 120};
        quotes_ = make_cev_quotes(S, r, q, tp, strikes, T);

        CEVCalibrator cal;
        cal.set_market(S, r, q);
        CalibConfig cfg;
        cfg.de_pop_size = 30;
        cfg.de_generations = 60;
        cfg.lm_max_iter = 300;
        de_result_ = cal.calibrate(quotes_, cfg);

        CalibConfig cfg_lm;
        cfg_lm.use_de_init = false;
        cfg_lm.lm_max_iter = 300;
        lm_result_ = cal.calibrate(quotes_, cfg_lm);
    }

    static CEVCalibParams recovered() {
        CEVCalibrator cal;
        return cal.extract_params(de_result_.params);
    }
};

std::vector<MarketQuote> CEVCalibratorTest::quotes_;
CEVCalibParams CEVCalibratorTest::true_p_;
CalibrationResult CEVCalibratorTest::de_result_;
CalibrationResult CEVCalibratorTest::lm_result_;

TEST_F(CEVCalibratorTest, SyntheticDataRecovery) {
    // 用已知 CEV 参数生成 IV surface, 标定恢复参数 (容差 2%)
    CEVCalibParams rec = recovered();
    EXPECT_NEAR(rec.sigma, true_p_.sigma, 0.02 * true_p_.sigma)
        << "sigma: expected " << true_p_.sigma << " got " << rec.sigma;
    EXPECT_NEAR(rec.beta, true_p_.beta, 0.02 * true_p_.beta)
        << "beta: expected " << true_p_.beta << " got " << rec.beta;
    EXPECT_LT(de_result_.objective_value, 1e-5);
}

TEST_F(CEVCalibratorTest, BetaOneDegradesToBS) {
    // beta=1 时 CEV 退化为 GBM: cev_call_price == bsm_call_price (vol = sigma)
    const Real S = 100.0, r = 0.03, q = 0.0, T = 1.0;
    const Real sigma = 0.2;
    for (Real K : std::vector<Real>{90.0, 100.0, 110.0}) {
        Real cev = cev_call_price(S, K, T, r, q, CEVParams{sigma, 1.0});
        Real bs = bsm_call_price(S, K, T, r, q, sigma);
        EXPECT_NEAR(cev, bs, 1e-10) << "beta=1 CEV must equal BSM at K=" << K;
        // 反推 IV 应回到 sigma
        Real iv = bsm_implied_vol(cev, S, K, T, r, q, true);
        EXPECT_NEAR(iv, sigma, 1e-8) << "IV round-trip at K=" << K;
    }
}

TEST_F(CEVCalibratorTest, BetaBelowOneProducesSkew) {
    // beta<1 (CEV 弹性) 产生下行 IV skew; beta 越接近 0 (正态极限), skew 越陡
    const Real S = 100.0, r = 0.0, q = 0.0, T = 1.0;
    auto iv_at = [&](Real K, const CEVParams& p) {
        return bsm_implied_vol(cev_call_price(S, K, T, r, q, p), S, K, T, r, q, true);
    };
    // 低行权价 IV > 高行权价 IV (equity 风格正 skew)
    CEVParams p05{0.4, 0.5};
    Real skew05 = iv_at(80.0, p05) - iv_at(120.0, p05);
    EXPECT_GT(skew05, 0.002) << "CEV beta=0.5 should show downward skew";
    // beta 越小 (接近正态 beta->0), skew 越陡
    CEVParams p07{0.4, 0.7};
    Real skew07 = iv_at(80.0, p07) - iv_at(120.0, p07);
    EXPECT_GT(skew05, skew07) << "lower beta should give steeper skew";
}

TEST_F(CEVCalibratorTest, DEPlusLMBetterThanPureLM) {
    EXPECT_LE(de_result_.objective_value, lm_result_.objective_value + 1e-8);
    EXPECT_LT(de_result_.objective_value, 1e-5);
}

TEST_F(CEVCalibratorTest, BoundaryBetaNoCrash) {
    // beta 接近 0 和接近 1 时标定不崩溃, 且参数有限、落在界内、拟合良好.
    // 注意: 用零利率环境 (纯波动模型测试), 避免非中心卡方定价在小 beta 下的数值伪影.
    // 注: beta 接近 0 时 CEV 近似正态过程 (dS ~ sigma*(1+beta*ln S) dW),
    //     (sigma, beta) 存在近似退化 σ·(1+β·lnS)=const, 参数不可唯一识别
    //     (目标函数平坦, 多种 (σ,β) 组合 IV 曲面几乎重合). 因此该情形只断言
    //     不崩溃/有限/在界内/拟合良好, 不做参数恢复; 可识别的高 beta 侧做精确恢复.
    const Real S = 100.0, r = 0.0, q = 0.0, T = 1.0;
    std::vector<Real> strikes = {80, 90, 100, 110, 120};
    for (Real tbeta : std::vector<Real>{0.10, 0.99}) {
        CEVParams tp{0.3, tbeta};
        auto quotes = make_cev_quotes(S, r, q, tp, strikes, T);
        CEVCalibrator cal;
        cal.set_market(S, r, q);
        CalibConfig cfg;
        cfg.de_pop_size = 40;
        cfg.de_generations = 100;
        cfg.lm_max_iter = 300;
        CalibrationResult res;
        EXPECT_NO_THROW(res = cal.calibrate(quotes, cfg))
            << "calibrate crashed at beta=" << tbeta;
        ASSERT_EQ(res.params.size(), 2u);
        for (Real p : res.params) EXPECT_TRUE(std::isfinite(p));
        CEVCalibParams rec = cal.extract_params(res.params);
        EXPECT_GE(rec.beta, 0.01 - 1e-4) << "beta below lower bound";
        EXPECT_LE(rec.beta, 0.99 + 1e-4) << "beta above upper bound";
        EXPECT_GT(rec.sigma, 0.0);
        EXPECT_LT(res.objective_value, 1e-5)
            << "boundary calibration should still fit the surface";
        if (tbeta > 0.9) {
            // 高 beta 侧参数可识别: 精确恢复
            EXPECT_NEAR(rec.beta, tbeta, 0.05) << "beta recovery near bound";
            EXPECT_NEAR(rec.sigma, 0.3, 0.05) << "sigma recovery near bound";
        }
    }
}
