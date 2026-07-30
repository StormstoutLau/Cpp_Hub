// Phase 4 LITE - G2 整改项: Python 交叉验证测试
//
// 验证 Cpp_Hub C++ 实现与 Python (numpy/scipy/statsmodels) 基准值一致:
//   1. VaR 模块 (Historical/Parametric/Weighted) - 算法逻辑一致性
//   2. SSVI 公式 (Power-law 参数化) - 位精确匹配 (1e-12)
//   3. SSVI 无套利条件 - 与 Python 一致
//
// 基准来源: tests/validation/python/cross_validate_var.py
//           tests/validation/python/cross_validate_calibration.py
// 生成的 JSON: tests/validation/benchmarks_var.json
//              tests/validation/benchmarks_calib.json
//
// 容差:
//   - SSVI 公式 (直接计算): 1e-12 (位精确)
//   - VaR 算法 (同 PnL 输入): 1e-12 (算法相同)
//   - Parametric VaR (z_alpha): 1e-12 (scipy.stats.norm.ppf 精确)
//   - 无套利条件: 完全一致 (布尔值)
//
// 参考:
//   - Gatheral-Jacquier (2014) "Arbitrage-free SVI volatility surfaces" arXiv:1204.0646
//   - Cpp_Hub AUDIT_CHECKLIST.md G2 项

#include <gtest/gtest.h>
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/risk/var/historical_var.hpp"
#include "cpphub/risk/var/parametric_var.hpp"
#include "cpphub/models/vol_surface/ssvi.hpp"
#include <vector>
#include <cmath>
#include <numeric>
#include <random>

using namespace cpphub::v1;

namespace {

// ===========================================================================
// Python 基准值 (从 benchmarks_var.json 提取,seed=42, n_samples=5000)
// 生成命令: python cross_validate_var.py --seed 42 --n-samples 5000
// ===========================================================================

// PnL 统计量 (Python numpy 计算)
constexpr double kPnlMeanPython = 0.00027092810659017487;
constexpr double kPnlStdPython  = 0.023828037914391854;

// Parametric VaR 基准 (Python scipy.stats.norm.ppf)
constexpr double kZAlpha99 = -2.3263478740408408;   // Φ⁻¹(0.01)
constexpr double kZAlpha95 = -1.6448536269514729;   // Φ⁻¹(0.05)
constexpr double kParametricVar99 = 0.05516137723811986;
constexpr double kParametricVar95 = 0.03892270648003446;

// ===========================================================================
// Python 基准值 (从 benchmarks_calib.json 提取)
// ===========================================================================

// SSVI Power-law 公式验证 (直接计算,应位精确)
// 参数: rho=-0.3, eta=1.0, gamma=0.25, theta=0.16, k=-0.2
constexpr double kSsviRho   = -0.3;
constexpr double kSsviEta   = 1.0;
constexpr double kSsviGamma = 0.25;
constexpr double kSsviTheta = 0.16;
constexpr double kSsviK     = -0.2;
constexpr double kSsviPhiExpected     = 1.5811388300841898;   // η * θ^(-γ)
constexpr double kSsviTotalVarExpected = 0.17844272511244674;
constexpr double kSsviImpliedVolExpected = 0.42242481592875997;

// ===========================================================================
// 辅助函数: 生成确定性 PnL (C++ MT19937,与 Python 不同 RNG 但相同分布)
// ===========================================================================

std::vector<Real> generate_test_pnl(Size n = 1000, uint64_t seed = 12345) {
    // 用 C++ 自己的 RNG 生成 PnL,用于验证 VaR 算法逻辑
    // (不验证 RNG 一致性,只验证算法实现)
    std::mt19937_64 gen(seed);
    std::normal_distribution<Real> normal(0.001, 0.02);
    std::bernoulli_distribution tail(0.05);
    std::normal_distribution<Real> tail_dist(-0.005, 0.06);

    std::vector<Real> pnl;
    pnl.reserve(n);
    for (Size i = 0; i < n; ++i) {
        if (tail(gen)) {
            pnl.push_back(tail_dist(gen));
        } else {
            pnl.push_back(normal(gen));
        }
    }
    return pnl;
}

// ===========================================================================
// 独立实现 Python VaR 算法 (作为 reference,镜像 cross_validate_var.py)
// ===========================================================================

double python_historical_var_linear(const std::vector<Real>& pnl, double confidence) {
    // 镜像 Python historical_var_linear:
    // losses = sort(-pnl), index = q*(n-1), linear interp
    std::vector<Real> losses;
    losses.reserve(pnl.size());
    for (auto p : pnl) losses.push_back(-p);
    std::sort(losses.begin(), losses.end());
    Size n = losses.size();
    double q = 1.0 - confidence;
    double index = q * (n - 1);
    Size lo = static_cast<Size>(std::floor(index));
    Size hi = static_cast<Size>(std::ceil(index));
    if (lo == hi) return -losses[lo];
    double frac = index - lo;
    double val = losses[lo] + frac * (losses[hi] - losses[lo]);
    return -val;
}

double python_historical_var_conservative(const std::vector<Real>& pnl, double confidence) {
    // 镜像 Python historical_var_conservative
    std::vector<Real> losses;
    losses.reserve(pnl.size());
    for (auto p : pnl) losses.push_back(-p);
    std::sort(losses.begin(), losses.end());
    Size n = losses.size();
    double q = 1.0 - confidence;
    double index = q * (n - 1);
    Size lo = static_cast<Size>(std::floor(index));
    Size hi = static_cast<Size>(std::ceil(index));
    if (lo >= n) lo = n - 1;
    if (hi >= n) hi = n - 1;
    if (lo == hi) return -losses[lo];
    if (std::abs(losses[lo]) >= std::abs(losses[hi])) return -losses[lo];
    return -losses[hi];
}

double python_historical_var_empirical(const std::vector<Real>& pnl, double confidence) {
    // 镜像 Python historical_var_empirical
    std::vector<Real> losses;
    losses.reserve(pnl.size());
    for (auto p : pnl) losses.push_back(-p);
    std::sort(losses.begin(), losses.end());
    Size n = losses.size();
    double q = 1.0 - confidence;
    Size idx = static_cast<Size>(q * n);
    if (idx >= n) idx = n - 1;
    return -losses[idx];
}

double python_parametric_var_normal(double mean, double std, double confidence) {
    // VaR = -(mean + std * z_alpha), z_alpha = Φ⁻¹(1-confidence)
    double z = inv_normal_cdf(1.0 - confidence);
    return -(mean + std * z);
}

}  // namespace

// ===========================================================================
// 1. Parametric VaR 交叉验证 (vs Python scipy.stats.norm.ppf)
// ===========================================================================

TEST(PythonCrossValidation, ParametricVarNormal99) {
    // Python: parametric_var_normal(mean, std, 0.99) = 0.05516137723811986
    double var = python_parametric_var_normal(kPnlMeanPython, kPnlStdPython, 0.99);
    EXPECT_NEAR(var, kParametricVar99, 1e-12)
        << "Parametric VaR 99% 应与 Python scipy 一致";

    // C++ ParametricVaR 实现
    PortfolioStats stats;
    stats.mean = kPnlMeanPython;
    stats.variance = kPnlStdPython * kPnlStdPython;
    ParametricVaR pv(stats, 0.99, 1);
    double cpp_var = pv.var(ParametricMethod::Normal);
    EXPECT_NEAR(cpp_var, kParametricVar99, 1e-12)
        << "C++ ParametricVaR 应与 Python scipy 一致";
}

TEST(PythonCrossValidation, ParametricVarNormal95) {
    // Python: parametric_var_normal(mean, std, 0.95) = 0.03892270648003446
    double var = python_parametric_var_normal(kPnlMeanPython, kPnlStdPython, 0.95);
    EXPECT_NEAR(var, kParametricVar95, 1e-12);

    PortfolioStats stats;
    stats.mean = kPnlMeanPython;
    stats.variance = kPnlStdPython * kPnlStdPython;
    ParametricVaR pv(stats, 0.95, 1);
    double cpp_var = pv.var(ParametricMethod::Normal);
    EXPECT_NEAR(cpp_var, kParametricVar95, 1e-12);
}

TEST(PythonCrossValidation, ZAlphaQuantile) {
    // 验证 C++ inv_normal_cdf 与 scipy.stats.norm.ppf 一致
    EXPECT_NEAR(inv_normal_cdf(0.01), kZAlpha99, 1e-12);
    EXPECT_NEAR(inv_normal_cdf(0.05), kZAlpha95, 1e-12);
}

// ===========================================================================
// 2. Historical VaR 算法交叉验证 (C++ HistoricalVaR vs Python reference)
// ===========================================================================

TEST(PythonCrossValidation, HistoricalVarLinearConsistency) {
    // 用 C++ 生成 PnL,验证 C++ HistoricalVaR 与 Python 算法实现一致
    auto pnl = generate_test_pnl(1000, 12345);

    // C++ HistoricalVaR
    HistoricalVaR hv(pnl, 0.99, 1);
    double cpp_var = hv.var(QuantileInterpolation::Linear);

    // Python reference (独立实现)
    double py_var = python_historical_var_linear(pnl, 0.99);

    EXPECT_NEAR(cpp_var, py_var, 1e-12)
        << "C++ HistoricalVaR (Linear) 应与 Python 算法一致";
}

TEST(PythonCrossValidation, HistoricalVarConservativeConsistency) {
    auto pnl = generate_test_pnl(1000, 12345);

    HistoricalVaR hv(pnl, 0.99, 1);
    double cpp_var = hv.var(QuantileInterpolation::Conservative);
    double py_var = python_historical_var_conservative(pnl, 0.99);

    EXPECT_NEAR(cpp_var, py_var, 1e-12)
        << "C++ HistoricalVaR (Conservative) 应与 Python 算法一致";
}

TEST(PythonCrossValidation, HistoricalVarEmpiricalConsistency) {
    auto pnl = generate_test_pnl(1000, 12345);

    HistoricalVaR hv(pnl, 0.99, 1);
    double cpp_var = hv.var(QuantileInterpolation::Empirical);
    double py_var = python_historical_var_empirical(pnl, 0.99);

    EXPECT_NEAR(cpp_var, py_var, 1e-12)
        << "C++ HistoricalVaR (Empirical) 应与 Python 算法一致";
}

TEST(PythonCrossValidation, HistoricalVarMultipleConfidenceLevels) {
    auto pnl = generate_test_pnl(2000, 54321);

    HistoricalVaR hv(pnl, 0.95, 1);
    EXPECT_NEAR(hv.var(QuantileInterpolation::Linear),
                python_historical_var_linear(pnl, 0.95), 1e-12);
    EXPECT_NEAR(hv.var(QuantileInterpolation::Conservative),
                python_historical_var_conservative(pnl, 0.95), 1e-12);
    EXPECT_NEAR(hv.var(QuantileInterpolation::Empirical),
                python_historical_var_empirical(pnl, 0.95), 1e-12);
}

// ===========================================================================
// 3. SSVI Power-law 公式交叉验证 (vs Python scipy, 位精确)
// ===========================================================================

TEST(PythonCrossValidation, SsviPowerLawPhi) {
    // Python: ssvi_power_law_phi(0.16, 1.0, 0.25) = 1.5811388300841898
    // 公式: φ(θ) = η * θ^(-γ) = 1.0 * 0.16^(-0.25)
    // 手算: 0.16^0.25 = sqrt(sqrt(0.16)) = sqrt(0.4) = 0.63246, 1/0.63246 = 1.58114
    auto params = SSVI::Power_law(kSsviRho, kSsviEta, kSsviGamma, {kSsviTheta});
    SSVI ssvi(params);

    double phi = ssvi.params().phi(kSsviTheta);
    EXPECT_NEAR(phi, kSsviPhiExpected, 1e-12)
        << "SSVI Power-law φ(θ) 应与 Python 一致";
}

TEST(PythonCrossValidation, SsviTotalVariance) {
    // Python: ssvi_total_variance(-0.2, 0.16, 1.5811, -0.3) = 0.17844272511244674
    auto params = SSVI::Power_law(kSsviRho, kSsviEta, kSsviGamma, {kSsviTheta});
    SSVI ssvi(params);

    double w = ssvi.total_variance(kSsviK, kSsviTheta);
    EXPECT_NEAR(w, kSsviTotalVarExpected, 1e-12)
        << "SSVI 总方差应与 Python 位精确一致";
}

TEST(PythonCrossValidation, SsviImpliedVol) {
    // Python: sqrt(0.17844/1.0) = 0.42242481592875997
    auto params = SSVI::Power_law(kSsviRho, kSsviEta, kSsviGamma, {kSsviTheta});
    SSVI ssvi(params);

    double iv = ssvi.implied_vol(kSsviK, 1.0, kSsviTheta);
    EXPECT_NEAR(iv, kSsviImpliedVolExpected, 1e-12)
        << "SSVI 隐含波动率应与 Python 位精确一致";
}

// ===========================================================================
// 4. SSVI 无套利条件交叉验证 (vs Python)
// ===========================================================================

TEST(PythonCrossValidation, SsviNoArbitragePowerLaw) {
    // Python: no_arbitrage_check 返回所有条件 true
    // 参数: rho=-0.3, eta=1.0, gamma=0.25, theta_slice=[0.04, 0.09, 0.16, 0.25]
    std::vector<Real> thetas = {0.04, 0.09, 0.16, 0.25};
    auto params = SSVI::Power_law(kSsviRho, kSsviEta, kSsviGamma, thetas);
    SSVI ssvi(params);

    EXPECT_TRUE(ssvi.check_no_arbitrage())
        << "SSVI Power-law (rho=-0.3, eta=1.0, gamma=0.25) 应无套利";

    // Python 验证的具体条件
    EXPECT_TRUE(std::abs(kSsviRho) < 1.0) << "|rho| < 1";
    Real phi = params.phi(kSsviTheta);
    EXPECT_GT(phi, 0.0) << "phi > 0";
    EXPECT_LT(phi * kSsviTheta * (1.0 + std::abs(kSsviRho)), 4.0)
        << "butterfly 充分条件: phi*theta*(1+|rho|) < 4";
    EXPECT_TRUE(kSsviGamma > 0.0 && kSsviGamma < 0.5)
        << "gamma 在 (0, 0.5) 内";
}

TEST(PythonCrossValidation, SsviCalendarArbitrage) {
    // Power-law: θ*φ(θ) = η*θ^(1-γ), d/dθ = η*(1-γ)*θ^(-γ) > 0 (γ<1)
    // 应满足 calendar 无套利
    std::vector<Real> thetas = {0.04, 0.09, 0.16, 0.25};
    auto params = SSVI::Power_law(kSsviRho, kSsviEta, kSsviGamma, thetas);
    SSVI ssvi(params);

    EXPECT_TRUE(ssvi.check_calendar_arbitrage())
        << "Power-law (gamma<0.5) 应满足 calendar 无套利";

    // 验证 θ*φ(θ) 单调递增 (Python 也验证此条件)
    for (Size i = 1; i < thetas.size(); ++i) {
        Real prod_prev = thetas[i-1] * params.phi(thetas[i-1]);
        Real prod_curr = thetas[i] * params.phi(thetas[i]);
        EXPECT_GT(prod_curr, prod_prev)
            << "θ*φ(θ) 应单调递增 (calendar 无套利)";
    }
}

TEST(PythonCrossValidation, SsviButterflyArbitrage) {
    // 验证充分 butterfly 无套利条件
    std::vector<Real> thetas = {0.04, 0.09, 0.16, 0.25};
    auto params = SSVI::Power_law(kSsviRho, kSsviEta, kSsviGamma, thetas);
    SSVI ssvi(params);

    EXPECT_TRUE(ssvi.check_butterfly_arbitrage())
        << "应满足 butterfly 充分无套利条件";

    // Python 验证: phi*theta*(1+|rho|) < 4 对所有 theta
    Real abs_rho = std::abs(kSsviRho);
    for (Real theta : thetas) {
        Real phi = params.phi(theta);
        EXPECT_LT(phi * theta * (1.0 + abs_rho), 4.0)
            << "theta=" << theta << " 应满足 butterfly 条件";
    }
}

// ===========================================================================
// 5. SSVI 参数化工厂验证 (Heston-like vs Power-law)
// ===========================================================================

TEST(PythonCrossValidation, SsviHestonLikePhi) {
    // Heston-like: φ(θ) = η * θ^(-λ)
    // 当 λ=γ 时,与 Power-law 相同
    Real theta = 0.16;
    Real eta = 1.0;
    Real lambda = 0.25;  // 等于 Power-law 的 gamma

    auto params_h = SSVI::Heston_like(kSsviRho, eta, lambda, {theta});
    auto params_p = SSVI::Power_law(kSsviRho, eta, lambda, {theta});

    EXPECT_NEAR(params_h.phi(theta), params_p.phi(theta), 1e-15)
        << "Heston-like (λ=γ) 应与 Power-law 相同";
}

TEST(PythonCrossValidation, SsviTotalVarianceAtATM) {
    // ATM (k=0): w(0, θ) = θ/2 * (1 + 0 + sqrt(ρ² + 1 - ρ²)) = θ/2 * (1+1) = θ
    auto params = SSVI::Power_law(kSsviRho, kSsviEta, kSsviGamma, {kSsviTheta});
    SSVI ssvi(params);

    double w_atm = ssvi.total_variance(0.0, kSsviTheta);
    EXPECT_NEAR(w_atm, kSsviTheta, 1e-15)
        << "ATM 总方差应等于 θ";
}

// ===========================================================================
// 6. VaR 回测交叉验证 (Kupiec POF 检验,验证统计量计算)
// ===========================================================================

TEST(PythonCrossValidation, VarBacktestLogic) {
    // 验证 Kupiec POF 统计量计算 (与 Python/标准实现一致)
    // Case 1: 100 天,3 天突破,99% 置信度 (期望 1 次)
    //   LR_POF ≈ 2.63 < 6.635 (不拒绝,突破偏多但不显著)
    {
        Size n = 100;
        Size violations = 3;
        double p_expected = 0.01;

        double num = std::pow(1 - p_expected, n - violations) * std::pow(p_expected, violations);
        double den = std::pow(1 - static_cast<double>(violations) / n, n - violations)
                     * std::pow(static_cast<double>(violations) / n, violations);
        double lr_pof = -2.0 * std::log(num / den);

        // 验证统计量数值 (Python: scipy.stats.chi2.logsf 一致)
        EXPECT_NEAR(lr_pof, 2.6324, 1e-3)
            << "Kupiec POF (3/100, 99%) 应约等于 2.63";
        EXPECT_LT(lr_pof, 6.635)
            << "3/100 突破在 99% 不应拒绝";
    }

    // Case 2: 250 天,10 天突破,99% 置信度 (期望 2.5 次)
    //   LR_POF ≈ 9.97 > 6.635 (拒绝,突破过多)
    {
        Size n = 250;
        Size violations = 10;
        double p_expected = 0.01;

        double num = std::pow(1 - p_expected, n - violations) * std::pow(p_expected, violations);
        double den = std::pow(1 - static_cast<double>(violations) / n, n - violations)
                     * std::pow(static_cast<double>(violations) / n, violations);
        double lr_pof = -2.0 * std::log(num / den);

        EXPECT_GT(lr_pof, 6.635)
            << "10/250 突破在 99% 应拒绝模型";
    }
}
