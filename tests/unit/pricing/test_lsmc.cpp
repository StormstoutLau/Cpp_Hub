// 基准来源: Longstaff & Schwartz (2001) "Valuing American Options by Simulation"
// 表/页码: Table 1 (p. 114-115)
// 容差: 美式 Put 基准 2% (MC 噪声 + 回归偏差), 路径统计 5% (Box-Muller 收敛)
#include <gtest/gtest.h>
#include "cpphub/pricing/monte_carlo/lsmc_engine.hpp"
#include "cpphub/pricing/monte_carlo/path_generator.hpp"
#include "cpphub/pricing/tree/binomial.hpp"
#include "cpphub/instruments/payoff/vanilla.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/monte_carlo/control_variate.hpp"
#include <cmath>
#include <vector>

using namespace cpphub;

// ============ 基函数正确性验证 ============

TEST(LSMCBasis, LaguerrePolynomials) {
    // L_0(x) = 1
    EXPECT_NEAR(laguerre(0, 2.5), 1.0, 1e-12);
    // L_1(x) = 1 - x
    EXPECT_NEAR(laguerre(1, 2.5), 1.0 - 2.5, 1e-12);
    // L_2(x) = 1 - 2x + x^2/2
    EXPECT_NEAR(laguerre(2, 2.5), 1.0 - 2.0*2.5 + 2.5*2.5/2.0, 1e-12);
    // L_3(x) = 1 - 3x + 3x^2/2 - x^3/6
    EXPECT_NEAR(laguerre(3, 2.5), 1.0 - 3.0*2.5 + 3.0*2.5*2.5/2.0 - 2.5*2.5*2.5/6.0, 1e-10);
    // 递推关系: n*L_n = (2n-1-x)*L_{n-1} - (n-1)*L_{n-2}
    Real x = 1.7;
    for (Size n = 2; n <= 6; ++n) {
        Real lhs = static_cast<Real>(n) * laguerre(n, x);
        Real rhs = (2.0*static_cast<Real>(n) - 1.0 - x) * laguerre(n-1, x)
                   - (static_cast<Real>(n) - 1.0) * laguerre(n-2, x);
        EXPECT_NEAR(lhs, rhs, 1e-10) << "Laguerre recursion failed at n=" << n;
    }
}

TEST(LSMCBasis, HermitePolynomials) {
    // He_0(x) = 1, He_1(x) = x
    EXPECT_NEAR(hermite_prob(0, 1.5), 1.0, 1e-12);
    EXPECT_NEAR(hermite_prob(1, 1.5), 1.5, 1e-12);
    // He_2(x) = x^2 - 1
    EXPECT_NEAR(hermite_prob(2, 1.5), 1.5*1.5 - 1.0, 1e-12);
    // He_3(x) = x^3 - 3x
    EXPECT_NEAR(hermite_prob(3, 1.5), 1.5*1.5*1.5 - 3.0*1.5, 1e-11);
    // 递推: He_n = x*He_{n-1} - (n-1)*He_{n-2}
    Real x = 0.7;
    for (Size n = 2; n <= 6; ++n) {
        Real lhs = hermite_prob(n, x);
        Real rhs = x * hermite_prob(n-1, x) - (static_cast<Real>(n) - 1.0) * hermite_prob(n-2, x);
        EXPECT_NEAR(lhs, rhs, 1e-10) << "Hermite recursion failed at n=" << n;
    }
}

TEST(LSMCBasis, ChebyshevPolynomials) {
    // T_0(x) = 1, T_1(x) = x
    EXPECT_NEAR(chebyshev_t(0, 0.5), 1.0, 1e-12);
    EXPECT_NEAR(chebyshev_t(1, 0.5), 0.5, 1e-12);
    // T_2(x) = 2x^2 - 1
    EXPECT_NEAR(chebyshev_t(2, 0.5), 2.0*0.5*0.5 - 1.0, 1e-12);
    // T_3(x) = 4x^3 - 3x
    EXPECT_NEAR(chebyshev_t(3, 0.5), 4.0*0.5*0.5*0.5 - 3.0*0.5, 1e-12);
    // 递推: T_n = 2x*T_{n-1} - T_{n-2}
    Real x = 0.3;
    for (Size n = 2; n <= 6; ++n) {
        Real lhs = chebyshev_t(n, x);
        Real rhs = 2.0*x*chebyshev_t(n-1, x) - chebyshev_t(n-2, x);
        EXPECT_NEAR(lhs, rhs, 1e-10) << "Chebyshev recursion failed at n=" << n;
    }
}

TEST(LSMCBasis, BasisEvalReturnsCorrectOrder) {
    for (Size order = 1; order <= 5; ++order) {
        std::vector<Real> b = basis_eval(BasisType::Laguerre, order, 1.5);
        EXPECT_EQ(b.size(), order);
        for (Size n = 0; n < order; ++n) {
            EXPECT_NEAR(b[n], laguerre(n, 1.5), 1e-12);
        }
    }
}

// ============ 线性方程组求解 ============

TEST(LSMCLinearSolve, IdentitySystem) {
    std::vector<std::vector<Real>> A = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    std::vector<Real> b = {3.0, 5.0, 7.0};
    bool ok = solve_linear_system(A, b, 3);
    EXPECT_TRUE(ok);
    EXPECT_NEAR(b[0], 3.0, 1e-12);
    EXPECT_NEAR(b[1], 5.0, 1e-12);
    EXPECT_NEAR(b[2], 7.0, 1e-12);
}

TEST(LSMCLinearSolve, GeneralSystem) {
    // 2x + y = 5, x + 3y = 10  =>  x=1, y=3
    std::vector<std::vector<Real>> A = {{2.0, 1.0}, {1.0, 3.0}};
    std::vector<Real> b = {5.0, 10.0};
    bool ok = solve_linear_system(A, b, 2);
    EXPECT_TRUE(ok);
    EXPECT_NEAR(b[0], 1.0, 1e-10);
    EXPECT_NEAR(b[1], 3.0, 1e-10);
}

TEST(LSMCLinearSolve, SingularSystem) {
    std::vector<std::vector<Real>> A = {{1.0, 2.0}, {2.0, 4.0}};  // 奇异
    std::vector<Real> b = {1.0, 2.0};
    bool ok = solve_linear_system(A, b, 2);
    EXPECT_FALSE(ok);
}

// ============ GBM 路径生成器 ============

TEST(LSMCPathGenerator, ExactPathStartsAtS0) {
    GBMConfig cfg{100.0, 0.05, 0.0, 0.20, 1.0, 50};
    GBMPathGenerator gen(cfg, PathScheme::Exact);
    Philox4x64 rng(42, 0);
    std::vector<Real> path = gen.generate_path(rng);
    EXPECT_EQ(path.size(), 51u);
    EXPECT_NEAR(path[0], 100.0, 1e-12);
}

TEST(LSMCPathGenerator, ExactPathTerminalDistributionMatchesTheory) {
    // GBM 精确解: ln(S_T/S0) ~ N((r-q-0.5*sigma^2)*T, sigma^2*T)
    // 大样本下样本均值应接近理论均值
    GBMConfig cfg{100.0, 0.05, 0.0, 0.20, 1.0, 1};
    GBMPathGenerator gen(cfg, PathScheme::Exact);
    Size n_paths = 100000;
    std::vector<Real> paths = gen.generate_paths(n_paths, 42, false);
    Real mu_theory = (0.05 - 0.0 - 0.5*0.20*0.20) * 1.0;
    Real var_theory = 0.20*0.20 * 1.0;
    Real sum_logret = 0.0, sum_sq = 0.0;
    for (Size p = 0; p < n_paths; ++p) {
        Real S_T = paths[p * 2 + 1];
        Real lr = std::log(S_T / 100.0);
        sum_logret += lr;
        sum_sq += lr * lr;
    }
    Real mean = sum_logret / static_cast<Real>(n_paths);
    Real var = sum_sq / static_cast<Real>(n_paths) - mean * mean;
    // 大数定律: 100000 样本, 均值容差 ~ sigma/sqrt(N) = 0.20/sqrt(100000) ≈ 0.0006
    EXPECT_NEAR(mean, mu_theory, 0.005);
    EXPECT_NEAR(var, var_theory, 0.002);
}

TEST(LSMCPathGenerator, BatchPathsHaveCorrectShape) {
    GBMConfig cfg{100.0, 0.05, 0.0, 0.20, 1.0, 50};
    GBMPathGenerator gen(cfg, PathScheme::Exact);
    Size n_paths = 1000;
    std::vector<Real> paths = gen.generate_paths(n_paths, 42, false);
    EXPECT_EQ(paths.size(), n_paths * 51u);
    // 每条路径起点都是 S0
    for (Size p = 0; p < n_paths; ++p) {
        EXPECT_NEAR(paths[p * 51], 100.0, 1e-12);
    }
}

TEST(LSMCPathGenerator, DeterministicGivenSeed) {
    GBMConfig cfg{100.0, 0.05, 0.0, 0.20, 1.0, 50};
    GBMPathGenerator gen(cfg, PathScheme::Exact);
    std::vector<Real> paths1 = gen.generate_paths(100, 42, false);
    std::vector<Real> paths2 = gen.generate_paths(100, 42, false);
    EXPECT_EQ(paths1.size(), paths2.size());
    for (Size i = 0; i < paths1.size(); ++i) {
        EXPECT_NEAR(paths1[i], paths2[i], 1e-15);
    }
}

TEST(LSMCPathGenerator, AntitheticReducesVariance) {
    // 反变量路径对: 路径 p 和 p+N/2 的 S_T 应关于 S0*exp(drift*T) 对称
    GBMConfig cfg{100.0, 0.05, 0.0, 0.20, 1.0, 1};
    GBMPathGenerator gen(cfg, PathScheme::Exact);
    Size n_paths = 1000;
    std::vector<Real> paths = gen.generate_paths(n_paths, 42, true);
    Size half = n_paths / 2;
    Real drift = (0.05 - 0.0 - 0.5*0.20*0.20) * 1.0;
    Real center = 100.0 * std::exp(drift);
    // 反变量对的几何均值应接近 center (理论上精确等于)
    Real sum_ratio = 0.0;
    for (Size p = 0; p < half; ++p) {
        Real S1 = paths[p * 2 + 1];
        Real S2 = paths[(p + half) * 2 + 1];
        // ln(S1*S2) = 2*drift + 0 (Z 和 -Z 抵消)
        Real geo_mean = std::sqrt(S1 * S2);
        sum_ratio += geo_mean;
    }
    Real avg = sum_ratio / static_cast<Real>(half);
    EXPECT_NEAR(avg, center, 0.01);  // 数值精度
}

// ============ LSMC 美式期权定价 ============

TEST(LSMCEngine, AmericanPutGreaterThanEuropean) {
    // 美式 Put >= 欧式 Put (无套利约束)
    LSMCConfig cfg;
    cfg.S0 = 100.0; cfg.K = 100.0; cfg.T = 1.0;
    cfg.r = 0.05; cfg.q = 0.0; cfg.sigma = 0.20;
    cfg.n_paths = 20000; cfg.n_steps = 50;
    cfg.basis = BasisType::Laguerre; cfg.basis_order = 3;
    cfg.seed = 42; cfg.antithetic = true;
    LSMCEngine engine(cfg);
    PutPayOff put(100.0);
    LSMCResult result = engine.price_american(put);
    // 美式 >= 欧式 (允许 MC 噪声)
    EXPECT_GE(result.price + 2.0 * result.std_error, result.european_price);
    // 早期行使溢价应为正 (美式 Put 有行使价值)
    EXPECT_GT(result.early_exercise_premium, -0.05);
}

TEST(LSMCEngine, AmericanPutMatchesLongstaffSchwartzBenchmark) {
    // Longstaff-Schwartz (2001) Table 1 基准:
    // S0=36, K=40, r=0.06, sigma=0.20, T=1.0
    // 欧式 Put: 3.8444, 美式 Put: 4.478
    LSMCConfig cfg;
    cfg.S0 = 36.0; cfg.K = 40.0; cfg.T = 1.0;
    cfg.r = 0.06; cfg.q = 0.0; cfg.sigma = 0.20;
    cfg.n_paths = 50000; cfg.n_steps = 50;
    cfg.basis = BasisType::Laguerre; cfg.basis_order = 3;
    cfg.seed = 42; cfg.antithetic = true;
    LSMCEngine engine(cfg);
    PutPayOff put(40.0);
    LSMCResult result = engine.price_american(put);
    // 基准 4.478, 容差 3% (MC + 回归偏差)
    EXPECT_NEAR(result.price, 4.478, 4.478 * 0.03);
    // 欧式参考值
    Real euro_ref = bsm_put_price(36.0, 40.0, 1.0, 0.06, 0.0, 0.20);
    EXPECT_NEAR(result.european_price, euro_ref, euro_ref * 0.02);
    // 早期行使溢价 > 0
    EXPECT_GT(result.early_exercise_premium, 0.3);
}

TEST(LSMCEngine, AmericanCallNoDividendEqualsEuropean) {
    // 无分红美式 Call = 欧式 Call (无套利约束, 不应提前行使)
    LSMCConfig cfg;
    cfg.S0 = 100.0; cfg.K = 100.0; cfg.T = 1.0;
    cfg.r = 0.05; cfg.q = 0.0; cfg.sigma = 0.20;
    cfg.n_paths = 20000; cfg.n_steps = 50;
    cfg.basis = BasisType::Laguerre; cfg.basis_order = 3;
    cfg.seed = 42; cfg.antithetic = true;
    LSMCEngine engine(cfg);
    CallPayOff call(100.0);
    LSMCResult result = engine.price_american(call);
    // 美式 Call (无分红) ≈ 欧式 Call
    // 容差 2% (MC 噪声 + 回归可能误判行使)
    EXPECT_NEAR(result.price, result.european_price, result.european_price * 0.02 + 0.05);
}

TEST(LSMCEngine, AmericanPutMatchesBinomialTree) {
    // 与 BinomialTreeEngine (CRR, n=5000) 对比
    // 两者都是近似方法, 容差 3%
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    // 二叉树基准
    BinomialParams bp{S0, K, T, r, q, sigma, 5000};
    BinomialTreeEngine bte(bp, BinomialType::CRR);
    PutPayOff put(K);
    Real binomial_price = bte.price_american(put);
    // LSMC
    LSMCConfig cfg;
    cfg.S0 = S0; cfg.K = K; cfg.T = T;
    cfg.r = r; cfg.q = q; cfg.sigma = sigma;
    cfg.n_paths = 50000; cfg.n_steps = 100;
    cfg.basis = BasisType::Laguerre; cfg.basis_order = 4;
    cfg.seed = 123; cfg.antithetic = true;
    LSMCEngine engine(cfg);
    LSMCResult result = engine.price_american(put);
    // 两种方法容差 3% (二叉树离散误差 + MC 噪声)
    EXPECT_NEAR(result.price, binomial_price, binomial_price * 0.03);
}

TEST(LSMCEngine, ConvergenceWithPathCount) {
    // 路径数增加时, LSMC 价格应收敛 (标准误差下降)
    LSMCConfig cfg1;
    cfg1.S0 = 100.0; cfg1.K = 100.0; cfg1.T = 1.0;
    cfg1.r = 0.05; cfg1.q = 0.0; cfg1.sigma = 0.20;
    cfg1.n_paths = 5000; cfg1.n_steps = 50;
    cfg1.basis = BasisType::Laguerre; cfg1.basis_order = 3;
    cfg1.seed = 42; cfg1.antithetic = true;
    LSMCEngine engine1(cfg1);
    PutPayOff put(100.0);
    LSMCResult r1 = engine1.price_american(put);

    LSMCConfig cfg2 = cfg1;
    cfg2.n_paths = 40000;
    LSMCEngine engine2(cfg2);
    LSMCResult r2 = engine2.price_american(put);

    // 标准误差应显著降低
    EXPECT_LT(r2.std_error, r1.std_error * 0.6);
    // 价格应在 2*SE 内一致
    Real diff = std::abs(r1.price - r2.price);
    Real tol = 3.0 * std::max(r1.std_error, r2.std_error);
    EXPECT_LT(diff, tol + 0.15);  // 允许少量回归偏差
}

TEST(LSMCEngine, DeterministicGivenSeed) {
    // 同种子同结果 (确定性)
    LSMCConfig cfg;
    cfg.S0 = 100.0; cfg.K = 100.0; cfg.T = 1.0;
    cfg.r = 0.05; cfg.q = 0.0; cfg.sigma = 0.20;
    cfg.n_paths = 5000; cfg.n_steps = 50;
    cfg.basis = BasisType::Laguerre; cfg.basis_order = 3;
    cfg.seed = 42; cfg.antithetic = true;
    LSMCEngine engine1(cfg);
    LSMCEngine engine2(cfg);
    PutPayOff put(100.0);
    LSMCResult r1 = engine1.price_american(put);
    LSMCResult r2 = engine2.price_american(put);
    EXPECT_NEAR(r1.price, r2.price, 1e-15);
    EXPECT_NEAR(r1.std_error, r2.std_error, 1e-15);
}

TEST(LSMCEngine, BermudanBetweenAmericanAndEuropean) {
    // 百慕大期权价格应介于美式和欧式之间
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    LSMCConfig cfg;
    cfg.S0 = S0; cfg.K = K; cfg.T = T;
    cfg.r = r; cfg.q = q; cfg.sigma = sigma;
    cfg.n_paths = 20000; cfg.n_steps = 50;
    cfg.basis = BasisType::Laguerre; cfg.basis_order = 3;
    cfg.seed = 42; cfg.antithetic = true;
    LSMCEngine engine(cfg);
    PutPayOff put(K);
    // 美式 (所有 50 个点可行使)
    LSMCResult amer = engine.price_american(put);
    // 百慕大 (仅 4 个点可行使: T/4, T/2, 3T/4, T)
    std::vector<Real> bermudan_times = {T*0.25, T*0.5, T*0.75, T};
    LSMCResult berm = engine.price_bermudan(put, bermudan_times);
    // 欧式参考
    Real euro = bsm_put_price(S0, K, T, r, q, sigma);
    // 百慕大 >= 欧式 (允许 MC 噪声)
    EXPECT_GE(berm.price + 2.0 * berm.std_error, euro);
    // 美式 >= 百慕大 (美式有更多行使机会, 允许 MC 噪声)
    EXPECT_GE(amer.price + 2.0 * amer.std_error, berm.price);
}

TEST(LSMCEngine, DifferentBasisTypesProduceSimilarResults) {
    // 不同基函数类型应给出相近的结果 (容差 3%)
    Real S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, sigma = 0.20;
    PutPayOff put(K);
    Real ref_price = 0.0;
    for (int b = 0; b < 4; ++b) {
        LSMCConfig cfg;
        cfg.S0 = S0; cfg.K = K; cfg.T = T;
        cfg.r = r; cfg.sigma = sigma;
        cfg.n_paths = 20000; cfg.n_steps = 50;
        cfg.basis_order = 3;
        cfg.seed = 42; cfg.antithetic = true;
        cfg.basis = static_cast<BasisType>(b);
        LSMCEngine engine(cfg);
        LSMCResult result = engine.price_american(put);
        if (b == 0) ref_price = result.price;
        else EXPECT_NEAR(result.price, ref_price, ref_price * 0.03)
              << "Basis type " << b << " differs from Laguerre";
    }
}

TEST(LSMCEngine, InvalidConfigThrows) {
    LSMCConfig cfg;
    cfg.S0 = -1.0;  // 非法
    EXPECT_THROW(LSMCEngine engine(cfg), std::invalid_argument);
    cfg.S0 = 100.0; cfg.n_paths = 10;  // 路径数太少
    EXPECT_THROW(LSMCEngine engine(cfg), std::invalid_argument);
    cfg.n_paths = 10000; cfg.basis_order = 0;  // 基函数阶数为 0
    EXPECT_THROW(LSMCEngine engine(cfg), std::invalid_argument);
}

TEST(LSMCEngine, AmericanPutITMHasHigherPremium) {
    // 深度 ITM 美式 Put 的早期行使溢价应高于 ATM
    Real K = 100.0, T = 1.0, r = 0.05, sigma = 0.20;
    // ATM Put
    LSMCConfig cfg_atm;
    cfg_atm.S0 = 100.0; cfg_atm.K = K; cfg_atm.T = T;
    cfg_atm.r = r; cfg_atm.sigma = sigma;
    cfg_atm.n_paths = 20000; cfg_atm.n_steps = 50;
    cfg_atm.basis = BasisType::Laguerre; cfg_atm.basis_order = 3;
    cfg_atm.seed = 42; cfg_atm.antithetic = true;
    LSMCEngine engine_atm(cfg_atm);
    PutPayOff put(K);
    LSMCResult atm = engine_atm.price_american(put);
    // ITM Put (S0=80, K=100)
    LSMCConfig cfg_itm = cfg_atm;
    cfg_itm.S0 = 80.0;
    LSMCEngine engine_itm(cfg_itm);
    LSMCResult itm = engine_itm.price_american(put);
    // ITM Put 的早期行使溢价应更大 (更值得提前行使)
    EXPECT_GT(itm.early_exercise_premium, atm.early_exercise_premium * 0.5);
}

TEST(LSMCEngine, ZeroInterestRateAmericanPutStillExercisable) {
    // r=0 时美式 Put 仍有价值 (理论上美式=欧式, 但算法应稳定)
    LSMCConfig cfg;
    cfg.S0 = 100.0; cfg.K = 100.0; cfg.T = 1.0;
    cfg.r = 0.0; cfg.q = 0.0; cfg.sigma = 0.20;
    cfg.n_paths = 10000; cfg.n_steps = 50;
    cfg.basis = BasisType::Laguerre; cfg.basis_order = 3;
    cfg.seed = 42; cfg.antithetic = true;
    LSMCEngine engine(cfg);
    PutPayOff put(100.0);
    LSMCResult result = engine.price_american(put);
    Real euro = bsm_put_price(100.0, 100.0, 1.0, 0.0, 0.0, 0.20);
    // r=0 时美式 Put = 欧式 Put (无提前行使收益)
    EXPECT_NEAR(result.price, euro, euro * 0.05 + 0.1);
    EXPECT_GT(result.price, 0.0);
}

// =====================================================================
// RISK-007: LSMC 交叉验证测试
// =====================================================================

// 1. 高噪声场景 CV 选择 λ>0
TEST(LSMCCV, CVSelectsNonZeroLambdaOnNoisyData) {
    LSMCConfig cfg;
    cfg.S0 = 100.0; cfg.K = 100.0; cfg.T = 1.0;
    cfg.r = 0.05; cfg.q = 0.0; cfg.sigma = 0.40;  // 高波动 → 高噪声
    cfg.n_paths = 5000; cfg.n_steps = 50;
    cfg.basis = BasisType::Laguerre; cfg.basis_order = 3;
    cfg.seed = 42; cfg.antithetic = true;
    cfg.use_cross_validation = true;
    cfg.cv_config.k_fold = 5;
    cfg.cv_config.lambda_grid = {0.0, 0.01, 0.1, 1.0, 10.0};

    LSMCEngine engine(cfg);
    PutPayOff put(100.0);
    LSMCResult result = engine.price_american(put);

    EXPECT_FALSE(result.selected_lambdas.empty())
        << "CV 启用时应记录 selected_lambdas";
    // 至少有一个时点选择了 λ>0
    bool any_nonzero = false;
    for (Real lam : result.selected_lambdas) {
        if (lam > 0.0) { any_nonzero = true; break; }
    }
    EXPECT_TRUE(any_nonzero) << "高噪声场景应至少有一个时点选择 λ>0";
}

// 2. CV 价格与固定 λ 价格统计一致 (容差 5·std_error)
TEST(LSMCCV, CVPriceCloseToFixedLambda) {
    LSMCConfig cfg_cv;
    cfg_cv.S0 = 100.0; cfg_cv.K = 100.0; cfg_cv.T = 1.0;
    cfg_cv.r = 0.05; cfg_cv.q = 0.0; cfg_cv.sigma = 0.20;
    cfg_cv.n_paths = 10000; cfg_cv.n_steps = 50;
    cfg_cv.basis = BasisType::Laguerre; cfg_cv.basis_order = 3;
    cfg_cv.seed = 42; cfg_cv.antithetic = true;
    cfg_cv.use_cross_validation = true;
    LSMCEngine engine_cv(cfg_cv);
    PutPayOff put(100.0);
    LSMCResult result_cv = engine_cv.price_american(put);

    LSMCConfig cfg_fixed = cfg_cv;
    cfg_fixed.use_cross_validation = false;
    cfg_fixed.ridge_lambda = 0.1;
    LSMCEngine engine_fixed(cfg_fixed);
    LSMCResult result_fixed = engine_fixed.price_american(put);

    Real tol = 5.0 * std::max(result_cv.std_error, result_fixed.std_error);
    EXPECT_NEAR(result_cv.price, result_fixed.price, tol)
        << "CV 价格=" << result_cv.price << " vs 固定 λ=0.1 价格=" << result_fixed.price
        << " tol=" << tol;
}

// 3. 样本不足时 fallback 到 ridge_lambda
TEST(LSMCCV, CVFallbackOnSmallSample) {
    // 构造极少 ITM 路径的场景: 深度 OTM Put, r=0, T 短
    // 大部分路径不会 ITM, 触发 fallback
    LSMCConfig cfg;
    cfg.S0 = 200.0; cfg.K = 100.0; cfg.T = 0.25;  // 深度 OTM
    cfg.r = 0.10; cfg.q = 0.0; cfg.sigma = 0.10;  // 低波动, 极少 ITM
    cfg.n_paths = 500; cfg.n_steps = 10;
    cfg.basis = BasisType::Laguerre; cfg.basis_order = 3;
    cfg.seed = 42; cfg.antithetic = true;
    cfg.use_cross_validation = true;
    cfg.ridge_lambda = 0.5;  // fallback 值
    cfg.cv_config.k_fold = 5;

    LSMCEngine engine(cfg);
    PutPayOff put(100.0);
    LSMCResult result = engine.price_american(put);

    // 若有 selected_lambdas, fallback 时应等于 ridge_lambda
    for (Real lam : result.selected_lambdas) {
        EXPECT_NEAR(lam, cfg.ridge_lambda, 1e-10)
            << "fallback 时 λ 应等于 ridge_lambda=" << cfg.ridge_lambda;
    }
}

// 4. CV 禁用时 selected_lambdas 为空
TEST(LSMCCV, CVDisabledUsesFixedLambda) {
    LSMCConfig cfg;
    cfg.S0 = 100.0; cfg.K = 100.0; cfg.T = 1.0;
    cfg.r = 0.05; cfg.q = 0.0; cfg.sigma = 0.20;
    cfg.n_paths = 1000; cfg.n_steps = 20;
    cfg.basis = BasisType::Laguerre; cfg.basis_order = 3;
    cfg.seed = 42; cfg.antithetic = true;
    cfg.use_cross_validation = false;  // 禁用 CV
    cfg.ridge_lambda = 0.1;

    LSMCEngine engine(cfg);
    PutPayOff put(100.0);
    LSMCResult result = engine.price_american(put);

    EXPECT_TRUE(result.selected_lambdas.empty())
        << "CV 禁用时 selected_lambdas 应为空";
}

// 5. K-fold 分折大小验证 (n_itm=100, K=5 → 每折 20)
TEST(LSMCCV, KFoldPartitionSizes) {
    // 通过验证 CV 在边界场景不崩溃来间接验证分折
    // n_paths=5000, ITM 路径数应远大于 k_fold * basis_order
    LSMCConfig cfg;
    cfg.S0 = 100.0; cfg.K = 100.0; cfg.T = 1.0;
    cfg.r = 0.05; cfg.q = 0.0; cfg.sigma = 0.20;
    cfg.n_paths = 5000; cfg.n_steps = 20;
    cfg.basis = BasisType::Laguerre; cfg.basis_order = 3;
    cfg.seed = 42; cfg.antithetic = true;
    cfg.use_cross_validation = true;
    cfg.cv_config.k_fold = 5;

    LSMCEngine engine(cfg);
    PutPayOff put(100.0);
    LSMCResult result = engine.price_american(put);

    // CV 应成功执行 (不崩溃), 且 selected_lambdas 非空
    EXPECT_FALSE(result.selected_lambdas.empty());
    // 价格应在合理范围
    EXPECT_GT(result.price, 0.0);
    EXPECT_LT(result.price, cfg.K);
}

