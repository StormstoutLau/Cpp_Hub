include(GoogleTest)

# 测试宏: 添加单元测试
function(cpphub_add_test name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE cpphub_core GTest::gtest_main)
    gtest_discover_tests(${name})
endfunction()

# core 测试
cpphub_add_test(test_core_math unit/core/test_core_math.cpp)
cpphub_add_test(test_linalg unit/core/test_linalg.cpp)
cpphub_add_test(test_rng unit/core/test_rng.cpp)

# instruments 测试
cpphub_add_test(test_payoff_bridge unit/instruments/test_payoff_bridge.cpp)
cpphub_add_test(test_payoff_factory unit/instruments/test_payoff_factory.cpp)

# pricing 测试
cpphub_add_test(test_bsm_analytic unit/pricing/test_bsm_analytic.cpp)

# v1.4: 路径相关解析解测试 (Kemna-Vorst / Reiner-Rubinstein / Goldman-Sosin-Gatto)
cpphub_add_test(test_path_dependent_analytic unit/pricing/test_path_dependent_analytic.cpp)

# Phase 2 M1: models 测试 (A 站实现)
cpphub_add_test(test_heston_process unit/models/test_heston_process.cpp)

# Phase 2 M1: models/pricing 测试 (B 站实现)
cpphub_add_test(test_gbm unit/models/test_gbm.cpp)
cpphub_add_test(test_heston_cf unit/pricing/test_heston_cf.cpp)

# v1.3 P5 (Discovery 002): Heston CF 极端参数稳定性测试
cpphub_add_test(test_heston_cf_edge unit/pricing/test_heston_cf_edge.cpp)

# v1.3 F: CGMY + Kou Levy 过程特征函数测试
cpphub_add_test(test_levy_cf unit/pricing/test_levy_cf.cpp)

# Phase 2 M2: PDE 引擎测试 (A 站实现)
cpphub_add_test(test_pde_engine unit/pricing/test_pde_engine.cpp)

# v1.4: 2D PDE ADI 引擎测试 (Craig-Sneyd / Hundsdorfer-Verwer for Heston)
cpphub_add_test(test_pde_adi_2d unit/pricing/test_pde_adi_2d.cpp)

# v1.4: Chebyshev 配点谱方法测试 (1D BSM PDE, spectral convergence)
cpphub_add_test(test_chebyshev_pde unit/pricing/test_chebyshev_pde.cpp)

# v1.4: RBF-FD 多资产 PDE 测试 (1D BSM 极限 + 2D Margrabe + Basket vs MC)
cpphub_add_test(test_rbf_fd unit/pricing/test_rbf_fd.cpp)

# Phase 2 M2: 树形引擎测试 (B 站实现)
cpphub_add_test(test_tree_engine unit/pricing/test_tree_engine.cpp)

# Phase 2 M3: Sobol QMC 测试 (A 站实现)
cpphub_add_test(test_sobol_qmc unit/monte_carlo/test_sobol_qmc.cpp)

# Phase 2 M3: 方差缩减测试 (B 站实现)
cpphub_add_test(test_variance_reduction unit/monte_carlo/test_variance_reduction.cpp)

# 高级方差缩减方法测试 (Importance Sampling / Conditional MC / Stratified Sampling)
cpphub_add_test(test_variance_reduction_advanced unit/monte_carlo/test_variance_reduction_advanced.cpp)

# v1.4 统一方差缩减框架测试 (Antithetic 独立类 + VREngine + IS Put + Z 复用)
cpphub_add_test(test_unified_vr unit/monte_carlo/test_unified_vr.cpp)

# Phase 2 M3: 端到端集成测试 (B 站实现)
cpphub_add_test(test_integration_phase2 unit/monte_carlo/test_integration_phase2.cpp)

# autodiff 库编译验证 (非测试, 仅验证头文件可用)
add_executable(test_autodiff_compile unit/test_autodiff_compile.cpp)
target_link_libraries(test_autodiff_compile PRIVATE cpphub_core)

# Phase 3 M1: AAD Greeks 测试 (A 站实现, 反向模式 autodiff::var)
cpphub_add_test(test_aad_greeks unit/risk/greeks/test_aad_greeks.cpp)

# Phase 3 M1: 解析/数值/前向AD Greeks 测试 (B 站实现, 前向模式 autodiff::dual)
cpphub_add_test(test_greeks_analytic_numerical unit/risk/greeks/test_greeks_analytic_numerical.cpp)

# Phase 3 M1: Pathwise Greeks 测试 (主控站实现, 光滑 payoff 路径法)
cpphub_add_test(test_pathwise_greeks unit/risk/greeks/test_pathwise_greeks.cpp)

# Phase 3 M1: LR Greeks 测试 (主控站实现, 不连续 payoff 似然比法)
cpphub_add_test(test_lr_greeks unit/risk/greeks/test_lr_greeks.cpp)

# Phase 3 M1: Greeks Factory 测试 (主控站实现, Auto/Analytic/Pathwise/LR/FD/AAD 分发)
cpphub_add_test(test_greeks_factory unit/risk/greeks/test_greeks_factory.cpp)

# Phase 3 M2: VaR/ES 引擎测试 (A 站实现, 历史/参数化/MC/ES)
cpphub_add_test(test_var_engine unit/risk/var/test_var_engine.cpp)

# Phase 3 M2: VaR 回测测试 (B 站实现, Kupiec/Christoffersen/Basel)
cpphub_add_test(test_backtesting unit/risk/var/test_backtesting.cpp)

# Phase 3 M2: 场景/压力测试 (B 站实现, StressTester/SensitivityAnalysis)
cpphub_add_test(test_stress_sensitivity unit/risk/scenario/test_stress_sensitivity.cpp)

# Phase 3 M3: 标定框架测试 (A 站实现, LM/NelderMead/DE + Heston/SABR Calibrator)
cpphub_add_test(test_calibration_framework unit/calibration/test_calibration_framework.cpp)

# E 任务: 多模型标定器测试 (Bates/VG/CEV Calibrator, DE + LM + COS/解析定价)
cpphub_add_test(test_multi_calibrator unit/calibration/test_multi_calibrator.cpp)

# Phase 3 M3: 波动率曲面测试 (B 站实现, SVI/SSVI/VolSurface)
cpphub_add_test(test_vol_surface unit/models/vol_surface/test_vol_surface.cpp)

# Phase 3 M3: Dupire 局部波动率 + 端到端集成测试 (主控站实现)
cpphub_add_test(test_dupire_local_vol unit/models/vol_surface/test_dupire_local_vol.cpp)
cpphub_add_test(test_integration_phase3 unit/integration/test_integration_phase3.cpp)

# Phase 4 LITE - D2 整改项: SSVI 跨期限模块测试 (主控站实现)
cpphub_add_test(test_ssvi unit/models/vol_surface/test_ssvi.cpp)

# 临时诊断测试 (RISK-013 SVI 标定失败排查)
cpphub_add_test(test_svi_diagnostic unit/calibration/test_svi_diagnostic.cpp)

# Phase 4 LITE - G2 整改项: Python 交叉验证测试 (主控站实现)
# 基准值由 tests/validation/python/ 下的脚本生成 (numpy/scipy/statsmodels)
cpphub_add_test(test_python_cross_validation unit/validation/test_python_cross_validation.cpp)

# v1.1: LSMC 美式期权引擎测试 (Longstaff-Schwartz 2001)
cpphub_add_test(test_lsmc unit/pricing/test_lsmc.cpp)

# v1.1: 利率模型测试 (Vasicek/CIR/Hull-White/G2++ 仿射短期利率模型)
cpphub_add_test(test_short_rate unit/models/test_short_rate.cpp)

# v1.1: 通用傅里叶引擎测试 (COS Fang-Oosterlee + Carr-Madan FFT + CharFn 接口)
cpphub_add_test(test_fourier_engine unit/pricing/test_fourier_engine.cpp)

# v1.1: 利率衍生品测试 (IRS / Cap-Floor / Swaption)
cpphub_add_test(test_ir_derivs unit/instruments/test_ir_derivs.cpp)

# v1.2: HJM 远期利率框架 + LMM (BGM) 市场模型测试
cpphub_add_test(test_hjm_lmm unit/models/test_hjm_lmm.cpp)

# v1.2: 多资产 MC + 路径相关/篮子期权测试
cpphub_add_test(test_multi_asset_mc unit/pricing/test_multi_asset_mc.cpp)

# v1.2 Batch 2: OIS / 隔夜利率衍生品 (ZeroCurve + OIS + FRA + Basis Swap + Dual-curve IRS)
cpphub_add_test(test_ois unit/instruments/test_ois.cpp)

# v1.2 Batch 3: xVA (CVA / DVA / FVA / BVA) 框架
cpphub_add_test(test_xva unit/risk/test_xva.cpp)

# v1.2 Batch 4: 信用违约互换 (CDS) + CDS curve bootstrap
cpphub_add_test(test_cds unit/instruments/test_cds.cpp)

# v1.2 Batch 5: PFE (潜在未来暴露) + SA-CCR (Basel III 交易对手信用风险)
cpphub_add_test(test_pfe_sa_ccr unit/risk/test_pfe_sa_ccr.cpp)

# v1.2 Batch 6: Total Return Swap (TRS) - 总回报互换
cpphub_add_test(test_trs unit/instruments/test_trs.cpp)

# v1.2 Batch 7: 信用利差期权 (Credit Spread Option, Black 1976)
cpphub_add_test(test_credit_spread_option unit/instruments/test_credit_spread_option.cpp)

# v1.2 Batch 8: Copula 模型 (Gaussian / One-Factor / t) + 合成 CDO + 篮子 CDS
cpphub_add_test(test_copula unit/instruments/test_copula.cpp)

# v1.3 D: CVA with Wrong-Way Risk (WWR) — Hull-White / Gaussian Copula / MC
cpphub_add_test(test_wrong_way_risk unit/risk/test_wrong_way_risk.cpp)
cpphub_add_test(test_cdo unit/instruments/test_cdo.cpp)
cpphub_add_test(test_basket_cds unit/instruments/test_basket_cds.cpp)

# v1.2 Batch 9: SABR 模型 (Hagan IV 闭式公式 + 随机过程)
cpphub_add_test(test_sabr unit/models/test_sabr.cpp)

# v1.2 Batch 10: CEV 模型 (非中心卡方解析定价 + Euler/log-Euler 随机过程)
cpphub_add_test(test_cev unit/models/test_cev.cpp)

# v1.2 Batch 11: Bates 模型 (Heston + Merton 跳跃, CF 乘积 + Full Truncation Euler + Poisson 跳跃)
cpphub_add_test(test_bates unit/models/test_bates.cpp)

# v1.2 Batch 12: Variance Gamma 模型 (Madan-Carr-Chang 1998, 纯跳跃 Levy 过程, CF + Gamma 时变 Brownian + COS 定价)
cpphub_add_test(test_vg unit/models/test_vg.cpp)

# v1.2 Batch 13: Rough Bergomi 模型 (rBergomi) — 粗糙波动率
# (a) 解析层: 近似特征函数 + 累积量展开 (Gatheral-Jaisson-Rosenbaum 2018)
cpphub_add_test(test_rough_bergomi_analytic unit/models/test_rough_bergomi_analytic.cpp)
# (b) 过程层: RL-fBm Hybrid Scheme 采样器 (Bennedsen-Lunde-Pakkanen 2017, O(N·b) per path)
cpphub_add_test(test_rbergomi_hybrid unit/models/test_rbergomi_hybrid.cpp)

# v1.3 G: Rough Heston CF (分数阶 Riccati) + 分数阶 Adams 求解器
cpphub_add_test(test_fractional_adams unit/numerics/test_fractional_adams.cpp)
cpphub_add_test(test_rough_heston_cf unit/pricing/test_rough_heston_cf.cpp)

# v1.3 H: Rough Heston 过程层 (分数阶 Euler + MC 定价, H=0.5 退化为 Heston)
cpphub_add_test(test_rough_heston_process unit/models/test_rough_heston_process.cpp)

# Phase 4 LITE - ADR-009: C ABI 测试 (条件编译, 需 CPPHUB_ENABLE_C_API=ON)
if(CPPHUB_ENABLE_C_API)
    add_executable(test_c_api unit/c_api/test_c_api.cpp)
    target_link_libraries(test_c_api PRIVATE cpphub_core cpphub_c_api GTest::gtest_main)
    gtest_discover_tests(test_c_api)
endif()

# v1.4.0 Phase 5: 高频计量经济学模块 (HFE) 测试
# 对标 R highfrequency 1.0.3 (Boudt, Kleen, Sjørup 2022, JSS doi:10.18637/jss.v104.i08)
# 基准数据: tests/fixtures/hfe/baselines.json (由 generate_r_baselines.R 生成)
cpphub_add_test(test_hfe_realized_measures unit/hfecon/test_realized_measures.cpp)

# v1.4.1 Phase 5: HFE Realized Kernel (BNS 2008 ECTA) + 噪声方差 + 最优带宽
# 对标 R highfrequency 1.0.3 rKernelCov (KK() + kernelEstimator() 源码实测)
# 基准数据: tests/fixtures/hfe/generate_v141_baselines.R 生成
cpphub_add_test(test_realized_kernel unit/hfecon/test_realized_kernel.cpp)

# v1.4.2 Phase 5 Wave A: HFE Hayashi-Yoshida 协方差 + makePsd 工具
# 对标 R highfrequency 1.0.3 rHYCov + pcovcc + makePsd
cpphub_add_test(test_hayashi_yoshida_cov unit/hfecon/test_hayashi_yoshida_cov.cpp)

# v1.4.2 Phase 5 Wave B: HFE Two-Scale 协方差 (rTSCov + TSRV + TSCov_bi)
# 对标 R highfrequency 1.0.3 rTSCov
cpphub_add_test(test_two_scale_cov unit/hfecon/test_two_scale_cov.cpp)

# v1.4.2 Phase 5 Wave B: HFE Modulated Realized 协方差 (rMRCov)
# 对标 R highfrequency 1.0.3 rMRCov
cpphub_add_test(test_modulated_realized_cov unit/hfecon/test_modulated_realized_cov.cpp)

# v1.4.2 Phase 5 Wave B: HFE Pre-averaged 协方差 (rAVGCov)
# 对标 R highfrequency 1.0.3 rAVGCov (preAveragingReturnsInternal + rAVGCov_bi)
cpphub_add_test(test_preaveraged_cov unit/hfecon/test_preaveraged_cov.cpp)

# v1.4.2 Phase 5 Wave B: HFE Robust Two-Scale 协方差 (rRTSCov + RTSRV)
# 对标 R highfrequency 1.0.3 rRTSCov (RTSRV + RTSCov_bi + inv_pchisq3)
cpphub_add_test(test_robust_two_scale_cov unit/hfecon/test_robust_two_scale_cov.cpp)

# v1.4.2 Phase 5 Wave C: HFE HAR 模型 (Corsi 2009 + BPQ 2016 RQ 变换)
# 对标 R highfrequency 1.0.3 HARmodel (HAR/HARJ/HARQ/CHAR + log/sqrt 变换)
cpphub_add_test(test_har_model unit/hfecon/test_har_model.cpp)

# v1.4.2 Phase 5 Wave C: HFE HEAVY 模型 (Shephard & Sheppard 2010)
# 对标 R highfrequency 1.0.3 HEAVYmodel (calcRecVarEq + heavyLLH + NelderMead MLE)
cpphub_add_test(test_heavy_model unit/hfecon/test_heavy_model.cpp)

# v1.4.3 Phase 5 Wave D: HFE Spread Cleaner (rmLargeSpread/rmNegativeSpread/spreadPrices)
# 对标 R highfrequency 1.0.3 dataHandling.R L1617/1670/3019
cpphub_add_test(test_spread_cleaner unit/hfecon/test_spread_cleaner.cpp)

# v1.4.3 Phase 5 Wave D: HFE Amihud Illiquidity (Amihud 2002, 无 R 对照)
cpphub_add_test(test_amihud unit/hfecon/test_amihud.cpp)

# v1.4.3 Phase 5 Wave D: HFE Liquidity Measures (23 种 + getTradeDirection)
# 对标 R highfrequency 1.0.3 liquidityMeasures.R L231/346 (排幻觉 D1-D3)
cpphub_add_test(test_liquidity_measures unit/hfecon/test_liquidity_measures.cpp)

# v1.4.3 Phase 5 Wave D: HFE AJ Jump Test (Aït-Sahalia & Jacod 2009)
# 对标 R highfrequency 1.0.3 jumpTests.R L106 + internalJumpTests.R (排幻觉 D4-D9)
cpphub_add_test(test_aj_jump_test unit/hfecon/test_aj_jump_test.cpp)

# v1.4.3 Phase 5 Wave D: HFE JO Jump Test (Jiang & Oomen 2008)
# 对标 R highfrequency 1.0.3 jumpTests.R L446 + internalJumpTests.R + internals.cpp L207 (排幻觉 D10-D13)
cpphub_add_test(test_jo_jump_test unit/hfecon/test_jo_jump_test.cpp)

# v1.4.3 Phase 5 Wave D: HFE Intraday Jump Test (Lee & Mykland 2008)
# 对标 R highfrequency 1.0.3 jumpTests.R L583 + internalSpotVolAndDrift.R L940 (排幻觉 D14-D16)
cpphub_add_test(test_intraday_jump_test unit/hfecon/test_intraday_jump_test.cpp)

# v1.4.3 Phase 5 Wave D: HFE Rank Jump Test (Bollerslev & Todorov 2011)
# 对标 R highfrequency 1.0.3 jumpTests.R L976 + internalJumpTests.R L115/125/149 (排幻觉 D17-D23)
cpphub_add_test(test_rank_jump_test unit/hfecon/test_rank_jump_test.cpp)

# === Phase 6 v1.5 M1: 经典参数计量模块 (ADR-013 双层 linalg) ===
# test_linalg_dynamic: Eigen3 封装层 (MatrixXD/VectorXD/SVD/LLT/LU/QR)
# 注: 此测试链接 cpphub_econometrics (含 Eigen3), 不污染 cpphub_core
add_executable(test_linalg_dynamic unit/econometrics/test_linalg_dynamic.cpp)
target_link_libraries(test_linalg_dynamic PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_linalg_dynamic)

# === Phase 6 v1.5 M1: econometrics/core 基础类型 (任务 1.2-1.5) ===
add_executable(test_econometrics_core unit/econometrics/test_econometrics_core.cpp)
target_link_libraries(test_econometrics_core PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_econometrics_core)

# === Phase 6 v1.5 M1: econometrics/inference HAC 内核 (任务 1.8, 排幻觉 E4/E5) ===
add_executable(test_hac_kernels unit/econometrics/test_hac_kernels.cpp)
target_link_libraries(test_hac_kernels PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_hac_kernels)

# === Phase 6 v1.5 M1: OLS + HC0-HC5 (任务 1.6-1.7, 排幻觉 E1/E2/E3) ===
# test_ols_hc: OLSEstimator + compute_hc_vcov (手算 N=4,K=2 解析值, 容差 1e-12)
add_executable(test_ols_hc unit/econometrics/test_ols_hc.cpp)
target_link_libraries(test_ols_hc PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_ols_hc)

# === Phase 6 v1.5 M1: Newey-West HAC 协方差 (任务 1.9, 排幻觉 E4/E5) ===
add_executable(test_newey_west unit/econometrics/test_newey_west.cpp)
target_link_libraries(test_newey_west PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_newey_west)

# === Phase 6 v1.5 M1: Cluster SE (任务 1.10, 排幻觉 E6) ===
# test_cluster_se: 单向聚类 (Liang-Zeger 1986) + 双向聚类 (CGM 2011)
# 手算 N=4,K=1,K=2 解析值, 容差 1e-12
add_executable(test_cluster_se unit/econometrics/test_cluster_se.cpp)
target_link_libraries(test_cluster_se PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_cluster_se)

# === Phase 6 v1.5 M1: 面板数据工具 (任务 1.11) ===
# test_panel_data: detect_balance / extract_cluster_ids / reshape_y_to_wide / reshape_x_to_wide
add_executable(test_panel_data unit/econometrics/test_panel_data.cpp)
target_link_libraries(test_panel_data PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_panel_data)

# === Phase 6 v1.5 M1: Cluster SE 集成测试 (任务 1.10 集成, vs R sandwich::vcovCL) ===
# test_cluster_se_integration: Grunfeld 数据 (N=200, K=3) vs R sandwich 3.1-3 vcovCL type="HC1"
# 数据来源: tests/fixtures/econometrics/grunfeld_data.hpp (auto-generated by generate_grunfeld_baseline.R)
# baseline: tests/fixtures/econometrics/grunfeld_cluster_baseline.json (digits=17 全精度)
# 容差: 1e-8 (spec §3.5 M1 验收标准)
add_executable(test_cluster_se_integration unit/econometrics/test_cluster_se_integration.cpp)
target_include_directories(test_cluster_se_integration PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/econometrics)
target_link_libraries(test_cluster_se_integration PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_cluster_se_integration)

# === Phase 6 v1.5 M1: 端到端集成测试 (任务 M1 集成) ===
# test_integration_m1: End-to-end pipeline (EconData → OLS → HC/HAC/Cluster → EstimationResult)
# 包含 HAR 风格时间序列回归 + Newey-West HAC (Corsi 2009 场景)
add_executable(test_integration_m1 unit/econometrics/test_integration_m1.cpp)
target_link_libraries(test_integration_m1 PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_integration_m1)

# === Phase 6 v1.5 M2: 信息准则 (任务 2.4, AIC/BIC/HQ/AICc) ===
# test_diagnostics: 手算解析值 + statsmodels 对照 (容差 1e-10)
add_executable(test_diagnostics unit/econometrics/test_diagnostics.cpp)
target_link_libraries(test_diagnostics PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_diagnostics)

# === Phase 6 v1.5 M2: 假设检验 (任务 2.3, Wald/LR/LM/J-test) ===
# test_hypothesis_tests: 手算解析值 + χ²/F 分布临界值对照 (容差 1e-6)
# 排幻觉点 E9 (Wald F/χ² 双模式) / H1 (LM=N·R²_aux) / H2 (J-test df=q-k)
add_executable(test_hypothesis_tests unit/econometrics/test_hypothesis_tests.cpp)
target_link_libraries(test_hypothesis_tests PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_hypothesis_tests)

# === Phase 6 v1.5 M2: MLE Gaussian (任务 2.1-2.2, 排幻觉点 G1-G5) ===
# test_mle_gaussian: 手算 N=4,K=2 数据集 + OLS 闭式解对齐 + 三种协方差公式验证
# 排幻觉点 G3 (Hessian V = σ²·(X'X)^{-1}, 非 (X'X)^{-1})
# 排幻觉点 G4 (Sandwich V = HC0 White 1980, meat 用 raw residual 无 1/σ²)
add_executable(test_mle_gaussian unit/econometrics/test_mle_gaussian.cpp)
target_link_libraries(test_mle_gaussian PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_mle_gaussian)

# === Phase 6 v1.5 M2: MLE Logistic (任务 2.1-2.2, 排幻觉点 E7/E8/L1-L3) ===
# test_mle_logistic: link/mean 互逆 + Newton-Raphson 一步手算 + 对称数据 MLE 解析解
# 排幻觉点 E7 (R glm IRLS vs C++ Newton-Raphson, canonical link 等价)
# 排幻觉点 E8 (GLM bread = (X'WX)^{-1}, meat = X' diag(ε²) X)
add_executable(test_mle_logistic unit/econometrics/test_mle_logistic.cpp)
target_link_libraries(test_mle_logistic PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_mle_logistic)

# === Phase 6 v1.5 M2: MLE Poisson (任务 2.1-2.2, 排幻觉点 P1-P5) ===
# test_mle_poisson: link/mean 互逆 + IRLS weight = μ + 含截距 K=2 手算解析解
# 排幻觉点 P2 (canonical log link: w = μ, V(μ) = μ)
# 排幻觉点 P3 (ℓ_i = y·η - μ - log(y!), 用 lgamma(y+1) 计算 log(y!))
add_executable(test_mle_poisson unit/econometrics/test_mle_poisson.cpp)
target_link_libraries(test_mle_poisson PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_mle_poisson)

# === Phase 6 v1.5 M2: MLE Probit/NegativeBinomial (任务 2.1-2.2, 排幻觉点 PB1-PB4, NB1-NB3) ===
# test_mle_probit_nb: Probit link 互逆 + IRLS weight = 2/π (μ=0.5) + 与 Logistic 系数关系 (√3/π)
# 排幻觉点 PB4 (β_probit ≈ β_logistic × √3/π, 因 logistic σ=π/√3 比正态宽)
# 排幻觉点 NB1 (w = μ/(1+αμ), α 增大 w 减小)
# 排幻觉点 NB2 (α→0 时 w→μ, 退化为 Poisson)
add_executable(test_mle_probit_nb unit/econometrics/test_mle_probit_nb.cpp)
target_link_libraries(test_mle_probit_nb PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_mle_probit_nb)

# === Phase 6 v1.5 M2: 端到端集成测试 (任务 2.5, 工厂 + Wald/LR + AIC/BIC) ===
# test_integration_m2: EconData → 工厂创建估计器 → 估计 → Wald/LR 检验 → AIC/BIC 模型选择
# 排幻觉点 F1 (工厂全新实例) / F2 (未注册抛异常) / F3 (unique_ptr 转移所有权) / F4 (Meyers Singleton)
# 排幻觉点 E7 (Newton-Raphson ≡ IRLS) / E8 (Sandwich bread/meat) / E9 (Wald χ²+F) / G3 (Gaussian σ² 因子)
add_executable(test_integration_m2 unit/econometrics/test_integration_m2.cpp)
target_link_libraries(test_integration_m2 PRIVATE cpphub_econometrics GTest::gtest_main)
gtest_discover_tests(test_integration_m2)
