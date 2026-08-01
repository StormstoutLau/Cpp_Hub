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
