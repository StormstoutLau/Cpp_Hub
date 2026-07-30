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

# Phase 2 M1: models 测试 (A 站实现)
cpphub_add_test(test_heston_process unit/models/test_heston_process.cpp)

# Phase 2 M1: models/pricing 测试 (B 站实现)
cpphub_add_test(test_gbm unit/models/test_gbm.cpp)
cpphub_add_test(test_heston_cf unit/pricing/test_heston_cf.cpp)

# Phase 2 M2: PDE 引擎测试 (A 站实现)
cpphub_add_test(test_pde_engine unit/pricing/test_pde_engine.cpp)

# Phase 2 M2: 树形引擎测试 (B 站实现)
cpphub_add_test(test_tree_engine unit/pricing/test_tree_engine.cpp)

# Phase 2 M3: Sobol QMC 测试 (A 站实现)
cpphub_add_test(test_sobol_qmc unit/monte_carlo/test_sobol_qmc.cpp)

# Phase 2 M3: 方差缩减测试 (B 站实现)
cpphub_add_test(test_variance_reduction unit/monte_carlo/test_variance_reduction.cpp)

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
