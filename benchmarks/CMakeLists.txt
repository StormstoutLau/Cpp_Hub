# Phase 4 LITE - 性能基准 (补齐 Phase 3 审计 E1/E2 整改项)
# 基准不纳入 ctest 默认运行,通过 `make bench` 或直接运行可执行文件触发
#
# E1: 万期权 AAD Greeks 单次扫描 (目标 < 100ms)
# E2: 10000 场景 MC VaR Full Revaluation (目标 < 500ms)

# 基准宏: 添加基准可执行文件 (不注册到 ctest)
function(cpphub_add_benchmark name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE cpphub_core)
    # 基准默认在 Release 模式下编译,继承顶层 CMake 的 Release 优化选项
endfunction()

# E1: 万期权 AAD Greeks 批量基准
cpphub_add_benchmark(bench_greeks_batch bench_greeks_batch.cpp)

# E2: 10000 场景 MC VaR 基准
cpphub_add_benchmark(bench_var_mc bench_var_mc.cpp)

# PDE 性能基准 (PHASE3_SPEC §6: PDE 400x2000 < 100ms)
cpphub_add_benchmark(bench_pde bench_pde.cpp)

# 便捷目标: 编译并运行所有基准
add_custom_target(bench
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    COMMENT "Running performance benchmarks (not ctest unit tests)"
    DEPENDS bench_greeks_batch bench_var_mc bench_pde
)
