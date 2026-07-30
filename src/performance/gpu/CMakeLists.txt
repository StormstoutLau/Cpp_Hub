# Phase 4 LITE - M3: GPU MC 子目录构建脚本
#
# 行为:
#   - CPPHUB_ENABLE_CUDA=ON 且 nvcc 可用 → 编译 gpu_mc.cu 为 cpphub_cuda 共享库
#   - CPPHUB_ENABLE_CUDA=OFF 或无 nvcc → 编译 gpu_mc_cpu_stub.cpp 为 cpphub_cuda 静态库
#
# 测试: 仅在 CUDA 可用时构建 test_gpu_mc (主控站验证, 不入 A/B 站回归)
# 基准: 仅在 CUDA 可用时构建 bench_gpu_mc

set(GPU_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
# src/performance/gpu/ → 3 层上到根目录
set(GPU_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../include)
set(GPU_THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../third_party)
set(GPU_TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../tests)

# ============================================================================
# 1. cpphub_cuda 库 (CUDA 或 CPU stub)
# ============================================================================
if(CPPHUB_CUDA_AVAILABLE)
    # === CUDA 实现路径 ===
    # STATIC 库避免 Windows DLL 导出/导入复杂性, 研究项目无需动态加载
    add_library(cpphub_cuda STATIC
        gpu_mc.cu
        gpu_config.cpp
    )

    target_sources(cpphub_cuda PRIVATE
        FILE_SET cuda_headers
        TYPE HEADERS
        BASE_DIRS ${GPU_INCLUDE_DIR}
        FILES
            ${GPU_INCLUDE_DIR}/cpphub/performance/gpu/gpu_config.hpp
            ${GPU_INCLUDE_DIR}/cpphub/performance/gpu/gpu_mc.hpp
    )

    target_include_directories(cpphub_cuda PUBLIC
        ${GPU_INCLUDE_DIR}
        ${GPU_THIRD_PARTY_DIR}
    )

    target_compile_features(cpphub_cuda PUBLIC cxx_std_20)

    # CUDA 编译选项 (默认禁用 fast-math 以保证位精确, 仅在 Release 启用 -O3)
    target_compile_options(cpphub_cuda PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:
            -Xptxas=-v           # PTX 优化日志
            -lineinfo             # 行号信息 (nsight 调试用)
            # 不加 --use_fast_math, 默认即禁用, 保证 FP64 位精确
        >
    )

    # MSVC 主机编译器 UTF-8 (CXX 和 CUDA host code 均需 /utf-8, 否则中文注释导致 C1020)
    if(MSVC)
        target_compile_options(cpphub_cuda PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:/utf-8 /W3>
            $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/utf-8 -Xcompiler=/W3>
        )
        target_compile_definitions(cpphub_cuda PRIVATE _USE_MATH_DEFINES)
    endif()

    target_compile_definitions(cpphub_cuda PUBLIC CPPHUB_HAS_CUDA=1)

    message(STATUS "  cpphub_cuda: SHARED library (CUDA implementation)")

else()
    # === CPU stub 实现路径 (A/B 站或无 CUDA 环境) ===
    add_library(cpphub_cuda STATIC
        gpu_mc_cpu_stub.cpp
        gpu_config.cpp
    )

    target_include_directories(cpphub_cuda PUBLIC
        ${GPU_INCLUDE_DIR}
        ${GPU_THIRD_PARTY_DIR}
    )

    target_compile_features(cpphub_cuda PUBLIC cxx_std_20)

    if(MSVC)
        target_compile_options(cpphub_cuda PRIVATE /utf-8 /W3)
        target_compile_definitions(cpphub_cuda PRIVATE _USE_MATH_DEFINES)
    endif()

    target_compile_definitions(cpphub_cuda PUBLIC CPPHUB_HAS_CUDA=0)

    message(STATUS "  cpphub_cuda: STATIC library (CPU stub, GPU calls will throw)")

endif()

# ============================================================================
# 2. GPU MC 测试 (仅 CUDA 可用时构建)
# ============================================================================
if(CPPHUB_CUDA_AVAILABLE)
    add_executable(test_gpu_mc
        ${GPU_TESTS_DIR}/unit/performance/test_gpu_mc.cpp
    )
    target_link_libraries(test_gpu_mc PRIVATE
        cpphub_cuda
        cpphub_core
        GTest::gtest_main
    )
    target_compile_features(test_gpu_mc PRIVATE cxx_std_20)
    if(MSVC)
        target_compile_options(test_gpu_mc PRIVATE /utf-8 /W3)
        target_compile_definitions(test_gpu_mc PRIVATE _USE_MATH_DEFINES)
    endif()

    # 不纳入默认 ctest (主控站单独验证, 避免在 A/B 站触发)
    # 手动运行: ctest -R gpu_mc --verbose
    add_test(NAME gpu_mc COMMAND test_gpu_mc)

    message(STATUS "  test_gpu_mc: ENABLED")
else()
    message(STATUS "  test_gpu_mc: SKIPPED (no CUDA)")
endif()
