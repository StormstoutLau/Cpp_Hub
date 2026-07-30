cmake_minimum_required(VERSION 3.25)
project(CppHub VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# MSVC UTF-8 (仅 CXX, CUDA 通过 -Xcompiler 传递需单独处理)
if(MSVC)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/utf-8> $<$<COMPILE_LANGUAGE:CXX>:/W3>)
    add_compile_definitions(_USE_MATH_DEFINES)
endif()

# Release 优化 (禁用 -ffast-math 保证位精确复现)
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(MSVC)
        add_compile_options(/O2 /arch:AVX2 /fp:precise)
    else()
        add_compile_options(-O3 -march=x86-64-v3 -ffp-contract=off)
    endif()
endif()

# 头文件库
add_library(cpphub_core INTERFACE)
target_include_directories(cpphub_core INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party  # autodiff header-only library
)
target_compile_features(cpphub_core INTERFACE cxx_std_20)

# GoogleTest (FetchContent)
include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
)
set(gtest_force_shared_runtime ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# 测试
enable_testing()
add_subdirectory(tests)

# Phase 4 LITE - 性能基准 (不纳入 ctest 默认运行, 通过 `make bench` 触发)
add_subdirectory(benchmarks)

# Phase 4 LITE - M3: GPU MC (CUDA, 主控站 RTX 4060, A/B 站无 CUDA 时自动跳过)
# 启用方式: cmake -DCPPHUB_ENABLE_CUDA=ON ..
# 默认 OFF, A/B 站 GCC 无 CUDA 编译器时会自动跳过
option(CPPHUB_ENABLE_CUDA "Build CUDA GPU MC kernel (requires NVIDIA GPU + nvcc)" OFF)
if(CPPHUB_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        # RTX 4060 = Ada Lovelace = compute capability 8.9
        # 允许用户覆盖: -DCMAKE_CUDA_ARCHITECTURES=89 (默认) 或 86/90 等
        if(NOT CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES 89)
        endif()
        # CUDA 标准与 C++ 标准对齐
        if(NOT CMAKE_CUDA_STANDARD)
            set(CMAKE_CUDA_STANDARD 20)
        endif()
        message(STATUS "CUDA GPU MC: ENABLED (nvcc=${CMAKE_CUDA_COMPILER}, arch=${CMAKE_CUDA_ARCHITECTURES})")
        set(CPPHUB_CUDA_AVAILABLE TRUE)
    else()
        message(WARNING
            "CPPHUB_ENABLE_CUDA=ON but no CUDA compiler found. "
            "Install CUDA Toolkit (https://developer.nvidia.com/cuda-toolkit). "
            "Falling back to CPU stub.")
        set(CPPHUB_CUDA_AVAILABLE FALSE)
    endif()
else()
    message(STATUS "CUDA GPU MC: DISABLED (set CPPHUB_ENABLE_CUDA=ON to enable)")
    set(CPPHUB_CUDA_AVAILABLE FALSE)
endif()

# GPU MC 子目录 (CUDA 可用时编译 .cu, 不可用时编译 CPU stub)
add_subdirectory(src/performance/gpu)

# Phase 4 LITE - ADR-009: C ABI 稳定边界 (可选, 默认 OFF)
# 启用方式: cmake -DCPPHUB_ENABLE_C_API=ON ..
# 编译静态库 cpphub_c_api, 暴露 cpphub_v1_* 版本化符号
option(CPPHUB_ENABLE_C_API "Build C ABI library for cross-language bindings" OFF)
if(CPPHUB_ENABLE_C_API)
    message(STATUS "C ABI library: ENABLED")
    add_subdirectory(src/c_api)
else()
    message(STATUS "C ABI library: DISABLED (set CPPHUB_ENABLE_C_API=ON to enable)")
endif()

# Phase 4 LITE - M2: Python 绑定 (nanobind, 可选)
# 启用方式: cmake -DCPPHUB_ENABLE_PYTHON=ON ..
#           或 pip install . (自动设置)
# 依赖: nanobind (pip install nanobind), Python 开发头文件
option(CPPHUB_ENABLE_PYTHON "Build Python bindings via nanobind" OFF)
if(CPPHUB_ENABLE_PYTHON)
    # 寻找 Python
    find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

    # 检查 nanobind 是否可用
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.cmake_dir())"
        OUTPUT_VARIABLE _nanobind_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _nanobind_result
    )
    if(_nanobind_result EQUAL 0)
        message(STATUS "Python bindings: ENABLED (nanobind found at ${_nanobind_dir})")
        add_subdirectory(python)
    else()
        message(WARNING
            "CPPHUB_ENABLE_PYTHON=ON but nanobind not found. "
            "Install with: pip install nanobind. Skipping Python bindings.")
    endif()
else()
    message(STATUS "Python bindings: DISABLED (set CPPHUB_ENABLE_PYTHON=ON to enable)")
endif()
