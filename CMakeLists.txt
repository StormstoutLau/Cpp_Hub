cmake_minimum_required(VERSION 3.25)
# RISK-032: 版本单一源 = 仓库根 VERSION
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" CPPHUB_LIB_VERSION LIMIT_COUNT 1)
project(CppHub VERSION ${CPPHUB_LIB_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# MSVC UTF-8 (仅 CXX, CUDA 通过 -Xcompiler 传递需单独处理)
if(MSVC)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/utf-8> $<$<COMPILE_LANGUAGE:CXX>:/W3>
        # C1128 根治: COFF .obj 默认 65,279 节上限; Debug (/Zi) 下每个 COMDAT
        # 函数额外带调试节, 本库 5 个巨型 TU (c_api + VECM + 2 集成 +
        # tool_registry) 实测超限。 /bigobj 仅把节容量扩到 2^32, 不改代码
        # 生成/ABI; VS v140+ 与 /Zi 及 /ZI 兼容。 全配置启用 (非仅 Debug):
        # $<CONFIG:> genex 在 VS 生成器下有静默丢弃先例, 且配置间格式一致
        # 无收益; COMPILE_LANGUAGE genex 与现有 /utf-8 同款, 避免传给 nvcc。
        $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
    add_compile_definitions(_USE_MATH_DEFINES)
endif()

# 浮点语义旗标 (禁用 -ffast-math 保证位精确复现; 与 Debug 的 /RTC1 / -O0 共存无害)
# 优化级别不在此注入, 交给 CMake per-config 默认 flags (MSVC Release=/O2 /Ob2,
# GCC Release=-O3): 全局 add_compile_options 会连带 Debug 配置 — MSVC 下 /O2 与
# /RTC1 冲突 (D8016 实测); CMAKE_BUILD_TYPE 判断对多配置生成器 (VS) 不适用
# (configure 期 -D CMAKE_BUILD_TYPE=Release 会让 if 对所有 per-config 全局 firing);
# $<CONFIG:> genex 在 add_compile_options 中也被 VS 生成器静默丢弃 (实测 0 落点)。
if(MSVC)
    add_compile_options(/fp:precise)
else()
    add_compile_options(-ffp-contract=off)
endif()

# RISK-029: x86-64-v3 / AVX2 代码生成改为显式 opt-in (默认 OFF, 基线 x86-64-v1)。
# 旧实现无条件在 Release 全局追加 /arch:AVX2 (MSVC) 或 -march=x86-64-v3 (GCC/Clang):
#   - 不支持 AVX2 的 CPU 上, 静态初始化器先于 main 执行 AVX2 指令 → SIGILL
#     (社区先例: TF#17411 / LanceDB / manylinux#1725 / cryptography#12069);
#   - aarch64 上 GCC 直接报 unknown -march, 配置即败。
# opt-in 后探测编译器支持, 拒绝即 FATAL_ERROR (显式选择不应静默降级);
# CheckCXXCompilerFlag 只证 "编译器能生成", 不证 "目标 CPU 能执行" —
# 运行时能力检测唯一真源仍是 config.hpp 的 CPPHUB_HAS_AVX2。
# 注: 旗标用普通 add_compile_options — COMPILE_LANGUAGE genex 在 VS 生成器
# 下被静默丢弃 (本仓库 cmake 4.4 实测); CUDA 开启时 nvcc 对该旗标会响亮报错
# 或透传宿主编译器, 两种形态都不产生静默错配。
option(CPPHUB_X86_V3
    "Enable x86-64-v3 (AVX2) codegen - the binary then requires an AVX2-capable CPU"
    OFF)
if(CPPHUB_X86_V3)
    string(TOUPPER "${CMAKE_SYSTEM_PROCESSOR}" _cpphub_arch)
    if(NOT _cpphub_arch MATCHES "^(X86_64|AMD64|X86|I686|I386)$")
        message(FATAL_ERROR
            "CPPHUB_X86_V3=ON requires an x86 host target "
            "(CMAKE_SYSTEM_PROCESSOR='${CMAKE_SYSTEM_PROCESSOR}'); "
            "for cross builds configure ISA flags via toolchain file instead.")
    endif()
    include(CheckCXXCompilerFlag)
    if(MSVC)
        check_cxx_compiler_flag("/arch:AVX2" CPPHUB_COMPILER_HAS_AVX2)
        if(NOT CPPHUB_COMPILER_HAS_AVX2)
            message(FATAL_ERROR "CPPHUB_X86_V3=ON but the compiler rejects /arch:AVX2")
        endif()
        add_compile_options(/arch:AVX2)
    else()
        check_cxx_compiler_flag("-march=x86-64-v3" CPPHUB_COMPILER_HAS_X86_64_V3)
        if(NOT CPPHUB_COMPILER_HAS_X86_64_V3)
            message(FATAL_ERROR "CPPHUB_X86_V3=ON but the compiler rejects -march=x86-64-v3")
        endif()
        add_compile_options(-march=x86-64-v3)
    endif()
    message(STATUS "x86-64-v3 (AVX2) codegen: ENABLED (opt-in; binary requires AVX2-capable CPU)")
endif()

# 头文件库
add_library(cpphub_core INTERFACE)
target_include_directories(cpphub_core INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party  # autodiff header-only library
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/random123/include  # RISK-027: Random123 Philox4x64
)
target_compile_features(cpphub_core INTERFACE cxx_std_20)

# === Phase 6 v1.5: Eigen3 + 计量模块 (ADR-013) ===
# Eigen3 3.4.0 header-only, MPL2.0 (非 MIT), 不编译 Eigen 自身
# 锁定 3.4.0 (C++14 兼容稳定版; 5.0.1 已发布但 v1.5 不升级)
set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen)
add_library(eigen3_interface INTERFACE)
target_include_directories(eigen3_interface INTERFACE ${EIGEN3_INCLUDE_DIR})
# 禁用 Eigen 对 Boost 等的可选依赖 (仅使用 MPL2.0 许可代码)
target_compile_definitions(eigen3_interface INTERFACE EIGEN_MPL2_ONLY)

# 计量模块 INTERFACE 库 (ADR-013: 双层 linalg, 与 cpphub_core 隔离)
# cpphub_core 不链接 eigen3_interface, 定价模块保持 header-only 不受 Eigen3 污染
add_library(cpphub_econometrics INTERFACE)
target_link_libraries(cpphub_econometrics INTERFACE cpphub_core eigen3_interface)
target_include_directories(cpphub_econometrics INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
# MSVC: 抑制 Eigen 内部警告 (C4505 未引用局部函数, C4714 强制 inline 警告)
if(MSVC)
    target_compile_options(cpphub_econometrics INTERFACE /wd4505 /wd4714)
endif()

# === Phase 7C v1.7 M2/M3: 多变量时序 Eigen 载体 (PHASE7C_SPEC §8.2, ADR-013) ===
# var/cointegration 头文件专用; 单变量 timeseries 头文件零 Eigen include
add_library(cpphub_timeseries_mat INTERFACE)
target_link_libraries(cpphub_timeseries_mat INTERFACE cpphub_core eigen3_interface)
target_compile_features(cpphub_timeseries_mat INTERFACE cxx_std_20)
if(MSVC)
    target_compile_options(cpphub_timeseries_mat INTERFACE /wd4505 /wd4714)
endif()

# GoogleTest (FetchContent)
include(FetchContent)
# 供应链 (RISK-033): 以 commit 为不可变身份钉住 googletest v1.14.0
# 不用 tag 名 / tag 归档 URL: tag 可被重指, tag 归档字节不承诺稳定。
# 不用 GIT_REPOSITORY: 本机 git 协议到 github.com:443 不通,
#   而 codeload HTTPS 可达 (URL 形态同时保证本地与 CI 都能全新 configure)。
FetchContent_Declare(
    googletest
    URL https://codeload.github.com/google/googletest/tar.gz/f8d7d77c06936315286eb55f8de22cd23c188571
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
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
