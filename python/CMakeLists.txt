# Phase 4 LITE - M2: Python 绑定构建脚本 (nanobind)
#
# 此 CMakeLists.txt 是自包含顶层文件,scikit-build-core 直接构建此目录。
# 它独立于根 CMakeLists.txt (后者用于 C++ 测试/基准)。
#
# 依赖: nanobind (pip install nanobind), Python 开发头文件
# 构建目标: cpphub._core (CPython 扩展模块)
# Python 包: src/cpphub/__init__.py (re-export)

cmake_minimum_required(VERSION 3.25)

# RISK-032: Python 发行版本单一源 = pyproject.toml。
# wheel 构建时 scikit-build-core 注入 SKBUILD_PROJECT_VERSION;
# 直接构建 python/ 目录时回退到仓库根 VERSION (库版本), 保证无硬编码。
if(DEFINED SKBUILD_PROJECT_VERSION)
    set(CPPHUB_PY_VERSION "${SKBUILD_PROJECT_VERSION}")
else()
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../VERSION" CPPHUB_PY_VERSION LIMIT_COUNT 1)
endif()
project(cpphub_python VERSION ${CPPHUB_PY_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# MSVC UTF-8 (Python 3.11+ 模块字符串默认 UTF-8)
if(MSVC)
    add_compile_options(/utf-8 /W3)
    add_compile_definitions(_USE_MATH_DEFINES)
endif()

# ----------------------------------------------------------------------------
# 1. Cpp_Hub 头文件接口库 (header-only, 仅包含路径)
# ----------------------------------------------------------------------------
add_library(cpphub_headers INTERFACE)
target_include_directories(cpphub_headers INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../third_party  # autodiff header-only
    ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/eigen  # Eigen3 3.4.0 (MPL2)
    ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/random123/include  # RISK-027
)
target_compile_features(cpphub_headers INTERFACE cxx_std_20)
# 与根 CMakeLists.txt 的 eigen3_interface 保持一致: 仅启用 MPL2.0 许可部分
target_compile_definitions(cpphub_headers INTERFACE EIGEN_MPL2_ONLY)

# ----------------------------------------------------------------------------
# 2. 寻找 Python + nanobind
#    (nanobind 要求先 find_package(Python ...))
# ----------------------------------------------------------------------------
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.cmake_dir())"
    OUTPUT_VARIABLE nanobind_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE nanobind_result
)

if(nanobind_result AND NOT nanobind_result EQUAL 0)
    message(FATAL_ERROR
        "nanobind not found. Install with: pip install nanobind\n"
        "Or disable Python bindings by not invoking this CMakeLists.txt.")
endif()

list(APPEND CMAKE_PREFIX_PATH "${nanobind_DIR}")
find_package(nanobind CONFIG REQUIRED)

# ----------------------------------------------------------------------------
# 3. 创建 cpphub._core 扩展模块
#    (NB_MODULE 名为 _core, 与 __init__.py 的 from ._core import 对齐)
# ----------------------------------------------------------------------------
nanobind_add_module(_core
    src/cpphub_python.cpp
    # MCP 工具面 (spec/mcp-server): JSON-through-C-ABI 实现作为独立 TU 编入,
    # 避免与 nanobind TU 的 calibrator.hpp 宏重命名冲突; 无 ctypes/无独立 .so。
    ../src/c_api/c_api.cpp
)

target_link_libraries(_core PRIVATE cpphub_headers)

# RISK-032: __version__ 由构建系统注入 (cpphub_python.cpp 消费 CPPHUB_PY_VERSION 宏)
target_compile_definitions(_core PRIVATE CPPHUB_PY_VERSION="${CPPHUB_PY_VERSION}")

# MSVC UTF-8 (扩展模块)
if(MSVC)
    target_compile_options(_core PRIVATE /utf-8 /W3)
    target_compile_definitions(_core PRIVATE _USE_MATH_DEFINES)
endif()

# ----------------------------------------------------------------------------
# 4. 安装到 Python 包目录
#    scikit-build-core 将 install prefix 中的文件按路径并入 wheel,
#    故 DESTINATION 应为 "cpphub" (wheel 中的包名),与 __init__.py 合并。
# ----------------------------------------------------------------------------
install(TARGETS _core
    LIBRARY DESTINATION cpphub
    RUNTIME DESTINATION cpphub
)
