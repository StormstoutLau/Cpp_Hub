// SOURCE: PHASE1_SPEC §2.4 - 编译配置、SIMD检测、版本
#pragma once

// 版本（单一源: 仓库根 VERSION → tools/sync_version.py 生成, RISK-032）
#include "cpphub/core/version.hpp"

// 编译器检测
#if defined(_MSC_VER)
    #define CPPHUB_COMPILER_MSVC 1
#elif defined(__GNUC__)
    #define CPPHUB_COMPILER_GCC 1
#elif defined(__clang__)
    #define CPPHUB_COMPILER_CLANG 1
#endif

// SIMD 检测 (编译期)
#if defined(__AVX2__) || (defined(_MSC_VER) && defined(__AVX2__))
    #define CPPHUB_HAS_AVX2 1
#endif
#if defined(__AVX512F__)
    #define CPPHUB_HAS_AVX512 1
#endif
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    #define CPPHUB_HAS_NEON 1
#endif

// 强制内联
#if defined(CPPHUB_COMPILER_MSVC)
    #define CPPHUB_FORCE_INLINE __forceinline
#else
    #define CPPHUB_FORCE_INLINE inline __attribute__((always_inline))
#endif

// 命名空间
namespace cpphub {
inline namespace v1 {
}  // namespace v1
}  // namespace cpphub
