// SOURCE: RISK-001 修复 - 运行时 CPU 特征检测
// 设计原则: 不修改现有编译期宏 (config.hpp/simd.hpp), 仅新增运行时检测能力
// 跨平台: MSVC 用 __cpuid/__cpuidex, GCC/Clang 用 __builtin_cpu_supports
#pragma once
#include "cpphub/core/config.hpp"

#if defined(CPPHUB_COMPILER_MSVC)
    #include <intrin.h>
#elif defined(CPPHUB_COMPILER_GCC) || defined(CPPHUB_COMPILER_CLANG)
    #include <cpuid.h>
#endif

namespace cpphub {
inline namespace v1 {

// 运行时检测的 CPU 特征位
struct CpuFeatures {
    bool has_sse2 = false;
    bool has_avx = false;
    bool has_avx2 = false;
    bool has_avx512f = false;
    bool has_fma = false;
};

// 单次调用后缓存 (Meyers Singleton, C++11+ 线程安全)
inline const CpuFeatures& runtime_cpu_features() {
    static const CpuFeatures features = []() {
        CpuFeatures f;
#if defined(CPPHUB_COMPILER_MSVC)
        int cpuinfo[4] = {0};
        __cpuid(cpuinfo, 0);
        int n_ids = cpuinfo[0];
        if (n_ids >= 1) {
            __cpuidex(cpuinfo, 1, 0);
            f.has_sse2 = (cpuinfo[3] & (1 << 26)) != 0;
            f.has_avx  = (cpuinfo[2] & (1 << 28)) != 0;
            f.has_fma  = (cpuinfo[2] & (1 << 12)) != 0;
        }
        if (n_ids >= 7) {
            __cpuidex(cpuinfo, 7, 0);
            f.has_avx2    = (cpuinfo[1] & (1 << 5))  != 0;
            f.has_avx512f = (cpuinfo[1] & (1 << 16)) != 0;
        }
#elif defined(CPPHUB_COMPILER_GCC) || defined(CPPHUB_COMPILER_CLANG)
        f.has_sse2    = __builtin_cpu_supports("sse2");
        f.has_avx     = __builtin_cpu_supports("avx");
        f.has_avx2    = __builtin_cpu_supports("avx2");
        f.has_avx512f = __builtin_cpu_supports("avx512f");
        f.has_fma     = __builtin_cpu_supports("fma");
#else
        // 未知编译器, 保守禁用所有特性
#endif
        return f;
    }();
    return features;
}

// 便捷查询函数
inline bool runtime_has_sse2()   { return runtime_cpu_features().has_sse2; }
inline bool runtime_has_avx()    { return runtime_cpu_features().has_avx; }
inline bool runtime_has_avx2()   { return runtime_cpu_features().has_avx2; }
inline bool runtime_has_avx512() { return runtime_cpu_features().has_avx512f; }
inline bool runtime_has_fma()    { return runtime_cpu_features().has_fma; }

}  // namespace v1
}  // namespace cpphub
