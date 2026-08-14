// RISK-001 验收测试: 运行时 CPU 特征检测
#include "cpphub/core/cpu_features.hpp"
#include "cpphub/core/config.hpp"
#include <gtest/gtest.h>
#include <thread>
#include <vector>

using namespace cpphub;

// 单向蕴含: 若编译期宏定义, 则运行时必须检测到
// (反向不成立: 运行时检测到但编译期未启用是合法的)
TEST(CpuFeatures, RuntimeConsistentWithCompileTime) {
    const CpuFeatures& f = runtime_cpu_features();
#ifdef CPPHUB_HAS_AVX2
    EXPECT_TRUE(f.has_avx2) << "编译期 CPPHUB_HAS_AVX2 已定义, 运行时必须检测到 AVX2";
#endif
#ifdef CPPHUB_HAS_AVX512
    EXPECT_TRUE(f.has_avx512f) << "编译期 CPPHUB_HAS_AVX512 已定义, 运行时必须检测到 AVX-512F";
#endif
    // SSE2 在 x86-64 上必定存在
#if defined(_M_X64) || defined(__x86_64__)
    EXPECT_TRUE(f.has_sse2) << "x86-64 必定支持 SSE2";
#endif
}

// Meyers Singleton 缓存: 两次调用返回同一引用
TEST(CpuFeatures, CachedAfterFirstCall) {
    const CpuFeatures& f1 = runtime_cpu_features();
    const CpuFeatures& f2 = runtime_cpu_features();
    EXPECT_EQ(&f1, &f2) << "runtime_cpu_features() 应返回同一缓存对象的引用";
}

// 线程安全: 多线程并发调用不崩溃
TEST(CpuFeatures, ThreadSafeConcurrentAccess) {
    constexpr int n_threads = 8;
    std::vector<std::thread> threads;
    std::vector<bool> results(n_threads, false);
    for (int t = 0; t < n_threads; ++t) {
        threads.emplace_back([&results, t]() {
            const CpuFeatures& f = runtime_cpu_features();
            results[t] = f.has_sse2;  // 读取字段验证不崩溃
        });
    }
    for (auto& th : threads) th.join();
    // 所有线程应成功完成并读取到值
    for (int t = 0; t < n_threads; ++t) {
        EXPECT_TRUE(results[t]) << "线程 " << t << " 未能正确读取 CPU 特征";
    }
}
