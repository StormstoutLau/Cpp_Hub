// RISK-002 验收: Philox 性能基准测试 (信息性, 无 pass/fail 断言)
// 对比 MSVC (_umul128) vs GCC (__uint128_t) 路径的吞吐量
// 不注册到 ctest, 作为独立可执行文件手动运行
#include "cpphub/core/rng.hpp"
#include "cpphub/core/config.hpp"
#include <chrono>
#include <iostream>
#include <vector>
#include <cstdint>

int main() {
    constexpr size_t N = 100'000'000;  // 10^8

    std::cout << "=== Philox4x64 Benchmark ===" << std::endl;
    std::cout << "平台: ";
#if defined(CPPHUB_COMPILER_MSVC)
    std::cout << "MSVC (_umul128 path)";
#elif defined(CPPHUB_COMPILER_GCC)
    std::cout << "GCC (__uint128_t path)";
#elif defined(CPPHUB_COMPILER_CLANG)
    std::cout << "Clang (__uint128_t path)";
#else
    std::cout << "Unknown";
#endif
    std::cout << std::endl;
    std::cout << "生成数量: " << N << std::endl;

    cpphub::Philox4x64 rng(42);

    // 预分配避免分配器干扰
    volatile uint64_t sink = 0;

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N; ++i) {
        sink ^= rng();
    }
    auto end = std::chrono::high_resolution_clock::now();

    auto duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    double duration_sec = static_cast<double>(duration_ns) / 1e9;
    double throughput = static_cast<double>(N) / duration_sec;

    std::cout << "总耗时: " << duration_sec << " sec" << std::endl;
    std::cout << "吞吐量: " << throughput << " numbers/sec" << std::endl;
    std::cout << "单次耗时: " << (duration_sec / N * 1e9) << " ns/number" << std::endl;
    std::cout << "(sink=" << sink << ")" << std::endl;

    return 0;
}
