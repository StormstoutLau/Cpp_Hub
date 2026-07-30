// Phase 4 LITE - M3: GPU MC 配置
//
// GPU MC 引擎配置: 设备选择、block size、流数量、内存策略。
// 主控站 RTX 4060 (Ada Lovelace, SM 8.9, 8GB VRAM)。
//
// 设计原则:
//   - 头文件纯 C++ (不依赖 CUDA),便于 CPU 侧测试与跨平台编译
//   - 实现细节 (.cu) 隔离在 src/performance/gpu/
//   - 无 CUDA 时回退到 CPU stub (抛出 std::runtime_error)

#pragma once
#include "cpphub/core/types.hpp"
#include <string>
#include <vector>

namespace cpphub {
inline namespace v1 {

/// GPU MC 引擎配置
struct GpuConfig {
    int device_id = 0;          ///< CUDA 设备 ID (默认 0)
    int block_size = 256;       ///< 每块线程数 (warp 32 的倍数, 256 为 Ada 最优)
    int stream_count = 1;       ///< CUDA 流数量 (1=单流, 多流用于异步重叠)
    Size memory_pool_bytes = 0; ///< 显存池 (0 = 按需分配, >0 = 预分配)
    bool sync_after_kernel = true;  ///< 内核后是否 cudaDeviceSynchronize
    bool use_managed_memory = false; ///< 是否使用统一内存 (调试用, 生产用 device malloc)

    /// 默认配置 (RTX 4060 优化: 256 线程/block, 单流, 按需分配)
    static GpuConfig Default() { return {}; }

    /// 大批量配置 (10000+ 期权, 预分配显存池)
    static GpuConfig Batch() {
        GpuConfig cfg;
        cfg.memory_pool_bytes = 512ULL * 1024 * 1024;  // 512 MB
        return cfg;
    }
};

/// GPU 设备信息 (运行时查询)
struct GpuDeviceInfo {
    std::string name;             ///< 设备名 (e.g. "NVIDIA GeForce RTX 4060")
    int device_id = 0;            ///< 设备 ID
    int compute_capability_major = 0;  ///< SM 主版本 (Ada=8)
    int compute_capability_minor = 0;  ///< SM 次版本 (Ada=9 → SM 8.9)
    Size total_global_memory = 0; ///< 全局显存 (字节)
    int multiprocessor_count = 0; ///< SM 数量
    int max_threads_per_block = 0; ///< 每块最大线程数
    int warp_size = 32;           ///< warp 大小 (NVIDIA 恒为 32)
    bool supports_double = true;  ///< 是否支持双精度 (Ada 支持)

    /// 格式化字符串 (调试用)
    std::string describe() const;
};

/// 查询指定设备 ID 的信息 (CPU stub 返回空)
/// @throws std::runtime_error 当 CUDA 不可用 (CPU stub 模式)
GpuDeviceInfo gpu_query_device(int device_id);

/// 列出所有可见 GPU 设备 (CPU stub 返回空 vector)
std::vector<GpuDeviceInfo> gpu_list_devices();

}  // namespace cpphub::v1
}  // namespace cpphub
