// Phase 4 LITE - M3: GPU 配置实现 (纯 C++, 不依赖 CUDA)
//
// gpu_query_device / gpu_list_devices 的实现位于:
//   - gpu_mc.cu (CUDA 可用时, 调用 cudaGetDeviceProperties)
//   - gpu_mc_cpu_stub.cpp (无 CUDA 时, 返回空/抛异常)

#include "cpphub/performance/gpu/gpu_config.hpp"
#include <sstream>
#include <iomanip>

namespace cpphub {
inline namespace v1 {

std::string GpuDeviceInfo::describe() const {
    std::ostringstream oss;
    oss << "GPU[" << device_id << "] " << name
        << " (SM " << compute_capability_major << "." << compute_capability_minor << ")"
        << " VRAM=" << std::fixed << std::setprecision(1)
        << static_cast<double>(total_global_memory) / (1024.0 * 1024.0 * 1024.0) << " GB"
        << " SMs=" << multiprocessor_count
        << " maxThreads/block=" << max_threads_per_block
        << " warp=" << warp_size
        << " double=" << (supports_double ? "yes" : "no");
    return oss.str();
}

}  // namespace cpphub::v1
}  // namespace cpphub
