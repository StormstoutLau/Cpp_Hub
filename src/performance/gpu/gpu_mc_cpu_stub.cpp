// Phase 4 LITE - M3: GPU MC CPU Stub (无 CUDA 时的回退实现)
//
// 当 CPPHUB_ENABLE_CUDA=OFF 或 nvcc 不可用时编译此文件。
// 所有 GPU 调用抛出 std::runtime_error, 但库仍可链接/编译,
// 保证 A/B 站 (无 CUDA) 的 C++ 核心测试不受影响。

#include "cpphub/performance/gpu/gpu_config.hpp"
#include "cpphub/performance/gpu/gpu_mc.hpp"
#include <stdexcept>
#include <vector>

namespace cpphub {
inline namespace v1 {

// ============================================================================
// 设备查询 (CPU stub: 返回空)
// ============================================================================

GpuDeviceInfo gpu_query_device(int /*device_id*/) {
    // CPU stub: 返回空设备信息, name 为空表示无 GPU
    return GpuDeviceInfo{};
}

std::vector<GpuDeviceInfo> gpu_list_devices() {
    // CPU stub: 无可见设备
    return {};
}

// ============================================================================
// GpuMCEngine (CPU stub: 所有调用抛异常)
// ============================================================================

struct GpuMCEngine::Impl {
    GpuConfig config;
    explicit Impl(GpuConfig c) : config(std::move(c)) {}
};

GpuMCEngine::GpuMCEngine(GpuConfig config)
    : config_(std::move(config)), impl_(new Impl(config_)) {}

GpuMCEngine::~GpuMCEngine() = default;

GpuMCEngine::GpuMCEngine(GpuMCEngine&&) noexcept = default;
GpuMCEngine& GpuMCEngine::operator=(GpuMCEngine&&) noexcept = default;

const GpuConfig& GpuMCEngine::config() const {
    return impl_->config;
}

GpuDeviceInfo GpuMCEngine::device_info() const {
    // CPU stub: 返回空
    return GpuDeviceInfo{};
}

GpuMCResult GpuMCEngine::price_european_gbm(const GpuBSMOption& /*opt*/,
                                             Size /*n_paths*/,
                                             uint64_t /*seed*/) const {
    throw std::runtime_error(
        "GpuMCEngine::price_european_gbm: CUDA not available. "
        "Rebuild with -DCPPHUB_ENABLE_CUDA=ON and NVIDIA GPU + nvcc installed.");
}

GpuMCBatchResult GpuMCEngine::price_european_gbm_batch(
    const std::vector<GpuBSMOption>& /*options*/,
    Size /*n_paths_per_option*/,
    uint64_t /*seed*/) const {
    throw std::runtime_error(
        "GpuMCEngine::price_european_gbm_batch: CUDA not available. "
        "Rebuild with -DCPPHUB_ENABLE_CUDA=ON and NVIDIA GPU + nvcc installed.");
}

GpuMCResult GpuMCEngine::price_european_heston() const {
    throw std::runtime_error(
        "GpuMCEngine::price_european_heston: not implemented in M3 (planned for v2.1).");
}

}  // namespace cpphub::v1
}  // namespace cpphub
