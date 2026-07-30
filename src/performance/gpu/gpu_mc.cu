// Phase 4 LITE - M3: GPU MC 引擎 (CUDA 实现)
//
// 设计要点:
//   1. Philox4x64-10 RNG 与 CPU (cpphub/core/rng.hpp) 算法完全一致:
//      - 同样常量 (PHILOX_M4x64_0/1, PHILOX_W32_0/1)
//      - 同样 10 轮
//      - 同样 mulhi/mullo 分离
//      → 同 seed + counter → 同 Z (位精确)
//   2. Box-Muller 与 CPU 一致: u=(r>>11)*(1/2^53), r=sqrt(-2 ln u1), theta=2π u2
//   3. 精确解一步到期 (欧式 GBM): S_T = S·exp((r-σ²/2)T + σ√T·Z)
//   4. 每线程 1 条路径, block 级 shared memory 归约 + global atomic
//   5. 双精度强制 (Real=double), RTX 4060 (SM 8.9) 支持原生 FP64
//
// CPU/GPU 数值差异来源:
//   - mulhi: CPU (MSVC _umul128 / GCC __uint128_t) vs GPU (__umul64hi) 实现等价
//   - 求和顺序: GPU warp/block 归约 vs CPU 顺序累加 → 价格相对差 < 1e-12
//   - 单测容差: 价格 1e-6 (MC 误差 O(1/√N) 主导), 标准误差 1e-4

#include "cpphub/performance/gpu/gpu_config.hpp"
#include "cpphub/performance/gpu/gpu_mc.hpp"
#include "cpphub/core/constants.hpp"

#include <cuda_runtime.h>
#include <vector>
#include <stdexcept>
#include <chrono>
#include <cmath>
#include <cstdio>

namespace cpphub {
inline namespace v1 {

// ============================================================================
// 1. 设备端 Philox4x64-10 (与 CPU cpphub/core/rng.hpp 算法一致)
// ============================================================================

namespace {

__device__ __constant__ uint64_t d_PHILOX_M4x64_0 = 0xD2B74407B1CE6E93ULL;
__device__ __constant__ uint64_t d_PHILOX_M4x64_1 = 0x9CCD9A8C5AE13B41ULL;
__device__ __constant__ uint64_t d_PHILOX_W32_0   = 0x9E3779B9;
__device__ __constant__ uint64_t d_PHILOX_W32_1   = 0xBB67AE85;

/// 64x64 高位乘法 (与 CPU mulhi 等价)
__device__ __forceinline__ uint64_t gpu_mulhi(uint64_t a, uint64_t b) {
    return __umul64hi(a, b);
}

/// 64x64 低位乘法
__device__ __forceinline__ uint64_t gpu_mullo(uint64_t a, uint64_t b) {
    return a * b;
}

/// Philox4x64 单次计算 (返回 2 个 64 位输出)
/// key=[seed, stream], counter=[path_id, 0]
__device__ __forceinline__ void philox4x64_10(
    uint64_t c0, uint64_t c1,
    uint64_t k0, uint64_t k1,
    uint64_t& out0, uint64_t& out1)
{
    #pragma unroll
    for (int r = 0; r < 10; ++r) {
        uint64_t hi0 = gpu_mulhi(c0, d_PHILOX_M4x64_0);
        uint64_t lo0 = gpu_mullo(c0, d_PHILOX_M4x64_0);
        uint64_t hi1 = gpu_mulhi(c1, d_PHILOX_M4x64_1);
        uint64_t lo1 = gpu_mullo(c1, d_PHILOX_M4x64_1);
        c0 = hi1 ^ lo0 ^ k0;
        c1 = hi0 ^ lo1 ^ k1;
        k0 += d_PHILOX_W32_0;
        k1 += d_PHILOX_W32_1;
    }
    out0 = c0;
    out1 = c1;
}

/// Box-Muller (与 CPU box_muller 一致)
/// u1, u2 ∈ (0, 1) (由 (r>>11)*(1/2^53) 生成)
__device__ __forceinline__ void gpu_box_muller(double u1, double u2,
                                                double& z1, double& z2) {
    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * 3.14159265358979323846 * u2;
    z1 = r * cos(theta);
    z2 = r * sin(theta);
}

/// 从 64 位随机数生成 (0, 1) double (与 CPU 一致: 53 位精度)
__device__ __forceinline__ double uint64_to_unit(uint64_t r) {
    return (r >> 11) * (1.0 / 9007199254740992.0);  // 1/2^53
}

}  // namespace (anonymous, device-only)

// ============================================================================
// 2. MC 内核: 欧式 GBM 单期权
// ============================================================================

/// 单期权 GBM MC 内核
/// 每线程 1 条路径, block 级归约 + global atomic
__global__ void mc_gbm_european_kernel(
    GpuBSMOption opt,
    Size n_paths,
    uint64_t seed,
    Real* __restrict__ g_sum,     // sum of discounted payoff
    Real* __restrict__ g_sum_sq,  // sum of squared discounted payoff
    Real* __restrict__ g_min,     // min payoff
    Real* __restrict__ g_max)     // max payoff
{
    extern __shared__ Real smem[];

    Real* s_sum    = smem;
    Real* s_sum_sq = smem + blockDim.x;
    Real* s_min    = smem + 2 * blockDim.x;
    Real* s_max    = smem + 3 * blockDim.x;

    Size tid = blockIdx.x * blockDim.x + threadIdx.x;
    Size stride = gridDim.x * blockDim.x;

    Real local_sum = 0.0;
    Real local_sum_sq = 0.0;
    Real local_min = 1e300;
    Real local_max = -1e300;

    // 预计算常量
    Real drift = (opt.r - opt.q - Real(0.5) * opt.sigma * opt.sigma) * opt.T;
    Real vol_sqrt_t = opt.sigma * sqrt(opt.T);
    Real discount = exp(-opt.r * opt.T);

    // 每线程处理多条路径 (grid stride loop)
    for (Size path = tid; path < n_paths; path += stride) {
        // Philox: key=[seed, 0], counter=[path, 0]
        uint64_t r0, r1;
        philox4x64_10(static_cast<uint64_t>(path), 0ULL,
                      seed, 0ULL, r0, r1);

        double u1 = uint64_to_unit(r0);
        double u2 = uint64_to_unit(r1);
        double z1, z2;
        gpu_box_muller(u1, u2, z1, z2);

        // 用 z1 (z2 留作扩展, 例如多因子或 antithetic)
        Real S_T = opt.S * exp(drift + vol_sqrt_t * static_cast<Real>(z1));

        Real payoff;
        if (opt.is_call) {
            payoff = fmax(S_T - opt.K, Real(0));
        } else {
            payoff = fmax(opt.K - S_T, Real(0));
        }

        Real discounted = discount * payoff;

        local_sum    += discounted;
        local_sum_sq += discounted * discounted;
        local_min = fmin(local_min, payoff);
        local_max = fmax(local_max, payoff);
    }

    // Block 级归约 (warp shuffle)
    int lane = threadIdx.x % 32;
    int warp_id = threadIdx.x / 32;

    // Warp 归约
    for (int offset = 16; offset > 0; offset >>= 1) {
        local_sum    += __shfl_down_sync(0xFFFFFFFF, local_sum, offset);
        local_sum_sq += __shfl_down_sync(0xFFFFFFFF, local_sum_sq, offset);
        local_min = fmin(local_min, __shfl_down_sync(0xFFFFFFFF, local_min, offset));
        local_max = fmax(local_max, __shfl_down_sync(0xFFFFFFFF, local_max, offset));
    }

    // 第一个 lane 写入 shared memory
    int n_warps = (blockDim.x + 31) / 32;
    if (lane == 0) {
        s_sum[warp_id]    = local_sum;
        s_sum_sq[warp_id] = local_sum_sq;
        s_min[warp_id]    = local_min;
        s_max[warp_id]    = local_max;
    }
    __syncthreads();

    // 第一个 warp 归约所有 warp 的结果
    if (warp_id == 0) {
        local_sum    = (lane < n_warps) ? s_sum[lane]    : Real(0);
        local_sum_sq = (lane < n_warps) ? s_sum_sq[lane] : Real(0);
        local_min    = (lane < n_warps) ? s_min[lane]    : Real(1e300);
        local_max    = (lane < n_warps) ? s_max[lane]    : Real(-1e300);

        for (int offset = 16; offset > 0; offset >>= 1) {
            local_sum    += __shfl_down_sync(0xFFFFFFFF, local_sum, offset);
            local_sum_sq += __shfl_down_sync(0xFFFFFFFF, local_sum_sq, offset);
            local_min = fmin(local_min, __shfl_down_sync(0xFFFFFFFF, local_min, offset));
            local_max = fmax(local_max, __shfl_down_sync(0xFFFFFFFF, local_max, offset));
        }

        // thread 0 写入 global (atomic)
        if (lane == 0) {
            atomicAdd(g_sum, local_sum);
            atomicAdd(g_sum_sq, local_sum_sq);
            // atomicMin/Max for double via CAS loop
            // 关键: CAS 返回 unsigned long long, 必须用 __longlong_as_double 转回,
            // 不能用隐式转换 (数值转换, 非位重解释), 否则 CAS 失败后无法正确重试
            unsigned long long old_min_ull = __double_as_longlong(*g_min);
            while (local_min < __longlong_as_double(old_min_ull)) {
                unsigned long long assumed = old_min_ull;
                old_min_ull = atomicCAS(reinterpret_cast<unsigned long long*>(g_min),
                                        assumed,
                                        __double_as_longlong(local_min));
                if (old_min_ull == assumed) break;
            }
            unsigned long long old_max_ull = __double_as_longlong(*g_max);
            while (local_max > __longlong_as_double(old_max_ull)) {
                unsigned long long assumed = old_max_ull;
                old_max_ull = atomicCAS(reinterpret_cast<unsigned long long*>(g_max),
                                        assumed,
                                        __double_as_longlong(local_max));
                if (old_max_ull == assumed) break;
            }
        }
    }
}

// ============================================================================
// 3. Host 侧: 设备查询 + GpuMCEngine 实现
// ============================================================================

GpuDeviceInfo gpu_query_device(int device_id) {
    GpuDeviceInfo info;
    info.device_id = device_id;

    cudaDeviceProp prop;
    cudaError_t err = cudaGetDeviceProperties(&prop, device_id);
    if (err != cudaSuccess) {
        return info;  // 空 info
    }

    info.name = prop.name;
    info.compute_capability_major = prop.major;
    info.compute_capability_minor = prop.minor;
    info.total_global_memory = static_cast<Size>(prop.totalGlobalMem);
    info.multiprocessor_count = prop.multiProcessorCount;
    info.max_threads_per_block = prop.maxThreadsPerBlock;
    info.warp_size = prop.warpSize;
    info.supports_double = (prop.major >= 1) &&  // 所有 SM 1.x+ 支持 FP64
                           ((prop.major > 1) || (prop.minor >= 3));
    return info;
}

std::vector<GpuDeviceInfo> gpu_list_devices() {
    std::vector<GpuDeviceInfo> devices;
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    if (err != cudaSuccess) return devices;

    devices.reserve(count);
    for (int i = 0; i < count; ++i) {
        devices.push_back(gpu_query_device(i));
    }
    return devices;
}

// ---------------------------------------------------------------------------
// GpuMCEngine::Impl
// ---------------------------------------------------------------------------

struct GpuMCEngine::Impl {
    GpuConfig config;
    GpuDeviceInfo device;
    cudaStream_t stream = nullptr;

    explicit Impl(GpuConfig c) : config(std::move(c)) {
        // 设置设备
        cudaSetDevice(config.device_id);
        device = gpu_query_device(config.device_id);

        // 创建流
        if (config.stream_count > 1) {
            // 多流模式 (M3 仅单流, 留扩展点)
            cudaStreamCreate(&stream);
        } else {
            cudaStreamCreate(&stream);
        }
    }

    ~Impl() {
        if (stream) {
            cudaStreamDestroy(stream);
        }
        // 重置设备 (同步, 调试用)
        if (config.sync_after_kernel) {
            cudaDeviceSynchronize();
        }
    }
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
    return impl_->device;
}

// ---------------------------------------------------------------------------
// 单期权欧式 GBM MC 定价
// ---------------------------------------------------------------------------

GpuMCResult GpuMCEngine::price_european_gbm(const GpuBSMOption& opt,
                                             Size n_paths,
                                             uint64_t seed) const {
    GpuMCResult result;
    result.n_paths = n_paths;
    result.used_gpu = true;

    auto t_start = std::chrono::steady_clock::now();

    // 分配 device 内存
    Real* d_sum = nullptr;
    Real* d_sum_sq = nullptr;
    Real* d_min = nullptr;
    Real* d_max = nullptr;

    cudaMalloc(&d_sum, sizeof(Real));
    cudaMalloc(&d_sum_sq, sizeof(Real));
    cudaMalloc(&d_min, sizeof(Real));
    cudaMalloc(&d_max, sizeof(Real));

    // 初始化
    Real h_zero = 0.0;
    Real h_max_init = 1e300;
    Real h_min_init = -1e300;
    cudaMemcpy(d_sum, &h_zero, sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_sq, &h_zero, sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d_min, &h_max_init, sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d_max, &h_min_init, sizeof(Real), cudaMemcpyHostToDevice);

    // 启动内核
    int block_size = impl_->config.block_size;
    // 限制 grid 大小避免过多 block (RTX 4060 SM=24, 每 SM 最多 32 block)
    int max_blocks = impl_->device.multiprocessor_count * 8;
    int grid_size = static_cast<int>((n_paths + block_size - 1) / block_size);
    if (grid_size > max_blocks) grid_size = max_blocks;
    if (grid_size < 1) grid_size = 1;

    // Shared memory: 4 arrays × block_size
    size_t smem_bytes = 4 * block_size * sizeof(Real);

    auto t_kernel_start = std::chrono::steady_clock::now();

    mc_gbm_european_kernel<<<grid_size, block_size, smem_bytes, impl_->stream>>>(
        opt, n_paths, seed, d_sum, d_sum_sq, d_min, d_max);

    if (impl_->config.sync_after_kernel) {
        cudaStreamSynchronize(impl_->stream);
    }

    auto t_kernel_end = std::chrono::steady_clock::now();

    // 拷回结果
    Real h_sum, h_sum_sq, h_min, h_max;
    cudaMemcpy(&h_sum, d_sum, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_sum_sq, d_sum_sq, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_min, d_min, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_max, d_max, sizeof(Real), cudaMemcpyDeviceToHost);

    // 计算 price 和 standard error
    result.price = h_sum / static_cast<Real>(n_paths);
    Real mean_sq = h_sum_sq / static_cast<Real>(n_paths);
    Real variance = mean_sq - result.price * result.price;
    if (variance < 0.0) variance = 0.0;  // 数值保护
    result.standard_error = sqrt(variance / static_cast<Real>(n_paths));
    result.min_payoff = h_min;
    result.max_payoff = h_max;

    // 释放
    cudaFree(d_sum);
    cudaFree(d_sum_sq);
    cudaFree(d_min);
    cudaFree(d_max);

    auto t_end = std::chrono::steady_clock::now();

    result.elapsed_ms = std::chrono::duration<Real, std::milli>(t_kernel_end - t_kernel_start).count();
    result.total_ms = std::chrono::duration<Real, std::milli>(t_end - t_start).count();

    return result;
}

// ---------------------------------------------------------------------------
// 批量欧式 GBM MC 定价 (多期权, 串行调用单期权内核)
// ---------------------------------------------------------------------------

GpuMCBatchResult GpuMCEngine::price_european_gbm_batch(
    const std::vector<GpuBSMOption>& options,
    Size n_paths_per_option,
    uint64_t seed) const
{
    GpuMCBatchResult result;
    result.n_paths_per_option = n_paths_per_option;
    result.prices.reserve(options.size());
    result.standard_errors.reserve(options.size());

    auto t_start = std::chrono::steady_clock::now();

    for (Size i = 0; i < options.size(); ++i) {
        // 每期权独立 seed (stream_id = i), 确保路径独立
        auto single = price_european_gbm(options[i], n_paths_per_option,
                                          seed + i);
        result.prices.push_back(single.price);
        result.standard_errors.push_back(single.standard_error);
    }

    auto t_end = std::chrono::steady_clock::now();
    result.total_ms = std::chrono::duration<Real, std::milli>(t_end - t_start).count();

    Size total_paths = options.size() * n_paths_per_option;
    result.throughput_paths_per_sec = static_cast<Real>(total_paths) / (result.total_ms * 1e-3);

    return result;
}

// ---------------------------------------------------------------------------
// Heston (M3 不实现, 留 v2.1)
// ---------------------------------------------------------------------------

GpuMCResult GpuMCEngine::price_european_heston() const {
    throw std::runtime_error(
        "GpuMCEngine::price_european_heston: not implemented in M3 (planned for v2.1).");
}

}  // namespace cpphub::v1
}  // namespace cpphub
