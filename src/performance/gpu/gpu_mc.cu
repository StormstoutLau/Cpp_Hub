// Phase 4 LITE - M3: GPU MC 引擎 (CUDA 实现)
//
// 设计要点:
//   1. 同一 Random123 Philox4x64-10 实现 (third_party/random123), 常量/轮数由上游接管:
//      - 设备端与 CPU (cpphub/core/rng.hpp) 共用同一 vendor 头
//      → 同 seed + counter → 同 Z (uniform 位精确)
//   2. Box-Muller 与 CPU 一致: u=(r>>11)*(1/2^53), r=sqrt(-2 ln u1), theta=2π u2
//   3. 精确解一步到期 (欧式 GBM): S_T = S·exp((r-σ²/2)T + σ√T·Z)
//   4. 每线程 1 条路径, block 级 shared memory 归约 + global atomic
//   5. 双精度强制 (Real=double), RTX 4060 (SM 8.9) 支持原生 FP64
//
// CPU/GPU 数值差异来源:
//   - mulhi: 上游 Random123 按平台分派 (MSVC _umul128 / GCC __uint128_t / CUDA __umul64hi) 等价
//   - 求和顺序: GPU warp/block 归约 vs CPU 顺序累加 → 价格相对差 < 1e-12
//   - 单测容差: 价格 1e-6 (MC 误差 O(1/√N) 主导), 标准误差 1e-4

#include "cpphub/performance/gpu/gpu_config.hpp"
#include "cpphub/performance/gpu/gpu_mc.hpp"
#include "cpphub/core/constants.hpp"

// RISK-027 / DESIGN §11 D-6: Windows+nvcc 下 msvcfeatures.h 与 nvccfeatures.h
// 同时置位 R123_USE_MULHILO64_MSVC_INTRIN 与 _CUDA_INTRIN, 而 philox.h 让前者优先
// → device 代码调用 host-only 的 _umul128 编译失败。必须在 include 前关闭。
#if defined(__CUDA_ARCH__)
#define R123_USE_MULHILO64_MSVC_INTRIN 0
#endif
#include <Random123/philox.h>
#include "cpphub/performance/gpu/rng_device_probe.hpp"

#include <cuda_runtime.h>
#include <atomic>
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>
#include <chrono>
#include <cmath>
#include <cstdio>

namespace cpphub {
inline namespace v1 {

// ============================================================================
// 1. 设备端工具 (Box-Muller / 单位化; Philox 由 Random123 提供)
// ============================================================================

namespace {

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
        // Philox: Random123 vendor (RISK-027), key=[seed, 0], counter=[path, 0, 0, 0]
        r123::Philox4x64::ctr_type ctr{};
        ctr.v[0] = static_cast<uint64_t>(path);
        r123::Philox4x64::key_type key{};
        key.v[0] = seed;
        const r123::Philox4x64::ctr_type r = r123::Philox4x64{}(ctr, key);

        double u1 = uint64_to_unit(r.v[0]);
        double u2 = uint64_to_unit(r.v[1]);
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
    // 全掩码 0xFFFFFFFF 的安全性: Impl 构造期已校验 block_size 为 32 的正整数倍
    // (RISK-026), 归约时所有 lane 必然活跃, 不存在部分掩码未定义行为。
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
// CUDA 错误检查 (RISK-026): 所有 CUDA 运行时调用必须显式检查
// ---------------------------------------------------------------------------

namespace {

/// 检查 CUDA 返回码; 失败抛 std::runtime_error (携带调用点与 cudaGetErrorString)
void cuda_check(cudaError_t err, const char* call, const char* file, int line) {
    if (err == cudaSuccess) return;
    throw std::runtime_error(
        std::string("GpuMCEngine: CUDA 调用失败: ") + cudaGetErrorString(err) +
        " (" + std::to_string(static_cast<int>(err)) + ") in " + call +
        " at " + file + ":" + std::to_string(line));
}

/// 设备显存 RAII (RISK-025): 作用域退出即 cudaFree, 覆盖异常与提前返回路径
class DeviceBuffer {
public:
    explicit DeviceBuffer(std::size_t bytes) {
        void* p = nullptr;
        cuda_check(cudaMalloc(&p, bytes), "cudaMalloc", __FILE__, __LINE__);
        ptr_ = static_cast<Real*>(p);
    }
    ~DeviceBuffer() {
        if (ptr_) cudaFree(ptr_);  // 析构不可抛异常, 忽略返回值
    }
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    DeviceBuffer(DeviceBuffer&& other) noexcept : ptr_(other.ptr_) {
        other.ptr_ = nullptr;
    }
    DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
        if (this != &other) {
            if (ptr_) cudaFree(ptr_);
            ptr_ = other.ptr_;
            other.ptr_ = nullptr;
        }
        return *this;
    }

    Real* get() const { return ptr_; }

private:
    Real* ptr_ = nullptr;
};

}  // namespace

/// 包装 CUDA 调用并检查返回码 (RISK-026)
#define CPPHUB_CUDA_CHECK(call) \
    ::cpphub::v1::cuda_check((call), #call, __FILE__, __LINE__)

// ---------------------------------------------------------------------------
// GpuMCEngine::Impl
// ---------------------------------------------------------------------------

namespace {
/// RISK-025 诊断钩子: 存活的 Impl 计数 (验证析构释放; 仅测试用)
std::atomic<std::size_t> g_live_impl_count{0};
}  // namespace

struct GpuMCEngine::Impl {
    GpuConfig config;
    GpuDeviceInfo device;
    cudaStream_t stream = nullptr;

    explicit Impl(GpuConfig c) : config(std::move(c)) {
        // RISK-026: 配置校验前置 —— 非法配置不创建任何资源
        // block_size 必须为 32 的正整数倍: warp 归约使用全掩码 __shfl_down_sync,
        // 非满 warp 的 lane 不活跃会触发未定义行为
        if (config.block_size <= 0 || config.block_size % 32 != 0) {
            throw std::invalid_argument(
                "GpuMCEngine: block_size 必须为 32 的正整数倍, 实际 " +
                std::to_string(config.block_size));
        }

        int device_count = 0;
        CPPHUB_CUDA_CHECK(cudaGetDeviceCount(&device_count));
        if (config.device_id < 0 || config.device_id >= device_count) {
            throw std::runtime_error(
                "GpuMCEngine: 非法 device_id " + std::to_string(config.device_id) +
                " (可见设备数: " + std::to_string(device_count) + ")");
        }

        device = gpu_query_device(config.device_id);
        // 属性查询失败时 gpu_query_device 返回空 info (max_threads_per_block=0), 需显式报错
        if (device.max_threads_per_block <= 0) {
            throw std::runtime_error(
                "GpuMCEngine: 无法查询设备属性 (device_id " +
                std::to_string(config.device_id) + ")");
        }
        if (config.block_size > device.max_threads_per_block) {
            throw std::invalid_argument(
                "GpuMCEngine: block_size " + std::to_string(config.block_size) +
                " 超过设备上限 max_threads_per_block " +
                std::to_string(device.max_threads_per_block));
        }

        CPPHUB_CUDA_CHECK(cudaSetDevice(config.device_id));
        // 多流模式 (stream_count > 1) 留待扩展; M3 仅单流
        CPPHUB_CUDA_CHECK(cudaStreamCreate(&stream));

        ++g_live_impl_count;  // 全部资源就绪后才计数: 构造抛异常不泄漏计数
    }

    ~Impl() {
        if (stream) {
            // 析构不可抛异常; cudaStreamDestroy 保证流上残留工作完成后自动回收资源
            cudaStreamDestroy(stream);
            stream = nullptr;
        }
        --g_live_impl_count;
    }
};

Size GpuMCEngine::live_impl_count() {
    return g_live_impl_count.load(std::memory_order_relaxed);
}

GpuMCEngine::GpuMCEngine(GpuConfig config)
    : impl_(std::make_unique<Impl>(std::move(config))) {}

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

    auto t_start = std::chrono::steady_clock::now();

    // 分配 device 内存 (RAII: 异常/提前返回路径亦释放)
    DeviceBuffer d_sum(sizeof(Real));
    DeviceBuffer d_sum_sq(sizeof(Real));
    DeviceBuffer d_min(sizeof(Real));
    DeviceBuffer d_max(sizeof(Real));

    // 初始化
    Real h_zero = 0.0;
    Real h_max_init = 1e300;
    Real h_min_init = -1e300;
    CPPHUB_CUDA_CHECK(cudaMemcpy(d_sum.get(), &h_zero, sizeof(Real), cudaMemcpyHostToDevice));
    CPPHUB_CUDA_CHECK(cudaMemcpy(d_sum_sq.get(), &h_zero, sizeof(Real), cudaMemcpyHostToDevice));
    CPPHUB_CUDA_CHECK(cudaMemcpy(d_min.get(), &h_max_init, sizeof(Real), cudaMemcpyHostToDevice));
    CPPHUB_CUDA_CHECK(cudaMemcpy(d_max.get(), &h_min_init, sizeof(Real), cudaMemcpyHostToDevice));

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
        opt, n_paths, seed, d_sum.get(), d_sum_sq.get(), d_min.get(), d_max.get());
    // 启动失败 (显存不足 / 配置非法) 必须显式报错, 不得静默返回 price=0
    CPPHUB_CUDA_CHECK(cudaGetLastError());

    if (impl_->config.sync_after_kernel) {
        // 捕获运行期错误 (如非法地址访问), 而非静默返回伪结果
        CPPHUB_CUDA_CHECK(cudaStreamSynchronize(impl_->stream));
    }

    auto t_kernel_end = std::chrono::steady_clock::now();

    // 拷回结果
    Real h_sum = 0.0, h_sum_sq = 0.0, h_min = 0.0, h_max = 0.0;
    CPPHUB_CUDA_CHECK(cudaMemcpy(&h_sum, d_sum.get(), sizeof(Real), cudaMemcpyDeviceToHost));
    CPPHUB_CUDA_CHECK(cudaMemcpy(&h_sum_sq, d_sum_sq.get(), sizeof(Real), cudaMemcpyDeviceToHost));
    CPPHUB_CUDA_CHECK(cudaMemcpy(&h_min, d_min.get(), sizeof(Real), cudaMemcpyDeviceToHost));
    CPPHUB_CUDA_CHECK(cudaMemcpy(&h_max, d_max.get(), sizeof(Real), cudaMemcpyDeviceToHost));

    // 计算 price 和 standard error
    result.price = h_sum / static_cast<Real>(n_paths);
    Real mean_sq = h_sum_sq / static_cast<Real>(n_paths);
    Real variance = mean_sq - result.price * result.price;
    if (variance < 0.0) variance = 0.0;  // 数值保护
    result.standard_error = sqrt(variance / static_cast<Real>(n_paths));
    result.min_payoff = h_min;
    result.max_payoff = h_max;

    // 全部 CUDA 调用成功后才标记 (RISK-026: 此前提前置 true, 失败时伪报成功)
    result.used_gpu = true;

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

// ============================================================================
// 4. 测试探针 (RISK-027): 设备端 Philox/KAT/Box-Muller 结果导出到 host
// ============================================================================

namespace {

constexpr int kProbeBlock = 128;

__global__ void probe_philox_block_kernel(uint64_t seed, uint64_t* out, uint64_t n) {
    const uint64_t path = blockIdx.x * blockDim.x + threadIdx.x;
    if (path >= n) return;  // grid 向上取整 → 必须守卫, 否则越界写
    r123::Philox4x64::ctr_type ctr{};
    ctr.v[0] = path;
    r123::Philox4x64::key_type key{};
    key.v[0] = seed;
    const r123::Philox4x64::ctr_type r = r123::Philox4x64{}(ctr, key);
    for (int i = 0; i < 4; ++i) out[path * 4 + i] = r.v[i];
}

__global__ void probe_philox_kat_kernel(const gpu_probe::KatVector* vecs,
                                        uint64_t* out, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    r123::Philox4x64::ctr_type ctr{};
    for (int j = 0; j < 4; ++j) ctr.v[j] = vecs[i].ctr[j];
    r123::Philox4x64::key_type key{};
    key.v[0] = vecs[i].key[0];
    key.v[1] = vecs[i].key[1];
    const r123::Philox4x64::ctr_type r = r123::Philox4x64{}(ctr, key);
    for (int j = 0; j < 4; ++j) out[i * 4 + j] = r.v[j];
}

__global__ void probe_box_muller_z1_kernel(uint64_t seed, double* out, uint64_t n) {
    const uint64_t path = blockIdx.x * blockDim.x + threadIdx.x;
    if (path >= n) return;
    r123::Philox4x64::ctr_type ctr{};
    ctr.v[0] = path;
    r123::Philox4x64::key_type key{};
    key.v[0] = seed;
    const r123::Philox4x64::ctr_type r = r123::Philox4x64{}(ctr, key);
    double z1, z2;
    gpu_box_muller(uint64_to_unit(r.v[0]), uint64_to_unit(r.v[1]), z1, z2);
    out[path] = z1;
}

unsigned probe_grid(std::size_t n) {
    return static_cast<unsigned>((n + kProbeBlock - 1) / kProbeBlock);
}

}  // namespace

namespace gpu_probe {  // 已在 cpphub::v1 内, 直接开子命名空间

std::vector<std::uint64_t> philox_blocks(std::uint64_t seed, std::size_t n_blocks) {
    if (n_blocks == 0) return {};
    std::vector<std::uint64_t> host(n_blocks * 4);
    DeviceBuffer dev(host.size() * sizeof(std::uint64_t));
    probe_philox_block_kernel<<<probe_grid(n_blocks), kProbeBlock>>>(
        seed, reinterpret_cast<std::uint64_t*>(dev.get()),
        static_cast<std::uint64_t>(n_blocks));
    CPPHUB_CUDA_CHECK(cudaGetLastError());
    CPPHUB_CUDA_CHECK(cudaMemcpy(host.data(), dev.get(),
                                 host.size() * sizeof(std::uint64_t),
                                 cudaMemcpyDeviceToHost));
    return host;
}

std::vector<std::uint64_t> philox_kat(const std::vector<KatVector>& vectors) {
    if (vectors.empty()) return {};
    std::vector<std::uint64_t> out(vectors.size() * 4);
    DeviceBuffer d_in(vectors.size() * sizeof(KatVector));
    DeviceBuffer d_out(out.size() * sizeof(std::uint64_t));
    CPPHUB_CUDA_CHECK(cudaMemcpy(d_in.get(), vectors.data(),
                                 vectors.size() * sizeof(KatVector),
                                 cudaMemcpyHostToDevice));
    probe_philox_kat_kernel<<<probe_grid(vectors.size()), kProbeBlock>>>(
        reinterpret_cast<const KatVector*>(d_in.get()),
        reinterpret_cast<std::uint64_t*>(d_out.get()),
        static_cast<int>(vectors.size()));
    CPPHUB_CUDA_CHECK(cudaGetLastError());
    CPPHUB_CUDA_CHECK(cudaMemcpy(out.data(), d_out.get(),
                                 out.size() * sizeof(std::uint64_t),
                                 cudaMemcpyDeviceToHost));
    return out;
}

std::vector<double> box_muller_z1_batch(std::uint64_t seed, std::size_t n) {
    if (n == 0) return {};
    std::vector<double> out(n);
    DeviceBuffer dev(n * sizeof(double));
    probe_box_muller_z1_kernel<<<probe_grid(n), kProbeBlock>>>(
        seed, dev.get(), static_cast<std::uint64_t>(n));
    CPPHUB_CUDA_CHECK(cudaGetLastError());
    CPPHUB_CUDA_CHECK(cudaMemcpy(out.data(), dev.get(), n * sizeof(double),
                                 cudaMemcpyDeviceToHost));
    return out;
}

}  // namespace gpu_probe

}  // namespace cpphub::v1
}  // namespace cpphub
