// Phase 4 LITE - M3: GPU MC 引擎 (Host API)
//
// 设计:
//   - 头文件纯 C++ (无 CUDA 依赖), CPU 侧可包含
//   - 实现位于 gpu_mc.cu (CUDA) 或 gpu_mc_cpu_stub.cpp (无 CUDA 回退)
//   - 内核: 每线程 1 条路径, Philox counter RNG, Box-Muller 法线
//   - 精确解一步到期 (欧式 GBM): S_T = S * exp((r-σ²/2)T + σ√T·Z)
//   - 折扣: price = E[exp(-rT) · payoff(S_T)]
//
// CPU/GPU 一致性:
//   - Philox4x64 相同种子 + counter 编号 → 相同 Z (位精确, CPU/GPU 实现等价)
//   - 但求和顺序不同 (GPU warp reduce vs CPU 顺序累加) → 价格相对差 < 1e-12
//   - 单测容差: 价格 1e-6 (MC 误差本身 O(1/√N)), Greeks 1e-4

#pragma once
#include "cpphub/core/types.hpp"
#include "cpphub/performance/gpu/gpu_config.hpp"
#include <vector>
#include <string>

namespace cpphub {
inline namespace v1 {

/// GPU MC 单期权定价结果
struct GpuMCResult {
    Real price = 0.0;        ///< 折现后期权价格
    Real standard_error = 0.0; ///< MC 标准误差 (1σ)
    Real min_payoff = 0.0;   ///< 路径最小 payoff
    Real max_payoff = 0.0;   ///< 路径最大 payoff
    Size n_paths = 0;        ///< 实际模拟路径数
    Real elapsed_ms = 0.0;   ///< 内核执行耗时 (不含 H2D/D2H)
    Real total_ms = 0.0;     ///< 总耗时 (含 H2D/D2H)
    bool used_gpu = false;   ///< 是否实际使用 GPU (false=CPU 回退)
};

/// GPU MC 批量结果 (多期权)
struct GpuMCBatchResult {
    std::vector<Real> prices;        ///< 各期权价格
    std::vector<Real> standard_errors; ///< 各期权标准误差
    Size n_paths_per_option = 0;     ///< 每期权路径数
    Real total_ms = 0.0;             ///< 总耗时
    Real throughput_paths_per_sec = 0.0; ///< 吞吐 (路径/秒)
};

/// GBM 欧式期权参数 (GPU MC 输入)
struct GpuBSMOption {
    Real S = 100.0;      ///< 标的现价
    Real K = 100.0;      ///< 行权价
    Real T = 1.0;        ///< 到期 (年)
    Real r = 0.05;       ///< 无风险利率
    Real q = 0.0;        ///< 股息率
    Real sigma = 0.20;   ///< 波动率
    bool is_call = true; ///< true=看涨, false=看跌
};

/// GPU MC 引擎
///
/// 使用方式:
/// ```cpp
/// GpuMCEngine engine(GpuConfig::Default());
/// GpuBSMOption opt{100, 100, 1.0, 0.05, 0.0, 0.20, true};
/// auto result = engine.price_european_gbm(opt, 1'000'000, /*seed=*/42);
/// ```
class GpuMCEngine {
public:
    /// 构造引擎 (CPU stub 模式下不抛异常, 在 price_* 时抛)
    explicit GpuMCEngine(GpuConfig config = GpuConfig::Default());

    /// 析构 (释放 CUDA 资源, CPU stub 为 no-op)
    ~GpuMCEngine();

    // 不可拷贝 (持有 CUDA stream/event 等资源)
    GpuMCEngine(const GpuMCEngine&) = delete;
    GpuMCEngine& operator=(const GpuMCEngine&) = delete;
    GpuMCEngine(GpuMCEngine&&) noexcept;
    GpuMCEngine& operator=(GpuMCEngine&&) noexcept;

    /// @return 引擎配置
    const GpuConfig& config() const;

    /// @return 设备信息 (CPU stub 返回空 name)
    GpuDeviceInfo device_info() const;

    /// 单期权欧式 GBM MC 定价
    /// @param opt 期权参数
    /// @param n_paths 路径数 (建议 2^k 以对齐 warp)
    /// @param seed RNG 种子 (Philox key)
    /// @throws std::runtime_error 当 CUDA 不可用
    GpuMCResult price_european_gbm(const GpuBSMOption& opt,
                                    Size n_paths,
                                    uint64_t seed = 42) const;

    /// 批量欧式 GBM MC 定价 (多期权并行, 每期权 n_paths 条路径)
    /// @param options 期权数组
    /// @param n_paths_per_option 每期权路径数
    /// @param seed RNG 种子 (每期权 stream_id = i, 确保独立)
    /// @throws std::runtime_error 当 CUDA 不可用
    GpuMCBatchResult price_european_gbm_batch(
        const std::vector<GpuBSMOption>& options,
        Size n_paths_per_option,
        uint64_t seed = 42) const;

    /// Heston MC 定价 (GPU QE 或 Broadie-Kaya, M3 仅 GBM, Heston 留待 v2.1)
    /// @throws std::runtime_error 当 CUDA 不可用 (M3 不实现, 留作占位)
    GpuMCResult price_european_heston(/* HestonOption, ... */) const;

private:
    GpuConfig config_;
    struct Impl;
    Impl* impl_ = nullptr;  ///< PIMPL, 隔离 CUDA 类型
};

}  // namespace cpphub::v1
}  // namespace cpphub
