// Phase 4 LITE - M3-3: GPU MC 完整 15 测试 (主控站 RTX 4060 验证)
//
// 测试矩阵:
//   1-3:  GpuConfig / GpuDeviceInfo 基础结构
//   4-6:  设备查询 + 引擎构造
//   7-10: MC 定价 vs BSM 解析解 (ATM/ITM/OTM call + ATM put)
//   11-12: 确定性 (同 seed 同结果, 不同 seed 接近)
//   13:   标准误差收敛 (4x 路径 → 2x SE)
//   14:   批量定价
//   15:   性能 (1M 路径耗时)
//   16-19: RISK-025/026 回归 (PIMPL 释放验证 + CUDA 错误检查 + 配置校验)
//
// 容差策略:
//   - MC 价格: 4*SE (99.99% 置信区间, 避免偶发失败)
//   - 确定性: 位精确 (同 seed)
//   - SE 收敛: 1.5x ~ 2.5x (理论 2x, 留数值余量)

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "cpphub/performance/gpu/gpu_config.hpp"
#include "cpphub/performance/gpu/gpu_mc.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/performance/gpu/rng_device_probe.hpp"
#include "kat_philox4x64_10.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"

using namespace cpphub;

// ============================================================================
// 辅助: BSM 解析价格
// ============================================================================
static double bsm_price(double S, double K, double T, double r, double q,
                         double sigma, bool is_call) {
    return AnalyticGreeksEngine::bsm_european(S, K, T, r, q, sigma, is_call).price;
}

// ============================================================================
// 1-3: GpuConfig / GpuDeviceInfo 基础结构
// ============================================================================

TEST(GpuMC, GpuConfigDefault) {
    GpuConfig cfg = GpuConfig::Default();
    EXPECT_EQ(cfg.device_id, 0);
    EXPECT_EQ(cfg.block_size, 256);
    EXPECT_EQ(cfg.stream_count, 1);
    EXPECT_EQ(cfg.memory_pool_bytes, 0u);
    EXPECT_TRUE(cfg.sync_after_kernel);
    EXPECT_FALSE(cfg.use_managed_memory);
}

TEST(GpuMC, GpuConfigBatch) {
    GpuConfig cfg = GpuConfig::Batch();
    EXPECT_EQ(cfg.device_id, 0);
    EXPECT_EQ(cfg.block_size, 256);
    EXPECT_EQ(cfg.memory_pool_bytes, static_cast<Size>(512ULL * 1024 * 1024));
}

TEST(GpuMC, GpuDeviceInfoDescribeFormat) {
    GpuDeviceInfo info;
    info.name = "Test GPU";
    info.device_id = 0;
    info.compute_capability_major = 8;
    info.compute_capability_minor = 9;
    info.total_global_memory = static_cast<Size>(8ULL * 1024 * 1024 * 1024);
    info.multiprocessor_count = 24;
    info.max_threads_per_block = 1024;
    info.warp_size = 32;
    info.supports_double = true;

    std::string desc = info.describe();
    EXPECT_NE(desc.find("Test GPU"), std::string::npos);
    EXPECT_NE(desc.find("SM 8.9"), std::string::npos);
    EXPECT_NE(desc.find("SMs=24"), std::string::npos);
    EXPECT_NE(desc.find("warp=32"), std::string::npos);
}

// ============================================================================
// 4-6: 设备查询 + 引擎构造
// ============================================================================

TEST(GpuMC, ListDevicesFindsAtLeastOne) {
    auto devices = gpu_list_devices();
    EXPECT_GE(devices.size(), 1u) << "主控站应有至少 1 个 GPU (RTX 4060)";
}

TEST(GpuMC, EngineConstruction) {
    GpuConfig cfg = GpuConfig::Default();
    GpuMCEngine engine(cfg);
    EXPECT_EQ(engine.config().device_id, 0);
    EXPECT_EQ(engine.config().block_size, 256);
}

TEST(GpuMC, EngineDeviceInfoRTX4060) {
    GpuMCEngine engine;
    auto info = engine.device_info();
    EXPECT_FALSE(info.name.empty()) << "设备名不应为空";
    EXPECT_EQ(info.compute_capability_major, 8) << "RTX 4060 = Ada = SM 8.x";
    EXPECT_EQ(info.compute_capability_minor, 9) << "RTX 4060 = SM 8.9";
    EXPECT_TRUE(info.supports_double) << "Ada 支持 FP64";
    EXPECT_EQ(info.warp_size, 32) << "NVIDIA warp 恒为 32";
    EXPECT_GT(info.multiprocessor_count, 0);
}

// ============================================================================
// 7-10: MC 定价 vs BSM 解析解
// ============================================================================

TEST(GpuMC, EuropeanCall_ATM_MatchesBSM) {
    GpuMCEngine engine;
    GpuBSMOption opt{100.0, 100.0, 1.0, 0.05, 0.0, 0.20, true};
    Size n_paths = 500000;
    auto result = engine.price_european_gbm(opt, n_paths, 42);

    double expected = bsm_price(100, 100, 1, 0.05, 0, 0.20, true);
    double tol = 4.0 * result.standard_error;

    EXPECT_NEAR(result.price, expected, tol)
        << "ATM call: MC=" << result.price << " BSM=" << expected
        << " SE=" << result.standard_error << " tol=" << tol;
    EXPECT_TRUE(result.used_gpu);
    EXPECT_EQ(result.n_paths, n_paths);
}

TEST(GpuMC, EuropeanPut_ATM_MatchesBSM) {
    GpuMCEngine engine;
    GpuBSMOption opt{100.0, 100.0, 1.0, 0.05, 0.0, 0.20, false};
    Size n_paths = 500000;
    auto result = engine.price_european_gbm(opt, n_paths, 42);

    double expected = bsm_price(100, 100, 1, 0.05, 0, 0.20, false);
    double tol = 4.0 * result.standard_error;

    EXPECT_NEAR(result.price, expected, tol)
        << "ATM put: MC=" << result.price << " BSM=" << expected
        << " SE=" << result.standard_error << " tol=" << tol;
}

TEST(GpuMC, EuropeanCall_ITM_MatchesBSM) {
    GpuMCEngine engine;
    // ITM call: S=120 > K=100
    GpuBSMOption opt{120.0, 100.0, 1.0, 0.05, 0.0, 0.20, true};
    Size n_paths = 500000;
    auto result = engine.price_european_gbm(opt, n_paths, 42);

    double expected = bsm_price(120, 100, 1, 0.05, 0, 0.20, true);
    double tol = 4.0 * result.standard_error;

    EXPECT_NEAR(result.price, expected, tol)
        << "ITM call: MC=" << result.price << " BSM=" << expected
        << " SE=" << result.standard_error << " tol=" << tol;
}

TEST(GpuMC, EuropeanCall_OTM_MatchesBSM) {
    GpuMCEngine engine;
    // OTM call: S=80 < K=100
    GpuBSMOption opt{80.0, 100.0, 1.0, 0.05, 0.0, 0.20, true};
    Size n_paths = 500000;
    auto result = engine.price_european_gbm(opt, n_paths, 42);

    double expected = bsm_price(80, 100, 1, 0.05, 0, 0.20, true);
    double tol = 5.0 * result.standard_error;  // OTM 方差大, 放宽

    EXPECT_NEAR(result.price, expected, tol)
        << "OTM call: MC=" << result.price << " BSM=" << expected
        << " SE=" << result.standard_error << " tol=" << tol;
}

// ============================================================================
// 11-12: 确定性
// ============================================================================

TEST(GpuMC, Determinism_SameSeed_SameResult) {
    GpuMCEngine engine;
    GpuBSMOption opt{100.0, 100.0, 1.0, 0.05, 0.0, 0.20, true};
    Size n_paths = 100000;

    auto r1 = engine.price_european_gbm(opt, n_paths, 42);
    auto r2 = engine.price_european_gbm(opt, n_paths, 42);

    EXPECT_DOUBLE_EQ(r1.price, r2.price)
        << "同 seed 应产生位精确相同结果";
    EXPECT_DOUBLE_EQ(r1.standard_error, r2.standard_error);
    EXPECT_DOUBLE_EQ(r1.min_payoff, r2.min_payoff);
    EXPECT_DOUBLE_EQ(r1.max_payoff, r2.max_payoff);
}

TEST(GpuMC, DifferentSeed_CloseResult) {
    GpuMCEngine engine;
    GpuBSMOption opt{100.0, 100.0, 1.0, 0.05, 0.0, 0.20, true};
    Size n_paths = 500000;

    auto r1 = engine.price_european_gbm(opt, n_paths, 42);
    auto r2 = engine.price_european_gbm(opt, n_paths, 999);

    // 不同 seed 应产生不同但接近的结果 (差异 < 5*SE_combined)
    double diff = std::abs(r1.price - r2.price);
    double se_combined = std::sqrt(r1.standard_error * r1.standard_error
                                   + r2.standard_error * r2.standard_error);
    EXPECT_LT(diff, 5.0 * se_combined)
        << "不同 seed: diff=" << diff << " 5*SE_combined=" << 5.0 * se_combined;
}

// ============================================================================
// 13: 标准误差收敛
// ============================================================================

TEST(GpuMC, StandardError_Decreases_WithMorePaths) {
    GpuMCEngine engine;
    GpuBSMOption opt{100.0, 100.0, 1.0, 0.05, 0.0, 0.20, true};

    Size n1 = 62500;   // N
    Size n2 = 250000;  // 4N → SE 应降为 ~1/2

    auto r1 = engine.price_european_gbm(opt, n1, 42);
    auto r2 = engine.price_european_gbm(opt, n2, 42);

    EXPECT_GT(r1.standard_error, 0.0);
    EXPECT_GT(r2.standard_error, 0.0);

    double ratio = r1.standard_error / r2.standard_error;
    // 理论 ratio = sqrt(4N/N) = 2, 容忍 [1.5, 2.8] 避免偶发
    EXPECT_GT(ratio, 1.5)
        << "SE 应随路径数减少: ratio=" << ratio;
    EXPECT_LT(ratio, 2.8)
        << "SE 收敛比不应过大: ratio=" << ratio;
}

// ============================================================================
// 14: 批量定价
// ============================================================================

TEST(GpuMC, BatchPricing_MultipleOptions) {
    GpuMCEngine engine;
    std::vector<GpuBSMOption> opts = {
        {100.0, 100.0, 1.0, 0.05, 0.0, 0.20, true},   // ATM call
        {100.0, 100.0, 1.0, 0.05, 0.0, 0.20, false},  // ATM put
        {120.0, 100.0, 0.5, 0.03, 0.0, 0.25, true},   // ITM call short T
        {80.0,  100.0, 2.0, 0.05, 0.0, 0.15, true},   // OTM call long T
        {100.0, 90.0,  0.25, 0.05, 0.02, 0.30, true}, // ITM call, dividend
    };
    Size n_paths = 200000;

    auto result = engine.price_european_gbm_batch(opts, n_paths, 42);

    ASSERT_EQ(result.prices.size(), opts.size());
    ASSERT_EQ(result.standard_errors.size(), opts.size());

    for (Size i = 0; i < opts.size(); ++i) {
        double expected = bsm_price(opts[i].S, opts[i].K, opts[i].T,
                                     opts[i].r, opts[i].q, opts[i].sigma,
                                     opts[i].is_call);
        double tol = 5.0 * result.standard_errors[i];
        EXPECT_NEAR(result.prices[i], expected, tol)
            << "Option " << i << ": MC=" << result.prices[i]
            << " BSM=" << expected << " SE=" << result.standard_errors[i];
    }

    EXPECT_GT(result.throughput_paths_per_sec, 0.0);
}

// ============================================================================
// 15: 性能 (1M 路径)
// ============================================================================

TEST(GpuMC, Performance_1MPaths_CompletesQuickly) {
    GpuMCEngine engine;
    GpuBSMOption opt{100.0, 100.0, 1.0, 0.05, 0.0, 0.20, true};
    Size n_paths = 1'000'000;

    auto result = engine.price_european_gbm(opt, n_paths, 42);

    // RTX 4060 应在合理时间内完成 1M 路径 (含 H2D/D2H)
    // 保守上界: 2 秒 (含首次 CUDA 初始化开销)
    EXPECT_LT(result.total_ms, 2000.0)
        << "1M 路径总耗时: " << result.total_ms << " ms";

    double expected = bsm_price(100, 100, 1, 0.05, 0, 0.20, true);
    double tol = 4.0 * result.standard_error;
    EXPECT_NEAR(result.price, expected, tol)
        << "1M paths: MC=" << result.price << " BSM=" << expected;
}

// ============================================================================
// 16-19: RISK-025/026 回归
// ============================================================================

namespace {
/// 构造即析构助手 (避开 EXPECT_THROW 内的声明语句怪癖)
void construct_engine(GpuConfig cfg) {
    GpuMCEngine engine(cfg);
    (void)engine;
}
}  // namespace

TEST(GpuMC, EngineRejectsInvalidBlockSize) {
    // RISK-026: 此前 block_size 未校验 → 非 32 倍数会触发 warp shuffle 全掩码 UB,
    // 超上限会 launch 失败但静默返回 price=0 且 used_gpu=true
    GpuConfig cfg = GpuConfig::Default();
    cfg.block_size = 100;  // 非 32 倍数
    EXPECT_THROW(construct_engine(cfg), std::invalid_argument);

    GpuConfig cfg_zero = GpuConfig::Default();
    cfg_zero.block_size = 0;
    EXPECT_THROW(construct_engine(cfg_zero), std::invalid_argument);

    GpuConfig cfg_neg = GpuConfig::Default();
    cfg_neg.block_size = -32;  // 负值: C++ 取模为 0, 必须显式拒绝
    EXPECT_THROW(construct_engine(cfg_neg), std::invalid_argument);

    GpuConfig cfg_big = GpuConfig::Default();
    cfg_big.block_size = 4096;  // 超过设备 maxThreadsPerBlock
    EXPECT_THROW(construct_engine(cfg_big), std::invalid_argument);
}

TEST(GpuMC, EngineRejectsUnknownDevice) {
    // RISK-026: 此前 cudaSetDevice 失败无检查 → 伪结果 (0 价格 + used_gpu=true)
    GpuConfig cfg = GpuConfig::Default();
    cfg.device_id = 9999;
    EXPECT_THROW(construct_engine(cfg), std::runtime_error);
}

TEST(GpuMC, EngineDestructorReleasesImpl) {
    // RISK-025: PIMPL 析构此前永不执行 (裸指针 + default 析构) → 每个引擎泄漏
    // 一个 Impl 与一个 CUDA stream
    EXPECT_EQ(GpuMCEngine::live_impl_count(), 0u)
        << "前序用例的引擎应已全部析构";

    {
        GpuMCEngine engine;
        EXPECT_EQ(GpuMCEngine::live_impl_count(), 1u);
    }
    EXPECT_EQ(GpuMCEngine::live_impl_count(), 0u)
        << "引擎离开作用域后 Impl 必须被释放";

    for (int i = 0; i < 64; ++i) {
        GpuMCEngine engine;
    }
    EXPECT_EQ(GpuMCEngine::live_impl_count(), 0u)
        << "循环构造 64 个引擎后不得残留 Impl";
}

TEST(GpuMC, EngineMoveTransfersOwnership) {
    // RISK-025: 默认移动此前只浅拷裸指针 → 移出对象与移入对象共享同一 Impl
    EXPECT_EQ(GpuMCEngine::live_impl_count(), 0u);

    {
        GpuMCEngine src;
        GpuMCEngine dst(std::move(src));
        EXPECT_EQ(GpuMCEngine::live_impl_count(), 1u)
            << "移动构造不产生新 Impl";

        GpuMCEngine dst2;
        dst2 = std::move(dst);
        EXPECT_EQ(GpuMCEngine::live_impl_count(), 1u)
            << "移动赋值应释放 dst2 原有 Impl";
        EXPECT_EQ(dst2.config().block_size, GpuConfig::Default().block_size);
    }
    EXPECT_EQ(GpuMCEngine::live_impl_count(), 0u)
        << "移动后的所有对象析构后不得残留";
}

// ============================================================================
// 20-23: RISK-027 RNG 对齐 (设备端 Random123 vs 官方 KAT / CPU 流 / 正态层级)
// ============================================================================

// P2a: 设备端 Random123 与官方 KAT 一致（证明 nvcc 编译出同一实现）
TEST(GpuMC, PhiloxMatchesKnownAnswers) {
    std::vector<gpu_probe::KatVector> vecs(rng_kat::kKatCount);
    for (std::size_t i = 0; i < rng_kat::kKatCount; ++i) {
        for (int j = 0; j < 4; ++j) vecs[i].ctr[j] = rng_kat::kKatVectors[i].ctr[j];
        vecs[i].key[0] = rng_kat::kKatVectors[i].key[0];
        vecs[i].key[1] = rng_kat::kKatVectors[i].key[1];
    }
    const auto got = gpu_probe::philox_kat(vecs);
    for (std::size_t i = 0; i < rng_kat::kKatCount; ++i) {
        for (int j = 0; j < 4; ++j) {
            ASSERT_EQ(got[i * 4 + j], rng_kat::kKatVectors[i].answer[j])
                << "vector " << i << " word " << j;
        }
    }
}

// P2b: uniform 流 CPU/GPU 逐字相同（把 E1-9 固化为回归; U-1 分层承诺第 1 层）
TEST(GpuMC, UniformStreamMatchesCpu) {
    constexpr std::size_t kN = 10000;
    Philox4x64 cpu_rng(12345);
    const auto gpu = gpu_probe::philox_blocks(12345, kN);
    for (std::size_t p = 0; p < kN; ++p) {
        const auto blk = cpu_rng.next4();
        for (int j = 0; j < 4; ++j) {
            ASSERT_EQ(blk[j], gpu[p * 4 + j]) << "block " << p << " word " << j;
        }
    }
}

// P3: 正态层级 1e-15 相对差（U-1 分层承诺第 2 层; 回归锁, 非红驱动）
TEST(GpuMC, NormalSamplesWithinTolerance) {
    constexpr std::size_t kN = 10000;
    Philox4x64 cpu_rng(12345);
    const auto gpu = gpu_probe::box_muller_z1_batch(12345, kN);
    for (std::size_t p = 0; p < kN; ++p) {
        auto u = cpu_rng.next4();
        double u1 = (u[0] >> 11) * (1.0 / 9007199254740992.0);
        double u2 = (u[1] >> 11) * (1.0 / 9007199254740992.0);
        const double zc1 = box_muller(u1, u2).first;  // 只用 z1, 与 GPU 一致
        const double diff = std::fabs(zc1 - gpu[p]);
        const double denom = std::fabs(zc1);
        ASSERT_TRUE(denom > 1.0 ? diff / denom <= 1e-15 : diff <= 1e-15)
            << "path " << p << " cpu=" << zc1 << " gpu=" << gpu[p];
    }
}
