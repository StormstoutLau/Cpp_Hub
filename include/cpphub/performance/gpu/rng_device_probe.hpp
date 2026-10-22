// GPU RNG 测试探针（RISK-027）: 把设备端 Philox 块/KAT/Box-Muller 结果导出到 host,
// 供 test_gpu_mc.cpp（.cpp, 不能含 __global__）做 CPU/GPU 对齐断言。
// 仅在 CPPHUB_CUDA_AVAILABLE 构建下有定义（gpu_mc.cu 末尾实现）;
// CPU stub 构建下未链接 —— 调用方（test_gpu_mc）本身也仅在 CUDA 分支构建。
#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>

namespace cpphub::v1::gpu_probe {

struct KatVector {
    std::uint64_t ctr[4];
    std::uint64_t key[2];
};

// 第 p 块 = Philox4x64(ctr={p,0,0,0}, key={seed,0}) 的 4 字;
// 返回 4*n_blocks 个字, 块优先顺序 (block0[0..3], block1[0..3], ...)。
std::vector<std::uint64_t> philox_blocks(std::uint64_t seed, std::size_t n_blocks);

// 逐 KAT 向量在设备端计算, 返回 4*vectors.size() 个字（每向量 4 字, 按输入序）。
std::vector<std::uint64_t> philox_kat(const std::vector<KatVector>& vectors);

// 与 MC 内核同一 uint64_to_unit + gpu_box_muller 路径的 z1 批量结果（第 p 路径）。
std::vector<double> box_muller_z1_batch(std::uint64_t seed, std::size_t n);

}  // namespace cpphub::v1::gpu_probe
