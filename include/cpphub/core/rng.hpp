#pragma once
#include <cstdint>
#include <array>
#include <utility>
#include <limits>
#include <Random123/philox.h>   // RISK-027: vendor Random123 (third_party/random123)
#include "cpphub/core/constants.hpp"
#include "cpphub/core/simd.hpp"
#include "cpphub/core/linalg.hpp"

namespace cpphub {
inline namespace v1 {

// Philox4x64-10 (Random123): counter 4xu64, key 2xu64, 每次调用产出 4 字。
// 本类为流式适配层: 逐字输出 + 块缓冲, 公共 API 不变;
// 输出流语义变更见 USAGE_GUIDE §3.1.3 (2026-10-22 R-027)。
// （operator() = 流顺序下一字; next4() = 恰 1 个完整块 = 与 GPU 每路径 1 块对齐）。
class Philox4x64 {
    r123::Philox4x64::ctr_type ctr_{};    // 4xu64 计数器（下一块的 counter）
    r123::Philox4x64::key_type key_{};    // 2xu64 密钥 = {seed, stream}
    r123::Philox4x64::ctr_type block_{};  // 当前缓冲块的 4 字
    int idx_ = 4;                          // 4 = 块已耗尽

public:
    using result_type = uint64_t;
    static constexpr uint64_t min() { return 0; }
    static constexpr uint64_t max() { return UINT64_MAX; }

    explicit Philox4x64(uint64_t seed = 0, uint64_t stream = 0) {
        key_.v[0] = seed;
        key_.v[1] = stream;
    }

    // 流顺序的下一字; 块耗尽时以当前 counter 计算整块并前移
    uint64_t operator()() {
        if (idx_ == 4) refill();
        return block_.v[idx_++];
    }

    // 恰 1 个完整块（4 字）。调用方若从块边界出发则保持块对齐
    // （normal_simd 即此用法 → CPU 第 p 步 = GPU 第 p 路径同一 counter 同一 4 字）。
    std::array<uint64_t, 4> next4() {
        return {operator()(), operator()(), operator()(), operator()()};
    }

    // 丢弃 n 个输出字（DESIGN §6.3: ctr += n/4, 再消费 n%4）
    void discard(uint64_t n) {
        if (n == 0) return;
        const uint64_t avail = (idx_ == 4) ? 0 : static_cast<uint64_t>(4 - idx_);
        if (n < avail) {
            idx_ += static_cast<int>(n);
            return;
        }
        n -= avail;
        idx_ = 4;
        const uint64_t blocks = n / 4;
        if (blocks > 0) add_to_ctr(blocks);
        const uint64_t rem = n % 4;
        if (rem > 0) {
            refill();
            idx_ = static_cast<int>(rem);
        }
        // rem == 0: 保持 idx_=4, 下次 operator() 触发 refill
    }

private:
    void refill() {
        block_ = r123::Philox4x64{}(ctr_, key_);
        bump_ctr();
        idx_ = 0;
    }
    void bump_ctr() {
        for (int i = 0; i < 4; ++i) {
            if (++ctr_.v[i] != 0) break;  // 256 位小端进位
        }
    }
    void add_to_ctr(uint64_t k) {
        for (int i = 0; i < 4 && k != 0; ++i) {
            const uint64_t old = ctr_.v[i];
            ctr_.v[i] = old + k;
            k = (ctr_.v[i] < old) ? 1 : 0;
        }
    }
};

inline std::pair<Real, Real> box_muller(Real u1, Real u2) noexcept {
    Real r = std::sqrt(-2 * std::log(u1));
    Real theta = 2 * PI * u2;
    Real z1 = r * std::cos(theta);
    Real z2 = r * std::sin(theta);
    return {z1, z2};
}

inline f64x4 normal_simd(Philox4x64& rng) noexcept {
    auto u = rng.next4();
    double u1 = (u[0] >> 11) * (1.0 / 9007199254740992.0);
    double u2 = (u[1] >> 11) * (1.0 / 9007199254740992.0);
    double u3 = (u[2] >> 11) * (1.0 / 9007199254740992.0);
    double u4 = (u[3] >> 11) * (1.0 / 9007199254740992.0);
    auto [z1, z2] = box_muller(u1, u2);
    auto [z3, z4] = box_muller(u3, u4);
    return f64x4(z1, z2, z3, z4);
}

template <Size N>
Vector<N> generate_correlated(Philox4x64& rng, const Matrix<N, N>& L) {
    Vector<N> z;
    for (Size i = 0; i < N; i += 2) {
        uint64_t r1 = rng();
        uint64_t r2 = rng();
        double u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
        double u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
        auto [z1, z2] = box_muller(u1, u2);
        z[i] = z1;
        if (i + 1 < N) z[i + 1] = z2;
    }
    return L * z;
}

}  // namespace v1
}  // namespace cpphub
