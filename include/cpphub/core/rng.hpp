#pragma once
#include <cstdint>
#include <array>
#include <utility>
#include <limits>
#include "cpphub/core/constants.hpp"
#include "cpphub/core/simd.hpp"
#include "cpphub/core/linalg.hpp"

// MSVC 128 位乘法内建
#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace cpphub {
inline namespace v1 {

static constexpr uint64_t PHILOX_M4x64_0 = 0xD2B74407B1CE6E93ULL;
static constexpr uint64_t PHILOX_M4x64_1 = 0x9CCD9A8C5AE13B41ULL;
static constexpr uint64_t PHILOX_W32_0 = 0x9E3779B9;
static constexpr uint64_t PHILOX_W32_1 = 0xBB67AE85;

class Philox4x64 {
    uint64_t counter_[2]{0, 0};
    uint64_t key_[2]{0, 0};
public:
    using result_type = uint64_t;
    static constexpr uint64_t min() { return 0; }
    static constexpr uint64_t max() { return UINT64_MAX; }

    explicit Philox4x64(uint64_t seed = 0, uint64_t stream = 0) {
        key_[0] = seed;
        key_[1] = stream;
    }

    uint64_t operator()() {
        auto r = compute();
        inc_counter();
        return r[0];
    }

    void discard(uint64_t n) {
        uint64_t lo = counter_[0] + n;
        if (lo < counter_[0]) counter_[1]++;
        counter_[0] = lo;
    }

    std::array<uint64_t, 4> next4() {
        auto r0 = compute();
        inc_counter();
        auto r1 = compute();
        inc_counter();
        auto r2 = compute();
        inc_counter();
        auto r3 = compute();
        inc_counter();
        return {r0[0], r1[0], r2[0], r3[0]};
    }

private:
    void inc_counter() {
        counter_[0]++;
        if (counter_[0] == 0) counter_[1]++;
    }

    static uint64_t mulhi(uint64_t a, uint64_t b) {
#if defined(CPPHUB_COMPILER_MSVC) || defined(_MSC_VER)
        // MSVC 无 __uint128_t, 使用 _umul128 内建
        uint64_t hi, lo;
        lo = _umul128(a, b, &hi);
        (void)lo;
        return hi;
#else
        return (uint64_t)(((__uint128_t)a * b) >> 64);
#endif
    }

    static uint64_t mullo(uint64_t a, uint64_t b) {
        // 64 位乘法的低 64 位, 所有平台都是普通乘法
        return a * b;
    }

    std::array<uint64_t, 2> compute() const {
        uint64_t c0 = counter_[0];
        uint64_t c1 = counter_[1];
        uint64_t k0 = key_[0];
        uint64_t k1 = key_[1];

        for (int r = 0; r < 10; ++r) {
            uint64_t hi0 = mulhi(c0, PHILOX_M4x64_0);
            uint64_t lo0 = mullo(c0, PHILOX_M4x64_0);
            uint64_t hi1 = mulhi(c1, PHILOX_M4x64_1);
            uint64_t lo1 = mullo(c1, PHILOX_M4x64_1);
            c0 = hi1 ^ lo0 ^ k0;
            c1 = hi0 ^ lo1 ^ k1;
            k0 += PHILOX_W32_0;
            k1 += PHILOX_W32_1;
        }

        return {c0, c1};
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
