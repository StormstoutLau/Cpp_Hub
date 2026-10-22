// 生成头: 由 tmp/extract_kat.py 从上游 tests/kat_vectors 抽取
// (Random123 commit 9545ff6413f258be2f04c1d319d99aaef7521150, philox4x64 10 轮全部行), 请勿手改。
#pragma once
#include <cstdint>
#include <cstddef>

namespace cpphub {
inline namespace v1 {
namespace rng_kat {

struct KatVector {
    std::uint64_t ctr[4];
    std::uint64_t key[2];
    std::uint64_t answer[4];
};

inline constexpr KatVector kKatVectors[] = {
    {{0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x16554d9eca36314cULL, 0xdb20fe9d672d0fdcULL, 0xd7e772cee186176bULL, 0x7e68b68aec7ba23bULL}},
    {{0xffffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffffULL},
     {0xffffffffffffffffULL, 0xffffffffffffffffULL},
     {0x87b092c3013fe90bULL, 0x438c3c67be8d0224ULL, 0x9cc7d7c69cd777b6ULL, 0xa09caebf594f0ba0ULL}},
    {{0x243f6a8885a308d3ULL, 0x13198a2e03707344ULL, 0xa4093822299f31d0ULL, 0x082efa98ec4e6c89ULL},
     {0x452821e638d01377ULL, 0xbe5466cf34e90c6cULL},
     {0xa528f45403e61d95ULL, 0x38c72dbd566e9788ULL, 0xa5a1610e72fd18b5ULL, 0x57bd43b5e52b7fe6ULL}},
};
inline constexpr std::size_t kKatCount =
    sizeof(kKatVectors) / sizeof(kKatVectors[0]);

}  // namespace rng_kat
}  // namespace v1
}  // namespace cpphub
