// SOURCE: PHASE1_SPEC §2.4 - 基础类型别名
#pragma once
#include <complex>
#include <cstddef>
#include <cstdint>

namespace cpphub {
inline namespace v1 {

using Real = double;
using Complex = std::complex<double>;
using Size = std::size_t;
using Index = std::ptrdiff_t;

}  // namespace v1
}  // namespace cpphub
