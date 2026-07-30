#pragma once
#include <cmath>
#include <cstring>
#include <array>
#include "cpphub/core/config.hpp"

// immintrin.h 必须在全局命名空间中 include，否则 __m128d/__m256d 等内建类型
// 会被拉入 cpphub::v1 命名空间，污染 GCC <random> 等系统头文件的类型查找
// (GCC bug: opt_random.h 引用 __m128d 时 ADL 找到 cpphub::v1::__m128d)
#if defined(CPPHUB_HAS_AVX2)
#include <immintrin.h>
#endif

namespace cpphub {
inline namespace v1 {

#if defined(CPPHUB_HAS_AVX2)

struct f64x4 {
    __m256d v;
    f64x4() = default;
    explicit f64x4(__m256d x) : v(x) {}
    f64x4(double x) : v(_mm256_set1_pd(x)) {}
    f64x4(double a, double b, double c, double d) : v(_mm256_set_pd(d, c, b, a)) {}
};

inline f64x4 load(const double* p) noexcept {
    return f64x4(_mm256_load_pd(p));
}

inline f64x4 loadu(const double* p) noexcept {
    return f64x4(_mm256_loadu_pd(p));
}

inline void store(double* p, f64x4 x) noexcept {
    _mm256_store_pd(p, x.v);
}

inline void storeu(double* p, f64x4 x) noexcept {
    _mm256_storeu_pd(p, x.v);
}

inline f64x4 operator+(f64x4 a, f64x4 b) noexcept {
    return f64x4(_mm256_add_pd(a.v, b.v));
}

inline f64x4 operator-(f64x4 a, f64x4 b) noexcept {
    return f64x4(_mm256_sub_pd(a.v, b.v));
}

inline f64x4 operator*(f64x4 a, f64x4 b) noexcept {
    return f64x4(_mm256_mul_pd(a.v, b.v));
}

inline f64x4 operator/(f64x4 a, f64x4 b) noexcept {
    return f64x4(_mm256_div_pd(a.v, b.v));
}

inline double hsum(f64x4 x) noexcept {
    __m128d lo = _mm256_castpd256_pd128(x.v);
    __m128d hi = _mm256_extractf128_pd(x.v, 1);
    __m128d sum = _mm_add_pd(lo, hi);
    __m128d hsum = _mm_hadd_pd(sum, sum);
    return _mm_cvtsd_f64(hsum);
}

inline f64x4 sqrt(f64x4 x) noexcept {
    return f64x4(_mm256_sqrt_pd(x.v));
}

inline f64x4 exp(f64x4 x) noexcept {
    alignas(32) double arr[4];
    store(arr, x);
    arr[0] = std::exp(arr[0]);
    arr[1] = std::exp(arr[1]);
    arr[2] = std::exp(arr[2]);
    arr[3] = std::exp(arr[3]);
    return load(arr);
}

inline f64x4 log(f64x4 x) noexcept {
    alignas(32) double arr[4];
    store(arr, x);
    arr[0] = std::log(arr[0]);
    arr[1] = std::log(arr[1]);
    arr[2] = std::log(arr[2]);
    arr[3] = std::log(arr[3]);
    return load(arr);
}

inline f64x4 blend(f64x4 a, f64x4 b, __m256d mask) noexcept {
    return f64x4(_mm256_blendv_pd(a.v, b.v, mask));
}

#else

struct f64x4 {
    alignas(32) double data_[4]{};
    f64x4() = default;
    explicit f64x4(double x) { data_[0] = data_[1] = data_[2] = data_[3] = x; }
    f64x4(double a, double b, double c, double d) {
        data_[0] = a; data_[1] = b; data_[2] = c; data_[3] = d;
    }
    double& operator[](int i) { return data_[i]; }
    const double& operator[](int i) const { return data_[i]; }
};

inline f64x4 load(const double* p) noexcept {
    f64x4 r;
    std::memcpy(r.data_, p, 32);
    return r;
}

inline f64x4 loadu(const double* p) noexcept {
    f64x4 r;
    std::memcpy(r.data_, p, 32);
    return r;
}

inline void store(double* p, f64x4 x) noexcept {
    std::memcpy(p, x.data_, 32);
}

inline void storeu(double* p, f64x4 x) noexcept {
    std::memcpy(p, x.data_, 32);
}

inline f64x4 operator+(f64x4 a, f64x4 b) noexcept {
    f64x4 r;
    for (int i = 0; i < 4; ++i) r.data_[i] = a.data_[i] + b.data_[i];
    return r;
}

inline f64x4 operator-(f64x4 a, f64x4 b) noexcept {
    f64x4 r;
    for (int i = 0; i < 4; ++i) r.data_[i] = a.data_[i] - b.data_[i];
    return r;
}

inline f64x4 operator*(f64x4 a, f64x4 b) noexcept {
    f64x4 r;
    for (int i = 0; i < 4; ++i) r.data_[i] = a.data_[i] * b.data_[i];
    return r;
}

inline f64x4 operator/(f64x4 a, f64x4 b) noexcept {
    f64x4 r;
    for (int i = 0; i < 4; ++i) r.data_[i] = a.data_[i] / b.data_[i];
    return r;
}

inline double hsum(f64x4 x) noexcept {
    return x.data_[0] + x.data_[1] + x.data_[2] + x.data_[3];
}

inline f64x4 sqrt(f64x4 x) noexcept {
    f64x4 r;
    for (int i = 0; i < 4; ++i) r.data_[i] = std::sqrt(x.data_[i]);
    return r;
}

inline f64x4 exp(f64x4 x) noexcept {
    f64x4 r;
    for (int i = 0; i < 4; ++i) r.data_[i] = std::exp(x.data_[i]);
    return r;
}

inline f64x4 log(f64x4 x) noexcept {
    f64x4 r;
    for (int i = 0; i < 4; ++i) r.data_[i] = std::log(x.data_[i]);
    return r;
}

// MSVC 不允许前向声明 __m256d 内建类型, 非 AVX2 分支用占位类型
struct f64x4_mask {};
inline f64x4 blend(f64x4 a, f64x4 b, f64x4_mask mask) noexcept {
    (void)a; (void)mask;
    f64x4 r;
    for (int i = 0; i < 4; ++i) r.data_[i] = b.data_[i];
    return r;
}

#endif

struct f64x8 {
    f64x4 lo, hi;
    f64x8() = default;
    f64x8(f64x4 a, f64x4 b) : lo(a), hi(b) {}
    f64x8(double x) : lo(x), hi(x) {}
};

}  // namespace v1
}  // namespace cpphub
