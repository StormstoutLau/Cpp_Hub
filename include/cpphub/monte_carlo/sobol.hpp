#pragma once
#include "cpphub/core/types.hpp"
#include <vector>
#include <array>
#include <cstdint>

namespace cpphub {
inline namespace v1 {

class SobolSequence {
public:
    explicit SobolSequence(Size dimension, uint64_t seed = 0)
        : dim_(dimension)
    {
        (void)seed;
        initialize_direction_vectors();
    }

    std::vector<Real> next()
    {
        auto result = operator()(count_);
        ++count_;
        return result;
    }

    std::vector<Real> operator()(uint64_t n) const
    {
        std::vector<Real> result(dim_);
        for (Size d = 0; d < dim_; ++d) {
            uint64_t x = 0;
            const auto& dir = direction_vectors_[d];
            uint64_t bits = n;
            for (int i = 0; i < 64; ++i) {
                if (bits & 1) {
                    x ^= dir[i];
                }
                bits >>= 1;
                if (bits == 0) break;
            }
            result[d] = Real(x) * 0x1p-64;
        }
        return result;
    }

    Size dimension() const noexcept { return dim_; }

private:
    Size dim_;
    uint64_t count_{0};
    std::vector<std::array<uint64_t, 64>> direction_vectors_;

    struct PolyData {
        Size degree;
        std::vector<uint64_t> a;     // a_1..a_{m-1}
        std::vector<uint64_t> init;  // initial m-values
    };

    static constexpr Size poly_table_size = 20;

    static const PolyData* poly_table()
    {
        static const PolyData table[poly_table_size] = {
            // dim 1: x + 1
            {1, {}, {1}},
            // dim 2: x^2 + x + 1
            {2, {1}, {1, 3}},
            // dim 3: x^3 + x + 1
            {3, {0, 1}, {1, 3, 7}},
            // dim 4: x^4 + x + 1
            {4, {0, 0, 1}, {1, 1, 5, 9}},
            // dim 5: x^4 + x^3 + x^2 + x + 1
            {4, {1, 1, 1}, {1, 3, 5, 11}},
            // dim 6: x^5 + x^2 + 1
            {5, {0, 1, 0, 0}, {1, 1, 3, 3, 5}},
            // dim 7: x^5 + x^4 + x^3 + x^2 + 1
            {5, {1, 1, 1, 0}, {1, 3, 5, 11, 3}},
            // dim 8: x^6 + x + 1
            {6, {0, 0, 0, 0, 1}, {1, 1, 1, 1, 1, 1}},
            // dim 9: x^6 + x^4 + x^3 + x + 1
            {6, {0, 1, 1, 0, 1}, {1, 3, 5, 7, 9, 11}},
            // dim 10: x^6 + x^5 + x^2 + x + 1
            {6, {1, 0, 0, 1, 1}, {1, 1, 7, 5, 3, 9}},
            // dim 11: x^7 + x + 1
            {7, {0, 0, 0, 0, 0, 1}, {1, 1, 1, 1, 1, 1, 1}},
            // dim 12: x^7 + x^3 + x^2 + x + 1
            {7, {0, 1, 1, 1, 0, 0}, {1, 3, 5, 7, 9, 11, 13}},
            // dim 13: x^7 + x^4 + x^3 + x^2 + 1
            {7, {0, 0, 1, 1, 1, 0}, {1, 1, 3, 3, 5, 5, 7}},
            // dim 14: x^7 + x^5 + x^4 + x^3 + 1
            {7, {1, 1, 1, 0, 0, 0}, {1, 3, 5, 11, 13, 7, 9}},
            // dim 15: x^8 + x^4 + x^3 + x^2 + 1
            {8, {0, 0, 1, 1, 1, 0, 0}, {1, 1, 3, 5, 7, 9, 11, 13}},
            // dim 16: x^8 + x^6 + x^5 + x^4 + 1
            {8, {1, 1, 1, 0, 0, 0, 0}, {1, 3, 5, 11, 3, 7, 9, 13}},
            // dim 17: x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + 1
            {8, {1, 1, 1, 1, 1, 1, 0}, {1, 1, 1, 1, 1, 1, 1, 1}},
            // dim 18: x^8 + x^7 + x^5 + x^3 + 1
            {8, {1, 0, 1, 0, 1, 0, 0}, {1, 3, 7, 5, 11, 9, 13, 15}},
            // dim 19: x^8 + x^7 + x^3 + x^2 + 1
            {8, {1, 0, 0, 1, 1, 0, 0}, {1, 3, 7, 13, 9, 11, 5, 15}},
            // dim 20: x^9 + x^4 + 1
            {9, {0, 0, 0, 1, 0, 0, 0, 0}, {1, 1, 1, 1, 1, 1, 1, 1, 1}},
        };
        return table;
    }

    void initialize_direction_vectors()
    {
        direction_vectors_.resize(dim_);
        for (Size d = 0; d < dim_; ++d) {
            auto& dir = direction_vectors_[d];
            if (d == 0) {
                for (int i = 0; i < 64; ++i) {
                    dir[i] = uint64_t(1) << (63 - i);
                }
            } else {
                Size idx = d - 1;
                if (idx >= poly_table_size) idx = poly_table_size - 1;
                const auto& p = poly_table()[idx];
                Size m = p.degree;
                // initial values: v[i] = init[i] << (63 - i)
                for (Size i = 0; i < m && i < 64; ++i) {
                    dir[i] = p.init[i] << (63 - i);
                }
                for (Size i = m; i < 64; ++i) {
                    uint64_t val = dir[i - m] ^ (dir[i - m] >> m);
                    for (Size k = 0; k < m - 1; ++k) {
                        if (p.a[k]) {
                            val ^= dir[i - 1 - k];
                        }
                    }
                    dir[i] = val;
                }
            }
        }
    }
};

}  // namespace v1
}  // namespace cpphub
