#pragma once
#include <cmath>
#include <stdexcept>
#include "cpphub/core/types.hpp"
#include "cpphub/core/error.hpp"

namespace cpphub {
inline namespace v1 {

template <Size N>
class Vector {
    Real data_[N]{};
public:
    constexpr Real& operator[](Size i) noexcept { return data_[i]; }
    constexpr const Real& operator[](Size i) const noexcept { return data_[i]; }
    constexpr Size size() const noexcept { return N; }

    Vector() = default;

    Vector operator+(const Vector& rhs) const {
        Vector r;
        for (Size i = 0; i < N; ++i) r[i] = data_[i] + rhs[i];
        return r;
    }

    Vector operator-(const Vector& rhs) const {
        Vector r;
        for (Size i = 0; i < N; ++i) r[i] = data_[i] - rhs[i];
        return r;
    }

    Vector operator*(Real s) const {
        Vector r;
        for (Size i = 0; i < N; ++i) r[i] = data_[i] * s;
        return r;
    }

    friend Vector operator*(Real s, const Vector& v) {
        return v * s;
    }

    Real dot(const Vector& rhs) const {
        Real s = 0;
        for (Size i = 0; i < N; ++i) s += data_[i] * rhs[i];
        return s;
    }

    Real norm() const {
        return std::sqrt(dot(*this));
    }
};

template <Size R, Size C>
class Matrix {
    Real data_[R][C]{};
public:
    constexpr Real& operator()(Size i, Size j) noexcept { return data_[i][j]; }
    constexpr const Real& operator()(Size i, Size j) const noexcept { return data_[i][j]; }
    constexpr Size rows() const noexcept { return R; }
    constexpr Size cols() const noexcept { return C; }

    Matrix() = default;

    Vector<R> row(Size i) const {
        Vector<R> v;
        for (Size j = 0; j < C; ++j) v[j] = data_[i][j];
        return v;
    }

    Vector<C> col(Size j) const {
        Vector<C> v;
        for (Size i = 0; i < R; ++i) v[i] = data_[i][j];
        return v;
    }

    template <Size C2>
    Matrix<R, C2> operator*(const Matrix<C, C2>& rhs) const {
        Matrix<R, C2> r;
        for (Size i = 0; i < R; ++i) {
            for (Size k = 0; k < C; ++k) {
                for (Size j = 0; j < C2; ++j) {
                    r(i, j) += data_[i][k] * rhs(k, j);
                }
            }
        }
        return r;
    }

    Vector<R> operator*(const Vector<C>& v) const {
        Vector<R> r;
        for (Size i = 0; i < R; ++i) {
            for (Size j = 0; j < C; ++j) {
                r[i] += data_[i][j] * v[j];
            }
        }
        return r;
    }

    Matrix<C, R> transpose() const {
        Matrix<C, R> r;
        for (Size i = 0; i < R; ++i)
            for (Size j = 0; j < C; ++j)
                r(j, i) = data_[i][j];
        return r;
    }

    Real determinant() const {
        static_assert(R == C, "determinant requires square matrix");
        if constexpr (R == 1) {
            return data_[0][0];
        } else if constexpr (R == 2) {
            return data_[0][0] * data_[1][1] - data_[0][1] * data_[1][0];
        } else if constexpr (R == 3) {
            return data_[0][0] * (data_[1][1] * data_[2][2] - data_[1][2] * data_[2][1])
                 - data_[0][1] * (data_[1][0] * data_[2][2] - data_[1][2] * data_[2][0])
                 + data_[0][2] * (data_[1][0] * data_[2][1] - data_[1][1] * data_[2][0]);
        } else {
            Real det = 1;
            Matrix<R, C> tmp = *this;
            for (Size i = 0; i < R; ++i) {
                Size pivot = i;
                for (Size k = i + 1; k < R; ++k) {
                    if (std::abs(tmp(k, i)) > std::abs(tmp(pivot, i)))
                        pivot = k;
                }
                if (std::abs(tmp(pivot, i)) < 1e-30)
                    return 0;
                if (pivot != i) {
                    for (Size j = i; j < R; ++j)
                        std::swap(tmp(i, j), tmp(pivot, j));
                    det = -det;
                }
                det *= tmp(i, i);
                for (Size k = i + 1; k < R; ++k) {
                    Real factor = tmp(k, i) / tmp(i, i);
                    for (Size j = i + 1; j < R; ++j)
                        tmp(k, j) -= factor * tmp(i, j);
                }
            }
            return det;
        }
    }
};

template <Size N>
Matrix<N, N> cholesky(const Matrix<N, N>& A) {
    Matrix<N, N> L;
    for (Size j = 0; j < N; ++j) {
        Real s = 0;
        for (Size k = 0; k < j; ++k) {
            s += L(j, k) * L(j, k);
        }
        Real val = A(j, j) - s;
        if (val <= 0) {
            throw CppHubException("Cholesky decomposition failed: matrix not positive definite", ErrorCode::NumericalError);
        }
        L(j, j) = std::sqrt(val);
        for (Size i = j + 1; i < N; ++i) {
            Real s2 = 0;
            for (Size k = 0; k < j; ++k) {
                s2 += L(i, k) * L(j, k);
            }
            L(i, j) = (A(i, j) - s2) / L(j, j);
        }
    }
    return L;
}

template <Size N>
Vector<N> thomas_algorithm(const Vector<N>& a, const Vector<N>& b,
                           const Vector<N>& c, const Vector<N>& d) {
    Vector<N> x;
    Vector<N> cp;
    Vector<N> dp;

    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    for (Size i = 1; i < N; ++i) {
        Real denom = b[i] - a[i] * cp[i - 1];
        if (i < N - 1) cp[i] = c[i] / denom;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }

    x[N - 1] = dp[N - 1];
    for (Size i = N - 1; i-- > 0; ) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    return x;
}

}  // namespace v1
}  // namespace cpphub
