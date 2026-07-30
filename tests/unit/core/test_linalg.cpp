#include <gtest/gtest.h>
#include "cpphub/core/linalg.hpp"
#include "cpphub/core/error.hpp"
#include <stdexcept>
#include <cmath>

using namespace cpphub;

TEST(Linalg, VectorBasic) {
    Vector<3> v;
    v[0] = 1.0; v[1] = 2.0; v[2] = 3.0;
    EXPECT_EQ(v.size(), 3u);
    EXPECT_DOUBLE_EQ(v[0], 1.0);
    EXPECT_DOUBLE_EQ(v[2], 3.0);
}

TEST(Linalg, VectorDotProduct) {
    Vector<3> a, b;
    a[0] = 1; a[1] = 2; a[2] = 3;
    b[0] = 4; b[1] = 5; b[2] = 6;
    EXPECT_DOUBLE_EQ(a.dot(b), 32.0);
}

TEST(Linalg, VectorNorm) {
    Vector<3> v;
    v[0] = 3; v[1] = 4; v[2] = 0;
    EXPECT_DOUBLE_EQ(v.norm(), 5.0);
}

TEST(Linalg, VectorAdd) {
    Vector<3> a, b;
    a[0] = 1; a[1] = 2; a[2] = 3;
    b[0] = 4; b[1] = 5; b[2] = 6;
    auto c = a + b;
    EXPECT_DOUBLE_EQ(c[0], 5);
    EXPECT_DOUBLE_EQ(c[1], 7);
    EXPECT_DOUBLE_EQ(c[2], 9);
}

TEST(Linalg, VectorSub) {
    Vector<3> a, b;
    a[0] = 4; a[1] = 5; a[2] = 6;
    b[0] = 1; b[1] = 2; b[2] = 3;
    auto c = a - b;
    EXPECT_DOUBLE_EQ(c[0], 3);
    EXPECT_DOUBLE_EQ(c[1], 3);
    EXPECT_DOUBLE_EQ(c[2], 3);
}

TEST(Linalg, VectorScalarMul) {
    Vector<3> v;
    v[0] = 1; v[1] = 2; v[2] = 3;
    auto r = v * 2.0;
    EXPECT_DOUBLE_EQ(r[0], 2);
    EXPECT_DOUBLE_EQ(r[1], 4);
    EXPECT_DOUBLE_EQ(r[2], 6);
}

TEST(Linalg, MatrixBasic) {
    Matrix<2, 3> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0; m(0, 2) = 3.0;
    m(1, 0) = 4.0; m(1, 1) = 5.0; m(1, 2) = 6.0;
    EXPECT_EQ(m.rows(), 2u);
    EXPECT_EQ(m.cols(), 3u);
    EXPECT_DOUBLE_EQ(m(1, 2), 6.0);
}

TEST(Linalg, MatrixMultiply) {
    Matrix<2, 2> a, b, c;
    a(0,0)=1; a(0,1)=2; a(1,0)=3; a(1,1)=4;
    b(0,0)=5; b(0,1)=6; b(1,0)=7; b(1,1)=8;
    c = a * b;
    EXPECT_DOUBLE_EQ(c(0,0), 19);  // 1*5+2*7
    EXPECT_DOUBLE_EQ(c(0,1), 22);  // 1*6+2*8
    EXPECT_DOUBLE_EQ(c(1,0), 43);  // 3*5+4*7
    EXPECT_DOUBLE_EQ(c(1,1), 50);  // 3*6+4*8
}

TEST(Linalg, MatrixTranspose) {
    Matrix<2, 3> m;
    m(0,0)=1; m(0,1)=2; m(0,2)=3;
    m(1,0)=4; m(1,1)=5; m(1,2)=6;
    auto t = m.transpose();
    EXPECT_EQ(t.rows(), 3u);
    EXPECT_EQ(t.cols(), 2u);
    EXPECT_DOUBLE_EQ(t(0,0), 1);
    EXPECT_DOUBLE_EQ(t(1,0), 2);
    EXPECT_DOUBLE_EQ(t(2,1), 6);
}

TEST(Linalg, MatrixVectorMul) {
    Matrix<2, 2> m;
    m(0,0)=1; m(0,1)=2; m(1,0)=3; m(1,1)=4;
    Vector<2> v;
    v[0]=1; v[1]=1;
    auto r = m * v;
    EXPECT_DOUBLE_EQ(r[0], 3);
    EXPECT_DOUBLE_EQ(r[1], 7);
}

TEST(Linalg, CholeskyIdentity) {
    Matrix<3, 3> I{};
    I(0,0)=1; I(1,1)=1; I(2,2)=1;
    auto L = cholesky(I);
    EXPECT_NEAR(L(0,0), 1.0, 1e-15);
    EXPECT_NEAR(L(1,0), 0.0, 1e-15);
    EXPECT_NEAR(L(1,1), 1.0, 1e-15);
}

TEST(Linalg, CholeskySPD) {
    Matrix<2,2> A;
    A(0,0)=4; A(0,1)=2; A(1,0)=2; A(1,1)=3;
    auto L = cholesky(A);
    EXPECT_NEAR(L(0,0), 2.0, 1e-14);
    EXPECT_NEAR(L(1,0), 1.0, 1e-14);
    EXPECT_NEAR(L(1,1), std::sqrt(2.0), 1e-14);
}

TEST(Linalg, CholeskyNonSPDThrows) {
    Matrix<2,2> A;
    A(0,0)=-1; A(0,1)=0; A(1,0)=0; A(1,1)=-1;
    EXPECT_THROW(cholesky(A), CppHubException);
}

TEST(Linalg, ThomasAlgorithm) {
    Vector<3> a, b, c, d;
    a[0]=0; a[1]=-1; a[2]=-1;
    b[0]=2; b[1]=2; b[2]=2;
    c[0]=-1; c[1]=-1; c[2]=0;
    d[0]=1; d[1]=0; d[2]=1;
    auto x = thomas_algorithm(a, b, c, d);
    EXPECT_NEAR(x[0], 1.0, 1e-14);
    EXPECT_NEAR(x[1], 1.0, 1e-14);
    EXPECT_NEAR(x[2], 1.0, 1e-14);
}

TEST(Linalg, MatrixDeterminant2x2) {
    Matrix<2,2> m;
    m(0,0)=1; m(0,1)=2; m(1,0)=3; m(1,1)=4;
    EXPECT_DOUBLE_EQ(m.determinant(), -2.0);
}

TEST(Linalg, MatrixDeterminant3x3) {
    Matrix<3,3> m;
    m(0,0)=6; m(0,1)=1; m(0,2)=1;
    m(1,0)=4; m(1,1)=-2; m(1,2)=5;
    m(2,0)=2; m(2,1)=8; m(2,2)=7;
    EXPECT_DOUBLE_EQ(m.determinant(), -306.0);
}

TEST(Linalg, CholeskyProduct) {
    Matrix<3,3> A;
    A(0,0)=25; A(0,1)=15; A(0,2)=-5;
    A(1,0)=15; A(1,1)=18; A(1,2)=0;
    A(2,0)=-5; A(2,1)=0;  A(2,2)=11;
    auto L = cholesky(A);
    auto LLT = L * L.transpose();
    for (Size i = 0; i < 3; ++i)
        for (Size j = 0; j < 3; ++j)
            EXPECT_NEAR(LLT(i,j), A(i,j), 1e-12);
}
