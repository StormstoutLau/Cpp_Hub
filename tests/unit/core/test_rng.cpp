#include "cpphub/core/rng.hpp"
#include "cpphub/core/linalg.hpp"
#include <gtest/gtest.h>
#include <vector>
#include <numeric>
#include <cmath>

using namespace cpphub;

TEST(Rng, PhiloxDeterministic) {
    Philox4x64 rng1(42);
    Philox4x64 rng2(42);
    for (int i = 0; i < 100; ++i) {
        EXPECT_EQ(rng1(), rng2());
    }
}

TEST(Rng, PhiloxDistinctSeeds) {
    Philox4x64 rng1(42);
    Philox4x64 rng2(43);
    int diff = 0;
    for (int i = 0; i < 100; ++i) {
        if (rng1() != rng2()) ++diff;
    }
    EXPECT_GT(diff, 95);
}

TEST(Rng, PhiloxNext4) {
    Philox4x64 rng(123);
    auto v1 = rng.next4();
    auto v2 = rng.next4();
    EXPECT_NE(v1[0], v1[1]);
    EXPECT_NE(v1[0], v1[2]);
    bool all_same = (v1[0]==v2[0] && v1[1]==v2[1] && v1[2]==v2[2] && v1[3]==v2[3]);
    EXPECT_FALSE(all_same);
}

TEST(Rng, BoxMuller) {
    auto [z1, z2] = box_muller(0.5, 0.5);
    EXPECT_NEAR(z1, -std::sqrt(2*std::log(2.0)), 1e-14);
    EXPECT_NEAR(z2, 0.0, 1e-14);
}

TEST(Rng, NormalSimdStats) {
    Philox4x64 rng(2024);
    std::vector<double> samples;
    samples.reserve(10000);
    for (int i = 0; i < 2500; ++i) {
        f64x4 z = normal_simd(rng);
        alignas(32) double arr[4];
        store(arr, z);
        for (int j = 0; j < 4; ++j) samples.push_back(arr[j]);
    }
    double mean = std::accumulate(samples.begin(), samples.end(), 0.0) / samples.size();
    double sq_sum = 0;
    for (double x : samples) sq_sum += (x - mean) * (x - mean);
    double var = sq_sum / samples.size();
    EXPECT_NEAR(mean, 0.0, 0.05);
    EXPECT_NEAR(var, 1.0, 0.1);
}

TEST(Rng, GenerateCorrelated) {
    Matrix<2,2> L;
    L(0,0)=1.0; L(0,1)=0.0;
    L(1,0)=0.5; L(1,1)=std::sqrt(3.0)/2.0;
    Philox4x64 rng(99);
    std::vector<Vector<2>> samples;
    for (int i = 0; i < 10000; ++i) {
        samples.push_back(generate_correlated<2>(rng, L));
    }
    double mean0 = 0, mean1 = 0;
    for (auto& v : samples) { mean0 += v[0]; mean1 += v[1]; }
    mean0 /= samples.size(); mean1 /= samples.size();
    double cov00=0, cov01=0, cov11=0;
    for (auto& v : samples) {
        cov00 += (v[0]-mean0)*(v[0]-mean0);
        cov01 += (v[0]-mean0)*(v[1]-mean1);
        cov11 += (v[1]-mean1)*(v[1]-mean1);
    }
    cov00 /= samples.size(); cov01 /= samples.size(); cov11 /= samples.size();
    EXPECT_NEAR(cov00, 1.0, 0.1);
    EXPECT_NEAR(cov01, 0.5, 0.1);
    EXPECT_NEAR(cov11, 1.0, 0.1);
}

TEST(Rng, PhiloxZeroSeed) {
    Philox4x64 rng(0);
    auto v = rng.next4();
    EXPECT_NE(v[0], 0);
    EXPECT_NE(v[1], 0);
}

TEST(Rng, PhiloxDiscard) {
    Philox4x64 rng1(42);
    Philox4x64 rng2(42);
    rng2.discard(50);
    for (int i = 0; i < 50; ++i) rng1();
    for (int i = 0; i < 50; ++i) {
        EXPECT_EQ(rng1(), rng2());
    }
}

TEST(Rng, PhiloxNext4Consistency) {
    Philox4x64 rng(77);
    uint64_t op1 = rng();
    uint64_t op2 = rng();
    uint64_t op3 = rng();
    uint64_t op4 = rng();
    Philox4x64 rng2(77);
    auto n4 = rng2.next4();
    EXPECT_EQ(op1, n4[0]);
    EXPECT_EQ(op2, n4[1]);
    EXPECT_EQ(op3, n4[2]);
    EXPECT_EQ(op4, n4[3]);
}
