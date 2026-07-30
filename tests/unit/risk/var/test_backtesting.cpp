#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/risk/var/backtesting.hpp"

namespace cpphub {
inline namespace v1 {

TEST(KupiecPOF, CorrectCoverage) {
    auto result = KupiecPOF::test(2, 250, 0.99);
    EXPECT_EQ(result.n_violations, 2u);
    EXPECT_EQ(result.n_observations, 250u);
    EXPECT_NEAR(result.expected_rate, 0.01, 1e-10);
    EXPECT_GT(result.p_value, 0.05);
    EXPECT_FALSE(result.reject_null);
}

TEST(KupiecPOF, UnderestimatedRisk) {
    auto result = KupiecPOF::test(10, 250, 0.99);
    EXPECT_EQ(result.n_violations, 10u);
    EXPECT_LT(result.p_value, 0.01);
    EXPECT_TRUE(result.reject_null);
}

TEST(KupiecPOF, OverestimatedRisk) {
    auto result = KupiecPOF::test(0, 250, 0.99);
    EXPECT_EQ(result.n_violations, 0u);
    EXPECT_EQ(result.violation_rate, 0.0);
    EXPECT_GT(result.p_value, 0.01);
    EXPECT_LT(result.p_value, 0.05);
}

TEST(KupiecPOF, FromSeries) {
    std::vector<Real> var_series(250, 100.0);
    std::vector<Real> losses(250, 50.0);
    losses[0] = 110.0;
    losses[1] = 120.0;

    auto result = KupiecPOF::test(var_series, losses, 0.99);
    EXPECT_EQ(result.n_violations, 2u);

    auto direct = KupiecPOF::test(2, 250, 0.99);
    EXPECT_NEAR(result.p_value, direct.p_value, 1e-10);
    EXPECT_EQ(result.reject_null, direct.reject_null);
}

TEST(ChristoffersenIID, IndependentBreaches) {
    std::vector<Real> var_series(500, 100.0);
    std::vector<Real> losses(500, 50.0);
    for (Size i = 0; i < 500; i += 50) {
        losses[i] = 110.0;
    }

    auto result = ChristoffersenIID::test(var_series, losses);
    EXPECT_GT(result.p_value, 0.05);
}

TEST(ChristoffersenIID, ClusteredBreaches) {
    std::vector<Real> var_series(500, 100.0);
    std::vector<Real> losses(500, 50.0);
    for (Size i = 200; i < 210; ++i) {
        losses[i] = 110.0;
    }

    auto result = ChristoffersenIID::test(var_series, losses);
    EXPECT_LT(result.p_value, 0.05);
}

TEST(ChristoffersenIID, JointTest) {
    std::vector<Real> var_series(500, 100.0);
    std::vector<Real> losses(500, 50.0);
    for (Size i = 200; i < 210; ++i) {
        losses[i] = 110.0;
    }

    auto ind = ChristoffersenIID::test(var_series, losses);
    auto joint = ChristoffersenIID::joint_test(var_series, losses, 0.99);

    EXPECT_LE(ind.p_value, joint.p_value + 1e-6);
}

TEST(BaselTrafficLight, GreenZone) {
    for (Size v = 0; v <= 4; ++v) {
        auto result = BaselTrafficLight::assess(v, 250);
        EXPECT_EQ(result.zone, BaselZone::Green);
        EXPECT_EQ(result.capital_multiplier, 3.0);
    }
}

TEST(BaselTrafficLight, YellowZone) {
    Real expected_mult[] = {3.4, 3.5, 3.65, 3.75, 3.85};
    for (Size v = 5; v <= 9; ++v) {
        auto result = BaselTrafficLight::assess(v, 250);
        EXPECT_EQ(result.zone, BaselZone::Yellow);
        EXPECT_NEAR(result.capital_multiplier, expected_mult[v - 5], 1e-10);
    }
}

TEST(BaselTrafficLight, RedZone) {
    for (Size v = 10; v <= 15; ++v) {
        auto result = BaselTrafficLight::assess(v, 250);
        EXPECT_EQ(result.zone, BaselZone::Red);
        EXPECT_EQ(result.capital_multiplier, 4.0);
    }
}

} // namespace v1
} // namespace cpphub
