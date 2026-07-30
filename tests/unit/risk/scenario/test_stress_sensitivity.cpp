#include <gtest/gtest.h>
#include <vector>
#include <map>
#include <string>
#include <cmath>

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/risk/scenario/stress_test.hpp"
#include "cpphub/risk/scenario/sensitivity.hpp"

namespace cpphub {
inline namespace v1 {

static Real linear_portfolio(const std::map<std::string, Real>& factors) {
    Real val = 0;
    auto it = factors.find("Equity");
    if (it != factors.end()) val += it->second * 1000;
    it = factors.find("Bond");
    if (it != factors.end()) val += it->second * 500;
    it = factors.find("IR");
    if (it != factors.end()) val += (1.0 - it->second) * 200;
    it = factors.find("Vol");
    if (it != factors.end()) val -= it->second * 100;
    it = factors.find("CreditSpread");
    if (it != factors.end()) val -= it->second * 300;
    it = factors.find("Basis");
    if (it != factors.end()) val -= it->second * 150;
    return val;
}

static Real quadratic_portfolio(const std::map<std::string, Real>& factors) {
    Real val = 0;
    auto it = factors.find("Equity");
    if (it != factors.end()) val += it->second * 1000 + it->second * it->second * 5;
    it = factors.find("IR");
    if (it != factors.end()) val += (1.0 - it->second) * 200;
    it = factors.find("Vol");
    if (it != factors.end()) val -= it->second * 100;
    return val;
}

TEST(StressTest, SingleScenario) {
    std::map<std::string, Real> current = {{"Equity", 100.0}, {"Bond", 100.0}, {"IR", 0.05}, {"Vol", 0.20}};
    StressTester tester(linear_portfolio, current);

    StressScenario scenario;
    scenario.name = "Equity Crash";
    scenario.spot_shocks["Equity"] = -0.30;
    scenario.probability = 1;

    auto result = tester.run(scenario);
    Real expected_loss = 100.0 * 1000 * (-0.30);
    EXPECT_NEAR(result.total_pnl, expected_loss, 1e-6);
    EXPECT_LT(result.total_pnl, 0);
}

TEST(StressTest, MultiFactorScenario) {
    std::map<std::string, Real> current = {{"Equity", 100.0}, {"IR", 0.05}, {"Vol", 0.20}};
    StressTester tester(linear_portfolio, current);

    StressScenario scenario;
    scenario.name = "Multi Shock";
    scenario.spot_shocks["Equity"] = -0.30;
    scenario.rate_shocks["IR"] = 200;
    scenario.vol_shocks["Vol"] = 0.50;
    scenario.probability = 1;

    auto result = tester.run(scenario);
    Real eq_loss = 100.0 * 1000 * (-0.30);
    Real ir_loss = (1.0 - (0.05 + 200.0/10000.0)) * 200 - (1.0 - 0.05) * 200;
    Real vol_loss = -(0.20 * 1.50) * 100 - (-0.20 * 100);
    Real expected = eq_loss + ir_loss + vol_loss;
    EXPECT_NEAR(result.total_pnl, expected, 1e-6);
    EXPECT_LT(result.total_pnl, 0);
}

TEST(StressTest, Crisis2008) {
    std::map<std::string, Real> current = {{"Equity", 100.0}, {"Bond", 100.0}, {"IR", 0.05}, {"Vol", 0.20}};
    StressTester tester(linear_portfolio, current);

    auto crisis = StressTester::crisis_2008();
    std::vector<HistoricalCrisis> crises = {crisis};
    auto results = tester.run_historical(crises);

    ASSERT_EQ(results.size(), 1u);
    EXPECT_LT(results[0].total_pnl, 0);
}

TEST(StressTest, Crisis2020) {
    std::map<std::string, Real> current = {{"Equity", 100.0}, {"Bond", 100.0}, {"IR", 0.05}, {"Vol", 0.20}};
    StressTester tester(linear_portfolio, current);

    auto crisis = StressTester::crisis_2020();
    std::vector<HistoricalCrisis> crises = {crisis};
    auto results = tester.run_historical(crises);

    ASSERT_EQ(results.size(), 1u);
    EXPECT_LT(results[0].total_pnl, 0);
}

TEST(StressTest, ScenarioWeightedES) {
    std::map<std::string, Real> current = {{"Equity", 100.0}, {"IR", 0.05}, {"Vol", 0.20}};
    StressTester tester(linear_portfolio, current);

    StressScenario mild;
    mild.name = "Mild Down";
    mild.spot_shocks["Equity"] = -0.05;
    mild.probability = 0.5;

    StressScenario severe;
    severe.name = "Severe Down";
    severe.spot_shocks["Equity"] = -0.30;
    severe.probability = 0.3;

    StressScenario crash;
    crash.name = "Crash";
    crash.spot_shocks["Equity"] = -0.50;
    crash.probability = 0.2;

    std::vector<StressScenario> scenarios = {mild, severe, crash};
    Real es = tester.scenario_weighted_es(scenarios, 0.95);

    EXPECT_LT(es, 0);
    EXPECT_GT(es, -60000);
}

TEST(StressTest, BaselFRTBScenarios) {
    auto scenarios = StressTester::basel_frtb_scenarios();
    EXPECT_GE(scenarios.size(), 7u);

    bool has_equity_up = false, has_equity_down = false;
    bool has_rate_up = false, has_rate_down = false;
    for (const auto& s : scenarios) {
        if (s.name == "Equity Up") has_equity_up = true;
        if (s.name == "Equity Down") has_equity_down = true;
        if (s.name == "Rate Up") has_rate_up = true;
        if (s.name == "Rate Down") has_rate_down = true;
    }
    EXPECT_TRUE(has_equity_up);
    EXPECT_TRUE(has_equity_down);
    EXPECT_TRUE(has_rate_up);
    EXPECT_TRUE(has_rate_down);
}

TEST(Sensitivity, SingleFactor) {
    std::map<std::string, Real> current = {{"Equity", 100.0}, {"IR", 0.05}};
    SensitivityAnalysis sa(linear_portfolio, current);

    auto results = sa.single_factor({"Equity"}, 0.01);

    ASSERT_EQ(results.size(), 1u);
    EXPECT_EQ(results[0].factor_name, "Equity");
    EXPECT_NEAR(results[0].shock, 0.01, 1e-10);
    Real expected_pnl = 100.0 * 1000 * 0.01;
    EXPECT_NEAR(results[0].pnl, expected_pnl, 1e-6);
}

TEST(Sensitivity, MultiFactor) {
    std::map<std::string, Real> current = {{"Equity", 100.0}, {"IR", 0.05}};
    SensitivityAnalysis sa(linear_portfolio, current);

    std::vector<Real> levels = {-0.05, 0.0, 0.05};
    auto results = sa.multi_factor({"Equity", "IR"}, levels);

    ASSERT_EQ(results.size(), 2u);
    ASSERT_EQ(results[0].size(), 3u);
    ASSERT_EQ(results[1].size(), 3u);

    EXPECT_EQ(results[0][0].shock, -0.05);
    EXPECT_NEAR(results[0][0].pnl, 100.0 * 1000 * (-0.05), 1e-6);
    EXPECT_EQ(results[0][1].shock, 0.0);
    EXPECT_NEAR(results[0][1].pnl, 0.0, 1e-6);
    EXPECT_EQ(results[0][2].shock, 0.05);
    EXPECT_NEAR(results[0][2].pnl, 100.0 * 1000 * 0.05, 1e-6);
}

TEST(Sensitivity, Directional) {
    std::map<std::string, Real> current = {{"Equity", 100.0}, {"IR", 0.05}};
    SensitivityAnalysis sa(linear_portfolio, current);

    std::map<std::string, Real> direction = {{"Equity", -1.0}};
    std::vector<Real> magnitudes = {0.01, 0.02, 0.03};
    auto results = sa.directional(direction, magnitudes);

    ASSERT_EQ(results.size(), 3u);
    EXPECT_LT(results[0].pnl, 0);
    EXPECT_LT(results[1].pnl, results[0].pnl);
    EXPECT_LT(results[2].pnl, results[1].pnl);
}

TEST(Sensitivity, WorstCase) {
    std::map<std::string, Real> current = {{"Equity", 100.0}, {"IR", 0.05}};
    SensitivityAnalysis sa(linear_portfolio, current);

    std::map<std::string, std::pair<Real, Real>> ranges;
    ranges["Equity"] = {-0.10, 0.10};
    ranges["IR"] = {-0.02, 0.02};

    Real worst = sa.worst_case(ranges, 1000);

    EXPECT_LT(worst, 0);

    Real avg = 0;
    {
        auto results = sa.single_factor({"Equity"}, 0.0);
        avg = results[0].pnl;
    }
    EXPECT_LT(worst, avg);
}

} // namespace v1
} // namespace cpphub
