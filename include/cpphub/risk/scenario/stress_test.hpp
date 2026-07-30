#pragma once
#include "cpphub/core/types.hpp"
#include <vector>
#include <map>
#include <string>
#include <functional>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace cpphub {
inline namespace v1 {

struct StressScenario {
    std::string name;
    std::map<std::string, Real> spot_shocks;
    std::map<std::string, Real> rate_shocks;
    std::map<std::string, Real> vol_shocks;
    std::map<std::pair<std::string, std::string>, Real> corr_shocks;
    Real probability;
    std::string description;
};

struct PortfolioPnL {
    Real total_pnl;
    std::map<std::string, Real> by_asset;
    std::map<std::string, Real> by_risk_factor;
};

struct HistoricalCrisis {
    std::string name;
    std::string start_date;
    std::string end_date;
    std::map<std::string, Real> asset_returns;
    std::map<std::string, Real> rate_changes;
    std::map<std::string, Real> vol_changes;
};

class StressTester {
public:
    using PortfolioValueFn = std::function<Real(const std::map<std::string, Real>&)>;

    StressTester(PortfolioValueFn value_fn, std::map<std::string, Real> current_factors)
        : value_fn_(std::move(value_fn)), current_factors_(std::move(current_factors)) {}

    PortfolioPnL run(const StressScenario& scenario) {
        auto shocked = apply_shocks(scenario);
        Real shocked_value = value_fn_(shocked);
        Real base_value = value_fn_(current_factors_);
        Real total_pnl = shocked_value - base_value;

        PortfolioPnL result;
        result.total_pnl = total_pnl;

        for (const auto& [asset, shock] : scenario.spot_shocks) {
            auto tmp = current_factors_;
            if (tmp.find(asset) != tmp.end()) {
                tmp[asset] = current_factors_.at(asset) * (Real(1) + shock);
            } else {
                tmp[asset] = shock;
            }
            Real val = value_fn_(tmp);
            result.by_asset[asset] = val - base_value;
        }

        for (const auto& [factor, shock] : scenario.rate_shocks) {
            auto tmp = current_factors_;
            if (tmp.find(factor) != tmp.end()) {
                tmp[factor] = current_factors_.at(factor) + shock / Real(10000);
            } else {
                tmp[factor] = shock;
            }
            Real val = value_fn_(tmp);
            result.by_risk_factor[factor + "_rate"] = val - base_value;
        }

        for (const auto& [factor, shock] : scenario.vol_shocks) {
            auto tmp = current_factors_;
            if (tmp.find(factor) != tmp.end()) {
                tmp[factor] = current_factors_.at(factor) * (Real(1) + shock);
            } else {
                tmp[factor] = shock;
            }
            Real val = value_fn_(tmp);
            result.by_risk_factor[factor + "_vol"] = val - base_value;
        }

        return result;
    }

    std::vector<PortfolioPnL> run_historical(const std::vector<HistoricalCrisis>& crises) {
        std::vector<PortfolioPnL> results;
        results.reserve(crises.size());
        for (const auto& crisis : crises) {
            StressScenario scenario;
            scenario.name = crisis.name;
            for (const auto& [asset, ret] : crisis.asset_returns) {
                scenario.spot_shocks[asset] = ret;
            }
            for (const auto& [rate, chg] : crisis.rate_changes) {
                scenario.rate_shocks[rate] = chg * Real(10000);
            }
            for (const auto& [vol, chg] : crisis.vol_changes) {
                scenario.vol_shocks[vol] = chg;
            }
            scenario.probability = Real(1);
            results.push_back(run(scenario));
        }
        return results;
    }

    static std::vector<StressScenario> basel_frtb_scenarios() {
        std::vector<StressScenario> scenarios;

        StressScenario equity_up;
        equity_up.name = "Equity Up";
        equity_up.spot_shocks["Equity"] = 0.30;
        equity_up.probability = 1;
        equity_up.description = "Equity spot +30%";
        scenarios.push_back(equity_up);

        StressScenario equity_down;
        equity_down.name = "Equity Down";
        equity_down.spot_shocks["Equity"] = -0.30;
        equity_down.probability = 1;
        equity_down.description = "Equity spot -30%";
        scenarios.push_back(equity_down);

        StressScenario rate_up;
        rate_up.name = "Rate Up";
        rate_up.rate_shocks["IR"] = 200;
        rate_up.probability = 1;
        rate_up.description = "IR +200bp";
        scenarios.push_back(rate_up);

        StressScenario rate_down;
        rate_down.name = "Rate Down";
        rate_down.rate_shocks["IR"] = -200;
        rate_down.probability = 1;
        rate_down.description = "IR -200bp";
        scenarios.push_back(rate_down);

        StressScenario credit_spread;
        credit_spread.name = "Credit Spread Widening";
        credit_spread.rate_shocks["CreditSpread"] = 150;
        credit_spread.probability = 1;
        credit_spread.description = "Credit spread +150bp";
        scenarios.push_back(credit_spread);

        StressScenario basis;
        basis.name = "Basis Risk";
        basis.rate_shocks["Basis"] = 100;
        basis.probability = 1;
        basis.description = "Basis +100bp";
        scenarios.push_back(basis);

        StressScenario vol_up;
        vol_up.name = "Volatility Up";
        vol_up.vol_shocks["Vol"] = 0.50;
        vol_up.probability = 1;
        vol_up.description = "Vol +50%";
        scenarios.push_back(vol_up);

        return scenarios;
    }

    static HistoricalCrisis crisis_2008() {
        HistoricalCrisis c;
        c.name = "2008 Financial Crisis";
        c.start_date = "2008-09-01";
        c.end_date = "2008-12-31";
        c.asset_returns["Equity"] = -0.40;
        c.asset_returns["Bond"] = 0.10;
        c.rate_changes["IR"] = -0.03;
        c.vol_changes["Vol"] = 0.80;
        return c;
    }

    static HistoricalCrisis crisis_2020() {
        HistoricalCrisis c;
        c.name = "2020 COVID-19 Shock";
        c.start_date = "2020-02-15";
        c.end_date = "2020-03-31";
        c.asset_returns["Equity"] = -0.35;
        c.asset_returns["Bond"] = 0.05;
        c.rate_changes["IR"] = -0.015;
        c.vol_changes["Vol"] = 0.60;
        return c;
    }

    Real scenario_weighted_es(const std::vector<StressScenario>& scenarios, Real confidence) {
        std::vector<Real> pnls;
        pnls.reserve(scenarios.size());

        for (const auto& scenario : scenarios) {
            auto pnl = run(scenario);
            pnls.push_back(pnl.total_pnl);
        }

        std::sort(pnls.begin(), pnls.end());

        Size n = pnls.size();
        Size es_idx = static_cast<Size>(std::floor(static_cast<Real>(n) * (Real(1) - confidence)));
        if (es_idx >= n) es_idx = n - 1;

        Real sum = Real(0);
        Size count = 0;
        for (Size i = 0; i <= es_idx; ++i) {
            sum += pnls[i];
            ++count;
        }

        return (count > 0) ? sum / static_cast<Real>(count) : Real(0);
    }

private:
    PortfolioValueFn value_fn_;
    std::map<std::string, Real> current_factors_;

    std::map<std::string, Real> apply_shocks(const StressScenario& scenario) const {
        auto factors = current_factors_;
        for (const auto& [asset, shock] : scenario.spot_shocks) {
            if (factors.find(asset) != factors.end()) {
                factors[asset] = current_factors_.at(asset) * (Real(1) + shock);
            } else {
                factors[asset] = shock;
            }
        }
        for (const auto& [rate, shock] : scenario.rate_shocks) {
            if (factors.find(rate) != factors.end()) {
                factors[rate] = current_factors_.at(rate) + shock / Real(10000);
            } else {
                factors[rate] = shock;
            }
        }
        for (const auto& [vol, shock] : scenario.vol_shocks) {
            if (factors.find(vol) != factors.end()) {
                factors[vol] = current_factors_.at(vol) * (Real(1) + shock);
            } else {
                factors[vol] = shock;
            }
        }
        return factors;
    }
};

} // namespace v1
} // namespace cpphub
