#pragma once
// SOURCE: PHASE3_SPEC §4.1 - Model calibrators (Heston/SABR/SVI)
// STUB: Interface declarations only, no implementation yet.
//       HestonCalibrator / SABRCalibrator are planned for v1.1+.
//       SVI calibration is available via SVI::calibrate() in svi.hpp.
// NOTE: CalibrationResult and CalibConfig are defined in optimizer.hpp
#include "cpphub/core/types.hpp"
#include "cpphub/calibration/optimizer.hpp"
#include "cpphub/calibration/objective.hpp"
#include <vector>
#include <string>

namespace cpphub {
inline namespace v1 {

class Calibrator {
public:
    virtual ~Calibrator() = default;
    virtual CalibrationResult calibrate(
        const std::vector<MarketQuote>& quotes,
        const CalibConfig& cfg = CalibConfig{}) = 0;
    virtual std::string name() const = 0;
};

struct HestonParams {
    Real v0, kappa, theta, sigma_v, rho;
};

class HestonCalibrator : public Calibrator {
public:
    CalibrationResult calibrate(
        const std::vector<MarketQuote>& quotes,
        const CalibConfig& cfg = CalibConfig{}) override;
    std::string name() const override { return "HestonCalibrator"; }

    static std::vector<Bounds> default_bounds();
    static bool check_feller(const HestonParams& p);
    HestonParams extract_params(const std::vector<Real>& x) const;
};

struct SABRParams {
    Real alpha, beta, nu, rho;
};

class SABRCalibrator : public Calibrator {
public:
    CalibrationResult calibrate(
        const std::vector<MarketQuote>& quotes,
        const CalibConfig& cfg = CalibConfig{}) override;
    std::string name() const override { return "SABRCalibrator"; }

    static std::vector<Bounds> default_bounds();
    SABRParams extract_params(const std::vector<Real>& x) const;
};

}  // inline namespace v1
}  // namespace cpphub
