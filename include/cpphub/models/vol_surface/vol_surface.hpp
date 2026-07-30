#pragma once
// SOURCE: PHASE3_SPEC §4.2 - Volatility surface (interpolation + SVI parameterization)
// STUB: Interface declarations only, no implementation yet.
//       VolSurface class is planned for v1.1+.
//       For Dupire local volatility, use DupireLocalVol (dupire_local_vol.hpp)
//       which has a self-contained IV grid path and does not depend on this class.
//       For SVI slice calibration, use SVI::calibrate() in svi.hpp.
#include "cpphub/core/types.hpp"
#include "cpphub/models/vol_surface/svi.hpp"
#include <vector>
#include <map>
#include <string>
#include <memory>

namespace cpphub {
inline namespace v1 {

enum class InterpolationMethod {
    Bilinear,
    Bicubic,
    CubicSpline,
    SVI
};

class VolSurface {
public:
    VolSurface();
    VolSurface(const std::vector<Real>& strikes,
                const std::vector<Real>& maturities,
                const std::vector<std::vector<Real>>& implied_vols);

    Real implied_vol(Real K, Real T, InterpolationMethod method = InterpolationMethod::Bilinear) const;
    Real total_variance(Real K, Real T) const;
    Real log_moneyness(Real K, Real T) const;

    // 微分 (用于 Dupire)
    Real dC_dT(Real K, Real T) const;
    Real dC_dK(Real K, Real T) const;
    Real d2C_dK2(Real K, Real T) const;

    // BSM 反推
    Real call_price(Real K, Real T, Real S, Real r, Real q) const;
    Real put_price(Real K, Real T, Real S, Real r, Real q) const;

    // 数据访问
    const std::vector<Real>& strikes() const { return strikes_; }
    const std::vector<Real>& maturities() const { return maturities_; }
    const std::vector<std::vector<Real>>& vols() const { return implied_vols_; }

    // SVI 参数化拟合
    void fit_svi(Real T);
    const SVI* svi_slice(Real T) const;

    // 套利检查
    bool check_calendar_arbitrage() const;
    bool check_butterfly_arbitrage() const;

private:
    std::vector<Real> strikes_;
    std::vector<Real> maturities_;
    std::vector<std::vector<Real>> implied_vols_;
    Real forward_ = 100.0;
    Real r_ = 0.0;
    Real q_ = 0.0;

    // 缓存的 SVI slice
    std::map<Real, std::unique_ptr<SVI>> svi_slices_;

    Real bilinear_interp(Real K, Real T) const;
    Real cubic_spline_interp(Real K, Real T) const;
    void build_splines() const;
};

}  // namespace v1
}  // namespace cpphub
