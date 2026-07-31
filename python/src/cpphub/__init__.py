"""Cpp_Hub: Quantitative finance library.

Modules:
    - core: normal_pdf/cdf/inv_cdf
    - bsm: Black-Scholes-Merton pricing and Greeks
    - heston: Heston characteristic function
    - greeks: Analytic + AAD Greeks
    - var: Historical/Parametric/Weighted VaR
    - rng: Philox4x64 counter-based RNG
    - cos: COS method pricing engine (GBM/Heston/Bates/VG)
    - bates: Bates model (CF + path simulation)
    - vg: Variance Gamma model (CF + path simulation)
    - cev: CEV model (analytic pricing + path simulation)
    - sabr: SABR Hagan implied volatility
    - calibration: Heston/SABR calibrators + BSM implied vol inversion

Example:
    >>> import cpphub
    >>> cpphub.bsm_price(100, 100, 1, 0.05, 0, 0.2, True)
    10.450583572185097
    >>> greeks = cpphub.bsm_greeks(100, 100, 1, 0.05, 0, 0.2, True)
    >>> greeks["delta"]
    0.6368306511756272
    >>> cpphub.cos_call_gbm(100, 100, 1, 0.05, 0, 0.2)  # COS method
    10.450583572185097
"""

from ._core import (
    # Core
    normal_pdf,
    normal_cdf,
    inv_normal_cdf,
    # BSM
    bsm_price,
    bsm_delta,
    bsm_gamma,
    bsm_vega,
    bsm_theta,
    bsm_rho,
    bsm_greeks,
    # AAD Greeks
    aad_greeks_bsm,
    # Heston
    heston_cf,
    # VaR
    historical_var,
    parametric_var_normal,
    weighted_var,
    # RNG
    Philox4x64,
    # COS method
    COSConfig,
    COSEngine,
    make_gbm_cf,
    make_heston_cf,
    make_bates_cf,
    make_vg_cf,
    make_vg_cf_direct,
    cos_call_gbm,
    cos_put_gbm,
    cos_call_heston,
    cos_put_heston,
    cos_call_bates,
    cos_put_bates,
    cos_call_vg,
    cos_put_vg,
    # Bates model
    BatesCFParams,
    BatesParams,
    BatesProcess,
    HestonScheme,
    bates_characteristic_function,
    bates_jump_compensation,
    merton_jump_cf,
    # Variance Gamma
    VGParams,
    VarianceGammaProcess,
    vg_omega,
    vg_characteristic_function,
    vg_cumulant_mean,
    vg_cumulant_variance,
    vg_cumulant_skewness,
    vg_cumulant_kurtosis_excess,
    # CEV
    CEVParams,
    CEVProcess,
    CEVScheme,
    cev_call_price,
    cev_put_price,
    validate_cev_params,
    # SABR
    SABRParams,
    sabr_implied_vol_hagan,
    # Calibration
    HestonCalibParams,
    MarketQuote,
    CalibConfig,
    CalibrationResult,
    WeightingScheme,
    HestonCalibrator,
    SABRCalibrator,
    bsm_implied_vol,
    # Metadata
    __version__,
    __author__,
)

__all__ = [
    # Core
    "normal_pdf", "normal_cdf", "inv_normal_cdf",
    # BSM
    "bsm_price", "bsm_delta", "bsm_gamma", "bsm_vega",
    "bsm_theta", "bsm_rho", "bsm_greeks",
    # AAD Greeks
    "aad_greeks_bsm",
    # Heston
    "heston_cf",
    # VaR
    "historical_var", "parametric_var_normal", "weighted_var",
    # RNG
    "Philox4x64",
    # COS method
    "COSConfig", "COSEngine",
    "make_gbm_cf", "make_heston_cf", "make_bates_cf",
    "make_vg_cf", "make_vg_cf_direct",
    "cos_call_gbm", "cos_put_gbm",
    "cos_call_heston", "cos_put_heston",
    "cos_call_bates", "cos_put_bates",
    "cos_call_vg", "cos_put_vg",
    # Bates
    "BatesCFParams", "BatesParams", "BatesProcess",
    "HestonScheme",
    "bates_characteristic_function",
    "bates_jump_compensation", "merton_jump_cf",
    # Variance Gamma
    "VGParams", "VarianceGammaProcess",
    "vg_omega", "vg_characteristic_function",
    "vg_cumulant_mean", "vg_cumulant_variance",
    "vg_cumulant_skewness", "vg_cumulant_kurtosis_excess",
    # CEV
    "CEVParams", "CEVProcess", "CEVScheme",
    "cev_call_price", "cev_put_price", "validate_cev_params",
    # SABR
    "SABRParams", "sabr_implied_vol_hagan",
    # Calibration
    "HestonCalibParams", "MarketQuote", "CalibConfig",
    "CalibrationResult", "WeightingScheme",
    "HestonCalibrator", "SABRCalibrator",
    "bsm_implied_vol",
    # Metadata
    "__version__", "__author__",
]
