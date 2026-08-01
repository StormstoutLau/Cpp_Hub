"""Cpp_Hub: Quantitative finance library.

Modules:
    - core: normal_pdf/cdf/inv_cdf
    - bsm: Black-Scholes-Merton pricing and Greeks
    - heston: Heston characteristic function
    - greeks: Analytic + AAD Greeks + unified GreeksFactory
    - var: Historical/Parametric/Weighted/MC VaR + ES + Backtesting
    - rng: Philox4x64 counter-based RNG
    - cos: COS method pricing engine (GBM/Heston/Bates/VG/CGMY/Kou/NIG)
    - bates: Bates model (CF + path simulation)
    - vg: Variance Gamma model (CF + path simulation)
    - cev: CEV model (analytic pricing + path simulation)
    - sabr: SABR Hagan implied volatility + SABR process
    - calibration: Heston/SABR calibrators + BSM implied vol inversion
    - payoff: PayOff class hierarchy (Call/Put/Digital/DoubleDigital)
    - pde: PDE finite difference engine (Explicit/Implicit/Crank-Nicolson)
    - tree: Binomial/Trinomial tree engines
    - mc: Monte Carlo engine (European + path-dependent Asian/Barrier/Lookback)
    - models: Heston/HestonQE/RoughHeston/RoughBergomi/SABR stochastic processes
    - risk: ES / MCVaR / Kupiec / Christoffersen / Basel traffic light

Example:
    >>> import cpphub
    >>> cpphub.bsm_price(100, 100, 1, 0.05, 0, 0.2, True)
    10.450583572185097
    >>> greeks = cpphub.bsm_greeks(100, 100, 1, 0.05, 0, 0.2, True)
    >>> greeks["delta"]
    0.6368306511756272
    >>> cpphub.cos_call_gbm(100, 100, 1, 0.05, 0, 0.2)  # COS method
    10.450583572185097
    >>> # Unified Greeks via factory (auto-dispatch)
    >>> ug = cpphub.GreeksFactory.compute_bsm(100, 100, 1, 0.05, 0, 0.2,
    ...                                       cpphub.PayoffType.VanillaCall)
    >>> ug.delta
    0.6368306511756272
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
    # VaR (basic)
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
    # v1.2: Levy CF extensions (CGMY / Kou / NIG)
    make_cgmy_cf,
    make_kou_cf,
    make_nig_cf,
    cos_call_cgmy,
    cos_put_cgmy,
    cos_call_kou,
    cos_put_kou,
    cos_call_nig,
    cos_put_nig,
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
    # SABR (calibration params + Hagan IV)
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
    # v1.3: PayOff class hierarchy
    PayOff,
    CallPayOff,
    PutPayOff,
    DigitalCallPayOff,
    DigitalPutPayOff,
    DoubleDigitalPayOff,
    # v1.3: PDE engine
    FDMSchemeType,
    PDEEngineConfig,
    PDEGreeks,
    PDEEngine,
    # v1.3: Tree engines
    BinomialType,
    BinomialParams,
    BinomialTreeEngine,
    TrinomialType,
    TrinomialParams,
    TrinomialTreeEngine,
    # v1.3: MC engine + path-dependent
    MCConfig,
    MCResult,
    mc_price_european,
    mc_price_asian_arithmetic,
    mc_price_barrier_up_out_call,
    mc_price_lookback_floating_call,
    # v1.4: Model processes
    HestonProcessParams,
    HestonProcess,
    HestonQEProcess,
    RoughHestonParams,
    RoughHestonProcess,
    validate_rough_heston_params,
    rough_heston_kernel,
    RoughBergomiParams,
    RoughBergomiProcess,
    validate_rough_bergomi_params,
    RLFbmSampler,
    SABRScheme,
    SABRProcess,
    # v1.4: Risk (ES / MCVaR / Backtesting)
    ExpectedShortfall,
    MCVarConfig,
    VaRApproximation,
    MCVaR,
    BacktestResult,
    KupiecPOF,
    ChristoffersenIID,
    BaselZone,
    BaselTrafficLightResult,
    BaselTrafficLight,
    # v1.4: GreeksFactory
    GreeksMethod,
    PayoffType,
    UnifiedGreeks,
    GreeksFactory,
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
    # VaR (basic)
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
    # v1.2: Levy CF extensions
    "make_cgmy_cf", "make_kou_cf", "make_nig_cf",
    "cos_call_cgmy", "cos_put_cgmy",
    "cos_call_kou", "cos_put_kou",
    "cos_call_nig", "cos_put_nig",
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
    # v1.3: PayOff
    "PayOff", "CallPayOff", "PutPayOff",
    "DigitalCallPayOff", "DigitalPutPayOff", "DoubleDigitalPayOff",
    # v1.3: PDE
    "FDMSchemeType", "PDEEngineConfig", "PDEGreeks", "PDEEngine",
    # v1.3: Tree
    "BinomialType", "BinomialParams", "BinomialTreeEngine",
    "TrinomialType", "TrinomialParams", "TrinomialTreeEngine",
    # v1.3: MC
    "MCConfig", "MCResult",
    "mc_price_european", "mc_price_asian_arithmetic",
    "mc_price_barrier_up_out_call", "mc_price_lookback_floating_call",
    # v1.4: Model processes
    "HestonProcessParams", "HestonProcess", "HestonQEProcess",
    "RoughHestonParams", "RoughHestonProcess",
    "validate_rough_heston_params", "rough_heston_kernel",
    "RoughBergomiParams", "RoughBergomiProcess",
    "validate_rough_bergomi_params", "RLFbmSampler",
    "SABRScheme", "SABRProcess",
    # v1.4: Risk
    "ExpectedShortfall",
    "MCVarConfig", "VaRApproximation", "MCVaR",
    "BacktestResult", "KupiecPOF", "ChristoffersenIID",
    "BaselZone", "BaselTrafficLightResult", "BaselTrafficLight",
    # v1.4: GreeksFactory
    "GreeksMethod", "PayoffType", "UnifiedGreeks", "GreeksFactory",
    # Metadata
    "__version__", "__author__",
]
