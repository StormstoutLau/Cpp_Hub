"""Cpp_Hub: Quantitative finance library.

Modules:
    - core: normal_pdf/cdf/inv_cdf
    - bsm: Black-Scholes-Merton pricing and Greeks
    - heston: Heston characteristic function
    - greeks: Analytic + AAD Greeks
    - var: Historical/Parametric/Weighted VaR

Example:
    >>> import cpphub
    >>> cpphub.bsm_price(100, 100, 1, 0.05, 0, 0.2, True)
    10.450583572185097
    >>> greeks = cpphub.bsm_greeks(100, 100, 1, 0.05, 0, 0.2, True)
    >>> greeks["delta"]
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
    # VaR
    historical_var,
    parametric_var_normal,
    weighted_var,
    # Metadata
    __version__,
    __author__,
)

__all__ = [
    "normal_pdf", "normal_cdf", "inv_normal_cdf",
    "bsm_price", "bsm_delta", "bsm_gamma", "bsm_vega",
    "bsm_theta", "bsm_rho", "bsm_greeks",
    "aad_greeks_bsm",
    "heston_cf",
    "historical_var", "parametric_var_normal", "weighted_var",
    "__version__", "__author__",
]
