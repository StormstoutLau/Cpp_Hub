# Cpp_Hub Python Bindings

Python bindings for the [Cpp_Hub](https://github.com/) quantitative finance library.

## Quick Start

```bash
pip install nanobind scikit-build-core cmake ninja
pip install .
```

## Example

```python
import cpphub

# Black-Scholes-Merton European call price
price = cpphub.bsm_price(S=100, K=100, T=1.0, r=0.05, q=0.0, sigma=0.20, is_call=True)
print(f"Price: {price:.6f}")

# Full Greeks (returned as dict)
greeks = cpphub.bsm_greeks(100, 100, 1.0, 0.05, 0.0, 0.20, True)
print(f"Delta: {greeks['delta']:.6f}")
print(f"Gamma: {greeks['gamma']:.6f}")
print(f"Vega:  {greeks['vega']:.6f}")

# Historical VaR
pnl = [0.001, -0.002, 0.0005, -0.003, 0.002, ...]
var = cpphub.historical_var(pnl, confidence=0.99, interpolation="linear")
```

## Modules

- **core**: `normal_pdf`, `normal_cdf`, `inv_normal_cdf`
- **bsm**: Black-Scholes-Merton European option pricing and Greeks
- **heston**: Heston characteristic function
- **greeks**: Analytic + AAD Greeks (reverse-mode autodiff)
- **var**: Historical / Parametric (Normal) / Weighted VaR

## Testing

```bash
pip install pytest numpy scipy
pytest python/tests -v
```
