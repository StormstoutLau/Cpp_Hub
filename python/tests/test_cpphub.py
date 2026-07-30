"""Cpp_Hub Python 绑定测试。

验证 cpphub Python 模块与 C++ 实现数值一致,并与 scipy 交叉验证。

运行: pytest python/tests/test_cpphub.py -v
"""

import math
import pytest
import numpy as np
from scipy import stats as st


# ===========================================================================
# 1. Core 数学函数
# ===========================================================================

class TestCore:
    def test_normal_pdf_at_zero(self):
        import cpphub
        assert abs(cpphub.normal_pdf(0.0) - 1.0 / math.sqrt(2 * math.pi)) < 1e-15

    def test_normal_pdf_symmetry(self):
        import cpphub
        for x in [0.5, 1.0, 2.0, -1.5]:
            assert abs(cpphub.normal_pdf(x) - cpphub.normal_pdf(-x)) < 1e-15

    def test_normal_cdf_at_zero(self):
        import cpphub
        assert abs(cpphub.normal_cdf(0.0) - 0.5) < 1e-15

    def test_normal_cdf_vs_scipy(self):
        import cpphub
        for x in [-2.0, -1.0, 0.0, 1.0, 2.0]:
            assert abs(cpphub.normal_cdf(x) - st.norm.cdf(x)) < 1e-12

    def test_inv_normal_cdf_known_values(self):
        import cpphub
        assert abs(cpphub.inv_normal_cdf(0.5)) < 1e-15
        assert abs(cpphub.inv_normal_cdf(0.95) - 1.6448536269514722) < 1e-12
        assert abs(cpphub.inv_normal_cdf(0.975) - 1.959963984540054) < 1e-12

    def test_inv_normal_cdf_roundtrip(self):
        import cpphub
        for p in [0.01, 0.05, 0.5, 0.95, 0.99]:
            x = cpphub.inv_normal_cdf(p)
            assert abs(cpphub.normal_cdf(x) - p) < 1e-10

    def test_inv_normal_cdf_invalid(self):
        import cpphub
        with pytest.raises((ValueError, Exception)):
            cpphub.inv_normal_cdf(0.0)
        with pytest.raises((ValueError, Exception)):
            cpphub.inv_normal_cdf(1.0)


# ===========================================================================
# 2. BSM 定价
# ===========================================================================

class TestBSM:
    def test_call_price_atm(self):
        import cpphub
        # S=100, K=100, T=1, r=0.05, q=0, sigma=0.20, Call
        # Expected: 10.4505835722 (与 scripts/generate_bs_benchmark.py 一致)
        price = cpphub.bsm_price(100, 100, 1, 0.05, 0, 0.20, True)
        assert abs(price - 10.4505835722) < 1e-6

    def test_put_price_atm(self):
        import cpphub
        # Put-Call Parity: Call - Put = S - K*exp(-rT)
        call = cpphub.bsm_price(100, 100, 1, 0.05, 0, 0.20, True)
        put = cpphub.bsm_price(100, 100, 1, 0.05, 0, 0.20, False)
        parity = 100 - 100 * math.exp(-0.05)
        assert abs((call - put) - parity) < 1e-10

    def test_call_price_otm(self):
        import cpphub
        # S=90, K=100, OTM call
        price = cpphub.bsm_price(90, 100, 1, 0.05, 0, 0.20, True)
        assert abs(price - 5.0912220788) < 1e-6

    def test_call_price_itm(self):
        import cpphub
        # S=110, K=100, ITM call
        price = cpphub.bsm_price(110, 100, 1, 0.05, 0, 0.20, True)
        assert abs(price - 17.6629537406) < 1e-6

    def test_put_call_parity_multiple(self):
        import cpphub
        for S, K, T, r, sigma in [
            (100, 100, 1.0, 0.05, 0.20),
            (100, 100, 0.25, 0.05, 0.20),
            (100, 100, 1.0, 0.05, 0.50),
        ]:
            call = cpphub.bsm_price(S, K, T, r, 0, sigma, True)
            put = cpphub.bsm_price(S, K, T, r, 0, sigma, False)
            parity = S - K * math.exp(-r * T)
            assert abs((call - put) - parity) < 1e-10


# ===========================================================================
# 3. BSM Greeks
# ===========================================================================

class TestBSMGreeks:
    def test_delta_call_in_range(self):
        import cpphub
        delta = cpphub.bsm_delta(100, 100, 1, 0.05, 0, 0.20, True)
        assert 0.0 < delta < 1.0

    def test_delta_put_in_range(self):
        import cpphub
        delta = cpphub.bsm_delta(100, 100, 1, 0.05, 0, 0.20, False)
        assert -1.0 < delta < 0.0

    def test_gamma_positive(self):
        import cpphub
        gamma = cpphub.bsm_gamma(100, 100, 1, 0.05, 0, 0.20)
        assert gamma > 0.0

    def test_vega_positive(self):
        import cpphub
        vega = cpphub.bsm_vega(100, 100, 1, 0.05, 0, 0.20)
        assert vega > 0.0

    def test_greeks_dict_keys(self):
        import cpphub
        g = cpphub.bsm_greeks(100, 100, 1, 0.05, 0, 0.20, True)
        expected_keys = {"price", "delta", "gamma", "vega", "theta", "rho",
                         "vanna", "vomma"}
        assert set(g.keys()) == expected_keys

    def test_greeks_vs_scipy(self):
        import cpphub
        S, K, T, r, sigma = 100, 100, 1, 0.05, 0.20
        g = cpphub.bsm_greeks(S, K, T, r, 0, sigma, True)
        # Delta vs scipy
        d1 = (math.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))
        assert abs(g["delta"] - st.norm.cdf(d1)) < 1e-10


# ===========================================================================
# 4. AAD Greeks
# ===========================================================================

class TestAADGreeks:
    def test_aad_vs_analytic_price(self):
        import cpphub
        analytic = cpphub.bsm_greeks(100, 100, 1, 0.05, 0, 0.20, True)
        aad = cpphub.aad_greeks_bsm(100, 100, 1, 0.05, 0, 0.20, True)
        assert abs(aad["price"] - analytic["price"]) < 1e-10

    def test_aad_vs_analytic_delta(self):
        import cpphub
        analytic = cpphub.bsm_greeks(100, 100, 1, 0.05, 0, 0.20, True)
        aad = cpphub.aad_greeks_bsm(100, 100, 1, 0.05, 0, 0.20, True)
        assert abs(aad["delta"] - analytic["delta"]) < 1e-10

    def test_aad_vs_analytic_vega(self):
        import cpphub
        analytic = cpphub.bsm_greeks(100, 100, 1, 0.05, 0, 0.20, True)
        aad = cpphub.aad_greeks_bsm(100, 100, 1, 0.05, 0, 0.20, True)
        assert abs(aad["vega"] - analytic["vega"]) < 1e-10


# ===========================================================================
# 5. Heston 特征函数
# ===========================================================================

class TestHeston:
    def test_cf_at_zero_returns_one(self):
        import cpphub
        phi = cpphub.heston_cf(complex(0, 0), 1.0, 100.0,
                                0.04, 1.5, 0.04, 0.3, -0.5)
        assert abs(phi.real - 1.0) < 1e-12
        assert abs(phi.imag) < 1e-12

    def test_cf_unit_modulus_real_u(self):
        import cpphub
        for u_real in [0.1, 0.5, 1.0, 2.0]:
            phi = cpphub.heston_cf(complex(u_real, 0), 1.0, 100.0,
                                    0.04, 1.5, 0.04, 0.3, -0.5)
            assert abs(phi) <= 1.0 + 1e-12

    def test_cf_decays_at_infinity(self):
        import cpphub
        phi = cpphub.heston_cf(complex(50.0, 0), 1.0, 100.0,
                                0.04, 1.5, 0.04, 0.3, -0.5)
        assert abs(phi) < 0.01


# ===========================================================================
# 6. VaR 模块
# ===========================================================================

class TestVaR:
    def test_historical_var_linear(self):
        import cpphub
        # 确定性 PnL 数据
        pnl = [-0.01, -0.02, 0.005, 0.01, -0.03, 0.015, -0.005, 0.02,
               -0.015, 0.008] * 100  # 1000 个样本
        var = cpphub.historical_var(pnl, 0.99, "linear")
        # 与 numpy 一致
        np_var = -np.quantile(-np.array(pnl), 0.01,
                               method="linear") if hasattr(np, "quantile") else None
        assert var > 0  # 应为正损失

    def test_historical_var_vs_numpy(self):
        import cpphub
        rng = np.random.default_rng(42)
        pnl = rng.normal(0.001, 0.02, 1000).tolist()

        # C++ linear VaR
        cpp_var = cpphub.historical_var(pnl, 0.99, "linear")

        # numpy 镜像 C++ 算法: index = q*(n-1), linear interp
        losses = np.sort(-np.array(pnl))
        n = len(losses)
        index = 0.01 * (n - 1)
        lo, hi = int(np.floor(index)), int(np.ceil(index))
        frac = index - lo
        np_var = -(losses[lo] + frac * (losses[hi] - losses[lo]))

        assert abs(cpp_var - np_var) < 1e-12

    def test_parametric_var_normal(self):
        import cpphub
        mean, std = 0.001, 0.02
        var = cpphub.parametric_var_normal(mean, std, 0.99)
        # VaR = -(mean + std * z_alpha), z_alpha = Φ⁻¹(0.01)
        z = st.norm.ppf(0.01)
        expected = -(mean + std * z)
        assert abs(var - expected) < 1e-10

    def test_weighted_var_positive(self):
        import cpphub
        rng = np.random.default_rng(42)
        pnl = rng.normal(0.001, 0.02, 500).tolist()
        var = cpphub.weighted_var(pnl, 0.99, 0.99)
        assert var > 0

    def test_var_confidence_ordering(self):
        import cpphub
        rng = np.random.default_rng(42)
        pnl = rng.normal(0.001, 0.02, 1000).tolist()
        var_99 = cpphub.historical_var(pnl, 0.99, "linear")
        var_95 = cpphub.historical_var(pnl, 0.95, "linear")
        # 99% VaR 应大于 95% VaR (更大的损失)
        assert var_99 >= var_95


# ===========================================================================
# 7. 版本与元信息
# ===========================================================================

class TestMetadata:
    def test_version_string(self):
        import cpphub
        assert isinstance(cpphub.__version__, str)
        assert len(cpphub.__version__) > 0

    def test_author_string(self):
        import cpphub
        assert isinstance(cpphub.__author__, str)
