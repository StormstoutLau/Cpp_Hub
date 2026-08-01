"""Cpp_Hub Python 绑定 v1.4 冒烟测试。

覆盖第二批绑定:
  - 模型过程: Heston / HestonQE / RoughHeston / RoughBergomi / SABR
  - Risk 模块: ExpectedShortfall / KupiecPOF / ChristoffersenIID / BaselTrafficLight
  - GreeksFactory: Auto 分派 (Vanilla→Analytic, Digital→LR)

运行: py -3.12 -m pytest python/tests/test_cpphub_v14.py -v
"""

import math
import pytest
import numpy as np
from scipy import stats as st

import cpphub


# ===========================================================================
# 1. 模型过程: Heston / HestonQE / RoughHeston / RoughBergomi / SABR
# ===========================================================================

class TestHestonProcess:
    def _make_params(self):
        p = cpphub.HestonProcessParams()
        p.S0 = 100.0
        p.v0 = 0.04
        p.kappa = 1.5
        p.theta = 0.04
        p.sigma = 0.3
        p.rho = -0.5
        p.r = 0.03
        p.q = 0.0
        return p

    def test_construct_and_metadata(self):
        p = self._make_params()
        proc = cpphub.HestonProcess(p)
        assert proc.dimension() == 2
        assert proc.spot() == pytest.approx(100.0)

    def test_characteristic_function_at_zero(self):
        p = self._make_params()
        proc = cpphub.HestonProcess(p)
        phi = proc.characteristic_function(complex(0.0, 0.0), 1.0)
        assert abs(phi.real - 1.0) < 1e-12
        assert abs(phi.imag) < 1e-12

    def test_cf_unit_modulus(self):
        p = self._make_params()
        proc = cpphub.HestonProcess(p)
        for u_real in [0.1, 0.5, 1.0, 2.0]:
            phi = proc.characteristic_function(complex(u_real, 0.0), 1.0)
            assert abs(phi) <= 1.0 + 1e-12

    def test_generate_path_length(self):
        p = self._make_params()
        proc = cpphub.HestonProcess(p)
        rng = cpphub.Philox4x64(seed=42)
        path = proc.generate_path(T=1.0, n_steps=252, path=[], rng=rng)
        assert len(path) == 253
        assert path[0] == pytest.approx(100.0, rel=1e-12)

    def test_qe_subclass(self):
        p = self._make_params()
        proc = cpphub.HestonQEProcess(p)
        rng = cpphub.Philox4x64(seed=7)
        path = proc.generate_path(T=0.5, n_steps=100, path=[], rng=rng)
        assert len(path) == 101
        assert path[0] == pytest.approx(100.0, rel=1e-12)


class TestRoughHestonProcess:
    def _make_params(self):
        p = cpphub.RoughHestonParams()
        p.H = 0.1
        p.kappa = 0.5
        p.theta = 0.02
        p.sigma = 0.3
        p.rho = -0.6
        p.v0 = 0.02
        p.S0 = 100.0
        p.r = 0.03
        p.q = 0.0
        return p

    def test_construct(self):
        p = self._make_params()
        proc = cpphub.RoughHestonProcess(p, T=0.25, n_steps=64)
        assert proc.dimension() == 2
        assert proc.spot() == pytest.approx(100.0)
        assert proc.alpha() == pytest.approx(0.6)  # H + 0.5

    def test_validate_params(self):
        p = self._make_params()
        cpphub.validate_rough_heston_params(p)  # 不抛异常

    def test_kernel_shape(self):
        K = cpphub.rough_heston_kernel(T=0.25, n_steps=8, alpha=0.6)
        assert len(K) == 8
        for row in K:
            assert len(row) == 8

    def test_generate_path(self):
        p = self._make_params()
        proc = cpphub.RoughHestonProcess(p, T=0.25, n_steps=64)
        rng = cpphub.Philox4x64(seed=11)
        path = proc.generate_path(T=0.25, n_steps=64, path=[], rng=rng)
        assert len(path) == 65
        assert path[0] == pytest.approx(100.0, rel=1e-12)

    def test_generate_path_sv(self):
        p = self._make_params()
        proc = cpphub.RoughHestonProcess(p, T=0.25, n_steps=64)
        rng = cpphub.Philox4x64(seed=13)
        result = proc.generate_path_sv(rng=rng)
        assert len(result) == 2  # [spot_path, var_path]
        assert len(result[0]) == 65
        assert len(result[1]) == 65


class TestRoughBergomiProcess:
    def _make_params(self):
        p = cpphub.RoughBergomiParams()
        p.H = 0.07
        p.eta = 1.9
        p.rho = -0.9
        p.xi0 = 0.235 ** 2
        p.S0 = 100.0
        p.r = 0.03
        p.q = 0.0
        return p

    def test_construct(self):
        p = self._make_params()
        proc = cpphub.RoughBergomiProcess(p)
        assert proc.dimension() == 2
        assert proc.spot() == pytest.approx(100.0)

    def test_validate(self):
        p = self._make_params()
        cpphub.validate_rough_bergomi_params(p)

    def test_generate_path(self):
        p = self._make_params()
        proc = cpphub.RoughBergomiProcess(p)
        rng = cpphub.Philox4x64(seed=23)
        path = proc.generate_path(T=0.25, n_steps=64, path=[], rng=rng)
        assert len(path) == 65
        assert path[0] == pytest.approx(100.0, rel=1e-12)


class TestSABRProcess:
    def test_construct_and_path(self):
        params = cpphub.SABRParams()
        params.alpha = 0.3
        params.beta = 0.5
        params.rho = -0.3
        params.nu = 0.4
        proc = cpphub.SABRProcess(params=params, F0=100.0,
                                  r=0.0, q=0.0,
                                  scheme=cpphub.SABRScheme.EulerAbsorbing)
        assert proc.dimension() == 2
        assert proc.forward0() == pytest.approx(100.0)
        rng = cpphub.Philox4x64(seed=31)
        path = proc.generate_path(T=1.0, n_steps=128, path=[], rng=rng)
        assert len(path) == 129
        assert path[0] == pytest.approx(100.0, rel=1e-12)


# ===========================================================================
# 2. Risk 模块: ExpectedShortfall / Kupiec / Basel
# ===========================================================================

class TestExpectedShortfall:
    def test_normal_es_vs_scipy(self):
        es = cpphub.ExpectedShortfall()
        mean, sigma, c = 0.0, 0.02, 0.99
        result = es.normal_es(mean, sigma, c)
        # Normal ES = -μ + σ·φ(z)/(1-c), z = Φ⁻¹(c) (loss convention)
        z = st.norm.ppf(c)
        expected = mean + sigma * st.norm.pdf(z) / (1 - c)
        assert abs(result - expected) < 1e-12

    def test_from_losses(self):
        es = cpphub.ExpectedShortfall()
        rng = np.random.default_rng(42)
        losses = rng.normal(0.0, 0.02, 5000).tolist()
        result = es.from_losses(losses, confidence=0.99)
        # ES = mean of worst 1% losses
        arr = np.sort(np.array(losses))
        n_tail = int(np.ceil(0.01 * len(arr)))
        expected = arr[-n_tail:].mean()
        assert abs(result - expected) < 1e-12

    def test_student_t_es_finite(self):
        es = cpphub.ExpectedShortfall()
        result = es.student_t_es(mean=0.0, sigma=0.02, dof=5.0,
                                 confidence=0.99)
        assert math.isfinite(result)
        assert result > 0

    def test_cornish_fisher_es_finite(self):
        es = cpphub.ExpectedShortfall()
        result = es.cornish_fisher_es(mean=0.0, sigma=0.02, skew=-0.5,
                                      kurt=4.0, confidence=0.99)
        assert math.isfinite(result)
        assert result > 0


class TestKupiecPOF:
    def test_well_calibrated_model(self):
        # 250 obs, 1% expected = 2.5 violations; p should be high
        result = cpphub.KupiecPOF.test(n_violations=2, n_observations=250,
                                       confidence=0.99)
        assert result.n_violations == 2
        assert result.n_observations == 250
        assert result.violation_rate == pytest.approx(2/250)
        assert result.expected_rate == pytest.approx(0.01)
        assert not result.reject_null  # 应不拒绝原假设

    def test_underestimated_risk(self):
        # 10 violations / 250 obs (4%) at 99% → 应拒绝
        result = cpphub.KupiecPOF.test(n_violations=10, n_observations=250,
                                       confidence=0.99)
        assert result.reject_null

    def test_test_series(self):
        rng = np.random.default_rng(7)
        losses = rng.normal(0.0, 0.02, 500).tolist()
        var_series = [0.05] * 500  # 固定 VaR
        result = cpphub.KupiecPOF.test_series(var_series, losses,
                                              confidence=0.99)
        assert 0 <= result.violation_rate


class TestBaselTrafficLight:
    def test_green_zone(self):
        result = cpphub.BaselTrafficLight.assess(n_violations=3)
        assert result.zone == cpphub.BaselZone.Green
        assert result.capital_multiplier == pytest.approx(3.0)

    def test_yellow_zone(self):
        result = cpphub.BaselTrafficLight.assess(n_violations=7)
        assert result.zone == cpphub.BaselZone.Yellow

    def test_red_zone(self):
        result = cpphub.BaselTrafficLight.assess(n_violations=12)
        assert result.zone == cpphub.BaselZone.Red
        assert result.capital_multiplier == pytest.approx(4.0)


# ===========================================================================
# 3. GreeksFactory: Auto 分派
# ===========================================================================

class TestGreeksFactory:
    def test_vanilla_call_auto_dispatches_to_analytic(self):
        g = cpphub.GreeksFactory.compute_bsm(
            S=100, K=100, T=1.0, r=0.05, q=0.0, sigma=0.20,
            payoff=cpphub.PayoffType.VanillaCall,
            method=cpphub.GreeksMethod.Auto)
        # 与解析 BSM 对照
        analytic = cpphub.bsm_greeks(100, 100, 1.0, 0.05, 0.0, 0.20, True)
        assert abs(g.price - analytic["price"]) < 1e-10
        assert abs(g.delta - analytic["delta"]) < 1e-10
        assert abs(g.gamma - analytic["gamma"]) < 1e-10
        assert g.method_used == cpphub.GreeksMethod.Analytic

    def test_vanilla_put_analytic(self):
        g = cpphub.GreeksFactory.compute_bsm(
            S=100, K=100, T=1.0, r=0.05, q=0.0, sigma=0.20,
            payoff=cpphub.PayoffType.VanillaPut,
            method=cpphub.GreeksMethod.Analytic)
        analytic = cpphub.bsm_greeks(100, 100, 1.0, 0.05, 0.0, 0.20, False)
        assert abs(g.price - analytic["price"]) < 1e-10
        assert abs(g.delta - analytic["delta"]) < 1e-10

    def test_digital_call_auto_dispatches_to_lr(self):
        g = cpphub.GreeksFactory.compute_bsm(
            S=100, K=100, T=1.0, r=0.05, q=0.0, sigma=0.20,
            payoff=cpphub.PayoffType.DigitalCall,
            method=cpphub.GreeksMethod.Auto,
            n_paths=200000, seed=42)
        # 数字看涨期权理论价格 = e^{-rT} N(d2)
        d2 = (math.log(100/100) + (0.05 - 0.5*0.20**2)*1.0) / (0.20*math.sqrt(1.0))
        expected_price = math.exp(-0.05) * st.norm.cdf(d2)
        assert abs(g.price - expected_price) < 0.01  # MC 误差容忍
        # 不连续 payoff 应分派到 LR
        assert g.method_used == cpphub.GreeksMethod.LR

    def test_fd_method_universal(self):
        g = cpphub.GreeksFactory.compute_bsm(
            S=100, K=100, T=1.0, r=0.05, q=0.0, sigma=0.20,
            payoff=cpphub.PayoffType.VanillaCall,
            method=cpphub.GreeksMethod.FD)
        assert g.method_used == cpphub.GreeksMethod.FD
        analytic = cpphub.bsm_price(100, 100, 1.0, 0.05, 0.0, 0.20, True)
        assert abs(g.price - analytic) < 0.05
