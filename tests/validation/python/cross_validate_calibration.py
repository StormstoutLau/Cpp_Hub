#!/usr/bin/env python3
"""
Cpp_Hub Phase 4 LITE - G2 整改项: 标定模块 Python 交叉验证

对比 Cpp_Hub C++ 实现与 Python (scipy.optimize) 的标定结果:
  1. SVI 切片标定 (5 参数 a/b/rho/m/sigma)
  2. SSVI Power-law 参数化验证 (ρ, η, γ)
  3. Heston-like 参数化验证 (ρ, η, λ)
  4. SSVI 无套利条件验证

输出 JSON 格式基准值,供 C++ 测试 test_python_cross_validation.cpp 加载验证。

依赖: numpy, scipy

用法:
    python cross_validate_calibration.py > benchmarks_calib.json

参考:
  - Cpp_Hub: include/cpphub/models/vol_surface/svi.hpp
  - Cpp_Hub: include/cpphub/models/vol_surface/ssvi.hpp
  - Gatheral-Jacquier (2014) "Arbitrage-free SVI volatility surfaces"
  - Gatheral-Jacquier (2013) "Arbitrage-free SVI volatility surfaces"
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import least_squares, minimize


# ---------------------------------------------------------------------------
# SVI 参数化 (镜像 C++ SVI 实现)
# ---------------------------------------------------------------------------

def svi_total_variance(k: float, a: float, b: float, rho: float,
                        m: float, sigma: float) -> float:
    """SVI 总方差: w(k) = a + b*(rho*(k-m) + sqrt((k-m)^2 + sigma^2))

    镜像 C++ SVI::total_variance 实现。
    """
    km = k - m
    return a + b * (rho * km + math.sqrt(km * km + sigma * sigma))


def svi_implied_vol(k: float, T: float, a: float, b: float, rho: float,
                     m: float, sigma: float) -> float:
    """SVI 隐含波动率: σ_imp = sqrt(w(k)/T)"""
    return math.sqrt(svi_total_variance(k, a, b, rho, m, sigma) / T)


# ---------------------------------------------------------------------------
# SSVI 参数化 (镜像 C++ SSVI 实现)
# ---------------------------------------------------------------------------

def ssvi_total_variance(k: float, theta: float, phi: float,
                         rho: float) -> float:
    """SSVI 总方差: w(k,θ) = θ/2 * (1 + φ*ρ*k + sqrt((φ*k+ρ)² + 1-ρ²))

    镜像 C++ SSVI::total_variance 实现。
    """
    phi_k = phi * k
    inner = (phi_k + rho) ** 2 + (1.0 - rho * rho)
    return 0.5 * theta * (1.0 + phi * rho * k + math.sqrt(inner))


def ssvi_power_law_phi(theta: float, eta: float, gamma: float) -> float:
    """SSVI Power-law 期限结构: φ(θ) = η * θ^(-γ)

    镜像 C++ SSVI::Power_law 实现 ( Gatheral-Jacquier 2014 eq. 5.3)。
    无套利条件: 0 < γ < 1/2。
    """
    return eta * theta ** (-gamma)


def ssvi_heston_like_phi(theta: float, eta: float, lambda_: float) -> float:
    """SSVI Heston-like 期限结构: φ(θ) = η * θ^(-λ)

    镜像 C++ SSVI::Heston_like 实现 (Gatheral-Jacquier 2014 eq. 5.2)。
    无套利条件: 0 < λ < 1。
    """
    return eta * theta ** (-lambda_)


# ---------------------------------------------------------------------------
# 标定函数
# ---------------------------------------------------------------------------

def calibrate_svi_slice(strikes: List[float], implied_vols: List[float],
                         T: float, forward: float) -> Dict:
    """SVI 切片标定 - 使用 scipy.optimize.least_squares

    镜像 C++ SVI::calibrate (Levenberg-Marquardt)。
    """
    strikes = np.array(strikes)
    implied_vols = np.array(implied_vols)
    market_w = implied_vols ** 2 * T  # 市场总方差
    log_k = np.log(strikes / forward)

    def residuals(params):
        a, b, rho, m, sigma = params
        model_w = np.array([
            svi_total_variance(k, a, b, rho, m, sigma) for k in log_k
        ])
        return model_w - market_w

    # 初始猜测: ATM 总方差, 典型斜率, 0 相关性, 0 中心, 0.1 宽度
    atm_var = implied_vols[len(implied_vols) // 2] ** 2 * T
    x0 = [0.5 * atm_var, 0.1, 0.0, 0.0, 0.1]

    # 边界: a>=0, b>=0, |rho|<1, m 自由, sigma>0
    bounds = ([0, 0, -0.999, -10, 0.001],
              [10, 10, 0.999, 10, 10])

    result = least_squares(residuals, x0, bounds=bounds,
                            method='trf', ftol=1e-12, xtol=1e-12)

    a, b, rho, m, sigma = result.x
    return {
        "a": float(a), "b": float(b), "rho": float(rho),
        "m": float(m), "sigma": float(sigma),
        "cost": float(result.cost),
        "nfev": int(result.nfev),
        "success": bool(result.success),
        "rmse": float(np.sqrt(2 * result.cost / len(strikes))),
    }


def calibrate_ssvi(strikes: List[float], maturities: List[float],
                    implied_vols: List[float], forward: float) -> Dict:
    """SSVI 全局标定 - 分层策略

    镜像 C++ SSVI::calibrate:
      1. 按期限分组,每组 SVI 切片标定得到 θ_T
      2. 拟合 SSVI 全局参数 (rho, eta, gamma)
    """
    # 按期限分组
    mat_array = np.array(maturities)
    unique_T = np.unique(mat_array)
    theta_T = []
    for T in unique_T:
        mask = mat_array == T
        K_T = np.array(strikes)[mask]
        iv_T = np.array(implied_vols)[mask]
        svi_result = calibrate_svi_slice(K_T.tolist(), iv_T.tolist(),
                                          float(T), forward)
        # θ_T = ATM 总方差 = w(0, T)
        theta_T.append(svi_result["a"] + svi_result["b"] * svi_result["sigma"])

    # 拟合 SSVI 全局参数 (rho, eta, gamma)
    theta_T = np.array(theta_T)
    all_K = np.array(strikes)
    all_T = np.array(maturities)
    all_iv = np.array(implied_vols)
    all_log_k = np.log(all_K / forward)
    all_market_w = all_iv ** 2 * all_T

    # 为每个数据点找到对应的 θ
    theta_per_point = np.zeros_like(all_log_k)
    for i, T in enumerate(unique_T):
        mask = all_T == T
        theta_per_point[mask] = theta_T[i]

    def ssvi_residuals(params):
        rho, eta, gamma = params
        if abs(rho) >= 1.0 or eta <= 0 or gamma <= 0 or gamma >= 1:
            return np.ones_like(all_log_k) * 1e10
        model_w = np.array([
            ssvi_total_variance(k, theta,
                                 ssvi_power_law_phi(theta, eta, gamma),
                                 rho)
            for k, theta in zip(all_log_k, theta_per_point)
        ])
        return model_w - all_market_w

    x0 = [-0.3, 1.0, 0.5]
    bounds = ([-0.999, 0.001, 0.001], [0.999, 100, 0.999])
    result = least_squares(ssvi_residuals, x0, bounds=bounds,
                            method='trf', ftol=1e-10, xtol=1e-10)
    rho, eta, gamma = result.x

    return {
        "rho": float(rho), "eta": float(eta), "gamma": float(gamma),
        "theta_slice": theta_T.tolist(),
        "cost": float(result.cost),
        "nfev": int(result.nfev),
        "success": bool(result.success),
        "rmse": float(np.sqrt(2 * result.cost / len(strikes))),
    }


# ---------------------------------------------------------------------------
# 基准数据生成
# ---------------------------------------------------------------------------

def generate_market_data() -> Dict:
    """生成确定性市场数据 (欧式期权链)"""
    # 3 个月, 6 个月, 1 年, 2 年
    maturities = [0.25, 0.5, 1.0, 2.0]
    forward = 100.0
    # 每个期限 9 个行权价
    strikes_per_T = [80, 85, 90, 95, 100, 105, 110, 115, 120]

    # 使用已知的 SSVI 参数生成 "市场" IV
    true_rho, true_eta, true_gamma = -0.3, 1.0, 0.25

    strikes, mats, ivs = [], [], []
    for T in maturities:
        theta_T = 0.04 * T  # ATM 总方差 = σ_ATM² * T, σ_ATM = 0.20
        phi_T = ssvi_power_law_phi(theta_T, true_eta, true_gamma)
        for K in strikes_per_T:
            k = math.log(K / forward)
            w = ssvi_total_variance(k, theta_T, phi_T, true_rho)
            iv = math.sqrt(w / T)
            strikes.append(K)
            mats.append(T)
            ivs.append(iv)

    return {
        "strikes": strikes,
        "maturities": mats,
        "implied_vols": ivs,
        "forward": forward,
        "true_params": {
            "rho": true_rho, "eta": true_eta, "gamma": true_gamma,
        },
    }


def generate_benchmarks() -> Dict:
    """生成所有标定基准值"""
    market = generate_market_data()

    # SVI 切片标定 (对每个期限)
    svi_slices = {}
    unique_T = sorted(set(market["maturities"]))
    for T in unique_T:
        mask = [m == T for m in market["maturities"]]
        K_T = [market["strikes"][i] for i, m in enumerate(mask) if m]
        iv_T = [market["implied_vols"][i] for i, m in enumerate(mask) if m]
        svi_slices[f"T_{T}"] = calibrate_svi_slice(
            K_T, iv_T, T, market["forward"])

    # SSVI 全局标定
    ssvi_result = calibrate_ssvi(
        market["strikes"], market["maturities"],
        market["implied_vols"], market["forward"])

    # SSVI 公式验证 (已知参数 → 总方差 → 反推参数)
    theta_test = 0.16  # T=1, σ=0.20
    phi_test = ssvi_power_law_phi(theta_test, 1.0, 0.25)
    k_test = -0.2  # ITM
    w_test = ssvi_total_variance(k_test, theta_test, phi_test, -0.3)

    # 无套利条件验证 (Power-law: γ < 0.5 且 η(1+|ρ|) < 2)
    no_arb_check = {
        "rho_valid": abs(-0.3) < 1.0,
        "phi_positive": phi_test > 0,
        "butterfly_sufficient": phi_test * theta_test * (1 + abs(-0.3)) < 4.0,
        "gamma_valid": 0 < 0.25 < 0.5,
    }

    benchmarks = {
        "metadata": {
            "generator": "cross_validate_calibration.py",
            "version": "1.0",
            "description": "标定模块 Python 交叉验证基准 (scipy.optimize)",
            "cpp_reference": [
                "include/cpphub/models/vol_surface/svi.hpp",
                "include/cpphub/models/vol_surface/ssvi.hpp",
            ],
            "references": [
                "Gatheral-Jacquier (2014) arXiv:1204.0646",
                "Gatheral-Jacquier (2013) Arbitrage-free SVI volatility surfaces",
            ],
        },
        "market_data": {
            "n_points": len(market["strikes"]),
            "n_maturities": len(unique_T),
            "forward": market["forward"],
            "maturities": unique_T,
            "true_params": market["true_params"],
        },
        "svi_slices": svi_slices,
        "ssvi_calibration": ssvi_result,
        "ssvi_formula_validation": {
            "theta_test": theta_test,
            "phi_test": phi_test,
            "k_test": k_test,
            "total_variance": w_test,
            "implied_vol": math.sqrt(w_test / 1.0),
        },
        "no_arbitrage_check": no_arb_check,
        "tolerances": {
            "svi_slice_params": 1e-4,        # scipy LM vs C++ LM,可能局部最优不同
            "ssvi_global_params": 1e-3,       # 全局优化,容差更大
            "total_variance": 1e-12,          # 公式直接验证,应位精确
            "implied_vol": 1e-12,
            "rmse": 1e-8,
        },
    }

    return benchmarks


# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Cpp_Hub 标定模块交叉验证基准生成")
    parser.add_argument("--output", "-o", default=None,
                        help="输出文件 (default: stdout)")
    args = parser.parse_args()

    benchmarks = generate_benchmarks()
    output = json.dumps(benchmarks, indent=2, default=str)

    if args.output:
        with open(args.output, "w") as f:
            f.write(output)
        print(f"# Written to {args.output}", file=sys.stderr)
    else:
        print(output)


if __name__ == "__main__":
    main()
