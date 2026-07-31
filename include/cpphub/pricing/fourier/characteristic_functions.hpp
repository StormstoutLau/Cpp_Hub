#pragma once
// SOURCE: Carr & Madan (1999) "Option valuation using the fast Fourier transform"
// SOURCE: Fang & Oosterlee (2009) "A novel pricing method for options based on Fourier-cosine series expansion"
// SOURCE: Madan, Carr & Chang (1998) "The Variance Gamma Process and Option Pricing"
// 模块: 通用特征函数 (Characteristic Function, CF) 接口与已知模型 CF 工厂
//
// 统一约定: 所有 CF 均为 ln(S_T) 的特征函数, 即
//   phi(u) = E[exp(iu * ln S_T)]
// 输入 Complex u (复数), 输出 Complex (复数).
//
// 调用方约定: CF 已绑定到期时间 T (即 phi 是 T-固定函数).
// 工厂函数 make_*_cf(...) 接受 T 参数并返回闭包.

#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/pricing/analytic/heston_cf.hpp"
#include <functional>
#include <complex>
#include <cmath>
#include <stdexcept>

namespace cpphub {
inline namespace v1 {

// 通用特征函数签名: 输入 u (复数), 输出 phi(u) (复数)
// 数学含义: phi(u) = E[exp(iu * ln S_T)]
using CharFn = std::function<Complex(Complex)>;

// ============ Black-Scholes-Merton (GBM) 特征函数 ============
// S_T = S_0 * exp((r - q - sigma^2/2) T + sigma * W_T)
// ln S_T ~ N(ln S_0 + (r-q-sigma^2/2)T, sigma^2 * T)
// phi(u) = exp(iu*ln S_0 + iu*(r-q-sigma^2/2)*T - sigma^2*u^2*T/2)
inline CharFn make_gbm_cf(Real S0, Real r, Real q, Real sigma, Real T) {
    if (S0 <= 0.0) throw std::invalid_argument("gbm_cf: S0 must be positive");
    if (sigma <= 0.0) throw std::invalid_argument("gbm_cf: sigma must be positive");
    if (T <= 0.0) throw std::invalid_argument("gbm_cf: T must be positive");
    Real ln_S0 = std::log(S0);
    Real mu = (r - q - 0.5 * sigma * sigma) * T;
    Real var = sigma * sigma * T;
    return [ln_S0, mu, var](Complex u) -> Complex {
        // phi(u) = exp(iu*ln_S0 + iu*mu - var*u^2/2)
        // 对复数 u (Carr-Madan 需要 u - (alpha+1)i): u^2 为复数平方, 不是 |u|^2
        Complex iu = Complex(0.0, 1.0) * u;  // i*u, 复数乘法
        Complex exponent = iu * (ln_S0 + mu) - Complex(var / 2.0, 0.0) * (u * u);
        return std::exp(exponent);
    };
}

// ============ Heston 特征函数 (复用 heston_cf.hpp) ============
// phi(u) = exp(C(u,T) + D(u,T)*v0 + iu*ln S_0) 已在 heston_characteristic_function 中实现
inline CharFn make_heston_cf(Real S0, Real r, Real q, const HestonCFParams& p, Real T) {
    if (S0 <= 0.0) throw std::invalid_argument("heston_cf: S0 must be positive");
    if (T <= 0.0) throw std::invalid_argument("heston_cf: T must be positive");
    // 拷贝参数 (闭包捕获 by value)
    HestonCFParams params = p;
    params.r = r;
    params.q = q;
    return [S0, T, params](Complex u) -> Complex {
        return heston_characteristic_function(u, T, S0, params);
    };
}

// ============ Variance Gamma (VG) 特征函数 ============
// Madan-Carr-Chang (1998): S_T = S_0 * exp((r - q + omega)T + X_T)
// X_T = theta * G_T + sigma * sqrt(G_T) * Z,  G_T ~ Gamma(T/nu, 1/nu) (速率参数化)
// omega = (1/nu) * ln(1 - theta*nu - sigma^2*nu/2)
// phi(u) = exp(iu*(ln S_0 + (r-q+omega)T)) * (1 - iu*theta*nu + sigma^2*nu*u^2/2)^(-T/nu)
inline CharFn make_vg_cf(Real S0, Real r, Real q,
                          Real sigma, Real nu, Real theta, Real T) {
    if (S0 <= 0.0) throw std::invalid_argument("vg_cf: S0 must be positive");
    if (sigma <= 0.0) throw std::invalid_argument("vg_cf: sigma must be positive");
    if (nu <= 0.0) throw std::invalid_argument("vg_cf: nu must be positive");
    if (T <= 0.0) throw std::invalid_argument("vg_cf: T must be positive");
    // Feller 条件: 1 - theta*nu - sigma^2*nu/2 > 0
    Real condition = 1.0 - theta * nu - 0.5 * sigma * sigma * nu;
    if (condition <= 0.0) {
        throw std::invalid_argument("vg_cf: Feller condition violated (1 - theta*nu - sigma^2*nu/2 > 0)");
    }
    Real omega = std::log(condition) / nu;
    Real ln_S0 = std::log(S0);
    Real drift = (r - q + omega) * T;
    Real T_over_nu = T / nu;
    Real half_sigma2_nu = 0.5 * sigma * sigma * nu;
    return [ln_S0, drift, T_over_nu, theta, nu, half_sigma2_nu](Complex u) -> Complex {
        // phi(u) = exp(iu*(ln_S0 + drift)) * (1 - iu*theta*nu + sigma^2*nu*u^2/2)^(-T/nu)
        Real u_re = u.real();
        Real u_im = u.imag();
        Complex iu_theta_nu = Complex(-u_im * theta * nu, u_re * theta * nu);
        // u^2 (复数平方) - 注意 VG CF 中 u^2 是复数平方, 不是 |u|^2
        Complex u_sq = u * u;
        Real u_sq_re = u_sq.real();
        Real u_sq_im = u_sq.imag();
        Complex base = Complex(1.0, 0.0) - iu_theta_nu
                     + Complex(half_sigma2_nu * u_sq_re, half_sigma2_nu * u_sq_im);
        // ln_S0 + drift 实数项
        Real real_shift = ln_S0 + drift;
        Complex iu_shift = Complex(-u_im * real_shift, u_re * real_shift);
        Complex phase = std::exp(iu_shift);
        Complex denom_pow = std::pow(base, -T_over_nu);
        return phase * denom_pow;
    };
}

// ============ Normal Inverse Gaussian (NIG) 特征函数 ============
// S_T = S_0 * exp((r - q)T + X_T), X_T ~ NIG(alpha, beta, delta*T)
// phi(u) = exp(iu*ln S_0 + iu*(r-q)T + delta*(sqrt(alpha^2 - beta^2) - sqrt(alpha^2 - (beta+iu)^2)))
inline CharFn make_nig_cf(Real S0, Real r, Real q,
                           Real alpha, Real beta, Real delta, Real T) {
    if (S0 <= 0.0) throw std::invalid_argument("nig_cf: S0 must be positive");
    if (alpha <= 0.0) throw std::invalid_argument("nig_cf: alpha must be positive");
    if (std::abs(beta) >= alpha) {
        throw std::invalid_argument("nig_cf: require |beta| < alpha");
    }
    if (delta <= 0.0) throw std::invalid_argument("nig_cf: delta must be positive");
    if (T <= 0.0) throw std::invalid_argument("nig_cf: T must be positive");
    Real ln_S0 = std::log(S0);
    Real drift = (r - q) * T;
    Real delta_T = delta * T;
    Real sqrt_a2_b2 = std::sqrt(alpha * alpha - beta * beta);
    Real a2 = alpha * alpha;
    return [ln_S0, drift, delta_T, sqrt_a2_b2, a2, beta](Complex u) -> Complex {
        Real u_re = u.real();
        Real u_im = u.imag();
        Real real_shift = ln_S0 + drift;
        Complex iu_shift = Complex(-u_im * real_shift, u_re * real_shift);
        Complex phase = std::exp(iu_shift);
        // (beta + iu)^2 = (beta - u_im)^2 - u_re^2 + 2i*(beta - u_im)*u_re
        Real bp_re = beta - u_im;
        Real bp_im = u_re;
        Real bp_sq_re = bp_re * bp_re - bp_im * bp_im;
        Real bp_sq_im = 2.0 * bp_re * bp_im;
        // sqrt(alpha^2 - (beta+iu)^2)
        Complex inside = Complex(a2 - bp_sq_re, -bp_sq_im);
        Complex sqrt_inside = std::sqrt(inside);
        Complex exponent = Complex(delta_T * sqrt_a2_b2, 0.0) - Complex(delta_T, 0.0) * sqrt_inside;
        return phase * std::exp(exponent);
    };
}

// ============ 辅助: 累积分布截断区间 (用于 COS 方法) ============
// 给定 CF, 数值估计 X = ln S_T 的 [a, b] 截断区间
// 基于 Fang-Oosterlee (2009) eq 27-30: 通过 CF 反演求 cdf, 找 1e-8 / 1-1e-8 分位数
// 简化版: 用 BSM 近似估计 (足够稳定, 通用性弱), 调用方可自行覆盖
inline std::pair<Real, Real> cos_truncation_range(const CharFn& phi, Real S0, Real T,
                                                   Real L = 10.0) {
    // 通过 CF 计算 ln S_T 的一阶/二阶矩:
    // E[ln S_T] = -Im[phi'(0)] / i (需数值微分)
    // Var[ln S_T] = -Re[phi''(0)] + (Re[phi'(0)])^2
    Real h = 1e-4;
    Complex phi0 = phi(Complex(0.0, 0.0));
    Complex phi_p = (phi(Complex(h, 0.0)) - phi(Complex(-h, 0.0))) / (2.0 * h);  // dphi/du
    Complex phi_pp = (phi(Complex(h, 0.0)) - 2.0 * phi0 + phi(Complex(-h, 0.0))) / (h * h);
    // E[X] = phi'(0) / (i) = -i * phi'(0)
    // E[X^2] = -phi''(0)
    Complex i_(0.0, 1.0);
    Real mean = (phi_p / i_).real();  // 应为实数
    Real second_moment = -phi_pp.real();
    Real variance = second_moment - mean * mean;
    if (variance < 0.0) variance = 0.01;  // 数值失败时回退
    Real std_dev = std::sqrt(variance);
    Real a = mean - L * std_dev;
    Real b = mean + L * std_dev;
    return {a, b};
}

}  // namespace v1
}  // namespace cpphub
