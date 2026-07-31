#pragma once
// SOURCE: Carr & Madan (1999) "Option valuation using the fast Fourier transform"
//         (J. Comp. Finance, 2(4), 61-73).
//
// Carr-Madan FFT 方法核心思路:
//   1. 通过乘以 e^{αk} (α>0 阻尼因子) 使 call 价格 c(k)=C(K) 在 k=ln K→-∞ 时可积
//   2. 对阻尼价格 c̃(k) = e^{αk} c(k) 做傅里叶变换, 得到解析的 ψ(v)
//   3. 用 FFT 在离散网格上数值反演 ψ(v), 得到一组 k 网格上的 c(k)
//
// 公式 (Carr-Madan 1999, eq. 6+9):
//   c(k) = (e^{-αk} / π) ∫_0^∞ Re[e^{-ivk} ψ(v)] dv
//   ψ(v) = e^{-rT} φ(v - (α+1)i) / (α² + α - v² + i(2α+1)v)
//
// FFT 离散化 (Simpson 法则加权):
//   c(k_u) ≈ (e^{-αk_u} / π) Σ_{j=0}^{N-1} e^{-i 2π j u / N} x_j
//   x_j = e^{-i v_j b} ψ(v_j) · (η/3) [3 + (-1)^j - δ_{j,0}]
//   其中 v_j = jη, k_u = -b + uλ, b = Nλ/2, λ = 2π/(Nη)
//
// 调用方约定: φ(u) = E[exp(iu ln S_T)] (与 characteristic_functions.hpp 一致)

#include "cpphub/core/types.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/pricing/fourier/characteristic_functions.hpp"
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ Radix-2 Cooley-Tukey FFT 内核 ============
// 迭代实现, 位反转排序 + 蝶形运算
// sign = -1: 正向 FFT (exp(-i 2π k n / N))
// sign = +1: 逆向 FFT (exp(+i 2π k n / N))
// 注意: 不做 1/N 归一化 (调用方自行处理)
inline void fft_radix2(std::vector<Complex>& a, int sign) {
    const Size n = a.size();
    if (n <= 1) return;
    if (n & (n - 1)) {
        throw std::invalid_argument("fft_radix2: size must be power of 2");
    }

    // 位反转排序 (Bit-reversal permutation)
    for (Size i = 1, j = 0; i < n; ++i) {
        Size bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    // 蝶形运算 (Cooley-Tukey)
    for (Size len = 2; len <= n; len <<= 1) {
        Real angle = sign * 2.0 * PI / static_cast<Real>(len);
        Complex wlen(std::cos(angle), std::sin(angle));
        for (Size i = 0; i < n; i += len) {
            Complex w(1.0, 0.0);
            for (Size j = 0; j < len / 2; ++j) {
                Complex u = a[i + j];
                Complex v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

// ============ Carr-Madan FFT 引擎 ============
class CarrMadanFFT {
public:
    struct Config {
        Real alpha = 1.5;        // 阻尼因子 (典型 1.5-5, 短期期权需较大值)
        Size n_fft = 4096;       // FFT 点数 N (必须为 2 的幂)
        Real eta = 0.25;         // v 空间步长 (典型 0.1-0.5)
        // 可选: 限制输出 strike 网格的中心 (通常设为 ln(S0))
        // 网格: k_u = -b + u*lambda, b = N*lambda/2, lambda = 2π/(N*eta)
        // 若不指定, 自动以 ln(S0) 为中心选择最近的网格点
    };

    explicit CarrMadanFFT(CharFn phi, Real S0, Real r, Real q, Real T, Config cfg = {})
        : phi_(std::move(phi)), S0_(S0), r_(r), q_(q), T_(T), cfg_(cfg) {
        if (S0 <= 0.0) throw std::invalid_argument("CarrMadanFFT: S0 must be positive");
        if (T <= 0.0) throw std::invalid_argument("CarrMadanFFT: T must be positive");
        if (cfg_.alpha <= 0.0) throw std::invalid_argument("CarrMadanFFT: alpha must be positive");
        // 验证 n_fft 是 2 的幂
        Size n = cfg_.n_fft;
        if (n < 8 || (n & (n - 1)) != 0) {
            throw std::invalid_argument("CarrMadanFFT: n_fft must be power of 2 and >= 8");
        }
        if (cfg_.eta <= 0.0) throw std::invalid_argument("CarrMadanFFT: eta must be positive");

        // 预计算 strike 网格
        lambda_ = 2.0 * PI / (static_cast<Real>(n) * cfg_.eta);
        b_grid_ = static_cast<Real>(n) * lambda_ / 2.0;
        strikes_.resize(n);
        log_strikes_.resize(n);
        for (Size u = 0; u < n; ++u) {
            log_strikes_[u] = -b_grid_ + static_cast<Real>(u) * lambda_;
            strikes_[u] = std::exp(log_strikes_[u]);
        }
    }

    // 计算整个 strike 网格上的 call 价格
    // 返回 (strikes, call_prices) 对, 长度均为 n_fft
    std::pair<std::vector<Real>, std::vector<Real>> price_call_grid() const {
        const Size n = cfg_.n_fft;
        const Real alpha = cfg_.alpha;
        const Real eta = cfg_.eta;
        const Real disc = std::exp(-r_ * T_);

        // 构建 FFT 输入 x_j = e^{-i v_j b} ψ(v_j) * Simpson 权重
        std::vector<Complex> x(n, Complex(0.0, 0.0));
        for (Size j = 0; j < n; ++j) {
            Real v_j = static_cast<Real>(j) * eta;
            // ψ(v) = e^{-rT} φ(v - (α+1)i) / (α² + α - v² + i(2α+1)v)
            Complex u_shifted = Complex(v_j, -(alpha + 1.0));  // v - (α+1)i
            Complex phi_val = phi_(u_shifted);
            Complex num = disc * phi_val;
            Complex denom = Complex(alpha * alpha + alpha - v_j * v_j,
                                     (2.0 * alpha + 1.0) * v_j);
            Complex psi_v = num / denom;

            // e^{-i v_j b}
            Complex phase = std::exp(Complex(0.0, -v_j * b_grid_));

            // Simpson 权重: (η/3) [3 - (-1)^j - δ_{j,0}]
            // = (η/3) [1, 4, 2, 4, 2, ..., 4, 2, ...]  (标准 Simpson 法则)
            // 注: j=0 时 δ_{j,0}=1, 权重 = (η/3)(3-1-1) = η/3
            //     j 奇数时 权重 = (η/3)(3+1) = 4η/3
            //     j 偶数 (j>0) 时 权重 = (η/3)(3-1) = 2η/3
            Real w_simpson = (eta / 3.0) * (3.0 - ((j % 2 == 0) ? 1.0 : -1.0)
                                             - (j == 0 ? 1.0 : 0.0));
            x[j] = phase * psi_v * w_simpson;
        }

        // 逆 FFT (sign = -1, 不归一化)
        // 注意: Carr-Madan 公式中的求和 Σ e^{-i 2π j u / N} x_j 对应正向 FFT (sign=-1)
        fft_radix2(x, -1);

        // c(k_u) = e^{-α k_u} / π * Re[FFT(x)_u]
        // 注意: FFT 输出未归一化, 但我们在 Simpson 权重中已隐含 η 因子
        std::vector<Real> calls(n);
        for (Size u = 0; u < n; ++u) {
            Real k_u = log_strikes_[u];
            calls[u] = (std::exp(-alpha * k_u) / PI) * x[u].real();
            if (calls[u] < 0.0) calls[u] = 0.0;
        }

        return {strikes_, calls};
    }

    // 通过线性插值获取指定行权价 K 的 call 价格
    Real price_call(Real K) const {
        if (K <= 0.0) throw std::invalid_argument("CarrMadanFFT: K must be positive");
        auto [ks, calls] = price_call_grid();
        return interpolate_log_strike(std::log(K), ks, calls);
    }

    // 通过 put-call parity 获取 put 价格
    // P = C - S0 * e^{-qT} + K * e^{-rT}
    Real price_put(Real K) const {
        if (K <= 0.0) throw std::invalid_argument("CarrMadanFFT: K must be positive");
        Real C = price_call(K);
        return C - S0_ * std::exp(-q_ * T_) + K * std::exp(-r_ * T_);
    }

    // 访问器
    const std::vector<Real>& strikes() const { return strikes_; }
    const std::vector<Real>& log_strikes() const { return log_strikes_; }
    Real S0() const { return S0_; }
    Real T() const { return T_; }
    Real alpha() const { return cfg_.alpha; }
    Size n_fft() const { return cfg_.n_fft; }
    Real eta() const { return cfg_.eta; }
    Real lambda() const { return lambda_; }

private:
    CharFn phi_;
    Real S0_, r_, q_, T_;
    Config cfg_;
    Real lambda_;       // k 空间步长
    Real b_grid_;       // k 空间半宽
    std::vector<Real> strikes_;
    std::vector<Real> log_strikes_;

    // 在 log-strike 网格上做线性插值
    // k 网格: k_u = -b + u*lambda, 等间距, 单调递增
    static Real interpolate_log_strike(Real k_target,
                                        const std::vector<Real>& strikes,
                                        const std::vector<Real>& calls) {
        const Size n = strikes.size();
        Real lambda = (n > 1) ? (std::log(strikes[1]) - std::log(strikes[0])) : 0.0;
        Real b = (n > 0) ? -std::log(strikes[0]) : 0.0;  // k_0 = -b
        // 找到 k_target 在网格中的位置: u = (k_target + b) / lambda
        Real u_real = (k_target + b) / lambda;
        if (u_real < 0.0) u_real = 0.0;
        if (u_real > static_cast<Real>(n - 1)) u_real = static_cast<Real>(n - 1);
        Size u_lo = static_cast<Size>(std::floor(u_real));
        Size u_hi = (u_lo + 1 < n) ? u_lo + 1 : u_lo;
        Real frac = u_real - static_cast<Real>(u_lo);
        return calls[u_lo] * (1.0 - frac) + calls[u_hi] * frac;
    }
};

// ============ 便捷工厂函数 ============
// BSM call via Carr-Madan FFT (用于验证)
inline Real fft_call_gbm(Real S0, Real K, Real T, Real r, Real q, Real sigma,
                          Real alpha = 1.5, Size n_fft = 4096, Real eta = 0.25) {
    auto phi = make_gbm_cf(S0, r, q, sigma, T);
    CarrMadanFFT::Config cfg;
    cfg.alpha = alpha;
    cfg.n_fft = n_fft;
    cfg.eta = eta;
    CarrMadanFFT engine(phi, S0, r, q, T, cfg);
    return engine.price_call(K);
}

// Heston call via Carr-Madan FFT
inline Real fft_call_heston(Real S0, Real K, Real T, Real r, Real q,
                             const HestonCFParams& hp,
                             Real alpha = 1.5, Size n_fft = 4096, Real eta = 0.25) {
    auto phi = make_heston_cf(S0, r, q, hp, T);
    CarrMadanFFT::Config cfg;
    cfg.alpha = alpha;
    cfg.n_fft = n_fft;
    cfg.eta = eta;
    CarrMadanFFT engine(phi, S0, r, q, T, cfg);
    return engine.price_call(K);
}

}  // namespace v1
}  // namespace cpphub
