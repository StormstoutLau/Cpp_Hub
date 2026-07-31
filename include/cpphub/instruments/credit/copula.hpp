#pragma once
// SOURCE: Li (2000) "On Default Correlation: A Copula Function Approach" J. Fixed Income 9(4), 43-54
// SOURCE: Nelsen (2006) "An Introduction to Copulas" 2nd ed., Springer
// SOURCE: Schönbucher (2003) "Credit Derivatives Pricing Models" Ch.10 (Copula)
// SOURCE: O'Kane (2008) "Modelling Single-Name and Multi-Name Credit Derivatives" Ch.9 (Copula)
// SOURCE: Hull & White (2004) "Valuation of a CDO and an n-th to Default CDS without Monte Carlo Simulation"
// SOURCE: Demarta & McNeil (2005) "The t Copula and Related Copulas" Int. Stat. Review 73(1), 111-129
// 模块: Copula 模型 (Gaussian / t) + 违约时间采样
//
// ==================== Copula 框架 ====================
//
// Sklar 定理: 任意联合分布 F(x_1,...,x_n) 可分解为边际 F_i 与一个 Copula C:
//   F(x_1,...,x_n) = C(F_1(x_1), ..., F_n(x_n))
// Copula 捕捉"依赖结构", 与边际分布解耦.
//
// 在信用衍生品中, 用于建模多个名字违约时间 τ_1, ..., τ_n 的联合分布:
//   P(τ_1 ≤ t_1, ..., τ_n ≤ t_n) = C(U_1, ..., U_n), U_i = F_i(t_i) = PD_i(0, t_i)
// 其中 F_i 为第 i 个名字违约时间的边际分布 (由 CreditCurve 给出).
//
// ==================== Gaussian Copula (Li 2000) ====================
//
// 1. 生成相关高斯向量 X = (X_1, ..., X_n) ~ N(0, Σ), Σ 为相关矩阵
//    通过 Cholesky 分解 Σ = LL^T, X = L * Z, Z ~ N(0, I)
// 2. 转换为均匀变量 U_i = Φ(X_i), Φ 为标准正态 CDF
// 3. 通过边际反函数得到违约时间: τ_i = F_i^{-1}(U_i) = -ln(1 - U_i) / h_i  (flat hazard 近似)
//
// 单因子 Gaussian Copula (Li 2000 / Vasicek 1987):
//   X_i = √ρ * M + √(1-ρ) * Z_i
//   M ~ N(0,1) 为系统因子 (共同市场风险), Z_i ~ N(0,1) 为特质因子
//   ρ 为资产相关性 (equity correlation proxy)
//   条件独立: 给定 M=m, 各名字违约独立, 违约概率为 PD_i|m = Φ((Φ^{-1}(PD_i) - √ρ * m) / √(1-ρ))
//
// ==================== t-Copula ====================
//
// t-Copula 捕捉尾部相关性 (extreme joint defaults), Gaussian Copula 无尾部相关:
//   1. 生成 (X, W): X ~ N(0, Σ), W ~ χ²(ν) 独立
//   2. Y_i = X_i / √(W / ν) ~ t_ν(0, Σ), 自由度 ν 控制尾部厚度
//   3. U_i = t_ν(Y_i), τ_i = F_i^{-1}(U_i)
// ν 越小, 尾部相关性越强; ν → ∞ 退化为 Gaussian Copula
//
// ==================== 违约时间反演 (flat hazard 近似) ====================
//
// 给定 CreditCurve, 生存概率 Q_i(t) = exp(-H_i(t))
//   U_i = F_i(τ_i) = 1 - Q_i(τ_i)  →  Q_i(τ_i) = 1 - U_i
//   H_i(τ_i) = -ln(1 - U_i)
// 若 hazard rate 平坦 h_i: τ_i = -ln(1 - U_i) / h_i
// 一般情况: 数值求解 H_i(τ_i) = -ln(1 - U_i), 即 τ_i = H_i^{-1}(-ln(1 - U_i))
//
// ==================== 约定 ====================
//
// - 输入: 独立均匀变量 [0,1]^n → 输出: 相关均匀变量 [0,1]^n (或违约时间)
// - 相关矩阵 Σ: n×n 对称正定 (或半正定), 对角线为 1
// - 自由度 ν: t-Copula 参数, ν > 2 通常 (保证协方差有限)
// - 违约时间: +∞ 表示不违约 (U_i = 0 时, 即 1-U_i = 1)

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/rng.hpp"
#include "cpphub/instruments/credit/credit_curve.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>

namespace cpphub {
inline namespace v1 {

// ============ 动态 Cholesky 分解 (半正定, 允许 ρ=±1) ============
// 复用 multi_asset_path_generator.hpp 的 cholesky_semi_definite 逻辑,
// 但独立实现以避免循环依赖.
inline std::vector<std::vector<Real>> cholesky_dynamic(
        const std::vector<std::vector<Real>>& A) {
    const Size n = A.size();
    if (n == 0) return {};
    for (Size i = 0; i < n; ++i) {
        if (A[i].size() != n) {
            throw std::invalid_argument("cholesky_dynamic: matrix must be square");
        }
    }
    std::vector<std::vector<Real>> L(n, std::vector<Real>(n, 0.0));
    for (Size i = 0; i < n; ++i) {
        for (Size j = 0; j <= i; ++j) {
            Real s = A[i][j];
            for (Size k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
            if (i == j) {
                if (s < -1e-12) {
                    throw std::invalid_argument(
                        "cholesky_dynamic: matrix not positive semidefinite");
                }
                L[i][j] = std::sqrt(std::max(s, 0.0));
            } else {
                L[i][j] = (L[j][j] > 1e-15) ? s / L[j][j] : 0.0;
            }
        }
    }
    return L;
}

// ============ Gaussian Copula ============
// 通过相关矩阵 Σ 参数化, 输出相关均匀变量 (或正态变量).
// 单因子模型可由 Σ_ij = ρ (i≠j), Σ_ii = 1 实现.
class GaussianCopula {
public:
    explicit GaussianCopula(std::vector<std::vector<Real>> correlation)
        : corr_(std::move(correlation)) {
        validate();
        L_ = cholesky_dynamic(corr_);
    }

    Size n_names() const noexcept { return corr_.size(); }
    const std::vector<std::vector<Real>>& correlation() const noexcept { return corr_; }
    const std::vector<std::vector<Real>>& cholesky_L() const noexcept { return L_; }

    // 验证相关矩阵: 对称, 对角线=1, 半正定 (由 cholesky 检查)
    void validate() const {
        const Size n = corr_.size();
        if (n == 0) throw std::invalid_argument("GaussianCopula: empty correlation");
        for (Size i = 0; i < n; ++i) {
            if (corr_[i].size() != n) {
                throw std::invalid_argument("GaussianCopula: matrix not square");
            }
            if (std::abs(corr_[i][i] - 1.0) > 1e-10) {
                throw std::invalid_argument("GaussianCopula: diagonal must be 1");
            }
            for (Size j = 0; j < i; ++j) {
                if (std::abs(corr_[i][j] - corr_[j][i]) > 1e-10) {
                    throw std::invalid_argument("GaussianCopula: matrix not symmetric");
                }
                if (std::abs(corr_[i][j]) > 1.0 + 1e-10) {
                    throw std::invalid_argument("GaussianCopula: |correlation| > 1");
                }
            }
        }
    }

    // 将 n 个独立均匀变量转换为相关均匀变量
    // 输入: u_indep[i] ∈ (0,1) 独立均匀
    // 输出: u_corr[i] = Φ((L * Φ^{-1}(u_indep))[i]), Φ 为标准正态 CDF
    std::vector<Real> transform_uniforms(const std::vector<Real>& u_indep) const {
        const Size n = n_names();
        if (u_indep.size() != n) {
            throw std::invalid_argument("GaussianCopula::transform_uniforms: size mismatch");
        }
        // 1. 转换为独立标准正态
        std::vector<Real> z(n);
        for (Size i = 0; i < n; ++i) {
            // 保护边界, 避免 inv_normal_cdf 返回 ±∞
            Real u = std::min(std::max(u_indep[i], 1e-15), 1.0 - 1e-15);
            z[i] = inv_normal_cdf(u);
        }
        // 2. 相关化: X = L * Z
        std::vector<Real> x(n, 0.0);
        for (Size i = 0; i < n; ++i) {
            for (Size j = 0; j <= i; ++j) {
                x[i] += L_[i][j] * z[j];
            }
        }
        // 3. 转换为均匀: U = Φ(X)
        std::vector<Real> u_corr(n);
        for (Size i = 0; i < n; ++i) {
            u_corr[i] = normal_cdf(x[i]);
        }
        return u_corr;
    }

    // 从 RNG 直接生成相关均匀变量 (一步完成)
    std::vector<Real> sample_uniforms(Philox4x64& rng) const {
        const Size n = n_names();
        std::vector<Real> u_indep(n);
        for (Size i = 0; i < n; ++i) {
            uint64_t r = rng();
            u_indep[i] = (r >> 11) * (1.0 / 9007199254740992.0);
            // 保护边界
            if (u_indep[i] < 1e-15) u_indep[i] = 1e-15;
            if (u_indep[i] > 1.0 - 1e-15) u_indep[i] = 1.0 - 1e-15;
        }
        return transform_uniforms(u_indep);
    }

    // 从相关均匀变量 + 信用曲线, 得到违约时间
    // τ_i = H_i^{-1}(-ln(1 - U_i)), H_i 为累积 hazard
    // 若 U_i ≤ 0 (不可能发生), 返回 +∞ (不违约)
    // 若 U_i → 1, 违约时间 → ∞
    // 输入: u_corr (相关均匀), credit_curves (各名字信用曲线)
    // 返回: default_times[i] = τ_i (可能为 +∞)
    std::vector<Real> uniforms_to_default_times(
            const std::vector<Real>& u_corr,
            const std::vector<CreditCurve>& credit_curves) const {
        const Size n = n_names();
        if (credit_curves.size() != n) {
            throw std::invalid_argument(
                "GaussianCopula::uniforms_to_default_times: credit_curves size mismatch");
        }
        if (u_corr.size() != n) {
            throw std::invalid_argument(
                "GaussianCopula::uniforms_to_default_times: u_corr size mismatch");
        }
        std::vector<Real> tau(n);
        for (Size i = 0; i < n; ++i) {
            Real u = u_corr[i];
            if (u <= 0.0) {
                tau[i] = std::numeric_limits<Real>::infinity();
                continue;
            }
            // 1 - U_i = Q_i(τ_i), 求解 τ_i = H_i^{-1}(-ln(1-U_i))
            Real target_h = -std::log(1.0 - std::min(u, 1.0 - 1e-15));
            tau[i] = solve_default_time_(credit_curves[i], target_h);
        }
        return tau;
    }

    // 一步采样违约时间 (从 RNG)
    std::vector<Real> sample_default_times(
            Philox4x64& rng,
            const std::vector<CreditCurve>& credit_curves) const {
        auto u = sample_uniforms(rng);
        return uniforms_to_default_times(u, credit_curves);
    }

private:
    std::vector<std::vector<Real>> corr_;
    std::vector<std::vector<Real>> L_;

    // 数值求解 H(τ) = target_h, 即 τ = H^{-1}(target_h)
    // H 单调递增, 用二分法
    // 搜索区间 [0, T_max], T_max = 100 年 (足够覆盖)
    static Real solve_default_time_(const CreditCurve& credit, Real target_h) {
        if (target_h <= 0.0) return std::numeric_limits<Real>::infinity();
        const Real T_max = 100.0;
        const Real h_tol = 1e-10;
        const Size max_iter = 100;

        // 检查 T_max 是否足够大
        Real h_max = credit.cumulative_hazard(T_max);
        if (h_max < target_h) {
            // hazard 太小, 违约时间超出 T_max → 返回 T_max (近似)
            // 或返回 +∞ (不违约). 这里返回 T_max 作为保守估计.
            return T_max;
        }

        Real t_lo = 0.0;
        Real t_hi = T_max;
        for (Size iter = 0; iter < max_iter; ++iter) {
            Real t_mid = 0.5 * (t_lo + t_hi);
            Real h_mid = credit.cumulative_hazard(t_mid);
            if (std::abs(h_mid - target_h) < h_tol) return t_mid;
            if (h_mid < target_h) t_lo = t_mid;
            else t_hi = t_mid;
            if (t_hi - t_lo < 1e-12) break;
        }
        return 0.5 * (t_lo + t_hi);
    }
};

// ============ 单因子 Gaussian Copula (Vasicek / Li 2000) ============
// 参数化: X_i = √ρ * M + √(1-ρ) * Z_i
// 优势: 计算高效, 条件独立框架支持半解析 tranche 定价 (LHA / recursive)
class OneFactorGaussianCopula {
public:
    OneFactorGaussianCopula(Real rho, Size n_names)
        : rho_(rho), n_(n_names) {
        if (rho_ < 0.0 || rho_ >= 1.0) {
            throw std::invalid_argument("OneFactorGaussianCopula: rho must be in [0, 1)");
        }
        if (n_ == 0) {
            throw std::invalid_argument("OneFactorGaussianCopula: n_names must be positive");
        }
        sqrt_rho_ = std::sqrt(rho_);
        sqrt_1mrho_ = std::sqrt(1.0 - rho_);
    }

    Real rho() const noexcept { return rho_; }
    Size n_names() const noexcept { return n_; }

    // 从两个独立均匀变量生成相关均匀向量
    // u_systemic: 共同因子 M 的均匀变量
    // u_idiosyncratic[i]: 第 i 个名字特质因子 Z_i 的均匀变量
    std::vector<Real> transform_uniforms(
            Real u_systemic,
            const std::vector<Real>& u_idiosyncratic) const {
        if (u_idiosyncratic.size() != n_) {
            throw std::invalid_argument(
                "OneFactorGaussianCopula::transform_uniforms: size mismatch");
        }
        Real u_sys = std::min(std::max(u_systemic, 1e-15), 1.0 - 1e-15);
        Real m = inv_normal_cdf(u_sys);
        std::vector<Real> u_corr(n_);
        for (Size i = 0; i < n_; ++i) {
            Real u_id = std::min(std::max(u_idiosyncratic[i], 1e-15), 1.0 - 1e-15);
            Real z = inv_normal_cdf(u_id);
            Real x = sqrt_rho_ * m + sqrt_1mrho_ * z;
            u_corr[i] = normal_cdf(x);
        }
        return u_corr;
    }

    // 从 RNG 一步采样相关均匀向量
    std::vector<Real> sample_uniforms(Philox4x64& rng) const {
        Real u_sys = uniform_from_rng_(rng);
        std::vector<Real> u_id(n_);
        for (Size i = 0; i < n_; ++i) u_id[i] = uniform_from_rng_(rng);
        return transform_uniforms(u_sys, u_id);
    }

    // 条件违约概率: 给定 M=m, 第 i 个名字在 t 前违约的概率
    // PD_i(t | m) = Φ((Φ^{-1}(PD_i(t)) - √ρ * m) / √(1-ρ))
    // 这是单因子模型的核心公式, 用于条件独立损失分布计算
    Real conditional_pd(Real pd_marginal, Real m) const {
        Real u = std::min(std::max(pd_marginal, 1e-15), 1.0 - 1e-15);
        Real z = inv_normal_cdf(u);
        Real x = (z - sqrt_rho_ * m) / sqrt_1mrho_;
        return normal_cdf(x);
    }

    // 生成相关矩阵 (n×n, Σ_ij = ρ for i≠j, Σ_ii = 1)
    // 用于等价的高斯 Copula 视角
    std::vector<std::vector<Real>> correlation_matrix() const {
        std::vector<std::vector<Real>> corr(n_, std::vector<Real>(n_, rho_));
        for (Size i = 0; i < n_; ++i) corr[i][i] = 1.0;
        return corr;
    }

    // 从相关均匀变量得到违约时间 (委托给 GaussianCopula 的等价逻辑)
    std::vector<Real> uniforms_to_default_times(
            const std::vector<Real>& u_corr,
            const std::vector<CreditCurve>& credit_curves) const {
        if (credit_curves.size() != n_ || u_corr.size() != n_) {
            throw std::invalid_argument(
                "OneFactorGaussianCopula::uniforms_to_default_times: size mismatch");
        }
        std::vector<Real> tau(n_);
        for (Size i = 0; i < n_; ++i) {
            Real u = u_corr[i];
            if (u <= 0.0) {
                tau[i] = std::numeric_limits<Real>::infinity();
                continue;
            }
            Real target_h = -std::log(1.0 - std::min(u, 1.0 - 1e-15));
            tau[i] = solve_default_time_(credit_curves[i], target_h);
        }
        return tau;
    }

    // 一步采样违约时间
    std::vector<Real> sample_default_times(
            Philox4x64& rng,
            const std::vector<CreditCurve>& credit_curves) const {
        auto u = sample_uniforms(rng);
        return uniforms_to_default_times(u, credit_curves);
    }

private:
    Real rho_;
    Size n_;
    Real sqrt_rho_;
    Real sqrt_1mrho_;

    static Real uniform_from_rng_(Philox4x64& rng) {
        uint64_t r = rng();
        Real u = (r >> 11) * (1.0 / 9007199254740992.0);
        if (u < 1e-15) u = 1e-15;
        if (u > 1.0 - 1e-15) u = 1.0 - 1e-15;
        return u;
    }

    static Real solve_default_time_(const CreditCurve& credit, Real target_h) {
        if (target_h <= 0.0) return std::numeric_limits<Real>::infinity();
        const Real T_max = 100.0;
        const Real h_tol = 1e-10;
        const Size max_iter = 100;

        Real h_max = credit.cumulative_hazard(T_max);
        if (h_max < target_h) return T_max;

        Real t_lo = 0.0, t_hi = T_max;
        for (Size iter = 0; iter < max_iter; ++iter) {
            Real t_mid = 0.5 * (t_lo + t_hi);
            Real h_mid = credit.cumulative_hazard(t_mid);
            if (std::abs(h_mid - target_h) < h_tol) return t_mid;
            if (h_mid < target_h) t_lo = t_mid;
            else t_hi = t_mid;
            if (t_hi - t_lo < 1e-12) break;
        }
        return 0.5 * (t_lo + t_hi);
    }
};

// ============ t-Copula ============
// 在 Gaussian Copula 基础上引入 χ²(ν) 缩放, 实现尾部相关性.
// Y_i = X_i / √(W/ν), W ~ χ²(ν), X ~ N(0, Σ) 独立
class TCopula {
public:
    TCopula(std::vector<std::vector<Real>> correlation, Real degrees_of_freedom)
        : gaussian_(std::move(correlation)),
          nu_(degrees_of_freedom) {
        if (nu_ <= 2.0) {
            throw std::invalid_argument("TCopula: degrees_of_freedom must be > 2");
        }
    }

    Size n_names() const noexcept { return gaussian_.n_names(); }
    Real degrees_of_freedom() const noexcept { return nu_; }
    const std::vector<std::vector<Real>>& correlation() const noexcept {
        return gaussian_.correlation();
    }

    // 将 n 个独立均匀变量 + 一个 χ² 样本转换为相关均匀变量
    // u_indep[i]: 独立均匀 (用于生成 Z_i)
    // w_sample: χ²(ν) 的样本 (用于缩放)
    std::vector<Real> transform_uniforms(
            const std::vector<Real>& u_indep, Real w_sample) const {
        const Size n = n_names();
        if (u_indep.size() != n) {
            throw std::invalid_argument("TCopula::transform_uniforms: size mismatch");
        }
        if (w_sample <= 0.0) {
            throw std::invalid_argument("TCopula::transform_uniforms: w_sample must be positive");
        }
        // 1. 独立正态 Z
        std::vector<Real> z(n);
        for (Size i = 0; i < n; ++i) {
            Real u = std::min(std::max(u_indep[i], 1e-15), 1.0 - 1e-15);
            z[i] = inv_normal_cdf(u);
        }
        // 2. 相关正态 X = L * Z
        const auto& L = gaussian_.cholesky_L();
        std::vector<Real> x(n, 0.0);
        for (Size i = 0; i < n; ++i) {
            for (Size j = 0; j <= i; ++j) x[i] += L[i][j] * z[j];
        }
        // 3. t 缩放: Y = X / √(W/ν)
        Real scale = std::sqrt(nu_ / w_sample);
        std::vector<Real> y(n);
        for (Size i = 0; i < n; ++i) y[i] = x[i] * scale;
        // 4. 均匀变量: U = t_ν(Y)
        std::vector<Real> u_corr(n);
        for (Size i = 0; i < n; ++i) u_corr[i] = student_t_cdf_(y[i], nu_);
        return u_corr;
    }

    // 从 RNG 一步采样相关均匀向量 (内部生成 χ² 样本)
    // χ²(ν) 样本: 通过 ν 个独立 N(0,1)^2 之和 (标准方法)
    std::vector<Real> sample_uniforms(Philox4x64& rng) const {
        const Size n = n_names();
        std::vector<Real> u_indep(n);
        for (Size i = 0; i < n; ++i) {
            uint64_t r = rng();
            u_indep[i] = (r >> 11) * (1.0 / 9007199254740992.0);
            if (u_indep[i] < 1e-15) u_indep[i] = 1e-15;
            if (u_indep[i] > 1.0 - 1e-15) u_indep[i] = 1.0 - 1e-15;
        }
        // 生成 χ²(ν) 样本: 若 ν 为整数, 用 Σ Z_k^2; 否则用 Gamma(ν/2, 2)
        Real w = sample_chi_squared_(rng, nu_);
        return transform_uniforms(u_indep, w);
    }

    // 从相关均匀变量得到违约时间
    std::vector<Real> uniforms_to_default_times(
            const std::vector<Real>& u_corr,
            const std::vector<CreditCurve>& credit_curves) const {
        return gaussian_.uniforms_to_default_times(u_corr, credit_curves);
    }

    // 一步采样违约时间
    std::vector<Real> sample_default_times(
            Philox4x64& rng,
            const std::vector<CreditCurve>& credit_curves) const {
        auto u = sample_uniforms(rng);
        return uniforms_to_default_times(u, credit_curves);
    }

private:
    GaussianCopula gaussian_;
    Real nu_;

    // χ²(ν) 采样: 若 ν 为整数, 用 Σ Z_k^2; 否则用 Marsaglia-Tsang Gamma 采样
    static Real sample_chi_squared_(Philox4x64& rng, Real nu) {
        // 整数自由度: 直接 Σ Z_k^2 (k = ν)
        Real nu_int = std::round(nu);
        if (std::abs(nu - nu_int) < 1e-10 && nu_int >= 1.0) {
            Real w = 0.0;
            Size k = static_cast<Size>(nu_int);
            Size remaining = k;
            while (remaining > 0) {
                // 每次 box_muller 产生两个独立 N(0,1)
                Real u1 = uniform_from_uint_(rng());
                Real u2 = uniform_from_uint_(rng());
                if (u1 < 1e-15) u1 = 1e-15;
                if (u2 < 1e-15) u2 = 1e-15;
                auto [z1, z2] = box_muller(u1, u2);
                w += z1 * z1;
                --remaining;
                if (remaining > 0) {
                    w += z2 * z2;
                    --remaining;
                }
            }
            return w;
        }
        // 非整数自由度: Gamma(nu/2, 2) = χ²(ν)
        return sample_gamma_(rng, nu * 0.5, 2.0);
    }

    static Real uniform_from_uint_(uint64_t r) {
        Real u = (r >> 11) * (1.0 / 9007199254740992.0);
        if (u < 1e-15) u = 1e-15;
        if (u > 1.0 - 1e-15) u = 1.0 - 1e-15;
        return u;
    }

    // Marsaglia-Tsang Gamma 采样 (shape > 1 推荐, shape ≤ 1 用 boost)
    static Real sample_gamma_(Philox4x64& rng, Real shape, Real scale) {
        if (shape < 1.0) {
            // Boost: Gamma(shape) = Gamma(shape+1) * U^(1/shape)
            Real g = sample_gamma_(rng, shape + 1.0, 1.0);
            uint64_t r = rng();
            Real u = (r >> 11) * (1.0 / 9007199254740992.0);
            if (u < 1e-15) u = 1e-15;
            return g * std::pow(u, 1.0 / shape) * scale;
        }
        // Marsaglia-Tsang
        Real d = shape - 1.0 / 3.0;
        Real c = 1.0 / std::sqrt(9.0 * d);
        while (true) {
            Real x, v;
            // 生成标准正态
            uint64_t r1 = rng();
            uint64_t r2 = rng();
            Real u1 = (r1 >> 11) * (1.0 / 9007199254740992.0);
            Real u2 = (r2 >> 11) * (1.0 / 9007199254740992.0);
            if (u1 < 1e-15) u1 = 1e-15;
            if (u2 < 1e-15) u2 = 1e-15;
            auto [z1, z2] = box_muller(u1, u2);
            x = z1;
            v = 1.0 + c * x;
            if (v <= 0.0) continue;
            v = v * v * v;
            Real u = (rng() >> 11) * (1.0 / 9007199254740992.0);
            if (u < 1e-300) continue;
            if (u < 1.0 - 0.0331 * x * x * x * x) return d * v * scale;
            if (std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v))) return d * v * scale;
        }
    }

    // 学生 t CDF (ν 自由度)
    // 使用不完全贝塔函数: t_cdf(x, ν) = 1 - 0.5 * I_{ν/(ν+x²)}(ν/2, 1/2)  (x>0)
    //                               = 0.5 * I_{ν/(ν+x²)}(ν/2, 1/2)        (x<0)
    static Real student_t_cdf_(Real x, Real nu) {
        if (nu <= 0.0) {
            throw std::invalid_argument("student_t_cdf: nu must be positive");
        }
        // ν → ∞ 退化为标准正态
        if (nu > 1e6) return normal_cdf(x);
        // ν = 1 为柯西分布
        if (nu == 1.0) return 0.5 + std::atan(x) / PI;

        Real x2 = x * x;
        Real z = nu / (nu + x2);
        Real ib = incomplete_beta_(z, nu * 0.5, 0.5);
        return (x >= 0.0) ? (1.0 - 0.5 * ib) : (0.5 * ib);
    }

    // 不完全贝塔函数 I_x(a, b) — 使用连分数展开 (Numerical Recipes)
    static Real incomplete_beta_(Real x, Real a, Real b) {
        if (x <= 0.0) return 0.0;
        if (x >= 1.0) return 1.0;
        // 对称性: I_x(a,b) = 1 - I_{1-x}(b,a)
        // 使用连分数
        const Real tiny = 1e-30;
        const Real eps = 1e-14;
        const Size max_iter = 300;

        Real bt = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b)
                            + a * std::log(x) + b * std::log(1.0 - x));
        if (x < (a + 1.0) / (a + b + 2.0)) {
            // 正向连分数
            return bt * betacf_(x, a, b, tiny, eps, max_iter) / a;
        } else {
            // 反向连分数
            return 1.0 - bt * betacf_(1.0 - x, b, a, tiny, eps, max_iter) / b;
        }
    }

    static Real betacf_(Real x, Real a, Real b, Real tiny, Real eps, Size max_iter) {
        Real qab = a + b;
        Real qap = a + 1.0;
        Real qam = a - 1.0;
        Real c = 1.0;
        Real d = 1.0 - qab * x / qap;
        if (std::abs(d) < tiny) d = tiny;
        d = 1.0 / d;
        Real h = d;
        for (Size m = 1; m <= max_iter; ++m) {
            Real m2 = 2.0 * static_cast<Real>(m);
            Real aa = static_cast<Real>(m) * (b - static_cast<Real>(m)) * x
                      / ((qam + m2) * (a + m2));
            d = 1.0 + aa * d;
            if (std::abs(d) < tiny) d = tiny;
            c = 1.0 + aa / c;
            if (std::abs(c) < tiny) c = tiny;
            d = 1.0 / d;
            h *= d * c;
            aa = -(a + static_cast<Real>(m)) * (qab + static_cast<Real>(m)) * x
                 / ((a + m2) * (qap + m2));
            d = 1.0 + aa * d;
            if (std::abs(d) < tiny) d = tiny;
            c = 1.0 + aa / c;
            if (std::abs(c) < tiny) c = tiny;
            d = 1.0 / d;
            Real del = d * c;
            h *= del;
            if (std::abs(del - 1.0) < eps) break;
        }
        return h;
    }
};

// ============ 便捷工厂 ============

// 构建等相关的单因子 Gaussian Copula 相关矩阵
inline std::vector<std::vector<Real>> make_equicorrelation(Size n, Real rho) {
    if (rho < -1.0 / static_cast<Real>(n - 1 > 0 ? n - 1 : 1) || rho >= 1.0) {
        throw std::invalid_argument("make_equicorrelation: invalid rho for given n");
    }
    std::vector<std::vector<Real>> corr(n, std::vector<Real>(n, rho));
    for (Size i = 0; i < n; ++i) corr[i][i] = 1.0;
    return corr;
}

// 从相关矩阵构造 Gaussian Copula (便捷)
inline GaussianCopula make_gaussian_copula(const std::vector<std::vector<Real>>& correlation) {
    return GaussianCopula(correlation);
}

// 构造 t-Copula
inline TCopula make_t_copula(const std::vector<std::vector<Real>>& correlation,
                               Real degrees_of_freedom) {
    return TCopula(correlation, degrees_of_freedom);
}

}  // namespace v1
}  // namespace cpphub
