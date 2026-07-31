#pragma once
// SOURCE: Pykhtin & Zhu (2007) "A Guide to Modeling Counterparty Credit Risk" (Gaussian Copula WWR)
// SOURCE: Hull & White (2012) "CVA and Wrong Way Risk" (近似公式)
// SOURCE: Gregory (2015) "XVA" Ch.15-17
// SOURCE: BCBS (2015) "Margin requirements for non-centrally cleared derivatives"
// 模块: CVA with Wrong-Way Risk (WWR) — 违约-暴露相关性
//
// ==================== 数学框架 ====================
//
// 基础 CVA (xva.hpp) 假设违约与暴露独立:
//   CVA_independent = -(1-R_c) * ∫ EPE_disc(t) * dPD_c(t)
//
// WWR 场景下, 交易对手违约往往发生在暴露上升时 (如卖出看跌期权 → 标的下跌 →
// 交易对手违约 → 暴露上升), 需对违约-暴露相关性建模.
//
// 单因子 Gaussian Copula 模型:
//   暴露驱动因子 M_V ~ N(0,1), 违约驱动因子 M_D ~ N(0,1), corr(M_V, M_D) = rho_ww
//   暴露 V(t) = f(M_V, t), 违约时间 τ = inf{t : PD(t) > Φ(M_D)}
//   (即 τ = PD^{-1}(Φ(M_D)), 若 Φ(M_D) ≤ PD(T) 则在期限内违约)
//
// 方法:
//   1. Hull-White 近似 (半解析): CVA_wwr ≈ CVA_ind * A(rho, sigma_V, PD_T)
//      A = exp(rho*sigma_V*sqrt(T)*Φ^{-1}(PD_T) + 0.5*(rho*sigma_V*sqrt(T))²)
//   2. Gaussian Copula (Pykhtin-Zhu 2007, 半解析):
//      对数正态暴露 ln V ~ N(mu_V(t), sigma_V(t)²), sigma_V(t) = sigma_V*sqrt(t)
//      条件化 ln V | M_D=m ~ N(mu_V + rho*sigma_V(t)*m, sigma_V(t)²*(1-rho²))
//      E[max(V,0)|M_D=m] = E[max(V,0)] * exp(rho*sigma_V(t)*m - 0.5*rho²*sigma_V(t)²) * Φ(d1)
//      CVA_wwr = -(1-R) * ∫_0^{PD_T} EPE_disc(τ(u)) * ratio(Φ^{-1}(u)) du
//      其中 u = Φ(M_D), τ(u) = PD^{-1}(u); 用 20 节点 Gauss-Hermite 求积
//   3. MC 模拟: 每路径生成相关高斯 [M_V, M_D], M_V 驱动资产价格,
//      M_D 决定违约时间, 若 τ ≤ T 记损失 (1-R)*max(V(τ),0)*P_d(0,τ)

#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"  // normal_cdf / normal_pdf / inv_normal_cdf
#include "cpphub/instruments/ir/ois_curve.hpp"       // ZeroCurve
#include "cpphub/instruments/credit/credit_curve.hpp" // PDCurve
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/risk/xva.hpp"  // ExposureProfile / XVAConfig / XVAResult / compute_xva
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <limits>

namespace cpphub {
inline namespace v1 {

// ============ WWR 配置 ============
// rho_ww > 0: WWR (违约时暴露高, CVA 损失放大)
// rho_ww < 0: RWR (Right-Way Risk, 违约时暴露低, CVA 损失缩小)
struct WWRConfig {
    Real rho_ww = 0.0;  // 暴露-违约相关性 ∈ [-0.99, 0.99]
};

namespace detail {

// ============ 20 节点 Gauss-Hermite (概率论者约定) ============
// ∫ f(x) φ(x) dx ≈ Σ_i w_i f(x_i),  x_i 为标准正态节点
// 节点/权重由 H_20 根 (Newton 迭代) + w = 2^{n-1} n! √π / (n² H_{n-1}(x)²) 得到,
// 再做 m = √2·x, w = w/√π 变换为标准正态约定.
inline const std::vector<std::pair<Real, Real>>& gauss_hermite_nodes_weights() {
    static const std::vector<std::pair<Real, Real>> gw = {
        {-7.6190485416797591e+00, 1.2578006724379317e-13},
        {-6.5105901570136551e+00, 2.4820623623151797e-10},
        {-5.5787388058932015e+00, 6.1274902599829531e-08},
        {-4.7345813340460552e+00, 4.4021210902308510e-06},
        {-3.9439673506573163e+00, 1.2882627996192944e-04},
        {-3.1890148165533900e+00, 1.8301031310804894e-03},
        {-2.4586636111723679e+00, 1.3997837447101012e-02},
        {-1.7452473208141270e+00, 6.1506372063976862e-02},
        {-1.0429453488027511e+00, 1.6173933398399989e-01},
        {-3.4696415708135597e-01, 2.6079306344955483e-01},
        { 3.4696415708135597e-01, 2.6079306344955483e-01},
        { 1.0429453488027511e+00, 1.6173933398399989e-01},
        { 1.7452473208141270e+00, 6.1506372063976862e-02},
        { 2.4586636111723679e+00, 1.3997837447101012e-02},
        { 3.1890148165533900e+00, 1.8301031310804894e-03},
        { 3.9439673506573163e+00, 1.2882627996192944e-04},
        { 4.7345813340460552e+00, 4.4021210902308510e-06},
        { 5.5787388058932015e+00, 6.1274902599829531e-08},
        { 6.5105901570136551e+00, 2.4820623623151797e-10},
        { 7.6190485416797591e+00, 1.2578006724379317e-13},
    };
    return gw;
}

// 分段线性插值 y(t), 网格外钳制到端点
inline Real interpolate_clamped(const std::vector<Real>& x,
                                const std::vector<Real>& y, Real t) {
    if (x.empty()) return 0.0;
    if (t <= x.front()) return y.front();
    if (t >= x.back()) return y.back();
    for (Size i = 0; i + 1 < x.size(); ++i) {
        if (t >= x[i] && t <= x[i + 1]) {
            Real w = (t - x[i]) / (x[i + 1] - x[i]);
            return y[i] + w * (y[i + 1] - y[i]);
        }
    }
    return y.back();
}

// 违约时间反演: 求 t ∈ [0, T] 使得 default_prob(t) = p
// (default_prob 单调递增 → 二分查找)
inline Real default_time_for_prob(const PDCurve& pd, Real p, Real T) {
    if (T <= 0.0) return 0.0;
    p = std::max(std::min(p, 1.0 - 1e-15), 0.0);
    if (p <= 0.0) return 0.0;
    const Real pd_T = pd.default_prob(T);
    if (p >= pd_T) return T;
    if (pd_T <= 0.0) return T;  // 永无违约
    Real lo = 0.0, hi = T;
    for (int i = 0; i < 100; ++i) {
        const Real mid = 0.5 * (lo + hi);
        if (pd.default_prob(mid) < p) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

// ============ 对数正态条件暴露比率 ============
// ln V(t) ~ N(mu_V(t), sigma_t²),  sigma_t = sigma_V*sqrt(t)
// 条件化 ln V | M_D=m ~ N(mu_V + rho*sigma_t*m, sigma_t²*(1-rho²)):
//   E[max(V,0)|M_D=m] = E[max(V,0)] * exp(mu_cond + 0.5*sigma_cond² - (mu_V + 0.5*sigma_t²))
//                      * Φ(d1_cond)
//   mu_cond    = mu_V + rho*sigma_t*m
//   sigma_cond = sigma_t*sqrt(1-rho²)
//   d1_cond    = (mu_cond + sigma_cond² - ln K_adj)/sigma_cond
// 暴露恒正 (K_adj → 0) 时 Φ(d1_cond) → 1, 得:
//   ratio(m) = exp(rho*sigma_t*m - 0.5*rho²*sigma_t²)
inline Real conditional_exposure_ratio(Real m, Real sigma_t, Real rho) {
    const Real a = rho * sigma_t;
    return std::exp(a * m - 0.5 * a * a);
}

}  // namespace detail

// ============ 方法 1: Hull-White WWR 近似 (半解析) ============
// CVA_wwr ≈ CVA_independent * exp(rho*sigma_V*sqrt(T)*Φ^{-1}(PD_T) + 0.5*(rho*sigma_V*sqrt(T))²)
inline XVAResult compute_cva_wwr_hw(
        const ExposureProfile& profile,
        const PDCurve& pd_counterparty,
        const PDCurve& pd_self,
        const XVAConfig& xva_cfg,
        const WWRConfig& wwr_cfg,
        Real sigma_V,
        Real risk_free_price = 0.0) {
    if (profile.size() < 2) {
        throw std::invalid_argument("compute_cva_wwr_hw: exposure profile needs >= 2 points");
    }
    if (sigma_V < 0.0) {
        throw std::invalid_argument("compute_cva_wwr_hw: sigma_V must be non-negative");
    }
    XVAResult base = compute_xva(profile, pd_counterparty, pd_self, xva_cfg, risk_free_price);
    if (base.cva == 0.0) {
        return base;  // 无信用损失 (PD=0 或 R=1), WWR 无影响
    }
    const Real T = profile.times.back();
    const Real PD_T = std::max(std::min(pd_counterparty.default_prob(T), 1.0 - 1e-12), 1e-12);
    const Real a = wwr_cfg.rho_ww * sigma_V * std::sqrt(T);
    const Real factor = std::exp(a * inv_normal_cdf(PD_T) + 0.5 * a * a);
    XVAResult r = base;
    r.cva = base.cva * factor;
    r.bva = r.cva + r.dva + r.fva;
    r.adjusted_price = risk_free_price + r.bva;
    return r;
}

// ============ 方法 2: Gaussian Copula WWR (Pykhtin-Zhu 2007, 半解析) ============
// CVA_wwr = -(1-R_c) * ∫_0^{PD_T} EPE_disc(τ(u)) * ratio(Φ^{-1}(u)) du
// 其中 u = Φ(M_D), τ(u) = PD^{-1}(u) 为违约时间, ratio 为对数正态条件暴露比率.
// 通过换元 u = PD_T * Φ(y) 映射到全实轴, 用 20 节点 Gauss-Hermite 求积.
inline XVAResult compute_cva_wwr_copula(
        const ExposureProfile& profile,
        const PDCurve& pd_counterparty,
        const PDCurve& pd_self,
        const XVAConfig& xva_cfg,
        const WWRConfig& wwr_cfg,
        Real sigma_V,
        Real risk_free_price = 0.0) {
    if (profile.size() < 2) {
        throw std::invalid_argument("compute_cva_wwr_copula: exposure profile needs >= 2 points");
    }
    if (sigma_V < 0.0) {
        throw std::invalid_argument("compute_cva_wwr_copula: sigma_V must be non-negative");
    }
    XVAResult base = compute_xva(profile, pd_counterparty, pd_self, xva_cfg, risk_free_price);
    const Real L = 1.0 - xva_cfg.recovery_counterparty;
    const Real T = profile.times.back();
    const Real PD_T = pd_counterparty.default_prob(T);

    Real cva = 0.0;
    if (PD_T > 0.0 && L > 0.0) {
        const Real rho = wwr_cfg.rho_ww;
        const auto& gw = detail::gauss_hermite_nodes_weights();
        Real integral = 0.0;
        for (const auto& [y, w] : gw) {
            // u = PD_T * Φ(y) ∈ [0, PD_T], 违约时间 τ(u) = PD^{-1}(u)
            const Real u = PD_T * normal_cdf(y);
            const Real tau = detail::default_time_for_prob(pd_counterparty, u, T);
            const Real epe_disc = detail::interpolate_clamped(profile.times, profile.epe, tau);
            const Real sigma_t = sigma_V * std::sqrt(tau);
            const Real m = inv_normal_cdf(u);  // Φ^{-1}(u): 违约时 M_D 阈值
            const Real ratio = detail::conditional_exposure_ratio(m, sigma_t, rho);
            integral += w * epe_disc * ratio;
        }
        cva = -L * PD_T * integral;
    }

    XVAResult r = base;
    r.cva = cva;
    r.bva = r.cva + r.dva + r.fva;
    r.adjusted_price = risk_free_price + r.bva;
    return r;
}

// ============ 方法 3: MC WWR 模拟 ============
// 单因子 Gaussian copula 违约时间模型:
//   1. 每路径生成 M_V, M_D ~ N(0,1), corr(M_V, M_D) = rho_ww
//   2. 资产路径由 M_V 驱动: S(t) = S0 * exp((r-q-0.5σ²)t + σ*sqrt(t)*M_V)
//   3. 违约时间 τ = inf{t : PD(t) > Φ(M_D)} = PD^{-1}(Φ(M_D)) (若 Φ(M_D) ≤ PD(T))
//   4. 若 τ ≤ T, 损失 = (1-R_c) * max(V(τ),0) * P_d(0,τ)
//   CVA_wwr = -E[损失]
// sigma_V 为暴露波动率参数 (本方法由 value_fn 直接计算暴露, 无需使用)
inline Real compute_cva_wwr_mc(
        const MultiAssetGBMPathGenerator& gen,
        std::function<Real(Real, const std::vector<Real>&)> value_fn,
        const std::vector<Real>& exposure_times,
        const PDCurve& pd_counterparty,
        const ZeroCurve& discount_curve,
        const XVAConfig& xva_cfg,
        const WWRConfig& wwr_cfg,
        Real sigma_V,
        Size n_paths,
        uint64_t seed) {
    (void)sigma_V;
    if (n_paths == 0) {
        throw std::invalid_argument("compute_cva_wwr_mc: n_paths must be positive");
    }
    if (exposure_times.empty()) {
        throw std::invalid_argument("compute_cva_wwr_mc: exposure_times empty");
    }
    const auto& cfg = gen.config();
    if (cfg.n_assets() < 1) {
        throw std::invalid_argument("compute_cva_wwr_mc: generator needs at least 1 asset");
    }
    const Real T = exposure_times.back();
    if (T > cfg.T + 1e-9) {
        throw std::invalid_argument("compute_cva_wwr_mc: exposure horizon exceeds generator T");
    }
    const Real S0 = cfg.S0[0];
    const Real sigma = cfg.sigma[0];
    const Real q = cfg.dividend(0);
    const Real r = cfg.r;
    const Real L = 1.0 - xva_cfg.recovery_counterparty;
    const Real rho = wwr_cfg.rho_ww;
    const Real sqrt_one_minus_rho2 = std::sqrt(std::max(0.0, 1.0 - rho * rho));

    Real sum_loss = 0.0;
    for (Size p = 0; p < n_paths; ++p) {
        Philox4x64 rng(seed, static_cast<uint64_t>(p));
        const Real M_V = next_normal(rng);
        const Real M_D = rho * M_V + sqrt_one_minus_rho2 * next_normal(rng);

        // 违约时间 τ = inf{t : PD(t) > Φ(M_D)}; 若期限内不违约则跳过
        const Real u = normal_cdf(M_D);
        if (u >= pd_counterparty.default_prob(T)) continue;
        const Real tau = detail::default_time_for_prob(pd_counterparty, u, T);

        // M_V 驱动资产价格 (单因子)
        const Real S_tau = S0 * std::exp((r - q - 0.5 * sigma * sigma) * tau
                                         + sigma * std::sqrt(tau) * M_V);
        const Real V_tau = value_fn(tau, std::vector<Real>{S_tau});
        sum_loss += L * std::max(V_tau, 0.0) * discount_curve.discount_factor(tau);
    }
    return -sum_loss / static_cast<Real>(n_paths);
}

// ============ 辅助: 从 MC 路径估计暴露波动率 ============
// 对每个时间 t_j, 计算 ln(max(V, eps)) 的横截面标准差并除以 sqrt(t_j) 年化,
// 最后对时间网格取平均. 若暴露为对数正态 ln V ~ N(·, σ²t), 结果 ≈ σ.
inline Real estimate_exposure_volatility(
        const std::vector<std::vector<Real>>& V_samples,
        const std::vector<Real>& times) {
    if (V_samples.empty() || V_samples[0].empty()) {
        throw std::invalid_argument("estimate_exposure_volatility: empty samples");
    }
    const Size n_paths = V_samples.size();
    const Size n_times = times.size();
    for (const auto& v : V_samples) {
        if (v.size() != n_times) {
            throw std::invalid_argument("estimate_exposure_volatility: column size mismatch");
        }
    }
    const Real eps = 1e-12;
    Real sum = 0.0;
    Size count = 0;
    for (Size j = 0; j < n_times; ++j) {
        if (times[j] <= 0.0) continue;
        Real mean = 0.0;
        for (Size p = 0; p < n_paths; ++p) {
            mean += std::log(std::max(V_samples[p][j], eps));
        }
        mean /= static_cast<Real>(n_paths);
        Real var = 0.0;
        for (Size p = 0; p < n_paths; ++p) {
            const Real d = std::log(std::max(V_samples[p][j], eps)) - mean;
            var += d * d;
        }
        var /= static_cast<Real>(n_paths - 1);
        sum += std::sqrt(std::max(var, 0.0)) / std::sqrt(times[j]);
        ++count;
    }
    return (count > 0) ? sum / static_cast<Real>(count) : 0.0;
}

}  // namespace v1
}  // namespace cpphub
