#pragma once
// SOURCE: PHASE3_SPEC §4.2 - Dupire local volatility recovery
// Implemented on main station (MSVC) - 2026-07-30
// Strategy: self-contained IV grid + BSM analytic call price + numerical differentiation
//   sigma^2_loc(K,T) = (dC/dT + qC + (r-q)K dC/dK) / (0.5 K^2 d^2C/dK^2)
// Ref: Dupire (1994), "Pricing with a smile"
//      Gatheral (2006), "The Volatility Surface: A Practitioner's Guide"
//
// NOTE: 也保留 VolSurface 引用构造接口 (B 站实现 VolSurface 后可对接),
//       但 IV grid 自包含路径独立可用, 不依赖 VolSurface 完整实现.
#include "cpphub/core/types.hpp"
#include "cpphub/core/math.hpp"
#include "cpphub/core/constants.hpp"
#include "cpphub/risk/greeks/greeks_analytic.hpp"
#include "cpphub/models/vol_surface/vol_surface.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace cpphub {
inline namespace v1 {

class DupireLocalVol {
public:
    // ---- 自包含构造 (推荐): IV grid + 市场参数 ----
    DupireLocalVol(const std::vector<Real>& strikes,
                   const std::vector<Real>& maturities,
                   const std::vector<std::vector<Real>>& implied_vols,
                   Real S, Real r, Real q)
        : strikes_(strikes), maturities_(maturities),
          implied_vols_(implied_vols), S_(S), r_(r), q_(q),
          has_surface_(false), surface_ptr_(nullptr) {
        validate_grid();
    }

    // ---- VolSurface 引用构造 (B 站实现后对接) ----
    DupireLocalVol(const VolSurface& surface, Real r, Real q)
        : r_(r), q_(q),
          has_surface_(true), surface_ptr_(&surface) {
        // 从 surface 提取 strikes/maturities (如果可用)
        // 当前 VolSurface 未完整实现, 此构造路径仅占位
    }

    // ---- 局部波动率 ----
    Real local_vol(Real K, Real T) const {
        Real lv2 = local_variance(K, T);
        if (lv2 < 0.0) {
            // 数值噪声导致的负方差, clamp 到 0
            return 0.0;
        }
        return std::sqrt(lv2);
    }

    // ---- 局部方差 (Dupire 公式) ----
    Real local_variance(Real K, Real T) const {
        if (K <= 0.0 || T <= 0.0) {
            throw std::invalid_argument("DupireLocalVol::local_variance: K and T must be positive");
        }

        // 选择差分步长 (相对步长, 避免数值误差)
        Real h_K = 1e-3 * K;        // K 步长
        Real h_T = 1e-4 * std::max(T, 0.1);  // T 步长

        // 中心差分计算 dC/dT, dC/dK, d²C/dK²
        Real C      = call_price_at(K, T);
        Real C_Kp   = call_price_at(K + h_K, T);
        Real C_Km   = call_price_at(K - h_K, T);
        Real C_Tp   = call_price_at(K, T + h_T);
        Real C_Tm   = call_price_at(K, T - h_T);

        Real dC_dK  = (C_Kp - C_Km) / (2.0 * h_K);
        Real d2C_dK2 = (C_Kp - 2.0 * C + C_Km) / (h_K * h_K);
        Real dC_dT  = (C_Tp - C_Tm) / (2.0 * h_T);

        // Dupire 公式: σ²_loc = (∂C/∂T + qC + (r-q)K ∂C/∂K) / (0.5 K² ∂²C/∂K²)
        Real numerator = dC_dT + q_ * C + (r_ - q_) * K * dC_dK;
        Real denominator = 0.5 * K * K * d2C_dK2;

        if (std::abs(denominator) < 1e-14) {
            // 退化情况 (e.g., 极 OTM), 返回 0
            return 0.0;
        }
        return numerator / denominator;
    }

    // ---- 网格化输出 (供 PDE/MC 使用) ----
    std::vector<std::vector<Real>> local_vol_grid(
        const std::vector<Real>& strikes,
        const std::vector<Real>& maturities) const {
        std::vector<std::vector<Real>> grid(
            maturities.size(), std::vector<Real>(strikes.size()));
        for (Size j = 0; j < maturities.size(); ++j) {
            for (Size i = 0; i < strikes.size(); ++i) {
                grid[j][i] = local_vol(strikes[i], maturities[j]);
            }
        }
        return grid;
    }

    // ---- 验证: 用恢复的局部波动率跑简化 MC, 对比 BSM 价格 ----
    // 误差 = |MC_price - BSM_price| (绝对误差)
    Real recovery_error(Real K, Real T, Real S, Size n_paths = 100000) const {
        // BSM 基准价格: 用 (K, T) 处的隐含波动率
        Real sigma_impl = implied_vol_at(K, T);
        auto g = AnalyticGreeksEngine::bsm_european(S, K, T, r_, q_, sigma_impl, true);
        Real bsm_price = g.price;

        // MC 价格: 用局部波动率模型跑 Euler SDE
        // S_{t+dt} = S_t * exp((r - 0.5*σ²_loc(S_t,t))*dt + σ_loc*sqrt(dt)*Z)
        Real mc_price = mc_local_vol_price(S, K, T, n_paths);

        return std::abs(mc_price - bsm_price);
    }

    // ---- 访问器 ----
    const std::vector<Real>& strikes() const { return strikes_; }
    const std::vector<Real>& maturities() const { return maturities_; }
    const std::vector<std::vector<Real>>& vols() const { return implied_vols_; }
    Real S() const { return S_; }
    Real r() const { return r_; }
    Real q() const { return q_; }
    const VolSurface& surface() const {
        if (!has_surface_) throw std::runtime_error("no VolSurface attached");
        return *surface_ptr_;
    }

private:
    // ---- 数据 ----
    std::vector<Real> strikes_;
    std::vector<Real> maturities_;
    std::vector<std::vector<Real>> implied_vols_;  // [T_idx][K_idx]
    Real S_ = 100.0;
    Real r_ = 0.0;
    Real q_ = 0.0;

    bool has_surface_ = false;
    const VolSurface* surface_ptr_ = nullptr;

    // ---- 验证 grid ----
    void validate_grid() const {
        if (strikes_.size() < 3) {
            throw std::invalid_argument("DupireLocalVol: need at least 3 strikes for 2nd-order diff");
        }
        if (maturities_.size() < 2) {
            throw std::invalid_argument("DupireLocalVol: need at least 2 maturities for d/dT");
        }
        if (implied_vols_.size() != maturities_.size()) {
            throw std::invalid_argument("DupireLocalVol: implied_vols rows != maturities");
        }
        for (Size j = 0; j < maturities_.size(); ++j) {
            if (implied_vols_[j].size() != strikes_.size()) {
                throw std::invalid_argument("DupireLocalVol: implied_vols cols != strikes");
            }
        }
        // 检查 strikes 严格递增
        for (Size i = 1; i < strikes_.size(); ++i) {
            if (strikes_[i] <= strikes_[i-1]) {
                throw std::invalid_argument("DupireLocalVol: strikes must be strictly increasing");
            }
        }
        for (Size j = 1; j < maturities_.size(); ++j) {
            if (maturities_[j] <= maturities_[j-1]) {
                throw std::invalid_argument("DupireLocalVol: maturities must be strictly increasing");
            }
        }
    }

    // ---- 双线性插值获取 (K, T) 处的隐含波动率 ----
    Real implied_vol_at(Real K, Real T) const {
        // 边界 clamp (避免外推)
        Real K_clamped = std::max(strikes_.front(), std::min(K, strikes_.back()));
        Real T_clamped = std::max(maturities_.front(), std::min(T, maturities_.back()));

        // 找 K 区间 [K_i, K_{i+1}]
        Size i = 0;
        while (i + 1 < strikes_.size() && strikes_[i + 1] < K_clamped) ++i;
        if (i + 1 >= strikes_.size()) i = strikes_.size() - 2;
        Real K0 = strikes_[i], K1 = strikes_[i + 1];
        Real alpha_K = (K_clamped - K0) / (K1 - K0);

        // 找 T 区间 [T_j, T_{j+1}]
        Size j = 0;
        while (j + 1 < maturities_.size() && maturities_[j + 1] < T_clamped) ++j;
        if (j + 1 >= maturities_.size()) j = maturities_.size() - 2;
        Real T0 = maturities_[j], T1 = maturities_[j + 1];
        Real alpha_T = (T_clamped - T0) / (T1 - T0);

        // 双线性插值
        Real v00 = implied_vols_[j][i];
        Real v01 = implied_vols_[j][i + 1];
        Real v10 = implied_vols_[j + 1][i];
        Real v11 = implied_vols_[j + 1][i + 1];

        Real v0 = v00 * (1.0 - alpha_K) + v01 * alpha_K;
        Real v1 = v10 * (1.0 - alpha_K) + v11 * alpha_K;
        return v0 * (1.0 - alpha_T) + v1 * alpha_T;
    }

    // ---- BSM call price at (K, T) using interpolated IV ----
    Real call_price_at(Real K, Real T) const {
        if (K <= 0.0 || T <= 0.0) return 0.0;
        Real sigma = implied_vol_at(K, T);
        if (sigma <= 0.0) return 0.0;
        auto g = AnalyticGreeksEngine::bsm_european(S_, K, T, r_, q_, sigma, true);
        return g.price;
    }

    // ---- MC local vol pricing (Euler SDE on log-spot) ----
    Real mc_local_vol_price(Real S0, Real K, Real T, Size n_paths) const {
        // 简化: 时间步固定 50 步 (足够平坦 IV 场景)
        Size n_steps = 50;
        Real dt = T / static_cast<Real>(n_steps);
        Real sqrt_dt = std::sqrt(dt);

        // xorshift64* PRNG (避免 <random> MSVC ICE 风险)
        uint64_t state = 0x9E3779B97F4A7C15ULL;
        auto next_gauss = [&]() -> Real {
            // Box-Muller
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            Real u1 = static_cast<Real>(state * 0x2545F4914F6CDD1DULL) /
                      static_cast<Real>(0xFFFFFFFFFFFFFFFFULL);
            if (u1 < 1e-12) u1 = 1e-12;
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            Real u2 = static_cast<Real>(state * 0x2545F4914F6CDD1DULL) /
                      static_cast<Real>(0xFFFFFFFFFFFFFFFFULL);
            if (u2 < 1e-12) u2 = 1e-12;
            Real mag = std::sqrt(-2.0 * std::log(u1));
            return mag * std::cos(2.0 * PI * u2);
        };

        Real sum_payoff = 0.0;
        for (Size p = 0; p < n_paths; ++p) {
            Real S = S0;
            Real t = 0.0;
            for (Size s = 0; s < n_steps; ++s) {
                Real lv = local_vol(S, t + 0.5 * dt);  // mid-step
                if (lv < 0.0) lv = 0.0;
                Real z = next_gauss();
                // log-spot Euler: dS/S = (r - q - 0.5*lv²)*dt + lv*sqrt(dt)*Z
                Real drift = (r_ - q_ - 0.5 * lv * lv) * dt;
                Real diff = lv * sqrt_dt * z;
                S = S * std::exp(drift + diff);
                t += dt;
            }
            Real payoff = S - K;
            if (payoff < 0.0) payoff = 0.0;
            sum_payoff += payoff;
        }
        Real mc_price = std::exp(-r_ * T) * sum_payoff / static_cast<Real>(n_paths);
        return mc_price;
    }
};

}  // namespace v1
}  // namespace cpphub
