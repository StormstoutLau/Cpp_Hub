#pragma once
// SOURCE: Pykhtin & Zhu (2007) "A Guide to Modeling Counterparty Credit Risk"
// SOURCE: Gregory (2015) "XVA: Credit, Funding and Capital Valuation Adjustments"
// SOURCE: Hull & White (2012) "The FVA Debate" + (2014) "Valuation Adjustments"
// SOURCE: Brigo, Pallavicini & Papatheodorou (2013) "FVA, CVA and the multiple curves"
// 模块: xVA (X-Value Adjustments) 框架
//
// ==================== 数学框架 ====================
//
// 设 V(t) 为衍生品组合在 t 时刻的市值 (从银行视角, 即 bank 收钱为正),
// 假设违约时间 τ 与暴露 V(t) 独立 (post-crisis 主流简化假设).
//
// 1. 暴露 (Exposure):
//    EPE(t) = E[ max(V(t), 0) ]    (Expected Positive Exposure)
//    ENE(t) = E[ max(-V(t), 0) ]   (Expected Negative Exposure)
//    折现: EPE_disc(t) = EPE(t) * P_d(0, t), 其中 P_d 为 OIS discount curve
//
// 2. CVA (Counterparty Credit VA):
//    CVA = -(1 - R_c) * ∫_0^T EPE_disc(t) * dPD_c(t)
//    离散 (trapezoidal): CVA = -(1-R_c) * Σ_i 0.5*(EPE_disc_i + EPE_disc_{i+1}) * ΔPD_c_i
//    CVA ≤ 0 (counterparty 违约导致 bank 损失)
//
// 3. DVA (Debit VA — 自身违约):
//    DVA = +(1 - R_self) * ∫_0^T ENE_disc(t) * dPD_self(t)
//    DVA ≥ 0 (bank 自身违约获益, by accounting standard IFRS 13)
//
// 4. FVA (Funding VA, Hull-White symmetric 框架):
//    FVA = -∫_0^T [EPE(t) - ENE(t)] * s_f(t) * P_d(0,t) * dt
//    其中 s_f(t) 为 funding spread (借入/借出对称假设),
//    EPE-ENE 为 net exposure (正暴露净额).
//    FVA ≤ 0 当 net funding cost (正暴露主导).
//
// 5. Bilateral VA:  BVA = CVA + DVA + FVA
//    调整后价格 = risk_free_price + BVA
//
// ============ PD 曲线 (违约概率曲线) ============
// 生存概率 Q(0,T) = exp(-∫_0^T h(s) ds), h(s) 为 hazard rate
// 分段常数 hazard rate: 给定 (T_i, h_i), 在 [T_i, T_{i+1}) 内 h=h_i
// PD(0, T) = 1 - Q(0, T);  PD(t1, t2) = Q(0, t1) - Q(0, t2)
//
// 从 CDS spread s 单点近似: h ≈ s / (1 - R)  (假设违约回收率 R 已知, flat CDS curve)
//
// 注: PDCurve 与 CreditCurve 定义于 instruments/credit/credit_curve.hpp,
//     本模块通过 include 复用, 避免定义重复.

#include "cpphub/core/types.hpp"
#include "cpphub/instruments/ir/ois_curve.hpp"
#include "cpphub/instruments/credit/credit_curve.hpp"  // PDCurve (信用模块基础设施)
#include "cpphub/pricing/monte_carlo/multi_asset_path_generator.hpp"
#include "cpphub/core/rng.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <functional>

namespace cpphub {
inline namespace v1 {

// ============ 暴露轮廓 (Exposure Profile) ============
// 沿时间网格 t_0=0 < t_1 < ... < t_n=T 计算的折现 EPE/ENE
struct ExposureProfile {
    std::vector<Real> times;   // 时间网格, n+1 点 (含 t_0=0)
    std::vector<Real> epe;     // EPE_disc(t_i) = E[max(V(t_i),0)] * P_d(0,t_i), size n+1
    std::vector<Real> ene;     // ENE_disc(t_i) = E[max(-V(t_i),0)] * P_d(0,t_i), size n+1

    Size size() const noexcept { return times.size(); }
    bool empty() const noexcept { return times.empty(); }
};

// ============ 从 MC 路径样本计算暴露轮廓 ============
// V_samples[path_idx][time_idx] = 第 path 条路径在 time_idx 处的衍生品市值 V(t)
// (已折现或未折现均可, 内部按 discount_curve 折现)
// 假设 V_samples 的列数 == times.size()
inline ExposureProfile compute_exposure(
        const std::vector<std::vector<Real>>& V_samples,
        const std::vector<Real>& times,
        const ZeroCurve& discount_curve) {
    if (V_samples.empty()) {
        throw std::invalid_argument("compute_exposure: V_samples empty");
    }
    const Size n_paths = V_samples.size();
    const Size n_times = times.size();
    for (const auto& v : V_samples) {
        if (v.size() != n_times) {
            throw std::invalid_argument("compute_exposure: V_samples column size mismatch");
        }
    }
    ExposureProfile prof;
    prof.times = times;
    prof.epe.assign(n_times, 0.0);
    prof.ene.assign(n_times, 0.0);
    for (Size j = 0; j < n_times; ++j) {
        Real sum_pos = 0.0;
        Real sum_neg = 0.0;
        for (Size p = 0; p < n_paths; ++p) {
            Real v = V_samples[p][j];
            if (v > 0.0) sum_pos += v;
            else sum_neg += -v;
        }
        Real P_d = discount_curve.discount_factor(times[j]);
        prof.epe[j] = (sum_pos / static_cast<Real>(n_paths)) * P_d;
        prof.ene[j] = (sum_neg / static_cast<Real>(n_paths)) * P_d;
    }
    return prof;
}

// ============ xVA 配置 ============
struct XVAConfig {
    Real recovery_counterparty = 0.40;  // 交易对手违约回收率 R_c
    Real recovery_self = 0.40;          // 银行自身违约回收率 R_self
    Real funding_spread = 0.005;        // s_f (年化, 50bp 默认), flat term structure
};

// ============ xVA 结果 ============
struct XVAResult {
    Real cva = 0.0;            // CVA (≤ 0)
    Real dva = 0.0;            // DVA (≥ 0)
    Real fva = 0.0;            // FVA (符号取决于净暴露)
    Real bva = 0.0;            // BVA = CVA + DVA + FVA
    Real risk_free_price = 0.0;  // 无违约/无 funding 调整的衍生品价格
    Real adjusted_price = 0.0;   // risk_free_price + BVA
};

// ============ 从暴露轮廓计算 xVA ============
// CVA = -(1-R_c) * Σ_i 0.5*(EPE_i + EPE_{i+1}) * ΔPD_c_i
// DVA = +(1-R_self) * Σ_i 0.5*(ENE_i + ENE_{i+1}) * ΔPD_self_i
// FVA = -s_f * Σ_i 0.5*((EPE_i-ENE_i) + (EPE_{i+1}-ENE_{i+1})) * Δt_i
// (EPE/ENE 已折现, FVA 中 s_f*Δt 为 funding cost over interval)
inline XVAResult compute_xva(const ExposureProfile& prof,
                                const PDCurve& pd_counterparty,
                                const PDCurve& pd_self,
                                const XVAConfig& cfg,
                                Real risk_free_price = 0.0) {
    if (prof.size() < 2) {
        throw std::invalid_argument("compute_xva: exposure profile needs >= 2 points");
    }
    XVAResult r;
    r.risk_free_price = risk_free_price;
    Real L_c = 1.0 - cfg.recovery_counterparty;   // LGD counterparty
    Real L_s = 1.0 - cfg.recovery_self;           // LGD self
    Real s_f = cfg.funding_spread;

    Real cva = 0.0, dva = 0.0, fva = 0.0;
    for (Size i = 0; i + 1 < prof.size(); ++i) {
        Real t_lo = prof.times[i];
        Real t_hi = prof.times[i + 1];
        Real dt = t_hi - t_lo;
        if (dt <= 0.0) continue;

        //梯形积分权重
        Real avg_epe = 0.5 * (prof.epe[i] + prof.epe[i + 1]);
        Real avg_ene = 0.5 * (prof.ene[i] + prof.ene[i + 1]);
        // 区间违约概率 (条件存活到 t_lo 后在 (t_lo, t_hi] 违约)
        Real dPD_c = pd_counterparty.default_prob(t_lo, t_hi);
        Real dPD_s = pd_self.default_prob(t_lo, t_hi);

        // CVA: 损失 = LGD * 暴露 * PD (从 bank 视角为负)
        cva -= L_c * avg_epe * dPD_c;
        // DVA: 银行自身违约 → 债务减记 → bank 收益
        dva += L_s * avg_ene * dPD_s;
        // FVA: net funding cost = s_f * (EPE - ENE) * dt, 银行付出 (符号为负)
        fva -= s_f * (avg_epe - avg_ene) * dt;
    }
    r.cva = cva;
    r.dva = dva;
    r.fva = fva;
    r.bva = cva + dva + fva;
    r.adjusted_price = risk_free_price + r.bva;
    return r;
}

// ============ 端到端: MC 模拟计算 xVA ============
// value_fn(t, S_t) 返回衍生品在 t 时刻的市值 V(t) (从 bank 视角).
// 注意: V(t) 应为 t 时刻未折现的市场价值, 内部按 discount_curve 折现.
//       typical 约定: V(t) = E_t^Q[payoff * P(t, T_pay)] (即 t 时刻的 PV)
//
// 实现流程:
//   1. 用 path generator 模拟 n_paths 条路径
//   2. 对每条路径在每个 t_i 计算 V(t_i) = value_fn(t_i, S(t_i))
//   3. 对 V_samples 计算 EPE/ENE profile
//   4. 用 PD 曲线 + funding spread 计算 xVA
inline XVAResult compute_xva_mc(
        const MultiAssetGBMPathGenerator& gen,
        std::function<Real(Real, const std::vector<Real>&)> value_fn,
        const std::vector<Real>& exposure_times,
        const PDCurve& pd_counterparty,
        const PDCurve& pd_self,
        const ZeroCurve& discount_curve,
        const XVAConfig& cfg,
        Size n_paths,
        uint64_t seed,
        Real risk_free_price = 0.0) {
    if (n_paths == 0) {
        throw std::invalid_argument("compute_xva_mc: n_paths must be positive");
    }
    if (exposure_times.empty()) {
        throw std::invalid_argument("compute_xva_mc: exposure_times empty");
    }
    // exposure_times 必须在 path generator 的网格内
    // path generator 生成 t = 0, dt, 2dt, ..., T = n_steps*dt
    // 通过 step_idx = round(t / dt) 查找
    const Real dt_gen = gen.dt();
    const Real T_gen = gen.config().T;
    const Size n_steps = gen.config().n_steps;

    // 校验: exposure_times 中的每个 t 必须可对应到 path generator 网格
    std::vector<Size> step_indices;
    step_indices.reserve(exposure_times.size());
    for (Real t : exposure_times) {
        if (t < 0.0 || t > T_gen + 1e-9) {
            throw std::invalid_argument("compute_xva_mc: exposure_time out of range");
        }
        Size idx = static_cast<Size>(std::round(t / dt_gen));
        if (idx > n_steps) idx = n_steps;
        step_indices.push_back(idx);
    }

    std::vector<std::vector<Real>> V_samples(n_paths,
                                               std::vector<Real>(exposure_times.size()));
    for (Size p = 0; p < n_paths; ++p) {
        Philox4x64 rng(seed, static_cast<uint64_t>(p));
        auto paths = gen.generate_path(rng, 1.0);
        // paths[asset_idx][step_idx]
        for (Size j = 0; j < exposure_times.size(); ++j) {
            Size step = step_indices[j];
            std::vector<Real> S_at_t(paths.size());
            for (Size a = 0; a < paths.size(); ++a) {
                S_at_t[a] = paths[a][step];
            }
            V_samples[p][j] = value_fn(exposure_times[j], S_at_t);
        }
    }

    auto profile = compute_exposure(V_samples, exposure_times, discount_curve);
    return compute_xva(profile, pd_counterparty, pd_self, cfg, risk_free_price);
}

// ============ 单期衍生品近似 (解析验证用) ============
// 对仅在 T 时刻产生一次性现金流的衍生品:
//   EPE_disc(T) = E[max(V(T), 0)] * P(0, T), ENE_disc(T) = E[max(-V(T), 0)] * P(0, T)
//   CVA = -(1-R_c) * EPE_disc(T) * PD_c(0, T)
//   DVA = +(1-R_self) * ENE_disc(T) * PD_self(0, T)
//   FVA = -s_f * (EPE_disc(T) - ENE_disc(T)) * T  (近似: funding cost over T)
//
// 此函数接收已折现的 epe/ene 单点值, 用于测试和快速估算
inline XVAResult compute_xva_single_period(
        Real epe_disc_T, Real ene_disc_T, Real T,
        Real pd_counterparty, Real pd_self,
        const XVAConfig& cfg, Real risk_free_price = 0.0) {
    XVAResult r;
    r.risk_free_price = risk_free_price;
    Real L_c = 1.0 - cfg.recovery_counterparty;
    Real L_s = 1.0 - cfg.recovery_self;
    r.cva = -L_c * epe_disc_T * pd_counterparty;
    r.dva = L_s * ene_disc_T * pd_self;
    // 单期 FVA 近似: funding spread * net exposure * T
    r.fva = -cfg.funding_spread * (epe_disc_T - ene_disc_T) * T;
    r.bva = r.cva + r.dva + r.fva;
    r.adjusted_price = risk_free_price + r.bva;
    return r;
}

}  // namespace v1
}  // namespace cpphub
