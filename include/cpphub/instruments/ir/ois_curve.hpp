#pragma once
// SOURCE: Hull (2018) "Options, Futures, and Other Derivatives" Ch.7, Ch.9
// SOURCE: Henrard (2014) "Interest Rate Modelling in the Multi-Curve Framework"
// SOURCE: Ametrano & Bianchetti (2013) "Everything You Always Wanted to Know About Multiple Interest Rate Curve Bootstrapping but Were Afraid to Ask"
// 模块: 零息债曲线 + OIS Bootstrap (post-crisis discounting 基础设施)
//
// 2008 金融危机后, LIBOR 不再是无风险利率代理, 隔夜指数 (SOFR/ESTR/SONIA) 成为贴现基准.
// OIS (Overnight Index Swap): 固定利率 K vs 隔夜指数复合平均 (几何平均, 实际按日复利).
//
// 零息债曲线 P(0,T) 三种插值方式:
//   - LinearZero: r(T) 线性插值, P(T) = exp(-r(T)*T), 远期利率阶梯状
//   - LogLinearDF: ln P(T) 线性插值, 远期利率在段内常数
//   - CubicSplineZero: r(T) 自然三次样条, 远期利率平滑
//
// OIS Bootstrap (单一曲线自举):
//   - OIS 浮端 (基于 OIS curve 自身): PV = P(0,T_start) - P(0,T_end)  (FRN par 性质)
//   - OIS 固端: PV = K * Σ τ_j * P(0, T_j)
//   - 令 PV_payer = 0 (Float = Fixed) 解出 P(0, T_i):
//     P(0,T_i) = [1 - K * Σ_{j<i} τ_j * P(0,T_j)] / (1 + K * τ_i)
//   - 逐点 bootstrap: 第 i 个 quote 仅依赖前 i-1 个已知的 P(0,T_j)

#include "cpphub/core/types.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cpphub {
inline namespace v1 {

// ============ 零息债曲线 ============
class ZeroCurve {
public:
    enum class InterpType {
        LinearZero,        // r(T) 线性插值
        LogLinearDF,       // ln P(T) 线性插值 (远期常数)
        CubicSplineZero    // r(T) 自然三次样条
    };

    ZeroCurve() = default;

    ZeroCurve(std::vector<Real> maturities,
              std::vector<Real> zero_rates,
              InterpType interp = InterpType::LinearZero)
        : T_(std::move(maturities)), r_(std::move(zero_rates)), interp_(interp) {
        validate();
        if (interp_ == InterpType::CubicSplineZero) build_cubic_spline();
    }

    // 折现因子 P(0,T) = exp(-r(T) * T)
    Real discount_factor(Real T) const {
        if (T <= 0.0) return 1.0;
        Real r = zero_rate(T);
        return std::exp(-r * T);
    }

    // 零息债收益率 r(T) (连续复利)
    Real zero_rate(Real T) const {
        if (T <= 0.0) return r_.empty() ? 0.0 : r_[0];
        if (T_.size() == 1) return r_[0];
        return interpolate_zero(T);
    }

    // 简单复利远期利率 f(T1, T2) = (P(T1)/P(T2) - 1) / (T2 - T1)
    Real forward_rate(Real T1, Real T2) const {
        if (T2 <= T1) {
            throw std::invalid_argument("ZeroCurve::forward_rate: require T2 > T1");
        }
        Real P1 = discount_factor(T1);
        Real P2 = discount_factor(T2);
        return (P1 / P2 - 1.0) / (T2 - T1);
    }

    // 连续复利远期 (瞬时) f(T) = -∂ ln P / ∂T = r(T) + T * r'(T)
    Real instantaneous_forward(Real T) const {
        if (T <= 0.0) return r_.empty() ? 0.0 : r_[0];
        Real r = zero_rate(T);
        Real dr = zero_rate_derivative(T);
        return r + T * dr;
    }

    // ============ 兼容 short_rate.hpp model 接口 (用于复用 IRS/CapFloor/Swaption) ============
    Real zero_coupon_bond(Real T) const { return discount_factor(T); }
    Real yield(Real T) const { return zero_rate(T); }

    // 访问器
    const std::vector<Real>& maturities() const noexcept { return T_; }
    const std::vector<Real>& zero_rates() const noexcept { return r_; }
    InterpType interp_type() const noexcept { return interp_; }
    Size size() const noexcept { return T_.size(); }
    bool empty() const noexcept { return T_.empty(); }

private:
    std::vector<Real> T_;
    std::vector<Real> r_;
    InterpType interp_ = InterpType::LinearZero;
    std::vector<Real> spline_y2_;  // cubic spline 二阶导数 (自然样条: y2[0]=y2[n-1]=0)

    void validate() const {
        if (T_.size() != r_.size()) {
            throw std::invalid_argument("ZeroCurve: maturities and zero_rates size mismatch");
        }
        if (T_.empty()) {
            throw std::invalid_argument("ZeroCurve: at least 1 point required");
        }
        for (Size i = 0; i < T_.size(); ++i) {
            if (T_[i] <= 0.0) {
                throw std::invalid_argument("ZeroCurve: maturities must be positive");
            }
            if (i > 0 && T_[i] <= T_[i - 1]) {
                throw std::invalid_argument("ZeroCurve: maturities must be strictly increasing");
            }
        }
    }

    // 三次样条预处理: 求二阶导数 y2[i] (Numerical Recipes 自然样条)
    void build_cubic_spline() {
        Size n = T_.size();
        spline_y2_.assign(n, 0.0);
        if (n < 3) return;  // 少于 3 点退化为线性

        std::vector<Real> u(n - 1, 0.0);
        spline_y2_[0] = 0.0;  // 自然样条下边界
        u[0] = 0.0;

        for (Size i = 1; i < n - 1; ++i) {
            Real h_prev = T_[i] - T_[i - 1];
            Real h_next = T_[i + 1] - T_[i];
            Real sig = h_prev / (T_[i + 1] - T_[i - 1]);
            Real p = sig * spline_y2_[i - 1] + 2.0;
            spline_y2_[i] = (sig - 1.0) / p;
            u[i] = (6.0 * ((r_[i + 1] - r_[i]) / h_next - (r_[i] - r_[i - 1]) / h_prev)
                    / (T_[i + 1] - T_[i - 1]) - sig * u[i - 1]) / p;
        }
        spline_y2_[n - 1] = 0.0;  // 自然样条上边界
        for (Size i = n - 1; i-- > 0; ) {
            spline_y2_[i] = spline_y2_[i] * spline_y2_[i + 1] + u[i];
        }
    }

    // 二分查找 T 所在区间 [T_[i], T_[i+1]]
    Size find_index(Real T) const {
        if (T <= T_[0]) return 0;
        if (T >= T_.back()) return T_.size() - 1;
        // 二分
        Size lo = 0, hi = T_.size() - 1;
        while (hi - lo > 1) {
            Size mid = (lo + hi) / 2;
            if (T_[mid] <= T) lo = mid;
            else hi = mid;
        }
        return lo;
    }

    // 插值求 r(T)
    Real interpolate_zero(Real T) const {
        if (T <= T_[0]) {
            // 外推 (线性)
            if (T_.size() < 2) return r_[0];
            Real slope = (r_[1] - r_[0]) / (T_[1] - T_[0]);
            return r_[0] + slope * (T - T_[0]);
        }
        if (T >= T_.back()) {
            // 外推 (线性)
            Size n = T_.size();
            Real slope = (r_[n - 1] - r_[n - 2]) / (T_[n - 1] - T_[n - 2]);
            return r_[n - 1] + slope * (T - T_[n - 1]);
        }
        Size i = find_index(T);
        Real h = T_[i + 1] - T_[i];
        Real t = (T - T_[i]) / h;

        switch (interp_) {
            case InterpType::LinearZero:
                return r_[i] + t * (r_[i + 1] - r_[i]);
            case InterpType::LogLinearDF: {
                // ln P 线性, P = exp(-r*T), ln P = -r*T
                Real lnP_i = -r_[i] * T_[i];
                Real lnP_i1 = -r_[i + 1] * T_[i + 1];
                Real lnP = lnP_i + t * (lnP_i1 - lnP_i);
                return -lnP / T;
            }
            case InterpType::CubicSplineZero: {
                if (spline_y2_.empty()) return r_[i] + t * (r_[i + 1] - r_[i]);
                Real A = (1.0 - t);
                Real B = t;
                Real C = (A * A * A - A) * h * h / 6.0;
                Real D = (B * B * B - B) * h * h / 6.0;
                return A * r_[i] + B * r_[i + 1] + C * spline_y2_[i] + D * spline_y2_[i + 1];
            }
        }
        return r_[i] + t * (r_[i + 1] - r_[i]);
    }

    // r'(T) 数值微分 (用于瞬时远期), 中心差分
    Real zero_rate_derivative(Real T) const {
        Real eps = 1e-6;
        Real r_plus = interpolate_zero(T + eps);
        Real r_minus = interpolate_zero(std::max(T - eps, T_.empty() ? 0.0 : T_[0] * 0.5));
        return (r_plus - r_minus) / (2.0 * eps);
    }
};

// ============ OIS Bootstrap ============
// 输入 OIS 报价, 逐点 bootstrap 出零息曲线.
// 假设: OIS 固定端按年支付 (annual), 浮动端基于 OIS curve 自身 (单曲线框架).
class OISCurveBootstrapper {
public:
    struct OISQuote {
        Real maturity;  // OIS 期限 (年)
        Real rate;      // OIS swap rate (年化, 简单复利)
    };

    // day_count_fraction: 每期年化分数 (1.0 = ACT/365 等价于按年)
    OISCurveBootstrapper(std::vector<OISQuote> quotes, Real day_count_fraction = 1.0)
        : quotes_(std::move(quotes)), dcf_(day_count_fraction) {
        validate();
    }

    // 从 OIS 报价 bootstrap 零息曲线
    // 假设: 第 i 个 OIS maturity = i * dcf_ (annual coupon), 形成连续 annual term structure
    ZeroCurve bootstrap(ZeroCurve::InterpType interp = ZeroCurve::InterpType::LinearZero) const {
        std::vector<Real> maturities;
        std::vector<Real> zero_rates;
        maturities.reserve(quotes_.size());
        zero_rates.reserve(quotes_.size());

        // 逐点 bootstrap: P(0,T_i) = [1 - K_i * Σ_{j<i} τ * P(0,T_j)] / (1 + K_i * τ)
        // 其中 Σ_{j<i} τ * P(0,T_j) = annuity_{i-1} (前 i-1 期的固定端 annuity)
        Real annuity = 0.0;
        Real prev_T = 0.0;
        for (Size i = 0; i < quotes_.size(); ++i) {
            Real T_i = quotes_[i].maturity;
            Real K_i = quotes_[i].rate;
            Real tau_i = (i == 0) ? T_i - prev_T : T_i - maturities.back();
            // 解出 P(0, T_i)
            Real P_i = (1.0 - K_i * annuity) / (1.0 + K_i * tau_i);
            if (P_i <= 0.0 || P_i > 1.0 + 1e-10) {
                throw std::runtime_error("OISCurveBootstrapper: bootstrap produced invalid P(0,T)");
            }
            Real r_i = -std::log(P_i) / T_i;  // 连续复利零息利率
            maturities.push_back(T_i);
            zero_rates.push_back(r_i);
            annuity += tau_i * P_i;
        }
        return ZeroCurve(std::move(maturities), std::move(zero_rates), interp);
    }

    const std::vector<OISQuote>& quotes() const noexcept { return quotes_; }
    Real day_count_fraction() const noexcept { return dcf_; }

private:
    std::vector<OISQuote> quotes_;
    Real dcf_;

    void validate() const {
        if (quotes_.empty()) {
            throw std::invalid_argument("OISCurveBootstrapper: at least 1 quote required");
        }
        for (Size i = 0; i < quotes_.size(); ++i) {
            if (quotes_[i].maturity <= 0.0) {
                throw std::invalid_argument("OISCurveBootstrapper: maturity must be positive");
            }
            if (quotes_[i].rate < 0.0) {
                throw std::invalid_argument("OISCurveBootstrapper: rate must be non-negative");
            }
            if (i > 0 && quotes_[i].maturity <= quotes_[i - 1].maturity) {
                throw std::invalid_argument("OISCurveBootstrapper: maturities must be strictly increasing");
            }
        }
    }
};

// ============ 工厂函数: 构造平坦零息曲线 ============
inline ZeroCurve make_flat_curve(Real rate, Real max_maturity = 30.0, Size n_points = 31) {
    std::vector<Real> maturities, zero_rates;
    maturities.reserve(n_points);
    zero_rates.reserve(n_points);
    for (Size i = 0; i < n_points; ++i) {
        Real T = (i == 0) ? 0.01 : static_cast<Real>(i);
        if (T > max_maturity) break;
        maturities.push_back(T);
        zero_rates.push_back(rate);
    }
    return ZeroCurve(std::move(maturities), std::move(zero_rates));
}

}  // namespace v1
}  // namespace cpphub
