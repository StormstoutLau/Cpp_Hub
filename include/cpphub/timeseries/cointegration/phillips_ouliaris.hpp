// =============================================================================
// phillips_ouliaris.hpp - Phillips-Ouliaris 协整检验 (Pu / Pz 双实现)
//
// Phase 7C v1.7 M3 (PHASE7C_SPEC.md §5.3; 决策 20)
//
// 教材锚点: Phillips-Ouliaris 1990 OBES 52(1):165-193;
//   对照: urca ca.po (demean/lag/type 参数配对, 1e-8)
//
// 算法 (复刻 urca ca.po 源码, ca_po_source.txt 逐行):
//   z = [y0, y1]; zl = z[1:], zr = z[:-1]  (水平对滞后水平, 非差分)
//   res = resid(zl ~ zr [+1/+trd])         (ari3: none/const/trend)
//   Ω = Bartlett LRV(res, lmax):
//     Ω = res'res/nobs + (1/nobs)·Σ_l (1−l/(lmax+1))·(Γ̂_l + Γ̂_l')
//     权重约定与 unit_root_common::long_run_variance 一致 (对角元相等)
//   Pu (残差基, 方向依赖): resu = resid(z₀ ~ z₁ [+det]);
//     w112 = Ω₁₁ − Ω₂₁'Ω₂₂⁻¹Ω₂₁; stat = nobs·w112/(Σresu²/nobs)
//   Pz (方向无关, 对协整向量归一化不变): Mzz = zl'zl/nobs;
//     stat = nobs·tr(Ω·Mzz⁻¹)
//   lmax: lag=0 → trunc(4·(nobs/100)^{1/4}) (urca "short");
//         lag≥1 → 显式带宽 (urca "long" = trunc(12·(nobs/100)^{1/4}) 由调用方
//         传入数值; 跨库对照 §1.4-4)
//   CV 表: urca ca.po 内嵌 PO1990 表 (m=2 行), 10/5/1%
//
// 排幻觉点 (spec §9.3): CI12 (Pu 方向依赖 / Pz 方向无关, Pz 优先)
// =============================================================================

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace cointegration {

/// PO 检验结果 (§5.3 v1.2 签名; 属性名镜像 urca)
struct POResult {
    Real statistic = 0.0;
    Real cv_5pct = 0.0;
    std::string type;           ///< "Pu"/"Pz" 回显 (Pz 优先, CI12)
    std::string demean = "none";  ///< "none"/"constant"/"trend"
    Size lag = 0;               ///< Bartlett 带宽 lmax 回显
    Size n_obs = 0;             ///< T−1
    bool reject_null = false;   ///< statistic > cv_5pct (右尾, H0: 无协整)
};

namespace detail {

// urca ca.po 内嵌 PO1990 临界值表 (m=2 行; 全表 m≤6, 本 API 双变量仅用首行)
// [ari3-1][pct]: ari3 1=none/2=const/3=trend; pct 0=10%/1=5%/2=1%
inline constexpr double PO_PU_CV[3][3] = {
    {20.3933, 25.9711, 38.3413},   // none
    {27.8536, 33.7130, 48.0021},   // constant
    {41.2488, 48.8439, 65.1714}};  // trend
inline constexpr double PO_PZ_CV[3][3] = {
    {33.9267, 40.8217, 55.1911},   // none
    {47.5877, 55.2202, 71.9273},   // constant
    {71.9586, 81.3812, 102.0167}}; // trend
static_assert(PO_PU_CV[0][1] == 25.9711, "PO anchor: Pu none 5pct");
static_assert(PO_PZ_CV[0][2] == 55.1911, "PO anchor: Pz none 1pct");
static_assert(PO_PZ_CV[2][1] == 81.3812, "PO anchor: Pz trend 5pct");

inline Size po_ari3(const std::string& demean) {
    if (demean == "none") return 1;
    if (demean == "constant") return 2;
    if (demean == "trend") return 3;
    throw std::invalid_argument(
        "phillips_ouliaris: demean must be none/constant/trend");
}

// Bartlett 核长期协方差矩阵 (urca ca.po Ω; 权重 1−l/(lmax+1))
// res: nobs×m (行对齐 zl); 返回 m×m
inline Eigen::MatrixXd po_omega(const Eigen::MatrixXd& res, Size lmax) {
    const Size nobs = static_cast<Size>(res.rows());
    const Eigen::Index m = res.cols();
    Eigen::MatrixXd omega = res.transpose() * res /
                            static_cast<Real>(nobs);
    for (Size l = 1; l <= lmax; ++l) {
        const Real w = 1.0 - static_cast<Real>(l) /
                                 static_cast<Real>(lmax + 1);
        const Eigen::Index li = static_cast<Eigen::Index>(l);
        const Eigen::MatrixXd g = res.bottomRows(nobs - l).transpose() *
                                  res.topRows(nobs - l);
        omega += w * (g + g.transpose()) / static_cast<Real>(nobs);
    }
    return omega;
}

// OLS 残差 (y − X·lstsq)
inline Eigen::MatrixXd po_resid(const Eigen::MatrixXd& y,
                                const Eigen::MatrixXd& x) {
    Eigen::MatrixXd beta = x.colPivHouseholderQr().solve(y);
    return y - x * beta;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Phillips-Ouliaris 协整检验 (双变量)
//
// @param y0 左侧变量 (Pu 方向依赖; Pz 对变量次序不变 — CI12 断言点)
// @param y1 右侧变量
// @param type "Pz" (方向无关, 优先) / "Pu" (残差基)
// @param demean "none"/"constant"/"trend" (urca 直译)
// @param lag 0 → urca "short" 带宽 trunc(4·(nobs/100)^0.25);
//   ≥1 → 显式 Bartlett 带宽 (urca "long" = trunc(12·(nobs/100)^0.25) 传数值)
// ---------------------------------------------------------------------------
inline POResult phillips_ouliaris(const std::vector<Real>& y0,
                                  const std::vector<Real>& y1,
                                  const std::string& type = "Pz",
                                  const std::string& demean = "none",
                                  Size lag = 0) {
    if (type != "Pu" && type != "Pz") {
        throw std::invalid_argument(
            "phillips_ouliaris: type must be Pu/Pz");
    }
    if (y0.size() != y1.size()) {
        throw std::invalid_argument("phillips_ouliaris: length mismatch");
    }
    const Size t = y0.size();
    if (t < 5) {
        throw std::invalid_argument("phillips_ouliaris: sample too small");
    }
    for (Size i = 0; i < t; ++i) {
        if (std::isnan(y0[i]) || std::isnan(y1[i])) {
            throw std::invalid_argument("phillips_ouliaris: NaN in input");
        }
    }
    const Size ari3 = detail::po_ari3(demean);

    Eigen::MatrixXd z(t, 2);
    for (Size i = 0; i < t; ++i) {
        z(i, 0) = y0[i];
        z(i, 1) = y1[i];
    }
    const Size nobs = t - 1;
    const Eigen::MatrixXd zl = z.bottomRows(
        static_cast<Eigen::Index>(nobs));
    const Eigen::MatrixXd zr = z.topRows(
        static_cast<Eigen::Index>(nobs));

    // res: zl ~ zr [+1 +trd]  (ari3: 1 无确定项, 2 常数, 3 常数+趋势)
    // 列数 = 2 (zr) + (ari3−1) 确定项: none→2 / const→3 / trend→4 [zr,1,trd]
    Eigen::MatrixXd x_reg(nobs, ari3 + 1);
    x_reg << zr, Eigen::MatrixXd::Ones(nobs, ari3 - 1);
    if (ari3 == 3) {
        for (Size i = 0; i < nobs; ++i) {
            x_reg(i, 3) = static_cast<Real>(i + 1);  // trd = 1..nobs
        }
    }
    const Eigen::MatrixXd res = detail::po_resid(zl, x_reg);

    // lmax
    Size lmax;
    if (lag == 0) {
        lmax = static_cast<Size>(std::trunc(
            4.0 * std::pow(static_cast<Real>(nobs) / 100.0, 0.25)));
    } else {
        lmax = lag;
    }

    const Eigen::MatrixXd omega = detail::po_omega(res, lmax);

    POResult out;
    out.type = type;
    out.demean = demean;
    out.lag = lmax;
    out.n_obs = nobs;

    if (type == "Pz") {
        Eigen::MatrixXd mzz = zl.transpose() * zl /
                              static_cast<Real>(nobs);
        out.statistic = static_cast<Real>(nobs) *
                        (omega * mzz.inverse()).trace();
        out.cv_5pct = detail::PO_PZ_CV[ari3 - 1][1];
    } else {
        // Pu: 首列对余列回归 (全样本 T; trend 情形 trd = 1..nobs+1 = 1..T)
        // 列数 = 1 (z₂) + (ari3−1) 确定项: none→1 / const→2 / trend→3 [z₂,1,trd]
        Eigen::MatrixXd xu(t, ari3);
        xu << z.col(1), Eigen::MatrixXd::Ones(t, ari3 - 1);
        if (ari3 == 3) {
            for (Size i = 0; i < t; ++i) {
                xu(i, 2) = static_cast<Real>(i + 1);
            }
        }
        const Eigen::MatrixXd resu = detail::po_resid(z.col(0), xu);
        const Real w11 = omega(0, 0);
        const Real w21 = omega(1, 0);
        const Real o22 = omega(1, 1);
        const Real w112 = w11 - w21 * (1.0 / o22) * w21;
        const Real ssqr = resu.squaredNorm() / static_cast<Real>(nobs);
        out.statistic = static_cast<Real>(nobs) * w112 / ssqr;
        out.cv_5pct = detail::PO_PU_CV[ari3 - 1][1];
    }
    out.reject_null = out.statistic > out.cv_5pct;
    return out;
}

}  // namespace cointegration
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
