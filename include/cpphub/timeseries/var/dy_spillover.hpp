// =============================================================================
// dy_spillover.hpp - Diebold-Yilmaz 溢出指数: TCI/TO/FROM/NET + 滚动 (spec §4.3)
//
// Phase 7C v1.7 M2 (PHASE7C_SPEC.md v1.2 §4.3, V10/§13-a)
//
// 输入: θ̃ = GFEVD DY 框架 (行归一化, H 步, 默认 H=10 可配置 —
//   V10: 论文日度惯例非数学常数; H=10 vs H=50 敏感性测试)
//
// 指数定义 (DY 2012 IJF; 对齐 R Spillover G.spillover standardized 表):
//   TCI    = 100 · Σ_{i≠j} θ̃_ij / K          (总溢出; Spillover 表尾
//            "C. to others (spillover)" 总和 = TCI, 实测 20.3467 对齐)
//   TO_j   = 100 · Σ_{i≠j} θ̃_ij / K          (j 发出 = 列和去对角/K)
//   FROM_j = 100 · Σ_{j'≠j} θ̃_j j' / K       (j 接收 = 行和去对角/K)
//   NET_j  = TO_j − FROM_j
//   (Spillover standardized 表 = θ̃_ij·100/K, 故 TO/FROM/TCI 均带 /K)
//
// 滚动窗口 (§13-a 裁决):
//   window 必填 (> 2·K 强制, 无硬编码默认); step=1; 每窗口完整重估
//   VAR + GFEVD; tci_path/net_path 按窗口尾时间索引 (T = window..T 全长)
//
// 主基准: R Spillover g.fevd/G.spillover (1e-8); 滚动对照 roll.spillover
// =============================================================================

#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include "cpphub/core/types.hpp"
#include "cpphub/timeseries/var/fevd.hpp"
#include "cpphub/timeseries/var/multivariate_data.hpp"
#include "cpphub/timeseries/var/var_model.hpp"

namespace cpphub {
inline namespace v1 {
namespace timeseries {
namespace var {

/// DY 溢出结果 (spec §4.3)
struct DYResult {
    Real tci = 0.0;                                    ///< 总溢出 (单窗口或全样本)
    std::vector<Real> to_spillover, from_spillover, net_spillover;  ///< 按变量序
    std::vector<Real> tci_path;                        ///< 逐窗口 TCI (step=1)
    std::vector<std::vector<Real>> net_path;           ///< 逐窗口 × 变量
    Size window = 0, horizon = 10;                     ///< 回显
    Size lag = 0;                                      ///< 实际滞后阶
    Eigen::MatrixXd fevd;                              ///< 末窗口 θ̃ (附)
};

namespace detail {

/// 从 θ̃ (行归一 GFEVD) 计算 TCI/TO/FROM/NET (全部 ×100/K 口径)
inline void dy_indices_from_theta(const Eigen::MatrixXd& theta, Real& tci,
                                  std::vector<Real>& to, std::vector<Real>& from,
                                  std::vector<Real>& net) {
    const Size K = theta.rows();
    to.assign(K, 0.0);
    from.assign(K, 0.0);
    net.assign(K, 0.0);
    Real total = 0.0;
    for (Size i = 0; i < K; ++i) {
        for (Size j = 0; j < K; ++j) {
            if (i == j) continue;
            to[j] += theta(i, j);
            from[i] += theta(i, j);
            total += theta(i, j);
        }
    }
    const Real Kd = static_cast<Real>(K);
    for (Size j = 0; j < K; ++j) {
        to[j] *= 100.0 / Kd;
        from[j] *= 100.0 / Kd;
        net[j] = to[j] - from[j];
    }
    tci = 100.0 * total / Kd;
}

}  // namespace detail

/// 静态 (全样本) DY 溢出指数 (无滚动)
inline DYResult dy_spillover_static(const MultivariateTSData& data,
                                    Size horizon = 10,
                                    const std::string& trend = "c",
                                    Size lag = 0) {
    data.validate();
    VARSpec spec;
    spec.lag = lag;
    spec.trend = trend;
    VARResult fit = var_fit(data, spec);
    FEVDResult f = var_fevd(fit, horizon, FevdFramework::GeneralizedDY);

    DYResult res;
    detail::dy_indices_from_theta(f.fevd, res.tci, res.to_spillover,
                                  res.from_spillover, res.net_spillover);
    res.horizon = horizon;
    res.lag = fit.lag;
    res.fevd = f.fevd;
    return res;
}

/// 滚动 DY 溢出指数 (window 必填 > 2·K, §13-a; step=1 每窗口全重估).
/// lag=0 ⇒ 每窗口 IC (aic) 自动选阶 (V10: H 可配; trend 透传).
inline DYResult dy_spillover(const MultivariateTSData& data, Size window,
                             Size horizon = 10, const std::string& trend = "c",
                             Size lag = 0, Size seed = 42) {
    data.validate();
    const Size T = data.T(), K = data.K();
    if (window <= 2 * K) {
        throw std::invalid_argument("dy_spillover: window must be > 2*K");
    }
    if (window > T) {
        throw std::invalid_argument("dy_spillover: window > T");
    }
    (void)seed;  // 确定性路径 (IC 选阶无随机性); 参数保留对齐 spec 签名

    DYResult res;
    res.window = window;
    res.horizon = horizon;

    const Size n_win = T - window + 1;
    res.tci_path.reserve(n_win);
    res.net_path.reserve(n_win);

    for (Size w = 0; w < n_win; ++w) {
        // 窗口 [w, w+window)
        MultivariateTSData sub;
        sub.columns.resize(K);
        for (Size j = 0; j < K; ++j) {
            sub.columns[j].assign(data.columns[j].begin() + w,
                                  data.columns[j].begin() + w + window);
        }
        VARSpec spec;
        spec.lag = lag;
        spec.trend = trend;
        VARResult fit = var_fit(sub, spec);
        FEVDResult f = var_fevd(fit, horizon, FevdFramework::GeneralizedDY);
        Real tci;
        std::vector<Real> to, from, net;
        detail::dy_indices_from_theta(f.fevd, tci, to, from, net);
        res.tci_path.push_back(tci);
        res.net_path.push_back(net);
        res.lag = fit.lag;
        // 末窗口覆盖标量输出
        res.tci = tci;
        res.to_spillover = to;
        res.from_spillover = from;
        res.net_spillover = net;
        res.fevd = f.fevd;
    }
    return res;
}

}  // namespace var
}  // namespace timeseries
}  // inline namespace v1
}  // namespace cpphub
